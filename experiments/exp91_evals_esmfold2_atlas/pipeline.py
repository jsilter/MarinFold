# Copyright The MarinFold Authors
# SPDX-License-Identifier: Apache-2.0

"""On-instance worker for the ESM Atlas distillation funnel (runs in AWS us-west-2).

This is the production realization of ``funnel.py`` + ``selection.py``: a sequence +
metadata pipeline that reads the Atlas *in region* (no egress, no ``structure_blob``
decode for the bulk) and emits a manifest of selected ``protein_hash``es. Launched by
``run_aws.py``; see that file and the README for the EC2 wiring.

Stages (``--stage all`` runs them in order; each can also run standalone):

1. ``scan``     — stream ``folds_1B.lance``; keep rows passing length / mean-pLDDT /
                  pTM; write ``survivors.fasta`` + ``survivors_meta.parquet``
                  (``protein_hash, seq_len, mean_plddt, ptm[, plddt_std]``). This is
                  the one pass over all 1.1B rows; shardable via ``--num-shards`` /
                  ``--shard-id`` for parallel instances.
2. ``novelty``  — ``mmseqs easy-search`` survivors vs the AFDB reference; drop any
                  survivor with a hit ≥ ``--max-afdb-seq-id`` identity (keeps the
                  sequence-novel set). Queries with no hit are novel → kept.
3. ``leakage``  — same vs the held-out eval sequences; drop ≥ ``--eval-max-seq-id``.
4. ``cluster``  — ``mmseqs easy-linclust`` the kept survivors at ``--cluster-id``
                  (40%; the Atlas is already 70%-clustered so only <70% collapses).
5. ``select``   — one representative per cluster (``selection.select_representatives``,
                  ESMFold2 rule: longest then lowest ``plddt_std``) →
                  ``selected_manifest.csv``.

All MMseqs identity thresholds are fraction (0–1). Outputs land under ``--work-dir``;
``run_aws.py`` syncs that to S3.
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import atlas_io
from selection import select_representatives

# MMseqs search output columns (fraction-identity ``fident`` is what we threshold on).
_SEARCH_FMT = "query,target,fident,alnlen,qcov,tcov"

_T0 = time.monotonic()


def _mem_avail_gb() -> float:
    """Available RAM in GB from /proc/meminfo (-1 if unreadable)."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable"):
                    return int(line.split()[1]) / 1e6  # kB -> GB
    except OSError:
        pass
    return -1.0


def _resources(work: Path) -> str:
    """One-line disk + RAM snapshot for the work volume (for progress logging)."""
    du = shutil.disk_usage(work)
    return (f"disk {du.used / 1e12:.2f}/{du.total / 1e12:.2f} TB used, "
            f"{du.free / 1e12:.2f} TB free; RAM avail {_mem_avail_gb():.0f} GB")


def _log(msg: str) -> None:
    """Timestamped, flushed progress line (elapsed seconds since process start).

    Everything the pipeline prints goes to a log that ``run_aws.py`` ships to S3
    every couple minutes, so verbose, frequently-flushed logging is what keeps a
    long run observable instead of a black box.
    """
    print(f"[{time.monotonic() - _T0:8.1f}s] {msg}", flush=True)


def run(cmd: list[str]) -> None:
    """Run a subprocess, streaming output, raising on non-zero exit."""
    _log(f"[run] {' '.join(str(c) for c in cmd)}")
    subprocess.run([str(c) for c in cmd], check=True)


def _threads_arg(threads: int) -> list[str]:
    """``--threads N`` only when N > 0 (MMseqs rejects ``--threads 0``)."""
    return ["--threads", str(threads)] if threads and threads > 0 else []


# --------------------------------------------------------------------------- scan

def stage_scan(
    work: Path, *, min_len: int, max_len: int, min_plddt: float, min_ptm: float,
    compute_plddt_std: bool, num_shards: int, shard_id: int, limit: int | None,
    batch_rows: int,
) -> None:
    """Stream ``folds_1B.lance`` and write survivors passing the metadata gates.

    Reads only metadata columns (``protein_hash, sequence, mean_plddt, ptm`` and,
    if ``compute_plddt_std``, ``per_residue_plddt``) — never ``structure_blob``.
    Shard ``shard_id`` of ``num_shards`` processes a contiguous row range so several
    instances can run the scan in parallel.
    """
    ds = atlas_io.open_folds()
    total = ds.count_rows()
    shard = total // num_shards
    start = shard_id * shard
    stop = total if shard_id == num_shards - 1 else start + shard
    if limit is not None:
        stop = min(stop, start + limit)
    print(f"[scan] rows {start:,}..{stop:,} of {total:,} "
          f"(shard {shard_id}/{num_shards})", flush=True)

    cols = ["protein_hash", "sequence", "mean_plddt", "ptm"]
    if compute_plddt_std:
        cols.append("per_residue_plddt")

    fasta_path = work / f"survivors_{shard_id}.fasta"
    meta_path = work / f"survivors_meta_{shard_id}.parquet"
    # Flush metadata to parquet in row-group batches rather than holding every
    # survivor's dict in RAM (at ~150M survivors a single list is 100+ GB). The
    # ParquetWriter is opened lazily on the first flush so the schema matches the
    # rows actually produced (with/without plddt_std).
    meta_rows: list[dict] = []
    writer: pq.ParquetWriter | None = None
    n_seen = n_kept = 0

    def flush_meta() -> None:
        nonlocal writer, meta_rows
        if not meta_rows:
            return
        tbl = pa.Table.from_pylist(meta_rows)
        if writer is None:
            writer = pq.ParquetWriter(meta_path, tbl.schema)
        writer.write_table(tbl)
        meta_rows = []

    scanner = ds.scanner(columns=cols, offset=start, limit=stop - start,
                         batch_size=batch_rows)
    with open(fasta_path, "w") as fasta:
        for batch in scanner.to_batches():
            d = batch.to_pydict()
            for i in range(len(d["protein_hash"])):
                n_seen += 1
                seq = d["sequence"][i]
                L = len(seq)
                mp = d["mean_plddt"][i]
                pt = d["ptm"][i]
                if not (min_len <= L <= max_len and mp >= min_plddt and pt >= min_ptm):
                    continue
                h = d["protein_hash"][i]
                row = {"protein_hash": h, "seq_len": L,
                       "mean_plddt": float(mp), "ptm": float(pt)}
                if compute_plddt_std:
                    arr = atlas_io.decode_npy_zip(d["per_residue_plddt"][i]).astype(float)
                    if arr.max() > 1.5:  # stored 0–100; normalise to 0–1
                        arr = arr / 100.0
                    row["plddt_std"] = float(arr.std())
                fasta.write(f">{h}\n{seq}\n")
                meta_rows.append(row)
                n_kept += 1
                if len(meta_rows) >= 1_000_000:
                    flush_meta()
            if n_seen % (batch_rows * 20) < batch_rows:
                _log(f"[scan] seen {n_seen:,} kept {n_kept:,} "
                     f"({n_kept / max(1, n_seen):.1%}); {_resources(work)}")
    flush_meta()
    if writer is not None:
        writer.close()
    else:  # no survivors at all — still emit an empty parquet with the schema
        schema_cols = {"protein_hash": pa.string(), "seq_len": pa.int64(),
                       "mean_plddt": pa.float64(), "ptm": pa.float64()}
        if compute_plddt_std:
            schema_cols["plddt_std"] = pa.float64()
        pq.write_table(pa.table({k: pa.array([], type=t)
                                 for k, t in schema_cols.items()}), meta_path)
    _log(f"[scan] done: kept {n_kept:,}/{n_seen:,} "
         f"({n_kept / max(1, n_seen):.1%}) -> {fasta_path}")


def _merge_shards(work: Path) -> tuple[Path, Path]:
    """Concatenate per-shard FASTA + meta into ``survivors.fasta`` / ``.parquet``."""
    fasta = work / "survivors.fasta"
    shards = sorted(work.glob("survivors_*.fasta"))
    with open(fasta, "wb") as out:
        for s in shards:
            out.write(s.read_bytes())
    metas = [pq.read_table(p) for p in sorted(work.glob("survivors_meta_*.parquet"))]
    meta = pa.concat_tables(metas)
    meta_path = work / "survivors_meta.parquet"
    pq.write_table(meta, meta_path)
    _log(f"[merge] {meta.num_rows:,} survivors from {len(shards)} shard(s)")
    return fasta, meta_path


# ------------------------------------------------------------ mmseqs drop helpers

def _split_fasta(src: Path, out_dir: Path, chunk_seqs: int) -> list[Path]:
    """Split ``src`` into FASTA files of at most ``chunk_seqs`` records each.

    Streaming, one pass, constant memory (never loads the FASTA). Returns the chunk
    paths in order.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    parts: list[Path] = []
    fout = None
    n_in_chunk = 0
    with open(src) as fin:
        for line in fin:
            if line.startswith(">"):
                if fout is None or n_in_chunk >= chunk_seqs:
                    if fout is not None:
                        fout.close()
                    p = out_dir / f"chunk_{len(parts):05d}.fasta"
                    parts.append(p)
                    fout = open(p, "w")
                    n_in_chunk = 0
                n_in_chunk += 1
            fout.write(line)
    if fout is not None:
        fout.close()
    return parts


def _drop_by_search(query_fasta: Path, target_fasta: Path, work: Path, tag: str,
                    max_seq_id: float, cov: float, threads: int, *,
                    chunk_seqs: int, split_memory_limit: str | None,
                    sensitivity: str | None = None) -> set[str]:
    """Return query ids with a target hit ≥ ``max_seq_id`` identity (to drop).

    The query set is searched in **chunks of ``chunk_seqs`` sequences**, not all at
    once: a single ``easy-search`` over ~150M queries drove MMseqs into split-
    prefilter mode and produced multiple TB of scratch that overran the disk and
    wedged the run. Chunking bounds peak RAM and disk to one chunk, and we delete
    each chunk's scratch + ``.m8`` immediately after pulling the matched query ids
    out of it (the ``.m8`` would otherwise reach hundreds of GB).

    The **target DB + prefilter index are built once and reused** for every chunk.
    ``easy-search`` re-``createdb``s (and re-indexes) the whole reference on every
    call; for one fixed reference searched by many query chunks that is redundant,
    so we drop to the lower-level ``createdb``/``createindex``/``search``/
    ``convertalis`` path and keep the indexed target between chunks. Results are
    identical to ``easy-search`` (the index is just the precomputed target k-mers).

    ``--min-seq-id`` makes MMseqs report only hits at/above the threshold, so any
    query appearing in the result has a too-similar match. ``sensitivity`` (``-s``)
    and ``--split-memory-limit`` are passed through to the search when set.
    """
    cols = _SEARCH_FMT.split(",")
    # createindex and search MUST use the same -s: createindex defaults to 7.5 but
    # search defaults to 5.7, so a default-built index would silently make the
    # search more sensitive than plain easy-search. Pin both to the same value
    # (mmseqs search default 5.7 unless overridden) so results are reproducible.
    s = sensitivity if sensitivity is not None else "5.7"
    search_extra: list = ["-s", s]
    if split_memory_limit:
        search_extra += ["--split-memory-limit", split_memory_limit]

    # Build the fixed reference DB once (createdb) and reuse it for every chunk;
    # easy-search would re-createdb the whole reference on each call. We do NOT
    # createindex: a precomputed index measurably changes prefilter seeding vs
    # plain easy-search (more hits), and the target createdb — not the in-search
    # index build — is the redundant work worth removing here.
    tdir = work / f"_targetdb_{tag}"
    tdir.mkdir(parents=True, exist_ok=True)
    target_db = tdir / "targetDB"
    run(["mmseqs", "createdb", target_fasta, target_db])

    chunk_dir = work / f"_chunks_{tag}"
    parts = _split_fasta(query_fasta, chunk_dir, chunk_seqs)
    _log(f"[{tag}] query split into {len(parts)} chunk(s) of <= {chunk_seqs:,} seqs; "
         f"target DB + index prebuilt once")

    drop: set[str] = set()
    for i, part in enumerate(parts):
        t0 = time.monotonic()
        cdir = work / f"_tmp_{tag}_{i:05d}"      # per-chunk query DB + search scratch
        cdir.mkdir(parents=True, exist_ok=True)
        query_db = cdir / "queryDB"
        result_db = cdir / "resultDB"
        search_tmp = cdir / "search_tmp"
        res = work / f"{tag}_{i:05d}.m8"
        run(["mmseqs", "createdb", part, query_db])
        # --alignment-mode 3 computes exact seq id + coverage so --min-seq-id / -c
        # filter identically to easy-search (which sets it internally); without it
        # the bare `search` default lets extra sub-threshold hits through.
        run(["mmseqs", "search", query_db, target_db, result_db, search_tmp,
             "--alignment-mode", 3, "--min-seq-id", max_seq_id, "-c", cov,
             "--cov-mode", 1, *search_extra, *_threads_arg(threads)])
        run(["mmseqs", "convertalis", query_db, target_db, result_db, res,
             "--format-output", _SEARCH_FMT, *_threads_arg(threads)])
        before = len(drop)
        if res.exists() and res.stat().st_size > 0:
            for ch in pd.read_csv(res, sep="\t", header=None, names=cols,
                                  usecols=[cols[0]], dtype=str, chunksize=2_000_000):
                drop.update(ch[cols[0]])
        # Free this chunk's scratch immediately so nothing accumulates on disk.
        shutil.rmtree(cdir, ignore_errors=True)
        res.unlink(missing_ok=True)
        part.unlink(missing_ok=True)
        _log(f"[{tag}] chunk {i + 1}/{len(parts)}: +{len(drop) - before:,} hit "
             f"({len(drop):,} dropped so far) in {time.monotonic() - t0:.0f}s; "
             f"{_resources(work)}")

    shutil.rmtree(chunk_dir, ignore_errors=True)
    shutil.rmtree(tdir, ignore_errors=True)
    _log(f"[{tag}] {len(drop):,} queries hit at >= {max_seq_id:.0%} id -> dropped")
    return drop


def _subset_fasta(src: Path, keep: set[str], dst: Path) -> int:
    """Write only records whose header id is in ``keep``. Returns count kept."""
    n = 0
    with open(src) as fin, open(dst, "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                write = line[1:].strip() in keep
                if write:
                    n += 1
            if write:
                fout.write(line)
    return n


def stage_filter(work: Path, *, ref: Path, tag: str, max_seq_id: float,
                 cov: float, threads: int, chunk_seqs: int,
                 split_memory_limit: str | None,
                 sensitivity: str | None = None) -> None:
    """Drop survivors too similar to ``ref`` (novelty or leakage)."""
    cur = work / "survivors.fasta"
    meta = pq.read_table(work / "survivors_meta.parquet").column("protein_hash")
    all_ids = set(meta.to_pylist())
    _log(f"[{tag}] {len(all_ids):,} survivors in; searching vs {ref}")
    drop = _drop_by_search(cur, ref, work, tag, max_seq_id, cov, threads,
                           chunk_seqs=chunk_seqs,
                           split_memory_limit=split_memory_limit,
                           sensitivity=sensitivity)
    keep = all_ids - drop
    kept = _subset_fasta(cur, keep, work / f"survivors_{tag}.fasta")
    (work / f"survivors_{tag}.fasta").replace(cur)  # advance the working fasta
    _log(f"[{tag}] kept {kept:,}/{len(all_ids):,}")


# ------------------------------------------------------------------- cluster/select

def stage_cluster(work: Path, *, cluster_id: float, cov: float, threads: int,
                  split_memory_limit: str | None) -> None:
    """linclust the surviving sequences; write ``cluster_map.csv`` (hash,cluster_id)."""
    pref = work / "clu"
    tmp = work / "_tmp_clu"
    cmd = ["mmseqs", "easy-linclust", work / "survivors.fasta", pref, tmp,
           "--min-seq-id", cluster_id, "-c", cov, "--cov-mode", 1,
           *_threads_arg(threads)]
    if split_memory_limit:
        cmd += ["--split-memory-limit", split_memory_limit]
    run(cmd)
    shutil.rmtree(tmp, ignore_errors=True)  # linclust scratch is large + useless
    tsv = pref.with_name("clu_cluster.tsv")  # mmseqs writes <pref>_cluster.tsv
    cm = pd.read_csv(tsv, sep="\t", header=None,
                     names=["cluster_id", "protein_hash"])
    cm.to_csv(work / "cluster_map.csv", index=False)
    _log(f"[cluster] {cm['cluster_id'].nunique():,} clusters over {len(cm):,} seqs")


def stage_select(work: Path, *, reps_per_cluster: int, select_by: str,
                 min_cluster_size: int) -> None:
    """One rep per cluster -> ``selected_manifest.csv``."""
    cm = pd.read_csv(work / "cluster_map.csv")
    meta = pq.read_table(work / "survivors_meta.parquet").to_pandas()
    frame = cm.merge(meta, on="protein_hash", how="left")
    if "plddt_std" not in frame and select_by == "esmfold2":
        print("[select] no plddt_std column; falling back to select_by=length")
        select_by = "length"
    sel = select_representatives(frame, reps_per_cluster=reps_per_cluster,
                                 select_by=select_by,
                                 min_cluster_size=min_cluster_size)
    out = work / "selected_manifest.csv"
    sel.to_csv(out, index=False)
    _log(f"[select] {len(sel):,} reps from {frame['cluster_id'].nunique():,} "
         f"clusters -> {out}")


# ------------------------------------------------------------------------- main

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--stage", default="all",
                    choices=["all", "scan", "merge", "novelty", "leakage",
                             "cluster", "select"])
    # scan / metadata gates
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--max-len", type=int, default=1000)
    ap.add_argument("--min-plddt", type=float, default=0.70)
    ap.add_argument("--min-ptm", type=float, default=0.50)
    ap.add_argument("--compute-plddt-std", action="store_true",
                    help="decode per_residue_plddt for the uniform-confidence "
                         "rep-selection criterion (extra read; otherwise select by length)")
    ap.add_argument("--num-shards", type=int, default=1)
    ap.add_argument("--shard-id", type=int, default=0)
    ap.add_argument("--limit", type=int, default=None, help="cap rows (testing)")
    ap.add_argument("--batch-rows", type=int, default=100_000)
    # references
    ap.add_argument("--afdb-ref", type=Path, help="AFDB/UniRef FASTA or mmseqs DB")
    ap.add_argument("--eval-ref", type=Path, help="held-out eval sequences FASTA")
    # thresholds
    ap.add_argument("--max-afdb-seq-id", type=float, default=0.40)
    ap.add_argument("--eval-max-seq-id", type=float, default=0.40)
    ap.add_argument("--search-cov", type=float, default=0.50)
    ap.add_argument("--cluster-id", type=float, default=0.40)
    ap.add_argument("--cluster-cov", type=float, default=0.80)
    # select
    ap.add_argument("--reps-per-cluster", type=int, default=1)
    ap.add_argument("--select-by", default="esmfold2")
    ap.add_argument("--min-cluster-size", type=int, default=1)
    ap.add_argument("--threads", type=int, default=0, help="0 = mmseqs default")
    # scaling knobs for the mmseqs stages (see _drop_by_search)
    ap.add_argument("--query-chunk-seqs", type=int, default=5_000_000,
                    help="novelty/leakage: sequences per query chunk (bounds RAM + "
                         "disk; one easy-search over all ~150M queries overruns both)")
    ap.add_argument("--split-memory-limit", default=None,
                    help="mmseqs --split-memory-limit (e.g. 100G); extra RAM cap for "
                         "search + linclust. Default: mmseqs auto")
    ap.add_argument("--search-sensitivity", default=None,
                    help="mmseqs -s for novelty/leakage (default mmseqs 5.7). Lower "
                         "(e.g. 4.0) is faster but may miss some near-threshold hits")
    args = ap.parse_args(argv)

    work = args.work_dir
    work.mkdir(parents=True, exist_ok=True)
    threads = args.threads or 0
    _log(f"pipeline start: stage={args.stage} work={work}; {_resources(work)}")

    def _run_stage(name: str, fn) -> None:
        _log(f"=== stage {name}: START — {_resources(work)}")
        t0 = time.monotonic()
        fn()
        _log(f"=== stage {name}: DONE in {time.monotonic() - t0:.0f}s — "
             f"{_resources(work)}")

    if args.stage in ("all", "scan"):
        _run_stage("scan", lambda: stage_scan(
            work, min_len=args.min_len, max_len=args.max_len,
            min_plddt=args.min_plddt, min_ptm=args.min_ptm,
            compute_plddt_std=args.compute_plddt_std,
            num_shards=args.num_shards, shard_id=args.shard_id,
            limit=args.limit, batch_rows=args.batch_rows))
    if args.stage in ("all", "merge"):
        _run_stage("merge", lambda: _merge_shards(work))
    if args.stage in ("all", "novelty"):
        if not args.afdb_ref:
            sys.exit("--afdb-ref required for the novelty stage")
        _run_stage("novelty", lambda: stage_filter(
            work, ref=args.afdb_ref, tag="novelty",
            max_seq_id=args.max_afdb_seq_id, cov=args.search_cov, threads=threads,
            chunk_seqs=args.query_chunk_seqs,
            split_memory_limit=args.split_memory_limit,
            sensitivity=args.search_sensitivity))
    if args.stage in ("all", "leakage"):
        if not args.eval_ref:
            sys.exit("--eval-ref required for the leakage stage")
        _run_stage("leakage", lambda: stage_filter(
            work, ref=args.eval_ref, tag="leakage",
            max_seq_id=args.eval_max_seq_id, cov=args.search_cov, threads=threads,
            chunk_seqs=args.query_chunk_seqs,
            split_memory_limit=args.split_memory_limit,
            sensitivity=args.search_sensitivity))
    if args.stage in ("all", "cluster"):
        _run_stage("cluster", lambda: stage_cluster(
            work, cluster_id=args.cluster_id, cov=args.cluster_cov,
            threads=threads, split_memory_limit=args.split_memory_limit))
    if args.stage in ("all", "select"):
        _run_stage("select", lambda: stage_select(
            work, reps_per_cluster=args.reps_per_cluster,
            select_by=args.select_by, min_cluster_size=args.min_cluster_size))

    _log("pipeline: ALL STAGES COMPLETE")


if __name__ == "__main__":
    main()
