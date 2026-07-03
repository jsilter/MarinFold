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
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import atlas_io
from selection import select_representatives

# MMseqs search output columns (fraction-identity ``fident`` is what we threshold on).
_SEARCH_FMT = "query,target,fident,alnlen,qcov,tcov"


def run(cmd: list[str]) -> None:
    """Run a subprocess, streaming output, raising on non-zero exit."""
    print(f"[run] {' '.join(str(c) for c in cmd)}", flush=True)
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
    meta_rows: list[dict] = []
    n_seen = n_kept = 0
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
            if n_seen % (batch_rows * 50) < batch_rows:
                print(f"[scan] seen {n_seen:,} kept {n_kept:,}", flush=True)

    pq.write_table(pa.Table.from_pylist(meta_rows),
                   work / f"survivors_meta_{shard_id}.parquet")
    print(f"[scan] done: kept {n_kept:,}/{n_seen:,} "
          f"({n_kept / max(1, n_seen):.1%}) -> {fasta_path}", flush=True)


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
    print(f"[merge] {meta.num_rows:,} survivors from {len(shards)} shard(s)")
    return fasta, meta_path


# ------------------------------------------------------------ mmseqs drop helpers

def _drop_by_search(query_fasta: Path, target: Path, work: Path, tag: str,
                    max_seq_id: float, cov: float, threads: int) -> set[str]:
    """Return query ids with a target hit ≥ ``max_seq_id`` identity (to be dropped).

    ``--min-seq-id`` makes MMseqs report only hits at/above the threshold, so any
    query appearing in the result has a too-similar match.
    """
    res = work / f"{tag}.m8"
    tmp = work / f"_tmp_{tag}"
    run(["mmseqs", "easy-search", query_fasta, target, res, tmp,
         "--min-seq-id", max_seq_id, "-c", cov, "--cov-mode", 1,
         "--format-output", _SEARCH_FMT, *_threads_arg(threads)])
    if res.stat().st_size == 0:
        return set()
    # The .m8 has one row per hit, so at full scale it can reach hundreds of GB —
    # never materialise it as a DataFrame. We only need the set of query ids that
    # matched, so stream the query column in chunks and accumulate.
    cols = _SEARCH_FMT.split(",")
    drop: set[str] = set()
    for chunk in pd.read_csv(res, sep="\t", header=None, names=cols,
                             usecols=[cols[0]], dtype=str, chunksize=5_000_000):
        drop.update(chunk[cols[0]])
    print(f"[{tag}] {len(drop):,} survivors hit at >= {max_seq_id:.0%} id -> dropped")
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
                 cov: float, threads: int) -> None:
    """Drop survivors too similar to ``ref`` (novelty or leakage)."""
    cur = work / "survivors.fasta"
    meta = pq.read_table(work / "survivors_meta.parquet").column("protein_hash")
    all_ids = set(meta.to_pylist())
    drop = _drop_by_search(cur, ref, work, tag, max_seq_id, cov, threads)
    keep = all_ids - drop
    kept = _subset_fasta(cur, keep, work / f"survivors_{tag}.fasta")
    (work / f"survivors_{tag}.fasta").replace(cur)  # advance the working fasta
    print(f"[{tag}] kept {kept:,}/{len(all_ids):,}")


# ------------------------------------------------------------------- cluster/select

def stage_cluster(work: Path, *, cluster_id: float, cov: float, threads: int) -> None:
    """linclust the surviving sequences; write ``cluster_map.csv`` (hash,cluster_id)."""
    pref = work / "clu"
    tmp = work / "_tmp_clu"
    run(["mmseqs", "easy-linclust", work / "survivors.fasta", pref, tmp,
         "--min-seq-id", cluster_id, "-c", cov, "--cov-mode", 1,
         *_threads_arg(threads)])
    tsv = pref.with_name("clu_cluster.tsv")  # mmseqs writes <pref>_cluster.tsv
    cm = pd.read_csv(tsv, sep="\t", header=None,
                     names=["cluster_id", "protein_hash"])
    cm.to_csv(work / "cluster_map.csv", index=False)
    print(f"[cluster] {cm['cluster_id'].nunique():,} clusters over {len(cm):,} seqs")


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
    print(f"[select] {len(sel):,} reps from {frame['cluster_id'].nunique():,} "
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
    args = ap.parse_args(argv)

    work = args.work_dir
    work.mkdir(parents=True, exist_ok=True)
    threads = args.threads or 0

    def do_scan() -> None:
        stage_scan(work, min_len=args.min_len, max_len=args.max_len,
                   min_plddt=args.min_plddt, min_ptm=args.min_ptm,
                   compute_plddt_std=args.compute_plddt_std,
                   num_shards=args.num_shards, shard_id=args.shard_id,
                   limit=args.limit, batch_rows=args.batch_rows)

    if args.stage in ("all", "scan"):
        do_scan()
    if args.stage in ("all", "merge"):
        _merge_shards(work)
    if args.stage in ("all", "novelty"):
        if not args.afdb_ref:
            sys.exit("--afdb-ref required for the novelty stage")
        stage_filter(work, ref=args.afdb_ref, tag="novelty",
                     max_seq_id=args.max_afdb_seq_id, cov=args.search_cov,
                     threads=threads)
    if args.stage in ("all", "leakage"):
        if not args.eval_ref:
            sys.exit("--eval-ref required for the leakage stage")
        stage_filter(work, ref=args.eval_ref, tag="leakage",
                     max_seq_id=args.eval_max_seq_id, cov=args.search_cov,
                     threads=threads)
    if args.stage in ("all", "cluster"):
        stage_cluster(work, cluster_id=args.cluster_id, cov=args.cluster_cov,
                      threads=threads)
    if args.stage in ("all", "select"):
        stage_select(work, reps_per_cluster=args.reps_per_cluster,
                     select_by=args.select_by, min_cluster_size=args.min_cluster_size)


if __name__ == "__main__":
    main()
