# Copyright The MarinFold Authors
# SPDX-License-Identifier: Apache-2.0

"""Materialize selected Atlas reps into an afdb-24M-layout parquet (stage 3).

``pipeline.py`` produces ``selected_manifest.csv`` — one row per chosen cluster
representative (``protein_hash`` + metadata), but **no structures**. This script
is the final step: it pulls those rows' ``structure_blob`` from ``folds_1B.lance``,
decodes each to mmCIF text, and writes a parquet whose schema matches the
``afdb-24M`` datasets the training pipeline already reads — each row carries the
raw mmCIF in ``cif_content`` and the id in ``entry_id`` (the
``contacts-v1`` ``generate`` defaults, see
``marinfold/document_structures/contacts_v1/parse.py``). The published parquet
then feeds ``contacts-v1 generate`` unchanged.

Like ``scan``, this is one pass over ``folds_1B.lance`` keyed on a membership set,
so it is **shardable** (``--num-shards`` / ``--shard-id``) across instances and
**must run in AWS us-west-2** — it decodes ~one ``structure_blob`` per selected
rep (the heavy column we deliberately skipped in the funnel), so the output is
multi-TB at the full ~100M scale. Each shard writes one parquet to ``--out-dir``;
``run_aws.py`` syncs that to S3. Publishing the result to HuggingFace is a
**separate, sign-off-gated** step (cross-cloud, >10 GB) — this script never does it.

Output parquet columns (per rep):
  entry_id        protein_hash (the afdb-24M id column)
  cif_content     decoded mmCIF text (per-residue pLDDT in the B-factor column)
  sequence        residue sequence (from the decoded structure)
  seq_len, mean_plddt, ptm, plddt_std, cluster_id, cluster_size   (from the manifest)
  source          constant "esm-atlas-v1" provenance tag
"""

import argparse
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import atlas_io

# Match the contacts-v1 afdb-24M defaults (parse.DEFAULT_ID_COLUMN / _CIF_COLUMN)
# so the published parquet feeds ``contacts-v1 generate`` with no column remapping.
ID_COLUMN = "entry_id"
CIF_COLUMN = "cif_content"
SOURCE_TAG = "esm-atlas-v1"

# Manifest metadata carried through onto each materialized row (everything except
# protein_hash, which becomes entry_id).
_CARRY_COLS = ["seq_len", "mean_plddt", "ptm", "plddt_std", "cluster_id",
               "cluster_size"]


def decode_to_cif(blob: bytes, name: str) -> tuple[str, str]:
    """Decode a ``structure_blob`` once to ``(mmCIF text, sequence)``.

    The mmCIF carries per-residue pLDDT in the B-factor column (``atom37_to_structure``).
    """
    struct = atlas_io.decode_structure_blob(blob)
    doc = atlas_io.atom37_to_structure(struct, name=name).make_mmcif_document()
    return doc.as_string(), struct["sequence"]


def stage_materialize(
    manifest: Path, out_dir: Path, *, num_shards: int, shard_id: int,
    batch_rows: int, limit: int | None, write_cifs: bool,
) -> None:
    """Decode the manifest's selected reps to mmCIF and write a parquet shard.

    Streams ``folds_1B.lance`` over this shard's contiguous row range and decodes
    only rows whose ``protein_hash`` is in the manifest. Metadata is taken from the
    manifest (not re-read), so only ``protein_hash`` + ``structure_blob`` are pulled
    from Lance.
    """
    sel = pd.read_csv(manifest)
    if "protein_hash" not in sel.columns:
        raise SystemExit(f"{manifest} has no 'protein_hash' column")
    meta = sel.set_index("protein_hash")
    wanted = set(meta.index.astype(str))
    print(f"[materialize] manifest: {len(wanted):,} reps from {manifest}", flush=True)

    ds = atlas_io.open_folds()
    total = ds.count_rows()
    shard = total // num_shards
    start = shard_id * shard
    stop = total if shard_id == num_shards - 1 else start + shard
    if limit is not None:
        stop = min(stop, start + limit)
    print(f"[materialize] scan rows {start:,}..{stop:,} of {total:,} "
          f"(shard {shard_id}/{num_shards})", flush=True)

    cif_dir = out_dir / "cifs"
    if write_cifs:
        cif_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    cols = ["protein_hash", "structure_blob"]
    rows: list[dict] = []
    remaining = set(wanted)  # shrinks as reps are found; lets us stop early
    n_seen = n_hit = 0
    scanner = ds.scanner(columns=cols, offset=start, limit=stop - start,
                         batch_size=batch_rows)
    for batch in scanner.to_batches():
        d = batch.to_pydict()
        for i in range(len(d["protein_hash"])):
            n_seen += 1
            h = d["protein_hash"][i]
            if h not in remaining:
                continue
            remaining.discard(h)
            cif, sequence = decode_to_cif(d["structure_blob"][i], name=h[:4])
            m = meta.loc[h]
            row = {ID_COLUMN: h, CIF_COLUMN: cif, "sequence": sequence,
                   "source": SOURCE_TAG}
            for c in _CARRY_COLS:
                if c in meta.columns:
                    row[c] = m[c]
            rows.append(row)
            if write_cifs:
                (cif_dir / f"{h}.cif").write_text(cif)
            n_hit += 1
        if n_seen % (batch_rows * 50) < batch_rows:
            print(f"[materialize] seen {n_seen:,} materialized {n_hit:,}", flush=True)
        if not remaining:  # every manifest rep found; no need to scan further
            print(f"[materialize] all {n_hit:,} reps found by row {n_seen:,}; "
                  f"stopping scan", flush=True)
            break

    out = out_dir / f"atlas_distill_{shard_id}.parquet"
    pq.write_table(pa.Table.from_pylist(rows), out)
    print(f"[materialize] wrote {n_hit:,} structures -> {out}", flush=True)
    missing = len(wanted) - n_hit
    if missing and num_shards == 1:
        print(f"[materialize] WARNING: {missing:,} manifest reps not found in the "
              f"scanned range (id mismatch or out-of-range shard?)", flush=True)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", type=Path, required=True,
                    help="selected_manifest.csv from pipeline.py")
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="directory for the parquet shard(s)")
    ap.add_argument("--num-shards", type=int, default=1)
    ap.add_argument("--shard-id", type=int, default=0)
    ap.add_argument("--batch-rows", type=int, default=50_000)
    ap.add_argument("--limit", type=int, default=None,
                    help="cap rows scanned (smoke test)")
    ap.add_argument("--write-cifs", action="store_true",
                    help="also dump one <hash>.cif per rep (debug / spot-check)")
    args = ap.parse_args(argv)

    stage_materialize(
        args.manifest, args.out_dir, num_shards=args.num_shards,
        shard_id=args.shard_id, batch_rows=args.batch_rows, limit=args.limit,
        write_cifs=args.write_cifs)


if __name__ == "__main__":
    main()
