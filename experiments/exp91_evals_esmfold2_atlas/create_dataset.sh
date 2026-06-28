#!/usr/bin/env bash
# Copyright The MarinFold Authors
# SPDX-License-Identifier: Apache-2.0
#
# End-to-end driver for building the ESM Atlas distillation subset.
#
# Stages (pass one as $1; "all" chains prep-ref -> run -> manifest):
#   local-smoke  run the WHOLE chain locally on a tiny Atlas slice (no AWS): funnel
#                (scan->novelty->leakage->cluster->select) + materialize to parquet.
#                Needs mmseqs on PATH; reads the Atlas anonymously. Sanity check only.
#   prep-ref  download our AFDB training set from HuggingFace, extract sequences
#             to FASTA, upload to S3 (the novelty reference; run once).
#   smoke     launch a capped EC2 run (--limit) to validate the whole chain cheaply.
#   run       launch the full funnel on EC2; writes selected_manifest.csv to S3.
#   manifest  download the finished manifest locally.
#   all       prep-ref + run + manifest.
#
# Everything funnel-side runs in AWS us-west-2, co-located with the Atlas (no
# egress). See README "Building the distillation subset" + run_aws.py for details.
set -euo pipefail

# ============================ CONFIG (edit these) ============================
BUCKET="YOURBUCKET"                       # S3 bucket you own, in us-west-2
REGION="us-west-2"
IAM_PROFILE="marinfold-exp91"             # instance profile w/ S3 read+write
INSTANCE_TYPE="r7i.16xlarge"              # big-RAM CPU; scale to survivor count
SMOKE_LIMIT="2000000"                     # rows scanned in the smoke run

# --- AFDB novelty reference (downloaded from the MarinFold HF bucket) ---------
HF_BUCKET="open-athena/MarinFold"         # the MarinFold HF bucket (namespace/name)
HF_AFDB_PREFIX="afdb-24m"                 # path within the bucket w/ the seqs (adjust;
                                          #   `hf buckets list hf://buckets/${HF_BUCKET}`)
HF_AFDB_INCLUDE="*.parquet"               # narrow to the sequence parquets
AFDB_ID_COL="id"                          # id column in those parquets (adjust)
AFDB_SEQ_COL="sequence"                   # sequence column (adjust)
# Alternative: pull from a dataset repo instead of the bucket (set to 1).
USE_DATASET_REPO="0"
HF_AFDB_REPO="timodonnell/afdb-24M"

# --- funnel thresholds (mirror funnel.py / pipeline.py defaults) --------------
# pLDDT is a QUALITY FLOOR, not the size dial. 0.70 = the paper's "high confidence"
# line (~320M survive here, before clustering). Keep it loose so we don't delete
# the sole member of a novel structural cluster; clustering + per-cluster rep
# selection (CLUSTER_ID / REPS_PER_CLUSTER below) sets the final published size and
# picks the best-pLDDT exemplar per cluster. Tighten CLUSTER_ID, not this, to shrink.
MIN_PLDDT="0.70"
MIN_PTM="0.50"
MAX_AFDB_SEQ_ID="0.40"                    # novelty: drop >=40% id to AFDB
EVAL_MAX_SEQ_ID="0.40"                    # leakage: drop >=40% id to eval set
CLUSTER_ID="0.40"                         # linclust identity (size dial)
REPS_PER_CLUSTER="1"

# --- local-smoke knobs --------------------------------------------------------
LOCAL_SMOKE_LIMIT="50000"                 # Atlas rows scanned in the local smoke
LOCAL_WORK="${LOCAL_WORK:-_localsmoke}"   # scratch dir (gitignored); wiped each run
EVAL_SEQS_CSV="../exp65_evals_low_msa_depth_proteins/data/candidate_sequences.csv"
# ============================================================================

# Derived locations
STAGING="s3://${BUCKET}/exp91/staging"
OUT="s3://${BUCKET}/exp91/out"
SMOKE_OUT="s3://${BUCKET}/exp91/smoke"
AFDB_REF_S3="s3://${BUCKET}/refs/afdb_ref.fasta"
AFDB_LOCAL_DIR="${AFDB_LOCAL_DIR:-data/_afdb_ref}"     # local scratch (gitignored)
AFDB_FASTA="${AFDB_LOCAL_DIR}/afdb_ref.fasta"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

local_smoke() {  # whole pipeline on a tiny slice, no AWS — does it run + emit valid output?
  command -v mmseqs >/dev/null || { echo "[local-smoke] mmseqs not on PATH" >&2; exit 1; }
  local w="$LOCAL_WORK"
  rm -rf "$w"; mkdir -p "$w"
  echo "[local-smoke] building eval + tiny AFDB references"
  EVAL_SEQS_CSV="$EVAL_SEQS_CSV" SMOKE_DIR="$w" uv run python - <<'PY'
import os, pandas as pd
w = os.environ["SMOKE_DIR"]
df = pd.read_csv(os.environ["EVAL_SEQS_CSV"])
with open(f"{w}/eval_seqs.fasta", "w") as f:
    for r in df.itertuples():
        f.write(f">{r.stem}\n{r.sequence}\n")
# Tiny stand-in AFDB novelty reference (exercises the mmseqs search path). The
# real run uses prep-ref's full afdb_ref.fasta; for a local sanity check two
# arbitrary sequences are enough.
with open(f"{w}/afdb_ref.fasta", "w") as f:
    f.write(">ref1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR\n")
    f.write(">ref2\nMSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSI\n")
print(f"[local-smoke] eval seqs: {len(df)}, afdb ref: 2")
PY

  echo "[local-smoke] funnel (limit ${LOCAL_SMOKE_LIMIT})"
  uv run python pipeline.py --work-dir "$w/work" --stage all --limit "$LOCAL_SMOKE_LIMIT" \
    --afdb-ref "$w/afdb_ref.fasta" --eval-ref "$w/eval_seqs.fasta" --compute-plddt-std \
    --min-plddt "$MIN_PLDDT" --min-ptm "$MIN_PTM" \
    --max-afdb-seq-id "$MAX_AFDB_SEQ_ID" --eval-max-seq-id "$EVAL_MAX_SEQ_ID" \
    --cluster-id "$CLUSTER_ID" --reps-per-cluster "$REPS_PER_CLUSTER"

  # Materialize only a small, row-ordered slice of the manifest so the (heavy,
  # cross-region) structure_blob reads finish quickly; early-exit stops the scan
  # once those reps are found.
  echo "[local-smoke] mini manifest (first 30 reps) + materialize"
  SMOKE_DIR="$w" uv run python - <<'PY'
import os, pandas as pd
w = os.environ["SMOKE_DIR"]
m = pd.read_parquet(f"{w}/work/survivors_meta.parquet").head(30).copy()
m["cluster_id"] = m["protein_hash"]; m["cluster_size"] = 1
m.to_csv(f"{w}/mini_manifest.csv", index=False)
PY
  uv run python materialize.py --manifest "$w/mini_manifest.csv" \
    --out-dir "$w/dataset" --limit "$LOCAL_SMOKE_LIMIT"

  SMOKE_DIR="$w" uv run python - <<'PY'
import os, pandas as pd, gemmi
w = os.environ["SMOKE_DIR"]
df = pd.read_parquet(f"{w}/dataset/atlas_distill_0.parquet")
ok = sum(1 for c in df.cif_content if len(gemmi.read_structure_string(c)[0][0]) > 0)
assert {"entry_id", "cif_content"} <= set(df.columns), "schema mismatch"
assert ok == len(df), f"only {ok}/{len(df)} CIFs parse"
print(f"[local-smoke] OK: {len(df)} structures, schema + all CIFs valid -> {w}/dataset/")
PY
}

prep_ref() {
  mkdir -p "$AFDB_LOCAL_DIR"
  if [ "$USE_DATASET_REPO" = "1" ]; then
    echo "[prep-ref] hf download ${HF_AFDB_REPO} (${HF_AFDB_INCLUDE})"
    hf download "$HF_AFDB_REPO" --type dataset \
        --include "$HF_AFDB_INCLUDE" --local-dir "$AFDB_LOCAL_DIR"
  else
    echo "[prep-ref] hf buckets sync hf://buckets/${HF_BUCKET}/${HF_AFDB_PREFIX}"
    hf buckets sync "hf://buckets/${HF_BUCKET}/${HF_AFDB_PREFIX}" "$AFDB_LOCAL_DIR" \
        --include "$HF_AFDB_INCLUDE"
  fi

  echo "[prep-ref] extracting sequences -> ${AFDB_FASTA}"
  AFDB_LOCAL_DIR="$AFDB_LOCAL_DIR" AFDB_FASTA="$AFDB_FASTA" \
  AFDB_ID_COL="$AFDB_ID_COL" AFDB_SEQ_COL="$AFDB_SEQ_COL" \
  uv run python - <<'PY'
import os, glob, pyarrow.parquet as pq
d, out = os.environ["AFDB_LOCAL_DIR"], os.environ["AFDB_FASTA"]
idc, sc = os.environ["AFDB_ID_COL"], os.environ["AFDB_SEQ_COL"]
files = sorted(glob.glob(os.path.join(d, "**", "*.parquet"), recursive=True))
assert files, f"no parquet files under {d} (check HF_AFDB_INCLUDE)"
n = 0
with open(out, "w") as f:
    for p in files:
        for batch in pq.ParquetFile(p).iter_batches(columns=[idc, sc], batch_size=50000):
            col = batch.to_pydict()
            for i, s in zip(col[idc], col[sc]):
                f.write(f">{i}\n{s}\n"); n += 1
print(f"[prep-ref] wrote {n:,} sequences from {len(files)} parquet file(s)")
PY

  echo "[prep-ref] uploading -> ${AFDB_REF_S3}"
  aws s3 cp "$AFDB_FASTA" "$AFDB_REF_S3" --region "$REGION"
  echo "[prep-ref] done. Novelty reference ready at ${AFDB_REF_S3}"
}

launch() {  # $1 = output prefix, remaining args appended to run_aws.py
  local out_prefix="$1"; shift
  uv run python run_aws.py \
    --s3-staging "$STAGING" \
    --s3-output  "$out_prefix" \
    --afdb-ref-uri "$AFDB_REF_S3" \
    --iam-instance-profile "$IAM_PROFILE" \
    --region "$REGION" \
    --instance-type "$INSTANCE_TYPE" \
    --min-plddt "$MIN_PLDDT" --min-ptm "$MIN_PTM" \
    --max-afdb-seq-id "$MAX_AFDB_SEQ_ID" --eval-max-seq-id "$EVAL_MAX_SEQ_ID" \
    --cluster-id "$CLUSTER_ID" --reps-per-cluster "$REPS_PER_CLUSTER" \
    --compute-plddt-std --watch "$@"
}

smoke() {
  echo "[smoke] capped run (--limit ${SMOKE_LIMIT}) -> ${SMOKE_OUT}"
  launch "$SMOKE_OUT" --limit "$SMOKE_LIMIT"
  echo "[smoke] check: aws s3 ls ${SMOKE_OUT}/"
}

run() {
  echo "[run] full funnel -> ${OUT}"
  launch "$OUT"
  echo "[run] manifest -> ${OUT}/selected_manifest.csv"
}

manifest() {
  aws s3 cp "${OUT}/selected_manifest.csv" ./selected_manifest.csv --region "$REGION"
  echo "[manifest] saved ./selected_manifest.csv ($(wc -l < selected_manifest.csv) lines)"
}

case "${1:-}" in
  local-smoke) local_smoke ;;
  prep-ref)    prep_ref ;;
  smoke)       smoke ;;
  run)         run ;;
  manifest)    manifest ;;
  all)         prep_ref; run; manifest ;;
  *) echo "usage: $0 {local-smoke|prep-ref|smoke|run|manifest|all}" >&2; exit 2 ;;
esac
