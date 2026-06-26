#!/usr/bin/env bash
# Copyright The MarinFold Authors
# SPDX-License-Identifier: Apache-2.0
#
# End-to-end driver for building the ESM Atlas distillation subset.
#
# Stages (pass one as $1; "all" chains prep-ref -> run -> manifest):
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
  prep-ref) prep_ref ;;
  smoke)    smoke ;;
  run)      run ;;
  manifest) manifest ;;
  all)      prep_ref; run; manifest ;;
  *) echo "usage: $0 {prep-ref|smoke|run|manifest|all}" >&2; exit 2 ;;
esac
