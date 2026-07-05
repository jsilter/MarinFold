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
#   prep-ref  extract the afdb-24M struct-cluster rep sequences from exp41's
#             MMseqs/Foldseek target DB to FASTA, upload to S3 (novelty ref; once).
#   smoke     launch a capped EC2 run (--limit) to validate the whole chain cheaply.
#   run       launch the full funnel on EC2; writes selected_manifest.csv to S3.
#   manifest  download the finished manifest locally.
#   all       prep-ref + run + manifest.
#
# Everything funnel-side runs in AWS us-west-2, co-located with the Atlas (no
# egress). See README "Building the distillation subset" + run_aws.py for details.
set -euo pipefail

# ============================ CONFIG (edit these) ============================
BUCKET="marinfold-exp91-usw2"             # S3 bucket you own in us-west-2 (name must start "marinfold")
REGION="us-west-2"
IAM_PROFILE="marinfold-exp91-instance-profile"   # EC2 instance profile wrapping role marinfold-exp91-instance-role
INSTANCE_TYPE="r7i.24xlarge"              # 96 vCPU / 768 GB: novelty is alignment-bound so 3x cores ~3x faster
                                          #   (~10h vs ~30h); 768 GB holds the parallel scan's ~400 GB at 64 workers.
SMOKE_LIMIT="2000000"                     # rows scanned in the smoke run
VOLUME_GB="2000"                          # root EBS: full-scale .m8 + linclust scratch (~1-2 TB)
SMOKE_VOLUME_GB="300"                     # smoke doesn't need the big scratch volume

# --- optional SSH access (both empty = no SSH; the .live log heartbeat is the
#     default way to watch a run) -----------------------------------------------
# Requires an EC2 key pair + a security group opening :22, which in turn need IAM
# perms the current marinfold-exp91 policy does NOT grant (ec2:CreateKeyPair or
# ImportKeyPair, ec2:CreateSecurityGroup, ec2:AuthorizeSecurityGroupIngress). Ask
# the account admin to grant those or to provision a key pair + SG and hand you
# the names, then fill these in.
KEY_NAME=""                               # EC2 key pair name (SSH login)
SECURITY_GROUP_ID=""                      # sg-... allowing inbound :22 from your IP

# --- AFDB novelty reference (afdb-24M struct-cluster reps, via exp41) ----------
# The novelty stage drops Atlas proteins already represented in our training set
# (afdb-24M). The reference is the 1.33M struct-cluster representatives exp41
# isolated; their amino-acid sequences already live in that experiment's MMseqs/
# Foldseek target DB (dbtype 0), so prep-ref extracts them with `convert2fasta`.
# We deliberately do NOT pull `timodonnell/afdb-24M` (1.2 TB, and it has no
# sequence column — only mmCIF text) just to recover sequences.
EXP41_REP_DB="../exp41_evals_foldseek_train_similarity/db_full/db/targetDB"

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

# --- mmseqs scaling (the novelty search is the heavy stage, not the scan) ------
# The first full run wedged because a single easy-search over ~150M survivors blew
# past RAM+disk. pipeline.py now searches the query set in chunks of this size and
# deletes each chunk's scratch immediately. SPLIT_MEMORY_LIMIT is an extra RAM cap.
QUERY_CHUNK_SEQS="5000000"                # survivor seqs per novelty/leakage chunk
SPLIT_MEMORY_LIMIT="180G"                 # cap mmseqs RAM (leaves ~76 GB on the 256 GB box for OS + id-sets)
SEARCH_SENSITIVITY="4.0"                  # mmseqs -s (5.7 default); 4.0 = faster, may miss a few near-40%-id hits
SCAN_WORKERS="64"                         # parallel scan procs (< 96 vCPU; ~400 GB of the 768 GB box)

# --- mid-scale test run (validate scaling before a full run) ------------------
TEST_LIMIT="40000000"                     # ~40M rows (~20x the smoke, ~3M survivors)
TEST_INSTANCE_TYPE="r7i.8xlarge"          # same box class as the full run; test just scans fewer rows
TEST_VOLUME_GB="600"

# --- local-smoke knobs --------------------------------------------------------
LOCAL_SMOKE_LIMIT="50000"                 # Atlas rows scanned in the local smoke
LOCAL_WORK="${LOCAL_WORK:-_localsmoke}"   # scratch dir (gitignored); wiped each run
EVAL_SEQS_CSV="../exp65_evals_low_msa_depth_proteins/data/candidate_sequences.csv"
# ============================================================================

# Derived locations
STAGING="s3://${BUCKET}/exp91/staging"
OUT="s3://${BUCKET}/exp91/out"
SMOKE_OUT="s3://${BUCKET}/exp91/smoke"
TEST_OUT="s3://${BUCKET}/exp91/test"
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
  # Extract the afdb-24M struct-cluster rep sequences (exp41's target DB, dbtype 0)
  # to FASTA and upload as the novelty reference. Run once; ~1.33M seqs, ~350 MB.
  command -v mmseqs >/dev/null || { echo "[prep-ref] mmseqs not on PATH" >&2; exit 1; }
  [ -e "$EXP41_REP_DB" ] || {
    echo "[prep-ref] rep DB not found: ${EXP41_REP_DB}" >&2
    echo "[prep-ref] sync exp41's db_full first (see exp41 README)" >&2
    exit 1; }
  mkdir -p "$AFDB_LOCAL_DIR"
  echo "[prep-ref] mmseqs convert2fasta ${EXP41_REP_DB} -> ${AFDB_FASTA}"
  mmseqs convert2fasta "$EXP41_REP_DB" "$AFDB_FASTA"
  echo "[prep-ref] extracted $(grep -c '^>' "$AFDB_FASTA") sequences"

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
    --volume-size-gb "$VOLUME_GB" \
    ${KEY_NAME:+--key-name "$KEY_NAME"} \
    ${SECURITY_GROUP_ID:+--security-group-id "$SECURITY_GROUP_ID"} \
    --min-plddt "$MIN_PLDDT" --min-ptm "$MIN_PTM" \
    --max-afdb-seq-id "$MAX_AFDB_SEQ_ID" --eval-max-seq-id "$EVAL_MAX_SEQ_ID" \
    --cluster-id "$CLUSTER_ID" --reps-per-cluster "$REPS_PER_CLUSTER" \
    --query-chunk-seqs "$QUERY_CHUNK_SEQS" \
    --scan-workers "$SCAN_WORKERS" \
    ${SPLIT_MEMORY_LIMIT:+--split-memory-limit "$SPLIT_MEMORY_LIMIT"} \
    ${SEARCH_SENSITIVITY:+--search-sensitivity "$SEARCH_SENSITIVITY"} \
    --compute-plddt-std --watch "$@"
}

smoke() {
  echo "[smoke] capped run (--limit ${SMOKE_LIMIT}) -> ${SMOKE_OUT}"
  # Trailing --volume-size-gb overrides the full-run VOLUME_GB (argparse: last wins).
  launch "$SMOKE_OUT" --limit "$SMOKE_LIMIT" --volume-size-gb "$SMOKE_VOLUME_GB"
  echo "[smoke] check: aws s3 ls ${SMOKE_OUT}/"
}

test-run() {
  # Mid-scale validation before committing to the full 1.1B run: ~40M rows on a
  # smaller box, so the chunked novelty search + logging get exercised at a scale
  # where the old code would already have thrashed. Watch pipeline.log.live.
  echo "[test] ${TEST_LIMIT} rows on ${TEST_INSTANCE_TYPE} -> ${TEST_OUT}"
  # Trailing overrides win (argparse: last wins) over the full-run defaults.
  launch "$TEST_OUT" --limit "$TEST_LIMIT" \
    --instance-type "$TEST_INSTANCE_TYPE" --volume-size-gb "$TEST_VOLUME_GB"
  echo "[test] live log: aws s3 cp ${TEST_OUT}/pipeline.log.live -"
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
  test)        test-run ;;
  run)         run ;;
  manifest)    manifest ;;
  all)         prep_ref; run; manifest ;;
  *) echo "usage: $0 {local-smoke|prep-ref|smoke|test|run|manifest|all}" >&2; exit 2 ;;
esac
