# exp91 distillation pipeline — how to use it

This is the operator guide for turning the ESM Atlas (1.1B predicted monomers,
`s3://esm-protein-atlas/`, AWS us-west-2) into a quality-filtered, novel,
de-duplicated parquet that feeds `contacts-v1 generate`. For *why* (the
characterization, the funnel design, the distillation survey), read
[`README.md`](README.md); this file is just the runbook.

## The pieces

| File | Role |
|---|---|
| `atlas_io.py` | anonymous Lance reads + `structure_blob` (brotli→msgpack atom37) decode |
| `funnel.py` | sample-tuning of the metadata gates + the pure `apply_funnel()` |
| `pipeline.py` | on-instance worker: scan → novelty → leakage → cluster → select → `selected_manifest.csv` |
| `selection.py` | one rep per cluster (ESMFold2 rule: longest, then lowest pLDDT std-dev) |
| `materialize.py` | manifest → decode each Atlas rep's structure → parquet in the schema `contacts-v1 generate` reads (`entry_id` + `cif_content`) |
| `run_aws.py` | boto3 launcher: stages the scripts, runs `pipeline.py` on a throwaway EC2 box, self-terminates |
| `create_dataset.sh` | end-to-end driver wrapping all of the above |

The funnel reads only metadata + sequence (cheap, one pass). `materialize.py` is
the only stage that decodes `structure_blob`, so it is the heavy, multi-TB step.

## 1. Sanity-check locally (no AWS)

Runs the *entire* chain (funnel + materialize) on a tiny Atlas slice, asserts the
output parquet's schema and that every CIF re-parses. Needs `mmseqs` on PATH; reads
the Atlas anonymously over the network.

```bash
./create_dataset.sh local-smoke
```

Use this after any change to a pipeline script. It is a correctness check, not a
size estimate (it scans 50k rows, not 1.1B).

## 2. Production run (AWS us-west-2)

Edit the `CONFIG` block at the top of `create_dataset.sh` first (see below), then:

```bash
./create_dataset.sh prep-ref    # one-time: pull afdb-24M seqs from the HF bucket
                                #   -> FASTA -> s3://$BUCKET/refs/afdb_ref.fasta
./create_dataset.sh smoke       # capped EC2 run (~$1) to validate the cloud path
./create_dataset.sh run         # full funnel on EC2 -> selected_manifest.csv in S3
./create_dataset.sh manifest    # download the manifest locally
```

Then materialize the chosen reps to structures (decoding the Atlas `structure_blob`;
runs on an in-region instance):

```bash
python materialize.py --manifest selected_manifest.csv --out-dir dataset/
# shardable across instances: --num-shards N --shard-id i
```

Publishing `dataset/` to HuggingFace is a **separate, deliberate** step (cross-cloud,
>10 GB → needs sign-off per the repo rules, and the Atlas is CC BY-SA). The pipeline
never pushes on its own.

## Config you must set (`create_dataset.sh`)

- `BUCKET` — an S3 bucket you own in `us-west-2`.
- `IAM_PROFILE` — an instance profile with S3 read/write on that bucket.
- `HF_AFDB_PREFIX`, `AFDB_ID_COL`, `AFDB_SEQ_COL` — where the afdb-24M sequences live
  in `hf://buckets/open-athena/MarinFold` and their parquet column names
  (`hf buckets ls hf://buckets/open-athena/MarinFold` to find them). These are
  placeholders until confirmed.

## The two knobs that set quality and size

- **`MIN_PLDDT` (default 0.70) is a quality floor, not the size dial.** It is the
  paper's "high confidence" line. Keep it loose so the sole member of a novel
  structural cluster is not deleted just because its single prediction landed at 0.78.
- **`CLUSTER_ID` (default 0.40) sets the final size.** Clustering + one best rep per
  cluster is what controls how many structures you publish (and guarantees
  diversity). Tighten `CLUSTER_ID` (coarser clusters) to shrink the set. Calibrate it
  on the `smoke` run; the cluster count is not predictable a priori.

## Running it in a locked-down AWS account

`run_aws.py` needs `ec2:RunInstances`, `iam:PassRole` (on the instance-profile role),
`s3:*` on `$BUCKET`, and read EC2/SSM for the AMI lookup. With those granted, no
static access keys are required — AWS CloudShell runs as your identity and the
launcher picks up its credentials automatically. The full least-privilege IAM policy
and the resources an account admin must provision are in the PR description for this
branch.
