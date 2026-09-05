# Runbook: fetching the val-eligible dimers on AWS

Phase 2 of [`CURATION_PLAN.md`](CURATION_PLAN.md). Fetches **881,010 val-eligible
complexes** (103 GB gzipped) from EBI into tar shards on S3.

No sequence collapse. PINDER clusters every system in its corpus rather than a
sequence-collapsed subset, and following PINDER is the decision recorded for this
dataset, so we fetch the whole val side. The collapse experiment and its numbers
stay in `val_cluster_summary.json` as a measured alternative, not as the path
taken.

## Shape of the job

| | |
|---|---|
| complexes | 881,010, from `dimer_split_assignment_id30_cov50.csv.gz` where `split == val` |
| transfer | 103 GB gzipped (117.0 kB/model measured over 8,012 models) |
| objects | 441 tars of 2,000 models each, ~234 MB per object, plus a sidecar JSON each |
| wall time | ~21 h at 11.5 files/s |
| cost | ~$2 EC2 + $0 ingress + ~$0.005 PUT + $2.37/mo storage |

**One instance, eight workers.** The pilot showed throughput linear in workers to
8 with no throttling, because `alphafold.ebi.ac.uk/files/` is GCS-backed. We do
not fan out across hosts: EBI publishes no rate limit, 11.5 files/s is the rate we
have demonstrated is safe, and 21 hours on a $0.09/h instance is not worth
multiplying our footprint against a public resource. The `--shard` / `--n-shards`
flags exist for the full-corpus case and are left at their defaults here.

## Launch

- **Instance**: `c7i.large` (2 vCPU, 4 GiB), on demand. The job is network-bound,
  not CPU-bound. One shard is assembled in memory before it is written, which is
  ~234 MB, so 4 GiB is ample.
- **Disk**: 20 GB root is enough. Structures stream to S3; nothing large is staged
  locally.
- **Region**: `us-west-2`, matching the account's existing
  `marinfold-exp91-usw2` bucket. Ingress from EBI is free from anywhere, so the
  region to pick is the one the data will be *read* from later, and exp91's data
  is already there.
- **IAM**: attach an instance role with `s3:PutObject`, `s3:GetObject` and
  `s3:ListBucket` on the destination prefix. **Do not put credentials in the
  environment.** exp91's short-lived CloudShell STS creds expire in ~15 minutes,
  which will not survive a 21-hour run; the instance role auto-refreshes and is
  the only thing that works here.

**Never commit the credential file.** `aws_keys.env` is gitignored at the repo
root; a `git add -A` in an experiment directory will otherwise sweep it into a
commit.

## Bootstrap

Do **not** `pip install boto3`. It reliably wedges these boxes about 60 seconds in
and the instance self-terminates; exp91 lost two days to this. `pyarrow` installs
cleanly and picks up the instance role on its own.

```bash
sudo apt-get update -qq && sudo apt-get install -y -qq python3-pip git
pip3 install --quiet pyarrow
git clone --depth 1 --branch exp145/multimer-data-survey \
    https://github.com/Open-Athena/MarinFold.git
cd MarinFold/experiments/exp145_data_multimer_data_survey
```

## Step 1: smoke test one shard

Run this on the instance role before committing to 21 hours. The procedure was
validated from a workstation on 2026-09-05 against `s3://marinfold-exp91-usw2/`:
23 models written as a 2.7 MB tar plus a 330-byte sidecar, read back to valid
mmCIF, resume correctly skipping the existing shard. What it has *not* been run
against is an EC2 instance role, which is the only credential source that
survives 21 hours.

```bash
python3 curate_fetch_structures.py --split val \
    --ids data/curation/dimer_split_assignment_id30_cov50.csv.gz \
    --out s3://marinfold-exp91-usw2/MarinFold/exp145-afdb-dimers/val/ \
    --shard 0 --n-shards 40000 --shard-size 25 --workers 8
```

Expect `n_ok: 23, n_missing_404: 0, n_failed: 0` and two objects at the
destination. Then confirm a member is real mmCIF:

```bash
python3 -c "
import io, tarfile, gzip, pyarrow.fs as pafs
fs = pafs.S3FileSystem(region=pafs.resolve_s3_region('marinfold-exp91-usw2'))
raw = fs.open_input_stream('marinfold-exp91-usw2/MarinFold/exp145-afdb-dimers/val/shard_000_00000.tar').read()
t = tarfile.open(fileobj=io.BytesIO(raw))
name = t.getnames()[0]
text = gzip.decompress(t.extractfile(name).read()).decode()
print(name, text.splitlines()[0])
# label_asym_id is field 6 and both chains must appear. Field 7 is
# label_entity_id, legitimately '1' for both chains of a homodimer, so checking
# that one makes a correct two-chain file look like it has only one.
print('chains:', sorted({l.split()[6] for l in text.splitlines() if l.startswith('ATOM')}))
"
```

Delete those two test objects before step 2, or the real run will skip
`shard_000_00000.tar` and silently omit 2,000 models.

## Step 2: the full run

```bash
tmux new -s fetch
python3 curate_fetch_structures.py --split val \
    --ids data/curation/dimer_split_assignment_id30_cov50.csv.gz \
    --out s3://marinfold-exp91-usw2/MarinFold/exp145-afdb-dimers/val/ \
    --workers 8 2>&1 | tee fetch.log
```

`tmux` because 21 hours will outlive any SSH session.

## Resume

Interrupt-safe at shard granularity. Rerunning the same command lists the
destination once, skips every tar already present, and continues. At most one
shard's work (2,000 models, ~3 minutes) is redone. A shard with a non-404 failure
after all five retries raises and is **not** written, so a partial tar can never
be mistaken for a complete one.

## What to check when it finishes

- **441 tars and 441 sidecar JSONs.** Fewer means a shard is missing, not that
  the run was short.
- **The 404 manifest.** Aggregate `missing_404` across the sidecars. The pilot saw
  28 in 2,000 in one ID band and 0 in 400 sampled across the whole confident set,
  so expect a small clustered count. This list ships with the dataset: the
  dataset is defined by what actually downloaded.
- **`fetch_timings.csv`** stays on the instance and is not written to S3. Copy it
  off before terminating; `AGENTS.md` requires per-input timings for every
  predictor and fetch run.
