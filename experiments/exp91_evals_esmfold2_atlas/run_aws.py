# Copyright The MarinFold Authors
# SPDX-License-Identifier: Apache-2.0

"""Launch the ESM Atlas distillation funnel on a throwaway EC2 instance.

The Atlas lives in AWS ``us-west-2``; ``pipeline.py`` reads it *in region* (no
egress) and is a CPU MMseqs + metadata job. This launcher provisions one
large-memory instance, hands it a cloud-init bootstrap that installs MMseqs +
the Python deps, downloads the pipeline scripts and the AFDB reference, runs the
funnel, uploads the result manifest to S3, and **terminates itself**. You bring
your own AWS account; this never runs implicitly.

Flow:
  1. tar the pipeline scripts (+ a generated eval FASTA) and upload to S3 staging.
  2. RunInstances in us-west-2 with a user-data bootstrap (terminate-on-shutdown).
  3. The instance runs ``pipeline.py --stage all`` and writes
     ``<s3-output>/selected_manifest.csv`` + ``pipeline.log`` + a ``_DONE`` marker.
  4. (optional) ``--watch`` polls S3 for the ``_DONE`` marker.

Prerequisites (all yours to provide):
  * AWS credentials in the environment (``aws configure`` / env vars / SSO).
  * An **S3 bucket you own**, in us-west-2, for staging + output
    (``--s3-staging``, ``--s3-output``).
  * An **IAM instance profile** (``--iam-instance-profile``) whose role can
    read/write that bucket (the instance uses it to fetch scripts + push results).
  * An **AFDB / training-set reference FASTA** (``--afdb-ref-uri``, ``s3://`` or
    ``https://``) — your ~10M afdb-24M training sequences — for the novelty stage.
  * Optionally a key pair (``--key-name``) + security group for SSH debugging.

``--dry-run`` renders the user-data and run config without touching AWS (and
without needing boto3 installed).
"""

import argparse
import io
import sys
import tarfile
import textwrap
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
EXP65_SEQS = (HERE.parent / "exp65_evals_low_msa_depth_proteins"
              / "data" / "candidate_sequences.csv")
# Scripts the instance needs (pipeline imports atlas_io + selection; funnel is
# carried for parity / post-selection QC).
SCRIPTS = ["pipeline.py", "atlas_io.py", "selection.py", "funnel.py",
           "materialize.py", "probe_throughput.py"]

# Cloud-init bootstrap. {placeholders} are filled by render_user_data().
USER_DATA_TMPL = r"""#!/bin/bash
exec > /var/log/marinfold-pipeline.log 2>&1

# ------------------------------- fire-and-forget safety -------------------------
# This box is meant to be launched and forgotten, so termination must be
# guaranteed on EVERY exit path — success, pipeline failure, OR a crash during
# bootstrap (apt/pip/mmseqs) that would otherwise leave it running idle forever.
#
#  1. Hard 48h cap via a scheduled shutdown, set before anything can fail. Even if
#     the whole script wedges, the instance powers off (and with
#     InstanceInitiatedShutdownBehavior=terminate, terminates -> billing stops,
#     the DeleteOnTermination volume is freed). 48h is well above the ~12-30h
#     expected runtime; lower it if you want a tighter ceiling.
#  2. An EXIT trap that uploads whatever /opt/work has plus the log, drops a
#     _DONE / _FAILED marker, and powers off. `shutdown` and coreutils always
#     exist, so this fires even if apt never installed awscli (the S3 copies just
#     no-op in that case, but the box still terminates).
shutdown -h +2880 "marinfold-exp91 48h safety cap" || true
# NB: the function braces are doubled ({{ }}) so str.format() leaves them intact
# and only fills the {output} placeholders.
finish() {{
  rc=$?
  trap - EXIT INT TERM   # de-register so finish runs exactly once
  # Save the small, critical artifacts FIRST. On a SIGTERM (manual terminate or
  # the 48h cap) systemd grants only ~90s before SIGKILL — not enough for the full
  # /opt/work upload — so prioritise the log + status marker, then best-effort the
  # bulk work dir. On a normal exit there's no time limit and all of it uploads.
  aws s3 cp /var/log/marinfold-pipeline.log {output}/pipeline.log 2>/dev/null || true
  if [ -f /opt/work/_pipeline_ok ]; then echo ok | aws s3 cp - {output}/_DONE 2>/dev/null || true
  else echo "rc=$rc" | aws s3 cp - {output}/_FAILED 2>/dev/null || true; fi
  # Skip mmseqs scratch dirs (_tmp_*, _chunks_*) — worthless, and can be multi-TB.
  aws s3 cp --recursive /opt/work {output}/ \
    --exclude "_tmp*" --exclude "_chunks*" 2>/dev/null || true
  shutdown -c 2>/dev/null || true   # cancel the 48h cap, then go now
  shutdown -h now
}}
# Catch signals too (EXIT alone does not fire on a SIGTERM from terminate/timeout),
# so a killed or capped run still uploads its log + partial results.
trap finish EXIT INT TERM
# --------------------------------------------------------------------------------

set -euxo pipefail
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y python3-pip awscli wget tar

# MMseqs2 static binary (AVX2 build).
cd /opt
wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
tar xzf mmseqs-linux-avx2.tar.gz
export PATH=/opt/mmseqs/bin:$PATH

# Python deps for pipeline.py / atlas_io.py.
pip3 install --quiet pylance pyarrow pandas numpy gemmi msgpack msgpack-numpy \
    brotli zstandard biopython

# Pipeline scripts + generated eval FASTA.
mkdir -p /opt/exp91 && cd /opt/exp91
aws s3 cp {staging}/scripts.tar.gz .
tar xzf scripts.tar.gz

# AFDB / training reference for the novelty stage.
REF={afdb_ref_uri}
if [[ "$REF" == s3://* ]]; then aws s3 cp "$REF" afdb_ref.fasta; else wget -q -O afdb_ref.fasta "$REF"; fi

# Heartbeat: publish the live log to S3 every 2 min so the run is observable
# mid-flight (there is otherwise no window in — no SSH/SSM). Detached so it can't
# block boot; the shutdown in the exit trap kills it. Writes to a distinct
# .live key so it never races the final pipeline.log the trap uploads.
( while true; do
    aws s3 cp /var/log/marinfold-pipeline.log {output}/pipeline.log.live 2>/dev/null || true
    aws s3 sync /opt/work {output}/ --exclude "*" --include "_progress_*/*" 2>/dev/null || true
    sleep 120
  done ) &
disown || true

# Resume: pull any prior novelty/leakage chunk checkpoints so a relaunch continues
# from the last finished chunk instead of restarting the multi-hour search. No-op
# on a fresh run. The ~37-min scan re-runs deterministically, regenerating the same
# survivors -> same chunk splits, so an earlier chunk's drop file stays valid.
mkdir -p /opt/work
aws s3 sync {output}/ /opt/work/ --exclude "*" --include "_progress_*/*" 2>/dev/null || true

# Run the funnel. On success, drop the marker the EXIT trap checks for; on any
# failure the trap still uploads the partial work + log + a _FAILED marker.
python3 pipeline.py --work-dir /opt/work --stage all \
    --afdb-ref afdb_ref.fasta --eval-ref eval_seqs.fasta \
    {pipeline_args}
touch /opt/work/_pipeline_ok
"""


# Materialize (stage 3): decode the selected reps' structures. Same fire-and-forget
# safety as the funnel, but (a) output is ~3.2 TB of parquet parts streamed to S3 as
# they are written (each part IS a checkpoint), and (b) the finish trap uses
# `aws s3 sync` — never `cp --recursive` — so it does not re-upload terabytes the
# heartbeat already pushed. Resume: a relaunch pulls the (small) plan/ back and
# builds a skip-list of already-uploaded parts, so decode continues where it stopped.
MATERIALIZE_USER_DATA_TMPL = r"""#!/bin/bash
exec > /var/log/marinfold-pipeline.log 2>&1
shutdown -h +2880 "marinfold-exp91 48h safety cap" || true
finish() {{
  rc=$?
  trap - EXIT INT TERM
  aws s3 cp /var/log/marinfold-pipeline.log {output}/pipeline{marker_suffix}.log 2>/dev/null || true
  if [ -f /opt/work/_pipeline_ok ]; then echo ok | aws s3 cp - {output}/_DONE{marker_suffix} 2>/dev/null || true
  else echo "rc=$rc" | aws s3 cp - {output}/_FAILED{marker_suffix} 2>/dev/null || true; fi
  # sync (not cp --recursive): parts the heartbeat already pushed are skipped, so
  # this only flushes stragglers instead of re-uploading the whole multi-TB output.
  # Exclude in-progress .tmp: only completed part_*.parquet are real checkpoints.
  aws s3 sync /opt/work/parts {output}/parts --exclude "*.tmp" 2>/dev/null || true
  shutdown -c 2>/dev/null || true
  shutdown -h now
}}
trap finish EXIT INT TERM

set -euxo pipefail
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y python3-pip awscli
pip3 install --quiet pylance pyarrow pandas numpy gemmi msgpack msgpack-numpy \
    brotli zstandard

mkdir -p /opt/exp91 && cd /opt/exp91
aws s3 cp {staging}/scripts.tar.gz .
tar xzf scripts.tar.gz
aws s3 cp {manifest_uri} selected_manifest.csv

# Heartbeat: stream the live log + finished parts (the checkpoints) to S3 every 2 min,
# so a crash loses at most the in-flight chunks and the run is observable mid-flight.
# ONLY completed part_*.parquet are synced (--exclude "*.tmp"): the atomic-write temp
# files are large, grow for the whole chunk, and re-uploading them every cycle wasted
# bandwidth and starved the decode workers. The plan already lives in S3 (pulled below),
# so it is never synced back up.
mkdir -p /opt/work/parts /opt/work/plan
( while true; do
    aws s3 cp /var/log/marinfold-pipeline.log {output}/pipeline{marker_suffix}.log.live 2>/dev/null || true
    aws s3 sync /opt/work/parts {output}/parts --exclude "*.tmp" 2>/dev/null || true
    sleep 120
  done ) &
disown || true

# Resume: pull the (small) plan back and list already-COMPLETED parts as chunk ids to
# skip. The `$` anchor matches only final part_NNNNN.parquet, never a partial
# .part_NNNNN.parquet.tmp (which contains the same substring), so a partial is never
# mistaken for done. The plan is deterministic, so on a fresh run these are simply empty.
aws s3 sync {output}/plan /opt/work/plan 2>/dev/null || true
aws s3 ls {output}/parts/ 2>/dev/null | grep -oE 'part_[0-9]+\.parquet$' \
    | grep -oE '[0-9]+' > /opt/work/skip.txt || true

python3 materialize.py full --manifest selected_manifest.csv --work-dir /opt/work \
    --skip-list /opt/work/skip.txt {materialize_args}
touch /opt/work/_pipeline_ok
"""


def build_materialize_args(args: argparse.Namespace) -> str:
    """Render the materialize.py passthrough flags from the launcher args."""
    parts = ["--chunk-size", args.chunk_size,
             "--scan-workers", args.scan_workers,
             "--workers", args.decode_workers]
    if args.take_batch is not None:
        parts += ["--take-batch", args.take_batch]
    if args.mat_num_shards > 1:
        parts += ["--num-shards", args.mat_num_shards, "--shard-id", args.mat_shard_id]
    if args.limit is not None:
        parts += ["--limit", args.limit]
    return " ".join(str(p) for p in parts)


def render_materialize_user_data(args: argparse.Namespace) -> str:
    # Per-shard done/failed markers so a sharded run's completion is unambiguous: with
    # N boxes all writing one parts/ prefix, a single _DONE would be written N times and
    # can't tell you all shards finished. Single-box (num_shards=1) keeps the bare _DONE.
    marker_suffix = "" if args.mat_num_shards <= 1 else f"_shard_{args.mat_shard_id}"
    return MATERIALIZE_USER_DATA_TMPL.format(
        staging=args.s3_staging.rstrip("/"),
        output=args.s3_output.rstrip("/"),
        manifest_uri=args.manifest_uri,
        marker_suffix=marker_suffix,
        materialize_args=build_materialize_args(args))


# Probe (throughput sweep): pull the plan from S3, sweep (workers x take_batch) on the
# real target box for a few minutes each, ship the results log to S3, self-terminate.
# No manifest / no output parts — this only measures, it does not materialize.
PROBE_USER_DATA_TMPL = r"""#!/bin/bash
exec > /var/log/marinfold-pipeline.log 2>&1
shutdown -h +60 "marinfold-exp91 probe 1h cap" || true
finish() {{
  rc=$?
  trap - EXIT INT TERM
  aws s3 cp /var/log/marinfold-pipeline.log {output}/probe.log 2>/dev/null || true
  shutdown -c 2>/dev/null || true
  shutdown -h now
}}
trap finish EXIT INT TERM

set -euxo pipefail
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y python3-pip awscli
pip3 install --quiet pylance pyarrow pandas numpy gemmi msgpack msgpack-numpy \
    brotli zstandard

mkdir -p /opt/exp91 && cd /opt/exp91
aws s3 cp {staging}/scripts.tar.gz .
tar xzf scripts.tar.gz

# Heartbeat the live probe log to S3 so we can read the sweep as it runs.
( while true; do
    aws s3 cp /var/log/marinfold-pipeline.log {output}/probe.log.live 2>/dev/null || true
    sleep 30
  done ) &
disown || true

mkdir -p /opt/work/plan
aws s3 sync {output}/plan /opt/work/plan 2>/dev/null || true
python3 probe_throughput.py --work-dir /opt/work {probe_args}
"""


def build_probe_args(args: argparse.Namespace) -> str:
    parts = ["--seconds", args.probe_seconds,
             "--workers", *args.probe_workers,
             "--take-batch", *args.probe_take_batch]
    if args.probe_blobs_only:
        parts.append("--blobs-only")
    return " ".join(str(p) for p in parts)


def render_probe_user_data(args: argparse.Namespace) -> str:
    return PROBE_USER_DATA_TMPL.format(
        staging=args.s3_staging.rstrip("/"),
        output=args.s3_output.rstrip("/"),
        probe_args=build_probe_args(args))


def build_pipeline_args(args: argparse.Namespace) -> str:
    """Render the pipeline.py passthrough flags from the launcher args."""
    parts = [
        "--min-plddt", args.min_plddt, "--min-ptm", args.min_ptm,
        "--max-afdb-seq-id", args.max_afdb_seq_id,
        "--eval-max-seq-id", args.eval_max_seq_id,
        "--cluster-id", args.cluster_id, "--reps-per-cluster", args.reps_per_cluster,
        "--num-shards", args.num_shards,
    ]
    if args.compute_plddt_std:
        parts.append("--compute-plddt-std")
    if args.limit is not None:
        parts += ["--limit", args.limit]
    parts += ["--query-chunk-seqs", args.query_chunk_seqs,
              "--scan-workers", args.scan_workers]
    if args.split_memory_limit:
        parts += ["--split-memory-limit", args.split_memory_limit]
    if args.search_sensitivity:
        parts += ["--search-sensitivity", args.search_sensitivity]
    return " ".join(str(p) for p in parts)


def render_user_data(args: argparse.Namespace) -> str:
    return USER_DATA_TMPL.format(
        staging=args.s3_staging.rstrip("/"),
        output=args.s3_output.rstrip("/"),
        afdb_ref_uri=args.afdb_ref_uri,
        pipeline_args=build_pipeline_args(args))


def build_scripts_tarball() -> bytes:
    """Tar the pipeline scripts + a freshly generated eval FASTA (in memory)."""
    df = pd.read_csv(EXP65_SEQS)
    eval_fasta = "".join(f">{r.stem}\n{r.sequence}\n" for r in df.itertuples())
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name in SCRIPTS:
            tar.add(HERE / name, arcname=name)
        data = eval_fasta.encode()
        info = tarfile.TarInfo("eval_seqs.fasta")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    print(f"[pack] scripts.tar.gz: {len(SCRIPTS)} scripts + eval_seqs.fasta "
          f"({len(df)} seqs), {buf.tell() / 1e3:.0f} KB")
    return buf.getvalue()


def resolve_ami(ec2, region: str) -> str:
    """Latest Canonical Ubuntu 22.04 amd64 server AMI via ``DescribeImages``.

    We look the image up directly (Canonical owner ``099720109477``) rather than
    via the public SSM parameter path: SSM's ``get-parameter`` is not reachable in
    every account (locked-down IAM returns ``ParameterNotFound``), whereas
    ``ec2:DescribeImages`` is already required for the launch. Sort the available
    ``hvm-ssd*`` (gp2 or gp3) images by creation date and take the newest.
    """
    resp = ec2.describe_images(
        Owners=["099720109477"],
        Filters=[
            {"Name": "name",
             "Values": ["ubuntu/images/hvm-ssd*/ubuntu-jammy-22.04-amd64-server-*"]},
            {"Name": "state", "Values": ["available"]},
        ],
    )
    images = sorted(resp["Images"], key=lambda im: im["CreationDate"])
    if not images:
        raise SystemExit(
            f"no Canonical Ubuntu 22.04 amd64 AMI found in {region}; pass --ami")
    return images[-1]["ImageId"]


def launch(args: argparse.Namespace) -> None:
    import boto3
    s3 = boto3.client("s3", region_name=args.region)
    ec2 = boto3.client("ec2", region_name=args.region)

    # 1. stage scripts
    bucket, _, key_prefix = args.s3_staging.replace("s3://", "").partition("/")
    s3.put_object(Bucket=bucket, Key=f"{key_prefix.rstrip('/')}/scripts.tar.gz",
                  Body=build_scripts_tarball())

    ami = args.ami or resolve_ami(ec2, args.region)
    bdm = [{"DeviceName": "/dev/sda1",
            "Ebs": {"VolumeSize": args.volume_size_gb, "VolumeType": "gp3",
                    "DeleteOnTermination": True}}]
    if args.task == "materialize":
        user_data = render_materialize_user_data(args)
        name_tag = "marinfold-exp91-materialize"
    elif args.task == "probe":
        user_data = render_probe_user_data(args)
        name_tag = "marinfold-exp91-probe"
    else:
        user_data = render_user_data(args)
        name_tag = "marinfold-exp91-funnel"
    run_kw = dict(
        ImageId=ami, InstanceType=args.instance_type, MinCount=1, MaxCount=1,
        UserData=user_data,
        IamInstanceProfile={"Name": args.iam_instance_profile},
        BlockDeviceMappings=bdm,
        InstanceInitiatedShutdownBehavior="terminate",
        TagSpecifications=[{"ResourceType": "instance",
                            "Tags": [{"Key": "Name", "Value": name_tag},
                                     {"Key": "project", "Value": "MarinFold"}]}])
    if args.key_name:
        run_kw["KeyName"] = args.key_name
    if args.security_group_id:
        run_kw["SecurityGroupIds"] = [args.security_group_id]
    if args.subnet_id:
        run_kw["SubnetId"] = args.subnet_id
    if args.spot:
        run_kw["InstanceMarketOptions"] = {"MarketType": "spot"}

    iid = ec2.run_instances(**run_kw)["Instances"][0]["InstanceId"]
    print(f"[launch] {iid} ({args.instance_type}, ami {ami}) in {args.region}")
    print(f"[launch] results -> {args.s3_output}  (manifest: selected_manifest.csv)")
    print(f"[launch] watch:  aws s3 ls {args.s3_output.rstrip('/')}/")
    print(f"[launch] logs:   aws s3 cp {args.s3_output.rstrip('/')}/pipeline.log -")
    if args.watch:
        _watch(s3, args.s3_output, args.region)


def _watch(s3, s3_output: str, region: str, poll_s: int = 60) -> None:
    bucket, _, prefix = s3_output.replace("s3://", "").partition("/")
    prefix = prefix.rstrip("/")
    print(f"[watch] polling {s3_output} for _DONE/_FAILED every {poll_s}s ...")
    while True:
        listing = s3.list_objects_v2(Bucket=bucket, Prefix=prefix).get("Contents", [])
        names = {o["Key"].rsplit("/", 1)[-1] for o in listing}
        if "_DONE" in names:
            print("[watch] _DONE — pipeline finished. Pull selected_manifest.csv.")
            return
        if "_FAILED" in names:
            print("[watch] _FAILED — check pipeline.log in the output prefix.")
            return
        time.sleep(poll_s)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    # task: funnel (scan->select) or materialize (decode selected reps to parquet)
    ap.add_argument("--task", choices=["funnel", "materialize", "probe"],
                    default="funnel")
    # AWS wiring
    ap.add_argument("--s3-staging", required=True, help="s3://bucket/prefix for scripts")
    ap.add_argument("--s3-output", required=True, help="s3://bucket/prefix for results")
    ap.add_argument("--afdb-ref-uri", default=None,
                    help="funnel: s3:// or https:// FASTA of your AFDB training seqs")
    # materialize task
    ap.add_argument("--manifest-uri", default=None,
                    help="materialize: s3:// selected_manifest.csv from the funnel")
    ap.add_argument("--chunk-size", type=int, default=20_000,
                    help="materialize: reps per plan chunk / output part")
    ap.add_argument("--decode-workers", type=int, default=1,
                    help="materialize: parallel decode processes (set ~vCPU count)")
    ap.add_argument("--take-batch", type=int, default=None,
                    help="materialize: reps per Lance take/decode sub-batch")
    ap.add_argument("--mat-num-shards", type=int, default=1,
                    help="materialize: split chunks across this many boxes")
    ap.add_argument("--mat-shard-id", type=int, default=0,
                    help="materialize: this box's shard index")
    # probe task (throughput sweep)
    ap.add_argument("--probe-workers", type=int, nargs="+",
                    default=[64, 128, 192, 256], help="probe: worker counts to sweep")
    ap.add_argument("--probe-take-batch", type=int, nargs="+", default=[512],
                    help="probe: take/decode sub-batch sizes to sweep")
    ap.add_argument("--probe-seconds", type=int, default=90,
                    help="probe: seconds per (workers, take_batch) config")
    ap.add_argument("--probe-blobs-only", action="store_true",
                    help="probe: skip decode to isolate pure S3 read throughput")
    ap.add_argument("--iam-instance-profile", required=True,
                    help="instance profile name with S3 read/write on those buckets")
    ap.add_argument("--region", default="us-west-2")
    ap.add_argument("--instance-type", default="r7i.8xlarge",
                    help="big-RAM CPU instance; scale to survivor count")
    ap.add_argument("--ami", default=None, help="override; else latest Ubuntu 22.04")
    ap.add_argument("--volume-size-gb", type=int, default=2000,
                    help="root EBS gp3 size; the full-scale novelty .m8 + linclust "
                         "scratch needs ~1-2 TB (500 GB is only enough for smoke)")
    ap.add_argument("--key-name", default=None, help="EC2 key pair for SSH (optional)")
    ap.add_argument("--security-group-id", default=None)
    ap.add_argument("--subnet-id", default=None)
    ap.add_argument("--spot", action="store_true", help="request a spot instance")
    ap.add_argument("--watch", action="store_true", help="poll S3 until done")
    # pipeline passthrough (defaults mirror funnel.py)
    ap.add_argument("--min-plddt", type=float, default=0.70)
    ap.add_argument("--min-ptm", type=float, default=0.50)
    ap.add_argument("--max-afdb-seq-id", type=float, default=0.40)
    ap.add_argument("--eval-max-seq-id", type=float, default=0.40)
    ap.add_argument("--cluster-id", type=float, default=0.40)
    ap.add_argument("--reps-per-cluster", type=int, default=1)
    ap.add_argument("--compute-plddt-std", action="store_true")
    ap.add_argument("--num-shards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=None, help="cap rows (smoke test)")
    ap.add_argument("--query-chunk-seqs", type=int, default=5_000_000,
                    help="novelty/leakage query-chunk size (bounds mmseqs RAM+disk)")
    ap.add_argument("--split-memory-limit", default=None,
                    help="mmseqs --split-memory-limit (e.g. 100G); extra RAM cap")
    ap.add_argument("--search-sensitivity", default=None,
                    help="mmseqs -s for novelty/leakage (default 5.7)")
    ap.add_argument("--scan-workers", type=int, default=1,
                    help="parallel scan processes (set ~vCPU count)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print user-data + config; touch nothing")
    args = ap.parse_args(argv)

    # Per-task required inputs (kept out of argparse `required` so each task only
    # demands what it uses).
    if args.task == "funnel" and not args.afdb_ref_uri:
        ap.error("--afdb-ref-uri is required for --task funnel")
    if args.task == "materialize" and not args.manifest_uri:
        ap.error("--manifest-uri is required for --task materialize")

    if args.dry_run:
        if args.task == "probe":
            print("=== probe args ===")
            print(build_probe_args(args))
            print("\n=== user-data (cloud-init) ===")
            print(render_probe_user_data(args))
        elif args.task == "materialize":
            print("=== materialize args ===")
            print(build_materialize_args(args))
            print("\n=== user-data (cloud-init) ===")
            print(render_materialize_user_data(args))
        else:
            print("=== pipeline args ===")
            print(build_pipeline_args(args))
            print("\n=== user-data (cloud-init) ===")
            print(render_user_data(args))
        print("=== run config ===")
        for k in ("task", "region", "instance_type", "volume_size_gb", "s3_staging",
                  "s3_output", "afdb_ref_uri", "manifest_uri", "iam_instance_profile",
                  "spot"):
            print(f"  {k}: {getattr(args, k)}")
        return

    launch(args)


if __name__ == "__main__":
    main()
