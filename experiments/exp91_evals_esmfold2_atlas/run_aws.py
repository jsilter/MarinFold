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
           "materialize.py"]

# Cloud-init bootstrap. {placeholders} are filled by render_user_data().
USER_DATA_TMPL = r"""#!/bin/bash
set -euxo pipefail
exec > /var/log/marinfold-pipeline.log 2>&1
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

# Run the funnel. Failure still uploads the log + a _FAILED marker for debugging.
set +e
python3 pipeline.py --work-dir /opt/work --stage all \
    --afdb-ref afdb_ref.fasta --eval-ref eval_seqs.fasta \
    {pipeline_args}
RC=$?
set -e

aws s3 cp --recursive /opt/work {output}/ || true
aws s3 cp /var/log/marinfold-pipeline.log {output}/pipeline.log || true
if [[ $RC -eq 0 ]]; then echo ok | aws s3 cp - {output}/_DONE; else echo "rc=$RC" | aws s3 cp - {output}/_FAILED; fi
shutdown -h now
"""


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
    run_kw = dict(
        ImageId=ami, InstanceType=args.instance_type, MinCount=1, MaxCount=1,
        UserData=render_user_data(args),
        IamInstanceProfile={"Name": args.iam_instance_profile},
        BlockDeviceMappings=bdm,
        InstanceInitiatedShutdownBehavior="terminate",
        TagSpecifications=[{"ResourceType": "instance",
                            "Tags": [{"Key": "Name", "Value": "marinfold-exp91-funnel"},
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
    # AWS wiring
    ap.add_argument("--s3-staging", required=True, help="s3://bucket/prefix for scripts")
    ap.add_argument("--s3-output", required=True, help="s3://bucket/prefix for results")
    ap.add_argument("--afdb-ref-uri", required=True,
                    help="s3:// or https:// FASTA of your AFDB training sequences")
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
    ap.add_argument("--dry-run", action="store_true",
                    help="print user-data + config; touch nothing")
    args = ap.parse_args(argv)

    if args.dry_run:
        print("=== pipeline args ===")
        print(build_pipeline_args(args))
        print("\n=== user-data (cloud-init) ===")
        print(render_user_data(args))
        print("=== run config ===")
        for k in ("region", "instance_type", "volume_size_gb", "s3_staging",
                  "s3_output", "afdb_ref_uri", "iam_instance_profile", "spot"):
            print(f"  {k}: {getattr(args, k)}")
        return

    launch(args)


if __name__ == "__main__":
    main()
