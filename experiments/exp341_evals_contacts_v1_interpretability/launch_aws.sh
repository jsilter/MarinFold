#!/usr/bin/env bash
# Launch one self-terminating EC2 instance for an exp341 job.
#
# Jobs:
#   harness    GPU. Harness checks 1-2 at the plan's 4 seeds on CUDA
#              (check_harness.py) and check 5, the timed clean readout over
#              eval-val and eval-denovo (time_readout.py).
#   panel      GPU. A3/A4 readouts (score_panel.py): exp232 and exp199 on the
#              1,000 removed + 1,000 kept structures of data/panel.parquet.
#   relatives  CPU. Panel stage 2 (relative_counts.py): decode exp232's
#              training sequences from the HF bucket and MMseqs2-search the
#              stage 1 candidates (CANDIDATES=path/to/candidates.parquet).
#
# The instance copies everything under out/ plus its log to S3, writes a
# _DONE_rc<N> marker, and shuts down (shutdown behaviour is terminate). A
# `shutdown -h +MAX_MINUTES` scheduled first caps the cost if anything hangs.
#
# No secrets reach the instance: code (and any inputs) come from S3 through
# the instance role, and model checkpoints and corpora from the public
# MarinFold HF bucket, read anonymously.
#
# Account facts (no key pairs; the instance profile name is not discoverable
# with this user's IAM permissions) are from the exp91/exp145 EC2 workers.
#
# Usage (CloudShell STS credentials expire after ~15 minutes; REGION and
# SUBNET must name the same region):
#   set -a; source experiments/exp341_evals_contacts_v1_interpretability/aws_keys.env; set +a
#   experiments/exp341_evals_contacts_v1_interpretability/launch_aws.sh harness
#   INSTANCE_TYPE=g6e.xlarge experiments/exp341_evals_contacts_v1_interpretability/launch_aws.sh harness
#   CANDIDATES=/path/candidates.parquet experiments/exp341_evals_contacts_v1_interpretability/launch_aws.sh relatives
set -euo pipefail

JOB=${1:?usage: launch_aws.sh harness|panel|relatives}
# The bucket stays in us-west-2; harness inputs and outputs are tens of MB.
BUCKET=marinfold-exp91-usw2
BUCKET_REGION=us-west-2
INSTANCE_PROFILE=marinfold-exp91-instance-profile
EXP=experiments/exp341_evals_contacts_v1_interpretability
EXP245=experiments/exp245_evals_foldbench_held_out_monomers/data
# Default subnets: us-east-2a subnet-03f9b7abd4284ad5e, 2b subnet-02c630af81a9bdc13,
# 2c subnet-091e8c3572c921332; us-west-2a subnet-0a393c46ff06d81c3,
# 2b subnet-0563fac5b8bee3f0d. Capacity varies by zone, so SUBNET picks one.
case "$JOB" in
    harness)
        # GPU quota (2026-10): us-east-2 has 192 P / 128 G on-demand vCPUs;
        # us-west-2 has 4 P / 0 G, too few for any GPU type.
        INSTANCE_TYPE=${INSTANCE_TYPE:-p5.4xlarge}  # 1x H100 80GB
        REGION=${REGION:-us-east-2}
        SUBNET=${SUBNET:-subnet-03f9b7abd4284ad5e}
        MAX_MINUTES=${MAX_MINUTES:-120}
        DISK_GB=150
        AMI_PARAMETER=/aws/service/deeplearning/ami/x86_64/base-oss-nvidia-driver-gpu-ubuntu-22.04/latest/ami-id
        ;;
    panel)
        # ~1.4 L40S-hours per model, extrapolated from check 5's eval-val
        # timings to the panel's lengths (mean 430, max 1,551).
        INSTANCE_TYPE=${INSTANCE_TYPE:-g6e.xlarge}  # 1x L40S 48GB
        REGION=${REGION:-us-east-2}
        SUBNET=${SUBNET:-subnet-03f9b7abd4284ad5e}
        MAX_MINUTES=${MAX_MINUTES:-360}
        DISK_GB=150
        AMI_PARAMETER=/aws/service/deeplearning/ami/x86_64/base-oss-nvidia-driver-gpu-ubuntu-22.04/latest/ami-id
        ;;
    relatives)
        # CPU, next to the bucket. 256 GB RAM keeps the MMseqs2 prefilter from
        # splitting heavily; the disk holds ~40 GB of FASTA, the target DB and
        # search scratch.
        INSTANCE_TYPE=${INSTANCE_TYPE:-r7i.8xlarge}
        REGION=${REGION:-us-west-2}
        SUBNET=${SUBNET:-subnet-0a393c46ff06d81c3}
        MAX_MINUTES=${MAX_MINUTES:-720}
        DISK_GB=1000
        AMI_PARAMETER=/aws/service/canonical/ubuntu/server/22.04/stable/current/amd64/hvm/ebs-gp2/ami-id
        : "${CANDIDATES:?set CANDIDATES to the stage 1 candidates.parquet}"
        ;;
    *)
        echo "unknown job $JOB (harness|panel|relatives)" >&2
        exit 1
        ;;
esac

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"
if ! git diff --quiet HEAD -- marinfold "$EXP"; then
    echo "uncommitted changes under marinfold/ or $EXP; the instance runs HEAD, commit first" >&2
    exit 1
fi
sha=$(git rev-parse --short HEAD)
prefix="exp341/$JOB/$(date -u +%Y%m%dT%H%M%SZ)-$sha"
s3="s3://$BUCKET/$prefix"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

git archive --format=tar.gz -o "$work/code.tgz" HEAD \
    marinfold "$EXP" "$EXP245/eval_sets.csv" "$EXP245/gt_universe_scored.jsonl"
aws s3 cp --region "$BUCKET_REGION" "$work/code.tgz" "$s3/code.tgz"

case "$JOB" in
    harness)
        steps=$(cat <<'STEPS'
nvidia-smi
uv run python check_harness.py --device cuda --seeds 4 --output out/harness_checks_cuda_4seed.json
rc_a=$?
uv run python time_readout.py --device cuda --runner-tag aws --output-dir out/readout
rc_b=$?
STEPS
)
        ;;
    panel)
        steps=$(cat <<'STEPS'
nvidia-smi
uv run python score_panel.py --device cuda --runner-tag aws --output-dir out/panel
rc_a=$?
rc_b=0
STEPS
)
        ;;
    relatives)
        aws s3 cp --region "$BUCKET_REGION" "$CANDIDATES" "$s3/candidates.parquet"
        steps=$(cat <<'STEPS'
aws s3 cp $S3/candidates.parquet /data/candidates.parquet
uv run python relative_counts.py --candidates /data/candidates.parquet --work /data/rel --workers 2
rc_a=$?
mkdir -p out && cp /data/rel/relative_counts.parquet out/ 2>/dev/null
rc_b=$?
STEPS
)
        ;;
esac

cat > "$work/user-data.sh" <<EOF
#!/bin/bash
set -uo pipefail
shutdown -h +$MAX_MINUTES
exec > /var/log/exp341.log 2>&1
S3=$s3
export HOME=/root
export AWS_DEFAULT_REGION=$BUCKET_REGION
command -v aws >/dev/null || snap install aws-cli --classic
mkdir -p /opt/run /data && cd /opt/run
aws s3 cp \$S3/code.tgz . && tar xzf code.tgz
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH=/root/.local/bin:\$PATH
cd $EXP
uv sync
$steps
aws s3 cp --recursive out \$S3/out
aws s3 cp /var/log/exp341.log \$S3/exp341.log
rc=\$(( rc_a != 0 || rc_b != 0 ))
echo "a=\$rc_a b=\$rc_b" | aws s3 cp - \$S3/_DONE_rc\$rc
shutdown -h now
EOF

ami=$(aws ssm get-parameter --region "$REGION" --name "$AMI_PARAMETER" --query Parameter.Value --output text)
instance=$(aws ec2 run-instances --region "$REGION" \
    --image-id "$ami" \
    --instance-type "$INSTANCE_TYPE" \
    --subnet-id "$SUBNET" \
    --iam-instance-profile "Name=$INSTANCE_PROFILE" \
    --instance-initiated-shutdown-behavior terminate \
    --block-device-mappings "[{\"DeviceName\":\"/dev/sda1\",\"Ebs\":{\"VolumeSize\":$DISK_GB,\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]" \
    --user-data "file://$work/user-data.sh" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=exp341-$JOB-$sha}]" \
    --query 'Instances[0].InstanceId' --output text)

echo "instance $instance ($INSTANCE_TYPE, $REGION, ami $ami)"
echo "results  $s3/"
echo "done     $s3/_DONE_rc<N>"
