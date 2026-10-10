#!/usr/bin/env bash
# Launch one self-terminating EC2 GPU instance for the exp341 harness checks.
#
# The instance runs harness checks 1-2 at the plan's 4 seeds on CUDA
# (check_harness.py) and check 5, the timed clean readout over eval-val and
# eval-denovo (time_readout.py), copies everything under out/ plus its log to
# S3, writes a _DONE_rc<N> marker, and shuts down (shutdown behaviour is
# terminate). A `shutdown -h +MAX_MINUTES` scheduled first caps the cost if
# anything hangs.
#
# No secrets reach the instance: code comes from S3 through the instance
# role, and the checkpoint from the public MarinFold HF bucket.
#
# Account facts (no key pairs; the instance profile name is not discoverable
# with this user's IAM permissions) are from the exp91/exp145 EC2 workers.
#
# Usage (CloudShell STS credentials expire after ~15 minutes; REGION and
# SUBNET must name the same region):
#   set -a; source experiments/exp341_evals_contacts_v1_interpretability/aws_keys.env; set +a
#   experiments/exp341_evals_contacts_v1_interpretability/launch_aws.sh
#   INSTANCE_TYPE=g6e.xlarge experiments/exp341_evals_contacts_v1_interpretability/launch_aws.sh
set -euo pipefail

INSTANCE_TYPE=${INSTANCE_TYPE:-p5.4xlarge}  # 1x H100 80GB
MAX_MINUTES=${MAX_MINUTES:-120}
# The instance's region. GPU quota (2026-10): us-east-2 has 192 P / 128 G
# on-demand vCPUs; us-west-2 has 4 P / 0 G, too few for any GPU type.
REGION=${REGION:-us-east-2}
# The bucket stays in us-west-2; code in and results out are tens of MB.
BUCKET=marinfold-exp91-usw2
BUCKET_REGION=us-west-2
# Default subnets: us-east-2a subnet-03f9b7abd4284ad5e, 2b subnet-02c630af81a9bdc13,
# 2c subnet-091e8c3572c921332; us-west-2a subnet-0a393c46ff06d81c3,
# 2b subnet-0563fac5b8bee3f0d. Capacity varies by zone, so SUBNET picks one.
SUBNET=${SUBNET:-subnet-03f9b7abd4284ad5e}
INSTANCE_PROFILE=marinfold-exp91-instance-profile
EXP=experiments/exp341_evals_contacts_v1_interpretability
EXP245=experiments/exp245_evals_foldbench_held_out_monomers/data

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"
if ! git diff --quiet HEAD -- marinfold "$EXP"; then
    echo "uncommitted changes under marinfold/ or $EXP; the instance runs HEAD, commit first" >&2
    exit 1
fi
sha=$(git rev-parse --short HEAD)
prefix="exp341/harness-checks/$(date -u +%Y%m%dT%H%M%SZ)-$sha"
s3="s3://$BUCKET/$prefix"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

git archive --format=tar.gz -o "$work/code.tgz" HEAD \
    marinfold "$EXP" "$EXP245/eval_sets.csv" "$EXP245/gt_universe_scored.jsonl"
aws s3 cp --region "$BUCKET_REGION" "$work/code.tgz" "$s3/code.tgz"

ami=$(aws ssm get-parameter --region "$REGION" \
    --name /aws/service/deeplearning/ami/x86_64/base-oss-nvidia-driver-gpu-ubuntu-22.04/latest/ami-id \
    --query Parameter.Value --output text)

cat > "$work/user-data.sh" <<EOF
#!/bin/bash
set -uo pipefail
shutdown -h +$MAX_MINUTES
exec > /var/log/exp341.log 2>&1
S3=$s3
export HOME=/root
export AWS_DEFAULT_REGION=$BUCKET_REGION
command -v aws >/dev/null || snap install aws-cli --classic
mkdir -p /opt/run && cd /opt/run
aws s3 cp \$S3/code.tgz . && tar xzf code.tgz
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH=/root/.local/bin:\$PATH
cd $EXP
nvidia-smi
uv sync
uv run python check_harness.py --device cuda --seeds 4 --output out/harness_checks_cuda_4seed.json
rc_checks=\$?
uv run python time_readout.py --device cuda --runner-tag aws --output-dir out/readout
rc_readout=\$?
aws s3 cp --recursive out \$S3/out
aws s3 cp /var/log/exp341.log \$S3/exp341.log
rc=\$(( rc_checks != 0 || rc_readout != 0 ))
echo "checks=\$rc_checks readout=\$rc_readout" | aws s3 cp - \$S3/_DONE_rc\$rc
shutdown -h now
EOF

instance=$(aws ec2 run-instances --region "$REGION" \
    --image-id "$ami" \
    --instance-type "$INSTANCE_TYPE" \
    --subnet-id "$SUBNET" \
    --iam-instance-profile "Name=$INSTANCE_PROFILE" \
    --instance-initiated-shutdown-behavior terminate \
    --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":150,"VolumeType":"gp3","DeleteOnTermination":true}}]' \
    --user-data "file://$work/user-data.sh" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=exp341-harness-$sha}]" \
    --query 'Instances[0].InstanceId' --output text)

echo "instance $instance ($INSTANCE_TYPE, ami $ami)"
echo "results  $s3/"
echo "done     $s3/_DONE_rc<N>"
