#!/bin/bash
# Launch the self-terminating Spot instance that runs bench.sh as user-data.
# Run from this directory (it references file://bench.sh). Edit bench.sh first
# (BUCKET + SNS_TOPIC_ARN), then: bash launch.sh
set -euo pipefail

############################ EDIT IF NEEDED ############################
REGION="us-east-1"
ROLE="rustar-bench-ec2"          # instance profile created by setup_aws.sh
INSTANCE_TYPE="m7i.4xlarge"      # 16 vCPU / 64 GB
DISK_GB=200
######################################################################

cd "$(dirname "$0")"
[ -f bench.sh ] || { echo "bench.sh not found in $(pwd)"; exit 1; }
if grep -q 'CHANGE-ME' bench.sh; then
  echo "Edit bench.sh first: set BUCKET (and SNS_TOPIC_ARN)."; exit 1
fi

ID=$(aws ec2 run-instances --region "$REGION" \
  --image-id resolve:ssm:/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 \
  --instance-type "$INSTANCE_TYPE" \
  --instance-market-options '{"MarketType":"spot"}' \
  --instance-initiated-shutdown-behavior terminate \
  --iam-instance-profile "Name=$ROLE" \
  --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":$DISK_GB,\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]" \
  --user-data file://bench.sh \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=rustar-bench}]' \
  --query 'Instances[0].InstanceId' --output text)

echo "Launched $ID. It self-terminates when done; you'll get an SNS email."
echo "Watch:  aws ec2 describe-instances --instance-ids $ID --query 'Reservations[0].Instances[0].State.Name' --output text"
