#!/bin/bash
# ---------------------------------------------------------------------------
# One-time AWS setup for the rustar/STARsolo/CellRanger benchmark:
#   - S3 bucket for inputs + results
#   - least-privilege IAM role + instance profile (S3 rw to that bucket only)
#   - a $5 monthly budget backstop
#   - prints the input-upload commands and the run-instances launch line
#
# Run locally with the AWS CLI configured. Idempotent-ish (ignores "already
# exists" errors). Edit BUCKET to a globally-unique name first.
# ---------------------------------------------------------------------------
set -uo pipefail
REGION="us-east-1"
BUCKET="rustar-bench"     # must be globally unique
ROLE="rustar-bench-ec2"
EMAIL="driver.ian@gmail.com"        # where the "benchmark done" email goes
TOPIC="rustar-bench"

# 1. bucket -----------------------------------------------------------------
aws s3 mb "s3://$BUCKET" --region "$REGION" || true

# 2. SNS topic + email subscription (one summary email per run) -------------
TOPIC_ARN=$(aws sns create-topic --name "$TOPIC" --region "$REGION" \
  --query TopicArn --output text)
aws sns subscribe --topic-arn "$TOPIC_ARN" --protocol email \
  --notification-endpoint "$EMAIL" --region "$REGION" || true
echo ">> Check $EMAIL and CLICK the SNS confirmation link, or no email will arrive."

# 3. IAM role + instance profile (trust EC2, allow S3 rw + SNS publish) ------
cat > /tmp/trust.json <<'J'
{"Version":"2012-10-17","Statement":[{"Effect":"Allow",
 "Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}
J
aws iam create-role --role-name "$ROLE" \
  --assume-role-policy-document file:///tmp/trust.json || true

cat > /tmp/policy.json <<J
{"Version":"2012-10-17","Statement":[
 {"Effect":"Allow","Action":["s3:GetObject","s3:PutObject","s3:ListBucket"],
  "Resource":["arn:aws:s3:::$BUCKET","arn:aws:s3:::$BUCKET/*"]},
 {"Effect":"Allow","Action":["sns:Publish"],"Resource":["$TOPIC_ARN"]}]}
J
aws iam put-role-policy --role-name "$ROLE" --policy-name s3rw \
  --policy-document file:///tmp/policy.json

aws iam create-instance-profile --instance-profile-name "$ROLE" || true
aws iam add-role-to-instance-profile \
  --instance-profile-name "$ROLE" --role-name "$ROLE" || true

# 4. $5 budget backstop (free; shows in Billing console) ---------------------
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
cat > /tmp/budget.json <<J
{"BudgetName":"rustar-bench-5usd","BudgetLimit":{"Amount":"5","Unit":"USD"},
 "TimeUnit":"MONTHLY","BudgetType":"COST"}
J
aws budgets create-budget --account-id "$ACCOUNT" \
  --budget file:///tmp/budget.json || true

# Emit the one value bench.sh needs; everything else is in README.md.
echo "Infra ready. SNS topic ARN (paste into bench.sh SNS_TOPIC_ARN):"
echo "  $TOPIC_ARN"
echo "Next steps: see README.md  (stage.sh -> edit bench.sh -> launch.sh)"
