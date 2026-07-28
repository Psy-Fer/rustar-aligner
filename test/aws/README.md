# rustar vs STARsolo vs CellRanger — AWS benchmark

A fair, native-x86_64 wall-time + memory comparison on a self-terminating EC2
Spot instance. Docker-on-Mac can't compare these fairly (CellRanger is x86-only,
STAR can't run native on Apple Silicon, and rustar's mmap'd index is penalized by
Rosetta+virtiofs) — a real Linux x86 box is the only clean substrate.

## Files
| file | what it is | run where |
|---|---|---|
| `setup_aws.sh` | one-time infra: S3 bucket, SNS topic+email, IAM role/profile, $5 budget | local |
| `stage.sh` | upload inputs (CellRanger tarball, ref, whitelist, FASTQs) to S3 | local |
| `launch.sh` | launch the self-terminating Spot instance | local |
| `bench.sh` | the benchmark itself (EC2 user-data); builds tools+indexes, runs 3× | on instance |
| `scrape_results.py` | turn the logs into a comparison table | local |
| `bench-user-policy.json` | least-privilege IAM policy for the launching user | local |

## Prerequisites
- AWS CLI configured (`aws configure` or `aws configure sso`); verify with
  `aws sts get-caller-identity`. Don't use the root account.
- Your CellRanger 10.0.0 tarball + mouse reference tarball locally.

## Run order
```bash
cd test/aws

# 1. one-time infra (prints the SNS topic ARN)
bash setup_aws.sh
#    -> click the SNS confirmation link in your email, or no "done" email arrives

# 2. edit stage.sh paths, then upload inputs
bash stage.sh

# 3. edit bench.sh:  BUCKET="s3://<bucket>"  and  SNS_TOPIC_ARN="<from step 1>"
#    (also check SAMPLE / REF_TGZ / OVERHANG)

# 4. launch (self-terminating Spot; ~$1-2; emails you when done)
bash launch.sh

# 5. when the email lands, pull + summarize
aws s3 cp s3://<bucket>/results/ ./results --recursive
python3 scrape_results.py ./results/<timestamp>
```

## Fairness controls (baked into bench.sh)
- All three native x86_64, same instance, same data, `--runThreadN 10`.
- No BAM for any tool (`--outSAMtype None` / `--create-bam=false --nosecondary`).
- STAR + rustar indexes built from CellRanger's own `genome.fa`/`genes.gtf`.
- Page cache dropped before each run; 3 reps for variance.
- Captures wall (`Elapsed`) + Max RSS via `/usr/bin/time -v`. rustar's RSS
  includes reclaimable mmap pages, so it overstates real memory pressure — noted
  in the output table.

## Cost / safety
- Spot `m7i.4xlarge` (16 vCPU / 64 GB), self-terminates on completion + 4 h
  dead-man's switch. ~$1-2 total. The `$5` budget alarm is a backstop.
- Least-privilege IAM throughout; delete the launching user's access key when done.
