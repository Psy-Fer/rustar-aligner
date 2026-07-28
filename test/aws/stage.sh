#!/bin/bash
# Upload the benchmark inputs to S3 (run once, locally). Edit the paths below
# to point at your local files, then: bash stage.sh
set -euo pipefail

############################ EDIT THESE ############################
BUCKET="rustar-bench"                                   # same name used in setup_aws.sh
SAMPLE="5k_Mouse_PBMCs_5p_gem-x_GEX"                    # FASTQ prefix before _S1_L001_*
CR_TGZ="$HOME/Downloads/cellranger-10.0.0.tar.gz"       # your CellRanger tarball
REF_TGZ="$HOME/Downloads/refdata-gex-GRCm39-2024-A.tar.gz" # your CellRanger mouse ref
WHITELIST="/Users/iandriver/solo_bench/bench/whitelist_5pgex.txt"
FASTQ_DIR="/Users/iandriver/solo_bench/sub10"
###################################################################

aws s3 cp "$CR_TGZ"    "s3://$BUCKET/"
aws s3 cp "$REF_TGZ"   "s3://$BUCKET/"
aws s3 cp "$WHITELIST" "s3://$BUCKET/whitelist.txt"
aws s3 cp "$FASTQ_DIR/" "s3://$BUCKET/fastq/" --recursive \
  --exclude "*" --include "${SAMPLE}_S1_L001_R*_001.fastq.gz"

echo "Uploaded inputs to s3://$BUCKET/"
