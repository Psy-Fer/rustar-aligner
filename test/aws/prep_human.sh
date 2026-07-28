#!/bin/bash
# Data-prep user-data: download a human 10x 3' v3 run from SRA, subsample to fixed
# read counts (paired, synced), and stage to S3 for the benchmark. Self-terminates.
# Source: GSM9842100 / SRP712988, run SRR39334... -> SRR39340842 (~106M reads).
# 10x 3' v3: R1 = 28 bp (CB16 + UMI12), R2 = 91 bp cDNA. soloStrand Forward.
set -uxo pipefail

BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
SRR="SRR39340842"
SAMPLE="human_pancreas_3pv3"
SUBSAMPLES="10000000 50000000"   # reads to emit (10M matches the mouse set; 50M = scaling)
SEED=11

cd /root
export HOME=/root
mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
mkdir -p /root/.ncbi
echo '/LIBS/GUID = "a1b2c3d4-e5f6-7890-abcd-ef1234567890"' > /root/.ncbi/user-settings.mkfg
TS=$(date -u +%Y%m%dT%H%M%SZ)
exec > >(tee /var/log/prep.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/prep.log "$BUCKET/prep/$TS.log" 2>/dev/null
  [ -n "$SNS_TOPIC_ARN" ] && aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
    --subject "human data prep done ($TS)" \
    --message "Subsampled $SAMPLE to S3: $BUCKET/fastq_human_*M/  (log: $BUCKET/prep/$TS.log)" 2>/dev/null
  shutdown -h now
}
trap finish EXIT
shutdown -h +180   # 3 h dead-man's switch

echo "=== START $TS ==="
set -e
dnf -y install wget tar gzip

# SRA toolkit (prefetch + fasterq-dump) + seqkit (paired subsample)
wget -q "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
tar xf sratoolkit.current-centos_linux64.tar.gz
export PATH="$PATH:$(ls -d /root/sratoolkit.*/bin | head -1)"
wget -q "https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz"
tar xf seqkit_linux_amd64.tar.gz
SEQKIT=/root/seqkit

# Download + extract. --include-technical is ESSENTIAL: 10x R1 (CB+UMI) is flagged
# technical and fasterq-dump would otherwise drop it.
prefetch "$SRR" -O /root/sra --max-size 100g
fasterq-dump "/root/sra/$SRR/$SRR.sra" --split-files --include-technical \
  -e "$(nproc)" -t /root/rtmp -O /root/fq

# Identify barcode (28 bp -> R1) vs cDNA (~91 bp -> R2) by read length.
len() { head -2 "$1" | tail -1 | tr -d '\n' | wc -c; }
L1=$(len /root/fq/${SRR}_1.fastq); L2=$(len /root/fq/${SRR}_2.fastq)
echo "read lengths: _1=$L1 _2=$L2"
if [ "$L1" -le 40 ]; then BC=/root/fq/${SRR}_1.fastq; CDNA=/root/fq/${SRR}_2.fastq;
else BC=/root/fq/${SRR}_2.fastq; CDNA=/root/fq/${SRR}_1.fastq; fi
echo "barcode read: $BC   cDNA read: $CDNA"

# seqkit `-n` (by number) loads ALL reads into RAM (OOM on 106M reads). Use `-p`
# (by proportion) which STREAMS; same seed + same input order => synced R1/R2 pairs.
# Counts come out ~N (±sqrt(N)), fine for a benchmark.
TOTAL=$(($(wc -l < "$CDNA") / 4))
echo "total reads: $TOTAL"
set +e
for N in $SUBSAMPLES; do
  M=$((N/1000000))
  PROP=$(awk "BEGIN{p=$N/$TOTAL; print (p>1)?1:p}")
  echo "=== subsampling to ~${M}M reads (p=$PROP, seed $SEED, synced pairs) ==="
  R1="${SAMPLE}_S1_L001_R1_001.fastq"; R2="${SAMPLE}_S1_L001_R2_001.fastq"
  $SEQKIT sample -p "$PROP" -s "$SEED" "$BC"   -o "$R1"
  $SEQKIT sample -p "$PROP" -s "$SEED" "$CDNA" -o "$R2"
  echo "emitted: $(($(wc -l < "$R1")/4)) reads"
  gzip -f "$R1" "$R2"
  aws s3 cp "${R1}.gz" "$BUCKET/fastq_human_${M}M/" --only-show-errors
  aws s3 cp "${R2}.gz" "$BUCKET/fastq_human_${M}M/" --only-show-errors
  echo "uploaded ${M}M: $(du -h ${R1}.gz ${R2}.gz | awk '{print $1}' | tr '\n' ' ')"
  rm -f "${R1}.gz" "${R2}.gz"
done

echo "=== DONE $(date -u) ==="
