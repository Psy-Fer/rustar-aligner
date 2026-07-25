#!/bin/bash
# Localize rustar's ~17 GB working heap: sample RssAnon/RssFile every 1s alongside
# the run phase (processed-N / matrix-build), so we can see WHEN it peaks (alignment
# vs the per-feature matrix build) and correlate with the logged per-feature record
# counts. Mouse solo, cached index, default purge (-1). Self-terminates.
set -uxo pipefail
BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
cd /root; export HOME=/root; mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
TS=$(date -u +%Y%m%dT%H%M%SZ)
exec > >(tee /var/log/traj.log) 2>&1
finish() { set +e
  aws s3 cp /var/log/traj.log "$BUCKET/mem/traj_$TS.log" 2>/dev/null
  aws s3 cp /root/traj.txt "$BUCKET/mem/traj_$TS.txt" 2>/dev/null
  aws s3 cp /root/run.log "$BUCKET/mem/trajrun_$TS.log" 2>/dev/null
  [ -n "$SNS_TOPIC_ARN" ] && aws sns publish --topic-arn "$SNS_TOPIC_ARN" --subject "rustar mem trajectory done ($TS)" --message "$BUCKET/mem/traj_$TS.*" 2>/dev/null
  shutdown -h now; }
trap finish EXIT
shutdown -h +90
set -e
dnf -y install git gcc make tar gzip
curl -sSf https://sh.rustup.rs | sh -s -- -y; source "$HOME/.cargo/env"
git clone -b perf-hashers https://github.com/iandriver/rustar-aligner
( cd rustar-aligner && cargo build --release )
RUSTAR=/root/rustar-aligner/target/release/rustar-aligner
aws s3 sync "$BUCKET/idx/rust_idx/" /root/rust_idx --only-show-errors
aws s3 cp "$BUCKET/mouse_genes.gtf.gz" . && gunzip -f mouse_genes.gtf.gz
aws s3 cp "$BUCKET/whitelist.txt" /root/whitelist.txt
aws s3 cp "$BUCKET/fastq/" /root/fastq --recursive --exclude "human*"
R1=$(ls /root/fastq/*R1_001.fastq.gz|head -1); R2=$(ls /root/fastq/*R2_001.fastq.gz|head -1)

set +e
mkdir -p /root/out
"$RUSTAR" --genomeDir rust_idx --sjdbGTFfile mouse_genes.gtf --readFilesIn "$R2" "$R1" \
  --soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Reverse \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000 \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 --runThreadN 10 \
  --soloOutGzip yes --outSAMtype None --outFileNamePrefix /root/out/ > /root/run.log 2>&1 &
PID=$!
T0=$(date +%s)
echo "elapsed RssAnon_G RssFile_G phase" > /root/traj.txt
while kill -0 "$PID" 2>/dev/null; do
  a=$(awk '/RssAnon:/{print $2}' /proc/$PID/status 2>/dev/null)
  f=$(awk '/RssFile:/{print $2}' /proc/$PID/status 2>/dev/null)
  ph=$(grep -oE 'processed [0-9]+|barcode stats|aligning|wrote (Gene|GeneFull|SJ|Velocyto)[^ ]*|collected [0-9]+|EmptyDrops_CR:' /root/run.log | tail -1 | tr ' ' '_')
  printf "%ss %s %s %s\n" "$(( $(date +%s)-T0 ))" "$(awk "BEGIN{printf \"%.1f\", ${a:-0}/1048576}")" "$(awk "BEGIN{printf \"%.1f\", ${f:-0}/1048576}")" "${ph:-start}" >> /root/traj.txt
  sleep 1
done
echo "=== per-feature record counts (the recorder memory) ==="
grep -E "collected [0-9]+ resolved|wrote .* matrix" /root/run.log
echo "=== trajectory peak ==="
sort -k2 -n /root/traj.txt | tail -5
echo "=== DONE $(date -u) ==="
