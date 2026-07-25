#!/bin/bash
# Memory-breakdown measurement on Linux (the platform where the 42 GB Max RSS was
# seen). Runs the mouse solo workload and samples /proc/<pid>/status to split
# RssFile (reclaimable mmap'd index) from RssAnon (real heap), under two mimalloc
# purge settings (-1 = never, default; 1000 ms = capped reserve) to quantify the
# purge/RSS/speed tradeoff. Uses the cached index. Self-terminates.
set -uxo pipefail

BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
THREADS=10

cd /root
export HOME=/root
mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
TS=$(date -u +%Y%m%dT%H%M%SZ)
exec > >(tee /var/log/mem.log) 2>&1
finish() {
  set +e
  aws s3 cp /var/log/mem.log "$BUCKET/mem/$TS.log" 2>/dev/null
  [ -n "$SNS_TOPIC_ARN" ] && aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
    --subject "rustar mem breakdown done ($TS)" --message "Log: $BUCKET/mem/$TS.log" 2>/dev/null
  shutdown -h now
}
trap finish EXIT
shutdown -h +120

echo "=== START $TS ==="
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
R1=$(ls /root/fastq/*R1_001.fastq.gz | head -1); R2=$(ls /root/fastq/*R2_001.fastq.gz | head -1)

SOLO=(--soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Reverse
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000)

# /proc sampler: track peak VmHWM (=Max RSS), RssAnon (real), RssFile (reclaimable).
sample_peak() { # $1 = pid
  local p=$1 a=0 f=0 h=0
  while kill -0 "$p" 2>/dev/null; do
    while read -r k v _; do
      case "$k" in
        VmHWM:)   [ "$v" -gt "$h" ] && h=$v ;;
        RssAnon:) [ "$v" -gt "$a" ] && a=$v ;;
        RssFile:) [ "$v" -gt "$f" ] && f=$v ;;
      esac
    done < <(grep -E "VmHWM:|RssAnon:|RssFile:" "/proc/$p/status" 2>/dev/null)
    sleep 1
  done
  echo "PEAK  MaxRSS(VmHWM)=$((h/1048576))G  RssAnon(real)=$((a/1048576))G  RssFile(reclaimable)=$((f/1048576))G"
}

set +e
for PURGE in -1 1000; do
  echo "=== mimalloc purge_delay=${PURGE}ms ==="
  rm -rf /root/out_$PURGE; mkdir -p /root/out_$PURGE
  START=$(date +%s)
  RUSTAR_PURGE_DELAY_MS=$PURGE "$RUSTAR" \
    --genomeDir rust_idx --sjdbGTFfile mouse_genes.gtf \
    --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
    --soloOutGzip yes --outSAMtype None --outFileNamePrefix "/root/out_$PURGE/" \
    > "/root/run_$PURGE.log" 2>&1 &
  PID=$!            # rustar directly (not wrapped) so /proc/$PID is the real process
  sample_peak "$PID"
  wait "$PID"
  echo "wall: $(( $(date +%s) - START ))s"
done

echo "=== DONE $(date -u) ==="
