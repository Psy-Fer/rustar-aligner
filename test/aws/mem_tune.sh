#!/bin/bash
# The ~17 GB RssAnon is up-front mimalloc committed reserve (live data ~3 GB), not
# solo data. Test mimalloc COMMIT knobs (via MIMALLOC_* env, read at init) to find a
# setting that shrinks committed RSS at acceptable speed cost. purge stays -1 (FFI).
# Mouse solo, cached index. Self-terminates.
set -uxo pipefail
BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
cd /root; export HOME=/root; mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
TS=$(date -u +%Y%m%dT%H%M%SZ)
exec > >(tee /var/log/tune.log) 2>&1
finish(){ set +e; aws s3 cp /var/log/tune.log "$BUCKET/mem/tune_$TS.log" 2>/dev/null
  [ -n "$SNS_TOPIC_ARN" ] && aws sns publish --topic-arn "$SNS_TOPIC_ARN" --subject "rustar mem tune done ($TS)" --message "$BUCKET/mem/tune_$TS.log" 2>/dev/null
  shutdown -h now; }
trap finish EXIT
shutdown -h +120
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
SOLO=(--soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt --soloCBstart 1 --soloCBlen 16
  --soloUMIstart 17 --soloUMIlen 12 --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Reverse
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000
  --clipAdapterType CellRanger4 --outFilterScoreMin 30)

peak() { local p=$1 a=0 h=0; while kill -0 "$p" 2>/dev/null; do
    v=$(awk '/RssAnon:/{print $2}' /proc/$p/status 2>/dev/null); [ "${v:-0}" -gt "$a" ] && a=$v
    w=$(awk '/VmHWM:/{print $2}' /proc/$p/status 2>/dev/null); [ "${w:-0}" -gt "$h" ] && h=$w
    sleep 1; done; echo "$((a/1048576)) $((h/1048576))"; }

run() { # $1=label  $2..=env assignments
  local label=$1; shift
  rm -rf /root/o; mkdir -p /root/o
  local s=$(date +%s)
  env RUSTAR_PURGE_DELAY_MS=-1 "$@" "$RUSTAR" --genomeDir rust_idx --sjdbGTFfile mouse_genes.gtf \
    --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN 10 --soloOutGzip yes --outSAMtype None \
    --outFileNamePrefix /root/o/ > /root/o.log 2>&1 &
  local pid=$!; read anon hwm < <(peak "$pid"); wait "$pid"
  echo "RESULT  $label  RssAnon=${anon}G  MaxRSS=${hwm}G  wall=$(( $(date +%s)-s ))s"
}

set +e
echo "###### MIMALLOC COMMIT TUNING ######"
run "baseline(eager)"
run "arena_eager_commit=0"        MIMALLOC_ARENA_EAGER_COMMIT=0
run "eager_commit_delay=16"       MIMALLOC_EAGER_COMMIT_DELAY=16
run "both"                        MIMALLOC_ARENA_EAGER_COMMIT=0 MIMALLOC_EAGER_COMMIT_DELAY=16
run "eager_commit=0"              MIMALLOC_EAGER_COMMIT=0
echo "=== DONE $(date -u) ==="
