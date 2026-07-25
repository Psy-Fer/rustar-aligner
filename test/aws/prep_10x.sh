#!/bin/bash
# Data-prep user-data: download 10x Genomics' OFFICIAL 5k Human PBMC 3' v3 FASTQs
# (CellRanger 3.0.2-era -> classic 6.8M 3M-february-2018 inclusion list, no barcode
# translation), VERIFY the cell barcodes match that whitelist (gate), subsample to
# fixed read counts (paired, synced), and stage to S3. Self-terminates.
#
# Why this dataset: the prior GEO pick (SRR39340842, "3pv3") matched NO 10x inclusion
# list (february-2018 ~7-12%, may-2023 0.1%, v2 0.1%) -> unusable for solo. A 10x-own
# v3 dataset is guaranteed to match february-2018.
# 10x 3' v3: R1 = 28 bp (CB16 + UMI12), R2 = 91 bp cDNA. soloStrand Forward.
set -uxo pipefail

BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
TARURL="https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar"
SAMPLE="human_pbmc5k_3pv3"
WL="$BUCKET/whitelists/3M-february-2018.txt.gz"
SUBSAMPLES="10000000 50000000"   # 10M matches the mouse set; 50M = scaling
SEED=11
MIN_MATCH=70                     # abort if <70% of sampled reads match the whitelist

export HOME=/root
TS=$(date -u +%Y%m%dT%H%M%SZ)

# --- mount instance-store NVMe at /root (27GB tar + extract + subsamples need space+speed)
dnf -y install nvme-cli wget tar gzip >/dev/null 2>&1 || true
INST=$(nvme list 2>/dev/null | awk '/Instance Storage/{print $1; exit}')
[ -z "$INST" ] && INST=/dev/nvme1n1
mkfs.xfs -f "$INST" >/dev/null 2>&1 && mount "$INST" /root && chmod 755 /root
cd /root
mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
exec > >(tee /var/log/prep.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/prep.log "$BUCKET/prep/${SAMPLE}_$TS.log" 2>/dev/null
  [ -n "$SNS_TOPIC_ARN" ] && aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
    --subject "10x human prep done ($TS)" \
    --message "$SAMPLE staged to $BUCKET/fastq_${SAMPLE}_*M/ (log: $BUCKET/prep/${SAMPLE}_$TS.log)" 2>/dev/null
  shutdown -h now
}
trap finish EXIT
shutdown -h +180   # 3 h dead-man's switch

echo "=== START $TS  ($SAMPLE) ==="
set -e
wget -q "https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz"
tar xf seqkit_linux_amd64.tar.gz; SEQKIT=/root/seqkit

# --- download + extract the official tarball (fast on AWS<->10x CDN)
echo "=== downloading $TARURL ==="
curl -sL --retry 3 "$TARURL" -o tenx.tar
echo "tar size: $(du -h tenx.tar | cut -f1)"
tar xf tenx.tar
DIR=$(find /root -maxdepth 1 -type d -name '*fastqs*' | head -1)
echo "extracted to $DIR"; ls -la "$DIR"

# --- concatenate lanes: R1 = barcode (28bp), R2 = cDNA (91bp)
cat "$DIR"/*_R1_*.fastq.gz > all_R1.fastq.gz
cat "$DIR"/*_R2_*.fastq.gz > all_R2.fastq.gz
L1=$(zcat all_R1.fastq.gz | sed -n '2p' | tr -d '\n' | wc -c)
L2=$(zcat all_R2.fastq.gz | sed -n '2p' | tr -d '\n' | wc -c)
echo "read lengths: R1=$L1 (expect 28)  R2=$L2 (expect 91)"

# From here on, `zcat | head` deliberately closes the pipe early -> zcat gets SIGPIPE
# (exit 141), which under `set -e` + `pipefail` would abort the whole script. Drop -e;
# the gate does its own explicit exit-on-low-match, and the subsample loop tolerates it.
set +e

# --- GATE: verify cell barcodes (R1[1:16]) match the february-2018 whitelist.
# Use an awk HASH set-membership (count EVERY read), NOT `comm`: comm pairs up
# duplicate lines, so a barcode seen 5000x but listed once matches only ONCE -> it
# collapses read-weighted match to ~distinct-matches/reads (~5%), a garbage metric.
aws s3 cp "$WL" wl.txt.gz; gzip -dc wl.txt.gz | awk '{print $1}' > wl.txt
zcat all_R1.fastq.gz | head -4000000 | awk 'NR%4==2{print substr($0,1,16)}' > cb.txt
read M T PCT <<<"$(awk 'NR==FNR{wl[$1]=1;next}{t++;if($1 in wl)m++}END{printf "%d %d %.1f",m,t,(t?100*m/t:0)}' wl.txt cb.txt)"
echo "=== WHITELIST MATCH: $PCT% ($M/$T reads, exact) ==="
if awk "BEGIN{exit !($PCT < $MIN_MATCH)}"; then
  echo "FATAL: match $PCT% < ${MIN_MATCH}% -> wrong whitelist/chemistry, aborting"; exit 1
fi
echo "match OK (>= ${MIN_MATCH}%) -> proceeding"

# --- subsample (seqkit -p STREAMS; same seed + same order => synced R1/R2 pairs)
TOTAL=$(($(zcat all_R2.fastq.gz | wc -l) / 4))
echo "total reads: $TOTAL"
set +e
for N in $SUBSAMPLES; do
  M=$((N/1000000))
  PROP=$(awk "BEGIN{p=$N/$TOTAL; print (p>1)?1:p}")
  echo "=== subsampling to ~${M}M reads (p=$PROP, seed $SEED) ==="
  R1="${SAMPLE}_S1_L001_R1_001.fastq"; R2="${SAMPLE}_S1_L001_R2_001.fastq"
  $SEQKIT sample -p "$PROP" -s "$SEED" all_R1.fastq.gz -o "$R1"
  $SEQKIT sample -p "$PROP" -s "$SEED" all_R2.fastq.gz -o "$R2"
  echo "emitted: $(($(wc -l < "$R1")/4)) reads"
  gzip -f "$R1" "$R2"
  aws s3 cp "${R1}.gz" "$BUCKET/fastq_${SAMPLE}_${M}M/" --only-show-errors
  aws s3 cp "${R2}.gz" "$BUCKET/fastq_${SAMPLE}_${M}M/" --only-show-errors
  echo "uploaded ${M}M: $(du -h ${R1}.gz ${R2}.gz | awk '{print $1}' | tr '\n' ' ')"
  rm -f "${R1}.gz" "${R2}.gz"
done
echo "=== PREP DONE ==="
