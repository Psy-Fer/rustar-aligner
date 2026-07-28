#!/bin/bash
# Profiling user-data: fresh WARM alignment profile of rustar on native x86_64,
# plus a baseline-vs-x86-64-v4 build A/B. Uses the S3-cached rust_idx (no rebuild).
# Profiles the CORE ALIGNER (plain SE align of the cDNA reads, --outSAMtype None,
# no solo) so the hotspots are seed/SA/compare/stitch, not solo quant.
# Uploads timing + perf reports to S3 and self-terminates.
set -uxo pipefail

BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
THREADS=10
RUSTAR_BRANCH="phase14-starsolo"

cd /root
export HOME=/root
mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
TS=$(date -u +%Y%m%dT%H%M%SZ)
RESULTS="$BUCKET/profile/$TS"
exec > >(tee /var/log/profile.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/profile.log "$RESULTS/" 2>/dev/null
  for f in /root/*.txt /root/time_*.log; do aws s3 cp "$f" "$RESULTS/" 2>/dev/null; done
  [ -n "$SNS_TOPIC_ARN" ] && aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
    --subject "rustar profile done ($TS)" --message "Profile results: $RESULTS/" 2>/dev/null
  shutdown -h now
}
trap finish EXIT
shutdown -h +120   # 2 h dead-man's switch (profiling is short)

echo "=== START $TS ==="
set -e
dnf -y install git gcc make tar gzip time wget perf
curl -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
git clone https://github.com/iandriver/rustar-aligner
cd /root/rustar-aligner
# Both unstripped + debug syms + frame pointers (for perf call graphs).
COMMON="-C force-frame-pointers=yes"
echo "=== building NEW (perf-hashers, with P1) ==="
git checkout perf-hashers
RUSTFLAGS="$COMMON" CARGO_PROFILE_RELEASE_STRIP=false CARGO_PROFILE_RELEASE_DEBUG=2 cargo build --release
cp target/release/rustar-aligner /root/rustar_new
echo "=== building OLD (phase14-starsolo, pre-P1) ==="
git checkout phase14-starsolo
RUSTFLAGS="$COMMON" CARGO_PROFILE_RELEASE_STRIP=false CARGO_PROFILE_RELEASE_DEBUG=2 cargo build --release
cp target/release/rustar-aligner /root/rustar_old
cd /root

# Cached index + reads from S3 (no rebuild).
aws s3 sync "$BUCKET/idx/rust_idx/" /root/rust_idx --only-show-errors
aws s3 cp "$BUCKET/fastq/5k_Mouse_PBMCs_5p_gem-x_GEX_S1_L001_R2_001.fastq.gz" /root/R2.fastq.gz --only-show-errors

ALIGN=(--genomeDir /root/rust_idx --readFilesIn /root/R2.fastq.gz --runThreadN "$THREADS" --outSAMtype None)

set +e
# Warm the page cache (untimed) so we profile compute, not cold-EBS I/O.
/root/rustar_new "${ALIGN[@]}" --outFileNamePrefix /root/warm/ >/dev/null 2>&1

# A/B: old (pre-P1) vs new (P1), 3 warm reps each.
for b in old new; do
  for r in 1 2 3; do
    /usr/bin/time -v /root/rustar_$b "${ALIGN[@]}" \
      --outFileNamePrefix /root/o_${b}_${r}/ &> /root/time_${b}_${r}.log
  done
done

# perf profile the NEW build (warm). FLAT self-time (`-g none` = no call trees, so
# the ranking isn't truncated) + a separate call-graph view.
perf record -g -F 999 -o /root/perf.data -- \
  /root/rustar_new "${ALIGN[@]}" --outFileNamePrefix /root/perfout/ 2> /root/perf_record.txt
perf report -i /root/perf.data --stdio --no-children -g none 2>/dev/null | head -45 > /root/perf_self.txt
perf report -i /root/perf.data --stdio 2>/dev/null | head -120 > /root/perf_callgraph.txt

echo "=== DONE $(date -u) ==="
