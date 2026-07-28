#!/bin/bash
# ---------------------------------------------------------------------------
# Detailed x86_64 profiling run to decide whether further alignment-speed
# work is worth it. Reuses the S3-cached dense mouse index (no genomeGenerate
# needed — on-disk format is unchanged this session, only content when
# --genomeSAsparseD>1 is requested) so this is fast (~20-30 min).
#
# Profiles the CORE ALIGNER (plain SE, --outSAMtype None, no solo) so hotspots
# are seed/SA/compare/stitch, not solo quant noise — same convention as the
# existing profile.sh.
#
# Three binaries:
#   OLD      — pinned pre-session commit ($OLD_COMMIT), plain release build
#   NEW      — this session's working tree (uncommitted), plain release build
#              [includes: zlib-rs, FxHash collation+stitch-finalization,
#               genome MADV_RANDOM, SA prefetch, get_base #[inline] fix]
#   NEW_PGO  — NEW, profile-guided-optimized via scripts/pgo-build.sh
#
# Deliverables:
#   - timing A/B: OLD vs NEW vs NEW_PGO (3 warm reps each)
#   - perf record -g on NEW: flat self-time + call-graph views
#   - perf stat on NEW: IPC, branch mispredict rate, cache miss rate
# ---------------------------------------------------------------------------
set -uxo pipefail

BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
THREADS=10
OLD_COMMIT="473ccc43a073df1aa70c25616eacd7bd3b5a545a"
NEW_SRC_S3="$BUCKET/src/rustar-new-src.tar.gz"

# --- Mount instance-store NVMe at /root (26 GB cached index + 3 target dirs
# + PGO training scratch would overflow the default small EBS root volume
# otherwise — see bench_oldnew.sh for the same pattern).
dnf -y install nvme-cli >/dev/null 2>&1 || true
INST=$(nvme list 2>/dev/null | awk '/Instance Storage/{print $1; exit}')
[ -z "$INST" ] && INST=/dev/nvme1n1
mkfs.xfs -f "$INST"
mount "$INST" /root
chmod 755 /root

cd /root
export HOME=/root
df -h /root
mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
TS=$(date -u +%Y%m%dT%H%M%SZ)
RESULTS="$BUCKET/profile/detailed-$TS"
exec > >(tee /var/log/profile2.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/profile2.log "$RESULTS/" 2>/dev/null
  for f in /root/*.log /root/*.txt; do aws s3 cp "$f" "$RESULTS/" 2>/dev/null; done
  if [ -n "$SNS_TOPIC_ARN" ]; then
    SUMMARY=$(for f in /root/time_*.log; do
      [ -f "$f" ] || continue
      w=$(grep -m1 'Elapsed (wall clock)' "$f" | sed 's/.*: //')
      echo "$(basename "$f" .log): wall=${w:-FAILED}"
    done)
    aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
      --subject "rustar detailed profile done ($TS)" \
      --message "Results: $RESULTS/

$SUMMARY

Pull:
  aws s3 cp $RESULTS/ ./results --recursive" 2>/dev/null
  fi
  shutdown -h now
}
trap finish EXIT
shutdown -h +120   # 2 h dead-man's switch

echo "=== START $TS ==="
set -e
dnf -y install git gcc make tar gzip time wget perf
curl -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
rustup component add llvm-tools-preview

COMMON="-C force-frame-pointers=yes"

# ---- OLD ----
git clone https://github.com/iandriver/rustar-aligner old_src
( cd old_src && git checkout "$OLD_COMMIT" && \
  RUSTFLAGS="$COMMON" CARGO_PROFILE_RELEASE_STRIP=false CARGO_PROFILE_RELEASE_DEBUG=2 \
  cargo build --release ) &> build_old.log
cp old_src/target/release/rustar-aligner /root/rustar_old

# ---- NEW ----
aws s3 cp "$NEW_SRC_S3" new_src.tar.gz
mkdir -p new_src && tar xzf new_src.tar.gz -C new_src --strip-components=1
( cd new_src && \
  RUSTFLAGS="$COMMON" CARGO_PROFILE_RELEASE_STRIP=false CARGO_PROFILE_RELEASE_DEBUG=2 \
  cargo build --release ) &> build_new.log
cp new_src/target/release/rustar-aligner /root/rustar_new

# ---- NEW_PGO (plain, no debug/frame-pointer flags — production-shaped build) ----
aws s3 cp "$BUCKET/idx/rust_idx/" /root/pgo_train_idx --recursive --only-show-errors
aws s3 cp "$BUCKET/fastq/5k_Mouse_PBMCs_5p_gem-x_GEX_S1_L001_R2_001.fastq.gz" /root/pgo_train.fastq.gz --only-show-errors
( cd new_src && bash scripts/pgo-build.sh /root/pgo_train_idx /root/pgo_train.fastq.gz 1500000 ) &> build_new_pgo.log
cp new_src/target/release/rustar-aligner /root/rustar_new_pgo

# ---- cached dense index + reads (no genomeGenerate — on-disk format unchanged) ----
aws s3 sync "$BUCKET/idx/rust_idx/" /root/rust_idx --only-show-errors
aws s3 cp "$BUCKET/fastq/5k_Mouse_PBMCs_5p_gem-x_GEX_S1_L001_R2_001.fastq.gz" /root/R2.fastq.gz --only-show-errors

ALIGN=(--genomeDir /root/rust_idx --readFilesIn /root/R2.fastq.gz --runThreadN "$THREADS" --outSAMtype None)

set +e
# Warm the page cache (untimed).
/root/rustar_new "${ALIGN[@]}" --outFileNamePrefix /root/warm/ >/dev/null 2>&1

for b in old new new_pgo; do
  BIN="/root/rustar_$b"
  for r in 1 2 3; do
    mkdir -p "/root/out_${b}_${r}"
    /usr/bin/time -v "$BIN" "${ALIGN[@]}" \
      --outFileNamePrefix "/root/out_${b}_${r}/" &> "time_${b}_${r}.log"
    echo "$b rep$r rc=$?"
  done
done

# ---- perf record: NEW, warm, flat self-time + call-graph ----
perf record -g -F 999 -o /root/perf.data -- \
  /root/rustar_new "${ALIGN[@]}" --outFileNamePrefix /root/perfout/ 2> /root/perf_record.txt
perf report -i /root/perf.data --stdio --no-children -g none 2>/dev/null | head -60 > /root/perf_self.txt
perf report -i /root/perf.data --stdio 2>/dev/null | head -150 > /root/perf_callgraph.txt

# ---- perf stat: IPC, branch mispredicts, cache misses (NEW, warm) ----
perf stat -e cycles,instructions,branches,branch-misses,cache-references,cache-misses \
  -o /root/perf_stat.txt -- \
  /root/rustar_new "${ALIGN[@]}" --outFileNamePrefix /root/perfstat_out/ 2>&1 | tee -a /root/perf_stat.txt

echo "=== DONE $(date -u) ==="
