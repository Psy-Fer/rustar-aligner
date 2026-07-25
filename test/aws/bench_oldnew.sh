#!/bin/bash
# ---------------------------------------------------------------------------
# rustar OLD (pre-session baseline) vs NEW (this session's changes) —
# native x86_64, self-terminating Spot instance, local NVMe for all I/O.
#
# Three variants, identical STARsolo-style mouse 5' GEM-X workload:
#   old_d1 — baseline commit $OLD_COMMIT, dense SA (D=1)
#   new_d1 — this session's changes, dense SA (D=1)
#            [decode/writer 3-stage pipeline, zlib-rs gzip backend, FxHash
#             solo collation, genome MADV_RANDOM, SA-search prefetch]
#   new_d2 — this session's changes, sparse SA (--genomeSAsparseD 2)
#            [expects ~half the SA file vs d1, STAR-faithful construction]
#
# "old" is built from a fresh git clone pinned to $OLD_COMMIT (the public
# fork tip at session start). "new" is built from a source tarball of the
# CURRENT WORKING TREE (uncommitted changes) staged to S3 beforehand —
# no commits are made to get it there.
#
# Skips STAR/CellRanger this round (this is an isolated rustar-vs-rustar
# regression check); re-run the full bench_nvme.sh 3-way if these numbers
# move enough to be worth re-establishing the STARsolo/CellRanger comparison.
# ---------------------------------------------------------------------------
set -uxo pipefail

############################ CONFIG — EDIT THESE ############################
BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
THREADS=16                                    # m6id.4xlarge = 16 vCPU
SAMPLE="5k_Mouse_PBMCs_5p_gem-x_GEX"
REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
OVERHANG=89                                   # R2 read length 90 - 1
OLD_COMMIT="473ccc43a073df1aa70c25616eacd7bd3b5a545a"   # fork/phase14-starsolo tip pre-session
NEW_SRC_S3="$BUCKET/src/rustar-new-src.tar.gz"          # uploaded working-tree tarball
############################################################################

# --- Mount instance-store NVMe at /root (see bench_nvme.sh for rationale).
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
RESULTS="$BUCKET/results/oldnew-$TS"
exec > >(tee /var/log/bench.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/bench.log "$RESULTS/" 2>/dev/null
  for f in /root/*.log /root/*.txt; do aws s3 cp "$f" "$RESULTS/" 2>/dev/null; done
  if [ -n "$SNS_TOPIC_ARN" ]; then
    SUMMARY=$(for f in /root/rustar_*.log; do
      [ -f "$f" ] || continue
      w=$(grep -m1 'Elapsed (wall clock)' "$f" | sed 's/.*: //')
      rss=$(grep -m1 'Maximum resident set size' "$f" | sed 's/.*: //')
      echo "$(basename "$f" .log): wall=${w:-FAILED} maxRSS_KB=${rss:-?}"
    done)
    SIZES=$(cat /root/index_sizes.txt 2>/dev/null)
    aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
      --subject "rustar OLD-vs-NEW benchmark done ($TS)" \
      --message "Results: $RESULTS/

Index sizes:
$SIZES

Wall time + Max RSS per run:
$SUMMARY

Pull:
  aws s3 cp $RESULTS/ ./results --recursive" 2>/dev/null
  fi
  shutdown -h now
}
trap finish EXIT
shutdown -h +240   # dead-man's switch

echo "=== START $TS  (threads=$THREADS) ==="

set -e
dnf -y install git gcc make tar gzip time wget unzip which
curl -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

# ---- OLD: fresh clone pinned to the pre-session commit ----
git clone https://github.com/iandriver/rustar-aligner old_src
( cd old_src && git checkout "$OLD_COMMIT" && cargo build --release ) &> build_old_bin.log
RUSTAR_OLD=/root/old_src/target/release/rustar-aligner

# ---- NEW: this session's uncommitted working tree, staged via S3 tarball ----
aws s3 cp "$NEW_SRC_S3" new_src.tar.gz
mkdir -p new_src && tar xzf new_src.tar.gz -C new_src --strip-components=1
( cd new_src && cargo build --release ) &> build_new_bin.log
RUSTAR_NEW=/root/new_src/target/release/rustar-aligner

# ---- reference + inputs ----
wget -q "$REF_URL"
tar xf "$(basename "$REF_URL")"
rm -f "$(basename "$REF_URL")"
REF="$(find /root -maxdepth 1 -type d -name 'refdata-gex-*' | head -1)"
[ -f "$REF/genes/genes.gtf.gz" ] && gunzip -kf "$REF/genes/genes.gtf.gz"
FASTA="$REF/fasta/genome.fa"
GTF="$REF/genes/genes.gtf"

aws s3 cp "$BUCKET/whitelist.txt" /root/whitelist.txt
aws s3 cp "$BUCKET/fastq/" /root/fastq --recursive
R1="/root/fastq/${SAMPLE}_S1_L001_R1_001.fastq.gz"   # barcode read
R2="/root/fastq/${SAMPLE}_S1_L001_R2_001.fastq.gz"   # cDNA read

SOLO=(--soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Reverse
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000)

# ---- build the 3 indexes under test ----
mkdir -p idx_old_d1 idx_new_d1 idx_new_d2
/usr/bin/time -v "$RUSTAR_OLD" --runMode genomeGenerate --genomeDir idx_old_d1 \
  --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
  --runThreadN "$THREADS" &> build_old_d1.log
/usr/bin/time -v "$RUSTAR_NEW" --runMode genomeGenerate --genomeDir idx_new_d1 \
  --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
  --runThreadN "$THREADS" &> build_new_d1.log
/usr/bin/time -v "$RUSTAR_NEW" --runMode genomeGenerate --genomeDir idx_new_d2 \
  --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
  --genomeSAsparseD 2 --runThreadN "$THREADS" &> build_new_d2.log

{
  echo "=== index directory sizes ==="
  du -sh idx_old_d1 idx_new_d1 idx_new_d2
  echo "=== SA file sizes (bytes) ==="
  ls -la idx_old_d1/SA idx_new_d1/SA idx_new_d2/SA
} > index_sizes.txt
cat index_sizes.txt

# ---- the comparison: solo align+count to matrix, 3 reps per variant ----
# Warm-cache grouping (see bench_nvme.sh): group each variant's reps together so
# the page cache holds that variant's index across its own reps; rep1 of each
# group is the cold/EBS-adjacent load, reps 2-3 are the true warm compute number.
set +e
for variant in old_d1 new_d1 new_d2; do
  case $variant in
    old_d1) BIN="$RUSTAR_OLD"; IDX=idx_old_d1 ;;
    new_d1) BIN="$RUSTAR_NEW"; IDX=idx_new_d1 ;;
    new_d2) BIN="$RUSTAR_NEW"; IDX=idx_new_d2 ;;
  esac
  for rep in 1 2 3; do
    mkdir -p "/root/out_${variant}_${rep}"
    /usr/bin/time -v "$BIN" --genomeDir "$IDX" --sjdbGTFfile "$GTF" \
      --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
      --soloOutGzip yes --outSAMtype None \
      --outFileNamePrefix "/root/out_${variant}_${rep}/" &> "rustar_${variant}_${rep}.log"
    echo "$variant rep$rep rc=$?"
    grep -H "Uniquely mapped reads %" "/root/out_${variant}_${rep}/Log.final.out" \
      >> mapping_summary.txt 2>/dev/null
  done
done

echo "=== DONE $(date -u) ==="
# trap finish() uploads logs + index_sizes.txt + mapping_summary.txt + powers off here.
