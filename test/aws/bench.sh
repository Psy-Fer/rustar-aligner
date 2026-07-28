#!/bin/bash
# ---------------------------------------------------------------------------
# rustar vs STARsolo vs CellRanger — fair native-x86_64 wall-time + Max-RSS
# benchmark. Runs as EC2 user-data on a self-terminating Spot instance.
#
# Inputs are pulled from S3 (you upload them once via setup_aws.sh notes):
#   $BUCKET/cellranger-10.0.0.tar.gz        your CellRanger tarball
#   $BUCKET/refdata-gex-*.tar.gz            your CellRanger mouse reference
#   $BUCKET/whitelist.txt                   the 5' GEM-X CB whitelist
#   $BUCKET/fastq/<SAMPLE>_S1_L001_R{1,2}_001.fastq.gz
#
# STAR + rustar indexes are BUILT here from the reference's genome.fa/genes.gtf
# (the same sequence/annotation CellRanger uses) so all three are apples-to-apples.
# Results (all logs) are uploaded to $BUCKET/results/<ts>/ and the box self-
# terminates. A 4 h dead-man's switch guarantees termination even if this hangs.
# ---------------------------------------------------------------------------
set -uxo pipefail

############################ CONFIG — EDIT THESE ############################
BUCKET="s3://rustar-bench"          # your bucket (no trailing slash)
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"  # summary email
THREADS=10                                    # fixed for ALL tools (fairness)
SAMPLE="5k_Mouse_PBMCs_5p_gem-x_GEX"          # FASTQ prefix before _S1_L001_*
CR_TGZ="cellranger-10.0.0.tar.xz"             # staged in S3 (xz; tar auto-detects)
REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz" # public 10x CDN — wget on instance, not uploaded (9.6 GB)
OVERHANG=89                                   # R2 read length 90 - 1
STAR_VER="2.7.11b"
RUSTAR_BRANCH="perf-hashers"   # optimized aligner (P1+P2): FxHash + mimalloc purge-off
############################################################################

cd /root
export HOME=/root  # cloud-init runs user-data with HOME unset; `set -u` + $HOME would abort
# AL2023 mounts /tmp as tmpfs (RAM, ~32 GB). rustar's caps-sa external-sort temp uses
# std::env::temp_dir() (=/tmp), so it fills RAM-backed /tmp regardless of root disk size.
# Point all temp at the big root volume instead.
mkdir -p /root/rtmp && export TMPDIR=/root/rtmp
TS=$(date -u +%Y%m%dT%H%M%SZ)
RESULTS="$BUCKET/results/$TS"
exec > >(tee /var/log/bench.log) 2>&1

# Always upload whatever we have, email a summary, and power off — even on error.
finish() {
  set +e
  aws s3 cp /var/log/bench.log "$RESULTS/" 2>/dev/null
  for f in /root/*.log; do aws s3 cp "$f" "$RESULTS/" 2>/dev/null; done
  if [ -n "$SNS_TOPIC_ARN" ]; then
    # one-line wall time per tool/rep, scraped straight from the time -v logs
    SUMMARY=$(for f in /root/rustar_*.log /root/star_*.log /root/cr_*.log; do
      [ -f "$f" ] || continue
      w=$(grep -m1 'Elapsed (wall clock)' "$f" | sed 's/.*: //')
      echo "$(basename "$f" .log): ${w:-FAILED}"
    done)
    aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
      --subject "rustar benchmark done ($TS)" \
      --message "Results: $RESULTS/

Wall time (Elapsed) per run:
$SUMMARY

Pull + summarize:
  aws s3 cp $RESULTS/ ./results --recursive
  python3 scrape_results.py ./results" 2>/dev/null
  fi
  shutdown -h now
}
trap finish EXIT
shutdown -h +240        # dead-man's switch: terminate after 4 h no matter what

echo "=== START $TS  (sample=$SAMPLE threads=$THREADS) ==="

# ---- dependencies (fail fast here) ----
set -e
dnf -y install git gcc make tar gzip time wget unzip which
curl -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

# rustar (built from source = exact version under test)
git clone -b "$RUSTAR_BRANCH" https://github.com/iandriver/rustar-aligner
( cd rustar-aligner && cargo build --release )
RUSTAR=/root/rustar-aligner/target/release/rustar-aligner

# STARsolo (official static x86_64 build = the version you validated against)
wget -q "https://github.com/alexdobin/STAR/releases/download/${STAR_VER}/STAR_${STAR_VER}.zip"
unzip -q "STAR_${STAR_VER}.zip"
STAR="/root/STAR_${STAR_VER}/Linux_x86_64_static/STAR"; chmod +x "$STAR"

# ---- inputs from S3 ----
aws s3 cp "$BUCKET/$CR_TGZ" . ; tar xf "$CR_TGZ"
CELLRANGER="$(ls -d /root/cellranger-*/cellranger | head -1)"
wget -q "$REF_URL" ; tar xf "$(basename "$REF_URL")" ; rm -f "$(basename "$REF_URL")"
REF="$(find /root -maxdepth 1 -type d -name 'refdata-gex-*' | head -1)"
# CellRanger 2024-A refs ship genes/genes.gtf GZIPPED; STAR + rustar need it plain
# (CellRanger reads the .gz itself, so leaving the .gz in place is fine).
[ -f "$REF/genes/genes.gtf.gz" ] && gunzip -kf "$REF/genes/genes.gtf.gz"
aws s3 cp "$BUCKET/whitelist.txt" /root/whitelist.txt
aws s3 cp "$BUCKET/fastq/" /root/fastq --recursive

FASTA="$REF/fasta/genome.fa"
GTF="$REF/genes/genes.gtf"
R1="/root/fastq/${SAMPLE}_S1_L001_R1_001.fastq.gz"   # barcode read
R2="/root/fastq/${SAMPLE}_S1_L001_R2_001.fastq.gz"   # cDNA read

# ---- STAR + rustar indexes: reuse cached copies in S3 if present, else build+cache.
# (Caching skips the ~75 min rebuild on re-runs; build timings come from the first run.)
mkdir -p star_idx rust_idx
echo "=== disk before builds ==="; df -h /
if aws s3 ls "$BUCKET/idx/rust_idx/SAindex" >/dev/null 2>&1; then
  echo "=== using cached indexes from S3 ==="
  aws s3 sync "$BUCKET/idx/star_idx/" star_idx --only-show-errors
  aws s3 sync "$BUCKET/idx/rust_idx/" rust_idx --only-show-errors
else
  /usr/bin/time -v "$STAR" --runMode genomeGenerate --genomeDir star_idx \
    --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
    --runThreadN 16 &> build_star.log
  echo "=== disk after STAR build (before rustar caps-sa) ==="; df -h /
  /usr/bin/time -v "$RUSTAR" --runMode genomeGenerate --genomeDir rust_idx \
    --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
    --runThreadN 16 &> build_rustar.log
  aws s3 sync star_idx "$BUCKET/idx/star_idx/" --only-show-errors || true
  aws s3 sync rust_idx "$BUCKET/idx/rust_idx/" --only-show-errors || true
fi

# Shared CellRanger-matching solo flags (rustar mirrors STAR's flag names).
SOLO=(--soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Reverse
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000)
  #              ^ STAR requires exactly 10 EmptyDrops_CR params (last = simN)

# ---- the comparison: align+count to matrix (no BAM for any tool) ----
# We measure WARM compute (no drop_caches): rustar mmaps its index and accesses the
# SA randomly, so a cold index on EBS is I/O-bound, not compute-bound (the same
# mmap-on-slow-storage effect as Rosetta/virtiofs). With 64 GB RAM the ~26 GB index
# page-caches, so we group each tool's 3 reps together (cross-tool would thrash the
# cache: rustar 26 GB + STAR 30 GB + CR-ref 30 GB > RAM). rep1 is cold (loads from
# EBS), reps 2-3 warm — the warm median is each tool's true compute speed.
set +e   # one tool failing must not abort the others
for rep in 1 2 3; do
  mkdir -p "/root/out_rustar_${rep}"   # rustar (like STAR) needs the prefix dir to exist
  /usr/bin/time -v "$RUSTAR" --genomeDir rust_idx --sjdbGTFfile "$GTF" \
    --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
    --soloOutGzip yes --outSAMtype None \
    --outFileNamePrefix "/root/out_rustar_${rep}/" &> "rustar_${rep}.log"
  echo "rustar rep$rep rc=$?"
done
for rep in 1 2 3; do
  mkdir -p "/root/out_star_${rep}"
  /usr/bin/time -v "$STAR" --genomeDir star_idx \
    --readFilesIn "$R2" "$R1" --readFilesCommand zcat "${SOLO[@]}" \
    --runThreadN "$THREADS" --outSAMtype None \
    --outFileNamePrefix "/root/out_star_${rep}/" &> "star_${rep}.log"
  echo "star rep$rep rc=$?"
done
for rep in 1 2 3; do
  rm -rf "/root/cr_${rep}"
  /usr/bin/time -v "$CELLRANGER" count --id="cr_${rep}" \
    --transcriptome="$REF" --fastqs=/root/fastq --sample="$SAMPLE" \
    --localcores="$THREADS" --localmem=56 --create-bam=false --nosecondary \
    &> "cr_${rep}.log"
  echo "cellranger rep$rep rc=$?"
done

echo "=== DONE $(date -u) ==="
# trap finish() uploads logs + powers off here.
