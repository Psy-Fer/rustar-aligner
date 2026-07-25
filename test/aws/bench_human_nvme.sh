#!/bin/bash
# ---------------------------------------------------------------------------
# HUMAN GRCh38 benchmark: rustar vs STARsolo vs CellRanger — fair native-x86_64
# wall-time + Max-RSS, PLUS index-build peak-RAM (the GRCh38 genomeGenerate is
# the memory-heavy step we want to measure). EC2 user-data, self-terminating,
# all I/O on local NVMe (no gp3 EBS index-reload variance).
#
# Dataset: 10x Genomics OFFICIAL 5k Human PBMC 3' v3 (staged by prep_10x.sh after
# a whitelist-match gate). 3' v3: R1=28bp (CB16+UMI12), R2=91bp cDNA, Strand Forward,
# inclusion list = full 6.8M 3M-february-2018 (in S3 whitelists/).
#
# Index-build peak RAM = the "Maximum resident set size" line from /usr/bin/time -v
# wrapping each genomeGenerate (build_star.log / build_rustar.log), surfaced in the
# SNS summary. STAR GRCh38 build ~32-40 GB; rustar caps-sa ~15-18 GB; 64 GB fits both.
# ---------------------------------------------------------------------------
set -uxo pipefail

############################ CONFIG ############################
BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
THREADS=10
SAMPLE="human_pbmc5k_3pv3"                     # FASTQ prefix (staged by prep_10x.sh)
FASTQ_S3="$BUCKET/fastq_${SAMPLE}_10M"         # 10M-read subsample
CR_TGZ="cellranger-10.0.0.tar.xz"
REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"  # public 10x CDN (~16 GB)
WL_S3="$BUCKET/whitelists/3M-february-2018.txt.gz"
OVERHANG=90                                    # R2 length 91 - 1
STAR_GENRAM=48000000000                        # 48 GB cap for STAR genomeGenerate (GRCh38 needs >31 GB default)
STAR_VER="2.7.11b"
RUSTAR_BRANCH="perf-hashers"
###############################################################

# --- mount instance-store NVMe at /root (all build/index/ref/fastq/output I/O local)
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
RESULTS="$BUCKET/results/$TS"
exec > >(tee /var/log/bench.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/bench.log "$RESULTS/" 2>/dev/null
  for f in /root/*.log; do aws s3 cp "$f" "$RESULTS/" 2>/dev/null; done
  # index-build peak RAM (kbytes -> GB) from the time -v build logs
  bpk() { grep -m1 'Maximum resident' "$1" 2>/dev/null | grep -oE '[0-9]+' | awk '{printf "%.1f GB", $1/1048576}'; }
  if [ -n "$SNS_TOPIC_ARN" ]; then
    WALL=$(for f in /root/rustar_*.log /root/star_*.log /root/cr_*.log; do
      [ -f "$f" ] || continue
      w=$(grep -m1 'Elapsed (wall clock)' "$f" | sed 's/.*: //')
      echo "$(basename "$f" .log): ${w:-FAILED}"
    done)
    aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
      --subject "rustar HUMAN benchmark done ($TS)" \
      --message "Results: $RESULTS/

=== GRCh38 index-build PEAK RAM ===
STAR genomeGenerate : $(bpk /root/build_star.log)
rustar genomeGenerate: $(bpk /root/build_rustar.log)

=== align+count wall time per run ===
$WALL

Pull + summarize:
  aws s3 cp $RESULTS/ ./results --recursive
  python3 scrape_results.py ./results" 2>/dev/null
  fi
  shutdown -h now
}
trap finish EXIT
shutdown -h +300        # 5 h dead-man's switch (human builds + CR are slower)

echo "=== START $TS  (HUMAN sample=$SAMPLE threads=$THREADS) ==="

set -e
dnf -y install git gcc make tar gzip time wget unzip which
curl -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

git clone -b "$RUSTAR_BRANCH" https://github.com/iandriver/rustar-aligner
( cd rustar-aligner && cargo build --release )
RUSTAR=/root/rustar-aligner/target/release/rustar-aligner

wget -q "https://github.com/alexdobin/STAR/releases/download/${STAR_VER}/STAR_${STAR_VER}.zip"
unzip -q "STAR_${STAR_VER}.zip"
STAR="/root/STAR_${STAR_VER}/Linux_x86_64_static/STAR"; chmod +x "$STAR"

# ---- inputs from S3 + CDN ----
aws s3 cp "$BUCKET/$CR_TGZ" . ; tar xf "$CR_TGZ"
CELLRANGER="$(ls -d /root/cellranger-*/cellranger | head -1)"
wget -q "$REF_URL" ; tar xf "$(basename "$REF_URL")" ; rm -f "$(basename "$REF_URL")"
REF="$(find /root -maxdepth 1 -type d -name 'refdata-gex-*' | head -1)"
[ -f "$REF/genes/genes.gtf.gz" ] && gunzip -kf "$REF/genes/genes.gtf.gz"
aws s3 cp "$WL_S3" wl.txt.gz ; gzip -df wl.txt.gz ; mv wl.txt /root/whitelist.txt
aws s3 cp "$FASTQ_S3/" /root/fastq --recursive

FASTA="$REF/fasta/genome.fa"
GTF="$REF/genes/genes.gtf"
R1="/root/fastq/${SAMPLE}_S1_L001_R1_001.fastq.gz"   # barcode (28bp)
R2="/root/fastq/${SAMPLE}_S1_L001_R2_001.fastq.gz"   # cDNA (91bp)

# ---- STAR + rustar indexes: build from GRCh38 (this is the index-build peak-RAM measurement).
# Cache to idx_human/ so re-runs skip the rebuild; build timings/peaks come from the first run.
mkdir -p star_idx rust_idx
echo "=== disk before builds ==="; df -h /root
if aws s3 ls "$BUCKET/idx_human/rust_idx/SAindex" >/dev/null 2>&1; then
  echo "=== using cached GRCh38 indexes from S3 (build peaks already measured) ==="
  aws s3 sync "$BUCKET/idx_human/star_idx/" star_idx --only-show-errors
  aws s3 sync "$BUCKET/idx_human/rust_idx/" rust_idx --only-show-errors
else
  echo "=== building STAR GRCh38 index (peak RAM via time -v) ==="
  /usr/bin/time -v "$STAR" --runMode genomeGenerate --genomeDir star_idx \
    --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
    --limitGenomeGenerateRAM "$STAR_GENRAM" --runThreadN 16 &> build_star.log
  echo "=== disk after STAR build ==="; df -h /root
  echo "=== building rustar GRCh38 index (peak RAM via time -v) ==="
  /usr/bin/time -v "$RUSTAR" --runMode genomeGenerate --genomeDir rust_idx \
    --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
    --runThreadN 16 &> build_rustar.log
  aws s3 sync star_idx "$BUCKET/idx_human/star_idx/" --only-show-errors || true
  aws s3 sync rust_idx "$BUCKET/idx_human/rust_idx/" --only-show-errors || true
fi
echo "=== STAR build peak: $(grep -m1 'Maximum resident' build_star.log)"
echo "=== rustar build peak: $(grep -m1 'Maximum resident' build_rustar.log)"

# 3' v3 solo flags (Strand FORWARD — differs from the mouse 5' run which used Reverse).
SOLO=(--soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Forward
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000)

set +e
for rep in 1 2 3; do
  mkdir -p "/root/out_rustar_${rep}"
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
