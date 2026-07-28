#!/bin/bash
# ---------------------------------------------------------------------------
# rustar (current code, dense + sparse-D2) vs STARsolo vs CellRanger —
# fresh, single-instance, apples-to-apples solo benchmark. Wall time + Max RSS
# for all four, same mouse 5' GEM-X workload, same threads, NVMe, warm cache.
#
# Differs from bench_nvme.sh:
#   - rustar is built from the CURRENT working tree (S3 src tarball), not a git
#     branch, so it reflects this session's changes.
#   - adds a 4th variant: rustar with --genomeSAsparseD 2 (the sparse-SA memory
#     win) so the memory comparison is current, not stitched from an old run.
#   - captures per-tool cell counts alongside wall/RSS.
# ---------------------------------------------------------------------------
set -uxo pipefail

############################ CONFIG ############################
BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
THREADS=10                                    # fixed for ALL tools (fairness)
SAMPLE="5k_Mouse_PBMCs_5p_gem-x_GEX"
CR_TGZ="cellranger-10.0.0.tar.xz"
REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
OVERHANG=89
STAR_VER="2.7.11b"
NEW_SRC_S3="$BUCKET/src/rustar-new-src.tar.gz"   # current working tree
###############################################################

# --- Mount instance-store NVMe at /root (see bench_nvme.sh rationale).
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
RESULTS="$BUCKET/results/3way-$TS"
exec > >(tee /var/log/bench.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/bench.log "$RESULTS/" 2>/dev/null
  for f in /root/*.log /root/*.csv /root/*.txt; do aws s3 cp "$f" "$RESULTS/" 2>/dev/null; done
  if [ -n "$SNS_TOPIC_ARN" ]; then
    SUMMARY=$(for f in /root/rustar_d1_*.log /root/rustar_d2_*.log /root/star_*.log /root/cr_*.log; do
      [ -f "$f" ] || continue
      w=$(grep -m1 'Elapsed (wall clock)' "$f" | sed 's/.*: //')
      r=$(grep -m1 'Maximum resident set size' "$f" | sed 's/.*: //')
      echo "$(basename "$f" .log): wall=${w:-FAILED} rssKB=${r:-?}"
    done)
    aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
      --subject "rustar 3-way benchmark done ($TS)" \
      --message "Results: $RESULTS/

$SUMMARY

Pull: aws s3 cp $RESULTS/ ./results --recursive" 2>/dev/null
  fi
  shutdown -h now
}
trap finish EXIT
shutdown -h +300        # 5 h dead-man's switch (4 variants x 3 reps + CR is slow)

echo "=== START $TS  (sample=$SAMPLE threads=$THREADS) ==="

set -e
dnf -y install git gcc make tar gzip time wget unzip which
curl -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

# rustar from the current working-tree tarball
aws s3 cp "$NEW_SRC_S3" new_src.tar.gz
mkdir -p new_src && tar xzf new_src.tar.gz -C new_src --strip-components=1
( cd new_src && cargo build --release ) &> build_rustar_bin.log
RUSTAR=/root/new_src/target/release/rustar-aligner

# STARsolo static x86_64
wget -q "https://github.com/alexdobin/STAR/releases/download/${STAR_VER}/STAR_${STAR_VER}.zip"
unzip -q "STAR_${STAR_VER}.zip"
STAR="/root/STAR_${STAR_VER}/Linux_x86_64_static/STAR"; chmod +x "$STAR"

# inputs
aws s3 cp "$BUCKET/$CR_TGZ" . ; tar xf "$CR_TGZ"
CELLRANGER="$(ls -d /root/cellranger-*/cellranger | head -1)"
wget -q "$REF_URL" ; tar xf "$(basename "$REF_URL")" ; rm -f "$(basename "$REF_URL")"
REF="$(find /root -maxdepth 1 -type d -name 'refdata-gex-*' | head -1)"
[ -f "$REF/genes/genes.gtf.gz" ] && gunzip -kf "$REF/genes/genes.gtf.gz"
aws s3 cp "$BUCKET/whitelist.txt" /root/whitelist.txt
aws s3 cp "$BUCKET/fastq/" /root/fastq --recursive

FASTA="$REF/fasta/genome.fa"
GTF="$REF/genes/genes.gtf"
R1="/root/fastq/${SAMPLE}_S1_L001_R1_001.fastq.gz"
R2="/root/fastq/${SAMPLE}_S1_L001_R2_001.fastq.gz"

# ---- indexes ----
# STAR: reuse cached (STAR unchanged). rustar D1 + D2: build FRESH from the
# current binary so the index (and its size) reflects this session's code.
mkdir -p star_idx
if aws s3 ls "$BUCKET/idx/star_idx/SAindex" >/dev/null 2>&1; then
  echo "=== using cached STAR index ==="
  aws s3 sync "$BUCKET/idx/star_idx/" star_idx --only-show-errors
else
  /usr/bin/time -v "$STAR" --runMode genomeGenerate --genomeDir star_idx \
    --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
    --runThreadN 16 &> build_star.log
  aws s3 sync star_idx "$BUCKET/idx/star_idx/" --only-show-errors || true
fi

/usr/bin/time -v "$RUSTAR" --runMode genomeGenerate --genomeDir rust_idx_d1 \
  --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
  --runThreadN 16 &> build_rustar_d1.log
/usr/bin/time -v "$RUSTAR" --runMode genomeGenerate --genomeDir rust_idx_d2 \
  --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
  --genomeSAsparseD 2 --runThreadN 16 &> build_rustar_d2.log

{
  echo "=== index directory sizes ==="
  du -sh star_idx rust_idx_d1 rust_idx_d2
  echo "=== SA file sizes (bytes) ==="
  ls -la star_idx/SA rust_idx_d1/SA rust_idx_d2/SA
} > index_sizes.txt
cat index_sizes.txt

# Shared CellRanger-matching solo flags.
SOLO=(--soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Reverse
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000)

# ---- runs: group each tool's reps for warm-cache fairness ----
set +e
for rep in 1 2 3; do
  mkdir -p "/root/out_rustar_d1_${rep}"
  /usr/bin/time -v "$RUSTAR" --genomeDir rust_idx_d1 --sjdbGTFfile "$GTF" \
    --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
    --soloOutGzip yes --outSAMtype None \
    --outFileNamePrefix "/root/out_rustar_d1_${rep}/" &> "rustar_d1_${rep}.log"
  echo "rustar_d1 rep$rep rc=$?"
done
for rep in 1 2 3; do
  mkdir -p "/root/out_rustar_d2_${rep}"
  /usr/bin/time -v "$RUSTAR" --genomeDir rust_idx_d2 --sjdbGTFfile "$GTF" \
    --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
    --soloOutGzip yes --outSAMtype None \
    --outFileNamePrefix "/root/out_rustar_d2_${rep}/" &> "rustar_d2_${rep}.log"
  echo "rustar_d2 rep$rep rc=$?"
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

# ---- collect comparable metrics (best-effort; failures are non-fatal) ----
{
  echo "=== cell counts / mapping (per tool, rep 2) ==="
  echo "-- rustar_d1 Gene Summary --"; cat /root/out_rustar_d1_2/Solo.out/Gene/Summary.csv 2>/dev/null
  echo "-- rustar_d2 Gene Summary --"; cat /root/out_rustar_d2_2/Solo.out/Gene/Summary.csv 2>/dev/null
  echo "-- star Gene Summary --";      cat /root/out_star_2/Solo.out/Gene/Summary.csv 2>/dev/null
  echo "-- cellranger metrics --";     cat /root/cr_2/outs/metrics_summary.csv 2>/dev/null
} > metrics_3way.txt
cat metrics_3way.txt

echo "=== DONE $(date -u) ==="
