#!/bin/bash
# ---------------------------------------------------------------------------
# Profile the SOLO run specifically (not the plain aligner) to locate the
# rustar-vs-STARsolo runtime gap. Same workload as bench_3way's rustar arm:
# mouse 5' GEM-X, 10M reads, GeneFull/Gene/SJ/Velocyto, EmptyDrops, GRCm39.
#
# Deliverables:
#   - baseline: 3 timed solo reps (warm) — confirms the ~120s figure
#   - perf record -g on a warm solo run: flat self-time + call graph, so we see
#     how much is aligner core vs solo classify/count vs collation/gzip/EmptyDrops
#   - phase split from the rustar log's own INFO timestamps: parallel align
#     phase (start -> "Processed 10000000 reads") vs the serial post-align tail
#     (matrix collation + EmptyDrops + gzip matrix write)
# ---------------------------------------------------------------------------
set -uxo pipefail

BUCKET="s3://rustar-bench"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:418696582915:rustar-bench"
THREADS=10
SAMPLE="5k_Mouse_PBMCs_5p_gem-x_GEX"
REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
OVERHANG=89
NEW_SRC_S3="$BUCKET/src/rustar-new-src.tar.gz"

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
RESULTS="$BUCKET/profile/solo-$TS"
exec > >(tee /var/log/profile_solo.log) 2>&1

finish() {
  set +e
  aws s3 cp /var/log/profile_solo.log "$RESULTS/" 2>/dev/null
  for f in /root/*.log /root/*.txt; do aws s3 cp "$f" "$RESULTS/" 2>/dev/null; done
  if [ -n "$SNS_TOPIC_ARN" ]; then
    SUMMARY=$(for f in /root/solo_*.log; do
      [ -f "$f" ] || continue
      w=$(grep -m1 'Elapsed (wall clock)' "$f" | sed 's/.*: //')
      echo "$(basename "$f" .log): wall=${w:-FAILED}"
    done)
    aws sns publish --topic-arn "$SNS_TOPIC_ARN" \
      --subject "rustar SOLO profile done ($TS)" \
      --message "Results: $RESULTS/

$SUMMARY

Pull: aws s3 cp $RESULTS/ ./results --recursive" 2>/dev/null
  fi
  shutdown -h now
}
trap finish EXIT
shutdown -h +180

echo "=== START $TS ==="
set -e
dnf -y install git gcc make tar gzip time wget unzip which perf
curl -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"

# rustar from current working tree, with debug syms + frame pointers for perf.
aws s3 cp "$NEW_SRC_S3" new_src.tar.gz
mkdir -p new_src && tar xzf new_src.tar.gz -C new_src --strip-components=1
( cd new_src && RUSTFLAGS="-C force-frame-pointers=yes" \
  CARGO_PROFILE_RELEASE_STRIP=false CARGO_PROFILE_RELEASE_DEBUG=2 \
  cargo build --release ) &> build_rustar.log
RUSTAR=/root/new_src/target/release/rustar-aligner

# reference (for the GTF), whitelist, fastq
wget -q "$REF_URL" ; tar xf "$(basename "$REF_URL")" ; rm -f "$(basename "$REF_URL")"
REF="$(find /root -maxdepth 1 -type d -name 'refdata-gex-*' | head -1)"
[ -f "$REF/genes/genes.gtf.gz" ] && gunzip -kf "$REF/genes/genes.gtf.gz"
FASTA="$REF/fasta/genome.fa"
GTF="$REF/genes/genes.gtf"
aws s3 cp "$BUCKET/whitelist.txt" /root/whitelist.txt
aws s3 cp "$BUCKET/fastq/" /root/fastq --recursive
R1="/root/fastq/${SAMPLE}_S1_L001_R1_001.fastq.gz"
R2="/root/fastq/${SAMPLE}_S1_L001_R2_001.fastq.gz"

# rustar index built fresh from THIS binary (with GTF), then cached for re-runs.
/usr/bin/time -v "$RUSTAR" --runMode genomeGenerate --genomeDir rust_idx \
  --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OVERHANG" \
  --runThreadN 16 &> build_index.log
aws s3 sync rust_idx "$BUCKET/idx/rust_idx_current/" --only-show-errors || true

SOLO=(--soloType CB_UMI_Simple --soloCBwhitelist /root/whitelist.txt
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloFeatures Gene GeneFull SJ Velocyto --soloStrand Reverse
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30
  --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 100 0.01 20000 0.01 10000)

set +e
# warm the index page cache (untimed)
mkdir -p /root/warm/
"$RUSTAR" --genomeDir rust_idx --sjdbGTFfile "$GTF" \
  --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
  --soloOutGzip yes --outSAMtype None --outFileNamePrefix /root/warm/ > warm.stdout.log 2>&1

# 3 timed reps — the log INFO timestamps give the align vs post-align phase split.
for r in 1 2 3; do
  mkdir -p "/root/out_${r}/"
  /usr/bin/time -v "$RUSTAR" --genomeDir rust_idx --sjdbGTFfile "$GTF" \
    --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
    --soloOutGzip yes --outSAMtype None --outFileNamePrefix "/root/out_${r}/" &> "solo_${r}.log"
  echo "solo rep$r rc=$?"
done

# perf record on one warm solo run (full run incl. collation tail).
mkdir -p /root/perfout/
perf record -g -F 999 -o /root/perf.data -- \
  "$RUSTAR" --genomeDir rust_idx --sjdbGTFfile "$GTF" \
  --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
  --soloOutGzip yes --outSAMtype None --outFileNamePrefix /root/perfout/ 2> perf_record.txt
perf report -i /root/perf.data --stdio --no-children -g none 2>/dev/null | head -70 > perf_self.txt
perf report -i /root/perf.data --stdio 2>/dev/null | head -170 > perf_callgraph.txt
perf stat -e cycles,instructions,branches,branch-misses,cache-references,cache-misses \
  -o perf_stat.txt -- "$RUSTAR" --genomeDir rust_idx --sjdbGTFfile "$GTF" \
  --readFilesIn "$R2" "$R1" "${SOLO[@]}" --runThreadN "$THREADS" \
  --soloOutGzip yes --outSAMtype None --outFileNamePrefix /root/perfstat_out/ 2>&1 | tee -a perf_stat.txt

# phase-split helper: pull the timeline markers from rep2's own log.
{
  echo "=== rustar solo phase timeline (rep2 INFO timestamps) ==="
  grep -E "Aligning reads|Processed .*reads|Writing splice junction|STARsolo|Collating|matrix|Log.final|Alignment complete|Solo" solo_2.log | head -60
} > phase_timeline.txt

echo "=== DONE $(date -u) ==="
