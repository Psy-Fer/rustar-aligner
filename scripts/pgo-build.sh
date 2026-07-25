#!/usr/bin/env bash
# Profile-Guided Optimization build for rustar-aligner.
#
# Two-phase: (1) build an instrumented binary and run it on a representative
# workload to collect branch/layout profiles, (2) rebuild using that profile.
# Output-identical to a normal build; PGO only reorders code / tunes inlining.
#
# Usage: scripts/pgo-build.sh <genomeDir> <readsFastq> [nTrainReads]
# Produces: target/release/rustar-aligner (PGO-optimized), or
# target/<PGO_TARGET>/release/rustar-aligner when PGO_TARGET is set.
#
# Optional env vars (used by the CI release build; empty = host default):
#   PGO_TARGET          — cargo --target triple, e.g. x86_64-unknown-linux-gnu
#   PGO_EXTRA_RUSTFLAGS — extra rustflags composed with the PGO flags on every
#                         build in this script, e.g. "-Ctarget-cpu=x86-64-v3"
#                         (must be applied to BOTH the instrumented build and
#                         the profile-use rebuild, or the training profile
#                         reflects a different instruction-selection profile
#                         than what ships).
set -euo pipefail
GENOME_DIR="${1:?usage: pgo-build.sh <genomeDir> <readsFastq> [nTrainReads]}"
READS="${2:?need a training FASTQ}"
NTRAIN="${3:-1000000}"
PGO_TARGET="${PGO_TARGET:-}"
PGO_EXTRA_RUSTFLAGS="${PGO_EXTRA_RUSTFLAGS:-}"
PGO_DIR="$(pwd)/target/pgo-data"
PROFDATA_BIN="$(find "$(rustc --print sysroot)" -name llvm-profdata | head -1)"
[ -x "$PROFDATA_BIN" ] || { echo "llvm-profdata not found; run: rustup component add llvm-tools-preview"; exit 1; }

TARGET_ARGS=()
BIN_DIR="target/release"
if [ -n "$PGO_TARGET" ]; then
  TARGET_ARGS=(--target "$PGO_TARGET")
  BIN_DIR="target/$PGO_TARGET/release"
fi
BIN="$BIN_DIR/rustar-aligner"

echo "== PGO phase 1: build instrumented binary =="
rm -rf "$PGO_DIR"
RUSTFLAGS="-Cprofile-generate=$PGO_DIR $PGO_EXTRA_RUSTFLAGS" \
  cargo build --release "${TARGET_ARGS[@]+"${TARGET_ARGS[@]}"}"

echo "== PGO phase 1: training run ($NTRAIN reads) =="
TRAIN_OUT="$(mktemp -d)"
"./$BIN" --genomeDir "$GENOME_DIR" --readFilesIn "$READS" \
  --runThreadN 8 --outSAMtype BAM Unsorted --outFileNamePrefix "$TRAIN_OUT/" \
  --readMapNumber "$NTRAIN" >/dev/null 2>&1
rm -rf "$TRAIN_OUT"

echo "== PGO: merge profiles =="
"$PROFDATA_BIN" merge -o "$PGO_DIR/merged.profdata" "$PGO_DIR"

echo "== PGO phase 2: rebuild with profile =="
RUSTFLAGS="-Cprofile-use=$PGO_DIR/merged.profdata -Cllvm-args=-pgo-warn-missing-function $PGO_EXTRA_RUSTFLAGS" \
  cargo build --release "${TARGET_ARGS[@]+"${TARGET_ARGS[@]}"}"
echo "== done: $BIN is PGO-optimized =="
