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
#   PGO_EXTRA_RUSTFLAGS — extra rustflags for the SHIPPED (profile-use) build,
#                         e.g. "-Ctarget-cpu=x86-64-v4".
#   PGO_TRAIN_EXTRA_RUSTFLAGS — extra rustflags for the INSTRUMENTED build that
#                         is executed to collect the profile. Defaults to
#                         PGO_EXTRA_RUSTFLAGS. Set this to a runner-executable
#                         ISA when the ship target uses instructions the build
#                         machine may lack: e.g. ship at x86-64-v4 (AVX-512) but
#                         train at x86-64-v3 (AVX2), since GitHub runners do not
#                         reliably have AVX-512 and would SIGILL (exit 132) on
#                         the training run. PGO profiles are behavioral (block/
#                         branch frequencies) and transfer across -Ctarget-cpu
#                         levels, so the shipped binary is still v4-optimized.
set -euo pipefail
GENOME_DIR="${1:?usage: pgo-build.sh <genomeDir> <readsFastq> [nTrainReads]}"
READS="${2:?need a training FASTQ}"
NTRAIN="${3:-1000000}"
PGO_TARGET="${PGO_TARGET:-}"
PGO_EXTRA_RUSTFLAGS="${PGO_EXTRA_RUSTFLAGS:-}"
PGO_TRAIN_EXTRA_RUSTFLAGS="${PGO_TRAIN_EXTRA_RUSTFLAGS:-$PGO_EXTRA_RUSTFLAGS}"
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

echo "== PGO phase 1: build instrumented binary (train flags: ${PGO_TRAIN_EXTRA_RUSTFLAGS:-none}) =="
rm -rf "$PGO_DIR"
RUSTFLAGS="-Cprofile-generate=$PGO_DIR $PGO_TRAIN_EXTRA_RUSTFLAGS" \
  cargo build --release "${TARGET_ARGS[@]+"${TARGET_ARGS[@]}"}"

echo "== PGO phase 1: training run ($NTRAIN reads) =="
TRAIN_OUT="$(mktemp -d)"
"./$BIN" --genomeDir "$GENOME_DIR" --readFilesIn "$READS" \
  --runThreadN 8 --outSAMtype BAM Unsorted --outFileNamePrefix "$TRAIN_OUT/" \
  --readMapNumber "$NTRAIN" >/dev/null 2>&1
rm -rf "$TRAIN_OUT"

echo "== PGO: merge profiles =="
"$PROFDATA_BIN" merge -o "$PGO_DIR/merged.profdata" "$PGO_DIR"

echo "== PGO phase 2: rebuild with profile (ship flags: ${PGO_EXTRA_RUSTFLAGS:-none}) =="
RUSTFLAGS="-Cprofile-use=$PGO_DIR/merged.profdata -Cllvm-args=-pgo-warn-missing-function $PGO_EXTRA_RUSTFLAGS" \
  cargo build --release "${TARGET_ARGS[@]+"${TARGET_ARGS[@]}"}"
echo "== done: $BIN is PGO-optimized =="
