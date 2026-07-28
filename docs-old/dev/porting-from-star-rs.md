# Porting features from STAR-rs into rustar-aligner

This note records the conventions for bringing features over from
[STAR-rs](https://github.com/BenjaminDEMAILLE/STAR-rs) (a clean-room Rust re-port of
STAR 2.7.11b, validated byte-identical against native STAR) into rustar-aligner.
It is the shared reference for the themed, dependency-stacked porting PRs.

## Why

rustar-aligner and STAR-rs are two independent Rust reimplementations of STAR.
rustar-aligner already covers single/paired-end alignment, splice junctions,
two-pass, 4-tier chimeric, GeneCounts, TranscriptomeSAM and sorted BAM. STAR-rs
additionally implements STARsolo, STARlong, WASP, genome transforms, signal tracks,
adapter clipping, PE mate-overlap, SuperTranscriptome, liftOver and more. These
notes cover moving that surface over, one theme per PR.

## Licensing

STAR-rs is dual-licensed `MIT OR Apache-2.0` (Copyright (c) 2026 Benjamin Demaille).
Ported code is taken under the **MIT** option, compatible with this repo's MIT
license. `LICENSES/STAR_RS_LICENSE` preserves the STAR-rs MIT notice for attribution.

**Never copy `vendor/cellranger-upstream` code from STAR-rs.** It is 10x Genomics
source-available (non-OSI) and kept in STAR-rs as algorithm *reference only*. When
porting STARsolo, reimplement the algorithms from the STAR C++ source, as STAR-rs did.

## Type / module mapping

STAR-rs is a 10-crate workspace (`star_core`, `star_index`, `star_seed`,
`star_align`, `star_cli`, ...); rustar-aligner is a single crate (`rustar_aligner`)
with a flat module tree. A ported file's `use star_*::...` header is repointed to
`crate::...`, and the core types are adapted field-by-field:

| STAR-rs | rustar-aligner (`src/align/transcript.rs`) |
| --- | --- |
| `star_align::star_model::Transcript` | `crate::align::transcript::Transcript` |
| `Exon { r, g, len, i_frag }` | `Exon { read_start, genome_start, genome_end, i_frag }` (start+end, not start+len) |
| `tr.str_: u8` (0 fwd / 1 rev) | `tr.is_reverse: bool` |
| `tr.chr` | `tr.chr_idx` |
| `tr.n_exons()` | `tr.exons.len()` |
| `star_core::Cigar` (derived from exons late) | `Vec<noodles ...::cigar::Op>` (materialized in `Transcript`) |
| `star_index::load::StarGenome` (unified) | split `genome::Genome` + `index::GenomeIndex` + `index::suffix_array::SuffixArray` |

Base encoding is identical in both: `A=0, C=1, G=2, T=3, N=4`, reverse complement in
the second half of the genome array.

## Determinism

rustar-aligner currently uses `rand` (`StdRng`, seeded by `--runRNGseed`) for
multimapper-order / tie-break randomness. STAR-rs instead uses a hand-rolled,
per-read-seeded `splitmix64` shuffle so output is **byte-identical regardless of
thread count** (a stronger correctness property than STAR's per-thread `mt19937`).

**Decision for the porting effort:** adopt STAR-rs's thread-invariant model. Ported
code should not introduce new `rand` usage; the migration of the existing tie-break
path lands in its own PR before STARsolo (whose UMI dedup depends on deterministic
order). This is a behavioural change from STAR's per-thread RNG and is flagged here
for team discussion.

## Validation

STAR-rs's `DIVERGENCES.md` (entries D1-D24) enumerates every intentional difference
from the native-STAR oracle, each locked by a named test. Treat it as the acceptance
spec: a ported feature should match native STAR except where DIVERGENCES.md documents
an intentional difference.

rustar-aligner already has a differential harness under `test/` (`compare_sam.py`,
`compare_pe.py`, `assess_faithfulness.py`, `compare_junctions.py`,
`compare_chimeric.py`, driven by `test/ci.sh`) that runs native STAR and diffs
outputs. Use it to validate each port against `STAR` on the yeast test set; it is a
local/manual harness and is not part of the GitHub `alls-green` CI gate.

Every PR must still pass the standard gate: `cargo build`, `cargo test`,
`cargo clippy --all-targets -- -D warnings`, `cargo fmt --check`, and MSRV
`cargo +1.89 check`.
