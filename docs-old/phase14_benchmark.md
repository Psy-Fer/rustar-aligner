[← Back to ROADMAP](../ROADMAP.md) · [Phase 14](phase14_starsolo.md)

# Phase 14 Benchmark: CellRanger vs STARsolo vs rustar-aligner

Runtime + output-stats comparison of the three single-cell quantifiers on a real
10x mouse dataset, run in one consistent Linux/x86_64 environment.

## Setup

- **Reference**: CellRanger mouse `refdata-gex-GRCm39-2024-A` (genome 2.79 Gb, 61
  contigs, 33,696 genes). STAR + rustar build their indexes from the refdata
  `fasta/genome.fa` + `genes/genes.gtf` (`--sjdbOverhang 89`); CellRanger uses
  the refdata directly.
- **Data**: 5k Mouse PBMCs, **5′ GEM-X** (SC5P-R2-v3); first **10,000,000 read
  pairs** of the GEX library — identical reads for all three tools.
- **Solo params** (CellRanger-matching, 5′): `--soloType CB_UMI_Simple`,
  CB 16 / UMI 12, `--soloStrand Reverse`, whitelist `3M-5pgex-jan-2023`,
  `--soloFeatures Gene`, `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`,
  `--soloUMIfiltering MultiGeneUMI_CR`, `--soloUMIdedup 1MM_CR`.
- **Environment**: Docker (colima) on Apple-Silicon macOS, **everything x86_64
  via Rosetta** (CellRanger is x86_64-only), 14 cores / 40 GB. All absolute
  times are inflated ~2–3× by emulation; the *relative* picture holds.
- **Tooling**: CellRanger 10.0.0, STAR 2.7.10b, rustar-aligner (this branch).
  Driver: [`test/solo_bench.py`](../test/solo_bench.py) (each step under
  `/usr/bin/time -v`), image [`test/Dockerfile.bench`](../test/Dockerfile.bench).

## Results

| Tool | Index build | Count (align+quant) | Peak RSS | Raw barcodes | Genes | Total UMIs |
|------|------------:|--------------------:|---------:|-------------:|------:|-----------:|
| **CellRanger 10.0.0** | (prebuilt) | 356 s | 12.5 GB | 161,465 | 17,258 | 4,843,682 |
| **STARsolo 2.7.10b** | 3,626 s | 152 s | 30 GB | 143,490 | 15,675 | 4,067,946 |
| **rustar-aligner** | 2,801 s | **670 s** | 37 GB | 156,258 | 16,278 | 4,219,582 |

CellRanger reported: 3,858 cells, 599 median genes/cell, 88.5 % valid barcodes,
58.5 % reads mapped to transcriptome.

### Correctness

On identical reads, rustar's raw matrix is in line with the references:
**4,219,582 UMIs** (exonic `Gene`), ~4 % above STARsolo's 4,067,946 (also exonic
`Gene`). CellRanger's 4,843,682 is higher because it counts **intronic** reads by
default (`include-introns`), whereas `--soloFeatures Gene` is exonic-only.
rustar's read-stage barcode match rate was **86 % exact** on this real data.

### The buffered-I/O fix

The first rustar count run took 1,774 s. A breakdown showed the raw-matrix write
dominated:

```
                  before        after
matrix write:    1,306 s   →       3 s     (~435×; byte-identical output)
align (10M):       402 s   →     627 s     (unchanged logic; emulation variance)
count total:     1,774 s   →     670 s
```

Cause: `write_barcodes` / `write_matrix_mtx` wrote to a raw `std::fs::File`
(unbuffered) — one `write(2)` syscall per line, so `barcodes.tsv` (the full
3,686,400-barcode whitelist) cost ~3.7M syscalls, amplified by Rosetta+virtiofs.
Fix: wrap the writers in `BufWriter` + a no-alloc barcode unpack
(`unpack_barcode_into`). The write dropped to ~3 s.

## Notes & limitations

- **Index build**: rustar (2,801 s) was *faster* than STARsolo (3,626 s) under
  emulation; CellRanger ships a prebuilt index (its 356 s "count" includes the
  internal STAR alignment + cell calling + full metrics).
- **Memory**: rustar's 37 GB peak is dominated by the **loaded index (~27 GB:
  5.4 B-entry SA for the 2.79 Gb genome)** plus the alignment working set — *not*
  the matrix build (Step 1 per-cell `build_matrix` already bounds that). Reducing
  the peak further is about the SA representation and alignment buffers, not the
  matrix.
- **Read count**: 10M (of ~200M total) keeps the run tractable and memory under
  the 40 GB cap. Stats scale with depth (CellRanger called 3,858 cells at this
  subsample vs the dataset's ~4,725).

## Reproduce

```bash
brew install colima docker && colima start --cpu 14 --memory 40 --vm-type vz --vz-rosetta
# build the amd64 image (colima can't build amd64 directly; run+commit a base):
docker run --platform linux/amd64 --name b rust:1-bookworm \
  bash -c "apt-get update -qq && apt-get install -y -qq rna-star python3 procps time"
docker commit b rustar-bench-amd64 && docker rm -f b
# then run test/solo_bench.py inside it with the ref/whitelist/fastqs mounted
# (see test/solo_bench.py header for the full argument list).
```
