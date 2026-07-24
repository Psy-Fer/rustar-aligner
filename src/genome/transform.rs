//! `--genomeTransformType Haploid` (STAR `Genome_transformGenome.cpp`): substitute VCF alleles into
//! the genome before the suffix array is built, so reads carrying the alternate allele align with
//! fewer mismatches.
//!
//! Ported from STAR-rs `crates/star-index/src/transform.rs`. Only **Haploid** (one allele per site,
//! including indels that shift coordinates and split the `transformGenomeBlocks.tsv` block map) is
//! implemented; **Diploid** (the genome duplicated into `_h1`/`_h2` haplotypes) is a follow-up. With
//! the only supported `--genomeTransformOutput` (`None`, the default) this is a pure `genomeGenerate`
//! transform: the aligner reports transformed-genome coordinates directly, so no align-time
//! back-transform is implemented either (that's a separate follow-up, gated in `Parameters::validate`).
//!
//! Unlike STAR-rs, which transforms an already-laid-out, bin-padded genome buffer, this operates on
//! rustar-aligner's `Vec<Chromosome>` (name + unpadded base-code sequence) directly, before
//! [`Genome::from_fasta`](super::Genome::from_fasta)'s own padding pass runs. Blocks are computed in
//! two stages: substitution happens in chromosome-local coordinates (no padding involved), then
//! globalized using [`compute_chr_starts`](super::compute_chr_starts) run once over the original
//! chromosome lengths and once over the transformed lengths — the same deterministic function
//! `Genome::from_fasta` itself uses, so the geometry is guaranteed to match the final on-disk index.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::io::fastq::encode_base;

use super::compute_chr_starts;
use super::fasta::Chromosome;

/// One VCF variant applied to a chromosome: 0-based chr-local `pos`, the reference length, and the
/// alternate allele bytes (raw ASCII letters; [`encode_base`] is applied at substitution time).
#[derive(Debug, Clone)]
pub struct Variant {
    pub pos: u64,
    pub ref_len: usize,
    pub alt: Vec<u8>,
}

impl Variant {
    /// STAR's `len = alt.size() - ref.size()` (0 for a SNV, `>0` insertion, `<0` deletion).
    fn len_delta(&self) -> i64 {
        self.alt.len() as i64 - self.ref_len as i64
    }
}

/// Parse a VCF for the Haploid transform (`Genome_transformGenome.cpp:41-107`): every record on a
/// known chromosome contributes its **first** alternate allele (no genotype filtering); `#` lines and
/// records on unknown chromosomes are skipped. Returns per-chromosome-index variants, each list sorted
/// by position with STAR's overlap filter applied (keep a variant only if it starts at or after the
/// end of the last kept one).
pub fn parse_vcf_haploid(text: &str, chr_name: &[String]) -> BTreeMap<usize, Vec<Variant>> {
    let index: BTreeMap<&str, usize> = chr_name
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    let mut per_chr: BTreeMap<usize, Vec<Variant>> = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 5 {
            continue;
        }
        let Some(&ci) = index.get(f[0]) else {
            continue;
        };
        let Ok(pos1): Result<u64, _> = f[1].parse() else {
            continue;
        };
        if pos1 == 0 {
            continue;
        }
        let ref_allele = f[3].as_bytes();
        // STAR takes the FIRST alternate allele for Haploid (`altV[0]`).
        let alt = f[4].split(',').next().unwrap_or(f[4]).as_bytes().to_vec();
        per_chr.entry(ci).or_default().push(Variant {
            pos: pos1 - 1,
            ref_len: ref_allele.len(),
            alt,
        });
    }

    for variants in per_chr.values_mut() {
        *variants = filter_sort_variants(std::mem::take(variants));
    }
    per_chr
}

/// Sort variants by position and apply STAR's overlap filter: keep a variant only if it starts at or
/// after the end of the previously kept variant's reference span.
fn filter_sort_variants(mut variants: Vec<Variant>) -> Vec<Variant> {
    variants.sort_by_key(|v| v.pos);
    let mut kept: Vec<Variant> = Vec::with_capacity(variants.len());
    let mut g0: u64 = 0;
    for v in variants {
        if v.pos >= g0 {
            g0 = v.pos + v.ref_len as u64;
            kept.push(v);
        }
    }
    kept
}

/// One chromosome's Haploid substitution, in chromosome-LOCAL coordinates (no padding): the
/// transformed sequence, and its blocks `[orig_local_start, len, new_local_start]` (globalized by the
/// caller). A chromosome with no variants yields one identity block spanning the whole sequence.
fn transform_one_chromosome(seq: &[u8], variants: &[Variant]) -> (Vec<u8>, Vec<[u64; 3]>) {
    if variants.is_empty() {
        return (seq.to_vec(), vec![[0, seq.len() as u64, 0]]);
    }

    let cl0 = seq.len() as u64;
    let mut gnew: Vec<u8> = Vec::with_capacity(seq.len());
    let mut blocks: Vec<[u64; 3]> = Vec::new();
    let mut iv = 0usize;
    let mut g0: u64 = 0;
    let mut g1: u64 = 0;
    blocks.push([g0, 0, g1]); // first block

    while g0 < cl0 {
        if g0 == variants[iv].pos {
            let v = &variants[iv];
            for &b in &v.alt {
                gnew.push(encode_base(b));
            }
            g0 += v.ref_len as u64;
            g1 += v.alt.len() as u64;
            if v.len_delta() != 0 {
                // Close the previous block; STAR's length formula, then open a new one.
                let last = blocks.last_mut().unwrap();
                last[1] = g0 - v.ref_len as u64 + v.ref_len.min(v.alt.len()) as u64 - last[0];
                blocks.push([g0, 0, g1]);
            }
            if iv < variants.len() - 1 {
                iv += 1;
            }
        } else {
            gnew.push(seq[g0 as usize]);
            g0 += 1;
            g1 += 1;
        }
    }
    if blocks.last().unwrap()[1] == 0 {
        let last = blocks.last_mut().unwrap();
        last[1] = g0 - last[0];
    }
    (gnew, blocks)
}

/// The result of transforming a genome's chromosome list: the substituted chromosomes, and the global
/// `[orig_start, length, new_start]` block map (STAR's `array<uint64,3>`), in chromosome order.
pub struct TransformedGenome {
    pub chromosomes: Vec<Chromosome>,
    pub blocks: Vec<[u64; 3]>,
}

/// Apply the Haploid VCF transform to a genome's chromosome list. `variants` maps a chromosome index
/// (matching `chromosomes`' order) to its filtered, sorted variants (from [`parse_vcf_haploid`]);
/// `chr_bin_nbits` is `--genomeChrBinNbits`, used only to globalize the block coordinates the same way
/// [`Genome::from_fasta`](super::Genome::from_fasta) will pad the transformed chromosomes.
pub fn transform_chromosomes(
    chromosomes: &[Chromosome],
    variants: &BTreeMap<usize, Vec<Variant>>,
    chr_bin_nbits: u32,
) -> TransformedGenome {
    let orig_lengths: Vec<u64> = chromosomes
        .iter()
        .map(|c| c.sequence.len() as u64)
        .collect();

    let empty: Vec<Variant> = Vec::new();
    let mut new_chromosomes = Vec::with_capacity(chromosomes.len());
    let mut local_blocks_per_chr: Vec<Vec<[u64; 3]>> = Vec::with_capacity(chromosomes.len());
    for (ci, chrom) in chromosomes.iter().enumerate() {
        let vs = variants.get(&ci).unwrap_or(&empty);
        let (seq, blocks) = transform_one_chromosome(&chrom.sequence, vs);
        new_chromosomes.push(Chromosome {
            name: chrom.name.clone(),
            sequence: seq,
        });
        local_blocks_per_chr.push(blocks);
    }

    let new_lengths: Vec<u64> = new_chromosomes
        .iter()
        .map(|c| c.sequence.len() as u64)
        .collect();
    let orig_chr_start = compute_chr_starts(&orig_lengths, chr_bin_nbits);
    let new_chr_start = compute_chr_starts(&new_lengths, chr_bin_nbits);

    let mut blocks = Vec::new();
    for (ci, local) in local_blocks_per_chr.into_iter().enumerate() {
        for [lo, len, ln] in local {
            blocks.push([orig_chr_start[ci] + lo, len, new_chr_start[ci] + ln]);
        }
    }

    TransformedGenome {
        chromosomes: new_chromosomes,
        blocks,
    }
}

/// Render `transformGenomeBlocks.tsv` (STAR's `transformBlocksWrite`): header `<nBlocks>\t-1`, then
/// one `new_start\tlength\torig_start` line per block (reverting the stored
/// `[orig_start, length, new_start]` order, for reverse conversion).
pub fn blocks_to_tsv(blocks: &[[u64; 3]]) -> String {
    let mut s = format!("{}\t-1\n", blocks.len());
    for b in blocks {
        let _ = writeln!(s, "{}\t{}\t{}", b[2], b[1], b[0]);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_first_alt_and_filters_overlaps() {
        let chr = vec!["chr1".to_string(), "chr2".to_string()];
        let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10\t.\tG\tA,T\t.\t.\t.
chrX\t5\t.\tG\tA\t.\t.\t.
chr2\t3\t.\tC\tT\t.\t.\t.
";
        let v = parse_vcf_haploid(vcf, &chr);
        assert_eq!(v[&0].len(), 1);
        assert_eq!(v[&0][0].pos, 9); // 0-based
        assert_eq!(v[&0][0].alt, b"A"); // first alt only
        assert_eq!(v[&1][0].pos, 2);
        assert_eq!(v[&1][0].alt, b"T");
        assert!(!v.contains_key(&2)); // chrX not in index; only 2 contigs
    }

    #[test]
    fn snv_substitution_and_identity_blocks() {
        let chromosomes = vec![
            Chromosome {
                name: "chr1".to_string(),
                sequence: vec![2, 2, 2, 2, 2], // GGGGG
            },
            Chromosome {
                name: "chr2".to_string(),
                sequence: vec![2, 2, 2, 2], // GGGG
            },
        ];
        let mut variants = BTreeMap::new();
        variants.insert(
            0usize,
            vec![Variant {
                pos: 2,
                ref_len: 1,
                alt: b"A".to_vec(),
            }],
        );
        let t = transform_chromosomes(&chromosomes, &variants, 18);
        assert_eq!(t.chromosomes[0].sequence[2], 0); // A code
        assert_eq!(t.chromosomes[0].sequence.len(), 5);
        assert_eq!(t.chromosomes[1].sequence.len(), 4);
        // One identity block per chromosome; lengths unchanged.
        let cs = compute_chr_starts(&[5, 4], 18);
        assert_eq!(t.blocks, vec![[cs[0], 5, cs[0]], [cs[1], 4, cs[1]]]);
    }

    #[test]
    fn indel_shifts_coordinates_and_splits_blocks() {
        // chr1 = 6 bases; a 1-base insertion at pos 2 (ref C -> alt CTT, +2) then a deletion at pos 4
        // (ref GA -> alt G, -1). Net length 6 + 2 - 1 = 7.
        let chromosomes = vec![Chromosome {
            name: "chr1".to_string(),
            sequence: vec![0, 0, 1, 2, 2, 0], // AACGGA
        }];
        let mut variants = BTreeMap::new();
        variants.insert(
            0usize,
            vec![
                Variant {
                    pos: 2,
                    ref_len: 1,
                    alt: b"CTT".to_vec(),
                },
                Variant {
                    pos: 4,
                    ref_len: 2,
                    alt: b"G".to_vec(),
                },
            ],
        );
        let t = transform_chromosomes(&chromosomes, &variants, 18);
        assert_eq!(t.chromosomes[0].sequence.len(), 7);
        let seq: Vec<u8> = t.chromosomes[0]
            .sequence
            .iter()
            .map(|&c| b"ACGT"[c as usize])
            .collect();
        assert_eq!(&seq, b"AACTTGG");
        // Blocks split at each indel: [orig_start, length, new_start] (global, chr1 starts at 0).
        assert_eq!(t.blocks[0], [0, 3, 0]);
    }

    #[test]
    fn blocks_to_tsv_reverts_column_order() {
        let tsv = blocks_to_tsv(&[[0, 5, 0], [5, 2, 7]]);
        assert_eq!(tsv, "2\t-1\n0\t5\t0\n7\t2\t5\n");
    }
}
