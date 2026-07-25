//! Deterministic synthetic genome + read generator for PGO training data.
//!
//! Produces a small, fully reproducible (LCG-seeded, no OS randomness) FASTA +
//! FASTQ pair that exercises the aligner's hot code paths — exact/mismatched
//! seed extension, GT-AG splice detection across multiple introns, and
//! multi-mapping window/cluster handling via a duplicated repeat — so a
//! profile-guided-optimization training run collects representative branch/
//! block frequencies without needing a real reference genome or FASTQ dataset
//! checked into the repo or fetched over the network.
//!
//! Usage:
//!   gen_pgo_fixture --out-dir <dir> [--n-reads 50000] [--read-len 100]
//!
//! Writes `<out-dir>/genome.fa` and `<out-dir>/reads.fastq`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// LCG pseudo-random sequence generator (same recurrence as
/// `tests/alignment_features.rs`, so genome layout stays comprehensible).
fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let bases: [u8; 4] = *b"ACGT";
    let mut state = seed;
    let mut seq = Vec::with_capacity(length);
    for _ in 0..length {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(bases[((state >> 16) & 3) as usize]);
    }
    seq
}

/// A small xorshift-based PRNG (deterministic, no OS entropy) for read
/// sampling decisions — separate from `lcg_seq` so genome layout and read
/// sampling don't interact.
struct Xorshift(u64);
impl Xorshift {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn next_range(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// chr1: background with 3 planted GT-AG introns (exercises splice detection
/// + stitch_recurse across multiple exon boundaries).
fn build_chr1() -> Vec<u8> {
    let mut g = lcg_seq(88_888, 20_000);
    // Three introns at well-separated positions, each with a distinct body
    // seed so no two introns/exons are identical sequences.
    for (i, &start) in [4_000usize, 9_000, 14_000].iter().enumerate() {
        let exon_seed = 11_111 + i as u32 * 1000;
        let intron_seed = 22_222 + i as u32 * 1000;
        let exon = lcg_seq(exon_seed, 50);
        let intron_body = lcg_seq(intron_seed, 196);
        g[start..start + 50].copy_from_slice(&exon);
        g[start + 50] = b'G';
        g[start + 51] = b'T';
        g[start + 52..start + 248].copy_from_slice(&intron_body);
        g[start + 248] = b'A';
        g[start + 249] = b'G';
    }
    g
}

/// chr2: background with a 1 kb repeat block duplicated at two positions
/// (exercises multi-mapping / cluster_seeds window + overlap-dedup logic).
fn build_chr2() -> Vec<u8> {
    let mut g = lcg_seq(55_555, 20_000);
    let repeat = lcg_seq(66_666, 1_000);
    g[2_000..3_000].copy_from_slice(&repeat);
    g[15_000..16_000].copy_from_slice(&repeat);
    g
}

fn write_fasta(path: &PathBuf, chrs: &[(&str, &[u8])]) -> std::io::Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    for (name, seq) in chrs {
        writeln!(f, ">{name}")?;
        // 70-column wrapping, matching typical FASTA convention.
        for chunk in seq.chunks(70) {
            f.write_all(chunk)?;
            f.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b,
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn write_reads(
    path: &PathBuf,
    chr1: &[u8],
    chr2: &[u8],
    n_reads: usize,
    read_len: usize,
) -> std::io::Result<()> {
    let mut f = BufWriter::new(File::create(path)?);
    let mut rng = Xorshift(0x9E3779B97F4A7C15);
    let intron_starts = [4_000usize, 9_000, 14_000];
    let intron_len = 200; // 2 (GT) + 196 (body) + 2 (AG)

    for i in 0..n_reads {
        let roll = rng.next_f64();
        let seq: Vec<u8> = if roll < 0.55 {
            // Exact/near-exact match sampled from chr1 background (0-2 random
            // mismatches to exercise extend_alignment's mismatch/early-exit path).
            let max_start = chr1.len().saturating_sub(read_len);
            let start = rng.next_range(max_start.max(1));
            let mut s = chr1[start..start + read_len].to_vec();
            let n_mm = rng.next_range(3); // 0, 1, or 2 mismatches
            for _ in 0..n_mm {
                let pos = rng.next_range(read_len);
                let bases = *b"ACGT";
                s[pos] = bases[rng.next_range(4)];
            }
            s
        } else if roll < 0.80 {
            // Spliced read crossing one of chr1's planted introns: half the
            // read from the exon before, half from the exon after — exercises
            // GT-AG detection + multi-exon stitching.
            let intron_idx = rng.next_range(intron_starts.len());
            let donor = intron_starts[intron_idx] + 50; // intron start (GT)
            let acceptor = donor + intron_len; // first base after AG
            let left_len = read_len / 2;
            let right_len = read_len - left_len;
            let mut s = Vec::with_capacity(read_len);
            s.extend_from_slice(&chr1[donor - left_len..donor]);
            s.extend_from_slice(&chr1[acceptor..acceptor + right_len]);
            s
        } else if roll < 0.95 {
            // Read from the chr2 duplicated repeat block — multi-maps to
            // both copies, exercising cluster_seeds window/overlap handling.
            let repeat_starts = [2_000usize, 15_000];
            let which = rng.next_range(2);
            let start = repeat_starts[which] + rng.next_range(1_000 - read_len);
            chr2[start..start + read_len].to_vec()
        } else {
            // Pure random noise — exercises the "no seed found" / unmapped path.
            let bases = *b"ACGT";
            (0..read_len).map(|_| bases[rng.next_range(4)]).collect()
        };

        // Reverse-complement roughly half the reads so both search directions
        // get exercised.
        let final_seq = if rng.next_f64() < 0.5 { rc(&seq) } else { seq };
        let qual = "I".repeat(read_len);

        writeln!(f, "@pgo_train.{i}")?;
        f.write_all(&final_seq)?;
        writeln!(f)?;
        writeln!(f, "+")?;
        writeln!(f, "{qual}")?;
    }
    Ok(())
}

fn main() -> std::io::Result<()> {
    let mut out_dir: Option<PathBuf> = None;
    let mut n_reads: usize = 50_000;
    let mut read_len: usize = 100;

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--out-dir" => out_dir = args.next().map(PathBuf::from),
            "--n-reads" => {
                n_reads = args.next().and_then(|s| s.parse().ok()).unwrap_or(n_reads);
            }
            "--read-len" => {
                read_len = args.next().and_then(|s| s.parse().ok()).unwrap_or(read_len);
            }
            other => {
                eprintln!("gen_pgo_fixture: unknown argument '{other}'");
                std::process::exit(2);
            }
        }
    }

    let out_dir = out_dir.unwrap_or_else(|| {
        eprintln!("usage: gen_pgo_fixture --out-dir <dir> [--n-reads 50000] [--read-len 100]");
        std::process::exit(2);
    });
    std::fs::create_dir_all(&out_dir)?;

    let chr1 = build_chr1();
    let chr2 = build_chr2();

    write_fasta(
        &out_dir.join("genome.fa"),
        &[("chr1", chr1.as_slice()), ("chr2", chr2.as_slice())],
    )?;
    write_reads(
        &out_dir.join("reads.fastq"),
        &chr1,
        &chr2,
        n_reads,
        read_len,
    )?;

    eprintln!(
        "gen_pgo_fixture: wrote {} reads (len={}) over a {}bp genome (chr1+chr2) to {}",
        n_reads,
        read_len,
        chr1.len() + chr2.len(),
        out_dir.display()
    );
    Ok(())
}
