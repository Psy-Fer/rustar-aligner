//! End-to-end checks for the STAR parameters closed in this change, so that
//! "accepted by the CLI" and "actually changes the output" stay different
//! claims.
//!
//! Genome: 20 kb of LCG(88888) background, the same generator the other
//! integration tests use, with one 100 bp segment planted at four extra
//! positions so a read from it multimaps.

use assert_cmd::cargo::cargo_bin_cmd;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const READ_LEN: usize = 100;
/// Where the repeated 100 bp block is copied to. Four copies plus the
/// original gives a read from it five alignments.
const REPEAT_SRC: usize = 1_000;
const REPEAT_COPIES: [usize; 4] = [4_000, 8_000, 12_000, 16_000];
/// A unique region used for the mismatch tests.
const UNIQUE_START: usize = 6_000;

fn lcg_seq(seed: u32, length: usize) -> Vec<u8> {
    let bases: [u8; 4] = *b"ACGT";
    let mut state = seed;
    let mut seq = Vec::with_capacity(length);
    for _ in 0..length {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        seq.push(bases[((state >> 16) & 3) as usize]);
    }
    seq
}

fn build_genome() -> Vec<u8> {
    let mut genome = lcg_seq(88888, 20_000);
    let block: Vec<u8> = genome[REPEAT_SRC..REPEAT_SRC + READ_LEN].to_vec();
    for &dst in &REPEAT_COPIES {
        genome[dst..dst + READ_LEN].copy_from_slice(&block);
    }
    genome
}

fn write_fasta(dir: &Path, genome: &[u8]) -> PathBuf {
    let path = dir.join("genome.fa");
    let mut f = fs::File::create(&path).unwrap();
    writeln!(f, ">chr1").unwrap();
    f.write_all(genome).unwrap();
    writeln!(f).unwrap();
    path
}

fn build_index(fasta: &Path, genome_dir: &Path) {
    fs::create_dir_all(genome_dir).unwrap();
    cargo_bin_cmd!("rustar-aligner")
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            genome_dir.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta.to_str().unwrap(),
            "--genomeSAindexNbases",
            "7",
        ])
        .assert()
        .success();
}

fn write_fastq(path: &Path, reads: &[(String, Vec<u8>)]) {
    let mut f = fs::File::create(path).unwrap();
    for (name, seq) in reads {
        writeln!(f, "@{name}").unwrap();
        f.write_all(seq).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(seq.len())).unwrap();
    }
}

/// Run the aligner and return the SAM records (header lines dropped).
fn align(genome_dir: &Path, fastq: &Path, prefix: &str, extra: &[&str]) -> Vec<String> {
    let mut cmd = cargo_bin_cmd!("rustar-aligner");
    cmd.args([
        "--runMode",
        "alignReads",
        "--genomeDir",
        genome_dir.to_str().unwrap(),
        "--readFilesIn",
        fastq.to_str().unwrap(),
        "--outFileNamePrefix",
        prefix,
    ]);
    cmd.args(extra);
    cmd.assert().success();

    fs::read_to_string(format!("{prefix}Aligned.out.sam"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(str::to_string)
        .collect()
}

fn mismatched_read(genome: &[u8], n_mismatches: usize) -> Vec<u8> {
    let mut seq = genome[UNIQUE_START..UNIQUE_START + READ_LEN].to_vec();
    // Spread the substitutions out so they cannot all be soft-clipped off one
    // end instead of being counted as mismatches.
    for i in 0..n_mismatches {
        let pos = 15 + i * 20;
        seq[pos] = match seq[pos] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    seq
}

// ── outFilterMismatchNoverReadLmax ──────────────────────────────────────────

#[test]
fn mismatch_nover_read_lmax_filters_by_read_length_ratio() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = build_genome();
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);

    let fq = root.join("mm.fq");
    write_fastq(&fq, &[("mm3".to_string(), mismatched_read(&genome, 3))]);

    // Control: the read maps with the default ratio of 1.0.
    let permissive = align(
        &genome_dir,
        &fq,
        &format!("{}/permissive_", root.display()),
        &["--outFilterMismatchNoverReadLmax", "1.0"],
    );
    assert_eq!(permissive.len(), 1, "read should map by default");
    let control_flag: u32 = permissive[0].split('\t').nth(1).unwrap().parse().unwrap();
    assert_eq!(
        control_flag & 0x4,
        0,
        "control read must be mapped: {}",
        permissive[0]
    );

    // 0.01 of 100 bases allows one mismatch; the read carries three.
    let strict = align(
        &genome_dir,
        &fq,
        &format!("{}/strict_", root.display()),
        &["--outFilterMismatchNoverReadLmax", "0.01"],
    );
    let mapped = strict
        .iter()
        .filter(|r| {
            let flag: u32 = r.split('\t').nth(1).unwrap().parse().unwrap();
            flag & 0x4 == 0
        })
        .count();
    assert_eq!(
        mapped, 0,
        "a 3-mismatch read must fail a 0.01 read-length ratio: {strict:?}"
    );
}

// ── alignTranscriptsPerReadNmax ─────────────────────────────────────────────

#[test]
fn transcripts_per_read_nmax_caps_the_alignments_kept() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = build_genome();
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);

    let fq = root.join("multi.fq");
    let repeat = genome[REPEAT_SRC..REPEAT_SRC + READ_LEN].to_vec();
    write_fastq(&fq, &[("multi".to_string(), repeat)]);

    // Control: all five copies are reported when nothing caps them.
    let uncapped = align(
        &genome_dir,
        &fq,
        &format!("{}/uncapped_", root.display()),
        &["--outFilterMultimapNmax", "20"],
    );
    assert!(
        uncapped.len() >= 3,
        "the planted repeat should multimap, got {} record(s)",
        uncapped.len()
    );

    let capped = align(
        &genome_dir,
        &fq,
        &format!("{}/capped_", root.display()),
        &[
            "--outFilterMultimapNmax",
            "20",
            "--alignTranscriptsPerReadNmax",
            "2",
        ],
    );
    assert!(
        capped.len() < uncapped.len(),
        "the cap must reduce the alignments kept: {} capped vs {} uncapped",
        capped.len(),
        uncapped.len()
    );
    assert!(
        capped.len() <= 2,
        "at most 2 alignments may survive --alignTranscriptsPerReadNmax 2, got {}",
        capped.len()
    );
}

// ── alignSoftClipAtReferenceEnds ────────────────────────────────────────────

#[test]
fn soft_clip_at_reference_ends_no_prohibits_clipping_past_the_chromosome() {
    let dir = TempDir::new().unwrap();
    let root = dir.path();
    let genome = build_genome();
    let genome_dir = root.join("genome");
    build_index(&write_fasta(root, &genome), &genome_dir);

    // A read whose last 20 bases fall past the end of chr1: 80 bases of the
    // genome tail, then 20 bases that exist nowhere, so the aligner has to
    // soft-clip them past the reference end.
    let mut seq = genome[genome.len() - 80..].to_vec();
    seq.extend_from_slice(b"ACACACACACGTGTGTGTGT");
    let fq = root.join("edge.fq");
    write_fastq(&fq, &[("edge".to_string(), seq)]);

    let allowed = align(
        &genome_dir,
        &fq,
        &format!("{}/allowed_", root.display()),
        &["--alignSoftClipAtReferenceEnds", "Yes"],
    );
    let mapped_allowed = allowed
        .iter()
        .filter(|r| {
            let flag: u32 = r.split('\t').nth(1).unwrap().parse().unwrap();
            flag & 0x4 == 0
        })
        .count();
    assert_eq!(
        mapped_allowed, 1,
        "control: the read maps with a soft clip past the end: {allowed:?}"
    );
    assert!(
        allowed[0].split('\t').nth(5).unwrap().contains('S'),
        "control alignment should carry a soft clip: {}",
        allowed[0]
    );

    let prohibited = align(
        &genome_dir,
        &fq,
        &format!("{}/prohibited_", root.display()),
        &["--alignSoftClipAtReferenceEnds", "No"],
    );
    let mapped_prohibited = prohibited
        .iter()
        .filter(|r| {
            let flag: u32 = r.split('\t').nth(1).unwrap().parse().unwrap();
            flag & 0x4 == 0
        })
        .count();
    assert_eq!(
        mapped_prohibited, 0,
        "--alignSoftClipAtReferenceEnds No must reject it: {prohibited:?}"
    );
}
