//! In-tree deterministic RNG (splitmix64), replacing the external `rand` crate.
//!
//! rustar-aligner's only randomness is tie-breaking / Monte-Carlo where the exact
//! stream is not STAR-faithful anyway (STAR uses mt19937; we already diverge and
//! report a tie-adjusted metric). What matters is **determinism** — reproducible
//! regardless of thread scheduling — which a fixed-seed splitmix64 gives, without
//! pulling in `rand` + its `getrandom`/`zerocopy`/`ppv-lite86` dependency chain.
//!
//! splitmix64 is Steele, Lea & Flood's generator (the same one `SplitMix64`/the
//! seeding of xoshiro use); constants are the reference values.

/// One splitmix64 step over `state`, returning the mixed output.
#[inline]
pub fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A minimal seedable deterministic RNG.
#[derive(Debug, Clone)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Seed the generator. Any seed (including 0) is valid.
    #[inline]
    pub fn seed(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next raw 64-bit value.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        splitmix64(&mut self.state)
    }

    /// Uniform `f64` in `[0, 1)` using the top 53 bits (full mantissa precision).
    #[inline]
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform integer in `[0, n)`. `n` must be > 0. Uses the float draw (`uniform`
    /// is strictly < 1, so the result is always in range); the tiny modulo bias is
    /// irrelevant for tie-breaking / Monte-Carlo use.
    #[inline]
    pub fn below(&mut self, n: usize) -> usize {
        debug_assert!(n > 0);
        ((self.uniform() * n as f64) as usize).min(n - 1)
    }
}

/// Deterministic in-place Fisher–Yates shuffle of `items`, seeded by `seed`.
pub fn shuffle_deterministic<T>(items: &mut [T], seed: u64) {
    let mut rng = SplitMix64::seed(seed);
    for i in (1..items.len()).rev() {
        let j = rng.below(i + 1);
        items.swap(i, j);
    }
}

/// Sample an index from ascending cumulative weights (prefix sums; `cumulative.last()`
/// is the total) via a uniform draw in `[0, total)`. Distribution-equivalent to
/// `rand::distr::weighted::WeightedIndex`. Returns an index in `[0, cumulative.len())`.
#[inline]
pub fn sample_cumulative(cumulative: &[f64], rng: &mut SplitMix64) -> usize {
    let total = *cumulative.last().unwrap_or(&0.0);
    let u = rng.uniform() * total;
    cumulative
        .partition_point(|&c| c <= u)
        .min(cumulative.len().saturating_sub(1))
}

/// Build ascending cumulative weights (prefix sums) from `weights`, for
/// [`sample_cumulative`].
pub fn cumulative_weights(weights: &[f64]) -> Vec<f64> {
    let mut c = Vec::with_capacity(weights.len());
    let mut acc = 0.0f64;
    for &w in weights {
        acc += w;
        c.push(acc);
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_across_calls() {
        let mut a = SplitMix64::seed(12345);
        let mut b = SplitMix64::seed(12345);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn uniform_in_range() {
        let mut r = SplitMix64::seed(7);
        for _ in 0..10_000 {
            let u = r.uniform();
            assert!((0.0..1.0).contains(&u));
        }
    }

    #[test]
    fn below_in_range() {
        let mut r = SplitMix64::seed(99);
        for n in 1..50usize {
            for _ in 0..200 {
                assert!(r.below(n) < n);
            }
        }
    }

    #[test]
    fn shuffle_is_a_permutation_and_deterministic() {
        let mut a: Vec<u32> = (0..50).collect();
        let mut b = a.clone();
        shuffle_deterministic(&mut a, 42);
        shuffle_deterministic(&mut b, 42);
        assert_eq!(a, b, "same seed -> same shuffle");
        let mut sorted = a.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..50).collect::<Vec<_>>(), "still a permutation");
    }

    #[test]
    fn sample_cumulative_respects_weights() {
        // Weight index 2 heavily; it should dominate the samples.
        let cum = cumulative_weights(&[1.0, 1.0, 100.0, 1.0]);
        let mut r = SplitMix64::seed(1);
        let mut counts = [0usize; 4];
        for _ in 0..10_000 {
            counts[sample_cumulative(&cum, &mut r)] += 1;
        }
        assert!(counts[2] > counts[0] + counts[1] + counts[3]);
    }
}
