use rustar_aligner::cpu;
use rustar_aligner::params::Parameters;

/// Global allocator override — see the comment in `Cargo.toml`
/// next to the `mimalloc` dependency for the rationale. The
/// `#[global_allocator]` attribute applies to the whole binary; no
/// further work needed at allocation sites.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> anyhow::Result<()> {
    // mimalloc returns freed pages to the OS aggressively by default (~10 ms purge
    // delay). On this high-throughput aligner that means constant page re-faulting +
    // kernel zeroing — ~12% of alignment time (clear_page_erms) in profiling. We
    // default to never-purge (-1) for speed (~13% faster). The tradeoff: mimalloc then
    // retains a resident free-page reserve that inflates Max RSS (reusable scratch, not
    // data — the actual heap footprint is ~3 GB, vs a reclaimable ~25 GB mmap'd index).
    // Set RUSTAR_PURGE_DELAY_MS to a finite value (e.g. 1000) to cap that reserve at a
    // small wall-time cost. `mi_option_set` overrides the cached default, so it takes
    // effect even though mimalloc initialized before main(). `mi_option_purge_delay` is
    // option 15 in mimalloc v2.x (not a named const in libmimalloc-sys 0.1.49).
    // SAFETY: mi_option_t is a plain c_int; trivial FFI call.
    const MI_OPTION_PURGE_DELAY: std::os::raw::c_int = 15;
    let purge_ms: std::os::raw::c_long = std::env::var("RUSTAR_PURGE_DELAY_MS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(-1);
    unsafe { libmimalloc_sys::mi_option_set(MI_OPTION_PURGE_DELAY, purge_ms) };

    // Don't eagerly commit whole arenas. On Linux (an overcommit OS) mimalloc's
    // default (arena_eager_commit=2) commits a large arena pool up front — measured
    // as ~15 GB of committed-but-never-touched reserve, inflating RssAnon to ~17 GB
    // despite only ~2 GB of live data. Committing arenas on demand instead drops
    // RssAnon 17 GB -> 2 GB and Max RSS 41 -> 26 GB at ~3% wall cost (the on-demand
    // commit faults). See test/aws/mem_tune.sh. `mi_option_arena_eager_commit` is
    // option 4 in mimalloc v2.x. SAFETY: trivial c_int FFI call.
    const MI_OPTION_ARENA_EAGER_COMMIT: std::os::raw::c_int = 4;
    unsafe { libmimalloc_sys::mi_option_set(MI_OPTION_ARENA_EAGER_COMMIT, 0) };

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    cpu::check_cpu_compat()?;

    rustar_aligner::run(&Parameters::parse())
}
