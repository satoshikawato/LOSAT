//! EXPERIMENT (LOSAT_X_QDUSTPAR, shadow LOSAT_X_QDUSTPARSHADOW): the DUST masks of a query batch
//! made on the search pool, one query per job, collected in query order.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/dust_filter.cpp:92-102,171-180
//! ```c++
//! void s_CombineDustMasksWithUserProvidedMasks(...)
//! {
//!     CSymDustMasker duster(level, window, linker);
//!     CRef<CPacked_seqint> masked_locations =
//!         duster.GetMaskedInts(*query_id, data);
//! ...
//!     for (size_t i = 0; i < queries.Size(); i++) {
//!         ...
//!         s_CombineDustMasksWithUserProvidedMasks(data, queries.GetQuerySeqLoc(i),
//!                                                 queries.GetScope(i), query_id,
//!                                                 masks, level, window, linker);
//! ```
//! NCBI makes a fresh `CSymDustMasker` for each query of the batch and masks the queries one
//! after the other. `apply_dust_masking` does the same, serially, in query order. For a batch of
//! thousands of queries that is most of the time before the lookup table is built, and no other
//! thread works meanwhile.
//!
//! With the switch on, and a LOSAT search pool of two or more threads current, the per-query map
//! `DustMasker::new(level, window, linker).mask_sequence(seq)` is an ordered parallel map
//! (`par_iter().map().collect()`, which returns the results in the order of the input) run with
//! `SearchPool::install`, so on the pool of the search that calls it: the whole search pool, or
//! the private pool of a LOSAT_X_QBATCHPAR slot with two or more threads. Rayon's global pool is
//! never used. A search without a pool (one thread, a slot of one thread, a build without the
//! `parallel` feature, WebAssembly without `wasm-threads`), and a batch of fewer than two
//! queries, run the serial map.
//!
//! Why the masks are the same: the mask of a query is a function of that query's letters and of
//! the three DUST parameters only. `DustMasker::new` builds `level`, `window`, `linker`, `low_k`
//! and a threshold table from the three parameters and keeps no other state; every call of
//! `mask_subsequence` makes its own `NcbiLfgRandom` (the random numbers that resolve `N`),
//! windows and perfect-interval lists as locals, which is the reset of NCBI's fresh
//! `CSymDustMasker` for each query. The statics of `utils/dust.rs` are the switches
//! `LOSAT_X_DUSTFAST`, `LOSAT_X_DUSTRING` and `LOSAT_X_DUSTSHADOW`, which are read once and never
//! written again. So no query's mask depends on another query or on the order in which they are
//! masked, and the results are collected in query order. The verbose lines, the totals and the
//! empty case stay in `apply_dust_masking`, in the same order. No floating-point operation
//! changes order.
//!
//! Programs that reach it: only BLASTN (megablast, dc-megablast, blastn, blastn-short with
//! `-dust yes`), through `coordination::query_masks` (every query batch of `search_query_batch`)
//! and `chunk_query_masks` (the parts of a split query). The `-db` search of BLASTN uses the same
//! functions, so it is reached too. BLASTP, BLASTX, TBLASTN and TBLASTX do not call
//! `apply_dust_masking`; they mask with SEG.
//!
//! The shadow also makes the masks serially and compares them element for element. On a
//! difference it prints `[X_SHADOW] QDUSTPAR ... DIFFERENT` and panics; otherwise `main` prints
//! `[X_SHADOW] QDUSTPAR batches= queries= intervals= identical` at exit. With `LOSAT_TIMING`, `main`
//! prints `[TIMING] query_dust: <seconds>s`, and with `LOSAT_X_STATS` `[X_STATS] QDUSTPAR
//! dust_s=<seconds>`: the wall time of the DUST step of all batches, with or without the switch.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::utils::dust::{DustMasker, MaskedInterval};
use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

// No NCBI counterpart: reads the switches once (any value turns them on; the shadow implies the
// switch; neither is in the strict set). They choose between two LOSAT code paths with the same
// result and change no value NCBI computes.
/// 0 = off, 1 = LOSAT_X_QDUSTPAR, 2 = LOSAT_X_QDUSTPARSHADOW.
fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_QDUSTPARSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_QDUSTPAR").is_some() {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: counters of the shadow's summary line.
static BATCHES: AtomicU64 = AtomicU64::new(0);
static QUERIES: AtomicU64 = AtomicU64::new(0);
static INTERVALS: AtomicU64 = AtomicU64::new(0);

/// The DUST masks of `seqs`, one list per sequence, in order: the serial map of NCBI's loop.
fn serial(seqs: &[&[u8]], level: u32, window: usize, linker: usize) -> Vec<Vec<MaskedInterval>> {
    seqs.iter()
        .map(|seq| DustMasker::new(level, window, linker).mask_sequence(seq))
        .collect()
}

/// Whether the masks of `seq_count` queries are made on `pool`: two or more queries and a LOSAT
/// search pool of two or more threads.
fn on_pool(pool: &SearchPool<'_>, seq_count: usize) -> bool {
    seq_count > 1 && pool.enabled() && pool.threads() > 1
}

/// The parallel map on `pool`; the same lists as `serial`. Only called when `on_pool`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn parallel(
    seqs: &[&[u8]],
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
) -> Vec<Vec<MaskedInterval>> {
    pool.install(|| {
        seqs.par_iter()
            .map(|seq| DustMasker::new(level, window, linker).mask_sequence(seq))
            .collect()
    })
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn parallel(
    seqs: &[&[u8]],
    level: u32,
    window: usize,
    linker: usize,
    _pool: &SearchPool<'_>,
) -> Vec<Vec<MaskedInterval>> {
    serial(seqs, level, window, linker)
}

/// The masks of `seqs` as the switch `mode` makes them (0 = off, 1 = on, 2 = on and compared with
/// the serial map). Separate from the environment so that tests choose the mode.
pub(crate) fn masks_with_mode(
    seqs: &[&[u8]],
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
    mode: u8,
) -> Vec<Vec<MaskedInterval>> {
    if mode == 0 || !on_pool(pool, seqs.len()) {
        return serial(seqs, level, window, linker);
    }
    let masks = parallel(seqs, level, window, linker, pool);
    if mode == 2 {
        let reference = serial(seqs, level, window, linker);
        let intervals: usize = masks.iter().map(Vec::len).sum();
        BATCHES.fetch_add(1, Ordering::Relaxed);
        QUERIES.fetch_add(seqs.len() as u64, Ordering::Relaxed);
        INTERVALS.fetch_add(intervals as u64, Ordering::Relaxed);
        if masks != reference {
            let first = masks
                .iter()
                .zip(&reference)
                .position(|(fast, slow)| fast != slow);
            eprintln!(
                "[X_SHADOW] QDUSTPAR batches={} queries={} intervals={} DIFFERENT (first query {:?}: {:?} / {:?})",
                BATCHES.load(Ordering::Relaxed),
                QUERIES.load(Ordering::Relaxed),
                INTERVALS.load(Ordering::Relaxed),
                first,
                first.map(|i| &masks[i]),
                first.map(|i| &reference[i])
            );
            panic!("LOSAT_X_QDUSTPARSHADOW: the DUST masks made on the pool differ from the serial ones");
        }
    }
    masks
}

// No NCBI counterpart: the wall time of the batches' DUST step for LOSAT_TIMING (any switch
// state), so that a run with the switch and one without can be compared; it changes no value NCBI
// computes.
static DUST_NANOS: AtomicU64 = AtomicU64::new(0);

/// Whether the DUST step is timed: with LOSAT_TIMING (the `[TIMING]` line) or LOSAT_X_STATS (the
/// `[X_STATS]` line, which costs no more than two clock reads per batch; LOSAT_TIMING also times
/// every ungapped extension on shared counters).
fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        std::env::var_os("LOSAT_TIMING").is_some() || std::env::var_os("LOSAT_X_STATS").is_some()
    })
}

/// The masks of `seqs` as `apply_dust_masking` makes them: the environment's switch decides.
pub(crate) fn masks(
    seqs: &[&[u8]],
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
) -> Vec<Vec<MaskedInterval>> {
    let start = timing_on().then(std::time::Instant::now);
    let masks = masks_with_mode(seqs, level, window, linker, pool, mode());
    if let Some(start) = start {
        DUST_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
    masks
}

/// The shadow's summary line (printed by `main` at exit), only with LOSAT_X_QDUSTPARSHADOW, and
/// with LOSAT_TIMING the wall time of the DUST step of all batches.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    if timing_on() {
        let seconds = DUST_NANOS.load(Ordering::Relaxed) as f64 / 1e9;
        if std::env::var_os("LOSAT_TIMING").is_some() {
            eprintln!("[TIMING] query_dust: {seconds:.3}s");
        }
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            eprintln!("[X_STATS] QDUSTPAR dust_s={seconds:.3}");
        }
    }
    if mode() == 2 {
        eprintln!(
            "[X_SHADOW] QDUSTPAR batches={} queries={} intervals={} identical",
            BATCHES.load(Ordering::Relaxed),
            QUERIES.load(Ordering::Relaxed),
            INTERVALS.load(Ordering::Relaxed)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::threading::with_search_pool;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// A random sequence with low-complexity runs, `N`s, lower case and other IUPAC letters.
    fn sequence(rng: &mut Rng) -> Vec<u8> {
        let len = match rng.below(6) {
            0 => rng.below(4) as usize,
            1 => 3 + rng.below(40) as usize,
            _ => 100 + rng.below(3000) as usize,
        };
        let mut seq = Vec::with_capacity(len);
        while seq.len() < len {
            match rng.below(8) {
                0 => {
                    let unit =
                        [b"A".as_slice(), b"AT", b"CAG", b"GGGA", b"TTTTC"][rng.below(5) as usize];
                    for _ in 0..(5 + rng.below(60)) {
                        seq.extend_from_slice(unit);
                    }
                }
                1 => seq.extend(std::iter::repeat_n(b'N', 1 + rng.below(12) as usize)),
                2 => seq.push(b"RYKMSWnacgt"[rng.below(11) as usize]),
                _ => seq.push(b"ACGT"[rng.below(4) as usize]),
            }
        }
        seq.truncate(len);
        seq
    }

    /// Parallel == serial on random sequences with low-complexity runs at 1, 2, 3 and 7 threads,
    /// with empty sequences and batches of one and of two.
    #[test]
    fn pool_masks_equal_serial_masks() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let mut seqs: Vec<Vec<u8>> = (0..90).map(|_| sequence(&mut rng)).collect();
        seqs.push(Vec::new());
        seqs.push(b"AC".to_vec());
        seqs.push(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec());
        let masked: usize = seqs
            .iter()
            .map(|s| DustMasker::new(20, 64, 1).mask_sequence(s).len())
            .sum();
        assert!(
            masked > 20,
            "the sequences must have masked regions ({masked})"
        );
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        for (level, window, linker) in [(20, 64, 1), (10, 32, 2), (30, 16, 1)] {
            let expect = serial(&refs, level, window, linker);
            for threads in [1usize, 2, 3, 7] {
                with_search_pool(threads, "blastn", |pool| {
                    for take in [refs.len(), 2, 1, 0] {
                        let part = &refs[..take];
                        let want = serial(part, level, window, linker);
                        for mode in [1u8, 2] {
                            let got = masks_with_mode(part, level, window, linker, pool, mode);
                            assert_eq!(got, want, "threads {threads} take {take} mode {mode}");
                        }
                    }
                    let got = masks_with_mode(&refs, level, window, linker, pool, 1);
                    assert_eq!(got, expect, "threads {threads}");
                    Ok(())
                })
                .unwrap();
            }
        }
    }

    /// The parallel path acts only where a search pool of two or more threads is current, and
    /// runs on that pool's threads.
    #[cfg(feature = "parallel")]
    #[test]
    fn only_a_search_pool_of_two_threads_acts() {
        with_search_pool(1, "blastn", |pool| {
            assert!(!on_pool(pool, 10));
            Ok(())
        })
        .unwrap();
        with_search_pool(3, "blastn", |pool| {
            assert!(on_pool(pool, 10));
            assert!(!on_pool(pool, 1));
            assert!(!on_pool(pool, 0));
            let seqs: Vec<Vec<u8>> = (0..30).map(|i| vec![b"ACGT"[i % 4]; 200]).collect();
            let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
            // `install` runs on the search pool (3 threads), not on rayon's global pool.
            assert_eq!(pool.install(rayon::current_num_threads), 3);
            let got = masks_with_mode(&refs, 20, 64, 1, pool, 1);
            assert_eq!(got, serial(&refs, 20, 64, 1));
            Ok(())
        })
        .unwrap();
    }
}
