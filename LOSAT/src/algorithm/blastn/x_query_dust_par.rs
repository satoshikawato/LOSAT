//! EXPERIMENT (LOSAT_X_QDUSTPAR): the DUST masks of a query batch
//! made on the search pool, the per-query map split over the pool's threads and collected in
//! query order.
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
//! (`par_iter().map().collect()`, which returns the results in the order of the input; rayon
//! splits the queries into ranges, so the masks are per query but the jobs are not) run with
//! `SearchPool::install`, so on the pool of the search that calls it. Rayon's global pool is never
//! used. A search without a pool (one thread, a build without the `parallel` feature, WebAssembly
//! without `wasm-threads`), and a batch of fewer than two queries, run the serial map. With `-dust no` there is no map at all (`apply_dust_masking`
//! makes no masks).
//!
//! Why the masks are the same: the mask of a query is a function of that query's letters and of
//! the three DUST parameters only. `DustMasker::new` builds `level`, `window`, `linker`, `low_k`
//! and a threshold table from the three parameters and keeps no other state; every call of
//! `mask_subsequence` makes its own `NcbiLfgRandom` (the random numbers that resolve `N`),
//! windows and perfect-interval lists as locals, which is the reset of NCBI's fresh
//! `CSymDustMasker` for each query. `utils/dust.rs` has no statics that a call writes. So no query's mask depends on another query or on the order in which they are
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
//! With `LOSAT_TIMING`, `main` prints `[TIMING] query_dust: <seconds>s`, and with `LOSAT_X_STATS`
//! `[X_STATS] QDUSTPAR dust_s=<seconds>`, only when at least one batch was masked (so not for
//! BLASTP, BLASTX, TBLASTN, TBLASTX or `-dust no`), with or without the switch. The seconds are
//! the sum of the wall times of every DUST call of the process (every batch's).

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use super::coordination::QueryLetters;
use crate::utils::dust::{DustMasker, MaskedInterval};
use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

// No NCBI counterpart: reads the switch once (LOSAT_X_QDUSTPAR is a strict-set switch since
// 2026-10-17: on unless it is `0`). It chooses between two LOSAT code paths with the same result
// and changes no value NCBI computes.
/// Whether LOSAT_X_QDUSTPAR is on.
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_QDUSTPAR"))
}

/// The DUST masks of `seqs`, one list per sequence, in order: the serial map of NCBI's loop.
fn serial<L: QueryLetters + ?Sized>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
) -> Vec<Vec<MaskedInterval>> {
    (0..seqs.count())
        .map(|index| {
            seqs.read(index, |seq| {
                DustMasker::new(level, window, linker).mask_sequence(seq)
            })
        })
        .collect()
}

/// Whether the masks of `seq_count` queries are made on `pool`: two or more queries and a LOSAT
/// search pool of two or more threads.
fn on_pool(pool: &SearchPool<'_>, seq_count: usize) -> bool {
    seq_count > 1 && pool.enabled() && pool.threads() > 1
}

/// The parallel map on `pool`; the same lists as `serial`. Only called when `on_pool`. `probe` is
/// called once per query, on the thread that masks it; production passes `|| {}`, which the
/// compiler removes (a unit test passes a closure that looks at the thread).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn parallel<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
    probe: &P,
) -> Vec<Vec<MaskedInterval>> {
    pool.install(|| {
        (0..seqs.count())
            .into_par_iter()
            .map(|index| {
                probe();
                seqs.read(index, |seq| {
                    DustMasker::new(level, window, linker).mask_sequence(seq)
                })
            })
            .collect()
    })
}

#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn parallel<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    _pool: &SearchPool<'_>,
    _probe: &P,
) -> Vec<Vec<MaskedInterval>> {
    serial(seqs, level, window, linker)
}

/// The masks of `seqs` as the switch `on` makes them. Separate from the environment so that tests
/// choose the switch.
pub(crate) fn masks_with_mode<L: QueryLetters + ?Sized>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
    on: bool,
) -> Vec<Vec<MaskedInterval>> {
    masks_with_probe(seqs, level, window, linker, pool, on, &|| {})
}

/// `masks_with_mode` with a `probe` that the parallel map calls once per query (see `parallel`).
fn masks_with_probe<L: QueryLetters + ?Sized, P: Fn() + Sync>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
    on: bool,
    probe: &P,
) -> Vec<Vec<MaskedInterval>> {
    if on {
        crate::utils::xreach::x_reach!("LOSAT_X_QDUSTPAR", on_pool(pool, seqs.count()));
    }
    if !on || !on_pool(pool, seqs.count()) {
        return serial(seqs, level, window, linker);
    }
    parallel(seqs, level, window, linker, pool, probe)
}

// No NCBI counterpart: the summed wall time of the DUST calls and their number, for LOSAT_TIMING
// (any switch state), so that a run with the switch and one without can be compared; they change
// no value NCBI computes.
static DUST_NANOS: AtomicU64 = AtomicU64::new(0);
static DUST_CALLS: AtomicU64 = AtomicU64::new(0);

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
pub(crate) fn masks<L: QueryLetters + ?Sized>(
    seqs: &L,
    level: u32,
    window: usize,
    linker: usize,
    pool: &SearchPool<'_>,
) -> Vec<Vec<MaskedInterval>> {
    let start = timing_on().then(std::time::Instant::now);
    let masks = masks_with_mode(seqs, level, window, linker, pool, on());
    if let Some(start) = start {
        DUST_NANOS.fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
        DUST_CALLS.fetch_add(1, Ordering::Relaxed);
    }
    masks
}

/// With LOSAT_TIMING / LOSAT_X_STATS the summed wall time of every DUST call (not the wall time of
/// the step: see the module doc), printed by `main` at exit, only when at least one batch was
/// masked.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    if timing_on() && DUST_CALLS.load(Ordering::Relaxed) > 0 {
        let seconds = DUST_NANOS.load(Ordering::Relaxed) as f64 / 1e9;
        if std::env::var_os("LOSAT_TIMING").is_some() {
            eprintln!("[TIMING] query_dust: {seconds:.3}s");
        }
        if std::env::var_os("LOSAT_X_STATS").is_some() {
            eprintln!("[X_STATS] QDUSTPAR dust_s={seconds:.3}");
        }
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
    /// with empty sequences and batches of one and of two. (That the masks are made on the pool
    /// is the business of `masks_are_made_on_the_pool_threads`; at 1 thread and for batches of
    /// one this test compares the serial map with itself.)
    #[cfg(feature = "parallel")]
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
            let expect = serial(&refs[..], level, window, linker);
            for threads in [1usize, 2, 3, 7] {
                with_search_pool(threads, "blastn", |pool| {
                    for take in [refs.len(), 2, 1, 0] {
                        let part = &refs[..take];
                        let want = serial(part, level, window, linker);
                        let got = masks_with_mode(part, level, window, linker, pool, true);
                        assert_eq!(got, want, "threads {threads} take {take}");
                    }
                    let got = masks_with_mode(&refs[..], level, window, linker, pool, true);
                    assert_eq!(got, expect, "threads {threads}");
                    Ok(())
                })
                .unwrap();
            }
        }
    }

    /// The parallel path acts only where a search pool of two or more threads is current.
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
            Ok(())
        })
        .unwrap();
    }

    /// With the switch on, the masks are made by the threads of the search pool, and by more than
    /// one of them: the probe that the parallel map calls once per query sees a rayon worker of
    /// a pool of three threads, and the queries are spread over at least two workers. (Each probe
    /// call waits until a second worker has called it (the first call that waits 10 s in vain
    /// ends all waiting), so that the spread does not depend on the scheduler; a serial map has
    /// no second worker, and fails after those 10 s.) With the
    /// switch off, or on a pool of one thread, or for one query, the probe is never called.
    #[cfg(feature = "parallel")]
    #[test]
    fn masks_are_made_on_the_pool_threads() {
        use std::collections::BTreeSet;
        use std::sync::atomic::{AtomicBool, AtomicUsize};
        use std::sync::Mutex;
        use std::time::{Duration, Instant};

        let seqs: Vec<Vec<u8>> = (0..30).map(|i| vec![b"ACGT"[i % 4]; 200]).collect();
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        let want = serial(&refs[..], 20, 64, 1);
        with_search_pool(3, "blastn", |pool| {
            let calls = AtomicUsize::new(0);
            let workers = Mutex::new(BTreeSet::new());
            let in_pool_of_three = AtomicUsize::new(0);
            let gave_up = AtomicBool::new(false);
            let probe = || {
                calls.fetch_add(1, Ordering::SeqCst);
                if let Some(index) = rayon::current_thread_index() {
                    workers.lock().unwrap().insert(index);
                    if rayon::current_num_threads() == 3 {
                        in_pool_of_three.fetch_add(1, Ordering::SeqCst);
                    }
                }
                let start = Instant::now();
                while !gave_up.load(Ordering::SeqCst) && workers.lock().unwrap().len() < 2 {
                    if start.elapsed() > Duration::from_secs(10) {
                        gave_up.store(true, Ordering::SeqCst);
                    }
                    std::thread::yield_now();
                }
            };
            calls.store(0, Ordering::SeqCst);
            in_pool_of_three.store(0, Ordering::SeqCst);
            workers.lock().unwrap().clear();
            gave_up.store(false, Ordering::SeqCst);
            let got = masks_with_probe(&refs[..], 20, 64, 1, pool, true, &probe);
            assert_eq!(got, want);
            assert_eq!(calls.load(Ordering::SeqCst), refs.len());
            assert_eq!(
                in_pool_of_three.load(Ordering::SeqCst),
                refs.len(),
                "every query is masked by a worker of the 3-thread pool"
            );
            {
                let workers = workers.lock().unwrap();
                assert!(
                    workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                    "workers {workers:?}"
                );
            }
            // Not acted on: switch off, one query.
            calls.store(0, Ordering::SeqCst);
            assert_eq!(
                masks_with_probe(&refs[..], 20, 64, 1, pool, false, &probe),
                want
            );
            assert_eq!(
                masks_with_probe(&refs[..1], 20, 64, 1, pool, true, &probe),
                want[..1]
            );
            assert_eq!(calls.load(Ordering::SeqCst), 0);
            Ok(())
        })
        .unwrap();
        with_search_pool(1, "blastn", |pool| {
            let calls = AtomicUsize::new(0);
            let probe = || {
                calls.fetch_add(1, Ordering::SeqCst);
            };
            assert_eq!(
                masks_with_probe(&refs[..], 20, 64, 1, pool, true, &probe),
                want
            );
            assert_eq!(calls.load(Ordering::SeqCst), 0);
            Ok(())
        })
        .unwrap();
    }
}
