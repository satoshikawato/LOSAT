//! EXPERIMENT (LOSAT_X_TBXQFRAMES): the SEG of a TBLASTX query
//! batch's frames made on the search pool, one task per frame, each task writing only its own
//! frame.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
//! ```c
//!     for (context = query_info->first_context;
//!          context <= query_info->last_context; ++context) {
//!
//!         BlastSeqLoc *filter_per_context = NULL;
//!         status = s_GetFilteringLocationsForOneContext(query_blk,
//!                                                       query_info,
//!                                                       context,
//!                                                       program_number,
//!                                                       filter_options,
//!                                                       &filter_per_context,
//!                                                       blast_message);
//! ```
//! NCBI filters the contexts of a query batch (for TBLASTX, the six frames of each query) one
//! after the other. `run_impl` did the same, serially, query by query and frame by frame. For a
//! genome-sized query that is most of the time outside the search at eight threads (AP027131 x
//! AP027133, code 4: SEG 0.21 of 0.34 s), and no other thread works meanwhile.
//!
//! With the switch on, and a LOSAT search pool of two or more threads current, the frames of the
//! batch that SEG reads (three or more bytes with the sentinels) are listed in query and frame
//! order and masked by `mask_frame` in a parallel loop (`into_par_iter().for_each()`) run with
//! `SearchPool::install`, so on the pool of the search that calls it. Rayon's global pool is never
//! used. A search without a pool (one thread, a build without the `parallel` feature, WebAssembly
//! without `wasm-threads`) and a batch with fewer than two such frames run the serial loop.
//! Without SEG (`-seg no`) there is no loop at all.
//!
//! Why the frames are the same: `mask_frame` reads and writes only the frame it is given (its
//! letters, its interval list and its unmasked copy) and the `SegMasker`, which it only reads (the
//! window, the cut-offs, the alphabet; no interior mutability). SEG keeps two caches, the
//! entropy memo of `entropy_from_state_vector` and the histogram memo of the fast SEG path, both
//! `thread_local!`, so each pool thread has its own; both store the value of a pure function of
//! their key (the sorted letter counts of a window), computed by the same reference function on a
//! miss, so a frame masked on any thread, with a warm or a cold cache, gets the same intervals as in
//! the serial loop. The statics of `utils/seg.rs` are those switches, read once and never written.
//! So no frame's mask depends on another frame or on the order in which the frames are masked.
//! No floating-point operation changes order: each window's entropy is the same expression on
//! the same counts.
//!
//! Programs that reach it: only TBLASTX (`-subject` and `-db`, every query batch of
//! `search_query_batch`). The translation of the frames (at most a few ms per batch) stays serial.
//! BLASTX already masks the contexts of a large batch on the pool (at 2+ threads,
//! `blastx/query_setup.rs`) and prepares the chunks of a split query in parallel
//! (`blastx/split.rs`).

use std::sync::OnceLock;

use crate::algorithm::tblastx::translation::QueryFrame;
use crate::utils::seg::SegMasker;
use crate::utils::threading::SearchPool;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

// No NCBI counterpart: reads the switch once (LOSAT_X_TBXQFRAMES is a strict-set switch: on unless
// it is `0`). It chooses between two LOSAT code paths with the same result and changes no value
// NCBI computes.
/// LOSAT_X_TBXQFRAMES.
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXQFRAMES"))
}

/// SEG of one query frame: the masked intervals into `seg_masks`, then the masked letters
/// overwritten with X, after keeping an unmasked copy. A frame of fewer than three bytes (the two
/// sentinels and at most one letter) is left as it is.
pub(super) fn mask_frame(seg: &SegMasker, frame: &mut QueryFrame) {
    if frame.aa_seq.len() >= 3 {
        for m in seg.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
            frame.seg_masks.push((m.start, m.end));
        }

        // NCBI BLAST query masking semantics (SEG, etc.): keep an unmasked copy
        // (`sequence_nomask`), then overwrite masked residues in the working
        // query sequence buffer with X.
        //
        // NCBI reference (verbatim):
        //   const Uint1 kProtMask = 21;     /* X in NCBISTDAA */
        //   query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
        //   query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
        //   buffer[index] = kMaskingLetter;
        //
        // LOSAT uses NCBISTDAA encoding where X = 21.
        if !frame.seg_masks.is_empty() {
            if frame.aa_seq_nomask.is_none() {
                frame.aa_seq_nomask = Some(frame.aa_seq.clone());
            }
            const X_MASK_NCBISTDAA: u8 = 21; // NCBI: kProtMask = 21
            let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1); // keep last sentinel untouched
            for &(s, e) in &frame.seg_masks {
                let raw_s = 1usize.saturating_add(s);
                let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                for pos in raw_s..raw_e {
                    frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                }
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
// ```c
//     if (!mask_at_hash) {
//         BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
//                              program_number);
//     }
// ```
/// As `mask_frame` with `-soft_masking true` (`mask_at_hash`): only the masked intervals into
/// `seg_masks`, which keep words out of the lookup table (`BLAST_ComplementMaskLocations`); the
/// frame keeps its letters and has no unmasked copy, so the extensions and the frame's Karlin block
/// read the translation.
fn mask_frame_soft(seg: &SegMasker, frame: &mut QueryFrame) {
    if frame.aa_seq.len() >= 3 {
        for m in seg.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
            frame.seg_masks.push((m.start, m.end));
        }
    }
}

/// `mask_frame_soft` with `MASK_AT_HASH`, `mask_frame` without. `MASK_AT_HASH` is a const parameter
/// so that the default (hard) instance of the loops below keeps the code it had before
/// `-soft_masking` (one instance per value).
#[inline(always)]
fn mask_frame_as<const MASK_AT_HASH: bool>(seg: &SegMasker, frame: &mut QueryFrame) {
    if MASK_AT_HASH {
        mask_frame_soft(seg, frame)
    } else {
        mask_frame(seg, frame)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
// ```c
//     for (context = query_info->first_context;
//          context <= query_info->last_context; ++context) {
// ```
/// The serial loop: every frame of every query, in query and frame order.
fn serial<const MASK_AT_HASH: bool>(query_frames: &mut [Vec<QueryFrame>], seg: &SegMasker) {
    for frames in query_frames.iter_mut() {
        for frame in frames.iter_mut() {
            mask_frame_as::<MASK_AT_HASH>(seg, frame);
        }
    }
}

// No NCBI counterpart: the number of frames SEG reads (three or more bytes), which decides
// between the pool and the serial loop; it computes no value NCBI computes.
fn masked_frame_count(query_frames: &[Vec<QueryFrame>]) -> usize {
    query_frames
        .iter()
        .flatten()
        .filter(|frame| frame.aa_seq.len() >= 3)
        .count()
}

/// Whether the frames are masked on `pool`: two or more frames that SEG reads and a LOSAT search
/// pool of two or more threads.
// No NCBI counterpart: scheduling only.
fn on_pool(pool: &SearchPool<'_>, frame_count: usize) -> bool {
    frame_count > 1 && pool.enabled() && pool.threads() > 1
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290 (the loop of the
// module header). The same `mask_frame` per frame, the frames spread over the pool's threads;
// each task writes only its own frame.
/// The parallel loop on `pool`; the same frames as `serial`. Only called when `on_pool`. `probe`
/// is called once per frame that SEG reads, on the thread that masks it; production passes
/// `|| {}`, which the compiler removes (a unit test passes a closure that looks at the thread).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn parallel<const MASK_AT_HASH: bool, P: Fn() + Sync>(
    query_frames: &mut [Vec<QueryFrame>],
    seg: &SegMasker,
    pool: &SearchPool<'_>,
    probe: &P,
) {
    let frames: Vec<&mut QueryFrame> = query_frames
        .iter_mut()
        .flatten()
        .filter(|frame| frame.aa_seq.len() >= 3)
        .collect();
    pool.install(|| {
        frames.into_par_iter().for_each(|frame| {
            probe();
            mask_frame_as::<MASK_AT_HASH>(seg, frame);
        })
    });
}

// No NCBI counterpart: a build without threads runs the serial loop.
#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn parallel<const MASK_AT_HASH: bool, P: Fn() + Sync>(
    query_frames: &mut [Vec<QueryFrame>],
    seg: &SegMasker,
    _pool: &SearchPool<'_>,
    _probe: &P,
) {
    serial::<MASK_AT_HASH>(query_frames, seg);
}

/// SEG of every frame of `query_frames` with the switch `on` or off. Separate from the environment
/// so that tests choose the switch.
pub(super) fn mask_with_mode<const MASK_AT_HASH: bool>(
    query_frames: &mut [Vec<QueryFrame>],
    seg: &SegMasker,
    pool: &SearchPool<'_>,
    on: bool,
) {
    mask_with_probe::<MASK_AT_HASH, _>(query_frames, seg, pool, on, &|| {})
}

/// `mask_with_mode` with a `probe` that the parallel loop calls once per frame (see `parallel`).
// No NCBI counterpart: the dispatch between the serial loop and the pool; it computes no value
// NCBI computes.
fn mask_with_probe<const MASK_AT_HASH: bool, P: Fn() + Sync>(
    query_frames: &mut [Vec<QueryFrame>],
    seg: &SegMasker,
    pool: &SearchPool<'_>,
    on: bool,
    probe: &P,
) {
    let frame_count = masked_frame_count(query_frames);
    let pooled = on_pool(pool, frame_count);
    if on {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXQFRAMES", pooled);
    }
    if !on || !pooled {
        serial::<MASK_AT_HASH>(query_frames, seg);
        return;
    }
    parallel::<MASK_AT_HASH, P>(query_frames, seg, pool, probe);
}

/// SEG of every frame of `query_frames` as `run_impl` makes it: the environment's switch decides.
pub(super) fn mask(query_frames: &mut [Vec<QueryFrame>], seg: &SegMasker, pool: &SearchPool<'_>) {
    mask_with_mode::<false>(query_frames, seg, pool, on())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
// ```c
//     if (!mask_at_hash) {
//         BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
//                              program_number);
//     }
// ```
/// As `mask`, for `-soft_masking true` (`mask_at_hash`): the frames get their SEG intervals in
/// `seg_masks` and keep their letters. Cold and out of line (no NCBI counterpart: code placement
/// only), so that the caller, `run_impl::search_query_batch`, keeps the default path's code.
#[cold]
#[inline(never)]
pub(super) fn mask_soft(
    query_frames: &mut [Vec<QueryFrame>],
    seg: &SegMasker,
    pool: &SearchPool<'_>,
) {
    mask_with_mode::<true>(query_frames, seg, pool, on())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastx::translation::generate_frames;
    use crate::utils::genetic_code::GeneticCode;
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

    /// A random DNA sequence with low-complexity runs (so that SEG masks some of every frame),
    /// `N`s and other IUPAC letters; short ones too (no frame, or frames of zero or one letter).
    fn sequence(rng: &mut Rng) -> Vec<u8> {
        let len = match rng.below(6) {
            0 => rng.below(8) as usize,
            1 => 8 + rng.below(60) as usize,
            _ => 300 + rng.below(4000) as usize,
        };
        let mut seq = Vec::with_capacity(len);
        while seq.len() < len {
            match rng.below(8) {
                0 => {
                    let unit = [b"CAG".as_slice(), b"GCA", b"AAAAAG", b"CCTCCA", b"GGT"]
                        [rng.below(5) as usize];
                    for _ in 0..(5 + rng.below(40)) {
                        seq.extend_from_slice(unit);
                    }
                }
                1 => seq.extend(std::iter::repeat_n(b'N', 1 + rng.below(9) as usize)),
                2 => seq.push(b"RYKMSWnacgt"[rng.below(11) as usize]),
                _ => seq.push(b"ACGT"[rng.below(4) as usize]),
            }
        }
        seq.truncate(len);
        seq
    }

    fn frames_of(seqs: &[Vec<u8>]) -> Vec<Vec<QueryFrame>> {
        let code = GeneticCode::from_id(1);
        seqs.iter().map(|s| generate_frames(s, &code)).collect()
    }

    /// The fields SEG writes, per frame, for comparison.
    fn masked(frames: &[Vec<QueryFrame>]) -> Vec<(Vec<u8>, Option<Vec<u8>>, Vec<(usize, usize)>)> {
        frames
            .iter()
            .flatten()
            .map(|f| {
                (
                    f.aa_seq.clone(),
                    f.aa_seq_nomask.clone(),
                    f.seg_masks.clone(),
                )
            })
            .collect()
    }

    /// Pool == serial on random queries with low-complexity runs at 1, 2, 3 and 7 threads, for
    /// several SEG parameter sets, for batches of many queries, of one query and of none. (That the
    /// frames are masked on the pool is the business of `frames_are_masked_on_the_pool_threads`.)
    #[cfg(feature = "parallel")]
    #[test]
    fn pool_frames_equal_serial_frames() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let mut seqs: Vec<Vec<u8>> = (0..40).map(|_| sequence(&mut rng)).collect();
        seqs.push(Vec::new());
        seqs.push(b"ACG".to_vec());
        seqs.push(b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG".repeat(4));
        for (window, locut, hicut) in [(12, 2.2, 2.5), (10, 1.8, 2.1), (20, 3.0, 3.3)] {
            let seg = SegMasker::new(window, locut, hicut);
            let mut want_all = frames_of(&seqs);
            serial::<false>(&mut want_all, &seg);
            let intervals: usize = want_all.iter().flatten().map(|f| f.seg_masks.len()).sum();
            assert!(
                intervals > 40,
                "the frames must have masked regions ({intervals})"
            );
            for threads in [1usize, 2, 3, 7] {
                with_search_pool(threads, "tblastx", |pool| {
                    for take in [seqs.len(), 2, 1, 0] {
                        let mut want = frames_of(&seqs[..take]);
                        serial::<false>(&mut want, &seg);
                        let mut got = frames_of(&seqs[..take]);
                        mask_with_mode::<false>(&mut got, &seg, pool, true);
                        assert_eq!(
                            masked(&got),
                            masked(&want),
                            "window {window} threads {threads} take {take}"
                        );
                    }
                    let mut got = frames_of(&seqs);
                    mask_with_mode::<false>(&mut got, &seg, pool, true);
                    assert_eq!(masked(&got), masked(&want_all), "threads {threads}");
                    Ok(())
                })
                .unwrap();
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
    // ```c
    //     if (!mask_at_hash) {
    //         BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
    //                              program_number);
    //     }
    // ```
    /// `-soft_masking true`: every frame gets the intervals of the hard frame and keeps its letters
    /// (no unmasked copy), in the serial loop and on pools of 1, 2 and 3 threads.
    #[test]
    fn soft_frames_keep_letters_and_hard_intervals() {
        let mut rng = Rng(0x2545_f491_4f6c_dd1d);
        let seqs: Vec<Vec<u8>> = (0..12).map(|_| sequence(&mut rng)).collect();
        let seg = SegMasker::new(12, 2.2, 2.5);
        let plain = frames_of(&seqs);
        let mut hard = frames_of(&seqs);
        serial::<false>(&mut hard, &seg);
        let mut soft = frames_of(&seqs);
        serial::<true>(&mut soft, &seg);
        let intervals: usize = hard.iter().flatten().map(|f| f.seg_masks.len()).sum();
        assert!(intervals > 0, "the frames must have masked regions");
        let frames = soft.iter().flatten().zip(hard.iter().flatten());
        for ((soft, hard), plain) in frames.zip(plain.iter().flatten()) {
            assert_eq!(soft.seg_masks, hard.seg_masks);
            assert_eq!(soft.aa_seq, plain.aa_seq);
            assert!(soft.aa_seq_nomask.is_none());
            assert_eq!(
                hard.aa_seq_nomask.as_ref().unwrap_or(&hard.aa_seq),
                &plain.aa_seq
            );
        }
        for threads in [1usize, 2, 3] {
            with_search_pool(threads, "tblastx", |pool| {
                let mut got = frames_of(&seqs);
                mask_with_mode::<true>(&mut got, &seg, pool, true);
                assert_eq!(masked(&got), masked(&soft), "threads {threads}");
                Ok(())
            })
            .unwrap();
        }
    }

    /// The parallel path acts only where a search pool of two or more threads is current and the
    /// batch has two or more frames that SEG reads.
    #[cfg(feature = "parallel")]
    #[test]
    fn only_a_search_pool_of_two_threads_acts() {
        with_search_pool(1, "tblastx", |pool| {
            assert!(!on_pool(pool, 6));
            Ok(())
        })
        .unwrap();
        with_search_pool(3, "tblastx", |pool| {
            assert!(on_pool(pool, 6));
            assert!(!on_pool(pool, 1));
            assert!(!on_pool(pool, 0));
            Ok(())
        })
        .unwrap();
        assert_eq!(masked_frame_count(&frames_of(&[b"ACGT".to_vec()])), 4);
        assert_eq!(masked_frame_count(&frames_of(&[b"ACGTAC".to_vec()])), 6);
        assert_eq!(masked_frame_count(&frames_of(&[b"AC".to_vec()])), 0);
    }

    /// With the switch on, the frames are masked by the threads of the search pool, and by more
    /// than one of them: the probe that the parallel loop calls once per frame sees a rayon worker
    /// of a pool of three threads, and the frames are spread over at least two workers. (Each probe
    /// call waits until a second worker has called it (the first call that waits 10 s in vain ends
    /// all waiting), so that the spread does not depend on the scheduler; a serial loop has no
    /// second worker, and fails after those 10 s.) With the switch off, or on a pool of one
    /// thread, the probe is never called.
    #[cfg(feature = "parallel")]
    #[test]
    fn frames_are_masked_on_the_pool_threads() {
        use std::collections::BTreeSet;
        use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
        use std::sync::Mutex;
        use std::time::{Duration, Instant};

        let seqs: Vec<Vec<u8>> = (0..5).map(|i| b"CAGCAGTTGACCA".repeat(20 + i)).collect();
        let seg = SegMasker::new(12, 2.2, 2.5);
        let mut want = frames_of(&seqs);
        serial::<false>(&mut want, &seg);
        let frame_count = masked_frame_count(&want);
        with_search_pool(3, "tblastx", |pool| {
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
            // The block drops the `workers` guard (it holds the lock the probe takes) before the
            // next run.
            {
                let mut got = frames_of(&seqs);
                mask_with_probe::<false, _>(&mut got, &seg, pool, true, &probe);
                assert_eq!(masked(&got), masked(&want));
                assert_eq!(calls.load(Ordering::SeqCst), frame_count);
                assert_eq!(
                    in_pool_of_three.load(Ordering::SeqCst),
                    frame_count,
                    "every frame is masked by a worker of the 3-thread pool"
                );
                let workers = workers.lock().unwrap();
                assert!(
                    workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                    "workers {workers:?}"
                );
            }
            // Not acted on: switch off.
            calls.store(0, Ordering::SeqCst);
            let mut got = frames_of(&seqs);
            mask_with_probe::<false, _>(&mut got, &seg, pool, false, &probe);
            assert_eq!(masked(&got), masked(&want));
            assert_eq!(calls.load(Ordering::SeqCst), 0);
            Ok(())
        })
        .unwrap();
        with_search_pool(1, "tblastx", |pool| {
            let calls = AtomicUsize::new(0);
            let probe = || {
                calls.fetch_add(1, Ordering::SeqCst);
            };
            let mut got = frames_of(&seqs);
            mask_with_probe::<false, _>(&mut got, &seg, pool, true, &probe);
            assert_eq!(masked(&got), masked(&want));
            assert_eq!(calls.load(Ordering::SeqCst), 0);
            Ok(())
        })
        .unwrap();
    }
}
