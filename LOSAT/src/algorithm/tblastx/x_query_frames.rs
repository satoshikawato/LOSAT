//! EXPERIMENT (LOSAT_X_TBXQFRAMES): the SEG of a TBLASTX query
//! batch's frames made on the search pool, one task per frame, each task writing only its own
//! frame.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
//! NCBI filters the contexts of a query batch (for TBLASTX, the six frames of each query) one
//! after the other. `run_impl` did the same, serially, query by query and frame by frame. For a
//! genome-sized query that is most of the time outside the search at eight threads (AP027131 x
//! AP027133, code 4: SEG 0.21 of 0.34 s), and no other thread works meanwhile.
//!
//! With the switch on, and a LOSAT search pool of two or more threads current, the frames of the
//! batch that SEG reads (three or more bytes with the sentinels) are listed in query and frame
//! order and masked by `mask_frame` (`mask_frame_soft` with `-soft_masking true`) in a parallel loop
//! (`into_par_iter().for_each()`) run with
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
//! `mask_frame_soft` (`-soft_masking true`) is the first half of `mask_frame`: it writes only the
//! interval list of its frame, so the same holds for it. So no frame's mask depends on another
//! frame or on the order in which the frames are masked.
//! No floating-point operation changes order: each window's entropy is the same expression on
//! the same counts.
//!
//! Programs that reach it: only TBLASTX (`-subject` and `-db`, every query batch of
//! `search_query_batch`, with either value of `-soft_masking`). The translation of the frames (at most a few ms per batch) stays serial.
//! BLASTX already masks the contexts of a large batch on the pool (at 2+ threads,
//! `blastx/query_setup.rs`) and prepares the chunks of a split query in parallel
//! (`blastx/split.rs`).

use std::sync::OnceLock;

use crate::algorithm::common::masks::{combine_masks, dna_to_protein_masks, lowercase_runs};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::blastinput::fasta_reader::FastaRecord;
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:39,1331,1381-1382
        //
        // LOSAT uses NCBISTDAA encoding where X = 21.
        if !frame.seg_masks.is_empty() {
            if frame.aa_seq_nomask.is_none() {
                frame.aa_seq_nomask = Some(frame.aa_seq.clone());
            }
            const X_MASK_NCBISTDAA: u8 = 21; // NCBI protein mask residue `kProtMask` (21)
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
/// As `mask_frame` with `-soft_masking true` (`mask_at_hash`): only the masked intervals into
/// `seg_masks`, which keep words out of the lookup table (`BLAST_ComplementMaskLocations`,
/// blast_setup.c:634-637); the
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
// No NCBI counterpart: chooses at compile time between the two bodies of blast_setup.c:621-624
// (`mask_frame_soft`, `mask_frame`); it computes nothing NCBI does not.
#[inline(always)]
fn mask_frame_as<const MASK_AT_HASH: bool>(seg: &SegMasker, frame: &mut QueryFrame) {
    if MASK_AT_HASH {
        mask_frame_soft(seg, frame)
    } else {
        mask_frame(seg, frame)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
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
/// As `mask`, for `-soft_masking true` (`mask_at_hash`): the frames get their SEG intervals in
/// `seg_masks` and keep their letters. Cold and out of line (no NCBI counterpart: code placement
/// only), so that the caller, `batch_setup::query_frames`, keeps the default path's code.
#[cold]
#[inline(never)]
pub(super) fn mask_soft(
    query_frames: &mut [Vec<QueryFrame>],
    seg: &SegMasker,
    pool: &SearchPool<'_>,
) {
    mask_with_mode::<true>(query_frames, seg, pool, on())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1239-1256
// The lower-case runs reach each frame in its protein coordinates (`s_AddMask`,
// c++/src/algo/blast/api/blast_setup_cxx.cpp:320-350, `UseProteinCoords`: `dna_to_protein_masks`).
/// `mask` / `mask_soft` for `-lcase_masking`: every frame's SEG intervals (none with `-seg no`) and
/// its query's lower-case runs in `seg_masks`, combined as NCBI combines them. Without `soft` the
/// frame keeps an unmasked copy and its masked letters become X (`mask_frame`). `query_frames[i]`
/// are the frames of `records[i]`. Cold and out of line (no NCBI counterpart: code placement only),
/// so that the caller, `batch_setup::query_frames`, keeps the default path's code.
#[cold]
#[inline(never)]
pub(super) fn mask_lcase(
    query_frames: &mut [Vec<QueryFrame>],
    records: &[FastaRecord],
    seg: Option<&SegMasker>,
    soft: bool,
    pool: &SearchPool<'_>,
) {
    if let Some(seg) = seg {
        mask_with_mode::<true>(query_frames, seg, pool, on());
    }
    for (frames, record) in query_frames.iter_mut().zip(records) {
        let letters = record.seq();
        let runs = lowercase_runs(letters);
        for frame in frames.iter_mut().filter(|frame| frame.aa_seq.len() >= 3) {
            // SEG's intervals first, as `BlastSeqLocAppend` of `lcase_mask_slp` to `filter_out`:
            // the sort of `BlastSeqLocCombine` compares only the left ends, and a one-letter run on
            // the first letter of a SEG interval is merged only in this order (the results' masks).
            let mut masks: Vec<(i32, i32)> = frame
                .seg_masks
                .iter()
                .map(|&(start, end)| (start as i32, end as i32 - 1))
                .collect();
            masks.extend(dna_to_protein_masks(&runs, letters.len(), frame.frame));
            let combined = combine_masks(masks);
            frame.seg_masks = combined
                .iter()
                .map(|&(from, to)| (from as usize, to as usize + 1))
                .collect();
            if !soft && !frame.seg_masks.is_empty() {
                frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                for &(from, to) in &frame.seg_masks {
                    // NCBI `kProtMask` is 21 (X in NCBISTDAA), past the leading sentinel.
                    frame.aa_seq[1 + from..1 + to].fill(21);
                }
            }
        }
    }
}

#[cfg(test)]
#[path = "x_query_frames_tests.rs"]
mod tests;
