//! EXPERIMENT (LOSAT_X_BNTBRANGE): the fence sentries that this thread's traceback alignments
//! met.
//!
//! No NCBI counterpart: NCBI's kernels set `fence_hit` for `Blast_TracebackFromHSPList`, which
//! breaks its loop, and the caller fetches the whole sequence and traces the HSP list again
//! (blast_traceback.c:516-520, 1651-1688). LOSAT's kernels return as before and count the fence
//! here; a caller compares the count before and after its call (`x_fenced`). Only a sequence
//! decoded in ranges (`database::traceback_listed`) has fences.

use std::cell::Cell;
use std::sync::atomic::{AtomicBool, Ordering};

thread_local! {
    static X_FENCE_HITS: Cell<u64> = const { Cell::new(0) };
}

/// A traceback alignment of this thread stopped at a fence.
#[cold]
#[inline(never)]
pub(crate) fn x_fence_met() {
    X_FENCE_HITS.with(|hits| hits.set(hits.get() + 1));
}

/// The fences this thread's traceback alignments met so far, when `fenced` is the flag of a
/// sequence decoded in ranges (`None` otherwise: nothing to note).
pub(crate) fn x_fence_mark(fenced: Option<&AtomicBool>) -> Option<u64> {
    fenced.map(|_| X_FENCE_HITS.with(Cell::get))
}

/// Sets `fenced` when this thread met a fence since `mark` (`x_fence_mark`): the call between
/// them read a fence of the sequence, and the sequence must be traced again whole.
pub(crate) fn x_fence_note(fenced: Option<&AtomicBool>, mark: Option<u64>) {
    if let (Some(fenced), Some(mark)) = (fenced, mark) {
        if X_FENCE_HITS.with(Cell::get) != mark {
            fenced.store(true, Ordering::Relaxed);
        }
    }
}

/// Whether the sequence of `fenced` met a fence: its traceback stops (`x_fence_note`).
pub(crate) fn x_fenced(fenced: Option<&AtomicBool>) -> bool {
    fenced.is_some_and(|fenced| fenced.load(Ordering::Relaxed))
}
