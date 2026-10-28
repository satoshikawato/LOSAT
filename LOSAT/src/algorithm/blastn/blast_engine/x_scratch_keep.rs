//! EXPERIMENT (LOSAT_X_BNSCRATCHKEEP; in the strict set: on unless it is `0`): the working memory
//! of a query batch's subject search (each thread's gap-alignment scratch and subject scratch:
//! the subject's ungapped and gapped hits, the diagonal arrays, the scan's offset pairs, the
//! partition tables of the folded strips, the interval tree, the greedy and DP scratch of the
//! look-ahead and speculative tracebacks, the traceback's BLASTNA letters) kept from one query
//! batch to the next of the same search, instead of made for every batch and freed at its end.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:991-1041
//! ```c
//!    Int4 offset_array_size = GetOffsetArraySize(lookup_wrap);
//! ...
//!    aux_struct->offset_pairs =
//!       (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
//! ```
//!
//! What it changes. A query batch's subject search makes one scratch per thread of its pool (the
//! state of `map_init`, or the one state of a one-thread search) and keeps it from one subject to
//! the next, as NCBI keeps a thread's `aux_struct` and `gap_align`; the batch drops it when it
//! ends. The next batch makes it again and grows every buffer again: at one thread the heap's top
//! (15-21 MiB of these buffers) is trimmed after a batch and faulted in by the next, at 16 threads
//! the slots' arenas are given back with `madvise` and the ungapped hits (mapped from 1 MiB on,
//! `LOSAT_X_BNMAPPED`, up to 32 MiB per batch) are unmapped and mapped again. Here the search
//! holds a `Keep`: a batch borrows its scratches from it (`lend`) and they go back to it when the
//! batch drops them (`Lent`'s `Drop`), so the next batch, of the batch loop or of a
//! LOSAT_X_QBATCHPAR slot, takes them with their buffers resident.
//!
//! Why the values are NCBI's. A borrowed scratch is renewed for its batch first (`SubjectScratch::
//! renew`): every part whose contents a batch reads before it writes them (the hit lists, the
//! diagonal arrays and their offset, the diagonal hash, the offset pairs) is set to what `new`
//! makes, keeping its allocation. The parts that every use writes before it reads them (the
//! interval tree, reset to each subject's bounds; the partition tables, copied from the diagonal
//! hash by every stage; the alignment scratch) are kept as they are, as each thread's scratch
//! keeps them today from one subject to the next. Which scratch a subject gets already depends on
//! how the pool's threads take the subjects, so no value depends on a scratch's history.
//!
//! Memory. The keep holds as many scratches as the search's batches had in use at once (one per
//! thread of the loop's pool, or of each slot's), each with the largest buffers its batches grew.
//! Between two batches they stay resident where they were freed before; within a batch they are
//! what the batch used anyway.
//!
//! Timing: with LOSAT_TIMING, `[TIMING] x_scratch_keep: renewed= fresh=` counts the scratches
//! batches took from a keep and those made new (with a keep or without one).

use std::ops::{Deref, DerefMut};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

static RENEWED: AtomicU64 = AtomicU64::new(0);
static FRESH: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: LOSAT_X_BNSCRATCHKEEP (read once).
fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNSCRATCHKEEP"))
}

/// The scratches one search keeps between its query batches (LOSAT_X_BNSCRATCHKEEP).
pub(super) struct Keep<T> {
    spare: Mutex<Vec<T>>,
}

impl<T> Keep<T> {
    /// A keep for one search, when LOSAT_X_BNSCRATCHKEEP is on.
    pub(super) fn for_search() -> Option<Self> {
        on().then(|| Keep {
            spare: Mutex::new(Vec::new()),
        })
    }
}

/// A scratch a batch uses: borrowed from the search's keep, which takes it back when the batch
/// drops it, or made for the batch alone.
pub(super) struct Lent<'a, T> {
    item: Option<T>,
    keep: Option<&'a Keep<T>>,
}

/// A scratch for a batch: a kept one renewed for it (`renew`) when `keep` has one, else a new one.
pub(super) fn lend<'a, T>(
    keep: Option<&'a Keep<T>>,
    renew: impl FnOnce(T) -> T,
    new: impl FnOnce() -> T,
) -> Lent<'a, T> {
    let kept = keep.and_then(|keep| {
        keep.spare
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .pop()
    });
    let counter = if kept.is_some() { &RENEWED } else { &FRESH };
    counter.fetch_add(1, Ordering::Relaxed);
    Lent {
        item: Some(match kept {
            Some(item) => renew(item),
            None => new(),
        }),
        keep,
    }
}

/// No NCBI counterpart: with LOSAT_TIMING, the `[TIMING] x_scratch_keep` line `main` prints at
/// exit (when a batch made or took a scratch).
pub(super) fn print_summary() {
    let renewed = RENEWED.load(Ordering::Relaxed);
    let fresh = FRESH.load(Ordering::Relaxed);
    if crate::utils::xtiming::on() && renewed + fresh > 0 {
        eprintln!("[TIMING] x_scratch_keep: renewed={renewed} fresh={fresh}");
    }
}

impl<T> Deref for Lent<'_, T> {
    type Target = T;
    fn deref(&self) -> &T {
        self.item.as_ref().expect("a lent scratch until its drop")
    }
}

impl<T> DerefMut for Lent<'_, T> {
    fn deref_mut(&mut self) -> &mut T {
        self.item.as_mut().expect("a lent scratch until its drop")
    }
}

impl<T> Drop for Lent<'_, T> {
    fn drop(&mut self) {
        if let (Some(keep), Some(item)) = (self.keep, self.item.take()) {
            // A panicking batch's scratch is not taken back (its state is not trusted).
            if !std::thread::panicking() {
                keep.spare
                    .lock()
                    .unwrap_or_else(PoisonError::into_inner)
                    .push(item);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_dropped_scratch_is_renewed_for_the_next_batch() {
        let keep = Keep {
            spare: Mutex::new(Vec::new()),
        };
        {
            let mut a = lend(Some(&keep), |_: Vec<u8>| unreachable!(), || vec![1, 2, 3]);
            a.push(4);
        }
        let b = lend(
            Some(&keep),
            |mut v: Vec<u8>| {
                assert_eq!(v, [1, 2, 3, 4]);
                v.clear();
                v
            },
            || unreachable!(),
        );
        assert!(b.is_empty() && b.capacity() >= 4);
        let c = lend(Some(&keep), |_| unreachable!(), Vec::new);
        assert!(c.is_empty());
        drop((b, c));
        assert_eq!(keep.spare.lock().unwrap().len(), 2);
        let d = lend(None, |_: Vec<u8>| unreachable!(), || vec![7]);
        assert_eq!(*d, [7]);
    }
}
