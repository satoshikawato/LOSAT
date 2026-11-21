//! The redo-stage SEG of a BLASTX subject computed once per subject instead of once per window.
//!
//! NCBI copies the whole subject and runs SEG on it every time the composition-based redo asks
//! for a subject range, i.e. for every window of every query context of a match:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1638
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1440-1451
//!
//! The SEG intervals are a pure function of the subject bytes (the masker parameters are the
//! constants of `kappa::get_range`). This module keeps, per thread, the bytes of the last
//! subject it masked and their interval list. A call with byte-equal input (the whole slice is
//! compared) returns the stored list; any other call evaluates the reference masker and stores
//! the result. The caller's control flow (when SEG runs, `subject_maybe_biased`, the masking
//! with residue 21, the range fit) is unchanged, so every value the redo computes is the
//! reference value.

use crate::utils::dust::MaskedInterval;
use std::cell::RefCell;

thread_local! {
    /// The last subject this thread masked and its SEG intervals.
    static X_LAST: RefCell<Option<(Vec<u8>, Vec<MaskedInterval>)>> = const { RefCell::new(None) };
}

/// The intervals `reference(subject)` returns, from the per-thread memo when `subject` equals
/// the last masked subject byte for byte. `reference` must be the pure masker call of the caller.
pub(crate) fn x_subject_seg(
    subject: &[u8],
    reference: impl Fn(&[u8]) -> Vec<MaskedInterval>,
) -> Vec<MaskedInterval> {
    X_LAST.with(|last| {
        let mut last = last.borrow_mut();
        if let Some((bytes, intervals)) = last.as_ref() {
            if bytes.as_slice() == subject {
                return intervals.clone();
            }
        }
        let fresh = reference(subject);
        *last = Some((subject.to_vec(), fresh.clone()));
        fresh
    })
}

#[cfg(test)]
#[path = "x_seg_memo_tests.rs"]
mod tests;
