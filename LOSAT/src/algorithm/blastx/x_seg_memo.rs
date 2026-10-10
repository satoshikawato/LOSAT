//! EXPERIMENT (LOSAT_X_BXSEGMEMO / LOSAT_X_BXSEGMEMOSHADOW): the redo-stage SEG of a BLASTX
//! subject computed once per subject instead of once per window.
//!
//! NCBI copies the whole subject and runs SEG on it every time the composition-based redo asks
//! for a subject range, i.e. for every window of every query context of a match:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1638
//! ```c
//!         if (compo_adjust_mode
//!             && (!subject_maybe_biased || *subject_maybe_biased)) {
//!
//!             if ( (!shouldTestIdentical)
//!                  || (shouldTestIdentical
//!                      && (!s_TestNearIdentical(seqData, 0, queryData,
//!                                               q_range->begin, query_words,
//!                                               align)))) {
//!
//!                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
//!                                              subject_maybe_biased);
//!             }
//!         }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1440-1451
//! ```c
//!         status = BlastSetUp_Filter(program_name, seqData->data,
//!                                    seqData->length, 0, filter_options,
//!                                    &mask_seqloc, NULL);
//!         filter_options = SBlastFilterOptionsFree(filter_options);
//!     }
//!     if (is_seq_biased) {
//!         *is_seq_biased = (mask_seqloc != NULL);
//!     }
//!     if (status == 0) {
//!         Blast_MaskTheResidues(seqData->data, seqData->length,
//!                               FALSE, mask_seqloc, FALSE, 0);
//!     }
//! ```
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

/// EXPERIMENT (LOSAT_X_BXSEGMEMO): 0 = off, 1 = memo, 2 = shadow
/// (`LOSAT_X_BXSEGMEMOSHADOW`: memo and reference on every call, compared).
pub(crate) fn x_bxsegmemo_mode() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BXSEGMEMOSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_BXSEGMEMO").is_some() {
            1
        } else {
            0
        }
    })
}

/// EXPERIMENT (LOSAT_X_BXSEGMEMOSHADOW): memo hits compared with the masker so far.
pub(crate) static X_SHADOW_CALLS: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

thread_local! {
    /// The last subject this thread masked and its SEG intervals.
    static X_LAST: RefCell<Option<(Vec<u8>, Vec<MaskedInterval>)>> = const { RefCell::new(None) };
}

/// EXPERIMENT (LOSAT_X_BXSEGMEMO): the intervals `reference(subject)` returns, from the
/// per-thread memo when `subject` equals the last masked subject byte for byte. `reference`
/// must be the pure masker call of the caller. With the shadow switch the reference is
/// evaluated on every call and compared with the memoised list.
pub(crate) fn x_subject_seg(
    subject: &[u8],
    reference: impl Fn(&[u8]) -> Vec<MaskedInterval>,
) -> Vec<MaskedInterval> {
    let shadow = x_bxsegmemo_mode() == 2;
    X_LAST.with(|last| {
        let mut last = last.borrow_mut();
        if let Some((bytes, intervals)) = last.as_ref() {
            if bytes.as_slice() == subject {
                if shadow {
                    let fresh = reference(subject);
                    assert!(
                        fresh == *intervals,
                        "LOSAT_X_BXSEGMEMOSHADOW: memoised SEG differs on a {}-residue subject",
                        subject.len()
                    );
                    X_SHADOW_CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                }
                return intervals.clone();
            }
        }
        let fresh = reference(subject);
        *last = Some((subject.to_vec(), fresh.clone()));
        fresh
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memo_returns_reference_result_for_equal_bytes_and_recomputes_otherwise() {
        use std::cell::Cell;
        let calls = Cell::new(0usize);
        let reference = |s: &[u8]| {
            calls.set(calls.get() + 1);
            vec![MaskedInterval::new(0, s.len())]
        };
        let a = [1u8, 2, 3, 4];
        let b = [1u8, 2, 3, 5];
        let first = x_subject_seg(&a, reference);
        let again = x_subject_seg(&a, reference);
        assert_eq!(first, again);
        assert_eq!(calls.get(), 1);
        let other = x_subject_seg(&b, reference);
        assert_eq!(other, vec![MaskedInterval::new(0, 4)]);
        assert_eq!(calls.get(), 2);
        // A different length with the same prefix is a different subject.
        let shorter = x_subject_seg(&a[..3], reference);
        assert_eq!(shorter, vec![MaskedInterval::new(0, 3)]);
        assert_eq!(calls.get(), 3);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
    // ```c
    //                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
    //                                              subject_maybe_biased);
    // ```
    #[test]
    fn memo_of_the_redo_masker_equals_the_masker() {
        use crate::utils::matrix::aa_char_to_ncbistdaa;
        use crate::utils::seg::{SegMasker, SegParams};
        let masker = |data: &[u8]| {
            SegMasker::with_params(&SegParams::new(10, 1.8, 2.1))
                .keeping_all_left_segments()
                .mask_sequence(data)
        };
        let subjects: Vec<Vec<u8>> = [
            "MKTAYIAKQRQISFVKSHFSRQAAAAAAAAAAAAAAAAAPPPPPPPGGGGGSLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ",
            "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYG",
        ]
        .iter()
        .map(|s| s.bytes().map(aa_char_to_ncbistdaa).collect())
        .collect();
        assert!(
            !masker(&subjects[0]).is_empty(),
            "the fixture must be masked"
        );
        for subject in subjects
            .iter()
            .chain(subjects.iter())
            .chain(subjects.iter().rev())
        {
            assert_eq!(x_subject_seg(subject, masker), masker(subject));
            assert_eq!(x_subject_seg(subject, masker), masker(subject));
        }
    }
}
