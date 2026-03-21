//! SEG filter implementation
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:1625-1638
// ```c
// if ( !(KAPPA_BLASTP_NO_SEG_SEQUENCE) ) {
//     if (compo_adjust_mode
//         && (!subject_maybe_biased || *subject_maybe_biased)) {
//         ...
//         status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
//                                      subject_maybe_biased);
//     }
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_seg.c:2282-2325
// ```c
// BlastSeqLocFromSeqLoc(..., SegParameters* sparamsp, BlastSeqLoc** seg_locs)
// {
//    ...
//    status = s_SegSeq(seqwin, sparamsp, &segs, 0);
//    ...
//    s_SegsToBlastSeqLoc(segs, offset, seg_locs);
// }
// ```
//
// Kappa subject SEG and query SEG must use the same `blast_seg.c` port.
// Re-export the canonical implementation from `utils::seg` so there is a
// single NCBI-faithful SEG code path in LOSAT.
pub use crate::utils::seg::{SegMasker, SegParams};
