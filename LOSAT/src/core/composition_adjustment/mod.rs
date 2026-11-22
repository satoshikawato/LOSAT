//! NCBI composition-adjustment support.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/

pub mod adjust_scores;
// The core of the position-based composition-based statistics (mode 1) of a PSSM (line PS,
// item PS4b).
pub(crate) mod compo_scale;
pub mod matrix_tables;
pub mod redo_alignment;
// NCBI's c++/src/algo/blast/composition_adjustment/smith_waterman.c: the kernels of the Smith-Waterman redo (`-use_sw_tback`; series PR,
// item PS11c-3).
pub(crate) mod smith_waterman;
pub(crate) mod smith_waterman_match;
// No NCBI counterpart: `LOSAT_X_REDOCONTAINIDX`, the redo's kept alignments with an index by
// subject start (`s_IsContained` / `s_WithDistinctEnds` answer from the alignments that can meet).
pub(crate) mod x_align_index;
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-758,762-768
// `x_newton_exact` reproduces `Blast_OptimizeTargetFrequencies` (optimize_target_freq.c:686-815) with
// the same floating-point operations on every value, in the same order. It solves the problems with
// 20 letters and the relative-entropy constraint; see the module comment.
pub(crate) mod x_newton_exact;
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// `x_newton_lanes` solves several independent calls of this function at once (one per vector lane,
// the operations of `x_newton_exact` per lane) ahead of the calls, and hands a result to a call whose
// input has the same bits. The callers that know their coming inputs (BLASTP's redo loops, TBLASTN's
// window loop) decide whether to prefetch; see the module comment.
pub(crate) mod x_newton_lanes;
