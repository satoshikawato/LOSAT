//! Shared BLASTN `-db` handling of HSPs whose Seq-align query start lies below 0 in record
//! coordinates (class N: plus-strand extension left of the query context or minus-strand
//! extension right of it, still below 0 after any `-query_loc` shift).
//!
//! For `-subject` searches, `s_RemapToSubjectLoc` stops the batch before any report runs
//! (`blast_engine::run::left_of_query_refusal`). For `-db` searches NCBI does not remap
//! subjects, so every formatter reads the Seq-align directly and acts according to its own
//! query-row access (`Dense_diag::GetSeqRange` in unsigned `TSeqPos`, `CAlnVec` throwing
//! `CAlnException::eInvalidDenseg` in `x_GetSeqRightSeg` / `x_GetSeqLeftSeg`, `CDense_seg`
//! throwing `CSeqalignException::eInvalidAlignment` in `GetSeqStart` / `GetSeqStop`, or
//! `CAlnMixMerger::x_Merge` throwing `CAlnException::eMergeFailure`).
//!
//! NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_diag.cpp:72-90
//! NCBI reference (598d8ae6): c++/include/objects/seqalign/Dense_diag.hpp:98-121
//! NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:165-238
//! NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:724-760
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-1095

use crate::blastinput::seq_range::Placements;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::hit_sort::HspScores;
use crate::report::tabular::{TabularField, TabularSpec};

use super::hsp::{BlastnHitList, BlastnHsp};

const fn dense_seg_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/objects/seq/../seqalign/Dense_seg.cpp",
        line,
        kind: "CSeqalignException::eInvalidAlignment",
        function: None,
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:197-199
/// `CDense_seg::GetSeqStart(): Row is empty` thrown when BLAST XML / XML2 / JSON reads the
/// query row start of a Dense-seg whose query row has no segment with start >= 0.
const DENSE_SEG_GET_SEQ_START: ThrowSite = dense_seg_site(197);

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:235-237
/// `CDense_seg::GetSeqStop(): Row is empty` thrown when `CDense_seg::GetSeqRange` (GCC
/// evaluates `GetSeqStop` before `GetSeqStart`), `SortHspByMasterStartAscending`,
/// `GetMasterCoverage`, `GetUniqSeqCoverage`, or `x_InitDeflineTable` reads row 0 of a gapped
/// Dense-seg whose query starts are all negative.
const DENSE_SEG_GET_SEQ_STOP: ThrowSite = dense_seg_site(235);

/// `CDense_seg::GetSeqStart(): Row is empty` (outfmt 5, 13, 14, 15, 16).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:197-199
#[cold]
#[inline(never)]
pub(crate) fn dense_seg_empty_row_start() -> NcbiException {
    NcbiException::new(
        DENSE_SEG_GET_SEQ_START,
        "CDense_seg::GetSeqStart(): Row is empty",
    )
}

/// `CDense_seg::GetSeqStop(): Row is empty` (gapped `Dense-seg` in `GetSeqRange`,
/// `SortHspByMasterStartAscending`, `GetMasterCoverage`, `GetUniqSeqCoverage`, and
/// `CShowBlastDefline::x_InitDeflineTable`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:235-237
/// NCBI reference (598d8ae6): c++/include/objects/seqalign/Dense_seg.hpp:234-237
#[cold]
#[inline(never)]
pub(crate) fn dense_seg_empty_row_stop() -> NcbiException {
    NcbiException::new(
        DENSE_SEG_GET_SEQ_STOP,
        "CDense_seg::GetSeqStop(): Row is empty",
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:811-825
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:572-660
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2969-2992,4838-4897
/// The exception thrown by `CBlastFormatUtil::InsertSubjectScores` before printing any tabular
/// row of `hit_list` when a gapped (`Dense-seg`) HSP starts below 0 in record coordinates and
/// `qcovs` or `qcovus` is requested (`CDense_seg::GetSeqRange` throws `Dense_seg.cpp:235`),
/// or `None` when `InsertSubjectScores` completes without throwing.
///
/// Returns `None` without allocating unless at least one HSP of `hit_list` lies outside its
/// query context (`outside_query_context`). When `needs_qcovus` is set, delegates to
/// `tabular::uniq_seq_coverage` so the subject-range subtraction and whole-range throw branch
/// share one owner with `GetUniqSeqCoverage`.
pub(crate) fn gapped_insert_subject_scores_throw(
    hit_list: &BlastnHitList,
    q_offset: usize,
    query_length: usize,
    needs_qcovs: bool,
    needs_qcovus: bool,
) -> Option<NcbiException> {
    if query_length == 0 || (!needs_qcovs && !needs_qcovus) {
        return None;
    }
    if !hit_list
        .hsplist_array
        .iter()
        .flat_map(|list| &list.hsps)
        .any(BlastnHsp::outside_query_context)
    {
        return None;
    }
    gapped_insert_subject_scores_throw_slow(hit_list, q_offset, needs_qcovs, needs_qcovus)
}

#[cold]
#[inline(never)]
fn gapped_insert_subject_scores_throw_slow(
    hit_list: &BlastnHitList,
    q_offset: usize,
    needs_qcovs: bool,
    needs_qcovus: bool,
) -> Option<NcbiException> {
    if needs_qcovs
        && hit_list
            .hsplist_array
            .iter()
            .flat_map(|list| &list.hsps)
            .any(|hsp| hsp.starts_left_of_record(q_offset, false))
    {
        return Some(dense_seg_empty_row_stop());
    }
    if needs_qcovus {
        for hsp_list in &hit_list.hsplist_array {
            if !hsp_list
                .hsps
                .iter()
                .any(|hsp| hsp.starts_left_of_record(q_offset, false))
            {
                continue;
            }
            if let Err(exc) = super::tabular::uniq_seq_coverage(&hsp_list.hsps, q_offset, false) {
                return Some(exc);
            }
        }
    }
    None
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-1095
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:724-760
/// The exception that NCBI's `CBlastTabularInfo::SetFields` throws when formatting a BLASTN
/// `-db` HSP whose Seq-align query start is below 0 in record coordinates, or `None` when the
/// requested fields do not inspect row 0 through `CAlnVec`:
/// - Outside `opens_counts`: returns `None` (no alignment block opened).
/// - Length/gap/mismatch/percent fields (`length`, `gaps`, `mismatch`, `ppos`, `pident`,
///   `gapopen`): `x_SetCounts` calls `GetAlnStop` -> `x_GetSeqRightSeg` (`alnmap.cpp:753`).
/// - Row/query-coord/frame fields (`qseq`, `sseq`, `positive`, `ppos`, `btop`, `qstart`, `qend`,
///   `qframe`, `frames`): calls `GetWholeAlnSeqString` or `GetSeqStart` on row 0, throwing
///   `x_GetSeqLeftSeg` (`alnmap.cpp:732`).
/// - Subject-only fields inside the block (`sstart`, `send`, `sstrand`, `sframe`, `nident`
///   without `frames`): reads row 1 only, returns `None`.
#[cold]
#[inline(never)]
pub(crate) fn tabular_exception(spec: &TabularSpec, num_ident: usize) -> Option<NcbiException> {
    use TabularField::*;
    if !spec.opens_counts(&num_ident) {
        return None;
    }
    if spec.has_any(&[
        AlignmentLength,
        Gaps,
        Mismatches,
        PercentPositives,
        PercentIdentical,
        GapOpenings,
    ]) {
        return Some(crate::report::serial::seqalign::aln_vec_row_gaps_only(
            0, false,
        ));
    }
    if spec.has_any(&[
        QuerySeq, SubjectSeq, Positives, Btop, QueryStart, QueryEnd, QueryFrame, Frames,
    ]) {
        return Some(crate::report::serial::seqalign::aln_vec_row_gaps_only(
            0, true,
        ));
    }
    None
}

#[inline]
fn record_query_from_len(hsp: &BlastnHsp, offset: usize, ungapped: bool) -> (u32, u32) {
    let from = hsp.record_query_start(offset, ungapped) as u32;
    let len = (hsp.internal_q_end_0 - hsp.internal_q_offset_0) as u32;
    (from, len)
}

// NCBI reference (598d8ae6): c++/include/objects/seqalign/Dense_diag.hpp:98-121
// NCBI reference (598d8ae6): c++/include/util/range.hpp:139-175
/// The 0-based `TSeqPos` query range `(from, to)` of `CDense_diag::GetSeqRange` on row 0 for
/// `hsp` after `RemapToQueryLoc` by `offset`, and its length in residues.
#[inline]
pub(crate) fn record_query_u32_range(hsp: &BlastnHsp, offset: usize, ungapped: bool) -> (u32, u32) {
    let (from, len) = record_query_from_len(hsp, offset, ungapped);
    (from, from.wrapping_add(len).wrapping_sub(1))
}

// NCBI reference (598d8ae6): c++/include/objects/seqalign/Dense_diag.hpp:118-121
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:162-169
/// The inclusive 0-based `TSeqPos` query range of `hsp` when non-empty in NCBI's
/// `CRange<TSeqPos>`, or `None` when `from + len` overflows `u32` (wrapping across 0 so
/// `to_open <= from`, which `CRange::Empty()` treats as empty and `CSeq_align::GetSeqRange`
/// drops via `CombineWith`).
#[inline]
pub(crate) fn non_empty_record_query_range(
    hsp: &BlastnHsp,
    offset: usize,
    ungapped: bool,
) -> Option<(u64, u64)> {
    let (from, len) = record_query_from_len(hsp, offset, ungapped);
    let to_open = from.checked_add(len)?;
    (to_open > from).then(|| (u64::from(from), u64::from(to_open - 1)))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1402-1430
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2969-2992
/// `HspScores` of an unshifted BLASTN HSP directly from `BlastnHsp`, matching `HspScores::of` on
/// the `PairwiseHit` that `shift_to_records` + `pairwise::Rows::list` will build for it with
/// query `-query_loc` start `offset` and run mode `ungapped` (including class N `TSeqPos` query
/// ranges, signed sort start, and ungapped negative-start HSPs lifted to `>= 0` by `offset`).
pub(crate) fn blastn_hsp_scores(hsp: &BlastnHsp, offset: usize, ungapped: bool) -> HspScores {
    if hsp.starts_left_of_record(offset, ungapped) {
        let (from, len) = record_query_from_len(hsp, offset, ungapped);
        return HspScores::of_left_of_record(
            hsp.bit_score,
            hsp.e_value,
            hsp.num_ident as i64,
            hsp.length,
            (from as i32, len),
            (hsp.s_start, hsp.s_end),
        );
    }
    let query = if hsp.outside_query_context() && hsp.seqalign_query_start() < 0 {
        let start = hsp.record_query_start(offset, ungapped);
        let len = i64::from(hsp.internal_q_end_0) - i64::from(hsp.internal_q_offset_0);
        ((start + 1).max(0) as usize, (start + len).max(0) as usize)
    } else {
        (hsp.q_start + offset, hsp.q_end + offset)
    };
    HspScores::of_values(
        hsp.bit_score,
        hsp.e_value,
        hsp.num_ident as i64,
        hsp.length,
        query,
        (hsp.s_start, hsp.s_end),
    )
}

/// Whether any HSP in `hit_lists` starts left of its query record (`starts_left_of_record`).
/// Used to route `x_tab_par` to the serial writer when a report will raise an NCBI exception
/// mid-stream.
pub(crate) fn any_left_of_record(
    hit_lists: &[Option<BlastnHitList>],
    placements: &Placements,
    ungapped: bool,
) -> bool {
    hit_lists.iter().enumerate().any(|(q_idx, list)| {
        let Some(list) = list else {
            return false;
        };
        let offset = placements.offset(q_idx);
        list.hsplist_array
            .iter()
            .flat_map(|hsp_list| &hsp_list.hsps)
            .any(|hsp| hsp.starts_left_of_record(offset, ungapped))
    })
}

#[cfg(test)]
#[path = "left_of_query_tests.rs"]
mod tests;
