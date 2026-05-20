use anyhow::Result;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use crate::report::{
    write_hit_fields, write_outfmt7, write_pairwise_simple, OutputConfig, OutputFormat,
    PairwiseConfig, ReportContext,
};

// =============================================================================
// Gap Edit Script for traceback trimming
// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c
// =============================================================================

/// Gap edit operation type for traceback
/// Reference: blast_hits.c eGapAlignOpType
///
/// ```c
/// // NCBI GapEditScript structure (gapinfo.h)
/// typedef enum {
///     eGapAlignSub = 0,  // Substitution (match or mismatch)
///     eGapAlignDel = 1,  // Deletion from query (gap in query)
///     eGapAlignIns = 2   // Insertion to query (gap in subject)
/// } EGapAlignOpType;
///
/// typedef struct GapEditScript {
///     Uint1 *op_type;  // Array of operation types
///     Int4 *num;       // Array of counts for each operation
///     Int4 size;       // Number of operations
/// } GapEditScript;
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapEditOp {
    /// Substitution (match or mismatch): consumes `num` query and `num` subject positions
    /// Reference: blast_hits.c:2404-2407
    /// ```c
    /// if (esp->op_type[index] == eGapAlignSub) {
    ///    qid++;
    ///    sid++;
    /// ```
    Sub(u32),
    /// Deletion: gap in query, consumes `num` subject positions only
    /// Reference: blast_hits.c:2408-2410
    /// ```c
    /// } else if (esp->op_type[index] == eGapAlignDel) {
    ///    sid+=esp->num[index];
    /// ```
    Del(u32),
    /// Insertion: gap in subject, consumes `num` query positions only
    /// Reference: blast_hits.c:2411-2413
    /// ```c
    /// } else if (esp->op_type[index] == eGapAlignIns) {
    ///    qid+=esp->num[index];
    /// ```
    Ins(u32),
}

impl GapEditOp {
    /// Get the count for this operation
    #[inline]
    pub fn num(&self) -> u32 {
        match self {
            GapEditOp::Sub(n) | GapEditOp::Del(n) | GapEditOp::Ins(n) => *n,
        }
    }

    /// Check if this is a substitution operation
    #[inline]
    pub fn is_sub(&self) -> bool {
        matches!(self, GapEditOp::Sub(_))
    }
}

#[derive(Debug, Clone)]
pub struct Hit {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    pub identity: f64,
    pub length: usize,
    pub mismatch: usize,
    pub gapopen: usize,
    pub q_start: usize,
    pub q_end: usize,
    pub s_start: usize,
    pub s_end: usize,
    pub e_value: f64,
    pub bit_score: f64,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;           /**< This HSP's raw score */
    //    Int4 num_ident;       /**< Number of identical base pairs in this HSP */
    //    double bit_score;     /**< Bit score, calculated from score */
    //    double evalue;        /**< This HSP's e-value */
    //    BlastSeg query;       /**< Query sequence info. */
    //    BlastSeg subject;     /**< Subject sequence info. */
    //    Int4     context;     /**< Context number of query */
    //    GapEditScript* gap_info;/**< ALL gapped alignment is here */
    //    ...
    //    Int4 num_positives;
    // } BlastHSP;
    // ```
    pub num_ident: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
    // ```c
    // if (hsp->query.frame != hsp->subject.frame) {
    //    *q_end = query_length - hsp->query.offset;
    //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
    //    *s_end = hsp->subject.offset + 1;
    //    *s_start = hsp->subject.end;
    // }
    // ```
    pub query_frame: i32,
    pub query_length: usize,
    // Fields for NCBI-style output ordering (not printed, used for sorting)
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    ...
    // } BlastHSPList;
    // ```
    /// Query index (input order) - NCBI uses query order for grouping
    pub q_idx: u32,
    /// Subject index (oid) - NCBI uses oid for tie-breaking in s_EvalueCompareHSPLists
    pub s_idx: u32,
    /// Raw alignment score - NCBI ScoreCompareHSPs uses raw score, not bit score
    pub raw_score: i32,
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
    // ```c
    // typedef struct BlastHSP {
    //    BlastSeg query;       /**< Query sequence info. */
    //    BlastSeg subject;     /**< Subject sequence info. */
    // } BlastHSP;
    // ```
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2384-2388
    // ```c
    // init_hsp->offsets.qs_offsets.q_off -= query_start;
    // if (init_hsp->ungapped_data) {
    //     init_hsp->ungapped_data->q_start -= query_start;
    // }
    // ```
    // Optional context/frame-relative BlastSeg offsets used by ScoreCompareHSPs
    // when formatted nucleotide coordinates cannot recover the translated HSP
    // coordinates.
    pub sort_query_offset: usize,
    pub sort_query_end: usize,
    pub sort_subject_offset: usize,
    pub sort_subject_end: usize,
    pub has_sort_offsets: bool,
    /// Edit script for traceback trimming (NCBI: gap_info in BlastHSP)
    /// Reference: blast_hits.c BlastHSP.gap_info (GapEditScript*)
    ///
    /// This is used by Blast_HSPListPurgeHSPsWithCommonEndpoints to trim
    /// overlapping HSPs instead of deleting them entirely.
    /// Reference: blast_hits.c:2490-2492, 2516-2518
    pub gap_info: Option<Vec<GapEditOp>>,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
    // ```c
    // typedef struct BlastHSP {
    //    ...
    //    Int4 num_positives;
    //
    //    BlastHSPMappingInfo* map_info;
    // } BlastHSP;
    // ```
    pub num_positives: usize,
}

impl Hit {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    // } BlastHSPList;
    // ```
    #[inline]
    pub fn resolve_ids<'a>(
        &self,
        query_ids: &'a [Arc<str>],
        subject_ids: &'a [Arc<str>],
    ) -> (&'a str, &'a str) {
        let query_id = query_ids
            .get(self.q_idx as usize)
            .map(|id| id.as_ref())
            .unwrap_or("unknown");
        let subject_id = subject_ids
            .get(self.s_idx as usize)
            .map(|id| id.as_ref())
            .unwrap_or("unknown");
        (query_id, subject_id)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
    // ```c
    // s_Blast_HSPGetNumIdentitiesAndPositives(..., Int4* num_ident_ptr, ...,
    //                                         Int4* num_pos_ptr)
    // {
    //    ...
    //    if (align_length_ptr) {
    //        *align_length_ptr = align_length;
    //    }
    //    *num_ident_ptr = num_ident;
    //
    //    if(NULL != matrix)
    //        *num_pos_ptr = num_pos + num_ident;
    // }
    // ```
    #[inline]
    pub fn gap_letters(&self) -> usize {
        self.length
            .saturating_sub(self.num_ident.saturating_add(self.mismatch))
    }
}

// =============================================================================
// NCBI-style comparators for outfmt6 output ordering
// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c
// =============================================================================

/// Compare two evalues, treating both as equal if they're close enough to zero.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c s_EvalueComp()
/// ```c
/// static int s_EvalueComp(double evalue1, double evalue2) {
///     const double epsilon = 1.0e-180;
///     if (evalue1 < epsilon && evalue2 < epsilon) { return 0; }
///     if (evalue1 < evalue2) { return -1; }
///     else if (evalue1 > evalue2) { return 1; }
///     else { return 0; }
/// }
/// ```
#[inline]
fn evalue_comp(evalue1: f64, evalue2: f64) -> Ordering {
    const EPSILON: f64 = 1.0e-180;
    if evalue1 < EPSILON && evalue2 < EPSILON {
        Ordering::Equal
    } else if evalue1 < evalue2 {
        Ordering::Less
    } else if evalue1 > evalue2 {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

/// NCBI ScoreCompareHSPs - compare two HSPs by score and coordinates.
/// Used to sort HSPs within a subject (BlastHSPList).
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c ScoreCompareHSPs()
/// ```c
/// int ScoreCompareHSPs(const void* h1, const void* h2) {
///    BlastHSP* hsp1,* hsp2;
///    int result = 0;
///    hsp1 = *((BlastHSP**) h1);
///    hsp2 = *((BlastHSP**) h2);
///    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
///        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
///        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
///        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
///        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
///    }
///    return result;
/// }
/// ```
/// Order: score DESC → s_start ASC → s_end DESC → q_start ASC → q_end DESC
pub fn score_compare_hsps(a: &Hit, b: &Hit) -> Ordering {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
    // ```c
    // if (hsp->query.frame != hsp->subject.frame) {
    //    *q_end = query_length - hsp->query.offset;
    //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
    // }
    // ```
    // Recover internal query offsets from output coordinates for blastn minus strand.
    let query_offsets = |hit: &Hit| {
        if hit.query_length > 0 && hit.query_frame < 0 {
            let q_offset = hit.query_length.saturating_sub(hit.q_end);
            let q_end = hit
                .query_length
                .saturating_sub(hit.q_start)
                .saturating_add(1);
            (q_offset, q_end)
        } else {
            (hit.q_start.saturating_sub(1), hit.q_end)
        }
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1347-1353
    // ```c
    // if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
    //     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
    //     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
    //     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
    //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
    // }
    // ```
    let (a_q_offset, a_q_end, a_s_offset, a_s_end) = if a.has_sort_offsets {
        (
            a.sort_query_offset,
            a.sort_query_end,
            a.sort_subject_offset,
            a.sort_subject_end,
        )
    } else {
        let (q_offset, q_end) = query_offsets(a);
        (
            q_offset,
            q_end,
            a.s_start.min(a.s_end).saturating_sub(1),
            a.s_start.max(a.s_end),
        )
    };
    let (b_q_offset, b_q_end, b_s_offset, b_s_end) = if b.has_sort_offsets {
        (
            b.sort_query_offset,
            b.sort_query_end,
            b.sort_subject_offset,
            b.sort_subject_end,
        )
    } else {
        let (q_offset, q_end) = query_offsets(b);
        (
            q_offset,
            q_end,
            b.s_start.min(b.s_end).saturating_sub(1),
            b.s_start.max(b.s_end),
        )
    };

    // score DESC (BLAST_CMP(hsp2->score, hsp1->score))
    match b.raw_score.cmp(&a.raw_score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // s_start ASC (BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset))
    match a_s_offset.cmp(&b_s_offset) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // s_end DESC (BLAST_CMP(hsp2->subject.end, hsp1->subject.end))
    match b_s_end.cmp(&a_s_end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // q_start ASC (BLAST_CMP(hsp1->query.offset, hsp2->query.offset))
    match a_q_offset.cmp(&b_q_offset) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // q_end DESC (BLAST_CMP(hsp2->query.end, hsp1->query.end))
    b_q_end.cmp(&a_q_end)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1415-1434
// ```c
// static int
// s_EvalueCompareHSPs(const void* v1, const void* v2)
// {
//    if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
//       return retval;
//    return ScoreCompareHSPs(v1, v2);
// }
// ```
fn evalue_compare_hsps(a: &Hit, b: &Hit) -> Ordering {
    match evalue_comp(a.e_value, b.e_value) {
        Ordering::Equal => score_compare_hsps(a, b),
        ord => ord,
    }
}

/// Subject group for NCBI-style output ordering.
/// Represents all HSPs for a single subject (BlastHSPList equivalent).
#[derive(Debug)]
struct SubjectGroup {
    s_idx: u32,
    best_evalue: f64,
    best_score: i32,
    hits: Vec<Hit>,
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1431-1454
// ```c
// if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
//    return retval;
// return ScoreCompareHSPs(v1, v2);
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       s_EvalueCompareHSPs);
// ```
#[derive(Clone, Copy)]
enum HspOutputOrder {
    ScoreCompare,
    EvalueCompare,
}

/// NCBI s_EvalueCompareHSPLists - compare two subjects by best e-value/score/oid.
/// Used to sort subjects within a query.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c s_EvalueCompareHSPLists()
/// ```c
/// static int s_EvalueCompareHSPLists(const void* v1, const void* v2) {
///    BlastHSPList* h1,* h2;
///    int retval = 0;
///    h1 = *(BlastHSPList**) v1;
///    h2 = *(BlastHSPList**) v2;
///    if (h1->hspcnt == 0 && h2->hspcnt == 0) return 0;
///    else if (h1->hspcnt == 0) return 1;
///    else if (h2->hspcnt == 0) return -1;
///    if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0)
///       return retval;
///    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
///    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
///    return BLAST_CMP(h2->oid, h1->oid);
/// }
/// ```
/// Order: best_evalue ASC → best_score DESC → oid DESC
fn compare_subject_groups(a: &SubjectGroup, b: &SubjectGroup) -> Ordering {
    // Empty groups go to the end
    if a.hits.is_empty() && b.hits.is_empty() {
        return Ordering::Equal;
    } else if a.hits.is_empty() {
        return Ordering::Greater;
    } else if b.hits.is_empty() {
        return Ordering::Less;
    }

    // best_evalue ASC (s_EvalueComp)
    match evalue_comp(a.best_evalue, b.best_evalue) {
        Ordering::Equal => {}
        ord => return ord,
    }

    // best_score DESC
    match b.best_score.cmp(&a.best_score) {
        Ordering::Equal => {}
        ord => return ord,
    }

    // oid DESC (BLAST_CMP(h2->oid, h1->oid))
    b.s_idx.cmp(&a.s_idx)
}

type HitCompare = fn(&Hit, &Hit) -> Ordering;
type SubjectGroupCompare = fn(&SubjectGroup, &SubjectGroup) -> Ordering;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       ScoreCompareHSPs);
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1452-1454
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       s_EvalueCompareHSPs);
// ```
// LOSAT owns Hit values in a Vec rather than a BlastHSP* array. Sort pointer-
// sized indices with the native C qsort, then replay that order onto the Vec so
// comparator-equal entries follow the same qsort behavior without LOSAT-only
// tie-break fields.
fn qsort_hits_by(hits: &mut Vec<Hit>, compare: HitCompare) {
    if hits.len() <= 1 {
        return;
    }

    let mut slots: Vec<Option<Hit>> = std::mem::take(hits).into_iter().map(Some).collect();
    let mut indices: Vec<usize> = (0..slots.len()).collect();
    qsort_hit_indices(&slots, &mut indices, compare);

    hits.extend(
        indices
            .into_iter()
            .map(|index| slots[index].take().expect("qsort replay index used once")),
    );
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3331-3337
// ```c
// if (hit_list && hit_list->hsplist_count > 1) {
//    qsort(hit_list->hsplist_array, hit_list->hsplist_count,
//             sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
// }
// ```
// Sort BlastHSPList-equivalent subject groups through the same pointer-sized
// index replay used for HSPs.
fn qsort_subject_groups_by(subject_groups: &mut Vec<SubjectGroup>, compare: SubjectGroupCompare) {
    if subject_groups.len() <= 1 {
        return;
    }

    let mut slots: Vec<Option<SubjectGroup>> = std::mem::take(subject_groups)
        .into_iter()
        .map(Some)
        .collect();
    let mut indices: Vec<usize> = (0..slots.len()).collect();
    qsort_subject_group_indices(&slots, &mut indices, compare);

    subject_groups.extend(
        indices
            .into_iter()
            .map(|index| slots[index].take().expect("qsort replay index used once")),
    );
}

#[cfg(target_arch = "wasm32")]
fn qsort_hit_indices(slots: &[Option<Hit>], indices: &mut [usize], compare: HitCompare) {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
    // ```c
    // qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //       ScoreCompareHSPs);
    // ```
    // Wasm has no C qsort in this build. Keep the same index-array replay
    // shape and avoid adding ordering criteria beyond NCBI's comparator.
    indices.sort_by(|&lhs, &rhs| {
        compare(
            slots[lhs].as_ref().expect("qsort lhs slot present"),
            slots[rhs].as_ref().expect("qsort rhs slot present"),
        )
    });
}

#[cfg(target_arch = "wasm32")]
fn qsort_subject_group_indices(
    slots: &[Option<SubjectGroup>],
    indices: &mut [usize],
    compare: SubjectGroupCompare,
) {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3331-3337
    // ```c
    // qsort(hit_list->hsplist_array, hit_list->hsplist_count,
    //          sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    // ```
    // Wasm replays sorted indices so subject grouping follows the same
    // comparator surface as the native qsort wrapper.
    indices.sort_by(|&lhs, &rhs| {
        compare(
            slots[lhs].as_ref().expect("qsort lhs slot present"),
            slots[rhs].as_ref().expect("qsort rhs slot present"),
        )
    });
}

#[cfg(not(target_arch = "wasm32"))]
fn qsort_hit_indices(slots: &[Option<Hit>], indices: &mut [usize], compare: HitCompare) {
    output_qsort_native::qsort_hit_indices(slots, indices, compare);
}

#[cfg(not(target_arch = "wasm32"))]
fn qsort_subject_group_indices(
    slots: &[Option<SubjectGroup>],
    indices: &mut [usize],
    compare: SubjectGroupCompare,
) {
    output_qsort_native::qsort_subject_group_indices(slots, indices, compare);
}

#[cfg(not(target_arch = "wasm32"))]
mod output_qsort_native {
    use super::{Hit, HitCompare, SubjectGroup, SubjectGroupCompare};
    use std::cell::Cell;
    use std::cmp::Ordering;
    use std::ffi::c_void;
    use std::mem;
    use std::os::raw::c_int;

    #[derive(Clone, Copy)]
    struct HitQsortContext {
        slots: *const Option<Hit>,
        len: usize,
        compare: HitCompare,
    }

    #[derive(Clone, Copy)]
    struct SubjectGroupQsortContext {
        slots: *const Option<SubjectGroup>,
        len: usize,
        compare: SubjectGroupCompare,
    }

    thread_local! {
        static HIT_QSORT_CONTEXT: Cell<Option<HitQsortContext>> = const { Cell::new(None) };
        static SUBJECT_GROUP_QSORT_CONTEXT: Cell<Option<SubjectGroupQsortContext>> = const { Cell::new(None) };
    }

    unsafe extern "C" {
        fn qsort(
            base: *mut c_void,
            nmemb: usize,
            size: usize,
            compar: extern "C" fn(*const c_void, *const c_void) -> c_int,
        );
    }

    extern "C" fn compare_hit_indices(lhs: *const c_void, rhs: *const c_void) -> c_int {
        // Safety: qsort passes pointers to elements of the usize index array
        // supplied in qsort_hit_indices.
        let lhs_index = unsafe { *(lhs as *const usize) };
        let rhs_index = unsafe { *(rhs as *const usize) };

        HIT_QSORT_CONTEXT.with(|cell| {
            let Some(ctx) = cell.get() else {
                return 0;
            };
            if lhs_index >= ctx.len || rhs_index >= ctx.len {
                return 0;
            }

            // Safety: ctx.slots points at the immutable slots slice that
            // outlives the qsort call; indices are checked against ctx.len.
            let lhs_hit = unsafe { &*ctx.slots.add(lhs_index) }
                .as_ref()
                .expect("qsort lhs slot present");
            let rhs_hit = unsafe { &*ctx.slots.add(rhs_index) }
                .as_ref()
                .expect("qsort rhs slot present");

            match (ctx.compare)(lhs_hit, rhs_hit) {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            }
        })
    }

    extern "C" fn compare_subject_group_indices(lhs: *const c_void, rhs: *const c_void) -> c_int {
        // Safety: qsort passes pointers to elements of the usize index array
        // supplied in qsort_subject_group_indices.
        let lhs_index = unsafe { *(lhs as *const usize) };
        let rhs_index = unsafe { *(rhs as *const usize) };

        SUBJECT_GROUP_QSORT_CONTEXT.with(|cell| {
            let Some(ctx) = cell.get() else {
                return 0;
            };
            if lhs_index >= ctx.len || rhs_index >= ctx.len {
                return 0;
            }

            // Safety: ctx.slots points at the immutable slots slice that
            // outlives the qsort call; indices are checked against ctx.len.
            let lhs_group = unsafe { &*ctx.slots.add(lhs_index) }
                .as_ref()
                .expect("qsort lhs slot present");
            let rhs_group = unsafe { &*ctx.slots.add(rhs_index) }
                .as_ref()
                .expect("qsort rhs slot present");

            match (ctx.compare)(lhs_group, rhs_group) {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            }
        })
    }

    pub(super) fn qsort_hit_indices(
        slots: &[Option<Hit>],
        indices: &mut [usize],
        compare: HitCompare,
    ) {
        let ctx = HitQsortContext {
            slots: slots.as_ptr(),
            len: slots.len(),
            compare,
        };

        HIT_QSORT_CONTEXT.with(|cell| {
            let previous = cell.replace(Some(ctx));
            // Safety: indices is a valid mutable array of pointer-sized usize
            // elements; the comparator only reads from `slots`, which outlives
            // this call.
            unsafe {
                qsort(
                    indices.as_mut_ptr().cast::<c_void>(),
                    indices.len(),
                    mem::size_of::<usize>(),
                    compare_hit_indices,
                );
            }
            cell.set(previous);
        });
    }

    pub(super) fn qsort_subject_group_indices(
        slots: &[Option<SubjectGroup>],
        indices: &mut [usize],
        compare: SubjectGroupCompare,
    ) {
        let ctx = SubjectGroupQsortContext {
            slots: slots.as_ptr(),
            len: slots.len(),
            compare,
        };

        SUBJECT_GROUP_QSORT_CONTEXT.with(|cell| {
            let previous = cell.replace(Some(ctx));
            // Safety: indices is a valid mutable array of pointer-sized usize
            // elements; the comparator only reads from `slots`, which outlives
            // this call.
            unsafe {
                qsort(
                    indices.as_mut_ptr().cast::<c_void>(),
                    indices.len(),
                    mem::size_of::<usize>(),
                    compare_subject_group_indices,
                );
            }
            cell.set(previous);
        });
    }
}

pub fn write_output(
    hits: &[Hit],
    out_path: Option<&PathBuf>,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
) -> Result<()> {
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };

    for hit in hits {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        // } BlastHSPList;
        // ```
        let (query_id, subject_id) = hit.resolve_ids(query_ids, subject_ids);
        writeln!(
            writer,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1e}\t{:.1}",
            query_id,
            subject_id,
            hit.identity,
            hit.length,
            hit.mismatch,
            hit.gapopen,
            hit.q_start,
            hit.q_end,
            hit.s_start,
            hit.s_end,
            hit.e_value,
            hit.bit_score
        )?;
    }
    Ok(())
}

// =============================================================================
// NCBI-compatible output with format selection
// =============================================================================

/// Write output with NCBI-compatible formatting and format selection
///
/// Supports:
/// - outfmt 0: Pairwise alignment view
/// - outfmt 6: Tabular output (default)
/// - outfmt 7: Tabular with comment lines
pub fn write_output_with_format(
    hits: &[Hit],
    out_path: Option<&PathBuf>,
    outfmt: OutputFormat,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
) -> Result<()> {
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };

    // NCBI-compatible output config
    let config = OutputConfig::ncbi_compat();

    match outfmt {
        OutputFormat::Pairwise => {
            // outfmt 0: Pairwise alignment view
            let pairwise_config = PairwiseConfig {
                program: context.program.clone(),
                ..Default::default()
            };
            write_pairwise_simple(
                hits,
                &mut writer,
                &pairwise_config,
                query_ids,
                subject_ids,
                context,
            )?;
        }
        OutputFormat::Tabular => {
            // outfmt 6: Tabular output
            // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
            // ```c
            // void CBlastTabularInfo::Print()
            // {
            //     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
            //         if (iter != m_FieldsToShow.begin())
            //             m_Ostream << m_FieldDelimiter;
            //         x_PrintField(*iter);
            //     }
            //     m_Ostream << "\n";
            // }
            // ```
            for hit in hits {
                // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
                // ```c
                // typedef struct BlastHSPList {
                //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
                //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
                //                       Set to 0 if not applicable */
                // } BlastHSPList;
                // ```
                let (query_id, subject_id) = hit.resolve_ids(query_ids, subject_ids);
                write_hit_fields(
                    &mut writer,
                    query_id,
                    subject_id,
                    hit.identity,
                    hit.num_ident,
                    hit.length,
                    hit.mismatch,
                    hit.gapopen,
                    hit.q_start,
                    hit.q_end,
                    hit.s_start,
                    hit.s_end,
                    hit.e_value,
                    hit.bit_score,
                    &config,
                )?;
            }
        }
        OutputFormat::TabularWithComments => {
            // outfmt 7: Tabular with comment lines
            write_outfmt7(hits, &mut writer, &config, query_ids, subject_ids, context)?;
        }
    }

    Ok(())
}

/// Write output with NCBI-style ordering:
/// 1. Group by query (input order)
/// 2. Within query, sort subjects by s_EvalueCompareHSPLists (best_evalue ASC → best_score DESC → oid DESC)
/// 3. Within subject, sort HSPs by ScoreCompareHSPs (score DESC → s_start ASC → s_end DESC → q_start ASC → q_end DESC)
///
/// Reference:
/// - BLAST_LinkHsps() calls Blast_HSPListSortByScore() after linking
/// - BlastHitList sorts subjects by s_EvalueCompareHSPLists
/// - Final output iterates query → subject → HSP
pub fn write_output_ncbi_order(
    hits: Vec<Hit>,
    out_path: Option<&PathBuf>,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
) -> Result<()> {
    write_output_ncbi_order_with_format(
        hits,
        out_path,
        OutputFormat::Tabular,
        query_ids,
        subject_ids,
        &ReportContext::default(),
    )
}

/// Write default tabular output with NCBI-style ordering to an existing writer.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
/// ```c
/// void CBlastTabularInfo::Print()
/// {
///     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
///         if (iter != m_FieldsToShow.begin())
///             m_Ostream << m_FieldDelimiter;
///         x_PrintField(*iter);
///     }
///     m_Ostream << "\n";
/// }
/// ```
pub fn write_output_ncbi_order_to_writer<W: Write>(
    hits: Vec<Hit>,
    writer: &mut W,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
) -> Result<()> {
    write_output_ncbi_order_with_format_to_writer(
        hits,
        writer,
        OutputFormat::Tabular,
        query_ids,
        subject_ids,
        &ReportContext::default(),
    )
}

/// Write default tabular output with NCBI subject ordering and formatter HSP
/// ordering by e-value.
///
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1569-1577
/// ```c
/// for (int index = 0; index < hit_list->hsplist_count; index++) {
///     BlastHSPList* hsp_list = hit_list->hsplist_array[index];
///     Blast_HSPListSortByEvalue(hsp_list);
/// }
/// ```
pub fn write_output_ncbi_order_evalue_hsp_order(
    hits: Vec<Hit>,
    out_path: Option<&PathBuf>,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
) -> Result<()> {
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };

    write_output_ncbi_order_with_format_to_writer_impl(
        hits,
        &mut writer,
        OutputFormat::Tabular,
        query_ids,
        subject_ids,
        &ReportContext::default(),
        HspOutputOrder::EvalueCompare,
    )
}

/// Write default tabular output with NCBI subject ordering to an existing writer,
/// using formatter HSP ordering by e-value.
///
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1437-1455
/// ```c
/// void Blast_HSPListSortByEvalue(BlastHSPList* hsp_list)
/// {
///    qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
///          s_EvalueCompareHSPs);
/// }
/// ```
pub fn write_output_ncbi_order_evalue_hsp_order_to_writer<W: Write>(
    hits: Vec<Hit>,
    writer: &mut W,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
) -> Result<()> {
    write_output_ncbi_order_with_format_to_writer_impl(
        hits,
        writer,
        OutputFormat::Tabular,
        query_ids,
        subject_ids,
        &ReportContext::default(),
        HspOutputOrder::EvalueCompare,
    )
}

/// Write output with NCBI-style ordering and format selection
///
/// Supports:
/// - outfmt 0: Pairwise alignment view  
/// - outfmt 6: Tabular output (default)
/// - outfmt 7: Tabular with comment lines
pub fn write_output_ncbi_order_with_format(
    hits: Vec<Hit>,
    out_path: Option<&PathBuf>,
    outfmt: OutputFormat,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
) -> Result<()> {
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };

    write_output_ncbi_order_with_format_to_writer(
        hits,
        &mut writer,
        outfmt,
        query_ids,
        subject_ids,
        context,
    )
}

/// Write output with NCBI-style ordering to an existing writer.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3383-3400
/// ```c
/// Int2 Blast_HSPResultsSortByEvalue(BlastHSPResults* results)
/// {
///    for (index = 0; index < results->num_queries; ++index) {
///       hit_list = results->hitlist_array[index];
///       if (hit_list != NULL && hit_list->hsplist_count > 1) {
///          qsort(hit_list->hsplist_array, hit_list->hsplist_count,
///                sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
///       }
///       s_BlastHitListPurge(hit_list);
///    }
/// }
/// ```
pub fn write_output_ncbi_order_with_format_to_writer<W: Write>(
    hits: Vec<Hit>,
    writer: &mut W,
    outfmt: OutputFormat,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
) -> Result<()> {
    write_output_ncbi_order_with_format_to_writer_impl(
        hits,
        writer,
        outfmt,
        query_ids,
        subject_ids,
        context,
        HspOutputOrder::ScoreCompare,
    )
}

fn write_output_ncbi_order_with_format_to_writer_impl<W: Write>(
    mut hits: Vec<Hit>,
    writer: &mut W,
    outfmt: OutputFormat,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
    hsp_order: HspOutputOrder,
) -> Result<()> {
    if hits.is_empty() {
        // For outfmt 7, still write header even with no hits
        if outfmt == OutputFormat::TabularWithComments {
            crate::report::write_outfmt7_header(writer, context, 0)?;
        } else if outfmt == OutputFormat::Pairwise {
            writeln!(writer, " ***** No hits found *****")?;
        }
        return Ok(());
    }

    // NCBI-compatible output config
    let config = OutputConfig::ncbi_compat();

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3391-3397
    // ```c
    // for (index = 0; index < results->num_queries; ++index) {
    //    hit_list = results->hitlist_array[index];
    //    if (hit_list != NULL
    //            && hit_list->hsplist_count > 1
    //            && hit_list->hsplist_array != NULL) {
    //       qsort(hit_list->hsplist_array, hit_list->hsplist_count,
    // ```
    // NCBI walks query slots by query index. Do the same instead of depending
    // on first-seen order from a flattened LOSAT Vec.
    let mut query_order: Vec<u32> = hits.iter().map(|h| h.q_idx).collect();
    query_order.sort_unstable();
    query_order.dedup();

    // Group hits by (q_idx, s_idx)
    let mut query_subject_hits: HashMap<(u32, u32), Vec<Hit>> = HashMap::new();
    for h in hits.drain(..) {
        query_subject_hits
            .entry((h.q_idx, h.s_idx))
            .or_default()
            .push(h);
    }

    // For outfmt 0 (pairwise), we need to collect all sorted hits first
    let mut all_sorted_hits: Vec<Hit> = Vec::new();

    // Step 2: For each query, build subject groups and sort
    for &q_idx in &query_order {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3331-3337
        // ```c
        // if (hit_list && hit_list->hsplist_count > 1) {
        //    qsort(hit_list->hsplist_array, hit_list->hsplist_count,
        //             sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
        // }
        // s_BlastHitListPurge(hit_list);
        // ```
        // Start from deterministic OID order, then apply NCBI's HSPList
        // comparator below. This avoids HashSet iteration order for groups
        // that will later be handed to a partial qsort-equivalent comparator.
        let mut subject_indices: Vec<u32> = query_subject_hits
            .keys()
            .filter(|(qidx, _)| *qidx == q_idx)
            .map(|(_, sidx)| *sidx)
            .collect();
        subject_indices.sort_unstable();
        subject_indices.dedup();

        // Build subject groups
        let mut subject_groups: Vec<SubjectGroup> = subject_indices
            .into_iter()
            .filter_map(|s_idx| {
                let key = (q_idx, s_idx);
                query_subject_hits.remove(&key).map(|mut hsp_hits| {
                    match hsp_order {
                        HspOutputOrder::ScoreCompare => {
                            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1381
                            // ```c
                            // if (hsp_list != NULL && hsp_list->hspcnt > 1
                            //         && hsp_list->hsp_array != NULL) {
                            //    qsort(hsp_list->hsp_array, hsp_list->hspcnt,
                            //          sizeof(BlastHSP*), ScoreCompareHSPs);
                            // }
                            // ```
                            qsort_hits_by(&mut hsp_hits, score_compare_hsps);
                        }
                        HspOutputOrder::EvalueCompare => {
                            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
                            // ```c
                            // // Sort HSPs with e-values as first priority and scores as
                            // // tie-breakers, since that is the order we want to see them in
                            // // in Seq-aligns.
                            // Blast_HSPListSortByEvalue(hsp_list);
                            // ```
                            qsort_hits_by(&mut hsp_hits, evalue_compare_hsps);
                        }
                    }

                    // Compute best_evalue and best_score
                    let best_evalue = hsp_hits
                        .iter()
                        .map(|h| h.e_value)
                        .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                        .unwrap_or(f64::MAX);
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3098-3104
                    // ```c
                    // if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0)
                    //    return retval;
                    // if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
                    // if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
                    // ```
                    let best_score = hsp_hits.first().map(|h| h.raw_score).unwrap_or(0);

                    SubjectGroup {
                        s_idx,
                        best_evalue,
                        best_score,
                        hits: hsp_hits,
                    }
                })
            })
            .collect();

        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3331-3337
        // ```c
        // if (hit_list && hit_list->hsplist_count > 1) {
        //    qsort(hit_list->hsplist_array, hit_list->hsplist_count,
        //             sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
        // }
        // ```
        qsort_subject_groups_by(&mut subject_groups, compare_subject_groups);

        // Step 3: Output in order based on format
        match outfmt {
            OutputFormat::Pairwise => {
                // Collect hits for later pairwise output
                for group in subject_groups {
                    all_sorted_hits.extend(group.hits);
                }
            }
            OutputFormat::Tabular => {
                // outfmt 6: Tabular output with NCBI formatting
                for group in subject_groups {
                    for hit in group.hits {
                        // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
                        // ```c
                        // void CBlastTabularInfo::Print()
                        // {
                        //     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
                        //         if (iter != m_FieldsToShow.begin())
                        //             m_Ostream << m_FieldDelimiter;
                        //         x_PrintField(*iter);
                        //     }
                        //     m_Ostream << "\n";
                        // }
                        // ```
                        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
                        // ```c
                        // typedef struct BlastHSPList {
                        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
                        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
                        //                       Set to 0 if not applicable */
                        // } BlastHSPList;
                        // ```
                        let (query_id, subject_id) = hit.resolve_ids(query_ids, subject_ids);
                        write_hit_fields(
                            writer,
                            query_id,
                            subject_id,
                            hit.identity,
                            hit.num_ident,
                            hit.length,
                            hit.mismatch,
                            hit.gapopen,
                            hit.q_start,
                            hit.q_end,
                            hit.s_start,
                            hit.s_end,
                            hit.e_value,
                            hit.bit_score,
                            &config,
                        )?;
                    }
                }
            }
            OutputFormat::TabularWithComments => {
                // outfmt 7: Collect for later output with headers
                for group in subject_groups {
                    all_sorted_hits.extend(group.hits);
                }
            }
        }
    }

    // Handle formats that need post-processing
    match outfmt {
        OutputFormat::Pairwise => {
            let pairwise_config = PairwiseConfig {
                program: context.program.clone(),
                ..Default::default()
            };
            write_pairwise_simple(
                &all_sorted_hits,
                writer,
                &pairwise_config,
                query_ids,
                subject_ids,
                context,
            )?;
        }
        OutputFormat::TabularWithComments => {
            write_outfmt7(
                &all_sorted_hits,
                writer,
                &config,
                query_ids,
                subject_ids,
                context,
            )?;
        }
        OutputFormat::Tabular => {
            // Already written above
        }
    }

    Ok(())
}
