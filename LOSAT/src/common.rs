use std::cmp::Ordering;
use std::sync::Arc;

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
    // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:476-480
    // ```c++
    // inline void CBlastTabularInfo::x_PrintMismatches(void)
    // {
    //     int num_mismatches = m_AlignLength - m_NumIdent - m_NumGaps;
    //     m_Ostream << num_mismatches;
    // }
    // ```
    // Fewer than none when a stored gap character of the subject (`-`) faces a gap of the
    // query: the letters are equal, an identity, and the column is a gap.
    pub mismatch: i64,
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
        usize::try_from(self.length as i64 - self.num_ident as i64 - self.mismatch).unwrap_or(0)
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
pub(crate) fn evalue_compare_hsps(a: &Hit, b: &Hit) -> Ordering {
    match evalue_comp(a.e_value, b.e_value) {
        Ordering::Equal => score_compare_hsps(a, b),
        ord => ord,
    }
}
