use std::cmp::Ordering;

// =============================================================================
// Gap Edit Script for traceback trimming
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c
// =============================================================================

/// Gap edit operation type for traceback
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c eGapAlignOpType
///
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapEditOp {
    /// Substitution (match or mismatch): consumes `num` query and `num` subject positions
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2404-2407
    Sub(u32),
    /// Deletion: gap in query, consumes `num` subject positions only
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2408-2410
    Del(u32),
    /// Insertion: gap in subject, consumes `num` query positions only
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2411-2413
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
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
    pub identity: f64,
    pub length: usize,
    // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:476-480
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
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
    pub num_ident: usize,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1132
    pub query_frame: i32,
    pub query_length: usize,
    // Fields for NCBI-style output ordering (not printed, used for sorting)
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
    /// Query index (input order) - NCBI uses query order for grouping
    pub q_idx: u32,
    /// Subject index (oid) - NCBI uses oid for tie-breaking in s_EvalueCompareHSPLists
    pub s_idx: u32,
    /// Raw alignment score - NCBI ScoreCompareHSPs uses raw score, not bit score
    pub raw_score: i32,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2388
    // Optional context/frame-relative BlastSeg offsets used by ScoreCompareHSPs
    // when formatted nucleotide coordinates cannot recover the translated HSP
    // coordinates.
    pub sort_query_offset: usize,
    pub sort_query_end: usize,
    pub sort_subject_offset: usize,
    pub sort_subject_end: usize,
    pub has_sort_offsets: bool,
    /// Edit script for traceback trimming (NCBI: gap_info in BlastHSP)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c BlastHSP.gap_info (GapEditScript*)
    ///
    /// This is used by Blast_HSPListPurgeHSPsWithCommonEndpoints to trim
    /// overlapping HSPs instead of deleting them entirely.
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2490-2492, 2516-2518
    pub gap_info: Option<Vec<GapEditOp>>,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
    pub num_positives: usize,
}

impl Hit {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-824
    #[inline]
    pub fn gap_letters(&self) -> usize {
        usize::try_from(self.length as i64 - self.num_ident as i64 - self.mismatch).unwrap_or(0)
    }
}

// =============================================================================
// NCBI-style comparators for outfmt6 output ordering
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c
// =============================================================================

/// Compare two evalues, treating both as equal if they're close enough to zero.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c s_EvalueComp()
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1335-1355 ScoreCompareHSPs()
/// Order: score DESC → s_start ASC → s_end DESC → q_start ASC → q_end DESC
pub fn score_compare_hsps(a: &Hit, b: &Hit) -> Ordering {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1132
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1353
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1414-1434
pub(crate) fn evalue_compare_hsps(a: &Hit, b: &Hit) -> Ordering {
    match evalue_comp(a.e_value, b.e_value) {
        Ordering::Equal => score_compare_hsps(a, b),
        ord => ord,
    }
}
