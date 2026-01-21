use std::cmp::Ordering;
use std::sync::Arc;

use crate::common::{GapEditOp, Hit};

#[derive(Debug, Clone)]
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
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
// } BlastHSP;
// ```
pub struct BlastnHsp {
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
    pub query_frame: i32,
    pub query_length: usize,
    pub q_idx: u32,
    pub s_idx: u32,
    pub raw_score: i32,
    pub gap_info: Option<Vec<GapEditOp>>,
}

#[derive(Debug)]
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
// ```c
// typedef struct BlastHSPList {
//    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
//    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
//    BlastHSP** hsp_array;
//    Int4 hspcnt;
//    double best_evalue;
// } BlastHSPList;
// ```
pub struct BlastnHspList {
    pub oid: u32,
    pub query_index: u32,
    pub hsps: Vec<BlastnHsp>,
    pub best_evalue: f64,
}

#[derive(Debug)]
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:168-180
// ```c
// typedef struct BlastHitList {
//    Int4 hsplist_count;
//    BlastHSPList** hsplist_array;
// } BlastHitList;
// ```
pub struct BlastnHitList {
    pub hsplist_array: Vec<BlastnHspList>,
}

impl BlastnHsp {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    Int4 num_ident;
    //    double bit_score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    //    GapEditScript* gap_info;
    // } BlastHSP;
    // ```
    pub fn from_hit(hit: Hit) -> Self {
        let Hit {
            query_id: _,
            subject_id: _,
            identity,
            length,
            mismatch,
            gapopen,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value,
            bit_score,
            query_frame,
            query_length,
            q_idx,
            s_idx,
            raw_score,
            gap_info,
        } = hit;
        Self {
            identity,
            length,
            mismatch,
            gapopen,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value,
            bit_score,
            query_frame,
            query_length,
            q_idx,
            s_idx,
            raw_score,
            gap_info,
        }
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
    // } BlastHSPList;
    // ```
    pub fn into_hit(self, query_ids: &[Arc<str>], subject_ids: &[Arc<str>]) -> Hit {
        let query_id = query_ids
            .get(self.q_idx as usize)
            .cloned()
            .unwrap_or_else(|| Arc::<str>::from("unknown"));
        let subject_id = subject_ids
            .get(self.s_idx as usize)
            .cloned()
            .unwrap_or_else(|| Arc::<str>::from("unknown"));
        Hit {
            query_id,
            subject_id,
            identity: self.identity,
            length: self.length,
            mismatch: self.mismatch,
            gapopen: self.gapopen,
            q_start: self.q_start,
            q_end: self.q_end,
            s_start: self.s_start,
            s_end: self.s_end,
            e_value: self.e_value,
            bit_score: self.bit_score,
            query_frame: self.query_frame,
            query_length: self.query_length,
            q_idx: self.q_idx,
            s_idx: self.s_idx,
            raw_score: self.raw_score,
            gap_info: self.gap_info,
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1389-1403
// ```c
// static int s_EvalueComp(double evalue1, double evalue2)
// {
//     const double epsilon = 1.0e-180;
//     if (evalue1 < epsilon && evalue2 < epsilon) { return 0; }
//     if (evalue1 < evalue2) return -1;
//     else if (evalue1 > evalue2) return 1;
//     else return 0;
// }
// ```
pub fn evalue_comp(evalue1: f64, evalue2: f64) -> Ordering {
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// int ScoreCompareHSPs(const void* h1, const void* h2) {
//    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
//        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
//        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
//        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
//        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
//    }
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
// ```c
// if (hsp->query.frame != hsp->subject.frame) {
//    *q_end = query_length - hsp->query.offset;
//    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
// }
// ```
pub fn score_compare_hsps(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    let query_offsets = |h: &BlastnHsp| {
        if h.query_length > 0 && h.query_frame < 0 {
            let q_offset = h.query_length.saturating_sub(h.q_end);
            let q_end = h
                .query_length
                .saturating_sub(h.q_start)
                .saturating_add(1);
            (q_offset, q_end)
        } else {
            (h.q_start.saturating_sub(1), h.q_end)
        }
    };
    let (a_q_offset, a_q_end) = query_offsets(a);
    let (b_q_offset, b_q_end) = query_offsets(b);
    let (a_s_offset, a_s_end) = (a.s_start.min(a.s_end).saturating_sub(1), a.s_start.max(a.s_end));
    let (b_s_offset, b_s_end) = (b.s_start.min(b.s_end).saturating_sub(1), b.s_start.max(b.s_end));

    match b.raw_score.cmp(&a.raw_score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match a_s_offset.cmp(&b_s_offset) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match b_s_end.cmp(&a_s_end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match a_q_offset.cmp(&b_q_offset) {
        Ordering::Equal => {}
        ord => return ord,
    }
    b_q_end.cmp(&a_q_end)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1739-1748
// ```c
// static double s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
// {
//     double best_evalue = (double) INT4_MAX;
//     for (index=0; index<hsp_list->hspcnt; index++)
//         best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
//     return best_evalue;
// }
// ```
pub fn update_best_evalue(list: &mut BlastnHspList) {
    if list.hsps.is_empty() {
        list.best_evalue = f64::MAX;
        return;
    }
    let mut best = f64::MAX;
    for hsp in &list.hsps {
        if hsp.e_value < best {
            best = hsp.e_value;
        }
    }
    list.best_evalue = best;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2049-2067
// ```c
// Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
//                                const BlastHitSavingOptions* hit_options)
// {
//    if ((hsp_list == NULL) ||
//        (hit_options->max_hsps_per_subject == 0) ||
//        (hsp_list->hspcnt <= hit_options->max_hsps_per_subject))
//       return 0;
//    hsp_list->hspcnt = hsp_max;
// }
// ```
pub fn trim_by_max_hsps(list: &mut BlastnHspList, max_hsps_per_subject: usize) {
    if max_hsps_per_subject == 0 || list.hsps.len() <= max_hsps_per_subject {
        return;
    }
    list.hsps.truncate(max_hsps_per_subject);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3071-3106
// ```c
// static int s_EvalueCompareHSPLists(const void* v1, const void* v2)
// {
//    if (h1->hspcnt == 0 && h2->hspcnt == 0) return 0;
//    else if (h1->hspcnt == 0) return 1;
//    else if (h2->hspcnt == 0) return -1;
//    if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0)
//       return retval;
//    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
//    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
//    return BLAST_CMP(h2->oid, h1->oid);
// }
// ```
pub fn compare_hsp_lists(a: &BlastnHspList, b: &BlastnHspList) -> Ordering {
    if a.hsps.is_empty() && b.hsps.is_empty() {
        return Ordering::Equal;
    } else if a.hsps.is_empty() {
        return Ordering::Greater;
    } else if b.hsps.is_empty() {
        return Ordering::Less;
    }

    let cmp = evalue_comp(a.best_evalue, b.best_evalue);
    if cmp != Ordering::Equal {
        return cmp;
    }
    let a_score = a.hsps.first().map(|h| h.raw_score).unwrap_or(0);
    let b_score = b.hsps.first().map(|h| h.raw_score).unwrap_or(0);
    match b_score.cmp(&a_score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    b.oid.cmp(&a.oid)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
// ```c
// static void s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
// {
//    for (subject_index = hitlist_size;
//         subject_index < hit_list->hsplist_count; ++subject_index) {
//       hit_list->hsplist_array[subject_index] =
//       Blast_HSPListFree(hit_list->hsplist_array[subject_index]);
//    }
//    hit_list->hsplist_count = MIN(hit_list->hsplist_count, hitlist_size);
// }
// ```
pub fn prune_hitlists_by_size(lists: &mut Vec<BlastnHspList>, hitlist_size: usize) {
    if hitlist_size == 0 || lists.len() <= hitlist_size {
        return;
    }
    lists.truncate(hitlist_size);
}
