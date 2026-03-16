//! BLASTP-specific HSP list and hit list helpers.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c

use std::cmp::Ordering;
use std::collections::HashMap;

use crate::common::{GapEditOp, Hit};

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
#[derive(Debug, Clone)]
pub(crate) struct BlastpHsp {
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
    pub num_ident: usize,
    pub query_frame: i32,
    pub query_length: usize,
    pub q_idx: u32,
    pub s_idx: u32,
    pub raw_score: i32,
    pub gap_info: Option<Vec<GapEditOp>>,
    pub num_positives: usize,
}

impl BlastpHsp {
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
    pub(crate) fn from_hit(hit: Hit) -> Self {
        let Hit {
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
            num_ident,
            query_frame,
            query_length,
            q_idx,
            s_idx,
            raw_score,
            gap_info,
            num_positives,
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
            num_ident,
            query_frame,
            query_length,
            q_idx,
            s_idx,
            raw_score,
            gap_info,
            num_positives,
        }
    }

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
    pub(crate) fn into_hit(self) -> Hit {
        Hit {
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
            num_ident: self.num_ident,
            query_frame: self.query_frame,
            query_length: self.query_length,
            q_idx: self.q_idx,
            s_idx: self.s_idx,
            raw_score: self.raw_score,
            gap_info: self.gap_info,
            num_positives: self.num_positives,
        }
    }
}

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
#[derive(Debug)]
#[allow(dead_code)]
pub(crate) struct BlastpHspList {
    pub oid: u32,
    pub query_index: u32,
    pub hsps: Vec<BlastpHsp>,
    pub best_evalue: f64,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:168-180
// ```c
// typedef struct BlastHitList {
//    Int4 hsplist_count;
//    Int4 hsplist_max;
//    double worst_evalue;
//    Int4 low_score;
//    Boolean heapified;
//    BlastHSPList** hsplist_array;
//    Int4 hsplist_current;
//    Int4 num_hits;
// } BlastHitList;
// ```
#[derive(Debug)]
#[allow(dead_code)]
pub(crate) struct BlastpHitList {
    pub hsplist_count: usize,
    pub hsplist_max: usize,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplist_array: Vec<BlastpHspList>,
    pub hsplist_current: usize,
    pub num_hits: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70
// ```c
// Int4
// GetPrelimHitlistSize(Int4 hitlist_size, Int4 compositionBasedStats, Boolean gapped_calculation)
// {
//     ...
// }
// ```
pub(crate) fn get_prelim_hitlist_size(
    hitlist_size: usize,
    composition_based_stats: bool,
    gapped_calculation: bool,
) -> usize {
    let mut prelim_hitlist_size = hitlist_size;
    let adaptive_cbs = std::env::var_os("ADAPTIVE_CBS").is_some();
    if composition_based_stats {
        if adaptive_cbs {
            if hitlist_size < 1000 {
                prelim_hitlist_size = std::cmp::max(prelim_hitlist_size + 1000, 1500);
            } else {
                prelim_hitlist_size = prelim_hitlist_size.saturating_mul(2).saturating_add(50);
            }
        } else if hitlist_size <= 500 {
            prelim_hitlist_size = 1050;
        } else {
            prelim_hitlist_size = prelim_hitlist_size.saturating_mul(2).saturating_add(50);
        }
    } else if gapped_calculation {
        prelim_hitlist_size = std::cmp::min(
            std::cmp::max(prelim_hitlist_size.saturating_mul(2), 10),
            prelim_hitlist_size.saturating_add(50),
        );
    }
    prelim_hitlist_size
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// int ScoreCompareHSPs(const void* h1, const void* h2) {
//    BlastHSP* hsp1,* hsp2;
//    int result = 0;
//    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
//        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
//        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
//        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
//        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
//    }
//    return result;
// }
// ```
pub(crate) fn score_compare_hits(a: &Hit, b: &Hit) -> Ordering {
    let (a_q_offset, a_q_end) = (a.q_start.saturating_sub(1), a.q_end);
    let (b_q_offset, b_q_end) = (b.q_start.saturating_sub(1), b.q_end);
    let (a_s_offset, a_s_end) = (
        a.s_start.min(a.s_end).saturating_sub(1),
        a.s_start.max(a.s_end),
    );
    let (b_s_offset, b_s_end) = (
        b.s_start.min(b.s_end).saturating_sub(1),
        b.s_start.max(b.s_end),
    );

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// int ScoreCompareHSPs(const void* h1, const void* h2) {
//    ...
// }
// ```
pub(crate) fn score_compare_hsps(a: &BlastpHsp, b: &BlastpHsp) -> Ordering {
    let (a_q_offset, a_q_end) = (a.q_start.saturating_sub(1), a.q_end);
    let (b_q_offset, b_q_end) = (b.q_start.saturating_sub(1), b.q_end);
    let (a_s_offset, a_s_end) = (
        a.s_start.min(a.s_end).saturating_sub(1),
        a.s_start.max(a.s_end),
    );
    let (b_s_offset, b_s_end) = (
        b.s_start.min(b.s_end).saturating_sub(1),
        b.s_start.max(b.s_end),
    );

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1414-1435
// ```c
// static int
// s_EvalueCompareHSPs(const void* v1, const void* v2)
// {
//    ...
//    if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
//       return retval;
//    return ScoreCompareHSPs(v1, v2);
// }
// ```
pub(crate) fn evalue_compare_hsps(a: &BlastpHsp, b: &BlastpHsp) -> Ordering {
    let cmp = evalue_comp(a.e_value, b.e_value);
    if cmp == Ordering::Equal {
        score_compare_hsps(a, b)
    } else {
        cmp
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1437-1455
// ```c
// void Blast_HSPListSortByEvalue(BlastHSPList* hsp_list)
// {
//     ...
// }
// ```
fn sort_hsplist_by_evalue(list: &mut BlastpHspList) {
    if list.hsps.len() > 1 {
        let mut index = 0usize;
        while index < list.hsps.len() - 1 {
            if evalue_compare_hsps(&list.hsps[index], &list.hsps[index + 1]) == Ordering::Greater {
                break;
            }
            index += 1;
        }
        if index < list.hsps.len() - 1 {
            list.hsps.sort_by(evalue_compare_hsps);
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1739-1748
// ```c
// static double s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
// {
//     ...
// }
// ```
pub(crate) fn update_best_evalue(list: &mut BlastpHspList) {
    let mut best = i32::MAX as f64;
    for hsp in &list.hsps {
        if hsp.e_value < best {
            best = hsp.e_value;
        }
    }
    list.best_evalue = best;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2049-2067
// ```c
// Int2 Blast_HSPListReapByEvalue(BlastHSPList* hsp_list,
//                                const BlastHitSavingOptions* hit_options)
// {
//    ...
// }
// ```
pub(crate) fn reap_hsplist_by_evalue(list: &mut BlastpHspList, expect_value: f64) {
    let mut write_index = 0usize;
    for read_index in 0..list.hsps.len() {
        if list.hsps[read_index].e_value > expect_value {
            continue;
        }
        if read_index > write_index {
            list.hsps.swap(write_index, read_index);
        }
        write_index += 1;
    }
    list.hsps.truncate(write_index);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2049-2067
// ```c
// Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
//                                const BlastHitSavingOptions* hit_options)
// {
//    ...
// }
// ```
pub(crate) fn trim_by_max_hsps(list: &mut BlastpHspList, max_hsps_per_subject: usize) {
    if max_hsps_per_subject == 0 || list.hsps.len() <= max_hsps_per_subject {
        return;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2059-2067
    // ```c
    // hsp_max = hit_options->max_hsps_per_subject;
    // for (index = hsp_max; index < hsp_list->hspcnt; index++) {
    //    hsp_array[index] = Blast_HSPFree(hsp_array[index]);
    // }
    // hsp_list->hspcnt = hsp_max;
    // ```
    list.hsps.truncate(max_hsps_per_subject);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3077-3106
// ```c
// static int s_EvalueCompareHSPLists(const void* v1, const void* v2)
// {
//    ...
// }
// ```
fn compare_hsp_lists(a: &BlastpHspList, b: &BlastpHspList) -> Ordering {
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1627-1676
// ```c
// static void s_Heapify(...);
// static void s_CreateHeap(...);
// ```
fn heapify_hsplist_array(lists: &mut [BlastpHspList], start: usize, end: usize) {
    let mut root = start;
    loop {
        let left = root.saturating_mul(2).saturating_add(1);
        if left > end {
            break;
        }
        let mut large = left;
        let right = left + 1;
        if right <= end && compare_hsp_lists(&lists[left], &lists[right]) == Ordering::Less {
            large = right;
        }
        if compare_hsp_lists(&lists[root], &lists[large]) == Ordering::Less {
            lists.swap(root, large);
            root = large;
        } else {
            break;
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1659-1676
// ```c
// static void s_CreateHeap(...);
// ```
fn create_hsplist_heap(lists: &mut [BlastpHspList]) {
    let nel = lists.len();
    if nel < 2 {
        return;
    }
    let mut i = nel / 2;
    while i > 0 {
        i -= 1;
        heapify_hsplist_array(lists, i, nel - 1);
    }
}

impl BlastpHitList {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3125-3133
    // ```c
    // BlastHitList* Blast_HitListNew(Int4 hitlist_size)
    // {
    //    ...
    // }
    // ```
    pub(crate) fn new(hitlist_size: usize) -> Self {
        Self {
            hsplist_count: 0,
            hsplist_max: hitlist_size,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
            hsplist_array: Vec::new(),
            hsplist_current: 0,
            num_hits: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3219-3240
    // ```c
    // static Int2 s_Blast_HitListGrowHSPListArray(BlastHitList* hit_list)
    // {
    //     ...
    // }
    // ```
    fn grow_hsplist_array(&mut self) -> bool {
        const K_START_VALUE: usize = 100;
        if self.hsplist_current >= self.hsplist_max {
            return false;
        }
        if self.hsplist_current == 0 {
            self.hsplist_current = K_START_VALUE.min(self.hsplist_max);
        } else {
            let next = self.hsplist_current.saturating_mul(2);
            self.hsplist_current = next.min(self.hsplist_max);
        }
        if self.hsplist_array.capacity() < self.hsplist_current {
            self.hsplist_array
                .reserve(self.hsplist_current - self.hsplist_array.capacity());
        }
        true
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3196-3209
    // ```c
    // static void s_BlastHitListInsertHSPListInHeap(BlastHitList* hit_list,
    //                                              BlastHSPList* hsp_list)
    // {
    //    ...
    // }
    // ```
    fn insert_hsplist_in_heap(&mut self, hsp_list: BlastpHspList) {
        if self.hsplist_array.is_empty() {
            self.hsplist_array.push(hsp_list);
            self.hsplist_count = 1;
        } else {
            self.hsplist_array[0] = hsp_list;
        }
        if self.hsplist_count >= 2 {
            heapify_hsplist_array(&mut self.hsplist_array, 0, self.hsplist_count - 1);
        }
        if let Some(root) = self.hsplist_array.first() {
            self.worst_evalue = root.best_evalue;
            if let Some(score) = root.hsps.first().map(|h| h.raw_score) {
                self.low_score = score;
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3243-3297
    // ```c
    // Int2 Blast_HitListUpdate(BlastHitList* hit_list, BlastHSPList* hsp_list)
    // {
    //    ...
    // }
    // ```
    pub(crate) fn update(&mut self, mut hsp_list: BlastpHspList) {
        update_best_evalue(&mut hsp_list);

        if self.hsplist_count < self.hsplist_max {
            if self.hsplist_current == self.hsplist_count && !self.grow_hsplist_array() {
                return;
            }
            self.hsplist_array.push(hsp_list);
            self.hsplist_count += 1;
            self.worst_evalue = self
                .worst_evalue
                .max(self.hsplist_array.last().unwrap().best_evalue);
            if let Some(score) = self
                .hsplist_array
                .last()
                .unwrap()
                .hsps
                .first()
                .map(|h| h.raw_score)
            {
                self.low_score = self.low_score.min(score);
            }
        } else {
            if !self.heapified {
                for list in &mut self.hsplist_array {
                    sort_hsplist_by_evalue(list);
                    update_best_evalue(list);
                }
                create_hsplist_heap(&mut self.hsplist_array);
                self.heapified = true;
            }

            sort_hsplist_by_evalue(&mut hsp_list);
            update_best_evalue(&mut hsp_list);
            let evalue_order = compare_hsp_lists(&self.hsplist_array[0], &hsp_list);
            if evalue_order == Ordering::Less {
                return;
            }
            self.insert_hsplist_in_heap(hsp_list);
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3172-3187
    // ```c
    // static void s_BlastHitListPurge(BlastHitList* hit_list)
    // {
    //    ...
    // }
    // ```
    fn purge(&mut self) {
        let mut index = 0usize;
        while index < self.hsplist_count {
            if self.hsplist_array[index].hsps.is_empty() {
                break;
            }
            index += 1;
        }
        self.hsplist_count = index;
        self.hsplist_array.truncate(index);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3331-3337
    // ```c
    // Int2 Blast_HitListSortByEvalue(BlastHitList* hit_list)
    // {
    //    ...
    // }
    // ```
    pub(crate) fn sort_by_evalue(&mut self) {
        if self.hsplist_count > 1 {
            self.hsplist_array.sort_by(compare_hsp_lists);
        }
        self.purge();
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
    // ```c
    // static void s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
    // {
    //    ...
    // }
    // ```
    pub(crate) fn prune_by_size(&mut self, hitlist_size: usize) {
        if hitlist_size == 0 {
            self.hsplist_array.clear();
            self.hsplist_count = 0;
            return;
        }
        if self.hsplist_count > hitlist_size {
            self.hsplist_array.truncate(hitlist_size);
            self.hsplist_count = hitlist_size;
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
// ```c
// Int4
// Blast_HSPListPurgeHSPsWithCommonEndpoints(EBlastProgramType program,
//                                           BlastHSPList* hsp_list,
//                                           Boolean purge)
// {
//    ...
// }
// ```
pub(crate) fn purge_hits_with_common_endpoints(hits: Vec<Hit>) -> Vec<Hit> {
    let (hits, _) = purge_hits_with_common_endpoints_ex(hits, true);
    hits
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
// ```c
// Int4
// Blast_HSPListPurgeHSPsWithCommonEndpoints(..., Boolean purge)
// {
//    ...
// }
// ```
#[allow(dead_code)]
pub(crate) fn purge_hits_with_common_endpoints_ex(
    hits: Vec<Hit>,
    purge: bool,
) -> (Vec<Hit>, usize) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len);
    }

    let mut subject_groups: HashMap<u32, Vec<Hit>> = HashMap::new();
    for hit in hits {
        subject_groups.entry(hit.s_idx).or_default().push(hit);
    }

    let mut result = Vec::new();
    let mut total_extra_start = 0usize;
    for (_subject_idx, group_hits) in subject_groups {
        let (purged, extra_start) = purge_hits_for_subject_ex(group_hits, purge);
        let offset = result.len();
        result.extend(purged);
        if extra_start > 0 {
            total_extra_start = offset + extra_start;
        }
    }
    if total_extra_start == 0 {
        total_extra_start = result.len();
    }
    (result, total_extra_start)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2392-2452
// ```c
// static void
// s_CutOffGapEditScript(BlastHSP* hsp, Int4 q_cut, Int4 s_cut, Boolean cut_begin)
// {
//    ...
// }
// ```
fn cut_off_gap_edit_script(hit: &mut Hit, q_cut: usize, s_cut: usize, cut_begin: bool) -> bool {
    let gap_info = match &hit.gap_info {
        Some(info) if !info.is_empty() => info.clone(),
        _ => return false,
    };

    let (mut q_offset_0, mut q_end_0) = (hit.q_start.saturating_sub(1), hit.q_end);
    let (mut s_offset_0, mut s_end_0) = (
        hit.s_start.min(hit.s_end).saturating_sub(1),
        hit.s_start.max(hit.s_end),
    );

    let q_cut_rel = q_cut.saturating_sub(q_offset_0);
    let s_cut_rel = s_cut.saturating_sub(s_offset_0);

    let mut qid = 0usize;
    let mut sid = 0usize;
    let mut found = false;
    let mut found_index = 0usize;
    let mut found_opid = 0usize;

    'outer: for (index, op) in gap_info.iter().enumerate() {
        let num = op.num() as usize;
        match op {
            GapEditOp::Sub(_) => {
                for opid in 0..num {
                    qid += 1;
                    sid += 1;
                    if qid >= q_cut_rel && sid >= s_cut_rel {
                        found = true;
                        found_index = index;
                        found_opid = opid + 1;
                        break 'outer;
                    }
                }
            }
            GapEditOp::Del(n) => {
                sid += *n as usize;
                if qid >= q_cut_rel && sid >= s_cut_rel {
                    found = true;
                    found_index = index;
                    found_opid = *n as usize;
                    break 'outer;
                }
            }
            GapEditOp::Ins(n) => {
                qid += *n as usize;
                if qid >= q_cut_rel && sid >= s_cut_rel {
                    found = true;
                    found_index = index;
                    found_opid = *n as usize;
                    break 'outer;
                }
            }
        }
    }

    if !found {
        return false;
    }

    let new_gap_info = if cut_begin {
        let mut new_ops = Vec::new();
        if let GapEditOp::Sub(n) = gap_info[found_index] {
            if found_opid < n as usize {
                let remaining = n as usize - found_opid;
                if remaining > 0 {
                    new_ops.push(GapEditOp::Sub(remaining as u32));
                }
            }
        }
        for op in gap_info.iter().skip(found_index + 1) {
            new_ops.push(*op);
        }
        q_offset_0 = q_offset_0.saturating_add(qid);
        s_offset_0 = s_offset_0.saturating_add(sid);
        new_ops
    } else {
        let mut new_ops = Vec::new();
        for (i, op) in gap_info.iter().enumerate() {
            if i < found_index {
                new_ops.push(*op);
            } else if i == found_index {
                if let GapEditOp::Sub(n) = op {
                    if found_opid < *n as usize {
                        if found_opid > 0 {
                            new_ops.push(GapEditOp::Sub(found_opid as u32));
                        }
                    } else {
                        new_ops.push(*op);
                    }
                } else {
                    new_ops.push(*op);
                }
                break;
            }
        }
        q_end_0 = q_offset_0.saturating_add(qid);
        s_end_0 = s_offset_0.saturating_add(sid);
        new_ops
    };

    hit.q_start = q_offset_0 + 1;
    hit.q_end = q_end_0;
    hit.s_start = s_offset_0 + 1;
    hit.s_end = s_end_0;
    hit.gap_info = if new_gap_info.is_empty() {
        None
    } else {
        Some(new_gap_info)
    };
    if let Some(ref ops) = hit.gap_info {
        hit.length = ops.iter().map(|op| op.num() as usize).sum();
    }
    true
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
// ```c
// Int4 Blast_HSPListPurgeHSPsWithCommonEndpoints(..., Boolean purge)
// {
//    ...
// }
// ```
fn purge_hits_for_subject_ex(mut hits: Vec<Hit>, purge: bool) -> (Vec<Hit>, usize) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len);
    }

    let context = |h: &Hit| (h.q_idx, h.query_frame);
    let q_offset = |h: &Hit| h.q_start.saturating_sub(1);
    let q_end = |h: &Hit| h.q_end;
    let s_offset = |h: &Hit| h.s_start.min(h.s_end).saturating_sub(1);
    let s_end = |h: &Hit| h.s_start.max(h.s_end);

    let mut trimmed_hits = Vec::new();

    hits.sort_by(|a, b| {
        context(a)
            .cmp(&context(b))
            .then_with(|| q_offset(a).cmp(&q_offset(b)))
            .then_with(|| s_offset(a).cmp(&s_offset(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score))
            .then_with(|| q_end(b).cmp(&q_end(a)))
            .then_with(|| s_end(b).cmp(&s_end(a)))
    });

    let mut i = 0usize;
    while i < hits.len() {
        let j = i + 1;
        if j >= hits.len() {
            i += 1;
            continue;
        }

        let same_context = context(&hits[i]) == context(&hits[j]);
        let same_q_start = q_offset(&hits[i]) == q_offset(&hits[j]);
        let same_s_start = s_offset(&hits[i]) == s_offset(&hits[j]);
        if same_context && same_q_start && same_s_start {
            let mut removed_hit = hits.remove(j);
            if !purge && q_end(&removed_hit) > q_end(&hits[i]) {
                let _ = cut_off_gap_edit_script(
                    &mut removed_hit,
                    q_end(&hits[i]),
                    s_end(&hits[i]),
                    true,
                );
                trimmed_hits.push(removed_hit);
            }
        } else {
            i += 1;
        }
    }

    hits.sort_by(|a, b| {
        context(a)
            .cmp(&context(b))
            .then_with(|| q_end(a).cmp(&q_end(b)))
            .then_with(|| s_end(a).cmp(&s_end(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score))
            .then_with(|| q_offset(b).cmp(&q_offset(a)))
            .then_with(|| s_offset(b).cmp(&s_offset(a)))
    });

    let mut i = 0usize;
    while i < hits.len() {
        let j = i + 1;
        if j >= hits.len() {
            i += 1;
            continue;
        }

        let same_context = context(&hits[i]) == context(&hits[j]);
        let same_q_end = q_end(&hits[i]) == q_end(&hits[j]);
        let same_s_end = s_end(&hits[i]) == s_end(&hits[j]);
        if same_context && same_q_end && same_s_end {
            let mut removed_hit = hits.remove(j);
            if !purge && q_offset(&removed_hit) < q_offset(&hits[i]) {
                let _ = cut_off_gap_edit_script(
                    &mut removed_hit,
                    q_offset(&hits[i]),
                    s_offset(&hits[i]),
                    false,
                );
                trimmed_hits.push(removed_hit);
            }
        } else {
            i += 1;
        }
    }

    let extra_start = hits.len();
    hits.extend(trimmed_hits);
    (hits, extra_start)
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-180
// ```c
// typedef struct BlastHSPList { ... } BlastHSPList;
// typedef struct BlastHitList { ... } BlastHitList;
// ```
pub(crate) fn collect_hits_from_hit_lists(hit_lists: &[Option<BlastpHitList>]) -> Vec<Hit> {
    let mut hits = Vec::new();
    for hit_list_opt in hit_lists {
        let Some(hit_list) = hit_list_opt else {
            continue;
        };
        for hsp_list in hit_list.hsplist_array.iter().take(hit_list.hsplist_count) {
            for hsp in &hsp_list.hsps {
                hits.push(hsp.clone().into_hit());
            }
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(score: i32, q_start: usize, q_end: usize, s_start: usize, s_end: usize) -> Hit {
        Hit {
            identity: 0.0,
            length: q_end.saturating_sub(q_start).saturating_add(1),
            mismatch: 0,
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: 0,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: score,
            gap_info: None,
            num_positives: 0,
        }
    }

    #[test]
    fn test_score_compare_hits_prefers_higher_score_then_offsets() {
        let left = make_hit(100, 10, 30, 20, 40);
        let right = make_hit(90, 5, 25, 10, 30);
        assert_eq!(score_compare_hits(&left, &right), Ordering::Less);
    }

    #[test]
    fn test_purge_hits_with_common_endpoints_keeps_longer_match() {
        let short = make_hit(100, 1, 50, 1, 50);
        let long = make_hit(100, 1, 80, 1, 80);
        let purged = purge_hits_with_common_endpoints(vec![short, long]);
        assert_eq!(purged.len(), 1);
        assert_eq!(purged[0].q_end, 80);
    }

    #[test]
    fn test_hit_list_update_and_prune_keep_best_subjects() {
        let mut hit_list = BlastpHitList::new(2);

        let mut first = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
            best_evalue: 0.0,
        };
        let mut second = BlastpHspList {
            oid: 1,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 2, 11))],
            best_evalue: 0.0,
        };
        let mut third = BlastpHspList {
            oid: 2,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(40, 1, 10, 3, 12))],
            best_evalue: 0.0,
        };

        first.hsps[0].e_value = 1e-10;
        second.hsps[0].e_value = 1e-20;
        third.hsps[0].e_value = 1e-5;

        hit_list.update(first);
        hit_list.update(second);
        hit_list.update(third);
        hit_list.sort_by_evalue();
        hit_list.prune_by_size(2);

        assert_eq!(hit_list.hsplist_count, 2);
        assert_eq!(hit_list.hsplist_array[0].oid, 1);
        assert_eq!(hit_list.hsplist_array[1].oid, 0);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1976-1999
    // ```c
    // for (index = 0; index < hsp_list->hspcnt; index++) {
    //    hsp = hsp_array[index];
    //    if (hsp->evalue > cutoff) {
    //       hsp_array[index] = Blast_HSPFree(hsp_array[index]);
    //    } else {
    //       if (index > hsp_cnt)
    //          hsp_array[hsp_cnt] = hsp_array[index];
    //       hsp_cnt++;
    //    }
    // }
    // hsp_list->hspcnt = hsp_cnt;
    // ```
    #[test]
    fn test_reap_hsplist_by_evalue_preserves_ncbi_compaction_order() {
        let mut hsp_list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![
                BlastpHsp::from_hit(make_hit(80, 1, 10, 1, 10)),
                BlastpHsp::from_hit(make_hit(70, 2, 11, 2, 11)),
                BlastpHsp::from_hit(make_hit(60, 3, 12, 3, 12)),
            ],
            best_evalue: 0.0,
        };
        hsp_list.hsps[0].e_value = 1e-20;
        hsp_list.hsps[1].e_value = 1e-2;
        hsp_list.hsps[2].e_value = 1e-10;

        reap_hsplist_by_evalue(&mut hsp_list, 1e-5);

        assert_eq!(hsp_list.hsps.len(), 2);
        assert_eq!(hsp_list.hsps[0].raw_score, 80);
        assert_eq!(hsp_list.hsps[1].raw_score, 60);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3243-3297
    // ```c
    // evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
    // if (evalue_order < 0) {
    //    Blast_HSPListFree(hsp_list);
    // } else {
    //    s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list);
    // }
    // ```
    #[test]
    fn test_hit_list_update_replaces_equal_worst_with_newer_oid() {
        let mut hit_list = BlastpHitList::new(1);

        let mut first = BlastpHspList {
            oid: 1,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
            best_evalue: 0.0,
        };
        let mut second = BlastpHspList {
            oid: 2,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 2, 11))],
            best_evalue: 0.0,
        };

        first.hsps[0].e_value = 1e-20;
        second.hsps[0].e_value = 1e-20;

        hit_list.update(first);
        hit_list.update(second);
        hit_list.sort_by_evalue();

        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 2);
    }
}
