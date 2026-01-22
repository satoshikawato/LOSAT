use std::cmp::Ordering;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use crate::common::{GapEditOp, Hit};
use crate::report::{write_hit_fields, OutputConfig};

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
//    Int4 hsplist_count; /**< Filled size of the HSP lists array */
//    Int4 hsplist_max; /**< Maximal allowed size of the HSP lists array */
//    double worst_evalue; /**< Highest of the best e-values among the HSP lists */
//    Int4 low_score; /**< The lowest of the best scores among the HSP lists */
//    Boolean heapified; /**< Is this hit list already heapified? */
//    BlastHSPList** hsplist_array; /**< Array of HSP lists for individual database hits */
//    Int4 hsplist_current; /**< Number of allocated HSP list arrays. */
//    Int4 num_hits; /**< Number of similar hits for the query (for mapping) */
// } BlastHitList;
// ```
pub struct BlastnHitList {
    pub hsplist_count: usize,
    pub hsplist_max: usize,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplist_array: Vec<BlastnHspList>,
    pub hsplist_current: usize,
    pub num_hits: usize,
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70 (GetPrelimHitlistSize)
// ```c
// Int4
// GetPrelimHitlistSize(Int4 hitlist_size, Int4 compositionBasedStats, Boolean gapped_calculation)
// {
//     Int4 prelim_hitlist_size = hitlist_size;
//     char * ADAPTIVE_CBS_ENV = getenv("ADAPTIVE_CBS");
//     if (compositionBasedStats) {
//         if(ADAPTIVE_CBS_ENV != NULL) {
//             if(hitlist_size < 1000) {
//                 prelim_hitlist_size = MAX(prelim_hitlist_size + 1000, 1500);
//             }
//             else {
//                 prelim_hitlist_size = prelim_hitlist_size*2 + 50;
//             }
//         }
//         else {
//             if(hitlist_size <= 500) {
//                 prelim_hitlist_size = 1050;
//             }
//             else {
//                 prelim_hitlist_size = prelim_hitlist_size*2 + 50;
//             }
//         }
//     }
//     else if (gapped_calculation) {
//          prelim_hitlist_size = MIN(MAX(2 * prelim_hitlist_size, 10),
//                                   prelim_hitlist_size + 50);
//     }
//     return prelim_hitlist_size;
// }
// ```
pub fn get_prelim_hitlist_size(
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1414-1435
// ```c
// static int
// s_EvalueCompareHSPs(const void* v1, const void* v2)
// {
//    BlastHSP* h1,* h2;
//    int retval = 0;
//    ...
//    if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
//       return retval;
//    return ScoreCompareHSPs(v1, v2);
// }
// ```
pub fn evalue_compare_hsps(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
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
//     if (hsp_list->hspcnt > 1) {
//         Int4 index;
//         BlastHSP** hsp_array = hsp_list->hsp_array;
//         for (index = 0; index < hsp_list->hspcnt - 1; ++index) {
//             if (s_EvalueCompareHSPs(&hsp_array[index],
//                                     &hsp_array[index+1]) > 0) {
//                 break;
//             }
//         }
//         if (index < hsp_list->hspcnt - 1) {
//             qsort(hsp_list->hsp_array, hsp_list->hspcnt,
//                   sizeof(BlastHSP*), s_EvalueCompareHSPs);
//         }
//     }
// }
// ```
pub fn sort_hsplist_by_evalue(list: &mut BlastnHspList) {
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
//     double best_evalue = (double) INT4_MAX;
//     for (index=0; index<hsp_list->hspcnt; index++)
//         best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
//     return best_evalue;
// }
// ```
pub fn update_best_evalue(list: &mut BlastnHspList) {
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3077-3106
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1627-1650
// ```c
// static void
// s_Heapify (char* base0, char* base, char* lim, char* last,
//            size_t width, int (*compar )(const void*, const void* ))
// {
//    ...
//    left_son = base0 + 2*(base-base0) + width;
//    while (base <= lim) {
//       ...
//       large_son = (*compar)(left_son, left_son+width) >= 0 ?
//          left_son : left_son+width;
//       if ((*compar)(base, large_son) < 0) {
//          ...
//       } else
//          break;
//    }
// }
// ```
fn heapify_hsplist_array(lists: &mut [BlastnHspList], start: usize, end: usize) {
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
// static void
// s_CreateHeap (void* b, size_t nel, size_t width,
//    int (*compar )(const void*, const void* ))
// {
//    if (nel < 2)
//       return;
//    ...
//    i = nel/2;
//    for (base = &base0[(i - 1)*width]; i > 0; base = base - width) {
//       s_Heapify(base0, base, lim, basef, width, compar);
//       i--;
//    }
// }
// ```
fn create_hsplist_heap(lists: &mut [BlastnHspList]) {
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

impl BlastnHitList {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3125-3133
    // ```c
    // BlastHitList* Blast_HitListNew(Int4 hitlist_size)
    // {
    //    BlastHitList* new_hitlist = (BlastHitList*) calloc(1, sizeof(BlastHitList));
    //    new_hitlist->hsplist_max = hitlist_size;
    //    new_hitlist->low_score = INT4_MAX;
    //    new_hitlist->hsplist_count = 0;
    //    new_hitlist->hsplist_current = 0;
    //    return new_hitlist;
    // }
    // ```
    pub fn new(hitlist_size: usize) -> Self {
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
    //     const int kStartValue = 100;
    //     if (hit_list->hsplist_current >= hit_list->hsplist_max)
    //        return 1;
    //     if (hit_list->hsplist_current <= 0)
    //        hit_list->hsplist_current = kStartValue;
    //     else
    //        hit_list->hsplist_current =
    //           MIN(2*hit_list->hsplist_current, hit_list->hsplist_max);
    //     hit_list->hsplist_array =
    //        (BlastHSPList**) realloc(hit_list->hsplist_array,
    //                                 hit_list->hsplist_current*sizeof(BlastHSPList*));
    //     if (hit_list->hsplist_array == NULL)
    //        return BLASTERR_MEMORY;
    //     return 0;
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
    // static void
    // s_BlastHitListInsertHSPListInHeap(BlastHitList* hit_list,
    //                                  BlastHSPList* hsp_list)
    // {
    //       Blast_HSPListFree(hit_list->hsplist_array[0]);
    //       hit_list->hsplist_array[0] = hsp_list;
    //       if (hit_list->hsplist_count >= 2) {
    //          s_Heapify((char*)hit_list->hsplist_array, (char*)hit_list->hsplist_array,
    //                  (char*)&hit_list->hsplist_array[hit_list->hsplist_count/2 - 1],
    //                  (char*)&hit_list->hsplist_array[hit_list->hsplist_count-1],
    //                  sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //       }
    //       hit_list->worst_evalue = hit_list->hsplist_array[0]->best_evalue;
    //       hit_list->low_score = hit_list->hsplist_array[0]->hsp_array[0]->score;
    // }
    // ```
    fn insert_hsplist_in_heap(&mut self, hsp_list: BlastnHspList) {
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
    //    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //    if (hit_list->hsplist_count < hit_list->hsplist_max) {
    //       if (hit_list->hsplist_current == hit_list->hsplist_count)
    //       {
    //          Int2 status = s_Blast_HitListGrowHSPListArray(hit_list);
    //          if (status)
    //            return status;
    //       }
    //       hit_list->hsplist_array[hit_list->hsplist_count++] = hsp_list;
    //       hit_list->worst_evalue =
    //          MAX(hsp_list->best_evalue, hit_list->worst_evalue);
    //       hit_list->low_score =
    //          MIN(hsp_list->hsp_array[0]->score, hit_list->low_score);
    //    } else {
    //       if (!hit_list->heapified) {
    //           for (index =0; index < hit_list->hsplist_count; index++) {
    //               Blast_HSPListSortByEvalue(hit_list->hsplist_array[index]);
    //               hit_list->hsplist_array[index]->best_evalue =
    //                   s_BlastGetBestEvalue(hit_list->hsplist_array[index]);
    //           }
    //           s_CreateHeap(hit_list->hsplist_array, hit_list->hsplist_count,
    //                        sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //           hit_list->heapified = TRUE;
    //       }
    //       Blast_HSPListSortByEvalue(hsp_list);
    //       hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
    //       if (evalue_order < 0) {
    //          Blast_HSPListFree(hsp_list);
    //       } else {
    //          s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list);
    //       }
    //    }
    //    return 0;
    // }
    // ```
    pub fn update(&mut self, mut hsp_list: BlastnHspList) {
        update_best_evalue(&mut hsp_list);

        if self.hsplist_count < self.hsplist_max {
            if self.hsplist_current == self.hsplist_count && !self.grow_hsplist_array() {
                return;
            }
            self.hsplist_array.push(hsp_list);
            self.hsplist_count += 1;
            self.worst_evalue = self.worst_evalue.max(self.hsplist_array.last().unwrap().best_evalue);
            if let Some(score) = self.hsplist_array.last().unwrap().hsps.first().map(|h| h.raw_score) {
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
    //    if (!hit_list) return;
    //    hsplist_count = hit_list->hsplist_count;
    //    for (index = 0; index < hsplist_count &&
    //            hit_list->hsplist_array[index]->hspcnt > 0; ++index);
    //    hit_list->hsplist_count = index;
    //    for ( ; index < hsplist_count; ++index) {
    //       Blast_HSPListFree(hit_list->hsplist_array[index]);
    //    }
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
    //    if (hit_list && hit_list->hsplist_count > 1) {
    //       qsort(hit_list->hsplist_array, hit_list->hsplist_count,
    //             sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //    }
    //    s_BlastHitListPurge(hit_list);
    //    return 0;
    // }
    // ```
    pub fn sort_by_evalue(&mut self) {
        if self.hsplist_count > 1 {
            self.hsplist_array.sort_by(compare_hsp_lists);
        }
        self.purge();
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
    pub fn prune_by_size(&mut self, hitlist_size: usize) {
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

/// Write BLASTN output in NCBI HSP list order without regrouping.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
/// ```c
/// int ScoreCompareHSPs(const void* h1, const void* h2) {
///    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
///        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
///        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
///        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
///        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
///    }
/// }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3077-3106
/// ```c
/// static int s_EvalueCompareHSPLists(const void* v1, const void* v2) {
///    if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0)
///       return retval;
///    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
///    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
///    return BLAST_CMP(h2->oid, h1->oid);
/// }
/// ```
pub fn write_output_blastn_hitlists(
    hit_lists: &[Option<BlastnHitList>],
    out_path: Option<&PathBuf>,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
) -> io::Result<()> {
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };

    let config = OutputConfig::ncbi_compat();

    for hit_list_opt in hit_lists.iter() {
        let hit_list = match hit_list_opt {
            Some(value) => value,
            None => continue,
        };
        for hsp_list in &hit_list.hsplist_array {
            for hsp in &hsp_list.hsps {
                let query_id = query_ids
                    .get(hsp.q_idx as usize)
                    .map(|id| id.as_ref())
                    .unwrap_or("unknown");
                let subject_id = subject_ids
                    .get(hsp.s_idx as usize)
                    .map(|id| id.as_ref())
                    .unwrap_or("unknown");
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
                write_hit_fields(
                    &mut writer,
                    query_id,
                    subject_id,
                    hsp.identity,
                    hsp.length,
                    hsp.mismatch,
                    hsp.gapopen,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.e_value,
                    hsp.bit_score,
                    &config,
                )?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    // } BlastHSP;
    // ```
    fn make_hsp(e_value: f64, raw_score: i32, s_idx: u32) -> BlastnHsp {
        BlastnHsp {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 1,
            s_start: 1,
            s_end: 1,
            e_value,
            bit_score: 0.0,
            query_frame: 1,
            query_length: 1,
            q_idx: 0,
            s_idx,
            raw_score,
            gap_info: None,
        }
    }

    #[test]
    fn test_prelim_hitlist_size_gapped() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70
        // ```c
        // else if (gapped_calculation) {
        //      prelim_hitlist_size = MIN(MAX(2 * prelim_hitlist_size, 10),
        //                               prelim_hitlist_size + 50);
        // }
        // ```
        assert_eq!(get_prelim_hitlist_size(1, false, true), 10);
        assert_eq!(get_prelim_hitlist_size(30, false, true), 60);
        assert_eq!(get_prelim_hitlist_size(1000, false, true), 1050);
    }

    #[test]
    fn test_hitlist_update_keeps_best() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3243-3297
        // ```c
        // if (hit_list->hsplist_count < hit_list->hsplist_max) { ... }
        // else {
        //    if (!hit_list->heapified) { ... }
        //    evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
        //    if (evalue_order < 0) { Blast_HSPListFree(hsp_list); }
        //    else { s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list); }
        // }
        // ```
        let mut hit_list = BlastnHitList::new(1);

        let hsp_list_a = BlastnHspList {
            oid: 10,
            query_index: 0,
            hsps: vec![make_hsp(5.0, 50, 10)],
            best_evalue: i32::MAX as f64,
        };
        hit_list.update(hsp_list_a);
        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 10);

        let hsp_list_b = BlastnHspList {
            oid: 20,
            query_index: 0,
            hsps: vec![make_hsp(1.0, 80, 20)],
            best_evalue: i32::MAX as f64,
        };
        hit_list.update(hsp_list_b);
        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 20);

        let hsp_list_c = BlastnHspList {
            oid: 30,
            query_index: 0,
            hsps: vec![make_hsp(10.0, 40, 30)],
            best_evalue: i32::MAX as f64,
        };
        hit_list.update(hsp_list_c);
        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 20);
    }
}
