//! BLASTP-specific HSP list and hit list helpers.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c

use std::cmp::Ordering;

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
    pub gapped_q_start: i32,
    pub gapped_s_start: i32,
    pub e_value: f64,
    pub bit_score: f64,
    pub num_ident: usize,
    pub query_context: i32,
    pub query_frame: i32,
    pub subject_frame: i32,
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
            gapped_q_start: i32::try_from(q_start.saturating_sub(1))
                .expect("NCBI BLAST HSP query.gapped_start must fit in Int4"),
            gapped_s_start: i32::try_from(s_start.saturating_sub(1))
                .expect("NCBI BLAST HSP subject.gapped_start must fit in Int4"),
            e_value,
            bit_score,
            num_ident,
            query_context: i32::try_from(q_idx).expect("NCBI BLAST HSP context must fit in Int4"),
            query_frame,
            subject_frame: 0,
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
#[derive(Debug, Clone)]
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:67-76
// ```c
// typedef struct BlastCompo_HeapRecord {
//     double        bestEvalue;
//     int           bestScore;
//     int           subject_index;
//     void *        theseAlignments;
// } BlastCompo_HeapRecord;
// ```
#[derive(Debug, Clone)]
pub(crate) struct BlastCompoHeapRecord {
    pub best_evalue: f64,
    pub best_score: i32,
    pub subject_index: u32,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:67-76
    // ```c
    // void * theseAlignments;
    // ```
    //
    // Rust stores the typed `BlastpHspList` directly instead of a raw
    // `void*`, preserving the same ownership and comparison behavior.
    pub hsplist: Option<BlastpHspList>,
}

impl Default for BlastCompoHeapRecord {
    fn default() -> Self {
        Self {
            best_evalue: 0.0,
            best_score: 0,
            subject_index: 0,
            hsplist: None,
        }
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/composition_adjustment/compo_heap.h:79-95
// ```c
// typedef struct BlastCompo_Heap {
//     int n;
//     int capacity;
//     int heapThreshold;
//     double ecutoff;
//     double worstEvalue;
//     struct BlastCompo_HeapRecord *array;
//     struct BlastCompo_HeapRecord *heapArray;
// } BlastCompo_Heap;
// ```
#[derive(Debug, Clone)]
pub(crate) struct BlastCompoHeap {
    pub n: usize,
    pub capacity: usize,
    pub heap_threshold: usize,
    pub ecutoff: f64,
    pub worst_evalue: f64,
    pub array: Option<Vec<BlastCompoHeapRecord>>,
    pub heap_array: Option<Vec<BlastCompoHeapRecord>>,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:50-57
// ```c
// #define HEAP_INITIAL_CAPACITY 100
// #define HEAP_RESIZE_FACTOR 1.5
// #define HEAP_MIN_RESIZE 100
// ```
const COMPO_HEAP_INITIAL_CAPACITY: usize = 100;
const COMPO_HEAP_MIN_RESIZE: usize = 100;
const COMPO_HEAP_RESIZE_FACTOR: f64 = 1.5;

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:53-54
// ```c
// /** by what factor might initially reported E-value exceed true E-value */
// #define EVALUE_STRETCH 5
// ```
const COMPO_HEAP_EVALUE_STRETCH: f64 = 5.0;

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:81-93
// ```c
// static int
// s_CompoHeapRecordCompare(BlastCompo_HeapRecord * place1,
//                          BlastCompo_HeapRecord * place2)
// {
//     int result;
//     if (0 == (result = CMP(place1->bestEvalue, place2->bestEvalue)) &&
//         0 == (result = CMP(place2->bestScore, place1->bestScore))) {
//         result = CMP(place2->subject_index, place1->subject_index);
//     }
//     return result > 0;
// }
// ```
fn compo_heap_record_compare(place1: &BlastCompoHeapRecord, place2: &BlastCompoHeapRecord) -> bool {
    let mut result = if place1.best_evalue > place2.best_evalue {
        1
    } else if place1.best_evalue < place2.best_evalue {
        -1
    } else {
        0
    };
    if result == 0 {
        result = match place2.best_score.cmp(&place1.best_score) {
            Ordering::Less => -1,
            Ordering::Equal => 0,
            Ordering::Greater => 1,
        };
    }
    if result == 0 {
        result = match place2.subject_index.cmp(&place1.subject_index) {
            Ordering::Less => -1,
            Ordering::Equal => 0,
            Ordering::Greater => 1,
        };
    }
    result > 0
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:154-186
// ```c
// static void
// s_CompoHeapifyDown(BlastCompo_HeapRecord * heapArray, int top, int n)
// {
//     ...
// }
// ```
fn compo_heapify_down(heap_array: &mut [BlastCompoHeapRecord], top: usize, n: usize) {
    let mut largest = top;
    loop {
        let i = largest;
        let left = 2 * i;
        let right = 2 * i + 1;
        largest = if left <= n && compo_heap_record_compare(&heap_array[left], &heap_array[i]) {
            left
        } else {
            i
        };
        if right <= n && compo_heap_record_compare(&heap_array[right], &heap_array[largest]) {
            largest = right;
        }
        if largest == i {
            break;
        }
        heap_array.swap(i, largest);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:197-214
// ```c
// static void
// s_CompoHeapifyUp(BlastCompo_HeapRecord * heapArray, int i)
// {
//     ...
// }
// ```
fn compo_heapify_up(heap_array: &mut [BlastCompoHeapRecord], mut i: usize) {
    let mut parent = i / 2;
    while parent >= 1 && compo_heap_record_compare(&heap_array[i], &heap_array[parent]) {
        heap_array.swap(i, parent);
        i = parent;
        parent /= 2;
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1355-1381
// ```c
// Boolean Blast_HSPListIsSortedByScore(const BlastHSPList* hsp_list)
// {
//     ...
//     if (ScoreCompareHSPs(&hsp_list->hsp_array[index],
//                          &hsp_list->hsp_array[index+1]) > 0) {
//         return FALSE;
//     }
// }
//
// void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
// {
//     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
//     }
// }
// ```
pub(crate) fn sort_hsplist_by_score(list: &mut BlastpHspList) {
    if list.hsps.len() <= 1 {
        return;
    }

    let mut index = 0usize;
    while index < list.hsps.len() - 1 {
        if score_compare_hsps(&list.hsps[index], &list.hsps[index + 1]) == Ordering::Greater {
            break;
        }
        index += 1;
    }
    if index < list.hsps.len() - 1 {
        list.hsps.sort_unstable_by(score_compare_hsps);
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
            list.hsps.sort_unstable_by(evalue_compare_hsps);
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

impl BlastCompoHeap {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:361-375
    // ```c
    // int
    // BlastCompo_HeapInitialize(BlastCompo_Heap * self, int heapThreshold,
    //                           double ecutoff)
    // {
    //     self->n             = 0;
    //     self->heapThreshold = heapThreshold;
    //     self->ecutoff       = ecutoff;
    //     self->heapArray     = NULL;
    //     self->capacity      = MIN(HEAP_INITIAL_CAPACITY, heapThreshold);
    //     self->worstEvalue   = 0;
    //     self->array = calloc(self->capacity + 1, sizeof(BlastCompo_HeapRecord));
    // }
    // ```
    pub(crate) fn new(heap_threshold: usize, ecutoff: f64) -> Self {
        let capacity = COMPO_HEAP_INITIAL_CAPACITY.min(heap_threshold);
        let mut array = Vec::with_capacity(capacity.saturating_add(1));
        array.push(BlastCompoHeapRecord::default());
        Self {
            n: 0,
            capacity,
            heap_threshold,
            ecutoff,
            worst_evalue: 0.0,
            array: Some(array),
            heap_array: None,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:221-238
    // ```c
    // static void
    // s_ConvertToHeap(BlastCompo_Heap * self)
    // {
    //     ...
    // }
    // ```
    fn convert_to_heap(&mut self) {
        let Some(mut heap_array) = self.array.take() else {
            return;
        };
        while heap_array.len() <= self.n {
            heap_array.push(BlastCompoHeapRecord::default());
        }
        if self.n >= 2 {
            for i in (1..=self.n / 2).rev() {
                compo_heapify_down(&mut heap_array, i, self.n);
            }
        }
        self.heap_array = Some(heap_array);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:257-284
    // ```c
    // int
    // BlastCompo_HeapWouldInsert(BlastCompo_Heap * self,
    //                            double eValue,
    //                            int score,
    //                            int subject_index)
    // {
    //     ...
    // }
    // ```
    pub(crate) fn would_insert(&mut self, evalue: f64, score: i32, subject_index: u32) -> bool {
        if self.n < self.heap_threshold || evalue <= self.ecutoff || evalue < self.worst_evalue {
            return true;
        }
        if self.heap_array.is_none() {
            self.convert_to_heap();
        }
        let mut record = BlastCompoHeapRecord::default();
        record.best_evalue = evalue;
        record.best_score = score;
        record.subject_index = subject_index;
        let heap_array = self
            .heap_array
            .as_ref()
            .expect("compo heap converted before root comparison");
        compo_heap_record_compare(&heap_array[1], &record)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:294-352
    // ```c
    // int
    // BlastCompo_HeapInsert(BlastCompo_Heap * self,
    //                       void * alignments,
    //                       double eValue,
    //                       int score,
    //                       int subject_index,
    //                       void ** discardedAlignments)
    // {
    //     ...
    // }
    // ```
    pub(crate) fn insert(
        &mut self,
        hsplist: BlastpHspList,
        evalue: f64,
        score: i32,
        subject_index: u32,
    ) {
        if self.array.is_some() && self.n >= self.heap_threshold {
            self.convert_to_heap();
        }
        if let Some(array) = self.array.as_mut() {
            self.n += 1;
            while array.len() <= self.n {
                array.push(BlastCompoHeapRecord::default());
            }
            array[self.n] = BlastCompoHeapRecord {
                best_evalue: evalue,
                best_score: score,
                subject_index,
                hsplist: Some(hsplist),
            };
            if self.worst_evalue < evalue {
                self.worst_evalue = evalue;
            }
            return;
        }

        let heap_array = self
            .heap_array
            .as_mut()
            .expect("compo heap storage exists when inserting");
        if self.n < self.heap_threshold
            || (evalue <= self.ecutoff && self.worst_evalue <= self.ecutoff)
        {
            self.n += 1;
            while heap_array.len() <= self.n {
                heap_array.push(BlastCompoHeapRecord::default());
            }
            heap_array[self.n] = BlastCompoHeapRecord {
                best_evalue: evalue,
                best_score: score,
                subject_index,
                hsplist: Some(hsplist),
            };
            compo_heapify_up(heap_array, self.n);
        } else {
            let candidate = BlastCompoHeapRecord {
                best_evalue: evalue,
                best_score: score,
                subject_index,
                hsplist: Some(hsplist),
            };
            if compo_heap_record_compare(&heap_array[1], &candidate) {
                heap_array[1] = candidate;
            }
            compo_heapify_down(heap_array, 1, self.n);
        }
        self.worst_evalue = heap_array[1].best_evalue;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:356-359
    // ```c
    // int
    // BlastCompo_HeapFilledToCutoff(const BlastCompo_Heap * self)
    // {
    //     return self->n >= self->heapThreshold &&
    //         self->worstEvalue <= self->ecutoff;
    // }
    // ```
    pub(crate) fn filled_to_cutoff(&self) -> bool {
        self.n >= self.heap_threshold && self.worst_evalue <= self.ecutoff
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:391-414
    // ```c
    // void *
    // BlastCompo_HeapPop(BlastCompo_Heap * self)
    // {
    //     ...
    // }
    // ```
    pub(crate) fn pop(&mut self) -> Option<BlastpHspList> {
        if self.heap_array.is_none() {
            self.convert_to_heap();
        }
        let heap_array = self
            .heap_array
            .as_mut()
            .expect("compo heap storage exists when popping");
        if self.n == 0 {
            return None;
        }
        let results = heap_array[1].hsplist.take();
        let last = std::mem::take(&mut heap_array[self.n]);
        self.n -= 1;
        if self.n > 0 {
            heap_array[1] = last;
            compo_heapify_down(heap_array, 1, self.n);
            self.worst_evalue = heap_array[1].best_evalue;
        } else {
            self.worst_evalue = 0.0;
        }
        results
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/redo_alignment.c:1561-1576
// ```c
// int
// BlastCompo_EarlyTermination(double evalue,
//                             BlastCompo_Heap significantMatches[],
//                             int numQueries)
// {
//     ...
// }
// ```
pub(crate) fn blast_compo_early_termination(
    evalue: f64,
    significant_matches: &[BlastCompoHeap],
) -> bool {
    for heap in significant_matches {
        if heap.filled_to_cutoff() {
            let ecutoff = heap.ecutoff;
            if evalue <= COMPO_HEAP_EVALUE_STRETCH * ecutoff {
                return false;
            }
        } else {
            return false;
        }
    }
    true
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2494-2514
// ```c
// static void
// s_FillResultsFromCompoHeaps(BlastHSPResults * results,
//                             BlastCompo_Heap heaps[],
//                             Int4 hitlist_size)
// {
//     ...
// }
// ```
pub(crate) fn fill_results_from_compo_heaps(
    hitlist_size: usize,
    heaps: &mut [BlastCompoHeap],
) -> Vec<Option<BlastpHitList>> {
    let mut hit_lists: Vec<Option<BlastpHitList>> = heaps
        .iter()
        .map(|_| Some(BlastpHitList::new(hitlist_size)))
        .collect();
    for (query_index, heap) in heaps.iter_mut().enumerate() {
        let hit_list = hit_lists[query_index]
            .as_mut()
            .expect("hit list initialized for each query");
        while let Some(hsp_list) = heap.pop() {
            hit_list.update(hsp_list);
        }
    }
    reverse_hit_list_order(&mut hit_lists);
    hit_lists
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3420-3438
// ```c
// Int2 Blast_HSPResultsReverseOrder(BlastHSPResults* results)
// {
//    for (index = 0; index < results->num_queries; ++index) {
//       hit_list = results->hitlist_array[index];
//       if (hit_list && hit_list->hsplist_count > 1) {
//          for (index1 = 0; index1 < hit_list->hsplist_count/2; ++index1) {
//             ...
//          }
//       }
//    }
// }
// ```
fn reverse_hit_list_order(hit_lists: &mut [Option<BlastpHitList>]) {
    for hit_list_opt in hit_lists {
        let Some(hit_list) = hit_list_opt.as_mut() else {
            continue;
        };
        if hit_list.hsplist_count > 1 {
            hit_list.hsplist_array[..hit_list.hsplist_count].reverse();
        }
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
            self.hsplist_array[..self.hsplist_count].sort_unstable_by(compare_hsp_lists);
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
//    hsp_array = hsp_list->hsp_array;
//    hsp_count = hsp_list->hspcnt;
//    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
//    ...
//    qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
//    ...
//    return hsp_count;
// }
// ```
#[allow(dead_code)]
pub(crate) fn purge_hits_with_common_endpoints_ex(
    hits: Vec<Hit>,
    purge: bool,
) -> (Vec<Hit>, usize) {
    purge_hits_for_subject_ex(hits, purge)
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2315
// ```c
// static int
// s_QueryOffsetCompareHSPs(const void* v1, const void* v2)
// {
//    ...
//    if (h1->query.offset < h2->query.offset) return -1;
//    if (h1->subject.offset < h2->subject.offset) return -1;
//    if (h1->score < h2->score) return 1;
//    if (h1->query.end < h2->query.end) return 1;
//    if (h1->subject.end < h2->subject.end) return 1;
// }
// ```
fn hit_query_offset_compare(a: &Hit, b: &Hit) -> Ordering {
    (a.q_idx, a.query_frame)
        .cmp(&(b.q_idx, b.query_frame))
        .then_with(|| {
            a.q_start
                .saturating_sub(1)
                .cmp(&b.q_start.saturating_sub(1))
        })
        .then_with(|| {
            a.s_start
                .min(a.s_end)
                .saturating_sub(1)
                .cmp(&b.s_start.min(b.s_end).saturating_sub(1))
        })
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.q_end.cmp(&a.q_end))
        .then_with(|| b.s_start.max(b.s_end).cmp(&a.s_start.max(a.s_end)))
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2333-2379
// ```c
// static int
// s_QueryEndCompareHSPs(const void* v1, const void* v2)
// {
//    ...
//    if (h1->query.end < h2->query.end) return -1;
//    if (h1->subject.end < h2->subject.end) return -1;
//    if (h1->score < h2->score) return 1;
//    if (h1->query.offset < h2->query.offset) return 1;
//    if (h1->subject.offset < h2->subject.offset) return 1;
// }
// ```
fn hit_query_end_compare(a: &Hit, b: &Hit) -> Ordering {
    (a.q_idx, a.query_frame)
        .cmp(&(b.q_idx, b.query_frame))
        .then_with(|| a.q_end.cmp(&b.q_end))
        .then_with(|| a.s_start.max(a.s_end).cmp(&b.s_start.max(b.s_end)))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| {
            b.q_start
                .saturating_sub(1)
                .cmp(&a.q_start.saturating_sub(1))
        })
        .then_with(|| {
            b.s_start
                .min(b.s_end)
                .saturating_sub(1)
                .cmp(&a.s_start.min(a.s_end).saturating_sub(1))
        })
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // ```
    hits.sort_unstable_by(hit_query_offset_compare);

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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2504
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
    // ```
    hits.sort_unstable_by(hit_query_end_compare);

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

    fn make_blastp_hsp(
        score: i32,
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
    ) -> BlastpHsp {
        BlastpHsp::from_hit(make_hit(score, q_start, q_end, s_start, s_end))
    }

    #[test]
    fn test_score_compare_hits_prefers_higher_score_then_offsets() {
        let left = make_hit(100, 10, 30, 20, 40);
        let right = make_hit(90, 5, 25, 10, 30);
        assert_eq!(score_compare_hits(&left, &right), Ordering::Less);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1355-1377
    // ```c
    // if (!Blast_HSPListIsSortedByScore(hsp_list)) {
    //     qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //           ScoreCompareHSPs);
    // }
    // ```
    #[test]
    fn test_sort_hsplist_by_score_matches_ncbi_tie_breakers() {
        let mut list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![
                make_blastp_hsp(100, 20, 40, 30, 50),
                make_blastp_hsp(120, 15, 35, 20, 40),
                make_blastp_hsp(100, 10, 25, 20, 35),
                make_blastp_hsp(100, 20, 30, 20, 30),
            ],
            best_evalue: 0.0,
        };

        sort_hsplist_by_score(&mut list);

        let ordered: Vec<(i32, usize, usize, usize, usize)> = list
            .hsps
            .iter()
            .map(|hsp| {
                (
                    hsp.raw_score,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.q_start,
                    hsp.q_end,
                )
            })
            .collect();
        assert_eq!(
            ordered,
            vec![
                (120, 20, 40, 15, 35),
                (100, 20, 35, 10, 25),
                (100, 20, 30, 20, 30),
                (100, 30, 50, 20, 40),
            ]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1414-1455
    // ```c
    // if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
    //    return retval;
    // return ScoreCompareHSPs(v1, v2);
    // ```
    #[test]
    fn test_sort_hsplist_by_evalue_uses_score_compare_as_tie_breaker() {
        let mut list = BlastpHspList {
            oid: 0,
            query_index: 0,
            hsps: vec![
                make_blastp_hsp(90, 20, 40, 30, 50),
                make_blastp_hsp(110, 15, 35, 20, 40),
                make_blastp_hsp(110, 10, 25, 20, 35),
            ],
            best_evalue: 0.0,
        };
        list.hsps[0].e_value = 1e-20;
        list.hsps[1].e_value = 1e-40;
        list.hsps[2].e_value = 1e-40;

        sort_hsplist_by_evalue(&mut list);

        let ordered: Vec<(f64, i32, usize, usize)> = list
            .hsps
            .iter()
            .map(|hsp| (hsp.e_value, hsp.raw_score, hsp.s_start, hsp.q_start))
            .collect();
        assert_eq!(
            ordered,
            vec![
                (1e-40, 110, 20, 10),
                (1e-40, 110, 20, 15),
                (1e-20, 90, 30, 20),
            ]
        );
    }

    #[test]
    fn test_purge_hits_with_common_endpoints_keeps_longer_match() {
        let short = make_hit(100, 1, 50, 1, 50);
        let long = make_hit(100, 1, 80, 1, 80);
        let purged = purge_hits_with_common_endpoints(vec![short, long]);
        assert_eq!(purged.len(), 1);
        assert_eq!(purged[0].q_end, 80);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2475-2534
    // ```c
    // hsp_array = hsp_list->hsp_array;
    // hsp_count = hsp_list->hspcnt;
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // ...
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
    // ```
    #[test]
    fn test_purge_hits_with_common_endpoints_ex_keeps_single_hsp_list_order_without_regrouping() {
        let left = make_hit(110, 5, 40, 50, 85);
        let right = make_hit(100, 1, 30, 10, 39);

        let (purged, extra_start) = purge_hits_with_common_endpoints_ex(vec![left, right], false);

        assert_eq!(purged.len(), 2);
        assert_eq!(extra_start, 2);
        assert_eq!(purged[0].raw_score, 100);
        assert_eq!(purged[1].raw_score, 110);
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:257-284
    // ```c
    // BlastCompo_HeapWouldInsert(...);
    // ```
    #[test]
    fn test_compo_heap_would_insert_uses_ncbi_score_and_subject_tie_breakers() {
        let mut heap = BlastCompoHeap::new(1, 1e-50);
        let mut first = BlastpHspList {
            oid: 1,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
            best_evalue: 1e-20,
        };
        first.hsps[0].e_value = 1e-20;
        heap.insert(first, 1e-20, 50, 1);

        assert!(heap.would_insert(1e-20, 50, 2));
        assert!(!heap.would_insert(1e-20, 50, 0));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:391-414
    // ```c
    // void * BlastCompo_HeapPop(BlastCompo_Heap * self)
    // {
    //     ...
    // }
    // ```
    #[test]
    fn test_compo_heap_pop_returns_worst_record_first() {
        let mut heap = BlastCompoHeap::new(2, 1e-50);

        let mut better = BlastpHspList {
            oid: 7,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(70, 1, 10, 1, 10))],
            best_evalue: 1e-20,
        };
        let mut worse = BlastpHspList {
            oid: 3,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 2, 11))],
            best_evalue: 1e-10,
        };
        better.hsps[0].e_value = 1e-20;
        worse.hsps[0].e_value = 1e-10;

        heap.insert(better, 1e-20, 70, 7);
        heap.insert(worse, 1e-10, 60, 3);

        let first = heap.pop().expect("worst record present");
        let second = heap.pop().expect("best record present");

        assert_eq!(first.oid, 3);
        assert_eq!(second.oid, 7);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:2494-2514
    // ```c
    // while (NULL != (hsp_list = BlastCompo_HeapPop(heap))) {
    //     Blast_HitListUpdate(hitlist, hsp_list);
    // }
    // ```
    #[test]
    fn test_fill_results_from_compo_heaps_builds_hitlist_via_hitlist_update() {
        let mut heaps = vec![BlastCompoHeap::new(2, 1e-50)];

        let mut first = BlastpHspList {
            oid: 1,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(50, 1, 10, 1, 10))],
            best_evalue: 1e-10,
        };
        let mut second = BlastpHspList {
            oid: 2,
            query_index: 0,
            hsps: vec![BlastpHsp::from_hit(make_hit(60, 1, 10, 2, 11))],
            best_evalue: 1e-20,
        };
        first.hsps[0].e_value = 1e-10;
        second.hsps[0].e_value = 1e-20;

        heaps[0].insert(first, 1e-10, 50, 1);
        heaps[0].insert(second, 1e-20, 60, 2);

        let mut hit_lists = fill_results_from_compo_heaps(2, &mut heaps);
        let hit_list = hit_lists[0].as_mut().expect("query hit list initialized");
        hit_list.sort_by_evalue();

        assert_eq!(hit_list.hsplist_count, 2);
        assert_eq!(hit_list.hsplist_array[0].oid, 2);
        assert_eq!(hit_list.hsplist_array[1].oid, 1);
    }
}
