//! BLASTP-specific HSP list and hit list helpers.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c

use std::cmp::Ordering;

use crate::common::{GapEditOp, Hit};

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
#[derive(Debug, Clone)]
pub(crate) struct BlastpHsp {
    pub identity: f64,
    pub length: usize,
    /// The columns less the identities and the gaps (`Hit::mismatch`), fewer than none when
    /// stored gap characters of the subject face gaps of the query.
    pub mismatch: i64,
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
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:139-141
    /// `BlastHSP::comp_adjustment_method`: 0 (no adjustment), 1 (composition-based
    /// statistics) or 2 (compositional matrix adjustment), set from the redone
    /// alignment's matrix adjustment rule (`kappa.rs`).
    pub comp_adjustment_method: u8,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:135-138
    /// The size of the HSP's linked set with sum statistics (psiblast's `-sum_stats true`,
    /// `core::link_hsps`), 1 otherwise; above 1 the report shows it in parentheses after `Expect`.
    pub num: i32,
}

impl BlastpHsp {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
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
            ..
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
            comp_adjustment_method: 0,
            num: 1,
        }
    }

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
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
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: self.gap_info,
            num_positives: self.num_positives,
        }
    }
}

// No NCBI counterpart: NCBI keeps one BlastHSP type in every HSP list. LOSAT's preliminary
// hit lists hold a smaller element (`x_prelim_compact.rs`);
// this trait gives the hit-list helpers the fields they read, in `BlastpHsp`'s convention.
/// The fields of an HSP that the HSP-list and hit-list helpers read: the E-value, the raw score
/// and the offsets of `ScoreCompareHSPs`.
pub(crate) trait HitListHsp {
    fn e_value(&self) -> f64;
    fn raw_score(&self) -> i32;
    /// `(q_start, q_end, s_start, s_end)` as `BlastpHsp` stores them: 1-based starts, ends
    /// inclusive (NCBI's 0-based offset + 1, and its exclusive end).
    fn offsets(&self) -> (usize, usize, usize, usize);
    /// The length of the query (or query part, or searched range) the HSP was found in:
    /// `contexts[hsp->context].query_length`.
    fn query_length(&self) -> usize;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:508-511 (cited in
    // full at `query_split::merge_query_chunk`)
    /// Moves the HSP to query `q_idx` (a blastp query has one context, its index), whose letters
    /// start `offset` letters before the HSP's query coordinates, of `query_length` letters
    /// (`query_split::rebase_hit_list`).
    fn move_to_query(&mut self, q_idx: u32, offset: i32, query_length: usize);
}

impl HitListHsp for BlastpHsp {
    #[inline]
    fn e_value(&self) -> f64 {
        self.e_value
    }
    #[inline]
    fn raw_score(&self) -> i32 {
        self.raw_score
    }
    #[inline]
    fn offsets(&self) -> (usize, usize, usize, usize) {
        (self.q_start, self.q_end, self.s_start, self.s_end)
    }
    #[inline]
    fn query_length(&self) -> usize {
        self.query_length
    }
    #[inline]
    fn move_to_query(&mut self, q_idx: u32, offset: i32, query_length: usize) {
        self.query_context = q_idx as i32;
        self.q_idx = q_idx;
        self.query_length = query_length;
        self.q_start += offset as usize;
        self.q_end += offset as usize;
        self.gapped_q_start += offset;
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(crate) struct BlastpHspList<H = BlastpHsp> {
    pub oid: u32,
    pub query_index: u32,
    pub hsps: Vec<H>,
    pub best_evalue: f64,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:168-180
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(crate) struct BlastpHitList<H = BlastpHsp> {
    pub hsplist_count: usize,
    pub hsplist_max: usize,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplist_array: Vec<BlastpHspList<H>>,
    pub hsplist_current: usize,
    pub num_hits: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:67-110
#[derive(Debug, Clone)]
pub(crate) struct BlastCompoHeapRecord {
    pub best_evalue: f64,
    pub best_score: i32,
    pub subject_index: u32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:110
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

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/compo_heap.h:79-100
#[derive(Debug, Clone)]
pub(crate) struct BlastCompoHeap {
    pub n: usize,
    pub heap_threshold: usize,
    pub ecutoff: f64,
    pub worst_evalue: f64,
    pub array: Option<Vec<BlastCompoHeapRecord>>,
    pub heap_array: Option<Vec<BlastCompoHeapRecord>>,
    /// No NCBI counterpart: `false` for a heap of a redo that ends no match early
    /// ([`BlastCompoHeap::without_early_end`]); `filled_to_cutoff` and `may_fill` are then `false`.
    pub x_ends_early: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:50-57
// (`HEAP_RESIZE_FACTOR` 1.5 and `HEAP_MIN_RESIZE` 100 size the C array's growth; the `Vec` grows
// by itself.)
const COMPO_HEAP_INITIAL_CAPACITY: usize = 100;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:63-64
pub(crate) const COMPO_HEAP_EVALUE_STRETCH: f64 = 5.0;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1389-1403
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:81-97
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:154-186
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:197-214
fn compo_heapify_up(heap_array: &mut [BlastCompoHeapRecord], mut i: usize) {
    let mut parent = i / 2;
    while parent >= 1 && compo_heap_record_compare(&heap_array[i], &heap_array[parent]) {
        heap_array.swap(i, parent);
        i = parent;
        parent /= 2;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1355
#[cfg(test)]
pub(crate) fn score_compare_hits(a: &Hit, b: &Hit) -> Ordering {
    let (a_q_offset, a_q_end) = (a.q_start.saturating_sub(1), a.q_end);
    let (b_q_offset, b_q_end) = (b.q_start.saturating_sub(1), b.q_end);
    let (a_s_offset, a_s_end) = subject_offset_end(a.s_start, a.s_end);
    let (b_s_offset, b_s_end) = subject_offset_end(b.s_start, b.s_end);

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

/// Zero-based subject offset and end of an HSP from its one-based start and end. A protein
/// subject has no reverse strand, so the start is at most one past the end; one past it is an
/// HSP of length 0, whose offset equals its end, as NCBI's HSP holds it.
#[inline(always)]
pub(crate) fn subject_offset_end(s_start: usize, s_end: usize) -> (usize, usize) {
    debug_assert!(
        s_start <= s_end + 1,
        "a protein HSP's subject start is at most one past its end"
    );
    (s_start.saturating_sub(1), s_end)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
pub(crate) fn score_compare_hsps<H: HitListHsp>(a: &H, b: &H) -> Ordering {
    let (a_q_start, a_q_end, a_s_start, a_s_end) = a.offsets();
    let (b_q_start, b_q_end, b_s_start, b_s_end) = b.offsets();
    let a_q_offset = a_q_start.saturating_sub(1);
    let b_q_offset = b_q_start.saturating_sub(1);
    let (a_s_offset, a_s_end) = subject_offset_end(a_s_start, a_s_end);
    let (b_s_offset, b_s_end) = subject_offset_end(b_s_start, b_s_end);

    match b.raw_score().cmp(&a.raw_score()) {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1414-1435
pub(crate) fn evalue_compare_hsps<H: HitListHsp>(a: &H, b: &H) -> Ordering {
    let cmp = evalue_comp(a.e_value(), b.e_value());
    if cmp == Ordering::Equal {
        score_compare_hsps(a, b)
    } else {
        cmp
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1355-1381
pub(crate) fn sort_hsplist_by_score<H: HitListHsp>(list: &mut BlastpHspList<H>) {
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1382
        // The pinned NCBI BLAST+ runs on glibc, whose qsort is a stable merge sort.
        // ScoreCompareHSPs ties HSPs with the same score and ends (other gapped starts
        // or edit scripts); they keep their order, as the containment test and the
        // endpoint purge keep the first.
        list.hsps.sort_by(score_compare_hsps);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1437-1455
// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1760-1811 (cited at
// `core::link_hsps`)
/// An HSP that takes its linked set's values from the linking (`BLAST_LinkHsps`).
pub(crate) trait LinkedHspValues: HitListHsp {
    fn take_linked(&mut self, linked: &crate::core::link_hsps::LinkedHsp);
}

impl LinkedHspValues for BlastpHsp {
    #[inline]
    fn take_linked(&mut self, linked: &crate::core::link_hsps::LinkedHsp) {
        self.e_value = linked.evalue;
        self.num = linked.num;
    }
}

impl LinkedHspValues for crate::algorithm::blastp::x_prelim_compact::XPrelimHsp {
    /// The preliminary HSP keeps no set size: every later linking starts from 1
    /// (`BLAST_LinkHsps`, link_hsps.c:1774-1776), and nothing reads it before.
    #[inline]
    fn take_linked(&mut self, linked: &crate::core::link_hsps::LinkedHsp) {
        self.e_value = linked.evalue;
    }
}

/// `BLAST_LinkHsps` on a subject's list of one query (psiblast's `-sum_stats true`): every HSP takes
/// its linked set's e-value and size (`LinkedHspValues`), in the list the linking builds, which is
/// then sorted by score, and the best e-value is updated. `context` is the query's
/// (`core::link_hsps::LinkContext`); `subject_length` the subject's length in letters.
pub(crate) fn link_hsp_list<H: LinkedHspValues>(
    list: &mut BlastpHspList<H>,
    context: crate::core::link_hsps::LinkContext,
    subject_length: i32,
    params: &crate::core::link_hsps::EvenGapLinkParams,
) {
    if list.hsps.is_empty() {
        return;
    }
    crate::core::link_hsps::link_in_place(
        &mut list.hsps,
        // `BlastHSP`'s 0-based offsets and exclusive ends (`HitListHsp::offsets`: 1-based starts).
        |hsp| {
            let (q_start, q_end, s_start, s_end) = hsp.offsets();
            crate::core::link_hsps::LinkHsp {
                context: 0,
                query_offset: q_start as i32 - 1,
                query_end: q_end as i32,
                subject_offset: subject_offset_end(s_start, s_end).0 as i32,
                subject_end: subject_offset_end(s_start, s_end).1 as i32,
                score: hsp.raw_score(),
            }
        },
        |_| context,
        subject_length,
        params,
        |hsp, linked| hsp.take_linked(linked),
    );
    sort_hsplist_by_score(list);
    update_best_evalue(list);
}

pub(crate) fn sort_hsplist_by_evalue<H: HitListHsp>(list: &mut BlastpHspList<H>) {
    if list.hsps.len() > 1 {
        let mut index = 0usize;
        while index < list.hsps.len() - 1 {
            if evalue_compare_hsps(&list.hsps[index], &list.hsps[index + 1]) == Ordering::Greater {
                break;
            }
            index += 1;
        }
        if index < list.hsps.len() - 1 {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1453-1454
            // Stable, as glibc's qsort under the pinned NCBI BLAST+ (ties of
            // ScoreCompareHSPs keep their order).
            list.hsps.sort_by(evalue_compare_hsps);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1578
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:672-674
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1485-1490
// BlastHitList2SeqAlign_OMF, which makes the Seq-aligns that the reports print (`-subject` is
// a db-scan CLocalDbAdapter, blast_app_util.cpp:207-210, so local_blast.cpp:287-290 keeps the
// database-search conversion). Each HSP list is sorted by E-value first. The redo sorted it by
// the scaled scores (blast_kappa.c:355-357) and s_HSPListNormalizeScores rounded them only
// afterwards (blast_kappa.c:108-114), so two HSPs whose scaled scores round to the same score
// and whose E-values are both below 1.0e-180 (s_EvalueComp, blast_hits.c:1390-1395) are put in
// ScoreCompareHSPs order here, the smaller subject offset first. An HSP of score 0 (kept by a
// very large -evalue) is in no report, and a subject with only such HSPs has none (blastp's hit
// list was already cut to -max_target_seqs; psiblast's is cut afterwards, `psiblast_hits_found`).
pub(crate) fn hit_list_to_seqalign_order(hit_list: &mut BlastpHitList) {
    let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
    let mut kept: Vec<_> = hit_list
        .hsplist_array
        .drain(..count)
        .filter_map(|mut hsp_list| {
            sort_hsplist_by_evalue(&mut hsp_list);
            hsp_list.hsps.retain(|hsp| hsp.raw_score != 0);
            (!hsp_list.hsps.is_empty()).then_some(hsp_list)
        })
        .collect();
    hit_list.hsplist_count = kept.len();
    kept.append(&mut hit_list.hsplist_array);
    hit_list.hsplist_array = kept;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1739-1748
pub(crate) fn update_best_evalue<H: HitListHsp>(list: &mut BlastpHspList<H>) {
    let mut best = i32::MAX as f64;
    for hsp in &list.hsps {
        if hsp.e_value() < best {
            best = hsp.e_value();
        }
    }
    list.best_evalue = best;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-1977,2049-2067
pub(crate) fn reap_hsplist_by_evalue<H: HitListHsp>(
    list: &mut BlastpHspList<H>,
    expect_value: f64,
) {
    let mut write_index = 0usize;
    for read_index in 0..list.hsps.len() {
        if list.hsps[read_index].e_value() > expect_value {
            continue;
        }
        if read_index > write_index {
            list.hsps.swap(write_index, read_index);
        }
        write_index += 1;
    }
    list.hsps.truncate(write_index);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2049-2067
// The filter of every program (`core::blast_hits::trim_by_max_hsps`).
pub(crate) fn trim_by_max_hsps(list: &mut BlastpHspList, max_hsps_per_subject: usize) {
    crate::core::blast_hits::trim_by_max_hsps(&mut list.hsps, max_hsps_per_subject);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3077-3106
fn compare_hsp_lists<H: HitListHsp>(a: &BlastpHspList<H>, b: &BlastpHspList<H>) -> Ordering {
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
    let a_score = a.hsps.first().map(|h| h.raw_score()).unwrap_or(0);
    let b_score = b.hsps.first().map(|h| h.raw_score()).unwrap_or(0);
    match b_score.cmp(&a_score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    b.oid.cmp(&a.oid)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1627-1676
fn heapify_hsplist_array<H: HitListHsp>(lists: &mut [BlastpHspList<H>], start: usize, end: usize) {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1659-1676
fn create_hsplist_heap<H: HitListHsp>(lists: &mut [BlastpHspList<H>]) {
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:414-424
    pub(crate) fn new(heap_threshold: usize, ecutoff: f64) -> Self {
        let capacity = COMPO_HEAP_INITIAL_CAPACITY.min(heap_threshold);
        let mut array = Vec::with_capacity(capacity.saturating_add(1));
        array.push(BlastCompoHeapRecord::default());
        Self {
            n: 0,
            heap_threshold,
            ecutoff,
            worst_evalue: 0.0,
            array: Some(array),
            heap_array: None,
            x_ends_early: true,
        }
    }

    /// The heap of a query whose redo ends no match early: the early termination test
    /// (`blast_compo_early_termination`) reads it as not filled to its cutoff, whatever it holds.
    ///
    /// For BLASTP with unified P-values (`ResolvedBlastpArgs::redone_match_heap`). NCBI's test can
    /// end only matches of the last query of a query batch (the heaps of the queries after the
    /// one in work are still empty), and that query's E-values are the ones NCBI adjusts with
    /// values read past its array. LOSAT's result for it is NCBI's report of the query when a
    /// copy of it follows in the batch (PD-LOSAT-NCBI-DEFECTS, the exception for unified
    /// P-values): there the query is not the last, and NCBI redoes every match of it.
    pub(crate) fn without_early_end(mut self) -> Self {
        self.x_ends_early = false;
        self
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:221-238
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-255,257-284
    pub(crate) fn would_insert(&mut self, evalue: f64, score: i32, subject_index: u32) -> bool {
        if self.would_insert_before_root(evalue) {
            return true;
        }
        if self.heap_array.is_none() {
            self.convert_to_heap();
        }
        self.would_replace_root(evalue, score, subject_index)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-255
    /// The clauses of `would_insert` that need no root: the heap holds fewer lists than its
    /// threshold, or the E-value is within the cutoff, or below the worst E-value kept.
    #[inline(always)]
    fn would_insert_before_root(&self, evalue: f64) -> bool {
        self.n < self.heap_threshold || evalue <= self.ecutoff || evalue < self.worst_evalue
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:257-284
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:230-247
    /// The last clause of `would_insert`: the root (`compo_heap_record_compare`'s worst record) is
    /// worse than the candidate. In the heap form the root is `heap_array[1]`; in the array form it
    /// is the record of `array[1..=n]` that the comparison puts above every other one, which is the
    /// root that the conversion (`convert_to_heap`) makes: the comparison reads only the E-value,
    /// the score and the subject index, and orders distinct keys strictly, so every maximum has the
    /// root's key.
    #[inline(always)]
    fn would_replace_root(&self, evalue: f64, score: i32, subject_index: u32) -> bool {
        let record = BlastCompoHeapRecord {
            best_evalue: evalue,
            best_score: score,
            subject_index,
            hsplist: None,
        };
        let root = match (&self.heap_array, &self.array) {
            (Some(heap_array), _) => &heap_array[1],
            (None, Some(array)) => array[2..=self.n].iter().fold(&array[1], |root, other| {
                if compo_heap_record_compare(other, root) {
                    other
                } else {
                    root
                }
            }),
            (None, None) => unreachable!("a compo heap holds its array or its heap array"),
        };
        compo_heap_record_compare(root, &record)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:251-276
    /// The answer of `would_insert` on the heap as it is, without its one side effect: an array
    /// that has reached the threshold is not turned into its heap form (`convert_to_heap`); the
    /// same clauses (`would_insert_before_root`, `would_replace_root`, which reads the root of the
    /// array form). The windows of the parallel redo ask it from many threads on the heaps as a
    /// window's start left them (`x_redo_par.rs`).
    // Only the parallel redo's windows ask it (the module `x_redo_par` is compiled where they run).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(crate) fn would_insert_view(&self, evalue: f64, score: i32, subject_index: u32) -> bool {
        self.would_insert_before_root(evalue)
            || self.would_replace_root(evalue, score, subject_index)
    }

    // No NCBI counterpart: a copy for looking ahead; NCBI's heap has none.
    /// A copy of the heap whose records keep their keys (E-value, score, subject index) but not
    /// their alignments. `would_insert` and `insert` read and order only the keys, so the copy
    /// answers and changes as the heap would under the same calls. The parallel redo's windows
    /// replay a query's carried matches ahead on it (`x_redo_par.rs`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(crate) fn keys(&self) -> Self {
        let keys = |records: &Vec<BlastCompoHeapRecord>| {
            records
                .iter()
                .map(|record| BlastCompoHeapRecord {
                    best_evalue: record.best_evalue,
                    best_score: record.best_score,
                    subject_index: record.subject_index,
                    hsplist: None,
                })
                .collect()
        };
        Self {
            n: self.n,
            heap_threshold: self.heap_threshold,
            ecutoff: self.ecutoff,
            worst_evalue: self.worst_evalue,
            array: self.array.as_ref().map(keys),
            heap_array: self.heap_array.as_ref().map(keys),
            x_ends_early: self.x_ends_early,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:294-352
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:405-408
    pub(crate) fn filled_to_cutoff(&self) -> bool {
        self.x_ends_early && self.n >= self.heap_threshold && self.worst_evalue <= self.ecutoff
    }

    /// Whether `more` insertions can fill the heap to its cutoff (`filled_to_cutoff`): the bound
    /// that the window redo's look ahead uses (`x_redo_par::x_end_above`). No NCBI counterpart.
    #[cfg_attr(
        all(
            not(test),
            not(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))
        ),
        expect(
            dead_code,
            reason = "used only with `parallel` threads (x_redo_par) and by unit tests"
        )
    )]
    pub(crate) fn may_fill(&self, more: usize) -> bool {
        self.x_ends_early && self.n.saturating_add(more) >= self.heap_threshold
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:443-444
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

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1560-1576
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2492-2514
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3420-3438
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

impl<H: HitListHsp> BlastpHitList<H> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3125-3133
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3219-3240
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3196-3209
    fn insert_hsplist_in_heap(&mut self, hsp_list: BlastpHspList<H>) {
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
            if let Some(score) = root.hsps.first().map(|h| h.raw_score()) {
                self.low_score = score;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3297
    pub(crate) fn update(&mut self, mut hsp_list: BlastpHspList<H>) {
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
                .map(|h| h.raw_score())
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3172-3187
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3331-3337
    pub(crate) fn sort_by_evalue(&mut self) {
        if self.hsplist_count > 1 {
            self.hsplist_array[..self.hsplist_count].sort_unstable_by(compare_hsp_lists);
        }
        self.purge();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:827-857
    // (`core::blast_hspstream`)
    /// `BlastHSPCBSStreamClose`'s trim of this query's preliminary hit list (`ADAPTIVE_CBS`;
    /// BLASTP's, PSIBLAST's and RPSBLAST's lists): sorted by E-value and cut where the E-values
    /// pass the limit, when it holds enough lists.
    pub(crate) fn cbs_close_trim(&mut self, hitlist_size: usize) {
        use crate::core::blast_hspstream::{cbs_trim_count, cbs_trim_sorts};
        let count = self.hsplist_count;
        if !cbs_trim_sorts(count, hitlist_size) {
            return;
        }
        self.sort_by_evalue();
        let lists = &self.hsplist_array;
        if let Some(keep) = cbs_trim_count(count, self.hsplist_count, hitlist_size, |i| {
            lists[i].best_evalue
        }) {
            self.hsplist_array.truncate(keep);
            self.hsplist_count = keep;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3314-3326
    /// The lists that a filter left without HSPs (NCBI frees them) removed, the others kept
    /// in their order.
    pub(crate) fn purge_empty_hsp_lists(&mut self) {
        self.hsplist_array.truncate(self.hsplist_count);
        self.hsplist_array.retain(|list| !list.hsps.is_empty());
        self.hsplist_count = self.hsplist_array.len();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:877-892
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2454-2535
#[cfg(test)]
pub(crate) fn purge_hits_with_common_endpoints(hits: Vec<Hit>) -> Vec<Hit> {
    let (hits, _) = purge_hits_with_common_endpoints_ex(hits, true);
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2537
#[allow(dead_code)]
pub(crate) fn purge_hits_with_common_endpoints_ex(
    hits: Vec<Hit>,
    purge: bool,
) -> (Vec<Hit>, usize) {
    purge_hits_for_subject_ex(hits, purge)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2391-2452
fn cut_off_gap_edit_script(hit: &mut Hit, q_cut: usize, s_cut: usize, cut_begin: bool) -> bool {
    let gap_info = match &hit.gap_info {
        Some(info) if !info.is_empty() => info.clone(),
        _ => return false,
    };

    let (mut q_offset_0, mut q_end_0) = (hit.q_start.saturating_sub(1), hit.q_end);
    let (mut s_offset_0, mut s_end_0) = subject_offset_end(hit.s_start, hit.s_end);

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2267-2315
fn hit_query_offset_compare(a: &Hit, b: &Hit) -> Ordering {
    (a.q_idx, a.query_frame)
        .cmp(&(b.q_idx, b.query_frame))
        .then_with(|| {
            a.q_start
                .saturating_sub(1)
                .cmp(&b.q_start.saturating_sub(1))
        })
        .then_with(|| {
            subject_offset_end(a.s_start, a.s_end)
                .0
                .cmp(&subject_offset_end(b.s_start, b.s_end).0)
        })
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.q_end.cmp(&a.q_end))
        .then_with(|| {
            subject_offset_end(b.s_start, b.s_end)
                .1
                .cmp(&subject_offset_end(a.s_start, a.s_end).1)
        })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2332-2379
fn hit_query_end_compare(a: &Hit, b: &Hit) -> Ordering {
    (a.q_idx, a.query_frame)
        .cmp(&(b.q_idx, b.query_frame))
        .then_with(|| a.q_end.cmp(&b.q_end))
        .then_with(|| {
            subject_offset_end(a.s_start, a.s_end)
                .1
                .cmp(&subject_offset_end(b.s_start, b.s_end).1)
        })
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| {
            b.q_start
                .saturating_sub(1)
                .cmp(&a.q_start.saturating_sub(1))
        })
        .then_with(|| {
            subject_offset_end(b.s_start, b.s_end)
                .0
                .cmp(&subject_offset_end(a.s_start, a.s_end).0)
        })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2500
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2514-2525
const ENDPOINT_REMOVAL_BITMAP_THRESHOLD_DIVISOR: usize = 8;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2500
//
// NCBI mutates the pointer array in place without allocating a side bitmap.
// Rust keeps the same active-prefix result, but only pays for a bitmap on large
// HSP lists where repeated removed-index checks would dominate.
const ENDPOINT_REMOVAL_BITMAP_MIN_HITS: usize = 256;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2500
//
// NCBI's pointer array removes common-endpoint HSPs from the active prefix and
// leaves freed or trimmed entries outside that prefix. Rust records the same
// removed active-prefix indices lazily, switching to a bitmap only when the
// duplicate set is dense enough that index scanning becomes wasteful.
enum EndpointRemovalMarks {
    Indices(Vec<usize>),
    Bitmap(Vec<bool>),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2500
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct EndpointPurgePassStats {
    removed_count: usize,
    trimmed_count: usize,
    allocated_removal_marks: bool,
    used_bitmap: bool,
    compacted: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2527
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct EndpointPurgeStats {
    start: EndpointPurgePassStats,
    end: EndpointPurgePassStats,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2500
fn mark_endpoint_removed<const TRACK_STATS: bool>(
    marks: &mut Option<EndpointRemovalMarks>,
    index: usize,
    hits_len: usize,
    stats: &mut EndpointPurgePassStats,
) {
    if marks.is_none() {
        if TRACK_STATS {
            stats.allocated_removal_marks = true;
        }
        *marks = Some(EndpointRemovalMarks::Indices(Vec::with_capacity(4)));
    }

    let dense_threshold = if hits_len >= ENDPOINT_REMOVAL_BITMAP_MIN_HITS {
        hits_len / ENDPOINT_REMOVAL_BITMAP_THRESHOLD_DIVISOR
    } else {
        usize::MAX
    };
    match marks
        .as_mut()
        .expect("removal marks are allocated before recording an index")
    {
        EndpointRemovalMarks::Indices(indices)
            if indices.len().saturating_add(1) > dense_threshold =>
        {
            let mut bitmap = vec![false; hits_len];
            for removed_index in indices.drain(..) {
                bitmap[removed_index] = true;
            }
            bitmap[index] = true;
            *marks = Some(EndpointRemovalMarks::Bitmap(bitmap));
            if TRACK_STATS {
                stats.used_bitmap = true;
            }
        }
        EndpointRemovalMarks::Indices(indices) => {
            indices.push(index);
        }
        EndpointRemovalMarks::Bitmap(bitmap) => {
            bitmap[index] = true;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2496-2500
//
// Rust performs the same active-prefix compaction once after the endpoint pass.
// The active HSP order is preserved; removed HSPs are dropped in one retain pass
// without building a replacement active Vec.
fn compact_hits_after_endpoint_purge(mut hits: Vec<Hit>, marks: EndpointRemovalMarks) -> Vec<Hit> {
    match marks {
        EndpointRemovalMarks::Indices(indices) => {
            let mut removed_iter = indices.into_iter();
            let mut next_removed = removed_iter.next();
            let mut read_index = 0usize;
            hits.retain(|_| {
                let keep = next_removed != Some(read_index);
                if !keep {
                    next_removed = removed_iter.next();
                }
                read_index += 1;
                keep
            });
        }
        EndpointRemovalMarks::Bitmap(bitmap) => {
            let mut read_index = 0usize;
            hits.retain(|_| {
                let keep = !bitmap[read_index];
                read_index += 1;
                keep
            });
        }
    }
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2502
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2504-2527
fn purge_sorted_hits_with_common_endpoints<const TRACK_STATS: bool, S, T>(
    hits: Vec<Hit>,
    purge: bool,
    same_endpoint: S,
    trim_if_kept: T,
) -> (Vec<Hit>, Vec<Hit>, EndpointPurgePassStats)
where
    S: Fn(&Hit, &Hit) -> bool,
    T: Fn(&Hit, &Hit, &mut Hit) -> bool,
{
    let hits_len = hits.len();
    let mut removal_marks: Option<EndpointRemovalMarks> = None;
    let mut trimmed_hits = Vec::new();
    let mut stats = EndpointPurgePassStats::default();

    let mut i = 0usize;
    while i < hits_len {
        let mut j = i + 1;
        while j < hits_len && same_endpoint(&hits[i], &hits[j]) {
            mark_endpoint_removed::<TRACK_STATS>(&mut removal_marks, j, hits_len, &mut stats);
            if TRACK_STATS {
                stats.removed_count += 1;
            }
            if !purge {
                let mut trimmed_hit = hits[j].clone();
                if trim_if_kept(&hits[i], &hits[j], &mut trimmed_hit) {
                    trimmed_hits.push(trimmed_hit);
                    if TRACK_STATS {
                        stats.trimmed_count += 1;
                    }
                }
            }
            j += 1;
        }
        i = j;
    }

    let hits = if let Some(removal_marks) = removal_marks {
        if TRACK_STATS {
            stats.compacted = true;
        }
        compact_hits_after_endpoint_purge(hits, removal_marks)
    } else {
        hits
    };

    (hits, trimmed_hits, stats)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535
fn purge_hits_for_subject_ex(hits: Vec<Hit>, purge: bool) -> (Vec<Hit>, usize) {
    let (hits, extra_start, _) = purge_hits_for_subject_ex_impl::<false>(hits, purge);
    (hits, extra_start)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535
#[cfg(test)]
fn purge_hits_for_subject_ex_with_stats(
    hits: Vec<Hit>,
    purge: bool,
) -> (Vec<Hit>, usize, EndpointPurgeStats) {
    purge_hits_for_subject_ex_impl::<true>(hits, purge)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535
fn purge_hits_for_subject_ex_impl<const TRACK_STATS: bool>(
    mut hits: Vec<Hit>,
    purge: bool,
) -> (Vec<Hit>, usize, EndpointPurgeStats) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len, EndpointPurgeStats::default());
    }

    let context = |h: &Hit| (h.q_idx, h.query_frame);
    let q_offset = |h: &Hit| h.q_start.saturating_sub(1);
    let q_end = |h: &Hit| h.q_end;
    let s_offset = |h: &Hit| subject_offset_end(h.s_start, h.s_end).0;
    let s_end = |h: &Hit| subject_offset_end(h.s_start, h.s_end).1;

    let mut trimmed_hits = Vec::new();
    let mut purged_count = 0usize;
    let mut stats = EndpointPurgeStats::default();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478
    // Stable, as glibc's qsort under the pinned NCBI BLAST+: HSPs with the same
    // context, ends and score tie, and the purge keeps the first.
    hits.sort_by(hit_query_offset_compare);
    let (active_hits, mut start_trimmed_hits, start_stats) =
        purge_sorted_hits_with_common_endpoints::<TRACK_STATS, _, _>(
            hits,
            purge,
            |best, candidate| {
                context(best) == context(candidate)
                    && q_offset(best) == q_offset(candidate)
                    && s_offset(best) == s_offset(candidate)
            },
            |best, candidate, trimmed_hit| {
                if q_end(candidate) > q_end(best) {
                    let _ = cut_off_gap_edit_script(trimmed_hit, q_end(best), s_end(best), true);
                    true
                } else {
                    false
                }
            },
        );
    hits = active_hits;
    if TRACK_STATS {
        purged_count += start_stats.removed_count;
    }
    trimmed_hits.append(&mut start_trimmed_hits);
    if TRACK_STATS {
        stats.start = start_stats;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2504
    // Stable, as glibc's qsort under the pinned NCBI BLAST+ (as above).
    hits.sort_by(hit_query_end_compare);
    let (active_hits, mut end_trimmed_hits, end_stats) =
        purge_sorted_hits_with_common_endpoints::<TRACK_STATS, _, _>(
            hits,
            purge,
            |best, candidate| {
                context(best) == context(candidate)
                    && q_end(best) == q_end(candidate)
                    && s_end(best) == s_end(candidate)
            },
            |best, candidate, trimmed_hit| {
                if q_offset(candidate) < q_offset(best) {
                    let _ =
                        cut_off_gap_edit_script(trimmed_hit, q_offset(best), s_offset(best), false);
                    true
                } else {
                    false
                }
            },
        );
    hits = active_hits;
    if TRACK_STATS {
        purged_count += end_stats.removed_count;
    }
    trimmed_hits.append(&mut end_trimmed_hits);
    if TRACK_STATS {
        stats.end = end_stats;
    }

    let extra_start = hits.len();
    hits.extend(trimmed_hits);
    debug_assert!(!TRACK_STATS || extra_start + purged_count == len);
    (hits, extra_start, stats)
}

/// `collect_hits_from_hit_lists` with each HSP's `comp_adjustment_method`, which the
/// outfmt 0 report prints (showalign.cpp:3600-3603), and its Seq-align's `num_ident`
/// (`PairwiseHit::seqalign_num_ident`): the engine's identities that `engine_counts` keeps
/// where the display recount replaced them, else the HSP's. `first_query` is the index of the
/// query of `hit_lists[0]` (the keys of `engine_counts` count the queries from the run's first).
pub(crate) fn collect_hits_and_methods_from_hit_lists(
    hit_lists: &[Option<BlastpHitList>],
    first_query: usize,
    engine_counts: &super::seqalign::EngineCounts,
) -> Vec<(Hit, u8, usize, i32)> {
    let reported = |hit_list: &BlastpHitList| {
        hit_list
            .hsplist_array
            .iter()
            .take(hit_list.hsplist_count)
            .map(|hsp_list| hsp_list.hsps.len())
            .sum::<usize>()
    };
    let mut hits = Vec::with_capacity(hit_lists.iter().flatten().map(reported).sum());
    for (q_offset, hit_list_opt) in hit_lists.iter().enumerate() {
        let Some(hit_list) = hit_list_opt else {
            continue;
        };
        let q_idx = (first_query + q_offset) as u32;
        for (list_index, hsp_list) in hit_list
            .hsplist_array
            .iter()
            .take(hit_list.hsplist_count)
            .enumerate()
        {
            for (hsp_index, hsp) in hsp_list.hsps.iter().enumerate() {
                let seqalign_num_ident = engine_counts
                    .num_ident((q_idx, list_index as u32, hsp_index as u32), hsp.num_ident);
                hits.push((
                    hsp.clone().into_hit(),
                    hsp.comp_adjustment_method,
                    seqalign_num_ident,
                    hsp.num,
                ));
            }
        }
    }
    hits
}

#[cfg(test)]
#[path = "hsp_tests.rs"]
mod tests;
