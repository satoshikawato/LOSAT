//! EXPERIMENT (LOSAT_X_TBXKMERGE / LOSAT_X_TBXKMERGESHADOW): the merge of a subject's unit lists
//! (`XPairPar::merge_subject`) without its full sorts. A child module of x_pair_par.rs.
//!
//! With the switch off, `merge_subject_listed` builds the subject's `combined_ungapped_hits` as the
//! reference frame loop does: for each subject frame and chunk it appends the lists of the chunk's
//! units in unit order and sorts the chunk's init list (`Blast_InitHitListSortByScore`), converts
//! it (`BLAST_GetUngappedHSPList`, which ends with `Blast_HSPListSortByScore`), adjusts the offsets
//! and merges the chunk list into the frame list; then it appends the frame list to the combined
//! list and sorts that (`Blast_HSPListAppend`). Each of these sorts is a full stable sort of a list
//! made of pieces that are already sorted, or nearly so.
//!
//! What the switch changes (read once per process, `mode`; in the strict set):
//! 1. Each unit sorts its own list, in its task on the pool, after its scan (`run_units`), by the
//!    comparator of the chunk's init sort (`XRecord::cmp`); a compact list with the task's scratch
//!    vector as the buffer (`x_stable_sort_with`, no allocation; `map_init` gives each split of the
//!    pool's work its own scratch).
//! 2. A chunk's init list is the k-way merge of its units' sorted lists (`x_kway`) in a Vec of their
//!    total length: the smallest head by the comparator, ties to the lower unit. No sort of the
//!    chunk list.
//! 3. The conversion runs the ported loop (`convert_init_hsps_unsorted`); its final sort is a stable
//!    sort of each maximal run of HSPs equal in K1 = (score desc, subject start asc, subject end
//!    desc), the first three keys of `ScoreCompareHSPs` (`x_sort_k1_runs`). A list that is not in
//!    K1 order (never expected) gets the final sort of the ported conversion instead.
//! 4. Each frame list is merged into the combined list by a backward two-way merge (`x_merge_frame`),
//!    ties to the combined list, instead of being appended and the whole list sorted. A frame list
//!    that is not sorted (never expected) is sorted first.
//! The merges of the chunks of one frame (`merge_tblastx_subject_chunk_hits`, subjects with more
//! than 5 M residues in a frame) are unchanged from the second chunk on; the first chunk's list
//! becomes the empty frame list, as `merge_subject_compact` does in both forms, which is what the
//! empty-list branch of `merge_tblastx_subject_chunk_hits` (`combined.append`) computes.
//!
//! Why `combined_ungapped_hits` is the list of the switch off:
//! a. Every sort involved is a stable sort by a lexicographic comparator over integers, so its
//!    result is unique: the elements in comparator order, equal elements in input order. Hence, for
//!    lists A and B, stable_sort(A ++ B) is the merge of stable_sort(A) and stable_sort(B) that takes
//!    from A on ties (the equal elements of A come before those of B in A ++ B, in A's order), and
//!    for k lists the merge with ties to the lower list.
//! b. Items 1 and 2: the chunk's init list of the switch off is stable_sort(U_1 ++ ... ++ U_k) of
//!    its units' lists in unit order; the units' sorts are stable_sort(U_i) (`XRecord::cmp` on a
//!    compact record is `score_compare_init_hsps_ncbi` on its `InitHSP`: the record's `len` is
//!    `s_end - s_start` and its `q_start` is `q_start_absolute`), and `x_kway` is the merge of (a).
//!    The records are expanded by the statements of `x_expand_unit` (`XExpand`).
//! c. Item 3: the conversion keeps the init order and assigns `hsp_list_order` in it. The init list
//!    is sorted by (score, s_start, length), so the converted list is sorted by K1: before the
//!    offset adjustment `raw_score` is the score, `s_aa_start` is `s_start` and `s_aa_end` is
//!    `s_start + length`. A stable sort of a K1-sorted list by a comparator that begins with K1
//!    moves no element across the boundary of a K1 run and orders each run stably by the rest of
//!    the comparator: the run-wise sort. Without K1 order the ported sort runs, so the result does
//!    not rest on this argument.
//! d. Item 4: the combined list of the switch off after frame f is stable_sort(C ++ F) where C,
//!    the list after the frames before f, is sorted (from the empty list on) and F, the frame list,
//!    is sorted (it leaves a sort: the conversion's, or that of `merge_tblastx_subject_chunk_hits`;
//!    `adjust_tblastx_chunk_subject_offsets` between them adds one offset to the subject start and
//!    end of every HSP, which keeps the order); the backward merge with ties to C is the merge of
//!    (a).
//!
//! Memory: the units' lists of a chunk live until its merge ends (with the switch off each is
//! dropped once it is expanded); the chunk's init list is made at its length; the frame merge holds
//! the combined list grown by the frame list's length and the frame list, as the append does, and
//! no index array or sort buffer.
//!
//! `LOSAT_X_TBXKMERGESHADOW=1` (opt-in, whatever the switch): the units return unsorted lists; the
//! merge copies them, runs `merge_subject_listed` on the copy and this file's merge on the lists
//! (sorting them first; a compact list both with `x_stable_sort_with` and with `slice::sort_by`,
//! which must agree), and compares every field of every HSP; the first difference stops the run.
//! With `LOSAT_X_STATS` the counters of the items above are printed (`print_summary`) whenever the
//! switch or its shadow is on, with zeros where no unit path ran.

use super::*;
use crate::algorithm::tblastx::blast_gapalign::{
    convert_init_hsps_unsorted, score_compare_init_hsps_ncbi, score_compare_ungapped_hits_ncbi,
    sort_ungapped_hits_by_score_ncbi, ungapped_hits_is_sorted_by_score_ncbi,
};
use std::cmp::Ordering;

// No NCBI counterpart: counters of LOSAT_X_TBXKMERGE (LOSAT_X_STATS) and of its shadow; they do not
// change any value NCBI computes.
static KM_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static KM_CHUNKS: AtomicU64 = AtomicU64::new(0);
static KM_TIE_RUNS: AtomicU64 = AtomicU64::new(0);
static KM_K1_FALLBACK: AtomicU64 = AtomicU64::new(0);
static KM_FRAME_UNSORTED: AtomicU64 = AtomicU64::new(0);
static KM_SHADOW_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static KM_SHADOW_HSPS: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: chooses how LOSAT orders the same lists; it does not change any value NCBI
// computes. Read once per process, as the other modes of x_pair_par.rs.
/// 0 = off, 1 = LOSAT_X_TBXKMERGE (a strict-set switch: on unless `LOSAT_X_TBXKMERGE` is `0`;
/// `utils/xswitch.rs`), 2 = LOSAT_X_TBXKMERGESHADOW.
pub(super) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_TBXKMERGESHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_TBXKMERGE") {
            1
        } else {
            0
        }
    })
}

/// A saved init HSP of a unit, in either form of the units' lists (`XUnitOutputs`).
pub(super) trait XRecord: Sized {
    /// The order of the chunk's init sort (`score_compare_init_hsps_ncbi` on the `InitHSP`).
    fn cmp(a: &Self, b: &Self) -> Ordering;
    /// The `InitHSP` the switch off puts into the chunk's init list for this record.
    fn init_hsp(self, unit: &XExpand<'_, '_>) -> InitHSP;
}

impl XRecord for InitHSP {
    fn cmp(a: &Self, b: &Self) -> Ordering {
        score_compare_init_hsps_ncbi(a, b)
    }

    fn init_hsp(self, _unit: &XExpand<'_, '_>) -> InitHSP {
        self
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-304
// ```c
//     if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
//                                  h1->ungapped_data->score)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
//                                  h2->ungapped_data->s_start)) &&
//         0 == (result = BLAST_CMP(h2->ungapped_data->length,
//                                  h1->ungapped_data->length)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
//                                  h2->ungapped_data->q_start))) {
// ```
// `score_compare_init_hsps_ncbi` on the record's own fields, which are NCBI's `BlastUngappedData`
// (section above `XUnitHsp` in x_pair_par.rs): `len` is `s_end - s_start` of the `InitHSP` and
// `q_start` its `q_start_absolute`.
impl XRecord for XUnitHsp {
    fn cmp(a: &Self, b: &Self) -> Ordering {
        b.score
            .cmp(&a.score)
            .then_with(|| a.s_start.cmp(&b.s_start))
            .then_with(|| b.len.cmp(&a.len))
            .then_with(|| a.q_start.cmp(&b.q_start))
    }

    fn init_hsp(self, unit: &XExpand<'_, '_>) -> InitHSP {
        unit.init_hsp(&self)
    }
}

/// A unit's list, sorted in its task (item 1).
pub(super) trait XUnitSort {
    /// The sort of the shadow mode and of the tests (`slice::sort_by`).
    fn x_sort(&mut self);
    /// The sort of the unit's task, which may use the task's scratch as its buffer.
    fn x_sort_in(&mut self, scratch: &mut XScratch);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
// ```c
// void Blast_InitHitListSortByScore(BlastInitHitList * init_hitlist)
// {
//     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//           sizeof(BlastInitHSP), score_compare_match);
// }
// ```
// The comparator of the chunk's sort on one unit's part of the chunk's list; `slice::sort_by` is
// stable, as `sort_init_hsps_by_score_ncbi` is.
impl XUnitSort for Vec<InitHSP> {
    fn x_sort(&mut self) {
        self.sort_by(<InitHSP as XRecord>::cmp);
    }

    fn x_sort_in(&mut self, _scratch: &mut XScratch) {
        self.x_sort();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
// ```c
//     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//           sizeof(BlastInitHSP), score_compare_match);
// ```
// The compact records of a unit are its list twice when its task ends: the exact copy the task
// returns and the scratch vector it was copied from (`run_unit_compact`; the scratch, one per split
// of the pool's work, is cleared at the start of its next unit). In the task the stable sort runs
// between the two (`x_stable_sort_with`), so it allocates nothing; `slice::sort_by` allocates a
// buffer of up to the list's length in each of the pool's tasks at once.
impl XUnitSort for Vec<XUnitHsp> {
    fn x_sort(&mut self) {
        self.sort_by(<XUnitHsp as XRecord>::cmp);
    }

    fn x_sort_in(&mut self, scratch: &mut XScratch) {
        x_stable_sort_with(self, &mut scratch.records, <XUnitHsp as XRecord>::cmp);
    }
}

/// Length of the runs `x_stable_sort_with` sorts by insertion before it merges.
const X_SORT_RUN: usize = 16;

/// A stable sort of `v` by `cmp` with `buf` as the merge buffer (its capacity is used, its contents
/// are overwritten): runs of `X_SORT_RUN` sorted by insertion, then merged bottom-up between `v`
/// and `buf`, taking from the left run unless the right element is smaller.
fn x_stable_sort_with<T: Copy>(v: &mut [T], buf: &mut Vec<T>, cmp: impl Fn(&T, &T) -> Ordering) {
    let n = v.len();
    if n <= 1 {
        return;
    }
    for run in v.chunks_mut(X_SORT_RUN) {
        for i in 1..run.len() {
            let x = run[i];
            let mut j = i;
            while j > 0 && cmp(&run[j - 1], &x) == Ordering::Greater {
                run[j] = run[j - 1];
                j -= 1;
            }
            run[j] = x;
        }
    }
    if n <= X_SORT_RUN {
        return;
    }
    buf.clear();
    buf.extend_from_slice(v);
    let mut in_v = true;
    let mut width = X_SORT_RUN;
    while width < n {
        if in_v {
            x_merge_pass(v, buf, width, &cmp);
        } else {
            x_merge_pass(buf, v, width, &cmp);
        }
        in_v = !in_v;
        width *= 2;
    }
    if !in_v {
        v.copy_from_slice(&buf[..n]);
    }
}

/// One bottom-up pass of `x_stable_sort_with`: merges the neighbouring runs of `width` of `src`
/// into `dst`.
fn x_merge_pass<T: Copy>(
    src: &[T],
    dst: &mut [T],
    width: usize,
    cmp: &impl Fn(&T, &T) -> Ordering,
) {
    let n = src.len();
    let mut start = 0;
    while start < n {
        let mid = (start + width).min(n);
        let end = (start + 2 * width).min(n);
        let (mut i, mut j, mut k) = (start, mid, start);
        while i < mid && j < end {
            if cmp(&src[j], &src[i]) == Ordering::Less {
                dst[k] = src[j];
                j += 1;
            } else {
                dst[k] = src[i];
                i += 1;
            }
            k += 1;
        }
        dst[k..k + (mid - i)].copy_from_slice(&src[i..mid]);
        k += mid - i;
        dst[k..k + (end - j)].copy_from_slice(&src[j..end]);
        start = end;
    }
}

/// The stable merge of `lists`, each sorted by `R::cmp`, ties to the lower list (item 2), each
/// record expanded by `unit`, in a Vec of their total length.
fn x_kway<R: XRecord>(lists: Vec<Vec<R>>, unit: &XExpand<'_, '_>) -> Vec<InitHSP> {
    let total = lists.iter().map(Vec::len).sum();
    let mut merged = Vec::with_capacity(total);
    // The empty lists are left out; the others keep their order, so a tie still goes to the lower
    // unit.
    let mut lists: Vec<std::vec::IntoIter<R>> = lists
        .into_iter()
        .filter(|list| !list.is_empty())
        .map(Vec::into_iter)
        .collect();
    let mut heads: Vec<Option<R>> = lists.iter_mut().map(Iterator::next).collect();
    let mut live = heads.len();
    while live > 1 {
        let mut best: Option<usize> = None;
        for (i, head) in heads.iter().enumerate() {
            let Some(r) = head else { continue };
            let smaller = match best {
                None => true,
                Some(b) => R::cmp(r, heads[b].as_ref().expect("a live head")) == Ordering::Less,
            };
            if smaller {
                best = Some(i);
            }
        }
        let b = best.expect("a live head");
        let next = lists[b].next();
        if next.is_none() {
            live -= 1;
        }
        let r = std::mem::replace(&mut heads[b], next).expect("a live head");
        merged.push(r.init_hsp(unit));
    }
    if let Some(b) = heads.iter().position(Option::is_some) {
        let r = heads[b].take().expect("the last live head");
        merged.push(r.init_hsp(unit));
        merged.extend(lists[b].by_ref().map(|r| r.init_hsp(unit)));
    }
    merged
}

/// K1: the first three keys of `score_compare_ungapped_hits_ncbi` (score, subject start, subject end).
fn x_k1(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.s_aa_start.cmp(&b.s_aa_start))
        .then_with(|| b.s_aa_end.cmp(&a.s_aa_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
// ```c
// void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
// {
//     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// The final sort of the conversion as the stable sort of each K1 run (item 3, argument c); a list
// out of K1 order gets the ported sort.
/// Sorts a converted chunk list as `Blast_HSPListSortByScore` does.
fn x_sort_k1_runs(hits: &mut [UngappedHit]) {
    if hits
        .windows(2)
        .any(|pair| x_k1(&pair[0], &pair[1]) == Ordering::Greater)
    {
        KM_K1_FALLBACK.fetch_add(1, AtomicOrdering::Relaxed);
        sort_ungapped_hits_by_score_ncbi(hits);
        return;
    }
    let mut runs = 0u64;
    let mut start = 0usize;
    while start < hits.len() {
        let mut end = start + 1;
        while end < hits.len() && x_k1(&hits[end - 1], &hits[end]) == Ordering::Equal {
            end += 1;
        }
        if end - start > 1 {
            hits[start..end].sort_by(score_compare_ungapped_hits_ncbi);
            runs += 1;
        }
        start = end;
    }
    KM_TIE_RUNS.fetch_add(runs, AtomicOrdering::Relaxed);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766
// ```c
// for (index=combined_hsp_list->hspcnt, index1=0;
//      index1<hsp_list->hspcnt; index1++) {
//    combined_hsp_list->hsp_array[index++] = hsp_list->hsp_array[index1];
// }
// combined_hsp_list->hspcnt = new_hspcnt;
// Blast_HSPListSortByScore(combined_hsp_list);
// ```
// The append and the sort as one backward merge of the sorted frame list into the sorted combined
// list, ties to the combined list (item 4, argument d). The frame list's elements are first copied
// onto the end of the combined list as place holders, so that every step is a swap.
/// Appends `frame` to the sorted `combined` and sorts the result as `Blast_HSPListAppend` does.
fn x_merge_frame(combined: &mut Vec<UngappedHit>, mut frame: Vec<UngappedHit>) {
    if frame.is_empty() {
        return;
    }
    if !ungapped_hits_is_sorted_by_score_ncbi(&frame) {
        KM_FRAME_UNSORTED.fetch_add(1, AtomicOrdering::Relaxed);
        sort_ungapped_hits_by_score_ncbi(&mut frame);
    }
    if combined.is_empty() {
        *combined = frame;
        return;
    }
    let (mut i, mut j) = (combined.len(), frame.len());
    combined.reserve_exact(j);
    combined.extend_from_slice(&frame);
    // combined[..i] are the combined elements not yet placed, combined[i + j..] the merged tail,
    // frame[..j] the frame elements not yet placed.
    while j > 0 {
        let k = i + j - 1;
        if i > 0
            && score_compare_ungapped_hits_ncbi(&combined[i - 1], &frame[j - 1])
                == Ordering::Greater
        {
            combined.swap(i - 1, k);
            i -= 1;
        } else {
            std::mem::swap(&mut combined[k], &mut frame[j - 1]);
            j -= 1;
        }
    }
}

/// `combined_ungapped_hits` of one subject from its units' sorted lists (items 2-4); `count` adds
/// the init HSPs to the `[X_STATS] PAIRPAR` counter.
fn x_merge<R: XRecord>(
    frame_chunks: &[Vec<SubjectChunk>],
    units: &[XUnit],
    lists: Vec<Vec<R>>,
    inp: &XSubject<'_>,
    count: bool,
) -> Vec<UngappedHit> {
    let mut lists = lists.into_iter();
    let mut next_unit = 0usize;
    let mut chunks_merged = 0u64;
    let mut combined_ungapped_hits: Vec<UngappedHit> = Vec::new();
    for (f, chunks) in frame_chunks.iter().enumerate() {
        let unit = XExpand::new(f, inp);
        let mut frame_ungapped_hits: Vec<UngappedHit> = Vec::new();
        for (k, chunk) in chunks.iter().enumerate() {
            let mut chunk_lists = Vec::new();
            while next_unit < units.len() && units[next_unit].f == f && units[next_unit].k == k {
                chunk_lists.push(lists.next().expect("one output per unit"));
                next_unit += 1;
            }
            let init_hsps = x_kway(chunk_lists, &unit);
            chunks_merged += 1;
            if count {
                INIT_HSPS.fetch_add(init_hsps.len() as u64, AtomicOrdering::Relaxed);
            }
            let mut hits = if init_hsps.is_empty() {
                Vec::new()
            } else {
                let exact = Vec::with_capacity(init_hsps.len());
                let mut hits =
                    convert_init_hsps_unsorted(init_hsps, inp.contexts, inp.s_frames, exact);
                x_sort_k1_runs(&mut hits);
                hits
            };
            adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
            if !hits.is_empty() && frame_ungapped_hits.is_empty() {
                frame_ungapped_hits = hits;
            } else if !hits.is_empty() {
                merge_tblastx_subject_chunk_hits(
                    &mut frame_ungapped_hits,
                    hits,
                    chunk.offset,
                    chunk.overlap,
                );
            }
        }
        x_merge_frame(&mut combined_ungapped_hits, frame_ungapped_hits);
    }
    KM_CHUNKS.fetch_add(chunks_merged, AtomicOrdering::Relaxed);
    combined_ungapped_hits
}

/// `XPairPar::merge_subject` with the switch on (`mode` 1: the units sorted their lists) or in
/// the shadow mode (`mode` 2: they did not).
pub(super) fn merge_subject(x_units: XSubjectUnits, inp: &XSubject<'_>) -> Vec<UngappedHit> {
    if mode() == 2 {
        return shadow(x_units, inp);
    }
    let t_merge = timing().then(Instant::now);
    let XSubjectUnits {
        frame_chunks,
        units,
        outputs,
    } = x_units;
    SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
    UNITS.fetch_add(units.len() as u64, AtomicOrdering::Relaxed);
    KM_SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
    let combined = match outputs {
        XUnitOutputs::Full(lists) => x_merge(&frame_chunks, &units, lists, inp, true),
        XUnitOutputs::Compact(records) => {
            // The LOSAT_X_SEEDCOMPACT counters of `merge_subject_compact`.
            SC_SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
            SC_UNITS.fetch_add(units.len() as u64, AtomicOrdering::Relaxed);
            SC_RECORDS.fetch_add(
                records.iter().map(|r| r.len() as u64).sum::<u64>(),
                AtomicOrdering::Relaxed,
            );
            x_merge(&frame_chunks, &units, records, inp, true)
        }
    };
    add_ns(&NS_MERGE, t_merge);
    combined
}

/// The shadow mode: the list of the switch off (on a copy of the units' lists, with its counters)
/// against this file's, every field; returns this file's.
fn shadow(x_units: XSubjectUnits, inp: &XSubject<'_>) -> Vec<UngappedHit> {
    let XSubjectUnits {
        frame_chunks,
        units,
        outputs,
    } = x_units;
    let copy = XSubjectUnits {
        frame_chunks: frame_chunks.clone(),
        units: units.clone(),
        outputs: match &outputs {
            XUnitOutputs::Full(lists) => XUnitOutputs::Full(lists.clone()),
            XUnitOutputs::Compact(records) => XUnitOutputs::Compact(records.clone()),
        },
    };
    let reference = XPairPar::merge_subject_listed(copy, inp);
    KM_SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
    // This file's merge, without the counters of x_pair_par.rs (the reference counted the subject).
    let merged = match outputs {
        XUnitOutputs::Full(mut lists) => {
            lists.iter_mut().for_each(XUnitSort::x_sort);
            x_merge(&frame_chunks, &units, lists, inp, false)
        }
        XUnitOutputs::Compact(mut records) => {
            // the task's sort (with a buffer) against `slice::sort_by`, unit by unit
            let mut buffer = Vec::new();
            for (k, list) in records.iter_mut().enumerate() {
                let mut by_std = list.clone();
                by_std.x_sort();
                x_stable_sort_with(list, &mut buffer, <XUnitHsp as XRecord>::cmp);
                assert!(
                    *list == by_std,
                    "LOSAT_X_TBXKMERGESHADOW: subject {} unit {k}: the buffered sort differs from slice::sort_by",
                    inp.s_idx
                );
            }
            x_merge(&frame_chunks, &units, records, inp, false)
        }
    };
    let s_idx = inp.s_idx;
    for k in 0..reference.len().max(merged.len()) {
        match (merged.get(k), reference.get(k)) {
            (Some(a), Some(b)) if x_same_hit(a, b) => {}
            (a, b) => panic!(
                "LOSAT_X_TBXKMERGESHADOW: subject {s_idx} index {k} differs (merged {} hits, reference {} hits): merged {a:?} vs reference {b:?}",
                merged.len(),
                reference.len()
            ),
        }
    }
    KM_SHADOW_SUBJECTS.fetch_add(1, AtomicOrdering::Relaxed);
    KM_SHADOW_HSPS.fetch_add(reference.len() as u64, AtomicOrdering::Relaxed);
    merged
}

/// The summary lines of the switch (printed by `x_pair_par::print_summary`).
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub(super) fn print_summary() {
    let load = |c: &AtomicU64| c.load(AtomicOrdering::Relaxed);
    let m = mode();
    if m == 2 {
        eprintln!(
            "[X_SHADOW] TBXKMERGE subjects_compared={} hsps_compared={} (combined_ungapped_hits identical)",
            load(&KM_SHADOW_SUBJECTS),
            load(&KM_SHADOW_HSPS)
        );
    }
    if m != 0 && std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] TBXKMERGE subjects={} chunks={} tie_runs={} k1_fallback={} frame_unsorted={}",
            load(&KM_SUBJECTS),
            load(&KM_CHUNKS),
            load(&KM_TIE_RUNS),
            load(&KM_K1_FALLBACK),
            load(&KM_FRAME_UNSORTED)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small deterministic generator (xorshift), so the tests need no crate.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn below(&mut self, n: u64) -> i32 {
            (self.next() % n) as i32
        }
    }

    fn record(rng: &mut Rng) -> XUnitHsp {
        // Few distinct values, so that ties on every key are common.
        XUnitHsp {
            q_start: rng.below(6),
            s_start: rng.below(4),
            len: 1 + rng.below(3),
            q_seed: rng.below(1000),
            s_seed: rng.below(1000),
            score: rng.below(3),
            ctx_idx: rng.below(2) as u32,
        }
    }

    fn init(r: &XUnitHsp) -> InitHSP {
        InitHSP {
            q_start_absolute: r.q_start,
            q_end_absolute: r.q_start + r.len,
            s_start: r.s_start,
            s_end: r.s_start + r.len,
            q_seed_absolute: r.q_seed,
            s_seed: r.s_seed,
            score: r.score,
            ctx_idx: r.ctx_idx as usize,
            s_f_idx: 0,
            q_idx: 0,
            s_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_orig_len: 100,
            s_orig_len: 100,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-304
    // ```c
    //     if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
    //                                  h1->ungapped_data->score)) &&
    // ```
    // The record's comparator is the init comparator on the expanded records.
    #[test]
    fn the_record_order_is_the_init_order() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        for _ in 0..20_000 {
            let (a, b) = (record(&mut rng), record(&mut rng));
            assert_eq!(
                <XUnitHsp as XRecord>::cmp(&a, &b),
                score_compare_init_hsps_ncbi(&init(&a), &init(&b)),
                "{a:?} {b:?}"
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
    // ```c
    //     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
    //           sizeof(BlastInitHSP), score_compare_match);
    // ```
    // Argument (a): the k-way merge of the units' sorted lists is the stable sort of their
    // concatenation, record for record (the seed offsets tell equal records apart), with empty
    // lists among them.
    #[test]
    fn the_kway_merge_is_the_stable_sort_of_the_concatenation() {
        let mut rng = Rng(0x2545_f491_4f6c_dd1d);
        for case in 0..2_000 {
            let k = 1 + (case % 6);
            let lists: Vec<Vec<InitHSP>> = (0..k)
                .map(|_| {
                    let n = rng.below(12) as usize;
                    (0..n).map(|_| init(&record(&mut rng))).collect()
                })
                .collect();
            let mut expected: Vec<InitHSP> = lists.concat();
            expected.sort_by(score_compare_init_hsps_ncbi);
            let mut sorted = lists.clone();
            sorted.iter_mut().for_each(XUnitSort::x_sort);
            let unit_free = |lists: Vec<Vec<InitHSP>>| {
                // `InitHSP::init_hsp` does not read the unit; any `XExpand` of an empty subject.
                let frames: [QueryFrame; 0] = [];
                let inp = XSubject {
                    s_frames_preliminary: &frames,
                    s_frames: &frames,
                    contexts: &[],
                    cutoff_scores: &[],
                    x_dropoff_per_context: &[],
                    window: 40,
                    wordsize: 3,
                    s_idx: 0,
                    s_len: 0,
                };
                let unit = XExpand {
                    f: 0,
                    s_frame: 1,
                    s_idx: 0,
                    inp: &inp,
                };
                x_kway(lists, &unit)
            };
            assert_eq!(unit_free(sorted), expected, "case {case}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
    // ```c
    //     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
    //           sizeof(BlastInitHSP), score_compare_match);
    // ```
    // The unit's sort in its task is `slice::sort_by` (stable), record for record, at every length
    // around the run and pass boundaries, with a buffer of any capacity and contents.
    #[test]
    fn the_buffered_sort_is_the_stable_sort() {
        let mut rng = Rng(0x0bad_5eed_1234_5678);
        for n in (0..70).chain([127, 128, 129, 255, 256, 257, 1000, 4099]) {
            for _ in 0..4 {
                // seed offsets tell equal records apart
                let list: Vec<XUnitHsp> = (0..n)
                    .map(|i| XUnitHsp {
                        q_seed: i as i32,
                        ..record(&mut rng)
                    })
                    .collect();
                let mut expected = list.clone();
                expected.sort_by(<XUnitHsp as XRecord>::cmp);
                let mut sorted = list.clone();
                let mut buf: Vec<XUnitHsp> = (0..rng.below(8)).map(|_| record(&mut rng)).collect();
                x_stable_sort_with(&mut sorted, &mut buf, <XUnitHsp as XRecord>::cmp);
                assert_eq!(sorted, expected, "n {n}");
            }
        }
    }

    fn hit(rng: &mut Rng, order: u32) -> UngappedHit {
        let s = rng.below(4) as u32;
        let q = rng.below(4) as u32;
        let len = 1 + rng.below(2) as u32;
        UngappedHit {
            q_idx: 0,
            s_idx: 0,
            ctx_idx: rng.below(3) as u32,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: q,
            q_aa_end: q + len,
            s_aa_start: s,
            s_aa_end: s + len,
            q_seed_off: q,
            s_seed_off: s,
            q_orig_len: 100,
            s_orig_len: 100,
            raw_score: rng.below(3),
            e_value: f64::INFINITY,
            num_ident: 0,
            hsp_list_order: order,
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
            link_id: 0,
            chain_next_link_id: None,
            hsp_link_num: 0,
            num: 0,
        }
    }

    fn same_lists(a: &[UngappedHit], b: &[UngappedHit]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x_same_hit(x, y))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
    // ```c
    //     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
    //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //               ScoreCompareHSPs);
    // ```
    // Argument (c): on a K1-sorted list the run-wise sort is the full stable sort; on a list out of
    // K1 order the full sort runs. `hsp_list_order` tells equal HSPs apart.
    #[test]
    fn the_run_wise_sort_is_the_full_sort() {
        let mut rng = Rng(0x0123_4567_89ab_cdef);
        for case in 0..2_000 {
            let n = rng.below(40) as u32;
            let mut list: Vec<UngappedHit> = (0..n).map(|i| hit(&mut rng, i)).collect();
            if case % 2 == 0 {
                list.sort_by(x_k1);
            }
            let mut expected = list.clone();
            sort_ungapped_hits_by_score_ncbi(&mut expected);
            x_sort_k1_runs(&mut list);
            assert!(same_lists(&list, &expected), "case {case}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766
    // ```c
    // combined_hsp_list->hspcnt = new_hspcnt;
    // Blast_HSPListSortByScore(combined_hsp_list);
    // ```
    // Argument (d): the backward merge is the append and the stable sort, also for an unsorted or
    // empty frame list and an empty combined list.
    #[test]
    fn the_frame_merge_is_the_append_and_the_sort() {
        let mut rng = Rng(0xdead_beef_cafe_f00d);
        for case in 0..2_000 {
            let mut combined: Vec<UngappedHit> = (0..rng.below(30) as u32)
                .map(|i| hit(&mut rng, i))
                .collect();
            sort_ungapped_hits_by_score_ncbi(&mut combined);
            let mut frame: Vec<UngappedHit> = (0..rng.below(30) as u32)
                .map(|i| hit(&mut rng, 100 + i))
                .collect();
            if case % 3 != 0 {
                sort_ungapped_hits_by_score_ncbi(&mut frame);
            }
            let mut expected = combined.clone();
            expected.extend(frame.iter().cloned());
            if !ungapped_hits_is_sorted_by_score_ncbi(&expected) {
                sort_ungapped_hits_by_score_ncbi(&mut expected);
            }
            x_merge_frame(&mut combined, frame);
            assert!(same_lists(&combined, &expected), "case {case}");
        }
    }
}
