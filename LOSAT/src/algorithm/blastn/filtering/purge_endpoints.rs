//! HSP endpoint purging for blastn - removes or trims HSPs with common start/end points
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535
//! Function: Blast_HSPListPurgeHSPsWithCommonEndpoints
//!
//! NCBI's endpoint purging has two modes:
//! - `purge` disabled: Trim overlapping HSPs using gap_info (edit script)
//! - `purge` enabled: Delete HSPs with common endpoints entirely
//!
//! The blastn traceback flow is:
//! 1. Call with `purge` disabled first (trim overlapping HSPs)
//! 2. Re-evaluate trimmed fragments
//! 3. Call with `purge` enabled for final cleanup
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:637-669

use super::super::alignment::ContextQuery;
use super::super::hsp::{adjust_blastn_offsets, BlastnHsp};
use crate::utils::glibc_fma_libm::GlibcLibm;
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
use super::super::tracing as blastn_trace;
use crate::common::GapEditOp;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::num::NonZeroUsize;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2484-2487
fn blastn_hsp_context(hsp: &BlastnHsp) -> u32 {
    hsp.q_idx * 2 + if hsp.query_frame < 0 { 1 } else { 0 }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2286-2320
fn compare_query_offset_hsps_for_common_endpoint(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    blastn_hsp_context(a)
        .cmp(&blastn_hsp_context(b))
        .then_with(|| a.internal_q_offset_0.cmp(&b.internal_q_offset_0))
        .then_with(|| a.internal_s_offset_0.cmp(&b.internal_s_offset_0))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.internal_q_end_0.cmp(&a.internal_q_end_0))
        .then_with(|| b.internal_s_end_0.cmp(&a.internal_s_end_0))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2351-2386
fn compare_query_end_hsps_for_common_endpoint(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    blastn_hsp_context(a)
        .cmp(&blastn_hsp_context(b))
        .then_with(|| a.internal_q_end_0.cmp(&b.internal_q_end_0))
        .then_with(|| a.internal_s_end_0.cmp(&b.internal_s_end_0))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.internal_q_offset_0.cmp(&a.internal_q_offset_0))
        .then_with(|| b.internal_s_offset_0.cmp(&a.internal_s_offset_0))
}

// =============================================================================
// Traceback Trimming Functions
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2392-2452 s_CutOffGapEditScript
// =============================================================================

/// Cut off gap edit script at specified position.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2392-2452 s_CutOffGapEditScript
///
/// This function trims an HSP's edit script at a specified cut point.
/// - If cut_begin=true: Keep the END portion, update offsets
/// - If cut_begin=false: Keep the START portion, update ends
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2391-2416
///
/// # Arguments
/// * `hit` - HSP to modify in place
/// * `q_cut` - Query cut position (absolute, 0-based internal offset)
/// * `s_cut` - Subject cut position (absolute, 0-based internal offset)
/// * `cut_begin` - If true, keep end portion; if false, keep start portion
///
/// # Returns
/// * `true` if trimming was successful
/// * `false` if no gap_info or cut point not found
pub fn cut_off_gap_edit_script(
    hit: &mut BlastnHsp,
    q_cut: i32,
    s_cut: usize,
    cut_begin: bool,
) -> bool {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2392-2396
    // NCBI mutates the HSP's existing edit script in place. Move the Vec out
    // while trimming so Rust does not clone the whole script in this hot path.
    let mut gap_info = match hit.gap_info.take() {
        Some(info) if !info.is_empty() => info,
        other => {
            hit.gap_info = other;
            return false;
        }
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2396-2401
    // Endpoint trimming operates on canonical internal offsets already stored
    // in the HSP. Do not reconstruct them from output coordinates.
    let mut q_offset_0 = hit.internal_q_offset_0;
    let mut q_end_0 = hit.internal_q_end_0;
    let mut s_offset_0 = hit.internal_s_offset_0;
    let mut s_end_0 = hit.internal_s_end_0;

    // Convert absolute cut positions to relative (within HSP)
    // Subtract the HSP's starting query and subject offsets from `q_cut` and `s_cut`.
    // NCBI's `Int4` difference; one below 0 tests as 0 (`qid` and `sid` are at least 0 when
    // they are tested).
    let q_cut_rel = (q_cut - q_offset_0).max(0) as usize;
    let s_cut_rel = s_cut.saturating_sub(s_offset_0);

    // Walk through edit script to find cut point
    // NCBI tracks: qid (query consumed), sid (subject consumed), opid (ops within current run)
    let mut qid: usize = 0;
    let mut sid: usize = 0;
    let mut found = false;
    let mut found_index: usize = 0;
    let mut found_opid: usize = 0;

    // Iterate over every run in the edit script.
    'outer: for (index, op) in gap_info.iter().enumerate() {
        let num = op.num() as usize;

        // Advance through the operations in the current edit-script run.
        match op {
            GapEditOp::Sub(_) => {
                // Sub: iterate position by position to find exact cut point
                // Increment `qid`, `sid`, and `opid` once per aligned substitution base.
                for opid in 0..num {
                    qid += 1;
                    sid += 1;

                    // Stop as soon as both `qid` and `sid` reach or pass their cut targets.
                    if qid >= q_cut_rel && sid >= s_cut_rel {
                        found = true;
                        found_index = index;
                        found_opid = opid + 1; // Position after the cut
                        break 'outer;
                    }
                }
            }
            GapEditOp::Del(n) => {
                // Del: gap in query, consumes only subject
                // Advance both `sid` and `opid` by the full deletion run length.
                sid += *n as usize;

                // Check if cut point reached after this run
                if qid >= q_cut_rel && sid >= s_cut_rel {
                    found = true;
                    found_index = index;
                    found_opid = *n as usize; // Entire run consumed
                    break 'outer;
                }
            }
            GapEditOp::Ins(n) => {
                // Ins: gap in subject, consumes only query
                // Advance both `qid` and `opid` by the full insertion run length.
                qid += *n as usize;

                // Check if cut point reached after this run
                if qid >= q_cut_rel && sid >= s_cut_rel {
                    found = true;
                    found_index = index;
                    found_opid = *n as usize; // Entire run consumed
                    break 'outer;
                }
            }
        }
    }

    if !found {
        hit.gap_info = Some(gap_info);
        return false;
    }

    if cut_begin {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2426-2439
        // Keep the tail of the same edit script buffer.
        let drain_end = match gap_info[found_index] {
            GapEditOp::Sub(n) if found_opid < n as usize => {
                gap_info[found_index] = GapEditOp::Sub((n as usize - found_opid) as u32);
                found_index
            }
            _ => found_index.saturating_add(1),
        };
        if drain_end > 0 {
            gap_info.drain(0..drain_end);
        }

        // Update HSP coordinates
        // Advance the HSP's query and subject starting offsets by `qid` and `sid`.
        q_offset_0 += qid as i32;
        s_offset_0 = s_offset_0.saturating_add(sid);
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2442-2449
        // Keep the prefix and truncate the same edit script buffer.
        if let GapEditOp::Sub(n) = gap_info[found_index] {
            if found_opid < n as usize {
                gap_info[found_index] = GapEditOp::Sub(found_opid as u32);
            }
        }
        gap_info.truncate(found_index.saturating_add(1));

        // Update HSP coordinates
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2449-2449
        q_end_0 = q_offset_0 + qid as i32;
        s_end_0 = s_offset_0.saturating_add(sid);
    }

    let (q_start, q_end, s_start, s_end) = adjust_blastn_offsets(
        q_offset_0,
        q_end_0,
        s_offset_0,
        s_end_0,
        hit.query_length,
        hit.query_frame,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2436-2449
    hit.internal_q_offset_0 = q_offset_0;
    hit.internal_q_end_0 = q_end_0;
    hit.internal_s_offset_0 = s_offset_0;
    hit.internal_s_end_0 = s_end_0;
    hit.q_start = q_start;
    hit.q_end = q_end;
    hit.s_start = s_start;
    hit.s_end = s_end;

    // Update gap_info
    hit.gap_info = if gap_info.is_empty() {
        None
    } else {
        Some(gap_info)
    };

    // Update alignment length based on new gap_info
    if let Some(ref ops) = hit.gap_info {
        let mut new_len = 0usize;
        for op in ops {
            new_len += op.num() as usize;
        }
        hit.length = new_len;
    }

    true
}

/// Purge HSPs with common start or end positions.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535 Blast_HSPListPurgeHSPsWithCommonEndpoints
///
/// CRITICAL: NCBI processes HSPs per-subject (BlastHSPList is per-subject).
/// This function groups hits by subject index (oid) before purging to match NCBI behavior.
///
/// # Parameters
/// * `purge` - If true, delete HSPs entirely (current behavior).
///   If false, trim overlapping HSPs using gap_info and move to end.
///
/// For blastn:
/// - context = query index + strand (q_idx * 2 + {0,1})
/// - subject.frame is always +1 (blastn subject only), so strand is encoded by context.
///
/// IMPORTANT: Use internal offsets (query.offset/subject.offset), not output coordinates.
///
/// # Returns
/// Returns the index of the first trimmed HSP (for re-evaluation in purge=false mode).
/// In purge=true mode, this is always equal to the final hit count.
pub fn purge_hsps_with_common_endpoints_ex(
    hits: Vec<BlastnHsp>,
    purge: bool,
) -> (Vec<BlastnHsp>, usize) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:637-668
    // Traceback calls this on one BlastHSPList, i.e. one subject, in the hot path.
    // Avoid rebuilding the same per-subject grouping when the input already matches
    // NCBI's per-subject HSP list shape.
    if hits.iter().all(|h| h.s_idx == hits[0].s_idx) {
        return purge_hsps_for_subject_ex(hits, purge);
    }

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535
    // Blast_HSPListPurgeHSPsWithCommonEndpoints operates on BlastHSPList which is per-subject.
    // We must group hits by subject oid (s_idx) and process each group separately.

    // Group hits by subject oid (s_idx)
    let mut subject_groups: FxHashMap<u32, Vec<BlastnHsp>> = FxHashMap::default();
    for hit in hits {
        subject_groups.entry(hit.s_idx).or_default().push(hit);
    }

    // Process each subject group independently and collect results
    let mut result: Vec<BlastnHsp> = Vec::new();
    let mut total_extra_start = 0usize;

    for (_subject_idx, group_hits) in subject_groups {
        let (purged, extra_start) = purge_hsps_for_subject_ex(group_hits, purge);
        let offset = result.len();
        result.extend(purged);
        // Adjust extra_start to account for position in combined result
        if extra_start > 0 {
            // extra_start is relative to the group, we want position in final result
            total_extra_start = offset + extra_start;
        }
    }

    // If no trimming occurred, extra_start equals final count
    if total_extra_start == 0 {
        total_extra_start = result.len();
    }

    (result, total_extra_start)
}

/// Backward-compatible wrapper that always uses purge=true mode.
pub fn purge_hsps_with_common_endpoints(hits: Vec<BlastnHsp>) -> Vec<BlastnHsp> {
    let (result, _) = purge_hsps_with_common_endpoints_ex(hits, true);
    result
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
// A nonzero slot+1 is the Rust equivalent of a nullable HSP pointer.
type HspHandle = Option<NonZeroUsize>;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2482-2487
// Both slots in `hsp_array` must be live and belong to the same query `context`.
#[inline]
fn hsp_for_handle(hsp_storage: &[Option<BlastnHsp>], handle: HspHandle) -> Option<&BlastnHsp> {
    handle.map(|index| {
        hsp_storage[index.get() - 1]
            .as_ref()
            .expect("live HSP handle must reference a populated arena slot")
    })
}

/// No NCBI counterpart (memory placement only): the HSPs of the live handles of `hsp_array`, in
/// its order, built in the allocation of `storage` (the input list's, through `map(Some)`):
/// slot `k` takes the HSP of the `k`-th live handle by swaps, with the position of each slot's HSP
/// tracked, and the slots past the live ones are dropped. `Option<BlastnHsp>` has the size of
/// `BlastnHsp` (asserted below), so the `map` and `collect` at both ends keep the allocation too,
/// and the endpoint purge holds one array of HSPs instead of the input and a copy of its survivors
/// (Sakai x EDL933 `-subject` blastn: 38,696 HSPs, 7 MiB). A missing or repeated handle is an
/// error, as in the copy it replaces.
fn x_in_handle_order(
    mut storage: Vec<Option<BlastnHsp>>,
    hsp_array: &[HspHandle],
) -> Vec<BlastnHsp> {
    let n = storage.len();
    // pos_of[s]: where the HSP of slot s is now; slot_at[p]: the slot whose HSP is at p
    let mut pos_of: Vec<usize> = (0..n).collect();
    let mut slot_at: Vec<usize> = (0..n).collect();
    let mut k = 0usize;
    for handle in hsp_array.iter().flatten() {
        let s = handle.get() - 1;
        let p = pos_of[s];
        assert!(
            p >= k && storage[p].is_some(),
            "live output handle requires a unique HSP"
        );
        if p != k {
            storage.swap(k, p);
            let t = slot_at[k];
            slot_at[p] = t;
            pos_of[t] = p;
            slot_at[k] = s;
            pos_of[s] = k;
        }
        k += 1;
    }
    storage.truncate(k);
    storage
        .into_iter()
        .map(|hsp| hsp.expect("live output handle requires a unique HSP"))
        .collect()
}

// No NCBI counterpart: `x_in_handle_order` keeps the allocation only when the two element types have
// one layout (the niche of `BlastnHsp`'s fields gives `Option` no extra tag).
const _: () = assert!(
    std::mem::size_of::<Option<BlastnHsp>>() == std::mem::size_of::<BlastnHsp>()
        && std::mem::align_of::<Option<BlastnHsp>>() == std::mem::align_of::<BlastnHsp>()
);

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2457,2478
/// Purge HSPs with common endpoints for a single subject.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535
/// This matches NCBI's per-BlastHSPList processing.
///
/// # Parameters
/// * `purge` - If true, delete HSPs entirely. If false, trim and move to end of array.
///
/// # Returns
/// Tuple of (result hits, index of first trimmed HSP for re-evaluation)
fn purge_hsps_for_subject_ex(hits: Vec<BlastnHsp>, purge: bool) -> (Vec<BlastnHsp>, usize) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len);
    }

    let initial_count = len;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2471-2478
    // Keep HSP payloads in fixed slots while sorting/shifting nullable handles.
    // This arena never grows or reorders; only final materialization moves live HSPs.
    let mut hsp_storage: Vec<Option<BlastnHsp>> = hits.into_iter().map(Some).collect();
    let mut hsp_array: Vec<HspHandle> =
        (0..len).map(|index| NonZeroUsize::new(index + 1)).collect();
    let mut hsp_count = initial_count;
    let mut start_removed = 0usize;
    let mut end_removed = 0usize;
    let mut start_trimmed = 0usize;
    let mut end_trimmed = 0usize;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2471-2478
    qsort_hsp_array_prefix_by(
        &hsp_storage,
        &mut hsp_array,
        hsp_count,
        compare_query_offset_hsps_for_common_endpoint,
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2480-2500
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2500
    // The C code keeps the sorted active prefix and moves trimmed/null HSPs to
    // the tail by decrementing hsp_count before shifting hsp_array left.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2502
    // No NCBI counterpart: LOSAT makes the same keep/trim/delete decisions on the same pairs as the C
    // loop above. It writes the survivors in order and the removed entries to the tail slots the C
    // shifts would leave them in, so the array is the same; it only avoids the repeated left shifts.
    // NCBI removes hsp_array[i+1], shifts the rest of the active prefix left and parks the removed
    // pointer at the shrinking end, so the k-th removal of a pass ends at index (count at pass
    // start - 1 - k) and the survivors stay in order. This loop compares the current keeper
    // with the next unprocessed entry exactly as NCBI does.
    let n0 = hsp_count;
    if n0 > 0 {
        let mut w = 0usize;
        let mut tails: Vec<HspHandle> = Vec::new();
        for r in 1..n0 {
            if same_common_start(&hsp_storage, &hsp_array, w, r) {
                let (keeper_q_end, keeper_s_end) = {
                    let keeper = hsp_for_handle(&hsp_storage, hsp_array[w])
                        .expect("same_common_start requires a keeper HSP");
                    (keeper.internal_q_end_0, keeper.internal_s_end_0)
                };
                let removed = hsp_array[r].take();
                let mut tail_hsp = None;

                if let Some(handle) = removed {
                    let removed_hit = hsp_for_handle(&hsp_storage, Some(handle))
                        .expect("removed handle requires a live HSP");
                    let keeper = hsp_for_handle(&hsp_storage, hsp_array[w])
                        .expect("common endpoint requires a keeper HSP");
                    if !purge && removed_hit.internal_q_end_0 > keeper_q_end {
                        trace_common_endpoint_purge_decision(
                            "common_start",
                            "trim_begin",
                            purge,
                            keeper,
                            removed_hit,
                            keeper_q_end,
                            keeper_s_end,
                        );
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2490-2492,2516-2518
                        // End immutable keeper/removed borrows before editing this one slot.
                        let removed_hit = hsp_storage[handle.get() - 1]
                            .as_mut()
                            .expect("trimmed handle requires a live HSP");
                        let _ =
                            cut_off_gap_edit_script(removed_hit, keeper_q_end, keeper_s_end, true);
                        tail_hsp = Some(handle);
                        start_trimmed += 1;
                    } else {
                        trace_common_endpoint_purge_decision(
                            "common_start",
                            "delete",
                            purge,
                            keeper,
                            removed_hit,
                            keeper_q_end,
                            keeper_s_end,
                        );
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2494,2520
                        // Release owned script storage now, before shifting handles.
                        drop(
                            hsp_storage[handle.get() - 1]
                                .take()
                                .expect("deleted handle requires a live HSP"),
                        );
                    }
                    start_removed += 1;
                }
                tails.push(tail_hsp);
            } else {
                w += 1;
                if w != r {
                    hsp_array[w] = hsp_array[r];
                }
            }
        }
        hsp_count = w + 1;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2496-2499,2522-2525
        // The removed handles go to the tail slots, the first removal to the last slot.
        for (k, tail) in tails.into_iter().enumerate() {
            hsp_array[n0 - 1 - k] = tail;
        }
    }

    // Pass 2: Remove HSPs with common END positions
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2504-2526
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2504
    qsort_hsp_array_prefix_by(
        &hsp_storage,
        &mut hsp_array,
        hsp_count,
        compare_query_end_hsps_for_common_endpoint,
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2499-2528
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2512-2528
    // Repeat the same active-prefix mutation NCBI uses for common endpoints.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2504-2528
    // No NCBI counterpart: LOSAT makes the same keep/trim/delete decisions on the same pairs as the C
    // loop above. It writes the survivors in order and the removed entries to the tail slots the C
    // shifts would leave them in, so the array is the same; it only avoids the repeated left shifts.
    // NCBI removes hsp_array[i+1], shifts the rest of the active prefix left and parks the removed
    // pointer at the shrinking end, so the k-th removal of a pass ends at index (count at pass
    // start - 1 - k) and the survivors stay in order. This loop compares the current keeper
    // with the next unprocessed entry exactly as NCBI does.
    let n0 = hsp_count;
    if n0 > 0 {
        let mut w = 0usize;
        let mut tails: Vec<HspHandle> = Vec::new();
        for r in 1..n0 {
            if same_common_end(&hsp_storage, &hsp_array, w, r) {
                let (keeper_q_offset, keeper_s_offset) = {
                    let keeper = hsp_for_handle(&hsp_storage, hsp_array[w])
                        .expect("same_common_end requires a keeper HSP");
                    (keeper.internal_q_offset_0, keeper.internal_s_offset_0)
                };
                let removed = hsp_array[r].take();
                let mut tail_hsp = None;

                if let Some(handle) = removed {
                    let removed_hit = hsp_for_handle(&hsp_storage, Some(handle))
                        .expect("removed handle requires a live HSP");
                    let keeper = hsp_for_handle(&hsp_storage, hsp_array[w])
                        .expect("common endpoint requires a keeper HSP");
                    if !purge && removed_hit.internal_q_offset_0 < keeper_q_offset {
                        trace_common_endpoint_purge_decision(
                            "common_end",
                            "trim_end",
                            purge,
                            keeper,
                            removed_hit,
                            keeper_q_offset,
                            keeper_s_offset,
                        );
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2490-2492,2516-2518
                        // End immutable keeper/removed borrows before editing this one slot.
                        let removed_hit = hsp_storage[handle.get() - 1]
                            .as_mut()
                            .expect("trimmed handle requires a live HSP");
                        let _ = cut_off_gap_edit_script(
                            removed_hit,
                            keeper_q_offset,
                            keeper_s_offset,
                            false,
                        );
                        tail_hsp = Some(handle);
                        end_trimmed += 1;
                    } else {
                        trace_common_endpoint_purge_decision(
                            "common_end",
                            "delete",
                            purge,
                            keeper,
                            removed_hit,
                            keeper_q_offset,
                            keeper_s_offset,
                        );
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2494,2520
                        // Release owned script storage now, before shifting handles.
                        drop(
                            hsp_storage[handle.get() - 1]
                                .take()
                                .expect("deleted handle requires a live HSP"),
                        );
                    }
                    end_removed += 1;
                }
                tails.push(tail_hsp);
            } else {
                w += 1;
                if w != r {
                    hsp_array[w] = hsp_array[r];
                }
            }
        }
        hsp_count = w + 1;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2496-2499,2522-2525
        // The removed handles go to the tail slots, the first removal to the last slot.
        for (k, tail) in tails.into_iter().enumerate() {
            hsp_array[n0 - 1 - k] = tail;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2530-2535
    // In purge=false mode NCBI leaves trimmed HSPs in the tail for
    // blast_traceback.c:647-665 to re-evaluate, while NULL tail entries are
    // skipped. Store the same live order and return the active-prefix boundary.
    let extra_start = hsp_count;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2530-2535
    // No NCBI counterpart (memory placement only): the live handle order is applied to the storage
    // in place, so the list is never copied into a second array (`x_in_handle_order`).
    let hits = x_in_handle_order(hsp_storage, &hsp_array);

    // Debug output for purge statistics
    if crate::utils::xenv::debug_blastn_is_ok() {
        eprintln!(
            "[DEBUG_PURGE] initial={}, start_removed={} (trimmed={}), end_removed={} (trimmed={}), active_prefix={}, final={}",
            initial_count, start_removed, start_trimmed, end_removed, end_trimmed, extra_start, hits.len()
        );
    }

    let _ = (
        initial_count,
        start_removed,
        end_removed,
        start_trimmed,
        end_trimmed,
    ); // Silence unused warnings

    (hits, extra_start)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2278-2283,2343-2348
// NCBI's qsort callbacks sort NULL HSP pointers after real HSPs. The active
// prefix normally contains only real HSPs, but preserve the NULL ordering when
// mirroring Blast_HSPListPurgeHSPsWithCommonEndpoints exactly.
fn qsort_hsp_array_prefix_by(
    hsp_storage: &[Option<BlastnHsp>],
    hsp_array: &mut [HspHandle],
    hsp_count: usize,
    compare: fn(&BlastnHsp, &BlastnHsp) -> Ordering,
) {
    let prefix_len = hsp_count.min(hsp_array.len());
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478,2504,2268-2387
    // Retain the existing Rust stable tie order, without adding an index tie-break.
    hsp_array[..prefix_len].sort_by(|&lhs, &rhs| {
        match (
            hsp_for_handle(hsp_storage, lhs),
            hsp_for_handle(hsp_storage, rhs),
        ) {
            (Some(a), Some(b)) => compare(a, b),
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            (None, None) => Ordering::Equal,
        }
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2482-2487
fn same_common_start(
    hsp_storage: &[Option<BlastnHsp>],
    hsp_array: &[HspHandle],
    lhs: usize,
    rhs: usize,
) -> bool {
    let (Some(a), Some(b)) = (
        hsp_array
            .get(lhs)
            .and_then(|&handle| hsp_for_handle(hsp_storage, handle)),
        hsp_array
            .get(rhs)
            .and_then(|&handle| hsp_for_handle(hsp_storage, handle)),
    ) else {
        return false;
    };
    blastn_hsp_context(a) == blastn_hsp_context(b)
        && a.internal_q_offset_0 == b.internal_q_offset_0
        && a.internal_s_offset_0 == b.internal_s_offset_0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2508-2513
fn same_common_end(
    hsp_storage: &[Option<BlastnHsp>],
    hsp_array: &[HspHandle],
    lhs: usize,
    rhs: usize,
) -> bool {
    let (Some(a), Some(b)) = (
        hsp_array
            .get(lhs)
            .and_then(|&handle| hsp_for_handle(hsp_storage, handle)),
        hsp_array
            .get(rhs)
            .and_then(|&handle| hsp_for_handle(hsp_storage, handle)),
    ) else {
        return false;
    };
    blastn_hsp_context(a) == blastn_hsp_context(b)
        && a.internal_q_end_0 == b.internal_q_end_0
        && a.internal_s_end_0 == b.internal_s_end_0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2480-2499
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2508-2525
fn trace_common_endpoint_purge_decision(
    pass: &str,
    action: &str,
    purge: bool,
    keeper: &BlastnHsp,
    removed: &BlastnHsp,
    cut_q: i32,
    cut_s: usize,
) {
    if !blastn_trace::enabled() {
        return;
    }

    let keeper_context = keeper.q_idx * 2 + if keeper.query_frame < 0 { 1 } else { 0 };
    let removed_context = removed.q_idx * 2 + if removed.query_frame < 0 { 1 } else { 0 };
    let subject_id_for_numeric_filter = "";
    let should_trace_keeper = blastn_trace::should_trace_range(
        "purge",
        keeper_context,
        keeper.s_idx as usize,
        subject_id_for_numeric_filter,
        keeper.internal_q_offset_0,
        keeper.internal_q_end_0,
        keeper.internal_s_offset_0,
        keeper.internal_s_end_0,
        keeper.query_length,
        keeper.query_frame,
    );
    let should_trace_removed = blastn_trace::should_trace_range(
        "purge",
        removed_context,
        removed.s_idx as usize,
        subject_id_for_numeric_filter,
        removed.internal_q_offset_0,
        removed.internal_q_end_0,
        removed.internal_s_offset_0,
        removed.internal_s_end_0,
        removed.query_length,
        removed.query_frame,
    );
    if !should_trace_keeper && !should_trace_removed {
        return;
    }

    blastn_trace::log(
        "purge",
        format!(
            "subject_idx={} context={} pass={} action={} purge={} keeper=q{}..{} s{}..{} raw_score={} removed=q{}..{} s{}..{} raw_score={} cut=({}, {})",
            removed.s_idx,
            removed_context,
            pass,
            action,
            purge,
            keeper.internal_q_offset_0,
            keeper.internal_q_end_0,
            keeper.internal_s_offset_0,
            keeper.internal_s_end_0,
            keeper.raw_score,
            removed.internal_q_offset_0,
            removed.internal_q_end_0,
            removed.internal_s_offset_0,
            removed.internal_s_end_0,
            removed.raw_score,
            cut_q,
            cut_s
        ),
    );
}

// =============================================================================
// HSP Re-evaluation with Ambiguities
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:479-647 Blast_HSPReevaluateWithAmbiguitiesGapped
// =============================================================================

/// BLASTNA alphabet size.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1060-1068 (BLASTNA_SIZE usage)
const BLASTNA_SIZE: usize = 16;

/// Compute alignment stats from a gap edit script.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-818 (identity counting)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1055-1074 (gap counts)
fn stats_from_edit_ops(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    edit_ops: &[GapEditOp],
) -> (usize, usize, usize, usize, usize) {
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_letters = 0usize;
    let mut align_length = 0usize;

    let mut qi = q_start;
    let mut si = s_start;

    for op in edit_ops {
        let num = op.num() as usize;
        align_length += num;
        match *op {
            GapEditOp::Sub(_) => {
                for _ in 0..num {
                    if qi < q_seq.len() && si < s_seq.len() {
                        if q_seq[qi] == s_seq[si] {
                            matches += 1;
                        } else {
                            mismatches += 1;
                        }
                    }
                    qi += 1;
                    si += 1;
                }
            }
            GapEditOp::Del(_) => {
                gap_opens += 1;
                gap_letters += num;
                si += num;
            }
            GapEditOp::Ins(_) => {
                gap_opens += 1;
                gap_letters += num;
                qi += num;
            }
        }
    }

    (matches, mismatches, gap_opens, gap_letters, align_length)
}

/// Update identity/length stats and apply Blast_HSPTestIdentityAndLength filtering.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:656-663
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-836,1003-1025
pub fn blast_hsp_test_identity_and_length(
    hit: &mut BlastnHsp,
    query_nomask: ContextQuery<'_>,
    s_seq: &[u8],
    percent_identity: f64,
) -> bool {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:758-765
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:675-676,678-692
    // BLASTN re-evaluation and identity counting operate on canonical internal
    // offsets; output coordinates are adjusted only after pruning. NCBI's `query` is the
    // context's start in the concatenated query (`query_nomask`), read at the HSP's `Int4`
    // offset, outside the context for an HSP outside it (`ContextQuery`).
    let q_seq = query_nomask.buf;
    let q_offset = query_nomask.at(hit.internal_q_offset_0);
    let s_offset = hit.internal_s_offset_0;

    let (matches, mismatches, gap_opens, _gap_letters, align_length) =
        if let Some(ref ops) = hit.gap_info {
            stats_from_edit_ops(q_seq, s_seq, q_offset, s_offset, ops)
        } else {
            // Ungapped fallback (N.B. not expected for BLASTN reevaluation path).
            let align_length = (hit.internal_q_end_0 - hit.internal_q_offset_0).max(0) as usize;
            let mut matches = 0usize;
            let mut mismatches = 0usize;
            for i in 0..align_length {
                let qi = q_offset.saturating_add(i);
                let si = s_offset.saturating_add(i);
                if qi < q_seq.len() && si < s_seq.len() {
                    if q_seq[qi] == s_seq[si] {
                        matches += 1;
                    } else {
                        mismatches += 1;
                    }
                }
            }
            (matches, mismatches, 0usize, 0usize, align_length)
        };

    hit.length = align_length;
    hit.mismatch = mismatches;
    hit.gapopen = gap_opens;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-836
    hit.num_ident = matches;
    hit.num_positives = matches;
    hit.identity = if align_length > 0 {
        (matches as f64 / align_length as f64) * 100.0
    } else {
        0.0
    };

    // `Blast_HSPTest`; NCBI's `min_hit_length` is always 0 (`core::blast_hits::hsp_test`).
    crate::core::blast_hits::hsp_test(matches as i64, align_length as i64, percent_identity, 0)
}

/// Parameters needed for E-value and bit score recalculation after re-evaluation.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-250 s_HSPListPostTracebackUpdate
pub struct ReevalParams {
    pub karlin: crate::stats::KarlinParams,
    pub eff_searchsp: i64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1864-1870
    pub round_down_evalue_score: bool,
}

/// Extended version with optional Karlin parameters for E-value recalculation.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-250 Blast_HSPListGetEvalues, Blast_HSPListGetBitScores
pub fn reevaluate_hsp_with_ambiguities_gapped_ex(
    hit: &mut BlastnHsp,
    query: ContextQuery<'_>,
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
    score_matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
    reeval_params: Option<&ReevalParams>,
) -> bool {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:539-540
    let gap_ops = match hit.gap_info.as_mut() {
        Some(info) if !info.is_empty() => info,
        _ => return true, // No gap_info = delete HSP
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:528-535
    // `q` is the context's start in the concatenated query and `qlen` the context's length
    // (`ContextQuery`): the positions below are in the concatenated query, from the HSP's `Int4`
    // offset, which lies outside the context for an HSP outside it.
    let q_seq = query.buf;
    let q_offset = query.at(hit.internal_q_offset_0);
    let s_offset = hit.internal_s_offset_0;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:508-526 (factor and gap penalties)
    let mut factor = 1;
    let (gap_open_eval, gap_extend_eval) = if gap_open == 0 && gap_extend == 0 {
        if reward % 2 == 1 {
            factor = 2;
        }
        (0, (reward - 2 * penalty) * factor / 2)
    } else {
        (gap_open, gap_extend)
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:507-514
    let matrix = score_matrix;

    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // NCBI uses raw pointers without bounds checks (c++/src/algo/blast/core/blast_hits.c:528-559).
    // Pre-validate the alignment span so we can safely use unchecked loads.
    let mut q_pos_check = q_offset;
    let mut s_pos_check = s_offset;
    for op in gap_ops.iter() {
        match *op {
            GapEditOp::Sub(n) => {
                let n = n as usize;
                q_pos_check = q_pos_check.saturating_add(n);
                s_pos_check = s_pos_check.saturating_add(n);
            }
            GapEditOp::Del(n) => {
                s_pos_check = s_pos_check.saturating_add(n as usize);
            }
            GapEditOp::Ins(n) => {
                q_pos_check = q_pos_check.saturating_add(n as usize);
            }
        }
    }
    let bounds_safe = q_pos_check <= q_seq.len() && s_pos_check <= s_seq.len();
    let q_ptr = q_seq.as_ptr();
    let s_ptr = s_seq.as_ptr();

    let mut query_pos = q_offset;
    let mut subject_pos = s_offset;

    let mut best_q_start = query_pos;
    let mut best_q_end = query_pos;
    let mut best_s_start = subject_pos;
    let mut best_s_end = subject_pos;

    let mut current_q_start = query_pos;
    let mut current_s_start = subject_pos;

    let mut best_start_esp_index = 0usize;
    let mut best_end_esp_index = 0usize;
    let mut current_start_esp_index = 0usize;
    let mut best_end_esp_num: i32 = -1;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:562-608
    macro_rules! update_after_op {
        ($index:ident, $op_index:ident) => {{
            if sum < 0 {
                if $op_index < gap_ops[$index].num() as usize {
                    if let GapEditOp::Sub(n) = gap_ops[$index] {
                        gap_ops[$index] = GapEditOp::Sub((n as usize - $op_index) as u32);
                    }
                    current_start_esp_index = $index;
                    $op_index = 0;
                } else {
                    current_start_esp_index = $index + 1;
                }
                sum = 0;
                current_q_start = query_pos;
                current_s_start = subject_pos;

                if score < cutoff_score {
                    best_q_start = query_pos;
                    best_s_start = subject_pos;
                    score = 0;
                    best_start_esp_index = current_start_esp_index;
                    best_end_esp_index = current_start_esp_index;
                }
            } else if sum > score {
                score = sum;
                best_q_start = current_q_start;
                best_s_start = current_s_start;
                best_q_end = query_pos;
                best_s_end = subject_pos;
                best_start_esp_index = current_start_esp_index;
                best_end_esp_index = $index;
                best_end_esp_num = $op_index as i32;
            }
        }};
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:541-610
    for index in 0..gap_ops.len() {
        let mut op_index = 0usize;
        while op_index < gap_ops[index].num() as usize {
            match gap_ops[index] {
                GapEditOp::Sub(_) => {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:547-551
                    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
                    {
                        use std::arch::wasm32::*;
                        let remaining = gap_ops[index].num() as usize - op_index;
                        if remaining >= 16
                            && query_pos + 16 <= q_seq.len()
                            && subject_pos + 16 <= s_seq.len()
                        {
                            // Safety: bounds are checked above; v128_load reads 16 bytes.
                            let (qv, sv) = unsafe {
                                (
                                    v128_load(q_seq.as_ptr().add(query_pos) as *const v128),
                                    v128_load(s_seq.as_ptr().add(subject_pos) as *const v128),
                                )
                            };
                            macro_rules! step_lane {
                                ($lane:expr, $qv:ident, $sv:ident) => {{
                                    let q = i8x16_extract_lane::<$lane>($qv) as u8;
                                    let s = i8x16_extract_lane::<$lane>($sv) as u8;
                                    let q_code = (q & 0x0f) as usize;
                                    let s_code = (s & 0x0f) as usize;
                                    // NCBI `kResidueMask` is `0x0f`, giving indices 0..15.
                                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:510-511
                                    // SAFETY: `q_code` and `s_code` are `& 0x0f`, so the index
                                    // `q_code * 16 + s_code` is at most 255, inside `matrix`
                                    // (`&[i32; BLASTNA_SIZE * BLASTNA_SIZE]`, 256 entries).
                                    let score = unsafe {
                                        *matrix.get_unchecked(q_code * BLASTNA_SIZE + s_code)
                                    };
                                    sum += factor * score;
                                    query_pos += 1;
                                    subject_pos += 1;
                                    op_index += 1;
                                    update_after_op!(index, op_index);
                                }};
                            }
                            step_lane!(0, qv, sv);
                            step_lane!(1, qv, sv);
                            step_lane!(2, qv, sv);
                            step_lane!(3, qv, sv);
                            step_lane!(4, qv, sv);
                            step_lane!(5, qv, sv);
                            step_lane!(6, qv, sv);
                            step_lane!(7, qv, sv);
                            step_lane!(8, qv, sv);
                            step_lane!(9, qv, sv);
                            step_lane!(10, qv, sv);
                            step_lane!(11, qv, sv);
                            step_lane!(12, qv, sv);
                            step_lane!(13, qv, sv);
                            step_lane!(14, qv, sv);
                            step_lane!(15, qv, sv);
                            continue;
                        }
                    }
                    if bounds_safe {
                        // Safety: bounds were validated above and kResidueMask keeps indices in range.
                        let q_code = unsafe { *q_ptr.add(query_pos) & 0x0f } as usize;
                        let s_code = unsafe { *s_ptr.add(subject_pos) & 0x0f } as usize;
                        let score =
                            unsafe { *matrix.get_unchecked(q_code * BLASTNA_SIZE + s_code) };
                        sum += factor * score;
                    } else if query_pos < q_seq.len() && subject_pos < s_seq.len() {
                        let q_code = (q_seq[query_pos] & 0x0f) as usize;
                        let s_code = (s_seq[subject_pos] & 0x0f) as usize;
                        sum += factor * matrix[q_code * BLASTNA_SIZE + s_code];
                    }
                    query_pos += 1;
                    subject_pos += 1;
                    op_index += 1;
                }
                GapEditOp::Del(n) => {
                    sum -= gap_open_eval + gap_extend_eval * (n as i32);
                    subject_pos += n as usize;
                    op_index += n as usize;
                }
                GapEditOp::Ins(n) => {
                    sum -= gap_open_eval + gap_extend_eval * (n as i32);
                    query_pos += n as usize;
                    op_index += n as usize;
                }
            }

            update_after_op!(index, op_index);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:612-616
    score /= factor;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:617-638
    if best_start_esp_index < gap_ops.len() && best_end_esp_index < gap_ops.len() {
        debug_assert!(matches!(gap_ops[best_start_esp_index], GapEditOp::Sub(_)));
        debug_assert!(matches!(gap_ops[best_end_esp_index], GapEditOp::Sub(_)));

        // NCBI's `qp > 0` and `qp < qlen` on the offsets in the context: an end outside the
        // context is not extended.
        let mut qp = best_q_start;
        let mut sp = best_s_start;
        let mut ext = 0usize;
        while qp > query.start && sp > 0 {
            let q_prev = qp - 1;
            let s_prev = sp - 1;
            if q_seq[q_prev] == s_seq[s_prev] && q_seq[q_prev] < 4 {
                qp -= 1;
                sp -= 1;
                ext += 1;
            } else {
                break;
            }
        }
        if ext > 0 {
            best_q_start -= ext;
            best_s_start -= ext;
            if let GapEditOp::Sub(n) = gap_ops[best_start_esp_index] {
                gap_ops[best_start_esp_index] = GapEditOp::Sub(n + ext as u32);
            }
            if best_end_esp_index == best_start_esp_index {
                best_end_esp_num += ext as i32;
            }
            score += (ext as i32) * reward;
        }

        let mut qp = best_q_end;
        let mut sp = best_s_end;
        let mut ext = 0usize;
        while qp < query.at(query.len) && sp < s_seq.len() && q_seq[qp] < 4 {
            if q_seq[qp] == s_seq[sp] {
                qp += 1;
                sp += 1;
                ext += 1;
            } else {
                break;
            }
        }
        if ext > 0 {
            best_q_end += ext;
            best_s_end += ext;
            if let GapEditOp::Sub(n) = gap_ops[best_end_esp_index] {
                gap_ops[best_end_esp_index] = GapEditOp::Sub(n + ext as u32);
            }
            best_end_esp_num += ext as i32;
            score += (ext as i32) * reward;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:640-646
    if score < cutoff_score {
        return true;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:414-470 s_UpdateReevaluatedHSP
    hit.raw_score = score;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:452-458
    // Back to NCBI's `Int4` offsets in the context.
    let best_q_start = (best_q_start as isize - query.start as isize) as i32;
    let best_q_end = (best_q_end as isize - query.start as isize) as i32;
    hit.internal_q_offset_0 = best_q_start;
    hit.internal_q_end_0 = best_q_end;
    hit.internal_s_offset_0 = best_s_start;
    hit.internal_s_end_0 = best_s_start.saturating_add(best_s_end.saturating_sub(best_s_start));
    let (q_start, q_end, s_start, s_end) = adjust_blastn_offsets(
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
        hit.query_length,
        hit.query_frame,
    );
    hit.q_start = q_start;
    hit.q_end = q_end;
    hit.s_start = s_start;
    hit.s_end = s_end;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:460-470
    let mut new_gap_info: Option<Vec<GapEditOp>> = None;
    if best_end_esp_index != gap_ops.len().saturating_sub(1) || best_start_esp_index > 0 {
        let mut subset = gap_ops[best_start_esp_index..=best_end_esp_index].to_vec();
        if let Some(last) = subset.last_mut() {
            let end_num = best_end_esp_num as u32;
            *last = match *last {
                GapEditOp::Sub(_) => GapEditOp::Sub(end_num),
                GapEditOp::Del(_) => GapEditOp::Del(end_num),
                GapEditOp::Ins(_) => GapEditOp::Ins(end_num),
            };
        }
        if !subset.is_empty() {
            new_gap_info = Some(subset);
        }
    } else if let Some(last) = gap_ops.last_mut() {
        let end_num = best_end_esp_num as u32;
        *last = match *last {
            GapEditOp::Sub(_) => GapEditOp::Sub(end_num),
            GapEditOp::Del(_) => GapEditOp::Del(end_num),
            GapEditOp::Ins(_) => GapEditOp::Ins(end_num),
        };
    }

    if let Some(new_gap_info) = new_gap_info {
        hit.gap_info = Some(new_gap_info);
    }

    // Recalculate bit_score and e_value if Karlin parameters are provided
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-250 Blast_HSPListGetEvalues, Blast_HSPListGetBitScores
    if let Some(params) = reeval_params {
        const LN2: f64 = std::f64::consts::LN_2; // NCBIMATH_LN2 (ncbi_math.h:161), the same f64

        // Bit score calculation
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1923-1928
        let log_k = params.karlin.k.glibc_ln();
        hit.bit_score = ((params.karlin.lambda * (score as f64)) - log_k) / LN2;

        // E-value calculation
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1864-1870
        let evalue_score = if params.round_down_evalue_score {
            score & !1
        } else {
            score
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
        if params.eff_searchsp > 0 {
            hit.e_value = crate::stats::karlin::evalue_from_raw_score_with_log_k(
                evalue_score,
                &params.karlin,
                log_k,
                params.eff_searchsp as f64,
            );
        }
    }

    false
}

#[cfg(test)]
#[path = "purge_endpoints_tests.rs"]
mod tests;
