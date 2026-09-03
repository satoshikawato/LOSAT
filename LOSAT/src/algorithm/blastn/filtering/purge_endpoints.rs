//! HSP endpoint purging for blastn - removes or trims HSPs with common start/end points
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535
//! Function: Blast_HSPListPurgeHSPsWithCommonEndpoints
//!
//! NCBI's endpoint purging has two modes:
//! - purge=FALSE: Trim overlapping HSPs using gap_info (edit script)
//! - purge=TRUE: Delete HSPs with common endpoints entirely
//!
//! The blastn traceback flow is:
//! 1. Call with purge=FALSE first (trim overlapping HSPs)
//! 2. Re-evaluate trimmed fragments
//! 3. Call with purge=TRUE for final cleanup
//!
//! Reference: blast_traceback.c:637-669

use super::super::hsp::{qsort_blastn_hsps_by, BlastnHsp};
use super::super::tracing as blastn_trace;
use crate::common::GapEditOp;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
// ```c
// if (hsp->query.frame != hsp->subject.frame) {
//    *q_end = query_length - hsp->query.offset;
//    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
//    *s_end = hsp->subject.offset + 1;
//    *s_start = hsp->subject.end;
// } else {
//    *q_start = hsp->query.offset + 1;
//    *q_end = hsp->query.end;
//    *s_start = hsp->subject.offset + 1;
//    *s_end = hsp->subject.end;
// }
// ```
fn adjust_blastn_offsets(
    query_offset: usize,
    query_end: usize,
    subject_offset: usize,
    subject_end: usize,
    query_length: usize,
    query_frame: i32,
) -> (usize, usize, usize, usize) {
    if query_frame < 0 {
        let q_end = query_length.saturating_sub(query_offset);
        let q_start = query_length.saturating_sub(query_end).saturating_add(1);
        let s_end = subject_offset + 1;
        let s_start = subject_end;
        (q_start, q_end, s_start, s_end)
    } else {
        let q_start = query_offset + 1;
        let q_end = query_end;
        let s_start = subject_offset + 1;
        let s_end = subject_end;
        (q_start, q_end, s_start, s_end)
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2484-2487
// ```c
// hsp_array[i]->context == hsp_array[i+j]->context &&
// hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
// hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
// hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame
// ```
fn blastn_hsp_context(hsp: &BlastnHsp) -> u32 {
    hsp.q_idx * 2 + if hsp.query_frame < 0 { 1 } else { 0 }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2285-2318
// ```c
// if (h1->context < h2->context) return -1;
// if (h1->context > h2->context) return 1;
// if (h1->query.offset < h2->query.offset) return -1;
// if (h1->query.offset > h2->query.offset) return 1;
// if (h1->subject.offset < h2->subject.offset) return -1;
// if (h1->subject.offset > h2->subject.offset) return 1;
// if (h1->score < h2->score) return 1;
// if (h1->score > h2->score) return -1;
// if (h1->query.end < h2->query.end) return 1;
// if (h1->query.end > h2->query.end) return -1;
// if (h1->subject.end < h2->subject.end) return 1;
// if (h1->subject.end > h2->subject.end) return -1;
// ```
fn compare_query_offset_hsps_for_common_endpoint(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    blastn_hsp_context(a)
        .cmp(&blastn_hsp_context(b))
        .then_with(|| a.internal_q_offset_0.cmp(&b.internal_q_offset_0))
        .then_with(|| a.internal_s_offset_0.cmp(&b.internal_s_offset_0))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.internal_q_end_0.cmp(&a.internal_q_end_0))
        .then_with(|| b.internal_s_end_0.cmp(&a.internal_s_end_0))
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2350-2384
// ```c
// if (h1->context < h2->context) return -1;
// if (h1->context > h2->context) return 1;
// if (h1->query.end < h2->query.end) return -1;
// if (h1->query.end > h2->query.end) return 1;
// if (h1->subject.end < h2->subject.end) return -1;
// if (h1->subject.end > h2->subject.end) return 1;
// if (h1->score < h2->score) return 1;
// if (h1->score > h2->score) return -1;
// if (h1->query.offset < h2->query.offset) return 1;
// if (h1->query.offset > h2->query.offset) return -1;
// if (h1->subject.offset < h2->subject.offset) return 1;
// if (h1->subject.offset > h2->subject.offset) return -1;
// ```
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
// Reference: blast_hits.c:2392-2452 s_CutOffGapEditScript
// =============================================================================

/// Cut off gap edit script at specified position.
///
/// NCBI reference: blast_hits.c:2392-2452 s_CutOffGapEditScript
///
/// This function trims an HSP's edit script at a specified cut point.
/// - If cut_begin=true: Keep the END portion, update offsets
/// - If cut_begin=false: Keep the START portion, update ends
///
/// ```c
/// // NCBI s_CutOffGapEditScript (blast_hits.c:2392-2452)
/// static void
/// s_CutOffGapEditScript(BlastHSP* hsp, Int4 q_cut, Int4 s_cut, Boolean cut_begin)
/// {
///    int index, opid, qid, sid;
///    Boolean found = FALSE;
///    GapEditScript *esp = hsp->gap_info;
///    qid = 0;
///    sid = 0;
///    opid = 0;
///    q_cut -= hsp->query.offset;
///    s_cut -= hsp->subject.offset;
///    for (index=0; index < esp->size; index++) {
///        for(opid=0; opid < esp->num[index];){
///           if (esp->op_type[index] == eGapAlignSub) {
///              qid++;
///              sid++;
///              opid++;
///           } else if (esp->op_type[index] == eGapAlignDel) {
///              sid+=esp->num[index];
///              opid+=esp->num[index];
///           } else if (esp->op_type[index] == eGapAlignIns) {
///              qid+=esp->num[index];
///              opid+=esp->num[index];
///           }
///           if (qid >= q_cut && sid >= s_cut) found = TRUE;
///           if (found) break;
///        }
///        if (found) break;
///    }
///    // ... (trimming logic follows)
/// }
/// ```
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
    q_cut: usize,
    s_cut: usize,
    cut_begin: bool,
) -> bool {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2392-2396
    // ```c
    // s_CutOffGapEditScript(BlastHSP* hsp, Int4 q_cut, Int4 s_cut, Boolean cut_begin)
    // {
    //    ...
    //    GapEditScript *esp = hsp->gap_info;
    // ```
    // NCBI mutates the HSP's existing edit script in place. Move the Vec out
    // while trimming so Rust does not clone the whole script in this hot path.
    let mut gap_info = match hit.gap_info.take() {
        Some(info) if !info.is_empty() => info,
        other => {
            hit.gap_info = other;
            return false;
        }
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2398-2401
    // ```c
    // GapEditScript *esp = hsp->gap_info;
    // q_cut -= hsp->query.offset;
    // s_cut -= hsp->subject.offset;
    // ```
    // Endpoint trimming operates on canonical internal offsets already stored
    // in the HSP. Do not reconstruct them from output coordinates.
    let mut q_offset_0 = hit.internal_q_offset_0;
    let mut q_end_0 = hit.internal_q_end_0;
    let mut s_offset_0 = hit.internal_s_offset_0;
    let mut s_end_0 = hit.internal_s_end_0;

    // Convert absolute cut positions to relative (within HSP)
    // NCBI: q_cut -= hsp->query.offset; s_cut -= hsp->subject.offset;
    let q_cut_rel = q_cut.saturating_sub(q_offset_0);
    let s_cut_rel = s_cut.saturating_sub(s_offset_0);

    // Walk through edit script to find cut point
    // NCBI tracks: qid (query consumed), sid (subject consumed), opid (ops within current run)
    let mut qid: usize = 0;
    let mut sid: usize = 0;
    let mut found = false;
    let mut found_index: usize = 0;
    let mut found_opid: usize = 0;

    // NCBI: for (index=0; index < esp->size; index++)
    'outer: for (index, op) in gap_info.iter().enumerate() {
        let num = op.num() as usize;

        // NCBI: for(opid=0; opid < esp->num[index];)
        match op {
            GapEditOp::Sub(_) => {
                // Sub: iterate position by position to find exact cut point
                // NCBI: qid++; sid++; opid++;
                for opid in 0..num {
                    qid += 1;
                    sid += 1;

                    // NCBI: if (qid >= q_cut && sid >= s_cut) found = TRUE;
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
                // NCBI: sid+=esp->num[index]; opid+=esp->num[index];
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
                // NCBI: qid+=esp->num[index]; opid+=esp->num[index];
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
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2426-2439
        // ```c
        // if (cut_begin) {
        //     int new_index = 0;
        //     if (opid < esp->num[index]) {
        //        esp->op_type[0] = esp->op_type[index];
        //        esp->num[0] = esp->num[index] - opid;
        //        new_index++;
        //     }
        //     ++index;
        //     for (; index < esp->size; index++, new_index++) {
        //        esp->op_type[new_index] = esp->op_type[index];
        //        esp->num[new_index] = esp->num[index];
        //     }
        //     esp->size = new_index;
        // ```
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
        // NCBI: hsp->query.offset += qid; hsp->subject.offset += sid;
        q_offset_0 = q_offset_0.saturating_add(qid);
        s_offset_0 = s_offset_0.saturating_add(sid);
    } else {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2442-2449
        // ```c
        // } else {
        //     if (opid < esp->num[index]) {
        //        ASSERT(esp->op_type[index] == eGapAlignSub);
        //        esp->num[index] = opid;
        //     }
        //     esp->size = index+1;
        //     hsp->query.end = hsp->query.offset + qid;
        //     hsp->subject.end = hsp->subject.offset + sid;
        // }
        // ```
        // Keep the prefix and truncate the same edit script buffer.
        if let GapEditOp::Sub(n) = gap_info[found_index] {
            if found_opid < n as usize {
                gap_info[found_index] = GapEditOp::Sub(found_opid as u32);
            }
        }
        gap_info.truncate(found_index.saturating_add(1));

        // Update HSP coordinates
        // NCBI: hsp->query.end = hsp->query.offset + qid;
        // hsp->subject.end = hsp->subject.offset + sid;
        q_end_0 = q_offset_0.saturating_add(qid);
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2436-2444
    // ```c
    // hsp->query.offset += qid;
    // hsp->subject.offset += sid;
    // ...
    // hsp->query.end = hsp->query.offset + qid;
    // hsp->subject.end = hsp->subject.offset + sid;
    // ```
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
/// NCBI reference: blast_hits.c:2455-2535 Blast_HSPListPurgeHSPsWithCommonEndpoints
///
/// CRITICAL: NCBI processes HSPs per-subject (BlastHSPList is per-subject).
/// This function groups hits by subject index (oid) before purging to match NCBI behavior.
///
/// # Parameters
/// * `purge` - If true, delete HSPs entirely (current behavior).
///             If false, trim overlapping HSPs using gap_info and move to end.
///
/// ```c
/// // NCBI s_QueryOffsetCompareHSPs (blast_hits.c:2268-2320)
/// // Sort by: context, query.offset, subject.offset, score DESC
/// if (h1->context < h2->context) return -1;
/// if (h1->context > h2->context) return 1;
/// if (h1->query.offset < h2->query.offset) return -1;
/// if (h1->query.offset > h2->query.offset) return 1;
/// if (h1->subject.offset < h2->subject.offset) return -1;
/// if (h1->subject.offset > h2->subject.offset) return 1;
/// if (h1->score < h2->score) return 1;  // DESC
/// if (h1->score > h2->score) return -1;
/// ```
///
/// ```c
/// // NCBI purge condition (blast_hits.c:2482-2487)
/// hsp_array[i]->context == hsp_array[i+j]->context &&
/// hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
/// hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
/// hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame
/// ```
///
/// For blastn:
/// - context = query index + strand (q_idx * 2 + {0,1})
/// - subject.frame is always +1 (blastn subject only), so strand is encoded by context.
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:637-668
    // ```c
    // Int4 extra_start =
    //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
    // ...
    // if(program_number == eBlastTypeBlastn) {
    //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // }
    // ```
    // Traceback calls this on one BlastHSPList, i.e. one subject, in the hot path.
    // Avoid rebuilding the same per-subject grouping when the input already matches
    // NCBI's per-subject HSP list shape.
    if hits.iter().all(|h| h.s_idx == hits[0].s_idx) {
        return purge_hsps_for_subject_ex(hits, purge);
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
    //    BlastHSP** hsp_array;
    //    Int4 hspcnt;
    //    ...
    // } BlastHSPList;
    // ```
    // NCBI reference: blast_hits.c:2455-2535
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

/// Purge HSPs with common endpoints for a single subject.
///
/// NCBI reference: blast_hits.c:2455-2535
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
    let mut hsp_array: Vec<Option<BlastnHsp>> = hits.into_iter().map(Some).collect();
    let mut hsp_count = initial_count;
    let mut start_removed = 0usize;
    let mut end_removed = 0usize;
    let mut start_trimmed = 0usize;
    let mut end_trimmed = 0usize;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2471-2478
    // ```c
    // hsp_array = hsp_list->hsp_array;
    // hsp_count = hsp_list->hspcnt;
    //
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // ```
    qsort_hsp_array_prefix_by(
        &mut hsp_array,
        hsp_count,
        compare_query_offset_hsps_for_common_endpoint,
    );

    // NCBI reference: blast_hits.c:2480-2500
    // if (!purge && (hsp->query.end > hsp_array[i]->query.end)) {
    //     s_CutOffGapEditScript(hsp, hsp_array[i]->query.end, hsp_array[i]->subject.end, TRUE);
    // } else {
    //     hsp = Blast_HSPFree(hsp);
    // }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2500
    // ```c
    // while (i+j < hsp_count &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
    //        hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
    //        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
    //    hsp_count--;
    //    hsp = hsp_array[i+j];
    //    ...
    //    for (k=i+j; k<hsp_count; k++) hsp_array[k] = hsp_array[k+1];
    //    hsp_array[hsp_count] = hsp;
    // }
    // ```
    // The C code keeps the sorted active prefix and moves trimmed/null HSPs to
    // the tail by decrementing hsp_count before shifting hsp_array left.
    let mut i = 0usize;
    while i < hsp_count {
        let j = 1usize;
        while i + j < hsp_count && same_common_start(&hsp_array, i, i + j) {
            let (keeper_q_end, keeper_s_end) = {
                let keeper = hsp_array[i]
                    .as_ref()
                    .expect("same_common_start requires a keeper HSP");
                (keeper.internal_q_end_0, keeper.internal_s_end_0)
            };
            hsp_count -= 1;
            let removed = hsp_array[i + j].take();
            let mut tail_hsp = None;

            if let Some(mut removed_hit) = removed {
                let keeper = hsp_array[i]
                    .as_ref()
                    .expect("same_common_start requires a keeper HSP");
                if !purge && removed_hit.internal_q_end_0 > keeper_q_end {
                    trace_common_endpoint_purge_decision(
                        "common_start",
                        "trim_begin",
                        purge,
                        keeper,
                        &removed_hit,
                        keeper_q_end,
                        keeper_s_end,
                    );
                    let _ =
                        cut_off_gap_edit_script(&mut removed_hit, keeper_q_end, keeper_s_end, true);
                    tail_hsp = Some(removed_hit);
                    start_trimmed += 1;
                } else {
                    trace_common_endpoint_purge_decision(
                        "common_start",
                        "delete",
                        purge,
                        keeper,
                        &removed_hit,
                        keeper_q_end,
                        keeper_s_end,
                    );
                }
                start_removed += 1;
            }

            for k in (i + j)..hsp_count {
                hsp_array[k] = hsp_array[k + 1].take();
            }
            hsp_array[hsp_count] = tail_hsp;
        }
        i += j;
    }

    // Pass 2: Remove HSPs with common END positions
    // NCBI reference: blast_hits.c:2504-2526
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2504
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
    // ```
    qsort_hsp_array_prefix_by(
        &mut hsp_array,
        hsp_count,
        compare_query_end_hsps_for_common_endpoint,
    );

    // NCBI reference: blast_hits.c:2516-2520
    // if (!purge && (hsp->query.offset < hsp_array[i]->query.offset)) {
    //     s_CutOffGapEditScript(hsp, hsp_array[i]->query.offset, hsp_array[i]->subject.offset, FALSE);
    // } else {
    //     hsp = Blast_HSPFree(hsp);
    // }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2512-2528
    // ```c
    // while (i+j < hsp_count &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.end == hsp_array[i+j]->query.end &&
    //        hsp_array[i]->subject.end == hsp_array[i+j]->subject.end &&
    //        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
    //    hsp_count--;
    //    hsp = hsp_array[i+j];
    //    ...
    //    for (k=i+j; k<hsp_count; k++) hsp_array[k] = hsp_array[k+1];
    //    hsp_array[hsp_count] = hsp;
    // }
    // ```
    // Repeat the same active-prefix mutation NCBI uses for common endpoints.
    let mut i = 0usize;
    while i < hsp_count {
        let j = 1usize;
        while i + j < hsp_count && same_common_end(&hsp_array, i, i + j) {
            let (keeper_q_offset, keeper_s_offset) = {
                let keeper = hsp_array[i]
                    .as_ref()
                    .expect("same_common_end requires a keeper HSP");
                (keeper.internal_q_offset_0, keeper.internal_s_offset_0)
            };
            hsp_count -= 1;
            let removed = hsp_array[i + j].take();
            let mut tail_hsp = None;

            if let Some(mut removed_hit) = removed {
                let keeper = hsp_array[i]
                    .as_ref()
                    .expect("same_common_end requires a keeper HSP");
                if !purge && removed_hit.internal_q_offset_0 < keeper_q_offset {
                    trace_common_endpoint_purge_decision(
                        "common_end",
                        "trim_end",
                        purge,
                        keeper,
                        &removed_hit,
                        keeper_q_offset,
                        keeper_s_offset,
                    );
                    let _ = cut_off_gap_edit_script(
                        &mut removed_hit,
                        keeper_q_offset,
                        keeper_s_offset,
                        false,
                    );
                    tail_hsp = Some(removed_hit);
                    end_trimmed += 1;
                } else {
                    trace_common_endpoint_purge_decision(
                        "common_end",
                        "delete",
                        purge,
                        keeper,
                        &removed_hit,
                        keeper_q_offset,
                        keeper_s_offset,
                    );
                }
                end_removed += 1;
            }

            for k in (i + j)..hsp_count {
                hsp_array[k] = hsp_array[k + 1].take();
            }
            hsp_array[hsp_count] = tail_hsp;
        }
        i += j;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2530-2535
    // ```c
    // if (purge) {
    //    Blast_HSPListPurgeNullHSPs(hsp_list);
    // }
    //
    // return hsp_count;
    // ```
    // In purge=false mode NCBI leaves trimmed HSPs in the tail for
    // blast_traceback.c:647-665 to re-evaluate, while NULL tail entries are
    // skipped. Store the same live order and return the active-prefix boundary.
    let extra_start = hsp_count;
    let mut hits = Vec::with_capacity(initial_count);
    for slot in hsp_array {
        if let Some(hit) = slot {
            hits.push(hit);
        }
    }

    // Debug output for purge statistics
    if std::env::var("LOSAT_DEBUG_BLASTN").is_ok() {
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2387
// ```c
// if (!h1 && !h2) return 0;
// else if (!h1) return 1;
// else if (!h2) return -1;
// ```
// NCBI's qsort callbacks sort NULL HSP pointers after real HSPs. The active
// prefix normally contains only real HSPs, but preserve the NULL ordering when
// mirroring Blast_HSPListPurgeHSPsWithCommonEndpoints exactly.
fn qsort_hsp_array_prefix_by(
    hsp_array: &mut [Option<BlastnHsp>],
    hsp_count: usize,
    compare: fn(&BlastnHsp, &BlastnHsp) -> Ordering,
) {
    if hsp_count <= 1 {
        return;
    }
    let prefix_len = hsp_count.min(hsp_array.len());
    let mut active = Vec::with_capacity(prefix_len);
    let mut null_count = 0usize;
    for slot in &mut hsp_array[..prefix_len] {
        if let Some(hit) = slot.take() {
            active.push(hit);
        } else {
            null_count += 1;
        }
    }

    qsort_blastn_hsps_by(&mut active, compare);

    let mut index = 0usize;
    for hit in active {
        hsp_array[index] = Some(hit);
        index += 1;
    }
    for slot in &mut hsp_array[index..index + null_count] {
        *slot = None;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2482-2487
// ```c
// hsp_array[i] && hsp_array[i+j] &&
// hsp_array[i]->context == hsp_array[i+j]->context &&
// hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
// hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
// hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame
// ```
fn same_common_start(hsp_array: &[Option<BlastnHsp>], lhs: usize, rhs: usize) -> bool {
    let (Some(a), Some(b)) = (
        hsp_array.get(lhs).and_then(Option::as_ref),
        hsp_array.get(rhs).and_then(Option::as_ref),
    ) else {
        return false;
    };
    blastn_hsp_context(a) == blastn_hsp_context(b)
        && a.internal_q_offset_0 == b.internal_q_offset_0
        && a.internal_s_offset_0 == b.internal_s_offset_0
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2508-2513
// ```c
// hsp_array[i] && hsp_array[i+j] &&
// hsp_array[i]->context == hsp_array[i+j]->context &&
// hsp_array[i]->query.end == hsp_array[i+j]->query.end &&
// hsp_array[i]->subject.end == hsp_array[i+j]->subject.end &&
// hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame
// ```
fn same_common_end(hsp_array: &[Option<BlastnHsp>], lhs: usize, rhs: usize) -> bool {
    let (Some(a), Some(b)) = (
        hsp_array.get(lhs).and_then(Option::as_ref),
        hsp_array.get(rhs).and_then(Option::as_ref),
    ) else {
        return false;
    };
    blastn_hsp_context(a) == blastn_hsp_context(b)
        && a.internal_q_end_0 == b.internal_q_end_0
        && a.internal_s_end_0 == b.internal_s_end_0
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2480-2499
// ```c
// while (i+j < hsp_count &&
//        hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
//        hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
//        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
//    hsp_count--;
//    hsp = hsp_array[i+j];
//    if (!purge && (hsp->query.end > hsp_array[i]->query.end)) {
//        s_CutOffGapEditScript(hsp, hsp_array[i]->query.end,
//                                  hsp_array[i]->subject.end, TRUE);
//    } else {
//        hsp = Blast_HSPFree(hsp);
//    }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2512-2525
// ```c
// while (i+j < hsp_count &&
//        hsp_array[i]->query.end == hsp_array[i+j]->query.end &&
//        hsp_array[i]->subject.end == hsp_array[i+j]->subject.end &&
//        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame) {
//    hsp_count--;
//    hsp = hsp_array[i+j];
//    if (!purge && (hsp->query.offset < hsp_array[i]->query.offset)) {
//        s_CutOffGapEditScript(hsp, hsp_array[i]->query.offset,
//                                  hsp_array[i]->subject.offset, FALSE);
//    } else {
//        hsp = Blast_HSPFree(hsp);
//    }
// ```
fn trace_common_endpoint_purge_decision(
    pass: &str,
    action: &str,
    purge: bool,
    keeper: &BlastnHsp,
    removed: &BlastnHsp,
    cut_q: usize,
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
// Reference: blast_hits.c:479-647 Blast_HSPReevaluateWithAmbiguitiesGapped
// =============================================================================

/// BLASTNA alphabet size.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1060-1068 (BLASTNA_SIZE usage)
const BLASTNA_SIZE: usize = 16;

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:45-53
    // ```c
    // typedef enum {
    //    eGapAlignDel = 0, /**< Deletion: a gap in query */
    //    eGapAlignSub = 3, /**< Substitution */
    //    eGapAlignIns = 6, /**< Insertion: a gap in subject */
    // } EGapAlignOpType;
    // ```
    fn edit_script_consumption(edit_ops: &[GapEditOp]) -> (usize, usize) {
        edit_ops.iter().fold(
            (0usize, 0usize),
            |(query_consumed, subject_consumed), op| match op {
                GapEditOp::Sub(count) => (
                    query_consumed + *count as usize,
                    subject_consumed + *count as usize,
                ),
                GapEditOp::Del(count) => (query_consumed, subject_consumed + *count as usize),
                GapEditOp::Ins(count) => (query_consumed + *count as usize, subject_consumed),
            },
        )
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    Int4 num_ident;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    //    GapEditScript* gap_info;
    // } BlastHSP;
    // ```
    fn opaque_equal_endpoint_hsps() -> (BlastnHsp, BlastnHsp) {
        let g_ops = vec![
            GapEditOp::Sub(10),
            GapEditOp::Ins(1),
            GapEditOp::Sub(1),
            GapEditOp::Del(1),
            GapEditOp::Sub(10),
        ];
        let s_ops = vec![GapEditOp::Sub(22)];
        let opaque_query = [1u8; 22];
        let opaque_subject = [1u8; 22];
        let g_stats = stats_from_edit_ops(&opaque_query, &opaque_subject, 0, 0, &g_ops);
        let s_stats = stats_from_edit_ops(&opaque_query, &opaque_subject, 0, 0, &s_ops);

        let make_hsp = |gap_info: Vec<GapEditOp>,
                        (matches, mismatches, gap_opens, _gap_letters, align_length): (
            usize,
            usize,
            usize,
            usize,
            usize,
        )| BlastnHsp {
            identity: matches as f64 * 100.0 / align_length as f64,
            length: align_length,
            mismatch: mismatches,
            gapopen: gap_opens,
            q_start: 6,
            q_end: 27,
            s_start: 12,
            s_end: 33,
            e_value: 1.0e-9,
            bit_score: 31.0,
            num_ident: matches,
            query_frame: 1,
            query_length: 64,
            q_idx: 2,
            s_idx: 3,
            raw_score: 42,
            internal_q_offset_0: 5,
            internal_q_end_0: 27,
            internal_s_offset_0: 11,
            internal_s_end_0: 33,
            internal_query_context_offset: 0,
            gap_info: Some(gap_info),
            num_positives: matches,
        };

        (make_hsp(g_ops, g_stats), make_hsp(s_ops, s_stats))
    }

    // Test-only post-sort injection seam for the NCBI first-survivor loop.
    // It intentionally omits qsort so tests can supply either permitted order
    // without depending on the host libc's ordering of comparator-equal items.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2500
    // ```c
    // qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
    // while (i+j < hsp_count &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
    //        hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset) {
    //    hsp_count--;
    //    hsp = hsp_array[i+j];
    //    hsp = Blast_HSPFree(hsp);
    //    for (k=i+j; k<hsp_count; k++) hsp_array[k] = hsp_array[k+1];
    // }
    // ```
    fn purge_after_injected_common_start_sort(post_sort_order: Vec<BlastnHsp>) -> Vec<BlastnHsp> {
        let mut hsp_array: Vec<Option<BlastnHsp>> = post_sort_order.into_iter().map(Some).collect();
        let mut hsp_count = hsp_array.len();
        let mut i = 0usize;
        while i < hsp_count {
            let j = 1usize;
            while i + j < hsp_count && same_common_start(&hsp_array, i, i + j) {
                hsp_count -= 1;
                let _removed = hsp_array[i + j].take();
                for k in (i + j)..hsp_count {
                    hsp_array[k] = hsp_array[k + 1].take();
                }
                hsp_array[hsp_count] = None;
            }
            i += j;
        }
        hsp_array.into_iter().take(hsp_count).flatten().collect()
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2387
    // Both common-endpoint comparators return zero after context, endpoints,
    // and score tie; neither comparator examines GapEditScript.
    #[test]
    fn equivalent_scripts_keep_endpoints_and_score() {
        let (g_path, s_path) = opaque_equal_endpoint_hsps();
        let g_ops = g_path
            .gap_info
            .as_deref()
            .expect("G path has an edit script");
        let s_ops = s_path
            .gap_info
            .as_deref()
            .expect("S path has an edit script");

        assert_eq!(edit_script_consumption(g_ops), (22, 22));
        assert_eq!(edit_script_consumption(s_ops), (22, 22));
        assert_eq!(g_path.internal_q_offset_0, s_path.internal_q_offset_0);
        assert_eq!(g_path.internal_q_end_0, s_path.internal_q_end_0);
        assert_eq!(g_path.internal_s_offset_0, s_path.internal_s_offset_0);
        assert_eq!(g_path.internal_s_end_0, s_path.internal_s_end_0);
        assert_eq!(g_path.raw_score, s_path.raw_score);
        assert_ne!(g_ops, s_ops);
        assert_eq!(
            (
                g_path.num_ident,
                g_path.length,
                g_path.mismatch,
                g_path.gapopen
            ),
            (21, 23, 0, 2)
        );
        assert_eq!(
            (
                s_path.num_ident,
                s_path.length,
                s_path.mismatch,
                s_path.gapopen
            ),
            (22, 22, 0, 0)
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2285-2387
    // ```c
    // if (h1->score < h2->score) return 1;
    // if (h1->score > h2->score) return -1;
    // ... compare the remaining endpoint ...
    // return 0;
    // ```
    #[test]
    fn common_endpoint_comparator_ignores_edit_script() {
        let (g_path, s_path) = opaque_equal_endpoint_hsps();

        assert_ne!(g_path.gap_info, s_path.gap_info);
        assert_eq!(
            compare_query_offset_hsps_for_common_endpoint(&g_path, &s_path),
            Ordering::Equal
        );
        assert_eq!(
            compare_query_offset_hsps_for_common_endpoint(&s_path, &g_path),
            Ordering::Equal
        );
        assert_eq!(
            compare_query_end_hsps_for_common_endpoint(&g_path, &s_path),
            Ordering::Equal
        );
        assert_eq!(
            compare_query_end_hsps_for_common_endpoint(&s_path, &g_path),
            Ordering::Equal
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478-2500
    // The common-start purge keeps hsp_array[i] and frees hsp_array[i+j].
    #[test]
    fn source_compatible_first_survivor_accepts_either_equal_order() {
        let (g_path, s_path) = opaque_equal_endpoint_hsps();
        let g_script = g_path.gap_info.clone();
        let s_script = s_path.gap_info.clone();

        let g_then_s = purge_after_injected_common_start_sort(vec![g_path.clone(), s_path.clone()]);
        let s_then_g = purge_after_injected_common_start_sort(vec![s_path, g_path]);

        assert_eq!(g_then_s.len(), 1);
        assert_eq!(s_then_g.len(), 1);
        assert_eq!(g_then_s[0].gap_info, g_script);
        assert_eq!(s_then_g[0].gap_info, s_script);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2488-2498
    // In purge mode the non-survivor is freed; the keeper's GapEditScript is
    // neither trimmed nor combined with the removed HSP.
    #[test]
    fn source_compatible_purge_never_synthesizes_edit_script() {
        let (g_path, s_path) = opaque_equal_endpoint_hsps();
        let exact_inputs = [
            (
                g_path.gap_info.clone(),
                g_path.num_ident,
                g_path.length,
                g_path.mismatch,
                g_path.gapopen,
            ),
            (
                s_path.gap_info.clone(),
                s_path.num_ident,
                s_path.length,
                s_path.mismatch,
                s_path.gapopen,
            ),
        ];

        let survivors =
            purge_after_injected_common_start_sort(vec![g_path.clone(), s_path.clone()])
                .into_iter()
                .chain(purge_after_injected_common_start_sort(vec![s_path, g_path]))
                .collect::<Vec<_>>();
        assert_eq!(survivors.len(), 2);

        for survivor in survivors {
            let observed = (
                survivor.gap_info,
                survivor.num_ident,
                survivor.length,
                survivor.mismatch,
                survivor.gapopen,
            );
            assert!(exact_inputs.contains(&observed));
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2285-2318
    // ```c
    // if (h1->query.offset < h2->query.offset) return -1;
    // if (h1->query.offset > h2->query.offset) return 1;
    // if (h1->subject.offset < h2->subject.offset) return -1;
    // if (h1->subject.offset > h2->subject.offset) return 1;
    // if (h1->score < h2->score) return 1;
    // if (h1->score > h2->score) return -1;
    // if (h1->query.end < h2->query.end) return 1;
    // if (h1->query.end > h2->query.end) return -1;
    // if (h1->subject.end < h2->subject.end) return 1;
    // if (h1->subject.end > h2->subject.end) return -1;
    // ```
    #[test]
    fn test_purge_common_start_keeps_longer_end_on_score_tie() {
        fn make_hit(q_end: usize, s_end: usize) -> BlastnHsp {
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
            BlastnHsp {
                identity: 100.0,
                length: q_end,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end,
                s_start: 1,
                s_end,
                e_value: 0.0,
                bit_score: 0.0,
                num_ident: q_end,
                query_frame: 1,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
                raw_score: 100,
                internal_q_offset_0: 0,
                internal_q_end_0: q_end,
                internal_s_offset_0: 0,
                internal_s_end_0: s_end,
                internal_query_context_offset: 0,
                gap_info: None,
                num_positives: q_end,
            }
        }

        let short = make_hit(50, 50);
        let long = make_hit(80, 80);

        let (purged, _) = purge_hsps_with_common_endpoints_ex(vec![short, long], true);
        assert_eq!(purged.len(), 1);
        assert_eq!(purged[0].q_end, 80);
        assert_eq!(purged[0].s_end, 80);
    }
}

/// Compute alignment stats from a gap edit script.
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:745-818 (identity counting)
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1055-1074 (gap counts)
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

/// Test if an HSP should be deleted based on percent identity and minimum hit length.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
/// ```c
/// return ((hsp->num_ident * 100.0 <
///         align_length * hit_options->percent_identity) ||
///         align_length < hit_options->min_hit_length) ;
/// ```
pub fn hsp_test(
    num_ident: usize,
    align_length: usize,
    percent_identity: f64,
    min_hit_length: usize,
) -> bool {
    let identity_check = if percent_identity > 0.0 {
        (num_ident as f64 * 100.0) < (align_length as f64 * percent_identity)
    } else {
        false
    };

    let length_check = if min_hit_length > 0 {
        align_length < min_hit_length
    } else {
        false
    };

    identity_check || length_check
}

/// Update identity/length stats and apply Blast_HSPTestIdentityAndLength filtering.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:656-663
/// ```c
/// delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
/// if (!delete_hsp)
///     delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp,
///                                                 query_nomask, subject,
///                                                 score_options, hit_options);
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:745-818
/// ```c
/// for (index=0; index<esp->size; index++) {
///     align_length += esp->num[index];
///     if (esp->op_type[index] == eGapAlignSub) {
///         if (*q == *s) num_ident++;
///         ...
///     }
/// }
/// ```
pub fn blast_hsp_test_identity_and_length(
    hit: &mut BlastnHsp,
    q_seq: &[u8],
    s_seq: &[u8],
    percent_identity: f64,
    min_hit_length: usize,
) -> bool {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:758-765
    // ```c
    // q_off = hsp->query.offset;
    // s_off = hsp->subject.offset;
    // q = (Uint1*) &query[q_off];
    // s = (Uint1*) &subject[s_off];
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:678-692
    // ```c
    // Since the list is sorted by score already, any HSP
    // contained by a previous HSP is guaranteed to have a
    // lower score, and may be purged.
    // ```
    // BLASTN re-evaluation and identity counting operate on canonical internal
    // offsets; output coordinates are adjusted only after pruning.
    let q_offset = hit.internal_q_offset_0;
    let s_offset = hit.internal_s_offset_0;

    let (matches, mismatches, gap_opens, _gap_letters, align_length) =
        if let Some(ref ops) = hit.gap_info {
            stats_from_edit_ops(q_seq, s_seq, q_offset, s_offset, ops)
        } else {
            // Ungapped fallback (N.B. not expected for BLASTN reevaluation path).
            let align_length = hit.internal_q_end_0.saturating_sub(hit.internal_q_offset_0);
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
    // ```c
    // s_Blast_HSPGetNumIdentitiesAndPositives(..., Int4* num_ident_ptr, ...,
    //                                         Int4* num_pos_ptr)
    // {
    //    ...
    //    *num_ident_ptr = num_ident;
    //    ...
    //    *num_pos_ptr = num_pos + num_ident;
    // }
    // ```
    hit.num_ident = matches;
    hit.num_positives = matches;
    hit.identity = if align_length > 0 {
        (matches as f64 / align_length as f64) * 100.0
    } else {
        0.0
    };

    hsp_test(matches, align_length, percent_identity, min_hit_length)
}

/// Parameters needed for E-value and bit score recalculation after re-evaluation.
/// NCBI reference: blast_traceback.c:234-250 s_HSPListPostTracebackUpdate
pub struct ReevalParams {
    pub lambda: f64,
    pub k: f64,
    pub eff_searchsp: i64,
    pub db_len: usize,
    pub db_num_seqs: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1864-1870
    // ```c
    // score = hsp->score;
    // if (hsp_list && hsp_list->hspcnt != 0
    //         && gapped_calculation && sbp->round_down) {
    //     score &= ~1;
    // }
    // ```
    pub round_down_evalue_score: bool,
}

/// Re-evaluate a gapped HSP after trimming, finding the best-scoring sub-alignment.
///
/// NCBI reference: blast_hits.c:479-647 Blast_HSPReevaluateWithAmbiguitiesGapped
///
/// This function walks through the edit script to find the best-scoring
/// sub-alignment, using the BLASTNA scoring matrix for ambiguous bases.
///
/// ```c
/// // NCBI Blast_HSPReevaluateWithAmbiguitiesGapped (blast_hits.c:479-647)
/// Boolean Blast_HSPReevaluateWithAmbiguitiesGapped(BlastHSP* hsp,
///            const Uint1* q, const Int4 qlen,
///            const Uint1* s, const Int4 slen,
///            const BlastHitSavingParameters* hit_params,
///            const BlastScoringParameters* score_params,
///            const BlastScoreBlk* sbp)
/// {
///    // Walk through edit script, find best scoring sub-alignment
///    esp = hsp->gap_info;
///    if (!esp) return TRUE;  // No gap_info = delete HSP
///    for (index=0; index<esp->size; index++) {
///        for (op_index=0; op_index<esp->num[index]; ) {
///           if (esp->op_type[index] == eGapAlignSub) {
///               sum += factor*matrix[*query & kResidueMask][*subject];
///               query++; subject++; op_index++;
///           } else if (esp->op_type[index] == eGapAlignDel) {
///               sum -= gap_open + gap_extend * esp->num[index];
///               subject += esp->num[index]; op_index += esp->num[index];
///           } else if (esp->op_type[index] == eGapAlignIns) {
///               sum -= gap_open + gap_extend * esp->num[index];
///               query += esp->num[index]; op_index += esp->num[index];
///           }
///           // Track best scoring region
///           if (sum < 0) { reset to current position }
///           else if (sum > score) { update best region }
///        }
///    }
///    // Update HSP with best region found
///    return s_UpdateReevaluatedHSP(hsp, ...);
/// }
/// ```
///
/// # Parameters
/// * `hit` - The HSP to re-evaluate (modified in place if best region found)
/// * `q_seq` - Query sequence in BLASTNA encoding (full, 0-indexed)
/// * `s_seq` - Subject sequence in BLASTNA encoding (full, 0-indexed, aligned orientation)
/// * `reward` - Match reward (positive)
/// * `penalty` - Mismatch penalty (negative)
/// * `gap_open` - Gap open penalty (positive)
/// * `gap_extend` - Gap extend penalty (positive)
/// * `cutoff_score` - Minimum score to keep HSP
/// * `score_matrix` - Prebuilt BLASTNA scoring matrix (sbp->matrix->data)
/// * `reeval_params` - Optional parameters for E-value/bit score recalculation
///
/// # Returns
/// * `true` if HSP should be DELETED (score < cutoff or no valid region found)
/// * `false` if HSP is OK (score >= cutoff, hit is updated in place)
pub fn reevaluate_hsp_with_ambiguities_gapped(
    hit: &mut BlastnHsp,
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
    score_matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> bool {
    let delete = reevaluate_hsp_with_ambiguities_gapped_ex(
        hit,
        q_seq,
        s_seq,
        reward,
        penalty,
        gap_open,
        gap_extend,
        cutoff_score,
        score_matrix,
        None,
    );
    if delete {
        return true;
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:656-663
    // ```c
    // delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(...);
    // if (!delete_hsp)
    //     delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp,
    //                                                 query_nomask, subject,
    //                                                 score_options, hit_options);
    // ```
    blast_hsp_test_identity_and_length(hit, q_seq, s_seq, 0.0, 0)
}

/// Extended version with optional Karlin parameters for E-value recalculation.
/// NCBI reference: blast_traceback.c:234-250 Blast_HSPListGetEvalues, Blast_HSPListGetBitScores
pub fn reevaluate_hsp_with_ambiguities_gapped_ex(
    hit: &mut BlastnHsp,
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
    score_matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
    reeval_params: Option<&ReevalParams>,
) -> bool {
    // NCBI reference: blast_hits.c:539-540
    // esp = hsp->gap_info;
    // if (!esp) return TRUE;
    let gap_ops = match hit.gap_info.as_mut() {
        Some(info) if !info.is_empty() => info,
        _ => return true, // No gap_info = delete HSP
    };

    // NCBI reference: blast_hits.c:528-535
    // query = q + hsp->query.offset;
    // subject = s + hsp->subject.offset;
    let q_offset = hit.internal_q_offset_0;
    let s_offset = hit.internal_s_offset_0;

    // NCBI reference: blast_hits.c:508-526 (factor and gap penalties)
    let mut factor = 1;
    let (gap_open_eval, gap_extend_eval) = if gap_open == 0 && gap_extend == 0 {
        if reward % 2 == 1 {
            factor = 2;
        }
        (0, (reward - 2 * penalty) * factor / 2)
    } else {
        (gap_open, gap_extend)
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:507-514
    // ```c
    // const Uint1* query,* subject;
    // Int4** matrix;
    // ...
    // matrix = sbp->matrix->data;
    // ```
    let matrix = score_matrix;

    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // NCBI uses raw pointers without bounds checks (blast_hits.c:528-559).
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:562-608
    // ```c
    // if (sum < 0) {
    //     if (op_index < esp->num[index]) {
    //         esp->num[index] -= op_index;
    //         current_start_esp_index = index;
    //         op_index = 0;
    //     } else {
    //         current_start_esp_index = index + 1;
    //     }
    //     sum = 0;
    //     current_q_start = query;
    //     current_s_start = subject;
    //     if (score < cutoff_score) {
    //         best_q_start = query;
    //         best_s_start = subject;
    //         score = 0;
    //         best_start_esp_index = current_start_esp_index;
    //         best_end_esp_index = current_start_esp_index;
    //     }
    // } else if (sum > score) {
    //     score = sum;
    //     best_q_start = current_q_start;
    //     best_s_start = current_s_start;
    //     best_q_end = query;
    //     best_s_end = subject;
    //     best_start_esp_index = current_start_esp_index;
    //     best_end_esp_index = index;
    //     best_end_esp_num = op_index;
    // }
    // ```
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

    // NCBI reference: blast_hits.c:541-610
    // ```c
    // if (op_index < esp->num[index]) {
    //     esp->num[index] -= op_index;
    //     current_start_esp_index = index;
    //     op_index = 0;
    // } else {
    //     current_start_esp_index = index + 1;
    // }
    // ```
    for index in 0..gap_ops.len() {
        let mut op_index = 0usize;
        while op_index < gap_ops[index].num() as usize {
            match gap_ops[index] {
                GapEditOp::Sub(_) => {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:547-551
                    // ```c
                    // sum += factor*matrix[*query & kResidueMask][*subject];
                    // query++;
                    // subject++;
                    // op_index++;
                    // ```
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
                                    // NCBI kResidueMask = 0x0f -> indices 0..15
                                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:510-511
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

    // NCBI reference: blast_hits.c:612-616
    score /= factor;

    // NCBI reference: blast_hits.c:617-638
    if best_start_esp_index < gap_ops.len() && best_end_esp_index < gap_ops.len() {
        debug_assert!(matches!(gap_ops[best_start_esp_index], GapEditOp::Sub(_)));
        debug_assert!(matches!(gap_ops[best_end_esp_index], GapEditOp::Sub(_)));

        let mut qp = best_q_start;
        let mut sp = best_s_start;
        let mut ext = 0usize;
        while qp > 0 && sp > 0 {
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
        while qp < q_seq.len() && sp < s_seq.len() && q_seq[qp] < 4 {
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:640-646
    // ```c
    // return s_UpdateReevaluatedHSP(hsp, TRUE, cutoff_score,
    //                               score, q, s, best_q_start,
    //                               best_q_end, best_s_start, best_s_end,
    //                               best_start_esp_index, best_end_esp_index,
    //                               best_end_esp_num);
    // ```
    if score < cutoff_score {
        return true;
    }

    // NCBI reference: blast_hits.c:414-470 s_UpdateReevaluatedHSP
    hit.raw_score = score;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:452-455
    // ```c
    // hsp->query.offset = (Int4)(best_q_start - query_start);
    // hsp->query.end = (Int4)(hsp->query.offset + best_q_end - best_q_start);
    // hsp->subject.offset = (Int4)(best_s_start - subject_start);
    // hsp->subject.end = (Int4)(hsp->subject.offset + best_s_end - best_s_start);
    // ```
    hit.internal_q_offset_0 = best_q_start;
    hit.internal_q_end_0 = best_q_start.saturating_add(best_q_end.saturating_sub(best_q_start));
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:460-470
    // ```c
    // if (best_end_esp_index != last_num || best_start_esp_index > 0) {
    //     GapEditScript* esp_temp = GapEditScriptNew(...);
    //     GapEditScriptPartialCopy(esp_temp, 0, hsp->gap_info, ...);
    //     hsp->gap_info = GapEditScriptDelete(hsp->gap_info);
    //     hsp->gap_info = esp_temp;
    // }
    // last_num = hsp->gap_info->size - 1;
    // hsp->gap_info->num[last_num] = best_end_esp_num;
    // ```
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
    // NCBI reference: blast_traceback.c:234-250 Blast_HSPListGetEvalues, Blast_HSPListGetBitScores
    if let Some(params) = reeval_params {
        const LN2: f64 = 0.69314718055994530941723212145818;

        // Bit score calculation
        // NCBI reference: blast_stat.c BLAST_KarlinStoE_simple
        let log_k = params.k.ln();
        hit.bit_score = ((params.lambda * (score as f64)) - log_k) / LN2;

        // E-value calculation
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1864-1870
        // ```c
        // /* Round score down to even number for E-value calculations only. */
        // score = hsp->score;
        // if (hsp_list && hsp_list->hspcnt != 0
        //         && gapped_calculation && sbp->round_down) {
        //     score &= ~1;
        // }
        // ```
        let evalue_score = if params.round_down_evalue_score {
            score & !1
        } else {
            score
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4157-4171
        // ```c
        // return (double) searchsp * exp((double)(-Lambda * S) + kbp->logK);
        // ```
        // E = K * searchsp * exp(-lambda * score)
        if params.eff_searchsp > 0 {
            hit.e_value = params.k
                * (params.eff_searchsp as f64)
                * (-params.lambda * (evalue_score as f64)).exp();
        }
    }

    false
}
