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

use crate::common::GapEditOp;
use super::super::hsp::BlastnHsp;
use rustc_hash::FxHashMap;

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
pub fn cut_off_gap_edit_script(hit: &mut BlastnHsp, q_cut: usize, s_cut: usize, cut_begin: bool) -> bool {
    // NCBI reference: blast_hits.c:2392-2452
    let gap_info = match &hit.gap_info {
        Some(info) if !info.is_empty() => info.clone(),
        _ => return false, // No gap_info, cannot trim
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
    // Convert output coordinates back to internal offsets for trimming.
    let (mut q_offset_0, mut q_end_0) = if hit.query_length > 0 && hit.query_frame < 0 {
        (
            hit.query_length.saturating_sub(hit.q_end),
            hit.query_length
                .saturating_sub(hit.q_start)
                .saturating_add(1),
        )
    } else {
        (hit.q_start.saturating_sub(1), hit.q_end)
    };
    let (mut s_offset_0, mut s_end_0) = (
        hit.s_start.min(hit.s_end).saturating_sub(1),
        hit.s_start.max(hit.s_end),
    );

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
        return false;
    }

    // Now trim the edit script based on cut_begin flag
    // NCBI: if (cut_begin) { ... } else { ... }
    let new_gap_info: Vec<GapEditOp>;

    if cut_begin {
        // Keep END portion, update offsets
        // NCBI reference: blast_hits.c:2426-2441
        let mut new_ops = Vec::new();

        // If cut is in middle of a Sub run, keep remaining part
        let current_op = &gap_info[found_index];
        if let GapEditOp::Sub(n) = current_op {
            if found_opid < *n as usize {
                // NCBI: ASSERT(esp->op_type[index] == eGapAlignSub);
                // esp->num[0] = esp->num[index] - opid;
                let remaining = *n as usize - found_opid;
                if remaining > 0 {
                    new_ops.push(GapEditOp::Sub(remaining as u32));
                }
            }
        }

        // Copy remaining operations
        // NCBI: for (; index < esp->size; index++, new_index++) { ... }
        for op in gap_info.iter().skip(found_index + 1) {
            new_ops.push(*op);
        }

        new_gap_info = new_ops;

        // Update HSP coordinates
        // NCBI: hsp->query.offset += qid; hsp->subject.offset += sid;
        q_offset_0 = q_offset_0.saturating_add(qid);
        s_offset_0 = s_offset_0.saturating_add(sid);
    } else {
        // Keep START portion, update ends
        // NCBI reference: blast_hits.c:2442-2450
        let mut new_ops = Vec::new();

        // Copy operations up to the cut point
        for (i, op) in gap_info.iter().enumerate() {
            if i < found_index {
                new_ops.push(*op);
            } else if i == found_index {
                // If cut is in middle of a Sub run, truncate it
                if let GapEditOp::Sub(n) = op {
                    if found_opid < *n as usize {
                        // NCBI: esp->num[index] = opid;
                        if found_opid > 0 {
                            new_ops.push(GapEditOp::Sub(found_opid as u32));
                        }
                    } else {
                        // Full run consumed
                        new_ops.push(*op);
                    }
                } else {
                    // Del/Ins: include full run (already consumed)
                    new_ops.push(*op);
                }
                break;
            }
        }

        new_gap_info = new_ops;

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
    hit.q_start = q_start;
    hit.q_end = q_end;
    hit.s_start = s_start;
    hit.s_end = s_end;

    // Update gap_info
    hit.gap_info = if new_gap_info.is_empty() {
        None
    } else {
        Some(new_gap_info)
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
pub fn purge_hsps_with_common_endpoints_ex(hits: Vec<BlastnHsp>, purge: bool) -> (Vec<BlastnHsp>, usize) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len);
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
fn purge_hsps_for_subject_ex(mut hits: Vec<BlastnHsp>, purge: bool) -> (Vec<BlastnHsp>, usize) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len);
    }

    let initial_count = hits.len();
    let mut start_purged = 0usize;
    let mut end_purged = 0usize;
    let mut start_trimmed = 0usize;
    let mut end_trimmed = 0usize;

    // Track trimmed HSPs that need re-evaluation (when purge=false)
    let mut trimmed_hits: Vec<BlastnHsp> = Vec::new();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
    // ```c
    // if (hsp->query.frame != hsp->subject.frame) {
    //    *q_end = query_length - hsp->query.offset;
    //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
    // }
    // ```
    let q_offsets = |h: &BlastnHsp| {
        if h.query_length > 0 && h.query_frame < 0 {
            (
                h.query_length.saturating_sub(h.q_end),
                h.query_length
                    .saturating_sub(h.q_start)
                    .saturating_add(1),
            )
        } else {
            (h.q_start.saturating_sub(1), h.q_end)
        }
    };
    let context = |h: &BlastnHsp| -> u32 {
        h.q_idx * 2 + if h.query_frame < 0 { 1 } else { 0 }
    };
    // NCBI uses CANONICAL coordinates: subject.offset < subject.end always
    // ASSERT(hsp->subject.offset < hsp->subject.end) at blast_engine.c:1312
    let s_offset = |h: &BlastnHsp| h.s_start.min(h.s_end).saturating_sub(1);  // NCBI subject.offset
    let s_end_canon = |h: &BlastnHsp| h.s_start.max(h.s_end);  // NCBI subject.end

    // Pass 1: Remove HSPs with common START positions
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2285-2318
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
    hits.sort_by(|a, b| {
        let (a_q_offset, a_q_end) = q_offsets(a);
        let (b_q_offset, b_q_end) = q_offsets(b);
        context(a)
            .cmp(&context(b)) // context ASC
            .then_with(|| a_q_offset.cmp(&b_q_offset)) // query.offset ASC
            .then_with(|| s_offset(a).cmp(&s_offset(b))) // subject.offset ASC
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // score DESC
            .then_with(|| b_q_end.cmp(&a_q_end)) // query.end DESC
            .then_with(|| s_end_canon(b).cmp(&s_end_canon(a))) // subject.end DESC
    });

    // NCBI reference: blast_hits.c:2480-2500
    // if (!purge && (hsp->query.end > hsp_array[i]->query.end)) {
    //     s_CutOffGapEditScript(hsp, hsp_array[i]->query.end, hsp_array[i]->subject.end, TRUE);
    // } else {
    //     hsp = Blast_HSPFree(hsp);
    // }
    let mut i = 0;
    while i < hits.len() {
        let j = i + 1;
        if j >= hits.len() {
            i += 1;
            continue;
        }

        let (i_q_offset, i_q_end) = q_offsets(&hits[i]);
        let (j_q_offset, j_q_end) = q_offsets(&hits[j]);
        let same_context = context(&hits[i]) == context(&hits[j]);
        let same_q_start = i_q_offset == j_q_offset;
        let same_s_offset = s_offset(&hits[i]) == s_offset(&hits[j]);

        if same_context && same_q_start && same_s_offset {
            // Found duplicate - either trim or delete
            // NCBI: hsp = hsp_array[i+j] is the lower-scoring one to remove/trim
            let mut removed_hit = hits.remove(j);

            if !purge && j_q_end > i_q_end {
                // NCBI: s_CutOffGapEditScript(hsp, hsp_array[i]->query.end, hsp_array[i]->subject.end, TRUE)
                // Trim the beginning, keep the end portion
                // NCBI reference: blast_hits.c:2480-2498
                // CRITICAL: NCBI's s_CutOffGapEditScript is a VOID function - it always moves
                // the HSP to the end for re-evaluation, even if trimming doesn't change anything.
                // LOSAT must match this behavior: always move to end when condition is met.
                let q_cut = i_q_end;
                let s_cut = s_end_canon(&hits[i]);
                let _ = cut_off_gap_edit_script(&mut removed_hit, q_cut, s_cut, true);
                // Always move to end for re-evaluation (matches NCBI void function behavior)
                trimmed_hits.push(removed_hit);
                start_trimmed += 1;
            }
            // else: purge=true or hsp doesn't extend beyond, just delete (already removed)

            start_purged += 1;
            // Don't increment i - check next element at same position
        } else {
            i += 1;
        }
    }

    // Pass 2: Remove HSPs with common END positions
    // NCBI reference: blast_hits.c:2504-2526
    hits.sort_by(|a, b| {
        let (a_q_offset, a_q_end) = q_offsets(a);
        let (b_q_offset, b_q_end) = q_offsets(b);
        context(a)
            .cmp(&context(b))
            .then_with(|| a_q_end.cmp(&b_q_end))
            .then_with(|| s_end_canon(a).cmp(&s_end_canon(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // score DESC
            .then_with(|| b_q_offset.cmp(&a_q_offset)) // query.offset DESC
            .then_with(|| s_offset(b).cmp(&s_offset(a))) // subject.offset DESC
    });

    // NCBI reference: blast_hits.c:2516-2520
    // if (!purge && (hsp->query.offset < hsp_array[i]->query.offset)) {
    //     s_CutOffGapEditScript(hsp, hsp_array[i]->query.offset, hsp_array[i]->subject.offset, FALSE);
    // } else {
    //     hsp = Blast_HSPFree(hsp);
    // }
    let mut i = 0;
    while i < hits.len() {
        let j = i + 1;
        if j >= hits.len() {
            i += 1;
            continue;
        }

        let (i_q_offset, i_q_end) = q_offsets(&hits[i]);
        let (j_q_offset, j_q_end) = q_offsets(&hits[j]);
        let same_context = context(&hits[i]) == context(&hits[j]);
        let same_q_end = i_q_end == j_q_end;
        let same_s_end = s_end_canon(&hits[i]) == s_end_canon(&hits[j]);

        if same_context && same_q_end && same_s_end {
            // Found duplicate - either trim or delete
            let mut removed_hit = hits.remove(j);

            if !purge && j_q_offset < i_q_offset {
                // NCBI: s_CutOffGapEditScript(hsp, hsp_array[i]->query.offset, hsp_array[i]->subject.offset, FALSE)
                // Trim the end, keep the start portion
                // NCBI reference: blast_hits.c:2516-2524
                // CRITICAL: NCBI's s_CutOffGapEditScript is a VOID function - always move to end
                let q_cut = i_q_offset;
                let s_cut = s_offset(&hits[i]);
                let _ = cut_off_gap_edit_script(&mut removed_hit, q_cut, s_cut, false);
                // Always move to end for re-evaluation (matches NCBI void function behavior)
                trimmed_hits.push(removed_hit);
                end_trimmed += 1;
            }
            // else: purge=true or hsp doesn't start before, just delete

            end_purged += 1;
            // Don't increment i - check next element at same position
        } else {
            i += 1;
        }
    }

    // Calculate the index where trimmed HSPs start
    let extra_start = hits.len();

    // Append trimmed HSPs to the end (for re-evaluation in purge=false mode)
    // NCBI reference: blast_hits.c:2496-2499, 2522-2525
    // for (k=i+j; k<hsp_count; k++) { hsp_array[k] = hsp_array[k+1]; }
    // hsp_array[hsp_count] = hsp;  // Move to end
    hits.extend(trimmed_hits);

    // Debug output for purge statistics
    if std::env::var("LOSAT_DEBUG_BLASTN").is_ok() {
        eprintln!(
            "[DEBUG_PURGE] initial={}, start_purged={} (trimmed={}), end_purged={} (trimmed={}), final={}",
            initial_count, start_purged, start_trimmed, end_purged, end_trimmed, hits.len()
        );
    }

    let _ = (initial_count, start_purged, end_purged, start_trimmed, end_trimmed); // Silence unused warnings

    (hits, extra_start)
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
                query_frame: 1,
                query_length: 100,
                q_idx: 0,
                s_idx: 0,
                raw_score: 100,
                gap_info: None,
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
    let query_is_minus = hit.query_frame < 0;
    let q_offset = if hit.query_length > 0 && query_is_minus {
        hit.query_length.saturating_sub(hit.q_end)
    } else {
        hit.q_start.saturating_sub(1)
    };
    let s_offset = hit.s_start.min(hit.s_end).saturating_sub(1);

    let (matches, mismatches, gap_opens, _gap_letters, align_length) =
        if let Some(ref ops) = hit.gap_info {
            stats_from_edit_ops(q_seq, s_seq, q_offset, s_offset, ops)
        } else {
            // Ungapped fallback (N.B. not expected for BLASTN reevaluation path).
            let align_length = if query_is_minus {
                hit.q_start
                    .saturating_sub(hit.q_end)
                    .saturating_add(1)
            } else {
                hit.q_end
                    .saturating_sub(hit.q_start)
                    .saturating_add(1)
            };
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

    let query_is_minus = hit.query_frame < 0;

    // NCBI reference: blast_hits.c:528-535
    // query = q + hsp->query.offset;
    // subject = s + hsp->subject.offset;
    let q_offset = if hit.query_length > 0 && query_is_minus {
        hit.query_length.saturating_sub(hit.q_end)
    } else {
        hit.q_start.saturating_sub(1)
    };
    let s_offset = hit.s_start.min(hit.s_end).saturating_sub(1);

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
                                    sum += factor * matrix[q_code * BLASTNA_SIZE + s_code];
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
                    if query_pos < q_seq.len() && subject_pos < s_seq.len() {
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
        // NCBI reference: blast_stat.c BLAST_KarlinStoE_simple
        // E = K * searchsp * exp(-lambda * score)
        if params.eff_searchsp > 0 {
            hit.e_value =
                params.k * (params.eff_searchsp as f64) * (-params.lambda * (score as f64)).exp();
        }
    }

    false
}

