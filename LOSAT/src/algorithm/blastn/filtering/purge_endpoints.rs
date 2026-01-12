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

use crate::common::{GapEditOp, Hit};
use rustc_hash::FxHashMap;

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
/// * `q_cut` - Query cut position (absolute, 1-based)
/// * `s_cut` - Subject cut position (absolute, 1-based canonical = min of s_start/s_end)
/// * `cut_begin` - If true, keep end portion; if false, keep start portion
///
/// # Returns
/// * `true` if trimming was successful
/// * `false` if no gap_info or cut point not found
pub fn cut_off_gap_edit_script(hit: &mut Hit, q_cut: usize, s_cut: usize, cut_begin: bool) -> bool {
    // NCBI reference: blast_hits.c:2392-2452
    let gap_info = match &hit.gap_info {
        Some(info) if !info.is_empty() => info.clone(),
        _ => return false, // No gap_info, cannot trim
    };

    // Get canonical offsets (NCBI uses subject.offset < subject.end)
    let q_offset = hit.q_start;
    let s_offset = hit.s_start.min(hit.s_end);

    // Convert absolute cut positions to relative (within HSP)
    // NCBI: q_cut -= hsp->query.offset; s_cut -= hsp->subject.offset;
    let q_cut_rel = q_cut.saturating_sub(q_offset);
    let s_cut_rel = s_cut.saturating_sub(s_offset);

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
        hit.q_start = q_offset + qid;
        // For subject, we need to handle strand correctly
        if hit.s_start <= hit.s_end {
            hit.s_start = s_offset + sid;
        } else {
            // Minus strand: s_start > s_end, so we update s_end (the lower coordinate)
            hit.s_end = s_offset + sid;
        }
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
        hit.q_end = q_offset + qid;
        // For subject, we need to handle strand correctly
        if hit.s_start <= hit.s_end {
            hit.s_end = s_offset + sid;
        } else {
            // Minus strand: s_start > s_end, so we update s_start (the higher coordinate)
            hit.s_start = s_offset + sid;
        }
    }

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
/// This function groups hits by subject_id before purging to match NCBI behavior.
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
/// - context = query_id + query strand (we use query_id as first group key)
/// - subject.frame = subject strand (+1 or -1), encoded by s_start > s_end
/// IMPORTANT: Use ACTUAL s_start/s_end, NOT canonical (min/max).
///
/// # Returns
/// Returns the index of the first trimmed HSP (for re-evaluation in purge=false mode).
/// In purge=true mode, this is always equal to the final hit count.
pub fn purge_hsps_with_common_endpoints_ex(hits: Vec<Hit>, purge: bool) -> (Vec<Hit>, usize) {
    let len = hits.len();
    if len <= 1 {
        return (hits, len);
    }

    // NCBI reference: blast_hits.c:2455-2535
    // Blast_HSPListPurgeHSPsWithCommonEndpoints operates on BlastHSPList which is per-subject.
    // We must group hits by subject_id and process each group separately.

    // Group hits by subject_id
    let mut subject_groups: FxHashMap<String, Vec<Hit>> = FxHashMap::default();
    for hit in hits {
        subject_groups.entry(hit.subject_id.clone()).or_default().push(hit);
    }

    // Process each subject group independently and collect results
    let mut result: Vec<Hit> = Vec::new();
    let mut total_extra_start = 0usize;

    for (_subject_id, group_hits) in subject_groups {
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
pub fn purge_hsps_with_common_endpoints(hits: Vec<Hit>) -> Vec<Hit> {
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
fn purge_hsps_for_subject_ex(mut hits: Vec<Hit>, purge: bool) -> (Vec<Hit>, usize) {
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
    let mut trimmed_hits: Vec<Hit> = Vec::new();

    // Helper to determine strand: minus if s_start > s_end
    // NCBI reference: For blastn, subject.frame encodes strand direction
    let is_minus_strand = |h: &Hit| h.s_start > h.s_end;

    // NCBI uses CANONICAL coordinates: subject.offset < subject.end always
    // ASSERT(hsp->subject.offset < hsp->subject.end) at blast_engine.c:1312
    // We need to use min/max for comparison to match NCBI behavior
    let s_offset = |h: &Hit| h.s_start.min(h.s_end);  // NCBI subject.offset
    let s_end_canon = |h: &Hit| h.s_start.max(h.s_end);  // NCBI subject.end

    // Pass 1: Remove HSPs with common START positions
    // NCBI reference: blast_hits.c:2268-2291 s_QueryOffsetCompareHSPs
    // Sort order: context ASC, query.offset ASC, subject.offset ASC,
    //             score DESC, query.end ASC, subject.end ASC
    hits.sort_by(|a, b| {
        a.query_id.cmp(&b.query_id)                       // context ASC
            .then_with(|| a.q_start.cmp(&b.q_start))      // query.offset ASC
            .then_with(|| s_offset(a).cmp(&s_offset(b)))  // subject.offset ASC
            .then_with(|| b.raw_score.cmp(&a.raw_score))  // score DESC
            .then_with(|| a.q_end.cmp(&b.q_end))          // query.end ASC (FIXED: was DESC)
            .then_with(|| s_end_canon(a).cmp(&s_end_canon(b))) // subject.end ASC (FIXED: was DESC)
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

        let same_query = hits[i].query_id == hits[j].query_id;
        let same_strand = is_minus_strand(&hits[i]) == is_minus_strand(&hits[j]);
        let same_q_start = hits[i].q_start == hits[j].q_start;
        let same_s_offset = s_offset(&hits[i]) == s_offset(&hits[j]);

        if same_query && same_strand && same_q_start && same_s_offset {
            // Found duplicate - either trim or delete
            // NCBI: hsp = hsp_array[i+j] is the lower-scoring one to remove/trim
            let mut removed_hit = hits.remove(j);

            if !purge && removed_hit.q_end > hits[i].q_end {
                // NCBI: s_CutOffGapEditScript(hsp, hsp_array[i]->query.end, hsp_array[i]->subject.end, TRUE)
                // Trim the beginning, keep the end portion
                // NCBI reference: blast_hits.c:2480-2498
                // CRITICAL: NCBI's s_CutOffGapEditScript is a VOID function - it always moves
                // the HSP to the end for re-evaluation, even if trimming doesn't change anything.
                // LOSAT must match this behavior: always move to end when condition is met.
                let q_cut = hits[i].q_end;
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
        a.query_id.cmp(&b.query_id)
            .then_with(|| a.q_end.cmp(&b.q_end))
            .then_with(|| s_end_canon(a).cmp(&s_end_canon(b)))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // score DESC
            .then_with(|| b.q_start.cmp(&a.q_start)) // query.offset DESC
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

        let same_query = hits[i].query_id == hits[j].query_id;
        let same_strand = is_minus_strand(&hits[i]) == is_minus_strand(&hits[j]);
        let same_q_end = hits[i].q_end == hits[j].q_end;
        let same_s_end = s_end_canon(&hits[i]) == s_end_canon(&hits[j]);

        if same_query && same_strand && same_q_end && same_s_end {
            // Found duplicate - either trim or delete
            let mut removed_hit = hits.remove(j);

            if !purge && removed_hit.q_start < hits[i].q_start {
                // NCBI: s_CutOffGapEditScript(hsp, hsp_array[i]->query.offset, hsp_array[i]->subject.offset, FALSE)
                // Trim the end, keep the start portion
                // NCBI reference: blast_hits.c:2516-2524
                // CRITICAL: NCBI's s_CutOffGapEditScript is a VOID function - always move to end
                let q_cut = hits[i].q_start;
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
/// sub-alignment. It handles ambiguous bases (N's) which score 0 and can
/// create low-scoring regions within an alignment.
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
/// * `q_seq` - Query sequence (full, 0-indexed)
/// * `s_seq` - Subject sequence (full, 0-indexed)
/// * `reward` - Match reward (positive)
/// * `penalty` - Mismatch penalty (negative)
/// * `gap_open` - Gap open penalty (positive)
/// * `gap_extend` - Gap extend penalty (positive)
/// * `cutoff_score` - Minimum score to keep HSP
/// * `reeval_params` - Optional parameters for E-value/bit score recalculation
///
/// # Returns
/// * `true` if HSP should be DELETED (score < cutoff or no valid region found)
/// * `false` if HSP is OK (score >= cutoff, hit is updated in place)
pub fn reevaluate_hsp_with_ambiguities_gapped(
    hit: &mut Hit,
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
) -> bool {
    reevaluate_hsp_with_ambiguities_gapped_ex(
        hit, q_seq, s_seq, reward, penalty, gap_open, gap_extend, cutoff_score, None
    )
}

/// Extended version with optional Karlin parameters for E-value recalculation.
/// NCBI reference: blast_traceback.c:234-250 Blast_HSPListGetEvalues, Blast_HSPListGetBitScores
pub fn reevaluate_hsp_with_ambiguities_gapped_ex(
    hit: &mut Hit,
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    cutoff_score: i32,
    reeval_params: Option<&ReevalParams>,
) -> bool {
    // NCBI reference: blast_hits.c:539-540
    // esp = hsp->gap_info;
    // if (!esp) return TRUE;
    let gap_info = match &hit.gap_info {
        Some(info) if !info.is_empty() => info,
        _ => return true, // No gap_info = delete HSP
    };

    // Get starting positions (0-indexed internal coordinates)
    // NCBI uses 0-indexed offsets for sequence access
    let q_offset = hit.q_start.saturating_sub(1); // Convert 1-based to 0-based
    let s_offset = hit.s_start.min(hit.s_end).saturating_sub(1); // Canonical, 0-based

    let is_minus = hit.s_start > hit.s_end;

    // NCBI reference: blast_hits.c:528-535
    // query = q + hsp->query.offset;
    // subject = s + hsp->subject.offset;
    // score = 0; sum = 0;
    // best_q_start = best_q_end = current_q_start = query;
    // best_s_start = best_s_end = current_s_start = subject;

    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // Track positions in sequences (0-indexed offsets from start of HSP)
    let mut query_pos = q_offset;
    let mut subject_pos = s_offset;

    // Best scoring region (as sequence positions)
    let mut best_q_start = query_pos;
    let mut best_q_end = query_pos;
    let mut best_s_start = subject_pos;
    let mut best_s_end = subject_pos;

    // Current scoring region start
    let mut current_q_start = query_pos;
    let mut current_s_start = subject_pos;

    // NCBI reference: blast_hits.c:541-610
    // for (index=0; index<esp->size; index++)
    for op in gap_info.iter() {
        let num = op.num() as usize;

        match op {
            GapEditOp::Sub(_) => {
                // NCBI reference: blast_hits.c:547-551
                // if (esp->op_type[index] == eGapAlignSub) {
                //     sum += factor*matrix[*query & kResidueMask][*subject];
                //     query++; subject++; op_index++;
                // }
                for _ in 0..num {
                    if query_pos < q_seq.len() && subject_pos < s_seq.len() {
                        // NCBI uses scoring matrix; we use simple match/mismatch
                        let match_score = if q_seq[query_pos] == s_seq[subject_pos] {
                            reward
                        } else {
                            penalty
                        };
                        sum += match_score;
                    }
                    query_pos += 1;
                    subject_pos += 1;

                    // NCBI reference: blast_hits.c:562-593
                    // Check if sum dropped below 0 (reset region)
                    if sum < 0 {
                        // NCBI reference: blast_hits.c:574-578
                        // current_q_start = query;
                        // current_s_start = subject;
                        current_q_start = query_pos;
                        current_s_start = subject_pos;

                        // NCBI reference: blast_hits.c:584-588
                        // if (score < cutoff_score) {
                        //     best_q_start = query;
                        //     best_s_start = subject;
                        //     score = 0;
                        // }
                        if score < cutoff_score {
                            best_q_start = query_pos;
                            best_s_start = subject_pos;
                            score = 0;
                        }
                        sum = 0;
                    } else if sum > score {
                        // NCBI reference: blast_hits.c:595-604
                        // Remember this as best scoring region
                        score = sum;
                        best_q_start = current_q_start;
                        best_s_start = current_s_start;
                        best_q_end = query_pos;
                        best_s_end = subject_pos;
                    }
                }
            }
            GapEditOp::Del(n) => {
                // NCBI reference: blast_hits.c:552-555
                // } else if (esp->op_type[index] == eGapAlignDel) {
                //     sum -= gap_open + gap_extend * esp->num[index];
                //     subject += esp->num[index];
                // }
                sum -= gap_open + gap_extend * (*n as i32);
                subject_pos += *n as usize;

                // Check for region reset
                if sum < 0 {
                    current_q_start = query_pos;
                    current_s_start = subject_pos;
                    if score < cutoff_score {
                        best_q_start = query_pos;
                        best_s_start = subject_pos;
                        score = 0;
                    }
                    sum = 0;
                } else if sum > score {
                    score = sum;
                    best_q_start = current_q_start;
                    best_s_start = current_s_start;
                    best_q_end = query_pos;
                    best_s_end = subject_pos;
                }
            }
            GapEditOp::Ins(n) => {
                // NCBI reference: blast_hits.c:556-559
                // } else if (esp->op_type[index] == eGapAlignIns) {
                //     sum -= gap_open + gap_extend * esp->num[index];
                //     query += esp->num[index];
                // }
                sum -= gap_open + gap_extend * (*n as i32);
                query_pos += *n as usize;

                // Check for region reset
                if sum < 0 {
                    current_q_start = query_pos;
                    current_s_start = subject_pos;
                    if score < cutoff_score {
                        best_q_start = query_pos;
                        best_s_start = subject_pos;
                        score = 0;
                    }
                    sum = 0;
                } else if sum > score {
                    score = sum;
                    best_q_start = current_q_start;
                    best_s_start = current_s_start;
                    best_q_end = query_pos;
                    best_s_end = subject_pos;
                }
            }
        }
    }

    // NCBI reference: blast_hits.c:614-638
    // Post-processing: try to extend further through exact matches
    // This is optional and primarily helps with ambiguous bases at boundaries
    if best_q_start < best_q_end && best_s_start < best_s_end {
        // Extend left through exact matches
        while best_q_start > 0 && best_s_start > 0 {
            let qi = best_q_start - 1;
            let si = best_s_start - 1;
            if qi < q_seq.len() && si < s_seq.len() && q_seq[qi] == s_seq[si] && q_seq[qi] < 4 {
                best_q_start -= 1;
                best_s_start -= 1;
                score += reward;
            } else {
                break;
            }
        }

        // Extend right through exact matches
        while best_q_end < q_seq.len() && best_s_end < s_seq.len() {
            if q_seq[best_q_end] == s_seq[best_s_end] && q_seq[best_q_end] < 4 {
                score += reward;
                best_q_end += 1;
                best_s_end += 1;
            } else {
                break;
            }
        }
    }

    // NCBI reference: blast_hits.c:641-646
    // return s_UpdateReevaluatedHSP(hsp, TRUE, cutoff_score, score, ...);
    // Returns TRUE if HSP should be deleted

    // Check if score passes cutoff
    if score < cutoff_score {
        return true; // Delete HSP
    }

    // Update HSP coordinates with best region (convert back to 1-based)
    hit.q_start = best_q_start + 1;
    hit.q_end = best_q_end;

    if is_minus {
        // For minus strand: s_start > s_end
        hit.s_start = best_s_end; // larger value
        hit.s_end = best_s_start + 1; // smaller value
    } else {
        hit.s_start = best_s_start + 1;
        hit.s_end = best_s_end;
    }

    // Update score
    hit.raw_score = score;

    // Update alignment length
    let q_len = best_q_end.saturating_sub(best_q_start);
    let s_len = best_s_end.saturating_sub(best_s_start);
    hit.length = q_len.max(s_len);

    // Recalculate identity, mismatch, and gapopen from coordinates
    // NCBI reference: blast_hits.c:745-790 s_Blast_HSPGetNumIdentitiesAndPositives
    if hit.length > 0 {
        // Estimate statistics from score and penalties
        let gap_letters = (q_len as i32 - s_len as i32).unsigned_abs() as usize;
        let non_gap_len = hit.length.saturating_sub(gap_letters);

        // Estimate matches from score: score = matches * reward + mismatches * penalty
        // mismatches = non_gap_len - matches
        // score = matches * reward + (non_gap_len - matches) * penalty
        // score = matches * (reward - penalty) + non_gap_len * penalty
        // matches = (score - non_gap_len * penalty) / (reward - penalty)
        let matches_est = if reward > 0 && non_gap_len > 0 && reward != penalty {
            let est = (score - (non_gap_len as i32) * penalty) / (reward - penalty);
            est.max(0).min(non_gap_len as i32) as usize
        } else {
            non_gap_len
        };

        let mismatches_est = non_gap_len.saturating_sub(matches_est);

        hit.identity = (matches_est as f64 / hit.length as f64) * 100.0;
        hit.mismatch = mismatches_est;

        // Estimate gap openings from gap_letters (rough estimate, at least 1 per direction if gaps exist)
        if gap_letters > 0 {
            // Rough estimate: 1 gap opening per gap region (very approximate)
            hit.gapopen = 1;
        } else {
            hit.gapopen = 0;
        }
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
            hit.e_value = params.k * (params.eff_searchsp as f64) * (-params.lambda * (score as f64)).exp();
        }
    }

    false // Keep HSP
}

