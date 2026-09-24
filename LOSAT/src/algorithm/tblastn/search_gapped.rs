//! Internal TBLASTN protein gapped-score boundary.
//! The public TBLASTN CLI remains unsupported until downstream stages agree.

use super::search_init::InitHsp;
use super::search_seed::resolve_local_subject_ncbi2na;
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::blastp::gapalign::{
    blast_gapped_alignment_with_traceback_with_scratch, blastp_get_start_for_gapped_alignment,
    blastp_score_only_gapped_alignment_with_scratch, BlastpGappedAlignmentMode, GapAlignScratch,
};
use crate::algorithm::tblastx::translation::generate_frames;
use crate::config::ScoringMatrix;
use crate::utils::genetic_code::GeneticCode;
use anyhow::{bail, Context, Result};

// NCBI c++/include/algo/blast/core/blast_hits.h:96-103,126-143:
// typedef struct BlastSeg { Int2 frame; Int4 offset, end, gapped_start; } BlastSeg;
// typedef struct BlastHSP { Int4 score; ... BlastSeg query, subject; } BlastHSP;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct GappedHsp {
    pub frame: i8,
    pub score: i32,
    pub q_start: i32,
    pub q_end: i32,
    pub q_gapped_start: i32,
    pub s_start: i32,
    pub s_end: i32,
    pub s_gapped_start: i32,
}

// NCBI c++/src/algo/blast/core/blast_engine.c:522-552:
// status = aux_struct->GetGappedScore(..., init_hitlist, &hsp_list, ...);
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// Blast_HSPListSortByScore(hsp_list);
// NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
// max_offset = BlastGetStartForGappedAlignment(query_tmp.sequence,
//     subject->sequence, gap_align->sbp, q_start, length, s_start, length);
// init_hsp->offsets.qs_offsets.s_off +=
//     max_offset - init_hsp->offsets.qs_offsets.q_off;
// init_hsp->offsets.qs_offsets.q_off = max_offset;
// status = s_BlastProtGappedAlignment(...);
// This internal entry currently covers the pinned BLOSUM62/word-3 score-only
// profile and preserves the frame-local coordinates returned by GetGappedScore.
#[allow(dead_code)]
pub(super) fn gapped_blosum62_word3_hsps(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    initial: &[InitHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    cutoff_score: i32,
) -> Result<Vec<GappedHsp>> {
    Ok(gapped_blosum62_word3_hsps_multi(
        &[query],
        subject,
        db_gencode,
        initial,
        gap_open,
        gap_extend,
        x_drop,
        &[cutoff_score],
    )?
    .into_iter()
    .map(|(_, hsp)| hsp)
    .collect())
}

// NCBI c++/src/algo/blast/core/blast_gapalign.c:2371-2468:
// context = BSearchContextInfo(init_hsp->offsets.qs_offsets.q_off, query_info);
// query_start = query_info->contexts[context].query_offset;
// init_hsp->offsets.qs_offsets.q_off -= query_start;
// init_hsp->ungapped_data->q_start -= query_start;
// NCBI c++/src/algo/blast/core/blast_itree.c:212-229:
// s_GetQueryStrandOffset(query_info, context) returns the protein context's
// query_offset, which is also used by containment and insertion.
#[allow(dead_code)]
pub(super) fn gapped_blosum62_word3_hsps_multi(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    initial: &[InitHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    cutoff_scores: &[i32],
) -> Result<Vec<(usize, GappedHsp)>> {
    Ok(gapped_blosum62_word3_hsps_multi_chunks(
        queries,
        subject,
        db_gencode,
        initial,
        gap_open,
        gap_extend,
        x_drop,
        cutoff_scores,
    )?
    .into_iter()
    .map(|(context, hsp, _)| (context, hsp))
    .collect())
}

// NCBI c++/src/algo/blast/core/blast_engine.c:478-552,572-586:
// each GetGappedScore output belongs to the current chunk; backup.offset
// is applied only after endpoint purge, before Blast_HSPListsMerge.
#[allow(dead_code)]
fn gapped_blosum62_word3_hsps_multi_chunks(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    initial: &[InitHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    cutoff_scores: &[i32],
) -> Result<Vec<(usize, GappedHsp, usize)>> {
    gapped_protein_hsps_multi_chunks(
        queries,
        subject,
        db_gencode,
        initial,
        gap_open,
        gap_extend,
        x_drop,
        cutoff_scores,
        ScoringMatrix::Blosum62,
    )
}

// NCBI c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4058:
// cutoff = hit_params->cutoffs[context].cutoff_score;
// BLAST_GetGappedScore uses gap_align->sbp->matrix for the chosen profile.
#[allow(dead_code)]
fn gapped_protein_hsps_multi_chunks(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    initial: &[InitHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    cutoff_scores: &[i32],
    matrix: ScoringMatrix,
) -> Result<Vec<(usize, GappedHsp, usize)>> {
    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3924-3927:
    // cutoff = hit_params->cutoffs[context].cutoff_score;
    // Every query context must supply its own measured cutoff.
    if cutoff_scores.len() != queries.len() {
        bail!("one preliminary gapped cutoff is required per query context");
    }
    let resolved = resolve_local_subject_ncbi2na(subject)?;
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    let frames = generate_frames(&resolved, &code);
    let query_frames: Vec<_> = queries
        .iter()
        .map(|query| encode_protein_query_frame_with_seg(query, None))
        .collect();
    let mut context_offsets = Vec::with_capacity(queries.len());
    let mut total_query_length = 0i32;
    for query in queries {
        context_offsets.push(total_query_length);
        total_query_length += i32::try_from(query.len())? + 1;
    }
    total_query_length -= 1;
    let mut scratch = GapAlignScratch::new();
    let mut results = Vec::new();

    for frame in frames {
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-552:
        // each chunk resets the initial-hit list, calls WordFinder, then
        // GetGappedScore before endpoint purge and offset adjustment.
        for chunk in super::search_seed::unmasked_translated_chunks(frame.aa_len) {
            let subject_sequence = &frame.aa_seq[1 + chunk.offset..1 + chunk.offset + chunk.length];
            // NCBI c++/src/algo/blast/core/aa_ungapped.c:223-235:
            // Blast_InitHitListSortByScore(init_hitlist);
            // NCBI c++/src/algo/blast/core/blast_extend.c:273-313:
            // compare score DESC, subject start ASC, length DESC, query start ASC.
            let mut frame_initial: Vec<_> = initial
                .iter()
                .filter(|hit| hit.frame == frame.frame && hit.chunk_offset == chunk.offset as u32)
                .collect();
            frame_initial.sort_unstable_by(|a, b| {
                b.score
                    .cmp(&a.score)
                    .then(a.s_start.cmp(&b.s_start))
                    .then(b.length.cmp(&a.length))
                    .then(a.q_start.cmp(&b.q_start))
            });
            // NCBI c++/src/algo/blast/core/blast_engine.c:489-493:
            // if (init_hitlist->total == 0) continue; GetGappedScore is skipped.
            if frame_initial.is_empty() {
                continue;
            }
            // NCBI c++/src/algo/blast/core/blast_gapalign.c:3827-3834,3908-3919,
            // 4065-4091:
            // tree = Blast_IntervalTreeInit(0, query->length+1, 0, subject->length+1);
            // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
            //                                   hit_options->min_diag_separation)) { ... }
            // BlastIntervalTreeAddHSP(new_hsp, tree, query_info, eQueryAndSubject);
            let mut tree = BlastIntervalTree::new(
                0,
                total_query_length + 1,
                0,
                i32::try_from(chunk.length)? + 1,
            );
            for hit in frame_initial {
                let context =
                    context_offsets.partition_point(|&offset| offset <= hit.q_seed as i32) - 1;
                let context_offset = context_offsets[context];
                let query_sequence = &query_frames[context].aa_seq[1..1 + queries[context].len()];
                let initial_tree_hsp = TreeHsp {
                    query_offset: hit.q_start - context_offset,
                    query_end: hit.q_start - context_offset + hit.length,
                    subject_offset: hit.s_start,
                    subject_end: hit.s_start + hit.length,
                    score: hit.score,
                    query_frame: 0,
                    query_length: i32::try_from(queries[context].len())?,
                    query_context_offset: context_offset,
                    subject_frame_sign: i32::from(frame.frame.signum()),
                };
                if tree.contains_hsp(&initial_tree_hsp, context_offset, 0) {
                    continue;
                }
                let q_start = usize::try_from(hit.q_start - context_offset)
                    .context("negative context-local initial query start")?;
                let s_start =
                    usize::try_from(hit.s_start).context("negative initial subject start")?;
                let length = usize::try_from(hit.length).context("negative initial HSP length")?;
                let q_gapped_start = blastp_get_start_for_gapped_alignment(
                    query_sequence,
                    subject_sequence,
                    q_start,
                    length,
                    s_start,
                    length,
                    matrix,
                );
                let s_gapped_start = i64::from(hit.s_seed) + i64::try_from(q_gapped_start)?
                    - i64::from(hit.q_seed - u32::try_from(context_offset)?);
                let s_gapped_start =
                    usize::try_from(s_gapped_start).context("negative gapped subject start")?;
                let gapped = blastp_score_only_gapped_alignment_with_scratch(
                    query_sequence,
                    subject_sequence,
                    q_gapped_start,
                    s_gapped_start,
                    matrix,
                    gap_open,
                    gap_extend,
                    x_drop,
                    BlastpGappedAlignmentMode::Exact,
                    &mut scratch,
                );
                // NCBI c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4058:
                // cutoff = hit_params->cutoffs[context].cutoff_score;
                // if (gap_align->score >= cutoff) {
                //     Blast_HSPInit(..., gap_align->score, ..., &new_hsp);
                //     Blast_HSPListSaveHSP(hsp_list, new_hsp);
                //     BlastIntervalTreeAddHSP(new_hsp, tree, query_info, eQueryAndSubject);
                // }
                if gapped.score >= cutoff_scores[context] {
                    let saved = GappedHsp {
                        frame: frame.frame,
                        score: gapped.score,
                        q_start: gapped.query_start,
                        q_end: gapped.query_stop,
                        q_gapped_start: i32::try_from(q_gapped_start)?,
                        s_start: gapped.subject_start,
                        s_end: gapped.subject_stop,
                        s_gapped_start: i32::try_from(s_gapped_start)?,
                    };
                    tree.add_hsp(
                        TreeHsp {
                            query_offset: saved.q_start,
                            query_end: saved.q_end,
                            subject_offset: saved.s_start,
                            subject_end: saved.s_end,
                            score: saved.score,
                            query_frame: 0,
                            query_length: i32::try_from(queries[context].len())?,
                            query_context_offset: context_offset,
                            subject_frame_sign: i32::from(frame.frame.signum()),
                        },
                        context_offset,
                        IndexMethod::QueryAndSubject,
                    );
                    results.push((context, saved, chunk.offset));
                }
            }
        }
    }
    Ok(results)
}

// NCBI c++/src/algo/blast/core/blast_hits.c:1330-1382:
// ScoreCompareHSPs sorts score DESC, then subject offset ASC,
// subject end DESC, query offset ASC, query end DESC.
fn sort_preliminary_by_score(hsps: &mut [(usize, GappedHsp)]) {
    if hsps
        .windows(2)
        .any(|pair| compare_gapped_score(&pair[0].1, &pair[1].1).is_gt())
    {
        hsps.sort_unstable_by(|a, b| compare_gapped_score(&a.1, &b.1));
    }
}

// NCBI c++/src/algo/blast/core/blast_hits.c:2268-2379,2455-2537:
// qsort(hsp_array, ..., s_QueryOffsetCompareHSPs); delete equal starts;
// qsort(hsp_array, ..., s_QueryEndCompareHSPs); delete equal ends.
// Context and subject frame must also agree before deletion.
fn purge_preliminary_common_endpoints(
    mut incoming: Vec<(usize, GappedHsp)>,
) -> Vec<(usize, GappedHsp)> {
    incoming.sort_unstable_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.q_start.cmp(&b.1.q_start))
            .then(a.1.s_start.cmp(&b.1.s_start))
            .then(b.1.score.cmp(&a.1.score))
            .then(b.1.q_end.cmp(&a.1.q_end))
            .then(b.1.s_end.cmp(&a.1.s_end))
    });
    let mut i = 0;
    while i + 1 < incoming.len() {
        let a = incoming[i];
        let b = incoming[i + 1];
        if a.0 == b.0
            && a.1.frame == b.1.frame
            && a.1.q_start == b.1.q_start
            && a.1.s_start == b.1.s_start
        {
            incoming.remove(i + 1);
        } else {
            i += 1;
        }
    }
    incoming.sort_unstable_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.q_end.cmp(&b.1.q_end))
            .then(a.1.s_end.cmp(&b.1.s_end))
            .then(b.1.score.cmp(&a.1.score))
            .then(b.1.q_start.cmp(&a.1.q_start))
            .then(b.1.s_start.cmp(&a.1.s_start))
    });
    i = 0;
    while i + 1 < incoming.len() {
        let a = incoming[i];
        let b = incoming[i + 1];
        if a.0 == b.0 && a.1.frame == b.1.frame && a.1.q_end == b.1.q_end && a.1.s_end == b.1.s_end
        {
            incoming.remove(i + 1);
        } else {
            i += 1;
        }
    }
    sort_preliminary_by_score(&mut incoming);
    incoming
}

// NCBI c++/src/algo/blast/core/blast_hits.c:2762-2864:
// if (!combined) adopt incoming; otherwise combine by ScoreCompareHSPs
// up to hsp_num_max, preferring the existing list on a comparator tie.
fn append_preliminary_hsps(
    mut combined: Vec<(usize, GappedHsp)>,
    mut incoming: Vec<(usize, GappedHsp)>,
    hsp_num_max: usize,
) -> Vec<(usize, GappedHsp)> {
    if incoming.is_empty() {
        return combined;
    }
    if combined.is_empty() {
        return incoming;
    }
    let retained = combined
        .len()
        .saturating_add(incoming.len())
        .min(hsp_num_max);
    if retained == combined.len() + incoming.len() {
        combined.extend(incoming);
        sort_preliminary_by_score(&mut combined);
        return combined;
    }
    sort_preliminary_by_score(&mut combined);
    sort_preliminary_by_score(&mut incoming);
    let mut selected = Vec::with_capacity(retained);
    let mut old = 0;
    let mut new = 0;
    while selected.len() < retained {
        if old < combined.len()
            && (new >= incoming.len()
                || compare_gapped_score(&combined[old].1, &incoming[new].1).is_le())
        {
            selected.push(combined[old]);
            old += 1;
        } else {
            selected.push(incoming[new]);
            new += 1;
        }
    }
    selected
}

// NCBI c++/src/algo/blast/core/blast_engine.c:522-552,804-850:
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// Blast_HSPListSortByScore(hsp_list);
// Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
#[allow(dead_code)]
fn preliminary_merged_hsps(
    gapped: &[(usize, GappedHsp)],
    hsp_num_max: usize,
) -> Vec<Vec<(usize, GappedHsp)>> {
    let mut combined = Vec::new();
    let mut snapshots = Vec::new();
    for frame in [1, 2, 3, -1, -2, -3] {
        let incoming: Vec<_> = gapped
            .iter()
            .copied()
            .filter(|(_, hsp)| hsp.frame == frame)
            .collect();
        let incoming = purge_preliminary_common_endpoints(incoming);
        combined = append_preliminary_hsps(combined, incoming, hsp_num_max);
        snapshots.push(combined.clone());
    }
    snapshots
}

// NCBI c++/src/algo/blast/core/blast_engine.c:478-586:
// for each chunk: GetGappedScore, endpoint purge, score sort,
// Blast_HSPListAdjustOffsets, Blast_HSPListsMerge; after all chunks,
// Blast_HSPListAppend retains the frame list under kHspNumMax.
// This covers one unmasked hard range per translated subject frame.
#[allow(dead_code)]
fn preliminary_chunked_hsps(
    gapped: &[(usize, GappedHsp, usize)],
    frame_lengths: &[(i8, usize)],
    hsp_num_max: usize,
) -> Result<Vec<Vec<(usize, GappedHsp)>>> {
    let mut combined = Vec::new();
    let mut snapshots = Vec::new();
    for &(frame, length) in frame_lengths {
        let mut per_frame = Vec::new();
        for chunk in super::search_seed::unmasked_translated_chunks(length) {
            let incoming: Vec<_> = gapped
                .iter()
                .copied()
                .filter(|(_, hsp, offset)| hsp.frame == frame && *offset == chunk.offset)
                .map(|(context, hsp, _)| (context, hsp))
                .collect();
            if incoming.is_empty() {
                continue;
            }
            let mut incoming = purge_preliminary_common_endpoints(incoming);
            // NCBI c++/src/algo/blast/core/blast_hits.c:3037-3050:
            // Blast_HSPListAdjustOffsets adds backup.offset to subject
            // offset, end, and gapped_start before overlap merge.
            let offset = i32::try_from(chunk.offset)?;
            for (_, hsp) in &mut incoming {
                hsp.s_start = hsp
                    .s_start
                    .checked_add(offset)
                    .context("subject offset overflow")?;
                hsp.s_end = hsp
                    .s_end
                    .checked_add(offset)
                    .context("subject end overflow")?;
                hsp.s_gapped_start = hsp
                    .s_gapped_start
                    .checked_add(offset)
                    .context("subject gapped start overflow")?;
            }
            let overlap = if chunk.offset == 0 { 0 } else { 100 };
            per_frame =
                merge_subject_chunk_hsps(per_frame, incoming, offset, overlap, hsp_num_max, true);
        }
        combined = append_preliminary_hsps(combined, per_frame, hsp_num_max);
        snapshots.push(combined.clone());
    }
    Ok(snapshots)
}

// NCBI c++/src/algo/blast/core/blast_hits.c:1465-1537:
// start_diag = query.offset - subject.offset;
// end_diag = query.end - subject.end;
// if (subject.frame differs) return FALSE;
// if either new endpoint is inside the old HSP, score_density is
// (old.score + new.score) / (old.query length + new.query length);
// merged score = MAX((int)(score_density * merged query length), old.score).
fn merge_two_chunk_hsps(old: &mut GappedHsp, new: &GappedHsp, allow_gap: bool) -> bool {
    if !allow_gap && old.s_start - new.s_start - old.q_start + new.q_start != 0 {
        return false;
    }
    if old.frame != new.frame {
        return false;
    }
    let contained =
        |q: i32, s: i32| old.q_start <= q && q <= old.q_end && old.s_start <= s && s <= old.s_end;
    if !contained(new.q_start, new.s_start) && !contained(new.q_end, new.s_end) {
        return false;
    }
    let density = f64::from(old.score + new.score)
        / f64::from((old.q_end - old.q_start) + (new.q_end - new.q_start));
    old.q_start = old.q_start.min(new.q_start);
    old.s_start = old.s_start.min(new.s_start);
    old.q_end = old.q_end.max(new.q_end);
    old.s_end = old.s_end.max(new.s_end);
    if new.score > old.score {
        old.q_gapped_start = new.q_gapped_start;
        old.s_gapped_start = new.s_gapped_start;
        old.score = new.score;
    }
    old.score = ((density * f64::from(old.q_end - old.q_start)) as i32).max(old.score);
    true
}

// NCBI c++/src/algo/blast/core/blast_engine.c:572-586:
// Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
// overlap = backup.offset == hard_range.left ? 0 : DBSEQ_CHUNK_OVERLAP;
// Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, kHspNumMax,
//                     &backup.offset, INT4_MIN, overlap, TRUE, FALSE);
// NCBI c++/src/algo/blast/core/blast_hits.c:2857-3035:
// move overlap HSPs to the front, match equal contexts whose diagonals
// differ by less than 10, merge, purge the consumed new HSPs, combine by score.
// Input HSP coordinates must already include the chunk offset.
#[allow(dead_code)]
fn merge_subject_chunk_hsps(
    mut old: Vec<(usize, GappedHsp)>,
    mut incoming: Vec<(usize, GappedHsp)>,
    split_offset: i32,
    overlap: i32,
    hsp_num_max: usize,
    allow_gap: bool,
) -> Vec<(usize, GappedHsp)> {
    if incoming.is_empty() {
        return old;
    }
    if old.is_empty() {
        return incoming;
    }
    let mut old_overlap = 0;
    for index in 0..old.len() {
        if old[index].1.s_end > split_offset {
            old.swap(old_overlap, index);
            old_overlap += 1;
        }
    }
    let mut new_overlap = 0;
    for index in 0..incoming.len() {
        if incoming[index].1.s_start < split_offset + overlap {
            incoming.swap(new_overlap, index);
            new_overlap += 1;
        }
    }
    let mut incoming: Vec<Option<(usize, GappedHsp)>> = incoming.into_iter().map(Some).collect();
    for old_index in 0..old_overlap {
        for new_index in 0..new_overlap {
            let Some((context, candidate)) = incoming[new_index] else {
                continue;
            };
            if old[old_index].0 != context {
                continue;
            }
            let old_end_diag = old[old_index].1.q_end - old[old_index].1.s_end;
            let new_start_diag = candidate.q_start - candidate.s_start;
            if (old_end_diag - new_start_diag).abs() < 10
                && merge_two_chunk_hsps(&mut old[old_index].1, &candidate, allow_gap)
            {
                incoming[new_index] = None;
            }
        }
    }
    let mut incoming: Vec<_> = incoming.into_iter().flatten().collect();
    let retained = old.len().saturating_add(incoming.len()).min(hsp_num_max);
    // NCBI blast_hits.c:2762-2807, 1330-1382:
    // on no cap, append then SortByScore if out of order;
    // on cap, sort both lists and merge score order, preferring old on ties.
    let sorted = |hsps: &mut Vec<(usize, GappedHsp)>| {
        if hsps
            .windows(2)
            .any(|pair| compare_gapped_score(&pair[0].1, &pair[1].1).is_gt())
        {
            hsps.sort_unstable_by(|a, b| compare_gapped_score(&a.1, &b.1));
        }
    };
    if retained == old.len() + incoming.len() {
        old.extend(incoming);
        sorted(&mut old);
        return old;
    }
    sorted(&mut old);
    sorted(&mut incoming);
    let mut selected = Vec::with_capacity(retained);
    let mut old_index = 0;
    let mut new_index = 0;
    while selected.len() < retained {
        if old_index < old.len()
            && (new_index >= incoming.len()
                || compare_gapped_score(&old[old_index].1, &incoming[new_index].1).is_le())
        {
            selected.push(old[old_index]);
            old_index += 1;
        } else {
            selected.push(incoming[new_index]);
            new_index += 1;
        }
    }
    selected
}

// NCBI c++/src/algo/blast/core/blast_util.c:1268-1317:
// retval->partial = !is_ooframe;
// retval->range = (Int4*) calloc(2*num_frames, sizeof(Int4));
// NCBI c++/src/algo/blast/core/blast_hits.c:1156-1228:
// start = target_t->range[2*context];
// stop = target_t->range[2*context+1];
// return target_t->translations[context] - target_t->range[2*context] + 1;
#[derive(Default)]
struct TargetFrameTranslation {
    start: usize,
    stop: usize,
    base: usize,
    sequence: Vec<u8>,
}

struct TargetTranslation<'a> {
    subject: &'a [u8],
    code: &'a GeneticCode,
    frames: [TargetFrameTranslation; 6],
}

impl<'a> TargetTranslation<'a> {
    fn new(subject: &'a [u8], code: &'a GeneticCode) -> Self {
        Self {
            subject,
            code,
            frames: std::array::from_fn(|_| TargetFrameTranslation::default()),
        }
    }

    // NCBI c++/src/algo/blast/core/blast_hits.c:1154-1228:
    // nucl_start = MAX(0, 3*hsp->subject.offset - kMaxTranslation);
    // nucl_end = MIN(subject->length, 3*hsp->subject.end + kMaxTranslation);
    // if (subject->length - nucl_end <= 21) nucl_end = subject->length;
    // translation_length = 1+nucl_length/CODON_LENGTH;
    // start_shift = nucl_start/CODON_LENGTH;
    // if (start_shift < start || start_shift+translation_length > stop) {
    //     length = BLAST_GetTranslation(...);
    //     target_t->range[2*context+1] = start_shift + length;
    //     translations[context][0] = FENCE_SENTRY;
    //     translations[context][length+1] = FENCE_SENTRY;
    // }
    // NCBI c++/src/algo/blast/core/blast_util.c:428-455:
    // prot_seq[0] = NULLB; ... prot_seq[index_prot] = NULLB;
    fn get(&mut self, frame: i8, offset: i32, end: i32) -> Result<(&[u8], usize, usize)> {
        let context = if frame > 0 { frame - 1 } else { 2 - frame };
        let context = usize::try_from(context).context("invalid subject frame")?;
        if context >= self.frames.len() || frame == 0 {
            bail!("invalid subject frame {frame}");
        }
        let target = &mut self.frames[context];
        let n = self.subject.len();
        if target.start != 0 || target.stop < (n / 3).saturating_sub(3) {
            let (nucl_start, mut nucl_end) = if offset < 0 {
                (0, n)
            } else {
                (
                    usize::try_from(offset)?
                        .saturating_mul(3)
                        .saturating_sub(99),
                    usize::try_from(end)?
                        .saturating_mul(3)
                        .saturating_add(99)
                        .min(n),
                )
            };
            if n - nucl_end <= 21 {
                nucl_end = n;
            }
            let nucl_length = nucl_end - nucl_start;
            let translation_length = 1 + nucl_length / 3;
            let start_shift = nucl_start / 3;
            if start_shift < target.start || start_shift + translation_length > target.stop {
                let nucl_shift = if frame < 0 { n - nucl_end } else { nucl_start };
                let translated = generate_frames(
                    &self.subject[nucl_shift..nucl_shift + nucl_length],
                    self.code,
                );
                let translated = translated
                    .iter()
                    .find(|candidate| candidate.frame == frame)
                    .context("missing translated subject frame")?;
                target.start = start_shift;
                target.stop = start_shift + translated.aa_len;
                // NCBI c++/src/algo/blast/core/blast_hits.c:1213-1228:
                // translations[context][0] = FENCE_SENTRY;
                // return translations[context] - range[2*context] + 1;
                // The Rust slice begins at the left fence and `base` restores
                // NCBI's absolute subject coordinates without allocating the
                // unused prefix preceding a late partial window.
                target.base = if offset >= 0 && start_shift > 0 {
                    start_shift - 1
                } else {
                    0
                };
                target.sequence = Vec::with_capacity(translated.aa_len + 2);
                if offset >= 0 && start_shift > 0 {
                    target.sequence.push(201);
                }
                target
                    .sequence
                    .extend_from_slice(&translated.aa_seq[1..1 + translated.aa_len]);
                target.sequence.push(if offset >= 0 { 201 } else { 0 });
            }
        }
        Ok((&target.sequence, target.stop, target.base))
    }
}

// NCBI c++/src/algo/blast/core/blast_hits.c:1330-1382:
// compare score DESC, subject.offset ASC, subject.end DESC,
// query.offset ASC, query.end DESC; sort only if an adjacent pair is out of order.
fn compare_gapped_score(a: &GappedHsp, b: &GappedHsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.s_start.cmp(&b.s_start))
        .then(b.s_end.cmp(&a.s_end))
        .then(a.q_start.cmp(&b.q_start))
        .then(b.q_end.cmp(&a.q_end))
}

fn sort_gapped_score_if_needed(hsps: &mut [GappedHsp]) {
    if hsps
        .windows(2)
        .any(|pair| compare_gapped_score(&pair[0], &pair[1]).is_gt())
    {
        hsps.sort_unstable_by(compare_gapped_score);
    }
}

// NCBI c++/src/algo/blast/core/blast_itree.c:817-838,931-995:
// query_start = s_GetQueryStrandOffset(query_info, hsp->context);
// return s_HSPIsContained(hsp, query_start, node->hsp, node->leftptr,
//                         min_diag_separation);
// NCBI c++/src/algo/blast/core/blast_traceback.c:352-360:
// Blast_IntervalTreeInit(0, query_blk->length + 1,
//   0, (subject_length > 0 ? subject_length : subject_blk->length / 3) + 1);
fn traceback_tree_hsp(hsp: &GappedHsp, query_length: i32) -> TreeHsp {
    TreeHsp {
        query_offset: hsp.q_start,
        query_end: hsp.q_end,
        subject_offset: hsp.s_start,
        subject_end: hsp.s_end,
        score: hsp.score,
        query_frame: 0,
        query_length,
        query_context_offset: 0,
        subject_frame_sign: i32::from(hsp.frame.signum()),
    }
}

// NCBI c++/src/algo/blast/core/blast_hits.c:2268-2379,2455-2537:
// purge |= (program != eBlastTypeBlastn);
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
// if (same context, query.offset, subject.offset, subject.frame) Blast_HSPFree(hsp);
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
// if (same context, query.end, subject.end, subject.frame) Blast_HSPFree(hsp);
// NCBI c++/src/algo/blast/core/blast_hits.c:1330-1382:
// BLAST_CMP(hsp2->score, hsp1->score), then subject.offset ASC,
// subject.end DESC, query.offset ASC, query.end DESC.
// This internal Stage C list has one protein query context; multi-query
// context offsets remain a separate Stage C boundary.
fn purge_traceback_common_endpoints(mut hsps: Vec<(GappedHsp, bool)>) -> Vec<(GappedHsp, bool)> {
    hsps.sort_unstable_by(|a, b| {
        a.0.q_start
            .cmp(&b.0.q_start)
            .then(a.0.s_start.cmp(&b.0.s_start))
            .then(b.0.score.cmp(&a.0.score))
            .then(b.0.q_end.cmp(&a.0.q_end))
            .then(b.0.s_end.cmp(&a.0.s_end))
    });
    let mut i = 0;
    while i + 1 < hsps.len() {
        let a = hsps[i].0;
        let b = hsps[i + 1].0;
        if a.frame == b.frame && a.q_start == b.q_start && a.s_start == b.s_start {
            hsps.remove(i + 1);
        } else {
            i += 1;
        }
    }
    hsps.sort_unstable_by(|a, b| {
        a.0.q_end
            .cmp(&b.0.q_end)
            .then(a.0.s_end.cmp(&b.0.s_end))
            .then(b.0.score.cmp(&a.0.score))
            .then(b.0.q_start.cmp(&a.0.q_start))
            .then(b.0.s_start.cmp(&a.0.s_start))
    });
    i = 0;
    while i + 1 < hsps.len() {
        let a = hsps[i].0;
        let b = hsps[i + 1].0;
        if a.frame == b.frame && a.q_end == b.q_end && a.s_end == b.s_end {
            hsps.remove(i + 1);
        } else {
            i += 1;
        }
    }
    // NCBI c++/src/algo/blast/core/blast_hits.c:1374-1382:
    // if (!Blast_HSPListIsSortedByScore(hsp_list)) qsort(..., ScoreCompareHSPs);
    if hsps
        .windows(2)
        .any(|pair| compare_gapped_score(&pair[0].0, &pair[1].0).is_gt())
    {
        hsps.sort_unstable_by(|a, b| compare_gapped_score(&a.0, &b.0));
    }
    hsps
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:303-312,408-440,503-536,
// 709-721,1644-1684:
// orig_hsplist = BlastHSPListDup(hsp_list);
// subject = Blast_HSPGetTargetTranslation(target_t, hsp, &subject_length);
// BLAST_GappedAlignmentWithTraceback(..., fence_hit);
// if (fence_error) break;
// if (fence_error) Blast_HSPListSwap(hsp_list, orig_hsplist);
// if (fence_hit) { /* refetch whole subject */
//     Blast_TracebackFromHSPList(..., &fence_hit);
// }
// The NCBI endpoint purge and interval-tree pass follow traceback; ambiguity
// reevaluation has no HSPs in this non-greedy, purge=TRUE profile.
#[allow(dead_code)]
pub(super) fn full_translation_traceback_blosum62(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    gapped: &[GappedHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
) -> Result<Vec<(GappedHsp, bool)>> {
    full_translation_traceback_with_matrix(
        query,
        subject,
        db_gencode,
        gapped,
        ScoringMatrix::Blosum62,
        gap_open,
        gap_extend,
        x_drop_final,
        0.0,
        0,
    )
}

// NCBI c++/src/algo/blast/core/blast_hits.c:993-1001:
// return ((hsp->num_ident * 100.0 < align_length * percent_identity) ||
//         align_length < hit_options->min_hit_length);
fn blast_hsp_test(
    num_ident: usize,
    align_length: usize,
    percent_identity: f64,
    min_hit_length: i32,
) -> bool {
    (num_ident as f64) * 100.0 < (align_length as f64) * percent_identity
        || (align_length as i64) < i64::from(min_hit_length)
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:508-513:
// BLAST_GappedAlignmentWithTraceback(program_number, query, adjusted_subject,
//     gap_align, score_params, q_start, s_start, query_length,
//     adjusted_s_length, fence_hit);
#[allow(dead_code)]
fn full_translation_traceback_with_matrix(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    gapped: &[GappedHsp],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    percent_identity: f64,
    min_hit_length: i32,
) -> Result<Vec<(GappedHsp, bool)>> {
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    let query_frame = encode_protein_query_frame_with_seg(query, None);
    let query_sequence = &query_frame.aa_seq[1..query_frame.aa_seq.len() - 1];
    // NCBI c++/src/algo/blast/core/blast_engine.c:539-552,840-850:
    // Blast_HSPListSortByScore(hsp_list);
    // Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
    let mut ordered = gapped.to_vec();
    sort_gapped_score_if_needed(&mut ordered);
    let query_length = i32::try_from(query.len())?;
    let subject_limit = i32::try_from(subject.len() / 3)? + 1;
    let mut target = TargetTranslation::new(subject, &code);
    for pass in 0..2 {
        let mut results = Vec::new();
        let mut retry = false;
        // NCBI c++/src/algo/blast/core/blast_traceback.c:352-360,401-405,
        // 583-605:
        // tree = Blast_IntervalTreeInit(...);
        // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info,
        //                                   hit_options->min_diag_separation)) { ... }
        // BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
        let mut tree = BlastIntervalTree::new(0, query_length + 1, 0, subject_limit);
        for hit in &ordered {
            if tree.contains_hsp(&traceback_tree_hsp(hit, query_length), 0, 0) {
                continue;
            }
            let (subject_sequence, _, subject_base) = target.get(
                hit.frame,
                if pass == 0 { hit.s_start } else { -1 },
                hit.s_end,
            )?;
            let q_start = usize::try_from(hit.q_gapped_start)?;
            let s_start = usize::try_from(hit.s_gapped_start)?
                .checked_sub(subject_base)
                .context("gapped start before translated window")?;
            let mut fence_hit = false;
            let mut scratch = GapAlignScratch::new();
            let alignment = blast_gapped_alignment_with_traceback_with_scratch(
                query_sequence,
                subject_sequence,
                q_start,
                s_start,
                matrix,
                None,
                gap_open,
                gap_extend,
                x_drop_final,
                &mut scratch,
                Some(&mut fence_hit),
            );
            if fence_hit {
                retry = true;
                break;
            }
            let alignment = alignment.context("NCBI protein traceback returned no alignment")?;
            // NCBI c++/src/algo/blast/core/blast_traceback.c:585-605:
            // Blast_HSPGetNumIdentitiesAndPositives(..., &align_length, sbp);
            // delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
            // if (!delete_hsp) BlastIntervalTreeAddHSP(...);
            // else hsp_array[index] = Blast_HSPFree(hsp);
            let align_length: usize = alignment
                .edit_script
                .iter()
                .map(|op| op.num() as usize)
                .sum();
            if blast_hsp_test(
                alignment.num_ident,
                align_length,
                percent_identity,
                min_hit_length,
            ) {
                continue;
            }
            let saved = GappedHsp {
                frame: hit.frame,
                score: alignment.score,
                q_start: i32::try_from(alignment.query_start)?,
                q_end: i32::try_from(alignment.query_stop)?,
                q_gapped_start: hit.q_gapped_start,
                s_start: i32::try_from(alignment.subject_start + subject_base)?,
                s_end: i32::try_from(alignment.subject_stop + subject_base)?,
                s_gapped_start: hit.s_gapped_start,
            };
            tree.add_hsp(
                traceback_tree_hsp(&saved, query_length),
                0,
                IndexMethod::QueryAndSubject,
            );
            results.push((saved, pass != 0));
        }
        if !retry {
            let mut results = purge_traceback_common_endpoints(results);
            // NCBI c++/src/algo/blast/core/blast_traceback.c:675-693:
            // Blast_HSPListSortByScore(hsp_list);
            // Blast_IntervalTreeReset(tree);
            // for (...) if (BlastIntervalTreeContainsHSP(tree, hsp, query_info,
            //        hit_options->min_diag_separation)) Blast_HSPFree(hsp);
            // else BlastIntervalTreeAddHSP(hsp, tree, query_info, eQueryAndSubject);
            let mut final_tree = BlastIntervalTree::new(0, query_length + 1, 0, subject_limit);
            results.retain(|(hsp, _)| {
                let tree_hsp = traceback_tree_hsp(hsp, query_length);
                if final_tree.contains_hsp(&tree_hsp, 0, 0) {
                    false
                } else {
                    final_tree.add_hsp(tree_hsp, 0, IndexMethod::QueryAndSubject);
                    true
                }
            });
            return Ok(results);
        }
        if pass == 1 {
            bail!("NCBI full-subject TBLASTN traceback reached a fence");
        }
        target = TargetTranslation::new(subject, &code);
    }
    unreachable!("NCBI traceback retries at most once")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::search_init::{
        find_blosum62_word3_init_hsps, find_blosum62_word3_init_hsps_multi,
        find_protein_init_hsps_multi,
    };
    use std::fs;

    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3924-3927:
    // cutoff = hit_params->cutoffs[context].cutoff_score;
    // Read each query context's cutoff from the comparison-only probe.
    fn multi_query_cutoffs() -> [i32; 3] {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/context_cutoffs_20260924/blosum62_word3.tsv"
        );
        let mut cutoffs = [None; 3];
        let mut count = 0;
        for line in fs::read_to_string(path).unwrap().lines() {
            let f: Vec<_> = line.split('\t').collect();
            let context: usize = f[2].parse().unwrap();
            let value: i32 = f[4].parse().unwrap();
            assert!(context < cutoffs.len());
            if let Some(previous) = cutoffs[context] {
                assert_eq!(previous, value);
            } else {
                cutoffs[context] = Some(value);
            }
            count += 1;
        }
        assert_eq!(count, 18);
        cutoffs.map(Option::unwrap)
    }

    fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
        let text = fs::read_to_string(path).unwrap();
        let mut records: Vec<(String, Vec<u8>)> = Vec::new();
        for line in text.lines() {
            if let Some(id) = line.strip_prefix('>') {
                records.push((id.to_string(), Vec::new()));
            } else {
                records.last_mut().unwrap().1.extend(line.as_bytes());
            }
        }
        records
    }

    // NCBI c++/src/algo/blast/core/blast_gapalign.c:2371-2468:
    // s_AdjustHspOffsetsAndGetQueryData converts global initial offsets to
    // local query coordinates and preserves the context on Blast_HSPInit.
    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3908-3919,4057-4091:
    // the frame-local tree checks containment before each saved gapped HSP.
    #[test]
    fn multi_query_gapped_hsps_match_ncbi_context_scores_offsets_and_order() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let initial =
            find_blosum62_word3_init_hsps_multi(&refs, subject, 1, None, 13, 40, 16, 0, false)
                .unwrap();
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("GAPPED_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                assert_eq!(f[5], "0");
                (
                    f[4].parse::<usize>().unwrap(),
                    GappedHsp {
                        frame: f[9].parse().unwrap(),
                        score: f[3].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        q_end: f[7].parse().unwrap(),
                        q_gapped_start: f[8].parse().unwrap(),
                        s_start: f[10].parse().unwrap(),
                        s_end: f[11].parse().unwrap(),
                        s_gapped_start: f[12].parse().unwrap(),
                    },
                )
            })
            .collect();
        let actual = gapped_blosum62_word3_hsps_multi(
            &refs,
            subject,
            1,
            &initial,
            11,
            1,
            38,
            &multi_query_cutoffs(),
        )
        .unwrap();
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!(
                "first multi-query gapped difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                actual.get(first),
                expected.get(first),
                actual.len(),
                expected.len()
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:539-552,840-850:
    // each frame is endpoint-purged, score-sorted, then appended to the
    // cumulative list at kHspNumMax. Compare every append snapshot.
    // NCBI c++/src/algo/blast/core/blast_hits.c:2809-2864:
    // the first list is adopted and later lists are combined by score.
    #[test]
    fn multi_query_preliminary_append_snapshots_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let initial =
            find_blosum62_word3_init_hsps_multi(&refs, subject, 1, None, 13, 40, 16, 0, false)
                .unwrap();
        let gapped = gapped_blosum62_word3_hsps_multi(
            &refs,
            subject,
            1,
            &initial,
            11,
            1,
            38,
            &multi_query_cutoffs(),
        )
        .unwrap();
        let actual = preliminary_merged_hsps(&gapped, i32::MAX as usize);
        let trace = fs::read_to_string(format!("{root}append_events.tsv")).unwrap();
        let mut expected: Vec<Vec<(usize, i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); 6];
        for line in trace
            .lines()
            .filter(|line| line.starts_with("APPEND_OUT_HSP\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            expected[f[1].parse::<usize>().unwrap()].push((
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
                f[9].parse().unwrap(),
            ));
        }
        let actual: Vec<Vec<_>> = actual
            .iter()
            .map(|snapshot| {
                snapshot
                    .iter()
                    .map(|(context, hsp)| {
                        (
                            *context,
                            hsp.score,
                            hsp.frame,
                            hsp.q_start,
                            hsp.q_end,
                            hsp.s_start,
                            hsp.s_end,
                        )
                    })
                    .collect()
            })
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(
            actual.iter().map(Vec::len).collect::<Vec<_>>(),
            [3, 5, 11, 19, 25, 31]
        );
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:259-312,635-721:
    // query-indexed preliminary lists enter traceback independently;
    // a fence retries the same original list on the full subject.
    // NCBI c++/src/algo/blast/core/blast_traceback.c:1644-1684:
    // only the second return with fence_hit=0 is retained here.
    #[test]
    fn multi_query_traceback_lists_match_ncbi_post_fence_hsps() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let initial =
            find_blosum62_word3_init_hsps_multi(&refs, subject, 1, None, 13, 40, 16, 0, false)
                .unwrap();
        let gapped = gapped_blosum62_word3_hsps_multi(
            &refs,
            subject,
            1,
            &initial,
            11,
            1,
            38,
            &multi_query_cutoffs(),
        )
        .unwrap();
        let appended = preliminary_merged_hsps(&gapped, i32::MAX as usize)
            .pop()
            .unwrap();

        let trace = fs::read_to_string(format!("{root}traceback_events.tsv")).unwrap();
        let mut expected: Vec<Vec<(i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); refs.len()];
        let mut context = 0usize;
        let mut successful_pass = false;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
                "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
                "TRACEBACK_OUT_HSP" if successful_pass => expected[context].push((
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                )),
                _ => {}
            }
        }
        for context in [1usize, 0] {
            let per_query: Vec<_> = appended
                .iter()
                .filter(|(index, _)| *index == context)
                .map(|(_, hsp)| *hsp)
                .collect();
            let actual = full_translation_traceback_blosum62(
                refs[context],
                subject,
                1,
                &per_query,
                11,
                1,
                64,
            )
            .unwrap();
            let actual: Vec<_> = actual
                .into_iter()
                .map(|(hsp, full_retry)| {
                    assert!(full_retry, "NCBI retried this list after a fence");
                    (
                        hsp.score,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    )
                })
                .collect();
            assert_eq!(actual, expected[context], "query context {context}");
        }
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:572-586:
    // adjust the new HSP's subject offsets, then merge with overlap 100.
    // NCBI c++/src/algo/blast/core/blast_hits.c:2857-3035:
    // only equal-context, close-diagonal overlap HSPs are merged.
    // The pinned oracle records old/new lists immediately before each call.
    #[test]
    fn long_chunk_overlap_merge_matches_ncbi_input_score_and_order() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/long_chunk_20260924/"
        );
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let parse_hsp = |line: &str| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse::<usize>().unwrap(),
                GappedHsp {
                    score: f[4].parse().unwrap(),
                    frame: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[9].parse().unwrap(),
                    s_end: f[10].parse().unwrap(),
                    s_gapped_start: f[11].parse().unwrap(),
                },
            )
        };
        let rows: Vec<_> = trace.lines().collect();
        let for_call = |event: &str, call: usize| {
            rows.iter()
                .copied()
                .filter(|line| {
                    let mut fields = line.split('\t');
                    fields.next() == Some(event)
                        && fields.next() == Some(if call == 0 { "0" } else { "1" })
                })
                .map(&parse_hsp)
                .collect::<Vec<_>>()
        };
        let first = merge_subject_chunk_hsps(
            Vec::new(),
            for_call("MERGE_IN_HSP", 0),
            0,
            0,
            i32::MAX as usize,
            true,
        );
        assert_eq!(first, for_call("MERGE_OUT_HSP", 0));
        let old = for_call("MERGE_OLD_HSP", 1);
        assert_eq!(first, old);
        let second = merge_subject_chunk_hsps(
            old,
            for_call("MERGE_IN_HSP", 1),
            4_999_900,
            100,
            i32::MAX as usize,
            true,
        );
        assert_eq!(second, for_call("MERGE_OUT_HSP", 1));
        assert_eq!(second[0].1.score, 660);
    }

    // NCBI c++/src/algo/blast/core/blast_hits.c:2762-2864:
    // with hsp_num_max=3, sort both lists and select by ScoreCompareHSPs;
    // when its five compared fields tie, choose the old list's HSP first.
    // The expected rows come from the direct comparison-only NCBI C API.
    #[test]
    fn preliminary_append_cap_and_tie_match_ncbi_api() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/append_cap_20260924/oracle.tsv"
        );
        let hsp = |context: usize, score: i32, frame: i8, start: i32| {
            (
                context,
                GappedHsp {
                    frame,
                    score,
                    q_start: 0,
                    q_end: 20,
                    q_gapped_start: 5,
                    s_start: start,
                    s_end: start + 20,
                    s_gapped_start: start + 5,
                },
            )
        };
        let input = vec![
            hsp(0, 100, 1, 100),
            hsp(0, 80, 1, 300),
            hsp(0, 60, 1, 500),
            hsp(1, 90, 2, 200),
            hsp(1, 80, 2, 300),
            hsp(1, 70, 2, 400),
        ];
        let snapshots = preliminary_merged_hsps(&input, 3);
        let actual: Vec<_> = snapshots[1]
            .iter()
            .enumerate()
            .map(|(index, (context, hsp))| {
                (
                    index,
                    *context,
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        let trace = fs::read_to_string(path).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("CAP_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].parse::<usize>().unwrap(),
                    f[2].parse::<usize>().unwrap(),
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i8>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                )
            })
            .collect();
        assert_eq!(actual, expected);
    }

    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953,4071-4091:
    // max_offset = BlastGetStartForGappedAlignment(...);
    // status = s_BlastProtGappedAlignment(...);
    // Blast_HSPInit(gap_align->query_start, gap_align->query_stop,
    //               gap_align->subject_start, gap_align->subject_stop, ...);
    fn compare_fixture(fixture_name: &str, oracle_name: &str, lowercase: bool) {
        let root = format!(
            "{}/../docs/evidence/tlosan_stage_c/{fixture_name}/",
            env!("CARGO_MANIFEST_DIR")
        );
        let oracle = format!(
            "{}/../docs/evidence/tlosan_stage_c/{oracle_name}/gapped_output.tsv",
            env!("CARGO_MANIFEST_DIR")
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subjects = read_fasta(&format!("{root}subjects.fna"));
        let mut actual = Vec::new();
        for (id, subject) in &subjects {
            let initial =
                find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, lowercase)
                    .unwrap();
            for hsp in
                gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38, 0).unwrap()
            {
                actual.push((id.clone(), hsp));
            }
        }
        let expected: Vec<_> = fs::read_to_string(oracle)
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].to_string(),
                    GappedHsp {
                        frame: f[2].parse().unwrap(),
                        score: f[4].parse().unwrap(),
                        q_start: f[7].parse().unwrap(),
                        q_end: f[8].parse().unwrap(),
                        q_gapped_start: f[9].parse().unwrap(),
                        s_start: f[11].parse().unwrap(),
                        s_end: f[12].parse().unwrap(),
                        s_gapped_start: f[13].parse().unwrap(),
                    },
                )
            })
            .collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn fixed_local_subject_gapped_scores_offsets_and_order_match_ncbi() {
        compare_fixture("run_20260923", "gapped_20260924", false);
    }

    #[test]
    fn lowercase_gapped_scores_offsets_and_order_match_ncbi() {
        compare_fixture("lowercase_20260923", "gapped_lowercase_20260924", true);
    }

    #[test]
    fn ambiguity_gapped_scores_offsets_and_order_match_ncbi() {
        compare_fixture("ambiguity_20260923", "gapped_ambiguity_20260924", false);
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:522-552:
    // status = aux_struct->GetGappedScore(..., init_hitlist, &hsp_list, ...);
    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3936-3953:
    // max_offset = BlastGetStartForGappedAlignment(...);
    // status = s_BlastProtGappedAlignment(...);
    #[test]
    fn code32_api_gapped_scores_offsets_and_order_match_ncbi() {
        let fixture = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures/"
        );
        let oracle = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/gapped_code32_20260924/gapped_code32_api.tsv"
        );
        let query = &read_fasta(&format!("{fixture}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{fixture}subject_code32.fna"))[0].1;
        let initial =
            find_blosum62_word3_init_hsps(query, subject, 32, None, 13, 40, 16, 13, false).unwrap();
        let actual =
            gapped_blosum62_word3_hsps(query, subject, 32, &initial, 11, 1, 38, 1).unwrap();
        let expected: Vec<GappedHsp> = fs::read_to_string(oracle)
            .unwrap()
            .lines()
            .skip(1)
            .filter_map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                if f[0] != "32" {
                    return None;
                }
                Some(GappedHsp {
                    frame: f[16].parse().unwrap(),
                    score: f[10].parse().unwrap(),
                    q_start: f[13].parse().unwrap(),
                    q_end: f[14].parse().unwrap(),
                    q_gapped_start: f[15].parse().unwrap(),
                    s_start: f[17].parse().unwrap(),
                    s_end: f[18].parse().unwrap(),
                    s_gapped_start: f[19].parse().unwrap(),
                })
            })
            .collect();
        assert_eq!(actual, expected);
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:1644-1677:
    // if (fence_hit) { /* refetch whole subject */ }
    // Blast_TracebackFromHSPList(..., &fence_hit);
    // NCBI c++/src/algo/blast/core/blast_traceback.c:503-589:
    // BLAST_GappedAlignmentWithTraceback(...);
    // Blast_HSPUpdateWithTraceback(gap_align, hsp);
    #[test]
    fn full_translation_traceback_matches_ncbi_post_fence_scores_and_offsets() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let oracle = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/gapped_20260924/traceback_hsps.tsv"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subjects = read_fasta(&format!("{root}subjects.fna"));
        let mut actual = Vec::new();
        for (id, subject) in &subjects {
            let initial =
                find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, false)
                    .unwrap();
            let gapped =
                gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38, 0).unwrap();
            for (hsp, fence_hit) in
                full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 64).unwrap()
            {
                actual.push((
                    id.clone(),
                    hsp.frame,
                    hsp.score,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    fence_hit,
                ));
            }
        }
        let expected: Vec<_> = fs::read_to_string(oracle)
            .unwrap()
            .lines()
            .skip(1)
            .filter_map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                if f[3] != "2" {
                    return None;
                }
                Some((
                    f[2].to_string(),
                    f[16].parse().unwrap(),
                    f[15].parse().unwrap(),
                    f[17].parse().unwrap(),
                    f[18].parse().unwrap(),
                    f[19].parse().unwrap(),
                    f[20].parse().unwrap(),
                    true,
                ))
            })
            .collect();
        assert_eq!(actual, expected);
        assert!(actual
            .iter()
            .any(|(id, _, score, _, _, _, _, _)| { id == "ambiguous" && *score == 646 }));
    }

    // NCBI c++/src/algo/blast/core/blast_hits.c:1154-1228:
    // start_shift = nucl_start/CODON_LENGTH;
    // target_t->range[2*context+1] = start_shift + length;
    // NCBI c++/src/algo/blast/core/blast_traceback.c:1644-1684:
    // if (fence_hit) { /* refetch whole subject and retry the HSP list */ }
    #[test]
    fn long_multi_hsp_partial_windows_and_retry_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_hsp_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let code = GeneticCode::try_from_id(1).unwrap();
        let mut target = TargetTranslation::new(subject, &code);
        for (index, line) in fs::read_to_string(format!("{root}translation_windows.tsv"))
            .unwrap()
            .lines()
            .skip(1)
            .enumerate()
        {
            if index == 1 {
                target = TargetTranslation::new(subject, &code);
            }
            let f: Vec<_> = line.split('\t').collect();
            let frame: i8 = f[1].parse().unwrap();
            let offset: i32 = f[2].parse().unwrap();
            let end: i32 = f[3].parse().unwrap();
            let (_, stop, _) = target.get(frame, offset, end).unwrap();
            let context = if frame > 0 { frame - 1 } else { 2 - frame };
            assert_eq!(
                target.frames[context as usize].start,
                f[4].parse::<usize>().unwrap()
            );
            assert_eq!(stop, f[5].parse::<usize>().unwrap());
        }
        let expected: Vec<GappedHsp> = fs::read_to_string(format!("{root}gapped_hsps.tsv"))
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                GappedHsp {
                    frame: f[9].parse().unwrap(),
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                }
            })
            .collect();
        let initial =
            find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, false).unwrap();
        let gapped = gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38, 0).unwrap();
        assert_eq!(gapped, expected);
        let actual =
            full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 64).unwrap();
        let mut pass = 0;
        let mut oracle = Vec::new();
        for line in fs::read_to_string(format!("{root}traceback_events.tsv"))
            .unwrap()
            .lines()
            .skip(1)
        {
            let f: Vec<_> = line.split('\t').collect();
            if f[0] == "TRACEBACK_INPUT" {
                pass += 1;
            }
            if pass == 2 && f[0] == "TRACEBACK_OUT_HSP" {
                oracle.push((
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i8>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                ));
            }
        }
        assert_eq!(pass, 2);
        let actual: Vec<_> = actual
            .iter()
            .map(|(hsp, retried)| {
                assert!(*retried);
                (
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        assert_eq!(actual, oracle);
    }

    // NCBI c++/src/algo/blast/core/blast_hits.c:1154-1228:
    // nucl_start = MAX(0, 3*hsp->subject.offset - kMaxTranslation);
    // return target_t->translations[context] - target_t->range[2*context] + 1;
    // NCBI c++/src/algo/blast/core/blast_traceback.c:503-536,1644-1684:
    // subject = Blast_HSPGetTargetTranslation(target_t, hsp, &subject_length);
    // BLAST_GappedAlignmentWithTraceback(..., fence_hit);
    // Only a fence causes full-subject retry.
    #[test]
    fn nonzero_base_partial_traceback_without_fence_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_hsp_no_fence_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let code = GeneticCode::try_from_id(1).unwrap();
        let mut target = TargetTranslation::new(subject, &code);
        for line in fs::read_to_string(format!("{root}translation_windows.tsv"))
            .unwrap()
            .lines()
            .skip(1)
        {
            let f: Vec<_> = line.split('\t').collect();
            let frame: i8 = f[1].parse().unwrap();
            let offset: i32 = f[2].parse().unwrap();
            let end: i32 = f[3].parse().unwrap();
            let (_, stop, base) = target.get(frame, offset, end).unwrap();
            let context = if frame > 0 { frame - 1 } else { 2 - frame };
            assert_eq!(
                target.frames[context as usize].start,
                f[4].parse::<usize>().unwrap()
            );
            assert_eq!(stop, f[5].parse::<usize>().unwrap());
            assert!(base > 0);
        }
        let initial =
            find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, false).unwrap();
        let gapped = gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 12, 0).unwrap();
        let oracle_gapped: Vec<_> = fs::read_to_string(format!("{root}gapped_hsps.tsv"))
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse::<i32>().unwrap(),
                    f[9].parse::<i8>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[10].parse::<i32>().unwrap(),
                    f[11].parse::<i32>().unwrap(),
                )
            })
            .collect();
        let actual_gapped: Vec<_> = gapped
            .iter()
            .map(|h| (h.score, h.frame, h.q_start, h.q_end, h.s_start, h.s_end))
            .collect();
        assert_eq!(actual_gapped, oracle_gapped);
        let actual =
            full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 12).unwrap();
        let oracle: Vec<_> = fs::read_to_string(format!("{root}traceback_events.tsv"))
            .unwrap()
            .lines()
            .skip(1)
            .filter(|line| line.starts_with("TRACEBACK_OUT_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i8>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                )
            })
            .collect();
        let actual: Vec<_> = actual
            .iter()
            .map(|(h, retried)| {
                assert!(!retried);
                (h.score, h.frame, h.q_start, h.q_end, h.s_start, h.s_end)
            })
            .collect();
        assert_eq!(actual, oracle);
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:804-844:
    // for (context=first_context; context<=last_context; context++) {
    //     Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
    // }
    // NCBI c++/src/algo/blast/core/blast_traceback.c:635-693:
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
    // Blast_HSPListSortByScore(hsp_list);
    #[test]
    fn six_frame_merged_traceback_endpoint_deletion_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/six_frame_merged_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let expected_gapped: Vec<GappedHsp> = fs::read_to_string(format!("{root}gapped_hsps.tsv"))
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                GappedHsp {
                    frame: f[9].parse().unwrap(),
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                }
            })
            .collect();
        let initial =
            find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, false).unwrap();
        let gapped = gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38, 0).unwrap();
        assert_eq!(gapped, expected_gapped);
        // NCBI c++/src/algo/blast/core/blast_traceback.c:503-589,635-666:
        // BLAST_GappedAlignmentWithTraceback(...);
        // Blast_HSPUpdateWithTraceback(gap_align, hsp);
        // Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
        // The extra -1 HSP grows from raw 38 to raw 236 before endpoint purge.
        let removed = gapped
            .iter()
            .find(|hsp| hsp.frame == -1 && hsp.score == 38)
            .unwrap();
        let code = GeneticCode::try_from_id(1).unwrap();
        let mut target = TargetTranslation::new(subject, &code);
        let (subject_sequence, _, _) = target.get(-1, -1, 0).unwrap();
        let encoded = encode_protein_query_frame_with_seg(query, None);
        let query_sequence = &encoded.aa_seq[1..encoded.aa_seq.len() - 1];
        let mut fence = false;
        let mut scratch = GapAlignScratch::new();
        let before_purge = blast_gapped_alignment_with_traceback_with_scratch(
            query_sequence,
            subject_sequence,
            removed.q_gapped_start as usize,
            removed.s_gapped_start as usize,
            ScoringMatrix::Blosum62,
            None,
            11,
            1,
            64,
            &mut scratch,
            Some(&mut fence),
        )
        .unwrap();
        assert!(!fence);
        let endpoint = fs::read_to_string(format!("{root}endpoint_events.tsv")).unwrap();
        let mut in_traceback_purge = false;
        let mut expected_removed = None;
        for line in endpoint.lines().skip(1) {
            if line == "ENDPOINT_INPUT\t0\t18" {
                in_traceback_purge = true;
            }
            if in_traceback_purge && line.starts_with("ENDPOINT_IN_HSP\t6\t") {
                let f: Vec<_> = line.split('\t').collect();
                expected_removed = Some((
                    f[2].parse::<i32>().unwrap(),
                    f[4].parse::<usize>().unwrap(),
                    f[5].parse::<usize>().unwrap(),
                    f[6].parse::<usize>().unwrap(),
                    f[7].parse::<usize>().unwrap(),
                ));
                break;
            }
        }
        assert_eq!(
            Some((
                before_purge.score,
                before_purge.query_start,
                before_purge.query_stop,
                before_purge.subject_start,
                before_purge.subject_stop
            )),
            expected_removed
        );
        let actual =
            full_translation_traceback_blosum62(query, subject, 1, &gapped, 11, 1, 64).unwrap();
        let mut pass = 0;
        let mut expected = Vec::new();
        for line in fs::read_to_string(format!("{root}traceback_events.tsv"))
            .unwrap()
            .lines()
            .skip(1)
        {
            let f: Vec<_> = line.split('\t').collect();
            if f[0] == "TRACEBACK_INPUT" {
                pass += 1;
            }
            if pass == 2 && f[0] == "TRACEBACK_OUT_HSP" {
                expected.push((
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i8>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                ));
            }
        }
        assert_eq!(pass, 2);
        assert_eq!(gapped.len(), 18);
        assert_eq!(expected.len(), 17);
        let actual: Vec<_> = actual
            .iter()
            .map(|(hsp, retried)| {
                assert!(*retried);
                (
                    hsp.score,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                )
            })
            .collect();
        assert_eq!(actual, expected);
    }
    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4057-4091:
    // cutoff = hit_params->cutoffs[context].cutoff_score;
    // BLAST_GetGappedScore saves each qualifying HSP in frame/chunk order.
    // NCBI c++/src/algo/blast/core/blast_engine.c:539-586,840-850:
    // purge, merge, and append follow each chunk's gapped extension.
    #[test]
    fn alternate_matrix_word2_gapped_and_append_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let initial = find_protein_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            16,
            60,
            21,
            0,
            false,
            ScoringMatrix::Blosum45,
            2,
        )
        .unwrap();
        let gapped = gapped_protein_hsps_multi_chunks(
            &refs,
            subject,
            1,
            &initial,
            14,
            2,
            53,
            &[0, 0, i32::MAX],
            ScoringMatrix::Blosum45,
        )
        .unwrap();
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("GAPPED_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                assert_eq!(f[5], "0");
                (
                    f[4].parse::<usize>().unwrap(),
                    GappedHsp {
                        score: f[3].parse().unwrap(),
                        frame: f[9].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        q_end: f[7].parse().unwrap(),
                        q_gapped_start: f[8].parse().unwrap(),
                        s_start: f[10].parse().unwrap(),
                        s_end: f[11].parse().unwrap(),
                        s_gapped_start: f[12].parse().unwrap(),
                    },
                    0usize,
                )
            })
            .collect();
        if gapped != expected {
            let first = gapped
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(gapped.len().min(expected.len()));
            panic!(
                "first BLOSUM45/word2 gapped difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                gapped.get(first), expected.get(first), gapped.len(), expected.len()
            );
        }

        let frames = generate_frames(
            &resolve_local_subject_ncbi2na(subject).unwrap(),
            &GeneticCode::try_from_id(1).unwrap(),
        );
        let frame_lengths: Vec<_> = frames
            .iter()
            .map(|frame| (frame.frame, frame.aa_len))
            .collect();
        let snapshots =
            preliminary_chunked_hsps(&gapped, &frame_lengths, i32::MAX as usize).unwrap();
        let append_trace = fs::read_to_string(format!("{root}append_events.tsv")).unwrap();
        let expected_snapshots: Vec<Vec<_>> = (0..6)
            .map(|call| {
                append_trace
                    .lines()
                    .filter(|line| line.starts_with(&format!("APPEND_OUT_HSP\t{call}\t")))
                    .map(|line| {
                        let f: Vec<_> = line.split('\t').collect();
                        (
                            f[3].parse::<usize>().unwrap(),
                            f[4].parse::<i32>().unwrap(),
                            f[5].parse::<i8>().unwrap(),
                            f[6].parse::<i32>().unwrap(),
                            f[7].parse::<i32>().unwrap(),
                            f[8].parse::<i32>().unwrap(),
                            f[9].parse::<i32>().unwrap(),
                        )
                    })
                    .collect()
            })
            .collect();
        let actual_snapshots: Vec<Vec<_>> = snapshots
            .iter()
            .map(|snapshot| {
                snapshot
                    .iter()
                    .map(|(context, hsp)| {
                        (
                            *context,
                            hsp.score,
                            hsp.frame,
                            hsp.q_start,
                            hsp.q_end,
                            hsp.s_start,
                            hsp.s_end,
                        )
                    })
                    .collect()
            })
            .collect();
        assert_eq!(actual_snapshots, expected_snapshots);
    }
    // NCBI c++/src/algo/blast/core/blast_traceback.c:401-405,583-605:
    // if (!BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...)) {
    //     BLAST_GappedAlignmentWithTraceback(...);
    // } else hsp_array[index] = Blast_HSPFree(hsp);
    // The pinned BLOSUM45/word-2 trace contains positive traceback-side
    // containment, including the full-translation retry.
    #[test]
    fn alternate_matrix_traceback_positive_containment_matches_ncbi() {
        use std::collections::HashMap;
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let append_trace = fs::read_to_string(format!("{root}append_events.tsv")).unwrap();
        let mut by_key = HashMap::new();
        let mut appended = Vec::new();
        for line in append_trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            if f[0] == "MERGE_OUT_HSP" {
                let context: usize = f[3].parse().unwrap();
                let hsp = GappedHsp {
                    score: f[4].parse().unwrap(),
                    frame: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    s_start: f[9].parse().unwrap(),
                    s_end: f[10].parse().unwrap(),
                    s_gapped_start: f[11].parse().unwrap(),
                };
                by_key.insert(
                    (
                        context,
                        hsp.score,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    ),
                    hsp,
                );
            } else if f[0] == "APPEND_OUT_HSP" && f[1] == "5" {
                appended.push((
                    f[3].parse::<usize>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[5].parse::<i8>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                    f[9].parse::<i32>().unwrap(),
                ));
            }
        }
        assert_eq!(appended.len(), 22);
        let trace = fs::read_to_string(format!("{root}traceback_events.tsv")).unwrap();
        let mut expected: Vec<Vec<(i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); queries.len()];
        let mut context = 0usize;
        let mut successful_pass = false;
        let mut in_traceback = false;
        let mut positive_containment = 0;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "TRACEBACK_INPUT" => in_traceback = true,
                "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
                "CONTAINS" if in_traceback && f[1] == "1" => positive_containment += 1,
                "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
                "TRACEBACK_OUT_HSP" if successful_pass => expected[context].push((
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                )),
                _ => {}
            }
        }
        assert!(positive_containment > 0);
        for context in [1usize, 0] {
            let per_query: Vec<_> = appended
                .iter()
                .filter(|entry| entry.0 == context)
                .map(|entry| *by_key.get(entry).unwrap())
                .collect();
            let actual = full_translation_traceback_with_matrix(
                &queries[context].1,
                subject,
                1,
                &per_query,
                ScoringMatrix::Blosum45,
                14,
                2,
                88,
                0.0,
                0,
            )
            .unwrap();
            let actual: Vec<_> = actual
                .into_iter()
                .map(|(hsp, _)| {
                    (
                        hsp.score,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    )
                })
                .collect();
            assert_eq!(actual, expected[context], "context {context}");
        }
    }
    // NCBI c++/src/algo/blast/core/blast_hits.c:993-1001,1027-1032:
    // return ((hsp->num_ident * 100.0 < align_length * percent_identity) ||
    //         align_length < hit_options->min_hit_length);
    // The direct API oracle checks equality, just-above equality, and length.
    #[test]
    fn blast_hsp_test_deletion_matches_ncbi_api() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/hsp_test_20260924/oracle.tsv"
        );
        for line in fs::read_to_string(path).unwrap().lines() {
            let f: Vec<_> = line.split('\t').collect();
            let observed = blast_hsp_test(
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
            );
            assert_eq!(observed, f[5] == "1", "{line}");
        }
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:478-552:
    // GetGappedScore receives the chunk-local translated subject before
    // Blast_HSPListAdjustOffsets and Blast_HSPListsMerge.
    #[test]
    fn long_subject_chunk_local_gapped_hsps_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let query = &read_fasta(&format!("{root}long_chunk_20260924/query.faa"))[0].1;
        let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
        let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
        let mut subject = Vec::with_capacity(15_000_962);
        for _ in 0..4_999_950 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 15_000_962);
        let initial =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 0, false).unwrap();
        let chunked = gapped_blosum62_word3_hsps_multi_chunks(
            &[query],
            &subject,
            1,
            &initial,
            11,
            1,
            38,
            &[28],
        )
        .unwrap();
        assert_eq!(
            chunked.iter().map(|row| row.2).collect::<Vec<_>>(),
            [0, 4_999_900]
        );
        let actual: Vec<_> = chunked.iter().map(|(_, hsp, _)| *hsp).collect();
        let trace =
            fs::read_to_string(format!("{root}long_chunk_20260924/gapped_events.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("GAPPED_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                GappedHsp {
                    score: f[3].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    q_end: f[7].parse().unwrap(),
                    q_gapped_start: f[8].parse().unwrap(),
                    frame: f[9].parse().unwrap(),
                    s_start: f[10].parse().unwrap(),
                    s_end: f[11].parse().unwrap(),
                    s_gapped_start: f[12].parse().unwrap(),
                }
            })
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(actual.len(), 2);
        // NCBI c++/src/algo/blast/core/blast_engine.c:572-586:
        // Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
        // Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
        //                     kHspNumMax, &(backup.offset), INT4_MIN, 100, ...);
        let mut second = actual[1];
        second.s_start += 4_999_900;
        second.s_end += 4_999_900;
        second.s_gapped_start += 4_999_900;
        let merged = merge_subject_chunk_hsps(
            vec![(0, actual[0])],
            vec![(0, second)],
            4_999_900,
            100,
            i32::MAX as usize,
            true,
        );
        let merged_expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("MERGE_OUT_HSP\t1\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse().unwrap(),
                    GappedHsp {
                        score: f[4].parse().unwrap(),
                        frame: f[5].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        q_end: f[7].parse().unwrap(),
                        q_gapped_start: f[8].parse().unwrap(),
                        s_start: f[9].parse().unwrap(),
                        s_end: f[10].parse().unwrap(),
                        s_gapped_start: f[11].parse().unwrap(),
                    },
                )
            })
            .collect();
        assert_eq!(merged, merged_expected);
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
        // per-chunk endpoint purge/offset/merge precedes one frame append.
        let frame_lengths = [
            (1, 5_000_320),
            (2, 5_000_320),
            (3, 5_000_320),
            (-1, 5_000_320),
            (-2, 5_000_320),
            (-3, 5_000_320),
        ];
        let snapshots =
            preliminary_chunked_hsps(&chunked, &frame_lengths, i32::MAX as usize).unwrap();
        assert_eq!(snapshots.iter().map(Vec::len).collect::<Vec<_>>(), [1; 6]);
        assert_eq!(snapshots.last().unwrap(), &merged_expected);
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:283-310,478-586:
    // the masked middle chunk returns SUBJECT_SPLIT_NO_RANGE; the later
    // positive-frame chunk yields one local HSP, then gets offset-adjusted.
    #[test]
    fn no_range_middle_chunk_gapped_and_merged_hsp_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let query = &read_fasta(&format!("{root}no_range_middle_20260924/query.faa"))[0].1;
        let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
        let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
        let mut subject = Vec::with_capacity(30_001_262);
        for _ in 0..10_000_050 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        subject[14_997_000..30_000_000].make_ascii_lowercase();
        let initial =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 0, true).unwrap();
        let chunked = gapped_blosum62_word3_hsps_multi_chunks(
            &[query],
            &subject,
            1,
            &initial,
            11,
            1,
            38,
            &[30],
        )
        .unwrap();
        assert_eq!(
            chunked.iter().map(|row| row.2).collect::<Vec<_>>(),
            [9_999_800]
        );
        let trace = fs::read_to_string(format!("{root}no_range_middle_20260924/gapped_events.tsv"))
            .unwrap();
        let local: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("GAPPED_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[4].parse().unwrap(),
                    GappedHsp {
                        score: f[3].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        q_end: f[7].parse().unwrap(),
                        q_gapped_start: f[8].parse().unwrap(),
                        frame: f[9].parse().unwrap(),
                        s_start: f[10].parse().unwrap(),
                        s_end: f[11].parse().unwrap(),
                        s_gapped_start: f[12].parse().unwrap(),
                    },
                    9_999_800,
                )
            })
            .collect();
        assert_eq!(chunked, local);
        let frame_lengths = [
            (1, 10_000_420),
            (2, 10_000_420),
            (3, 10_000_420),
            (-1, 10_000_420),
            (-2, 10_000_420),
            (-3, 10_000_420),
        ];
        let snapshots =
            preliminary_chunked_hsps(&chunked, &frame_lengths, i32::MAX as usize).unwrap();
        assert_eq!(snapshots.iter().map(Vec::len).collect::<Vec<_>>(), [1; 6]);
        let merged: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("MERGE_OUT_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse().unwrap(),
                    GappedHsp {
                        score: f[4].parse().unwrap(),
                        frame: f[5].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        q_end: f[7].parse().unwrap(),
                        q_gapped_start: f[8].parse().unwrap(),
                        s_start: f[9].parse().unwrap(),
                        s_end: f[10].parse().unwrap(),
                        s_gapped_start: f[11].parse().unwrap(),
                    },
                )
            })
            .collect();
        assert_eq!(snapshots.last().unwrap(), &merged);
    }

    // NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
    // skip the fully masked middle chunk; retain the first chunk's list,
    // then merge the later chunk and append that frame before the next frame.
    // NCBI c++/src/algo/blast/core/blast_hits.c:2857-3035:
    // Blast_HSPListsMerge retains both nonoverlapping equal-score HSPs.
    #[test]
    fn skipped_middle_chunk_preserves_both_hsps_and_all_append_rows() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let case = format!("{root}no_range_two_hits_20260924/");
        let query = &read_fasta(&format!("{case}query.faa"))[0].1;
        let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
        let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
        let mut subject = Vec::with_capacity(30_001_262);
        for _ in 0..4_998_450 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        subject.push(b'A');
        for _ in 0..(10_000_050 - 4_998_450 - 121) {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 30_001_262);
        subject[14_997_000..30_000_000].make_ascii_lowercase();
        let initial =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 0, true).unwrap();
        let word_trace = fs::read_to_string(format!("{case}frame_chunks.tsv")).unwrap();
        let expected_init: Vec<_> = word_trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                assert!(call < 2);
                super::super::search_init::InitHsp {
                    frame: 1,
                    chunk_offset: if call == 0 { 0 } else { 9_999_800 },
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        assert_eq!(initial, expected_init);
        let chunked = gapped_blosum62_word3_hsps_multi_chunks(
            &[query],
            &subject,
            1,
            &initial,
            11,
            1,
            38,
            &[30],
        )
        .unwrap();
        let trace = fs::read_to_string(format!("{case}gapped_events.tsv")).unwrap();
        let expected_local: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("GAPPED_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                assert!(call < 2);
                (
                    f[4].parse().unwrap(),
                    GappedHsp {
                        score: f[3].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        q_end: f[7].parse().unwrap(),
                        q_gapped_start: f[8].parse().unwrap(),
                        frame: f[9].parse().unwrap(),
                        s_start: f[10].parse().unwrap(),
                        s_end: f[11].parse().unwrap(),
                        s_gapped_start: f[12].parse().unwrap(),
                    },
                    if call == 0 { 0 } else { 9_999_800 },
                )
            })
            .collect();
        assert_eq!(chunked, expected_local);
        let frame_lengths = [
            (1, 10_000_420),
            (2, 10_000_420),
            (3, 10_000_420),
            (-1, 10_000_420),
            (-2, 10_000_420),
            (-3, 10_000_420),
        ];
        let snapshots =
            preliminary_chunked_hsps(&chunked, &frame_lengths, i32::MAX as usize).unwrap();
        let expected_append: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("APPEND_OUT_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].parse::<usize>().unwrap(),
                    f[3].parse::<usize>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[5].parse::<i8>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                    f[9].parse::<i32>().unwrap(),
                )
            })
            .collect();
        let actual_append: Vec<_> = snapshots
            .iter()
            .enumerate()
            .flat_map(|(call, hsps)| {
                hsps.iter().map(move |(context, hsp)| {
                    (
                        call,
                        *context,
                        hsp.score,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    )
                })
            })
            .collect();
        assert_eq!(actual_append, expected_append);
        assert_eq!(snapshots.iter().map(Vec::len).collect::<Vec<_>>(), [2; 6]);
    }
}
