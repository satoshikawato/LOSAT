//! Internal TBLASTN protein gapped-score boundary.
//! The public TBLASTN CLI remains unsupported until downstream stages agree.

use super::search_init::{find_protein_init_hsps_by_chunk_with_mask_mode, InitHsp, InitStageEvent};
use super::search_seed::{
    encode_tblastn_lookup_query, resolve_local_subject_ncbi2na, Seed, TranslatedChunk,
};
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::blastp::gapalign::{
    blast_gapped_alignment_with_traceback_with_scratch, blastp_get_start_for_gapped_alignment,
    blastp_score_only_gapped_alignment_with_scratch, protein_identities_from_edit_ops,
    BlastpGappedAlignmentMode, GapAlignScratch,
};
use crate::algorithm::tblastx::translation::{generate_frames, QueryFrame};
use crate::config::ScoringMatrix;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::protein_score;
use crate::utils::seg::SegParams;
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
            let chunk_hsps = score_gapped_chunk(
                queries,
                &query_frames,
                &context_offsets,
                total_query_length,
                &frame,
                chunk,
                initial,
                gap_open,
                gap_extend,
                x_drop,
                cutoff_scores,
                matrix,
                &mut scratch,
            )?;
            results.extend(
                chunk_hsps
                    .into_iter()
                    .map(|(context, hsp)| (context, hsp, chunk.offset)),
            );
        }
    }
    Ok(results)
}

// NCBI c++/src/algo/blast/core/blast_engine.c:478-552:
// BlastInitHitListReset(init_hitlist); WordFinder(..., init_hitlist, ...);
// GetGappedScore(..., init_hitlist, &hsp_list, ...);
// Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
// Score only this chunk's sorted initial HSPs before the next chunk scan.
fn score_gapped_chunk(
    queries: &[&[u8]],
    query_frames: &[QueryFrame],
    context_offsets: &[i32],
    total_query_length: i32,
    frame: &QueryFrame,
    chunk: TranslatedChunk,
    initial: &[InitHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    cutoff_scores: &[i32],
    matrix: ScoringMatrix,
    scratch: &mut GapAlignScratch,
) -> Result<Vec<(usize, GappedHsp)>> {
    let mut results = Vec::new();
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
        return Ok(Vec::new());
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
        let context = context_offsets.partition_point(|&offset| offset <= hit.q_seed as i32) - 1;
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
        let s_start = usize::try_from(hit.s_start).context("negative initial subject start")?;
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
            scratch,
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
            results.push((context, saved));
        }
    }
    Ok(results)
}

// NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
// WordFinder -> GetGappedScore -> endpoint purge -> offset adjustment ->
// Blast_HSPListsMerge runs within each chunk; Blast_HSPListAppend ends a frame.
// This diagnostic keeps that function and input order while the public CLI
// remains explicitly unimplemented through Stages C-E.
#[derive(Clone, Copy)]
struct PreliminaryProfile<'a> {
    seg: Option<&'a SegParams>,
    soft_masking: bool,
    threshold: i32,
    window: i32,
    word_xdrop: &'a [i32],
    word_cutoff: &'a [i32],
    mask_lowercase: bool,
    matrix: ScoringMatrix,
    word_size: usize,
    gap_open: i32,
    gap_extend: i32,
    gap_xdrop: i32,
    gapped_cutoff: &'a [i32],
    hsp_num_max: usize,
}

// NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
// one chunk's WordFinder and gapped result precedes the next chunk's scan;
// each frame's append precedes the next frame's scan.
#[derive(Clone, Debug, PartialEq, Eq)]
enum PreliminaryEvent {
    Candidates(i8, usize, Vec<Seed>),
    Initial(i8, usize, Vec<InitHsp>),
    Gapped(i8, usize, Vec<(usize, GappedHsp)>),
    Purged(i8, usize, Vec<(usize, GappedHsp)>),
    Merged(i8, usize, Vec<(usize, GappedHsp)>),
    Appended(i8, Vec<(usize, GappedHsp)>),
}

// NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
// s_BlastSearchEngineOneContext handles each chunk's WordFinder, gapped
// score, purge and merge before Blast_HSPListAppend for that frame.
#[allow(dead_code)]
fn preliminary_protein_hsps_in_ncbi_order(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
) -> Result<(Vec<(usize, GappedHsp)>, Vec<PreliminaryEvent>)> {
    preliminary_protein_hsps_with_comparison_input(
        queries,
        subject,
        db_gencode,
        profile,
        |_, _, _| Ok(()),
    )
}

// NCBI c++/src/algo/blast/core/blast_engine.c:525-552:
// GetGappedScore returns hsp_list; immediately afterward the same list is
// passed to Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE).
// The callback is for comparison-only artificial NCBI input states; the
// ordinary integrated path above passes an identity callback.
#[allow(dead_code)]
fn preliminary_protein_hsps_with_comparison_input(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
    mut at_purge_input: impl FnMut(i8, usize, &mut Vec<(usize, GappedHsp)>) -> Result<()>,
) -> Result<(Vec<(usize, GappedHsp)>, Vec<PreliminaryEvent>)> {
    if profile.gapped_cutoff.len() != queries.len() {
        bail!("one preliminary gapped cutoff is required per query context");
    }
    // NCBI c++/src/algo/blast/core/blast_setup.c:614-625;
    // c++/src/algo/blast/core/blast_engine.c:484-525;
    // c++/src/algo/blast/core/blast_gapalign.c:2410-2442:
    // if (!mask_at_hash) BlastSetUp_MaskQuery(query_blk, ...);
    // WordFinder(..., query, ...); GetGappedScore(..., query, ...);
    // single_query->sequence = concatenated_query->sequence + *query_start;
    // Keep the hard-SEG working query identical at both function calls.
    let query_frames: Vec<_> = queries
        .iter()
        // NCBI c++/src/algo/blast/core/blast_setup.c:614-625:
        // if (!mask_at_hash) BlastSetUp_MaskQuery(query_blk, ...);
        .map(|query| {
            if profile.soft_masking {
                encode_protein_query_frame_with_seg(query, None)
            } else {
                encode_tblastn_lookup_query(query, profile.seg, profile.mask_lowercase)
            }
        })
        .collect();
    let mut context_offsets = Vec::with_capacity(queries.len());
    let mut total_query_length = 0i32;
    for query in queries {
        context_offsets.push(total_query_length);
        total_query_length += i32::try_from(query.len())? + 1;
    }
    total_query_length -= 1;
    let mut scratch = GapAlignScratch::new();
    let mut per_frame = Vec::new();
    let mut combined = Vec::new();
    let mut events = Vec::new();
    find_protein_init_hsps_by_chunk_with_mask_mode(
        queries,
        subject,
        db_gencode,
        profile.seg,
        profile.threshold,
        profile.window,
        profile.word_xdrop,
        profile.word_cutoff,
        profile.mask_lowercase,
        profile.matrix,
        profile.word_size,
        profile.soft_masking,
        |event| {
            match event {
                InitStageEvent::Chunk(frame, chunk, seeds, initial) => {
                    events.push(PreliminaryEvent::Candidates(
                        frame.frame,
                        chunk.offset,
                        seeds.to_vec(),
                    ));
                    events.push(PreliminaryEvent::Initial(
                        frame.frame,
                        chunk.offset,
                        initial.to_vec(),
                    ));
                    // NCBI c++/src/algo/blast/core/blast_engine.c:489-493:
                    // if (init_hitlist->total == 0) continue;
                    if initial.is_empty() {
                        return Ok(());
                    }
                    let gapped = score_gapped_chunk(
                        queries,
                        &query_frames,
                        &context_offsets,
                        total_query_length,
                        frame,
                        chunk,
                        initial,
                        profile.gap_open,
                        profile.gap_extend,
                        profile.gap_xdrop,
                        profile.gapped_cutoff,
                        profile.matrix,
                        &mut scratch,
                    )?;
                    events.push(PreliminaryEvent::Gapped(
                        frame.frame,
                        chunk.offset,
                        gapped.clone(),
                    ));
                    let mut gapped = gapped;
                    at_purge_input(frame.frame, chunk.offset, &mut gapped)?;
                    // NCBI c++/src/algo/blast/core/blast_engine.c:539-552:
                    // Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE);
                    // Blast_HSPListSortByScore(hsp_list);
                    let purged = purge_preliminary_common_endpoints(gapped);
                    events.push(PreliminaryEvent::Purged(
                        frame.frame,
                        chunk.offset,
                        purged.clone(),
                    ));
                    if purged.is_empty() {
                        return Ok(());
                    }
                    per_frame = merge_adjusted_chunk(
                        std::mem::take(&mut per_frame),
                        purged,
                        chunk,
                        profile.hsp_num_max,
                    )?;
                    events.push(PreliminaryEvent::Merged(
                        frame.frame,
                        chunk.offset,
                        per_frame.clone(),
                    ));
                }
                InitStageEvent::FrameEnd(frame) => {
                    // NCBI c++/src/algo/blast/core/blast_engine.c:840-850:
                    // Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out,
                    //                     kHspNumMax);
                    combined = append_preliminary_hsps(
                        std::mem::take(&mut combined),
                        std::mem::take(&mut per_frame),
                        profile.hsp_num_max,
                    );
                    events.push(PreliminaryEvent::Appended(frame.frame, combined.clone()));
                }
            }
            Ok(())
        },
    )?;
    Ok((combined, events))
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
            per_frame = merge_adjusted_chunk(per_frame, incoming, chunk, hsp_num_max)?;
        }
        combined = append_preliminary_hsps(combined, per_frame, hsp_num_max);
        snapshots.push(combined.clone());
    }
    Ok(snapshots)
}

// NCBI c++/src/algo/blast/core/blast_engine.c:572-586:
// Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
// Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, kHspNumMax,
//                     &(backup.offset), INT4_MIN, overlap, TRUE, FALSE);
// NCBI c++/src/algo/blast/core/blast_hits.c:3037-3050:
// Adjust subject offset, end and gapped_start before merging.
fn merge_adjusted_chunk(
    per_frame: Vec<(usize, GappedHsp)>,
    mut incoming: Vec<(usize, GappedHsp)>,
    chunk: TranslatedChunk,
    hsp_num_max: usize,
) -> Result<Vec<(usize, GappedHsp)>> {
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
    Ok(merge_subject_chunk_hsps(
        per_frame,
        incoming,
        offset,
        overlap,
        hsp_num_max,
        true,
    ))
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

pub(super) struct TargetTranslation<'a> {
    subject: &'a [u8],
    code: &'a GeneticCode,
    frames: [TargetFrameTranslation; 6],
}

impl<'a> TargetTranslation<'a> {
    pub(super) fn new(subject: &'a [u8], code: &'a GeneticCode) -> Self {
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
    pub(super) fn get(
        &mut self,
        frame: i8,
        offset: i32,
        end: i32,
    ) -> Result<(&[u8], usize, usize)> {
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

// NCBI c++/src/algo/blast/core/blast_gapalign.c:3248-3321:
// if (q_length <= HSP_MAX_WINDOW) { *q_retval = q_start + q_length/2;
//     *s_retval = s_start + q_length/2; return TRUE; }
// score += sbp->matrix->data[*query_var][*subject_var];
// if (score > max_score) { max_score = score; max_offset = index1; }
// if (max_score > 0) { *q_retval = max_offset;
//     *s_retval = max_offset - q_start + s_start; return TRUE; }
// if (score > 0) { *q_retval = hsp->query.end - HSP_MAX_WINDOW/2;
//     *s_retval = hsp->subject.end - HSP_MAX_WINDOW/2; return TRUE; }
// return FALSE;
fn blast_get_offsets_for_gapped_alignment_protein(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    matrix: ScoringMatrix,
) -> Option<(usize, usize)> {
    const HSP_MAX_WINDOW: usize = 11;
    let q_length = q_end.checked_sub(q_start)?;
    let s_length = s_end.checked_sub(s_start)?;
    if q_end > query.len() || s_end > subject.len() {
        return None;
    }
    if q_length <= HSP_MAX_WINDOW {
        return Some((q_start + q_length / 2, s_start + q_length / 2));
    }
    // NCBI uses raw pointers; these are Rust bounds checks for the same reads.
    if s_start.checked_add(HSP_MAX_WINDOW)? > subject.len() || s_end < HSP_MAX_WINDOW {
        return None;
    }
    let mut score = 0i32;
    for index in 0..HSP_MAX_WINDOW {
        score += protein_score(matrix, query[q_start + index], subject[s_start + index]);
    }
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;
    for index in (q_start + HSP_MAX_WINDOW)..(q_start + q_length.min(s_length)) {
        let subject_index = s_start + index - q_start;
        score -= protein_score(
            matrix,
            query[index - HSP_MAX_WINDOW],
            subject[subject_index - HSP_MAX_WINDOW],
        );
        score += protein_score(matrix, query[index], subject[subject_index]);
        if score > max_score {
            max_score = score;
            max_offset = index;
        }
    }
    if max_score > 0 {
        return Some((max_offset, max_offset - q_start + s_start));
    }
    score = 0;
    for index in 0..HSP_MAX_WINDOW {
        score += protein_score(
            matrix,
            query[q_end - HSP_MAX_WINDOW + index],
            subject[s_end - HSP_MAX_WINDOW + index],
        );
    }
    if score > 0 {
        Some((q_end - HSP_MAX_WINDOW / 2, s_end - HSP_MAX_WINDOW / 2))
    } else {
        None
    }
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
    full_translation_traceback_with_matrix_and_events(
        query,
        subject,
        db_gencode,
        gapped,
        matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        percent_identity,
        min_hit_length,
        None,
        None,
        None,
        None,
        None,
    )
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:585-605:
// Blast_HSPUpdateWithTraceback(gap_align, hsp);
// delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
// if (delete_hsp) hsp_array[index] = Blast_HSPFree(hsp);
// This diagnostic observer records the real call order and input offsets.
#[allow(dead_code)]
fn full_translation_traceback_with_matrix_and_events(
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
    seg: Option<&SegParams>,
    test_events: Option<&mut Vec<(bool, usize, usize, usize, usize, usize)>>,
    start_events: Option<&mut Vec<(bool, i32, i32, i32, i32, i32, i32)>>,
    containment_events: Option<&mut Vec<(bool, i8, i32, i32, i32, i32)>>,
    identity_events: Option<&mut Vec<(usize, usize)>>,
) -> Result<Vec<(GappedHsp, bool)>> {
    full_translation_traceback_with_matrix_and_events_with_mask_mode(
        query,
        subject,
        db_gencode,
        gapped,
        matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        percent_identity,
        min_hit_length,
        seg,
        false,
        false,
        test_events,
        start_events,
        containment_events,
        identity_events,
    )
}

// NCBI c++/src/algo/blast/core/blast_setup.c:614-625;
// c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
// if (!mask_at_hash) BlastSetUp_MaskQuery(query_blk, ...);
// query = query_blk->sequence + context_offset;
// query_nomask = query_blk->sequence_nomask + context_offset;
// Soft masking keeps both traceback and identity inputs unmasked.
#[allow(dead_code)]
fn full_translation_traceback_with_matrix_and_events_with_mask_mode(
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
    seg: Option<&SegParams>,
    soft_masking: bool,
    mask_lowercase: bool,
    mut test_events: Option<&mut Vec<(bool, usize, usize, usize, usize, usize)>>,
    mut start_events: Option<&mut Vec<(bool, i32, i32, i32, i32, i32, i32)>>,
    mut containment_events: Option<&mut Vec<(bool, i8, i32, i32, i32, i32)>>,
    mut identity_events: Option<&mut Vec<(usize, usize)>>,
) -> Result<Vec<(GappedHsp, bool)>> {
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    // NCBI c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
    // query = query_blk->sequence + context_offset;
    // query_nomask = query_blk->sequence_nomask + context_offset;
    // Blast_HSPGetNumIdentitiesAndPositives(query_nomask, ...);
    // Trace the masked working query, then count identities on its
    // separately retained unmasked sequence.
    let query_frame = if soft_masking {
        encode_protein_query_frame_with_seg(query, None)
    } else {
        encode_tblastn_lookup_query(query, seg, mask_lowercase)
    };
    let query_sequence = &query_frame.aa_seq[1..query_frame.aa_seq.len() - 1];
    let query_nomask = query_frame
        .aa_seq_nomask
        .as_ref()
        .map_or(query_sequence, |bytes| &bytes[1..bytes.len() - 1]);
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
            // NCBI c++/src/algo/blast/core/blast_traceback.c:403-405,607-609:
            // if (... || !BlastIntervalTreeContainsHSP(tree, hsp, ...)) { ... }
            // else { hsp_array[index] = Blast_HSPFree(hsp); }
            // Test containment before
            // obtaining target translation or attempting traceback.
            let contained = tree.contains_hsp(&traceback_tree_hsp(hit, query_length), 0, 0);
            if let Some(events) = containment_events.as_deref_mut() {
                events.push((
                    contained,
                    hit.frame,
                    hit.q_start,
                    hit.q_end,
                    hit.s_start,
                    hit.s_end,
                ));
            }
            if contained {
                continue;
            }
            let (subject_sequence, _, subject_base) = target.get(
                hit.frame,
                if pass == 0 { hit.s_start } else { -1 },
                hit.s_end,
            )?;
            // NCBI c++/src/algo/blast/core/blast_traceback.c:436-449:
            // if (hsp->query.gapped_start == 0 && hsp->subject.gapped_start == 0) {
            //   retval = BlastGetOffsetsForGappedAlignment(...);
            //   if (!retval) { hsp_array[index] = Blast_HSPFree(hsp); continue; }
            // } else { q_start = hsp->query.gapped_start;
            //          s_start = hsp->subject.gapped_start; }
            let (q_start, s_start) = if hit.q_gapped_start == 0 && hit.s_gapped_start == 0 {
                let q_offset = usize::try_from(hit.q_start)?;
                let q_end = usize::try_from(hit.q_end)?;
                let s_offset = usize::try_from(hit.s_start)?
                    .checked_sub(subject_base)
                    .context("HSP subject start before translated window")?;
                let s_end = usize::try_from(hit.s_end)?
                    .checked_sub(subject_base)
                    .context("HSP subject end before translated window")?;
                let start = blast_get_offsets_for_gapped_alignment_protein(
                    query_sequence,
                    subject_sequence,
                    q_offset,
                    q_end,
                    s_offset,
                    s_end,
                    matrix,
                );
                if let Some(events) = start_events.as_deref_mut() {
                    events.push((
                        start.is_some(),
                        hit.q_start,
                        hit.q_end,
                        hit.s_start,
                        hit.s_end,
                        start.map_or(-1, |row| row.0 as i32),
                        start.map_or(-1, |row| (row.1 + subject_base) as i32),
                    ));
                }
                let Some(start) = start else { continue };
                start
            } else {
                (
                    usize::try_from(hit.q_gapped_start)?,
                    usize::try_from(hit.s_gapped_start)?
                        .checked_sub(subject_base)
                        .context("gapped start before translated window")?,
                )
            };
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
            // NCBI c++/src/algo/blast/core/blast_traceback.c:583-596:
            // Blast_HSPGetNumIdentitiesAndPositives(query_nomask,
            //     adjusted_subject, hsp, score_options, &align_length, sbp);
            // delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
            let num_ident = protein_identities_from_edit_ops(
                query_nomask,
                subject_sequence,
                alignment.query_start,
                alignment.subject_start,
                &alignment.edit_script,
                matrix,
            );
            if let Some(events) = identity_events.as_deref_mut() {
                events.push((num_ident, align_length));
            }
            let delete_hsp =
                blast_hsp_test(num_ident, align_length, percent_identity, min_hit_length);
            // NCBI blast_traceback.c:585-605 records updated HSP offsets
            // before Blast_HSPTest and before Blast_HSPAdjustSubjectOffset.
            if let Some(events) = test_events.as_deref_mut() {
                events.push((
                    delete_hsp,
                    align_length,
                    alignment.query_start,
                    alignment.query_stop,
                    alignment.subject_start + subject_base,
                    alignment.subject_stop + subject_base,
                ));
            }
            if delete_hsp {
                continue;
            }
            // NCBI c++/src/algo/blast/core/blast_traceback.c:446-448:
            // hsp->query.gapped_start = q_start;
            // hsp->subject.gapped_start = s_start;
            // The Rust subject slice is rebased; restore its global offset.
            let saved = GappedHsp {
                frame: hit.frame,
                score: alignment.score,
                q_start: i32::try_from(alignment.query_start)?,
                q_end: i32::try_from(alignment.query_stop)?,
                q_gapped_start: i32::try_from(q_start)?,
                s_start: i32::try_from(alignment.subject_start + subject_base)?,
                s_end: i32::try_from(alignment.subject_stop + subject_base)?,
                s_gapped_start: i32::try_from(s_start + subject_base)?,
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
                let contained = final_tree.contains_hsp(&tree_hsp, 0, 0);
                // NCBI c++/src/algo/blast/core/blast_traceback.c:675-693:
                // if (BlastIntervalTreeContainsHSP(tree, hsp, query_info, ...))
                //     hsp_array[index] = Blast_HSPFree(hsp);
                // After score sort, test
                // interval-tree containment again before retaining HSPs.
                if let Some(events) = containment_events.as_deref_mut() {
                    events.push((
                        contained,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    ));
                }
                if contained {
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
    use crate::algorithm::tblastn::stage_d_linking::{
        link_preliminary_hsps, reap_by_evalue, LinkedHspList,
    };
    use crate::algorithm::tblastn::stage_d_stats::{
        local_parameters_for_call, LocalParameterCall, LocalParameterOptions,
    };
    use crate::config::ProteinScoringSpec;
    use crate::stats::spouge::lookup_protein_gumbel_params;
    use crate::stats::tables::{lookup_protein_params_gapped, lookup_protein_params_ungapped};
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

    // NCBI c++/src/algo/blast/core/link_hsps.c:1765-1810:
    // Blast_HSPListSortByScore(hsp_list);
    // hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
    // compare the saved ordered NCBI link output after a Rust C-to-D call.
    fn assert_stage_d_link_list_matches_trace(list: &LinkedHspList, trace: &str, event: usize) {
        let prefix = format!("D_HSP\t{event}\tlink_after\t");
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with(&prefix))
            .collect();
        assert_eq!(list.hsps.len(), expected.len());
        for (actual, line) in list.hsps.iter().zip(expected) {
            let f: Vec<_> = line.split('\t').collect();
            assert_eq!(actual.context, f[4].parse::<usize>().unwrap());
            assert_eq!(actual.hsp.frame, f[5].parse::<i8>().unwrap());
            assert_eq!(actual.hsp.q_start, f[6].parse::<i32>().unwrap());
            assert_eq!(actual.hsp.q_end, f[7].parse::<i32>().unwrap());
            assert_eq!(actual.hsp.s_start, f[8].parse::<i32>().unwrap());
            assert_eq!(actual.hsp.s_end, f[9].parse::<i32>().unwrap());
            assert_eq!(actual.hsp.score, f[10].parse::<i32>().unwrap());
            assert_eq!(actual.num, f[12].parse::<i32>().unwrap());
            assert_eq!(
                actual.evalue.to_bits(),
                f[13].parse::<f64>().unwrap().to_bits()
            );
        }
        let list_prefix = format!("D_LIST\t{event}\tlink_after\t");
        let list_row = trace
            .lines()
            .find(|line| line.starts_with(&list_prefix))
            .unwrap();
        assert_eq!(
            list.best_evalue.to_bits(),
            list_row
                .split('\t')
                .nth(6)
                .unwrap()
                .parse::<f64>()
                .unwrap()
                .to_bits()
        );
    }

    // NCBI c++/src/algo/blast/core/blast_setup.c:964-985;
    // c++/src/algo/blast/core/blast_parameters.c:774-815,902-999,302-383;
    // c++/src/algo/blast/core/blast_engine.c:1390-1445,870-885:
    // initial hit/word/link cutoff state is computed before the preliminary
    // search; the complete appended list is passed to BLAST_LinkHsps.
    #[test]
    fn computed_default_parameters_preserve_ncbi_preliminary_link_input() {
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
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            (subject.len() / 3) as i64,
        )
        .unwrap();
        let parameters = local_parameters_for_call(
            &[(120, true), (70, true), (120, false)],
            subject.len(),
            &[gapped; 3],
            &[ungapped; 3],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: (subject.len() / 3) as i32,
                composition_based_stats: 2,
            },
        );
        let word_xdrop: Vec<_> = parameters.cutoffs.iter().map(|p| p.word_xdrop).collect();
        let word_cutoff: Vec<_> = parameters.cutoffs.iter().map(|p| p.word_cutoff).collect();
        let hit_cutoff: Vec<_> = parameters.cutoffs.iter().map(|p| p.hit_cutoff).collect();
        assert_eq!(hit_cutoff, [19, 17, i32::MAX]);
        assert_eq!(word_cutoff, [19, 17, i32::MAX]);
        assert_eq!(word_xdrop, [16, 16, 0]);
        let seg = SegParams::default();
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &word_xdrop,
            word_cutoff: &word_cutoff,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &hit_cutoff,
            hsp_num_max: i32::MAX as usize,
        };
        let (preliminary, _) =
            preliminary_protein_hsps_in_ncbi_order(&refs, subject, 1, profile).unwrap();
        let trace = fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.trace"
        ))
        .unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("D_HSP\t0\tlink_before\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[4].parse::<usize>().unwrap(),
                    f[5].parse::<i8>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                    f[9].parse::<i32>().unwrap(),
                    f[10].parse::<i32>().unwrap(),
                )
            })
            .collect();
        let actual: Vec<_> = preliminary
            .iter()
            .map(|(context, hsp)| {
                (
                    *context,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.score,
                )
            })
            .collect();
        assert_eq!(actual, expected);
        // NCBI link_hsps.c:1765-1810 and blast_engine.c:643-676:
        // the scored Stage C list is linked before preliminary E-value reap.
        let mut linked = link_preliminary_hsps(
            &preliminary,
            &[120, 70, 120],
            &parameters.lengths,
            subject.len() as i32,
            &[gapped; 3],
            &gumbel,
            &parameters.link.unwrap(),
        )
        .unwrap();
        assert_stage_d_link_list_matches_trace(&linked, &trace, 0);
        reap_by_evalue(&mut linked, parameters.prelim_evalue);
        assert_eq!(linked.hsps.len(), 20);
    }

    // NCBI c++/src/algo/blast/core/blast_setup.c:964-985;
    // c++/src/algo/blast/core/blast_engine.c:478-586,870-885:
    // the initial computed hit/word cutoffs feed the chunked preliminary search;
    // the appended output is the link input before NCBI modifies num/E-value.
    #[test]
    fn computed_positive_uneven_gap_inputs_match_ncbi_link_before() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/uneven_gap_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let gumbel = lookup_protein_gumbel_params(
            &ProteinScoringSpec {
                matrix: ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            (subject.len() / 3) as i64,
        )
        .unwrap();
        let parameters = local_parameters_for_call(
            &[(query.len(), true)],
            subject.len(),
            &[gapped],
            &[ungapped],
            LocalParameterOptions {
                expect_value: 10.0,
                do_sum_stats: true,
                max_intron_length: 0,
                gap_trigger_bits: 22.0,
                word_xdrop_bits: 7.0,
                scale_factor: 1.0,
                gumbel: Some(&gumbel),
            },
            LocalParameterCall::Initial {
                min_subject_length: (subject.len() / 3) as i32,
                composition_based_stats: 2,
            },
        );
        let word_xdrop = [parameters.cutoffs[0].word_xdrop];
        let word_cutoff = [parameters.cutoffs[0].word_cutoff];
        let gapped_cutoff = [parameters.cutoffs[0].hit_cutoff];
        assert_eq!(word_cutoff, [10]);
        let seg = SegParams::default();
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &word_xdrop,
            word_cutoff: &word_cutoff,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &gapped_cutoff,
            hsp_num_max: i32::MAX as usize,
        };
        let (preliminary, _) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let trace = fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/uneven_gap_run_20260924/uneven_gap_20260924_default.trace"
        )).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("D_HSP\t0\tlink_before\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[4].parse::<usize>().unwrap(),
                    f[5].parse::<i8>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                    f[9].parse::<i32>().unwrap(),
                    f[10].parse::<i32>().unwrap(),
                )
            })
            .collect();
        let actual: Vec<_> = preliminary
            .iter()
            .map(|(context, hsp)| {
                (
                    *context,
                    hsp.frame,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.score,
                )
            })
            .collect();
        assert_eq!(actual, expected);
        // NCBI link_hsps.c:1765-1810 and blast_engine.c:643-676:
        // the scored Stage C list is linked before preliminary E-value reap.
        let mut linked = link_preliminary_hsps(
            &preliminary,
            &[query.len() as i32],
            &parameters.lengths,
            subject.len() as i32,
            &[gapped],
            &gumbel,
            &parameters.link.unwrap(),
        )
        .unwrap();
        assert_stage_d_link_list_matches_trace(&linked, &trace, 0);
        reap_by_evalue(&mut linked, parameters.prelim_evalue);
        assert_eq!(linked.hsps.len(), 2);
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
        let initial = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
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
        let initial = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
        )
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

    // NCBI c++/src/algo/blast/core/blast_engine.c:539-552,840-850:
    // Blast_HSPListPurgeHSPsWithCommonEndpoints(..., TRUE);
    // Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
    // NCBI c++/src/algo/blast/core/blast_hits.c:2455-2537,2809-2864:
    // Remove a lower-scoring
    // common-start HSP, then cap the cumulative score-sorted list.
    #[test]
    fn real_chunk_endpoint_deletion_and_append_cap_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
        let initial = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
        )
        .unwrap();
        let mut gapped = gapped_blosum62_word3_hsps_multi(
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
        // Comparison-only NCBI input intervention, before the first
        // chunk-local endpoint purge call in ncbi_chunk_purge_cap_inject.c.
        let targets: Vec<_> = gapped
            .iter_mut()
            .filter(|(context, hsp)| {
                *context == 0
                    && hsp.frame == 1
                    && hsp.score == 20
                    && hsp.q_start == 17
                    && hsp.s_start == 1703
            })
            .collect();
        assert_eq!(targets.len(), 1);
        targets.into_iter().next().unwrap().1.q_start = 0;
        gapped
            .iter_mut()
            .find(|(context, hsp)| {
                *context == 0
                    && hsp.frame == 1
                    && hsp.score == 20
                    && hsp.q_start == 0
                    && hsp.s_start == 1703
            })
            .unwrap()
            .1
            .s_start = 200;
        let trace = fs::read_to_string(format!(
            "{root}chunk_purge_cap_real_path_20260924/trace.tsv"
        ))
        .unwrap();
        let mut expected: Vec<Vec<(usize, i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); 6];
        let mut first_input = Vec::new();
        let mut first_output = Vec::new();
        let mut purge_call = 0;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "ENDPOINT_INPUT" => purge_call += 1,
                "ENDPOINT_IN_HSP" if purge_call == 1 => first_input.push((
                    f[2].parse::<i32>().unwrap(),
                    f[3].parse::<i8>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                )),
                "ENDPOINT_OUT_HSP" if purge_call == 1 => first_output.push((
                    f[2].parse::<i32>().unwrap(),
                    f[3].parse::<i8>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                )),
                "APPEND_OUT_HSP" => expected[f[1].parse::<usize>().unwrap()].push((
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                )),
                _ => {}
            }
        }
        assert_eq!(
            first_input,
            [(656, 1, 0, 200), (393, 1, 0, 225), (20, 1, 0, 200)]
        );
        assert_eq!(first_output, [(656, 1, 0, 200), (393, 1, 0, 225)]);
        let actual = preliminary_merged_hsps(&gapped, 3);
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
        // NCBI c++/src/algo/blast/core/blast_engine.c:525-552,840-850:
        // the comparison-only C probe changes this one HSP after the first
        // GetGappedScore return and before the actual endpoint purge call;
        // all six real append calls receive hsp_num_max=3.
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16, 16, 0],
            word_cutoff: &[0, 0, i32::MAX],
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0, 0, i32::MAX],
            hsp_num_max: 3,
        };
        let mut changes = 0;
        let (_, events) = preliminary_protein_hsps_with_comparison_input(
            &refs,
            subject,
            1,
            profile,
            |frame, offset, hsps| {
                if frame == 1 && offset == 0 {
                    let matches: Vec<_> = hsps
                        .iter_mut()
                        .filter(|(context, hsp)| {
                            *context == 0
                                && hsp.score == 20
                                && hsp.q_start == 17
                                && hsp.s_start == 1703
                        })
                        .collect();
                    assert_eq!(matches.len(), 1);
                    matches.into_iter().next().unwrap().1.q_start = 0;
                    hsps.iter_mut()
                        .find(|(context, hsp)| {
                            *context == 0
                                && hsp.score == 20
                                && hsp.q_start == 0
                                && hsp.s_start == 1703
                        })
                        .unwrap()
                        .1
                        .s_start = 200;
                    changes += 1;
                }
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(changes, 1);
        let first_purged = events
            .iter()
            .find_map(|event| match event {
                PreliminaryEvent::Purged(1, 0, hsps) => Some(hsps),
                _ => None,
            })
            .unwrap();
        let actual_purged: Vec<_> = first_purged
            .iter()
            .map(|(_, hsp)| (hsp.score, hsp.frame, hsp.q_start, hsp.s_start))
            .collect();
        assert_eq!(actual_purged, first_output);
        let mut appended = Vec::new();
        for event in &events {
            if let PreliminaryEvent::Appended(_, hsps) = event {
                appended.push(
                    hsps.iter()
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
                        .collect::<Vec<_>>(),
                );
            }
        }
        assert_eq!(appended, expected);
        assert_eq!(
            appended.iter().map(Vec::len).collect::<Vec<_>>(),
            [2, 3, 3, 3, 3, 3]
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
        let initial = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
        )
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

    // NCBI c++/src/algo/blast/core/blast_traceback.c:585-605:
    // Blast_HSPUpdateWithTraceback(gap_align, hsp);
    // delete_hsp = Blast_HSPTest(hsp, hit_options, align_length);
    // if (delete_hsp) hsp_array[index] = Blast_HSPFree(hsp);
    // The comparison-only C API input sets percent_identity=100.0 during
    // each real NCBI HSPTest call, yielding positive deletions.
    #[test]
    fn real_path_hsp_test_deletion_order_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
        let initial = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
        )
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
        let trace =
            fs::read_to_string(format!("{root}hsp_test_real_path_20260924/trace.tsv")).unwrap();
        let mut expected_tests = vec![Vec::new(); refs.len()];
        let mut expected_output = vec![Vec::new(); refs.len()];
        let mut context = 0usize;
        let mut successful_pass = false;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
                "HSP_TEST" => expected_tests[context].push((
                    f[1] == "1",
                    f[2].parse::<usize>().unwrap(),
                    f[3].parse::<usize>().unwrap(),
                    f[4].parse::<usize>().unwrap(),
                    f[5].parse::<usize>().unwrap(),
                    f[6].parse::<usize>().unwrap(),
                )),
                "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
                "TRACEBACK_OUT_HSP" if successful_pass => expected_output[context].push((
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i8>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                )),
                _ => {}
            }
        }
        assert_eq!(expected_tests.iter().map(Vec::len).sum::<usize>(), 31);
        assert_eq!(
            expected_tests.iter().flatten().filter(|row| row.0).count(),
            17
        );
        for context in [1usize, 0] {
            let per_query: Vec<_> = appended
                .iter()
                .filter(|(index, _)| *index == context)
                .map(|(_, hsp)| *hsp)
                .collect();
            let mut events = Vec::new();
            let actual = full_translation_traceback_with_matrix_and_events(
                refs[context],
                subject,
                1,
                &per_query,
                ScoringMatrix::Blosum62,
                11,
                1,
                64,
                100.0,
                0,
                None,
                Some(&mut events),
                None,
                None,
                None,
            )
            .unwrap();
            assert_eq!(
                events, expected_tests[context],
                "query context {context} HSPTest order"
            );
            let actual: Vec<_> = actual
                .into_iter()
                .map(|(hsp, retried)| {
                    assert!(retried);
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
            assert_eq!(
                actual, expected_output[context],
                "query context {context} survivors"
            );
        }
    }

    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3259-3321:
    // q_length <= 11 returns the midpoint; a strictly positive sliding
    // maximum wins; otherwise the terminal 11-residue window is tested.
    #[test]
    fn protein_start_offset_window_boundaries_follow_ncbi() {
        let encode = |residues: &[u8]| {
            let frame = encode_protein_query_frame_with_seg(residues, None);
            frame.aa_seq[1..1 + residues.len()].to_vec()
        };
        let all_m = encode(&b"M".repeat(20));
        assert_eq!(
            blast_get_offsets_for_gapped_alignment_protein(
                &all_m[..11],
                &all_m[..11],
                0,
                11,
                0,
                11,
                ScoringMatrix::Blosum62,
            ),
            Some((5, 5))
        );
        assert_eq!(
            blast_get_offsets_for_gapped_alignment_protein(
                &all_m,
                &all_m,
                0,
                20,
                0,
                20,
                ScoringMatrix::Blosum62,
            ),
            Some((10, 10))
        );
        let mut end_only = b"D".repeat(9);
        end_only.extend_from_slice(&b"M".repeat(11));
        let end_only = encode(&end_only);
        assert_eq!(
            blast_get_offsets_for_gapped_alignment_protein(
                &end_only,
                &all_m[..12],
                0,
                20,
                0,
                12,
                ScoringMatrix::Blosum62,
            ),
            Some((15, 7))
        );
        let all_k = encode(&b"K".repeat(20));
        assert_eq!(
            blast_get_offsets_for_gapped_alignment_protein(
                &all_k,
                &all_m,
                0,
                20,
                0,
                20,
                ScoringMatrix::Blosum62,
            ),
            None
        );
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:436-445:
    // retval = BlastGetOffsetsForGappedAlignment(query, subject, sbp, hsp, ...);
    // if (!retval) { hsp_array[index] = Blast_HSPFree(hsp); continue; }
    // NCBI blast_gapalign.c:3259-3321: positive 11-residue window scores
    // choose a start; two nonpositive endpoint windows return FALSE.
    // NCBI blast_traceback.c:446-448: on TRUE, write q_start and s_start
    // into the HSP before alignment; the real-path positive probe tests this.
    fn compare_real_path_start_offset_with_ncbi(positive: bool) {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
        let initial = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
        )
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
        let case = if positive {
            "start_success_real_path_20260924"
        } else {
            "start_failure_real_path_20260924"
        };
        let trace = fs::read_to_string(format!("{root}{case}/trace.tsv")).unwrap();
        let mut expected_output = vec![Vec::new(); refs.len()];
        let mut expected_starts = Vec::new();
        let mut expected_after = Vec::new();
        let mut context = 0usize;
        let mut successful_pass = false;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
                "START_RESULT" => expected_starts.push((
                    f[1] == "1",
                    f[2].parse::<i32>().unwrap(),
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                )),
                "AFTER_START_HSP" => expected_after.push((
                    f[2].parse::<i32>().unwrap(),
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                )),
                "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
                "TRACEBACK_OUT_HSP" if successful_pass => expected_output[context].push((
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i8>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                )),
                _ => {}
            }
        }
        assert_eq!(expected_starts.len(), 2);
        assert!(expected_starts.iter().all(|row| row.0 == positive));
        for context in [1usize, 0] {
            let mut per_query: Vec<_> = appended
                .iter()
                .filter(|(index, _)| *index == context)
                .map(|(_, hsp)| *hsp)
                .collect();
            if context == 0 {
                // Comparison-only NCBI input-state intervention in
                // ncbi_start_failure_inject.c, before the traceback call.
                per_query[0].q_start = 0;
                per_query[0].q_end = if positive { 120 } else { 20 };
                per_query[0].q_gapped_start = 0;
                per_query[0].frame = 1;
                per_query[0].s_start = if positive { 200 } else { 0 };
                per_query[0].s_end = if positive { 320 } else { 20 };
                per_query[0].s_gapped_start = 0;
            }
            let mut starts = Vec::new();
            let actual = full_translation_traceback_with_matrix_and_events(
                refs[context],
                subject,
                1,
                &per_query,
                ScoringMatrix::Blosum62,
                11,
                1,
                64,
                0.0,
                0,
                None,
                None,
                Some(&mut starts),
                None,
                None,
            )
            .unwrap();
            if context == 0 {
                assert_eq!(starts, expected_starts);
            } else {
                assert!(starts.is_empty());
            }
            if context == 0 && positive {
                let hsp = &actual[0].0;
                let &(score, q_start, q_end, q_gap, s_start, s_end, s_gap) =
                    expected_after.last().unwrap();
                assert_eq!(
                    (
                        hsp.score,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.q_gapped_start,
                        hsp.s_start,
                        hsp.s_end,
                        hsp.s_gapped_start
                    ),
                    (score, q_start, q_end, q_gap, s_start, s_end, s_gap),
                    "NCBI successful acquired start writeback"
                );
            }
            let actual: Vec<_> = actual
                .into_iter()
                .map(|(hsp, retried)| {
                    assert!(retried);
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
            assert_eq!(actual, expected_output[context], "query context {context}");
        }
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:436-445:
    // if (!retval) { hsp_array[index] = Blast_HSPFree(hsp); continue; }
    #[test]
    fn real_path_start_offset_failure_deletion_matches_ncbi() {
        compare_real_path_start_offset_with_ncbi(false);
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:436-448:
    // if (retval) { hsp->query.gapped_start = q_start;
    //               hsp->subject.gapped_start = s_start; }
    #[test]
    fn real_path_start_offset_success_writeback_matches_ncbi() {
        compare_real_path_start_offset_with_ncbi(true);
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
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
        // same function order and code-32 subject translation enter the
        // comparison-only NCBI C++ API gapped score oracle.
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[13],
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[1],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, _) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 32, profile).unwrap();
        assert_eq!(
            integrated,
            expected.into_iter().map(|hsp| (0, hsp)).collect::<Vec<_>>()
        );
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
            &[21, 21, 0],
            &[0, 0, i32::MAX],
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
        let mut appended = Vec::new();
        let mut merged_by_key = HashMap::new();
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
                merged_by_key.insert(
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
        // NCBI c++/src/algo/blast/core/blast_engine.c:804-850:
        // The same natural query/subject
        // inputs pass through WordFinder, gapped extension, chunk purge,
        // merge and append before this positive containment call.
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let initial = find_protein_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            16,
            60,
            &[21, 21, 0],
            &[0, 0, i32::MAX],
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
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
        // preserve the same BLOSUM45/word-2 chunk-to-append call path before
        // testing positive traceback interval-tree containment.
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 16,
            window: 60,
            word_xdrop: &[21, 21, 0],
            word_cutoff: &[0, 0, i32::MAX],
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum45,
            word_size: 2,
            gap_open: 14,
            gap_extend: 2,
            gap_xdrop: 53,
            gapped_cutoff: &[0, 0, i32::MAX],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, _) =
            preliminary_protein_hsps_in_ncbi_order(&refs, subject, 1, profile).unwrap();
        assert_eq!(integrated, *snapshots.last().unwrap());
        let complete_path = &integrated;
        let actual_appended: Vec<_> = complete_path
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
            .collect();
        assert_eq!(actual_appended, appended);
        // NCBI c++/src/algo/blast/core/blast_engine.c:572-586,840-850;
        // c++/src/algo/blast/core/blast_hits.c:2809-2864:
        // append keeps each merged HSP's gapped starts. Check the exact
        // traceback input state, which APPEND_OUT_HSP does not print.
        for (context, hsp) in complete_path {
            let key = (
                *context,
                hsp.score,
                hsp.frame,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
            );
            assert_eq!(*hsp, *merged_by_key.get(&key).unwrap());
        }
        let trace = fs::read_to_string(format!("{root}traceback_events.tsv")).unwrap();
        let mut expected: Vec<Vec<(i32, i8, i32, i32, i32, i32)>> = vec![Vec::new(); queries.len()];
        let mut expected_containment = vec![Vec::new(); queries.len()];
        let mut context = 0usize;
        let mut successful_pass = false;
        let mut in_traceback = false;
        let mut positive_containment = 0;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "TRACEBACK_INPUT" => in_traceback = true,
                "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
                "CONTAINS" if in_traceback => {
                    let result = f[1] == "1";
                    positive_containment += usize::from(result);
                    expected_containment[context].push((
                        result,
                        f[2].parse::<i8>().unwrap(),
                        f[3].parse::<i32>().unwrap(),
                        f[4].parse::<i32>().unwrap(),
                        f[5].parse::<i32>().unwrap(),
                        f[6].parse::<i32>().unwrap(),
                    ));
                }
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
            let per_query: Vec<_> = complete_path
                .iter()
                .filter(|(index, _)| *index == context)
                .map(|(_, hsp)| *hsp)
                .collect();
            let mut containment_events = Vec::new();
            let actual = full_translation_traceback_with_matrix_and_events(
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
                None,
                None,
                None,
                Some(&mut containment_events),
                None,
            )
            .unwrap();
            assert_eq!(
                containment_events, expected_containment[context],
                "NCBI containment predicate and call order, context {context}"
            );
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
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 28.
        let initial =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 28, false).unwrap();
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
    // NCBI c++/src/algo/blast/core/blast_gapalign.c:3924-3927:
    // each long chunk uses its query-context cutoff during gapped scoring.
    // NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
    // Both chunks merge, then the
    // frame list is appended in score order for both valid query contexts.
    #[test]
    fn long_multi_query_chunk_hsps_and_append_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let fixture = format!("{root}long_multi_query_20260924/");
        let queries = read_fasta(&format!("{fixture}query.faa"));
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let plus1 = read_fasta(&format!("{root}run_20260923/subjects.fna"))
            .into_iter()
            .find(|(id, _)| id == "plus1")
            .unwrap()
            .1;
        let mut subject = Vec::with_capacity(15_000_962);
        for _ in 0..4_999_950 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(&plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 15_000_962);
        let cutoffs = fs::read_to_string(format!("{fixture}context_cutoffs.tsv")).unwrap();
        let gapped_cutoff_rows: Vec<_> = cutoffs
            .lines()
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                assert_eq!(f[0], "GAPPED_CONTEXT_CUTOFF");
                (
                    f[1].parse::<usize>().unwrap(),
                    f[2].parse::<usize>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                )
            })
            .collect();
        assert_eq!(gapped_cutoff_rows.len(), 6);
        for call in 0..2 {
            for (context, value) in [28, 25, i32::MAX].into_iter().enumerate() {
                assert_eq!(
                    gapped_cutoff_rows[call * 3 + context],
                    (call, context, value)
                );
            }
        }
        let expected_cutoffs = [28, 25, i32::MAX];
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:547-583:
        // cutoffs = word_params->cutoffs + curr_context;
        // the saved probe records 28, 25, INT4_MAX in all 12 chunk calls.
        let word_trace = fs::read_to_string(format!("{fixture}word_context_cutoffs.tsv")).unwrap();
        let word_rows: Vec<_> = word_trace
            .lines()
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                assert_eq!(f[0], "WORD_CONTEXT_CUTOFF");
                (
                    f[1].parse::<usize>().unwrap(),
                    f[2].parse::<usize>().unwrap(),
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                )
            })
            .collect();
        assert_eq!(word_rows.len(), 36);
        for call in 0..12 {
            for (context, (init, drop, cutoff)) in [(16, 16, 28), (16, 16, 25), (0, 0, i32::MAX)]
                .into_iter()
                .enumerate()
            {
                assert_eq!(
                    word_rows[call * 3 + context],
                    (call, context, init, drop, cutoff)
                );
            }
        }
        let word_cutoffs = [28, 25, i32::MAX];
        let initial = find_blosum62_word3_init_hsps_multi(
            &refs,
            &subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &word_cutoffs,
            false,
        )
        .unwrap();
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:200-234,547-614:
        // after per-context extension/cutoff, sorted chunk-local init HSPs
        // enter GetGappedScore; compare all rows before any gapped deletion.
        let wordfinder_trace = fs::read_to_string(format!("{fixture}frame_chunks.tsv")).unwrap();
        let expected_initial: Vec<_> = wordfinder_trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                InitHsp {
                    frame: [1, 2, 3, -1, -2, -3][call / 2],
                    chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        assert_eq!(initial, expected_initial);
        assert_eq!(initial.len(), 4);
        let chunked = gapped_blosum62_word3_hsps_multi_chunks(
            &refs,
            &subject,
            1,
            &initial,
            11,
            1,
            38,
            &expected_cutoffs,
        )
        .unwrap();
        let trace = fs::read_to_string(format!("{fixture}gapped_events.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("GAPPED_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                (
                    f[4].parse::<usize>().unwrap(),
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
                    if call == 0 { 0 } else { 4_999_900 },
                )
            })
            .collect();
        assert_eq!(chunked, expected);
        assert_eq!(chunked.len(), 4);
        let frames = generate_frames(
            &resolve_local_subject_ncbi2na(&subject).unwrap(),
            &GeneticCode::try_from_id(1).unwrap(),
        );
        let frame_lengths: Vec<_> = frames
            .iter()
            .map(|frame| (frame.frame, frame.aa_len))
            .collect();
        let snapshots =
            preliminary_chunked_hsps(&chunked, &frame_lengths, i32::MAX as usize).unwrap();
        let expected_snapshots: Vec<Vec<_>> = (0..6)
            .map(|call| {
                trace
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
            .map(|list| {
                list.iter()
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
        // NCBI c++/src/algo/blast/core/blast_engine.c:572-586,840-850:
        // The final merged HSPs
        // retain exact gapped starts through append into traceback.
        let expected_merged: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("MERGE_OUT_HSP\t1\t"))
            .map(|line| {
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
            })
            .collect();
        assert_eq!(snapshots.last().unwrap(), &expected_merged);
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
        // each chunk's scansub -> WordFinder -> GetGappedScore -> purge ->
        // offset adjustment -> merge precedes the next chunk; append ends
        // each frame. Compare this one-call Rust path with every saved NCBI
        // function input and output, including all 394 candidate pairs.
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16, 16, 0],
            word_cutoff: &word_cutoffs,
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &expected_cutoffs,
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&refs, &subject, 1, profile).unwrap();
        assert_eq!(integrated, expected_merged);
        let frame_order = [1, 2, 3, -1, -2, -3];
        let candidate_trace = fs::read_to_string(format!("{fixture}candidate.stderr")).unwrap();
        let mut candidate_by_call = vec![Vec::new(); 12];
        for line in candidate_trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            candidate_by_call[call].push(Seed {
                frame: frame_order[call / 2],
                chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            });
        }
        assert_eq!(candidate_by_call.iter().map(Vec::len).sum::<usize>(), 394);
        let mut tags = Vec::new();
        let mut call = 0usize;
        let mut append = 0usize;
        for event in &events {
            match event {
                PreliminaryEvent::Candidates(frame, offset, seeds) => {
                    tags.push("C");
                    assert_eq!(
                        (*frame, *offset),
                        (
                            frame_order[call / 2],
                            if call % 2 == 0 { 0 } else { 4_999_900 }
                        )
                    );
                    assert_eq!(seeds, &candidate_by_call[call], "candidate call {call}");
                    call += 1;
                }
                PreliminaryEvent::Initial(frame, offset, hsps) => {
                    tags.push("I");
                    let expected: Vec<_> = expected_initial
                        .iter()
                        .copied()
                        .filter(|hit| hit.frame == *frame && hit.chunk_offset as usize == *offset)
                        .collect();
                    assert_eq!(*hsps, expected, "WordFinder call {}", call - 1);
                }
                PreliminaryEvent::Gapped(frame, offset, hsps) => {
                    tags.push("G");
                    let expected: Vec<_> = expected
                        .iter()
                        .filter(|(_, hsp, chunk_offset)| {
                            hsp.frame == *frame && chunk_offset == offset
                        })
                        .map(|(context, hsp, _)| (*context, *hsp))
                        .collect();
                    assert_eq!(*hsps, expected, "GetGappedScore call {}", call - 1);
                }
                PreliminaryEvent::Purged(frame, offset, hsps) => {
                    tags.push("P");
                    let expected: Vec<_> = expected
                        .iter()
                        .filter(|(_, hsp, chunk_offset)| {
                            hsp.frame == *frame && chunk_offset == offset
                        })
                        .map(|(context, hsp, _)| (*context, *hsp))
                        .collect();
                    assert_eq!(*hsps, expected, "natural endpoint purge call {}", call - 1);
                }
                PreliminaryEvent::Merged(frame, offset, hsps) => {
                    tags.push("M");
                    let merge_call = if *offset == 0 { 0 } else { 1 };
                    assert_eq!(*frame, 1);
                    let expected: Vec<_> = trace
                        .lines()
                        .filter(|line| line.starts_with(&format!("MERGE_OUT_HSP\t{merge_call}\t")))
                        .map(|line| {
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
                        })
                        .collect();
                    assert_eq!(*hsps, expected, "merge call {merge_call}");
                }
                PreliminaryEvent::Appended(frame, hsps) => {
                    tags.push("A");
                    assert_eq!(*frame, frame_order[append]);
                    let actual: Vec<_> = hsps
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
                        .collect();
                    assert_eq!(actual, expected_snapshots[append], "append call {append}");
                    append += 1;
                }
            }
        }
        let mut expected_tags = Vec::new();
        for call in 0..12 {
            expected_tags.extend(["C", "I"]);
            if call < 2 {
                expected_tags.extend(["G", "P", "M"]);
            }
            if call % 2 == 1 {
                expected_tags.push("A");
            }
        }
        assert_eq!(tags, expected_tags, "NCBI function execution order");
        assert_eq!((call, append), (12, 6));
        // NCBI c++/src/algo/blast/core/blast_traceback.c:259-312,401-449,583-721:
        // Process each
        // query-indexed appended list; a fence retries on full translation.
        let mut expected_final = vec![Vec::new(); refs.len()];
        let mut context = 0usize;
        let mut successful_pass = false;
        for line in trace.lines() {
            let f: Vec<_> = line.split('\t').collect();
            match f[0] {
                "TRACEBACK_CONTEXT" => context = f[1].parse().unwrap(),
                "TRACEBACK_OUTPUT" => successful_pass = f[4] == "0",
                "TRACEBACK_OUT_HSP" if successful_pass => expected_final[context].push((
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i8>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                )),
                _ => {}
            }
        }
        for context in [1usize, 0] {
            // NCBI c++/src/algo/blast/core/blast_traceback.c:259-312:
            // each query-indexed HSP list enters traceback from the completed
            // frame append stream; the saved trace visits context 1 then 0.
            let per_query: Vec<_> = integrated
                .iter()
                .filter(|(index, _)| *index == context)
                .map(|(_, hsp)| *hsp)
                .collect();
            let actual = full_translation_traceback_blosum62(
                refs[context],
                &subject,
                1,
                &per_query,
                11,
                1,
                64,
            )
            .unwrap();
            let actual: Vec<_> = actual
                .into_iter()
                .map(|(hsp, retry)| {
                    assert!(retry);
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
            assert_eq!(actual, expected_final[context], "query context {context}");
        }
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
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 30.
        let initial =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 30, true).unwrap();
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
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
        // SUBJECT_SPLIT_NO_RANGE skips WordFinder/GetGappedScore in the
        // masked middle chunk; append still occurs once for every frame.
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[30],
            mask_lowercase: true,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[30],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], &subject, 1, profile).unwrap();
        assert_eq!(integrated, merged);
        assert_eq!(
            events
                .iter()
                .filter(|event| matches!(event, PreliminaryEvent::Appended(..)))
                .count(),
            6
        );
        assert_eq!(
            events
                .iter()
                .filter(|event| matches!(event, PreliminaryEvent::Gapped(..)))
                .count(),
            1
        );
    }

    // NCBI c++/src/algo/blast/core/blast_setup.c:614-625;
    // c++/src/algo/blast/core/blast_engine.c:484-525;
    // c++/src/algo/blast/core/blast_gapalign.c:2410-2442:
    // the hard-SEG working query is shared by WordFinder and GetGappedScore;
    // the unmasked suffix still produces a real HSP and traceback.
    #[test]
    fn hard_seg_query_bytes_and_integrated_hsp_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/seg_hard_query_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let seg = SegParams::default();
        let encoded = encode_protein_query_frame_with_seg(query, Some(&seg));
        assert!(!encoded.seg_masks.is_empty());
        assert!(encoded.aa_seq[1..41].iter().all(|&b| b == 21));
        let hex: String = encoded.aa_seq[1..1 + query.len()]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
        let word: Vec<_> = state
            .lines()
            .filter(|line| line.starts_with("WORD_QUERY_BYTES\t"))
            .collect();
        let gap: Vec<_> = state
            .lines()
            .filter(|line| line.starts_with("GAPPED_QUERY_BYTES\t"))
            .collect();
        assert_eq!((word.len(), gap.len()), (6, 1));
        for (call, row) in word.iter().enumerate() {
            assert_eq!(*row, format!("WORD_QUERY_BYTES\t{call}\t160\t{hex}"));
        }
        assert_eq!(gap[0], format!("GAPPED_QUERY_BYTES\t0\t160\t{hex}"));
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let mut candidates = vec![Vec::new(); 6];
        for line in candidate_trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            candidates[call].push(Seed {
                frame: [1, 2, 3, -1, -2, -3][call],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            });
        }
        assert_eq!(candidates.iter().map(Vec::len).sum::<usize>(), 185);
        let gapped_trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let fields: Vec<_> = gapped_trace
            .lines()
            .find(|line| line.starts_with("GAPPED_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        let ncbi_hsp = GappedHsp {
            score: fields[3].parse().unwrap(),
            frame: fields[9].parse().unwrap(),
            q_start: fields[6].parse().unwrap(),
            q_end: fields[7].parse().unwrap(),
            q_gapped_start: fields[8].parse().unwrap(),
            s_start: fields[10].parse().unwrap(),
            s_end: fields[11].parse().unwrap(),
            s_gapped_start: fields[12].parse().unwrap(),
        };
        assert_eq!(integrated, vec![(0, ncbi_hsp)]);
        let mut call = 0usize;
        let mut append = 0usize;
        for event in &events {
            match event {
                PreliminaryEvent::Candidates(frame, offset, seeds) => {
                    assert_eq!(
                        (*frame, *offset, seeds),
                        ([1, 2, 3, -1, -2, -3][call], 0, &candidates[call])
                    );
                    call += 1;
                }
                PreliminaryEvent::Initial(frame, _, hsps) => {
                    assert_eq!(*frame, [1, 2, 3, -1, -2, -3][call - 1]);
                    assert_eq!(hsps.len(), usize::from(call == 1));
                    if call == 1 {
                        assert_eq!(
                            (
                                hsps[0].q_seed,
                                hsps[0].s_seed,
                                hsps[0].q_start,
                                hsps[0].s_start,
                                hsps[0].length,
                                hsps[0].score
                            ),
                            (43, 3, 40, 0, 120, 656)
                        );
                    }
                }
                PreliminaryEvent::Gapped(_, _, hsps) | PreliminaryEvent::Purged(_, _, hsps) => {
                    assert_eq!(*hsps, integrated)
                }
                PreliminaryEvent::Merged(_, _, hsps) => assert_eq!(*hsps, integrated),
                PreliminaryEvent::Appended(frame, hsps) => {
                    assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                    assert_eq!(*hsps, integrated);
                    append += 1;
                }
            }
        }
        assert_eq!((call, append), (6, 6));
        let traceback =
            full_translation_traceback_blosum62(query, subject, 1, &[integrated[0].1], 11, 1, 64)
                .unwrap();
        assert_eq!(traceback, vec![(ncbi_hsp, true)]);
    }

    // NCBI c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
    // query_blk->sequence (masked) drives the traceback score, while
    // query_blk->sequence_nomask drives identities before Blast_HSPTest.
    // The retained NCBI call has a real HSP spanning masked residues 60..78.
    #[test]
    fn internal_seg_crossing_traceback_uses_masked_score_and_unmasked_identity() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/seg_cross_traceback_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let seg = SegParams::default();
        let encoded = encode_protein_query_frame_with_seg(query, Some(&seg));
        assert_eq!(query.len(), 138);
        assert!(encoded.aa_seq[61..79].iter().all(|&b| b == 21));
        let unmasked = encoded.aa_seq_nomask.as_ref().unwrap();
        assert!(unmasked[61..79].iter().all(|&b| b == 10));
        let masked_hex: String = encoded.aa_seq[1..139]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let unmasked_hex: String = unmasked[1..139]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
        let word: Vec<_> = state
            .lines()
            .filter(|line| line.starts_with("WORD_QUERY_BYTES\t"))
            .collect();
        let gap: Vec<_> = state
            .lines()
            .filter(|line| line.starts_with("GAPPED_QUERY_BYTES\t"))
            .collect();
        let identity: Vec<_> = state
            .lines()
            .filter(|line| line.starts_with("IDENTITY_STATE\t"))
            .collect();
        assert_eq!((word.len(), gap.len(), identity.len()), (6, 1, 1));
        for (call, row) in word.iter().enumerate() {
            assert_eq!(*row, format!("WORD_QUERY_BYTES\t{call}\t138\t{masked_hex}"));
        }
        assert_eq!(gap[0], format!("GAPPED_QUERY_BYTES\t0\t138\t{masked_hex}"));
        assert_eq!(
            identity[0],
            format!("IDENTITY_STATE\t0\t638\t138\t138\t0\t138\t{unmasked_hex}")
        );
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let f: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("GAPPED_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected = GappedHsp {
            score: f[3].parse().unwrap(),
            frame: f[9].parse().unwrap(),
            q_start: f[6].parse().unwrap(),
            q_end: f[7].parse().unwrap(),
            q_gapped_start: f[8].parse().unwrap(),
            s_start: f[10].parse().unwrap(),
            s_end: f[11].parse().unwrap(),
            s_gapped_start: f[12].parse().unwrap(),
        };
        assert_eq!(expected.score, 638);
        assert_eq!((expected.q_start, expected.q_end), (0, 138));
        assert_eq!(integrated, vec![(0, expected)]);
        assert_eq!(
            events
                .iter()
                .filter(|event| matches!(event, PreliminaryEvent::Appended(..)))
                .count(),
            6
        );
        let mut identities = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events(
            query,
            subject,
            1,
            &[integrated[0].1],
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            0,
            Some(&seg),
            None,
            None,
            None,
            Some(&mut identities),
        )
        .unwrap();
        assert_eq!(actual, vec![(expected, true)]);
        assert_eq!(identities, [(138, 138)]);
        let unmasked_result = full_translation_traceback_with_matrix(
            query,
            subject,
            1,
            &[integrated[0].1],
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            0,
        )
        .unwrap();
        assert_eq!(unmasked_result[0].0.score, 746);
        assert_ne!(unmasked_result[0].0.score, expected.score);
    }

    // NCBI c++/src/algo/blast/core/blast_setup.c:614-638;
    // c++/src/algo/blast/core/blast_engine.c:484-525;
    // c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
    // if (!mask_at_hash) BlastSetUp_MaskQuery(query_blk, ...);
    // BLAST_ComplementMaskLocations(..., filter_maskloc, lookup_segments);
    // Soft masking excludes the same seed intervals but keeps K in every
    // WordFinder, GetGappedScore and traceback query input.
    #[test]
    fn internal_soft_seg_crossing_matches_ncbi_query_hsp_and_traceback() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/seg_soft_traceback_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let seg = SegParams::default();
        let unmasked = encode_protein_query_frame_with_seg(query, None);
        let hex: String = unmasked.aa_seq[1..139]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
        let word: Vec<_> = state
            .lines()
            .filter(|line| line.starts_with("WORD_QUERY_BYTES\t"))
            .collect();
        assert_eq!(word.len(), 6);
        for (call, row) in word.iter().enumerate() {
            assert_eq!(*row, format!("WORD_QUERY_BYTES\t{call}\t138\t{hex}"));
        }
        assert!(state
            .lines()
            .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t138\t{hex}")));
        assert!(state
            .lines()
            .any(|row| row == format!("IDENTITY_STATE\t0\t746\t138\t138\t0\t138\t{hex}")));
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: true,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: false,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let mut candidates = vec![Vec::new(); 6];
        for line in candidate_trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            candidates[call].push(Seed {
                frame: [1, 2, 3, -1, -2, -3][call],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            });
        }
        assert_eq!(candidates.iter().map(Vec::len).sum::<usize>(), 183);
        let mut call = 0;
        let mut append = 0;
        for event in &events {
            match event {
                PreliminaryEvent::Candidates(frame, offset, seeds) => {
                    assert_eq!(
                        (*frame, *offset, seeds),
                        ([1, 2, 3, -1, -2, -3][call], 0, &candidates[call])
                    );
                    call += 1;
                }
                PreliminaryEvent::Initial(_, _, hsps) => {
                    assert_eq!(hsps.len(), usize::from(call == 1));
                    if call == 1 {
                        assert_eq!(
                            (
                                hsps[0].q_seed,
                                hsps[0].s_seed,
                                hsps[0].q_start,
                                hsps[0].s_start,
                                hsps[0].length,
                                hsps[0].score
                            ),
                            (3, 3, 0, 0, 138, 746)
                        );
                    }
                }
                PreliminaryEvent::Appended(frame, hsps) => {
                    assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                    assert_eq!(*hsps, integrated);
                    append += 1;
                }
                _ => {}
            }
        }
        assert_eq!((call, append), (6, 6));
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let f: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("GAPPED_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected = GappedHsp {
            score: f[3].parse().unwrap(),
            frame: f[9].parse().unwrap(),
            q_start: f[6].parse().unwrap(),
            q_end: f[7].parse().unwrap(),
            q_gapped_start: f[8].parse().unwrap(),
            s_start: f[10].parse().unwrap(),
            s_end: f[11].parse().unwrap(),
            s_gapped_start: f[12].parse().unwrap(),
        };
        assert_eq!(expected.score, 746);
        assert_eq!(integrated, vec![(0, expected)]);
        let mut identities = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
            query,
            subject,
            1,
            &[expected],
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            0,
            Some(&seg),
            true,
            false,
            None,
            None,
            None,
            Some(&mut identities),
        )
        .unwrap();
        assert_eq!(actual, vec![(expected, true)]);
        assert_eq!(identities, [(138, 138)]);
    }

    // NCBI c++/src/algo/blast/blastinput/blast_fasta_input.cpp:486-502;
    // c++/src/algo/blast/core/blast_filter.c:1241-1255;
    // c++/src/algo/blast/core/blast_setup.c:614-638:
    // -lcase_masking includes query lowercase spans in the filter locations.
    // The hard-mask working query has X in positions 40..60; identities use
    // the original sequence after the traceback score is computed.
    #[test]
    fn lowercase_query_mask_matches_ncbi_candidates_hsps_and_traceback() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/lcase_query_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let masked = encode_tblastn_lookup_query(query, None, true);
        assert_eq!(query.len(), 120);
        assert!(masked.aa_seq[41..61].iter().all(|&residue| residue == 21));
        let nomask = masked.aa_seq_nomask.as_ref().unwrap();
        let masked_hex: String = masked.aa_seq[1..121]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let nomask_hex: String = nomask[1..121]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
        for call in 0..6 {
            assert!(state
                .lines()
                .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t120\t{masked_hex}")));
        }
        assert!(state
            .lines()
            .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t120\t{masked_hex}")));
        assert!(state
            .lines()
            .any(|row| row == format!("IDENTITY_STATE\t0\t528\t120\t120\t0\t120\t{nomask_hex}")));
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: true,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let mut candidates = vec![Vec::new(); 6];
        for line in candidate_trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            candidates[call].push(Seed {
                frame: [1, 2, 3, -1, -2, -3][call],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            });
        }
        assert_eq!(candidates.iter().map(Vec::len).sum::<usize>(), 154);
        let initial_trace = fs::read_to_string(format!("{root}wordfinder.stderr")).unwrap();
        let initial: Vec<_> = initial_trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[3].parse::<u32>().unwrap(),
                    f[4].parse::<u32>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                )
            })
            .collect();
        assert_eq!(initial, [(63, 63, 60, 60, 60, 325), (3, 3, 0, 0, 40, 223)]);
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let f: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("GAPPED_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected = GappedHsp {
            score: f[3].parse().unwrap(),
            frame: f[9].parse().unwrap(),
            q_start: f[6].parse().unwrap(),
            q_end: f[7].parse().unwrap(),
            q_gapped_start: f[8].parse().unwrap(),
            s_start: f[10].parse().unwrap(),
            s_end: f[11].parse().unwrap(),
            s_gapped_start: f[12].parse().unwrap(),
        };
        assert_eq!(expected.score, 528);
        assert_eq!(integrated, vec![(0, expected)]);
        let mut call = 0;
        let mut append = 0;
        for event in &events {
            match event {
                PreliminaryEvent::Candidates(frame, offset, seeds) => {
                    assert_eq!(
                        (*frame, *offset, seeds),
                        ([1, 2, 3, -1, -2, -3][call], 0, &candidates[call])
                    );
                    call += 1;
                }
                PreliminaryEvent::Initial(_, _, hsps) => {
                    assert_eq!(
                        hsps.iter()
                            .map(|h| (h.q_seed, h.s_seed, h.q_start, h.s_start, h.length, h.score))
                            .collect::<Vec<_>>(),
                        if call == 1 {
                            initial.clone()
                        } else {
                            Vec::new()
                        }
                    );
                }
                PreliminaryEvent::Appended(frame, hsps) => {
                    assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                    assert_eq!(*hsps, integrated);
                    append += 1;
                }
                _ => {}
            }
        }
        assert_eq!((call, append), (6, 6));
        let mut identities = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
            query,
            subject,
            1,
            &[expected],
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            0,
            None,
            false,
            true,
            None,
            None,
            None,
            Some(&mut identities),
        )
        .unwrap();
        assert_eq!(actual, vec![(expected, true)]);
        assert_eq!(identities, [(120, 120)]);
    }

    // NCBI c++/src/algo/blast/core/blast_setup.c:614-638;
    // c++/src/algo/blast/core/blast_filter.c:1241-1255:
    // mask_at_hash leaves the lowercase query unmasked during extension and
    // traceback, while its complement still restricts lookup positions.
    #[test]
    fn lowercase_query_soft_mask_matches_ncbi_lookup_and_traceback() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/lcase_soft_query_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        // NCBI blast_setup.c:614-638 keeps this query unmasked when
        // mask_at_hash is true; all six WordFinder calls see the same bytes.
        let unmasked = encode_protein_query_frame_with_seg(query, None);
        let hex: String = unmasked.aa_seq[1..121]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
        for call in 0..6 {
            assert!(state
                .lines()
                .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t120\t{hex}")));
        }
        assert!(state
            .lines()
            .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t120\t{hex}")));
        assert!(state
            .lines()
            .any(|row| row == format!("IDENTITY_STATE\t0\t656\t120\t120\t0\t120\t{hex}")));
        let profile = PreliminaryProfile {
            seg: None,
            soft_masking: true,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: true,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let reference = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let f: Vec<_> = reference
            .lines()
            .find(|line| line.starts_with("GAPPED_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected = GappedHsp {
            score: f[3].parse().unwrap(),
            frame: f[9].parse().unwrap(),
            q_start: f[6].parse().unwrap(),
            q_end: f[7].parse().unwrap(),
            q_gapped_start: f[8].parse().unwrap(),
            s_start: f[10].parse().unwrap(),
            s_end: f[11].parse().unwrap(),
            s_gapped_start: f[12].parse().unwrap(),
        };
        assert_eq!(expected.score, 656);
        assert_eq!(integrated, vec![(0, expected)]);
        let hard = fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/lcase_query_20260924/candidate.stderr"
        ))
        .unwrap();
        assert_eq!(
            fs::read_to_string(format!("{root}candidate.stderr")).unwrap(),
            hard
        );
        let mut saved_candidates = vec![Vec::new(); 6];
        for line in hard.lines().filter(|line| line.starts_with("CAND\t")) {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            saved_candidates[call].push(Seed {
                frame: [1, 2, 3, -1, -2, -3][call],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            });
        }
        let mut candidate_call = 0;
        let mut initial_call = 0;
        let mut append_call = 0;
        for event in &events {
            match event {
                PreliminaryEvent::Candidates(frame, offset, seeds) => {
                    assert_eq!(
                        (*frame, *offset, seeds),
                        (
                            [1, 2, 3, -1, -2, -3][candidate_call],
                            0,
                            &saved_candidates[candidate_call]
                        )
                    );
                    candidate_call += 1;
                }
                PreliminaryEvent::Initial(_, _, hsps) => {
                    assert_eq!(hsps.len(), usize::from(initial_call == 0));
                    if initial_call == 0 {
                        assert_eq!(
                            (
                                hsps[0].q_seed,
                                hsps[0].s_seed,
                                hsps[0].q_start,
                                hsps[0].s_start,
                                hsps[0].length,
                                hsps[0].score
                            ),
                            (3, 3, 0, 0, 120, 656)
                        );
                    }
                    initial_call += 1;
                }
                PreliminaryEvent::Appended(frame, hsps) => {
                    assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append_call]);
                    assert_eq!(*hsps, integrated);
                    append_call += 1;
                }
                _ => {}
            }
        }
        assert_eq!((candidate_call, initial_call, append_call), (6, 6, 6));
        let mut identities = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
            query,
            subject,
            1,
            &[expected],
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            0,
            None,
            true,
            true,
            None,
            None,
            None,
            Some(&mut identities),
        )
        .unwrap();
        assert_eq!(actual, vec![(expected, true)]);
        assert_eq!(identities, [(120, 120)]);
    }

    // NCBI c++/src/algo/blast/core/blast_filter.c:1241-1255;
    // c++/src/algo/blast/core/blast_setup.c:614-638:
    // BlastSeqLocAppend(filter_out, lcase_mask_slp);
    // BlastSeqLocCombine(filter_out, 0);
    // The SEG prefix and overlapping lowercase span combine into [0,45).
    #[test]
    fn overlapping_seg_and_lowercase_query_mask_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/seg_lcase_overlap_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let seg = SegParams::default();
        let masked = encode_tblastn_lookup_query(query, Some(&seg), true);
        assert_eq!(masked.seg_masks, [(0, 45)]);
        assert!(masked.aa_seq[1..46].iter().all(|&residue| residue == 21));
        let hex: String = masked.aa_seq[1..161]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
        for call in 0..6 {
            assert!(state
                .lines()
                .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t160\t{hex}")));
        }
        assert!(state
            .lines()
            .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t160\t{hex}")));
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: false,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: true,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let candidate_trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let mut saved_candidates = vec![Vec::new(); 6];
        for line in candidate_trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
        {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            saved_candidates[call].push(Seed {
                frame: [1, 2, 3, -1, -2, -3][call],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            });
        }
        assert_eq!(saved_candidates.iter().map(Vec::len).sum::<usize>(), 168);
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let f: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("GAPPED_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected = GappedHsp {
            score: f[3].parse().unwrap(),
            frame: f[9].parse().unwrap(),
            q_start: f[6].parse().unwrap(),
            q_end: f[7].parse().unwrap(),
            q_gapped_start: f[8].parse().unwrap(),
            s_start: f[10].parse().unwrap(),
            s_end: f[11].parse().unwrap(),
            s_gapped_start: f[12].parse().unwrap(),
        };
        assert_eq!(expected.score, 623);
        assert_eq!(integrated, vec![(0, expected)]);
        let mut call = 0;
        let mut append = 0;
        for event in &events {
            match event {
                PreliminaryEvent::Candidates(frame, offset, seeds) => {
                    assert_eq!(
                        (*frame, *offset, seeds),
                        ([1, 2, 3, -1, -2, -3][call], 0, &saved_candidates[call])
                    );
                    call += 1;
                }
                PreliminaryEvent::Initial(_, _, hsps) => {
                    assert_eq!(hsps.len(), usize::from(call == 1));
                    if call == 1 {
                        assert_eq!(
                            (
                                hsps[0].q_seed,
                                hsps[0].s_seed,
                                hsps[0].q_start,
                                hsps[0].s_start,
                                hsps[0].length,
                                hsps[0].score
                            ),
                            (48, 8, 45, 5, 115, 623)
                        );
                    }
                }
                PreliminaryEvent::Appended(frame, hsps) => {
                    assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                    assert_eq!(*hsps, integrated);
                    append += 1;
                }
                _ => {}
            }
        }
        assert_eq!((call, append), (6, 6));
        let mut identities = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
            query,
            subject,
            1,
            &[expected],
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            0,
            Some(&seg),
            false,
            true,
            None,
            None,
            None,
            Some(&mut identities),
        )
        .unwrap();
        assert_eq!(actual, vec![(expected, true)]);
        assert_eq!(identities, [(115, 115)]);
    }

    // NCBI c++/src/algo/blast/core/blast_filter.c:1241-1255;
    // c++/src/algo/blast/core/blast_setup.c:614-638:
    // BlastSeqLocCombine(filter_out, 0);
    // if (!mask_at_hash) BlastSetUp_MaskQuery(query_blk, ...);
    // The overlapping SEG/lowercase intervals still restrict lookup, while
    // WordFinder, GetGappedScore and traceback use the unmasked query.
    #[test]
    fn overlapping_seg_lowercase_soft_mask_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/seg_lcase_overlap_soft_20260924/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let seg = SegParams::default();
        let unmasked = encode_protein_query_frame_with_seg(query, None);
        let hex: String = unmasked.aa_seq[1..161]
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let state = fs::read_to_string(format!("{root}seg_state.stderr")).unwrap();
        for call in 0..6 {
            assert!(state
                .lines()
                .any(|row| row == format!("WORD_QUERY_BYTES\t{call}\t160\t{hex}")));
        }
        assert!(state
            .lines()
            .any(|row| row == format!("GAPPED_QUERY_BYTES\t0\t160\t{hex}")));
        assert!(state
            .lines()
            .any(|row| row == format!("IDENTITY_STATE\t0\t656\t120\t120\t40\t160\t{hex}")));
        let profile = PreliminaryProfile {
            seg: Some(&seg),
            soft_masking: true,
            threshold: 13,
            window: 40,
            word_xdrop: &[16],
            word_cutoff: &[0],
            mask_lowercase: true,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            gap_open: 11,
            gap_extend: 1,
            gap_xdrop: 38,
            gapped_cutoff: &[0],
            hsp_num_max: i32::MAX as usize,
        };
        let (integrated, events) =
            preliminary_protein_hsps_in_ncbi_order(&[query], subject, 1, profile).unwrap();
        let trace = fs::read_to_string(format!("{root}gapped_events.tsv")).unwrap();
        let f: Vec<_> = trace
            .lines()
            .find(|line| line.starts_with("GAPPED_HSP\t"))
            .unwrap()
            .split('\t')
            .collect();
        let expected = GappedHsp {
            score: f[3].parse().unwrap(),
            frame: f[9].parse().unwrap(),
            q_start: f[6].parse().unwrap(),
            q_end: f[7].parse().unwrap(),
            q_gapped_start: f[8].parse().unwrap(),
            s_start: f[10].parse().unwrap(),
            s_end: f[11].parse().unwrap(),
            s_gapped_start: f[12].parse().unwrap(),
        };
        assert_eq!(expected.score, 656);
        assert_eq!(integrated, vec![(0, expected)]);
        let candidates = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let hard_candidates = fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/seg_lcase_overlap_20260924/candidate.stderr"
        ))
        .unwrap();
        assert_eq!(candidates, hard_candidates);
        let mut saved_candidates = vec![Vec::new(); 6];
        for line in candidates.lines().filter(|line| line.starts_with("CAND\t")) {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            saved_candidates[call].push(Seed {
                frame: [1, 2, 3, -1, -2, -3][call],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            });
        }
        assert_eq!(saved_candidates.iter().map(Vec::len).sum::<usize>(), 168);
        let mut call = 0;
        let mut append = 0;
        for event in &events {
            match event {
                PreliminaryEvent::Candidates(frame, offset, seeds) => {
                    assert_eq!(
                        (*frame, *offset, seeds),
                        ([1, 2, 3, -1, -2, -3][call], 0, &saved_candidates[call])
                    );
                    call += 1;
                }
                PreliminaryEvent::Initial(_, _, hsps) => {
                    assert_eq!(hsps.len(), usize::from(call == 1));
                    if call == 1 {
                        assert_eq!(
                            (
                                hsps[0].q_seed,
                                hsps[0].s_seed,
                                hsps[0].q_start,
                                hsps[0].s_start,
                                hsps[0].length,
                                hsps[0].score
                            ),
                            (48, 8, 40, 0, 120, 656)
                        );
                    }
                }
                PreliminaryEvent::Appended(frame, hsps) => {
                    assert_eq!(*frame, [1, 2, 3, -1, -2, -3][append]);
                    assert_eq!(*hsps, integrated);
                    append += 1;
                }
                _ => {}
            }
        }
        assert_eq!((call, append), (6, 6));
        let mut identities = Vec::new();
        let actual = full_translation_traceback_with_matrix_and_events_with_mask_mode(
            query,
            subject,
            1,
            &[expected],
            ScoringMatrix::Blosum62,
            11,
            1,
            64,
            0.0,
            0,
            Some(&seg),
            true,
            true,
            None,
            None,
            None,
            Some(&mut identities),
        )
        .unwrap();
        assert_eq!(actual, vec![(expected, true)]);
        assert_eq!(identities, [(120, 120)]);
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
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 30.
        let initial =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 30, true).unwrap();
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
