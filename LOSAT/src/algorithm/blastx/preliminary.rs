//! Context-local BLASTX preliminary HSPs before linking/traceback.
use super::{
    args::ResolvedOptions, parameters::ContextParameters, query_setup::PreparedQueryBatch,
    seed::InitHsp,
};
use crate::algorithm::blastn::interval_tree::{
    x_grid_cells, BlastIntervalTree, IndexMethod, TreeHsp,
};
use crate::algorithm::blastp::gapalign::{
    blastp_get_start_for_gapped_alignment, blastp_score_only_gapped_alignment_with_scratch,
    gapped_start_int4, gapped_start_window, BlastpGappedAlignmentMode, GapAlignScratch,
};
use crate::config::ScoringMatrix;
use anyhow::Result;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-103
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:126-143
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PreliminaryHsp {
    pub context: usize,
    pub frame: i8,
    pub score: i32,
    pub q_start: i32,
    pub q_end: i32,
    pub q_gapped_start: i32,
    pub s_start: i32,
    pub s_end: i32,
    pub s_gapped_start: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-226
fn tree_hsp(h: &PreliminaryHsp, batch: &PreparedQueryBatch) -> TreeHsp {
    let strand_context = h.context / 6 * 6 + if h.context % 6 < 3 { 0 } else { 3 };
    TreeHsp {
        query_offset: h.q_start,
        query_end: h.q_end,
        subject_offset: h.s_start,
        subject_end: h.s_end,
        score: h.score,
        query_frame: h.frame as i32,
        query_length: batch.contexts[h.context].length as i32,
        query_context_offset: batch.contexts[strand_context].offset as i32,
        subject_frame_sign: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1635-1675
fn heapify(hits: &mut [PreliminaryHsp], mut base: usize) {
    while hits.len() >= 2 && base <= (hits.len() - 2) / 2 {
        let left = 2 * base + 1;
        let worst = if left == hits.len() - 1 || compare_score(&hits[left], &hits[left + 1]).is_ge()
        {
            left
        } else {
            left + 1
        };
        if compare_score(&hits[base], &hits[worst]).is_lt() {
            hits.swap(base, worst);
            base = worst;
        } else {
            break;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1687-1706
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1751-1809
fn save_hsp(hits: &mut Vec<PreliminaryHsp>, h: PreliminaryHsp, max: usize, heap: &mut bool) {
    if hits.len() < max {
        hits.push(h);
        return;
    }
    if !*heap {
        for i in (0..hits.len() / 2).rev() {
            heapify(hits, i);
        }
        *heap = true;
    }
    if compare_score(&h, &hits[0]).is_gt() {
        return;
    }
    hits[0] = h;
    heapify(hits, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3860-3869
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3900-3953
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4058-4087
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3909-3921
pub(crate) fn gapped_observed(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    initial: &[InitHsp],
    xdrop: i32,
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<PreliminaryHsp>> {
    gapped_observed_with(
        super::seed::Blosum62Bx,
        batch,
        parameters,
        options,
        subject,
        initial,
        xdrop,
        observe,
    )
}

/// [`gapped_observed`] with another matrix than BLOSUM62 (item SC6 of line NP). Not inlined, so that
/// BLOSUM62's keeps its code.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn gapped_observed_of<S: super::seed::BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    initial: &[InitHsp],
    xdrop: i32,
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<PreliminaryHsp>> {
    gapped_observed_with(
        scoring, batch, parameters, options, subject, initial, xdrop, observe,
    )
}

/// The body of [`gapped_observed`], with the matrix of `scoring`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn gapped_observed_with<S: super::seed::BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    initial: &[InitHsp],
    xdrop: i32,
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<PreliminaryHsp>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3707-3708
    // NCBI returns at once when the initial hit list is empty (and the engine skips the gapped
    // stage, blast_engine.c:500). This port returns the empty list before it builds the interval
    // tree and the DP scratch; the loop below would have done nothing. Most (chunk, subject)
    // pairs have no initial HSP.
    if initial.is_empty() {
        return Ok(Vec::new());
    }
    let last = batch.contexts.last().expect("BLASTX contexts");
    let (q_max, s_max) = ((last.offset + last.length + 1) as i32, subject.len() as i32);
    let mut tree = BlastIntervalTree::new(0, q_max, 0, s_max);
    tree.reset_with_bounds_cells(0, q_max, 0, s_max, x_grid_cells(initial.len()));
    let mut scratch = GapAlignScratch::new();
    let mut hits = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:213-227
    // hsp_num_max defaults to zero; max_hsps_per_subject is a later collector filter.
    let max = i32::MAX as usize;
    let mut heap = false;
    for init in initial {
        let context = batch
            .contexts
            .partition_point(|c| c.offset as i32 <= init.q_seed)
            .saturating_sub(1);
        let c = &batch.contexts[context];
        let mut h = PreliminaryHsp {
            context,
            frame: c.frame,
            score: init.score,
            q_start: init.q_start - c.offset as i32,
            q_end: init.q_start - c.offset as i32 + init.length,
            q_gapped_start: init.q_seed - c.offset as i32,
            s_start: init.s_start,
            s_end: init.s_start + init.length,
            s_gapped_start: init.s_seed,
        };
        let candidate = tree_hsp(&h, batch);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3909-3921
        let contained = tree.contains_hsp(&candidate, candidate.query_context_offset, 0);
        observe(&h, contained);
        if contained {
            continue;
        }
        let query = &batch.sequence_start[1 + c.offset..1 + c.offset + c.length];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3935-3942
        // A negative length comes only from the one-hit extension (`-window_size 0`,
        // aa_ungapped.c:1054,1083), which NCBI passes as `Uint4` (`gapped_start_negative`).
        let start = if init.length < 0 {
            gapped_start_negative(
                query,
                &subject[..subject.len() - 1],
                h.q_start as usize,
                h.s_start as usize,
                scoring.matrix(),
            )
        } else {
            blastp_get_start_for_gapped_alignment(
                query,
                &subject[..subject.len() - 1],
                h.q_start as usize,
                init.length as usize,
                h.s_start as usize,
                init.length as usize,
                scoring.matrix(),
            )
        };
        h.s_gapped_start += start as i32 - h.q_gapped_start;
        h.q_gapped_start = start as i32;
        let a = blastp_score_only_gapped_alignment_with_scratch(
            query,
            &subject[..subject.len() - 1],
            start,
            h.s_gapped_start as usize,
            scoring.matrix(),
            options.gap_open,
            options.gap_extend,
            xdrop,
            BlastpGappedAlignmentMode::Exact,
            &mut scratch,
        );
        if a.score < parameters[context].hit_cutoff {
            continue;
        }
        h.score = a.score;
        h.q_start = a.query_start;
        h.q_end = a.query_stop;
        h.s_start = a.subject_start;
        h.s_end = a.subject_stop;
        let node = tree_hsp(&h, batch);
        save_hsp(&mut hits, h, max, &mut heap);
        tree.add_hsp(
            node,
            node.query_context_offset,
            IndexMethod::QueryAndSubject,
        );
    }
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1054,1083 and
// c++/src/algo/blast/core/blast_gapalign.c:3394,3401,3405-3413,3419-3420
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501,885-893,
// c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1802 and seqdbatlas.cpp:114-127
/// The gapped start of a one-hit HSP of negative length at `q_start` of its context's letters
/// (`query`) and at `s_start` of the subject's (`subject`). NCBI reads the length as `Uint4`, so it
/// scores only the first window of `HSP_MAX_WINDOW` letters (`gapped_start_int4`). A window that
/// passes the end of its context or of the subject holds the NULLB after it (`BLAST_SCORE_MIN`), so
/// its sum is below 0 and the start is `q_start` whatever bytes follow: inside NCBI's buffers the
/// next contexts' letters, the next database sequences or the zeros of the mapped file's last page.
/// Past those buffers NCBI reads memory that is not its own (the byte after the query batch's
/// `calloc` block of `buflen` plus 1 bytes, after a `-subject` protein's `malloc` block of `len`
/// plus 2 bytes, after the last memory page of a database volume's sequence file) and can crash;
/// LOSAT reads every letter past the sequence as the sentinel (`gapped_start_window`;
/// PD-LOSAT-NCBI-DEFECTS exception 3, which version 1.48 extends from BLASTP to BLASTX). Out of
/// line, so that the gapped stage of the other HSPs keeps its code.
#[cold]
#[inline(never)]
fn gapped_start_negative(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    matrix: ScoringMatrix,
) -> usize {
    gapped_start_int4(
        &gapped_start_window(query, q_start),
        &gapped_start_window(subject, s_start),
        q_start,
        matrix,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4745-4776
pub fn ungapped(batch: &PreparedQueryBatch, initial: &[InitHsp]) -> Vec<PreliminaryHsp> {
    let mut hits = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:213-227
    // hsp_num_max defaults to zero; max_hsps_per_subject is a later collector filter.
    let max = i32::MAX as usize;
    let mut heap = false;
    for init in initial {
        let context = batch
            .contexts
            .partition_point(|c| c.offset as i32 <= init.q_seed)
            .saturating_sub(1);
        let c = &batch.contexts[context];
        save_hsp(
            &mut hits,
            PreliminaryHsp {
                context,
                frame: c.frame,
                score: init.score,
                q_start: init.q_start - c.offset as i32,
                q_end: init.q_start - c.offset as i32 + init.length,
                q_gapped_start: init.q_seed - c.offset as i32,
                s_start: init.s_start,
                s_end: init.s_start + init.length,
                s_gapped_start: init.s_seed,
            },
            max,
            &mut heap,
        );
    }
    sort_by_score(&mut hits);
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1358
pub fn compare_score(a: &PreliminaryHsp, b: &PreliminaryHsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.s_start.cmp(&b.s_start))
        .then(b.s_end.cmp(&a.s_end))
        .then(a.q_start.cmp(&b.q_start))
        .then(b.q_end.cmp(&a.q_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1383
pub fn sort_by_score(hits: &mut [PreliminaryHsp]) {
    hits.sort_by(compare_score);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2295-2321
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2362-2387
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2505
pub fn purge_endpoints(hits: &mut Vec<PreliminaryHsp>) {
    hits.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then(a.q_start.cmp(&b.q_start))
            .then(a.s_start.cmp(&b.s_start))
            .then(b.score.cmp(&a.score))
            .then(b.q_end.cmp(&a.q_end))
            .then(b.s_end.cmp(&a.s_end))
    });
    hits.dedup_by(|b, a| {
        a.context == b.context && a.q_start == b.q_start && a.s_start == b.s_start
    });
    hits.sort_by(|a, b| {
        a.context
            .cmp(&b.context)
            .then(a.q_end.cmp(&b.q_end))
            .then(a.s_end.cmp(&b.s_end))
            .then(b.score.cmp(&a.score))
            .then(b.q_start.cmp(&a.q_start))
            .then(b.s_start.cmp(&a.s_start))
    });
    hits.dedup_by(|b, a| a.context == b.context && a.q_end == b.q_end && a.s_end == b.s_end);
    sort_by_score(hits);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-226
#[cfg(test)]
#[path = "preliminary_tests.rs"]
mod tests;
