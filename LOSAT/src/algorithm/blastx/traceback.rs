//! Ordinary BLASTX composition-0 traceback in context-local amino-acid coordinates.
use super::{
    args::ResolvedOptions,
    kappa::RedoneHsp,
    parameters::ContextParameters,
    preliminary::{purge_endpoints, sort_by_score, PreliminaryHsp},
    query_setup::PreparedQueryBatch,
};
use crate::{
    algorithm::{
        blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp},
        blastp::gapalign::{
            adjust_subject_range, blast_gapped_alignment_with_traceback_with_scratch,
            GapAlignScratch,
        },
    },
    config::{ProteinScoringSpec, ScoringMatrix},
    stats::lookup_protein_params,
};
use anyhow::{ensure, Result};
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:367-403
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-407
pub(crate) fn ordinary_traceback_observed(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    _parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    observe: &mut dyn FnMut(&PreliminaryHsp, bool),
) -> Result<Vec<RedoneHsp>> {
    ensure!(
        options.gapped && options.composition == 0,
        "BLASTX ordinary traceback requires gapped composition 0"
    );
    let ka = lookup_protein_params(&ProteinScoringSpec {
        matrix: options.scoring,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    // blast_parameters.c:454-463: the final X-drop (`parameters::gapped_x_drops`).
    let xdrop = super::parameters::gapped_x_drops(options, &ka).1;
    let mut tree = BlastIntervalTree::new(
        0,
        batch.sequence_start.len() as i32,
        0,
        subject.len() as i32 + 1,
    );
    let mut scratch = GapAlignScratch::new();
    let mut output = Vec::new();
    for original in input {
        let mut h = original.clone();
        let node = tree_hsp(&h, batch);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:402-407
        let contained = tree.contains_hsp(&node, node.query_context_offset, 0);
        observe(&h, contained);
        if contained {
            continue;
        }
        let c = &batch.contexts[h.context];
        let query = &batch.sequence_start[1 + c.offset..1 + c.offset + c.length];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:437-475
        if h.q_gapped_start == 0 && h.s_gapped_start == 0 {
            let Some((q, s)) =
                get_offsets_for_gapped_alignment(query, subject, &h, options.scoring)
            else {
                continue;
            };
            h.q_gapped_start = q;
            h.s_gapped_start = s;
        }
        let (s_start, s_length, shift) = adjust_subject_range(
            h.s_gapped_start as usize,
            subject.len(),
            h.q_gapped_start as usize,
            c.length,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:479-513
        let Some(aligned) = blast_gapped_alignment_with_traceback_with_scratch(
            query,
            &subject[shift..shift + s_length],
            h.q_gapped_start as usize,
            s_start,
            options.scoring,
            None,
            options.gap_open,
            options.gap_extend,
            xdrop,
            &mut scratch,
            None,
        ) else {
            continue;
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:581-605
        h.score = aligned.score;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583,599
        h.q_start = aligned.query_start;
        h.q_end = aligned.query_stop;
        h.s_start = aligned.subject_start + i32::try_from(shift)?;
        h.s_end = aligned.subject_stop + i32::try_from(shift)?;
        let node = tree_hsp(&h, batch);
        tree.add_hsp(
            node,
            node.query_context_offset,
            IndexMethod::QueryAndSubject,
        );
        output.push(RedoneHsp {
            hsp: h,
            edit_script: aligned.edit_script,
            composition_method: 0,
            identity_matrix: None,
        });
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2467
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:634-693
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:672-693
    // Preserve each edit script with its original HSP while the source comparator
    // and endpoint purge rearrange the HSP values. Duplicate keys retain input order.
    let key = |h: &PreliminaryHsp| {
        (
            h.context,
            h.frame,
            h.score,
            h.q_start,
            h.q_end,
            h.q_gapped_start,
            h.s_start,
            h.s_end,
            h.s_gapped_start,
        )
    };
    let mut payloads = std::collections::HashMap::<_, std::collections::VecDeque<_>>::new();
    for p in &output {
        payloads
            .entry(key(&p.hsp))
            .or_default()
            .push_back(p.clone());
    }
    let mut hsps: Vec<_> = output.iter().map(|h| h.hsp.clone()).collect();
    purge_endpoints(&mut hsps);
    sort_by_score(&mut hsps);
    tree.reset();
    let mut retained = Vec::new();
    for h in hsps {
        let node = tree_hsp(&h, batch);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:677-690
        let contained = tree.contains_hsp(&node, node.query_context_offset, 0);
        observe(&h, contained);
        if contained {
            continue;
        }
        tree.add_hsp(
            node,
            node.query_context_offset,
            IndexMethod::QueryAndSubject,
        );
        retained.push(
            payloads
                .get_mut(&key(&h))
                .and_then(|queue| queue.pop_front())
                .expect("traceback payload remains attached"),
        );
    }
    Ok(retained)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3330
fn get_offsets_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    h: &PreliminaryHsp,
    matrix: ScoringMatrix,
) -> Option<(i32, i32)> {
    let q_start = h.q_start as usize;
    let s_start = h.s_start as usize;
    let q_length = (h.q_end - h.q_start) as usize;
    let s_length = (h.s_end - h.s_start) as usize;
    if q_length <= 11 {
        return Some((
            h.q_start + (q_length / 2) as i32,
            h.s_start + (q_length / 2) as i32,
        ));
    }
    let score_at =
        |q: usize, s: usize| crate::utils::matrix::protein_score(matrix, query[q], subject[s]);
    let mut score: i32 = (0..11).map(|i| score_at(q_start + i, s_start + i)).sum();
    let mut max_score = score;
    let mut max_offset = q_start + 10;
    for i in 11..q_length.min(s_length) {
        score -= score_at(q_start + i - 11, s_start + i - 11);
        score += score_at(q_start + i, s_start + i);
        if score > max_score {
            max_score = score;
            max_offset = q_start + i;
        }
    }
    if max_score > 0 {
        return Some((max_offset as i32, (max_offset - q_start + s_start) as i32));
    }
    let end_score: i32 = (0..11)
        .map(|i| score_at(q_start + q_length - 11 + i, s_start + s_length - 11 + i))
        .sum();
    if end_score > 0 {
        Some((h.q_end - 5, h.s_end - 5))
    } else {
        None
    }
}

#[cfg(test)]
#[path = "traceback_tests.rs"]
mod tests;
