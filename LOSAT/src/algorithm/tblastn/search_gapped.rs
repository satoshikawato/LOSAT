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
) -> Result<Vec<GappedHsp>> {
    let resolved = resolve_local_subject_ncbi2na(subject)?;
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    let frames = generate_frames(&resolved, &code);
    let query_frame = encode_protein_query_frame_with_seg(query, None);
    let query_sequence = &query_frame.aa_seq[1..query_frame.aa_seq.len() - 1];
    let mut scratch = GapAlignScratch::new();
    let mut results = Vec::new();

    for frame in frames {
        let subject_sequence = &frame.aa_seq[1..frame.aa_seq.len() - 1];
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:223-235:
        // Blast_InitHitListSortByScore(init_hitlist);
        // NCBI c++/src/algo/blast/core/blast_extend.c:273-313:
        // compare score DESC, subject start ASC, length DESC, query start ASC.
        let mut frame_initial: Vec<_> = initial
            .iter()
            .filter(|hit| hit.frame == frame.frame)
            .collect();
        frame_initial.sort_unstable_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.s_start.cmp(&b.s_start))
                .then(b.length.cmp(&a.length))
                .then(a.q_start.cmp(&b.q_start))
        });
        // NCBI c++/src/algo/blast/core/blast_gapalign.c:3827-3834,3908-3919,
        // 4065-4091:
        // tree = Blast_IntervalTreeInit(0, query->length+1, 0, subject->length+1);
        // if (!BlastIntervalTreeContainsHSP(tree, &tmp_hsp, query_info,
        //                                   hit_options->min_diag_separation)) { ... }
        // BlastIntervalTreeAddHSP(new_hsp, tree, query_info, eQueryAndSubject);
        let mut tree = BlastIntervalTree::new(
            0,
            i32::try_from(query.len())? + 1,
            0,
            i32::try_from(frame.aa_len)? + 1,
        );
        for hit in frame_initial {
            let initial_tree_hsp = TreeHsp {
                query_offset: hit.q_start,
                query_end: hit.q_start + hit.length,
                subject_offset: hit.s_start,
                subject_end: hit.s_start + hit.length,
                score: hit.score,
                query_frame: 0,
                query_length: i32::try_from(query.len())?,
                query_context_offset: 0,
                subject_frame_sign: i32::from(frame.frame.signum()),
            };
            if tree.contains_hsp(&initial_tree_hsp, 0, 0) {
                continue;
            }
            let q_start = usize::try_from(hit.q_start).context("negative initial query start")?;
            let s_start = usize::try_from(hit.s_start).context("negative initial subject start")?;
            let length = usize::try_from(hit.length).context("negative initial HSP length")?;
            let q_gapped_start = blastp_get_start_for_gapped_alignment(
                query_sequence,
                subject_sequence,
                q_start,
                length,
                s_start,
                length,
                ScoringMatrix::Blosum62,
            );
            let s_gapped_start =
                i64::from(hit.s_seed) + i64::try_from(q_gapped_start)? - i64::from(hit.q_seed);
            let s_gapped_start =
                usize::try_from(s_gapped_start).context("negative gapped subject start")?;
            let gapped = blastp_score_only_gapped_alignment_with_scratch(
                query_sequence,
                subject_sequence,
                q_gapped_start,
                s_gapped_start,
                ScoringMatrix::Blosum62,
                gap_open,
                gap_extend,
                x_drop,
                BlastpGappedAlignmentMode::Exact,
                &mut scratch,
            );
            // NCBI c++/src/algo/blast/core/blast_gapalign.c:4057-4091:
            // if (gap_align->score >= cutoff) {
            //     Blast_HSPInit(..., gap_align->score, ..., &new_hsp);
            //     Blast_HSPListSaveHSP(hsp_list, new_hsp);
            //     BlastIntervalTreeAddHSP(new_hsp, tree, query_info, eQueryAndSubject);
            // }
            if gapped.score >= 0 {
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
                        query_length: i32::try_from(query.len())?,
                        query_context_offset: 0,
                        subject_frame_sign: i32::from(frame.frame.signum()),
                    },
                    0,
                    IndexMethod::QueryAndSubject,
                );
                results.push(saved);
            }
        }
    }
    Ok(results)
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
                ScoringMatrix::Blosum62,
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
    use crate::algorithm::tblastn::search_init::find_blosum62_word3_init_hsps;
    use std::fs;

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
            for hsp in gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38).unwrap() {
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
        let actual = gapped_blosum62_word3_hsps(query, subject, 32, &initial, 11, 1, 38).unwrap();
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
                gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38).unwrap();
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
        let gapped = gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38).unwrap();
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
        let gapped = gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 12).unwrap();
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
        let gapped = gapped_blosum62_word3_hsps(query, subject, 1, &initial, 11, 1, 38).unwrap();
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
}
