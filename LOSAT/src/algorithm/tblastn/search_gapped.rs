//! Internal TBLASTN protein gapped-score boundary.
//! The public TBLASTN CLI remains unsupported until downstream stages agree.

use super::search_init::InitHsp;
use super::search_seed::resolve_local_subject_ncbi2na;
use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::blastp::gapalign::{
    blast_gapped_alignment_with_traceback, blast_gapped_alignment_with_traceback_with_scratch,
    blastp_get_start_for_gapped_alignment, blastp_score_only_gapped_alignment_with_scratch,
    BlastpGappedAlignmentMode, GapAlignScratch,
};
use crate::algorithm::tblastx::translation::generate_frames;
use crate::config::ScoringMatrix;
use crate::utils::genetic_code::GeneticCode;
use anyhow::{Context, Result};

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
        for hit in initial.iter().filter(|hit| hit.frame == frame.frame) {
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
            results.push(GappedHsp {
                frame: frame.frame,
                score: gapped.score,
                q_start: gapped.query_start,
                q_end: gapped.query_stop,
                q_gapped_start: i32::try_from(q_gapped_start)?,
                s_start: gapped.subject_start,
                s_end: gapped.subject_stop,
                s_gapped_start: i32::try_from(s_gapped_start)?,
            });
        }
    }
    Ok(results)
}

// NCBI c++/src/algo/blast/core/blast_traceback.c:408-440,503-509,
// 1644-1677:
// subject = Blast_HSPGetTargetTranslation(target_t, hsp, &subject_length);
// BLAST_GappedAlignmentWithTraceback(program_number, query, subject, ...);
// if (fence_hit) { /* refetch whole subject */ }
//     Blast_TracebackFromHSPList(..., &fence_hit);
// This internal diagnostic models the first-pass right fence and second-pass
// full ncbi4na-derived translation for the fixed small full-frame HSPs. General
// target windows and the public search lifecycle remain unimplemented.
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
    let frames = generate_frames(subject, &code);
    let query_frame = encode_protein_query_frame_with_seg(query, None);
    let query_sequence = &query_frame.aa_seq[1..query_frame.aa_seq.len() - 1];
    let mut results = Vec::new();
    for hit in gapped {
        let frame = frames
            .iter()
            .find(|frame| frame.frame == hit.frame)
            .context("missing translated subject frame")?;
        let subject_sequence = &frame.aa_seq[1..frame.aa_seq.len() - 1];
        let q_start = usize::try_from(hit.q_gapped_start)?;
        let s_start = usize::try_from(hit.s_gapped_start)?;
        // NCBI c++/src/algo/blast/core/blast_hits.c:1160-1222:
        // target_t->translations[context][length+1] = FENCE_SENTRY;
        // NCBI c++/src/algo/blast/core/blast_traceback.c:529-536,1644-1684:
        // BLAST_GappedAlignmentWithTraceback(..., fence_hit);
        // if (fence_hit) { /* refetch whole subject; retry */ }
        // The fixed small-fixture HSPs span the translated subject. The
        // range translation contains the complete frame plus its right fence.
        let mut partial_sequence = subject_sequence.to_vec();
        partial_sequence.push(201); // FENCE_SENTRY, blast_util.h:364
        let mut fence_hit = false;
        let mut scratch = GapAlignScratch::new();
        let first = blast_gapped_alignment_with_traceback_with_scratch(
            query_sequence,
            &partial_sequence,
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
        let alignment = if fence_hit {
            blast_gapped_alignment_with_traceback(
                query_sequence,
                subject_sequence,
                q_start,
                s_start,
                ScoringMatrix::Blosum62,
                None,
                gap_open,
                gap_extend,
                x_drop_final,
            )
        } else {
            first
        }
        .context("NCBI protein traceback returned no alignment")?;
        results.push((
            GappedHsp {
                frame: hit.frame,
                score: alignment.score,
                q_start: i32::try_from(alignment.query_start)?,
                q_end: i32::try_from(alignment.query_stop)?,
                q_gapped_start: hit.q_gapped_start,
                s_start: i32::try_from(alignment.subject_start)?,
                s_end: i32::try_from(alignment.subject_stop)?,
                s_gapped_start: hit.s_gapped_start,
            },
            fence_hit,
        ));
    }
    Ok(results)
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
}
