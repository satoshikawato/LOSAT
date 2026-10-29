//! The steps of NCBI's ungapped blastn (`-ungapped`) that run on a subject's HSP list after its
//! word finding and ungapped extension: the re-evaluation with ambiguity codes, the identity count
//! and `Blast_HSPTest`, and the word cutoff of the batch. Pure functions; the engine
//! (`blast_engine/run.rs`) calls them for its `UNGAPPED` instance of `search_query_batch`.
//!
//! NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1482-1541
//! ```c
//!     if (!score_params->options->gapped_calculation) {
//!         status = Blast_HSPListReevaluateUngapped(program_number, hsp_list, query,
//!                     subject, word_params, hit_params, query_info, sbp, score_params,
//!                     seq_src, gen_code_string);
//!         ...
//!         Blast_HSPListGetEvalues(program_number, query_info, hsp_list, gapped_calculation,
//!                     sbp, 0, 1.0);
//!         s_Blast_HSPListReapByPrelimEvalue(...)
//!         Blast_HSPListReapByQueryCoverage(...)
//!         Blast_HSPListGetBitScores(hsp_list, gapped_calculation, sbp);
//!     }
//! ```

use super::filtering::hsp_test;
use super::ncbi_cutoffs::compute_blastn_cutoff_score_ungapped;
use crate::stats::KarlinParams;

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1060-1068 (BLASTNA_SIZE)
const BLASTNA_SIZE: usize = 16;

/// The range of an ungapped HSP after its re-evaluation: the best-scoring stretch of the
/// HSP's diagonal. Offsets are those of the sequences passed in (context-relative for the query).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReevaluatedUngapped {
    pub score: i32,
    pub q_offset: usize,
    pub q_end: usize,
    pub s_offset: usize,
    pub s_end: usize,
}

/// `Blast_HSPReevaluateWithAmbiguitiesUngapped` for a nucleotide HSP (`translated = FALSE`, so the
/// query letters are masked with 0x0f), with `s_UpdateReevaluatedHSPUngapped`. `query` starts at the
/// first letter of the HSP's context (NCBI's `query_start`); `subject` is the subject in BLASTNA.
/// Returns `None` when the HSP is to be deleted (its best score is below `cutoff_score`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:676-735
/// ```c
///    Int4 hsp_length = hsp->query.end - hsp->query.offset;
///    Int4 cutoff_score = word_params->cutoffs[hsp->context].cutoff_score;
///    ...
///    for (index = 0; index < hsp_length; ++index) {
///       sum += matrix[*query & kResidueMask][*subject];
///       query++;
///       subject++;
///       if (sum < 0) {
///           sum = 0;
///           current_q_start = query;
///           current_s_start = subject;
///          if (score < cutoff_score) {
///             best_q_start = best_q_end = query;
///             best_s_start = best_s_end = subject;
///             score = 0;
///          }
///       } else if (sum > score) {
///          score = sum;
///          best_q_end = query;
///          best_s_end = subject;
///          best_q_start = current_q_start;
///          best_s_start = current_s_start;
///       }
///    }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:440-484 (`s_UpdateReevaluatedHSP`
/// with `gapped = FALSE`)
/// ```c
///     hsp->score = score;
///     if (hsp->score >= cutoff_score) {
///         hsp->query.offset = (Int4)(best_q_start - query_start);
///         hsp->query.end = (Int4)(hsp->query.offset + best_q_end - best_q_start);
///         hsp->subject.offset = (Int4)(best_s_start - subject_start);
///         hsp->subject.end = (Int4)(hsp->subject.offset + best_s_end - best_s_start);
///         delete_hsp = FALSE;
///     }
/// ```
pub fn reevaluate_with_ambiguities_ungapped(
    query: &[u8],
    subject: &[u8],
    q_offset: usize,
    s_offset: usize,
    hsp_length: usize,
    cutoff_score: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> Option<ReevaluatedUngapped> {
    let mut score: i32 = 0;
    let mut sum: i32 = 0;
    let (mut best_q_start, mut best_q_end, mut current_q_start) = (q_offset, q_offset, q_offset);
    let (mut best_s_start, mut best_s_end, mut current_s_start) = (s_offset, s_offset, s_offset);
    let mut q = q_offset;
    let mut s = s_offset;
    for _ in 0..hsp_length {
        sum += matrix[(usize::from(query[q]) & 0x0f) * BLASTNA_SIZE + usize::from(subject[s])];
        q += 1;
        s += 1;
        if sum < 0 {
            // Start from new offset.
            sum = 0;
            current_q_start = q;
            current_s_start = s;
            // If previous top score never reached the cutoff, discard the front part of the
            // alignment completely.
            if score < cutoff_score {
                best_q_start = q;
                best_q_end = q;
                best_s_start = s;
                best_s_end = s;
                score = 0;
            }
        } else if sum > score {
            score = sum;
            best_q_end = q;
            best_s_end = s;
            best_q_start = current_q_start;
            best_s_start = current_s_start;
        }
    }
    if score >= cutoff_score {
        Some(ReevaluatedUngapped {
            score,
            q_offset: best_q_start,
            q_end: best_q_end,
            s_offset: best_s_start,
            s_end: best_s_end,
        })
    } else {
        None
    }
}

/// The number of identities of an ungapped HSP: `*q == *s` in BLASTNA (so `N` and `N` are a
/// match), over `query.end - query.offset` letters. `query` is the unmasked query of the context.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-790
/// (`s_Blast_HSPGetNumIdentitiesAndPositives`, the `!hsp->gap_info` branch)
/// ```c
///       align_length = q_length;
///       for (i=0; i<align_length; i++) {
///          if (*q == *s)
///             num_ident++;
///          ...
///          q++;
///          s++;
///       }
/// ```
pub fn count_identities_ungapped(
    query: &[u8],
    subject: &[u8],
    q_offset: usize,
    s_offset: usize,
    align_length: usize,
) -> usize {
    query[q_offset..q_offset + align_length]
        .iter()
        .zip(&subject[s_offset..s_offset + align_length])
        .filter(|(q, s)| q == s)
        .count()
}

/// `Blast_HSPTest`: true when the HSP is to be deleted (`-perc_identity`, minimum hit length).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1000, 1027-1032
/// ```c
/// return ((hsp->num_ident * 100.0 < align_length * hit_options->percent_identity) ||
///         align_length < hit_options->min_hit_length) ;
/// ```
pub fn blast_hsp_test(
    num_ident: usize,
    align_length: usize,
    percent_identity: f64,
    min_hit_length: usize,
) -> bool {
    hsp_test(num_ident, align_length, percent_identity, min_hit_length)
}

/// A limit of the word cutoff's arithmetic that LOSAT refuses instead of reproducing NCBI's
/// overflow (memory "no NCBI UB / overflow"): `2 * query_length` must fit NCBI's `Int4`, and
/// `MIN(subject, 2 * query) * subject` must fit the signed 64 bits that NCBI's `Uint8` is cast to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CutoffLimit {
    /// `2 * context length` does not fit an `Int4`.
    QueryLength,
    /// The search space `MIN(subject, 2 * query) * subject` does not fit an `Int8`.
    SearchSpace,
}

/// The ungapped word cutoff of one valid context: `new_cutoff` of
/// `BlastInitialWordParametersUpdate` for an ungapped blastn search.
///
/// `avg_subject_length` is the subject length that NCBI passes (the average length of the
/// subjects, see [`average_subject_length_subject_set`] and [`average_subject_length_database`]),
/// `kbp_std` the ungapped Karlin block of the context, and `cutoff_score_max` the context's
/// `hit_params->cutoffs[context].cutoff_score_max`.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:346-374
/// ```c
///          double cutoff_e = s_GetCutoffEvalue(program_number);
///          Int4 query_length = query_info->contexts[context].query_length;
///          if (program_number == eBlastTypeBlastn ||
///              program_number == eBlastTypeMapping)
///             query_length *= 2;
///          kbp = kbp_array[context];
///          BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
///                        MIN((Uint8)subj_length,
///                            (Uint8)query_length)*((Uint8)subj_length),
///                        TRUE, gap_decay_rate);
///       ...
///       new_cutoff *= (Int4)sbp->scale_factor;
///       new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
/// ```
pub fn word_cutoff_score(
    context_length: usize,
    avg_subject_length: u64,
    kbp_std: &KarlinParams,
    cutoff_score_max: i32,
) -> Result<i32, CutoffLimit> {
    let query_length = context_length
        .checked_mul(2)
        .filter(|n| i32::try_from(*n).is_ok())
        .ok_or(CutoffLimit::QueryLength)? as u64;
    let search_space = query_length
        .min(avg_subject_length)
        .checked_mul(avg_subject_length)
        .filter(|n| i64::try_from(*n).is_ok())
        .ok_or(CutoffLimit::SearchSpace)?;
    debug_assert!(search_space <= i64::MAX as u64);
    // `gap_decay_rate` is 0 (the link parameters are NULL), `sbp->scale_factor` is 1.
    Ok(compute_blastn_cutoff_score_ungapped(
        query_length as i64,
        avg_subject_length as i64,
        i32::MAX,
        kbp_std,
        1.0,
        0.0,
        cutoff_score_max,
    ))
}

/// The average subject length of a `-subject` set: `(Uint4)(total / count)`.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/seqsrc_multiseq.cpp:245-266
/// ```c
///     retval = (Uint4) (total / num_seqs)   (s_MultiSeqGetAvgLength)
/// ```
pub fn average_subject_length_subject_set(total_length: u64, num_seqs: u64) -> u64 {
    // `Int8 / Uint4`, then the cast to `Uint4`.
    (total_length / num_seqs.max(1)) as u32 as u64
}

/// The average subject length of a database: `(Int4)(total / MAX(1, count))`.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/seqsrc_seqdb.cpp:203-209
/// ```c
///     return (Int4) (seqdb->GetTotalLength() / MAX(1, seqdb->GetNumSeqs()));  (s_SeqDbGetAvgLength)
/// ```
pub fn average_subject_length_database(total_length: u64, num_seqs: u64) -> u64 {
    (total_length / num_seqs.max(1)) as i32 as u32 as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastn::alignment::gapped::build_blastna_matrix;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:471-478
    // ```c
    // const Uint1 kQuerySeq[kLength+2] = { 15, 0, 1, 2, 3, 0, 0, 1, 1, 2, 2, 3,
    //                                     3, 3, 3, 0, 1, 2, 2, 1, 1, 0, 0, 15 };
    // const Uint1 kSubjectSeq[kLength] = { 0, 1, 2, 2, 14, 14, 1, 1, 2, 2, 3,
    //                                     3, 3, 3, 0, 1, 8, 2, 1, 1, 14, 0 };
    // ```
    const NCBI_QUERY: [u8; 22] = [
        0, 1, 2, 3, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 0, 1, 2, 2, 1, 1, 0, 0,
    ];
    const NCBI_SUBJECT: [u8; 22] = [
        0, 1, 2, 2, 14, 14, 1, 1, 2, 2, 3, 3, 3, 3, 0, 1, 8, 2, 1, 1, 14, 0,
    ];

    /// `testReevaluateWithAmbiguitiesUngappedNucl` (blasthits_unit_test.cpp:500-600): the HSP
    /// 0..20 on 0..20 (score 20) of a 1/-2 search with a word cutoff of 3 is cut to 6..20 on
    /// 6..20 with score 11, and not deleted.
    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:500-585 testReevaluateWithAmbiguitiesUngappedNucl
    #[test]
    fn ncbi_test_reevaluate_with_ambiguities_ungapped_nucl() {
        let matrix = build_blastna_matrix(1, -2);
        let got =
            reevaluate_with_ambiguities_ungapped(&NCBI_QUERY, &NCBI_SUBJECT, 0, 0, 20, 3, &matrix);
        assert_eq!(
            got,
            Some(ReevaluatedUngapped {
                score: 11,
                q_offset: 6,
                q_end: 20,
                s_offset: 6,
                s_end: 20
            })
        );
    }

    /// The same HSP with the cutoff of the other NCBI constant (`kHitCutoff` 11) is kept, with 12
    /// it is deleted (the best stretch scores 11).
    #[test]
    fn reevaluate_ungapped_cutoff_boundary() {
        let matrix = build_blastna_matrix(1, -2);
        let kept =
            reevaluate_with_ambiguities_ungapped(&NCBI_QUERY, &NCBI_SUBJECT, 0, 0, 20, 11, &matrix);
        assert_eq!(
            kept.map(|r| (r.score, r.q_offset, r.q_end)),
            Some((11, 6, 20))
        );
        assert_eq!(
            reevaluate_with_ambiguities_ungapped(&NCBI_QUERY, &NCBI_SUBJECT, 0, 0, 20, 12, &matrix),
            None
        );
    }

    /// A front part that reached the cutoff before the sum went negative stays the best part
    /// unless the later stretch scores more; below the cutoff it is discarded at the break.
    /// (Expected values worked out outside with the same recurrence; an `N` scores -1 at 1/-2.)
    #[test]
    fn reevaluate_ungapped_front_part_kept_only_above_cutoff() {
        let matrix = build_blastna_matrix(1, -2);
        // Two Ns inside: the sum stays positive, the whole stretch is one HSP of score 6.
        let q = [0u8, 1, 2, 3, 0, 0, 0, 1, 2, 3];
        let s = [0u8, 1, 2, 3, 0, 14, 14, 1, 2, 3];
        let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 10, 4, &matrix).unwrap();
        assert_eq!(
            (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
            (6, 0, 10, 0, 10)
        );
        assert_eq!(
            reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 10, 7, &matrix),
            None
        );
        // Three mismatches take the sum below zero (5 - 6): front 5, back 7.
        let q = [0u8, 1, 2, 3, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 3];
        let s = [0u8, 1, 2, 3, 0, 2, 2, 2, 1, 2, 3, 0, 1, 2, 3];
        let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 15, 4, &matrix).unwrap();
        assert_eq!(
            (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
            (7, 8, 15, 8, 15)
        );
        // Cutoff 6: the front (5) is discarded at the break, the back (7) is kept; cutoff 8: both fail.
        let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 15, 6, &matrix).unwrap();
        assert_eq!((got.score, got.q_offset, got.q_end), (7, 8, 15));
        assert_eq!(
            reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 15, 8, &matrix),
            None
        );
        // Front 5 kept at cutoff 4 when the back stretch (3) scores less; at cutoff 6 nothing is left.
        let q = [0u8, 1, 2, 3, 0, 0, 0, 0, 1, 2, 3];
        let s = [0u8, 1, 2, 3, 0, 2, 2, 2, 1, 2, 3];
        let got = reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 11, 4, &matrix).unwrap();
        assert_eq!(
            (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
            (5, 0, 5, 0, 5)
        );
        assert_eq!(
            reevaluate_with_ambiguities_ungapped(&q, &s, 0, 0, 11, 6, &matrix),
            None
        );
    }

    /// Offsets are those of the passed sequences (the HSP need not start at 0 on either).
    #[test]
    fn reevaluate_ungapped_with_offsets_and_sentinel_mask() {
        let matrix = build_blastna_matrix(2, -3);
        let q = [3u8, 3, 0, 1, 2, 3, 0, 1];
        let s = [9u8, 0, 1, 2, 3, 0, 1, 1, 1];
        // HSP: q 2..8 against s 1..7 (6 matches except the last: q1 vs s1 equals, 6 letters).
        let got = reevaluate_with_ambiguities_ungapped(&q, &s, 2, 1, 6, 5, &matrix).unwrap();
        assert_eq!(
            (got.score, got.q_offset, got.q_end, got.s_offset, got.s_end),
            (12, 2, 8, 1, 7)
        );
        // A query letter of 0x10 or more is masked with 0x0f (the translated mask is 0xff).
        let qm = [3u8, 3, 0x10, 1, 2, 3, 0, 1];
        let got2 = reevaluate_with_ambiguities_ungapped(&qm, &s, 2, 1, 6, 5, &matrix).unwrap();
        assert_eq!(got, got2);
    }

    #[test]
    fn identities_count_blastna_codes_equal() {
        // N (14) against N (14) is an identity; N against A is not.
        let q = [0u8, 14, 14, 3, 2];
        let s = [0u8, 14, 0, 3, 1];
        assert_eq!(count_identities_ungapped(&q, &s, 0, 0, 5), 3);
        assert_eq!(count_identities_ungapped(&q, &s, 3, 3, 2), 1);
        assert_eq!(count_identities_ungapped(&q, &s, 2, 2, 0), 0);
    }

    /// `Blast_HSPTest`: `num_ident * 100.0 < align_length * percent_identity` deletes; equality keeps.
    #[test]
    fn hsp_test_boundary() {
        // 95 identities in 100 letters at -perc_identity 95: 9500.0 < 9500.0 is false, kept.
        assert!(!blast_hsp_test(95, 100, 95.0, 0));
        assert!(blast_hsp_test(94, 100, 95.0, 0));
        assert!(!blast_hsp_test(0, 100, 0.0, 0));
        assert!(blast_hsp_test(90, 90, 0.0, 91));
        assert!(!blast_hsp_test(91, 91, 0.0, 91));
    }

    /// The word cutoff of NCBI's `BLAST_Cutoffs(.., MIN(avg, 2q) * avg, TRUE, 0)`: worked by hand
    /// as `ceil(ln(K * searchsp / 0.05) / Lambda)` (values from the formula, computed outside).
    #[test]
    fn word_cutoff_is_ceil_of_etos_with_min_of_avg_and_two_q() {
        let kbp = KarlinParams {
            lambda: 1.33,
            k: 0.621,
            h: 1.12,
            alpha: 0.0,
            beta: 0.0,
        };
        // q 377 (2q 754), avg 373: searchsp 373 * 373 -> ceil(ln(0.621*139129/0.05)/1.33) = 11.
        assert_eq!(word_cutoff_score(377, 373, &kbp, i32::MAX), Ok(11));
        // A short query: 2q 600 < avg 20000 -> searchsp 600 * 20000 (not 300 * 20000).
        assert_eq!(word_cutoff_score(300, 20000, &kbp, i32::MAX), Ok(15));
        // The cutoff is capped by cutoff_score_max.
        assert_eq!(word_cutoff_score(377, 373, &kbp, 7), Ok(7));
        // A very large subject average: 100000 * 100000 -> 20.
        assert_eq!(word_cutoff_score(100_000, 100_000, &kbp, i32::MAX), Ok(20));
        let kbp2 = KarlinParams {
            lambda: 0.634,
            k: 0.408,
            h: 0.912,
            alpha: 0.0,
            beta: -2.0,
        };
        assert_eq!(word_cutoff_score(377, 373, &kbp2, i32::MAX), Ok(22));
    }

    #[test]
    fn word_cutoff_refuses_overflowing_arithmetic() {
        let kbp = KarlinParams {
            lambda: 1.33,
            k: 0.621,
            h: 1.12,
            alpha: 0.0,
            beta: 0.0,
        };
        assert_eq!(
            word_cutoff_score(1 << 30, 1000, &kbp, i32::MAX),
            Err(CutoffLimit::QueryLength)
        );
        assert!(word_cutoff_score((1 << 30) - 1, 1000, &kbp, i32::MAX).is_ok());
        // min(avg, 2q) * avg > i64::MAX (not reachable with NCBI's 32-bit averages; the function
        // takes the average as a u64).
        assert_eq!(
            word_cutoff_score((1 << 30) - 1, 5_000_000_000, &kbp, i32::MAX),
            Err(CutoffLimit::SearchSpace)
        );
    }

    #[test]
    fn average_subject_lengths_follow_each_source() {
        assert_eq!(average_subject_length_subject_set(1000, 3), 333);
        assert_eq!(average_subject_length_database(1000, 0), 1000);
        assert_eq!(average_subject_length_database(10, 3), 3);
        // Uint4 and Int4 casts of a huge average.
        assert_eq!(
            average_subject_length_subject_set(5_000_000_000, 1),
            705_032_704
        );
        // A database average above `Int4` becomes negative and then a huge `Uint4` (3,000,000,000 here).
        assert_eq!(
            average_subject_length_database(3_000_000_000, 1),
            3_000_000_000
        );
    }
}
