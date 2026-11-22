//! The steps of NCBI's ungapped blastn (`-ungapped`) that run on a subject's HSP list after its
//! word finding and ungapped extension: the re-evaluation with ambiguity codes, the identity count
//! and the word cutoff of the batch (`Blast_HSPTest` is `core::blast_hits::hsp_test`). Pure
//! functions; the engine (`blast_engine/run.rs`) calls them for its `UNGAPPED` instance of
//! `search_query_batch`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1482-1541

use super::ncbi_cutoffs::compute_blastn_cutoff_score_ungapped;
use crate::stats::KarlinParams;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1060-1068 (BLASTNA_SIZE)
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

/// `Blast_HSPReevaluateWithAmbiguitiesUngapped` for a nucleotide HSP (`translated` false, so the
/// query letters are masked with 0x0f), with `s_UpdateReevaluatedHSPUngapped`. `query` starts at the
/// first letter of the HSP's context (NCBI's `query_start`); `subject` is the subject in BLASTNA.
/// Returns `None` when the HSP is to be deleted (its best score is below `cutoff_score`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:676-735
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:440-484 (`s_UpdateReevaluatedHSP`
/// with `gapped` false)
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-790
/// (`s_Blast_HSPGetNumIdentitiesAndPositives`, the branch without `hsp->gap_info`)
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

/// A limit of the word cutoff's arithmetic that LOSAT refuses instead of reproducing NCBI's
/// overflow (memory "no NCBI UB / overflow"): `2 * query_length` must fit NCBI's `Int4`, and
/// the product of `subject` and the smaller of `subject` and twice `query` must fit the signed 64
/// bits that NCBI's `Uint8` is cast to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CutoffLimit {
    /// `2 * context length` does not fit an `Int4`.
    QueryLength,
    /// The search space, `subject` times the smaller of `subject` and twice `query`, does not fit
    /// an `Int8`.
    SearchSpace,
}

/// Whether `2 * context length` fits an `Int4`, as NCBI's `query_length * 2` in the word cutoff
/// needs; a longer context (1 Gb or more) is refused ([`CutoffLimit::QueryLength`]), before any
/// output (`x_rpt_stream::check` asks it of every query of the search).
pub fn cutoff_context_length_supported(context_length: usize) -> bool {
    context_length
        .checked_mul(2)
        .is_some_and(|n| i32::try_from(n).is_ok())
}

/// The ungapped word cutoff of one valid context: `new_cutoff` of
/// `BlastInitialWordParametersUpdate` for an ungapped blastn search.
///
/// `avg_subject_length` is the subject length that NCBI passes (the average length of the
/// subjects, see [`average_subject_length_subject_set`] and
/// [`SeqDb::average_length`](crate::seqdb::SeqDb::average_length)),
/// `kbp_std` the ungapped Karlin block of the context, and `cutoff_score_max` the context's
/// `hit_params->cutoffs[context].cutoff_score_max`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:346-374
pub fn word_cutoff_score(
    context_length: usize,
    avg_subject_length: u64,
    kbp_std: &KarlinParams,
    cutoff_score_max: i32,
) -> Result<i32, CutoffLimit> {
    if !cutoff_context_length_supported(context_length) {
        return Err(CutoffLimit::QueryLength);
    }
    let query_length = context_length as u64 * 2;
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

/// The average subject length of a `-subject` set: the total over the count, as a `Uint4`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:245-266
pub fn average_subject_length_subject_set(total_length: u64, num_seqs: u64) -> u64 {
    // `Int8 / Uint4`, then the cast to `Uint4`.
    (total_length / num_seqs.max(1)) as u32 as u64
}

#[cfg(test)]
#[path = "ungapped_tests.rs"]
mod tests;
