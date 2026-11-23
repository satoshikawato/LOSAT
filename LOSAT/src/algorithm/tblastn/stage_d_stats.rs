//! TBLASTN local-subject statistical inputs before preliminary search.

use crate::algorithm::tblastx::ncbi_cutoffs::{
    cutoff_score_from_evalue, cutoff_score_sum_stats, gap_trigger_raw_score, x_drop_raw_score,
};
use crate::stats::setup::SearchSpaceOverflow;
use crate::stats::spouge::{blast_spouge_etos, BlastGumbelBlk};
use crate::stats::tables::KarlinParams;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_query_info.h:60-76:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1446-1467:
// Uses the subject sequence length for statistics, dividing by 3 for translated subjects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LocalContextLength {
    pub length_adjustment: i64,
    pub eff_searchsp: i64,
    pub subject_stat_length: i64,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:69-71;
// c++/src/algo/blast/core/blast_parameters.c:604-625,793-805:
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct LocalLinkParameters {
    pub gap_decay_rate: f64,
    pub gap_size: i32,
    pub overlap_size: i32,
    pub longest_intron: i32,
    pub cutoff_small_gap: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:902-999,302-383:
// Clamps each context's cutoff score to `cutoff_score_max`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LocalContextCutoffs {
    pub hit_cutoff: i32,
    pub hit_cutoff_max: i32,
    pub word_cutoff: i32,
    pub word_xdrop: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1011-1024:
#[derive(Clone, Debug, PartialEq)]
pub(super) struct LocalSubjectParameters {
    pub lengths: Vec<LocalContextLength>,
    pub cutoffs: Vec<LocalContextCutoffs>,
    pub link: Option<LocalLinkParameters>,
    pub hit_cutoff_min: i32,
    pub word_cutoff_min: i32,
    pub prelim_evalue: f64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:774-815,831-856;
// c++/src/algo/blast/core/blast_setup.c:1016-1024:
// The latter call passes zero compositionBasedStats even when the search uses
// composition mode 2; the score block's Gumbel state is still consulted.
#[derive(Clone, Copy)]
pub(super) struct LocalParameterOptions<'a> {
    pub expect_value: f64,
    pub do_sum_stats: bool,
    pub max_intron_length: i32,
    pub gap_trigger_bits: f64,
    pub word_xdrop_bits: f64,
    pub scale_factor: f64,
    pub gumbel: Option<&'a BlastGumbelBlk>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:676-697: the options' effective search
    // spaces (`searchsp_eff`), which NCBI sets from the whole batch for a query chunk.
    pub eff_searchsp_override: Option<&'a [i64]>,
    /// `-ungapped`: the search has no gapped Karlin-Altschul blocks (`LocalUngappedSearch`);
    /// `None` for a gapped search.
    pub ungapped: Option<LocalUngappedSearch>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:517-527
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1372-1374
// NCBI computes the gapped Karlin blocks only for a gapped search and drops the Gumbel block of
// an ungapped one; the word parameters take the sequence source's average subject length.
/// What an ungapped search (`-ungapped`) changes in `local_parameters_for_blocks`: its score
/// block has no gapped blocks and no Gumbel block, so the cutoffs and the linking read the
/// ungapped `kbp` (the caller's `LocalKarlinBlocks::gapped`), the sum statistics do not lower the
/// hit cutoffs, the linking decays at `BLAST_GAP_DECAY_RATE` and keeps NCBI's ungapped
/// intron length, and the word cutoff is an e-value cutoff below the gap trigger. That cutoff
/// reads the average subject length of the sequence source (`BlastSeqSrcGetAvgSeqLen`, in
/// nucleotides, `average_subject_length`) for the query set, a subject's own length where NCBI
/// computes the statistics of each subject (`LocalParameterCall::Update`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LocalUngappedSearch {
    pub average_subject_length: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:964-985,1001-1024:
// The update is NCBI's statistics of each subject (`db_length == 0`, blast_engine.c:1433-1443):
// the subject's length in nucleotides and a composition-based statistics of 0 (no stretch of the
// e-value).
#[derive(Clone, Copy)]
pub(super) enum LocalParameterCall {
    Initial {
        min_subject_length: i32,
        composition_based_stats: i32,
    },
    Update {
        subject_length: i32,
    },
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-735,770-847:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1442,1446-1467:
#[allow(dead_code)] // Stage D's public path remains gated until C/D/E are complete.
pub(super) fn local_subject_effective_lengths(
    query_contexts: &[(usize, bool)],
    subject_nt_length: i64,
    db_num_seqs: i32,
    gapped_params: &[KarlinParams],
) -> Result<Vec<LocalContextLength>, SearchSpaceOverflow> {
    local_subject_effective_lengths_with_search_spaces(
        query_contexts,
        subject_nt_length,
        db_num_seqs,
        gapped_params,
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:676-697
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:778-847
// A query chunk reads the search spaces of the whole batch by its own context index (RP-4);
// its length adjustment is that of the chunk's part. The warning is not written: a chunk
// with one search space has one context.
// The lengths are `stats::setup::calc_eff_lengths` (one owner for the six programs, line RG
// D1a-3): `subject_nt_length` is NCBI's `Int8` `real_db_length` (blast_parameters.h:217; a
// database's total, which an alias file's `LENGTH` declares, `-dbsize` or a subject's length),
// `db_num_seqs` is NCBI's `Int4` (negative for an alias file's negative `NSEQ`), and an effective
// search space past NCBI's `Int8` is refused (`SearchSpaceOverflow`).
pub(super) fn local_subject_effective_lengths_with_search_spaces(
    query_contexts: &[(usize, bool)],
    subject_nt_length: i64,
    db_num_seqs: i32,
    gapped_params: &[KarlinParams],
    search_spaces: Option<&[i64]>,
) -> Result<Vec<LocalContextLength>, SearchSpaceOverflow> {
    assert_eq!(query_contexts.len(), gapped_params.len());
    let search_space_set =
        search_spaces.is_some_and(|spaces| spaces.iter().any(|&space| space != 0));
    let Some((db_length, db_num_seqs)) = crate::stats::setup::eff_lengths_database(
        subject_nt_length,
        db_num_seqs,
        0,
        search_space_set,
        true,
    ) else {
        return Ok(query_contexts
            .iter()
            .map(|_| LocalContextLength {
                length_adjustment: 0,
                eff_searchsp: 0,
                subject_stat_length: 0,
            })
            .collect());
    };
    query_contexts
        .iter()
        .zip(gapped_params)
        .enumerate()
        .map(|(index, (&(query_length, is_valid), params))| {
            let search_space = match search_spaces {
                None | Some([]) => 0,
                Some([only]) => *only,
                Some(spaces) => spaces[index],
            };
            if !is_valid || query_length == 0 {
                return Ok(LocalContextLength {
                    length_adjustment: 0,
                    eff_searchsp: search_space,
                    subject_stat_length: db_length,
                });
            }
            let lengths = crate::stats::setup::calc_eff_lengths(
                query_length as i64,
                db_length,
                db_num_seqs,
                params,
                search_space,
                "TBLASTN",
            )?;
            Ok(LocalContextLength {
                length_adjustment: i64::from(lengths.length_adjustment),
                eff_searchsp: lengths.eff_searchsp,
                subject_stat_length: db_length,
            })
        })
        .collect()
}

#[allow(dead_code)] // The public TBLASTN path remains gated until C/D/E pass.
pub(super) fn local_parameters_for_call(
    query_contexts: &[(usize, bool)],
    subject_nt_length: i64,
    db_num_seqs: i32,
    gapped_params: &[KarlinParams],
    ungapped_params: &[KarlinParams],
    options: LocalParameterOptions<'_>,
    call: LocalParameterCall,
) -> Result<LocalSubjectParameters, SearchSpaceOverflow> {
    local_parameters_for_blocks(
        query_contexts,
        subject_nt_length,
        db_num_seqs,
        LocalKarlinBlocks {
            gapped: gapped_params,
            lengths: gapped_params,
            ungapped: ungapped_params,
            x_dropoff: ungapped_params,
        },
        options,
        call,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2826
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:120-126,766-768
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:204-221,338-345
/// The Karlin-Altschul blocks of each context that `local_parameters_for_blocks` reads: the
/// gapped `kbp_gap` of the cutoffs and the linking, the gapped `kbp_gap_std` of the effective
/// lengths, the ungapped `kbp_std` of the gap trigger and the ungapped `kbp` of the word X-drop.
/// A matrix query's are two blocks (`kbp_gap` is `kbp_gap_std`, `kbp` is `kbp_std`); a PSSM
/// query's `kbp_gap` and `kbp` are its `kbp_gap_psi` and `kbp_psi` (tblastn `-in_pssm`).
#[derive(Clone, Copy)]
pub(super) struct LocalKarlinBlocks<'a> {
    pub gapped: &'a [KarlinParams],
    pub lengths: &'a [KarlinParams],
    pub ungapped: &'a [KarlinParams],
    pub x_dropoff: &'a [KarlinParams],
}

/// `local_parameters_for_call` with each kind of Karlin-Altschul block of its own
/// (`LocalKarlinBlocks`).
pub(super) fn local_parameters_for_blocks(
    query_contexts: &[(usize, bool)],
    subject_nt_length: i64,
    db_num_seqs: i32,
    blocks: LocalKarlinBlocks<'_>,
    options: LocalParameterOptions<'_>,
    call: LocalParameterCall,
) -> Result<LocalSubjectParameters, SearchSpaceOverflow> {
    let gapped_params = blocks.gapped;
    assert!(!query_contexts.is_empty());
    assert_eq!(query_contexts.len(), gapped_params.len());
    assert_eq!(query_contexts.len(), blocks.lengths.len());
    assert_eq!(query_contexts.len(), blocks.ungapped.len());
    assert_eq!(query_contexts.len(), blocks.x_dropoff.len());
    let lengths = local_subject_effective_lengths_with_search_spaces(
        query_contexts,
        subject_nt_length,
        db_num_seqs,
        blocks.lengths,
        options.eff_searchsp_override,
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:964-985,1011-1024:
    // initial creation passes min_subject_length and compositionBasedStats;
    // the conditional update passes subject_length and literal zero (`LocalParameterCall`).
    let (hit_subject_length, composition_based_stats) = match call {
        LocalParameterCall::Initial {
            min_subject_length,
            composition_based_stats,
        } => (min_subject_length, composition_based_stats),
        LocalParameterCall::Update { subject_length } => (subject_length, 0),
    };
    let gapped = options.ungapped.is_none();
    let mut link = if options.do_sum_stats {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:774-815: the conversion of
        // `-max_intron_length` has one owner, shared with BLASTX
        // (`translated_longest_intron`, with NCBI's lines): a gapped search's value is positive
        // or the linking is off; an ungapped search's is at least 0, and 0 links with even gaps
        // (`stage_d_linking::link_even_gap_translated`).
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:609-615
        // An ungapped search links with the ungapped gap probability and decay rate (0.5), a
        // gapped one with the gapped values (0.1).
        crate::core::link_hsps::translated_longest_intron(options.max_intron_length, gapped).map(
            |intron| LocalLinkParameters {
                gap_decay_rate: if gapped {
                    crate::core::link_hsps::BLAST_GAP_DECAY_RATE_GAPPED
                } else {
                    crate::core::link_hsps::BLAST_GAP_DECAY_RATE
                },
                gap_size: 40,
                overlap_size: 9,
                longest_intron: intron,
                cutoff_small_gap: 0,
            },
        )
    } else {
        None
    };
    let mut cutoffs = Vec::with_capacity(query_contexts.len());
    let mut hit_cutoff_min = i32::MAX;
    let mut word_cutoff_min = i32::MAX;
    let mut prelim_evalue = options.expect_value;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:902-999:
    // Computes per-context hit-saving cutoffs from `query_length` and `avg_subject_length`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:931-946:
    let mut hit = Vec::with_capacity(query_contexts.len());
    for (context, (&(query_length, valid), params)) in
        query_contexts.iter().zip(gapped_params).enumerate()
    {
        if !valid || query_length == 0 {
            hit.push((i32::MAX, 0));
            continue;
        }
        let max_cutoff = if let Some(gumbel) = options.gumbel {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:931-945:
            let cbs_stretch = if composition_based_stats > 1 {
                5.0
            } else {
                1.0
            };
            prelim_evalue = cbs_stretch * options.expect_value;
            blast_spouge_etos(
                prelim_evalue,
                params,
                gumbel,
                query_length as i32,
                hit_subject_length,
            )
        } else {
            cutoff_score_from_evalue(options.expect_value, lengths[context].eff_searchsp, params)
                .max(1)
        };
        hit.push((max_cutoff, max_cutoff));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:952
    // NCBI lowers the hit cutoffs for sum statistics in a gapped search only: an ungapped search
    // keeps the cutoffs of its e-value.
    if let Some(link_params) = link.filter(|_| gapped) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:950-976:
        // Concatenated query length equals the last context's query offset plus its query length.
        let concat_qlen: usize = query_contexts
            .iter()
            .map(|(length, _)| length + 1)
            .sum::<usize>()
            - 1;
        let avg_qlen = (concat_qlen / query_contexts.len()) as i32;
        for (context, (&(_, valid), params)) in query_contexts.iter().zip(gapped_params).enumerate()
        {
            if valid {
                let cutoff = cutoff_score_sum_stats(
                    avg_qlen,
                    hit_subject_length,
                    link_params.gap_decay_rate,
                    params,
                );
                hit[context].0 = hit[context].0.min(cutoff);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:983-999,302-383:
    // hit cutoff and maximum are scaled before the word cutoff is computed;
    // gap_trigger uses kbp_std and the word cutoff is capped by hit cutoff_max.
    for (context, ((&(query_length, valid), ungapped), x_dropoff_kbp)) in query_contexts
        .iter()
        .zip(blocks.ungapped)
        .zip(blocks.x_dropoff)
        .enumerate()
    {
        if !valid || query_length == 0 {
            cutoffs.push(LocalContextCutoffs {
                hit_cutoff: i32::MAX,
                hit_cutoff_max: 0,
                word_cutoff: i32::MAX,
                word_xdrop: 0,
            });
            continue;
        }
        let hit_cutoff = hit[context].0 * options.scale_factor as i32;
        let hit_cutoff_max = hit[context].1 * options.scale_factor as i32;
        hit_cutoff_min = hit_cutoff_min.min(hit_cutoff);
        let gap_trigger = gap_trigger_raw_score(options.gap_trigger_bits, ungapped);
        let new_cutoff = match options.ungapped {
            None => gap_trigger,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:348-367
            // An ungapped search's word cutoff is the score of e-value 1 (`CUTOFF_E_TBLASTN`) in
            // the space of the shorter of subject and query times the subject, with the linking's
            // gap decay, at most the
            // gap trigger (`ungapped_word_cutoff`).
            // `kbp_array` is `kbp_std` without gapped blocks (blast_parameters.c:298-302),
            // `gap_decay_rate` the linking's, 0 without sum statistics (blast_parameters.c:288,
            // 310-311).
            Some(search) => {
                let subj_length = match call {
                    LocalParameterCall::Initial { .. } => search.average_subject_length,
                    LocalParameterCall::Update { subject_length } => subject_length,
                };
                ungapped_word_cutoff(
                    query_length,
                    subj_length,
                    link.map_or(0.0, |link| link.gap_decay_rate),
                    ungapped,
                )
                .min(gap_trigger)
            }
        };
        let word_cutoff = ((new_cutoff as f64 * options.scale_factor) as i32).min(hit_cutoff_max);
        let xdrop_init =
            x_drop_raw_score(options.word_xdrop_bits, x_dropoff_kbp, options.scale_factor);
        let word_xdrop = if xdrop_init == 0 {
            word_cutoff
        } else {
            xdrop_init
        };
        word_cutoff_min = word_cutoff_min.min(word_cutoff);
        cutoffs.push(LocalContextCutoffs {
            hit_cutoff,
            hit_cutoff_max,
            word_cutoff,
            word_xdrop,
        });
    }
    if let Some(link_params) = &mut link {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:625-649:
        link_params.cutoff_small_gap = word_cutoff_min;
    }
    Ok(LocalSubjectParameters {
        lengths,
        cutoffs,
        link,
        hit_cutoff_min,
        word_cutoff_min,
        prelim_evalue,
    })
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:79
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4108-4129
/// The `BLAST_Cutoffs` of an ungapped TBLASTN search's word cutoff (`local_parameters_for_blocks`,
/// which cites the call): the score of e-value `CUTOFF_E_TBLASTN` (1.0) with the gap decay
/// `gap_decay_rate` in the search space (the shorter of subject and query times the subject) of
/// NCBI's `Uint8` (the subject length is NCBI's `Uint4` argument, in nucleotides), at least the
/// cutoff of 1 that NCBI passes in. `cutoff_score_from_evalue` is NCBI's `BlastKarlinEtoS_simple`.
fn ungapped_word_cutoff(
    query_length: usize,
    subj_length: i32,
    gap_decay_rate: f64,
    kbp: &KarlinParams,
) -> i32 {
    const CUTOFF_E_TBLASTN: f64 = 1.0;
    let subj_length = u64::from(subj_length as u32);
    // C's conversion of the `Uint8` product to the `Int8` argument (two's complement).
    let searchsp = (subj_length.min(query_length as u64)).wrapping_mul(subj_length) as i64;
    let e = if gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
        CUTOFF_E_TBLASTN * crate::stats::sum_statistics::gap_decay_divisor(gap_decay_rate, 1)
    } else {
        CUTOFF_E_TBLASTN
    };
    cutoff_score_from_evalue(e, searchsp, kbp).max(1)
}

#[cfg(test)]
#[path = "stage_d_stats_tests.rs"]
mod tests;
