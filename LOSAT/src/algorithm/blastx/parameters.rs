//! BLASTX score-block and preliminary parameter setup, using Rust calculations.
use super::{args::ResolvedOptions, query_setup::PreparedQueryBatch};
use crate::algorithm::tblastx::ncbi_cutoffs::{
    cutoff_score_from_evalue, cutoff_score_sum_stats, gap_trigger_raw_score, x_drop_raw_score,
};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::karlin_calc::{
    apply_check_ideal, compute_aa_composition, compute_karlin_params_ungapped,
    compute_score_freq_profile_for_matrix, compute_std_aa_composition,
};
use crate::stats::length_adjustment::compute_length_adjustment_ncbi;
use crate::stats::spouge::{blast_spouge_etos, lookup_protein_gumbel_params};
use crate::stats::{lookup_protein_params, KarlinParams};
use anyhow::{ensure, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2743-2754
// ```c++
//    Int4 context; /* loop variable. */
//    Blast_ResFreq* rfp,* stdrfp;
//    BlastContextInfo* contexts = query_info->contexts;
//    Boolean check_ideal =
//       (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
//        program == eBlastTypeRpsTblastn);
//    Boolean valid_context = FALSE;
//
//    ASSERT(contexts);
//
//    /* Ideal Karlin block is filled unconditionally. */
//    status = Blast_ScoreBlkKbpIdealCalc(sbp);
// ```
#[derive(Clone, Debug)]
pub struct ContextParameters {
    pub valid: bool,
    pub ungapped: KarlinParams,
    pub search_space: i64,
    pub length_adjustment: i32,
    pub word_cutoff: i32,
    pub word_xdrop_init: i32,
    pub word_xdrop: i32,
    pub hit_cutoff: i32,
    pub hit_cutoff_max: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2776-2801
// ```c++
//       buffer = &query[context_offset];
//
//       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
//       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
//       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
//       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
//       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
//       if (loop_status) {
//           contexts[context].is_valid = FALSE;
//           sbp->sfp[context] = Blast_ScoreFreqFree(sbp->sfp[context]);
//           sbp->kbp_std[context] = Blast_KarlinBlkFree(sbp->kbp_std[context]);
//           if (!Blast_QueryIsTranslated(program) ) {
//              Blast_MessageWrite(blast_message, eBlastSevWarning, context,
//              kBlastErrMsg_CantCalculateUngappedKAParams);
//           }
//           continue;
//       }
//       /* For searches with translated queries, check whether ideal values
//          should be substituted instead of calculated values, so a more
//          conservative (smaller) Lambda is used. */
//       if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
//          Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
//
//       sbp->kbp_psi[context] = Blast_KarlinBlkNew();
//       loop_status = Blast_KarlinBlkUngappedCalc(sbp->kbp_psi[context],
//                                            sbp->sfp[context]);
// ```
pub fn score_block(batch: &PreparedQueryBatch) -> Vec<ContextParameters> {
    let std = compute_std_aa_composition();
    let ideal_profile =
        compute_score_freq_profile_for_matrix(&std, &std, -4, 11, ScoringMatrix::Blosum62);
    let ideal = compute_karlin_params_ungapped(&ideal_profile).expect("BLOSUM62 ideal KA block");
    batch
        .contexts
        .iter()
        .map(|c| {
            let mut seq = vec![0];
            seq.extend_from_slice(&batch.sequence_start[1 + c.offset..1 + c.offset + c.length]);
            seq.push(0);
            let comp = compute_aa_composition(&seq, c.length);
            let profile =
                compute_score_freq_profile_for_matrix(&comp, &std, -4, 11, ScoringMatrix::Blosum62);
            let ka = compute_karlin_params_ungapped(&profile)
                .ok()
                .map(|p| apply_check_ideal(p, ideal));
            ContextParameters {
                valid: c.is_valid && ka.is_some(),
                ungapped: ka.unwrap_or_default(),
                search_space: 0,
                length_adjustment: 0,
                word_cutoff: i32::MAX,
                word_xdrop_init: 0,
                word_xdrop: 0,
                hit_cutoff: i32::MAX,
                hit_cutoff_max: 0,
            }
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:784-824
// ```c++
//       if (query_info->contexts[index].is_valid &&
//           ((query_length = query_info->contexts[index].query_length) > 0) ) {
//
//          /* Use the correct Karlin block. For blastn, two identical Karlin
//           * blocks are allocated for each sequence (one per strand), but we
//           * only need one of them.
//           */
//          if (program_number == eBlastTypeBlastn) {
//              /* Setting reward and penalty to zero is being used to indicate
//               * that matrix scoring should be used for ungapped and gapped
//               * alignment.  For now reward/penalty are being reset to the
//               * default blastn values to not disturb the KA calcs  -RMH- */
//              if ( scoring_options->reward == 0 && scoring_options->penalty == 0 )
//              {
//                  Blast_GetNuclAlphaBeta(BLAST_REWARD,
//                                     BLAST_PENALTY,
//                                     scoring_options->gap_open,
//                                     scoring_options->gap_extend,
//                                     sbp->kbp_std[index],
//                                     scoring_options->gapped_calculation,
//                                     &alpha, &beta);
//              }else {
//                  Blast_GetNuclAlphaBeta(scoring_options->reward,
//                                     scoring_options->penalty,
//                                     scoring_options->gap_open,
//                                     scoring_options->gap_extend,
//                                     sbp->kbp_std[index],
//                                     scoring_options->gapped_calculation,
//                                     &alpha, &beta);
//              }
//          } else {
//              BLAST_GetAlphaBeta(sbp->name, &alpha, &beta,
//                                 scoring_options->gapped_calculation,
//                                 scoring_options->gap_open,
//                                 scoring_options->gap_extend,
//                                 sbp->kbp_std[index]);
//          }
//          BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK,
//                                        alpha/kbp->Lambda, beta,
//                                        query_length, db_length,
//                                        db_num_seqs, &length_adjustment);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:836-847
// ```c++
//         	 Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
//
//         	 // Just in case effective_db_length < 0
//         	 if (effective_db_length <= 0)
//         		 effective_db_length = 1;
//
//              effective_search_space = effective_db_length *
//                              (query_length - length_adjustment);
//          }
//       }
//       query_info->contexts[index].eff_searchsp = effective_search_space;
//       query_info->contexts[index].length_adjustment = length_adjustment;
// ```
pub fn effective_lengths(
    batch: &PreparedQueryBatch,
    parameters: &mut [ContextParameters],
    options: &ResolvedOptions,
    db_length: i64,
    db_count: i64,
    fixed_spaces: Option<&[i64]>,
) -> Result<()> {
    ensure!(
        parameters.len() == batch.contexts.len(),
        "one statistical block per BLASTX context"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-732
    // ```c
    //    if (db_length == 0 &&
    //        !BlastEffectiveLengthsOptions_IsSearchSpaceSet(eff_len_options)) {
    //       return 0;
    //    }
    // ```
    if db_length == 0 && fixed_spaces.is_none() {
        return Ok(());
    }
    let gapped = lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    for (i, (c, p)) in batch.contexts.iter().zip(parameters).enumerate() {
        if !p.valid {
            p.length_adjustment = 0;
            p.search_space = fixed_spaces.map_or(0, |spaces| spaces[i]);
            continue;
        }
        let ka = if options.gapped { gapped } else { p.ungapped };
        let adj = compute_length_adjustment_ncbi(c.length as i64, db_length, db_count, &ka)
            .length_adjustment;
        p.length_adjustment = i32::try_from(adj)?;
        p.search_space = match fixed_spaces {
            Some(spaces) => spaces[i],
            None => (db_length - db_count * adj)
                .max(1)
                .checked_mul(c.length as i64 - adj)
                .ok_or_else(|| anyhow::anyhow!("BLASTX effective search space exceeds Int8"))?,
        };
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-944
// ```c++
//    if (seq_src) {
//       total_length = BlastSeqSrcGetTotLenStats(seq_src);
//       if (total_length <= 0)
//           total_length = BlastSeqSrcGetTotLen(seq_src);
//
//       /* Set the database length for new FSC */
//       if (sbp->gbp) {
//           Int8 dbl = total_length;
//           /* Override the database length or set one (e.g., blast2seq). */
//           if (eff_len_options->db_length) {
//               dbl = eff_len_options->db_length;
//           }
//           sbp->gbp->db_length =
//               (Blast_SubjectIsTranslated(program_number))?
//               dbl/3 : dbl;
//       }
//
//       if (total_length > 0) {
//           num_seqs = BlastSeqSrcGetNumSeqsStats(seq_src);
//           if (num_seqs <= 0)
//               num_seqs = BlastSeqSrcGetNumSeqs(seq_src);
//       } else {
//           /* Not a database search; each subject sequence is considered
//              individually */
//           Int4 oid = 0;  /* Get length of first sequence. */
//           if ( (total_length = BlastSeqSrcGetSeqLen(seq_src, (void*) &oid)) < 0) {
//               total_length = -1;
//               num_seqs = -1;
//           }
//           num_seqs = 1;
//       }
//    }
//
//    /* Initialize the effective length parameters with real values of
//       database length and number of sequences */
//    BlastEffectiveLengthsParametersNew(eff_len_options, total_length, num_seqs,
//                                       eff_len_params);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-987
// ```c++
//    if (sbp->gbp) {
//        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
//        if (Blast_SubjectIsTranslated(program_number)) {
//            min_subject_length/=3;
//        }
//    } else {
//        min_subject_length = (Int4) (total_length/num_seqs);
//    }
//
//    if(min_subject_length <=0) {
// 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
//    }
//
//    if ((status = BlastHitSavingParametersNew(program_number, hit_options, sbp,
// 		                                     query_info, min_subject_length,
// 		                                     (*ext_params)->options->compositionBasedStats,
// 		                                     hit_params)) != 0){
//        return status;
//    }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1372-1374
// ```c++
//     BlastInitialWordParametersNew(program_number, word_options,
//       hit_params, lookup_wrap, sbp, query_info,
//       BlastSeqSrcGetAvgSeqLen(seq_src), &word_params);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:935-975
// ```c++
//          if (sbp->gbp && sbp->gbp->filled) {
// 	     /* If cbs greater than 1 (2 or 3), then increase expect value by 5 for preliminary search. */
// 	     int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
// 	     params->prelim_evalue = cbs_stretch*evalue;
//              new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
//                          query_info->contexts[context].query_length,
//                          avg_subject_length);
//          } else {
//              BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
//          }
//          params->cutoffs[context].cutoff_score = new_cutoff;
//          params->cutoffs[context].cutoff_score_max = new_cutoff;
//       }
//
//       /* If using sum statistics, use a modified cutoff score
//          if that turns out smaller */
//       if (params->link_hsp_params && gapped_calculation) {
//
//          double evalue_hsp = 1.0;
//          Int4 concat_qlen =
//              query_info->contexts[query_info->last_context].query_offset +
//              query_info->contexts[query_info->last_context].query_length;
//          Int4 avg_qlen = concat_qlen / (query_info->last_context + 1);
//          Int8 searchsp = (Int8)MIN(avg_qlen, avg_subject_length) *
//                          (Int8)avg_subject_length;
//
//          ASSERT(params->link_hsp_params);
//
//          for (context = query_info->first_context;
//                              context <= query_info->last_context; ++context) {
//             Int4 new_cutoff = 1;
//
//             if (!(query_info->contexts[context].is_valid))
//                 continue;
//
//             kbp = kbp_array[context];
//             ASSERT(s_BlastKarlinBlkIsValid(kbp));
//             BLAST_Cutoffs(&new_cutoff, &evalue_hsp, kbp, searchsp,
//                        TRUE, params->link_hsp_params->gap_decay_rate);
//             params->cutoffs[context].cutoff_score = MIN(new_cutoff,
//                                     params->cutoffs[context].cutoff_score);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:340-383
// ```c++
//       if (sbp->kbp_std) {     /* this may not be set for gapped blastn */
//          kbp = sbp->kbp_std[context];
//          if (s_BlastKarlinBlkIsValid(kbp)) {
//             gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 +
//                                      kbp->logK) / kbp->Lambda);
//          }
//       }
//
//       if (!gapped_calculation || sbp->matrix_only_scoring) {
//          double cutoff_e = s_GetCutoffEvalue(program_number);
//          Int4 query_length = query_info->contexts[context].query_length;
//
//          /* include the length of reverse complement for blastn searchs. */
//          ASSERT(query_length > 0);
//          if (program_number == eBlastTypeBlastn ||
//              program_number == eBlastTypeMapping)
//             query_length *= 2;
//
//          kbp = kbp_array[context];
//          ASSERT(s_BlastKarlinBlkIsValid(kbp));
//          BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
//                        MIN((Uint8)subj_length,
//                            (Uint8)query_length)*((Uint8)subj_length),
//                        TRUE, gap_decay_rate);
//
//          /* Perform this check for compatibility with the old code */
//          if (program_number != eBlastTypeBlastn)
//             new_cutoff = MIN(new_cutoff, gap_trigger);
//       } else {
//          new_cutoff = gap_trigger;
//       }
//       new_cutoff *= (Int4)sbp->scale_factor;
//       new_cutoff = MIN(new_cutoff,
//                        hit_params->cutoffs[context].cutoff_score_max);
//       curr_cutoffs->cutoff_score = new_cutoff;
//
//       /* Note that x_dropoff_init stays constant throughout the search,
//          but the cutoff_score and x_dropoff parameters may be updated
//          multiple times, if every subject sequence is treated individually */
//
//       if (curr_cutoffs->x_dropoff_init == 0)
//          curr_cutoffs->x_dropoff = new_cutoff;
//       else
//          curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
// ```
pub fn subject_parameters(
    batch: &PreparedQueryBatch,
    parameters: &mut [ContextParameters],
    options: &ResolvedOptions,
    db_length: i64,
    db_count: i64,
    min_subject_length: usize,
    fixed_spaces: Option<&[i64]>,
) -> Result<i32> {
    let avg_subject_length = i32::try_from(db_length / db_count)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:233-239
    // ```c++
    //     for (index=0; index<(*seq_info)->GetNumSeqs(); ++index)
    //         retval = MIN(retval, (*seq_info)->GetSeqBlk(index)->length);
    //
    //     if(retval < BLAST_SEQSRC_MINLENGTH)
    // 	retval = BLAST_SEQSRC_MINLENGTH;
    //
    //     return retval;
    // ```
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205-205
    // ```c++
    // #define BLAST_SEQSRC_MINLENGTH  10    /**< Default minimal sequence length */
    // ```
    let n = if options.gapped {
        i32::try_from(min_subject_length.max(10))?
    } else {
        avg_subject_length
    };
    effective_lengths(
        batch,
        parameters,
        options,
        db_length,
        db_count,
        fixed_spaces,
    )?;
    let spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    };
    let gapped = lookup_protein_params(&spec);
    let gumbel = if options.gapped {
        lookup_protein_gumbel_params(&spec, db_length)
    } else {
        None
    };
    let decay = if options.sum_stats
        && (!options.gapped
            || options.max_intron_length == 0
            || (options.max_intron_length - 2) / 3 > 0)
    {
        if options.gapped {
            0.1
        } else {
            0.5
        }
    } else {
        0.0
    };
    let last = batch.contexts.last().expect("nonempty contexts");
    let avg_q = i32::try_from((last.offset + last.length) / batch.contexts.len())?;
    for (c, p) in batch.contexts.iter().zip(parameters) {
        if !p.valid {
            continue;
        }
        let ka = if options.gapped { gapped } else { p.ungapped };
        p.hit_cutoff_max = if let Some(gbp) = gumbel.as_ref() {
            blast_spouge_etos(
                if options.composition > 1 {
                    5.0 * options.evalue
                } else {
                    options.evalue
                },
                &ka,
                gbp,
                c.length as i32,
                n,
            )
        } else {
            cutoff_score_from_evalue(options.evalue, p.search_space, &ka).max(1)
        };
        p.hit_cutoff = p.hit_cutoff_max;
        if options.gapped && decay != 0.0 {
            p.hit_cutoff = p
                .hit_cutoff
                .min(cutoff_score_sum_stats(avg_q, n, decay, &ka).max(1));
        }
        let trigger = gap_trigger_raw_score(options.gap_trigger, &p.ungapped);
        p.word_cutoff = if options.gapped {
            trigger
        } else {
            // NCBI blast_parameters.c:143-144 uses CUTOFF_E_BLASTX;
            // c++/include/algo/blast/core/blast_parameters.h:78:
            // #define CUTOFF_E_BLASTX 1.0
            cutoff_score_sum_stats(c.length as i32, avg_subject_length, decay, &ka)
                .max(1)
                .min(trigger)
        }
        .min(p.hit_cutoff_max);
        p.word_xdrop_init = x_drop_raw_score(options.x_dropoff, &p.ungapped, 1.0);
        p.word_xdrop = if p.word_xdrop_init == 0 {
            p.word_cutoff
        } else {
            p.word_xdrop_init
        };
    }
    Ok(if options.gapped {
        (options.gap_x_dropoff * std::f64::consts::LN_2 / gapped.lambda) as i32
    } else {
        0
    })
}
