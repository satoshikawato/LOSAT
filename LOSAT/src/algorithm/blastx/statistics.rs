//! BLASTX HSP statistics; protein-subject lengths remain in amino acids.
use super::{
    args::ResolvedOptions,
    linking::{link_parameters, link_uneven, LinkedHsp, LinkedHspList},
    parameters::ContextParameters,
    preliminary::PreliminaryHsp,
    query_setup::PreparedQueryBatch,
};
use crate::{
    config::{ProteinScoringSpec, ScoringMatrix},
    stats::{
        lookup_protein_params,
        spouge::{blast_spouge_stoe, lookup_protein_gumbel_params},
    },
};
use anyhow::{ensure, Result};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

// No NCBI counterpart: reads the LOSAT_X_BXEVALEMPTY / LOSAT_X_BXEVALEMPTYSHADOW switches once; it
// does not change any value NCBI computes.
/// LOSAT_X_BXEVALEMPTY (a strict-set switch, on unless `0`)
/// : `Some(shadow)` when `evaluate_gapped` returns the empty list of a subject without HSPs at
/// once, `shadow` with `LOSAT_X_BXEVALEMPTYSHADOW` set (the full evaluation is run too and
/// compared).
fn x_eval_empty() -> Option<bool> {
    static MODE: OnceLock<Option<bool>> = OnceLock::new();
    *MODE.get_or_init(|| {
        crate::utils::xswitch::on("LOSAT_X_BXEVALEMPTY")
            .then(|| std::env::var_os("LOSAT_X_BXEVALEMPTYSHADOW").is_some())
    })
}

/// LOSAT_X_BXEVALEMPTYSHADOW: the empty lists compared with the full evaluation.
static X_EVAL_EMPTY_COMPARED: AtomicU64 = AtomicU64::new(0);

/// LOSAT_X_BXEVALEMPTYSHADOW: the full evaluation of no input, which must give `empty`.
#[cold]
#[inline(never)]
fn x_eval_empty_shadow(
    empty: &LinkedHspList,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject_length: i32,
    db_length: i64,
    scale: f64,
) {
    let full = x_evaluate_gapped_full(
        &[],
        batch,
        parameters,
        options,
        subject_length,
        db_length,
        scale,
    )
    .expect("the full evaluation of no input");
    if !(full.hsps.is_empty() && full.best_evalue.to_bits() == empty.best_evalue.to_bits()) {
        eprintln!("{}", x_eval_empty_line());
        panic!(
            "LOSAT_X_BXEVALEMPTYSHADOW: no input evaluates to {} HSPs, best e-value {}",
            full.hsps.len(),
            full.best_evalue
        );
    }
    X_EVAL_EMPTY_COMPARED.fetch_add(1, Ordering::Relaxed);
}

/// LOSAT_X_BXEVALEMPTYSHADOW: the shadow line.
fn x_eval_empty_line() -> String {
    format!(
        "[X_SHADOW] BXEVALEMPTY lists_compared={} identical",
        X_EVAL_EMPTY_COMPARED.load(Ordering::Relaxed)
    )
}

/// LOSAT_X_BXEVALEMPTYSHADOW: the shadow line at exit (`x_seed_scan::print_shadow_stats`), after
/// a comparison.
pub(crate) fn x_print_eval_empty_summary() {
    if X_EVAL_EMPTY_COMPARED.load(Ordering::Relaxed) > 0 {
        eprintln!("{}", x_eval_empty_line());
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1811-1839
// ```c++
// Int2 Blast_HSPListGetEvalues(EBlastProgramType program_number,
//                              const BlastQueryInfo* query_info,
//                              Int4 subject_length,
//                              BlastHSPList* hsp_list,
//                              Boolean gapped_calculation,
//                              Boolean RPS_prelim,
//                              const BlastScoreBlk* sbp, double gap_decay_rate,
//                              double scaling_factor)
// {
//    BlastHSP* hsp;
//    BlastHSP** hsp_array;
//    Blast_KarlinBlk** kbp;
//    Int4 hsp_cnt;
//    Int4 index;
//    Int4 kbp_context;
//    Int4 score;
//    double gap_decay_divisor = 1.;
//    Boolean isRPS = Blast_ProgramIsRpsBlast(program_number);
//
//    if (hsp_list == NULL || hsp_list->hspcnt == 0)
//       return 0;
//
//    kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
//    hsp_cnt = hsp_list->hspcnt;
//    hsp_array = hsp_list->hsp_array;
//
//    if (gap_decay_rate != 0.)
//       gap_decay_divisor = BLAST_GapDecayDivisor(gap_decay_rate, 1);
//
// ```
pub fn evaluate_gapped(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject_length: i32,
    db_length: i64,
    scale: f64,
) -> Result<LinkedHspList> {
    ensure!(
        options.gapped,
        "BLASTX even-gap ungapped linking is not implemented"
    );
    // Speed (LOSAT_X_BXEVALEMPTY, in the strict set): a subject without HSPs gets at once the list
    // that both ways below give for no input (`link_uneven`'s first lines, and the loop below: no
    // HSP, best e-value 0), without their tables; the sink evaluates every subject, and most
    // subjects of a database have no HSP. No NCBI counterpart: it does not change any value NCBI
    // computes.
    if input.is_empty() {
        if let Some(shadow) = x_eval_empty() {
            let empty = LinkedHspList {
                hsps: Vec::new(),
                best_evalue: 0.0,
            };
            if shadow {
                x_eval_empty_shadow(
                    &empty,
                    batch,
                    parameters,
                    options,
                    subject_length,
                    db_length,
                    scale,
                );
            }
            return Ok(empty);
        }
    }
    x_evaluate_gapped_full(
        input,
        batch,
        parameters,
        options,
        subject_length,
        db_length,
        scale,
    )
}

/// `evaluate_gapped` after its check, for any input (LOSAT_X_BXEVALEMPTY's shadow evaluates no
/// input here).
#[inline]
fn x_evaluate_gapped_full(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject_length: i32,
    db_length: i64,
    scale: f64,
) -> Result<LinkedHspList> {
    let spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    };
    let mut ka = lookup_protein_params(&spec);
    ka.lambda /= scale;
    let gumbel = super::parameters::gumbel_block(options, db_length);
    let query_lengths: Vec<_> = batch.contexts.iter().map(|c| c.length as i32).collect();
    if let Some(link) = link_parameters(options) {
        return link_uneven(
            input,
            &query_lengths,
            parameters,
            subject_length,
            &vec![ka; parameters.len()],
            gumbel.as_ref(),
            &link,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1875-1906
    // ```c++
    //               hsp->evalue =
    //                   BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
    //                                query_info->contexts[hsp->context].query_length,
    //                                subject_length);
    //           } else {
    //               /* for RPS blast, query and subject is swapped */
    //               hsp->evalue =
    //                   BLAST_SpougeStoE(score, kbp[kbp_context], sbp->gbp,
    //                                subject_length,
    //                                query_info->contexts[hsp->context].query_length);
    //           }
    //       } else {
    //           /* Get effective search space from the query information block */
    //           hsp->evalue =
    //               BLAST_KarlinStoE_simple(score, kbp[kbp_context],
    //                                query_info->contexts[hsp->context].eff_searchsp);
    //       }
    //
    //       hsp->evalue /= gap_decay_divisor;
    //       /* Put back the unscaled value of Lambda. */
    //       kbp[kbp_context]->Lambda *= scaling_factor;
    //    }
    //
    //    /* Assign the best e-value field. Here the best e-value will always be
    //       attained for the first HSP in the list. Check that the incoming
    //       HSP list is properly sorted by score. */
    //    ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    //    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //
    //    return 0;
    // }
    //
    // ```
    let hsps: Vec<_> = input
        .iter()
        .enumerate()
        .map(|(source_index, h)| LinkedHsp {
            hsp: h.clone(),
            num: 0,
            evalue: match &gumbel {
                Some(gumbel) => blast_spouge_stoe(
                    h.score,
                    &ka,
                    gumbel,
                    query_lengths[h.context],
                    subject_length,
                ),
                None => {
                    parameters[h.context].search_space as f64
                        * (-ka.lambda * h.score as f64 + ka.k.ln()).exp()
                }
            },
            source_index,
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1737-1749
    // ```c++
    //  * @return TRUE if OK, FALSE otherwise.
    //  */
    // static double
    // s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
    // {
    //     int index = 0;
    //     double best_evalue = (double) INT4_MAX;
    //
    //     for (index=0; index<hsp_list->hspcnt; index++)
    //        best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
    //
    //     return best_evalue;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1900-1904
    // ```c++
    //       HSP list is properly sorted by score. */
    //    ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    //    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //
    //    return 0;
    // ```
    let best_evalue = if hsps.is_empty() {
        0.0
    } else {
        hsps.iter()
            .fold(i32::MAX as f64, |best, h| h.evalue.min(best))
    };
    Ok(LinkedHspList { hsps, best_evalue })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-669
// ```c++
// s_Blast_HSPListReapByPrelimEvalue(BlastHSPList* hsp_list, const BlastHitSavingParameters* hit_params)
// {
//    BlastHSP* hsp;
//    BlastHSP** hsp_array;
//    Int4 hsp_cnt = 0;
//    Int4 index;
//    double cutoff;
//
//    if (hsp_list == NULL)
//       return 0;
//
//    cutoff = hit_params->prelim_evalue;
//
//    hsp_array = hsp_list->hsp_array;
//    for (index = 0; index < hsp_list->hspcnt; index++) {
//       hsp = hsp_array[index];
//
//       ASSERT(hsp != NULL);
//
//       if (hsp->evalue > cutoff) {
//          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
//       } else {
//          if (index > hsp_cnt)
//             hsp_array[hsp_cnt] = hsp_array[index];
//          hsp_cnt++;
//       }
//    }
// ```
pub fn reap_preliminary(list: &mut LinkedHspList, options: &ResolvedOptions) -> Vec<LinkedHsp> {
    // `params->prelim_evalue = evalue` without a Gumbel block (blast_parameters.c:931-935).
    let cutoff = options.evalue
        * if options.composition > 0 && !options.old_fsc {
            5.0
        } else {
            1.0
        };
    let mut deleted = Vec::new();
    list.hsps.retain(|h| {
        if h.evalue > cutoff {
            deleted.push(h.clone());
            false
        } else {
            true
        }
    });
    deleted
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-2008
// ```c++
// Int2 Blast_HSPListReapByEvalue(BlastHSPList* hsp_list,
//         const BlastHitSavingOptions* hit_options)
// {
//    BlastHSP* hsp;
//    BlastHSP** hsp_array;
//    Int4 hsp_cnt = 0;
//    Int4 index;
//    double cutoff;
//
//    if (hsp_list == NULL)
//       return 0;
//
//    cutoff = hit_options->expect_value;
//
//    hsp_array = hsp_list->hsp_array;
//    for (index = 0; index < hsp_list->hspcnt; index++) {
//       hsp = hsp_array[index];
//
//       ASSERT(hsp != NULL);
//
//       if (hsp->evalue > cutoff) {
//          hsp_array[index] = Blast_HSPFree(hsp_array[index]);
//       } else {
//          if (index > hsp_cnt)
//             hsp_array[hsp_cnt] = hsp_array[index];
//          hsp_cnt++;
//       }
//    }
//
//    hsp_list->hspcnt = hsp_cnt;
//
//    return 0;
// }
// ```
pub fn reap_final(list: &mut LinkedHspList, options: &ResolvedOptions) -> Vec<LinkedHsp> {
    let mut deleted = Vec::new();
    list.hsps.retain(|h| {
        if h.evalue > options.evalue {
            deleted.push(h.clone());
            false
        } else {
            true
        }
    });
    deleted
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:102-114
// ```c++
// s_HSPListNormalizeScores(BlastHSPList * hsp_list,
//                          double lambda,
//                          double logK,
//                          double scoreDivisor)
// {
//     int hsp_index;
//     for(hsp_index = 0; hsp_index < hsp_list->hspcnt; hsp_index++) {
//         BlastHSP * hsp = hsp_list->hsp_array[hsp_index];
//
//         hsp->score = (Int4)BLAST_Nint(((double) hsp->score) / scoreDivisor);
//         /* Compute the bit score using the newly computed scaled score. */
//         hsp->bit_score = (hsp->score*lambda*scoreDivisor - logK)/NCBIMATH_LN2;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1917-1934
// ```c++
//
//    kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
//
//    for (index=0; index<hsp_list->hspcnt; index++) {
//       hsp = hsp_list->hsp_array[index];
//       ASSERT(hsp != NULL);
// #if 0
//       ASSERT(sbp->round_down == FALSE || (hsp->score & 1) == 0);
// #endif
//       hsp->bit_score =
//          (hsp->score*kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) /
//          NCBIMATH_LN2;
//    }
//
//    return 0;
// }
//
// void Blast_HSPListPHIGetBitScores(BlastHSPList* hsp_list, BlastScoreBlk* sbp)
// ```
pub fn normalize_scores(list: &mut LinkedHspList, options: &ResolvedOptions) -> Vec<f64> {
    let ka = lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    let scale = if options.composition == 2 { 32.0 } else { 1.0 };
    list.hsps
        .iter_mut()
        .map(|h| {
            if scale != 1.0 {
                let value = h.hsp.score as f64 / scale;
                h.hsp.score = (value + if value >= 0.0 { 0.5 } else { -0.5 }) as i32;
            }
            (h.hsp.score as f64 * (ka.lambda / scale) * scale - ka.k.ln()) / std::f64::consts::LN_2
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:511-532
// ```c++
//
//         /* Initialize the query */
//         if (program_number == eBlastTypeBlastx && kIsOutOfFrame) {
//             Int4 context = hsp->context - hsp->context % CODON_LENGTH;
//             Int4 context_offset = query_info->contexts[context].query_offset;
//             query = query_blk->oof_sequence + CODON_LENGTH + context_offset;
//             query_nomask = query_blk->oof_sequence + CODON_LENGTH + context_offset;
//         } else {
//             query = query_blk->sequence +
//                 query_info->contexts[hsp->context].query_offset;
//             query_nomask = query_blk->sequence_nomask +
//                 query_info->contexts[hsp->context].query_offset;
//         }
//
//         /* Translate subject if needed. */
//         if (program_number == eBlastTypeTblastn) {
//             const Uint1* target_sequence = Blast_HSPGetTargetTranslation(target_t, hsp, NULL);
//             status = Blast_HSPGetNumIdentitiesAndPositives(query, target_sequence, hsp, scoring_options, 0, sbp);
//         }
//         else
//             status = Blast_HSPGetNumIdentitiesAndPositives(query_nomask, subject, hsp, scoring_options, 0, sbp);
//
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:767-833
// ```c++
//    num_ident = 0;
//    align_length = 0;
//
//    if(NULL != sbp)
//    {
// 	   if(sbp->protein_alphabet)
// 		   matrix = sbp->matrix->data;
//    }
//
//    if (!hsp->gap_info) {
//       /* Ungapped case. Check that lengths are the same in query and subject,
//          then count number of matches. */
//       if (q_length != s_length)
//          return -1;
//       align_length = q_length;
//       for (i=0; i<align_length; i++) {
//          if (*q == *s)
//             num_ident++;
//          else if (NULL != matrix) {
//         	 if (matrix[*q][*s] > 0)
//         		 num_pos ++;
//              }
//          q++;
//          s++;
//       }
//    	}
//     else {
//       Int4 index;
//       GapEditScript* esp = hsp->gap_info;
//       for (index=0; index<esp->size; index++)
//       {
//          align_length += esp->num[index];
//          switch (esp->op_type[index]) {
//          case eGapAlignSub:
//             for (i=0; i<esp->num[index]; i++) {
//                if (*q == *s) {
//                   num_ident++;
//                }
//                else if (NULL != matrix) {
//             	   if (matrix[*q][*s] > 0)
//             		   num_pos ++;
//                }
//                q++;
//                s++;
//             }
//             break;
//          case eGapAlignDel:
//             s += esp->num[index];
//             break;
//          case eGapAlignIns:
//             q += esp->num[index];
//             break;
//          default:
//             s += esp->num[index];
//             q += esp->num[index];
//             break;
//          }
//       }
//    }
//
//    if (align_length_ptr) {
//        *align_length_ptr = align_length;
//    }
//    *num_ident_ptr = num_ident;
//
//    if(NULL != matrix)
// 	   *num_pos_ptr = num_pos + num_ident;
// ```
pub fn identity_and_positive(
    h: &super::kappa::RedoneHsp,
    batch: &PreparedQueryBatch,
    subject: &[u8],
) -> (usize, usize) {
    let p = &h.hsp;
    let c = &batch.contexts[p.context];
    let query = &batch.sequence_start_nomask[1 + c.offset..1 + c.offset + c.length];
    let mut qi = p.q_start as usize;
    let mut si = p.s_start as usize;
    let mut identity = 0;
    let mut positive = 0;
    for op in &h.edit_script {
        match *op {
            crate::common::GapEditOp::Sub(n) => {
                for _ in 0..n {
                    let q = query[qi];
                    let s = subject[si];
                    if q == s {
                        identity += 1;
                        positive += 1;
                    } else if h.identity_matrix.as_ref().map_or_else(
                        || crate::utils::matrix::blosum62_score(q, s) > 0,
                        |matrix| matrix.score(q, s) > 0,
                    ) {
                        positive += 1;
                    }
                    qi += 1;
                    si += 1;
                }
            }
            crate::common::GapEditOp::Del(n) => si += n as usize,
            crate::common::GapEditOp::Ins(n) => qi += n as usize,
            _ => {
                unreachable!("BLASTX in-frame edit script")
            }
        }
    }
    (identity, positive)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
// ```c++
// BLAST_KarlinStoE_simple(Int4 S,
//       Blast_KarlinBlk* kbp,
//       Int8  searchsp)   /* size of search space. */
// {
//    double   Lambda, K, H; /* parameters for Karlin statistics */
//
//    Lambda = kbp->Lambda;
//    K = kbp->K;
//    H = kbp->H;
//    if (Lambda < 0. || K < 0. || H < 0.) {
//       return -1.;
//    }
//
//    return (double) searchsp * exp((double)(-Lambda * S) + kbp->logK);
// }
// ```
pub(crate) fn evaluate_ungapped(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject_length: i32,
    prepared: Option<([i32; 2], f64)>,
) -> Result<LinkedHspList> {
    if let Some(link) = link_parameters(options) {
        if link.longest_intron <= 0 {
            return Ok(super::even_gap::link_even(
                input,
                batch,
                parameters,
                subject_length,
                prepared.expect("per-subject even-gap cutoffs"),
            ));
        }
        let lengths: Vec<_> = batch.contexts.iter().map(|c| c.length as i32).collect();
        let ka: Vec<_> = parameters.iter().map(|p| p.ungapped).collect();
        return link_uneven(
            input,
            &lengths,
            parameters,
            subject_length,
            &ka,
            None,
            &link,
        );
    }
    let hsps: Vec<_> = input
        .iter()
        .enumerate()
        .map(|(source_index, h)| {
            let p = &parameters[h.context];
            LinkedHsp {
                hsp: h.clone(),
                num: 0,
                evalue: p.search_space as f64
                    * (-p.ungapped.lambda * h.score as f64 + p.ungapped.k.ln()).exp(),
                source_index,
            }
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1900-1904
    // ```c++
    //       HSP list is properly sorted by score. */
    //    ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
    //    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //
    //    return 0;
    // ```
    let best_evalue = if hsps.is_empty() {
        0.0
    } else {
        hsps.iter()
            .fold(i32::MAX as f64, |best, h| h.evalue.min(best))
    };
    Ok(LinkedHspList { hsps, best_evalue })
}
