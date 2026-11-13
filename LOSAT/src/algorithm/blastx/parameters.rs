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
use crate::stats::spouge::{blast_spouge_etos, lookup_protein_gumbel_params, BlastGumbelBlk};
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

/// No NCBI counterpart: whether `LOSAT_X_BXKAONCE` is on (a strict-set switch, on unless `=0`):
/// a query batch's Karlin blocks computed once when the batch is not split (`search.rs`), and the
/// BLOSUM62 ideal block once per process (`ideal_block`).
pub(crate) fn ka_once_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXKAONCE"))
}

/// The BLOSUM62 ideal Karlin block of the standard composition against itself, as `score_block`
/// computes it for each call without `LOSAT_X_BXKAONCE`.
fn ideal_block_computed(std: &[f64; crate::utils::matrix::BLASTAA_SIZE]) -> KarlinParams {
    let ideal_profile =
        compute_score_freq_profile_for_matrix(std, std, -4, 11, ScoringMatrix::Blosum62);
    compute_karlin_params_ungapped(&ideal_profile).expect("BLOSUM62 ideal KA block")
}

/// No NCBI counterpart (`LOSAT_X_BXKAONCE`): the ideal block is a constant of the program (no
/// input), so it is computed on first use and shared; the value `ideal_block_computed` gives.
/// BLOSUM62's only: another matrix's searches compute theirs in `score_block_of`.
fn ideal_block() -> KarlinParams {
    static IDEAL: std::sync::OnceLock<KarlinParams> = std::sync::OnceLock::new();
    *IDEAL.get_or_init(|| ideal_block_computed(&compute_std_aa_composition()))
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
/// `score_block` with the search's matrix: BLOSUM62's is `score_block`, another matrix's
/// `score_block_of`.
pub fn score_block_for(
    batch: &PreparedQueryBatch,
    matrix: ScoringMatrix,
) -> Vec<ContextParameters> {
    if matrix == ScoringMatrix::Blosum62 {
        score_block(batch)
    } else {
        score_block_of(batch, matrix)
    }
}

/// `score_block` with another matrix than BLOSUM62: the contexts' score profiles over the
/// matrix's score range (`sbp->loscore`, `sbp->hiscore`) and its ideal block
/// (`ideal_karlin_params_of`, the block TBLASTN's and BLASTP's searches with that matrix take).
/// Not inlined, so that BLOSUM62's `score_block` keeps its code.
#[cold]
#[inline(never)]
fn score_block_of(batch: &PreparedQueryBatch, matrix: ScoringMatrix) -> Vec<ContextParameters> {
    use crate::algorithm::tblastx::lookup::{ideal_karlin_params_of, matrix_score_bounds};
    let std = compute_std_aa_composition();
    let ideal = ideal_karlin_params_of(matrix);
    let (low, high) = matrix_score_bounds(matrix);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3127-3130 (the alpha and
    // beta of an ungapped search, which the length adjustment reads, blast_setup.c:815)
    // ```c
    //    else if (num_values > 0) {
    //      (*alpha) = alpha_arr[0];
    //      (*beta) = beta_arr[0];
    // ```
    // The matrix's ungapped row (BLOSUM62's are the computed block's defaults), as TBLASTX's
    // contexts take them (`tblastx/blast_engine/run_impl.rs`).
    let row = crate::stats::tables::lookup_protein_params_ungapped(matrix);
    batch
        .contexts
        .iter()
        .map(|c| {
            let mut seq = vec![0];
            seq.extend_from_slice(&batch.sequence_start[1 + c.offset..1 + c.offset + c.length]);
            seq.push(0);
            let comp = compute_aa_composition(&seq, c.length);
            let profile = compute_score_freq_profile_for_matrix(&comp, &std, low, high, matrix);
            let ka = compute_karlin_params_ungapped(&profile).ok().map(|p| {
                let mut ka = apply_check_ideal(p, ideal);
                ka.alpha = row.alpha;
                ka.beta = row.beta;
                ka
            });
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

pub fn score_block(batch: &PreparedQueryBatch) -> Vec<ContextParameters> {
    let std = compute_std_aa_composition();
    let ideal = if ka_once_on() {
        ideal_block()
    } else {
        ideal_block_computed(&std)
    };
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:913-920
// ```c
//       /* Set the database length for new FSC */
//       if (sbp->gbp) {
//           Int8 dbl = total_length;
//           /* Override the database length or set one (e.g., blast2seq). */
//           if (eff_len_options->db_length) {
//               dbl = eff_len_options->db_length;
//           }
//           sbp->gbp->db_length =
// ```
/// The database length of the Gumbel block: `-dbsize` when it is not 0 (also below 0, as
/// NCBI), else the real total `db_length`; not divided by 3, since BLASTX's subject is not
/// translated (blast_setup.c:921-922).
pub fn gumbel_db_length(options: &ResolvedOptions, db_length: i64) -> i64 {
    if options.db_length != 0 {
        options.db_length
    } else {
        db_length
    }
}

/// The gapped Karlin block of the search's matrix and gap costs (`kbp_gap`, from NCBI's tables);
/// an ungapped search reads none (its gap costs need not be in the matrix's table, `args.rs`), and
/// gets the default block, which it never reads.
pub(crate) fn gapped_karlin(options: &ResolvedOptions) -> KarlinParams {
    if options.gapped {
        lookup_protein_params(&ProteinScoringSpec {
            matrix: options.scoring,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        })
    } else {
        KarlinParams::default()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
// ```c++
//     if (args.Exist(kArgEffSearchSpace) && args[kArgEffSearchSpace]) {
//         CNcbiEnvironment env;
//         env.Set("OLD_FSC", "true");
//         opt.SetEffectiveSearchSpace(args[kArgEffSearchSpace].AsInt8());
//     }
// ```
/// The Gumbel block of a gapped search, none with `-searchsp` (at any value) or `OLD_FSC`, nor in
/// an ungapped search (`ResolvedOptions::has_gumbel_block`).
pub fn gumbel_block(options: &ResolvedOptions, db_length: i64) -> Option<BlastGumbelBlk> {
    if !options.has_gumbel_block() {
        return None;
    }
    lookup_protein_gumbel_params(
        &ProteinScoringSpec {
            matrix: options.scoring,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        },
        gumbel_db_length(options, db_length),
    )
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
/// LOSAT's refusal of a search whose effective search space overflows NCBI's `Int8`
/// (`db_length` letters of the database, or of `-dbsize`, against a context of
/// `query_length` residues): NCBI's arithmetic overflows there (as BLASTP's `-dbsize`).
#[derive(Debug)]
pub struct SearchSpaceOverflow {
    pub db_length: i64,
    pub query_length: usize,
    /// The program of the refusal (`args::ProgramName::upper`).
    pub program: &'static str,
}

impl std::fmt::Display for SearchSpaceOverflow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "a database length of {} letters with a query of {} letters overflows NCBI BLAST+'s 64-bit effective search space, which is not supported by LOSAT's {}",
            self.db_length, self.query_length, self.program
        )
    }
}

impl std::error::Error for SearchSpaceOverflow {}

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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-720
    // ```c
    //    /* use overriding value from effective lengths options or the real value
    //       from effective lengths parameters. */
    //    if (eff_len_options->db_length > 0)
    //       db_length = eff_len_options->db_length;
    //    else
    // ```
    // `db_length` is the real total; `-dbsize` above 0 replaces it here (never divided by
    // 3: BLASTX's subject is not translated).
    let db_length = if options.db_length > 0 {
        options.db_length
    } else {
        db_length
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:729-732
    // ```c
    //    if (db_length == 0 &&
    //        !BlastEffectiveLengthsOptions_IsSearchSpaceSet(eff_len_options)) {
    //       return 0;
    //    }
    // ```
    let options_space = options.effective_search_space;
    if db_length == 0 && fixed_spaces.is_none() && options_space == 0 {
        return Ok(());
    }
    let gapped = gapped_karlin(options);
    for (i, (c, p)) in batch.contexts.iter().zip(parameters).enumerate() {
        if !p.valid {
            p.length_adjustment = 0;
            p.search_space = fixed_spaces.map_or(options_space, |spaces| spaces[i]);
            continue;
        }
        let ka = if options.gapped { gapped } else { p.ungapped };
        let adj = compute_length_adjustment_ncbi(c.length as i64, db_length, db_count, &ka)
            .length_adjustment;
        p.length_adjustment = i32::try_from(adj)?;
        p.search_space = match fixed_spaces {
            Some(spaces) => spaces[i],
            None if options_space != 0 => options_space,
            None => (db_length - db_count * adj)
                .max(1)
                .checked_mul(c.length as i64 - adj)
                .ok_or_else(|| SearchSpaceOverflow {
                    db_length,
                    query_length: c.length,
                    program: options.program().upper,
                })?,
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
/// The numbers of the subjects that the parameters of a search read (`Subjects::numbers`), as
/// NCBI's sequence source answers them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubjectNumbers {
    /// The total length of the statistics (blast_setup.c:908-911): `STATS_TOTLEN` or
    /// `BlastSeqSrcGetTotLen` (`SeqDb::statistics_numbers`), the sum of a FASTA set's records.
    pub total: i64,
    /// The number of sequences of the statistics (blast_setup.c:925-928): `STATS_NSEQ` or
    /// `GetNumSeqs`, 0 for an alias file's `NSEQ 0`; a FASTA set's records.
    pub count: i64,
    /// `BlastSeqSrcGetMinSeqLen` (blast_setup.c:969-973): a database's `GetMinLength`
    /// (`SeqDb::min_length`), a FASTA set's shortest record raised to `BLAST_SEQSRC_MINLENGTH`.
    pub shortest: usize,
    /// `BlastSeqSrcGetAvgSeqLen` (blast_engine.c:1372-1374, the word parameters): a database's
    /// `GetTotalLength / MAX(1, GetNumSeqs)` (`SeqDb::average_length`, no `STATS_*` key), a
    /// FASTA set's total over its records.
    pub average: i32,
}

impl SubjectNumbers {
    /// The numbers of one FASTA subject of `length` residues.
    #[cfg(test)]
    pub fn one(length: usize) -> Self {
        Self {
            total: length as i64,
            count: 1,
            shortest: length.max(super::subjects::BLAST_SEQSRC_MINLENGTH),
            average: length as i32,
        }
    }
}

pub fn subject_parameters(
    batch: &PreparedQueryBatch,
    parameters: &mut [ContextParameters],
    options: &ResolvedOptions,
    numbers: SubjectNumbers,
    fixed_spaces: Option<&[i64]>,
) -> Result<i32> {
    let SubjectNumbers {
        total: db_length,
        count: db_count,
        shortest,
        average,
    } = numbers;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // ```c
    //    if (sbp->gbp) {
    //        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
    //        ...
    //    } else {
    //        min_subject_length = (Int4) (total_length/num_seqs);
    //    }
    //
    //    if(min_subject_length <=0) {
    // 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
    //    }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-217
    // ```c
    //     case BLASTERR_SUBJECT_LENGTH_INVALID:
    //         new_msg->message = strdup("The average subject length is too short");
    //         new_msg->severity = eBlastSevFatal;
    // ```
    // BLASTX's subject is not translated, so the shortest length is not divided by 3. Without the
    // Gumbel block (`ResolvedOptions::has_gumbel_block`: `-ungapped`, `-searchsp`, `OLD_FSC`) the
    // statistics' average divides by their number of sequences, which an alias file's `NSEQ 0`
    // makes 0: NCBI crashes there (SIGFPE), LOSAT refuses (`common::database::check_statistics_count`,
    // DECISIONS 14 of line RG). A set of subjects whose average is below one letter fails the
    // batch's setup: the engine error after the start of the outfmt 0 report, exit status 3, as
    // every other error of this setup (`native::database_search_error`; NCBI 2.17.0 blastx
    // `-subject` of three empty records, or of `M`, `` and ``, with `-searchsp`, `OLD_FSC` or
    // `-ungapped -comp_based_stats 0`, env parity probe f_probe of 2026-11-12).
    let gumbel_block_used = options.has_gumbel_block();
    if !gumbel_block_used && db_count == 0 && db_length > 0 {
        return Err(crate::algorithm::common::database::refusal(
            "the statistics count 0 sequences of the database (an alias file's NSEQ 0), by which NCBI BLAST+ divides its length",
            "BLASTX",
        ));
    }
    let n = if gumbel_block_used {
        i32::try_from(shortest)?
    } else {
        i32::try_from(db_length / db_count.max(1))?
    };
    ensure!(n > 0, "The average subject length is too short");
    effective_lengths(
        batch,
        parameters,
        options,
        db_length,
        db_count,
        fixed_spaces,
    )?;
    let gapped = gapped_karlin(options);
    let gumbel = gumbel_block(options, db_length);
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
    // The average context length is that of the query batch (a reader batch searched in several
    // units has one: `reader_batch.rs`, c++/src/algo/blast/core/blast_parameters.c:954-957).
    let avg_q = i32::try_from(batch.batch_concatenated_length() / batch.batch_context_count())?;
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
            cutoff_score_sum_stats(c.length as i32, average, decay, &ka)
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
        gapped_x_drops(options, &gapped).0
    } else {
        0
    })
}

/// The preliminary and final gapped X-drops in raw scores (`gap_x_dropoff`,
/// `gap_x_dropoff_final`) of the search's X-drops in bits (`-xdrop_gap`, `-xdrop_gap_final` or
/// BLASTX's defaults) at the gapped Lambda `gapped` (`gapped_karlin`; the contexts share it, the
/// smallest one), converted as the protein programs convert them (`core::x_drop`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:454-463 (quoted at
/// `core::x_drop::gapped_x_dropoffs`)
pub(crate) fn gapped_x_drops(options: &ResolvedOptions, gapped: &KarlinParams) -> (i32, i32) {
    crate::core::x_drop::gapped_x_dropoffs(
        options.gap_x_dropoff,
        options.gap_x_dropoff_final,
        gapped.lambda,
    )
}

/// The refusal, before the search, of a gapped X-drop that NCBI's kernels cannot size
/// (`core::x_drop::check_gapped_x_dropoffs`); an ungapped search has none.
pub(crate) fn check_gapped_x_drops(options: &ResolvedOptions) -> Result<()> {
    if !options.gapped {
        return Ok(());
    }
    // The final X-drop reaches the unscaled `ALIGN_EX` only on the ordinary traceback
    // (composition 0, `runtime.rs`).
    crate::core::x_drop::check_gapped_x_dropoffs(
        "BLASTX",
        (options.gap_x_dropoff, options.gap_x_dropoff_final),
        gapped_x_drops(options, &gapped_karlin(options)),
        (options.gap_open, options.gap_extend),
        options.composition == 0,
    )
}
