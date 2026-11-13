//! The database numbers of NCBI's statistics setup (`BLAST_GapAlignSetUp` and
//! `BLAST_CalcEffLengths`, c++/src/algo/blast/core/blast_setup.c), with NCBI's types: the one
//! owner of the effective lengths of a query context and of the "average subject length is too
//! short" test for the six search programs (line RG D1a-3).

use super::length_adjustment::compute_length_adjustment_ncbi;
use super::tables::KarlinParams;

/// The effective lengths of one query context (`BlastContextInfo::length_adjustment`,
/// `eff_searchsp`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ContextEffLengths {
    /// `length_adjustment`: NCBI's `Int4`.
    pub(crate) length_adjustment: i32,
    /// `eff_searchsp`: NCBI's `Int8`.
    pub(crate) eff_searchsp: i64,
}

/// LOSAT's refusal of an effective search space that overflows NCBI's `Int8` (signed overflow,
/// undefined in NCBI's C): `db_length` letters (after `-dbsize` and the division of a translated
/// subject) with a query context of `query_length` letters, in the search of `program`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SearchSpaceOverflow {
    pub(crate) db_length: i64,
    pub(crate) query_length: i64,
    pub(crate) program: String,
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-740
// ```c
//    /* use overriding value from effective lengths options or the real value
//       from effective lengths parameters. */
//    if (eff_len_options->db_length > 0)
//       db_length = eff_len_options->db_length;
//    else
//       db_length = eff_len_params->real_db_length;
//
//    /* If database (subject) length is not available at this stage, and
//     * overriding value of effective search space is not provided by user,
//     * do nothing.
//     ...
//    if (db_length == 0 &&
//        !BlastEffectiveLengthsOptions_IsSearchSpaceSet(eff_len_options)) {
//       return 0;
//    }
//
//    if (Blast_SubjectIsTranslated(program_number))
//       db_length = db_length/3;
//
//    if (eff_len_options->dbseq_num > 0)
//       db_num_seqs = eff_len_options->dbseq_num;
//    else
//       db_num_seqs = eff_len_params->real_num_seqs;
// ```
// The command lines have no `dbseq_num` option (it is 0): the number is the real one.
/// The database length and number of sequences that `BLAST_CalcEffLengths` computes the
/// contexts with: `dbsize` (`-dbsize`) when above 0, else `real_db_length` (the setup's total,
/// or a subject's length), divided by 3 for a translated subject, and `real_num_seqs`. `None`
/// where NCBI returns before the contexts and leaves their values as they are: a length of 0
/// and no search space option (`searchsp_set`, `-searchsp`, or a query chunk's spaces).
pub(crate) fn eff_lengths_database(
    real_db_length: i64,
    real_num_seqs: i32,
    dbsize: i64,
    searchsp_set: bool,
    translated_subject: bool,
) -> Option<(i64, i32)> {
    let db_length = if dbsize > 0 { dbsize } else { real_db_length };
    if db_length == 0 && !searchsp_set {
        return None;
    }
    let db_length = if translated_subject {
        db_length / 3
    } else {
        db_length
    };
    Some((db_length, real_num_seqs))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:778-848
// ```c
//       Int8 effective_search_space =
//           s_GetEffectiveSearchSpaceForContext(eff_len_options, index,
//                                               blast_message);
//       ...
//       if (query_info->contexts[index].is_valid &&
//           ((query_length = query_info->contexts[index].query_length) > 0) ) {
//          ...
//          BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK,
//                                        alpha/kbp->Lambda, beta,
//                                        query_length, db_length,
//                                        db_num_seqs, &length_adjustment);
//
//          if (effective_search_space == 0) {
//             ...
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
// `db_num_seqs` is NCBI's `Int4`, negative for an alias file's negative `NSEQ` (or one of 2^31
// or more): the length adjustment is computed with it as a double, and `(Int8)db_num_seqs *
// length_adjustment` cannot overflow. The subtraction (a negative product with a length near
// 2^63) and the product (a database of about 10^17 letters and more) can: that signed overflow
// is LOSAT's refusal (`SearchSpaceOverflow`; NCBI BLAST+ 2.17.0 prints a wrapped space, probe
// p5_big of line RG D1a-3).
/// `BLAST_CalcEffLengths` for one valid query context of `query_length` (> 0) letters, after
/// `eff_lengths_database`: the length adjustment for the Karlin-Altschul block `ka` (`K`,
/// `Lambda` and the context's `alpha` and `beta`) and the effective search space, which is
/// `options_searchsp` when that is not 0 (`-searchsp`, a query chunk's space).
pub(crate) fn calc_eff_lengths(
    query_length: i64,
    db_length: i64,
    db_num_seqs: i32,
    ka: &KarlinParams,
    options_searchsp: i64,
    program: &str,
) -> Result<ContextEffLengths, SearchSpaceOverflow> {
    let adjustment =
        compute_length_adjustment_ncbi(query_length, db_length, i64::from(db_num_seqs), ka)
            .length_adjustment;
    // `BLAST_ComputeLengthAdjustment` writes an `Int4`.
    let length_adjustment = adjustment as i32;
    if options_searchsp != 0 {
        return Ok(ContextEffLengths {
            length_adjustment,
            eff_searchsp: options_searchsp,
        });
    }
    let overflow = || SearchSpaceOverflow {
        db_length,
        query_length,
        program: program.to_string(),
    };
    let mut effective_db_length = db_length
        .checked_sub(i64::from(db_num_seqs) * i64::from(length_adjustment))
        .ok_or_else(overflow)?;
    if effective_db_length <= 0 {
        effective_db_length = 1;
    }
    let eff_searchsp = effective_db_length
        .checked_mul(query_length - i64::from(length_adjustment))
        .ok_or_else(overflow)?;
    Ok(ContextEffLengths {
        length_adjustment,
        eff_searchsp,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-220
// ```c++
//     CEffectiveSearchSpacesMemento eff_memento(m_Options);
//     SplitQuery_SetEffectiveSearchSpace(m_Options, m_QueryFactory,
//                                        m_InternalData);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:155-181
// ```c++
//     // If the effective search options have been set, we don't need to
//     // recompute those...
//     if (options->GetEffectiveSearchSpace() != 0) {
//         return;
//     }
//
//     const BlastSeqSrc* seqsrc = full_data->m_SeqSrc->GetPointer();
//     Int8 total_length = BlastSeqSrcGetTotLenStats(seqsrc);
//     if (total_length <= 0)
//         total_length = BlastSeqSrcGetTotLen(seqsrc);
//     Int4 num_seqs = BlastSeqSrcGetNumSeqsStats(seqsrc);
//     if (num_seqs <= 0)
//         num_seqs = BlastSeqSrcGetNumSeqs(seqsrc);
//
//     CEffectiveSearchSpaceCalculator calc(full_query_fact, *options,
//                                          num_seqs, total_length,
//                                          full_data->m_ScoreBlk->GetPointer());
//     BlastQueryInfo* qinfo = full_data->m_QueryInfo;
//     _ASSERT(qinfo);
//
//     vector<Int8> eff_searchsp;
//     for (size_t index = 0; index <= (size_t)qinfo->last_context; index++) {
//         eff_searchsp.push_back(calc.GetEffSearchSpaceForContext(index));
//     }
//     options->SetEffectiveSearchSpace(eff_searchsp);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/effsearchspace_calc.cpp:115-121
// ```c++
//     BlastEffectiveLengthsParametersNew(opts_memento->m_EffLenOpts,
//                                        db_num_bases, db_num_seqs,
//                                        &eff_len_params);
//
//     Int2 status =
//         BLAST_CalcEffLengths(m_Program, opts_memento->m_ScoringOpts,
//                              eff_len_params, sbp, m_QueryInfo, NULL);
// ```
// The calculator's query information is a new one of the whole query (its spaces start at 0),
// and its setup has no fallback to the length of OID 0: with a total length of 0 (an alias
// file's `LENGTH 0`, a list whose scan totals nothing) `BLAST_CalcEffLengths` returns before the
// contexts (`eff_lengths_database`) unless `-dbsize` is given, and every option space is 0 (the
// options then count as unset, `BlastEffectiveLengthsOptions_IsSearchSpaceSet`); with a negative
// total (a `LENGTH` of 2^63 or more) the effective database length is 1 and a context's space is
// its query length. The setups of the preliminary stage read these spaces in place of their own
// (`s_GetEffectiveSearchSpaceForContext`, blast_setup.c:676-697 and 778-780), the memento gives
// the traceback and the reports the options' own spaces back. With a positive total the spaces
// are the setup's own (the same total, number and `-dbsize`).
/// The search space that NCBI's preliminary stage sets in the options for a context of
/// `query_length` letters of the whole query (`SplitQuery_SetEffectiveSearchSpace`): `searchsp`
/// (`-searchsp`, the options' own) when it is not 0, else the effective search space of the
/// statistics' total `statistics_total` and number `statistics_count` before their fallback to OID
/// 0 (`DbStatistics::gumbel_db_length`, `request_count`) with `dbsize` (`-dbsize`), or 0 where
/// `BLAST_CalcEffLengths` returns before the contexts. The setups of the preliminary stage and the
/// statistics of each subject (`subject_eff_lengths`) take it as the context's option space
/// (`calc_eff_lengths`' `options_searchsp`); an invalid or empty context has 0.
#[allow(clippy::too_many_arguments)]
pub(crate) fn preliminary_option_space(
    query_length: i64,
    statistics_total: i64,
    statistics_count: i32,
    dbsize: i64,
    searchsp: i64,
    translated_subject: bool,
    ka: &KarlinParams,
    program: &str,
) -> Result<i64, SearchSpaceOverflow> {
    if searchsp != 0 {
        return Ok(searchsp);
    }
    if query_length <= 0 {
        return Ok(0);
    }
    match eff_lengths_database(
        statistics_total,
        statistics_count,
        dbsize,
        false,
        translated_subject,
    ) {
        None => Ok(0),
        Some((db_length, num_seqs)) => {
            Ok(calc_eff_lengths(query_length, db_length, num_seqs, ka, 0, program)?.eff_searchsp)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1001-1014
// ```c
// Int2 BLAST_OneSubjectUpdateParameters(EBlastProgramType program_number,
//                     Uint4 subject_length,
//                     const BlastScoringOptions* scoring_options,
//                     BlastQueryInfo* query_info,
//                     const BlastScoreBlk* sbp,
//                     BlastHitSavingParameters* hit_params,
//                     BlastInitialWordParameters* word_params,
//                     BlastEffectiveLengthsParameters* eff_len_params)
// {
//    Int2 status = 0;
//    eff_len_params->real_db_length = subject_length;
//    if ((status = BLAST_CalcEffLengths(program_number, scoring_options,
//                                 eff_len_params, sbp, query_info, NULL)) != 0)
//       return status;
// ```
// `real_num_seqs` stays the setup's number (`DbStatistics::setup_count`: 1 after the fallback to
// OID 0, `STATS_NSEQ` with `STATS_TOTLEN`), and the options' spaces are the preliminary stage's
// (`preliminary_option_space`). A subject that NCBI searches has letters (blast_engine.c:1429).
/// The effective lengths of a valid context of `query_length` letters for the subject of
/// `subject_length` letters (nucleotides for a translated subject) that NCBI's preliminary search
/// computes the statistics of on its own (`BLAST_OneSubjectUpdateParameters`, a database length
/// of 0: `DbStatistics::per_subject`): `setup_count` sequences, `dbsize` (`-dbsize`) in place of
/// the subject's length, and the context's preliminary option space `option_space`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn subject_eff_lengths(
    query_length: i64,
    subject_length: i64,
    setup_count: i32,
    dbsize: i64,
    option_space: i64,
    translated_subject: bool,
    ka: &KarlinParams,
    program: &str,
) -> Result<ContextEffLengths, SearchSpaceOverflow> {
    let (db_length, num_seqs) = eff_lengths_database(
        subject_length,
        setup_count,
        dbsize,
        true,
        translated_subject,
    )
    .expect("eff_lengths_database returns the lengths when the option spaces count as set");
    calc_eff_lengths(query_length, db_length, num_seqs, ka, option_space, program)
}

/// The first overflow of NCBI's `Int8` search space that the statistics of a subject met in a
/// subject loop that cannot return an error (`subject_eff_lengths`): the loop notes it, and the
/// search returns it (LOSAT's refusal, `SearchSpaceOverflow`) after the loop.
#[derive(Debug, Default)]
pub(crate) struct SubjectSpaceOverflow(std::sync::OnceLock<SearchSpaceOverflow>);

impl SubjectSpaceOverflow {
    /// Notes `overflow` unless one was noted before.
    pub(crate) fn note(&self, overflow: SearchSpaceOverflow) {
        let _ = self.0.set(overflow);
    }

    /// The overflow that the loop noted, as the search's error.
    pub(crate) fn check(&self) -> Result<(), SearchSpaceOverflow> {
        match self.0.get() {
            Some(overflow) => Err(overflow.clone()),
            None => Ok(()),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:900-904, 969-980
// ```c
//    Uint4 max_subject_length;
//    Uint4 min_subject_length;
//    Int8 total_length = -1;
//    Int4 num_seqs = -1;
//    ...
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
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-217
// ```c
//     case BLASTERR_SUBJECT_LENGTH_INVALID:
//         new_msg->message = strdup("The average subject length is too short");
// ```
// `min_subject_length` is a `Uint4`, so `<= 0` is `== 0`: a negative quotient (a negative
// `NSEQ`) passes, as does one of 2^31, while one that is a multiple of 2^32 fails (an `NSEQ` of
// 2^31 as NCBI's `int`, or a `LENGTH` of 2^32 times the count; probes p5_ns_*, p5_avg_* of line
// RG D1a-3). The Gumbel block's shortest length is divided as a `Uint4` too.
/// The subject length of `BLAST_GapAlignSetUp` that the hit-saving parameters read (their
/// `Int4 avg_subj_length`), or `None` where the setup fails with `The average subject length is
/// too short`: with the Gumbel block (`gumbel_min`: `GetMinLength`, and whether the subject is
/// translated) the shortest length (divided by 3 for a translated subject), else the statistics'
/// total `total_length` over their number `num_seqs` as NCBI's `Int4`. `num_seqs` is not 0 there
/// (NCBI divides by 0: LOSAT refuses, `common::database::check_statistics_count`).
pub(crate) fn setup_subject_length(
    total_length: i64,
    num_seqs: i32,
    gumbel_min: Option<(i32, bool)>,
) -> Option<i32> {
    let min_subject_length: u32 = match gumbel_min {
        Some((min_length, translated)) => {
            let min_length = min_length as u32;
            if translated {
                min_length / 3
            } else {
                min_length
            }
        }
        None => {
            let average = total_length
                .checked_div(i64::from(num_seqs))
                .expect("a statistics count of 0 is refused before the setup");
            average as i32 as u32
        }
    };
    (min_subject_length != 0).then_some(min_subject_length as i32)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:755-760
// ```c
//    if (!sbp->kbp_gap)
//       gapped_calculation = FALSE;
//
//    if (options->do_sum_stats && gapped_calculation && avg_subj_length <= 0)
//        return 1;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:242-248
// ```c
//     /* Unknown error */
//     default:
//         {
//             char buf[512];
//             snprintf(buf, sizeof(buf) - 1, "Unknown error code %d", error_code);
// ```
// `BLAST_GapAlignSetUp` returns the status after the "too short" test (blast_setup.c:982-987):
// NCBI 2.17.0 blastx and tblastn with `-searchsp` over an alias file's `NSEQ 4294967295` (-1 as
// NCBI's `int`) end with `BLAST engine error: Unknown error code 1`, exit status 3 (probe run1 of
// line RG D1a-3, impl1).
/// Whether `BlastHitSavingParametersNew` fails in the setup: a gapped search with sum statistics
/// (`sum_stats`, `gapped`) and a setup subject length (`setup_subject_length`) that is not
/// positive (a negative average of a negative `NSEQ`, without the Gumbel block). The setup's
/// error is then `HIT_SAVING_PARAMETERS_ERROR`.
pub(crate) fn hit_saving_parameters_fail(
    sum_stats: bool,
    gapped: bool,
    avg_subj_length: i32,
) -> bool {
    sum_stats && gapped && avg_subj_length <= 0
}

/// The engine error of `hit_saving_parameters_fail` (`Blast_PerrorEx` of status 1).
pub(crate) const HIT_SAVING_PARAMETERS_ERROR: &str = "Unknown error code 1";

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI BLAST+ 2.17.0 probes of line RG D1a-3 (p5_ns_*, p5_avg_*, p5_min_*): which setups
    // fail with `The average subject length is too short`.
    #[test]
    fn too_short_is_the_uint4_of_ncbis_quotient() {
        let tfn = 10_937;
        assert_eq!(setup_subject_length(tfn, 16, None), Some(683));
        // NSEQ -1, -5, 2^32 - 1: negative quotients pass.
        assert_eq!(setup_subject_length(tfn, -1, None), Some(-10_937));
        assert_eq!(setup_subject_length(tfn, -5, None), Some(-2_187));
        // NSEQ 2^31 - 1, 2^31, 2^31 + 1, 3e9: a quotient of 0.
        for n in [
            i32::MAX,
            i32::MIN,
            (2_147_483_649u32) as i32,
            (3_000_000_000u32) as i32,
        ] {
            assert_eq!(setup_subject_length(tfn, n, None), None, "{n}");
        }
        // LENGTH 16 * 2^32 over 16 sequences: the `Int4` of 2^32 is 0; 2^31 and 2^32 + 5 run.
        assert_eq!(setup_subject_length(16 << 32, 16, None), None);
        assert_eq!(setup_subject_length(16 << 31, 16, None), Some(i32::MIN));
        assert_eq!(
            setup_subject_length(16 * ((1 << 32) + 5), 16, None),
            Some(5)
        );
        // blastn LENGTH 4294967296000000000 over 16: 268435456000000000 = 2^32 * 62500000.
        assert_eq!(
            setup_subject_length(4_294_967_296_000_000_000, 16, None),
            None
        );
        // The Gumbel block's MIN_SEQ_LENGTH: 1 and 2 fail for a translated subject, 3 runs.
        assert_eq!(setup_subject_length(tfn, 16, Some((1, true))), None);
        assert_eq!(setup_subject_length(tfn, 16, Some((2, true))), None);
        assert_eq!(setup_subject_length(tfn, 16, Some((3, true))), Some(1));
        assert_eq!(setup_subject_length(tfn, -1, Some((10, false))), Some(10));
        // A negative `int` is a large `Uint4`: 2^32 - 1 over 3.
        assert_eq!(
            setup_subject_length(tfn, 16, Some((-1, true))),
            Some(1_431_655_765)
        );
    }

    // probe p5_ns_2p32m1#10 of line RG D1a-3 (NSEQ -1 over tfn): the adjustment grows with the
    // length the negative count adds; `-searchsp` keeps the option's space.
    #[test]
    fn effective_lengths_take_ncbis_int4_count() {
        let ka = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        let positive = calc_eff_lengths(142, 3263, 23, &ka, 0, "BLASTP").unwrap();
        let negative = calc_eff_lengths(142, 3263, -1, &ka, 0, "BLASTP").unwrap();
        let adj = i64::from(negative.length_adjustment);
        assert_eq!(negative.eff_searchsp, (3263 + adj) * (142 - adj));
        assert_eq!(
            positive.eff_searchsp,
            (3263 - 23 * i64::from(positive.length_adjustment))
                * (142 - i64::from(positive.length_adjustment))
        );
        let fixed = calc_eff_lengths(142, 3263, -1, &ka, 1_000_000, "BLASTP").unwrap();
        assert_eq!(
            (fixed.length_adjustment, fixed.eff_searchsp),
            (negative.length_adjustment, 1_000_000)
        );
        let overflow = calc_eff_lengths(142, i64::MAX, -1, &ka, 0, "BLASTP").unwrap_err();
        assert!(overflow
            .to_string()
            .contains("64-bit effective search space"));
        assert_eq!(eff_lengths_database(0, 1, 0, false, false), None);
        assert_eq!(eff_lengths_database(0, 1, 0, true, true), Some((0, 1)));
        assert_eq!(eff_lengths_database(30, 2, 0, false, true), Some((10, 2)));
        assert_eq!(eff_lengths_database(30, 2, 99, false, false), Some((99, 2)));
    }

    // SURVEY_D1a3 §0.3 / impl2 of line RG D1a-3: the preliminary stage's option spaces
    // (`SplitQuery_SetEffectiveSearchSpace`) and the statistics of one subject.
    #[test]
    fn preliminary_option_spaces_follow_ncbis_split_query_setup() {
        let ka = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        // -searchsp is the options' own.
        assert_eq!(
            preliminary_option_space(142, 0, 23, 0, 777, false, &ka, "BLASTP").unwrap(),
            777
        );
        // A total of 0 without -dbsize: CalcEffLengths returns before the contexts.
        assert_eq!(
            preliminary_option_space(142, 0, 23, 0, 0, false, &ka, "BLASTP").unwrap(),
            0
        );
        // With -dbsize the walk's count, not the setup's fallback of 1.
        let dbsize =
            preliminary_option_space(142, 0, 23, 100_000, 0, false, &ka, "BLASTP").unwrap();
        assert_eq!(
            dbsize,
            calc_eff_lengths(142, 100_000, 23, &ka, 0, "BLASTP")
                .unwrap()
                .eff_searchsp
        );
        // A negative total: an effective length of 1, the context's query length.
        assert_eq!(
            preliminary_option_space(
                142,
                -9_223_372_036_854_775_807,
                23,
                0,
                0,
                false,
                &ka,
                "BLASTP"
            )
            .unwrap(),
            142
        );
        // A positive total: the setup's own space; an empty context has none.
        assert_eq!(
            preliminary_option_space(142, 3263, 23, 0, 0, false, &ka, "BLASTP").unwrap(),
            calc_eff_lengths(142, 3263, 23, &ka, 0, "BLASTP")
                .unwrap()
                .eff_searchsp
        );
        assert_eq!(
            preliminary_option_space(0, 3263, 23, 0, 0, false, &ka, "BLASTP").unwrap(),
            0
        );
        // One subject: its length (/3 translated) with the setup's count, an option space kept.
        let subject = subject_eff_lengths(142, 300, 1, 0, 0, true, &ka, "TBLASTN").unwrap();
        assert_eq!(
            subject,
            calc_eff_lengths(142, 100, 1, &ka, 0, "TBLASTN").unwrap()
        );
        let kept = subject_eff_lengths(142, 300, 7, 0, 481_440, false, &ka, "BLASTP").unwrap();
        assert_eq!(kept.eff_searchsp, 481_440);
        let fixed = subject_eff_lengths(142, 300, 1, 100_000, 0, false, &ka, "BLASTP").unwrap();
        assert_eq!(
            fixed,
            calc_eff_lengths(142, 100_000, 1, &ka, 0, "BLASTP").unwrap()
        );
    }
}
