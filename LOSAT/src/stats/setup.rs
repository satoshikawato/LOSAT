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
// `db_num_seqs` is NCBI's `Int4`, negative for an alias file's negative `NSEQ` (or one of 2^31
// or more): the length adjustment is computed with it as a double, and the product of `db_num_seqs`
// (cast to `Int8`) and `length_adjustment` cannot overflow. The subtraction (a negative product
// with a length near 2^63) and the product (a database of about 10^17 letters and more) can: that
// signed overflow is LOSAT's refusal (`SearchSpaceOverflow`; NCBI BLAST+ 2.17.0 prints a wrapped
// space, probe p5_big of line RG D1a-3).
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
    let length_adjustment = length_adjustment(query_length, db_length, db_num_seqs, ka);
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-824
/// The length adjustment of `BLAST_CalcEffLengths` for a valid context of `query_length` letters
/// (`calc_eff_lengths`' own, which the reports that keep only the search space read again):
/// `BLAST_ComputeLengthAdjustment` with NCBI's `Int4` count, written as an `Int4`.
pub(crate) fn length_adjustment(
    query_length: i64,
    db_length: i64,
    db_num_seqs: i32,
    ka: &KarlinParams,
) -> i32 {
    compute_length_adjustment_ncbi(query_length, db_length, i64::from(db_num_seqs), ka)
        .length_adjustment as i32
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-221
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:155-181
// NCBI reference (598d8ae6): c++/src/algo/blast/api/effsearchspace_calc.cpp:115-121
// NCBI reference (598d8ae6): c++/src/algo/blast/api/query_data.cpp:52-58 and
// blast_aux_priv.cpp:150-152 (`MakeLocalQueryData` caches the query data; the setup's
// `m_QueryInfo` is its query information).
// The calculator's query information is the setup's own of the whole query, whose `is_valid` the
// score block has set (`CSetupFactory::CreateScoreBlock`, before the search) and whose spaces no
// setup has computed yet (they start at 0): a context that is not valid keeps 0, and the
// calculator has no fallback to the length of OID 0: with a total length of 0 (an alias
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

// No NCBI counterpart: a call boundary; it computes no value NCBI computes.
/// Calls `f`, the statistics of each subject of a subject loop (or another call whose statistics
/// are not the setup's), in a function of its own that is never inlined, so that the loop's default
/// path keeps its code and pays one test of a flag per subject.
#[cold]
#[inline(never)]
pub(crate) fn apart<R>(f: impl FnOnce() -> R) -> R {
    f()
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:903-906,969-983
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-217
/// NCBI's error of a setup whose `min_subject_length` is 0 (`setup_subject_length`): the one
/// text of blastp's, psiblast's, rpsblast's, rpstblastn's, blastx's, tblastn's, blastn's and
/// tblastx's error (review C of line RG D1a-4, N1; reviews D and E, N1).
pub(crate) const SUBJECT_LENGTH_TOO_SHORT: &str = "The average subject length is too short";

// `min_subject_length` is a `Uint4`, so `<= 0` is `== 0`: a negative quotient (a negative
// `NSEQ`) passes, as does one of 2^31, while one that is a multiple of 2^32 fails (an `NSEQ` of
// 2^31 as NCBI's `int`, or a `LENGTH` of 2^32 times the count; probes p5_ns_*, p5_avg_* of line
// RG D1a-3). The Gumbel block's shortest length is divided as a `Uint4` too.
/// The subject length of `BLAST_GapAlignSetUp` that the hit-saving parameters read (their
/// `avg_subj_length`, an `Int4`), or `None` where the setup fails with `The average subject length is
/// too short`: with the Gumbel block (`gumbel_min`: `GetMinLength`, and whether the subject is
/// translated) the shortest length (divided by 3 for a translated subject), else the statistics'
/// total `total_length` over their number `num_seqs` as NCBI's `Int4`. A `num_seqs` of 0 (NCBI
/// divides by 0: LOSAT refuses it before, `common::database::check_setup_count`) gives `None`
/// here too, so that no input reaches a division by 0.
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
            let average = total_length.checked_div(i64::from(num_seqs))?;
            average as i32 as u32
        }
    };
    (min_subject_length != 0).then_some(min_subject_length as i32)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:755-760
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:242-248
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
#[path = "setup_tests.rs"]
mod tests;
