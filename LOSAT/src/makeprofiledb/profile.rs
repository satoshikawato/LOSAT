//! One profile of makeprofiledb: the checks of `x_CheckInputScoremat`, the parameters each
//! profile shares with its volume (`x_UpdateRPSDbInfo`, `x_FillInRPSDbParameters`), and the
//! scores of a profile of frequency ratios only, made by NCBI's PSSM engine
//! (`x_RPSUpdateStatistics`).

use crate::algorithm::pssm_search::in_pssm::{engine_score_block, query_ungapped_block};
use crate::config::ScoringMatrix;
use crate::core::blast_util::BLASTAA_SIZE;
use crate::report::serial::pssm::{
    FormatRpsDbParameters, PssmFinalData, PssmWithParameters, SeqData,
};
use crate::writedb::convert;

use super::rps_files::{freq_ratio_matrix, RpsError, VolumeParams};

/// What `x_CheckInputScoremat` found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Scoremat {
    /// `sm_valid_has_pssm`: the profile has its scores.
    HasPssm,
    /// `sm_valid_freq_only`: the profile has frequency ratios only.
    FreqOnly,
}

/// Why a profile stops the build, as makeprofiledb's handler reports it.
#[derive(Debug)]
pub(super) enum ProfileError {
    /// `CInputException`: "INPUT ERROR: " and the message, exit status 1.
    Input(Vec<u8>),
    /// The RPS files' errors (`CBlastException`, ...).
    Rps(RpsError),
    /// `CBlastException` with a message of the engine (`CPssmEngineException` too, which derives
    /// from it): exit status 1.
    Blast(String),
    /// Another `CException` (`CSeqUtilException` of the query's conversion): exit status 255.
    Exception(crate::ncbi_exception::NcbiException),
}

impl From<RpsError> for ProfileError {
    fn from(error: RpsError) -> Self {
        ProfileError::Rps(error)
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:601-668
// NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:62-74
// `GetQueryLength` is 0 for a query that is a set, or without a length. The length is compared
// with `numColumns` as an `int` ([`query_length`]).
/// `x_CheckInputScoremat` of the profile of the file `filename` (a `-dbtype rps` build).
pub(super) fn check_input_scoremat(
    pssm_w_parameters: &PssmWithParameters,
    filename: &[u8],
) -> Result<Scoremat, ProfileError> {
    let input = |what: &[u8]| ProfileError::Input([filename, what].concat());
    let pssm = &pssm_w_parameters.pssm;
    let length = query_length(pssm_w_parameters);
    if length == 0 {
        return Err(input(b" contains no bioseq data"));
    }
    if length != pssm.num_columns {
        return Err(input(b" 's num of columns does not match size of sequence"));
    }
    if pssm.num_rows <= 0 || pssm.num_rows > BLASTAA_SIZE as i32 {
        return Err(input(b" has invalid alphabet size"));
    }
    if pssm.intermediate_data.is_none() {
        return Err(input(
            b" contains no frequency ratios.\nPlease use a recent version of psiblast to \
              regenerate PSSM files\n",
        ));
    }
    Ok(if pssm.final_data.is_some() {
        Scoremat::HasPssm
    } else {
        Scoremat::FreqOnly
    })
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:624,1406
// The `size_t` of a `TSeqPos` (`Uint4`) converted to `int`, which GCC defines as the value modulo
// 2^32: a length past `INT_MAX` is negative.
/// `GetQueryLength()` cast to `int`: the length of the query when it is a Bioseq with one, else 0,
/// as an `int`.
pub(super) fn query_length(pssm_w_parameters: &PssmWithParameters) -> i32 {
    pssm_w_parameters
        .pssm
        .query
        .as_ref()
        .and_then(|query| query.seq())
        .and_then(|seq| seq.inst.length)
        .unwrap_or(0) as i32
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:786-803
/// `x_UpdateRPSDbInfo`: the first profile's `rpsdbparams` replace the volume's gap costs and
/// matrix.
pub(super) fn update_rps_db_info(
    params: &mut VolumeParams,
    pssm_w_parameters: &PssmWithParameters,
) {
    let Some(rps) = rpsdbparams(pssm_w_parameters) else {
        return;
    };
    if let Some(gap_extend) = rps.gap_extend {
        params.gap_extend = gap_extend;
    }
    if let Some(gap_open) = rps.gap_open {
        params.gap_open = gap_open;
    }
    params.matrix.clone_from(&rps.matrix_name);
}

fn rpsdbparams(pssm_w_parameters: &PssmWithParameters) -> Option<&FormatRpsDbParameters> {
    pssm_w_parameters
        .params
        .as_ref()
        .and_then(|params| params.rpsdbparams.as_ref())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:805-829
/// `x_FillInRPSDbParameters`: the profile's gap costs and matrix agree with the volume's
/// (the ones it does not set are the volume's).
pub(super) fn fill_in_rps_db_parameters(
    params: &VolumeParams,
    pssm_w_parameters: &PssmWithParameters,
) -> Result<(), RpsError> {
    let Some(rps) = rpsdbparams(pssm_w_parameters) else {
        return Ok(());
    };
    if rps
        .gap_extend
        .is_some_and(|gap_extend| gap_extend != params.gap_extend)
    {
        return Err(RpsError::Blast("Gap extend penalties do not match"));
    }
    if rps
        .gap_open
        .is_some_and(|gap_open| gap_open != params.gap_open)
    {
        return Err(RpsError::Blast("Gap open penalties do not match"));
    }
    if rps.matrix_name != params.matrix {
        return Err(RpsError::Blast("Score matrix does not match"));
    }
    Ok(())
}

fn query_seq_data(pssm_w_parameters: &PssmWithParameters) -> Option<&SeqData> {
    pssm_w_parameters
        .pssm
        .query
        .as_ref()
        .and_then(|query| query.seq())
        .and_then(|seq| seq.inst.seq_data.as_ref())
}

// NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:76-99
// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:2009-2066
// `CSeqportUtil::Convert` hands NCBIeaa and IUPACaa to `CSeqConvert`, whose tables are the ones
// the database writer converts with (`writedb::convert`): every byte has a code, 255 for one
// that is no residue (NCBIeaa: other than a capital letter, `*` or `-`; IUPACaa: other than a
// letter of either case). NCBIstdaa is copied.
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.hpp:50-70
// `CSeqConvert` has no conversion from a nucleotide coding to NCBIstdaa: data of IUPACna,
// NCBI2na or NCBI4na is its exception (NCBI BLAST+ 2.17.0: `Unknown conversion.`, exit status
// 255, row PS10d survey `FA_*`), but empty data returns before the conversion, an empty query.
/// `GetQuerySequenceData`: the query's NCBIstdaa codes (none without data).
fn query_ncbistdaa(pssm_w_parameters: &PssmWithParameters) -> Result<Vec<u8>, ProfileError> {
    match query_seq_data(pssm_w_parameters) {
        None => Ok(Vec::new()),
        Some(SeqData::Ncbistdaa(codes)) => Ok(codes.clone()),
        Some(SeqData::Ncbieaa(data)) => Ok(convert::ncbieaa_codes(data)),
        Some(SeqData::Iupacaa(data)) => Ok(convert::iupacaa_codes(data)),
        Some(SeqData::Iupacna(data) | SeqData::Ncbi2na(data) | SeqData::Ncbi4na(data)) => {
            if data.is_empty() {
                Ok(Vec::new())
            } else {
                Err(ProfileError::Exception(convert::unknown_conversion()))
            }
        }
        Some(SeqData::Other(_)) => Err(ProfileError::Rps(RpsError::Refused(
            "a profile of frequency ratios whose query has data in another coding (NCBI's \
             sequence conversion of it is not followed)"
                .to_string(),
        ))),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:498-515
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1962-1971,1984-2007
// The engine's score block reads the matrix first (a name it does not know is its error, the
// engine's `Unknown error code -1`), then counts the query's composition in an array of the 28
// codes: a code past them (255 for a byte that is no residue) is counted outside the array, the
// first access of NCBI's run outside its data (valgrind's first invalid access). What NCBI BLAST+
// 2.17.0 does after it depends on the code (row PS10c probes `PR/ps10c/probe*`): for 28, 29, 64,
// 128 and 255 it ends with SIGSEGV where the PSSM conversion reads the matrix's row of the code,
// for 31 malloc aborts on its corrupted heap (also where a negative frequency ratio is otherwise
// reported), for 32 it builds the database. Before the count, the lengths, the frequency ratios
// and the matrix are NCBI's.
/// The refusal of a query with a code past the 28 residues of NCBIstdaa, for a matrix that
/// NCBI's engine reads; `None` otherwise.
fn query_code_past_residues(
    pssm_w_parameters: &PssmWithParameters,
    query: &[u8],
    matrix_name: &str,
) -> Option<ProfileError> {
    matrix_name.parse::<ScoringMatrix>().ok()?;
    let position = query
        .iter()
        .position(|&code| usize::from(code) >= BLASTAA_SIZE)?;
    let what = match query_seq_data(pssm_w_parameters) {
        Some(SeqData::Ncbieaa(data) | SeqData::Iupacaa(data)) => {
            format!("the letter 0x{:02x}", data[position])
        }
        _ => format!("the NCBIstdaa code {}", query[position]),
    };
    Some(ProfileError::Rps(RpsError::Refused(format!(
        "a profile of frequency ratios whose query has {what} at position {}, no residue of \
         NCBIstdaa (NCBI's PSSM engine writes outside its table of residues as it counts the \
         query's composition)",
        position + 1
    ))))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:841-882
// NCBI reference (598d8ae6): c++/include/algo/blast/api/psi_pssm_input.hpp:78-128
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:174-186,336-361,512-583
// The engine takes the gap costs 0 as 11 and 1 (`IPssmInputFreqRatios`'s defaults) and the
// IMPALA scale factor 0 as 1 (`kPSSM_NoImpalaScaling`). Its score block's errors come first
// (`x_InitializeScoreBlock` in the constructor), then the negative ratios (`s_Validate`), then
// the status of `PSICreatePssmFromFrequencyRatios`, which reads the matrix's frequency ratios
// (a null table for IDENTITY) and the query's ungapped Karlin-Altschul block: one that cannot be
// computed leaves NCBI's `kbp_psi[0]` null, which the scaling reads (a crash).
/// `x_RPSUpdateStatistics`: the profile's `finalData` (the scores, lambda, kappa and h of the
/// engine, the scale factor of the volume) made from its frequency ratios.
pub(super) fn update_statistics(
    params: &VolumeParams,
    pssm_w_parameters: &mut PssmWithParameters,
    seq_size: i32,
) -> Result<(), ProfileError> {
    let query = query_ncbistdaa(pssm_w_parameters)?;
    // The size of `query_v` cast to `Int4`, compared with `seq_size`: a negative length (one past
    // `INT_MAX`) never matches a query of fewer than 2^31 letters.
    if query.len() as i32 != seq_size {
        return Err(ProfileError::Blast(
            "Query sequence lengths mismatch".to_string(),
        ));
    }
    let seq_size = usize::try_from(seq_size).map_err(|_| {
        ProfileError::Rps(RpsError::Refused(
            "a profile of frequency ratios whose query has 2^31 letters or more".to_string(),
        ))
    })?;
    let pssm = &pssm_w_parameters.pssm;
    let ratios = pssm
        .intermediate_data
        .as_ref()
        .map_or(&[][..], |data| &data.freq_ratios[..]);
    let matrix = freq_ratio_matrix(
        ratios,
        pssm.num_rows as usize,
        seq_size,
        pssm.by_row.unwrap_or(false),
    )?;
    let gap_open = if params.gap_open != 0 {
        params.gap_open
    } else {
        11
    };
    let gap_extend = if params.gap_extend != 0 {
        params.gap_extend
    } else {
        1
    };
    let matrix_name = String::from_utf8_lossy(&params.matrix).into_owned();
    if let Some(refusal) = query_code_past_residues(pssm_w_parameters, &query, &matrix_name) {
        return Err(refusal);
    }
    let refuse = |what: String| ProfileError::Rps(RpsError::Refused(what));
    let (scoring, sbp) = match engine_score_block(&matrix_name, gap_open, gap_extend, &query)
        .map_err(|error| refuse(format!("the scoring system of a profile ({error:#})")))?
    {
        Ok(block) => block,
        Err(message) => return Err(ProfileError::Blast(message)),
    };
    if matrix.iter().any(|&ratio| ratio < 0.0) {
        return Err(ProfileError::Blast(
            "PSSM frequency ratios cannot have negative values".to_string(),
        ));
    }
    let Some(sbp) = sbp else {
        return Err(refuse(format!(
            "a profile of frequency ratios for the matrix {matrix_name}, which has no frequency \
             ratios of its own (NCBI's PSSM engine reads a null table there)"
        )));
    };
    let Some(kbp_psi) = query_ungapped_block(scoring, &query) else {
        return Err(refuse(
            "a profile of frequency ratios whose query has no ungapped Karlin-Altschul \
             parameters (NCBI's PSSM engine reads a null block there)"
                .to_string(),
        ));
    };
    let rows: Vec<[f64; BLASTAA_SIZE]> = (0..seq_size)
        .map(|c| std::array::from_fn(|r| matrix[r * seq_size + c]))
        .collect();
    let impala = if params.scale_factor != 0 {
        f64::from(params.scale_factor)
    } else {
        1.0
    };
    let (scores, karlin) =
        crate::core::psi::pssm_from_freq_ratios_engine(&query, &rows, &sbp, kbp_psi, impala)
            .map_err(|error| match error {
                crate::core::psi::PsiError::ScalingUndefined => refuse(
                    "a profile of frequency ratios whose scaling leaves C's defined behaviour in \
                     NCBI's PSSM engine (a score past int, or a scale that grows without end)"
                        .to_string(),
                ),
                error => ProfileError::Blast(
                    crate::algorithm::psiblast::pssm::error_code_to_string(error).to_string(),
                ),
            })?;
    let final_data = PssmFinalData {
        scores: scores.iter().flatten().copied().collect(),
        lambda: karlin.lambda,
        kappa: karlin.k,
        h: karlin.h,
        scaling_factor: Some(params.scale_factor),
        lambda_ungapped: None,
        kappa_ungapped: None,
        h_ungapped: None,
        word_score_threshold: None,
    };
    pssm_w_parameters.pssm.final_data = Some(final_data);
    Ok(())
}

#[cfg(test)]
#[path = "profile_tests.rs"]
mod tests;
