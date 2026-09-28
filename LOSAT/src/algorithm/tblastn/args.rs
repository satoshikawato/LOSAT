//! Pinned NCBI TBLASTN local-subject search and Stage E report options.
use anyhow::{bail, Context, Result};
use bio::io::fasta;
use clap::Args;
use std::io::Write;
use std::path::PathBuf;

use crate::blastinput::value_parsers::*;
use crate::utils::genetic_code::GeneticCode;

use super::scoring::{matrix_params, suggested_threshold, suggested_window_size};
use super::stage_d_pipeline::{
    run_local_for_report_threads, LocalStageDProfile, LocalStageDScoring,
};
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::render;
use crate::config::ScoringMatrix;

// NCBI reference: c++/src/algo/blast/blastinput/tblastn_args.cpp:45-125
// ```c++
// static const char kDefaultTask[] = "tblastn";
// tasks.insert(kDefaultTask); tasks.insert("tblastn-fast");
// arg.Reset(new CGenericSearchArgs(kQueryIsProtein));
// arg.Reset(new CGeneticCodeArgs(CGeneticCodeArgs::eDatabase));
// arg.Reset(new CCompositionBasedStatsArgs);
// m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eNucleotideDb));
// ```
// NCBI reference: c++/src/algo/blast/api/tblastn_options.cpp:54-85
// ```c++
// SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
// m_Opts->SetSumStatisticsMode();
// m_Opts->SetCompositionBasedStats(eCompositionMatrixAdjust);
// SetDbGeneticCode(BLAST_GENETIC_CODE);
// ```
#[derive(Args, Debug)]
#[command(rename_all = "snake_case")]
pub struct TblastnArgs {
    // NCBI blast_args.cpp:3425-3427:
    // arg_desc.AddDefaultKey(kArgQuery, "input_file", "Input file name",
    //                        CArgDescriptions::eInputFile, kDfltArgQuery);
    // The NCBI stdin default is an explicitly unsupported LOSAT input path.
    #[arg(long, value_parser = file_path(), value_name = "PATH")]
    pub query: Option<PathBuf>,
    #[arg(long, value_parser = file_path(), value_name = "PATH", conflicts_with = "db")]
    pub subject: Option<PathBuf>,
    // NCBI blast_args.cpp:2377-2380:
    // arg_desc.SetDependency(kArgSubjectLocation,
    //                        CArgDescriptions::eExcludes, *dbarg);
    #[arg(long, value_name = "DB", conflicts_with_all = ["subject", "subject_loc"])]
    pub db: Option<String>,
    #[arg(long, default_value = "tblastn", value_parser = ["tblastn", "tblastn-fast"])]
    pub task: String,
    #[arg(long, value_name = "PATH")]
    pub out: Option<PathBuf>,
    #[arg(long, default_value_t = 10.0, value_parser = tblastn_evalue)]
    pub evalue: f64,
    #[arg(long, default_value_t = 3, value_parser = tblastn_word_size)]
    pub word_size: usize,
    // NCBI blast_args.cpp:258-273: when omitted, each cost keeps the
    // selected matrix's BLAST_MATRIX_BEST value.
    // if (args.Exist(kArgGapOpen) && args[kArgGapOpen]) opt.SetGapOpeningCost(...);
    // else if (args.Exist(kArgMatrixName) && args[kArgMatrixName]) opt.SetGapOpeningCost(gap_open);
    #[arg(long = "gapopen", value_parser = nonnegative_i32)]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = nonnegative_i32)]
    pub gap_extend: Option<i32>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1029-1056
    // ```c++
    // arg_desc.AddDefaultKey(kArgDbGeneticCode, "int_value",
    //                        "Genetic code to use to translate "
    //                        "database/subjects (see user manual for details)\n",
    //                        CArgDescriptions::eInteger,
    //                        NStr::IntToString(BLAST_GENETIC_CODE));
    // opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
    // ```
    // PD-TLOSAN-LOCAL-GENCODE-32 adds gc.prt ID 32 only for TBLASTN.
    #[arg(long, default_value_t = 1, value_parser = tblastn_genetic_code)]
    pub db_gencode: u8,
    #[arg(long, default_value_t = 0, value_parser = nonnegative_usize)]
    pub max_intron_length: usize,
    #[arg(long, default_value = "BLOSUM62")]
    pub matrix: String,
    #[arg(long, value_parser = tblastn_threshold)]
    pub threshold: Option<f64>,
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:220-229,280-286:
    // arg_desc.AddOptionalKey(kArgGappedXDropoff, "float_value", ...);
    // arg_desc.AddOptionalKey(kArgFinalGappedXDropoff, "float_value", ...);
    // opt.SetGapXDropoff(args[kArgGappedXDropoff].AsDouble());
    // opt.SetGapXDropoffFinal(args[kArgFinalGappedXDropoff].AsDouble());
    #[arg(long, value_parser = nonnegative_f64)]
    pub xdrop_gap: Option<f64>,
    #[arg(long, value_parser = nonnegative_f64)]
    pub xdrop_gap_final: Option<f64>,
    #[arg(long, default_value = "2", value_parser = tblastn_composition)]
    pub comp_based_stats: String,
    #[arg(long, default_value = "0", value_parser = tblastn_outfmt)]
    pub outfmt: String,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:332-342
    // ```c++
    // arg_desc.AddDefaultKey(kArgSegFiltering, "SEG_options",
    //                        "Filter query sequence with SEG "
    //                        "(Format: '" + kDfltArgApplyFiltering + "', " +
    //                        "'window locut hicut', or '" + kDfltArgNoFiltering +
    //                        "' to disable)",
    //                        CArgDescriptions::eString, m_FilterByDefault
    //                        ? kDfltArgSegFiltering : kDfltArgNoFiltering);
    // arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                        "Apply filtering locations as soft masks",
    //                        CArgDescriptions::eBoolean,
    //                        kDfltArgLookupTableMaskingOnlyProt);
    // ```
    #[arg(long, default_value = "12 2.2 2.5", value_parser = parse_seg_filtering)]
    pub seg: SegSpec,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942,2549-2557
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //     "Use lower case filtering in query and subject sequence(s)?", true);
    // ReadSequencesToBlast(..., use_lcase_masks, subjects, ...);
    #[arg(long)]
    pub lcase_masking: bool,
    #[arg(long, default_value = "false", value_parser = ncbi_bool, action = clap::ArgAction::Set, num_args = 1)]
    pub soft_masking: bool,
    #[arg(long, default_value = "true", value_parser = ncbi_bool, action = clap::ArgAction::Set, num_args = 1)]
    pub sum_stats: bool,
    #[arg(long, value_parser = nonnegative_usize)]
    pub window_size: Option<usize>,
    #[arg(long, default_value_t = 500, value_parser = positive_usize)]
    pub max_target_seqs: usize,
    #[arg(long, default_value_t = 1, value_parser = positive_usize)]
    pub num_threads: usize,
    #[arg(long)]
    pub ungapped: bool,
    #[arg(long, value_name = "PATH", conflicts_with = "query")]
    pub in_pssm: Option<PathBuf>,
    // NCBI blast_args.cpp:3164-3166,1139-1142,2384-2385:
    // arg_desc.SetDependency(kArgNumThreads, CArgDescriptions::eExcludes, kArgRemote);
    // arg_desc.SetDependency(kArgPSIInputChkPntFile, CArgDescriptions::eExcludes, kArgRemote);
    // arg_desc.SetDependency(kArgSubjectLocation, CArgDescriptions::eExcludes, kArgRemote);
    #[arg(long, conflicts_with_all = ["subject_loc", "num_threads", "in_pssm"])]
    pub remote: bool,
    #[arg(long, conflicts_with = "remote")]
    pub subject_loc: Option<String>,
}

// NCBI reference: c++/src/algo/blast/core/blast_options.c:1518-1523
// ```c
// if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
// {
//     Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
//         "expect value or cutoff score must be greater than zero");
//     return BLASTERR_OPTION_VALUE_INVALID;
// }
// ```
// No cutoff-score option is exposed on this local TBLASTN profile.
fn tblastn_evalue(value: &str) -> std::result::Result<f64, String> {
    let evalue = nonnegative_f64(value)?;
    if evalue == 0.0 {
        return Err("expect value or cutoff score must be greater than zero".into());
    }
    Ok(evalue)
}

// NCBI reference: c++/src/algo/blast/core/blast_options.c:1335-1348
// ```c
// if (program_number == eBlastTypeBlastp ||
//     program_number == eBlastTypeTblastn ||
//     program_number == eBlastTypeBlastx)
// {
//     if (options->word_size > 7) {
//         Blast_MessageWrite(blast_msg, eBlastSevError,
//                            kBlastMessageNoContext,
//                            "Word-size must be less than "
//                            "8 for a tblastn, blastp or blastx search");
//         return BLASTERR_OPTION_VALUE_INVALID;
//     }
// }
// ```
fn tblastn_word_size(value: &str) -> std::result::Result<usize, String> {
    let word_size = blastp_word_size(value)?;
    if word_size > 7 {
        return Err("Word-size must be less than 8 for a tblastn, blastp or blastx search".into());
    }
    Ok(word_size)
}

// NCBI reference: c++/src/algo/blast/core/blast_options.c:1303-1310
// ```c
// if (program_number != eBlastTypeBlastn &&
//     program_number != eBlastTypeMapping &&
//     (!Blast_ProgramIsRpsBlast(program_number)) &&
//     options->threshold <= 0)
// {
//     Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
//                        "Non-zero threshold required");
//     return BLASTERR_OPTION_VALUE_INVALID;
// }
// ```
fn tblastn_threshold(value: &str) -> std::result::Result<f64, String> {
    let threshold = nonnegative_f64(value)?;
    if threshold == 0.0 {
        return Err("Non-zero threshold required".into());
    }
    Ok(threshold)
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:997-1056
// ```c++
// static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,
//                     21,22,23,24,25,26,27,28,29,30,31,33};
// return (genetic_codes.find(val) != genetic_codes.end());
// ```
// Stage A product decision admits 32 from gc.prt:340-347 for TBLASTN only.
pub fn tblastn_genetic_code(value: &str) -> std::result::Result<u8, String> {
    let id = value.parse::<u8>().map_err(|_| "invalid genetic code ID")?;
    GeneticCode::try_from_id(id).map(|_| id)
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:838-866
// ```c++
// ECompoAdjustModes compo_mode = eNoCompositionBasedStats;
// switch (comp_stat_string[0]) {
//     case '0': case 'F': case 'f':
//         compo_mode = eNoCompositionBasedStats;
//         break;
//     case '2':
//         compo_mode = eCompositionMatrixAdjust;
//         break;
// }
// ```
pub fn tblastn_composition(value: &str) -> std::result::Result<String, String> {
    // The NCBI switch has no default error: an unknown initial character
    // leaves eNoCompositionBasedStats selected. The pinned CLI accepts it.
    Ok(value.to_string())
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
// ```c++
// arg_desc.AddDefaultKey(kArgOutputFormat, "format",
//                        kOutputFormatDescription,
//                        CArgDescriptions::eString,
//                        NStr::IntToString(dft_outfmt));
// ```
pub fn tblastn_outfmt(value: &str) -> std::result::Result<String, String> {
    if matches!(value, "0" | "6" | "7") {
        Ok(value.to_string())
    } else {
        Err("unsupported TBLASTN outfmt: available formats are 0, 6 and 7".into())
    }
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:339-342,389-390
// ```c++
// arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
//                        "Apply filtering locations as soft masks",
//                        CArgDescriptions::eBoolean,
//                        kDfltArgLookupTableMaskingOnlyProt);
// opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
// ```
fn ncbi_bool(value: &str) -> std::result::Result<bool, String> {
    match value.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" | "1" => Ok(true),
        "false" | "f" | "no" | "n" | "0" => Ok(false),
        _ => Err(format!("invalid Boolean '{value}'")),
    }
}

impl TblastnArgs {
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:289-304,599-618
    // ```c++
    // if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
    //    opt.SetWordThreshold(19.3);
    //    if (args[kArgWordSize].AsInteger() > 5) opt.SetWordThreshold(21.0);
    //    if (args[kArgWordSize].AsInteger() > 6) opt.SetWordThreshold(20.25);
    // }
    // if (args[kArgWordScoreThreshold]) opt.SetWordThreshold(args[kArgWordScoreThreshold].AsDouble());
    // else if (s_IsDefaultWordThreshold(...)) BLAST_GetSuggestedThreshold(...);
    // ```
    pub fn effective_threshold(&self) -> f64 {
        self.threshold.unwrap_or_else(|| {
            if self.word_size == 5 {
                19.3
            } else if self.word_size == 6 {
                21.0
            } else if self.word_size == 7 {
                20.25
            } else {
                suggested_threshold(&self.matrix)
            }
        })
    }

    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:486-497
    // ```c++
    // if (args[kArgWindowSize]) opt.SetWindowSize(args[kArgWindowSize].AsInteger());
    // else BLAST_GetSuggestedWindowSize(opt.GetProgramType(), opt.GetMatrixName(), &window);
    // ```
    pub fn effective_window_size(&self) -> usize {
        self.window_size
            .unwrap_or_else(|| suggested_window_size(&self.matrix))
    }

    // NCBI reference: c++/src/algo/blast/blastinput/tblastn_args.cpp:55-62,134-152
    // ```c++
    // tasks.insert("tblastn-fast");
    // if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile])
    //     rv->SetPSITblastnDefaults();
    // else return x_CreateOptionsHandleWithTask(locality, args[kTask].AsString());
    // ```
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2364-2385,872-876
    // ```c++
    // arg_desc.SetDependency(kArgSubject, CArgDescriptions::eExcludes, *dbarg);
    // if (ungapped && *ungapped && compo_mode != eNoCompositionBasedStats)
    //     NCBI_THROW(CInputException, eInvalidInput, "Composition-adjusted searched are not supported ...");
    // ```
    pub fn validate(&self) -> Result<()> {
        // NCBI reference: c++/src/algo/blast/core/blast_options.c:913-943,1778-1796
        // ```c
        // Blast_KarlinBlkGappedLoadFromTables(NULL, options->gap_open,
        //     options->gap_extend, options->matrix, std_matrix_only);
        // if (lookup_options->word_size > 5 && is_identity) { return BLASTERR_OPTION_VALUE_INVALID; }
        // ```
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:258-273
        // ```c++
        // BLAST_GetProteinGapExistenceExtendParams(matrix, &gap_open, &gap_extend);
        // if (args[kArgGapOpen]) opt.SetGapOpeningCost(args[kArgGapOpen].AsInteger());
        // else if (args[kArgMatrixName]) opt.SetGapOpeningCost(gap_open);
        // ```
        let params = matrix_params(&self.matrix)
            .ok_or_else(|| anyhow::anyhow!("{} is not a supported matrix", self.matrix))?;
        let (preferred_open, preferred_extend) = params.preferred;
        let gap_pair = (
            self.gap_open.unwrap_or(preferred_open),
            self.gap_extend.unwrap_or(preferred_extend),
        );
        if !self.ungapped && !params.allowed.contains(&gap_pair) {
            bail!(
                "Gap existence and extension values of {} and {} not supported for {}",
                gap_pair.0,
                gap_pair.1,
                self.matrix
            );
        }
        if self.matrix.eq_ignore_ascii_case("IDENTITY") && self.word_size > 5 {
            bail!("Word size larger than 5 is not supported for the identity scoring matrix");
        }
        if self.task != "tblastn" {
            bail!("unsupported TBLASTN task: {}", self.task);
        }
        if self.in_pssm.is_some() {
            bail!("unsupported PSI-TBLASTN checkpoint search");
        }
        if self.db.is_some() {
            bail!("unsupported TBLASTN database search");
        }
        if self.remote {
            bail!("unsupported TBLASTN remote search");
        }
        if self.subject.is_none() {
            bail!("TBLASTN requires -subject for the local path");
        }
        if self.query.is_none() {
            bail!("TBLASTN requires -query file input; stdin is unimplemented");
        }
        if self.subject_loc.is_some() {
            bail!("unsupported TBLASTN -subject_loc until local range handling is implemented");
        }
        if self.ungapped
            && matches!(
                self.comp_based_stats.chars().next(),
                Some('1' | 'D' | 'd' | '2' | '3' | 'T' | 't')
            )
        {
            bail!("Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a gapped search");
        }
        Ok(())
    }

    // NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
    // results = lcl_blast.Run();
    // formatter.PrintOneResultSet(**result, query);
    // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1458:
    // PrintOneResultSet dispatches the complete result to outfmt 0/6/7.
    pub fn run(self) -> Result<()> {
        self.validate()?;
        // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:838-866;
        // core/blast_traceback.c:1481-1501:
        // compo_mode selects ordinary traceback or composition redo.
        // The Stage D port covers only these two choices.
        let composition_mode2 = match self.comp_based_stats.chars().next() {
            Some('0' | 'F' | 'f') => false,
            Some('2') => true,
            _ => bail!(
                "unsupported TBLASTN -comp_based_stats value: {}",
                self.comp_based_stats
            ),
        };
        // NCBI c++/src/algo/blast/core/blast_options.c:903-936;
        // core/blast_parameters.c:774-815:
        // Unsupported scoring, intron, and ungapped branches fail before
        // opening the output stream.
        if self.max_intron_length != 0 {
            bail!("unsupported TBLASTN -max_intron_length");
        }
        if self.ungapped {
            bail!("unsupported TBLASTN -ungapped search");
        }
        // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187:
        // arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
        // NCBI c++/src/algo/blast/api/prelim_stage.cpp:145-147:
        // TBlastThreads the_threads(GetNumberOfThreads());
        crate::utils::threading::validate_threads(self.num_threads)?;
        let matrix: ScoringMatrix = self.matrix.parse().map_err(anyhow::Error::msg)?;
        let params = matrix_params(&self.matrix).context("missing TBLASTN matrix parameters")?;
        let (gap_open, gap_extend) = (
            self.gap_open.unwrap_or(params.preferred.0),
            self.gap_extend.unwrap_or(params.preferred.1),
        );
        let threshold = self.effective_threshold();
        let window = self.effective_window_size();
        let profile_supported = (matrix == ScoringMatrix::Blosum62
            && self.word_size == 3
            && gap_open == 11
            && gap_extend == 1
            && threshold == 13.0
            && window == 40)
            || (matrix == ScoringMatrix::Blosum45
                && !composition_mode2
                && self.word_size == 2
                && gap_open == 14
                && gap_extend == 2
                && threshold == 16.0
                && window == 60);
        if !profile_supported {
            bail!("unsupported TBLASTN scoring and lookup option combination");
        }
        let query_path = self.query.as_ref().context("TBLASTN query file missing")?;
        let subject_path = self
            .subject
            .as_ref()
            .context("TBLASTN subject file missing")?;
        // NCBI c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246:
        // sequences = input->GetAllSeqs(*scope);
        let read = |path: &std::path::Path| -> Result<Vec<fasta::Record>> {
            fasta::Reader::from_file(path)
                .with_context(|| format!("failed to open FASTA {}", path.display()))?
                .records()
                .collect::<std::result::Result<Vec<_>, _>>()
                .with_context(|| format!("failed to read FASTA {}", path.display()))
        };
        let queries = read(query_path)?;
        let subjects = read(subject_path)?;
        let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
        let subject_seqs: Vec<_> = subjects
            .iter()
            .map(|record| record.seq().to_vec())
            .collect();
        let seg = self.seg.params();
        // NCBI c++/src/algo/blast/api/prelim_stage.cpp:172-188:
        // (*thread)->Run(); (*thread)->Join(&result);
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1417:
        // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
        //                         CConstRef<blast::CBlastQueryVector> queries,
        // NCBI c++/src/app/blast/tblastn_app.cpp:251-252,275-301:
        // input.Reset(new CBlastInput(&*fasta, GetQueryBatchSize()));
        // for (; !input->End(); ...) {
        //     query = input->GetNextSeqBatch(*scope);
        //     CLocalBlast lcl_blast(query_factory, ..., db_adapter);
        //     results = lcl_blast.Run();
        // }
        // NCBI c++/src/algo/blast/blastinput/blast_input_aux.cpp:70-127:
        // case eTblastn: retval = 20000;
        // Search each accumulated batch independently; preliminary linking
        // and sum-statistics depend on its concatenated query contexts.
        let mut results = Vec::with_capacity(queries.len());
        let mut lengths: Option<LocalSubjectParameters> = None;
        let mut ungapped_karlin = Vec::with_capacity(queries.len());
        let mut query_validity = Vec::with_capacity(queries.len());
        let mut query_batch_skipped = Vec::with_capacity(queries.len());
        for range in tblastn_query_batches(&query_seqs) {
            let (mut batch_results, batch_lengths, batch_karlin, batch_validity) =
                run_local_for_report_threads(
                    &query_seqs[range],
                    &subject_seqs,
                    LocalStageDProfile {
                        seg: seg.as_ref(),
                        soft_masking: self.soft_masking,
                        mask_lowercase: self.lcase_masking,
                        genetic_code: self.db_gencode,
                        expect_value: self.evalue,
                        max_target_seqs: self.max_target_seqs,
                    },
                    composition_mode2,
                    self.sum_stats,
                    LocalStageDScoring {
                        matrix,
                        gap_open,
                        gap_extend,
                        word_size: self.word_size,
                        threshold: threshold as i32,
                        window: i32::try_from(window)?,
                        gap_xdrop_bits: self.xdrop_gap.unwrap_or(15.0),
                        final_xdrop_bits: self.xdrop_gap_final.unwrap_or(25.0),
                    },
                    self.num_threads,
                )?;
            results.append(&mut batch_results);
            if let Some(all_lengths) = &mut lengths {
                // NCBI tblastn_app.cpp:288-301 formats each batch's own
                // query-context lengths; retain those values in input order.
                all_lengths.lengths.extend(batch_lengths.lengths);
                all_lengths.cutoffs.extend(batch_lengths.cutoffs);
            } else {
                lengths = Some(batch_lengths);
            }
            ungapped_karlin.extend(batch_karlin);
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-224:
            // if (CheckInternalData() != 0) each result in this Run() batch
            // receives an unsearched ancillary state and null align set.
            let skipped = batch_validity.iter().all(|&valid| !valid);
            query_batch_skipped.extend(std::iter::repeat(skipped).take(batch_validity.len()));
            query_validity.extend(batch_validity);
        }
        let lengths = lengths.context("TBLASTN query file is empty")?;
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1458:
        // formatter.PrintOneResultSet(...) writes only after a valid result.
        // Buffer the complete report so an error cannot leave partial success output.
        let mut bytes = Vec::new();
        render(
            &mut bytes,
            &self.outfmt,
            &queries,
            &subjects,
            subject_path,
            &mut results,
            &lengths,
            &ungapped_karlin,
            &query_validity,
            &query_batch_skipped,
            LocalStageDScoring {
                matrix,
                gap_open,
                gap_extend,
                word_size: self.word_size,
                threshold: threshold as i32,
                window: i32::try_from(window)?,
                gap_xdrop_bits: self.xdrop_gap.unwrap_or(15.0),
                final_xdrop_bits: self.xdrop_gap_final.unwrap_or(25.0),
            },
            self.db_gencode,
            seg.as_ref(),
            self.lcase_masking,
        )?;
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1443-1451:
        // if (results.HasWarnings()) ERR_POST(Warning << results.GetWarningStrings());
        // Print one setup warning per invalid query in formatter query order.
        for (index, (query, valid)) in queries.iter().zip(&query_validity).enumerate() {
            if !valid {
                std::io::stderr().write_all(&ncbi_invalid_query_warning(index, query))?;
            }
        }
        if let Some(path) = &self.out {
            std::fs::write(path, bytes)
                .with_context(|| format!("failed to write {}", path.display()))?;
        } else {
            std::io::stdout().write_all(&bytes)?;
        }
        Ok(())
    }
}

// NCBI c++/src/algo/blast/blastinput/blast_input_aux.cpp:105-127:
// case eTblastn: retval = 20000;
// NCBI c++/src/algo/blast/blastinput/blast_input.cpp:135-166:
// while (size_read < GetBatchSize()) { size_read += sequence::GetLength(...);
//                                 retval->AddQuery(q); }
// The query that reaches the threshold remains in the current batch.
fn tblastn_query_batches(queries: &[Vec<u8>]) -> Vec<std::ops::Range<usize>> {
    let mut batches = Vec::new();
    let mut start = 0;
    while start < queries.len() {
        let mut end = start;
        let mut residues = 0usize;
        while end < queries.len() && residues < 20_000 {
            residues += queries[end].len();
            end += 1;
        }
        batches.push(start..end);
        start = end;
    }
    batches
}

// NCBI c++/src/algo/blast/core/blast_stat.c:2780-2792:
// if (loop_status && !Blast_QueryIsTranslated(program))
//     Blast_MessageWrite(..., eBlastSevWarning, context,
//                        kBlastErrMsg_CantCalculateUngappedKAParams);
// NCBI c++/src/algo/blast/core/blast_message.c:37-40:
// kBlastErrMsg_CantCalculateUngappedKAParams = "Could not calculate ...".
// NCBI c++/src/algo/blast/api/blast_setup_cxx.cpp:535-543:
// query_id = id->GetSeqIdString() + " " + kTitle;
// if (query_id.size() > 35) query_id = query_id.substr(0, 25) + ".. ";
// NCBI c++/src/algo/blast/api/blast_results.cpp:277-293:
// retval = m_Errors.GetQueryId() + ": " + warning + " ";
fn ncbi_invalid_query_warning(index: usize, query: &fasta::Record) -> Vec<u8> {
    let mut query_id = format!("Query_{} {}", index + 1, query.id()).into_bytes();
    if let Some(desc) = query.desc() {
        query_id.extend_from_slice(b" ");
        query_id.extend_from_slice(desc.as_bytes());
    }
    if query_id.len() > 35 {
        query_id.truncate(25);
        query_id.extend_from_slice(b".. ");
    }
    let mut warning = b"Warning: [tblastn] ".to_vec();
    warning.extend_from_slice(&query_id);
    warning.extend_from_slice(b": Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n");
    warning
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI c++/src/algo/blast/blastinput/blast_input.cpp:137-165:
    // while (size_read < GetBatchSize()) { ... size_read += length;
    //                                  retval->AddQuery(q); }
    #[test]
    fn tblastn_batch_keeps_threshold_crossing_query() {
        let queries = vec![vec![b'A'; 19_999], vec![b'A'; 2], vec![b'A'; 1]];
        assert_eq!(tblastn_query_batches(&queries), vec![0..2, 2..3]);
    }
    use crate::cli::{try_parse_from, Cli, Commands};

    // NCBI c++/src/algo/blast/api/blast_setup_cxx.cpp:535-543;
    // c++/src/algo/blast/api/blast_results.cpp:277-293:
    // Query_1 + FASTA title is shortened after 35 bytes, then the
    // invalid-Karlin warning retains NCBI's trailing space and newline.
    #[test]
    fn invalid_query_warning_matches_ncbi_bytes() {
        let short = fasta::Record::with_attrs("nohit_query", None, b"W");
        assert_eq!(
            ncbi_invalid_query_warning(0, &short),
            b"Warning: [tblastn] Query_1 nohit_query: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options \n"
        );
        let long =
            fasta::Record::with_attrs("long_header", Some("abcdefghijklmnopqrstuvwxyz"), b"W");
        assert!(ncbi_invalid_query_warning(1, &long)
            .starts_with(b"Warning: [tblastn] Query_2 long_header abcde.. : "));
    }

    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //     "Use lower case filtering in query and subject sequence(s)?", true);
    #[test]
    fn lcase_masking_flag_is_explicit() {
        assert!(!parse(&[]).unwrap().lcase_masking);
        assert!(parse(&["-lcase_masking"]).unwrap().lcase_masking);
    }

    fn parse(extra: &[&str]) -> std::result::Result<TblastnArgs, clap::Error> {
        let mut argv = vec!["losat", "tblastn", "-query", "q.faa", "-subject", "s.fna"];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv)?;
        let Commands::Tblastn(args) = cli.command else {
            unreachable!()
        };
        Ok(args)
    }

    // NCBI tblastn_options.cpp:54-85 and blast_prot_options.cpp:85-141:
    // SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
    // m_Opts->SetSumStatisticsMode(); SetDbGeneticCode(BLAST_GENETIC_CODE);
    #[test]
    fn ordinary_defaults_and_code32_are_accepted() {
        let args = parse(&[]).unwrap();
        assert_eq!(args.task, "tblastn");
        assert_eq!(args.db_gencode, 1);
        assert_eq!(args.word_size, 3);
        assert_eq!(args.effective_threshold(), 13.0);
        assert_eq!(args.effective_window_size(), 40);
        assert_eq!(args.gap_open, None);
        assert_eq!(args.gap_extend, None);
        assert_eq!(args.evalue, 10.0);
        assert_eq!(args.comp_based_stats, "2");
        assert_eq!(args.matrix, "BLOSUM62");
        assert!(args.sum_stats);
        assert!(!args.soft_masking);
        assert_eq!(args.max_intron_length, 0);
        assert_eq!(args.outfmt, "0");
        args.validate().unwrap();
        assert!(parse(&["-soft_masking", "true", "-sum_stats", "false"]).is_ok());
        for id in [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            30, 31, 32, 33,
        ] {
            let id_string = id.to_string();
            assert_eq!(parse(&["-db_gencode", &id_string]).unwrap().db_gencode, id);
        }
    }

    // NCBI blast_args.cpp:997-1056 has 26 CLI values; gc.prt:340-347
    // defines ID 32 for the TBLASTN-only product decision.
    #[test]
    fn invalid_codes_are_explicit_errors() {
        for id in ["0", "7", "8", "17", "20", "34", "255", "256", "abc", "-1"] {
            let error = parse(&["-db_gencode", id]).unwrap_err().to_string();
            assert!(error.contains("genetic code"), "{id}: {error}");
        }
    }

    // NCBI tblastn_args.cpp:55-62,125-152 registers fast and PSI paths;
    // blast_args.cpp:2364-2385 excludes DB and subject combination.
    #[test]
    fn unsupported_paths_and_ncbi_exclusions_are_rejected() {
        assert!(parse(&["-task", "tblastn-fast"])
            .unwrap()
            .validate()
            .unwrap_err()
            .to_string()
            .contains("unsupported"));
        assert!(parse(&["-remote"])
            .unwrap()
            .validate()
            .unwrap_err()
            .to_string()
            .contains("unsupported"));
        let db_cli: Cli =
            try_parse_from(["losat", "tblastn", "-query", "q.faa", "-db", "db"]).unwrap();
        let Commands::Tblastn(db_args) = db_cli.command else {
            unreachable!()
        };
        assert!(db_args
            .validate()
            .unwrap_err()
            .to_string()
            .contains("database"));
        assert!(parse(&["-db", "db"]).is_err());
        let db_range: std::result::Result<Cli, clap::Error> = try_parse_from([
            "losat",
            "tblastn",
            "-query",
            "q.faa",
            "-db",
            "db",
            "-subject_loc",
            "1-30",
        ]);
        assert!(db_range.is_err());
        let psi_cli: Cli =
            try_parse_from(["losat", "tblastn", "-subject", "s.fna", "-in_pssm", "p.asn"]).unwrap();
        let Commands::Tblastn(psi_args) = psi_cli.command else {
            unreachable!()
        };
        assert!(psi_args
            .validate()
            .unwrap_err()
            .to_string()
            .contains("PSI-TBLASTN"));
        assert!(parse(&["-subject_loc", "1-30"])
            .unwrap()
            .validate()
            .unwrap_err()
            .to_string()
            .contains("subject_loc"));
        assert!(parse(&["-ungapped"])
            .unwrap()
            .validate()
            .unwrap_err()
            .to_string()
            .contains("Composition-adjusted"));
        assert!(parse(&["-ungapped", "-comp_based_stats", "F"])
            .unwrap()
            .validate()
            .is_ok());
        assert!(parse(&["-subject_loc", "1-30", "-remote"]).is_err());
        assert!(parse(&["-remote", "-num_threads", "2"]).is_err());
        let psi_remote: std::result::Result<Cli, clap::Error> = try_parse_from([
            "losat", "tblastn", "-subject", "s.fna", "-in_pssm", "p.asn", "-remote",
        ]);
        assert!(psi_remote.is_err());
        assert!(parse(&["-matrix", "BAD"])
            .unwrap()
            .validate()
            .unwrap_err()
            .to_string()
            .contains("not a supported matrix"));
        assert!(parse(&["-gapopen", "1", "-gapextend", "1"])
            .unwrap()
            .validate()
            .unwrap_err()
            .to_string()
            .contains("not supported for BLOSUM62"));
        assert!(parse(&["-matrix", "PAM30"]).unwrap().validate().is_ok());
        // NCBI blast_options.c:903-936 accepts PAM30, while the current
        // Rust Stage D search profile explicitly rejects its unported path.
        assert!(parse(&["-matrix", "PAM30"])
            .unwrap()
            .run()
            .unwrap_err()
            .to_string()
            .contains("unsupported"));
        assert!(
            parse(&["-matrix", "PAM30", "-gapopen", "9", "-gapextend", "1"])
                .unwrap()
                .validate()
                .is_ok()
        );
        assert!(parse(&["-matrix", "IDENTITY", "-word_size", "6"])
            .unwrap()
            .validate()
            .unwrap_err()
            .to_string()
            .contains("Word size larger than 5"));
        assert!(parse(&["-evalue", "0"])
            .unwrap_err()
            .to_string()
            .contains("expect value or cutoff score"));
        assert!(parse(&["-word_size", "1"]).is_err());
        assert!(parse(&["-word_size", "7"]).is_ok());
        assert!(parse(&["-word_size", "8"])
            .unwrap_err()
            .to_string()
            .contains("Word-size must be less than 8"));
        assert!(parse(&["-threshold", "0"])
            .unwrap_err()
            .to_string()
            .contains("Non-zero threshold required"));
        // NCBI blast_format.cpp:1411-1458 dispatches these three
        // implemented output paths from a complete local result set.
        for format in ["0", "6", "7"] {
            assert_eq!(parse(&["-outfmt", format]).unwrap().outfmt, format);
        }
        assert!(parse(&["-outfmt", "5"]).is_err());
        assert!(parse(&["-outfmt", "6 qseq sseq"]).is_err());
        assert!(parse(&["-comp_based_stats", "bogus"])
            .unwrap()
            .validate()
            .is_ok());
        assert!(parse(&["-comp_based_stats", "2u"]).is_ok());
        let unported = parse(&["-max_hsps", "2"]).unwrap_err().to_string();
        assert!(unported.contains("unsupported TBLASTN option"));
    }
}
