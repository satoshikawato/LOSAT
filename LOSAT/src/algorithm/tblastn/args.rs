//! Pinned NCBI TBLASTN option profile; no search engine is exposed yet.
use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

use crate::blastinput::value_parsers::*;
use crate::utils::genetic_code::GeneticCode;

use super::scoring::{matrix_params, suggested_threshold, suggested_window_size};

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
// No cutoff-score option is exposed on this Stage B TBLASTN profile.
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
        Err("unsupported TBLASTN outfmt: Stage E will implement 0, 6 and 7 output".into())
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

    // NCBI reference: c++/src/app/blast/tblastn_app.cpp:288-301
    // ```c++
    // results = lcl_blast.Run();
    // formatter.PrintOneResultSet(**result, query);
    // ```
    // Search and reporting are absent until Stages C-E, so no successful result is possible.
    pub fn run(self) -> Result<()> {
        self.validate()?;
        bail!("TBLASTN local search is unimplemented (Stages C-E)")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{try_parse_from, Cli, Commands};

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
