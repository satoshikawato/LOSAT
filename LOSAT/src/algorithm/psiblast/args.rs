//! Command-line arguments for PSIBLAST

use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
use crate::blastinput::value_parsers::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:58-116
// ```c++
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
//     m_BlastDbArgs->SetIPGFilteringSupport(true);
//     arg.Reset(m_BlastDbArgs);
//     m_Args.push_back(arg);
//
//     m_StdCmdLineArgs.Reset(new CStdCmdLineArgs);
//     arg.Reset(m_StdCmdLineArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CGenericSearchArgs(kQueryIsProtein));
//     m_Args.push_back(arg);
//
//     arg.Reset(new CFilteringArgs(kQueryIsProtein, kFilterByDefault));
//     m_Args.push_back(arg);
//
//     arg.Reset(new CMatrixNameArg);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CWordThresholdArg);
//     m_Args.push_back(arg);
//
//     m_HspFilteringArgs.Reset(new CHspFilteringArgs);
//     arg.Reset(m_HspFilteringArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CWindowSizeArg);
//     m_Args.push_back(arg);
//
//     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
//     arg.Reset(m_QueryOptsArgs);
//     m_Args.push_back(arg);
//
//     m_FormattingArgs.Reset(new CFormattingArgs);
//     arg.Reset(m_FormattingArgs);
//     m_Args.push_back(arg);
//
//     m_MTArgs.Reset(new CMTArgs);
//     arg.Reset(m_MTArgs);
//     m_Args.push_back(arg);
//
//     m_RemoteArgs.Reset(new CRemoteArgs);
//     arg.Reset(m_RemoteArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CCompositionBasedStatsArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CGapTriggerArgs(kQueryIsProtein));
//     m_Args.push_back(arg);
//
//     m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eProteinDb));
//     arg.Reset(m_PsiBlastArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CPssmEngineArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CPhiBlastArgs);
//     m_Args.push_back(arg);
// ```
// psiblast's arguments are blastp's without `-task` (its task is "psiblast"), `-ungapped`
// (no `CGappedArgs`), the database masks and `-mt_mode`, with the PSI-BLAST, PSSM engine and
// PHI-BLAST arguments. LOSAT runs psiblast's first round (`-num_iterations 1`), which is blastp's
// search (`run`); the PSSM files, `-in_msa`, `-in_pssm`, `-phi_pattern`, `-gap_trigger` and
// `-sum_stats` are refused by the parser (`cli.rs` `is_unported_psiblast_arg`). The values are
// read as NCBI's argument parser reads them (`value_parsers.rs`); the checks that NCBI makes
// after reading them are blastp's (`BlastpArgs::check_options`, with psiblast's program).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct PsiblastArgs {
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // ```c
    //     arg_desc.AddDefaultKey(kArgQuery, "input_file",
    //                      "Input file name",
    //                      CArgDescriptions::eInputFile, kDfltArgQuery);
    // ```
    // `-` (the default) is standard input.
    #[arg(long, default_value = "-", value_parser = ncbi_input_path(), value_name = "PATH")]
    pub query: PathBuf,
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastinput/app.rs` `missing_subject_error`).
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // ```c++
    //         arg_desc.AddOptionalKey(kArgSubject, "subject_input_file",
    //                                 "Subject sequence(s) to search",
    //                                 CArgDescriptions::eInputFile);
    //         ITERATE(vector<string>, dbarg, database_args) {
    //             arg_desc.SetDependency(kArgSubject, CArgDescriptions::eExcludes,
    //                                    *dbarg);
    //         }
    // ```
    // The exclusion is an error of the argument parser: LOSAT's parser reports it with its own
    // text (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1).
    #[arg(
        long,
        value_name = "NAME",
        conflicts_with_all = ["subject", "subject_loc"],
        help = "BLAST database name (a protein database without parsed Seq-ids)"
    )]
    pub db: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // ```c++
    //     if (!m_IsMapper) {
    //         arg_desc.SetCurrentGroup("Statistical options");
    //         arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
    //                                 "Effective length of the database ",
    //                                 CArgDescriptions::eInt8);
    //     }
    // ```
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS", help = "Effective length of the database")]
    pub dbsize: Option<i64>,
    #[arg(long, value_parser = psiblast_real, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = psiblast_threshold, help = "Minimum word score such that the word is added to the BLAST lookup table [default: 11, or the matrix's suggestion]")]
    pub threshold: Option<f64>,
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3]")]
    pub word_size: Option<i32>,
    #[arg(long, default_value_t = 1, value_parser = blastn_count)]
    pub num_threads: usize,
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]", conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    // NCBI blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    #[arg(long = "max_hsps", value_parser = blastn_count)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long = "qcov_hsp_perc", value_parser = psiblast_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40, or the matrix's suggestion]")]
    pub window_size: Option<i32>,
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11, or the matrix's best value]")]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1, or the matrix's best value]")]
    pub gap_extend: Option<i32>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // ```c
    //     arg_desc.AddDefaultKey(kArgCompBasedStats, "compo", legend,
    //                            CArgDescriptions::eString, m_DefaultOpt);
    // ```
    #[arg(
        long,
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
    )]
    pub comp_based_stats: String,
    // NCBI's filtering handler (`blastinput/app.rs` `parse_seg_option`); no filtering when
    // omitted (`CPSIBlastOptionsHandle::SetQueryOptionDefaults`).
    #[arg(long, help = "SEG: no, yes, or \"WINDOW LOCUT HICUT\" [default: no]")]
    pub seg: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // ```c++
    //         arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                         "Apply filtering locations as soft masks",
    //                         CArgDescriptions::eBoolean,
    //                         kDfltArgLookupTableMaskingOnlyProt);
    // ```
    // psiblast's filtering handler is blastp's (`CFilteringArgs`, psiblast_args.cpp:70).
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: false]")]
    pub soft_masking: Option<bool>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:803-806
    // ```c
    //     arg_desc.AddFlag(kArgUseSWTraceback,
    //                      "Compute locally optimal Smith-Waterman alignments?",
    //                      true);
    // ```
    #[arg(long)]
    pub use_sw_tback: bool,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // ```c++
    //     arg_desc.AddOptionalKey(kArgQueryLocation, "range",
    //                             "Location on the query sequence in 1-based offsets "
    //                             "(Format: start-stop)",
    //                             CArgDescriptions::eString);
    // ```
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // ```c++
    //         arg_desc.AddOptionalKey(kArgSubjectLocation, "range",
    //                         "Location on the subject sequence in 1-based offsets "
    //                         "(Format: start-stop)",
    //                         CArgDescriptions::eString);
    // ```
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    // NCBI's argument is a string (blast_args.cpp:2657-2660), parsed before the option
    // handlers (`blastinput/app.rs` `parse_formatting_string`).
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1147-1152
    // ```c++
    //         arg_desc.AddDefaultKey(kArgPSINumIterations, "int_value",
    //                                "Number of iterations to perform (0 means run "
    //                                "until convergence)", CArgDescriptions::eInteger,
    //                                NStr::IntToString(kDfltArgPSINumIterations));
    //         arg_desc.SetConstraint(kArgPSINumIterations,
    //                                new CArgAllowValuesGreaterThanOrEqual(0));
    // ```
    #[arg(
        long = "num_iterations",
        default_value = "1",
        value_name = "int_value",
        value_parser = nonnegative_ncbi_integer,
        help = "Number of iterations to perform (0 means run until convergence); LOSAT runs 1"
    )]
    pub num_iterations: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1086-1090
    // ```c++
    //     // Pseudo count
    //     arg_desc.AddDefaultKey(kArgPSIPseudocount, "pseudocount",
    //                            "Pseudo-count value used when constructing PSSM",
    //                            CArgDescriptions::eInteger,
    //                            NStr::IntToString(PSI_PSEUDO_COUNT_CONST));
    // ```
    // Only the PSSM engine reads it (`PSIBlastOptions::pseudo_count`), which one round without a
    // PSSM file does not run: the value changes nothing (NCBI 2.17.0: the same reports with 5 and
    // -1).
    #[arg(
        long,
        default_value = "0",
        value_parser = ncbi_integer,
        help = "Pseudo-count value used when constructing PSSM"
    )]
    pub pseudocount: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1100-1104
    // ```c++
    //     // Evalue inclusion threshold
    //     arg_desc.AddDefaultKey(kArgPSIInclusionEThreshold, "ethresh",
    //                    "E-value inclusion threshold for pairwise alignments",
    //                    CArgDescriptions::eDouble,
    //                    NStr::DoubleToString(PSI_INCLUSION_ETHRESH));
    // ```
    // The first round reads it too: it is the cutoff of the composition-based redo's heaps
    // (`ResolvedBlastpArgs::inclusion_ethresh`).
    #[arg(
        long,
        default_value = "0.002",
        value_parser = psiblast_real,
        help = "E-value inclusion threshold for pairwise alignments"
    )]
    pub inclusion_ethresh: f64,
}

impl PsiblastArgs {
    /// blastp's arguments for psiblast's first round, after LOSAT's refusals of what it does not
    /// run (more rounds, and the output formats of the later items of the series).
    pub fn into_blastp_args(self) -> Result<BlastpArgs> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1277-1281
        // ```c++
        //     if (m_DbTarget == eProteinDb) {
        //         if (args[kArgPSINumIterations]) {
        //         	if(m_NumIterations == 1)
        //         		m_NumIterations = args[kArgPSINumIterations].AsInteger();
        //         }
        // ```
        if self.num_iterations != 1 {
            bail!(
                "-num_iterations {} (rounds after the first, which search with a PSSM) is not supported by LOSAT's PSIBLAST",
                self.num_iterations
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2744-2748
        // ```c++
        // CFormattingArgs::ArchiveFormatRequested(const CArgs& args) const
        // {
        //     EOutputFormat output_fmt;
        //     string ignore1, ignore2;
        //     ParseFormattingString(args, output_fmt, ignore1, ignore2);
        // ```
        // NCBI parses `-outfmt` before the option handlers; its own errors come from that parse
        // and from blastp's (`run`). The formats that NCBI runs and LOSAT's psiblast does not
        // write yet are refused here.
        let choice = crate::blastinput::app::parse_formatting_string(&self.outfmt)?;
        if let Some(format) = crate::blastinput::app::ReportFormat::from_number(choice.number) {
            use crate::blastinput::app::ReportFormat;
            if !matches!(
                format,
                ReportFormat::Pairwise
                    | ReportFormat::Tabular
                    | ReportFormat::TabularWithComments
                    | ReportFormat::Sam
            ) {
                bail!(
                    "output format {} is not supported by LOSAT's PSIBLAST",
                    choice.number
                );
            }
        }
        Ok(BlastpArgs {
            query: self.query,
            subject: self.subject,
            db: self.db,
            dbsize: self.dbsize,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:46-50
            // ```c++
            // CPSIBlastOptionsHandle::CPSIBlastOptionsHandle(EAPILocality locality)
            //     : CBlastAdvancedProteinOptionsHandle(locality)
            // {
            //     SetDefaults();
            //     m_Opts->SetProgram(ePSIBlast);
            // ```
            // psiblast's options are those of blastp's task (`CBlastAdvancedProteinOptionsHandle`)
            // with no filtering (`SetQueryOptionDefaults`, as blastp's).
            task: "blastp".to_string(),
            evalue: self.evalue,
            threshold: self.threshold,
            word_size: self.word_size,
            num_threads: self.num_threads,
            out: self.out,
            max_target_seqs: self.max_target_seqs,
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            max_hsps_per_subject: self.max_hsps_per_subject,
            qcov_hsp_perc: self.qcov_hsp_perc,
            // psiblast has no -culling_limit (rejected above).
            culling_limit: None,
            ungapped: false,
            window_size: self.window_size,
            matrix: self.matrix,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            comp_based_stats: self.comp_based_stats,
            seg: self.seg,
            soft_masking: self.soft_masking,
            use_sw_tback: self.use_sw_tback,
            query_loc: self.query_loc,
            subject_loc: self.subject_loc,
            outfmt: self.outfmt,
            archive: None,
            app: BlastpApp::Psiblast,
            inclusion_ethresh: Some(self.inclusion_ethresh),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};

    fn parse(extra: &[&str]) -> Result<PsiblastArgs, String> {
        let mut argv = vec!["losat", "psiblast", "-query", "q.faa", "-subject", "s.faa"];
        argv.extend_from_slice(extra);
        match crate::cli::try_parse_from::<Cli, _, _>(argv) {
            Ok(Cli {
                command: Commands::Psiblast(args),
            }) => Ok(args),
            Ok(other) => panic!("not psiblast: {other:?}"),
            Err(error) => Err(crate::cli::render_message(&error)),
        }
    }

    // psiblast_args.cpp:44-118 and psiblast_options.cpp:46-60: blastp's task options, the
    // program psiblast, and the PSSM engine's defaults (inclusion threshold 0.002, pseudo-count 0,
    // one round).
    #[test]
    fn the_first_round_is_blastps_search_with_psiblasts_program() {
        let args = parse(&[]).unwrap();
        assert_eq!(args.num_iterations, 1);
        assert_eq!(args.pseudocount, 0);
        assert_eq!(args.inclusion_ethresh, 0.002);
        let blastp = args.into_blastp_args().unwrap();
        assert_eq!(blastp.app, BlastpApp::Psiblast);
        assert_eq!(blastp.task, "blastp");
        assert!(!blastp.ungapped);
        assert_eq!(blastp.inclusion_ethresh, Some(0.002));
        let resolved = blastp.resolve().unwrap();
        assert_eq!(resolved.app, BlastpApp::Psiblast);
        assert_eq!(resolved.inclusion_ethresh, 0.002);
        assert_eq!(resolved.query_batch_size().unwrap(), 1);
        // traceback_stage.cpp:216-225: the preliminary hit list size of CBS (1050 for 500).
        assert_eq!(resolved.max_target_seqs, 500);
        assert_eq!(resolved.traceback_hitlist_size(), 1050);
        let resolved = parse(&["-max_target_seqs", "600", "-inclusion_ethresh", "1e-300"])
            .unwrap()
            .into_blastp_args()
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!(resolved.traceback_hitlist_size(), 1250);
        assert_eq!(resolved.inclusion_ethresh, 1e-300);
        // blastp keeps its hit list size and its threshold in the traceback.
        let mut blastp = parse(&[]).unwrap().into_blastp_args().unwrap();
        blastp.app = BlastpApp::Blastp;
        blastp.inclusion_ethresh = None;
        let resolved = blastp.resolve().unwrap();
        assert_eq!(resolved.traceback_hitlist_size(), 500);
        assert_eq!(resolved.inclusion_ethresh, 0.002);
    }

    #[test]
    fn more_rounds_and_the_later_output_formats_are_refused() {
        for rounds in ["0", "2"] {
            let error = parse(&["-num_iterations", rounds])
                .unwrap()
                .into_blastp_args()
                .unwrap_err()
                .to_string();
            assert!(
                error.contains(&format!("-num_iterations {rounds} (rounds after the first")),
                "{error}"
            );
        }
        for format in ["1", "3", "5", "8", "10", "11", "13", "15", "16", "18"] {
            let error = parse(&["-outfmt", format])
                .unwrap()
                .into_blastp_args()
                .unwrap_err()
                .to_string();
            assert_eq!(
                error,
                format!("output format {format} is not supported by LOSAT's PSIBLAST")
            );
        }
        for format in ["0", "6", "7", "6 qseqid sseqid evalue", "17"] {
            assert!(parse(&["-outfmt", format])
                .unwrap()
                .into_blastp_args()
                .is_ok());
        }
        // NCBI's own error of the format string comes first (NCBI 2.17.0 `psiblast -outfmt 99`:
        // "Error: Formatting choice is out of range", exit 255).
        assert!(parse(&["-outfmt", "99"])
            .unwrap()
            .into_blastp_args()
            .unwrap_err()
            .to_string()
            .contains("Formatting choice is out of range"));
    }

    #[test]
    fn psiblasts_arguments_that_losat_does_not_run_are_named() {
        for (option, value) in [
            ("-in_pssm", Some("x")),
            ("-out_pssm", Some("x")),
            ("-out_ascii_pssm", Some("x")),
            ("-in_msa", Some("x")),
            ("-msa_master_idx", Some("2")),
            ("-ignore_msa_master", None),
            ("-save_pssm_after_last_round", None),
            ("-save_each_pssm", None),
            ("-phi_pattern", Some("x")),
            ("-gap_trigger", Some("22")),
            ("-sum_stats", Some("false")),
            ("-culling_limit", Some("1")),
            ("-html", None),
        ] {
            let mut extra = vec![option];
            extra.extend(value);
            let error = parse(&extra).unwrap_err();
            assert!(
                error.contains(&format!(
                    "the NCBI BLAST+ option {option} is not supported by LOSAT's PSIBLAST"
                )),
                "{option}: {error}"
            );
        }
        // Arguments that psiblast does not have are unknown, as in NCBI's parser.
        for option in ["-task", "-ungapped", "-db_soft_mask", "-mt_mode"] {
            let error = parse(&[option, "1"]).unwrap_err();
            assert!(
                error.contains("unknown option or argument"),
                "{option}: {error}"
            );
        }
    }
}
