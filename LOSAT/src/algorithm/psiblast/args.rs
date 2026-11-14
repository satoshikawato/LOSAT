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
// PHI-BLAST arguments. LOSAT runs psiblast's first round as blastp's search (`run`), its later
// rounds (`-num_iterations`) with the PSSM made from each round (`rounds.rs`), and writes the
// ASN.1 and ASCII PSSMs (`pssm.rs`); with `-in_pssm` the
// query is a PSSM, searched with the PSSM's stages (`in_pssm.rs`, `algorithm::pssm_search`), in
// rounds and with the PSSM files too (`in_pssm_files.rs`);
// `-in_msa` and `-phi_pattern` are refused by the parser (`cli.rs` `is_unported_psiblast_arg`).
// The values are
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
        help = "BLAST database name (a protein database, with or without parsed Seq-ids)"
    )]
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist`,
    /// `-negative_seqidlist` and `-no_taxid_expansion`: blastp's (line DI-8), which hold in every
    /// round (NCBI's `CLocalDbAdapter` of the first round searches every later one too).
    #[command(flatten)]
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663
    // ```c++
    //     arg_desc.AddFlag(kArgShowGIs, "Show NCBI GIs in deflines?", true);
    // ```
    #[arg(
        long = "show_gis",
        default_value_t = false,
        help = "Show NCBI GIs in deflines?"
    )]
    pub show_gis: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1961-1962
    // (`CQueryOptionsArgs`, psiblast_args.cpp:86)
    // ```c++
    //     arg_desc.AddFlag(kArgParseDeflines,
    //                  "Should the query and subject defline(s) be parsed?", true);
    // ```
    #[arg(
        long = "parse_deflines",
        default_value_t = false,
        help = "Should the query and subject defline(s) be parsed?"
    )]
    pub parse_deflines: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    // ```c++
    //     arg_desc.AddOptionalKey(kArgEffSearchSpace, "int_value",
    //                             "Effective length of the search space",
    //                             CArgDescriptions::eInt8);
    //     arg_desc.SetConstraint(kArgEffSearchSpace,
    //                            new CArgAllowValuesGreaterThanOrEqual(0));
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
    // ```c++
    //     if (args.Exist(kArgEffSearchSpace) && args[kArgEffSearchSpace]) {
    //         CNcbiEnvironment env;
    //         env.Set("OLD_FSC", "true");
    //         opt.SetEffectiveSearchSpace(args[kArgEffSearchSpace].AsInt8());
    //     }
    // ```
    // psiblast's arguments include `CGenericSearchArgs` (psiblast_args.cpp:66-67). The value is
    // blastp's (`BlastpArgs::searchsp`): every round, those of `-in_pssm` too, then runs without
    // the Gumbel block (`ResolvedBlastpArgs::old_fsc`).
    #[arg(long, value_parser = nonnegative_ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-245, 321-323
    // ```c++
    //     if (!m_SuppressSumStats) {
    //         arg_desc.AddOptionalKey(kArgSumStats, "bool_value",
    //                      	 	"Use sum statistics",
    //                      	 	CArgDescriptions::eBoolean);
    //     }
    //     ...
    //     if (args.Exist(kArgSumStats) && args[kArgSumStats]) {
    //         opt.SetSumStatisticsMode(args[kArgSumStats].AsBoolean());
    //     }
    // ```
    // psiblast's `CGenericSearchArgs` keeps it (psiblast_args.cpp:67; blastp suppresses it). A
    // gapped protein search does not link HSPs by default (blast_options.c:1460-1471), so `false`
    // is psiblast's search; `true` links the HSPs of each subject in every round
    // (`BlastpArgs::sum_stats`, `core::link_hsps`).
    #[arg(long, value_parser = ncbi_boolean, help = "Use sum statistics")]
    pub sum_stats: Option<bool>,
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
    /// The search strategy file (`-export_search_strategy`, quoted at
    /// `BlastpArgs::export_search_strategy`), created by the search's `run` before the other options
    /// are read (with `-in_pssm` refused: item F1b-2 of line RG).
    #[arg(long, value_name = "filename", value_parser = ncbi_output_path(), help = "File name to record the search strategy used")]
    pub export_search_strategy: Option<PathBuf>,
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]", conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    // NCBI blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    #[command(flatten)]
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    #[arg(long = "max_hsps", value_parser = blastn_count)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long = "qcov_hsp_perc", value_parser = psiblast_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3301, 3338-3340
    // ```c++
    //     arg_desc.AddOptionalKey(kArgCullingLimit, "int_value",
    //                      "If the query range of a hit is enveloped by that of at "
    //                      "least this many higher-scoring hits, delete the hit",
    //                      CArgDescriptions::eInteger);
    //     arg_desc.SetConstraint(kArgCullingLimit,
    //                new CArgAllowValuesGreaterThanOrEqual(kDfltArgCullingLimit));
    //     ...
    //     if (args[kArgCullingLimit]) {
    //         opts.SetCullingLimit(args[kArgCullingLimit].AsInteger());
    //     }
    // ```
    // psiblast's `CHspFilteringArgs` (psiblast_args.cpp:79), blastp's (`BlastpArgs::culling_limit`):
    // every round culls, the PSSM rounds and `-in_pssm` too (`pssm_hooks::pssm_prelim`).
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit")]
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229, 276-286
    // (quoted at `BlastpArgs::xdrop_ungap`)
    // psiblast's `CGenericSearchArgs` (psiblast_args.cpp:67): blastp's X-drops in bits, which
    // every round converts with its own Karlin-Altschul blocks (`pssm_hooks`).
    #[arg(long = "xdrop_ungap", value_name = "float_value", value_parser = psiblast_real, help = "X-dropoff value (in bits) for ungapped extensions [default: 7]")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long = "xdrop_gap", value_name = "float_value", value_parser = psiblast_real, help = "X-dropoff value (in bits) for preliminary gapped extensions [default: 15]")]
    pub xdrop_gap: Option<f64>,
    #[arg(long = "xdrop_gap_final", value_name = "float_value", value_parser = psiblast_real, help = "X-dropoff value (in bits) for final gapped alignment [default: 25]")]
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1058-1079
    // ```c++
    //     const double default_value = m_QueryIsProtein
    //         ? BLAST_GAP_TRIGGER_PROT : BLAST_GAP_TRIGGER_NUCL;
    //     arg_desc.AddDefaultKey(kArgGapTrigger, "float_value",
    //                            "Number of bits to trigger gapping",
    //                            CArgDescriptions::eDouble,
    //                            NStr::DoubleToString(default_value));
    //     ...
    //     if (args[kArgGapTrigger]) {
    //         opt.SetGapTrigger(args[kArgGapTrigger].AsDouble());
    //     }
    // ```
    // psiblast's `CGapTriggerArgs` (psiblast_args.cpp:105; blastp has none): the cutoff of the
    // ungapped alignments that a gapped extension starts from, in every round
    // (`ResolvedBlastpArgs::gap_trigger_bits`).
    #[arg(long = "gap_trigger", value_name = "float_value", default_value = "22", value_parser = psiblast_real, help = "Number of bits to trigger gapping")]
    pub gap_trigger: f64,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1940-1942
    // ```c++
    //     // lowercase masking
    //     arg_desc.AddFlag(kArgUseLCaseMasking,
    //          "Use lower case filtering in query and subject sequence(s)?", true);
    // ```
    #[arg(long)]
    pub lcase_masking: bool,
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
        conflicts_with_all = crate::blastinput::db_filter::DB_FILTER_OPTIONS,
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
        help = "Number of iterations to perform (0 means run until convergence)"
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
    // Only the PSSM engine reads it (`PSIBlastOptions::pseudo_count`), which one round runs only
    // for `-save_pssm_after_last_round` (`pssm.rs`); without it the value changes nothing (NCBI
    // 2.17.0: the same reports with 5 and -1).
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1156-1160
    // ```c++
    //         // checkpoint file
    //         arg_desc.AddOptionalKey(kArgPSIOutputChkPntFile, "checkpoint_file",
    //
    //                                 "File name to store checkpoint file",
    //                                 CArgDescriptions::eOutputFile);
    // ```
    // The file is opened when a PSSM is written (`pssm.rs` `AsnPssmWriter`, as
    // `CAutoOutputFileReset`), not by the parser.
    #[arg(
        long = "out_pssm",
        value_name = "checkpoint_file",
        help = "File name to store checkpoint file"
    )]
    pub out_pssm: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1161-1170
    // ```c++
    //         // ASCII matrix file
    //         arg_desc.AddOptionalKey(kArgAsciiPssmOutputFile, "ascii_mtx_file",
    //                                 "File name to store ASCII version of PSSM",
    //                                 CArgDescriptions::eOutputFile);
    //
    //         arg_desc.AddFlag(kArgSaveLastPssm, "Save PSSM after the last database "
    //                          "search");
    //         arg_desc.AddFlag(kArgSaveAllPssms, "Save PSSM after each iteration "
    //                          "(file name is given in -save_pssm or "
    //                          "-save_ascii_pssm options)");
    // ```
    // The file is opened when a PSSM is written (`pssm.rs` `write_ascii_pssm_files`, as
    // `CAutoOutputFileReset`), not by the parser.
    #[arg(
        long = "out_ascii_pssm",
        value_name = "ascii_mtx_file",
        help = "File name to store ASCII version of PSSM"
    )]
    pub out_ascii_pssm: Option<String>,
    #[arg(
        long = "save_pssm_after_last_round",
        help = "Save PSSM after the last database search"
    )]
    pub save_pssm_after_last_round: bool,
    #[arg(
        long = "save_each_pssm",
        help = "Save PSSM after each iteration (file name is given in -save_pssm or -save_ascii_pssm options)"
    )]
    pub save_each_pssm: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1231-1247
    // ```c++
    //             // PSI-BLAST checkpoint
    //             arg_desc.AddOptionalKey(kArgPSIInputChkPntFile, "psi_chkpt_file",
    //                                     "PSI-BLAST checkpoint file",
    //                                     CArgDescriptions::eInputFile);
    //             ...
    //     if (!m_IsDeltaBlast) {
    //         arg_desc.SetDependency(kArgPSIInputChkPntFile,
    //                                CArgDescriptions::eExcludes,
    //                                kArgQuery);
    //         arg_desc.SetDependency(kArgPSIInputChkPntFile,
    //                                CArgDescriptions::eExcludes,
    //                                kArgQueryLocation);
    //     }
    // ```
    // The exclusions are errors of the argument parser: LOSAT's parser reports them with its own
    // text (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1). The file is opened and read with
    // the options (`algorithm::psiblast::in_pssm`); `-in_msa` and `-phi_pattern` are refused by
    // the parser (`cli.rs` `is_unported_psiblast_arg`).
    #[arg(
        long = "in_pssm",
        value_name = "psi_chkpt_file",
        value_parser = ncbi_input_path(),
        conflicts_with_all = ["query", "query_loc"],
        help = "PSI-BLAST checkpoint file"
    )]
    pub in_pssm: Option<PathBuf>,

    /// `BlastpArgs::first_query_ordinal`: the number of the first query of the records of a call
    /// that starts there (`api::run_batches`). 0 on the command line.
    #[arg(skip)]
    pub first_query_ordinal: usize,
}

impl PsiblastArgs {
    /// blastp's arguments for psiblast's search, after LOSAT's refusals of what it does not run
    /// (the rounds of a range, and the output formats of the later items of the series).
    pub fn into_blastp_args(self) -> Result<BlastpArgs> {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:281-287
        // ```c++
        //         if (pssm.NotEmpty())
        //         {
        //             _ASSERT(queries.Empty());
        //             if(num_iters != 0)
        //             	retval.Reset(new blast::CExportStrategy(pssm, opts_hndl, search_db, client_id, num_iters));
        // ```
        // A PSSM query's request holds the PSSM as read (`queries pssm`), which item F1b-2 of line
        // RG writes; refused before the file is made.
        if self.export_search_strategy.is_some() && self.in_pssm.is_some() {
            bail!(
                "-export_search_strategy with -in_pssm is not supported by LOSAT's PSIBLAST yet (NCBI writes a PSI-BLAST request that holds the PSSM with -db and stops with \"Remote PSI-BL2SEQ is not supported\" with -subject)"
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1277-1281
        // ```c++
        //     if (m_DbTarget == eProteinDb) {
        //         if (args[kArgPSINumIterations]) {
        //         	if(m_NumIterations == 1)
        //         		m_NumIterations = args[kArgPSINumIterations].AsInteger();
        //         }
        // ```
        // A query given as a PSSM (`-in_pssm`) runs its rounds as a FASTA query's, round 1 with the
        // input PSSM (item PS5b).
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:147-165
        // ```c++
        // s_GetQueryBioseq(CConstRef<CBlastQueryVector> query, CRef<CScope> scope,
        //                  CRef<CPssmWithParameters> pssm)
        // {
        //     CConstRef<CBioseq> retval;
        //
        //     if (pssm) {
        //         retval.Reset(&pssm->SetQuery().SetSeq());
        //     } else {
        //         _ASSERT(query.NotEmpty());
        //         if (query->GetQuerySeqLoc(0)->GetId()) {
        //             CBioseq_Handle bh =
        //                 scope->GetBioseqHandle(*query->GetQuerySeqLoc(0)->GetId());
        // ```
        // The PSSM of the next round is made of the whole record's Bioseq while the search keeps
        // the range: NCBI's later rounds align a PSSM of the record to the range (NCBI 2.17.0
        // `-query_loc 20-300 -num_iterations 6`: 4.494 % identities in round 2), an output that
        // LOSAT does not reproduce.
        if self.num_iterations != 1 && self.query_loc.is_some() {
            bail!(
                "-query_loc with -num_iterations {} (NCBI BLAST+'s later rounds align the PSSM of the whole record to the range) is not supported by LOSAT's PSIBLAST",
                self.num_iterations
            );
        }
        Self::check_outfmt(&self.outfmt)?;
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-248 (quoted at
        // `in_pssm_files`: the loop writes the PSSM of each round at its start, the input PSSM of
        // `-in_pssm` in round 1)
        // The round driver runs `-num_iterations` other than 1, and with `-in_pssm` one round too
        // when a PSSM file is written (the input PSSM before the search).
        let rounds = (self.num_iterations != 1
            || (self.in_pssm.is_some()
                && (self.out_pssm.is_some() || self.out_ascii_pssm.is_some())))
        .then(|| u32::try_from(self.num_iterations).unwrap_or(0));
        let request_inputs = self
            .export_search_strategy
            .is_some()
            .then(|| Box::new(self.option_inputs()));
        Ok(BlastpArgs {
            query: self.query,
            subject: self.subject,
            db: self.db,
            // blastp's database filters (line DI-8), which every round searches with.
            db_filter: self.db_filter,
            dbsize: self.dbsize,
            show_gis: self.show_gis,
            parse_deflines: self.parse_deflines,
            id_scope: None,
            searchsp: self.searchsp,
            lcase_masking: self.lcase_masking,
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
            // psiblast has no `-mt_mode` (`cli.rs` `is_unported_psiblast_arg`).
            mt_mode: None,
            out: self.out,
            max_target_seqs: self.max_target_seqs,
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            max_hsps_per_subject: self.max_hsps_per_subject,
            qcov_hsp_perc: self.qcov_hsp_perc,
            culling_limit: self.culling_limit,
            ungapped: false,
            xdrop_ungap: self.xdrop_ungap,
            xdrop_gap: self.xdrop_gap,
            xdrop_gap_final: self.xdrop_gap_final,
            gap_trigger: Some(self.gap_trigger),
            sum_stats: self.sum_stats == Some(true),
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
            export_search_strategy: self.export_search_strategy,
            strategy: None,
            request_inputs,
            app: BlastpApp::Psiblast,
            subject_taxdb: None,
            inclusion_ethresh: Some(self.inclusion_ethresh),
            first_query_ordinal: self.first_query_ordinal,
            psiblast_pssm: Some(crate::algorithm::psiblast::pssm::PsiblastPssmOptions {
                save_pssm_after_last_round: self.save_pssm_after_last_round,
                save_each_pssm: self.save_each_pssm,
                out_pssm: self.out_pssm,
                out_ascii_pssm: self.out_ascii_pssm,
                pseudocount: self.pseudocount,
                rounds,
                in_pssm: self.in_pssm,
            }),
            // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:502
            // ```c++
            //         formatter.SetLineLength(fmt_args->GetLineLength());
            // ```
            // The line length only: `-sorthits` and `-sorthsps` are read with their warnings
            // (`BlastpArgs::check_options`) and sort nothing (`ResolvedBlastpArgs::sort`).
            formatting: self.formatting,
            rps: None,
        })
    }

    /// The options that psiblast's search strategy requests record (`options::algorithm_options`
    /// of `RequestKind::Psiblast`: NCBI's argument classes of `psiblast_args.cpp` and their
    /// conditions), from the arguments as given. `HitlistSize` is the request's
    /// (`StrategyExport::for_run`); the database lists are the database's
    /// (`OptionsHandle::set_database_filter`).
    pub(crate) fn option_inputs(&self) -> crate::report::archive::options::OptionInputs {
        use crate::report::archive::options::{seg_filter, OptionInputs, RequestKind};
        let mut inputs = OptionInputs::new(RequestKind::Psiblast, "blastp");
        inputs.dbsize = self.dbsize;
        inputs.searchsp = self.searchsp;
        inputs.sum_stats = self.sum_stats;
        inputs.evalue = self.evalue;
        inputs.gapopen = self.gap_open;
        inputs.gapextend = self.gap_extend;
        inputs.matrix = self.matrix.clone();
        // `-xdrop_ungap` is never recorded (`SetXDropoff` has no remote option).
        inputs.xdrop_gap = self.xdrop_gap;
        inputs.xdrop_gap_final = self.xdrop_gap_final;
        inputs.word_size = self.word_size;
        inputs.qcov_hsp_perc = self.qcov_hsp_perc;
        inputs.max_hsps = self.max_hsps_per_subject.map(|count| count as i32);
        inputs.culling_limit = self.culling_limit;
        // `static_cast<int>` of the double, as BLASTP's (`blastp::archive::option_inputs`).
        inputs.threshold = self
            .threshold
            .map(crate::core::blast_util::ncbi_int4_from_double);
        inputs.window_size = self.window_size;
        inputs.soft_masking = self.soft_masking;
        inputs.seg = self
            .seg
            .as_deref()
            .and_then(|value| seg_filter(value, "PSIBLAST"));
        if let Ok((mode, unified_p)) =
            crate::blastinput::app::parse_comp_based_stats(&self.comp_based_stats, true, false)
        {
            inputs.comp_based_stats = mode as i32;
            inputs.unified_p = unified_p;
        }
        inputs.use_sw_tback = self.use_sw_tback;
        inputs.gap_trigger = self.gap_trigger;
        inputs.pseudocount = self.pseudocount;
        inputs.inclusion_ethresh = self.inclusion_ethresh;
        inputs.num_iterations = Some(self.num_iterations);
        inputs
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2744-2748
    // ```c++
    // CFormattingArgs::ArchiveFormatRequested(const CArgs& args) const
    // {
    //     EOutputFormat output_fmt;
    //     string ignore1, ignore2;
    //     ParseFormattingString(args, output_fmt, ignore1, ignore2);
    // ```
    /// NCBI parses `-outfmt` before the option handlers; its own errors come from that parse
    /// and from blastp's (`run`). The formats that NCBI runs and LOSAT's psiblast does not
    /// write yet are refused here. The tabular reports with a field list, and outfmt 10 and
    /// 20, are blastp's (`blastp::tabular`), without a custom delimiter, which psiblast's
    /// formatter does not get (psiblast_app.cpp:482-505, `blast_engine.rs` `run`).
    pub fn check_outfmt(outfmt: &str) -> Result<()> {
        let choice = crate::blastinput::app::parse_formatting_string(outfmt)?;
        if let Some(format) = crate::blastinput::app::ReportFormat::from_number(choice.number) {
            use crate::blastinput::app::ReportFormat;
            if !matches!(
                format,
                ReportFormat::Pairwise
                    | ReportFormat::Tabular
                    | ReportFormat::TabularWithComments
                    | ReportFormat::CommaSeparatedValues
                    | ReportFormat::CommaSeparatedValuesWithHeader
                    | ReportFormat::Sam
            ) {
                bail!(
                    "output format {} is not supported by LOSAT's PSIBLAST",
                    choice.number
                );
            }
        }
        Ok(())
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

    // blast_args.cpp:234-238 and 303-307, ncbiargs.cpp:118-129: `-searchsp` takes NCBI's `Int8`
    // forms (a sign, hex, leading zeros) at 0 or above, and at any value (0 too) the search runs
    // without the Gumbel block (NCBI 2.17.0 probes `x.sp.{plus,hex,lead0}`, `err.sp.*` and
    // `x.sp.{space,trail}` of line PS, item PS5c).
    #[test]
    fn searchsp_takes_ncbis_int8_forms_and_drops_the_gumbel_block() {
        for (value, space) in [
            ("+1000000000", 1_000_000_000),
            ("0x3B9ACA00", 1_000_000_000),
            ("0001000000000", 1_000_000_000),
            ("9223372036854775807", i64::MAX),
            ("0", 0),
        ] {
            let args = parse(&["-searchsp", value]).unwrap();
            assert_eq!(args.searchsp, Some(space), "{value}");
            let resolved = args.into_blastp_args().unwrap().resolve().unwrap();
            assert_eq!(resolved.searchsp, Some(space), "{value}");
            assert!(resolved.old_fsc, "{value}");
        }
        let blastp = parse(&[]).unwrap().into_blastp_args().unwrap();
        assert_eq!(blastp.searchsp, None);
        for value in [" 1", "1 ", "1e6", "1.5", "", "-1", "9223372036854775808"] {
            assert!(
                parse(&["-searchsp", value]).is_err(),
                "-searchsp {value:?} is refused"
            );
        }
    }

    // psiblast_app.cpp:246-321 (`rounds.rs`): the rounds of a FASTA query run (0: until
    // convergence), and those of `-in_pssm` (PS5b, with one round too when a PSSM file is
    // written); with `-query_loc` they are refused.
    #[test]
    fn rounds_run_and_the_later_output_formats_are_refused() {
        for (rounds, plan) in [("0", Some(0)), ("2", Some(2)), ("1", None)] {
            let blastp = parse(&["-num_iterations", rounds])
                .unwrap()
                .into_blastp_args()
                .unwrap();
            assert_eq!(blastp.psiblast_pssm.as_ref().unwrap().rounds, plan);
            let resolved = blastp.resolve().unwrap();
            assert_eq!(
                resolved.psi_rounds.as_ref().map(|plan| plan.num_iterations),
                plan
            );
        }
        let in_pssm = |extra: &[&str]| {
            let mut argv = vec![
                "losat", "psiblast", "-in_pssm", "x.asn", "-subject", "s.faa",
            ];
            argv.extend_from_slice(extra);
            let Ok(Cli {
                command: Commands::Psiblast(args),
            }) = crate::cli::try_parse_from::<Cli, _, _>(argv)
            else {
                panic!("psiblast -in_pssm {extra:?}");
            };
            args.into_blastp_args()
                .unwrap()
                .psiblast_pssm
                .unwrap()
                .rounds
        };
        assert_eq!(in_pssm(&["-num_iterations", "3"]), Some(3));
        assert_eq!(in_pssm(&[]), None);
        assert_eq!(in_pssm(&["-out_pssm", "f"]), Some(1));
        assert_eq!(in_pssm(&["-out_ascii_pssm", "f"]), Some(1));
        let error = parse(&["-num_iterations", "2", "-query_loc", "1-20"])
            .unwrap()
            .into_blastp_args()
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("-query_loc with -num_iterations 2"),
            "{error}"
        );
        let error = parse(&["-num_iterations", "3", "-seg", "yes"])
            .unwrap()
            .into_blastp_args()
            .unwrap()
            .resolve()
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("-seg and -soft_masking true with -num_iterations 3"),
            "{error}"
        );
        for format in ["1", "3", "5", "8", "11", "13", "15", "16", "18"] {
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
        for format in [
            "0",
            "6",
            "7",
            "6 qseqid sseqid evalue",
            "6 delim=| std",
            "10",
            "20 qseqid sseqid",
            "17",
        ] {
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
    fn blastps_database_filters_and_id_options_reach_the_search() {
        let db = |extra: &[&str]| {
            let mut argv = vec!["losat", "psiblast", "-query", "q.faa", "-db", "d"];
            argv.extend_from_slice(extra);
            match crate::cli::try_parse_from::<Cli, _, _>(argv) {
                Ok(Cli {
                    command: Commands::Psiblast(args),
                }) => Ok(args.into_blastp_args().unwrap()),
                Ok(other) => panic!("not psiblast: {other:?}"),
                Err(error) => Err(crate::cli::render_message(&error)),
            }
        };
        let blastp = db(&[
            "-taxids",
            "9606",
            "-no_taxid_expansion",
            "-show_gis",
            "-parse_deflines",
            "-lcase_masking",
            "-num_iterations",
            "3",
        ])
        .unwrap();
        assert_eq!(blastp.db_filter.taxids.as_deref(), Some("9606"));
        assert!(blastp.db_filter.no_taxid_expansion);
        assert!(blastp.show_gis && blastp.parse_deflines && blastp.lcase_masking);
        let blastp = db(&["-negative_seqidlist", "l.txt"]).unwrap();
        assert_eq!(
            blastp.db_filter.negative_seqidlist.as_deref(),
            Some("l.txt")
        );
        let blastp = db(&[]).unwrap();
        assert!(!blastp.show_gis && !blastp.parse_deflines && !blastp.lcase_masking);
        assert!(blastp.db_filter.taxids.is_none() && blastp.db_filter.seqidlist.is_none());
        // The filters exclude `-subject` and `-subject_loc`, as blastp's do.
        for extra in [
            &["-subject", "s.faa", "-taxids", "9606"][..],
            &["-subject_loc", "1-10", "-seqidlist", "l.txt"][..],
        ] {
            let mut argv = vec!["losat", "psiblast", "-query", "q.faa"];
            argv.extend_from_slice(extra);
            let error = crate::cli::try_parse_from::<Cli, _, _>(argv).unwrap_err();
            assert!(
                crate::cli::render_message(&error).contains("cannot be used with"),
                "{extra:?}"
            );
        }
    }

    #[test]
    fn the_x_drops_gap_trigger_culling_and_sum_stats_reach_blastps_search() {
        let blastp = |extra: &[&str]| parse(extra).unwrap().into_blastp_args();
        let args = blastp(&[
            "-culling_limit",
            "2",
            "-xdrop_ungap",
            "3.5",
            "-xdrop_gap",
            "1e9",
            "-xdrop_gap_final",
            "-1",
            "-gap_trigger",
            "-5",
            "-sum_stats",
            "false",
        ])
        .unwrap();
        assert_eq!(args.culling_limit, Some(2));
        assert_eq!(
            (args.xdrop_ungap, args.xdrop_gap, args.xdrop_gap_final),
            (Some(3.5), Some(1e9), Some(-1.0))
        );
        assert_eq!(args.gap_trigger, Some(-5.0));
        // NCBI's default gap trigger of a protein search (BLAST_GAP_TRIGGER_PROT).
        let args = blastp(&[]).unwrap();
        assert_eq!(args.gap_trigger, Some(22.0));
        assert_eq!(args.culling_limit, None);
        assert!(!args.sum_stats);
        // Sum statistics link the HSPs, in every spelling of NCBI's boolean.
        for spelling in ["true", "T", "t", "1", "yes"] {
            assert!(
                blastp(&["-sum_stats", spelling]).unwrap().sum_stats,
                "{spelling}"
            );
        }
        assert!(!blastp(&["-sum_stats", "false"]).unwrap().sum_stats);
        // NCBI's constraint of -culling_limit (0 or more).
        assert!(parse(&["-culling_limit", "-1"]).is_err());
    }

    #[test]
    fn psiblasts_arguments_that_losat_does_not_run_are_named() {
        for (option, value) in [
            ("-in_msa", Some("x")),
            ("-msa_master_idx", Some("2")),
            ("-ignore_msa_master", None),
            ("-phi_pattern", Some("x")),
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
        // `-in_pssm` (PS4a) is an argument: with `-query` it is the parser's exclusion
        // (blast_args.cpp:1240-1247).
        let error = parse(&["-in_pssm", "x"]).unwrap_err();
        assert!(error.contains("cannot be used with"), "{error}");
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
