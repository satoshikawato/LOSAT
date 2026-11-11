//! BLASTX options and native serial local-subject search/report dispatch.
use crate::blastinput::value_parsers::{
    blastx_query_coverage, blastx_real, blastx_threshold, genetic_code, ncbi_boolean,
    ncbi_input_path, ncbi_int8, ncbi_integer, ncbi_output_path, ncbi_positive_integer,
    nonnegative_ncbi_integer, protein_word_size, report_count, SegSpec,
};
use crate::stats::protein_options::{
    protein_gap_existence_extend_params, suggested_threshold, suggested_window_size,
    SuggestionProgram,
};
use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

/// `-max_target_seqs` as the hit list count of `blastinput::report_counts`: a value below 1, which a
/// caller of `run_local` can set (the command line's is at least 1), counts 0 (NCBI's hit list size
/// below 1, "No hits are being saved").
pub(crate) fn max_target_seqs_count(max_target_seqs: Option<i32>) -> Option<usize> {
    max_target_seqs.map(|n| usize::try_from(n).unwrap_or(0))
}

/// BLASTX's `-task`: the command line's parser, and the check of `run_local` (`TaskParser::check`).
pub(crate) const BLASTX_TASK: crate::blastinput::value_parsers::TaskParser =
    crate::blastinput::value_parsers::task_parser("BLASTX", &["blastx", "blastx-fast"], &[]);

#[derive(Args, Debug)]
#[command(rename_all = "snake_case")]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:49-62
// ```c++
//                                   "Translated Query-Protein Subject BLAST"));
//     const bool kQueryIsProtein = false;
//     m_Args.push_back(arg);
//     m_ClientId = kProgram + " " + CBlastVersion().Print();
//
//     static const char kDefaultTask[] = "blastx";
//     SetTask(kDefaultTask);
//     set<string> tasks;
//     tasks.insert(kDefaultTask);
//     tasks.insert("blastx-fast");
//     arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
//     m_Args.push_back(arg);
//
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:72-78
// ```c++
//     // N.B.: query is not protein because the options are applied on the
//     // translated query
//     arg.Reset(new CGenericSearchArgs( !kQueryIsProtein ));
//     m_Args.push_back(arg);
//
//     //Disable until OOF is supported in align manager
//     //SB-1043
// ```
pub struct BlastxArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // ```c++
    //     arg_desc.AddDefaultKey(kArgQuery, "input_file",
    //                      "Input file name",
    //                      CArgDescriptions::eInputFile, kDfltArgQuery);
    // ```
    // The default `-` is standard input (`input_files::open_input`).
    #[arg(long, value_name = "PATH", default_value = "-", value_parser = ncbi_input_path())]
    pub query: PathBuf,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2366
    // ```c++
    //         arg_desc.AddOptionalKey(kArgSubject, "subject_input_file",
    //                                 "Subject sequence(s) to search",
    //                                 CArgDescriptions::eInputFile);
    // ```
    // Optional for the parser: without `-subject` and `-db` the run fails with NCBI's error
    // after the parsing (`blastinput/app.rs` `missing_subject_error`).
    #[arg(long, value_name = "PATH", value_parser = ncbi_input_path())]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // ```c++
    //         // subject location
    //         arg_desc.AddOptionalKey(kArgSubjectLocation, "range",
    //                         "Location on the subject sequence in 1-based offsets "
    //                         "(Format: start-stop)",
    //                         CArgDescriptions::eString);
    // ```
    // Read by the database arguments handler when it reads the subjects (`native::run`).
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        conflicts_with_all = crate::blastinput::db_filter::DB_FILTER_OPTIONS,
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
    // ```c++
    //     arg_desc.SetCurrentGroup("General search options");
    //     // database filename
    //     if (m_IsIgBlast){
    //         arg_desc.AddOptionalKey(kArgDb, "database_name", "Optional additional database name",
    //                                 CArgDescriptions::eString);
    //     } else {
    //         arg_desc.AddOptionalKey(kArgDb, "database_name", "BLAST database name",
    //                                 CArgDescriptions::eString);
    //     }
    // ```
    // NCBI excludes `-subject` and `-subject_loc` with `-db` (blast_args.cpp:2364-2381). The
    // exclusion is an error of the argument parser: LOSAT's parser reports it with its own text
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1).
    #[arg(
        long,
        value_name = "NAME",
        conflicts_with_all = ["subject", "subject_loc"],
        help = "BLAST database name (a protein database, with or without parsed Seq-ids)"
    )]
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    #[command(flatten)]
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // ```c++
    //     // DB size
    //     if (!m_IsMapper) {
    //         arg_desc.SetCurrentGroup("Statistical options");
    //         arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
    //                                 "Effective length of the database ",
    //                                 CArgDescriptions::eInt8);
    //     }
    // ```
    // No constraint and no dependency: read with `-db` and with `-subject`.
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2050-2052
    // ```c++
    //     arg_desc.AddDefaultKey(kArgParseDeflines, "TF", "Should the query and "
    // ```
    // `-parse_deflines` is a flag of the programs (blast_args.cpp:1961).
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
    // ```
    #[arg(long, value_parser = crate::blastinput::value_parsers::nonnegative_ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3441-3443
    // ```c++
    //     arg_desc.AddDefaultKey(kArgOutput, "output_file",
    //                    "Output file name",
    //                    CArgDescriptions::eOutputFile, "-");
    // ```
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:99-100
    // ```c++
    //         arg_desc.AddDefaultKey(kTask, "task_name", "Task to execute",
    //                                CArgDescriptions::eString, m_DefaultTask);
    // ```
    #[arg(long, default_value = "blastx", value_parser = BLASTX_TASK)]
    pub task: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:146-146
    // ```c++
    //         arg_desc.AddOptionalKey(kArgEvalue, "evalue", des, CArgDescriptions::eDouble);
    // ```
    #[arg(long, default_value_t = 10.0, value_parser = blastx_real)]
    pub evalue: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-167
    // ```c++
    //         arg_desc.AddOptionalKey(kArgWordSize, "int_value", description,
    //                                 CArgDescriptions::eInteger);
    // ```
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3; 5 for blastx-fast]")]
    pub word_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:203-205
    // ```c++
    //         arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value",
    //                            "Set maximum number of HSPs per subject sequence to save for each query",
    //                            CArgDescriptions::eInteger);
    // ```
    #[arg(long, value_parser = ncbi_positive_integer)]
    pub max_hsps: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201
    // ```c++
    //         arg_desc.AddOptionalKey(kArgQueryCovHspPerc, "float_value",
    //                             "Percent query coverage per hsp",
    //                             CArgDescriptions::eDouble);
    //         arg_desc.SetConstraint(kArgQueryCovHspPerc,
    //                            new CArgAllow_Doubles(0.0, 100.0));
    // ```
    #[arg(long, value_parser = blastx_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:241-243
    // ```c++
    //         arg_desc.AddOptionalKey(kArgSumStats, "bool_value",
    //                      	 	"Use sum statistics",
    //                      	 	CArgDescriptions::eBoolean);
    // ```
    #[arg(long, value_parser = ncbi_boolean, action = clap::ArgAction::Set, help = "Use sum statistics [default: true]")]
    pub sum_stats: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:175-177
    // ```c++
    //         arg_desc.AddOptionalKey(kArgGapOpen, "open_penalty",
    //                                 "Cost to open a gap",
    //                                 CArgDescriptions::eInteger);
    // ```
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11]")]
    pub gap_open: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:180-182
    // ```c++
    //         arg_desc.AddOptionalKey(kArgGapExtend, "extend_penalty",
    //                                "Cost to extend a gap",
    //                                CArgDescriptions::eInteger);
    // ```
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1]")]
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:937-942
    // ```c++
    //     arg_desc.AddDefaultKey(kArgMaxIntronLength, "length",
    //                     "Length of the largest intron allowed in a translated "
    //                     "nucleotide sequence when linking multiple distinct "
    //                     "alignments",
    //                     CArgDescriptions::eInteger,
    //                     NStr::IntToString(kDfltArgMaxIntronLength));
    // ```
    #[arg(long, default_value_t = 0, value_parser = nonnegative_ncbi_integer)]
    pub max_intron_length: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-338
    // ```c++
    //         arg_desc.AddDefaultKey(kArgSegFiltering, "SEG_options",
    //                         "Filter query sequence with SEG "
    //                         "(Format: '" + kDfltArgApplyFiltering + "', " +
    //                         "'window locut hicut', or '" + kDfltArgNoFiltering +
    //                         "' to disable)",
    //                         CArgDescriptions::eString, m_FilterByDefault
    //                         ? kDfltArgSegFiltering : kDfltArgNoFiltering);
    // ```
    // Read when NCBI's filtering handler reads it (`seg_spec`).
    #[arg(long, default_value = "12 2.2 2.5")]
    pub seg: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // ```c++
    //         arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                         "Apply filtering locations as soft masks",
    //                         CArgDescriptions::eBoolean,
    //                         kDfltArgLookupTableMaskingOnlyProt);
    // ```
    #[arg(long, default_value = "false", value_parser = ncbi_boolean, action = clap::ArgAction::Set)]
    pub soft_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:531-534
    // ```c++
    //     arg_desc.AddDefaultKey(kArgMatrixName, "matrix_name",
    //                            "Scoring matrix name",
    //                            CArgDescriptions::eString,
    //                            string(""));
    // ```
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:578-583
    // ```c++
    //     arg_desc.AddOptionalKey(kArgWordScoreThreshold, "float_value",
    //                  "Minimum word score such that the word is added to the "
    //                  "BLAST lookup table",
    //                  CArgDescriptions::eDouble);
    //     arg_desc.SetConstraint(kArgWordScoreThreshold,
    //                            new CArgAllowValuesGreaterThanOrEqual(0));
    // ```
    #[arg(
        long,
        value_parser = blastx_threshold,
        help = "Minimum word score such that the word is added to the BLAST lookup table [default: 12; 19.3, 21 or 20.25 for word sizes 5, 6 and 7; 20 for blastx-fast]"
    )]
    pub threshold: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:475-478
    // ```c++
    //     arg_desc.AddOptionalKey(kArgWindowSize, "int_value",
    //                             "Multiple hits window size, use 0 to specify "
    //                             "1-hit algorithm",
    //                             CArgDescriptions::eInteger);
    // ```
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40]")]
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3299
    // ```c++
    //     arg_desc.AddOptionalKey(kArgCullingLimit, "int_value",
    //                      "If the query range of a hit is enveloped by that of at "
    //                      "least this many higher-scoring hits, delete the hit",
    //                      CArgDescriptions::eInteger);
    // ```
    #[arg(
        long,
        value_parser = nonnegative_ncbi_integer,
        help = "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit [default: 0, no culling]"
    )]
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329-3329
    // ```c++
    //     arg_desc.AddFlag(kArgSubjectBestHit, "Return only the best HSP for each non overlapping query region", true);
    // ```
    #[arg(long)]
    pub subject_besthit: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917
    // ```c++
    //     arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
    // ```
    #[arg(long)]
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1941-1942
    // ```c++
    //     arg_desc.AddFlag(kArgUseLCaseMasking,
    //          "Use lower case filtering in query and subject sequence(s)?", true);
    // ```
    #[arg(long)]
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
    // ```c++
    //         arg_desc.AddDefaultKey(kArgStrand, "strand",
    //                      "Query strand(s) to search against database/subject",
    //                                CArgDescriptions::eString, kDfltArgStrand);
    // ```
    #[arg(long, default_value = "both", value_parser = ["both", "plus", "minus"])]
    pub strand: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // ```c++
    //     // query location
    //     arg_desc.AddOptionalKey(kArgQueryLocation, "range",
    //                             "Location on the query sequence in 1-based offsets "
    //                             "(Format: start-stop)",
    //                             CArgDescriptions::eString);
    // ```
    // Read by the query options handler (`native::run`).
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1021-1024
    // ```c++
    //         arg_desc.AddDefaultKey(kArgQueryGeneticCode, "int_value",
    //                                "Genetic code to use to translate query (see https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes for details)\n",
    //                                CArgDescriptions::eInteger,
    //                                NStr::IntToString(BLAST_GENETIC_CODE));
    // ```
    #[arg(long, default_value_t = 1, value_parser = genetic_code)]
    pub query_gencode: u8,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
    // ```c++
    //     arg_desc.AddDefaultKey(kArgOutputFormat, "format",
    //                            kOutputFormatDescription,
    //                            CArgDescriptions::eString,
    //                            NStr::IntToString(dft_outfmt));
    // ```
    // NCBI reads the value when it sets the options (`native::run`; `blastinput/app.rs`
    // `parse_formatting_string`).
    #[arg(long, default_value = "0")]
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2726-2730
    // ```c++
    //         arg_desc.AddOptionalKey(kArgMaxTargetSequences, "num_sequences",
    //                             "Maximum number of aligned sequences to keep \n"
    //     						"(value of 5 or more is recommended)\n"
    //     						"Default = `" + NStr::IntToString(BLAST_HITLIST_SIZE) + "'",
    //                             CArgDescriptions::eInteger);
    // ```
    #[arg(
        long,
        value_parser = ncbi_positive_integer,
        help = "Maximum number of aligned sequences to keep [default: 500]",
        conflicts_with_all = ["num_descriptions", "num_alignments"]
    )]
    pub max_target_seqs: Option<i32>,
    // NCBI blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3158-3161
    // ```c++
    //     arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
    //                            "Number of threads (CPUs) to use in the BLAST search",
    //                            CArgDescriptions::eInteger,
    //                            NStr::IntToString(kDfltValue));
    // ```
    #[arg(long, default_value_t = 1, value_parser = ncbi_positive_integer, help = "Number of threads (CPUs); serial WASI accepts only 1")]
    pub num_threads: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3169-3181
    // ```c++
    //     if (m_MTMode >= 0) {
    //     	arg_desc.AddDefaultKey(kArgMTMode, "int_value",
    //                                "Multi-thread mode to use in BLAST search:\n "
    //                                "0 auto split by database or queries \n "
    //                                "1 split by queries\n "
    //                                "2 split by database",
    //                                CArgDescriptions::eInteger,
    //                                NStr::IntToString(0));
    //     	arg_desc.SetConstraint(kArgMTMode,
    //                                new CArgAllowValuesBetween(0, 2, true));
    //     	arg_desc.SetDependency(kArgMTMode,
    //                                CArgDescriptions::eRequires,
    //                                kArgNumThreads);
    //     }
    // ```
    /// `-mt_mode` is read and has no effect: LOSAT's output is the same at every thread count
    /// and for every value, and equals NCBI's `-num_threads 1` output (Owner 2026-10-13;
    /// CLI DECISIONS 5). NCBI's `-mt_mode 1` (and `0` when it chooses to split by queries)
    /// runs the queries in batches and writes the end of an outfmt 0 report or each XML
    /// document per batch (`blastn_app.cpp:134-170`, `blastx_app.cpp` `x_RunMTBySplitQuery`);
    /// that framing is not reproduced. The value must be 0 to 2 and `-num_threads` must be
    /// given with it, as in NCBI (clap's `requires` ignores a default value).
    #[arg(long, value_name = "int_value", value_parser = crate::blastinput::value_parsers::ncbi_mt_mode, requires = "num_threads", help = "Multi-thread mode: 0 auto, 1 split by queries, 2 split by database (accepted; the output does not depend on it)")]
    pub mt_mode: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // ```c++
    //     arg_desc.AddDefaultKey(kArgCompBasedStats, "compo", legend,
    //                            CArgDescriptions::eString, m_DefaultOpt);
    // ```
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`: the first
    // character chooses the mode, and the 'u' of unified P is BLASTP's only).
    #[arg(
        long,
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
    )]
    pub comp_based_stats: String,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:838-871
// ```c++
//         switch (comp_stat_string[0]) {
//             case '0': case 'F': case 'f':
//                 compo_mode = eNoCompositionBasedStats;
//                 break;
//             case '1':
//                 compo_mode = eCompositionBasedStats;
//                 break;
//             case 'D': case 'd':
//                 if ((program == eRPSBlast) || (program == eRPSTblastn)) {
//                     compo_mode = eNoCompositionBasedStats;
//                 }
//                 else if (program == eDeltaBlast) {
//                     compo_mode = eCompositionBasedStats;
//                 }
//                 else {
//                     compo_mode = eCompositionMatrixAdjust;
//                 }
//                 break;
//             case '2':
//                 compo_mode = eCompositionMatrixAdjust;
//                 break;
//             case '3':
//                 compo_mode = eCompoForceFullMatrixAdjust;
//                 break;
//             case 'T': case 't':
//                 compo_mode = (program == eRPSBlast || program == eRPSTblastn || program == eDeltaBlast) ?
//                     eCompositionBasedStats : eCompositionMatrixAdjust;
//                 break;
//         }
//
//         if(program == ePSITblastn) {
//             compo_mode = eNoCompositionBasedStats;
//         }
//
// ```
/// The matrix of a `-matrix` name, as NCBI's matrix loader finds it (`None`: a name it does not
/// know).
///
/// NCBI reference (598d8ae6): c++/src/util/tables/raw_scoremat.c:131-166
/// ```c
/// const SNCBIPackedScoreMatrix* NCBISM_GetStandardMatrix(const char* name)
/// {
///     switch (name[0]) {
///     case 'B': case 'b':
///         if ( !NCBI_strncasecmp(name, "BLOSUM", 6) ) {
/// ...
///     case 'I': case 'i':
///         if ( !NCBI_strncasecmp(name, "IDENTITY", 8) ) {
///             return &NCBISM_Identity;
/// ```
/// (a name that starts with IDENTITY but is not IDENTITY is IDENTITY's matrix without its
/// statistics table: `identity_prefix`.)
pub(crate) fn blastx_matrix(name: &str) -> Option<crate::config::ScoringMatrix> {
    name.parse().ok()
}

/// Whether NCBI's matrix loader reads `name` as IDENTITY (`blastx_matrix`).
fn identity_prefix(name: &str) -> bool {
    name.get(..8)
        .is_some_and(|start| start.eq_ignore_ascii_case("identity"))
}

/// The composition mode of `-comp_based_stats` as NCBI reads it for BLASTX: 0 to 3 from the first
/// character (`blastinput/app.rs` `parse_comp_based_stats`, the switch above), any other first
/// character 0, and a 'u' after it ignored (unified P is BLASTP's).
pub(crate) fn composition_mode(value: &str) -> u8 {
    use crate::blastinput::app::{parse_comp_based_stats, CompositionMode};
    match parse_comp_based_stats(value, false, false) {
        Ok((CompositionMode::CompositionBasedStats, _)) => 1,
        Ok((CompositionMode::CompositionMatrixAdjust, _)) => 2,
        Ok((CompositionMode::CompoForceFullMatrixAdjust, _)) => 3,
        Ok((CompositionMode::NoCompositionBasedStats, _)) | Err(_) => 0,
    }
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:43-49
// ```c++
// const size_t kNumTabularOutputFormatSpecifiers = 50;
// const SFormatSpec sc_FormatSpecifiers[kNumTabularOutputFormatSpecifiers] = {
//     SFormatSpec("qseqid",
//                 "Query Seq-id",
//                 eQuerySeqId),
//     SFormatSpec("qgi",
//                 "Query GI",
// ```
/// The fields that ABI v1 accepts (`web::v1_outfmt`, frozen, plan TD-1): the custom fields
/// of BLASTX's report before it wrote them through the shared tabular writer.
pub const FIELDS: &[&str] = &[
    "qseqid", "qacc", "qaccver", "qlen", "sseqid", "sacc", "saccver", "slen", "qstart", "qend",
    "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident",
    "mismatch", "positive", "gapopen", "gaps", "ppos", "qframe", "sframe", "frames", "btop",
    "stitle",
];

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blastx_options.cpp:57-70
// ```c++
//     CBlastProteinOptionsHandle::SetLookupTableDefaults();
//     m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTX);
// }
//
// void
// CBlastxOptionsHandle::SetQueryOptionDefaults()
// {
//     CBlastProteinOptionsHandle::SetQueryOptionDefaults();
//     m_Opts->SetStrandOption(objects::eNa_strand_both);
//     m_Opts->SetQueryGeneticCode(BLAST_GENETIC_CODE);
//     SetSegFiltering(false); // disable SEG filtering because of eCompositionMatrixAdjust mode
// }
//
// void
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:131-146
// ```c++
// }
//
// void
// CBlastProteinOptionsHandle::SetHitSavingOptionsDefaults()
// {
//     SetHitlistSize(500);
//     SetEvalueThreshold(BLAST_EXPECT_VALUE);
//     SetMinDiagSeparation(0);
//     SetPercentIdentity(0);
//     // set some default here, allow INT4MAX to mean infinity
//     SetMaxNumHspPerSequence(0);
//     SetMaxHspsPerSubject(0);
//
//     SetCutoffScore(0); // will be calculated based on evalue threshold,
//     // effective lengths and Karlin-Altschul params in BLAST_Cutoffs_simple
//     // and passed to the engine in the params structure
// ```
pub struct ResolvedOptions {
    pub word_size: i32,
    pub threshold: f64,
    pub lookup_type: i32,
    /// The `-matrix` name as typed (BLOSUM62 when omitted), which the reports print.
    pub matrix: String,
    /// The matrix of `matrix` (BLOSUM62 for a name that NCBI's matrices do not know, which only an
    /// ungapped search takes: `unknown_matrix`).
    pub scoring: crate::config::ScoringMatrix,
    /// An ungapped search's `-matrix` that NCBI's matrix loader does not know: NCBI's set-up of the
    /// first query batch fails with "Unknown error code -1" (`query_setup.rs`).
    pub unknown_matrix: bool,
    /// The file that NCBI BLAST+ reads for `unknown_matrix`'s name
    /// (`ncbi_environment::ApplicationSettings::matrix_file`), found by the command line after
    /// the options are resolved; `None` elsewhere. The first query batch is refused then
    /// (`query_setup::matrix_error_message`).
    pub matrix_file: Option<std::path::PathBuf>,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gapped: bool,
    pub composition: u8,
    pub seg: SegSpec,
    pub soft_masking: bool,
    pub strand: String,
    pub query_gencode: u8,
    pub window_size: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
    // ```c++
    //
    // void
    // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
    // {
    //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
    // }
    //
    // void
    // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
    // {
    //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
    //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
    //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:123-145
    // ```c++
    //                                              searches except blastn */
    // #define BLAST_UNGAPPED_X_DROPOFF_NUCL 20 /**< ungapped dropoff score for
    //                                               blastn (and megablast) */
    //
    // /** default dropoff for preliminary gapped extensions */
    // #define BLAST_GAP_X_DROPOFF_PROT 15 /**< default dropoff (all protein-
    //                                          based gapped extensions) */
    // #define BLAST_GAP_X_DROPOFF_NUCL 30 /**< default dropoff for non-greedy
    //                                          nucleotide gapped extensions */
    // #define BLAST_GAP_X_DROPOFF_GREEDY 25 /**< default dropoff for greedy
    //                                          nucleotide gapped extensions */
    // #define BLAST_GAP_X_DROPOFF_TBLASTX 0 /**< default dropoff for tblastx */
    //
    // /** default bit score that will trigger gapped extension */
    // #define BLAST_GAP_TRIGGER_PROT 22.0 /**< default bit score that will trigger
    //                                          a gapped extension for all protein-
    //                                          based searches */
    // #define BLAST_GAP_TRIGGER_NUCL 27.0  /**< default bit score that will trigger
    //                                          a gapped extension for blastn */
    //
    // /** default dropoff for the final gapped extension with traceback */
    // #define BLAST_GAP_X_DROPOFF_FINAL_PROT 25 /**< default dropoff (all protein-
    //                                                based gapped extensions) */
    // ```
    pub x_dropoff: f64,
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub gap_trigger: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
    // ```c++
    // void
    // CBlastProteinOptionsHandle::SetEffectiveLengthsOptionsDefaults()
    // {
    //     SetDbLength(0);
    //     SetDbSeqNum(0);
    //     SetEffectiveSearchSpace(0);
    // }
    // ```
    pub db_length: i64,
    pub effective_search_space: i64,
    /// Whether `-searchsp` was given, at any value, or the environment variable `OLD_FSC` is set
    /// (`app::old_fsc_environment`): NCBI then has no Gumbel block (`gumbel_block`).
    pub old_fsc: bool,

    pub evalue: f64,
    pub sum_stats: bool,
    pub max_intron_length: i32,
    pub max_hsps: i32,
    /// The -qcov_hsp_perc value, or 0 (no filter).
    pub query_cov_hsp_perc: f64,
    pub culling_limit: i32,
    pub subject_besthit: bool,
    pub hitlist_size: i32,
    pub num_descriptions: i32,
    pub num_alignments: i32,
    pub outfmt: u8,
    /// The fields of the default tabular report (outfmt 6 and 7 without `tabular`: the `std`
    /// fields), by name; empty for outfmt 0.
    pub fields: Vec<String>,
    /// The tabular report with a field list or a delimiter (outfmt 6 and 7), or outfmt 10 or
    /// 20 (`TabularSpec::of`), which the shared writer writes (`tabular.rs`); `None` for the
    /// default outfmt 6 and 7 and the other formats.
    pub tabular: Option<&'static crate::report::tabular::TabularSpec>,
    /// Whether `taxdb.bti` is on NCBI's search path, for a `tabular` field list that shows a
    /// plural name of the taxonomy database (its formatters warn without it); set by
    /// `native::run` (`-db`: the database's search path, `-subject`: NCBI's). False otherwise,
    /// as Web reads it (DECISIONS 7 of line CLI C2).
    pub taxdb_found: bool,
    /// A list kept no sequence of the database (`SeqDb::filtered_empty`): the queries of the
    /// searched batches have `EMPTY_DATABASE_MESSAGE`.
    pub filtered_empty: bool,
    pub num_threads: i32,
    /// The `-query_loc` range as given, which the tabular formats print as `qlen` and divide
    /// `qcovs` by; read by `native::run` where NCBI's query options handler reads it (after the
    /// subjects), `None` otherwise.
    pub query_range: Option<crate::blastinput::seq_range::SequenceRange>,
}

/// NCBI's `eAaLookupTable` (`ELookupTableType`, include/algo/blast/core/blast_options.h): the
/// standard lookup table, of word sizes 2 to 4.
pub(crate) const AA_LOOKUP_TABLE: i32 = 3;
/// NCBI's `eCompressedAaLookupTable`: the compressed lookup table, of word sizes 5 to 7.
pub(crate) const COMPRESSED_AA_LOOKUP_TABLE: i32 = 4;

impl ResolvedOptions {
    /// Whether the search's lookup table is the compressed one (`-word_size` 5 to 7, blastx-fast).
    pub(crate) fn compressed_lookup(&self) -> bool {
        self.lookup_type == COMPRESSED_AA_LOOKUP_TABLE
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3629-3642
// ```c++
//     CRef<CBlastOptionsHandle> retval(x_CreateOptionsHandle(locality, args));
//     CBlastOptions& opts = retval->SetOptions();
//     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
//         (*arg)->ExtractAlgorithmOptions(args, opts);
//     }
//
//     m_IsUngapped = !opts.GetGappedMode();
//     try { retval->Validate(); }
//     catch (const CBlastException& e) {
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
//     }
//     return retval;
// }
//
// ```
impl BlastxArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3629-3642
    // ```c++
    //     CRef<CBlastOptionsHandle> retval(x_CreateOptionsHandle(locality, args));
    //     CBlastOptions& opts = retval->SetOptions();
    //     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    //         (*arg)->ExtractAlgorithmOptions(args, opts);
    //     }
    //
    //     m_IsUngapped = !opts.GetGappedMode();
    //     try { retval->Validate(); }
    //     catch (const CBlastException& e) {
    //         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    //     }
    //     return retval;
    // }
    //
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2940-2947
    // ```c++
    //     	    m_LineLength = args[kArgLineLength].AsInteger();
    //     	}
    //         if(args.Exist(kArgSortHits) && args[kArgSortHits])
    //         {
    //        	    m_HitsSortOption = args[kArgSortHits].AsInteger();
    //         }
    //     }
    //     else
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1310
    // ```c++
    //     if (program_number != eBlastTypeBlastn &&
    //         program_number != eBlastTypeMapping &&
    //         (!Blast_ProgramIsRpsBlast(program_number)) &&
    //         options->threshold <= 0)
    //     {
    //         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //                          "Non-zero threshold required");
    //         return BLASTERR_OPTION_VALUE_INVALID;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
    // ```c++
    // 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
    // 	{
    // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //          "expect value or cutoff score must be greater than zero");
    // 		return BLASTERR_OPTION_VALUE_INVALID;
    // 	}
    // ```
    pub fn resolve(&self) -> Result<ResolvedOptions> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2844-2851
        // ```c++
        //         fmt_type = static_cast<EOutputFormat>(val);
        //         if ( !(fmt_type == eTabular ||
        //                fmt_type == eTabularWithComments ||
        //                fmt_type == eCommaSeparatedValues ||
        //                fmt_type == eCommaSeparatedValuesWithHeader ||
        //                fmt_type == eSAM) ) {
        //                custom_fmt_spec.clear();
        //         }
        // ```
        // The format as NCBI reads it (`native::run` reads it first, with NCBI's errors).
        let choice = crate::blastinput::app::parse_formatting_string(&self.outfmt)?;
        self.resolve_with_format(&choice)
    }

    /// `resolve` with the format already read (Web ABI v1 reads it its own way, `web.rs`).
    pub(crate) fn resolve_with_format(
        &self,
        choice: &crate::blastinput::app::FormatChoice,
    ) -> Result<ResolvedOptions> {
        // The command line's parser has checked the values below, a caller of `run_local` has
        // not (the task: `run_local`, first).
        crate::algorithm::blastn::args::QueryStrand::parse(&self.strand)?;
        let outfmt = u8::try_from(choice.number)?;
        if self
            .subject
            .as_deref()
            .is_some_and(|subject| subject.as_os_str() == "-")
        {
            bail!("unsupported BLASTX stdin subject");
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:872-877
        // ```c++
        //         if (ungapped && *ungapped && compo_mode != eNoCompositionBasedStats) {
        //             NCBI_THROW(CInputException, eInvalidInput,
        //                        "Composition-adjusted searched are not supported with "
        //                        "an ungapped search, please add -comp_based_stats F or "
        //                        "do a gapped search");
        //         }
        // ```
        // The filtering handler (`CFilteringArgs`, before the composition handler,
        // blastx_args.cpp:44-82): NCBI's error of `-seg` here, LOSAT's refusal of a value that
        // NCBI reads (a locut or hicut that is not finite) after NCBI's checks.
        let (seg, seg_refusal) = match self.seg_spec() {
            Ok(seg) => (seg, None),
            Err(crate::blastinput::app::SegOptionError::Ncbi(message)) => {
                return Err(SemanticOptionsError::from(message).into())
            }
            Err(crate::blastinput::app::SegOptionError::Losat(refusal)) => {
                (SegSpec::No, Some(refusal))
            }
        };
        let composition = composition_mode(&self.comp_based_stats);
        if self.ungapped && composition != 0 {
            return Err(SemanticOptionsError::from("Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a gapped search").into());
        }
        // The matrix as typed (BLOSUM62 when omitted) and its gap costs: with `-matrix`, a cost
        // not given is the matrix's BEST one.
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:255-275
        // ```c++
        //     if (args[kArgMatrixName]) {
        //         opt.SetMatrixName(args[kArgMatrixName].AsString().c_str());
        //         int gap_open=0, gap_extend=0;
        //         BLAST_GetProteinGapExistenceExtendParams(
        //                     args[kArgMatrixName].AsString().c_str(),
        //                     &gap_open, &gap_extend);
        //         if (!args[kArgGapOpen])
        //             opt.SetGapOpeningCost(gap_open);
        //         if (!args[kArgGapExtend])
        //             opt.SetGapExtensionCost(gap_extend);
        //     }
        // ```
        let matrix_name = self.matrix.clone().unwrap_or_else(|| "BLOSUM62".into());
        let (best_open, best_extend) = match &self.matrix {
            Some(name) => protein_gap_existence_extend_params(name).unwrap_or((0, 0)),
            None => (11, 1),
        };
        let gap_open = self.gap_open.unwrap_or(best_open);
        let gap_extend = self.gap_extend.unwrap_or(best_extend);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:913-924
        // ```c
        //     if (program_number != eBlastTypeBlastn && options->gapped_calculation)
        //     {
        //        ...
        //        Boolean std_matrix_only =
        //                   (program_number != eBlastTypeBlastp &&
        //                    program_number != eBlastTypeTblastn);
        //        Int2 status = Blast_KarlinBlkGappedLoadFromTables(kbp, options->gap_open,
        //                    options->gap_extend, options->matrix, std_matrix_only);
        // ```
        // A gapped search takes a matrix of NCBI's tables but IDENTITY (BLASTX reads the standard
        // matrices only), with one of the matrix's gap costs, NCBI's errors otherwise (after the
        // composition check above and before the threshold's); an ungapped one takes any name,
        // and any gap costs.
        if !self.ungapped {
            use crate::stats::protein_options::{
                karlin_blk_gapped_load_from_tables, print_allowed_values, print_matrix_message,
                TableStatus,
            };
            match karlin_blk_gapped_load_from_tables(gap_open, gap_extend, &matrix_name, true) {
                Ok(_) => {}
                Err(TableStatus::UnknownMatrix) => {
                    return Err(SemanticOptionsError::from(print_matrix_message(
                        &matrix_name,
                        true,
                    ))
                    .into())
                }
                Err(TableStatus::UnknownGapCosts) => {
                    return Err(SemanticOptionsError::from(print_allowed_values(
                        &matrix_name,
                        gap_open,
                        gap_extend,
                    ))
                    .into())
                }
            }
        }
        let scoring = blastx_matrix(&matrix_name);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1765-1773
        // ```c++
        //    if ((status = LookupTableOptionsValidate(program_number,
        //                     lookup_options, blast_msg)) != 0)
        //        return status;
        //    if ((status = BlastInitialWordOptionsValidate(program_number,
        //                     word_options, blast_msg)) != 0)
        //        return status;
        //    if ((status = BlastHitSavingOptionsValidate(program_number, hit_options,
        //                                                blast_msg)) != 0)
        //        return status;
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1310
        // ```c++
        //     if (program_number != eBlastTypeBlastn &&
        //         program_number != eBlastTypeMapping &&
        //         (!Blast_ProgramIsRpsBlast(program_number)) &&
        //         options->threshold <= 0)
        //     {
        //         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
        //                          "Non-zero threshold required");
        //         return BLASTERR_OPTION_VALUE_INVALID;
        // ```
        if self.threshold.is_some_and(|x| x <= 0.0) {
            return Err(SemanticOptionsError::from("Non-zero threshold required").into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1312-1348
        // ```c++
        //     if (options->word_size <= 0)
        //     {
        //         if ( !Blast_ProgramIsRpsBlast(program_number)) {
        //             Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
        //                                      "Word-size must be greater than zero");
        // ...
        //         if (program_number == eBlastTypeBlastp ||
        //             program_number == eBlastTypeTblastn ||
        //             program_number == eBlastTypeBlastx)
        //         {
        //             if (options->word_size > 7) {
        //                 Blast_MessageWrite(blast_msg, eBlastSevError,
        //                                    kBlastMessageNoContext,
        //                                    "Word-size must be less than "
        //                                    "8 for a tblastn, blastp or blastx search");
        // ```
        // The command line's `-word_size` is at least 2 (`protein_word_size`); a caller of
        // `run_local` sets the field itself.
        if let Some(word) = self.word_size {
            if word <= 0 {
                return Err(
                    SemanticOptionsError::from("Word-size must be greater than zero").into(),
                );
            }
            if word > 7 {
                return Err(SemanticOptionsError::from(
                    "Word-size must be less than 8 for a tblastn, blastp or blastx search",
                )
                .into());
            }
        }
        // The hit list size and the counts of the pairwise report
        // (`blastinput::report_counts`). A `-max_target_seqs` below 1 (a caller of `run_local`
        // sets it; the command line's is at least 1) is NCBI's hit list size below 1.
        let counts = crate::blastinput::report_counts::report_counts(
            i32::from(outfmt),
            max_target_seqs_count(self.max_target_seqs),
            self.num_descriptions,
            self.num_alignments,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1511-1517
        // (the hit list size is checked before the e-value)
        // ```c++
        // 	if (options->hitlist_size < 1)
        // 	{
        // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
        //                          "No hits are being saved");
        // 		return BLASTERR_OPTION_VALUE_INVALID;
        // 	}
        // ```
        if counts.hitlist_size < 1 {
            return Err(SemanticOptionsError::from("No hits are being saved").into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
        // ```c++
        // 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
        // 	{
        // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
        //          "expect value or cutoff score must be greater than zero");
        // 		return BLASTERR_OPTION_VALUE_INVALID;
        // 	}
        // ```
        if self.evalue <= 0.0 {
            return Err(SemanticOptionsError::from(
                "expect value or cutoff score must be greater than zero",
            )
            .into());
        }
        if let Some(refusal) = seg_refusal {
            return Err(refusal);
        }
        // LOSAT's limit, for values that NCBI's checks above accept: an e-value that is not
        // finite (`+inf`, `1e999`, `-nan`; NCBI's argument reads `inf` and `nan` without a sign as
        // no number). Every word size of NCBI's, 2 to 7, runs (item SC9 of line NP).
        if !self.evalue.is_finite() {
            bail!(
                "an -evalue that is not finite ({}) is not supported by LOSAT's BLASTX",
                self.evalue
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:429-433
        // ```c++
        //         if (task == "blastx-fast") {
        //             opts->SetWordSize(5);
        //             opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
        //             opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST);
        //         }
        // ```
        let fast = self.task == "blastx-fast";
        let word_size = self.word_size.unwrap_or(if fast { 5 } else { 3 });
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:288-300
        // ```c++
        //     if ( args.Exist(kArgWordSize) && args[kArgWordSize]) {
        //         if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
        //            opt.SetLookupTableType(eCompressedAaLookupTable);
        //            opt.SetWordThreshold(19.3);
        //            if (args[kArgWordSize].AsInteger() > 5) {
        //                opt.SetWordThreshold(21.0);
        //            }
        //            if (args[kArgWordSize].AsInteger() > 6) {
        //                opt.SetWordThreshold(20.25);
        //            }
        //         }
        //         opt.SetWordSize(args[kArgWordSize].AsInteger());
        // ```
        // (BLASTX's query is protein for these arguments: blastx_args.cpp:74 passes
        // `!kQueryIsProtein`.) A word size given on the command line above 4 is the compressed
        // table with its own threshold (19.3, 21 or 20.25), and one of 4 or less leaves the task's:
        // `-task blastx-fast -word_size 3` is the standard table with threshold 20
        // (api/blast_options_local_priv.hpp:612-619).
        let task_threshold = if fast { 20.0 } else { 12.0 };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:604-622
        // ```c++
        //     if (args[kArgWordScoreThreshold]) {
        //         opt.SetWordThreshold(args[kArgWordScoreThreshold].AsDouble());
        //     } else if (s_IsDefaultWordThreshold(opt.GetProgram(),
        //                                         opt.GetWordThreshold())) {
        //         double threshold = -1;
        //         BLAST_GetSuggestedThreshold(opt.GetProgramType(),
        //                                     opt.GetMatrixName(),
        //                                     &threshold);
        // ```
        // The task's threshold 12 (`BLAST_WORD_THRESHOLD_BLASTX`) is the matrix's suggestion
        // (BLOSUM62's 11 + 1 for a translated query); blastx-fast's 20 and `-word_size 5`'s 19.3
        // are not the default, and stay.
        let threshold = match self.threshold {
            Some(threshold) => threshold,
            None if self.word_size == Some(5) => 19.3,
            None if self.word_size == Some(6) => 21.0,
            None if self.word_size == Some(7) => 20.25,
            None if fast => task_threshold,
            None => suggested_threshold(SuggestionProgram::TranslatedQuery, &matrix_name),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:485-497
        // ```c++
        //     if (args[kArgWindowSize]) {
        //         opt.SetWindowSize(args[kArgWindowSize].AsInteger());
        //     } else {
        //         int window = -1;
        //         BLAST_GetSuggestedWindowSize(opt.GetProgramType(),
        //                                      opt.GetMatrixName(),
        //                                      &window);
        // ```
        let window_size = self
            .window_size
            .unwrap_or_else(|| suggested_window_size(&matrix_name));
        // NCBI 2.17.0 crashes on a name that its matrix loader reads as IDENTITY without IDENTITY's
        // statistics table (`IDENTITYX`), and on IDENTITY with the compressed lookup table
        // (`-word_size` 5 to 7, blastx-fast): LOSAT refuses both (ungapped searches only: a gapped
        // one stops at NCBI's matrix error first).
        if scoring.is_none() && identity_prefix(&matrix_name) {
            bail!(
                "-matrix {matrix_name} is read as IDENTITY without IDENTITY's statistics table (NCBI BLAST+ 2.17.0 aborts); this is not supported by LOSAT's BLASTX"
            );
        }
        crate::blastinput::protein_options_check::check_identity_compressed_lookup(
            &matrix_name,
            word_size > 4,
            word_size,
            "blastx",
            "BLASTX",
        )?;
        // The counts are at most `i32::MAX` (the options are `Int4`).
        let hitlist_size = counts.hitlist_size as i32;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:76-99
        // ```c++
        //     vector<string> format_tokens;
        //     NStr::Split(format, " ", format_tokens);
        //
        //     if (format_tokens.empty())
        //         x_AddDefaultFieldsToShow();
        //
        //     ITERATE (vector<string>, iter, format_tokens) {
        //         if (*iter == kDfltArgTabularOutputFmtTag)
        //             x_AddDefaultFieldsToShow();
        //         else if ((*iter)[0] == '-') {
        // ...
        //     if (m_FieldsToShow.empty()) {
        //         x_AddDefaultFieldsToShow();
        //     }
        // ```
        // The fields as NCBI reads them (`report::tabular::parse_fields`: unknown words are
        // dropped, `-name` removes a field). A list that is not the `std` fields with a tab, and
        // outfmt 10 and 20, are written by the shared writer (`tabular`); the default report
        // writes the `std` fields (`fields`).
        let fields = if outfmt != 0 {
            crate::report::tabular::parse_fields(&choice.spec)
                .into_iter()
                .map(|field| field.name().to_string())
                .collect()
        } else {
            Vec::new()
        };
        let tabular = crate::blastinput::app::ReportFormat::from_number(choice.number)
            .and_then(|format| crate::report::tabular::TabularSpec::of(format, choice))
            .map(crate::report::tabular::TabularSpec::intern);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:43-70 (the
        // preliminary hit list size through NCBI's `int`; `blastn/hsp.rs`)
        // A hit list size whose preliminary size is not positive crashes NCBI (-max_target_seqs
        // 1073741799 to 2147483622 with composition-based statistics: SIGSEGV).
        crate::algorithm::blastn::hsp::check_prelim_hitlist_size(
            counts.hitlist_size,
            self.max_target_seqs.is_some(),
            composition > 0,
            !self.ungapped,
            "BLASTX",
        )?;
        Ok(ResolvedOptions {
            word_size,
            threshold,
            lookup_type: if word_size > 4 {
                COMPRESSED_AA_LOOKUP_TABLE
            } else {
                AA_LOOKUP_TABLE
            },
            matrix: matrix_name,
            scoring: scoring.unwrap_or(crate::config::ScoringMatrix::Blosum62),
            unknown_matrix: scoring.is_none(),
            matrix_file: None,
            gap_open,
            gap_extend,
            gapped: !self.ungapped,
            composition,
            seg,
            soft_masking: self.soft_masking,
            strand: self.strand.clone(),
            query_gencode: self.query_gencode,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
            // ```c++
            //
            // void
            // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
            // {
            //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
            //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
            // }
            //
            // void
            // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
            // {
            //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
            //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
            //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
            //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
            //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
            // }
            // ```
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:123-145
            // ```c++
            //                                              searches except blastn */
            // #define BLAST_UNGAPPED_X_DROPOFF_NUCL 20 /**< ungapped dropoff score for
            //                                               blastn (and megablast) */
            //
            // /** default dropoff for preliminary gapped extensions */
            // #define BLAST_GAP_X_DROPOFF_PROT 15 /**< default dropoff (all protein-
            //                                          based gapped extensions) */
            // #define BLAST_GAP_X_DROPOFF_NUCL 30 /**< default dropoff for non-greedy
            //                                          nucleotide gapped extensions */
            // #define BLAST_GAP_X_DROPOFF_GREEDY 25 /**< default dropoff for greedy
            //                                          nucleotide gapped extensions */
            // #define BLAST_GAP_X_DROPOFF_TBLASTX 0 /**< default dropoff for tblastx */
            //
            // /** default bit score that will trigger gapped extension */
            // #define BLAST_GAP_TRIGGER_PROT 22.0 /**< default bit score that will trigger
            //                                          a gapped extension for all protein-
            //                                          based searches */
            // #define BLAST_GAP_TRIGGER_NUCL 27.0  /**< default bit score that will trigger
            //                                          a gapped extension for blastn */
            //
            // /** default dropoff for the final gapped extension with traceback */
            // #define BLAST_GAP_X_DROPOFF_FINAL_PROT 25 /**< default dropoff (all protein-
            //                                                based gapped extensions) */
            // ```
            x_dropoff: 7.0,
            gap_x_dropoff: 15.0,
            gap_x_dropoff_final: 25.0,
            gap_trigger: 22.0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
            // ```c++
            // void
            // CBlastProteinOptionsHandle::SetEffectiveLengthsOptionsDefaults()
            // {
            //     SetDbLength(0);
            //     SetDbSeqNum(0);
            //     SetEffectiveSearchSpace(0);
            // }
            // ```
            // `-dbsize` (`CBlastOptions::SetDbLength`, blast_args.cpp:2570-2572).
            db_length: self.dbsize.unwrap_or(0),
            effective_search_space: self.searchsp.unwrap_or(0),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:922-923
            // ```c
            //     use_old_fsc = getenv("OLD_FSC");
            //     if (!use_old_fsc) sbp->gbp = s_BlastGumbelBlkNew();
            // ```
            // `-searchsp` sets OLD_FSC (blast_args.cpp:303-307, `parameters::gumbel_block`).
            old_fsc: self.searchsp.is_some() || crate::blastinput::app::old_fsc_environment(),
            window_size,
            evalue: self.evalue,
            sum_stats: self.sum_stats.unwrap_or(true),
            max_intron_length: self.max_intron_length,
            max_hsps: self.max_hsps.unwrap_or(0),
            query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
            culling_limit: self.culling_limit.unwrap_or(0),
            subject_besthit: self.subject_besthit,
            hitlist_size,
            num_descriptions: counts.num_descriptions as i32,
            num_alignments: counts.num_alignments as i32,
            outfmt,
            fields,
            tabular,
            taxdb_found: false,
            filtered_empty: false,
            num_threads: self.num_threads,
            query_range: None,
        })
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-220
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    //         _ASSERT(db_adapter && scope);
    //
    //         /*** Get the query sequence(s) ***/
    //         CRef<CQueryOptionsArgs> query_opts =
    //             m_CmdLineArgs->GetQueryOptionsArgs();
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    //         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                      query_opts->UseLowercaseMasks(),
    //                                      query_opts->GetParseDeflines(),
    //                                      query_opts->GetRange());
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    //
    //         /*** Get the formatting options ***/
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:276-290
    // ```c++
    // 	        BLAST_PROF_START( APP.LOOP.BLAST );
    //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
    //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
    //                 results = lcl_blast.Run();
    // 	        BLAST_PROF_STOP( APP.LOOP.BLAST );
    //             }
    // 	    BLAST_PROF_START( APP.LOOP.FMT );
    //             if (fmt_args->ArchiveFormatRequested(args)) {
    //                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
    //                 m_Bah.ResetMessages();
    //             } else {
    //                 BlastFormatter_PreFetchSequenceData(*results, scope,
    //                 		                            fmt_args->GetFormattedOutputChoice());
    //             	ITERATE(CSearchResultSet, result, *results) {
    //                	    formatter.PrintOneResultSet(**result, query_batch);
    // ```
    pub fn run(self) -> Result<()> {
        super::native::run(&self)
    }

    /// `-seg` as NCBI's filtering handler reads it (`app::read_seg_option`).
    pub(crate) fn seg_spec(
        &self,
    ) -> std::result::Result<SegSpec, crate::blastinput::app::SegOptionError> {
        crate::blastinput::app::read_seg_option(&self.seg, "BLASTX")
    }

    /// The `-db` value as NCBI keeps it (`seqdb::db_argument_name`), when given.
    pub fn database_name(&self) -> Option<&str> {
        self.db.as_deref().map(crate::seqdb::db_argument_name)
    }
}

#[cfg(test)]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3630-3642
// ```c++
//     CBlastOptions& opts = retval->SetOptions();
//     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
//         (*arg)->ExtractAlgorithmOptions(args, opts);
//     }
//
//     m_IsUngapped = !opts.GetGappedMode();
//     try { retval->Validate(); }
//     catch (const CBlastException& e) {
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
//     }
//     return retval;
// }
//
// ```
mod tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3630-3642
    // ```c++
    //     CBlastOptions& opts = retval->SetOptions();
    //     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    //         (*arg)->ExtractAlgorithmOptions(args, opts);
    //     }
    //
    //     m_IsUngapped = !opts.GetGappedMode();
    //     try { retval->Validate(); }
    //     catch (const CBlastException& e) {
    //         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    //     }
    //     return retval;
    // }
    //
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:288-300
    // ```c++
    //     if ( args.Exist(kArgWordSize) && args[kArgWordSize]) {
    //         if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
    //            opt.SetLookupTableType(eCompressedAaLookupTable);
    //            opt.SetWordThreshold(19.3);
    //            if (args[kArgWordSize].AsInteger() > 5) {
    //                opt.SetWordThreshold(21.0);
    //            }
    //            if (args[kArgWordSize].AsInteger() > 6) {
    //                opt.SetWordThreshold(20.25);
    //            }
    //         }
    //         opt.SetWordSize(args[kArgWordSize].AsInteger());
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
    // ```c++
    // CBlastOptionsLocal::SetWordSize(int ws)
    // {
    //     m_LutOpts->word_size = ws;
    //     if (m_LutOpts->lut_type == eCompressedAaLookupTable && ws <= 4)
    // 	m_LutOpts->lut_type = eAaLookupTable;
    //     else if (m_LutOpts->lut_type == eAaLookupTable && ws > 4)
    // 	m_LutOpts->lut_type = eCompressedAaLookupTable;
    // }
    // ```
    // Item SC9: every word size of NCBI's runs; one above 4 is the compressed table with its own
    // threshold, one of 4 or less the lookup table with the task's (blastx-fast's 20 too), and a
    // given -threshold stays. IDENTITY is refused with the compressed table, where NCBI crashes.
    #[test]
    fn every_word_size_takes_ncbi_s_table_and_threshold() {
        for (word, threshold, compressed) in [
            ("2", 12.0, false),
            ("3", 12.0, false),
            ("4", 12.0, false),
            ("5", 19.3, true),
            ("6", 21.0, true),
            ("7", 20.25, true),
        ] {
            let resolved = parse(&["-word_size", word]).resolve().unwrap();
            assert_eq!(resolved.word_size.to_string(), word);
            assert_eq!(resolved.threshold, threshold, "{word}");
            assert_eq!(resolved.compressed_lookup(), compressed, "{word}");
            let fast = parse(&["-task", "blastx-fast", "-word_size", word])
                .resolve()
                .unwrap();
            let fast_threshold = if compressed { threshold } else { 20.0 };
            assert_eq!(fast.threshold, fast_threshold, "blastx-fast {word}");
            assert_eq!(fast.compressed_lookup(), compressed, "blastx-fast {word}");
            let given = parse(&["-word_size", word, "-threshold", "17"])
                .resolve()
                .unwrap();
            assert_eq!(given.threshold, 17.0, "{word}");
            let identity = parse(&[
                "-word_size",
                word,
                "-matrix",
                "IDENTITY",
                "-ungapped",
                "-comp_based_stats",
                "0",
            ])
            .resolve();
            if compressed {
                assert_eq!(
                    identity.err().unwrap().to_string(),
                    format!("-matrix IDENTITY with the compressed lookup table (-word_size {word}), with which NCBI BLAST+'s blastx crashes, is not supported by LOSAT's BLASTX")
                );
            } else {
                assert!(identity.is_ok(), "{word}");
            }
        }
        // The task's default word sizes: 3 is the lookup table, blastx-fast's 5 the compressed one.
        assert!(!parse(&[]).resolve().unwrap().compressed_lookup());
        assert!(parse(&["-task", "blastx-fast"])
            .resolve()
            .unwrap()
            .compressed_lookup());
    }
    fn parse(options: &[&str]) -> BlastxArgs {
        let mut argv = vec!["losat", "blastx", "-query", "q", "-subject", "s"];
        argv.extend_from_slice(options);
        let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).unwrap();
        match cli.command {
            crate::cli::Commands::Blastx(a) => a,
            _ => unreachable!(),
        }
    }
    #[test]
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1300-1306
    // ```c++
    //     if (kPhiBlast)
    //         return 0;
    //
    //     if (program_number != eBlastTypeBlastn &&
    //         program_number != eBlastTypeMapping &&
    //         (!Blast_ProgramIsRpsBlast(program_number)) &&
    //         options->threshold <= 0)
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
    // ```c++
    // 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
    // 	{
    // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //          "expect value or cutoff score must be greater than zero");
    // 		return BLASTERR_OPTION_VALUE_INVALID;
    // 	}
    // ```
    fn cli_acceptance_precedes_core_validation() {
        assert_eq!(parse(&["-threshold", "0"]).threshold, Some(0.0));
        assert!(parse(&["-threshold", "0"]).resolve().is_err());
        assert_eq!(parse(&["-evalue", "0"]).evalue, 0.0);
        assert!(parse(&["-evalue", "0"]).resolve().is_err());
    }
    #[test]
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2910-2928
    // ```c++
    //     if(m_OutputFormat <= eFlatQueryAnchoredNoIdentities) {
    //
    //
    //     	 m_NumDescriptions = m_DfltNumDescriptions;
    //     	 m_NumAlignments = m_DfltNumAlignments;
    //
    //     	 if (args.Exist(kArgNumDescriptions) && args[kArgNumDescriptions]) {
    //     	    m_NumDescriptions = args[kArgNumDescriptions].AsInteger();
    //     	 }
    //
    //          if (args.Exist(kArgNumAlignments) && args[kArgNumAlignments]) {
    //     		m_NumAlignments = args[kArgNumAlignments].AsInteger();
    //     	}
    //
    //          if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //     	    m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //     		m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    //     		hitlist_size = m_NumAlignments;
    //          }
    // ```
    fn omission_and_explicit_target_count_remain_distinct() {
        let omitted = parse(&[]).resolve().unwrap();
        let explicit = parse(&["-max_target_seqs", "500"]).resolve().unwrap();
        assert_eq!(
            (
                omitted.hitlist_size,
                omitted.num_descriptions,
                omitted.num_alignments
            ),
            (500, 500, 250)
        );
        assert_eq!(
            (
                explicit.hitlist_size,
                explicit.num_descriptions,
                explicit.num_alignments
            ),
            (500, 500, 500)
        );
    }

    // NCBI blast_args.cpp:2894-2973 (`blastinput::report_counts`).
    #[test]
    fn num_descriptions_and_num_alignments_set_the_counts() {
        let counts = |extra: &[&str]| {
            let resolved = parse(extra).resolve().unwrap();
            (
                resolved.hitlist_size,
                resolved.num_descriptions,
                resolved.num_alignments,
            )
        };
        assert_eq!(counts(&["-num_descriptions", "3"]), (250, 3, 250));
        assert_eq!(counts(&["-num_alignments", "2"]), (500, 500, 2));
        assert_eq!(counts(&["-outfmt", "6", "-num_alignments", "2"]), (2, 2, 2));
        assert_eq!(
            counts(&["-outfmt", "6", "-num_descriptions", "2"]),
            (500, 500, 500)
        );
        let error = parse(&["-num_descriptions", "0", "-num_alignments", "0"])
            .resolve()
            .unwrap_err();
        assert_eq!(error.to_string(), "No hits are being saved");
    }

    // NCBI blast_args.cpp:2801-2855 and tabular.cpp:76-99: BLASTX reads -outfmt and its
    // fields as NCBI (a specification counts only for the tabular formats; unknown words are
    // dropped, `-name` removes a field); a list other than the `std` fields with a tab, and
    // outfmt 10 and 20, go to the shared writer (`tabular`).
    #[test]
    fn outfmt_and_fields_are_read_as_ncbi() {
        use crate::report::tabular::{TabularKind, STD_FIELDS};
        let fields = |outfmt: &str| {
            let resolved = parse(&["-outfmt", outfmt]).resolve().unwrap();
            (resolved.outfmt, resolved.fields)
        };
        let std: Vec<String> = STD_FIELDS.iter().map(|f| f.name().to_string()).collect();
        assert_eq!(fields("0 std"), (0, Vec::new()));
        assert_eq!(fields("0 foo"), (0, Vec::new()));
        assert_eq!(fields(" 6 delim="), (6, std.clone()));
        assert_eq!(fields("7"), (7, std.clone()));
        assert_eq!(fields("6 foo"), (6, std.clone()));
        // Words between single spaces; an empty word names no field.
        let mut listed = vec!["qseqid".to_string(), "sseqid".to_string()];
        listed.extend(std.iter().cloned());
        assert_eq!(fields("7 qseqid  sseqid std"), (7, listed));
        // The XInclude formats to standard output fail later (`native::run`).
        assert_eq!(fields("13").0, 13);
        let tabular = |outfmt: &str| {
            parse(&["-outfmt", outfmt])
                .resolve()
                .unwrap()
                .tabular
                .map(|spec| (spec.kind, spec.fields.len(), spec.delimiter.clone()))
        };
        for outfmt in ["6", "7 std", "6 foo", "6 -qseqid", "0 qseqid", "5"] {
            assert_eq!(tabular(outfmt), None, "{outfmt:?}");
        }
        assert_eq!(
            tabular("6 std -evalue staxids"),
            Some((TabularKind::Tabular, 12, "\t".to_string()))
        );
        assert_eq!(
            tabular("7 delim=| std"),
            Some((TabularKind::Commented, 12, "|".to_string()))
        );
        assert_eq!(tabular("10"), Some((TabularKind::Csv, 12, ",".to_string())));
        assert_eq!(
            tabular("20 qseqid"),
            Some((TabularKind::CsvWithHeader, 1, ",".to_string()))
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3636-3639
// ```c++
//     try { retval->Validate(); }
//     catch (const CBlastException& e) {
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
//     }
// ```
// The native app catches this category after stream extraction. Product scope
// refusals retain the separate registered CLI-v2 boundary.
#[derive(Debug)]
pub(crate) struct SemanticOptionsError(pub std::borrow::Cow<'static, str>);

impl From<&'static str> for SemanticOptionsError {
    fn from(message: &'static str) -> Self {
        Self(std::borrow::Cow::Borrowed(message))
    }
}

impl From<String> for SemanticOptionsError {
    fn from(message: String) -> Self {
        Self(std::borrow::Cow::Owned(message))
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
// ```c++
//     catch (const blast::CInputException& e) {                               \
//         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
//         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
//         exit_code = BLAST_INPUT_ERROR;                                      \
// ```
impl std::fmt::Display for SemanticOptionsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3638-3638
// ```c++
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
// ```
impl std::error::Error for SemanticOptionsError {}
