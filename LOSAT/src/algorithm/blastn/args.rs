use crate::blastinput::value_parsers::*;
use clap::Args;
use std::path::PathBuf;

/// The strand(s) of the query that a search covers (`-strand`): the strand of the query's
/// Seq-loc. A strand that is not searched has a context of length 0 and, in a split search, a
/// chunk context of `None` (`kInvalidContext`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:394-410
/// ```c
///                     if (ctx % NUM_STRANDS == 0) {
///                         if (kStrand == eNa_strand_minus) {
///                             m_SplitBlk->AddContextToChunk(chunk_num,
///                                                           kInvalidContext);
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QueryStrand {
    Both,
    Plus,
    Minus,
}

impl QueryStrand {
    /// The value of `-strand`: `both`, `plus` or `minus`, anything else is NCBI's argument
    /// error (the command line checks it at parsing, `run_local` with this function).
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
    /// ```c
    ///         arg_desc.AddDefaultKey(kArgStrand, "strand",
    ///                      "Query strand(s) to search against database/subject",
    ///                                CArgDescriptions::eString, kDfltArgStrand);
    /// ```
    pub fn parse(value: &str) -> anyhow::Result<Self> {
        match value {
            "both" => Ok(Self::Both),
            "plus" => Ok(Self::Plus),
            "minus" => Ok(Self::Minus),
            _ => anyhow::bail!(
                "Argument \"strand\". Illegal value, expected `both', `minus', `plus':  `{value}'"
            ),
        }
    }

    /// The value of `-strand` of arguments that have passed `parse` (clap's value parser, or
    /// `run_local`'s check): a value that is none of the three would be a bug of the caller.
    pub fn from_option(value: &str) -> Self {
        Self::parse(value).expect("-strand was checked when the arguments were read")
    }

    /// The first valid context of a query (1 for the minus strand).
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:146-149
    /// ```c
    ///     // Adjust the first context if the requested strand is the minus strand
    ///     if (strand == eNa_strand_minus) {
    ///         query_info->first_context = translate ? 3 : 1;
    ///     }
    /// ```
    pub fn first_context(self) -> usize {
        usize::from(self == Self::Minus)
    }

    /// Whether the plus (context 0) and minus (context 1) strand of a query is searched.
    pub fn searched(self) -> (bool, bool) {
        (self != Self::Minus, self != Self::Plus)
    }
}

// CLI v2 names/defaults: NCBI c++/src/algo/blast/blastinput/blast_args.cpp:
// 166-170,203-207,332-349,2657-2660,3158-3163; cmdline_flags.cpp:46-94.
// arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value", ..., eInteger);
// arg_desc.SetConstraint(kArgMaxHSPsPerSubject, new CArgAllowValuesGreaterThanOrEqual(1));
// arg_desc.AddDefaultKey(kArgNumThreads, "int_value", ..., NStr::IntToString(kDfltValue));
// The single-dash lexical translation is owned by crate::cli.
#[derive(Args, Debug)]
#[command(rename_all = "snake_case")]
pub struct BlastnArgs {
    // NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:47
    // ```c
    // const string kDfltArgQuery("-");
    // ```
    // `-` is standard input (`blastn/blast_engine/run.rs`).
    #[arg(long, default_value = "-", value_parser = blastn_input_path(), value_name = "PATH")]
    pub query: PathBuf,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
    // ```c
    //     } else if (!m_IsIgBlast){
    //         // IgBlast permits use of germline database
    //         NCBI_THROW(CInputException, eInvalidInput,
    //            "Either a BLAST database or subject sequence(s) must be specified");
    //     }
    // ```
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastn/blast_engine/run.rs`).
    #[arg(long, value_parser = blastn_input_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
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
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). The exclusion
    // is an error of the argument parser: LOSAT's parser reports it with its own text
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1). The database is searched by
    // `blastn/blast_engine/database.rs`.
    #[arg(
        long,
        value_name = "NAME",
        conflicts_with_all = ["subject", "subject_loc"],
        help = "BLAST database name (a nucleotide database without parsed Seq-ids)"
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
    // No constraint and no dependency: read with `-db` and with `-subject`
    // (`BlastnArgs::db_length`).
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
    #[arg(long, default_value = "megablast", help = "Task to execute: megablast, blastn, dc-megablast or blastn-short", long_help = "Implemented tasks: megablast, blastn, dc-megablast and blastn-short. Task defaults: megablast uses word size 28, reward 1, penalty -2, gaps 0/0; blastn uses word size 11, reward 2, penalty -3, gaps 5/2; dc-megablast uses the 11-of-18 coding template, reward 2, penalty -3, gaps 5/2 and two hits in a window of 40; blastn-short uses word size 7, reward 1, penalty -3, gaps 5/2, e-value 1000 and no DUST. An omitted option takes the default of the task.", value_parser = blastn_task)]
    pub task: String,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:708-730
    // ```c
    //     arg_desc.AddOptionalKey(kArgDMBTemplateType, "type",
    //                  "Discontiguous MegaBLAST template type",
    //                  CArgDescriptions::eString);
    //     ...
    //     arg_desc.SetDependency(kArgDMBTemplateType,
    //                            CArgDescriptions::eRequires,
    //                            kArgDMBTemplateLength);
    //
    //     arg_desc.AddOptionalKey(kArgDMBTemplateLength, "int_value",
    //                  "Discontiguous MegaBLAST template length",
    //                  CArgDescriptions::eInteger);
    //     ...
    //     arg_desc.SetDependency(kArgDMBTemplateLength,
    //                            CArgDescriptions::eRequires,
    //                            kArgDMBTemplateType);
    // ```
    // NCBI applies them to any task (`coordination.rs` `determine_template`).
    #[arg(long, value_name = "TYPE", value_parser = blastn_template_type, requires = "template_length", help = "Discontiguous MegaBLAST template type: coding, optimal or coding_and_optimal")]
    pub template_type: Option<crate::algorithm::blastn::disc_lookup::DiscWordType>,
    #[arg(long, value_name = "INT", value_parser = blastn_template_length, requires = "template_type", help = "Discontiguous MegaBLAST template length: 16, 18 or 21")]
    pub template_length: Option<u8>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:166-170,288-300
    // ```c
    //         arg_desc.AddOptionalKey(kArgWordSize, "int_value", description,
    //                                 CArgDescriptions::eInteger);
    // ...
    //     if ( args.Exist(kArgWordSize) && args[kArgWordSize]) {
    // ...
    //         opt.SetWordSize(args[kArgWordSize].AsInteger());
    // ```
    // An omitted value keeps the default of the task (`coordination.rs`).
    #[arg(long, value_parser = blastn_word_size, help = "Word size for wordfinder algorithm (default: 28 for megablast, 11 for blastn and dc-megablast, 7 for blastn-short)")]
    pub word_size: Option<usize>,
    #[arg(long, default_value_t = 1, value_parser = blastn_count)]
    pub num_threads: usize,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3354-3363
    // ```c++
    //     arg_desc.AddDefaultKey(
    //             kArgUseIndex, "boolean",
    //             "Use MegaBLAST database index",
    //             CArgDescriptions::eBoolean, NStr::BoolToString(kDfltArgUseIndex));
    //     arg_desc.AddOptionalKey(
    //             kArgIndexName, "string",
    //             "MegaBLAST database index name (deprecated; use only for old style indices)",
    //             CArgDescriptions::eString );
    // ```
    /// `-use_index false` (NCBI's default) is accepted and changes nothing; `true` is refused
    /// by the value parser: LOSAT has no MegaBLAST database index. `-index_name` stays unported.
    #[arg(long, value_name = "boolean", value_parser = crate::blastinput::value_parsers::blastn_use_index, help = "Use MegaBLAST database index (only false is supported) [default: false]")]
    pub use_index: Option<bool>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:142-146
    // ```c
    //     	string des = "Expectation value (E) threshold for saving hits. Default = 10";
    //     	if(m_IsBlastn) {
    //     		des += " (1000 for blastn-short)";
    //     	}
    //         arg_desc.AddOptionalKey(kArgEvalue, "evalue", des, CArgDescriptions::eDouble);
    // ```
    // An omitted value keeps the default of the task (`coordination.rs` `determine_evalue`).
    #[arg(long, value_parser = blastn_evalue, help = "Expectation value (E) threshold for saving hits (default: 10; 1000 for blastn-short)")]
    pub evalue: Option<f64>,
    /// Percent identity threshold for filtering HSPs (Blast_HSPTest).
    ///
    /// Default: 0.0 (disabled, matches NCBI BLAST default from calloc initialization).
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    /// ```c
    /// return ((hsp->num_ident * 100.0 <
    ///         align_length * hit_options->percent_identity) ||
    ///         align_length < hit_options->min_hit_length) ;
    /// ```
    #[arg(long = "perc_identity", value_parser = blastn_percentage)]
    pub percent_identity: Option<f64>,
    /// Minimum hit length for filtering HSPs (Blast_HSPTest).
    ///
    /// Default: 0 (disabled, matches NCBI BLAST default from calloc initialization).
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    /// ```c
    /// return ((hsp->num_ident * 100.0 <
    ///         align_length * hit_options->percent_identity) ||
    ///         align_length < hit_options->min_hit_length) ;
    /// ```
    // NCBI blastn has no option for it (the LOSAT-only option was removed, AGENTS.md
    // rule 5).
    #[arg(skip = 0usize)]
    pub min_hit_length: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960-2968
    // ```c
    // if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //    m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //    m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    //    hitlist_size = m_NumAlignments;
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:219,221
    // ```c
    // const size_t kDfltArgNumDescriptions = 500;
    // const size_t kDfltArgNumAlignments = 250;
    // ```
    // An omitted value keeps the default hit list size (500, `hitlist_size`) but the
    // pairwise report then shows 250 alignments, so the option has no clap default.
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep (default: 500)", conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    // NCBI blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// The hit list size (NCBI BLAST hitlist_size) and the numbers of descriptions and
    /// alignments of the pairwise report, which `process_options` sets from
    /// `-max_target_seqs`, `-num_descriptions`, `-num_alignments` and the output format
    /// (`blastinput::report_counts`). Not a public option.
    /// Reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
    #[arg(skip)]
    pub report_counts: crate::blastinput::report_counts::ReportCounts,
    /// Remove word seeds with high frequency in the searched database (off).
    /// Reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:257 (limit_lookup)
    // NCBI blastn has no such option: it is magicblast's (blast_args.cpp:1492-1501), so
    // the LOSAT option was removed (AGENTS.md rule 5).
    #[arg(skip = false)]
    pub limit_lookup: bool,
    /// Maximum database word count for lookup filtering.
    /// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:172-174
    /// #define MAX_DB_WORD_COUNT_MAPPER 30
    #[arg(skip = 30u8)]
    pub max_db_word_count: u8,
    /// Maximum number of HSPs per subject (unlimited when omitted)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    #[arg(long = "max_hsps", value_parser = blastn_count)]
    pub max_hsps_per_subject: Option<usize>,
    // NCBI blast_args.cpp:196-201,313-315: no filter when omitted or 0.
    #[arg(long = "qcov_hsp_perc", value_parser = blastn_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    /// Minimum diagonal separation between HSPs on the same subject (0 = auto, task-specific)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    /// Default: 50 for blastn, 6 for megablast
    // NCBI api/blast_nucl_options.cpp:239,259:
    // SetMinDiagSeparation(50); SetMinDiagSeparation(6);
    // Task resolution owns this internal value; no CLI override.
    #[arg(skip = 0usize)]
    pub min_diag_separation: usize,
    #[arg(long, value_name = "PATH", value_parser = blastn_output_path())]
    pub out: Option<PathBuf>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:648-660,674-679
    // ```c
    //     arg_desc.AddOptionalKey(kArgMismatch, "penalty",
    //                            "Penalty for a nucleotide mismatch",
    //                            CArgDescriptions::eInteger);
    //     arg_desc.SetConstraint(kArgMismatch,
    //                            new CArgAllowValuesLessThanOrEqual(0));
    // ...
    //     arg_desc.SetConstraint(kArgMatch,
    //                            new CArgAllowValuesGreaterThanOrEqual(0));
    // ...
    //     if (cmd_line_args.Exist(kArgMismatch) && cmd_line_args[kArgMismatch]) {
    //         options.SetMismatchPenalty(cmd_line_args[kArgMismatch].AsInteger());
    //     }
    //     if (cmd_line_args.Exist(kArgMatch) && cmd_line_args[kArgMatch]) {
    //         options.SetMatchReward(cmd_line_args[kArgMatch].AsInteger());
    //     }
    // ```
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:175-182,262-271
    // ```c
    //         arg_desc.AddOptionalKey(kArgGapOpen, "open_penalty",
    //                                 "Cost to open a gap",
    //                                 CArgDescriptions::eInteger);
    // ...
    //     if (args.Exist(kArgGapOpen) && args[kArgGapOpen]) {
    //         opt.SetGapOpeningCost(args[kArgGapOpen].AsInteger());
    //     }
    // ```
    // An omitted option keeps the default of the task (`coordination.rs`); a given one is
    // checked as NCBI checks it (`scoring.rs`). NCBI keeps the reward and the penalty in 16
    // bits; a reward of 0 or less (NCBI's rmblastn matrix scoring, or no valid query) is
    // not implemented and is rejected after NCBI's checks (`scoring.rs`).
    #[arg(long, value_parser = blastn_reward, help = "Reward for a nucleotide match (default: 1 for megablast and blastn-short, 2 for blastn and dc-megablast; LOSAT does not support 0)")]
    pub reward: Option<i32>,
    #[arg(long, value_parser = blastn_penalty, help = "Penalty for a nucleotide mismatch (default: -2 for megablast, -3 for the other tasks)")]
    pub penalty: Option<i32>,
    #[arg(
        long = "gapopen",
        value_parser = blastn_gap_cost,
        help = "Cost to open a gap (default: 0 for megablast, 5 for the other tasks)"
    )]
    pub gap_open: Option<i32>,
    #[arg(
        long = "gapextend",
        value_parser = blastn_gap_cost,
        help = "Cost to extend a gap (default: 0 for megablast, 2 for the other tasks)"
    )]
    pub gap_extend: Option<i32>,
    // NCBI blast_args.cpp:410-420: opt.SetDustFiltering(false/true);
    // opt.SetDustFilteringLevel(...); opt.SetDustFilteringWindow(...);
    // opt.SetDustFilteringLinker(...);
    // NCBI reads -dust in its filtering handler, after the output is opened
    // (`resolve_dust`); an omitted value keeps the default (20 64 1).
    #[arg(
        long = "dust",
        value_name = "DUST",
        help = "DUST: no, yes, or LEVEL WINDOW LINKER (default: 20 64 1; no for blastn-short)"
    )]
    pub dust_filtering: Option<String>,
    #[arg(skip = DustSpec::Yes)]
    pub dust: DustSpec,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:1939-1942
    // ```c
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //      "Use lower case filtering in query and subject sequence(s)?", true);
    // ```
    #[arg(long = "lcase_masking", default_value_t = false)]
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:358-363
    // ```c++
    //         arg_desc.AddOptionalKey(kArgWindowMaskerDatabase, "window_masker_db",
    //                 "Enable WindowMasker filtering using this repeats database.",
    //                 CArgDescriptions::eString);
    //         arg_desc.SetDependency(kArgWindowMaskerDatabase,
    //                                CArgDescriptions::eExcludes,
    //                                kArgRemote);
    // ```
    // Read by the filtering handler after `-dust` (`resolve_window_masker`).
    #[arg(
        long = "window_masker_db",
        value_name = "window_masker_db",
        help = "Enable WindowMasker filtering using this repeats database."
    )]
    pub window_masker_db: Option<String>,
    /// The masker of `-window_masker_db` (`resolve_window_masker`); `None` without the option.
    #[arg(skip)]
    pub window_masker: Option<std::sync::Arc<super::window_masker::WindowMasker>>,
    /// Where `-window_masker_db` is read from: the file system when `None`; LOSAT Web's adapter
    /// gives the counts files it holds under their names.
    #[arg(skip)]
    pub window_masker_counts: Option<std::sync::Arc<dyn crate::utils::winmask::CountsSource>>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath("taxdb.bti")`), asked
    /// by the command line before the search where a tabular field list shows a name of the
    /// taxonomy database (`algorithm::blastn::tabular`); `None` elsewhere (`run_local`'s callers
    /// give no search path: not found). Not a public option.
    #[arg(skip)]
    pub tabular_taxdb: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
    // ```c++
    //         arg_desc.AddDefaultKey(kArgStrand, "strand",
    //                      "Query strand(s) to search against database/subject",
    //                                CArgDescriptions::eString, kDfltArgStrand);
    // ```
    #[arg(long, default_value = "both", value_parser = ["both", "plus", "minus"])]
    pub strand: String,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // ```c++
    //     // query location
    //     arg_desc.AddOptionalKey(kArgQueryLocation, "range",
    //                             "Location on the query sequence in 1-based offsets "
    //                             "(Format: start-stop)",
    //                             CArgDescriptions::eString);
    // ```
    // Read by the query options handler (`search_cli`), after the filtering handler.
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // ```c++
    //         // subject location
    //         arg_desc.AddOptionalKey(kArgSubjectLocation, "range",
    //                         "Location on the subject sequence in 1-based offsets "
    //                         "(Format: start-stop)",
    //                         CArgDescriptions::eString);
    // ```
    // Read by the database arguments handler when it reads the subjects (`run`).
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        conflicts_with_all = crate::blastinput::db_filter::DB_FILTER_OPTIONS,
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    /// Apply subject best hit filtering (disabled by default).
    /// Reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:135
    #[arg(long = "subject_besthit", default_value_t = false)]
    pub subject_besthit: bool,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340
    // ```c
    // CHspFilteringArgs::ExtractAlgorithmOptions(const CArgs& args,
    //                                            CBlastOptions& opts)
    // {
    //     if (args[kArgCullingLimit]) {
    //         opts.SetCullingLimit(args[kArgCullingLimit].AsInteger());
    //     }
    // ```
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    #[arg(long, value_parser = nonnegative_ncbi_integer)]
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917,929-929
    // ```c++
    //     arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
    // ...
    //     opt.SetGappedMode(!args[kArgUngapped]);
    // ```
    /// Perform ungapped alignment only.
    #[arg(long)]
    pub ungapped: bool,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:212-229,276-286
    // ```c++
    //         arg_desc.AddOptionalKey(kArgUngappedXDropoff, "float_value",
    //                             "X-dropoff value (in bits) for ungapped extensions",
    //                             CArgDescriptions::eDouble);
    // ...
    //             arg_desc.AddOptionalKey(kArgGappedXDropoff, "float_value",
    //                  "X-dropoff value (in bits) for preliminary gapped extensions",
    //                  CArgDescriptions::eDouble);
    // ...
    //             arg_desc.AddOptionalKey(kArgFinalGappedXDropoff, "float_value",
    //                          "X-dropoff value (in bits) for final gapped alignment",
    //                          CArgDescriptions::eDouble);
    // ...
    //     if (args.Exist(kArgUngappedXDropoff) && args[kArgUngappedXDropoff]) {
    //         opt.SetXDropoff(args[kArgUngappedXDropoff].AsDouble());
    //     }
    //
    //     if (args.Exist(kArgGappedXDropoff) && args[kArgGappedXDropoff]) {
    //         opt.SetGapXDropoff(args[kArgGappedXDropoff].AsDouble());
    //     }
    //
    //     if (args.Exist(kArgFinalGappedXDropoff) && args[kArgFinalGappedXDropoff]) {
    //         opt.SetGapXDropoffFinal(args[kArgFinalGappedXDropoff].AsDouble());
    //     }
    // ```
    // An omitted value keeps the default of the task (`coordination.rs` `task_defaults`: the
    // nucleotide handles of api/blast_nucl_options.cpp:163-193 and disc_nucl_options.cpp:66-84).
    // The ungapped X-drop of megablast is 0 (its handle sets the window only, and the options
    // are calloc'd, api/blast_options_local_priv.cpp:54): the raw X-drop is then the word cutoff
    // (blast_parameters.c:378-381). The values are bits; NCBI converts them to raw scores once,
    // with the smallest Lambda of the batch (`scoring.rs` `gap_x_dropoffs`, `run.rs`
    // `x_dropoff_init`).
    #[arg(long, value_parser = blastn_real, help = "X-dropoff value (in bits) for ungapped extensions (default: 20 for blastn, dc-megablast and blastn-short; megablast uses the ungapped cutoff score)")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long, value_parser = blastn_real, help = "X-dropoff value (in bits) for preliminary gapped extensions (default: 25 for megablast, 30 for the other tasks)")]
    pub xdrop_gap: Option<f64>,
    #[arg(long, value_parser = blastn_real, help = "X-dropoff value (in bits) for final gapped alignment (default: 100)")]
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:663-665,681-685
    // ```c++
    //     arg_desc.AddFlag(kArgNoGreedyExtension,
    //                      "Use non-greedy dynamic programming extension",
    //                      true);
    // ...
    //     if (cmd_line_args.Exist(kArgNoGreedyExtension) &&
    //         cmd_line_args[kArgNoGreedyExtension]) {
    //         options.SetGapExtnAlgorithm(eDynProgScoreOnly);
    //         options.SetGapTracebackAlgorithm(eDynProgTbck);
    //     }
    // ```
    // The X-drops keep the task's defaults (megablast's greedy 25 bits). With the gap costs
    // 0/0 of megablast NCBI refuses the search (`scoring.rs` `check_scoring_options`).
    /// Use non-greedy dynamic programming extension.
    #[arg(long)]
    pub no_greedy: bool,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:701-704,739-741
    // ```c++
    //     arg_desc.AddOptionalKey(kArgMinRawGappedScore, "int_value",
    //                             "Minimum raw gapped score to keep an alignment "
    //                             "in the preliminary gapped and traceback stages",
    //                             CArgDescriptions::eInteger);
    // ...
    //     if (args[kArgMinRawGappedScore]) {
    //         options.SetCutoffScore(args[kArgMinRawGappedScore].AsInteger());
    //     }
    // ```
    // No constraint on the value: 0 and below keep the cutoff of the e-value
    // (`coordination.rs` `configure_task`).
    #[arg(long, value_parser = ncbi_integer, help = "Minimum raw gapped score to keep an alignment in the preliminary gapped and traceback stages")]
    pub min_raw_gapped_score: Option<i32>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:475-481,487-497
    // ```c++
    //     arg_desc.AddOptionalKey(kArgWindowSize, "int_value",
    //                             "Multiple hits window size, use 0 to specify "
    //                             "1-hit algorithm",
    //                             CArgDescriptions::eInteger);
    //     arg_desc.SetConstraint(kArgWindowSize,
    //                            new CArgAllowValuesGreaterThanOrEqual(0));
    // ...
    //     if (args[kArgWindowSize]) {
    //         opt.SetWindowSize(args[kArgWindowSize].AsInteger());
    //     } else {
    //         int window = -1;
    //         BLAST_GetSuggestedWindowSize(opt.GetProgramType(),
    //                                      opt.GetMatrixName(),
    //                                      &window);
    // ```
    // BLASTN has no matrix (no suggested window), so an omitted value keeps the default of the
    // task (`coordination.rs` `determine_window_size`).
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm (default: 40 for dc-megablast, 0 for the other tasks)")]
    pub window_size: Option<i32>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:505-511,518-522
    // ```c++
    //     arg_desc.AddDefaultKey(kArgOffDiagonalRange, "int_value",
    //                             "Number of off-diagonals to search for the 2nd hit, "
    //                             "use 0 to turn off",
    //                             CArgDescriptions::eInteger,
    //                             NStr::IntToString(kDfltOffDiagonalRange));
    //     arg_desc.SetConstraint(kArgOffDiagonalRange,
    //                            new CArgAllowValuesGreaterThanOrEqual(0));
    // ...
    //     if (args[kArgOffDiagonalRange]) {
    //         opt.SetOffDiagonalRange(args[kArgOffDiagonalRange].AsInteger());
    //     } else {
    //         opt.SetOffDiagonalRange(0);
    //     }
    // ```
    // NCBI refuses a range without a window (`scoring.rs` `check_scoring_options`).
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Number of off-diagonals to search for the 2nd hit, use 0 to turn off (default: 0)")]
    pub off_diagonal_range: Option<i32>,
    // NCBI blastn has no such option (the LOSAT-only option was removed, AGENTS.md rule 5).
    #[arg(skip = false)]
    pub verbose: bool,
    /// Scan stride for subject sequence scanning (NCBI BLAST optimization).
    /// Higher values skip more positions, reducing k-mer lookups but potentially missing some seeds.
    /// Default: 0 (auto-calculate based on word_size: 1 for word_size < 16, 4 for word_size >= 16).
    /// For megablast (word_size=28), scan_step=4 reduces lookups by ~4x with minimal sensitivity loss.
    // NCBI core/blast_nalookup.c:397,1334:
    // lookup->scan_step = lookup->word_length - lookup->lut_word_length + 1;
    // Lookup construction owns this internal value; no CLI override.
    #[arg(skip = 0usize)]
    pub scan_step: usize,

    /// Output format (NCBI BLAST compatible).
    ///
    /// Supported formats:
    ///   0 = Pairwise (the NCBI report; the default)
    ///   6 = Tabular (tab-separated values)
    ///   7 = Tabular with comment lines (headers)
    ///
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-782
    // ```c
    // case 6: case 7:
    //     formatter = new CBlastTabularFormatter(...);
    //     break;
    // ```
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
    // ```c
    // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    //     x_PrintField(*iter);
    // }
    // ```
    /// Custom field specifications are not yet ported and fail explicitly.
    ///
    /// Default fields: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
    // NCBI's argument is a string (blast_args.cpp:2657-2660), parsed when the options
    // are set (`blastinput/app.rs` `parse_formatting_string` and `report_format`).
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,

    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    #[arg(skip)]
    pub archive: Option<crate::report::archive::ArchiveHandle>,

    /// The command line of the run as the SAM `@PG` line shows it: every argument followed by
    /// one space, raw bytes (`GetCmdlineArgs`); `None` where there is no command line (the
    /// library, ABI v1). Set by `run` from the process's arguments.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:878-885
    /// ```c++
    /// string
    /// GetCmdlineArgs(const CNcbiArguments & a)
    /// {
    /// 	string cmd = kEmptyStr;
    /// 	for(unsigned int i=0; i < a.Size(); i++) {
    /// 		cmd += a[i] + " ";
    /// 	}
    /// 	return cmd;
    /// }
    /// ```
    #[arg(skip)]
    pub command_line: Option<Vec<u8>>,
}

impl BlastnArgs {
    /// The `-subject` path (empty before NCBI's check that it is given).
    pub fn subject_path(&self) -> &std::path::Path {
        self.subject.as_deref().unwrap_or(std::path::Path::new(""))
    }

    /// The `-db` value as NCBI keeps it: without the white space at its ends; `None` for a
    /// search of `-subject`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    /// ```c++
    ///     if (args.Exist(kArgDb) && args[kArgDb]) {
    /// 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
    ///
    ///         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
    ///                                              mol_type));
    /// ```
    pub fn database_name(&self) -> Option<&str> {
        self.db.as_deref().map(crate::seqdb::db_argument_name)
    }

    /// The database length of the options: the `-dbsize` value, 0 without it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2570-2572
    /// ```c++
    ///     if (args.Exist(kArgDbSize) && args[kArgDbSize]) {
    ///         opts.SetDbLength(args[kArgDbSize].AsInt8());
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:275
    /// ```c++
    ///     SetDbLength(0);
    /// ```
    /// `-perc_identity`, 0 (disabled) when it is not given. The BLAST archive records whether
    /// it was given (`percent_identity`).
    pub fn percent_identity_value(&self) -> f64 {
        self.percent_identity.unwrap_or(0.0)
    }

    /// `-culling_limit`, 0 (no culling) when it is not given. The BLAST archive records
    /// whether it was given (`culling_limit`).
    pub fn culling_limit_value(&self) -> i32 {
        self.culling_limit.unwrap_or(0)
    }

    /// `-searchsp`, 0 without it (`CBlastOptions::SetEffectiveSearchSpace`).
    pub fn options_searchsp(&self) -> i64 {
        self.searchsp.unwrap_or(0)
    }

    pub fn db_length(&self) -> i64 {
        self.dbsize.unwrap_or(0)
    }

    /// Reads the `-dust` value (`parse_dust_filtering`) where NCBI's filtering handler reads
    /// it, with NCBI's error.
    ///
    /// NCBI reference: c++/src/app/blast/blast_app_util.hpp:172-175
    /// ```c
    ///     catch (const blast::CInputException& e) {                               \
    ///         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
    ///         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
    ///         exit_code = BLAST_INPUT_ERROR;                                      \
    /// ```
    ///
    /// Without `-dust`, blastn-short does not filter: its task clears the filtering
    /// options (`coordination.rs` `task_dust_by_default`), and the argument has no default.
    ///
    /// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:344-349,410
    /// ```c
    ///         arg_desc.AddOptionalKey(kArgDustFiltering, "DUST_options",
    ///                         "Filter query sequence with DUST "
    ///                         "(Format: '" + kDfltArgApplyFiltering + "', " +
    ///                         "'level window linker', or '" + kDfltArgNoFiltering +
    ///                         "' to disable) Default = '20 64 1' ('" + kDfltArgNoFiltering + "' for blastn-short)",
    ///                         CArgDescriptions::eString);
    /// ...
    ///         if ( !m_QueryIsProtein && args[kArgDustFiltering]) {
    /// ```
    /// The filtering handler's `-window_masker_db`, after `-dust`: NCBI's check of the counts
    /// file (`window_masker::open`), and the masker made from it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-458
    /// ```c++
    ///     if (args.Exist(kArgWindowMaskerDatabase) &&
    ///         args[kArgWindowMaskerDatabase]) {
    ///         const string& stat_file = args[kArgWindowMaskerDatabase].AsString();
    ///         const CSeqMaskerIstatFactory::EStatType type =
    ///             CSeqMaskerIstatFactory::DiscoverStatType(stat_file);
    /// ```
    pub fn resolve_window_masker(&mut self) -> anyhow::Result<()> {
        if let Some(name) = &self.window_masker_db {
            let source = super::window_masker::source_of(self.window_masker_counts.as_ref());
            self.window_masker = Some(std::sync::Arc::new(super::window_masker::open(
                name, &*source,
            )?));
        }
        Ok(())
    }

    pub fn resolve_dust(&mut self) -> anyhow::Result<()> {
        if self.dust_filtering.is_none()
            && !crate::algorithm::blastn::coordination::task_dust_by_default(&self.task)
        {
            self.dust = DustSpec::No;
        }
        if let Some(value) = self.dust_filtering.take() {
            self.dust = parse_dust_filtering(&value).map_err(|message| crate::cli::NativeError {
                exit: 1,
                message: format!(
                    "BLAST query/options error: {message}\nPlease refer to the BLAST+ user manual.\n"
                ),
            })?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::QueryStrand;

    #[test]
    fn the_strand_option_is_read_once() {
        assert_eq!(QueryStrand::from_option("both").searched(), (true, true));
        assert_eq!(QueryStrand::from_option("plus").searched(), (true, false));
        assert_eq!(QueryStrand::from_option("minus").searched(), (false, true));
    }

    #[test]
    fn an_unknown_strand_is_an_argument_error() {
        for value in ["foo", "", "Plus", "PLUS", "both "] {
            let error = QueryStrand::parse(value).unwrap_err().to_string();
            assert_eq!(
                error,
                format!(
                    "Argument \"strand\". Illegal value, expected `both', `minus', `plus':  `{value}'"
                )
            );
        }
    }
}
