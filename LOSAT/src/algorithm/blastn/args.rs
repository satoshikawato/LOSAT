use super::disc_lookup::DiscWordType;
use crate::blastinput::value_parsers::*;
use crate::wording::QUERY_OPTIONS_ERROR;
use clap::Args;
use std::path::PathBuf;

/// BLASTN's `-task` (NCBI's nucleotide-nucleotide tasks, vecscreen aside: a program of its own),
/// for the check of `run_local` (`TaskParser::check`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:211-222
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:57-59
pub(crate) const BLASTN_TASK: TaskParser = task_parser("BLASTN", "blastn");

/// The strand(s) of the query that a search covers (`-strand`): the strand of the query's
/// Seq-loc. A strand that is not searched has a context of length 0 and, in a split search, a
/// chunk context of `None` (`kInvalidContext`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:394-410
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:146-149
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
// Constrains `kArgMaxHSPsPerSubject` to integers >= 1 via `CArgAllowValuesGreaterThanOrEqual`,
// and registers `kArgNumThreads` with its default integer value.
// The single-dash lexical translation is owned by crate::cli.
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct BlastnArgs {
    /// The arguments that the command line gave, for the import of a search strategy
    /// (`cli::GivenArguments`; read only with `-import_search_strategy`). The first field: the
    /// parser takes each later field's values out of the matches as it reads them
    /// (`from_arg_matches_mut`), so this one reads the matches whole.
    #[command(flatten)]
    pub(crate) given: crate::cli::GivenArguments<BlastnArgs>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:47
    // `-` is standard input (`blastn/blast_engine/run.rs`).
    #[arg(long, default_value = "-", value_parser = ncbi_path(), value_name = "PATH")]
    pub query: PathBuf,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastn/blast_engine/run.rs`).
    #[arg(long, value_parser = ncbi_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). NCBI's argument reader reports the
    // exclusion (`cli::reader`). The database is searched by `blastn/blast_engine/database.rs`.
    #[arg(
        long,
        value_name = "NAME",
        help = "BLAST database name (a nucleotide database, with or without parsed Seq-ids)"
    )]
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    #[command(flatten)]
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    /// `-db_soft_mask` and `-db_hard_mask` (`blastinput::db_mask`; series NP, M4a-2).
    #[command(flatten)]
    pub db_mask: crate::blastinput::db_mask::DbMaskArgs,
    /// The database masking of the search as the formatter resolved it
    /// (`DbMaskOption::resolve`, in `run::search`); `None` before, and for a search without it.
    #[arg(skip)]
    pub db_mask_resolved: Option<crate::blastinput::db_mask::ResolvedMask>,
    /// The MegaBLAST database index that the search uses (`-use_index true` that NCBI's loader
    /// kept, `blast_engine::mb_index`, set in `run::search`); `None` before, and for a search
    /// without one.
    #[arg(skip)]
    pub(crate) megablast_index_loaded: Option<std::sync::Arc<crate::dbindex::search::LoadedIndex>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // No constraint and no dependency: read with `-db` and with `-subject`
    // (`BlastnArgs::db_length`).
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS", help = "Effective length of the database")]
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663
    #[arg(
        long = "show_gis",
        default_value_t = false,
        help = "Show NCBI GIs in deflines?"
    )]
    pub show_gis: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2692-2693
    // NCBI's `CBlastFormat::m_IsHTML` (blast_format.cpp:90,224); with `-outfmt` 1 to 4 LOSAT does
    // not write it yet (`blastinput::app::html_format_check`).
    #[arg(long = "html", default_value_t = false, help = "Produce HTML output?")]
    pub html: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2050-2052
    // `-parse_deflines` is a flag of the programs (blast_args.cpp:1961).
    #[arg(
        long = "parse_deflines",
        default_value_t = false,
        help = "Should the query and subject defline(s) be parsed?"
    )]
    pub parse_deflines: bool,
    /// The ids of the search's records under `-parse_deflines` (`run`).
    #[arg(skip)]
    pub id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    #[arg(long, value_parser = ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // The help texts name the tasks (the Web form shows the short one); `describe` lists NCBI's
    // set (`ncbi_app::args::choices`, in the order of NCBI's BLASTN page).
    #[arg(
        long,
        default_value = "megablast",
        help = "Task to execute: megablast, blastn, dc-megablast, blastn-short or rmblastn",
        long_help = "Implemented tasks: megablast, blastn, dc-megablast, blastn-short and rmblastn. Task defaults: megablast uses word size 28, reward 1, penalty -2, gaps 0/0; blastn uses word size 11, reward 2, penalty -3, gaps 5/2; dc-megablast uses the 11-of-18 coding template, reward 2, penalty -3, gaps 5/2 and two hits in a window of 40; blastn-short uses word size 7, reward 1, penalty -3, gaps 5/2, e-value 1000 and no DUST; rmblastn is blastn with reward 0 and penalty 0 (matrix-only scoring: X-drops as raw scores, a preliminary cutoff by raw score, and a pairwise report with raw scores only). An omitted option takes the default of the task."
    )]
    pub task: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:708-730
    // NCBI applies them to any task (`coordination.rs` `determine_template`).
    #[arg(long, value_name = "TYPE", value_parser = blastn_template_type, help = "Discontiguous MegaBLAST template type: coding, optimal or coding_and_optimal")]
    pub template_type: Option<DiscWordType>,
    #[arg(long, value_name = "INT", value_parser = ncbi_integer_as::<u8>, help = "Discontiguous MegaBLAST template length: 16, 18 or 21")]
    pub template_length: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-170,288-300
    // An omitted value keeps the default of the task (`coordination.rs`).
    #[arg(long, value_parser = ncbi_integer_as::<usize>, help = "Word size for wordfinder algorithm (default: 28 for megablast, 11 for blastn and dc-megablast, 7 for blastn-short)")]
    pub word_size: Option<usize>,
    #[arg(long, default_value_t = 1, value_parser = ncbi_integer_as::<usize>)]
    pub num_threads: usize,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3168-3181
    // `-mt_mode` is read and has no effect: LOSAT's output is the same at every thread count
    // and for every value, and equals NCBI's `-num_threads 1` output (Owner 2026-10-13;
    // CLI DECISIONS 5). NCBI's `-mt_mode 1` (and `0` when it chooses to split by queries)
    // runs the queries in batches and writes the end of an outfmt 0 report or each XML
    // document per batch (`blastn_app.cpp:134-170`, `blastx_app.cpp` `x_RunMTBySplitQuery`);
    // that framing is not reproduced. The value must be 0 to 2 and `-num_threads` must be
    // given with it, as in NCBI (clap's `requires` ignores a default value).
    #[arg(long, value_name = "int_value", value_parser = ncbi_integer_as::<u8>, help = "Multi-thread mode: 0 auto, 1 split by queries, 2 split by database (accepted; the output does not depend on it)")]
    pub mt_mode: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3354-3363
    /// `-use_index` as NCBI's handler reads it (`blast_engine::mb_index`): ignored with
    /// `-subject` and for a task other than megablast; with `-db` and megablast, NCBI loads the
    /// database's index. When it loads none, or the word size is below the index's, NCBI warns
    /// and searches without it; otherwise it searches with the index. LOSAT does both; with
    /// `-db_hard_mask` its search reads the masked letters as matching nothing (approved exception
    /// 24 of PD-LOSAT-NCBI-DEFECTS, `dbindex::search`), and it refuses the later searches of an
    /// index of only some database volumes where NCBI frees its list of initial hits
    /// (`blast_engine::mb_index_seeds`).
    #[arg(long, value_name = "boolean", value_parser = ncbi_boolean, help = "Use MegaBLAST database index [default: false]")]
    pub use_index: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3359-3362
    /// `-index_name`: the old-style index of `-use_index true` (`megablast_index_name`); without
    /// `-use_index true` NCBI ignores it.
    #[arg(
        long,
        value_name = "string",
        help = "MegaBLAST database index name (deprecated; use only for old style indices)"
    )]
    pub index_name: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:142-146
    // An omitted value keeps the default of the task (`coordination.rs` `determine_evalue`).
    #[arg(long, value_parser = ncbi_double, help = "Expectation value (E) threshold for saving hits (default: 10; 1000 for blastn-short)")]
    pub evalue: Option<f64>,
    /// Percent identity threshold for filtering HSPs (default: 0, no filtering)
    // Blast_HSPTest. The default 0.0 disables it, as NCBI's calloc'ed option does.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    #[arg(long = "perc_identity", value_parser = ncbi_double)]
    pub percent_identity: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2924-2927,2960-2968
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:219,221
    // An omitted value keeps the default hit list size (500, `hitlist_size`) but the
    // pairwise report then shows 250 alignments, so the option has no clap default.
    #[arg(long, value_parser = ncbi_integer_as::<usize>, help = "Maximum number of aligned sequences to keep (default: 500)")]
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = ncbi_integer_as::<usize>, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = ncbi_integer_as::<usize>, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    #[command(flatten)]
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    // The hit list size (NCBI BLAST hitlist_size) and the numbers of descriptions and
    // alignments of the pairwise report, which `process_options` sets from
    // `-max_target_seqs`, `-num_descriptions`, `-num_alignments` and the output format
    // (`blastinput::report_counts`). Not a public option.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
    #[arg(skip)]
    pub report_counts: crate::blastinput::report_counts::ReportCounts,
    // The `-html` links of the run, which `process_options` sets
    // (`blastinput::app::html_links`). Not a public option.
    #[arg(skip)]
    pub html_links: Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
    // Remove word seeds with high frequency in the searched database (off).
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:257 (limit_lookup)
    // NCBI blastn has no such option: it is magicblast's (c++/src/algo/blast/blastinput/blast_args.cpp:1492-1501), so
    // the LOSAT option was removed (AGENTS.md rule 5).
    #[arg(skip = false)]
    pub limit_lookup: bool,
    // Maximum database word count for lookup filtering.
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:172-174
    #[arg(skip = 30u8)]
    pub max_db_word_count: u8,
    /// Maximum number of HSPs per subject (unlimited when omitted)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    #[arg(long = "max_hsps", value_parser = ncbi_integer_as::<usize>)]
    pub max_hsps_per_subject: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201,313-315: no filter when omitted or 0.
    #[arg(long = "qcov_hsp_perc", value_parser = ncbi_double, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // Minimum diagonal separation between HSPs on the same subject (0 = auto, task-specific)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    // Default: 50 for blastn, 6 for megablast
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:239,259:
    // Task resolution owns this internal value; no CLI override.
    #[arg(skip = 0usize)]
    pub min_diag_separation: usize,
    #[arg(long, value_name = "PATH", value_parser = ncbi_path())]
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    // The search strategy file (`report::archive::strategy`), created by `run` before the other
    // options are read.
    #[arg(long, value_name = "filename", value_parser = ncbi_path(), help = "File name to record the search strategy used")]
    pub export_search_strategy: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3527-3538
    // The search strategy to import (`report::archive::read::import`), read by `run` before the
    // other options (`import::apply`). NCBI's argument reader reports its exclusion of
    // `-export_search_strategy` (`cli::reader`).
    #[arg(
        long,
        value_name = "filename",
        value_parser = ncbi_path(),
        help = "Search strategy to use"
    )]
    pub import_search_strategy: Option<PathBuf>,
    /// The run state of an imported search strategy (`import::apply`).
    #[arg(skip)]
    pub(crate) import: Option<crate::report::archive::read::import::ImportRun>,
    /// A search strategy's gapped extension and traceback algorithms (`GapExtnAlgorithm`,
    /// `GapTracebackAlgorithm`, the task's where one is missing) other than the pairs that LOSAT
    /// runs (dynamic programming, and megablast's greedy pair): NCBI's options check reads them
    /// where it comes (`scoring::check_scoring_options`), and LOSAT refuses them after NCBI's
    /// checks. `None` from the command line.
    #[arg(skip)]
    pub(crate) strategy_algorithms: Option<(i32, i32)>,
    /// A search strategy's discontiguous template type above 2 (`MBTemplateType`, NCBI's
    /// `Uint1`), which NCBI's options check refuses (`scoring::check_scoring_options`).
    #[arg(skip)]
    pub(crate) strategy_template_type: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:648-660,674-679
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:175-182,262-271
    // An omitted option keeps the default of the task (`coordination.rs`); a given one is
    // checked as NCBI checks it (`scoring.rs`). NCBI keeps the reward and the penalty in 16
    // bits; reward 0 with penalty 0 is NCBI's matrix-only scoring (`-task rmblastn`), a reward
    // of 0 or less with another penalty (no valid query) is rejected after NCBI's checks
    // (`scoring.rs`).
    #[arg(long, value_parser = ncbi_integer, help = "Reward for a nucleotide match (default: 1 for megablast and blastn-short, 2 for blastn and dc-megablast, 0 for rmblastn; 0 needs penalty 0)")]
    pub reward: Option<i32>,
    #[arg(long, value_parser = ncbi_integer, help = "Penalty for a nucleotide mismatch (default: -2 for megablast, 0 for rmblastn, -3 for the other tasks)")]
    pub penalty: Option<i32>,
    #[arg(
        long = "gapopen",
        value_parser = ncbi_integer,
        help = "Cost to open a gap (default: 0 for megablast, 5 for the other tasks)"
    )]
    pub gap_open: Option<i32>,
    #[arg(
        long = "gapextend",
        value_parser = ncbi_integer,
        help = "Cost to extend a gap (default: 0 for megablast, 2 for the other tasks)"
    )]
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:410-420
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1939-1942
    #[arg(long = "lcase_masking", default_value_t = false)]
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:365-368
    // Read by `soft_masking`.
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: true]")]
    pub soft_masking: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:350-352
    // Read by the filtering handler after `-dust` and before the WindowMasker options, which
    // count it (`resolve_window_masker`); the repeat search runs where the masks of each query
    // batch are made (`repeat_filter`).
    #[arg(
        long = "filtering_db",
        value_name = "filtering_database",
        help = "BLAST database containing filtering elements (i.e.: repeats)"
    )]
    pub filtering_db: Option<String>,
    /// The repeat search of `-filtering_db` (`repeat_filter::RepeatFilter`): the database's name
    /// and the settings that look it up, given by the caller that has them (`resolve_repeat_filter`);
    /// `None` without the option.
    #[arg(skip)]
    pub repeat_filter: Option<std::sync::Arc<super::repeat_filter::RepeatFilter>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:354-356
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2285
    // Read by the filtering handler after `-dust` (`resolve_window_masker`); 0 is no
    // WindowMasker, but it is one of the filtering databases given. The database arguments
    // exclude it with `-no_taxid_expansion` (`DbFilterArgs`, shared with programs that have no
    // `-window_masker_taxid`): an argument error, also with `-dryrun` (NCBI BLAST+ 2.17.0, line RG
    // item E1b probe `blastn.x.no_taxid_expansion.window_masker_taxid`: USAGE, rc 1).
    #[arg(
        long = "window_masker_taxid",
        value_name = "window_masker_taxid",
        value_parser = ncbi_integer,
        help = "Enable WindowMasker filtering using a Taxonomic ID"
    )]
    pub window_masker_taxid: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:358-363
    // Read by the filtering handler after `-dust` (`resolve_window_masker`).
    #[arg(
        long = "window_masker_db",
        value_name = "window_masker_db",
        help = "Enable WindowMasker filtering using this repeats database."
    )]
    pub window_masker_db: Option<String>,
    /// The masker of `-window_masker_db` or of `-window_masker_taxid` (`resolve_window_masker`);
    /// `None` without them (or with the taxid 0).
    #[arg(skip)]
    pub window_masker: Option<std::sync::Arc<super::window_masker::WindowMasker>>,
    /// Where `-window_masker_db` and the file of `-window_masker_taxid` are read from: the file
    /// system when `None`; LOSAT Web's adapter gives the counts files it holds under their names.
    #[arg(skip)]
    pub window_masker_counts: Option<std::sync::Arc<dyn crate::utils::winmask::CountsSource>>,
    /// The registry of `-window_masker_taxid` as the caller read it (the directory of its files
    /// and what each load writes, `ncbi_environment::ApplicationSettings::window_masker_registry`),
    /// or LOSAT Web's registry of the counts files it holds (`window_masker_counts`: its directory,
    /// no posts); `None` when the caller gives neither, which `-window_masker_taxid` other than 0
    /// then refuses (`run_local_formats`).
    #[arg(skip)]
    pub window_masker_registry:
        Option<Result<crate::blastinput::ncbi_environment::WindowMaskerRegistry, String>>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath` of `taxdb.bti`), asked
    /// before the search where a tabular field list looks every taxid up (the plural names of the
    /// taxonomy database, the organism report; `algorithm::blastn::tabular`), by the command line
    /// and the library (`common::database::subject_taxonomy`); `None` elsewhere (LOSAT Web has no
    /// files: not found). Not a public option.
    #[arg(skip)]
    pub tabular_taxdb: Option<bool>,
    /// The taxonomy database of a `-subject` search's organism report (outfmt 18), made with
    /// `tabular_taxdb` (`common::database::SubjectTaxdb`); `None` elsewhere. Not a public option.
    #[arg(skip)]
    pub subject_taxdb: Option<crate::algorithm::common::database::SubjectTaxdb>,
    /// The `[BLAST] BLASTDB` setting through which a `-subject` search's tabular reports look for
    /// the taxonomy database (`common::database::subject_taxonomy`), which a search of records of
    /// the library (`api::run`) reads from the call's application settings when a report looks
    /// every taxid up; `None` elsewhere (the command line passes its settings, LOSAT Web has no
    /// files). Not a public option.
    #[arg(skip)]
    pub taxdb_setting: Option<crate::algorithm::common::database::TaxdbSetting>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
    #[arg(long, default_value = "both")]
    pub strand: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // Read by the query options handler (`search_cli`), after the filtering handler.
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // Read by the database arguments handler when it reads the subjects (`run`).
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    /// Apply subject best hit filtering (disabled by default)
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:135
    #[arg(long = "subject_besthit", default_value_t = false)]
    pub subject_besthit: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    #[arg(long, value_parser = ncbi_integer)]
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
    // With `-culling_limit` NCBI refuses the command line, LOSAT's parser too (with its own
    // text, as for the other exclusions).
    /// Best Hit algorithm overhang value (recommended value: 0.1)
    #[arg(long, value_parser = ncbi_double, value_name = "float_value")]
    pub best_hit_overhang: Option<f64>,
    /// Best Hit algorithm score edge value (recommended value: 0.1)
    #[arg(long, value_parser = ncbi_double, value_name = "float_value")]
    pub best_hit_score_edge: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917,929-929
    /// Perform ungapped alignment only.
    #[arg(long)]
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:212-229,276-286
    // An omitted value keeps the default of the task (`coordination.rs` `task_defaults`: the
    // nucleotide handles of api/blast_nucl_options.cpp:163-193 and disc_nucl_options.cpp:66-84).
    // The ungapped X-drop of megablast is 0 (its handle sets the window only, and the options
    // are calloc'd, api/blast_options_local_priv.cpp:54): the raw X-drop is then the word cutoff
    // (blast_parameters.c:378-381). The values are bits; NCBI converts them to raw scores once,
    // with the smallest Lambda of the batch (`scoring.rs` `gap_x_dropoffs`, `run.rs`
    // `x_dropoff_init`).
    #[arg(long, value_parser = ncbi_double, help = "X-dropoff value (in bits) for ungapped extensions (default: 20 for blastn, dc-megablast and blastn-short; megablast uses the ungapped cutoff score)")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long, value_parser = ncbi_double, help = "X-dropoff value (in bits) for preliminary gapped extensions (default: 25 for megablast, 30 for the other tasks)")]
    pub xdrop_gap: Option<f64>,
    #[arg(long, value_parser = ncbi_double, help = "X-dropoff value (in bits) for final gapped alignment (default: 100)")]
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:663-665,681-685
    // The X-drops keep the task's defaults (megablast's greedy 25 bits). With the gap costs
    // 0/0 of megablast NCBI refuses the search (`scoring.rs` `check_scoring_options`).
    /// Use non-greedy dynamic programming extension.
    #[arg(long)]
    pub no_greedy: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:701-704,739-741
    // No constraint on the value: 0 and below keep the cutoff of the e-value
    // (`coordination.rs` `configure_task`).
    #[arg(long, value_parser = ncbi_integer, help = "Minimum raw gapped score to keep an alignment in the preliminary gapped and traceback stages")]
    pub min_raw_gapped_score: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:475-481,487-497
    // BLASTN has no matrix (no suggested window), so an omitted value keeps the default of the
    // task (`coordination.rs` `determine_window_size`).
    #[arg(long, value_parser = ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm (default: 40 for dc-megablast, 0 for the other tasks)")]
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:505-511,518-522
    // NCBI refuses a range without a window (`scoring.rs` `check_scoring_options`).
    #[arg(long, value_parser = ncbi_integer, help = "Number of off-diagonals to search for the 2nd hit, use 0 to turn off (default: 0)")]
    pub off_diagonal_range: Option<i32>,
    // NCBI blastn has no such option (the LOSAT-only option was removed, AGENTS.md rule 5).
    #[arg(skip = false)]
    pub verbose: bool,
    /// Scan stride for subject sequence scanning (NCBI BLAST optimization).
    /// Higher values skip more positions, reducing k-mer lookups but potentially missing some seeds.
    /// Default: 0 (auto-calculate based on word_size: 1 for word_size < 16, 4 for word_size >= 16).
    /// For megablast (word_size=28), scan_step=4 reduces lookups by ~4x with minimal sensitivity loss.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:397,1334:
    // Lookup construction owns this internal value; no CLI override.
    #[arg(skip = 0usize)]
    pub scan_step: usize,

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-782
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1098-1108
    // NCBI's argument is a string (c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660), parsed when the options
    // are set (`blastinput/app.rs` `parse_formatting_string` and `report_format`).
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,

    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    #[arg(skip)]
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The run's search strategy export (`-export_search_strategy`), made by `run`
    /// (`strategy::StrategyExport::for_run`).
    #[arg(skip)]
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,

    // The command line of the run as the SAM `@PG` line shows it: every argument followed by
    // one space, raw bytes (`GetCmdlineArgs`); `None` where there is no command line (the
    // library, ABI v1). Set by `run` from the process's arguments.
    //
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:877-885
    #[arg(skip)]
    pub command_line: Option<Vec<u8>>,
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
    // NCBI's application writes nothing to its standard error while it searches. LOSAT's
    // progress bar over the subjects of a parallel search (no NCBI counterpart) belongs to the
    // command line (`run`), which sets this; the library API and the web entries (`run_local`)
    // leave it off, so a search from memory writes nothing to the process's standard error.
    #[arg(skip)]
    pub progress_bar: bool,
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-264
    // `api::run_batches` hands over one of NCBI's query batches at a time, sized as the loop
    // sizes them (`QueryBatchSizes`); these records are that batch, whatever its size. No
    // argument of the command line sets this.
    #[arg(skip)]
    pub one_query_batch: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:394-399
    // The search of the repeats of `-filtering_db` (`repeat_filter`) is a `CLocalBlast` whose
    // alignments its caller reads: no formatter reads its results, so nothing that only the
    // formatter reads (`PRE_FETCH_SEQS_LIMIT`) stops it. Its caller takes the hits
    // (`blast_engine::repeat_search_hits`). No argument of the command line sets this.
    #[arg(skip)]
    pub hits_only: bool,
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input.hpp:180-182
    // The position in the whole input of the first of these records, as a reader's counter
    // starts there (`api::run_batches`): a query's `Query_N` in the messages. 0 on the command
    // line, whose records are the whole input.
    #[arg(skip)]
    pub first_query_ordinal: usize,
    /// LOSAT's refusal of a record of these records with the id of a query of an earlier call of
    /// the run (`query_batch::ScopeRefusal`), where the library searches each query batch as a call
    /// of its own (`api::run_batches`): the call's add of its batch gives it, so that the batch
    /// stops in the command line's form (line WQA D10 round 5, 2026-11-20). `None` on the command
    /// line. Not a public option.
    #[arg(skip)]
    pub(crate) scope_refusal: Option<crate::blastinput::query_batch::ScopeRefusal>,
}

impl BlastnArgs {
    /// The running application: the program of the command line, the library call or the Web run,
    /// which is not this program's when it runs another program's strategy
    /// (`ImportRun::app`). Its name prefixes the run's diagnostics and its formatting handler
    /// checks the formats.
    pub(crate) fn app(&self) -> crate::blastinput::app::Program {
        crate::report::archive::read::import::ImportRun::running(
            self.import.as_ref(),
            crate::blastinput::app::Program::Blastn,
        )
    }

    /// Whether the pairwise report (outfmt 0, `CBlastFormat`) cites the megablast paper: the
    /// formatter is made by the running application, and BLASTN's makes it with its `-task`
    /// megablast, TBLASTX's never does. So tblastx running a megablast strategy cites Gapped
    /// BLAST here (`import::apply`, probe `C.tblastx.from.blastn.db`).
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:244
    /// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:158-161
    pub(crate) fn pairwise_megablast_reference(&self) -> bool {
        self.task == "megablast" && self.app() == crate::blastinput::app::Program::Blastn
    }

    /// Whether the XML and JSON reports (outfmt 5, 15 and 16) cite the megablast paper: their
    /// formatters read the program of the options, which is the strategy's task, whichever
    /// application runs it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:653-668
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:637-653
    pub(crate) fn structured_megablast_reference(&self) -> bool {
        self.task == "megablast"
    }

    /// The `-subject` path (empty before NCBI's check that it is given).
    pub fn subject_path(&self) -> &std::path::Path {
        self.subject.as_deref().unwrap_or(std::path::Path::new(""))
    }

    /// The `-db` value as NCBI keeps it: without the white space at its ends; `None` for a
    /// search of `-subject`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    pub fn database_name(&self) -> Option<&str> {
        self.db.as_deref().map(crate::seqdb::db_argument_name)
    }

    /// The name of the MegaBLAST index that NCBI's handler asks for, and whether it is an
    /// old-style one: with `-use_index true`, task megablast and a database (with `-subject` or
    /// another task NCBI ignores `-use_index`; `blast_engine::mb_index`), the `-index_name` value
    /// (old style) or else the `-db` value as given (new style first). A run of an imported
    /// strategy asks for the import's index (`ImportRun::mb_index`: the strategy's, or the one
    /// that the command line names over it).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3381-3414
    pub(crate) fn megablast_index_name(&self) -> Option<(&str, bool)> {
        if let Some(run) = &self.import {
            return run
                .mb_index
                .as_ref()
                .map(|(name, old_style)| (name.as_str(), *old_style));
        }
        if self.use_index != Some(true) || self.task != "megablast" {
            return None;
        }
        match (&self.index_name, &self.db) {
            (Some(name), Some(_)) => Some((name, true)),
            (None, Some(db)) => Some((db, false)),
            (_, None) => None,
        }
    }

    /// The database length of the options: the `-dbsize` value, 0 without it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2570-2572
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:275
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

    /// The best hit filter of `-best_hit_overhang` and `-best_hit_score_edge` (the one not given
    /// takes its default), `None` when neither is given. The BLAST archive records the options
    /// as given.
    pub(crate) fn best_hit(&self) -> Option<crate::core::hspfilter_besthit::BestHitParams> {
        crate::core::hspfilter_besthit::BestHitParams::from_options(
            self.best_hit_overhang,
            self.best_hit_score_edge,
        )
    }

    /// `-searchsp`, 0 without it (`CBlastOptions::SetEffectiveSearchSpace`).
    pub fn options_searchsp(&self) -> i64 {
        self.searchsp.unwrap_or(0)
    }

    pub fn db_length(&self) -> i64 {
        self.dbsize.unwrap_or(0)
    }

    /// `-soft_masking` (`mask_at_hash`): `false` writes the query's masks into the searched copy
    /// of the query (`mask_query_letters`). The task's `ClearFilterOptions` does not change the
    /// value: blastn-short takes the default as well.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    pub fn soft_masking(&self) -> bool {
        self.soft_masking.unwrap_or(true)
    }

    /// Reads the `-dust` value (`parse_dust_filtering`) where NCBI's filtering handler reads
    /// it, with NCBI's error.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
    ///
    /// Without `-dust`, blastn-short does not filter: its task clears the filtering
    /// options (`coordination.rs` `task_dust_by_default`), and the argument has no default.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:344-349,410
    pub fn resolve_dust(&mut self) -> anyhow::Result<()> {
        if self.dust_filtering.is_none()
            && !crate::algorithm::blastn::coordination::task_dust_by_default(&self.task)
        {
            self.dust = DustSpec::No;
        }
        if let Some(value) = self.dust_filtering.take() {
            // See docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md exception 12.
            let footer = crate::wording::USER_MANUAL_FOOTER;
            self.dust =
                parse_dust_filtering(&value).map_err(|message| crate::cli::NativeError {
                    exit: 1,
                    message: format!("{QUERY_OPTIONS_ERROR} {message}\n{footer}\n").into(),
                })?;
        }
        Ok(())
    }

    /// The filtering handler's `-window_masker_taxid` and `-window_masker_db`, after `-dust`:
    /// NCBI's check of the counts file of `-window_masker_db` (`window_masker::open`), then
    /// that only one of them is given, and the masker made from the file (that of the taxid:
    /// `window_masker::open_taxid`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-458
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:430-466
    /// `-window_masker_taxid` counts also with the taxid 0, which uses no WindowMasker
    /// (`Blast_FindWindowMaskerLoc_Fwd` reads a taxid only when it is not 0); its file is
    /// looked for in the set-up of the query batches (`window_masker::open_taxid`).
    /// `-filtering_db` counts with any value, an empty one too (its database is opened by the
    /// repeat search of each query batch, `repeat_filter`).
    pub fn resolve_window_masker(&mut self) -> anyhow::Result<()> {
        if self.window_masker_db.is_none() && self.window_masker_taxid.is_none() {
            return Ok(());
        }
        let source = super::window_masker::source_of(self.window_masker_counts.as_ref());
        let database = match &self.window_masker_db {
            Some(name) => Some(super::window_masker::open(
                name,
                &*source,
                self.app().ncbi_name(),
            )?),
            None => None,
        };
        let filter_dbs = usize::from(self.filtering_db.is_some())
            + usize::from(self.window_masker_taxid.is_some())
            + usize::from(database.is_some());
        if filter_dbs > 1 {
            // See docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md exception 12.
            return Err(crate::cli::NativeError {
                exit: 1,
                message: format!(
                    "{QUERY_OPTIONS_ERROR} Please specify at most one of filtering_db, window_masker_taxid, or window_masker_db.\n{}\n",
                    crate::wording::USER_MANUAL_FOOTER
                ).into(),
            }
            .into());
        }
        let masker = match self.window_masker_taxid.filter(|&taxid| taxid != 0) {
            Some(taxid) => {
                // The command line and the library give the registry, LOSAT Web the registry of
                // its counts files; a caller that gives none is refused before
                // (`run_local_formats`).
                let registry = self
                    .window_masker_registry
                    .as_ref()
                    .ok_or_else(|| {
                        anyhow::anyhow!(
                            "-window_masker_taxid: the caller gave no directory of WindowMasker counts files"
                        )
                    })?
                    .as_ref()
                    .map_err(|reason| anyhow::anyhow!("{reason}"))?;
                Some(super::window_masker::open_taxid(
                    taxid,
                    registry,
                    &*source,
                    self.app().ncbi_name(),
                )?)
            }
            None => database,
        };
        self.window_masker = masker.map(std::sync::Arc::new);
        Ok(())
    }
}

/// A discontiguous megablast template type: `coding`, `optimal` or `coding_and_optimal`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:688-693,708-714
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:747-752
pub fn blastn_template_type(value: &str) -> Result<DiscWordType, String> {
    match value {
        "coding" => Ok(DiscWordType::Coding),
        "optimal" => Ok(DiscWordType::Optimal),
        "coding_and_optimal" => Ok(DiscWordType::TwoTemplates),
        _ => Err("expected one of coding, coding_and_optimal, optimal".into()),
    }
}

#[cfg(test)]
#[path = "args_inline_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "args_tests.rs"]
mod args_tests;
