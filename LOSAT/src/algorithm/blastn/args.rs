use super::disc_lookup::DiscWordType;
use crate::blastinput::value_parsers::*;
use crate::utils::os_bytes::os_str_bytes;
use crate::wording::QUERY_OPTIONS_ERROR;
use std::ffi::OsString;
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

    /// The value of `-strand` of arguments that have passed `parse` (NCBI's argument reader, or
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
#[derive(Debug, Clone)]
pub struct BlastnArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:47
    // `-` is standard input (`blastn/blast_engine/run.rs`).
    pub query: PathBuf,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastn/blast_engine/run.rs`).
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). NCBI's argument reader reports the
    // exclusion (`cli::reader`). The database is searched by `blastn/blast_engine/database.rs`.
    pub db: Option<OsString>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    /// `-db_soft_mask` and `-db_hard_mask` (`blastinput::db_mask`; series NP, M4a-2).
    pub db_mask: crate::blastinput::db_mask::DbMaskArgs,
    /// The database masking of the search as the formatter resolved it
    /// (`DbMaskOption::resolve`, in `run::search`); `None` before, and for a search without it.
    pub db_mask_resolved: Option<crate::blastinput::db_mask::ResolvedMask>,
    /// The MegaBLAST database index that the search uses (`-use_index true` that NCBI's loader
    /// kept, `blast_engine::mb_index`, set in `run::search`); `None` before, and for a search
    /// without one.
    pub(crate) megablast_index_loaded: Option<std::sync::Arc<crate::dbindex::search::LoadedIndex>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // No constraint and no dependency: read with `-db` and with `-subject`
    // (`BlastnArgs::db_length`).
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663
    pub show_gis: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2692-2693
    // NCBI's `CBlastFormat::m_IsHTML` (blast_format.cpp:90,224), for `-outfmt` 0 to 4 and 18.
    pub html: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2050-2052
    // `-parse_deflines` is a flag of the programs (blast_args.cpp:1961).
    pub parse_deflines: bool,
    /// The ids of the search's records under `-parse_deflines` (`run`).
    pub id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    pub searchsp: Option<i64>,
    // The help texts name the tasks (the Web form shows the short one); `describe` lists NCBI's
    // set (`ncbi_app::args::choices`, in the order of NCBI's BLASTN page).
    pub task: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:708-730
    // NCBI applies them to any task (`coordination.rs` `determine_template`).
    pub template_type: Option<DiscWordType>,
    pub template_length: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-170,288-300
    // An omitted value keeps the default of the task (`coordination.rs`).
    pub word_size: Option<usize>,
    pub num_threads: usize,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3168-3181
    // `-mt_mode` is read and has no effect: LOSAT's output is the same at every thread count
    // and for every value, and equals NCBI's `-num_threads 1` output (Owner 2026-10-13;
    // CLI DECISIONS 5). NCBI's `-mt_mode 1` (and `0` when it chooses to split by queries)
    // runs the queries in batches and writes the end of an outfmt 0 report or each XML
    // document per batch (`blastn_app.cpp:134-170`, `blastx_app.cpp` `x_RunMTBySplitQuery`);
    // that framing is not reproduced. The value must be 0 to 2 and `-num_threads` must be
    // given with it, as in NCBI (NCBI's argument reader checks the dependency, `cli::reader`).
    pub mt_mode: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3354-3363
    /// `-use_index` as NCBI's handler reads it (`blast_engine::mb_index`): ignored with
    /// `-subject` and for a task other than megablast; with `-db` and megablast, NCBI loads the
    /// database's index. When it loads none, or the word size is below the index's, NCBI warns
    /// and searches without it; otherwise it searches with the index. LOSAT does both; with
    /// `-db_hard_mask` its search reads the masked letters as matching nothing (approved exception
    /// 24 of PD-LOSAT-NCBI-DEFECTS, `dbindex::search`), and with an index of only some database
    /// volumes every query batch and query chunk searches the index anew (approved exception 30,
    /// `blast_engine::mb_index_seeds`).
    pub use_index: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3359-3362
    /// `-index_name`: the old-style index of `-use_index true` (`megablast_index_name`); without
    /// `-use_index true` NCBI ignores it.
    pub index_name: Option<OsString>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:142-146
    // An omitted value keeps the default of the task (`coordination.rs` `determine_evalue`).
    pub evalue: Option<f64>,
    /// Percent identity threshold for filtering HSPs (default: 0, no filtering)
    // Blast_HSPTest. The default 0.0 disables it, as NCBI's calloc'ed option does.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    pub percent_identity: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2924-2927,2960-2968
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:219,221
    // An omitted value keeps the default hit list size (500, `hitlist_size`) but the
    // pairwise report then shows 250 alignments, so the option has no default
    // (`Args::given_option`).
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    // The hit list size (NCBI BLAST hitlist_size) and the numbers of descriptions and
    // alignments of the pairwise report, which `process_options` sets from
    // `-max_target_seqs`, `-num_descriptions`, `-num_alignments` and the output format
    // (`blastinput::report_counts`). Not a public option.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
    pub report_counts: crate::blastinput::report_counts::ReportCounts,
    // The `-html` links of the run, which `process_options` sets
    // (`blastinput::app::html_links`). Not a public option.
    pub html_links: Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
    // Remove word seeds with high frequency in the searched database (off).
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:257 (limit_lookup)
    // NCBI blastn has no such option: it is magicblast's (c++/src/algo/blast/blastinput/blast_args.cpp:1492-1501), so
    // the LOSAT option was removed (AGENTS.md rule 5).
    pub limit_lookup: bool,
    // Maximum database word count for lookup filtering.
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:172-174
    pub max_db_word_count: u8,
    /// Maximum number of HSPs per subject (unlimited when omitted)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    pub max_hsps_per_subject: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201,313-315: no filter when omitted or 0.
    pub qcov_hsp_perc: Option<f64>,
    // Minimum diagonal separation between HSPs on the same subject (0 = auto, task-specific)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    // Default: 50 for blastn, 6 for megablast
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:239,259:
    // Task resolution owns this internal value; no CLI override.
    pub min_diag_separation: usize,
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    // The search strategy file (`report::archive::strategy`), created by `run` before the other
    // options are read.
    pub export_search_strategy: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:648-660,674-679
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:175-182,262-271
    // An omitted option keeps the default of the task (`coordination.rs`); a given one is
    // checked as NCBI checks it (`scoring.rs`). NCBI keeps the reward and the penalty in 16
    // bits; reward 0 with penalty 0 is NCBI's matrix-only scoring (`-task rmblastn`), a reward
    // of 0 or less with another penalty (no valid query) is rejected after NCBI's checks
    // (`scoring.rs`).
    pub reward: Option<i32>,
    pub penalty: Option<i32>,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:410-420
    // NCBI reads -dust in its filtering handler, after the output is opened
    // (`resolve_dust`); an omitted value keeps the default (20 64 1).
    pub dust_filtering: Option<OsString>,
    pub dust: DustSpec,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1939-1942
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:365-368
    // Read by `soft_masking`.
    pub soft_masking: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:350-352
    // Read by the filtering handler after `-dust` and before the WindowMasker options, which
    // count it (`resolve_window_masker`); the repeat search runs where the masks of each query
    // batch are made (`repeat_filter`).
    pub filtering_db: Option<OsString>,
    /// The repeat search of `-filtering_db` (`repeat_filter::RepeatFilter`): the database's name
    /// and the settings that look it up, given by the caller that has them (`resolve_repeat_filter`);
    /// `None` without the option.
    pub repeat_filter: Option<std::sync::Arc<super::repeat_filter::RepeatFilter>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:354-356
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2285
    // Read by the filtering handler after `-dust` (`resolve_window_masker`); 0 is no
    // WindowMasker, but it is one of the filtering databases given. The database arguments
    // exclude it with `-no_taxid_expansion` (`DbFilterArgs`, shared with programs that have no
    // `-window_masker_taxid`): an argument error, also with `-dryrun` (NCBI BLAST+ 2.17.0, line RG
    // item E1b probe `blastn.x.no_taxid_expansion.window_masker_taxid`: USAGE, rc 1).
    pub window_masker_taxid: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:358-363
    // Read by the filtering handler after `-dust` (`resolve_window_masker`).
    pub window_masker_db: Option<OsString>,
    /// The masker of `-window_masker_db` or of `-window_masker_taxid` (`resolve_window_masker`);
    /// `None` without them (or with the taxid 0).
    pub window_masker: Option<std::sync::Arc<super::window_masker::WindowMasker>>,
    /// Where `-window_masker_db` and the file of `-window_masker_taxid` are read from: the file
    /// system when `None`; LOSAT Web's adapter gives the counts files it holds under their names.
    pub window_masker_counts: Option<std::sync::Arc<dyn crate::utils::winmask::CountsSource>>,
    /// The registry of `-window_masker_taxid` as the caller read it (the directory of its files
    /// and what each load writes, `ncbi_environment::ApplicationSettings::window_masker_registry`),
    /// or LOSAT Web's registry of the counts files it holds (`window_masker_counts`: its directory,
    /// no posts); `None` when the caller gives neither, which `-window_masker_taxid` other than 0
    /// then refuses (`run_local_formats`).
    pub window_masker_registry:
        Option<Result<crate::blastinput::ncbi_environment::WindowMaskerRegistry, String>>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath` of `taxdb.bti`), asked
    /// before the search where a tabular field list looks every taxid up (the plural names of the
    /// taxonomy database, the organism report; `algorithm::blastn::tabular`), by the command line
    /// and the library (`common::database::subject_taxonomy`); `None` elsewhere (LOSAT Web has no
    /// files: not found). Not a public option.
    pub tabular_taxdb: Option<bool>,
    /// The taxonomy database of a `-subject` search's organism report (outfmt 18), made with
    /// `tabular_taxdb` (`common::database::SubjectTaxdb`); `None` elsewhere. Not a public option.
    pub subject_taxdb: Option<crate::algorithm::common::database::SubjectTaxdb>,
    /// The `[BLAST] BLASTDB` setting through which a `-subject` search's tabular reports look for
    /// the taxonomy database (`common::database::subject_taxonomy`), which a search of records of
    /// the library (`api::run`) reads from the call's application settings when a report looks
    /// every taxid up; `None` elsewhere (the command line passes its settings, LOSAT Web has no
    /// files). Not a public option.
    pub taxdb_setting: Option<crate::algorithm::common::database::TaxdbSetting>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
    pub strand: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // Read by the query options handler (`search_cli`), after the filtering handler.
    pub query_loc: Option<OsString>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // Read by the database arguments handler when it reads the subjects (`run`).
    pub subject_loc: Option<OsString>,
    /// Apply subject best hit filtering (disabled by default)
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:135
    pub subject_besthit: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
    // With `-culling_limit` NCBI refuses the command line, LOSAT's parser too (with its own
    // text, as for the other exclusions).
    /// Best Hit algorithm overhang value (recommended value: 0.1)
    pub best_hit_overhang: Option<f64>,
    /// Best Hit algorithm score edge value (recommended value: 0.1)
    pub best_hit_score_edge: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917,929-929
    /// Perform ungapped alignment only.
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:212-229,276-286
    // An omitted value keeps the default of the task (`coordination.rs` `task_defaults`: the
    // nucleotide handles of api/blast_nucl_options.cpp:163-193 and disc_nucl_options.cpp:66-84).
    // The ungapped X-drop of megablast is 0 (its handle sets the window only, and the options
    // are calloc'd, api/blast_options_local_priv.cpp:54): the raw X-drop is then the word cutoff
    // (blast_parameters.c:378-381). The values are bits; NCBI converts them to raw scores once,
    // with the smallest Lambda of the batch (`scoring.rs` `gap_x_dropoffs`, `run.rs`
    // `x_dropoff_init`).
    pub xdrop_ungap: Option<f64>,
    pub xdrop_gap: Option<f64>,
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:663-665,681-685
    // The X-drops keep the task's defaults (megablast's greedy 25 bits). With the gap costs
    // 0/0 of megablast NCBI refuses the search (`scoring.rs` `check_scoring_options`).
    /// Use non-greedy dynamic programming extension.
    pub no_greedy: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:701-704,739-741
    // No constraint on the value: 0 and below keep the cutoff of the e-value
    // (`coordination.rs` `configure_task`).
    pub min_raw_gapped_score: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:475-481,487-497
    // BLASTN has no matrix (no suggested window), so an omitted value keeps the default of the
    // task (`coordination.rs` `determine_window_size`).
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:505-511,518-522
    // NCBI refuses a range without a window (`scoring.rs` `check_scoring_options`).
    pub off_diagonal_range: Option<i32>,
    // NCBI blastn has no such option (the LOSAT-only option was removed, AGENTS.md rule 5).
    pub verbose: bool,
    /// Scan stride for subject sequence scanning (NCBI BLAST optimization).
    /// Higher values skip more positions, reducing k-mer lookups but potentially missing some seeds.
    /// Default: 0 (auto-calculate based on word_size: 1 for word_size < 16, 4 for word_size >= 16).
    /// For megablast (word_size=28), scan_step=4 reduces lookups by ~4x with minimal sensitivity loss.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:397,1334:
    // Lookup construction owns this internal value; no CLI override.
    pub scan_step: usize,

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:772-782
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1098-1108
    // NCBI's argument is a string (c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660), parsed when the options
    // are set (`blastinput/app.rs` `parse_formatting_string` and `report_format`).
    pub outfmt: OsString,

    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The run's search strategy export (`-export_search_strategy`), made by `run` and by LOSAT
    /// Web's `web_adapter::give_strategy` (`strategy::StrategyExport::for_args`).
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,

    // The command line of the run as the SAM `@PG` line shows it: every argument followed by
    // one space, raw bytes (`GetCmdlineArgs`); `None` where there is no command line (the
    // library, ABI v1). Set by `run` from the process's arguments.
    //
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:877-885
    pub command_line: Option<Vec<u8>>,
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
    // NCBI's application writes nothing to its standard error while it searches. LOSAT's
    // progress bar over the subjects of a parallel search (no NCBI counterpart) belongs to the
    // command line (`run`), which sets this; the library API and the web entries (`run_local`)
    // leave it off, so a search from memory writes nothing to the process's standard error.
    pub progress_bar: bool,
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-264
    // `api::run_batches` hands over one of NCBI's query batches at a time, sized as the loop
    // sizes them (`QueryBatchSizes`); these records are that batch, whatever its size. No
    // argument of the command line sets this.
    pub one_query_batch: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:394-399
    // The search of the repeats of `-filtering_db` (`repeat_filter`) is a `CLocalBlast` whose
    // alignments its caller reads: no formatter reads its results, so nothing that only the
    // formatter reads (`PRE_FETCH_SEQS_LIMIT`) stops it. Its caller takes the hits
    // (`blast_engine::repeat_search_hits`). No argument of the command line sets this.
    pub hits_only: bool,
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input.hpp:180-182
    // The position in the whole input of the first of these records, as a reader's counter
    // starts there (`api::run_batches`): a query's `Query_N` in the messages. 0 on the command
    // line, whose records are the whole input.
    pub first_query_ordinal: usize,
    /// LOSAT's refusal of a record of these records with the id of a query of an earlier call of
    /// the run (`query_batch::ScopeRefusal`), where the library searches each query batch as a call
    /// of its own (`api::run_batches`): the call's add of its batch gives it, so that the batch
    /// stops in the command line's form (line WQA D10 round 5, 2026-11-20). `None` on the command
    /// line. Not a public option.
    pub(crate) scope_refusal: Option<crate::blastinput::query_batch::ScopeRefusal>,
}

impl BlastnArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, blastn's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(BlastnArgs {
            query: args.option_or_default("query")?,
            subject: args.given_option("subject")?,
            db: args.given_option("db")?,
            db_filter: crate::blastinput::db_filter::DbFilterArgs::from_args(args)?,
            db_mask: crate::blastinput::db_mask::DbMaskArgs::from_args(args)?,
            db_mask_resolved: Default::default(),
            megablast_index_loaded: Default::default(),
            dbsize: args.given_option("dbsize")?,
            show_gis: args.flag("show_gis"),
            html: args.flag("html"),
            parse_deflines: args.flag("parse_deflines"),
            id_scope: Default::default(),
            searchsp: args.given_option("searchsp")?,
            task: args.option_or_default("task")?,
            template_type: args.given_option("template_type")?,
            template_length: args.given_option("template_length")?,
            word_size: args.given_option("word_size")?,
            num_threads: args.option_or_default("num_threads")?,
            mt_mode: args.given_option("mt_mode")?,
            use_index: args.given_option("use_index")?,
            index_name: args.given_option("index_name")?,
            evalue: args.given_option("evalue")?,
            percent_identity: args.given_option("perc_identity")?,
            max_target_seqs: args.given_option("max_target_seqs")?,
            num_descriptions: args.given_option("num_descriptions")?,
            num_alignments: args.given_option("num_alignments")?,
            formatting: crate::blastinput::report_counts::FormattingArgs::from_args(args)?,
            report_counts: Default::default(),
            html_links: Default::default(),
            limit_lookup: false,
            max_db_word_count: 30u8,
            max_hsps_per_subject: args.given_option("max_hsps")?,
            qcov_hsp_perc: args.given_option("qcov_hsp_perc")?,
            min_diag_separation: 0usize,
            out: args.given_option("out")?,
            export_search_strategy: args.given_option("export_search_strategy")?,
            reward: args.given_option("reward")?,
            penalty: args.given_option("penalty")?,
            gap_open: args.given_option("gapopen")?,
            gap_extend: args.given_option("gapextend")?,
            dust_filtering: args.given_option("dust")?,
            dust: DustSpec::Yes,
            lcase_masking: args.flag("lcase_masking"),
            soft_masking: args.given_option("soft_masking")?,
            filtering_db: args.given_option("filtering_db")?,
            repeat_filter: Default::default(),
            window_masker_taxid: args.given_option("window_masker_taxid")?,
            window_masker_db: args.given_option("window_masker_db")?,
            window_masker: Default::default(),
            window_masker_counts: Default::default(),
            window_masker_registry: Default::default(),
            tabular_taxdb: Default::default(),
            subject_taxdb: Default::default(),
            taxdb_setting: Default::default(),
            strand: args.option_or_default("strand")?,
            query_loc: args.given_option("query_loc")?,
            subject_loc: args.given_option("subject_loc")?,
            subject_besthit: args.flag("subject_besthit"),
            culling_limit: args.given_option("culling_limit")?,
            best_hit_overhang: args.given_option("best_hit_overhang")?,
            best_hit_score_edge: args.given_option("best_hit_score_edge")?,
            ungapped: args.flag("ungapped"),
            xdrop_ungap: args.given_option("xdrop_ungap")?,
            xdrop_gap: args.given_option("xdrop_gap")?,
            xdrop_gap_final: args.given_option("xdrop_gap_final")?,
            no_greedy: args.flag("no_greedy"),
            min_raw_gapped_score: args.given_option("min_raw_gapped_score")?,
            window_size: args.given_option("window_size")?,
            off_diagonal_range: args.given_option("off_diagonal_range")?,
            verbose: false,
            scan_step: 0usize,
            outfmt: args.option_or_default("outfmt")?,
            archive: Default::default(),
            strategy: Default::default(),
            command_line: Default::default(),
            progress_bar: Default::default(),
            one_query_batch: Default::default(),
            hits_only: Default::default(),
            first_query_ordinal: Default::default(),
            scope_refusal: Default::default(),
        })
    }
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
    pub fn database_name(&self) -> Option<Vec<u8>> {
        self.db
            .as_deref()
            .map(|db| crate::seqdb::db_argument_name_bytes(&os_str_bytes(db)).to_vec())
    }

    /// The name of the MegaBLAST index that NCBI's handler asks for, and whether it is an
    /// old-style one: with `-use_index true`, task megablast and a database (with `-subject` or
    /// another task NCBI ignores `-use_index`; `blast_engine::mb_index`), the `-index_name` value
    /// (old style) or else the `-db` value as given (new style first).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3381-3414
    pub(crate) fn megablast_index_name(&self) -> Option<(Vec<u8>, bool)> {
        if self.use_index != Some(true) || self.task != "megablast" {
            return None;
        }
        match (&self.index_name, &self.db) {
            (Some(name), Some(_)) => Some((os_str_bytes(name).into_owned(), true)),
            (None, Some(db)) => Some((os_str_bytes(db).into_owned(), false)),
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
            self.dust = parse_dust_filtering(os_str_bytes(&value)).map_err(|message| {
                crate::cli::NativeError {
                    exit: 1,
                    message: format!("{QUERY_OPTIONS_ERROR} {message}\n{footer}\n").into(),
                }
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
            Some(name) => Some(super::window_masker::open(&os_str_bytes(name), &*source)?),
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
                Some(super::window_masker::open_taxid(taxid, registry, &*source)?)
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

/// `-template_type`'s templates (`blastn_template_type` of the value, which NCBI's reader checked
/// against its set).
impl crate::ncbi_app::args::OptionValue for DiscWordType {
    const KINDS: &'static [crate::ncbi_app::args::Kind] = &[crate::ncbi_app::args::Kind::Text];
    fn of(value: &std::ffi::OsStr, _: Option<crate::ncbi_app::args::Number>) -> Option<Self> {
        blastn_template_type(value.to_str()?).ok()
    }
}

#[cfg(test)]
#[path = "args_inline_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "args_tests.rs"]
mod args_tests;
