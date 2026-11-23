//! BLASTX options and native serial local-subject search/report dispatch.
use crate::blastinput::value_parsers::{
    best_hit_value, blastx_query_coverage, blastx_real, blastx_threshold, genetic_code,
    ncbi_boolean, ncbi_input_path, ncbi_int8, ncbi_integer, ncbi_output_path,
    ncbi_positive_integer, nonnegative_ncbi_integer, protein_word_size, report_count, SegSpec,
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

#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:49-62
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:72-78
pub struct BlastxArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // The default `-` is standard input (`input_files::open_input`).
    #[arg(long, value_name = "PATH", default_value = "-", value_parser = ncbi_input_path())]
    pub query: PathBuf,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2366
    // Optional for the parser: without `-subject` and `-db` the run fails with NCBI's error
    // after the parsing (`blastinput/app.rs` `missing_subject_error`).
    #[arg(long, value_name = "PATH", value_parser = ncbi_input_path())]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // Read by the database arguments handler when it reads the subjects (`native::run`).
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        conflicts_with_all = crate::blastinput::db_filter::DB_FILTER_OPTIONS,
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
    // NCBI excludes `-subject` and `-subject_loc` with `-db` (c++/src/algo/blast/blastinput/blast_args.cpp:2364-2381). The
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
    /// `-db_soft_mask` and `-db_hard_mask` (`blastinput::db_mask`; series NP, M4b).
    #[command(flatten)]
    pub db_mask: crate::blastinput::db_mask::DbMaskArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // No constraint and no dependency: read with `-db` and with `-subject`.
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS", help = "Effective length of the database")]
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663
    #[arg(
        long = "show_gis",
        default_value_t = false,
        help = "Show NCBI GIs in deflines?"
    )]
    pub show_gis: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2050-2052
    // `-parse_deflines` is a flag of the programs (blast_args.cpp:1961).
    #[arg(
        long = "parse_deflines",
        default_value_t = false,
        help = "Should the query and subject defline(s) be parsed?"
    )]
    pub parse_deflines: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    #[arg(long, value_parser = crate::blastinput::value_parsers::nonnegative_ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3441-3443
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    // The search strategy file (`report::archive::strategy`), created by `run` before the other
    // options are read.
    #[arg(long, value_name = "filename", value_parser = ncbi_output_path(), help = "File name to record the search strategy used")]
    pub export_search_strategy: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:99-100
    #[arg(long, default_value = "blastx", value_parser = BLASTX_TASK)]
    pub task: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:146-146
    #[arg(long, default_value_t = 10.0, value_parser = blastx_real)]
    pub evalue: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-167
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3; 5 for blastx-fast]")]
    pub word_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:203-205
    #[arg(long, value_parser = ncbi_positive_integer)]
    pub max_hsps: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201
    #[arg(long, value_parser = blastx_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:241-243
    #[arg(long, value_parser = ncbi_boolean, action = clap::ArgAction::Set, help = "Use sum statistics [default: true]")]
    pub sum_stats: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:175-177
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11]")]
    pub gap_open: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:180-182
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1]")]
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:937-942
    #[arg(long, default_value_t = 0, value_parser = nonnegative_ncbi_integer)]
    pub max_intron_length: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-338
    // Read when NCBI's filtering handler reads it (`seg_spec`).
    #[arg(long, default_value = "12 2.2 2.5")]
    pub seg: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    #[arg(long, default_value = "false", value_parser = ncbi_boolean, action = clap::ArgAction::Set)]
    pub soft_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:531-534
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:578-583
    #[arg(
        long,
        value_parser = blastx_threshold,
        help = "Minimum word score such that the word is added to the BLAST lookup table [default: 12; 19.3, 21 or 20.25 for word sizes 5, 6 and 7; 20 for blastx-fast]"
    )]
    pub threshold: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:475-478
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40]")]
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229, 276-286
    // (cited at `blastp::args::BlastpArgs::xdrop_ungap`): `-xdrop_ungap`, `-xdrop_gap` and
    // `-xdrop_gap_final` in bits, given to the options when given. `None`: BLASTX's defaults
    // (`ResolvedOptions`).
    #[arg(long = "xdrop_ungap", value_name = "float_value", value_parser = blastx_real, help = "X-dropoff value (in bits) for ungapped extensions [default: 7]")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long = "xdrop_gap", value_name = "float_value", value_parser = blastx_real, help = "X-dropoff value (in bits) for preliminary gapped extensions [default: 15]")]
    pub xdrop_gap: Option<f64>,
    #[arg(long = "xdrop_gap_final", value_name = "float_value", value_parser = blastx_real, help = "X-dropoff value (in bits) for final gapped alignment [default: 25]")]
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3299
    #[arg(
        long,
        value_parser = nonnegative_ncbi_integer,
        help = "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit [default: 0, no culling]"
    )]
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
    // With `-culling_limit` NCBI refuses the command line, LOSAT's parser too (with its own
    // text, as for the other exclusions; BLASTP's `best_hit_overhang`).
    /// Best Hit algorithm overhang value (recommended value: 0.1)
    #[arg(long, value_parser = best_hit_value, value_name = "float_value", conflicts_with = "culling_limit")]
    pub best_hit_overhang: Option<f64>,
    /// Best Hit algorithm score edge value (recommended value: 0.1)
    #[arg(long, value_parser = best_hit_value, value_name = "float_value", conflicts_with = "culling_limit")]
    pub best_hit_score_edge: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329-3329
    #[arg(long)]
    pub subject_besthit: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917
    #[arg(long)]
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1941-1942
    #[arg(long)]
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    // NCBI's flag, off by default: the Smith-Waterman redo of each match (`kappa::sw_redo_context`,
    // item X4b-2), with every `-comp_based_stats`. RPSTBLASTN's own (`rpstblastn::args`) never
    // sets it here.
    #[arg(long)]
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
    #[arg(long, default_value = "both", value_parser = ["both", "plus", "minus"])]
    pub strand: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // Read by the query options handler (`native::run`).
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1021-1024
    #[arg(long, default_value_t = 1, value_parser = genetic_code)]
    pub query_gencode: u8,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
    // NCBI reads the value when it sets the options (`native::run`; `blastinput/app.rs`
    // `parse_formatting_string`).
    #[arg(long, default_value = "0", help = crate::cli::outfmt_help("0 to 16, 18 and 20", "6, 7, 10 and 20", ""))]
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2726-2730
    #[arg(
        long,
        value_parser = ncbi_positive_integer,
        help = "Maximum number of aligned sequences to keep [default: 500]",
        conflicts_with_all = ["num_descriptions", "num_alignments"]
    )]
    pub max_target_seqs: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    #[command(flatten)]
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3158-3161
    #[arg(long, default_value_t = 1, value_parser = ncbi_positive_integer, help = "Number of threads (CPUs); serial WASI accepts only 1")]
    pub num_threads: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3168-3181
    // `-mt_mode` is read and has no effect: LOSAT's output is the same at every thread count
    // and for every value, and equals NCBI's `-num_threads 1` output (Owner 2026-10-13;
    // CLI DECISIONS 5). NCBI's `-mt_mode 1` (and `0` when it chooses to split by queries)
    // runs the queries in batches and writes the end of an outfmt 0 report or each XML
    // document per batch (`blastn_app.cpp:134-170`, `blastx_app.cpp` `x_RunMTBySplitQuery`);
    // that framing is not reproduced. The value must be 0 to 2 and `-num_threads` must be
    // given with it, as in NCBI (clap's `requires` ignores a default value).
    #[arg(long, value_name = "int_value", value_parser = crate::blastinput::value_parsers::ncbi_mt_mode, requires = "num_threads", help = "Multi-thread mode: 0 auto, 1 split by queries, 2 split by database (accepted; the output does not depend on it)")]
    pub mt_mode: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`: the first
    // character chooses the mode, and the 'u' of unified P is BLASTP's only).
    #[arg(
        long,
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
    )]
    pub comp_based_stats: String,
    /// The RPS search of an RPSTBLASTN run (`algorithm::rpstblastn`), which takes the place of
    /// BLASTX's search of the query batches; `None` in every BLASTX run.
    #[arg(skip)]
    pub(crate) rpstblastn: Option<std::sync::Arc<dyn RpsTblastnSearch>>,
    /// The recorded options of RPSTBLASTN's search strategy requests (`archive::option_inputs`),
    /// made by its command line from its own arguments as given, with `-export_search_strategy`;
    /// `None` in every BLASTX run.
    #[arg(skip)]
    pub(crate) request_inputs: Option<Box<crate::report::archive::options::OptionInputs>>,
    /// The `[BLAST] BLASTDB` setting through which a `-subject` search's tabular reports look for
    /// the taxonomy database (`common::database::subject_taxonomy`), which a search of records of
    /// the library (`api::run`) reads from the call's application settings when a report looks
    /// every taxid up; `None` elsewhere (the command line passes its settings, LOSAT Web has no
    /// files). Not a public option.
    #[arg(skip)]
    pub taxdb_setting: Option<crate::algorithm::common::database::TaxdbSetting>,
    /// The settings by which a search of the library looks up the file of an unknown matrix
    /// (`ResolvedOptions::find_matrix_file`), from its application start: set by `api::run` for an
    /// ungapped search whose `-matrix` NCBI's tables do not have (`unknown_matrix_name`) and by
    /// `run_local_database`; `None` for the command line, which looks it up with its own
    /// settings (`native::run`), and LOSAT Web, which has no file system. Not a public option.
    #[arg(skip)]
    pub matrix_search: Option<crate::blastinput::ncbi_environment::MatrixSearch>,
    /// The position in the whole input of the first of these records, where the library searches
    /// each query batch as a call of its own (`api::run_batches`); 0 on the command line, whose
    /// records are the whole input. A later call's prefetch has no record of the loads of the
    /// calls before (`local::DatabaseReports::unrecorded_calls`). Not a public option.
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

/// The search of RPSTBLASTN's query batches (`rps_search::translated`, series PR item PS12a,
/// DECISIONS 38 D1): the user's nucleotide queries, in BLASTX's query batches with BLASTX's frames,
/// masks and statistics, against a profile database, whose results BLASTX's reports write. A run
/// that has one (`ResolvedOptions::rpstblastn`) is RPSTBLASTN's: the command line's loop of the
/// query batches (`native::run_query_batches`) is this search's, once per run; nothing of BLASTX's
/// own search runs. It is a trait of BLASTX's so that BLASTX does not depend on the RPS engine.
pub(crate) trait RpsTblastnSearch: std::fmt::Debug + Send + Sync {
    /// The run's search over `database`, once per run before the first byte of the reports
    /// (`native::run_inner`): LOSAT's refusals of the database's files and of the options come
    /// there, NCBI's own errors of the files where NCBI's search of a batch raises them.
    /// `requests` is whether the run writes an archive or a search strategy, whose requests hold
    /// what the search's set-up writes into the options handle (`handle_sets`).
    fn for_database(
        &self,
        options: &ResolvedOptions,
        database: Option<&super::database::BlastxDatabase>,
        requests: bool,
    ) -> Result<std::sync::Arc<dyn RpsTblastnSearch>>;
    /// `native::run_query_batches` of RPSTBLASTN: the query batches of `read`, searched against
    /// the profile volumes of `database`, each handed to `take` (the search strategy's request,
    /// before its search, and NCBI's posts of the bytes replaced in it to the diagnostics, after
    /// the batch's reader messages) and to `render` in the reader's order. Returns the number of
    /// queries reported.
    #[allow(clippy::too_many_arguments)]
    fn run_query_batches(
        &self,
        options: &ResolvedOptions,
        database: Option<&super::database::BlastxDatabase>,
        read: &mut super::native::QueryBatchReader<'_>,
        diagnostics: &mut dyn std::io::Write,
        search_error: &dyn Fn(anyhow::Error) -> anyhow::Error,
        render: &mut super::native::BatchRender<'_>,
        take: &mut dyn FnMut(
            &mut dyn std::io::Write,
            &[super::input::FastaRecord],
            usize,
        ) -> Result<()>,
    ) -> Result<usize>;
    /// What each query batch's set-up of the search writes into NCBI's options handle beyond the
    /// batch's search space (`rps_search::volumes::search_sets`), which the requests of the archive
    /// and of the search strategy hold (the archive's of a batch after its set-up, the strategy's
    /// after the first batch's): the search that `for_database` returned for a run that writes
    /// one knows it.
    fn handle_sets(&self) -> Vec<crate::report::archive::options::Param>;
    /// Whether the batches' results are merged over a database of several volumes, as NCBI's
    /// `CLocalRPSBlast` merges them: the merge keeps only NCBI's messages of error severity (no
    /// query warning, `report::query_warning`) and gives a batch that was not searched an empty
    /// alignment set (outfmt 7's `# 0 hits found`), while its statistics stay those of a batch
    /// that was not searched.
    fn merges_volumes(&self) -> bool;
}

/// The program of a run of BLASTX's driver as NCBI's and LOSAT's texts name it: BLASTX, or
/// RPSTBLASTN for an RPSTBLASTN run (`BlastxArgs::rpstblastn`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ProgramName {
    /// `BLASTX`: LOSAT's refusals ("not supported by LOSAT's BLASTX") and the reports' program
    /// lines.
    pub(crate) upper: &'static str,
    /// `blastx`: NCBI's messages (`Warning: [blastx] ...`) and its settings.
    pub(crate) lower: &'static str,
}

impl ProgramName {
    pub(crate) const BLASTX: Self = Self {
        upper: "BLASTX",
        lower: "blastx",
    };
    pub(crate) const RPSTBLASTN: Self = Self {
        upper: "RPSTBLASTN",
        lower: "rpstblastn",
    };

    fn of(rpstblastn: &Option<std::sync::Arc<dyn RpsTblastnSearch>>) -> Self {
        match rpstblastn {
            Some(_) => Self::RPSTBLASTN,
            None => Self::BLASTX,
        }
    }
}

impl ResolvedOptions {
    /// The program of the run's texts (`ProgramName`).
    pub(crate) fn program(&self) -> ProgramName {
        ProgramName::of(&self.rpstblastn)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:78-81 (cited at
    // `report::tabular::TabularSpec::std`)
    /// The report the tabular writer writes (`tabular::write_blastx_tabular`): `tabular`, or for
    /// the default outfmt 6 and 7 the `std` fields and a tab of the same formatter
    /// (`TabularSpec::std`); `None` for a format that is not tabular. `tabular` itself stays
    /// `None` for the default report, as the shared rules of the default fields read it
    /// (`report::tabular::check_failed_lookups`, `TabularIds::of`).
    pub(crate) fn tabular_report(&self) -> Option<&'static crate::report::tabular::TabularSpec> {
        match self.tabular {
            Some(spec) => Some(spec),
            None if matches!(self.outfmt, 6 | 7) => {
                Some(crate::report::tabular::TabularSpec::std(self.outfmt == 7))
            }
            None => None,
        }
    }

    /// Whether the run's results are merged over several volumes (`RpsTblastnSearch::merges_volumes`).
    pub(crate) fn merged_volumes(&self) -> bool {
        self.rpstblastn
            .as_ref()
            .is_some_and(|search| search.merges_volumes())
    }

    /// The number of a query's lists that the reports NCBI's formatter cuts write: the merge of
    /// several volumes keeps each volume's lists (`s_MergeAlignSet`), which the formatter cuts to
    /// the hit list size (`BatchResults::shown_lists`): `PruneSeqalign` before the tabular
    /// reports (outfmt 7 counts their hits before it, as rpsblast's `psiblast_hits_found`) and
    /// the Seq-annots (8, 9, 12), `TrimSeqAlign` before BLAST XML, XML2 and JSON (5, 13 to 16);
    /// the archive (11) and the organism report (18) take every list. All of them otherwise.
    pub(crate) fn shown_lists(&self) -> usize {
        if self.merged_volumes() {
            usize::try_from(self.hitlist_size).unwrap_or(0)
        } else {
            usize::MAX
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:184 (c++/src/app/blast/blastx_app.cpp:247-249
    // hands the sorts over too, cited at `formatting`)
    /// The sorts of the pairwise report (`-sorthits`, `-sorthsps`): BLASTX's; none for an
    /// RPSTBLASTN run, whose application hands its formatter the line length only, as rpsblast's
    /// and psiblast's (`blastp::args`, `hit_sort::ReportSort`): its sorts are read with NCBI's
    /// warnings (`native::warnings_of`) and sort nothing.
    pub(crate) fn report_sort(&self) -> crate::report::hit_sort::ReportSort {
        match self.rpstblastn {
            Some(_) => Default::default(),
            None => self.formatting.report_sort(),
        }
    }

    /// The view of the pairwise report (outfmt 0 to 4, `PairwiseView::of`), `None` for the other
    /// formats: formats 1 to 4 are outfmt 0's report with the query-anchored section in place of
    /// its alignments.
    pub(crate) fn pairwise_view(&self) -> Option<crate::report::pairwise::PairwiseView> {
        crate::blastinput::app::ReportFormat::from_number(i32::from(self.outfmt))
            .and_then(crate::report::pairwise::PairwiseView::of)
    }

    /// How the pairwise report reads each query's subjects, as its writer reads them
    /// (`PairwiseView::reading` of `-sorthits`, `-sorthsps`, `-num_descriptions` and
    /// `-num_alignments`), `None` for the other formats.
    pub(crate) fn pairwise_reading(&self) -> Option<crate::report::pairwise::PairwiseReading> {
        self.pairwise_view().map(|view| {
            view.reading(
                self.report_sort(),
                crate::report::db_subject::Shown {
                    descriptions: usize::try_from(self.num_descriptions).unwrap_or(0),
                    alignments: usize::try_from(self.num_alignments).unwrap_or(0),
                },
            )
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:838-871
/// The matrix of a `-matrix` name, as NCBI's matrix loader finds it (`None`: a name it does not
/// know).
///
/// NCBI reference (598d8ae6): c++/src/util/tables/raw_scoremat.c:131-166
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
/// The fields that ABI v1 accepts (`web::v1_outfmt`, frozen, plan TD-1): the custom fields
/// of BLASTX's report before it wrote them through the shared tabular writer.
#[cfg(any(target_arch = "wasm32", test))]
pub const FIELDS: &[&str] = &[
    "qseqid", "qacc", "qaccver", "qlen", "sseqid", "sacc", "saccver", "slen", "qstart", "qend",
    "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident",
    "mismatch", "positive", "gapopen", "gaps", "ppos", "qframe", "sframe", "frames", "btop",
    "stitle",
];

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blastx_options.cpp:57-70
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:131-146
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
    /// (`ncbi_environment::ApplicationSettings::matrix_file`), found by the command line and the
    /// library API after the options are resolved (`find_matrix_file`); `None` in LOSAT Web. The
    /// first query batch is refused then
    /// (`query_setup::matrix_error_message`).
    pub matrix_file: Option<std::path::PathBuf>,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gapped: bool,
    pub composition: u8,
    /// `-use_sw_tback` (NCBI's `SetSmithWatermanMode`, blast_args.cpp:885): the redo of each
    /// match is NCBI's Smith-Waterman redo, with or without composition-based statistics.
    pub use_sw_tback: bool,
    pub seg: SegSpec,
    pub soft_masking: bool,
    pub strand: String,
    pub query_gencode: u8,
    pub window_size: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:123-145
    pub x_dropoff: f64,
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub gap_trigger: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
    pub db_length: i64,
    pub effective_search_space: i64,
    /// Whether `-searchsp` was given, at any value, or the environment variable `OLD_FSC` is set
    /// (`app::old_fsc_environment`): NCBI then has no Gumbel block (`has_gumbel_block`).
    pub old_fsc: bool,

    pub evalue: f64,
    pub sum_stats: bool,
    pub max_intron_length: i32,
    pub max_hsps: i32,
    /// The -qcov_hsp_perc value, or 0 (no filter).
    pub query_cov_hsp_perc: f64,
    pub culling_limit: i32,
    /// The best hit filter of `-best_hit_overhang` and `-best_hit_score_edge`, `None` without
    /// them (`BlastxArgs::best_hit`).
    pub(crate) best_hit: Option<crate::core::hspfilter_besthit::BestHitParams>,
    pub subject_besthit: bool,
    pub hitlist_size: i32,
    pub num_descriptions: i32,
    pub num_alignments: i32,
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:247-249
    /// `-line_length`, `-sorthits` and `-sorthsps` of the pairwise report.
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    pub outfmt: u8,
    /// `-show_gis`: the labels of the query-anchored sections (outfmt 1 to 4) show a sequence's GI.
    pub show_gis: bool,
    /// The tabular report with a field list or a delimiter (outfmt 6 and 7), or outfmt 10 or
    /// 20 (`TabularSpec::of_app`); `None` for the default outfmt 6 and 7 (the `std` fields and a
    /// tab, `TabularSpec::std`, which the reports write with the same shared writer, `tabular.rs`)
    /// and the other formats.
    pub tabular: Option<&'static crate::report::tabular::TabularSpec>,
    /// Whether `taxdb.bti` is on NCBI's search path, for a `tabular` field list that shows a
    /// plural name of the taxonomy database (its formatters warn without it) or the organism
    /// report; set where the command line asks (`native::run`; `-db`: the database's search
    /// path, `-subject`: NCBI's, `common::database::subject_taxonomy`), by the library's searches
    /// too. False otherwise, and in LOSAT Web, which has no files.
    pub taxdb_found: bool,
    /// The taxonomy database of a `-subject` search's organism report (outfmt 18), set with
    /// `taxdb_found` (`common::database::SubjectTaxdb`); `None` otherwise.
    pub subject_taxdb: Option<crate::algorithm::common::database::SubjectTaxdb>,
    /// A list kept no sequence of the database (`SeqDb::filtered_empty`): the queries of the
    /// searched batches have `EMPTY_DATABASE_MESSAGE`.
    pub filtered_empty: bool,
    pub num_threads: i32,
    /// The `-query_loc` range as given, which the tabular formats print as `qlen` and divide
    /// `qcovs` by; read by `native::run` where NCBI's query options handler reads it (after the
    /// subjects), `None` otherwise.
    pub query_range: Option<crate::blastinput::seq_range::SequenceRange>,
    /// RPSTBLASTN's search (`BlastxArgs::rpstblastn`): the run is RPSTBLASTN's, whose reports name
    /// it; `None` in every BLASTX run.
    pub(crate) rpstblastn: Option<std::sync::Arc<dyn RpsTblastnSearch>>,
    /// The database masking as the formatter resolved it (`DbMaskOption::resolve`, in
    /// `native::run` and `local::run_local_database`; series NP, M4b); `None` before, and for a
    /// search without it.
    pub db_mask: Option<crate::blastinput::db_mask::ResolvedMask>,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:185-193
/// NCBI's `eAaLookupTable` (`ELookupTableType`): the standard lookup table, of word sizes 2 to 4.
pub(crate) const AA_LOOKUP_TABLE: i32 = 3;
/// NCBI's `eCompressedAaLookupTable`: the compressed lookup table, of word sizes 5 to 7.
pub(crate) const COMPRESSED_AA_LOOKUP_TABLE: i32 = 4;

impl ResolvedOptions {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1306 (cited at
    // `ncbi_environment::ApplicationSettings::matrix_file`: `BlastFindMatrixPath`)
    /// Looks up the file that NCBI reads for `unknown_matrix`'s name with `search`
    /// (`MatrixSearch::file`) into `matrix_file`, after the options are resolved: the command
    /// line's run (`native::run`) and the library's searches (`local::read_format_options`).
    pub(crate) fn find_matrix_file(
        &mut self,
        search: &crate::blastinput::ncbi_environment::MatrixSearch,
    ) -> Result<()> {
        if self.unknown_matrix {
            self.matrix_file = search.file(&self.matrix).map_err(anyhow::Error::msg)?;
        }
        Ok(())
    }

    /// Whether the search's lookup table is the compressed one (`-word_size` 5 to 7, blastx-fast).
    pub(crate) fn compressed_lookup(&self) -> bool {
        self.lookup_type == COMPRESSED_AA_LOOKUP_TABLE
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:921-923
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:517-527
    /// Whether NCBI's score block has its Gumbel block (`sbp->gbp`): a gapped search without
    /// `-searchsp` and `OLD_FSC` (`old_fsc`). Every choice between the Gumbel block's rules and
    /// those without it (the block itself, `parameters::gumbel_block`; the length of the
    /// hit-saving cutoff and its average-length error, `parameters::subject_parameters`; the
    /// report's block) asks this.
    pub(crate) fn has_gumbel_block(&self) -> bool {
        self.gapped && !self.old_fsc
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3629-3642
impl BlastxArgs {
    /// Whether the search is ungapped with a `-matrix` that NCBI's tables do not have: the name
    /// that NCBI looks for as a file (`ResolvedOptions::unknown_matrix`, which the resolution
    /// sets for such a search; a gapped one refuses the name first).
    pub(crate) fn unknown_matrix_name(&self) -> bool {
        self.ungapped && blastx_matrix(self.matrix.as_deref().unwrap_or("BLOSUM62")).is_none()
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3629-3642
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2940-2947
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1310
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
    pub fn resolve(&self) -> Result<ResolvedOptions> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2844-2851
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:872-877
        // The filtering handler (`CFilteringArgs`, before the composition handler,
        // blastx_args.cpp:44-82): NCBI's error of `-seg` here.
        let seg = match self.seg_spec() {
            Ok(seg) => seg,
            Err(message) => return Err(SemanticOptionsError::from(message).into()),
        };
        let composition = composition_mode(&self.comp_based_stats);
        if self.ungapped && composition != 0 {
            return Err(SemanticOptionsError::from("Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a gapped search").into());
        }
        // The matrix as typed (BLOSUM62 when omitted) and its gap costs: with `-matrix`, a cost
        // not given is the matrix's BEST one.
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:255-275
        let matrix_name = self.matrix.clone().unwrap_or_else(|| "BLOSUM62".into());
        let (best_open, best_extend) = match &self.matrix {
            Some(name) => protein_gap_existence_extend_params(name).unwrap_or((0, 0)),
            None => (11, 1),
        };
        let gap_open = self.gap_open.unwrap_or(best_open);
        let gap_extend = self.gap_extend.unwrap_or(best_extend);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:913-924
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1310
        if self.threshold.is_some_and(|x| x <= 0.0) {
            return Err(SemanticOptionsError::from("Non-zero threshold required").into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1312-1348
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:608-618 (the initial
        // word options' check, after the lookup table's and before the hit saving options')
        // (`<= 0.0` is false for NaN, as in C.)
        if self.xdrop_ungap.is_some_and(|x_dropoff| x_dropoff <= 0.0) {
            return Err(SemanticOptionsError::from("x_dropoff must be greater than zero").into());
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
        if counts.hitlist_size < 1 {
            return Err(SemanticOptionsError::from("No hits are being saved").into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
        if self.evalue <= 0.0 {
            return Err(SemanticOptionsError::from(
                "expect value or cutoff score must be greater than zero",
            )
            .into());
        }
        // LOSAT's limit, for values that NCBI's checks above accept: an e-value that is not
        // finite (`+inf`, `1e999`, `-nan`; NCBI's argument reads `inf` and `nan` without a sign as
        // no number). Every word size of NCBI's, 2 to 7, runs (item SC9 of line NP).
        if !self.evalue.is_finite() {
            bail!(
                "an -evalue that is not finite ({}) is not supported by LOSAT's {}",
                self.evalue,
                self.program().upper
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:429-433
        let fast = self.task == "blastx-fast";
        let word_size = self.word_size.unwrap_or(if fast { 5 } else { 3 });
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:288-300
        // (BLASTX's query is protein for these arguments: blastx_args.cpp:74 passes
        // the negation of `kQueryIsProtein`.) A word size given on the command line above 4 is the compressed
        // table with its own threshold (19.3, 21 or 20.25), and one of 4 or less leaves the task's:
        // `-task blastx-fast -word_size 3` is the standard table with threshold 20
        // (api/blast_options_local_priv.hpp:612-619).
        let task_threshold = if fast { 20.0 } else { 12.0 };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:604-622
        // The task's threshold 12 (`BLAST_WORD_THRESHOLD_BLASTX`) is the matrix's suggestion
        // (BLOSUM62's 11 + 1 for a translated query); blastx-fast's 20 and the 19.3, 21 and 20.25
        // of `-word_size` 5, 6 and 7 are not the default, and stay.
        let threshold = match self.threshold {
            Some(threshold) => threshold,
            None if self.word_size == Some(5) => 19.3,
            None if self.word_size == Some(6) => 21.0,
            None if self.word_size == Some(7) => 20.25,
            None if fast => task_threshold,
            None => suggested_threshold(SuggestionProgram::TranslatedQuery, &matrix_name),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:485-497
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
            self.program().lower,
            self.program().upper,
        )?;
        // The counts are at most `i32::MAX` (the options are `Int4`).
        let hitlist_size = counts.hitlist_size as i32;
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:250
        // rpstblastn_app.cpp does not set it (as rpsblast's and psiblast's apps): its tabular
        // formats keep their own delimiter, the tab or the comma (`TabularSpec::of_app`;
        // `delim=` is still checked).
        let app = if self.rpstblastn.is_some() {
            crate::report::archive::options::RequestKind::Rpstblastn
        } else {
            crate::report::archive::options::RequestKind::Blastx
        };
        let tabular = crate::blastinput::app::ReportFormat::from_number(choice.number)
            .and_then(|format| crate::report::tabular::TabularSpec::of_app(format, choice, app))
            .map(crate::report::tabular::TabularSpec::intern);
        let options = ResolvedOptions {
            word_size,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:60-67
            // RPSTBLASTN's threshold is 0: its reports print none.
            threshold: if self.rpstblastn.is_some() {
                0.0
            } else {
                threshold
            },
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
            use_sw_tback: self.use_sw_tback,
            seg,
            soft_masking: self.soft_masking,
            strand: self.strand.clone(),
            query_gencode: self.query_gencode,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:123-145
            // `-xdrop_ungap`, `-xdrop_gap` and `-xdrop_gap_final` replace them
            // (blast_args.cpp:276-286).
            x_dropoff: self.xdrop_ungap.unwrap_or(7.0),
            gap_x_dropoff: self.xdrop_gap.unwrap_or(15.0),
            gap_x_dropoff_final: self.xdrop_gap_final.unwrap_or(25.0),
            gap_trigger: 22.0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
            // `-dbsize` (`CBlastOptions::SetDbLength`, blast_args.cpp:2570-2572).
            db_length: self.dbsize.unwrap_or(0),
            effective_search_space: self.searchsp.unwrap_or(0),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:922-923
            // `-searchsp` sets OLD_FSC (blast_args.cpp:303-307, `parameters::gumbel_block`).
            old_fsc: self.searchsp.is_some() || crate::blastinput::app::old_fsc_environment(),
            window_size,
            evalue: self.evalue,
            sum_stats: self.sum_stats.unwrap_or(true),
            max_intron_length: self.max_intron_length,
            max_hsps: self.max_hsps.unwrap_or(0),
            query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
            culling_limit: self.culling_limit.unwrap_or(0),
            best_hit: self.best_hit(),
            subject_besthit: self.subject_besthit,
            hitlist_size,
            num_descriptions: counts.num_descriptions as i32,
            num_alignments: counts.num_alignments as i32,
            formatting: self.formatting,
            outfmt,
            show_gis: self.show_gis,
            tabular,
            taxdb_found: false,
            subject_taxdb: None,
            filtered_empty: false,
            num_threads: self.num_threads,
            query_range: None,
            rpstblastn: self.rpstblastn.clone(),
            db_mask: None,
        };
        // The gapped X-drops that NCBI's kernels cannot size are refused before the search.
        super::parameters::check_gapped_x_drops(&options)?;
        Ok(options)
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-220
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:276-290
    /// Runs the search with the settings of the program's application start
    /// (`ApplicationStart::settings`).
    pub(crate) fn run(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    ) -> Result<()> {
        super::native::run(&self, settings)
    }

    /// The program of the run's texts (`ProgramName`).
    pub(crate) fn program(&self) -> ProgramName {
        ProgramName::of(&self.rpstblastn)
    }

    /// `-seg` as NCBI's filtering handler reads it (`app::read_seg_option`).
    pub(crate) fn seg_spec(&self) -> std::result::Result<SegSpec, &'static str> {
        crate::blastinput::app::read_seg_option(&self.seg)
    }

    /// The `-db` value as NCBI keeps it (`seqdb::db_argument_name`), when given.
    pub fn database_name(&self) -> Option<&str> {
        self.db.as_deref().map(crate::seqdb::db_argument_name)
    }
}

#[cfg(test)]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3630-3642
#[path = "args_tests.rs"]
mod tests;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3636-3639
// The native app catches this category after stream extraction. Product scope
// refusals retain the separate registered CLI-v2 boundary.
#[derive(Debug)]
pub(crate) struct SemanticOptionsError(pub std::borrow::Cow<'static, [u8]>);

impl From<&'static str> for SemanticOptionsError {
    fn from(message: &'static str) -> Self {
        Self(std::borrow::Cow::Borrowed(message.as_bytes()))
    }
}

impl From<String> for SemanticOptionsError {
    fn from(message: String) -> Self {
        Self(std::borrow::Cow::Owned(message.into_bytes()))
    }
}

/// A message's bytes (`stats::protein_options::print_matrix_message`).
impl From<Vec<u8>> for SemanticOptionsError {
    fn from(message: Vec<u8>) -> Self {
        Self(std::borrow::Cow::Owned(message))
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
/// The text of the message; a byte that is not UTF-8 reads as U+FFFD (the error keeps it).
impl std::fmt::Display for SemanticOptionsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.0))
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3638-3638
impl std::error::Error for SemanticOptionsError {}
