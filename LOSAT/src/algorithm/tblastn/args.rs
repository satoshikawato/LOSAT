//! TBLASTN's options (NCBI c++/src/algo/blast/blastinput/tblastn_args.cpp and blast_args.cpp)
//! and LOSAT's limits.
use anyhow::Result;
use clap::Args;
use std::io::Write;
use std::path::PathBuf;

use crate::blastinput::value_parsers::*;

use crate::utils::genetic_code::GeneticCode;

use super::stage_d_pipeline::LocalStageDScoring;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:54-61
/// TBLASTN's `-task`: the command line's parser, and the check of `run_local` (`TaskParser::check`).
pub(crate) const TBLASTN_TASK: TaskParser =
    task_parser("TBLASTN", &["tblastn", "tblastn-fast"], &[]);

// NCBI reference (598d8ae6): c++/src/algo/blast/api/tblastn_options.cpp:54-85
// The arguments are read with NCBI's grammar (`blastinput/value_parsers.rs`); the values
// that NCBI checks after parsing are checked by `check_options`, in NCBI's order. The
// NCBI options that LOSAT's TBLASTN does not implement are rejected by the parser
// (`cli.rs`).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct TblastnArgs {
    /// The arguments that the command line gave, for the import of a search strategy
    /// (`cli::GivenArguments`; read only with `-import_search_strategy`). The first field: the
    /// parser takes each later field's values out of the matches as it reads them
    /// (`from_arg_matches_mut`), so this one reads the matches whole.
    #[command(flatten)]
    pub(crate) given: crate::cli::GivenArguments<TblastnArgs>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427:
    // The default `-` is standard input.
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH", default_value = "-")]
    pub query: PathBuf,
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). The exclusion
    // is an error of the argument parser: LOSAT's parser reports it with its own text
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1). The database is searched by
    // `tblastn/database.rs`.
    #[arg(
        long,
        value_name = "NAME",
        conflicts_with_all = ["subject", "subject_loc"],
        help = "BLAST database name (a nucleotide database, with or without parsed Seq-ids)"
    )]
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    #[command(flatten)]
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    /// `-db_soft_mask` and `-db_hard_mask` (`blastinput::db_mask`; series NP, M4b PR 2).
    #[command(flatten)]
    pub db_mask: crate::blastinput::db_mask::DbMaskArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // No constraint and no dependency: read with `-db` and with `-subject`
    // (`stage_d_pipeline::SubjectSetNumbers`).
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
    /// The ids of the search's records under `-parse_deflines` (`run`).
    #[arg(skip)]
    pub id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    #[arg(long, value_parser = crate::blastinput::value_parsers::nonnegative_ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // `tblastn-fast`: word size 5, the compressed lookup table and threshold 20 (`resolve`).
    #[arg(long, default_value = "tblastn", value_parser = TBLASTN_TASK)]
    pub task: String,
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    // (cited at `BlastpArgs::export_search_strategy`). The search strategy file
    // (`report::archive::strategy`), created by `run` before the other options are read.
    #[arg(long, value_name = "filename", value_parser = ncbi_output_path(), help = "File name to record the search strategy used")]
    pub export_search_strategy: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3527-3538
    // The search strategy to import (`report::archive::read::import`), read by `run` before the
    // other options (`import::apply`).
    #[arg(
        long,
        value_name = "filename",
        conflicts_with = "export_search_strategy",
        help = "Search strategy to use"
    )]
    pub import_search_strategy: Option<PathBuf>,
    /// The run state of an imported search strategy (`import::apply`).
    #[arg(skip)]
    pub(crate) import: Option<crate::report::archive::read::import::ImportRun>,
    // An omitted value keeps the default of the task (10).
    #[arg(long, value_parser = tblastn_real, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3; 5 for tblastn-fast]")]
    pub word_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:258-273: when omitted, each cost keeps the selected matrix's
    // BLAST_MATRIX_BEST value.
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11, or the matrix's best value]")]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1, or the matrix's best value]")]
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1029-1056
    // PD-TLOSAN-LOCAL-GENCODE-32 adds gc.prt ID 32 only for TBLASTN.
    #[arg(long, default_value_t = 1, value_parser = tblastn_genetic_code)]
    pub db_gencode: u8,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Length of the largest intron allowed in a translated nucleotide sequence when linking multiple distinct alignments [default: 0]")]
    pub max_intron_length: Option<i32>,
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    #[arg(long, value_parser = tblastn_threshold_value, help = "Minimum word score such that the word is added to the BLAST lookup table [default: 13, or the matrix's suggestion; 19.3, 21 or 20.25 for word sizes 5, 6 and 7; 20 for tblastn-fast]")]
    pub threshold: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:220-229,280-286:
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-214,276-278:
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for ungapped extensions [default: 7]")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for preliminary gapped extensions [default: 15]")]
    pub xdrop_gap: Option<f64>,
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for final gapped alignment [default: 25]")]
    pub xdrop_gap_final: Option<f64>,
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`).
    #[arg(
        long,
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
    )]
    pub comp_based_stats: String,
    #[arg(long, default_value = "0", value_name = "SPEC", help = crate::cli::outfmt_help("0 to 16, 18 and 20", "6, 7, 10 and 20", "; with -in_pssm: 0 to 4, 6, 7 and 10"))]
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // Read by the query options handler (`check_options`).
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
        conflicts_with_all = crate::blastinput::db_filter::DB_FILTER_OPTIONS,
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`);
    // an omitted value is NCBI's default for tblastn, "12 2.2 2.5".
    #[arg(
        long,
        help = "SEG: no, yes, or \"WINDOW LOCUT HICUT\" [default: 12 2.2 2.5]"
    )]
    pub seg: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942,2549-2557
    #[arg(long)]
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: false]")]
    pub soft_masking: Option<bool>,
    #[arg(long, value_parser = ncbi_boolean, help = "Use sum statistics [default: true]")]
    pub sum_stats: Option<bool>,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40, or the matrix's suggestion]")]
    pub window_size: Option<i32>,
    // An omitted value keeps the default hit list size (500).
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]", conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    #[command(flatten)]
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-206: -qcov_hsp_perc in 0..100, -max_hsps >= 1.
    #[arg(long = "max_hsps", value_parser = blastn_count, help = "Set maximum number of HSPs per subject sequence to save for each query")]
    pub max_hsps: Option<usize>,
    #[arg(long = "qcov_hsp_perc", value_parser = tblastn_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340 (CHspFilteringArgs::ExtractAlgorithmOptions):
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    #[arg(long, value_parser = nonnegative_ncbi_integer)]
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
    // With `-culling_limit` NCBI refuses the command line, LOSAT's parser too (with its own
    // text, as for the other exclusions).
    /// Best Hit algorithm overhang value (recommended value: 0.1)
    #[arg(long, value_parser = best_hit_value, value_name = "float_value", conflicts_with = "culling_limit")]
    pub best_hit_overhang: Option<f64>,
    /// Best Hit algorithm score edge value (recommended value: 0.1)
    #[arg(long, value_parser = best_hit_value, value_name = "float_value", conflicts_with = "culling_limit")]
    pub best_hit_score_edge: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329,3347-3349
    /// Return only the best HSP for each non overlapping query region
    #[arg(long)]
    pub subject_besthit: bool,
    #[arg(long, default_value_t = 1, value_parser = blastn_count)]
    pub num_threads: usize,
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
    #[arg(long)]
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1133-1143,1241-1247
    // PSI-TBLASTN (line PS, item PS7, DECISIONS 35): the query is a PSSM. The exclusions are
    // errors of the argument parser, reported with LOSAT's text
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1; `-remote` is refused by the parser,
    // `cli.rs`). The file is read with the options (`check_options`, `in_pssm.rs`); the default
    // `-query -` is not read.
    #[arg(
        long = "in_pssm",
        value_name = "psi_chkpt_file",
        value_parser = ncbi_input_path(),
        conflicts_with_all = ["query", "query_loc"],
        help = "PSI-TBLASTN checkpoint file"
    )]
    pub in_pssm: Option<PathBuf>,
    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    #[arg(skip)]
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The run's search strategy export (`-export_search_strategy`), made by `run`
    /// (`strategy::StrategyExport::for_run`).
    #[arg(skip)]
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath` of `taxdb.bti`), asked
    /// before the search where a tabular field list looks every taxid up (the plural names of the
    /// taxonomy database, the organism report; `algorithm::tblastn::tabular`), by the command line
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:997-1056
// Stage A product decision admits 32 from gc.prt:340-347 for TBLASTN only.
/// `NStr::StringToInt` of `value` in `CArgAllowGeneticCodeInteger` is
/// `value_parsers::genetic_code_number` (LOSAT's "invalid" message for the values outside 0..=255).
pub fn tblastn_genetic_code(value: &str) -> std::result::Result<u8, String> {
    let id = crate::blastinput::value_parsers::genetic_code_number(value)
        .ok_or("invalid genetic code ID")?;
    GeneticCode::try_from_id(id).map(|_| id)
}

/// The options of a TBLASTN search as NCBI's option handlers set them.
#[derive(Debug, Clone)]
pub struct ResolvedTblastnArgs {
    pub task: String,
    pub evalue: f64,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub xdrop_gap: f64,
    pub xdrop_gap_final: f64,
    /// `-xdrop_ungap` or 7 bits (`CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults`,
    /// blast_prot_options.cpp:100-104).
    pub xdrop_ungap: f64,
    pub word_size: i32,
    /// Whether the compressed-alphabet lookup table is used (word sizes 5-7, tblastn-fast).
    pub compressed_lookup: bool,
    pub sum_stats: bool,
    pub db_gencode: u8,
    pub gapped: bool,
    pub max_intron_length: i32,
    pub seg: SegSpec,
    pub soft_masking: bool,
    pub lcase_masking: bool,
    /// The matrix name as typed (the report prints it so).
    pub matrix_name: String,
    pub threshold: f64,
    pub window_size: i32,
    /// The hit list size (`ReportCounts::hitlist_size`): the -max_target_seqs value, the
    /// larger of the numbers of descriptions and alignments of the pairwise report when
    /// one of them is given, the -num_alignments value of a tabular format, or 500.
    pub hitlist_size: usize,
    /// The numbers of descriptions and alignments of the pairwise report: the
    /// -max_target_seqs value, -num_descriptions and -num_alignments, or 500 and 250.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:263-265
    /// `-line_length`, `-sorthits` and `-sorthsps` of the pairwise report.
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    /// The -max_hsps value, or 0 (every HSP).
    pub max_hsps: usize,
    /// The -qcov_hsp_perc value, or 0 (no filter).
    pub query_cov_hsp_perc: f64,
    /// The -culling_limit value (0: no culling).
    pub culling_limit: i32,
    /// The best hit filter of `-best_hit_overhang` and `-best_hit_score_edge`, `None` without
    /// them (`TblastnArgs::best_hit`).
    pub(crate) best_hit: Option<crate::core::hspfilter_besthit::BestHitParams>,
    /// Whether -subject_besthit was given.
    pub subject_besthit: bool,
    pub composition_mode: crate::blastinput::app::CompositionMode,
    /// The `-comp_based_stats` value as typed.
    pub comp_based_stats: String,
    pub num_threads: usize,
    /// The `-query_loc` range, as given.
    pub query_range: Option<crate::blastinput::seq_range::SequenceRange>,
    /// `-in_pssm`: the query's PSSM (PSI-TBLASTN, `in_pssm.rs`).
    pub pssm: Option<std::sync::Arc<crate::algorithm::pssm_search::PsiSearch>>,
    /// The database masking as the formatter resolved it (`tblastn_app::resolve_db_mask`; series
    /// NP, M4b PR 2); `None` before, and for a search without it.
    pub db_mask: Option<crate::blastinput::db_mask::ResolvedMask>,
}

impl TblastnArgs {
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

    /// The `-db` value as NCBI keeps it (white space at its ends removed), `None` for a
    /// search of `-subject`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    pub fn database_name(&self) -> Option<&str> {
        self.db.as_deref().map(crate::seqdb::db_argument_name)
    }

    /// The options as `check_options` sets them, without its warnings.
    pub fn resolve(&self) -> Result<ResolvedTblastnArgs> {
        self.check_options(&mut std::io::sink())
    }

    /// `check_options` for the output formats of a `run_local` call (their numbers).
    pub fn check_options_for_formats(
        &self,
        outfmts: &[i32],
        diagnostics: &mut dyn Write,
    ) -> Result<ResolvedTblastnArgs> {
        self.check_options_with(Some(outfmts), diagnostics)
    }

    /// NCBI's `CBlastAppArgs::SetOptions` for tblastn after the files are opened: the task's
    /// options, each argument handler in the order of `CTblastnAppArgs` (with the warnings
    /// that they post to `diagnostics`), and `Validate`, with NCBI's errors.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
    pub fn check_options(&self, diagnostics: &mut dyn Write) -> Result<ResolvedTblastnArgs> {
        self.check_options_with(None, diagnostics)
    }

    fn check_options_with(
        &self,
        outfmts: Option<&[i32]>,
        diagnostics: &mut dyn Write,
    ) -> Result<ResolvedTblastnArgs> {
        use crate::blastinput::app::{
            formatting_handler_check, parse_comp_based_stats, parse_formatting_string,
            parse_seg_option,
        };
        use crate::blastinput::protein_options_check::{
            validate_protein_options, ProteinOptionsCheck,
        };
        use crate::stats::protein_options::{
            protein_gap_existence_extend_params, suggested_threshold, suggested_window_size,
            SuggestionProgram,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:436-447
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/tblastn_options.cpp:54-85 (the handle's
        // defaults: threshold 13, word size 3, window 40, X-drops 15 and 25, BLOSUM62
        // 11/1, e-value 10, sum statistics, composition mode 2, genetic code 1)
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:138-151
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:69-75
        // With `-in_pssm` (PSI-TBLASTN, DECISIONS 35 of line PS) the options are psiblast's handle
        // with tblastn's program, threshold 13 and genetic code: the defaults of blastp's protein
        // handle (word size 3, window 40, X-drops 15 and 25, BLOSUM62 11/1, e-value 10, the same
        // values as tblastn's) without sum statistics (the handle never calls
        // `SetSumStatisticsMode`). `-task` is checked and has no effect (the test names
        // `blastx-fast`, which is not a tblastn task).
        let psi = self.in_pssm.is_some();
        let fast = !psi && self.task == "tblastn-fast";
        let mut word_size: i32 = if fast { 5 } else { 3 };
        let mut compressed_lookup = fast;
        let mut threshold: f64 = if fast { 20.0 } else { 13.0 };
        let mut evalue = 10.0;
        let (mut gap_open, mut gap_extend) = (11, 1);
        let (mut xdrop_gap, mut xdrop_gap_final) = (15.0, 25.0);
        let mut sum_stats = !psi;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:251-322 (the
        // generic search handler)
        if let Some(value) = self.evalue {
            evalue = value;
        }
        let (matrix_gap_open, matrix_gap_extend) = self
            .matrix
            .as_deref()
            .and_then(protein_gap_existence_extend_params)
            .unwrap_or((0, 0));
        match (self.gap_open, &self.matrix) {
            (Some(value), _) => gap_open = value,
            (None, Some(_)) => gap_open = matrix_gap_open,
            (None, None) => {}
        }
        match (self.gap_extend, &self.matrix) {
            (Some(value), _) => gap_extend = value,
            (None, Some(_)) => gap_extend = matrix_gap_extend,
            (None, None) => {}
        }
        if let Some(value) = self.xdrop_gap {
            xdrop_gap = value;
        }
        if let Some(value) = self.xdrop_gap_final {
            xdrop_gap_final = value;
        }
        if let Some(value) = self.word_size {
            if value > 4 {
                compressed_lookup = true;
                threshold = 19.3;
                if value > 5 {
                    threshold = 21.0;
                }
                if value > 6 {
                    threshold = 20.25;
                }
            }
            word_size = value;
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
            // (tblastn-fast with -word_size 3 or 4 uses the plain lookup table.)
            if compressed_lookup && value <= 4 {
                compressed_lookup = false;
            }
        }
        if let Some(value) = self.sum_stats {
            sum_stats = value;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-408 (the
        // filtering handler; the default -seg of tblastn is "12 2.2 2.5")
        let soft_masking = self.soft_masking.unwrap_or(false);
        let seg = parse_seg_option(self.seg.as_deref().unwrap_or("12 2.2 2.5"))?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:586-623 (the word
        // threshold handler: the matrix's suggestion replaces only tblastn's default 13)
        let matrix_name = self
            .matrix
            .clone()
            .unwrap_or_else(|| "BLOSUM62".to_string());
        if let Some(value) = self.threshold {
            threshold = value;
        } else if threshold as i32 == 13 {
            threshold = suggested_threshold(SuggestionProgram::TranslatedSubject, &matrix_name);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:485-497
        let window_size = self
            .window_size
            .unwrap_or_else(|| suggested_window_size(&matrix_name));
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
        // The query options handler comes after the window size handler and before the
        // formatting handler (tblastn_args.cpp:44-135).
        // An imported strategy's query range stands for `-query_loc`'s (`import::apply`).
        let query_range = match self.import.as_ref().and_then(|import| import.query_range) {
            Some(range) => Some(range),
            None => crate::blastinput::seq_range::parse_optional_range(
                self.query_loc.as_deref(),
                crate::blastinput::seq_range::RangeRole::Query,
            )?,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886 and
        // 2960-2977 (the formatting handler: other programs' formats, the hit list size
        // and its warning)
        let choice = parse_formatting_string(&self.outfmt)?;
        formatting_handler_check(&choice, crate::blastinput::app::Program::Tblastn)?;
        // The hit list size and the counts of the pairwise report
        // (`blastinput::report_counts`), with the handler's warnings. The formats of a
        // `run_local` call are those of `outfmts`; the command line has the one `-outfmt`.
        let own_outfmt = [choice.number];
        let (counts, warnings) = crate::blastinput::report_counts::run_report_counts_with(
            "tblastn",
            outfmts.unwrap_or(&own_outfmt),
            self.max_target_seqs,
            self.num_descriptions,
            self.num_alignments,
            &self.formatting,
        )?;
        diagnostics.write_all(&warnings)?;
        // An imported strategy's warnings for the options it replaces
        // (`x_IssueWarningsForIgnoredOptions`, after the formatting handler's).
        if let Some(import) = &self.import {
            diagnostics.write_all(&import.option_warning_lines())?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:893-903 (the
        // composition-based statistics handler; `blastinput/app.rs`)
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:866-875
        // PSI-TBLASTN takes any value as mode 0, before the test of an ungapped search.
        let (composition_mode, _) =
            parse_comp_based_stats(&self.comp_based_stats, false, self.ungapped && !psi)?;
        let composition_mode = if psi {
            crate::blastinput::app::CompositionMode::NoCompositionBasedStats
        } else {
            composition_mode
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1340
        // The PSSM is read by the last handler of tblastn's options (`CPsiBlastArgs`, after the
        // composition-based statistics, tblastn_args.cpp:128-130), before `Validate`.
        let pssm = match self.in_pssm.as_deref() {
            Some(path) => Some(std::sync::Arc::new(super::in_pssm::read(
                path,
                &matrix_name,
                self.export_search_strategy.is_some(),
                self.db.is_some(),
                self.archive.as_ref(),
                diagnostics,
            )?)),
            None => None,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1050-1055
        // PSI-TBLASTN keeps the genetic code 1 of its handle; the value is only checked.
        let db_gencode = if psi {
            super::in_pssm::check_db_gencode(self.db_gencode)?;
            1
        } else {
            self.db_gencode
        };
        validate_protein_options(&ProteinOptionsCheck {
            gapped: !self.ungapped,
            gap_open,
            gap_extend,
            matrix_name: &matrix_name,
            threshold,
            word_size,
            compressed_lookup,
            hitlist_size: counts.hitlist_size,
            evalue,
            psiblast: false,
            psi_tblastn: psi,
            x_dropoff: self.xdrop_ungap,
        })?;
        Ok(ResolvedTblastnArgs {
            task: self.task.clone(),
            evalue,
            gap_open,
            gap_extend,
            xdrop_gap,
            xdrop_gap_final,
            xdrop_ungap: self.xdrop_ungap.unwrap_or(7.0),
            word_size,
            compressed_lookup,
            sum_stats,
            db_gencode,
            gapped: !self.ungapped,
            max_intron_length: self.max_intron_length.unwrap_or(0),
            seg,
            soft_masking,
            lcase_masking: self.lcase_masking,
            matrix_name,
            threshold,
            window_size,
            hitlist_size: counts.hitlist_size,
            num_descriptions: counts.num_descriptions,
            formatting: self.formatting,
            num_alignments: counts.num_alignments,
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:313-319: SetQueryCovHspPerc, SetMaxHspsPerSubject.
            max_hsps: self.max_hsps.unwrap_or(0),
            query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
            culling_limit: self.culling_limit_value(),
            best_hit: self.best_hit(),
            subject_besthit: self.subject_besthit,
            composition_mode,
            comp_based_stats: self.comp_based_stats.clone(),
            num_threads: self.num_threads,
            query_range,
            pssm,
            db_mask: None,
        })
    }
}

/// LOSAT's limits of the search: the options that NCBI runs and LOSAT's TBLASTN does not, refused
/// where NCBI searches the first query batch (with IDENTITY's compressed lookup table,
/// `check_identity_compressed_lookup_of`): a run that NCBI stops before (a first batch without
/// letters, a PSSM that its checks reject, the PSSM engine's error, the object manager's stop) gives
/// NCBI's error (NCBI BLAST+ 2.17.0 `-xdrop_gap 5e8` or `-ungapped`: line RG D3-1 review 4c N1, as
/// psiblast's `validate_requested_blastp_support`).
pub(super) fn check_search_limits(args: &ResolvedTblastnArgs) -> Result<()> {
    use crate::blastinput::app::CompositionMode;
    if !args.gapped {
        anyhow::bail!("{}", crate::cli::run_refusal("TBLASTN", "ungapped"));
    }
    // The gapped X-drops that NCBI's kernels cannot size are refused where the search starts
    // (`core::x_drop::check_gapped_x_dropoffs`), at the search's gapped Lambda: the matrix's, or a
    // PSSM query's `kbp_gap_psi` (`PsiSearch::gapped_karlin`, whose Lambda NCBI's X-drops read,
    // blast_parameters.c:455-458). A PSI-TBLASTN run that NCBI's object manager ends as its query
    // joins the scope never gets here (NCBI BLAST+ 2.17.0 writes that stop with any -xdrop_gap,
    // line RG D3-1 fix 4). (A name that LOSAT has no matrix for is refused below.)
    if let Ok(matrix) = args.matrix_name.parse::<ScoringMatrix>() {
        let gapped = crate::stats::lookup_protein_params(&crate::config::ProteinScoringSpec {
            matrix,
            gap_open: args.gap_open,
            gap_extend: args.gap_extend,
        });
        // The final X-drop reaches the unscaled `ALIGN_EX` only without the composition-based
        // redo: mode 0, or IDENTITY's reset (below).
        crate::core::x_drop::check_gapped_x_dropoffs(
            "TBLASTN",
            (args.xdrop_gap, args.xdrop_gap_final),
            crate::core::x_drop::gapped_x_dropoffs(
                args.xdrop_gap,
                args.xdrop_gap_final,
                args.pssm
                    .as_deref()
                    .map_or(gapped.lambda, |psi| psi.gapped_karlin(&gapped).lambda),
            ),
            (args.gap_open, args.gap_extend),
            args.composition_mode == CompositionMode::NoCompositionBasedStats
                || args.matrix_name.eq_ignore_ascii_case("identity"),
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:409,3687
    // A list whose HSPs all fail keeps the best e-value DBL_MAX, so an -evalue of DBL_MAX
    // or more (+inf, 1e999, 1.7976931348623157e308) puts an empty list into the hit list,
    // and NCBI's tblastn dies of SIGSEGV reading its first HSP (blast_hits.c:3266) when
    // the search keeps enough HSPs, and runs otherwise; LOSAT cannot tell beforehand, so
    // it rejects the value (decision D12 of docs/evidence/losat_web_e2e/AUTHORITY.md). A
    // smaller -evalue such as 1e308 runs.
    if args.evalue >= f64::MAX {
        anyhow::bail!(
            "an -evalue of DBL_MAX or more ({}), with which NCBI BLAST+'s tblastn crashes on some inputs, is not supported by LOSAT's TBLASTN",
            args.evalue
        );
    }
    Ok(())
}

/// The search's composition mode, and IDENTITY's reset (IDENTITY with a mode above 0, which NCBI
/// searches in mode 0).
fn composition(args: &ResolvedTblastnArgs) -> (BlastCompoAdjustMode, bool) {
    use crate::blastinput::app::CompositionMode;
    let identity_cbs_reset = args.matrix_name.eq_ignore_ascii_case("identity")
        && args.composition_mode != CompositionMode::NoCompositionBasedStats;
    let mode = match args.composition_mode {
        _ if identity_cbs_reset => BlastCompoAdjustMode::NoCompositionBasedStats,
        CompositionMode::NoCompositionBasedStats => BlastCompoAdjustMode::NoCompositionBasedStats,
        CompositionMode::CompositionBasedStats => BlastCompoAdjustMode::CompositionBasedStats,
        CompositionMode::CompositionMatrixAdjust => BlastCompoAdjustMode::CompositionMatrixAdjust,
        CompositionMode::CompoForceFullMatrixAdjust => BlastCompoAdjustMode::ForceFullMatrixAdjust,
    };
    (mode, identity_cbs_reset)
}

/// The options of the search: the composition mode, IDENTITY's reset and the scoring; a matrix and
/// word size that LOSAT's TBLASTN does not implement is refused (LOSAT's limits of the search
/// itself are `check_search_limits`).
pub(super) fn search_settings(args: &ResolvedTblastnArgs) -> Result<SearchSettings> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1499
    // The composition mode selects the ordinary traceback (mode 0) or the composition
    // redo, whose matrix adjustment rule modes 1 to 3 pick (`blast_adjust_scores`; item SC5 of
    // line NP).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:143-155
    // IDENTITY has no composition data: NCBI searches it in mode 0, whatever mode was asked,
    // with a warning for each query of the first query batch (`IDENTITY_CBS_RESET_MESSAGE`; the
    // later batches start in mode 0, `tblastn_app::search`). The reset is made with the score
    // block, before the HSP writer takes the preliminary hit list size
    // (`blast_aux_priv.cpp:175-229`), so that size is mode 0's too.
    let (compo_adjust_mode, identity_cbs_reset) = composition(args);
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187:
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-147:
    crate::utils::threading::validate_threads(args.num_threads)?;
    let matrix: Option<ScoringMatrix> = args.matrix_name.parse().ok();
    let threshold = args.threshold;
    let window = args.window_size;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:137 and
    // c++/src/algo/blast/core/aa_ungapped.c:214-224
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:910-936
    // Every matrix of NCBI's tables (the nine of blastp and tblastn, IDENTITY with them) and
    // every gap cost pair of its table, as `validate_protein_options` has checked them; the
    // threshold and window are the options' (the matrix's suggested values by default). The
    // standard table takes word sizes 2 to 4 and the compressed one 5 to 7; the word finder is
    // the two-hit one (a window above 0) or the one-hit one (a window of 0,
    // `XTwoHit::x_pairs_one_hit`).
    let profile_supported = (2..=7).contains(&args.word_size);
    let Some(matrix) = matrix.filter(|_| profile_supported) else {
        anyhow::bail!(
            "the matrix {} with word size {} is not supported by LOSAT's TBLASTN; it implements word sizes 2 to 7",
            args.matrix_name,
            args.word_size
        );
    };
    Ok(SearchSettings {
        compo_adjust_mode,
        identity_cbs_reset,
        scoring: LocalStageDScoring {
            matrix,
            gap_open: args.gap_open,
            gap_extend: args.gap_extend,
            word_size: args.word_size as usize,
            threshold,
            window,
            gap_xdrop_bits: args.xdrop_gap,
            final_xdrop_bits: args.xdrop_gap_final,
            word_xdrop_bits: args.xdrop_ungap,
        },
    })
}

// No NCBI counterpart: LOSAT Web's host has no file system.
/// The browser's refusal of `-in_pssm` (a PSSM query): the host has no PSSM files (line PS, item
/// PSW1). The Web's checks of the options (`check_options`) and its runs (the adapter's `run`,
/// before `run_local`) make it; the library API's records take the option (`api::Inputs::pssm`).
pub fn refuse_browser_pssm_query(args: &TblastnArgs) -> Result<()> {
    if args.in_pssm.is_some() {
        anyhow::bail!("-in_pssm (a PSSM query) is not supported in the browser yet");
    }
    Ok(())
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options (`TblastnArgs::check_options`) and LOSAT's limits.
pub fn check_options(args: &TblastnArgs) -> Result<()> {
    // The host has no database files (`run_local`).
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported in the browser yet");
    }
    refuse_browser_pssm_query(args)?;
    let resolved = args.resolve()?;
    search_settings(&resolved)?;
    check_search_limits(&resolved)?;
    // Without the inputs: IDENTITY with the compressed lookup table, which the search refuses where
    // it searches its first query batch.
    check_identity_compressed_lookup_of(&resolved)
}

/// `protein_options_check::check_identity_compressed_lookup` of a TBLASTN search.
pub(super) fn check_identity_compressed_lookup_of(args: &ResolvedTblastnArgs) -> Result<()> {
    crate::blastinput::protein_options_check::check_identity_compressed_lookup(
        &args.matrix_name,
        args.compressed_lookup,
        args.word_size,
        "tblastn",
        "TBLASTN",
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/tblastn_options.cpp:54-58
// The options handle keeps the composition mode and the scoring choices that the
// search reads; these are the ones the Stage D port uses.
/// The search settings that the options select.
pub(super) struct SearchSettings {
    /// The composition mode of the search (mode 0 for IDENTITY, `identity_cbs_reset`).
    pub(super) compo_adjust_mode: BlastCompoAdjustMode,
    /// IDENTITY with a composition mode above 0, which NCBI searches in mode 0 with a warning
    /// for each query of the first query batch (`IDENTITY_CBS_RESET_MESSAGE`).
    pub(super) identity_cbs_reset: bool,
    pub(super) scoring: LocalStageDScoring,
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "tblastn_culling_tests.rs"]
mod tblastn_culling_tests;
