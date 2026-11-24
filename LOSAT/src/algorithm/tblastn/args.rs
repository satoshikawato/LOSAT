//! TBLASTN's options (NCBI c++/src/algo/blast/blastinput/tblastn_args.cpp and blast_args.cpp)
//! and LOSAT's limits.
use anyhow::Result;
use std::io::Write;
use std::path::PathBuf;

use crate::blastinput::value_parsers::*;

use super::stage_d_pipeline::LocalStageDScoring;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:54-61
/// TBLASTN's `-task`: the command line's parser, and the check of `run_local` (`TaskParser::check`).
pub(crate) const TBLASTN_TASK: TaskParser = task_parser("TBLASTN", "tblastn");

// NCBI reference (598d8ae6): c++/src/algo/blast/api/tblastn_options.cpp:54-85
// The arguments are those that NCBI's argument reader read (`ncbi_app::args`, `from_args`); the
// values that NCBI checks after parsing are checked by `check_options`, in NCBI's order. The
// NCBI options that LOSAT's TBLASTN does not implement are refused by name after the reader
// (`ncbi_app::args::refusal`).
#[derive(Debug, Clone)]
pub struct TblastnArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427:
    // The default `-` is standard input.
    pub query: PathBuf,
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). NCBI's argument reader reports the
    // exclusion (`cli::reader`). The database is searched by `tblastn/database.rs`.
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    /// `-db_soft_mask` and `-db_hard_mask` (`blastinput::db_mask`; series NP, M4b PR 2).
    pub db_mask: crate::blastinput::db_mask::DbMaskArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // No constraint and no dependency: read with `-db` and with `-subject`
    // (`stage_d_pipeline::SubjectSetNumbers`).
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
    // `tblastn-fast`: word size 5, the compressed lookup table and threshold 20 (`resolve`).
    pub task: String,
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    // (cited at `BlastpArgs::export_search_strategy`). The search strategy file
    // (`report::archive::strategy`), created by `run` before the other options are read.
    pub export_search_strategy: Option<PathBuf>,
    // An omitted value keeps the default of the task (10).
    pub evalue: Option<f64>,
    pub word_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:258-273: when omitted, each cost keeps the selected matrix's
    // BLAST_MATRIX_BEST value.
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1029-1056
    // PD-TLOSAN-LOCAL-GENCODE-32 adds gc.prt ID 32 only for TBLASTN.
    pub db_gencode: u8,
    pub max_intron_length: Option<i32>,
    pub matrix: Option<String>,
    pub threshold: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:220-229,280-286:
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-214,276-278:
    pub xdrop_ungap: Option<f64>,
    pub xdrop_gap: Option<f64>,
    pub xdrop_gap_final: Option<f64>,
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`).
    pub comp_based_stats: String,
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // Read by the query options handler (`check_options`).
    pub query_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // Read by the database arguments handler when it reads the subjects (`run`).
    pub subject_loc: Option<String>,
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`);
    // an omitted value is NCBI's default for tblastn, "12 2.2 2.5".
    pub seg: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942,2549-2557
    pub lcase_masking: bool,
    // NCBI's flag, off by default: the Smith-Waterman redo of each match
    // (`match_redo::SmithWatermanRedo`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    pub soft_masking: Option<bool>,
    pub sum_stats: Option<bool>,
    pub window_size: Option<i32>,
    // An omitted value keeps the default hit list size (500).
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-206: -qcov_hsp_perc in 0..100, -max_hsps >= 1.
    pub max_hsps: Option<usize>,
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340 (CHspFilteringArgs::ExtractAlgorithmOptions):
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329,3347-3349
    /// Return only the best HSP for each non overlapping query region
    pub subject_besthit: bool,
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
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1133-1143,1241-1247
    // PSI-TBLASTN (line PS, item PS7, DECISIONS 35): the query is a PSSM. NCBI's argument reader reports the exclusions
    // (`cli::reader`; `-remote` is refused by name after it). The file is read with the options (`check_options`,
    // `in_pssm.rs`); the default `-query -` is not read.
    pub in_pssm: Option<PathBuf>,
    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// `SERIAL_VERIFY_DATA_WRITE` of the run's application settings when it changes what NCBI
    /// writes (`ncbi_environment::WriteVerification`), set with the settings by the command line's
    /// `run` and the library's entries; `-in_pssm`'s reader refuses a PSSM that the run would write
    /// otherwise (`in_pssm::PssmReader::write_verification`). Not a public option.
    pub write_verification: Option<crate::blastinput::ncbi_environment::WriteVerification>,
    /// The run's search strategy export (`-export_search_strategy`), made by `run` and by LOSAT
    /// Web's `web_adapter::give_strategy` (`strategy::StrategyExport::for_args`).
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath` of `taxdb.bti`), asked
    /// before the search where a tabular field list looks every taxid up (the plural names of the
    /// taxonomy database, the organism report; `algorithm::tblastn::tabular`), by the command line
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
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input.hpp:180-182
    // The position in the whole input of the first of these records, as a reader's counter
    // starts there (`api::run_batches`): a query's `Query_N` in the messages. 0 on the command
    // line, whose records are the whole input.
    pub first_query_ordinal: usize,
    /// The file that NCBI BLAST+ reads for an ungapped search's `-matrix` that its tables do not
    /// have (`ncbi_environment::ApplicationSettings::matrix_file`), found before the search by the
    /// command line and the library (`find_matrix_file`); `None` in LOSAT Web. Not a public
    /// option.
    pub matrix_file: Option<PathBuf>,
    /// LOSAT's refusal of a record of these records with the id of a query of an earlier call of
    /// the run (`query_batch::ScopeRefusal`), where the library searches each query batch as a call
    /// of its own (`api::run_batches`): the call's add of its batch gives it, so that the batch
    /// stops in the command line's form (line WQA D10 round 5, 2026-11-20). `None` on the command
    /// line. Not a public option.
    pub(crate) scope_refusal: Option<crate::blastinput::query_batch::ScopeRefusal>,
}

impl TblastnArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, tblastn's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(TblastnArgs {
            query: args.option_or_default("query")?,
            subject: args.given_option("subject")?,
            db: args.given_option("db")?,
            db_filter: crate::blastinput::db_filter::DbFilterArgs::from_args(args)?,
            db_mask: crate::blastinput::db_mask::DbMaskArgs::from_args(args)?,
            dbsize: args.given_option("dbsize")?,
            show_gis: args.flag("show_gis"),
            html: args.flag("html"),
            parse_deflines: args.flag("parse_deflines"),
            id_scope: Default::default(),
            searchsp: args.given_option("searchsp")?,
            task: args.option_or_default("task")?,
            out: args.given_option("out")?,
            export_search_strategy: args.given_option("export_search_strategy")?,
            evalue: args.given_option("evalue")?,
            word_size: args.given_option("word_size")?,
            gap_open: args.given_option("gapopen")?,
            gap_extend: args.given_option("gapextend")?,
            db_gencode: args.option_or_default("db_gencode")?,
            max_intron_length: args.given_option("max_intron_length")?,
            matrix: args.given_option("matrix")?,
            threshold: args.given_option("threshold")?,
            xdrop_ungap: args.given_option("xdrop_ungap")?,
            xdrop_gap: args.given_option("xdrop_gap")?,
            xdrop_gap_final: args.given_option("xdrop_gap_final")?,
            comp_based_stats: args.option_or_default("comp_based_stats")?,
            outfmt: args.option_or_default("outfmt")?,
            query_loc: args.given_option("query_loc")?,
            subject_loc: args.given_option("subject_loc")?,
            seg: args.given_option("seg")?,
            lcase_masking: args.flag("lcase_masking"),
            use_sw_tback: args.flag("use_sw_tback"),
            soft_masking: args.given_option("soft_masking")?,
            sum_stats: args.given_option("sum_stats")?,
            window_size: args.given_option("window_size")?,
            max_target_seqs: args.given_option("max_target_seqs")?,
            num_descriptions: args.given_option("num_descriptions")?,
            num_alignments: args.given_option("num_alignments")?,
            formatting: crate::blastinput::report_counts::FormattingArgs::from_args(args)?,
            max_hsps: args.given_option("max_hsps")?,
            qcov_hsp_perc: args.given_option("qcov_hsp_perc")?,
            culling_limit: args.given_option("culling_limit")?,
            best_hit_overhang: args.given_option("best_hit_overhang")?,
            best_hit_score_edge: args.given_option("best_hit_score_edge")?,
            subject_besthit: args.flag("subject_besthit"),
            num_threads: args.option_or_default("num_threads")?,
            mt_mode: args.given_option("mt_mode")?,
            ungapped: args.flag("ungapped"),
            in_pssm: args.given_option("in_pssm")?,
            archive: Default::default(),
            write_verification: Default::default(),
            strategy: Default::default(),
            tabular_taxdb: Default::default(),
            subject_taxdb: Default::default(),
            taxdb_setting: Default::default(),
            first_query_ordinal: Default::default(),
            matrix_file: Default::default(),
            scope_refusal: Default::default(),
        })
    }
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
    /// `-use_sw_tback` (NCBI's `SetSmithWatermanMode`, blast_args.cpp:885).
    pub use_sw_tback: bool,
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
    /// The `-html` links of the run (`blastinput::app::html_links`); `None` without `-html`.
    pub html: Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
    /// What NCBI's score block makes of an ungapped search's `-matrix` that its tables do not
    /// have (`unknown_matrix_name`); `None` for a matrix of the tables and for a gapped search,
    /// which refuses such a name first (`validate_protein_options`).
    pub unknown_matrix: Option<crate::algorithm::tblastx::args::TblastxMatrix>,
    /// The file that NCBI reads for an `Unknown` matrix (`TblastnArgs::matrix_file`).
    pub matrix_file: Option<PathBuf>,
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

    /// What NCBI's score block makes of an ungapped search's `-matrix` name that its tables do
    /// not have (`tblastx::args::ungapped_matrix`: IDENTITY without a table, or a name that it
    /// looks for as a file); `None` for a gapped search or a matrix of the tables.
    pub(crate) fn unknown_matrix_name(
        &self,
    ) -> Option<crate::algorithm::tblastx::args::TblastxMatrix> {
        use crate::algorithm::tblastx::args::{ungapped_matrix, TblastxMatrix};
        if !self.ungapped {
            return None;
        }
        match ungapped_matrix(self.matrix.as_deref().unwrap_or("BLOSUM62")) {
            TblastxMatrix::Table(_) => None,
            other => Some(other),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1306
    /// Looks up the file that NCBI reads for an `Unknown` matrix name (`unknown_matrix_name`)
    /// with `search` (`MatrixSearch::file`) into `matrix_file`: the first step of the command
    /// line's run and of the library's searches, as TBLASTX's (`TblastxArgs::find_matrix_file`);
    /// `check_search_limits` refuses a search with such a file where NCBI sets up its first batch.
    pub(crate) fn find_matrix_file(
        &mut self,
        search: &crate::blastinput::ncbi_environment::MatrixSearch,
    ) -> Result<()> {
        if self.unknown_matrix_name()
            == Some(crate::algorithm::tblastx::args::TblastxMatrix::Unknown)
        {
            self.matrix_file = search
                .file(self.matrix.as_deref().unwrap_or("BLOSUM62"))
                .map_err(anyhow::Error::msg)?;
        }
        Ok(())
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
        let query_range = crate::blastinput::seq_range::parse_optional_range(
            self.query_loc.as_deref(),
            crate::blastinput::seq_range::RangeRole::Query,
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2976 and
        // 2960-2977 (the formatting handler: other programs' formats, the hit list size
        // and its warning)
        let choice = parse_formatting_string(&self.outfmt)?;
        formatting_handler_check(&choice, crate::blastinput::app::Program::Tblastn)?;
        // The hit list size and the counts of the pairwise report
        // (`blastinput::report_counts`), with the handler's warnings. The formats of a
        // `run_local` call are those of `outfmts`; the command line has the one `-outfmt`.
        let own_outfmt = [choice.number];
        let html = crate::blastinput::app::html_links(
            self.html,
            crate::blastinput::app::Program::Tblastn,
            self.database_name(),
            outfmts.unwrap_or(&own_outfmt),
        )?;
        let (counts, warnings) = crate::blastinput::report_counts::run_report_counts_with(
            "tblastn",
            outfmts.unwrap_or(&own_outfmt),
            self.max_target_seqs,
            self.num_descriptions,
            self.num_alignments,
            &self.formatting,
        )?;
        diagnostics.write_all(&warnings)?;
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
                self.write_verification.as_ref(),
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
            phi: false,
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
            use_sw_tback: self.use_sw_tback,
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
            html,
            unknown_matrix: self.unknown_matrix_name(),
            matrix_file: self.matrix_file.clone(),
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
    // An ungapped search of a PSSM query: NCBI BLAST+ 2.17.0's tblastn dies with SIGSEGV
    // (fixture `db.tbn.inpssm.rej.ungapped`, line PS7), so no result is defined.
    if !args.gapped && args.pssm.is_some() {
        anyhow::bail!("{}", ungapped_pssm_refusal());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3103-3105,3136
    // An ungapped search's matrix name that NCBI's loader reads as IDENTITY without a table, or
    // whose file it reads: NCBI 2.17.0 aborts when it sets up the first batch (as TBLASTX's,
    // `tblastx::args::ungapped_matrix`). A name without a file fails the batch with NCBI's
    // engine error (`stage_d_pipeline::SetupEngineError::unknown_matrix`).
    if args.unknown_matrix == Some(crate::algorithm::tblastx::args::TblastxMatrix::Aborts) {
        return Err(
            crate::algorithm::tblastx::args::identity_named_matrix_refusal(
                &args.matrix_name,
                "TBLASTN",
            ),
        );
    }
    if let Some(file) = args.matrix_file.as_deref() {
        anyhow::bail!(
            "{}",
            crate::blastinput::ncbi_environment::matrix_file_refusal(
                &args.matrix_name,
                file,
                "TBLASTN"
            )
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:449-464
    // NCBI converts the gapped X-drops to raw scores only when the search has gapped Karlin
    // blocks (with the smallest gapped Lambda). An ungapped search has no gapped Karlin block (blast_setup.c:517-527), so it reads neither
    // gapped X-drop nor the gap costs.
    // The gapped X-drops that NCBI's kernels cannot size are refused where the search starts
    // (`core::x_drop::check_gapped_x_dropoffs`), at the search's gapped Lambda: the matrix's, or a
    // PSSM query's `kbp_gap_psi` (`PsiSearch::gapped_karlin`, whose Lambda NCBI's X-drops read,
    // blast_parameters.c:455-458). A PSI-TBLASTN run that NCBI's object manager ends as its query
    // joins the scope never gets here (NCBI BLAST+ 2.17.0 writes that stop with any -xdrop_gap,
    // line RG D3-1 fix 4). (A name that LOSAT has no matrix for is refused below.)
    if let Some(matrix) = args
        .matrix_name
        .parse::<ScoringMatrix>()
        .ok()
        .filter(|_| args.gapped)
    {
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
    let (compo_adjust_mode, _) = composition(args);
    // An infinite -evalue, or DBL_MAX with the composition-based redo (decision D12;
    // `protein_options_check::check_evalue_limit`). A smaller -evalue such as 1e308 runs.
    crate::blastinput::protein_options_check::check_evalue_limit(
        args.evalue,
        compo_adjust_mode.is_enabled(),
        "TBLASTN",
    )?;
    Ok(())
}

/// The refusal of `-ungapped` with `-in_pssm` (`check_search_limits`).
pub(super) fn ungapped_pssm_refusal() -> &'static str {
    "-ungapped with -in_pssm (an ungapped search of a PSSM query, which crashes NCBI BLAST+ 2.17.0's tblastn) is not supported by LOSAT's TBLASTN"
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
    // An ungapped search's name that NCBI's tables do not have scores with BLOSUM62 until the
    // setup of its first batch fails (`LocalStageDScoring::unknown_matrix`).
    let matrix: Option<ScoringMatrix> = match args.unknown_matrix {
        Some(_) => Some(ScoringMatrix::Blosum62),
        None => args.matrix_name.parse().ok(),
    };
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
            gapped: args.gapped,
            unknown_matrix: args.unknown_matrix.is_some(),
            // NCBI's traceback runs the redo when composition-based statistics are on or the
            // traceback is Smith-Waterman: `-use_sw_tback` runs the redo with or without them.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1487
            smith_waterman: args.use_sw_tback,
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
    // `-db` is checked as the command line checks it: the database is opened by the run
    // (`run_local_database`, which LOSAT Web's host gives the run's files; line WQA D12a).
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
