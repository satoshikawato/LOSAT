//! Command-line arguments for TBLASTX

use crate::blastinput::value_parsers::*;
use crate::utils::os_bytes::os_str_bytes;
use std::ffi::OsString;
use std::path::PathBuf;

// Command-line arguments for TBLASTX (translated DNA vs translated DNA search)
//
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastx_args.cpp:66-118 (CTblastxAppArgs argument set)
// CLI v2 names/defaults: NCBI c++/src/algo/blast/blastinput/blast_args.cpp:
// 166-170,203-207,332-349,2657-2660,3158-3163; cmdline_flags.cpp:46-94.
// The arguments are those that NCBI's argument reader read (`ncbi_app::args`, `from_args`); the
// values that NCBI checks after parsing (-seg, -threshold, -word_size, -evalue, -outfmt) are
// checked in NCBI's order by `blast_engine::check_options`.
#[derive(Debug, Clone)]
pub struct TblastxArgs {
    /// The arguments that the command line gave, for the import of a search strategy
    /// (`cli::GivenArguments`; read only with `-import_search_strategy`).
    pub(crate) given: crate::cli::GivenArguments,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427:
    // The default `-` is standard input.
    pub query: PathBuf,
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). NCBI's argument reader reports the
    // exclusion (`cli::reader`). The database is searched by `tblastx/database.rs`.
    pub db: Option<OsString>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    /// `-db_soft_mask` and `-db_hard_mask` (`blastinput::db_mask`; series NP, M4b PR 2).
    pub db_mask: crate::blastinput::db_mask::DbMaskArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // No constraint and no dependency: read with `-db` and with `-subject`
    // (`subjects::SubjectSetNumbers`).
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
    pub evalue: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:578-583: a double of at least 0; the lookup table takes its
    // (Int4) value and the report prints the double. Omitted, the matrix's suggestion
    // (`word_threshold`).
    pub threshold: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:624-629 (c++/src/algo/blast/blastinput/tblastx_args.cpp:81
    // adds it)
    pub matrix: Option<OsString>,
    // Option, not a value with a default: whether it was given is what a BLAST archive (outfmt 11)
    // records (`WordSize`, NCBI blast_args.cpp:290-292); the search uses `word_size()` (default 3).
    pub word_size: Option<i32>,
    pub num_threads: usize,

    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    // (cited at `BlastpArgs::export_search_strategy`). The search strategy file
    // (`report::archive::strategy`), created by `run` before the other options are read.
    pub export_search_strategy: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3527-3538
    // The search strategy to import (`report::archive::read::import`), read by `run` before the
    // other options (`import::apply`). NCBI's argument reader reports its exclusion of
    // `-export_search_strategy` (`cli::reader`).
    pub import_search_strategy: Option<PathBuf>,
    /// The run state of an imported search strategy (`import::apply`).
    pub(crate) import: Option<crate::report::archive::read::import::ImportRun>,
    pub query_gencode: u8,
    pub db_gencode: u8,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastx_args.cpp:74-75
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:937-944
    // A value other than 0 is NCBI's options error (`check_ncbi_options`); 0 changes nothing.
    pub max_intron_length: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2910-2927
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:219,221
    // An omitted value keeps the default hit list size (500) but the pairwise report then
    // shows 250 alignments, so the option has no default (`Args::given_option`).
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    // The hit list size and the numbers of descriptions and alignments of the pairwise
    // report, which `check_ncbi_options` sets from `-max_target_seqs`, `-num_descriptions`,
    // `-num_alignments` and the output format (`blastinput::report_counts`). Not a public
    // option.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
    pub report_counts: crate::blastinput::report_counts::ReportCounts,
    // The `-html` links of the run, which `check_ncbi_options` makes
    // (`blastinput::app::html_links`). Not a public option.
    pub html_links: Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-207
    pub max_hsps: Option<usize>,
    pub qcov_hsp_perc: Option<f64>,
    // NCBI low-complexity filtering selection:
    // - dust is used only for blastn (and mapping)
    // - otherwise seg is used
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:572-580
    //
    // Therefore, for tblastx we do NOT apply nucleotide-level DUST masking.

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-406
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`).
    pub seg: OsString,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // (`CFilteringArgs` of the negation of `kQueryIsProtein`, tblastx_args.cpp:78, takes this
    // protein branch.) Read by `soft_masking`.
    pub soft_masking: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-244
    // (`CGenericSearchArgs` of the negation of `kQueryIsProtein`, false, false and true,
    // tblastx_args.cpp:68, does not suppress it.) Option, not a value with a default: whether it
    // was given is what a BLAST archive (outfmt 11) records (`SumStatistics`); the search reads
    // `sum_stats()`.
    pub sum_stats: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942,2549-2557
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1951-1957 (the
    // query strand of the programs whose query is nucleotide: blastn, blastx, tblastx)
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:128 (the query reader takes it)
    // The search translates the query's frames of the strands it names (`QueryStrand`,
    // `blast_engine/batch_setup.rs` `query_frames`); the subject keeps its six frames.
    pub strand: String,

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-214, 276-278
    // TBLASTX has no gapped X-drops (blast_args.cpp:216-229, gated on `m_IsTblastx` being false). Read by
    // `x_drop_ungapped_bits`.
    pub xdrop_ungap: Option<f64>,

    // The window of the two-hit word finder (`window`); 0 takes the one-hit word finder. The help is
    // NCBI's, as tblastn's. The command line takes 0 or more (NCBI's argument reader checks it); an
    // imported search strategy carries any `Int4`, which NCBI runs: a negative window takes the
    // one-hit word finder too, on a table that keeps its state across the batch's calls
    // (`blast_engine::word_finder::WideOffset`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:479-480
    pub window_size: Option<i32>,

    // Read when NCBI reads it (`blastinput/app.rs` `parse_formatting_string`).
    pub outfmt: OsString,

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    // Option, not a value with a default: whether it was given (also 0) is what a BLAST archive
    // (outfmt 11) records (`Culling`); the search uses `culling_limit()`.
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // Read by the query options handler (`check_ncbi_options`).
    pub query_loc: Option<OsString>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // Read by the database arguments handler when it reads the subjects (`run`).
    pub subject_loc: Option<OsString>,
    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The run's search strategy export (`-export_search_strategy`), made by `run` and by LOSAT
    /// Web's `web_adapter::give_strategy` (`strategy::StrategyExport::for_args`).
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath` of `taxdb.bti`), asked
    /// before the search where a tabular field list looks every taxid up (the plural names of the
    /// taxonomy database, the organism report; `algorithm::tblastx::tabular`), by the command line
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
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-264
    // `api::run_batches` hands over one of NCBI's query batches at a time, sized as the running
    // application's loop sizes them (BLASTN's mixer where blastn runs a TBLASTX strategy); these
    // records are that batch, whatever its size. No argument of the command line sets this.
    pub(crate) one_query_batch: bool,
    /// LOSAT's refusal of a record of these records with the id of a query of an earlier call of
    /// the run (`query_batch::ScopeRefusal`), where the library searches each query batch as a call
    /// of its own (`api::run_batches`): the call's add of its batch gives it, so that the batch
    /// stops in the command line's form (line WQA D10 round 5, 2026-11-20). `None` on the command
    /// line. Not a public option.
    pub(crate) scope_refusal: Option<crate::blastinput::query_batch::ScopeRefusal>,
    /// The file that NCBI BLAST+ reads for a `-matrix` that its tables do not have
    /// (`ncbi_environment::ApplicationSettings::matrix_file`), found before the search by the
    /// command line and the library API (`find_matrix_file`); `None` in LOSAT Web. Not a public
    /// option.
    pub matrix_file: Option<std::path::PathBuf>,
    // The command line of the run as the SAM `@PG` line shows it: every argument followed by one
    // space, raw bytes (`GetCmdlineArgs`); `None` where there is no command line (the library, the
    // web ABIs). Set by `run` from the process's arguments; only a TBLASTX strategy that BLASTN's
    // application runs writes SAM (`tblastx::sam`).
    //
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:877-885
    pub command_line: Option<Vec<u8>>,
}

impl TblastxArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, tblastx's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(TblastxArgs {
            given: crate::cli::GivenArguments::from_args(args),
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
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:157-158
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1454
            // NCBI's XML gives `-evalue` no default (its search sets `BLAST_EXPECT_VALUE`, 10,
            // where the option is not given): the one default that LOSAT writes itself.
            evalue: args.given_option("evalue")?.unwrap_or(10.0),
            threshold: args.given_option("threshold")?,
            matrix: args.given_option("matrix")?,
            word_size: args.given_option("word_size")?,
            num_threads: args.option_or_default("num_threads")?,
            out: args.given_option("out")?,
            export_search_strategy: args.given_option("export_search_strategy")?,
            import_search_strategy: args.given_option("import_search_strategy")?,
            import: Default::default(),
            one_query_batch: Default::default(),
            command_line: Default::default(),
            query_gencode: args.option_or_default("query_gencode")?,
            db_gencode: args.option_or_default("db_gencode")?,
            max_intron_length: args.given_option("max_intron_length")?,
            max_target_seqs: args.given_option("max_target_seqs")?,
            num_descriptions: args.given_option("num_descriptions")?,
            num_alignments: args.given_option("num_alignments")?,
            formatting: crate::blastinput::report_counts::FormattingArgs::from_args(args)?,
            report_counts: Default::default(),
            html_links: Default::default(),
            max_hsps: args.given_option("max_hsps")?,
            qcov_hsp_perc: args.given_option("qcov_hsp_perc")?,
            seg: args.option_or_default("seg")?,
            soft_masking: args.given_option("soft_masking")?,
            sum_stats: args.given_option("sum_stats")?,
            lcase_masking: args.flag("lcase_masking"),
            strand: args.option_or_default("strand")?,
            xdrop_ungap: args.given_option("xdrop_ungap")?,
            window_size: args.given_option("window_size")?,
            outfmt: args.option_or_default("outfmt")?,
            culling_limit: args.given_option("culling_limit")?,
            best_hit_overhang: args.given_option("best_hit_overhang")?,
            best_hit_score_edge: args.given_option("best_hit_score_edge")?,
            subject_besthit: args.flag("subject_besthit"),
            query_loc: args.given_option("query_loc")?,
            subject_loc: args.given_option("subject_loc")?,
            archive: Default::default(),
            strategy: Default::default(),
            tabular_taxdb: Default::default(),
            subject_taxdb: Default::default(),
            taxdb_setting: Default::default(),
            first_query_ordinal: Default::default(),
            scope_refusal: Default::default(),
            matrix_file: Default::default(),
        })
    }
}

/// What NCBI's score block makes of the `-matrix` name of an ungapped search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TblastxMatrix {
    /// A matrix of NCBI's tables.
    Table(crate::config::ScoringMatrix),
    /// A name that NCBI's matrices and tables do not know: NCBI's search stops with "Unknown
    /// error code -1" (`check_matrix`).
    Unknown,
    /// A name that NCBI's matrix loader reads as IDENTITY but its tables do not know: NCBI
    /// 2.17.0 aborts (`free(): invalid size`), LOSAT refuses the name.
    Aborts,
}

// NCBI reference (598d8ae6): c++/src/util/tables/raw_scoremat.c:157-161
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3103-3105,3136
/// What NCBI's score block makes of the `-matrix` name `name` of an ungapped search (TBLASTX's,
/// and TBLASTN's with `-ungapped`; `TblastxArgs::scoring_matrix` tells the rules): a matrix of the
/// tables, a name that the matrix loader reads as IDENTITY without a table (NCBI aborts), or a
/// name that it does not know (looked for as a file).
pub(crate) fn ungapped_matrix(name: impl AsRef<[u8]>) -> TblastxMatrix {
    let name = name.as_ref();
    if let Ok(matrix) = std::str::from_utf8(name)
        .ok()
        .and_then(|s| s.parse::<crate::config::ScoringMatrix>().ok())
        .ok_or(())
    {
        return TblastxMatrix::Table(matrix);
    }
    let identity = name
        .get(..8)
        .is_some_and(|start| start.eq_ignore_ascii_case(b"identity"));
    if identity {
        TblastxMatrix::Aborts
    } else {
        TblastxMatrix::Unknown
    }
}

/// LOSAT's refusal of a search with a batch whose `-matrix` name `name` the matrix loader reads
/// as IDENTITY without IDENTITY's statistics table (`TblastxMatrix::Aborts`), where NCBI BLAST+
/// 2.17.0 aborts when it sets up the batch; `program` is the program's name in capitals.
pub(crate) fn identity_named_matrix_refusal(
    name: impl AsRef<[u8]>,
    program: &str,
) -> anyhow::Error {
    anyhow::anyhow!(
        "-matrix {} is read as IDENTITY without IDENTITY's statistics table (NCBI BLAST+ 2.17.0 aborts); this is not supported by LOSAT's {program}",
        String::from_utf8_lossy(name.as_ref())
    )
}

impl TblastxArgs {
    /// The running application: the program of the command line, the library call or the Web run,
    /// which is not this program's when it runs another program's strategy
    /// (`ImportRun::app`). Its name prefixes the run's diagnostics and its formatting handler
    /// checks the formats.
    pub(crate) fn app(&self) -> crate::blastinput::app::Program {
        crate::report::archive::read::import::ImportRun::running(
            self.import.as_ref(),
            crate::blastinput::app::Program::Tblastx,
        )
    }

    /// The `-matrix` name as typed, BLOSUM62 when omitted (the pairwise report prints it so).
    pub fn matrix_name(&self) -> std::borrow::Cow<'_, str> {
        match self.matrix.as_deref() {
            Some(name) => name.to_string_lossy(),
            None => std::borrow::Cow::Borrowed("BLOSUM62"),
        }
    }

    /// Raw bytes of `-matrix` (BLOSUM62 when omitted).
    pub fn matrix_name_bytes(&self) -> std::borrow::Cow<'_, [u8]> {
        match self.matrix.as_deref() {
            Some(name) => os_str_bytes(name),
            None => std::borrow::Cow::Borrowed(b"BLOSUM62"),
        }
    }

    /// The matrix of `-matrix`.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/tables/raw_scoremat.c:157-161 (the matrix that the score
    /// block loads, `BlastScoreBlkProteinMatrixLoad`, for the upper-cased name; blast_setup.c:385)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3103-3105, 3136 (the tables of
    /// the length adjustment, found by `strcasecmp`; a name without a table leaves the
    /// arrays unset, and `sfree` of them aborts)
    /// The other names that the loader knows (`blosum` or `pam` and the exact digits) are the
    /// tables' names ignoring the case. A name that the loader does not know is looked for as
    /// a file (`BlastFindMatrixPath`, blast_stat.c:1631-1660); LOSAT reads no matrix file, and
    /// NCBI 2.17.0 aborts for a file that it finds (`SURVEY_scoring.md` §1.3), which LOSAT
    /// refuses (`matrix_file`).
    pub fn scoring_matrix(&self) -> TblastxMatrix {
        ungapped_matrix(self.matrix_name_bytes())
    }

    /// Looks up the file that NCBI reads for an `Unknown` matrix (`scoring_matrix`) with `search`
    /// (`MatrixSearch::file`) into `matrix_file`: the first step of the command line's run
    /// (`run`) and of the library's searches (`run_local_database`, `api::run`'s search of
    /// records), before the options and the outputs, as NCBI's score block, which looks for it
    /// when it sets up a batch, refuses the run there (`run_in_pool`).
    pub(crate) fn find_matrix_file(
        &mut self,
        search: &crate::blastinput::ncbi_environment::MatrixSearch,
    ) -> anyhow::Result<()> {
        if self.scoring_matrix() == TblastxMatrix::Unknown {
            self.matrix_file = search
                .file(self.matrix_name_bytes())
                .map_err(anyhow::Error::msg)?;
        }
        Ok(())
    }

    /// The word threshold: `-threshold`, or the matrix's suggestion for a translated subject; a
    /// search strategy without one keeps the handle's 13 (`BLAST_WORD_THRESHOLD_TBLASTX`): after
    /// an import NCBI runs only the handlers that `SetOptionsForSavedStrategy` names.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:586-623 (TBLASTX's default
    /// threshold is always the suggestion's)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3801-3837
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:115
    pub fn word_threshold(&self) -> f64 {
        self.threshold.unwrap_or_else(|| {
            if self.import.is_some() {
                return 13.0;
            }
            crate::stats::protein_options::suggested_threshold(
                crate::stats::protein_options::SuggestionProgram::TranslatedSubject,
                self.matrix_name_bytes(),
            )
        })
    }

    /// The two-hit window: `-window_size`, or the matrix's suggestion; a search strategy without
    /// one keeps the handle's 40 (`BLAST_WINDOW_SIZE_PROT`, as `word_threshold`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:485-497
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:57
    pub fn window(&self) -> i32 {
        self.window_size.unwrap_or_else(|| {
            if self.import.is_some() {
                return 40;
            }
            crate::stats::protein_options::suggested_window_size(self.matrix_name_bytes())
        })
    }

    /// The ungapped X-drop in bits (`word_options->x_dropoff`): `-xdrop_ungap`, else TBLASTX's
    /// default.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:100-104 (the
    /// base handle of `CTBlastxOptionsHandle`, tblastx_options.cpp:47-49)
    pub fn x_drop_ungapped_bits(&self) -> f64 {
        self.xdrop_ungap
            .unwrap_or(super::constants::X_DROP_UNGAPPED_BITS)
    }

    /// The word size the search uses: the argument, else TBLASTX's default (3).
    pub fn word_size(&self) -> i32 {
        self.word_size.unwrap_or(3)
    }

    /// The culling limit the search uses: the argument, else 0 (no culling).
    pub fn culling_limit(&self) -> i32 {
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

    /// The -seg value as NCBI's filtering handler reads it (`check_options` reports its
    /// errors first).
    pub fn seg_spec(&self) -> anyhow::Result<SegSpec> {
        crate::blastinput::app::parse_seg_option(os_str_bytes(&self.seg))
    }

    /// `-soft_masking` (`mask_at_hash`): the SEG intervals keep words out of the lookup table but
    /// the searched frames keep their letters (`x_query_frames.rs` `mask_soft`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    /// The argument always has its default (`kDfltArgLookupTableMaskingOnlyProt`, "false"), and no
    /// other handler of tblastx sets `mask_at_hash`.
    pub fn soft_masking(&self) -> bool {
        self.soft_masking.unwrap_or(false)
    }

    /// Whether the search links its HSPs with sum statistics (`do_sum_stats`): `-sum_stats`, else
    /// TBLASTX's default, true.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:321-323
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/tblastx_options.cpp:76-81
    pub fn sum_stats(&self) -> bool {
        self.sum_stats.unwrap_or(true)
    }

    /// The `-db` value as NCBI keeps it (white space at its ends removed), `None` for a
    /// search of `-subject`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    pub fn database_name(&self) -> Option<Vec<u8>> {
        self.db
            .as_deref()
            .map(|db| crate::seqdb::db_argument_name_bytes(&os_str_bytes(db)).to_vec())
    }

    /// The subject file name as the reports show it.
    pub fn subject_label(&self) -> std::path::Display<'_> {
        self.subject
            .as_deref()
            .unwrap_or(std::path::Path::new(""))
            .display()
    }
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod args_tests;
