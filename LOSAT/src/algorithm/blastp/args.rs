//! Command-line arguments for BLASTP

use anyhow::{bail, Result};
use std::path::PathBuf;

use crate::config::{ProteinScoringSpec, ScoringMatrix};

pub use crate::blastinput::value_parsers::SegSpec as BlastpSegSpec;
use crate::blastinput::value_parsers::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:825-876
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastpCompositionMode {
    NoCompositionBasedStats,
    CompositionBasedStats,
    CompositionMatrixAdjust,
    ForceFullMatrixAdjust,
}

/// How a format's formatter cuts an app's traceback hit lists to the hit list size
/// (`BlastpApp::list_cut`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ListCut {
    /// The lists as the search made them.
    Uncut,
    /// A copy of the lists cut after the round's PSSM and inclusion set are made from them.
    AfterPssm,
    /// The round's lists cut in place before its PSSM and inclusion set are made (psiblast's XML).
    BeforePssm,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:877-886
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:44-56
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:418-421
// NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:117-120
/// The NCBI application whose search the blastp engine runs: blastp, the first round of
/// psiblast (`algorithm/psiblast`), whose program is `ePSIBlast` (`eBlastTypePsiBlast` in the
/// core) with blastp's search until a PSSM exists except blastp's restricted gapped alignment
/// (blast_parameters.c:852-855, `blast_engine::blastp_restricted_align_enabled`), or rpsblast
/// (`algorithm/rpsblast`), whose search blastp's driver hands to the RPS engine
/// (`algorithm::rps_search`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BlastpApp {
    #[default]
    Blastp,
    Psiblast,
    Rpsblast,
}

impl BlastpApp {
    /// The program's name in NCBI's messages (`Warning: [psiblast] ...`) and reports
    /// (`# PSIBLAST 2.17.0+` in upper case): its application's name
    /// (`options::RequestKind::application`, the one table of the names).
    pub fn ncbi_name(self) -> &'static str {
        self.request_kind().application()
    }

    /// The program type of the app's Seq-aligns where BLASTP does not name it
    /// (`RunRecord::program_type`): psiblast's (every round) and rpsblast's.
    pub fn program_type(self) -> Option<crate::report::record::ProgramType> {
        use crate::report::record::ProgramType;
        match self {
            BlastpApp::Blastp => None,
            BlastpApp::Psiblast => Some(ProgramType::PsiBlast),
            BlastpApp::Rpsblast => Some(ProgramType::RpsBlast),
        }
    }

    /// The app's NCBI application (`options::RequestKind`).
    pub fn request_kind(self) -> crate::report::archive::options::RequestKind {
        use crate::report::archive::options::RequestKind;
        match self {
            BlastpApp::Blastp => RequestKind::Blastp,
            BlastpApp::Psiblast => RequestKind::Psiblast,
            BlastpApp::Rpsblast => RequestKind::Rpsblast,
        }
    }

    /// Whether the app's Seq-aligns come from traceback hit lists longer than `max_target_seqs`
    /// that the formatter then cuts (`PruneSeqalign`): psiblast's larger lists, and rpsblast's
    /// merged volumes. `# N hits found` of outfmt 7 counts before the cut.
    pub fn prunes_after_traceback(self) -> bool {
        matches!(self, BlastpApp::Psiblast | BlastpApp::Rpsblast)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1066-1133
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:707-711
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2358-2362
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:873-901
    /// How the formatter of the output format `outfmt` cuts the app's hit lists: the one owner of
    /// this rule for the search (`report_prep::report_values` reads it of the first format) and
    /// for the library's grouping of formats into searches (`api::run::search_groups`), since a
    /// search's lists are cut once for all its reports. `WriteArchive` and the organism report
    /// (`x_PrintTaxReport`, outfmt 18) cut nothing; psiblast's XML formatter and its XML2 and JSON
    /// formatters (outfmt 13 to 16) cut the round's results themselves before the round's PSSM and
    /// inclusion set are made (`TrimSeqAlign`); the other formats of psiblast and rpsblast cut a
    /// copy after them (`PruneSeqalign`); the lists of the other apps are never longer than the
    /// formats show.
    pub fn list_cut(self, outfmt: i32) -> ListCut {
        use crate::blastinput::app::output_format::{
            ARCHIVE_FORMAT, JSON, TAX_FORMAT, XML, XML2_S,
        };
        match (self, outfmt) {
            _ if !self.prunes_after_traceback() => ListCut::Uncut,
            (_, ARCHIVE_FORMAT | TAX_FORMAT) => ListCut::Uncut,
            (BlastpApp::Psiblast, XML | JSON..=XML2_S) => ListCut::BeforePssm,
            _ => ListCut::AfterPssm,
        }
    }

    /// The name in LOSAT's rejections ("... is not supported by LOSAT's PSIBLAST").
    pub fn name(self) -> &'static str {
        match self {
            BlastpApp::Blastp => "BLASTP",
            BlastpApp::Psiblast => "PSIBLAST",
            BlastpApp::Rpsblast => "RPSBLAST",
        }
    }

    /// The application's name in NCBI's posts (`[psiblast]`).
    pub fn post_name(self) -> &'static str {
        match self {
            BlastpApp::Blastp => "blastp",
            BlastpApp::Psiblast => "psiblast",
            BlastpApp::Rpsblast => "rpsblast",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastpCompBasedStats {
    pub mode: BlastpCompositionMode,
    pub unified_p: bool,
}

impl BlastpCompBasedStats {
    #[inline]
    pub fn is_enabled(&self) -> bool {
        self.mode != BlastpCompositionMode::NoCompositionBasedStats
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:825-886
    pub fn to_ncbi_cli_string(self) -> String {
        let mode = match self.mode {
            BlastpCompositionMode::NoCompositionBasedStats => '0',
            BlastpCompositionMode::CompositionBasedStats => '1',
            BlastpCompositionMode::CompositionMatrixAdjust => '2',
            BlastpCompositionMode::ForceFullMatrixAdjust => '3',
        };
        if self.unified_p && self.mode != BlastpCompositionMode::NoCompositionBasedStats {
            format!("{mode}u")
        } else {
            mode.to_string()
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:838-886
#[cfg(any(losat_web_abi_v1, test))]
pub fn parse_comp_based_stats(value: &str) -> Result<BlastpCompBasedStats, String> {
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return Err("composition-based statistics option cannot be empty".to_string());
    };

    let mode = match first {
        '0' | 'F' | 'f' => BlastpCompositionMode::NoCompositionBasedStats,
        '1' => BlastpCompositionMode::CompositionBasedStats,
        '2' | 'D' | 'd' | 'T' | 't' => BlastpCompositionMode::CompositionMatrixAdjust,
        '3' => BlastpCompositionMode::ForceFullMatrixAdjust,
        _ => {
            return Err(format!(
                "invalid composition-based statistics mode '{value}'"
            ))
        }
    };

    // CLI v2 requires the complete mode token. NCBI blast_args.cpp:879-883:
    // for BLASTP with composition-based statistics enabled, NCBI accepts 'u' or 'U'
    // (it lowercases first) to enable unified p-values (`SetUnifiedP`); LOSAT accepts
    // only lowercase 'u', and only on enabled modes.
    let unified_p = match (chars.next(), chars.next()) {
        (None, None) => false,
        (Some('u'), None) if mode != BlastpCompositionMode::NoCompositionBasedStats => true,
        _ => {
            return Err(format!(
                "invalid composition-based statistics token '{value}'"
            ))
        }
    };

    Ok(BlastpCompBasedStats { mode, unified_p })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:55-83
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastpLookupTableType {
    AaLookupTable,
    CompressedAaLookupTable,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:85-148
//
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_advprot_options.cpp:48-66
// CLI v2 names/defaults: NCBI c++/src/algo/blast/blastinput/blast_args.cpp:
// 166-170,203-207,332-349,2657-2660,3158-3163; cmdline_flags.cpp:46-94.
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:44-60
/// BLASTP's `-task`: the command line's parser, and the check of `run_local` (`TaskParser::check`).
pub(crate) const BLASTP_TASK: TaskParser = task_parser("BLASTP", "blastp");

// The values are those that NCBI's argument reader read (`ncbi_app::args`, `from_args`), with
// NCBI's declared constraints only; the checks that NCBI makes after reading the
// arguments are in `check_options`, in NCBI's order.
#[derive(Debug, Clone)]
pub struct BlastpArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // `-` (the default) is standard input.
    pub query: PathBuf,
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastinput/app.rs` `missing_subject_error`).
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). NCBI's argument reader reports the
    // exclusion (`cli::reader`).
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    /// `-db_soft_mask` and `-db_hard_mask` (`blastinput::db_mask`; series NP, M4b).
    pub db_mask: crate::blastinput::db_mask::DbMaskArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // No constraint and no dependency: read with `-db` and with `-subject`.
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:225-228
    // `blastp-short` (word size 2, PAM30 9/1, `BlastpTaskOptions::create`) runs the lookup table
    // of word size 2 (item SC3).
    pub task: String,
    // An omitted value keeps the default of the task (10).
    pub evalue: Option<f64>,
    pub threshold: Option<f64>,
    pub word_size: Option<i32>,
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
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    // The search strategy file (`report::archive::strategy`), created by `run` before the other
    // options are read.
    pub export_search_strategy: Option<PathBuf>,
    // An omitted value keeps the default hit list size (500) and NCBI's numbers of
    // descriptions and alignments of the pairwise report (500, 250).
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`; psiblast's and rpsblast's from their own
    /// parsers, `into_blastp_args`). Only blastp's report sorts (`ResolvedBlastpArgs::sort`).
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    pub max_hsps_per_subject: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201,313-315: no filter when omitted or 0.
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
    // With `-culling_limit` NCBI refuses the command line, LOSAT's parser too (with its own
    // text, as for the other exclusions). psiblast's are its parser's
    // (`PsiblastArgs::best_hit_overhang`); rpsblast has neither.
    /// Best Hit algorithm overhang value (recommended value: 0.1)
    pub best_hit_overhang: Option<f64>,
    /// Best Hit algorithm score edge value (recommended value: 0.1)
    pub best_hit_score_edge: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329,3347-3349
    /// Return only the best HSP for each non overlapping query region
    pub subject_besthit: bool,
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229, 276-286
    // BLASTP's X-drops in bits; `None`: the defaults of every task (`ResolvedBlastpArgs`).
    // psiblast's are its parser's (`PsiblastArgs::xdrop_ungap`); rpsblast's are its RPS engine's
    // (`rps`).
    pub xdrop_ungap: Option<f64>,
    pub xdrop_gap: Option<f64>,
    pub xdrop_gap_final: Option<f64>,
    pub window_size: Option<i32>,
    pub matrix: Option<String>,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`).
    pub comp_based_stats: String,
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`);
    // an omitted value keeps the default of the task (no filtering).
    pub seg: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    pub soft_masking: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942,2549-2557
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // Read by the query options handler (`check_options`).
    pub query_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // Read by the database arguments handler when it reads the subjects (`run`).
    pub subject_loc: Option<String>,
    // NCBI's argument is a string (c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660), parsed before the option
    // handlers (`blastinput/app.rs` `parse_formatting_string`).
    pub outfmt: String,
    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The run's search strategy export (`-export_search_strategy`), made by `run` and by LOSAT
    /// Web's `web_adapter::give_strategy` (`strategy::StrategyExport::for_args`).
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,
    /// The recorded options of the search strategy's requests of an application whose search is
    /// BLASTP's (psiblast, rpsblast: `archive::option_inputs`), made by its command line from its
    /// own arguments as given, with `-export_search_strategy`; `None` for BLASTP's.
    pub request_inputs: Option<Box<crate::report::archive::options::OptionInputs>>,
    /// The application whose search this is: psiblast's command line sets it
    /// (`algorithm/psiblast`).
    pub app: BlastpApp,
    /// The `[BLAST] BLASTDB` setting through which a `-subject` search's tabular reports look for
    /// the taxonomy database (`common::database::subject_taxonomy`), which a search of records of
    /// the library (`api::run`) reads from the call's application settings when a report looks
    /// every taxid up; `None` elsewhere (the command line passes its settings, LOSAT Web has no
    /// files). Not a public option.
    pub taxdb_setting: Option<crate::algorithm::common::database::TaxdbSetting>,
    /// The settings by which the options' resolution looks up the file of an ungapped search's
    /// `-matrix` that NCBI's tables do not have (`ResolvedBlastpArgs::matrix_file`): the command
    /// line's own settings (`search_cli`, `search_database`) and, for a search of the library with
    /// `-subject`, the call's (`api::run`, from its application start); `None` in LOSAT Web, which
    /// has no file system. Not a public option.
    pub matrix_search: Option<crate::blastinput::ncbi_environment::MatrixSearch>,
    /// `SERIAL_VERIFY_DATA_WRITE` of the run's application settings when it changes what NCBI
    /// writes (`ncbi_environment::WriteVerification`), set with the settings by the command line's
    /// `run` and the library's entries; `-in_pssm`'s reader refuses a PSSM that the run would write
    /// otherwise (`in_pssm::PssmReader::write_verification`). Not a public option.
    pub write_verification: Option<crate::blastinput::ncbi_environment::WriteVerification>,
    /// `SERIAL_VERIFY_DATA_READ` and `SERIAL_VERIFY_DATA_GET` of the run's application settings
    /// (`ncbi_environment::ReadVerification`), set with them as `write_verification` is;
    /// `-in_pssm`'s reader refuses a PSSM that fails a check they turn off
    /// (`in_pssm::PssmReader::read_verification`). Not a public option.
    pub read_verification: crate::blastinput::ncbi_environment::ReadVerification,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1117-1119
    /// psiblast's `-inclusion_ethresh`; `None` keeps the options' default
    /// (`PSI_INCLUSION_ETHRESH`, which blastp's options have too).
    pub inclusion_ethresh: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:584-591
    /// psiblast's `-gap_trigger` in bits (`CGapTriggerArgs`); `None` keeps the options' default
    /// (`BLAST_GAP_TRIGGER_PROT`, 22 bits), the only value of blastp and rpsblast, which have no
    /// such argument.
    pub gap_trigger: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1460-1471 (cited at
    // `ResolvedBlastpArgs::sum_stats`)
    /// psiblast's `-sum_stats true` (`CGenericSearchArgs`): every round links the HSPs of each
    /// subject with sum statistics (`core::link_hsps`). blastp and rpsblast have no such argument
    /// (NCBI's rpsblast's `-sum_stats true` is refused, `rps_hooks`).
    pub sum_stats: bool,
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
    /// psiblast's PSSM options (`-save_pssm_after_last_round`, `-save_each_pssm`,
    /// `-out_ascii_pssm`, `-pseudocount`); `None` for blastp.
    pub psiblast_pssm: Option<crate::algorithm::psiblast::pssm::PsiblastPssmOptions>,
    /// rpsblast's options that its RPS engine applies (`algorithm::rps_search`); `None` for
    /// blastp and psiblast.
    pub rps: Option<crate::algorithm::rps_search::engine::RpsSearch>,
}

impl BlastpArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, blastp's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(BlastpArgs {
            query: args.option_or_default("query")?,
            subject: args.given_option("subject")?,
            db: args.given_option("db")?,
            db_filter: crate::blastinput::db_filter::DbFilterArgs::from_args_with_ipg(args)?,
            db_mask: crate::blastinput::db_mask::DbMaskArgs::from_args(args)?,
            dbsize: args.given_option("dbsize")?,
            show_gis: args.flag("show_gis"),
            html: args.flag("html"),
            parse_deflines: args.flag("parse_deflines"),
            id_scope: Default::default(),
            searchsp: args.given_option("searchsp")?,
            task: args.option_or_default("task")?,
            evalue: args.given_option("evalue")?,
            threshold: args.given_option("threshold")?,
            word_size: args.given_option("word_size")?,
            num_threads: args.option_or_default("num_threads")?,
            mt_mode: args.given_option("mt_mode")?,
            out: args.given_option("out")?,
            export_search_strategy: args.given_option("export_search_strategy")?,
            max_target_seqs: args.given_option("max_target_seqs")?,
            num_descriptions: args.given_option("num_descriptions")?,
            num_alignments: args.given_option("num_alignments")?,
            formatting: crate::blastinput::report_counts::FormattingArgs::from_args(args)?,
            max_hsps_per_subject: args.given_option("max_hsps")?,
            qcov_hsp_perc: args.given_option("qcov_hsp_perc")?,
            culling_limit: args.given_option("culling_limit")?,
            best_hit_overhang: args.given_option("best_hit_overhang")?,
            best_hit_score_edge: args.given_option("best_hit_score_edge")?,
            subject_besthit: args.flag("subject_besthit"),
            ungapped: args.flag("ungapped"),
            xdrop_ungap: args.given_option("xdrop_ungap")?,
            xdrop_gap: args.given_option("xdrop_gap")?,
            xdrop_gap_final: args.given_option("xdrop_gap_final")?,
            window_size: args.given_option("window_size")?,
            matrix: args.given_option("matrix")?,
            gap_open: args.given_option("gapopen")?,
            gap_extend: args.given_option("gapextend")?,
            comp_based_stats: args.option_or_default("comp_based_stats")?,
            seg: args.given_option("seg")?,
            soft_masking: args.given_option("soft_masking")?,
            lcase_masking: args.flag("lcase_masking"),
            use_sw_tback: args.flag("use_sw_tback"),
            query_loc: args.given_option("query_loc")?,
            subject_loc: args.given_option("subject_loc")?,
            outfmt: args.option_or_default("outfmt")?,
            archive: Default::default(),
            strategy: Default::default(),
            request_inputs: Default::default(),
            app: Default::default(),
            taxdb_setting: Default::default(),
            matrix_search: Default::default(),
            write_verification: Default::default(),
            read_verification: Default::default(),
            inclusion_ethresh: Default::default(),
            gap_trigger: Default::default(),
            sum_stats: Default::default(),
            first_query_ordinal: Default::default(),
            scope_refusal: Default::default(),
            psiblast_pssm: Default::default(),
            rps: Default::default(),
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:381-402
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:56-82
// The options of a task before the command line's values (BLAST_EXPECT_VALUE 10,
// BLAST_WORDSIZE_PROT 3, BLAST_WORD_THRESHOLD_BLASTP 11, BLOSUM62 11/1, window 40, SEG
// off for the advanced protein handle).
#[derive(Debug, Clone)]
struct BlastpTaskOptions {
    evalue: f64,
    matrix_name: String,
    gap_open: i32,
    gap_extend: i32,
    word_size: i32,
    threshold: f64,
    lookup_table_type: BlastpLookupTableType,
    window_size: i32,
    seg: BlastpSegSpec,
    chaining: bool,
}

impl BlastpTaskOptions {
    fn create(task: &str) -> Result<Self> {
        let mut options = Self {
            evalue: 10.0,
            matrix_name: "BLOSUM62".to_string(),
            gap_open: 11,
            gap_extend: 1,
            word_size: 3,
            threshold: 11.0,
            lookup_table_type: BlastpLookupTableType::AaLookupTable,
            window_size: 40,
            seg: BlastpSegSpec::No,
            chaining: false,
        };
        match task {
            "blastp" => {}
            "blastp-short" => {
                options.matrix_name = "PAM30".to_string();
                options.gap_open = 9;
                options.gap_extend = 1;
                options.evalue = 20000.0;
                options.word_size = 2;
                options.threshold = 11.0;
                options.lookup_table_type = BlastpLookupTableType::AaLookupTable;
            }
            "blastp-fast" => {
                options.word_size = 5;
                options.threshold = 20.0;
                options.lookup_table_type = BlastpLookupTableType::CompressedAaLookupTable;
                options.chaining = true;
            }
            _ => bail!("unsupported blastp task '{task}': expected one of blastp, blastp-fast, blastp-short"),
        }
        Ok(options)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
    fn set_word_size(&mut self, ws: i32) {
        self.word_size = ws;
        if self.lookup_table_type == BlastpLookupTableType::CompressedAaLookupTable && ws <= 4 {
            self.lookup_table_type = BlastpLookupTableType::AaLookupTable;
        } else if self.lookup_table_type == BlastpLookupTableType::AaLookupTable && ws > 4 {
            self.lookup_table_type = BlastpLookupTableType::CompressedAaLookupTable;
        }
    }
}

#[derive(Debug, Clone)]
pub struct ResolvedBlastpArgs {
    pub query: PathBuf,
    pub subject: PathBuf,
    /// `-parse_deflines`, for the reports of the `-subject` records whose ids it read.
    pub parse_deflines: bool,
    /// `-show_gis`, for the same reports.
    pub show_gis: bool,
    /// The `-html` links of the outfmt 0 to 4 reports (`blastinput::app::html_links`); None
    /// without `-html` or any of those formats.
    pub html: Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
    /// `-html` (NCBI's `CBlastFormat::m_IsHTML`, blast_format.cpp:90,224): the HTML form of the
    /// outfmt 0 to 4 and 18 reports (the organism report takes no links).
    pub is_html: bool,
    /// The `-db` value as NCBI keeps it (`db_name`), when the search is of a database.
    pub db: Option<String>,
    /// `-db_soft_mask` or `-db_hard_mask` as the command line keeps it (`DbMaskArgs::option`;
    /// series NP, M4b).
    pub db_mask: Option<crate::blastinput::db_mask::DbMaskOption>,
    /// The database masking as the formatter resolved it (`DbMaskOption::resolve`, in
    /// `search`); `None` before, and for a search without it.
    pub db_mask_resolved: Option<crate::blastinput::db_mask::ResolvedMask>,
    /// `-dbsize` (`eff_len_options->db_length`; 0 when it is not given).
    pub dbsize: i64,
    /// `-searchsp`: the effective search space when above 0, and at any value no Gumbel block
    /// (blast_args.cpp:303-307, blast_stat.c:921-923).
    pub searchsp: Option<i64>,
    /// No Gumbel block (finite-size correction): `-searchsp` at any value, or BLASTP,
    /// PSIBLAST and RPSBLAST with the environment variable `OLD_FSC` (`app::old_fsc_environment`).
    /// An ungapped search has none either, whatever this says (`has_gumbel`).
    pub old_fsc: bool,
    pub task: String,
    pub evalue: f64,
    pub threshold: f64,
    pub word_size: usize,
    pub lookup_table_type: BlastpLookupTableType,
    pub num_threads: usize,
    pub out: Option<PathBuf>,
    /// The hit list size (`ReportCounts::hitlist_size`): the -max_target_seqs value, the
    /// larger of the numbers of descriptions and alignments of the pairwise report when
    /// one of them is given, the -num_alignments value of a tabular format, or 500.
    pub max_target_seqs: usize,
    /// The numbers of descriptions and alignments of the pairwise report: the
    /// -max_target_seqs value, -num_descriptions and -num_alignments, or 500 and 250.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    /// The width of the pairwise alignment rows (`m_LineLength`): `-line_length`, or 60.
    pub line_length: usize,
    /// `-sorthits` and `-sorthsps` of the pairwise report (`FormattingArgs::report_sort`):
    /// blastp's; none for psiblast and rpsblast, which read them and sort nothing.
    pub sort: crate::report::hit_sort::ReportSort,
    pub max_hsps_per_subject: usize,
    /// The -qcov_hsp_perc value, or 0 (no filter).
    pub query_cov_hsp_perc: f64,
    /// The -culling_limit value (0: no culling, as NCBI's `SetCullingLimit` of 0).
    pub culling_limit: i32,
    /// The best hit filter of `-best_hit_overhang` and `-best_hit_score_edge`, `None` without
    /// them (`BlastpArgs::best_hit`).
    pub(crate) best_hit: Option<crate::core::hspfilter_besthit::BestHitParams>,
    /// Whether -subject_besthit was given.
    pub subject_besthit: bool,
    pub ungapped: bool,
    /// The X-drops in bits (`word_options->x_dropoff`, `ext_options->gap_x_dropoff` and
    /// `gap_x_dropoff_final`): `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final` or the defaults
    /// of every BLASTP task, converted to raw scores by the search (`core::x_drop`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:100-111 (no
    /// task changes them, blast_options_handle.cpp:381-402)
    pub x_drop_ungapped_bits: f64,
    pub x_drop_gapped_bits: f64,
    pub x_drop_final_bits: f64,
    /// The gap trigger in bits (`kOptions->gap_trigger`): psiblast's `-gap_trigger`, or
    /// `BLAST_GAP_TRIGGER_PROT`. Every round converts it with the query's ungapped block
    /// (`kbp_std`) into the cutoff of the ungapped alignments that a gapped extension starts from.
    pub gap_trigger_bits: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:774-784
    /// `hit_params->do_sum_stats`: psiblast's `-sum_stats true`. The search links the HSPs of each
    /// subject (`core::link_hsps`) in the preliminary stage, the traceback and the
    /// composition-based redo, and lowers the cutoff scores of the HSPs it keeps
    /// (`BlastHitSavingParametersUpdate`).
    pub sum_stats: bool,
    pub window_size: usize,
    /// The matrix as the command line names it (NCBI keeps the spelling); `scoring.matrix`
    /// is its LOSAT matrix when LOSAT has it, otherwise BLOSUM62 (and the search rejects
    /// the name, `blast_engine.rs` `validate_requested_blastp_support`).
    pub matrix_name: String,
    /// An ungapped search whose `-matrix` NCBI's tables do not have (`scoring.matrix` is then
    /// BLOSUM62, which no search scores with): NCBI's set-up of the first query batch fails with
    /// "Unknown error code -1" (`blastx/query_setup.rs` `unknown_matrix_message`,
    /// `preliminary_search_error`).
    pub unknown_matrix: bool,
    /// The file that NCBI BLAST+ reads for `unknown_matrix`'s name
    /// (`ncbi_environment::ApplicationSettings::matrix_file`), found with the settings' search of
    /// the call (`BlastpArgs::matrix_search`) where the options are resolved, so for every entry of
    /// a search alike; `None` without such a file and in LOSAT Web. The first query batch is
    /// refused then (`blastx::query_setup::matrix_error_message`).
    pub matrix_file: Option<PathBuf>,
    pub scoring: ProteinScoringSpec,
    /// The composition mode that the search runs: the `-comp_based_stats` value, or mode 0 when
    /// IDENTITY resets it (`identity_cbs_reset`; the archive records the value as given).
    pub comp_based_stats: BlastpCompBasedStats,
    /// Whether IDENTITY reset a composition mode above 0 to 0 (NCBI's `CreateScoreBlock`): the
    /// queries of the first query batch then have NCBI's warning
    /// (`report::query_warnings::IDENTITY_CBS_RESET_MESSAGE`).
    pub identity_cbs_reset: bool,
    pub seg: BlastpSegSpec,
    /// `-soft_masking` (`mask_at_hash`): the SEG intervals keep words out of the lookup table
    /// but the searched query keeps its letters (`encoding.rs` `encode_query_frame`).
    pub soft_masking: bool,
    /// `-lcase_masking`: the query's lower-case letters join the SEG intervals
    /// (`encoding.rs` `encode_protein_query_frame_lcase`), and those of a `-subject` leave its
    /// word scan (`x_soft_subject.rs`).
    pub lcase_masking: bool,
    pub use_sw_tback: bool,
    pub chaining: bool,
    /// The `-query_loc` range, as given.
    pub query_range: Option<crate::blastinput::seq_range::SequenceRange>,
    pub outfmt: String,
    /// The run's BLAST archive (outfmt 11), shared by the blocks and batches of the run.
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The run's search strategy export (`BlastpArgs::strategy`), shared as the archive is.
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,
    /// The application whose search this is (`BlastpArgs::app`).
    pub app: BlastpApp,
    /// The inclusion threshold (`psi_options->inclusion_ethresh`): psiblast's
    /// `-inclusion_ethresh`, or `PSI_INCLUSION_ETHRESH`.
    pub inclusion_ethresh: f64,
    /// `BlastpArgs::first_query_ordinal`.
    pub first_query_ordinal: usize,
    /// `BlastpArgs::scope_refusal`.
    pub(crate) scope_refusal: Option<crate::blastinput::query_batch::ScopeRefusal>,
    /// The ASCII PSSMs that psiblast writes after its round (`algorithm::psiblast::pssm`).
    pub psi_pssm: Option<crate::algorithm::psiblast::pssm::PssmRequest>,
    /// psiblast `-in_pssm`: the query given as a PSSM, read with the options
    /// (`algorithm::psiblast::in_pssm`); the search then runs the PSSM's stages
    /// (`algorithm::pssm_search`).
    pub psi_search: Option<std::sync::Arc<crate::algorithm::pssm_search::PsiSearch>>,
    /// rpsblast: the RPS engine's options (`BlastpArgs::rps`); the search then runs the RPS
    /// engine (`algorithm::rps_search`) in place of blastp's stages.
    pub rps_search: Option<std::sync::Arc<crate::algorithm::rps_search::engine::RpsSearch>>,
    /// psiblast `-num_iterations` other than 1 for FASTA queries: the rounds that the round
    /// driver runs (`algorithm::psiblast::rounds`), with the PSSM files it writes.
    pub psi_rounds: Option<std::sync::Arc<crate::algorithm::psiblast::rounds::PsiRoundPlan>>,
    /// `traceback_hitlist_size`, computed once (psiblast's reads `ADAPTIVE_CBS`).
    traceback_hitlist: usize,
    /// Whether `taxdb.bti` is on NCBI's search path for a `-subject` search whose custom tabular
    /// report looks every taxid up (`x_CheckTaxDB`; `common::database::subject_taxonomy`, asked
    /// where the command line asks, `blast_engine.rs` `open_search_database`); false otherwise,
    /// and in LOSAT Web, which has no files. A `-db` search asks its database
    /// (`BlastpDatabase::taxdb_found`).
    pub subject_taxdb_found: bool,
    /// The taxonomy database of a `-subject` search's organism report (outfmt 18), made with
    /// `subject_taxdb_found` (`common::database::SubjectTaxdb`); `None` otherwise.
    pub subject_taxdb: Option<crate::algorithm::common::database::SubjectTaxdb>,
}

impl ResolvedBlastpArgs {
    /// Whether the search has a Gumbel block (the finite-size correction): not with `-searchsp`
    /// or `OLD_FSC` (`old_fsc`), and not in an ungapped search, whose score block frees it.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:517-527 (cited at
    /// `blast_engine::call_setup::gapped_block`) and c++/src/algo/blast/core/blast_stat.c:921-923
    /// (cited at `ResolvedBlastpArgs::old_fsc`'s resolution, `BlastpArgs::check_options`)
    pub fn has_gumbel(&self) -> bool {
        !self.ungapped && !self.old_fsc
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:213-224
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:92-93
    /// The hit list size of the traceback (`hitParams->options->hitlist_size`: the heaps and the
    /// hit lists of the composition-based redo): the hit list size, or psiblast's preliminary hit
    /// list size. psiblast's formatter then keeps the first `max_target_seqs` subjects
    /// (`PruneSeqalign` with `m_HitlistSize`, blast_format.cpp:813). The preliminary size is
    /// positive (`validate_requested_blastp_support`).
    pub fn traceback_hitlist_size(&self) -> usize {
        self.traceback_hitlist
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3288-3292
    /// The heap of one query's redone matches as the composition-based redo starts.
    ///
    /// With unified P-values (BLASTP's `-comp_based_stats 1u`, `2u`, `3u`) the redo ends no match
    /// early ([`super::hsp::BlastCompoHeap::without_early_end`]): the only query whose matches NCBI's test can
    /// end is the last one of a query batch, whose result LOSAT defines by NCBI's search with a
    /// copy of it after it, where NCBI redoes every match of it.
    pub(crate) fn redone_match_heap(&self) -> super::hsp::BlastCompoHeap {
        let heap =
            super::hsp::BlastCompoHeap::new(self.traceback_hitlist_size(), self.inclusion_ethresh);
        if self.comp_based_stats.unified_p {
            heap.without_early_end()
        } else {
            heap
        }
    }

    fn compute_traceback_hitlist_size(&self) -> usize {
        match self.app {
            // RPS-BLAST's traceback keeps `hitlist_size` lists per volume (`s_RPSComputeTraceback`,
            // `Blast_HSPResultsInsertHSPList` of `hsp_list` into `results` with `hitlist_size`).
            BlastpApp::Blastp | BlastpApp::Rpsblast => self.max_target_seqs,
            BlastpApp::Psiblast => crate::algorithm::blastn::hsp::get_prelim_hitlist_size(
                self.max_target_seqs,
                self.comp_based_stats.is_enabled(),
                true,
            ),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:465-468
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:215-216
    /// The query batch size: blastp's (`BATCH_SIZE` or 10000 letters), or psiblast's 1, which
    /// reads no variable (a batch of one query, and of the queries without letters before it).
    pub fn query_batch_size(&self) -> Result<u32> {
        match self.app {
            BlastpApp::Blastp => crate::blastinput::app::query_batch_size(10000),
            BlastpApp::Psiblast => Ok(1),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:141-144
            // (`eRPSBlast` takes the default of `GetQueryBatchSize`, 10000, or `BATCH_SIZE`.)
            BlastpApp::Rpsblast => crate::blastinput::app::query_batch_size(10000),
        }
    }
}

impl BlastpArgs {
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

    /// Whether the search is ungapped with a `-matrix` that NCBI's tables do not have: the name
    /// that NCBI looks for as a file (`ResolvedBlastpArgs::unknown_matrix`, which the resolution
    /// sets for such a search; a gapped one refuses the name first).
    pub(crate) fn unknown_matrix_name(&self) -> bool {
        self.ungapped
            && self
                .matrix
                .as_deref()
                .is_some_and(|name| name.parse::<ScoringMatrix>().is_err())
    }

    /// The options as `check_options` resolves them, without NCBI's warnings.
    pub fn resolve(&self) -> Result<ResolvedBlastpArgs> {
        self.check_options(&mut std::io::sink())
    }

    /// `resolve` without opening `-in_pssm` or `-in_msa` from disk (for web ABI v2's `validate`,
    /// which checks options alone without the run's working directory).
    pub(crate) fn resolve_without_files(&self) -> Result<ResolvedBlastpArgs> {
        self.check_options_with(false, &mut std::io::sink())
    }

    /// NCBI's `CBlastAppArgs::SetOptions` for blastp after the files are opened: the task's
    /// options, each argument handler in the order of `CBlastpAppArgs` (with the warnings
    /// that they post to `diagnostics`), and `Validate`, with NCBI's errors.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
    pub fn check_options(
        &self,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<ResolvedBlastpArgs> {
        self.check_options_with(true, diagnostics)
    }

    fn check_options_with(
        &self,
        read_files: bool,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<ResolvedBlastpArgs> {
        use crate::blastinput::app::CompositionMode;
        use crate::blastinput::protein_options_check::{
            validate_protein_options, ProteinOptionsCheck,
        };
        use crate::stats::protein_options::{
            protein_gap_existence_extend_params, suggested_threshold, suggested_window_size,
            SuggestionProgram,
        };
        let mut options = BlastpTaskOptions::create(&self.task)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:251-301
        if let Some(evalue) = self.evalue {
            options.evalue = evalue;
        }
        let (matrix_gap_open, matrix_gap_extend) = self
            .matrix
            .as_deref()
            .and_then(protein_gap_existence_extend_params)
            .unwrap_or((0, 0));
        match (self.gap_open, &self.matrix) {
            (Some(gap_open), _) => options.gap_open = gap_open,
            (None, Some(_)) => options.gap_open = matrix_gap_open,
            (None, None) => {}
        }
        match (self.gap_extend, &self.matrix) {
            (Some(gap_extend), _) => options.gap_extend = gap_extend,
            (None, Some(_)) => options.gap_extend = matrix_gap_extend,
            (None, None) => {}
        }
        if let Some(word_size) = self.word_size {
            if word_size > 4 {
                options.lookup_table_type = BlastpLookupTableType::CompressedAaLookupTable;
                options.threshold = 19.3;
                if word_size > 5 {
                    options.threshold = 21.0;
                }
                if word_size > 6 {
                    options.threshold = 20.25;
                }
            }
            options.set_word_size(word_size);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-406 (the
        // filtering handler; `blastinput/app.rs` `parse_seg_option`)
        if let Some(seg) = &self.seg {
            options.seg = crate::blastinput::app::parse_seg_option(seg)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
        // The argument always has its default (`kDfltArgLookupTableMaskingOnlyProt`, "false"),
        // and no other handler of blastp sets `mask_at_hash`.
        let soft_masking = self.soft_masking.unwrap_or(false);
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:634-639
        if let Some(matrix) = &self.matrix {
            options.matrix_name = matrix.clone();
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:586-623
        if let Some(threshold) = self.threshold {
            options.threshold = threshold;
        } else if options.threshold as i32 == 11 {
            options.threshold =
                suggested_threshold(SuggestionProgram::Protein, &options.matrix_name);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:485-497
        options.window_size = match self.window_size {
            Some(window) => window,
            None => suggested_window_size(&options.matrix_name),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
        // The query options handler comes after the window size handler and before the
        // formatting handler (blastp_args.cpp:44-120).
        let query_range = crate::blastinput::seq_range::parse_optional_range(
            self.query_loc.as_deref(),
            crate::blastinput::seq_range::RangeRole::Query,
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2976 and
        // 2960-2977 (the formatting handler: other programs' formats, the hit list size
        // and its warning)
        let choice = crate::blastinput::app::parse_formatting_string(&self.outfmt)?;
        crate::blastinput::app::formatting_handler_check(
            &choice,
            crate::blastinput::app::Program::Blastp,
        )?;
        let html = crate::blastinput::app::html_links(
            self.html,
            self.app.request_kind(),
            self.db.as_deref().map(crate::seqdb::db_argument_name),
            &[choice.number],
        )?;
        // The hit list size and the counts of the pairwise report
        // (`blastinput::report_counts`), with the handler's warnings.
        let (counts, warnings) = crate::blastinput::report_counts::run_report_counts_with(
            self.app.ncbi_name(),
            &[choice.number],
            self.max_target_seqs,
            self.num_descriptions,
            self.num_alignments,
            &self.formatting,
        )?;
        diagnostics.write_all(&warnings)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:893-903 (the
        // composition-based statistics handler; `blastinput/app.rs`)
        // The unified P-values are blastp's only (`program` is `eBlastp`, blast_args.cpp:880-883).
        let (mode, unified_p) = crate::blastinput::app::parse_comp_based_stats(
            &self.comp_based_stats,
            self.app == BlastpApp::Blastp,
            self.ungapped,
        )?;
        if let Some(pssm) = &self.psiblast_pssm {
            pssm.check()?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1340
        // psiblast's `-in_pssm` is read by the same handler, after the checks of its PSSM
        // files (`algorithm::psiblast::in_pssm`); the file as read goes to the rounds' plan
        // (`in_pssm_files`), and to the search strategy's request with `-db` (`PssmReader::strategy`).
        let mut in_pssm_file = None;
        let pssm_reader = crate::algorithm::pssm_search::in_pssm::PssmReader {
            // The PSSM engine makes the scores of a PSSM of frequency ratios only with
            // the options' matrix, as tblastn's (`opts.GetMatrixName()`; NCBI BLAST+
            // 2.17.0 `-matrix BLOSUM45 -gapopen 15 -gapextend 2`: the engine's error,
            // line RG D3-1 review 4c `out_m/` m13-m15).
            matrix_name: &options.matrix_name,
            strategy: self.export_search_strategy.is_some(),
            database: self.db.is_some(),
            archive: self.archive.as_ref(),
            write_verification: self.write_verification.as_ref(),
            read_verification: Some(&self.read_verification),
            out_pssm: self
                .psiblast_pssm
                .as_ref()
                .is_some_and(|pssm| pssm.out_pssm.is_some()),
            ..crate::algorithm::pssm_search::in_pssm::PssmReader::PSIBLAST
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1319-1328
        // `-in_msa`'s PSSM is made first, with the matrix and gap costs that the generic search
        // handler has set by then (`algorithm::psiblast::in_msa`); the search of it is `-in_pssm`'s.
        let msa_read = match self
            .psiblast_pssm
            .as_ref()
            .filter(|_| read_files)
            .and_then(|pssm| pssm.in_msa.as_ref())
        {
            Some(input) => Some(crate::algorithm::psiblast::in_msa::read_in_msa(
                input,
                &crate::algorithm::psiblast::in_msa::MsaScoring {
                    matrix_name: &options.matrix_name,
                    gap_open: options.gap_open,
                    gap_extend: options.gap_extend,
                    save_ascii_pssm: self
                        .psiblast_pssm
                        .as_ref()
                        .is_some_and(|pssm| pssm.out_ascii_pssm.is_some()),
                },
                pssm_reader,
                diagnostics,
            )?),
            None => None,
        };
        let pssm_read = match msa_read {
            Some(read) => Some(read),
            None => match self
                .psiblast_pssm
                .as_ref()
                .filter(|_| read_files)
                .and_then(|pssm| pssm.in_pssm.as_deref())
            {
                Some(path) => Some(crate::algorithm::pssm_search::in_pssm::read_in_pssm(
                    path,
                    pssm_reader,
                    diagnostics,
                )?),
                None => None,
            },
        };
        let psi_search = match pssm_read {
            Some((mut search, mut file)) => {
                // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:340
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:366-369
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:557
                // The search strategy's request holds the object itself and is written in ASN.1
                // text before the rounds: that write reads each delay buffer of another format
                // into the object (`PssmWithParameters::update_delayed`), so round 1's files and
                // archive write the values read, without the messages of that read again (line
                // RG-F3, PR 2: probe AR4 K2). The request itself was taken as read
                // (`PsiSearch::request`).
                if file.valid
                    && crate::algorithm::pssm_search::pssm_request_written(
                        self.export_search_strategy.is_some(),
                        self.db.is_some(),
                    )
                {
                    file.object
                        .update_delayed(crate::report::serial::asn_input::Encoding::Text);
                }
                // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:265-270
                // psiblast's archive of each round of the PSSM holds the object that the search
                // holds (the `CPsiBlast` constructor keeps the application's reference to it): after the
                // constructor's scores, with the round's K/A written back
                // (`blastp::archive::PsiArchive`).
                if self.app == BlastpApp::Psiblast && self.archive.is_some() && file.valid {
                    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-248
                    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-396
                    // With `-out_pssm` NCBI writes the object in ASN.1 text as round 1 starts
                    // (`SavePssmToFile`), before the round's archive, which then writes the
                    // values read in every format (line RG-F3, PR 2: probe AR3 B, an XML PSSM's
                    // XML archive with `-out_pssm`).
                    let mut object = file.object.clone();
                    if self
                        .psiblast_pssm
                        .as_ref()
                        .is_some_and(|pssm| pssm.out_pssm.is_some())
                    {
                        object.update_delayed(crate::report::serial::asn_input::Encoding::Text);
                    }
                    search.archive_pssm = Some(std::sync::Arc::new(object));
                    search.archive_matrix_name_unwritable = file.matrix_name_unwritable();
                }
                in_pssm_file = Some(file);
                search.setup_messages = crate::algorithm::psiblast::in_pssm::setup_messages(
                    &search,
                    mode,
                    &options.matrix_name,
                );
                Some(std::sync::Arc::new(search))
            }
            None => None,
        };
        // A PSSM search runs modes 2 and 3 as 1 (`PsiBlastSetupScoreBlock`), which every later
        // reader of the options sees (`algorithm::psiblast::in_pssm::position_based_mode`).
        let has_pssm_query = psi_search.is_some()
            || self
                .psiblast_pssm
                .as_ref()
                .is_some_and(|pssm| pssm.in_pssm.is_some() || pssm.in_msa.is_some());
        let mode = if has_pssm_query {
            crate::algorithm::psiblast::in_pssm::position_based_mode(mode)
        } else {
            mode
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:143-155
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3043
        // IDENTITY has no composition data: NCBI searches it in mode 0 whatever mode was asked,
        // with a warning for each query of the first query batch (`identity_cbs_reset`; the reset
        // stays in the options of the later batches, `blastp_app.cpp:274`), as TBLASTN
        // (`tblastn::args`). The reset is made with the score block, before the HSP writer takes
        // the preliminary hit list size (`blast_aux_priv.cpp:175-229`), so that size is mode 0's
        // too; the unified P-values are read by the composition-based redo only, so a reset
        // search has none. The BLAST archive records the mode as given (`archive.rs` reads
        // `-comp_based_stats`).
        let identity_cbs_reset = options.matrix_name.eq_ignore_ascii_case("identity")
            && mode != CompositionMode::NoCompositionBasedStats;
        let (mode, unified_p) = if identity_cbs_reset {
            (CompositionMode::NoCompositionBasedStats, false)
        } else {
            (mode, unified_p)
        };
        // NCBI's `Validate` (`BLAST_ValidateOptions`), shared with TBLASTN.
        validate_protein_options(&ProteinOptionsCheck {
            gapped: !self.ungapped,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
            matrix_name: &options.matrix_name,
            threshold: options.threshold,
            word_size: options.word_size,
            compressed_lookup: options.lookup_table_type
                == BlastpLookupTableType::CompressedAaLookupTable,
            hitlist_size: counts.hitlist_size,
            evalue: options.evalue,
            psiblast: self.app == BlastpApp::Psiblast,
            psi_tblastn: false,
            // blastp's and psiblast's `-xdrop_ungap`, or rpsblast's (`RpsSearch`).
            x_dropoff: self
                .xdrop_ungap
                .or_else(|| self.rps.as_ref().and_then(|rps| rps.x_drop_ungapped_bits)),
        })?;
        // NCBI reference (598d8ae6): c++/src/util/tables/raw_scoremat.c:131-166 (cited at
        // `blastx::args::blastx_matrix`)
        // `validate_protein_options` has refused a matrix that NCBI's tables do not know in a
        // gapped search (blast_options.c:907-940). An ungapped search takes any name, and NCBI
        // fails at the score block of the first query batch ("BLAST engine error: Error: Unknown
        // error code -1 " for each query, exit 3; NCBI BLAST+ 2.17.0, line LIB P-S3): the search
        // reads `unknown_matrix` (`blast_engine.rs` `run_call`), as BLASTX does
        // (`blastx/query_setup.rs` `unknown_matrix_message`), so no search scores with the
        // BLOSUM62 below. NCBI 2.17.0 crashes on a name that its matrix loader reads as IDENTITY
        // without IDENTITY's statistics table (`IDENTITYX`, exit 139): LOSAT refuses it, as BLASTX
        // does (`blastx::args`, ungapped searches only: a gapped one stops at NCBI's matrix error
        // first).
        let matrix = options.matrix_name.parse::<ScoringMatrix>().ok();
        let unknown_matrix = matrix.is_none();
        if unknown_matrix {
            debug_assert!(
                self.ungapped,
                "validate_protein_options refuses an unknown matrix of a gapped search"
            );
            if options
                .matrix_name
                .get(..8)
                .is_some_and(|start| start.eq_ignore_ascii_case("identity"))
            {
                bail!(
                    "-matrix {} is read as IDENTITY without IDENTITY's statistics table (NCBI BLAST+ 2.17.0 aborts); this is not supported by LOSAT's BLASTP",
                    options.matrix_name
                );
            }
        }
        let matrix = matrix.unwrap_or(ScoringMatrix::Blosum62);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1306 (cited at
        // `ncbi_environment::ApplicationSettings::matrix_file`: `BlastFindMatrixPath`)
        // The file that NCBI reads for that name, looked up with the settings' search where every
        // entry of a search resolves the options (`matrix_search`).
        let matrix_file = match (&self.matrix_search, unknown_matrix) {
            (Some(search), true) => search
                .file(&options.matrix_name)
                .map_err(anyhow::Error::msg)?,
            _ => None,
        };
        let psi_pssm = self
            .psiblast_pssm
            .as_ref()
            .and_then(crate::algorithm::psiblast::pssm::PsiblastPssmOptions::request);
        let psi_rounds = match self.psiblast_pssm.as_ref() {
            Some(options) => {
                crate::algorithm::psiblast::rounds::PsiRoundPlan::of(options, in_pssm_file)
                    .map(std::sync::Arc::new)
            }
            None => None,
        };
        let mut resolved = ResolvedBlastpArgs {
            query: self.query.clone(),
            subject: self.subject.clone().unwrap_or_default(),
            parse_deflines: self.parse_deflines,
            show_gis: self.show_gis,
            html,
            is_html: self.html,
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
            // The name without the white space at its ends: what the reader resolves and
            // what the reports show.
            db: self
                .db
                .as_deref()
                .map(crate::seqdb::db_argument_name)
                .map(str::to_string),
            db_mask: self.db_mask.option(),
            db_mask_resolved: None,
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2570-2572
            // The options start with a database length of 0
            // (blast_prot_options.cpp:150-156, `SetDbLength` of 0).
            dbsize: self.dbsize.unwrap_or(0),
            searchsp: self.searchsp,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:922-923
            // `-searchsp` sets OLD_FSC (blast_args.cpp:303-307); every app reads the variable.
            old_fsc: self.searchsp.is_some() || crate::blastinput::app::old_fsc_environment(),
            task: self.task.clone(),
            evalue: options.evalue,
            threshold: options.threshold,
            word_size: options.word_size as usize,
            lookup_table_type: options.lookup_table_type,
            num_threads: self.num_threads,
            out: self.out.clone(),
            max_target_seqs: counts.hitlist_size,
            num_descriptions: counts.num_descriptions,
            num_alignments: counts.num_alignments,
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2939-2941
            // (Read for the output formats up to 4, the others warn, `warnings_with`.)
            line_length: self.formatting.line_length(),
            // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:246-248
            // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:502
            // (rpsblast_app.cpp:187 sets the line length only, as psiblast_app.cpp does.)
            sort: match self.app {
                BlastpApp::Blastp => self.formatting.report_sort(),
                BlastpApp::Psiblast | BlastpApp::Rpsblast => Default::default(),
            },
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:317-318
            max_hsps_per_subject: self.max_hsps_per_subject.unwrap_or(0),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:313-315
            query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
            culling_limit: self.culling_limit_value(),
            best_hit: self.best_hit(),
            subject_besthit: self.subject_besthit,
            ungapped: self.ungapped,
            x_drop_ungapped_bits: self
                .xdrop_ungap
                .unwrap_or(crate::algorithm::tblastx::constants::X_DROP_UNGAPPED_BITS),
            x_drop_gapped_bits: self.xdrop_gap.unwrap_or(f64::from(
                crate::algorithm::tblastx::constants::X_DROP_GAPPED_PRELIM,
            )),
            x_drop_final_bits: self.xdrop_gap_final.unwrap_or(f64::from(
                crate::algorithm::tblastx::constants::X_DROP_GAPPED_FINAL,
            )),
            gap_trigger_bits: self
                .gap_trigger
                .unwrap_or(crate::algorithm::tblastx::constants::GAP_TRIGGER_BIT_SCORE),
            sum_stats: self.sum_stats,
            window_size: usize::try_from(options.window_size).unwrap_or(0),
            matrix_name: options.matrix_name,
            unknown_matrix,
            matrix_file,
            scoring: ProteinScoringSpec {
                matrix,
                gap_open: options.gap_open,
                gap_extend: options.gap_extend,
            },
            comp_based_stats: BlastpCompBasedStats {
                mode: match mode {
                    CompositionMode::NoCompositionBasedStats => {
                        BlastpCompositionMode::NoCompositionBasedStats
                    }
                    CompositionMode::CompositionBasedStats => {
                        BlastpCompositionMode::CompositionBasedStats
                    }
                    CompositionMode::CompositionMatrixAdjust => {
                        BlastpCompositionMode::CompositionMatrixAdjust
                    }
                    CompositionMode::CompoForceFullMatrixAdjust => {
                        BlastpCompositionMode::ForceFullMatrixAdjust
                    }
                },
                unified_p,
            },
            identity_cbs_reset,
            seg: options.seg,
            soft_masking,
            lcase_masking: self.lcase_masking,
            use_sw_tback: self.use_sw_tback,
            chaining: options.chaining,
            query_range,
            outfmt: self.outfmt.clone(),
            archive: self.archive.clone(),
            strategy: self.strategy.clone(),
            first_query_ordinal: self.first_query_ordinal,
            scope_refusal: self.scope_refusal.clone(),
            psi_pssm,
            psi_search,
            rps_search: self.rps.clone().map(std::sync::Arc::new),
            psi_rounds,
            app: self.app,
            inclusion_ethresh: self
                .inclusion_ethresh
                .unwrap_or(super::blast_engine::PSI_INCLUSION_ETHRESH),
            traceback_hitlist: 0,
            subject_taxdb_found: false,
            subject_taxdb: None,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:60-66 (no word
        // threshold for RPS-BLAST: the reports print none, blast_format.cpp:2281-2284)
        if resolved.rps_search.is_some() {
            resolved.threshold = 0.0;
        }
        resolved.traceback_hitlist = resolved.compute_traceback_hitlist_size();
        Ok(resolved)
    }
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
