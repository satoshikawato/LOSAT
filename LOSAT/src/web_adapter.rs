//! The engine items that LOSAT Web's adapter (`web/adapter`, ABI v2: `docs/web/abi_v2.md`)
//! names, and nothing else. No NCBI counterpart: this is the boundary between this crate and
//! the adapter crate.
//!
//! Hidden from the documentation and outside the compatibility promise of the version number,
//! as the hidden modules are; `api` is the public surface. An item here can move to `api` later,
//! one at a time.
//!
//! Rules:
//! - The adapter names the engine through `losat::api` (public items) and this module only.
//! - Each item is named on its own (no glob, no module), so a new item shows in the pull
//!   request's diff of this file, and every item here is one the adapter uses.
//! - The adapter is code outside the crate: a function here must be sound for every value a
//!   safe caller can pass (no out-of-range read or write behind `unsafe`, whatever the input).
//!   A function whose soundness needs a narrower input gets a wrapper that fixes it
//!   ([`try_parse_cli`]).
//!
//! `web/adapter/tests/engine_boundary.rs` checks the first two rules; the third is the review's.

use std::ffi::OsString;

// The command line's arguments (`cli.rs`): the clap types of every program, the parser's
// messages, the text that the command line prints where its reader stopped at a help word, and
// the command line's own errors.
pub use crate::cli::{
    command_line_message, printed_text, Cli, Commands, NativeError, WEB_V2_PROGRAMS,
};

/// `cli::try_parse_with_printing_start` for the command line's own arguments (`Cli`): the argv read
/// as the command line reads it after its application start, NCBI's single-dash options included.
/// A run of ABI v2 makes its own start (`start_run`), so the toolkit's words that the start reads
/// (`-logfile`, `-conffile`, `-dryrun`) are taken, as the library API takes them, and the help
/// words stop the reader as on the command line, `-xmlhelp` included: the run prints their text
/// (`printed_text`; line WQA's D14, 2026-11-20). The versions are printed by the start
/// (`CallerStart::version`) and never reach the reader.
///
/// Only `Cli`: the adapter parses nothing else, so the façade does not offer the parser for any
/// other clap type.
pub fn try_parse_cli<I, S>(argv: I) -> Result<Cli, clap::Error>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    crate::cli::try_parse_with_printing_start(argv)
}

// The NCBI application start of a run (`api::start`, line RG item E1w-3a): the words that it reads
// before the arguments (`preparse`, NCBI's scan of `-logfile`), the registry files that the host
// registered (`RegistryFiles`), and where the run's files go (`RunFiles`: the file of `-logfile`).
pub use crate::api::start::{CallerSink, CallerStart, ErrorSource};
pub use crate::blastinput::ncbi_environment::{
    log_file, preparse, ApplicationSettings, Preparsed, RegistryFiles, TextPlace,
};
pub use crate::report::run_files::{FileReceiver, RunFiles};

/// The NCBI application start of a run of `program` with the arguments `words` (after the
/// program's name), as the command line's (`cli::main_entry`) in an environment of the
/// `variables` alone, with the registry files of `files` (`ApplicationStart::with_files`) and the
/// file of `-logfile` going to `run_files` (`api::start::begin_with`): the messages that go to
/// `cerr` are written to `diagnostics`.
pub fn start_run(
    program: &str,
    words: &[&str],
    variables: &[(&str, &str)],
    files: &dyn RegistryFiles,
    run_files: &RunFiles,
    diagnostics: &mut dyn std::io::Write,
) -> Result<CallerStart, String> {
    let start = crate::blastinput::ncbi_environment::ApplicationStart::with_files(
        program, words, variables, files,
    );
    crate::api::start::begin_with(start, words, run_files, diagnostics)
}

// The checks of the options alone, without the inputs (the `validate` of ABI v2), and BLASTN's
// effective values of the options (the defaults that `describe` shows).
pub use crate::algorithm::blastn::coordination::{
    configure_task, determine_effective_word_size, determine_evalue, determine_scoring_params,
    determine_window_size,
};
pub use crate::algorithm::blastn::scoring::{check_losat_limits, check_scoring_options};
pub use crate::algorithm::blastn::BlastnArgs;
pub use crate::algorithm::blastp::blast_engine::check_options as check_blastp_options;
pub use crate::algorithm::blastx::check_options as check_blastx_options;
pub use crate::algorithm::blastx::BlastxArgs;
pub use crate::algorithm::tblastn::check_options as check_tblastn_options;
pub use crate::algorithm::tblastn::refuse_browser_pssm_query as refuse_tblastn_pssm_query;
pub use crate::algorithm::tblastx::blast_engine::check_options as check_tblastx_options;
pub use crate::utils::threading::validate_threads;

// The searches: each program's `run_local` over records in memory, the pool that the searches
// of one argv share, and the routing of their outputs (`api/local_blast.rs`).
pub use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
    with_run_pool, FormatObserver, FormatOutput, FormatRefusal, HspIndex, OutputSink,
    ReportOutputs,
};
pub use crate::report::blast_xml2::SideFileSink;
pub use crate::report::PairwiseHit;

// LOSAT's words of NCBI's messages (`wording.rs`, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 12):
// the adapter's tests compare the engine's messages with NCBI's texts through it (AGENTS.md
// rule 8).
pub use crate::wording::{
    losat_bytes, losat_output_bytes, losat_text, LOSAT_VERSION, REPORT_VERSION,
};

// The output formats of the application layer (`blastinput/app.rs`, NCBI's CFormattingArgs):
// `-outfmt` parsed and checked as the command line does. `EngineProgram` is the program of
// these checks (the adapter has its own `Program`, and `api::Program` is another type).
pub use crate::blastinput::app::output_format::{ARCHIVE_FORMAT, SAM as SAM_FORMAT};
pub use crate::blastinput::app::Program as EngineProgram;
pub use crate::blastinput::app::{
    engine_error, formatting_handler_check, missing_subject_error, parse_formatting_string,
    report_format,
};

// The hit list sizes of the formats (`blastinput/report_counts.rs`), the scoring matrices that
// `describe` lists, and the `-query_loc`/`-subject_loc` ranges (`blastinput/seq_range.rs`).
pub use crate::blastinput::report_counts::{report_counts, run_report_counts};
pub use crate::blastinput::seq_range::{
    parse_optional_range, record_interval, RangeRole, RecordInterval,
};
pub use crate::config::ScoringMatrix;

// The FASTA reader (the port of CBlastFastaInputSource, `blastinput/fasta_reader`) that reads
// the registered inputs; its records are `api::FastaRecord`. The inputs of a run whose reading
// needs its options (`-parse_deflines`: the ids of the deflines, read in one scope of ids for the
// query and the subjects) are read together, by the reading of the library API's `read_inputs`.
pub use crate::api::input::{parse_deflines_of, read_search_inputs};
pub use crate::blastinput::fasta_reader::{
    read_all, FastaInputSource, ParseErrorCode, ReadError, ReaderConfig, RecordSource,
};
// The reader's own counts for LOSAT Web's index scan (`web/adapter/src/scan/ncbi.rs`), which reads
// the inputs as the reader does without holding their letters: the letters of a `>?` gap line and
// the refusal of a record longer than LOSAT reads.
pub use crate::blastinput::fasta_reader::{gap_line_letters, record_too_long, MAX_RECORD_LETTERS};
pub use crate::blastinput::input_files::is_blank;

/// The options of `program`'s command line that the browser's searches of records do not take:
/// the search strategy file, which the library API refuses too (`api::search::file_options`), and
/// tblastn's `-in_pssm` (the query as a PSSM file), which the browser refuses until line PS's PSW1
/// (`tblastn::refuse_browser_pssm_query`; the library API's records take it, `api::Inputs::pssm`).
/// The adapter's `describe` leaves them out.
pub fn records_refused_options(program: crate::api::Program) -> &'static [&'static str] {
    match program {
        crate::api::Program::Tblastn => &["in_pssm", "export_search_strategy"],
        program => crate::api::file_options(program),
    }
}

// A search strategy file given as bytes (`-import_search_strategy`, `file_begin` kind 3): the
// engine's `apply_strategy` over them, the file that the parsed arguments name, and the warnings
// of the import, which the run writes first (`api::strategy`; the same path as the library's
// `Search::read_strategy`).
pub use crate::api::strategy::{
    apply as import_strategy, import_name as strategy_file_name,
    take_warnings as take_strategy_warnings,
};

// The WindowMasker unit counts of `-window_masker_db` that the adapter holds in memory under their
// names (`file_begin`): NCBI's check of the format and the counts readers on the bytes, and the
// source through which BLASTN opens them (`BlastnArgs::window_masker_counts`).
pub use crate::utils::winmask::{
    discover_stat_type_in, CountsSource, MaskerParams, StatFile, StatType, UnitCounts, WinMaskError,
};

// The decoders of compressed inputs (`blastinput/compressed.rs`): the name rule (`.gz`, `.bz2`,
// `.zst`, in any case), the compression with its sentence for data that is not that format, and
// the push form of the decoder that the command line's reader runs, which inflates a file the
// application gives in pieces (`inflate_begin`/`inflate_chunk`/`inflate_end` of ABI v2).
pub use crate::blastinput::compressed::compression_of_name;
pub use crate::blastinput::compressed::Compression;
pub use crate::blastinput::compressed::PushInflater;
