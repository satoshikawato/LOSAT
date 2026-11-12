//! The engine items that LOSAT Web's adapter (`web/adapter`, ABI v2: `docs/web/abi_v2.md`)
//! names, and nothing else. No NCBI counterpart: this is the boundary between this crate and
//! the adapter crate.
//!
//! Hidden from the documentation and outside the compatibility promise of the version number,
//! as the hidden modules are; `api` is the public surface. An item here can move to `api` later,
//! one at a time.
//!
//! Rules:
//! - The adapter names the engine through `LOSAT::api` (public items) and this module only.
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
// messages, and the command line's own errors.
pub use crate::cli::{render_message, Cli, Commands, NativeError};

/// `cli::try_parse_from` for the command line's own arguments (`Cli`): the argv read as the
/// command line reads it, NCBI's single-dash options included.
///
/// Only `Cli`: the adapter parses nothing else, so the façade does not offer the parser for any
/// other clap type.
pub fn try_parse_cli<I, S>(argv: I) -> Result<Cli, clap::Error>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    crate::cli::try_parse_from(argv)
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
pub use crate::algorithm::tblastx::blast_engine::check_options as check_tblastx_options;
pub use crate::utils::threading::validate_threads;

// The searches: each program's `run_local` over records in memory, the pool that the searches
// of one argv share, and the routing of their outputs (`api/local_blast.rs`).
pub use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
    with_run_pool, FormatObserver, FormatOutput, HspIndex, OutputSink, ReportOutputs,
};
pub use crate::report::blast_xml2::{SideFileReceiver, SideFileSink};
pub use crate::report::PairwiseHit;

// The output formats of the application layer (`blastinput/app.rs`, NCBI's CFormattingArgs):
// `-outfmt` parsed and checked as the command line does. `EngineProgram` is the program of
// these checks (the adapter has its own `Program`, and `api::Program` is another type).
pub use crate::blastinput::app::output_format::{ARCHIVE_FORMAT, SAM as SAM_FORMAT};
pub use crate::blastinput::app::Program as EngineProgram;
pub use crate::blastinput::app::{
    engine_error, formatting_handler_check, parse_formatting_string, report_format,
};

// The hit list sizes of the formats (`blastinput/report_counts.rs`), the scoring matrices that
// `describe` lists, and the `-query_loc`/`-subject_loc` ranges (`blastinput/seq_range.rs`).
pub use crate::blastinput::report_counts::{report_counts, run_report_counts};
pub use crate::blastinput::seq_range::{
    parse_optional_range, record_interval, RangeRole, RecordInterval,
};
pub use crate::config::ScoringMatrix;

// The FASTA reader (the port of CBlastFastaInputSource, `blastinput/fasta_reader`) that reads
// the registered inputs; its records are `api::FastaRecord`.
pub use crate::blastinput::fasta_reader::{
    read_all, FastaInputSource, ParseErrorCode, ReadError, ReaderConfig, RecordSource,
};
pub use crate::blastinput::input_files::is_blank;

// The decoders of compressed inputs (`blastinput/compressed.rs`): the name rule (`.gz`, `.bz2`,
// `.zst`, in any case), the compression with its sentence for data that is not that format, and
// the push form of the decoder that the command line's reader runs, which inflates a file the
// application gives in pieces (`inflate_begin`/`inflate_chunk`/`inflate_end` of ABI v2).
pub use crate::blastinput::compressed::compression_of_name;
pub use crate::blastinput::compressed::Compression;
pub use crate::blastinput::compressed::PushInflater;
