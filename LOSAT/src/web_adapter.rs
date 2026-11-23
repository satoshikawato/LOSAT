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

// The values that NCBI's arguments allow, as LOSAT's argument tables give them (`ncbi_app::args`,
// made from each program's NCBI XML): the genetic codes that `describe` lists and the tasks whose
// defaults it resolves.
pub use crate::ncbi_app::args::{genetic_codes, string_values};

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

// The search strategy file of a run (`-export_search_strategy`, line RG item F1w): made through the
// run's files where the command line makes it, and given to each search of the run.
pub use crate::report::archive::strategy::StrategyFile;

/// The `-export_search_strategy` file of the arguments `command` (`words`: the arguments after the
/// program's name) of a run whose files go to `run_files` (LOSAT Web: `RunFiles::Receiver`), made
/// where the command line makes it (`strategy::open_strategy_file_in`, as the command line's
/// `open_strategy_file`): after the settings and a dry run's end, before `-outfmt` is read. `None`
/// without the option. The error is the command line's (NCBI's `File is not accessible` for a name
/// that the command line cannot open); the run's other files that the name may be are its `-query`,
/// `-subject`, `-logfile`, BLASTN's `-window_masker_db`, and the files of the run's working directory
/// that the search reads (`-db`, the lists of its filters, BLASTN's `-filtering_db`).
pub fn open_strategy_file(
    command: &Commands,
    words: &[&str],
    run_files: &RunFiles,
) -> anyhow::Result<Option<StrategyFile>> {
    use crate::report::run_files::RunFileName;
    use std::path::Path;
    let (path, label, query, subject, masker): (Option<&Path>, &str, &Path, _, Option<&str>) =
        match command {
            Commands::Blastp(args) => (
                args.export_search_strategy.as_deref(),
                args.app.name(),
                args.query.as_path(),
                args.subject.as_deref(),
                None,
            ),
            Commands::Blastn(args) => (
                args.export_search_strategy.as_deref(),
                "BLASTN",
                args.query.as_path(),
                args.subject.as_deref(),
                args.window_masker_db.as_deref(),
            ),
            Commands::Tblastn(args) => (
                args.export_search_strategy.as_deref(),
                "TBLASTN",
                args.query.as_path(),
                args.subject.as_deref(),
                None,
            ),
            Commands::Tblastx(args) => (
                args.export_search_strategy.as_deref(),
                "TBLASTX",
                args.query.as_path(),
                args.subject.as_deref(),
                None,
            ),
            Commands::Blastx(args) => (
                args.export_search_strategy.as_deref(),
                "BLASTX",
                args.query.as_path(),
                args.subject.as_deref(),
                None,
            ),
            _ => return Ok(None),
        };
    let mut others = vec![RunFileName {
        option: "query",
        name: query.as_os_str(),
    }];
    others.extend(subject.map(|subject| RunFileName {
        option: "subject",
        name: subject.as_os_str(),
    }));
    others.extend(log_file(words).map(|name| RunFileName {
        option: "logfile",
        name,
    }));
    others.extend(masker.map(|masker| RunFileName {
        option: "window_masker_db",
        name: std::ffi::OsStr::new(masker),
    }));
    // The files of the run's working directory that the search reads (docs/web/abi_v2.md §3), which
    // the host gives the engine there: with `-db` (line WQA D12a) the database and the lists of its
    // filters, with BLASTN's `-filtering_db` (line RG item D4) the database of repeats. The lists are
    // compared under their options; any other file of that directory that a relative name opens is
    // a database's (the taxonomy database's too): `-filtering_db`'s for a name that starts with its
    // value or in a run without `-db`, else `-db`'s.
    if let Some(reads) = directory_reads(command) {
        if let Some(filter) = reads.filter {
            for (option, list) in [
                ("seqidlist", &filter.seqidlist),
                ("negative_seqidlist", &filter.negative_seqidlist),
                ("taxidlist", &filter.taxidlist),
                ("negative_taxidlist", &filter.negative_taxidlist),
            ] {
                others.extend(list.as_deref().map(|name| RunFileName {
                    option,
                    name: std::ffi::OsStr::new(name),
                }));
            }
        }
        if let Some(path) = path.filter(|path| path.is_relative() && path.is_file()) {
            let repeats = reads.filter.is_none()
                || reads.repeats.is_some_and(|repeats| {
                    path.as_os_str()
                        .as_encoded_bytes()
                        .starts_with(repeats.as_bytes())
                });
            others.push(RunFileName {
                option: if repeats { "filtering_db" } else { "db" },
                name: path.as_os_str(),
            });
        }
    }
    crate::report::archive::strategy::open_strategy_file_in(run_files, path, None, label, &others)
}

/// What a search reads from the run's working directory (`open_strategy_file`): the filters of the
/// database of `-db` (`filter`, `None` without `-db`) and BLASTN's `-filtering_db` (`repeats`).
struct DirectoryReads<'a> {
    filter: Option<&'a crate::blastinput::db_filter::DbFilterArgs>,
    repeats: Option<&'a str>,
}

/// What the search of `command` reads from the run's working directory; `None` for a search that
/// reads no file there.
fn directory_reads(command: &Commands) -> Option<DirectoryReads<'_>> {
    let (db, filter, repeats) = match command {
        Commands::Blastp(args) => (&args.db, &args.db_filter, None),
        Commands::Blastn(args) => (&args.db, &args.db_filter, args.filtering_db.as_deref()),
        Commands::Tblastn(args) => (&args.db, &args.db_filter, None),
        Commands::Tblastx(args) => (&args.db, &args.db_filter, None),
        Commands::Blastx(args) => (&args.db, &args.db_filter, None),
        _ => return None,
    };
    let filter = db.is_some().then_some(filter);
    (filter.is_some() || repeats.is_some()).then_some(DirectoryReads { filter, repeats })
}

/// Gives the search of `command` (a search of LOSAT Web's run of the argv, `run_local`) the export
/// of its requests into `file` (`open_strategy_file`), as the command line's run makes it
/// (`StrategyExport::for_run`): the recorded options of the arguments, with the `HitlistSize` of
/// the argv's format numbered `outfmt` and the database of `-db` as the command line names it
/// (the archive's request shares the export's options handle, `run_local` and
/// `run_local_database`).
pub fn give_strategy(command: &mut Commands, file: StrategyFile, outfmt: i32) {
    use crate::report::archive::strategy::{HitlistArgs, StrategyExport};
    let options = crate::report::archive::OptionsHandle::default();
    let export = |inputs, counts, database: Option<&str>| {
        Some(StrategyExport::for_run(
            file,
            inputs,
            Some(outfmt),
            counts,
            database,
            &options,
        ))
    };
    match command {
        Commands::Blastp(args) => {
            let counts = HitlistArgs {
                max_target_seqs: args.max_target_seqs,
                num_descriptions: args.num_descriptions,
                num_alignments: args.num_alignments,
            };
            args.strategy = export(
                crate::algorithm::blastp::archive::option_inputs(args),
                counts,
                args.db.as_deref().map(crate::seqdb::db_argument_name),
            );
        }
        Commands::Blastn(args) => {
            let counts = HitlistArgs {
                max_target_seqs: args.max_target_seqs,
                num_descriptions: args.num_descriptions,
                num_alignments: args.num_alignments,
            };
            let inputs = crate::algorithm::blastn::archive::option_inputs(
                args,
                args.dust_filtering.as_deref(),
            );
            args.strategy = export(inputs, counts, args.database_name());
        }
        Commands::Tblastn(args) => {
            let counts = HitlistArgs {
                max_target_seqs: args.max_target_seqs,
                num_descriptions: args.num_descriptions,
                num_alignments: args.num_alignments,
            };
            args.strategy = export(
                crate::algorithm::tblastn::archive::option_inputs(args),
                counts,
                args.database_name(),
            );
        }
        Commands::Tblastx(args) => {
            let counts = HitlistArgs {
                max_target_seqs: args.max_target_seqs,
                num_descriptions: args.num_descriptions,
                num_alignments: args.num_alignments,
            };
            args.strategy = export(
                crate::algorithm::tblastx::option_inputs(args),
                counts,
                args.database_name(),
            );
        }
        Commands::Blastx(args) => {
            let counts = HitlistArgs {
                max_target_seqs: args.max_target_seqs.map(|count| count as usize),
                num_descriptions: args.num_descriptions,
                num_alignments: args.num_alignments,
            };
            args.strategy = export(
                crate::algorithm::blastx::option_inputs(args),
                counts,
                args.database_name(),
            );
        }
        _ => {}
    }
}

/// The `-export_search_strategy` name of the arguments `command`, if any: a run that writes the file
/// reads its records as the command line's reader then reads them (their deflines kept, LOSAT
/// Web's `store`).
pub fn export_search_strategy_of(command: &Commands) -> Option<&std::path::Path> {
    match command {
        Commands::Blastp(args) => args.export_search_strategy.as_deref(),
        Commands::Blastn(args) => args.export_search_strategy.as_deref(),
        Commands::Tblastn(args) => args.export_search_strategy.as_deref(),
        Commands::Tblastx(args) => args.export_search_strategy.as_deref(),
        Commands::Blastx(args) => args.export_search_strategy.as_deref(),
        _ => None,
    }
}

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

/// Gives BLASTN's `-filtering_db` its repeat search with the run's `settings` (the lookup of the
/// database as of `-db`: the run's working directory and registry), as the command line does
/// before a search; nothing without the option (`BlastnArgs::resolve_repeat_filter`).
pub fn resolve_blastn_repeat_filter(args: &mut BlastnArgs, settings: &ApplicationSettings) {
    args.resolve_repeat_filter(settings);
}
pub use crate::algorithm::blastp::blast_engine::check_options as check_blastp_options;
pub use crate::algorithm::blastx::check_options as check_blastx_options;
pub use crate::algorithm::blastx::BlastxArgs;
pub use crate::algorithm::tblastn::check_options as check_tblastn_options;
pub use crate::algorithm::tblastn::refuse_browser_pssm_query as refuse_tblastn_pssm_query;
pub use crate::algorithm::tblastx::blast_engine::check_options as check_tblastx_options;
pub use crate::utils::threading::validate_threads;

// The searches: each program's `run_local` over records in memory and its `run_local_database`
// over the database of `-db` (opened by name with the run's application settings, as the command
// line opens it: LOSAT Web's host gives the engine the run's working directory, line WQA D12a), the
// pool that the searches of one argv share, and the routing of their outputs (`api/local_blast.rs`).
pub use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_database_blastn,
    run_local_database_blastp, run_local_database_blastx, run_local_database_tblastn,
    run_local_database_tblastx, run_local_tblastn, run_local_tblastx, with_run_pool,
    FormatObserver, FormatOutput, FormatRefusal, HspIndex, OutputSink, ReportOutputs,
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
// A query is read as the command line's query batches read it (`read_query_records`): the records
// before the reader's error, and the error, which the run meets in its batch.
pub use crate::api::input::{parse_deflines_of, read_query_records, read_search_inputs};
pub use crate::blastinput::fasta_reader::{
    read_all, FastaInputSource, ParseErrorCode, ReadError, ReaderConfig, RecordSource,
};
// The reader's own counts for LOSAT Web's index scan (`web/adapter/src/scan/ncbi.rs`), which reads
// the inputs as the reader does without holding their letters: the letters of a `>?` gap line and
// the refusal of a record longer than LOSAT reads.
pub use crate::blastinput::fasta_reader::{gap_line_letters, record_too_long, MAX_RECORD_LETTERS};

/// The options of `program`'s command line that the browser's searches of records do not take:
/// tblastn's `-in_pssm` (the query as a PSSM file), which the browser refuses until line PS's PSW1
/// (`tblastn::refuse_browser_pssm_query`; the library API's records take it, `api::Inputs::pssm`).
/// The adapter's `describe` leaves them out. The search strategy file, which the library API
/// refuses (`api::search::file_options`), is a file of the run in the browser
/// (`open_strategy_file`, line RG item F1w).
pub fn records_refused_options(program: crate::api::Program) -> &'static [&'static str] {
    match program {
        crate::api::Program::Tblastn => &["in_pssm"],
        _ => &[],
    }
}

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
