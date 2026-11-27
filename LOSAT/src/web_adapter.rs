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

// The command line's arguments (`cli.rs`): the option types of every program, how NCBI's reader
// and LOSAT's refusals end an argv (the text that the command line prints where its reader
// stopped at a help word, an argument error's report, a refusal), and the command line's own
// errors.
pub use crate::cli::{ArgvError, Cli, Commands, NativeError, WEB_V2_PROGRAMS};

// The values that NCBI's arguments allow, as LOSAT's argument tables give them (`ncbi_app::args`,
// made from each program's NCBI XML): the genetic codes that `describe` lists and the tasks whose
// defaults it resolves (the defaults and the choices come with `described_parameters`).
pub use crate::ncbi_app::args::{genetic_codes, string_values};

/// `cli::try_parse_with_printing_start` for the command line's own arguments (`Cli`): the argv read
/// as the command line reads it after its application start (NCBI's reader, LOSAT's refusals,
/// the program's options: `cli::reader`). A run of ABI v2 makes its own start (`start_run`), so the
/// toolkit's words that the start reads (`-logfile`, `-conffile`, `-dryrun`) are taken, as the
/// library API takes them, and the help words stop the reader as on the command line, `-xmlhelp`
/// included: the run prints their text (`ArgvError::Printed`; line WQA's D14, 2026-11-20). The
/// versions are printed by the start (`CallerStart::version`) and never reach the reader.
///
/// The posts of an argument error go to the `-logfile` target of the run's start `started`
/// (`CallerStart::target`: the file of `-logfile`, which the run hands over with
/// `CallerStart::finish`), and the error's report is what the command line writes to its standard
/// error then: NCBI's usage, its line that points to the error logs, the line of `=`. Without a
/// start (an argv read on its own: `None`) the posts are the report's, as without `-logfile`.
///
/// Only `Cli`: the adapter parses nothing else, so the façade does not offer the parser for any
/// other type.
pub fn try_parse_cli<I, S>(argv: I, started: Option<&CallerStart>) -> Result<Cli, ArgvError>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let target = started.map_or(
        &crate::blastinput::ncbi_environment::DiagTarget::Stderr,
        CallerStart::target,
    );
    crate::cli::try_parse_with_printing_start(argv, target)
}

/// The message of an argv that the reader or LOSAT's refusals end (`try_parse_cli`), as LOSAT Web
/// gives it to the host (docs/web/abi_v2.md §4): an argument error's report as the command line
/// writes it to its standard error, ending with its line end; LOSAT's refusal as the library API
/// gives it (`api::Error::Argv`: without the command line's `Error: ` and newline); the text that
/// the command line prints where the reader stopped at a help word.
pub fn argv_message(error: &ArgvError) -> String {
    error.to_string()
}

/// Scans `words` (the command line after the program's name) using `program`'s NCBI argument
/// definitions (`ncbi_app::args::scan_argv`), calling `visit(canonical_name, value)` for each
/// recognized option before `--` (`None` for a flag, `Some(value)` for an option taking a value).
pub fn scan_program_argv(
    program: &str,
    words: &[&str],
    mut visit: impl FnMut(&'static str, Option<&str>),
) {
    crate::ncbi_app::args::scan_argv(program, words, |name, value| {
        visit(name, value.and_then(std::ffi::OsStr::to_str));
    });
}

/// The input file names that a `windowmasker` run with `words` opens for reading
/// (`windowmasker::read_files`): `(for_logfile, for_out)`.
pub fn windowmasker_read_files(words: &[&str]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    crate::windowmasker::read_files(words)
}

/// File names that `seqdb` checks or opens for a search's `-db` (`dbs`) and BLASTN's
/// `-filtering_db` (`filtering_dbs`), plus the taxonomy lookup files
/// (`seqdb::taxdb::lookup_files_for_use`) when `outfmt` consults `taxdb`
/// (`algorithm::common::database::outfmt_taxdb_use`) and `taxonomy4blast.sqlite3` when
/// `expands_taxids` is true and a taxonomy file exists on the search path (line WQA D19b).
pub fn search_database_read_files(
    program: crate::api::Program,
    arguments: &[&str],
    registry: &dyn RegistryFiles,
    dbs: &[&str],
    filtering_dbs: &[&str],
    outfmt: &str,
    expands_taxids: bool,
) -> Vec<Vec<u8>> {
    use crate::algorithm::common::database::{
        outfmt_prefetches_taxdb, outfmt_taxdb_use, OutfmtTaxdbUse,
    };
    use crate::report::archive::options::RequestKind;
    let request_kind = match program {
        crate::api::Program::Blastp => RequestKind::Blastp,
        crate::api::Program::Blastn => RequestKind::Blastn,
        crate::api::Program::Blastx => RequestKind::Blastx,
        crate::api::Program::Tblastn => RequestKind::Tblastn,
        crate::api::Program::Tblastx => RequestKind::Tblastx,
        crate::api::Program::Psiblast => RequestKind::Psiblast,
    };
    let has_db = !dbs.is_empty();
    let taxdb_use = outfmt_taxdb_use(outfmt, request_kind, has_db);
    if !has_db && filtering_dbs.is_empty() && taxdb_use == OutfmtTaxdbUse::None {
        return Vec::new();
    }
    let blastdb = crate::blastinput::ncbi_environment::ApplicationStart::with_files(
        program.name(),
        arguments,
        &[],
        registry,
    )
    .settings()
    .ok()
    .and_then(|settings| settings.blastdb);
    let seqtype = match program {
        crate::api::Program::Blastp
        | crate::api::Program::Blastx
        | crate::api::Program::Psiblast => crate::seqdb::SeqType::Protein,
        crate::api::Program::Blastn
        | crate::api::Program::Tblastn
        | crate::api::Program::Tblastx => crate::seqdb::SeqType::Nucleotide,
    };
    let mut out = Vec::new();
    let push_db = |out: &mut Vec<Vec<u8>>, db: &str, db_type: crate::seqdb::SeqType| {
        let trimmed = crate::seqdb::db_argument_name_bytes(db.as_bytes());
        out.extend(crate::seqdb::alias::database_lookup_files(
            trimmed,
            db_type,
            blastdb.as_deref(),
        ));
    };
    for db in dbs {
        push_db(&mut out, db, seqtype);
    }
    let prefetches = has_db && outfmt_prefetches_taxdb(outfmt);
    out.extend(crate::seqdb::taxdb::lookup_files_for_use(
        blastdb.as_deref(),
        taxdb_use,
        prefetches,
    ));
    if has_db && expands_taxids {
        let tax_file = crate::blastinput::db_filter::TAXONOMY_FILE.as_bytes();
        let search = crate::seqdb::SearchPath::from_environment(blastdb.as_deref());
        if search.resolve_file(tax_file).is_some() {
            out.push(tax_file.to_vec());
        }
    }
    for fdb in filtering_dbs {
        push_db(&mut out, fdb, crate::seqdb::SeqType::Nucleotide);
    }
    out
}

// The NCBI application start of a run (`api::start`, line RG item E1w-3a): the words that it reads
// before the arguments (`preparse`, NCBI's scan of `-logfile`), the registry files that the host
// registered (`RegistryFiles`), and where the run's files go (`RunFiles`: the file of `-logfile`).
pub use crate::api::start::{CallerSink, CallerStart, ErrorSource};
pub use crate::blastinput::ncbi_environment::{
    log_file, preparse, ApplicationSettings, Preparsed, RegistryFiles, TextPlace,
};
pub use crate::report::run_files::{FileReceiver, FileStream, RunFiles};
/// The names of the command line's programs (`losat <program>`, and the names that run one without
/// its word): LOSAT Web's tool run tells a program that it does not run from a name that is none.
pub use crate::utils::entry_names::ENTRY_NAMES as PROGRAMS;

// The search strategy file of a run (`-export_search_strategy`, line RG item F1w): made through the
// run's files where the command line makes it, and given to each search of the run.
pub use crate::report::archive::strategy::StrategyFile;

/// What LOSAT Web reads of the arguments of a search that may write `-export_search_strategy`
/// (`search_files`): the name of the strategy file, the program's label in NCBI's messages, the
/// files of `-query`, `-subject` and BLASTN's `-window_masker_db`, and what the search reads from
/// the run's working directory: the filters of the database of `-db` (`filter`, `None` without
/// `-db`) and BLASTN's `-filtering_db` (`repeats`).
struct SearchFiles<'a> {
    strategy: Option<&'a std::path::Path>,
    label: &'a str,
    query: &'a std::path::Path,
    subject: Option<&'a std::path::Path>,
    in_pssm: Option<&'a std::path::Path>,
    in_msa: Option<&'a std::path::Path>,
    out_pssm: Option<&'a std::path::Path>,
    out_ascii_pssm: Option<&'a std::path::Path>,
    masker: Option<&'a std::ffi::OsStr>,
    filter: Option<&'a crate::blastinput::db_filter::DbFilterArgs>,
    repeats: Option<&'a std::ffi::OsStr>,
}

impl<'a> SearchFiles<'a> {
    fn new(
        label: &'static str,
        strategy: Option<&'a std::path::Path>,
        query: &'a std::path::Path,
        subject: Option<&'a std::path::Path>,
        filter: Option<&'a crate::blastinput::db_filter::DbFilterArgs>,
    ) -> Self {
        Self {
            strategy,
            label,
            query,
            subject,
            in_pssm: None,
            in_msa: None,
            out_pssm: None,
            out_ascii_pssm: None,
            masker: None,
            filter,
            repeats: None,
        }
    }
}

/// `SearchFiles` of the arguments `command`: `None` for a program that writes no strategy file.
fn search_files(command: &Commands) -> Option<SearchFiles<'_>> {
    macro_rules! base_files {
        ($label:expr, $args:expr) => {
            SearchFiles::new(
                $label,
                $args.export_search_strategy.as_deref(),
                $args.query.as_path(),
                $args.subject.as_deref(),
                $args.db.is_some().then_some(&$args.db_filter),
            )
        };
    }
    match command {
        Commands::Blastp(args) => {
            let mut files = base_files!(args.app.name(), args);
            if let Some(pssm) = args.psiblast_pssm.as_ref() {
                files.in_pssm = pssm.in_pssm.as_deref();
                files.in_msa = pssm.in_msa.as_ref().map(|m| m.path.as_path());
                files.out_pssm = pssm.out_pssm.as_deref();
                files.out_ascii_pssm = pssm.out_ascii_pssm.as_deref();
            }
            Some(files)
        }
        Commands::Psiblast(args) => {
            let mut files = base_files!("PSIBLAST", args);
            files.in_pssm = args.in_pssm.as_deref();
            files.in_msa = args.in_msa.as_deref();
            files.out_pssm = args.out_pssm.as_deref();
            files.out_ascii_pssm = args.out_ascii_pssm.as_deref();
            Some(files)
        }
        Commands::Blastn(args) => {
            let mut files = base_files!("BLASTN", args);
            files.masker = args.window_masker_db.as_deref();
            files.repeats = args.filtering_db.as_deref();
            Some(files)
        }
        Commands::Tblastn(args) => Some(SearchFiles {
            in_pssm: args.in_pssm.as_deref(),
            ..base_files!("TBLASTN", args)
        }),
        Commands::Tblastx(args) => Some(base_files!("TBLASTX", args)),
        Commands::Blastx(args) => Some(base_files!("BLASTX", args)),
        _ => None,
    }
}

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
    let Some(files) = search_files(command) else {
        return Ok(None);
    };
    let mut others = vec![RunFileName {
        option: "query",
        name: files.query.as_os_str(),
    }];
    others.extend(files.subject.map(|subject| RunFileName {
        option: "subject",
        name: subject.as_os_str(),
    }));
    others.extend(files.in_pssm.map(|in_pssm| RunFileName {
        option: "in_pssm",
        name: in_pssm.as_os_str(),
    }));
    others.extend(files.in_msa.map(|in_msa| RunFileName {
        option: "in_msa",
        name: in_msa.as_os_str(),
    }));
    others.extend(files.out_pssm.map(|out_pssm| RunFileName {
        option: "out_pssm",
        name: out_pssm.as_os_str(),
    }));
    others.extend(files.out_ascii_pssm.map(|out_ascii_pssm| RunFileName {
        option: "out_ascii_pssm",
        name: out_ascii_pssm.as_os_str(),
    }));
    others.extend(log_file(words).map(|name| RunFileName {
        option: "logfile",
        name,
    }));
    others.extend(files.masker.map(|masker| RunFileName {
        option: "window_masker_db",
        name: masker,
    }));
    // The files of the run's working directory that the search reads (docs/web/abi_v2.md §3), which
    // the host gives the engine there: with `-db` (line WQA D12a) the database and the lists of its
    // filters, with BLASTN's `-filtering_db` (line RG item D4) the database of repeats. The lists are
    // compared under their options; any other file of that directory that a relative name opens is
    // a database's (the taxonomy database's too): `-filtering_db`'s for a name that starts with its
    // value or in a run without `-db`, else `-db`'s.
    if files.filter.is_some() || files.repeats.is_some() {
        if let Some(filter) = files.filter {
            for (option, list) in [
                ("gilist", &filter.gilist),
                ("negative_gilist", &filter.negative_gilist),
                ("seqidlist", &filter.seqidlist),
                ("negative_seqidlist", &filter.negative_seqidlist),
                ("taxidlist", &filter.taxidlist),
                ("negative_taxidlist", &filter.negative_taxidlist),
                ("ipglist", &filter.ipglist),
                ("negative_ipglist", &filter.negative_ipglist),
            ] {
                others.extend(list.as_deref().map(|name| RunFileName {
                    option,
                    name: std::ffi::OsStr::new(name),
                }));
            }
        }
        let read_there = files
            .strategy
            .filter(|path| path.is_relative() && path.is_file());
        if let Some(path) = read_there {
            let repeats = files.filter.is_none()
                || files.repeats.is_some_and(|repeats| {
                    crate::utils::os_bytes::os_str_bytes(path.as_os_str())
                        .starts_with(&crate::utils::os_bytes::os_str_bytes(repeats))
                });
            others.push(RunFileName {
                option: if repeats { "filtering_db" } else { "db" },
                name: path.as_os_str(),
            });
        }
    }
    crate::report::archive::strategy::open_strategy_file_in(
        run_files,
        files.strategy,
        files.label,
        &others,
    )
}

/// Gives the search of `command` (a search of LOSAT Web's run of the argv, `run_local`) the export
/// of its requests into `file` (`open_strategy_file`), made as the command line's run makes it
/// (`StrategyExport::for_args`, the program's `StrategyArgs`), with the `HitlistSize` of the argv's
/// format numbered `outfmt`, none for a value that the early reading left to the formatting
/// handler (`early_format_choice`), as the command line's (the archive's request shares the
/// export's options handle, `run_local` and `run_local_database`).
pub fn give_strategy(command: &mut Commands, file: StrategyFile, outfmt: Option<i32>) {
    use crate::report::archive::strategy::{StrategyArgs, StrategyExport};
    let options = crate::report::archive::OptionsHandle::default();
    let export =
        |args: &dyn StrategyArgs| Some(StrategyExport::for_args(args, file, outfmt, &options));
    match command {
        Commands::Blastp(args) => args.strategy = export(&*args),
        Commands::Blastn(args) => args.strategy = export(&*args),
        Commands::Tblastn(args) => args.strategy = export(&*args),
        Commands::Tblastx(args) => args.strategy = export(&*args),
        Commands::Blastx(args) => args.strategy = export(&*args),
        _ => {}
    }
}

/// The `-export_search_strategy` name of the arguments `command`, if any: a run that writes the file
/// reads its records as the command line's reader then reads them (their deflines kept, LOSAT
/// Web's `store`).
pub fn export_search_strategy_of(command: &Commands) -> Option<&std::path::Path> {
    search_files(command)?.strategy
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

/// Runs `windowmasker` with `words` (the command line after `losat windowmasker`) as the command
/// line runs it in its working directory, and its exit status (the tool runs of ABI v2,
/// `losat_web2_tool`; line NP WW3a, 2026-11-22): it opens its input files by name in the run's
/// working directory, which LOSAT Web's host gives the engine (docs/web/abi_v2.md §3), reads `input`
/// where the command line reads its standard input, writes its standard output to `out` and its
/// standard error to `err`, and gives the files that it writes (`-out FILE`, `-logfile FILE`) to
/// `run_files` (`RunFiles::Receiver`: each file as it is written when the receiver opens it,
/// `FileReceiver::open`, else whole when the run ends; one file per path, also where the run opens
/// it twice, `-out x -logfile x`, as the system shares it). The application start is the command
/// line's (its environment and its registry files: the run's `.ncbirc`, the file of `-conffile`),
/// and `-xmlhelp` prints NCBI's XML, as there.
pub fn run_windowmasker(
    words: &[&str],
    run_files: &RunFiles,
    input: &mut dyn std::io::Read,
    out: &mut dyn std::io::Write,
    err: &mut dyn std::io::Write,
) -> i32 {
    let words: Vec<OsString> = words.iter().map(OsString::from).collect();
    crate::windowmasker::run_in(&words, run_files, input, out, err)
}

/// Checks the words of a `windowmasker` run as `run_windowmasker` runs them, up to where the tool
/// reads its input or opens a file: the application start (which reads the registry files of the
/// run's working directory, `.ncbirc`, the file of `-conffile`, as the command line does), the
/// arguments (NCBI's reader and its checks), LOSAT's settings and the checks of the mode before its
/// files (`-mk_counts` with `-ids` and `-exclude_ids`, or with `-infmt seqids`; `-convert -in -`)
/// (LOSAT Web's validation of a tool's words, `losat_web2_validate`; line NP WW3a, D-WW3-7 and
/// D-WW3-10, 2026-11-23). `Err` with what the command line writes for the words that it refuses:
/// its standard error, then the file of `-logfile` where the argv has one, without the last line
/// end.
pub fn check_windowmasker(words: &[&str]) -> Result<(), String> {
    use std::sync::{Arc, Mutex};
    /// The file of `-logfile`, kept for the message, carrying the tool's read input names so
    /// `-logfile` and `-out` naming a missing input that the CLI reads back are refused.
    struct Log {
        bytes: Mutex<Vec<u8>>,
        read: Vec<Vec<u8>>,
        out_read: Vec<Vec<u8>>,
    }
    impl FileReceiver for Log {
        fn receive(&self, _name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
            self.bytes
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .extend_from_slice(bytes);
            Ok(())
        }

        fn read_files(&self) -> &[Vec<u8>] {
            &self.read
        }

        fn output_read_files(&self) -> &[Vec<u8>] {
            &self.out_read
        }
    }
    let words: Vec<OsString> = words.iter().map(OsString::from).collect();
    let (read, out_read) = crate::windowmasker::read_files(&words);
    let log = Arc::new(Log {
        bytes: Mutex::new(Vec::new()),
        read,
        out_read,
    });
    let mut out = Vec::new();
    let mut err = Vec::new();
    let status =
        crate::windowmasker::check_in(&words, &RunFiles::Receiver(log.clone()), &mut out, &mut err);
    if status == 0 {
        return Ok(());
    }
    err.extend_from_slice(&out);
    err.extend_from_slice(
        &log.bytes
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner),
    );
    if err.last() == Some(&b'\n') {
        err.pop();
    }
    Err(String::from_utf8_lossy(&err).into_owned())
}

/// Runs `makeblastdb` with `words` (the command line after `losat makeblastdb`) as the command
/// line runs it in its working directory, and its exit status (the tool runs of ABI v2,
/// `losat_web2_tool`; line NP WW4b1, D-WW4b-1, 2026-11-27): it reads and writes its files by name
/// in the run's working directory as the command line does (the file system, which LOSAT Web's
/// host gives the engine for a tool run: docs/web/abi_v2.md §3), reads `input` where the command
/// line reads its standard input, writes its standard output to `out` and its standard error to
/// `err`, and its dates are those of the zone of `TZ` (`writedb::time::LocalTime::now`). The
/// application start is the command line's (its environment and its registry files: the run's
/// `.ncbirc`, `makeblastdb.ini`, the file of `-conffile`), and `-xmlhelp` prints NCBI's XML, as
/// there.
pub fn run_makeblastdb(
    words: &[&str],
    input: &mut dyn std::io::Read,
    out: &mut dyn std::io::Write,
    err: &mut dyn std::io::Write,
) -> i32 {
    crate::api::tools::makeblastdb_call(
        words,
        input,
        out,
        err,
        crate::api::tools::MakeblastdbCall::WebRun,
    )
    .status
}

/// Checks the words of a `makeblastdb` run as `run_makeblastdb` runs them, up to where the build
/// begins (`api::run_makeblastdb` before `Run`): the application start (which reads the registry
/// files of the run's working directory), NCBI's reader, LOSAT's refusals and settings, and
/// makeblastdb's checks of the environment and the registry (LOSAT Web's validation of a tool's
/// words, `losat_web2_validate`; line NP WW4b1, D-WW4b-3, 2026-11-27, as D-WW3-7). It opens no
/// file for writing: the file of `-logfile` is kept in memory. `Err` with what the command line
/// writes for the words that it refuses: its standard error, then its standard output and the file
/// of `-logfile`, without the last line end.
pub fn check_makeblastdb(words: &[&str]) -> Result<(), String> {
    let mut out = Vec::new();
    let mut err = Vec::new();
    let end = crate::api::tools::makeblastdb_call(
        words,
        &mut std::io::empty(),
        &mut out,
        &mut err,
        crate::api::tools::MakeblastdbCall::WebCheck,
    );
    if end.status == 0 {
        return Ok(());
    }
    err.extend_from_slice(&out);
    err.extend_from_slice(&end.log);
    if err.last() == Some(&b'\n') {
        err.pop();
    }
    Err(String::from_utf8_lossy(&err).into_owned())
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
pub use crate::algorithm::blastp::args::ListCut;
pub use crate::algorithm::blastp::blast_engine::check_options as check_blastp_options;
pub use crate::algorithm::blastx::check_options as check_blastx_options;
pub use crate::algorithm::blastx::BlastxArgs;
pub use crate::algorithm::tblastn::check_options as check_tblastn_options;
pub use crate::algorithm::tblastx::blast_engine::check_options as check_tblastx_options;
pub use crate::utils::threading::validate_threads;

/// Routes a `psiblast` run's `-out_pssm` and `-out_ascii_pssm` writes through `run_files`
/// (`RunFiles::Receiver` in LOSAT Web). Nothing if `args` is not a `psiblast` search.
pub fn give_pssm_run_files(args: &mut crate::algorithm::blastp::BlastpArgs, run_files: RunFiles) {
    if let Some(pssm) = args.psiblast_pssm.as_mut() {
        pssm.run_files = run_files;
    }
}

/// Checks `outfmt` against the output formats that LOSAT's `psiblast` supports
/// (`PsiblastArgs::check_outfmt`).
pub fn check_psiblast_outfmt(outfmt: &str, strategy: bool) -> anyhow::Result<()> {
    crate::algorithm::psiblast::PsiblastArgs::check_outfmt(outfmt, None, strategy, false)
}

/// Refuses `-phi_pattern` in LOSAT Web's `psiblast` (PHI-BLAST pattern files are not supported in
/// the browser host yet; `records_refused_options`).
pub fn refuse_psiblast_phi_pattern(
    args: &crate::algorithm::blastp::BlastpArgs,
) -> anyhow::Result<()> {
    if args
        .psiblast_pssm
        .as_ref()
        .is_some_and(|pssm| pssm.phi_pattern.is_some())
    {
        anyhow::bail!("option -phi_pattern is not supported by LOSAT Web's PSIBLAST yet");
    }
    Ok(())
}

/// How `psiblast`'s formatter of `outfmt` cuts hit lists (`BlastpApp::Psiblast.list_cut`), used
/// by `web/adapter` to group formats into searches.
pub fn psiblast_list_cut(outfmt: i32) -> ListCut {
    crate::algorithm::blastp::args::BlastpApp::Psiblast.list_cut(outfmt)
}

// The searches: each program's `run_local` over records in memory and its `run_local_database`
// over the database of `-db` (opened by name with the run's application settings, as the command
// line opens it: LOSAT Web's host gives the engine the run's working directory, line WQA D12a), the
// pool that the searches of one argv share, and the routing of their outputs (`api/local_blast.rs`).
pub use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_database_blastn,
    run_local_database_blastp, run_local_database_blastx, run_local_database_tblastn,
    run_local_database_tblastx, run_local_tblastn, run_local_tblastx, with_run_pool,
    FormatObserver, FormatOutput, FormatRefusal, HspIndex, OutputSink, ReportOutputs,
    ReportedSubject,
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
    early_format_choice, engine_error, formatting_handler_check, missing_subject_error,
    parse_formatting_string, report_format, EarlyFormat,
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
// the inputs as the reader does without holding their letters: the letters of a `>?` gap line,
// the refusal of a record longer than LOSAT reads, and the `-parse_deflines` molecule and letter
// rules (`defline_named_molecule`, `DataLineAlphabet`).
pub use crate::blastinput::fasta_reader::believed::defline_named_molecule;
pub use crate::blastinput::fasta_reader::reader::DataLineAlphabet;
pub use crate::blastinput::fasta_reader::{gap_line_letters, record_too_long, MAX_RECORD_LETTERS};

/// Every option of `program` that LOSAT's descriptions list, in the order of the usage line of the
/// command line's `-help`, as the library API's `describe` gives each (NCBI's argument tables,
/// line RG DECISIONS K4i-3-1): the rows of the adapter's `describe`, which leaves out the options
/// that the adapter sets itself. `api::describe` lists the same rows but those that the library
/// API sets.
pub fn described_parameters(program: crate::api::Program) -> Vec<crate::api::Parameter> {
    crate::api::described_parameters(program).collect()
}

/// The options of `program`'s command line that the browser's searches of records do not take.
/// Since line PS's PSW1, tblastn's `-in_pssm` and psiblast's `-in_pssm`/`-in_msa`/`-out_pssm`/
/// `-out_ascii_pssm` are supported in the browser; psiblast's `-phi_pattern` (item PS8a) is
/// refused in LOSAT Web until pattern files are wired.
pub fn records_refused_options(program: crate::api::Program) -> &'static [&'static str] {
    match program {
        crate::api::Program::Psiblast => &["phi_pattern"],
        _ => &[],
    }
}

// A search strategy file given as bytes (`-import_search_strategy`, `file_begin` kind 3): the
// engine's import over them (the strategy's program's arguments, `program_of`), the file that the
// parsed arguments name, and the warnings of the import, which the run writes first
// (`api::strategy`; the same path as the command line's and the library's
// `Search::read_strategy`).
pub use crate::api::strategy::{
    apply as import_strategy, import_name as strategy_file_name, program_of,
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
