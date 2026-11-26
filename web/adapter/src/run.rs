//! The search: argv parsing (`validate`) and `run` (docs/web/abi_v2.md §4, §7-8).

use std::fmt::Write as _;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use losat::api::{FastaRecord, QueryEnd};
use losat::web_adapter::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_database_blastn,
    run_local_database_blastp, run_local_database_blastx, run_local_database_tblastn,
    run_local_database_tblastx, run_local_tblastn, run_local_tblastx, ApplicationSettings,
    CallerSink, CallerStart, Cli, Commands, CountsSource, ErrorSource, FormatObserver,
    FormatOutput, FormatRefusal, HspIndex, OutputSink, PairwiseHit, ReportOutputs, RunFiles,
    StrategyFile,
};

use crate::emit::{
    emit, StreamWriter, STREAM_CUSTOM, STREAM_DIAGNOSTICS, STREAM_HITS, STREAM_PRINTED,
    STREAM_RESPONSE, STREAM_RUN_FILES,
};
use crate::imported::Supplied;
use crate::json;

/// The programs that ABI v2 runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
    Blastp,
    Tblastn,
    Blastn,
    Tblastx,
    Blastx,
}

impl Program {
    pub fn parse(name: &str) -> Result<Self, String> {
        match name {
            "blastp" => Ok(Self::Blastp),
            "tblastn" => Ok(Self::Tblastn),
            "blastn" => Ok(Self::Blastn),
            "tblastx" => Ok(Self::Tblastx),
            "blastx" => Ok(Self::Blastx),
            "psiblast" => Err("psiblast is not available in the browser (ABI v2)".into()),
            "rpsblast" => Err("rpsblast is not available in the browser (ABI v2)".into()),
            "rpstblastn" => Err("rpstblastn is not available in the browser (ABI v2)".into()),
            other => Err(format!("unknown program '{other}'")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Blastp => "blastp",
            Self::Tblastn => "tblastn",
            Self::Blastn => "blastn",
            Self::Tblastx => "tblastx",
            Self::Blastx => "blastx",
        }
    }

    /// The program in the library API's terms.
    pub fn api(self) -> losat::api::Program {
        match self {
            Self::Blastp => losat::api::Program::Blastp,
            Self::Tblastn => losat::api::Program::Tblastn,
            Self::Blastn => losat::api::Program::Blastn,
            Self::Tblastx => losat::api::Program::Tblastx,
            Self::Blastx => losat::api::Program::Blastx,
        }
    }

    /// The output formats that the engine implements for this program; a run writes
    /// exactly these, each on the stream of the same number.
    pub fn formats(self) -> &'static [u32] {
        match self {
            Self::Blastp | Self::Tblastn | Self::Blastn | Self::Tblastx | Self::Blastx => {
                &[0, 6, 7]
            }
        }
    }
}

/// One output of a run: a format that every run of the program writes (`Program::formats`), on
/// the stream of its number, or the output of the argv's own `-outfmt` (any other format that the
/// CLI writes, or 6 and 7 with fields), on `STREAM_CUSTOM` (docs/web/abi_v2.md §5, §7).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Output {
    pub stream: u32,
    /// The format's number, which sets its hit list size.
    pub number: u32,
    /// The `-outfmt` value that the engine writes it for.
    pub outfmt: String,
    /// Whether this is the argv's own format (`Requested`): LOSAT's refusal of its report fails
    /// the run, as the CLI's run of the argv fails. A refusal of another output's report fails
    /// that output alone (`FormatRefusal::FailsAlone`, `search`).
    pub requested: bool,
}

impl Output {
    fn format(number: u32) -> Self {
        Self {
            stream: number,
            number,
            outfmt: number.to_string(),
            requested: false,
        }
    }
}

/// The argv's `-outfmt` as a run reads it (`parse_run`): the `-outfmt` value of the argv's own
/// output (`STREAM_CUSTOM`), if it has one, and the number of the argv's format. Without an output
/// of its own the argv's format is the program's format of that number (0, 6 or 7). No number: a
/// value that does not parse, of an argv that names a search strategy file, which the command
/// line's early reading leaves to the program's formatting handler
/// (`losat::web_adapter::early_format_choice`); it is the run's one output and the arguments'
/// `-outfmt`, which the engine's handlers read again and reject after the arguments' own errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Requested {
    pub custom: Option<String>,
    pub number: Option<u32>,
}

/// The outputs of a run: the output of the argv's `-outfmt` (`parse_run`) if it has one, then
/// the program's formats; the argv's format among them is `Output::requested`. The first format
/// of a search writes the warnings that the formats share (the formatting handler's, the
/// queries'), and the output of the `-outfmt` adds those of its fields' writer (NCBI's taxdb
/// warning): the argv's format comes first in its search (`search_groups`), so that stream 3
/// starts with the CLI's standard error for the argv (BLASTX writes the warnings of its first
/// format only).
///
/// A value left to the formatting handler (`Requested::number` none) is the run's one output,
/// numbered 0 (the CLI's default; its search stops at the formatting handler at the latest, before
/// a reader of the number), as the CLI's run of the argv has that format alone: no search of the
/// program's formats writes before its error.
pub fn outputs(program: Program, requested: &Requested) -> Vec<Output> {
    let mut outputs = Vec::new();
    if let Some(outfmt) = &requested.custom {
        outputs.push(Output {
            stream: STREAM_CUSTOM,
            number: requested.number.unwrap_or(0),
            outfmt: outfmt.clone(),
            requested: true,
        });
    }
    let Some(requested_number) = requested.number else {
        return outputs;
    };
    outputs.extend(program.formats().iter().map(|&number| Output {
        requested: requested.custom.is_none() && number == requested_number,
        ..Output::format(number)
    }));
    outputs
}

/// The `-outfmt` value of the program's arguments.
fn outfmt_of(command: &mut Commands) -> &mut std::ffi::OsString {
    match command {
        Commands::Blastp(a) => &mut a.outfmt,
        Commands::Tblastn(a) => &mut a.outfmt,
        Commands::Blastn(a) => &mut a.outfmt,
        Commands::Tblastx(a) => &mut a.outfmt,
        Commands::Blastx(a) => &mut a.outfmt,
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {
            unreachable!("rejected by Program::parse")
        }
    }
}

/// The engine's name of the program, for its reading of `-outfmt`.
fn engine_program(program: Program) -> losat::web_adapter::EngineProgram {
    use losat::web_adapter::EngineProgram as Engine;
    match program {
        Program::Blastp => Engine::Blastp,
        Program::Tblastn => Engine::Tblastn,
        Program::Blastn => Engine::Blastn,
        Program::Tblastx => Engine::Tblastx,
        Program::Blastx => Engine::Blastx,
    }
}

/// Parses an argv (program name first) with the engine's CLI parser (`parse_run`), without the
/// output of its `-outfmt`.
pub fn parse(words: &[&str]) -> Result<(Program, Commands), String> {
    parse_run(words).map(|(program, command, _)| (program, command))
}

/// Parses an argv (program name first) with the engine's CLI parser, and reads its `-outfmt` as
/// the CLI does (`parse_formatting_string`, `report_format`, with the CLI's messages; the
/// formatting handler's check of a custom delimiter comes later, `custom_handler_check`). A run
/// writes the program's formats 0, 6 and 7 whatever the value; every other format that the CLI
/// writes (5, 8 to 16, 18, 20, BLASTN's 17), and 6 and 7 with words after the number (fields,
/// `delim=`), add an output of their own, the third value returned (docs/web/abi_v2.md §7). The
/// rule reads the words, not the fields they make, so that the application reads it from the argv
/// alike: `6 std` adds an output with the bytes of outfmt 6. 0 with words adds none (the CLI writes
/// outfmt 0); 1 to 4 add one without the words. The output of outfmt 13 and 14 is a
/// ZIP archive of the master file and the side files (`search`). The arguments come back with
/// `-outfmt` at the CLI's default, `0`: the engine is given the formats of each search
/// (`ReportOutputs`), and `validate` gives the checks of the arguments alone (`check_options`)
/// the argv's own format.
pub fn parse_run(words: &[&str]) -> Result<(Program, Commands, Requested), String> {
    let (program, command, supplied) = parse_imported(words)?;
    read_outfmt(words, program, command, supplied.is_some())
}

/// The first part of `parse_run`: the adapter's own rules and the engine's CLI parser, without the
/// reading of `-outfmt`, which NCBI makes in `Run`, after the settings and a dry run's end
/// (`started_command`). An argv whose reader stopped at a help word fails with the text of the
/// help (`parse_printing` tells it apart).
fn parse_words(words: &[&str]) -> Result<(Program, Commands), String> {
    match parse_printing(words, None)? {
        Parsed::Command(program, command, _) => Ok((program, *command)),
        Parsed::Prints(text) => Err(text),
    }
}

/// What the engine's CLI parser makes of an argv (`parse_printing`).
enum Parsed {
    /// The program's arguments (boxed: they are far larger than a text), and with an
    /// `-import_search_strategy NAME` the run's import of the registered strategy file, which has
    /// what it carries as inputs (`Supplied::take`).
    Command(Program, Box<Commands>, Option<Supplied>),
    /// The text that the command line prints to its standard output in place of a run, where
    /// NCBI's reader stopped at a help word (`-h`, `-help`, `-help-full`, `-xmlhelp`, the words
    /// before it read: `losat::web_adapter::ArgvError::Printed`).
    Prints(String),
}

/// `parse_words`, with the help that stopped the reader as its own outcome (`Parsed::Prints`; line
/// WQA's D14, 2026-11-20). The help is printed as the command line prints it, whatever the words
/// around it (`-out` too: line WQA's D14 review, 2026-11-20); any other argv with `-out` is refused.
/// With an `-import_search_strategy NAME` the command line's import of the registered strategy file
/// `NAME` (kind 3, `files::strategy`): the arguments are the import's
/// (`losat::web_adapter::import_strategy`: the engine's import, the path of the command line),
/// and what the strategy carries as inputs comes with them; the program is the strategy's where
/// the argv's program runs another program's strategy as the command line does. A file that is
/// not registered is the command line's file that is not accessible.
/// A run gives its start (`started_command`): the posts of an argument error go to the run's
/// `-logfile` target, as the command line's (`try_parse_cli`).
fn parse_printing(words: &[&str], started: Option<&CallerStart>) -> Result<Parsed, String> {
    // The formats come from `Program::formats`, and every program accepts the CLI's
    // default `-outfmt 0`, so the host's words are parsed exactly as on the command line,
    // with the CLI's binary name: an error (an unknown program too) is the CLI's message.
    let argv = std::iter::once("losat").chain(words.iter().copied());
    let parsed = losat::web_adapter::try_parse_cli(argv, started);
    if let Err(losat::web_adapter::ArgvError::Printed(text)) = &parsed {
        return Program::parse(words.first().copied().unwrap_or(""))
            .map(|_| Parsed::Prints(text.clone()));
    }
    for word in words.iter().skip(1) {
        let key = word.split('=').next().unwrap_or("");
        if matches!(key, "-out" | "--out") {
            return Err(format!(
                "{key} is not accepted: a run writes every output format of the program"
            ));
        }
    }
    let cli: Cli = parsed.map_err(|error| losat::web_adapter::argv_message(&error))?;
    let program = Program::parse(words.first().copied().unwrap_or(""))?;
    let command = cli.command;
    let Some(name) =
        losat::web_adapter::strategy_file_name(&command).map(std::path::Path::to_path_buf)
    else {
        return Ok(Parsed::Command(program, Box::new(command), None));
    };
    let bytes = crate::files::strategy(&name).ok_or_else(|| {
        format!(
            "Command line argument error: Argument \"import_search_strategy\". File is not accessible:  `{}'",
            name.display()
        )
    })?;
    // The import is the run's, made once (`imported`): a reading of the same argv over the same
    // file takes a copy of its arguments. The strategy's program reads the inputs and runs
    // (`losat::web_adapter::program_of`).
    let supplied = Supplied::of(words, program, command, &bytes)?;
    Ok(Parsed::Command(
        supplied.program(),
        Box::new(supplied.command()),
        Some(supplied),
    ))
}

/// `parse` with the inputs that the argv's search strategy supplies (`parse_printing`).
pub fn parse_imported(words: &[&str]) -> Result<(Program, Commands, Option<Supplied>), String> {
    match parse_printing(words, None)? {
        Parsed::Command(program, command, supplied) => Ok((program, *command, supplied)),
        Parsed::Prints(text) => Err(text),
    }
}

/// What a run of `words` (program name first) prints in place of a search (`STREAM_PRINTED`), as the
/// command line does: `"version"` for a version word (`-version` or a full version), which NCBI's
/// application reads before the arguments, `"help"` where the arguments' reader stops at a help
/// word (`-xmlhelp` included); `None` for an argv that runs (or fails). Read without the
/// configuration, which does not change it.
pub fn prints(words: &[&str]) -> Option<&'static str> {
    use losat::web_adapter::Preparsed;
    match losat::web_adapter::preparse(words.get(1..).unwrap_or(&[])) {
        Preparsed::Version(_) => Some("version"),
        Preparsed::Words { .. } => {
            matches!(parse_printing(words, None), Ok(Parsed::Prints(_))).then_some("help")
        }
    }
}

/// Whether the argv of `command` names its subjects: `-subject` or `-db`.
fn names_subjects(command: &Commands) -> bool {
    match command {
        Commands::Blastp(a) => a.subject.is_some() || a.db.is_some(),
        Commands::Tblastn(a) => a.subject.is_some() || a.db.is_some(),
        Commands::Blastn(a) => a.subject.is_some() || a.db.is_some(),
        Commands::Tblastx(a) => a.subject.is_some() || a.db.is_some(),
        Commands::Blastx(a) => a.subject.is_some() || a.db.is_some(),
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {
            unreachable!("rejected by Program::parse")
        }
    }
}

/// The second part of `parse_run`: the argv's `-outfmt`, read as the CLI reads it, through the
/// CLI's early reading (`losat::web_adapter::early_format_choice`), with the search strategy file
/// that the argv names (`-export_search_strategy`, which the run makes before it reads `-outfmt`,
/// `started_command`): with one, a value that does not parse is left to the formatting handler
/// (`Requested`).
fn read_outfmt(
    words: &[&str],
    program: Program,
    mut command: Commands,
    imported: bool,
) -> Result<(Program, Commands, Requested), String> {
    use losat::web_adapter::EarlyFormat;
    let strategy = losat::web_adapter::export_search_strategy_of(&command).is_some();
    let outfmt = std::mem::replace(outfmt_of(&mut command), std::ffi::OsString::from("0"))
        .into_string()
        .map_err(|_| "invalid UTF-8 in -outfmt".to_string())?;
    let early = losat::web_adapter::early_format_choice(&outfmt, strategy)
        .map_err(|error| format!("{error:#}"))?;
    // NCBI's error for an argv without `-subject` (and without `-db`), which the command line
    // raises after it read `-outfmt`'s number and before its other checks, as it does where the
    // run has inputs (the registered records are the files of `-query` and `-subject`; line
    // WQA's D14, 2026-11-21: `blastp -query -h`, whose `-h` is the query's name). A search strategy
    // gives the subjects or the database itself (its import has refused one that gives none).
    if !imported && !names_subjects(&command) {
        return Err(format!("{:#}", losat::web_adapter::missing_subject_error()));
    }
    let choice = match early {
        EarlyFormat::Parsed(choice) => choice,
        // The value stays the arguments' and is the run's one output, which the engine's
        // formatting handler reads again (line RG F1w review 2, S-a).
        EarlyFormat::Deferred => {
            *outfmt_of(&mut command) = std::ffi::OsString::from(&outfmt);
            let requested = Requested {
                custom: Some(outfmt),
                number: None,
            };
            return Ok((program, command, requested));
        }
    };
    // The `@PG` line of SAM (outfmt 17) shows the command line, which `main` gives BLASTN: the
    // run's, `losat` and its argv, each word followed by one space (`cli::command_line_bytes`). Only
    // a run that writes SAM holds it.
    if let (Commands::Blastn(args), losat::web_adapter::SAM_FORMAT) = (&mut command, choice.number)
    {
        let mut line = b"losat ".to_vec();
        for word in words {
            line.extend_from_slice(word.as_bytes());
            line.push(b' ');
        }
        args.command_line = Some(line);
    }
    // What the CLI rejects, with its message: a format that LOSAT's program does not write, and the
    // formats that NCBI's formatting handler rejects (`report_format` gives none: SAM of another
    // program, AIRR, FASTA), which a run never gives the engine; the CLI rejects those after the
    // ranges and the query options, so for an argv with an error of both the CLI reports the
    // other. The CLI is taken with `-out`, which the formats 13 and 14 need there (a run gives
    // them a receiver of their side files, `search`). The handler's
    // check of a custom delimiter comes with the program's checks of the options (`validate`) and
    // in the engine's `run_local`.
    let engine = engine_program(program);
    losat::web_adapter::report_format(&choice, engine, Some(std::path::Path::new("out")))
        .and_then(|format| match format {
            Some(_) => Ok(()),
            None => losat::web_adapter::formatting_handler_check(&choice, engine),
        })
        .map_err(|error| format!("{error:#}"))?;
    // NCBI's white space around the value (`NStr::TruncateSpaces`), as `parse_formatting_string`
    // reads it: a space inside is followed by words.
    let words_follow = outfmt
        .trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r'))
        .contains(' ');
    let custom = match choice.number {
        0 => None,
        6 | 7 => words_follow.then_some(outfmt),
        // The query-anchored formats (every program's, line RG T3) without the words after the
        // number, which NCBI and the CLI drop (`FormatChoice::normalized`). The engine's library
        // drops them as well (`parse_formatting_string`; BLASTP's since line RG F/d36), so this
        // only gives the output one spelling of its number.
        1..=4 => Some(
            String::from_utf8(choice.normalized())
                .map_err(|_| "invalid UTF-8 in -outfmt".to_string())?,
        ),
        _ => Some(outfmt),
    };
    let number = Some(choice.number as u32);
    Ok((program, command, Requested { custom, number }))
}

/// Whether a run of `words` (program name first) takes its query and subject from
/// `store::register_inputs`, which reads them together for its argv, and not from `register`: an
/// argv with `-import_search_strategy`, whose inputs may come from the strategy, or with
/// `-parse_deflines` (`parse_deflines_of`), whose reading of the inputs needs the argv
/// (the ids of the deflines, read in one scope of ids for the query and the subjects), that is not
/// a dry run, which reads no input. False for an argv that the CLI's parser rejects.
pub fn inputs_together(words: &[&str]) -> bool {
    let dry_run = matches!(
        losat::web_adapter::preparse(words.get(1..).unwrap_or(&[])),
        losat::web_adapter::Preparsed::Words { dryrun: true, .. }
    );
    // A strategy's queries and subjects are read here too, with the argv's options.
    !dry_run
        && parse_words(words).is_ok_and(|(_, command)| {
            losat::web_adapter::parse_deflines_of(&command)
                || losat::web_adapter::strategy_file_name(&command).is_some()
        })
}

/// Whether the run of `words` writes the BLAST archive (outfmt 11) or a search strategy file
/// (`-export_search_strategy`), whose Bioseqs show the whole records, as the command line's reader
/// keeps them for either (`store::with_inputs`).
pub fn keeps_whole_records(words: &[&str]) -> bool {
    parse_run(words).is_ok_and(|(_, command, requested)| {
        losat::web_adapter::export_search_strategy_of(&command).is_some()
            || requested.custom.is_some_and(|outfmt| {
                losat::web_adapter::parse_formatting_string(&outfmt)
                    .is_ok_and(|choice| choice.number == losat::web_adapter::ARCHIVE_FORMAT)
            })
    })
}

/// The run's working directory: the name under which the host preopens it (docs/web/abi_v2.md §3),
/// and the module's current directory. A relative name of the argv (`-db db/nt`) is a file of it, as
/// on the command line run from that directory; an absolute one (`/db/nt`) names a file outside it,
/// where the module has none, as the command line has none at its root. (Under a preopen `.` and the
/// current directory `/`, WASI's C library gives both names to the same file; review of D12a, low 4.)
pub const RUN_DIRECTORY: &str = "/work";

/// Makes [`RUN_DIRECTORY`] the module's current directory, once per instance, before its first call
/// that may read a file. An instance whose host preopened no run directory keeps `/`: its runs find
/// no file, as the command line finds none that is not there.
pub fn enter_run_directory() {
    #[cfg(target_os = "wasi")]
    {
        static ENTERED: std::sync::Once = std::sync::Once::new();
        ENTERED.call_once(|| {
            let _ = std::env::set_current_dir(RUN_DIRECTORY);
        });
    }
}

/// The database of the arguments `command` (`-db NAME`), which the run opens by name as the command
/// line opens it, from the files of the run's working directory that the host gives the engine
/// (docs/web/abi_v2.md §3, §4; line WQA D12a, 2026-11-20). `None` for a search of `-subject`.
fn database_of(command: &Commands) -> Option<&std::ffi::OsStr> {
    match command {
        Commands::Blastp(args) => args.db.as_deref(),
        Commands::Tblastn(args) => args.db.as_deref(),
        Commands::Blastn(args) => args.db.as_deref(),
        Commands::Tblastx(args) => args.db.as_deref(),
        Commands::Blastx(args) => args.db.as_deref(),
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => None,
    }
}

/// Whether a run of `words` (program name first) searches the database of `-db` (`database_of`):
/// its `run` takes no subject handle (0), and `register_inputs` reads its query alone. False for
/// an argv that the CLI's parser rejects.
pub fn searches_database(words: &[&str]) -> bool {
    parse_words(words).is_ok_and(|(_, command)| database_of(&command).is_some())
}

/// The environment of the application start of a run in LOSAT Web (`start`): no variable but the
/// directory of BLASTN's `-window_masker_taxid` (`files::WINDOW_MASKER_PATH`), set in the layer of
/// NCBI's registry that comes before its files, so that the registry of `-window_masker_taxid`
/// (`.ncbirc` and the environment) gives that directory whatever a registered `.ncbirc` sets.
const START_VARIABLES: &[(&str, &str)] = &[(
    "NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH",
    crate::files::WINDOW_MASKER_PATH,
)];

/// The NCBI application start of a run of `words` (program name first), as the command line's
/// before it reads the arguments (`losat::web_adapter::start_run`): in the environment
/// `START_VARIABLES`, with the registry files registered now (kind 2, `files::registry_files`) and
/// the file of `-logfile` going to `run_files`; the messages that the command line writes to its
/// standard error go to `diagnostics`. `Err` when the start ends the run (a usage report or a
/// load that LOSAT refuses, the error of NCBI's load:
/// `-conffile` of a file that is not registered, a registry file that NCBI cannot read), with the
/// message that the command line writes to its standard error without `-logfile`, and the start
/// when its file of `-logfile` is open (it has the messages); and for a program that ABI v2 does
/// not run, with the error of its argv (`parse`: an unknown program's is the command line's).
fn start(
    words: &[&str],
    registry: &crate::files::RegisteredRegistry,
    run_files: &RunFiles,
    diagnostics: &mut dyn std::io::Write,
) -> Result<CallerStart, (String, Option<Box<CallerStart>>)> {
    let (name, arguments) = words.split_first().unwrap_or((&"", &[]));
    if let Err(error) = Program::parse(name) {
        return Err((parse(words).err().unwrap_or(error), None));
    }
    let start = losat::web_adapter::start_run(
        name,
        arguments,
        START_VARIABLES,
        registry,
        run_files,
        diagnostics,
    )
    .map_err(|message| (message, None))?;
    match start.ended() {
        Some(message) => Err((message.to_string(), Some(Box::new(start)))),
        None => Ok(start),
    }
}

/// The file names of a search run (`search_run_files`): those the host holds outside `/work`
/// (`registered`: `FileReceiver::registered_files`) and those in `/work` the run opens for reading
/// whether present or missing (`read`: `FileReceiver::read_files`, a `-db` search's list files).
#[derive(Default)]
struct SearchRunFiles {
    registered: Vec<Vec<u8>>,
    read: Vec<Vec<u8>>,
}

/// The names of the search's files that the host holds outside `/work` and the files of a `-db`
/// or `-filtering_db` search that the run opens or checks for reading (line WQA D19 rounds R1 and
/// R2, D19b):
/// - `registered`: `-query` and `-subject` (unless `-`), BLASTN's `-window_masker_db` or
///   `-window_masker_taxid` (`windowmasker/<taxid>/wmasker.obinary`) when registered as kind 1,
///   and the kind-2 registry files of `registry` that this run's start reads (`.ncbirc`, and
///   without a version word either the registered file of `-conffile` or `<program>.ini`);
/// - `read`: with `-db`, the list files of its filters (`-gilist`, `-negative_gilist`,
///   `-ipglist`, `-negative_ipglist`, `-seqidlist`, `-negative_seqidlist`, `-taxidlist`,
///   `-negative_taxidlist`), the database lookup files (`search_database_read_files`, including
///   BLASTN's `-filtering_db` even with `-subject`) and taxonomy lookup files, kept apart from
///   `registered` so a missing file does not make its parent folders present.
fn search_run_files(words: &[&str], registry: &crate::files::RegisteredRegistry) -> SearchRunFiles {
    use losat::web_adapter::{
        preparse, scan_program_argv, search_database_read_files, Preparsed, RegistryFiles as _,
    };
    use std::path::Path;

    let (program_name, arguments) = words.split_first().unwrap_or((&"", &[]));
    let Ok(program) = Program::parse(program_name) else {
        return SearchRunFiles::default();
    };

    fn push_input(names: &mut Vec<Vec<u8>>, value: &str) {
        if !value.is_empty() && value != "-" {
            names.push(value.as_bytes().to_vec());
        }
    }
    fn push_counts(names: &mut Vec<Vec<u8>>, value: &str) {
        if !value.is_empty() && crate::files::registered(value) {
            names.push(value.as_bytes().to_vec());
        }
    }
    fn push_taxid(names: &mut Vec<Vec<u8>>, value: &str) {
        if let Ok(taxid) = value.parse::<i32>() {
            if taxid != 0 {
                let path = format!(
                    "{}/{taxid}/wmasker.obinary",
                    crate::files::WINDOW_MASKER_PATH
                );
                push_counts(names, &path);
            }
        }
    }

    // The command line opens `-logfile` before it reads the arguments, so the run's input files do
    // not depend on whether the argv parses. Walk `arguments` once through the program's own
    // argument definitions (`scan_program_argv`), skipping the value of every option that takes
    // one (even when it starts with `-`).
    let mut registered: Vec<Vec<u8>> = Vec::new();
    let mut list_files: Vec<Vec<u8>> = Vec::new();
    let mut dbs: Vec<String> = Vec::new();
    let mut filtering_dbs: Vec<String> = Vec::new();
    let mut outfmt: Option<String> = None;
    let mut has_taxids = false;
    let mut no_taxid_expansion = false;
    let is_blastn = program == Program::Blastn;
    scan_program_argv(program_name, arguments, |name, value| {
        if name == "no_taxid_expansion" && value.is_none() {
            no_taxid_expansion = true;
            return;
        }
        let Some(value) = value else {
            return;
        };
        match name {
            "query" | "subject" => push_input(&mut registered, value),
            "window_masker_db" if is_blastn => push_counts(&mut registered, value),
            "window_masker_taxid" if is_blastn => push_taxid(&mut registered, value),
            "db" => dbs.push(value.to_string()),
            "filtering_db" if is_blastn => filtering_dbs.push(value.to_string()),
            "outfmt" => outfmt = Some(value.to_string()),
            "taxids" | "negative_taxids" => has_taxids = true,
            "taxidlist" | "negative_taxidlist" => {
                has_taxids = true;
                if !value.is_empty() {
                    list_files.push(value.as_bytes().to_vec());
                }
            }
            "gilist" | "negative_gilist" | "ipglist" | "negative_ipglist" | "seqidlist"
            | "negative_seqidlist" => {
                if !value.is_empty() {
                    list_files.push(value.as_bytes().to_vec());
                }
            }
            _ => {}
        }
    });
    if registry.is_file(Path::new(".ncbirc")) {
        registered.push(b".ncbirc".to_vec());
    }
    if let Preparsed::Words { conffile, .. } = preparse(arguments) {
        match conffile {
            Some(name) => {
                if registry.is_file(Path::new(&name)) {
                    registered.push(name.as_encoded_bytes().to_vec());
                }
            }
            None => {
                let ini = format!("{program_name}.ini");
                if registry.is_file(Path::new(&ini)) {
                    registered.push(ini.into_bytes());
                }
            }
        }
    }
    let mut read = if !dbs.is_empty() {
        list_files
    } else {
        Vec::new()
    };
    let db_refs: Vec<&str> = dbs.iter().map(String::as_str).collect();
    let filtering_db_refs: Vec<&str> = filtering_dbs.iter().map(String::as_str).collect();
    read.extend(search_database_read_files(
        program.api(),
        arguments,
        registry,
        &db_refs,
        &filtering_db_refs,
        outfmt.as_deref().unwrap_or("0"),
        has_taxids && !no_taxid_expansion,
    ));
    SearchRunFiles { registered, read }
}

/// The names of the search's files that the host holds outside `/work` (`search_run_files`).
#[cfg(test)]
fn registered_run_files(
    words: &[&str],
    registry: &crate::files::RegisteredRegistry,
) -> Vec<Vec<u8>> {
    search_run_files(words, registry).registered
}

/// What a run does after its start (`started_command`).
enum Next {
    /// It prints this text to the command line's standard output and ends, exit status 0.
    Prints(String),
    /// It ends, exit status 0: a dry run.
    Ends,
    /// It searches (boxed: the arguments are far larger than a text), with the file of
    /// `-export_search_strategy` that the run made (`losat::web_adapter::open_strategy_file`).
    Search(Box<(Program, Commands, Requested)>, Option<StrategyFile>),
}

/// The argv `words` after the run's start `started`, in the command line's order (`cli::main_entry`):
/// the version that the start printed (`CallerStart::version`), the arguments parsed
/// (`parse_printing`: the help that stops the reader is printed), the settings that LOSAT refuses
/// (`start.settings()?`), the end of a dry run (rc 0), the file of `-export_search_strategy`, made
/// through the run's files `run_files` (`losat::web_adapter::open_strategy_file`, NCBI's
/// `GetExportSearchStrategyStream` in `Run`), then `-outfmt`, which NCBI reads after it
/// (`read_outfmt`: a dry run does not read it), with the error of an argv without subjects.
/// `validate` and `run` take the argv here. With `to_log` (`run`) the posts of an argument error
/// go to the run's `-logfile` target, which the run hands over (`CallerStart::finish`), and the
/// error is the command line's standard error, NCBI's line that points to the error logs included,
/// as on the command line; without (`validate`, whose start keeps no file) they stay in the error,
/// the report as the command line writes it without `-logfile`, so that the form names the
/// refused fields.
fn started_command(
    words: &[&str],
    started: &CallerStart,
    run_files: &RunFiles,
    to_log: bool,
) -> Result<Next, String> {
    if let Some(text) = started.version() {
        return Ok(Next::Prints(text.to_string()));
    }
    let (program, command, imported) = match parse_printing(words, to_log.then_some(started))? {
        Parsed::Prints(text) => return Ok(Next::Prints(text)),
        Parsed::Command(program, command, supplied) => (program, *command, supplied.is_some()),
    };
    started.settings()?;
    if started.dry_run() {
        return Ok(Next::Ends);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627 (cited at
    // `strategy::open_strategy_file`): the file is made before `-outfmt` is read; an error after it
    // leaves it empty (the receiver gets the file when it is dropped).
    let strategy = losat::web_adapter::open_strategy_file(&command, &words[1..], run_files)
        .map_err(|error| format!("{error:#}"))?;
    read_outfmt(words, program, command, imported)
        .map(|parsed| Next::Search(Box::new(parsed), strategy))
}

/// A receiver of the files of a run that keeps none (`validate`, which emits nothing), carrying
/// the search's registered and read file names so `-logfile` is checked against the run's files.
struct NoFiles(SearchRunFiles);

impl losat::web_adapter::FileReceiver for NoFiles {
    fn receive(&self, _name: &[u8], _bytes: &[u8]) -> std::io::Result<()> {
        Ok(())
    }

    fn registered_files(&self) -> &[Vec<u8>] {
        &self.0.registered
    }

    fn read_files(&self) -> &[Vec<u8>] {
        &self.0.read
    }
}

/// The *start* JSON of `losat_web2_start` (docs/web/abi_v2.md §4, §9): what NCBI's application
/// reads of the argv `words` (program name first) before its arguments (`preparse`): the file of
/// the last `-conffile` (which the host registers as kind 2 before it validates; none for a
/// version word, which the command line prints before it loads the file), the file of `-logfile`
/// (NCBI's scan of it, `log_file`, and the command line's normalized file name under which a run
/// gives it on `STREAM_RUN_FILES`, `log_file_name`, `RunFiles::log_host_name`), the file of
/// `-export_search_strategy` (the name as the argv's parser reads it, `strategy_file`: a run gives
/// it on `STREAM_RUN_FILES` too; none for an argv that the parser rejects), whether the run is a dry
/// run, whether the run takes its inputs from `losat_web2_register_inputs` (`inputs_together`: an
/// argv with `-parse_deflines`), whether it prints a text in place of a search (`prints`:
/// `"help"` or `"version"`, `STREAM_PRINTED`), which needs no inputs, and whether it searches the
/// database of `-db` (`searches_database`: its `run` takes no subject handle, and the engine opens
/// the database from the run's files, line WQA D12a). An unknown program is an error.
pub fn start_json(words: &[&str]) -> Result<String, String> {
    use losat::web_adapter::Preparsed;
    Program::parse(words.first().copied().unwrap_or(""))?;
    let arguments = &words[1..];
    let (conffile, dryrun) = match losat::web_adapter::preparse(arguments) {
        Preparsed::Words {
            conffile, dryrun, ..
        } => (conffile, dryrun),
        Preparsed::Version(_) => (None, false),
    };
    let prints = prints(words);
    let log_file = losat::web_adapter::log_file(arguments);
    let log_file_name = log_file
        .and_then(RunFiles::log_host_name)
        .and_then(|bytes| String::from_utf8(bytes).ok());
    let mut out = String::from("{\"conffile\":");
    json::optional_string(&mut out, conffile.as_deref().and_then(|name| name.to_str()));
    out.push_str(",\"log_file\":");
    json::optional_string(&mut out, log_file.and_then(|name| name.to_str()));
    out.push_str(",\"log_file_name\":");
    json::optional_string(&mut out, log_file_name.as_deref());
    out.push_str(",\"strategy_file\":");
    json::optional_string(&mut out, strategy_file(words).as_deref());
    // A run that prints its help ends before the end of a dry run (the reader stops first).
    let _ = write!(
        out,
        ",\"dry_run\":{},\"register_inputs\":{},\"prints\":",
        dryrun && prints.is_none(),
        inputs_together(words)
    );
    json::optional_string(&mut out, prints);
    let _ = write!(out, ",\"database\":{}}}", searches_database(words));
    Ok(out)
}

/// The name of the file of `-export_search_strategy` of the argv `words` (program name first), as
/// the run gives the file (`losat::web_adapter::open_strategy_file`): the name as the argv's parser
/// reads it. `None` without the option, for a name that is not UTF-8 (which the run refuses) and
/// for an argv that the parser rejects or that prints a text.
fn strategy_file(words: &[&str]) -> Option<String> {
    let (_, command) = parse_words(words).ok()?;
    losat::web_adapter::export_search_strategy_of(&command)
        .and_then(|path| path.to_str())
        .map(str::to_string)
}

/// `losat_web2_start`: emits the *start* JSON (`start_json`) on stream 2.
pub fn emit_start(words: &[&str]) -> Result<(), String> {
    let response = start_json(words)?;
    emit(STREAM_RESPONSE, response.as_bytes());
    Ok(())
}

/// Whether the run of `words` reads its FASTA inputs with NCBI's data loaders on
/// (`ApplicationSettings::data_loaders`: `[BLAST] DATA_LOADERS` of the registry, which `-conffile`
/// names): with them on, the reader refuses a first line that NCBI would fetch as a Seq-id; with
/// `DATA_LOADERS=none` it reads that line as data, as the CLI does (line WQA's D16, 2026-11-20).
/// `Err` with the start's error where the start, or LOSAT's check of its settings, ends the run.
pub fn data_loaders(words: &[&str]) -> Result<bool, String> {
    let registry = crate::files::registry_files();
    let run_names = search_run_files(words, &registry);
    let started = start(
        words,
        &registry,
        &RunFiles::Receiver(Arc::new(NoFiles(run_names))),
        &mut std::io::sink(),
    )
    .map_err(|(message, _)| message)?;
    Ok(started.settings()?.data_loaders)
}

/// `validate`: the application start (`start`: its errors are the command line's first), then
/// `parse`, the settings that LOSAT refuses, and nothing more for a dry run (`-dryrun`); then the
/// syntax of `-subject_loc` and `-query_loc`, BLASTN's
/// `-dust` value and the scoring options that NCBI rejects or that its Karlin-Altschul
/// tables do not support, with NCBI's message (for a batch of one query), and TBLASTX's
/// checks of its options alone (NCBI's `-evalue` check, the options that LOSAT's TBLASTX
/// rejects).
/// The checks include the hit list size, which differs between the searches of a run
/// (`search_groups`): they are those of the search of the argv's own format (`own_search`; the
/// format itself for BLASTP, TBLASTN and TBLASTX), whose outcome is the run's (`run`,
/// `search_all`), as the CLI checks the argv's `-outfmt` alone.
/// The other search's rejection fails only its formats, which the argv does not ask for.
/// `run` leaves the other checks to the engine, which reports them as the CLI does.
/// Returns whether the run reads its inputs with NCBI's data loaders on (`data_loaders`), which
/// the host gives `register` and `scan_begin` (line WQA's D16, 2026-11-20).
pub fn validate(words: &[&str]) -> Result<bool, String> {
    validate_with(words, &mut |_| ())
}

/// `validate`, which also gives `known` the data loaders of the argv's registry as soon as the
/// start has read it and LOSAT takes its settings, before the checks of the arguments and the
/// options: the setting does not depend on them, so it is known while another option is refused
/// (line WQA's D13 round 3, L1, 2026-11-20). Not called when the start or the settings fail.
pub fn validate_with(words: &[&str], known: &mut dyn FnMut(bool)) -> Result<bool, String> {
    use losat::web_adapter::{parse_optional_range, RangeRole};
    let registry = crate::files::registry_files();
    let run_names = search_run_files(words, &registry);
    let no_files = RunFiles::Receiver(Arc::new(NoFiles(run_names)));
    let started =
        start(words, &registry, &no_files, &mut std::io::sink()).map_err(|(message, _)| message)?;
    if let Ok(settings) = started.settings() {
        known(settings.data_loaders);
    }
    let (program, command, requested) = match started_command(words, &started, &no_files, false)? {
        Next::Search(parsed, _) => *parsed,
        // A run that prints reads no input; the command line prints before it refuses settings.
        Next::Prints(_) => {
            return Ok(started
                .settings()
                .is_ok_and(|settings| settings.data_loaders))
        }
        Next::Ends => return Ok(started.settings()?.data_loaders),
    };
    let custom = requested.custom.as_deref();
    let settings = started.settings()?;
    // The ranges' syntax, in NCBI's order: the subject range where the subjects are read
    // (before the other option handlers), the query range with the query options (the
    // engine's option checks read it). Whether a range fits the records is checked by
    // `run`, which has them.
    let subject_range = |subject_loc: Option<&std::ffi::OsStr>| {
        parse_optional_range(subject_loc, RangeRole::Subject).map(|_| ())
    };
    let groups = search_groups(&outputs(program, &requested), &command);
    let own = &groups[own_search(&groups)];
    let custom_check = || custom_handler_check(program, custom);
    // A custom delimiter that a field uses makes the CLI reject the argv, whatever the hit list
    // sizes of the searches: its error is that of the CLI's checks with the argv's `-outfmt`, which
    // make the handler's check where the CLI does (the error may come before it).
    let custom_rejected = custom_check().is_err();
    // BLASTP, TBLASTN and TBLASTX check their options with one `-outfmt`, which sets the hit list
    // size: the argv's own format, first in its search (`search_groups`), as the engine's search
    // checks it and as the CLI does, whose hit list check (`No hits are being saved`) comes before
    // the e-value's (line WQA's sweep, D3 review r2, 2026-11-19: with format 0's size, `-outfmt 6
    // -num_alignments 0 -evalue 0` gave the e-value's error).
    let with_own = |outfmt: &mut std::ffi::OsString| {
        *outfmt = std::ffi::OsString::from(&own[0].outfmt);
    };
    // The CLI checks the thread count before the options (`validate_threads`), the checks of the
    // options alone last; with the argv's `-outfmt` they stop at its delimiter first.
    if custom_rejected {
        losat::web_adapter::validate_threads(num_threads(&command))
            .map_err(|error| format!("{error:#}"))?;
    }
    match command {
        Commands::Blastn(mut args) => subject_range(args.subject_loc.as_deref())
            .and_then(|()| args.resolve_dust())
            .and_then(|()| {
                if validates_window_masker(&args) {
                    give_blastn_counts(&mut args, &crate::files::source(), settings);
                    args.resolve_window_masker()
                } else {
                    Ok(())
                }
            })
            .and_then(|()| {
                parse_optional_range(args.query_loc.as_deref(), RangeRole::Query).map(|_| ())
            })
            .map_err(|error| format!("{error:#}"))
            // The formatting handler comes after the query options (`process_options`).
            .and_then(|()| custom_check())
            .and_then(|()| blastn_group_check(&mut args, own))?,
        // The checks of the options read the argv's own format (`with_own`): with an `-outfmt`
        // whose custom delimiter the CLI rejects they fail as the CLI's.
        Commands::Tblastx(mut args) => subject_range(args.subject_loc.as_deref())
            .and_then(|()| {
                with_own(&mut args.outfmt);
                losat::web_adapter::check_tblastx_options(&args)
            })
            .map_err(|error| format!("{error:#}"))
            .and_then(|()| custom_check())?,
        Commands::Blastp(mut args) => subject_range(args.subject_loc.as_deref())
            .and_then(|()| {
                with_own(&mut args.outfmt);
                losat::web_adapter::check_blastp_options(&args)
            })
            .map_err(|error| format!("{error:#}"))
            .and_then(|()| custom_check())?,
        Commands::Tblastn(mut args) => subject_range(args.subject_loc.as_deref())
            .and_then(|()| {
                with_own(&mut args.outfmt);
                losat::web_adapter::check_tblastn_options(&args)
            })
            .map_err(|error| format!("{error:#}"))
            .and_then(|()| custom_check())?,
        // The checks of the argv's own search's formats, as for BLASTN: BLASTX resolves its
        // options for each format (the hit list size). Its ranges are read among them, where its
        // run reads them: after LOSAT's limits of the options (`blastx::check_options`); the
        // handler's check of a custom delimiter follows them.
        Commands::Blastx(args) => {
            blastx_group_check(&args, own).and_then(|()| custom_check())?;
        }
        _ => {}
    }
    Ok(settings.data_loaders)
}

/// NCBI's formatting handler's check of the argv's own `-outfmt` (a custom delimiter that a field
/// uses to separate its entries), with the CLI's message. `validate` asks it whether the CLI
/// rejects the argv, and then reports the error of the program's checks of the options with that
/// `-outfmt`, which make it where the CLI does (BLASTN's right after the query range). A run's
/// engine makes it (`run_local`) in the argv's own search, whose error is the run's (`search_all`).
fn custom_handler_check(program: Program, custom: Option<&str>) -> Result<(), String> {
    let Some(outfmt) = custom else {
        return Ok(());
    };
    losat::web_adapter::parse_formatting_string(outfmt)
        .and_then(|choice| {
            losat::web_adapter::formatting_handler_check(&choice, engine_program(program))
        })
        .map_err(|error| format!("{error:#}"))
}

/// BLASTN's checks of the options for the search of the output formats `group`, whose hit list
/// size is the one the CLI derives for them: the checks that `run_local` makes before it reads
/// a query (`process_options`: the counts of the report, `check_scoring_options`; then
/// `check_losat_limits`). `validate` runs it for the argv's own search (`own_search`).
///
/// It does not check the Karlin-Altschul tables (scores and gap costs without a row, NCBI's
/// `Error: Substitution scores ...`) nor LOSAT's limit of greedy gap costs: a CLI run meets them
/// only when some query context has an ungapped Karlin-Altschul block (`context_blocks`,
/// `any_valid_context` in `blast_engine/run.rs`), which depends on the query's letters (a query
/// of N only has none, and the same scores give rc 0 with a warning). `validate` has no query,
/// so `run` makes these checks and stops with the CLI's message and exit code.
/// Making them for a query of any letters would reject what the CLI accepts (WA-FM-01).
fn blastn_group_check(
    args: &mut losat::web_adapter::BlastnArgs,
    group: &[Output],
) -> Result<(), String> {
    use losat::web_adapter::{check_losat_limits, check_scoring_options};
    let outfmts: Vec<i32> = group.iter().map(|output| output.number as i32).collect();
    // `process_options` sets the counts that the checks read.
    let (counts, _) = losat::web_adapter::run_report_counts(
        "blastn",
        &outfmts,
        args.max_target_seqs,
        args.num_descriptions,
        args.num_alignments,
    )
    .map_err(|error| format!("{error:#}"))?;
    args.report_counts = counts;
    check_scoring_options(args)
        .and_then(|()| check_losat_limits(args))
        .map_err(|error| format!("{error:#}"))
}

/// BLASTX's checks (`blastx::check_options`) for the search of the output formats `group`, as
/// `run_local` makes them before it searches. `validate` runs it for the argv's own search
/// (`own_search`).
fn blastx_group_check(
    args: &losat::web_adapter::BlastxArgs,
    group: &[Output],
) -> Result<(), String> {
    let outfmts: Vec<&str> = group.iter().map(|output| output.outfmt.as_str()).collect();
    losat::web_adapter::check_blastx_options(args, &outfmts).map_err(|error| format!("{error:#}"))
}

/// The search of a run (`search_groups`) that writes the argv's own format (`Output::requested`).
/// Its outcome is the run's (`search_all`): the CLI runs that format alone, so the run completes
/// when that search does, and fails with its error and the warnings that it wrote before it when it
/// fails (line WQA's sweep, D3 review M1, 2026-11-19). The other search writes only formats that the
/// argv does not ask for (the Web writes outfmt 0, 6 and 7 for every run, line NP): its failure
/// fails those formats alone, as the CLI's run of each of them alone fails.
fn own_search(groups: &[Vec<Output>]) -> usize {
    groups
        .iter()
        .position(|group| group.iter().any(|output| output.requested))
        .expect("the argv's format is one of the run's outputs (`outputs`)")
}

/// Records the byte range of every HSP row (stream 6) and section (stream 0), and of
/// the subject heading (stream 0) that each section belongs to.
struct RangeRecorder {
    streams: Vec<u32>,
    positions: Vec<Arc<AtomicU64>>,
    open: Vec<Option<(HspIndex, u64)>>,
    out6: Vec<Option<(u64, u64)>>,
    out0: Vec<Option<(u64, u64)>>,
    out0_subject: Vec<Option<(u64, u64)>>,
    /// The start of the heading being written, and the last complete heading.
    subject_start: Option<u64>,
    subject: Option<(u64, u64)>,
    /// The first observer event that does not fit the formatter contract; the run fails.
    fault: Option<String>,
}

impl RangeRecorder {
    fn position(&self, format: usize) -> u64 {
        self.positions[format].load(Ordering::SeqCst)
    }

    /// Forgets the events of the format `format`, which failed alone (`search`): its stream is
    /// dropped, so no HSP record points into it (stream 6: `out6`; stream 0: `out0` and
    /// `out0_subject`), and an HSP that its report had begun when it stopped is not a fault.
    fn forget(&mut self, format: usize) {
        self.open[format] = None;
        match self.streams[format] {
            6 => self.out6.clear(),
            0 => {
                self.out0.clear();
                self.out0_subject.clear();
                self.subject_start = None;
                self.subject = None;
            }
            _ => {}
        }
    }
}

/// Stores the range of `hsp`; a second range for the same HSP is a fault.
fn store(
    ranges: &mut Vec<Option<(u64, u64)>>,
    hsp: HspIndex,
    range: (u64, u64),
) -> Result<(), String> {
    if ranges.len() <= hsp {
        ranges.resize(hsp + 1, None);
    }
    match ranges[hsp].replace(range) {
        None => Ok(()),
        Some(_) => Err(format!("HSP {hsp} was reported twice in one format")),
    }
}

impl FormatObserver for RangeRecorder {
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex) {
        let start = self.position(format);
        if let Some((open, _)) = self.open[format].replace((hsp, start)) {
            self.fault
                .get_or_insert(format!("HSP {hsp} started before HSP {open} ended"));
        }
    }

    fn hsp_end(&mut self, format: usize, hsp: HspIndex) {
        let Some((begun, start)) = self.open[format].take() else {
            self.fault
                .get_or_insert(format!("HSP {hsp} ended without a start"));
            return;
        };
        if begun != hsp {
            self.fault
                .get_or_insert(format!("HSP {begun} started but HSP {hsp} ended"));
            return;
        }
        let range = (start, self.position(format));
        let stored = match self.streams[format] {
            6 => store(&mut self.out6, hsp, range),
            0 => store(&mut self.out0, hsp, range).and_then(|()| match self.subject {
                Some(subject) => store(&mut self.out0_subject, hsp, subject),
                None => Ok(()),
            }),
            _ => Ok(()),
        };
        if let Err(fault) = stored {
            self.fault.get_or_insert(fault);
        }
    }

    fn subject_begin(&mut self, format: usize, _first_hsp: HspIndex) {
        if self.streams[format] == 0 {
            self.subject_start = Some(self.position(format));
        }
    }

    fn subject_end(&mut self, format: usize, _first_hsp: HspIndex) {
        if self.streams[format] == 0 {
            if let Some(start) = self.subject_start.take() {
                self.subject = Some((start, self.position(format)));
            }
        }
    }
}

/// The registered index of each query that the engine searches: NCBI's batch reader skips
/// a query whose `-query_loc` interval starts more than one letter past its end
/// (`seq_range::cut_queries`), so the engine's `q_idx` counts the searched queries; the
/// HSP records give the registered index (docs/web/abi_v2.md §8, port plan Q8). `None`
/// without `-query_loc` (every query is searched); a range that does not parse fails the
/// run before any hit.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-155
fn searched_queries(
    query_loc: Option<&std::ffi::OsStr>,
    records: &[FastaRecord],
) -> Option<Vec<usize>> {
    use losat::web_adapter::{parse_optional_range, record_interval, RangeRole, RecordInterval};
    let range = parse_optional_range(query_loc, RangeRole::Query).ok()??;
    Some(
        records
            .iter()
            .enumerate()
            .filter(|(_, record)| {
                record_interval(&range, record.seq().len()) != RecordInterval::PastEnd
            })
            .map(|(index, _)| index)
            .collect(),
    )
}

/// The options that set the hit list size: `-max_target_seqs`, `-num_descriptions` and
/// `-num_alignments` of the program's arguments.
fn count_options(command: &Commands) -> [Option<usize>; 3] {
    match command {
        Commands::Blastp(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Tblastn(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Blastn(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        Commands::Tblastx(a) => [a.max_target_seqs, a.num_descriptions, a.num_alignments],
        // A positive `i32` (the parser's `positive`).
        Commands::Blastx(a) => [
            a.max_target_seqs.map(|n| n as usize),
            a.num_descriptions,
            a.num_alignments,
        ],
        // Rejected by `Program::parse`.
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => [None; 3],
    }
}

/// `-num_threads` of the program's arguments.
fn num_threads(command: &Commands) -> usize {
    match command {
        Commands::Blastp(a) => a.num_threads,
        Commands::Tblastn(a) => a.num_threads,
        Commands::Blastn(a) => a.num_threads,
        Commands::Tblastx(a) => a.num_threads,
        // A positive `i32` (the parser's `positive`).
        Commands::Blastx(a) => a.num_threads as usize,
        // Rejected by `Program::parse`.
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => 1,
    }
}

/// The hit list size that the CLI derives for the output format `format` from the options
/// of `command` (`report_counts`, NCBI `blast_args.cpp:2894-2973`).
pub fn hit_list_size(command: &Commands, format: u32) -> usize {
    let [max_target_seqs, descriptions, alignments] = count_options(command);
    losat::web_adapter::report_counts(format as i32, max_target_seqs, descriptions, alignments)
        .hitlist_size
}

/// The searches of a run: the output formats grouped by hit list size, the group of the tabular
/// formats first. One group is one `run_local` call, and every format of it is written from
/// that search, as the CLI writes it alone: the CLI's search depends on the hit list size,
/// which `-num_descriptions` and `-num_alignments` (without `-max_target_seqs`) can make
/// differ between the pairwise format 0 and the tabular formats. Without those options the
/// sizes are equal and the run is one search of every format. The output of the argv's `-outfmt`
/// has the hit list size of its number's formats (that of the tabular formats above 4, that of
/// format 0 for 1 to 4), so it joins their group: it adds no search.
///
/// In each group the argv's own format (`Output::requested`) comes first, then the tabular
/// formats: the first format of a search writes the warnings that its formats share (the
/// formatting handler's, the queries'), and the argv's are the ones that the CLI writes for the
/// argv (the handler warns of `-num_descriptions`, `-line_length` and `-sorthits` for formats
/// above 4 only and of `-sorthsps` for formats other than 0, blast_args.cpp:2947-2989), whatever
/// the outcome of the run (line WQA's sweep, D3 review r2 M1, 2026-11-19: with the tabular
/// formats first, a failed run whose argv asks for format 0 to 4 gave their lines). `search_all`
/// adds the lines that the other formats' handler writes where the run completes.
pub fn search_groups(formats: &[Output], command: &Commands) -> Vec<Vec<Output>> {
    let [_, descriptions, alignments] = count_options(command);
    let mut ordered: Vec<Output> = formats.to_vec();
    ordered.sort_by_key(|format| format.number <= 4);
    let mut groups: Vec<(usize, Vec<Output>)> = Vec::new();
    for format in ordered {
        let size = match descriptions.is_none() && alignments.is_none() {
            true => 0,
            false => hit_list_size(command, format.number),
        };
        match groups.iter_mut().find(|(group, _)| *group == size) {
            Some((_, group)) => group.push(format),
            None => groups.push((size, vec![format])),
        }
    }
    groups
        .into_iter()
        .map(|(_, mut group)| {
            group.sort_by_key(|format| !format.requested);
            group
        })
        .collect()
}

/// The formatting handler's warnings that the CLI writes for the output format `number` with
/// the options of `command` (`ReportCounts::warnings_with`, which each program's option handler
/// writes for the first format of a search): `search_all` adds those of the formats after the
/// first where the run completes.
fn handler_warnings(program: Program, command: &Commands, number: u32) -> Vec<u8> {
    let [max_target_seqs, descriptions, alignments] = count_options(command);
    let counts =
        losat::web_adapter::report_counts(number as i32, max_target_seqs, descriptions, alignments);
    let warn = |formatting| {
        counts.warnings_with(
            program.name(),
            number as i32,
            descriptions.is_some(),
            formatting,
        )
    };
    match command {
        Commands::Blastp(a) => warn(&a.formatting),
        Commands::Tblastn(a) => warn(&a.formatting),
        Commands::Blastn(a) => warn(&a.formatting),
        Commands::Tblastx(a) => warn(&a.formatting),
        Commands::Blastx(a) => warn(&a.formatting),
        // Rejected by `Program::parse`.
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => Vec::new(),
    }
}

/// A search that failed.
struct Failed {
    /// The error of the run, when this failure fails it.
    message: String,
    /// What the CLI writes to standard error for the error: NCBI's message (a
    /// `NativeError`), or `Error: ` and the error, as `main` returns it.
    cli_stderr: Vec<u8>,
    /// NCBI's message (a `NativeError`), which the CLI also writes to the file of `-logfile` when
    /// the run's messages go there (`CallerStart::post`); `None` for LOSAT's own errors, which the
    /// CLI writes to the standard error whatever `-logfile` is (`Failed::source`).
    native: Option<Vec<u8>>,
    /// An error of the adapter itself (the output ranges, the archive of outfmt 13 or 14), which
    /// no CLI run has: it ends the run, whichever search it comes from (`search_all`).
    fault: bool,
    /// The writers of the argv's own outputs (`Output::requested`) of the search that failed, with what
    /// they hold: the bytes that the CLI's run of the argv leaves in its output before its error, which
    /// the run sends (`search_all`; line WQA's sweep, D15, 2026-11-20). None for a fault.
    kept: Vec<StreamWriter>,
}

impl Failed {
    /// `error` is the engine's `anyhow::Error`, which this crate does not name.
    fn new<E>(error: &E) -> Self
    where
        E: std::ops::Deref<Target = dyn std::error::Error + Send + Sync + 'static>
            + std::fmt::Display
            + std::fmt::Debug,
    {
        let native = error
            .downcast_ref::<losat::web_adapter::NativeError>()
            .map(|native| native.message.clone());
        let cli_stderr = match &native {
            Some(native) => native.clone(),
            None => format!("Error: {error:?}\n").into_bytes(),
        };
        Self {
            message: format!("{error:#}"),
            cli_stderr,
            native,
            fault: false,
            kept: Vec::new(),
        }
    }

    /// Whose error this is, which decides where the CLI writes it (`CallerSink::error`).
    fn source(&self) -> ErrorSource {
        match self.native {
            Some(_) => ErrorSource::Ncbi,
            None => ErrorSource::Losat,
        }
    }
}

/// The name that the side files of outfmt 13 and 14 follow, as the CLI's `-out` (no directory and
/// no extension, so that the side files have the names that the master file lists: `results_1.json`
/// for 13 too, where an `-out` ending in `.json` lists other names than its side files have,
/// `report::blast_xml2::files`).
pub const SIDE_FILE_BASE: &str = "results";

/// The name of the master file in the archive of outfmt 13 or 14.
pub fn master_name(number: u32) -> String {
    let extension = if number == 13 { "json" } else { "xml" };
    format!("{SIDE_FILE_BASE}.{extension}")
}

/// The output of outfmt 13 or 14 (stream 4): a ZIP archive (`zip::ZipStream`) of the side files,
/// added as the engine hands each over, and of the master file, added after the run.
struct ArchivedFiles(std::sync::Mutex<Option<crate::zip::ZipStream<StreamWriter>>>);

impl ArchivedFiles {
    fn take(&self) -> crate::zip::ZipStream<StreamWriter> {
        self.0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .take()
            .expect("the archive is taken once, after the run")
    }
}

impl losat::web_adapter::FileReceiver for ArchivedFiles {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
        match self
            .0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .as_mut()
        {
            Some(zip) => zip.add(name, bytes),
            None => Err(std::io::Error::other("the archive was closed")),
        }
    }
}

/// What one search wrote.
struct Searched {
    /// The streams of the formats that were written (not those of `refused`).
    writers: Vec<StreamWriter>,
    recorder: RangeRecorder,
    hits: Option<std::sync::Arc<Vec<PairwiseHit>>>,
    ordinals: Option<Vec<usize>>,
    /// The formats that failed alone (`FormatRefusal::FailsAlone`): each one's stream and
    /// refusal, in the order of the formats. Their streams stay empty.
    refused: Vec<(u32, Failed)>,
    /// By format: what the CLI's standard error of a format after the first has of its queries'
    /// reports where the first format, a structured one, does not write the queries' warnings
    /// (`ReportOutputs::later_warnings`): each format's own text, up to where its report stopped
    /// if it failed alone, for stream 3 where the run completes (`search_all`, which keeps each
    /// line once).
    later: Vec<Vec<u8>>,
}

impl Searched {
    /// Whether the format of `stream` failed alone.
    fn refused(&self, stream: u32) -> bool {
        self.refused.iter().any(|(refused, _)| *refused == stream)
    }
}

/// One `run_local` call that writes the output formats `formats` (a group of
/// `search_groups`), with the warnings of the CLI to `diagnostics`. A format other than the
/// argv's own (`Output::requested`) whose report fails, at a check before its rows or in its
/// middle (`FormatRefusal`), fails alone (`FormatRefusal::FailsAlone`, `Searched::refused`): the
/// call writes the others, as the CLI's run of each format alone writes it; an error of the argv's
/// own format fails the call, as it fails the CLI. The bytes that a format that failed alone had
/// written stay with its writer, which is dropped, and its observer events with the recorder,
/// which forgets them (`RangeRecorder::forget`); `search_all` tells the host that the stream
/// failed (`emit::discard`), and the host drops what it received of it where a chunk reached it
/// already (a report that failed after 1 MiB, BLASTX's later query batch).
fn search(
    mut command: Commands,
    formats: &[Output],
    inputs: RunInputs,
    settings: &ApplicationSettings,
    diagnostics: &mut (dyn std::io::Write + Send),
) -> Result<Searched, Failed> {
    // An imported strategy's warnings (`RecoverSearchStrategy`, `x_IssueWarningsForIgnoredOptions`)
    // come first, before the program's own, as the CLI posts them at its start.
    let warnings = losat::web_adapter::take_strategy_warnings(&mut command);
    if let Err(error) = diagnostics.write_all(&warnings) {
        let message = format!("the diagnostics could not be written: {error}");
        return Err(Failed {
            cli_stderr: message.clone().into_bytes(),
            message,
            native: None,
            fault: true,
            kept: Vec::new(),
        });
    }
    let RunInputs {
        query,
        query_end,
        subject,
    } = inputs;
    let mut writers: Vec<StreamWriter> = formats
        .iter()
        .map(|format| StreamWriter::new(format.stream))
        .collect();
    let mut recorder = RangeRecorder {
        streams: formats.iter().map(|format| format.stream).collect(),
        positions: writers.iter().map(StreamWriter::position).collect(),
        open: vec![None; formats.len()],
        out6: Vec::new(),
        out0: Vec::new(),
        out0_subject: Vec::new(),
        subject_start: None,
        subject: None,
        fault: None,
    };
    let mut collected: Option<std::sync::Arc<Vec<PairwiseHit>>> = None;
    let mut later = vec![Vec::new(); formats.len()];
    // Outfmt 13 and 14: the master file goes to memory and the side files to a ZIP archive on the
    // format's stream, to which the master file is added after the run (`ArchivedFiles`).
    let mut archives: Vec<Option<(Vec<u8>, Arc<ArchivedFiles>)>> = formats
        .iter()
        .zip(writers.iter_mut())
        .map(|(format, writer)| {
            matches!(format.number, 13 | 14).then(|| {
                let stream = std::mem::replace(writer, StreamWriter::new(format.stream));
                let zip = crate::zip::ZipStream::new(stream);
                (
                    Vec::new(),
                    Arc::new(ArchivedFiles(std::sync::Mutex::new(Some(zip)))),
                )
            })
        })
        .collect();
    let ordinals;
    let outcome = {
        let mut on_hits = |hits| collected = Some(hits);
        let mut outputs = ReportOutputs {
            // The run's formats besides the argv's fail alone, each as the CLI's run of it alone
            // (`FormatRefusal::FailsAlone`): no search answers for another's formats. (A subject
            // whose lookup throws where the prefetch is off, which `other_outfmts` is for, needs
            // `PRE_FETCH_SEQS_LIMIT`, which LOSAT Web does not set.)
            other_outfmts: &[],
            html_posts: None,
            formats: formats
                .iter()
                .zip(writers.iter_mut())
                .zip(archives.iter_mut())
                .map(|((format, writer), archive)| FormatOutput {
                    outfmt: format.outfmt.as_bytes(),
                    refusal: if format.requested {
                        FormatRefusal::FailsCall
                    } else {
                        FormatRefusal::FailsAlone
                    },
                    sink: match archive {
                        Some((master, files)) => OutputSink::Files(
                            master,
                            losat::web_adapter::SideFileSink::new(SIDE_FILE_BASE, files.clone()),
                        ),
                        None => OutputSink::Writer(writer),
                    },
                })
                .collect(),
            diagnostics,
            later_warnings: Some(&mut later),
            hits: Some(&mut on_hits),
            observer: Some(&mut recorder),
        };
        // A search of `-db` (line WQA D12a): the program's database entry, as the library's
        // `api::run` runs it, with the run's application settings, which find the database (the
        // run's working directory, then `BLASTDB` of a registered `.ncbirc` or `-conffile`, as on
        // the command line). The engine opens its files by name, from the files that the host
        // gives it (docs/web/abi_v2.md §3); the subject records are none.
        let database = database_of(&command).is_some();
        let result = match command {
            Commands::Blastp(args) if database => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_database_blastp(args, query, query_end, settings.clone(), &mut outputs)
            }
            Commands::Tblastn(args) if database => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                losat::web_adapter::refuse_tblastn_pssm_query(&args).and_then(|()| {
                    run_local_database_tblastn(
                        args,
                        query,
                        query_end,
                        settings.clone(),
                        &mut outputs,
                    )
                })
            }
            Commands::Blastn(args) if database => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_database_blastn(args, query, query_end, settings.clone(), &mut outputs)
            }
            Commands::Tblastx(args) if database => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_database_tblastx(args, query, query_end, settings.clone(), &mut outputs)
            }
            Commands::Blastx(args) if database => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_database_blastx(args, query, query_end, settings.clone(), &mut outputs)
            }
            Commands::Blastp(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_blastp(args, query, query_end, subject, "", "", &mut outputs)
            }
            // The host has no PSSM files: `-in_pssm` is refused as by `validate`
            // (`check_tblastn_options`), before the search.
            Commands::Tblastn(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                losat::web_adapter::refuse_tblastn_pssm_query(&args)
                    .and_then(|()| run_local_tblastn(args, query, query_end, subject, &mut outputs))
            }
            Commands::Blastn(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_blastn(args, query, query_end, subject, &mut outputs)
            }
            Commands::Tblastx(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_tblastx(args, query, query_end, subject, &mut outputs)
            }
            Commands::Blastx(args) => {
                ordinals = searched_queries(args.query_loc.as_deref(), query);
                run_local_blastx(args, query, query_end, subject, &mut outputs)
            }
            Commands::Psiblast(_)
            | Commands::Rpsblast(_)
            | Commands::Rpstblastn(_)
            | Commands::Makedb(_)
            | Commands::Rpsbproc(_) => {
                unreachable!("rejected by Program::parse")
            }
        };
        // The refusals of the formats that failed alone, by format.
        let refusals: Vec<_> = outputs
            .formats
            .iter_mut()
            .map(|format| format.refusal.take_refused())
            .collect();
        (result, refusals)
    };
    let (outcome, refusals) = outcome;
    // The archives: the master file after the side files and the central directory, on the
    // format's stream again. A run that failed leaves the side files written (the host discards
    // them).
    for ((format, writer), archive) in formats.iter().zip(writers.iter_mut()).zip(&mut archives) {
        let Some((master, files)) = archive.take() else {
            continue;
        };
        let zip = files.take();
        *writer = if outcome.is_ok() {
            let mut zip = zip;
            match zip
                .add(master_name(format.number).as_bytes(), &master)
                .and_then(|()| zip.finish())
            {
                Ok(stream) => stream,
                Err(error) => {
                    return Err(Failed {
                        cli_stderr: format!("{error}\n").into_bytes(),
                        message: format!("the archive of outfmt {}: {error}", format.number),
                        native: None,
                        fault: true,
                        kept: Vec::new(),
                    })
                }
            }
        } else {
            // The side files written before the error and the master file as the run holds it (the CLI
            // leaves its `-out` file empty when it fails before the master file is written), with the
            // archive's central directory: a run that failed keeps them (line WQA's sweep, D15).
            let mut zip = zip;
            if let Err(error) = zip.add(master_name(format.number).as_bytes(), &master) {
                return Err(Failed {
                    cli_stderr: format!("{error}\n").into_bytes(),
                    message: format!("the archive of outfmt {}: {error}", format.number),
                    native: None,
                    fault: true,
                    kept: Vec::new(),
                });
            }
            match zip.finish() {
                Ok(stream) => stream,
                Err(error) => {
                    return Err(Failed {
                        cli_stderr: format!("{error}\n").into_bytes(),
                        message: format!("the archive of outfmt {}: {error}", format.number),
                        native: None,
                        fault: true,
                        kept: Vec::new(),
                    })
                }
            }
        };
    }
    if let Err(error) = outcome {
        // The writers of the argv's own outputs keep what they hold: the CLI's run of the argv leaves it
        // in its output before the error (`Failed::kept`).
        let mut failed = Failed::new(&error);
        failed.kept = formats
            .iter()
            .zip(writers)
            .filter(|(format, _)| format.requested)
            .map(|(_, writer)| writer)
            .collect();
        return Err(failed);
    }
    // The formats that failed alone: their writers are dropped with what they hold, and the
    // recorder forgets their events (`search_all` tells the host that each failed).
    let mut written = Vec::with_capacity(writers.len());
    let mut refused = Vec::new();
    for (index, ((format, writer), refusal)) in
        formats.iter().zip(writers).zip(refusals).enumerate()
    {
        match refusal {
            None => written.push(writer),
            Some(error) => {
                recorder.forget(index);
                refused.push((format.stream, Failed::new(&error)));
            }
        }
    }
    if let Some((hsp, _)) = recorder.open.iter().flatten().next() {
        recorder
            .fault
            .get_or_insert(format!("HSP {hsp} started but did not end"));
    }
    if let Some(fault) = recorder.fault.take() {
        let message = format!("internal error in the output ranges: {fault}");
        return Err(Failed {
            cli_stderr: message.clone().into_bytes(),
            message,
            native: None,
            fault: true,
            kept: Vec::new(),
        });
    }
    Ok(Searched {
        writers: written,
        recorder,
        hits: collected,
        ordinals,
        refused,
        later,
    })
}

/// The HSPs of one (query, subject) pair, in the list's order.
type Pairs = std::collections::HashMap<(u32, u32), Vec<usize>>;

fn pairs(hits: &[PairwiseHit], only: impl Fn(usize) -> bool) -> Pairs {
    let mut pairs = Pairs::new();
    for (index, record) in hits.iter().enumerate() {
        if only(index) {
            pairs
                .entry((record.hit.q_idx, record.hit.s_idx))
                .or_default()
                .push(index);
        }
    }
    pairs
}

/// Whether two searches found the same HSP: every field of the `PairwiseHit`, the numbers
/// of `f64` by their bits.
fn same_hsp(a: &PairwiseHit, b: &PairwiseHit) -> bool {
    let (x, y) = (&a.hit, &b.hit);
    x.identity.to_bits() == y.identity.to_bits()
        && x.length == y.length
        && x.mismatch == y.mismatch
        && x.gapopen == y.gapopen
        && (x.q_start, x.q_end, x.s_start, x.s_end) == (y.q_start, y.q_end, y.s_start, y.s_end)
        && x.e_value.to_bits() == y.e_value.to_bits()
        && x.bit_score.to_bits() == y.bit_score.to_bits()
        && x.num_ident == y.num_ident
        && x.query_frame == y.query_frame
        && x.query_length == y.query_length
        && (x.q_idx, x.s_idx, x.raw_score) == (y.q_idx, y.s_idx, y.raw_score)
        && (x.sort_query_offset, x.sort_query_end) == (y.sort_query_offset, y.sort_query_end)
        && (x.sort_subject_offset, x.sort_subject_end)
            == (y.sort_subject_offset, y.sort_subject_end)
        && x.has_sort_offsets == y.has_sort_offsets
        && x.gap_info == y.gap_info
        && x.num_positives == y.num_positives
        && a.query_seq == b.query_seq
        && a.subject_seq == b.subject_seq
        && a.query_frame == b.query_frame
        && a.subject_frame == b.subject_frame
        && a.positives == b.positives
        && a.gaps == b.gaps
        && a.subject_length == b.subject_length
        && a.subject_title == b.subject_title
        && a.comp_adjust_method == b.comp_adjust_method
        && a.sum_n == b.sum_n
}

/// Gives the HSPs of the tabular search `first` the sections of format 0 from the other
/// search, `second`, for the pairs (query, subject) whose HSPs in `second`'s format 0 are
/// exactly those of `first`: the same count and, in order, every field of the HSP equal
/// (`same_hsp`). A format 0 hit list
/// of another size can keep other subjects or other HSPs of a subject (the preliminary
/// list evicts a subject before its alignment is extended), and such a pair stays without
/// a section. Records, rows and sections that this leaves unlinked are still correct.
fn link_format0(
    first: &[PairwiseHit],
    first_recorder: &mut RangeRecorder,
    second: &[PairwiseHit],
    second_recorder: &RangeRecorder,
) {
    let headed = pairs(second, |index| {
        second_recorder.out0.get(index).copied().flatten().is_some()
    });
    let ours = pairs(first, |_| true);
    let same = same_hsp;
    for (pair, theirs) in headed {
        let Some(mine) = ours.get(&pair) else {
            continue;
        };
        if mine.len() != theirs.len()
            || !mine
                .iter()
                .zip(&theirs)
                .all(|(&a, &b)| same(&first[a], &second[b]))
        {
            continue;
        }
        for (&a, &b) in mine.iter().zip(&theirs) {
            for (target, source) in [
                (&mut first_recorder.out0, &second_recorder.out0),
                (
                    &mut first_recorder.out0_subject,
                    &second_recorder.out0_subject,
                ),
            ] {
                if target.len() <= a {
                    target.resize(a + 1, None);
                }
                target[a] = source.get(b).copied().flatten();
            }
        }
    }
}

/// Whether `validate` reads the WindowMasker options where the command line reads them: after
/// `-dust`, before the query range and the other options (`resolve_window_masker`: the check of
/// the file of `-window_masker_db`, with the CLI's message for a format other than the two binary
/// ones, then "at most one of" the options, then the file of `-window_masker_taxid`). Not when
/// `-window_masker_db` names a file that is not registered: a host that validates an argv without
/// its files (the application's control reactor) leaves the file to `run`, where a name that is
/// not registered is the CLI's file that cannot be opened. A taxid whose file is not registered
/// is the CLI's taxid without a file, whose error comes with the search (`batch_error`). An argv
/// without the options is validated as before named input files.
fn validates_window_masker(args: &losat::web_adapter::BlastnArgs) -> bool {
    match args.window_masker_db.as_deref() {
        Some(name) => name.to_str().is_some_and(crate::files::registered),
        None => args.window_masker_taxid.is_some(),
    }
}

/// The counts source of BLASTN's `-window_masker_db` and `-window_masker_taxid`: the named input
/// files registered when the run starts (`files::source`), which the engine opens by the option's
/// value, or by the taxid's name in `files::WINDOW_MASKER_PATH`, as the command line opens the
/// file of that name. `None` without the options, so a run without them gives the engine what it
/// gave before named input files.
fn window_masker_counts(command: &Commands) -> Option<Arc<dyn CountsSource>> {
    match command {
        Commands::Blastn(args)
            if args.window_masker_db.is_some() || args.window_masker_taxid.is_some() =>
        {
            Some(crate::files::source())
        }
        _ => None,
    }
}

/// Gives a search's command what the run's start and host provide: the counts source
/// (`window_masker_counts`) and, for BLASTN's `-filtering_db`, the repeat search that looks its
/// database up with the run's settings.
fn give_run_resources(
    command: &mut Commands,
    counts: Option<&Arc<dyn CountsSource>>,
    settings: &ApplicationSettings,
) {
    if let (Commands::Blastn(args), Some(counts)) = (&mut *command, counts) {
        give_blastn_counts(args, counts, settings);
    }
    // The repeat search of `-filtering_db`, for a `-db` and a `-subject` search alike (a `-db`
    // search's entry makes the same call).
    if let Commands::Blastn(args) = command {
        losat::web_adapter::resolve_blastn_repeat_filter(args, settings);
    }
}

/// Gives BLASTN's arguments the counts source, and with `-window_masker_taxid` the registry of
/// `-window_masker_taxid` of the run's start (`ApplicationSettings::window_masker_registry`), as
/// the command line gives its own: the directory `files::WINDOW_MASKER_PATH` (`START_VARIABLES`),
/// and the posts of a registered `.ncbirc`.
fn give_blastn_counts(
    args: &mut losat::web_adapter::BlastnArgs,
    counts: &Arc<dyn CountsSource>,
    settings: &ApplicationSettings,
) {
    args.window_masker_counts = Some(counts.clone());
    if args.window_masker_taxid.is_some() {
        args.window_masker_registry = Some(settings.window_masker_registry.clone());
    }
}

/// Runs the searches of the run and emits its streams. On failure the output streams stop: the
/// output of the argv's own `-outfmt` sends what the CLI's run of the argv leaves in its output
/// before the error, and the other outputs are dropped (`keep_failed_outputs`; line WQA's sweep,
/// D15, 2026-11-20); the run still sends stream 3 (what the CLI writes to its standard error before
/// the error that ends the run, without that error, which is the run's) and stream 8 (the file of
/// `-logfile`), which the host keeps with the failed run (docs/web/abi_v2.md §5).
/// The query's records end with `query_end`, how its reading ended (`store::Registered::end`): a
/// reader's error ends the run in the query batch that reads the record after them, after the
/// batches before, as on the command line (line WQA's D10, 2026-11-19).
///
/// Without `-num_descriptions` and `-num_alignments` this is one search that writes every
/// format. With them it is one search per distinct hit list size of the formats
/// (`search_groups`), so that each format's stream is the CLI's output for that format
/// alone. The output of the argv's `-outfmt` (stream 4) is written by the search of its number's
/// formats, first, as the CLI writes that `-outfmt`. The HSP records follow the tabular search;
/// a record gets the section of format 0 only where the search of format 0 has exactly its HSPs
/// (`link_format0`) (docs/web/abi_v2.md §5, §7, §8).
///
/// The run's outcome is that of the search of the argv's own format (`own_search`), as the CLI
/// runs the argv: the run completes when that search completes, and fails with its error when it
/// fails, with stream 3 the CLI's standard error for the argv up to that error (that search's
/// warnings). The failure of the other search of a run of two searches (`-num_alignments 0` for the
/// tabular formats) fails its formats alone, as the CLI fails those formats alone: their streams stay empty, what the host received of them is discarded, and
/// without the tabular search there are no records. A fault of the adapter itself (`Failed::fault`)
/// fails the run in any search.
///
/// Stream 3 of a run that completes is the CLI's standard error for the argv (the warnings of the
/// argv's own search), followed by what the CLI's standard error for the run's other formats adds
/// to it (`search_all`): the formatting handler's lines of the formats besides the argv's, the
/// lines of the other search's diagnostics that stream 3 does not hold yet, then the error of a
/// search that failed, or the refusals of the formats that failed alone, so that each text comes
/// once (`Diagnostics`). Where the CLI writes the run's messages to the file of `-logfile`
/// (TBLASTX), the file has what the CLI's file has for the argv, and stream 3 the rest.
pub fn run(
    words: &[&str],
    query: &[FastaRecord],
    query_end: &QueryEnd,
    subject: &[FastaRecord],
) -> Result<(), String> {
    let registry = crate::files::registry_files();
    let files = Arc::new(RunFileArchive::new(search_run_files(words, &registry)));
    let run_files = RunFiles::Receiver(files.clone());
    let mut diagnostics = StreamWriter::new(STREAM_DIAGNOSTICS);
    let (result, started) = match start(words, &registry, &run_files, &mut diagnostics) {
        Ok(mut started) => (
            run_started(
                words,
                RunInputs {
                    query,
                    query_end,
                    subject,
                },
                &mut started,
                &run_files,
                &mut diagnostics,
            ),
            Some(started),
        ),
        Err((message, started)) => (Err(message), started.map(|started| *started)),
    };
    // Stream 3 of a run that failed: what the command line wrote to its standard error before the
    // error that ends it (its warnings: NCBI's `Subject sequence contains no data` before `The
    // average subject length is too short`), which the host keeps with the failed run as it keeps
    // stream 8 (docs/web/abi_v2.md §5; line WQA's sweep, D3, 2026-11-19). The error itself is the
    // run's error. A run that succeeded sent it already.
    diagnostics.finish();
    // The run's files (the file of `-logfile`), after the run, whether it failed or not.
    let handed = started.map_or(Ok(()), CallerStart::finish);
    // The archive closes here, also where the start's files failed: a file handed over later is
    // refused (`RunFileArchive::receive`).
    let kept = files.close();
    let emitted = handed.and_then(|()| RunFileArchive::emit(kept));
    // The run's import of its search strategy ends with it (`imported`).
    crate::imported::end_run();
    result.and(emitted.map_err(|error| format!("the files of the run: {error}")))
}

/// Entries queued in [`RunFileArchive`] before the run emits them on `STREAM_RUN_FILES`.
type ArchivedEntries = Vec<(Vec<u8>, Vec<u8>)>;

/// The files of a run besides its outputs (the file of `-logfile`), which the engine hands over
/// when the run ends (`RunFiles::Receiver`), emitted as a ZIP archive on `STREAM_RUN_FILES` after
/// the run, whether it succeeded or failed (docs/web/abi_v2.md §5), and only when there are any.
/// The archive closes once, after the run (`close`): a file handed over after it (`None`) is
/// refused, so that a file that the engine hands over late is its error (LOSAT's internal error of
/// `-export_search_strategy`, `strategy::end_run`) instead of a file lost without one.
struct RunFileArchive {
    files: std::sync::Mutex<Option<ArchivedEntries>>,
    run_names: SearchRunFiles,
}

impl Default for RunFileArchive {
    fn default() -> Self {
        Self::new(SearchRunFiles::default())
    }
}

impl RunFileArchive {
    fn new(run_names: SearchRunFiles) -> Self {
        Self {
            files: std::sync::Mutex::new(Some(Vec::new())),
            run_names,
        }
    }

    /// The files handed over so far; any later one is refused.
    fn close(&self) -> ArchivedEntries {
        self.files
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .take()
            .unwrap_or_default()
    }

    /// Emits the archive of the files `files` (`close`; stored entries, as stream 4's,
    /// `zip::ZipStream`). `Err` when the archive cannot hold a file (a name of 65,536 bytes or
    /// more).
    fn emit(files: ArchivedEntries) -> std::io::Result<()> {
        if files.is_empty() {
            return Ok(());
        }
        let mut zip = crate::zip::ZipStream::new(StreamWriter::new(STREAM_RUN_FILES));
        for (name, bytes) in &files {
            zip.add(name, bytes)?;
        }
        zip.finish()?.finish();
        Ok(())
    }
}

impl losat::web_adapter::FileReceiver for RunFileArchive {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
        match &mut *self
            .files
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
        {
            Some(files) => {
                files.push((name.to_vec(), bytes.to_vec()));
                Ok(())
            }
            None => Err(std::io::Error::other(format!(
                "{} came after the run's files were emitted",
                String::from_utf8_lossy(name)
            ))),
        }
    }

    fn registered_files(&self) -> &[Vec<u8>] {
        &self.run_names.registered
    }

    fn read_files(&self) -> &[Vec<u8>] {
        &self.run_names.read
    }
}

/// The inputs of a run (`run`): the query records, how the reading of the query ended after them,
/// and the subject records.
#[derive(Clone, Copy)]
struct RunInputs<'a> {
    query: &'a [FastaRecord],
    query_end: &'a QueryEnd,
    subject: &'a [FastaRecord],
}

/// `run` after the application start `started` (`started_command`: the argv, the settings that
/// LOSAT refuses, a dry run's end), then the searches. The messages that the command line writes
/// to its standard error for the argv go where the run's messages go (`CallerStart::sink`: the
/// file of `-logfile` in TBLASTX, else stream 3), as they come, or when the searches end
/// (`search_all`), and the error that ends the run with an NCBI message is posted there after them
/// (`CallerStart::post`), as on the command line. What the other formats' standard error adds to
/// them goes to stream 3 (`search_all`, `CallerSink::add`), which no file of the argv's CLI run
/// holds (line WQA's sweep, D3 review r2, 2026-11-19: TBLASTX's file of `-logfile` had it), and
/// nowhere where the run's messages go nowhere. The host is told first where the run's messages
/// go (`emit::messages`, stream 10); the error of each output that the run could not write is
/// where the sink wrote it (`CallerSink::error`, stream 9).
fn run_started(
    words: &[&str],
    inputs: RunInputs,
    started: &mut CallerStart,
    run_files: &RunFiles,
    diagnostics: &mut StreamWriter,
) -> Result<(), String> {
    // Where the messages written to the sink end up, which the host is told first (stream 10; line
    // WQA's sweep, D7, 2026-11-20).
    crate::emit::messages(crate::emit::place(started.messages_place()));
    let (parsed, strategy) = match started_command(words, started, run_files, true)? {
        Next::Search(parsed, strategy) => (*parsed, strategy),
        // The command line's text on the stream of its standard output (line WQA's D14,
        // 2026-11-20); no search.
        Next::Prints(text) => {
            emit(STREAM_PRINTED, text.as_bytes());
            diagnostics.finish();
            return Ok(());
        }
        Next::Ends => {
            diagnostics.finish();
            return Ok(());
        }
    };
    let settings = started.settings()?.clone();
    let mut posted = None;
    let searched = {
        let mut sink = started.sink(diagnostics);
        search_all(
            parsed,
            strategy,
            words,
            inputs,
            &settings,
            &mut sink,
            &mut posted,
        )
    };
    let Done {
        mut writers,
        records,
    } = match searched {
        Ok(done) => done,
        Err(message) => {
            if let Some(text) = posted {
                started.post(&text);
            }
            return Err(message);
        }
    };
    for writer in &mut writers {
        writer.finish();
    }
    diagnostics.finish();
    if let Some((hits, recorder, ordinals)) = records {
        emit_hit_records(&hits, &recorder, ordinals.as_deref());
    }
    Ok(())
}

/// The HSP records of the tabular search (`emit_hit_records`): its hits, the recorder of their
/// ranges and the hits' ordinals.
type HitRecords = (Arc<Vec<PairwiseHit>>, RangeRecorder, Option<Vec<usize>>);

/// What the searches of a run wrote (`search_all`): the streams of the formats and the HSP records
/// of the tabular search.
struct Done {
    writers: Vec<StreamWriter>,
    records: Option<HitRecords>,
}

/// The searches of a run of the parsed argv `parsed` (`run`), whose file of
/// `-export_search_strategy` is `strategy` (`started_command`): the argv's own search writes it as
/// the command line's run of the argv does (`losat::web_adapter::give_strategy`), and each other
/// search writes one of the same name whose bytes no one keeps, so that it reads its records and
/// refuses the file as the argv's search does. The messages that the command line
/// writes to its standard error for the argv go to `diagnostics`, and what the other formats'
/// standard error adds to them goes to stream 3 after them (`CallerSink::add`); for an error that
/// ends the run, `posted` gets the NCBI message that the command line writes for it
/// (`Failed::native`). Each output stream that fails while the run goes on (a format that failed
/// alone, the formats of the other search when it failed) is signalled to the host with the place
/// where the sink wrote its error (`CallerSink::error`, `emit::discard`).
fn search_all(
    parsed: (Program, Commands, Requested),
    mut strategy: Option<StrategyFile>,
    words: &[&str],
    inputs: RunInputs,
    settings: &ApplicationSettings,
    diagnostics: &mut CallerSink<'_>,
    posted: &mut Option<Vec<u8>>,
) -> Result<Done, String> {
    use std::io::Write as _;
    let (program, command, requested) = parsed;
    let groups = search_groups(&outputs(program, &requested), &command);
    let own = own_search(&groups);
    let handlers: Vec<Vec<Vec<u8>>> = groups
        .iter()
        .map(|group| {
            group
                .iter()
                .map(|format| handler_warnings(program, &command, format.number))
                .collect()
        })
        .collect();
    let threads = num_threads(&command);
    let counts = window_masker_counts(&command);
    let mut command = Some(command);
    // The warnings of the argv's own search go to `diagnostics` whatever the outcome: as the first
    // search it writes them as it goes (with a copy, whose lines the other formats' repeat,
    // `Diagnostics`), else they are held until the searches end. Those of the other search are
    // stream 3's only when the run completes: they are held until the searches end.
    let streamed = own == 0;
    // Each search: what it wrote or why it failed, and the warnings that it wrote to the standard
    // error.
    let mut outcomes: Vec<(Result<Searched, Failed>, Vec<u8>)> = Vec::new();
    // The search whose error ends the run: the argv's own search (`own_search`), or one with a
    // fault of the adapter. No search follows it.
    let mut ending = None;
    // The searches of the run share one thread pool (`with_run_pool`): a browser holds a thread
    // worker for every thread that may start, and the threads of a pool may not have returned when
    // the next pool starts, so each pool of its own needed thread workers of its own.
    losat::web_adapter::with_run_pool(threads, program.name(), || -> Result<(), String> {
        for (position, group) in groups.iter().enumerate() {
            let mut command = match command.take() {
                Some(command) => command,
                None => parse(words)?.1,
            };
            give_run_resources(&mut command, counts.as_ref(), settings);
            let file = match position == own {
                true => strategy.take(),
                false => losat::web_adapter::open_strategy_file(
                    &command,
                    &words[1..],
                    &RunFiles::Receiver(Arc::new(NoFiles(SearchRunFiles::default()))),
                )
                .map_err(|error| format!("{error:#}"))?,
            };
            if let Some(file) = file {
                losat::web_adapter::give_strategy(
                    &mut command,
                    file,
                    requested.number.map(|number| number as i32),
                );
            }
            let mut warnings = Vec::new();
            let result = match streamed && position == 0 {
                true => search(
                    command,
                    group,
                    inputs,
                    settings,
                    &mut Kept {
                        sink: &mut *diagnostics,
                        copy: &mut warnings,
                    },
                ),
                false => search(command, group, inputs, settings, &mut warnings),
            };
            let ends = matches!(&result, Err(failed) if position == own || failed.fault);
            outcomes.push((result, warnings));
            if ends {
                ending = Some(position);
                break;
            }
        }
        Ok(())
    })?;
    // A run that fails: stream 3 is the CLI's standard error for the argv up to its error, which
    // is the run's error and not stream 3's: the warnings of the argv's own search alone, which the
    // CLI's run of the argv writes, not those of a search before it, which writes formats that the
    // argv does not ask for (line WQA's sweep, D3 review M1, 2026-11-19).
    if let Some(position) = ending {
        let (result, warnings) = outcomes.swap_remove(position);
        if position == own && !streamed {
            let _ = diagnostics.write_all(&warnings);
        }
        return Err(match result {
            Err(failed) => {
                keep_failed_outputs(&groups, failed.kept, failed.fault);
                *posted = failed.native;
                failed.message
            }
            Ok(_) => unreachable!("the search that ends the run failed"),
        });
    }
    // A run that completes: where the run's messages go (stream 3, or TBLASTX's file of `-logfile`)
    // has the CLI's standard error for the argv, the warnings of the argv's own search. What the
    // standard error of the run's other formats adds to it goes to stream 3 after it, or nowhere
    // where the run's messages go nowhere (`CallerSink::add`), each text once (`Diagnostics`), for
    // the argv's own search and then the other: the lines of the other search's warnings that stream
    // 3 does not hold yet, then, format by format in the search's order, what the CLI's run of
    // that format adds: the formatting handler's warnings of a format after the first of the
    // search (which the search does not write: `handler_warnings`), the queries' warnings that it
    // writes where the search's first format, a structured one, writes none (`Searched::later`),
    // and its refusal if it failed alone; a search that failed alone (its formats are empty) adds
    // the handler's warnings and its error (line WQA's sweep, D3 review r2, 2026-11-19; the order by
    // format: review r1 of the later formats' warnings, L1, 2026-11-20). Each stream that failed,
    // alone or with its search, is signalled with the place where the sink wrote its error (stream
    // 9, also where no chunk of it was sent: the host tells which outputs a run could not write,
    // and where their errors are, from these signals alone; line WQA's sweep, D7, 2026-11-19 and
    // 2026-11-20).
    if !streamed {
        let _ = diagnostics.write_all(&outcomes[own].1);
    }
    // The lines of the argv's own search, which the other formats' repeat.
    let mut texts = Diagnostics::after(&outcomes[own].1);
    let others = (0..outcomes.len()).filter(|&position| position != own);
    for position in std::iter::once(own).chain(others) {
        let (result, warnings) = &outcomes[position];
        if position != own {
            texts.warnings(warnings, diagnostics);
        }
        match result {
            Err(failed) => {
                for handler in &handlers[position][1..] {
                    texts.warnings(handler, diagnostics);
                }
                let place = texts.error(failed, diagnostics);
                for output in &groups[position] {
                    crate::emit::discard(output.stream, place);
                }
            }
            Ok(done) => {
                for (index, output) in groups[position].iter().enumerate() {
                    if index > 0 {
                        texts.warnings(&handlers[position][index], diagnostics);
                    }
                    // The queries' warnings that a structured first format does not write and
                    // this format's run writes, as far as its report went (line WQA's sweep, D3
                    // review r3, 2026-11-19: they were lost).
                    texts.warnings(&done.later[index], diagnostics);
                    let refused = done
                        .refused
                        .iter()
                        .find(|(stream, _)| *stream == output.stream);
                    if let Some((stream, refused)) = refused {
                        let place = texts.error(refused, diagnostics);
                        crate::emit::discard(*stream, place);
                    }
                }
            }
        }
    }
    let mut searched: Vec<Searched> = Vec::new();
    let mut tabular_failed = false;
    for (position, (result, _)) in outcomes.into_iter().enumerate() {
        match result {
            Ok(done) => searched.push(done),
            Err(_) => tabular_failed |= position == 0,
        }
    }
    let mut writers = Vec::new();
    let mut records = None;
    if !tabular_failed {
        let mut first = searched.remove(0);
        if let Some(hits) = &first.hits {
            for other in searched.iter().filter(|other| !other.refused(0)) {
                if let Some(second) = &other.hits {
                    link_format0(hits, &mut first.recorder, second, &other.recorder);
                }
            }
        }
        writers.append(&mut first.writers);
        // The HSP records point into stream 6 (§8): none where it failed alone.
        if !first.refused(6) {
            records = first
                .hits
                .map(|hits| (hits, first.recorder, first.ordinals));
        }
    }
    for other in &mut searched {
        writers.append(&mut other.writers);
    }
    Ok(Done { writers, records })
}

/// The outputs of a run that fails (line WQA's sweep, D15, 2026-11-20; docs/web/abi_v2.md §5): the
/// output of the argv's own `-outfmt` (`Output::requested`) is what the CLI's run of the argv writes,
/// and its file holds the bytes written before the error, so its writer (`Failed::kept`) sends what it
/// holds, and the host keeps what it received. The other outputs, which the Web writes besides
/// (outfmt 0, 6 and 7), are not the CLI's run of the argv, and a format of the search that failed
/// stopped where that search did: each is signalled as dropped, with the error that ended the run as
/// the place of its error (`emit::PLACE_RUN`), also where no search made it. An adapter fault
/// (`Failed::fault`) is no CLI error: it keeps no output.
fn keep_failed_outputs(groups: &[Vec<Output>], kept: Vec<StreamWriter>, fault: bool) {
    for mut writer in kept {
        writer.finish();
    }
    for output in groups.iter().flatten() {
        if fault || !output.requested {
            crate::emit::discard(output.stream, crate::emit::PLACE_RUN);
        }
    }
}

/// A writer to `sink` that keeps a copy of what it writes in `copy` (`search_all`: the warnings of
/// the argv's own search as the first, whose lines the other formats' repeat).
struct Kept<'w> {
    sink: &'w mut (dyn std::io::Write + Send),
    copy: &'w mut Vec<u8>,
}

impl std::io::Write for Kept<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.sink.write(buf)?;
        self.copy.extend_from_slice(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.sink.flush()
    }
}

/// Stream 3 of a run that completes (`search_all`, docs/web/abi_v2.md §5): the CLI's standard
/// error for the argv, then what that of the other formats adds, where each text comes once. The
/// other search writes only the lines of its warnings that stream 3 does not hold yet: the reader's
/// and the queries' warnings, which every search of the run writes (up to the place where its
/// formats end), come with the argv's search, and so does a formatting handler's line that both
/// searches write (`Examining 5 or more matches`). Each format after the argv's own one then
/// adds (in the other search, its first format too), in the search's order, its formatting handler's lines, what its run writes with its
/// queries' reports where the argv's format is a structured one (`Searched::later`: the queries'
/// warnings; BLASTX: every line of its reports' standard error), and its refusal where it failed
/// alone (every format of the other search fails alone), each where stream 3 does not hold it;
/// the error of a search that failed follows that
/// search's warnings. An error or a refusal comes once as a whole text (streams 6 and 7 of a
/// search stop at the same row). What is added goes where the sink writes it (`CallerSink::add`,
/// `CallerSink::error`).
struct Diagnostics<'t> {
    /// The lines that stream 3 holds, with the number of times it holds each.
    lines: std::collections::HashMap<&'t [u8], usize>,
    /// The errors that the run holds, each with the place where the sink wrote it (`emit::place`).
    errors: Vec<(&'t [u8], u32)>,
}

impl<'t> Diagnostics<'t> {
    /// Stream 3 after the lines `first` of the argv's own search (`search_all`).
    fn after(first: &'t [u8]) -> Self {
        let mut texts = Diagnostics {
            lines: std::collections::HashMap::new(),
            errors: Vec::new(),
        };
        texts.hold(first);
        texts
    }

    /// Stream 3 holds `text` (one more time for each of its lines).
    fn hold(&mut self, text: &'t [u8]) {
        for line in text.split_inclusive(|&byte| byte == b'\n') {
            *self.lines.entry(line).or_default() += 1;
        }
    }

    /// The warnings `text` of the other search, or of the formatting handler for a format after
    /// the first of its search: a line is added by `sink` (`CallerSink::add`) where `text` has it
    /// more times than stream 3 holds it (its n-th copy of a line repeats the n-th copy that
    /// stream 3 holds).
    fn warnings(&mut self, text: &'t [u8], sink: &mut CallerSink<'_>) {
        let mut copies: std::collections::HashMap<&'t [u8], usize> =
            std::collections::HashMap::new();
        for line in text.split_inclusive(|&byte| byte == b'\n') {
            let copy = copies.entry(line).or_default();
            *copy += 1;
            if *copy > self.lines.get(line).copied().unwrap_or(0) {
                sink.add(line);
            }
        }
        for (line, copies) in copies {
            let held = self.lines.entry(line).or_default();
            *held = (*held).max(copies);
        }
    }

    /// An error of a search, or a refusal of a format that failed alone, `failed`: written by
    /// `sink` where the CLI writes it (`CallerSink::error`) unless the run holds it, and the place
    /// where it is, which the sink gave when it wrote it (`emit::place`).
    fn error(&mut self, failed: &'t Failed, sink: &mut CallerSink<'_>) -> u32 {
        let text = failed.cli_stderr.as_slice();
        if let Some((_, place)) = self.errors.iter().find(|(held, _)| *held == text) {
            return *place;
        }
        let place = crate::emit::place(sink.error(text, failed.source()));
        self.errors.push((text, place));
        self.hold(text);
        place
    }
}

/// Emits one JSON object per HSP on stream 1 (docs/web/abi_v2.md §8); `ordinals` maps the
/// engine's `q_idx` (searched queries) to the registered index (`searched_queries`).
fn emit_hit_records(hits: &[PairwiseHit], recorder: &RangeRecorder, ordinals: Option<&[usize]>) {
    let mut lines = String::new();
    let mut ranks: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for (index, record) in hits.iter().enumerate() {
        let hit = &record.hit;
        let q_idx = ordinals
            .and_then(|ordinals| ordinals.get(hit.q_idx as usize).copied())
            .unwrap_or(hit.q_idx as usize);
        let rank = ranks.entry(q_idx).or_default();
        let _ = write!(
            lines,
            "{{\"index\":{index},\"q_idx\":{q_idx},\"s_idx\":{},\"rank\":{rank},\"raw_score\":{},\"bit_score\":",
            hit.s_idx, hit.raw_score
        );
        *rank += 1;
        json::number(&mut lines, hit.bit_score);
        lines.push_str(",\"e_value\":");
        json::number(&mut lines, hit.e_value);
        let _ = write!(
            lines,
            ",\"q_start\":{},\"q_end\":{},\"s_start\":{},\"s_end\":{},\"query_frame\":",
            hit.q_start, hit.q_end, hit.s_start, hit.s_end
        );
        json::optional(&mut lines, record.query_frame);
        lines.push_str(",\"subject_frame\":");
        json::optional(&mut lines, record.subject_frame);
        lines.push_str(",\"subject_length\":");
        json::optional(&mut lines, record.subject_length);
        lines.push_str(",\"query_aligned\":");
        json::optional_string(&mut lines, record.query_seq.as_deref());
        lines.push_str(",\"subject_aligned\":");
        json::optional_string(&mut lines, record.subject_seq.as_deref());
        lines.push_str(",\"out6\":");
        json::range(&mut lines, recorder.out6.get(index).copied().flatten());
        lines.push_str(",\"out0\":");
        json::range(&mut lines, recorder.out0.get(index).copied().flatten());
        lines.push_str(",\"out0_subject\":");
        json::range(
            &mut lines,
            recorder.out0_subject.get(index).copied().flatten(),
        );
        lines.push_str("}\n");
        if lines.len() >= crate::emit::CHUNK {
            emit(STREAM_HITS, lines.as_bytes());
            lines.clear();
        }
    }
    emit(STREAM_HITS, lines.as_bytes());
}

#[cfg(test)]
#[path = "run_start_tests.rs"]
mod start_tests;

#[cfg(test)]
#[path = "run_import_tests.rs"]
mod import_tests;

#[cfg(test)]
#[path = "run_strategy_tests.rs"]
mod strategy_tests;

#[cfg(test)]
mod tests {
    use super::*;

    /// NCBI's text in LOSAT's words (`losat::web_adapter::losat_text`): the tests state what NCBI
    /// writes and compare it with what LOSAT writes.
    fn ncbi_in_losat_words(ncbi: &str) -> String {
        losat::web_adapter::losat_text(ncbi).into_owned()
    }

    /// A frozen NCBI output or standard error in LOSAT's words (`losat::web_adapter::losat_output_bytes`).
    fn ncbi_output_in_losat_words(ncbi: &[u8]) -> Vec<u8> {
        losat::web_adapter::losat_output_bytes(ncbi).into_owned()
    }

    #[test]
    fn argv_is_parsed_as_on_the_command_line() {
        let words = ["blastp", "-query", "q", "-subject", "s", "-evalue", "0.1"];
        assert!(matches!(
            parse(&words),
            Ok((Program::Blastp, Commands::Blastp(_)))
        ));
        for owned in ["-out", "-out=x"] {
            let error = parse(&["blastp", "-query", "q", "-subject", "s", owned, "6"]).unwrap_err();
            assert!(error.contains("is not accepted"), "{owned}: {error}");
        }
        // A help word prints the help, as on the command line (line WQA's D14 review).
        for help in [["blastp", "-out", "x", "-h"], ["blastp", "-h", "-out", "x"]] {
            assert!(
                matches!(parse_printing(&help, None), Ok(Parsed::Prints(_))),
                "{help:?}"
            );
        }
        assert!(matches!(
            parse_printing(&["blastp", "-out=x", "-help"], None),
            Ok(Parsed::Prints(_))
        ));
        // -num_descriptions and -num_alignments are accepted (also as `-x=v`),
        // parsed as on the command line; -max_target_seqs excludes both (NCBI's reader's
        // exclusion).
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            for option in [
                &["-num_descriptions", "3"][..],
                &["-num_alignments", "0"],
                &["-num_alignments=3"],
            ] {
                let words = [&[program, "-query", "q", "-subject", "s"][..], option].concat();
                assert!(parse(&words).is_ok(), "{program} {option:?}");
                assert!(validate(&words).is_ok(), "{program} {option:?}");
            }
            let words = [
                program,
                "-query",
                "q",
                "-subject",
                "s",
                "-max_target_seqs",
                "3",
            ];
            assert!(parse(&words).is_ok(), "{program}");
            for option in ["-num_descriptions", "-num_alignments"] {
                let both = [&words[..], &[option, "3"]].concat();
                let cli = losat::web_adapter::try_parse_cli(
                    std::iter::once("losat").chain(both.iter().copied()),
                    None,
                )
                .unwrap_err();
                assert_eq!(
                    parse(&both).unwrap_err(),
                    losat::web_adapter::argv_message(&cli),
                    "{program} {option}"
                );
                assert!(parse(&both)
                    .unwrap_err()
                    .contains("Incompatible with argument:  `"));
            }
        }
        assert!(matches!(
            parse(&["blastx", "-query", "q", "-subject", "s"]),
            Ok((Program::Blastx, Commands::Blastx(_)))
        ));
        assert!(parse(&["psiblast", "-query", "q", "-subject", "s"]).is_err());
        // An unknown program and a parser error are the CLI's messages.
        for words in [
            &["nosuch"][..],
            &["tblastx", "-query"][..],
            &["tblastx", "-query", "q", "-subject", "s"][..],
        ] {
            let cli = losat::web_adapter::try_parse_cli(
                std::iter::once("losat").chain(words.iter().copied()),
                None,
            )
            .map(|_| ())
            .map_err(|error| losat::web_adapter::argv_message(&error));
            assert_eq!(parse(words).map(|_| ()), cli, "{words:?}");
        }
        // An argv without subjects is NCBI's error after the parse (line WQA's D14, 2026-11-21).
        assert_eq!(
            parse(&["tblastx", "-query", "q"]).map(|_| ()),
            Err(format!("{:#}", losat::web_adapter::missing_subject_error()))
        );
        assert_eq!(Program::Tblastx.formats(), &[0, 6, 7]);
        assert_eq!(Program::Blastx.formats(), &[0, 6, 7]);
        // validate rejects BLASTN scoring that NCBI rejects, with NCBI's message.
        // validate reads BLASTN's -dust as the run does.
        let error =
            validate(&["blastn", "-query", "q", "-subject", "s", "-dust", "20 64"]).unwrap_err();
        assert!(
            error.starts_with(&ncbi_in_losat_words(
                "BLAST query/options error: Invalid number of arguments"
            )),
            "{error}"
        );
        let blastn = ["blastn", "-query", "q", "-subject", "s", "-task", "blastn"];
        assert!(validate(&blastn).is_ok());
        let unsupported = [&blastn[..], &["-reward", "1", "-penalty", "-6"]].concat();
        assert!(
            parse(&unsupported).is_ok(),
            "run leaves the check to the engine"
        );
        let error =
            validate(&[&blastn[..], &["-gapopen", "0", "-gapextend", "0"]].concat()).unwrap_err();
        assert!(
            error.starts_with(&ncbi_in_losat_words(
                "BLAST query/options error: Greedy extension"
            )),
            "{error}"
        );
        // The tables depend on the query's letters, so `run` checks them (WA-FM-01,
        // `blastn_scores_without_table_row_are_left_to_the_run`).
        assert!(validate(&unsupported).is_ok());
        let letters = b"ACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGA";
        let error = run(
            &unsupported,
            &[FastaRecord::new("Query_1", b"q", letters)],
            &losat::api::QueryEnd::Input,
            &[FastaRecord::new("Subject_1", b"s", letters)],
        )
        .unwrap_err();
        assert!(
            error.starts_with(&ncbi_in_losat_words(
                "BLAST engine error: Error: Substitution scores 1 and -6"
            )),
            "{error}"
        );
        // dc-megablast and blastn-short (Session SD): the task defaults and NCBI's template
        // checks, with NCBI's messages (blast_options.c:1247-1261,1399-1413).
        for task in ["dc-megablast", "blastn-short"] {
            assert!(validate(&["blastn", "-query", "q", "-subject", "s", "-task", task]).is_ok());
        }
        let dc = [
            "blastn",
            "-query",
            "q",
            "-subject",
            "s",
            "-task",
            "dc-megablast",
        ];
        let template = [
            "-template_type",
            "coding_and_optimal",
            "-template_length",
            "21",
        ];
        assert!(validate(&[&dc[..], &template[..], &["-word_size", "12"]].concat()).is_ok());
        let error = validate(&[&dc[..], &["-word_size", "13"]].concat()).unwrap_err();
        assert!(
            error.starts_with(
                &ncbi_in_losat_words("BLAST query/options error: Invalid discontiguous template parameters: word size must be either 11 or 12")
            ),
            "{error}"
        );
        let error =
            validate(&[&blastn[..], &template[..], &["-word_size", "11"]].concat()).unwrap_err();
        assert!(
            error.starts_with(&ncbi_in_losat_words(
                "BLAST query/options error: Invalid lookup table type for discontiguous Mega BLAST"
            )),
            "{error}"
        );
        // The options require each other: a parser error, the CLI's message.
        assert!(validate(&[&dc[..], &["-template_type", "coding"]].concat()).is_err());
        // validate reports TBLASTX's checks of its options alone, as the run does; the one-hit
        // word finder (`-window_size 0`) runs as the CLI runs it (line RG, X2).
        let tblastx = ["tblastx", "-query", "q", "-subject", "s"];
        assert!(validate(&tblastx).is_ok());
        assert!(validate(&[&tblastx[..], &["-window_size", "0"]].concat()).is_ok());
        let error = validate(&[&tblastx[..], &["-evalue", "0"]].concat()).unwrap_err();
        assert!(
            error.starts_with(&ncbi_in_losat_words(
                "BLAST query/options error: expect value"
            )),
            "{error}"
        );
        // validate reports BLASTX's checks of the argv's own search's options, as the run does
        // (NCBI's errors and LOSAT's refusals); the tabular search that -num_alignments 0 rejects
        // alone fails only its formats where the argv's own is outfmt 0, as in `run`.
        let blastx = ["blastx", "-query", "q", "-subject", "s"];
        assert!(validate(&blastx).is_ok());
        assert!(validate(&[&blastx[..], &["-num_alignments", "0"]].concat()).is_ok());
        assert!(
            validate(&[&blastx[..], &["-num_alignments", "0", "-outfmt", "7"]].concat())
                .unwrap_err()
                .starts_with(&ncbi_in_losat_words(
                    "BLAST query/options error: No hits are being saved"
                ))
        );
        for (extra, start) in [
            (
                &["-evalue", "0"][..],
                ncbi_in_losat_words("BLAST query/options error: expect value"),
            ),
            (
                &["-matrix", "IDENTITY"],
                ncbi_in_losat_words(
                    "BLAST query/options error: IDENTITY is not a supported matrix",
                ),
            ),
            (
                &["-num_descriptions", "0", "-num_alignments", "0"],
                ncbi_in_losat_words("BLAST query/options error: No hits are being saved"),
            ),
        ] {
            let error = validate(&[&blastx[..], extra].concat()).unwrap_err();
            assert!(error.starts_with(&start), "{extra:?}: {error}");
        }
        // An error of the CLI parser keeps the CLI's message.
        for words in [
            &["blastp", "-query", "q", "-subject", "s", "-evalue"][..],
            &[
                "blastn",
                "-query",
                "q",
                "-subject",
                "s",
                "-nosuchoption",
                "1",
            ][..],
            &[
                "tblastn",
                "-query",
                "q",
                "-subject",
                "s",
                "-db_gencode",
                "7",
            ][..],
        ] {
            let cli = losat::web_adapter::try_parse_cli(
                ["losat"].into_iter().chain(words.iter().copied()),
                None,
            )
            .unwrap_err();
            assert_eq!(
                parse(words).unwrap_err(),
                losat::web_adapter::argv_message(&cli),
                "{words:?}"
            );
        }
    }

    /// The text a run wrote on `stream`.
    fn emitted(stream: u32) -> String {
        let bytes = crate::emit::EMITTED.with(|emitted| {
            emitted
                .borrow()
                .iter()
                .filter(|(number, _)| *number == stream)
                .flat_map(|(_, bytes)| bytes.clone())
                .collect::<Vec<u8>>()
        });
        String::from_utf8(bytes).unwrap()
    }

    // WA-FM-01. The CLI meets the Karlin-Altschul tables (scores and gap costs without a row,
    // LOSAT's limit of greedy gap costs) only when some query context has an ungapped
    // Karlin-Altschul block (`context_blocks`), which depends on the query's letters: a query
    // of N only has none, and the same scores give a warning and rc 0. `validate` has no query,
    // so it makes the checks of the options only, and `run` rejects with the CLI's message and
    // exit code (`validate` accepts a subset of what the CLI accepts).
    #[test]
    fn blastn_scores_without_table_row_are_left_to_the_run() {
        let megablast = ["blastn", "-query", "q", "-subject", "s"];
        let blastn = [&megablast[..], &["-task", "blastn"]].concat();
        let with = |base: &[&str], extra: &[&str]| validate(&[base, extra].concat());
        // The scores of the symptom: no table row, and a pair that the CLI rejects for some
        // queries.
        for base in [&megablast[..], &blastn[..]] {
            for pair in [
                &["-reward", "7"][..],
                &["-reward", "7", "-penalty", "-2"],
                &["-reward", "3", "-penalty", "-1"],
                &["-reward", "4", "-penalty", "-1"],
                &["-reward", "5", "-penalty", "-1"],
                &["-reward", "2", "-penalty", "-1"],
            ] {
                assert_eq!(with(base, pair), Ok(true), "{base:?} {pair:?}");
            }
        }
        // The checks of the options stand: NCBI's and LOSAT's refusals.
        let error = with(&megablast, &["-reward", "7", "-evalue", "0"]).unwrap_err();
        assert!(
            error.starts_with(&ncbi_in_losat_words(
                "BLAST query/options error: expect value"
            )),
            "{error}"
        );
        let gaps = [
            "-reward",
            "3",
            "-penalty",
            "-1",
            "-gapopen",
            "5",
            "-gapextend",
            "0",
        ];
        let error = with(&blastn, &gaps).unwrap_err();
        assert!(
            error.contains("gap extension penalty cannot be 0"),
            "{error}"
        );
        let error = with(&megablast, &["-reward", "0"]).unwrap_err();
        assert!(error.contains("a reward of 0"), "{error}");

        // `run` meets the tables with the query's letters, as the CLI does.
        let search = |base: &[&str], extra: &[&str], letters: &[u8]| {
            crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
            let query = [FastaRecord::new("Query_1", b"mq1 query", letters)];
            let subject = [FastaRecord::new("Subject_1", b"ms1", ACGT)];
            let result = run(
                &[base, extra].concat(),
                &query,
                &losat::api::QueryEnd::Input,
                &subject,
            );
            (
                result,
                emitted(crate::emit::STREAM_DIAGNOSTICS),
                emitted(crate::emit::STREAM_HITS),
            )
        };
        const ACGT: &[u8] = b"ACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGA";
        const NNNN: &[u8] = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        const WARNING: &str =
            "Warning: [blastn] Query_1 mq1 query: Could not calculate ungapped Karlin-Altschul";
        // An expected score that is not negative (no block for ACGT): a warning, no hit.
        for base in [&megablast[..], &blastn[..]] {
            for pair in [
                &["-reward", "3", "-penalty", "-1"][..],
                &["-reward", "5", "-penalty", "-1"],
            ] {
                let (result, diagnostics, hits) = search(base, pair, ACGT);
                assert_eq!(result, Ok(()), "{base:?} {pair:?}");
                assert!(diagnostics.starts_with(WARNING), "{diagnostics}");
                assert_eq!(hits, "");
            }
        }
        let (result, diagnostics, _) = search(&megablast, &["-reward", "7"], ACGT);
        assert_eq!(result, Ok(()));
        assert!(diagnostics.starts_with(WARNING), "{diagnostics}");
        // A query of N only has no block for any scores: `-task blastn -reward 7` is a warning.
        let (result, diagnostics, hits) = search(&blastn, &["-reward", "7"], NNNN);
        assert_eq!(result, Ok(()));
        assert!(diagnostics.starts_with(WARNING), "{diagnostics}");
        assert_eq!(hits, "");
        // A query of A, C, G and T has a block for 7 and -3, which has no row: the CLI's error.
        let (result, _, _) = search(&blastn, &["-reward", "7"], ACGT);
        let error = result.unwrap_err();
        assert_eq!(
            error.trim_end(),
            ncbi_in_losat_words(
                "BLAST engine error: Error: Substitution scores 7 and -3 are not supported"
            )
        );
        // Letters of ambiguity put a block under scores whose expected score is 0: (3, -1) with
        // a B (LOSAT's own limit of the gap costs, or NCBI's table error, as the CLI).
        let mut ambiguous = ACGT.to_vec();
        ambiguous.extend_from_slice(b"B");
        let (result, _, _) = search(&megablast, &["-reward", "3", "-penalty", "-1"], &ambiguous);
        assert!(result.is_err());
    }

    fn recorder() -> RangeRecorder {
        RangeRecorder {
            streams: vec![0, 6],
            positions: vec![Arc::new(AtomicU64::new(0)), Arc::new(AtomicU64::new(0))],
            open: vec![None, None],
            out6: Vec::new(),
            out0: Vec::new(),
            out0_subject: Vec::new(),
            subject_start: None,
            subject: None,
            fault: None,
        }
    }

    // `-parse_deflines` is the CLI's option (line WQA's sweep, D7, 2026-11-19): `validate` takes it,
    // and the host gives a run of it the inputs of `register_inputs` (`inputs_together`; the
    // *start* JSON's `register_inputs`), but not a dry run, which reads no input.
    #[test]
    fn parse_deflines_reads_the_inputs_together() {
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            let words = [program, "-query", "q", "-subject", "s", "-parse_deflines"];
            assert_eq!(validate(&words), Ok(true));
            assert!(inputs_together(&words), "{program}");
            assert!(start_json(&words)
                .unwrap()
                .ends_with(",\"register_inputs\":true,\"prints\":null,\"database\":false}"));
            let dry = [&words[..], &["-dryrun"]].concat();
            assert!(!inputs_together(&dry), "{program}");
            assert!(start_json(&dry)
                .unwrap()
                .ends_with(",\"register_inputs\":false,\"prints\":null,\"database\":false}"));
            assert!(!inputs_together(&words[..5]), "{program}");
            assert!(!inputs_together(&[&words[..], &["-nosuchoption"]].concat()));
        }
    }

    /// The `shown_id` of each record of the *register* response `response` (no escape in them).
    fn shown_ids(response: &str) -> Vec<String> {
        response
            .split("\"shown_id\":\"")
            .skip(1)
            .map(|rest| rest[..rest.find('"').unwrap()].to_string())
            .collect()
    }

    // A run with `-parse_deflines` writes NCBI's bytes (line WQA's sweep, D7, 2026-11-19): every row
    // of NCBI's frozen outputs of the five programs with the option (LOSAT/tests/outfmt0_manifest.tsv,
    // rows pd.* and the others), in its format, its inputs registered together for its argv
    // (`registered_streams`), and stream 3 NCBI's standard error. The IDs that the *inputs*
    // response shows for the records are the first two columns of outfmt 6 (`qaccver`, `saccver`),
    // which the application's results show.
    #[test]
    fn parse_deflines_runs_write_ncbis_outputs() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        let rows = manifest_rows_with(|program, _, files, extra| {
            Program::parse(program).is_ok()
                && files.is_empty()
                && extra
                    .split_whitespace()
                    .any(|word| word == "-parse_deflines")
        });
        // Not a PSSM query (a query column `pssm-gz:`, `-in_pssm`), which LOSAT Web does not take.
        let rows: Vec<_> = rows
            .into_iter()
            .filter(|(_, argv)| !argv[2].starts_with("pssm"))
            .collect();
        let programs: std::collections::BTreeSet<&str> =
            rows.iter().map(|(_, argv)| argv[0].as_str()).collect();
        assert_eq!(programs.len(), 5, "{programs:?}");
        assert!(rows.len() >= 130, "{}", rows.len());
        let (mut shown, mut gaps, mut alone) = (0, 0, 0);
        for (id, argv) in &rows {
            let words: Vec<&str> = argv.iter().map(String::as_str).collect();
            let at = |flag: &str| words[words.iter().position(|w| *w == flag).unwrap() + 1];
            let (query, subject) = (read(at("-query")), read(at("-subject")));
            // A `>?` gap line runs as in the CLI (line WQA's D13 (a), 2026-11-20; the store refused
            // it before): its rows are counted, to show that they run.
            let gap_line = [&query, &subject].iter().any(|bytes| {
                bytes
                    .split(|b| *b == b'\n')
                    .any(|line| line.starts_with(b">?"))
            });
            gaps += usize::from(gap_line);
            // The pd.self rows other than outfmt 0: a record that is both a query and a subject and
            // has masks, whose outfmt 0 LOSAT refuses, that format's alone
            // (`believed::refuse_shared_masks_in`; line WQA's sweep, D7, 2026-11-19): the run
            // writes the argv's format and outfmt 0 fails alone, its refusal on stream 3 after the
            // CLI's warnings.
            let shared_masks = id.starts_with("pd.self.") && !id.ends_with(".0");
            let streams = try_registered_streams(&words, &query, &subject)
                .unwrap_or_else(|error| panic!("{id}: {error}"));
            let outfmt = at("-outfmt");
            let stream = match outfmt {
                "0" => 0,
                "6" => 6,
                "7" => 7,
                _ => STREAM_CUSTOM,
            };
            let ncbi =
                ncbi_output_in_losat_words(&read(&format!("tests/fixtures/outfmt0/{id}.out")));
            let ours = streams.get(&stream).cloned().unwrap_or_default();
            assert!(ours == ncbi, "{id}: stream {stream}\n{}", text(&ours));
            let stderr = ncbi_output_in_losat_words(
                &std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                    .unwrap_or_default(),
            );
            let diagnostics = streams
                .get(&STREAM_DIAGNOSTICS)
                .cloned()
                .unwrap_or_default();
            if shared_masks {
                let program = match words[0] {
                    "blastn" => "BLASTN",
                    "blastp" => "BLASTP",
                    _ => "TBLASTX",
                };
                let refusal = format!(
                    "is both a query and a subject under -parse_deflines and has masks (NCBI \
                     BLAST+ draws them at the other record's places in outfmt 0); this is not \
                     supported by LOSAT's {program}\n"
                );
                let rest = text(&diagnostics)
                    .strip_prefix(text(&stderr).as_str())
                    .map(str::to_string)
                    .unwrap_or_else(|| panic!("{id}: stream 3 {}", text(&diagnostics)));
                assert!(
                    rest.starts_with("Error: the record ") && rest.ends_with(&refusal),
                    "{id}: stream 3 after the CLI's: {rest}"
                );
                assert!(
                    streams.get(&0).is_none_or(Vec::is_empty),
                    "{id}: outfmt 0 failed alone"
                );
                alone += 1;
            } else {
                assert_eq!(text(&diagnostics), text(&stderr), "{id}: stream 3");
            }
            if outfmt != "6" {
                continue;
            }
            let (query_handle, subject_handle, response) =
                crate::store::register_inputs(&words, &query, &subject).unwrap();
            crate::store::release(query_handle).unwrap();
            crate::store::release(subject_handle).unwrap();
            let split = response.find(",\"subject\":").unwrap();
            let (query_ids, subject_ids) =
                (shown_ids(&response[..split]), shown_ids(&response[split..]));
            let hits = text(streams.get(&STREAM_HITS).map_or(&[][..], Vec::as_slice));
            let index = |line: &str, field: &str| {
                let at = line.find(field).unwrap() + field.len();
                line[at..]
                    .split(',')
                    .next()
                    .unwrap()
                    .parse::<usize>()
                    .unwrap()
            };
            for (hit, row) in hits.lines().zip(text(&ncbi).lines()) {
                let columns: Vec<&str> = row.split('\t').collect();
                assert_eq!(
                    query_ids[index(hit, "\"q_idx\":")],
                    columns[0],
                    "{id}: {hit}"
                );
                assert_eq!(
                    subject_ids[index(hit, "\"s_idx\":")],
                    columns[1],
                    "{id}: {hit}"
                );
                shown += 1;
            }
        }
        assert!(shown > 0);
        assert!(gaps > 0, "the rows of a gap line run");
        assert_eq!(alone, 33, "the pd.self rows other than outfmt 0");
        eprintln!(
            "{} rows, {gaps} with a gap line, {alone} with outfmt 0 failed alone, {shown} outfmt 6 rows' ids",
            rows.len()
        );
    }

    // The errors of a run with `-parse_deflines` are the CLI's, in its order (line WQA's sweep, D7,
    // 2026-11-19; the texts are those of the native CLI of this branch for the same inputs, less its
    // `Error: ` and line end): the subjects are read before the query, so a subject's error comes
    // first; a query with a gap line runs (line WQA's D13 (a), 2026-11-20); a subject input without
    // records ends the run before the query is read; a ranged id fails a run that writes the BLAST
    // archive alone. A dry run reads no input.
    #[test]
    fn parse_deflines_errors_are_the_clis() {
        let words = |outfmt: &'static str| {
            vec![
                "blastn",
                "-query",
                "q",
                "-subject",
                "s",
                "-parse_deflines",
                "-outfmt",
                outfmt,
            ]
        };
        let good: &[u8] = b">s1 x\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTT\n";
        let dup: &[u8] = b">a x\nACGTACGTACGTAAACCCGGGTTT\n>a y\nACGTACGTAAAACCCGGGTTTAAA\n";
        let no_defline: &[u8] = b"ACGTACGTACGT\n";
        let gap: &[u8] = b">q1 x\nACGT\n>?10\nACGT\n";
        let ranged: &[u8] = b">rs:3-20 ranged\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTT\n";
        let error = |outfmt, query: &[u8], subject: &[u8]| {
            try_registered_streams(&words(outfmt), query, subject).unwrap_err()
        };
        let dup_subject = "two subject records have the Seq-id lcl|a under -parse_deflines (NCBI \
                           BLAST+ shows the first for both); this is not supported by LOSAT's BLASTN";
        assert_eq!(
            error("6", dup, good),
            "two query records have the Seq-id lcl|a under -parse_deflines in one query batch, or \
             in any with -subject or outfmt 5 (NCBI BLAST+'s scope holds both and searches the \
             letters of the first for both); this is not supported by LOSAT's BLASTN"
        );
        assert_eq!(error("6", good, dup), dup_subject);
        assert_eq!(error("6", gap, dup), dup_subject);
        assert!(try_registered_streams(&words("6"), gap, good).is_ok());
        // LOSAT's refusal of ids that the query's gap line does not change: here a query and a
        // subject record of one id and other letters.
        let shared_gap: &[u8] = b">s1 x\nACGT\n>?10\nACGT\n";
        assert_eq!(
            error("6", shared_gap, good),
            "a query and a subject record have the Seq-id lcl|s1 under -parse_deflines (NCBI \
             BLAST+ shows them as one sequence); this is not supported by LOSAT's BLASTN"
        );
        // The query's reader error ends the run, which gives NCBI's errors with the CLI's line end
        // (line WQA's D10, 2026-11-19).
        assert_eq!(
            error("6", no_defline, good),
            ncbi_in_losat_words("BLAST query error: CFastaReader: Input doesn't start with a defline or comment around line 1\n")
        );
        assert_eq!(
            error("6", good, no_defline),
            ncbi_in_losat_words("BLAST query/options error: Query input doesn't start with a defline or comment, line \
             1\0\nPlease refer to the BLAST+ user manual.")
        );
        assert_eq!(
            error("6", no_defline, b"\n"),
            ncbi_in_losat_words("BLAST engine error: Empty CBlastQueryVector")
        );
        assert_eq!(
            error("11", good, ranged),
            "a subject Seq-id with a range under -parse_deflines in the BLAST archive is not \
             supported by LOSAT's BLASTN"
        );
        assert_eq!(
            error("11", ranged, good),
            "a query Seq-id with a range under -parse_deflines in the BLAST archive is not \
             supported by LOSAT's BLASTN"
        );
        // Without the archive the ranged record runs, and shows the id that NCBI generates for it.
        let streams = registered_streams(&words("6"), good, ranged);
        assert!(text(&streams[&6]).starts_with("s1\tranged\t100.000\t48\t"));
        let (query, subject, response) =
            crate::store::register_inputs(&words("6"), good, ranged).unwrap();
        crate::store::release(query).unwrap();
        crate::store::release(subject).unwrap();
        assert_eq!(shown_ids(&response), ["s1", "ranged"]);
        // A dry run registers its inputs one at a time, as any run without the option.
        let dry = [&words("6")[..], &["-dryrun"]].concat();
        assert!(!inputs_together(&dry));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-152
    // A reader's error on a query record ends the run in the query batch that reads the record,
    // after the batches before, as the CLI's batches do (line WQA's D10, 2026-11-19; before,
    // `register` refused the input and `register_inputs` failed before the run, with the reader's
    // error where the CLI ends with an error of an earlier batch): a line that is not plausible data
    // (`register`), and under `-parse_deflines` an id of the other molecule (`register_inputs`, whose
    // response lists the scan's records past the error without a `shown_id`). BLASTP reads 10000
    // letters a batch: 10 records of 1100 fill the first. A failed run emits no report (the host
    // discards it); the library's tests compare the reports before the error with the CLI's
    // (`LOSAT/src/api/reader_end_tests.rs`).
    #[test]
    fn a_query_reader_error_ends_the_run_in_its_batch() {
        let protein: Vec<u8> = (0..1100u32)
            .map(|i| b"ACDEFGHIKLMNPQRSTVWY"[(i * 7 + i / 20) as usize % 20])
            .collect();
        let fasta = |names: &[&str]| {
            let mut text = Vec::new();
            for name in names {
                text.extend_from_slice(format!(">{name}\n").as_bytes());
                for line in protein.chunks(60) {
                    text.extend_from_slice(line);
                    text.push(b'\n');
                }
            }
            text
        };
        let names: Vec<String> = (1..=11).map(|n| format!("q{n}")).collect();
        let names: Vec<&str> = names.iter().map(String::as_str).collect();
        let subject = fasta(&["s1"]);
        let words = ["blastp", "-query", "q", "-subject", "s", "-outfmt", "6"];

        let mut bad = fasta(&names);
        bad.extend_from_slice(b">bad\nMKVL%%%%\n");
        let line = bad.iter().filter(|&&b| b == b'\n').count();
        let cli = ncbi_in_losat_words(&format!("BLAST query error: CFastaReader: Near line {line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."));
        let (handle, response) = crate::store::register("blastp", 0, true, &bad).unwrap();
        crate::store::release(handle).unwrap();
        assert!(
            response.contains(r#"{"index":10,"id":"q11","length":1100}],"reader_error":""#),
            "{response}"
        );
        assert!(response.ends_with(&format!(r#""reader_error":"{cli}"}}"#)));
        let (result, _) = registered_run(&words, &bad, &subject).unwrap();
        assert_eq!(result.unwrap_err(), format!("{cli}\n"));

        // `-parse_deflines`: the record `ref|NM_000001.1|` names a nucleotide sequence, in batch 2.
        let parse = [&words[..], &["-parse_deflines"]].concat();
        let mut other = fasta(&names);
        other.extend_from_slice(&fasta(&["ref|NM_000001.1|", "q12"]));
        let (result, _) = registered_run(&parse, &other, &subject).unwrap();
        let error = result.unwrap_err();
        assert!(
            error.contains("names the other molecule under -parse_deflines"),
            "{error}"
        );
        let (query, subject_handle, response) =
            crate::store::register_inputs(&parse, &other, &subject).unwrap();
        crate::store::release(query).unwrap();
        crate::store::release(subject_handle).unwrap();
        let split = response.find(",\"subject\":").unwrap();
        let listed = &response[..split];
        assert!(
            listed.contains(r#"{"index":11,"id":"ref|NM_000001.1|","length":1100,"shown_id":null},{"index":12,"id":"q12","length":1100,"shown_id":null}],"reader_error":""#),
            "{listed}"
        );
        assert!(
            listed.contains(&format!(
                r#""reader_error":"{}""#,
                error.trim_end().replace('"', "\\\"")
            )),
            "{listed}"
        );
        // Two records of one id in batch 1 end the run there, before the reading of batch 2 meets
        // the record of the other molecule, as on the CLI.
        let mut dup_names = names.clone();
        dup_names[1] = "q1";
        let mut dup = fasta(&dup_names);
        dup.extend_from_slice(&fasta(&["ref|NM_000001.1|"]));
        let (result, _) = registered_run(&parse, &dup, &subject).unwrap();
        let error = result.unwrap_err();
        assert!(
            error.starts_with("two query records have the Seq-id lcl|q1 under -parse_deflines"),
            "{error}"
        );
    }

    /// The repository's `LOSAT/` directory, where the paths of the `-db` fixtures start
    /// (`LOSAT/tests/db_search_manifest`).
    fn losat_dir() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../LOSAT")
    }

    /// The streams of a run of `words`, whose paths name files of `LOSAT/` (the run's working
    /// directory, as the host's preopened `.` gives it, docs/web/abi_v2.md §3), with the query of
    /// `register` and no subject (`-db`, line WQA D12a).
    fn database_streams(
        words: &[&str],
        query: &str,
    ) -> Result<std::collections::BTreeMap<u32, Vec<u8>>, String> {
        let dir = losat_dir();
        let absolute = |word: &str| {
            if word.starts_with("tests/") {
                dir.join(word).to_string_lossy().into_owned()
            } else {
                word.to_string()
            }
        };
        let words: Vec<String> = words.iter().map(|word| absolute(word)).collect();
        let words: Vec<&str> = words.iter().map(String::as_str).collect();
        let bytes = std::fs::read(dir.join(query)).unwrap();
        // The query is read with the data loaders of the argv, as the host reads it (abi_v2.md §9).
        let loaders = data_loaders(&words)?;
        let (handle, _) = crate::store::register(words[0], 0, loaders, &bytes)?;
        crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
        let program = Program::parse(words[0]).unwrap();
        let result = crate::store::with_inputs(
            program,
            handle,
            0,
            &words,
            Some(loaders),
            |queries, end, subjects| {
                assert!(subjects.is_empty());
                run(&words, queries, end, subjects)
            },
        );
        crate::store::release(handle).unwrap();
        result??;
        Ok(emitted_streams())
    }

    // A search of `-db` (line WQA D12a, 2026-11-20): `validate` takes it, `run` searches the
    // database with the library's database entry (`run_local_database_<program>`), as the command
    // line does, and each format is NCBI BLAST+ 2.17.0's frozen bytes of the same argv in LOSAT's
    // words (LOSAT/tests/fixtures/db_search, the rows of LOSAT/tests/db_search_manifest), the
    // filters too. The run takes no subject handle: `with_inputs` refuses one, and
    // `register_inputs` reads the query alone.
    #[test]
    fn a_database_search_runs_as_the_command_line() {
        let frozen = |row: &str| {
            ncbi_output_in_losat_words(
                &std::fs::read(losat_dir().join(format!("tests/fixtures/db_search/{row}.out")))
                    .unwrap(),
            )
        };
        for (row, words, query) in [
            (
                "db.d1.6",
                &["blastp", "-db", "tests/db_search/paja", "-query", "q"][..],
                "tests/db_search/d1_q.faa",
            ),
            (
                "db.blastn.mega.6",
                &["blastn", "-db", "tests/db_search/nslice", "-query", "q"],
                "tests/db_search/nslice_q.fna",
            ),
            (
                "db.sl.bp.bsl.6",
                &[
                    "blastp",
                    "-db",
                    "tests/db_search/taxfilter/tfp",
                    "-query",
                    "q",
                    "-seqidlist",
                    "tests/db_search/seqidlist/sl_tfp.bsl",
                ],
                "tests/db_search/taxfilter/q.faa",
            ),
            (
                "db.tf.bp.nosql.6",
                &[
                    "blastp",
                    "-db",
                    "tests/db_search/taxfilter/tfp",
                    "-query",
                    "q",
                    "-taxids",
                    "40674",
                ],
                "tests/db_search/taxfilter/q.faa",
            ),
            (
                "db.mb.bp.mb1.6",
                &[
                    "blastp",
                    "-db",
                    "tests/seqdb/memb_bit/db/mb_1",
                    "-query",
                    "q",
                ],
                "tests/seqdb/memb_bit/in/qp.fa",
            ),
            (
                "db.mb.bn.mb1.6",
                &[
                    "blastn",
                    "-db",
                    "tests/seqdb/memb_bit/db/nmb_gm_1",
                    "-query",
                    "q",
                    "-dust",
                    "no",
                ],
                "tests/seqdb/memb_bit/in/qn.fa",
            ),
        ] {
            validate(words).unwrap_or_else(|error| panic!("{row}: {error}"));
            assert!(searches_database(words), "{row}");
            assert!(start_json(words).unwrap().ends_with(",\"database\":true}"));
            let streams = database_streams(words, query).unwrap_or_else(|e| panic!("{row}: {e}"));
            assert_eq!(text(&streams[&6]), text(&frozen(row)), "{row}");
        }
        // A run of -db takes no subject handle.
        let words = ["blastp", "-db", "d", "-query", "q"];
        // Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1).
        let protein = b">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ\n";
        let (query, _) = crate::store::register("blastp", 0, true, protein).unwrap();
        let (subject, _) = crate::store::register("blastp", 1, true, protein).unwrap();
        let refused = crate::store::with_inputs(
            Program::Blastp,
            query,
            subject,
            &words,
            Some(true),
            |_, _, _| (),
        );
        assert_eq!(
            refused.unwrap_err(),
            format!("a run of -db takes no subject handle: give 0, not {subject}")
        );
        crate::store::release(subject).unwrap();
        // Its query handle is checked as any run's: 0 is no query (review of D12a, low 1: the
        // check skipped every handle 0 of a -db run, and the run stopped the instance).
        let no_query =
            crate::store::with_inputs(Program::Blastp, 0, 0, &words, Some(true), |_, _, _| ());
        assert_eq!(no_query.unwrap_err(), "0 is not a registered query handle");
        // A database that is not there: the CLI's error, from the run.
        let error = crate::store::with_inputs(
            Program::Blastp,
            query,
            0,
            &words,
            Some(true),
            |queries, end, _| run(&words, queries, end, &[]),
        )
        .unwrap()
        .unwrap_err();
        assert!(
            error.starts_with(&ncbi_in_losat_words(
                "BLAST Database error: No alias or index file found for protein database [d]"
            )),
            "{error}"
        );
        crate::store::release(query).unwrap();
        // -parse_deflines: `register_inputs` reads the query alone, with the ids of its deflines;
        // the subject handle is 0 and the JSON's subject is null.
        let parsed = ["blastp", "-db", "d", "-query", "q", "-parse_deflines"];
        assert!(inputs_together(&parsed));
        let (query, subject, response) =
            crate::store::register_inputs(&parsed, b">sp|P1|X q\nMKTAYIAKQ\n", b"").unwrap();
        assert_eq!(subject, 0);
        assert!(response.ends_with(",\"subject\":null}"), "{response}");
        crate::store::with_inputs(
            Program::Blastp,
            query,
            0,
            &parsed,
            Some(true),
            |queries, _, subjects| {
                assert_eq!(queries.len(), 1);
                assert!(subjects.is_empty());
            },
        )
        .unwrap();
        crate::store::release(query).unwrap();
        assert!(crate::store::register_inputs(&parsed, b">q\nMK\n", b">s\nMK\n").is_err());
        // The query of a -db run is read with the data loaders of its registry, as every run's
        // inputs are (line WQA's D13 and D16 with D12a, 2026-11-21): with `[BLAST]
        // DATA_LOADERS=none` of `-conffile`, a first line that may be a Seq-id is data, and a run
        // with the data loaders on takes no query read with them off.
        let registry = b"[BLAST]\nDATA_LOADERS=none\n";
        let file = crate::files::begin(crate::files::KIND_REGISTRY, "d12a-db.ini", registry.len())
            .unwrap();
        crate::files::chunk(file, registry).unwrap();
        crate::files::end(file).unwrap();
        let off = [
            "blastp",
            "-conffile",
            "d12a-db.ini",
            "-db",
            "d",
            "-query",
            "q",
        ];
        let on = [&off[..1], &off[3..]].concat();
        assert_eq!(data_loaders(&off), Ok(false));
        assert_eq!(data_loaders(&on), Ok(true));
        let seq_id = b"2HHB\nMKTAYIAKQR\n";
        // The host registers the query with the data loaders of each argv's `validate`.
        let (query, response) = crate::store::register("blastp", 0, false, seq_id).unwrap();
        assert!(response.contains(r#""length":13"#), "{response}");
        let error = crate::store::register("blastp", 0, true, seq_id).unwrap_err();
        assert!(
            error.contains("the first line of the query (\"2HHB\")"),
            "{error}"
        );
        let refused =
            crate::store::with_inputs(Program::Blastp, query, 0, &on, Some(true), |_, _, _| ());
        assert!(
            refused
                .unwrap_err()
                .contains("read with the data loaders off"),
            "a -db run with the data loaders on"
        );
        assert!(crate::store::with_inputs(
            Program::Blastp,
            query,
            0,
            &off,
            Some(false),
            |_, _, _| ()
        )
        .is_ok());
        crate::store::release(query).unwrap();
        // With -parse_deflines the reader of the ids wants a defline, with the data loaders off
        // too. `register_inputs` keeps the reader's error with the query (`reader_error`, line
        // WQA's D10), and the run meets it in its first query batch, after it opens the database,
        // as the CLI does (2026-11-22): a database that is not there gives the CLI's database
        // error, and the db_search database `paja` the CLI's query error.
        let parsed = [&off[..], &["-parse_deflines"]].concat();
        let cli = ncbi_in_losat_words(
            "BLAST query error: CFastaReader: Input doesn't start with a defline or comment around line 1",
        );
        let paja = losat_dir()
            .join("tests/db_search/paja")
            .to_string_lossy()
            .into_owned();
        let in_paja: Vec<&str> = parsed
            .iter()
            .map(|&word| if word == "d" { paja.as_str() } else { word })
            .collect();
        for words in [&parsed[..], &in_paja[..]] {
            let (query, subject, response) =
                crate::store::register_inputs(words, seq_id, b"").unwrap();
            assert_eq!(subject, 0);
            assert!(
                response.ends_with(&format!(r#","reader_error":"{cli}"}},"subject":null}}"#)),
                "{response}"
            );
            let error = crate::store::with_inputs(
                Program::Blastp,
                query,
                0,
                words,
                Some(false),
                |queries, end, _| run(words, queries, end, &[]),
            )
            .unwrap()
            .unwrap_err();
            crate::store::release(query).unwrap();
            if words == parsed {
                assert!(
                    error.starts_with(&ncbi_in_losat_words(
                        "BLAST Database error: No alias or index file found for protein database [d]"
                    )),
                    "{error}"
                );
            } else {
                assert_eq!(error, format!("{cli}\n"));
            }
        }
        crate::files::release(file).unwrap();
        // BLASTN's database of repeats is accepted as on the command line (line RG item D4).
        let repeats = [
            "blastn",
            "-query",
            "q",
            "-subject",
            "s",
            "-filtering_db",
            "repeats",
        ];
        validate(&repeats).unwrap();
        // A repeat database whose files are absent: NCBI posts the missing-database error as a
        // batch warning after the engine prefix (as on the CLI; under WebAssembly its search path
        // names `/work::`, abi_v2.md §3).
        let dna = b">q\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
        let (q_nuc, _) = crate::store::register("blastn", 0, true, dna).unwrap();
        let (s_nuc, _) = crate::store::register("blastn", 1, true, dna).unwrap();
        let missing_repeats = crate::store::with_inputs(
            Program::Blastn,
            q_nuc,
            s_nuc,
            &repeats,
            Some(true),
            |queries, end, subjects| run(&repeats, queries, end, subjects),
        )
        .unwrap()
        .unwrap_err();
        assert!(
            missing_repeats.starts_with(&ncbi_in_losat_words(
                "BLAST engine error: Warning: No alias or index file found for nucleotide database [repeats]"
            )),
            "{missing_repeats}"
        );
        crate::store::release(q_nuc).unwrap();
        crate::store::release(s_nuc).unwrap();
        // `-db` with `-subject` is NCBI's reader's exclusion, as on the command line.
        let both = ["blastp", "-query", "q", "-subject", "s", "-db", "swissprot"];
        assert!(validate(&both)
            .unwrap_err()
            .contains("Incompatible with argument:  `"));
        assert!(!searches_database(&both));
        // `-dbsize` is an option of the search: accepted with the subject records.
        validate(&[
            "blastp", "-query", "q", "-subject", "s", "-dbsize", "1000000",
        ])
        .unwrap();
    }

    // The HSP records give the registered index of a query where `-query_loc` skips
    // an earlier record (its interval starts past its end), and records without residues
    // are registered records too.
    #[test]
    fn q_idx_is_the_registered_index() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [
            FastaRecord::new("Query_1", b"short", b"MK"),
            FastaRecord::new("Query_2", b"empty", b""),
            FastaRecord::new("Query_3", b"q", protein),
        ];
        let subject = [FastaRecord::new("Subject_1", b"s", protein)];
        let q_idx = |words: &[&str]| {
            crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
            run(words, &query[1..], &losat::api::QueryEnd::Input, &subject).unwrap();
            let hits = crate::emit::EMITTED.with(|emitted| {
                emitted
                    .borrow()
                    .iter()
                    .filter(|(stream, _)| *stream == crate::emit::STREAM_HITS)
                    .flat_map(|(_, bytes)| bytes.clone())
                    .collect::<Vec<u8>>()
            });
            let hits = String::from_utf8(hits).unwrap();
            assert!(!hits.is_empty(), "{words:?}");
            hits.lines()
                .map(|line| {
                    let at = line.find("\"q_idx\":").unwrap() + 8;
                    line[at..]
                        .split(',')
                        .next()
                        .unwrap()
                        .parse::<usize>()
                        .unwrap()
                })
                .collect::<Vec<_>>()
        };
        // The query without residues is searched as a query without data (engine q_idx 1);
        // with the range it is skipped (engine q_idx 0), and the record stays the second.
        let words = ["blastp", "-query", "q", "-subject", "s"];
        assert!(q_idx(&words).iter().all(|&index| index == 1));
        let ranged = [&words[..], &["-query_loc", "5-40"]].concat();
        assert!(q_idx(&ranged).iter().all(|&index| index == 1));
        let ordinals = searched_queries(Some(std::ffi::OsStr::new("5-40")), &query).unwrap();
        assert_eq!(ordinals, vec![2]);
        assert_eq!(searched_queries(None, &query), None);
    }

    fn fasta(text: &[u8]) -> Vec<FastaRecord> {
        let text = std::str::from_utf8(text).unwrap();
        text.split('>')
            .filter(|entry| !entry.is_empty())
            .map(|entry| {
                let (name, residues) = entry.split_once('\n').unwrap();
                FastaRecord::new(name, name.as_bytes(), residues.trim_end().as_bytes())
            })
            .collect()
    }

    /// Bytes as text, so that a failed comparison prints readable text.
    fn text(bytes: &[u8]) -> String {
        String::from_utf8_lossy(bytes).into_owned()
    }

    /// The bytes of each stream that `run` emitted, joined.
    /// The streams of a run as a host keeps them (docs/web/abi_v2.md §5): each stream's chunks in
    /// order, and for each signal of stream 9, what came of the stream it names dropped (stream 9
    /// itself holds the signals).
    fn streams(
        words: &[&str],
        query: &[FastaRecord],
        subject: &[FastaRecord],
    ) -> Result<std::collections::BTreeMap<u32, Vec<u8>>, String> {
        crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
        run(words, query, &losat::api::QueryEnd::Input, subject)?;
        Ok(emitted_streams())
    }

    /// The streams emitted since `EMITTED` was cleared, as a host keeps them (`streams`).
    fn emitted_streams() -> std::collections::BTreeMap<u32, Vec<u8>> {
        let mut streams = std::collections::BTreeMap::new();
        crate::emit::EMITTED.with(|emitted| {
            for (stream, bytes) in emitted.borrow().iter() {
                if *stream == crate::emit::STREAM_DISCARD {
                    let discarded = u32::from_le_bytes(bytes[..4].try_into().unwrap());
                    streams.insert(discarded, Vec::new());
                }
                streams
                    .entry(*stream)
                    .or_insert_with(Vec::new)
                    .extend_from_slice(bytes);
            }
        });
        streams
    }

    /// The streams that the signals of stream 9 of `streams` discard, in order. Each signal says
    /// that the error of its stream is on stream 3: these runs have no `-logfile`, so their
    /// messages go to stream 3 (`run_start_tests` has the places of `-logfile`).
    fn discarded(streams: &std::collections::BTreeMap<u32, Vec<u8>>) -> Vec<u32> {
        let signals = streams
            .get(&crate::emit::STREAM_DISCARD)
            .map_or(&[][..], Vec::as_slice);
        assert_eq!(signals.len() % 8, 0, "signals of 8 bytes");
        signals
            .chunks(8)
            .map(|signal| {
                let place = u32::from_le_bytes(signal[4..].try_into().unwrap());
                assert_eq!(place, STREAM_DIAGNOSTICS, "the place of the error");
                u32::from_le_bytes(signal[..4].try_into().unwrap())
            })
            .collect()
    }

    // Line WQA's sweep, D15 (2026-11-20): a run that fails sends the bytes that the CLI's run of the argv leaves in
    // its output before the error (NCBI's `The average subject length is too short` comes after the outfmt 0 prolog:
    // an empty file in the other formats), and drops every other output with the error that ended the run as the
    // place of its error; stream 3 is the warnings before the error.
    #[test]
    fn a_failed_run_sends_the_output_of_its_argv_and_drops_the_others() {
        let query = fasta(
            b">q\nACGTACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGACTGACTAGCTAGCTAGCATGCATCGATCGATC\n",
        );
        let subject = fasta(b">s1\n>s2\n");
        for (outfmt, kept) in [
            (None, 0),
            (Some("0"), 0),
            (Some("6"), 6),
            (Some("7"), 7),
            (Some("6 qseqid"), crate::emit::STREAM_CUSTOM),
            (Some("5"), crate::emit::STREAM_CUSTOM),
        ] {
            let mut words = vec!["blastn", "-query", "q", "-subject", "s"];
            if let Some(outfmt) = outfmt {
                words.extend(["-outfmt", outfmt]);
            }
            crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
            let error = run(&words, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
            assert!(
                error.contains("The average subject length is too short"),
                "{outfmt:?}: {error}"
            );
            let emitted = crate::emit::EMITTED.with(|emitted| emitted.borrow().clone());
            let mut sent = std::collections::BTreeMap::<u32, Vec<u8>>::new();
            let mut dropped = Vec::new();
            for (stream, bytes) in &emitted {
                if *stream == crate::emit::STREAM_DISCARD {
                    let place = u32::from_le_bytes(bytes[4..].try_into().unwrap());
                    assert_eq!(
                        place,
                        crate::emit::PLACE_RUN,
                        "{outfmt:?}: the place of the error"
                    );
                    dropped.push(u32::from_le_bytes(bytes[..4].try_into().unwrap()));
                } else {
                    sent.entry(*stream).or_default().extend_from_slice(bytes);
                }
            }
            dropped.sort_unstable();
            // The program's formats besides the argv's own, which a custom -outfmt adds to.
            let others: Vec<u32> = [0, 6, 7]
                .into_iter()
                .filter(|&stream| stream != kept)
                .collect();
            assert_eq!(dropped, others, "{outfmt:?}: the outputs dropped");
            for stream in &dropped {
                assert!(
                    !sent.contains_key(stream),
                    "{outfmt:?}: nothing of dropped stream {stream}"
                );
            }
            // The outfmt 0 prolog is the file that the CLI's run of outfmt 0 leaves; the others leave none.
            let output = sent.get(&kept).map(|bytes| text(bytes));
            if kept == 0 {
                let prolog = output.expect("the prolog of outfmt 0");
                assert!(
                    prolog.starts_with(&format!("BLASTN {}", losat::web_adapter::REPORT_VERSION)),
                    "{outfmt:?}: {prolog}"
                );
                assert!(
                    prolog.contains("2 sequences; 0 total letters"),
                    "{outfmt:?}: {prolog}"
                );
            } else {
                assert_eq!(output, None, "{outfmt:?}: an empty output is no chunk");
            }
            assert!(
                text(&sent[&crate::emit::STREAM_DIAGNOSTICS])
                    .contains("Subject sequence contains no data"),
                "{outfmt:?}: stream 3"
            );
        }
    }

    // The searches of a run: one per distinct hit list size of the formats 0, 6 and 7 (NCBI
    // `blast_args.cpp:2894-2973`), the tabular formats first; without -num_descriptions and
    // -num_alignments the one search of the program's formats, the tabular ones first.
    #[test]
    fn a_run_searches_once_for_each_hit_list_size() {
        let groups = |program: &str, options: &[&str]| {
            let words = [&[program, "-query", "q", "-subject", "s"][..], options].concat();
            let (program, command, requested) = parse_run(&words).unwrap();
            search_groups(&outputs(program, &requested), &command)
                .iter()
                .map(|group| group.iter().map(|output| output.stream).collect())
                .collect::<Vec<Vec<u32>>>()
        };
        // The argv's own format first in its search (line WQA's sweep, D3 review r2 M1,
        // 2026-11-19), then the tabular formats; the tabular formats' search first.
        let one = vec![vec![0, 6, 7]];
        let one_tabular_first = vec![vec![0, 6, 7]];
        let two = vec![vec![6, 7], vec![0]];
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            for (options, expected) in [
                (&[][..], &one),
                (&["-max_target_seqs", "3"], &one),
                // The hit list of format 0 is max(-num_descriptions, -num_alignments), with
                // 500 and 250 for the one that is not given; the tabular formats' is
                // -num_alignments, or 500.
                (&["-num_descriptions", "500"], &one_tabular_first),
                (&["-num_descriptions", "3"], &two),
                (&["-num_alignments", "600"], &one_tabular_first),
                (&["-num_alignments", "3"], &two),
                (&["-num_alignments", "0"], &two),
                (
                    &["-num_descriptions", "2", "-num_alignments", "4"],
                    &one_tabular_first,
                ),
                (&["-num_descriptions", "4", "-num_alignments", "2"], &two),
                (
                    &["-num_descriptions", "0", "-num_alignments", "0"],
                    &one_tabular_first,
                ),
            ] {
                assert_eq!(&groups(program, options), expected, "{program} {options:?}");
            }
            // The output of the argv's -outfmt (stream 4) has the hit list size of its number's
            // formats: it joins their search, first, and adds none.
            for (options, expected) in [
                (&["-outfmt", "6 qseqid"][..], vec![vec![4, 6, 7, 0]]),
                (
                    &["-outfmt", "7 std qcovs", "-num_alignments", "3"],
                    vec![vec![4, 6, 7], vec![0]],
                ),
                (
                    &["-outfmt=6 sseqid", "-num_descriptions", "500"],
                    vec![vec![4, 6, 7, 0]],
                ),
                (&["-outfmt", "7"], vec![vec![7, 6, 0]]),
                (&["-outfmt", "6"], vec![vec![6, 7, 0]]),
                (&["-outfmt", "0"], vec![vec![0, 6, 7]]),
                (&["-outfmt", "1"], vec![vec![4, 6, 7, 0]]),
                (
                    &["-outfmt", "2", "-num_alignments", "3"],
                    vec![vec![6, 7], vec![4, 0]],
                ),
            ] {
                assert_eq!(groups(program, options), expected, "{program} {options:?}");
            }
        }
    }

    // -outfmt is read as the CLI reads it, with the CLI's messages; 0, 6 and 7 add no output, 6 and
    // 7 with words after the number add one (WA2-OP-03, Owner 2026-11-05), and so does every other
    // format that the CLI writes (W11), 13 and 14 without the CLI's need of -out.
    #[test]
    fn outfmt_is_read_as_on_the_command_line() {
        let custom = |program: &str, outfmt: &str| {
            parse_run(&[program, "-query", "q", "-subject", "s", "-outfmt", outfmt])
                .map(|(_, _, requested)| requested.custom)
        };
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            for outfmt in ["0", "6", "7", " 7 ", "0 qseqid", "+6"] {
                assert_eq!(custom(program, outfmt), Ok(None), "{program} {outfmt:?}");
            }
            for outfmt in [
                "6 qseqid sseqid pident",
                "7 std qcovs",
                "6 std",
                "6 delim=,",
                " 6 -qseqid ",
                "+7 sseqid",
            ] {
                assert_eq!(
                    custom(program, outfmt),
                    Ok(Some(outfmt.to_string())),
                    "{program} {outfmt:?}"
                );
            }
            // The -flag=value form, and the arguments come back with the CLI's default.
            let words = [program, "-query", "q", "-subject", "s", "-outfmt=6 qseqid"];
            let (_, mut command, given) = parse_run(&words).unwrap();
            assert_eq!(given.custom.as_deref(), Some("6 qseqid"));
            assert_eq!(outfmt_of(&mut command), "0");
            // Every other format that the CLI writes adds an output, with or without words (W11).
            for outfmt in [
                "5",
                "8",
                "9",
                "10",
                "12",
                "15",
                "16",
                "20",
                " 5 ",
                "+010",
                "10 qseqid",
                "20 delim=| std",
                "5 qseqid",
                "11",
                "13",
                "14",
                "14 qseqid",
            ] {
                assert_eq!(
                    custom(program, outfmt),
                    Ok(Some(outfmt.to_string())),
                    "{program} {outfmt:?}"
                );
                let words = [program, "-query", "q", "-subject", "s", "-outfmt", outfmt];
                assert_eq!(validate(&words), Ok(true), "{program} {outfmt:?}");
            }
            // Every program writes the query-anchored formats 1 to 4 (line RG T3: BLASTP in T3a,
            // BLASTN and TBLASTN in T3b, BLASTX and TBLASTX in T3c), each an output of its own
            // without the words after the number, which the CLI drops.
            let anchored: &[(&str, &str)] = &[
                ("1", "1"),
                ("2", "2"),
                ("3", "3"),
                ("4", "4"),
                (" 3 ", "3"),
                ("+02", "2"),
                ("1 qseqid", "1"),
                ("4 std", "4"),
            ];
            for &(outfmt, number) in anchored {
                assert_eq!(
                    custom(program, outfmt),
                    Ok(Some(number.to_string())),
                    "{program} {outfmt:?}"
                );
                let words = [program, "-query", "q", "-subject", "s", "-outfmt", outfmt];
                assert_eq!(validate(&words), Ok(true), "{program} {outfmt:?}");
            }
            // What the CLI rejects, with its message.
            let engine = engine_program(Program::parse(program).unwrap());
            let rejected: &[&str] = &["abc", "99", "6\tqseqid"];
            for &outfmt in rejected {
                let cli = losat::web_adapter::parse_formatting_string(outfmt)
                    .and_then(|choice| {
                        losat::web_adapter::report_format(&choice, engine, None)?;
                        losat::web_adapter::formatting_handler_check(&choice, engine)
                    })
                    .unwrap_err();
                assert_eq!(
                    custom(program, outfmt),
                    Err(format!("{cli:#}")),
                    "{program} {outfmt:?}"
                );
            }
            // The formatting handler's check of a custom delimiter comes where the CLI makes it:
            // after the ranges (and BLASTN's query options), in validate and in the engine's run.
            let conflict = "6 delim=; qseqid staxids";
            assert_eq!(custom(program, conflict), Ok(Some(conflict.to_string())));
            let words = [program, "-query", "q", "-subject", "s", "-outfmt", conflict];
            let handler = losat::web_adapter::parse_formatting_string(conflict)
                .and_then(|choice| losat::web_adapter::formatting_handler_check(&choice, engine))
                .unwrap_err();
            assert_eq!(validate(&words), Err(format!("{handler:#}")), "{program}");
            assert!(validate(&words)
                .unwrap_err()
                .contains("Your custom record separator (;)"));
            let ranged = [&words[..], &["-subject_loc", "5-2"]].concat();
            assert!(
                validate(&ranged)
                    .unwrap_err()
                    .starts_with(&ncbi_in_losat_words(
                        "BLAST engine error: Invalid specification of subject location"
                    )),
                "{program}: {:?}",
                validate(&ranged)
            );
            let query = [FastaRecord::new("q", b"q", b"ACGTACGTACGT")];
            assert_eq!(
                run(&words, &query, &losat::api::QueryEnd::Input, &query),
                Err(format!("{handler:#}")),
                "{program}"
            );
            // The CLI rejects the argv whatever the hit list sizes: a run of two searches does not
            // write the other one, and validate gives the error of the CLI's checks with the argv's
            // -outfmt (review 2).
            for extra in [
                &["-num_descriptions", "3"][..],
                &["-num_alignments", "3"],
                &["-num_alignments", "0"],
                &["-evalue", "0"],
            ] {
                let both = [&words[..], extra].concat();
                assert_eq!(
                    validate(&both),
                    Err(format!("{handler:#}")),
                    "{program} {extra:?}"
                );
                assert_eq!(
                    run(&both, &query, &losat::api::QueryEnd::Input, &query),
                    Err(format!("{handler:#}")),
                    "{program} {extra:?}"
                );
            }
            let twice = [
                program, "-query", "q", "-subject", "s", "-outfmt", "6", "-outfmt", "7",
            ];
            assert!(parse(&twice)
                .unwrap_err()
                .contains("Argument with this name is defined already: outfmt"));
            assert!(validate(&[
                program,
                "-query",
                "q",
                "-subject",
                "s",
                "-outfmt",
                "7 std qcovs"
            ])
            .is_ok());
        }
        assert_eq!(
            custom("blastn", "abc").unwrap_err(),
            ncbi_in_losat_words("BLAST query/options error: 'abc' is not a valid output format\nPlease refer to the BLAST+ user manual.\n")
        );
        assert_eq!(
            custom("blastn", "99").unwrap_err(),
            "Error: Formatting choice is out of range\n"
        );
        assert_eq!(custom("blastp", "17").unwrap_err(), ncbi_in_losat_words("BLAST query/options error: SAM format is only applicable to blastn\nPlease refer to the BLAST+ user manual.\n"));
        assert_eq!(custom("blastn", "17 SQ"), Ok(Some("17 SQ".to_string())));
    }

    /// The streams of `run` over inputs that the program's reader reads, as a host registers them
    /// for the argv (`store::register`, or `store::register_inputs` where `inputs_together`).
    fn registered_streams(
        words: &[&str],
        query: &[u8],
        subject: &[u8],
    ) -> std::collections::BTreeMap<u32, Vec<u8>> {
        try_registered_streams(words, query, subject).unwrap()
    }

    /// `registered_streams`, or the error of the registration or the run.
    fn try_registered_streams(
        words: &[&str],
        query: &[u8],
        subject: &[u8],
    ) -> Result<std::collections::BTreeMap<u32, Vec<u8>>, String> {
        let (result, streams) = registered_run(words, query, subject)?;
        result.map(|()| streams)
    }

    /// `try_registered_streams` that keeps the streams of a run that fails: what the run emitted
    /// before its error. `Err` when the inputs are not registered.
    #[allow(clippy::type_complexity)]
    fn registered_run(
        words: &[&str],
        query: &[u8],
        subject: &[u8],
    ) -> Result<(Result<(), String>, std::collections::BTreeMap<u32, Vec<u8>>), String> {
        let program = Program::parse(words[0]).unwrap();
        let (query, subject) = if inputs_together(words) {
            let (query, subject, _) = crate::store::register_inputs(words, query, subject)?;
            (query, subject)
        } else {
            let (query, _) = crate::store::register(words[0], 0, true, query)?;
            let (subject, _) =
                crate::store::register(words[0], 1, true, subject).inspect_err(|_| {
                    crate::store::release(query).unwrap();
                })?;
            (query, subject)
        };
        crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
        let result = crate::store::with_inputs(
            program,
            query,
            subject,
            words,
            data_loaders(words).ok(),
            |queries, end, subjects| run(words, queries, end, subjects),
        );
        crate::store::release(query).unwrap();
        crate::store::release(subject).unwrap();
        let mut streams = std::collections::BTreeMap::new();
        crate::emit::EMITTED.with(|emitted| {
            for (stream, bytes) in emitted.borrow().iter() {
                streams
                    .entry(*stream)
                    .or_insert_with(Vec::new)
                    .extend_from_slice(bytes);
            }
        });
        Ok((result?, streams))
    }

    /// NCBI BLAST+ 2.17.0's frozen outputs of `-outfmt` with fields, for every program
    /// (`LOSAT/tests/outfmt0_manifest.tsv`, rows `tab.*`: id, program, query, subject, -outfmt).
    const CUSTOM_ROWS: [(&str, &str, &str, &str, &str); 12] = [
        ("tab.bn.delim_pipe.6", "blastn", "tests/db_search/nslice_q.fna", "tests/db_search/nslice.fna", "6 delim=| std sstrand"),
        ("tab.bn.align.7", "blastn", "tests/db_search/nslice_q.fna", "tests/db_search/nslice.fna", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        ("tab.bp.cov.6", "blastp", "tests/db_search/ids/ids_q.faa", "tests/fasta/outfmt11/arch_bp_s.faa", "6 std qcovs qcovhsp qcovus"),
        ("tab.bp.ids.7", "blastp", "tests/db_search/ids/ids_q.faa", "tests/db_search/ids/prot.faa", "7 std qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen"),
        ("tab.bx.qframe.6", "blastx", "tests/db_search/ids/ids_q.fna", "tests/fasta/outfmt11/arch_bp_s.faa", "6 qseqid qstart frames sframe qframe"),
        ("tab.bx.align.7", "blastx", "tests/db_search/ids/ids_q.fna", "tests/fasta/outfmt11/arch_bp_s.faa", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        ("tab.tbn.delim_pipe.6", "tblastn", "tests/db_search/tab_tbn_q.faa", "tests/db_search/tslice.fna", "6 delim=| std sframe sstrand"),
        ("tab.tbn.align.7", "tblastn", "tests/db_search/tab_tbn_q.faa", "tests/db_search/tslice.fna", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        ("tab.tbx.removed.6", "tblastx", "tests/db_search/tab_tbx_q.fna", "tests/db_search/xslice.fna", "6 -qseqid"),
        ("tab.tbx.align.7", "tblastx", "tests/db_search/tab_tbx_q.fna", "tests/db_search/xslice.fna", "7 qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"),
        // Fields whose writer warns (taxdb, once per query block).
        ("tab.bn.tax.7", "blastn", "tests/db_search/nslice_q.fna", "tests/db_search/nslice.fna", "7 staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp qcovus"),
        ("tab.bx.delim_comma.6", "blastx", "tests/db_search/ids/ids_q.fna", "tests/fasta/outfmt11/arch_bp_s.faa", "6 delim=, qseqid sseqid scomnames"),
    ];

    // The output of the argv's -outfmt (6 or 7 with fields here) is stream 4, NCBI's bytes of the same
    // command line, in the five programs, and stream 3 NCBI's standard error of it (with the
    // warnings of its fields' writer); the formats 0, 6 and 7 and the HSP records are those of
    // the run without -outfmt.
    #[test]
    fn a_custom_outfmt_is_written_as_the_cli_writes_it() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        for (id, program, query, subject, outfmt) in CUSTOM_ROWS {
            let words = [program, "-query", query, "-subject", subject];
            let plain = registered_streams(&words, &read(query), &read(subject));
            let with = [&words[..], &["-outfmt", outfmt]].concat();
            let custom = registered_streams(&with, &read(query), &read(subject));
            let ncbi =
                ncbi_output_in_losat_words(&read(&format!("tests/fixtures/outfmt0/{id}.out")));
            assert_eq!(text(&custom[&STREAM_CUSTOM]), text(&ncbi), "{id}");
            let stderr = ncbi_output_in_losat_words(
                &std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                    .unwrap_or_default(),
            );
            let diagnostics = custom.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(text(&diagnostics), text(&stderr), "{id}: stream 3");
            for stream in [0, 6, 7, STREAM_HITS] {
                assert_eq!(
                    plain.get(&stream),
                    custom.get(&stream),
                    "{id}: stream {stream}"
                );
            }
            assert!(!plain.contains_key(&STREAM_CUSTOM), "{id}");
        }
    }

    /// NCBI BLAST+ 2.17.0's frozen outputs of the formats besides 0, 6 and 7 (rows of
    /// `LOSAT/tests/outfmt0_manifest.tsv`: id and argv, run from LOSAT/): every program in 5, 8, 9,
    /// 12, 15 and 16, 10 and 20 with and without fields, the BLAST archive (11) of four programs
    /// with `-query_loc`, `-subject_loc` and deflines that the title does not give, SAM (17), and
    /// BLASTP's ungapped search in 8 and 11 (Dense-diag Seq-aligns, line NP item O4b) and with an HSP
    /// of length 0 in 5, 8 and 15 (`-window_size 0`, also item O4b), and BLASTN's ungapped search in
    /// SAM (one record per subject, item O4a-2).
    const FORMAT_ROWS: [(&str, &[&str]); 62] = [
        (
            "sm1.blastp.x16.soft.5",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/x16_seg_q.faa",
                "-subject",
                "tests/fasta/outfmt0/x16_seg_s.faa",
                "-seg",
                "yes",
                "-soft_masking",
                "true",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sm1.blastp.x16.soft.11",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/x16_seg_q.faa",
                "-subject",
                "tests/fasta/outfmt0/x16_seg_s.faa",
                "-seg",
                "yes",
                "-soft_masking",
                "true",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sm1.blastp.x16.soft.16",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/x16_seg_q.faa",
                "-subject",
                "tests/fasta/outfmt0/x16_seg_s.faa",
                "-seg",
                "yes",
                "-soft_masking",
                "true",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sm2.tblastx.mask.soft.5",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-soft_masking",
                "true",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sm2.tblastx.mask.soft.16",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-soft_masking",
                "true",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sm3.megablast.mask.hard.5",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-task",
                "megablast",
                "-soft_masking",
                "false",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sm3.megablast.mask.hard.16",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-task",
                "megablast",
                "-soft_masking",
                "false",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sm3.megablast.mask.hard.11",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/mask_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/mask_subject.fasta",
                "-task",
                "megablast",
                "-soft_masking",
                "false",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.multi.megablast.8",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/multi_subject.fasta",
                "-task",
                "megablast",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.multi.megablast.12",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/multi_subject.fasta",
                "-task",
                "megablast",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.8",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.12",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.e2e.tblastn.xdrop_final.8",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-xdrop_gap_final",
                "1e9",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.e2d.tblastn.sloc.12",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-subject",
                "tests/fasta/outfmt0/e2d_n_sub8k.fa",
                "-subject_loc",
                "1001-6000",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.tblastx.nohit.8",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/tblastx_multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/tblastx_multi_subject.fasta",
                "-evalue",
                "1e-100",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.tblastx.nohit.12",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/tblastx_multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/tblastx_multi_subject.fasta",
                "-evalue",
                "1e-100",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.blastx.multi.8",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-outfmt",
                "8",
            ],
        ),
        (
            "sa.blastx.multi.12",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-outfmt",
                "12",
            ],
        ),
        (
            "sa.multi.megablast.9",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/multi_subject.fasta",
                "-task",
                "megablast",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.9",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.e2d.tblastn.sloc.9",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-subject",
                "tests/fasta/outfmt0/e2d_n_sub8k.fa",
                "-subject_loc",
                "1001-6000",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.tblastx.nohit.9",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/tblastx_multi_query.fasta",
                "-subject",
                "tests/fasta/outfmt0/tblastx_multi_subject.fasta",
                "-evalue",
                "1e-100",
                "-outfmt",
                "9",
            ],
        ),
        (
            "sa.blastx.multi.9",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject",
                "tests/fasta/outfmt0/e2d_t_qp2.faa",
                "-outfmt",
                "9",
            ],
        ),
        (
            "x5.tblastn.sloc2_470",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject_loc",
                "2-470",
                "-outfmt",
                "5",
            ],
        ),
        (
            "x5.tblastn.ungapped_cbs0",
            &[
                "tblastn",
                "-query",
                "tests/db_search/of1/tnq.fa",
                "-subject",
                "tests/db_search/of1/tns.fa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "5",
            ],
        ),
        (
            "x5.blastx.ungapped_cbs0",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/x5_bx_q_mix.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "1",
                "-outfmt",
                "5",
            ],
        ),
        (
            "sa.e2e.blastp.invalid_o.15",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/e2e_invalid_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_protein_subject.faa",
                "-outfmt",
                "15",
            ],
        ),
        (
            "sa.x16.blastn.minus.15",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/x16n_q.fna",
                "-subject",
                "tests/fasta/outfmt0/x16n_s.fna",
                "-outfmt",
                "15",
            ],
        ),
        (
            "x16.blastx.ungapped_cbs0.16",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/x5_bx_q_mix.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "1",
                "-outfmt",
                "16",
            ],
        ),
        (
            "x16.blastx.ungapped_cbs0.15",
            &[
                "blastx",
                "-query",
                "tests/fasta/outfmt0/x5_bx_q_mix.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "1",
                "-outfmt",
                "15",
            ],
        ),
        (
            "sa.x16.tblastn.sloc.16",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject_loc",
                "2-470",
                "-outfmt",
                "16",
            ],
        ),
        (
            "sa.x16.tblastn.sloc.15",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/e2e_protein_query.faa",
                "-subject",
                "tests/fasta/outfmt0/e2e_tblastn_subject.fna",
                "-subject_loc",
                "2-470",
                "-outfmt",
                "15",
            ],
        ),
        (
            "sa.x16.tblastx.sub.15",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt0/x16x_q.fna",
                "-subject",
                "tests/fasta/outfmt0/x16x_s.fna",
                "-evalue",
                "0.01",
                "-max_target_seqs",
                "2",
                "-outfmt",
                "15",
            ],
        ),
        (
            "culling.tblastn.dom_cl2.11",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt0/culling_dom_q.faa",
                "-subject",
                "tests/fasta/outfmt0/culling_dom_s.fna",
                "-culling_limit",
                "2",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.sam.multi.17p",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/sam_multi_q.fna",
                "-subject",
                "tests/fasta/outfmt0/sam_multi_s.fna",
                "-task",
                "megablast",
                "-outfmt",
                "17",
            ],
        ),
        (
            "sa.sam.multi.17q",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt0/sam_multi_q.fna",
                "-subject",
                "tests/fasta/outfmt0/sam_multi_s.fna",
                "-task",
                "megablast",
                "-outfmt",
                "17 SQ",
            ],
        ),
        (
            "sa.arch.bn.g4.11",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt11/arch_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_s.fna",
                "-max_target_seqs",
                "1",
                "-dust",
                "20 64 1",
                "-query_loc",
                "10-300",
                "-subject_loc",
                "10-300",
                "-dbsize",
                "1000000",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.bn.qloc.11",
            &[
                "blastn",
                "-query",
                "tests/fasta/outfmt11/arch_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_s.fna",
                "-query_loc",
                "10-300",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.bp.sloc.11",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt11/arch_bp_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-subject_loc",
                "10-100",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.bp.titles.11",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt11/arch_bp_titles.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbn.qloc.11",
            &[
                "tblastn",
                "-query",
                "tests/fasta/outfmt11/arch_tbn_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_tbn_s.fna",
                "-query_loc",
                "10-100",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbx.def.11",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt11/arch_tbx_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_tbx_s.fna",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbx.qloc.11",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt11/arch_tbx_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_tbx_s.fna",
                "-query_loc",
                "10-100",
                "-outfmt",
                "11",
            ],
        ),
        (
            "sa.arch.tbx.stitles.11",
            &[
                "tblastx",
                "-query",
                "tests/fasta/outfmt11/arch_tbx_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_tbx_s_titles.fna",
                "-outfmt",
                "11",
            ],
        ),
        (
            "tab.bn.csv.10",
            &[
                "blastn",
                "-query",
                "tests/db_search/nslice_q.fna",
                "-subject",
                "tests/db_search/nslice.fna",
                "-outfmt",
                "10 std sstrand",
            ],
        ),
        (
            "tab.bn.csvhdr.20",
            &[
                "blastn",
                "-query",
                "tests/db_search/nslice_q.fna",
                "-subject",
                "tests/db_search/nslice.fna",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.bp.csv.10",
            &[
                "blastp",
                "-query",
                "tests/db_search/ids/ids_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "10 std",
            ],
        ),
        (
            "tab.bp.csvhdr.20",
            &[
                "blastp",
                "-query",
                "tests/db_search/ids/ids_q.faa",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.bx.csv.10",
            &[
                "blastx",
                "-query",
                "tests/db_search/ids/ids_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "10 std",
            ],
        ),
        (
            "tab.bx.csvhdr.20",
            &[
                "blastx",
                "-query",
                "tests/db_search/ids/ids_q.fna",
                "-subject",
                "tests/fasta/outfmt11/arch_bp_s.faa",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.tbn.csv.10",
            &[
                "tblastn",
                "-query",
                "tests/db_search/tab_tbn_q.faa",
                "-subject",
                "tests/db_search/tslice.fna",
                "-outfmt",
                "10 std sstrand",
            ],
        ),
        (
            "tab.tbn.csvhdr.20",
            &[
                "tblastn",
                "-query",
                "tests/db_search/tab_tbn_q.faa",
                "-subject",
                "tests/db_search/tslice.fna",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "tab.tbx.csv.10",
            &[
                "tblastx",
                "-query",
                "tests/db_search/tab_tbx_q.fna",
                "-subject",
                "tests/db_search/xslice.fna",
                "-outfmt",
                "10 std sstrand",
            ],
        ),
        (
            "tab.tbx.csvhdr.20",
            &[
                "tblastx",
                "-query",
                "tests/db_search/tab_tbx_q.fna",
                "-subject",
                "tests/db_search/xslice.fna",
                "-outfmt",
                "20 qseqid sseqid evalue",
            ],
        ),
        (
            "sxd.blastx.qloc.10",
            &[
                "blastx",
                "-query",
                "tests/db_search/bslice_q.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-query_loc",
                "1001-3000",
                "-outfmt",
                "10",
            ],
        ),
        (
            "sxd.blastx.qframe.20",
            &[
                "blastx",
                "-query",
                "tests/db_search/bslice_q.fna",
                "-subject",
                "tests/db_search/bslice.faa",
                "-query_loc",
                "1002-3000",
                "-outfmt",
                "20",
            ],
        ),
        (
            "ungapped.blastp.sc4a.8",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/sc4a_q.faa",
                "-subject",
                "tests/fasta/PajaWSV.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "5",
                "-outfmt",
                "8",
            ],
        ),
        (
            "ungapped.blastp.sc4a.11",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/sc4a_q.faa",
                "-subject",
                "tests/fasta/PajaWSV.faa",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-max_target_seqs",
                "5",
                "-outfmt",
                "11",
            ],
        ),
        (
            "ungapped.blastp.len0.5",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/len0_l0015_q.faa",
                "-subject",
                "tests/fasta/outfmt0/len0_l0015_s.faa",
                "-window_size",
                "0",
                "-word_size",
                "7",
                "-threshold",
                "9",
                "-matrix",
                "PAM30",
                "-evalue",
                "1e5",
                "-seg",
                "yes",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "5",
            ],
        ),
        (
            "ungapped.blastp.len0.8",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/len0_l0015_q.faa",
                "-subject",
                "tests/fasta/outfmt0/len0_l0015_s.faa",
                "-window_size",
                "0",
                "-word_size",
                "7",
                "-threshold",
                "9",
                "-matrix",
                "PAM30",
                "-evalue",
                "1e5",
                "-seg",
                "yes",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "8",
            ],
        ),
        (
            "ungapped.blastp.len0.15",
            &[
                "blastp",
                "-query",
                "tests/fasta/outfmt0/len0_l0015_q.faa",
                "-subject",
                "tests/fasta/outfmt0/len0_l0015_s.faa",
                "-window_size",
                "0",
                "-word_size",
                "7",
                "-threshold",
                "9",
                "-matrix",
                "PAM30",
                "-evalue",
                "1e5",
                "-seg",
                "yes",
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "15",
            ],
        ),
        (
            "ungapped.sam.mb.17q",
            &[
                "blastn",
                "-query",
                "tests/fixtures/blastn_regression/inputs/sam_ug_q.fa",
                "-subject",
                "tests/fixtures/blastn_regression/inputs/sam_ug_s.fa",
                "-ungapped",
                "-outfmt",
                "17 SQ",
            ],
        ),
    ];

    /// Stream 3 of a run whose own -outfmt is structured: NCBI's standard error `own` of its command
    /// line, then each line of `other` (stream 3 of the run without -outfmt: outfmt 0, 6 and 7) past
    /// the copies of it that `own` holds. The structured formats write none of the queries' warnings
    /// that outfmt 0, 6 and 7 write (line WQA's sweep, D3 review r3 M2, 2026-11-19). The lines
    /// come in the order of the run without -outfmt, which is stream 3's order (format by format,
    /// §5) for the rows that use this: none has a format that failed alone before another
    /// format's warning (review r2 of 2026-11-20, L2).
    fn and_the_lines_of(own: &[u8], other: &[u8]) -> String {
        let mut held: std::collections::HashMap<&[u8], usize> = std::collections::HashMap::new();
        for line in own.split_inclusive(|&byte| byte == b'\n') {
            *held.entry(line).or_default() += 1;
        }
        let mut stream = own.to_vec();
        for line in other.split_inclusive(|&byte| byte == b'\n') {
            match held.get_mut(line) {
                Some(copies) if *copies > 0 => *copies -= 1,
                _ => stream.extend_from_slice(line),
            }
        }
        text(&stream)
    }

    /// A SAM document without the CL tag of its `@PG` line (`db_search_fixtures.strip_sam_cl`, the
    /// rule of the SAM rows: the CLI writes its own command line there).
    fn without_cl(document: &[u8]) -> Vec<u8> {
        let text = String::from_utf8(document.to_vec()).unwrap();
        let head = format!("\n@PG\tID:0\tVN:{}", losat::web_adapter::LOSAT_VERSION);
        let Some(at) = text.find(&format!("{head}\tCL:")) else {
            return document.to_vec();
        };
        let start = at + head.len();
        let end = start + text[start..].find("\tPN:blastn\n").unwrap();
        format!("{}{}", &text[..start], &text[end..]).into_bytes()
    }

    // The output of an -outfmt besides 0, 6 and 7 is stream 4, NCBI's bytes of the same command
    // line (W11), and stream 3 NCBI's standard error of it, then the lines that the run's outfmt 0,
    // 6 and 7 add (`and_the_lines_of`); the formats 0, 6 and 7 and the HSP records are those of the
    // run without -outfmt. The archive's Bioseqs are the whole records
    // and the deflines as read; SAM's @PG line shows the run's command line.
    #[test]
    fn every_format_is_written_as_the_cli_writes_it() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        for (id, words) in FORMAT_ROWS {
            let at = |flag: &str| words[words.iter().position(|w| *w == flag).unwrap() + 1];
            let (query, subject) = (read(at("-query")), read(at("-subject")));
            let with = registered_streams(words, &query, &subject);
            let plain = registered_streams(&words[..words.len() - 2], &query, &subject);
            let mut ncbi =
                ncbi_output_in_losat_words(&read(&format!("tests/fixtures/outfmt0/{id}.out")));
            let mut ours = with.get(&STREAM_CUSTOM).cloned().unwrap_or_default();
            if at("-outfmt").starts_with("17") {
                let cl = format!("\tCL:losat {} \tPN:blastn\n", words.join(" "));
                assert!(text(&ours).contains(&cl), "{id}: {cl:?}");
                (ncbi, ours) = (without_cl(&ncbi), without_cl(&ours));
            }
            assert!(ours == ncbi, "{id}: stream 4");
            let stderr = ncbi_output_in_losat_words(
                &std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                    .unwrap_or_default(),
            );
            let diagnostics = with.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            let other = plain.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(
                text(&diagnostics),
                and_the_lines_of(&stderr, &other),
                "{id}: stream 3"
            );
            for stream in [0, 6, 7, STREAM_HITS] {
                assert!(
                    plain.get(&stream) == with.get(&stream),
                    "{id}: stream {stream}"
                );
            }
        }
    }

    /// The rows of NCBI's frozen outputs (`LOSAT/tests/outfmt0_manifest.tsv`) for which `keep(program,
    /// format number, files column)` holds, as (fixture id, argv): the arguments of
    /// `run_oracle.search_argv` (none of the rows taken here quotes a value).
    fn manifest_rows(keep: impl Fn(&str, &str, &str) -> bool) -> Vec<(String, Vec<String>)> {
        manifest_rows_with(|program, number, files, _| keep(program, number, files))
    }

    /// `manifest_rows` for which `keep(program, format number, files column, extra arguments)`
    /// holds.
    fn manifest_rows_with(
        keep: impl Fn(&str, &str, &str, &str) -> bool,
    ) -> Vec<(String, Vec<String>)> {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let manifest =
            std::fs::read_to_string(format!("{root}/tests/outfmt0_manifest.tsv")).unwrap();
        let mut lines = manifest.lines().filter(|line| !line.starts_with('#'));
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let column = |name: &str| header.iter().position(|h| *h == name).unwrap();
        let mut rows = Vec::new();
        for line in lines {
            let row: Vec<&str> = line.split('\t').collect();
            let field = |name: &str| row.get(column(name)).copied().unwrap_or("");
            let outfmt = match field("outfmt") {
                "" => "0",
                outfmt => outfmt,
            };
            let number = outfmt
                .trim()
                .split(' ')
                .next()
                .unwrap()
                .trim_start_matches('+');
            let number = number.trim_start_matches('0');
            if !keep(
                field("program"),
                number,
                field("files"),
                field("extra_args"),
            ) {
                continue;
            }
            let extra = field("extra_args");
            assert!(
                !extra.contains(['"', '\'', '\\']),
                "{}",
                field("fixture_id")
            );
            let mut argv = vec![
                field("program").to_string(),
                "-query".to_string(),
                field("query").to_string(),
                "-subject".to_string(),
                field("subject").to_string(),
            ];
            if !field("task").is_empty() {
                argv.extend(["-task".to_string(), field("task").to_string()]);
            }
            argv.extend(extra.split_whitespace().map(str::to_string));
            argv.extend(["-outfmt".to_string(), outfmt.to_string()]);
            rows.push((field("fixture_id").to_string(), argv));
        }
        rows
    }

    // Outfmt 13 and 14 are a ZIP archive on stream 4: NCBI's side files, named as the CLI names them
    // for `-out results`, and NCBI's master file (`results.json` / `results.xml`), which lists them by
    // those names; stream 3 is NCBI's standard error, then the lines that the run's outfmt 0, 6
    // and 7 add (`and_the_lines_of`), and the other streams are those of the run without -outfmt
    // (the 23 rows of the manifest with side files; `ungapped.blastp.sc4a.13` is BLASTP's ungapped
    // search, line NP item O4b).
    #[test]
    fn outfmt_13_and_14_are_an_archive_of_ncbis_files() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        let rows = manifest_rows(|_, _, files| files == "xinclude");
        assert_eq!(rows.len(), 23);
        let programs: std::collections::BTreeSet<&str> =
            rows.iter().map(|(_, argv)| argv[0].as_str()).collect();
        assert_eq!(programs.len(), 5);
        for (id, argv) in &rows {
            let words: Vec<&str> = argv.iter().map(String::as_str).collect();
            let at = |flag: &str| words[words.iter().position(|w| *w == flag).unwrap() + 1];
            let number: u32 = at("-outfmt").trim().parse().unwrap();
            let extension = if number == 13 { "json" } else { "xml" };
            let (query, subject) = (read(at("-query")), read(at("-subject")));
            let with = registered_streams(&words, &query, &subject);
            let plain = registered_streams(&words[..words.len() - 2], &query, &subject);
            let mut expected = Vec::new();
            for n in 1.. {
                let path = format!("{root}/tests/fixtures/outfmt0/{id}.out_{n}.{extension}");
                let Ok(side) = std::fs::read(path) else {
                    break;
                };
                expected.push((
                    format!("results_{n}.{extension}").into_bytes(),
                    ncbi_output_in_losat_words(&side),
                ));
            }
            let master = text(&ncbi_output_in_losat_words(&read(&format!(
                "tests/fixtures/outfmt0/{id}.out"
            ))))
            .replace(&format!("{id}.out_"), "results_");
            expected.push((master_name(number).into_bytes(), master.into_bytes()));
            let archive = with.get(&STREAM_CUSTOM).cloned().unwrap_or_default();
            let files = crate::zip::read(&archive);
            assert_eq!(files.len(), expected.len(), "{id}: entries");
            for ((name, bytes), (want_name, want)) in files.iter().zip(&expected) {
                assert_eq!(text(name), text(want_name), "{id}");
                assert!(bytes == want, "{id}: {}", text(name));
            }
            let stderr = ncbi_output_in_losat_words(
                &std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                    .unwrap_or_default(),
            );
            let diagnostics = with.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            let other = plain.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(
                text(&diagnostics),
                and_the_lines_of(&stderr, &other),
                "{id}: stream 3"
            );
            for stream in [0, 6, 7, STREAM_HITS] {
                assert!(
                    plain.get(&stream) == with.get(&stream),
                    "{id}: stream {stream}"
                );
            }
        }
    }

    // BLASTX writes the BLAST archive (outfmt 11) as NCBI does: the 60 rows of the manifest that a
    // run takes (not the one with a `>?` gap line, which the store refuses: maintainer decision 3),
    // stream 4 and stream 3 (NCBI's standard error, then the lines that the run's outfmt 0, 6 and 7
    // add: `and_the_lines_of`).
    #[test]
    fn blastx_writes_ncbis_archive() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../LOSAT");
        let read = |path: &str| std::fs::read(format!("{root}/{path}")).unwrap();
        let gap_line = |path: &str| {
            read(path)
                .split(|b| *b == b'\n')
                .any(|l| l.starts_with(b">?"))
        };
        let rows = manifest_rows(|program, number, _| program == "blastx" && number == "11");
        let rows: Vec<_> = rows
            .into_iter()
            .filter(|(_, argv)| {
                let file = |flag: &str| &argv[argv.iter().position(|w| w == flag).unwrap() + 1];
                !gap_line(file("-query")) && !gap_line(file("-subject"))
            })
            .collect();
        assert_eq!(rows.len(), 60);
        for (id, argv) in &rows {
            let words: Vec<&str> = argv.iter().map(String::as_str).collect();
            let at = |flag: &str| words[words.iter().position(|w| *w == flag).unwrap() + 1];
            let (query, subject) = (read(at("-query")), read(at("-subject")));
            let with = registered_streams(&words, &query, &subject);
            let plain = registered_streams(&words[..words.len() - 2], &query, &subject);
            let ncbi =
                ncbi_output_in_losat_words(&read(&format!("tests/fixtures/outfmt0/{id}.out")));
            let ours = with.get(&STREAM_CUSTOM).cloned().unwrap_or_default();
            assert!(ours == ncbi, "{id}: stream 4");
            let stderr = ncbi_output_in_losat_words(
                &std::fs::read(format!("{root}/tests/fixtures/outfmt0/{id}.err"))
                    .unwrap_or_default(),
            );
            let diagnostics = with.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            let other = plain.get(&STREAM_DIAGNOSTICS).cloned().unwrap_or_default();
            assert_eq!(
                text(&diagnostics),
                and_the_lines_of(&stderr, &other),
                "{id}: stream 3"
            );
        }
    }

    /// A BLASTN argv that writes the BLAST archive (`whole`) or not.
    fn argv_writing_archive(whole: bool) -> &'static [&'static str] {
        if whole {
            &["blastn", "-query", "q", "-subject", "s", "-outfmt", "11"]
        } else {
            &["blastn", "-query", "q", "-subject", "s"]
        }
    }

    // The store keeps the deflines that the archive cannot make from the title (a control byte, white
    // space after the `>` or at the end; not a CRLF line end), aside from the records: only a run that
    // writes the archive gives the records a source, which they hold no longer after it.
    #[test]
    fn the_store_keeps_the_deflines_that_the_title_does_not_give() {
        let fasta: &[u8] = b">a one\r\nACGT\n>b\twith tab\nGG\n>  c\nTT\n>d \nAA\n>e\r\nCC\n";
        let (handle, _) = crate::store::register("blastn", 0, true, fasta).unwrap();
        let (subject, _) = crate::store::register("blastn", 1, true, b">s\nACGTACGT\n").unwrap();
        let sources = |whole: bool| {
            crate::store::with_inputs(
                Program::Blastn,
                handle,
                subject,
                argv_writing_archive(whole),
                Some(true),
                |queries, _, subjects| {
                    let source = |record: &FastaRecord| {
                        record.source.as_ref().map(|source| source.defline.clone())
                    };
                    (
                        queries.iter().map(source).collect::<Vec<_>>(),
                        subjects.iter().map(source).collect::<Vec<_>>(),
                    )
                },
            )
            .unwrap()
        };
        assert_eq!(sources(false), (vec![None; 5], vec![None]));
        let kept = |line: &[u8]| Some(line.to_vec());
        assert_eq!(
            sources(true),
            (
                vec![
                    kept(b""),
                    kept(b">b\twith tab"),
                    kept(b">  c"),
                    kept(b">d "),
                    kept(b"")
                ],
                vec![kept(b"")]
            )
        );
        assert_eq!(sources(false), (vec![None; 5], vec![None]));
        crate::store::release(handle).unwrap();
        crate::store::release(subject).unwrap();
        // A first record without a defline, a kept defline at the end of the input without a line
        // end, and CR line ends.
        for (fasta, kept_lines) in [
            (
                &b"ACGT\n>b\tx\nGG\n"[..],
                vec![b"".to_vec(), b">b\tx".to_vec()],
            ),
            (&b">a\nAC\n>b \t"[..], vec![b"".to_vec(), b">b \t".to_vec()]),
            (
                &b">a \rAC\r>b\rGG\r"[..],
                vec![b">a ".to_vec(), b"".to_vec()],
            ),
        ] {
            let (handle, _) = crate::store::register("blastn", 0, true, fasta).unwrap();
            let (subject, _) =
                crate::store::register("blastn", 1, true, b">s\nACGTACGT\n").unwrap();
            let defline = |whole: bool| {
                crate::store::with_inputs(
                    Program::Blastn,
                    handle,
                    subject,
                    argv_writing_archive(whole),
                    Some(true),
                    |queries, _, _| {
                        queries
                            .iter()
                            .map(|record| {
                                record.source.as_ref().map(|source| source.defline.clone())
                            })
                            .collect::<Vec<_>>()
                    },
                )
                .unwrap()
            };
            assert_eq!(defline(false), vec![None; kept_lines.len()], "{fasta:?}");
            assert_eq!(
                defline(true),
                kept_lines.into_iter().map(Some).collect::<Vec<_>>(),
                "{fasta:?}"
            );
            crate::store::release(handle).unwrap();
            crate::store::release(subject).unwrap();
        }
    }

    const SPLIT_QUERY: &[u8] =
        include_bytes!("../../../LOSAT/tests/fasta/outfmt0/numdesc_split_q.fna");
    const SPLIT_SUBJECT: &[u8] =
        include_bytes!("../../../LOSAT/tests/fasta/outfmt0/numdesc_split_s.fna");
    const SPLIT_OUT0: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_na3.out");
    const SPLIT_OUT6: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_na3.6.out");
    const SPLIT_ERR6: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_na3.6.err");

    // The one input where outfmt 0 and outfmt 6/7 keep other subjects with -num_alignments 3
    // (the hit list of outfmt 0 is 500, the tabular formats' 3, and the subject A_split, whose
    // gap the preliminary X-drop cannot cross, is evicted from the small list): two searches.
    // Each stream is NCBI's frozen bytes of its format alone; the HSP records are the
    // tabular list, and only the HSPs that outfmt 0 shows identically have its sections.
    #[test]
    fn two_searches_write_each_format_as_the_cli_does() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        let words = [
            "blastn",
            "-query",
            "tests/fasta/outfmt0/numdesc_split_q.fna",
            "-subject",
            "tests/fasta/outfmt0/numdesc_split_s.fna",
            "-task",
            "blastn",
            "-num_alignments",
            "3",
        ];
        let streams = streams(&words, &query, &subject).unwrap();
        assert_eq!(streams[&0], ncbi_output_in_losat_words(SPLIT_OUT0));
        assert_eq!(streams[&6], ncbi_output_in_losat_words(SPLIT_OUT6));
        // Neither search failed: no signal of stream 9 (line WQA's sweep, D7, 2026-11-19).
        assert_eq!(discarded(&streams), [] as [u32; 0]);
        // The warning of the tabular formats, once (outfmt 0 has none).
        assert_eq!(streams[&3], ncbi_output_in_losat_words(SPLIT_ERR6));
        let out0 = String::from_utf8(streams[&0].clone()).unwrap();
        let headings: Vec<&str> = out0.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(headings, ["> A_split", "> D00", "> D01"]);
        let records = String::from_utf8(streams[&1].clone()).unwrap();
        let field = |line: &str, name: &str| {
            let at = line.find(&format!("\"{name}\":")).unwrap() + name.len() + 3;
            line[at..]
                .split(&[',', '}'][..])
                .next()
                .unwrap()
                .to_string()
        };
        let lines: Vec<&str> = records.lines().collect();
        // D00 D01 D02 are the subjects 1, 2 and 3; D02 is not in outfmt 0's list.
        assert_eq!(
            lines.iter().map(|l| field(l, "s_idx")).collect::<Vec<_>>(),
            ["1", "2", "3"]
        );
        let pointed: Vec<bool> = lines.iter().map(|l| !l.contains("\"out0\":null")).collect();
        assert_eq!(pointed, [true, true, false]);
        assert!(lines[2].contains("\"out0_subject\":null"));
        let out0_subject = |line: &str| {
            let at = line.find("\"out0_subject\":[").unwrap() + 16;
            let range: Vec<usize> = line[at..]
                .split(']')
                .next()
                .unwrap()
                .split(',')
                .map(|n| n.parse().unwrap())
                .collect();
            out0[range[0]..range[1]].lines().next().unwrap().to_string()
        };
        assert_eq!(out0_subject(lines[0]), "> D00");
        assert_eq!(out0_subject(lines[1]), "> D01");
        // The output of -outfmt "6 std" comes from the tabular search: NCBI's outfmt 6 with
        // -num_alignments 3, without A_split, which the hit list of outfmt 0 keeps. The other
        // streams do not change.
        let custom = [&words[..], &["-outfmt", "6 std"]].concat();
        let with = self::streams(&custom, &query, &subject).unwrap();
        assert_eq!(with[&STREAM_CUSTOM], ncbi_output_in_losat_words(SPLIT_OUT6));
        for stream in [0, 6, 7, STREAM_HITS, STREAM_DIAGNOSTICS] {
            assert_eq!(with.get(&stream), streams.get(&stream), "stream {stream}");
        }
    }

    // The same command line without a difference of the hit list sizes is one search.
    #[test]
    fn equal_hit_list_sizes_are_one_search() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        fn words<'a>(extra: &[&'a str]) -> Vec<&'a str> {
            [
                &["blastn", "-query", "q", "-subject", "s", "-task", "blastn"][..],
                extra,
            ]
            .concat()
        }
        let plain = streams(&words(&[]), &query, &subject).unwrap();
        // -num_alignments 500: format 0 keeps the default hit list of 500 as the tabular
        // formats do (one search, and no warning).
        let same = streams(&words(&["-num_alignments", "500"]), &query, &subject).unwrap();
        for stream in [0, 6, 7, 1] {
            assert_eq!(plain[&stream], same[&stream], "stream {stream}");
        }
        assert!(plain.get(&3).is_none_or(|bytes| bytes.is_empty()));
        assert!(same.get(&3).is_none_or(|bytes| bytes.is_empty()));
        // -num_descriptions 500 is the other default that makes one search; outfmt 6 and 7
        // ignore it with NCBI's warning, which stream 3 carries once.
        let ignored = streams(&words(&["-num_descriptions", "500"]), &query, &subject).unwrap();
        for stream in [0, 6, 7, 1] {
            assert_eq!(plain[&stream], ignored[&stream], "stream {stream}");
        }
        assert_eq!(
            String::from_utf8(ignored[&3].clone()).unwrap(),
            "Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n"
        );
    }

    // The CLI's standard error for `-outfmt 6` with `-num_descriptions 3`: the formatting
    // handler's two warnings, then NCBI's error for a hit list of 0 (`-num_alignments 0`).
    const ERR6_ND3_NA0: &str = "Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\nWarning: [blastn] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n";
    const SPLIT_ND3_NA0_OUT0: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_nd3_na0.out");
    const SPLIT_ND3_NA0_ERR0: &[u8] =
        include_bytes!("../../../LOSAT/tests/fixtures/outfmt0/numdesc.blastn.split_nd3_na0.err");

    // -num_alignments 0 gives the tabular formats a hit list of 0, which the CLI rejects with
    // NCBI's error for -outfmt 6 and 7; format 0 is written as the CLI writes it and the run
    // does not fail. Stream 0 is NCBI's frozen bytes, streams 6, 7 and 1 are empty, and stream 3
    // is the CLI's standard error for the argv (outfmt 0: the handler's `Examining 5 or more
    // matches`, for the hit list size 3) followed by the lines of that for outfmt 6 that it does not
    // hold (the handler's `-num_descriptions` line), then its error, each text once (line WQA's
    // sweep, D3 review r2, 2026-11-19: the argv's lines first).
    #[test]
    fn a_hit_list_size_of_zero_fails_the_tabular_formats_only() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        let words = [
            "blastn",
            "-query",
            "tests/fasta/outfmt0/numdesc_split_q.fna",
            "-subject",
            "tests/fasta/outfmt0/numdesc_split_s.fna",
            "-task",
            "blastn",
            "-num_descriptions",
            "3",
            "-num_alignments",
            "0",
        ];
        let streams = streams(&words, &query, &subject).unwrap();
        assert_eq!(
            text(&streams[&0]),
            text(&ncbi_output_in_losat_words(SPLIT_ND3_NA0_OUT0))
        );
        for stream in [1, 6, 7] {
            assert!(
                streams.get(&stream).is_none_or(|bytes| bytes.is_empty()),
                "stream {stream}"
            );
        }
        // The host is told that streams 6 and 7 failed (stream 9; line WQA's sweep, D7, 2026-11-19).
        assert_eq!(discarded(&streams), [6, 7]);
        // NCBI's frozen bytes and texts in LOSAT's words (AGENTS.md rule 8).
        assert!(ncbi_in_losat_words(ERR6_ND3_NA0)
            .contains(&text(&losat::web_adapter::losat_bytes(SPLIT_ND3_NA0_ERR0))));
        assert_eq!(
            text(&streams[&3]),
            losat::web_adapter::losat_text(&format!(
                "{}{}",
                text(SPLIT_ND3_NA0_ERR0),
                "Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n"
            ))
        );
        // Both formats' hit lists of 0 fail every format, as the CLI does.
        let both = [
            &words[..7],
            &["-num_descriptions", "0", "-num_alignments", "0"],
        ]
        .concat();
        let error = run(&both, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
        assert!(error.contains("No hits are being saved"), "{error}");
        assert!(validate(&both)
            .unwrap_err()
            .contains("No hits are being saved"));
        // validate accepts -num_alignments 0: the format 0 search, the argv's own, is valid.
        validate(&words).unwrap();
        // The argv's own -outfmt 6, 7 or a tabular -outfmt with fields is the tabular search's: the
        // run fails as the CLI's run of that argv (exit 1), whose standard error is stream 3 and
        // then the run's error, and validate rejects it with that error (native LOSAT 16bef7acd,
        // 2026-11-19; line WQA's sweep, D3 review, 2026-11-19).
        for outfmt in ["6", "7", "6 qseqid sseqid"] {
            let tabular = [&words[..], &["-outfmt", outfmt]].concat();
            let (error, failed) = failed_streams(&tabular, &query, &subject);
            assert_eq!(
                format!("{}{error}", diagnostics_of(&failed)),
                ncbi_in_losat_words(ERR6_ND3_NA0),
                "-outfmt {outfmt}"
            );
            assert_outputs_empty(&failed, &format!("-outfmt {outfmt}"));
            assert_eq!(validate(&tabular).unwrap_err(), error, "-outfmt {outfmt}");
        }
    }

    // A -num_descriptions of 2,147,483,598 or more makes BLASTN's preliminary hit list size of
    // format 0 overflow NCBI's 32-bit `int` (N + 50; NCBI crashes), while the tabular formats ignore
    // -num_descriptions. LOSAT computes the size without overflow (PD-LOSAT-NCBI-DEFECTS, hit list
    // sizes past `Int4`), so the run writes every format, with the argv's `-outfmt 6` and without
    // `-outfmt`; these inputs have fewer hits than either size, so each format is the plain run's.
    // Stream 3 is the formatting handler's warning of the tabular formats, as the CLI's standard
    // error of `-outfmt 6`.
    #[test]
    fn a_hit_list_size_past_ncbi_s_int_writes_every_format() {
        let query = fasta(SPLIT_QUERY);
        let subject = fasta(SPLIT_SUBJECT);
        let plain_words = ["blastn", "-query", "q", "-subject", "s", "-task", "blastn"];
        let plain = streams(&plain_words, &query, &subject).unwrap();
        assert_eq!(discarded(&plain), [] as [u32; 0]);
        let words = [&plain_words[..], &["-num_descriptions", "2147483647"]].concat();
        let tabular = [&words[..], &["-outfmt", "6"]].concat();
        for argv in [&words, &tabular] {
            validate(argv).unwrap();
            let streams = streams(argv, &query, &subject).unwrap();
            assert_eq!(discarded(&streams), [] as [u32; 0], "{argv:?}");
            for stream in [0, 1, 6, 7] {
                assert_eq!(streams[&stream], plain[&stream], "{argv:?} stream {stream}");
            }
            assert_eq!(
                String::from_utf8(streams[&3].clone()).unwrap(),
                "Warning: [blastn] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n",
                "{argv:?}"
            );
        }
    }

    // The argv's own search decides the run (`own_search`; line WQA's sweep, D3 review, 2026-11-19).
    // `-num_alignments 0` gives the tabular formats a hit list of 0, which the CLI rejects. With the
    // argv's `-outfmt 6`, `7` or a tabular `-outfmt` with fields (the first search) the run fails as
    // the CLI fails (exit 1): stream 3 and then the error are the CLI's standard error, the
    // handler's warning and NCBI's error, and no output is kept; validate rejects the argv with that
    // error (line WQA's sweep, D3 review r2, 2026-11-19: BLASTP's, TBLASTN's and TBLASTX's left the
    // hit list to the run). Without `-outfmt` (outfmt 0, the second search) the run completes:
    // format 0 is written, the tabular formats fail alone, and stream 3 is the CLI's standard error
    // for the argv (none) followed by what that for `-outfmt 6` adds (native LOSAT 16bef7acd,
    // 2026-11-19, every program).
    #[test]
    fn the_argv_s_own_search_decides_the_run() {
        let protein: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let nucleotides = random_nucleotides(300, 11);
        let cases: [(&str, &[u8], &[u8]); 5] = [
            ("blastp", protein, protein),
            ("tblastn", protein, CODING),
            ("blastn", &nucleotides, &nucleotides),
            ("tblastx", CODING, CODING),
            ("blastx", CODING, protein),
        ];
        for (program, query, subject) in cases {
            let query = [FastaRecord::new("q", b"q", query)];
            let subject = [FastaRecord::new("s", b"s", subject)];
            let words = [
                program,
                "-query",
                "q",
                "-subject",
                "s",
                "-num_alignments",
                "0",
            ];
            let cli = ncbi_in_losat_words(&format!("Warning: [{program}] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n"));
            for outfmt in ["6", "7", "6 qseqid sseqid"] {
                let argv = [&words[..], &["-outfmt", outfmt]].concat();
                let label = format!("{program} -outfmt {outfmt}");
                let (error, failed) = failed_streams(&argv, &query, &subject);
                assert_eq!(
                    format!("{}{error}", diagnostics_of(&failed)),
                    cli,
                    "{label}"
                );
                assert_outputs_empty(&failed, &label);
                assert_eq!(validate(&argv).unwrap_err(), error, "{label}");
            }
            validate(&words).unwrap();
            let completed = streams(&words, &query, &subject).unwrap();
            assert!(
                text(&completed[&0]).contains("Sequences producing significant alignments"),
                "{program}"
            );
            for stream in [1, 6, 7] {
                assert!(
                    completed.get(&stream).is_none_or(Vec::is_empty),
                    "{program} {stream}"
                );
            }
            assert_eq!(text(&completed[&3]), cli, "{program}");
        }
    }

    // A run of two searches that fails gives the CLI's standard error for the argv up to its error
    // (line WQA's sweep, D3 review M1, 2026-11-19): the warnings of the argv's own search, not
    // those of the other search, which writes formats that the argv does not ask for. A subject
    // whose one record has no residues fails BLASTN and TBLASTX with `-num_descriptions 5` (format
    // 0 has a hit list size of its own): the CLI's outfmt 0 writes the subject's warning and the
    // error, its `-outfmt 6` the handler's `-num_descriptions` line first. A BLASTP subject titled
    // `, ,` fails the CLI's outfmt 0 with LOSAT's error and no warning; its `-outfmt 6` completes,
    // where stream 0 fails alone (native LOSAT 16bef7acd, 2026-11-19).
    #[test]
    fn a_failed_run_of_two_searches_gives_the_warnings_of_the_argv_s_own_search() {
        let letters = random_nucleotides(300, 11);
        let query = [FastaRecord::new("Query_1", b"q1", &letters)];
        let empty = [FastaRecord::new("Subject_1", b"empty record", b"")];
        let ignored = |program: &str| {
            format!("Warning: [{program}] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n")
        };
        for program in ["blastn", "tblastx"] {
            let words = [
                program,
                "-query",
                "q",
                "-subject",
                "s",
                "-num_descriptions",
                "5",
            ];
            let warning = format!(
                "Warning: [{program}] Subject_1 empty record: Subject sequence contains no data\n"
            );
            for (outfmt, expected) in [
                (&[][..], warning.clone()),
                (
                    &["-outfmt", "6"][..],
                    format!("{}{warning}", ignored(program)),
                ),
            ] {
                let label = format!("{program} {outfmt:?}");
                let argv = [&words[..], outfmt].concat();
                let (error, failed) = failed_streams(&argv, &query, &empty);
                assert_eq!(
                    error.trim_end(),
                    ncbi_in_losat_words(
                        "BLAST engine error: The average subject length is too short"
                    ),
                    "{label}"
                );
                assert_eq!(diagnostics_of(&failed), expected, "{label}");
                assert_outputs_empty(&failed, &label);
            }
        }
        let protein: &[u8] = b"MKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMA";
        let query = [FastaRecord::new(
            "Query_1",
            b"p",
            b"MKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMAWRSTY",
        )];
        let subject = [
            FastaRecord::new("Subject_1", b"s1 ok", protein),
            FastaRecord::new("Subject_2", b", ,", protein),
        ];
        let words = [
            "blastp",
            "-query",
            "q",
            "-subject",
            "s",
            "-num_descriptions",
            "5",
        ];
        let (error, failed) = failed_streams(&words, &query, &subject);
        assert!(
            error.starts_with("subject record 2 has a title of punctuation"),
            "{error}"
        );
        assert_eq!(diagnostics_of(&failed), "");
        assert_outputs_empty(&failed, "blastp");
        let tabular = streams(&[&words[..], &["-outfmt", "6"]].concat(), &query, &subject).unwrap();
        assert!(!tabular[&6].is_empty());
        assert!(tabular.get(&0).is_none_or(Vec::is_empty));
        assert_eq!(
            text(&tabular[&3]),
            format!("{}Error: {error}\n", ignored("blastp"))
        );
    }

    // A coding sequence of the protein of `q_idx_is_the_registered_index`.
    const CODING: &[u8] = b"ATGAAAACCGCGTATATTGCGAAACAGCGTCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGTCAGCTGGAAGAACGTCTGGGCCTGATTGAAGTGCAGGCGCCGATTCTGAGCCGTGTGGGCGATGGCACCCAGGATAACCTGAGCGGCGCGGAAAAAGCGGTGCAGGTGAAAGTGAAAGCGCTGCCGGATGCGCAG";

    // Line RG item T7 (PD-LOSAT-NCBI-DEFECTS 1.37): a subject titled `Subject_9 bad x\x81y title`,
    // whose title NCBI's `GenerateDefline` cannot decode (a `Subject_` first word and a byte in no
    // encoding NCBI guesses). Each stream is the CLI's run of its format: the default tabular
    // formats stop where NCBI stops (at the row's `saccver`, a null pointer), outfmt 0 writes
    // NCBI's line in place of the alignment, a field list without `sacc` and `saccver` writes
    // NCBI's catch rows. Streams 6 and 7 that the argv does not ask for fail alone in the run's
    // one search (`FormatRefusal::FailsAlone`): they and the HSP records (stream 1, which point
    // into stream 6) stay empty, stream 3 has the refusal once after the warnings, and the other
    // formats are written. The argv's own `-outfmt 6`, `7` or `sacc` fails the run, as the CLI
    // fails (exit 1). BLASTX reads the title as its bytes as the other programs do.
    #[test]
    fn a_subject_title_ncbi_cannot_decode_fails_the_tabular_streams_the_argv_does_not_ask_for() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let title = b"Subject_9 bad x\x81y title";
        let refusal = "subject record 1 (Subject_9)";
        for (program, query, subject) in [
            ("blastn", CODING, CODING),
            ("tblastx", CODING, CODING),
            ("blastx", CODING, &protein[..]),
        ] {
            let query = [FastaRecord::new("Query_1", b"q", query)];
            let subject = [FastaRecord::new("Subject_1", title, subject)];
            let words = [program, "-query", "q", "-subject", "s"];
            let plain = streams(&words, &query, &subject).unwrap();
            assert!(
                text(&plain[&0]).contains(
                    "Sequence with id Subject_1 no longer exists in database...alignment skipped"
                ),
                "{program}"
            );
            for stream in [1, 6, 7] {
                assert!(
                    plain.get(&stream).is_none_or(|bytes| bytes.is_empty()),
                    "{program}: stream {stream}"
                );
            }
            assert_eq!(discarded(&plain), [6, 7], "{program}");
            let error = text(&plain[&3]);
            assert_eq!(error.matches(refusal).count(), 1, "{program}: {error}");
            assert!(
                error.contains("null pointer exception (exit 255)"),
                "{program}: {error}"
            );
            let fields = [&words[..], &["-outfmt", "6 qseqid sseqid sallseqid"]].concat();
            let custom = streams(&fields, &query, &subject).unwrap();
            assert!(
                text(&custom[&STREAM_CUSTOM]).starts_with("q\tUnknown\tSubject_1\n"),
                "{program}: {}",
                text(&custom[&STREAM_CUSTOM])
            );
            assert!(!custom[&0].is_empty(), "{program}");
            assert!(text(&custom[&3]).contains(refusal), "{program}");
            for argv in [
                &["-outfmt", "6"][..],
                &["-outfmt", "7"],
                &["-outfmt", "6 qseqid sacc"],
            ] {
                let stops = [&words[..], argv].concat();
                let error =
                    run(&stops, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
                assert!(
                    error.contains(refusal) && error.contains("null pointer exception"),
                    "{program} {argv:?}: {error}"
                );
            }
        }
    }

    // NCBI BLAST+ 2.17.0 writes BLASTP's `-outfmt 6` of a subject titled `, ,` (exit 0: `,` is its
    // `sseqid`, line RG item T7 review 2 `web/n6.out`) and crashes in outfmt 0, where
    // `x_CleanAndCompress` reads past the title's end; the CLI refuses outfmt 0 (rule 3). A run
    // whose argv asks for `-outfmt 6` or a field list writes those, and stream 0 fails alone:
    // empty, its refusal on stream 3 after the warnings of the CLI's `-outfmt 6` (those of a second
    // query, whose O letters are replaced by X: stream 0 fails before it takes them), the HSP
    // records without `out0`. The default `-outfmt 0` fails the run, as the CLI fails.
    #[test]
    fn a_blastp_title_of_punctuation_fails_stream_0_alone_where_the_argv_asks_for_another_format() {
        let protein = b"CFDYQGKGFISMYSLRLTSNPFPEIFSRVERPMRGRLALIDVRHRVPYACGDESHIVAGRPARPTLTRSVRHTISPDPFFVCRCRRDYSIQYMLQPRLLRPSENVESGTCANLVMRAVVKVSAYDFAVKTSNDFDGGVCERDKSYVTLIPSGLSTPFVNLVRHAMQLKVPSRLAACYTHLAALKGLGSPRLTEYNRCSTHSGRGLFEPESTLGDPSRVQIPLTETRCNSRLVTPITNYCTSKSREPRHSCNLVRVNVILICNSARGYLCKRSRGTEYHKFTSSEFSSPPWMNLTLVNYQN";
        let query = [
            FastaRecord::new("Query_1", b"q1 one", protein),
            FastaRecord::new("Query_2", b"q2 two", b"OOOOOOOOOOOOOOOOOOOO"),
        ];
        let subject = [FastaRecord::new("Subject_1", b", ,", protein)];
        let words = ["blastp", "-query", "q", "-subject", "s"];
        let refusal = "subject record 1 has a title of punctuation";
        let ncbi = "q1\t,\t100.000\t300\t0\t0\t1\t300\t1\t300\t0.0\t618\n";
        let six = [&words[..], &["-outfmt", "6"]].concat();
        // The warnings of the same run with a title that NCBI writes in every format.
        let good = [FastaRecord::new("Subject_1", b"s1 good", protein)];
        let warnings = text(&streams(&six, &query, &good).unwrap()[&3]);
        assert!(
            warnings.contains("O characters replaced by X"),
            "{warnings}"
        );
        let written = streams(&six, &query, &subject).unwrap();
        assert_eq!(text(&written[&6]), ncbi);
        assert!(text(&written[&7]).contains(ncbi));
        assert!(written.get(&0).is_none_or(|bytes| bytes.is_empty()));
        let error = text(&written[&3]);
        assert!(error.starts_with(&warnings), "{error}");
        assert_eq!(error.matches(refusal).count(), 1, "{error}");
        let records = text(&written[&STREAM_HITS]);
        assert!(
            records.contains("\"out6\":[0,") && records.contains("\"out0\":null"),
            "{records}"
        );
        let fields = [&words[..], &["-outfmt", "6 qseqid sseqid evalue"]].concat();
        let custom = streams(&fields, &query, &subject).unwrap();
        assert_eq!(text(&custom[&STREAM_CUSTOM]), "q1\t,\t0.0\n");
        assert_eq!(text(&custom[&6]), ncbi);
        assert!(custom.get(&0).is_none_or(|bytes| bytes.is_empty()));
        let error = run(&words, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
        assert!(error.contains(refusal), "{error}");
    }

    // BLASTX searches a nucleotide query against protein subjects and writes every format, the
    // HSP records and the warnings as the other programs do, one search for each hit list size:
    // with -num_alignments 0 the tabular search fails alone with NCBI's error on stream 3, and
    // format 0 is written.
    #[test]
    fn blastx_runs_as_the_other_programs_do() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [FastaRecord::new("q", b"q", CODING)];
        let subject = [FastaRecord::new("s", b"s", protein)];
        let words = ["blastx", "-query", "q.fna", "-subject", "s.faa"];
        let plain = streams(&words, &query, &subject).unwrap();
        let row = text(&plain[&6]);
        assert!(
            row.starts_with("q\ts\t100.000\t66\t0\t0\t1\t198\t1\t66\t"),
            "{row}"
        );
        assert_eq!(row.lines().count(), 1);
        assert!(text(&plain[&7]).contains("# BLASTX "));
        assert!(text(&plain[&0]).starts_with("BLASTX "));
        let records = text(&plain[&1]);
        assert_eq!(records.lines().count(), 1);
        assert!(records.contains("\"q_idx\":0"), "{records}");
        assert!(!records.contains("\"out0\":null"), "{records}");
        assert!(plain.get(&3).is_none_or(|bytes| bytes.is_empty()));
        let zero = streams(
            &[&words[..], &["-num_alignments", "0"]].concat(),
            &query,
            &subject,
        )
        .unwrap();
        assert!(text(&zero[&0]).contains("Sequences producing significant alignments"));
        for stream in [1, 6, 7] {
            assert!(
                zero.get(&stream).is_none_or(|bytes| bytes.is_empty()),
                "stream {stream}"
            );
        }
        // The CLI's standard error for `-outfmt 6`: the hit list size's warning, then NCBI's error.
        assert_eq!(
            text(&zero[&3]),
            ncbi_in_losat_words("Warning: [blastx] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n")
        );
        // Both searches fail, each for its own reason (the tabular formats' hit list is checked
        // before the e-value): the run fails with the error of the argv's own search, as validate,
        // which is the CLI's for the argv: the e-value's for outfmt 0, NCBI's of the hit list
        // after the handler's warning for `-outfmt 6` (native LOSAT 16bef7acd, 2026-11-19).
        let both = [&words[..], &["-num_alignments", "0", "-evalue", "0"]].concat();
        let (error, failed) = failed_streams(&both, &query, &subject);
        assert_eq!(
            format!("{}{error}", diagnostics_of(&failed)),
            ncbi_in_losat_words("BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the BLAST+ user manual.\n")
        );
        assert_eq!(validate(&both).unwrap_err(), error);
        let tabular = [&both[..], &["-outfmt", "6"]].concat();
        let (error, failed) = failed_streams(&tabular, &query, &subject);
        assert_eq!(
            format!("{}{error}", diagnostics_of(&failed)),
            ncbi_in_losat_words("Warning: [blastx] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n")
        );
        assert_eq!(validate(&tabular).unwrap_err(), error);
    }

    // BLASTX reads -query_loc and -subject_loc as the CLI does (DW-39): the rows give record
    // coordinates, a query that the range skips keeps its registered index in the HSP records,
    // and validate reads the ranges where the run reads them, after LOSAT's limits of the
    // options and the subject range first.
    #[test]
    fn blastx_ranges_as_the_other_programs() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let query = [
            FastaRecord::new("short", b"short", b"ACGTACGT"),
            FastaRecord::new("q", b"q", CODING),
        ];
        let subject = [FastaRecord::new("s", b"s", protein)];
        let words = [
            "blastx",
            "-query",
            "q.fna",
            "-subject",
            "s.faa",
            "-query_loc",
            "31-198",
            "-subject_loc",
            "5-66",
        ];
        assert!(validate(&words).is_ok());
        let ranged = streams(&words, &query, &subject).unwrap();
        let row = text(&ranged[&6]);
        assert!(
            row.starts_with("q\ts\t100.000\t56\t0\t0\t31\t198\t11\t66\t"),
            "{row}"
        );
        let records = text(&ranged[&1]);
        assert_eq!(records.lines().count(), 1);
        assert!(records.contains("\"q_idx\":1"), "{records}");
        assert_eq!(
            searched_queries(Some(std::ffi::OsStr::new("31-198")), &query),
            Some(vec![1])
        );
        // The errors of the argv, as the run reports them.
        for extra in [
            &["-subject_loc", "5-2", "-query_loc", "9-1"][..],
            &["-query_loc", "9-1", "-evalue", "0"],
            &["-subject_loc", "5-2", "-matrix", "IDENTITY"],
        ] {
            let argv = [&words[..5], extra].concat();
            let error = run(&argv, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
            assert_eq!(validate(&argv).unwrap_err(), error, "{extra:?}");
        }
        let error =
            validate(&[&words[..5], &["-subject_loc", "5-2", "-query_loc", "9-1"]].concat())
                .unwrap_err();
        assert!(error.contains("subject location"), "{error}");
        // The subject range comes before NCBI's errors of the options (IDENTITY, which BLASTX
        // does not take).
        let error =
            validate(&[&words[..5], &["-subject_loc", "5-2", "-matrix", "IDENTITY"]].concat())
                .unwrap_err();
        assert!(error.contains("subject location"), "{error}");
        // A subject range that starts past the subject's end is an error of the records, which
        // validate does not see.
        let past = [&words[..5], &["-subject_loc", "70-80"]].concat();
        assert!(validate(&past).is_ok());
        let error = run(&past, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
        assert!(error.contains("Invalid from coordinate"), "{error}");
    }

    // Observer events that break the formatter contract fail the run instead of
    // producing wrong ranges.
    #[test]
    fn unmatched_or_repeated_hsp_events_are_faults() {
        let mut unmatched = recorder();
        unmatched.hsp_end(1, 0);
        assert!(unmatched.fault.is_some());

        let mut crossed = recorder();
        crossed.hsp_begin(1, 0);
        crossed.hsp_end(1, 1);
        assert!(crossed.fault.is_some());

        let mut repeated = recorder();
        for _ in 0..2 {
            repeated.hsp_begin(1, 3);
            repeated.hsp_end(1, 3);
        }
        assert!(repeated.fault.is_some());

        let mut doubled = recorder();
        doubled.hsp_begin(1, 0);
        doubled.hsp_begin(1, 1);
        assert!(doubled.fault.is_some());

        let mut valid = recorder();
        valid.hsp_begin(1, 0);
        valid.hsp_end(1, 0);
        valid.subject_begin(0, 0);
        valid.subject_end(0, 0);
        valid.hsp_begin(0, 0);
        valid.hsp_end(0, 0);
        assert!(valid.fault.is_none());
        assert_eq!(
            (valid.out6.len(), valid.out0.len(), valid.out0_subject.len()),
            (1, 1, 1)
        );
    }

    // A format that failed alone leaves no trace in the HSP records (review 3 N5 of line RG T7):
    // the recorder forgets the ranges of its stream and an HSP that its report had begun.
    #[test]
    fn a_format_that_failed_alone_is_forgotten_by_the_recorder() {
        let mut recorder = recorder();
        recorder.hsp_begin(1, 0);
        recorder.hsp_end(1, 0);
        recorder.subject_begin(0, 0);
        recorder.subject_end(0, 0);
        recorder.hsp_begin(0, 0);
        recorder.hsp_end(0, 0);
        recorder.hsp_begin(0, 1);
        recorder.forget(0);
        assert!(recorder.open.iter().all(Option::is_none));
        assert!(recorder.out0.is_empty() && recorder.out0_subject.is_empty());
        assert_eq!(recorder.out6.len(), 1);
        recorder.forget(1);
        assert!(recorder.out6.is_empty());
        assert!(recorder.fault.is_none());
    }

    /// `letters` pseudo-random nucleotides (a fixed linear congruential generator, `seed`).
    fn random_nucleotides(letters: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..letters)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect()
    }

    /// `letters` pseudo-random residues of the 20 amino acids (the generator of
    /// `random_nucleotides`, `seed`).
    fn random_protein(letters: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..letters)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACDEFGHIKLMNPQRSTVWY"[((state >> 33) % 20) as usize]
            })
            .collect()
    }

    /// A coding sequence of the protein `letters` (one codon for each of the 20 amino acids).
    fn codons(letters: &[u8]) -> Vec<u8> {
        letters
            .iter()
            .flat_map(|&residue| {
                let codon: &[u8] = match residue {
                    b'A' => b"GCG",
                    b'C' => b"TGC",
                    b'D' => b"GAT",
                    b'E' => b"GAA",
                    b'F' => b"TTT",
                    b'G' => b"GGC",
                    b'H' => b"CAT",
                    b'I' => b"ATT",
                    b'K' => b"AAA",
                    b'L' => b"CTG",
                    b'M' => b"ATG",
                    b'N' => b"AAC",
                    b'P' => b"CCG",
                    b'Q' => b"CAG",
                    b'R' => b"CGT",
                    b'S' => b"AGC",
                    b'T' => b"ACC",
                    b'V' => b"GTG",
                    b'W' => b"TGG",
                    b'Y' => b"TAT",
                    _ => unreachable!("not one of the 20 amino acids"),
                };
                codon.iter().copied()
            })
            .collect()
    }

    // Line RG item T7-W (review 3 W3 of line RG T7, dev too): with `-sorthits`, NCBI BLAST+ 2.17.0
    // pads a description table whose last total is above 99999 bits without end
    // (showdefline.cpp:940-944), which the CLI refuses in outfmt 0 (`TotalScoreRefusal`, in the
    // middle of the report), while its `-outfmt 6`, field list and 5 are NCBI's (the 60 kb self hit
    // here: `1.108e+05` bits, NCBI's row). The queries before it write more than 1 MiB of outfmt
    // 0, so a chunk of stream 0 has reached the host when it fails: a run whose argv asks for
    // another format writes it, and stream 0 fails alone, discarded by the signal of stream 9,
    // its refusal once on stream 3 after the warnings, the HSP records without `out0`. The
    // default outfmt 0 fails the run with the refusal, as the CLI fails.
    #[test]
    fn a_sorthits_total_ncbi_pads_without_end_fails_stream_0_alone_after_its_chunks() {
        let mut query = Vec::new();
        let mut subject = Vec::new();
        for n in 0..7u64 {
            let letters = random_nucleotides(50_000, n + 1);
            query.push(FastaRecord::new(
                format!("Query_{}", n + 1),
                format!("q{n}").as_bytes(),
                &letters,
            ));
            subject.push(FastaRecord::new(
                format!("Subject_{}", n + 1),
                format!("s{n}").as_bytes(),
                &letters,
            ));
        }
        let big = random_nucleotides(60_000, 99);
        query.push(FastaRecord::new("Query_8", b"bq", &big));
        subject.push(FastaRecord::new("Subject_8", b"bs", &big));
        let words = ["blastn", "-query", "q", "-subject", "s", "-sorthits", "0"];
        let refusal = "-sorthits: the total score of the description table's last subject is wider than its column";
        let row = "bq\tbs\t100.000\t60000\t0\t0\t1\t60000\t1\t60000\t0.0\t1.108e+05\n";
        let six = streams(&[&words[..], &["-outfmt", "6"]].concat(), &query, &subject).unwrap();
        assert!(text(&six[&6]).ends_with(row), "{}", text(&six[&6]));
        assert_eq!(text(&six[&6]).lines().count(), 8);
        assert_eq!(discarded(&six), [0]);
        assert!(six[&0].is_empty());
        let error = text(&six[&3]);
        assert!(
            error.starts_with(
                "Warning: [blastn] The parameter -sorthits is ignored for output formats > 4.\n"
            ),
            "{error}"
        );
        assert_eq!(error.matches(refusal).count(), 1, "{error}");
        let records = text(&six[&STREAM_HITS]);
        assert_eq!(records.lines().count(), 8);
        assert!(
            records.contains("\"out6\":[0,")
                && !records.contains("\"out0\":[")
                && !records.contains("\"out0_subject\":["),
            "{records}"
        );
        for argv in ["6 qseqid sseqid evalue", "5"] {
            let custom =
                streams(&[&words[..], &["-outfmt", argv]].concat(), &query, &subject).unwrap();
            assert!(!custom[&STREAM_CUSTOM].is_empty(), "{argv}");
            assert_eq!(text(&custom[&6]), text(&six[&6]), "{argv}");
            assert!(custom[&0].is_empty(), "{argv}");
            assert_eq!(discarded(&custom), [0], "{argv}");
            assert_eq!(text(&custom[&3]).matches(refusal).count(), 1, "{argv}");
        }
        assert_eq!(
            text(&streams(&[&words[..], &["-outfmt", "6 qseqid sseqid evalue"]].concat(), &query, &subject).unwrap()[&STREAM_CUSTOM]),
            "q0\ts0\t0.0\nq1\ts1\t0.0\nq2\ts2\t0.0\nq3\ts3\t0.0\nq4\ts4\t0.0\nq5\ts5\t0.0\nq6\ts6\t0.0\nbq\tbs\t0.0\n"
        );
        let error = run(&words, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
        assert!(error.contains(refusal), "{error}");
    }

    // Line RG item T7-W (review 3 W4 of line RG T7, dev too): an NCBI error of outfmt 0 alone.
    // NCBI BLAST+ 2.17.0's BLASTX writes outfmt 0 of this ungapped one-hit run up to
    // `CSeq_align::GetSeqRange(): Row is empty` (Seq_align.cpp:226-229, exit 255), while its
    // `-outfmt 6` and field list complete (exit 0), and so does the CLI. A run whose argv asks
    // for them writes NCBI's rows, stream 0 fails alone (no chunk of it reached the host) with
    // NCBI's exception on stream 3 after the prefetch warning, each once; the default outfmt 0
    // fails the run with NCBI's exception, as the CLI fails.
    #[test]
    fn an_ncbi_error_of_outfmt_0_alone_fails_stream_0_alone_where_the_argv_asks_for_another_format()
    {
        let query = [FastaRecord::new(
            "Query_1",
            b"c00556_q",
            b"GACTCCCTCGAGCTTC",
        )];
        let subject = [FastaRecord::new("Subject_1", b"c00556_s", b"WRTIVII")];
        let words = [
            "blastx",
            "-query",
            "q",
            "-subject",
            "s",
            "-ungapped",
            "-window_size",
            "0",
            "-comp_based_stats",
            "0",
            "-seg",
            "no",
            "-evalue",
            "1e5",
            "-word_size",
            "3",
            "-threshold",
            "2",
        ];
        let exception = "CSeq_align::GetSeqRange(): Row is empty";
        let warning = "Warning: [blastx] Error pre-fetching sequence data\n";
        for (argv, stream, ncbi) in [
            (
                "6",
                6,
                "Query_1\tc00556_s\t0.000\t0\t0\t0\t0\t0\t0\t0\t3.7\t 4.2\n",
            ),
            (
                "6 qseqid sseqid evalue",
                STREAM_CUSTOM,
                "Query_1\tc00556_s\t3.7\n",
            ),
        ] {
            let written =
                streams(&[&words[..], &["-outfmt", argv]].concat(), &query, &subject).unwrap();
            assert_eq!(text(&written[&stream]), ncbi, "{argv}");
            assert!(written.get(&0).is_none_or(Vec::is_empty), "{argv}");
            // The host is told that stream 0 failed, also where no chunk of it was sent (line WQA's
            // sweep, D7, 2026-11-19).
            assert_eq!(discarded(&written), [0], "{argv}");
            let error = text(&written[&3]);
            assert!(error.starts_with(warning), "{argv}: {error}");
            assert_eq!(error.matches(warning).count(), 1, "{argv}: {error}");
            assert_eq!(error.matches(exception).count(), 1, "{argv}: {error}");
        }
        let error = run(&words, &query, &losat::api::QueryEnd::Input, &subject).unwrap_err();
        assert!(error.contains(exception), "{error}");
    }

    // Line RG item T7-W (decision T7-4, review 3 N5 of line RG T7): BLASTX writes its formats
    // query batch by query batch (10,002 query letters, NCBI's). The first batch (103 queries of
    // 99 letters against 320 subjects) writes more than 1 MiB of streams 6 and 7; the last query,
    // in the second batch, alone hits a subject whose title NCBI cannot decode, where NCBI's
    // tabular `saccver` stops (PD-LOSAT-NCBI-DEFECTS 1.37). With the default outfmt 0, or a field
    // list without `saccver`, the CLI completes; streams 6 and 7, which the argv does not ask
    // for, fail alone at that batch with a chunk sent: the signals of stream 9 discard them, the
    // run writes the others, the HSP records (which point into stream 6) are not written, and
    // stream 3 has the refusal once. The argv's own `-outfmt 6` fails the run, as the CLI fails.
    #[test]
    fn a_blastx_stream_that_fails_alone_at_a_later_batch_is_discarded() {
        let protein = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ";
        let other = b"CFDYQGKGFISMYSLRLTSNPFPEIFSRVERPMRGRLALIDVRHRVPYACGDESHIVAGRPARPTLTRS";
        let mut query: Vec<FastaRecord> = (0..103)
            .map(|n| {
                FastaRecord::new(
                    format!("Query_{}", n + 1),
                    format!("q{n}").as_bytes(),
                    &CODING[..99],
                )
            })
            .collect();
        query.push(FastaRecord::new("Query_104", b"qlast", &codons(other)));
        let mut subject: Vec<FastaRecord> = (0..320)
            .map(|n| {
                FastaRecord::new(
                    format!("Subject_{}", n + 1),
                    format!("s{n}").as_bytes(),
                    protein,
                )
            })
            .collect();
        subject.push(FastaRecord::new(
            "Subject_321",
            b"Subject_9 bad x\x81y title",
            other,
        ));
        let words = ["blastx", "-query", "q", "-subject", "s", "-evalue", "1e-10"];
        let refusal = "subject record 321 (Subject_9)";
        let plain = streams(&words, &query, &subject).unwrap();
        assert_eq!(discarded(&plain), [6, 7]);
        for stream in [6, 7, STREAM_HITS] {
            assert!(
                plain.get(&stream).is_none_or(Vec::is_empty),
                "stream {stream}"
            );
        }
        let report = text(&plain[&0]);
        assert!(report.len() > 1 << 20, "{}", report.len());
        assert!(
            report.contains("Query= qlast")
                && report.contains("Sequence with id Subject_321 no longer exists"),
            "{}",
            &report[report.len() - 2000..]
        );
        let last = &report[report.find("Query= qlast").unwrap()..];
        assert!(last.contains("Effective search space used:"), "{last}");
        let error = text(&plain[&3]);
        assert_eq!(error.matches(refusal).count(), 1, "{error}");
        let fields = [&words[..], &["-outfmt", "6 qseqid evalue"]].concat();
        let custom = streams(&fields, &query, &subject).unwrap();
        assert_eq!(discarded(&custom), [6, 7]);
        let rows = text(&custom[&STREAM_CUSTOM]);
        assert!(
            rows.starts_with("q0\t") && rows.contains("\nqlast\t"),
            "{}",
            &rows[rows.len().saturating_sub(300)..]
        );
        assert_eq!(rows.lines().count(), 103 * 320 + 1);
        assert_eq!(text(&custom[&0]), report);
        let error = run(
            &[&words[..], &["-outfmt", "6"]].concat(),
            &query,
            &losat::api::QueryEnd::Input,
            &subject,
        )
        .unwrap_err();
        assert!(error.contains(refusal), "{error}");
    }

    // Line RG item T7-W (review N5): a run of two searches (`-num_descriptions 5` gives format 0
    // a hit list size of its own, §7) in which stream 0, which the argv's `-outfmt 6` does not ask
    // for, fails alone in its search, after the report of a query with a warning (`warning`, which
    // the format 0 search writes before it fails, as the CLI's outfmt 0 does). Stream 3 holds
    // each text once: the CLI's standard error for `-outfmt 6` (the first search's), then the
    // refusal of stream 0 (`refusal`), without the warnings that the format 0 search repeats. So
    // it is that of the run's one search without `-num_descriptions` after the formatting
    // handler's line of `-num_descriptions`, which the CLI's `-outfmt 6` writes first.
    fn two_searches_write_each_text_once(
        program: &str,
        words: &[&str],
        query: &[FastaRecord],
        subject: &[FastaRecord],
        warning: &str,
        refusal: &str,
    ) {
        let one = streams(&[words, &["-outfmt", "6"]].concat(), query, subject).unwrap();
        let alone = text(&one[&3]);
        assert_eq!(alone.matches(warning).count(), 1, "{program}: {alone}");
        assert_eq!(alone.matches(refusal).count(), 1, "{program}: {alone}");
        assert!(
            alone.find(warning) < alone.find(refusal),
            "{program}: {alone}"
        );
        let two = streams(
            &[words, &["-outfmt", "6", "-num_descriptions", "5"]].concat(),
            query,
            subject,
        )
        .unwrap();
        assert_eq!(
            text(&two[&3]),
            format!(
                "Warning: [{program}] The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output\n{alone}"
            ),
            "{program}"
        );
        assert!(!two[&6].is_empty(), "{program}");
        assert_eq!(two[&6], one[&6], "{program}");
        for streams in [&one, &two] {
            assert!(streams.get(&0).is_none_or(Vec::is_empty), "{program}");
            assert_eq!(discarded(streams), [0], "{program}");
        }
    }

    /// The refusal of a `-sorthits` total above 99999 bits (`TotalScoreRefusal`).
    const SORTHITS_REFUSAL: &str =
        "-sorthits: the total score of the description table's last subject is wider than its column";

    // N5 in BLASTN: a query of N (no Karlin-Altschul parameters, NCBI's warning), then a 60 kb
    // self hit whose `-sorthits` total is above 99999 bits (review `bnw_nd5_f6`).
    #[test]
    fn blastn_two_searches_write_each_text_once() {
        let big = random_nucleotides(60_000, 5);
        let query = [
            FastaRecord::new("Query_1", b"nq1 ns", &[b'N'; 60]),
            FastaRecord::new("Query_2", b"bq big", &big),
        ];
        let subject = [FastaRecord::new("Subject_1", b"bs big", &big)];
        two_searches_write_each_text_once(
            "blastn",
            &["blastn", "-query", "q", "-subject", "s", "-sorthits", "0"],
            &query,
            &subject,
            "Warning: [blastn] Query_1 nq1 ns: Could not calculate ungapped Karlin-Altschul parameters",
            SORTHITS_REFUSAL,
        );
    }

    // N5 in TBLASTX: a query without letters (NCBI's warning; a query of N has none in TBLASTX's
    // search of `-subject`), then a 30 kb self hit whose total over the frames is above 99999 bits
    // (review `tx_*`).
    #[test]
    fn tblastx_two_searches_write_each_text_once() {
        let big = random_nucleotides(30_000, 6);
        let query = [
            FastaRecord::new("Query_1", b"eq1 empty", b""),
            FastaRecord::new("Query_2", b"bq big", &big),
        ];
        let subject = [FastaRecord::new("Subject_1", b"bs big", &big)];
        two_searches_write_each_text_once(
            "tblastx",
            &["tblastx", "-query", "q", "-subject", "s", "-sorthits", "0"],
            &query,
            &subject,
            "Warning: [tblastx] Query_1 eq1 empty: Sequence contains no data",
            SORTHITS_REFUSAL,
        );
    }

    // N5 in BLASTP: a query of X (no Karlin-Altschul parameters), then a 60,000-residue self hit
    // whose total is above 99999 bits (review `bp_nd5_f6`).
    #[test]
    fn blastp_two_searches_write_each_text_once() {
        let big = random_protein(60_000, 7);
        let query = [
            FastaRecord::new("Query_1", b"xq1 xs", &[b'X'; 60]),
            FastaRecord::new("Query_2", b"bq big", &big),
        ];
        let subject = [FastaRecord::new("Subject_1", b"bs big", &big)];
        two_searches_write_each_text_once(
            "blastp",
            &["blastp", "-query", "q", "-subject", "s", "-sorthits", "0"],
            &query,
            &subject,
            "Warning: [blastp] Query_1 xq1 xs: Could not calculate ungapped Karlin-Altschul parameters",
            SORTHITS_REFUSAL,
        );
    }

    // N5 in TBLASTN: the queries of BLASTP's case against a coding sequence of the big protein.
    #[test]
    fn tblastn_two_searches_write_each_text_once() {
        let big = random_protein(60_000, 8);
        let query = [
            FastaRecord::new("Query_1", b"xq1 xs", &[b'X'; 60]),
            FastaRecord::new("Query_2", b"bq big", &big),
        ];
        let subject = [FastaRecord::new("Subject_1", b"bs big", &codons(&big))];
        two_searches_write_each_text_once(
            "tblastn",
            &["tblastn", "-query", "q", "-subject", "s", "-sorthits", "0"],
            &query,
            &subject,
            "Warning: [tblastn] Query_1 xq1 xs: Could not calculate ungapped Karlin-Altschul parameters",
            SORTHITS_REFUSAL,
        );
    }

    // N5 in BLASTX: the ungapped one-hit run of NCBI's `Row is empty` in outfmt 0 alone (W4,
    // `an_ncbi_error_of_outfmt_0_alone_...`), whose prefetch warning every search writes.
    #[test]
    fn blastx_two_searches_write_each_text_once() {
        let query = [FastaRecord::new(
            "Query_1",
            b"c00556_q",
            b"GACTCCCTCGAGCTTC",
        )];
        let subject = [FastaRecord::new("Subject_1", b"c00556_s", b"WRTIVII")];
        two_searches_write_each_text_once(
            "blastx",
            &[
                "blastx",
                "-query",
                "q",
                "-subject",
                "s",
                "-ungapped",
                "-window_size",
                "0",
                "-comp_based_stats",
                "0",
                "-seg",
                "no",
                "-evalue",
                "1e5",
                "-word_size",
                "3",
                "-threshold",
                "2",
            ],
            &query,
            &subject,
            "Warning: [blastx] Error pre-fetching sequence data\n",
            "CSeq_align::GetSeqRange(): Row is empty",
        );
    }

    /// The streams of a run that failed, as a host keeps them (docs/web/abi_v2.md §5: stream 3
    /// and stream 8), and its error.
    fn failed_streams(
        words: &[&str],
        query: &[FastaRecord],
        subject: &[FastaRecord],
    ) -> (String, std::collections::BTreeMap<u32, Vec<u8>>) {
        crate::emit::EMITTED.with(|emitted| emitted.borrow_mut().clear());
        let error = run(words, query, &losat::api::QueryEnd::Input, subject).unwrap_err();
        let mut streams = std::collections::BTreeMap::new();
        crate::emit::EMITTED.with(|emitted| {
            for (stream, bytes) in emitted.borrow().iter() {
                streams
                    .entry(*stream)
                    .or_insert_with(Vec::new)
                    .extend_from_slice(bytes);
            }
        });
        (error, streams)
    }

    /// Stream 3 of a run as text.
    fn diagnostics_of(streams: &std::collections::BTreeMap<u32, Vec<u8>>) -> String {
        text(
            streams
                .get(&STREAM_DIAGNOSTICS)
                .map_or(&[][..], Vec::as_slice),
        )
    }

    /// Whether the bytes are the CLI's prolog of a pairwise report (outfmt 0 to 4, 11, and the program's name and
    /// version on its first line) and nothing after it but blank lines and the reference.
    fn prolog_only(bytes: &[u8]) -> bool {
        let text = String::from_utf8_lossy(bytes);
        let first = text.lines().next().unwrap_or("");
        first.contains("BLAST")
            && first.contains(losat::web_adapter::REPORT_VERSION)
            && !text.contains("Query=")
            && !text.contains('\t')
    }

    /// A failed run sent no hits and no output with rows: what it keeps of its outputs is the CLI's file before the
    /// error (outfmt 0's prolog, `a_failed_run_sends_the_output_of_its_argv_and_drops_the_others`), and the records
    /// (each fits in one chunk) are none.
    fn assert_outputs_empty(streams: &std::collections::BTreeMap<u32, Vec<u8>>, label: &str) {
        for stream in [0, 6, 7, STREAM_CUSTOM, STREAM_HITS] {
            let bytes = streams.get(&stream).map_or(&[][..], Vec::as_slice);
            assert!(
                bytes.is_empty() || ([0, STREAM_CUSTOM].contains(&stream) && prolog_only(bytes)),
                "{label}: stream {stream}"
            );
        }
    }

    // A format 0 search that ends a run of two searches (§7) gives its warnings, and only its: the
    // run's stream 3 and error are the CLI's standard error for the argv. BLASTP with
    // `-num_alignments 0` and a subject whose first record has no residues and whose second is
    // titled `, ,`: the tabular search fails alone before it reads the subject, and the format 0
    // search, the argv's own `-outfmt 0`, reads the subject, warns, and fails the run on the title,
    // where the CLI writes the same warning before its error (native LOSAT and NCBI BLAST+ 2.17.0
    // 2026-11-19, `fix-diagnostics/native/twosearch`; NCBI crashes after the warning). The tabular
    // search's warning and error, which the CLI writes for `-outfmt 6` only, are not stream 3's
    // (line WQA's sweep, D3 review M1, 2026-11-19).
    #[test]
    fn a_format_0_search_that_fails_the_run_gives_its_warnings() {
        let query = [FastaRecord::new(
            "Query_1",
            b"p",
            b"MKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMAWRSTY",
        )];
        let subject = [
            FastaRecord::new("Subject_1", b"empty record", b""),
            FastaRecord::new(
                "Subject_2",
                b", ,",
                b"MKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMA",
            ),
        ];
        let words = [
            "blastp",
            "-query",
            "q",
            "-subject",
            "s",
            "-num_alignments",
            "0",
        ];
        let (error, streams) = failed_streams(&words, &query, &subject);
        assert!(
            error.contains("subject record 2 has a title of punctuation"),
            "{error}"
        );
        assert_eq!(
            diagnostics_of(&streams),
            "Warning: [blastp] Subject_1 empty record: Subject sequence contains no data\n"
        );
        assert_outputs_empty(&streams, "blastp");
    }

    // A run that fails gives on stream 3 what the command line writes to its standard error before
    // its error: NCBI BLAST+ 2.17.0 ($NCBI_BIN, 2026-11-19, line WQA's sweep D3) writes the subject's
    // warning and then the error, exit status 3, for a subject whose one record has no residues
    // (BLASTN and TBLASTX; the other programs search it with the warning), and the error alone for
    // a query whose one record has none. TBLASTX's `-logfile` file holds both, and stream 3 nothing.
    #[test]
    fn a_failed_run_gives_the_warnings_before_its_error() {
        let letters = random_nucleotides(300, 11);
        let query = [FastaRecord::new("Query_1", b"q1", &letters)];
        let empty = [FastaRecord::new("Subject_1", b"empty record", b"")];
        for program in ["blastn", "tblastx"] {
            let words = [program, "-query", "q", "-subject", "s", "-outfmt", "6"];
            let (error, streams) = failed_streams(&words, &query, &empty);
            assert_eq!(
                error.trim_end(),
                ncbi_in_losat_words("BLAST engine error: The average subject length is too short"),
                "{program}"
            );
            assert_eq!(
                text(streams.get(&STREAM_DIAGNOSTICS).map_or(&[][..], Vec::as_slice)),
                format!("Warning: [{program}] Subject_1 empty record: Subject sequence contains no data\n"),
                "{program}"
            );
            for stream in [0, 6, 7, STREAM_HITS] {
                assert!(
                    streams.get(&stream).is_none_or(Vec::is_empty),
                    "{program} {stream}"
                );
            }
            let (error, streams) = failed_streams(&words, &empty, &query);
            assert_eq!(
                error.trim_end(),
                ncbi_in_losat_words("BLAST engine error: Warning: Sequence contains no data"),
                "{program}"
            );
            assert!(
                streams.get(&STREAM_DIAGNOSTICS).is_none_or(Vec::is_empty),
                "{program}: {:?}",
                streams.get(&STREAM_DIAGNOSTICS).map(|bytes| text(bytes))
            );
        }
        let words = [
            "tblastx", "-query", "q", "-subject", "s", "-outfmt", "6", "-logfile", "run.log",
        ];
        let (_, streams) = failed_streams(&words, &query, &empty);
        assert!(streams.get(&STREAM_DIAGNOSTICS).is_none_or(Vec::is_empty));
        let archive = crate::zip::read(&streams[&crate::emit::STREAM_RUN_FILES]);
        assert_eq!(
            archive
                .iter()
                .map(|(name, bytes)| (text(name), text(bytes)))
                .collect::<Vec<_>>(),
            [(
                "run.log".to_string(),
                ncbi_in_losat_words("Warning: [tblastx] Subject_1 empty record: Subject sequence contains no data\nBLAST engine error: The average subject length is too short\n").to_string()
            )]
        );
    }

    /// The inputs of a run of `program` that fails with `-query_loc 500-600` (a query of 300
    /// letters, `BLAST engine error: Empty CBlastQueryVector`) and completes without it: the
    /// query again as the subject where both have the same letters, so that the search has hits.
    fn ranged_inputs(program: &str) -> (Vec<FastaRecord>, Vec<FastaRecord>) {
        let protein_query = matches!(program, "blastp" | "tblastn");
        let protein_subject = matches!(program, "blastp" | "blastx");
        let letters = |protein: bool| match protein {
            true => random_protein(300, 21),
            false => random_nucleotides(300, 21),
        };
        (
            vec![FastaRecord::new("Query_1", b"q1", &letters(protein_query))],
            vec![FastaRecord::new(
                "Subject_1",
                b"s1",
                &letters(protein_subject),
            )],
        )
    }

    // The formatting handler warns of `-num_descriptions`, `-line_length` and `-sorthits` for
    // output formats above 4 and of `-sorthsps` for formats other than 0 (blast_args.cpp:2947-2989):
    // a failed run whose argv asks for format 0 (no `-outfmt` too) to 4 gives none of the lines of
    // the formats that it does not ask for, as the CLI's standard error for the argv has none
    // (native LOSAT 16bef7acd and NCBI BLAST+ 2.17.0, review r2 probe `result5.txt`; line WQA's sweep,
    // D3 review r2 M1, 2026-11-19: the tabular formats came first in the search, and their lines
    // were stream 3's).
    #[test]
    fn a_failed_run_gives_the_handler_warnings_of_the_argv_s_format() {
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            let (query, subject) = ranged_inputs(program);
            let warning = |text: &str| format!("Warning: [{program}] {text}\n");
            let line_length =
                warning("The parameter -line_length is not applicable for output formats > 4 .");
            let sorthits = warning("The parameter -sorthits is ignored for output formats > 4.");
            let sorthsps = warning("The parameter -sorthsps is ignored for output formats != 0.");
            let descriptions = warning("The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output");
            let cases: [(&[&str], &[&str], String); 13] = [
                (&["-line_length", "80"], &[], String::new()),
                (&["-line_length", "80"], &["-outfmt", "0"], String::new()),
                (&["-line_length", "80"], &["-outfmt", "1"], String::new()),
                (
                    &["-line_length", "80"],
                    &["-outfmt", "6"],
                    line_length.clone(),
                ),
                (
                    &["-line_length", "80"],
                    &["-outfmt", "6 qseqid"],
                    line_length,
                ),
                (&["-sorthits", "1"], &[], String::new()),
                (&["-sorthits", "1"], &["-outfmt", "7"], sorthits),
                (&["-sorthsps", "1"], &[], String::new()),
                (&["-sorthsps", "1"], &["-outfmt", "1"], sorthsps.clone()),
                (&["-sorthsps", "1"], &["-outfmt", "6"], sorthsps),
                (
                    &["-num_descriptions", "5", "-num_alignments", "5"],
                    &[],
                    String::new(),
                ),
                (
                    &["-num_descriptions", "5", "-num_alignments", "5"],
                    &["-outfmt", "4"],
                    String::new(),
                ),
                (
                    &["-num_descriptions", "5", "-num_alignments", "5"],
                    &["-outfmt", "6"],
                    descriptions,
                ),
            ];
            for (options, outfmt, expected) in cases {
                let words = [
                    &[
                        program,
                        "-query",
                        "q",
                        "-subject",
                        "s",
                        "-query_loc",
                        "500-600",
                    ][..],
                    options,
                    outfmt,
                ]
                .concat();
                let label = format!("{words:?}");
                let (error, streams) = failed_streams(&words, &query, &subject);
                assert_eq!(
                    error,
                    ncbi_in_losat_words("BLAST engine error: Empty CBlastQueryVector\n"),
                    "{label}"
                );
                assert_eq!(diagnostics_of(&streams), expected, "{label}");
                assert_outputs_empty(&streams, &label);
            }
        }
    }

    // A run that completes gives the CLI's standard error for the argv first, then what the
    // standard error of the other formats adds (abi_v2.md §5): the formatting handler's lines that
    // the CLI writes for the tabular formats only come after the subject's warning where the argv
    // asks for format 0, and before it where it asks for a tabular format, as the CLI writes them
    // (line WQA's sweep, D3 review r2 M1, 2026-11-19).
    #[test]
    fn a_completed_run_gives_the_argv_s_standard_error_first() {
        let letters = random_nucleotides(300, 23);
        let query = [FastaRecord::new("Query_1", b"q1", &letters)];
        let subject = [
            FastaRecord::new("Subject_1", b"empty", b""),
            FastaRecord::new("Subject_2", b"s2", &letters[..200]),
        ];
        let empty = "Warning: [blastn] Subject_1 empty: Subject sequence contains no data\n";
        let line_length = "Warning: [blastn] The parameter -line_length is not applicable for output formats > 4 .\n";
        let sorthsps =
            "Warning: [blastn] The parameter -sorthsps is ignored for output formats != 0.\n";
        let cases: [(&[&str], String); 5] = [
            (&[], format!("{empty}{line_length}{sorthsps}")),
            (&["-outfmt", "0"], format!("{empty}{line_length}{sorthsps}")),
            (&["-outfmt", "1"], format!("{sorthsps}{empty}{line_length}")),
            (&["-outfmt", "6"], format!("{line_length}{sorthsps}{empty}")),
            (&["-outfmt", "7"], format!("{line_length}{sorthsps}{empty}")),
        ];
        for (outfmt, expected) in cases {
            let words = [
                &[
                    "blastn",
                    "-query",
                    "q",
                    "-subject",
                    "s",
                    "-line_length",
                    "80",
                    "-sorthsps",
                    "1",
                ][..],
                outfmt,
            ]
            .concat();
            let streams = streams(&words, &query, &subject).unwrap();
            assert_eq!(diagnostics_of(&streams), expected, "{words:?}");
            assert!(
                !streams[&0].is_empty() && !streams[&6].is_empty(),
                "{words:?}"
            );
        }
    }

    // TBLASTX writes the run's messages to the file of `-logfile`, which holds what the CLI's file
    // holds for the argv, in a run that completes and in one that fails (native LOSAT 16bef7acd,
    // review r2 probe `result6-base.txt`): not the formatting handler's lines of the formats that
    // the argv does not ask for, nor the warnings and the error of the other search of a run of two
    // searches (`-num_alignments 0`: the tabular formats' search fails alone, and the CLI's file for
    // the argv, outfmt 0, is empty). Those go to stream 3, which the file does not hold for the argv
    // (line WQA's sweep, D3 review r2, 2026-11-19).
    #[test]
    fn tblastx_s_file_of_logfile_holds_what_the_cli_s_file_holds_for_the_argv() {
        let (query, subject) = ranged_inputs("tblastx");
        let file_of = |streams: &std::collections::BTreeMap<u32, Vec<u8>>| {
            crate::zip::read(&streams[&crate::emit::STREAM_RUN_FILES])
                .into_iter()
                .map(|(name, bytes)| (text(&name), text(&bytes)))
                .collect::<Vec<_>>()
        };
        let log = |text: &str| vec![("run.log".to_string(), text.to_string())];
        let line_length = "Warning: [tblastx] The parameter -line_length is not applicable for output formats > 4 .\n";
        let words = [
            "tblastx",
            "-query",
            "q",
            "-subject",
            "s",
            "-line_length",
            "80",
            "-logfile",
            "run.log",
        ];
        // Completed: the file of outfmt 0 has no line; stream 3 has the tabular formats' line.
        let completed = streams(&words, &query, &subject).unwrap();
        assert_eq!(file_of(&completed), log(""));
        assert_eq!(diagnostics_of(&completed), line_length);
        // With `-outfmt 6` the line is the file's.
        let tabular = [&words[..], &["-outfmt", "6"]].concat();
        let completed = streams(&tabular, &query, &subject).unwrap();
        assert_eq!(file_of(&completed), log(line_length));
        assert_eq!(diagnostics_of(&completed), "");
        // Failed: the file has the error alone (with `-outfmt 6`, the line and the error).
        let ranged = [&words[..], &["-query_loc", "500-600"]].concat();
        let (error, failed) = failed_streams(&ranged, &query, &subject);
        assert_eq!(
            error,
            ncbi_in_losat_words("BLAST engine error: Empty CBlastQueryVector\n")
        );
        assert_eq!(file_of(&failed), log(&error));
        assert_eq!(diagnostics_of(&failed), "");
        let ranged = [&ranged[..], &["-outfmt", "6"]].concat();
        let (_, failed) = failed_streams(&ranged, &query, &subject);
        assert_eq!(file_of(&failed), log(&format!("{line_length}{error}")));
        assert_eq!(diagnostics_of(&failed), "");
        // Two searches: the tabular formats' search fails alone; the run completes as the CLI's run
        // of outfmt 0, whose file is empty, and stream 3 has the other search's warning and error.
        let words = [
            "tblastx",
            "-query",
            "q",
            "-subject",
            "s",
            "-num_alignments",
            "0",
            "-logfile",
            "run.log",
        ];
        let completed = streams(&words, &query, &subject).unwrap();
        assert_eq!(file_of(&completed), log(""));
        assert_eq!(
            diagnostics_of(&completed),
            ncbi_in_losat_words("Warning: [tblastx] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n")
        );
        assert!(!completed[&0].is_empty());
        assert!(completed.get(&6).is_none_or(Vec::is_empty));
        // Where the run's messages go nowhere (`-logfile ""`, `-logfile /dev/null`), so do the
        // other formats' lines: the CLI writes none of them (line WQA's sweep, D3 review r3,
        // 2026-11-19: stream 3 had them).
        for nowhere in ["", "/dev/null"] {
            for options in [&["-line_length", "80"][..], &["-num_alignments", "0"]] {
                let words = [
                    &["tblastx", "-query", "q", "-subject", "s"][..],
                    options,
                    &["-logfile", nowhere],
                ]
                .concat();
                let completed = streams(&words, &query, &subject).unwrap();
                assert_eq!(diagnostics_of(&completed), "", "{words:?}");
                assert!(!completed[&0].is_empty(), "{words:?}");
            }
        }
    }

    // `validate` checks the options of BLASTP, TBLASTN and TBLASTX with the argv's own format,
    // whose hit list size the CLI checks before the e-value (`No hits are being saved` for
    // `-num_alignments 0` with a tabular format), and gives the run's error, the CLI's (line WQA's
    // sweep, D3 review r2, 2026-11-19: it gave the e-value's error, as for format 0).
    #[test]
    fn validate_checks_the_options_with_the_argv_s_format() {
        let no_hits =
            ncbi_in_losat_words("BLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n");
        for program in ["blastp", "tblastn", "tblastx", "blastn", "blastx"] {
            let (query, subject) = ranged_inputs(program);
            for outfmt in [
                &[][..],
                &["-outfmt", "0"],
                &["-outfmt", "6"],
                &["-outfmt", "7"],
                &["-outfmt", "6 qseqid"],
                &["-outfmt", "5"],
            ] {
                let words = [
                    &[
                        program,
                        "-query",
                        "q",
                        "-subject",
                        "s",
                        "-num_alignments",
                        "0",
                        "-evalue",
                        "0",
                    ][..],
                    outfmt,
                ]
                .concat();
                let validated = validate(&words).unwrap_err();
                let (error, _) = failed_streams(&words, &query, &subject);
                assert_eq!(validated, error, "{words:?}");
                if !outfmt.is_empty() && outfmt[1] != "0" {
                    assert_eq!(validated, no_hits, "{words:?}");
                } else {
                    assert_ne!(validated, no_hits, "{words:?}");
                }
            }
        }
    }

    // A run whose own -outfmt is structured (5, 11, 13; BLASTN's 17) writes outfmt 0, 6 and 7 too,
    // whose CLI runs write the queries' result warnings, which the structured format's run does not
    // (blast_format.cpp:1431-1452): stream 3 is the CLI's standard error for the argv, then those
    // warnings, which the CLI's standard error for the other formats adds (line WQA's sweep, D3
    // review r3 M2, 2026-11-19: they were lost; native LOSAT 16bef7acd, review r3 probe
    // `result-warn`).
    #[test]
    fn a_structured_outfmt_s_run_adds_the_query_warnings_of_the_other_formats() {
        for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
            let (ranged, subject) = ranged_inputs(program);
            let mut query = vec![FastaRecord::new("Query_1", b"qempty no residues", b"")];
            query.push(FastaRecord::new("Query_2", b"q1", ranged[0].seq()));
            let warning = format!(
                "Warning: [{program}] Query_1 qempty no residues: Sequence contains no data \n"
            );
            let words = [program, "-query", "q", "-subject", "s", "-sorthsps", "1"];
            let tabular =
                streams(&[&words[..], &["-outfmt", "6"]].concat(), &query, &subject).unwrap();
            let tabular = diagnostics_of(&tabular);
            assert_eq!(
                tabular.matches(&warning).count(),
                1,
                "{program}: {tabular:?}"
            );
            // The CLI's standard error for a structured format: that of outfmt 6 (the formatting
            // handler warns alike for both) without the query's warning.
            let argv_s = tabular.replace(&warning, "");
            let structured: &[&str] = match program {
                "blastn" => &["5", "11", "13", "17"],
                _ => &["5", "11", "13"],
            };
            for outfmt in structured {
                let completed = streams(
                    &[&words[..], &["-outfmt", outfmt]].concat(),
                    &query,
                    &subject,
                )
                .unwrap();
                assert_eq!(
                    diagnostics_of(&completed),
                    format!("{argv_s}{warning}"),
                    "{program} -outfmt {outfmt}"
                );
            }
        }
    }

    // A format that the argv does not ask for and whose report stops after a query's warning gave
    // that warning before it stopped, and the formats after it give their own: stream 3 keeps it
    // once. BLASTX with HSPs of length 0 (`-ungapped -window_size 0`, the manifest's
    // `x2g.blastx.c0054` inputs): NCBI's outfmt 6, 7 and 0 throw at the second query's alignments,
    // after the first query's warning, and outfmt 5 completes (native LOSAT f80c94dd3). Stream 3 is
    // the CLI's standard error for `-outfmt 5`, then the warning once and the three refusals
    // (review r1 of the later formats' warnings, M2, 2026-11-20: one receiver took the warning of
    // each format that stopped, and stream 3 showed it three times).
    #[test]
    fn formats_that_stop_after_a_query_warning_add_it_once() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../LOSAT/tests/fasta/outfmt0"
        );
        let sequence = |name: &str| -> Vec<u8> {
            let text = std::fs::read_to_string(format!("{root}/{name}")).unwrap();
            text.lines().skip(1).flat_map(str::bytes).collect()
        };
        let query = vec![
            FastaRecord::new("Query_1", b"qempty no residues", b""),
            FastaRecord::new("Query_2", b"q1", &sequence("x2g_c0054_q.fna")),
        ];
        let subject = vec![FastaRecord::new(
            "Subject_1",
            b"s1",
            &sequence("x2g_c0054_s.faa"),
        )];
        let options = [
            "-ungapped",
            "-comp_based_stats",
            "0",
            "-window_size",
            "0",
            "-word_size",
            "5",
            "-threshold",
            "1",
            "-evalue",
            "1e5",
            "-seg",
            "no",
        ];
        let words = [&["blastx", "-query", "q", "-subject", "s"][..], &options].concat();
        let warning = "Warning: [blastx] Query_1 qempty no residues: Sequence contains no data \n";
        let written = streams(&[&words[..], &["-outfmt", "5"]].concat(), &query, &subject).unwrap();
        assert_eq!(discarded(&written), [6, 7, 0]);
        // The CLI's standard error for `-outfmt 5` (NCBI's post of the Iteration that its XML
        // leaves out, in LOSAT's words), then what outfmt 6 adds (the warning and its refusal;
        // outfmt 7 adds the same two lines, which stream 3 holds) and what outfmt 0 adds (its
        // refusal).
        let own = "Error: [blastx] Failed s_BlastXMLAddIteration Q(0/1: CreateDensegFromStdseg(): \
                   Std-seg segment lengths not accurate!\n\n";
        let tabular = "Error: CreateDensegFromStdseg(): Std-seg segment lengths not accurate!\n\n";
        let pairwise = "Error: CSeq_align::GetSeqRange(): Row is empty\n\n";
        assert_eq!(
            diagnostics_of(&written),
            format!("{own}{warning}{tabular}{pairwise}")
        );
        // The argv's own outfmt 6 fails the run after the warning, as the CLI fails.
        let error = run(
            &[&words[..], &["-outfmt", "6"]].concat(),
            &query,
            &losat::api::QueryEnd::Input,
            &subject,
        )
        .unwrap_err();
        assert!(error.contains(tabular.trim_end()), "{error}");
        // A run of two searches (`-num_descriptions 1` gives format 0 a hit list of 250 subjects
        // and leaves the tabular formats theirs of 500, §7): 250 subjects score above the one
        // with the HSP of length 0, so format 0 does not show it and completes; the tabular
        // search's formats, which fail alone, each stop at it, after the warning. They share
        // their search's diagnostics: outfmt 6 wrote the warning there before it stopped, and
        // outfmt 7 does not write it again (`blastx::local::write_run`; line WQA, 2026-11-20:
        // stream 3 showed it twice).
        let mut subjects: Vec<FastaRecord> = (1..=250)
            .map(|number| FastaRecord::new(format!("Subject_{number}"), b"sbest", b"CVFNR"))
            .collect();
        subjects.push(FastaRecord::new(
            "Subject_251",
            b"s1",
            &sequence("x2g_c0054_s.faa"),
        ));
        let one = ["-num_descriptions", "1"];
        let written = streams(&[&words[..], &one].concat(), &query, &subjects).unwrap();
        assert_eq!(discarded(&written), [6, 7]);
        // The CLI's standard error for the argv, then what outfmt 6 adds: its handler's line
        // and its refusal (outfmt 7 adds the same lines).
        let ignored = "Warning: [blastx] The parameter -num_descriptions is ignored for output \
                       formats > 4 . Use -max_target_seqs to control output\n";
        assert_eq!(
            diagnostics_of(&written),
            format!("{warning}{ignored}{tabular}")
        );
    }
}
