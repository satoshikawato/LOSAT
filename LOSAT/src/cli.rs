//! Public CLI v2 boundary: LOSAT's command line, `losat <program> <words>`, whose programs read
//! their arguments with NCBI's argument reader (`reader`) into their own options (`Commands`).

use crate::blastinput::ncbi_environment::DiagTarget;
use crate::ncbi_app::args::Refusal;
use crate::wording::{DATABASE_ERROR, ENGINE_ERROR, OPTIONS_ERROR, QUERY_OPTIONS_ERROR};
use std::{ffi::OsString, fmt};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
use crate::algorithm::{
    blastn, blastp, blastx, deltablast, psiblast, rpsblast, rpstblastn, tblastn, tblastx,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46,51,75
// LOSAT retains its program subcommand; search parameter names are NCBI names.
pub use crate::utils::entry_names::{entry_argv, ENTRY_NAMES};

// What the `losat` binary runs; the crate's root re-exports `main_entry` for it.
pub(crate) mod main_entry;

// NCBI's argument reader of the programs that read their arguments into their options.
pub(crate) mod reader;
pub use reader::ArgvError;

// NCBI's buffers of a report stream (`NcbiStdio`, `NcbiBuffer`).
mod stdio_buffer;
pub use stdio_buffer::{NcbiBuffer, NcbiStdio, NCBI_BUFSIZ};

// The process's standard output with no buffer of Rust's, under the command line's writers.
mod standard_output;
pub(crate) use standard_output::standard_output;

/// The words that ask for `install_shims`: `-install-shims` is how the top level's `-help` shows
/// it.
pub const INSTALL_SHIMS_WORDS: [&str; 2] = ["--install-shims", "-install-shims"];

/// `losat --install-shims [DIR]`: makes the `ENTRY_NAMES` in `DIR` (default: the directory of
/// this executable) and returns their paths. Each is a symbolic link to this executable, relative
/// when it lies in `DIR` so that the directory can move; Windows, where a link needs a privilege,
/// gets a copy. A name made by an earlier LOSAT (a link to a `losat`, or a copy of one) is
/// replaced. Any other file of that name is an error and nothing is made: it may be NCBI's
/// program. No NCBI counterpart.
pub fn install_shims(dir: Option<&std::ffi::OsStr>) -> Result<Vec<std::path::PathBuf>, String> {
    let exe = std::env::current_exe()
        .and_then(std::fs::canonicalize)
        .map_err(|error| format!("cannot find this executable: {error}"))?;
    let exe_dir = exe
        .parent()
        .map(std::path::PathBuf::from)
        .unwrap_or_default();
    let dir = dir.map_or_else(|| exe_dir.clone(), std::path::PathBuf::from);
    std::fs::create_dir_all(&dir)
        .map_err(|error| format!("cannot create {}: {error}", dir.display()))?;
    let dir = std::fs::canonicalize(&dir)
        .map_err(|error| format!("cannot read {}: {error}", dir.display()))?;
    let mut missing = Vec::new();
    let mut stale = Vec::new();
    let mut paths = Vec::new();
    let mut taken = Vec::new();
    for name in ENTRY_NAMES {
        let path = dir.join(format!("{name}{}", std::env::consts::EXE_SUFFIX));
        if std::fs::symlink_metadata(&path).is_err() {
            missing.push(path.clone());
        } else if is_this_executable(&path, &exe) {
            // Already made by this executable: kept.
        } else if is_earlier_shim(&path) {
            stale.push(path.clone());
        } else {
            taken.push(path.display().to_string());
        }
        paths.push(path);
    }
    if !taken.is_empty() {
        return Err(format!(
            "these files exist and are not this executable; nothing was installed: {}",
            taken.join(" ")
        ));
    }
    for path in &stale {
        std::fs::remove_file(path)
            .map_err(|error| format!("cannot replace {}: {error}", path.display()))?;
    }
    for path in missing.into_iter().chain(stale) {
        make_shim(&exe, exe_dir == dir, &path)
            .map_err(|error| format!("cannot make {}: {error}", path.display()))?;
    }
    Ok(paths)
}

/// Whether `path` was made by `install_shims` of another `losat` (an upgrade, a move): a link
/// whose target is named `losat`, which may be gone, or a file with LOSAT's `-version` text (a
/// Windows copy, a copy of an older build).
fn is_earlier_shim(path: &std::path::Path) -> bool {
    let losat = format!("losat{}", std::env::consts::EXE_SUFFIX);
    if let Ok(target) = std::fs::read_link(path) {
        return target
            .file_name()
            .is_some_and(|name| name == losat.as_str());
    }
    const MARKER: &[u8] = b" Package: LOSAT ";
    std::fs::read(path)
        .is_ok_and(|bytes| bytes.windows(MARKER.len()).any(|window| window == MARKER))
}

/// Whether `path` is a link to `exe`, or `exe` itself.
#[cfg(not(windows))]
fn is_this_executable(path: &std::path::Path, exe: &std::path::Path) -> bool {
    std::fs::canonicalize(path).is_ok_and(|target| target == exe)
}

/// Whether `path` is a copy of `exe`.
#[cfg(windows)]
fn is_this_executable(path: &std::path::Path, exe: &std::path::Path) -> bool {
    match (std::fs::read(path), std::fs::read(exe)) {
        (Ok(shim), Ok(this)) => shim == this,
        _ => false,
    }
}

#[cfg(unix)]
fn make_shim(exe: &std::path::Path, same_dir: bool, path: &std::path::Path) -> std::io::Result<()> {
    match exe.file_name().filter(|_| same_dir) {
        Some(file_name) => std::os::unix::fs::symlink(file_name, path),
        None => std::os::unix::fs::symlink(exe, path),
    }
}

#[cfg(windows)]
fn make_shim(exe: &std::path::Path, _: bool, path: &std::path::Path) -> std::io::Result<()> {
    std::fs::copy(exe, path).map(|_| ())
}

#[cfg(not(any(unix, windows)))]
fn make_shim(_: &std::path::Path, _: bool, _: &std::path::Path) -> std::io::Result<()> {
    Err(std::io::ErrorKind::Unsupported.into())
}

/// What a program that `main_entry` runs before `cli::reader` reads the command line (the
/// programs of `ENTRY_NAMES` that are not a `Commands` variant: they read their words with NCBI's
/// argument reader, `ncbi_app`) says about itself in the top level's `-help`: NCBI's DESCRIPTION
/// line of the program (its first line; `convert2blastmask`'s continues in NCBI's text).
fn entry_program_description(name: &str) -> &'static str {
    match name {
        "blastdbcmd" => "BLAST database client",
        "dustmasker" => "Low complexity region masker based on Symmetric DUST algorithm",
        "segmasker" => "Low complexity region masker based on the SEG algorithm",
        "windowmasker" => "Window based sequence masker",
        "convert2blastmask" => {
            "Convert masking information in lower-case masked FASTA input to file formats suitable for makeblastdb"
        }
        "makembindex" => "Create a BLAST database index.",
        "makeprofiledb" => "Application to create databases for rpsblast, cobalt and deltablast",
        "blastdb_aliastool" => "Application to create BLAST database aliases",
        "blastdbcheck" => "BLAST database integrity and validity checking application",
        "cleanup-blastdb-volumes.py" => "Remove needless BLAST database volumes.",
        _ => "",
    }
}

/// The `after_help` of the top level: the programs that `losat NAME` runs and `Commands` does not
/// list (derived from `ENTRY_NAMES` and `is_program`, so that a new name shows up), and the shims.
/// No NCBI counterpart.
fn top_level_after_help() -> &'static str {
    static TEXT: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    TEXT.get_or_init(|| {
        let others: Vec<&str> = ENTRY_NAMES
            .into_iter()
            .filter(|name| !is_program(name))
            .collect();
        let width = others.iter().map(|name| name.len()).max().unwrap_or(0);
        let mut text = String::new();
        if !others.is_empty() {
            text.push_str("Other programs (losat <program> -help shows each one's options):\n");
            for name in others {
                text.push_str(&format!(
                    "  {name:<width$}  {}\n",
                    entry_program_description(name)
                ));
            }
            text.push('\n');
        }
        text.push_str("To call the programs by NCBI's names (blastn, makeblastdb, ...), link them with:\n  losat -install-shims [DIR]   (default: the directory of this executable)");
        text
    })
}

// No NCBI counterpart: `losat <program>` names the program (NCBI has one executable per program).
/// LOSAT's command line, `losat <program> <words>`: the options of its program.
#[derive(Debug)]
pub struct Cli {
    pub command: Commands,
}

#[cfg(test)]
impl Cli {
    /// BLASTN's arguments of a command line that names blastn (the tests' argv).
    pub(crate) fn blastn(self) -> blastn::BlastnArgs {
        match self.command {
            Commands::Blastn(args) => args,
            command => panic!("not blastn: {command:?}"),
        }
    }
}

/// The options of the program that a command line runs, one variant for each program on the
/// reader's table (`reader::tool`), in its order: each read from NCBI's arguments by its type's
/// `from_args` (`reader::options`, line RG DECISIONS K4i-3-1).
#[derive(Debug, Clone)]
pub enum Commands {
    Blastx(blastx::BlastxArgs),
    Blastn(blastn::BlastnArgs),
    Blastp(blastp::BlastpArgs),
    Tblastx(tblastx::TblastxArgs),
    Tblastn(tblastn::TblastnArgs),
    Psiblast(psiblast::PsiblastArgs),
    Deltablast(deltablast::DeltablastArgs),
    Rpsblast(rpsblast::RpsblastArgs),
    Rpstblastn(rpstblastn::RpstblastnArgs),
    Makedb(crate::writedb::makeblastdb::MakedbArgs),
    Rpsbproc(crate::rpsbproc::RpsbprocArgs),
}

/// Reads `argv` (`losat <program> <words>`) into the program's arguments for a caller that runs no
/// NCBI application start (the engine's own argv: BLASTN's repeat search; the tests): NCBI's
/// reader of the program, LOSAT's refusals, the program's options (`reader::read_argv`). The
/// toolkit's words that the start reads (`-logfile`, `-conffile`, `-dryrun`) are refused
/// (`toolkit_file_word_refusal`; makeblastdb keeps its own `-logfile` argument), and so are the
/// versions (`version_word_refusal`) and `-xmlhelp` (`xml_help_refusal`).
pub fn try_parse_from<I, S>(argv: I) -> Result<Cli, ArgvError>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    read_cli(argv, reader::Caller::NoStart, &DiagTarget::Stderr)
}

/// `try_parse_from` for a caller whose calls run an NCBI application start of their own and return
/// reports: the library API (`api::search`). `-logfile`, `-conffile` and `-dryrun` are the
/// start's, whose `-logfile` target (`target`) takes the posts of an argument error, as the
/// command line's; the versions and `-xmlhelp`, whose texts the library has no place for, are
/// refused (`version_word_refusal`, `xml_help_refusal`).
pub(crate) fn try_parse_with_own_start<I, S>(argv: I, target: &DiagTarget) -> Result<Cli, ArgvError>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    read_cli(argv, reader::Caller::OwnStart, target)
}

/// `try_parse_from` for LOSAT Web's ABI v2, whose runs make their own NCBI application start
/// (`web_adapter::start_run`) and give the host the text that the command line prints in place of
/// a run (stream 11; line WQA's D14, 2026-11-20): `-logfile`, `-conffile` and `-dryrun` are the
/// start's, whose `-logfile` target (`target`) takes the posts of an argument error, as the
/// command line's, and `-xmlhelp` prints the XML, as on the command line (`ArgvError::Printed`).
/// The versions are printed before the arguments are read (`ncbi_version_request`,
/// `api::start::CallerStart::version`), so a version word that reaches this reader is refused
/// (`version_word_refusal`).
pub(crate) fn try_parse_with_printing_start<I, S>(
    argv: I,
    target: &DiagTarget,
) -> Result<Cli, ArgvError>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    read_cli(argv, reader::Caller::PrintingStart, target)
}

/// `reader::read_argv` of `argv` for `caller`, with the posts of an argument error going to
/// `target`, as `Cli`.
fn read_cli<I, S>(argv: I, caller: reader::Caller, target: &DiagTarget) -> Result<Cli, ArgvError>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let argv: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    reader::read_argv(&argv, caller, target).map(|command| Cli { command })
}

/// What LOSAT's top level writes for `argv` (`losat <words>`, whose first word names no program on
/// a reader, `ENTRY_NAMES`): its help or version on the standard output (exit status 0), or its
/// error on the standard error (exit status 2). The command line's entry (`main_entry`) writes it,
/// and the readers of an argv give it (`reader::read_argv`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TopLevel {
    /// The text, which ends with a line end.
    pub(crate) text: String,
    /// Whether it is an error, for the standard error.
    pub(crate) error: bool,
    /// The exit status.
    pub(crate) status: i32,
}

// No NCBI counterpart: NCBI has one executable per program. LOSAT's top level has no option of its
// own: its first word is its help (`-h`, `-help`, `--help`), its version (`--version`), or an
// error; the words after it are not read. Its texts are those of LOSAT 0.2.0, frozen
// (`top_level_tests`, line RG DECISIONS K4i-3-1 Q5).
/// What the top level writes for `argv` (`TopLevel`): without a word, its help as an error.
pub(crate) fn top_level(argv: &[OsString]) -> TopLevel {
    let Some(first) = argv.get(1) else {
        return TopLevel {
            text: top_level_help(usage_name(argv)),
            error: true,
            status: 2,
        };
    };
    let (text, error) = match first.to_str() {
        Some("-h" | "-help" | "--help") => (top_level_help(usage_name(argv)), false),
        Some("--version") => (format!("losat {}\n", env!("CARGO_PKG_VERSION")), false),
        _ => (
            format!(
                "error: unknown option or argument '{}'; use -help to list the options\n",
                first.to_string_lossy()
            ),
            true,
        ),
    };
    TopLevel {
        text,
        error,
        status: if error { 2 } else { 0 },
    }
}

/// The arguments that a program's command line gave, by NCBI's names, each with the text of its
/// value (empty for a flag), in the order of the command line: what the import of a search
/// strategy compares with NCBI's defaults (`report::archive::read::import::NcbiArguments`, which
/// keeps them by name, as NCBI's `CArgs` holds them). Read from NCBI's argument reader
/// (`ncbi_app::args::Args::given_values`) only for a command line that gives
/// `-import_search_strategy`: the only reader. Every read of an argv makes its own (a default
/// argument list read later does not replace them).
#[derive(Debug, Clone, Default)]
pub struct GivenArguments {
    list: Vec<(String, String)>,
}

impl GivenArguments {
    /// The arguments that `args` gave, when they give `-import_search_strategy`; else none.
    pub(crate) fn from_args(args: &crate::ncbi_app::args::Args) -> Self {
        if !args.given("import_search_strategy") {
            return Self::default();
        }
        let list = args
            .given_values()
            .map(|(name, value)| {
                let text = value.map(|value| value.to_string_lossy().into_owned());
                (name.to_string(), text.unwrap_or_default())
            })
            .collect();
        Self { list }
    }

    /// The arguments given, as `(name, text)`.
    pub(crate) fn list(&self) -> &[(String, String)] {
        &self.list
    }
}

/// What the top level's help says of each program whose options are LOSAT's types: the programs on
/// the reader's table (`reader::tool`), in its order
/// (`reader::options_tests::one_list_of_programs`).
const PROGRAMS: [(&str, &str); 11] = [
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
    (
        "blastx",
        "Translated nucleotide query vs protein subject (blastx)",
    ),
    (
        "blastn",
        "Pairwise nucleotide alignment (megablast [default], blastn, dc-megablast, blastn-short)",
    ),
    ("blastp", "Pairwise protein alignment (blastp)"),
    (
        "tblastx",
        "Pairwise 6-frame translated nucleotide alignment (tblastx)",
    ),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:47-52
    (
        "tblastn",
        "Protein query vs translated nucleotide subject (tblastn)",
    ),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:47-51
    ("psiblast", "Position-Specific Iterated BLAST (psiblast)"),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/deltablast_args.cpp:44-52
    (
        "deltablast",
        "Domain Enhanced Lookup Time Accelerated BLAST (deltablast)",
    ),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:77-79
    ("rpsblast", "Reverse Position Specific BLAST (rpsblast)"),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:48-51
    (
        "rpstblastn",
        "Translated Reverse Position Specific BLAST (rpstblastn)",
    ),
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:217-219
    (
        "makeblastdb",
        "makeblastdb: create a BLAST database from a FASTA file (format version 5 by default, 4 \
         with -blastdb_version 4)",
    ),
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:64-65
    (
        "rpsbproc",
        "Process output from offline rpsblast to annotate sequence with conserved domain \
         information (rpsbproc; its domain hits, sites and motifs of protein and nucleotide \
         queries)",
    ),
];

/// Whether `name` is a program whose options are LOSAT's types (`Commands`): one on the reader's
/// table (`reader::tool`).
pub(crate) fn is_program(name: &str) -> bool {
    reader::tool(name).is_some()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3815-3826
// The help opens with NCBI's two headings, `USAGE` and `DESCRIPTION`: tools that call BLAST+ look
// for them (BRIG takes a `blastn -help` that prints `USAGE` as BLAST+ found).
/// The top level's help, whose usage names the executable `executable` (`usage_name`): its usage,
/// its description, the programs (`PROGRAMS`, and a `help` that the top level of LOSAT 0.2.0
/// listed), its two flags, and the other programs and the shims (`top_level_after_help`).
fn top_level_help(executable: &str) -> String {
    const HELP: (&str, &str) = (
        "help",
        "Print this message or the help of the given subcommand(s)",
    );
    let listed = || PROGRAMS.iter().chain([&HELP]);
    let width = listed().map(|(name, _)| name.len()).max().unwrap_or(0);
    let mut text = format!(
        "USAGE\n  {executable} <COMMAND>\n\nDESCRIPTION\n   Pure-Rust pairwise alignment that \
         matches NCBI BLAST+ on LOSAT's test fixtures, with genetic-code exceptions\n\nCommands:\n",
    );
    for (name, about) in listed() {
        text += &format!("  {name:<width$}  {about}\n");
    }
    text += "\nOptions:\n  -help     Print help\n  --version  Print version\n\n\n\n";
    text += top_level_after_help();
    text.push('\n');
    text
}

// No NCBI counterpart: the name that the top level of LOSAT 0.2.0 wrote, frozen in
// `top_level_tests` for each form of `argv[0]`.
/// The name that the top level's usage gives the executable: the file name of `argv[0]`, as it was
/// called (`losat.exe`, the name of a link), when `argv[0]` has one that is UTF-8; else `losat`
/// (no `argv[0]`, an empty one, `/`, a file name that is not UTF-8).
fn usage_name(argv: &[OsString]) -> &str {
    argv.first()
        .and_then(|argv0| std::path::Path::new(argv0).file_name())
        .and_then(std::ffi::OsStr::to_str)
        .unwrap_or("losat")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:266-270
// NCBI keeps the bytes of an option's value.
/// LOSAT's refusal of the value of `-name` that is not UTF-8 in `program` (its name in LOSAT's
/// refusals, `MAKEBLASTDB`): NCBI's arguments keep the value's bytes. The one text of it, for the
/// programs on `cli::reader` (`reader::refusal_after_read`).
pub(crate) fn non_utf8_value_refusal(program: &str, name: &str) -> String {
    format!(
        "the value of -{name} is not UTF-8; NCBI BLAST+ reads the bytes of an option's value, \
         which is not supported by LOSAT's {program}"
    )
}

/// LOSAT's refusal of the run of `-name` in `program` (its name in LOSAT's refusals, `BLASTP`):
/// what the argument is, from its row of the refusals (`ncbi_app::args::refusal`, `Refusal::Run`).
pub(crate) fn run_refusal(program: &str, name: &str) -> String {
    let what = match crate::ncbi_app::args::refusal(&program.to_ascii_lowercase(), name) {
        Some(Refusal::Run(what)) => Some(what),
        _ => None,
    };
    debug_assert!(what.is_some(), "{program} -{name} is no run refusal");
    match what {
        Some(what) => format!("-{name} ({what}) is not supported by LOSAT's {program}"),
        None => format!("-{name} is not supported by LOSAT's {program}"),
    }
}

/// The text that `-version` or a full version (`-version-full`, `-version-full-xml`,
/// `-version-full-json`) prints for `argv` (`losat <program> ... -version ...`), or `None` when
/// `argv` has no such request.
///
/// NCBI's application reads these words before its arguments (`ncbi_environment::preparse`,
/// ncbiapp.cpp:926-1001): anywhere before `--`, in an option's value place too
/// (`blastn -out -version` prints the version) but not as the value of `-logfile` or
/// `-conffile`, and also when other arguments are invalid; the first one is printed.
/// `-version` prints the two lines that the flags `fVersionInfo` and `fPackageShort` select: the
/// program name and version, then the package. NCBI's second line names its own package and build
/// date; LOSAT's names LOSAT and the BLAST+ release its output follows. The full versions are
/// NCBI's forms with LOSAT's values (`ncbi_version::print`, line RG item E3). `makeblastdb` prints
/// `makeblastdb`, as NCBI's program does. `main` prints the text on stdout and exits 0, before any
/// argument is parsed. The top level's `losat --version` is not read here (`top_level`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:972-977
pub fn ncbi_version_request<I, S>(argv: I) -> Option<String>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let tokens: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    // A program's name (`is_program`).
    let name = tokens.get(1)?.to_str().filter(|word| is_program(word))?;
    // NCBI's reading of the words before the arguments (`ncbi_environment::preparse`): not
    // after `--`, nor as the value of `-logfile` or `-conffile` (`blastn -logfile -version`
    // logs to the file `-version`).
    match crate::blastinput::ncbi_environment::preparse(&tokens[2..]) {
        crate::blastinput::ncbi_environment::Preparsed::Version(request) => {
            // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:388
            // Every program prints LOSAT's version (NCBI's rpsbproc its own 0.5.1, the BLAST+
            // programs 2.17.0+; the Owner 2026-11-21 「Q1=B」).
            Some(crate::ncbi_version::print(request, name))
        }
        _ => None,
    }
}

/// The first `-logfile`, `-conffile` or `-dryrun` that NCBI's application reads in `tokens` (with
/// its value; `ncbi_environment::preparse`), as written: the words that a caller without an
/// application start refuses (`reader::Caller::NoStart`). makeblastdb's `-logfile` is its own
/// argument (`MakedbArgs::logfile`), which its run opens.
pub(crate) fn toolkit_file_word<'a>(tokens: &'a [OsString], program: &str) -> Option<&'a str> {
    use crate::blastinput::ncbi_environment::{log_file, preparse, Preparsed};
    let (conffile, dryrun) = match preparse(tokens) {
        Preparsed::Words {
            conffile, dryrun, ..
        } => (conffile.is_some(), dryrun),
        _ => (false, false),
    };
    if program != "makeblastdb" && log_file(tokens).is_some() {
        return toolkit_word_as_written(tokens, "-logfile");
    }
    if conffile {
        return toolkit_word_as_written(tokens, "-conffile");
    }
    if dryrun {
        return toolkit_word_as_written(tokens, "-dryrun");
    }
    None
}

/// The first word of `tokens` that is the toolkit option `name` (`-conffile` also as
/// `-conffile=FILE`), as written.
fn toolkit_word_as_written<'a>(tokens: &'a [OsString], name: &str) -> Option<&'a str> {
    tokens
        .iter()
        .filter_map(|token| token.to_str())
        .find(|text| match text.strip_prefix(name) {
            Some(rest) => rest.is_empty() || (name != "-dryrun" && rest.starts_with('=')),
            None => false,
        })
}

// No NCBI counterpart: NCBI has no library API or Web engine. The application start of the
// command line (`cli::main_entry`, `ncbi_environment::ApplicationStart`), of each library call
// (`api::start`) and of each run of LOSAT Web's ABI v2 (`web_adapter::start_run`) reads `-logfile`
// and `-conffile` before the arguments, and each ends a `-dryrun` after them; ABI v1 runs none, so
// it does not take them (line RG items E1w, E1w-3a).
/// LOSAT's refusal of the toolkit option `word` (`-logfile`, `-conffile`, `-dryrun`) in the
/// program of `tool` (named as its other refusals name it, `Tool::refusal_name`) when the reader's
/// caller runs no application start (`reader::Caller::NoStart`).
pub(crate) fn toolkit_file_word_refusal(tool: &crate::ncbi_app::Tool, word: &str) -> String {
    let name = toolkit_option_name(word);
    let (program, shown) = (tool.program, tool.refusal_name);
    let places = match library_part(program) {
        LibraryPart::Search if WEB_V2_PROGRAMS.contains(&program) => {
            "the command line, in the library API and in the Web engine's ABI v2"
        }
        LibraryPart::Search | LibraryPart::Tool => "the command line and in the library API",
        LibraryPart::Absent => "the command line",
    };
    format!(
        "the option -{name} runs on {places} of LOSAT's {shown} only: the Web engine's ABI v1 runs no NCBI application start, which reads it"
    )
}

/// The programs that LOSAT Web's ABI v2 runs (`web/adapter/src/run.rs` `Program`, whose tests
/// check this list), whose runs take the toolkit's words (`toolkit_file_word_refusal`).
pub const WEB_V2_PROGRAMS: &[&str] = &["blastn", "blastp", "blastx", "tblastn", "tblastx"];

/// What LOSAT's library API runs of the command line's `program`.
enum LibraryPart {
    /// A search (`api::Program`), whose results are reports: it prints no version.
    Search,
    /// A tool (`api::run_makeblastdb`), which prints the versions to its output.
    Tool,
    /// Nothing: rpsblast, rpstblastn and rpsbproc.
    Absent,
}

/// The part of the library API that runs `program` (`LibraryPart`). The other tools of the
/// library (`api::run_blastdbcmd` and the maskers) have a parser of their own (`ncbi_app`).
fn library_part(program: &str) -> LibraryPart {
    if program == "makeblastdb" {
        LibraryPart::Tool
    } else if crate::api::Program::parse(program).is_ok() {
        LibraryPart::Search
    } else {
        LibraryPart::Absent
    }
}

/// The name of the toolkit option `word` (`-conffile=FILE` is `conffile`).
fn toolkit_option_name(word: &str) -> &str {
    let name = word.trim_start_matches('-');
    name.split_once('=').map_or(name, |(key, _)| key)
}

// No NCBI counterpart: NCBI's application prints the versions before its arguments
// (ncbiapp.cpp:926-1001, cited at `ncbi_environment::preparse`). The library API of the search
// programs returns reports and has no place for the text (line RG DECISIONS E3-2; the library's
// tools print it to their output, `ncbi_app::run_with`, `api::run_makeblastdb`), nor has the Web
// engine's ABI v1. A run of the Web engine's ABI v2 prints it on its stream of the command line's
// standard output (`api::start::CallerStart::version`; line WQA's D14, 2026-11-20).
/// LOSAT's refusal of the version word `word` (`-version` or a full version) in the program of
/// `tool` (named as its other refusals name it, `Tool::refusal_name`) where the reader's caller
/// prints no version: the library API of the search programs and the Web engine's ABI v1
/// (`reader::read_argv`).
pub(crate) fn version_word_refusal(tool: &crate::ncbi_app::Tool, word: &str) -> String {
    let name = toolkit_option_name(word);
    let (program, shown) = (tool.program, tool.refusal_name);
    let (places, refusers) = match library_part(program) {
        LibraryPart::Search if WEB_V2_PROGRAMS.contains(&program) => {
            ("the command line and in LOSAT Web", "the library API does")
        }
        LibraryPart::Search => ("the command line", "the library API and the Web engine do"),
        LibraryPart::Tool => (
            "the command line and in the library API",
            "the Web engine does",
        ),
        LibraryPart::Absent => ("the command line", "the Web engine does"),
    };
    format!(
        "the option -{name} runs on {places} of LOSAT's {shown} only: {refusers} not write NCBI's version text"
    )
}

// No NCBI counterpart: NCBI has no library API. NCBI's application prints the XML of `-xmlhelp`
// (ncbiapp.cpp:685-690, cited at `ncbi_xmlhelp`); LOSAT's library refuses it (line RG DECISIONS
// E2x-1, the Owner's choice of 2026-11-13). The command line's entries reach the XML
// (`xml_help_text`, `ncbi_app::run`), and so does a run of the Web engine's ABI v2, which prints it
// on its stream of the command line's standard output (`try_parse_with_printing_start`; line WQA's D14 and Q-XML=A,
// the Owner's choice of 2026-11-20, which takes the XML into the Web engine where E2x-1 and E2x-2
// left it out); the library and the Web engine's ABI v1 refuse it.
/// LOSAT's refusal of `-xmlhelp` where no command line's entry prints it: the library API and the
/// Web engine's ABI v1 of `shown` (`ncbi_app::Tool::refusal_name`).
pub(crate) fn xml_help_refusal(shown: &str) -> String {
    let web = WEB_V2_PROGRAMS
        .iter()
        .any(|program| program.eq_ignore_ascii_case(shown));
    match web {
        true => format!(
            "the option -xmlhelp runs on the command line and in LOSAT Web of LOSAT's {shown} only: the library API does not write NCBI's XML description of the arguments"
        ),
        false => format!(
            "the option -xmlhelp runs on the command line of LOSAT's {shown} only: the library API and the Web engine do not write NCBI's XML description of the arguments"
        ),
    }
}

/// What `losat <program> -xmlhelp` prints: NCBI's XML of `program` without the arguments that
/// LOSAT refuses (`ncbi_app::args::refusal`, by name or at the run; `ncbi_xmlhelp::print`).
pub(crate) fn xml_help_text(program: &str) -> Option<String> {
    let refused = |name: &str| crate::ncbi_app::args::refusal(program, name).is_some();
    crate::ncbi_xmlhelp::print(program, refused)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:78-88
// The standard options of an NCBI C++ Toolkit application, which NCBI's BLAST+ programs
// accept without listing them in -help. `-logfile` and `-conffile` run in the BLAST programs and
// makeblastdb (line RG item E1), `-dryrun` in them and blastdbcmd (item E1b), the three in
// rpsbproc (item E1c, its own `DryRun`: `rpsbproc::dry_run`), `-help-full` in every program (item
// E2, NCBI's reader: `reader`), the full versions in every program (item E3,
// `ncbi_version_request`), `-xmlhelp` in every program on the command line (item E2x,
// `xml_help_text`).

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-230
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
#[derive(Debug, Clone)]
pub struct NativeError {
    pub exit: i32,
    /// The text NCBI's `CATCH_ALL` writes, as bytes: a name of the command line or of a file
    /// keeps the bytes NCBI's `std::string` keeps (a name that is not UTF-8: probes R1, R3 and R8
    /// of line LIB-B).
    pub message: Vec<u8>,
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:178-180
/// The text of `message`; a byte that is not UTF-8 reads as U+FFFD (`message` keeps it).
impl fmt::Display for NativeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.message))
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
impl std::error::Error for NativeError {}

/// Whether `error` is `NativeError::silent_end` (including a writer's `StrategyWriteFault`).
pub(crate) fn is_silent_end(error: &anyhow::Error) -> bool {
    NativeError::caught(error).is_some_and(|native| native.exit == 0 && native.message.is_empty())
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:946
// NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:342
/// The C-string prefix of `bytes` up to its first NUL byte (`const char *` output of
/// `CException::what()`).
pub(crate) fn c_string(bytes: &[u8]) -> &[u8] {
    bytes.split(|&byte| byte == 0).next().unwrap_or_default()
}

impl NativeError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
    /// The end of a run by a `CIOException` other than `eFlush`, which `CATCH_ALL` ends without a
    /// post and with the run's status, 0 (a `CEofException` of a PSSM file that ends after its
    /// header; the failed write of the object manager's message after the diagnostic stream
    /// failed). The programs that write an error archive write it of the messages saved before
    /// (`is_silent_end`).
    pub(crate) fn silent_end() -> Self {
        Self {
            exit: 0,
            message: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
    // NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:946
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:342
    /// `CATCH_ALL`'s report of an exception of NCBI's C++ Toolkit (`CException`) that ends a
    /// BLAST program: `Error: `, the exception's `what()` as a C string (ending at the first NUL
    /// byte, dropping any bytes and trailing line feed after it; without a NUL, `what()` ends
    /// with a line feed) and the post's line feed, exit status 255 (`BLAST_UNKNOWN_ERROR`); a
    /// `CSeqDBException`, which `CATCH_ALL` catches before (its branch comes first),
    /// `BLAST Database error: ` and its message, exit status 2 (`BLAST_DATABASE_ERROR`). Every
    /// place of a `CSeqDBException` that LOSAT raises reports the code `eInvalid`
    /// (`seqdb::throw_site`), so the texts of that branch for `eTooManyOpenFiles` and
    /// `eMemoryMappingFailure` do not arise here (the reader's errors of those codes are
    /// `common::database::database_error`'s). The one producer of these reports: an exception of
    /// NCBI's that LOSAT raises where NCBI throws it becomes this error.
    pub(crate) fn exception(e: &crate::ncbi_exception::NcbiException) -> NativeError {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-229
        // NCBI reference (598d8ae6): c++/include/algo/blast/api/remote_blast.hpp:63-64
        // `CATCH_ALL`'s clause of `CBlastException`, which also catches its subclass
        // `CRemoteBlastException`.
        if e.kind().starts_with("CBlastException::")
            || e.kind().starts_with("CRemoteBlastException::")
        {
            let msg = e.get_msg();
            let has = |text: &[u8]| msg.windows(text.len()).any(|part| part == text);
            return if e.kind() == "CBlastException::eInvalidOptions" {
                NativeError {
                    exit: 1,
                    message: [OPTIONS_ERROR.as_bytes(), b" ", msg, b"\n"].concat(),
                }
            } else if has(b"Out of memory") || has(b"Failed to allocate") {
                // The message names an allocation that failed: exit status 4
                // (`BLAST_OUT_OF_MEMORY`), the report of `out_of_memory` with the message.
                NativeError {
                    exit: 4,
                    message: [crate::wording::OUT_OF_MEMORY.as_bytes(), b": ", msg, b"\n"].concat(),
                }
            } else {
                NativeError {
                    exit: 3,
                    message: [ENGINE_ERROR.as_bytes(), b" ", msg, b"\n"].concat(),
                }
            };
        }
        if e.kind().starts_with("CSeqDBException::") {
            return NativeError {
                exit: 2,
                message: [DATABASE_ERROR.as_bytes(), b" ", e.get_msg(), b"\n"].concat(),
            };
        }
        let what = e.what();
        let what = c_string(&what);
        NativeError {
            exit: 255,
            message: [b"Error: ".as_slice(), what, b"\n"].concat(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-384 (cited at
    // `report::serial::encoder::NullElement`)
    /// The end of a run whose ASN.1 text stream threw at an element that holds no object
    /// (`encoder::NullElement`, the application `application`): the stream's post of its first
    /// failure, then `CATCH_ALL`'s report of the exception (`exception`, exit status 255).
    pub(crate) fn null_element(
        application: &str,
        thrown: &crate::report::serial::encoder::NullElement,
    ) -> NativeError {
        let mut error = Self::exception(&thrown.exception());
        let mut message = thrown.post(Some(application));
        message.append(&mut error.message);
        error.message = message;
        error
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:260-263
    /// `CATCH_ALL`'s report of a `std::exception` (one that is neither a `CException` nor
    /// `bad_alloc` nor `ios::failure`) whose `what()` is `what`: exit status 255.
    pub(crate) fn std_exception(what: &str) -> NativeError {
        NativeError {
            exit: 255,
            message: format!("Error: {what}\n").into(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:256-259
    /// `CATCH_ALL`'s report of a `std::bad_alloc` (a string or byte block of a PSSM file in ASN.1
    /// binary whose length passes every allocation, `report::serial::asn_input::reserve`): exit
    /// status 4 (`BLAST_OUT_OF_MEMORY`).
    pub(crate) fn out_of_memory() -> NativeError {
        NativeError {
            exit: 4,
            message: [crate::wording::OUT_OF_MEMORY.as_bytes(), b"\n"].concat(),
        }
    }

    /// The error of a report's writer as the run ends with it: an exception of NCBI's that the
    /// writer raised mid-report is `CATCH_ALL`'s report of it (`of_writer_exception`), after what
    /// the report wrote before it; any other error is the writer's.
    pub(crate) fn of_writer(error: std::io::Error) -> anyhow::Error {
        match NativeError::of_writer_exception::<dyn std::error::Error>(&error) {
            Some(native) => native.into(),
            None => error.into(),
        }
    }

    /// `CATCH_ALL`'s report (`exception`) of the exception of NCBI's that the error `error` of a
    /// report's writer carries, where it carries one (`NcbiException::of_io_error`: raised
    /// mid-report by the `-html` registries, `report::formatter_registry`, by a subject's load
    /// that ends the run, `report::x_pairwise_par::end_here`, or by the merge of a query-anchored
    /// section, `report::query_anchored::write_section` and
    /// `report::query_anchored::write_unreadable_section`), whichever program wrote the report. The
    /// one mapping of such an error: at the programs' report boundaries (`of_writer`, BLASTX's
    /// `native::report_error`) and, through `caught`, at the command line's end
    /// (`exit_on_native_error`) and in the library's message (`api::error::chain_bytes`), so the
    /// library's bytes are the command line's.
    pub(crate) fn of_writer_exception<E: CaughtError + ?Sized>(error: &E) -> Option<NativeError> {
        let io = error.find::<std::io::Error>()?;
        if crate::report::ostream::is_strategy_write_fault(io) {
            return Some(NativeError::silent_end());
        }
        crate::ncbi_exception::NcbiException::of_io_error(io).map(NativeError::exception)
    }

    /// The report with which `CATCH_ALL` ends a run whose error (or one error of whose chain) is
    /// `error`, where it has one: the report of an exception of NCBI's that a report's writer
    /// raised (`of_writer_exception`), or the `NativeError` itself.
    pub(crate) fn caught<E: CaughtError + ?Sized>(
        error: &E,
    ) -> Option<std::borrow::Cow<'_, NativeError>> {
        match NativeError::of_writer_exception(error) {
            Some(native) => Some(std::borrow::Cow::Owned(native)),
            None => error.find::<NativeError>().map(std::borrow::Cow::Borrowed),
        }
    }
}

/// An error that `NativeError::caught` reads: a run's error (anyhow's `downcast_ref`, which also
/// finds the error under a context) or one error of a chain (`dyn Error`'s).
pub(crate) trait CaughtError {
    /// The error as a `T`, where it is one.
    fn find<T: std::error::Error + Send + Sync + 'static>(&self) -> Option<&T>;
}

impl CaughtError for anyhow::Error {
    fn find<T: std::error::Error + Send + Sync + 'static>(&self) -> Option<&T> {
        self.downcast_ref()
    }
}

impl CaughtError for dyn std::error::Error + 'static {
    fn find<T: std::error::Error + Send + Sync + 'static>(&self) -> Option<&T> {
        self.downcast_ref()
    }
}

pub use crate::report::ostream::{insert, insert_streamed, CliStderr, NcbiOstream, ReportStream};

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
pub fn inaccessible(name: &str, path: &std::path::Path) -> anyhow::Error {
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:95-99 (`s_ArgExptMsg`, cited in
    // `utils::os_bytes`): the name's bytes (probe R1 of line LIB-B: BLASTX's `-query`,
    // `-subject` and `-out` of a name that is not UTF-8).
    let mut message =
        format!("Command line argument error: Argument \"{name}\". File is not accessible")
            .into_bytes();
    if !path.as_os_str().is_empty() {
        message.extend_from_slice(b":  `");
        message.extend_from_slice(&crate::utils::os_bytes::path_bytes(path));
        message.push(b'\'');
    }
    message.push(b'\n');
    NativeError { exit: 1, message }.into()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3478-3480
/// The `-out` file `path` of a search, created where NCBI's argument handling opens it (after the
/// subjects are read and the query is opened, before the options are processed): `inaccessible`
/// when it cannot be made. `program` names the program whose file names must be UTF-8
/// (`input_files::check_utf8_file_name`); `None` for BLASTX, which has no such check. The command
/// line's programs (`open_out_file`) and the library's `run_to_file` (`api::local_blast::OutFile`)
/// make the file here.
pub(crate) fn create_out_file(
    path: &std::path::Path,
    program: Option<&str>,
) -> anyhow::Result<std::fs::File> {
    if let Some(program) = program {
        crate::blastinput::input_files::check_utf8_file_name(path, "out", program)?;
    }
    std::fs::File::create(path).map_err(|_| inaccessible("out", path))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
/// Writes the text of the run's `NativeError` to standard error and exits with its code; a text
/// that the run already wrote (`post_native_error`) is not written again. A failed write is
/// dropped, as NCBI's post drops it (`eprint!` would abort: tblastx `-evalue 0 -logfile /dev/full`
/// exits with 1 and writes nothing in NCBI BLAST+ 2.17.0; line WQA's D11 round 2, 2026-11-20). An
/// error of a report's writer that carries an exception of NCBI's is that exception's report
/// (`NativeError::caught`), whichever program wrote it.
///
/// The text of `CATCH_ALL` is a post on `cerr`, which a failed standard error does not take
/// (`report::ostream::cerr_post`); the C++ runtime's text before an abort (exit -6, `std::terminate`)
/// is written to C's `stderr` whatever `cerr`'s state.
pub fn exit_on_native_error(error: &anyhow::Error) {
    if let Some(e) = NativeError::caught(error) {
        if posted(&e.message) {
        } else if e.exit < 0 {
            let mut stderr = std::io::stderr();
            let _ = std::io::Write::write_all(&mut stderr, catch_all_text(&e.message));
            let _ = std::io::Write::flush(&mut stderr);
        } else {
            let _ = crate::report::ostream::cerr_post(catch_all_text(&e.message));
        }
        if e.exit < 0 {
            std::process::abort();
        }
        std::process::exit(e.exit);
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-176 (`CATCH_ALL`'s one clause
// of two posts; every other clause posts once)
/// The text that `CATCH_ALL` writes for the `NativeError` text `message`: the text, or, after the
/// main thread's diagnostic stream failed (`ncbi_exception::main_diag_failbit`), the line feed of
/// each of its `LOG_POST`s (two for the query/options error of
/// `blastinput::app::options_error`, one for the others; NCBI BLAST+ 2.17.0, line RG D3-1 `om/`
/// cases k02, k03, k06 and `om/mute/` e05, e11). A text without a post (the end of the PSSM file
/// directly in its root object, exit status 0) stays empty.
pub(crate) fn catch_all_text(message: &[u8]) -> &[u8] {
    if message.is_empty() || !crate::ncbi_exception::main_diag_failbit() {
        return message;
    }
    if message.starts_with(QUERY_OPTIONS_ERROR.as_bytes())
        && message
            .strip_suffix(b"\n")
            .and_then(|s| s.strip_suffix(crate::wording::USER_MANUAL_FOOTER.as_bytes()))
            .is_some_and(|s| s.ends_with(b"\n"))
    {
        b"\n\n"
    } else {
        b"\n"
    }
}

/// Set by `main_entry` (`post_errors_in_run`): the command line's run writes the text of an error that
/// gets an error archive itself, before the archive (`post_native_error`). A library caller of the
/// command line entries (`api::run_blastn` and the others) does not set it, and the text stays in
/// the error it gets back.
static POST_ERRORS_IN_RUN: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

/// The text that the run wrote (`post_native_error`), which `exit_on_native_error` does not write
/// again.
static POSTED_TEXT: std::sync::Mutex<Option<Vec<u8>>> = std::sync::Mutex::new(None);

/// Whether the run already wrote `text` (`post_native_error`).
fn posted(text: &[u8]) -> bool {
    POSTED_TEXT
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .as_deref()
        == Some(text)
}

/// Makes the command line's run write an error's text before its error archive
/// (`post_native_error`); called by `main_entry` only.
pub fn post_errors_in_run() {
    POST_ERRORS_IN_RUN.store(true, std::sync::atomic::Ordering::Relaxed);
}

/// Whether this is the command line's process (`post_errors_in_run`, which `main_entry` alone
/// calls): its streams post NCBI's texts themselves (`report::ostream::write_side_file`).
pub fn is_command_line_run() -> bool {
    POST_ERRORS_IN_RUN.load(std::sync::atomic::Ordering::Relaxed)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:327-332
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1160
/// Writes the text of `error` (a `NativeError`) to standard error now. NCBI's `CATCH_ALL` posts
/// the error inside `Run` (its handler writes the post at once and saves it for the archive)
/// before `PrintErrorArchive` writes the error archive; LOSAT's `main` writes the text after the
/// run returns (`exit_on_native_error`). A command line run that writes an error archive calls this
/// first (`report::archive::write_error_archive_after_post`) and returns that error, whose text
/// `exit_on_native_error` then does not write again. Only when `main` asked for it
/// (`post_errors_in_run`). A standard error that cannot be written loses the text, not the archive
/// after it.
pub fn post_native_error(error: &anyhow::Error) {
    if !POST_ERRORS_IN_RUN.load(std::sync::atomic::Ordering::Relaxed) {
        return;
    }
    if let Some(native) = error.downcast_ref::<NativeError>() {
        let _ = crate::report::ostream::cerr_post(catch_all_text(&native.message));
        *POSTED_TEXT
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(native.message.clone());
    }
}

#[cfg(test)]
mod manual_refusals_tests;

#[cfg(test)]
mod top_level_tests;

#[cfg(test)]
mod syntax_cases_tests;

#[cfg(all(test, unix))]
#[path = "cli_tests.rs"]
mod tests;

/// The process's arguments as NCBI's `GetCmdlineArgs` joins them for the `@PG` line of a SAM
/// report: every argument (the program name first) followed by one space, as raw bytes where
/// the system keeps them (Unix), else lossily decoded.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:877-885
pub fn command_line_bytes() -> Vec<u8> {
    let mut line = Vec::new();
    for argument in std::env::args_os() {
        line.extend_from_slice(&crate::utils::os_bytes::os_str_bytes(&argument));
        line.push(b' ');
    }
    line
}

#[cfg(test)]
#[path = "cli_v2_tests.rs"]
mod cli_v2_tests;
