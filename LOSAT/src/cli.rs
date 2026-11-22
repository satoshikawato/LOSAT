//! Public CLI v2 boundary. Clap's double-dash representation is internal only.

use crate::wording::{DATABASE_ERROR, QUERY_OPTIONS_ERROR};
use std::{ffi::OsString, fmt};

use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
use crate::algorithm::{blastn, blastp, blastx, psiblast, rpsblast, rpstblastn, tblastn, tblastx};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46,51,75
// LOSAT retains its program subcommand; search parameter names are NCBI names.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3815-3826
// `-help` and `-h` open with NCBI's two headings, `USAGE` and `DESCRIPTION`: tools that call
// BLAST+ look for them (BRIG takes a `blastn -help` that prints `USAGE` as BLAST+ found).
// The options follow in clap's layout.
const NCBI_HELP_TEMPLATE: &str = "USAGE\n  {usage}\n\nDESCRIPTION\n   {about}\n\n{all-args}";

/// The top level's `-help`: NCBI's headings, and what the entry names need.
const LOSAT_HELP_TEMPLATE: &str =
    "USAGE\n  {usage}\n\nDESCRIPTION\n   {about}\n\n{all-args}\n\n{after-help}";

pub use crate::utils::entry_names::{entry_argv, ENTRY_NAMES};

// What the `losat` binary runs; the crate's root re-exports `main_entry` for it.
pub(crate) mod main_entry;

// NCBI's buffers of a report stream (`NcbiStdio`, `NcbiBuffer`).
mod stdio_buffer;
pub use stdio_buffer::{NcbiBuffer, NcbiStdio, NCBI_BUFSIZ};

/// The words that ask for `install_shims`: `-install-shims` is how the top level's `-help` shows
/// it (`render_message` writes long options with one dash).
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

/// The `-outfmt` help that every search program shows (`blastn`, `blastp`, `blastx`, `tblastn`,
/// `tblastx`, `psiblast`, `rpsblast`): the formats the program writes (`written`, from the
/// README's "Output formats" tables), the ones that take a field list (`fields`), and the
/// default fields. `note` closes the sentence with what differs for `-in_pssm`, or is empty.
/// `rpstblastn` shows NCBI's frozen text instead.
pub(crate) fn outfmt_help(written: &str, fields: &str, note: &str) -> String {
    format!(
        "Output format. Written: {written}{note}. Formats {fields} take a field list, for example \
         -outfmt \"6 qseqid sseqid evalue\"; without one they write qaccver saccver pident length \
         mismatch gapopen qstart qend sstart send evalue bitscore"
    )
}

/// What a program that `main_entry` runs before clap reads the command line (the programs of
/// `ENTRY_NAMES` that are not a `Commands` variant: they read their words with NCBI's argument
/// reader, `ncbi_app`) says about itself in the top level's `-help`: NCBI's DESCRIPTION line of
/// the program (its first line; `convert2blastmask`'s continues in NCBI's text).
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
        _ => "",
    }
}

/// The `after_help` of the top level: the programs that `losat NAME` runs and `Commands` does not
/// list (derived from `ENTRY_NAMES` and the subcommands of `Commands`, so that a new name shows
/// up), and the shims. No NCBI counterpart.
fn top_level_after_help() -> &'static str {
    static TEXT: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    TEXT.get_or_init(|| {
        let others: Vec<&str> = ENTRY_NAMES
            .into_iter()
            .filter(|name| !<Commands as clap::Subcommand>::has_subcommand(name))
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

// No NCBI counterpart: `losat <program>` names the program, and a run builds the clap command of
// that program only (`ProgramParser`).
// LOSAT's command line, `losat <program> <words>` (`Commands`). Its command (`CommandFactory`) is
// `cli_command`'s: `Cli::command()` has every program's arguments, which the readers of the
// definitions walk (`api::describe`, LOSAT Web's `describe`, `render_message`), and the parser
// reads a command line with only its program's arguments (`ProgramParser::program_command`). A
// doc comment here would be clap's description of the top level (`long_about`).
#[derive(clap::Args, Debug)]
#[command(
    version,
    about = "Pure-Rust pairwise alignment that matches NCBI BLAST+ on LOSAT's test fixtures, with genetic-code exceptions",
    help_template = LOSAT_HELP_TEMPLATE,
    after_help = top_level_after_help()
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

// The `Parser` that `#[derive(Parser)]` writes (clap_derive 4.6.7: src/derives/parser.rs:85-107,
// src/derives/into_app.rs:21-64): `Command::new("losat")` with `Cli`'s `Args`, its program
// arguments added by `cli_command`.
impl CommandFactory for Cli {
    fn command() -> clap::Command {
        cli_command::<false>(|_| true)
    }

    fn command_for_update() -> clap::Command {
        cli_command::<true>(|_| true)
    }
}

impl Parser for Cli {}

/// The type that LOSAT's parser reads a command line into (`try_parse_from` and the readers
/// beside it): a `Parser`, read with its whole command unless it builds less for one command line.
pub trait ProgramParser: Parser {
    /// The command that reads `words` (`losat <program> <words>`): `Self::command()`.
    fn program_command(_words: &[OsString]) -> clap::Command {
        Self::command()
    }
}

/// `Cli` reads a command line with the arguments of its program only, the program that the reader
/// reads (`read_program`), so that each option of another program costs a run nothing.
impl ProgramParser for Cli {
    fn program_command(words: &[OsString]) -> clap::Command {
        let program = read_program(words, <Commands as Subcommand>::has_subcommand);
        cli_command::<false>(|name| Some(name) == program)
    }
}

/// `Cli`'s command with the arguments of `program` only (`cli_command`).
pub(crate) fn command_of(program: &str) -> clap::Command {
    cli_command::<false>(|name| name == program)
}

/// The words that LOSAT's reader takes at the top level (`parse_command`), before the program's
/// name: its help and version flags. The top level has no option of its own
/// (`Cli::command().get_arguments()` is empty before clap builds it).
const TOP_LEVEL_FLAGS: [&str; 4] = ["-help", "--help", "-h", "--version"];

/// The program whose arguments LOSAT's reader reads in `words` (`losat <words>`,
/// `parse_command`): the first word after the top level's flags (`TOP_LEVEL_FLAGS`), when
/// `is_program` takes it. The reader rejects any other word there, before it reads a program, so
/// such a command line reads no program's arguments. The command that `Cli` is read with has the
/// arguments of this program only (`program_command`), and the reader checks that the program it
/// enters is this one, with its arguments.
fn read_program(words: &[OsString], is_program: impl Fn(&str) -> bool) -> Option<&str> {
    let word = words
        .iter()
        .skip(1)
        .find(|word| !TOP_LEVEL_FLAGS.iter().any(|flag| word.as_os_str() == *flag))?;
    word.to_str().filter(|word| is_program(word))
}

/// `Cli`'s command (`augment_args`, or `augment_args_for_update` with `UPDATE`) where the programs
/// whose names `selected` takes have their arguments (`program_args`). The others have their
/// names, descriptions and help template, which the top level's help lists, and no arguments to
/// LOSAT's reader, which reads a program's arguments before clap builds anything
/// (`parse_command`): it reads only a selected program. clap adds the others' arguments where it
/// builds them (`#[command(defer = true)]` on `Commands`). `UPDATE` is a constant so that a
/// program that never calls `command_for_update` carries none of its code.
fn cli_command<const UPDATE: bool>(selected: impl Fn(&str) -> bool) -> clap::Command {
    let command = clap::Command::new("losat");
    let (mut command, program_args) = if UPDATE {
        (
            <Cli as clap::Args>::augment_args_for_update(command),
            PROGRAM_ARGS_FOR_UPDATE,
        )
    } else {
        (<Cli as clap::Args>::augment_args(command), PROGRAM_ARGS)
    };
    debug_assert_eq!(command.get_subcommands().count(), program_args.len());
    for (program, &args) in command.get_subcommands_mut().zip(program_args) {
        let deferred = std::mem::take(program);
        *program = if selected(deferred.get_name()) {
            args(deferred.defer(no_args))
        } else {
            deferred.defer(args)
        };
    }
    command
}

/// What clap builds of a program whose arguments `cli_command` has added: nothing more.
fn no_args(command: clap::Command) -> clap::Command {
    command
}

/// What `Commands` defers for a program (`#[command(defer = true)]`): the arguments that
/// `cli_command` adds for `CommandFactory::command` (`PROGRAM_ARGS`) or for `command_for_update`
/// (`PROGRAM_ARGS_FOR_UPDATE`).
type ProgramArgs = fn(clap::Command) -> clap::Command;

// clap_derive 4.6.7: src/derives/subcommand.rs:288-358 (a variant's command: its description and
// help template, `final_from_attrs`, after its `Args`, or before them when deferred)
/// The arguments of the program whose command is `command` (`A::augment_args`, or
/// `augment_args_for_update` with `UPDATE`), as clap's derive adds them where it does not defer
/// them: the description stays the variant's. A deferred `Args` comes after the variant's
/// description, and the doc comment of an `Args` (`RpsbprocArgs`, `MakedbArgs` and the groups that
/// the programs flatten) is a description to clap. Every variant has a description and no `Args`
/// sets a help template (`cli_tests`), so the rest of the variant's command stays as the derive
/// writes it.
fn program_args<A: clap::Args, const UPDATE: bool>(command: clap::Command) -> clap::Command {
    use clap::builder::Resettable;
    let about = Resettable::from(command.get_about().cloned());
    let long_about = Resettable::from(command.get_long_about().cloned());
    let command = if UPDATE {
        A::augment_args_for_update(command)
    } else {
        A::augment_args(command)
    };
    command.about(about).long_about(long_about)
}

/// `Commands`, `PROGRAM_ARGS` and `PROGRAM_ARGS_FOR_UPDATE` from one list of the programs (a
/// program is added in one place).
macro_rules! programs {
    (
        $(#[$($attr:tt)*])*
        $vis:vis enum $name:ident {
            $( $(#[$($variant_attr:tt)*])* $variant:ident($($args:tt)+), )+
        }
    ) => {
        $(#[$($attr)*])*
        $vis enum $name {
            $( $(#[$($variant_attr)*])* $variant($($args)+), )+
        }

        /// Each program's arguments (`program_args` of its variant's type), in the order of
        /// `Commands`' variants: the order of `Cli`'s subcommands.
        const PROGRAM_ARGS: &[ProgramArgs] = &[ $( program_args::<$($args)+, false>, )+ ];

        /// `PROGRAM_ARGS` for `command_for_update` (`augment_args_for_update`).
        const PROGRAM_ARGS_FOR_UPDATE: &[ProgramArgs] = &[ $( program_args::<$($args)+, true>, )+ ];
    };
}

programs! {
    #[derive(Subcommand, Debug, Clone)]
    #[command(defer = true)]
    pub enum Commands {
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
        /// Translated nucleotide query vs protein subject (blastx)
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Blastx(blastx::BlastxArgs),
        /// Pairwise nucleotide alignment (megablast [default], blastn, dc-megablast, blastn-short)
        #[allow(rustdoc::broken_intra_doc_links, reason = "clap prints it as help")]
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Blastn(blastn::BlastnArgs),
        /// Pairwise protein alignment (blastp)
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Blastp(blastp::BlastpArgs),
        /// Pairwise 6-frame translated nucleotide alignment (tblastx)
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Tblastx(tblastx::TblastxArgs),
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:47-52:
        /// Protein query vs translated nucleotide subject (tblastn)
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Tblastn(tblastn::TblastnArgs),
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:47-51
        /// Position-Specific Iterated BLAST (psiblast)
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Psiblast(psiblast::PsiblastArgs),
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:77-79
        /// Reverse Position Specific BLAST (rpsblast)
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Rpsblast(rpsblast::RpsblastArgs),
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:48-51
        // `-h`, `-help` and `-help-full` write NCBI's frozen texts (`RpstblastnArgs::ncbi_help`; the
        // reader stops at them, `try_parse_from_command`).
        /// Translated Reverse Position Specific BLAST (rpstblastn)
        #[command(help_template = NCBI_HELP_TEMPLATE)]
        Rpstblastn(rpstblastn::RpstblastnArgs),
        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:217-219
        /// makeblastdb: create a BLAST database from a FASTA file (format version 5 by default, 4 with -blastdb_version 4)
        #[command(name = "makeblastdb", help_template = NCBI_HELP_TEMPLATE)]
        Makedb(crate::writedb::makeblastdb::MakedbArgs),
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:64-65
        /// Process output from offline rpsblast to annotate sequence with conserved domain information (rpsbproc; its domain hits, sites and motifs of protein and nucleotide queries)
        #[command(name = "rpsbproc", help_template = NCBI_HELP_TEMPLATE)]
        Rpsbproc(crate::rpsbproc::RpsbprocArgs),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-170,203-207,332-349
// Translate registered canonical keys only. Consume each value once so a path
// or a quoted filtering specification can never be interpreted as another key.
/// Parses `argv` (`losat <program> <words>`) for a caller that runs no NCBI application start
/// (LOSAT Web's ABI v1, and the engine's own parses of an argv): the toolkit's words that the
/// start reads (`-logfile`, `-conffile`, `-dryrun`) are refused (`toolkit_file_word_error`;
/// makeblastdb keeps its own `-logfile` argument), and so are the versions (`version_word_error`)
/// and `-xmlhelp` (`xml_help_refusal`).
pub fn try_parse_from<T, I, S>(argv: I) -> Result<T, clap::Error>
where
    T: ProgramParser,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let words: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    try_parse_from_command(T::program_command(&words), words)
}

/// `try_parse_from` for a caller whose calls run an NCBI application start of their own and return
/// reports: the library API (`api::start`). `-logfile`, `-conffile` and `-dryrun` are taken; the
/// versions and `-xmlhelp`, whose texts the library has no place for, are refused
/// (`version_word_error`, `xml_help_refusal`).
pub(crate) fn try_parse_with_own_start<T, I, S>(argv: I) -> Result<T, clap::Error>
where
    T: ProgramParser,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let words: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    parse_command(T::program_command(&words), words, Reader::OwnStart)
}

/// `try_parse_from` for LOSAT Web's ABI v2, whose runs make their own NCBI application start
/// (`web_adapter::start_run`) and give the host the text that the command line prints in place of
/// a run (stream 11; line WQA's D14, 2026-11-20): `-logfile`, `-conffile` and `-dryrun` are taken,
/// and `-xmlhelp` stops the reader with the request whose XML the run prints (`printed_text`), as
/// on the command line. The versions are printed before the arguments are read
/// (`ncbi_version_request`, `api::start::CallerStart::version`), so a version word that reaches
/// this reader is refused (`version_word_error`).
pub(crate) fn try_parse_with_printing_start<T, I, S>(argv: I) -> Result<T, clap::Error>
where
    T: ProgramParser,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let words: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    parse_command(T::program_command(&words), words, Reader::PrintingStart)
}

/// `try_parse_from` after the command line's application start `start` of the program that
/// `argv` names (`cli::main_entry`): the start has read `-logfile` and `-conffile`, the entry
/// ends a `-dryrun` after the arguments, and `-xmlhelp` stops the reader with the request whose
/// XML the entry prints (`xml_help_program`).
pub(crate) fn try_parse_from_start<T, I, S>(
    argv: I,
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> Result<T, clap::Error>
where
    T: ProgramParser,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let words: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    parse_command(
        T::program_command(&words),
        words,
        Reader::CommandLine(start),
    )
}

/// `try_parse_from` with `command`, one of `T`'s commands whose arguments a caller changed
/// (Web ABI v1's BLASTX keeps its own `-outfmt` values, `blastx/web.rs`), read as clap's
/// `Parser::try_parse_from` reads `T::command()`. The program that `argv` names has its arguments
/// in `command` (`program_command`).
pub(crate) fn try_parse_from_command<T, I, S>(
    command: clap::Command,
    argv: I,
) -> Result<T, clap::Error>
where
    T: clap::FromArgMatches,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    parse_command(command, argv, Reader::NoStart)
}

/// Who reads an argv (`parse_command`), which decides what NCBI's toolkit words do there.
#[derive(Clone, Copy)]
enum Reader<'a> {
    /// A caller that runs no NCBI application start (LOSAT Web's ABI v1): the toolkit's words,
    /// the versions and `-xmlhelp` are refused.
    NoStart,
    /// The library API, whose calls run a start of their own (`api::start`): `-logfile`,
    /// `-conffile` and `-dryrun` are taken; the versions and `-xmlhelp` are refused.
    OwnStart,
    /// LOSAT Web's ABI v2, whose runs make a start of their own (`web_adapter::start_run`) and
    /// print the command line's texts: the toolkit's words are taken, and `-xmlhelp` is the
    /// request whose XML the run prints, as on the command line.
    PrintingStart,
    /// The command line's entry with its application start: `-xmlhelp` is the request whose XML
    /// the entry prints.
    CommandLine(&'a crate::blastinput::ncbi_environment::ApplicationStart),
}

/// The reader of `try_parse_from`, `try_parse_with_own_start`, `try_parse_with_printing_start`,
/// `try_parse_from_start` and `try_parse_from_command`.
fn parse_command<T, I, S>(
    command: clap::Command,
    argv: I,
    reader: Reader<'_>,
) -> Result<T, clap::Error>
where
    T: clap::FromArgMatches,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let mut tokens: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    if let Reader::CommandLine(start) = reader {
        debug_assert_eq!(
            tokens.get(1).and_then(|token| token.to_str()),
            Some(start.program()),
            "the start of another program"
        );
    }
    // The program that the reader reads (`read_program`), the only one whose arguments it reads,
    // with the toolkit's file arguments.
    let read = read_program(&tokens, |word| command.find_subcommand(word).is_some())
        .and_then(|word| command.find_subcommand(word))
        .map(|program| program.get_name().to_string());
    let command = with_toolkit_file_args(command, |name| Some(name) == read.as_deref());
    // NCBI's application reads its own toolkit words before the program's arguments; the
    // argument reader gets the words that are left (`ncbi_environment::preparse`).
    if let Some(program) = tokens
        .get(1)
        .and_then(|token| token.to_str())
        .and_then(|text| command.find_subcommand(text))
        .filter(|program| TOOLKIT_PROGRAMS.contains(&program.get_name()))
    {
        if let Some(word) = ncbi_version_word(&tokens[2..]) {
            return Err(version_word_error(program.get_name(), word));
        }
        if matches!(reader, Reader::NoStart) {
            if let Some(word) = toolkit_file_word(&tokens[2..]) {
                return Err(toolkit_file_word_error(program.get_name(), word));
            }
        }
        use crate::blastinput::ncbi_environment::{preparse, Preparsed};
        if let Preparsed::Words { words, .. } = preparse(&tokens[2..]) {
            tokens.truncate(2);
            tokens.extend(words);
        }
    } else if tokens.get(1).is_some_and(|token| token == "makeblastdb") {
        // makeblastdb's application reads the same words (its own refusals of the others come
        // with its arguments): a last `-logfile` or `-conffile` without a value is dropped,
        // `-dryrun` goes to the entry (`main_entry`) or the library's call
        // (`api::run_makeblastdb`), and a caller without a start refuses `-conffile` and
        // `-dryrun` (NCBI BLAST+ 2.17.0, probes `makeblastdb.log_noval`,
        // `makeblastdb.conf_noval`: rc 0).
        use crate::blastinput::ncbi_environment::{preparse, Preparsed};
        if let Some(word) = ncbi_version_word(&tokens[2..]) {
            return Err(version_word_error("makeblastdb", word));
        }
        if let Preparsed::Words {
            words,
            conffile,
            dryrun,
        } = preparse(&tokens[2..])
        {
            if matches!(reader, Reader::NoStart) && (conffile.is_some() || dryrun) {
                // makeblastdb takes `-logfile` (its own argument): the word named is the one
                // refused.
                let name = if conffile.is_some() {
                    "-conffile"
                } else {
                    "-dryrun"
                };
                let word = toolkit_word_as_written(&tokens[2..], name).unwrap_or(name);
                return Err(toolkit_file_word_error("makeblastdb", word));
            }
            tokens.truncate(2);
            tokens.extend(words);
        }
    }
    let mut input = tokens.into_iter();
    let mut translated = vec![input.next().unwrap_or_else(|| "losat".into())];
    let mut scope = &command;
    // rpstblastn's help flag: the words up to its program name and the value of its
    // `RpstblastnArgs::ncbi_help`.
    let mut ncbi_help: Option<(usize, &str)> = None;
    // The first option of the NCBI search programs whose value is not UTF-8, rejected once the
    // command line parses (`non_utf8_value_error`).
    let mut non_utf8_value: Option<(&clap::Command, String)> = None;
    // The program whose `-xmlhelp` stopped the reader.
    let mut xml_help: Option<String> = None;
    while let Some(token) = input.next() {
        // In the NCBI search programs a word that is not UTF-8 keeps its key readable
        // (`-name=value` with a value that is not UTF-8 is read below).
        let lossy = token.to_string_lossy();
        let text: &str = if is_ncbi_search_program(scope) {
            &lossy
        } else {
            token.to_str().unwrap_or("")
        };
        if let Some(subcommand) = scope.find_subcommand(text) {
            // The program that `read_program` names, whose arguments the command has
            // (`ProgramParser::program_command`): the reader looks its options up itself.
            debug_assert_eq!(
                Some(subcommand.get_name()),
                read.as_deref(),
                "{translated:?}"
            );
            debug_assert!(
                subcommand.get_arguments().next().is_some(),
                "{} read without its arguments",
                subcommand.get_name()
            );
            scope = subcommand;
            translated.push(token);
            if scope.get_name() == "rpstblastn" {
                ncbi_help = Some((translated.len(), ""));
            }
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3020-3023, 2840-2852
        // NCBI's reader stops at its help flags, the words before them read first: a word after
        // one is never read (NCBI BLAST+ 2.17.0: `-help -bogus`, `-help ... -out` with no value
        // and `-h -xmlhelp` print the help, exit status 0; line RG item E2). Every program of
        // this parser follows NCBI's reader; rpstblastn writes NCBI's texts
        // (`RpstblastnArgs::ncbi_help`). An option's value is not a flag (`-out -help` names a
        // file). `--help` is LOSAT's spelling of `-help` (clap's) and stops the reader too.
        // `-help-full` (`s_AutoHelpShowAll`) is `-help` with the arguments defined `fHidden`
        // (`ShowAllArguments`). Of NCBI's programs only makeblastdb has one, `-randomize`, which
        // LOSAT's help (its own texts, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 1) does not
        // list, so the text is `-help`'s, as NCBI's other programs write the bytes of `-help` (in
        // every place: line RG, `E/SURVEY_E2.md`).
        let stops_at_help = !std::ptr::eq(scope, &command);
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:685-690 (cited at
        // `ncbi_xmlhelp`)
        // `-xmlhelp` (`s_AutoHelpXml`) stops the reader as `-help` does, the words before it
        // read: the command line's entry prints the program's arguments in NCBI's XML
        // (`xml_help_program`, `ncbi_xmlhelp::print`; line RG item E2x), and so does a run of
        // LOSAT Web's ABI v2 (`printed_text`; line WQA's D14, 2026-11-20). The library API and the
        // Web engine's ABI v1, which write no help of NCBI's, refuse it there, after the words
        // before it (`xml_help_refusal`).
        if stops_at_help && text == "-xmlhelp" {
            translated.push("--help".into());
            xml_help = Some(scope.get_name().to_string());
            break;
        }
        if matches!(text, "-help" | "--help") || (stops_at_help && text == "-help-full") {
            translated.push("--help".into());
            if stops_at_help {
                ncbi_help = ncbi_help.map(|(start, _)| (start, "full"));
                break;
            }
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:76-77
        // `-h` is the short usage of every NCBI program (rc 0); clap's `-h` is its short help.
        if text == "-h" {
            translated.push("-h".into());
            if stops_at_help {
                ncbi_help = ncbi_help.map(|(start, _)| (start, "short"));
                break;
            }
            continue;
        }
        if text == "--version" && std::ptr::eq(scope, &command) {
            translated.push(token);
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2866-2872
        // `--` starts NCBI's positional arguments, which the BLAST programs do not have:
        // a last `--` changes nothing, and a word after it is an extra positional
        // argument (NCBI's USAGE error; LOSAT's parser error, exception 1).
        // makeblastdb reads its arguments with the same parser.
        if text == "--"
            && input.len() == 0
            && matches!(
                scope.get_name(),
                "blastn"
                    | "blastp"
                    | "tblastn"
                    | "tblastx"
                    | "psiblast"
                    | "rpsblast"
                    | "rpstblastn"
                    | "makeblastdb"
                    | "rpsbproc"
            )
        {
            continue;
        }
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/argwrapper.cpp:42-43
        // rpsbproc's aliases (`-infile` of `i`, ...) are its `--infile`, ... (`rpsbproc::args`).
        let alias = (scope.get_name() == "rpsbproc")
            .then(|| crate::rpsbproc::alias_word(text))
            .flatten();
        let text: &str = alias.as_deref().unwrap_or(text);
        let key = text.strip_prefix('-').filter(|key| !key.starts_with('-'));
        let (name, inline) = key
            .unwrap_or("")
            .split_once('=')
            .map_or((key.unwrap_or(""), None), |(k, v)| (k, Some(v)));
        // rpstblastn's `RpstblastnArgs::ncbi_help` is the reader's own, named by no word; the
        // other hidden arguments (the toolkit's `-logfile`, `-conffile`) are the user's.
        let Some(arg) = scope
            .get_arguments()
            .find(|arg| arg.get_long() == Some(name))
            .filter(|_| name != rpstblastn::args::NCBI_HELP_ARG)
        else {
            return Err(unknown_option_error(scope, name, text));
        };
        if arg.get_action().takes_values() {
            let value = match inline {
                Some(value) if token.to_str().is_some() => OsString::from(value),
                Some(_) => inline_value(&token),
                None => input.next().ok_or_else(|| {
                    clap::Error::raw(
                        ErrorKind::InvalidValue,
                        format!("-{name} requires one value"),
                    )
                })?,
            };
            // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:266-270
            // NCBI keeps the bytes of an option's value. LOSAT's options are UTF-8 strings,
            // so a value that is not UTF-8 is rejected explicitly, after the parser: NCBI's
            // USAGE errors and `-help` come first, as for any value. The parser reads the
            // value with U+FFFD in place of the bytes. The file names keep their bytes and
            // have their own check where the files are opened (`check_utf8_file_name`).
            let value = if value.to_str().is_none() && refuses_non_utf8_value(scope, name) {
                non_utf8_value.get_or_insert_with(|| (scope, name.to_string()));
                OsString::from(value.to_string_lossy().into_owned())
            } else {
                value
            };
            let mut internal = OsString::from(format!("--{name}="));
            internal.push(value);
            translated.push(internal);
        } else {
            if inline.is_some() {
                return Err(clap::Error::raw(
                    ErrorKind::TooManyValues,
                    format!("-{name} is a flag and does not take a value"),
                ));
            }
            translated.push(format!("--{name}").into());
        }
    }
    let mut matches = match command.clone().try_get_matches_from(&translated) {
        // The words before `-xmlhelp` read, as before `-help`: the entry prints the XML.
        Err(error) if error.kind() == ErrorKind::DisplayHelp && xml_help.is_some() => {
            let program = xml_help.unwrap_or_default();
            return Err(match reader {
                Reader::CommandLine(_) | Reader::PrintingStart => xml_help_request(program),
                Reader::NoStart | Reader::OwnStart => clap::Error::raw(
                    ErrorKind::InvalidValue,
                    xml_help_refusal(&refusal_program_name(&program)),
                ),
            });
        }
        // The parser stopped at rpstblastn's help flag, the words before it read (NCBI's
        // order): the program writes NCBI's text (`RpstblastnArgs::ncbi_help_text`), and the
        // words before the flag, read, are not used (NCBI checks their dependencies after).
        Err(error) if error.kind() == ErrorKind::DisplayHelp => match ncbi_help {
            Some((start, value)) if !value.is_empty() => {
                translated.truncate(start);
                translated.push(format!("--ncbi_help={value}").into());
                let mut matches = command.clone().try_get_matches_from(translated)?;
                return T::from_arg_matches_mut(&mut matches)
                    .map_err(|error| error.format(&mut command.clone()));
            }
            _ => return Err(error),
        },
        matches => matches?,
    };
    let parsed = T::from_arg_matches_mut(&mut matches)
        .map_err(|error| error.format(&mut command.clone()))?;
    match non_utf8_value {
        Some((scope, name)) => Err(non_utf8_value_error(scope, &name)),
        None => Ok(parsed),
    }
}

/// The value of a word `-name=value` that is not UTF-8: its bytes after the first `=`
/// (`try_parse_from`, which found the name before the first `=` of the word's lossy text; a
/// replaced byte is never `=`, so the first `=` of the bytes is that one, whatever the name's
/// letters are).
fn inline_value(token: &std::ffi::OsStr) -> OsString {
    let bytes = token.as_encoded_bytes();
    let start = bytes
        .iter()
        .position(|&b| b == b'=')
        .map_or(bytes.len(), |at| at + 1);
    // SAFETY: `start` is just after an ASCII `=` or at the end, split points that
    // `OsStr::from_encoded_bytes_unchecked` allows (after a valid non-empty UTF-8 substring).
    unsafe { std::ffi::OsStr::from_encoded_bytes_unchecked(&bytes[start..]) }.to_os_string()
}

/// Whether `scope` is one of the programs whose command lines follow NCBI's.
fn is_ncbi_search_program(scope: &clap::Command) -> bool {
    matches!(
        scope.get_name(),
        "blastn"
            | "blastp"
            | "blastx"
            | "tblastn"
            | "tblastx"
            | "psiblast"
            | "rpsblast"
            | "rpstblastn"
    )
}

/// Whether an option's value that is not UTF-8 is rejected after the parser
/// ([`non_utf8_value_error`]): each option of makeblastdb, whose values are all UTF-8 strings,
/// and each option of the NCBI search programs but the file names that keep their bytes
/// (`check_utf8_file_name`).
fn refuses_non_utf8_value(scope: &clap::Command, name: &str) -> bool {
    scope.get_name() == "makeblastdb"
        || (is_ncbi_search_program(scope) && !matches!(name, "query" | "subject" | "out"))
}

/// The explicit rejection of an option's value that is not UTF-8 (`try_parse_from`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:266-270
fn non_utf8_value_error(scope: &clap::Command, name: &str) -> clap::Error {
    clap::Error::raw(
        ErrorKind::InvalidUtf8,
        format!(
            "the value of -{name} is not UTF-8; NCBI BLAST+ reads the bytes of an option's value, which is not supported by LOSAT's {}",
            scope.get_name().to_uppercase()
        ),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:107,143
// Render the same canonical names in help, usage and parser diagnostics.
/// `text` with each `--{name}` that is a whole word (not followed by a letter, a digit, `-` or `_`)
/// written `-{name}` (`render_rpsbproc_message`).
fn replace_whole_option(text: &str, name: &str) -> String {
    let pattern = format!("--{name}");
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find(&pattern) {
        let after = &rest[at + pattern.len()..];
        let whole = !after
            .chars()
            .next()
            .is_some_and(|c| c.is_alphanumeric() || c == '-' || c == '_');
        out.push_str(&rest[..at]);
        out.push_str(if whole { &pattern[1..] } else { &pattern });
        rest = after;
    }
    out.push_str(rest);
    out
}

pub fn render_message(error: &clap::Error) -> String {
    let mut message = error
        .to_string()
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version");
    // Every program's option names, whatever the program of the message: a word of the command
    // line that names another program's option is written with one dash too.
    let command = with_toolkit_file_args(Cli::command(), |_| true);
    for scope in std::iter::once(&command).chain(command.get_subcommands()) {
        for arg in scope.get_arguments() {
            // rpsbproc's one-letter arguments are left to `render_rpsbproc_message`: replacing
            // `--i` here would change other words (`--infile`) in every program's messages.
            if let Some(name) = arg
                .get_long()
                .filter(|name| !matches!(*name, "help" | "version") && name.len() > 1)
            {
                message = message.replace(&format!("--{name}"), &format!("-{name}"));
            }
        }
    }
    // clap's hints at the foot of the help name the other form: `-h` prints the summary and
    // `-help` the full text, as in NCBI, so the hint of the full text points back to `-h` and
    // that of the summary to `-help`.
    message
        .replace("'-h'", "'-help'")
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version")
        .replace("see a summary with '-help'", "see a summary with '-h'")
        .replace("see more with '--help'", "see more with '-help'")
}

/// `render_message` for `losat rpsbproc`, whose arguments are NCBI's one-letter keys (`-i`)
/// with long forms of two dashes (`--infile`, rpsbproc's `AddAlias`): only the one-letter
/// arguments are written with one dash, and the long forms keep their two.
pub fn render_rpsbproc_message(error: &clap::Error) -> String {
    let text = error.to_string();
    // The word that `unknown_option_error` quotes is the one on the command line (`--i` stays
    // `--i`); the rest is rendered.
    const UNKNOWN: &str = "unknown option or argument '";
    if let Some(start) = text.find(UNKNOWN).map(|at| at + UNKNOWN.len()) {
        if let Some(len) = text[start..].find('\'') {
            return format!(
                "{}{}{}",
                render_rpsbproc_text(&text[..start]),
                &text[start..start + len],
                render_rpsbproc_text(&text[start + len..])
            );
        }
    }
    render_rpsbproc_text(&text)
}

/// The text of a parser error as the command line writes it, which LOSAT Web's ABI v2 returns as
/// it is (`web/adapter/tests/v_abi.js` compares the two): `render_message`, or
/// `render_rpsbproc_message` for `losat rpsbproc`, ending with a line end. clap's own errors have
/// one; LOSAT's own refusals and syntax errors (`clap::Error::raw`: an option LOSAT refuses by
/// name, an unknown option, a missing value) are made without one, and the shell prompt would
/// follow on the message's line. ABI v1 (frozen) and the library return `render_message`.
pub fn command_line_message(error: &clap::Error, rpsbproc: bool) -> String {
    let mut message = if rpsbproc {
        render_rpsbproc_message(error)
    } else {
        render_message(error)
    };
    if error.use_stderr() && !message.ends_with('\n') {
        message.push('\n');
    }
    message
}

/// `render_rpsbproc_message` of a part of a message.
fn render_rpsbproc_text(text: &str) -> String {
    let mut message = text
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version");
    let command = command_of("rpsbproc");
    if let Some(scope) = command.find_subcommand("rpsbproc") {
        for name in scope.get_arguments().filter_map(clap::Arg::get_long) {
            if name.len() == 1 {
                message = replace_whole_option(&message, name);
            }
        }
    }
    message.replace("'-h'", "'-help'")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:63-125
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:63-70
// The options of NCBI blastn 2.17.0+ (-help) that LOSAT's BLASTN does not implement
// (AGENTS.md rule 2: explicit unsupported errors). `-db` and `-dbsize` are implemented
// (`blastn/args.rs`), as are the X-drops, `-window_size`, `-off_diagonal_range`, `-no_greedy`
// (DI-6), `-min_raw_gapped_score` (DI-6b), the database masks (`blastinput::db_mask`, series
// NP M4a-2), `-filtering_db` (`algorithm::blastn::repeat_filter`, line RG item D4) and
// `-best_hit_overhang` / `-best_hit_score_edge` (series NP O6a).
const UNPORTED_BLASTN_ARGS: &[&str] = &[
    "entrez_query",
    "gilist",
    "html",
    "import_search_strategy",
    "negative_gilist",
    "remote",
];

fn is_unported_blastn_arg(name: &str) -> bool {
    UNPORTED_BLASTN_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:44-50
// The options of NCBI blastp 2.17.0+ (-help) that LOSAT's BLASTP does not implement
// (AGENTS.md rule 2: explicit unsupported errors). `-best_hit_overhang` / `-best_hit_score_edge`
// are implemented (series NP O6a-2, psiblast's too).
const UNPORTED_BLASTP_ARGS: &[&str] = &[
    "entrez_query",
    "gilist",
    "html",
    "import_search_strategy",
    "ipglist",
    "negative_gilist",
    "negative_ipglist",
    "remote",
];

fn is_unported_blastp_arg(name: &str) -> bool {
    UNPORTED_BLASTP_ARGS.contains(&name)
}

// The rejection of an option that a program does not have: NCBI's options that LOSAT
// has not ported, the NCBI C++ Toolkit options, and the rest (unknown to NCBI too).
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1156-1164
/// psiblast's NCBI arguments that LOSAT does not run: PHI-BLAST (PS8; the ASCII and ASN.1 PSSM
/// files are PS2's and PS3's, `algorithm::psiblast::pssm`, `-in_pssm` PS4's,
/// `algorithm::psiblast::in_pssm`, the MSA inputs PS6's, `algorithm::psiblast::in_msa`), and the arguments that blastp's engine does
/// not run either (`is_unported_blastp_arg`; psiblast has no database masks, NCBI's unknown
/// arguments of its parser, and no `-mt_mode`).
/// The search strategy file is written, of a FASTA query and of `-in_pssm` (items F1b-1 and F1b-2
/// of line RG).
const UNPORTED_PSIBLAST_OWN_ARGS: &[&str] = &["phi_pattern"];

fn is_unported_psiblast_arg(name: &str) -> bool {
    UNPORTED_PSIBLAST_OWN_ARGS.contains(&name) || is_unported_blastp_arg(name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:85-127 (cited at
// `algorithm::rpsblast::args`)
/// rpsblast's NCBI arguments that LOSAT does not run: `-html` (PS15a) and the remote options and
/// `-import_search_strategy` of every program.
const UNPORTED_RPSBLAST_ARGS: &[&str] =
    &["entrez_query", "html", "import_search_strategy", "remote"];

fn is_unported_rpsblast_arg(name: &str) -> bool {
    UNPORTED_RPSBLAST_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:44-110 (cited at
// `algorithm::rpstblastn::args`)
/// rpstblastn's NCBI arguments that LOSAT does not run: `-html` (PS15a) and the remote options
/// and `-import_search_strategy` of every program. The composition-based statistics of a later
/// item are refused once per run (`rpstblastn::search::TranslatedRpsSearch::check_options`).
const UNPORTED_RPSTBLASTN_ARGS: &[&str] =
    &["entrez_query", "html", "import_search_strategy", "remote"];

fn is_unported_rpstblastn_arg(name: &str) -> bool {
    UNPORTED_RPSTBLASTN_ARGS.contains(&name)
}

fn unknown_option_error(scope: &clap::Command, name: &str, text: &str) -> clap::Error {
    named_refusal(scope.get_name(), name).unwrap_or_else(|| {
        clap::Error::raw(
            ErrorKind::UnknownArgument,
            format!("unknown option or argument '{text}'; use -help to list the options"),
        )
    })
}

/// Whether the reader of the program `scope` (`with_toolkit_file_args` of its command)
/// refuses NCBI's argument `-name` by name (`named_refusal`, for a name the program does not
/// define): the arguments that its `-xmlhelp` leaves out of NCBI's XML (`xml_help_text`).
fn refuses_option(scope: &clap::Command, name: &str) -> bool {
    let program = scope.get_name();
    (!scope
        .get_arguments()
        .any(|arg| arg.get_long() == Some(name))
        && named_refusal(program, name).is_some())
        || RUN_REFUSALS
            .iter()
            .any(|refusal| refusal.0 == program && refusal.1 == name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-918 (cited at
// blastp's check, `validate_requested_blastp_support`)
/// NCBI's arguments that a program's parser reads, so that NCBI's own errors about them and the
/// other arguments come first (tblastn's `-ungapped` with composition-based statistics:
/// `Composition-adjusted searched are not supported with an ungapped search`), and whose run LOSAT
/// refuses whatever their value, by name, where NCBI starts the search: the program, the name and
/// what the refusal says the argument is (`run_refusal`). `-xmlhelp` leaves them out of NCBI's XML
/// as it leaves the parser's refusals (`refuses_option`).
pub(crate) const RUN_REFUSALS: &[(&str, &str, &str)] = &[
    ("psiblast", "use_sw_tback", "the Smith-Waterman traceback"),
    ("tblastn", "ungapped", "an ungapped search"),
    (
        "rpstblastn",
        "ungapped",
        "NCBI BLAST+ 2.17.0's rpstblastn crashes after an ungapped search",
    ),
];

/// LOSAT's refusal of the run of `-name` in `program` (its name in LOSAT's refusals, `BLASTP`):
/// the entry of `RUN_REFUSALS`.
pub(crate) fn run_refusal(program: &str, name: &str) -> String {
    let entry = RUN_REFUSALS
        .iter()
        .find(|refusal| refusal.0.eq_ignore_ascii_case(program) && refusal.1 == name);
    debug_assert!(entry.is_some(), "{program} -{name} is not in RUN_REFUSALS");
    match entry {
        Some((_, _, what)) => {
            format!("-{name} ({what}) is not supported by LOSAT's {program}")
        }
        None => format!("-{name} is not supported by LOSAT's {program}"),
    }
}

/// LOSAT's refusal of NCBI's argument `-name`, which `program` does not define, by name; `None`
/// for a name that LOSAT does not refuse by name (an unknown option, as NCBI's).
fn named_refusal(program: &str, name: &str) -> Option<clap::Error> {
    let refusal = |shown: &str| {
        Some(clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's {shown}"),
        ))
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:64-129:
    // m_Args.push_back(arg) registers shared, formatting, DB, and PSI groups.
    // The named NCBI options below have no Stage B Rust behavior.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-66
    // Product scope explicitly refuses unported NCBI capabilities.
    if program == "blastx" && is_unported_blastx_arg(name) {
        return refusal("BLASTX");
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:246-250
    // makeblastdb's options that LOSAT's makeblastdb does not serve are named in
    // its refusal (AGENTS.md rule 2: explicit unsupported errors).
    if program == "makeblastdb" && crate::writedb::makeblastdb::UNSUPPORTED_OPTIONS.contains(&name)
    {
        // `-gi_mask` and `-gi_mask_name` make GI-indexed mask files (makeblastdb.cpp:279-288),
        // which LOSAT's writer does not write yet (line NP, item M4c); `-mask_data` itself works.
        if matches!(name, "gi_mask" | "gi_mask_name") {
            return refusal(
                "makeblastdb (it writes GI-indexed mask files, which LOSAT does not write yet; -mask_data, -mask_id and -mask_desc work)",
            );
        }
        return refusal("makeblastdb");
    }
    if program == "blastn" && is_unported_blastn_arg(name) {
        return refusal("BLASTN");
    }
    if program == "blastp" && is_unported_blastp_arg(name) {
        return refusal("BLASTP");
    }
    if program == "psiblast" && is_unported_psiblast_arg(name) {
        return refusal("PSIBLAST");
    }
    if program == "rpsblast" && is_unported_rpsblast_arg(name) {
        return refusal("RPSBLAST");
    }
    if program == "rpstblastn" && is_unported_rpstblastn_arg(name) {
        return refusal("RPSTBLASTN");
    }
    if program == "tblastx" && is_unported_tblastx_arg(name) {
        return refusal("TBLASTX");
    }
    if program == "tblastn" && is_unported_tblastn_arg(name) {
        return refusal("TBLASTN");
    }
    None
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
/// argument is parsed. The top level `losat --version` is clap's and is not read here.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:972-977
pub fn ncbi_version_request<I, S>(argv: I) -> Option<String>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let tokens: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    // A program's name (`Commands`), with no command built.
    let name = tokens
        .get(1)?
        .to_str()
        .filter(|word| <Commands as Subcommand>::has_subcommand(word))?;
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

/// The NCBI applications whose words before the arguments the parser reads (`parse_command`,
/// `ncbi_environment::preparse`).
const TOOLKIT_PROGRAMS: &[&str] = &[
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "rpsblast",
    "rpstblastn",
    "rpsbproc",
];

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001 (cited at
// `ncbi_environment::preparse`: `-version` and the full versions are read before the
// arguments, anywhere before `--`, also where another option expects its value, but not as the
// value of `-logfile` or `-conffile`)
/// The version word (`-version` or a full version) that NCBI's application prints for `tokens`
/// (`ncbi_environment::preparse`), in an option's value too. The command line never brings one
/// here: `main` prints the version before parsing (`ncbi_version_request`). The library API of
/// the search programs and the Web engine, which parse their options here, refuse it
/// (`version_word_error`).
fn ncbi_version_word(tokens: &[OsString]) -> Option<&'static str> {
    use crate::blastinput::ncbi_environment::{preparse, Preparsed};
    match preparse(tokens) {
        Preparsed::Version(request) => Some(request.word()),
        Preparsed::Words { .. } => None,
    }
}

/// The first `-logfile`, `-conffile` or `-dryrun` that NCBI's application reads in `tokens` (with
/// its value; `ncbi_environment::preparse`), as written.
fn toolkit_file_word(tokens: &[OsString]) -> Option<&str> {
    use crate::blastinput::ncbi_environment::{log_file, preparse, Preparsed};
    let (conffile, dryrun) = match preparse(tokens) {
        Preparsed::Words {
            conffile, dryrun, ..
        } => (conffile.is_some(), dryrun),
        _ => (false, false),
    };
    if log_file(tokens).is_some() {
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
/// LOSAT's refusal of the toolkit option `word` (`-logfile`, `-conffile`, `-dryrun`) in
/// `program` when the parser's caller runs no application start: LOSAT Web's ABI v1.
fn toolkit_file_word_error(program: &str, word: &str) -> clap::Error {
    clap::Error::raw(
        ErrorKind::InvalidValue,
        toolkit_file_word_refusal(program, word),
    )
}

/// The text of `toolkit_file_word_error`.
pub(crate) fn toolkit_file_word_refusal(program: &str, word: &str) -> String {
    let name = toolkit_option_name(word);
    let shown = refusal_program_name(program);
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
/// LOSAT's refusal of the version word `word` (`-version` or a full version) in `program` where
/// the parser's caller prints no version: the library API of the search programs and the Web
/// engine's ABI v1.
fn version_word_error(program: &str, word: &str) -> clap::Error {
    clap::Error::raw(ErrorKind::InvalidValue, version_word_refusal(program, word))
}

/// The text of `version_word_error`.
pub(crate) fn version_word_refusal(program: &str, word: &str) -> String {
    let name = toolkit_option_name(word);
    let shown = refusal_program_name(program);
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

/// `program` as the refusals of the toolkit's words name it: the BLAST programs in capitals,
/// makeblastdb as its other refusals do.
fn refusal_program_name(program: &str) -> String {
    match program {
        "makeblastdb" => program.to_string(),
        _ => program.to_ascii_uppercase(),
    }
}

// No NCBI counterpart: NCBI has no library API. NCBI's application prints the XML of `-xmlhelp`
// (ncbiapp.cpp:685-690, cited at `ncbi_xmlhelp`); LOSAT's library refuses it (line RG DECISIONS
// E2x-1, the Owner's choice of 2026-11-13). The command line's entries reach the XML
// (`xml_help_text`, `ncbi_app::run`), and so does a run of the Web engine's ABI v2, which prints it
// on its stream of the command line's standard output (`printed_text`; line WQA's D14 and Q-XML=A,
// the Owner's choice of 2026-11-20, which takes the XML into the Web engine where E2x-1 and E2x-2
// left it out); the library and the Web engine's ABI v1 refuse it.
/// LOSAT's refusal of `-xmlhelp` where no command line's entry prints it: the library API and the
/// Web engine's ABI v1 of `shown` (`refusal_program_name`, `ncbi_app::Tool::refusal_name`).
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

/// The help request that `try_parse_from_start` returns where `-xmlhelp` stopped the reader of
/// `program`: a `DisplayHelp` error that names the program (`xml_help_program`).
fn xml_help_request(program: String) -> clap::Error {
    use clap::error::{ContextKind, ContextValue};
    let mut error = clap::Error::new(ErrorKind::DisplayHelp);
    error.insert(ContextKind::Custom, ContextValue::String(program));
    error
}

/// The program whose `-xmlhelp` the parser's `error` is (`try_parse_from_start`): the command
/// line's entry prints its XML (`xml_help_text`) in place of the error, exit status 0.
pub(crate) fn xml_help_program(error: &clap::Error) -> Option<&str> {
    use clap::error::{ContextKind, ContextValue};
    match error.get(ContextKind::Custom) {
        Some(ContextValue::String(program)) if error.kind() == ErrorKind::DisplayHelp => {
            Some(program)
        }
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:685-697
// NCBI's application writes the help and the XML of `-xmlhelp` to its standard output and exits
// with status 0.
/// What the command line prints to its standard output, exit status 0, where its argument reader
/// stopped with `error` (`cli::main_entry`): the XML of `-xmlhelp` (`xml_help_program`,
/// `xml_help_text`), or the help that stopped the reader (`-h`, `-help`, `-help-full`, `--help`;
/// `command_line_message`, `rpsbproc` as there). `None` for an error, which the command line writes
/// to its standard error. LOSAT Web's ABI v2 prints the same text (stream 11; line WQA's D14,
/// 2026-11-20): one producer of the text for both.
pub fn printed_text(error: &clap::Error, rpsbproc: bool) -> Option<String> {
    xml_help_program(error)
        .and_then(xml_help_text)
        .or_else(|| (!error.use_stderr()).then(|| command_line_message(error, rpsbproc)))
}

/// What `losat <program> -xmlhelp` prints: NCBI's XML of `program` without the arguments that its
/// reader refuses by name (`refuses_option`, `ncbi_xmlhelp::print`).
pub(crate) fn xml_help_text(program: &str) -> Option<String> {
    let command = with_toolkit_file_args(command_of(program), |name| name == program);
    let scope = command.find_subcommand(program)?;
    crate::ncbi_xmlhelp::print(program, |name| refuses_option(scope, name))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2229-2242 (`AddDefaultFileArguments`:
// every NCBI application's arguments have `-logfile` and `-conffile`; the BLAST programs hide
// them from their `-help`, through `HideStdArgs` with `fHideLogfile`, `fHideConffile` and other
// flags, blastn_app.cpp:92)
/// `command` with NCBI's two toolkit file arguments, `-logfile` and `-conffile`, in each NCBI
/// application that runs them (the BLAST programs and makeblastdb) and `selected` takes, hidden as
/// NCBI hides them. The program's own arguments are in `command` (`cli_command`).
/// Their values are read by the application start (`ncbi_environment::log_file` and
/// `ApplicationStart::of_process`), before the arguments: the argument reader only checks them
/// (one of each, as NCBI's `Argument with this name is defined already`).
fn with_toolkit_file_args(
    mut command: clap::Command,
    selected: impl Fn(&str) -> bool,
) -> clap::Command {
    let names: Vec<String> = command
        .get_subcommands()
        .map(|sub| sub.get_name().to_string())
        .filter(|name| TOOLKIT_PROGRAMS.contains(&name.as_str()) || name == "makeblastdb")
        .filter(|name| selected(name.as_str()))
        .collect();
    for name in names {
        command = command.mut_subcommand(name, |mut sub| {
            for (key, id) in [
                ("logfile", "ncbi_toolkit_logfile"),
                ("conffile", "ncbi_toolkit_conffile"),
            ] {
                if !sub.get_arguments().any(|arg| arg.get_long() == Some(key)) {
                    sub = sub.arg(
                        clap::Arg::new(id)
                            .long(key)
                            .value_name("File_Name")
                            .value_parser(clap::value_parser!(OsString))
                            .action(clap::ArgAction::Set)
                            .hide(true),
                    );
                }
            }
            sub
        });
    }
    command
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:78-88
// The standard options of an NCBI C++ Toolkit application, which NCBI's BLAST+ programs
// accept without listing them in -help. `-logfile` and `-conffile` run in the BLAST programs and
// makeblastdb (line RG item E1), `-dryrun` in them and blastdbcmd (item E1b), the three in
// rpsbproc (item E1c, its own `DryRun`: `rpsbproc::dry_run`), `-help-full` in every program (item
// E2, `try_parse_from`), the full versions in every program (item E3, `ncbi_version_request`),
// `-xmlhelp` in every program on the command line (item E2x, `xml_help_program`).

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:44-50
// The options of NCBI tblastn 2.17.0+ (-help) that LOSAT's TBLASTN does not implement
// (AGENTS.md rule 2: explicit unsupported errors). `-best_hit_overhang` / `-best_hit_score_edge`
// are implemented (series NP O6a-3).
const UNPORTED_TBLASTN_ARGS: &[&str] = &[
    "entrez_query",
    "gilist",
    "html",
    "import_search_strategy",
    "negative_gilist",
    "remote",
];

fn is_unported_tblastn_arg(name: &str) -> bool {
    UNPORTED_TBLASTN_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastx_args.cpp:44-50
// The options of NCBI tblastx 2.17.0+ (-help) that LOSAT's TBLASTX does not implement
// (AGENTS.md rule 2: explicit unsupported errors). `-best_hit_overhang` / `-best_hit_score_edge`
// are implemented (series NP O6a-3).
const UNPORTED_TBLASTX_ARGS: &[&str] = &[
    "entrez_query",
    "gilist",
    "html",
    "import_search_strategy",
    "negative_gilist",
    "remote",
];

fn is_unported_tblastx_arg(name: &str) -> bool {
    UNPORTED_TBLASTX_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-66
// Named capabilities exist in NCBI but are explicitly outside the BLASTX product scope.
// `-best_hit_overhang` / `-best_hit_score_edge` are implemented (series NP O6a-3).
const UNPORTED_BLASTX_ARGS: &[&str] = &[
    "html",
    "gilist",
    "negative_gilist",
    "entrez_query",
    "ipglist",
    "negative_ipglist",
    "import_search_strategy",
    "remote",
];

fn is_unported_blastx_arg(name: &str) -> bool {
    UNPORTED_BLASTX_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-230
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
#[derive(Debug)]
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

/// Whether `error` is `NativeError::silent_end`.
pub(crate) fn is_silent_end(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<NativeError>()
        .is_some_and(|native| native.exit == 0 && native.message.is_empty())
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
    /// `CATCH_ALL`'s report of an exception of NCBI's C++ Toolkit (`CException`) that ends a
    /// BLAST program: `Error: `, the exception's `what()` (which ends with a line feed) and the
    /// post's line feed, exit status 255 (`BLAST_UNKNOWN_ERROR`); a `CSeqDBException`, which
    /// `CATCH_ALL` catches before (its branch comes first), `BLAST Database error: ` and its
    /// message, exit status 2 (`BLAST_DATABASE_ERROR`). Every place of a `CSeqDBException` that
    /// LOSAT raises reports the code `eInvalid` (`seqdb::throw_site`), so the texts of that
    /// branch for `eTooManyOpenFiles` and `eMemoryMappingFailure` do not arise here (the reader's
    /// errors of those codes are `common::database::database_error`'s). The one producer of these
    /// reports: an exception of NCBI's that LOSAT raises where NCBI throws it becomes this error.
    pub(crate) fn exception(e: &crate::ncbi_exception::NcbiException) -> NativeError {
        if e.kind().starts_with("CSeqDBException::") {
            return NativeError {
                exit: 2,
                message: [DATABASE_ERROR.as_bytes(), b" ", e.get_msg(), b"\n"].concat(),
            };
        }
        NativeError {
            exit: 255,
            message: [b"Error: ".as_slice(), &e.what(), b"\n"].concat(),
        }
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

    /// The error of a report's writer as the run ends with it: an exception of NCBI's that the
    /// writer raised mid-report (`NcbiException::of_io_error`: TBLASTX's query-anchored section,
    /// `report::query_anchored::write_section`) is `CATCH_ALL`'s report of it (`exception`), after
    /// what the report wrote before it; any other error is the writer's.
    pub(crate) fn of_writer(error: std::io::Error) -> anyhow::Error {
        match crate::ncbi_exception::NcbiException::of_io_error(&error) {
            Some(exception) => NativeError::exception(exception).into(),
            None => error.into(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
/// The report's output stream (the `-out` file or standard output), which records whether
/// a write to it failed: NCBI's outfmt 0 stream throws there, which the application
/// reports as "BLAST failed to write output" (blast_app_util.hpp:252-255).
pub struct ReportStream {
    pub inner: Box<dyn std::io::Write + Send>,
    pub failed: bool,
}

impl ReportStream {
    /// The command line's report stream: the `-out` file `out_file` (`None`: standard output),
    /// behind NCBI's buffer for the organism report (`tax_report`, `report_destination`) and the
    /// hold of a run whose strategy file `strategy` LOSAT refuses (`report_destination`).
    pub fn of_cli(
        out_file: Option<std::io::BufWriter<std::fs::File>>,
        tax_report: bool,
        strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    ) -> Self {
        ReportStream {
            inner: report_destination(
                out_file.map(|file| Box::new(file) as Box<dyn std::io::Write + Send>),
                tax_report,
                strategy,
            ),
            failed: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:252-255
    /// The end of a run whose report went to this stream, `result` being the run's: what the
    /// report wrote stays (a flush, also after an error, such as the outfmt 0 prolog). Where the
    /// report's stream throws (`throws`: the stream of the outfmt 0 to 4 formatter,
    /// blast_format.cpp:118-119, and the organism report's, `NcbiReportBuffer`), a write that
    /// failed (`failed`) or an insertion that set the stream's `badbit` (`BadbitInsertion`) is
    /// NCBI's "BLAST failed to write output", exit status 6. For the other formats NCBI aborts and
    /// LOSAT reports the run's error, then the flush's (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES
    /// exception 3). The command line's runs of blastn, blastp (psiblast), tblastn and tblastx end
    /// their report here, and so does the library's `-out` file (`api::local_blast::OutFile`).
    pub(crate) fn finish(
        &mut self,
        result: anyhow::Result<()>,
        throws: bool,
    ) -> anyhow::Result<()> {
        use anyhow::Context;
        let flushed = std::io::Write::flush(self);
        let badbit = result.as_ref().is_err_and(BadbitInsertion::is);
        if (self.failed || badbit) && throws {
            return Err(NativeError {
                exit: 6,
                message: "BLAST failed to write output\n".into(),
            }
            .into());
        }
        result?;
        flushed.context("failed to write the output")
    }
}

impl std::io::Write for ReportStream {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf);
        self.failed |= written.is_err();
        written
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let flushed = self.inner.flush();
        self.failed |= flushed.is_err();
        flushed
    }
}

// C++ standard, [iostate.flags]: `setstate` calls `clear`, which throws `ios_base::failure` when the bit is in `exceptions()`.
/// NCBI's `ios::failure` of the outfmt 0 to 4 formatter's stream that an insertion of the report
/// raises itself, not a write that the file or device refused (`ReportStream::failed`): an
/// insertion that sets `badbit` (`setstate`, which calls `clear`) throws on a stream whose
/// `exceptions` mask holds `NcbiBadbit` (blast_format.cpp:118-119), and the application reports it
/// as "BLAST failed to write output", exit status 6 (blast_app_util.hpp:252-255), after the bytes
/// the report wrote before it. It travels as the `io::Error` of the report's writes (`error`),
/// which a program reports as a failed write of the stream (`is`); the epilog's null matrix name
/// raises it (`report::pairwise::write_matrix_line`).
#[derive(Debug)]
pub(crate) struct BadbitInsertion;

impl BadbitInsertion {
    /// The failure as the error of the report's writes.
    pub(crate) fn error() -> std::io::Error {
        std::io::Error::other(BadbitInsertion)
    }

    /// Whether `error` is this failure (an `io::Error` that holds it).
    pub(crate) fn is(error: &anyhow::Error) -> bool {
        error
            .downcast_ref::<std::io::Error>()
            .and_then(std::io::Error::get_ref)
            .is_some_and(|inner| inner.is::<BadbitInsertion>())
    }
}

impl std::fmt::Display for BadbitInsertion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("BLAST failed to write output")
    }
}

impl std::error::Error for BadbitInsertion {}

/// The stream of the organism report (outfmt 18) as NCBI's buffer writes it: NCBI's formatter
/// never flushes it, a write to the file or device comes only where the buffer writes
/// (`NcbiBuffer`), and `~CBlastFormat` turns the stream's exceptions off before its last flush
/// (the exit flushes C's `stdout` without a check). So a write that fails where NCBI's buffer
/// writes is "BLAST failed to write output" (`ReportStream::failed`; a reader that closed the
/// pipe: NCBI is ended by SIGPIPE there), and the bytes still held are written when the run ends,
/// a failure then unreported (NCBI BLAST+ 2.17.0, probe `T/t2/fix1` `f1.d5.*`: `-out /dev/full`
/// exit 0 for 2 KB, 6 for 90 KB). One type for both streams: the `-out` file's `basic_filebuf`
/// (`file`) and standard output's glibc buffer (`standard_output`, line-buffered on a terminal),
/// and `cerr` is tied to `cout`: a write to standard error writes what standard output holds first
/// (`tie_standard_output`, `CliStderr`), where NCBI aborts if that write fails (LOSAT: its error,
/// non-zero exit). LOSAT's flushes of the report hold. Each `write` is one insertion of NCBI's
/// (`insert`, `insert_streamed`: the writer of the report makes NCBI's insertions, into the
/// `NcbiOstream` that the run opened, `sink`'s, `standard_output`'s or `file`'s, with nothing
/// between).
pub struct NcbiReportBuffer(std::sync::Arc<std::sync::Mutex<HeldReport>>);

/// What `NcbiReportBuffer` holds.
struct HeldReport {
    inner: Box<dyn std::io::Write + Send>,
    /// The bytes that NCBI's buffer holds.
    held: Vec<u8>,
    /// NCBI's buffer of the stream.
    buffer: NcbiBuffer,
}

impl HeldReport {
    /// NCBI's buffer writes what it holds (a flush).
    fn write_held(&mut self) -> std::io::Result<()> {
        self.buffer.flush();
        let held = std::mem::take(&mut self.held);
        self.inner.write_all(&held)?;
        self.inner.flush()
    }
}

impl Drop for HeldReport {
    fn drop(&mut self) {
        let _ = self.inner.write_all(&self.held);
        let _ = self.inner.flush();
    }
}

/// The organism report's standard output of the run, which a write to standard error writes
/// first (`tie_standard_output`): standard output is the process's.
static TIED_STANDARD_OUTPUT: std::sync::Mutex<
    Option<std::sync::Weak<std::sync::Mutex<HeldReport>>>,
> = std::sync::Mutex::new(None);

fn lock<T>(mutex: &std::sync::Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:317-323
impl NcbiReportBuffer {
    fn new(inner: Box<dyn std::io::Write + Send>, buffer: NcbiBuffer) -> Self {
        NcbiReportBuffer(std::sync::Arc::new(std::sync::Mutex::new(HeldReport {
            inner,
            held: Vec::new(),
            buffer,
        })))
    }

    /// The `-out` file `inner` (`NcbiStdio::File`: `BUFSIZ - 2` bytes held, the 8191st written
    /// with them).
    pub fn file(inner: Box<dyn std::io::Write + Send>) -> Self {
        Self::new(inner, NcbiBuffer::new(NcbiStdio::File))
    }

    /// Standard output, glibc's buffer (`NcbiStdio::cout`, line-buffered on a terminal), tied to
    /// standard error (`tie_standard_output`), writing to `inner` (the process's standard output).
    pub fn standard_output(inner: Box<dyn std::io::Write + Send>) -> Self {
        // NCBI's Windows (MSVC) and macOS (libc++, BSD stdio) builds buffer otherwise: there
        // the exit status of a failing write is unverified (this model is their Linux builds').
        let report = Self::new(inner, NcbiBuffer::new(NcbiStdio::cout()));
        *lock(&TIED_STANDARD_OUTPUT) = Some(std::sync::Arc::downgrade(&report.0));
        report
    }
}

impl std::io::Write for NcbiReportBuffer {
    /// One insertion of NCBI's (`NcbiBuffer::sputn`), whose stream is the bytes held and `buf`:
    /// what NCBI's buffer writes of them is written, the first of them (all that it held, or
    /// none), and it holds the rest.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let mut report = lock(&self.0);
        let report = &mut *report;
        report.held.extend_from_slice(buf);
        report.buffer.sputn(buf.len(), &report.held);
        let written = report.held.len() - report.buffer.held();
        if written > 0 {
            report.inner.write_all(&report.held[..written])?;
            report.inner.flush()?;
            report.held.drain(..written);
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// A stream insertion (`operator<<`) in NCBI's report: one insertion (a string or a number) of
/// `bytes`, one write call, which the organism report's stream takes as one (`NcbiReportBuffer`);
/// other writers just take the bytes. A string that LOSAT streams is `insert_streamed`.
pub fn insert(out: &mut (impl std::io::Write + ?Sized), bytes: &[u8]) -> std::io::Result<()> {
    out.write_all(bytes)
}

/// One insertion of NCBI's of `len` bytes that LOSAT streams instead of holding (the organism
/// report, whose runs of about 4 GiB of spaces `tax_report` keeps as counts): `fill` gives the bytes to
/// the `Insertion`, which writes them in pieces that NCBI's buffers (`NcbiBuffer`) write as they
/// write the whole. Each piece is longer than a buffer (`basic_filebuf` writes it with what it
/// holds and keeps nothing; glibc keeps `(held + n) % block`, which adds up over the pieces as
/// over the whole), and each but the last ends in `NCBI_BUFSIZ` bytes without a newline (so a
/// line-buffered `cout` keeps its whole remainder there, as within the whole); the last is the
/// rest. A string without such a place is written in one piece.
pub fn insert_streamed<W: std::io::Write + ?Sized>(
    out: &mut W,
    len: u64,
    fill: impl FnOnce(&mut Insertion<'_, W>) -> std::io::Result<()>,
) -> std::io::Result<()> {
    let mut insertion = Insertion {
        out,
        piece: Vec::new(),
        left: len,
        since_newline: 0,
    };
    fill(&mut insertion)?;
    debug_assert_eq!(insertion.left, insertion.piece.len() as u64);
    insertion.out.write_all(&insertion.piece)
}

/// The bytes of one insertion that LOSAT streams (`insert_streamed`).
pub struct Insertion<'w, W: std::io::Write + ?Sized> {
    out: &'w mut W,
    /// The piece being made.
    piece: Vec<u8>,
    /// The bytes of the insertion not written yet, `piece`'s among them.
    left: u64,
    /// How many bytes follow the last newline.
    since_newline: u64,
}

impl<W: std::io::Write + ?Sized> Insertion<'_, W> {
    /// The length from which a piece is written where it may end.
    const PIECE: usize = 1 << 16;

    /// The next bytes of the insertion.
    pub fn bytes(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        self.piece.extend_from_slice(bytes);
        self.since_newline = match bytes.iter().rposition(|&byte| byte == b'\n') {
            Some(at) => (bytes.len() - at - 1) as u64,
            None => self.since_newline + bytes.len() as u64,
        };
        self.end_piece()
    }

    /// The next `count` bytes of the insertion, spaces.
    pub fn spaces(&mut self, mut count: u64) -> std::io::Result<()> {
        while count > 0 {
            let n = count.min(Self::PIECE as u64);
            self.piece.resize(self.piece.len() + n as usize, b' ');
            self.since_newline += n;
            count -= n;
            self.end_piece()?;
        }
        Ok(())
    }

    /// Writes the piece where it may end (`insert_streamed`) once it is `PIECE` bytes long. A
    /// piece ends only if more than `NCBI_BUFSIZ` bytes follow it, which is at least a block (a
    /// block never exceeds `NCBI_BUFSIZ`), so the last piece always takes the buffer's branch of
    /// an insertion longer than the room, as the whole insertion does; a final piece that exactly
    /// filled a block would differ from the whole.
    fn end_piece(&mut self) -> std::io::Result<()> {
        let n = self.piece.len() as u64;
        if n >= Self::PIECE as u64
            && self.since_newline >= NCBI_BUFSIZ as u64
            && self.left - n > NCBI_BUFSIZ as u64
        {
            self.out.write_all(&self.piece)?;
            self.left -= n;
            self.piece.clear();
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68-93
/// The stream of a report as NCBI's formatter writes to it (`CBlastFormat::m_Outfile`), with
/// nothing of a program's before it: each write is one insertion (`operator<<`) of NCBI's
/// (`insert`, `insert_streamed`), passed on as one write. The organism report (`tax_report`)
/// writes only to this type, and a writer becomes one only through the constructors below, one per
/// destination, each saying what comes after it: they are the only places where a buffer or hold
/// could come between the report's insertions and NCBI's buffer (`NcbiReportBuffer`, which counts
/// each write as one insertion). A format's sink (`sink`: the command line's `ReportStream` of
/// blastn, blastp, tblastn and tblastx, the library's and LOSAT Web's writers), blastx's command
/// line (`standard_output` and `file` for the organism report, `buffered` for every other format),
/// and the hold of a library call (`held`).
pub struct NcbiOstream<'a>(Box<dyn std::io::Write + 'a>);

impl<'a> NcbiOstream<'a> {
    fn new(stream: Box<dyn std::io::Write + 'a>) -> Self {
        NcbiOstream(stream)
    }

    /// A format's sink opened as its report's stream (`OutputSink::open`): the caller's writer as
    /// it is, which the command line of blastn, blastp, tblastn and tblastx makes its
    /// `ReportStream` (the organism report's is NCBI's buffer, `report_destination`), which the
    /// library's `run_to_file` makes the same (`api::local_blast::OutFile`) and which a library
    /// caller or LOSAT Web gives (no stream of NCBI's to model); or, for a library call's `Stdout`
    /// and `File` sinks, LOSAT's buffered standard output or file (no stream of NCBI's either: the
    /// report's bytes are NCBI's, the times at which they reach the device are not).
    pub(crate) fn sink(
        sink: &'a mut crate::api::local_blast::OutputSink<'_>,
    ) -> std::io::Result<Self> {
        sink.open().map(Self::new)
    }

    /// The organism report's standard output on blastx's command line: NCBI's buffer of `cout`
    /// (`NcbiReportBuffer::standard_output`) before the process's standard output
    /// (`report_standard_output`). Nothing comes before the buffer or between it and the device.
    pub(crate) fn standard_output(
        strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    ) -> Self {
        Self::new(report_destination(None, true, strategy))
    }

    /// The organism report's `-out` file `file` on blastx's command line: NCBI's `basic_filebuf`
    /// (`NcbiReportBuffer::file`) before the file, as `standard_output`'s.
    pub(crate) fn file(
        file: std::fs::File,
        strategy: Option<&crate::report::archive::strategy::StrategyExport>,
    ) -> Self {
        Self::new(report_destination(Some(Box::new(file)), true, strategy))
    }

    /// The stream of a format other than the organism report on blastx's command line
    /// (`native::run_inner`): LOSAT's buffer of 8 KiB (`BufWriter`) before `destination` (the file
    /// or standard output, behind the held document, `blast_xml::HeldOutput`). NCBI's buffers of these formats
    /// are not modelled (the order of their texts against the report is item T8 of line RG); the
    /// organism report never gets this stream (`run_inner` chooses by the format).
    pub(crate) fn buffered(destination: Box<dyn std::io::Write + 'a>) -> Self {
        Self::new(Box::new(std::io::BufWriter::new(destination)))
    }

    /// A format's stream (`sink`'s) behind the hold of a library call (`local::report_batches`,
    /// `blast_xml::HeldFormats`), which passes the writes and flushes on in their order (at the
    /// release, while it holds). The caller's writer has no stream of NCBI's (`sink`), so no buffer
    /// of NCBI's comes after the hold.
    pub(crate) fn held(
        stream: crate::algorithm::blastx::HeldFormatWriter<NcbiOstream<'a>>,
    ) -> Self {
        Self::new(Box::new(stream))
    }
}

impl std::io::Write for NcbiOstream<'_> {
    /// One insertion, written with one write of the stream.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.write_all(buf)?;
        Ok(buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.0.write_all(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

/// `cerr.tie()`: before a write to standard error, the organism report's standard output writes
/// what it holds (`NcbiReportBuffer::standard_output`); nothing for every other report.
pub fn tie_standard_output() -> std::io::Result<()> {
    let tied = lock(&TIED_STANDARD_OUTPUT)
        .as_ref()
        .and_then(std::sync::Weak::upgrade);
    match tied {
        Some(report) => lock(&report).write_held(),
        None => Ok(()),
    }
}

/// The command line's standard error, tied to the organism report's standard output
/// (`tie_standard_output`): the process's (`process`), so that what standard output writes before
/// a text is written before it where the text is made. A whole message is written with one
/// `write_all` of the stream, which `std::io::Stderr` writes under one lock of standard error, so
/// that no other thread's message comes inside it. An empty write is no write of NCBI's to `cerr`
/// (a query without warnings, whose `ERR_POST` NCBI does not make): it writes nothing of the
/// report. Where the standard error is the file of `-logfile` (tblastx's, rpsblast's and
/// rpstblastn's run, `DiagTarget::into_standard_error`), a failed write is dropped, as NCBI's
/// handler of that file drops it: `CFileHandleDiagHandler::Post` ignores what the write returns
/// (NCBI BLAST+ 2.17.0 tblastx `-logfile /dev/full` runs and writes its report; line WQA's D11
/// round 2, 2026-11-20).
/// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:6695-6696
pub struct CliStderr(std::io::Stderr);

impl CliStderr {
    /// `result` of a write, dropped where the standard error is the file of `-logfile`.
    fn logged<T: Default>(result: std::io::Result<T>) -> std::io::Result<T> {
        match result {
            Err(_) if crate::blastinput::ncbi_environment::terminal_stderr().is_some() => {
                Ok(T::default())
            }
            result => result,
        }
    }
}

impl CliStderr {
    /// The process's standard error.
    pub fn process() -> Self {
        CliStderr(std::io::stderr())
    }
}

impl std::io::Write for CliStderr {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        tie_standard_output()?;
        Self::logged(self.0.write_all(buf)).map(|()| buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        if buf.is_empty() {
            return Ok(());
        }
        tie_standard_output()?;
        Self::logged(self.0.write_all(buf))
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        self.write_all(args.to_string().as_bytes())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

/// The report stream's writer of a run: the `-out` file `file` (`None`: standard output), behind
/// NCBI's buffer for the organism report (`tax_report`, `NcbiReportBuffer`) and, nearest the
/// destination, the hold of a run whose strategy file `strategy` LOSAT refuses
/// (`strategy::hold_report_send`).
pub fn report_destination(
    file: Option<Box<dyn std::io::Write + Send>>,
    tax_report: bool,
    strategy: Option<&crate::report::archive::strategy::StrategyExport>,
) -> Box<dyn std::io::Write + Send> {
    use crate::report::archive::strategy::hold_report_send;
    match (file, tax_report) {
        (Some(file), true) => Box::new(NcbiReportBuffer::file(hold_report_send(strategy, file))),
        (Some(file), false) => hold_report_send(strategy, file),
        (None, true) => Box::new(NcbiReportBuffer::standard_output(hold_report_send(
            strategy,
            report_standard_output(),
        ))),
        (None, false) => hold_report_send(strategy, report_standard_output()),
    }
}

/// Standard output as the report's stream (`ReportStream::inner`).
///
/// NCBI writes the report to `cout`; when the process starts with its standard output
/// closed (`>&-`), its first write fails (outfmt 0: "BLAST failed to write output", exit
/// 6; outfmt 6 and 7 abort). Rust's runtime opens `/dev/null` on a closed standard
/// descriptor before `main`, which LOSAT cannot tell from a `/dev/null` that the caller
/// opened (read and write, as Python's `subprocess.DEVNULL`), so the report is written
/// there and the run succeeds (approved exception 6 of
/// `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, accepted by the maintainer in session S08b).
pub fn report_standard_output() -> Box<dyn std::io::Write + Send> {
    Box::new(std::io::BufWriter::new(std::io::stdout()))
}

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
/// exits with 1 and writes nothing in NCBI BLAST+ 2.17.0; line WQA's D11 round 2, 2026-11-20).
pub fn exit_on_native_error(error: &anyhow::Error) {
    if let Some(e) = error.downcast_ref::<NativeError>() {
        if !posted(&e.message) {
            let mut stderr = std::io::stderr();
            let _ = std::io::Write::write_all(&mut stderr, catch_all_text(&e.message));
            let _ = std::io::Write::flush(&mut stderr);
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
        let _ = std::io::Write::write_all(&mut std::io::stderr(), catch_all_text(&native.message));
        *POSTED_TEXT
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(native.message.clone());
    }
}

#[cfg(test)]
mod manual_refusals_tests;

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
