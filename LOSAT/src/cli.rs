//! Public CLI v2 boundary. Clap's double-dash representation is internal only.

use std::{ffi::OsString, fmt};

use clap::{error::ErrorKind, CommandFactory, Parser, Subcommand};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
// ```c++
//     static const string kProgram("blastx");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                   "Translated Query-Protein Subject BLAST"));
// ```
use crate::algorithm::{blastn, blastp, blastx, psiblast, rpsblast, rpstblastn, tblastn, tblastx};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46,51,75
// ```c++
// const string kArgQuery("query");
// const string kArgSubject("subject");
// const string kArgNumThreads("num_threads");
// ```
// LOSAT retains its program subcommand; search parameter names are NCBI names.
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3815-3826
// ```c++
// string& CArgDescriptions::PrintUsage(string& str, bool detailed) const
// {
//     CPrintUsage x(*this);
//     list<string> arr;
//
//     // SYNOPSIS
//     arr.push_back("USAGE");
//     x.AddSynopsis(arr, m_UsageName, "    ");
//
//     // DESCRIPTION
//     arr.push_back(kEmptyStr);
//     x.AddDescription(arr, detailed);
// ```
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

#[derive(Parser, Debug)]
#[command(
    name = "losat",
    version,
    about = "Pure-Rust pairwise alignment that matches NCBI BLAST+ on LOSAT's test fixtures, with genetic-code exceptions",
    help_template = LOSAT_HELP_TEMPLATE,
    after_help = "To call the programs by NCBI's names (blastn, makeblastdb, ...), link them with:\n  losat -install-shims [DIR]   (default: the directory of this executable)"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:47-49
    // ```c++
    //     static const string kProgram("blastx");
    //     arg.Reset(new CProgramDescriptionArgs(kProgram,
    //                                   "Translated Query-Protein Subject BLAST"));
    // ```
    /// Translated nucleotide query vs protein subject (local blastx)
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Blastx(blastx::BlastxArgs),
    /// Pairwise nucleotide alignment (megablast [default], blastn)
    #[allow(rustdoc::broken_intra_doc_links, reason = "clap prints it as help")]
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Blastn(blastn::BlastnArgs),
    /// Pairwise protein alignment (blastp)
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Blastp(blastp::BlastpArgs),
    /// Pairwise 6-frame translated nucleotide alignment (tblastx)
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Tblastx(tblastx::TblastxArgs),
    // NCBI c++/src/algo/blast/blastinput/tblastn_args.cpp:47-52:
    // static const string kProgram("tblastn"); SetTask("tblastn");
    /// Protein query vs translated nucleotide subject (local tblastn)
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Tblastn(tblastn::TblastnArgs),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:47-51
    // ```c++
    //     static const string kProgram("psiblast");
    //     CRef<IBlastCmdLineArgs> arg;
    //     arg.Reset(new CProgramDescriptionArgs(kProgram,
    //                                           "Position-Specific Iterated BLAST"));
    // ```
    /// Position-Specific Iterated BLAST (psiblast; its first round)
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Psiblast(psiblast::PsiblastArgs),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:77-79
    // ```c++
    //     static const string kProgram("rpsblast");
    //     arg.Reset(new CProgramDescriptionArgs(kProgram,
    //                                           "Reverse Position Specific BLAST"));
    // ```
    /// Reverse Position Specific BLAST (rpsblast)
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Rpsblast(rpsblast::RpsblastArgs),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:48-50
    // ```c++
    //     static const string kProgram("rpstblastn");
    //     arg.Reset(new CProgramDescriptionArgs(kProgram,
    //                                "Translated Reverse Position Specific BLAST"));
    // ```
    // `-h`, `-help` and `-help-full` write NCBI's frozen texts (`RpstblastnArgs::ncbi_help`; the
    // reader stops at them, `try_parse_from_command`).
    /// Translated Reverse Position Specific BLAST (rpstblastn)
    #[command(help_template = NCBI_HELP_TEMPLATE)]
    Rpstblastn(rpstblastn::RpstblastnArgs),
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:217-219
    // ```c++
    //     arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
    //                   "Application to create BLAST databases, version "
    //                   + CBlastVersion().Print());
    // ```
    /// makeblastdb: create a BLAST database from a FASTA file, format version 4
    #[command(name = "makeblastdb", help_template = NCBI_HELP_TEMPLATE)]
    Makedb(crate::writedb::makeblastdb::MakedbArgs),
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:64-65
    // ```c++
    // static const char *PROGRAM_USAGE =	//multiline
    // "Process output from offline rpsblast to annotate sequence with conserved domain information";
    // ```
    /// Process output from offline rpsblast to annotate sequence with conserved domain information (rpsbproc; its domain hits, sites and motifs of protein and nucleotide queries)
    #[command(name = "rpsbproc", help_template = NCBI_HELP_TEMPLATE)]
    Rpsbproc(crate::rpsbproc::RpsbprocArgs),
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:166-170,203-207,332-349
// ```c++
// arg_desc.AddOptionalKey(kArgWordSize, "int_value", description, CArgDescriptions::eInteger);
// arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value", ...);
// arg_desc.AddDefaultKey(kArgSegFiltering, "SEG_options", ..., CArgDescriptions::eString, ...);
// ```
// Translate registered canonical keys only. Consume each value once so a path
// or a quoted filtering specification can never be interpreted as another key.
/// Parses `argv` (`losat <program> <words>`) for a caller that runs no NCBI application start
/// (the Web engine): the toolkit's words that the start reads (`-logfile`, `-conffile`,
/// `-dryrun`) are refused (`toolkit_file_word_error`; makeblastdb keeps its own `-logfile`
/// argument), and so are the versions (`version_word_error`) and `-xmlhelp`
/// (`xml_help_refusal`).
pub fn try_parse_from<T, I, S>(argv: I) -> Result<T, clap::Error>
where
    T: Parser + CommandFactory,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    try_parse_from_command(T::command(), argv)
}

/// `try_parse_from` for the library API, whose calls run an NCBI application start of their own
/// (`api::start`): `-logfile`, `-conffile` and `-dryrun` are taken; the versions and `-xmlhelp`,
/// whose texts only the command line writes, are refused (`version_word_error`,
/// `xml_help_refusal`).
pub(crate) fn try_parse_from_library<T, I, S>(argv: I) -> Result<T, clap::Error>
where
    T: Parser + CommandFactory,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    parse_command(T::command(), argv, Reader::Library)
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
    T: Parser + CommandFactory,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    parse_command(T::command(), argv, Reader::CommandLine(start))
}

/// `try_parse_from` with `command`, one of `T`'s commands whose arguments a caller changed
/// (Web ABI v1's BLASTX keeps its own `-outfmt` values, `blastx/web.rs`), read as clap's
/// `Parser::try_parse_from` reads `T::command()`.
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
    /// A caller that runs no NCBI application start (the Web engine): the toolkit's words, the
    /// versions and `-xmlhelp` are refused.
    NoStart,
    /// The library API, whose calls run a start of their own (`api::start`): `-logfile`,
    /// `-conffile` and `-dryrun` are taken; the versions and `-xmlhelp` are refused.
    Library,
    /// The command line's entry with its application start: `-xmlhelp` is the request whose XML
    /// the entry prints.
    CommandLine(&'a crate::blastinput::ncbi_environment::ApplicationStart),
}

/// The reader of `try_parse_from`, `try_parse_from_library`, `try_parse_from_start` and
/// `try_parse_from_command`.
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
    let command = with_toolkit_file_args(command);
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
            scope = subcommand;
            translated.push(token);
            if scope.get_name() == "rpstblastn" {
                ncbi_help = Some((translated.len(), ""));
            }
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3020-3023, 2840-2852
        // ```c++
        //     const CArgDesc& arg = **it;
        //
        //     if ( s_IsFlag(arg) ) {
        //         x_CheckAutoHelp(arg1);
        // ...
        // void CArgDescriptions::x_CheckAutoHelp(const string& arg) const
        // {
        //     if (arg.compare(string("-") + s_AutoHelp) == 0) {
        //         if (m_AutoHelp) {
        //             NCBI_THROW(CArgHelpException,eHelp,kEmptyStr);
        //         }
        //     } else if (arg.compare(string("-") + s_AutoHelpFull) == 0) {
        //         NCBI_THROW(CArgHelpException,eHelpFull,kEmptyStr);
        //     } else if (arg.compare(string("-") + s_AutoHelpXml) == 0) {
        //         NCBI_THROW(CArgHelpException,eHelpXml,kEmptyStr);
        //     } else if (arg.compare(string("-") + s_AutoHelpShowAll) == 0) {
        //         NCBI_THROW(CArgHelpException,eHelpShowAll,kEmptyStr);
        //     }
        // ```
        // NCBI's reader stops at its help flags, the words before them read first: a word after
        // one is never read (NCBI BLAST+ 2.17.0: `-help -bogus`, `-help ... -out` with no value
        // and `-h -xmlhelp` print the help, exit status 0; line RG item E2). Every program of
        // this parser follows NCBI's reader; rpstblastn writes NCBI's texts
        // (`RpstblastnArgs::ncbi_help`). An option's value is not a flag (`-out -help` names a
        // file). `--help` is LOSAT's spelling of `-help` (clap's) and stops the reader too.
        // `-help-full` (`s_AutoHelpShowAll`) is `-help` with the arguments defined `fHidden`
        // (`ShowAllArguments`). Of NCBI's programs only makeblastdb has one, `-randomize`, which
        // LOSAT refuses; LOSAT's help lists the arguments LOSAT runs (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES
        // exception 1), so its text is `-help`'s, as NCBI's other programs write the bytes of
        // `-help` (in every place: line RG, `E/SURVEY_E2.md`).
        let stops_at_help = !std::ptr::eq(scope, &command);
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:685-690 (quoted at
        // `ncbi_xmlhelp`)
        // `-xmlhelp` (`s_AutoHelpXml`) stops the reader as `-help` does, the words before it
        // read: the command line's entry prints the program's arguments in NCBI's XML
        // (`xml_help_program`, `ncbi_xmlhelp::print`; line RG item E2x). The library API and the
        // Web engine, which write no help of NCBI's, refuse it there, after the words before it
        // (`xml_help_refusal`).
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
        // ```c++
        // static const char* s_AutoHelp     = "h";
        // static const char* s_AutoHelpFull = "help";
        // ```
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
        // ```c++
        //     if (*n_plain == kMax_UInt || m_PositionalMode == ePositionalMode_Loose) {
        //         // Check for the s_ArgDelimiter delimiter
        //         if (arg1.compare(s_ArgDelimiter) == 0) {
        //             if (*n_plain == kMax_UInt) {
        //                 *n_plain = 0;  // pos.args started
        //             }
        //             return false;
        // ```
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
        // ```c++
        // 				if (TArgDefinition::EMPTYSTR != argdef.alias)
        // 					pMyArgs->AddAlias(argdef.alias, argdef.name);
        // ```
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
            // NCBI reference: c++/src/corelib/ncbiargs.cpp:266-270
            // ```c++
            // inline CArg_String::CArg_String(const string& name, const string& value)
            //     : CArgValue(name)
            // {
            //     m_StringList.push_back(value);
            // }
            // ```
            // NCBI keeps the bytes of an option's value. LOSAT's options are UTF-8 strings,
            // so a value that is not UTF-8 is rejected explicitly, after the parser: NCBI's
            // USAGE errors and `-help` come first, as for any value. The parser reads the
            // value with U+FFFD in place of the bytes. The file names keep their bytes and
            // have their own check where the files are opened (`check_utf8_file_name`).
            let value = if value.to_str().is_none()
                && !matches!(name, "query" | "subject" | "out")
                && is_ncbi_search_program(scope)
            {
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
                Reader::CommandLine(_) => xml_help_request(program),
                Reader::NoStart | Reader::Library => clap::Error::raw(
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

/// The explicit rejection of an option's value that is not UTF-8 (`try_parse_from`).
///
/// NCBI reference: c++/src/corelib/ncbiargs.cpp:266-270
/// ```c++
/// inline CArg_String::CArg_String(const string& name, const string& value)
///     : CArgValue(name)
/// {
///     m_StringList.push_back(value);
/// }
/// ```
fn non_utf8_value_error(scope: &clap::Command, name: &str) -> clap::Error {
    clap::Error::raw(
        ErrorKind::InvalidUtf8,
        format!(
            "the value of -{name} is not UTF-8; NCBI BLAST+ reads the bytes of an option's value, which is not supported by LOSAT's {}",
            scope.get_name().to_uppercase()
        ),
    )
}

// NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:107,143
// ```c++
// const string kArgWordSize("word_size");
// const string kArgCompBasedStats("comp_based_stats");
// ```
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
    let command = with_toolkit_file_args(Cli::command());
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
    message
        .replace("'-h'", "'-help'")
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version")
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

/// `render_rpsbproc_message` of a part of a message.
fn render_rpsbproc_text(text: &str) -> String {
    let mut message = text
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version");
    let command = Cli::command();
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
// ```c++
// m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
// arg.Reset(new CGenericSearchArgs(kQueryIsProtein));
// m_HspFilteringArgs.Reset(new CHspFilteringArgs);
// m_FormattingArgs.Reset(new CFormattingArgs);
// m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eNucleotideDb));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:63-70
// ```c++
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
//     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
//     arg.Reset(m_BlastDbArgs);
//     m_Args.push_back(arg);
//
//     m_StdCmdLineArgs.Reset(new CStdCmdLineArgs);
//     arg.Reset(m_StdCmdLineArgs);
//     m_Args.push_back(arg);
// ```
// The options of NCBI blastn 2.17.0+ (-help) that LOSAT's BLASTN does not implement
// (AGENTS.md rule 2: explicit unsupported errors). `-db` and `-dbsize` are implemented
// (`blastn/args.rs`), as are the X-drops, `-window_size`, `-off_diagonal_range`, `-no_greedy`
// (DI-6), `-min_raw_gapped_score` (DI-6b) and the database masks (`blastinput::db_mask`, series
// NP M4a-2).
const UNPORTED_BLASTN_ARGS: &[&str] = &[
    "best_hit_overhang",
    "best_hit_score_edge",
    "entrez_query",
    "filtering_db",
    "gilist",
    "html",
    "import_search_strategy",
    "index_name",
    "negative_gilist",
    "remote",
];

fn is_unported_blastn_arg(name: &str) -> bool {
    UNPORTED_BLASTN_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:44-50
// ```c++
// CBlastpAppArgs::CBlastpAppArgs()
// {
//     CRef<IBlastCmdLineArgs> arg;
//     static const string kProgram("blastp");
//     arg.Reset(new CProgramDescriptionArgs(kProgram, "Protein-Protein BLAST"));
//     const bool kQueryIsProtein = true;
//     bool const kFilterByDefault = false;
// ```
// The options of NCBI blastp 2.17.0+ (-help) that LOSAT's BLASTP does not implement
// (AGENTS.md rule 2: explicit unsupported errors).
const UNPORTED_BLASTP_ARGS: &[&str] = &[
    "best_hit_overhang",
    "best_hit_score_edge",
    "db_hard_mask",
    "db_soft_mask",
    "entrez_query",
    "gilist",
    "html",
    "import_search_strategy",
    "ipglist",
    "negative_gilist",
    "negative_ipglist",
    "remote",
    "subject_besthit",
];

fn is_unported_blastp_arg(name: &str) -> bool {
    UNPORTED_BLASTP_ARGS.contains(&name)
}

// The rejection of an option that a program does not have: NCBI's options that LOSAT
// has not ported, the NCBI C++ Toolkit options, and the rest (unknown to NCBI too).
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1156-1164
// ```c++
//         // checkpoint file
//         arg_desc.AddOptionalKey(kArgPSIOutputChkPntFile, "checkpoint_file",
//
//                                 "File name to store checkpoint file",
//                                 CArgDescriptions::eOutputFile);
//         // ASCII matrix file
//         arg_desc.AddOptionalKey(kArgAsciiPssmOutputFile, "ascii_mtx_file",
//                                 "File name to store ASCII version of PSSM",
//                                 CArgDescriptions::eOutputFile);
// ```
/// psiblast's NCBI arguments that LOSAT does not run: the MSA inputs (series PS6; the ASCII and
/// ASN.1 PSSM files are PS2's and PS3's, `algorithm::psiblast::pssm`, `-in_pssm` PS4's,
/// `algorithm::psiblast::in_pssm`), PHI-BLAST (PS8), and the arguments that blastp's engine does
/// not run either (`is_unported_blastp_arg`; psiblast has no database masks and no `-mt_mode`).
/// The search strategy file of a FASTA query is written (item F1b-1 of line RG; with `-in_pssm`
/// it is refused, `PsiblastArgs::into_blastp_args`).
const UNPORTED_PSIBLAST_OWN_ARGS: &[&str] = &[
    "ignore_msa_master",
    "in_msa",
    "msa_master_idx",
    "phi_pattern",
];

fn is_unported_psiblast_arg(name: &str) -> bool {
    UNPORTED_PSIBLAST_OWN_ARGS.contains(&name)
        || (is_unported_blastp_arg(name) && !matches!(name, "db_hard_mask" | "db_soft_mask"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:85-127 (quoted at
// `algorithm::rpsblast::args`)
/// rpsblast's NCBI arguments that LOSAT does not run: `-use_sw_tback` (the Smith-Waterman
/// redo, PS11c-3), `-html` (PS15a) and the remote options and `-import_search_strategy` of every
/// program.
const UNPORTED_RPSBLAST_ARGS: &[&str] = &[
    "entrez_query",
    "html",
    "import_search_strategy",
    "remote",
    "use_sw_tback",
];

fn is_unported_rpsblast_arg(name: &str) -> bool {
    UNPORTED_RPSBLAST_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:44-110 (quoted at
// `algorithm::rpstblastn::args`)
/// rpstblastn's NCBI arguments that LOSAT does not run: `-html` (PS15a) and the remote options
/// and `-import_search_strategy` of every program. The composition-based statistics of a later
/// item are refused once per run (`rpstblastn::search::TranslatedRpsSearch::check_options`).
fn is_unported_rpstblastn_arg(name: &str) -> bool {
    matches!(
        name,
        "entrez_query" | "html" | "import_search_strategy" | "remote"
    )
}

fn unknown_option_error(scope: &clap::Command, name: &str, text: &str) -> clap::Error {
    named_refusal(scope.get_name(), name).unwrap_or_else(|| {
        clap::Error::raw(
            ErrorKind::UnknownArgument,
            format!("unknown option or argument '{text}'; use -help for CLI v2 syntax"),
        )
    })
}

/// Whether the reader of the program `scope` (of `with_toolkit_file_args(Cli::command())`)
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-918 (quoted at
// blastp's check, `validate_requested_blastp_support`)
// ```c++
//     arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
// ```
/// NCBI's arguments that a program's parser reads, so that NCBI's own errors about them and the
/// other arguments come first (blastp's `-ungapped` with composition-based statistics:
/// `Composition-adjusted searched are not supported with an ungapped search`), and whose run LOSAT
/// refuses whatever their value, by name, where NCBI starts the search: the program, the name and
/// what the refusal says the argument is (`run_refusal`). `-xmlhelp` leaves them out of NCBI's XML
/// as it leaves the parser's refusals (`refuses_option`).
pub(crate) const RUN_REFUSALS: &[(&str, &str, &str)] = &[
    ("blastp", "ungapped", "an ungapped search"),
    ("blastp", "use_sw_tback", "the Smith-Waterman traceback"),
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
    // NCBI c++/src/algo/blast/blastinput/tblastn_args.cpp:64-129:
    // m_Args.push_back(arg) registers shared, formatting, DB, and PSI groups.
    // The named NCBI options below have no Stage B Rust behavior.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-66
    // ```c++
    //     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
    //     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
    //     m_BlastDbArgs->SetIPGFilteringSupport(true);
    //     arg.Reset(m_BlastDbArgs);
    //     m_Args.push_back(arg);
    // ```
    // Product scope explicitly refuses unported NCBI capabilities.
    if program == "blastx" && is_unported_blastx_arg(name) {
        return Some(clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("unsupported BLASTX option '-{name}': outside the declared local FASTA scope"),
        ));
    }
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:246-250
    // ```c++
    //     arg_desc->AddFlag("parse_seqids",
    //                       "Option to parse seqid for FASTA input if set, for all other input types seqids are parsed automatically", true);
    //     arg_desc->AddFlag("hash_index",
    //                       "Create index of sequence hash values.",
    // ```
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
/// `-version` prints the two lines of `fVersionInfo | fPackageShort`: the program name and
/// version, then the package. NCBI's second line names its own package and build date; LOSAT's
/// names LOSAT and the BLAST+ release its output follows. The full versions are NCBI's forms
/// with LOSAT's values (`ncbi_version::print`, line RG item E3).
/// `makeblastdb` prints `makeblastdb`, as NCBI's program does. `main` prints the
/// text on stdout and exits 0, before any argument is parsed. The top level `losat --version`
/// is clap's and is not read here.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:972-977
/// ```c++
///             } else if ( NStr::strcmp(argv[i], s_ArgVersion) == 0 ) {
///                 delete[] v;
///                 // Print VERSION
///                 cout << GetFullVersion().Print( appname, CVersionAPI::fVersionInfo | CVersionAPI::fPackageShort );
///                 diag_context.DiscardMessages();
///                 return 0;
/// ```
pub fn ncbi_version_request<I, S>(argv: I) -> Option<String>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let tokens: Vec<OsString> = argv.into_iter().map(Into::into).collect();
    let command = Cli::command();
    let program = command.find_subcommand(tokens.get(1)?.to_str()?)?;
    let name = match program.get_name() {
        "blastn" => "blastn",
        "blastp" => "blastp",
        "blastx" => "blastx",
        "tblastn" => "tblastn",
        "tblastx" => "tblastx",
        "psiblast" => "psiblast",
        "rpsblast" => "rpsblast",
        "rpstblastn" => "rpstblastn",
        "makeblastdb" => "makeblastdb",
        "rpsbproc" => "rpsbproc",
        _ => return None,
    };
    // NCBI's reading of the words before the arguments (`ncbi_environment::preparse`): not
    // after `--`, nor as the value of `-logfile` or `-conffile` (`blastn -logfile -version`
    // logs to the file `-version`).
    match crate::blastinput::ncbi_environment::preparse(&tokens[2..]) {
        crate::blastinput::ncbi_environment::Preparsed::Version(request) => {
            // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:388
            // ```c++
            // 	SetVersion(CVersionInfo(VER.major, VER.minor, VER.patch_level));
            // ```
            // rpsbproc's own version (0.5.1), not BLAST's; the second line names the BLAST+
            // toolkit that bioconda builds it on (orchestrator's choice, Owner-delegated
            // 2026-11-11).
            let version = if name == "rpsbproc" {
                crate::rpsbproc::VERSION
            } else {
                crate::NCBI_BLAST_VERSION
            };
            Some(crate::ncbi_version::print(request, name, version))
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

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1001 (quoted at
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
// command line (`cli::main_entry`, `ncbi_environment::ApplicationStart`) and of each library
// call (`api::start`) reads `-logfile` and `-conffile` before the arguments, and both end a
// `-dryrun` after them; the Web engine runs none, so it does not take them (line RG item E1w).
/// LOSAT's refusal of the toolkit option `word` (`-logfile`, `-conffile`, `-dryrun`) in
/// `program` when the parser's caller runs no application start: the Web engine.
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
        LibraryPart::Search | LibraryPart::Tool => "the command line and in the library API",
        LibraryPart::Absent => "the command line",
    };
    format!(
        "the NCBI C++ Toolkit option -{name} runs on {places} of LOSAT's {shown} only: the Web engine runs no NCBI application start, which reads it"
    )
}

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
// (ncbiapp.cpp:926-1001, quoted at `ncbi_environment::preparse`). The library API of the search
// programs returns reports, and the Web engine shows them: neither has a place for the text
// (line RG DECISIONS E3-2; the library's tools print it to their output, `ncbi_app::run_with`,
// `api::run_makeblastdb`, so their refusal names the Web engine alone).
/// LOSAT's refusal of the version word `word` (`-version` or a full version) in `program` where
/// the parser's caller prints no version: the library API of the search programs and the Web
/// engine.
fn version_word_error(program: &str, word: &str) -> clap::Error {
    clap::Error::raw(ErrorKind::InvalidValue, version_word_refusal(program, word))
}

/// The text of `version_word_error`.
pub(crate) fn version_word_refusal(program: &str, word: &str) -> String {
    let name = toolkit_option_name(word);
    let shown = refusal_program_name(program);
    let (places, refusers) = match library_part(program) {
        LibraryPart::Search => ("the command line", "the library API and the Web engine do"),
        LibraryPart::Tool => (
            "the command line and in the library API",
            "the Web engine does",
        ),
        LibraryPart::Absent => ("the command line", "the Web engine does"),
    };
    format!(
        "the NCBI C++ Toolkit option -{name} runs on {places} of LOSAT's {shown} only: {refusers} not write NCBI's version text"
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
// (ncbiapp.cpp:685-690, quoted at `ncbi_xmlhelp`); LOSAT's library and the Web engine refuse it
// (line RG DECISIONS E2x-1, the Owner's choice of 2026-11-13). Only the command line's entries
// reach the XML (`xml_help_text`, `ncbi_app::run`), so the library and the Web engine do not carry
// it (E2x-2).
/// LOSAT's refusal of `-xmlhelp` where no command line's entry prints it: the library API and the
/// Web engine of `shown` (`refusal_program_name`, `ncbi_app::Tool::refusal_name`).
pub(crate) fn xml_help_refusal(shown: &str) -> String {
    format!(
        "the NCBI C++ Toolkit option -xmlhelp runs on the command line of LOSAT's {shown} only: the library API and the Web engine do not write NCBI's XML description of the arguments"
    )
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

/// What `losat <program> -xmlhelp` prints: NCBI's XML of `program` without the arguments that its
/// reader refuses by name (`refuses_option`, `ncbi_xmlhelp::print`).
pub(crate) fn xml_help_text(program: &str) -> Option<String> {
    let command = with_toolkit_file_args(Cli::command());
    let scope = command.find_subcommand(program)?;
    crate::ncbi_xmlhelp::print(program, |name| refuses_option(scope, name))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2337-2357 (`AddDefaultFileArguments`:
// every NCBI application's arguments have `-logfile` and `-conffile`; the BLAST programs hide
// them from their `-help`, `HideStdArgs(fHideLogfile | fHideConffile | ...)`, blastn_app.cpp:92)
// ```c++
// void CArgDescriptions::AddDefaultFileArguments(const string& default_config)
// {
//     if (!Exist(s_ArgLogFile + 1) ) {
//         AddOptionalKey
//             (s_ArgLogFile+1, "File_Name",
//                 "File to which the program log should be redirected",
//                 CArgDescriptions::eOutputFile);
//     }
//     if (!Exist(s_ArgCfgFile + 1) ) {
//         if (default_config.empty()) {
//             AddOptionalKey
//                 (s_ArgCfgFile + 1, "File_Name",
//                     "Program's configuration (registry) data file",
//                     CArgDescriptions::eInputFile);
// ```
/// `command` with NCBI's two toolkit file arguments, `-logfile` and `-conffile`, in each NCBI
/// application that runs them (the BLAST programs and makeblastdb), hidden as NCBI hides them.
/// Their values are read by the application start (`ncbi_environment::log_file` and
/// `ApplicationStart::of_process`), before the arguments: the argument reader only checks them
/// (one of each, as NCBI's `Argument with this name is defined already`).
fn with_toolkit_file_args(mut command: clap::Command) -> clap::Command {
    let names: Vec<String> = command
        .get_subcommands()
        .map(|sub| sub.get_name().to_string())
        .filter(|name| TOOLKIT_PROGRAMS.contains(&name.as_str()) || name == "makeblastdb")
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
// ```c++
// static const char* s_AutoHelpShowAll  = "help-full";
// static const char* s_AutoHelpXml  = "xmlhelp";
// static const char* s_ExtraName    = "....";
//
// const char* s_ArgLogFile         = "-logfile";
// const char* s_ArgCfgFile         = "-conffile";
// const char* s_ArgVersion         = "-version";
// const char* s_ArgFullVersion     = "-version-full";
// const char* s_ArgFullVersionXml  = "-version-full-xml";
// const char* s_ArgFullVersionJson = "-version-full-json";
// const char* s_ArgDryRun          = "-dryrun";
// ```
// The standard options of an NCBI C++ Toolkit application, which NCBI's BLAST+ programs
// accept without listing them in -help. `-logfile` and `-conffile` run in the BLAST programs and
// makeblastdb (line RG item E1), `-dryrun` in them and blastdbcmd (item E1b), the three in
// rpsbproc (item E1c, its own `DryRun`: `rpsbproc::dry_run`), `-help-full` in every program (item
// E2, `try_parse_from`), the full versions in every program (item E3, `ncbi_version_request`),
// `-xmlhelp` in every program on the command line (item E2x, `xml_help_program`).

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:44-50
// ```c++
// CTblastnAppArgs::CTblastnAppArgs()
// {
//     CRef<IBlastCmdLineArgs> arg;
//     static const string kProgram("tblastn");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                   "Protein Query-Translated Subject BLAST"));
// ```
// The options of NCBI tblastn 2.17.0+ (-help) that LOSAT's TBLASTN does not implement
// (AGENTS.md rule 2: explicit unsupported errors).
const UNPORTED_TBLASTN_ARGS: &[&str] = &[
    "best_hit_overhang",
    "best_hit_score_edge",
    "db_hard_mask",
    "db_soft_mask",
    "entrez_query",
    "gilist",
    "html",
    "import_search_strategy",
    "negative_gilist",
    "remote",
    "subject_besthit",
    "use_sw_tback",
];

fn is_unported_tblastn_arg(name: &str) -> bool {
    UNPORTED_TBLASTN_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastx_args.cpp:44-50
// ```c++
// CTblastxAppArgs::CTblastxAppArgs()
// {
//     CRef<IBlastCmdLineArgs> arg;
//     static const string kProgram("tblastx");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                   "Translated Query-Translated Subject BLAST"));
//     const bool kQueryIsProtein = false;
// ```
// The options of NCBI tblastx 2.17.0+ (-help) that LOSAT's TBLASTX does not implement
// (AGENTS.md rule 2: explicit unsupported errors).
const UNPORTED_TBLASTX_ARGS: &[&str] = &[
    "best_hit_overhang",
    "best_hit_score_edge",
    "db_hard_mask",
    "db_soft_mask",
    "entrez_query",
    "gilist",
    "html",
    "import_search_strategy",
    "max_intron_length",
    "negative_gilist",
    "remote",
    "strand",
    "subject_besthit",
    "sum_stats",
];

fn is_unported_tblastx_arg(name: &str) -> bool {
    UNPORTED_TBLASTX_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-66
// ```c++
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
//     m_BlastDbArgs->SetDatabaseMaskingSupport(true);
//     m_BlastDbArgs->SetIPGFilteringSupport(true);
//     arg.Reset(m_BlastDbArgs);
//     m_Args.push_back(arg);
// ```
// Named capabilities exist in NCBI but are explicitly outside the BLASTX product scope.
const UNPORTED_BLASTX_ARGS: &[&str] = &[
    "html",
    "gilist",
    "negative_gilist",
    "entrez_query",
    "db_soft_mask",
    "db_hard_mask",
    "ipglist",
    "negative_ipglist",
    "best_hit_overhang",
    "best_hit_score_edge",
    "import_search_strategy",
    "remote",
    "use_sw_tback",
];

fn is_unported_blastx_arg(name: &str) -> bool {
    UNPORTED_BLASTX_ARGS.contains(&name)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
// ```c++
//     catch (const CArgException& e) {                                        \
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-230
// ```c++
//     catch (const blast::CBlastException& e) {                               \
//         const string& msg = e.GetMsg();                                     \
//         if (e.GetErrCode() == CBlastException::eInvalidOptions) {           \
//             LOG_POST(Error << "BLAST options error: " << e.GetMsg());       \
//             exit_code = BLAST_INPUT_ERROR;                                  \
//         } else if ((NStr::Find(msg, "Out of memory") != NPOS) ||            \
//             (NStr::Find(msg, "Failed to allocate") != NPOS)) {              \
//             LOG_POST(Error << "BLAST ran out of memory: " << e.GetMsg());   \
//             exit_code = BLAST_OUT_OF_MEMORY;                                \
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
//         }                                                                   \
//     }                                                                       \
//     catch (const blast::CBlastSystemException& e) {                         \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
// ```c++
//     }                                                                       \
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
// ```
#[derive(Debug)]
pub struct NativeError {
    pub exit: i32,
    pub message: String,
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:178-180
// ```c++
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
impl fmt::Display for NativeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-184
// ```c++
//     catch (const CArgException& e) {                                        \
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
impl std::error::Error for NativeError {}

impl NativeError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    /// `CATCH_ALL`'s report of an exception of NCBI's C++ Toolkit (`CException`) that ends a
    /// BLAST program: `Error: `, the exception's `what()` (which ends with a line feed) and the
    /// post's line feed, exit status 255 (`BLAST_UNKNOWN_ERROR`). The one producer of this
    /// report: an exception of NCBI's that LOSAT raises where NCBI throws it becomes this error.
    pub(crate) fn exception(e: &crate::ncbi_exception::NcbiException) -> NativeError {
        NativeError {
            exit: 255,
            message: format!("Error: {}\n", e.what()),
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:118-119
// ```c
// {
//     m_Outfile.exceptions(NcbiBadbit);
// ```
/// The report's output stream (the `-out` file or standard output), which records whether
/// a write to it failed: NCBI's outfmt 0 stream throws there, which the application
/// reports as "BLAST failed to write output" (blast_app_util.hpp:252-255).
pub struct ReportStream {
    pub inner: Box<dyn std::io::Write + Send>,
    pub failed: bool,
}

impl ReportStream {
    /// The command line's report stream: the `-out` file `out_file` (`None`: standard output),
    /// behind NCBI's buffer for the organism report (`tax_report`, `report_destination`).
    pub fn of_cli(out_file: Option<std::io::BufWriter<std::fs::File>>, tax_report: bool) -> Self {
        ReportStream {
            inner: report_destination(
                out_file.map(|file| Box::new(file) as Box<dyn std::io::Write + Send>),
                tax_report,
            ),
            failed: false,
        }
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

// NCBI reference (glibc 2.17 to 2.39): libio/filedoalloc.c `_IO_file_doallocate`
// ```c
//   size = BUFSIZ;
//   if (fp->_fileno >= 0 && __builtin_expect (_IO_SYSSTAT (fp, &st), 0) >= 0)
//     {
//       if (S_ISCHR (st.st_mode))
// 	{
// 	  /* Possibly a tty.  */
// 	  if (
// #ifdef DEV_TTY_P
// 	      DEV_TTY_P (&st) ||
// #endif
// 	      local_isatty (fp->_fileno))
// 	    fp->_flags |= _IO_LINE_BUF;
// 	}
// #if defined _STATBUF_ST_BLKSIZE
//       if (st.st_blksize > 0 && st.st_blksize < BUFSIZ)
// 	size = st.st_blksize;
// #endif
//     }
// ```
/// The size of glibc's buffer of standard output (`cout` writes into it, synchronized with C's
/// `stdout`): `st_blksize` of the standard output where it is below `BUFSIZ`, else `BUFSIZ` (4096
/// for a pipe and for most files, 1024 for a terminal); 4096 where LOSAT cannot ask.
pub fn stdio_buffer_size() -> usize {
    #[cfg(unix)]
    let block = rustix::fs::fstat(std::io::stdout())
        .ok()
        .and_then(|stat| usize::try_from(stat.st_blksize).ok())
        .filter(|&size| size > 0 && size < NCBI_BUFSIZ)
        .unwrap_or(NCBI_BUFSIZ);
    #[cfg(not(unix))]
    let block = 4096;
    block
}

/// glibc's and libstdc++'s `BUFSIZ`.
pub const NCBI_BUFSIZ: usize = 8192;

/// NCBI's buffer of a report stream, which decides when the bytes inserted into the stream reach
/// the file or device (`NcbiBuffer`). The buffers are those of NCBI's builds of BLAST+ and
/// rpsbproc, which are Linux builds (glibc, GCC 13's libstdc++); LOSAT keeps their counts on every
/// platform.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NcbiStdio {
    /// `cout`, synchronized with C's `stdout` (NCBI does not call `sync_with_stdio(false)`):
    /// each insertion is an `fwrite` into glibc's buffer of `block` bytes (`stdio_buffer_size`),
    /// which writes at each newline too where `lines` is set (a terminal: `_IO_LINE_BUF`).
    Cout { block: usize, lines: bool },
    /// An `ofstream` (`-out`, rpsbproc's `-o`): libstdc++'s `basic_filebuf` with a buffer of
    /// `BUFSIZ` bytes, of which an insertion that does not fit in what is left is written
    /// directly with it.
    File,
}

impl NcbiStdio {
    /// The stream of a report on standard output: glibc's buffer of it, line-buffered on a
    /// terminal (`_IO_file_doallocate`, quoted at `stdio_buffer_size`: `local_isatty`, as
    /// `is_terminal`).
    pub fn cout() -> Self {
        use std::io::IsTerminal;
        NcbiStdio::Cout {
            block: stdio_buffer_size(),
            lines: std::io::stdout().is_terminal(),
        }
    }
}

/// What NCBI's buffer of a report stream (`NcbiStdio`) holds of the bytes inserted into it: the
/// last `held` bytes; the others are written to the file or device. One model for the organism
/// report's streams (`NcbiReportBuffer`) and rpsbproc's report (`rpsbproc::Report`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NcbiBuffer {
    stream: NcbiStdio,
    held: usize,
    /// Whether glibc has made `stdout`'s buffer (at the first insertion of a byte).
    allocated: bool,
}

impl NcbiBuffer {
    /// The empty buffer of `stream`.
    pub fn new(stream: NcbiStdio) -> Self {
        NcbiBuffer {
            stream,
            held: 0,
            allocated: false,
        }
    }

    /// The buffer's stream.
    pub fn stream(&self) -> NcbiStdio {
        self.stream
    }

    /// How many of the bytes inserted the buffer holds (the last ones).
    pub fn held(&self) -> usize {
        self.held
    }

    /// An insertion of `n` bytes (`sputn`: a string or a number), the last bytes of `stream`: the
    /// stream's bytes up to the insertion's end, at least the last `NCBI_BUFSIZ` of them or all
    /// since the buffer last wrote (what a line-buffered `cout` holds is read for newlines).
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_xsputn`
    // ```c
    //   /* First figure out how much space is available in the buffer. */
    //   ...
    //   else if (f->_IO_write_end > f->_IO_write_ptr)
    //     count = f->_IO_write_end - f->_IO_write_ptr; /* Space available. */
    //   ...
    //   if (to_do + must_flush > 0)
    //     {
    //       size_t block_size, do_write;
    //       /* Next flush the (full) buffer. */
    //       if (_IO_OVERFLOW (f, EOF) == EOF)
    //     ...
    //       block_size = f->_IO_buf_end - f->_IO_buf_base;
    //       do_write = to_do - (block_size >= 128 ? to_do % block_size : 0);
    // ```
    // Before the first insertion `stdout` has no buffer (no space available): `_IO_OVERFLOW`
    // makes it, and the first insertion keeps `n % block` (glibc 2.39, `T/t2/fix2/
    // glibc_first_write.c` on a file of `st_blksize` 4096: an `fwrite` of 4,096 bytes is written at
    // once, 4,095 and then 1 are held, and 10 more write 4,096).
    // NCBI reference (GCC 13 libstdc++): include/c++/13/bits/fstream.tcc:766-778
    // (`basic_filebuf::xsputn`)
    // ```c++
    // 	  streamsize __bufavail = this->epptr() - this->pptr();
    //
    // 	  // Don't mistake 'uncommitted' mode buffered with unbuffered.
    // 	  if (!_M_writing && _M_buf_size > 1)
    // 	    __bufavail = _M_buf_size - 1;
    //
    // 	  if (__n >= __bufavail)
    // 	    {
    // 	      const streamsize __buffill = this->pptr() - this->pbase();
    // 	      const char* __buf = reinterpret_cast<const char*>(this->pbase());
    // 	      __ret = _M_file.xsputn_2(__buf, __buffill,
    // 				       reinterpret_cast<const char*>(__s),
    // 				       __n);
    // ```
    // A line-buffered `stdout` (a terminal, `_IO_LINE_BUF` from `_IO_file_doallocate`): the same
    // function's branch for it writes, of an insertion that fits, up to its last newline, and the
    // rest of one that does not fit is put byte by byte (`_IO_default_xsputn`), where
    // `_IO_new_file_overflow` writes at each newline; so the buffer keeps what follows the last
    // newline of what it would keep (glibc 2.39, `T/t2/fix3/glibc_tty_lines.c` on a pty of
    // `st_blksize` 1024: 100 then 50 ending in a newline write 150; 1000 then 30 with a newline at
    // 27 write 1,028; 500 then 2,000 with newlines at 100 and 1,900 write 2,401).
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_xsputn`
    // ```c
    //   if ((f->_flags & _IO_LINE_BUF) && (f->_flags & _IO_CURRENTLY_PUTTING))
    //     {
    //       count = f->_IO_buf_end - f->_IO_write_ptr;
    //       if (count >= n)
    //         {
    //           const char *p;
    //           for (p = s + n; p > s; )
    //             {
    //               if (*--p == '\n')
    //                 {
    //                   count = p - s + 1;
    //                   must_flush = 1;
    //                   break;
    //                 }
    //             }
    //         }
    //     }
    // ```
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_overflow`
    // ```c
    //   if ((f->_flags & _IO_UNBUFFERED)
    //       || ((f->_flags & _IO_LINE_BUF) && ch == '\n'))
    //     if (_IO_do_write (f, f->_IO_write_base,
    // 		      f->_IO_write_ptr - f->_IO_write_base) == EOF)
    // ```
    pub fn sputn(&mut self, n: usize, stream: &[u8]) {
        if n == 0 {
            return;
        }
        self.held = match self.stream {
            NcbiStdio::Cout { block, lines } => {
                let room = if self.allocated { block - self.held } else { 0 };
                let held = if n <= room {
                    self.held + n
                } else if block >= 128 {
                    (self.held + n) % block
                } else {
                    // `do_write` is the whole rest; the buffer keeps nothing.
                    0
                };
                if lines {
                    let kept = &stream[stream.len() - held..];
                    kept.iter()
                        .rposition(|&byte| byte == b'\n')
                        .map_or(held, |at| held - at - 1)
                } else {
                    held
                }
            }
            NcbiStdio::File => {
                if n >= NCBI_BUFSIZ - 1 - self.held {
                    0
                } else {
                    self.held + n
                }
            }
        };
        self.allocated = true;
    }

    /// An insertion of one `char` (`put`, `sputc`), the last byte of `stream` (as `sputn`'s): the
    /// buffer takes it while it has room; a full buffer is written with it
    /// (`basic_filebuf::overflow`). NCBI rpsbproc 0.5.1, both builds, with 8,190 bytes held
    /// before a tab of the `SESSION` line (line RG E1c review 2, `s_same_mat8165`,
    /// `s_same_ev8156`). On `cout` the `char` is glibc's `putc` (`stdio_sync_filebuf::overflow`),
    /// which a full or unallocated buffer passes to `__overflow`: the buffer is made or written,
    /// then it keeps the byte (an `fwrite` of one byte keeps none where the buffer is under 128
    /// bytes, `sputn`), and a line-buffered `stdout` writes it with what it holds at a newline.
    // NCBI reference (GCC 13 libstdc++): include/c++/13/ostream:570-576 (a `char` is put, not
    // written with `sputn`)
    // ```c++
    //   template<typename _Traits>
    //     inline basic_ostream<char, _Traits>&
    //     operator<<(basic_ostream<char, _Traits>& __out, char __c)
    //     {
    //       if (__out.width() != 0)
    // 	return __ostream_insert(__out, &__c, 1);
    //       __out.put(__c);
    //       return __out;
    //     }
    // ```
    // NCBI reference (GCC 13 libstdc++): include/ext/stdio_sync_filebuf.h `overflow`
    // ```c++
    //       virtual int_type
    //       overflow(int_type __c = traits_type::eof())
    //       {
    // 	int_type __ret;
    // 	if (traits_type::eq_int_type(__c, traits_type::eof()))
    // 	  {
    // 	    if (std::fflush(_M_file))
    // 	      __ret = traits_type::eof();
    // 	    else
    // 	      __ret = traits_type::not_eof(__c);
    // 	  }
    // 	else
    // 	  __ret = std::putc(__c, _M_file);
    // 	return __ret;
    //       }
    // ```
    // NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_overflow` (`putc` calls it
    // when `_IO_write_ptr >= _IO_write_end`: a full or unallocated buffer, and every `char` of a
    // line-buffered stream, whose `_IO_write_end` stays at `_IO_write_base`)
    // ```c
    //   if ((f->_flags & _IO_CURRENTLY_PUTTING) == 0 || f->_IO_write_base == NULL)
    //     {
    //       /* Allocate a buffer if needed. */
    //       if (f->_IO_write_base == NULL)
    // 	{
    // 	  _IO_doallocbuf (f);
    //     ...
    //       if (f->_mode <= 0 && f->_flags & (_IO_LINE_BUF | _IO_UNBUFFERED))
    // 	f->_IO_write_end = f->_IO_write_ptr;
    //     }
    //   ...
    //   if (f->_IO_write_ptr == f->_IO_buf_end ) /* Buffer is really full */
    //     if (_IO_do_flush (f) == EOF)
    //       return EOF;
    //   *f->_IO_write_ptr++ = ch;
    //   if ((f->_flags & _IO_UNBUFFERED)
    //       || ((f->_flags & _IO_LINE_BUF) && ch == '\n'))
    //     if (_IO_do_write (f, f->_IO_write_base,
    // 		      f->_IO_write_ptr - f->_IO_write_base) == EOF)
    //       return EOF;
    // ```
    pub fn sputc(&mut self, stream: &[u8]) {
        match self.stream {
            NcbiStdio::Cout { block, lines } => {
                self.held = if lines && stream.last() == Some(&b'\n') {
                    0
                } else if !self.allocated || self.held == block {
                    1
                } else {
                    self.held + 1
                };
                self.allocated = true;
            }
            NcbiStdio::File => {
                self.held = if self.held < NCBI_BUFSIZ - 1 {
                    self.held + 1
                } else {
                    0
                };
            }
        }
    }

    /// The stream flushed (`flush`, `endl`, a write to `cerr`, which is tied to `cout`): the
    /// buffer writes what it holds.
    pub fn flush(&mut self) {
        self.held = 0;
    }
}

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
/// `NcbiOstream` that the run opened, with nothing between).
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
// ```c++
// CBlastFormat::~CBlastFormat()
// {
//     try {
//         m_Outfile.exceptions(m_OrigExceptionMask);
//     } catch (...) {/*ignore exceptions*/}
//     m_Outfile.flush();
// }
// ```
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
    /// standard error (`tie_standard_output`), writing to `inner` (the process's standard output,
    /// or a hold of it, `report_destination_through`).
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

/// `out << value` in NCBI's report: one insertion (a string or a number) of `bytes`, one write
/// call, which the organism report's stream takes as one (`NcbiReportBuffer`); other writers just
/// take the bytes. A string that LOSAT streams is `insert_streamed`.
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

    /// Writes the piece where it may end (`insert_streamed`) once it is `PIECE` bytes long.
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
// ```c
// CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
//     : m_FormatType(format_type), ..., m_Outfile(outfile),
// ```
/// The stream of a report as NCBI's formatter writes to it (`CBlastFormat::m_Outfile`): the writer
/// where the run's report goes, as the run opened it (a format's sink, `OutputSink::ostream`: the
/// command line's `ReportStream`, whose organism report is NCBI's buffer, `report_destination`, or
/// a library caller's writer, which has no stream of NCBI's; blastx's command-line destination),
/// with nothing of a program's before it. Each write is one insertion of NCBI's (`out << value`,
/// `insert`, `insert_streamed`), passed on as one write. The organism report (`tax_report`) writes
/// only to this type, so no buffer or hold of a program's can come between its insertions and
/// NCBI's buffer (`NcbiReportBuffer`), which counts each write as one insertion.
pub struct NcbiOstream<'a>(Box<dyn std::io::Write + 'a>);

impl<'a> NcbiOstream<'a> {
    /// The stream `stream` where the run's report goes, as the run opened it.
    pub fn new(stream: Box<dyn std::io::Write + 'a>) -> Self {
        NcbiOstream(stream)
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
/// (`tie_standard_output`): the process's (`process`), or a stream of the run's standard error
/// (`through`: blastx's hold of a run whose later batch may be refused, `archive::ArchiveHold`),
/// so that what standard output writes before a text is written before it where the text is
/// made. A whole message is written with one `write_all` of the stream, which `std::io::Stderr`
/// writes under one lock of standard error, so that no other thread's message comes inside it.
/// An empty write is no write of NCBI's to `cerr` (a query without warnings, whose `ERR_POST`
/// NCBI does not make): it writes nothing of the report.
pub struct CliStderr<W: std::io::Write = std::io::Stderr>(W);

impl CliStderr {
    /// The process's standard error.
    pub fn process() -> Self {
        CliStderr(std::io::stderr())
    }
}

impl<W: std::io::Write> CliStderr<W> {
    /// The run's standard error `stream`.
    pub fn through(stream: W) -> Self {
        CliStderr(stream)
    }
}

impl<W: std::io::Write> std::io::Write for CliStderr<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        tie_standard_output()?;
        self.0.write(buf)
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        if buf.is_empty() {
            return Ok(());
        }
        tie_standard_output()?;
        self.0.write_all(buf)
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        self.write_all(args.to_string().as_bytes())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

/// The report stream's writer of a run: the `-out` file `file` (`None`: standard output), behind
/// NCBI's buffer for the organism report (`tax_report`, `NcbiReportBuffer`).
pub fn report_destination(
    file: Option<Box<dyn std::io::Write + Send>>,
    tax_report: bool,
) -> Box<dyn std::io::Write + Send> {
    report_destination_through(file, tax_report, |device| device)
}

/// `report_destination` with the file or standard output behind `through` (blastx's hold of a run
/// whose later batch may be refused, `archive::ArchiveHold`): NCBI's buffer comes before it, so
/// that what it holds is what the buffer writes, where the buffer writes it.
pub fn report_destination_through(
    file: Option<Box<dyn std::io::Write + Send>>,
    tax_report: bool,
    through: impl FnOnce(Box<dyn std::io::Write + Send>) -> Box<dyn std::io::Write + Send>,
) -> Box<dyn std::io::Write + Send> {
    match (file, tax_report) {
        (Some(file), true) => Box::new(NcbiReportBuffer::file(through(file))),
        (Some(file), false) => through(file),
        (None, true) => Box::new(NcbiReportBuffer::standard_output(through(
            report_standard_output(),
        ))),
        (None, false) => through(report_standard_output()),
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
// ```c++
// string s_ArgExptMsg(const string& name, const string& what, const string& attr)
// {
//     return string("Argument \"") + (name.empty() ? s_ExtraName : name) +
//         "\". " + what + (attr.empty() ? attr : ":  `" + attr + "'");
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:615-619
// ```c++
// void CArg_Ios::x_Open(CArgValue::TFileFlags /*flags*/) const
// {
//     if ( !m_Ios ) {
//         NCBI_THROW(CArgException,eNoFile, s_ArgExptMsg(GetName(),
//             "File is not accessible",AsString()));
// ```
pub fn inaccessible(name: &str, path: &std::path::Path) -> anyhow::Error {
    let value = if path.as_os_str().is_empty() {
        String::new()
    } else {
        format!(":  `{}'", path.display())
    };
    NativeError {
        exit: 1,
        message: format!(
            "Command line argument error: Argument \"{name}\". File is not accessible{value}\n"
        ),
    }
    .into()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
// ```c++
//     catch (const CArgException& e) {                                        \
//         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// ```c++
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:251-255
// ```c++
//     }                                                                       \
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163
// ```c++
// CBlastTabularInfo::~CBlastTabularInfo()
// {
//     m_Ostream.flush();
// }
// ```
/// Writes the text of the run's `NativeError` to standard error and exits with its code; a text
/// that the run already wrote (`post_native_error`) is not written again.
pub fn exit_on_native_error(error: &anyhow::Error) {
    if let Some(e) = error.downcast_ref::<NativeError>() {
        if !posted(&e.message) {
            eprint!("{}", e.message);
        }
        if e.exit < 0 {
            std::process::abort();
        }
        std::process::exit(e.exit);
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
static POSTED_TEXT: std::sync::Mutex<Option<String>> = std::sync::Mutex::new(None);

/// Whether the run already wrote `text` (`post_native_error`).
fn posted(text: &str) -> bool {
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
// ```c++
//     } CATCH_ALL(status)
//
//     if(!m_Bah.GetMessages().empty()) {
//     	const CArgs & a = GetArgs();
//     	PrintErrorArchive(a, m_Bah.GetMessages());
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1160
// ```c++
// void CBlastAppDiagHandler::Post(const SDiagMessage & mess)
// {
// 	if(m_handler != NULL) {
// 		m_handler->Post(mess);
// 	}
// 	if(m_save) {
// ```
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
        let _ = std::io::Write::write_all(&mut std::io::stderr(), native.message.as_bytes());
        *POSTED_TEXT
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(native.message.clone());
    }
}

#[cfg(test)]
mod readme_refusals_tests;

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    use std::os::unix::ffi::OsStringExt;

    /// A device that takes nothing (`/dev/full`), counting the bytes offered to it.
    struct Full(std::sync::Arc<std::sync::atomic::AtomicUsize>);

    impl std::io::Write for Full {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0
                .fetch_add(buf.len(), std::sync::atomic::Ordering::Relaxed);
            Err(std::io::Error::other("full"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn the_organism_reports_file_buffer_writes_at_its_8191st_byte() {
        // probe `T/t2/fix1` `f1.d5.*` and `stdout_probe.sh`: NCBI exits 0 for a report that its
        // buffer holds to the end, 6 where the buffer writes.
        use std::io::Write;
        let offered = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let mut buffer = NcbiReportBuffer::file(Box::new(Full(offered.clone())));
        buffer.write_all(&[b'x'; 8190]).unwrap();
        buffer.flush().unwrap();
        assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 0);
        assert!(buffer.write_all(b"x").is_err());
        assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 8191);
        // What a buffer holds at the end is written, a failure unreported.
        let quiet = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let mut buffer = NcbiReportBuffer::file(Box::new(Full(quiet.clone())));
        buffer.write_all(b"held").unwrap();
        drop(buffer);
        assert_eq!(quiet.load(std::sync::atomic::Ordering::Relaxed), 4);
    }

    /// A device that takes everything, keeping the size of each write.
    struct Writes(std::sync::Arc<std::sync::Mutex<Vec<usize>>>);

    impl std::io::Write for Writes {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            lock(&self.0).push(buf.len());
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn the_organism_reports_standard_output_writes_where_glibc_writes() {
        // glibc 2.39 `_IO_new_file_xsputn` on a file of `st_blksize` 4096
        // (`T/t2/fix2/glibc_first_write.log`): the writes end at multiples of the block and the
        // rest stays held; the first `fwrite` has no buffer yet; a flush empties the buffer.
        use std::io::Write;
        let run = |pieces: &[usize]| {
            let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
            let mut after = Vec::new();
            let mut buffer = NcbiReportBuffer::new(
                Box::new(Writes(writes.clone())),
                NcbiBuffer::new(NcbiStdio::Cout {
                    block: 4096,
                    lines: false,
                }),
            );
            for &piece in pieces {
                match piece {
                    0 => lock(&buffer.0).write_held().unwrap(),
                    n => buffer.write_all(&vec![b'x'; n]).unwrap(),
                }
                after.push(lock(&writes).iter().sum::<usize>());
            }
            drop(buffer);
            let total = lock(&writes).iter().sum::<usize>();
            (after, total)
        };
        // `4000 200 4000`: 4,096 written at the second, 8,192 at the third, 8 held to the end.
        assert_eq!(run(&[4000, 200, 4000]), (vec![0, 4096, 8192], 8200));
        assert_eq!(run(&[4096, 10]), (vec![4096, 4096], 4106));
        assert_eq!(run(&[4095, 1, 10]), (vec![0, 0, 4096], 4106));
        assert_eq!(run(&[8192, 1]), (vec![8192, 8192], 8193));
        assert_eq!(run(&[100, 4096]), (vec![0, 4096], 4196));
        // `100 f 4096 1` (0: the flush of a write to `cerr`, which is tied to `cout`).
        assert_eq!(run(&[100, 0, 4096, 1]), (vec![0, 100, 100, 4196], 4197));
        assert_eq!(run(&[8000, 300]), (vec![4096, 8192], 8300));
    }

    /// The bytes that `buffer` has written to its device after each step of `steps`
    /// (`Ok(bytes)`: one insertion; `Err(())`: a flush, as before a write to `cerr`), and in all.
    fn written_after(
        stream: NcbiStdio,
        steps: &[Result<Vec<u8>, ()>],
    ) -> (Vec<usize>, usize, usize) {
        use std::io::Write;
        let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut after = Vec::new();
        let mut buffer =
            NcbiReportBuffer::new(Box::new(Writes(writes.clone())), NcbiBuffer::new(stream));
        for step in steps {
            match step {
                Ok(bytes) => buffer.write_all(bytes).unwrap(),
                Err(()) => lock(&buffer.0).write_held().unwrap(),
            }
            after.push(lock(&writes).iter().sum::<usize>());
        }
        drop(buffer);
        let writes = lock(&writes);
        (after, writes.iter().sum(), writes.len())
    }

    #[test]
    fn a_terminals_standard_output_writes_at_its_newlines_as_glibc_does() {
        // glibc 2.39 on a pty of `st_blksize` 1024 (`T/t2/fix3/glibc_tty_lines.log`): `N` is an
        // `fwrite` of N bytes, `N:K` with a newline at K, `f` a flush.
        let tty = NcbiStdio::Cout {
            block: 1024,
            lines: true,
        };
        let run = |steps: &str| {
            let steps = steps
                .split(' ')
                .map(|step| match step {
                    "f" => Err(()),
                    _ => {
                        let mut parts = step.split(':').map(|n| n.parse::<usize>().unwrap());
                        let mut bytes = vec![b'x'; parts.next().unwrap()];
                        parts.for_each(|at| bytes[at] = b'\n');
                        Ok(bytes)
                    }
                })
                .collect::<Vec<_>>();
            written_after(tty, &steps).0
        };
        assert_eq!(run("100"), [0]);
        assert_eq!(run("100 50:49"), [0, 150]);
        assert_eq!(run("100 50:10"), [0, 111]);
        assert_eq!(run("3000"), [2048]);
        assert_eq!(run("3000:2500"), [2501]);
        assert_eq!(run("1000 30"), [0, 1024]);
        assert_eq!(run("1000 24 1"), [0, 0, 1024]);
        assert_eq!(run("1000 30:5"), [0, 1024]);
        assert_eq!(run("1000 30:27"), [0, 1028]);
        assert_eq!(run("500 2000:100"), [0, 2048]);
        assert_eq!(run("500 2000:100:1900"), [0, 2401]);
        assert_eq!(run("f 100 f 50"), [0, 0, 100, 100]);
    }

    #[test]
    fn a_char_on_standard_output_is_kept_as_glibcs_putc_keeps_it() {
        // glibc 2.39 (`T/t2/fix4/glibc_putc.log`): `stdout` with a buffer of `block` bytes
        // (`setvbuf`, full or line buffering) on a file; `cN` is N `putc` of 'x', `n` a `putc` of a
        // newline, `wN` an `fwrite` of N bytes, `W` an `fwrite` of a newline; the file's size
        // after each step. `putc` keeps its byte after a full or unallocated buffer, where an
        // `fwrite` of one byte into a buffer under 128 bytes keeps none.
        let run = |lines: bool, block: usize, steps: &str| -> Vec<usize> {
            let mut buffer = NcbiBuffer::new(NcbiStdio::Cout { block, lines });
            let mut stream = Vec::new();
            steps
                .split(' ')
                .map(|step| {
                    let count = step[1..].parse::<usize>().unwrap_or(1);
                    match step.as_bytes()[0] {
                        b'c' => (0..count).for_each(|_| {
                            stream.push(b'x');
                            buffer.sputc(&stream);
                        }),
                        b'n' => {
                            stream.push(b'\n');
                            buffer.sputc(&stream);
                        }
                        b'w' => {
                            stream.resize(stream.len() + count, b'x');
                            buffer.sputn(count, &stream);
                        }
                        _ => {
                            stream.push(b'\n');
                            buffer.sputn(1, &stream);
                        }
                    }
                    stream.len() - buffer.held()
                })
                .collect()
        };
        assert_eq!(run(false, 64, "c1"), [0]);
        assert_eq!(run(false, 64, "w1"), [1]);
        assert_eq!(run(false, 64, "w63 c64 c1 w1"), [63, 63, 127, 127]);
        assert_eq!(run(false, 64, "w63 c64 w1"), [63, 63, 128]);
        assert_eq!(run(false, 64, "w64 c1"), [64, 64]);
        assert_eq!(run(false, 64, "w64 w1"), [64, 64]);
        assert_eq!(run(false, 4096, "w4096 c1"), [4096, 4096]);
        assert_eq!(run(false, 4096, "w4095 c1 c1"), [0, 0, 4096]);
        assert_eq!(run(false, 4096, "w4095 w1 w1"), [0, 0, 4096]);
        assert_eq!(run(true, 64, "c1"), [0]);
        assert_eq!(run(true, 64, "n"), [1]);
        assert_eq!(run(true, 64, "W"), [1]);
        assert_eq!(run(true, 64, "c10 n"), [0, 11]);
        assert_eq!(run(true, 64, "w63 c1 c1"), [63, 63, 63]);
        assert_eq!(run(true, 64, "w63 w1 w1"), [63, 63, 63]);
        assert_eq!(run(true, 64, "c64 c1"), [0, 64]);
        assert_eq!(run(true, 64, "c63 n"), [0, 64]);
        assert_eq!(run(true, 64, "c64 n"), [0, 65]);
        assert_eq!(run(true, 1024, "c100 n c5"), [0, 101, 101]);
    }

    #[test]
    fn a_streamed_insertion_is_written_as_the_whole_insertion() {
        // `insert_streamed` writes pieces; every model of NCBI's buffers must write what it writes
        // for the whole insertion as one write, before it and after it. A row of 70,000 bytes
        // ends 99 bytes after a newline: no piece may end there (a terminal would keep only 99).
        let row = [&[b'r'; 69_900][..], b"\n", &[b'r'; 99]].concat();
        let segments: &[(&[u8], u64)] = &[
            (b"\n  Organism Report\nrow\n", 200_000),
            (&row, 3_000),
            (b"|x\ny", 5_000),
            (b"z", 70_000),
            (b"tail\nend", 0),
        ];
        let mut whole = Vec::new();
        for (bytes, spaces) in segments {
            whole.extend_from_slice(bytes);
            whole.resize(whole.len() + *spaces as usize, b' ');
        }
        let streamed = |out: &mut dyn std::io::Write| {
            insert_streamed(out, whole.len() as u64, |insertion| {
                for (bytes, spaces) in segments {
                    insertion.bytes(bytes)?;
                    insertion.spaces(*spaces)?;
                }
                Ok(())
            })
            .unwrap()
        };
        let mut flat = Vec::new();
        streamed(&mut flat);
        assert_eq!(flat, whole);
        for stream in [
            NcbiStdio::File,
            NcbiStdio::Cout {
                block: 4096,
                lines: false,
            },
            NcbiStdio::Cout {
                block: 1024,
                lines: true,
            },
            NcbiStdio::Cout {
                block: 100,
                lines: false,
            },
        ] {
            for before in [0, 1, 1000, 4095, 4096, 8000] {
                for after in [&b"ab"[..], b"ab\ncd", &[b'y'; 9000]] {
                    let lead = vec![b'w'; before];
                    // The whole as one insertion, between the same insertions.
                    let one = written_after(
                        stream,
                        &[Ok(lead.clone()), Ok(whole.clone()), Ok(after.to_vec())],
                    );
                    // The whole streamed.
                    let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
                    let mut buffer = NcbiReportBuffer::new(
                        Box::new(Writes(writes.clone())),
                        NcbiBuffer::new(stream),
                    );
                    let mut pieces = Vec::new();
                    std::io::Write::write_all(&mut buffer, &lead).unwrap();
                    pieces.push(lock(&writes).iter().sum::<usize>());
                    streamed(&mut buffer);
                    pieces.push(lock(&writes).iter().sum::<usize>());
                    std::io::Write::write_all(&mut buffer, after).unwrap();
                    pieces.push(lock(&writes).iter().sum::<usize>());
                    drop(buffer);
                    let total = lock(&writes).iter().sum::<usize>();
                    assert_eq!((pieces, total), (one.0, one.1), "{stream:?} {before}");
                }
            }
        }
        // In pieces: the run of 200,000 spaces ends three, the rest is the last.
        let pieces = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        streamed(&mut Writes(pieces.clone()));
        assert_eq!(*lock(&pieces), [65_559, 65_536, 65_536, 151_405]);
    }

    fn blastn_argv(words: Vec<OsString>) -> Vec<OsString> {
        let mut argv = vec![
            OsString::from("losat"),
            OsString::from("blastn"),
            OsString::from("-query"),
            OsString::from("q.fa"),
            OsString::from("-subject"),
            OsString::from("s.fa"),
        ];
        argv.extend(words);
        argv
    }

    // Without `main`'s `post_errors_in_run` (a library caller of the command line entries such as
    // `api::run_blastn`), the run writes no error text of its own: the text stays in the error.
    // No unit test calls `post_errors_in_run`.
    #[test]
    fn a_run_outside_main_posts_no_error_text() {
        let text = "BLAST query error: Could not construct seq-id from 'gi|abc'\n";
        let error: anyhow::Error = NativeError {
            exit: 1,
            message: text.to_string(),
        }
        .into();
        post_native_error(&error);
        assert!(!posted(text));
    }

    // NCBI c++/src/corelib/ncbiargs.cpp:266-270 keeps the bytes of a string value; LOSAT
    // rejects a value that is not UTF-8 explicitly, in both spellings of an option.
    #[test]
    fn option_values_that_are_not_utf8_are_rejected_explicitly() {
        for (option, value) in [
            ("-query_loc", vec![0xff]),
            ("-subject_loc", b"1-\xff".to_vec()),
            ("-dust", vec![0xfe]),
        ] {
            let separate = blastn_argv(vec![
                OsString::from(option),
                OsString::from_vec(value.clone()),
            ]);
            let mut inline = format!("{option}=").into_bytes();
            inline.extend(&value);
            let joined = blastn_argv(vec![OsString::from_vec(inline)]);
            for argv in [separate, joined] {
                let error = try_parse_from::<Cli, _, _>(argv).unwrap_err();
                let message = render_message(&error);
                assert!(
                    message.contains(&format!("the value of {option} is not UTF-8"))
                        && message.contains("not supported by LOSAT's BLASTN"),
                    "{message}"
                );
            }
        }
    }

    // BLASTX reads both ranges (line SX, SXd1), excludes `-subject_loc` with `-db` as the other
    // programs do, and rejects an option's value that is not UTF-8 with their text (DECISIONS 10
    // (c) of line SX).
    #[test]
    fn blastx_reads_the_ranges_and_rejects_values_that_are_not_utf8() {
        let parsed = try_parse_from::<Cli, _, _>([
            "losat",
            "blastx",
            "-query",
            "q.fa",
            "-subject",
            "s.faa",
            "-query_loc",
            "1-10",
            "-subject_loc",
            "2-5",
        ])
        .unwrap();
        let Commands::Blastx(args) = parsed.command else {
            panic!("a BLASTX command line");
        };
        assert_eq!(args.query_loc.as_deref(), Some("1-10"));
        assert_eq!(args.subject_loc.as_deref(), Some("2-5"));
        let error = try_parse_from::<Cli, _, _>([
            "losat",
            "blastx",
            "-query",
            "q.fa",
            "-db",
            "d",
            "-subject_loc",
            "2-5",
        ])
        .unwrap_err();
        assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
        for (option, value) in [
            ("-query_loc", vec![0xff]),
            ("-subject_loc", b"1-\xff".to_vec()),
            ("-outfmt", vec![0xfe]),
        ] {
            let mut argv = blastn_argv(vec![OsString::from(option), OsString::from_vec(value)]);
            argv[1] = OsString::from("blastx");
            let message = render_message(&try_parse_from::<Cli, _, _>(argv).unwrap_err());
            assert!(
                message.contains(&format!("the value of {option} is not UTF-8"))
                    && message.contains("not supported by LOSAT's BLASTX"),
                "{message}"
            );
        }
    }

    // The rejection comes after the parser, as NCBI reads the value only after its
    // arguments parse: `-help` and the parser's errors come first.
    #[test]
    fn help_and_parser_errors_come_before_a_value_that_is_not_utf8() {
        let mut argv = vec![OsString::from("losat"), OsString::from("blastn")];
        argv.push(OsString::from("-help"));
        argv.push(OsString::from("-query_loc"));
        argv.push(OsString::from_vec(vec![0xff]));
        let error = try_parse_from::<Cli, _, _>(argv).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::DisplayHelp);
        let argv = blastn_argv(vec![
            OsString::from("-query_loc"),
            OsString::from_vec(vec![0xff]),
            OsString::from("-no_such_option"),
        ]);
        let message = render_message(&try_parse_from::<Cli, _, _>(argv).unwrap_err());
        assert!(!message.contains("not UTF-8"), "{message}");
    }

    // `try_parse_from` is generic: a caller's own type can name an option with letters that are
    // not ASCII, here U+FFFD (3 bytes), which a byte that is not UTF-8 also reads as in the
    // word's lossy text. The value of `-name=value` starts after the first `=` of the bytes,
    // not `name.len()` bytes on (which skipped `\xfe\xfd` and read `abc`).
    #[test]
    fn inline_value_starts_after_the_first_equals_sign() {
        #[derive(clap::Parser)]
        enum Own {
            Blastn(OwnArgs),
        }
        #[derive(clap::Args)]
        struct OwnArgs {
            #[arg(long = "q\u{fffd}")]
            q: Option<String>,
        }
        let argv = vec![
            OsString::from("losat"),
            OsString::from("blastn"),
            OsString::from_vec(b"-q\xff=\xfe\xfdabc".to_vec()),
        ];
        let error = try_parse_from::<Own, _, _>(argv)
            .err()
            .expect("a value that is not UTF-8");
        assert_eq!(error.kind(), ErrorKind::InvalidUtf8);
    }

    // A file name keeps its bytes in both spellings (`check_utf8_file_name` reads it).
    #[test]
    fn file_names_keep_their_bytes_in_both_spellings() {
        let name = b"q\xff.fa".to_vec();
        let mut inline = b"-query=".to_vec();
        inline.extend(&name);
        let argv = vec![
            OsString::from("losat"),
            OsString::from("blastn"),
            OsString::from_vec(inline),
            OsString::from("-subject"),
            OsString::from("s.fa"),
        ];
        let Commands::Blastn(args) = try_parse_from::<Cli, _, _>(argv).unwrap().command else {
            panic!("not blastn");
        };
        assert_eq!(args.query.into_os_string(), OsString::from_vec(name));
    }
}

/// The process's arguments as NCBI's `GetCmdlineArgs` joins them for the `@PG` line of a SAM
/// report: every argument (the program name first) followed by one space, as raw bytes where
/// the system keeps them (Unix), else lossily decoded.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:878-885
/// ```c++
/// string
/// GetCmdlineArgs(const CNcbiArguments & a)
/// {
/// 	string cmd = kEmptyStr;
/// 	for(unsigned int i=0; i < a.Size(); i++) {
/// 		cmd += a[i] + " ";
/// 	}
/// 	return cmd;
/// }
/// ```
pub fn command_line_bytes() -> Vec<u8> {
    let mut line = Vec::new();
    for argument in std::env::args_os() {
        #[cfg(unix)]
        {
            use std::os::unix::ffi::OsStrExt;
            line.extend_from_slice(argument.as_bytes());
        }
        #[cfg(not(unix))]
        line.extend_from_slice(argument.to_string_lossy().as_bytes());
        line.push(b' ');
    }
    line
}

#[cfg(test)]
#[path = "cli_v2_tests.rs"]
mod cli_v2_tests;
