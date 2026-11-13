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
pub fn try_parse_from<T, I, S>(argv: I) -> Result<T, clap::Error>
where
    T: Parser + CommandFactory,
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    try_parse_from_command(T::command(), argv)
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
    let mut tokens: Vec<OsString> = argv.into_iter().map(Into::into).collect();
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
            return Err(toolkit_file_word_error(program.get_name(), word));
        }
        if crate::blastinput::ncbi_environment::application_start_of_process().is_none() {
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
        // `-dryrun` goes to the entry (`main_entry`), and the library API refuses `-conffile`
        // and `-dryrun` (NCBI BLAST+ 2.17.0, probes `makeblastdb.log_noval`,
        // `makeblastdb.conf_noval`: rc 0).
        use crate::blastinput::ncbi_environment::{preparse, Preparsed};
        if let Some(word) = ncbi_version_word(&tokens[2..]) {
            return Err(toolkit_file_word_error("makeblastdb", word));
        }
        if let Preparsed::Words {
            words,
            conffile,
            dryrun,
        } = preparse(&tokens[2..])
        {
            let no_start =
                crate::blastinput::ncbi_environment::application_start_of_process().is_none();
            if no_start && (conffile.is_some() || dryrun) {
                // The library's makeblastdb takes `-logfile` (its own argument): the word named is
                // the one refused.
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
// (DI-6) and `-min_raw_gapped_score` (DI-6b); the database masks are not.
const UNPORTED_BLASTN_ARGS: &[&str] = &[
    "best_hit_overhang",
    "best_hit_score_edge",
    "db_hard_mask",
    "db_soft_mask",
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
/// `algorithm::psiblast::in_pssm`), PHI-BLAST (PS8), the search strategy file (item F1b of line
/// RG), and the arguments that blastp's engine does not run either (`is_unported_blastp_arg`;
/// psiblast has no database masks and no `-mt_mode`). `-sum_stats true` is refused with the
/// options (`PsiblastArgs::into_blastp_args`).
const UNPORTED_PSIBLAST_OWN_ARGS: &[&str] = &[
    // Line RG-F1 ported it to BLASTP's search, not to PSIBLAST's (item F1b).
    "export_search_strategy",
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
/// redo, PS11c-3), `-html` (PS15a) and the remote and search strategy options of every program.
const UNPORTED_RPSBLAST_ARGS: &[&str] = &[
    "entrez_query",
    "export_search_strategy",
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
/// rpstblastn's NCBI arguments that LOSAT does not run: `-html` (PS15a) and the remote and search
/// strategy options of every program. The composition-based statistics of a later item are
/// refused once per run (`rpstblastn::search::TranslatedRpsSearch::check_options`).
fn is_unported_rpstblastn_arg(name: &str) -> bool {
    matches!(
        name,
        "entrez_query" | "export_search_strategy" | "html" | "import_search_strategy" | "remote"
    )
}

/// LOSAT's refusal of the NCBI C++ Toolkit option `-name` in `program` (decision D9).
fn toolkit_option_error(program: &str, name: &str) -> clap::Error {
    // The BLAST programs are named in capitals; makeblastdb and rpsbproc as in their other
    // refusals.
    let shown = match program {
        "makeblastdb" | "rpsbproc" => program.to_string(),
        _ => program.to_ascii_uppercase(),
    };
    clap::Error::raw(
        ErrorKind::InvalidValue,
        format!("the NCBI C++ Toolkit option -{name} is not supported by LOSAT's {shown}"),
    )
}

fn unknown_option_error(scope: &clap::Command, name: &str, text: &str) -> clap::Error {
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
    if scope.get_name() == "blastx" && is_unported_blastx_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("unsupported BLASTX option '-{name}': outside the declared local FASTA scope"),
        );
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
    if scope.get_name() == "makeblastdb"
        && crate::writedb::makeblastdb::UNSUPPORTED_OPTIONS.contains(&name)
    {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's makeblastdb"),
        );
    }
    if scope.get_name() == "blastn" && is_unported_blastn_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's BLASTN"),
        );
    }
    if scope.get_name() == "blastp" && is_unported_blastp_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's BLASTP"),
        );
    }
    if scope.get_name() == "psiblast" && is_unported_psiblast_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's PSIBLAST"),
        );
    }
    if scope.get_name() == "rpsblast" && is_unported_rpsblast_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's RPSBLAST"),
        );
    }
    if scope.get_name() == "rpstblastn" && is_unported_rpstblastn_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's RPSTBLASTN"),
        );
    }
    // makeblastdb reads the toolkit's words with its own branch of `try_parse_from`, and refuses
    // the word left as the others do (`-xmlhelp`).
    if let Some(program) = TOOLKIT_PROGRAMS
        .iter()
        .chain(&["makeblastdb"])
        .copied()
        .find(|program| scope.get_name() == *program)
        .filter(|_| is_ncbi_toolkit_arg(name))
    {
        return toolkit_option_error(program, name);
    }
    if scope.get_name() == "tblastx" && is_unported_tblastx_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's TBLASTX"),
        );
    }
    if scope.get_name() == "tblastn" && is_unported_tblastn_arg(name) {
        return clap::Error::raw(
            ErrorKind::InvalidValue,
            format!("the NCBI BLAST+ option -{name} is not supported by LOSAT's TBLASTN"),
        );
    }
    clap::Error::raw(
        ErrorKind::UnknownArgument,
        format!("unknown option or argument '{text}'; use -help for CLI v2 syntax"),
    )
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

/// The NCBI applications whose words before the arguments `try_parse_from` reads
/// (`ncbi_environment::preparse`).
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
/// here: `main` prints the version before parsing (`ncbi_version_request`). The library API and
/// the Web engine, which parse their options here and run no application start, refuse it
/// (`toolkit_file_word_error`).
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

// No NCBI counterpart: NCBI has no library API. The application start of the command line
// (`cli::main_entry`, `ncbi_environment::start_application`) reads `-logfile` and `-conffile`
// before the arguments and prints the versions, and the entry ends a `-dryrun` after them;
// LOSAT's library and the Web engine run none, so they do not take them (as the library takes
// no PSSM file, `api::search`).
/// LOSAT's refusal of the toolkit option `word` (`-logfile`, `-conffile`, `-dryrun`, `-version`
/// or a full version) in `program` when the parser's caller ran no application start: the
/// library API and the Web engine.
fn toolkit_file_word_error(program: &str, word: &str) -> clap::Error {
    clap::Error::raw(
        ErrorKind::InvalidValue,
        toolkit_file_word_refusal(program, word),
    )
}

/// The text of `toolkit_file_word_error`, which blastdbcmd's library run writes too.
pub(crate) fn toolkit_file_word_refusal(program: &str, word: &str) -> String {
    let name = word.trim_start_matches('-');
    let name = name.split_once('=').map_or(name, |(key, _)| key);
    let shown = match program {
        "makeblastdb" => program.to_string(),
        _ => program.to_ascii_uppercase(),
    };
    format!(
        "the NCBI C++ Toolkit option -{name} runs on the command line of LOSAT's {shown} only: the library API and the Web engine run no NCBI application start, which reads it"
    )
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
/// `start_application`), before the arguments: the argument reader only checks them (one of
/// each, as NCBI's `Argument with this name is defined already`).
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
// E2, `try_parse_from`), the full versions in every program (item E3, `ncbi_version_request`).
// The word left is refused by name, in makeblastdb too.
const NCBI_TOOLKIT_ARGS: &[&str] = &["xmlhelp"];

fn is_ncbi_toolkit_arg(name: &str) -> bool {
    NCBI_TOOLKIT_ARGS.contains(&name)
}

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
