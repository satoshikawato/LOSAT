//! The configuration that NCBI's application layer (`CNcbiApplication`) reads before a BLAST
//! program runs: environment variables of the NCBI C++ Toolkit and the registry files
//! (`.ncbirc` and the program's `.ini`). LOSAT does not reproduce the settings that change
//! a program's output or exit status, and rejects them explicitly (plan DW-13); a registry
//! file with only entries that change no output (such as `[BLAST] BLASTDB`) is accepted.

use std::ffi::{OsStr, OsString};
use std::io::Write;
use std::path::{Path, PathBuf};

use super::ncbi_registry::{
    error_reading_post, last_value, read_registry, registry_text, repeated_entry_post,
    running_error_reading_post, running_syntax_error_notice, syntax_error_notice, Layer,
    RegistryEntry, RegistryRead,
};
use crate::ncbi_str::string_to_bool;
use crate::report::diag_log::{Posts, CRITICAL, ERROR, WARNING};

/// Registry entries (section, name; case-insensitive) that change no output of a search with
/// `-subject`: the database and data-loader settings, the data directory, the directory of
/// the WindowMasker files of blastn's `-window_masker_taxid` (`window_masker_path_of`), and
/// `[NCBI] DONT_USE_NCBIRC`, which only decides whether `.ncbirc` is read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1436-1439
const HARMLESS_ENTRIES: &[(&str, &str)] = &[
    ("BLAST", "BLASTDB"),
    ("BLAST", "BLASTMAT"),
    ("BLAST", "DATA_LOADERS"),
    ("BLAST", "BLASTDB_NUCL_DATA_LOADER"),
    ("BLAST", "BLASTDB_PROT_DATA_LOADER"),
    ("BLAST", "IGDATA"),
    ("BLAST", "MAX_SEQID_LENGTH"),
    ("NCBI", "DATA"),
    ("NCBI", "DONT_USE_NCBIRC"),
    ("WINDOW_MASKER", "WINDOW_MASKER_PATH"),
];

/// Environment variables, besides the `DIAG_` and `NCBI_CONFIG__` families, that the NCBI
/// C++ Toolkit reads for every application and that change its error output or exit status.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1031-1044
const REJECTED_VARIABLES: &[&str] = &[
    "ABORT_ON_THROW",
    "DEBUG_STACK_TRACE_LEVEL",
    "EXCEPTION_STACK_TRACE_LEVEL",
    "NCBI_ABORT_ON_COBJECT_THROW",
    "NCBI_ABORT_ON_NULL",
    "LOG_TRUNCATE",
    "LOG_NOCREATE",
    "LOG_PERFLOGGING",
    "THREAD_STACK_SIZE",
    // Not diagnostics: they change how NCBI's `CSeq_id_Mapper` matches Seq-ids, and with it which
    // sequence a mask of `makeblastdb -mask_data` reaches (`seqid::mask_handle`); LOSAT
    // implements the default (both on) and has its own message for them (`rejected_variable`).
    //
    // NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:692-707
    // NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:692-703
    // With `OBJECTS_PACK_TEXTID=0` a GenBank, EMBL or DDBJ accession no longer matches the same
    // accession of another of the three; with `OBJECTS_PACK_GENERAL=0` the order of the ranges of
    // two `general` ids changes (and both change the order of the messages of a PSSM query's
    // annotations, `algorithm::pssm_search::annot_index`). Refused at any value, NCBI's default 1 included: accepting it
    // would mean porting how NCBI reads a bool parameter, for a setting nobody gives. The registry
    // forms (`[OBJECTS] PACK_TEXTID` in a registry file, `NCBI_CONFIG__OBJECTS__PACK_TEXTID`) are
    // not harmless entries, so they are refused like every entry LOSAT does not know.
    "OBJECTS_PACK_TEXTID",
    "OBJECTS_PACK_GENERAL",
];

fn is_harmless_entry(section: &str, name: &str) -> bool {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:209-211
    // A Boolean only decides whether the usage is reported. Only the usage report's registry
    // reads the entry, and another value aborts NCBI there, which `application_start` refuses
    // (`usage_report_entry_refusal`); the application's registry never reads it (NCBI BLAST+
    // 2.17.0 with `[BLAST] BLAST_USAGE_REPORT=maybe` in `blastn.ini`, or in `.ncbirc` with
    // `BLAST_USAGE_REPORT=false`: rc 0, g_reg probes `l1.ini`, `l1.rc_ur0`).
    if section.eq_ignore_ascii_case("BLAST") && name.eq_ignore_ascii_case("BLAST_USAGE_REPORT") {
        return true;
    }
    HARMLESS_ENTRIES
        .iter()
        .any(|(s, n)| section.eq_ignore_ascii_case(s) && name.eq_ignore_ascii_case(n))
}

/// Why LOSAT rejects an environment variable, or `None`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/env_reg.cpp:418-469
/// `kPrefix` is "NCBI_CONFIG_" (env_reg.cpp:362): `NCBI_CONFIG__<SECTION>__<NAME>` sets a
/// registry entry, `NCBI_CONFIG_<NAME>__<SECTION>` a special one, and
/// `NCBI_CONFIG_OVERRIDES` names a file of entries (ncbireg.cpp:1585-1602), while
/// `NCBI_CONFIG_PATH` is the search path of the files (`registry_search_path`). The
/// `DIAG_` variables are the parameters of the diagnostics (`DIAG_POST_LEVEL`,
/// `DIAG_OLD_POST_FORMAT`, `DIAG_TEE_TO_STDERR`, ...; ncbidiag.cpp), each of which changes
/// stderr or, with a value its parser rejects, stops NCBI.
fn rejected_variable(name: &OsStr) -> Option<String> {
    let name = name.to_str()?;
    if name == "NCBI_CONFIG_PATH" {
        return None;
    }
    if let Some(entry) = name.strip_prefix("NCBI_CONFIG__") {
        if let Some((section, key)) = entry.split_once("__") {
            if !section.is_empty() && !key.is_empty() && is_harmless_entry(section, key) {
                return None;
            }
        }
        return Some(format!(
            "the environment variable {name}, which sets an entry of NCBI BLAST+'s registry that LOSAT does not know to change no output (it accepts only the entries listed for registry files), is not supported by LOSAT"
        ));
    }
    if name.starts_with("NCBI_CONFIG_") {
        return Some(format!(
            "the environment variable {name}, which sets entries of NCBI BLAST+'s registry, is not supported by LOSAT"
        ));
    }
    if name.starts_with("OBJECTS_PACK_") && REJECTED_VARIABLES.contains(&name) {
        return Some(format!(
            "the environment variable {name} is not supported by LOSAT: NCBI BLAST+ reads it to change how Seq-ids are matched (which sequences a mask of makeblastdb reaches); LOSAT does not read it, whatever its value"
        ));
    }
    if name.starts_with("DIAG_") || REJECTED_VARIABLES.contains(&name) {
        return Some(format!(
            "the environment variable {name}, a parameter of NCBI BLAST+'s diagnostics or error handling (LOSAT rejects the whole family: some members change NCBI's messages or exit status, some only on errors), is not supported by LOSAT"
        ));
    }
    // `BLAST_USAGE_REPORT` and `NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT` that are not Booleans
    // abort NCBI's usage report, which `application_start` refuses where NCBI aborts.
    None
}

/// The directories where NCBI looks for its registry files, in order.
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:331-396
/// `kConfigPathDelim` is ":;" outside Windows (metareg.cpp:54-56). An empty `NCBI_CONFIG_PATH`
/// splits into no directory (`NStr::Split` of an empty string adds nothing,
/// include/corelib/ncbistr_util.hpp:271-272), so it has no empty element and the search path
/// is empty. The home directory is `ncbi_home`'s, the program's directories are
/// `process.program_dirs`' (`program_directories`) for the time the path is made.
///
/// NCBI makes the path once, the first time a registry file is looked for (`CMetaRegistry`'s
/// one instance, metareg.hpp:258-261), and that is when the program's usage report object is
/// made, before `AppMain` sets the program's arguments, unless that object looks for no file
/// (`BLAST_USAGE_REPORT` with a false value, or `.ncbirc` turned off in the environment): then
/// it is when the application loads its configuration, after they are set.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:84-86, 195-209
/// (`CBlastUsageReport`, a member of every BLAST+ application, e.g. blastn_app.cpp)
/// NCBI BLAST+ 2.17.0 (`blastn` run by a link on `PATH`, a malformed `.ncbirc` or
/// `blastn.ini`): one in the link's directory is read only with `BLAST_USAGE_REPORT=false` or
/// `=0`, or (`blastn.ini`) with `NCBI_DONT_USE_NCBIRC` or `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`;
/// one in the directory of the executable in every case.
fn registry_search_path(
    program: &str,
    env: &dyn Fn(&str) -> Option<OsString>,
    process: &ProcessFacts,
) -> Result<Vec<PathBuf>, String> {
    let delimiter = |byte: &u8| *byte == b';' || (!cfg!(windows) && *byte == b':');
    let mut path: Vec<PathBuf> = Vec::new();
    let mut tail: Vec<PathBuf> = Vec::new();
    if let Some(config_path) = env("NCBI_CONFIG_PATH") {
        // The variable's bytes (`NcbiSys_getenv`), as the directories' names: a directory whose
        // name is not UTF-8 is searched as NCBI searches it (probe R7 of line LIB-B).
        let config_path = crate::utils::os_bytes::os_str_bytes(&config_path);
        if config_path.is_empty() {
            return Ok(Vec::new());
        }
        let directory = |part: &&[u8]| crate::utils::os_bytes::path_from_bytes(*part);
        let parts: Vec<&[u8]> = config_path.split(delimiter).collect();
        match parts.iter().position(|part| part.is_empty()) {
            None => return Ok(parts.iter().map(directory).collect()),
            Some(empty) => {
                path.extend(parts[..empty].iter().map(directory));
                tail.extend(
                    parts[empty + 1..]
                        .iter()
                        .filter(|part| !part.is_empty())
                        .map(directory),
                );
            }
        }
    }
    if env("NCBI_DONT_USE_LOCAL_CONFIG").is_none() {
        path.push(PathBuf::from("."));
        let home = ncbi_home(env, process.home_by_uid)?;
        if !home.is_empty() {
            path.push(PathBuf::from(home));
        }
    }
    if let Some(ncbi) = env("NCBI") {
        if !ncbi.is_empty() {
            path.push(PathBuf::from(ncbi));
        }
    }
    if cfg!(windows) {
        if let Some(root) = env("SYSTEMROOT") {
            if !root.is_empty() {
                path.push(PathBuf::from(root));
            }
        }
    } else {
        path.push(PathBuf::from("/etc"));
    }
    // The path is made when the usage report's registry looks for `.ncbirc`, or else later
    // (`usage_report_reads_ncbirc`).
    let arguments_set = !usage_report_reads_ncbirc(program, env);
    path.extend((process.program_dirs)(arguments_set)?);
    path.extend(tail);
    Ok(path)
}

/// What NCBI's registry search path reads of the process besides its environment
/// (`registry_search_path`), given as functions so that a test can set them.
struct ProcessFacts<'a> {
    /// The home directory of the user ID's entry in the user database (`s_GetHomeByUID`):
    /// `None` when the user ID has no entry. Asked only on Unix when `HOME` is not set.
    home_by_uid: &'a dyn Fn() -> Option<OsString>,
    /// The program's directories, in order, before or after the program's arguments are set
    /// (`program_directories`), or why LOSAT cannot tell them.
    program_dirs: &'a dyn Fn(bool) -> Result<Vec<PathBuf>, String>,
    /// The files that the registry reads look at.
    files: &'a dyn RegistryFiles,
}

/// The files that NCBI's registry reads look at (`find_registry`, `conf_file_path`,
/// `read_registry_file`): the file system of the process (`FileSystem`), or the files that LOSAT
/// Web's host registered (`ApplicationStart::with_files`). No NCBI counterpart: NCBI reads its file
/// system.
pub trait RegistryFiles {
    /// Whether `path` is a file (`CFile::Exists`: a file, links followed).
    fn is_file(&self, path: &Path) -> bool;
    /// Whether the file at `path` can be opened for reading.
    fn opens(&self, path: &Path) -> bool {
        self.read(path).is_some()
    }
    /// The bytes of the file at `path`: `None` when it cannot be opened, `Some(Err)` when it was
    /// opened and its read failed.
    fn read(&self, path: &Path) -> Option<std::io::Result<Vec<u8>>>;
}

/// The file system of this process, which the command line and the library read
/// (`ApplicationStart::of_process`).
pub struct FileSystem;

impl RegistryFiles for FileSystem {
    fn is_file(&self, path: &Path) -> bool {
        path.is_file()
    }

    fn opens(&self, path: &Path) -> bool {
        std::fs::File::open(path).is_ok()
    }

    fn read(&self, path: &Path) -> Option<std::io::Result<Vec<u8>>> {
        let mut file = std::fs::File::open(path).ok()?;
        let mut bytes = Vec::new();
        Some(std::io::Read::read_to_end(&mut file, &mut bytes).map(|_| bytes))
    }
}

impl ProcessFacts<'static> {
    /// The facts of this process.
    fn of_this_process() -> Self {
        Self {
            // Rust's `home_dir` on Unix reads `HOME` and, when it is not set, the `pw_dir` of
            // `getpwuid_r(getuid())`; it is asked only when `HOME` is not set (`ncbi_home`).
            home_by_uid: &|| std::env::home_dir().map(PathBuf::into_os_string),
            program_dirs: &program_directories,
            files: &FileSystem,
        }
    }
}

/// `CDir::GetHome()`, the directory that `registry_search_path` adds after ".".
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3624-3657
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3586-3618
/// A variable that is set is the home directory, an empty one too (the search path skips an
/// empty home, and `AddTrailingPathSeparator` changes no file name that LOSAT joins to it).
/// The login name's entry (`getpwnam`, `getlogin`) is not reachable from Rust's standard
/// library: when `HOME` is not set and the user ID has no entry, LOSAT refuses to run. A
/// target that is neither Unix nor Windows (WASI) has no user database: there the home is
/// `HOME` or none.
fn ncbi_home(
    env: &dyn Fn(&str) -> Option<OsString>,
    home_by_uid: &dyn Fn() -> Option<OsString>,
) -> Result<OsString, String> {
    if cfg!(windows) {
        return Ok(env("APPDATA")
            .or_else(|| env("USERPROFILE"))
            .unwrap_or_default());
    }
    if let Some(home) = env("HOME") {
        return Ok(home);
    }
    if !cfg!(unix) {
        return Ok(OsString::new());
    }
    home_by_uid().ok_or_else(|| {
        "the environment variable HOME is not set and the user ID has no entry in the user database; NCBI BLAST+ then looks up the home directory, where it reads its registry files, by the login name (USER, LOGNAME or getlogin), which LOSAT cannot do; this is not supported by LOSAT".to_string()
    })
}

/// The program's directories on NCBI's registry search path (`registry_search_path`): the
/// directory of the program's name, then that of its resolved name when it is another one,
/// as they are before the program's arguments are set or after (`arguments_set`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:378-385
/// NCBI reference (598d8ae6): c++/src/corelib/ncbienv.cpp:373-415
/// Before the arguments are set (`CNcbiArguments(0, 0)`, ncbiapp.cpp:190) the name is "ncbi",
/// which has no directory, and the resolved name is on Linux the target of `/proc/<pid>/exe`
/// (Rust's `current_exe`, the same `readlink` of `/proc/self/exe`: a path the system made, with
/// no link in it to follow), elsewhere "ncbi" normalized (on Unix a link "ncbi" in the working
/// directory is followed).
///
/// After they are set (`m_Arguments->Reset(argc, argv, exepath, m_RealExePath)`,
/// ncbiapp.cpp:882, 1026) the name is the path that `FindProgramExecutablePath` found
/// (`unix_program_name` on Unix: for a link `blastn` on `PATH` that runs LOSAT, the link's
/// path) and the resolved name on Linux the target of `/proc/<pid>/exe` (ncbiapp.cpp:1511-1524),
/// elsewhere on Unix the name normalized with its links followed. On Windows both come from
/// `GetModuleFileName` (Rust's `current_exe`), the resolved name lexically normalized (no link
/// is read outside Unix, ncbifile.cpp:928-946), which names the same directory as the full
/// name that the loader gives, so LOSAT adds that one directory. A target that is neither
/// (WASI) has no program path: no directory.
fn program_directories(arguments_set: bool) -> Result<Vec<PathBuf>, String> {
    #[cfg(windows)]
    {
        if !arguments_set {
            return Ok(Vec::new());
        }
        Ok(std::env::current_exe()
            .ok()
            .and_then(|exe| exe.parent().map(Path::to_path_buf))
            .into_iter()
            .collect())
    }
    #[cfg(unix)]
    {
        use std::os::unix::ffi::{OsStrExt, OsStringExt};
        let argv0 = std::env::args_os().next();
        let name = match (&argv0, arguments_set) {
            (Some(argv0), true) => {
                let path_env = std::env::var_os("PATH").unwrap_or_default();
                let cwd = std::env::current_dir()
                    .map(PathBuf::into_os_string)
                    .unwrap_or_default();
                let is_file = |name: &[u8]| Path::new(OsStr::from_bytes(name)).is_file();
                unix_program_name(
                    Some(argv0.as_bytes()),
                    path_env.as_bytes(),
                    cwd.as_bytes(),
                    &is_file,
                )
            }
            // No argument at all is also "ncbi" after they are set (`m_Args` is empty).
            _ => b"ncbi".to_vec(),
        };
        let read_link = |name: &[u8]| {
            std::fs::read_link(OsStr::from_bytes(name))
                .ok()
                .map(|target| target.into_os_string().into_vec())
        };
        let following = || {
            normalize_path(&name, Some(&read_link)).ok_or_else(|| {
                format!(
                    "NCBI BLAST+ follows the links of the program's name {} 1024 times (a loop of links) and then posts a warning, which LOSAT does not write; this is not supported by LOSAT",
                    String::from_utf8_lossy(&name)
                )
            })
        };
        // Linux: `/proc/<pid>/exe`'s target; when it cannot be read, NCBI's resolved name is
        // after the arguments are set the name with its links followed (ncbiapp.cpp:1603-1607),
        // before `/proc/<pid>/exe` itself, a directory with no file.
        let resolved = match (cfg!(target_os = "linux"), std::env::current_exe()) {
            (true, Ok(exe)) => Some(exe.into_os_string().into_vec()),
            (true, Err(_)) if !arguments_set => None,
            _ => Some(following()?),
        };
        Ok(program_dirs_of(&name, resolved.as_deref())
            .into_iter()
            .map(|dir| PathBuf::from(OsString::from_vec(dir)))
            .collect())
    }
    #[cfg(not(any(unix, windows)))]
    {
        let _ = arguments_set;
        Ok(Vec::new())
    }
}

/// The directories of the program's `name` and of its `resolved` name, in order, the second
/// only when it differs (`program_directories`).
#[cfg(any(unix, test))]
fn program_dirs_of(name: &[u8], resolved: Option<&[u8]>) -> Vec<Vec<u8>> {
    let dirname = |name: &[u8]| -> Vec<u8> {
        match name.iter().rposition(|byte| b"/\\:".contains(byte)) {
            Some(at) => name[..=at].to_vec(),
            None => Vec::new(),
        }
    };
    let dir = dirname(name);
    let dir2 = resolved.map(dirname).unwrap_or_default();
    let mut dirs = Vec::new();
    if !dir.is_empty() {
        dirs.push(dir.clone());
    }
    if !dir2.is_empty() && dir2 != dir {
        dirs.push(dir2);
    }
    dirs
}

/// The program's path that NCBI finds from `argv[0]` on Unix (`m_ProgramName`): `argv0` as
/// it is when absolute; else the file of that name in the working directory `cwd`; else the
/// first file of the name's base on `PATH` (`path_env`); then normalized without reading
/// links (`normalize_path`). `is_file` is `CFile(name).Exists()`: a file, links followed.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1452-1457 (`argv[0]`, else the
/// application's arguments' name, which before they are set is "ncbi", ncbienv.cpp:384-394)
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1549-1595
/// `CDir::GetCwd` is "" when `getcwd` fails (ncbifile.cpp:3708-3716). `GetBase` is the file
/// name after the last '/' without the part from its last '.' (`SplitPath`,
/// ncbifile.cpp:358-377), of the name without its trailing '/'s (`CDirEntry::Reset`, a lone
/// "/" kept); `MakePath` puts a '/' between the directory and the base unless the directory
/// ends with one (ncbifile.cpp:413-429, 455-462).
#[cfg(any(unix, test))]
fn unix_program_name(
    argv0: Option<&[u8]>,
    path_env: &[u8],
    cwd: &[u8],
    is_file: &dyn Fn(&[u8]) -> bool,
) -> Vec<u8> {
    let given: &[u8] = match argv0 {
        Some(name) if !name.is_empty() => name,
        _ => b"ncbi",
    };
    let mut app_path = given.to_vec();
    if !app_path.starts_with(b"/") {
        if is_file(&app_path) {
            app_path = [cwd, b"/", given].concat();
            if !is_file(&app_path) {
                app_path.clear();
            }
        } else {
            // `given` does not start with '/', so it is not the lone "/" that `Reset` keeps.
            let entry = trim_trailing_slashes(given);
            let file = match entry.iter().rposition(|&byte| byte == b'/') {
                Some(at) => &entry[at + 1..],
                None => entry,
            };
            let base = match file.iter().rposition(|&byte| byte == b'.') {
                Some(at) => &file[..at],
                None => file,
            };
            app_path = path_env
                .split(|&byte| byte == b':')
                .filter(|dir| !dir.is_empty())
                .map(|dir| {
                    let mut candidate = dir.to_vec();
                    if !candidate.ends_with(b"/") {
                        candidate.push(b'/');
                    }
                    candidate.extend_from_slice(base);
                    candidate
                })
                .find(|candidate| is_file(candidate))
                .unwrap_or_default();
        }
    }
    let name = if app_path.is_empty() {
        argv0.unwrap_or_default()
    } else {
        &app_path
    };
    normalize_path(name, None).expect("a path normalized without links has no link limit")
}

/// `path` without its trailing '/'s (`CDirEntry::DeleteTrailingPathSeparator` on Unix,
/// ncbifile.cpp:465-472): empty when it has only '/'s.
fn trim_trailing_slashes(path: &[u8]) -> &[u8] {
    match path.iter().rposition(|&byte| byte != b'/') {
        Some(at) => &path[..=at],
        None => &[],
    }
}

/// `CDirEntry::NormalizePath` on Unix: `path` without "." and empty parts, a ".." taking off
/// the part before it (none at the root, kept after another ".." or at the start of a relative
/// path), and with `read_link` (`eFollowLinks`) each link replaced by its target, which
/// restarts from the root when absolute.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-971
/// `NStr::Split` with no flag keeps the empty parts. At a depth of 1024 links NCBI posts a
/// warning and stops following them: `None` then (a loop of links; the path that the system
/// resolved to run the program, within its own limit of 40 links on Linux, never reaches it,
/// but the name "ncbi" of a program before its arguments are set can, `program_directories`).
/// Without `read_link`, never `None`. rpsbproc's launch path (`rpsbproc::program_paths`,
/// `eIgnoreLinks`) is normalized here on every target.
pub(crate) fn normalize_path(
    path: &[u8],
    read_link: Option<&dyn Fn(&[u8]) -> Option<Vec<u8>>>,
) -> Option<Vec<u8>> {
    if path.is_empty() {
        return Some(Vec::new());
    }
    let mut current = trim_trailing_slashes(path).to_vec();
    if current.is_empty() {
        return Some(b"/".to_vec());
    }
    let follow = read_link;
    let mut link_depth = 0;
    let mut head: Vec<Vec<u8>> = Vec::new();
    let mut tail: std::collections::VecDeque<Vec<u8>> = std::collections::VecDeque::new();
    while !current.is_empty() || !tail.is_empty() {
        if !current.is_empty() {
            let pretail: Vec<Vec<u8>> = current
                .split(|&byte| byte == b'/')
                .map(<[u8]>::to_vec)
                .collect();
            current.clear();
            if pretail[0].is_empty() {
                head.clear();
            }
            for part in pretail.into_iter().rev() {
                tail.push_front(part);
            }
        }
        let next = tail.pop_front().unwrap_or_default();
        if let Some(last) = head.last() {
            if last == b"." {
                if !next.is_empty() {
                    head.pop();
                }
            } else if next == b"." || next.is_empty() {
                continue;
            } else if next == b".." {
                if last.is_empty() {
                    continue;
                }
                if last != b".." {
                    head.pop();
                    continue;
                }
            }
        }
        if let Some(read) = follow {
            let link = if head.is_empty() {
                next.clone()
            } else {
                [head.join(&b'/').as_slice(), b"/", &next].concat()
            };
            if let Some(target) = read(&link).filter(|target| !target.is_empty()) {
                current = target;
                link_depth += 1;
                if link_depth >= 1024 {
                    return None;
                }
                continue;
            }
        }
        head.push(next);
    }
    if head.is_empty() || (head.len() == 2 && head[0] == b"." && head[1].is_empty()) {
        return Some(b".".to_vec());
    }
    if head.len() == 1 && head[0].is_empty() {
        return Some(b"/".to_vec());
    }
    Some(head.join(&b'/'))
}

/// The registry file that NCBI loads under `file_name` (the first one on the search path that
/// is a file, `CFile::Exists`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:256-265
fn find_registry(
    search_path: &[PathBuf],
    file_name: &str,
    files: &dyn RegistryFiles,
) -> Option<PathBuf> {
    search_path
        .iter()
        .map(|dir| dir.join(file_name))
        .find(|path| files.is_file(path))
}

/// Why the file of `-conffile` is not loaded (`conf_file_path`).
#[derive(Debug, PartialEq, Eq)]
enum ConfFile {
    /// NCBI's post that it is not a file (`ERR_POST_X(23, Critical ...)`), after which the
    /// application aborts.
    NotFound(String),
    /// NCBI's report of a name with a directory whose file cannot be opened (`CAppException`
    /// `eNoRegistry`, thrown by `LoadConfig`: the configuration load fails, exit status 2).
    CannotOpen(String),
    /// LOSAT's refusal: a name that is not UTF-8, whose bytes NCBI's post prints.
    Refused(String),
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1014-1022 (after the words before the
// arguments are read; `-version` there ends the program first)
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1262-1264,1267-1270 (`LoadConfig` with that file)
// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:256-269:x_FindRegistry (a name without a
// directory is looked for in the search path; a name with one is that file)
/// The registry file that the application loads in place of `<program>.ini` for `-conffile
/// name` (`Ok(None)`: none, the load goes on without one), or why it does not load one. `name`
/// must be a file where the program runs (`CFile::Exists`: a directory, an empty name or a
/// missing file is NCBI's post and abort, rc 255). A name without a directory is then looked for
/// in the registry search path, as `<program>.ini` is (`find_registry`): with
/// `NCBI_DONT_USE_LOCAL_CONFIG` the current directory is not in it, and NCBI loads no file (NCBI
/// BLAST+ 2.17.0, probe `conf_nolocal`, line RG `E/SURVEY_E1.md`); a name with one (`./c.ini`) is
/// that file.
fn conf_file_path(
    name: &OsStr,
    search_path: &[PathBuf],
    files: &dyn RegistryFiles,
) -> Result<Option<PathBuf>, ConfFile> {
    let Some(text) = name.to_str() else {
        return Err(ConfFile::Refused(format!(
            "the -conffile name {:?} is not UTF-8; NCBI BLAST+ reads its bytes, which is not supported by LOSAT",
            name
        )));
    };
    let path = Path::new(name);
    if !files.is_file(path) {
        return Err(ConfFile::NotFound(format!(
            "Critical: (106.23) Registry file \"{text}\" not found, aborting.\n"
        )));
    }
    let has_directory = path
        .parent()
        .is_some_and(|parent| !parent.as_os_str().is_empty());
    if has_directory {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1266-1281 (`LoadConfig`: a
        // name with a directory that gives no registry throws; one without is an `Info` post)
        // NCBI BLAST+ 2.17.0 (probes `conf_unread_dot`, `conf_unread_bare`, mode 000): `./c.ini`
        // ends the program with this report (exit status 2), `c.ini` loads no file.
        if !files.opens(path) {
            let message = format!("Registry file \"{text}\" cannot be opened");
            // NCBI's `Error: (CAppException::eNoRegistry) <message>` and `... initialization failed
            // (CAppException::eNoRegistry) <message>`, in LOSAT's words (`crate::wording`: without
            // the class and code, `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12).
            return Err(ConfFile::CannotOpen(format!(
                "Error: {message}\nError: (106.15) Application's initialization failed: {message}\n"
            )));
        }
        return Ok(Some(path.to_path_buf()));
    }
    Ok(find_registry(search_path, text, files))
}

// NCBI reference (598d8ae6): c++/include/corelib/metareg.hpp:64-69
// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:284-303 (`eName_Ini` adds `.ini`,
// `eName_DotRc` makes `.<name>rc`)
/// The name of NCBI's system-wide registry file, which `CMetaRegistry::Load("ncbi",
/// eName_RcOrIni)` looks for: `.ncbirc`, and `ncbi.ini` on Windows.
const SYSTEM_WIDE_REGISTRY: &str = if cfg!(windows) { "ncbi.ini" } else { ".ncbirc" };

/// One read of the registry file at `path` into `layer` (`ncbi_registry::read_registry`):
/// `Ok(None)` when NCBI cannot open it, which NCBI takes for no file (and looks no further).
/// `Err` holds LOSAT's refusal.
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:60-90 (`SEntry::Reload`, which
/// `x_Load` calls for the file found)
/// NCBI BLAST+ 2.17.0 (g_reg probes `unr`, `unr2`, 2026-11-11): a malformed `.ncbirc` that it
/// cannot open (mode 000) is not reported, and a `.ncbirc` in the home directory is then not
/// read either; a `blastn.ini` that it cannot open is no file. The file is opened in binary mode:
/// the line ends are `NcbiGetlineEOL`'s alone.
fn read_registry_file(
    path: &Path,
    layer: Layer,
    files: &dyn RegistryFiles,
) -> Result<Option<RegistryRead>, String> {
    let Some(read) = files.read(path) else {
        return Ok(None);
    };
    let bytes = read.map_err(|error| {
        format!(
            "NCBI BLAST+ would read the registry file {} ({error}), which LOSAT cannot read; this is not supported by LOSAT",
            path.display()
        )
    })?;
    read_registry(&bytes, layer).map(Some).map_err(|reason| {
        format!(
            "the registry file {} {reason}; this is not supported by LOSAT",
            path.display()
        )
    })
}

/// The first entry of `entries` that LOSAT does not know to change no output
/// (`is_harmless_entry`), as LOSAT's refusal of the file at `path`.
fn check_registry_entries(path: &Path, entries: &[RegistryEntry]) -> Result<(), String> {
    for RegistryEntry { section, name, .. } in entries {
        if !is_harmless_entry(section, name) {
            return Err(format!(
                "the registry file {} sets [{section}] {name}, which LOSAT does not know to change no output of NCBI BLAST+ (LOSAT accepts only [BLAST] BLASTDB, BLASTMAT, DATA_LOADERS, BLASTDB_NUCL_DATA_LOADER, BLASTDB_PROT_DATA_LOADER, IGDATA, MAX_SEQID_LENGTH, BLAST_USAGE_REPORT, [NCBI] DATA and DONT_USE_NCBIRC, [WINDOW_MASKER] WINDOW_MASKER_PATH); this is not supported by LOSAT",
                path.display()
            ));
        }
    }
    Ok(())
}

/// `.ncbirc` as the registry of `-window_masker_taxid` reads it (`window_masker_path_of`): its
/// `[WINDOW_MASKER] WINDOW_MASKER_PATH`, and what each load of the file writes. It is the read of
/// the file that the application start made (`read_registry_file` into `Layer::SystemWide`, which
/// is how that registry loads it, also when the program's `.ini` turned `.ncbirc` off for the
/// application). A file that NCBI reports as malformed is read up to the line it reports, as at
/// the start. NCBI loads that registry for each query factory of the search (`WindowMasker`'s
/// load posts), and `CMetaRegistry` keeps no read that threw: each load reads the file again and
/// posts the syntax error, now with the program's prefix (`running_syntax_error_notice`,
/// `Critical`), after the warning of each entry above it that sets an entry again
/// (`repeated_entry_post`, which the start's post level hides). It keeps a read that did not
/// throw (a stream that failed, `Error`, or entries set again): a load then posts only when no
/// read of the file was made before with the registry's flags (`kept`: those of the usage
/// report's registry, and of the application's load without a `.ini`), and only the first
/// (NCBI BLAST+ 2.17.0 on 2026-11-12, line NP M3c: 3 query batches give 3 `Critical` posts of
/// a syntax error, and one `Error reading the registry` of a one-byte `.ncbirc` with a
/// `blastn.ini` and `BLAST_USAGE_REPORT=0`, none without the `.ini` or with the usage report;
/// its review's `d1` and `d3`: one `Found multiple` warning of a well-formed file so read, and
/// one before each `Critical` of a malformed file).
/// `Err` is LOSAT's refusal, which only `-window_masker_taxid` gives: a file LOSAT cannot read,
/// a value that is not UTF-8, and `[NCBI] .Inherits`, whose files NCBI loads into the same layer.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1650-1662
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:200-209 (cited at
/// `application_start`: a read kept for the same name, style and flags is not read again)
fn window_masker_ncbirc(
    ncbirc: Option<&(Result<Option<RegistryRead>, String>, PathBuf)>,
    program: &str,
    kept: bool,
) -> Result<NcbircWindowMasker, String> {
    let Some((read, path)) = ncbirc else {
        return Ok(NcbircWindowMasker::default());
    };
    let Some(read) = read.as_ref().map_err(String::clone)? else {
        return Ok(NcbircWindowMasker::default());
    };
    if read
        .entries
        .iter()
        .any(|RegistryEntry { section, name, .. }| {
            section.eq_ignore_ascii_case("NCBI") && name.eq_ignore_ascii_case(".Inherits")
        })
    {
        return Err(format!(
            "the registry file {} sets [NCBI] .Inherits, whose files NCBI BLAST+ reads for -window_masker_taxid; this is not supported by LOSAT",
            path.display()
        ));
    }
    // What a read of the file posts: a warning of each entry that sets an entry again, then the
    // error that ended the read, or that of a stream that failed.
    let mut load_posts = Posts::new();
    if read.error.is_some() || !kept {
        for repeated in &read.repeated {
            load_posts.post(WARNING, repeated_entry_post(repeated, program));
        }
    }
    let posts_once = match &read.error {
        Some(error) => {
            load_posts.post(CRITICAL, running_syntax_error_notice(error, program));
            false
        }
        None if kept => false,
        None => {
            if read.stream_failed {
                load_posts.post(ERROR, running_error_reading_post(program));
            }
            true
        }
    };
    Ok(NcbircWindowMasker {
        value: file_text_entry(path, &read.entries, "WINDOW_MASKER", "WINDOW_MASKER_PATH")?,
        load_posts,
        posts_once,
    })
}

/// What `window_masker_ncbirc` reads of `.ncbirc` (no file: no value, no posts).
#[derive(Default)]
struct NcbircWindowMasker {
    value: Option<String>,
    load_posts: Posts,
    posts_once: bool,
}

/// The registry that blastn's `-window_masker_taxid` builds to find its files
/// (`window_masker_path_of`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct WindowMaskerRegistry {
    /// The directory where `-window_masker_taxid` looks for `<taxid>/wmasker.obinary`: `None`
    /// when no layer has the entry; NCBI uses the current directory for both `None` and an
    /// empty value.
    pub path: Option<String>,
    /// What NCBI posts each time it loads the registry (`window_masker_ncbirc`): the posts of a
    /// `.ncbirc` that is malformed or sets an entry again, or nothing.
    pub load_posts: Posts,
    /// Whether only the first load writes `load_posts` (a read that NCBI keeps).
    pub posts_once: bool,
}

/// The NCBI application settings that LOSAT reproduces: the command line's from its application
/// start (`ApplicationStart::settings`), the library's from the start of its call (`api::start`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ApplicationSettings {
    /// `SDataLoaderConfig::UseDataLoaders()` of the query and subject readers: whether
    /// NCBI's reader tries a record's first line as a sequence identifier
    /// (`fasta_reader::ReaderConfig::data_loaders`).
    pub data_loaders: bool,
    /// The registry entry `[BLAST] BLASTDB`, the last part of the path where a database
    /// name is looked for (`seqdb::SearchPath::new`): `Ok(None)` when the registry has no
    /// such entry, `Ok(Some(value))` with the value NCBI's registry gives. `Err(reason)`
    /// when a registry file sets the entry to a value that LOSAT does not read as NCBI
    /// reads it; only a program that opens a database has to refuse to run then, since
    /// the entry changes nothing else.
    pub blastdb: Result<Option<String>, String>,
    /// The first of `SEQ_ID_RANK_VARIABLES` that is set (with any value). NCBI's choice of the
    /// Seq-ids of a database with parsed Seq-ids depends on it; only a program that opens such
    /// a database has to refuse to run then (`algorithm::common::database::
    /// check_seqids`), since the ids of the other searches are local or ordinal.
    pub seq_id_rank_variable: Option<&'static str>,
    /// The value of the registry entry `[BLAST] MAX_SEQID_LENGTH` in the first layer that has
    /// the entry (the environment, then the program's `.ini` file, then `.ncbirc`), an empty
    /// value too: `registry.Get` answers from that layer, and an empty value there hides the
    /// layers below (makeblastdb 2.17.0 with an empty `NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH`
    /// or an empty `makeblastdb.ini` entry ignores a `.ncbirc` value). `Err` holds why LOSAT
    /// does not read a file's value as NCBI stores it. The entry changes no output of a
    /// search (`HARMLESS_ENTRIES`); makeblastdb gives a value it can convert to its FASTA
    /// reader as the limit of a Seq-id's length, which `losat makeblastdb -parse_seqids`
    /// refuses (`writedb::makeblastdb::check_registry_settings`).
    pub max_seqid_length: Result<Option<String>, String>,
    /// Where NCBI looks for the file of a matrix that its tables do not have
    /// (`ApplicationSettings::matrix_file`).
    pub matrix_search: MatrixSearch,
    /// The registry of blastn's `-window_masker_taxid` (`window_masker_path_of`): the directory
    /// of its files and what each load of the registry writes. `Err(reason)` when LOSAT does
    /// not read it as NCBI does; only `-window_masker_taxid` has to refuse to run then.
    pub window_masker_registry: Result<WindowMaskerRegistry, String>,
}

/// The parameters that move the GIs to the end of NCBI's rank of Seq-ids (the Seq-align id of a
/// database sequence and every choice by `Score`): `SEQ_ID_PREFER_ACCESSION_OVER_GI` and
/// `SEQ_ID_AVOID_GI` (PLAN of line DBID §9.2, DECISIONS 3 Q8: LOSAT follows their defaults).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3393-3413
/// (The registry entries `[SeqId] PreferAccessionOverGi` and `AvoidGi` and the variables
/// `NCBI_CONFIG__SEQID__...` set the same parameters; LOSAT rejects both as entries it does not
/// read.)
pub const SEQ_ID_RANK_VARIABLES: [&str; 2] = ["SEQ_ID_PREFER_ACCESSION_OVER_GI", "SEQ_ID_AVOID_GI"];

impl Default for ApplicationSettings {
    /// No registry entry and no variable: both data loaders on, no `[BLAST] BLASTDB`.
    ///
    /// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_scope_src.hpp:73,81
    fn default() -> Self {
        Self {
            data_loaders: true,
            blastdb: Ok(None),
            seq_id_rank_variable: None,
            max_seqid_length: Ok(None),
            matrix_search: MatrixSearch::default(),
            window_masker_registry: Ok(WindowMaskerRegistry::default()),
        }
    }
}

/// The settings by which NCBI looks for the file of a protein matrix that its tables do not have
/// (`ApplicationSettings::matrix_file`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MatrixSearch {
    /// `NCBI_DATA_PATH`, the parameter `[NCBI] DataPath` (a registry file's or
    /// `NCBI_CONFIG__NCBI__DATAPATH`'s value is refused, `HARMLESS_ENTRIES`).
    data_path: Option<OsString>,
    /// The parameter `[NCBI] Data`: `NCBI_CONFIG__NCBI__DATA`, else the registry's entry
    /// (`NCBI_CONFIG__NCBI__Data`, then the program's `.ini` file, then `.ncbirc`, as
    /// `blastdb_of`). `Err` holds why LOSAT does not read a file's value as NCBI stores it.
    data: Result<Option<OsString>, String>,
    /// `BLASTMAT`.
    blastmat: Option<OsString>,
}

impl Default for MatrixSearch {
    fn default() -> Self {
        Self {
            data_path: None,
            data: Ok(None),
            blastmat: None,
        }
    }
}

impl ApplicationSettings {
    /// The file that NCBI BLAST+ reads for `-matrix name` when the name is not one of its
    /// tables' (an ungapped search's: a gapped one refuses the name first), or `None` when it
    /// finds none or cannot open the one it finds (its search then fails with "Unknown error
    /// code -1"). `Err` when a registry file's `[NCBI] Data` decides and LOSAT cannot read it.
    /// NCBI reads the file as a matrix: 2.17.0 aborts after a matrix that it reads (the length
    /// adjustment's tables of the name are unset, `TblastxMatrix::Aborts`) and fails with the
    /// reader's error code on another file; LOSAT reads no matrix file, and a program refuses
    /// the run in both cases (`matrix_file_refusal`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:381-387 (the score
    /// block's name of a protein matrix is upper-cased)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1631-1652
    /// (`full_matrix_path` is the path found: `s_GetCStringOfMatrixPath` cuts the name off it,
    /// blast_setup_cxx.cpp:1288-1298.)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1381
    /// (`BlastFindMatrixPath`, `get_path` above: the name and its upper case, which are the
    /// same here; `is_prot` is FALSE, so the subdirectory is "nt")
    /// NCBI reference (598d8ae6): c++/src/util/util_misc.cpp:129-179 (`g_FindDataFile`; a
    /// parameter's variable, when set, is its value, an empty one too, ncbi_param.cpp:400-437)
    pub fn matrix_file(&self, name: &str) -> Result<Option<PathBuf>, String> {
        self.matrix_search.file(name)
    }
}

impl MatrixSearch {
    /// `ApplicationSettings::matrix_file` with these settings: the file that NCBI reads for the
    /// matrix `name`, which a search of the library looks up with the settings of its start
    /// (`api::start::CallStart::settings`), as the command line does with its own.
    pub fn file(&self, name: &str) -> Result<Option<PathBuf>, String> {
        Ok(find_matrix_file(name, self, &Path::is_file, &Path::is_dir)?
            .filter(|file| std::fs::File::open(file).is_ok()))
    }
}

/// `ApplicationSettings::matrix_file` with `is_file` (`CFile::Exists`: a file, links followed)
/// and `is_dir` (`CDir::Exists`).
fn find_matrix_file(
    name: &str,
    search: &MatrixSearch,
    is_file: &dyn Fn(&Path) -> bool,
    is_dir: &dyn Fn(&Path) -> bool,
) -> Result<Option<PathBuf>, String> {
    let upper = name.to_ascii_uppercase();
    // `kDelim` (util_misc.cpp:141-145), `CFile::GetPathSeparator` and the separators of
    // `AddTrailingPathSeparator` (ncbifile.cpp:105-120, 455-462).
    let (delimiter, separator, separators): (char, &str, &[u8]) = if cfg!(windows) {
        (';', "\\", b":/\\")
    } else {
        (':', "/", b"/")
    };
    let concat = |parts: &[&OsStr]| -> PathBuf {
        let mut path = OsString::new();
        parts.iter().for_each(|part| path.push(part));
        PathBuf::from(path)
    };
    let found = |path: PathBuf| is_file(&path).then_some(path);
    // `IsAbsolutePath` (ncbifile.cpp:508-523, `s_Win_IsDiskPath` and `s_Win_IsNetworkPath` on
    // Windows).
    let bytes = upper.as_bytes();
    let slash = |byte: u8| byte == b'/' || byte == b'\\';
    let absolute = if cfg!(windows) {
        (bytes.len() >= 3 && bytes[0].is_ascii_alphabetic() && bytes[1] == b':' && slash(bytes[2]))
            || (bytes.len() >= 2 && slash(bytes[0]) && slash(bytes[1]))
    } else {
        bytes.first() == Some(&b'/')
    };
    let mut dirs: Vec<OsString> = Vec::new();
    if absolute {
        dirs.push(OsString::new());
    } else {
        if let Some(data_path) = &search.data_path {
            dirs.extend(
                split_os(data_path, delimiter)
                    .into_iter()
                    .filter(|dir| !dir.is_empty()),
            );
        }
        if let Some(data) = search.data.clone()? {
            if !data.is_empty() {
                dirs.push(data);
            }
        }
    }
    for dir in &dirs {
        // `MakePath(dir, name)`: a separator between them unless the directory ends with one.
        let between = match dir.as_encoded_bytes().last() {
            Some(last) if !separators.contains(last) => separator,
            _ => "",
        };
        let candidate = concat(&[dir, OsStr::new(between), OsStr::new(&upper)]);
        if let Some(file) = found(candidate) {
            return Ok(Some(file));
        }
    }
    let separator = OsStr::new(separator);
    let upper = OsStr::new(&upper);
    if let Some(blastmat) = search.blastmat.as_deref() {
        if is_dir(Path::new(blastmat)) {
            for candidate in [
                concat(&[blastmat, separator, upper]),
                concat(&[blastmat, separator, OsStr::new("nt"), separator, upper]),
            ] {
                if let Some(file) = found(candidate) {
                    return Ok(Some(file));
                }
            }
        }
    }
    Ok(found(concat(&[OsStr::new("data"), separator, upper])))
}

/// `value` split at each `delimiter`, empty parts kept.
fn split_os(value: &OsStr, delimiter: char) -> Vec<OsString> {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        value
            .as_bytes()
            .split(|&byte| char::from(byte) == delimiter)
            .map(|part| OsStr::from_bytes(part).to_os_string())
            .collect()
    }
    #[cfg(not(unix))]
    {
        value
            .to_string_lossy()
            .split(delimiter)
            .map(OsString::from)
            .collect()
    }
}

/// The refusal of `program`'s run whose `-matrix name` NCBI BLAST+ reads from `file`
/// (`ApplicationSettings::matrix_file`).
pub fn matrix_file_refusal(name: &str, file: &Path, program: &str) -> String {
    format!(
        "-matrix {name} is not a matrix of NCBI BLAST+'s tables, and NCBI BLAST+ reads it from the file {} (2.17.0 then aborts without the matrix's statistics table, or fails on a file that is not a matrix); this is not supported by LOSAT's {program}",
        file.display()
    )
}

/// `NStr::FindNoCase(value, word) != NPOS` for an ASCII `word` (case folded in the C
/// locale, where bytes above 0x7f fold to themselves).
fn contains_no_case(value: &[u8], word: &[u8]) -> bool {
    value
        .windows(word.len())
        .any(|window| window.eq_ignore_ascii_case(word))
}

/// Whether the data loaders are used, from the values of `[BLAST] DATA_LOADERS` in the
/// layers of NCBI's registry, highest priority first: the environment
/// (`NCBI_CONFIG__BLAST__DATA_LOADERS`, `env`), then the registry files (`files`: the
/// program's `.ini` file, then `.ncbirc`); `None`: the layer has no such entry.
///
/// A variable that is set is an entry of the environment layer even when its value is
/// empty, so an empty `NCBI_CONFIG__BLAST__DATA_LOADERS=` answers first and turns both
/// loaders off (it contains neither `blastdb` nor `genbank`; audit finding B-1 of session
/// SFc, `AUTHORITY.md` §G1). In a registry file an empty value is no entry
/// (`CMemoryRegistry::x_HasEntry`, ncbireg.cpp:984-991; confirmed with the oracle in
/// session SFb), so the next layer decides.
///
/// NCBI reference (598d8ae6): c++/src/corelib/env_reg.cpp:157-167
/// NCBI reference (598d8ae6): c++/src/corelib/ncbienv.cpp:97-104,117-124
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:984-991
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:75-92
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_scope_src.hpp:111
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1577-1583
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1235-1246
fn data_loaders_of(env: Option<&[u8]>, files: [Option<&[u8]>; 2]) -> bool {
    let Some(value) = env.or_else(|| files.into_iter().flatten().find(|value| !value.is_empty()))
    else {
        return true;
    };
    let mut use_blast_dbs = true;
    let mut use_genbank = true;
    if !contains_no_case(value, b"blastdb") {
        use_blast_dbs = false;
    }
    if !contains_no_case(value, b"genbank") {
        use_genbank = false;
    }
    if contains_no_case(value, b"none") {
        use_blast_dbs = false;
        use_genbank = false;
    }
    use_blast_dbs || use_genbank
}

/// `last_value` of `[section] entry` as text, for a setting that LOSAT reads as text: `Err`
/// (LOSAT's refusal) for a value that is not UTF-8, as for the environment's variable.
fn file_text_entry(
    path: &Path,
    entries: &[RegistryEntry],
    section: &str,
    entry: &str,
) -> Result<Option<String>, String> {
    last_value(entries, section, entry)
        .map(|value| {
            String::from_utf8(value.to_vec()).map_err(|_| {
                format!(
                    "the registry file {} sets [{section}] {entry} to a value that is not valid UTF-8; this is not supported by LOSAT",
                    path.display()
                )
            })
        })
        .transpose()
}

/// A registry file's `[BLAST] MAX_SEQID_LENGTH` when the file has the entry (an empty value
/// too; `Err` for a value LOSAT cannot read, `file_text_entry`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:739-747
fn max_seqid_length_entry(
    path: &Path,
    entries: &[RegistryEntry],
) -> Option<Result<String, String>> {
    file_text_entry(path, entries, "BLAST", "MAX_SEQID_LENGTH").transpose()
}

/// The registry entry `[BLAST] BLASTDB` from the layers of the application's registry,
/// highest priority first: the environment (`NCBI_CONFIG__BLAST__BLASTDB`), then the
/// registry files (the program's `.ini` file, then `.ncbirc`). The first layer that has the
/// entry answers, also with an empty value: unlike `[BLAST] DATA_LOADERS`
/// (`data_loaders_of`), an entry of a registry file whose value is empty is an entry here
/// and hides the layers below it. A file's value that LOSAT cannot read (`Err`) ends the
/// search only when no layer above it has the entry.
///
/// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:527-533
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1290-1304
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1239-1244
/// `FindByContents` asks every layer with `fCountCleared`, under which a registry file's
/// entry with an empty value counts (`CMemoryRegistry::x_HasEntry`, ncbireg.cpp:984-991,
/// cited at `data_loaders_of`); the file reader stores such an entry
/// (`CMemoryRegistry::x_Set`, ncbireg.cpp:1066-1070: removing it is "historic behavior;
/// could inappropriately expose entries in lower layers").
fn blastdb_of(
    env: Option<Result<String, String>>,
    files: [Result<Option<String>, String>; 2],
) -> Result<Option<String>, String> {
    if let Some(value) = env {
        return value.map(Some);
    }
    for file in files {
        if let Some(value) = file? {
            return Ok(Some(value));
        }
    }
    Ok(None)
}

/// Where NCBI's diagnostics go once an application has set them up with its `-logfile`
/// (`SetLogFile`).
pub enum DiagTarget {
    /// The standard stream `cerr`: no `-logfile`, `-logfile -`, or a file that cannot be
    /// opened.
    Stderr,
    /// The file of `-logfile`, appended to.
    File(std::fs::File),
    /// Nowhere: `-logfile ""` and `-logfile /dev/null` (a stream handler without a stream).
    Discard,
    /// The file of `-logfile` of a run whose files LOSAT Web gives its host
    /// (`report::run_files::RunFiles::Receiver`): its bytes, kept until the run ends. Behind a
    /// lock, so that every target is written through a shared reference, as a file is
    /// (`Write for &File`), and a target stays `Sync` and `RefUnwindSafe` (`api::Batches` holds
    /// one; a `RefCell` would take both away from that public type).
    Memory(std::sync::Mutex<Vec<u8>>),
}

/// Where the text written to a `DiagTarget` ends up (`DiagTarget::place`, which `write_text`
/// follows). No NCBI counterpart: LOSAT Web tells its host where the error of an output that a run
/// could not write is (line WQA's sweep, D7, 2026-11-20).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TextPlace {
    /// The standard error (`cerr`): the caller's diagnostics.
    StandardError,
    /// The file of `-logfile`, on disk or in memory.
    LogFile,
    /// Nowhere (`-logfile ""` and `-logfile /dev/null`).
    Nowhere,
}

impl DiagTarget {
    /// Where the text written to this target ends up (`write_text`).
    pub fn place(&self) -> TextPlace {
        match self {
            DiagTarget::Stderr => TextPlace::StandardError,
            DiagTarget::File(_) | DiagTarget::Memory(_) => TextPlace::LogFile,
            DiagTarget::Discard => TextPlace::Nowhere,
        }
    }

    /// Writes `text` where the diagnostics go (`place`; `stderr` for `cerr`). A shared reference is
    /// enough: makeblastdb writes its log through the target's file while its diagnostics go to the
    /// target (`writedb::makeblastdb::run_with`).
    pub fn write_text(&self, text: &[u8], stderr: &mut dyn Write) {
        match self.place() {
            TextPlace::LogFile => self.write_to_file(text),
            TextPlace::StandardError => {
                let _ = stderr.write_all(text);
                let _ = stderr.flush();
            }
            TextPlace::Nowhere => {}
        }
    }

    /// Writes `text` to the file of `-logfile` (on disk or in memory); nothing for `cerr` and for
    /// no output.
    pub fn write_to_file(&self, text: &[u8]) {
        match self {
            DiagTarget::File(file) => {
                // `Write for &File`: the one open file, written through a shared reference.
                let mut file: &std::fs::File = file;
                let _ = file.write_all(text);
            }
            DiagTarget::Memory(bytes) => bytes
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .extend_from_slice(text),
            DiagTarget::Stderr | DiagTarget::Discard => {}
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:91-95 (likewise
// rpsblast_app.cpp:112-114: the handler is made in `Run`, from the handler that `-logfile` set up)
// NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:111-114 (and :241-244 for
// `x_RunMTBySplitQuery`; both are called from `Run`)
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_aux.hpp:252
/// Whether the messages of `program`'s run go where its diagnostics go (the file of
/// `-logfile`): tblastx, rpsblast and rpstblastn, whose application makes its
/// `CBlastAppDiagHandler` in `Run`, after `-logfile` set up the diagnostics. blastn, blastp, blastx, tblastn and psiblast
/// hold theirs as a member, made with the application before `AppMain` from the handler of
/// `cerr` (`replays_to_stderr`): their run's messages go to `cerr` whatever the `-logfile`.
/// NCBI BLAST+ 2.17.0 (line RG `E/SURVEY_E1.md`): `-max_target_seqs 1`'s warning, a missing
/// `-query` and `-evalue 0` go to `log.txt` in tblastx and rpsblast and to stderr in the others;
/// rpstblastn as tblastx (PS12a's probes of 2026-11-13, `PR/ws/ps12a/probe_e1_rpstbn.py`).
pub fn run_messages_follow_log(program: &str) -> bool {
    matches!(program, "tblastx" | "rpsblast" | "rpstblastn")
}

impl DiagTarget {
    /// Makes this process's standard error the target, for the rest of the process: every
    /// message of the run then goes where the diagnostics go (`run_messages_follow_log`). A
    /// target of `cerr` changes nothing; `Discard` writes to the null device. `Err` with LOSAT's
    /// refusal on a platform where the standard error cannot be replaced (`wasi`).
    pub fn into_standard_error(self, program: &str) -> Result<(), String> {
        let refusal = |reason: String| {
            format!(
                "{reason}; NCBI BLAST+ {program} writes the messages of its run to the -logfile file, which is not supported by LOSAT here"
            )
        };
        let file = match self {
            DiagTarget::Stderr => return Ok(()),
            DiagTarget::File(file) => file,
            DiagTarget::Discard => {
                let null = if cfg!(windows) { "NUL" } else { "/dev/null" };
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(null)
                    .map_err(|error| refusal(format!("{null}: {error}")))?
            }
            DiagTarget::Memory(_) => {
                return Err(refusal("the file is kept in memory".to_string()));
            }
        };
        replace_standard_error(file).map_err(refusal)
    }
}

/// The standard error that the process had before `replace_standard_error` (the terminal),
/// where LOSAT's own errors still go (`terminal_stderr`).
static TERMINAL_STDERR: std::sync::OnceLock<std::fs::File> = std::sync::OnceLock::new();

/// The standard error that the process had before its run's messages went to the `-logfile`
/// target (`DiagTarget::into_standard_error`), or `None` when it was not replaced.
pub fn terminal_stderr() -> Option<&'static std::fs::File> {
    TERMINAL_STDERR.get()
}

/// No NCBI counterpart: the standard error (fd 2) becomes a copy of `file`'s descriptor; a copy
/// of the old one is kept (`terminal_stderr`).
#[cfg(unix)]
fn replace_standard_error(file: std::fs::File) -> Result<(), String> {
    let old = rustix::io::dup(std::io::stderr())
        .map_err(|error| format!("the standard error: {error}"))?;
    rustix::stdio::dup2_stderr(&file).map_err(|error| format!("the standard error: {error}"))?;
    let _ = TERMINAL_STDERR.set(std::fs::File::from(old));
    Ok(())
}

/// No NCBI counterpart: the process's standard error handle becomes `file`'s, which is kept
/// open for the rest of the process (Rust's standard error reads the handle at each write).
#[cfg(windows)]
fn replace_standard_error(file: std::fs::File) -> Result<(), String> {
    use std::os::windows::io::{AsRawHandle, FromRawHandle, IntoRawHandle};
    let old = std::io::stderr().as_raw_handle();
    let handle = file.into_raw_handle();
    // SAFETY: `handle` is an open file handle that this process owns and never closes.
    let done = unsafe {
        windows_sys::Win32::System::Console::SetStdHandle(
            windows_sys::Win32::System::Console::STD_ERROR_HANDLE,
            handle as windows_sys::Win32::Foundation::HANDLE,
        )
    };
    if done == 0 {
        return Err(format!(
            "the standard error: {}",
            std::io::Error::last_os_error()
        ));
    }
    if !old.is_null() {
        // SAFETY: `old` is the process's former standard error handle, which nothing else
        // closes; the static keeps it for the rest of the process.
        let _ = TERMINAL_STDERR.set(unsafe { std::fs::File::from_raw_handle(old) });
    }
    Ok(())
}

#[cfg(not(any(unix, windows)))]
fn replace_standard_error(_file: std::fs::File) -> Result<(), String> {
    Err("this platform cannot replace the standard error of a process".to_string())
}

/// No NCBI counterpart: a `DiagTarget` as a stream, each write going where the diagnostics go
/// (`DiagTarget::write_text`; `stderr` stands for `cerr`). It holds the target by a shared
/// reference, so that a log can write through the target's file at the same time
/// (`writedb::makeblastdb::run_with`).
pub struct DiagStream<'a> {
    pub target: &'a DiagTarget,
    pub stderr: &'a mut dyn Write,
}

impl Write for DiagStream<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.target.write_text(buf, self.stderr);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7570-7608
/// The target of the `-logfile` value `name` (a file NCBI fails to open leaves `cerr`). The file
/// is opened (made when it does not exist) at once, as NCBI's first `SetupDiag` does, before
/// `-version` (NCBI BLAST+ 2.17.0 makeblastdb and blastdbcmd `-version -logfile F` make `F`:
/// g_reg probes `ver_log_ok`, `ver_log_ok2`).
pub fn diag_target(name: Option<&OsStr>) -> DiagTarget {
    match (special_log_target(name), name) {
        (Some(target), _) => target,
        (None, None) => DiagTarget::Stderr,
        (None, Some(name)) => match std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(name)
        {
            Ok(file) => DiagTarget::File(file),
            Err(_) => DiagTarget::Stderr,
        },
    }
}

/// The target of NCBI's special `-logfile` names (cited at `diag_target`): `cerr` without a name
/// and for `-`, no output for `""` and `/dev/null`; `None` for the name of a file. The one reading
/// of these names: `diag_target` (the file system) and `report::run_files::RunFiles::log_target`
/// (LOSAT Web).
pub fn special_log_target(name: Option<&OsStr>) -> Option<DiagTarget> {
    match name {
        None => Some(DiagTarget::Stderr),
        Some(name) if name == "-" => Some(DiagTarget::Stderr),
        Some(name) if name.is_empty() || name == "/dev/null" => Some(DiagTarget::Discard),
        Some(_) => None,
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-866
/// The file of `-logfile` (`m_LogFile`) in the arguments `words` (after the program's name):
/// the value of its last `-logfile FILE` or `-logfile=FILE` before `--`.
pub fn log_file<S: AsRef<OsStr>>(words: &[S]) -> Option<&OsStr> {
    let mut log = None;
    let mut at = 0;
    while at < words.len() {
        let word = words[at].as_ref();
        if word == "--" {
            break;
        }
        if word == "-logfile" {
            at += 1;
            if let Some(name) = words.get(at) {
                log = Some(name.as_ref());
            }
        } else if let Some(name) = word.as_encoded_bytes().strip_prefix(b"-logfile=") {
            // SAFETY: `name` is the rest of `word` after an ASCII prefix, so it starts at a
            // boundary of the platform's encoding and ends where `word` ends.
            log = Some(unsafe { OsStr::from_encoded_bytes_unchecked(name) });
        }
        at += 1;
    }
    log
}

/// What NCBI's application reads of its arguments before its argument reader
/// (`CNcbiApplicationAPI::AppMain`, `preparse`).
#[derive(Debug, PartialEq, Eq)]
pub enum Preparsed<S> {
    /// `-version` or a full version (`-version-full`, `-version-full-xml`,
    /// `-version-full-json`): print it, exit status 0 (`ncbi_version::print`).
    Version(crate::ncbi_version::VersionRequest),
    /// The words left for the argument reader, the registry file of the last `-conffile FILE`
    /// or `-conffile=FILE` (`conf_arg`, which NCBI loads in place of `<program>.ini`), and
    /// whether there was `-dryrun` (NCBI reads the arguments, then runs nothing: the command
    /// line's entry, `cli::main_entry`; rpsbproc runs its own `DryRun` after its `Init`,
    /// `rpsbproc::dry_run`; a call of the library API ends after its start,
    /// `api::start::CallStart::dry_run`; the Web engine refuses it).
    Words {
        words: Vec<S>,
        conffile: Option<OsString>,
        dryrun: bool,
    },
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1015
/// The words NCBI's application takes before its argument reader: `-version` and the full
/// versions anywhere before `--` (but not as the value of `-logfile` or `-conffile`; the first
/// one is printed, line RG item E3), `-dryrun` and `-conffile`. A `-logfile` or `-conffile`
/// without a value is dropped. `-dryrun` ends the run after the arguments (`cli::main_entry`;
/// rpsbproc's after its `Init`, `rpsbproc::dry_run`); the file of the last `-conffile` is the
/// application start's configuration (`ApplicationStart::of_process`, line RG item E1).
/// `words` are the arguments after the program's name. This is the one reading of these words:
/// the command line's entry (the versions, before the configuration load, and the words that
/// the argument reader gets, `cli::try_parse_from_start`), `ApplicationStart::of_process`,
/// `cli::ncbi_version_request` and `ncbi_app` read them here, and `log_file` is NCBI's first scan
/// of `-logfile`.
pub fn preparse<S: AsRef<OsStr> + Clone>(words: &[S]) -> Preparsed<S> {
    let mut kept = Vec::with_capacity(words.len());
    let mut conffile = None;
    let mut dryrun = false;
    let mut skip = false;
    let mut at = 0;
    while at < words.len() {
        let word = words[at].as_ref().as_encoded_bytes();
        if word == b"--" {
            skip = true;
        }
        if skip {
            kept.push(words[at].clone());
            at += 1;
            continue;
        }
        match word {
            b"-logfile" | b"-conffile" => {
                at += 1;
                if let Some(value) = words.get(at) {
                    kept.push(words[at - 1].clone());
                    kept.push(value.clone());
                    if word == b"-conffile" {
                        conffile = Some(value.as_ref().to_os_string());
                    }
                }
            }
            b"-dryrun" => dryrun = true,
            _ => {
                if let Some(request) = crate::ncbi_version::VersionRequest::of_word(word) {
                    return Preparsed::Version(request);
                }
                if let Some(name) = word.strip_prefix(b"-conffile=") {
                    // SAFETY: `name` is the rest of the word after an ASCII prefix, so it starts
                    // at a boundary of the platform's encoding and ends where the word ends.
                    let name = unsafe { OsStr::from_encoded_bytes_unchecked(name) };
                    conffile = Some(name.to_os_string());
                }
                kept.push(words[at].clone());
            }
        }
        at += 1;
    }
    Preparsed::Words {
        words: kept,
        conffile,
        dryrun,
    }
}

/// The file of the last `-conffile` of `words` (`preparse`): none before a version word.
fn conffile_of<S: AsRef<OsStr> + Clone>(words: &[S]) -> Option<OsString> {
    match preparse(words) {
        Preparsed::Words { conffile, .. } => conffile,
        Preparsed::Version(_) => None,
    }
}

/// What NCBI's application start writes about the registry files, and the settings it leaves
/// (`ApplicationStart::of_process`). NCBI reads the system-wide file (`.ncbirc`) twice: when the
/// program's usage report is made (`CBlastUsageReport`, a member of every BLAST+ application,
/// before `AppMain`; `usage_report`) and when the application loads its configuration
/// (`<program>.ini`, then `.ncbirc`, after the `-logfile`, `-version` and `-conffile` words and
/// before the program's arguments; `config`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:208-210
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (`x_TryInit`, which the
/// `-version` word ends before, ncbiapp.cpp:972-977)
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1650-1663 (`IncludeNcbircIfAllowed`)
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:669-674 (`x_TryMain`: an exception of
/// the configuration load ends the application with exit status 2)
/// Where NCBI's messages go (NCBI BLAST+ 2.17.0, the official Linux x64 release and bioconda's;
/// g_reg probes of 2026-11-11, env parity DECISIONS 13): the usage report's are written to
/// `cerr` at once (`usage_report`). An application with a `-logfile` file writes the messages
/// collected so far, the usage report's, into it when it opens it (`write_log_setup`), then
/// posts the configuration's messages to its diagnostics (`write_config_load`): `cerr`, the
/// file, or nowhere. A configuration whose `<program>.ini` throws ends the application with
/// exit status 2, after the exception's report and then every message collected (the usage
/// report's and the report again).
#[derive(Clone, Debug)]
pub struct ApplicationStart {
    /// The program whose start this is (the name of its `.ini`).
    program: String,
    /// The usage report's messages, or LOSAT's refusal of the run there (before `-version`).
    usage_report: Result<String, String>,
    /// The configuration load's messages, or LOSAT's refusal there.
    config: Result<ConfigLoad, String>,
    /// The settings after the start (`settings`).
    settings: Result<ApplicationSettings, String>,
}

/// The configuration load of an application start.
#[derive(Clone, Debug, PartialEq, Eq)]
enum ConfigLoad {
    /// The messages of a load (an empty text when there are none).
    Loaded(String),
    /// The report of the exception that ended the load (exit status 2).
    Failed(String),
    /// NCBI's post that the file of `-conffile` is not a file, after which the application
    /// aborts (exit status 255) before the load (`conf_file_path`).
    NotFound(String),
}

impl ApplicationStart {
    /// The application start of `program` in this process (its environment, its registry files
    /// and its facts, `ProcessFacts`): what the command line's entry writes before the program's
    /// arguments are read, and the settings that it leaves (`settings`). The entry makes one
    /// and gives it to what reads after it (the argument reader, `cli::try_parse_from_start`; the
    /// program's run), as NCBI's application starts once.
    ///
    /// `words` are the arguments after the program's name: the file of their last `-conffile`
    /// (`preparse`) is loaded in place of `<program>.ini`.
    pub fn of_process<S: AsRef<OsStr> + Clone>(program: &str, words: &[S]) -> ApplicationStart {
        application_start(
            program,
            conffile_of(words).as_deref(),
            std::env::vars_os(),
            &|name: &str| std::env::var_os(name),
            &ProcessFacts::of_this_process(),
        )
    }

    /// The application start of `program` in an environment of its own: the variables
    /// `variables` and no other, no home directory, no directory of the program, and the registry
    /// files of `files` (LOSAT Web, whose module has neither an environment nor a file system: its
    /// host registers the files, `web/adapter/src/files.rs`). Otherwise as `of_process`.
    pub fn with_files<S: AsRef<OsStr> + Clone>(
        program: &str,
        words: &[S],
        variables: &[(&str, &str)],
        files: &dyn RegistryFiles,
    ) -> ApplicationStart {
        let env = |name: &str| {
            variables
                .iter()
                .find(|(key, _)| *key == name)
                .map(|(_, value)| OsString::from(value))
        };
        application_start(
            program,
            conffile_of(words).as_deref(),
            variables
                .iter()
                .map(|(name, value)| (OsString::from(name), OsString::from(value))),
            &env,
            &ProcessFacts {
                // An empty home, which the search path leaves out (`registry_search_path`), as a
                // target without a user database has (`ncbi_home`).
                home_by_uid: &|| Some(OsString::new()),
                program_dirs: &|_| Ok(Vec::new()),
                files,
            },
        )
    }

    /// The program whose start this is.
    pub fn program(&self) -> &str {
        &self.program
    }

    /// The NCBI application settings that LOSAT reproduces after the start, or LOSAT's refusal
    /// of them (`settings_of` describes them).
    pub fn settings(&self) -> Result<ApplicationSettings, String> {
        self.settings.clone()
    }

    /// The usage report's messages, which go to `cerr` at once, or LOSAT's refusal there: a
    /// registry file that LOSAT does not read as NCBI does, a `[BLAST] BLAST_USAGE_REPORT` that
    /// is not a Boolean (NCBI aborts), or a refused environment variable when there are
    /// messages to write (LOSAT does not write them under it).
    pub fn usage_report(&self) -> Result<&str, &str> {
        self.usage_report.as_deref().map_err(String::as_str)
    }

    /// NCBI's first `SetupDiag` with a `-logfile` file: the file gets the messages collected
    /// so far, the usage report's (`CDiagContext::SetupDiag` with `eDCM_NoChange` flushes them
    /// to a new log; ncbiapp.cpp:868-871, ncbidiag.cpp:3863-3874).
    pub fn write_log_setup(&self, target: &mut DiagTarget) {
        if let Ok(text) = &self.usage_report {
            target.write_to_file(text.as_bytes());
        }
    }

    /// The configuration load's messages on `target` (`stderr` for `cerr`): `Ok(None)` after a
    /// load, `Ok(Some(2))` when the program's `.ini` (or the file of `-conffile`) threw (the
    /// application then exits with that status, having posted the report and then every
    /// collected message, the replay), `Ok(Some(255))` when the file of `-conffile` is not a
    /// file (NCBI's post, then its abort), `Err` with LOSAT's refusal. The replay goes to `target`, except in the programs whose
    /// application holds a `CBlastAppDiagHandler` (`replays_to_stderr`).
    pub fn write_config_load(
        &self,
        target: &mut DiagTarget,
        stderr: &mut dyn Write,
    ) -> Result<Option<i32>, String> {
        match &self.config {
            Err(reason) => Err(reason.clone()),
            Ok(ConfigLoad::Loaded(posts)) => {
                target.write_text(posts.as_bytes(), stderr);
                Ok(None)
            }
            Ok(ConfigLoad::NotFound(post)) => {
                target.write_text(post.as_bytes(), stderr);
                Ok(Some(255))
            }
            Ok(ConfigLoad::Failed(report)) => {
                target.write_text(report.as_bytes(), stderr);
                if !replays_collected(&self.program) {
                    return Ok(Some(2));
                }
                let replay = format!("{}{report}", self.usage_report.as_deref().unwrap_or(""));
                if replays_to_stderr(&self.program) {
                    DiagTarget::Stderr.write_text(replay.as_bytes(), stderr);
                } else {
                    target.write_text(replay.as_bytes(), stderr);
                }
                Ok(Some(2))
            }
        }
    }
}

/// Whether `program`'s application gives the messages it collected before its diagnostics were
/// set up to `cerr` when it ends after a failed configuration load, whatever its `-logfile`:
/// those whose application holds a `CBlastAppDiagHandler` as a member (blastn, blastp, blastx,
/// tblastn, psiblast). The handler takes the default handler (`cerr`) when the application is
/// made, before `AppMain`, and puts it back when the application ends, before the collected
/// messages are posted; tblastx, rpsblast and rpstblastn make theirs in `Run`, and makeblastdb and
/// blastdbcmd have none, so their messages go to the `-logfile` target. NCBI BLAST+ 2.17.0, the
/// official Linux x64 release (g_reg probes `probe_h`, 2026-11-12, a malformed `.ncbirc` and
/// `<program>.ini`): blastn, blastp, blastx, tblastn and psiblast with `-logfile log.txt` write
/// U U P on stderr and U P in the file, with `-logfile /dev/null` U U P on stderr; tblastx,
/// rpsblast, makeblastdb and blastdbcmd U on stderr and U P U P in the file.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:80 (likewise c++/src/app/blast/blastp_app.cpp:80,
/// blastx_app.cpp:80, tblastn_app.cpp:79, psiblast_app.cpp:122 `m_bah`)
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_aux.hpp:252
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1181-1187
/// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:94 (in `Run`; likewise
/// rpsblast_app.cpp:113)
// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:258-261
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:750-753 (`x_TryMain`, after a failed
// `x_TryInit` too)
/// Whether `program`'s application posts again, when it ends after a failed configuration load,
/// the messages it collected before its diagnostics were set up (the replay of
/// `write_config_load`): not segmasker, whose `Exit` takes the diagnostic stream away before
/// (NCBI BLAST+ 2.17.0 segmasker, the official release and the conda build, `-conffile` of a
/// malformed file: the report once, exit status 2; dustmasker writes it twice; fixtures
/// `conf_malformed`).
fn replays_collected(program: &str) -> bool {
    !matches!(program, "segmasker" | "convert2blastmask")
}

fn replays_to_stderr(program: &str) -> bool {
    matches!(
        program,
        "blastn" | "blastp" | "blastx" | "tblastn" | "psiblast"
    )
}

/// Whether the usage report's registry reads `.ncbirc`: not with `BLAST_USAGE_REPORT` set to a
/// false value (the report then makes no registry), nor with `.ncbirc` turned off in the
/// environment (`NCBI_DONT_USE_NCBIRC`, or the entry `[NCBI] DONT_USE_NCBIRC`, whose only layer
/// with entries there is the environment's: `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`, an empty
/// value too).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:195-206
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1634-1648 (`IncludeNcbircIfAllowed`,
/// cited at `settings_of`)
fn usage_report_reads_ncbirc(program: &str, env: &dyn Fn(&str) -> Option<OsString>) -> bool {
    has_usage_report(program)
        && usage_report_variable(env) != Some(Some(false))
        && !ncbirc_turned_off(env)
}

/// Whether `program`'s NCBI application has BLAST's usage report (`CBlastUsageReport`, a member
/// of the BLAST+ applications, which reads `BLAST_USAGE_REPORT` and `.ncbirc` before `AppMain`).
/// rpsbproc's has none (no `CBlastUsageReport` in RpsbProc 0.5.1's source): NCBI's rpsbproc
/// (the 2025-05-13 package and bioconda's 0.5.1) with a `.ncbirc` of a malformed line posts
/// nothing for `-version` and its `Critical` once for `-h` and for a run (env parity probe of
/// 2026-11-12 on ws), the configuration load's only. dustmasker's has none either: its
/// application is a `CNcbiApplication` with no member of BLAST's.
///
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.hpp:47-55
fn has_usage_report(program: &str) -> bool {
    !matches!(
        program,
        "rpsbproc" | "dustmasker" | "segmasker" | "windowmasker"
    )
}

/// `BLAST_USAGE_REPORT` as `NStr::StringToBool` reads it (`usage_report_reads_ncbirc`): `None`
/// when it is not set, `Some(None)` for a value that is not a Boolean (NCBI aborts).
fn usage_report_variable(env: &dyn Fn(&str) -> Option<OsString>) -> Option<Option<bool>> {
    env("BLAST_USAGE_REPORT").map(|value| string_to_bool(&value.to_string_lossy()))
}

/// Whether the environment turns `.ncbirc` off for every registry (`IncludeNcbircIfAllowed`,
/// cited at `settings_of`): `NCBI_DONT_USE_NCBIRC`, or
/// `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` (the environment's layer of `[NCBI] DONT_USE_NCBIRC`),
/// set with any value.
fn ncbirc_turned_off(env: &dyn Fn(&str) -> Option<OsString>) -> bool {
    env("NCBI_DONT_USE_NCBIRC").is_some() || env("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC").is_some()
}

/// The messages of `program`'s read of the system-wide file (`ncbi_registry::syntax_error_notice`,
/// `ncbi_registry::error_reading_post`).
fn system_wide_posts(read: &RegistryRead, program: &str) -> String {
    let mut posts = String::new();
    if read.stream_failed {
        posts.push_str(error_reading_post());
    }
    if let Some(error) = &read.error {
        posts.push_str(&syntax_error_notice(error, program));
    }
    posts
}

/// `[WINDOW_MASKER] WINDOW_MASKER_PATH` from the registry that blastn builds for
/// `-window_masker_taxid`, highest priority first: the environment layer
/// (`NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH`, `config`), the variable
/// `WINDOW_MASKER_PATH` (`plain`, mapped in at `ePriority_MaxUser`, below the environment
/// layer's `ePriority_Environment`), then `.ncbirc` (`ncbirc`). As for `[BLAST] BLASTDB`
/// (`blastdb_of`), a variable set to the empty string and an empty value in the file are
/// entries that hide the layers below; NCBI then uses the current directory (`Ok(None)` too).
/// The file is read before the entry is asked for (the registry's constructor loads it): its
/// posts are written also when a variable gives the directory, and a file that LOSAT does not
/// read as NCBI does (`window_masker_ncbirc`) stops `-window_masker_taxid` then too.
/// This registry is not the application's: it has no `.ini` file, so the program's `.ini`
/// neither sets the entry nor turns `.ncbirc` off (`tests/winmask/expected/taxid_path.tsv`,
/// NCBI 2.17.0 on 2026-10-19).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:374-391
/// NCBI reference (598d8ae6): c++/include/corelib/ncbireg.hpp:808-811,978-981
/// `WindowMaskerPathGet()` is empty in blastn (only `WindowMaskerPathInit` sets it, and the
/// applications do not call it).
fn window_masker_path_of(
    config: Option<Result<String, String>>,
    plain: Option<Result<String, String>>,
    ncbirc: Result<NcbircWindowMasker, String>,
) -> Result<WindowMaskerRegistry, String> {
    let NcbircWindowMasker {
        value,
        load_posts,
        posts_once,
    } = ncbirc?;
    let path = match config.or(plain) {
        Some(variable) => Some(variable?),
        None => value,
    };
    Ok(WindowMaskerRegistry {
        path,
        load_posts,
        posts_once,
    })
}

/// The settings of `program`'s application start without `-conffile` (`application_start`), for
/// the tests of the settings.
#[cfg(test)]
fn application_settings(
    program: &str,
    variables: impl Iterator<Item = (OsString, OsString)>,
    env: &dyn Fn(&str) -> Option<OsString>,
    process: &ProcessFacts,
) -> Result<ApplicationSettings, String> {
    application_start(program, None, variables, env, process).settings()
}

/// The application start of `program` (`ApplicationStart`) over the environment `variables`
/// (all of them) and `env` (one of them by name), and the facts of the process
/// (`ProcessFacts`). The search path is made once for both reads of `.ncbirc`
/// (`registry_search_path`), and the file found is read once: its two reads give the same.
///
/// The configuration load reads `.ncbirc` again unless `CMetaRegistry` keeps the usage report's
/// read: it does when that read threw nothing and the application's load asks for it with the
/// same flags, which it does when no `<program>.ini` was loaded (`LoadConfig` calls
/// `IncludeNcbircIfAllowed(fWithNcbirc)`, as the usage report's registry does; after a `.ini` it
/// is called with `fJustCore` added, ncbireg.cpp:1687-1692 and metareg.cpp:90). Only a read that
/// posts without throwing shows it: a one-byte `.ncbirc` of a byte order mark's first byte is
/// reported once, twice with a `blastn.ini` (NCBI BLAST+ 2.17.0, g_reg probes `b1full`,
/// `b1ini`); an exception is reported by both reads (`CMetaRegistry::x_Load` keeps no failed
/// read).
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:200-209 (`x_Load`: a read kept for
/// the same name, style and flags is not read again unless the file changed)
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1247-1252, 1278-1279
/// (`LoadConfig`: `.ncbirc` without a `.ini`)
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:210-211 (a
/// `[BLAST] BLAST_USAGE_REPORT` of the usage report's registry that is not a Boolean throws
/// in the constructor of the application, before `AppMain`: NCBI aborts)
/// The environment layer of `HasEntry("NCBI", "DONT_USE_NCBIRC")` (`IncludeNcbircIfAllowed`,
/// cited at `settings_of`) is `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`, set with
/// any value.
fn application_start(
    program: &str,
    conffile: Option<&OsStr>,
    mut variables: impl Iterator<Item = (OsString, OsString)>,
    env: &dyn Fn(&str) -> Option<OsString>,
    process: &ProcessFacts,
) -> ApplicationStart {
    let environment_refusal = variables.find_map(|(name, _)| rejected_variable(&name));
    let mut start = ApplicationStart {
        program: program.to_string(),
        usage_report: Ok(String::new()),
        config: Ok(ConfigLoad::Loaded(String::new())),
        settings: Err(String::new()),
    };
    // A start whose messages LOSAT would write under a variable that it refuses is refused
    // there.
    let unless_refused = |posts: String| match &environment_refusal {
        Some(reason) if !posts.is_empty() => Err(reason.clone()),
        _ => Ok(posts),
    };
    // The usage report reads `BLAST_USAGE_REPORT` first, before it makes its registry: a value
    // that is not a Boolean aborts NCBI, `-version` and `-h` included (NCBI BLAST+ 2.17.0 with
    // `abc` or an empty value: rc 134, g_reg probes `m2.ver`, `m2.h`, `m2.empty`).
    if has_usage_report(program) && usage_report_variable(env) == Some(None) {
        start.usage_report = Err(format!(
            "the environment variable BLAST_USAGE_REPORT has the value {:?}, which is not a Boolean of NCBI BLAST+ (it aborts); this is not supported by LOSAT",
            env("BLAST_USAGE_REPORT").unwrap_or_default().to_string_lossy()
        ));
        return ended(start);
    }
    let search_path = match registry_search_path(program, env, process) {
        Ok(path) => path,
        Err(reason) => {
            start.settings = Err(environment_refusal.clone().unwrap_or(reason));
            // NCBI's check of the `-conffile` file comes before its load and needs no search
            // path (`conf_file_path`).
            if let Some(Err(ConfFile::NotFound(post))) =
                conffile.map(|name| conf_file_path(name, &[], process.files))
            {
                start.config = Ok(ConfigLoad::NotFound(post));
            }
            return start;
        }
    };
    let dont_use_ncbirc = ncbirc_turned_off(env);
    let ncbirc = if dont_use_ncbirc {
        None
    } else {
        find_registry(&search_path, SYSTEM_WIDE_REGISTRY, process.files).map(|path| {
            (
                read_registry_file(&path, Layer::SystemWide, process.files),
                path,
            )
        })
    };
    // The usage report.
    let usage_report_reads = usage_report_reads_ncbirc(program, env);
    let mut usage_read = None;
    if usage_report_reads {
        match &ncbirc {
            Some((Err(reason), _)) => {
                start.usage_report = Err(reason.clone());
                return ended(start);
            }
            Some((Ok(Some(read)), path)) => usage_read = Some((read, path.as_path())),
            Some((Ok(None), _)) | None => {}
        }
    }
    if has_usage_report(program) && usage_report_variable(env) != Some(Some(false)) {
        if let Some(reason) = usage_report_entry_refusal(env, usage_read) {
            start.usage_report = Err(reason);
            return ended(start);
        }
    }
    if let Some((read, _)) = usage_read {
        start.usage_report = unless_refused(system_wide_posts(read, program));
        if start.usage_report.is_err() {
            return ended(start);
        }
    }
    // The configuration: `<program>.ini` (or the file of `-conffile`), then `.ncbirc`.
    let ini_path = match conffile {
        None => find_registry(&search_path, &format!("{program}.ini"), process.files),
        Some(name) => match conf_file_path(name, &search_path, process.files) {
            Ok(path) => path,
            Err(ConfFile::NotFound(post)) => {
                start.config = Ok(ConfigLoad::NotFound(post));
                return ended(start);
            }
            Err(ConfFile::CannotOpen(report)) => {
                start.config = match &environment_refusal {
                    Some(reason) => Err(reason.clone()),
                    None => Ok(ConfigLoad::Failed(report)),
                };
                return ended(start);
            }
            Err(ConfFile::Refused(reason)) => {
                start.config = Err(reason);
                return ended(start);
            }
        },
    };
    let ini = match ini_path {
        None => None,
        Some(path) => match read_registry_file(&path, Layer::Application, process.files) {
            Err(reason) => {
                start.config = Err(reason);
                return ended(start);
            }
            Ok(read) => read.map(|read| (path, read)),
        },
    };
    let mut posts = String::new();
    if let Some((_, read)) = &ini {
        if let Some(error) = &read.error {
            let report = error.report_default_ex(106, 15, "Application's initialization failed");
            start.config = match &environment_refusal {
                Some(reason) => Err(reason.clone()),
                None => Ok(ConfigLoad::Failed(registry_text(report))),
            };
            return ended(start);
        }
        if read.stream_failed {
            posts.push_str(error_reading_post());
        }
    }
    let ini_turns_it_off = ini
        .as_ref()
        .is_some_and(|(_, read)| read.last_value("NCBI", "DONT_USE_NCBIRC").is_some());
    let use_ncbirc = !dont_use_ncbirc && !ini_turns_it_off;
    let mut ncbirc_entries = None;
    if use_ncbirc {
        match &ncbirc {
            Some((Err(reason), _)) => {
                start.config = Err(reason.clone());
                return ended(start);
            }
            Some((Ok(Some(read)), path)) => {
                let kept = usage_report_reads && read.error.is_none() && ini.is_none();
                if !kept {
                    posts.push_str(&system_wide_posts(read, program));
                }
                ncbirc_entries = Some((path.as_path(), read.entries.as_slice()));
            }
            Some((Ok(None), _)) | None => {}
        }
    }
    start.config = unless_refused(posts).map(ConfigLoad::Loaded);
    let ini_entries = ini
        .as_ref()
        .map(|(path, read)| (path.as_path(), read.entries.as_slice()));
    // The registry of `-window_masker_taxid` reads `.ncbirc` also when the program's `.ini`
    // turned it off for the application (`window_masker_path_of`). Its flags are those of the
    // usage report's registry, and of the application's load without a `.ini`: a read of
    // either is kept for it.
    let window_masker_kept = usage_report_reads || (use_ncbirc && ini.is_none());
    let ncbirc_window_masker = window_masker_ncbirc(ncbirc.as_ref(), program, window_masker_kept);
    start.settings = match &environment_refusal {
        Some(reason) => Err(reason.clone()),
        None => settings_of(env, ini_entries, ncbirc_entries, ncbirc_window_masker),
    };
    start
}

/// LOSAT's refusal of the usage report's `[BLAST] BLAST_USAGE_REPORT` when it is not a Boolean
/// (NCBI aborts, before `AppMain`): the value of the environment's layer
/// (`NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT`), else that of `.ncbirc` when the report reads it
/// (`usage_read`). NCBI BLAST+ 2.17.0 (g_reg probes of 2026-11-12): `maybe` in the variable
/// aborts (rc 134) also with `NCBI_DONT_USE_NCBIRC`, not with `BLAST_USAGE_REPORT=0`; in
/// `.ncbirc` it aborts also with `[NCBI] DONT_USE_NCBIRC` in `blastn.ini`, not with
/// `NCBI_DONT_USE_NCBIRC`, nor under a variable `yes` (probes `l1.env`, `l1.env_dont`,
/// `l1.env_ur0`, `l1.rc_inidont`, `l1.rc_dont`, `l1.env_ok`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:208-211
fn usage_report_entry_refusal(
    env: &dyn Fn(&str) -> Option<OsString>,
    usage_read: Option<(&RegistryRead, &Path)>,
) -> Option<String> {
    let not_a_boolean = |value: &str| string_to_bool(value).is_none();
    if let Some(value) = env("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT") {
        let value = value.to_string_lossy();
        return not_a_boolean(&value).then(|| {
            format!(
                "the environment variable NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT has the value {value:?}, which is not a Boolean of NCBI BLAST+ (it aborts); this is not supported by LOSAT"
            )
        });
    }
    let (read, path) = usage_read?;
    let value = String::from_utf8_lossy(read.last_value("BLAST", "BLAST_USAGE_REPORT")?);
    not_a_boolean(&value).then(|| {
        format!(
            "the registry file {} sets [BLAST] BLAST_USAGE_REPORT to {value:?}, which is not a Boolean of NCBI BLAST+ (it aborts); this is not supported by LOSAT",
            path.display()
        )
    })
}

/// A start that a stage refused or ended (`application_start`): its settings are that refusal,
/// or the end (which the command line's entry and the library's calls reach first, at
/// `ApplicationStart::write_config_load`).
fn ended(mut start: ApplicationStart) -> ApplicationStart {
    let reason = match (&start.usage_report, &start.config) {
        (Err(reason), _) | (_, Err(reason)) => reason.clone(),
        (_, Ok(ConfigLoad::NotFound(post))) => format!(
            "{}; NCBI BLAST+ aborts there (exit status 255), and this is not supported by LOSAT",
            post.trim_end()
        ),
        _ => "a registry file that NCBI BLAST+ reports as malformed ends the program (exit status 2), and this is not supported by LOSAT".to_string(),
    };
    start.settings = Err(reason);
    start
}

/// Rejects the NCBI application settings that change `program`'s output: the environment
/// variables of the NCBI C++ Toolkit and the entries of the registry files that NCBI loads
/// (`<program>.ini`, then `.ncbirc` unless it is turned off). Returns the settings that
/// LOSAT reproduces: whether the data loaders are used (`[BLAST] DATA_LOADERS`, which
/// decides whether a record's first line can be a sequence identifier), the registry
/// entry `[BLAST] BLASTDB` of the database search path, and the directory of
/// `-window_masker_taxid` (`window_masker_path_of`).
///
/// An application start has these (`ApplicationStart::settings`) after it has written what
/// NCBI writes about the files: the command line's entry, and each call of LOSAT's library that
/// reads NCBI's settings (`api::start`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1254-1259
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1636-1648
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:299-303
fn settings_of(
    env: &dyn Fn(&str) -> Option<OsString>,
    ini: Option<(&Path, &[RegistryEntry])>,
    ncbirc: Option<(&Path, &[RegistryEntry])>,
    ncbirc_window_masker: Result<NcbircWindowMasker, String>,
) -> Result<ApplicationSettings, String> {
    for (path, entries) in ini.iter().chain(ncbirc.iter()) {
        check_registry_entries(path, entries)?;
    }
    let text = |file: Option<(&Path, &[RegistryEntry])>, section: &str, name: &str| {
        file.map_or(Ok(None), |(path, entries)| {
            file_text_entry(path, entries, section, name)
        })
    };
    let mut max_seqid_length = env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH")
        .map(|value| Ok(value.to_string_lossy().into_owned()));
    for (path, entries) in ini.iter().chain(ncbirc.iter()) {
        if max_seqid_length.is_none() {
            max_seqid_length = max_seqid_length_entry(path, entries);
        }
    }
    let utf8_variable = |name: &str| {
        env(name).map(|value| {
            value.into_string().map_err(|_| {
                format!("the environment variable {name} is not valid UTF-8; this is not supported by LOSAT")
            })
        })
    };
    let env_data_loaders = env("NCBI_CONFIG__BLAST__DATA_LOADERS");
    let env_blastdb = utf8_variable("NCBI_CONFIG__BLAST__BLASTDB");
    let window_masker_registry = window_masker_path_of(
        utf8_variable("NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH"),
        utf8_variable("WINDOW_MASKER_PATH"),
        ncbirc_window_masker,
    );
    Ok(ApplicationSettings {
        data_loaders: data_loaders_of(
            env_data_loaders.as_deref().map(OsStr::as_encoded_bytes),
            [ini, ncbirc].map(|file| {
                file.and_then(|(_, entries)| last_value(entries, "BLAST", "DATA_LOADERS"))
            }),
        ),
        blastdb: blastdb_of(
            env_blastdb,
            [
                text(ini, "BLAST", "BLASTDB"),
                text(ncbirc, "BLAST", "BLASTDB"),
            ],
        ),
        seq_id_rank_variable: SEQ_ID_RANK_VARIABLES
            .into_iter()
            .find(|name| env(name).is_some()),
        max_seqid_length: max_seqid_length.transpose(),
        matrix_search: MatrixSearch {
            data_path: env("NCBI_DATA_PATH"),
            data: match env("NCBI_CONFIG__NCBI__DATA").or_else(|| env("NCBI_CONFIG__NCBI__Data")) {
                Some(value) => Ok(Some(value)),
                None => blastdb_of(
                    None,
                    [text(ini, "NCBI", "Data"), text(ncbirc, "NCBI", "Data")],
                )
                .map(|data| data.map(OsString::from)),
            },
            blastmat: env("BLASTMAT"),
        },
        window_masker_registry,
    })
}

#[cfg(test)]
#[path = "ncbi_environment_tests.rs"]
mod tests;
