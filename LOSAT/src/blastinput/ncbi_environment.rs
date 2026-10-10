//! The configuration that NCBI's application layer (`CNcbiApplication`) reads before a BLAST
//! program runs: environment variables of the NCBI C++ Toolkit and the registry files
//! (`.ncbirc` and the program's `.ini`). LOSAT does not reproduce the settings that change
//! a program's output or exit status, and rejects them explicitly (plan DW-13); a registry
//! file with only entries that change no output (such as `[BLAST] BLASTDB`) is accepted.

use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

use crate::report::outfmt6::truncate_spaces;

/// Registry entries (section, name; case-insensitive) that change no output of a search with
/// `-subject`: the database and data-loader settings and the data directory, and
/// `[NCBI] DONT_USE_NCBIRC`, which only decides whether `.ncbirc` is read.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:1436-1439
/// ```c
///     if (app) {
///         const CNcbiRegistry& registry = app->GetConfig();
///         if (registry.HasEntry("BLAST", "BLASTDB"))
///             CDirEntry::NormalizePath(registry.Get("BLAST", "BLASTDB"), eFollowLinks);
/// ```
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
];

/// Environment variables, besides the `DIAG_` and `NCBI_CONFIG__` families, that the NCBI
/// C++ Toolkit reads for every application and that change its error output or exit status.
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbiapp.cpp:1031-1044
/// ```c
///     // Setup some debugging features from environment variables.
///     if ( !m_Environ->Get(DIAG_TRACE).empty() ) {
///         SetDiagTrace(eDT_Enable, eDT_Enable);
///     }
///     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
///     if ( !post_level.empty() ) {
///         EDiagSev sev;
///         if (CNcbiDiag::StrToSeverityLevel(post_level.c_str(), sev)) {
///             SetDiagFixedPostLevel(sev);
///         }
///     }
///     if ( !m_Environ->Get(ABORT_ON_THROW).empty() ) {
///         SetThrowTraceAbort(true);
///     }
/// ```
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
];

fn is_harmless_entry(section: &str, name: &str, value: &str) -> bool {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_usage_report.cpp:209-211
    // ```c
    // 	CRef<CNcbiRegistry> registry(new CNcbiRegistry(empty_stream, IRegistry::fWithNcbirc));
    // 	if (registry->HasEntry("BLAST", "BLAST_USAGE_REPORT")) {
    // 		bool enable = NStr::StringToBool(registry->Get("BLAST", "BLAST_USAGE_REPORT"));
    // ```
    // A Boolean only decides whether the usage is reported; another value aborts NCBI.
    if section.eq_ignore_ascii_case("BLAST") && name.eq_ignore_ascii_case("BLAST_USAGE_REPORT") {
        return ncbi_string_to_bool(value).is_some();
    }
    HARMLESS_ENTRIES
        .iter()
        .any(|(s, n)| section.eq_ignore_ascii_case(s) && name.eq_ignore_ascii_case(n))
}

/// `NStr::StringToBool`.
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbistr.cpp:2817-2838
/// ```c
/// bool NStr::StringToBool(const CTempString str)
/// {
///     if ( str == "1"  ||
///          AStrEquiv(str, s_kTrueString,  PNocase())  ||
///          AStrEquiv(str, s_kTString,     PNocase())  ||
///          AStrEquiv(str, s_kYesString,   PNocase())  ||
///          AStrEquiv(str, s_kYString,     PNocase())  ||
///          AStrEquiv(str, s_kOnString,    PNocase()) ) {
///         errno = 0;
///         return true;
///     }
///     if ( str == "0"  ||
///          AStrEquiv(str, s_kFalseString, PNocase())  ||
///          AStrEquiv(str, s_kFString,     PNocase())  ||
///          AStrEquiv(str, s_kNoString,    PNocase())  ||
///          AStrEquiv(str, s_kNString,     PNocase())  ||
///          AStrEquiv(str, s_kOffString,   PNocase()) ) {
///         errno = 0;
///         return false;
///     }
///     NCBI_THROW2(CStringException, eConvert,
///                 "String cannot be converted to bool", 0);
/// }
/// ```
pub fn ncbi_string_to_bool(value: &str) -> Option<bool> {
    let is = |word: &str| value.eq_ignore_ascii_case(word);
    if value == "1" || is("true") || is("t") || is("yes") || is("y") || is("on") {
        Some(true)
    } else if value == "0" || is("false") || is("f") || is("no") || is("n") || is("off") {
        Some(false)
    } else {
        None
    }
}

/// Why LOSAT rejects an environment variable, or `None`.
///
/// NCBI reference: ncbi-blast/c++/src/corelib/env_reg.cpp:418-469
/// ```c
/// bool CNcbiEnvRegMapper::EnvToReg(const string& env_in, string& section,
///                                  string& name) const
/// {
///     if (env_in.size() <= kPrefixLen  ||  !NStr::StartsWith(env_in, kPrefix) ) {
///         return false;
///     }
///     ...
///     SIZE_TYPE uu_pos = env.find("__", section_start_pos + 1);
///     if (uu_pos == NPOS  ||  uu_pos == env.size() - 2) {
///         return false;
///     }
///     /* Parse section and entry names from the variable */
///     if (env[kPrefixLen] == '_') { // regular entry
///         section = env.substr(kPrefixLen + 1, uu_pos - kPrefixLen - 1);
///         name    = env.substr(uu_pos + 2);
///     } else {
///         name    = env.substr(kPrefixLen - 1, uu_pos - kPrefixLen + 1);
///         _ASSERT(name[0] == '_');
///         name[0] = '.';
///         section = env.substr(uu_pos + 2);
///     }
/// ```
/// `kPrefix` is "NCBI_CONFIG_" (env_reg.cpp:362): `NCBI_CONFIG__<SECTION>__<NAME>` sets a
/// registry entry, `NCBI_CONFIG_<NAME>__<SECTION>` a special one, and
/// `NCBI_CONFIG_OVERRIDES` names a file of entries (ncbireg.cpp:1585-1602), while
/// `NCBI_CONFIG_PATH` is the search path of the files (`registry_search_path`). The
/// `DIAG_` variables are the parameters of the diagnostics (`DIAG_POST_LEVEL`,
/// `DIAG_OLD_POST_FORMAT`, `DIAG_TEE_TO_STDERR`, ...; ncbidiag.cpp), each of which changes
/// stderr or, with a value its parser rejects, stops NCBI.
fn rejected_variable(name: &OsStr, value: &OsStr) -> Option<String> {
    let name = name.to_str()?;
    if name == "NCBI_CONFIG_PATH" {
        return None;
    }
    if let Some(entry) = name.strip_prefix("NCBI_CONFIG__") {
        if let Some((section, key)) = entry.split_once("__") {
            if !section.is_empty()
                && !key.is_empty()
                && is_harmless_entry(section, key, &value.to_string_lossy())
            {
                return None;
            }
        }
        return Some(format!(
            "the environment variable {name}, which sets an entry of NCBI BLAST+'s registry that LOSAT does not know to change no output (it accepts only the entries listed for registry files), is not supported by LOSAT"
        ));
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1585-1586
    // ```c++
    //     const TXChar* xoverride_path = NcbiSys_getenv(_TX("NCBI_CONFIG_OVERRIDES"));
    //     if (xoverride_path  &&  *xoverride_path) {
    // ```
    // A set but empty NCBI_CONFIG_OVERRIDES names no file and changes nothing (re-audit round 2
    // finding R-7 of session SFd).
    if name == "NCBI_CONFIG_OVERRIDES" && value.is_empty() {
        return None;
    }
    if name.starts_with("NCBI_CONFIG_") {
        return Some(format!(
            "the environment variable {name}, which sets entries of NCBI BLAST+'s registry, is not supported by LOSAT"
        ));
    }
    if name.starts_with("DIAG_") || REJECTED_VARIABLES.contains(&name) {
        return Some(format!(
            "the environment variable {name}, a parameter of NCBI BLAST+'s diagnostics or error handling (LOSAT rejects the whole family: some members change NCBI's messages or exit status, some only on errors), is not supported by LOSAT"
        ));
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_usage_report.cpp:197-199
    // ```c
    // 	char * blast_usage_env = getenv("BLAST_USAGE_REPORT");
    // 	if(blast_usage_env != NULL){
    // 		bool enable = NStr::StringToBool(blast_usage_env);
    // ```
    if name == "BLAST_USAGE_REPORT" && ncbi_string_to_bool(&value.to_string_lossy()).is_none() {
        return Some(format!(
            "the environment variable BLAST_USAGE_REPORT has the value {:?}, which is not a Boolean of NCBI BLAST+ (it aborts); this is not supported by LOSAT",
            value.to_string_lossy()
        ));
    }
    None
}

/// The bytes of an `OsString` built from `OsStr::as_encoded_bytes` parts joined at ASCII
/// bytes (path separators): the bytes as they are on Unix and WASI, and a lossy conversion
/// elsewhere (Windows paths are not cut or joined by NCBI's Unix rules below).
fn os_string_from_bytes(bytes: Vec<u8>) -> OsString {
    #[cfg(unix)]
    {
        std::os::unix::ffi::OsStringExt::from_vec(bytes)
    }
    #[cfg(target_os = "wasi")]
    {
        std::os::wasi::ffi::OsStringExt::from_vec(bytes)
    }
    #[cfg(not(any(unix, target_os = "wasi")))]
    {
        OsString::from(String::from_utf8_lossy(&bytes).into_owned())
    }
}

fn path_from_bytes(bytes: &[u8]) -> PathBuf {
    PathBuf::from(os_string_from_bytes(bytes.to_vec()))
}

/// What NCBI's `CNcbiApplication` knows of the running program and its user, read by
/// `registry_search_path`.
struct ProgramContext {
    /// `argv[0]` as the program was started.
    argv0: Option<OsString>,
    /// The executable with its links resolved (`std::env::current_exe`, on Linux
    /// `readlink("/proc/self/exe")`, as NCBI's `/proc/<pid>/exe`).
    current_exe: Option<PathBuf>,
    /// The working directory (`CDir::GetCwd`, `getcwd`).
    cwd: Option<PathBuf>,
    /// The home directory of the user's passwd entry (`getpwuid(getuid())->pw_dir`), read only
    /// when `HOME` is not set; `None` when the user has no entry.
    passwd_home: Option<OsString>,
}

impl ProgramContext {
    fn of_this_process() -> Self {
        Self {
            argv0: std::env::args_os().next(),
            current_exe: std::env::current_exe().ok(),
            cwd: std::env::current_dir().ok(),
            // `std::env::home_dir` returns `HOME` when it is set and not empty, and otherwise
            // the passwd entry's directory (`getpwuid_r(getuid())`); it is asked only when
            // `HOME` is not set at all (`home_directory`).
            passwd_home: if std::env::var_os("HOME").is_none() {
                std::env::home_dir().map(PathBuf::into_os_string)
            } else {
                None
            },
        }
    }
}

/// A directory of NCBI's registry search path.
#[derive(Clone, Debug, PartialEq, Eq)]
enum SearchDir {
    Path(PathBuf),
    /// The home directory when `HOME` is not set and the user has no passwd entry: NCBI then
    /// looks the login name up (`USER`, `LOGNAME`, `getlogin()`, then `getpwnam`), which LOSAT
    /// does not port (`AUTHORITY.md` §J-8 of `docs/evidence/losat_web_e2h/`).
    UnknownHome,
}

/// `CDir::GetHome`: `HOME` when it is set (an empty value gives no directory), otherwise the
/// passwd entry's directory; on Windows `APPDATA`, then `USERPROFILE` (re-audit round 2
/// finding R-2 of session SFd).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3586-3597
/// ```c++
/// static bool s_GetHomeByUID(string& home)
/// {
///     // Get the info using user ID
///     struct passwd* pwd;
///
///     if ((pwd = getpwuid(getuid())) == 0) {
///         LOG_ERROR_ERRNO(48, "s_GetHomeByUID(): getpwuid() failed");
///         return false;
///     }
///     home = pwd->pw_dir;
///     return true;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3624-3657
/// ```c++
/// string CDir::GetHome(void)
/// {
///     string home;
///
/// #if defined(NCBI_OS_MSWIN)
///     // Get home dir from environment variables
///     // like - C:\Documents and Settings\user\Application Data
///     const TXChar* str = NcbiSys_getenv(_TX("APPDATA"));
///     if ( str ) {
///         home = _T_CSTRING(str);
///     } else {
///         // like - C:\Documents and Settings\user
///         str = NcbiSys_getenv(_TX("USERPROFILE"));
///         if ( str ) {
///             home = _T_CSTRING(str);
///         }
///     }
/// #elif defined(NCBI_OS_UNIX)
///     // Try get home dir from environment variable
///     char* str = NcbiSys_getenv(_TX("HOME"));
///     if ( str ) {
///         home = str;
///     } else {
///         // Try to retrieve the home dir -- first use user's ID,
///         // and if failed, then use user's login name.
///         if ( !s_GetHomeByUID(home) ) {
///             s_GetHomeByLOGIN(home);
///         }
///     }
/// #endif
///
///     // Add trailing separator if needed
///     return AddTrailingPathSeparator(home);
/// }
/// ```
/// The messages of `LOG_ERROR_ERRNO` are written only with `[NCBI] FileAPILogging`
/// (ncbifile.cpp:178-180, default false; its variable `NCBI_CONFIG__FILEAPILOGGING` is
/// rejected).
fn home_directory(env: &dyn Fn(&str) -> Option<OsString>, context: &ProgramContext) -> SearchDir {
    let home = if cfg!(windows) {
        env("APPDATA").or_else(|| env("USERPROFILE"))
    } else {
        match env("HOME") {
            Some(home) => Some(home),
            None => match &context.passwd_home {
                Some(home) => Some(home.clone()),
                None => return SearchDir::UnknownHome,
            },
        }
    };
    SearchDir::Path(PathBuf::from(home.unwrap_or_default()))
}

/// `CNcbiArguments::GetProgramDirname`: the name up to its last `/`, `\` or `:`, that byte
/// included, or nothing.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbienv.cpp:408-415
/// ```c++
/// string CNcbiArguments::GetProgramDirname(EFollowLinks follow_links) const
/// {
///     const string& name = GetProgramName(follow_links);
///     SIZE_TYPE base_pos = name.find_last_of("/\\:");
///     if (base_pos == NPOS)
///         return NcbiEmptyString;
///     return name.substr(0, base_pos + 1);
/// }
/// ```
fn program_dirname(name: &[u8]) -> &[u8] {
    match name
        .iter()
        .rposition(|&byte| matches!(byte, b'/' | b'\\' | b':'))
    {
        Some(pos) => &name[..=pos],
        None => &[],
    }
}

/// `CDirEntry::NormalizePath` on Unix: `.`, `..` and empty components removed by the text of
/// the path, and with `follow_links` every component that is a symbolic link replaced by its
/// target as it is reached.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-971
/// ```c++
///     if ( path.empty() ) {
///         return path;
///     }
///     ...
///     current = DeleteTrailingPathSeparator(path);
///     if ( current.empty() ) {
///         // root dir
///         return string(1, DIR_SEPARATOR);
///     }
///     while ( !current.empty()  ||  !tail.empty() ) {
///         std::list<string> pretail;
///         if ( !current.empty() ) {
///             NStr::Split(current, kSep, pretail);
///             current.erase();
///             if (pretail.front().empty()
///                 ) {
///                 // Absolute path
///                 head.clear();
///             }
///             tail.splice(tail.begin(), pretail);
///         }
///         string next;
///         if (!tail.empty()) {
///             next = tail.front();
///             tail.pop_front();
///         }
///         if ( !head.empty() ) { // empty heads should accept anything
///             string& last = head.back();
///             if (last == DIR_CURRENT) {
///                 if (!next.empty()) {
///                     head.pop_back();
///                 }
///             } else if (next == DIR_CURRENT) {
///                 // Leave out, since we already have content
///                 continue;
///             } else if (next.empty()) {
///                 continue; // leave out empty components in most cases
///             } else if (next == DIR_PARENT) {
///                 // Back up if possible, assuming existing path to be "physical"
///                 if (last.empty()) {
///                     // Already at the root; .. is a no-op
///                     continue;
///                 } else if (last != DIR_PARENT) {
///                     head.pop_back();
///                     continue;
///                 }
///             }
///         }
/// #ifdef NCBI_OS_UNIX
///         // Is there a Windows equivalent for readlink?
///         if ( follow_links ) {
///             string s(head.empty() ? next : NStr::Join(head, string(1, DIR_SEPARATOR)) + DIR_SEPARATOR + next);
///             char buf[PATH_MAX];
///             int  length = (int)readlink(s.c_str(), buf, sizeof(buf));
///             if (length > 0) {
///                 current.assign(buf, length);
///                 if (++link_depth >= 1024) {
///                     ...
///                     follow_links = eIgnoreLinks;
///                 }
///                 continue;
///             }
///         }
/// #endif
///         // Normal case: just append the next element to head
///         head.push_back(next);
///     }
///
///     // Special cases
///     if ( (head.size() == 0)  ||
///          (head.size() == 2  &&  head.front() == DIR_CURRENT  &&  head.back().empty()) ) {
///         // current dir
///         return DIR_CURRENT;
///     }
///     if (head.size() == 1  &&  head.front().empty()) {
///         // root dir
///         return string(1, DIR_SEPARATOR);
///     }
///     ...
///     // Compose path
///     return NStr::Join(head, string(1, DIR_SEPARATOR));
/// ```
/// The symlink depth warning (1024 links) is not reproduced.
fn normalize_path(path: &[u8], follow_links: bool) -> Vec<u8> {
    if path.is_empty() {
        return Vec::new();
    }
    // `DeleteTrailingPathSeparator` (ncbifile.cpp:465-472).
    let end = path
        .iter()
        .rposition(|&byte| byte != b'/')
        .map_or(0, |pos| pos + 1);
    let mut current = path[..end].to_vec();
    if current.is_empty() {
        return b"/".to_vec();
    }
    let mut follow_links = follow_links;
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
            if last.as_slice() == b"." {
                if !next.is_empty() {
                    head.pop();
                }
            } else if next.as_slice() == b"." || next.is_empty() {
                continue;
            } else if next.as_slice() == b".." {
                if last.is_empty() {
                    continue;
                } else if last.as_slice() != b".." {
                    head.pop();
                    continue;
                }
            }
        }
        if follow_links {
            let mut link = head.join(&b'/');
            if !head.is_empty() {
                link.push(b'/');
            }
            link.extend_from_slice(&next);
            if let Ok(target) = std::fs::read_link(path_from_bytes(&link)) {
                let target = target.as_os_str().as_encoded_bytes().to_vec();
                if !target.is_empty() {
                    current = target;
                    link_depth += 1;
                    if link_depth >= 1024 {
                        follow_links = false;
                    }
                    continue;
                }
            }
        }
        head.push(next);
    }
    if head.is_empty() || (head.len() == 2 && head[0].as_slice() == b"." && head[1].is_empty()) {
        return b".".to_vec();
    }
    if head.len() == 1 && head[0].is_empty() {
        return b"/".to_vec();
    }
    head.join(&b'/')
}

/// `CNcbiApplicationAPI::FindProgramExecutablePath` outside Windows, as `AppMain` calls it
/// (the program's name, as `CNcbiArguments` keeps it): `argv[0]` made absolute from the
/// working directory when it names a file there, otherwise searched in `PATH` by its base
/// name, then normalized without following links.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1451-1456
/// ```c++
///     if (argc > 0  &&  argv[0] != NULL  &&  argv[0][0] != '\0') {
///         ret_val = argv[0];
///     } else if (instance) {
///         ret_val = instance->GetArguments().GetProgramName();
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1544-1597
/// ```c++
///     string app_path = ret_val;
///
///     if ( !CDirEntry::IsAbsolutePath(app_path) ) {
///         ...
///         if ( CFile(app_path).Exists() ) {
///             // Relative path from the the current directory
///             app_path = CDir::GetCwd() + CDirEntry::GetPathSeparator() + app_path;
///             if ( !CFile(app_path).Exists() ) {
///                 app_path = kEmptyStr;
///             }
///         } else {
///             // Running from some path from PATH environment variable.
///             // Try to determine that path.
///             string env_path;
///             if (instance) {
///                 env_path = instance->GetEnvironment().Get("PATH");
///             } else {
///                 env_path = _T_STDSTRING(NcbiSys_getenv(_TX("PATH")));
///             }
///             list<string> split_path;
///             ...
///             NStr::Split(env_path, ":", split_path,
///                 NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
///             ...
///             string base_name = CDirEntry(app_path).GetBase();
///             ITERATE(list<string>, it, split_path) {
///                 app_path = CDirEntry::MakePath(*it, base_name);
///                 if ( CFile(app_path).Exists() ) {
///                     break;
///                 }
///                 app_path = kEmptyStr;
///             }
///         }
///     }
///     ret_val = CDirEntry::NormalizePath(
///         (app_path.empty() && argv != NULL && argv[0] != NULL) ? argv[0] : app_path);
/// ```
/// `CFile::Exists` is `IsFile` (ncbifile.hpp:4039-4042: `stat`, links followed, a regular
/// file). `GetBase` is the file name without its last extension (`SplitPath`,
/// ncbifile.cpp:358-377). An empty `argv[0]` is looked up as `ncbi`, the name
/// `CNcbiArguments` gives before its arguments are set (ncbienv.cpp:388-392); NCBI's warning
/// for a name it cannot find then (ncbiapp.cpp:901-908) is not reproduced.
fn find_program_executable_path(argv0: &[u8], cwd: Option<&Path>, path: Option<&OsStr>) -> Vec<u8> {
    let is_file = |bytes: &[u8]| path_from_bytes(bytes).is_file();
    let name: &[u8] = if argv0.is_empty() { b"ncbi" } else { argv0 };
    let mut app_path = name.to_vec();
    if !app_path.starts_with(b"/") {
        if is_file(&app_path) {
            let mut absolute = cwd
                .map(|cwd| cwd.as_os_str().as_encoded_bytes().to_vec())
                .unwrap_or_default();
            absolute.push(b'/');
            absolute.extend_from_slice(&app_path);
            app_path = if is_file(&absolute) {
                absolute
            } else {
                Vec::new()
            };
        } else {
            // `CDirEntry(app_path)` drops trailing separators (ncbifile.cpp:298-313).
            let entry = match app_path.iter().rposition(|&byte| byte != b'/') {
                Some(pos) if app_path.len() > 1 => &app_path[..=pos],
                _ => &app_path[..],
            };
            let file_name = match entry.iter().rposition(|&byte| byte == b'/') {
                Some(pos) => &entry[pos + 1..],
                None => entry,
            };
            let base_name = match file_name.iter().rposition(|&byte| byte == b'.') {
                Some(pos) => &file_name[..pos],
                None => file_name,
            }
            .to_vec();
            app_path = Vec::new();
            let path = path.map(OsStr::as_encoded_bytes).unwrap_or_default();
            for dir in path
                .split(|&byte| byte == b':')
                .filter(|dir| !dir.is_empty())
            {
                // `MakePath`: the directory, a separator unless it ends with one, the name.
                let mut candidate = dir.to_vec();
                if !candidate.ends_with(b"/") {
                    candidate.push(b'/');
                }
                candidate.extend_from_slice(&base_name);
                if is_file(&candidate) {
                    app_path = candidate;
                    break;
                }
            }
        }
    }
    let chosen: &[u8] = if app_path.is_empty() {
        argv0
    } else {
        &app_path
    };
    normalize_path(chosen, false)
}

/// The program's directories in NCBI's registry search path (the `args.GetProgramDirname`
/// part of `GetDefaultSearchPath`, quoted at `registry_search_path`).
///
/// The search path is built once, by the first use of `CMetaRegistry` (its constructor calls
/// `GetDefaultSearchPath`). Normally that is the usage report's load of `.ncbirc` while the
/// application object is constructed, before `AppMain` gives `CNcbiArguments` the program's
/// name: `GetProgramDirname(eIgnoreLinks)` is then empty (the name `ncbi`) and only the
/// resolved directory of `/proc/<pid>/exe` is searched on Linux (nothing elsewhere). When
/// the usage report does not load `.ncbirc` (`args_known`: `BLAST_USAGE_REPORT` false,
/// `NCBI_DONT_USE_NCBIRC` or `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` set), the first use is
/// `LoadConfig`, after `AppMain` has set the name: the directory of the name as invoked
/// (`FindProgramExecutablePath`, links not followed) comes first, then the resolved one when
/// it differs (re-audit round 2 finding R-1 of session SFd).
///
/// NCBI reference (598d8ae6): c++/include/corelib/metareg.hpp:257-261
/// ```c++
/// inline
/// CMetaRegistry::CMetaRegistry()
/// {
///     GetDefaultSearchPath(x_SetSearchPath());
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbienv.cpp:373-393
/// ```c++
/// const string& CNcbiArguments::GetProgramName(EFollowLinks follow_links) const
/// {
///     if (follow_links) {
///         CFastMutexGuard LOCK(m_ResolvedNameMutex);
///         if ( !m_ResolvedName.size() ) {
/// #ifdef NCBI_OS_LINUX
///             string proc_link = "/proc/" + NStr::IntToString(getpid()) + "/exe";
///             m_ResolvedName = CDirEntry::NormalizePath(proc_link, follow_links);
/// #else
///             m_ResolvedName = CDirEntry::NormalizePath(GetProgramName(eIgnoreLinks), follow_links);
/// #endif
///         }
///         return m_ResolvedName;
///     } else if ( !m_ProgramName.empty() ) {
///         return m_ProgramName;
///     } else if ( m_Args.size() ) {
///         return m_Args[0];
///     } else {
///         static CSafeStatic<string> kDefProgramName;
///         kDefProgramName->assign("ncbi");
///         return kDefProgramName.Get();
///     }
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:881-882,1025-1026
/// ```c++
///     string exepath = FindProgramExecutablePath(argc, argv, &m_RealExePath);
/// ...
///     // Reset command-line args and application name
///     m_Arguments->Reset(argc, argv, exepath, m_RealExePath);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1510-1528
/// ```c++
///     if (real_path) {
///         char buf[PATH_MAX + 1];
///         string procfile = "/proc/" + NStr::IntToString(getpid()) + "/exe";
///         int    ncount   = (int)readlink((procfile).c_str(), buf, PATH_MAX);
///         if (ncount > 0) {
///             real_path->assign(buf, ncount);
///             ...
///             real_path = 0;
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1603-1608
/// ```c++
///     // Save to cache and return
///     *s_Path = ret_val;
///     *s_RealPath = CDirEntry::NormalizePath(ret_val, eFollowLinks);
///     if (real_path) {
///         *real_path = *s_RealPath;
///     }
/// ```
/// On Windows NCBI takes the module's file name (`GetModuleFileName`, `current_exe`) for
/// both names (`NormalizePath` follows no links there). Only Linux was run against NCBI.
fn program_directories(
    context: &ProgramContext,
    path: Option<&OsStr>,
    args_known: bool,
) -> Vec<PathBuf> {
    let exe = context
        .current_exe
        .as_ref()
        .map(|exe| exe.as_os_str().as_encoded_bytes().to_vec());
    let (name, resolved): (Vec<u8>, Vec<u8>) = if !args_known {
        if cfg!(target_os = "linux") {
            (Vec::new(), exe.unwrap_or_default())
        } else {
            (Vec::new(), Vec::new())
        }
    } else if cfg!(windows) {
        let exe = exe.unwrap_or_default();
        (exe.clone(), exe)
    } else {
        let argv0 = context
            .argv0
            .as_ref()
            .map(|argv0| argv0.as_encoded_bytes().to_vec())
            .unwrap_or_default();
        let name = find_program_executable_path(&argv0, context.cwd.as_deref(), path);
        let resolved = match exe.filter(|_| cfg!(target_os = "linux")) {
            Some(exe) => exe,
            None => normalize_path(&name, true),
        };
        (name, resolved)
    };
    let dir = program_dirname(&name);
    let dir2 = program_dirname(&resolved);
    let mut dirs = Vec::new();
    if !dir.is_empty() {
        dirs.push(path_from_bytes(dir));
    }
    if !dir2.is_empty() && dir2 != dir {
        dirs.push(path_from_bytes(dir2));
    }
    dirs
}

/// The directories where NCBI looks for its registry files, in order.
///
/// NCBI reference: ncbi-blast/c++/src/corelib/metareg.cpp:331-396
/// ```c
/// void CMetaRegistry::GetDefaultSearchPath(CMetaRegistry::TSearchPath& path)
/// {
///     path.clear();
///
///     const TXChar* cfg_path = NcbiSys_getenv(_TX("NCBI_CONFIG_PATH"));
///     TSearchPath   path_tail;
///     if (cfg_path) {
///         NStr::Split(_T_STDSTRING(cfg_path), kConfigPathDelim, path);
///         TSearchPath::iterator it = find(path.begin(), path.end(), kEmptyStr);
///         if (it == path.end()) {
///             return;
///         } else {
///             path_tail.assign(it + 1, path.end());
///             path.erase(it, path.end());
///         }
///     }
///
///     if (NcbiSys_getenv(_TX("NCBI_DONT_USE_LOCAL_CONFIG")) == NULL) {
///         path.push_back(".");
///         string home = CDir::GetHome();
///         if ( !home.empty() ) {
///             path.push_back(home);
///         }
///     }
///
///     {{
///         const TXChar* ncbi = NcbiSys_getenv(_TX("NCBI"));
///         if (ncbi  &&  *ncbi) {
///             path.push_back(_T_STDSTRING(ncbi));
///         }
///     }}
///     ...
///     path.push_back("/etc");
///     ...
///             string                dir  = args.GetProgramDirname(eIgnoreLinks);
///             string                dir2 = args.GetProgramDirname(eFollowLinks);
///             if (dir.size()) {
///                 path.push_back(dir);
///             }
///             if (dir2.size() && dir2 != dir) {
///                 path.push_back(dir2);
///             }
///     ...
///     if ( !path_tail.empty() ) {
///         ITERATE (TSearchPath, it, path_tail) {
///             if ( !it->empty() ) {
///                 path.push_back(*it);
///             }
///         }
///     }
/// }
/// ```
/// `kConfigPathDelim` is ":;" outside Windows (metareg.cpp:54-56). The home directory is
/// `home_directory`, the program's directories `program_directories`. `NCBI_CONFIG_PATH` is
/// split as bytes, so a directory whose name is not UTF-8 is kept (finding R-3 of SFd).
///
/// `NStr::Split` adds no token for an empty string, so a set but empty `NCBI_CONFIG_PATH`
/// has no empty element and the search path stays empty: NCBI reads no `<program>.ini` and
/// no `.ncbirc` (audit finding A-3/B-4 of session SFd; Rust's `"".split` yields one empty
/// part, which would splice in the default directories).
///
/// NCBI reference (598d8ae6): c++/include/corelib/ncbistr_util.hpp:268-273
/// ```c++
///         auto target_initial_size = target.size();
///
///         // Special cases
///         if (m_Str.empty()) {
///             return;
///         } else if (m_Delim.empty()) {
/// ```
fn registry_search_path(
    env: &dyn Fn(&str) -> Option<OsString>,
    context: &ProgramContext,
    args_known: bool,
) -> Vec<SearchDir> {
    let delimiters: &[u8] = if cfg!(windows) { b";" } else { b":;" };
    let dir = |bytes: &[u8]| SearchDir::Path(path_from_bytes(bytes));
    let mut path: Vec<SearchDir> = Vec::new();
    let mut tail: Vec<SearchDir> = Vec::new();
    if let Some(config_path) = env("NCBI_CONFIG_PATH") {
        let config_path = config_path.as_encoded_bytes();
        if config_path.is_empty() {
            return path;
        }
        let parts: Vec<&[u8]> = config_path
            .split(|byte| delimiters.contains(byte))
            .collect();
        match parts.iter().position(|part| part.is_empty()) {
            None => return parts.iter().map(|part| dir(part)).collect(),
            Some(empty) => {
                path.extend(parts[..empty].iter().map(|part| dir(part)));
                tail.extend(
                    parts[empty + 1..]
                        .iter()
                        .filter(|part| !part.is_empty())
                        .map(|part| dir(part)),
                );
            }
        }
    }
    if env("NCBI_DONT_USE_LOCAL_CONFIG").is_none() {
        path.push(SearchDir::Path(PathBuf::from(".")));
        match home_directory(env, context) {
            SearchDir::Path(home) if home.as_os_str().is_empty() => {}
            home => path.push(home),
        }
    }
    if let Some(ncbi) = env("NCBI") {
        if !ncbi.is_empty() {
            path.push(SearchDir::Path(PathBuf::from(ncbi)));
        }
    }
    if cfg!(windows) {
        if let Some(root) = env("SYSTEMROOT") {
            if !root.is_empty() {
                path.push(SearchDir::Path(PathBuf::from(root)));
            }
        }
    } else {
        path.push(SearchDir::Path(PathBuf::from("/etc")));
    }
    path.extend(
        program_directories(context, env("PATH").as_deref(), args_known)
            .into_iter()
            .map(SearchDir::Path),
    );
    path.extend(tail);
    path
}

/// One entry of a registry file as NCBI's reader stores it: the section and the name as
/// written (NCBI compares both without case) and the value after its quotes and escapes.
#[derive(Clone, Debug, PartialEq, Eq)]
struct RegistryEntry {
    section: String,
    name: String,
    value: Vec<u8>,
}

/// A registry file that LOSAT does not read as NCBI does.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RegistryError {
    /// NCBI's reader throws `CRegistryException` on this line (1-based, continuation lines
    /// counted): a syntax error.
    Syntax { line: usize, reason: &'static str },
    /// A UTF-16 byte-order mark: NCBI converts the file to UTF-8 first (`ReadIntoUtf8`),
    /// which LOSAT does not port.
    Utf16,
}

/// `isspace` of an `unsigned char` in the C locale (NCBI's programs do not call
/// `setlocale`): the space and 0x09-0x0D, never a byte above 0x7f.
fn is_c_space(byte: u8) -> bool {
    matches!(byte, 0x09..=0x0D | b' ')
}

/// `IRegistry::IsNameSection` (and `IsNameEntry`) without `fInternalSpaces` and
/// `fSectionlessEntries`, the flags of NCBI's reads of `<program>.ini` and `.ncbirc`:
/// a nonempty name of ASCII letters and digits, `_`, `-`, `.` and `/`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:59-81
/// ```c++
/// // Valid symbols for a section/entry name
/// inline bool s_IsNameSectionSymbol(char ch, IRegistry::TFlags flags)
/// {
///     return (isalnum((unsigned char) ch)
///             ||  ch == '_'  ||  ch == '-' ||  ch == '.'  ||  ch == '/'
///             ||  ((flags & IRegistry::fInternalSpaces)  &&  ch == ' '));
/// }
///
///
/// bool IRegistry::IsNameSection(const string& str, TFlags flags)
/// {
///     // Allow empty section name in case of fSectionlessEntries set
///     if (str.empty() && !(flags & IRegistry::fSectionlessEntries) ) {
///         return false;
///     }
///
///     ITERATE (string, it, str) {
///         if (!s_IsNameSectionSymbol(*it, flags)) {
///             return false;
///         }
///     }
///     return true;
/// }
/// ```
fn is_registry_name(name: &[u8]) -> bool {
    !name.is_empty()
        && name
            .iter()
            .all(|&byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.' | b'/'))
}

/// Whether `pos` in `text` follows an odd number of backslashes.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:140-148
/// ```c++
/// // Does pos follow an odd number of backslashes?
/// inline bool s_Backslashed(const string& s, SIZE_TYPE pos)
/// {
///     if (pos == 0) {
///         return false;
///     }
///     SIZE_TYPE last_non_bs = s.find_last_not_of("\\", pos - 1);
///     return (pos - last_non_bs) % 2 == 0;
/// }
/// ```
fn backslashed(text: &[u8], pos: usize) -> bool {
    text[..pos]
        .iter()
        .rev()
        .take_while(|&&byte| byte == b'\\')
        .count()
        % 2
        == 1
}

/// `NStr::ParseEscapes(text)` (`eEscSeqRange_Standard`), or `None` where it throws
/// `CStringException`: a backslash at the end, `\x` without a hexadecimal digit, or a `\x`
/// number above `kMax_UInt` (`NStr::StringToUInt`; a smaller one keeps its low byte).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4796-4884
/// ```c++
///     while (pos < str.size()  ||  !is_error) {
///         SIZE_TYPE pos2 = str.find('\\', pos);
///         if (pos2 == NPOS) {
///             //~ out += str.substr(pos);
///             CTempString sub(str, pos);
///             out += sub;
///             break;
///         }
///         //~ out += str.substr(pos, pos2 - pos);
///         CTempString sub(str, pos, pos2-pos);
///         out += sub;
///         if (++pos2 == str.size()) {
///             NCBI_THROW2(CStringException, eFormat,
///                         "Unterminated escape sequence", pos2);
///         }
///         switch (str[pos2]) {
///         case 'a':  out += '\a';  break;
///         case 'b':  out += '\b';  break;
///         case 'f':  out += '\f';  break;
///         case 'n':  out += '\n';  break;
///         case 'r':  out += '\r';  break;
///         case 't':  out += '\t';  break;
///         case 'v':  out += '\v';  break;
///         case 'x':
///             {{
///                 pos = ++pos2;
///                 while (pos < str.size()
///                        &&  isxdigit((unsigned char) str[pos])) {
///                     pos++;
///                 }
///                 if (pos > pos2) {
///                     SIZE_TYPE len = pos-pos2;
///                     ...
///                     unsigned int value =
///                         StringToUInt(CTempString(str, pos2, len), 0, 16);
///                     ...
///                     out += static_cast<char>(value);
///                 } else {
///                     NCBI_THROW2(CStringException, eFormat,
///                                 "\\x followed by no hexadecimal digits", pos);
///                 }
///             }}
///             continue;
///         case '0':  case '1':  case '2':  case '3':
///         case '4':  case '5':  case '6':  case '7':
///             {{
///                 pos = pos2;
///                 unsigned char c = (unsigned char)(str[pos++] - '0');
///                 while (pos < pos2 + 3  &&  pos < str.size()
///                        &&  str[pos] >= '0'  &&  str[pos] <= '7') {
///                     c = (unsigned char)((c << 3) | (str[pos++] - '0'));
///                 }
///                 out += c;
///             }}
///             continue;
///         case '\n':
///             // quoted EOL means no EOL
///             break;
///         default:
///             out += str[pos2];
///             break;
///         }
///         pos = pos2 + 1;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:647-654
/// ```c++
/// NStr::StringToUInt(const CTempString str, TStringToNumFlags flags, int base)
/// {
///     S2N_CONVERT_GUARD_EX(flags);
///     Uint8 value = StringToUInt8(str, flags, base);
///     if ( value > kMax_UInt ) {
///         S2N_CONVERT_ERROR(unsigned int, "overflow", ERANGE, 0);
///     }
///     return (unsigned int) value;
/// ```
fn parse_escapes(text: &[u8]) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(text.len());
    let mut pos = 0;
    while let Some(found) = text[pos..].iter().position(|&byte| byte == b'\\') {
        out.extend_from_slice(&text[pos..pos + found]);
        let at = pos + found + 1;
        let &escaped = text.get(at)?;
        pos = at + 1;
        match escaped {
            b'a' => out.push(0x07),
            b'b' => out.push(0x08),
            b'f' => out.push(0x0c),
            b'n' => out.push(b'\n'),
            b'r' => out.push(b'\r'),
            b't' => out.push(b'\t'),
            b'v' => out.push(0x0b),
            b'x' => {
                let digits = text[pos..]
                    .iter()
                    .take_while(|byte| byte.is_ascii_hexdigit())
                    .count();
                if digits == 0 {
                    return None;
                }
                let mut value: u64 = 0;
                for &digit in &text[pos..pos + digits] {
                    value = value * 16 + u64::from((digit as char).to_digit(16)?);
                    if value > u64::from(u32::MAX) {
                        return None;
                    }
                }
                out.push(value as u8);
                pos += digits;
            }
            b'0'..=b'7' => {
                let mut value = escaped - b'0';
                let end = (at + 3).min(text.len());
                while pos < end && (b'0'..=b'7').contains(&text[pos]) {
                    value = (value << 3) | (text[pos] - b'0');
                    pos += 1;
                }
                out.push(value);
            }
            b'\n' => {}
            other => out.push(other),
        }
    }
    out.extend_from_slice(&text[pos..]);
    Some(out)
}

/// The lines of a registry file as `NcbiGetlineEOL` reads them: split at `\n` (on Windows
/// without one `\r` before it; on macOS at `\r`, `\n` or `\r\n`); a delimiter at the end of
/// the file starts no line.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:196-208
/// ```c++
/// CNcbiIstream& NcbiGetlineEOL(CNcbiIstream& is, string& str, SIZE_TYPE* count)
/// {
/// #if   defined(NCBI_OS_MSWIN)
///     NcbiGetline(is, str, '\n', count);
///     if (!str.empty()  &&  str[str.length() - 1] == '\r')
///         str.resize(str.length() - 1);
/// #elif defined(NCBI_OS_DARWIN)
///     NcbiGetline(is, str, "\r\n", count);
/// #else /* assume UNIX-like EOLs */
///     NcbiGetline(is, str, '\n', count);
/// #endif //NCBI_OS_...
///     return is;
/// }
/// ```
fn registry_lines(data: &[u8]) -> Vec<&[u8]> {
    fn trim(line: &[u8]) -> &[u8] {
        if cfg!(windows) {
            line.strip_suffix(b"\r").unwrap_or(line)
        } else {
            line
        }
    }
    let mac = cfg!(target_os = "macos");
    let mut lines = Vec::new();
    let mut start = 0;
    let mut pos = 0;
    while pos < data.len() {
        let byte = data[pos];
        if byte == b'\n' || (mac && byte == b'\r') {
            lines.push(trim(&data[start..pos]));
            pos += 1;
            if mac && byte == b'\r' && data.get(pos) == Some(&b'\n') {
                pos += 1;
            }
            start = pos;
        } else {
            pos += 1;
        }
    }
    if start < data.len() {
        lines.push(trim(&data[start..]));
    }
    lines
}

/// What `GetTextEncodingForm(is, eBOM_Discard)` leaves for the reader of a registry file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TextForm<'a> {
    /// The text after a UTF-8 byte-order mark, or all of it.
    Text(&'a [u8]),
    /// A UTF-16 byte-order mark: NCBI converts the rest to UTF-8 (`ReadIntoUtf8`), which
    /// LOSAT does not port.
    Utf16,
    /// A file of one byte 0xEF, 0xFE or 0xFF: the second `get` meets the end of the file and
    /// `unget()` then fails, so the stream is left failed without its end-of-file flag.
    FailedStream,
}

/// `GetTextEncodingForm(is, eBOM_Discard)` on a registry file's `bytes` (re-audit round 2
/// findings R-4 and R-5 of session SFd). A lead byte 0xEF, 0xFE or 0xFF followed by `BB BF`
/// is a UTF-8 mark (the lead byte is not checked again); bytes read without a mark are put
/// back (`Pushback` installs a new buffer through `rdbuf`, which clears the stream's state).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:782-826
/// ```c++
/// EEncodingForm GetTextEncodingForm(CNcbiIstream& input,
///                                   EBOMDiscard   discard_bom)
/// {
///     EEncodingForm ef = eEncodingForm_Unknown;
///     if (input.good()) {
///         const int bom_max = 4;
///         char tmp[bom_max];
///         memset(tmp, 0, bom_max);
///         Uint2* us = reinterpret_cast<Uint2*>(tmp);
///         Uchar* uc = reinterpret_cast<Uchar*>(tmp);
///         input.get(tmp[0]);
///         int n = (int) input.gcount();
///         if (n == 1  &&  (uc[0] == 0xEF  ||  uc[0] == 0xFE  ||  uc[0] == 0xFF)){
///             input.get(tmp[1]);
///             if (input.gcount() == 1) {
///                 ++n;
///                 if (us[0] == 0xFEFF) {
///                     ef = eEncodingForm_Utf16Native;
///                 } else if (us[0] == 0xFFFE) {
///                     ef = eEncodingForm_Utf16Foreign;
///                 } else if (uc[1] == 0xBB) {
///                     input.get(tmp[2]);
///                     if (input.gcount() == 1) {
///                         ++n;
///                         if (uc[2] == 0xBF) {
///                             ef = eEncodingForm_Utf8;
///                         }
///                     }
///                 }
///             }
///         }
///         if (ef == eEncodingForm_Unknown) {
///             if (n > 1) {
///                 CStreamUtils::Pushback(input, tmp, n);
///             } else if (n == 1) {
///                 input.unget();
///             }
///         } else {
///             if (discard_bom == eBOM_Keep) {
///                 CStreamUtils::Pushback(input, tmp, n);
///             }
///         }
///     }
///     return ef;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:143
/// ```c++
///     m_Sb = m_Is.rdbuf(this);
/// ```
/// When the stream already reads from a pushback buffer that holds the bytes, `Pushback`
/// steps back in that buffer instead (no `rdbuf`, the state stays; `read_registry_text`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:411-428
/// ```c++
///         // 2/ [make] equal to the [reasonably-sized, if copy] adjacent part
///         //    of the internal buffer with the data that have just been read?
///         if (how == ePushback_Stepback
///             ||  (how == ePushback_Copy
///                  &&  buf_size <= (del_ptr
///                                   ? CPushback_Streambuf::kMinBufSize
///                                   : CPushback_Streambuf::kMinBufSize >> 4))) {
///             CT_CHAR_TYPE* bp = sb->gptr();
///             size_t avail = bp - sb->m_Buf;
///             size_t take  = avail < buf_size ? avail : buf_size;
///             if (take) {
///                 bp -= take;
///                 buf_size -= take;
///                 if (how != ePushback_Stepback  &&  bp != buf + buf_size) {
///                     memmove(bp, buf + buf_size, take);
///                 }
///                 sb->setg(bp, bp, sb->egptr());
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:323-327
/// ```c++
///         streamsize n = m_Sb->sgetn(bp ? bp : (CT_CHAR_TYPE*) m_DelPtr, r);
///         if (n <= 0) {
///             // NB: For unknown reasons WorkShop6 can return -1 from sgetn :-/
///             delete[] bp;
///             return;
/// ```
fn text_encoding_form(bytes: &[u8]) -> TextForm<'_> {
    match bytes {
        [0xFF, 0xFE, ..] | [0xFE, 0xFF, ..] => TextForm::Utf16,
        [0xEF | 0xFE | 0xFF, 0xBB, 0xBF, rest @ ..] => TextForm::Text(rest),
        [0xEF | 0xFE | 0xFF] => TextForm::FailedStream,
        _ => TextForm::Text(bytes),
    }
}

/// What NCBI writes to stderr when `IRWRegistry::x_Read` ends on a failed stream before its
/// first line (`read_registry_text`): line 1, no path (`SEntry::Reload` and
/// `CNcbiRegistry::x_Read` pass none) and an empty line, at the severity `Error` with the
/// error code of the registry (110) and the subcode 4 (observed with NCBI BLAST+ 2.17.0 for
/// `.ncbirc` and `<program>.ini`, before and after `AppMain` sets up the diagnostics).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:813-816
/// ```c++
///     if ( !is.eof() ) {
///         ERR_POST_X(4, "Error reading the registry after line " << line
///                    << in_path << ": " << str);
///     }
/// ```
/// NCBI reference (598d8ae6): c++/include/corelib/error_codes.hpp:53
/// ```c++
/// NCBI_DEFINE_ERRCODE_X(Corelib_Reg,        110,  8);
/// ```
const REGISTRY_READ_ERROR: &str = "Error: (110.4) Error reading the registry after line 1: \n";

/// A registry file as NCBI's reader leaves it: its entries, and whether `x_Read` reported a
/// failed stream (`Error reading the registry after line 1: `, `REGISTRY_READ_ERROR`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct RegistryText {
    entries: Vec<RegistryEntry>,
    read_error: bool,
}

/// Which reader reads a registry file: the encoding form is checked once for `.ncbirc` and
/// twice for `<program>.ini`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RegistryKind {
    /// `.ncbirc`: `IRWRegistry::Read` on a `CCompoundRWRegistry` (`SEntry::Reload`), whose
    /// `x_Read` is `IRWRegistry::x_Read`.
    Ncbirc,
    /// `<program>.ini`: `IRWRegistry::Read` on the application's `CNcbiRegistry`, whose
    /// `x_Read` passes the stream to `m_FileRegistry->Read`: a second `IRWRegistry::Read`,
    /// which returns at once on a failed stream.
    ProgramIni,
}

/// The registry file `bytes` read by the reader of `kind`: the encoding-form checks, then
/// `registry_entries`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:605-627
/// ```c++
/// IRWRegistry* IRWRegistry::Read(CNcbiIstream& is, TFlags flags,
///                                const string& path)
/// {
///     ...
///     if ( !is ) {
///         return NULL;
///     }
///
///     // Ensure that x_Read gets a stream it can handle.
///     EEncodingForm ef = GetTextEncodingForm(is, eBOM_Discard);
///     if (ef == eEncodingForm_Utf16Native  ||  ef == eEncodingForm_Utf16Foreign) {
///         CStringUTF8 s;
///         ReadIntoUtf8(is, &s, ef);
///         CNcbiIstrstream iss(s);
///         return x_Read(iss, flags, path);
///     } else {
///         return x_Read(is, flags, path);
///     }
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1681-1692
/// ```c++
/// IRWRegistry* CNcbiRegistry::x_Read(CNcbiIstream& is, TFlags flags,
///                                    const string& path)
/// {
///     // Normally, all settings should go to the main portion.  However,
///     // loading an initial configuration file should instead go to the
///     // file portion so that environment settings can take priority.
///     CConstRef<IRegistry> main_reg(FindByName(sm_MainRegName));
///     if (main_reg->Empty()  &&  m_FileRegistry->Empty()) {
///         m_FileRegistry->Read(is, flags & ~fWithNcbirc);
///         LoadBaseRegistries(flags, 0, path);
///         IncludeNcbircIfAllowed(flags);
///         return NULL;
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:2080-2089
/// ```c++
/// IRWRegistry* CCompoundRWRegistry::x_Read(CNcbiIstream& in, TFlags flags,
///                                          const string& path)
/// {
///     ...
///     IRWRegistry::x_Read(in, flags, path);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:813-816
/// ```c++
///     if ( !is.eof() ) {
///         ERR_POST_X(4, "Error reading the registry after line " << line
///                    << in_path << ": " << str);
///     }
/// ```
/// So `<program>.ini` reads after two UTF-8 marks, and a UTF-8 mark followed by a UTF-16 one
/// is converted (LOSAT rejects it as UTF-16); a one-byte 0xEF/0xFE/0xFF `<program>.ini` is
/// empty without a message (the second `Read` sees the failed stream), while a UTF-8 mark and
/// one such byte, or a one-byte `.ncbirc`, is empty with the message (`x_Read` reads no line
/// and the stream has no end-of-file flag). An empty file reads no line and reaches the end
/// of the file: no message.
fn read_registry_text(bytes: &[u8], kind: RegistryKind) -> Result<RegistryText, RegistryError> {
    let failed = |read_error| RegistryText {
        entries: Vec::new(),
        read_error,
    };
    let text = match (text_encoding_form(bytes), kind) {
        (TextForm::Utf16, _) => return Err(RegistryError::Utf16),
        (TextForm::FailedStream, RegistryKind::Ncbirc) => return Ok(failed(true)),
        (TextForm::FailedStream, RegistryKind::ProgramIni) => return Ok(failed(false)),
        (TextForm::Text(text), RegistryKind::Ncbirc) => text,
        // The first check put the two bytes back after meeting the end of the file (a new
        // pushback buffer, state cleared); the second reads them from that buffer, meets the
        // end again and puts them back into the same buffer, which leaves the stream's
        // end-of-file and fail flags set: `x_Read` reads no line and posts nothing.
        (TextForm::Text(_), RegistryKind::ProgramIni)
            if matches!(bytes, [0xEF | 0xFE | 0xFF, 0xBB]) =>
        {
            return Ok(failed(false))
        }
        (TextForm::Text(text), RegistryKind::ProgramIni) => match text_encoding_form(text) {
            TextForm::Utf16 => return Err(RegistryError::Utf16),
            TextForm::FailedStream => return Ok(failed(true)),
            TextForm::Text(text) => text,
        },
    };
    Ok(RegistryText {
        entries: registry_entries(text)?,
        read_error: false,
    })
}

/// The entries of the text of a registry file (after `text_encoding_form`) as NCBI's reader
/// (`IRWRegistry::x_Read`) stores them, in order, or the line of a syntax error (NCBI stops
/// at that line). An entry before the first section is read and checked but not stored
/// (`IRWRegistry::Set` refuses the empty section name without a message).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:652-662
/// ```c++
///     for (line = 1;  NcbiGetlineEOL(is, str);  ++line) {
///         try {
///             SIZE_TYPE len = str.length();
///             SIZE_TYPE beg = 0;
///
///             while (beg < len  &&  isspace((unsigned char) str[beg])) {
///                 ++beg;
///             }
///             // If this line is empty, all comments
///             // that have just been read go to the current section.
///             if (beg == len) {
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:672-733
/// ```c++
///             switch (str[beg]) {
///
///             case '#':  { // file comment
///             ...
///             case ';':  { // section or entry comment
///             ...
///             case '[':  { // section name
///                 ++beg;
///                 SIZE_TYPE end = str.find_first_of(']', beg + 1);
///                 if (end == NPOS) {
///                     NCBI_THROW2(CRegistryException, eSection,
///                                 "Invalid registry section" + in_path
///                                 + " (']' is missing): `" + str + "'", line);
///                 }
///                 section = NStr::TruncateSpaces(str.substr(beg, end - beg));
///                 if (section.empty()) {
///                     NCBI_THROW2(CRegistryException, eSection,
///                                 "Unnamed registry section" + in_path + ": `"
///                                 + str + "'", line);
///                 } else if ( !IsNameSection(section, flags) ) {
///                     NCBI_THROW2(CRegistryException, eSection,
///                                 "Invalid registry section name" + in_path
///                                 + ": `" + str + "'", line);
///                 }
///             ...
///             default:  { // regular entry
///                 string name, value;
///                 if ( !NStr::SplitInTwo(str, "=", name, value) ) {
///                     NCBI_THROW2(CRegistryException, eEntry,
///                                 "Invalid registry entry format" + in_path
///                                 + ": '" + str + "'", line);
///                 }
///                 NStr::TruncateSpacesInPlace(name);
///                 if ( !IsNameEntry(name, flags) ) {
///                     NCBI_THROW2(CRegistryException, eEntry,
///                                 "Invalid registry entry name" + in_path + ": '"
///                                 + str + "'", line);
///                 }
///
///                 NStr::TruncateSpacesInPlace(value);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:743-781
/// ```c++
///                 // read continuation lines, if any
///                 string cont;
///                 while (s_Backslashed(value, value.size())
///                        &&  NcbiGetlineEOL(is, cont)) {
///                     ++line;
///                     value[value.size() - 1] = '\n';
///                     value += NStr::TruncateSpaces(cont);
///                     str   += 'n' + cont; // for presentation in exceptions
///                 }
///
///                 // Historically, " may appear unescaped at the beginning,
///                 // end, both, or neither.
///                 beg = 0;
///                 SIZE_TYPE end = value.size();
///                 for (SIZE_TYPE pos = value.find('\"');
///                      pos < end  &&  pos != NPOS;
///                      pos = value.find('\"', pos + 1)) {
///                     if (s_Backslashed(value, pos)) {
///                         continue;
///                     } else if (pos == beg) {
///                         ++beg;
///                     } else if (pos == end - 1) {
///                         --end;
///                     } else {
///                         NCBI_THROW2(CRegistryException, eValue,
///                                     "Single(unescaped) '\"' in the middle "
///                                     "of registry value" + in_path + ": '"
///                                     + str + "'", line);
///                     }
///                 }
///
///                 try {
///                     value = NStr::ParseEscapes(value.substr(beg, end - beg));
///                 } catch (CStringException&) {
///                     NCBI_THROW2(CRegistryException, eValue,
///                                 "Badly placed '\\' in the registry value"
///                                 + in_path + ": '" + str + "'", line);
///
///                 }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:833-838
/// ```c++
///     string clean_section = NStr::TruncateSpaces(section);
///     if ( !IsNameSection(clean_section, flags) ) {
///         _TRACE("IRWRegistry::Set: bad section name \""
///                << NStr::PrintableString(section) << '\"');
///         return false;
///     }
/// ```
/// The `CRegistryException` leaves `x_Read`, since NCBI's reads of `<program>.ini` and
/// `.ncbirc` do not pass `fIgnoreErrors`.
fn registry_entries(data: &[u8]) -> Result<Vec<RegistryEntry>, RegistryError> {
    let syntax = |line: usize, reason: &'static str| RegistryError::Syntax { line, reason };
    let mut entries = Vec::new();
    let mut section = String::new();
    let mut lines = registry_lines(data).into_iter();
    let mut number = 0;
    while let Some(line) = lines.next() {
        number += 1;
        let Some(beg) = line.iter().position(|&byte| !is_c_space(byte)) else {
            continue;
        };
        match line[beg] {
            b'#' | b';' => {}
            b'[' => {
                let beg = beg + 1;
                let end = line
                    .get(beg + 1..)
                    .and_then(|rest| rest.iter().position(|&byte| byte == b']'))
                    .map(|end| beg + 1 + end)
                    .ok_or_else(|| syntax(number, "a section line without ']'"))?;
                let name = truncate_spaces(&line[beg..end]);
                if name.is_empty() {
                    return Err(syntax(number, "a section without a name"));
                }
                if !is_registry_name(name) {
                    return Err(syntax(number, "an invalid section name"));
                }
                section = String::from_utf8_lossy(name).into_owned();
            }
            _ => {
                let equals = line
                    .iter()
                    .position(|&byte| byte == b'=')
                    .ok_or_else(|| syntax(number, "a line that is not an entry (no '=')"))?;
                let name = truncate_spaces(&line[..equals]);
                if !is_registry_name(name) {
                    return Err(syntax(number, "an invalid entry name"));
                }
                let mut value = truncate_spaces(&line[equals + 1..]).to_vec();
                while backslashed(&value, value.len()) {
                    let Some(next) = lines.next() else { break };
                    number += 1;
                    let last = value.len() - 1;
                    value[last] = b'\n';
                    value.extend_from_slice(truncate_spaces(next));
                }
                let mut beg = 0;
                let mut end = value.len();
                let mut from = 0;
                while let Some(found) = value[from..].iter().position(|&byte| byte == b'"') {
                    let pos = from + found;
                    if pos >= end {
                        break;
                    }
                    if !backslashed(&value, pos) {
                        if pos == beg {
                            beg += 1;
                        } else if pos == end - 1 {
                            end -= 1;
                        } else {
                            return Err(syntax(
                                number,
                                "an unescaped '\"' in the middle of a value",
                            ));
                        }
                    }
                    from = pos + 1;
                }
                let value = parse_escapes(&value[beg..end])
                    .ok_or_else(|| syntax(number, "a badly placed '\\' in a value"))?;
                if !section.is_empty() {
                    entries.push(RegistryEntry {
                        section: section.clone(),
                        name: String::from_utf8_lossy(name).into_owned(),
                        value,
                    });
                }
            }
        }
    }
    Ok(entries)
}

/// The registry file that NCBI loads under `file_name`: the first one on the search path that
/// is a regular file (links followed), or LOSAT's rejection when the search reaches a home
/// directory that LOSAT cannot know (`SearchDir::UnknownHome`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:258-269
/// ```c++
///     if ( dir.empty() ) {
///         ITERATE (TSearchPath, it, m_SearchPath) {
///             const string& result
///                 = x_FindRegistry(CDirEntry::MakePath(*it, name), style);
///             if ( !result.empty() ) {
///                 return result;
///             }
///         }
///     } else {
///         switch (style) {
///         case eName_AsIs:
///             if (CFile(name).Exists()) {
/// ```
fn find_registry(search_path: &[SearchDir], file_name: &str) -> Result<Option<PathBuf>, String> {
    for dir in search_path {
        match dir {
            SearchDir::Path(dir) => {
                let path = dir.join(file_name);
                if path.is_file() {
                    return Ok(Some(path));
                }
            }
            SearchDir::UnknownHome => {
                return Err(format!(
                    "HOME is not set and the user has no passwd entry, so NCBI BLAST+ would look for its registry file {file_name} in the home directory of the login name (USER, LOGNAME or getlogin), which LOSAT does not look up; this is not supported by LOSAT"
                ))
            }
        }
    }
    Ok(None)
}

/// The registry file at `path` as the reader of `kind` reads it, `None` when it cannot be
/// opened (NCBI then loads no registry from it and does not look further: `SEntry::Reload`
/// fails after `x_FindRegistry` has chosen the file; re-audit round 2 finding R-6 of session
/// SFd), or LOSAT's rejection of a file that NCBI reports as a syntax error, that is UTF-16,
/// or that LOSAT could open but not read. On a syntax error NCBI writes a message with the
/// text of the exception (for `.ncbirc` `Critical: ... Syntax error in system-wide
/// configuration file: NCBI C++ Exception:` with the path and line of NCBI's own source file,
/// once for each reader of the file, and goes on with the entries before the line; for
/// `<program>.ini` `Error: (CRegistryException::...) ...` and exit code 2), which LOSAT does
/// not reproduce (`AUTHORITY.md` §J of `docs/evidence/losat_web_e2h/`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:62-84
/// ```c++
/// bool CMetaRegistry::SEntry::Reload(CMetaRegistry::TFlags reload_flags)
/// {
///     CFile file(actual_name);
///     if ( !file.Exists() ) {
///         _TRACE("No such registry file " << actual_name);
///         return false;
///     }
///     ...
///     CNcbiIfstream ifs(actual_name.c_str(), IOS_BASE::in | IOS_BASE::binary);
///     if ( !ifs.good() ) {
///         _TRACE("Unable to (re)open registry file " << actual_name);
///         return false;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:228-231
/// ```c++
///     if (scratch_entry.actual_name.empty()
///         ||  !scratch_entry.Reload(flags | fAlwaysReload | fKeepContents) ) {
///         scratch_entry.registry.Reset();
///         return scratch_entry;
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1659-1664
/// ```c++
///     } catch (CRegistryException& e) {
///         ERR_POST_X(6, Critical << "CNcbiRegistry: "
///                       "Syntax error in system-wide configuration file: "
///                       << e.what());
///         return false;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:669-673
/// ```c++
///             catch (const CException& e) {
///                 NCBI_REPORT_EXCEPTION_X(15,
///                                         "Application's initialization failed", e);
///                 *got_exception = true;
///                 *exit_code = 2;
/// ```
fn read_registry_file(path: &Path, kind: RegistryKind) -> Result<Option<RegistryText>, String> {
    use std::io::Read;
    let Ok(mut file) = std::fs::File::open(path) else {
        return Ok(None);
    };
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes).map_err(|error| {
        format!(
            "NCBI BLAST+ would read the registry file {}, which LOSAT could open but not read ({error}); this is not supported by LOSAT",
            path.display()
        )
    })?;
    read_registry_text(&bytes, kind)
        .map(Some)
        .map_err(|error| match error {
            RegistryError::Syntax { line, reason } => format!(
                "the registry file {} has {reason} on line {line}, which NCBI BLAST+ reports as a syntax error; this is not supported by LOSAT",
                path.display()
            ),
            RegistryError::Utf16 => format!(
                "the registry file {} is UTF-16, which NCBI BLAST+ converts to UTF-8 before reading it; this is not supported by LOSAT",
                path.display()
            ),
        })
}

/// Rejects an entry of a registry file that LOSAT does not know to change no output.
fn check_registry_entries(path: &Path, entries: &[RegistryEntry]) -> Result<(), String> {
    for entry in entries {
        if !is_harmless_entry(
            &entry.section,
            &entry.name,
            &String::from_utf8_lossy(&entry.value),
        ) {
            return Err(format!(
                "the registry file {} sets [{}] {}, which LOSAT does not know to change no output of NCBI BLAST+ (LOSAT accepts only [BLAST] BLASTDB, BLASTMAT, DATA_LOADERS, BLASTDB_NUCL_DATA_LOADER, BLASTDB_PROT_DATA_LOADER, IGDATA, MAX_SEQID_LENGTH, BLAST_USAGE_REPORT with a Boolean, [NCBI] DATA and DONT_USE_NCBIRC); this is not supported by LOSAT",
                path.display(),
                entry.section,
                entry.name
            ));
        }
    }
    Ok(())
}

/// The NCBI application settings that LOSAT reproduces (`check_ncbi_application_settings`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ApplicationSettings {
    /// `SDataLoaderConfig::UseDataLoaders()` of the query and subject readers: whether
    /// NCBI's reader tries a record's first line as a sequence identifier
    /// (`fasta_reader::ReaderConfig::data_loaders`).
    pub data_loaders: bool,
}

impl Default for ApplicationSettings {
    /// No registry entry and no variable: both data loaders on.
    ///
    /// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_scope_src.hpp:73,81
    /// ```c++
    ///         eDefault = (eUseBlastDbDataLoader | eUseGenbankDataLoader)
    /// ...
    ///     SDataLoaderConfig(bool load_proteins, EConfigOpts options = eDefault)
    /// ```
    fn default() -> Self {
        Self { data_loaders: true }
    }
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
/// The first layer that has the entry answers, also with an empty value:
/// `CCompoundRegistry::FindByContents` asks each layer with `fCountCleared`, with which a
/// variable set to an empty value (`CEnvironmentRegistry::x_HasEntry`), an empty value of the
/// program's `.ini` (`CMemoryRegistry::x_HasEntry`) and an empty value that `.ncbirc` wrote
/// into the application's registry (recorded as cleared, `CCompoundRWRegistry::x_HasEntry`)
/// all count. An empty value contains neither `blastdb` nor `genbank`, so it turns both
/// loaders off (audit findings B-1 of session SFc, A-1/B-3 of SFd; `AUTHORITY.md` §G1).
/// Whether an empty value of `.ncbirc` reaches the application's registry at all is the
/// caller's question (`application_settings`, finding A-2 of SFd).
///
/// NCBI reference (598d8ae6): c++/src/corelib/env_reg.cpp:157-167
/// ```c++
/// bool CEnvironmentRegistry::x_HasEntry(const string& section,
///                                       const string& name,
///                                       TFlags flags) const
/// {
///     if (name.empty()) {
///         return x_HasSection(section, flags);
///     }
///     bool found = false;
///     x_Get(section, name, flags, found);
///     return found;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbienv.cpp:97-104,117-124
/// ```c++
///     for ( ;  *envp;  envp++) {
///         const char* s = *envp;
///         const char* eq = strchr(s, '=');
///         ...
///         m_Cache[string(s, (size_t)(eq - s))] = SEnvValue(eq + 1, kEmptyXCStr);
/// ...
///     if ( i != m_Cache.end() ) {
///         if (i->second.ptr == NULL  &&  i->second.value.empty()) {
///             *found = false;
///             return kEmptyStr;
///         } else {
///             *found = true;
///             return i->second.value;
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:984-991
/// ```c++
///     TEntries::const_iterator eit = entries.find(name);
///     if (eit == entries.end()) {
///         return false;
///     } else if ((flags & fCountCleared) != 0) {
///         return true;
///     } else {
///         return !eit->second.value.empty();
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1939-1947
/// ```c++
/// bool CCompoundRWRegistry::x_HasEntry(const string& section, const string& name,
///                                      TFlags flags) const
/// {
///     TClearedEntries::const_iterator it
///         = m_ClearedEntries.find(s_FlatKey(section, name));
///     if (it != m_ClearedEntries.end()) {
///         if ((flags & fCountCleared)  &&  (flags & it->second)) {
///             return true;
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:2041-2045
/// ```c++
///     if (value.empty()) {
///         bool was_empty = Get(section, name, flags).empty();
///         m_MainRegistry->Set(section, name, value, flags, comment);
///         m_ClearedEntries[s_FlatKey(section, name)] |= flags2;
///         return !was_empty;
/// ```
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:75-92
/// ```c++
/// void
/// SDataLoaderConfig::x_LoadDataLoadersConfig(const CNcbiRegistry& registry)
/// {
///     static const string kDataLoadersConfig("DATA_LOADERS");
///
///     if (registry.HasEntry("BLAST", kDataLoadersConfig)) {
///         const string& kLoaders = registry.Get("BLAST", kDataLoadersConfig);
///         if (NStr::FindNoCase(kLoaders, "blastdb") == NPOS) {
///             m_UseBlastDbs = false;
///         }
///         if (NStr::FindNoCase(kLoaders, "genbank") == NPOS) {
///             m_UseGenbank = false;
///         }
///         if (NStr::FindNoCase(kLoaders, "none") != NPOS) {
///             m_UseBlastDbs = false;
///             m_UseGenbank = false;
///         }
///     }
/// ```
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_scope_src.hpp:111
/// ```c++
///     bool UseDataLoaders() const { return m_UseBlastDbs || m_UseGenbank; }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1577-1583
/// ```c++
///     x_Add(*m_EnvRegistry, ePriority_Environment, sm_EnvRegName);
///
///     m_FileRegistry.Reset(new CTwoLayerRegistry(NULL, cf));
///     x_Add(*m_FileRegistry, ePriority_File, sm_FileRegName);
///
///     m_SysRegistry.Reset(new CCompoundRWRegistry(cf));
///     x_Add(*m_SysRegistry, ePriority_Default - 1, sm_SysRegName);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1235-1246
/// ```c++
/// CConstRef<IRegistry> CCompoundRegistry::FindByContents(const string& section,
///                                                        const string& entry,
///                                                        TFlags flags) const
/// {
///     TFlags has_entry_flags = (flags | fCountCleared) & ~fJustCore;
///     REVERSE_ITERATE(TPriorityMap, it, m_PriorityMap) {
///         if (it->second->HasEntry(section, entry, has_entry_flags)) {
///             return it->second;
///         }
///     }
///     return null;
/// }
/// ```
fn data_loaders_of(env: Option<&[u8]>, files: [Option<&[u8]>; 2]) -> bool {
    let Some(value) = env.or_else(|| files.into_iter().flatten().next()) else {
        return ApplicationSettings::default().data_loaders;
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

/// The value of `[BLAST] DATA_LOADERS` that a registry file stores: that of its last
/// entry (a later entry replaces an earlier one, `Set` without `fNoOverride`; an empty
/// value is stored too).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:794-801
/// ```c++
///                 } else if (was_empty  &&  HasEntry(section, name, flags)) {
///                     ERR_POST_X(8, Warning
///                                << "Found multiple [" << section << "] "
///                                << name << " settings" << in_path
///                                << "; using the one from line " << line);
///                 }
///                 Set(section, name, value, set_flags, comment);
///                 comment.erase();
/// ```
fn file_data_loaders(entries: &[RegistryEntry]) -> Option<&[u8]> {
    entries
        .iter()
        .rev()
        .find(|entry| {
            entry.section.eq_ignore_ascii_case("BLAST")
                && entry.name.eq_ignore_ascii_case("DATA_LOADERS")
        })
        .map(|entry| entry.value.as_slice())
}

/// Rejects the NCBI application settings that change `program`'s output: the environment
/// variables of the NCBI C++ Toolkit and the entries of the registry files that NCBI loads
/// (`<program>.ini`, then `.ncbirc` unless it is turned off). Returns the settings that
/// LOSAT reproduces: whether the data loaders are used (`[BLAST] DATA_LOADERS`, which
/// decides whether a record's first line can be a sequence identifier).
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbiapp.cpp:1254-1259
/// ```c
///     } else if (conf->empty()) {
///         entry = CMetaRegistry::Load(basename, CMetaRegistry::eName_Ini, 0,
///                                     reg_flags, &reg);
/// ```
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbireg.cpp:1636-1648
/// ```c
///     if (flags & fWithNcbirc) {
///         flags &= ~fWithNcbirc;
///     } else {
///         return false;
///     }
///
///     if (getenv("NCBI_DONT_USE_NCBIRC")) {
///         return false;
///     }
///
///     if (HasEntry("NCBI", "DONT_USE_NCBIRC")) {
///         return false;
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/corelib/metareg.cpp:299-303
/// ```c
///         case eName_DotRc: {
///             string base, ext;
///             CDirEntry::SplitPath(name, 0, &base, &ext);
///             return x_FindRegistry(CDirEntry::MakePath(dir, '.' + base, ext)
///                                   + "rc", eName_AsIs);
/// ```
///
/// NCBI's readers of the registry files write `REGISTRY_READ_ERROR` to stderr once for each
/// read of a file that leaves the stream failed (`read_registry_text`); LOSAT writes them
/// before the program's own output, as NCBI does while it starts.
pub fn check_ncbi_application_settings(program: &str) -> Result<ApplicationSettings, String> {
    let (settings, read_errors) = application_settings(
        program,
        std::env::vars_os(),
        &|name: &str| std::env::var_os(name),
        &ProgramContext::of_this_process(),
    )?;
    if read_errors > 0 {
        use std::io::Write;
        let mut stderr = std::io::stderr().lock();
        for _ in 0..read_errors {
            let _ = stderr.write_all(REGISTRY_READ_ERROR.as_bytes());
        }
        let _ = stderr.flush();
    }
    Ok(settings)
}

/// `check_ncbi_application_settings` over the environment `variables` (all of them), `env`
/// (one of them by name) and the program's `context`; returns the settings and the number of
/// `REGISTRY_READ_ERROR` lines that NCBI's readers of the registry files write.
///
/// NCBI reference (598d8ae6): c++/src/corelib/env_reg.cpp:366-375
/// ```c++
/// string CNcbiEnvRegMapper::RegToEnv(const string& section, const string& name)
///     const
/// {
///     string result;
///     result.assign(kPrefix, kPrefixLen);
///     if (NStr::StartsWith(name, ".")) {
///         result += name.substr(1) + "__" + section;
///     } else {
///         result += "_" + section + "__" + name;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/env_reg.cpp:142-153
/// ```c++
///     REVERSE_ITERATE (TPriorityMap, it, m_PriorityMap) {
///         string        var_name = it->second->RegToEnv(section, name);
///         const string* resultp  = &m_Env->Get(var_name, &found);
///         if ((m_Flags & fCaseFlags) == 0  &&  !found) {
///             // try capitalizing the name
///             resultp = &m_Env->Get(NStr::ToUpper(var_name), &found);
///         }
///         if (found) {
///             return *resultp;
///         }
///     }
/// ```
/// The environment layer of `HasEntry("NCBI", "DONT_USE_NCBIRC")` (`IncludeNcbircIfAllowed`,
/// quoted above) is `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`, set with any value.
///
/// Who reads `.ncbirc`, and whether its empty values reach the application's registry
/// (audit finding A-2 of session SFd). Every BLAST program has a `CBlastUsageReport` member,
/// built before `AppMain` runs. Unless the variable `BLAST_USAGE_REPORT` (the environment
/// only) is a false Boolean, it builds a registry with `fWithNcbirc`, which loads `.ncbirc`
/// through `CMetaRegistry::Load` with no flags and keeps it in `CMetaRegistry`'s cache. The
/// application then loads `<program>.ini` (`LoadConfig`). Without one it loads `.ncbirc`
/// under the same key, finds the cached registry and copies it into its own through `Write`
/// and `Read`; `Write` enumerates the entries without `fCountCleared`, so an empty value
/// does not reach the application. With a `<program>.ini`, `.ncbirc` is loaded from
/// `CNcbiRegistry::x_Read` with `fJustCore` (which `SEntry::Reload` adds), a key the cache
/// does not have, so the file is read into the application's registry itself and an empty
/// value counts (`data_loaders_of`). A syntax error makes the usage report's load fail
/// and leaves nothing in the cache; LOSAT rejects such a file anyway. The usage report reads
/// `.ncbirc` also when `<program>.ini` turns it off for the application (`[NCBI]
/// DONT_USE_NCBIRC`): its syntax errors are still reported and its `[BLAST]
/// BLAST_USAGE_REPORT` is still converted to a Boolean (another value aborts NCBI). The copy
/// changes no other value: `Printable` writes octal escapes that `ParseEscapes` reads back.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:81
/// ```c++
///     CBlastUsageReport m_UsageReport;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:195-211
/// ```c++
/// void CBlastUsageReport::x_CheckBlastUsageEnv()
/// {
///     char * blast_usage_env = getenv("BLAST_USAGE_REPORT");
///     if(blast_usage_env != NULL){
///         bool enable = NStr::StringToBool(blast_usage_env);
///         if (!enable) {
///             SetEnabled(false);
///             CUsageReportAPI::SetEnabled(false);
///             LOG_POST(Info <<"Phone home disabled");
///             return ;
///         }
///     }
///
///     CNcbiIstrstream empty_stream(kEmptyStr);
///     CRef<CNcbiRegistry> registry(new CNcbiRegistry(empty_stream, IRegistry::fWithNcbirc));
///     if (registry->HasEntry("BLAST", "BLAST_USAGE_REPORT")) {
///         bool enable = NStr::StringToBool(registry->Get("BLAST", "BLAST_USAGE_REPORT"));
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1622-1625
/// ```c++
///     x_Init();
///     m_FileRegistry->Read(is, flags & ~(fWithNcbirc | fCaseFlags));
///     LoadBaseRegistries(flags, 0, path);
///     IncludeNcbircIfAllowed(flags & ~fCaseFlags);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1650-1653
/// ```c++
///     try {
///         CMetaRegistry::SEntry entry
///             = CMetaRegistry::Load("ncbi", CMetaRegistry::eName_RcOrIni,
///                                   0, flags, m_SysRegistry.GetPointer());
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1254-1256,1278-1279
/// ```c++
///     } else if (conf->empty()) {
///         entry = CMetaRegistry::Load(basename, CMetaRegistry::eName_Ini, 0,
///                                     reg_flags, &reg);
///         // still consider pulling in defaults from .ncbirc
///         if (reg.IncludeNcbircIfAllowed(reg_flags)) {
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:89-90
/// ```c++
///         if ((reload_flags & fKeepContents)  ||  registry->Empty(rflags)) {
///             dest = registry->Read(ifs, reg_flags | IRegistry::fJustCore);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1688-1691
/// ```c++
///     if (main_reg->Empty()  &&  m_FileRegistry->Empty()) {
///         m_FileRegistry->Read(is, flags & ~fWithNcbirc);
///         LoadBaseRegistries(flags, 0, path);
///         IncludeNcbircIfAllowed(flags);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:200-208
/// ```c++
///     else { // see if we already have it
///         TIndex::const_iterator iit
///             = m_Index.find(SKey(name, style, flags, reg_flags));
///         if (iit != m_Index.end()) {
///             _TRACE("found in cache");
///             _ASSERT(iit->second < m_Contents.size());
///             SEntry& result = m_Contents[iit->second];
///             result.Reload(flags);
///             return result;
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:152-171
/// ```c++
///     if (reg  &&  entry.registry  &&  reg != entry.registry.GetPointer()) {
///         ...
///         TStrStream str;
///         entry.registry->Write(str, rflags);
///         str.seekg(0);
///         ...
///         reg->Read(str, reg_flags | IRegistry::fJustCore);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:226-234
/// ```c++
///         list<string> entries;
///         EnumerateEntries(*section, &entries, flags);
///         ITERATE (list<string>, entry, entries) {
///             ...
///             os << *entry << " = \""
///                << Printable(Get(*section, *entry, flags)) << "\""
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1997-2004
/// ```c++
///         ITERATE (list<string>, it2, tmp) {
///             // avoid reporting cleared entries
///             TClearedEntries::const_iterator ceci
///                 = (flags & fCountCleared) ? m_ClearedEntries.end()
///                 : m_ClearedEntries.find(s_FlatKey(section, *it2));
///             if (ceci == m_ClearedEntries.end()
///                 ||  (flags & ~fJustCore & ~ceci->second)) {
///                 accum.insert(*it2);
/// ```
fn application_settings(
    program: &str,
    variables: impl Iterator<Item = (OsString, OsString)>,
    env: &dyn Fn(&str) -> Option<OsString>,
    context: &ProgramContext,
) -> Result<(ApplicationSettings, usize), String> {
    for (name, value) in variables {
        if let Some(reason) = rejected_variable(&name, &value) {
            return Err(reason);
        }
    }
    let ncbirc_allowed = env("NCBI_DONT_USE_NCBIRC").is_none()
        && env("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC").is_none();
    // `rejected_variable` has rejected a BLAST_USAGE_REPORT that is not a Boolean.
    let usage_report_reads_ncbirc = ncbirc_allowed
        && env("BLAST_USAGE_REPORT")
            .is_none_or(|value| ncbi_string_to_bool(&value.to_string_lossy()) != Some(false));
    // The search path is built by the first use of CMetaRegistry: the usage report's load of
    // .ncbirc before the program's name is known, or else LoadConfig (`program_directories`).
    let search_path = registry_search_path(env, context, !usage_report_reads_ncbirc);
    let mut read_errors = 0;
    let mut application_reads_ncbirc = ncbirc_allowed;
    // A <program>.ini that cannot be opened is no registry: LoadConfig goes on as without one.
    let ini = match find_registry(&search_path, &format!("{program}.ini"))? {
        Some(path) => read_registry_file(&path, RegistryKind::ProgramIni)?.map(|text| (path, text)),
        None => None,
    };
    let ini_entries = match &ini {
        Some((path, text)) => {
            check_registry_entries(path, &text.entries)?;
            read_errors += usize::from(text.read_error);
            if text.entries.iter().any(|entry| {
                entry.section.eq_ignore_ascii_case("NCBI")
                    && entry.name.eq_ignore_ascii_case("DONT_USE_NCBIRC")
            }) {
                application_reads_ncbirc = false;
            }
            text.entries.clone()
        }
        None => Vec::new(),
    };
    // The application reads .ncbirc itself unless it copies the usage report's cached
    // registry (no <program>.ini and the usage report loaded the file).
    let ncbirc_read_directly = ini.is_some() || !usage_report_reads_ncbirc;
    let mut ncbirc_entries = Vec::new();
    if usage_report_reads_ncbirc || application_reads_ncbirc {
        let ncbirc = match find_registry(&search_path, ".ncbirc")? {
            Some(path) => read_registry_file(&path, RegistryKind::Ncbirc)?.map(|text| (path, text)),
            None => None,
        };
        if let Some((path, text)) = ncbirc {
            if text.read_error {
                read_errors += usize::from(usage_report_reads_ncbirc)
                    + usize::from(application_reads_ncbirc && ncbirc_read_directly);
            }
            if application_reads_ncbirc {
                check_registry_entries(&path, &text.entries)?;
                ncbirc_entries = text.entries;
            } else if env("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT").is_none() {
                // Only the usage report reads the file: its last [BLAST] BLAST_USAGE_REPORT,
                // empty included (read into the usage report's own registry).
                let usage = text.entries.iter().rev().find(|entry| {
                    entry.section.eq_ignore_ascii_case("BLAST")
                        && entry.name.eq_ignore_ascii_case("BLAST_USAGE_REPORT")
                });
                if let Some(entry) = usage {
                    check_registry_entries(&path, std::slice::from_ref(entry))?;
                }
            }
        }
    }
    let env_data_loaders = env("NCBI_CONFIG__BLAST__DATA_LOADERS");
    Ok((
        ApplicationSettings {
            data_loaders: data_loaders_of(
                env_data_loaders.as_deref().map(OsStr::as_encoded_bytes),
                [
                    file_data_loaders(&ini_entries),
                    file_data_loaders(&ncbirc_entries)
                        .filter(|value| ncbirc_read_directly || !value.is_empty()),
                ],
            ),
        },
        read_errors,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toolkit_variables_are_rejected_and_harmless_entries_accepted() {
        let reject =
            |name: &str, value: &str| rejected_variable(OsStr::new(name), OsStr::new(value));
        assert!(reject("DIAG_POST_LEVEL", "Error").is_some());
        assert!(reject("DIAG_OLD_POST_FORMAT", "false").is_some());
        assert!(reject("NCBI_CONFIG__DIAG__POST_FILTER", "!BLAST").is_some());
        assert!(reject("NCBI_CONFIG__BLAST__LONG_SEQID", "1").is_some());
        assert!(reject("NCBI_CONFIG__NCBI__MemorySizeLimit", "1").is_some());
        assert!(reject("NCBI_CONFIG_OVERRIDES", "/tmp/x").is_some());
        assert!(reject("NCBI_CONFIG_LONG_SEQID__BLAST", "1").is_some());
        assert!(reject("NCBI_CONFIG_PATH", "/etc").is_none());
        assert!(reject("ABORT_ON_THROW", "1").is_some());
        assert!(reject("BLAST_USAGE_REPORT", "bogus").is_some());
        assert!(reject("BLAST_USAGE_REPORT", "false").is_none());
        assert!(reject("NCBI_CONFIG__BLAST__BLASTDB", "/db").is_none());
        assert!(reject("NCBI_CONFIG__blast__blastdb", "/db").is_none());
        assert!(reject("BLASTDB", "/db").is_none());
        assert!(reject("NCBI", "/opt/ncbi").is_none());
        assert!(reject("HOME", "/home/x").is_none());
    }

    fn entry(section: &str, name: &str, value: &[u8]) -> RegistryEntry {
        RegistryEntry {
            section: section.into(),
            name: name.into(),
            value: value.to_vec(),
        }
    }

    fn syntax_line(text: &[u8]) -> Option<usize> {
        match registry_entries(text) {
            Err(RegistryError::Syntax { line, .. }) => Some(line),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:652-781 (IRWRegistry::x_Read:
    // C-locale white space, comments, sections, entries, continuation lines, quotes and
    // escapes) and :833-838 (IRWRegistry::Set ignores an entry without a section).
    #[test]
    fn registry_files_are_read_like_ncbis() {
        let entries = registry_entries(
            b"# file comment\n; entry comment\n[BLAST]\nBLASTDB = /db \\\n  /more\n  \n[ncbi]\ndata=/data\n",
        )
        .unwrap();
        assert_eq!(
            entries,
            vec![
                entry("BLAST", "BLASTDB", b"/db \n/more"),
                entry("ncbi", "data", b"/data")
            ]
        );
        assert!(entries
            .iter()
            .all(|e| is_harmless_entry(&e.section, &e.name, "")));
        assert!(!is_harmless_entry("BLAST", "LONG_SEQID", "1"));
        assert!(!is_harmless_entry("DEBUG", "DIAG_POST_LEVEL", "Error"));
        assert!(is_harmless_entry("blast", "blast_usage_report", "Off"));
        assert!(!is_harmless_entry("BLAST", "BLAST_USAGE_REPORT", "bogus"));
        // Syntax errors (CRegistryException), with NCBI's line numbers.
        for (text, line) in [
            (&b"garbage line without equals\n"[..], 1),
            (b"[BLAST\n", 1),
            (b"[]\n", 1),
            (b"[ ]\n", 1),
            (b"[]]\n", 1),
            (b"[B@D]\n", 1),
            (b"[BL AST]\n", 1),
            (b"[\xc2\xa0BLAST]\n", 1),
            (b"\xc2\xa0[BLAST]\n", 1),
            (b"[BLAST]\n\xc2\xa0DATA_LOADERS = none\n", 2),
            (b"[BLAST]\n\xe3\x80\x80DATA_LOADERS = none\n", 2),
            (b"[BLAST]\n\xc2\x85DATA_LOADERS = none\n", 2),
            (b"[BLAST]\nDATA LOADERS=none\n", 2),
            (b"[BLAST]\nDATA_LOADERS\x00=none\n", 2),
            (b"[BLAST]\n=none\n", 2),
            (b"\x00[BLAST]\n", 1),
            (b"[BLAST]\nBLASTDB = a\"b\n", 2),
            (b"[BLAST]\nBLASTDB = \"a\"b\"\n", 2),
            (b"[BLAST]\nBLASTDB = a\"\"\n", 2),
            (b"[BLAST]\nBLASTDB = a\\x\n", 2),
            (b"[BLAST]\nBLASTDB = \\x100000000\n", 2),
            (b"[BLAST]\nBLASTDB = a\\", 2),
            (b"[BLAST]\nX = a\\\nb\\\nc\\", 4),
            (b"DATA LOADERS=none\n", 1),
        ] {
            assert_eq!(syntax_line(text), Some(line), "{text:?}");
        }
        // Accepted: C-locale white space only (VT and FF too), the text after ']', sections
        // of the name symbols, a UTF-8 byte-order mark, an entry before the first section
        // (read, not stored), a line end CR (white space), a value of the other bytes as is.
        let value = |text: &[u8]| {
            let entries = registry_entries(text).unwrap();
            entries
                .iter()
                .rev()
                .find(|e| e.name.eq_ignore_ascii_case("DATA_LOADERS"))
                .map(|e| e.value.clone())
        };
        assert_eq!(
            value(b"[BLAST]\n\x0bDATA_LOADERS=\x0cnone\x0c\n"),
            Some(b"none".to_vec())
        );
        assert_eq!(
            value(b"[BLAST]] junk\nDATA_LOADERS=none\n"),
            Some(b"none".to_vec())
        );
        assert_eq!(
            value(b"[x.y-z_1/w]\n[BLAST]\nDATA_LOADERS=none\n"),
            Some(b"none".to_vec())
        );
        assert_eq!(value(b"DATA_LOADERS=none\n"), None);
        assert_eq!(registry_entries(b"LONG_SEQID=1\n").unwrap(), vec![]);
        assert_eq!(
            value(b"[BLAST]\r\nDATA_LOADERS = none\r\n"),
            Some(b"none".to_vec())
        );
        assert_eq!(
            value(b"[BLAST]\nDATA_LOADERS = \xc2\xa0\n"),
            Some(b"\xc2\xa0".to_vec())
        );
        assert_eq!(
            value(b"[BLAST]\nDATA_LOADERS==genbank=\n"),
            Some(b"=genbank=".to_vec())
        );
        assert_eq!(
            value(b"[BLAST]\nDATA_LOADERS=genbank # none\n"),
            Some(b"genbank # none".to_vec())
        );
        assert_eq!(value(b"[BLAST]\nDATA_LOADERS=\n"), Some(Vec::new()));
        // A lone CR is no line end outside macOS.
        if !cfg!(target_os = "macos") {
            assert_eq!(value(b"[BLAST]\rDATA_LOADERS=none\r"), None);
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:782-826 (GetTextEncodingForm),
    // c++/src/corelib/ncbireg.cpp:605-627, 1681-1692, 2080-2089 (one check for .ncbirc, two
    // for <program>.ini) and :813-816 (the message after a failed stream); re-audit round 2
    // findings R-4 and R-5 of session SFd.
    #[test]
    fn byte_order_marks_are_checked_as_ncbis_readers_do() {
        use RegistryKind::{Ncbirc, ProgramIni};
        let loaders = |text: &[u8], kind| {
            read_registry_text(text, kind).map(|read| {
                (
                    read.entries
                        .iter()
                        .rev()
                        .find(|e| e.name.eq_ignore_ascii_case("DATA_LOADERS"))
                        .map(|e| e.value.clone()),
                    read.read_error,
                )
            })
        };
        let none = Some(b"none".to_vec());
        let body = b"[BLAST]\nDATA_LOADERS=none\n";
        let with = |prefix: &[u8]| [prefix, &body[..]].concat();
        for kind in [Ncbirc, ProgramIni] {
            assert_eq!(loaders(body, kind), Ok((none.clone(), false)));
            // A UTF-8 mark is skipped, also after the lead bytes 0xFE and 0xFF.
            for lead in [0xEF, 0xFE, 0xFF] {
                assert_eq!(
                    loaders(&with(&[lead, 0xBB, 0xBF]), kind),
                    Ok((none.clone(), false))
                );
            }
            // An incomplete mark is put back and read as text.
            assert!(matches!(
                read_registry_text(&with(b"\xef\xbb"), kind),
                Err(RegistryError::Syntax { line: 1, .. })
            ));
            assert!(matches!(
                read_registry_text(b"\xfe", kind),
                Ok(RegistryText { ref entries, .. }) if entries.is_empty()
            ));
            assert_eq!(
                read_registry_text(&with(b"\xff\xfe"), kind),
                Err(RegistryError::Utf16)
            );
            assert_eq!(
                read_registry_text(&with(b"\xfe\xff"), kind),
                Err(RegistryError::Utf16)
            );
            assert_eq!(loaders(b"", kind), Ok((None, false)));
            assert_eq!(loaders(b"\xef\xbb\xbf", kind), Ok((None, false)));
        }
        // One byte 0xEF/0xFE/0xFF: empty; the message for .ncbirc only.
        for byte in [0xEF, 0xFE, 0xFF] {
            assert_eq!(loaders(&[byte], Ncbirc), Ok((None, true)));
            assert_eq!(loaders(&[byte], ProgramIni), Ok((None, false)));
            // After a UTF-8 mark: .ncbirc reads the byte as a line, <program>.ini checks again.
            assert!(matches!(
                read_registry_text(&[0xEF, 0xBB, 0xBF, byte], Ncbirc),
                Err(RegistryError::Syntax { line: 1, .. })
            ));
            assert_eq!(
                loaders(&[0xEF, 0xBB, 0xBF, byte], ProgramIni),
                Ok((None, true))
            );
        }
        // [lead, BB] alone: .ncbirc reads it as a line; <program>.ini's second check leaves
        // the stream failed at its end (empty, no message); after a UTF-8 mark it is a line.
        for lead in [0xEF, 0xFE, 0xFF] {
            assert!(matches!(
                read_registry_text(&[lead, 0xBB], Ncbirc),
                Err(RegistryError::Syntax { line: 1, .. })
            ));
            assert_eq!(loaders(&[lead, 0xBB], ProgramIni), Ok((None, false)));
            assert!(matches!(
                read_registry_text(&[0xEF, 0xBB, 0xBF, lead, 0xBB], ProgramIni),
                Err(RegistryError::Syntax { line: 1, .. })
            ));
        }
        // Two UTF-8 marks, or a UTF-8 mark and a UTF-16 one: <program>.ini checks twice.
        let twice = with(b"\xef\xbb\xbf\xef\xbb\xbf");
        assert_eq!(loaders(&twice, ProgramIni), Ok((none.clone(), false)));
        assert!(matches!(
            read_registry_text(&twice, Ncbirc),
            Err(RegistryError::Syntax { line: 1, .. })
        ));
        let utf16: Vec<u8> = [
            &b"\xef\xbb\xbf\xff\xfe"[..],
            &"[BLAST]\n"
                .encode_utf16()
                .flat_map(u16::to_le_bytes)
                .collect::<Vec<u8>>(),
        ]
        .concat();
        assert_eq!(
            read_registry_text(&utf16, ProgramIni),
            Err(RegistryError::Utf16)
        );
        assert!(matches!(
            read_registry_text(&utf16, Ncbirc),
            Err(RegistryError::Syntax { line: 1, .. })
        ));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:743-781 (continuation lines
    // after an odd number of backslashes, the quotes at the ends) and
    // c++/src/corelib/ncbistr.cpp:4796-4884 (NStr::ParseEscapes).
    #[test]
    fn registry_values_lose_their_quotes_and_escapes() {
        let value = |raw: &[u8]| {
            let mut text = b"[BLAST]\nDATA_LOADERS = ".to_vec();
            text.extend_from_slice(raw);
            text.push(b'\n');
            registry_entries(&text).map(|entries| entries[0].value.clone())
        };
        for (raw, stored) in [
            (&b"none"[..], &b"none"[..]),
            (b"\"none\"", b"none"),
            (b"\"none", b"none"),
            (b"none\"", b"none"),
            (b"\"\"", b""),
            (b"\"", b""),
            (b"\"\"\"", b""),
            (b"", b""),
            (b"\" none \"", b" none "),
            (b"\\\"genbank\\\"", b"\"genbank\""),
            (b"\\x6eone", b"none"),
            (b"gen\\x62ank", b"gen*nk"),
            (b"gen\\x62\\x61nk", b"genbank"),
            (b"\\x16eone", b"none"),
            (b"\\x000000006eone", b"none"),
            (b"\\xFFFFFFFFg", b"\xffg"),
            (b"\\156one", b"none"),
            (b"\\556one", b"none"),
            (b"\\0156", b"\x0d6"),
            (b"gen\\bank", b"gen\x08ank"),
            (b"gen\\\\bank", b"gen\\bank"),
            (b"\\genbank", b"genbank"),
            (b"\\a\\f\\n\\r\\t\\v", b"\x07\x0c\n\r\t\x0b"),
            (b"genbank\\\\", b"genbank\\"),
        ] {
            assert_eq!(value(raw), Ok(stored.to_vec()), "{raw:?}");
        }
        for raw in [
            &b"no\"ne"[..],
            b"\"a\"b\"",
            b"\\x",
            b"a\\xg",
            b"\\x100000000",
        ] {
            assert!(value(raw).is_err(), "{raw:?}");
        }
        // Continuation lines: an odd number of backslashes joins the next line (trimmed)
        // after a line feed; an even number does not.
        let entries =
            registry_entries(b"[BLAST]\nA=x\\\\\\\n  genbank  \nB=none\\\\\n[x]\nC=gen\\\n# c\n")
                .unwrap();
        assert_eq!(
            entries,
            vec![
                entry("BLAST", "A", b"x\\\ngenbank"),
                entry("BLAST", "B", b"none\\"),
                entry("x", "C", b"gen\n# c"),
            ]
        );
        // A backslash at the end of the file is a badly placed '\'.
        assert!(registry_entries(b"[BLAST]\nDATA_LOADERS=genbank\\").is_err());
    }

    fn context(
        argv0: Option<&str>,
        current_exe: Option<&str>,
        passwd_home: Option<&str>,
    ) -> ProgramContext {
        ProgramContext {
            argv0: argv0.map(OsString::from),
            current_exe: current_exe.map(PathBuf::from),
            cwd: std::env::current_dir().ok(),
            passwd_home: passwd_home.map(OsString::from),
        }
    }

    fn dirs(path: &[SearchDir]) -> Vec<PathBuf> {
        path.iter()
            .map(|dir| match dir {
                SearchDir::Path(path) => path.clone(),
                SearchDir::UnknownHome => PathBuf::from("<unknown home>"),
            })
            .collect()
    }

    #[test]
    fn the_search_path_follows_ncbi_config_path() {
        let env = |pairs: &'static [(&'static str, &'static str)]| {
            move |name: &str| {
                pairs
                    .iter()
                    .find(|(key, _)| *key == name)
                    .map(|(_, value)| OsString::from(value))
            }
        };
        let none = context(None, None, None);
        let search = |pairs: &'static [(&'static str, &'static str)]| {
            dirs(&registry_search_path(&env(pairs), &none, false))
        };
        let only = search(&[("NCBI_CONFIG_PATH", "/a:/b")]);
        assert_eq!(only, vec![PathBuf::from("/a"), PathBuf::from("/b")]);
        let spliced = search(&[("NCBI_CONFIG_PATH", "/a::/z"), ("HOME", "/h")]);
        assert_eq!(
            &spliced[..3],
            &[PathBuf::from("/a"), PathBuf::from("."), PathBuf::from("/h")]
        );
        assert_eq!(spliced.last(), Some(&PathBuf::from("/z")));
        let no_local = search(&[("NCBI_DONT_USE_LOCAL_CONFIG", ""), ("NCBI", "/n")]);
        assert_eq!(&no_local[..1], &[PathBuf::from("/n")]);
        // A set but empty NCBI_CONFIG_PATH has no token, so no directory is searched
        // (ncbistr_util.hpp:268-273; audit finding A-3/B-4 of session SFd); ":" has two empty
        // tokens and splices in the default directories.
        assert!(search(&[("NCBI_CONFIG_PATH", ""), ("HOME", "/h")]).is_empty());
        let colon = search(&[("NCBI_CONFIG_PATH", ":"), ("HOME", "/h")]);
        assert_eq!(&colon[..2], &[PathBuf::from("."), PathBuf::from("/h")]);
        // NCBI_CONFIG_PATH is split as bytes (finding R-3 of SFd).
        #[cfg(unix)]
        {
            use std::os::unix::ffi::OsStrExt;
            let bytes = |name: &str| {
                (name == "NCBI_CONFIG_PATH")
                    .then(|| OsStr::from_bytes(b"/x\xffy:/z").to_os_string())
            };
            assert_eq!(
                dirs(&registry_search_path(&bytes, &none, false)),
                vec![
                    PathBuf::from(OsStr::from_bytes(b"/x\xffy")),
                    PathBuf::from("/z")
                ]
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3586-3657 (CDir::GetHome: HOME
    // when set, also empty; otherwise the passwd entry's directory, then the login name);
    // re-audit round 2 finding R-2 of session SFd.
    #[test]
    fn the_home_directory_falls_back_to_the_passwd_entry() {
        if cfg!(windows) {
            return;
        }
        let env = |home: Option<&'static str>| {
            move |name: &str| {
                (name == "HOME")
                    .then_some(home)
                    .flatten()
                    .map(OsString::from)
            }
        };
        let home = |home: Option<&'static str>, passwd: Option<&str>| {
            let path = registry_search_path(&env(home), &context(None, None, passwd), false);
            path.get(1).cloned()
        };
        let at = |dir: &str| Some(SearchDir::Path(PathBuf::from(dir)));
        assert_eq!(home(Some("/h"), Some("/p")), at("/h"));
        assert_eq!(home(None, Some("/p")), at("/p"));
        // An empty HOME, or an empty passwd directory, gives no home directory.
        assert_eq!(home(Some(""), Some("/p")), at("/etc"));
        assert_eq!(home(None, Some("")), at("/etc"));
        // No passwd entry: NCBI looks the login name up, which LOSAT does not.
        assert_eq!(home(None, None), Some(SearchDir::UnknownHome));
        let path = registry_search_path(&env(None), &context(None, None, None), false);
        assert!(find_registry(&path[..1], ".ncbirc").unwrap().is_none());
        assert!(find_registry(&path, ".ncbirc")
            .unwrap_err()
            .contains("passwd"));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbienv.cpp:408-415 (GetProgramDirname) and
    // c++/src/corelib/ncbifile.cpp:820-971 (NormalizePath).
    #[test]
    fn program_names_are_cut_and_normalized_like_ncbis() {
        assert_eq!(program_dirname(b"/a/b/LOSAT"), b"/a/b/");
        assert_eq!(program_dirname(b"LOSAT"), b"");
        assert_eq!(program_dirname(b"/a/my:dir/LOSAT"), b"/a/my:dir/");
        assert_eq!(program_dirname(b"/a/my:LOSAT"), b"/a/my:");
        assert_eq!(program_dirname(b"/a/b\\LOSAT"), b"/a/b\\");
        for (path, normal) in [
            (&b""[..], &b""[..]),
            (b"/", b"/"),
            (b"///", b"/"),
            (b"/a/b/", b"/a/b"),
            (b"/a//b", b"/a/b"),
            (b"/a/./b", b"/a/b"),
            (b"/a/../b", b"/b"),
            (b"/../a", b"/a"),
            (b"/a/b/..", b"/a"),
            (b"./a", b"a"),
            (b"a/./b/", b"a/b"),
            (b".", b"."),
            (b"./", b"."),
            (b"a/..", b"."),
            (b"../a", b"../a"),
            (b"../../a", b"../../a"),
            (b"a/../../b", b"../b"),
        ] {
            assert_eq!(normalize_path(path, false), normal, "{path:?}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1426-1609
    // (FindProgramExecutablePath), metareg.cpp:375-388 and ncbienv.cpp:373-393 (the program's
    // directories before and after AppMain sets the name); re-audit round 2 finding R-1.
    #[cfg(unix)]
    #[test]
    fn the_program_directories_depend_on_when_the_search_path_is_built() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-program-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let dir = std::fs::canonicalize(&dir).unwrap();
        let real = dir.join("real");
        let sym = dir.join("sym");
        std::fs::create_dir_all(&real).unwrap();
        std::fs::create_dir_all(&sym).unwrap();
        let exe = real.join("LOSAT");
        std::fs::write(&exe, b"").unwrap();
        let link = sym.join("LOSAT");
        let _ = std::fs::remove_file(&link);
        std::os::unix::fs::symlink(&exe, &link).unwrap();
        let text = |path: &Path| path.as_os_str().as_encoded_bytes().to_vec();
        let argv0 = text(&link);
        let path_var = OsString::from(format!("/nonexistent::{}", sym.display()));
        // An absolute argv[0] is kept, links not followed.
        assert_eq!(find_program_executable_path(&argv0, None, None), argv0);
        // A bare name is searched in PATH (empty elements skipped).
        assert_eq!(
            find_program_executable_path(b"LOSAT", None, Some(&path_var)),
            argv0
        );
        // A name found nowhere stays as it is.
        assert_eq!(
            find_program_executable_path(b"nowhere/x", None, None),
            b"nowhere/x"
        );
        let cwd = std::env::current_dir().unwrap();
        let ctx = ProgramContext {
            argv0: Some(link.clone().into_os_string()),
            current_exe: Some(exe.clone()),
            cwd: Some(cwd),
            passwd_home: None,
        };
        let with_slash = |path: &Path| {
            let mut bytes = text(path);
            bytes.push(b'/');
            path_from_bytes(&bytes)
        };
        if cfg!(target_os = "linux") {
            assert_eq!(
                program_directories(&ctx, None, false),
                vec![with_slash(&real)]
            );
        } else {
            assert!(program_directories(&ctx, None, false).is_empty());
        }
        assert_eq!(
            program_directories(&ctx, None, true),
            vec![with_slash(&sym), with_slash(&real)]
        );
        // A direct start: one directory.
        let direct = ProgramContext {
            argv0: Some(exe.clone().into_os_string()),
            ..ctx
        };
        assert_eq!(
            program_directories(&direct, None, true),
            vec![with_slash(&real)]
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:75-92
    // (case-insensitive substrings `blastdb`, `genbank` and `none`) and
    // c++/src/corelib/ncbireg.cpp:1235-1246, 1577-1583 (the environment, then the program's
    // `.ini`, then `.ncbirc`); the first layer with the entry answers, an empty value included
    // (fCountCleared: env_reg.cpp:157-167, ncbireg.cpp:984-991, 1939-1947; audit findings B-1
    // of session SFc and A-1/B-3 of SFd).
    #[test]
    fn data_loaders_follow_the_registry_layers_and_the_substring_rules() {
        let only = |value: &str| data_loaders_of(None, [None, Some(value.as_bytes())]);
        for (value, used) in [
            ("blastdb", true),
            ("genbank", true),
            ("BlastDB", true),
            ("GENBANK", true),
            ("blastdb,genbank", true),
            ("genbank blastdb", true),
            ("myblastdbs", true),
            ("none", false),
            ("NONE", false),
            ("blastdb,none", false),
            ("genbank none blastdb", false),
            ("nonexistent", false),
            ("x", false),
            ("blast db", false),
            ("gen bank", false),
            ("", false),
        ] {
            assert_eq!(only(value), used, "{value:?}");
        }
        let layers = |env: Option<&str>, ini: Option<&str>, ncbirc: Option<&str>| {
            data_loaders_of(
                env.map(str::as_bytes),
                [ini.map(str::as_bytes), ncbirc.map(str::as_bytes)],
            )
        };
        assert!(layers(None, None, None));
        assert!(!layers(Some("none"), Some("blastdb"), Some("genbank")));
        assert!(layers(Some("genbank"), Some("none"), Some("none")));
        assert!(!layers(None, Some("none"), Some("blastdb")));
        assert!(layers(None, Some("blastdb"), Some("none")));
        assert!(!layers(None, None, Some("none")));
        // An empty value is an entry in every layer: it answers first and turns both
        // loaders off.
        assert!(!layers(None, Some(""), Some("genbank")));
        assert!(!layers(None, None, Some("")));
        assert!(!layers(Some(""), None, None));
        assert!(!layers(Some(""), Some("blastdb"), Some("genbank")));
        assert!(!layers(Some(""), None, Some("genbank")));
        assert!(layers(Some("genbank"), Some(""), None));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1636-1651 (`.ncbirc` is read
    // unless NCBI_DONT_USE_NCBIRC is set or the registry has [NCBI] DONT_USE_NCBIRC, also
    // from the environment), c++/src/corelib/env_reg.cpp:366-375 (the variable of
    // [BLAST] DATA_LOADERS is NCBI_CONFIG__BLAST__DATA_LOADERS) and
    // c++/src/algo/blast/api/blast_usage_report.cpp:195-211, metareg.cpp:152-171,200-208
    // (an empty value of `.ncbirc` counts unless the application copied the usage report's
    // cached registry: no `<program>.ini` and BLAST_USAGE_REPORT not false).
    #[test]
    fn application_settings_read_the_environment_and_the_registry_files() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-data-loaders-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let ini = dir.join("blastn.ini");
        let ncbirc = dir.join(".ncbirc");
        let settings = |pairs: &[(&str, &str)]| {
            let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
            all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_settings(
                "blastn",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &context(None, None, None),
            )
            .map(|(settings, _)| settings.data_loaders)
        };
        let _ = std::fs::remove_file(&ini);
        let _ = std::fs::remove_file(&ncbirc);
        assert_eq!(settings(&[]), Ok(true));
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        // The environment comes first; an empty variable is an entry that turns both
        // loaders off (env_reg.cpp:157-167).
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "blastdb")]),
            Ok(true)
        );
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
            Ok(false)
        );
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = genbank\n").unwrap();
        assert_eq!(settings(&[]), Ok(true));
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
            Ok(false)
        );
        // An empty value of .ncbirc: no entry when the application copies the usage
        // report's cached registry, an entry when it reads the file itself
        // (BLAST_USAGE_REPORT false, or a <program>.ini).
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS =\n").unwrap();
        assert_eq!(settings(&[]), Ok(true));
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "true")]), Ok(true));
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "Off")]), Ok(false));
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "0")]), Ok(false));
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT", "false")]),
            Ok(true)
        );
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
            Ok(false)
        );
        std::fs::write(&ini, "[BLAST]\nBLASTDB = /nonexistent\n").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        std::fs::write(&ini, "").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        std::fs::remove_file(&ini).unwrap();
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
        // Only the exact variable name is the entry.
        assert_eq!(
            settings(&[("NCBI_CONFIG__blast__data_loaders", "blastdb")]),
            Ok(false)
        );
        // The program's .ini comes before .ncbirc; its last entry counts, without quotes,
        // and an empty value is an entry.
        std::fs::write(
            &ini,
            "[blast]\ndata_loaders = none\n[BLAST]\nDATA_LOADERS = \"GenBank\"\n",
        )
        .unwrap();
        assert_eq!(settings(&[]), Ok(true));
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = genbank\n").unwrap();
        std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = \"\"\n").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = genbank\nDATA_LOADERS =\n").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = no\"ne\n").unwrap();
        assert!(settings(&[]).unwrap_err().contains("syntax error"));
        // [NCBI] DONT_USE_NCBIRC in the .ini or the environment turns .ncbirc off for the
        // application.
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
        std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
        assert_eq!(settings(&[]), Ok(true));
        // The usage report still reads it: a syntax error or a BLAST_USAGE_REPORT that is
        // not a Boolean is rejected unless BLAST_USAGE_REPORT is false; other entries are not
        // read.
        std::fs::write(&ncbirc, "[B@D]\n").unwrap();
        assert!(settings(&[]).unwrap_err().contains("syntax error"));
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "no")]), Ok(true));
        std::fs::write(&ncbirc, "[BLAST]\nBLAST_USAGE_REPORT =\n").unwrap();
        assert!(settings(&[]).is_err());
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "no")]), Ok(true));
        std::fs::write(&ncbirc, "[BLAST]\nLONG_SEQID = 1\n").unwrap();
        assert_eq!(settings(&[]), Ok(true));
        std::fs::remove_file(&ini).unwrap();
        assert!(settings(&[]).is_err());
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
        assert_eq!(
            settings(&[("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC", "")]),
            Ok(true)
        );
        assert_eq!(settings(&[("NCBI_DONT_USE_NCBIRC", "1")]), Ok(true));
        std::fs::remove_file(&ncbirc).unwrap();
        std::fs::remove_dir(&dir).unwrap();
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:813-816 (the message of a failed
    // stream, once for each read of the file: the usage report's, and the application's
    // unless it copies the usage report's cached registry) and metareg.cpp:62-84 (a file that
    // cannot be opened is no registry, and the search does not go on); re-audit round 2
    // findings R-5, R-6 and R-7 of session SFd.
    #[test]
    fn registry_read_errors_and_unreadable_files_follow_ncbi() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-read-errors-{}",
            std::process::id()
        ));
        let home = dir.join("home");
        std::fs::create_dir_all(&home).unwrap();
        let ini = dir.join("blastn.ini");
        let ncbirc = dir.join(".ncbirc");
        let settings = |pairs: &[(&str, &str)]| {
            let mut all = vec![(
                "NCBI_CONFIG_PATH".to_string(),
                format!("{}:{}", dir.display(), home.display()),
            )];
            all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_settings(
                "blastn",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &context(None, None, None),
            )
            .map(|(settings, read_errors)| (settings.data_loaders, read_errors))
        };
        let _ = std::fs::remove_file(&ini);
        std::fs::write(&ncbirc, [0xFF]).unwrap();
        assert_eq!(settings(&[]), Ok((true, 1)));
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "false")]), Ok((true, 1)));
        assert_eq!(settings(&[("NCBI_DONT_USE_NCBIRC", "")]), Ok((true, 0)));
        std::fs::write(&ini, b"[BLAST]\nBLASTDB=/x\n").unwrap();
        assert_eq!(settings(&[]), Ok((true, 2)));
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "no")]), Ok((true, 1)));
        std::fs::write(&ini, b"[NCBI]\nDONT_USE_NCBIRC=1\n").unwrap();
        assert_eq!(settings(&[]), Ok((true, 1)));
        assert_eq!(settings(&[("BLAST_USAGE_REPORT", "no")]), Ok((true, 0)));
        // A one-byte <program>.ini is empty and silent; after a UTF-8 mark it has the message.
        std::fs::write(&ini, [0xFE]).unwrap();
        assert_eq!(settings(&[]), Ok((true, 2)));
        std::fs::write(&ini, [0xEF, 0xBB, 0xBF, 0xFE]).unwrap();
        assert_eq!(settings(&[]), Ok((true, 3)));
        // Still a <program>.ini: .ncbirc is read into the application's registry itself.
        std::fs::write(&ini, [0xFE]).unwrap();
        std::fs::write(&ncbirc, b"[BLAST]\nDATA_LOADERS=\n").unwrap();
        assert_eq!(settings(&[]), Ok((false, 0)));
        std::fs::remove_file(&ini).unwrap();
        assert_eq!(settings(&[]), Ok((true, 0)));
        // A file that cannot be opened is no registry and hides the next one.
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::write(home.join(".ncbirc"), b"[BLAST]\nDATA_LOADERS=none\n").unwrap();
            std::fs::write(&ncbirc, b"[BLAST]\nLONG_SEQID=1\n").unwrap();
            std::fs::set_permissions(&ncbirc, std::fs::Permissions::from_mode(0o000)).unwrap();
            if std::fs::File::open(&ncbirc).is_err() {
                assert_eq!(settings(&[]), Ok((true, 0)));
                std::fs::write(&ini, b"[BLAST]\nLONG_SEQID=1\n").ok();
                std::fs::set_permissions(&ini, std::fs::Permissions::from_mode(0o000)).unwrap();
                assert_eq!(settings(&[]), Ok((true, 0)));
                std::fs::set_permissions(&ini, std::fs::Permissions::from_mode(0o644)).unwrap();
                std::fs::remove_file(&ini).unwrap();
            }
            std::fs::set_permissions(&ncbirc, std::fs::Permissions::from_mode(0o644)).unwrap();
            std::fs::remove_file(home.join(".ncbirc")).unwrap();
        }
        // An empty NCBI_CONFIG_OVERRIDES is not set.
        std::fs::remove_file(&ncbirc).unwrap();
        assert_eq!(settings(&[("NCBI_CONFIG_OVERRIDES", "")]), Ok((true, 0)));
        assert!(settings(&[("NCBI_CONFIG_OVERRIDES", "/x")]).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn ncbis_booleans() {
        for (text, value) in [
            ("1", Some(true)),
            ("TRUE", Some(true)),
            ("t", Some(true)),
            ("On", Some(true)),
            ("0", Some(false)),
            ("No", Some(false)),
            ("off", Some(false)),
            ("", None),
            ("2", None),
        ] {
            assert_eq!(ncbi_string_to_bool(text), value, "{text:?}");
        }
    }
}
