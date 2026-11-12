//! The configuration that NCBI's application layer (`CNcbiApplication`) reads before a BLAST
//! program runs: environment variables of the NCBI C++ Toolkit and the registry files
//! (`.ncbirc` and the program's `.ini`). LOSAT does not reproduce the settings that change
//! a program's output or exit status, and rejects them explicitly (plan DW-13); a registry
//! file with only entries that change no output (such as `[BLAST] BLASTDB`) is accepted.

use std::ffi::{OsStr, OsString};
use std::io::Write;
use std::path::{Path, PathBuf};

use super::ncbi_registry::{
    error_reading_post, last_value, read_registry, repeated_entry_post, running_error_reading_post,
    running_syntax_error_notice, syntax_error_notice, Layer, RegistryEntry, RegistryRead,
};
use crate::report::diag_log::{Posts, CRITICAL, ERROR, WARNING};

/// Registry entries (section, name; case-insensitive) that change no output of a search with
/// `-subject`: the database and data-loader settings, the data directory, the directory of
/// the WindowMasker files of blastn's `-window_masker_taxid` (`window_masker_path_of`), and
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
    ("WINDOW_MASKER", "WINDOW_MASKER_PATH"),
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

fn is_harmless_entry(section: &str, name: &str) -> bool {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_usage_report.cpp:209-211
    // ```c
    // 	CRef<CNcbiRegistry> registry(new CNcbiRegistry(empty_stream, IRegistry::fWithNcbirc));
    // 	if (registry->HasEntry("BLAST", "BLAST_USAGE_REPORT")) {
    // 		bool enable = NStr::StringToBool(registry->Get("BLAST", "BLAST_USAGE_REPORT"));
    // ```
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
/// ```c++
/// CBlastUsageReport::CBlastUsageReport()
/// {
/// 	x_CheckBlastUsageEnv();
/// ...
/// void CBlastUsageReport::x_CheckBlastUsageEnv()
/// {
/// 	char * blast_usage_env = getenv("BLAST_USAGE_REPORT");
/// 	if(blast_usage_env != NULL){
/// 		bool enable = NStr::StringToBool(blast_usage_env);
/// 		if (!enable) {
/// ...
/// 			return ;
/// 		}
/// 	}
///
/// 	CNcbiIstrstream empty_stream(kEmptyStr);
/// 	CRef<CNcbiRegistry> registry(new CNcbiRegistry(empty_stream, IRegistry::fWithNcbirc));
/// ```
/// NCBI BLAST+ 2.17.0 (`blastn` run by a link on `PATH`, a malformed `.ncbirc` or
/// `blastn.ini`): one in the link's directory is read only with `BLAST_USAGE_REPORT=false` or
/// `=0`, or (`blastn.ini`) with `NCBI_DONT_USE_NCBIRC` or `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`;
/// one in the directory of the executable in every case.
fn registry_search_path(
    program: &str,
    env: &dyn Fn(&str) -> Option<OsString>,
    process: &ProcessFacts,
) -> Result<Vec<PathBuf>, String> {
    let delimiters: &[char] = if cfg!(windows) { &[';'] } else { &[':', ';'] };
    let mut path: Vec<PathBuf> = Vec::new();
    let mut tail: Vec<PathBuf> = Vec::new();
    if let Some(config_path) = env("NCBI_CONFIG_PATH") {
        let config_path = config_path.to_string_lossy().into_owned();
        if config_path.is_empty() {
            return Ok(Vec::new());
        }
        let parts: Vec<&str> = config_path.split(delimiters).collect();
        match parts.iter().position(|part| part.is_empty()) {
            None => return Ok(parts.iter().map(PathBuf::from).collect()),
            Some(empty) => {
                path.extend(parts[..empty].iter().map(PathBuf::from));
                tail.extend(
                    parts[empty + 1..]
                        .iter()
                        .filter(|part| !part.is_empty())
                        .map(PathBuf::from),
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
}

impl ProcessFacts<'static> {
    /// The facts of this process.
    fn of_this_process() -> Self {
        Self {
            // Rust's `home_dir` on Unix reads `HOME` and, when it is not set, the `pw_dir` of
            // `getpwuid_r(getuid())`; it is asked only when `HOME` is not set (`ncbi_home`).
            home_by_uid: &|| std::env::home_dir().map(PathBuf::into_os_string),
            program_dirs: &program_directories,
        }
    }
}

/// `CDir::GetHome()`, the directory that `registry_search_path` adds after ".".
///
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
/// NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3586-3618
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
///
/// static bool s_GetHomeByLOGIN(string& home)
/// {
///     const TXChar* ptr = 0;
///     // Get user name
///     if ( !(ptr = NcbiSys_getenv(_TX("USER"))) ) {
///         if ( !(ptr = NcbiSys_getenv(_TX("LOGNAME"))) ) {
///             if ( !(ptr = getlogin()) ) {
/// ...
///     struct passwd* pwd = getpwnam(ptr);
/// ```
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
/// ```c++
///             string                dir  = args.GetProgramDirname(eIgnoreLinks);
///             string                dir2 = args.GetProgramDirname(eFollowLinks);
///             if (dir.size()) {
///                 path.push_back(dir);
///             }
///             if (dir2.size() && dir2 != dir) {
///                 path.push_back(dir2);
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbienv.cpp:373-415
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
/// ...
/// string CNcbiArguments::GetProgramDirname(EFollowLinks follow_links) const
/// {
///     const string& name = GetProgramName(follow_links);
///     SIZE_TYPE base_pos = name.find_last_of("/\\:");
///     if (base_pos == NPOS)
///         return NcbiEmptyString;
///     return name.substr(0, base_pos + 1);
/// }
/// ```
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
/// ```c++
///     if (argc > 0  &&  argv[0] != NULL  &&  argv[0][0] != '\0') {
///         ret_val = argv[0];
///     } else if (instance) {
///         ret_val = instance->GetArguments().GetProgramName();
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1549-1595
/// ```c++
///     string app_path = ret_val;
///
///     if ( !CDirEntry::IsAbsolutePath(app_path) ) {
/// ...
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
/// ...
///             NStr::Split(env_path, ":", split_path,
///                 NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
/// #  endif
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
/// ```c++
///     if ( path.empty() ) {
///         return path;
///     }
/// ...
///     current = DeleteTrailingPathSeparator(path);
/// #  endif
///
///     if ( current.empty() ) {
///         // root dir
///         return string(1, DIR_SEPARATOR);
///     }
/// ...
///     while ( !current.empty()  ||  !tail.empty() ) {
///         std::list<string> pretail;
///         if ( !current.empty() ) {
///             NStr::Split(current, kSep, pretail);
///             current.erase();
///             if (pretail.front().empty()
/// ...
///                 ) {
///                 // Absolute path
///                 head.clear();
/// ...
///             tail.splice(tail.begin(), pretail);
///         }
///
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
/// ...
///             } else if (next == DIR_CURRENT) {
///                 // Leave out, since we already have content
///                 continue;
///             } else if (next.empty()) {
///                 continue; // leave out empty components in most cases
///             } else if (next == DIR_PARENT) {
/// ...
///                 // Back up if possible, assuming existing path to be "physical"
///                 if (last.empty()) {
///                     // Already at the root; .. is a no-op
///                     continue;
/// ...
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
/// ...
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
/// ...
///     // Compose path
///     return NStr::Join(head, string(1, DIR_SEPARATOR));
/// ```
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
/// ```c++
///     string dir;
///     CDirEntry::SplitPath(name, &dir, 0, 0);
///     if ( dir.empty() ) {
///         ITERATE (TSearchPath, it, m_SearchPath) {
///             const string& result
///                 = x_FindRegistry(CDirEntry::MakePath(*it, name), style);
///             if ( !result.empty() ) {
///                 return result;
///             }
///         }
/// ```
fn find_registry(search_path: &[PathBuf], file_name: &str) -> Option<PathBuf> {
    search_path
        .iter()
        .map(|dir| dir.join(file_name))
        .find(|path| path.is_file())
}

// NCBI reference (598d8ae6): c++/include/corelib/metareg.hpp:64-69
// ```c++
//         /// C Toolkit style; mostly useful with name = "ncbi"
// #ifdef NCBI_OS_MSWIN
//         eName_RcOrIni = eName_Ini
// #else
//         eName_RcOrIni = eName_DotRc
// #endif
// ```
// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:284-303 (`eName_Ini` adds `.ini`,
// `eName_DotRc` makes `.<name>rc`)
// ```c++
//         case eName_Ini:
//             for (string name2(name); ; ) {
//                 string result = x_FindRegistry(name2 + ".ini", eName_AsIs);
//         ...
//         case eName_DotRc: {
//             string base, ext;
//             CDirEntry::SplitPath(name, 0, &base, &ext);
//             return x_FindRegistry(CDirEntry::MakePath(dir, '.' + base, ext)
//                                   + "rc", eName_AsIs);
// ```
/// The name of NCBI's system-wide registry file, which `CMetaRegistry::Load("ncbi",
/// eName_RcOrIni)` looks for: `.ncbirc`, and `ncbi.ini` on Windows.
const SYSTEM_WIDE_REGISTRY: &str = if cfg!(windows) { "ncbi.ini" } else { ".ncbirc" };

/// One read of the registry file at `path` into `layer` (`ncbi_registry::read_registry`):
/// `Ok(None)` when NCBI cannot open it, which NCBI takes for no file (and looks no further).
/// `Err` holds LOSAT's refusal.
///
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:60-90 (`SEntry::Reload`, which
/// `x_Load` calls for the file found)
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
///     ...
///         if ((reload_flags & fKeepContents)  ||  registry->Empty(rflags)) {
///             dest = registry->Read(ifs, reg_flags | IRegistry::fJustCore);
/// ```
/// NCBI BLAST+ 2.17.0 (g_reg probes `unr`, `unr2`, 2026-11-11): a malformed `.ncbirc` that it
/// cannot open (mode 000) is not reported, and a `.ncbirc` in the home directory is then not
/// read either; a `blastn.ini` that it cannot open is no file. The file is opened in binary mode:
/// the line ends are `NcbiGetlineEOL`'s alone.
fn read_registry_file(path: &Path, layer: Layer) -> Result<Option<RegistryRead>, String> {
    let Ok(mut file) = std::fs::File::open(path) else {
        return Ok(None);
    };
    let mut bytes = Vec::new();
    std::io::Read::read_to_end(&mut file, &mut bytes).map_err(|error| {
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
/// ```c++
///     try {
///         CMetaRegistry::SEntry entry
///             = CMetaRegistry::Load("ncbi", CMetaRegistry::eName_RcOrIni,
///                                   0, flags, m_SysRegistry.GetPointer());
///     ...
///     } catch (CRegistryException& e) {
///         ERR_POST_X(6, Critical << "CNcbiRegistry: "
///                       "Syntax error in system-wide configuration file: "
///                       << e.what());
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:200-209 (quoted at
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

/// The NCBI application settings that LOSAT reproduces (`check_ncbi_application_settings`).
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
/// ```c++
/// NCBI_PARAM_DECL(bool, SeqId, PreferAccessionOverGi);
/// NCBI_PARAM_DEF_EX(bool, SeqId, PreferAccessionOverGi, false, eParam_NoThread,
///     SEQ_ID_PREFER_ACCESSION_OVER_GI);
/// typedef NCBI_PARAM_TYPE(SeqId, PreferAccessionOverGi) TPreferAccessionOverGi;
///
/// NCBI_PARAM_DECL(bool, SeqId, AvoidGi);
/// NCBI_PARAM_DEF_EX(bool, SeqId, AvoidGi, false, eParam_NoThread,
///     SEQ_ID_AVOID_GI);
/// typedef NCBI_PARAM_TYPE(SeqId, AvoidGi) TAvoidGi;
///
///
/// bool CSeq_id::PreferAccessionOverGi(void)
/// {
///     return TPreferAccessionOverGi::GetDefault() || AvoidGi();
/// }
/// ```
/// (The registry entries `[SeqId] PreferAccessionOverGi` and `AvoidGi` and the variables
/// `NCBI_CONFIG__SEQID__...` set the same parameters; LOSAT rejects both as entries it does not
/// read.)
pub const SEQ_ID_RANK_VARIABLES: [&str; 2] = ["SEQ_ID_PREFER_ACCESSION_OVER_GI", "SEQ_ID_AVOID_GI"];

impl Default for ApplicationSettings {
    /// No registry entry and no variable: both data loaders on, no `[BLAST] BLASTDB`.
    ///
    /// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_scope_src.hpp:73,81
    /// ```c++
    ///         eDefault = (eUseBlastDbDataLoader | eUseGenbankDataLoader)
    /// ...
    ///     SDataLoaderConfig(bool load_proteins, EConfigOpts options = eDefault)
    /// ```
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:381-385 (the score
    /// block's name of a protein matrix is upper-cased)
    /// ```c
    ///         sbp->read_in_matrix = TRUE;
    ///         BLAST_ScoreSetAmbigRes(sbp, 'X');
    ///         sbp->name = BLAST_StrToUpper(scoring_options->matrix);
    ///     }
    ///     status = Blast_ScoreBlkMatrixFill(sbp, get_path);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1631-1652
    /// ```c
    ///     if (matrix_found == FALSE && sbp->read_in_matrix && get_path) {
    ///             // -RMH- Changed to FALSE
    ///             char* matrix_path = get_path(sbp->name, FALSE);
    ///             if (matrix_path) {
    /// ...
    ///                 if ( (fp=fopen(full_matrix_path, "r")) == NULL) {
    ///                    return -1;
    ///                 }
    /// ```
    /// (`full_matrix_path` is the path found: `s_GetCStringOfMatrixPath` cuts the name off it,
    /// blast_setup_cxx.cpp:1288-1298.)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1381
    /// (`BlastFindMatrixPath`, `get_path` above: the name and its upper case, which are the
    /// same here; `is_prot` is FALSE, so the subdirectory is "nt")
    /// ```c++
    ///        string mtx(matrix_name);
    ///        mtx = NStr::ToUpper(mtx);
    ///
    ///        // Try all the default directories
    ///        string full_path = g_FindDataFile(mtx);
    /// ...
    ///        const string& blastmat_env = app->GetEnvironment().Get("BLASTMAT");
    ///        if (CDir(blastmat_env).Exists()) {
    ///            full_path = blastmat_env;
    ///            full_path += CFile::GetPathSeparator();
    ///            full_path += mtx;
    ///            if (CFile(full_path).Exists()) {
    /// ...
    ///            full_path = blastmat_env;
    ///            full_path += CFile::GetPathSeparator();
    ///            full_path += is_prot ? "aa" : "nt";
    ///            full_path += CFile::GetPathSeparator();
    ///            full_path += mtx;
    /// ...
    ///        // Try local "data" directory
    ///        full_path = "data";
    ///        full_path += CFile::GetPathSeparator();
    ///        full_path += mtx;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/util/util_misc.cpp:131-179 (`g_FindDataFile`; a
    /// parameter's variable, when set, is its value, an empty one too, ncbi_param.cpp:400-437)
    /// ```c++
    /// NCBI_PARAM_DEF_EX(string, NCBI, DataPath, "", 0, NCBI_DATA_PATH);
    /// ...
    /// NCBI_PARAM_DEF (string, NCBI, Data, "");
    /// ...
    ///     if (CDirEntry::IsAbsolutePath(name)) {
    ///         dirs.push_back(kEmptyStr);
    ///     } else {
    ///         TNCBIDataPath path;
    ///         TNCBIDataDir dir;
    ///
    ///         if ( !path.Get().empty() ) {
    ///             NStr::Split(path.Get(), kDelim, dirs,
    ///                     NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
    ///         }
    ///         if ( !dir.Get().empty() ) {
    ///             dirs.push_back(dir.Get());
    ///         }
    ///     }
    /// ...
    ///         candidate.Reset(CDirEntry::MakePath(*dir, name));
    ///         if (candidate.Exists() &&  candidate.GetType(fl) == type) {
    ///             return candidate.GetPath();
    /// ```
    pub fn matrix_file(&self, name: &str) -> Result<Option<PathBuf>, String> {
        Ok(
            find_matrix_file(name, &self.matrix_search, &Path::is_file, &Path::is_dir)?
                .filter(|file| std::fs::File::open(file).is_ok()),
        )
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
/// ```c++
///     CNcbiApplication* app = CNcbiApplication::Instance();
///     if (app) {
///         const CNcbiRegistry& registry = app->GetConfig();
///         const string& value = registry.Get("BLAST", "MAX_SEQID_LENGTH");
///         if (!value.empty()) {
///             try { m_FastaReader->SetMaxIDLength(NStr::StringToUInt(value)); }
///             catch (const exception&) {}  // Ignore errors
///         }
///     }
/// ```
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
/// ```c++
///             const CNcbiRegistry& registry = app->GetConfig();
///             if (registry.HasEntry("BLAST", "BLASTDB")) {
///                 path += CDirEntry::NormalizePath(registry.Get("BLAST", "BLASTDB"),eFollowLinks);
///                 path += splitter;
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1290-1304
/// ```c++
/// const string& CCompoundRegistry::x_Get(const string& section,
///                                        const string& name,
///                                        TFlags flags) const
/// {
///     CConstRef<IRegistry> reg = FindByContents(section, name,
///                                               flags & ~fJustCore);
///     return reg ? reg->Get(section, name, flags & ~fJustCore) : kEmptyStr;
/// }
///
///
/// bool CCompoundRegistry::x_HasEntry(const string& section, const string& name,
///                                    TFlags flags) const
/// {
///     return FindByContents(section, name, flags).NotEmpty();
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1239-1244
/// ```c++
///     TFlags has_entry_flags = (flags | fCountCleared) & ~fJustCore;
///     REVERSE_ITERATE(TPriorityMap, it, m_PriorityMap) {
///         if (it->second->HasEntry(section, entry, has_entry_flags)) {
///             return it->second;
///         }
///     }
/// ```
/// `FindByContents` asks every layer with `fCountCleared`, under which a registry file's
/// entry with an empty value counts (`CMemoryRegistry::x_HasEntry`, ncbireg.cpp:984-991,
/// quoted at `data_loaders_of`); the file reader stores such an entry
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
}

impl DiagTarget {
    /// Writes `text` where the diagnostics go (`stderr` for `cerr`).
    pub fn write_text(&mut self, text: &[u8], stderr: &mut dyn Write) {
        match self {
            DiagTarget::File(file) => {
                let _ = file.write_all(text);
            }
            DiagTarget::Stderr => {
                let _ = stderr.write_all(text);
                let _ = stderr.flush();
            }
            DiagTarget::Discard => {}
        }
    }
}

/// No NCBI counterpart: a `DiagTarget` as a stream, each write going where the diagnostics go
/// (`DiagTarget::write_text`; `stderr` stands for `cerr`).
pub struct DiagStream<'a> {
    pub target: &'a mut DiagTarget,
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
// ```c++
//     // Check if a non-existing dir is specified
//     if ( !s_IsSpecialLogName(file_name) ) {
//         string dir = CFile(file_name).GetDir();
//         if ( !dir.empty()  &&  !CDir(dir).Exists() ) {
//             return false;
//         }
//     }
//     ...
//         // Check special filenames
//         if ( file_name.empty()  ||  file_name == "/dev/null" ) {
//             // no output
//             SetDiagStream(0, quick_flush, 0, 0, kLogName_None);
//         }
//         else if (file_name == "-") {
//             // output to stderr
//             SetDiagStream(&NcbiCerr, quick_flush, 0, 0, kLogName_Stderr);
//         }
//         else {
//             // output to file
//             unique_ptr<CFileDiagHandler> fhandler(new CFileDiagHandler());
//             if ( !fhandler->SetLogFile(file_name, eDiagFile_All, quick_flush) ) {
//                 ERR_POST_X(9, Info << "Failed to initialize log: " << file_name);
//                 return false;
//             }
// ```
/// The target of the `-logfile` value `name` (a file NCBI fails to open leaves `cerr`). The file
/// is opened (made when it does not exist) at once, as NCBI's first `SetupDiag` does, before
/// `-version` (NCBI BLAST+ 2.17.0 makeblastdb and blastdbcmd `-version -logfile F` make `F`:
/// g_reg probes `ver_log_ok`, `ver_log_ok2`).
pub fn diag_target(name: Option<&OsStr>) -> DiagTarget {
    match name {
        None => DiagTarget::Stderr,
        Some(name) if name == "-" => DiagTarget::Stderr,
        Some(name) if name.is_empty() || name == "/dev/null" => DiagTarget::Discard,
        Some(name) => match std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(name)
        {
            Ok(file) => DiagTarget::File(file),
            Err(_) => DiagTarget::Stderr,
        },
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-866
// ```c++
//     // Check if logfile is set in the args.
//     m_LogFile = 0;
//     if (!m_DisableArgDesc && argc > 1  &&  argv  &&  diag != eDS_User) {
//         for (int i = 1;  i < argc;  i++) {
//             if ( !argv[i] ) {
//                 continue;
//             }
//             if ( NStr::strcmp(argv[i], s_ArgDelimiter) == 0 ) {
//                 break;
//             }
//             if ( NStr::strcmp(argv[i], s_ArgLogFile) == 0 ) {
//                 if (!argv[++i]) {
//                     continue;
//                 }
//                 m_LogFile = argv[i];
//             } else if (NStr::StartsWith(argv[i], s_ArgLogFile)) {
//                 const char *a = argv[i] + strlen(s_ArgLogFile);
//                 if (*a == '=') {
//                     m_LogFile = ++a;
//                 }
//             }
//         }
//     }
// ```
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
    /// `-version`: print the version, exit status 0.
    Version,
    /// The words left for the argument reader, and whether there were words that LOSAT
    /// refuses although NCBI goes on reading: `-conffile FILE` (NCBI reads that registry file
    /// before the arguments) and `-dryrun` (NCBI reads the arguments, then runs nothing).
    Words {
        words: Vec<S>,
        conffile: bool,
        dryrun: bool,
    },
    /// A word LOSAT refuses where NCBI stops reading (`-version-full`, ...).
    Refused(&'static str),
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1015
// ```c++
//         for (int i = 1;  i < argc;  i++) {
//             if ( !argv[i] ) {
//                 continue;
//             }
//             if ( NStr::strcmp(argv[i], s_ArgDelimiter) == 0 ) {
//                 skip = true;
//             }
//             if (skip) {
//                 v[real_arg_index++] = argv[i];
//                 continue;
//             }
//             // Log file - ignore if diag is eDS_User - the user wants to
//             // take care about logging.
//             if ( diag != eDS_User  &&
//                 NStr::strcmp(argv[i], s_ArgLogFile) == 0 ) {
//                 if (!argv[++i]) {
//                     continue;
//                 }
//                 v[real_arg_index++] = argv[i - 1];
//                 v[real_arg_index++] = argv[i];
//                 // Configuration file
//             } else if ( NStr::strcmp(argv[i], s_ArgCfgFile) == 0 ) {
//                 if (!argv[++i]) {
//                     continue;
//                 }
//                 v[real_arg_index++] = argv[i - 1];
//                 v[real_arg_index++] = argv[i];
//                 conf_arg = argv[i];
//
//             }
//             else if (NStr::StartsWith(argv[i], s_ArgCfgFile)) {
//                 v[real_arg_index++] = argv[i];
//                 const char* a = argv[i] + strlen(s_ArgCfgFile);
//                 if (*a == '=') {
//                     conf_arg = ++a;
//                 }
//
//                 // Version
//             } else if ( NStr::strcmp(argv[i], s_ArgVersion) == 0 ) {
//                 delete[] v;
//                 // Print VERSION
//                 cout << GetFullVersion().Print( appname, CVersionAPI::fVersionInfo | CVersionAPI::fPackageShort );
//                 diag_context.DiscardMessages();
//                 return 0;
//
//                 // Full version
//             } else if ( NStr::strcmp(argv[i], s_ArgFullVersion) == 0 ) {
//                 ...
//             } else if ( NStr::strcmp(argv[i], s_ArgDryRun) == 0 ) {
//                 m_DryRun = true;
//
//                 // Save real argument
//             } else {
//                 v[real_arg_index++] = argv[i];
//             }
//         }
// ```
/// The words NCBI's application takes before its argument reader: `-version` anywhere before
/// `--` (but not as the value of `-logfile` or `-conffile`), the full versions, `-dryrun` and
/// `-conffile`. A `-logfile` or `-conffile` without a value is dropped. LOSAT refuses
/// `-conffile` (another registry file) and `-dryrun` (NCBI then runs nothing), as the search
/// programs do (line CLI DECISIONS 4), and the full versions (NCBI's build details). `words` are
/// the arguments after the program's name. This is the one reading of these words: the command
/// line's entry (`-version`, and `-conffile` and the full versions before the configuration
/// load), `cli::ncbi_version_request` and blastdbcmd read them here, and `log_file` is NCBI's
/// first scan of `-logfile`.
pub fn preparse<S: AsRef<OsStr> + Clone>(words: &[S]) -> Preparsed<S> {
    let mut kept = Vec::with_capacity(words.len());
    let mut conffile = false;
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
                    conffile |= word == b"-conffile";
                }
            }
            b"-version" => return Preparsed::Version,
            b"-version-full" => return Preparsed::Refused("-version-full"),
            b"-version-full-xml" => return Preparsed::Refused("-version-full-xml"),
            b"-version-full-json" => return Preparsed::Refused("-version-full-json"),
            b"-dryrun" => dryrun = true,
            _ => {
                conffile |= word.starts_with(b"-conffile=");
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

/// What NCBI's application start writes about the registry files, and the settings it leaves
/// (`start_application`). NCBI reads the system-wide file (`.ncbirc`) twice: when the
/// program's usage report is made (`CBlastUsageReport`, a member of every BLAST+ application,
/// before `AppMain`; `usage_report`) and when the application loads its configuration
/// (`<program>.ini`, then `.ncbirc`, after the `-logfile`, `-version` and `-conffile` words and
/// before the program's arguments; `config`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:208-210
/// ```c++
/// 	CNcbiIstrstream empty_stream(kEmptyStr);
/// 	CRef<CNcbiRegistry> registry(new CNcbiRegistry(empty_stream, IRegistry::fWithNcbirc));
/// 	if (registry->HasEntry("BLAST", "BLAST_USAGE_REPORT")) {
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (`x_TryInit`, which the
/// `-version` word ends before, ncbiapp.cpp:972-977)
/// ```c++
/// void CNcbiApplicationAPI::x_TryInit(EAppDiagStream diag, const char* conf)
/// {
///     // Load registry from the config file
///     if ( conf ) {
///         string x_conf(conf);
///         LoadConfig(*m_Config, &x_conf);
///     } else {
///         LoadConfig(*m_Config, NULL);
///     }
///     m_ConfigLoaded = true;
///
///     CDiagContext::SetupDiag(diag, m_Config, eDCM_Flush, m_LogFile);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1650-1663 (`IncludeNcbircIfAllowed`)
/// ```c++
///     try {
///         CMetaRegistry::SEntry entry
///             = CMetaRegistry::Load("ncbi", CMetaRegistry::eName_RcOrIni,
///                                   0, flags, m_SysRegistry.GetPointer());
///         ...
///     } catch (CRegistryException& e) {
///         ERR_POST_X(6, Critical << "CNcbiRegistry: "
///                       "Syntax error in system-wide configuration file: "
///                       << e.what());
///         return false;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:669-674 (`x_TryMain`: an exception of
/// the configuration load ends the application with exit status 2)
/// ```c++
///             catch (const CException& e) {
///                 NCBI_REPORT_EXCEPTION_X(15,
///                                         "Application's initialization failed", e);
///                 *got_exception = true;
///                 *exit_code = 2;
///             }
/// ```
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
    /// The first environment variable that LOSAT refuses (`rejected_variable`).
    environment_refusal: Option<String>,
    /// The usage report's messages, or LOSAT's refusal of the run there (before `-version`).
    usage_report: Result<String, String>,
    /// The configuration load's messages, or LOSAT's refusal there.
    config: Result<ConfigLoad, String>,
    /// A registry file that the start writes a message about, for LOSAT's library, which runs
    /// no application start and refuses it (`library_settings`).
    library_note: Option<String>,
    /// What `check_ncbi_application_settings` gives after the start.
    settings: Result<ApplicationSettings, String>,
}

/// The configuration load of an application start.
#[derive(Clone, Debug, PartialEq, Eq)]
enum ConfigLoad {
    /// The messages of a load (an empty text when there are none).
    Loaded(String),
    /// The report of the exception that ended the load (exit status 2).
    Failed(String),
}

/// The application start of this process: made once by the command line's entry for an NCBI
/// program (`start_application`), as NCBI's application starts once.
static APPLICATION_START: std::sync::OnceLock<ApplicationStart> = std::sync::OnceLock::new();

/// The application start of `program` in this process, made the first time and then kept: what
/// the command line's entry writes before the program's arguments are read
/// (`ApplicationStart`), and the settings that `check_ncbi_application_settings` gives after
/// it.
pub fn start_application(program: &str) -> &'static ApplicationStart {
    APPLICATION_START.get_or_init(|| {
        application_start(
            program,
            std::env::vars_os(),
            &|name: &str| std::env::var_os(name),
            &ProcessFacts::of_this_process(),
        )
    })
}

/// The application start of this process, when the command line's entry made one
/// (`start_application`).
pub fn application_start_of_process() -> Option<&'static ApplicationStart> {
    APPLICATION_START.get()
}

impl ApplicationStart {
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
        if let (DiagTarget::File(file), Ok(text)) = (target, &self.usage_report) {
            let _ = file.write_all(text.as_bytes());
        }
    }

    /// The configuration load's messages on `target` (`stderr` for `cerr`): `Ok(None)` after a
    /// load, `Ok(Some(2))` when the program's `.ini` threw (the application then exits with
    /// that status, having posted the report and then every collected message, the replay),
    /// `Err` with LOSAT's refusal. The replay goes to `target`, except in the programs whose
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
            Ok(ConfigLoad::Failed(report)) => {
                target.write_text(report.as_bytes(), stderr);
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

    /// The settings for LOSAT's library, which runs no application start: a refused variable,
    /// a refusal of the start, and a registry file that the start writes a message about are
    /// refusals there.
    fn library_settings(self) -> Result<ApplicationSettings, String> {
        if let Some(reason) = self.environment_refusal {
            return Err(reason);
        }
        self.usage_report?;
        self.config?;
        if let Some(note) = self.library_note {
            return Err(format!(
                "{note}; LOSAT's library does not run NCBI BLAST+'s application start that writes it, and this is not supported by LOSAT"
            ));
        }
        self.settings
    }
}

/// Whether `program`'s application gives the messages it collected before its diagnostics were
/// set up to `cerr` when it ends after a failed configuration load, whatever its `-logfile`:
/// those whose application holds a `CBlastAppDiagHandler` as a member (blastn, blastp, blastx,
/// tblastn, psiblast). The handler takes the default handler (`cerr`) when the application is
/// made, before `AppMain`, and puts it back when the application ends, before the collected
/// messages are posted; tblastx and rpsblast make theirs in `Run`, and makeblastdb and
/// blastdbcmd have none, so their messages go to the `-logfile` target. NCBI BLAST+ 2.17.0, the
/// official Linux x64 release (g_reg probes `probe_h`, 2026-11-12, a malformed `.ncbirc` and
/// `<program>.ini`): blastn, blastp, blastx, tblastn and psiblast with `-logfile log.txt` write
/// U U P on stderr and U P in the file, with `-logfile /dev/null` U U P on stderr; tblastx,
/// rpsblast, makeblastdb and blastdbcmd U on stderr and U P U P in the file.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:80 (likewise blastp_app.cpp:80,
/// blastx_app.cpp:80, tblastn_app.cpp:79, psiblast_app.cpp:122 `m_bah`)
/// ```c++
///     CBlastAppDiagHandler m_Bah;
/// ```
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_aux.hpp:252
/// ```c++
/// 	CBlastAppDiagHandler():m_handler(GetDiagHandler(true)), m_save (true) {}
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1181-1187
/// ```c++
/// CBlastAppDiagHandler::~CBlastAppDiagHandler()
/// {
/// 	if(m_handler) {
/// 		SetDiagHandler(m_handler);
/// 		m_handler = NULL;
/// 	}
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:94 (in `Run`; likewise
/// rpsblast_app.cpp:113)
/// ```c++
///     CBlastAppDiagHandler bah;
/// ```
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
/// ```c++
/// void CBlastUsageReport::x_CheckBlastUsageEnv()
/// {
/// 	char * blast_usage_env = getenv("BLAST_USAGE_REPORT");
/// 	if(blast_usage_env != NULL){
/// 		bool enable = NStr::StringToBool(blast_usage_env);
/// 		if (!enable) {
/// 			SetEnabled(false);
/// 			CUsageReportAPI::SetEnabled(false);
/// 			LOG_POST(Info <<"Phone home disabled");
/// 			return ;
/// 		}
/// 	}
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1634-1648 (`IncludeNcbircIfAllowed`,
/// quoted at `check_ncbi_application_settings`)
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
/// ```c++
/// class CDustMaskApplication : public CNcbiApplication
/// {
/// public:
///     /// Application constructor
///     CDustMaskApplication() {
///         CRef<CVersion> version(new CVersion());
///         version->SetVersionInfo(1, 0, 0);
///         SetFullVersion(version);
///     }
/// ```
fn has_usage_report(program: &str) -> bool {
    !matches!(program, "rpsbproc" | "dustmasker")
}

/// Whether NCBI's `program` is built apart from BLAST+ 2.17.0 (rpsbproc), so that the reports of
/// its registry files name its own toolkit build (`rpsbproc_registry_refusal`).
fn has_own_build(program: &str) -> bool {
    program == "rpsbproc"
}

/// LOSAT's refusal of a registry file that NCBI's rpsbproc reports when it loads its
/// configuration (`application_start`): its exception reports name the place of its own
/// toolkit build, which LOSAT does not write.
fn rpsbproc_registry_refusal(path: &Path) -> String {
    format!(
        "the registry file {} has an error that NCBI's rpsbproc reports when it starts; this is not supported by LOSAT's rpsbproc",
        path.display()
    )
}

/// `BLAST_USAGE_REPORT` as `NStr::StringToBool` reads it (`usage_report_reads_ncbirc`): `None`
/// when it is not set, `Some(None)` for a value that is not a Boolean (NCBI aborts).
fn usage_report_variable(env: &dyn Fn(&str) -> Option<OsString>) -> Option<Option<bool>> {
    env("BLAST_USAGE_REPORT").map(|value| ncbi_string_to_bool(&value.to_string_lossy()))
}

/// Whether the environment turns `.ncbirc` off for every registry (`IncludeNcbircIfAllowed`,
/// quoted at `check_ncbi_application_settings`): `NCBI_DONT_USE_NCBIRC`, or
/// `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` (the environment's layer of `[NCBI] DONT_USE_NCBIRC`),
/// set with any value.
fn ncbirc_turned_off(env: &dyn Fn(&str) -> Option<OsString>) -> bool {
    env("NCBI_DONT_USE_NCBIRC").is_some() || env("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC").is_some()
}

/// The messages of one read of the system-wide file (`ncbi_registry::syntax_error_notice`,
/// `ncbi_registry::error_reading_post`).
fn system_wide_posts(read: &RegistryRead) -> String {
    let mut posts = String::new();
    if read.stream_failed {
        posts.push_str(error_reading_post());
    }
    if let Some(error) = &read.error {
        posts.push_str(&syntax_error_notice(error));
    }
    posts
}

/// What a read of the registry file at `path` makes the start write, for LOSAT's library
/// (`ApplicationStart::library_note`).
fn library_note(path: &Path, read: &RegistryRead) -> Option<String> {
    if let Some(error) = &read.error {
        Some(format!(
            "the registry file {} has a syntax error that NCBI BLAST+ reports when its program starts ({})",
            path.display(),
            error.get_msg()
        ))
    } else if read.stream_failed {
        Some(format!(
            "the registry file {} is one byte that starts a byte order mark, which NCBI BLAST+ reports as an error reading it when its program starts",
            path.display()
        ))
    } else {
        None
    }
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
/// ```c++
///     string retval = WindowMaskerPathGet();
///     if ( !retval.empty() ) {
///         return retval;
///     }
///     const string kEnvVar("WINDOW_MASKER_PATH");
///     const string kSection("WINDOW_MASKER");
///     CNcbiIstrstream empty_stream(kEmptyStr);
///     CRef<CNcbiRegistry> reg(new CNcbiRegistry(empty_stream,
///                                               IRegistry::fWithNcbirc));
///     CRef<CSimpleEnvRegMapper> mapper(new CSimpleEnvRegMapper(kSection,
///                                                              kEmptyStr));
///     CRef<CEnvironmentRegistry> env_reg(new CEnvironmentRegistry);
///     env_reg->AddMapper(*mapper, CEnvironmentRegistry::ePriority_Max);
///     reg->Add(*env_reg, CNcbiRegistry::ePriority_MaxUser);
///     retval = reg->Get(kSection, kEnvVar);
///     if (retval == kEmptyStr) {
///         retval = CDir::GetCwd();
///     }
/// ```
/// NCBI reference (598d8ae6): c++/include/corelib/ncbireg.hpp:808-811,978-981
/// ```c++
///         ePriority_MaxUser  = CCompoundRegistry::ePriority_Max - 0x10000,
///         ePriority_Reserved ///< Everything greater is for internal use.
/// ...
///         ePriority_File = ePriority_Reserved,
///         ePriority_Overrides,
///         ePriority_Environment,
/// ```
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

/// Rejects the NCBI application settings that change `program`'s output: the environment
/// variables of the NCBI C++ Toolkit and the entries of the registry files that NCBI loads
/// (`<program>.ini`, then `.ncbirc` unless it is turned off). Returns the settings that
/// LOSAT reproduces: whether the data loaders are used (`[BLAST] DATA_LOADERS`, which
/// decides whether a record's first line can be a sequence identifier), the registry
/// entry `[BLAST] BLASTDB` of the database search path, and the directory of
/// `-window_masker_taxid` (`window_masker_path_of`).
///
/// After the command line's entry has started `program`'s application (`start_application`),
/// these are the settings of that start, which has written what NCBI writes about the files.
/// Without it (LOSAT's library: `api`), a registry file that NCBI's start writes a message
/// about (a syntax error) is refused as well (`ApplicationStart::library_settings`).
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
pub fn check_ncbi_application_settings(program: &str) -> Result<ApplicationSettings, String> {
    match APPLICATION_START.get() {
        Some(start) if start.program == program => start.settings.clone(),
        _ => application_settings(
            program,
            std::env::vars_os(),
            &|name: &str| std::env::var_os(name),
            &ProcessFacts::of_this_process(),
        ),
    }
}

/// `check_ncbi_application_settings` without an application start (LOSAT's library), over the
/// environment `variables` (all of them) and `env` (one of them by name), and the facts of the
/// process (`ProcessFacts`).
fn application_settings(
    program: &str,
    variables: impl Iterator<Item = (OsString, OsString)>,
    env: &dyn Fn(&str) -> Option<OsString>,
    process: &ProcessFacts,
) -> Result<ApplicationSettings, String> {
    application_start(program, variables, env, process).library_settings()
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
/// ```c++
///         TIndex::const_iterator iit
///             = m_Index.find(SKey(name, style, flags, reg_flags));
///         if (iit != m_Index.end()) {
///             _TRACE("found in cache");
///             _ASSERT(iit->second < m_Contents.size());
///             SEntry& result = m_Contents[iit->second];
///             result.Reload(flags);
///             return result;
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1247-1252, 1278-1279
/// (`LoadConfig`: `.ncbirc` without a `.ini`)
/// ```c++
///     if ( !conf ) {
///         if (reg.IncludeNcbircIfAllowed(reg_flags)) {
///     ...
///         // still consider pulling in defaults from .ncbirc
///         if (reg.IncludeNcbircIfAllowed(reg_flags)) {
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:210-211 (a
/// `[BLAST] BLAST_USAGE_REPORT` of the usage report's registry that is not a Boolean throws
/// in the constructor of the application, before `AppMain`: NCBI aborts)
/// ```c++
/// 	if (registry->HasEntry("BLAST", "BLAST_USAGE_REPORT")) {
/// 		bool enable = NStr::StringToBool(registry->Get("BLAST", "BLAST_USAGE_REPORT"));
/// ```
/// The environment layer of `HasEntry("NCBI", "DONT_USE_NCBIRC")` (`IncludeNcbircIfAllowed`,
/// quoted at `check_ncbi_application_settings`) is `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`, set with
/// any value.
fn application_start(
    program: &str,
    mut variables: impl Iterator<Item = (OsString, OsString)>,
    env: &dyn Fn(&str) -> Option<OsString>,
    process: &ProcessFacts,
) -> ApplicationStart {
    let environment_refusal = variables.find_map(|(name, _)| rejected_variable(&name));
    let mut start = ApplicationStart {
        program: program.to_string(),
        environment_refusal: environment_refusal.clone(),
        usage_report: Ok(String::new()),
        config: Ok(ConfigLoad::Loaded(String::new())),
        library_note: None,
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
            return start;
        }
    };
    let dont_use_ncbirc = ncbirc_turned_off(env);
    let ncbirc = if dont_use_ncbirc {
        None
    } else {
        find_registry(&search_path, SYSTEM_WIDE_REGISTRY)
            .map(|path| (read_registry_file(&path, Layer::SystemWide), path))
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
    if let Some((read, path)) = usage_read {
        start.library_note = library_note(path, read);
        start.usage_report = unless_refused(system_wide_posts(read));
        if start.usage_report.is_err() {
            return ended(start);
        }
    }
    // The configuration: `<program>.ini`, then `.ncbirc`.
    let ini = match find_registry(&search_path, &format!("{program}.ini")) {
        None => None,
        Some(path) => match read_registry_file(&path, Layer::Application) {
            Err(reason) => {
                start.config = Err(reason);
                return ended(start);
            }
            Ok(read) => read.map(|read| (path, read)),
        },
    };
    let mut posts = String::new();
    if let Some((path, read)) = &ini {
        start.library_note = start.library_note.take().or(library_note(path, read));
        if (read.error.is_some() || read.stream_failed) && has_own_build(program) {
            start.config = Err(rpsbproc_registry_refusal(path));
            return ended(start);
        }
        if let Some(error) = &read.error {
            let report = error.report_default_ex(106, 15, "Application's initialization failed");
            start.config = match &environment_refusal {
                Some(reason) => Err(reason.clone()),
                None => Ok(ConfigLoad::Failed(report)),
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
                    posts.push_str(&system_wide_posts(read));
                }
                if has_own_build(program) && !system_wide_posts(read).is_empty() {
                    start.config = Err(rpsbproc_registry_refusal(path));
                    return ended(start);
                }
                start.library_note = start.library_note.take().or(library_note(path, read));
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
/// ```c++
/// 	CNcbiIstrstream empty_stream(kEmptyStr);
/// 	CRef<CNcbiRegistry> registry(new CNcbiRegistry(empty_stream, IRegistry::fWithNcbirc));
/// 	if (registry->HasEntry("BLAST", "BLAST_USAGE_REPORT")) {
/// 		bool enable = NStr::StringToBool(registry->Get("BLAST", "BLAST_USAGE_REPORT"));
/// ```
fn usage_report_entry_refusal(
    env: &dyn Fn(&str) -> Option<OsString>,
    usage_read: Option<(&RegistryRead, &Path)>,
) -> Option<String> {
    let not_a_boolean = |value: &str| ncbi_string_to_bool(value).is_none();
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
/// or, after a `<program>.ini` that ends the program, the refusal of LOSAT's library.
fn ended(mut start: ApplicationStart) -> ApplicationStart {
    let reason = match (&start.usage_report, &start.config) {
        (Err(reason), _) | (_, Err(reason)) => reason.clone(),
        _ => format!(
            "{}; NCBI BLAST+ ends the program there (exit status 2), and this is not supported by LOSAT",
            start.library_note.as_deref().unwrap_or("a registry file")
        ),
    };
    start.settings = Err(reason);
    start
}

/// The settings from the environment and the entries of the registry files that the
/// application loaded (`<program>.ini`, then `.ncbirc`), after LOSAT's checks of the entries
/// (`check_registry_entries`).
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
mod tests {
    use super::*;

    #[test]
    fn toolkit_variables_are_rejected_and_harmless_entries_accepted() {
        let reject = |name: &str, _value: &str| rejected_variable(OsStr::new(name));
        assert!(reject("DIAG_POST_LEVEL", "Error").is_some());
        assert!(reject("DIAG_OLD_POST_FORMAT", "false").is_some());
        assert!(reject("NCBI_CONFIG__DIAG__POST_FILTER", "!BLAST").is_some());
        assert!(reject("NCBI_CONFIG__BLAST__LONG_SEQID", "1").is_some());
        assert!(reject("NCBI_CONFIG__NCBI__MemorySizeLimit", "1").is_some());
        assert!(reject("NCBI_CONFIG_OVERRIDES", "/tmp/x").is_some());
        assert!(reject("NCBI_CONFIG_LONG_SEQID__BLAST", "1").is_some());
        assert!(reject("NCBI_CONFIG_PATH", "/etc").is_none());
        assert!(reject("ABORT_ON_THROW", "1").is_some());
        // `application_start` refuses a `BLAST_USAGE_REPORT` that is not a Boolean.
        assert!(reject("BLAST_USAGE_REPORT", "bogus").is_none());
        assert!(reject("BLAST_USAGE_REPORT", "false").is_none());
        assert!(reject("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT", "bogus").is_none());
        assert!(reject("NCBI_CONFIG__BLAST__BLASTDB", "/db").is_none());
        assert!(reject("NCBI_CONFIG__blast__blastdb", "/db").is_none());
        assert!(reject("BLASTDB", "/db").is_none());
        assert!(reject("NCBI", "/opt/ncbi").is_none());
        assert!(reject("HOME", "/home/x").is_none());
    }

    #[test]
    fn harmless_entries_are_the_listed_ones() {
        let read = read_registry(
            b"# file comment\n; entry comment\n[BLAST]\nBLASTDB = /db \\\n  /more\n  \n[ncbi]\ndata=/data\n",
            Layer::SystemWide,
        )
        .unwrap();
        assert_eq!(read.error, None);
        assert_eq!(read.entries[0].value, b"/db \n/more");
        assert_eq!(
            check_registry_entries(Path::new("f"), &read.entries),
            Ok(())
        );
        assert!(!is_harmless_entry("BLAST", "LONG_SEQID"));
        assert!(!is_harmless_entry("DEBUG", "DIAG_POST_LEVEL"));
        assert!(!is_harmless_entry("NCBI", ".Inherits"));
        assert!(is_harmless_entry("blast", "blast_usage_report"));
        let other = read_registry(b"[BLAST]\nLONG_SEQID=1\n", Layer::SystemWide).unwrap();
        let refusal = check_registry_entries(Path::new("f"), &other.entries).unwrap_err();
        assert!(
            refusal.contains("the registry file f sets [BLAST] LONG_SEQID"),
            "{refusal}"
        );
    }

    /// No NCBI counterpart: a process with no user database entry and no program directory,
    /// for the tests whose search path does not reach them.
    fn no_home() -> Option<OsString> {
        None
    }
    fn no_dirs(_: bool) -> Result<Vec<PathBuf>, String> {
        Ok(Vec::new())
    }
    const NO_PROCESS: ProcessFacts<'static> = ProcessFacts {
        home_by_uid: &no_home,
        program_dirs: &no_dirs,
    };

    // NCBI reference (598d8ae6): c++/src/corelib/metareg.cpp:54-56 (`kConfigPathDelim` is ";" on
    // Windows and ":;" elsewhere) and :331-396 (`CMetaRegistry::GetDefaultSearchPath`, quoted at
    // `registry_search_path`); the home directory comes from `CDir::GetHome` (`ncbi_home`). The
    // test uses the delimiter and the variable of the platform it runs on, so it checks the same
    // rule everywhere.
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
        // The program's directories after the arguments are set, and before.
        let dirs = |arguments_set: bool| {
            Ok(match arguments_set {
                true => vec![PathBuf::from("/p/"), PathBuf::from("/q/")],
                false => vec![PathBuf::from("/q/")],
            })
        };
        let process = ProcessFacts {
            home_by_uid: &no_home,
            program_dirs: &dirs,
        };
        let search = |pairs| registry_search_path("blastn", &env(pairs), &process).unwrap();
        let only = search(&[("NCBI_CONFIG_PATH", "/a;/b")]);
        assert_eq!(only, vec![PathBuf::from("/a"), PathBuf::from("/b")]);
        // ':' splits the path outside Windows only (metareg.cpp:54-56).
        let colon = search(&[("NCBI_CONFIG_PATH", "/a:/b")]);
        if cfg!(windows) {
            assert_eq!(colon, vec![PathBuf::from("/a:/b")]);
        } else {
            assert_eq!(colon, vec![PathBuf::from("/a"), PathBuf::from("/b")]);
        }
        // An empty NCBI_CONFIG_PATH splits into no directory (ncbistr_util.hpp:271-272): no
        // empty element, so nothing is searched.
        assert_eq!(search(&[("NCBI_CONFIG_PATH", "")]), Vec::<PathBuf>::new());
        // No NCBI counterpart: test plumbing that picks the platform's empty-entry path and home
        // variable (the slices must be constants to be 'static).
        let spliced = search(if cfg!(windows) {
            &[("NCBI_CONFIG_PATH", "/a;;/z"), ("APPDATA", "/h")]
        } else {
            &[("NCBI_CONFIG_PATH", "/a::/z"), ("HOME", "/h")]
        });
        assert_eq!(
            &spliced[..3],
            &[PathBuf::from("/a"), PathBuf::from("."), PathBuf::from("/h")]
        );
        assert_eq!(
            &spliced[spliced.len() - 2..],
            &[PathBuf::from("/q/"), PathBuf::from("/z")]
        );
        let no_local = search(&[("NCBI_DONT_USE_LOCAL_CONFIG", ""), ("NCBI", "/n")]);
        assert_eq!(&no_local[..1], &[PathBuf::from("/n")]);
        assert_eq!(no_local.last(), Some(&PathBuf::from("/q/")));
        // A false BLAST_USAGE_REPORT, or .ncbirc turned off in the environment (an empty value
        // too), makes the path after the arguments are set (blast_usage_report.cpp:195-209,
        // ncbireg.cpp:1634-1648); a true BLAST_USAGE_REPORT, before.
        for (variable, value, set) in [
            ("BLAST_USAGE_REPORT", "false", true),
            ("BLAST_USAGE_REPORT", "0", true),
            ("BLAST_USAGE_REPORT", "on", false),
            ("NCBI_DONT_USE_NCBIRC", "", true),
            ("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC", "", true),
        ] {
            let env = |name: &str| match name {
                "HOME" => Some(OsString::from("/h")),
                _ if name == variable => Some(OsString::from(value)),
                _ => None,
            };
            let path = registry_search_path("blastn", &env, &process).unwrap();
            let last = dirs(set).unwrap();
            assert_eq!(
                &path[path.len() - last.len()..],
                &last[..],
                "{variable}={value}"
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:3586-3657 (`CDir::GetHome` and
    // `s_GetHomeByUID`, quoted at `ncbi_home`).
    #[test]
    fn the_home_directory_is_ncbis() {
        let env = |pairs: &'static [(&'static str, &'static str)]| {
            move |name: &str| {
                pairs
                    .iter()
                    .find(|(key, _)| *key == name)
                    .map(|(_, value)| OsString::from(value))
            }
        };
        let uid_home = || Some(OsString::from("/pw"));
        let home = |pairs, by_uid: &dyn Fn() -> Option<OsString>| ncbi_home(&env(pairs), by_uid);
        if cfg!(windows) {
            // APPDATA when it is set, an empty one too; then USERPROFILE.
            let both = &[("APPDATA", "/app"), ("USERPROFILE", "/user")];
            assert_eq!(home(both, &no_home), Ok(OsString::from("/app")));
            let empty = &[("APPDATA", ""), ("USERPROFILE", "/user")];
            assert_eq!(home(empty, &no_home), Ok(OsString::new()));
            assert_eq!(
                home(&[("USERPROFILE", "/user")], &no_home),
                Ok(OsString::from("/user"))
            );
            assert_eq!(home(&[("HOME", "/h")], &no_home), Ok(OsString::new()));
        } else if cfg!(unix) {
            // HOME when it is set, an empty one too; then the user ID's entry, an empty
            // directory too; without an entry the login name's, which LOSAT refuses.
            assert_eq!(home(&[("HOME", "/h")], &uid_home), Ok(OsString::from("/h")));
            assert_eq!(home(&[("HOME", "")], &uid_home), Ok(OsString::new()));
            assert_eq!(home(&[], &uid_home), Ok(OsString::from("/pw")));
            assert_eq!(home(&[], &|| Some(OsString::new())), Ok(OsString::new()));
            assert!(home(&[("USER", "u")], &no_home)
                .unwrap_err()
                .contains("by the login name"));
            assert_eq!(
                home(&[("APPDATA", "/app")], &uid_home),
                Ok(OsString::from("/pw"))
            );
        }
        // The search path leaves out an empty home (metareg.cpp:349-352).
        let process = ProcessFacts {
            home_by_uid: &|| Some(OsString::new()),
            program_dirs: &no_dirs,
        };
        let path = registry_search_path("blastn", &env(&[]), &process).unwrap();
        assert_eq!(path[0], PathBuf::from("."));
        assert!(!path.contains(&PathBuf::new()), "{path:?}");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1301-1381
    // (`BlastFindMatrixPath`) and c++/src/util/util_misc.cpp:139-179 (`g_FindDataFile`), quoted
    // at `ApplicationSettings::matrix_file`.
    #[test]
    fn matrix_files_are_found_where_ncbi_looks() {
        let sep = if cfg!(windows) { "\\" } else { "/" };
        let delimiter = if cfg!(windows) { ";" } else { ":" };
        let path = |parts: &[&str]| PathBuf::from(parts.join(sep));
        let files = [
            path(&["dp2", "MYMAT"]),
            path(&["data_dir", "MYMAT"]),
            path(&["bm", "MYMAT"]),
            path(&["bm", "nt", "NTMAT"]),
            path(&["data", "LOCAL"]),
            path(&["bm", "mixed"]),
        ];
        let dirs = [PathBuf::from("bm")];
        let is_file = |p: &Path| files.iter().any(|file| file == p);
        let is_dir = |p: &Path| dirs.iter().any(|dir| dir == p);
        let search =
            |data_path: Option<&str>, data: Option<&str>, blastmat: Option<&str>| MatrixSearch {
                data_path: data_path.map(OsString::from),
                data: Ok(data.map(OsString::from)),
                blastmat: blastmat.map(OsString::from),
            };
        let find = |name: &str, search: &MatrixSearch| {
            find_matrix_file(name, search, &is_file, &is_dir).unwrap()
        };
        let none = search(None, None, None);
        // NCBI_DATA_PATH's directories (empty parts dropped), then [NCBI] Data, each with a
        // separator added, for the upper-cased name.
        let data_path = format!("{delimiter}dp1{delimiter}{delimiter}dp2{sep}");
        assert_eq!(
            find("mymat", &search(Some(&data_path), None, None)),
            Some(path(&["dp2", "MYMAT"]))
        );
        assert_eq!(
            find("MyMat", &search(None, Some("data_dir"), Some("bm"))),
            Some(path(&["data_dir", "MYMAT"]))
        );
        // BLASTMAT, a directory: its file, then its "nt" subdirectory's.
        assert_eq!(
            find("mymat", &search(None, None, Some("bm"))),
            Some(path(&["bm", "MYMAT"]))
        );
        assert_eq!(
            find("ntmat", &search(None, None, Some("bm"))),
            Some(path(&["bm", "nt", "NTMAT"]))
        );
        assert_eq!(find("mymat", &search(None, None, Some("nobm"))), None);
        // The name is upper-cased before the search (blast_setup.c:381-385).
        assert_eq!(find("mixed", &search(None, None, Some("bm"))), None);
        // The working directory's "data".
        assert_eq!(find("local", &none), Some(path(&["data", "LOCAL"])));
        assert_eq!(find("mymat", &none), None);
        // An empty [NCBI] Data is no directory; an unreadable registry value is the reason.
        assert_eq!(find("mymat", &search(None, Some(""), None)), None);
        let unreadable = MatrixSearch {
            data: Err("reason".to_string()),
            ..none.clone()
        };
        assert_eq!(
            find_matrix_file("mymat", &unreadable, &is_file, &is_dir),
            Err("reason".to_string())
        );
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:820-971 (`NormalizePath`, quoted
    // at `normalize_path`).
    #[test]
    fn paths_are_normalized_like_ncbis() {
        let lexical =
            |path: &str| String::from_utf8(normalize_path(path.as_bytes(), None).unwrap()).unwrap();
        for (path, normal) in [
            ("", ""),
            ("/", "/"),
            ("///", "/"),
            ("/usr/bin/", "/usr/bin"),
            ("/usr//./bin", "/usr/bin"),
            ("/../a", "/a"),
            ("/a/b/../../..", "/"),
            ("/a/b/../c", "/a/c"),
            ("a/../..", ".."),
            ("../../a", "../../a"),
            ("./a", "a"),
            ("./", "."),
            (".", "."),
            ("a/.", "a"),
            ("a/..", "."),
            ("./..", ".."),
            ("bin/blastn", "bin/blastn"),
        ] {
            assert_eq!(lexical(path), normal, "{path:?}");
        }
        // With links: /l -> /real/dir, /real/dir/x -> ../y (relative to its directory).
        let read_link = |name: &[u8]| match name {
            b"/l" => Some(b"/real/dir".to_vec()),
            b"/real/dir/x" => Some(b"../y".to_vec()),
            _ => None,
        };
        let follow = |path: &str| {
            String::from_utf8(normalize_path(path.as_bytes(), Some(&read_link)).unwrap()).unwrap()
        };
        assert_eq!(follow("/l/z"), "/real/dir/z");
        assert_eq!(follow("/l/x"), "/real/y");
        assert_eq!(follow("/l/../z"), "/real/z");
        // A loop of links: NCBI stops at 1024 with a warning.
        let looping = |name: &[u8]| (name == b"ncbi").then(|| b"ncbi".to_vec());
        assert_eq!(normalize_path(b"ncbi", Some(&looping)), None);
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1452-1457, 1549-1595
    // (`FindProgramExecutablePath`, quoted at `unix_program_name`) and metareg.cpp:378-385,
    // ncbienv.cpp:408-415 (`GetProgramDirname`, quoted at `program_directories`).
    #[test]
    fn the_program_directories_are_ncbis() {
        let files = [
            "/cwd/bin/blastn",
            "/cwd/blastn",
            "/home/u/bin/blastn",
            "/opt/l/losat",
            "/cwd/ncbi",
        ];
        let is_file = |name: &[u8]| files.iter().any(|file| file.as_bytes() == name);
        // The working directory is /cwd: a relative name is tried there first.
        let in_cwd = |name: &[u8]| {
            let full = if name.starts_with(b"/") {
                name.to_vec()
            } else {
                [b"/cwd/".as_slice(), name].concat()
            };
            is_file(&normalize_path(&full, None).unwrap())
        };
        let name = |argv0: Option<&str>, path: &str| {
            String::from_utf8(unix_program_name(
                argv0.map(str::as_bytes),
                path.as_bytes(),
                b"/cwd",
                &in_cwd,
            ))
            .unwrap()
        };
        // An absolute name stays (normalized).
        assert_eq!(
            name(Some("/home/u/bin/../bin/blastn"), ""),
            "/home/u/bin/blastn"
        );
        // A relative name of a file in the working directory.
        assert_eq!(name(Some("./bin/blastn"), ""), "/cwd/bin/blastn");
        assert_eq!(name(Some("blastn"), "/home/u/bin"), "/cwd/blastn");
        // Else the first directory of PATH (empty parts dropped) with a file of the base name.
        assert_eq!(
            name(Some("bin2/blastn.x"), "::/nowhere::/home/u/bin:"),
            "/home/u/bin/blastn"
        );
        assert_eq!(name(Some("losat"), "/x:/opt/l/"), "/opt/l/losat");
        // Not found: argv[0] as it is.
        assert_eq!(name(Some("tblastn"), "/home/u/bin"), "tblastn");
        // No argv[0] or an empty one: the name "ncbi" is looked for, else none.
        assert_eq!(name(None, ""), "/cwd/ncbi");
        assert_eq!(name(Some(""), ""), "/cwd/ncbi");
        let none = |name: &[u8]| name == b"x";
        assert_eq!(unix_program_name(Some(b""), b"", b"/cwd", &none), b"");
        assert_eq!(unix_program_name(None, b"/bin", b"/cwd", &none), b"");
        // getcwd failed: "" + "/" + the name, and the name as it is when that is no file.
        let rooted = |name: &[u8]| name == b"x" || name == b"/x";
        assert_eq!(unix_program_name(Some(b"x"), b"", b"", &rooted), b"/x");
        assert_eq!(unix_program_name(Some(b"x"), b"", b"", &none), b"x");
        // The directories: the name's, then the resolved name's when it is another one (the
        // last of '/', '\\' and ':' ends a directory on every platform).
        let dirs = |name: &str, resolved: Option<&str>| {
            program_dirs_of(name.as_bytes(), resolved.map(str::as_bytes))
                .into_iter()
                .map(|dir| String::from_utf8(dir).unwrap())
                .collect::<Vec<_>>()
        };
        assert_eq!(
            dirs("/home/u/bin/blastn", Some("/opt/l/losat")),
            ["/home/u/bin/", "/opt/l/"]
        );
        assert_eq!(dirs("/opt/l/losat", Some("/opt/l/losat")), ["/opt/l/"]);
        assert_eq!(dirs("tblastn", Some("/opt/l/losat")), ["/opt/l/"]);
        assert_eq!(dirs("", None), Vec::<String>::new());
        assert_eq!(dirs("/a/b:c", None), ["/a/b:"]);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:75-92
    // (case-insensitive substrings `blastdb`, `genbank` and `none`) and
    // c++/src/corelib/ncbireg.cpp:1235-1246, 1577-1583 (the environment, then the program's
    // `.ini`, then `.ncbirc`); an empty value in a registry file counts as no entry
    // (ncbireg.cpp:984-991), while a variable set to an empty value is an entry
    // (env_reg.cpp:157-167; audit finding B-1 of session SFc).
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
        // An empty value in a registry file is no entry: the next layer decides, and with
        // none left the data loaders stay on.
        assert!(!layers(None, Some(""), Some("x")));
        assert!(layers(None, Some(""), Some("genbank")));
        assert!(layers(None, None, Some("")));
        assert!(layers(None, Some(""), Some("")));
        // An empty variable is an entry: it answers first and turns both loaders off, also
        // over a registry file that turns them on.
        assert!(!layers(Some(""), None, None));
        assert!(!layers(Some(""), Some("blastdb"), Some("genbank")));
        assert!(!layers(Some(""), None, Some("genbank")));
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1636-1651 (`.ncbirc` is read
    // unless NCBI_DONT_USE_NCBIRC is set or the registry has [NCBI] DONT_USE_NCBIRC, also
    // from the environment) and c++/src/corelib/env_reg.cpp:366-375 (the variable of
    // [BLAST] DATA_LOADERS is NCBI_CONFIG__BLAST__DATA_LOADERS).
    #[test]
    fn application_settings_read_the_environment_and_the_registry_files() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-data-loaders-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let ini = dir.join("blastn.ini");
        let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
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
                &NO_PROCESS,
            )
            .map(|settings| settings.data_loaders)
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
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS =\n").unwrap();
        assert_eq!(settings(&[]), Ok(true));
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
            Ok(false)
        );
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
        // Only the exact variable name is the entry.
        assert_eq!(
            settings(&[("NCBI_CONFIG__blast__data_loaders", "blastdb")]),
            Ok(false)
        );
        // The program's .ini comes before .ncbirc; its last entry counts, without quotes.
        std::fs::write(
            &ini,
            "[blast]\ndata_loaders = none\n[BLAST]\nDATA_LOADERS = \"GenBank\"\n",
        )
        .unwrap();
        assert_eq!(settings(&[]), Ok(true));
        std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = \"\"\n").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        // A value NCBI reports as malformed is refused by the library (no application start).
        std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = no\"ne\n").unwrap();
        assert!(settings(&[]).unwrap_err().contains("syntax error"));
        // The escapes are read: `\x6e` is `n`.
        std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = \\x6eone\n").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        // [NCBI] DONT_USE_NCBIRC in the .ini or the environment turns .ncbirc off.
        std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
        assert_eq!(settings(&[]), Ok(true));
        std::fs::remove_file(&ini).unwrap();
        assert_eq!(
            settings(&[("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC", "")]),
            Ok(true)
        );
        assert_eq!(settings(&[("NCBI_DONT_USE_NCBIRC", "1")]), Ok(true));
        std::fs::remove_file(&ncbirc).unwrap();
        std::fs::remove_dir(&dir).unwrap();
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:739-747 (the
    // registry's `[BLAST] MAX_SEQID_LENGTH`, used when it is not empty): the value of the first
    // layer that has the entry (the environment, the program's .ini, .ncbirc; an empty value
    // hides the layers below, as makeblastdb 2.17.0 reads them: pr2/fix/n4probe), for
    // makeblastdb to refuse with -parse_seqids when NCBI acts on it; the entry stays harmless
    // for the other programs.
    #[test]
    fn max_seqid_length_is_reported_from_every_layer() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-max-seqid-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let ini = dir.join("makeblastdb.ini");
        let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
        let _ = std::fs::remove_file(&ini);
        let _ = std::fs::remove_file(&ncbirc);
        let settings = |pairs: &[(&str, &str)]| {
            let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
            all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_settings(
                "makeblastdb",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
            .map(|settings| settings.max_seqid_length)
        };
        let value = |v: &str| Ok(Ok(Some(v.to_string())));
        assert_eq!(settings(&[]), Ok(Ok(None)));
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "60")]),
            value("60")
        );
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "")]),
            value("")
        );
        std::fs::write(&ncbirc, "[BLAST]\nMAX_SEQID_LENGTH = 70\n").unwrap();
        assert_eq!(settings(&[]), value("70"));
        // An empty environment value hides the files' values.
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "")]),
            value("")
        );
        std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH = 80\n").unwrap();
        assert_eq!(settings(&[]), value("80"));
        std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH =\n").unwrap();
        assert_eq!(settings(&[]), value(""));
        std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH = \\x37\\x30\n").unwrap();
        assert_eq!(settings(&[]), value("70"));
        std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH = 1\\xff\n").unwrap();
        assert!(matches!(settings(&[]), Ok(Err(_))));
        std::fs::write(&ini, "[BLAST]\nBLASTDB = /tmp\n").unwrap();
        assert_eq!(settings(&[]), value("70"));
        std::fs::write(&ncbirc, "[BLAST]\nLONG_SEQID = 1\n").unwrap();
        assert!(settings(&[]).is_err());
        std::fs::remove_file(&ini).unwrap();
        std::fs::remove_file(&ncbirc).unwrap();
        std::fs::remove_dir(&dir).unwrap();
    }

    /// A field of `registry.tsv`: backslash escapes, `\s` for a space at either end and
    /// `\e` for an empty field.
    fn unescape(field: &str) -> String {
        let mut out = String::new();
        let mut chars = field.chars();
        while let Some(c) = chars.next() {
            if c != '\\' {
                out.push(c);
                continue;
            }
            match chars.next() {
                Some('t') => out.push('\t'),
                Some('n') => out.push('\n'),
                Some('r') => out.push('\r'),
                Some('s') => out.push(' '),
                Some('e') => {}
                Some(other) => out.push(other),
                None => out.push('\\'),
            }
        }
        out
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:527-533
    // (`registry.HasEntry("BLAST", "BLASTDB")` and `registry.Get`, quoted at `blastdb_of`) and
    // c++/src/corelib/ncbireg.cpp:1235-1246, 1577-1583 (the environment, then the program's
    // `.ini`, then `.ncbirc`).
    /// `tests/seqdb/expected/search_path/registry.tsv`: the value of `[BLAST] BLASTDB` that
    /// NCBI's blastp put into the database search path (it prints the path when a database
    /// is not found) for a variable, a `blastp.ini` and a `.ncbirc`; `{R}` stands for the
    /// directory of the run, `<unset>` for a variable or file that is absent and `<none>`
    /// for a registry without the entry.
    #[test]
    fn the_registry_value_of_blastdb_is_ncbis() {
        let table = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/seqdb/expected/search_path/registry.tsv");
        let text = std::fs::read_to_string(&table).unwrap();
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-blastdb-{}",
            std::process::id()
        ));
        let place = |value: &str| value.replace("{R}", "/sandbox");
        let mut cases = 0;
        for line in text.lines().skip(1) {
            let row: Vec<String> = line.split('\t').map(unescape).collect();
            assert_eq!(row.len(), 8, "{line:?}");
            let (id, variable, value, ini, ncbirc, registry) =
                (&row[0], &row[1], &row[2], &row[3], &row[4], &row[5]);
            std::fs::create_dir_all(&dir).unwrap();
            for (name, content) in [("blastp.ini", ini), (SYSTEM_WIDE_REGISTRY, ncbirc)] {
                let file = dir.join(name);
                let _ = std::fs::remove_file(&file);
                if content != "<unset>" {
                    std::fs::write(&file, place(content)).unwrap();
                }
            }
            let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
            if value != "<unset>" {
                all.push((variable.clone(), place(value)));
            }
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            let settings = application_settings(
                "blastp",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
            .unwrap_or_else(|reason| panic!("{id}: {reason}"));
            let expected = (registry != "<none>").then(|| place(registry));
            assert_eq!(settings.blastdb, Ok(expected), "{id}");
            cases += 1;
        }
        assert!(cases >= 15, "{cases} cases");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// `tests/winmask/expected/taxid_path.tsv`: the directory whose `7227/wmasker.obinary`
    /// NCBI's blastn `-window_masker_taxid 7227` used (told apart by the masks of three
    /// different counts files) for the two variables, another variable, a `blastn.ini` and
    /// a `.ncbirc` in the current directory and in the home directory; `{R}` stands for the
    /// directory of the run and `<cwd>` for the current directory. The search path here is
    /// the current directory, then the home directory, as NCBI's (`registry_search_path`).
    #[test]
    fn the_window_masker_directory_is_ncbis() {
        let table =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/winmask/expected/taxid_path.tsv");
        let text = std::fs::read_to_string(&table).unwrap();
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-window-masker-{}",
            std::process::id()
        ));
        let (cwd, home) = (dir.join("cwd"), dir.join("home"));
        let place = |value: &str| value.replace("{R}", "/sandbox");
        let mut cases = 0;
        for line in text.lines().skip(1) {
            let row: Vec<String> = line.split('\t').map(unescape).collect();
            assert_eq!(row.len(), 9, "{line:?}");
            let id = &row[0];
            std::fs::create_dir_all(&cwd).unwrap();
            std::fs::create_dir_all(&home).unwrap();
            for (file, content) in [
                (cwd.join("blastn.ini"), &row[4]),
                (cwd.join(SYSTEM_WIDE_REGISTRY), &row[5]),
                (home.join(SYSTEM_WIDE_REGISTRY), &row[6]),
            ] {
                let _ = std::fs::remove_file(&file);
                if content != "<unset>" {
                    std::fs::write(&file, place(content)).unwrap();
                }
            }
            let mut all = vec![(
                "NCBI_CONFIG_PATH".to_string(),
                format!("{};{}", cwd.display(), home.display()),
            )];
            for (name, value) in [
                ("NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH", &row[1]),
                ("WINDOW_MASKER_PATH", &row[2]),
            ] {
                if value != "<unset>" {
                    all.push((name.to_string(), place(value)));
                }
            }
            if row[3] != "<unset>" {
                let (name, value) = row[3].split_once('=').unwrap();
                all.push((name.to_string(), place(value)));
            }
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            let settings = application_settings(
                "blastn",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
            .unwrap_or_else(|reason| panic!("{id}: {reason}"));
            let used = settings
                .window_masker_registry
                .unwrap_or_else(|reason| panic!("{id}: {reason}"))
                .path
                .filter(|value| !value.is_empty());
            let expected = (row[7] != "<cwd>").then(|| place(&row[7]));
            assert_eq!(used, expected, "{id}");
            cases += 1;
        }
        assert_eq!(cases, 20);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// The entry is accepted in every program's registry files, read by NCBI's reader
    /// (`read_registry_file`). A `.ncbirc` that NCBI reports as malformed (a syntax error, a bad
    /// value or name in any entry) gets the application start's messages and gives the entries
    /// above the line, to the registry of `-window_masker_taxid` too, which posts the error again
    /// at each load with the program's prefix, also when the program's `.ini` turned `.ncbirc`
    /// off for the application and when a variable gives the directory (NCBI reads the file
    /// before it asks for the entry). A one-byte file whose stream fails is posted by that
    /// registry's first load only where no read with its flags came before (NCBI BLAST+ 2.17.0 on
    /// 2026-11-12, line NP M3c). `[NCBI] .Inherits` and a value that is not UTF-8 are refused for
    /// `-window_masker_taxid`. LOSAT's library refuses a malformed file (`library_settings`).
    #[test]
    fn window_masker_registry_values_are_read_or_refused_as_ncbi_stores_them() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-window-masker-kept-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let config_path = dir.display().to_string();
        let start_of = |program: &str, extra: &[(&str, &str)]| {
            let mut all = vec![("NCBI_CONFIG_PATH", config_path.as_str())];
            all.extend_from_slice(extra);
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| *key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_start(
                program,
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
        };
        let registry_of = |program: &str, extra: &[(&str, &str)]| {
            start_of(program, extra)
                .settings
                .unwrap()
                .window_masker_registry
        };
        let registry = |extra: &[(&str, &str)]| registry_of("blastn", extra).unwrap();
        let path_of = |path: &str| WindowMaskerRegistry {
            path: Some(path.to_string()),
            ..WindowMaskerRegistry::default()
        };
        let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
        std::fs::write(&ncbirc, "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n").unwrap();
        for program in ["blastn", "blastp", "blastx", "tblastn", "tblastx"] {
            assert_eq!(registry_of(program, &[]), Ok(path_of("/w")), "{program}");
        }
        // NCBI's reader takes the escapes and the quotes at the ends.
        std::fs::write(
            &ncbirc,
            "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = \"/w\\\\x\"\n",
        )
        .unwrap();
        assert_eq!(registry(&[]), path_of("/w\\x"));
        // A malformed line: the application goes on with the entries above it, with NCBI's
        // message at the start, and so does the registry of `-window_masker_taxid`, which posts
        // it again at each load, with the program's prefix.
        std::fs::write(
            &ncbirc,
            "[BLAST]\nBLASTDB = /dbs\n[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\nWINDOW_MASKER_PATH = a\"b\n",
        )
        .unwrap();
        for program in ["blastn", "blastp"] {
            let start = start_of(program, &[]);
            let at_start = start.usage_report().unwrap().to_string();
            assert!(
                at_start.starts_with("Critical: (110.6) CNcbiRegistry: "),
                "{program}"
            );
            let settings = start.settings.unwrap();
            assert_eq!(settings.blastdb, Ok(Some("/dbs".to_string())), "{program}");
            let registry = settings.window_masker_registry.unwrap();
            assert_eq!(registry.path.as_deref(), Some("/w"), "{program}");
            let mut posts = Posts::new();
            posts.post(
                CRITICAL,
                at_start.replacen("(110.6) ", &format!("(110.6) [{program}] "), 1),
            );
            assert_eq!(registry.load_posts, posts, "{program}");
            assert!(!registry.posts_once, "{program}");
        }
        assert!(application_settings(
            "blastn",
            [(
                OsString::from("NCBI_CONFIG_PATH"),
                OsString::from(&config_path)
            )]
            .into_iter(),
            &|name: &str| (name == "NCBI_CONFIG_PATH").then(|| OsString::from(&config_path)),
            &NO_PROCESS,
        )
        .unwrap_err()
        .contains("LOSAT's library"));
        // A one-byte file whose stream fails, read before with the registry's flags (by the
        // usage report, or by the application without a `.ini`): its loads post nothing.
        std::fs::write(&ncbirc, b"\xef").unwrap();
        assert_eq!(registry(&[]), WindowMaskerRegistry::default());
        assert_eq!(
            registry(&[("BLAST_USAGE_REPORT", "0")]),
            WindowMaskerRegistry::default()
        );
        let ini = dir.join("blastn.ini");
        std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
        assert_eq!(registry(&[]), WindowMaskerRegistry::default());
        // Without such a read (a `.ini`, and no usage report): its first load posts it.
        let mut posts = Posts::new();
        posts.post(
            ERROR,
            "Error: (110.4) [blastn] Error reading the registry after line 1: \n",
        );
        assert_eq!(
            registry(&[("BLAST_USAGE_REPORT", "0")]),
            WindowMaskerRegistry {
                path: None,
                load_posts: posts,
                posts_once: true,
            }
        );
        // An entry set again while its value is not empty: a warning, kept as the stream that
        // failed is (NCBI BLAST+ 2.17.0 on 2026-11-12, line NP M3c review `d1`, `d2`).
        let twice = "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /u\nwindow_masker_path = /w\n";
        std::fs::write(&ncbirc, twice).unwrap();
        let mut posts = Posts::new();
        posts.post(
            WARNING,
            "Warning: (110.8) [blastn] Found multiple [WINDOW_MASKER] window_masker_path settings; using the one from line 3\n",
        );
        assert_eq!(
            registry(&[("BLAST_USAGE_REPORT", "0")]),
            WindowMaskerRegistry {
                path: Some("/w".to_string()),
                load_posts: posts.clone(),
                posts_once: true,
            }
        );
        assert_eq!(registry(&[]), path_of("/w"));
        // Above a malformed line: the warning, then the error, at each load (`d3`).
        std::fs::write(&ncbirc, format!("{twice}[my blast]\n")).unwrap();
        for extra in [&[][..], &[("BLAST_USAGE_REPORT", "0")][..]] {
            let registry = registry(extra);
            assert_eq!(registry.path.as_deref(), Some("/w"));
            assert!(!registry.posts_once);
            let mut expected = posts.clone();
            expected.post(
                CRITICAL,
                "Critical: (110.6) [blastn] CNcbiRegistry: Syntax error in system-wide configuration file: ",
            );
            assert!(
                registry.load_posts.text().starts_with(&expected.text()),
                "{}",
                String::from_utf8_lossy(&registry.load_posts.text())
            );
        }
        // An earlier entry with the empty value is no entry to NCBI's `HasEntry`.
        std::fs::write(
            &ncbirc,
            "[WINDOW_MASKER]\nWINDOW_MASKER_PATH =\nWINDOW_MASKER_PATH = /w\n",
        )
        .unwrap();
        assert_eq!(
            registry(&[("BLAST_USAGE_REPORT", "0")]).load_posts,
            Posts::new()
        );
        std::fs::write(
            &ncbirc,
            "[BLAST]\nLONG_SEQID = 1\n[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n",
        )
        .unwrap();
        assert_eq!(registry(&[]), path_of("/w"));
        assert_eq!(registry(&[("WINDOW_MASKER_PATH", "/v")]), path_of("/v"));
        // An entry before the first section is dropped by NCBI without an error, but its name
        // is checked (NCBI BLAST+ 2.17.0 on 2026-11-12: `x y = 1` there is a syntax error).
        std::fs::write(
            &ncbirc,
            "xy = 1\n[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n",
        )
        .unwrap();
        assert_eq!(registry(&[]), path_of("/w"));
        // Malformed files, with the `.ini` that turns `.ncbirc` off for the application: the
        // entries above the line, the posts at each load (also under a variable).
        let wm = "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\n";
        for (text, path) in [
            (
                "[WINDOW_MASKER\nWINDOW_MASKER_PATH = /w\n".to_string(),
                None,
            ),
            (
                "[WINDOW_MASKER]\nWINDOW_MASKER_PATH = a\"b\n".to_string(),
                None,
            ),
            (format!("[BLAST]\nLONG_SEQID = a\"b\n{wm}"), None),
            (format!("[my blast]\nx = 1\n{wm}"), None),
            (format!("[BLAST]\nlong seqid = 1\n{wm}"), None),
            (format!("x y = 1\n{wm}"), None),
            (format!("{wm}[my blast]\n"), Some("/w")),
        ] {
            std::fs::write(&ncbirc, &text).unwrap();
            for (extra, path) in [
                (&[][..], path),
                (&[("WINDOW_MASKER_PATH", "/v")][..], Some("/v")),
            ] {
                let registry = registry(extra);
                assert_eq!(registry.path.as_deref(), path, "{text:?} {extra:?}");
                assert!(
                    registry.load_posts.text().starts_with(
                        b"Critical: (110.6) [blastn] CNcbiRegistry: Syntax error in system-wide configuration file: NCBI C++ Exception:\n"
                    ),
                    "{text:?} {extra:?}: {}",
                    String::from_utf8_lossy(&registry.load_posts.text())
                );
                assert!(!registry.posts_once);
            }
        }
        std::fs::write(&ncbirc, "[NCBI]\n.Inherits = other\n").unwrap();
        for extra in [&[][..], &[("WINDOW_MASKER_PATH", "/v")][..]] {
            let reason = registry_of("blastn", extra).unwrap_err();
            assert!(reason.contains(".Inherits"), "{reason}");
        }
        // A byte that is not UTF-8 in the value.
        std::fs::write(&ncbirc, b"[WINDOW_MASKER]\nWINDOW_MASKER_PATH = /w\xff\n").unwrap();
        let reason = registry_of("blastn", &[]).unwrap_err();
        assert!(reason.contains("UTF-8"), "{reason}");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:753-781 (quoted at
    // `ncbi_registry::read_lines`): the escapes of a value are read, and a value that is not
    // UTF-8 is not a path LOSAT reads.
    /// A value of `[BLAST] BLASTDB` that LOSAT cannot read is no reason to reject a program
    /// that opens no database: it is kept as the reason for the program that does.
    #[test]
    fn an_unreadable_blastdb_value_is_kept_for_the_program_that_opens_a_database() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-blastdb-unreadable-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let settings = |pairs: &[(&str, &str)]| {
            let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
            all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_settings(
                "blastp",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
            .unwrap()
        };
        let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
        let ini = dir.join("blastp.ini");
        let _ = std::fs::remove_file(&ini);
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB = C:\\db\n").unwrap();
        assert_eq!(settings(&[]).blastdb, Ok(Some("C:db".to_string())));
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB = /d\\xffb\n").unwrap();
        let unreadable = settings(&[]);
        assert!(unreadable.data_loaders);
        let reason = unreadable.blastdb.unwrap_err();
        assert!(reason.contains("[BLAST] BLASTDB") && reason.contains("UTF-8"));
        // A layer above the file answers first.
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__BLASTDB", "/db")]).blastdb,
            Ok(Some("/db".to_string()))
        );
        std::fs::write(&ini, "[BLAST]\nBLASTDB = /ini\n").unwrap();
        assert_eq!(settings(&[]).blastdb, Ok(Some("/ini".to_string())));
        std::fs::write(&ini, "[BLAST]\nBLASTDB = /a\\\"b\n").unwrap();
        assert_eq!(settings(&[]).blastdb, Ok(Some("/a\"b".to_string())));
        assert_eq!(ApplicationSettings::default().blastdb, Ok(None));
        std::fs::remove_file(&ini).unwrap();
        std::fs::remove_file(&ncbirc).unwrap();
        std::fs::remove_dir(&dir).unwrap();
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

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:844-866 (quoted at `log_file`).
    #[test]
    fn the_log_file_is_the_last_one_before_the_delimiter() {
        let words = |line: &[&str]| line.iter().map(OsString::from).collect::<Vec<_>>();
        assert_eq!(
            log_file(&words(&["-logfile", "a", "-logfile=b"])),
            Some(OsStr::new("b"))
        );
        assert_eq!(
            log_file(&words(&["-in", "x", "-logfile", "a", "-logfile"])),
            Some(OsStr::new("a"))
        );
        assert_eq!(log_file(&words(&["-logfile"])), None);
        assert_eq!(log_file(&words(&["--", "-logfile", "a"])), None);
        assert_eq!(log_file(&words(&["-logfilex", "a"])), None);
        assert_eq!(log_file(&words(&["-logfile="])), Some(OsStr::new("")));
        assert!(matches!(diag_target(None), DiagTarget::Stderr));
        assert!(matches!(
            diag_target(Some(OsStr::new("-"))),
            DiagTarget::Stderr
        ));
        for name in ["", "/dev/null"] {
            assert!(matches!(
                diag_target(Some(OsStr::new(name))),
                DiagTarget::Discard
            ));
        }
        assert!(matches!(
            diag_target(Some(OsStr::new("/nonexistent-dir-of-losat/log.txt"))),
            DiagTarget::Stderr
        ));
    }

    /// The writes of an application start on the standard error stream and on the file of a
    /// `-logfile` for the three targets, as the command line's entry and blastdbcmd make them:
    /// the usage report's messages on stderr at once, the first `SetupDiag`, then the
    /// configuration's.
    fn writes_of(start: &ApplicationStart, target: &str) -> (String, String, Option<i32>) {
        let mut stderr = Vec::new();
        stderr.extend_from_slice(start.usage_report().unwrap().as_bytes());
        let file = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-start-{target}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&file);
        let mut diag = match target {
            "stderr" => DiagTarget::Stderr,
            "file" => diag_target(Some(file.as_os_str())),
            _ => diag_target(Some(OsStr::new("/dev/null"))),
        };
        start.write_log_setup(&mut diag);
        let status = start.write_config_load(&mut diag, &mut stderr).unwrap();
        drop(diag);
        let logged = std::fs::read_to_string(&file).unwrap_or_default();
        let _ = std::fs::remove_file(&file);
        (String::from_utf8(stderr).unwrap(), logged, status)
    }

    /// A start with the usage report's messages `u` and the configuration `config`.
    fn start_with(u: &str, config: ConfigLoad) -> ApplicationStart {
        start_of("tblastx", u, config)
    }

    /// `start_with` of `program`.
    fn start_of(program: &str, u: &str, config: ConfigLoad) -> ApplicationStart {
        ApplicationStart {
            program: program.to_string(),
            environment_refusal: None,
            usage_report: Ok(u.to_string()),
            config: Ok(config),
            library_note: None,
            settings: Ok(ApplicationSettings::default()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:3863-3874 (quoted at
    // `ApplicationStart::write_log_setup`), c++/src/corelib/ncbiapp.cpp:669-674 (quoted at
    // `ApplicationStart`) and c++/src/algo/blast/api/blast_aux.cpp:1181-1187 (quoted at
    // `replays_to_stderr`). NCBI BLAST+ 2.17.0 (g_reg probes: x_blastp, x_ini, x_iniUR,
    // x_inibothbad, o1, o2, l, h_log, inibad_null, inibad_empty; probe_h of 2026-11-12 for every
    // program): where the messages go. tblastx, rpsblast, makeblastdb and blastdbcmd replay the
    // collected messages to the target; blastn, blastp, blastx, tblastn and psiblast to stderr.
    #[test]
    fn the_start_writes_where_ncbi_writes() {
        let (u, a, p) = ("U\n", "A\n", "P\n");
        let s = |text: &str| text.to_string();
        let loaded = start_with(u, ConfigLoad::Loaded(s(a)));
        assert_eq!(writes_of(&loaded, "stderr"), (s("U\nA\n"), s(""), None));
        assert_eq!(writes_of(&loaded, "file"), (s(u), s("U\nA\n"), None));
        assert_eq!(writes_of(&loaded, "nowhere"), (s(u), s(""), None));
        let failed = start_with(u, ConfigLoad::Failed(s(p)));
        assert_eq!(
            writes_of(&failed, "stderr"),
            (s("U\nP\nU\nP\n"), s(""), Some(2))
        );
        assert_eq!(
            writes_of(&failed, "file"),
            (s(u), s("U\nP\nU\nP\n"), Some(2))
        );
        assert_eq!(writes_of(&failed, "nowhere"), (s(u), s(""), Some(2)));
        // Without the usage report's messages.
        let failed = start_with("", ConfigLoad::Failed(s(p)));
        assert_eq!(writes_of(&failed, "stderr"), (s("P\nP\n"), s(""), Some(2)));
        assert_eq!(writes_of(&failed, "file"), (s(""), s("P\nP\n"), Some(2)));
        let loaded = start_with("", ConfigLoad::Loaded(s(a)));
        assert_eq!(writes_of(&loaded, "stderr"), (s(a), s(""), None));
        assert_eq!(writes_of(&loaded, "file"), (s(""), s(a), None));
        assert_eq!(writes_of(&loaded, "nowhere"), (s(""), s(""), None));
        // The programs with a `CBlastAppDiagHandler`: the report to the target, the replay to
        // stderr.
        for program in ["blastn", "blastp", "blastx", "tblastn", "psiblast"] {
            let failed = start_of(program, u, ConfigLoad::Failed(s(p)));
            assert_eq!(
                writes_of(&failed, "stderr"),
                (s("U\nP\nU\nP\n"), s(""), Some(2))
            );
            assert_eq!(
                writes_of(&failed, "file"),
                (s("U\nU\nP\n"), s("U\nP\n"), Some(2))
            );
            assert_eq!(
                writes_of(&failed, "nowhere"),
                (s("U\nU\nP\n"), s(""), Some(2))
            );
            let failed = start_of(program, "", ConfigLoad::Failed(s(p)));
            assert_eq!(writes_of(&failed, "file"), (s(p), s(p), Some(2)));
            let loaded = start_of(program, u, ConfigLoad::Loaded(s(a)));
            assert_eq!(writes_of(&loaded, "file"), (s(u), s("U\nA\n"), None));
        }
        for program in ["tblastx", "rpsblast", "makeblastdb", "blastdbcmd"] {
            let failed = start_of(program, u, ConfigLoad::Failed(s(p)));
            assert_eq!(
                writes_of(&failed, "file"),
                (s(u), s("U\nP\nU\nP\n"), Some(2))
            );
        }
        // A refusal at the configuration.
        let mut refused = start_with(u, ConfigLoad::Loaded(s(a)));
        refused.config = Err(s("no"));
        let mut stderr = Vec::new();
        assert_eq!(
            refused.write_config_load(&mut DiagTarget::Stderr, &mut stderr),
            Err(s("no"))
        );
        assert!(stderr.is_empty());
    }

    /// NCBI's post of a `.ncbirc` with `garbage` on line `line` (the official Linux x64 build;
    /// g_reg probe `partial`).
    fn critical(line: usize) -> String {
        format!(
            "Critical: (110.6) CNcbiRegistry: Syntax error in system-wide configuration file: NCBI C++ Exception:\n    T0 \"{}/src/corelib/ncbireg.cpp\", line 722: Error: (CRegistryException::eEntry) ncbi::IRWRegistry::x_Read() - Invalid registry entry format: 'garbage' (m_Pos = {line})\n\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        )
    }

    // NCBI's rpsbproc 0.5.1 (the 2025-05-13 package and bioconda's, env parity probe of
    // 2026-11-12 on ws; quoted at `has_usage_report`): no usage report, so neither its read of
    // `.ncbirc` nor its `BLAST_USAGE_REPORT`; a file that rpsbproc reports is LOSAT's refusal.
    #[test]
    fn rpsbproc_starts_without_the_usage_report() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-rpsbproc-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
        let start = |pairs: &[(&str, &str)]| {
            let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
            all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_start(
                "rpsbproc",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
        };
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\n").unwrap();
        let read = start(&[("BLAST_USAGE_REPORT", "abc")]);
        assert_eq!(read.usage_report, Ok(String::new()));
        assert_eq!(read.config, Ok(ConfigLoad::Loaded(String::new())));
        assert_eq!(read.settings.unwrap().blastdb, Ok(Some("/db".to_string())));
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\ngarbage\n").unwrap();
        let malformed = start(&[]);
        assert_eq!(malformed.usage_report, Ok(String::new()));
        assert!(malformed.config.unwrap_err().contains("LOSAT's rpsbproc"));
        assert!(malformed.settings.unwrap_err().contains("LOSAT's rpsbproc"));
        std::fs::remove_file(&ncbirc).unwrap();
        std::fs::write(dir.join("rpsbproc.ini"), "garbage\n").unwrap();
        assert!(start(&[]).config.unwrap_err().contains("LOSAT's rpsbproc"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    // dustmasker has no usage report (dust_mask_app.hpp:47-55, quoted at `has_usage_report`):
    // `BLAST_USAGE_REPORT` is not read, and `.ncbirc` is read once, by the configuration load,
    // whose posts go to the diagnostics as blastdbcmd's do.
    #[test]
    fn dustmasker_starts_without_the_usage_report() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-dustmasker-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
        let start = |pairs: &[(&str, &str)]| {
            let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
            all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_start(
                "dustmasker",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
        };
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\n").unwrap();
        let read = start(&[("BLAST_USAGE_REPORT", "abc")]);
        assert_eq!(read.usage_report, Ok(String::new()));
        assert_eq!(read.config, Ok(ConfigLoad::Loaded(String::new())));
        assert_eq!(read.settings.unwrap().blastdb, Ok(Some("/db".to_string())));
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\ngarbage\n").unwrap();
        let malformed = start(&[]);
        assert_eq!(malformed.usage_report, Ok(String::new()));
        assert_eq!(malformed.config, Ok(ConfigLoad::Loaded(critical(3))));
        let _ = std::fs::remove_dir_all(&dir);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_usage_report.cpp:195-211,
    // c++/src/corelib/ncbiapp.cpp:1247-1279 and c++/src/corelib/metareg.cpp:200-209 (quoted at
    // `application_start`). NCBI BLAST+ 2.17.0 (g_reg probes of 2026-11-11: run.sh `partial`
    // default and BLAST_USAGE_REPORT=false, x_dontuse, x_inidont, x_ini, x_inibothbad, b1full,
    // b1ini, inib1r): what each read writes.
    #[test]
    fn the_start_reads_the_files_as_ncbi_reads_them() {
        let dir = std::env::temp_dir().join(format!(
            "losat-ncbi-environment-start-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ncbirc = dir.join(SYSTEM_WIDE_REGISTRY);
        let ini = dir.join("blastn.ini");
        let start = |pairs: &[(&str, &str)]| {
            let mut all = vec![("NCBI_CONFIG_PATH".to_string(), dir.display().to_string())];
            all.extend(pairs.iter().map(|(n, v)| (n.to_string(), v.to_string())));
            let env = |name: &str| {
                all.iter()
                    .find(|(key, _)| key == name)
                    .map(|(_, value)| OsString::from(value))
            };
            application_start(
                "blastn",
                all.iter()
                    .map(|(n, v)| (OsString::from(n), OsString::from(v))),
                &env,
                &NO_PROCESS,
            )
        };
        let loaded = |text: &str| Ok(ConfigLoad::Loaded(text.to_string()));
        let empty = start(&[]);
        assert_eq!(empty.usage_report, Ok(String::new()));
        assert_eq!(empty.config, loaded(""));
        assert_eq!(empty.library_note, None);
        // A syntax error: both reads report it, and the entries before it are used.
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB=/db\ngarbage\n").unwrap();
        let partial = start(&[]);
        assert_eq!(partial.usage_report, Ok(critical(3)));
        assert_eq!(partial.config, loaded(&critical(3)));
        assert_eq!(
            partial.settings.as_ref().unwrap().blastdb,
            Ok(Some("/db".to_string()))
        );
        assert!(partial
            .library_note
            .as_ref()
            .unwrap()
            .contains("syntax error"));
        let late = start(&[("BLAST_USAGE_REPORT", "false")]);
        assert_eq!(late.usage_report, Ok(String::new()));
        assert_eq!(late.config, loaded(&critical(3)));
        for variable in ["NCBI_DONT_USE_NCBIRC", "NCBI_CONFIG__NCBI__DONT_USE_NCBIRC"] {
            let off = start(&[(variable, "")]);
            assert_eq!(off.usage_report, Ok(String::new()));
            assert_eq!(off.config, loaded(""));
            assert_eq!(off.settings.unwrap().blastdb, Ok(None));
        }
        // `[NCBI] DONT_USE_NCBIRC` in the program's file: the usage report's read alone.
        std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC=1\n").unwrap();
        let ini_off = start(&[]);
        assert_eq!(ini_off.usage_report, Ok(critical(3)));
        assert_eq!(ini_off.config, loaded(""));
        // A malformed program's file: exit status 2 with NCBI's report, `.ncbirc` not read.
        std::fs::write(&ini, "garbage\n").unwrap();
        let report = "Error: (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\n";
        let failed = start(&[]);
        assert_eq!(failed.usage_report, Ok(critical(3)));
        assert_eq!(failed.config, Ok(ConfigLoad::Failed(report.to_string())));
        std::fs::remove_file(&ncbirc).unwrap();
        let failed = start(&[]);
        assert_eq!(failed.usage_report, Ok(String::new()));
        assert_eq!(failed.config, Ok(ConfigLoad::Failed(report.to_string())));
        assert!(application_settings(
            "blastn",
            [(
                OsString::from("NCBI_CONFIG_PATH"),
                dir.as_os_str().to_os_string()
            )]
            .into_iter(),
            &|name: &str| (name == "NCBI_CONFIG_PATH").then(|| dir.as_os_str().to_os_string()),
            &NO_PROCESS,
        )
        .unwrap_err()
        .contains("LOSAT's library"));
        std::fs::remove_file(&ini).unwrap();
        // A `.ncbirc` of one byte that starts a byte order mark: reported by each read that is
        // not kept from the usage report's.
        std::fs::write(&ncbirc, b"\xEF").unwrap();
        let one = start(&[]);
        assert_eq!(one.usage_report, Ok(error_reading_post().to_string()));
        assert_eq!(one.config, loaded(""));
        assert_eq!(
            start(&[("BLAST_USAGE_REPORT", "0")]).config,
            loaded(error_reading_post())
        );
        std::fs::write(&ini, "[BLAST]\n").unwrap();
        assert_eq!(start(&[]).config, loaded(error_reading_post()));
        // A program's file of that byte is read as empty, without a message.
        std::fs::remove_file(&ncbirc).unwrap();
        std::fs::write(&ini, b"\xFF").unwrap();
        let silent = start(&[]);
        assert_eq!(silent.config, loaded(""));
        assert_eq!(silent.library_note, None);
        std::fs::remove_file(&ini).unwrap();
        // LOSAT's refusals: a variable it refuses when there is a message to write; a UTF-16
        // file where it is read; a usage report's entry that NCBI aborts on.
        std::fs::write(&ncbirc, "garbage\n").unwrap();
        let refused = start(&[("DIAG_POST_LEVEL", "Error")]);
        assert!(refused
            .usage_report
            .unwrap_err()
            .contains("DIAG_POST_LEVEL"));
        let refused = start(&[("DIAG_POST_LEVEL", "Error"), ("BLAST_USAGE_REPORT", "off")]);
        assert_eq!(refused.usage_report, Ok(String::new()));
        assert!(refused.config.unwrap_err().contains("DIAG_POST_LEVEL"));
        std::fs::remove_file(&ncbirc).unwrap();
        let quiet = start(&[("DIAG_POST_LEVEL", "Error")]);
        assert_eq!(quiet.usage_report, Ok(String::new()));
        assert_eq!(quiet.config, loaded(""));
        assert!(quiet.settings.unwrap_err().contains("DIAG_POST_LEVEL"));
        std::fs::write(&ncbirc, b"\xFF\xFE[\0B\0]\0").unwrap();
        assert!(start(&[]).usage_report.unwrap_err().contains("UTF-16"));
        assert!(start(&[("BLAST_USAGE_REPORT", "no")])
            .config
            .unwrap_err()
            .contains("UTF-16"));
        assert_eq!(start(&[("NCBI_DONT_USE_NCBIRC", "1")]).config, loaded(""));
        std::fs::write(&ncbirc, "[BLAST]\nBLAST_USAGE_REPORT = maybe\n").unwrap();
        assert!(start(&[])
            .usage_report
            .unwrap_err()
            .contains("[BLAST] BLAST_USAGE_REPORT"));
        std::fs::write(&ncbirc, "[BLAST]\nBLAST_USAGE_REPORT = off\n").unwrap();
        assert_eq!(start(&[]).usage_report, Ok(String::new()));
        // The usage report alone reads `[BLAST] BLAST_USAGE_REPORT`: not with
        // `BLAST_USAGE_REPORT` false, nor from the program's file (probes `l1.*`).
        std::fs::write(&ncbirc, "[BLAST]\nBLAST_USAGE_REPORT = maybe\n").unwrap();
        let off = start(&[("BLAST_USAGE_REPORT", "false")]);
        assert_eq!(off.usage_report, Ok(String::new()));
        assert!(off.settings.is_ok());
        assert!(start(&[("NCBI_DONT_USE_NCBIRC", "")]).settings.is_ok());
        assert!(start(&[("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT", "yes")])
            .settings
            .is_ok());
        std::fs::write(&ini, "[NCBI]\nDONT_USE_NCBIRC = 1\n").unwrap();
        assert!(start(&[]).usage_report.is_err());
        std::fs::remove_file(&ncbirc).unwrap();
        std::fs::write(&ini, "[BLAST]\nBLAST_USAGE_REPORT = maybe\n").unwrap();
        let from_ini = start(&[]);
        assert_eq!(from_ini.usage_report, Ok(String::new()));
        assert!(from_ini.settings.is_ok());
        std::fs::remove_file(&ini).unwrap();
        // The variables: `BLAST_USAGE_REPORT` before any file (probes `m2.*`), the environment's
        // layer of the entry where the report makes its registry.
        for value in ["abc", ""] {
            let aborts = start(&[("BLAST_USAGE_REPORT", value)]);
            assert!(aborts
                .usage_report
                .unwrap_err()
                .contains("the environment variable BLAST_USAGE_REPORT"));
        }
        let layer = [("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT", "maybe")];
        assert!(start(&layer)
            .usage_report
            .unwrap_err()
            .contains("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT"));
        assert!(start(&[layer[0], ("NCBI_DONT_USE_NCBIRC", "1")])
            .usage_report
            .is_err());
        let off = start(&[layer[0], ("BLAST_USAGE_REPORT", "0")]);
        assert_eq!(off.usage_report, Ok(String::new()));
        assert!(off.settings.is_ok());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:926-1015 (quoted at `preparse`).
    // NCBI BLAST+ 2.17.0 (env parity review G, F4): `blastn -logfile -version` logs to a file
    // named `-version` and prints no version.
    #[test]
    fn preparse_takes_ncbi_s_words() {
        let words = |line: &[&str]| line.iter().map(|word| word.to_string()).collect::<Vec<_>>();
        let kept = |line: &[&str]| match preparse(&words(line)) {
            Preparsed::Words {
                words,
                conffile,
                dryrun,
            } => (words, conffile, dryrun),
            Preparsed::Version => panic!("version"),
            Preparsed::Refused(word) => panic!("refused {word}"),
        };
        assert_eq!(preparse(&words(&["-db", "-version"])), Preparsed::Version);
        // The value of `-logfile` or `-conffile` is not read as `-version`; a last `-logfile`
        // is dropped.
        assert_eq!(
            kept(&["-info", "-logfile", "-version"]).0,
            words(&["-info", "-logfile", "-version"])
        );
        assert!(kept(&["-conffile", "-version"]).1);
        assert_eq!(kept(&["-info", "-logfile"]).0, words(&["-info"]));
        assert_eq!(kept(&["--", "-version"]).0, words(&["--", "-version"]));
        assert_eq!(
            kept(&["-info", "-dryrun"]),
            (words(&["-info"]), false, true)
        );
        assert!(kept(&["-dryrun", "-conffile", "x"]).1);
        assert!(kept(&["-conffile=x"]).1);
        assert!(!kept(&["-conffilex"]).1);
        assert!(!kept(&["-info", "-conffile"]).1);
        assert_eq!(
            preparse(&words(&["-version-full", "-version"])),
            Preparsed::Refused("-version-full")
        );
        // `-version` comes before a `-conffile` NCBI would check after its scan.
        assert_eq!(
            preparse(&words(&["-conffile", "x", "-version"])),
            Preparsed::Version
        );
        // The words of the command line's entry.
        let os = [OsString::from("-logfile"), OsString::from("-version")];
        assert!(matches!(preparse(&os), Preparsed::Words { .. }));
        assert_eq!(log_file(&os), Some(OsStr::new("-version")));
    }
}
