//! The configuration that NCBI's application layer (`CNcbiApplication`) reads before a BLAST
//! program runs: environment variables of the NCBI C++ Toolkit and the registry files
//! (`.ncbirc` and the program's `.ini`). LOSAT does not reproduce the settings that change
//! a program's output or exit status, and rejects them explicitly (plan DW-13); a registry
//! file with only entries that change no output (such as `[BLAST] BLASTDB`) is accepted.

use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

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
    // The usage report's registry looks for no `.ncbirc` either when `NCBI_DONT_USE_NCBIRC` is
    // set or when it has the entry `[NCBI] DONT_USE_NCBIRC`, whose only layer with entries there
    // is the environment's (`NCBI_CONFIG__NCBI__DONT_USE_NCBIRC`, an empty value too;
    // `IncludeNcbircIfAllowed`, ncbireg.cpp:1634-1648).
    let arguments_set = env("BLAST_USAGE_REPORT")
        .is_some_and(|value| ncbi_string_to_bool(&value.to_string_lossy()) == Some(false))
        || env("NCBI_DONT_USE_NCBIRC").is_some()
        || env("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC").is_some();
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
#[cfg(any(unix, test))]
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
/// Without `read_link`, never `None`.
#[cfg(any(unix, test))]
fn normalize_path(
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

/// The entries of a registry file, or the reason NCBI reports it as malformed.
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbireg.cpp:652-749
/// ```c
///     for (line = 1;  NcbiGetlineEOL(is, str);  ++line) {
///         try {
///             SIZE_TYPE len = str.length();
///             SIZE_TYPE beg = 0;
///
///             while (beg < len  &&  isspace((unsigned char) str[beg])) {
///                 ++beg;
///             }
///             if (beg == len) {
///             ...
///                 continue;
///             }
///
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
///             ...
///             default:  { // regular entry
///                 string name, value;
///                 if ( !NStr::SplitInTwo(str, "=", name, value) ) {
///                     NCBI_THROW2(CRegistryException, eEntry,
///                                 "Invalid registry entry format" + in_path
///                                 + ": '" + str + "'", line);
///                 }
///                 NStr::TruncateSpacesInPlace(name);
///             ...
///                 NStr::TruncateSpacesInPlace(value);
///             ...
///                 string cont;
///                 while (s_Backslashed(value, value.size())
///                        &&  NcbiGetlineEOL(is, cont)) {
///                     ++line;
///                     value[value.size() - 1] = '\n';
///                     value += NStr::TruncateSpaces(cont);
/// ```
/// LOSAT reads the entries only to decide whether to reject the file, so a line it does
/// not read as an entry, a comment or a section makes the file rejected.
fn registry_entries(text: &str) -> Result<Vec<(String, String, String)>, String> {
    let mut entries = Vec::new();
    let mut section = String::new();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        let trimmed = line.trim_start();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with(';') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('[') {
            match rest.find(']') {
                Some(end) if !rest[..end].trim().is_empty() => {
                    section = rest[..end].trim().to_string();
                }
                _ => return Err(format!("a section line it cannot read: {line:?}")),
            }
            continue;
        }
        let Some((name, value)) = trimmed.split_once('=') else {
            return Err(format!("a line it cannot read: {line:?}"));
        };
        let mut value = value.trim().to_string();
        while value.ends_with('\\') {
            let Some(next) = lines.next() else { break };
            value.pop();
            value.push('\n');
            value.push_str(next.trim());
        }
        entries.push((section.clone(), name.trim().to_string(), value));
    }
    Ok(entries)
}

/// The registry file that NCBI loads under `file_name` (the first one on the search path).
fn find_registry(search_path: &[PathBuf], file_name: &str) -> Option<PathBuf> {
    search_path
        .iter()
        .map(|dir| dir.join(file_name))
        .find(|path| path.is_file())
}

fn check_registry_file(path: &Path) -> Result<Vec<(String, String, String)>, String> {
    let bytes = std::fs::read(path)
        .map_err(|error| format!("NCBI BLAST+ would read the registry file {} ({error}), which LOSAT cannot check; this is not supported by LOSAT", path.display()))?;
    let text = String::from_utf8_lossy(&bytes);
    let entries = registry_entries(&text).map_err(|reason| {
        format!(
            "the registry file {} has {reason}, which NCBI BLAST+ reports as a syntax error; this is not supported by LOSAT",
            path.display()
        )
    })?;
    for (section, name, value) in &entries {
        if !is_harmless_entry(section, name, value) {
            return Err(format!(
                "the registry file {} sets [{section}] {name}, which LOSAT does not know to change no output of NCBI BLAST+ (LOSAT accepts only [BLAST] BLASTDB, BLASTMAT, DATA_LOADERS, BLASTDB_NUCL_DATA_LOADER, BLASTDB_PROT_DATA_LOADER, IGDATA, MAX_SEQID_LENGTH, BLAST_USAGE_REPORT with a Boolean, [NCBI] DATA and DONT_USE_NCBIRC); this is not supported by LOSAT",
                path.display()
            ));
        }
    }
    Ok(entries)
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

/// A registry file's value of an entry as NCBI stores it: an unescaped `"` at either end
/// is removed, and one elsewhere is NCBI's error. LOSAT does not port `NStr::ParseEscapes`
/// and rejects a value with a backslash (which also marks an escaped `"`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:753-781
/// ```c++
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
fn stored_registry_value(value: &str) -> Result<&str, &'static str> {
    if value.contains('\\') {
        return Err("a backslash escape, which LOSAT does not read");
    }
    let bytes = value.as_bytes();
    let mut beg = 0;
    let mut end = bytes.len();
    let mut pos = 0;
    while let Some(found) = bytes[pos..].iter().position(|&byte| byte == b'"') {
        let at = pos + found;
        if at >= end {
            break;
        }
        if at == beg {
            beg += 1;
        } else if at == end - 1 {
            end -= 1;
        } else {
            return Err("an unescaped '\"' in the middle, which NCBI BLAST+ reports as an error");
        }
        pos = at + 1;
    }
    Ok(&value[beg..end.max(beg)])
}

/// The value of `[BLAST] <name>` that a registry file stores: that of its last entry (a
/// later entry replaces an earlier one, `Set` without `fNoOverride`). `Err` holds why
/// LOSAT does not read the value as NCBI reads it (`stored_registry_value`).
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
fn file_blast_entry(
    path: &Path,
    entries: &[(String, String, String)],
    entry: &str,
) -> Result<Option<String>, String> {
    file_entry(path, entries, "BLAST", entry)
}

/// `file_blast_entry` of the entry `[section] entry`.
fn file_entry(
    path: &Path,
    entries: &[(String, String, String)],
    section: &str,
    entry: &str,
) -> Result<Option<String>, String> {
    let Some((_, _, value)) = entries.iter().rev().find(|(in_section, name, _)| {
        in_section.eq_ignore_ascii_case(section) && name.eq_ignore_ascii_case(entry)
    }) else {
        return Ok(None);
    };
    stored_registry_value(value)
        .map(|value| Some(value.to_string()))
        .map_err(|reason| {
            format!(
                "the registry file {} sets [{section}] {entry} to a value with {reason}; this is not supported by LOSAT",
                path.display()
            )
        })
}

/// A registry file's `[BLAST] MAX_SEQID_LENGTH` when the file has the entry (an empty value
/// too; `Err` for a value LOSAT cannot read as NCBI stores it).
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
    entries: &[(String, String, String)],
) -> Option<Result<String, String>> {
    file_blast_entry(path, entries, "MAX_SEQID_LENGTH").transpose()
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

/// Rejects the NCBI application settings that change `program`'s output: the environment
/// variables of the NCBI C++ Toolkit and the entries of the registry files that NCBI loads
/// (`<program>.ini`, then `.ncbirc` unless it is turned off). Returns the settings that
/// LOSAT reproduces: whether the data loaders are used (`[BLAST] DATA_LOADERS`, which
/// decides whether a record's first line can be a sequence identifier), and the registry
/// entry `[BLAST] BLASTDB` of the database search path.
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
    application_settings(
        program,
        std::env::vars_os(),
        &|name: &str| std::env::var_os(name),
        &ProcessFacts::of_this_process(),
    )
}

/// `check_ncbi_application_settings` over the environment `variables` (all of them) and
/// `env` (one of them by name), and the facts of the process (`ProcessFacts`).
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
fn application_settings(
    program: &str,
    variables: impl Iterator<Item = (OsString, OsString)>,
    env: &dyn Fn(&str) -> Option<OsString>,
    process: &ProcessFacts,
) -> Result<ApplicationSettings, String> {
    for (name, value) in variables {
        if let Some(reason) = rejected_variable(&name, &value) {
            return Err(reason);
        }
    }
    let search_path = registry_search_path(env, process)?;
    let mut use_ncbirc = env("NCBI_DONT_USE_NCBIRC").is_none()
        && env("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC").is_none();
    let mut ini_data_loaders = None;
    let mut ini_data = Ok(None);
    let mut ini_blastdb = Ok(None);
    let mut max_seqid_length = env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH")
        .map(|value| Ok(value.to_string_lossy().into_owned()));
    if let Some(path) = find_registry(&search_path, &format!("{program}.ini")) {
        let entries = check_registry_file(&path)?;
        if entries.iter().any(|(section, name, _)| {
            section.eq_ignore_ascii_case("NCBI") && name.eq_ignore_ascii_case("DONT_USE_NCBIRC")
        }) {
            use_ncbirc = false;
        }
        ini_data_loaders = file_blast_entry(&path, &entries, "DATA_LOADERS")?;
        ini_blastdb = file_blast_entry(&path, &entries, "BLASTDB");
        ini_data = file_entry(&path, &entries, "NCBI", "Data");
        if max_seqid_length.is_none() {
            max_seqid_length = max_seqid_length_entry(&path, &entries);
        }
    }
    let mut ncbirc_data_loaders = None;
    let mut ncbirc_blastdb = Ok(None);
    let mut ncbirc_data = Ok(None);
    if use_ncbirc {
        if let Some(path) = find_registry(&search_path, ".ncbirc") {
            let entries = check_registry_file(&path)?;
            ncbirc_data_loaders = file_blast_entry(&path, &entries, "DATA_LOADERS")?;
            ncbirc_blastdb = file_blast_entry(&path, &entries, "BLASTDB");
            ncbirc_data = file_entry(&path, &entries, "NCBI", "Data");
            if max_seqid_length.is_none() {
                max_seqid_length = max_seqid_length_entry(&path, &entries);
            }
        }
    }
    let env_data_loaders = env("NCBI_CONFIG__BLAST__DATA_LOADERS");
    let env_blastdb = env("NCBI_CONFIG__BLAST__BLASTDB").map(|value| {
        value.into_string().map_err(|_| {
            "the environment variable NCBI_CONFIG__BLAST__BLASTDB is not valid UTF-8; this is not supported by LOSAT".to_string()
        })
    });
    Ok(ApplicationSettings {
        data_loaders: data_loaders_of(
            env_data_loaders.as_deref().map(OsStr::as_encoded_bytes),
            [
                ini_data_loaders.as_deref().map(str::as_bytes),
                ncbirc_data_loaders.as_deref().map(str::as_bytes),
            ],
        ),
        blastdb: blastdb_of(env_blastdb, [ini_blastdb, ncbirc_blastdb]),
        seq_id_rank_variable: SEQ_ID_RANK_VARIABLES
            .into_iter()
            .find(|name| env(name).is_some()),
        max_seqid_length: max_seqid_length.transpose(),
        matrix_search: MatrixSearch {
            data_path: env("NCBI_DATA_PATH"),
            data: match env("NCBI_CONFIG__NCBI__DATA").or_else(|| env("NCBI_CONFIG__NCBI__Data")) {
                Some(value) => Ok(Some(value)),
                None => {
                    blastdb_of(None, [ini_data, ncbirc_data]).map(|data| data.map(OsString::from))
                }
            },
            blastmat: env("BLASTMAT"),
        },
    })
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

    #[test]
    fn registry_files_are_read_like_ncbis() {
        let entries = registry_entries(
            "# file comment\n; entry comment\n[BLAST]\nBLASTDB = /db \\\n  /more\n  \n[ncbi]\ndata=/data\n",
        )
        .unwrap();
        assert_eq!(
            entries,
            vec![
                ("BLAST".into(), "BLASTDB".into(), "/db \n/more".into()),
                ("ncbi".into(), "data".into(), "/data".into())
            ]
        );
        assert!(registry_entries("garbage line without bracket\n").is_err());
        assert!(registry_entries("[BLAST\n").is_err());
        assert!(registry_entries("[ ]\n").is_err());
        assert!(entries.iter().all(|(s, n, v)| is_harmless_entry(s, n, v)));
        assert!(!is_harmless_entry("BLAST", "LONG_SEQID", "1"));
        assert!(!is_harmless_entry("DEBUG", "DIAG_POST_LEVEL", "Error"));
        assert!(is_harmless_entry("blast", "blast_usage_report", "Off"));
        assert!(!is_harmless_entry("BLAST", "BLAST_USAGE_REPORT", "bogus"));
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
        let search = |pairs| registry_search_path(&env(pairs), &process).unwrap();
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
            let path = registry_search_path(&env, &process).unwrap();
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
        let path = registry_search_path(&env(&[]), &process).unwrap();
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

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:753-781 (an unescaped `"` at
    // either end is removed; one elsewhere is an error; escapes are not ported).
    #[test]
    fn registry_values_lose_the_quotes_at_their_ends() {
        for (raw, stored) in [
            ("none", "none"),
            ("\"none\"", "none"),
            ("\"none", "none"),
            ("none\"", "none"),
            ("\"\"", ""),
            ("\"", ""),
            ("\"\"\"", ""),
            ("", ""),
        ] {
            assert_eq!(stored_registry_value(raw), Ok(stored), "{raw:?}");
        }
        for raw in ["no\"ne", "\"a\"b\"", "none\\x", "\\\"none\""] {
            assert!(stored_registry_value(raw).is_err(), "{raw:?}");
        }
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
        std::fs::write(&ini, "[BLAST]\nDATA_LOADERS = no\"ne\n").unwrap();
        assert!(settings(&[]).unwrap_err().contains("[BLAST] DATA_LOADERS"));
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
        let ncbirc = dir.join(".ncbirc");
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
        std::fs::write(&ini, "[BLAST]\nMAX_SEQID_LENGTH = 1\\2\n").unwrap();
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
            for (name, content) in [("blastp.ini", ini), (".ncbirc", ncbirc)] {
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

    // NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:753-781 (quoted at
    // `stored_registry_value`): LOSAT does not read a value with a backslash or with a `"`
    // in the middle as NCBI reads it.
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
        let ncbirc = dir.join(".ncbirc");
        let ini = dir.join("blastp.ini");
        let _ = std::fs::remove_file(&ini);
        std::fs::write(&ncbirc, "[BLAST]\nBLASTDB = C:\\db\n").unwrap();
        let unreadable = settings(&[]);
        assert!(unreadable.data_loaders);
        let reason = unreadable.blastdb.unwrap_err();
        assert!(reason.contains("[BLAST] BLASTDB") && reason.contains("backslash"));
        // A layer above the file answers first.
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__BLASTDB", "/db")]).blastdb,
            Ok(Some("/db".to_string()))
        );
        std::fs::write(&ini, "[BLAST]\nBLASTDB = /ini\n").unwrap();
        assert_eq!(settings(&[]).blastdb, Ok(Some("/ini".to_string())));
        std::fs::write(&ini, "[BLAST]\nBLASTDB = /a\"b\n").unwrap();
        assert!(settings(&[]).blastdb.is_err());
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
}
