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
/// `kConfigPathDelim` is ":;" outside Windows (metareg.cpp:54-56). The program's directory
/// is LOSAT's.
fn registry_search_path(env: &dyn Fn(&str) -> Option<OsString>) -> Vec<PathBuf> {
    let delimiters: &[char] = if cfg!(windows) { &[';'] } else { &[':', ';'] };
    let mut path: Vec<PathBuf> = Vec::new();
    let mut tail: Vec<PathBuf> = Vec::new();
    if let Some(config_path) = env("NCBI_CONFIG_PATH") {
        let config_path = config_path.to_string_lossy().into_owned();
        let parts: Vec<&str> = config_path.split(delimiters).collect();
        match parts.iter().position(|part| part.is_empty()) {
            None => return parts.iter().map(PathBuf::from).collect(),
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
        if let Some(home) = env(if cfg!(windows) { "USERPROFILE" } else { "HOME" }) {
            if !home.is_empty() {
                path.push(PathBuf::from(home));
            }
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
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            path.push(dir.to_path_buf());
            if let Ok(resolved) = exe.canonicalize() {
                if let Some(dir2) = resolved.parent() {
                    if dir2 != dir {
                        path.push(dir2.to_path_buf());
                    }
                }
            }
        }
    }
    path.extend(tail);
    path
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

/// Rejects the NCBI application settings that change `program`'s output: the environment
/// variables of the NCBI C++ Toolkit and the entries of the registry files that NCBI loads
/// (`<program>.ini`, then `.ncbirc` unless it is turned off).
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
pub fn check_ncbi_application_settings(program: &str) -> Result<(), String> {
    for (name, value) in std::env::vars_os() {
        if let Some(reason) = rejected_variable(&name, &value) {
            return Err(reason);
        }
    }
    let env = |name: &str| std::env::var_os(name);
    let search_path = registry_search_path(&env);
    let mut use_ncbirc = env("NCBI_DONT_USE_NCBIRC").is_none();
    if let Some(path) = find_registry(&search_path, &format!("{program}.ini")) {
        let entries = check_registry_file(&path)?;
        if entries.iter().any(|(section, name, _)| {
            section.eq_ignore_ascii_case("NCBI") && name.eq_ignore_ascii_case("DONT_USE_NCBIRC")
        }) {
            use_ncbirc = false;
        }
    }
    if use_ncbirc {
        if let Some(path) = find_registry(&search_path, ".ncbirc") {
            check_registry_file(&path)?;
        }
    }
    Ok(())
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
        let only = registry_search_path(&env(&[("NCBI_CONFIG_PATH", "/a:/b")]));
        assert_eq!(only, vec![PathBuf::from("/a"), PathBuf::from("/b")]);
        let spliced = registry_search_path(&env(&[("NCBI_CONFIG_PATH", "/a::/z"), ("HOME", "/h")]));
        assert_eq!(
            &spliced[..3],
            &[PathBuf::from("/a"), PathBuf::from("."), PathBuf::from("/h")]
        );
        assert_eq!(spliced.last(), Some(&PathBuf::from("/z")));
        let no_local =
            registry_search_path(&env(&[("NCBI_DONT_USE_LOCAL_CONFIG", ""), ("NCBI", "/n")]));
        assert_eq!(&no_local[..1], &[PathBuf::from("/n")]);
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
