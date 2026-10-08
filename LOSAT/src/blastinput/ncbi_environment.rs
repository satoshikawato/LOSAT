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
/// (`NCBI_CONFIG__BLAST__DATA_LOADERS`), the program's `.ini` file, then `.ncbirc` (`None`:
/// the layer has no such entry). An empty value counts as no entry (`DECISIONS.md`
/// 2026-10-08, `AUTHORITY.md` §G1), so a lower layer decides.
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
fn data_loaders_of(layers: [Option<&[u8]>; 3]) -> bool {
    let Some(value) = layers.into_iter().flatten().find(|value| !value.is_empty()) else {
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

/// The value of `[BLAST] DATA_LOADERS` that a registry file stores: that of its last
/// entry (a later entry replaces an earlier one, `Set` without `fNoOverride`).
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
fn file_data_loaders(
    path: &Path,
    entries: &[(String, String, String)],
) -> Result<Option<String>, String> {
    let Some((_, _, value)) = entries.iter().rev().find(|(section, name, _)| {
        section.eq_ignore_ascii_case("BLAST") && name.eq_ignore_ascii_case("DATA_LOADERS")
    }) else {
        return Ok(None);
    };
    stored_registry_value(value)
        .map(|value| Some(value.to_string()))
        .map_err(|reason| {
            format!(
                "the registry file {} sets [BLAST] DATA_LOADERS to a value with {reason}; this is not supported by LOSAT",
                path.display()
            )
        })
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
pub fn check_ncbi_application_settings(program: &str) -> Result<ApplicationSettings, String> {
    application_settings(program, std::env::vars_os(), &|name: &str| {
        std::env::var_os(name)
    })
}

/// `check_ncbi_application_settings` over the environment `variables` (all of them) and
/// `env` (one of them by name).
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
) -> Result<ApplicationSettings, String> {
    for (name, value) in variables {
        if let Some(reason) = rejected_variable(&name, &value) {
            return Err(reason);
        }
    }
    let search_path = registry_search_path(env);
    let mut use_ncbirc = env("NCBI_DONT_USE_NCBIRC").is_none()
        && env("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC").is_none();
    let mut ini_data_loaders = None;
    if let Some(path) = find_registry(&search_path, &format!("{program}.ini")) {
        let entries = check_registry_file(&path)?;
        if entries.iter().any(|(section, name, _)| {
            section.eq_ignore_ascii_case("NCBI") && name.eq_ignore_ascii_case("DONT_USE_NCBIRC")
        }) {
            use_ncbirc = false;
        }
        ini_data_loaders = file_data_loaders(&path, &entries)?;
    }
    let mut ncbirc_data_loaders = None;
    if use_ncbirc {
        if let Some(path) = find_registry(&search_path, ".ncbirc") {
            let entries = check_registry_file(&path)?;
            ncbirc_data_loaders = file_data_loaders(&path, &entries)?;
        }
    }
    let env_data_loaders = env("NCBI_CONFIG__BLAST__DATA_LOADERS");
    Ok(ApplicationSettings {
        data_loaders: data_loaders_of([
            env_data_loaders.as_deref().map(OsStr::as_encoded_bytes),
            ini_data_loaders.as_deref().map(str::as_bytes),
            ncbirc_data_loaders.as_deref().map(str::as_bytes),
        ]),
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:75-92
    // (case-insensitive substrings `blastdb`, `genbank` and `none`) and
    // c++/src/corelib/ncbireg.cpp:1235-1246, 1577-1583 (the environment, then the program's
    // `.ini`, then `.ncbirc`); an empty value counts as no entry (DECISIONS.md 2026-10-08).
    #[test]
    fn data_loaders_follow_the_registry_layers_and_the_substring_rules() {
        let only = |value: &str| data_loaders_of([None, None, Some(value.as_bytes())]);
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
            data_loaders_of([
                env.map(str::as_bytes),
                ini.map(str::as_bytes),
                ncbirc.map(str::as_bytes),
            ])
        };
        assert!(layers(None, None, None));
        assert!(!layers(Some("none"), Some("blastdb"), Some("genbank")));
        assert!(layers(Some("genbank"), Some("none"), Some("none")));
        assert!(!layers(None, Some("none"), Some("blastdb")));
        assert!(layers(None, Some("blastdb"), Some("none")));
        assert!(!layers(None, None, Some("none")));
        // An empty value is no entry: the next layer decides, and with none left the
        // data loaders stay on.
        assert!(!layers(Some(""), Some("none"), None));
        assert!(!layers(Some(""), Some(""), Some("x")));
        assert!(layers(Some(""), None, Some("")));
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
            )
            .map(|settings| settings.data_loaders)
        };
        let _ = std::fs::remove_file(&ini);
        let _ = std::fs::remove_file(&ncbirc);
        assert_eq!(settings(&[]), Ok(true));
        std::fs::write(&ncbirc, "[BLAST]\nDATA_LOADERS = none\n").unwrap();
        assert_eq!(settings(&[]), Ok(false));
        // The environment comes first; an empty variable leaves the decision to the files.
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "blastdb")]),
            Ok(true)
        );
        assert_eq!(
            settings(&[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
            Ok(false)
        );
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
