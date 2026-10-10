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
fn registry_search_path(env: &dyn Fn(&str) -> Option<OsString>) -> Vec<PathBuf> {
    let delimiters: &[char] = if cfg!(windows) { &[';'] } else { &[':', ';'] };
    let mut path: Vec<PathBuf> = Vec::new();
    let mut tail: Vec<PathBuf> = Vec::new();
    if let Some(config_path) = env("NCBI_CONFIG_PATH") {
        let config_path = config_path.to_string_lossy().into_owned();
        if config_path.is_empty() {
            return path;
        }
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

/// The entries of a registry file as NCBI's reader (`IRWRegistry::Read`, `x_Read`) stores
/// them, in order, or the reason it does not read the file: a syntax error (NCBI stops at
/// that line) or a UTF-16 file. A UTF-8 byte-order mark is skipped. An entry before the
/// first section is read and checked but not stored (`IRWRegistry::Set` refuses the empty
/// section name without a message).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:617-623
/// ```c++
///     // Ensure that x_Read gets a stream it can handle.
///     EEncodingForm ef = GetTextEncodingForm(is, eBOM_Discard);
///     if (ef == eEncodingForm_Utf16Native  ||  ef == eEncodingForm_Utf16Foreign) {
///         CStringUTF8 s;
///         ReadIntoUtf8(is, &s, ef);
///         CNcbiIstrstream iss(s);
///         return x_Read(iss, flags, path);
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:794-808
/// ```c++
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
/// ```
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
fn registry_entries(bytes: &[u8]) -> Result<Vec<RegistryEntry>, RegistryError> {
    if bytes.starts_with(&[0xFF, 0xFE]) || bytes.starts_with(&[0xFE, 0xFF]) {
        return Err(RegistryError::Utf16);
    }
    let data = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(bytes);
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

/// The registry file that NCBI loads under `file_name` (the first one on the search path).
fn find_registry(search_path: &[PathBuf], file_name: &str) -> Option<PathBuf> {
    search_path
        .iter()
        .map(|dir| dir.join(file_name))
        .find(|path| path.is_file())
}

/// The entries of the registry file at `path`, or LOSAT's rejection of a file that NCBI
/// reports as a syntax error or that is UTF-16. On a syntax error NCBI writes a message with
/// the text of the exception (for `.ncbirc` `Critical: ... Syntax error in system-wide
/// configuration file: NCBI C++ Exception:` with the path and line of NCBI's own source file,
/// once for each reader of the file, and goes on with the entries before the line; for
/// `<program>.ini` `Error: (CRegistryException::...) ...` and exit code 2), which LOSAT does
/// not reproduce (`AUTHORITY.md` §J of `docs/evidence/losat_web_e2h/`).
///
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
fn read_registry_file(path: &Path) -> Result<Vec<RegistryEntry>, String> {
    let bytes = std::fs::read(path)
        .map_err(|error| format!("NCBI BLAST+ would read the registry file {} ({error}), which LOSAT cannot check; this is not supported by LOSAT", path.display()))?;
    registry_entries(&bytes).map_err(|error| match error {
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
) -> Result<ApplicationSettings, String> {
    for (name, value) in variables {
        if let Some(reason) = rejected_variable(&name, &value) {
            return Err(reason);
        }
    }
    let search_path = registry_search_path(env);
    let ncbirc_allowed = env("NCBI_DONT_USE_NCBIRC").is_none()
        && env("NCBI_CONFIG__NCBI__DONT_USE_NCBIRC").is_none();
    // `rejected_variable` has rejected a BLAST_USAGE_REPORT that is not a Boolean.
    let usage_report_reads_ncbirc = ncbirc_allowed
        && env("BLAST_USAGE_REPORT")
            .is_none_or(|value| ncbi_string_to_bool(&value.to_string_lossy()) != Some(false));
    let mut application_reads_ncbirc = ncbirc_allowed;
    let ini = find_registry(&search_path, &format!("{program}.ini"));
    let ini_entries = match &ini {
        Some(path) => {
            let entries = read_registry_file(path)?;
            check_registry_entries(path, &entries)?;
            if entries.iter().any(|entry| {
                entry.section.eq_ignore_ascii_case("NCBI")
                    && entry.name.eq_ignore_ascii_case("DONT_USE_NCBIRC")
            }) {
                application_reads_ncbirc = false;
            }
            entries
        }
        None => Vec::new(),
    };
    let mut ncbirc_entries = Vec::new();
    if usage_report_reads_ncbirc || application_reads_ncbirc {
        if let Some(path) = find_registry(&search_path, ".ncbirc") {
            let entries = read_registry_file(&path)?;
            if application_reads_ncbirc {
                check_registry_entries(&path, &entries)?;
                ncbirc_entries = entries;
            } else if env("NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT").is_none() {
                // Only the usage report reads the file: its last [BLAST] BLAST_USAGE_REPORT,
                // empty included (read into the usage report's own registry).
                let usage = entries.iter().rev().find(|entry| {
                    entry.section.eq_ignore_ascii_case("BLAST")
                        && entry.name.eq_ignore_ascii_case("BLAST_USAGE_REPORT")
                });
                if let Some(entry) = usage {
                    check_registry_entries(&path, std::slice::from_ref(entry))?;
                }
            }
        }
    }
    let ncbirc_read_directly = ini.is_some() || !usage_report_reads_ncbirc;
    let env_data_loaders = env("NCBI_CONFIG__BLAST__DATA_LOADERS");
    Ok(ApplicationSettings {
        data_loaders: data_loaders_of(
            env_data_loaders.as_deref().map(OsStr::as_encoded_bytes),
            [
                file_data_loaders(&ini_entries),
                file_data_loaders(&ncbirc_entries)
                    .filter(|value| ncbirc_read_directly || !value.is_empty()),
            ],
        ),
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
        assert_eq!(
            value(b"\xef\xbb\xbf[BLAST]\nDATA_LOADERS=none\n"),
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
        // A UTF-16 file is not read.
        assert_eq!(
            registry_entries(b"\xff\xfe[\0B\0]\0"),
            Err(RegistryError::Utf16)
        );
        assert_eq!(
            registry_entries(b"\xfe\xff\0[\0B\0]"),
            Err(RegistryError::Utf16)
        );
        assert_eq!(
            syntax_line(b"\xef\xbb[BLAST]\nDATA_LOADERS=none\n"),
            Some(1)
        );
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
        // A set but empty NCBI_CONFIG_PATH has no token, so no directory is searched
        // (ncbistr_util.hpp:268-273; audit finding A-3/B-4 of session SFd); ":" has two empty
        // tokens and splices in the default directories.
        assert!(registry_search_path(&env(&[("NCBI_CONFIG_PATH", ""), ("HOME", "/h")])).is_empty());
        let colon = registry_search_path(&env(&[("NCBI_CONFIG_PATH", ":"), ("HOME", "/h")]));
        assert_eq!(&colon[..2], &[PathBuf::from("."), PathBuf::from("/h")]);
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
