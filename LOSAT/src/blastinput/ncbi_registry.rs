//! NCBI's reader of registry files (`.ncbirc`, `<program>.ini`): `IRWRegistry::Read` and
//! `IRWRegistry::x_Read` (`c++/src/corelib/ncbireg.cpp`), with the byte order mark check that
//! comes before the lines (`GetTextEncodingForm`, `ncbistre.cpp`) and the line reader
//! (`NcbiGetlineEOL`, `ncbistre.cpp`). It returns the entries that NCBI's registry keeps, in the
//! order NCBI sets them, and the `CRegistryException` that ended the read, if one did: the
//! entries before it stay in the registry, as in NCBI's. What a program does with them is
//! `ncbi_environment`'s.

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::{is_space, parse_escapes, truncate_spaces};

/// An entry that NCBI's registry keeps: its section and name (names of the symbols of
/// `s_IsNameSectionSymbol`, so ASCII; NCBI compares them without case) and its value as NCBI
/// stores it (quotes at the ends removed, escapes read).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RegistryEntry {
    pub(crate) section: String,
    pub(crate) name: String,
    pub(crate) value: Vec<u8>,
}

/// What one read of a registry file gives.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RegistryRead {
    /// The entries, in the order of the file (a later entry of the same section and name
    /// replaces an earlier one).
    pub(crate) entries: Vec<RegistryEntry>,
    /// The exception that ended the read (`CRegistryException`, thrown with the line number as
    /// its position), or `None`.
    pub(crate) error: Option<NcbiException>,
    /// The stream failed before the first line without reaching its end (a file of one byte that
    /// starts a byte order mark, `text_encoding_form`): NCBI then posts
    /// `ERROR_READING_THE_REGISTRY` (`error_reading_post`).
    pub(crate) stream_failed: bool,
    /// The entries that set a `[section] name` again while its value is not empty, in the order
    /// of the file, which NCBI warns of as it reads them (`RepeatedEntry`).
    pub(crate) repeated: Vec<RepeatedEntry>,
}

/// An entry of a registry file that sets a `[section] name` again while the value it has is not
/// empty: NCBI posts `Found multiple [section] name settings; using the one from line <line>`
/// (`repeated_entry_post`) when it reads it into an empty registry (`was_empty`, as each read of
/// a file into its layer is), with the section and the name as the entry writes them and the
/// entry's last line. `[NCBI] .Inherits` joins its values instead.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:783-799
/// ```c++
///                 if (NStr::EqualNocase(section, "NCBI")
///                     &&  NStr::EqualNocase(name, ".Inherits")
///                     &&  HasEntry(section, name, flags)) {
///                     ...
///                 } else if (was_empty  &&  HasEntry(section, name, flags)) {
///                     ERR_POST_X(8, Warning
///                                << "Found multiple [" << section << "] "
///                                << name << " settings" << in_path
///                                << "; using the one from line " << line);
///                 }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:983-990 (`HasEntry` without
/// `fCountCleared` counts an entry whose value is not empty)
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
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RepeatedEntry {
    pub(crate) section: String,
    pub(crate) name: String,
    pub(crate) line: usize,
}

/// NCBI's warning of a `RepeatedEntry`, posted while `program` runs (its diagnostics' prefix,
/// `running_syntax_error_notice`). NCBI BLAST+ 2.17.0 on 2026-11-12 (line NP M3c review, `d1`):
/// `Warning: (110.8) [blastn] Found multiple [WINDOW_MASKER] WINDOW_MASKER_PATH settings; using
/// the one from line 3`. `CMetaRegistry` reads a file without its path (`in_path` empty,
/// `error_reading_post`).
pub(crate) fn repeated_entry_post(entry: &RepeatedEntry, program: &str) -> String {
    format!(
        "Warning: (110.8) [{program}] Found multiple [{}] {} settings; using the one from line {}\n",
        entry.section, entry.name, entry.line
    )
}

impl RegistryRead {
    /// The value of `[section] name` that the read stores: that of its last entry (a later
    /// entry replaces an earlier one, `Set` without `fNoOverride`; sections and names without
    /// case), as NCBI stores it.
    ///
    /// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:794-800
    /// ```c++
    ///                 } else if (was_empty  &&  HasEntry(section, name, flags)) {
    ///                     ERR_POST_X(8, Warning
    ///                                << "Found multiple [" << section << "] "
    ///                                << name << " settings" << in_path
    ///                                << "; using the one from line " << line);
    ///                 }
    ///                 Set(section, name, value, set_flags, comment);
    /// ```
    pub(crate) fn last_value(&self, section: &str, name: &str) -> Option<&[u8]> {
        last_value(&self.entries, section, name)
    }
}

/// `RegistryRead::last_value` over `entries`.
pub(crate) fn last_value<'a>(
    entries: &'a [RegistryEntry],
    section: &str,
    name: &str,
) -> Option<&'a [u8]> {
    entries
        .iter()
        .rev()
        .find(|entry| {
            entry.section.eq_ignore_ascii_case(section) && entry.name.eq_ignore_ascii_case(name)
        })
        .map(|entry| entry.value.as_slice())
}

/// The post of `IncludeNcbircIfAllowed` for an exception of a read of the system-wide file
/// (`.ncbirc`), the whole text of one load, ending with its line end.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1650-1663
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
/// `<< e.what()` writes a C string, which ends at the first NUL byte (NCBI BLAST+ 2.17.0, g_reg
/// probe `nul`: the line ends there, and no empty line follows). The post is a `Critical` of
/// `Corelib_Reg` (110, error_codes.hpp:53), which the old post format writes as
/// `Critical: (110.6) ` before the text, and ends with a line feed (what() ends with one too,
/// so an empty line follows).
pub(crate) fn syntax_error_notice(error: &NcbiException) -> String {
    syntax_error_post(error, "")
}

/// `syntax_error_notice` posted while `program` runs, after its application has set its
/// diagnostics' prefix (`SetDiagPostPrefix`), which the old post format writes as `[program] `
/// after the error code: the post of a registry that the program loads in its search, such as
/// that of blastn's `-window_masker_taxid` (NCBI BLAST+ 2.17.0, 2026-11-12:
/// `Critical: (110.6) [blastn] CNcbiRegistry: ...`, the rest as at the start).
pub(crate) fn running_syntax_error_notice(error: &NcbiException, program: &str) -> String {
    syntax_error_post(error, &format!("[{program}] "))
}

/// The post of `syntax_error_notice` with the diagnostics' `prefix`.
fn syntax_error_post(error: &NcbiException, prefix: &str) -> String {
    let what = error.what();
    let what = what.split('\0').next().unwrap_or_default();
    format!(
        "Critical: (110.6) {prefix}CNcbiRegistry: Syntax error in system-wide configuration file: {what}\n"
    )
}

/// The registry that NCBI reads a file into, which decides how often the byte order mark is
/// checked before the lines.
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
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1687-1692 (`CNcbiRegistry::x_Read`, the
/// first read of an application's registry)
/// ```c++
///     CConstRef<IRegistry> main_reg(FindByName(sm_MainRegName));
///     if (main_reg->Empty()  &&  m_FileRegistry->Empty()) {
///         m_FileRegistry->Read(is, flags & ~fWithNcbirc);
///         LoadBaseRegistries(flags, 0, path);
///         IncludeNcbircIfAllowed(flags);
///         return NULL;
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:2080-2091 (`m_SysRegistry`'s class)
/// ```c++
/// IRWRegistry* CCompoundRWRegistry::x_Read(CNcbiIstream& in, TFlags flags,
///                                          const string& path)
/// {
///     ...
///     IRWRegistry::x_Read(in, flags, path);
///     LoadBaseRegistries(lbr_flags, 0, path);
///     return NULL;
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Layer {
    /// `.ncbirc`, read into `m_SysRegistry` (a `CCompoundRWRegistry`): `Read` checks the mark
    /// once, then its `x_Read` reads the lines.
    SystemWide,
    /// `<program>.ini`, read into the application's `CNcbiRegistry`: its `Read` checks the mark,
    /// its `x_Read` hands the stream to `m_FileRegistry->Read`, which checks it again (and
    /// returns at once on a stream that has failed) before the lines.
    Application,
}

/// The end of a line in NCBI's build for a platform (`NcbiGetlineEOL`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LineEnds {
    /// Unix (and WASI): `\n`.
    Unix,
    /// Windows: `\n`, and one `\r` before it is dropped.
    Windows,
    /// macOS (`NCBI_OS_DARWIN`): `\r`, `\n`, or `\r\n` as one end.
    Darwin,
}

impl LineEnds {
    /// The line ends of the platform LOSAT is built for.
    pub(crate) const HOST: LineEnds = if cfg!(windows) {
        LineEnds::Windows
    } else if cfg!(target_vendor = "apple") {
        LineEnds::Darwin
    } else {
        LineEnds::Unix
    };
}

/// The message of the post `ERR_POST_X(4, ...)` after a read whose stream failed before its
/// end; LOSAT's reads fail so only before the first line, with an empty line read.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:813-816
/// ```c++
///     if ( !is.eof() ) {
///         ERR_POST_X(4, "Error reading the registry after line " << line
///                    << in_path << ": " << str);
///     }
/// ```
/// `in_path` is empty: `CMetaRegistry` reads a file with `Read(ifs, flags)`, without a path
/// (metareg.cpp:90). The post is an `Error` of `Corelib_Reg` (110, error_codes.hpp:53), which the
/// old post format writes as `Error: (110.4) ` before the text (`SDiagMessage::x_OldWrite`,
/// ncbidiag.cpp:5612-5700).
pub(crate) fn error_reading_post() -> &'static str {
    "Error: (110.4) Error reading the registry after line 1: \n"
}

/// `error_reading_post` posted while `program` runs, with its diagnostics' prefix
/// (`running_syntax_error_notice`; NCBI BLAST+ 2.17.0 on 2026-11-12:
/// `Error: (110.4) [blastn] Error reading the registry after line 1: `).
pub(crate) fn running_error_reading_post(program: &str) -> String {
    error_reading_post().replacen("(110.4) ", &format!("(110.4) [{program}] "), 1)
}

/// What one `GetTextEncodingForm(is, eBOM_Discard)` leaves of the stream.
#[derive(Debug, PartialEq, Eq)]
enum Encoding<'a> {
    /// The bytes that the next read starts at (a UTF-8 byte order mark dropped).
    Text(&'a [u8]),
    /// The stream failed without reaching its end: a file of one byte, 0xEF, 0xFE or 0xFF.
    Failed,
    /// A UTF-16 byte order mark, after which NCBI converts the text to UTF-8 (`ReadIntoUtf8`).
    Utf16,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:782-826
// ```c++
// EEncodingForm GetTextEncodingForm(CNcbiIstream& input,
//                                   EBOMDiscard   discard_bom)
// {
//     EEncodingForm ef = eEncodingForm_Unknown;
//     if (input.good()) {
//         const int bom_max = 4;
//         char tmp[bom_max];
//         memset(tmp, 0, bom_max);
//         Uint2* us = reinterpret_cast<Uint2*>(tmp);
//         Uchar* uc = reinterpret_cast<Uchar*>(tmp);
//         input.get(tmp[0]);
//         int n = (int) input.gcount();
//         if (n == 1  &&  (uc[0] == 0xEF  ||  uc[0] == 0xFE  ||  uc[0] == 0xFF)){
//             input.get(tmp[1]);
//             if (input.gcount() == 1) {
//                 ++n;
//                 if (us[0] == 0xFEFF) {
//                     ef = eEncodingForm_Utf16Native;
//                 } else if (us[0] == 0xFFFE) {
//                     ef = eEncodingForm_Utf16Foreign;
//                 } else if (uc[1] == 0xBB) {
//                     input.get(tmp[2]);
//                     if (input.gcount() == 1) {
//                         ++n;
//                         if (uc[2] == 0xBF) {
//                             ef = eEncodingForm_Utf8;
//                         }
//                     }
//                 }
//             }
//         }
//         if (ef == eEncodingForm_Unknown) {
//             if (n > 1) {
//                 CStreamUtils::Pushback(input, tmp, n);
//             } else if (n == 1) {
//                 input.unget();
//             }
//         } else {
//             if (discard_bom == eBOM_Keep) {
//                 CStreamUtils::Pushback(input, tmp, n);
//             }
//         }
//     }
//     return ef;
// }
// ```
// `us[0]` is 0xFEFF or 0xFFFE for the bytes FF FE and FE FF on either byte order. The second
// `get` of a one-byte file fails at the end (`eofbit | failbit`), and `unget` of that stream
// clears `eofbit` and fails again: the stream has failed and is not at its end. NCBI 2.17.0 reads
// such a `.ncbirc` with the post "Error reading the registry after line 1" (g_reg probe, files
// `\xEF` and `\xFE`) and a `blastn.ini` of the byte `\xFF` with none (its second `Read` returns at
// the failed stream); a file of the two bytes EF BB is read as a line of them (the push-back
// leaves a stream that reads them). Any first byte followed by BB BF is taken for the UTF-8 mark,
// as the code reads (`uc[1] == 0xBB` alone).
/// One `GetTextEncodingForm(is, eBOM_Discard)` at the start of `data`.
fn text_encoding_form(data: &[u8]) -> Encoding<'_> {
    let [first, rest @ ..] = data else {
        return Encoding::Text(data);
    };
    if !matches!(first, 0xEF | 0xFE | 0xFF) {
        return Encoding::Text(data);
    }
    let Some(&second) = rest.first() else {
        return Encoding::Failed;
    };
    if [*first, second] == [0xFF, 0xFE] || [*first, second] == [0xFE, 0xFF] {
        Encoding::Utf16
    } else if second == 0xBB && data.get(2) == Some(&0xBF) {
        Encoding::Text(&data[3..])
    } else {
        Encoding::Text(data)
    }
}

/// The lines of a stream as `NcbiGetlineEOL` reads them.
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
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:90-103 (the delimiters of a string:
/// `\r\n` back to back is one end)
/// ```c++
///         SIZE_TYPE delim_pos = delims.find(CT_TO_CHAR_TYPE(ch));
///         if (delim_pos != NPOS) {
///             // Special case -- if two different delimiters are back to
///             // back and in the same order as in delims, treat them as
///             // a single delimiter (necessary for correct handling of
///             // DOS/MAC-style CR/LF endings).
///             ch = is.rdbuf()->sgetc();
///             if (!CT_EQ_INT_TYPE(ch, CT_EOF)
///                 &&  delims.find(CT_TO_CHAR_TYPE(ch), delim_pos + 1) != NPOS) {
///                 is.rdbuf()->sbumpc();
///                 delim_count = 2;
///             } else {
///                 delim_count = 1;
///             }
///             break;
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:153-176 (one delimiter: a read at the
/// end with nothing read fails, which ends the reading)
/// ```c++
///     do {
///         CT_INT_TYPE nextc = is.get();
///         if (CT_EQ_INT_TYPE(nextc, CT_EOF)
///             ||  CT_EQ_INT_TYPE(nextc, CT_TO_INT_TYPE(delim))) {
///             ++size;
///             break;
///         }
///         ...
///     } while ( is.good() );
/// #endif
///
///     if (is.rdstate() == NcbiEofbit  &&  str.empty())
///         is.setstate(NcbiFailbit);
/// ```
/// A last line without an end is read (its read reaches the end with text); the read after it,
/// like the read after a last end, fails. A NUL byte is part of a line.
struct Lines<'a> {
    data: &'a [u8],
    pos: usize,
    ends: LineEnds,
}

impl<'a> Lines<'a> {
    fn next_line(&mut self) -> Option<&'a [u8]> {
        if self.pos >= self.data.len() {
            return None;
        }
        let rest = &self.data[self.pos..];
        let found = match self.ends {
            LineEnds::Unix | LineEnds::Windows => rest.iter().position(|&c| c == b'\n'),
            LineEnds::Darwin => rest.iter().position(|&c| c == b'\r' || c == b'\n'),
        };
        let mut line = match found {
            Some(at) => {
                let pair = self.ends == LineEnds::Darwin
                    && rest[at] == b'\r'
                    && rest.get(at + 1) == Some(&b'\n');
                self.pos += at + if pair { 2 } else { 1 };
                &rest[..at]
            }
            None => {
                self.pos = self.data.len();
                rest
            }
        };
        if self.ends == LineEnds::Windows {
            if let [head @ .., b'\r'] = line {
                line = head;
            }
        }
        Some(line)
    }
}

/// `s_IsNameSectionSymbol` over a whole non-empty name (`IRegistry::IsNameSection` and
/// `IsNameEntry` without `fInternalSpaces` and `fSectionlessEntries`, which no read of LOSAT's
/// programs passes: `CMetaRegistry` reads `.ncbirc` with the flags 0 or `fJustCore`, and the
/// application's file with `fWithNcbirc`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:59-81
/// ```c++
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
fn is_name(name: &[u8]) -> bool {
    !name.is_empty()
        && name
            .iter()
            .all(|&c| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'-' | b'.' | b'/'))
}

/// `s_Backslashed(s, pos)`: an odd number of backslashes ends just before `pos`.
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
/// (With backslashes only before `pos`, `last_non_bs` is `NPOS`, and `pos - NPOS` is `pos + 1`.)
fn backslashed(s: &[u8], pos: usize) -> bool {
    s[..pos].iter().rev().take_while(|&&c| c == b'\\').count() % 2 == 1
}

/// Where `IRWRegistry::x_Read` throws `CRegistryException::<code>` at `line` of ncbireg.cpp.
const fn throw_site(line: u32, kind: &'static str) -> ThrowSite {
    ThrowSite {
        file: "src/corelib/ncbireg.cpp",
        line,
        kind,
        function: Some("ncbi::IRWRegistry::x_Read"),
    }
}

/// `NCBI_THROW2(CRegistryException, <code>, message, line)` with the line's text in `message`
/// (`before`, the text, `after`). `Err` when the text is not UTF-8: LOSAT's exception reports are
/// text, and NCBI quotes the line's bytes as they are, so LOSAT refuses such a file.
fn registry_exception(
    site: ThrowSite,
    before: &str,
    text: &[u8],
    after: &str,
    line: usize,
) -> Result<NcbiException, String> {
    let text = std::str::from_utf8(text).map_err(|_| {
        "has a syntax error on a line that is not UTF-8 (NCBI BLAST+ quotes the line's bytes in its message)".to_string()
    })?;
    Ok(NcbiException::with_pos(
        site,
        format!("{before}{text}{after}"),
        line,
    ))
}

/// Reads a registry file of `bytes` into `layer` as NCBI's build for this platform does.
/// `Err` holds why LOSAT does not read the file as NCBI does (completing "the registry file
/// NAME ..."): a UTF-16 byte order mark (NCBI converts the text first, `ReadIntoUtf8`, which
/// LOSAT does not port), or a syntax error on a line that is not UTF-8.
pub(crate) fn read_registry(bytes: &[u8], layer: Layer) -> Result<RegistryRead, String> {
    read_registry_with(bytes, layer, LineEnds::HOST)
}

/// `read_registry` with the line ends of `ends`.
fn read_registry_with(bytes: &[u8], layer: Layer, ends: LineEnds) -> Result<RegistryRead, String> {
    let checks = match layer {
        Layer::SystemWide => 1,
        Layer::Application => 2,
    };
    let mut data = bytes;
    for check in 1..=checks {
        match text_encoding_form(data) {
            Encoding::Text(rest) => data = rest,
            Encoding::Utf16 => {
                return Err(
                    "starts with a UTF-16 byte order mark (NCBI BLAST+ converts the text to UTF-8 before it reads it, which LOSAT does not do)"
                        .to_string(),
                );
            }
            // Only the last check is followed by `x_Read`; an earlier one by a `Read` that
            // returns at the failed stream.
            Encoding::Failed => {
                return Ok(RegistryRead {
                    entries: Vec::new(),
                    error: None,
                    stream_failed: check == checks,
                    repeated: Vec::new(),
                });
            }
        }
    }
    read_lines(data, ends)
}

/// `IRWRegistry::x_Read` of the lines of `data`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:652-680
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
///                 ...
///                 continue;
///             }
///
///             switch (str[beg]) {
///
///             case '#':  { // file comment
///                 SetComment(GetComment() + str + '\n');
///                 break;
///             }
///
///             case ';':  { // section or entry comment
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:694-711
/// ```c++
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
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:719-751
/// ```c++
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
///                 ...
///                 // read continuation lines, if any
///                 string cont;
///                 while (s_Backslashed(value, value.size())
///                        &&  NcbiGetlineEOL(is, cont)) {
///                     ++line;
///                     value[value.size() - 1] = '\n';
///                     value += NStr::TruncateSpaces(cont);
///                     str   += 'n' + cont; // for presentation in exceptions
///                 }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:753-801
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
///                 TFlags set_flags = flags & ~fJustCore;
///                 if (NStr::EqualNocase(section, "NCBI")
///                     &&  NStr::EqualNocase(name, ".Inherits")
///                     &&  HasEntry(section, name, flags)) {
///                     ...
///                 } else if (was_empty  &&  HasEntry(section, name, flags)) {
///                     ERR_POST_X(8, Warning
///                                << "Found multiple [" << section << "] "
///                                << name << " settings" << in_path
///                                << "; using the one from line " << line);
///                 }
///                 Set(section, name, value, set_flags, comment);
/// ```
/// `in_path` is empty (`error_reading_post`). `SplitInTwo` with no flag cuts at the first `=`. A
/// comment changes no entry. `Set` keeps no entry before the first section (`IRWRegistry::Set`,
/// ncbireg.cpp:833-837: the section "" is not a name). The warning of a repeated entry is not
/// written: a release build posts nothing below `Error` (`CDiagBuffer::sm_PostSeverity`,
/// ncbidiag.cpp:3975-3979), and NCBI 2.17.0 writes nothing for a repeated `[BLAST] BLASTDB` in
/// `.ncbirc`, in its `-logfile` neither (g_reg probes `multi`, `multi_log`). The values of
/// `[NCBI] .Inherits` that NCBI joins are kept as separate entries: no program of LOSAT reads the
/// entry, which makes NCBI read further files (`ncbi_environment` refuses it).
fn read_lines(data: &[u8], ends: LineEnds) -> Result<RegistryRead, String> {
    let mut lines = Lines { data, pos: 0, ends };
    let mut entries: Vec<RegistryEntry> = Vec::new();
    let mut repeated = Vec::new();
    let mut section = String::new();
    let mut line = 0usize;
    let done = |entries: Vec<RegistryEntry>, repeated, error: NcbiException| RegistryRead {
        entries,
        error: Some(error),
        stream_failed: false,
        repeated,
    };
    while let Some(str) = lines.next_line() {
        line += 1;
        let Some(beg) = str.iter().position(|&c| !is_space(c)) else {
            continue;
        };
        match str[beg] {
            b'#' | b';' => {}
            b'[' => {
                let beg = beg + 1;
                let end = str
                    .get(beg + 1..)
                    .and_then(|rest| rest.iter().position(|&c| c == b']'))
                    .map(|at| beg + 1 + at);
                let Some(end) = end else {
                    let error = registry_exception(
                        throw_site(698, "CRegistryException::eSection"),
                        "Invalid registry section (']' is missing): `",
                        str,
                        "'",
                        line,
                    )?;
                    return Ok(done(entries, repeated, error));
                };
                let name = truncate_spaces(&str[beg..end]);
                if name.is_empty() {
                    let error = registry_exception(
                        throw_site(704, "CRegistryException::eSection"),
                        "Unnamed registry section: `",
                        str,
                        "'",
                        line,
                    )?;
                    return Ok(done(entries, repeated, error));
                }
                if !is_name(name) {
                    let error = registry_exception(
                        throw_site(708, "CRegistryException::eSection"),
                        "Invalid registry section name: `",
                        str,
                        "'",
                        line,
                    )?;
                    return Ok(done(entries, repeated, error));
                }
                section = String::from_utf8_lossy(name).into_owned();
            }
            _ => {
                let Some(equals) = str.iter().position(|&c| c == b'=') else {
                    let error = registry_exception(
                        throw_site(722, "CRegistryException::eEntry"),
                        "Invalid registry entry format: '",
                        str,
                        "'",
                        line,
                    )?;
                    return Ok(done(entries, repeated, error));
                };
                let name = truncate_spaces(&str[..equals]);
                if !is_name(name) {
                    let error = registry_exception(
                        throw_site(728, "CRegistryException::eEntry"),
                        "Invalid registry entry name: '",
                        str,
                        "'",
                        line,
                    )?;
                    return Ok(done(entries, repeated, error));
                }
                let mut value = truncate_spaces(&str[equals + 1..]).to_vec();
                let mut text = str.to_vec();
                while backslashed(&value, value.len()) {
                    let Some(cont) = lines.next_line() else {
                        break;
                    };
                    line += 1;
                    let last = value.len() - 1;
                    value[last] = b'\n';
                    value.extend_from_slice(truncate_spaces(cont));
                    text.push(b'n');
                    text.extend_from_slice(cont);
                }
                let mut beg = 0;
                let mut end = value.len();
                let mut pos = value.iter().position(|&c| c == b'"');
                while let Some(at) = pos.filter(|&at| at < end) {
                    if backslashed(&value, at) {
                    } else if at == beg {
                        beg += 1;
                    } else if at == end - 1 {
                        end -= 1;
                    } else {
                        let error = registry_exception(
                            throw_site(767, "CRegistryException::eValue"),
                            "Single(unescaped) '\"' in the middle of registry value: '",
                            &text,
                            "'",
                            line,
                        )?;
                        return Ok(done(entries, repeated, error));
                    }
                    pos = value[at + 1..]
                        .iter()
                        .position(|&c| c == b'"')
                        .map(|next| at + 1 + next);
                }
                let Ok(value) = parse_escapes(&value[beg..end]) else {
                    let error = registry_exception(
                        throw_site(777, "CRegistryException::eValue"),
                        "Badly placed '\\' in the registry value: '",
                        &text,
                        "'",
                        line,
                    )?;
                    return Ok(done(entries, repeated, error));
                };
                if !section.is_empty() {
                    let name = String::from_utf8_lossy(name).into_owned();
                    let inherits = section.eq_ignore_ascii_case("NCBI")
                        && name.eq_ignore_ascii_case(".Inherits");
                    if !inherits
                        && last_value(&entries, &section, &name).is_some_and(|old| !old.is_empty())
                    {
                        repeated.push(RepeatedEntry {
                            section: section.clone(),
                            name: name.clone(),
                            line,
                        });
                    }
                    entries.push(RegistryEntry {
                        section: section.clone(),
                        name,
                        value,
                    });
                }
            }
        }
    }
    Ok(RegistryRead {
        entries,
        error: None,
        stream_failed: false,
        repeated,
    })
}

#[cfg(test)]
#[path = "ncbi_registry_tests.rs"]
mod ncbi_registry_tests;
