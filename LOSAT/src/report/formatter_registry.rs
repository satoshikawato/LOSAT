//! The registry that NCBI's formatter reads by itself: `./.ncbirc` of the current directory, which
//! the formatting code opens with a registry of its own each time it needs a key, apart from the
//! application's registry (`blastinput::ncbi_environment`, which also finds `~/.ncbirc`,
//! `$NCBI/.ncbirc` and `NCBI_CONFIG_PATH`, and which `NCBI_DONT_USE_NCBIRC` turns off; neither
//! touches this read).
//!
//! The organism report (outfmt 18, `report::tax_report`) reads `[TAX_BROWSER] BLASTFMTUTIL` and
//! `[BLASTFMTUTIL] PROTOCOL` here. NCBI's `-html` forms read the same file for the references'
//! scheme (`CBlastFormatUtil::BlastPrintReference`, `report::reference::write_reference`); the
//! `[LOCAL_ID]` links read another one, through `CAlignFormatUtil::InitConfig`: the file that
//! `FMTCFG` names (default `.ncbirc`) in the current directory, else the one under `$NCBI`
//! (line RG, item T4a).

use std::io::Write;

use crate::blastinput::ncbi_registry::RegistryRead;

/// `./.ncbirc` as NCBI's formatter reads it: `Ok(None)` without the file, `Err(what)` for a file
/// that LOSAT refuses: one whose read throws (a syntax error: NCBI's organism report then fails,
/// NCBI BLAST+ 2.17.0 exit 255 after the first query's `Query=` block, probe `T/t2/probe_rc`) and
/// one that reads other files (`[NCBI] .Inherits`, which LOSAT's registry reads refuse elsewhere
/// too). The caller words the refusal (`what` completes "the registry file ./.ncbirc ..."). The
/// messages of a read that does not throw are `registry_posts`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:171-177
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1615-1625
/// One `Read` of the stream into a file layer (the byte order mark checked once), as the
/// system-wide `.ncbirc` is read; flags 0 include no other `.ncbirc`.
pub(crate) fn read_cwd_ncbirc() -> Result<Option<RegistryRead>, String> {
    use crate::blastinput::ncbi_registry::{read_registry, Layer};
    let Ok(bytes) = std::fs::read(".ncbirc") else {
        return Ok(None);
    };
    let read = read_registry(&bytes, Layer::SystemWide)?;
    if read.error.is_some() {
        return Err("cannot be read as a registry".to_string());
    }
    if read.last_value("NCBI", ".Inherits").is_some() {
        return Err("sets [NCBI] .Inherits".to_string());
    }
    Ok(Some(read))
}

/// What one read of `./.ncbirc` (`read`) posts while `program` (in lower case) runs: a warning of
/// each entry that sets an entry again, then the error of a stream that failed
/// (`IRWRegistry::Read`, as the registry of `-window_masker_taxid` posts them,
/// `ncbi_environment::window_masker_ncbirc`). NCBI BLAST+ 2.17.0 (probe `T/t2/probe_rc`, `rep` and
/// `bom`): two of each per query of outfmt 18, `Warning: (110.8) [blastp] Found multiple [FOO] a
/// settings; using the one from line 3` and `Error: (110.4) [blastp] Error reading the registry
/// after line 1: `.
pub(crate) fn registry_posts(read: &RegistryRead, program: &str) -> Vec<u8> {
    use crate::blastinput::ncbi_registry::{repeated_entry_post, running_error_reading_post};
    let mut posts = String::new();
    for repeated in &read.repeated {
        posts.push_str(&repeated_entry_post(repeated, program));
    }
    if read.stream_failed {
        posts.push_str(&running_error_reading_post(program));
    }
    posts.into_bytes()
}

/// `Get(section, name)` of the read `registry` (`None`: no file) for a caller that takes an empty
/// value as none (`Get(...).empty()`, as `[TAX_BROWSER] BLASTFMTUTIL` is read): the last value of
/// the entry, `None` where it is not set or empty. (A `HasEntry` of NCBI's `CNcbiRegistry` counts
/// an empty value as set: `protocol`.)
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:983-991
pub(crate) fn value<'a>(
    registry: Option<&'a RegistryRead>,
    section: &str,
    name: &str,
) -> Option<&'a [u8]> {
    registry
        .and_then(|read| read.last_value(section, name))
        .filter(|value| !value.is_empty())
}

/// The scheme of NCBI's URLs: the last value of `[BLASTFMTUTIL] PROTOCOL`, also an empty one, else
/// `https:`. `CNcbiRegistry` is a compound registry, whose `HasEntry` asks its layers with
/// `fCountCleared`, so an entry set to an empty value is found: NCBI BLAST+ 2.17.0 with
/// `PROTOCOL=` writes `//www.ncbi.nlm.nih.gov/...` in the organism report's `<@taxBrowserURL@>`
/// and in the `-html` references (line RG T4 probes `rc_p_empty`, review of T4a0 `pe`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3643-3655
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1236-1245, 1300-1303
pub(crate) fn protocol(registry: Option<&RegistryRead>) -> &[u8] {
    registry
        .and_then(|read| read.last_value("BLASTFMTUTIL", "PROTOCOL"))
        .unwrap_or(super::reference::DEFAULT_PROTOCOL)
}

/// The organism report's taxonomy browser URL before its `<@protocol@>` is replaced
/// (`m_TaxBrowserURL`): `[TAX_BROWSER] BLASTFMTUTIL` (`value`), else NCBI's `TAX_BROWSER_URL`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:175-177
pub(crate) fn tax_browser_template(registry: Option<&RegistryRead>) -> &[u8] {
    const TAX_BROWSER_URL: &[u8] =
        b"<@protocol@>//www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi";
    value(registry, "TAX_BROWSER", "BLASTFMTUTIL").unwrap_or(TAX_BROWSER_URL)
}

/// `CAlignFormatUtil::MapTemplate` on bytes: every `<@name@>` of `input` replaced by `value`,
/// left to right; the text put in is not searched again by this call.
///
/// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:3593-3599
/// ```c++
/// string CAlignFormatUtil::MapTemplate(string inpString,string tmplParamName,string templParamVal)
/// {
///     string outString;
///     string tmplParam = "<@" + tmplParamName + "@>";
///     NStr::Replace(inpString,tmplParam,templParamVal,outString);
///     return outString;
/// }
/// ```
pub(crate) fn map_template(input: &[u8], name: &str, value: &[u8]) -> Vec<u8> {
    let token = [b"<@", name.as_bytes(), b"@>"].concat();
    let mut out = Vec::with_capacity(input.len());
    let mut rest = input;
    while let Some(at) = rest
        .windows(token.len())
        .position(|w| w == token.as_slice())
    {
        out.extend_from_slice(&rest[..at]);
        out.extend_from_slice(value);
        rest = &rest[at + token.len()..];
    }
    out.extend_from_slice(rest);
    out
}

/// `NStr::URLEncode(text)` (`eUrlEnc_SkipMarkChars`, table `s_Encode`): letters, digits and
/// `!$'()*,-._` stay, a space becomes `+`, every other byte `%XX` (upper-case hex).
///
/// NCBI reference (2.17.0): c++/src/corelib/ncbistr.cpp:5614-5619, 6066-6067 (`s_Encode`, rows of 0x20-0x3F)
/// ```c++
///     "+",   "!",   "%22", "%23", "$",   "%25", "%26", "'",
///     "(",   ")",   "*",   "%2B", ",",   "-",   ".",   "%2F",
///     "0",   "1",   "2",   "3",   "4",   "5",   "6",   "7",
///     "8",   "9",   "%3A", "%3B", "%3C", "%3D", "%3E", "%3F",
/// ...
///     case eUrlEnc_SkipMarkChars:
///         encode_table = s_Encode;
/// ```
/// The other rows keep `A`-`Z`, `a`-`z` and `_` only (`%40`, `%5B`-`%5E`, `%60`, `%7B`-`%7F`).
pub(crate) fn url_encode(text: &[u8]) -> Vec<u8> {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut out = Vec::with_capacity(text.len());
    for &byte in text {
        match byte {
            b' ' => out.push(b'+'),
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' => out.push(byte),
            b'!' | b'$' | b'\'' | b'(' | b')' | b'*' | b',' | b'-' | b'.' | b'_' => out.push(byte),
            _ => out.extend_from_slice(&[
                b'%',
                HEX[usize::from(byte >> 4)],
                HEX[usize::from(byte & 15)],
            ]),
        }
    }
    out
}

/// The file of `CAlignFormatUtil::InitConfig`: `FMTCFG` (unset or empty: `.ncbirc`) where it
/// is a file (`CFile::Exists`, links followed), else `$NCBI` (not empty) with `/` added where it
/// does not end in one, then that name, where it is a file; `None`: no file.
///
/// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:3383-3418
/// ```c++
/// void CAlignFormatUtil::InitConfig()
/// {
///     string l_cfg_file_name;
///     bool   l_dbg = CAlignFormatUtil::m_geturl_debug_flag;
///     if( getenv("GETURL_DEBUG") ) CAlignFormatUtil::m_geturl_debug_flag = l_dbg = true;
///     if( !m_Reg ) {
///         ...
///         if( NULL !=  getenv("NCBI")   ) l_ncbi_env = getenv("NCBI");
///         if( NULL !=  getenv("FMTCFG") ) l_fmtcfg_env = getenv("FMTCFG");
///         // config file name: value of FMTCFG or  default ( .ncbirc )
///         if( l_fmtcfg_env.empty()  )
///             l_cfg_file_name = ".ncbirc";
///         else
///             l_cfg_file_name = l_fmtcfg_env;
///         // checkinf existance of configuration file
///         CFile  l_fchecker( l_cfg_file_name );
///         cfgExists = l_fchecker.Exists();
///         if( (!cfgExists) && (!l_ncbi_env.empty()) ) {
///             if( l_ncbi_env.rfind("/") != (l_ncbi_env.length() -1 ))
///             l_ncbi_env.append("/");
///             l_cfg_file_name = l_ncbi_env + l_cfg_file_name;
///             CFile  l_fchecker2( l_cfg_file_name );
///             cfgExists = l_fchecker2.Exists();
///         }
///         if(cfgExists) {
///             CNcbiIfstream l_ConfigFile(l_cfg_file_name.c_str() );
///             m_Reg.reset(new CNcbiRegistry(l_ConfigFile));
///             if( l_dbg ) fprintf(stderr,"REGISTRY: %s\n",l_cfg_file_name.c_str());
///         }
///     }
/// ```
fn init_config_path(
    fmtcfg: Option<&std::ffi::OsStr>,
    ncbi: Option<&std::ffi::OsStr>,
    is_file: &dyn Fn(&std::path::Path) -> bool,
) -> Option<std::ffi::OsString> {
    let name = fmtcfg
        .filter(|name| !name.is_empty())
        .unwrap_or(std::ffi::OsStr::new(".ncbirc"));
    if is_file(name.as_ref()) {
        return Some(name.to_owned());
    }
    let ncbi = ncbi.filter(|ncbi| !ncbi.is_empty())?;
    let mut path = ncbi.to_owned();
    if !ncbi.as_encoded_bytes().ends_with(b"/") {
        path.push("/");
    }
    path.push(name);
    is_file(path.as_ref()).then_some(path)
}

/// What `InitConfig` found (`CAlignFormatUtil::m_Reg`, static in NCBI: one per run).
#[derive(Debug)]
enum ConfigFile {
    /// No file: `m_Reg` stays null and NCBI looks again at each call (a run's files and
    /// environment do not change, so each look finds none again).
    Absent,
    /// The file `path` and its read (`m_Reg`), with the posts of the read (`registry_posts`).
    Read {
        path: std::ffi::OsString,
        registry: RegistryRead,
        posts: Vec<u8>,
    },
    /// A file whose read throws (`CRegistryException`, `ncbireg.cpp` line 722): `m_Reg` stays
    /// null, so each call reads it again, posts again and throws again.
    Throws {
        posts: Vec<u8>,
        exception: crate::ncbi_exception::NcbiException,
    },
}

/// The registry of NCBI's id links with `-html` (`CAlignFormatUtil::InitConfig`, `m_Reg`,
/// `m_Protocol`, `m_geturl_debug_flag`): one per run, made only with `-html` (`InitConfig::read`,
/// which reads the file at once so that LOSAT's refusals come before the first byte), and used
/// where NCBI first calls `InitConfig`: at the first id URL of the report (`map_protocol`,
/// `local_id_url`), never by the organism report. The effects of a read, the posts, the
/// `GETURL_DEBUG` line and the exception of a malformed file, happen at that call, as in NCBI
/// (line RG T4 probes `env_geturl`, `cfg.*`: a malformed file fails only with hits, rc 255).
#[derive(Debug)]
pub(crate) struct InitConfig {
    file: ConfigFile,
    /// `GETURL_DEBUG` is set (to any value, the empty one too).
    debug: bool,
    /// Whether the first call has written the posts and the `GETURL_DEBUG` line of a `Read`.
    reported: std::sync::atomic::AtomicBool,
}

impl InitConfig {
    /// The registry of the run of `program` (in lower case, as NCBI's posts name it), from the
    /// environment (`FMTCFG`, `NCBI`, `GETURL_DEBUG`) and the current directory. LOSAT refuses a
    /// file that cannot be read and one with `[NCBI] .Inherits` (`CNcbiRegistry` reads the files
    /// it names, which LOSAT's registry reads refuse elsewhere too); `what` completes "the
    /// registry file <path> ...".
    pub(crate) fn read(program: &str) -> Result<Self, String> {
        let path = init_config_path(
            std::env::var_os("FMTCFG").as_deref(),
            std::env::var_os("NCBI").as_deref(),
            &|path| path.is_file(),
        );
        let file = match path {
            None => None,
            Some(path) => {
                let bytes = std::fs::read(&path).map_err(|error| {
                    format!(
                        "the registry file {} ({error}), which NCBI BLAST+'s -html links read, cannot be read; this is not supported by LOSAT",
                        std::path::Path::new(&path).display()
                    )
                })?;
                Some((path, bytes))
            }
        };
        Self::of(program, file, std::env::var_os("GETURL_DEBUG").is_some())
    }

    /// `read` of the file `file` (its path as NCBI names it and its bytes; `None`: no file),
    /// `debug`: `GETURL_DEBUG` is set.
    pub(super) fn of(
        program: &str,
        file: Option<(std::ffi::OsString, Vec<u8>)>,
        debug: bool,
    ) -> Result<Self, String> {
        use crate::blastinput::ncbi_registry::{read_registry, Layer};
        let refused = |path: &std::ffi::OsStr, what: &str| {
            format!(
                "the registry file {} {what}, which NCBI BLAST+'s -html links read (FMTCFG, NCBI); this is not supported by LOSAT",
                std::path::Path::new(path).display()
            )
        };
        let file = match file {
            None => ConfigFile::Absent,
            Some((path, bytes)) => {
                let mut registry = read_registry(&bytes, Layer::SystemWide)
                    .map_err(|reason| refused(&path, &reason))?;
                if registry.last_value("NCBI", ".Inherits").is_some() {
                    return Err(refused(&path, "sets [NCBI] .Inherits"));
                }
                let posts = registry_posts(&registry, program);
                match registry.error.take() {
                    Some(exception) => ConfigFile::Throws { posts, exception },
                    None => ConfigFile::Read {
                        path,
                        registry,
                        posts,
                    },
                }
            }
        };
        Ok(InitConfig {
            file,
            debug,
            reported: Default::default(),
        })
    }

    /// `InitConfig()` where `m_Reg` is null, then `m_Reg`: the registry (`None`: no file). The
    /// first call writes the posts of the read and, with `GETURL_DEBUG`, `REGISTRY: <path>` to
    /// `diagnostics`; a file whose read throws posts and throws at every call (the exception as
    /// the report's I/O error, `NcbiException::into_io_error`: NCBI's application then ends with
    /// it, exit status 255).
    fn registry(&self, diagnostics: &mut dyn Write) -> std::io::Result<Option<&RegistryRead>> {
        match &self.file {
            ConfigFile::Absent => Ok(None),
            ConfigFile::Read {
                path,
                registry,
                posts,
            } => {
                if !self
                    .reported
                    .swap(true, std::sync::atomic::Ordering::Relaxed)
                {
                    diagnostics.write_all(posts)?;
                    if self.debug {
                        diagnostics.write_all(b"REGISTRY: ")?;
                        diagnostics.write_all(path.as_encoded_bytes())?;
                        diagnostics.write_all(b"\n")?;
                    }
                    diagnostics.flush()?;
                }
                Ok(Some(registry))
            }
            ConfigFile::Throws { posts, exception } => {
                diagnostics.write_all(posts)?;
                diagnostics.flush()?;
                Err(exception.clone().into_io_error())
            }
        }
    }

    /// `CAlignFormatUtil::MapProtocol(url)`: `<@protocol@>` replaced by `[BLASTFMTUTIL] PROTOCOL`
    /// of the registry (also an empty value: `HasEntry`), else `https:`. NCBI applies it to the URL
    /// of every id (`GetIDUrlGen`), so a malformed file fails at the first id of any kind (probes
    /// `cfg.*.{L,A,D}.h`).
    ///
    /// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:3661-3671
    /// ```c++
    /// string CAlignFormatUtil::MapProtocol(string url_link)
    /// {
    ///     if(m_Protocol.empty()){
    ///         if(!m_Reg) {
    ///             InitConfig();
    ///         }
    ///         m_Protocol = (m_Reg && m_Reg->HasEntry("BLASTFMTUTIL","PROTOCOL")) ? m_Protocol = m_Reg->Get("BLASTFMTUTIL","PROTOCOL") : "https:";
    ///     }
    ///     url_link = CAlignFormatUtil::MapTemplate(url_link,"protocol",m_Protocol);
    ///     return url_link;
    /// }
    /// ```
    pub(crate) fn map_protocol(
        &self,
        url: &[u8],
        diagnostics: &mut dyn Write,
    ) -> std::io::Result<Vec<u8>> {
        let registry = self.registry(diagnostics)?;
        Ok(map_template(url, "protocol", protocol(registry)))
    }

    /// The URL of a local id (`GetIDUrlGen`, `CSeq_id::e_Local`), after `MapProtocol`: the
    /// registry's `[LOCAL_ID] TOOL_URL` (empty without it or without a file: the caller then has no
    /// URL) with its tags filled in NCBI's order: `seq_id` (`NStr::URLEncode` of the id's
    /// `eContent` label `label`, `Subject_1` for `-subject` records too), `db_name` (that of the
    /// `-db` text `database`, empty for `-subject`), `taxid`, `cssInf` (empty), `title` (`label`),
    /// `target` (empty), then `protocol`. `GetURLFromRegistry("LOCAL_ID")` comes first and reads
    /// `[BLASTFMTUTIL] INCLUDE_BASE_DIR`, `LOCAL_ID_HOST_PORT` and `LOCAL_ID_FORMAT` (and a file
    /// they name), whose result NCBI never uses: LOSAT reads none of them. `TOOL_URL_ALIGN` is read
    /// only with `addCssInfo`, which no NCBI application sets.
    ///
    /// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:3816-3838
    /// ```c++
    ///         } else if (wid->Which() == CSeq_id::e_Local){
    ///
    ///             string url_holder = CAlignFormatUtil::GetURLFromRegistry("LOCAL_ID");
    ///
    ///             string user_url;
    ///             if (m_Reg) {
    ///                 user_url = (seqUrlInfo->addCssInfo) ? m_Reg->Get("LOCAL_ID","TOOL_URL_ALIGN") : m_Reg->Get("LOCAL_ID","TOOL_URL");
    ///             }
    ///             string id_string;
    ///             wid->GetLabel(&id_string, CSeq_id::eContent);
    ///             url_link = CAlignFormatUtil::MapTemplate(user_url,"seq_id", NStr::URLEncode(id_string));
    ///             url_link = CAlignFormatUtil::MapTemplate(url_link,"db_name", NStr::URLEncode(seqUrlInfo->database));
    ///             url_link = CAlignFormatUtil::MapTemplate(url_link,"taxid", TAX_ID_TO(int, seqUrlInfo->taxid));
    ///             temp_class_info = (!seqUrlInfo->defline.empty())? CAlignFormatUtil::MapTemplate(temp_class_info,"defline",seqUrlInfo->defline):temp_class_info;
    ///             url_link = CAlignFormatUtil::MapTemplate(url_link,"cssInf",(seqUrlInfo->addCssInfo) ? temp_class_info.c_str() : "");
    ///             url_link = CAlignFormatUtil::MapTemplate(url_link,"title", id_string);
    ///             url_link = CAlignFormatUtil::MapTemplate(url_link,"target",seqUrlInfo->new_win ? "TARGET=\"EntrezView\"" : "");
    ///         }
    ///     }
    ///     url_link = CAlignFormatUtil::MapProtocol(url_link);
    /// ```
    pub(crate) fn local_id_url(
        &self,
        label: &[u8],
        database: &[u8],
        taxid: i32,
        diagnostics: &mut dyn Write,
    ) -> std::io::Result<Vec<u8>> {
        let registry = self.registry(diagnostics)?;
        let tool_url = registry
            .and_then(|read| read.last_value("LOCAL_ID", "TOOL_URL"))
            .unwrap_or_default();
        let url = map_template(tool_url, "seq_id", &url_encode(label));
        let url = map_template(&url, "db_name", &url_encode(database));
        let url = map_template(&url, "taxid", taxid.to_string().as_bytes());
        let url = map_template(&url, "cssInf", b"");
        let url = map_template(&url, "title", label);
        let url = map_template(&url, "target", b"");
        self.map_protocol(&url, diagnostics)
    }
}

#[cfg(test)]
#[path = "formatter_registry_tests.rs"]
mod tests;
