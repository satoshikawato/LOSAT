//! The registry that NCBI's formatter reads by itself: `./.ncbirc` of the current directory, which
//! the formatting code opens with a registry of its own each time it needs a key, apart from the
//! application's registry (`blastinput::ncbi_environment`, which also finds `~/.ncbirc`,
//! `$NCBI/.ncbirc` and `NCBI_CONFIG_PATH`, and which `NCBI_DONT_USE_NCBIRC` turns off; neither
//! touches this read).
//!
//! The organism report (outfmt 18, `report::tax_report`) reads `[TAX_BROWSER] BLASTFMTUTIL` and
//! `[BLASTFMTUTIL] PROTOCOL` here. NCBI's `-html` pairwise report reads the same file for the
//! references' scheme (`CBlastFormatUtil::BlastPrintReference`, `report::reference::write_reference`)
//! and again at the start of each query's description table and of its alignments
//! (`FormatterNcbirc`); the `[LOCAL_ID]` links read another one, through
//! `CAlignFormatUtil::InitConfig`: the file that `FMTCFG` names (default `.ncbirc`) in the current
//! directory, else the one under `$NCBI` (line RG, item T4a). The reads of `-html` have no effect
//! here: what each writes to the run's diagnostics is data (`FormatterNcbirc::posts`,
//! `InitConfig::call_posts`), which the report's writer writes where NCBI's read comes in the report
//! (`report::html_ids::HtmlPost`).

use crate::blastinput::ncbi_registry::RegistryRead;

/// The name of the file in the current directory that NCBI's formatter opens by itself.
const CWD_NCBIRC: &str = ".ncbirc";

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1615-1625
// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:605-614 (a stream that has failed gives
// no read)
/// A registry file of NCBI's formatter read as NCBI reads it (`./.ncbirc` of the references, the
/// description tables, the alignments and the organism report, and `InitConfig`'s file): a
/// `CNcbiIfstream` of the file `path`, which `CNcbiRegistry` reads with one `Read` into a file
/// layer (the byte order mark checked once; flags 0 include no other `.ncbirc`). `Ok(None)` where
/// the stream does not open (no file, a broken link, a file the process may not read): the read
/// returns at the failed stream, an empty registry and no post. A file that opens and whose first
/// read fails (a directory) is a stream that fails before its first line
/// (`ncbi_registry::read_failed_stream`), whose read posts `Error reading the registry after line
/// 1` (`registry_posts`). NCBI BLAST+ 2.17.0 (line RG T4a2 fix round 2, `T/t4/fix2` probes `c1.*`,
/// `./.ncbirc` a directory, mode 000, a broken link, an empty file) agrees at each of these reads.
/// `Err(what)` (completing "the registry file `PATH` ...") for what LOSAT refuses: a read that
/// fails after some bytes (the lines that NCBI reads before it follow its stream's buffer), and
/// what `read_registry` refuses.
fn read_stream(path: &std::path::Path) -> Result<Option<RegistryRead>, String> {
    use crate::blastinput::ncbi_registry::{read_failed_stream, read_registry, Layer};
    let Ok(mut file) = std::fs::File::open(path) else {
        return Ok(None);
    };
    let mut bytes = Vec::new();
    match std::io::Read::read_to_end(&mut file, &mut bytes) {
        Ok(_) => read_registry(&bytes, Layer::SystemWide).map(Some),
        Err(_) if bytes.is_empty() => Ok(Some(read_failed_stream(Layer::SystemWide))),
        Err(error) => Err(format!(
            "cannot be read after its first {} bytes ({error})",
            bytes.len()
        )),
    }
}

/// `./.ncbirc` as NCBI's organism report reads it (`read_stream`): `Ok(None)` where it does not
/// open, `Err(what)` for a file that LOSAT refuses: one whose read throws (a syntax error: NCBI's
/// organism report then fails, NCBI BLAST+ 2.17.0 exit 255 after the first query's `Query=`
/// block, probe `T/t2/probe_rc`) and one that reads other files (`[NCBI] .Inherits`, which LOSAT's
/// registry reads refuse elsewhere too). The caller words the refusal (`what` completes "the
/// registry file ./.ncbirc ..."). The messages of a read that does not throw are
/// `registry_posts`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:171-177
pub(crate) fn read_cwd_ncbirc() -> Result<Option<RegistryRead>, String> {
    let Some(read) = read_stream(CWD_NCBIRC.as_ref())? else {
        return Ok(None);
    };
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

/// `Get` of `section` and `name` in the read `registry` (`None`: no file) for a caller that takes
/// an empty value as none (it tests whether `Get` gives an empty string, as
/// `[TAX_BROWSER] BLASTFMTUTIL` is read): the last value of the entry, `None` where it is not set
/// or empty. (A `HasEntry` of NCBI's `CNcbiRegistry` counts an empty value as set: `protocol`.)
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
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3593-3599
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

/// `NStr::URLEncode` of `text` (`eUrlEnc_SkipMarkChars`, table `s_Encode`): letters, digits and
/// `!$'()*,-._` stay, a space becomes `+`, every other byte `%XX` (upper-case hex).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5614-5619, 6066-6067 (`s_Encode`, rows of 0x20-0x3F)
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
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3383-3418
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

/// What a read of a registry file of NCBI's formatter found: `InitConfig`'s file
/// (`CAlignFormatUtil::m_Reg`, static in NCBI: one per run) or the `./.ncbirc` of the `-html`
/// pairwise report (`FormatterNcbirc`).
#[derive(Debug, Default)]
enum ConfigFile {
    /// No file: for `InitConfig` none is found, `m_Reg` stays null and NCBI looks again at each
    /// call (a run's files and environment do not change, so each look finds none again); for the
    /// report's `./.ncbirc`, none that opens (each read makes an empty registry and posts nothing).
    #[default]
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

impl ConfigFile {
    /// What each read of the file writes to the run's diagnostics (`registry_posts`; nothing
    /// without a file).
    fn posts(&self) -> &[u8] {
        match self {
            ConfigFile::Absent => &[],
            ConfigFile::Read { posts, .. } | ConfigFile::Throws { posts, .. } => posts,
        }
    }

    /// The registry of a read (`None`: no file), or for a file whose read throws its exception as
    /// the report's I/O error (`NcbiException::into_io_error`: NCBI's application ends with it,
    /// exit status 255).
    fn registry(&self) -> std::io::Result<Option<&RegistryRead>> {
        match self {
            ConfigFile::Absent => Ok(None),
            ConfigFile::Read { registry, .. } => Ok(Some(registry)),
            ConfigFile::Throws { exception, .. } => Err(exception.clone().into_io_error()),
        }
    }

    /// The file `path` of the read `registry` (`read_stream`) while `program` (in lower case)
    /// runs. `Err(what)` for a file that LOSAT refuses (completing `the registry file <path>
    /// ...`): one with `[NCBI] .Inherits` (NCBI then reads the files it names, which LOSAT's
    /// registry reads refuse elsewhere too).
    fn of(
        program: &str,
        path: std::ffi::OsString,
        mut registry: RegistryRead,
    ) -> Result<Self, String> {
        if registry.last_value("NCBI", ".Inherits").is_some() {
            return Err("sets [NCBI] .Inherits".to_string());
        }
        let posts = registry_posts(&registry, program);
        Ok(match registry.error.take() {
            Some(exception) => ConfigFile::Throws { posts, exception },
            None => ConfigFile::Read {
                path,
                registry,
                posts,
            },
        })
    }
}

/// `./.ncbirc` as NCBI's `-html` pairwise report reads it: `CBlastFormatUtil::BlastPrintReference`
/// makes a registry of the file for each reference of the prolog and takes `[BLASTFMTUTIL]
/// PROTOCOL` from it (`protocol`), and each query with hits reads it again where its description
/// table starts (`CShowBlastDefline::x_InitDeflineTable`, also with `-num_descriptions 0`) and where
/// its alignments start (`CDisplaySeqalign::x_InitAlignParams`, where an alignment is shown). Read
/// once at the program's setup, with the run's other `-html` data (`html_ids::HtmlLinks`,
/// `blastinput::app::html_links`), so that LOSAT's refusals of the file come before the first byte;
/// each of NCBI's reads then has the effects of its own: the posts of the read (`posts`, which the
/// report's writer writes at the read's place, `html_ids::HtmlPost::Ncbirc`), and for a file whose
/// read throws the exception after them (`read`), which ends the prolog after the version line
/// (line RG T4 probes `bp.di.f0.rc_bad.h`, `rc_bad_dontuse.h`: exit status 255). The default is the
/// run without the file. NCBI BLAST+ 2.17.0 with a file that sets an entry twice (line RG T4a2 fix
/// round 1, `T/t4/fix1` probe `bp2`): two posts for blastp's references, then per query with hits
/// one at its table and one at its alignments.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:96-105
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1068-1070
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1882-1885
#[derive(Debug, Default)]
pub(crate) struct FormatterNcbirc(ConfigFile);

impl FormatterNcbirc {
    /// The read of `./.ncbirc` of the current directory for a run of `program` (in lower case, as
    /// NCBI's posts name it, `read_stream`); a file that does not open is the default (an empty
    /// registry without posts). `Err(what)` for a file that LOSAT refuses (`read_stream`,
    /// `ConfigFile::of`): the message of the refusal, which the caller ends with the program's
    /// name.
    pub(crate) fn read(program: &str) -> Result<Self, String> {
        Self::of_read(program, read_stream(CWD_NCBIRC.as_ref()))
    }

    /// `read` of the file's read `read` (`read_stream`).
    fn of_read(program: &str, read: Result<Option<RegistryRead>, String>) -> Result<Self, String> {
        let file = read.and_then(|read| match read {
            None => Ok(ConfigFile::Absent),
            Some(registry) => ConfigFile::of(program, CWD_NCBIRC.into(), registry),
        });
        file.map(FormatterNcbirc).map_err(|what| {
            format!(
                "the registry file ./.ncbirc {what}, which NCBI BLAST+'s -html reports read; this is not supported by LOSAT"
            )
        })
    }

    /// `read` of the file's bytes `bytes` (`None`: a file that does not open).
    #[cfg(test)]
    pub(crate) fn of(program: &str, bytes: Option<Vec<u8>>) -> Result<Self, String> {
        use crate::blastinput::ncbi_registry::{read_registry, Layer};
        let read = bytes
            .map(|bytes| read_registry(&bytes, Layer::SystemWide))
            .transpose();
        Self::of_read(program, read)
    }

    /// What each of NCBI's reads of the file writes to the run's diagnostics (standard error)
    /// before it gives its registry or throws (`registry`).
    pub(crate) fn posts(&self) -> &[u8] {
        self.0.posts()
    }

    /// The registry that each read makes (`None` without the file), or the exception of a file
    /// whose read throws, after its posts (`posts`).
    pub(crate) fn registry(&self) -> std::io::Result<Option<&RegistryRead>> {
        self.0.registry()
    }

    /// The scheme of one reference's PubMed link (`httpProt`) after the reference's read
    /// (`registry`): the file's `[BLASTFMTUTIL] PROTOCOL` (`protocol`: an empty value counts as
    /// set), `https:` without the file.
    pub(crate) fn protocol(&self) -> std::io::Result<&[u8]> {
        Ok(protocol(self.registry()?))
    }
}

/// The registry of NCBI's id links with `-html` (`CAlignFormatUtil::InitConfig`, `m_Reg`,
/// `m_Protocol`, `m_geturl_debug_flag`): one per run, made only with `-html` (`InitConfig::read`,
/// which reads the file at once so that LOSAT's refusals come before the first byte), and used
/// where NCBI calls `InitConfig`: at each id URL of the report (`map_protocol`, `local_id_url`),
/// never by the organism report. The effects of a call, the posts and the `GETURL_DEBUG` line of
/// the run's first (`call_posts`) and the exception of a malformed file at every call, happen at
/// the id URL's place in the report, as in NCBI (line RG T4 probes `env_geturl`, `cfg.*`: a
/// malformed file fails only with hits, rc 255): the exception is the call's result, the text is
/// written by the report's writer where the call comes in the report (`html_ids::HtmlPost`).
#[derive(Debug)]
pub(crate) struct InitConfig {
    file: ConfigFile,
    /// `GETURL_DEBUG` is set (to any value, the empty one too).
    debug: bool,
    /// By the report's place among the search's formats, whether it has written the posts and the
    /// `GETURL_DEBUG` line of a `Read`, those of its first call (`call_posts`).
    reported: std::sync::Mutex<Vec<bool>>,
}

impl InitConfig {
    /// The registry of the run of `program` (in lower case, as NCBI's posts name it), from the
    /// environment (`FMTCFG`, `NCBI`, `GETURL_DEBUG`) and the current directory: the file that
    /// `init_config_path` finds, read as NCBI's formatter reads its files (`read_stream`). A file
    /// that does not open is read as an empty registry, which NCBI keeps (`m_Reg`) and names in the
    /// `GETURL_DEBUG` line (NCBI BLAST+ 2.17.0, line RG T4a2 fix round 2, probes
    /// `T/t4/fix2/out/c1.noperm.*`: mode 000, rc 0, `REGISTRY: .ncbirc`, the default links). LOSAT
    /// refuses what `read_stream` refuses and a file with `[NCBI] .Inherits` (`CNcbiRegistry`
    /// reads the files it names, which LOSAT's registry reads refuse elsewhere too); `what`
    /// completes `the registry file <path> ...`.
    pub(crate) fn read(program: &str) -> Result<Self, String> {
        let path = init_config_path(
            std::env::var_os("FMTCFG").as_deref(),
            std::env::var_os("NCBI").as_deref(),
            &|path| path.is_file(),
        );
        Self::read_file(program, path, std::env::var_os("GETURL_DEBUG").is_some())
    }

    /// `read` of the file `path` that `init_config_path` found (`None`: none), `debug`:
    /// `GETURL_DEBUG` is set.
    fn read_file(
        program: &str,
        path: Option<std::ffi::OsString>,
        debug: bool,
    ) -> Result<Self, String> {
        let file = path.map(|path| {
            let read = read_stream(path.as_ref()).map(Option::unwrap_or_default);
            (path, read)
        });
        Self::of_read(program, file, debug)
    }

    /// `read` of the file `file` (its path as NCBI names it and its read, `read_stream`; `None`:
    /// no file), `debug`: `GETURL_DEBUG` is set.
    fn of_read(
        program: &str,
        file: Option<(std::ffi::OsString, Result<RegistryRead, String>)>,
        debug: bool,
    ) -> Result<Self, String> {
        let file = match file {
            None => ConfigFile::Absent,
            Some((path, read)) => {
                let shown = std::path::Path::new(&path).display().to_string();
                read.and_then(|registry| ConfigFile::of(program, path, registry))
                    .map_err(|what| {
                        format!(
                            "the registry file {shown} {what}, which NCBI BLAST+'s -html links read (FMTCFG, NCBI); this is not supported by LOSAT"
                        )
                    })?
            }
        };
        Ok(InitConfig {
            file,
            debug,
            reported: Default::default(),
        })
    }

    /// `read` of the file `file` (its path as NCBI names it and its bytes; `None`: no file).
    #[cfg(test)]
    pub(super) fn of(
        program: &str,
        file: Option<(std::ffi::OsString, Vec<u8>)>,
        debug: bool,
    ) -> Result<Self, String> {
        use crate::blastinput::ncbi_registry::{read_registry, Layer};
        let file = file.map(|(path, bytes)| (path, read_registry(&bytes, Layer::SystemWide)));
        Self::of_read(program, file, debug)
    }

    /// `InitConfig()` where `m_Reg` is null, then `m_Reg`: the registry (`None`: no file), or for a
    /// file whose read throws its exception at every call (as the report's I/O error,
    /// `NcbiException::into_io_error`: NCBI's application then ends with it, exit status 255).
    fn registry(&self) -> std::io::Result<Option<&RegistryRead>> {
        self.file.registry()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3410-3412
    /// What one call of `InitConfig` (an id URL, `map_protocol`) by the report of the search's
    /// `report`th format writes, where the report reaches it: the posts of the read, to the run's
    /// diagnostics, and with `GETURL_DEBUG` the line `REGISTRY: <path>`, which NCBI writes to the
    /// process's standard error itself, not as a post (`QueryWarnings::html_stderr_line`: with
    /// tblastx's `-logfile` the posts go to the file and the line stays on the terminal). The first
    /// call of a report that reads a file writes both, its later calls nothing (`m_Reg` is set); a
    /// file whose read throws posts at every call (before its exception, `registry`). NCBI's
    /// `m_Reg` is a static of the process, and each report writes what a run of its format alone
    /// writes, so each report of a library call of several formats has its first call (line RG
    /// T4a2 fix round 2, review c C3 = review d D5: formats 0 and 0); the command line has one. The
    /// report's writer asks it in the report's order, so the first call is the report's first id
    /// URL also where the parts of the report are formatted on threads
    /// (`x_pairwise_par::write_parallel`).
    pub(crate) fn call_posts(&self, report: usize) -> (Vec<u8>, Vec<u8>) {
        match &self.file {
            ConfigFile::Absent => Default::default(),
            ConfigFile::Read { path, posts, .. } => {
                let mut reported = self
                    .reported
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                if reported.len() <= report {
                    reported.resize(report + 1, false);
                }
                if std::mem::replace(&mut reported[report], true) {
                    return Default::default();
                }
                let mut line = Vec::new();
                if self.debug {
                    line.extend_from_slice(b"REGISTRY: ");
                    line.extend_from_slice(path.as_encoded_bytes());
                    line.push(b'\n');
                }
                (posts.clone(), line)
            }
            ConfigFile::Throws { posts, .. } => (posts.clone(), Vec::new()),
        }
    }

    /// `CAlignFormatUtil::MapProtocol` of `url`: `<@protocol@>` replaced by `[BLASTFMTUTIL]
    /// PROTOCOL` of the registry (also an empty value: `HasEntry`), else `https:`. NCBI applies it
    /// to the URL of every id (`GetIDUrlGen`), so a malformed file fails at the first id of any
    /// kind (probes `cfg.*.{L,A,D}.h`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3661-3671
    pub(crate) fn map_protocol(&self, url: &[u8]) -> std::io::Result<Vec<u8>> {
        let registry = self.registry()?;
        Ok(map_template(url, "protocol", protocol(registry)))
    }

    /// The URL of a local id (`GetIDUrlGen`, `CSeq_id::e_Local`), after `MapProtocol`: the
    /// registry's `[LOCAL_ID] TOOL_URL` (empty without it or without a file: the caller then has no
    /// URL) with its tags filled in NCBI's order: `seq_id` (`NStr::URLEncode` of the id's
    /// `eContent` label `label`, `Subject_1` for `-subject` records too), `db_name` (that of the
    /// `-db` text `database`, empty for `-subject`), `taxid`, `cssInf` (empty), `title` (`label`),
    /// `target` (empty), then `protocol`. `GetURLFromRegistry` of `LOCAL_ID` comes first and reads
    /// `[BLASTFMTUTIL] INCLUDE_BASE_DIR`, `LOCAL_ID_HOST_PORT` and `LOCAL_ID_FORMAT` (and a file
    /// they name), whose result NCBI never uses: LOSAT reads none of them. `TOOL_URL_ALIGN` is read
    /// only with `addCssInfo`, which no NCBI application sets.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3816-3838
    pub(crate) fn local_id_url(
        &self,
        label: &[u8],
        database: &[u8],
        taxid: i32,
    ) -> std::io::Result<Vec<u8>> {
        let registry = self.registry()?;
        let tool_url = registry
            .and_then(|read| read.last_value("LOCAL_ID", "TOOL_URL"))
            .unwrap_or_default();
        let url = map_template(tool_url, "seq_id", &url_encode(label));
        let url = map_template(&url, "db_name", &url_encode(database));
        let url = map_template(&url, "taxid", taxid.to_string().as_bytes());
        let url = map_template(&url, "cssInf", b"");
        let url = map_template(&url, "title", label);
        let url = map_template(&url, "target", b"");
        self.map_protocol(&url)
    }
}

#[cfg(test)]
#[path = "formatter_registry_tests.rs"]
mod tests;
