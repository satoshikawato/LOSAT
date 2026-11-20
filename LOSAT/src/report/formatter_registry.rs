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

use crate::blastinput::ncbi_registry::RegistryRead;

/// `./.ncbirc` as NCBI's formatter reads it: `Ok(None)` without the file, `Err(what)` for a file
/// that LOSAT refuses: one whose read throws (a syntax error: NCBI's organism report then fails,
/// NCBI BLAST+ 2.17.0 exit 255 after the first query's `Query=` block, probe `T/t2/probe_rc`) and
/// one that reads other files (`[NCBI] .Inherits`, which LOSAT's registry reads refuse elsewhere
/// too). The caller words the refusal (`what` completes "the registry file ./.ncbirc ..."). The
/// messages of a read that does not throw are `registry_posts`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:171-177
/// ```c++
///     //set config file
///     m_ConfigFile = new CNcbiIfstream(".ncbirc");
///     m_Reg = new CNcbiRegistry(*m_ConfigFile);
///     if(m_Reg) {
///         m_TaxBrowserURL = m_Reg->Get("TAX_BROWSER","BLASTFMTUTIL");
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1615-1625
/// ```c++
/// CNcbiRegistry::CNcbiRegistry(CNcbiIstream& is, TFlags flags,
///                              const string& path)
/// ...
///     x_Init();
///     m_FileRegistry->Read(is, flags & ~(fWithNcbirc | fCaseFlags));
///     LoadBaseRegistries(flags, 0, path);
///     IncludeNcbircIfAllowed(flags & ~fCaseFlags);
/// ```
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
/// ```c++
///     const TEntries& entries = sit->second.entries;
///     TEntries::const_iterator eit = entries.find(name);
///     if (eit == entries.end()) {
///         return false;
///     } else if ((flags & fCountCleared) != 0) {
///         return true;
///     } else {
///         return !eit->second.value.empty();
///     }
/// ```
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
/// ```c++
/// string CAlignFormatUtil::GetProtocol()
/// {
///     CNcbiIfstream  config_file(".ncbirc");
///     CNcbiRegistry config_reg(config_file);
///     string httpProt = "https:";
///     if(!config_reg.Empty()) {
///         if(config_reg.HasEntry("BLASTFMTUTIL","PROTOCOL")) {
///             httpProt = config_reg.Get("BLASTFMTUTIL","PROTOCOL");
///         }
///     }
///     return httpProt;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:1236-1245, 1300-1303
/// ```c++
/// CConstRef<IRegistry> CCompoundRegistry::FindByContents(const string& section,
///                                                        const string& entry,
///                                                        TFlags flags) const
/// {
///     TFlags has_entry_flags = (flags | fCountCleared) & ~fJustCore;
///     REVERSE_ITERATE(TPriorityMap, it, m_PriorityMap) {
///         if (it->second->HasEntry(section, entry, has_entry_flags)) {
///             return it->second;
/// ...
/// bool CCompoundRegistry::x_HasEntry(const string& section, const string& name,
///                                    TFlags flags) const
/// {
///     return FindByContents(section, name, flags).NotEmpty();
/// }
/// ```
pub(crate) fn protocol(registry: Option<&RegistryRead>) -> &[u8] {
    registry
        .and_then(|read| read.last_value("BLASTFMTUTIL", "PROTOCOL"))
        .unwrap_or(super::reference::DEFAULT_PROTOCOL)
}

/// The organism report's taxonomy browser URL before its `<@protocol@>` is replaced
/// (`m_TaxBrowserURL`): `[TAX_BROWSER] BLASTFMTUTIL` (`value`), else NCBI's `TAX_BROWSER_URL`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:175-177
/// ```c++
///         m_TaxBrowserURL = m_Reg->Get("TAX_BROWSER","BLASTFMTUTIL");
///     }
///     m_TaxBrowserURL = (m_TaxBrowserURL.empty()) ? s_TagToConstString("TAX_BROWSER_URL") : m_TaxBrowserURL;
/// ```
pub(crate) fn tax_browser_template(registry: Option<&RegistryRead>) -> &[u8] {
    const TAX_BROWSER_URL: &[u8] =
        b"<@protocol@>//www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi";
    value(registry, "TAX_BROWSER", "BLASTFMTUTIL").unwrap_or(TAX_BROWSER_URL)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::ncbi_registry::{read_registry, Layer};

    #[test]
    fn the_protocol_and_the_browser_url_read_ncbirc() {
        assert_eq!(protocol(None), b"https:");
        assert_eq!(
            tax_browser_template(None),
            b"<@protocol@>//www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"
        );
        let read = read_registry(
            b"[BLASTFMTUTIL]\nPROTOCOL=http:\n[TAX_BROWSER]\nBLASTFMTUTIL=<@protocol@>//x\n",
            Layer::SystemWide,
        )
        .unwrap();
        assert_eq!(protocol(Some(&read)), b"http:");
        assert_eq!(tax_browser_template(Some(&read)), b"<@protocol@>//x");
        let read = read_registry(
            b"[BLASTFMTUTIL]\nPROTOCOL=\n[TAX_BROWSER]\nBLASTFMTUTIL=\n",
            Layer::SystemWide,
        )
        .unwrap();
        // An empty PROTOCOL is set (`HasEntry` with `fCountCleared`); an empty browser URL is not
        // (`Get(...).empty()`).
        assert_eq!(protocol(Some(&read)), b"");
        assert_eq!(value(Some(&read), "TAX_BROWSER", "BLASTFMTUTIL"), None);
    }

    #[test]
    fn a_registry_read_posts_as_ncbi() {
        // NCBI BLAST+ 2.17.0 (probe `T/t2/probe_rc`, `rep` and `bom`): the posts of one read; each
        // query's report reads `./.ncbirc` twice.
        let read = read_registry(b"[FOO]\na=1\na=2\n", Layer::SystemWide).unwrap();
        assert_eq!(
            registry_posts(&read, "blastp"),
            b"Warning: (110.8) [blastp] Found multiple [FOO] a settings; using the one from line 3\n"
                .to_vec()
        );
        let read = read_registry(b"\xef", Layer::SystemWide).unwrap();
        assert_eq!(
            registry_posts(&read, "blastp"),
            b"Error: (110.4) [blastp] Error reading the registry after line 1: \n".to_vec()
        );
    }
}
