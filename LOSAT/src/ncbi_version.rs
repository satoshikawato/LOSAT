//! NCBI corelib/version.cpp: what an NCBI application prints for `-version`, `-version-full`,
//! `-version-full-xml` and `-version-full-json` (`CVersionAPI::Print`, `PrintXml`, `PrintJson`),
//! the words that its application start reads before the arguments (`ncbi_environment::preparse`).
//!
//! The forms are NCBI's; the values are LOSAT's (line RG DECISIONS E3-1). The first line (the
//! program and its version, `blastn: 2.17.0+`) is NCBI's. NCBI's package is the BLAST+ build,
//! LOSAT's is LOSAT: the text's `Package:` line names LOSAT and the BLAST+ release its output
//! follows, and the XML and JSON name LOSAT and its version. What NCBI's build writes about
//! itself (`Package-Config`'s configure line, `Build-Signature`, the build dates, `Revision`,
//! `Stable-Components-Version`, `Production-Version`, `Built-As`, the toolkit component of
//! NCBI's rpsbproc) has no LOSAT value and is not written: the keys that NCBI writes whenever
//! there is a package (`Package-Config`, the package's and the application's build information)
//! stay, empty, as NCBI's bioconda build writes an empty `Package-Config`. `GI-64bit` is
//! `TRUE`: LOSAT reads GIs as NCBI's default 64-bit build (`ParsedSeqId::Gi(i64)`).

/// The word of NCBI's application start that prints a version and ends the program with exit
/// status 0 (`ncbi_environment::preparse`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:84-87
/// ```c++
/// const char* s_ArgVersion         = "-version";
/// const char* s_ArgFullVersion     = "-version-full";
/// const char* s_ArgFullVersionXml  = "-version-full-xml";
/// const char* s_ArgFullVersionJson = "-version-full-json";
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VersionRequest {
    /// `-version`: `Print(appname, fVersionInfo | fPackageShort)`.
    Version,
    /// `-version-full`: `Print(appname)`, every part (`fPrintAll`).
    Full,
    /// `-version-full-xml`: `PrintXml(appname)`.
    FullXml,
    /// `-version-full-json`: `PrintJson(appname)`.
    FullJson,
}

impl VersionRequest {
    /// The request of the command-line word `word`, which NCBI compares whole (`strcmp`).
    pub fn of_word(word: &[u8]) -> Option<Self> {
        match word {
            b"-version" => Some(Self::Version),
            b"-version-full" => Some(Self::Full),
            b"-version-full-xml" => Some(Self::FullXml),
            b"-version-full-json" => Some(Self::FullJson),
            _ => None,
        }
    }

    /// The command-line word of the request.
    pub fn word(self) -> &'static str {
        match self {
            Self::Version => "-version",
            Self::Full => "-version-full",
            Self::FullXml => "-version-full-xml",
            Self::FullJson => "-version-full-json",
        }
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:972-1000
// ```c++
//             } else if ( NStr::strcmp(argv[i], s_ArgVersion) == 0 ) {
//                 delete[] v;
//                 // Print VERSION
//                 cout << GetFullVersion().Print( appname, CVersionAPI::fVersionInfo | CVersionAPI::fPackageShort );
//                 diag_context.DiscardMessages();
//                 return 0;
//
//                 // Full version
//             } else if ( NStr::strcmp(argv[i], s_ArgFullVersion) == 0 ) {
//                 delete[] v;
//                 // Print full VERSION
//                 cout << GetFullVersion().Print( appname );
//                 diag_context.DiscardMessages();
//                 return 0;
//             } else if ( NStr::strcmp(argv[i], s_ArgFullVersionXml) == 0 ) {
//                 delete[] v;
//                 // Print full VERSION in XML format
//                 cout << GetFullVersion().PrintXml( appname );
//                 diag_context.DiscardMessages();
//                 return 0;
//             } else if ( NStr::strcmp(argv[i], s_ArgFullVersionJson) == 0 ) {
//                 delete[] v;
//                 // Print full VERSION in JSON format
//                 cout << GetFullVersion().PrintJson( appname );
//                 diag_context.DiscardMessages();
//                 return 0;
// ```
/// The text that `request` prints for the program `appname` whose version prints as `version`
/// (`CVersionInfo::Print`: `2.17.0+` for the BLAST+ programs, whose `CBlastVersion::Print` adds
/// the `+`, `1.0.0` for dustmasker and segmasker, `0.5.1` for rpsbproc).
pub fn print(request: VersionRequest, appname: &str, version: &str) -> String {
    match request {
        VersionRequest::Version => format!("{appname}: {version}\n{}", package_line()),
        VersionRequest::Full => print_text(appname, version),
        VersionRequest::FullXml => print_xml(appname, version),
        VersionRequest::FullJson => print_json(appname, version),
    }
}

/// The package line of the text forms (`fPackageShort`). NCBI's names its package and the
/// package's build date; LOSAT's names LOSAT and the BLAST+ release its output follows.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:737-742
/// ```c++
///     if (flags & ( fPackageShort | fPackageFull )) {
///         os << " Package: " << GetPackageName() << ' '
///            << GetPackageVersion().Print() << ", build "
///            << NCBI_SBUILDINFO_DEFAULT().date
///            << endl;
///     }
/// ```
fn package_line() -> String {
    format!(
        " Package: LOSAT {}, output compatible with blast {}\n",
        env!("CARGO_PKG_VERSION"),
        crate::NCBI_BLAST_RELEASE
    )
}

/// `CVersionAPI::Print(appname, fPrintAll)` with LOSAT's values: no components, an empty
/// `Package-Config`, no `Build-Signature` (NCBI writes it only for a build that has one), 64-bit
/// GIs and an empty build information (`SBuildInfo::Print` writes a line per value).
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:743-769
/// ```c++
///     if (flags & fPackageFull) {
///         os << " Package-Config: " << ' ' << GetPackageConfig() << endl;
///     }
/// #endif
///
/// #ifdef NCBI_SIGNATURE
///     if (flags & fBuildSignature) {
///         os << " Build-Signature: " << ' ' << NCBI_SIGNATURE << endl;
///     }
/// #endif
///     if (flags & fGI64bit) {
/// #ifdef NCBI_INT8_GI
///         os << " GI-64bit:  TRUE" << endl;
/// #else
///         os << " GI-64bit:  FALSE" << endl;
/// #endif
///     }
///
///     if (flags & fBuildInfo) {
///         os << m_BuildInfo.Print();
///     }
///
///     if (flags & fComponents) {
///         for (const auto& c : m_Components) {
///             os << endl << ' ' <<  c->Print() << endl;
///         }
///     }
/// ```
fn print_text(appname: &str, version: &str) -> String {
    format!(
        "{appname}: {version}\n{} Package-Config:  \n GI-64bit:  TRUE\n",
        package_line()
    )
}

/// `CVersionAPI::PrintXml(appname, fPrintAll)` with LOSAT's values (`print_text`): the package
/// is LOSAT, with its version, an empty build information and an empty `config`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:778-820
/// ```c++
///     os << "<?xml version=\"1.0\"?>\n"
///         "<ncbi_version xmlns=\"ncbi:version\"\n"
///         "  xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
///         "  xs:schemaLocation=\"ncbi:version ncbi_version.xsd\">\n";
///
///     if (flags & fVersionInfo) {
///         if ( !appname.empty() ) {
///             os << "<appname>" << NStr::XmlEncode(appname) << "</appname>\n";
///         }
///         os << m_VersionInfo->PrintXml();
///     }
///     ...
///     if (flags & ( fPackageShort | fPackageFull )) {
///         os << "<package name=\"" << NStr::XmlEncode(GetPackageName()) << "\">\n" <<
///             GetPackageVersion().PrintXml() <<
///             NCBI_SBUILDINFO_DEFAULT().PrintXml();
///         if (flags & fPackageFull) {
///             os << "<config>" << NStr::XmlEncode(GetPackageConfig()) << "</config>\n";
///         }
///         os << "</package>\n";
///     }
///     ...
///     if (flags & fBuildInfo) {
///         os << m_BuildInfo.PrintXml();
///     }
///
///     os << "</ncbi_version>\n";
/// ```
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:567-579 (`SBuildInfo::PrintXml`,
/// without a date, tag or other value)
/// ```c++
///     os << "<build_info";
///     ...
///     os << ">" << endl;
///     ...
///     os << "</build_info>" << endl;
/// ```
fn print_xml(appname: &str, version: &str) -> String {
    format!(
        "<?xml version=\"1.0\"?>\n\
         <ncbi_version xmlns=\"ncbi:version\"\n  \
         xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n  \
         xs:schemaLocation=\"ncbi:version ncbi_version.xsd\">\n\
         <appname>{appname}</appname>\n\
         {}\
         <package name=\"LOSAT\">\n\
         {}\
         <build_info>\n\
         </build_info>\n\
         <config></config>\n\
         </package>\n\
         <build_info>\n\
         </build_info>\n\
         </ncbi_version>\n",
        version_info_xml(&VersionNumbers::of_print(version)),
        version_info_xml(&VersionNumbers::losat()),
    )
}

/// `CVersionAPI::PrintJson(appname, fPrintAll)` with LOSAT's values (`print_xml`): an empty
/// list of components, as NCBI's BLAST+ programs write.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:829-884
/// ```c++
///     os << "{\n  \"ncbi_version\": {\n";
///
///     if (flags & fVersionInfo) {
///         if ( !appname.empty() ) {
///             os << "    \"appname\": " << NStr::JsonEncode(appname, NStr::eJsonEnc_Quoted) << ",\n";
///         }
///         os << "    \"version_info\": " << m_VersionInfo->PrintJson();
///         need_separator = true;
///     }
///
///     if (flags & fComponents) {
///         if ( need_separator ) os << ",\n";
///         os << "    \"component\": [";
///         ...
///         os << "]";
///         need_separator = true;
///     }
///
/// #if NCBI_PACKAGE
///     if (flags & ( fPackageShort | fPackageFull )) {
///         if ( need_separator ) os << ",\n";
///         os << "    \"package\": {\n" <<
///             "      \"name\": " << NStr::JsonEncode(GetPackageName(), NStr::eJsonEnc_Quoted) << ",\n" <<
///             "      \"version_info\": " << GetPackageVersion().PrintJson() << ",\n" <<
///             "      \"build_info\": " << NCBI_SBUILDINFO_DEFAULT().PrintJson();
///         if (flags & fPackageFull) {
///             os << ",\n      \"config\": " << NStr::JsonEncode(GetPackageConfig(), NStr::eJsonEnc_Quoted);
///         }
///         os << "}";
///         need_separator = true;
///     }
/// #endif
///     ...
///     if (flags & fBuildInfo) {
///         if ( need_separator ) os << ",\n";
///         os << "    \"build_info\": " << m_BuildInfo.PrintJson();
///         need_separator = true;
///     }
///
///     os << "\n  }\n}\n";
/// ```
fn print_json(appname: &str, version: &str) -> String {
    format!(
        "{{\n  \"ncbi_version\": {{\n    \
         \"appname\": \"{appname}\",\n    \
         \"version_info\": {},\n    \
         \"component\": [],\n    \
         \"package\": {{\n      \
         \"name\": \"LOSAT\",\n      \
         \"version_info\": {},\n      \
         \"build_info\": {{}},\n      \
         \"config\": \"\"}},\n    \
         \"build_info\": {{}}\n  }}\n}}\n",
        version_info_json(&VersionNumbers::of_print(version)),
        version_info_json(&VersionNumbers::losat()),
    )
}

/// The parts of a `CVersionInfo` that its XML and JSON write: the numbers, and LOSAT's
/// pre-release text (`CARGO_PKG_VERSION_PRE`) as the version's name.
struct VersionNumbers<'a> {
    major: &'a str,
    minor: &'a str,
    patch: &'a str,
    name: &'a str,
}

impl<'a> VersionNumbers<'a> {
    /// The numbers of `version` as `CVersionInfo::Print` writes it (`2.17.0+`: `CBlastVersion`
    /// adds the `+` to the text, not to the numbers).
    fn of_print(version: &'a str) -> Self {
        let numbers = version.strip_suffix('+').unwrap_or(version);
        let mut parts = numbers.splitn(3, '.');
        let mut next = || parts.next().unwrap_or("0");
        Self {
            major: next(),
            minor: next(),
            patch: next(),
            name: "",
        }
    }

    /// LOSAT's version, the package's.
    fn losat() -> VersionNumbers<'static> {
        VersionNumbers {
            major: env!("CARGO_PKG_VERSION_MAJOR"),
            minor: env!("CARGO_PKG_VERSION_MINOR"),
            patch: env!("CARGO_PKG_VERSION_PATCH"),
            name: env!("CARGO_PKG_VERSION_PRE"),
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:137-153
/// ```c++
/// string CVersionInfo::PrintXml(void) const
/// {
///     CNcbiOstrstream os;
///     os << "<version_info";
///     if (m_Major >= 0) {
///         os << " major=\"" << m_Major <<
///             "\" minor=\"" << (m_Minor >= 0 ? m_Minor : 0) << "\"";
///         if (m_PatchLevel >= 0) {
///             os << " patch_level=\"" << m_PatchLevel << "\"";
///         }
///     }
///     if ( !m_Name.empty() ) {
///         os << " name=\"" << NStr::XmlEncode(m_Name) << "\"";
///     }
///     os << "/>\n";
/// ```
fn version_info_xml(version: &VersionNumbers) -> String {
    let VersionNumbers {
        major,
        minor,
        patch,
        name,
    } = version;
    let name = if name.is_empty() {
        String::new()
    } else {
        format!(" name=\"{name}\"")
    };
    format!("<version_info major=\"{major}\" minor=\"{minor}\" patch_level=\"{patch}\"{name}/>\n")
}

/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:156-175
/// ```c++
/// string CVersionInfo::PrintJson(void) const
/// {
///     CNcbiOstrstream os;
///     os << "{";
///     bool need_separator = false;
///     if (m_Major >= 0) {
///         os << "\"major\": " << m_Major <<
///             ", \"minor\": " << (m_Minor >= 0 ? m_Minor : 0);
///         if (m_PatchLevel >= 0) {
///             os << ", \"patch_level\": " << m_PatchLevel;
///         }
///         need_separator = true;
///     }
///     if ( !m_Name.empty() ) {
///         if ( need_separator ) os << ", ";
///         os << "\"name\": " << NStr::JsonEncode(m_Name, NStr::eJsonEnc_Quoted);
///     }
///     os << "}";
/// ```
fn version_info_json(version: &VersionNumbers) -> String {
    let VersionNumbers {
        major,
        minor,
        patch,
        name,
    } = version;
    let name = if name.is_empty() {
        String::new()
    } else {
        format!(", \"name\": \"{name}\"")
    };
    format!("{{\"major\": {major}, \"minor\": {minor}, \"patch_level\": {patch}{name}}}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_request_is_the_whole_word() {
        for request in [
            VersionRequest::Version,
            VersionRequest::Full,
            VersionRequest::FullXml,
            VersionRequest::FullJson,
        ] {
            assert_eq!(
                VersionRequest::of_word(request.word().as_bytes()),
                Some(request)
            );
        }
        for word in ["-version-full=1", "-version-fullx", "--version", "-Version"] {
            assert_eq!(VersionRequest::of_word(word.as_bytes()), None, "{word}");
        }
    }

    // NCBI BLAST+ 2.17.0 (the x64-linux release) `blastn -version-full`, with LOSAT's package and
    // without the lines of NCBI's build (line RG item E3, `E/e1/ncbi/blastn.vfull.out`).
    #[test]
    fn the_full_version_is_ncbi_s_form_with_losat_s_values() {
        let package = format!(
            " Package: LOSAT {}, output compatible with blast 2.17.0\n",
            env!("CARGO_PKG_VERSION")
        );
        assert_eq!(
            print(VersionRequest::Version, "blastn", "2.17.0+"),
            format!("blastn: 2.17.0+\n{package}")
        );
        assert_eq!(
            print(VersionRequest::Full, "blastn", "2.17.0+"),
            format!("blastn: 2.17.0+\n{package} Package-Config:  \n GI-64bit:  TRUE\n")
        );
        let pre = env!("CARGO_PKG_VERSION_PRE");
        let (xml_name, json_name) = if pre.is_empty() {
            (String::new(), String::new())
        } else {
            (format!(" name=\"{pre}\""), format!(", \"name\": \"{pre}\""))
        };
        let (major, minor, patch) = (
            env!("CARGO_PKG_VERSION_MAJOR"),
            env!("CARGO_PKG_VERSION_MINOR"),
            env!("CARGO_PKG_VERSION_PATCH"),
        );
        assert_eq!(
            print(VersionRequest::FullXml, "rpsbproc", "0.5.1"),
            format!(
                "<?xml version=\"1.0\"?>
<ncbi_version xmlns=\"ncbi:version\"
  xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"
  xs:schemaLocation=\"ncbi:version ncbi_version.xsd\">
<appname>rpsbproc</appname>
<version_info major=\"0\" minor=\"5\" patch_level=\"1\"/>
<package name=\"LOSAT\">
<version_info major=\"{major}\" minor=\"{minor}\" patch_level=\"{patch}\"{xml_name}/>
<build_info>
</build_info>
<config></config>
</package>
<build_info>
</build_info>
</ncbi_version>
"
            )
        );
        assert_eq!(
            print(VersionRequest::FullJson, "dustmasker", "1.0.0"),
            format!(
                "{{
  \"ncbi_version\": {{
    \"appname\": \"dustmasker\",
    \"version_info\": {{\"major\": 1, \"minor\": 0, \"patch_level\": 0}},
    \"component\": [],
    \"package\": {{
      \"name\": \"LOSAT\",
      \"version_info\": {{\"major\": {major}, \"minor\": {minor}, \"patch_level\": {patch}{json_name}}},
      \"build_info\": {{}},
      \"config\": \"\"}},
    \"build_info\": {{}}
  }}
}}
"
            )
        );
    }

    // The BLAST+ programs' version has the `+` of `CBlastVersion::Print` in the text only.
    #[test]
    fn the_numbers_leave_out_the_plus() {
        let numbers = VersionNumbers::of_print(crate::NCBI_BLAST_VERSION);
        assert_eq!(
            (numbers.major, numbers.minor, numbers.patch),
            ("2", "17", "0")
        );
        assert!(print(
            VersionRequest::FullXml,
            "makeblastdb",
            crate::NCBI_BLAST_VERSION
        )
        .contains("<version_info major=\"2\" minor=\"17\" patch_level=\"0\"/>\n<package"));
    }
}
