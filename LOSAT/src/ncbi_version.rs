//! NCBI reference (598d8ae6): c++/src/corelib/version.cpp: what an NCBI application prints for `-version`, `-version-full`,
//! `-version-full-xml` and `-version-full-json` (`CVersionAPI::Print`, `PrintXml`, `PrintJson`),
//! the words that its application start reads before the arguments (`ncbi_environment::preparse`).
//!
//! The forms are NCBI's; the values are LOSAT's (line RG DECISIONS E3-1). The first line names the
//! program and LOSAT's version (`blastn: LOSAT 0.2.0`, `wording::LOSAT_VERSION`; NCBI's
//! `blastn: 2.17.0+`, the Owner 2026-11-21 「Q1=B」), as do the XML and JSON version of the
//! program. NCBI's package is the BLAST+ build,
//! LOSAT's is LOSAT: the text's `Package:` line names LOSAT and the BLAST+ release its output
//! follows, and the XML and JSON name LOSAT and its version. What NCBI's build writes about
//! itself (`Package-Config`'s configure line, `Build-Signature`, the build dates, `Revision`,
//! `Stable-Components-Version`, `Production-Version`, `Built-As`, the toolkit component of
//! NCBI's rpsbproc) has no LOSAT value and is not written: the keys that NCBI writes whenever
//! there is a package (`Package-Config`, the package's and the application's build information)
//! stay, empty, as NCBI's bioconda build writes an empty `Package-Config`. `GI-64bit` is
//! `TRUE`: LOSAT reads GIs as NCBI's default 64-bit build (`ParsedSeqId::Gi(i64)`).
//! windowmasker's components (`Component`: its algorithm and its counts formats, which LOSAT
//! implements) are written with their versions and an empty build information.

/// The word of NCBI's application start that prints a version and ends the program with exit
/// status 0 (`ncbi_environment::preparse`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:84-87
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
/// The text that `request` prints for the program `appname`: LOSAT's version
/// (`wording::LOSAT_VERSION`) where NCBI's `CVersionInfo::Print` writes the program's own
/// (`2.17.0+` for the BLAST+ programs, `1.0.0` for dustmasker and segmasker, `0.5.1` for rpsbproc).
pub fn print(request: VersionRequest, appname: &str) -> String {
    print_with_components(request, appname, &[])
}

/// A component of an application's version (`CVersionAPI::AddComponentVersion`): its name,
/// the text before its version in the text form (`CSeqMaskerVersion`'s `ver_pfx`) and its
/// version. Its build information is NCBI's build's, which LOSAT does not write (empty, as the
/// application's).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Component {
    pub name: &'static str,
    pub prefix: &'static str,
    pub major: u32,
    pub minor: u32,
    pub patch: u32,
}

/// `print` of an application that has `components` (windowmasker's six).
pub fn print_with_components(
    request: VersionRequest,
    appname: &str,
    components: &[Component],
) -> String {
    let version = crate::wording::LOSAT_VERSION;
    match request {
        VersionRequest::Version => format!("{appname}: {version}\n{}", package_line()),
        VersionRequest::Full => {
            let mut text = print_text(appname, version);
            for component in components {
                text.push_str(&component_text(component));
            }
            text
        }
        VersionRequest::FullXml => print_xml(appname, components),
        VersionRequest::FullJson => print_json(appname, components),
    }
}

impl Component {
    fn numbers(&self) -> (String, String, String) {
        (
            self.major.to_string(),
            self.minor.to_string(),
            self.patch.to_string(),
        )
    }
}

/// A component in the text form: `CSeqMaskerVersion::Print` after an empty line.
///
/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_version.hpp:63-67
/// (`CVersionAPI::Print` writes a line break, a space, the component's own text and a line break;
/// `print_text`.)
fn component_text(component: &Component) -> String {
    format!(
        "\n {}:{}{}.{}.{}\n",
        component.name, component.prefix, component.major, component.minor, component.patch
    )
}

/// A component in the XML form, with LOSAT's empty build information.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:432-440
fn component_xml(component: &Component) -> String {
    let (major, minor, patch) = component.numbers();
    format!(
        "<component name=\"{}\">\n{}\n<build_info>\n</build_info>\n</component>\n",
        component.name,
        version_info_xml(&VersionNumbers {
            major: &major,
            minor: &minor,
            patch: &patch,
            name: "",
        })
    )
}

/// The JSON list of components, with LOSAT's empty build information.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:442-449,840-849
fn components_json(components: &[Component]) -> String {
    let items: Vec<String> = components
        .iter()
        .map(|component| {
            let (major, minor, patch) = component.numbers();
            format!(
                "\n      {{ \"name\": \"{}\", \"version_info\": {},\n        \"build_info\": {{}}}}",
                component.name,
                version_info_json(&VersionNumbers {
                    major: &major,
                    minor: &minor,
                    patch: &patch,
                    name: "",
                })
            )
        })
        .collect();
    format!("[{}]", items.join(","))
}

/// The package line of the text forms (`fPackageShort`). NCBI's names its package and the
/// package's build date; LOSAT's names LOSAT and the BLAST+ release its output follows.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:737-742
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
fn print_text(appname: &str, version: &str) -> String {
    format!(
        "{appname}: {version}\n{} Package-Config:  \n GI-64bit:  TRUE\n",
        package_line()
    )
}

/// `CVersionAPI::PrintXml(appname, fPrintAll)` with LOSAT's values (`print_text`): the program's
/// version and the package's are LOSAT's, with an empty build information and an empty `config`.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:778-820
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:567-579 (`SBuildInfo::PrintXml`,
/// without a date, tag or other value)
fn print_xml(appname: &str, components: &[Component]) -> String {
    format!(
        "<?xml version=\"1.0\"?>\n\
         <ncbi_version xmlns=\"ncbi:version\"\n  \
         xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n  \
         xs:schemaLocation=\"ncbi:version ncbi_version.xsd\">\n\
         <appname>{appname}</appname>\n\
         {}\
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
        version_info_xml(&VersionNumbers::losat()),
        components.iter().map(component_xml).collect::<String>(),
        version_info_xml(&VersionNumbers::losat()),
    )
}

/// `CVersionAPI::PrintJson(appname, fPrintAll)` with LOSAT's values (`print_xml`): the list of
/// components is empty, as NCBI's BLAST+ programs write it, except for windowmasker.
///
/// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:829-884
fn print_json(appname: &str, components: &[Component]) -> String {
    format!(
        "{{\n  \"ncbi_version\": {{\n    \
         \"appname\": \"{appname}\",\n    \
         \"version_info\": {},\n    \
         \"component\": {},\n    \
         \"package\": {{\n      \
         \"name\": \"LOSAT\",\n      \
         \"version_info\": {},\n      \
         \"build_info\": {{}},\n      \
         \"config\": \"\"}},\n    \
         \"build_info\": {{}}\n  }}\n}}\n",
        version_info_json(&VersionNumbers::losat()),
        components_json(components),
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

impl VersionNumbers<'_> {
    /// LOSAT's version, the program's and the package's.
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
#[path = "ncbi_version_tests.rs"]
mod tests;
