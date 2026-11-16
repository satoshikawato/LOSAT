//! NCBI corelib/ncbiargs.cpp: what an NCBI application prints for `-xmlhelp`
//! (`CArgDescriptions::PrintUsageXml`), the description of its arguments in XML.
//!
//! LOSAT prints NCBI BLAST+ 2.17.0's XML of each program, frozen (`ncbi_xmlhelp/<program>.xml`),
//! without the arguments that LOSAT refuses by name and the dependencies that name them (line RG
//! DECISIONS E2x-1, the Owner's choice of 2026-11-13): the arguments that LOSAT runs keep NCBI's
//! definitions, byte for byte (their types, defaults, constraints, groups and dependencies), and
//! an argument that LOSAT ports is written again as it stops being refused. The names to leave
//! out are the program reader's own refusals (`cli::refuses_option`; the readers of
//! `ncbi_app` refuse none), so the XML describes what the reader takes.

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:685-690
// ```c++
//     catch (const CArgHelpException& e) {
//         x_AddDefaultArgs();
//         // Print USAGE
//         if (e.GetErrCode() == CArgHelpException::eHelpXml) {
//             m_ArgDesc->PrintUsageXml(cout);
//         } else {
// ```
// The XML of each program as NCBI BLAST+ 2.17.0 prints it (the official x64-linux build; the
// bioconda build prints the same bytes; rpsbproc 0.5.1 of bioconda; line RG item E2, probe `xml`).
const PROGRAMS: &[(&str, &str)] = &[
    ("blastn", include_str!("ncbi_xmlhelp/blastn.xml")),
    ("blastp", include_str!("ncbi_xmlhelp/blastp.xml")),
    ("blastx", include_str!("ncbi_xmlhelp/blastx.xml")),
    ("tblastn", include_str!("ncbi_xmlhelp/tblastn.xml")),
    ("tblastx", include_str!("ncbi_xmlhelp/tblastx.xml")),
    ("psiblast", include_str!("ncbi_xmlhelp/psiblast.xml")),
    ("rpsblast", include_str!("ncbi_xmlhelp/rpsblast.xml")),
    ("rpstblastn", include_str!("ncbi_xmlhelp/rpstblastn.xml")),
    ("makeblastdb", include_str!("ncbi_xmlhelp/makeblastdb.xml")),
    ("rpsbproc", include_str!("ncbi_xmlhelp/rpsbproc.xml")),
    ("blastdbcmd", include_str!("ncbi_xmlhelp/blastdbcmd.xml")),
    ("dustmasker", include_str!("ncbi_xmlhelp/dustmasker.xml")),
    ("segmasker", include_str!("ncbi_xmlhelp/segmasker.xml")),
    (
        "convert2blastmask",
        include_str!("ncbi_xmlhelp/convert2blastmask.xml"),
    ),
];

/// NCBI 2.17.0's `-xmlhelp` of `program`, unchanged.
fn ncbi_xml(program: &str) -> Option<&'static str> {
    PROGRAMS
        .iter()
        .find(|(name, _)| *name == program)
        .map(|(_, xml)| *xml)
}

/// What `losat <program> -xmlhelp` prints: NCBI's XML of `program` without the arguments for
/// which `refused` is true (`-name` given as `name`) and the dependencies that name one; `None`
/// for a program that NCBI does not have.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3900-3943
/// ```c++
/// // keys
///     ITERATE (TArgs, a, desc.m_Args) {
///         if (s_IsKey(**a)) {
///             tag = (*a)->PrintXml(m_out);
///             desc.x_PrintAliasesAsXml(m_out, (*a)->GetName());
///             m_out << "</" << tag << ">" << endl;
///         }
///     }
/// ...
///     if (!desc.m_Dependencies.empty()) {
///         m_out << "<" << "dependencies" << ">" << endl;
///         ITERATE(TDependencies, dep, desc.m_Dependencies) {
///             if (dep->second.m_Dep == eRequires) {
///                 m_out << "<" << "first_requires_second" << ">" << endl;
///                 s_WriteXmlLine(m_out, "arg1", dep->first);
///                 s_WriteXmlLine(m_out, "arg2", dep->second.m_Arg);
///                 m_out << "</" << "first_requires_second" << ">" << endl;
///             }
///         }
/// ...
///         m_out << "</" << "dependencies" << ">" << endl;
///     }
/// ```
/// Each argument is an element of its own lines, from `<key name="NAME" ...>` (or `<flag`) to
/// `</key>` (`</flag>`), and each dependency the four lines of its two names; NCBI writes
/// `<dependencies>` only when there is one, so the block goes when none is left. The names are
/// written as NCBI's reader takes them (no character of a name needs escaping).
pub(crate) fn print(program: &str, refused: impl Fn(&str) -> bool) -> Option<String> {
    let xml = ncbi_xml(program)?;
    let mut out = String::with_capacity(xml.len());
    let mut lines = xml.split_inclusive('\n');
    while let Some(line) = lines.next() {
        if let Some((tag, name)) = argument_start(line) {
            let end = format!("</{tag}>\n");
            let mut element = vec![line];
            for line in lines.by_ref() {
                element.push(line);
                if line == end {
                    break;
                }
            }
            if !refused(name) {
                out.extend(element);
            }
            continue;
        }
        if let Some(tag) = line
            .strip_prefix('<')
            .and_then(|rest| rest.strip_suffix(">\n"))
            .filter(|tag| matches!(*tag, "first_requires_second" | "first_excludes_second"))
        {
            let element: Vec<&str> = std::iter::once(line)
                .chain(lines.by_ref().take(3))
                .collect();
            let names_refused = element[1..3].iter().any(|line| {
                xml_line_value(line)
                    .map(|(_, name)| refused(name))
                    .unwrap_or(false)
            });
            debug_assert_eq!(element.get(3).copied(), Some(&*format!("</{tag}>\n")));
            if !names_refused {
                out.extend(element);
            }
            continue;
        }
        out.push_str(line);
    }
    Some(out.replace("<dependencies>\n</dependencies>\n", ""))
}

/// The tag (`key` or `flag`) and the name of the argument whose element starts at `line`
/// (`CArgDesc::PrintXml`: `<key name="NAME" type=...>`).
fn argument_start(line: &str) -> Option<(&str, &str)> {
    let rest = line.strip_prefix('<')?;
    let (tag, rest) = rest.split_once(' ')?;
    if !matches!(tag, "key" | "flag") {
        return None;
    }
    let name = rest.strip_prefix("name=\"")?.split_once('"')?.0;
    Some((tag, name))
}

/// The tag and the text of a line `<tag>text</tag>` (`s_WriteXmlLine`).
fn xml_line_value(line: &str) -> Option<(&str, &str)> {
    let (tag, rest) = line.strip_prefix('<')?.split_once('>')?;
    let text = rest
        .strip_suffix('\n')?
        .strip_suffix(&*format!("</{tag}>"))?;
    Some((tag, text))
}

/// The names of the arguments of NCBI's XML of `program` (its keys and flags, in NCBI's order).
#[cfg(test)]
pub(crate) fn ncbi_argument_names(program: &str) -> Vec<(&'static str, bool)> {
    ncbi_xml(program)
        .into_iter()
        .flat_map(|xml| xml.split_inclusive('\n'))
        .filter_map(argument_start)
        .map(|(tag, name)| (name, tag == "flag"))
        .collect()
}

/// The programs whose XML LOSAT has.
#[cfg(test)]
pub(crate) fn programs() -> impl Iterator<Item = &'static str> {
    PROGRAMS.iter().map(|(name, _)| *name)
}

#[cfg(test)]
#[path = "ncbi_xmlhelp_tests.rs"]
mod tests;
