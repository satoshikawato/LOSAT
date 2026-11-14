//! The paths in NCBI's messages as NCBI's code writes them on Windows, for the tests that compare
//! LOSAT with NCBI's bytes frozen on Linux (`blastdbcmd_info.rs`, `masker_fixtures.rs`,
//! `search_environment.rs`).

/// `text`, NCBI's bytes from Linux with the placeholders `roots` for the directories a test
/// places there, as NCBI's code writes it on Windows, where the paths of its database messages
/// are the ones `SeqDB_ConvertOSPath` makes (the separator of the platform) and the search path
/// is split at `;` (`s_GetPathSplitter`): the rest of each path after a placeholder is written
/// with `\`, so is the name in `database [name]`, and the `:` of a search path `[<root>::]` is
/// `;`.
pub(crate) fn ncbi_paths_on_windows(text: &str, roots: &[&str]) -> String {
    let mut text = text.to_string();
    for root in roots {
        text = text.replace(&format!("{root}::]"), &format!("{root};;]"));
    }
    // The database name of "No alias or index file found for ... database [name]" goes through
    // the same conversion (`SeqDB_ConvertOSPath`).
    let mut converted = String::new();
    let mut rest = text.as_str();
    while let Some(at) = rest.find("database [") {
        let name_at = at + "database [".len();
        converted.push_str(&rest[..name_at]);
        let end = rest[name_at..]
            .find(']')
            .map_or(rest.len(), |n| name_at + n);
        converted.push_str(&rest[name_at..end].replace('/', "\\"));
        rest = &rest[end..];
    }
    converted.push_str(rest);
    let mut out = String::new();
    let mut rest = converted.as_str();
    while let Some((at, root)) = roots
        .iter()
        .filter_map(|root| rest.find(root).map(|at| (at, root)))
        .min()
    {
        let path_at = at + root.len();
        out.push_str(&rest[..path_at]);
        rest = &rest[path_at..];
        let end = rest
            .find(|c: char| !(c.is_ascii_alphanumeric() || "_./-".contains(c)))
            .unwrap_or(rest.len());
        out.push_str(&rest[..end].replace('/', "\\"));
        rest = &rest[end..];
    }
    out.push_str(rest);
    out
}

#[test]
fn ncbi_s_paths_are_written_as_on_windows() {
    let linux = "Error: (CSeqDBException::eInvalid) No alias or index file found for protein \
                 database [db/p 5] in search path [{R}::]. Please verify\n\
                 {R}/db/n5\t{T}/seq/a.fa {R}";
    assert_eq!(
        ncbi_paths_on_windows(linux, &["{R}", "{T}"]),
        "Error: (CSeqDBException::eInvalid) No alias or index file found for protein \
         database [db\\p 5] in search path [{R};;]. Please verify\n\
         {R}\\db\\n5\t{T}\\seq\\a.fa {R}"
    );
}
