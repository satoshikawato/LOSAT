//! The paths in NCBI's messages as NCBI's code writes them on Windows, for the tests that compare
//! LOSAT with NCBI's bytes frozen on Linux (the `blastdbcmd_*.rs` tests, `masker_fixtures.rs`,
//! the makembindex tests), and the directories those tests run in.

use std::fs;
use std::path::{Path, PathBuf};

/// `fs::canonicalize` without the `\\?\` prefix that Windows puts on its result: a verbatim path
/// takes `/` as a character of a name (`{O}/B23.out` is not a file of the directory `{O}` then),
/// and the product prints the paths it lists without the prefix.
pub(crate) fn canonical(path: &Path) -> PathBuf {
    let path = fs::canonicalize(path).unwrap();
    #[cfg(windows)]
    if let Some(plain) = path.to_str().and_then(|text| text.strip_prefix(r"\\?\")) {
        if plain.as_bytes().get(1) == Some(&b':') {
            return PathBuf::from(plain);
        }
    }
    path
}

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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1712-1722
/// `text`, NCBI's `-metadata` JSON from Linux for the `-db` value `db`, as NCBI's code writes it
/// on Windows: `"dbname"` holds each name of the list (`db` split at spaces) after its last
/// separator of the platform (`CFile::GetPathSeparator()`, `\` there), so a name given with `/`
/// stays whole; the names of `"files"` are `CFile::GetName()` on both platforms and do not
/// change. The value is a JSON string (a `"` or a `\` of a name is escaped).
pub(crate) fn ncbi_metadata_on_windows(text: &str, db: &str) -> String {
    let names: Vec<&str> = db
        .split(' ')
        .filter(|name| !name.is_empty())
        .map(|name| name.rsplit('\\').next().unwrap_or(name))
        .collect();
    let escaped = names.join(" ").replace('\\', "\\\\").replace('"', "\\\"");
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        if line.starts_with("  \"dbname\": \"") {
            out.push_str(&format!("  \"dbname\": \"{escaped}\",\n"));
        } else {
            out.push_str(line);
        }
    }
    out
}

/// `text`, NCBI's bytes from Linux for the blastdbcmd command line `words` with the placeholders
/// `roots`, as NCBI's code writes them on the platform of the run: as they are off Windows; on
/// Windows with its paths (`ncbi_paths_on_windows`) and, for `-metadata`, its database names
/// (`ncbi_metadata_on_windows`). A text that is not UTF-8 holds neither.
pub(crate) fn blastdbcmd_bytes_on_platform(text: &[u8], roots: &[&str], words: &[&str]) -> Vec<u8> {
    let text = match std::str::from_utf8(text) {
        Ok(text) if cfg!(windows) => text,
        _ => return text.to_vec(),
    };
    let text = ncbi_paths_on_windows(text, roots);
    let db = words
        .iter()
        .position(|word| *word == "-db")
        .and_then(|at| words.get(at + 1));
    match db {
        Some(db) if words.contains(&"-metadata") => ncbi_metadata_on_windows(&text, db),
        _ => text,
    }
    .into_bytes()
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

#[test]
fn a_metadata_dbname_keeps_a_name_given_with_slashes_on_windows() {
    let linux = "{\n  \"version\": \"1.2\",\n  \"dbname\": \"p5\\\" p5n\\\"\",\n  \"files\": [\n    \"p5.pdb\"\n";
    assert_eq!(
        ncbi_metadata_on_windows(linux, "\"j/lst/p5\" \"j/lst/p5n\""),
        "{\n  \"version\": \"1.2\",\n  \"dbname\": \"\\\"j/lst/p5\\\" \\\"j/lst/p5n\\\"\",\n  \"files\": [\n    \"p5.pdb\"\n"
    );
    assert_eq!(
        ncbi_metadata_on_windows("  \"dbname\": \"p5\",\n", r"a\b\p5  db/p5n"),
        "  \"dbname\": \"p5 db/p5n\",\n"
    );
}
