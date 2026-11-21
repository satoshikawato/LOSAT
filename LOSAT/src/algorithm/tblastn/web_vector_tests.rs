use std::path::PathBuf;

pub(crate) struct Case {
    pub(crate) name: String,
    /// The vector's files of the query and the subject.
    pub(crate) query: String,
    pub(crate) subject: String,
    /// The subject's name (the CLI's `-subject` path, the handle's label).
    pub(crate) label: String,
    pub(crate) outfmt: String,
    pub(crate) extra: Vec<String>,
    /// The report (`Ok`), or the error message (`Err`): the CLI's standard output or error.
    pub(crate) expected: Result<Vec<u8>, Vec<u8>>,
}

/// A file of the vector (read at run time: the wasm32-wasip1 test runner opens `/`).
pub(crate) fn read(name: &str) -> Vec<u8> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/web_v1_tblastn")
        .join(name);
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

pub(crate) fn cases() -> Vec<Case> {
    String::from_utf8(read("cases.tsv"))
        .expect("UTF-8 cases")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let [name, query, subject, label, outfmt, extra, kind] =
                line.split('\t').collect::<Vec<_>>()[..]
            else {
                panic!("seven columns: {line:?}");
            };
            // The files are NCBI's bytes (the vector was frozen identical to NCBI's); LOSAT
            // writes them in its own words (`crate::wording::losat_output_bytes`).
            let expected = read(&format!("{name}.{kind}"));
            let expected = crate::wording::losat_output_bytes(&expected).into_owned();
            Case {
                name: name.to_string(),
                query: query.to_string(),
                subject: subject.to_string(),
                label: label.to_string(),
                outfmt: outfmt.to_string(),
                extra: extra.split_whitespace().map(str::to_string).collect(),
                expected: if kind == "out" {
                    Ok(expected)
                } else {
                    Err(expected)
                },
            }
        })
        .collect()
}
