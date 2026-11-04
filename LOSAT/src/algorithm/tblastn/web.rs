//! Web ABI v1's TBLASTN (`web_api::run_pair`, `run_pair_handles`), added on the Owner's request
//! of 2026-11-04 (gbdraw's serial module searches a protein query against its records); the
//! programs ABI v1 had before stay frozen (plan TD-1). The CLI's parser reads the words, and
//! the CLI's run (`TblastnArgs::run_in_memory`) reads the call's bytes as the files with the
//! same bytes: the report, the diagnostics and the errors are the command line's.
use super::TblastnArgs;
use anyhow::{bail, Result};

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:114-117
// ```c++
// int CTblastnApp::Run(void)
// {
//     int status = BLAST_EXIT_SUCCESS;
// 	const CArgs& args = GetArgs();
// ```
/// The words of a request (`extra_args`, then `-outfmt` and the `outfmt` parameter) as the CLI's
/// parser reads them after `-query query -subject subject`: its defaults, value checks and errors.
/// The inputs and the output are the call's, so `-query`, `-subject`, `-out`, `-db` and the
/// help and version words are refused.
pub(crate) fn parse_args(extra: &[&str]) -> Result<TblastnArgs> {
    // The host has no database files; LOSAT Web searches the records it gives (`run_local`).
    if extra
        .iter()
        .any(|key| key.split('=').next().unwrap_or(key) == "-db")
    {
        bail!("-db (a database search) is not supported in the browser yet");
    }
    if extra.iter().any(|key| {
        matches!(
            key.split('=').next().unwrap_or(key),
            "-query" | "-subject" | "-out" | "-h" | "-help" | "--help" | "-version" | "--version"
        )
    }) {
        bail!("file/help arguments are unsupported for the TBLASTN web API");
    }
    let mut argv = vec!["losat", "tblastn", "-query", "query", "-subject", "subject"];
    argv.extend_from_slice(extra);
    match crate::cli::try_parse_from::<crate::cli::Cli, _, _>(argv) {
        Ok(cli) => match cli.command {
            crate::cli::Commands::Tblastn(args) => Ok(args),
            _ => unreachable!("fixed TBLASTN subcommand"),
        },
        // The CLI's message for the words (`main` writes it to standard error).
        Err(error) => Err(anyhow::Error::msg(crate::cli::render_message(&error))),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:213-219
// ```c++
//             if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
//                	ERR_POST(Warning << "Query is Empty!");
//                	return BLAST_EXIT_SUCCESS;
//             }
//             fasta.Reset(new CBlastFastaInputSource(
//                                          m_CmdLineArgs->GetInputStream(),
//                                          iconfig));
// ```
/// One search of `query` (protein FASTA) against `subject` (nucleotide FASTA) with `args`
/// (`parse_args`), as `losat tblastn` searches the files `query_label` and `subject_label` with
/// the same bytes, without NCBI's registry (`.ncbirc`, `tblastn.ini`, `NCBI_CONFIG__*`); the other
/// variables the CLI reads are read alike. The labels are only the names of the inputs (outfmt 0
/// and 7 show `Input: <subject_label>`; no TBLASTN report shows the query's).
pub(crate) fn run_pair(
    mut args: TblastnArgs,
    query: &[u8],
    subject: &[u8],
    query_label: &str,
    subject_label: &str,
) -> Result<Vec<u8>> {
    args.query = query_label.into();
    args.subject = Some(subject_label.into());
    // The report is returned, never written to a file (`parse_args` refuses `-out`).
    args.out = None;
    args.run_in_memory(Default::default(), query, subject)
}

/// The shared vector of ABI v1's TBLASTN (`tests/fixtures/web_v1_tblastn/cases.tsv`): the bytes
/// `losat tblastn` writes for each case (`tests/web_v1_tblastn.rs` checks them against the CLI),
/// which `run_pair` (the tests below) and the ABI's exports (`web_api.rs`) must return.
#[cfg(test)]
pub(crate) mod vector {
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
                let expected = read(&format!("{name}.{kind}"));
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn web_v1_run_pair_writes_the_cli_bytes_of_the_shared_vector() {
        let cases = vector::cases();
        assert_eq!(cases.len(), 11);
        for case in cases {
            let (query, subject) = (vector::read(&case.query), vector::read(&case.subject));
            let mut words: Vec<&str> = case.extra.iter().map(String::as_str).collect();
            words.extend(["-outfmt", case.outfmt.as_str()]);
            let result = parse_args(&words)
                .and_then(|args| run_pair(args, &query, &subject, "query", &case.label))
                .map_err(|error| format!("{error:#}").into_bytes());
            let shown = |bytes: &Vec<u8>| String::from_utf8_lossy(bytes).into_owned();
            assert!(
                result == case.expected,
                "{}: {:?}",
                case.name,
                result.as_ref().map(shown).map_err(shown)
            );
        }
    }

    #[test]
    fn web_v1_refuses_the_words_of_files_and_databases() {
        let error = |words: &[&str]| format!("{:#}", parse_args(words).unwrap_err());
        assert_eq!(
            error(&["-db", "nt"]),
            "-db (a database search) is not supported in the browser yet"
        );
        assert_eq!(
            error(&["-db=nt"]),
            "-db (a database search) is not supported in the browser yet"
        );
        for word in [
            "-query",
            "-subject",
            "-out",
            "-help",
            "-version",
            "--version",
        ] {
            assert_eq!(
                error(&[word, "x"]),
                "file/help arguments are unsupported for the TBLASTN web API",
                "{word}"
            );
        }
        // Any other word is the CLI parser's, with its message (the vector's `unknown` case).
        assert!(parse_args(&["-evalue", "1e-5", "-outfmt", "6"]).is_ok());
    }
}
