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
    assert_eq!(
        error(&["-in_pssm", "p.asn"]),
        "-in_pssm (a PSSM query) is not supported in the browser yet"
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
