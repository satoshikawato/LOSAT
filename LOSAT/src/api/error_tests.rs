use super::Error;

#[test]
fn the_message_of_an_error_is_its_bytes_and_its_text_a_view_of_them() {
    let bytes = b"Error: [blastdbcmd] Unrecognized format specification: '%\xC3'".to_vec();
    let error = Error::Search(bytes.clone());
    assert_eq!(error.message(), &bytes[..]);
    let text = "Error: [blastdbcmd] Unrecognized format specification: '%\u{FFFD}'";
    assert_eq!(error.to_string_lossy(), text);
    assert_eq!(error.to_string(), text);
    // One U+FFFD per maximal ill-formed subpart: a truncated `€` (E2 82) is one subpart, a
    // lone continuation byte (80) and FF are one each.
    let error = Error::Input(b"a\xE2\x82b\x80\xFF".to_vec());
    assert_eq!(error.to_string(), "a\u{FFFD}b\u{FFFD}\u{FFFD}");
    // A message that is UTF-8 is its own text, borrowed.
    let error = Error::Argv("-query is required: é".into());
    assert_eq!(error.message(), "-query is required: é".as_bytes());
    assert!(matches!(
        error.to_string_lossy(),
        std::borrow::Cow::Borrowed("-query is required: é")
    ));
    assert_eq!(error.to_string(), "-query is required: é");
    assert_eq!(
        Error::Argv(String::from("x").into()),
        Error::Argv(b"x".to_vec())
    );
    // A tool's diagnostics lose their last line feeds only.
    assert_eq!(
        Error::written(b"a\n\xE9\n\n".to_vec()),
        Error::Search(b"a\n\xE9".to_vec())
    );
}

#[test]
fn the_debug_form_of_an_error_is_its_kind_and_its_escaped_message() {
    assert_eq!(
        format!("{:?}", Error::Argv("-query is required".into())),
        r#"Argv("-query is required")"#
    );
    assert_eq!(
        format!("{:?}", Error::Search(b"q\xE9 \"t\"\n".to_vec())),
        r#"Search("q\xE9 \"t\"\n")"#
    );
    assert_eq!(
        format!("{:?}", Error::Input("é\0".into())),
        format!("Input({:?})", "é\0")
    );
    assert_eq!(
        format!("{:#?}", Error::Argv("x".into())),
        "Argv(\n    \"x\",\n)"
    );
}

// An exception of NCBI's that a report's writer raises mid-report (`-html` with a `./.ncbirc` that
// does not parse: the formatter's registry, `report::formatter_registry`) is the command line's
// report of it, `Error: `, the exception's text and its newlines, whichever program wrote the
// report (line RG T4b: the messages of blastn, blastx and tblastx lacked the `Error: ` and the
// last newline). The expected bytes are NCBI's: the `CATCH_ALL` part of the frozen standard error
// of the fixture row `db.html.bp.di.f0.rc_bad.use` (NCBI BLAST+ 2.17.0 in the run directory of that
// `./.ncbirc`), in LOSAT's words, as the fixtures' checker reads it (`wording::losat_bytes`).
#[test]
fn an_exception_of_a_report_s_writer_is_the_command_line_s_report_of_it() {
    use crate::blastinput::ncbi_registry::{read_registry, Layer};
    let frozen = crate::wording::losat_bytes(include_bytes!(
        "../../tests/fixtures/db_search/db.html.bp.di.f0.rc_bad.use.err"
    ));
    // The application's `Critical` post of its own read of the file comes first, then a blank line.
    let at = frozen
        .windows(8)
        .position(|window| window == b"\nError: ")
        .expect("the CATCH_ALL report")
        + 1;
    let expected = &frozen[at..];
    let exception = read_registry(
        include_bytes!("../../tests/db_search/html/rc/bad/.ncbirc"),
        Layer::SystemWide,
    )
    .expect("a read")
    .error
    .expect("the syntax error");
    let error = anyhow::Error::from(exception.into_io_error());
    assert_eq!(
        String::from_utf8_lossy(Error::search(&error).message()),
        String::from_utf8_lossy(expected)
    );
    // The command line's report is the same function's (`cli::NativeError::caught`).
    let native = crate::cli::NativeError::caught(&error).expect("CATCH_ALL's report");
    assert_eq!((native.exit, &native.message[..]), (255, expected));
}
