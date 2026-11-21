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
