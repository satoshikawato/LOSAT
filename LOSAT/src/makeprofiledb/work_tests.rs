use super::*;
use crate::report::serial::asn_types::Post;
use crate::report::serial::encoder::{ReplacedChar, StreamPosition};

#[test]
fn read_error_to_run_error_maps_all_variants() {
    assert!(matches!(
        read_error_to_run_error(ReadError::Eof, "PssmWithParameters", |w| w),
        RunError::Exception(msg) if msg == b" Frame type= eFrameNamed, Object type= PssmWithParameters"
    ));
    assert!(matches!(
        read_error_to_run_error(ReadError::Eof, "Blast-def-line-set", |w| w),
        RunError::Exception(msg) if msg == b" Frame type= eFrameNamed, Object type= Blast-def-line-set"
    ));
    assert!(matches!(
        read_error_to_run_error(
            ReadError::Std("basic_string::_M_create"),
            "PssmWithParameters",
            |w| w
        ),
        RunError::Unknown
    ));
    assert!(matches!(
        read_error_to_run_error(ReadError::OutOfMemory, "Blast-def-line-set", |w| w),
        RunError::Unknown
    ));
    assert!(matches!(
        read_error_to_run_error(ReadError::Unsupported("too big".to_string()), "PssmWithParameters", |w| format!("ctx: {w}")),
        RunError::Refused(msg) if msg == "ctx: too big"
    ));
}

#[test]
fn post_asn_posts_suppresses_utf8_warning_and_writes_errors() {
    let target = crate::blastinput::ncbi_environment::DiagTarget::Stderr;
    let mut stderr = Vec::new();
    let mut diag = Diag::new(&target, &mut stderr, &super::super::TOOL);
    post_asn_posts(
        &mut diag,
        &[
            Post::Utf8Tag("PssmWithParameters.pssm.query.seq.descr.E.title".to_string()),
            Post::UnendedByteBlock,
            Post::BadChar(ReplacedChar {
                trace: "Blast-def-line-set.E.title".to_string(),
                byte: 0x01,
                position: StreamPosition::Byte(20),
                string: b"a\x01b".as_slice().into(),
            }),
        ],
    );
    let text = String::from_utf8(stderr).unwrap();
    assert!(
        !text.contains("803.10"),
        "UTF8 warning must be suppressed: {text}"
    );
    assert!(
        text.contains("Error: (803.4) unended byte block\n"),
        "{text}"
    );
    assert!(
        text.contains("Error: (803.7) Blast-def-line-set.E.title\nBad char [0x1] in string at byte 20\na\x01b\n"),
        "{text}"
    );
}
