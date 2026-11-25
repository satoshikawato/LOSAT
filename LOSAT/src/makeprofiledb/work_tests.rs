use super::*;
use crate::report::serial::asn_types::Post;
use crate::report::serial::encoder::{NullElement, ReplacedChar, Stream, StreamPosition};

#[test]
fn read_error_to_run_error_maps_all_variants() {
    let default_rv = ReadVerification::default();
    assert!(matches!(
        read_error_to_run_error(
            ReadError::Eof,
            "PssmWithParameters",
            "a profile (-in)",
            &default_rv,
            |w| w,
        ),
        RunError::Exception(msg) if msg == b" Frame type= eFrameNamed, Object type= PssmWithParameters"
    ));
    assert!(matches!(
        read_error_to_run_error(
            ReadError::Eof,
            "Blast-def-line-set",
            "binary definition lines (ASN1_BlastDefLine)",
            &default_rv,
            |w| w,
        ),
        RunError::Exception(msg) if msg == b" Frame type= eFrameNamed, Object type= Blast-def-line-set"
    ));
    assert!(matches!(
        read_error_to_run_error(
            ReadError::Std("basic_string::_M_create"),
            "PssmWithParameters",
            "a profile (-in)",
            &default_rv,
            |w| w,
        ),
        RunError::Unknown
    ));
    assert!(matches!(
        read_error_to_run_error(
            ReadError::OutOfMemory,
            "Blast-def-line-set",
            "binary definition lines (ASN1_BlastDefLine)",
            &default_rv,
            |w| w,
        ),
        RunError::Unknown
    ));
    assert!(matches!(
        read_error_to_run_error(
            ReadError::Unsupported("too big".to_string()),
            "PssmWithParameters",
            "a profile (-in)",
            &default_rv,
            |w| format!("ctx: {w}"),
        ),
        RunError::Refused(msg) if msg == "ctx: too big"
    ));
    let read_off = ReadVerification::of(Some(std::ffi::OsStr::new("NO")), None);
    let expected = crate::ncbi_exception::NcbiException::new(
        crate::report::serial::object_frame::EXPECTED_MEMBER,
        "member seqid expected",
    );
    assert!(matches!(
        read_error_to_run_error(
            ReadError::Serial(expected),
            "Blast-def-line-set",
            "binary definition lines (ASN1_BlastDefLine)",
            &read_off,
            |w| w,
        ),
        RunError::Refused(msg) if msg.contains("SERIAL_VERIFY_DATA_READ=NO")
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

#[test]
fn aborts_strips_trailing_newline_from_exception_message() {
    let assert_refused =
        |first: RunError, then: &WriteDbError, want: &str| match aborts(first, then) {
            RunError::Refused(msg) => {
                assert!(!msg.contains('\n'), "refusal split across lines: {msg:?}");
                assert_eq!(msg, want);
            }
            other => panic!("expected Refused, got {other:?}"),
        };

    let null = NullElement {
        stream: Stream::AsnBinary,
        element_frame: true,
        trace: "Blast-def-line-set.E.seqid.E".to_string(),
        position: StreamPosition::Byte(15),
        frames: Vec::new(),
    };
    let then = WriteDbError::Exception(null.exception());
    assert!(then.message().ends_with(b"\n"));
    assert_refused(
        RunError::Input(b"in/nope.smp does not exists".to_vec()),
        &then,
        "an error (INPUT ERROR: in/nope.smp does not exists) after which the database's \
         writer stops with another (Error: byte 15: cannot write unassigned member NULL \
         element while writing container ; NCBI's makeprofiledb aborts there)",
    );

    let plain = WriteDbError::Other("No tax info for any oid".to_string());
    assert!(!plain.message().ends_with(b"\n"));
    assert_refused(
        RunError::Exception(b"boom".to_vec()),
        &plain,
        "an error (ERROR: boom) after which the database's writer stops with another \
         (Error: No tax info for any oid; NCBI's makeprofiledb aborts there)",
    );
}
