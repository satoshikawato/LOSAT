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
    let mut out = Vec::new();
    let mut stdout = cout(
        &mut out,
        crate::cli::NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
    );
    let mut stderr = Vec::new();
    let mut diag = Diag::new(&target, &mut stderr, &super::super::TOOL);
    post_asn_posts(
        &mut stdout,
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
    )
    .unwrap();
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
fn post_diag_flushes_cout_only_for_cerr_and_respects_empty_and_suppressed_posts() {
    let block_stream = crate::cli::NcbiStdio::Cout {
        block: 4096,
        lines: false,
    };

    // Empty cout on a closed pipe does not stop a cerr post (matches NCBI D25's first Bad char).
    let mut pipe = FailingWriter(std::io::ErrorKind::BrokenPipe);
    let mut stdout = cout(&mut pipe, block_stream);
    let target = crate::blastinput::ncbi_environment::DiagTarget::Stderr;
    let mut stderr = Vec::new();
    {
        let mut diag = Diag::new(&target, &mut stderr, &super::super::TOOL);
        post_diag(&mut stdout, &mut diag, |d| d.post("first")).unwrap();
    }
    assert_eq!(stderr, b"first\n");

    // Suppressed Utf8Tag warning does not flush buffered cout on a closed pipe.
    stdout.insert(b"buffered\n");
    {
        let mut diag = Diag::new(&target, &mut stderr, &super::super::TOOL);
        post_asn_posts(
            &mut stdout,
            &mut diag,
            &[Post::Utf8Tag(
                "PssmWithParameters.pssm.query.seq.descr.E.title".to_string(),
            )],
        )
        .unwrap();
    }
    assert!(!stdout.rdbuf().closed_pipe());

    // Non-empty cout on a closed pipe stops before writing to cerr (off Windows).
    let res = {
        let mut diag = Diag::new(&target, &mut stderr, &super::super::TOOL);
        post_diag(&mut stdout, &mut diag, |d| d.post("second"))
    };
    if cfg!(not(windows)) {
        assert!(matches!(res, Err(RunError::ClosedPipe)));
        assert_eq!(stderr, b"first\n");
    } else {
        assert!(res.is_ok());
        assert_eq!(stderr, b"first\nsecond\n");
    }

    // Non-cerr target (e.g. -logfile) is not tied to cout and does not flush buffered stdout.
    let mut pipe2 = FailingWriter(std::io::ErrorKind::BrokenPipe);
    let mut stdout2 = cout(&mut pipe2, block_stream);
    stdout2.insert(b"still buffered\n");
    let mem_target = crate::blastinput::ncbi_environment::DiagTarget::Memory(Default::default());
    let mut stderr2 = Vec::new();
    let mut diag2 = Diag::new(&mem_target, &mut stderr2, &super::super::TOOL);
    post_diag(&mut stdout2, &mut diag2, |d| d.post("logged")).unwrap();
    assert!(!stdout2.rdbuf().closed_pipe());
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

struct FailingWriter(std::io::ErrorKind);

impl Write for FailingWriter {
    fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::from(self.0))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Err(std::io::Error::from(self.0))
    }
}

#[test]
fn cout_log_buffers_exclusion_lines_until_flush_or_overflow() {
    let mut pipe = FailingWriter(std::io::ErrorKind::BrokenPipe);
    let mut stream = cout(
        &mut pipe,
        crate::cli::NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
    );
    // An insertion within the 4096-byte pipe buffer is held and does not observe EPIPE yet.
    let line = [b'x'; 70];
    for _ in 0..58 {
        stream.insert(&line);
        assert!(!stream.rdbuf().closed_pipe());
    }
    // The 59th line (4130 > 4096) flushes a block and observes BrokenPipe off Windows.
    stream.insert(&line);
    assert_eq!(stream.rdbuf().closed_pipe(), cfg!(not(windows)));

    let mut full = FailingWriter(std::io::ErrorKind::WriteZero);
    let mut full_stream = cout(
        &mut full,
        crate::cli::NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
    );
    full_stream.insert(b"Deleted");
    full_stream.endl();
    assert!(!full_stream.rdbuf().closed_pipe());
}
