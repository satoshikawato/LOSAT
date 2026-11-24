//! blast_formatter's object loop: the errors of NCBI's `CATCH_ALL` by the class of the exception,
//! the cut of `PruneSeqalign`, and runs over the fixture archive `tests/blast_formatter/in/bn.asn`
//! (a `blastn -subject` archive of NCBI BLAST+ 2.17.0; the frozen rows of the formats are in
//! `LOSAT/tests/blast_formatter/`).

use std::ffi::OsString;

use super::*;
use crate::ncbi_exception::ThrowSite;

/// An exception of NCBI's of the class and code `kind` with the message `msg`.
fn exception(kind: &'static str, msg: &str) -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: "src/x.cpp",
            line: 1,
            kind,
            function: None,
        },
        msg,
    )
}

/// The exit status and text of an error as `catch_all` ends the run with it.
fn native(error: anyhow::Error) -> (i32, String) {
    let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
    (
        native.exit,
        String::from_utf8(native.message.clone()).unwrap(),
    )
}

#[test]
fn exceptions_end_the_run_as_catch_all_s_branches() {
    let invalid = (
        1,
        format!(
            "LOSAT query/options error: Invalid input format for BLAST Archive.\n{}\n",
            crate::wording::USER_MANUAL_FOOTER
        ),
    );
    assert_eq!(
        native(exception_error(&exception(
            "CSerialException::eFormatError",
            "x"
        ))),
        invalid
    );
    assert_eq!(
        native(exception_error(&exception("CUnassignedMember::eGet", "x"))),
        invalid
    );
    assert_eq!(
        native(stop_error(Stop::Unassigned {
            type_name: "Seq-inst",
            member: "length",
        })),
        invalid
    );
    assert_eq!(
        native(exception_error(&exception(
            "CRemoteBlastException::eServiceNotAvailable",
            "no"
        ))),
        (3, "LOSAT engine error: no\n".to_string())
    );
    assert_eq!(
        native(exception_error(&exception(
            "CBlastException::eInvalidOptions",
            "bad"
        ))),
        (1, "LOSAT options error: bad\n".to_string())
    );
    assert_eq!(
        native(exception_error(&exception(
            "CBlastException::eCoreBlastError",
            "Out of memory"
        ))),
        (4, "BLAST ran out of memory: Out of memory\n".to_string())
    );
    assert_eq!(
        native(exception_error(&exception(
            "CSeqDBException::eFileErr",
            "gone"
        ))),
        (2, "LOSAT database error: gone\n".to_string())
    );
    assert_eq!(
        native(exception_error(&exception(
            "CSeqalignException::eUnsupported",
            "u"
        )))
        .0,
        255
    );
    // What NCBI does not decide is refused.
    let refused = stop_error(Stop::Undefined("front() of an empty list"));
    assert!(refused.downcast_ref::<Refusal>().is_some());
    // An element stated as `NULL` in a list: NCBI's null-pointer report, with its stack's
    // addresses (`read::Stop::NullElement`).
    assert_eq!(
        stop_error(Stop::NullElement).downcast_ref::<Refusal>(),
        Some(&Refusal::Undecided("an element stated as NULL in a list"))
    );
}

/// The run of `losat blast_formatter` with `words` over the fixture archive on the standard input.
fn run_archive(words: &[&str]) -> (i32, Vec<u8>, String) {
    let archive = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/blast_formatter/in/bn.asn"
    ))
    .unwrap();
    let mut line: Vec<OsString> = vec!["-archive".into(), "-".into()];
    line.extend(words.iter().map(OsString::from));
    let mut out = Vec::new();
    let mut err = Vec::new();
    let status = super::super::run_with(&line, &mut &archive[..], &mut out, &mut err);
    (status, out, String::from_utf8(err).unwrap())
}

#[test]
fn the_seq_annot_formats_write_the_archive_s_seq_aligns() {
    let (status, out, err) = run_archive(&["-outfmt", "8"]);
    assert_eq!((status, err.as_str()), (0, ""));
    let text = String::from_utf8(out).unwrap();
    assert!(text.starts_with("Seq-annot ::= {\n  desc {\n"), "{text}");
    // A `-subject` archive: no database name.
    assert!(text.contains("str \"n/a\""), "{text}");
    assert!(text.contains("local str \"Subject_1\""), "{text}");
    let (status, out, err) = run_archive(&["-outfmt", "12"]);
    assert_eq!((status, err.as_str()), (0, ""));
    assert!(out.starts_with(b"{\n"), "{}", String::from_utf8_lossy(&out));
    let (status, out, err) = run_archive(&["-outfmt", "9"]);
    assert_eq!((status, err.as_str()), (0, ""));
    // ASN.1 binary: a SEQUENCE (0x30) of indefinite length.
    assert_eq!(&out[..2], &[0x30, 0x80]);
}

#[test]
fn the_other_formats_are_refused_after_ncbi_s_checks() {
    let (status, out, err) = run_archive(&["-outfmt", "6"]);
    assert_eq!(
        (status, out.as_slice(), err),
        (
            1,
            &b""[..],
            format!("Error: {}\n", Refusal::OutputFormat(6).reason())
        )
    );
    // NCBI's check of SAM comes first; this archive is blastn's, so SAM is refused too.
    let (status, _, err) = run_archive(&["-outfmt", "17"]);
    assert_eq!(
        (status, err),
        (
            1,
            format!("Error: {}\n", Refusal::OutputFormat(17).reason())
        )
    );
    // 13 without `-out`: NCBI's error before the refusal.
    let (status, _, err) = run_archive(&["-outfmt", "13"]);
    assert_eq!(
        (status, err),
        (
            1,
            format!(
                "LOSAT query/options error: Please provide a file name for outfmt 13.\n{}\n",
                crate::wording::USER_MANUAL_FOOTER
            )
        )
    );
}

/// A standard output on which every write fails, as on a full device.
struct FullDevice;

impl std::io::Write for FullDevice {
    fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::from(std::io::ErrorKind::StorageFull))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// `run_archive` with its standard output on `FullDevice`.
fn run_archive_full(words: &[&str]) -> (i32, String) {
    let archive = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/blast_formatter/in/bn.asn"
    ))
    .unwrap();
    let mut line: Vec<OsString> = vec!["-archive".into(), "-".into()];
    line.extend(words.iter().map(OsString::from));
    let mut err = Vec::new();
    let status = super::super::run_with(&line, &mut &archive[..], &mut FullDevice, &mut err);
    (status, String::from_utf8(err).unwrap())
}

#[test]
fn a_failed_standard_output_ends_the_run_as_ncbi_s_runtime_does() {
    // The diagnostics on `cerr`: the close's post flushes `cout` through the tie, which throws
    // inside the destructor (`std::terminate`, exit status 134 in a shell; the library returns
    // -6 and the command line aborts).
    for outfmt in ["8", "9", "12"] {
        assert_eq!(
            run_archive_full(&["-outfmt", outfmt]),
            (
                -6,
                "terminate called after throwing an instance of 'std::__ios_failure'\n  what():  \
                 basic_ios::clear: iostream error\n"
                    .to_string()
            ),
            "outfmt {outfmt}"
        );
    }
    // Without the tie (`-logfile`): the one Seq-annot fits in the C library's buffer, its close
    // fails and posts where the diagnostics go; nothing writes to the stream after it.
    assert_eq!(
        run_archive_full(&["-outfmt", "8", "-logfile", "/dev/null"]),
        (0, String::new())
    );
}

#[test]
fn the_formatter_s_constructor_throws_on_a_stream_that_failed() {
    let files = crate::report::run_files::OpenFiles::disk();
    let mut output = Output {
        stdout: Some(Box::leak(Box::new(FullDevice))),
        files: &files,
        stream: None,
        formatting: false,
    };
    output.open(OsStr::new("-")).unwrap();
    output.begin_object().unwrap();
    let (stream, is_stdout) = output.opened();
    assert!(is_stdout);
    stream.write(b"x");
    stream.flush();
    assert!(stream.take_thrown());
    output.end_object();
    // The next object's formatter: `BLAST failed to write output`, exit status 6.
    let (status, message) = native(output.begin_object().unwrap_err());
    assert_eq!(
        (status, message.as_str()),
        (6, "BLAST failed to write output\n")
    );
}
