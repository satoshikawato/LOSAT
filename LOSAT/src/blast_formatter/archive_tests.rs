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
    let (status, out, err) = run_archive(&["-outfmt", "5"]);
    assert_eq!(
        (status, out.as_slice(), err),
        (
            1,
            &b""[..],
            format!("Error: {}\n", Refusal::OutputFormat(5).reason())
        )
    );
    // The tabular formats are written.
    let (status, out, err) = run_archive(&["-outfmt", "6"]);
    assert_eq!((status, err.as_str()), (0, ""));
    assert!(
        out.starts_with(b"q1\ts1\t"),
        "{}",
        String::from_utf8_lossy(&out)
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
fn a_failed_writer_of_the_caller_ends_the_run_with_the_write_failure() {
    // A library run's standard output is the caller's writer, no stream of NCBI's (as the
    // tools' and the search programs' library writers): its failed write ends the run with
    // `CATCH_ALL`'s text, wherever the diagnostics go. NCBI's runtime on a failed standard output
    // of the command line (`std::terminate`): `blast_formatter_standard_output_failure_aborts` of
    // tests/masker_fixtures.rs.
    for outfmt in ["8", "9", "11", "12"] {
        assert_eq!(
            run_archive_full(&["-outfmt", outfmt]),
            (6, "BLAST failed to write output\n".to_string()),
            "outfmt {outfmt}"
        );
    }
    assert_eq!(
        run_archive_full(&["-outfmt", "8", "-logfile", "/dev/null"]),
        (6, String::new())
    );
}

/// A file on which every write fails, as on a full device.
struct FullFile;

impl std::io::Write for FullFile {
    fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::from(std::io::ErrorKind::StorageFull))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn the_formatter_s_constructor_throws_on_a_stream_that_failed() {
    let files = crate::report::run_files::OpenFiles::disk();
    let posts = Posts::default();
    let mut output = Output {
        stdout: None,
        files: &files,
        stream: Some(OutStream::Report(ReportStream::formatter_output(
            Some(Box::new(FullFile)),
            Box::new(posts.clone()),
            false,
        ))),
        formatting: false,
        posts,
    };
    let target = crate::blastinput::ncbi_environment::DiagTarget::Stderr;
    let mut stderr = Vec::new();
    let mut diag = Diag::new(&target, &mut stderr, &super::super::TOOL);
    output.begin_object(Site::Report, &diag).unwrap();
    // The flush of what the file buffer holds fails (NCBI's write, made by the object's end at
    // the latest): the formatter's exception, "BLAST failed to write output" at its end.
    let writer = output.writer();
    writer.write_all(b"x").unwrap();
    let _ = writer.ncbi_flush();
    let (status, message) = native(output.end_object(Ok(()), &mut diag).unwrap_err());
    assert_eq!(
        (status, message.as_str()),
        (6, "BLAST failed to write output\n")
    );
    // The next object's formatter throws at its start: the same end.
    let started = output.begin_object(Site::Report, &diag);
    assert!(started.is_err());
    let (status, message) = native(output.end_object(started, &mut diag).unwrap_err());
    assert_eq!(
        (status, message.as_str()),
        (6, "BLAST failed to write output\n")
    );
    drop(diag);
    assert!(stderr.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2801-2851
#[cfg(unix)]
#[test]
fn non_utf8_outfmt_argument_preserves_raw_bytes_in_error() {
    use std::os::unix::ffi::OsStringExt;

    let archive = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/blast_formatter/in/bn.asn"
    ))
    .unwrap();
    for bad in [&b"8\xff"[..], &b"\xff"[..]] {
        let line = vec![
            OsString::from("-archive"),
            OsString::from("-"),
            OsString::from("-outfmt"),
            OsString::from_vec(bad.to_vec()),
        ];
        let mut out = Vec::new();
        let mut err = Vec::new();
        let status = super::super::run_with(&line, &mut &archive[..], &mut out, &mut err);
        let expected = [
            b"LOSAT query/options error: '",
            bad,
            b"' is not a valid output format\n",
            crate::wording::USER_MANUAL_FOOTER.as_bytes(),
            b"\n",
        ]
        .concat();
        assert_eq!((status, out, err), (1, Vec::new(), expected), "{bad:?}");
    }
}

/// `psi_round`, the one check of the round of the formats: a `Web_StepNumber` that is no integer
/// throws at the first result set (NCBI rows of7_psiblast_step_string and of7/of6/of20/of8
/// _ps_none_step_str, a query without hits), and not where an object has no result set.
#[test]
fn a_step_that_is_no_integer_throws_at_the_first_result_set() {
    assert_eq!(psi_round(Ok(Some(2)), 1).unwrap(), Some(2));
    assert_eq!(psi_round(Ok(None), 3).unwrap(), None);
    assert_eq!(psi_round(Err(()), 0).unwrap(), None);
    let error = psi_round(Err(()), 1).unwrap_err();
    assert_eq!(
        format!("{error}"),
        format!("{}", invalid_archive()),
        "the archive's error"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:796-818
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:599-608
#[test]
fn binary_and_xml_archives_format_identical_to_text_and_match_truncation_quirks() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/blast_formatter/in/");
    let text = std::fs::read(format!("{dir}bn.asn")).unwrap();
    let binary = std::fs::read(format!("{dir}bn.asnb")).unwrap();
    let xml = std::fs::read(format!("{dir}bn.xml")).unwrap();
    let run_bytes = |bytes: &[u8], outfmt: &str| {
        let line: Vec<OsString> = ["-archive", "-", "-outfmt", outfmt]
            .into_iter()
            .map(OsString::from)
            .collect();
        let mut out = Vec::new();
        let mut err = Vec::new();
        let status = super::super::run_with(&line, &mut &bytes[..], &mut out, &mut err);
        (status, out, String::from_utf8(err).unwrap())
    };
    for outfmt in ["8", "9", "12", "6", "7", "10", "20", "6 std qseq sseq"] {
        let expected = run_bytes(&text, outfmt);
        assert_eq!(expected.0, 0, "text outfmt {outfmt}");
        assert_eq!(run_bytes(&binary, outfmt), expected, "bin outfmt {outfmt}");
        assert_eq!(run_bytes(&xml, outfmt), expected, "xml outfmt {outfmt}");
    }
    // Binary archive missing only its last byte is skipped silently with rc 0 (F3-6 (3), P3.104).
    assert_eq!(
        run_bytes(&binary[..binary.len() - 1], "8"),
        (0, Vec::new(), String::new())
    );
    // Binary archive missing its last 2 bytes is NCBI's invalid input (rc 1).
    let invalid = format!(
        "LOSAT query/options error: Invalid input format for BLAST Archive.\n{}\n",
        crate::wording::USER_MANUAL_FOOTER
    );
    assert_eq!(
        run_bytes(&binary[..binary.len() - 2], "8"),
        (1, Vec::new(), invalid.clone())
    );
    // XML archive with trailing `<junk/>` formats object 1 then reports invalid input (P3.092).
    let (status, out, err) = run_bytes(&[&xml[..], b"<junk/>"].concat(), "8");
    assert_eq!((status, out, err), (1, run_bytes(&text, "8").1, invalid));
}
