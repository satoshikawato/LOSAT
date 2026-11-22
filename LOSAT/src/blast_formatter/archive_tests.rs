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

#[test]
fn bl2seq_legacy_is_refused_for_subject_archives_as_by_the_search_programs() {
    // The variable is read through `check_search_environment`; its refusal is the search
    // programs' text with blast_formatter's name.
    let error = crate::blastinput::app::check_search_environment(PROGRAM, true);
    if std::env::var_os("BL2SEQ_LEGACY").is_some() {
        assert!(error.unwrap_err().to_string().contains("BLAST_FORMATTER"));
    } else {
        assert!(error.is_ok());
    }
}
