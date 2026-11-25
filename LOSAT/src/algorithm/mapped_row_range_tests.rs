//! Shared `-subject` writer parity table for BLASTX and TBLASTN when a gapped edit script still
//! starts with a gap after `-subject`'s Seq-loc mapper drops at most one end segment (`[Del, Ins]`
//! -> `[Ins]`, `[Ins, Del]` -> `[Del]`).
//!
//! NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1568-1580
//! NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-191
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:243-245
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:816-824

use crate::common::GapEditOp::{self, Del, Ins, Sub};

const ROW_SEQIDS_MSG: &str = "CSeq_align::GetSeqRange(): Row seqids not consistent.";

enum Expected {
    /// Report succeeds (exit 0) and writes the HSP / row with no exception.
    WritesHit,
    /// XML / XML2 catches `Row seqids not consistent`, posts it to diagnostics, leaves `<Hit>`
    /// out of the document, and returns exit 0.
    XmlCaughtException,
    /// Tabular outfmt 7 with `qcovs` writes its comment header lines, then throws
    /// `Row seqids not consistent` (exit 255) with no data row.
    TabularThrowsAfterHeader,
}

#[derive(Clone, Copy, Debug)]
enum Format {
    Xml5,
    Xml16,
    Tabular7Qcovs,
}

#[test]
fn blastx_and_tblastn_subject_writers_handle_mapped_leading_gaps_identically() {
    let cases: &[(&str, Vec<GapEditOp>, Format, Expected)] = &[
        // `[Del(2), Ins(1), Sub(10)]` -> mapped `[Ins(1), Sub(10)]`: subject starts with gap.
        // XML (5 & 16) calls GetSeqStart(0) then GetSeqStart(1) -> caught exception, exit 0.
        // `qcovs` reads GetSeqRange(0) only -> query row has letters, so row is written.
        (
            "del_ins_xml5",
            vec![Del(2), Ins(1), Sub(10)],
            Format::Xml5,
            Expected::XmlCaughtException,
        ),
        (
            "del_ins_xml16",
            vec![Del(2), Ins(1), Sub(10)],
            Format::Xml16,
            Expected::XmlCaughtException,
        ),
        (
            "del_ins_7_qcovs",
            vec![Del(2), Ins(1), Sub(10)],
            Format::Tabular7Qcovs,
            Expected::WritesHit,
        ),
        // `[Ins(1), Del(2), Sub(10)]` -> mapped `[Del(2), Sub(10)]`: query starts with gap.
        // XML (5 & 16) catches exception, exit 0; outfmt 7 `qcovs` throws after its header lines.
        (
            "ins_del_xml5",
            vec![Ins(1), Del(2), Sub(10)],
            Format::Xml5,
            Expected::XmlCaughtException,
        ),
        (
            "ins_del_xml16",
            vec![Ins(1), Del(2), Sub(10)],
            Format::Xml16,
            Expected::XmlCaughtException,
        ),
        (
            "ins_del_7_qcovs",
            vec![Ins(1), Del(2), Sub(10)],
            Format::Tabular7Qcovs,
            Expected::TabularThrowsAfterHeader,
        ),
        // Controls: single leading gap `[Del(2), Sub(10)]` maps cleanly to `[Sub(10)]`.
        (
            "del_sub_xml5",
            vec![Del(2), Sub(10)],
            Format::Xml5,
            Expected::WritesHit,
        ),
        (
            "del_sub_xml16",
            vec![Del(2), Sub(10)],
            Format::Xml16,
            Expected::WritesHit,
        ),
        (
            "del_sub_7_qcovs",
            vec![Del(2), Sub(10)],
            Format::Tabular7Qcovs,
            Expected::WritesHit,
        ),
        ("sub_xml5", vec![Sub(10)], Format::Xml5, Expected::WritesHit),
        (
            "sub_xml16",
            vec![Sub(10)],
            Format::Xml16,
            Expected::WritesHit,
        ),
        (
            "sub_7_qcovs",
            vec![Sub(10)],
            Format::Tabular7Qcovs,
            Expected::WritesHit,
        ),
    ];

    for (name, ops, format, expected) in cases {
        let bx = {
            use super::blastx::mapped_row_range_tests::{run_subject_report, SubjectReportKind};
            let kind = match format {
                Format::Xml5 => SubjectReportKind::Xml5,
                Format::Xml16 => SubjectReportKind::Xml16,
                Format::Tabular7Qcovs => SubjectReportKind::Tabular7Qcovs,
            };
            run_subject_report(ops.clone(), kind)
        };
        let tbn = {
            use super::tblastn::mapped_row_range_tests::{run_subject_report, SubjectReportKind};
            let kind = match format {
                Format::Xml5 => SubjectReportKind::Xml5,
                Format::Xml16 => SubjectReportKind::Xml16,
                Format::Tabular7Qcovs => SubjectReportKind::Tabular7Qcovs,
            };
            run_subject_report(ops.clone(), kind)
        };

        for (prog, header_tag, stdout, stderr, error) in [
            ("blastx", "# BLASTX ", &bx.stdout, &bx.stderr, &bx.error),
            (
                "tblastn",
                "# TBLASTN ",
                &tbn.stdout,
                &tbn.stderr,
                &tbn.error,
            ),
        ] {
            let ctx = format!("{prog} {name}");
            match expected {
                Expected::WritesHit => {
                    assert!(error.is_none(), "{ctx}: unexpected error {error:?}");
                    assert!(stderr.is_empty(), "{ctx}: unexpected stderr {stderr:?}");
                    match format {
                        Format::Xml5 | Format::Xml16 => {
                            assert!(stdout.contains("<Hit>"), "{ctx}: missing <Hit> in {stdout}");
                        }
                        Format::Tabular7Qcovs => {
                            assert!(
                                stdout.contains(header_tag) && stdout.contains("# 1 hits found\n"),
                                "{ctx}: missing outfmt 7 header in {stdout}"
                            );
                            let data_rows = stdout.lines().filter(|l| !l.starts_with('#')).count();
                            assert_eq!(data_rows, 1, "{ctx}: expected 1 data row in {stdout}");
                        }
                    }
                }
                Expected::XmlCaughtException => {
                    assert!(
                        error.is_none(),
                        "{ctx}: XML must return exit 0, got {error:?}"
                    );
                    assert!(
                        stderr.contains(ROW_SEQIDS_MSG),
                        "{ctx}: missing caught exception in stderr {stderr:?}"
                    );
                    assert!(
                        !stdout.contains("<Hit>"),
                        "{ctx}: failing HSP must not write <Hit> in {stdout}"
                    );
                    match format {
                        Format::Xml5 => {
                            assert!(
                                stdout.contains("<BlastOutput>")
                                    && stdout.contains("</BlastOutput>"),
                                "{ctx}: incomplete XML5 document {stdout}"
                            );
                        }
                        Format::Xml16 => {
                            assert!(
                                stdout.contains("<BlastXML2") && stdout.contains("</BlastXML2>"),
                                "{ctx}: incomplete XML2 document {stdout}"
                            );
                        }
                        Format::Tabular7Qcovs => unreachable!(),
                    }
                }
                Expected::TabularThrowsAfterHeader => {
                    let native = error
                        .as_ref()
                        .unwrap_or_else(|| panic!("{ctx}: expected NativeError"));
                    assert_eq!(native.exit, 255, "{ctx}: expected exit 255");
                    let msg = String::from_utf8_lossy(&native.message);
                    assert!(
                        msg.contains(ROW_SEQIDS_MSG),
                        "{ctx}: unexpected error message {msg:?}"
                    );
                    assert!(
                        stdout.contains(header_tag)
                            && stdout.contains("# Query: q1\n")
                            && stdout.contains("# 1 hits found\n"),
                        "{ctx}: outfmt 7 header lines must precede exception: {stdout:?}"
                    );
                    let data_rows = stdout.lines().filter(|l| !l.starts_with('#')).count();
                    assert_eq!(
                        data_rows, 0,
                        "{ctx}: no data row may be written after qcovs throws: {stdout:?}"
                    );
                }
            }
        }
    }
}
