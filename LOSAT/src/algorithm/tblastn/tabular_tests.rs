use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786: every tabular kind goes to `CBlastTabularInfo`; LOSAT's
// default 6 and 7 (the `std` fields with a tab) keep `stage_e_report::write_tabular`.
#[test]
fn only_a_field_list_a_delimiter_and_10_and_20_take_the_shared_writer() {
    for outfmt in [
        "6",
        "7",
        "6 std",
        " 7 std ",
        "6 -qseqid",
        "6 delim=",
        "0 qseqid",
        "5",
    ] {
        assert!(tabular_spec(outfmt).is_none(), "{outfmt}");
    }
    for (outfmt, kind) in [
        ("6 qseqid", TabularKind::Tabular),
        ("6 delim=, std", TabularKind::Tabular),
        ("7 qaccver foo", TabularKind::Commented),
        ("10", TabularKind::Csv),
        ("20 qseqid", TabularKind::CsvWithHeader),
    ] {
        assert_eq!(
            tabular_spec(outfmt).map(|spec| spec.kind),
            Some(kind),
            "{outfmt}"
        );
    }
}

#[test]
fn outfmt_20_writes_its_line_and_its_formatter_s_warning() {
    let spec = tabular_spec("20 qseqid sscinames").unwrap();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    write_tabular_prolog(&spec, false, &mut out, &mut err).unwrap();
    assert_eq!(out, b"qseqid,sscinames\n");
    assert_eq!(err, taxdb_warning_line("tblastn").as_bytes());
    let (mut out, mut err) = (Vec::new(), Vec::new());
    write_tabular_prolog(&spec, true, &mut out, &mut err).unwrap();
    assert_eq!(
        (out.as_slice(), err.as_slice()),
        (&b"qseqid,sscinames\n"[..], &b""[..])
    );
    let (mut out, mut err) = (Vec::new(), Vec::new());
    write_tabular_prolog(
        &tabular_spec("10 sscinames").unwrap(),
        false,
        &mut out,
        &mut err,
    )
    .unwrap();
    assert!(out.is_empty() && err.is_empty());
}
