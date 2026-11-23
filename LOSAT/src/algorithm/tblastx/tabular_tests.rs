use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357: outfmt 20's line of the field names, after the warning of its
// formatter when `taxdb.bti` is not found and a plural name is shown.
#[test]
fn outfmt_20_writes_its_line_and_its_formatter_s_warning() {
    let spec = |outfmt: &str| match crate::algorithm::tblastx::report::output_format(outfmt) {
        crate::algorithm::tblastx::report::TblastxOutputFormat::CustomTabular(spec) => spec,
        other => panic!("{outfmt}: {other:?}"),
    };
    let (mut out, mut err) = (Vec::new(), Vec::new());
    write_tabular_prolog(
        spec("20 qseqid sscinames"),
        false,
        "tblastx",
        &mut out,
        &mut err,
    )
    .unwrap();
    assert_eq!(out, b"qseqid,sscinames\n");
    assert_eq!(err, taxdb_warning_line("tblastx").as_bytes());
    let (mut out, mut err) = (Vec::new(), Vec::new());
    write_tabular_prolog(
        spec("20 delim=| qseqid"),
        true,
        "tblastx",
        &mut out,
        &mut err,
    )
    .unwrap();
    assert_eq!(
        (out.as_slice(), err.as_slice()),
        (&b"qseqid\n"[..], &b""[..])
    );
    let (mut out, mut err) = (Vec::new(), Vec::new());
    write_tabular_prolog(spec("10 sscinames"), false, "tblastx", &mut out, &mut err).unwrap();
    assert!(out.is_empty() && err.is_empty());
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908 and 970-986: `nident` alone enters the block (and reads the rows)
// with a positive engine count only; the other count fields always, the row fields read the
// rows.
#[test]
fn the_counts_block_and_the_row_fields_are_ncbis() {
    for field in ROW_FIELDS {
        assert!(
            field == TabularField::NumIdentical || COUNTS_BLOCK.contains(&field),
            "{field:?}"
        );
    }
    assert!(!COUNTS_BLOCK.contains(&TabularField::NumIdentical));
    assert!(!COUNTS_BLOCK.contains(&TabularField::Frames));
    assert!(!COUNTS_BLOCK.contains(&TabularField::SubjectFrame));
}
