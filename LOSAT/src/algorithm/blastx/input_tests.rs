use super::*;

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:214-223 (cited at
// the reader): the title is the defline's bytes up to the first byte below 0x20 after its first,
// as they are. BLASTX's record and the shared reader's record of the same title have the same
// names (`fasta_reader::RecordNames`, their one owner).
#[test]
fn titles_keep_their_bytes_and_their_names_are_the_shared_records() {
    let titles: [&[u8]; 13] = [
        b"q1 a\xc3\xa9b z",
        b"q1 a\xe9b z",
        b"q1 a\xa9b z",
        b"q1 a\x81b z",
        b"q1 a\xc0\x80b z",
        b"q1 a\xe0\x81\x81b z",
        b"q1 a\xed\xa0\x80b z",
        b"q1 a\xf5\x80\x80\x80b z",
        b"q1 one\xe9",
        b"\x01q1 ab z",
        b"q1\xe9 rest",
        b"Subject_1\x81 rest",
        b"q1 a title longer than thirty-five bytes \xe9 end",
    ];
    for title in titles {
        let mut data = b">".to_vec();
        data.extend_from_slice(title);
        data.extend_from_slice(b"\nACGT\n");
        for protein in [false, true] {
            let records = parse_fasta(&data, protein, false).unwrap();
            let record = &records[0];
            assert_eq!(record.title, title, "{title:?}");
            let shared = crate::blastinput::fasta_reader::FastaRecord::new(
                record.internal_id.clone(),
                title,
                b"ACGT",
            );
            assert_eq!(record.query_line(), shared.query_line(), "{title:?}");
            assert_eq!(record.shown_id(), shared.shown_id(), "{title:?}");
            assert_eq!(record.warning_id(), shared.warning_id(), "{title:?}");
            assert_eq!(
                record.names().warning_query_id(),
                shared.names().warning_query_id(),
                "{title:?}"
            );
            assert_eq!(
                record.subject_tabular_id(),
                shared.subject_tabular_id(true).unwrap(),
                "{title:?}"
            );
        }
    }
    // The rules themselves, as NCBI writes them: a warning's id and title cut to 25 bytes and
    // `.. ` past 35 (blast_setup_cxx.cpp:533-543), the tabular id the title's first word
    // (tabular.cpp:474-504), a `Subject_` word's id the first word of the title
    // `GenerateDefline` makes, none where it throws (showdefline.cpp:225-231).
    let record = |title: &[u8], protein: bool| {
        let mut data = b">".to_vec();
        data.extend_from_slice(title);
        data.extend_from_slice(b"\nACGT\n");
        parse_fasta(&data, protein, false).unwrap().remove(0)
    };
    let long = b"q1 a title longer than thirty-five bytes \xe9 end";
    assert_eq!(
        record(long, false).names().warning_query_id(),
        b"Query_1 q1 a title longer.. "
    );
    assert_eq!(
        record(long, true).names().warning_query_id(),
        b"Subject_1 q1 a title long.. "
    );
    assert_eq!(
        record(b"q1 one\xe9", false).names().warning_query_id(),
        b"Query_1 q1 one\xe9"
    );
    assert_eq!(&*record(b"q1 a\x81b z", false).query_line(), b"q1 a\x81b z");
    assert_eq!(record(b"q1\xe9 rest", false).shown_id(), b"q1\xe9");
    assert_eq!(record(b"\x01q1 ab z", false).shown_id(), b"\x01q1");
    assert_eq!(
        record(b"q1 a\xe9b z", true).subject_tabular_id(),
        Some(b"q1".to_vec())
    );
    assert_eq!(
        record(b"Subject_1\x81 rest", true).subject_tabular_id(),
        None
    );
    // A byte below 0x20 after the first ends the title.
    for (data, title) in [
        (&b">q1 a\x01b z\nACGT\n"[..], &b"q1 a"[..]),
        (b">q1 a\tb z\nACGT\n", b"q1 a"),
        (b">q1\x01 rest\nACGT\n", b"q1"),
    ] {
        assert_eq!(parse_fasta(data, false, false).unwrap()[0].title, title);
    }
}
#[test]
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:736-751
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:974-990
fn ncbi_reader_preserves_normalized_positions_and_warnings() {
    let r = parse_fasta(b">  q description   \r\natg1-2-3ACT;rest\r\n", false, true).unwrap();
    assert_eq!(r[0].title, b"q description");
    assert_eq!(r[0].internal_id, "Query_1");
    assert_eq!(r[0].sequence, b"ATGACT");
    assert_eq!(r[0].lowercase_masks, vec![(0, 2)]);
    assert!(parse_fasta(b">q\nATG?!ACT", false, false).is_err());
}
#[test]
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:382-390
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:437-442
fn empty_file_record_and_optional_defline_differ() {
    assert!(parse_fasta(b"", false, false).unwrap().is_empty());
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:339-340
    assert!(parse_fasta(b">q\n", false, false).unwrap()[0]
        .sequence
        .is_empty());
    assert_eq!(
        parse_fasta(b"ATGACT", false, false).unwrap()[0].sequence,
        b"ATGACT"
    );
    let r = parse_fasta(b">dup\nATG\n>dup\nACT", false, false).unwrap();
    assert_eq!(
        r.iter().map(|r| r.internal_id.as_str()).collect::<Vec<_>>(),
        vec!["Query_1", "Query_2"]
    );
}
