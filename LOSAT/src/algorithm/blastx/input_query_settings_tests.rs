use super::*;

/// The directory of this test's files (removed at its end).
fn dir() -> std::path::PathBuf {
    crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastx-query-settings-{}",
        crate::utils::test_dir::process_id()
    ))
}

/// A file of this test with the bytes `text`, open for reading.
fn file(tag: &str, text: &[u8]) -> std::fs::File {
    std::fs::create_dir_all(dir()).unwrap();
    let path = dir().join(tag);
    std::fs::write(&path, text).unwrap();
    std::fs::File::open(&path).unwrap()
}

/// The records of `stream`, of the query (`protein` false) or the subject role.
fn read(
    stream: &mut FastaStream<crate::blastinput::input_files::InputFile>,
    protein: bool,
) -> Result<Vec<FastaRecord>> {
    parse_fasta_stream_keeping(
        stream,
        protein,
        false,
        100,
        false,
        None,
        &mut |_| Ok(()),
        &mut |_| Ok(()),
    )
}

// The read before the search reads as the search does (under `-parse_deflines` a first line
// that is not a defline is the reader's error) and adds its ids to a copy of the search's: a
// subject with the id of a query that only that read has read is not refused, and that read
// stops at a query with the id of a subject the search has read, as the search does.
#[test]
fn a_read_before_the_search_reads_as_the_search_and_keeps_the_search_ids() {
    let settings = QueryStreamSettings {
        program: crate::algorithm::blastx::args::ProgramName::BLASTX,
        parse_deflines: true,
        id_scope: Some(Default::default()),
        strand: crate::report::archive::bioseq::SeqLocStrand::Both,
        data_loaders: false,
    };
    let error = read(
        &mut settings.stream_before_search(file("nodef", b"ACGT\n>q1 x\nACGT\n")),
        false,
    )
    .unwrap_err();
    assert!(error.is::<NoDeflineError>(), "{error}");
    let queries = read(
        &mut settings.stream_before_search(file("query", b">q1 x\nACGTACGT\n")),
        false,
    )
    .unwrap();
    assert_eq!(queries.len(), 1);
    read(
        &mut settings.stream(file("subject", b">q1 y\nMKVL\n")),
        true,
    )
    .unwrap();
    let error = read(
        &mut settings.stream_before_search(file("again", b">q1 x\nACGTACGT\n")),
        false,
    )
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("a query and a subject record have the Seq-id lcl|q1"),
        "{error}"
    );
    std::fs::remove_dir_all(dir()).unwrap();
}
