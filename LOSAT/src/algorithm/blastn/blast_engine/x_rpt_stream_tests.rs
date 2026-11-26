use super::*;
use crate::api::local_blast::FormatOutput;
use crate::blastinput::seq_range::RecordPlacement;

// EXPERIMENT (LOSAT_X_RPTSTREAM): a batch's part of the query input numbers its queries from
// 0, as the batch's lists do; without `-query_loc` it has no placements either.
#[test]
fn a_batch_input_is_the_batch_part_numbered_from_zero() {
    let placements = (0..5)
        .map(|index| RecordPlacement {
            offset: 10 * index,
            length: 100 + index,
        })
        .collect();
    let input = QueryInput {
        placements: Placements(placements),
        ordinals: vec![0, 1, 3, 4, 6],
        input_lengths: vec![100; 7],
        skipped: vec![false, false, true, false, false, true, false],
        range_length: Some(50),
    };
    let batch = batch_input(&input, 2..4);
    assert!(batch.placements.ranged());
    assert_eq!(
        (batch.placements.offset(0), batch.placements.length(1, 7)),
        (20, 103)
    );
    assert_eq!(batch.placements.0.len(), 2);
    assert_eq!(batch.ordinals, [3, 4]);
    assert_eq!(batch.range_length, Some(50));
    let whole = QueryInput {
        placements: Placements::default(),
        ordinals: vec![0, 1, 2],
        input_lengths: vec![10; 3],
        skipped: vec![false; 3],
        range_length: None,
    };
    let batch = batch_input(&whole, 1..3);
    assert!(!batch.placements.ranged());
    assert_eq!(batch.placements.offset(1), 0);
    assert_eq!(batch.ordinals, [1, 2]);
}

/// The arguments of a `blastn` command line (`extra` after `-query q -subject s`).
fn args(extra: &[&str]) -> BlastnArgs {
    let mut argv = vec!["losat", "blastn", "-query", "q", "-subject", "s"];
    argv.extend_from_slice(extra);
    crate::cli::try_parse_from(argv).unwrap().blastn()
}

fn record(
    index: usize,
    title: &[u8],
    letters: usize,
) -> crate::blastinput::fasta_reader::FastaRecord {
    crate::blastinput::fasta_reader::FastaRecord {
        local_id: format!("Subject_{}", index + 1),
        title: title.to_vec(),
        sequence: b"ACGTTGCA".iter().copied().cycle().take(letters).collect(),
        warnings: Vec::new(),
        source: Default::default(),
    }
}

/// `check` for queries of 1,000 letters each (`queries`), first batches of `batch_size`
/// letters, the subjects `subjects` and the outputs `outputs`.
fn checked(
    on: bool,
    extra: &[&str],
    subjects: &[crate::blastinput::fasta_reader::FastaRecord],
    queries: usize,
    batch_size: u32,
    split: crate::algorithm::blastn::query_split::SplitSizes,
    outputs: &mut ReportOutputs<'_>,
) -> std::result::Result<Ready, Off> {
    checked_as(
        &[BlastnOutputFormat::Tabular],
        on,
        extra,
        subjects,
        queries,
        batch_size,
        split,
        outputs,
    )
}

/// `checked` for the output formats `formats`.
#[allow(clippy::too_many_arguments)]
fn checked_as(
    formats: &[BlastnOutputFormat],
    on: bool,
    extra: &[&str],
    subjects: &[crate::blastinput::fasta_reader::FastaRecord],
    queries: usize,
    batch_size: u32,
    split: crate::algorithm::blastn::query_split::SplitSizes,
    outputs: &mut ReportOutputs<'_>,
) -> std::result::Result<Ready, Off> {
    let args = args(extra);
    let query_records: Vec<_> = (0..queries)
        .map(|index| record(index, b"q", 1000))
        .collect();
    let query_input = QueryInput::whole(&query_records);
    let lengths: Vec<usize> = query_records.iter().map(|q| q.sequence.len()).collect();
    let prepared = super::super::PreparedSubjects {
        records: subjects,
        stores: None,
        metadata: super::super::subject_metadata_from_records(subjects, |index| {
            subjects[index].sequence.len()
        }),
        packed: Vec::new(),
    };
    let batching = QueryBatching {
        fixed: Some(batch_size),
        batch_size,
        split,
    };
    check(
        on,
        &args,
        &*outputs,
        formats,
        SubjectSource::Prepared(&prepared),
        SubjectRanges {
            placements: &Placements::default(),
            records: None,
        },
        &query_input,
        &lengths,
        batching,
    )
}

// EXPERIMENT (LOSAT_X_RPTSTREAM): where the switch acts, and why not elsewhere: one format of
// the command line's kind only; with one batch no refusal is checked (and nothing is made
// ahead); with more, the refusals of the options and the query lengths that could come after
// the first batch are, and a subject set's refused ids are left to the batches' reports.
#[test]
fn the_switch_acts_where_the_outputs_and_options_allow() {
    use crate::algorithm::blastn::query_split::SplitSizes;
    let plain = SplitSizes::new(false, None, None);
    let good = vec![
        record(0, b"s1 a subject", 4000),
        record(1, b"s2 another", 3000),
    ];
    let refused = vec![
        record(0, b"s1 a subject", 4000),
        record(1, b"Subject_\x81 x", 3000),
    ];
    let (mut one, mut two) = (Vec::<u8>::new(), Vec::<u8>::new());
    let mut errors = Vec::<u8>::new();
    let single = |sink: &mut Vec<u8>, errors: &mut Vec<u8>| -> std::result::Result<Ready, Off> {
        let mut outputs = ReportOutputs::single("6", OutputSink::Writer(sink), errors);
        checked(true, &[], &good, 10, 4000, plain, &mut outputs)
    };
    // Not set.
    let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors);
    assert_eq!(
        checked(false, &[], &good, 10, 4000, plain, &mut outputs).err(),
        Some(Off::NotSet)
    );
    drop(outputs);
    // Hit records, an observer-like second output: two formats; a file sink.
    let mut on_hits = |_: Arc<Vec<crate::report::PairwiseHit>>| {};
    let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors);
    outputs.hits = Some(&mut on_hits);
    assert_eq!(
        checked(true, &[], &good, 10, 4000, plain, &mut outputs).err(),
        Some(Off::Outputs)
    );
    drop(outputs);
    let mut errors2 = Vec::<u8>::new();
    let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors2);
    outputs.formats.push(FormatOutput {
        outfmt: "7",
        refusal: crate::api::local_blast::FormatRefusal::FailsCall,
        sink: OutputSink::Writer(&mut two),
    });
    assert_eq!(
        checked(true, &[], &good, 10, 4000, plain, &mut outputs).err(),
        Some(Off::Outputs)
    );
    drop(outputs);
    let path = std::path::Path::new("never-created");
    let mut outputs = ReportOutputs::single("6", OutputSink::File(path), &mut errors);
    assert_eq!(
        checked(true, &[], &good, 10, 4000, plain, &mut outputs).err(),
        Some(Off::Outputs)
    );
    drop(outputs);
    // Several batches (4,000 letters of 10,000): the subject values are made once.
    let ready = single(&mut one, &mut errors).expect("acts");
    assert!(!ready.one_batch);
    assert_eq!(
        ready
            .cache
            .as_ref()
            .and_then(|cache| cache.ids.as_ref())
            .map(Vec::len),
        Some(2)
    );
    // One batch (20,000 letters): no check, nothing made ahead, even with a refused id.
    let mut outputs = ReportOutputs::single("6", OutputSink::Writer(&mut one), &mut errors);
    let ready = checked(true, &[], &refused, 10, 20_000, plain, &mut outputs).expect("acts");
    assert!(ready.one_batch && ready.cache.is_none());
    // Several batches with a refused id: the switch acts, and each batch's report makes the
    // ids (and refuses the subjects that it reports, as without the switch).
    let ready = checked(true, &[], &refused, 10, 4000, plain, &mut outputs).expect("acts");
    assert!(!ready.one_batch);
    assert!(ready.cache.is_some_and(|cache| cache.ids.is_none()));
    // Several batches: the refusals of the options and the query lengths that could come
    // after the first batch.
    let scoring = [
        "-task",
        "blastn",
        "-reward",
        "2",
        "-penalty",
        "-3",
        "-gapopen",
        "1",
        "-gapextend",
        "1",
    ];
    assert_eq!(
        checked(true, &scoring, &good, 10, 4000, plain, &mut outputs).err(),
        Some(Off::Scoring)
    );
    // A `-dbsize` whose product could overflow 64-bit search space still streams: search space
    // overflow fails the batch after earlier reports (`setup_failure` -> `ErrorAfterBatchReports`)
    // in both stream on and stream off paths.
    assert!(checked(
        true,
        &["-dbsize", "9000000000000000000"],
        &good,
        10,
        4000,
        plain,
        &mut outputs
    )
    .is_ok_and(|ready| !ready.one_batch));
    // A negative CHUNK_SIZE streams: its refusals come after the reports of the batches before in
    // both paths (line RG item X-short-bat-2).
    assert!(checked(
        true,
        &[],
        &good,
        10,
        4000,
        SplitSizes::new(false, Some(-5), Some(-10)),
        &mut outputs
    )
    .is_ok());
    // BLAST XML, XML2 and JSON (outfmt 5, 16 and 15) are never streamed, with one batch
    // (20,000 letters) or several.
    for batch_size in [20_000, 4000] {
        for xml in [
            [BlastnOutputFormat::Xml],
            [BlastnOutputFormat::Xml2(
                crate::report::blast_xml2::Xml2Kind::Single16,
            )],
            [BlastnOutputFormat::Xml2(
                crate::report::blast_xml2::Xml2Kind::Single15,
            )],
        ] {
            assert_eq!(
                checked_as(&xml, true, &[], &good, 10, batch_size, plain, &mut outputs).err(),
                Some(Off::Xml)
            );
        }
    }
    // The BLAST archive (outfmt 11) is streamed, with several batches (4,000 letters) or one
    // (20,000 letters).
    let archive = [BlastnOutputFormat::Archive];
    assert!(
        checked_as(&archive, true, &[], &good, 10, 4000, plain, &mut outputs)
            .is_ok_and(|ready| !ready.one_batch)
    );
    assert!(
        checked_as(&archive, true, &[], &good, 10, 20_000, plain, &mut outputs)
            .is_ok_and(|ready| ready.one_batch)
    );
}
