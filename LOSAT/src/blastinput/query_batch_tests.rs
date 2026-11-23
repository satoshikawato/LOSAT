use super::{get_query_batch_size, BatchSizeMixer};
use crate::blastinput::seq_range::next_ranged_batch_end;

/// The end of the batch that starts at `start` when no record is skipped. The batch size
/// is a `TSeqPos` (32-bit unsigned): an `Int4` batch size converts to it.
fn next_query_batch_end(lengths: &[usize], start: usize, batch_size: u32) -> usize {
    next_ranged_batch_end(lengths, &vec![false; lengths.len()], start, batch_size)
}

// blast_input.cpp:134-171: a batch takes records until their letters reach the batch
// size; each searched query belongs to the batch that read it.
#[test]
fn every_searched_query_knows_where_its_batch_starts() {
    use super::QueryReading;
    use crate::blastinput::fasta_reader::FastaRecord;
    use crate::blastinput::seq_range::QueryInput;
    let records: Vec<FastaRecord> = [4usize, 7, 3, 12, 1]
        .iter()
        .enumerate()
        .map(|(i, &len)| FastaRecord::new(format!("q{i}"), b"", &vec![b'A'; len]))
        .collect();
    let input = QueryInput::whole(&records);
    let batches = super::read_query_batches(
        &records,
        &records,
        &input,
        QueryReading::Records,
        10,
        true,
        None,
    );
    assert!(batches.stop.is_none());
    assert_eq!(batches.starts, vec![(0, 0), (2, 2), (4, 4)]);
    let starts: Vec<(usize, usize)> = (0..5).map(|q| batches.start_of(q)).collect();
    assert_eq!(starts, vec![(0, 0), (0, 0), (2, 2), (2, 2), (4, 4)]);
}

// Probe K4f of line RG (NCBI 2.17.0, `-matrix blosum62`): the strategy's last request is that
// of the batch of the run's first searched query, or of its second where the report is held;
// a held run of one searched query (also one that an empty batch stops) writes every request.
#[test]
fn the_xml_matrix_failure_ends_the_requests_at_its_batch() {
    use super::QueryReading;
    use crate::blastinput::fasta_reader::FastaRecord;
    use crate::blastinput::seq_range::QueryInput;
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let batches_of = |lengths: &[usize], batch_size: u32| {
        let records: Vec<FastaRecord> = lengths
            .iter()
            .enumerate()
            .map(|(i, &len)| FastaRecord::new(format!("q{i}"), b"", &vec![b'A'; len]))
            .collect();
        let input = QueryInput::whole(&records);
        super::read_query_batches(
            &records,
            &records,
            &input,
            QueryReading::Records,
            batch_size,
            true,
            None,
        )
    };
    let (first, held) = (
        XmlMatrixFailure::FirstQuery,
        XmlMatrixFailure::SecondQueryOrEnd { second_id: Some(1) },
    );
    // Three batches of one query; one batch of three.
    let three = batches_of(&[4, 7, 3], 1);
    assert_eq!(three.xml_matrix_request_end(first), Some(1));
    assert_eq!(three.xml_matrix_request_end(held), Some(2));
    let one_batch = batches_of(&[4, 7, 3], 100);
    assert_eq!(one_batch.xml_matrix_request_end(first), Some(3));
    assert_eq!(one_batch.xml_matrix_request_end(held), Some(3));
    // One query; one query and a batch without letters, which stops the run.
    let alone = batches_of(&[4], 1);
    assert_eq!(alone.xml_matrix_request_end(first), Some(1));
    assert_eq!(alone.xml_matrix_request_end(held), None);
    let stopped = batches_of(&[4, 0], 1);
    assert!(stopped.stop.is_some() && stopped.searched_done == 1);
    assert_eq!(stopped.xml_matrix_request_end(first), Some(1));
    assert_eq!(stopped.xml_matrix_request_end(held), None);
}

/// The records of `fasta` read with `-parse_deflines` (BLASTP's reader).
fn parsed_records(fasta: &[u8]) -> Vec<crate::blastinput::fasta_reader::FastaRecord> {
    use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
    let config = ReaderConfig::query("BLASTP", true, false).with_parse_deflines(true);
    let mut source = FastaInputSource::from_bytes(fasta, config);
    read_all(&mut source, &mut |_| Ok(())).unwrap()
}

/// The scope keys of the local id `id`.
fn local_keys(id: &str) -> Vec<u64> {
    use crate::seqid::model::{ObjectId, ParsedSeqId};
    let mut keys = Vec::new();
    crate::seqid::defline_set::scope_keys(
        &ParsedSeqId::Local(ObjectId::Str(id.as_bytes().to_vec())),
        &mut keys,
    );
    keys
}

// blast_fasta_input.cpp:364-367 and 395, blast_input.cpp:146-155 (review 2 of line DBID
// loadord, H2 and M2): the scope at a batch's prefetch holds the batch's input records with
// their `Query_N` of the run, the record that `-query_loc` drops included, and with
// `-parse_deflines` the ids of the deflines alone.
#[test]
fn the_scope_of_a_batch_holds_its_input_records() {
    use super::QueryScope;
    let config = crate::blastinput::fasta_reader::ReaderConfig::query("BLASTP", true, false);
    let mut source = crate::blastinput::fasta_reader::FastaInputSource::from_bytes(
        b">qa\nACDEFGHIKL\n>qb\nACD\n>qc\nACDEFGHIKL\n",
        config,
    );
    let records = crate::blastinput::fasta_reader::read_all(&mut source, &mut |_| Ok(())).unwrap();
    let cut = crate::blastinput::seq_range::cut_queries(
        &records,
        &crate::blastinput::seq_range::SequenceRange { from: 5, to: 8 },
    );
    // `qb` (3 letters) is dropped; batches of 10 letters: `qa`, then `qb` and `qc`.
    let batches = cut.input.batches(10);
    assert_eq!(batches.len(), 2);
    assert_eq!(
        (batches[1].input.clone(), batches[1].searched.clone()),
        (1..3, 1..2)
    );
    assert_eq!(cut.records[1].local_id, "Query_3");
    let second = QueryScope::of_shown_records(&records[batches[1].input.clone()]);
    assert!(second.has_key(&local_keys("Query_2")));
    assert!(second.has_key(&local_keys("Query_3")));
    assert!(!second.has_key(&local_keys("Query_1")));
    // outfmt 5: the batches before too.
    let run = QueryScope::of_shown_records(&records[..batches[1].input.end]);
    assert!(run.has_key(&local_keys("Query_1")));
    let believed = parsed_records(b">qa x\nAC\n>qb x\nAC\n");
    let scope = QueryScope::of_shown_records(&believed);
    assert!(scope.has_key(&local_keys("qa")) && scope.has_key(&local_keys("qb")));
    assert!(!scope.has_key(&local_keys("Query_1")));
}

// blast_format.cpp:2295-2313 (series Q of line RG, Q1c): a batch's scope holds its own queries;
// a call that searches from a later query numbers them from there.
#[test]
fn query_scopes_number_the_searched_queries_of_each_batch() {
    use super::{QueryScopeSpan, QueryScopes};
    let records = parsed_records(b">a x\nAC\n>b x\nAC\n>a y\nAC\n>c x\nAC\n");
    let mut scopes = QueryScopes::new(QueryScopeSpan::Batch, "BLASTP");
    scopes.add_batch(&records[..2], &[false; 2], 0).unwrap();
    scopes.add_batch(&records[2..], &[false; 2], 2).unwrap();
    assert!(scopes.of_query(0).has_key(&local_keys("b")));
    assert!(scopes.of_query(1).has_key(&local_keys("b")));
    assert!(!scopes.of_query(2).has_key(&local_keys("b")));
    assert!(scopes.of_query(3).has_key(&local_keys("c")));
    assert!(!scopes.of_query(1).has_key(&local_keys("c")));
    let from = scopes.numbered_from(2);
    assert!(from.of_query(0).has_key(&local_keys("c")));
    assert!(!from.of_query(1).has_key(&local_keys("b")));
    // `-subject` and outfmt 5: one scope, in which the second `a` is refused.
    let mut run = QueryScopes::new(QueryScopeSpan::Run, "BLASTP");
    run.add_batch(&records[..2], &[false; 2], 0).unwrap();
    let error = run
        .add_batch(&records[2..], &[false; 2], 2)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("two query records have the Seq-id lcl|a"),
        "{error}"
    );
    // The ids of every query before the batches add theirs; a read ahead adds nothing again.
    let whole = QueryScopes::new(QueryScopeSpan::Run, "BLASTP").with_whole_run_of(&records);
    assert!(whole.of_query(0).has_key(&local_keys("c")));
    let mut ahead = QueryScopes::of_records(QueryScopeSpan::Run, "BLASTP", &records[1..])
        .unwrap()
        .read_ahead(Vec::new());
    ahead.add_batch(&records[1..], &[false; 3], 0).unwrap();
    // A refusal that the read met comes with the search's add of its batch (the second here).
    let mut ahead = QueryScopes::of_records(QueryScopeSpan::Run, "BLASTP", &records[..2])
        .unwrap()
        .read_ahead(vec![(1, anyhow::anyhow!("refused"))]);
    ahead.add_batch(&records[..2], &[false; 2], 0).unwrap();
    let refused = ahead.add_batch(&records[2..3], &[false], 2).unwrap_err();
    assert_eq!(refused.to_string(), "refused");
    ahead.add_batch(&records[3..], &[false], 3).unwrap();
    // A batch's scope has no ids of the whole run.
    let mut batch = QueryScopes::new(QueryScopeSpan::Batch, "BLASTP").with_whole_run_of(&records);
    batch.add_batch(&records[..2], &[false; 2], 0).unwrap();
    assert!(!batch.of_query(0).has_key(&local_keys("c")));
}

// The refusal of two queries of one scope stops the run at the batch of the second, after the
// batches before, which NCBI searches first (`QueryBatches::scope_refusal`; line WQA D10
// round 4, 2026-11-20).
#[test]
fn two_queries_of_one_scope_stop_the_run_at_the_second() {
    use super::{QueryReading, QueryScopeSpan, QueryScopes};
    use crate::blastinput::seq_range::QueryInput;
    let records = parsed_records(b">a x\nAC\n>a y\nAC\n");
    let input = QueryInput::whole(&records);
    for (span, refused) in [(QueryScopeSpan::Batch, false), (QueryScopeSpan::Run, true)] {
        let batches = super::read_query_batches(
            &records,
            &records,
            &input,
            QueryReading::Records,
            1,
            true,
            Some(QueryScopes::new(span, "BLASTP")),
        );
        assert_eq!(batches.stop.is_some(), refused);
        assert_eq!(batches.scope_refusal, refused);
        assert!(batches.stop_batch.is_none());
        assert_eq!(batches.starts.len(), if refused { 1 } else { 2 });
        assert_eq!(batches.searched_done, if refused { 1 } else { 2 });
    }
    let batches = super::read_query_batches(
        &records,
        &records,
        &input,
        QueryReading::Records,
        10,
        true,
        Some(QueryScopes::new(QueryScopeSpan::Batch, "BLASTP")),
    );
    assert!(batches.scope_refusal && batches.starts.is_empty());
    let (_, refusal) = batches.stop.expect("the refusal stops the first batch");
    assert!(refusal
        .to_string()
        .contains("two query records have the Seq-id lcl|a"));
}

/// The refusal of a batch with the id of an earlier record (line WQA D10 round 5, 2026-11-20;
/// round 6): the batches `batches` added to one scope of the run.
fn batch_refusal(
    program: &'static str,
    batches: &[&[u8]],
) -> Option<crate::blastinput::query_batch::ScopeRefusal> {
    use super::{QueryScopeSpan, QueryScopes, ScopeRefusal};
    let mut scopes = QueryScopes::new(QueryScopeSpan::Run, program);
    let mut first = 0;
    let mut refusal = None;
    for fasta in batches {
        let records = parsed_records(fasta);
        if let Err(error) = scopes.add_batch(&records, &vec![false; records.len()], first) {
            refusal = ScopeRefusal::of(&error).cloned();
            assert!(refusal.is_some(), "{error}");
            break;
        }
        first += records.len();
    }
    refusal
}

/// `batch_refusal`'s records that NCBI's set-up reads with their own letters.
fn setup_records(program: &'static str, batches: &[&[u8]]) -> Option<Vec<bool>> {
    batch_refusal(program, batches)
        .expect("a refusal")
        .setup_records()
        .map(<[bool]>::to_vec)
}

#[test]
fn a_refused_batch_knows_the_records_that_ncbis_set_up_reads_as_their_own() {
    let own = |flags: &[u8]| Some(flags.iter().map(|&flag| flag == 1).collect::<Vec<_>>());
    // NCBI's batch is LOSAT's (`BatchCount`): its set-up reads q1 and q2 with their own
    // letters, the repeat in the letters of the first q1.
    let batch: &[u8] = b">q1\nMKV\n>q2\nMKV\n>q1\nMK\n";
    assert_eq!(setup_records("BLASTP", &[batch]), own(&[1, 1, 0]));
    let refusal = batch_refusal("BLASTP", &[batch]).unwrap();
    assert!(refusal.to_string().contains("the Seq-id lcl|q1"));
    let batches: [&[u8]; 2] = [b">q0\nMKV\n", b">q3\nMK\n>q1\nMKV\n>q1\nMK\n"];
    assert_eq!(setup_records("BLASTP", &batches), own(&[1, 1, 0]));
    // A record after the repeat is set up with its own letters too
    // (`db.pd.q1c.bp.avgmid.subj.6`, review 4's `g2a`).
    let batch: &[u8] = b">q1\nMKV\n>q1\nMK\n>q2\nMKV\n";
    assert_eq!(setup_records("BLASTP", &[batch]), own(&[1, 0, 1]));
    // The repeat of a batch before, shorter or longer than its first record.
    let batches: [&[u8]; 2] = [b">q1\nMKV\n", b">q3\nMK\n>q1\nMKV\n"];
    assert_eq!(setup_records("BLASTP", &batches), own(&[1, 0]));
    assert_eq!(setup_records("BLASTN", &batches), own(&[1, 0]));
    let batches: [&[u8]; 2] = [b">q1\nMK\n", b">q3\nMK\n>q1\nMKV\n"];
    assert_eq!(setup_records("BLASTP", &batches), own(&[1, 0]));
    assert_eq!(setup_records("TBLASTX", &batches), own(&[1, 0]));
    let batch: &[u8] = b">q1\nACGT\n>q2\nACGT\n>q1\nACGTACGT\n";
    assert_eq!(setup_records("TBLASTX", &[batch]), own(&[1, 1, 0]));
    // BLASTN's set-up reads the repeat's own interval in the first record's letters (NCBI
    // 2.17.0: past their end, `CSeqVector_CI: cannot locate segment`): refused where those
    // are fewer.
    assert_eq!(setup_records("BLASTN", &[batch]), None);
    assert_eq!(setup_records("BLASTN", &batches), None);
    // Nothing with its own letters: no set-up (`db.pd.q1c.bp.avgdup2.subj.6`).
    assert_eq!(
        setup_records("BLASTP", &[b">q1\nMKV\n", b">q1\nMKV\n"]),
        None
    );
    // A record with a new id and an earlier one, or an earlier one's name with an accession.
    let batch: &[u8] = b">q1\nMKV\n>q2\nMKV\n>gi|5|lcl|q1\nMK\n";
    assert_eq!(setup_records("BLASTP", &[batch]), None);
    let batch: &[u8] = b">gb||NAMEA\nMKV\n>gb|AAA00007.1|NAMEA\nMK\n";
    assert_eq!(setup_records("BLASTP", &[batch]), None);
}

// NCBI BLAST+ 2.17.0 (`fix-reader-error/r6/probe` of line WQA, 2026-11-20): a repeat counts
// the letters of its first record in the scope, which the formatter empties after each batch
// but with `-subject` and outfmt 5; a record that `-query_loc` skips joins the scope and adds
// nothing.
#[test]
fn ncbis_reader_counts_a_repeat_by_its_first_records_letters() {
    use super::{BatchCount, QueryScopeSpan};
    let records = parsed_records(b">q1\nMKVLA\n>q1\nMK\n>q2\nMKV\n>q1\nMKVLAMKVLA\n");
    let lengths: Vec<usize> = records.iter().map(|record| record.seq().len()).collect();
    let none = vec![false; records.len()];
    let batches = |span: QueryScopeSpan, size: u32, skipped: &[bool]| {
        let mut count = Some(BatchCount::new(span, "BLASTP"));
        let mut batches = Vec::new();
        let mut start = 0;
        while start < records.len() {
            let batch = BatchCount::next_batch(
                count.as_mut(),
                &|record| records[record].believed_shared(),
                &lengths,
                skipped,
                start,
                size,
            );
            batches.push(batch);
            start = batch.0;
        }
        batches
    };
    // 5, then the second q1 by 5 (not 2): a batch of 10 letters ends there. The third q1
    // counts the first's 5 in the run's scope (3 + 5, the input ends first), its own 10 after
    // the scope was emptied (3 + 10).
    assert_eq!(
        batches(QueryScopeSpan::Run, 10, &none),
        [(2, true), (4, false)]
    );
    assert_eq!(
        batches(QueryScopeSpan::Batch, 10, &none),
        [(2, true), (4, true)]
    );
    assert_eq!(
        batches(QueryScopeSpan::Run, 11, &none),
        [(3, true), (4, false)]
    );
    assert_eq!(batches(QueryScopeSpan::Run, 14, &none), [(4, true)]);
    // The first q1 skipped: the second counts by its 5 letters all the same, and the skipped
    // record adds nothing.
    let skipped = [true, false, false, false];
    assert_eq!(
        batches(QueryScopeSpan::Run, 8, &skipped),
        [(3, true), (4, false)]
    );
    // Without -parse_deflines each record counts its own letters.
    let batch = BatchCount::next_batch(None, &|_| None, &lengths, &none, 0, 10);
    assert_eq!(batch, (3, true));
}

/// What `BatchCount`'s scope holds beside `QueryScopes`' with the run's scope (`-subject`,
/// outfmt 5), where both keep every query's keys (review 5 of line WQA D10, L2): run with
/// `--ignored --nocapture`, it prints the entries and the table bytes of each for 1,000,000
/// queries `>qN` of 2 letters.
#[test]
#[ignore]
fn batch_count_scope_of_a_million_queries() {
    use super::{BatchCount, QueryScopeSpan, QueryScopes};
    let mut text = Vec::new();
    for query in 0..1_000_000 {
        text.extend_from_slice(format!(">q{query}\nMK\n").as_bytes());
    }
    let records = parsed_records(&text);
    let mut count = BatchCount::new(QueryScopeSpan::Run, "BLASTP");
    for record in &records {
        count.count(record.believed(), record.seq().len(), false);
    }
    let mut scopes = QueryScopes::new(QueryScopeSpan::Run, "BLASTP");
    scopes
        .add_batch(&records, &vec![false; records.len()], 0)
        .expect("no repeated id");
    // A hashbrown table holds its capacity's entries and a control byte each.
    let bytes = |scope: &super::QueryScope| {
        scope.keys.capacity() * (8 + 1) + scope.first_records.capacity() * (24 + 1)
    };
    for (name, scope) in [
        ("BatchCount", &count.scope),
        ("QueryScopes", &scopes.scopes[0]),
    ] {
        println!(
            "{name}: keys {} first_records {} table bytes {}",
            scope.keys.len(),
            scope.first_records.len(),
            bytes(scope)
        );
    }
}

#[test]
fn a_refused_batch_gives_ncbis_set_up_error_or_the_refusal() {
    use super::{refused_batch_error, BatchSetUp, ScopeRefusal};
    let refusal = || -> anyhow::Error {
        batch_refusal("BLASTP", &[b">q1\nMKV\n>q2\nMKV\n>q1\nMK\n"])
            .unwrap()
            .into()
    };
    let native = || -> anyhow::Error {
        crate::blastinput::app::engine_error("The average subject length is too short")
    };
    let is_refusal = |error: &anyhow::Error| ScopeRefusal::of(error).is_some();
    // With -export_search_strategy the refusal comes before the batch's request: no set-up.
    let error = refused_batch_error(refusal(), true, |_| panic!("no set-up")).unwrap_err();
    assert!(is_refusal(&error));
    // NCBI's error of the set-up of the first records, which get it.
    let error = refused_batch_error(refusal(), false, |first| {
        assert_eq!(first, [true, true, false]);
        BatchSetUp::Fails(native())
    })
    .unwrap_err();
    assert!(error.is::<crate::cli::NativeError>(), "{error}");
    // A set-up that passes, or an error that is not NCBI's: the refusal (no XML failure takes it).
    let setups: [fn() -> BatchSetUp; 3] = [
        || BatchSetUp::Searched,
        || BatchSetUp::Unknown,
        || BatchSetUp::Fails(anyhow::anyhow!("LOSAT's")),
    ];
    for setup in setups {
        let error = refused_batch_error(refusal(), false, |_| setup()).unwrap_err();
        assert!(is_refusal(&error));
    }
    // No record of the batch to set up: the refusal without a set-up.
    let alone: anyhow::Error = batch_refusal("BLASTP", &[b">q1\nMKV\n", b">q1\nMKV\n"])
        .unwrap()
        .into();
    assert!(is_refusal(
        &refused_batch_error(alone, false, |_| panic!("no set-up")).unwrap_err()
    ));
}

/// The refusal of the first batch of `batches` (each with its records that a range skips) that
/// has one, under the run's XML failure `failure` with `span`, as the readers that read on past a
/// refusal that the failure takes add them: whether the failure takes it (`None`: no refusal), and
/// the failure with the query whose batch fails.
fn xml_refusal(
    failure: Option<crate::report::blast_xml::matrix::XmlMatrixFailure>,
    span: super::QueryScopeSpan,
    batches: &[(&[u8], &[bool])],
) -> (
    Option<bool>,
    Option<crate::report::blast_xml::matrix::XmlMatrixFailure>,
) {
    use super::{QueryScopes, ScopeRefusal};
    let mut scopes = QueryScopes::new(span, "BLASTP").with_xml_matrix_failure(failure);
    let mut first = 0;
    let mut taken = None;
    for &(fasta, skipped) in batches {
        let records = parsed_records(fasta);
        if let Err(error) = scopes.add_batch(&records, skipped, first) {
            let refusal = ScopeRefusal::of(&error).expect("a refusal");
            taken.get_or_insert(refusal.xml_takes());
            if !refusal.xml_takes() {
                break;
            }
        }
        first += skipped.iter().filter(|&&skipped| !skipped).count();
    }
    let at = failure.map(|failure| QueryScopes::held_xml_failure(Some(&scopes), failure));
    (taken, at)
}

/// Probe 2 of line RG-KB item 7 (NCBI BLAST+ 2.17.0, `-parse_deflines -matrix blosum62`): the
/// refusal of a repeat gives way to the run's XML failure where NCBI fails before it shows anything
/// that the repeated record's letters decide, and stays where NCBI reads them (`CSeqVector_CI`).
#[test]
fn the_xml_matrix_failure_takes_a_refusal_where_ncbi_reads_no_repeated_letters() {
    use super::QueryScopeSpan::{Batch, Run};
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let first = Some(XmlMatrixFailure::FirstQuery);
    let held = Some(XmlMatrixFailure::SecondQueryOrEnd { second_id: Some(1) });
    let at = |second_id| Some(XmlMatrixFailure::SecondQueryOrEnd { second_id });
    let (one, two, three) = (&[false][..], &[false; 2][..], &[false; 3][..]);
    let dup: &[u8] = b">q1\nMKV\n>q1\nMKV\n";
    let dup3: &[u8] = b">q1\nMKV\n>q2\nMKV\n>q1\nMKV\n";
    let (q1, q2): (&[u8], &[u8]) = (b">q1\nMKV\n", b">q2\nMKV\n");
    // FirstQuery: a repeat in the batch of the run's first searched query (`-db`: a scope per
    // batch; outfmt 5: one scope), not in a later batch (never read by NCBI).
    for span in [Batch, Run] {
        assert_eq!(xml_refusal(first, span, &[(dup, two)]).0, Some(true));
        assert_eq!(xml_refusal(first, span, &[(dup3, three)]).0, Some(true));
    }
    assert_eq!(
        xml_refusal(first, Run, &[(q1, one), (q1, one)]).0,
        Some(false)
    );
    // The first searched query is the run's, past the records that a range skips.
    let skipped_first: &[u8] = b">q0\nMKV\n>q1\nMKV\n>q1\nMKV\n";
    let skips = &[true, false, false][..];
    assert_eq!(
        xml_refusal(first, Run, &[(skipped_first, skips)]).0,
        Some(true)
    );
    // SecondQueryOrEnd (`-subject`: one scope): dup in one batch fails at the end; dup3 at q2.
    assert_eq!(
        xml_refusal(held, Run, &[(dup, two)]),
        (Some(true), at(None))
    );
    assert_eq!(
        xml_refusal(held, Run, &[(dup3, three)]),
        (Some(true), at(Some(1)))
    );
    // `q1 q1 q2`: the results of q1 are accumulated up to q2's batch.
    let dup_then: &[u8] = b">q1\nMKV\n>q1\nMKV\n>q2\nMKV\n";
    let accumulated = xml_refusal(held, Run, &[(dup_then, three)]);
    assert_eq!(accumulated, (Some(true), at(Some(2))));
    // A repeat after q2 in q2's batch.
    let after: &[u8] = b">q2\nMKV\n>q1\nMKV\n";
    let later = xml_refusal(held, Run, &[(q1, one), (after, two)]);
    assert_eq!(later, (Some(true), at(Some(1))));
    let repeat_of_second: &[u8] = b">q1\nMKV\n>q2\nMKV\n>q2\nMKV\n";
    assert_eq!(
        xml_refusal(held, Run, &[(repeat_of_second, three)]).0,
        Some(true)
    );
    // NCBI reads the repeat's letters: the first's id again in a later batch before q2 (dup b1),
    // or any repeat in a batch after q2's (dup3 b1, `q1 q2 | q1`).
    assert_eq!(
        xml_refusal(held, Run, &[(q1, one), (q1, one)]),
        (Some(false), at(None))
    );
    let dup3_b1 = xml_refusal(held, Run, &[(q1, one), (q2, one), (q1, one)]);
    assert_eq!(dup3_b1, (Some(false), at(Some(1))));
    let pair: &[u8] = b">q1\nMKV\n>q2\nMKV\n";
    assert_eq!(
        xml_refusal(held, Run, &[(pair, two), (q1, one)]).0,
        Some(false)
    );
    // A repeat that a range skips, a second Seq-id that the scope does not tell apart exactly, no
    // XML failure: the refusal.
    assert_eq!(
        xml_refusal(held, Run, &[(dup, &[false, true])]).0,
        Some(false)
    );
    assert_eq!(
        xml_refusal(first, Run, &[(dup, &[false, true])]).0,
        Some(false)
    );
    let inexact: &[u8] = b">gi|5 a\nMKV\n>gi|5|ref|NM_000001.1| b\nMKV\n";
    assert_eq!(xml_refusal(held, Run, &[(inexact, two)]).0, Some(false));
    assert_eq!(xml_refusal(None, Run, &[(dup, two)]), (Some(false), None));
    // No repeat: the held report fails with the second query of another Seq-id, as without the
    // scopes.
    assert_eq!(xml_refusal(held, Run, &[(pair, two)]), (None, held));
    let held = XmlMatrixFailure::SecondQueryOrEnd { second_id: Some(1) };
    assert_eq!(super::QueryScopes::held_xml_failure(None, held), held);
}

/// A refusal that the run's XML failure takes is searched where the batch's set-up passes with a
/// valid context; elsewhere it stops the run as other refusals do (`refused_batch_error`,
/// `QueryBatches::resolve_refusals`).
#[test]
fn a_refusal_that_the_xml_failure_takes_waits_for_the_set_up_of_its_batch() {
    use super::ScopeRefusal;
    use super::{refused_batch_error, BatchSetUp, QueryReading, QueryScopeSpan, QueryScopes};
    use crate::blastinput::seq_range::QueryInput;
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let taken = || -> anyhow::Error {
        let mut scopes = QueryScopes::new(QueryScopeSpan::Run, "BLASTP")
            .with_xml_matrix_failure(Some(XmlMatrixFailure::FirstQuery));
        let records = parsed_records(b">q1\nMKV\n>q2\nMKV\n>q1\nMK\n");
        scopes.add_batch(&records, &[false; 3], 0).unwrap_err()
    };
    let is_refusal = |error: &anyhow::Error| ScopeRefusal::of(error).is_some();
    let searched = |first: &[bool]| {
        assert_eq!(first, [true, true, false]);
        BatchSetUp::Searched
    };
    refused_batch_error(taken(), false, searched).unwrap();
    let error = refused_batch_error(taken(), true, |_| panic!("no set-up")).unwrap_err();
    assert!(is_refusal(&error));
    let error = refused_batch_error(taken(), false, |_| BatchSetUp::Unknown).unwrap_err();
    assert!(is_refusal(&error));
    let native = || crate::blastinput::app::engine_error("The average subject length is too short");
    let error = refused_batch_error(taken(), false, |_| BatchSetUp::Fails(native())).unwrap_err();
    assert!(error.is::<crate::cli::NativeError>(), "{error}");
    // Two batches of a scope each (`-db`, batches of 6 letters): the first's repeat waits for its
    // set-up, the second's (not the first searched query's) stops the reading.
    let records = parsed_records(b">q1\nMKV\n>q1\nMKV\n>q2\nMKV\n>q2\nMKV\n");
    let input = QueryInput::whole(&records);
    let read = || {
        let scopes = QueryScopes::new(QueryScopeSpan::Batch, "BLASTP")
            .with_xml_matrix_failure(Some(XmlMatrixFailure::FirstQuery));
        super::read_query_batches(
            &records,
            &records,
            &input,
            QueryReading::Records,
            6,
            true,
            Some(scopes),
        )
    };
    let mut batches = read();
    assert_eq!((batches.starts.len(), batches.xml_refusals.len()), (1, 1));
    assert!(batches.scope_refusal && batches.stop.is_some());
    let mut set_ups = Vec::new();
    batches.resolve_refusals(false, |start, _| {
        set_ups.push(start);
        BatchSetUp::Searched
    });
    assert_eq!(set_ups, [(0, 0), (2, 2)]);
    assert_eq!(batches.starts, [(0, 0)]);
    assert!(batches.scope_refusal && is_refusal(&batches.stop.as_ref().unwrap().1));
    // A first batch whose set-up does not pass: the run stops at it.
    let mut batches = read();
    batches.resolve_refusals(false, |_, _| BatchSetUp::Unknown);
    assert!(batches.starts.is_empty() && batches.xml_refusals.is_empty());
    assert_eq!((batches.input_done, batches.searched_done), (0, 0));
    assert!(batches.scope_refusal && is_refusal(&batches.stop.as_ref().unwrap().1));
    let mut batches = read();
    batches.resolve_refusals(false, |_, _| BatchSetUp::Fails(native()));
    assert!(batches.starts.is_empty() && !batches.scope_refusal);
    assert!(batches
        .stop
        .as_ref()
        .unwrap()
        .1
        .is::<crate::cli::NativeError>());
}

#[test]
fn a_call_gets_the_refusal_of_its_batch_where_it_adds_the_batch() {
    use super::{QueryScopeSpan, QueryScopes, ScopeRefusal};
    let refusal = batch_refusal("BLASTP", &[b">q1\nMKV\n", b">q1\nMKV\n"]).unwrap();
    let records = parsed_records(b">q1\nMKV\n");
    let mut scopes = QueryScopes::new(QueryScopeSpan::Run, "BLASTP").with_refusal(Some(refusal));
    let error = scopes.add_batch(&records, &[false], 0).unwrap_err();
    assert!(ScopeRefusal::of(&error).is_some());
    let mut scopes = QueryScopes::new(QueryScopeSpan::Run, "BLASTP").with_refusal(None);
    scopes.add_batch(&records, &[false], 0).unwrap();
}

// psiblast_app.cpp:468 and psiblast_aux_priv.cpp:636-649: a batch of one letter takes one
// query with letters and the queries without letters before it; `CPsiBlast` takes one query
// only (NCBI 2.17.0: `tk, empty, pol` stops with this error after tk's report; a last record
// without residues crashes it).
#[test]
fn a_psiblast_batch_holds_one_query() {
    use super::QueryReading;
    use crate::blastinput::fasta_reader::FastaRecord;
    use crate::blastinput::seq_range::QueryInput;
    let records = |lengths: &[usize]| -> Vec<FastaRecord> {
        lengths
            .iter()
            .enumerate()
            .map(|(i, &len)| FastaRecord::new(format!("q{i}"), b"", &vec![b'A'; len]))
            .collect()
    };
    let read = |records: &[FastaRecord]| {
        let input = QueryInput::whole(records);
        super::read_psiblast_query_batches(records, records, &input, QueryReading::Records, None)
    };
    let whole = records(&[4, 7, 3]);
    let batches = read(&whole);
    assert!(batches.stop.is_none());
    assert_eq!(batches.starts, vec![(0, 0), (1, 1), (2, 2)]);
    let joined = records(&[4, 0, 3]);
    let batches = read(&joined);
    assert_eq!(batches.starts, vec![(0, 0)]);
    assert_eq!((batches.input_done, batches.searched_done), (1, 1));
    let error = batches.stop.expect("two queries in one batch").1;
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .expect("an NCBI error");
    assert_eq!(native.exit, 3);
    assert_eq!(
        native.message,
        b"LOSAT engine error: PSI-BLAST only accepts one protein sequence as query\n"
    );
    for lengths in [&[4, 0][..], &[0][..]] {
        let alone = records(lengths);
        let batches = read(&alone);
        let error = batches.stop.expect("a query without letters alone").1;
        assert!(error.downcast_ref::<crate::cli::NativeError>().is_none());
        assert!(error
            .to_string()
            .contains("a query without letters alone in its batch"));
        assert_eq!(batches.starts.len(), lengths.len() - 1);
    }
    // blastp's batches of the same records: the record without letters is searched with the
    // next one, without an error.
    let input = QueryInput::whole(&joined);
    let blastp = super::read_query_batches(
        &joined,
        &joined,
        &input,
        QueryReading::Records,
        1,
        true,
        None,
    );
    assert!(blastp.stop.is_none());
    assert_eq!(blastp.starts, vec![(0, 0), (1, 1)]);
}

#[test]
fn the_query_that_reaches_the_batch_size_stays_in_the_batch() {
    assert_eq!(next_query_batch_end(&[19_999, 2, 1], 0, 20_000), 2);
    assert_eq!(next_query_batch_end(&[19_999, 2, 1], 2, 20_000), 3);
}

// NCBI BLAST+ 2.17.0 with `BATCH_SIZE` set (env parity probe h_probe of 2026-11-12, cases
// `BATCH_SIZE.*`): `NStr::StringToInt` with flags 0 takes a sign and no spaces, and stops at
// the line and position of each error.
#[test]
fn batch_size_strings_are_read_as_ncbis_string_to_int() {
    use std::ffi::OsStr;
    for (text, value) in [
        ("100", 100),
        ("+5", 5),
        ("-1", -1),
        ("0", 0),
        ("007", 7),
        ("2147483647", i32::MAX),
        ("-2147483648", i32::MIN),
        // NCBI reads up to a NUL (a loop over the `char` at each `str[pos]`,
        // c++/src/corelib/ncbistr.cpp:841).
        ("12\u{0}x", 12),
    ] {
        assert_eq!(
            get_query_batch_size(Some(OsStr::new(text))),
            Ok(value),
            "{text:?}"
        );
    }
    for (text, line, pos) in [
        ("abc", 861, 0),
        ("", 861, 0),
        (" 5", 861, 0),
        ("5 ", 870, 1),
        ("1.5", 870, 1),
        ("5x", 870, 1),
        ("+", 861, 1),
        ("1,000", 870, 1),
        ("0x10", 870, 1),
        ("99999999999999999999", 852, 18),
        ("-99999999999999999999", 852, 19),
        ("2147483648", 640, 0),
        ("-2147483649", 640, 0),
        ("\u{0}12", 861, 0),
    ] {
        let error = get_query_batch_size(Some(OsStr::new(text))).unwrap_err();
        assert_eq!(error.site().line, line, "{text:?}");
        let what = String::from_utf8(error.what()).unwrap();
        assert!(
            what.ends_with(&format!("(position {pos})\n")),
            "{text:?}: {what}"
        );
    }
    assert_eq!(get_query_batch_size(None), Ok(0));
}

#[test]
fn a_negative_batch_size_is_a_huge_tseqpos() {
    assert_eq!(
        next_query_batch_end(&[3000, 1999, 2, 7], 0, -1i32 as u32),
        4
    );
    assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 0, 1), 1);
}

#[test]
fn the_next_batch_takes_queries_until_the_batch_size_is_reached() {
    assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 0, 5000), 3);
    assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 3, 5000), 4);
    assert_eq!(next_query_batch_end(&[3000, 2000, 2], 0, 5000), 2);
}

// The values follow NCBI's CBatchSizeMixer (blast_app_util.cpp:65-81).
#[test]
fn the_mixer_is_ncbis() {
    let mut mixer = BatchSizeMixer::new(1_000_000 - 1000);
    assert_eq!(mixer.batch_size(None), 5000);
    // 4999 hits: ratio 1.0, 1e6 residues, above the maximum.
    assert_eq!(mixer.batch_size(Some(4999)), 999_000);
    // 999 hits in 999000 residues: ratio 0.001 (reset after the clamp), 1e9, clamped.
    assert_eq!(mixer.batch_size(Some(998)), 999_000);
    // 9999 hits in 5000 residues: ratio 2.0, 500000 residues.
    let mut mixer = BatchSizeMixer::new(4_999_000);
    mixer.batch_size(None);
    assert_eq!(mixer.batch_size(Some(9999)), 500_000);
    // Mixed: 0.3 * (20000 / 500000) + 0.7 * 2.0 = 1.412 -> 708215.
    assert_eq!(mixer.batch_size(Some(19_999)), 708_215);
    // No hits: 1e6 / (1 / 5000) overflows Int4, INT_MIN on x86-64, so the smallest batch.
    let mut mixer = BatchSizeMixer::new(999_000);
    mixer.batch_size(None);
    assert_eq!(mixer.batch_size(Some(0)), 100);
    // A large subject set raises the target and the first batch, within the maximum.
    let mut mixer = BatchSizeMixer::new(4_999_000);
    mixer.set_target_hits(3_000_000_000_i64.min(i64::from(i32::MAX)) as i32);
    assert_eq!(mixer.batch_size(None), 4_999_000);
    let mut mixer = BatchSizeMixer::new(999_000);
    mixer.set_target_hits(60_000_000);
    assert_eq!(mixer.batch_size(None), 300_000);
}
