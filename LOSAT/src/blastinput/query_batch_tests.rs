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
        XmlMatrixFailure::SecondQueryOrEnd,
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
/// the query whose batch fails (`QueryScopes::xml_failing_query`).
fn xml_refusal(
    failure: Option<crate::report::blast_xml::matrix::XmlMatrixFailure>,
    span: super::QueryScopeSpan,
    batches: &[(&[u8], &[bool])],
) -> (Option<bool>, Option<super::XmlFailingQuery>) {
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
    let at = failure.map(|failure| QueryScopes::xml_failing_query(Some(&scopes), failure));
    (taken, at)
}

/// Probe 2 of line RG-KB item 7 and its review (NCBI BLAST+ 2.17.0, `-parse_deflines -matrix
/// blosum62`): the refusal of a repeat gives way to the run's XML failure in every batch that
/// NCBI formats before the failure or fails in (`HeldIds::reaches`), wherever the repeat stands
/// in it, a repeat that a range skips too; a batch after the failing one keeps the refusal (never
/// read by NCBI, whose failure comes first).
#[test]
fn the_xml_matrix_failure_takes_a_refusal_in_every_batch_up_to_the_failing_one() {
    use super::QueryScopeSpan::{Batch, Run};
    use super::XmlFailingQuery;
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let first = Some(XmlMatrixFailure::FirstQuery);
    let held = Some(XmlMatrixFailure::SecondQueryOrEnd);
    let at = |query| Some(XmlFailingQuery { query });
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
    // The first's id again in a later batch before q2 (probe 2's dup b1, the review's `short`,
    // `equal` and `mixed`): NCBI sets it up in the first's letters and accumulates its results.
    assert_eq!(
        xml_refusal(held, Run, &[(q1, one), (q1, one)]),
        (Some(true), at(None))
    );
    let mixed = xml_refusal(held, Run, &[(dup, two), (b">q1\nMKVL\n>q3\nMK\n", two)]);
    assert_eq!(mixed, (Some(true), at(Some(3))));
    // Any repeat in a batch after q2's (dup3 b1, `q1 q2 | q1`): never read.
    let dup3_b1 = xml_refusal(held, Run, &[(q1, one), (q2, one), (q1, one)]);
    assert_eq!(dup3_b1, (Some(false), at(Some(1))));
    let pair: &[u8] = b">q1\nMKV\n>q2\nMKV\n";
    assert_eq!(
        xml_refusal(held, Run, &[(pair, two), (q1, one)]).0,
        Some(false)
    );
    // A repeat that a range skips is never set up, searched or formatted (review F3).
    assert_eq!(
        xml_refusal(held, Run, &[(dup, &[false, true])]),
        (Some(true), at(None))
    );
    assert_eq!(
        xml_refusal(first, Run, &[(dup, &[false, true])]).0,
        Some(true)
    );
    assert_eq!(xml_refusal(None, Run, &[(dup, two)]), (Some(false), None));
    // No repeat: the held report fails with the second query of another Seq-id, as without the
    // scopes.
    assert_eq!(xml_refusal(held, Run, &[(pair, two)]), (None, at(Some(1))));
    let held = XmlMatrixFailure::SecondQueryOrEnd;
    assert_eq!(
        super::QueryScopes::xml_failing_query(None, held),
        XmlFailingQuery { query: Some(1) }
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2392
/// The held report's Seq-id is the best of the record's ids, matched as `CSeq_id::Match` matches
/// them (review F2 of line RG-KB item 7, probe `subset_empty`): `ref|NP_000001.1|` after
/// `emb|CAA00001.1|ref|NP_000001.1|` is a repeat in the scope and of another Seq-id (its best id),
/// so the held report fails at it, not at the next record.
#[test]
fn the_held_report_matches_the_queries_best_ids() {
    use super::QueryScopeSpan::Run;
    use super::XmlFailingQuery;
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let held = Some(XmlMatrixFailure::SecondQueryOrEnd);
    let at = |query| Some(XmlFailingQuery { query });
    let subset: &[u8] = b">emb|CAA00001.1|ref|NP_000001.1| a\nMKV\n>ref|NP_000001.1| b\nMKV\n";
    let records = parsed_records(subset);
    let best: Vec<String> = records
        .iter()
        .map(|record| {
            let mut label = Vec::new();
            crate::seqid::label::fasta(&record.believed().unwrap().best, &mut label);
            String::from_utf8(label).unwrap()
        })
        .collect();
    assert_eq!(best, ["emb|CAA00001.1|", "ref|NP_000001.1|"]);
    let empty_after = [
        (subset, &[false; 2][..]),
        (b">q2\n".as_slice(), &[false][..]),
    ];
    assert_eq!(
        xml_refusal(held, Run, &empty_after),
        (Some(true), at(Some(1)))
    );
    // An accession with and without its version match.
    let versions: &[u8] = b">ref|NP_000001.1| a\nMKV\n>ref|NP_000001| b\nMKV\n>q2\nMKV\n";
    assert_eq!(
        xml_refusal(held, Run, &[(versions, &[false; 3])]).1,
        at(Some(2))
    );
}

/// Where the run stops (probe K/k4f, NCBI 2.17.0): the batch of the first searched query, or
/// of the second when the report is held (`-subject` 13-16: 3 batches of one query wrote
/// the requests of Query_1 and Query_2), or the end of a run of one query. With queries of one
/// Seq-id (probe 2 of line RG-KB item 7): the batch of the first query of another Seq-id, or the
/// end of a run without one.
#[test]
fn a_held_report_fails_with_the_second_query_or_at_the_end() {
    use super::{QueryScopes, XmlFailingQuery};
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let first = QueryScopes::xml_failing_query(None, XmlMatrixFailure::FirstQuery);
    let held = QueryScopes::xml_failing_query(None, XmlMatrixFailure::SecondQueryOrEnd);
    assert!(first.fails_in(0..1) && first.fails_in(0..3) && !first.fails_in(1..2));
    assert!(!held.fails_in(0..1) && held.fails_in(1..2) && held.fails_in(0..3));
    assert!(!held.fails_in(2..3) && !held.fails_in(1..1));
    assert!(held.fails_at_end(1) && !held.fails_at_end(0) && !held.fails_at_end(2));
    assert!(!first.fails_at_end(1) && !first.fails_at_end(0));
    assert_eq!((first.query(), held.query()), (Some(0), Some(1)));
    // `q1, q1, q2`: the third query's batch; `q1, q1`: the end, whatever was searched.
    let third = XmlFailingQuery { query: Some(2) };
    assert!(!third.fails_in(0..2) && third.fails_in(2..3) && third.fails_in(0..3));
    assert!(third.fails_at_end(2) && !third.fails_at_end(3));
    let end = XmlFailingQuery { query: None };
    assert!(!end.fails_in(0..3) && end.query().is_none());
    assert!(end.fails_at_end(1) && end.fails_at_end(3) && !end.fails_at_end(0));
}

/// A refused batch of repeats alone (review F1 of line RG-KB item 7): where the run's XML failure
/// takes its refusal, the program sets up the repeats within their earlier records' letters with
/// their own letters (`SetUpRead::InFirst`: NCBI's contexts of them have those lengths), whose
/// `Searched` gives way to the failure; where that set-up does not pass, NCBI's result would
/// depend on the earlier records' letters, and where every repeat has more letters than its
/// earlier record NCBI's set-up fails (`CSeqVector_CI: cannot locate segment`, probe `long_new`):
/// the refusal.
#[test]
fn a_batch_of_repeats_is_set_up_with_their_contexts_in_their_earlier_records() {
    use super::{refused_batch_error, BatchSetUp, QueryScopeSpan, QueryScopes, ScopeRefusal};
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let refusal_of = |failure: Option<XmlMatrixFailure>, batches: &[&[u8]]| -> anyhow::Error {
        let mut scopes =
            QueryScopes::new(QueryScopeSpan::Run, "BLASTP").with_xml_matrix_failure(failure);
        let mut first = 0;
        for fasta in batches {
            let records = parsed_records(fasta);
            if let Err(error) = scopes.add_batch(&records, &vec![false; records.len()], first) {
                return error;
            }
            first += records.len();
        }
        panic!("no refusal")
    };
    let is_refusal = |error: &anyhow::Error| ScopeRefusal::of(error).is_some();
    let held = Some(XmlMatrixFailure::SecondQueryOrEnd);
    let native = || crate::blastinput::app::engine_error("The average subject length is too short");
    // The review's `short` and `equal`: the repeat in a batch of its own.
    let short: [&[u8]; 2] = [b">q1\nMKVLA\n", b">q1\nMKV\n"];
    let equal: [&[u8]; 2] = [b">q1\nMKV\n", b">q1\nAKV\n"];
    for batches in [short, equal] {
        let searched = |records: &[bool]| {
            assert_eq!(records, [true]);
            BatchSetUp::Searched
        };
        refused_batch_error(refusal_of(held, &batches), false, searched).unwrap();
        for fails in [false, true] {
            let setup = |_: &[bool]| match fails {
                true => BatchSetUp::Fails(native()),
                false => BatchSetUp::Unknown,
            };
            let error = refused_batch_error(refusal_of(held, &batches), false, setup).unwrap_err();
            assert!(is_refusal(&error), "{error}");
        }
        let no_set_up = |_: &[bool]| -> BatchSetUp { panic!("no set-up") };
        let error = refused_batch_error(refusal_of(held, &batches), true, no_set_up).unwrap_err();
        assert!(is_refusal(&error));
        let error = refused_batch_error(refusal_of(None, &batches), false, no_set_up).unwrap_err();
        assert!(is_refusal(&error));
    }
    // `long_new`: longer than the earlier record and alone in its batch: NCBI's `SetupQueries_OMF`
    // fails with `CSeqVector_CI: cannot locate segment at 2` (exit 3) regardless of XML matrix
    // failure, except under `-export_search_strategy` (which stays refused).
    let long: [&[u8]; 2] = [b">q1\nMK\n", b">q1\nMKV\n"];
    let error =
        refused_batch_error(refusal_of(held, &long), false, |_| panic!("no set-up")).unwrap_err();
    let native_err = error.downcast_ref::<crate::cli::NativeError>().unwrap();
    assert_eq!(native_err.exit, 3);
    assert_eq!(
        native_err.message.as_slice(),
        format!(
            "{} Warning: CSeqVector_CI: cannot locate segment at 2 \n",
            crate::wording::ENGINE_ERROR
        )
        .as_bytes()
    );
    let error =
        refused_batch_error(refusal_of(None, &long), false, |_| panic!("no set-up")).unwrap_err();
    assert_eq!(
        error
            .downcast_ref::<crate::cli::NativeError>()
            .unwrap()
            .exit,
        3
    );
    let error =
        refused_batch_error(refusal_of(held, &long), true, |_| panic!("no set-up")).unwrap_err();
    assert!(is_refusal(&error));
    // `mixed`: a record of its own beside it, whose set-up decides as in the first batch.
    let mixed: [&[u8]; 2] = [b">q1\nMK\n", b">q1\nMKV\n>q3\nMK\n"];
    let own = |records: &[bool]| {
        assert_eq!(records, [false, true]);
        BatchSetUp::Fails(native())
    };
    let error = refused_batch_error(refusal_of(held, &mixed), false, own).unwrap_err();
    assert!(error.is::<crate::cli::NativeError>(), "{error}");
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
        super::read_psiblast_query_batches(
            records,
            records,
            &input,
            QueryReading::Records,
            None,
            false,
        )
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

// psiblast `-phi_pattern` with one round (CLocalBlast): a record without letters joins the next
// query's batch, and a last record without letters stops the run with blastp's error (NCBI 2.17.0
// `PS/tmp/rev8a/ncbi3/er.first.0`, rc 0, and `er.last.0`, `Sequence contains no data`, rc 3).
#[test]
fn a_phi_batch_holds_the_queries_without_letters_before_the_next() {
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
        super::read_psiblast_query_batches(
            records,
            records,
            &input,
            QueryReading::Records,
            None,
            true,
        )
    };
    let joined = records(&[4, 0, 0, 3]);
    let batches = read(&joined);
    assert!(batches.stop.is_none());
    assert_eq!(batches.starts, vec![(0, 0), (1, 1)]);
    let last = records(&[4, 0]);
    let batches = read(&last);
    assert_eq!(batches.starts, vec![(0, 0)]);
    let error = batches.stop.expect("a last query without letters").1;
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .expect("an NCBI error");
    assert_eq!(native.exit, 3);
    assert_eq!(
        native.message,
        b"LOSAT engine error: Warning: Sequence contains no data \n"
    );
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

/// `QueryScopes::xml_show` of every query of `fasta`, added in batches of `batch` records to the
/// scopes of `span`, as (label, joins_held), and `held_start_of` of the whole run.
fn xml_shows(
    span: super::QueryScopeSpan,
    fasta: &[u8],
    batch: usize,
) -> (Vec<(Option<usize>, bool)>, usize) {
    let records = parsed_records(fasta);
    let mut scopes = super::QueryScopes::new(span, "BLASTP");
    for (n, chunk) in records.chunks(batch).enumerate() {
        scopes
            .add_batch(chunk, &vec![false; chunk.len()], n * batch)
            .unwrap();
    }
    let shows = (0..records.len())
        .map(|query| {
            let show = super::QueryScopes::xml_show(Some(&scopes), query);
            (show.label, show.joins_held)
        })
        .collect();
    let start = super::held_start_of(Some(scopes.numbered_from(0)), records.len());
    (shows, start)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2413
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:713-719
/// Probes of line RG-KB item 7c (NCBI BLAST+ 2.17.0, `-parse_deflines`, blastp, tblastn, blastx and
/// tblastx alike): a result shows the first query of its batch whose best Seq-id matches its own,
/// with `-db` and `-subject`, in every XML family format; with `-subject` the XML2 and JSON
/// documents join every result whose Seq-id matches that of the held report's first, across
/// batches too (`BATCH_SIZE` 1), and a result of another Seq-id leads the next report.
#[test]
fn xml_results_show_and_join_the_queries_whose_seq_ids_match() {
    use super::QueryScopeSpan::{Batch, Run};
    const NO: (Option<usize>, bool) = (None, false);
    let vpair: &[u8] = b">ref|NP_000001.1| first\nMKV\n>ref|NP_000001| second\nMKVL\n";
    for span in [Run, Batch] {
        assert_eq!(xml_shows(span, vpair, 2), (vec![NO, (Some(0), true)], 0));
        // Batches of one record: the label is the batch's own, the report is still joined.
        assert_eq!(xml_shows(span, vpair, 1), (vec![NO, (None, true)], 0));
    }
    // Both orders, PDB chains, accession case with a version, GenBank and EMBL.
    for pair in [
        &b">ref|NP_000001| second\nMKVL\n>ref|NP_000001.1| first\nMKV\n"[..],
        b">pdb|1ABC|A first\nMKV\n>pdb|1ABC| second\nMKVL\n",
        b">ref|NP_000001.1| first\nMKV\n>ref|np_000001| second\nMKVL\n",
        b">gb|AAB00001.1| first\nMKV\n>emb|AAB00001| second\nMKVL\n",
    ] {
        assert_eq!(xml_shows(Run, pair, 2), (vec![NO, (Some(0), true)], 0));
    }
    // Three records: the third, of another Seq-id, leads the last report (`vpair3a`); one of the
    // first's after another's shows the first and leads a report of its own (`vpair3b`).
    let vpair3a: &[u8] =
        b">ref|NP_000001.1| first\nMKV\n>ref|NP_000001| second\nMKVL\n>ref|NP_000002.1| third\nMKV\n";
    assert_eq!(
        xml_shows(Run, vpair3a, 3),
        (vec![NO, (Some(0), true), NO], 2)
    );
    assert_eq!(xml_shows(Run, vpair3a, 1), (vec![NO, (None, true), NO], 2));
    let vpair3b: &[u8] =
        b">ref|NP_000001.1| first\nMKV\n>ref|NP_000002.1| third\nMKV\n>ref|NP_000001| second\nMKVL\n";
    assert_eq!(
        xml_shows(Run, vpair3b, 3),
        (vec![NO, NO, (Some(0), false)], 2)
    );
    assert_eq!(xml_shows(Run, vpair3b, 1), (vec![NO, NO, NO], 2));
    // The third matches the second, not the first (`vpair3c`): it shows the second, the first
    // query of its batch that it matches, and leads a report of its own.
    let vpair3c: &[u8] =
        b">ref|NP_000001.1| first\nMKV\n>ref|NP_000001| second\nMKVL\n>ref|NP_000001.2| third\nMKV\n";
    for span in [Run, Batch] {
        assert_eq!(
            xml_shows(span, vpair3c, 3),
            (vec![NO, (Some(0), true), (Some(1), false)], 2)
        );
        assert_eq!(xml_shows(span, vpair3c, 1), (vec![NO, (None, true), NO], 2));
    }
    // A shown query keeps its own names, the second of `vpair3c` too when the third is written
    // (NCBI: `Iteration_query-def` `second`); when a stopping run drops the third (`reports = 2`),
    // or when a second query joins across batches under `-subject` XML2/JSON (`holds_back = true`),
    // its names are unshown.
    let unshown = |fasta: &[u8], batch: usize, reports: usize, holds_back: bool| {
        let records = parsed_records(fasta);
        let mut scopes = super::QueryScopes::new(Run, "BLASTP");
        for (n, chunk) in records.chunks(batch).enumerate() {
            scopes
                .add_batch(chunk, &vec![false; chunk.len()], n * batch)
                .unwrap();
        }
        let call = Some(scopes.numbered_from(0));
        (0..reports)
            .map(|query| super::names_unshown_in(call, reports, holds_back, query))
            .collect::<Vec<_>>()
    };
    assert_eq!(unshown(vpair3c, 3, 3, false), [false, false, true]);
    assert_eq!(unshown(vpair3c, 3, 3, true), [false, false, true]);
    assert_eq!(unshown(vpair3c, 3, 2, true), [false, true]);
    assert_eq!(unshown(vpair3a, 3, 3, false), [false, true, false]);
    assert_eq!(unshown(vpair, 1, 2, true), [false, true]);
    assert_eq!(unshown(vpair, 1, 2, false), [false, false]);
    // Cross-batch join where a third query in the second batch labels itself with the second:
    // the second is shown when the third is written, unshown when the third is held unwritten.
    let records = parsed_records(vpair3c);
    let mut split = super::QueryScopes::new(Run, "BLASTP");
    split.add_batch(&records[..1], &[false], 0).unwrap();
    split.add_batch(&records[1..], &[false, false], 1).unwrap();
    let call = Some(split.numbered_from(0));
    assert!(!super::names_unshown_in(call, 3, true, 1));
    assert!(super::names_unshown_in(call, 2, true, 1));
    // A call that starts at a later batch numbers the label and the held report's start as it
    // numbers its queries.
    let later: &[u8] = b">ref|NP_000002.1| a\nMKV\n>ref|NP_000003.1| b\nMKV\n\
        >ref|NP_000001.1| c\nMKV\n>ref|NP_000001| d\nMKVL\n";
    let records = parsed_records(later);
    let mut scopes = super::QueryScopes::new(Run, "BLASTP");
    scopes.add_batch(&records[..2], &[false; 2], 0).unwrap();
    scopes.add_batch(&records[2..], &[false; 2], 2).unwrap();
    let call = scopes.numbered_from(2);
    let show = call.xml_show(1);
    assert_eq!((show.label, show.joins_held), (Some(0), true));
    assert_eq!(super::held_start_of(Some(call), 2), 0);
    assert_eq!(super::held_start_of(Some(scopes.numbered_from(0)), 4), 2);
    assert_eq!(super::held_start_of(None, 4), 3);
    assert_eq!(super::xml_show_of(None, 3), super::XmlShow::default());
    // The pairs that the scope holds as one sequence stay refused: NCBI searches the first
    // record's letters for the second (accession, `gnl` and `lcl` case, GenBank and EMBL of one
    // version).
    for pair in [
        &b">ref|NP_000001.1| first\nMKV\n>ref|np_000001.1| second\nMKVL\n"[..],
        b">gnl|db|tag1 first\nMKV\n>gnl|DB|TAG1 second\nMKVL\n",
        b">lcl|q1 first\nMKV\n>lcl|Q1 second\nMKVL\n",
        b">gb|AAB00001.1| first\nMKV\n>emb|AAB00001.1| second\nMKVL\n",
    ] {
        let records = parsed_records(pair);
        let mut scopes = super::QueryScopes::new(Run, "BLASTP");
        assert!(scopes.add_batch(&records, &[false; 2], 0).is_err());
    }
}

/// Only an accession with and without its version, in any case, or a PDB molecule with and
/// without its chain may be an earlier id of the scope in another form
/// (`QueryScope::add_version_form`): two chains of one molecule, two versions or two accessions
/// never are, nor ids of a scope that has none of the other form. A scope that takes another's
/// ids takes their forms (`QueryScope::extend`).
#[test]
fn only_another_form_of_a_scope_s_id_may_match_it() {
    let fasta: &[u8] = b">ref|NP_000001.1| a\nMKV\n>ref|NP_000002.1| b\nMKV\n>ref|np_000001| c\nMKV\n>pdb|1ABC|A d\nMKV\n>pdb|1ABC|B e\nMKV\n>pdb|1ABC| f\nMKV\n>ref|NP_000002.2| g\nMKV\n>h x\nMKV\n";
    let records = parsed_records(fasta);
    let ids = |n: usize| &records[n].believed().unwrap().ids;
    let mut scope = super::QueryScope::default();
    let others: Vec<bool> = (0..records.len())
        .map(|n| scope.add(ids(n), (n, 3), "BLASTP").unwrap())
        .collect();
    assert_eq!(
        others,
        [false, false, true, false, false, true, false, false]
    );
    let mut first = super::QueryScope::default();
    assert!(!first.add(ids(2), (0, 3), "BLASTP").unwrap());
    let mut first = super::QueryScope::default();
    first.add(ids(0), (0, 3), "BLASTP").unwrap();
    let mut whole = super::QueryScope::default();
    whole.extend(&first);
    assert!(whole.add_after_refusal(ids(2), (1, 3)));
}

/// A run whose queries' Seq-ids all differ keeps nothing per query for the XML family formats
/// (`HeldIds`): no label and no join, and the second query is NCBI's `F`.
#[test]
fn queries_of_different_seq_ids_keep_no_xml_lists() {
    use super::QueryScopeSpan::{Batch, Run};
    let fasta: &[u8] = b">ref|NP_000001.1| a\nMKV\n>ref|NP_000002.1| b\nMKV\n>c x\nMKV\n>pdb|1ABC|A d\nMKV\n>pdb|1ABC|B e\nMKV\n";
    let records = parsed_records(fasta);
    for (span, batch) in [(Run, 5), (Run, 1), (Batch, 2)] {
        let mut scopes = super::QueryScopes::new(span, "BLASTP");
        for (n, chunk) in records.chunks(batch).enumerate() {
            scopes
                .add_batch(chunk, &vec![false; chunk.len()], n * batch)
                .unwrap();
        }
        assert!(scopes.held.joins.is_empty() && scopes.held.labels.is_empty());
        assert_eq!(scopes.held.second, Some(1));
        assert_eq!(xml_shows(span, fasta, batch).1, 4);
    }
}

/// `SetupQueries_OMF` failures on refused batches (`PastFirst` / `Empty` / `< 3` nt in translated
/// queries) emit NCBI's exit-3 `CSeqVector_CI` / `Sequence contains no data` message across
/// programs, `-query_loc` and `-strand minus`, while `PSIBLAST` and mixed `BLASTN` batches stay
/// refused (line RG-KB item 7b).
#[test]
fn setup_queries_omf_failures_match_ncbi_across_ranges_strands_and_programs() {
    use super::{refused_batch_error, BatchSetUp, QueryScopeSpan, QueryScopes, ScopeRefusal};
    use crate::blastinput::seq_range::{record_interval, RecordInterval, SequenceRange};
    use crate::wording::ENGINE_ERROR;

    let run = |program: &'static str,
               range: Option<SequenceRange>,
               minus: bool,
               batches: &[&[u8]]|
     -> anyhow::Error {
        let mut scopes = QueryScopes::new(QueryScopeSpan::Run, program)
            .with_query_range(range)
            .with_minus_strand(minus);
        let mut first = 0;
        for fasta in batches {
            let records = parsed_records(fasta);
            let skipped: Vec<bool> = records
                .iter()
                .map(|r| {
                    range.is_some_and(|qr| {
                        record_interval(&qr, r.seq().len()) == RecordInterval::PastEnd
                    })
                })
                .collect();
            if let Err(error) = scopes.add_batch(&records, &skipped, first) {
                return error;
            }
            first += skipped.iter().filter(|&&s| !s).count();
        }
        panic!("no refusal")
    };

    let expect_native = |err: anyhow::Error, expected: &str| {
        let native = err
            .downcast_ref::<crate::cli::NativeError>()
            .unwrap_or_else(|| panic!("expected NativeError, got: {err}"));
        assert_eq!(native.exit, 3);
        assert_eq!(native.message.as_slice(), expected.as_bytes());
    };

    // `-query_loc`: mid-range (earlier record length minus range start), past-end (0), exact-end
    // (`Sequence contains no data`), and in-bounds (`InFirst`, stays refused without XML failure).
    let range_11_50 = Some(SequenceRange { from: 10, to: 49 });
    let err = refused_batch_error(
        run(
            "BLASTP",
            range_11_50,
            false,
            &[
                b">q1\nAAAAAAAAAAAAAAAAAAAA\n",
                b">q1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
            ],
        ),
        false,
        |_| panic!("no set-up"),
    )
    .unwrap_err();
    expect_native(
        err,
        &format!("{ENGINE_ERROR} Warning: CSeqVector_CI: cannot locate segment at 10 \n"),
    );

    let err = refused_batch_error(
        run(
            "BLASTX",
            Some(SequenceRange { from: 14, to: 49 }),
            false,
            &[
                b">q1\nACGTACGTAC\n",
                b">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n",
            ],
        ),
        false,
        |_| panic!("no set-up"),
    )
    .unwrap_err();
    expect_native(
        err,
        &format!("{ENGINE_ERROR} Warning: CSeqVector_CI: cannot locate segment at 0 \n"),
    );

    // `-strand minus` locates the reverse end first, failing at position 0.
    let err = refused_batch_error(
        run(
            "BLASTX",
            None,
            true,
            &[b">q1\nACGTACGTAC\n", b">q1\nACGTACGTACGTACGTACGT\n"],
        ),
        false,
        |_| panic!("no set-up"),
    )
    .unwrap_err();
    expect_native(
        err,
        &format!("{ENGINE_ERROR} Warning: CSeqVector_CI: cannot locate segment at 0 \n"),
    );

    // Multiple failing records (Empty + PastFirst + translated `< 3` nt).
    let err = refused_batch_error(
        run(
            "TBLASTX",
            None,
            false,
            &[
                b">q1\nACGTAC\n>q2\nACG\n",
                b">q3\n\n>q1\nACGTACGTAC\n>q4\nAT\n>q2\nACGTAC\n",
            ],
        ),
        false,
        |_| panic!("no set-up"),
    )
    .unwrap_err();
    expect_native(
        err,
        &format!(
            "{ENGINE_ERROR} Warning: Sequence contains no data \
             Warning: CSeqVector_CI: cannot locate segment at 6 \
             Warning: CSeqVector_CI: cannot locate segment at 3 \n"
        ),
    );

    // `BLASTN`: all-`PastFirst` emits `CSeqVector_CI`; mixed `Own` + `PastFirst` stays refused.
    let err = refused_batch_error(
        run(
            "BLASTN",
            None,
            false,
            &[b">q1\nACGTAC\n", b">q1\nACGTACGTAC\n"],
        ),
        false,
        |_| panic!("no set-up"),
    )
    .unwrap_err();
    expect_native(
        err,
        &format!("{ENGINE_ERROR} Warning: CSeqVector_CI: cannot locate segment at 6 \n"),
    );
    let err = refused_batch_error(
        run(
            "BLASTN",
            None,
            false,
            &[b">q1\nACGTAC\n", b">q1\nACGTACGTAC\n>q2\nACGTACGTAC\n"],
        ),
        false,
        |_| panic!("no set-up"),
    )
    .unwrap_err();
    assert!(ScopeRefusal::of(&err).is_some(), "{err}");

    // `PSIBLAST` crashes in NCBI 2.17.0 on `PastFirst`: stays refused.
    let err = refused_batch_error(
        run("PSIBLAST", None, false, &[b">q1\nMK\n", b">q1\nMKV\n"]),
        false,
        |_| panic!("no set-up"),
    )
    .unwrap_err();
    assert!(ScopeRefusal::of(&err).is_some(), "{err}");

    // `BatchSetUp::NoValidContext` (review 2 B1/B2): proceeds under XML matrix failure only when
    // `subjects_pass` is true.
    use crate::report::blast_xml::matrix::XmlMatrixFailure;
    let xml_refusal = |batches: &[&[u8]]| -> anyhow::Error {
        let mut scopes = QueryScopes::new(QueryScopeSpan::Run, "BLASTP")
            .with_xml_matrix_failure(Some(XmlMatrixFailure::SecondQueryOrEnd));
        let mut first = 0;
        for fasta in batches {
            let records = parsed_records(fasta);
            if let Err(error) = scopes.add_batch(&records, &vec![false; records.len()], first) {
                return error;
            }
            first += records.len();
        }
        panic!("no refusal")
    };
    // B2 (`Own` has no valid KA context):
    let b2: [&[u8]; 1] = [b">q1\nXXXXXX\n>q1\nMKVLAA\n"];
    refused_batch_error(xml_refusal(&b2), false, |_| BatchSetUp::NoValidContext {
        subjects_pass: true,
    })
    .unwrap();
    let err = refused_batch_error(xml_refusal(&b2), false, |_| BatchSetUp::NoValidContext {
        subjects_pass: false,
    })
    .unwrap_err();
    assert!(ScopeRefusal::of(&err).is_some(), "{err}");

    // B1 (`InFirst` repeat alone in later batch has no valid KA context):
    let b1: [&[u8]; 2] = [b">q1\nMKVLAA\n", b">q1\nXXXXXX\n"];
    refused_batch_error(xml_refusal(&b1), false, |_| BatchSetUp::NoValidContext {
        subjects_pass: true,
    })
    .unwrap();
    let err = refused_batch_error(xml_refusal(&b1), false, |_| BatchSetUp::NoValidContext {
        subjects_pass: false,
    })
    .unwrap_err();
    assert!(ScopeRefusal::of(&err).is_some(), "{err}");
}
