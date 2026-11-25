use super::*;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::archive::options::OptionInputs;
use crate::report::archive::OptionsHandle;
use std::sync::{Arc, Mutex};

const PROTEIN: RequestPolicy = RequestPolicy {
    query_protein: true,
    subject_protein: true,
    lcase: None,
};

#[derive(Clone, Default)]
struct Sink(Arc<Mutex<Vec<u8>>>);

impl crate::report::archive::strategy::StrategySink for Sink {
    fn write_request(
        &mut self,
        bytes: &[u8],
        _pieces: &[crate::report::serial::encoder::BufferPiece],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn std::io::Write,
        archive: Option<&crate::report::archive::ArchiveHandle>,
    ) -> std::io::Result<()> {
        for (_, text) in posts {
            match archive {
                Some(handle) => {
                    handle.write_post(diagnostics, crate::report::diag_log::ERROR, text)?
                }
                None => diagnostics.write_all(text)?,
            }
            diagnostics.flush()?;
        }
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(())
    }

    fn write_partial(
        &mut self,
        bytes: &[u8],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn std::io::Write,
        archive: Option<&crate::report::archive::ArchiveHandle>,
    ) -> std::io::Result<()> {
        for (_, text) in posts {
            match archive {
                Some(handle) => {
                    handle.write_post(diagnostics, crate::report::diag_log::ERROR, text)?
                }
                None => diagnostics.write_all(text)?,
            }
            diagnostics.flush()?;
        }
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(())
    }
}

fn export() -> (StrategyExport, Sink) {
    let sink = Sink::default();
    let export = StrategyExport::new(
        Box::new(sink.clone()),
        &OptionInputs::new(Program::Blastp, "blastp"),
        None,
        OptionsHandle::default(),
    );
    (export, sink)
}

/// A title that ASN.1 text does not write as it is: the request holds `#` and the byte is
/// posted (NCBI 2.17.0, probe RG K/k4b/probe case blastp_strat: `Error: (803.7) [blastp]
/// Blast4-request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.E.title`), one
/// error post that the run's archive saves as a message; a subject's title the same; a
/// request after the run's stop is not made and posts nothing.
#[test]
fn a_title_with_a_replaced_byte_is_written_with_a_post() {
    let good = [FastaRecord::new("Query_1", b"q", b"MKV")];
    let bad = [FastaRecord::new("Query_2", b"b\x7fta", b"MKV")];
    let requests = RunRequests::new::<FastaRecord>(PROTEIN, None);
    let (strategy, sink) = export();
    requests.write(&strategy, 0, &good, None).unwrap();
    strategy.note_searched(0..1, None).unwrap();
    // A request without a replaced byte: its mark alone, which shows nothing.
    let marks = strategy.take_posts(usize::MAX);
    assert_eq!(marks.len(), 1);
    assert!(!marks[0].1.has_text());
    let mut nothing = Vec::new();
    marks[0].1.write(&mut nothing, None).unwrap();
    assert!(nothing.is_empty());
    let handle = crate::report::archive::ArchiveHandle::new(
        crate::report::archive::ArchiveEncoding::Text,
        crate::report::archive::ArchiveRequest {
            program: "blastp",
            application: "blastp",
            service: "plain",
            algorithm: Vec::new(),
            database: None,
            step_number: None,
        },
        OptionsHandle::default(),
    );
    requests.write(&strategy, 1, &bad, None).unwrap();
    strategy.note_searched(1..2, None).unwrap();
    // The posts wait for the batch's start, by its first query.
    assert!(strategy.take_posts(1).is_empty());
    let taken = strategy.take_posts(2);
    assert_eq!(taken.len(), 1);
    assert_eq!(taken[0].0, 1);
    let mut diagnostics = Vec::new();
    taken[0].1.write(&mut diagnostics, Some(&handle)).unwrap();
    let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
    assert!(text.contains("title \"b#ta\""), "{text}");
    let post = "Error: (803.7) [blastp] Blast4-request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.E.title\nBad char [0x7F] in string at line 21\nb\x7fta\n";
    assert_eq!(String::from_utf8(diagnostics.clone()).unwrap(), post);
    assert_eq!(
        handle.take_messages(),
        vec![(
            crate::report::diag_log::ERROR,
            b"Error: (803.7) [blastp] Blast4-request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.E.title Bad char [0x7F] in string at line 21 b ta".to_vec()
        )]
    );
    // After a stop no request is made.
    strategy.stop_after(2);
    let written = sink.0.lock().unwrap().len();
    requests.write(&strategy, 2, &bad, None).unwrap();
    assert_eq!(sink.0.lock().unwrap().len(), written);
    assert!(strategy.take_all_posts().is_empty());
    // A subject's title: posted with each request that holds it.
    let subjects = [FastaRecord::new("Subject_1", b"s\x7f", b"MKV")];
    let requests = RunRequests::new(PROTEIN, Some(&subjects[..]));
    let (strategy, sink) = export();
    requests.write(&strategy, 0, &good, None).unwrap();
    let mut diagnostics = Vec::new();
    strategy
        .take_all_posts()
        .write(&mut diagnostics, None)
        .unwrap();
    let text = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
    assert!(text.contains("title \"s#\""), "{text}");
    let posted = String::from_utf8(diagnostics).unwrap();
    assert!(
        posted.starts_with(
            "Error: (803.7) [blastp] Blast4-request.body.queue-search.subject.sequences.E.descr.E.title\n"
        ),
        "{posted}"
    );
}

#[test]
fn the_lower_case_runs_take_the_frames_of_the_program() {
    let records = [
        FastaRecord::new("Query_1", b"a", b"acGTtt"),
        FastaRecord::new("Query_2", b"b", b"ACGT"),
    ];
    let frames = |lcase: LcaseFrames| -> Vec<(Frame, Vec<(i32, i32)>)> {
        lcase
            .masks(&records)
            .into_iter()
            .map(|mask| {
                assert_eq!(mask.id, "Query_1");
                (mask.frame, mask.intervals)
            })
            .collect()
    };
    assert_eq!(
        frames(LcaseFrames::NotSet),
        [(Frame::NotSet, vec![(0, 1), (4, 5)])]
    );
    assert_eq!(
        frames(LcaseFrames::PlusTwice),
        [(Frame::Plus1, vec![(0, 1), (0, 1), (4, 5), (4, 5)])]
    );
    assert_eq!(
        frames(LcaseFrames::PlusAndMinus),
        [
            (Frame::Plus1, vec![(0, 1)]),
            (Frame::Minus1, vec![(0, 1)]),
            (Frame::Plus1, vec![(4, 5)]),
            (Frame::Minus1, vec![(4, 5)]),
        ]
    );
}

/// A PSSM query's request holds the PSSM in place of the batch's queries, without a range or
/// `LCaseMask` parameters (`x_Process_Pssm`), whatever the policy and the records.
#[test]
fn a_pssm_query_replaces_the_queries_without_range_or_masks() {
    let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(
        &std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/db_search/db.ps.asn.rnd.6.asn"
        ))
        .unwrap(),
    )
    .unwrap();
    let requests = RunRequests::new::<FastaRecord>(
        RequestPolicy {
            lcase: Some(LcaseFrames::NotSet),
            ..PROTEIN
        },
        None,
    )
    .with_pssm_query(Some(std::sync::Arc::new(pssm.clone())));
    let queries = [FastaRecord::new("Query_1", b"q", b"mkvMKV")];
    let cut = RecordPlacement {
        offset: 4,
        length: 10,
    };
    requests.with_batch(&queries, Some(cut), |request| {
        assert!(matches!(request.queries, RequestQueries::Pssm(held) if *held == pssm));
        assert_eq!(request.required, None);
        assert!(request.lcase_masks.is_empty());
    });
}

#[test]
fn a_cut_first_query_gives_the_required_interval() {
    let requests = RunRequests::new::<FastaRecord>(PROTEIN, None);
    let queries = [FastaRecord::new("Query_1", b"q", b"MKV")];
    let required = |placement| requests.with_batch(&queries, placement, |r| r.required);
    assert_eq!(required(None), None);
    let whole = RecordPlacement {
        offset: 0,
        length: 3,
    };
    assert_eq!(required(Some(whole)), None);
    let cut = RecordPlacement {
        offset: 4,
        length: 10,
    };
    assert_eq!(required(Some(cut)), Some((4, 6)));
}
