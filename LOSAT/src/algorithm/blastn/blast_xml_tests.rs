use super::*;
use crate::blastinput::seq_range::QueryInput;
use crate::report::blast_xml::{BlastXmlWriter, REFERENCE_MEGABLAST};
use crate::report::blast_xml2::framing::{prolog_bytes, TAIL15, TAIL16};
use crate::report::blast_xml2::{Xml2Carry, Xml2Kind, Xml2Run, Xml2Target, Xml2Writer};

// The refusal of a title that NCBI BLAST+ 2.17.0 crashes on (`names::Refusal`) where NCBI makes the
// report that shows it (remaining-gaps DECISIONS D-SO-6): queries without hits against one subject,
// so that the reports show their titles alone, written by one report call (LOSAT_X_RPTSTREAM's off
// path), by a call for each query batch and the epilog's (the stream), or by one call that an error
// of NCBI's stops.

const LETTERS: &[u8] = b"ACGTTGCAACGTTGCAACGTTGCA";
/// The title that NCBI's `x_CleanAndCompress` reads past the end of: the one title LOSAT refuses.
const PUNCTUATION: &[u8] = b", ,";
const REFUSAL: &str = "not supported by LOSAT's BLASTN";

fn records(prefix: &str, titles: &[&[u8]]) -> Vec<FastaRecord> {
    titles
        .iter()
        .enumerate()
        .map(|(index, title)| FastaRecord::new(format!("{prefix}_{}", index + 1), title, LETTERS))
        .collect()
}

fn run() -> XmlRun<'static> {
    XmlRun {
        program: "blastn",
        app: "blastn",
        version: "BLASTN 2.17.0+",
        reference: REFERENCE_MEGABLAST,
        db_name: b"",
        matrix: None,
        expect: 10.0,
        gap_open: 0,
        gap_extend: 0,
        reward_penalty: Some((1, -2)),
        filter: b"L;m;",
        db_num: 0,
        db_len: 0,
        db_len_overflow: false,
        positives: None,
        mask_char: b'N',
    }
}

/// How the report calls of a run take its queries.
#[derive(Clone, Copy, Debug)]
enum Calls<'c> {
    /// One call with the epilog (the reports after the batch loop).
    One,
    /// A call for each query batch of these sizes, then the epilog's call (`x_rpt_stream`).
    Batches(&'c [usize]),
    /// One call without the epilog: an error of NCBI's stops the run after these batches.
    Stops,
}

impl Calls<'_> {
    /// Each call's queries of `count`, whether it has the epilog and whether it stops the run.
    fn plan(self, count: usize) -> Vec<(std::ops::Range<usize>, bool, bool)> {
        match self {
            Calls::One => vec![(0..count, true, false)],
            Calls::Stops => vec![(0..count, false, true)],
            Calls::Batches(sizes) => {
                assert_eq!(sizes.iter().sum::<usize>(), count);
                let mut at = 0;
                let mut plan: Vec<_> = sizes
                    .iter()
                    .map(|&size| {
                        at += size;
                        (at - size..at, false, false)
                    })
                    .collect();
                plan.push((count..count, true, false));
                plan
            }
        }
    }
}

/// One report call of `queries` (numbered from `first` in the run) into `out`.
fn call<W: StructuredXml>(
    document: &mut W,
    queries: &[FastaRecord],
    (epilog, stops): (bool, bool),
    out: &mut Vec<u8>,
) -> Result<()> {
    let subjects = records("Subject", &[b"subject"]);
    let query_input = QueryInput::whole(queries);
    let lengths = vec![LETTERS.len(); queries.len()];
    let length = |_: usize| LETTERS.len();
    let xml = BlastnXml {
        run: run(),
        queries,
        query_records: crate::algorithm::blastn::pairwise::QueryRecords::new(
            queries,
            &query_input,
            queries,
            0,
        ),
        searched_lengths: &lengths,
        stats: vec![XmlStats::NONE; queries.len()],
        messages: vec![Vec::new(); queries.len()],
        masks: Vec::new(),
        searched: Vec::new(),
        subjects: XmlSubjects::Records {
            records: &subjects,
            length: &length,
            protein: false,
        },
        taxdb_notes: Default::default(),
        letters: SubjectLetters::new(&subjects, None),
        ungapped: false,
        epilog,
        stops,
        only_xml2: true,
        query_scopes: None,
    };
    let hit_lists: Vec<Option<BlastnHitList>> = queries.iter().map(|_| None).collect();
    write_blastn_xml(&xml, document, &hit_lists, out, None)
}

/// The `-subject` XML2 or JSON document (`kind`) of queries titled `titles`, written by `calls`
/// with the document's state handed from call to call (`Xml2Carry`): the output and the error.
fn xml2(kind: Xml2Kind, titles: &[&[u8]], calls: Calls<'_>) -> (Vec<u8>, Option<String>) {
    let queries = records("Query", titles);
    let mut carry = Xml2Carry::default();
    let mut out = Vec::new();
    for (range, epilog, stops) in calls.plan(queries.len()) {
        let params = Xml2Run {
            target: Xml2Target::Subjects(1),
            subject_names: None,
            gapped: true,
            cbs: 0,
            query_gencode: None,
            db_gencode: None,
            program_type: None,
        };
        let mut document = Xml2Writer::for_run(kind, params, None, None);
        document.resume(std::mem::take(&mut carry));
        let result = call(&mut document, &queries[range], (epilog, stops), &mut out);
        carry = document.take_carry();
        if let Err(error) = result {
            return (out, Some(error.to_string()));
        }
    }
    (out, None)
}

fn tail(kind: Xml2Kind) -> &'static [u8] {
    match kind {
        Xml2Kind::Single15 => TAIL15,
        _ => TAIL16,
    }
}

/// What an unrefused run of `titles` writes before its document's end: the reports that NCBI
/// writes before it makes a later report.
fn before_end(kind: Xml2Kind, titles: &[&[u8]]) -> Vec<u8> {
    let (mut bytes, error) = xml2(kind, titles, Calls::One);
    assert_eq!(error, None, "{titles:?}");
    assert!(bytes.ends_with(tail(kind)), "{kind:?}");
    bytes.truncate(bytes.len() - tail(kind).len());
    bytes
}

const KINDS: [Xml2Kind; 2] = [Xml2Kind::Single16, Xml2Kind::Single15];

#[test]
fn every_layout_of_report_calls_writes_the_document_of_one() {
    let titles: [&[u8]; 4] = [b"qa", b"qb", b"qc", b"qd"];
    for kind in KINDS {
        let (one, error) = xml2(kind, &titles, Calls::One);
        assert_eq!(error, None);
        assert!(one.starts_with(&prolog_bytes(kind)) && one.ends_with(tail(kind)));
        for sizes in [&[1, 1, 1, 1][..], &[2, 2], &[3, 1], &[1, 3], &[4]] {
            assert_eq!(
                xml2(kind, &titles, Calls::Batches(sizes)),
                (one.clone(), None),
                "{kind:?} {sizes:?}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2403
// A `-subject` report is made when the result of another query arrives, or in the epilog
// (blast_format.cpp:2204-2230): the refusal of its title comes there, after the reports before
// it, in every layout of the calls.
#[test]
fn a_refused_report_comes_where_ncbi_makes_it_after_the_reports_before_it() {
    for kind in KINDS {
        // In the middle: qc's report is made when qd's result arrives.
        let titles: [&[u8]; 4] = [b"qa", b"qb", PUNCTUATION, b"qd"];
        let expected = before_end(kind, &[b"qa", b"qb"]);
        for calls in [
            Calls::One,
            Calls::Batches(&[1, 1, 1, 1]),
            Calls::Batches(&[2, 2]),
            Calls::Batches(&[3, 1]),
            Calls::Batches(&[1, 3]),
            Calls::Batches(&[2, 1, 1]),
        ] {
            let (bytes, error) = xml2(kind, &titles, calls);
            assert_eq!(bytes, expected, "{kind:?} {calls:?}");
            assert!(error.is_some_and(|error| error.contains(REFUSAL)));
        }
        // The last: made in the epilog.
        let titles: [&[u8]; 2] = [b"qa", PUNCTUATION];
        let expected = before_end(kind, &[b"qa"]);
        for calls in [Calls::One, Calls::Batches(&[1, 1]), Calls::Batches(&[2])] {
            let (bytes, error) = xml2(kind, &titles, calls);
            assert_eq!(bytes, expected, "{kind:?} {calls:?}");
            assert!(error.is_some_and(|error| error.contains(REFUSAL)));
        }
        // The first: NCBI's start waits in its buffer, which the crash loses, so nothing.
        for titles in [&[PUNCTUATION, b"qb"][..], &[PUNCTUATION]] {
            for calls in [Calls::One, Calls::Batches(&[1; 2][..titles.len()])] {
                let (bytes, error) = xml2(kind, titles, calls);
                assert!(
                    bytes.is_empty(),
                    "{kind:?} {calls:?}: {} bytes",
                    bytes.len()
                );
                assert!(error.is_some_and(|error| error.contains(REFUSAL)));
            }
        }
        // Two refused reports: the first one's refusal comes, where it is made.
        let titles: [&[u8]; 3] = [b"qa", PUNCTUATION, PUNCTUATION];
        let expected = before_end(kind, &[b"qa"]);
        for calls in [
            Calls::One,
            Calls::Batches(&[1, 1, 1]),
            Calls::Batches(&[2, 1]),
        ] {
            let (bytes, error) = xml2(kind, &titles, calls);
            assert_eq!(bytes, expected, "{kind:?} {calls:?}");
            assert!(error.is_some_and(|error| error.contains(REFUSAL)));
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:318-327
// A run that an error of NCBI's stops never reaches the epilog: the report of its last query is
// never made, nor refused, and NCBI's exit flushes what it wrote (the start, from `PrintProlog`).
#[test]
fn a_run_that_an_error_stops_never_makes_its_last_report_nor_refuses_it() {
    for kind in KINDS {
        let (bytes, error) = xml2(kind, &[b"qa", PUNCTUATION], Calls::Stops);
        assert_eq!((bytes, error), (before_end(kind, &[b"qa"]), None));
        let (bytes, error) = xml2(kind, &[PUNCTUATION], Calls::Stops);
        assert_eq!((bytes, error), (prolog_bytes(kind), None), "{kind:?}");
        // A refused report before the last is made at the last one's arrival.
        let (bytes, error) = xml2(kind, &[PUNCTUATION, b"qb"], Calls::Stops);
        assert!(bytes.is_empty() && error.is_some_and(|error| error.contains(REFUSAL)));
    }
}

/// A report call that only holds its report writes nothing, not even the start of the document:
/// the next call writes the start with the report (`StartFirst`, `start_not_written`).
#[test]
fn the_start_comes_with_the_first_written_report() {
    for kind in KINDS {
        let queries = records("Query", &[b"qa", b"qb"]);
        let params = || Xml2Run {
            target: Xml2Target::Subjects(1),
            subject_names: None,
            gapped: true,
            cbs: 0,
            query_gencode: None,
            db_gencode: None,
            program_type: None,
        };
        let mut out = Vec::new();
        let mut document = Xml2Writer::for_run(kind, params(), None, None);
        call(&mut document, &queries[..1], (false, false), &mut out).expect("qa's call");
        assert!(out.is_empty(), "{kind:?}: {} bytes", out.len());
        let carry = document.take_carry();
        let mut document = Xml2Writer::for_run(kind, params(), None, None);
        document.resume(carry);
        call(&mut document, &queries[1..], (false, false), &mut out).expect("qb's call");
        assert_eq!(out, before_end(kind, &[b"qa", b"qb"])[..out.len()]);
        assert!(out.starts_with(&prolog_bytes(kind)) && out.len() > prolog_bytes(kind).len());
    }
}

/// BLAST XML (outfmt 5) of queries titled `titles` written by `calls`: each call's writer numbers
/// its queries after those of the calls before (`BlastXmlWriter::new`).
fn xml5(titles: &[&[u8]], calls: Calls<'_>) -> (Vec<u8>, Option<String>) {
    let queries = records("Query", titles);
    let mut out = Vec::new();
    for (range, epilog, stops) in calls.plan(queries.len()) {
        let mut document = BlastXmlWriter::new(range.start);
        if let Err(error) = call(&mut document, &queries[range], (epilog, stops), &mut out) {
            return (out, Some(error.to_string()));
        }
    }
    (out, None)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:724-766
// BLAST XML writes each report at once: a call refuses its reports before its first byte, after the
// reports of the calls before (the Owner's Q-X2c, 2026-11-16), without the document's end.
#[test]
fn blast_xml_refuses_a_batch_before_its_first_byte() {
    let (one, error) = xml5(&[b"qa"], Calls::One);
    assert_eq!(error, None);
    let end = b"\n</BlastOutput_iterations>\n</BlastOutput>\n\n";
    assert!(one.ends_with(end));
    let titles: [&[u8]; 2] = [b"qa", PUNCTUATION];
    let (bytes, error) = xml5(&titles, Calls::Batches(&[1, 1]));
    assert_eq!(bytes, one[..one.len() - end.len()]);
    assert!(error.is_some_and(|error| error.contains(REFUSAL)));
    for calls in [Calls::One, Calls::Batches(&[2])] {
        let (bytes, error) = xml5(&titles, calls);
        assert!(bytes.is_empty(), "{calls:?}");
        assert!(error.is_some_and(|error| error.contains(REFUSAL)));
    }
}
