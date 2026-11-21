//! `HeldFormats`: what a run that held its output writes when it ends (`end_document_with`:
//! the XML, XML2 and JSON documents), and the side files of outfmt 13 and 14 that follow it
//! (`side_files`).
use super::*;
use crate::cli::NativeError;
use crate::report::blast_xml2::SideFiles;
use crate::seqdb::SeqDbError;
use std::sync::{Arc, Mutex};

/// A stream that the test reads after the hold wrote to it.
#[derive(Clone, Default)]
struct Shared(Arc<Mutex<Vec<u8>>>);

impl Write for Shared {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl Shared {
    fn text(&self) -> String {
        String::from_utf8(self.0.lock().unwrap().clone()).unwrap()
    }
}

/// Side files that record each file as `<n>:<bytes>;`.
struct Recorded(Shared);

impl SideFiles for Recorded {
    fn create(&mut self, _n: u32) -> std::io::Result<Box<dyn Write + '_>> {
        Err(std::io::Error::other("a run writes through write_file"))
    }

    fn listed(&self, n: u32) -> Vec<u8> {
        format!("f{n}").into_bytes()
    }

    fn write_file(&mut self, n: u32, report: &[u8]) -> Result<()> {
        write!(self.0, "{n}:{};", String::from_utf8_lossy(report))?;
        Ok(())
    }
}

/// A run that holds two streams and its side files, wrote `a1 b1 a2` and two side files, and
/// has not ended.
fn held_run() -> (HeldFormats<Shared>, [Shared; 3]) {
    let (a, b, files) = (Shared::default(), Shared::default(), Shared::default());
    let hold = HeldFormats::new();
    let (mut first, mut second) = (hold.writer(a.clone()), hold.writer(b.clone()));
    hold.hold();
    let mut side = hold.side_files(Box::new(Recorded(files.clone())));
    first.write_all(b"a1 ").unwrap();
    side.write_file(1, b"one").unwrap();
    second.write_all(b"b1 ").unwrap();
    first.write_all(b"a2 ").unwrap();
    side.write_file(2, b"two").unwrap();
    assert_eq!(side.listed(2), b"f2");
    for stream in [&a, &b, &files] {
        assert_eq!(stream.text(), "", "nothing passes a hold");
    }
    (hold, [a, b, files])
}

fn written(streams: &[Shared; 3]) -> [String; 3] {
    [streams[0].text(), streams[1].text(), streams[2].text()]
}

const ALL: [&str; 3] = ["a1 a2 ", "b1 ", "1:one;2:two;"];
const NOTHING: [&str; 3] = ["", "", ""];

fn ncbis_error() -> anyhow::Error {
    NativeError {
        exit: 3,
        message: "BLAST engine error: Warning: Sequence contains no data\n".into(),
    }
    .into()
}

/// A refusal of the taxonomy database as the XML2 and JSON documents meet it (`XmlRefusal`,
/// made of the reader's `TaxonomyRefusal`, which outfmt 5's prefetch returns as it is).
fn document_refusal() -> anyhow::Error {
    super::super::blast_xml::XmlRefusal::of_taxonomy(
        crate::algorithm::common::database::taxdb_error_for_tests(SeqDbError::Unsupported {
            path: "taxdb.btd".into(),
            message: "the taxonomy record of taxid 9606 reaches past the mapped end".into(),
        }),
    )
}

#[test]
fn a_release_writes_the_side_files_and_every_stream_in_order() {
    let (hold, streams) = held_run();
    hold.release(&|_| false).unwrap();
    assert_eq!(written(&streams), ALL);
}

/// Review 3 N5 of line RG T7: a format that failed alone (`FormatRefusal`) gets nothing of
/// what was held, so that its writer (a LOSAT Web stream) never had a byte of it.
#[test]
fn a_release_drops_the_writes_of_a_stream_that_failed_alone() {
    let (hold, streams) = held_run();
    hold.release(&|sink| sink == 0).unwrap();
    assert_eq!(written(&streams), ["", "b1 ", "1:one;2:two;"]);
}

/// NCBI's error of a side file that it cannot open (`CArgException`, exit status 1): no
/// error archive follows it (`error_lines`), and NCBI wrote the bytes before it all the same.
fn side_file_error() -> anyhow::Error {
    crate::report::blast_xml2::framing::side_file_open_error()
}

#[test]
fn the_documents_hold_drops_after_their_refusal_and_writes_after_any_other_error() {
    let (hold, streams) = held_run();
    hold.end_document_with(&document_refusal(), &|_| false)
        .unwrap();
    assert_eq!(written(&streams), NOTHING);
    // The taxonomy database's refusal of the prefetch (outfmt 5) is one of them.
    let (hold, streams) = held_run();
    let prefetch =
        crate::algorithm::common::database::taxdb_error_for_tests(SeqDbError::Unsupported {
            path: "taxdb.btd".into(),
            message: "the taxonomy record of taxid 9606 reaches past the mapped end".into(),
        });
    hold.end_document_with(&prefetch, &|_| false).unwrap();
    assert_eq!(written(&streams), NOTHING);
    for other in [
        ncbis_error(),
        side_file_error(),
        anyhow::anyhow!("taxdb.btd: Input/output error; this database is not supported"),
    ] {
        let (hold, streams) = held_run();
        hold.end_document_with(&other, &|_| false).unwrap();
        assert_eq!(written(&streams), ALL, "{other}");
    }
}

/// Side files whose every file cannot be opened.
struct Blocked;

impl SideFiles for Blocked {
    fn create(&mut self, _n: u32) -> std::io::Result<Box<dyn Write + '_>> {
        Err(std::io::Error::other("blocked"))
    }

    fn listed(&self, _n: u32) -> Vec<u8> {
        Vec::new()
    }
}

/// The side files of one document cannot be opened at the release (they could when their
/// reports were made): NCBI's error ends the run, and neither the streams nor the side files
/// of the document after it get what was held, then or when the run's writers are dropped.
#[test]
fn a_side_file_error_at_the_release_writes_nothing_after_it() {
    let (a, files) = (Shared::default(), Shared::default());
    let hold = HeldFormats::new();
    let mut stream = hold.writer(a.clone());
    hold.hold();
    let mut first = hold.side_files(Box::new(Recorded(files.clone())));
    let mut blocked = hold.side_files(Box::new(Blocked));
    let mut last = hold.side_files(Box::new(Recorded(files.clone())));
    stream.write_all(b"a1 ").unwrap();
    first.write_file(1, b"one").unwrap();
    blocked.write_file(1, b"x").unwrap();
    last.write_file(1, b"late").unwrap();
    let error = hold.release(&|_| false).unwrap_err();
    let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
    assert_eq!(
        native.message,
        b"Command line argument error: Cannot open output file\n"
    );
    drop((first, blocked, last, stream, hold));
    assert_eq!([a.text(), files.text()], ["", "1:one;"].map(String::from));
}

/// A stream that fails every write (a full device).
struct Full;

impl Write for Full {
    fn write(&mut self, _bytes: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("no space left"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Two formats' streams, one behind a buffer (`std::io::BufWriter`) and one unbuffered, both
/// read from one pipe (`log`), given `report 1` and a flush, the text `text 1`, `report 2`
/// without a flush and `text 2`, through a hold that holds (`held`) or not: the pipe after the
/// last text and at the end.
fn through_one_pipe(held: bool) -> [String; 2] {
    let log = Shared::default();
    let hold: HeldFormats<Box<dyn Write + Send>> = HeldFormats::new();
    let mut buffered = hold.writer(Box::new(std::io::BufWriter::new(log.clone())));
    let mut unbuffered = hold.writer(Box::new(log.clone()));
    if held {
        hold.hold();
    }
    buffered.write_all(b"report 1 ").unwrap();
    buffered.flush().unwrap();
    unbuffered.write_all(b"text 1 ").unwrap();
    buffered.write_all(b"report 2 ").unwrap();
    unbuffered.write_all(b"text 2 ").unwrap();
    if held {
        hold.release(&|_| false).unwrap();
    }
    let after_texts = log.text();
    drop((buffered, unbuffered, hold));
    [after_texts, log.text()]
}

/// The release makes the writes and flushes of the run in their order (line RG T2 review 5
/// D1): a buffered stream's bytes flushed before a text on another stream reach the pipe
/// before it, those not flushed stay in its buffer, as without the hold.
#[test]
fn a_release_writes_and_flushes_as_the_run_without_the_hold() {
    let unheld = through_one_pipe(false);
    assert_eq!(
        unheld,
        [
            "report 1 text 1 text 2 ",
            "report 1 text 1 text 2 report 2 "
        ]
    );
    assert_eq!(through_one_pipe(true), unheld);
}

/// A write of the release that fails stops it there: nothing held after it is written, then
/// or later, and the run ends with it in place of its own error (`end_document_with`: the
/// writer's error, as every program's library call ends; the command line's stream ends as
/// NCBI's, `ReportStream::finish`).
#[test]
fn a_failing_write_of_the_release_ends_the_run() {
    let log = Shared::default();
    let hold: HeldFormats<Box<dyn Write + Send>> = HeldFormats::new();
    let mut full = hold.writer(Box::new(Full));
    let mut other = hold.writer(Box::new(log.clone()));
    hold.hold();
    other.write_all(b"text 0 ").unwrap();
    full.write_all(b"report 1 ").unwrap();
    full.flush().unwrap();
    other.write_all(b"text 1 ").unwrap();
    let error = hold
        .end_document_with(&ncbis_error(), &|_| false)
        .unwrap_err();
    assert!(
        error.downcast_ref::<std::io::Error>().is_some(),
        "{error:?}"
    );
    other.write_all(b"text 2 ").unwrap();
    assert_eq!(log.text(), "text 0 ");
}

#[test]
fn after_a_discard_nothing_more_is_written() {
    let (hold, streams) = held_run();
    let mut late = hold.writer(streams[0].clone());
    let mut side = hold.side_files(Box::new(Recorded(streams[2].clone())));
    hold.discard();
    late.write_all(b"late").unwrap();
    side.write_file(3, b"three").unwrap();
    hold.release(&|_| false).unwrap();
    assert_eq!(written(&streams), NOTHING);
}
