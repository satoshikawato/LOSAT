//! The hold of the side files of outfmt 13 and 14 (`HeldSideFiles`, `HeldOutput`): nothing on
//! the file system while the run holds its output, the files in report order at `release`, none
//! after `discard`, and the files made so far when the run ends with another error.
use super::*;
use crate::cli::NativeError;
use crate::report::blast_xml2::{FsSideFiles, Xml2Kind};
use std::path::{Path, PathBuf};

fn work(stem: &str) -> PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat_held_side_{stem}_{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::remove_dir_all(&dir).ok();
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A held output over `<dir>/o.out` and its side files `<dir>/o.out_<n>.xml`.
fn run(dir: &Path) -> (HeldOutput, HeldSideFiles) {
    let master = std::fs::File::create(dir.join("o.out")).unwrap();
    let side = HeldSideFiles::new(Box::new(FsSideFiles::new(
        &dir.join("o.out"),
        Xml2Kind::Files14,
    )));
    (HeldOutput::new(Box::new(master), Some(side.clone())), side)
}

/// A report as a run writes it (`SideFiles::write_file`).
fn report(side: &mut HeldSideFiles, n: u32, bytes: &[u8]) {
    side.write_file(n, bytes).unwrap();
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// One log of the output (`o:`) and of standard error (`e:`) in the order written.
#[derive(Clone, Default)]
struct Log(Rc<RefCell<Vec<String>>>, &'static str);

impl Write for Log {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if !bytes.is_empty() {
            self.0
                .borrow_mut()
                .push(format!("{}{}", self.1, String::from_utf8_lossy(bytes)));
        }
        Ok(bytes.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn a_held_output_writes_the_standard_error_texts_where_they_came() {
    let log = Rc::new(RefCell::new(Vec::new()));
    let held = HeldOutput::new(Box::new(Log(Rc::clone(&log), "o:")), None);
    let mut stderr = held.stderr(Box::new(Log(Rc::clone(&log), "e:")));
    let mut output = held.clone();
    // Not held: each goes at once.
    output.write_all(b"head").unwrap();
    stderr.write_all(b"w1").unwrap();
    held.hold();
    output.write_all(b"hit1").unwrap();
    stderr.write_all(b"t2").unwrap();
    output.write_all(b"hit2").unwrap();
    assert_eq!(*log.borrow(), ["o:head", "e:w1"]);
    held.release().unwrap();
    assert_eq!(
        *log.borrow(),
        ["o:head", "e:w1", "o:hit1", "e:t2", "o:hit2"]
    );
    // A discarded output keeps its standard error texts and drops the output.
    let log = Rc::new(RefCell::new(Vec::new()));
    let held = HeldOutput::new(Box::new(Log(Rc::clone(&log), "o:")), None);
    let mut stderr = held.stderr(Box::new(Log(Rc::clone(&log), "e:")));
    let mut output = held.clone();
    held.hold();
    output.write_all(b"hit1").unwrap();
    stderr.write_all(b"t1").unwrap();
    held.discard();
    stderr.write_all(b"refusal").unwrap();
    assert_eq!(*log.borrow(), ["e:t1", "e:refusal"]);
}

#[test]
fn without_a_hold_the_files_are_written_at_once() {
    let dir = work("nohold");
    let (held, mut side) = run(&dir);
    report(&mut side, 1, b"one");
    assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"one");
    held.release().unwrap();
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn a_hold_keeps_the_files_in_memory_until_the_release_and_writes_them_in_order() {
    let dir = work("hold");
    let (held, mut side) = run(&dir);
    held.hold();
    report(&mut side, 1, b"one");
    report(&mut side, 2, b"two");
    let mut master = held.clone();
    master.write_all(b"master").unwrap();
    assert_eq!(names(&dir), ["o.out"]);
    assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"");
    held.release().unwrap();
    assert_eq!(names(&dir), ["o.out", "o.out_1.xml", "o.out_2.xml"]);
    assert_eq!(std::fs::read(dir.join("o.out_2.xml")).unwrap(), b"two");
    assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"master");
    // after the release the files are written at once again
    report(&mut side, 3, b"three");
    assert_eq!(std::fs::read(dir.join("o.out_3.xml")).unwrap(), b"three");
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn a_discard_leaves_no_side_file_and_an_empty_master() {
    let dir = work("discard");
    let (held, mut side) = run(&dir);
    held.hold();
    report(&mut side, 1, b"one");
    let mut master = held.clone();
    master.write_all(b"master").unwrap();
    held.discard();
    report(&mut side, 2, b"two");
    drop((held, side, master));
    assert_eq!(names(&dir), ["o.out"]);
    assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"");
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn any_other_end_of_the_run_writes_the_files_made_so_far() {
    let dir = work("dropped");
    let (held, mut side) = run(&dir);
    held.hold();
    report(&mut side, 1, b"one");
    drop((held, side));
    assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"one");
    std::fs::remove_dir_all(dir).ok();
}

/// DECISIONS 35 S1: the run that holds its side files opens each one when its report is made, as
/// NCBI does: a file that cannot be opened stops the run there (the earlier reports are written
/// by the drop, the master stays empty), and the check leaves no file of its own.
#[test]
fn a_side_file_that_cannot_be_opened_stops_at_its_report_while_holding() {
    let dir = work("blocked");
    std::fs::create_dir(dir.join("o.out_2.xml")).unwrap();
    let (held, mut side) = run(&dir);
    held.hold();
    report(&mut side, 1, b"one");
    // nothing but the master file yet, and the check of the first report left no mark
    assert_eq!(names(&dir), ["o.out", "o.out_2.xml"]);
    let error = side.create(2).err().expect("a directory cannot be opened");
    assert!(error.raw_os_error().is_some(), "{error:?}");
    // The run's route: NCBI's error of the file, and nothing of the report is kept.
    let error = side.write_file(2, b"two").unwrap_err();
    let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
    assert_eq!(
        native.message,
        b"Command line argument error: Cannot open output file\n"
    );
    let mut master = held.clone();
    master.write_all(b"master").unwrap();
    drop((held, side, master));
    // the file before the blocked one stays; the master has nothing of the held bytes of the
    // failed report
    assert_eq!(names(&dir), ["o.out", "o.out_1.xml", "o.out_2.xml"]);
    assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"one");
    std::fs::remove_dir_all(dir).ok();
}

/// The check opens an existing file for writing without truncating it, and a refusal after it
/// (`discard`) leaves the file as it was.
#[test]
fn the_check_of_an_existing_side_file_does_not_truncate_it() {
    let dir = work("existing");
    std::fs::write(dir.join("o.out_1.xml"), b"old").unwrap();
    let (held, mut side) = run(&dir);
    held.hold();
    report(&mut side, 1, b"new");
    assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"old");
    held.discard();
    drop((held, side));
    assert_eq!(std::fs::read(dir.join("o.out_1.xml")).unwrap(), b"old");
    std::fs::remove_dir_all(dir).ok();
}

#[cfg(unix)]
#[test]
fn a_read_only_side_file_stops_at_its_report_while_holding() {
    use std::os::unix::fs::PermissionsExt;
    let dir = work("readonly");
    std::fs::write(dir.join("o.out_1.xml"), b"old").unwrap();
    std::fs::set_permissions(
        dir.join("o.out_1.xml"),
        std::fs::Permissions::from_mode(0o444),
    )
    .unwrap();
    let (held, mut side) = run(&dir);
    held.hold();
    // a process with the right to write anyway (root) opens it, as it opens it for NCBI
    let writable = std::fs::OpenOptions::new()
        .write(true)
        .open(dir.join("o.out_1.xml"))
        .is_ok();
    assert_eq!(side.create(1).is_ok(), writable);
    drop((held, side));
    std::fs::remove_dir_all(dir).ok();
}

/// A write that fails after the open (a link to `/dev/full`) is the write error of every
/// program, not NCBI's "Cannot open output file" (DECISIONS 35 S2).
#[cfg(target_os = "linux")]
#[test]
fn a_side_file_write_error_at_the_release_is_the_write_error() {
    let dir = work("full");
    std::os::unix::fs::symlink("/dev/full", dir.join("o.out_1.xml")).unwrap();
    let (held, mut side) = run(&dir);
    held.hold();
    report(&mut side, 1, b"one");
    let mut master = held.clone();
    master.write_all(b"master").unwrap();
    let error = held.release().unwrap_err();
    let native = error.downcast_ref::<NativeError>().expect("a NativeError");
    assert_eq!(native.exit, 6);
    assert_eq!(native.message, b"BLAST failed to write output\n");
    drop((held, side, master));
    assert_eq!(std::fs::read(dir.join("o.out")).unwrap(), b"");
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn a_side_file_that_cannot_be_created_at_the_release_is_the_open_error() {
    // The check at the report passed; a directory put in its place before the release is
    // what `File::create` then refuses.
    let dir = work("late");
    let (held, mut side) = run(&dir);
    held.hold();
    report(&mut side, 1, b"one");
    std::fs::create_dir(dir.join("o.out_1.xml")).unwrap();
    let error = held.release().unwrap_err();
    let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
    assert_eq!(native.exit, 1);
    assert_eq!(
        native.message,
        b"Command line argument error: Cannot open output file\n"
    );
    drop((held, side));
    std::fs::remove_dir_all(dir).ok();
}

/// Discarding a held single-stream document (`HeldOutput::discard` on a LOSAT refusal) drops the
/// held report bytes while still writing both queued stderr texts and queued `HeldDiag::Posts`
/// (including a batch's `803.7` replaced-byte diagnostic from `-export_search_strategy`) to stderr
/// in order.
#[test]
fn discard_emits_queued_strategy_posts_to_stderr() {
    use crate::blastinput::app::Program;
    use crate::report::archive::bioseq::BioseqInput;
    use crate::report::archive::options::OptionInputs;
    use crate::report::archive::strategy::{StrategyExport, StrategySink};
    use crate::report::archive::writer::{RequestInputs, RequestQueries};
    use crate::report::archive::{ArchiveHandle, OptionsHandle};
    use crate::report::ostream::emit_strategy_posts;
    use crate::report::serial::encoder::BufferPiece;

    #[derive(Clone, Default)]
    struct Sink(Arc<Mutex<Vec<u8>>>);
    impl Sink {
        fn record(
            &mut self,
            bytes: &[u8],
            posts: &[(u64, Vec<u8>)],
            diagnostics: &mut dyn Write,
            archive: Option<&ArchiveHandle>,
        ) -> io::Result<()> {
            let _ = emit_strategy_posts(posts, None, diagnostics, archive)?;
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(())
        }
    }
    impl StrategySink for Sink {
        fn write_request(
            &mut self,
            bytes: &[u8],
            _pieces: &[BufferPiece],
            posts: &[(u64, Vec<u8>)],
            diagnostics: &mut dyn Write,
            archive: Option<&ArchiveHandle>,
        ) -> io::Result<()> {
            self.record(bytes, posts, diagnostics, archive)
        }
        fn write_partial(
            &mut self,
            bytes: &[u8],
            posts: &[(u64, Vec<u8>)],
            diag: &mut dyn Write,
            archive: Option<&ArchiveHandle>,
        ) -> io::Result<()> {
            self.record(bytes, posts, diag, archive)
        }
    }

    let log = Rc::new(RefCell::new(Vec::<String>::new()));
    let held = HeldOutput::new(Box::new(Log(Rc::clone(&log), "o:")), None);
    let mut err = held.stderr(Box::new(Log(Rc::clone(&log), "e:")));
    held.hold();

    let sink = Sink::default();
    let export = StrategyExport::new(
        Box::new(sink.clone()),
        &OptionInputs::new(Program::Blastx, "blastx"),
        None,
        OptionsHandle::default(),
    );
    let bad = [BioseqInput {
        id: "Query_1",
        ids: None,
        seq_loc: None,
        defline: std::borrow::Cow::Borrowed(b">q caf\xe9"),
        title: b"q caf\xe9",
        letters: b"ACGTACGTACGT",
        protein: false,
    }];
    let req = RequestInputs {
        queries: RequestQueries::Bioseqs(&bad),
        subjects: &[],
        required: None,
        lcase_masks: &[],
    };
    export.write_request(0..1, &req).unwrap();
    let posts = export.take_all_posts();

    let mut report = held.clone();
    report.write_all(b"<BlastOutput>").unwrap();
    err.write_all(b"reader warning\n").unwrap();
    assert!(held.push_posts(posts));
    err.write_all(b"later warning\n").unwrap();

    held.discard();
    drop((report, err, held));

    let events = log.borrow();
    // Held report bytes were dropped; stderr kept all warnings and the 803.7 post in order.
    assert!(
        events.iter().all(|e| e.starts_with("e:")),
        "unexpected report output after discard: {events:?}"
    );
    let joined: String = events.iter().map(|e| &e[2..]).collect();
    let w1 = joined.find("reader warning\n").expect("reader warning");
    let p803 = joined
        .find("Error: (803.7) [blastx]")
        .expect("803.7 strategy post");
    let w2 = joined.find("later warning\n").expect("later warning");
    assert!(w1 < p803 && p803 < w2, "unexpected stderr order:\n{joined}");
    assert!(
        !sink.0.lock().unwrap().is_empty(),
        "strategy request should be written when Posts are emitted on discard"
    );
}
