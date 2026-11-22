use super::*;

/// A device that takes nothing (`/dev/full`), counting the bytes offered to it.
struct Full(std::sync::Arc<std::sync::atomic::AtomicUsize>);

impl std::io::Write for Full {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0
            .fetch_add(buf.len(), std::sync::atomic::Ordering::Relaxed);
        Err(std::io::Error::other("full"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn the_organism_reports_file_buffer_writes_at_its_8191st_byte() {
    // probe `T/t2/fix1` `f1.d5.*` and `stdout_probe.sh`: NCBI exits 0 for a report that its
    // buffer holds to the end, 6 where the buffer writes.
    use std::io::Write;
    let offered = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let mut buffer = ReportStream::file(Box::new(Full(offered.clone())), "blastn", Site::Report);
    buffer.write_all(&[b'x'; 8190]).unwrap();
    buffer.flush().unwrap();
    assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 0);
    // The 8191st byte makes NCBI's write, which LOSAT makes before anything else is written
    // (`Core::sync`): all 8,191 bytes at once.
    buffer.write_all(b"x").unwrap();
    assert_eq!(lock(&buffer.0).ncbi_written(), 8191);
    assert!(lock(&buffer.0).sync().is_err());
    assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 8191);
    // What a buffer holds at the end is written, a failure unreported.
    let quiet = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let mut buffer = ReportStream::file(Box::new(Full(quiet.clone())), "blastn", Site::Report);
    buffer.write_all(b"held").unwrap();
    drop(buffer);
    assert_eq!(quiet.load(std::sync::atomic::Ordering::Relaxed), 4);
}

/// A device that takes everything, keeping the size of each write.
struct Writes(std::sync::Arc<std::sync::Mutex<Vec<usize>>>);

impl std::io::Write for Writes {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        lock(&self.0).push(buf.len());
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn the_organism_reports_standard_output_writes_where_ncbis_writes() {
    // `std::cout` on a file of `st_blksize` 4096 (observed, glibc 2.35 and 2.39,
    // `tests/fixtures/stdio_buffer/cout_file.tsv`): the writes end at multiples of the block
    // and the rest stays held; the first insertion has no buffer yet; a flush empties the
    // buffer.
    use std::io::Write;
    let run = |pieces: &[usize]| {
        let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut after = Vec::new();
        let mut buffer = ReportStream::new(
            Box::new(Writes(writes.clone())),
            Some(NcbiBuffer::new(NcbiStdio::Cout {
                block: 4096,
                lines: false,
            })),
            false,
            "blastn",
            Site::Report,
        );
        for &piece in pieces {
            match piece {
                0 => lock(&buffer.0).flush().unwrap(),
                n => buffer.write_all(&vec![b'x'; n]).unwrap(),
            }
            after.push(lock(&buffer.0).ncbi_written() as usize);
        }
        drop(buffer);
        let total = lock(&writes).iter().sum::<usize>();
        (after, total)
    };
    // `4000 200 4000`: 4,096 written at the second, 8,192 at the third, 8 held to the end.
    assert_eq!(run(&[4000, 200, 4000]), (vec![0, 4096, 8192], 8200));
    assert_eq!(run(&[4096, 10]), (vec![4096, 4096], 4106));
    assert_eq!(run(&[4095, 1, 10]), (vec![0, 0, 4096], 4106));
    assert_eq!(run(&[8192, 1]), (vec![8192, 8192], 8193));
    assert_eq!(run(&[100, 4096]), (vec![0, 4096], 4196));
    // `100 f 4096 1` (0: the flush of a write to `cerr`, which is tied to `cout`).
    assert_eq!(run(&[100, 0, 4096, 1]), (vec![0, 100, 100, 4196], 4197));
    assert_eq!(run(&[8000, 300]), (vec![4096, 8192], 8300));
}

/// The bytes that `buffer` has written to its device after each step of `steps`
/// (`Ok(bytes)`: one insertion; `Err(())`: a flush, as before a write to `cerr`), and in all.
fn written_after(stream: NcbiStdio, steps: &[Result<Vec<u8>, ()>]) -> (Vec<usize>, usize, usize) {
    use std::io::Write;
    let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut after = Vec::new();
    let mut buffer = ReportStream::new(
        Box::new(Writes(writes.clone())),
        Some(NcbiBuffer::new(stream)),
        false,
        "blastn",
        Site::Report,
    );
    for step in steps {
        match step {
            Ok(bytes) => buffer.write_all(bytes).unwrap(),
            Err(()) => lock(&buffer.0).flush().unwrap(),
        }
        after.push(lock(&buffer.0).ncbi_written() as usize);
    }
    drop(buffer);
    let writes = lock(&writes);
    (after, writes.iter().sum(), writes.len())
}

#[test]
fn a_terminals_standard_output_writes_at_its_newlines_as_ncbis_does() {
    // `std::cout` on a pty of `st_blksize` 1024 (observed, glibc 2.35 and 2.39,
    // `tests/fixtures/stdio_buffer/cout_tty.tsv`): `N` is an insertion of N bytes, `N:K` with
    // a newline at K, `f` a flush.
    let tty = NcbiStdio::Cout {
        block: 1024,
        lines: true,
    };
    let run = |steps: &str| {
        let steps = steps
            .split(' ')
            .map(|step| match step {
                "f" => Err(()),
                _ => {
                    let mut parts = step.split(':').map(|n| n.parse::<usize>().unwrap());
                    let mut bytes = vec![b'x'; parts.next().unwrap()];
                    parts.for_each(|at| bytes[at] = b'\n');
                    Ok(bytes)
                }
            })
            .collect::<Vec<_>>();
        written_after(tty, &steps).0
    };
    assert_eq!(run("100"), [0]);
    assert_eq!(run("100 50:49"), [0, 150]);
    assert_eq!(run("100 50:10"), [0, 111]);
    assert_eq!(run("3000"), [2048]);
    assert_eq!(run("3000:2500"), [2501]);
    assert_eq!(run("1000 30"), [0, 1024]);
    assert_eq!(run("1000 24 1"), [0, 0, 1024]);
    assert_eq!(run("1000 30:5"), [0, 1024]);
    assert_eq!(run("1000 30:27"), [0, 1028]);
    assert_eq!(run("500 2000:100"), [0, 2048]);
    assert_eq!(run("500 2000:100:1900"), [0, 2401]);
    assert_eq!(run("f 100 f 50"), [0, 0, 100, 100]);
}

#[test]
fn a_streamed_insertion_is_written_as_the_whole_insertion() {
    // `insert_streamed` writes pieces; every model of NCBI's buffers must write what it writes
    // for the whole insertion as one write, before it and after it. A row of 70,000 bytes
    // ends 99 bytes after a newline: no piece may end there (a terminal would keep only 99).
    let row = [&[b'r'; 69_900][..], b"\n", &[b'r'; 99]].concat();
    let segments: &[(&[u8], u64)] = &[
        (b"\n  Organism Report\nrow\n", 200_000),
        (&row, 3_000),
        (b"|x\ny", 5_000),
        (b"z", 70_000),
        (b"tail\nend", 0),
    ];
    let mut whole = Vec::new();
    for (bytes, spaces) in segments {
        whole.extend_from_slice(bytes);
        whole.resize(whole.len() + *spaces as usize, b' ');
    }
    let streamed = |out: &mut dyn std::io::Write| {
        insert_streamed(out, whole.len() as u64, |insertion| {
            for (bytes, spaces) in segments {
                insertion.bytes(bytes)?;
                insertion.spaces(*spaces)?;
            }
            Ok(())
        })
        .unwrap()
    };
    let mut flat = Vec::new();
    streamed(&mut flat);
    assert_eq!(flat, whole);
    for stream in [
        NcbiStdio::File,
        NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
        NcbiStdio::Cout {
            block: 1024,
            lines: true,
        },
        NcbiStdio::Cout {
            block: 100,
            lines: false,
        },
    ] {
        for before in [0, 1, 1000, 4095, 4096, 8000] {
            for after in [&b"ab"[..], b"ab\ncd", &[b'y'; 9000]] {
                let lead = vec![b'w'; before];
                // The whole as one insertion, between the same insertions.
                let one = written_after(
                    stream,
                    &[Ok(lead.clone()), Ok(whole.clone()), Ok(after.to_vec())],
                );
                // The whole streamed.
                let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
                let mut buffer = ReportStream::new(
                    Box::new(Writes(writes.clone())),
                    Some(NcbiBuffer::new(stream)),
                    false,
                    "blastn",
                    Site::Report,
                );
                let mut pieces = Vec::new();
                std::io::Write::write_all(&mut buffer, &lead).unwrap();
                pieces.push(lock(&buffer.0).ncbi_written() as usize);
                streamed(&mut buffer);
                pieces.push(lock(&buffer.0).ncbi_written() as usize);
                std::io::Write::write_all(&mut buffer, after).unwrap();
                pieces.push(lock(&buffer.0).ncbi_written() as usize);
                drop(buffer);
                let total = lock(&writes).iter().sum::<usize>();
                assert_eq!((pieces, total), (one.0, one.1), "{stream:?} {before}");
            }
        }
    }
    // In pieces: the run of 200,000 spaces ends three, the rest is the last.
    let pieces = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    streamed(&mut Writes(pieces.clone()));
    assert_eq!(*lock(&pieces), [65_559, 65_536, 65_536, 151_405]);
}

/// A device that takes `0` more bytes and then fails (a file at its size limit,
/// `prlimit --fsize`: the write that reaches the limit is cut there, the next fails).
struct FailAt(usize);

impl std::io::Write for FailAt {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if self.0 == 0 {
            return Err(std::io::Error::other("File too large"));
        }
        let n = buf.len().min(self.0);
        self.0 -= n;
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn native(result: anyhow::Result<()>) -> (i32, String) {
    let error = result.expect_err("the run fails");
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .expect("NCBI's end");
    (
        native.exit,
        String::from_utf8(native.message.clone()).unwrap(),
    )
}

#[test]
fn a_tabular_report_that_fails_at_byte_n_terminates_as_ncbis() {
    use std::io::Write;
    // Held to the end: the flush of the last `CBlastTabularInfo` fails while it throws.
    let mut stream = ReportStream::file(Box::new(FailAt(10)), "blastn", Site::Tabular);
    stream.write_all(&[b'x'; 100]).unwrap();
    assert_eq!(
        native(stream.finish(Ok(()))),
        (-6, TERMINATE_TEXT.to_string())
    );
    // Where NCBI's buffer writes: the insertion fails, and every later one.
    let mut stream = ReportStream::file(Box::new(FailAt(10)), "blastn", Site::Tabular);
    stream.write_all(&[b'x'; 8000]).unwrap();
    let failed = stream.write_all(&[b'x'; 500]).unwrap_err();
    assert!(is_report_write_failure(&failed));
    assert!(stream.write_all(b"x").is_err());
    assert_eq!(native(stream.finish(Err(failed.into()))).0, -6);
}

#[test]
fn the_alignments_post_802_4_in_a_file_and_on_standard_output_at_their_flush() {
    use std::io::Write;
    let post = "Error: (802.4) [blastp] [Cannot close serializing output stream] \
                    Exception: basic_ios::clear: iostream error\nBLAST failed to write output\n";
    let mut stream = ReportStream::file(Box::new(FailAt(5)), "blastp", Site::Report);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    ostream.write_all(b"prolog\n").unwrap();
    let mut display = DisplaySeqalign::default();
    display.enter(&mut ostream).unwrap();
    ostream.write_all(b"> subject\n").unwrap();
    // `out2`'s flush is NCBI's write; LOSAT makes it before anything else is written.
    let left = display
        .leave(&mut ostream)
        .and_then(|()| lock(&stream.0).sync());
    assert!(left.is_err());
    drop(ostream);
    assert_eq!(native(stream.finish(Ok(()))), (6, post.to_string()));
    // Standard output: a failed insertion leaves the buffer empty (no post); the flush
    // at the end of the alignments posts, and the post's tie terminates.
    let mut stream = ReportStream::new(
        Box::new(FailAt(5)),
        Some(NcbiBuffer::new(NcbiStdio::Cout {
            block: 4096,
            lines: false,
        })),
        true,
        "blastp",
        Site::Alignments,
    );
    assert!(stream.write_all(&[b'x'; 5000]).is_err());
    let message = "BLAST failed to write output\n".to_string();
    assert_eq!(native(stream.finish(Ok(()))), (6, message));
    let mut stream = ReportStream::new(
        Box::new(FailAt(5)),
        Some(NcbiBuffer::new(NcbiStdio::Cout {
            block: 4096,
            lines: false,
        })),
        true,
        "blastp",
        Site::Alignments,
    );
    stream.write_all(b"> subject\n").unwrap();
    lock(&stream.0).flush().unwrap();
    assert!(lock(&stream.0).sync().is_err());
    assert_eq!(native(stream.finish(Ok(()))).0, -6);
}

/// A second stream on the device of `stream` (the test's `NcbiOstream::report` owns it).
fn stream_of(stream: &mut ReportStream) -> ReportStream {
    ReportStream(stream.0.clone())
}

#[test]
fn a_failure_after_the_formatter_or_of_sam_is_not_reported() {
    use std::io::Write;
    // Held to the end, written after `~CBlastFormat`: unreported.
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Report);
    stream.write_all(b"held").unwrap();
    assert!(stream.finish(Ok(())).is_ok());
    assert!(stream.failed());
    // SAM writes after `~CBlastFormat`.
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Silent);
    stream.write_all(&[b'x'; 9000]).unwrap();
    stream.write_all(b"more").unwrap();
    assert!(stream.finish(Ok(())).is_ok());
}

/// The posts of `stream` go to a buffer, which is returned.
fn posts_of(stream: &ReportStream) -> std::sync::Arc<std::sync::Mutex<Vec<u8>>> {
    let posts = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    lock(&stream.0).diagnostics = Some((Box::new(Writes2(posts.clone())), Cerr::new()));
    posts
}

/// A device that keeps what it is given.
struct Writes2(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

impl std::io::Write for Writes2 {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        lock(&self.0).extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn an_xml_iteration_that_fails_is_posted_once_and_the_run_goes_on_to_the_epilog() {
    use std::io::Write;
    // probe `T/t6/survey` blastn.f5.qm full2: one post per query from the failing one, then
    // the epilog's insertion throws: "BLAST failed to write output", exit 6.
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Report);
    let posts = posts_of(&stream);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    let mut cuts = Cuts::default();
    cuts.push(0, Cut::Site(Site::Report));
    cuts.push(100, Cut::Site(Site::XmlIteration));
    cuts.push(8200, Cut::Insertion);
    cuts.push(8300, Cut::Site(Site::Report));
    let report = vec![b'x'; 8300];
    let failed = cuts
        .cursor()
        .write(&mut ostream, &report, 0, report.len())
        .unwrap_err();
    assert!(is_caught_report_write_failure(&failed));
    // The next query's `Iteration` fails at its first insertion.
    let mut next = Cuts::default();
    next.push(0, Cut::Site(Site::XmlIteration));
    next.push(10, Cut::Insertion);
    let failed = next
        .cursor()
        .write(&mut ostream, &report[..20], 0, 20)
        .unwrap_err();
    assert!(is_caught_report_write_failure(&failed));
    let post = "Error: [blastn] Failed s_BlastXMLAddIteration Q(0/1\n";
    assert_eq!(
        String::from_utf8(lock(&posts).clone()).unwrap(),
        post.repeat(2)
    );
    // The epilog, outside the `try`.
    ostream.site(Site::Report).unwrap();
    let failed = ostream.write_all(b"</BlastOutput>").unwrap_err();
    assert!(is_report_write_failure(&failed) && !is_caught_report_write_failure(&failed));
    drop(ostream);
    assert_eq!(
        native(stream.finish(Err(failed.into()))),
        (6, "BLAST failed to write output\n".to_string())
    );
    // On standard output the post's tie throws: `std::terminate`, no post.
    let mut stream = ReportStream::new(
        Box::new(FailAt(0)),
        Some(NcbiBuffer::new(NcbiStdio::Cout {
            block: 4096,
            lines: false,
        })),
        true,
        "blastn",
        Site::XmlIteration,
    );
    let posts = posts_of(&stream);
    let failed = stream.write_all(&[b'x'; 5000]).unwrap_err();
    assert!(!is_caught_report_write_failure(&failed));
    assert!(lock(&posts).is_empty());
    assert_eq!(native(stream.finish(Err(failed.into()))).0, -6);
}

#[test]
fn a_report_held_to_its_epilogs_flush_fails_there() {
    // A small BLAST XML report that NCBI's file buffer holds: `endl` at the epilog flushes it
    // while the formatter lives (exit 6); without that flush the bytes would reach the device
    // after `~CBlastFormat`, unreported.
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastp", Site::Report);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    let mut cuts = Cuts::default();
    cuts.push(0, Cut::Site(Site::Report));
    cuts.push(5, Cut::Insertion);
    cuts.push(6, Cut::Flush);
    // The epilog's flush is NCBI's write, made at the stream's end at the latest.
    let written = cuts.cursor().write(&mut ostream, b"<x/>\n\n", 0, 6);
    drop(ostream);
    assert_eq!(native(stream.finish(written.map_err(Into::into))).0, 6);
}

/// One Seq-annot written as `seqalign::write_query_to` writes it: whole, `Close`'s flush, and
/// NCBI's pieces where the whole insertion fails.
fn serial_object(
    ostream: &mut NcbiOstream<'_>,
    len: usize,
    pieces: &[usize],
) -> std::io::Result<()> {
    use std::io::Write;
    ostream.site(Site::Serial { whole: true })?;
    match ostream
        .write_all(&vec![b'x'; len])
        .and_then(|()| ostream.ncbi_flush())
    {
        Err(error) if is_undecided_report_write_failure(&error) => {
            ostream.resolve_object(pieces)?
        }
        written => written?,
    }
    ostream.site(Site::Report)
}

#[test]
fn a_serial_object_that_fails_inside_ends_the_run_and_one_that_fails_at_close_goes_on() {
    let post = "Error: (802.4) [blastn] [Cannot close serializing output stream] Exception: \
                    basic_ios::clear: iostream error\n";
    let failed = "BLAST failed to write output\n";
    let pieces = [4096, 4096, 1808];
    // The file buffer holds the first piece and writes 8,192 bytes at the second: past 8,000
    // bytes that write fails inside the object (802.4 at `Close`, then `CATCH_ALL`).
    let mut stream = ReportStream::file(Box::new(FailAt(8000)), "blastn", Site::Report);
    let posts = posts_of(&stream);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    let error = serial_object(&mut ostream, 10_000, &pieces).unwrap_err();
    drop(ostream);
    assert!(lock(&posts).is_empty());
    assert_eq!(
        native(stream.finish(Err(error.into()))),
        (6, format!("{post}{failed}"))
    );
    // Up to 9,000 bytes that write succeeds and only `Close`'s flush fails: posted, the run
    // goes on, and the next object's first insertion throws.
    let mut stream = ReportStream::file(Box::new(FailAt(9000)), "blastn", Site::Report);
    let posts = posts_of(&stream);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    serial_object(&mut ostream, 10_000, &pieces).unwrap();
    assert_eq!(String::from_utf8(lock(&posts).clone()).unwrap(), post);
    let error = serial_object(&mut ostream, 50, &[50]).unwrap_err();
    drop(ostream);
    assert_eq!(
        native(stream.finish(Err(error.into()))),
        (6, format!("{post}{failed}"))
    );
    // A report that the file's buffer holds to `Close`: its flush fails, the run ends well
    // (NCBI 2.17.0, `blastp.f9.q1`: exit 0 after one 802.4).
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastp", Site::Report);
    let posts = posts_of(&stream);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    serial_object(&mut ostream, 3000, &[3000]).unwrap();
    drop(ostream);
    assert!(stream.finish(Ok(())).is_ok());
    assert_eq!(lock(&posts).len(), post.len());
    // On standard output the post's tie throws.
    let mut stream = ReportStream::new(
        Box::new(FailAt(0)),
        Some(NcbiBuffer::new(NcbiStdio::Cout {
            block: 4096,
            lines: false,
        })),
        true,
        "blastn",
        Site::Report,
    );
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    let error = serial_object(&mut ostream, 3000, &[3000]).unwrap_err();
    drop(ostream);
    assert_eq!(native(stream.finish(Err(error.into()))).0, -6);
}

#[test]
fn an_error_archive_into_the_failed_stream_is_a_write_fault() {
    use std::io::Write;
    // `PrintErrorArchive` after `CATCH_ALL` (probe `T/t6/survey` blastn.f11.q1 full2: 802.5
    // at byte 356, the error archive's length).
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Report);
    let posts = posts_of(&stream);
    let failed = stream.write_all(&[b'x'; 9000]).unwrap_err();
    assert_eq!(native(stream.finish(Err(failed.into()))).0, 6);
    stream.site(Site::Serial { whole: true }).unwrap();
    let error = stream.write_all(&[b'a'; 356]).unwrap_err();
    assert!(is_undecided_report_write_failure(&error));
    stream.write_fault(356).unwrap();
    stream.site(Site::Report).unwrap();
    assert_eq!(
        String::from_utf8(lock(&posts).clone()).unwrap(),
        "Error: (802.5) [blastn] CObjectOStream: error at byte 356: stack is empty: write \
             fault\n"
    );
}

#[test]
fn an_xml2_report_that_fails_is_caught_after_its_close() {
    // probe blastn.f16: 802.4 at `Close`, then the catch's post; the run goes on.
    let mut stream = ReportStream::file(Box::new(FailAt(100)), "blastn", Site::Report);
    let posts = posts_of(&stream);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    let report = [
        &b"<BlastOutput2>\n"[..],
        &b"  <x>12345</x>\n".repeat(700),
        b"</BlastOutput2>\n",
    ]
    .concat();
    let mut cuts = Cuts::default();
    cuts.push(0, Cut::Site(Site::Xml2Report { json: false }));
    cuts.push(report.len(), Cut::Flush);
    cuts.push(report.len(), Cut::Site(Site::Report));
    let error = cuts
        .cursor()
        .write(&mut ostream, &report, 0, report.len())
        .unwrap_err();
    assert!(is_caught_report_write_failure(&error));
    assert_eq!(
        String::from_utf8(lock(&posts).clone()).unwrap(),
        "Error: (802.4) [blastn] [Cannot close serializing output stream] Exception: \
             basic_ios::clear: iostream error\nError: [blastn] XML format failed\n"
    );
}

#[test]
fn a_held_xml2_report_caught_at_its_close_ends_the_run_at_what_follows_it() {
    // review D4 (outfmt 16 with `-subject`): the report that `PrintEpilog` writes is caught
    // (802.4, `XML format failed`), and the end of the document after it throws.
    let mut stream = ReportStream::file(Box::new(FailAt(100)), "blastn", Site::Report);
    let posts = posts_of(&stream);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    let report = [&b"  <x>12345</x>\n".repeat(700)[..], b"</BlastXML2>\n"].concat();
    let tail = report.len() - b"</BlastXML2>\n".len();
    let mut cuts = Cuts::default();
    cuts.push(0, Cut::Site(Site::Xml2Report { json: false }));
    cuts.push(tail, Cut::Flush);
    cuts.push(tail, Cut::Site(Site::Report));
    let error = cuts
        .cursor()
        .write(&mut ostream, &report, 0, report.len())
        .unwrap_err();
    assert!(is_report_write_failure(&error) && !is_caught_report_write_failure(&error));
    drop(ostream);
    assert_eq!(
        native(stream.finish(Err(error.into()))),
        (6, "BLAST failed to write output\n".to_string())
    );
    assert_eq!(
        String::from_utf8(lock(&posts).clone()).unwrap(),
        "Error: (802.4) [blastn] [Cannot close serializing output stream] Exception: \
             basic_ios::clear: iostream error\nError: [blastn] XML format failed\n"
    );
}

#[test]
fn a_tabular_report_flushed_at_each_query_fails_before_the_next_querys_post() {
    // review D2: `~CBlastTabularInfo` flushes after each query; on `-out /dev/full` the first
    // query's flush terminates, before the second query's warning.
    use std::io::Write;
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Tabular);
    let posts = posts_of(&stream);
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    ostream.write_all(b"w1\tsubject\t100.000\n").unwrap();
    ostream.ncbi_flush().unwrap();
    // The warning of the next query: NCBI's writes before it are made first, and fail.
    let error = lock(&stream.0)
        .post("Warning: [blastn] wN: ...\n")
        .unwrap_err();
    assert!(is_report_write_failure(&error));
    drop(ostream);
    assert!(lock(&posts).is_empty());
    assert_eq!(native(stream.finish(Err(error.into()))).0, -6);
}

/// A device that takes `limit` bytes and fails after them, keeping what it took.
struct Limited {
    took: std::sync::Arc<std::sync::Mutex<Vec<u8>>>,
    limit: usize,
}

impl std::io::Write for Limited {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let mut took = lock(&self.took);
        let room = self.limit - took.len();
        if room == 0 {
            return Err(std::io::Error::other("File too large"));
        }
        let n = buf.len().min(room);
        took.extend_from_slice(&buf[..n]);
        Ok(n)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// What a run of `steps` on a stream of `stream` (`standard_output`) over a device of `limit`
/// bytes shows: the bytes the device took, the posts, and the run's end; the formatter stops at
/// an error, and after a caught one goes on in the report's region. `immediate`: every write
/// of NCBI's made where NCBI makes it.
fn shown(
    steps: &[Step],
    stream: NcbiStdio,
    standard_output: bool,
    limit: usize,
    immediate: bool,
) -> (Vec<u8>, String, Result<(), (i32, String)>) {
    use std::io::Write;
    let took = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut report = ReportStream::new(
        Box::new(Limited {
            took: took.clone(),
            limit,
        }),
        Some(NcbiBuffer::new(stream)),
        standard_output,
        "blastn",
        Site::Report,
    );
    lock(&report.0).immediate = immediate;
    let posts = posts_of(&report);
    let mut ostream = NcbiOstream::report(stream_of(&mut report));
    let mut result = Ok(());
    let mut skipping = false;
    for step in steps {
        let done = match step {
            Step::Insert(_) if skipping => Ok(()),
            Step::Insert(bytes) => ostream.write_all(bytes),
            Step::Run(..) if skipping => Ok(()),
            Step::Run(bytes, ends) => {
                let mut cuts = Cuts::default();
                for &end in ends {
                    cuts.push(end, Cut::Insertion);
                }
                cuts.cursor().write(&mut ostream, bytes, 0, bytes.len())
            }
            Step::Flush if skipping => Ok(()),
            Step::Flush => ostream.ncbi_flush(),
            Step::Sync => lock(&report.0).sync(),
            Step::Site(site) => {
                skipping = false;
                match ostream.site(*site) {
                    // The region before ended caught; the new one has started.
                    Err(error) if is_caught_report_write_failure(&error) => Ok(()),
                    other => other,
                }
            }
            // A post in a region whose failure was caught belongs to what NCBI no longer
            // makes (the reads of the `Iteration`'s later `Hit`s).
            Step::Post(_) if skipping => Ok(()),
            Step::Post(text) => lock(&report.0).post(text),
        };
        match done {
            Err(error) if is_caught_report_write_failure(&error) => {
                skipping = true;
                ostream.site(Site::Report).unwrap();
            }
            Err(error) => {
                result = Err(error);
                break;
            }
            Ok(()) => {}
        }
    }
    drop(ostream);
    let end = report.finish(result.map_err(Into::into)).map_err(|error| {
        let native = error
            .downcast_ref::<crate::cli::NativeError>()
            .expect("NCBI's end");
        (
            native.exit,
            String::from_utf8(native.message.clone()).unwrap(),
        )
    });
    drop(report);
    let took = lock(&took).clone();
    let posts = String::from_utf8(lock(&posts).clone()).unwrap();
    (took, posts, end)
}

#[derive(Clone, Debug)]
enum Step {
    Insert(Vec<u8>),
    /// Insertions of the bytes up to each end (`Cuts`), in one call.
    Run(Vec<u8>, Vec<usize>),
    Flush,
    Sync,
    Site(Site),
    Post(String),
}

#[test]
fn ncbis_writes_made_later_together_end_a_run_as_where_ncbi_makes_them() {
    // A small deterministic generator (xorshift) of runs of insertions, flushes, regions and
    // posts; every limit of the device from 0 to past the report.
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = |n: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed % n
    };
    let sites = [
        Site::Report,
        Site::Alignments,
        Site::Tabular,
        Site::XmlIteration,
    ];
    for case in 0..120 {
        let mut steps = Vec::new();
        let mut total = 0usize;
        for _ in 0..(5 + next(40)) {
            let step = match next(12) {
                0..=4 => {
                    let large = next(4) == 0;
                    let n = 1 + next(if large { 9000 } else { 300 }) as usize;
                    let mut bytes = vec![b'x'; n];
                    if next(2) == 0 {
                        bytes[next(n as u64) as usize] = b'\n';
                    }
                    Step::Insert(bytes)
                }
                5 | 6 => {
                    let n = 1 + next(20_000) as usize;
                    let mut ends: Vec<usize> =
                        (0..next(30)).map(|_| next(n as u64) as usize).collect();
                    ends.sort_unstable();
                    Step::Run(vec![b'r'; n], ends)
                }
                7 => Step::Flush,
                8 => Step::Sync,
                9 | 10 => Step::Site(sites[next(4) as usize]),
                _ => Step::Post(format!("Warning: post {case}\n")),
            };
            if let Step::Insert(bytes) | Step::Run(bytes, _) = &step {
                total += bytes.len();
            }
            steps.push(step);
        }
        for (stream, standard_output) in [
            (NcbiStdio::File, false),
            (
                NcbiStdio::Cout {
                    block: 4096,
                    lines: false,
                },
                true,
            ),
        ] {
            let limits = (0..=total + 1).step_by(1 + total / 97);
            for limit in limits.chain([0, 1, 8190, 8191, 8192, 4096]) {
                let deferred = shown(&steps, stream, standard_output, limit, false);
                let reference = shown(&steps, stream, standard_output, limit, true);
                let brief = |shown: &(Vec<u8>, String, Result<(), (i32, String)>)| {
                    (shown.0.len(), shown.1.clone(), shown.2.clone())
                };
                let steps: Vec<String> = steps
                    .iter()
                    .map(|step| match step {
                        Step::Insert(bytes) => format!("I{}", bytes.len()),
                        Step::Run(bytes, ends) => format!("R{}{ends:?}", bytes.len()),
                        other => format!("{other:?}"),
                    })
                    .collect();
                assert!(
                    deferred == reference,
                    "case {case} {stream:?} limit {limit}: {:?} vs reference {:?}: {steps:?}",
                    brief(&deferred),
                    brief(&reference)
                );
            }
        }
    }
}

#[test]
fn a_post_while_standard_output_failed_terminates() {
    use std::io::Write;
    let mut stream = ReportStream::new(
        Box::new(FailAt(0)),
        Some(NcbiBuffer::new(NcbiStdio::Cout {
            block: 4096,
            lines: false,
        })),
        true,
        "blastn",
        Site::Report,
    );
    stream.write_all(b"held").unwrap();
    // `cerr`'s tie writes what `cout` holds and fails inside the post's destructor.
    assert!(lock(&stream.0).tie().is_err());
    assert_eq!(native(stream.finish(Ok(()))).0, -6);
}

#[test]
fn a_failed_standard_error_is_bad_and_the_run_goes_on_without_its_tie() {
    // review 2 D2: a post that standard error does not take makes `cerr` bad; the run goes on,
    // a later post writes nothing and makes none of standard output's writes before it.
    use std::io::Write;
    let took = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut stream = ReportStream::new(
        Box::new(Limited {
            took: took.clone(),
            limit: usize::MAX,
        }),
        Some(NcbiBuffer::new(NcbiStdio::Cout {
            block: 4096,
            lines: false,
        })),
        true,
        "blastn",
        Site::Report,
    );
    lock(&stream.0).diagnostics = Some((Box::new(FailAt(0)), Cerr::new()));
    stream.write_all(b"first").unwrap();
    // The first post's tie writes what standard output holds; standard error fails.
    lock(&stream.0).post("Warning: [blastn] w1\n").unwrap();
    assert_eq!(lock(&took).as_slice(), b"first");
    stream.write_all(b"second").unwrap();
    lock(&stream.0).post("Warning: [blastn] w2\n").unwrap();
    assert_eq!(lock(&took).as_slice(), b"first");
    assert!(stream.finish(Ok(())).is_ok());
    assert_eq!(lock(&took).as_slice(), b"firstsecond");
}

#[test]
fn the_request_point_of_a_refused_strategy_file_makes_ncbis_writes_first() {
    // merge 4's edge (fix 3): `-out /dev/full -export_search_strategy /dev/full` outfmt 0. NCBI's
    // prolog flush fails before its first request: `BLAST failed to write output`, exit 6.
    use std::io::Write;
    let mut stream = ReportStream::file(Box::new(FailAt(0)), "blastn", Site::Report);
    let point = RequestPoint(std::sync::Arc::downgrade(&stream.0));
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    ostream.write_all(b"BLASTN (LOSAT 0.2.0)\n\n").unwrap();
    ostream.ncbi_flush().unwrap();
    let failed = point.make_writes().unwrap_err();
    assert!(is_report_write_failure(&failed));
    drop(ostream);
    assert_eq!(
        native(stream.finish(Err(failed.into()))),
        (6, "BLAST failed to write output\n".to_string())
    );
    // A device that takes them: NCBI's writes so far are made, and the refusal follows them.
    let took = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut stream = ReportStream::file(
        Box::new(Limited {
            took: took.clone(),
            limit: usize::MAX,
        }),
        "blastn",
        Site::Report,
    );
    let point = RequestPoint(std::sync::Arc::downgrade(&stream.0));
    let mut ostream = NcbiOstream::report(stream_of(&mut stream));
    ostream.write_all(b"BLASTN (LOSAT 0.2.0)\n\n").unwrap();
    ostream.ncbi_flush().unwrap();
    ostream.write_all(b"Reference: ...").unwrap();
    assert!(lock(&took).is_empty());
    point.make_writes().unwrap();
    assert_eq!(lock(&took).as_slice(), b"BLASTN (LOSAT 0.2.0)\n\n");
    drop(ostream);
    let refusal = stream.finish(Err(anyhow::anyhow!("refused"))).unwrap_err();
    assert_eq!(refusal.to_string(), "refused");
    assert_eq!(
        lock(&took).as_slice(),
        b"BLASTN (LOSAT 0.2.0)\n\nReference: ..."
    );
}
