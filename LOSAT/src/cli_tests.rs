use super::*;
use std::os::unix::ffi::OsStringExt;

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

// Line DBID pfoff-1 (DECISIONS 62): `CATCH_ALL` reports a `CSeqDBException` (blastp 2.17.0,
// outfmt 13 of `pP pU` with the prefetch off: `BLAST Database error: DB contains no accession
// info.`, exit 2) before any other `CException` (`Error: ` and `what()`, exit 255).
#[test]
fn catch_all_reports_a_seqdb_exception_as_a_database_error() {
    use crate::ncbi_exception::{NcbiException, ThrowSite};
    let site = crate::seqdb::lookup_throw_site(b"DB contains no accession info.").unwrap();
    let error = NativeError::exception(&NcbiException::new(site, "DB contains no accession info."));
    assert_eq!(error.exit, 2);
    assert_eq!(
        error.message,
        b"BLAST Database error: DB contains no accession info.\n"
    );
    let other = NcbiException::new(
        ThrowSite {
            file: "src/corelib/ncbistr.cpp",
            line: 852,
            kind: "CStringException::eConvert",
            function: Some("ncbi::NStr::StringToInt8"),
        },
        "x",
    );
    let error = NativeError::exception(&other);
    assert_eq!(error.exit, 255);
    assert_eq!(
        error.message,
        [b"Error: ".as_slice(), &other.what(), b"\n"].concat()
    );
}

#[test]
fn the_organism_reports_file_buffer_writes_at_its_8191st_byte() {
    // probe `T/t2/fix1` `f1.d5.*` and `stdout_probe.sh`: NCBI exits 0 for a report that its
    // buffer holds to the end, 6 where the buffer writes.
    use std::io::Write;
    let offered = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let mut buffer = NcbiReportBuffer::file(Box::new(Full(offered.clone())));
    buffer.write_all(&[b'x'; 8190]).unwrap();
    buffer.flush().unwrap();
    assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 0);
    assert!(buffer.write_all(b"x").is_err());
    assert_eq!(offered.load(std::sync::atomic::Ordering::Relaxed), 8191);
    // What a buffer holds at the end is written, a failure unreported.
    let quiet = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let mut buffer = NcbiReportBuffer::file(Box::new(Full(quiet.clone())));
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
fn the_organism_reports_standard_output_writes_where_glibc_writes() {
    // `std::cout` on a file of `st_blksize` 4096 (observed, glibc 2.35 and 2.39,
    // `tests/fixtures/stdio_buffer/cout_file.tsv`): the writes end at multiples of the block
    // and the rest stays held; the first insertion has no buffer yet; a flush empties the
    // buffer.
    use std::io::Write;
    let run = |pieces: &[usize]| {
        let writes = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut after = Vec::new();
        let mut buffer = NcbiReportBuffer::new(
            Box::new(Writes(writes.clone())),
            NcbiBuffer::new(NcbiStdio::Cout {
                block: 4096,
                lines: false,
            }),
        );
        for &piece in pieces {
            match piece {
                0 => lock(&buffer.0).write_held().unwrap(),
                n => buffer.write_all(&vec![b'x'; n]).unwrap(),
            }
            after.push(lock(&writes).iter().sum::<usize>());
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
    let mut buffer =
        NcbiReportBuffer::new(Box::new(Writes(writes.clone())), NcbiBuffer::new(stream));
    for step in steps {
        match step {
            Ok(bytes) => buffer.write_all(bytes).unwrap(),
            Err(()) => lock(&buffer.0).write_held().unwrap(),
        }
        after.push(lock(&writes).iter().sum::<usize>());
    }
    drop(buffer);
    let writes = lock(&writes);
    (after, writes.iter().sum(), writes.len())
}

#[test]
fn a_terminals_standard_output_writes_at_its_newlines_as_glibc_does() {
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
                let mut buffer = NcbiReportBuffer::new(
                    Box::new(Writes(writes.clone())),
                    NcbiBuffer::new(stream),
                );
                let mut pieces = Vec::new();
                std::io::Write::write_all(&mut buffer, &lead).unwrap();
                pieces.push(lock(&writes).iter().sum::<usize>());
                streamed(&mut buffer);
                pieces.push(lock(&writes).iter().sum::<usize>());
                std::io::Write::write_all(&mut buffer, after).unwrap();
                pieces.push(lock(&writes).iter().sum::<usize>());
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

fn blastn_argv(words: Vec<OsString>) -> Vec<OsString> {
    let mut argv = vec![
        OsString::from("losat"),
        OsString::from("blastn"),
        OsString::from("-query"),
        OsString::from("q.fa"),
        OsString::from("-subject"),
        OsString::from("s.fa"),
    ];
    argv.extend(words);
    argv
}

// Without `main`'s `post_errors_in_run` (a library caller of the command line entries such as
// `api::run_blastn`), the run writes no error text of its own: the text stays in the error.
// No unit test calls `post_errors_in_run`.
#[test]
fn a_run_outside_main_posts_no_error_text() {
    let text = "BLAST query error: Could not construct seq-id from 'gi|abc'\n";
    let error: anyhow::Error = NativeError {
        exit: 1,
        message: text.to_string().into(),
    }
    .into();
    post_native_error(&error);
    assert!(!posted(text.as_bytes()));
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:266-270 keeps the bytes of a string value; LOSAT
// rejects a value that is not UTF-8 explicitly, in both spellings of an option.
#[test]
fn option_values_that_are_not_utf8_are_rejected_explicitly() {
    for (option, value) in [
        ("-query_loc", vec![0xff]),
        ("-subject_loc", b"1-\xff".to_vec()),
        ("-dust", vec![0xfe]),
    ] {
        let separate = blastn_argv(vec![
            OsString::from(option),
            OsString::from_vec(value.clone()),
        ]);
        let mut inline = format!("{option}=").into_bytes();
        inline.extend(&value);
        let joined = blastn_argv(vec![OsString::from_vec(inline)]);
        for argv in [separate, joined] {
            let error = try_parse_from::<Cli, _, _>(argv).unwrap_err();
            let message = render_message(&error);
            assert!(
                message.contains(&format!("the value of {option} is not UTF-8"))
                    && message.contains("not supported by LOSAT's BLASTN"),
                "{message}"
            );
        }
    }
}

// makeblastdb keeps its values as UTF-8 strings: a value that is not UTF-8 (a file name of
// `-in`, `-out`, `-mask_data`, `-taxid_map`, `-logfile`, a `-title`) is LOSAT's explicit
// rejection, where the parser's own error named no option (line LIB-B, probe 2026-11-21).
#[test]
fn makeblastdb_rejects_values_that_are_not_utf8_explicitly() {
    for option in [
        "-in",
        "-out",
        "-mask_data",
        "-taxid_map",
        "-logfile",
        "-title",
    ] {
        let argv = vec![
            OsString::from("losat"),
            OsString::from("makeblastdb"),
            OsString::from("-dbtype"),
            OsString::from("nucl"),
            // `-taxid_map` requires it: the parser's errors come before the rejection.
            OsString::from("-parse_seqids"),
            OsString::from(option),
            OsString::from_vec(b"n\xe9.x".to_vec()),
        ];
        let error = try_parse_from::<Cli, _, _>(argv).unwrap_err();
        let message = render_message(&error);
        assert_eq!(error.kind(), ErrorKind::InvalidUtf8, "{option}: {message}");
        assert!(
            message.contains(&format!("the value of {option} is not UTF-8"))
                && message.contains("not supported by LOSAT's MAKEBLASTDB"),
            "{message}"
        );
    }
}

// BLASTX reads both ranges (line SX, SXd1), excludes `-subject_loc` with `-db` as the other
// programs do, and rejects an option's value that is not UTF-8 with their text (DECISIONS 10
// (c) of line SX).
#[test]
fn blastx_reads_the_ranges_and_rejects_values_that_are_not_utf8() {
    let parsed = try_parse_from::<Cli, _, _>([
        "losat",
        "blastx",
        "-query",
        "q.fa",
        "-subject",
        "s.faa",
        "-query_loc",
        "1-10",
        "-subject_loc",
        "2-5",
    ])
    .unwrap();
    let Commands::Blastx(args) = parsed.command else {
        panic!("a BLASTX command line");
    };
    assert_eq!(args.query_loc.as_deref(), Some("1-10"));
    assert_eq!(args.subject_loc.as_deref(), Some("2-5"));
    let error = try_parse_from::<Cli, _, _>([
        "losat",
        "blastx",
        "-query",
        "q.fa",
        "-db",
        "d",
        "-subject_loc",
        "2-5",
    ])
    .unwrap_err();
    assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
    for (option, value) in [
        ("-query_loc", vec![0xff]),
        ("-subject_loc", b"1-\xff".to_vec()),
        ("-outfmt", vec![0xfe]),
    ] {
        let mut argv = blastn_argv(vec![OsString::from(option), OsString::from_vec(value)]);
        argv[1] = OsString::from("blastx");
        let message = render_message(&try_parse_from::<Cli, _, _>(argv).unwrap_err());
        assert!(
            message.contains(&format!("the value of {option} is not UTF-8"))
                && message.contains("not supported by LOSAT's BLASTX"),
            "{message}"
        );
    }
}

// The rejection comes after the parser, as NCBI reads the value only after its
// arguments parse: `-help` and the parser's errors come first.
#[test]
fn help_and_parser_errors_come_before_a_value_that_is_not_utf8() {
    let mut argv = vec![OsString::from("losat"), OsString::from("blastn")];
    argv.push(OsString::from("-help"));
    argv.push(OsString::from("-query_loc"));
    argv.push(OsString::from_vec(vec![0xff]));
    let error = try_parse_from::<Cli, _, _>(argv).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::DisplayHelp);
    let argv = blastn_argv(vec![
        OsString::from("-query_loc"),
        OsString::from_vec(vec![0xff]),
        OsString::from("-no_such_option"),
    ]);
    let message = render_message(&try_parse_from::<Cli, _, _>(argv).unwrap_err());
    assert!(!message.contains("not UTF-8"), "{message}");
}

// `try_parse_from` is generic: a caller's own type can name an option with letters that are
// not ASCII, here U+FFFD (3 bytes), which a byte that is not UTF-8 also reads as in the
// word's lossy text. The value of `-name=value` starts after the first `=` of the bytes,
// not `name.len()` bytes on (which skipped `\xfe\xfd` and read `abc`).
#[test]
fn inline_value_starts_after_the_first_equals_sign() {
    #[derive(clap::Parser)]
    enum Own {
        Blastn(OwnArgs),
    }
    #[derive(clap::Args)]
    struct OwnArgs {
        #[arg(long = "q\u{fffd}")]
        q: Option<String>,
    }
    let argv = vec![
        OsString::from("losat"),
        OsString::from("blastn"),
        OsString::from_vec(b"-q\xff=\xfe\xfdabc".to_vec()),
    ];
    let error = try_parse_from::<Own, _, _>(argv)
        .err()
        .expect("a value that is not UTF-8");
    assert_eq!(error.kind(), ErrorKind::InvalidUtf8);
}

// A file name keeps its bytes in both spellings (`check_utf8_file_name` reads it).
#[test]
fn file_names_keep_their_bytes_in_both_spellings() {
    let name = b"q\xff.fa".to_vec();
    let mut inline = b"-query=".to_vec();
    inline.extend(&name);
    let argv = vec![
        OsString::from("losat"),
        OsString::from("blastn"),
        OsString::from_vec(inline),
        OsString::from("-subject"),
        OsString::from("s.fa"),
    ];
    let Commands::Blastn(args) = try_parse_from::<Cli, _, _>(argv).unwrap().command else {
        panic!("not blastn");
    };
    assert_eq!(args.query.into_os_string(), OsString::from_vec(name));
}
