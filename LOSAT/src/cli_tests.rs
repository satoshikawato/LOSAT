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
        b"LOSAT database error: DB contains no accession info.\n"
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
    let text = "LOSAT query error: Could not construct seq-id from 'gi|abc'\n";
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
            let error = try_parse_from(argv).unwrap_err();
            let message = error.to_string();
            assert!(matches!(error, ArgvError::Refused(_)), "{message}");
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
// rejection after NCBI's reader, in the words of makeblastdb's other refusals (line LIB-B,
// probe 2026-11-21; line RG DECISIONS K4i-2b-2 (1)).
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
        let error = try_parse_from(argv).unwrap_err();
        let message = error.to_string();
        assert!(
            matches!(error, ArgvError::Refused(_)),
            "{option}: {message}"
        );
        assert!(
            message.contains(&format!("the value of {option} is not UTF-8"))
                && message.contains("not supported by LOSAT's makeblastdb"),
            "{message}"
        );
    }
}

// BLASTX reads both ranges (line SX, SXd1), excludes `-subject_loc` with `-db` as the other
// programs do, and rejects an option's value that is not UTF-8 with their text (DECISIONS 10
// (c) of line SX).
#[test]
fn blastx_reads_the_ranges_and_rejects_values_that_are_not_utf8() {
    let parsed = try_parse_from([
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
    let error = try_parse_from([
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
    // NCBI's reader: `-subject_loc` excludes `-db` (`Incompatible with argument`).
    assert!(
        matches!(&error, ArgvError::Report(text) if String::from_utf8_lossy(text).contains("Incompatible with argument")),
        "{error}"
    );
    for (option, value) in [
        ("-query_loc", vec![0xff]),
        ("-subject_loc", b"1-\xff".to_vec()),
        ("-outfmt", vec![0xfe]),
    ] {
        let mut argv = blastn_argv(vec![OsString::from(option), OsString::from_vec(value)]);
        argv[1] = OsString::from("blastx");
        let message = try_parse_from(argv).unwrap_err().to_string();
        assert!(
            message.contains(&format!("the value of {option} is not UTF-8"))
                && message.contains("not supported by LOSAT's BLASTX"),
            "{message}"
        );
    }
}

// The rejection comes after NCBI's reader, which reads the value's bytes: `-help` and the
// reader's errors come first.
#[test]
fn help_and_parser_errors_come_before_a_value_that_is_not_utf8() {
    let mut argv = vec![OsString::from("losat"), OsString::from("blastn")];
    argv.push(OsString::from("-help"));
    argv.push(OsString::from("-query_loc"));
    argv.push(OsString::from_vec(vec![0xff]));
    let error = try_parse_from(argv).unwrap_err();
    assert!(matches!(error, ArgvError::Printed(_)), "{error}");
    let argv = blastn_argv(vec![
        OsString::from("-query_loc"),
        OsString::from_vec(vec![0xff]),
        OsString::from("-no_such_option"),
    ]);
    let error = try_parse_from(argv).unwrap_err();
    assert!(matches!(error, ArgvError::Report(_)), "{error}");
    assert!(!error.to_string().contains("not UTF-8"), "{error}");
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
    let Commands::Blastn(args) = try_parse_from(argv).unwrap().command else {
        panic!("not blastn");
    };
    assert_eq!(args.query.into_os_string(), OsString::from_vec(name));
}

// No NCBI counterpart: a run builds the clap command of its program only (`ProgramParser`).
/// The programs of `Cli`, in the order of its subcommands.
fn program_names() -> Vec<String> {
    Cli::command()
        .get_subcommands()
        .map(|program| program.get_name().to_string())
        .collect()
}

/// What clap writes where it stops at `words` with `command` (a help, or an error), or the
/// program that it parsed.
fn stopped_text(command: &clap::Command, words: &[&str]) -> String {
    match command.clone().try_get_matches_from(words) {
        Err(error) => format!("{:?} {error}", error.kind()),
        Ok(matches) => format!("parsed {:?}", matches.subcommand_name()),
    }
}

// The command that the parser reads a program's command line with has that program's arguments
// and no other's: each option of another program costs a run nothing.
#[test]
fn a_command_line_builds_the_arguments_of_its_program_only() {
    let names = program_names();
    for name in &names {
        for words in [
            vec!["losat", name.as_str()],
            vec!["losat", "-help", name.as_str(), "-h"],
            vec!["losat", "--version", name.as_str()],
        ] {
            let words: Vec<OsString> = words.into_iter().map(OsString::from).collect();
            let command = Cli::program_command(&words);
            for program in command.get_subcommands() {
                let built = program.get_arguments().next().is_some();
                assert_eq!(built, program.get_name() == name, "{words:?}: {program:?}");
            }
        }
    }
    // A word that the reader rejects at the top level, before any program (`read_program`).
    for words in [
        &["losat"][..],
        &["losat", "-help"],
        &["losat", "nosuch", "blastp"],
        &["losat", "-x", "blastn", "blastp"],
        &["losat", "-help-full", "blastp"],
        &["losat", "--", "blastp"],
    ] {
        let words: Vec<OsString> = words.iter().map(OsString::from).collect();
        let program = Cli::program_command(&words)
            .get_subcommands()
            .find(|program| program.get_arguments().next().is_some())
            .map(|program| program.get_name().to_string());
        assert_eq!(program, None, "{words:?}");
    }
}

// The reader takes only its help and version flags before the program's name, so the program
// it reads is the first word after them (`read_program`): the top level has no option of its
// own, which another word there could be the value of.
#[test]
fn the_top_level_has_no_option_of_its_own() {
    assert!(Cli::command().get_arguments().next().is_none());
    for flag in TOP_LEVEL_FLAGS {
        let words = [
            OsString::from("losat"),
            OsString::from(flag),
            "blastn".into(),
        ];
        assert_eq!(
            read_program(&words, |word| word == "blastn"),
            Some("blastn")
        );
    }
}

// The program's command is the one that `Cli::command()` has, every byte (its arguments with the
// description and help template), and clap
// writes the same help and errors with it: the parser's texts do not depend on the other
// programs' arguments.
#[test]
fn the_command_of_one_program_is_its_command_among_every_program() {
    let every = Cli::command();
    for name in program_names() {
        let words = [OsString::from("losat"), OsString::from(&name)];
        let one = Cli::program_command(&words);
        assert_eq!(
            format!("{:?}", one.find_subcommand(&name)),
            format!("{:?}", every.find_subcommand(&name)),
            "{name}"
        );
        for words in [
            vec!["losat", name.as_str(), "-h"],
            vec!["losat", name.as_str(), "--help"],
            vec!["losat", "--help", name.as_str(), "-h"],
            vec!["losat", name.as_str(), "--no_such_option"],
            vec!["losat", name.as_str()],
        ] {
            assert_eq!(stopped_text(&one, &words), stopped_text(&every, &words));
        }
    }
    let top = Cli::program_command(&[OsString::from("losat")]);
    for words in [
        &["losat"][..],
        &["losat", "--help"],
        &["losat", "-h"],
        &["losat", "nosuch"],
    ] {
        assert_eq!(
            stopped_text(&top, words),
            stopped_text(&Cli::command(), words)
        );
    }
}

// `program_args` adds the arguments that clap's derive defers for each variant of `Commands`
// (`PROGRAM_ARGS` is in the variants' order), and keeps the variant's description, which clap's
// deferred `Args` replaces where a doc comment of an `Args` gives one.
#[test]
fn program_args_are_the_arguments_that_commands_defers() {
    let mut deferred = Commands::augment_subcommands(clap::Command::new("losat"));
    let descriptions: Vec<_> = deferred
        .get_subcommands()
        .map(|program| {
            (
                program.get_about().cloned(),
                program.get_long_about().cloned(),
            )
        })
        .collect();
    deferred.build();
    let mut ours = Cli::command();
    ours.build();
    assert_eq!(
        deferred.get_subcommands().count(),
        ours.get_subcommands().count()
    );
    for ((by_clap, by_us), (about, long_about)) in deferred
        .get_subcommands()
        .zip(ours.get_subcommands())
        .zip(descriptions)
    {
        assert_eq!(by_clap.get_name(), by_us.get_name());
        let arguments = |program: &clap::Command| -> Vec<String> {
            program
                .get_arguments()
                .map(|arg| format!("{arg:?}"))
                .collect()
        };
        assert_eq!(arguments(by_clap), arguments(by_us), "{}", by_us.get_name());
        assert_eq!(by_us.get_about(), about.as_ref(), "{}", by_us.get_name());
        assert_eq!(
            by_us.get_long_about(),
            long_about.as_ref(),
            "{}",
            by_us.get_name()
        );
    }
}

// Every variant of `Commands` has a description, which `program_args` keeps as the derive's
// final attributes do, and no program's `Args` sets the help template: the variant's command is
// the derive's where clap defers the arguments too.
#[test]
fn every_program_keeps_its_description_and_help_template() {
    // One empty command with the help template "T": the parent of the variants here, and the
    // command that each program's `Args` alone is added to.
    let blank = clap::Command::new("alone").help_template("T");
    let programs = Commands::augment_subcommands(blank.clone());
    assert_eq!(programs.get_subcommands().count(), PROGRAM_ARGS.len());
    for (program, args) in programs.get_subcommands().zip(PROGRAM_ARGS) {
        assert!(program.get_about().is_some(), "{}", program.get_name());
        let mut alone = args(blank.clone());
        assert_eq!(
            alone.render_help().to_string().trim_end(),
            "T",
            "{}",
            program.get_name()
        );
    }
}

// A program that is not selected is deferred (`cli_command`); where clap builds it, it has every
// argument and help text of the program's command among every program.
#[test]
fn a_deferred_program_is_built_with_its_arguments() {
    let mut deferred = cli_command::<false>(|_| false);
    deferred.build();
    let mut every = Cli::command();
    every.build();
    assert_eq!(
        deferred.get_subcommands().count(),
        every.get_subcommands().count()
    );
    for (by_clap, ours) in deferred
        .get_subcommands_mut()
        .zip(every.get_subcommands_mut())
    {
        let name = ours.get_name().to_string();
        assert_eq!(by_clap.get_name(), name);
        let arguments = |program: &clap::Command| -> Vec<String> {
            program
                .get_arguments()
                .map(|arg| format!("{arg:?}"))
                .collect()
        };
        assert_eq!(arguments(by_clap), arguments(ours), "{name}");
        assert_eq!(
            by_clap.render_help().to_string(),
            ours.render_help().to_string(),
            "{name}"
        );
        assert_eq!(
            by_clap.render_long_help().to_string(),
            ours.render_long_help().to_string(),
            "{name}"
        );
    }
}
