use super::*;

/// The emitted calls, in order, taken.
fn emitted() -> Vec<(u32, Vec<u8>)> {
    crate::emit::EMITTED.with(|emitted| std::mem::take(&mut *emitted.borrow_mut()))
}

/// The bytes of `stream` in `calls`, joined.
fn joined(calls: &[(u32, Vec<u8>)], stream: u32) -> Vec<u8> {
    calls
        .iter()
        .filter(|(known, _)| *known == stream)
        .flat_map(|(_, bytes)| bytes.iter().copied())
        .collect()
}

/// A frame of stream 12: the file's number, the kind, the payload.
fn frames(calls: &[(u32, Vec<u8>)]) -> Vec<(u32, u32, &[u8])> {
    calls
        .iter()
        .filter(|(stream, _)| *stream == STREAM_TOOL_FILES)
        .map(|(_, frame)| {
            assert!(frame.len() >= 8 && frame.len() <= CHUNK, "one frame a call");
            let number = u32::from_le_bytes(frame[..4].try_into().unwrap());
            let kind = u32::from_le_bytes(frame[4..8].try_into().unwrap());
            (number, kind, &frame[8..])
        })
        .collect()
}

/// The files of the frames, in the order of their beginnings: the name and the bytes, each frame
/// of bytes written at its offset and each truncation applied, as a host applies them; every file
/// begins once, then has bytes and truncations, then ends once, the ends last, in order.
fn files(frames: &[(u32, u32, &[u8])]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut files: Vec<(Vec<u8>, Vec<u8>, bool)> = Vec::new();
    let number64 = |bytes: &[u8]| u64::from_le_bytes(bytes[..8].try_into().unwrap()) as usize;
    for &(number, kind, payload) in frames {
        let number = number as usize;
        assert!(
            kind == FRAME_BEGIN || !files[number].2,
            "a frame after the end"
        );
        match kind {
            FRAME_BEGIN => {
                assert_eq!(number, files.len(), "numbered in order");
                assert!(files.iter().all(|file| !file.2), "a file after an end");
                files.push((payload.to_vec(), Vec::new(), false));
            }
            FRAME_BYTES => {
                let (offset, bytes) = (number64(payload), &payload[8..]);
                assert!(!bytes.is_empty(), "a frame of no bytes");
                let file = &mut files[number].1;
                if file.len() < offset + bytes.len() {
                    file.resize(offset + bytes.len(), 0);
                }
                file[offset..offset + bytes.len()].copy_from_slice(bytes);
            }
            FRAME_TRUNCATE => {
                assert_eq!(payload.len(), 8);
                files[number].1.resize(number64(payload), 0);
            }
            FRAME_END => {
                assert!(payload.is_empty(), "an end has no payload");
                assert!(files[..number].iter().all(|file| file.2), "ends in order");
                files[number].2 = true;
            }
            other => panic!("kind {other}"),
        }
    }
    assert!(files.iter().all(|(_, _, ended)| *ended), "every file ends");
    files
        .into_iter()
        .map(|(name, bytes, _)| (name, bytes))
        .collect()
}

/// A fixture of LOSAT's windowmasker tests, by its absolute path (the adapter's tests run in
/// web/adapter).
fn fixture(path: &str) -> String {
    format!("{}/../../LOSAT/tests/{path}", env!("CARGO_MANIFEST_DIR"))
}

/// A folder of the test's own, made empty.
fn folder(name: &str) -> std::path::PathBuf {
    let folder = std::env::temp_dir().join(format!("ww3a-adapter-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&folder);
    std::fs::create_dir_all(&folder).unwrap();
    folder
}

/// The library's run of `words` (after the name) on the disk, whose `-out` and `-logfile` are the
/// files of their names in `folder`: the standard output and error, and the bytes of each file of
/// `names` in `folder` after the run.
fn in_folder(
    words: &[&str],
    folder: &std::path::Path,
    names: &[&str],
) -> (Vec<u8>, Vec<u8>, Vec<Vec<u8>>) {
    let paths: Vec<String> = words
        .iter()
        .map(|word| folder.join(word).to_str().unwrap().to_string())
        .collect();
    let mut replaced = Vec::new();
    for (at, &word) in words.iter().enumerate() {
        let named = at > 0 && matches!(words[at - 1], "-out" | "-logfile");
        replaced.push(if named { paths[at].as_str() } else { word });
    }
    let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
    let _ = losat::api::run_windowmasker(&replaced, &mut io::empty(), &mut stdout, &mut stderr);
    let files = names
        .iter()
        .map(|name| std::fs::read(folder.join(name)).unwrap())
        .collect();
    (stdout, stderr, files)
}

/// The library's run of `words` (after the name) on the disk, whose `-out` and `-logfile` are the
/// files `out` and `log` of `folder`: the standard output and error, and the two files.
fn on_disk(words: &[&str], folder: &std::path::Path) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let out = folder.join("out");
    let log = folder.join("log");
    let (out_name, log_name) = (out.to_str().unwrap(), log.to_str().unwrap());
    let mut replaced = Vec::new();
    let mut previous = "";
    for &word in words {
        replaced.push(match previous {
            "-out" => out_name,
            "-logfile" => log_name,
            _ => word,
        });
        previous = word;
    }
    let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
    let _ = losat::api::run_windowmasker(&replaced, &mut io::empty(), &mut stdout, &mut stderr);
    (
        stdout,
        stderr,
        std::fs::read(&out).unwrap(),
        std::fs::read(&log).unwrap(),
    )
}

/// A tool run of `-mk_counts` with `-out FILE` and `-logfile FILE` gives the files on stream 12 as
/// the library writes them on the disk, under the names of the argv (the file of `-logfile`, which
/// the run opens first, is file 0), no standard output, and the exit status 0 on stream 2 last;
/// nothing on the disk.
#[test]
fn a_tool_run_gives_the_command_line_s_files() {
    emitted();
    let fasta = fixture("winmask/fasta/dm2L_3k.fa");
    let words = [
        "windowmasker",
        "-mk_counts",
        "-in",
        &fasta,
        "-unit",
        "5",
        "-sformat",
        "obinary",
        "-out",
        "ww3a-tool.counts",
        "-logfile",
        "ww3a-tool.log",
    ];
    run(&words).unwrap();
    let calls = emitted();
    assert_eq!(
        calls.last().unwrap(),
        &(STREAM_RESPONSE, b"{\"status\":0}".to_vec())
    );
    assert!(joined(&calls, STREAM_PRINTED).is_empty());
    assert!(calls
        .iter()
        .all(|(stream, _)| *stream != crate::emit::STREAM_RUN_FILES));
    let files = files(&frames(&calls));
    let names: Vec<&[u8]> = files.iter().map(|(name, _)| name.as_slice()).collect();
    assert_eq!(names, [b"ww3a-tool.log".as_slice(), b"ww3a-tool.counts"]);
    assert!(
        !std::path::Path::new("ww3a-tool.counts").exists(),
        "no file on the disk"
    );
    let folder = folder("counts");
    let (stdout, stderr, out, log) = on_disk(&words[1..], &folder);
    assert!(stdout.is_empty());
    assert_eq!(joined(&calls, STREAM_DIAGNOSTICS), stderr);
    assert_eq!(files[1].1, out, "the counts");
    assert_eq!(files[0].1, log, "the file of -logfile");
    std::fs::remove_dir_all(&folder).unwrap();
}

/// A file of a tool run goes out as it is written, in frames of at most a chunk (D-WW3-6): the
/// masked FASTA of a sequence of 1.2 Mbp comes in two frames of bytes or more, and the files end
/// when the run ends, in the order of their numbers (the file of `-logfile` first); its bytes are
/// the library's on the disk.
#[test]
fn a_tool_run_streams_its_files_while_it_runs() {
    emitted();
    let folder = folder("stream");
    // A sequence of 1.2 Mbp of the four bases (a linear congruential generator).
    let mut state: u32 = 1;
    let mut fasta = b">big\n".to_vec();
    for line in 0..20_000 {
        for _ in 0..60 {
            state = state
                .wrapping_mul(1_103_515_245)
                .wrapping_add(12_345 + line % 7);
            fasta.push(b"ACGT"[(state >> 16) as usize % 4]);
        }
        fasta.push(b'\n');
    }
    let input = folder.join("big.fa");
    std::fs::write(&input, &fasta).unwrap();
    let counts = fixture("winmask/counts/mg1655_u9.obinary");
    let input = input.to_str().unwrap().to_string();
    let words = [
        "windowmasker",
        "-ustat",
        &counts,
        "-in",
        &input,
        "-outfmt",
        "fasta",
        "-out",
        "ww3a-big.masked",
        "-logfile",
        "ww3a-big.log",
    ];
    run(&words).unwrap();
    let calls = emitted();
    assert_eq!(
        calls.last().unwrap(),
        &(STREAM_RESPONSE, b"{\"status\":0}".to_vec())
    );
    let frames = frames(&calls);
    // The files end when the run ends, in the order of their numbers (the log's file first); the
    // output's last bytes, which it writes as it closes, may come after the log's end.
    let ends: Vec<u32> = frames
        .iter()
        .filter(|frame| frame.1 == FRAME_END)
        .map(|frame| frame.0)
        .collect();
    assert_eq!(ends, [0, 1], "the files end with the run");
    assert_eq!(
        frames.last().map(|frame| (frame.0, frame.1)),
        Some((1, FRAME_END))
    );
    let bytes = frames
        .iter()
        .filter(|&&(number, kind, _)| number == 1 && kind == FRAME_BYTES)
        .count();
    assert!(bytes >= 2, "{bytes} frames of bytes");
    let files = files(&frames);
    assert_eq!(files[1].0, b"ww3a-big.masked");
    let (_, _, out, log) = on_disk(&words[1..], &folder);
    assert!(out.len() > CHUNK);
    assert_eq!(files[1].1, out, "the masked FASTA");
    assert_eq!(files[0].1, log);
    std::fs::remove_dir_all(&folder).unwrap();
}

/// A file's frames go out as the run writes it (review of WW3a, L9): a frame of bytes as soon as
/// one fills, before the file ends, each with its offset; writes that follow each other are
/// joined, a write elsewhere or a truncation sends what waits first, and the end sends the rest.
#[test]
fn a_file_s_frames_go_out_before_it_ends() {
    emitted();
    let files = ToolFiles::default();
    let mut file = files.open(b"x").unwrap();
    let data = vec![7u8; FRAME_PAYLOAD + 10];
    file.write_at(0, &data[..FRAME_PAYLOAD - 5]).unwrap();
    assert_eq!(frames(&emitted()).len(), 1, "the beginning alone");
    file.write_at((FRAME_PAYLOAD - 5) as u64, &data[FRAME_PAYLOAD - 5..])
        .unwrap();
    let calls = emitted();
    let sent = frames(&calls);
    assert_eq!(sent.len(), 1, "one full frame, before the end");
    assert_eq!((sent[0].0, sent[0].1), (0, FRAME_BYTES));
    assert_eq!(&sent[0].2[..8], &0u64.to_le_bytes());
    assert_eq!(sent[0].2.len(), 8 + FRAME_PAYLOAD);
    file.write_at(3, b"abc").unwrap();
    file.set_len(2).unwrap();
    drop(file);
    let calls = emitted();
    let sent = frames(&calls);
    let kinds: Vec<u32> = sent.iter().map(|frame| frame.1).collect();
    assert_eq!(
        kinds,
        [FRAME_BYTES, FRAME_BYTES, FRAME_TRUNCATE, FRAME_END],
        "the rest at its offset, then 3 at 3, the truncation, the end"
    );
    assert_eq!(
        sent[0].2,
        [&(FRAME_PAYLOAD as u64).to_le_bytes()[..], &[7u8; 10]].concat()
    );
    assert_eq!(sent[1].2, [&3u64.to_le_bytes()[..], b"abc"].concat());
    assert_eq!(sent[2].2, 2u64.to_le_bytes());
}

/// A run that opens one file twice (`-out x -logfile x`, `-out ./x -logfile x`) gives one file with
/// the bytes of the command line's, on stream 12 (review of WW3a, M2, D-WW3-9): the log's messages
/// go to the end of the file, the output's opening empties it and writes from its start; the
/// frames come as the run writes (`-convert` logs a message before it opens its output, which
/// truncates it; `-mk_counts` opens its output first), and the output's bytes are written over a
/// message that the log added at the end before them.
#[test]
fn a_file_opened_twice_is_one_file() {
    emitted();
    let counts = fixture("windowmasker/counts/u4.obinary");
    let fasta = fixture("winmask/fasta/dm2L_3k.fa");
    for (case, words, first) in [
        (
            "convert",
            vec![
                "-convert", "-in", &counts, "-sformat", "oascii", "-out", "x", "-logfile", "x",
            ],
            [FRAME_BEGIN, FRAME_BYTES, FRAME_TRUNCATE],
        ),
        (
            "mk_counts",
            vec![
                "-mk_counts",
                "-in",
                &fasta,
                "-unit",
                "5",
                "-sformat",
                "ascii",
                "-out",
                "./x",
                "-logfile",
                "x",
            ],
            [FRAME_BEGIN, FRAME_TRUNCATE, FRAME_BYTES],
        ),
    ] {
        let mut argv = vec!["windowmasker"];
        argv.extend(&words);
        run(&argv).unwrap();
        let calls = emitted();
        let frames = frames(&calls);
        let kinds: Vec<u32> = frames.iter().take(3).map(|frame| frame.1).collect();
        assert_eq!(kinds, first, "{case}");
        // A frame of bytes before the file's end: the output's write over the log's message.
        let mut length = 0;
        let mut before_end = false;
        for &(_, kind, payload) in &frames {
            let at = || u64::from_le_bytes(payload[..8].try_into().unwrap());
            if kind == FRAME_BYTES {
                before_end |= at() < length;
                length = length.max(at() + (payload.len() - 8) as u64);
            } else if kind == FRAME_TRUNCATE {
                length = at();
            }
        }
        assert!(before_end, "{case}: a write before the end");
        let files = files(&frames);
        let folder = folder(case);
        let (stdout, stderr, on_disk) = in_folder(&words, &folder, &["x"]);
        assert_eq!(files.len(), 1, "{case}: one file");
        assert_eq!(files[0].0, b"x");
        assert_eq!(files[0].1, on_disk[0], "{case}: the command line's bytes");
        assert_eq!(joined(&calls, STREAM_PRINTED), stdout);
        assert_eq!(joined(&calls, STREAM_DIAGNOSTICS), stderr);
        std::fs::remove_dir_all(&folder).unwrap();
    }
}

/// A tool run's standard error goes out as it is written (review of WW3a, L6): a message at a
/// time, not when the run ends.
#[test]
fn a_tool_run_s_standard_error_goes_out_as_it_is_written() {
    emitted();
    let fasta = fixture("winmask/fasta/dm2L_3k.fa");
    run(&["windowmasker", "-mk_counts", "-in", &fasta, "-unit", "5"]).unwrap();
    let calls = emitted();
    let messages: Vec<&[u8]> = calls
        .iter()
        .filter(|(stream, _)| *stream == STREAM_DIAGNOSTICS)
        .map(|(_, bytes)| bytes.as_slice())
        .collect();
    assert_eq!(messages, [b"pass 1\n".as_slice(), b"pass 2\n"]);
}

/// A tool run's standard output is stream 11 (`-out -`, the help, NCBI's XML of `-xmlhelp`, the
/// version), its exit status that of the command line, also for a run that fails (stream 3 then
/// has the error), and a folder that the run's directory does not have opens no counts file.
#[test]
fn a_tool_run_prints_and_fails_as_the_command_line() {
    emitted();
    let fasta = fixture("winmask/fasta/dm2L_3k.fa");
    run(&["windowmasker", "-mk_counts", "-in", &fasta, "-unit", "5"]).unwrap();
    let calls = emitted();
    let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
    losat::api::run_windowmasker(
        &["-mk_counts", "-in", &fasta, "-unit", "5"],
        &mut io::empty(),
        &mut stdout,
        &mut stderr,
    )
    .unwrap();
    assert_eq!(joined(&calls, STREAM_PRINTED), stdout);
    assert_eq!(joined(&calls, STREAM_DIAGNOSTICS), stderr);
    assert!(frames(&calls).is_empty());
    run(&[
        "windowmasker",
        "-mk_counts",
        "-in",
        &fasta,
        "-out",
        "ww3a-no-such-folder/c.counts",
    ])
    .unwrap();
    let calls = emitted();
    assert!(frames(&calls).is_empty(), "no file");
    run(&["windowmasker", "-in", &fasta]).unwrap();
    let calls = emitted();
    assert_eq!(
        calls.last().unwrap(),
        &(STREAM_RESPONSE, b"{\"status\":3}".to_vec())
    );
    assert!(String::from_utf8_lossy(&joined(&calls, STREAM_DIAGNOSTICS))
        .contains("one of '-mk_counts', '-convert' or '-ustat <stat_file>' must be specified"));
    run(&["windowmasker", "-help"]).unwrap();
    assert!(joined(&emitted(), STREAM_PRINTED).starts_with(b"USAGE\n  windowmasker"));
    run(&["windowmasker", "-xmlhelp"]).unwrap();
    let calls = emitted();
    assert!(joined(&calls, STREAM_PRINTED).starts_with(b"<?xml"));
    assert_eq!(
        calls.last().unwrap(),
        &(STREAM_RESPONSE, b"{\"status\":0}".to_vec())
    );
    run(&["windowmasker", "-version"]).unwrap();
    let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
    losat::api::run_windowmasker(&["-version"], &mut io::empty(), &mut stdout, &mut stderr)
        .unwrap();
    assert!(!stdout.is_empty());
    assert_eq!(joined(&emitted(), STREAM_PRINTED), stdout);
}

/// The validation of a tool's argv (D-WW3-7, D-WW3-10): the command line's checks before the tool
/// reads its input or opens a file pass a good argv (also one of a missing input, which the run
/// reads, and one whose error comes after a file: the masking's `-ids` with `-exclude_ids`, after
/// `-in`; `-infmt seqids` after the id list), and refuse the others with the command line's error:
/// its standard error, without the last line end; nothing emitted.
#[test]
fn a_tool_s_argv_is_checked_as_the_command_line_checks_it() {
    emitted();
    assert_eq!(
        validate(&["windowmasker", "-mk_counts", "-in", "missing.fa"]),
        Ok(())
    );
    for passed in [
        &[
            "-ustat",
            "c",
            "-in",
            "missing.fa",
            "-ids",
            "a",
            "-exclude_ids",
            "b",
        ][..],
        &["-mk_counts", "-infmt", "seqids", "-ids", "missing.txt"],
        &["-convert", "-in", "missing.counts"],
    ] {
        let words = [&["windowmasker"][..], passed].concat();
        assert_eq!(validate(&words), Ok(()), "{passed:?}");
    }
    for refused in [
        &["-mk_counts", "-ids", "a", "-exclude_ids", "b"][..],
        &["-mk_counts", "-infmt", "seqids"],
        &["-convert", "-in", "-"],
        &["-convert", "-out", "-"],
    ] {
        let words = [&["windowmasker"][..], refused].concat();
        let error = validate(&words).unwrap_err();
        let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
        let failed =
            losat::api::run_windowmasker(refused, &mut io::empty(), &mut stdout, &mut stderr);
        assert!(failed.is_err(), "{refused:?}");
        assert_eq!(format!("{error}\n").as_bytes(), stderr, "{refused:?}");
    }
    assert_eq!(validate(&["windowmasker", "-help"]), Ok(()));
    let error = validate(&["windowmasker", "-in", "x.fa"]).unwrap_err();
    assert!(error.contains("must be specified"), "{error}");
    let error = validate(&["windowmasker", "-mk_counts", "-nosuchoption"]).unwrap_err();
    assert!(error.starts_with("USAGE\n"), "{error}");
    assert!(!error.ends_with('\n'));
    let words = ["-mk_counts", "-nosuchoption"];
    let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
    let _ = losat::api::run_windowmasker(&words, &mut io::empty(), &mut stdout, &mut stderr);
    assert_eq!(format!("{error}\n").as_bytes(), stderr);
    assert!(emitted().is_empty());
}

/// The names that are not a tool of LOSAT Web: the searches and the tools that it does not run yet
/// (`not available in the browser`), and an unknown name; nothing is emitted.
#[test]
fn a_tool_run_runs_windowmasker_alone() {
    emitted();
    let refused = |name: &str| Err(format!("{name} is not available in the browser (ABI v2)"));
    assert_eq!(run(&["dustmasker", "-in", "x"]), refused("dustmasker"));
    assert_eq!(run(&["blastn"]), refused("blastn"));
    assert_eq!(validate(&["blastn"]), refused("blastn"));
    assert_eq!(
        run(&["nosuch"]),
        Err("unknown program 'nosuch'".to_string())
    );
    assert_eq!(run(&[]), Err("unknown program ''".to_string()));
    for &program in PROGRAMS.iter().filter(|program| !TOOLS.contains(program)) {
        assert_eq!(run(&[program]), refused(program));
    }
    assert!(emitted().is_empty());
}

/// Tool runs and tool validation read `-logfile` against the run's real files (line WQA D19):
/// missing/existing folders keep diagnostics on stream 3 with no log file; existing files are
/// refused (exit status 1 in `run`, `Err` in `validate`); `/dev/stdout` writes diagnostics into
/// stream 11 in chronological order; and stream 12 carries normalized host file names for both
/// `-logfile` and `-out`.
#[test]
fn a_tool_run_reads_logfile_against_the_runs_files() {
    emitted();
    let fasta = fixture("winmask/fasta/dm2L_3k.fa");
    // Case 0 & 1: missing folder and existing folder (`src`) keep diagnostics on stream 3.
    for log in ["d19-no-such-folder/run.log", "src", "src/"] {
        run(&[
            "windowmasker",
            "-mk_counts",
            "-in",
            &fasta,
            "-unit",
            "5",
            "-out",
            "d19-tool.counts",
            "-logfile",
            log,
        ])
        .unwrap();
        let calls = emitted();
        assert_eq!(
            calls.last().unwrap(),
            &(STREAM_RESPONSE, b"{\"status\":0}".to_vec()),
            "{log}"
        );
        assert_eq!(
            joined(&calls, STREAM_DIAGNOSTICS),
            b"pass 1\npass 2\n",
            "{log}"
        );
        assert_eq!(
            files(&frames(&calls))
                .iter()
                .map(|(n, _)| n.as_slice())
                .collect::<Vec<_>>(),
            vec![b"d19-tool.counts".as_slice()],
            "{log}"
        );
    }
    // Case 2: existing file of the run's directory (`Cargo.toml` in `web/adapter`) is refused by
    // both `run` (status 1, stream 3) and `validate`.
    let expected_refusal = "Error: a -logfile file 'Cargo.toml' that is a file of the run's \
         directory (LOSAT Web reads the run's files and writes none of them) is not supported by \
         LOSAT's WINDOWMASKER";
    assert_eq!(
        validate(&[
            "windowmasker",
            "-mk_counts",
            "-in",
            &fasta,
            "-logfile",
            "Cargo.toml"
        ]),
        Err(expected_refusal.to_string())
    );
    run(&[
        "windowmasker",
        "-mk_counts",
        "-in",
        &fasta,
        "-logfile",
        "Cargo.toml",
    ])
    .unwrap();
    let calls = emitted();
    assert_eq!(
        calls.last().unwrap(),
        &(STREAM_RESPONSE, b"{\"status\":1}".to_vec())
    );
    assert_eq!(
        joined(&calls, STREAM_DIAGNOSTICS),
        format!("{expected_refusal}\n").as_bytes()
    );
    assert!(joined(&calls, STREAM_PRINTED).is_empty());
    assert!(frames(&calls).is_empty());
    // Case 3: `-logfile /dev/stdout` writes diagnostic messages into stream 11 before counts.
    run(&["windowmasker", "-mk_counts", "-in", &fasta, "-unit", "5"]).unwrap();
    let plain_calls = emitted();
    let plain_out = joined(&plain_calls, STREAM_PRINTED);
    assert_eq!(
        joined(&plain_calls, STREAM_DIAGNOSTICS),
        b"pass 1\npass 2\n"
    );
    for extra in [&[][..], &["-out", "/dev/stdout"][..]] {
        let mut words = vec!["windowmasker", "-mk_counts", "-in", &fasta, "-unit", "5"];
        words.extend_from_slice(extra);
        words.extend(["-logfile", "/dev/stdout"]);
        run(&words).unwrap();
        let calls = emitted();
        assert_eq!(
            calls.last().unwrap(),
            &(STREAM_RESPONSE, b"{\"status\":0}".to_vec())
        );
        assert!(joined(&calls, STREAM_DIAGNOSTICS).is_empty(), "{extra:?}");
        assert!(frames(&calls).is_empty(), "{extra:?}");
        let mut want = b"pass 1\npass 2\n".to_vec();
        want.extend_from_slice(&plain_out);
        assert_eq!(joined(&calls, STREAM_PRINTED), want, "{extra:?}");
    }
    // Case 3b: binary `-out - -logfile /dev/stdout` (`d19.t3_stdout_binary`): NCBI's `cout`
    // buffer holds the 4-byte binary stream prefix while `-logfile /dev/stdout` writes `pass 1`
    // and `pass 2`, so the log lines precede the binary prefix just as on the CLI.
    run(&[
        "windowmasker",
        "-mk_counts",
        "-in",
        &fasta,
        "-unit",
        "5",
        "-sformat",
        "binary",
        "-out",
        "-",
    ])
    .unwrap();
    let plain_bin = joined(&emitted(), STREAM_PRINTED);
    assert!(plain_bin.starts_with(&[0, 0, 0, 0]));
    run(&[
        "windowmasker",
        "-mk_counts",
        "-in",
        &fasta,
        "-unit",
        "5",
        "-sformat",
        "binary",
        "-out",
        "-",
        "-logfile",
        "/dev/stdout",
    ])
    .unwrap();
    let bin_calls = emitted();
    assert_eq!(
        bin_calls.last().unwrap(),
        &(STREAM_RESPONSE, b"{\"status\":0}".to_vec())
    );
    assert!(joined(&bin_calls, STREAM_DIAGNOSTICS).is_empty());
    let mut want_bin = b"pass 1\npass 2\n".to_vec();
    want_bin.extend_from_slice(&plain_bin);
    assert_eq!(joined(&bin_calls, STREAM_PRINTED), want_bin);
    // On a missing-input error with `-sformat binary -out - -logfile /dev/stdout`, the buffered
    // 4-byte prefix flushes at exit after the error posts, matching CLI/NCBI.
    run(&[
        "windowmasker",
        "-mk_counts",
        "-in",
        "no-such-input-d19-bin.fa",
        "-unit",
        "5",
        "-sformat",
        "binary",
        "-out",
        "-",
        "-logfile",
        "/dev/stdout",
    ])
    .unwrap();
    let err_bin_calls = emitted();
    let err_bin_out = joined(&err_bin_calls, STREAM_PRINTED);
    assert!(err_bin_out.starts_with(b"pass 1\nError:"));
    assert!(err_bin_out.ends_with(&[0, 0, 0, 0]));
    // Case 4: normalized file names on stream 12 (` run.log` -> `run.log`, `./x.counts` -> `x.counts`,
    // `d19-slash/` -> `d19-slash`).
    let file_names = |calls: &[(u32, Vec<u8>)]| {
        files(&frames(calls))
            .into_iter()
            .map(|(n, _)| n)
            .collect::<Vec<_>>()
    };
    run(&[
        "windowmasker",
        "-mk_counts",
        "-in",
        &fasta,
        "-unit",
        "5",
        "-out",
        "./x.counts",
        "-logfile",
        " run.log",
    ])
    .unwrap();
    let calls = emitted();
    assert!(joined(&calls, STREAM_DIAGNOSTICS).is_empty(), "run.log");
    let got = files(&frames(&calls));
    assert_eq!(
        file_names(&calls),
        vec![b"run.log".to_vec(), b"x.counts".to_vec()]
    );
    assert_eq!(got[0].1, b"pass 1\npass 2\n");
    run(&[
        "windowmasker",
        "-mk_counts",
        "-in",
        &fasta,
        "-unit",
        "5",
        "-out",
        "./x.counts",
        "-logfile",
        "d19-slash/",
    ])
    .unwrap();
    let calls = emitted();
    assert_eq!(
        file_names(&calls),
        vec![b"d19-slash".to_vec(), b"x.counts".to_vec()]
    );

    // Round R2 finding 2: `-logfile` or `-out` naming a missing input the CLI reads back is
    // refused by both `validate` and `run` (status 1, same sentence, no file).
    let log_refusal = |name: &str| {
        format!(
            "Error: a -logfile file '{name}' that is a file of the run's directory (LOSAT Web \
             reads the run's files and writes none of them) is not supported by LOSAT's \
             WINDOWMASKER"
        )
    };
    let out_refusal = |name: &str| {
        format!(
            "Error: an output file '{name}' that is a file of the run's directory (the command \
             line empties it; LOSAT Web reads the run's files and writes none of them) is not \
             supported by LOSAT's WINDOWMASKER"
        )
    };
    for (label, words, want) in [
        (
            "mk_counts -in missing -logfile same",
            vec![
                "windowmasker",
                "-mk_counts",
                "-in",
                "d19-r2-missing.fa",
                "-logfile",
                "d19-r2-missing.fa",
                "-out",
                "d19-r2.counts",
            ],
            log_refusal("d19-r2-missing.fa"),
        ),
        (
            "ustat missing -logfile same",
            vec![
                "windowmasker",
                "-ustat",
                "d19-r2-missing.counts",
                "-in",
                &fasta,
                "-logfile",
                "d19-r2-missing.counts",
                "-out",
                "d19-r2.mask",
            ],
            log_refusal("d19-r2-missing.counts"),
        ),
        (
            "mk_counts -in missing -out same",
            vec![
                "windowmasker",
                "-mk_counts",
                "-in",
                "d19-r2-missing-x.fa",
                "-out",
                "d19-r2-missing-x.fa",
            ],
            out_refusal("d19-r2-missing-x.fa"),
        ),
        (
            "ustat missing -out same",
            vec![
                "windowmasker",
                "-ustat",
                "d19-r2-missing.counts",
                "-in",
                &fasta,
                "-out",
                "d19-r2-missing.counts",
            ],
            out_refusal("d19-r2-missing.counts"),
        ),
    ] {
        assert_eq!(validate(&words), Err(want.clone()), "{label}: validate");
        run(&words).unwrap();
        let calls = emitted();
        assert_eq!(
            calls.last().unwrap(),
            &(STREAM_RESPONSE, b"{\"status\":1}".to_vec()),
            "{label}: status"
        );
        assert_eq!(
            joined(&calls, STREAM_DIAGNOSTICS),
            format!("{want}\n").as_bytes(),
            "{label}: diag"
        );
        assert!(joined(&calls, STREAM_PRINTED).is_empty(), "{label}: stdout");
        assert!(frames(&calls).is_empty(), "{label}: no files");
    }
    // Controls: `-mk_counts -in x.fa` (missing, no write name) stays status 3 with `Read error`,
    // and a missing input's parent folder must NOT become present (`-in sub/x.fa -logfile
    // sub/run.log` without `sub` keeps diagnostics on stream 3 with no file).
    for (label, words) in [
        (
            "missing input control",
            &["windowmasker", "-mk_counts", "-in", "d19-r2-missing-x.fa"][..],
        ),
        (
            "missing input parent folder",
            &[
                "windowmasker",
                "-mk_counts",
                "-in",
                "d19-r2-nosub/x.fa",
                "-logfile",
                "d19-r2-nosub/run.log",
            ][..],
        ),
    ] {
        assert_eq!(validate(words), Ok(()), "{label}: validate");
        run(words).unwrap();
        let calls = emitted();
        assert_eq!(
            calls.last().unwrap(),
            &(STREAM_RESPONSE, b"{\"status\":3}".to_vec()),
            "{label}: status"
        );
        assert!(
            String::from_utf8_lossy(&joined(&calls, STREAM_DIAGNOSTICS)).contains("Read error"),
            "{label}: Read error"
        );
        assert!(frames(&calls).is_empty(), "{label}: no stream 12 files");
    }
}
