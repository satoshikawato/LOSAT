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

/// The files of the frames, in the order of their beginnings: the name and the bytes; every file
/// begins once, then has bytes, then ends once.
fn files(frames: &[(u32, u32, &[u8])]) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut files: Vec<(Vec<u8>, Vec<u8>, bool)> = Vec::new();
    for &(number, kind, payload) in frames {
        let number = number as usize;
        match kind {
            FRAME_BEGIN => {
                assert_eq!(number, files.len(), "numbered in order");
                files.push((payload.to_vec(), Vec::new(), false));
            }
            FRAME_BYTES => {
                assert!(!files[number].2, "bytes after the end");
                files[number].1.extend_from_slice(payload);
            }
            FRAME_END => {
                assert!(payload.is_empty() && !files[number].2, "one end");
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
/// masked FASTA of a sequence of 1.2 Mbp comes in two frames of bytes or more, all of them before
/// the end of the file of `-logfile`, which the run closes when it ends; its bytes are the
/// library's on the disk.
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
    let log_end = frames
        .iter()
        .position(|&(number, kind, _)| number == 0 && kind == FRAME_END)
        .unwrap();
    assert_eq!(log_end, frames.len() - 1, "the log ends with the run");
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

/// The validation of a tool's argv (D-WW3-7): the command line's checks before the tool opens a
/// file pass a good argv (also one of a missing input, which the run reads), and refuse the others
/// with the command line's error: its standard error, without the last line end; nothing emitted.
#[test]
fn a_tool_s_argv_is_checked_as_the_command_line_checks_it() {
    emitted();
    assert_eq!(
        validate(&["windowmasker", "-mk_counts", "-in", "missing.fa"]),
        Ok(())
    );
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
    assert!(emitted().is_empty());
}
