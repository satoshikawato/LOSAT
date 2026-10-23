#![allow(warnings, clippy::all)]

//! The blastp search in blocks of whole query batches (`LOSAT_X_BPQBLOCK`) for a `run_local` call
//! with several outputs, the hit records and an observer (`LOSAT_X_BPQBLOCKWEB`, ABI v2 of the Web
//! build; `src/algorithm/blastp/x_query_blocks.rs`).
//!
//! The switches are read once per process, so every search runs in a child process: this test
//! binary again, running only `child_search` with the variables of its arm (`BATCH_SIZE` makes
//! small batches and `LOSAT_X_BPQBLOCK=1` closes a block at every batch end, so the fixture pair has
//! three blocks or more). The child writes what the call wrote (each output, the warnings, the hit
//! records as `Debug` text and the observer's ranges) to files and checks the observer's ranges
//! against the outputs and the records (`assert_observer_ranges`, with each record's `q_idx`
//! against the `qseqid` of its row). The parent compares the block search with the search of all
//! queries (`LOSAT_X_BPQBLOCKWEB=0`) and with the CLI run of each format, and counts the
//! `[TIMING] blastp_total` lines (one per search of a query set) to tell whether the blocks ran.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-295
//! ```c++
//!         formatter.PrintProlog();
//!         ...
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!             ...
//!                     formatter.PrintOneResultSet(**result, query_batch);
//!         ...
//!         formatter.PrintEpilog(opt);
//! ```

mod run_local_support;

use std::path::{Path, PathBuf};
use std::process::Command;

use run_local_support::{assert_observer_ranges, field, reader_records, run_formats, temp_path};
use LOSAT::api::local_blast::{
    run_local_blastp, FormatObserver, FormatOutput, HspIndex, OutputSink, ReportOutputs,
};
use LOSAT::blastinput::fasta_reader::ReaderConfig;
use LOSAT::cli::{try_parse_from, Cli, Commands};

const CUSTOM: &str =
    "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq btop";
/// The outputs of ABI v2 of the Web build.
const WEB: [&str; 3] = ["0", "6", "7"];
/// Every format the block search takes, outfmt 0, 6 and 7 first.
const EIGHT: [&str; 8] = ["0", "6", "7", CUSTOM, "8", "9", "12", "5"];
const BATCH: (&str, &str) = ("BATCH_SIZE", "3000");
const BLOCK_PER_BATCH: (&str, &str) = ("LOSAT_X_BPQBLOCK", "1");
const WEB_OFF: (&str, &str) = ("LOSAT_X_BPQBLOCKWEB", "0");

/// The child's instructions (set only in a child process).
const CHILD: &str = "QBLOCK_WEB_CHILD";

/// The first `count` records of a fixture, byte for byte, in a temporary file.
fn fixture_prefix(name: &str, count: usize) -> PathBuf {
    let source = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    let text = std::fs::read_to_string(source).expect("fixture FASTA");
    let mut records = 0;
    let mut kept = String::new();
    for line in text.split_inclusive('\n') {
        if line.starts_with('>') {
            records += 1;
            if records > count {
                break;
            }
        }
        kept.push_str(line);
    }
    let path = temp_path(&format!("qblockweb_{count}_{name}"));
    std::fs::write(&path, kept).expect("write fixture prefix");
    path
}

/// A query and a subject file: 20 SicyWSV proteins (11,114 letters, four batches of `BATCH` and so
/// four blocks; 13 of them with hits, from the second to the last) against the 86 PajaWSV proteins.
struct Inputs {
    query: PathBuf,
    subject: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        Self {
            query: fixture_prefix("SicyWSV.faa", 20),
            subject: Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/PajaWSV.faa"),
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.query);
    }
}

/// What the call gives the caller besides the outputs.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Observe {
    /// No hit records, no observer.
    Nothing,
    /// The hit records and the observer.
    Both,
    /// The observer alone.
    ObserverOnly,
    /// Nothing, each output a file the call creates.
    Files,
}

impl Observe {
    fn name(self) -> &'static str {
        match self {
            Observe::Nothing => "nothing",
            Observe::Both => "both",
            Observe::ObserverOnly => "observer",
            Observe::Files => "files",
        }
    }

    fn parse(name: &str) -> Self {
        match name {
            "nothing" => Observe::Nothing,
            "both" => Observe::Both,
            "observer" => Observe::ObserverOnly,
            "files" => Observe::Files,
            _ => panic!("observe {name}"),
        }
    }
}

/// One child's search.
struct Child {
    dir: PathBuf,
    /// The `[TIMING] blastp_total` lines: one per search of a query set.
    searches: usize,
    /// The `[X_SHADOW] BPQBLOCKFMT` lines.
    shadow: Vec<String>,
}

impl Child {
    fn file(&self, name: &str) -> Vec<u8> {
        std::fs::read(self.dir.join(name)).unwrap_or_else(|error| panic!("{name}: {error}"))
    }
}

impl Drop for Child {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The environment of a child or CLI process: no `LOSAT_*` variable that chooses a path and no
/// batch size of the test's environment, `LOSAT_BLASTX_WORKER_LOG` kept or set (a binary built with
/// the `blastx-worker-probe` feature, `cargo test --all-features` as in CI, needs it at 2 or more
/// threads; it selects nothing), then `vars`.
fn clean_env(command: &mut Command, vars: &[(&str, &str)]) {
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(
                &*text,
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | CHILD
            )
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join(format!("losat-qblockweb-worker-{}.log", std::process::id())),
        );
    }
    command.envs(vars.iter().copied());
}

/// Runs `child_search` in a child process with `vars` (and `LOSAT_TIMING`).
fn child(
    inputs: &Inputs,
    formats: &[&str],
    observe: Observe,
    threads: u32,
    vars: &[(&str, &str)],
    tag: &str,
) -> Child {
    let dir = temp_path(&format!("qblockweb_{tag}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("child directory");
    let spec = [
        dir.display().to_string(),
        inputs.query.display().to_string(),
        inputs.subject.display().to_string(),
        threads.to_string(),
        observe.name().to_string(),
    ]
    .into_iter()
    .chain(formats.iter().map(|format| format.to_string()))
    .collect::<Vec<_>>()
    .join("\n");
    let mut command = Command::new(std::env::current_exe().expect("test binary"));
    clean_env(&mut command, vars);
    let output = command
        .args([
            "child_search",
            "--exact",
            "--nocapture",
            "--test-threads",
            "1",
        ])
        .env(CHILD, spec)
        .env("LOSAT_TIMING", "1")
        .output()
        .expect("run the child");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "child {tag} failed: {}\n{stderr}",
        String::from_utf8_lossy(&output.stdout)
    );
    assert!(
        dir.join("done").exists(),
        "child {tag} did not search: {}",
        String::from_utf8_lossy(&output.stdout)
    );
    Child {
        dir,
        searches: stderr
            .lines()
            .filter(|line| line.starts_with("[TIMING] blastp_total:"))
            .count(),
        shadow: stderr
            .lines()
            .filter(|line| line.starts_with("[X_SHADOW] BPQBLOCKFMT"))
            .map(str::to_string)
            .collect(),
    }
}

/// The CLI output file and standard error of one single-format run with `vars`.
fn cli(inputs: &Inputs, outfmt: &str, threads: u32, vars: &[(&str, &str)]) -> (Vec<u8>, Vec<u8>) {
    let out = temp_path(&format!("qblockweb_cli_{}_{threads}.out", outfmt.len()));
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    clean_env(&mut command, vars);
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(&inputs.query)
        .arg("-subject")
        .arg(&inputs.subject)
        .args(["-evalue", "0.001", "-num_threads", &threads.to_string()])
        .arg("-outfmt")
        .arg(outfmt)
        .arg("-out")
        .arg(&out)
        .output()
        .expect("run LOSAT CLI");
    assert!(output.status.success(), "CLI failed: {output:?}");
    let bytes = std::fs::read(&out).expect("CLI output");
    let _ = std::fs::remove_file(&out);
    (bytes, output.stderr)
}

/// Lines that do not start with `[TIMING]` (the child's warnings sink gets none; kept defensive).
fn without_timing(bytes: &[u8]) -> Vec<u8> {
    bytes
        .split_inclusive(|&byte| byte == b'\n')
        .filter(|line| !line.starts_with(b"[TIMING]"))
        .flatten()
        .copied()
        .collect()
}

/// Counts the events of an observer that comes without the hit records.
struct CountingObserver(usize);

impl FormatObserver for CountingObserver {
    fn hsp_begin(&mut self, _format: usize, _hsp: HspIndex) {
        self.0 += 1;
    }
    fn hsp_end(&mut self, _format: usize, _hsp: HspIndex) {
        self.0 += 1;
    }
}

/// The search of one child process (`CHILD` holds the directory to write to, the query and subject
/// files, the threads, what to observe and the formats, one per line); without `CHILD` it does
/// nothing.
#[test]
fn child_search() {
    let Ok(spec) = std::env::var(CHILD) else {
        return;
    };
    let lines: Vec<&str> = spec.lines().collect();
    let [dir, query, subject, threads, observe, formats @ ..] = &lines[..] else {
        panic!("child spec {spec:?}");
    };
    let dir = Path::new(dir);
    let observe = Observe::parse(observe);
    let queries = reader_records(Path::new(query), ReaderConfig::query("BLASTP", true, true));
    let subjects = reader_records(
        Path::new(subject),
        ReaderConfig::subject("BLASTP", true, true),
    );
    let args = || {
        let argv = [
            "LOSAT",
            "blastp",
            "-query",
            query,
            "-subject",
            subject,
            "-evalue",
            "0.001",
            "-num_threads",
            threads,
        ];
        let cli: Cli = try_parse_from(argv).expect("blastp argv");
        match cli.command {
            Commands::Blastp(args) => args,
            _ => unreachable!("blastp subcommand"),
        }
    };
    let write = |name: &str, bytes: &[u8]| std::fs::write(dir.join(name), bytes).expect(name);
    if observe == Observe::Files {
        let paths: Vec<PathBuf> = (0..formats.len())
            .map(|index| dir.join(format!("out{index}")))
            .collect();
        let mut diagnostics = Vec::new();
        {
            let mut outputs = ReportOutputs {
                formats: formats
                    .iter()
                    .zip(&paths)
                    .map(|(&outfmt, path)| FormatOutput {
                        outfmt,
                        sink: OutputSink::File(path),
                    })
                    .collect(),
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            run_local_blastp(args(), &queries, &subjects, "", "", &mut outputs).expect("run_local");
        }
        write("diagnostics", &diagnostics);
        write("done", b"");
        return;
    }
    if observe == Observe::ObserverOnly {
        let mut sinks: Vec<Vec<u8>> = formats.iter().map(|_| Vec::new()).collect();
        let mut diagnostics = Vec::new();
        let mut observer = CountingObserver(0);
        {
            let mut outputs = ReportOutputs {
                formats: formats
                    .iter()
                    .zip(sinks.iter_mut())
                    .map(|(&outfmt, sink)| FormatOutput {
                        outfmt,
                        sink: OutputSink::Writer(sink),
                    })
                    .collect(),
                diagnostics: &mut diagnostics,
                hits: None,
                observer: Some(&mut observer),
            };
            run_local_blastp(args(), &queries, &subjects, "", "", &mut outputs).expect("run_local");
        }
        assert!(observer.0 > 0, "the observer saw no HSP");
        for (index, bytes) in sinks.iter().enumerate() {
            write(&format!("out{index}"), bytes);
        }
        write("diagnostics", &diagnostics);
        write("done", b"");
        return;
    }
    let run = run_formats(formats, observe == Observe::Both, |outputs| {
        run_local_blastp(args(), &queries, &subjects, "", "", outputs)
    });
    if observe == Observe::Both && formats.starts_with(&WEB) {
        assert_observer_ranges(&run, [0, 1, 2], &format!("blastp {formats:?} n{threads}"));
        // Each record's query is the query of its outfmt 6 row (the first word of the fixture's
        // definition lines is the `qseqid`).
        for &(hsp, start, end) in &run.ranges[1] {
            let q_idx = run.hits[hsp].hit.q_idx as usize;
            let title = String::from_utf8_lossy(&queries[q_idx].title).into_owned();
            assert_eq!(
                field(&run.outputs[1][start..end], 0),
                title.split_whitespace().next().expect("an id"),
                "HSP {hsp}: q_idx {q_idx}"
            );
        }
        let last = run.hits.last().expect("hits").hit.q_idx as usize;
        assert!(
            last >= queries.len() / 2,
            "the hits must reach the later blocks (last q_idx {last} of {})",
            queries.len()
        );
    }
    for (index, bytes) in run.outputs.iter().enumerate() {
        write(&format!("out{index}"), bytes);
    }
    write("diagnostics", &run.diagnostics);
    write(
        "hits",
        run.hits
            .iter()
            .map(|hit| format!("{hit:?}\n"))
            .collect::<String>()
            .as_bytes(),
    );
    write(
        "ranges",
        format!("{:?}\n{:?}\n", run.ranges, run.subjects).as_bytes(),
    );
    write("done", b"");
}

/// The block search of outfmt 0, 6 and 7 (ABI v2 of the Web build) and of eight formats, with the
/// hit records and the observer, at 1 and 4 threads: the same outputs, warnings, records and
/// ranges as the search of all queries (`LOSAT_X_BPQBLOCKWEB=0`), and each output the CLI run of
/// its format.
#[test]
fn outputs_records_and_ranges_of_the_blocks_are_those_of_all_queries() {
    let inputs = Inputs::new();
    for threads in [1, 4] {
        for formats in [&WEB[..], &EIGHT[..]] {
            let tag = format!("{}_{threads}", formats.len());
            let blocks = child(
                &inputs,
                formats,
                Observe::Both,
                threads,
                &[BATCH, BLOCK_PER_BATCH],
                &format!("blocks_{tag}"),
            );
            let all = child(
                &inputs,
                formats,
                Observe::Both,
                threads,
                &[BATCH, BLOCK_PER_BATCH, WEB_OFF],
                &format!("all_{tag}"),
            );
            assert!(
                blocks.searches >= 3,
                "{tag}: three blocks or more, got {} searches",
                blocks.searches
            );
            assert_eq!(all.searches, 1, "{tag}: LOSAT_X_BPQBLOCKWEB=0");
            for name in ["diagnostics", "hits", "ranges"]
                .into_iter()
                .map(str::to_string)
                .chain((0..formats.len()).map(|index| format!("out{index}")))
            {
                assert!(
                    blocks.file(&name) == all.file(&name),
                    "{tag}: {name} of the blocks differs from the search of all queries"
                );
            }
            let mut first_stderr = None;
            for (index, outfmt) in formats.iter().enumerate() {
                let (output, stderr) = cli(
                    &inputs,
                    outfmt,
                    threads,
                    &[BATCH, ("LOSAT_X_BPQBLOCK", "0")],
                );
                assert!(
                    blocks.file(&format!("out{index}")) == output,
                    "{tag}: outfmt {outfmt} differs from the CLI"
                );
                first_stderr.get_or_insert(stderr);
            }
            assert_eq!(
                String::from_utf8_lossy(&without_timing(&blocks.file("diagnostics"))),
                String::from_utf8_lossy(&first_stderr.expect("a format")),
                "{tag}: the warnings are the CLI's of the first format"
            );
        }
    }
}

/// The shadow (`LOSAT_X_BPQBLOCKFMTSHADOW`) compares the outputs, the warnings, the hit records and
/// each format's observer events of the blocks with the search of all queries inside the run.
#[test]
fn the_shadow_compares_outputs_records_and_events() {
    let inputs = Inputs::new();
    for threads in [1, 4] {
        let run = child(
            &inputs,
            &WEB,
            Observe::Both,
            threads,
            &[BATCH, BLOCK_PER_BATCH, ("LOSAT_X_BPQBLOCKFMTSHADOW", "1")],
            &format!("shadow_{threads}"),
        );
        assert_eq!(run.shadow.len(), 1, "n{threads}: {:?}", run.shadow);
        let line = &run.shadow[0];
        let value = |key: &str| -> usize {
            line.split(' ')
                .find_map(|word| word.strip_prefix(key))
                .unwrap_or_else(|| panic!("{key} in {line}"))
                .parse()
                .expect("a count")
        };
        assert!(line.ends_with(" identical"), "n{threads}: {line}");
        assert!(value("blocks=") >= 3, "n{threads}: {line}");
        assert_eq!(value("formats="), 3, "n{threads}: {line}");
        let hits = value("hits=");
        assert!(hits > 0, "n{threads}: {line}");
        // Each HSP opens and closes once in each of the three formats, and outfmt 0 heads each
        // subject of each query once.
        assert!(value("events=") > 6 * hits, "n{threads}: {line}");
    }
}

/// The calls the block search does not take (the search of all queries, one search) and those it
/// takes with one output.
#[test]
fn calls_the_block_search_takes_and_does_not() {
    let inputs = Inputs::new();
    let vars = [BATCH, BLOCK_PER_BATCH];
    // One output with the hit records and the observer, and several outputs alone.
    for (formats, observe) in [(&WEB[..1], Observe::Both), (&WEB[..], Observe::Nothing)] {
        let run = child(&inputs, formats, observe, 1, &vars, "takes");
        assert!(
            run.searches >= 3,
            "{formats:?} {observe:?}: {}",
            run.searches
        );
    }
    // An observer without the hit records; several outputs that are not the caller's writers; an
    // output the block search does not write (outfmt 15 and 16); the switch off
    // (`LOSAT_X_BPQBLOCKWEB=0`; one outfmt 6 output with the records and the observer, which only
    // `LOSAT_X_BPQBLOCKWEB` takes); the block search off.
    for (formats, observe, extra) in [
        (&WEB[..], Observe::ObserverOnly, None),
        (&WEB[..], Observe::Files, None),
        (&["0", "6", "7", "15", "16"][..], Observe::Both, None),
        (&WEB[1..2], Observe::Both, Some(WEB_OFF)),
        (&WEB[..], Observe::Both, Some(("LOSAT_X_BPQBLOCK", "0"))),
    ] {
        let mut vars = vec![BATCH];
        if extra != Some(("LOSAT_X_BPQBLOCK", "0")) {
            vars.push(BLOCK_PER_BATCH);
        }
        vars.extend(extra);
        let run = child(&inputs, formats, observe, 1, &vars, "declines");
        assert_eq!(run.searches, 1, "{formats:?} {observe:?} {extra:?}");
    }
}

/// A call of two outputs (outfmt 6 and 16, both writers of the caller) over queries that make four
/// blocks does not take the block search, whichever of `LOSAT_X_BPQBLOCKXML2` and
/// `LOSAT_X_BPQBLOCKWEB` is on: one search, and the bytes are those of the run with the block
/// search off. (The outfmt 16 run-level frame was checked for a call of one output only.) The CLI's
/// one-output outfmt 16 run takes the blocks and writes the same document.
#[test]
fn outfmt_16_in_a_call_of_several_outputs_is_not_blocked() {
    let inputs = Inputs::new();
    let formats = ["6", "16"];
    let unblocked = child(
        &inputs,
        &formats,
        Observe::Nothing,
        1,
        &[BATCH, ("LOSAT_X_BPQBLOCK", "0")],
        "x16_off",
    );
    let run = child(
        &inputs,
        &formats,
        Observe::Nothing,
        1,
        &[BATCH, BLOCK_PER_BATCH],
        "x16_on",
    );
    assert_eq!(run.searches, 1, "outfmt 6 and 16: one search");
    assert_eq!(unblocked.searches, 1);
    for name in ["out0", "out1", "diagnostics"] {
        assert!(
            run.file(name) == unblocked.file(name),
            "{name} differs from the unblocked run"
        );
    }
    // The CLI's one-output outfmt 16 run is blocked (four batches, one block each) and writes the
    // same document as the second output.
    let (document, _) = cli(&inputs, "16", 1, &[BATCH, BLOCK_PER_BATCH]);
    assert!(
        run.file("out1") == document,
        "outfmt 16 differs from the CLI"
    );
}
