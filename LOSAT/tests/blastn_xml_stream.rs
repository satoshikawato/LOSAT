//! BLAST XML (outfmt 5) and the XML2 and JSON documents (16 and 15 on the output, 14 and 13 in side
//! files) of blastn written batch by batch (`LOSAT_X_RPTSTREAM`,
//! `src/algorithm/blastn/blast_engine/x_rpt_stream.rs`; line RG series D, remaining-gaps DECISIONS
//! D-SO-2), as NCBI's formatter writes each batch's reports before it reads the next, the header
//! or start with the first report and the end in the epilog.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-748
//!
//! Every run that LOSAT does not refuse writes the same bytes with the switch (on by default) and
//! without it (`LOSAT_X_RPTSTREAM=0`, the reports after the loop): the exit status, standard output
//! or the `-out` file, the side files and standard error, for megablast and blastn, `-subject` and
//! `-db`, at 1 and 4 threads, with small batches (`BATCH_SIZE`), also where NCBI's reader error
//! stops a later batch; the `[X_STATS] RPTSTREAM` line shows that the loop wrote the batches. A
//! title that LOSAT refuses (one NCBI crashes on, PD-LOSAT-NCBI-DEFECTS) comes after the output
//! that NCBI has written by then (the Owner's decision of 2026-11-16 on line RG's Q-X2c, "the
//! results of the query batches before are written"; remaining-gaps DECISIONS D-SO-6): the start of
//! the document of the same search without the refused title, without its end. BLAST XML refuses
//! the refused title's batch before its first byte, after the batches before; the held report of a
//! `-subject` XML2 or JSON search is refused where NCBI makes it (the arrival of the next query's
//! result, or the epilog), after the reports before it, in any layout of the batches and with the
//! switch or without it alike, and is dropped with its refusal where an error of NCBI's stops the
//! run first, as NCBI never makes it. A refused first report writes nothing. That the bytes equal
//! NCBI's is the business of the blastn fixtures (`tests/db_search_manifest/blastn.dxs.tsv`).

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

mod temp_support;
use temp_support::{worker_log_if_unset, TempFile};

const OFF: (&str, &str) = ("LOSAT_X_RPTSTREAM", "0");
const STATS: (&str, &str) = ("LOSAT_X_STATS", "1");

/// The refusal of a title NCBI crashes on (`blast_xml::names::Refusal`), for any title.
const REFUSAL: &str =
    "with punctuation that NCBI BLAST+ 2.17.0 crashes on when it makes the definition line is not supported by LOSAT's BLASTN";

/// The title of punctuation that NCBI's definition line generator reads past the end of.
const PUNCTUATION: &[u8] = b", ,";

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The letters of the first record of a FASTA file.
fn letters(name: &str) -> Vec<u8> {
    std::fs::read(fasta(name))
        .expect("read FASTA")
        .split(|&byte| byte == b'\n')
        .skip(1)
        .flat_map(|line| {
            line.iter()
                .copied()
                .filter(|byte| !byte.is_ascii_whitespace())
        })
        .collect()
}

/// A directory of this test in the temporary directory, removed when dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat_xmlstream_{}_{}_{tag}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        TempDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A FASTA record with `title` and `seq` in lines of 61 letters.
fn record(text: &mut Vec<u8>, title: &[u8], seq: &[u8]) {
    text.push(b'>');
    text.extend_from_slice(title);
    text.push(b'\n');
    for line in seq.chunks(61) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
}

/// 16 pieces of LC738874 of 900 to 4,500 letters (every fourth with a run of `N` and lower case
/// letters), one of them with a line of digits (the reader's warning), then a record without
/// letters, one of `N`s only (the warning of an invalid query) and one of two letters: with a batch
/// size of 4,000 letters about a dozen batches.
fn queries(dir: &TempDir, name: &str) -> PathBuf {
    let genome = letters("LC738874.fasta");
    let mut text = Vec::new();
    for i in 0..16usize {
        let start = i * 11_000 + 137 * i;
        let len = 900 + (i * 1_531) % 3_600;
        let mut seq = genome[start..start + len].to_vec();
        if i % 4 == 0 {
            let at = len / 2;
            seq.splice(at..at, std::iter::repeat_n(b'N', 25));
            for byte in &mut seq[at + 40..at + 140] {
                *byte = byte.to_ascii_lowercase();
            }
        }
        record(
            &mut text,
            format!("piece{i} title words {i}").as_bytes(),
            &seq,
        );
        if i == 6 {
            // A line the reader keeps with a warning for its digits.
            text.extend_from_slice(b"ACGT12ACGT\n");
        }
        if i == 9 {
            text.extend_from_slice(b">without letters\n");
            record(&mut text, b"onlyN", &[b'N'; 600]);
        }
    }
    record(&mut text, b"two", b"AC");
    let path = dir.path(name);
    std::fs::write(&path, text).expect("write the queries");
    path
}

/// Six pieces of 400 letters of the subject LC738868 (each hits it), titled `p<n> words` except
/// the one at `refused` (the title NCBI crashes on), then `tail`: with a batch size of 400 letters
/// one query per batch.
fn pieces(dir: &TempDir, name: &str, refused: Option<usize>, tail: &[u8]) -> PathBuf {
    let genome = letters("LC738868.fasta");
    let mut text = Vec::new();
    for i in 0..6usize {
        let start = 3_000 + i * 17_000;
        let title = match refused {
            Some(at) if at == i => PUNCTUATION.to_vec(),
            _ => format!("p{i} words").into_bytes(),
        };
        record(&mut text, &title, &genome[start..start + 400]);
    }
    text.extend_from_slice(tail);
    let path = dir.path(name);
    std::fs::write(&path, text).expect("write the pieces");
    path
}

/// One run's exit status, report (standard output, or the `-out` file), side files (by name) and
/// standard error.
#[derive(Debug)]
struct Run {
    /// The run's worker log, when the test made one: kept until the caller drops the run.
    _log: Option<TempFile>,
    status: Option<i32>,
    report: Vec<u8>,
    files: BTreeMap<String, Vec<u8>>,
    stderr: Vec<u8>,
}

impl PartialEq for Run {
    fn eq(&self, other: &Run) -> bool {
        self.status == other.status
            && self.report == other.report
            && self.files == other.files
            && self.stderr == other.stderr
    }
}

/// One search with `args`, the given variables as its only `LOSAT_*` variables, `BATCH_SIZE`
/// `batch` (`CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` removed). `LOSAT_BLASTX_WORKER_LOG` stays (a
/// binary built with the `blastx-worker-probe` feature needs it once a search has more than one
/// thread; it selects nothing). With `out`, a file in a directory of its own (emptied first), the
/// report is that file's bytes (standard output must then be empty) and the other files of the
/// directory are the side files.
fn blastn(
    args: &[&str],
    threads: u32,
    batch: &str,
    switches: &[(&str, &str)],
    out: Option<&Path>,
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
            || text.starts_with("BLEMIR_")
        {
            command.env_remove(name);
        }
    }
    let log = worker_log_if_unset(&mut command, "blastn_xml_stream_worker");
    command
        .arg("blastn")
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .env("BATCH_SIZE", batch)
        .envs(switches.iter().copied());
    let out_dir = out.map(|out| out.parent().expect("the -out directory"));
    if let (Some(out), Some(dir)) = (out, out_dir) {
        let _ = std::fs::remove_dir_all(dir);
        std::fs::create_dir_all(dir).expect("create the -out directory");
        command.arg("-out").arg(out);
    }
    let output = command.output().expect("run LOSAT");
    let mut files = BTreeMap::new();
    let report = match (out, out_dir) {
        (Some(out), Some(dir)) => {
            assert!(output.stdout.is_empty(), "{args:?}: output beside -out");
            for entry in std::fs::read_dir(dir).expect("list the -out directory") {
                let path = entry.expect("a directory entry").path();
                if path != out {
                    let name = path.file_name().unwrap().to_string_lossy().into_owned();
                    files.insert(name, std::fs::read(&path).expect("read a side file"));
                }
            }
            std::fs::read(out).unwrap_or_default()
        }
        _ => output.stdout,
    };
    Run {
        _log: log,
        status: output.status.code(),
        report,
        files,
        stderr: output.stderr,
    }
}

/// The `[X_STATS] RPTSTREAM` line taken out of `run`'s standard error (with the other `[X_STATS]`
/// lines): its `batches=` and `off=`.
fn take_stats_lines(run: &mut Run) -> (u64, String) {
    let text = String::from_utf8_lossy(&run.stderr).into_owned();
    let lines: Vec<&str> = text
        .split_inclusive('\n')
        .filter(|line| line.starts_with("[X_STATS] RPTSTREAM "))
        .collect();
    assert_eq!(lines.len(), 1, "one RPTSTREAM line:\n{text}");
    let line = lines[0].trim_end();
    let field = |name: &str| -> String {
        line.split_whitespace()
            .find_map(|word| word.strip_prefix(name))
            .unwrap_or_else(|| panic!("no {name} in: {line}"))
            .to_string()
    };
    let numbers = (field("batches=").parse().unwrap(), field("off="));
    let rest: Vec<&str> = text
        .split_inclusive('\n')
        .filter(|line| !line.starts_with("[X_STATS] "))
        .collect();
    run.stderr = rest.concat().into_bytes();
    numbers
}

/// The run without the switch, with it (the default), and with it and `LOSAT_X_STATS`: the same
/// status, report, side files and standard error (the counters' lines taken out of the last).
/// Returns the off run and the counters (batches written by the loop, why the switch did not act:
/// `none` when it did).
fn same_two(
    what: &str,
    args: &[&str],
    threads: u32,
    batch: &str,
    out: Option<&Path>,
) -> (Run, (u64, String)) {
    let off = blastn(args, threads, batch, &[OFF], out);
    let on = blastn(args, threads, batch, &[], out);
    assert!(
        on == off,
        "{what} threads {threads}: on differs\n{on:?}\n{off:?}"
    );
    let mut counted = blastn(args, threads, batch, &[STATS], out);
    let numbers = take_stats_lines(&mut counted);
    assert!(
        counted == off,
        "{what} threads {threads}: on with statistics differs\n{counted:?}\n{off:?}"
    );
    (off, numbers)
}

/// The number of reports in a document of outfmt `outfmt` (5: its `Iteration`s).
fn reports_in(outfmt: &str, document: &[u8]) -> usize {
    let mark: &[u8] = match outfmt {
        "5" => b"</Iteration>",
        "16" => b"<Report>",
        "15" => b"\"report\": {",
        _ => unreachable!("{outfmt}"),
    };
    document
        .windows(mark.len())
        .filter(|window| *window == mark)
        .count()
}

#[test]
fn every_xml_format_gives_the_off_bytes() {
    let dir = TempDir::new("formats");
    let query = queries(&dir, "q.fa");
    let subject = fasta("LC738868.fasta");
    let db = dir.path("db").to_str().unwrap().to_string();
    let made = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makeblastdb")
        .arg("-in")
        .arg(&subject)
        .args(["-dbtype", "nucl", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(made.status.success(), "makeblastdb");
    let out = dir.path("out").join("o.out");
    let mut threads = [1, 4].into_iter().cycle();
    for task in ["megablast", "blastn"] {
        for (target, outfmt, with_out) in [
            (["-subject", subject.to_str().unwrap()], "5", false),
            (["-subject", subject.to_str().unwrap()], "16", false),
            (["-subject", subject.to_str().unwrap()], "15", false),
            (["-subject", subject.to_str().unwrap()], "5", true),
            (["-subject", subject.to_str().unwrap()], "14", true),
            (["-subject", subject.to_str().unwrap()], "13", true),
            (["-db", &db], "5", false),
            (["-db", &db], "16", false),
            (["-db", &db], "15", false),
            (["-db", &db], "14", true),
        ] {
            let mut args = vec![
                "-task",
                task,
                "-query",
                query.to_str().unwrap(),
                target[0],
                target[1],
                "-outfmt",
                outfmt,
            ];
            if task == "blastn" && outfmt == "5" && target[0] == "-subject" && !with_out {
                // The labels of the queries' parsed ids (`QueryScopes`, numbered from each
                // batch's first query).
                args.push("-parse_deflines");
            }
            let threads = threads.next().unwrap();
            let what = format!(
                "{task} {} outfmt {outfmt} -out {with_out} {args:?}",
                target[0]
            );
            let (off, (batches, reason)) =
                same_two(&what, &args, threads, "4000", with_out.then_some(&*out));
            assert_eq!(off.status, Some(0), "{what}: {off:?}");
            assert!(!off.report.is_empty(), "{what}");
            if matches!(outfmt, "13" | "14") {
                assert!(off.files.len() >= 8, "{what}: {:?}", off.files.keys());
            }
            assert!(
                reason == "none" && batches >= 8,
                "{what}: {batches} {reason}"
            );
        }
    }
}

#[test]
fn a_refused_title_comes_after_what_ncbi_wrote_before_it_in_any_layout() {
    let dir = TempDir::new("refused");
    let subject = fasta("LC738868.fasta");
    for outfmt in ["5", "16", "15"] {
        let args_of = |query: &Path| -> Vec<String> {
            [
                "-task",
                "megablast",
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                outfmt,
            ]
            .map(String::from)
            .to_vec()
        };
        let whole_query = pieces(&dir, "whole.fa", None, b"");
        let whole_args = args_of(&whole_query);
        let whole_args: Vec<&str> = whole_args.iter().map(String::as_str).collect();
        let whole = blastn(&whole_args, 1, "400", &[], None);
        assert_eq!(whole.status, Some(0), "outfmt {outfmt}: {whole:?}");
        assert_eq!(reports_in(outfmt, &whole.report), 6, "outfmt {outfmt}");
        // The title on the first query, in the third, the fourth (the second of a batch of two)
        // and the last (the held report of 16 and 15 is then refused by the epilog): a batch of
        // each query (400 letters), of two (800) and one batch.
        for (refused, threads) in [(0usize, 1), (2, 1), (3, 4), (5, 4)] {
            let query = pieces(&dir, "refused.fa", Some(refused), b"");
            let args = args_of(&query);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            for (batch, per_batch) in [("400", 1), ("800", 2), ("100000", 6)] {
                let what = format!("outfmt {outfmt}, refused query {refused}, batch {batch}");
                // BLAST XML: the batches before the refused title's; XML2 and JSON: the reports
                // of the queries before it, each made at the next query's arrival.
                let written = match outfmt {
                    "5" => refused / per_batch * per_batch,
                    _ => refused,
                };
                for switches in [&[][..], &[STATS], &[OFF]] {
                    let mut run = blastn(&args, threads, batch, switches, None);
                    if switches == [STATS] {
                        let (batches, reason) = take_stats_lines(&mut run);
                        assert_eq!(reason, "none", "{what}");
                        assert!(batches >= (refused / per_batch) as u64, "{what}: {batches}");
                    }
                    let stderr = String::from_utf8_lossy(&run.stderr);
                    assert_eq!(run.status, Some(1), "{what} {switches:?}: {stderr}");
                    assert!(
                        stderr.contains(&format!("a query title {REFUSAL}")),
                        "{what}: {stderr}"
                    );
                    // The reports after the loop (`OFF`): BLAST XML refuses before its first byte.
                    let written = if switches == [OFF] && outfmt == "5" {
                        0
                    } else {
                        written
                    };
                    assert!(
                        whole.report.starts_with(&run.report)
                            && run.report.len() < whole.report.len()
                            && reports_in(outfmt, &run.report) == written,
                        "{what} {switches:?}: {} reports, not {written}, of\n{}",
                        reports_in(outfmt, &run.report),
                        String::from_utf8_lossy(&run.report)
                    );
                    if written == 0 {
                        assert!(run.report.is_empty(), "{what} {switches:?}: nothing");
                    } else if outfmt != "5" {
                        // Where the next report would start (outfmt 15's separator is that
                        // report's).
                        let rest = &whole.report[run.report.len()..];
                        let next: &[&[u8]] = if outfmt == "16" {
                            &[b"<BlastOutput2>", b"\n<BlastOutput2>"]
                        } else {
                            &[b",\n{\n"]
                        };
                        assert!(
                            next.iter().any(|next| rest.starts_with(next)),
                            "{what} {switches:?}: ends inside a report"
                        );
                    }
                }
            }
        }
    }
}

/// Six records of 400 letters of LC738868 with the given deflines (`-parse_deflines` reads the
/// ids), one per batch of 400 letters.
fn records_titled(dir: &TempDir, name: &str, deflines: [&[u8]; 6]) -> PathBuf {
    let genome = letters("LC738868.fasta");
    let mut text = Vec::new();
    for (i, defline) in deflines.iter().enumerate() {
        let start = 3_000 + i * 17_000;
        record(&mut text, defline, &genome[start..start + 400]);
    }
    let path = dir.path(name);
    std::fs::write(&path, text).expect("write the records");
    path
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2403
/// Under `-parse_deflines`, the results of consecutive queries of one Seq-id join one held
/// `-subject` report of 16 and 15, across batches too (`XmlShow::joins_held`), and the report shows
/// its first query's id, title and length: a title that LOSAT refuses on a joined query is never
/// shown, so never refused (the bytes of the run without it), and one on the first query of a join
/// is refused where NCBI makes the joined report, after the reports before it, with the switch
/// and without it, a batch of each query or one batch.
#[test]
fn joined_results_of_one_seq_id_are_refused_where_ncbi_makes_their_report() {
    let dir = TempDir::new("joined");
    let subject = fasta("LC738868.fasta");
    let plain = records_titled(
        &dir,
        "plain.fa",
        [
            b"lcl|qa first a",
            b"lcl|qa second a",
            b"lcl|qb first b",
            b"lcl|qb second b",
            b"lcl|qc first c",
            b"lcl|qc second c",
        ],
    );
    let joined_refused = records_titled(
        &dir,
        "joined.fa",
        [
            b"lcl|qa first a",
            b"lcl|qa second a",
            b"lcl|qb first b",
            b"lcl|qb , ,",
            b"lcl|qc first c",
            b"lcl|qc second c",
        ],
    );
    let lead_refused = records_titled(
        &dir,
        "lead.fa",
        [
            b"lcl|qa first a",
            b"lcl|qa second a",
            b"lcl|qb , ,",
            b"lcl|qb second b",
            b"lcl|qc first c",
            b"lcl|qc second c",
        ],
    );
    for (outfmt, threads) in [("16", 1), ("15", 4)] {
        for batch in ["400", "100000"] {
            let args_of = |query: &Path| -> Vec<String> {
                [
                    "-task",
                    "megablast",
                    "-parse_deflines",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-outfmt",
                    outfmt,
                ]
                .map(String::from)
                .to_vec()
            };
            let what = format!("outfmt {outfmt} batch {batch}");
            let args = args_of(&plain);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let (whole, _) = same_two(&what, &args, threads, batch, None);
            assert_eq!(whole.status, Some(0), "{what}: {whole:?}");
            assert_eq!(reports_in(outfmt, &whole.report), 3, "{what}: one each id");
            let args = args_of(&joined_refused);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let (joined, _) = same_two(&what, &args, threads, batch, None);
            assert!(joined == whole, "{what}: the joined title is not shown\n{joined:?}");
            let args = args_of(&lead_refused);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let (lead, _) = same_two(&what, &args, threads, batch, None);
            let stderr = String::from_utf8_lossy(&lead.stderr);
            assert_eq!(lead.status, Some(1), "{what}: {stderr}");
            assert!(stderr.contains(REFUSAL), "{what}: {stderr}");
            assert!(
                whole.report.starts_with(&lead.report) && reports_in(outfmt, &lead.report) == 1,
                "{what}: qa's report alone\n{}",
                String::from_utf8_lossy(&lead.report)
            );
        }
    }
}

#[test]
fn a_held_report_that_an_error_stops_first_is_dropped_with_its_refusal() {
    // The last query's `-subject` report of 16 and 15, held when NCBI's reader error stops the
    // next batch: NCBI never makes it, so the title it would crash on is never met, in both
    // paths: the reports of the batches before it, then the reader's error.
    let dir = TempDir::new("stopped");
    let subject = fasta("LC738868.fasta");
    let query = pieces(&dir, "q.fa", Some(5), b">bad\n@@@@\n");
    for (outfmt, threads) in [("16", 1), ("15", 4)] {
        let args = [
            "-task",
            "megablast",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
            "-outfmt",
            outfmt,
        ];
        let what = format!("outfmt {outfmt}");
        let (off, (batches, reason)) = same_two(&what, &args, threads, "400", None);
        let stderr = String::from_utf8_lossy(&off.stderr);
        assert_eq!(off.status, Some(1), "{what}: {stderr}");
        assert!(
            stderr.contains("CFastaReader") && !stderr.contains(REFUSAL),
            "{what}: {stderr}"
        );
        assert_eq!(reports_in(outfmt, &off.report), 5, "{what}");
        assert!(
            reason == "none" && batches >= 5,
            "{what}: {batches} {reason}"
        );
    }
}
