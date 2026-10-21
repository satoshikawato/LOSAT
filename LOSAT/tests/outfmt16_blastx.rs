#![allow(warnings, clippy::all)]

//! BLAST XML2 and JSON (outfmt 16 and 15) of BLASTX: what the frozen-fixture rows cannot show.
//! The `-subject` report that NCBI holds back one query (a failed batch loses it), the run that stops
//! at its first batch (the prolog alone), and the searches LOSAT refuses before the first byte (a byte
//! >= 0x80 in a title or the `-db` name, DECISIONS 19 Q2; a taxonomy database scientific name, DECISIONS
//! 30; a query that is not a regular file, DECISIONS 22), with the output held in memory when a later
//! batch can still refuse (DECISIONS 21), and the taxonomy database of a `-db` taxid (DECISIONS 28).
//! BLASTX reads its queries in batches of 10,002 bases and not `BATCH_SIZE`, so the queries of the
//! failed-batch tests are padded to 10,003 bases (`tests/db_search/make_x16_blastx_inputs.py`). The bytes
//! of NCBI BLAST+ 2.17.0 are in `tests/fixtures/outfmt16` (`badbatch_blastx_{sub,db,dbtax}_{15,16}`:
//! `blastx -query qnx_bad.fna -subject sprot.faa | -db db/sprot | -db db/sprot_tax (BLASTDB=the `good`
//! taxdb) -num_threads 1 -outfmt N`, run in that directory; `firstbad_blastx_*` with `qnx_firstbad.fna`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2363-2404, :2206-2219
//! (the `-subject` result is written when the next query's arrives, or by `PrintEpilog`, which a
//! failed batch does not reach).

use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/outfmt16")
}

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt16_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn command(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_X_DEFAULT" {
            command.env_remove(name);
        }
    }
    // NCBI's variables change a report (and the search path): none leaks in, and `~/.ncbirc`
    // (which can name a taxonomy database or BLASTDB) is not read: HOME is an empty directory.
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("NCBI") {
            command.env_remove(name);
        }
    }
    for name in ["BATCH_SIZE", "CHUNK_SIZE", "BLASTDB"] {
        command.env_remove(name);
    }
    let home = std::env::temp_dir().join("losat_outfmt16_home");
    fs::create_dir_all(&home).expect("empty home");
    command.env("HOME", &home).env("USERPROFILE", &home);
    command.env(
        "LOSAT_BLASTX_WORKER_LOG",
        std::env::temp_dir().join("losat_outfmt16_worker.log"),
    );
    command.envs(env.iter().copied());
    command
}

fn run_in(dir: &Path, args: &[&str], env: &[(&str, &str)]) -> Output {
    command(env)
        .current_dir(dir)
        .args(args)
        .output()
        .expect("run LOSAT")
}

fn run_with_stdin(args: &[&str], stdin: &[u8]) -> Output {
    use std::io::Write;
    let mut child = command(&[])
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run LOSAT");
    child.stdin.take().unwrap().write_all(stdin).unwrap();
    child.wait_with_output().expect("wait for LOSAT")
}

const PROTEIN: &str = "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
const OTHER_PROTEIN: &str = "WQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNM";

/// One codon per amino acid.
fn codon(letter: char) -> &'static str {
    match letter {
        'A' => "GCT",
        'C' => "TGC",
        'D' => "GAT",
        'E' => "GAA",
        'F' => "TTC",
        'G' => "GGT",
        'H' => "CAC",
        'I' => "ATC",
        'K' => "AAA",
        'L' => "CTG",
        'M' => "ATG",
        'N' => "AAC",
        'P' => "CCT",
        'Q' => "CAG",
        'R' => "CGT",
        'S' => "TCT",
        'T' => "ACC",
        'V' => "GTG",
        'W' => "TGG",
        'Y' => "TAC",
        other => panic!("no codon for {other}"),
    }
}

fn back_translate(protein: &str) -> String {
    protein.chars().map(codon).collect()
}

/// `core` followed by repeats of a short sequence up to 10,003 bases: NCBI's query batch of 10,002
/// bases is then closed by this query alone.
fn batch_of_its_own(core: &str) -> String {
    let mut query = core.to_string();
    while query.len() < 10_003 {
        query.push_str("ACGTTGCAAGCTTGCA");
    }
    query
}

fn fasta(path: &Path, records: &[(&[u8], &[u8])]) {
    let mut text = Vec::new();
    for (title, sequence) in records {
        text.extend_from_slice(b">");
        text.extend_from_slice(title);
        text.push(b'\n');
        text.extend_from_slice(sequence);
        text.push(b'\n');
    }
    fs::write(path, text).expect("write FASTA");
}

/// Standard error of the failed batches (`qnx_bad.fna`, `qnx_firstbad.fna`).
const READER_ERROR: &str = "BLAST query error: CFastaReader: Near line";

/// `blastx -query QUERY (-subject sprot.faa | -db db/sprot | -db db/sprot_tax with the `good` taxdb in
/// BLASTDB)` in `tests/fixtures/outfmt16`: what NCBI BLAST+ 2.17.0 wrote in `<stem>_{sub,db,dbtax}_{15,16}`
/// (`.out`, and `.err` when stderr was not empty), with the exit status `rc`. LOSAT, at one and four
/// threads, writes the same stdout and stderr, which starts as `stderr_start`.
fn equals_frozen(stem: &str, query: &str, modes: &[&str], rc: i32, stderr_start: &str) {
    let dir = fixtures();
    let taxdb = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/good");
    for format in ["15", "16"] {
        for &mode in modes {
            let name = format!("{stem}_{mode}_{format}");
            let out = fs::read(dir.join(format!("{name}.out"))).unwrap();
            let err = fs::read(dir.join(format!("{name}.err"))).unwrap_or_default();
            for threads in ["1", "4"] {
                let mut args = vec!["blastx", "-query", query];
                args.extend(match mode {
                    "sub" => ["-subject", "sprot.faa"],
                    "db" => ["-db", "db/sprot"],
                    _ => ["-db", "db/sprot_tax"],
                });
                args.extend(["-num_threads", threads, "-outfmt", format]);
                let env: Vec<(&str, &str)> = if mode == "dbtax" {
                    vec![("BLASTDB", taxdb.to_str().unwrap())]
                } else {
                    vec![]
                };
                let output = run_in(&dir, &args, &env);
                let case = format!("{name} threads {threads}");
                assert_eq!(output.status.code(), Some(rc), "{case}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&out),
                    "{case}: stdout"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr),
                    String::from_utf8_lossy(&err),
                    "{case}: stderr"
                );
                assert!(
                    String::from_utf8_lossy(&output.stderr).starts_with(stderr_start),
                    "{case}"
                );
            }
        }
    }
}

/// Design §5.2 row 22 for BLASTX (`blastx_{sub,db}_badbatch_fmt{15,16}`): the third query's batch fails to
/// read. `-db` has written the reports of the two queries before it; `-subject` only the first one's (the
/// second is held back and lost), and neither has an epilog. `dbtax` is a `-db` search with the taxonomy
/// database on the search path: the output is written batch by batch there too (the taxonomy database holds
/// no name LOSAT would refuse, DECISIONS 32), and the run that ends with the reader error has written the
/// batches before it.
#[test]
fn a_failed_batch_keeps_the_db_report_and_loses_the_held_subject_report() {
    equals_frozen(
        "badbatch_blastx",
        "qnx_bad.fna",
        &["sub", "db", "dbtax"],
        1,
        READER_ERROR,
    );
    let dbtax = fs::read_to_string(fixtures().join("badbatch_blastx_dbtax_16.out")).unwrap();
    assert!(dbtax.contains("<sciname>Homo sapiens</sciname>"), "{dbtax}");
}

/// A first batch that fails to read (the whole file is one batch): NCBI has written the prolog and nothing
/// else (`PrintProlog` comes before the first batch is read); the error and the exit status 1 are NCBI's.
#[test]
fn a_first_batch_that_fails_to_read_leaves_the_prolog_alone() {
    equals_frozen(
        "firstbad_blastx",
        "qnx_firstbad.fna",
        &["sub", "db"],
        1,
        READER_ERROR,
    );
}

/// A query that is a header and no sequence stops the search with `Sequence contains no data` (exit
/// status 3) after the prolog (`emptyseq_blastx_*`): the document is the prolog alone, as NCBI's.
#[test]
fn a_query_without_a_sequence_leaves_the_prolog_alone() {
    equals_frozen(
        "emptyseq_blastx",
        "qnx_empty.fna",
        &["sub", "db"],
        3,
        "BLAST engine error: Warning: Sequence contains no data",
    );
}

fn refused(output: &Output, what: &str, label: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{label}: {stderr}");
    assert!(output.stdout.is_empty(), "{label}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTX"
        )),
        "{label}: {stderr}"
    );
}

/// Q2: a byte >= 0x80 in a title or the `-db` name is refused for 15 and 16 as for 5, before the first byte,
/// also when the title is that of a later query batch or of a subject that only a later batch hits.
#[test]
fn a_title_or_database_name_outside_us_ascii_is_refused_in_15_and_16() {
    let dir = temp_dir("q2");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    let nucleotide = back_translate(PROTEIN);
    let padded = batch_of_its_own("");
    for format in ["15", "16"] {
        for threads in ["1", "4"] {
            let run = || {
                run_in(
                    &dir,
                    &[
                        "blastx",
                        "-query",
                        query.to_str().unwrap(),
                        "-subject",
                        subject.to_str().unwrap(),
                        "-num_threads",
                        threads,
                        "-outfmt",
                        format,
                    ],
                    &[],
                )
            };
            // a query title, first or second, in the same batch or in a later one; UTF-8
            for position in 0..2 {
                for later_batch in [false, true] {
                    let filler: &[u8] = if later_batch {
                        padded.as_bytes()
                    } else {
                        nucleotide.as_bytes()
                    };
                    let mut records: Vec<(&[u8], &[u8])> =
                        vec![(b"q1 ok", filler), (b"q2 ok", nucleotide.as_bytes())];
                    records[position].0 = "caf\u{e9} ".as_bytes();
                    fasta(&query, &records);
                    fasta(&subject, &[(b"s1 subject", PROTEIN.as_bytes())]);
                    refused(
                        &run(),
                        "a query title",
                        &format!("{format} {threads} {position} {later_batch}"),
                    );
                }
            }
            // a subject title with a hit, in a batch after one without hits
            fasta(
                &query,
                &[
                    (b"q1 no hit", padded.as_bytes()),
                    (b"q2 ok", nucleotide.as_bytes()),
                ],
            );
            fasta(
                &subject,
                &[
                    (b"s1 ok", OTHER_PROTEIN.as_bytes()),
                    ("s2 caf\u{e9}".as_bytes(), PROTEIN.as_bytes()),
                ],
            );
            refused(
                &run(),
                "a subject title",
                &format!("{format} {threads} subject"),
            );
        }
    }
    // the -db name
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let copy = dir.join("bslice_caf\u{e9}");
    for extension in ["phr", "pin", "psq"] {
        fs::copy(
            tests.join(format!("bslice.{extension}")),
            format!("{}.{extension}", copy.display()),
        )
        .expect("copy database");
    }
    for format in ["15", "16"] {
        for threads in ["1", "4"] {
            let output = run_in(
                &dir,
                &[
                    "blastx",
                    "-query",
                    tests.join("bslice_q.fna").to_str().unwrap(),
                    "-db",
                    copy.to_str().unwrap(),
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                ],
                &[],
            );
            refused(
                &output,
                "a database name",
                &format!("{format} {threads} db"),
            );
        }
    }
    fs::remove_dir_all(dir).ok();
}

/// DECISIONS 22 for 15 and 16: a query that is not a regular file cannot be read twice (the pre-scan of the
/// titles, DECISIONS 21), so it is refused before the search with nothing written; the other formats read
/// it once.
#[cfg(unix)]
#[test]
fn a_query_that_is_not_a_regular_file_is_refused_in_15_and_16() {
    use std::io::Write;
    let dir = temp_dir("pipe");
    let subject = dir.join("s.fa");
    fasta(&subject, &[(b"s1 subject", PROTEIN.as_bytes())]);
    let text = format!(">q1 plain\n{}\n", back_translate(PROTEIN)).into_bytes();
    for (format, expected) in [("15", true), ("16", true), ("6", false)] {
        let mut child = command(&[])
            .args([
                "blastx",
                "-query",
                "/dev/stdin",
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                format,
                "-evalue",
                "1e-10",
            ])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("run LOSAT");
        let mut stdin = child.stdin.take().expect("stdin");
        let text = text.clone();
        let writer = std::thread::spawn(move || {
            let _ = stdin.write_all(&text);
        });
        let output = child.wait_with_output().expect("wait");
        writer.join().expect("writer thread");
        let stderr = String::from_utf8_lossy(&output.stderr);
        if expected {
            assert_eq!(output.status.code(), Some(1), "{format}: {stderr}");
            assert!(output.stdout.is_empty(), "{format}: nothing is written");
            assert!(
                stderr.contains(&format!(
                    "unsupported BLASTX query input for -outfmt {format}: use a regular file"
                )),
                "{format}: {stderr}"
            );
        } else {
            assert_eq!(output.status.code(), Some(0), "{format}: {stderr}");
            assert!(String::from_utf8_lossy(&output.stdout).starts_with("q1\ts1\t"));
        }
    }
    fs::remove_dir_all(dir).ok();
}

/// DECISIONS 28 and 30: a `-db` search with the taxonomy database on the search path can still refuse in a
/// later batch (here the scientific name of a taxid has a byte >= 0x80, and the sequence that carries it is
/// first hit by the second batch). The output of the first batch is held in memory and dropped with the
/// refusal: nothing is written. Without the taxonomy database the same search writes both reports.
#[test]
fn a_taxdb_refusal_in_a_later_batch_leaves_nothing_written() {
    let dir = fixtures();
    let taxdb = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/high_bytes_utf8");
    for format in ["15", "16"] {
        for threads in ["1", "4"] {
            let args = [
                "blastx",
                "-query",
                "qnx_late.fna",
                "-db",
                "db/sprot_tax",
                "-num_threads",
                threads,
                "-outfmt",
                format,
            ];
            let without = run_in(&dir, &args, &[]);
            assert_eq!(without.status.code(), Some(0), "{format} {threads}");
            let text = String::from_utf8_lossy(&without.stdout);
            // two sequences are hit, both by the second query
            assert_eq!(text.matches("gnl|BL_ORD_ID|").count(), 2, "{text}");
            assert!(text.contains("9606") && !text.contains("sciname"), "{text}");
            let with = run_in(&dir, &args, &[("BLASTDB", taxdb.to_str().unwrap())]);
            let stderr = String::from_utf8_lossy(&with.stderr);
            assert_eq!(with.status.code(), Some(1), "{format} {threads}: {stderr}");
            assert!(
                with.stdout.is_empty(),
                "{format} {threads}: nothing is written"
            );
            assert!(
                stderr.contains("a taxonomy database scientific name with a character outside the US-ASCII range"),
                "{stderr}"
            );
        }
    }
}

/// DECISIONS 32: the output of a `-db` search with the taxonomy database on the search path is held in memory
/// only when some sequence of the database has a taxid whose lookup would be refused (the pre-pass,
/// `algorithm::common::database::taxdb_may_refuse`, unit-tested in `taxdb_prepass_tests`). `db/sprot_mixtax` is
/// `sprot.faa` with taxid 10090 and a third protein, never hit, with taxid 9606; `qnx_late.fna` is two batches
/// (the second hits the first two sequences). With the `good` taxonomy database nothing is held; with
/// `high_bytes_utf8` (the name of 9606 has bytes >= 0x80) the pre-pass finds the unreported taxid and holds, and
/// the run writes what NCBI BLAST+ 2.17.0 writes (`mixtax_blastx_{good,badname}_{15,16}`: it never looks 9606 up).
/// The same taxonomy database with `db/sprot_tax`, where 9606 is reported in the second batch, is refused with
/// nothing written (`a_taxdb_refusal_in_a_later_batch_leaves_nothing_written`).
#[test]
fn a_taxid_that_is_in_the_database_but_never_reported_does_not_refuse_the_run() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for (case, taxdb) in [("good", "good"), ("badname", "high_bytes_utf8")] {
            let name = format!("mixtax_blastx_{case}_{format}");
            let out = fs::read(dir.join(format!("{name}.out"))).unwrap();
            let err = fs::read(dir.join(format!("{name}.err"))).unwrap_or_default();
            let taxdb = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/seqdb/taxdb")
                .join(taxdb);
            for threads in ["1", "4"] {
                let args = [
                    "blastx",
                    "-query",
                    "qnx_late.fna",
                    "-db",
                    "db/sprot_mixtax",
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                ];
                let output = run_in(&dir, &args, &[("BLASTDB", taxdb.to_str().unwrap())]);
                let case = format!("{name} threads {threads}");
                assert_eq!(output.status.code(), Some(0), "{case}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&out),
                    "{case}: stdout"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr),
                    String::from_utf8_lossy(&err),
                    "{case}: stderr"
                );
            }
        }
    }
}

/// The rows `db.x16_bx_taxid.*` of `db_search_manifest.tsv` (outfmt 15 and 16 of BLASTX with a taxonomy
/// database in `BLASTDB`, the offset error of `[blastx]` on standard error), at one and four threads: the
/// bytes of NCBI BLAST+ 2.17.0 (`blastx -num_threads 1`) for stdout, stderr and the exit status; the rows
/// of kind `losat-rejection` write nothing and name their reason.
#[test]
fn the_taxdb_rows_match_the_frozen_ncbi_bytes_at_one_and_four_threads() {
    let engine = Path::new(env!("CARGO_MANIFEST_DIR"));
    let manifest = fs::read_to_string(engine.join("tests/db_search_manifest.tsv")).unwrap();
    let mut rows = 0;
    for line in manifest
        .lines()
        .filter(|line| line.starts_with("db.x16_bx_taxid."))
    {
        let cells: Vec<&str> = line.split('\t').collect();
        let (id, kind, db, query, extra, outfmt, rc, covers) = (
            cells[0], cells[2], cells[3], cells[4], cells[5], cells[6], cells[7], cells[11],
        );
        let mut env: Vec<(&str, &str)> = Vec::new();
        for token in extra.split_whitespace() {
            let (name, value) = token
                .split_once('=')
                .unwrap_or_else(|| panic!("{id}: the extra argument {token:?} is not NAME=VALUE"));
            assert!(
                ["BLASTDB", "BATCH_SIZE", "CHUNK_SIZE"].contains(&name),
                "{id}: {token}"
            );
            env.push((name, value));
        }
        for threads in ["1", "4"] {
            let output = run_in(
                engine,
                &[
                    "blastx",
                    "-db",
                    db,
                    "-query",
                    query,
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    threads,
                ],
                &env,
            );
            let rc: i32 = rc.parse().unwrap();
            assert_eq!(output.status.code(), Some(rc), "{id} threads {threads}");
            if kind == "losat-rejection" {
                let expect = covers.split("expect: ").nth(1).unwrap();
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert!(output.stdout.is_empty(), "{id}: nothing is written");
                assert!(stderr.contains(expect), "{id}: {stderr}");
                continue;
            }
            let frozen = engine.join("tests/fixtures/db_search");
            let out = fs::read(frozen.join(format!("{id}.out"))).unwrap();
            let err = fs::read(frozen.join(format!("{id}.err"))).unwrap_or_default();
            assert!(
                output.stdout == out,
                "{id} threads {threads}: stdout differs"
            );
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                String::from_utf8_lossy(&err),
                "{id} threads {threads}: stderr"
            );
        }
        rows += 1;
    }
    assert!(rows >= 8, "{rows} rows");
}

/// The taxonomy database is read for outfmt 15 and 16 only (DECISIONS 30): a `taxdb.bti` of the search path
/// that cannot be read (here a damaged one) changes no byte of outfmt 0, 5, 6, 7 or 8 and adds no line to
/// standard error, while 15 and 16 refuse it. The database has taxids, so the lookup would run if any other
/// format asked for it.
#[test]
fn the_taxdb_is_not_opened_outside_outfmt_15_and_16() {
    let engine = Path::new(env!("CARGO_MANIFEST_DIR"));
    let (db, query) = (
        "tests/db_search/x16_taxid",
        "tests/db_search/x16_bx_taxid_q.fna",
    );
    let damaged = "tests/seqdb/taxdb/btd_empty";
    for format in ["0", "5", "6", "7", "8"] {
        for threads in ["1", "4"] {
            let args = [
                "blastx",
                "-db",
                db,
                "-query",
                query,
                "-outfmt",
                format,
                "-num_threads",
                threads,
            ];
            let without = run_in(engine, &args, &[]);
            let with = run_in(engine, &args, &[("BLASTDB", damaged)]);
            let case = format!("outfmt {format} threads {threads}");
            assert_eq!(without.status.code(), Some(0), "{case}");
            assert_eq!(with.status.code(), Some(0), "{case}");
            assert!(with.stderr.is_empty(), "{case}: {:?}", with.stderr);
            assert!(!with.stdout.is_empty(), "{case}");
            assert!(with.stdout == without.stdout, "{case}: stdout");
        }
    }
    for format in ["15", "16"] {
        let output = run_in(
            engine,
            &["blastx", "-db", db, "-query", query, "-outfmt", format],
            &[("BLASTDB", damaged)],
        );
        assert_eq!(output.status.code(), Some(1), "outfmt {format}");
        assert!(output.stdout.is_empty(), "outfmt {format}");
    }
}
