#![allow(warnings, clippy::all)]

//! BLAST XML2 and JSON (outfmt 16 and 15) of TBLASTN: what the frozen-fixture rows cannot show.
//! The `-subject` report that NCBI holds back one query (a failed batch loses it), the documents
//! of a run with no query, the `-matrix blosum62` error that NCBI raises after the prolog, and
//! the searches LOSAT refuses before the first byte (a byte >= 0x80 in a title or the `-db`
//! name, DECISIONS 19 Q2; a taxonomy database scientific name, DECISIONS 30: the rows
//! `db.sa.x16_ttaxid.nonascii.*` of `db_search_manifest.tsv`), the taxonomy database that only
//! outfmt 15 and 16 read, and the sink of a file. The bytes of NCBI BLAST+ 2.17.0 are in
//! `tests/fixtures/outfmt16` (`*_tblastn_{sub,db}_{15,16}.out`: `tblastn -query QUERY -subject
//! snucl.fna | -db db/snucl -num_threads 1 -outfmt N`, run in that directory; `badbatch` and
//! `unsearched` with `BATCH_SIZE=100`).
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

const SUBJECT: &str = "snucl.fna";
const DB: &str = "db/snucl";

/// The command of a frozen case: `tblastn -query QUERY (-subject snucl.fna | -db db/snucl)`, run
/// in `tests/fixtures/outfmt16` (the db name is part of the bytes).
fn args<'a>(query: &'a str, mode: &str, threads: &'a str, format: &'a str) -> Vec<&'a str> {
    let mut args = vec!["tblastn", "-query", query];
    args.extend(if mode == "sub" {
        ["-subject", SUBJECT]
    } else {
        ["-db", DB]
    });
    args.extend(["-num_threads", threads, "-outfmt", format]);
    args
}

fn frozen(name: &str) -> String {
    String::from_utf8(fs::read(fixtures().join(name)).unwrap()).unwrap()
}

/// What NCBI BLAST+ 2.17.0 wrote for `tblastn -query QUERY (-subject snucl.fna | -db db/snucl)
/// EXTRA -outfmt N -num_threads 1` in `tests/fixtures/outfmt16` (`<stem>_{sub,db}_{15,16}.out`, and
/// `.err` when stderr was not empty), with `env` and the exit status `rc`: LOSAT at one and four
/// threads writes the same stdout, stderr and exit status.
fn equals_frozen(stem: &str, query: &str, extra: &[&str], env: &[(&str, &str)], rc: i32) {
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let name = format!("{stem}_{mode}_{format}");
            let out = fs::read(dir.join(format!("{name}.out"))).unwrap();
            let err = fs::read(dir.join(format!("{name}.err"))).unwrap_or_default();
            for threads in ["1", "4"] {
                let mut args = args(query, mode, threads, format);
                args.extend(extra);
                let output = run_in(&dir, &args, env);
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
            }
        }
    }
}

/// Design §5.2 row 22 for TBLASTN (`badbatch_tblastn_*`): with `BATCH_SIZE=100` every query is a
/// batch and the third one fails to read. `-db` has written the reports of the two queries
/// before it; `-subject` the first one's only (the second is held back and lost), and neither
/// has an epilog. The document of a failed run is the same in a `-out` file.
#[test]
fn a_failed_batch_keeps_the_db_report_and_loses_the_held_subject_report() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("badbatch_tblastn_{mode}_{format}.out"));
            for threads in ["1", "4"] {
                let output = run_in(
                    &dir,
                    &args("qprot_bad.faa", mode, threads, format),
                    &[("BATCH_SIZE", "100")],
                );
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{mode} {format}: {stderr}");
                assert!(
                    stderr.starts_with("BLAST query error: CFastaReader: Near line 10"),
                    "{mode} {format}: {stderr}"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{mode} {format} threads {threads}"
                );
            }
            // a file sink: the same document (its start is written with the first report, or
            // by the error path when no report was written)
            let scratch = temp_dir("failed_out");
            let out = scratch.join("report.out");
            let mut args = args("qprot_bad.faa", mode, "1", format);
            args.extend(["-out", out.to_str().unwrap()]);
            let output = run_in(&dir, &args, &[("BATCH_SIZE", "100")]);
            assert_eq!(output.status.code(), Some(1), "{mode} {format} -out");
            assert!(output.stdout.is_empty());
            assert_eq!(
                String::from_utf8_lossy(&fs::read(&out).unwrap()),
                expected,
                "{mode} {format} -out"
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
}

/// A first batch that fails to read (the whole file is one batch, `firstbad_tblastn_*`): NCBI has
/// written the prolog (`PrintProlog` comes before the first batch is read) and nothing else, in
/// both modes and with or without `-out`.
#[test]
fn a_first_batch_that_fails_to_read_leaves_the_prolog_alone() {
    let dir = fixtures();
    equals_frozen("firstbad_tblastn", "qprot_bad.faa", &[], &[], 1);
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("firstbad_tblastn_{mode}_{format}.out"));
            let scratch = temp_dir("firstbad_out");
            let out = scratch.join("report.out");
            let mut args = args("qprot_bad.faa", mode, "1", format);
            args.extend(["-out", out.to_str().unwrap()]);
            let output = run_in(&dir, &args, &[]);
            assert_eq!(output.status.code(), Some(1), "{mode} {format} -out");
            assert_eq!(
                String::from_utf8_lossy(&fs::read(&out).unwrap()),
                expected,
                "{mode} {format} -out"
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
}

/// `-matrix blosum62` in lower case (design 4.5, probe `tblastn_sub_lcmatrix_*`): TBLASTN runs
/// the search, and NCBI's report data constructor throws `unsupported score matrix` when it makes
/// the first report, outside the stream's catch. Standard error has `BLAST engine error:
/// unsupported score matrix`, the exit status is 3, and the prolog is all that was written (16:
/// the header; 15: `{\n"BlastOutput2": [\n`). Without a query no report data is made, and the
/// document is the prolog and the epilog (`matrixempty_tblastn_*`).
#[test]
fn an_unsupported_matrix_leaves_the_prolog_alone() {
    equals_frozen(
        "matrix_tblastn",
        "qprot_unsearched.faa",
        &["-matrix", "blosum62"],
        &[],
        3,
    );
    let out = frozen("matrix_tblastn_sub_15.out");
    assert_eq!(out, "{\n\"BlastOutput2\": [\n");
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("matrixempty_tblastn_{mode}_{format}.out"));
            let (subject, db) = (
                dir.join(SUBJECT).to_str().unwrap().to_string(),
                dir.join(DB).to_str().unwrap().to_string(),
            );
            let mut args = vec!["tblastn", "-query", "-"];
            args.extend(if mode == "sub" {
                ["-subject", subject.as_str()]
            } else {
                ["-db", db.as_str()]
            });
            args.extend(["-matrix", "blosum62", "-outfmt", format]);
            let output = run_with_stdin(&args, b"");
            assert_eq!(output.status.code(), Some(0), "{mode} {format}");
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                expected,
                "{mode} {format}"
            );
        }
    }
}

/// A query of X only is a batch NCBI does not search when `BATCH_SIZE` puts it alone: its `Search`
/// has no `hits` member and ` \nNo hits found` after the message (`unsearched_tblastn_*`:
/// `BATCH_SIZE=100`, queries `qp1`, `inv1`, `qp2`), as for BLASTP.
#[test]
fn a_query_of_an_unsearched_batch_has_no_hits_member_and_the_no_hits_suffix() {
    equals_frozen(
        "unsearched_tblastn",
        "qprot_unsearched.faa",
        &[],
        &[("BATCH_SIZE", "100")],
        0,
    );
    let out = fs::read_to_string(fixtures().join("unsearched_tblastn_db_16.out")).unwrap();
    assert!(out.contains("No hits found &#xa;No hits found"), "{out}");
}

/// A query of an unsearched batch has no `query-masking` either, whatever its lower case
/// (`unsearchedmask_tblastn_*`: `BATCH_SIZE=10 -lcase_masking`, queries `qp1`, `inv_lc` (X with a
/// lower-case run) and `qp1b` (a hit with a lower-case run, which is listed)).
#[test]
fn a_query_of_an_unsearched_batch_lists_no_masks() {
    equals_frozen(
        "unsearchedmask_tblastn",
        "qprot_unsearched_lc.faa",
        &["-lcase_masking"],
        &[("BATCH_SIZE", "10")],
        0,
    );
    let out = fs::read_to_string(fixtures().join("unsearchedmask_tblastn_sub_16.out")).unwrap();
    assert_eq!(out.matches("<query-masking>").count(), 1, "{out}");
    assert!(out.contains("No hits found &#xa;No hits found"), "{out}");
}

/// A run with no query formatted writes the prolog and the epilog (`PrintProlog`, `PrintEpilog`;
/// `emptypipe_tblastn_*`); an empty query file is `Query is Empty!` and no bytes
/// (`emptyfile_tblastn_*`).
#[test]
fn an_empty_query_pipe_writes_the_prolog_and_the_epilog() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for mode in ["sub", "db"] {
            let expected = frozen(&format!("emptypipe_tblastn_{mode}_{format}.out"));
            for threads in ["1", "4"] {
                let (subject, db) = (
                    dir.join(SUBJECT).to_str().unwrap().to_string(),
                    dir.join(DB).to_str().unwrap().to_string(),
                );
                let mut args = vec!["tblastn", "-query", "-"];
                args.extend(if mode == "sub" {
                    ["-subject", subject.as_str()]
                } else {
                    ["-db", db.as_str()]
                });
                args.extend(["-num_threads", threads, "-outfmt", format]);
                let output = run_with_stdin(&args, b"");
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(0), "{mode} {format}: {stderr}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{mode} {format}"
                );
            }
            // an empty query file
            let scratch = temp_dir("emptyfile");
            let empty = scratch.join("empty.faa");
            fs::write(&empty, b"").unwrap();
            let output = run_in(&dir, &args(empty.to_str().unwrap(), mode, "1", format), &[]);
            assert_eq!(output.status.code(), Some(0));
            assert!(output.stdout.is_empty(), "Query is Empty!: no bytes");
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                frozen(&format!("emptyfile_tblastn_{mode}_{format}.err"))
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
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

/// A protein and a nucleotide sequence that it matches (`sn1` of `snucl.fna` translates to the
/// first part of `qp1`), and a nucleotide sequence that it does not match.
const PROTEIN: &[u8] =
    b"MWMEHSVLEAIFAVLRKVMDDPKIRDGPMAAMKGWENYCWPANISWKVLMKFEFSSDIHPCNLCICCLRIWQTQDCRRRIRCIGQQ";
/// The first record of `snucl.fna` (`sn1`).
fn nucleotide() -> Vec<u8> {
    let text = fs::read_to_string(fixtures().join("snucl.fna")).unwrap();
    text.lines()
        .skip(1)
        .take_while(|line| !line.starts_with('>'))
        .flat_map(|line| line.bytes())
        .collect()
}
const OTHER_NUCLEOTIDE: &[u8] = b"TCATGCTATCTATGGCACAAATGGACATTATTTTTGGAAGGTCGTTAGTGTCCGAGGGTCTCATAGATCCAGTTTAACGCCTTGGATTGATCATGCTATCTATGGCACAAATGGACATTATTTT";

fn refused(output: &Output, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    assert!(output.stdout.is_empty(), "{what}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's TBLASTN"
        )),
        "{stderr}"
    );
}

/// Q2: a byte >= 0x80 in a title or the `-db` name is refused for 15 and 16 as for 5, before the
/// first byte (NCBI cuts the 16 document at the byte, and writes the 15 one in UTF-8 after
/// guessing the encoding), at one and four threads.
#[test]
fn a_title_or_database_name_outside_us_ascii_is_refused_in_15_and_16() {
    let dir = temp_dir("q2");
    let sn1 = nucleotide();
    let (query, subject) = (dir.join("q.faa"), dir.join("s.fna"));
    for format in ["15", "16"] {
        for threads in ["1", "4"] {
            for title in ["caf\u{e9} ".as_bytes(), b"caf\xe9"] {
                for position in 0..2 {
                    let mut records: Vec<(&[u8], &[u8])> =
                        vec![(b"q1 ok", PROTEIN), (b"q2 ok", PROTEIN)];
                    records[position] = (title, PROTEIN);
                    fasta(&query, &records);
                    fasta(&subject, &[(b"s1 subject", &sn1[..])]);
                    let output = run_in(
                        &dir,
                        &[
                            "tblastn",
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
                    );
                    refused(&output, "a query title");
                }
            }
            // a subject title with a hit
            fasta(&query, &[(b"q1 ok", PROTEIN)]);
            fasta(
                &subject,
                &[
                    (b"s1 ok", OTHER_NUCLEOTIDE),
                    ("s2 caf\u{e9}".as_bytes(), &sn1[..]),
                ],
            );
            let output = run_in(
                &dir,
                &[
                    "tblastn",
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
            );
            refused(&output, "a subject title");
        }
    }
    // the -db name
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let copy = dir.join("x16_ntaxid_caf\u{e9}");
    for extension in ["nhr", "nin", "nsq"] {
        fs::copy(
            tests.join(format!("x16_ntaxid.{extension}")),
            format!("{}.{extension}", copy.display()),
        )
        .expect("copy database");
    }
    for format in ["15", "16"] {
        for threads in ["1", "4"] {
            let output = run_in(
                &dir,
                &[
                    "tblastn",
                    "-query",
                    tests.join("x16_ttaxid_q.faa").to_str().unwrap(),
                    "-db",
                    copy.to_str().unwrap(),
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                ],
                &[],
            );
            refused(&output, "a database name");
        }
    }
    fs::remove_dir_all(dir).ok();
}

/// A `-out` file holds the bytes of standard output for 15 and 16, in a search that succeeds
/// (the prolog comes with the first report, `Xml2Document`) at every thread count; the frozen
/// row `sa.x16.tblastn.sloc` is the same search on standard output.
#[test]
fn a_file_sink_holds_the_document_of_standard_output() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let query = tests.join("fasta/outfmt0/e2e_protein_query.faa");
    let subject = tests.join("fasta/outfmt0/e2e_tblastn_subject.fna");
    for format in ["15", "16"] {
        let expected = fs::read(
            tests
                .join("fixtures/outfmt0")
                .join(format!("sa.x16.tblastn.sloc.{format}.out")),
        )
        .unwrap();
        for threads in ["1", "4"] {
            let scratch = temp_dir("file_sink");
            let out = scratch.join("report.out");
            let output = run_in(
                &scratch,
                &[
                    "tblastn",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-subject_loc",
                    "2-470",
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                    "-out",
                    out.to_str().unwrap(),
                ],
                &[],
            );
            assert_eq!(output.status.code(), Some(0), "{format}");
            assert!(output.stdout.is_empty(), "{format}");
            assert_eq!(
                String::from_utf8_lossy(&fs::read(&out).unwrap()),
                String::from_utf8_lossy(&expected),
                "{format} threads {threads}"
            );
            fs::remove_dir_all(scratch).ok();
        }
    }
}

/// The rows `db.sa.x16_ttaxid.*` of `db_search_manifest.tsv` (OF5b part 2: outfmt 15 and 16 with
/// a taxonomy database in `BLASTDB`; the rows with `warn` in their name also have a query file with
/// a reader warning, and some a `BATCH_SIZE`), at one and four threads: the bytes of NCBI BLAST+
/// 2.17.0 (`tblastn -num_threads 1`) for stdout, stderr and the exit status; the rows of kind
/// `losat-rejection` write nothing and name their reason. Every row with a `NAME=VALUE` token
/// (`BLASTDB`, `BATCH_SIZE`, `CHUNK_SIZE`) runs with it in the environment, as the harness does.
#[test]
fn the_taxdb_rows_match_the_frozen_ncbi_bytes_at_one_and_four_threads() {
    let engine = Path::new(env!("CARGO_MANIFEST_DIR"));
    let manifest = fs::read_to_string(engine.join("tests/db_search_manifest.tsv")).unwrap();
    let mut rows = 0;
    let mut warn_rows = 0;
    for line in manifest
        .lines()
        .filter(|line| line.starts_with("db.sa.x16_ttaxid."))
    {
        let cells: Vec<&str> = line.split('\t').collect();
        // A row of outfmt 13 or 14 (the `files` column) writes files next to `-out`:
        // `tests/outfmt13_14_files.rs` and the harness check it.
        if cells.get(12) == Some(&"xinclude") {
            continue;
        }
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
                    "tblastn",
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
        warn_rows += id.contains(".warn.") as usize;
    }
    assert!(rows >= 10, "{rows} rows");
    assert!(warn_rows >= 3, "{warn_rows} rows with a reader warning");
}

/// The taxonomy database is read for outfmt 15 and 16 only (DECISIONS 30): a `taxdb.bti` of the
/// search path that cannot be read (here a damaged one) changes no byte of outfmt 0, 5, 6, 7 or 8
/// and adds no line to standard error, while 15 and 16 refuse it. The database has taxids, so the
/// lookup would run if any other format asked for it.
#[test]
fn the_taxdb_is_not_opened_outside_outfmt_15_and_16() {
    let engine = Path::new(env!("CARGO_MANIFEST_DIR"));
    let (db, query) = (
        "tests/db_search/x16_ntaxid",
        "tests/db_search/x16_ttaxid_q.faa",
    );
    let damaged = "tests/seqdb/taxdb/btd_empty";
    for format in ["0", "5", "6", "7", "8"] {
        for threads in ["1", "4"] {
            let args = [
                "tblastn",
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
            &["tblastn", "-db", db, "-query", query, "-outfmt", format],
            &[("BLASTDB", damaged)],
        );
        assert_eq!(output.status.code(), Some(1), "outfmt {format}");
        assert!(output.stdout.is_empty(), "outfmt {format}");
    }
}
