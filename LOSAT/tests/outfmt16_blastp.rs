#![allow(warnings, clippy::all)]

//! BLAST XML2 and JSON (outfmt 16 and 15) of BLASTP: what the frozen-fixture rows cannot show.
//! The `-subject` report that NCBI holds back one query (a failed batch loses it), the documents
//! of a run with no query, and the searches LOSAT refuses before the first byte (a byte >= 0x80 in
//! a title or the `-db` name, DECISIONS 19 Q2; a taxonomy database scientific name, DECISIONS 30)
//! and the taxonomy database of a `-db` taxid (DECISIONS 28). The bytes of NCBI BLAST+ 2.17.0 are in `tests/fixtures/outfmt16`
//! (`badbatch_{sub,db}_{15,16}.out`: `BATCH_SIZE=100 blastp -query qprot_bad.faa -subject
//! sprot.faa | -db db/sprot -num_threads 1 -outfmt N`, run in that directory).
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
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

/// Design §5.2 row 22 (`blastp_{sub,db}_badbatch_fmt{15,16}`): with `BATCH_SIZE=100` the third
/// query's batch fails to read. `-db` has written the reports of the two queries before it;
/// `-subject` only the first one's (the second is held back and lost), and neither has an epilog.
#[test]
fn a_failed_batch_keeps_the_db_report_and_loses_the_held_subject_report() {
    let dir = fixtures();
    for format in ["15", "16"] {
        for (mode, subject_args) in [
            ("sub", ["-subject", "sprot.faa"]),
            ("db", ["-db", "db/sprot"]),
        ] {
            let mut args = vec!["blastp", "-query", "qprot_bad.faa"];
            args.extend(subject_args);
            args.extend(["-num_threads", "1", "-outfmt", format]);
            let expected = fs::read(dir.join(format!("badbatch_{mode}_{format}.out"))).unwrap();
            for threads in ["1", "4"] {
                let mut args = args.clone();
                let at = args.iter().position(|a| *a == "-num_threads").unwrap();
                args[at + 1] = threads;
                let output = run_in(&dir, &args, &[("BATCH_SIZE", "100")]);
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{mode} {format}: {stderr}");
                assert!(
                    stderr.starts_with("BLAST query error: CFastaReader: Near line 10"),
                    "{mode} {format}: {stderr}"
                );
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&expected),
                    "{mode} {format} threads {threads}"
                );
            }
        }
    }
}

/// A run with no query formatted writes the prolog and the epilog (`PrintProlog`, `PrintEpilog`;
/// the probes `blastp_sub_emptypipe_fmt{15,16}`); an empty query file is `Query is Empty!` and no
/// bytes (the probe `blastp_sub_emptyfile_fmt{15,16}`).
#[test]
fn an_empty_query_pipe_writes_the_prolog_and_the_epilog() {
    let subject = fixtures().join("sprot.faa");
    let expected16 = "<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<BlastXML2\n    xmlns=\"http://www.ncbi.nlm.nih.gov\"\n    xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\n    xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\"\n></BlastXML2>\n";
    let expected15 = "{\n\"BlastOutput2\": [\n]\n}\n";
    let db = fixtures().join("db/sprot");
    for (format, expected) in [("16", expected16), ("15", expected15)] {
        for subject_args in [
            ["-subject", subject.to_str().unwrap()],
            ["-db", db.to_str().unwrap()],
        ] {
            for threads in ["1", "4"] {
                let mut args = vec!["blastp", "-query", "-"];
                args.extend(subject_args);
                args.extend(["-num_threads", threads, "-outfmt", format]);
                let output = run_with_stdin(&args, b"");
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(0), "{args:?}: {stderr}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout),
                    expected,
                    "{args:?}"
                );
            }
        }
        let dir = temp_dir("emptyfile");
        let empty = dir.join("empty.fa");
        fs::write(&empty, b"").unwrap();
        let output = run_in(
            &dir,
            &[
                "blastp",
                "-query",
                empty.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                format,
            ],
            &[],
        );
        assert_eq!(output.status.code(), Some(0));
        assert!(output.stdout.is_empty(), "Query is Empty!: no bytes");
        fs::remove_dir_all(dir).ok();
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

const PROTEIN: &[u8] = b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
const OTHER_PROTEIN: &[u8] = b"WQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNM";

fn refused(output: &Output, what: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{what}: {stderr}");
    assert!(output.stdout.is_empty(), "{what}: nothing is written");
    assert!(
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTP"
        )),
        "{stderr}"
    );
}

/// Q2: a byte >= 0x80 in a title or the `-db` name is refused for 15 and 16 as for 5, before the
/// first byte (NCBI cuts the 16 document at the byte, and writes the 15 one in UTF-8 after
/// guessing the encoding).
#[test]
fn a_title_or_database_name_outside_us_ascii_is_refused_in_15_and_16() {
    let dir = temp_dir("q2");
    let (query, subject) = (dir.join("q.fa"), dir.join("s.fa"));
    for format in ["15", "16"] {
        // a query title, first or second; UTF-8 and Latin-1
        for title in ["caf\u{e9} ".as_bytes(), b"caf\xe9"] {
            for position in 0..2 {
                let mut records: Vec<(&[u8], &[u8])> =
                    vec![(b"q1 ok", PROTEIN), (b"q2 ok", PROTEIN)];
                records[position] = (title, PROTEIN);
                fasta(&query, &records);
                fasta(&subject, &[(b"s1 subject", PROTEIN)]);
                let output = run_in(
                    &dir,
                    &[
                        "blastp",
                        "-query",
                        query.to_str().unwrap(),
                        "-subject",
                        subject.to_str().unwrap(),
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
                (b"s1 ok", OTHER_PROTEIN),
                ("s2 caf\u{e9}".as_bytes(), PROTEIN),
            ],
        );
        let output = run_in(
            &dir,
            &[
                "blastp",
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                format,
            ],
            &[],
        );
        refused(&output, "a subject title");
    }
    // the -db name
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let copy = dir.join("residues_caf\u{e9}");
    for extension in ["phr", "pin", "psq"] {
        fs::copy(
            tests.join(format!("residues.{extension}")),
            format!("{}.{extension}", copy.display()),
        )
        .expect("copy database");
    }
    for format in ["15", "16"] {
        let output = run_in(
            &dir,
            &[
                "blastp",
                "-query",
                tests.join("residues_q.faa").to_str().unwrap(),
                "-db",
                copy.to_str().unwrap(),
                "-outfmt",
                format,
            ],
            &[],
        );
        refused(&output, "a database name");
    }
    fs::remove_dir_all(dir).ok();
}

/// DECISIONS 28: a reported `-db` sequence with a taxid is written (`HitDescr_taxid`, the frozen
/// rows `db.sa.x16_taxid.{15,16}`) and, when `taxdb.bti` is on the search path (the current
/// directory, `BLASTDB`), NCBI asks the taxonomy database: whatever the file holds, a taxid it
/// cannot name leaves the `HitDescr` with the id alone (the empty `taxdb.bti` is a damaged
/// database). A database without taxids is not touched by the taxdb.
#[test]
fn a_taxdb_that_names_nothing_leaves_the_id_alone() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let taxid_db = tests.join("x16_taxid");
    let query = tests.join("x16_taxid_q.faa");
    let plain_db = tests.join("residues");
    let plain_query = tests.join("residues_q.faa");
    let dir = temp_dir("taxdb");
    fs::write(dir.join("taxdb.bti"), b"").unwrap();
    fn as_str(args: &[String]) -> Vec<&str> {
        args.iter().map(String::as_str).collect()
    }
    for format in ["15", "16"] {
        let args = |db: &Path, query: &Path| {
            vec![
                "blastp".to_string(),
                "-query".to_string(),
                query.to_str().unwrap().to_string(),
                "-db".to_string(),
                db.to_str().unwrap().to_string(),
                "-outfmt".to_string(),
                format.to_string(),
            ]
        };
        // without a taxdb: the taxid is written, no scientific name
        let work = temp_dir("nottaxdb");
        let output = run_in(&work, &as_str(&args(&taxid_db, &query)), &[]);
        assert_eq!(output.status.code(), Some(0));
        let text = String::from_utf8_lossy(&output.stdout);
        assert!(text.contains("9606"), "{text}");
        assert!(text.contains("accession"), "{text}");
        assert!(!text.contains("sciname"), "{text}");
        fs::remove_dir_all(work).ok();
        // taxdb.bti in BLASTDB, or in the current directory
        for (cwd, env) in [
            (
                temp_dir("nottaxdb2"),
                vec![("BLASTDB", dir.to_str().unwrap())],
            ),
            (dir.clone(), vec![]),
        ] {
            let output = run_in(&cwd, &as_str(&args(&taxid_db, &query)), &env);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(output.status.code(), Some(0), "{format} {env:?}: {stderr}");
            let text = String::from_utf8_lossy(&output.stdout);
            assert!(text.contains("gnl|BL_ORD_ID|"), "{text}");
            for gone in ["accession", "title", "9606", "sciname"] {
                assert!(
                    !text.contains(&format!("\"{gone}\"")) && !text.contains(&format!("<{gone}>")),
                    "{gone} in {text}"
                );
            }
            // the hits keep their lengths and alignments
            assert!(text.contains("hsps"), "{text}");
            assert!(stderr.is_empty(), "{stderr}");
            if cwd != dir {
                fs::remove_dir_all(cwd).ok();
            }
            // a database without taxids is not touched by the taxdb
            let output = run_in(&dir, &as_str(&args(&plain_db, &plain_query)), &[]);
            assert_eq!(
                output.status.code(),
                Some(0),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert!(!String::from_utf8_lossy(&output.stdout).contains("sciname"));
        }
    }
    fs::remove_dir_all(dir).ok();
}

/// The rows `db.sa.x16_taxid.*` of `db_search_manifest.tsv` (OF5a part 1b: outfmt 15 and 16 with
/// a taxonomy database in `BLASTDB`; the rows with `warn` in their name also have a query file with
/// a reader warning, and some a `BATCH_SIZE`), at one and four threads: the bytes of NCBI BLAST+
/// 2.17.0 (`blastp -num_threads 1`) for stdout, stderr and the exit status; the rows of kind
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
        .filter(|line| line.starts_with("db.sa.x16_taxid."))
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
        if env.is_empty() {
            continue;
        }
        for threads in ["1", "4"] {
            let output = run_in(
                engine,
                &[
                    "blastp",
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
    assert!(rows >= 20, "{rows} rows");
    assert!(warn_rows >= 4, "{warn_rows} rows with a reader warning");
}

/// What NCBI BLAST+ 2.17.0 wrote for `blastp -query QUERY (-subject sprot.faa | -db db/sprot)
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
                let mut args = vec!["blastp", "-query", query];
                args.extend(if mode == "sub" {
                    ["-subject", "sprot.faa"]
                } else {
                    ["-db", "db/sprot"]
                });
                args.extend(extra);
                args.extend(["-num_threads", threads, "-outfmt", format]);
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

/// A first batch that fails to read (the whole file is one batch, `firstbad_blastp_*`): NCBI has
/// written the prolog (`PrintProlog` comes before the first batch is read) and nothing else; the
/// error and the exit status 1 are NCBI's.
#[test]
fn a_first_batch_that_fails_to_read_leaves_the_prolog_alone() {
    equals_frozen("firstbad_blastp", "qprot_bad.faa", &[], &[], 1);
}

/// A score matrix that the engine does not support (`-matrix blosum62`, in lower case) stops the
/// search after the prolog (`matrix_blastp_*`): exit status 3 and `unsupported score matrix`.
#[test]
fn an_unsupported_matrix_leaves_the_prolog_alone() {
    equals_frozen(
        "matrix_blastp",
        "qprot_unsearched.faa",
        &["-matrix", "blosum62"],
        &[],
        3,
    );
}

/// A query of X only is a batch NCBI does not search when `BATCH_SIZE` puts it alone: its `Search`
/// has no `hits` member and ` \nNo hits found` after the message (`unsearched_blastp_*`:
/// `BATCH_SIZE=100`, queries `qp1`, `inv1`, `qp2`, `inv2`), as for BLASTN.
#[test]
fn a_query_of_an_unsearched_batch_has_no_hits_member_and_the_no_hits_suffix() {
    equals_frozen(
        "unsearched_blastp",
        "qprot_unsearched.faa",
        &[],
        &[("BATCH_SIZE", "100")],
        0,
    );
    let out = fs::read_to_string(fixtures().join("unsearched_blastp_db_16.out")).unwrap();
    assert!(out.contains("No hits found &#xa;No hits found"), "{out}");
}

/// A query of an unsearched batch lists no `query-masking`, and a query of a searched batch lists
/// its SEG region (`unsearchedmask_blastp_*`: `BATCH_SIZE=10 blastp -seg yes`, queries `qp1` with a
/// run of Q, `inv1` of X, a query of Q that SEG masks whole, `qp2`): the region of `qp1` is the
/// only one in the document.
#[test]
fn a_query_of_an_unsearched_batch_lists_no_masks() {
    equals_frozen(
        "unsearchedmask_blastp",
        "qprot_unsearched_seg.faa",
        &["-seg", "yes"],
        &[("BATCH_SIZE", "10")],
        0,
    );
    let out = fs::read_to_string(fixtures().join("unsearchedmask_blastp_sub_16.out")).unwrap();
    assert_eq!(out.matches("<query-masking>").count(), 1, "{out}");
}

/// The taxonomy database is read for outfmt 15 and 16 only (DECISIONS 30): a `taxdb.bti` of the
/// search path that cannot be read (here a damaged one) changes no byte of outfmt 0, 5, 6, 7 or 8
/// and adds no line to standard error, while 15 and 16 refuse it. The database has taxids, so the
/// lookup would run if any other format asked for it.
#[test]
fn the_taxdb_is_not_opened_outside_outfmt_15_and_16() {
    let engine = Path::new(env!("CARGO_MANIFEST_DIR"));
    let (db, query) = (
        "tests/db_search/x16_taxid",
        "tests/db_search/x16_taxid_q.faa",
    );
    let damaged = "tests/seqdb/taxdb/btd_empty";
    for format in ["0", "5", "6", "7", "8"] {
        for threads in ["1", "4"] {
            let args = [
                "blastp",
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
            &["blastp", "-db", db, "-query", query, "-outfmt", format],
            &[("BLASTDB", damaged)],
        );
        assert_eq!(output.status.code(), Some(1), "outfmt {format}");
        assert!(output.stdout.is_empty(), "outfmt {format}");
    }
}
