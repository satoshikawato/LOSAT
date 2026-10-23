#![allow(warnings, clippy::all)]

//! BLAST JSON and XML2 as files (outfmt 13 and 14): what the frozen-fixture rows of the manifests
//! (`files` column) cannot show. BLASTP and BLASTN first (OF5c part 1), then BLASTX, TBLASTN and TBLASTX
//! (part 2: the last section of this file). The names of the side files and of the entries of
//! the master file for the `-out` spellings of NCBI BLAST+ 2.17.0 (DESIGN_OF5 1.7, probes `naming/`
//! of `outfmt-20261019/S/probe5`), a side file that already exists (`stale14`) or cannot be opened
//! (`unwritable13`, `unwritable14`), a missing directory, the master of a run with no query and of
//! an empty query file, the `-subject` report that NCBI holds back one query (a failed batch loses
//! it, DESIGN_OF5 5.2 row 22), and the callers that have no output file.
//!
//! The bytes of NCBI are in `tests/fixtures/outfmt13_14`: the side files of the search of
//! `qprot.faa` against `outfmt16/sprot.faa` (`naming_<n>.xml|json`; the same for every spelling,
//! the file name does not enter a report) and of `outfmt16/qnucl_valid.fna` against
//! `outfmt16/snucl.fna` (`blastn_<n>`), the master of each spelling (`master_<case>`: they list the
//! names of the side files), and the files of the failed batch (`badbatch_{sub,db}_{13,14}.out*`:
//! `BATCH_SIZE=100 blastp -query qprot_bad.faa -subject sprot.faa | -db db/sprot -num_threads 1 -outfmt N
//! -out badbatch_{sub,db}_N.out`, run in `tests/fixtures/outfmt16`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2316-2353, :2425-2457.

mod run_local_support;

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn engine() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn inputs() -> PathBuf {
    engine().join("tests/fixtures/outfmt16")
}

fn frozen() -> PathBuf {
    engine().join("tests/fixtures/outfmt13_14")
}

fn read_frozen(name: &str) -> Vec<u8> {
    fs::read(frozen().join(name)).unwrap_or_else(|error| panic!("{name}: {error}"))
}

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt1314_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn command(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") || name_text.starts_with("NCBI") {
            command.env_remove(name);
        }
    }
    for name in ["BATCH_SIZE", "CHUNK_SIZE", "BLASTDB"] {
        command.env_remove(name);
    }
    let home = std::env::temp_dir().join("losat_outfmt1314_home");
    fs::create_dir_all(&home).expect("empty home");
    command.env("HOME", &home).env("USERPROFILE", &home);
    command.env(
        "LOSAT_BLASTX_WORKER_LOG",
        std::env::temp_dir().join("losat_outfmt1314_worker.log"),
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

fn run_with_stdin(dir: &Path, args: &[&str], stdin: &[u8]) -> Output {
    use std::io::Write;
    let mut child = command(&[])
        .current_dir(dir)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run LOSAT");
    child.stdin.take().unwrap().write_all(stdin).unwrap();
    child.wait_with_output().expect("wait for LOSAT")
}

/// The names of the entries of a directory, sorted.
fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

fn assert_ok(output: &Output, what: &str) {
    assert_eq!(
        output.status.code(),
        Some(0),
        "{what}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stdout.is_empty(), "{what}: standard output");
    assert!(output.stderr.is_empty(), "{what}: {:?}", output.stderr);
}

/// `S/probe5/out/naming` (DESIGN_OF5 1.7): the `-out` spelling, the format, the files that appear
/// (the master first) and the frozen master of NCBI, which lists `X_1.xml` for `d/X.XML`,
/// `Y.Json_1.json` for `d/Y.Json` (the XML rule for the master's names: DECISIONS 8), `k_1.json` for
/// `d/k.xml`, `m.a.json_1.json` for `d/m.a.json`, `j.json_1.xml`, `z_1.xml`, `_1.xml` and `rel_1.xml`.
const NAMING: [(&str, &str, [&str; 3], &str); 8] = [
    ("upper14", "14", ["X.XML", "X_1.xml", "X_2.xml"], "d/X.XML"),
    (
        "upper13",
        "13",
        ["Y.Json", "Y_1.json", "Y_2.json"],
        "d/Y.Json",
    ),
    ("noext14", "14", ["z", "z_1.xml", "z_2.xml"], "d/z"),
    (
        "json14",
        "14",
        ["j.json", "j.json_1.xml", "j.json_2.xml"],
        "d/j.json",
    ),
    (
        "xml13",
        "13",
        ["k.xml", "k.xml_1.json", "k.xml_2.json"],
        "d/k.xml",
    ),
    (
        "twodot13",
        "13",
        ["m.a.json", "m.a_1.json", "m.a_2.json"],
        "d/m.a.json",
    ),
    ("dotonly14", "14", [".xml", "_1.xml", "_2.xml"], "d/.xml"),
    (
        "rel14",
        "14",
        ["rel.xml", "rel_1.xml", "rel_2.xml"],
        "rel.xml",
    ),
];

fn blastp_args<'a>(
    query: &'a str,
    subject: &'a str,
    format: &'a str,
    out: &'a str,
) -> Vec<&'a str> {
    vec![
        "blastp",
        "-query",
        query,
        "-subject",
        subject,
        "-num_threads",
        "1",
        "-outfmt",
        format,
        "-out",
        out,
    ]
}

#[test]
fn blastp_names_the_side_files_and_the_master_entries_as_ncbi() {
    let query = frozen().join("qprot.faa");
    let subject = inputs().join("sprot.faa");
    let (query, subject) = (query.to_str().unwrap(), subject.to_str().unwrap());
    for (case, format, expected, out) in NAMING {
        for threads in ["1", "4"] {
            let root = temp_dir(case);
            fs::create_dir_all(root.join("d")).unwrap();
            // `rel.xml` is relative to the working directory (probe `naming/rel14`).
            let (cwd, out) = if case == "rel14" {
                (root.join("d"), out.to_string())
            } else {
                (root.clone(), out.to_string())
            };
            let mut args = blastp_args(query, subject, format, &out);
            let at = args.iter().position(|a| *a == "-num_threads").unwrap();
            args[at + 1] = threads;
            let output = run_in(&cwd, &args, &[]);
            assert_ok(&output, case);
            let dir = root.join("d");
            assert_eq!(
                names(&dir).into_iter().collect::<BTreeSet<_>>(),
                expected
                    .iter()
                    .map(|name| name.to_string())
                    .collect::<BTreeSet<_>>(),
                "{case}"
            );
            let extension = if format == "14" { "xml" } else { "json" };
            assert_eq!(
                fs::read(dir.join(expected[0])).unwrap(),
                read_frozen(&format!("master_{case}")),
                "{case}: the master file"
            );
            for n in 1..=2 {
                assert_eq!(
                    fs::read(dir.join(expected[n])).unwrap(),
                    read_frozen(&format!("naming_{n}.{extension}")),
                    "{case}: {}",
                    expected[n]
                );
            }
            fs::remove_dir_all(root).ok();
        }
    }
}

#[test]
fn blastn_writes_the_same_files() {
    let query = inputs().join("qnucl_valid.fna");
    let subject = inputs().join("snucl.fna");
    for (format, out, expected, extension, master) in [
        (
            "13",
            "d/Y.Json",
            ["Y.Json", "Y_1.json", "Y_2.json"],
            "json",
            "master_blastn13",
        ),
        (
            "14",
            "d/z",
            ["z", "z_1.xml", "z_2.xml"],
            "xml",
            "master_blastn14",
        ),
    ] {
        for threads in ["1", "4"] {
            let root = temp_dir("blastn");
            fs::create_dir_all(root.join("d")).unwrap();
            let output = run_in(
                &root,
                &[
                    "blastn",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                    "-out",
                    out,
                ],
                &[],
            );
            assert_ok(&output, out);
            let dir = root.join("d");
            assert_eq!(names(&dir), {
                let mut sorted = expected.map(String::from).to_vec();
                sorted.sort();
                sorted
            });
            assert_eq!(
                fs::read(dir.join(expected[0])).unwrap(),
                read_frozen(master)
            );
            for n in 1..=2 {
                assert_eq!(
                    fs::read(dir.join(expected[n])).unwrap(),
                    read_frozen(&format!("blastn_{n}.{extension}")),
                    "{out} {}",
                    expected[n]
                );
            }
            fs::remove_dir_all(root).ok();
        }
    }
}

/// Probe `naming/stale14`: a side file is opened for writing and cut; a file of another number,
/// left by an earlier run, stays as it was.
#[test]
fn an_existing_side_file_is_replaced_and_one_of_another_number_stays() {
    let root = temp_dir("stale");
    let dir = root.join("d");
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("s_1.xml"), "old\n").unwrap();
    fs::write(dir.join("s_3.xml"), "stale\n").unwrap();
    let query = frozen().join("qprot.faa");
    let subject = inputs().join("sprot.faa");
    let output = run_in(
        &root,
        &blastp_args(
            query.to_str().unwrap(),
            subject.to_str().unwrap(),
            "14",
            "d/s.xml",
        ),
        &[],
    );
    assert_ok(&output, "stale14");
    assert_eq!(
        fs::read(dir.join("s_1.xml")).unwrap(),
        read_frozen("naming_1.xml")
    );
    assert_eq!(
        fs::read(dir.join("s_2.xml")).unwrap(),
        read_frozen("naming_2.xml")
    );
    assert_eq!(fs::read(dir.join("s_3.xml")).unwrap(), b"stale\n");
    assert_eq!(
        fs::read(dir.join("s.xml")).unwrap(),
        read_frozen("master_stale14")
    );
    fs::remove_dir_all(root).ok();
}

/// Probes `naming/unwritable14`, `unwritable13`: a directory has the name of the second (14) or
/// the first (13) side file. NCBI: `Command line argument error: Cannot open output file`, exit
/// status 1; the side files before it stay, the master file exists and is empty.
#[test]
fn a_side_file_that_cannot_be_opened_stops_the_run_with_ncbis_message() {
    let query = frozen().join("qprot.faa");
    let subject = inputs().join("sprot.faa");
    for (format, out, blocked, written) in [
        (
            "14",
            "w/u.xml",
            "w/u_2.xml",
            Some(("w/u_1.xml", "naming_1.xml")),
        ),
        ("13", "w/v.json", "w/v_1.json", None),
    ] {
        let root = temp_dir("unwritable");
        fs::create_dir_all(root.join(blocked)).unwrap();
        let output = run_in(
            &root,
            &blastp_args(
                query.to_str().unwrap(),
                subject.to_str().unwrap(),
                format,
                out,
            ),
            &[],
        );
        assert_eq!(output.status.code(), Some(1), "{out}");
        assert!(output.stdout.is_empty(), "{out}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            "Command line argument error: Cannot open output file\n",
            "{out}"
        );
        assert_eq!(
            fs::read(root.join(out)).unwrap(),
            b"",
            "{out}: the master stays empty"
        );
        if let Some((file, frozen_name)) = written {
            assert_eq!(fs::read(root.join(file)).unwrap(), read_frozen(frozen_name));
        }
        assert!(root.join(blocked).is_dir(), "{out}");
        let expected_count = if written.is_some() { 3 } else { 2 };
        assert_eq!(names(&root.join("w")).len(), expected_count, "{out}");
        fs::remove_dir_all(root).ok();
    }
}

/// Probes `naming/missingdir14`, `missingdir16`: `-out` in a directory that does not exist is the
/// argument error of every format, before the search.
#[test]
fn an_output_in_a_missing_directory_is_an_argument_error_before_the_search() {
    let query = frozen().join("qprot.faa");
    let subject = inputs().join("sprot.faa");
    for format in ["13", "14", "15", "16"] {
        let root = temp_dir("missing");
        let output = run_in(
            &root,
            &blastp_args(
                query.to_str().unwrap(),
                subject.to_str().unwrap(),
                format,
                "nodir/x.xml",
            ),
            &[],
        );
        assert_eq!(output.status.code(), Some(1), "{format}");
        assert!(output.stdout.is_empty(), "{format}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            "Command line argument error: Argument \"out\". File is not accessible:  `nodir/x.xml'\n",
            "{format}"
        );
        assert_eq!(names(&root).len(), 0, "{format}: nothing is written");
        fs::remove_dir_all(root).ok();
    }
}

/// `UseXInclude`: outfmt 13 and 14 to standard output (no `-out`, or `-out -`) are NCBI's error
/// after `Query is Empty!`, before any search.
#[test]
fn outfmt_13_and_14_to_standard_output_ask_for_a_file_name() {
    let query = frozen().join("qprot.faa");
    let subject = inputs().join("sprot.faa");
    let root = temp_dir("stdout");
    for (format, number) in [("13", "13"), ("14", "14")] {
        for out in [None, Some("-")] {
            let mut args = vec![
                "blastp",
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                format,
            ];
            if let Some(out) = out {
                args.extend(["-out", out]);
            }
            let output = run_in(&root, &args, &[]);
            assert_eq!(output.status.code(), Some(1), "{format} {out:?}");
            assert!(output.stdout.is_empty());
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                format!("BLAST query/options error: Please provide a file name for outfmt {number}.\nPlease refer to the BLAST+ user manual.\n"),
                "{format} {out:?}"
            );
        }
    }
    assert_eq!(names(&root).len(), 0);
    fs::remove_dir_all(root).ok();
}

/// Probes `blastp_sub_emptypipe_fmt1[34]`: a run with no query formatted (an empty pipe) writes a
/// master file with no entry and no side file; `blastp_sub_emptyfile_fmt1[34]`: an empty query
/// file is `Query is Empty!` and an empty master file. (Bytes of NCBI, `-subject` and `-db`.)
#[test]
fn a_run_with_no_query_writes_a_master_with_no_entry() {
    let subject = inputs().join("sprot.faa");
    let db = inputs().join("db/sprot");
    let master14 = "<?xml version=\"1.0\"?>\n<BlastXML2\nxmlns=\"http://www.ncbi.nlm.nih.gov\"\nxmlns:xi=\"http://www.w3.org/2003/XInclude\"\nxmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\nxs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\">\n</BlastXML2>\n";
    let master13 = "{\n\t\"BlastJSON\": [\n\t]\n}";
    for (format, expected) in [("14", master14), ("13", master13)] {
        for target in [
            ["-subject", subject.to_str().unwrap()],
            ["-db", db.to_str().unwrap()],
        ] {
            for threads in ["1", "4"] {
                let root = temp_dir("emptypipe");
                let mut args = vec!["blastp", "-query", "-"];
                args.extend(target);
                args.extend(["-num_threads", threads, "-outfmt", format, "-out", "o.out"]);
                let output = run_with_stdin(&root, &args, b"");
                assert_ok(&output, &format!("{format} {target:?}"));
                assert_eq!(names(&root), ["o.out"], "{format} {target:?}");
                assert_eq!(
                    String::from_utf8_lossy(&fs::read(root.join("o.out")).unwrap()),
                    expected,
                    "{format} {target:?}"
                );
                fs::remove_dir_all(root).ok();
            }
        }
        // an empty query file
        let root = temp_dir("emptyfile");
        fs::write(root.join("empty.fa"), b"").unwrap();
        let output = run_in(
            &root,
            &[
                "blastp",
                "-query",
                "empty.fa",
                "-subject",
                subject.to_str().unwrap(),
                "-outfmt",
                format,
                "-out",
                "o.out",
            ],
            &[],
        );
        assert_eq!(output.status.code(), Some(0));
        assert!(output.stdout.is_empty());
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            "Warning: [blastp] Query is Empty!\n"
        );
        assert_eq!(names(&root), ["empty.fa", "o.out"]);
        assert_eq!(fs::read(root.join("o.out")).unwrap(), b"");
        fs::remove_dir_all(root).ok();
    }
}

/// Design §5.2 row 22 (`blastp_{sub,db}_badbatch_fmt1[34]`): with `BATCH_SIZE=100` the third
/// query's batch fails to read. `-db` has written the files of the two queries before it;
/// `-subject` only the first one's (the second is held back and lost); neither run has an epilog,
/// so the master file stays empty.
#[test]
fn a_failed_batch_keeps_the_db_files_and_loses_the_held_subject_file() {
    let dir = inputs();
    for format in ["13", "14"] {
        let extension = if format == "14" { "xml" } else { "json" };
        for (mode, target, files) in [
            ("sub", ["-subject", "sprot.faa"], 1),
            ("db", ["-db", "db/sprot"], 2),
        ] {
            for threads in ["1", "4"] {
                let work = temp_dir("badbatch");
                let out = work.join(format!("badbatch_{mode}_{format}.out"));
                let mut args = vec!["blastp", "-query", "qprot_bad.faa"];
                args.extend(target);
                args.extend([
                    "-num_threads",
                    threads,
                    "-outfmt",
                    format,
                    "-out",
                    out.to_str().unwrap(),
                ]);
                let output = run_in(&dir, &args, &[("BATCH_SIZE", "100")]);
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{mode} {format}: {stderr}");
                assert!(
                    stderr.starts_with("BLAST query error: CFastaReader: Near line 10"),
                    "{mode} {format}: {stderr}"
                );
                assert!(output.stdout.is_empty(), "{mode} {format}");
                let master = out.file_name().unwrap().to_string_lossy().into_owned();
                assert_eq!(
                    fs::read(&out).unwrap(),
                    b"",
                    "{mode} {format} threads {threads}: the master stays empty"
                );
                let mut expected_names = vec![master.clone()];
                for n in 1..=files {
                    let name = format!("{master}_{n}.{extension}");
                    assert_eq!(
                        fs::read(work.join(&name)).unwrap(),
                        read_frozen(&name),
                        "{mode} {format} threads {threads}: {name}"
                    );
                    expected_names.push(name);
                }
                expected_names.sort();
                assert_eq!(names(&work), expected_names, "{mode} {format}");
                fs::remove_dir_all(work).ok();
            }
        }
    }
}

/// 13 and 14 are written by the command line only (DECISIONS 26 Q10: `run_local` has no output
/// file, ABI v1 has no such format): the shared entries of BLASTP and BLASTN stop at NCBI's
/// XInclude error.
#[test]
fn run_local_does_not_write_outfmt_13_and_14() {
    use run_local_support::reader_records as records;
    use LOSAT::api::local_blast::{
        run_local_blastn, run_local_blastp, FormatOutput, OutputSink, ReportOutputs,
    };
    use LOSAT::blastinput::fasta_reader::ReaderConfig;
    use LOSAT::cli::{try_parse_from, Cli, Commands};
    let (query, subject) = (frozen().join("qprot.faa"), inputs().join("sprot.faa"));
    let (nquery, nsubject) = (inputs().join("qnucl_valid.fna"), inputs().join("snucl.fna"));
    for format in ["13", "14"] {
        let mut diagnostics = Vec::new();
        let mut sink = Vec::new();
        let mut outputs = ReportOutputs {
            formats: vec![FormatOutput {
                outfmt: format,
                sink: OutputSink::Writer(&mut sink),
            }],
            diagnostics: &mut diagnostics,
            hits: None,
            observer: None,
        };
        let cli: Cli = try_parse_from([
            "LOSAT",
            "blastp",
            "-query",
            query.to_str().unwrap(),
            "-subject",
            subject.to_str().unwrap(),
        ])
        .unwrap();
        let Commands::Blastp(args) = cli.command else {
            unreachable!()
        };
        let error = run_local_blastp(
            args,
            &records(&query, ReaderConfig::query("BLASTP", true, true)),
            &records(&subject, ReaderConfig::subject("BLASTP", true, true)),
            "",
            "",
            &mut outputs,
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains(&format!("Please provide a file name for outfmt {format}.")),
            "blastp {format}: {error}"
        );
        let cli: Cli = try_parse_from([
            "LOSAT",
            "blastn",
            "-query",
            nquery.to_str().unwrap(),
            "-subject",
            nsubject.to_str().unwrap(),
        ])
        .unwrap();
        let Commands::Blastn(args) = cli.command else {
            unreachable!()
        };
        let error = run_local_blastn(
            args,
            &records(&nquery, ReaderConfig::query("BLASTN", false, true)),
            &records(&nsubject, ReaderConfig::subject("BLASTN", false, true)),
            &mut outputs,
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains(&format!("Please provide a file name for outfmt {format}.")),
            "blastn {format}: {error}"
        );
        drop(outputs);
        assert!(sink.is_empty() && diagnostics.is_empty());
    }
}

// ---------------------------------------------------------------------------------------------
// OF5c part 2: BLASTX, TBLASTN and TBLASTX.
//
// The bytes of NCBI BLAST+ 2.17.0 are in `tests/fixtures/outfmt13_14/<case>.out[_<n>.xml|json]`; each case
// was run in `tests/fixtures/outfmt16` as `[BATCH_SIZE=100] PROGRAM -query QUERY -subject ... | -db ...
// -num_threads 1 -outfmt N -out <case>.out`:
//  * `badbatch_blastx_{sub,db}_N`: `qnx_bad.fna` (two queries of 10,003 bases, one batch each, then a
//    third whose batch fails to read; BLASTX reads batches of 10,002 bases and not `BATCH_SIZE`);
//  * `multibatch_blastx_{sub,db}_N`: `qnx_late.fna` (two batches, both read);
//  * `badbatch_tblastn_{sub,db}_N`: `BATCH_SIZE=100`, `qprot_bad.faa` against `snucl.fna` | `db/snucl`;
//  * `badbatch_tblastx_{sub,db}_N`: `BATCH_SIZE=100`, `qnucl_bad.fna`, `-evalue 1e-20 -max_hsps 2`;
//  * `holdmix_blastx_db_N`: `qnx_late.fna` against `db/sprot_mixtax` with `BLASTDB` the taxonomy
//    database `high_bytes_utf8` (DECISIONS 32: LOSAT holds the output of this run in memory).
// ---------------------------------------------------------------------------------------------

const BAD_BATCH: &str = "BLAST query error: CFastaReader: Near line";

/// The names of the frozen files of `case`: its master file and its side files, sorted.
fn frozen_names(case: &str) -> Vec<String> {
    let master = format!("{case}.out");
    let mut found: Vec<String> = names(&frozen())
        .into_iter()
        .filter(|name| *name == master || name.starts_with(&format!("{master}_")))
        .collect();
    found.sort();
    assert!(!found.is_empty(), "no frozen files for {case}");
    found
}

/// Runs `args` (with `-num_threads`, `-outfmt` and `-out <work>/<case>.out` added) in
/// `tests/fixtures/outfmt16` at one and four threads and checks that exactly the files of NCBI's run
/// appear, with its bytes, with `rc` and a standard error that starts as `stderr_start` (and is
/// empty when that is empty). Returns nothing; every run is in its own directory.
fn equals_frozen_files(
    case: &str,
    format: &str,
    args: &[&str],
    env: &[(&str, &str)],
    rc: i32,
    stderr_start: &str,
) {
    let expected = frozen_names(case);
    for threads in ["1", "4"] {
        let work = temp_dir(case);
        let out = work.join(format!("{case}.out"));
        let mut argv = args.to_vec();
        argv.extend([
            "-num_threads",
            threads,
            "-outfmt",
            format,
            "-out",
            out.to_str().unwrap(),
        ]);
        let output = run_in(&inputs(), &argv, env);
        let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
        let label = format!("{case} threads {threads}");
        assert_eq!(output.status.code(), Some(rc), "{label}: {stderr}");
        assert!(output.stdout.is_empty(), "{label}: standard output");
        if stderr_start.is_empty() {
            assert_eq!(stderr, "", "{label}");
        } else {
            assert!(stderr.starts_with(stderr_start), "{label}: {stderr}");
        }
        assert_eq!(names(&work), expected, "{label}: the files");
        for name in &expected {
            assert_eq!(
                fs::read(work.join(name)).unwrap(),
                read_frozen(name),
                "{label}: {name}"
            );
        }
        fs::remove_dir_all(work).ok();
    }
}

fn blastx_args<'a>(query: &'a str, target: [&'a str; 2]) -> Vec<&'a str> {
    let mut args = vec!["blastx", "-query", query];
    args.extend(target);
    args
}

/// BLASTX over several query batches (it reads batches of 10,002 bases): the master lists the side files of
/// every batch's queries, in query order, at every thread count.
#[test]
fn blastx_writes_the_files_of_several_query_batches() {
    for format in ["13", "14"] {
        for (mode, target) in [
            ("sub", ["-subject", "sprot.faa"]),
            ("db", ["-db", "db/sprot"]),
        ] {
            equals_frozen_files(
                &format!("multibatch_blastx_{mode}_{format}"),
                format,
                &blastx_args("qnx_late.fna", target),
                &[],
                0,
                "",
            );
        }
    }
}

/// Design §5.2 row 22 for BLASTX (`badbatch_blastx_*`): the third query's batch fails to read. `-db` has written
/// the files of the two queries before it; `-subject` only the first one's (the second is held back and lost);
/// neither run has an epilog, so the master stays empty.
#[test]
fn blastx_failed_batch_keeps_the_db_files_and_loses_the_held_subject_file() {
    for format in ["13", "14"] {
        for (mode, target) in [
            ("sub", ["-subject", "sprot.faa"]),
            ("db", ["-db", "db/sprot"]),
        ] {
            equals_frozen_files(
                &format!("badbatch_blastx_{mode}_{format}"),
                format,
                &blastx_args("qnx_bad.fna", target),
                &[],
                1,
                BAD_BATCH,
            );
        }
    }
}

/// The same for TBLASTN (`BATCH_SIZE=100`: the third query's batch fails).
#[test]
fn tblastn_failed_batch_keeps_the_db_files_and_loses_the_held_subject_file() {
    for format in ["13", "14"] {
        for (mode, target) in [
            ("sub", ["-subject", "snucl.fna"]),
            ("db", ["-db", "db/snucl"]),
        ] {
            let mut args = vec!["tblastn", "-query", "qprot_bad.faa"];
            args.extend(target);
            equals_frozen_files(
                &format!("badbatch_tblastn_{mode}_{format}"),
                format,
                &args,
                &[("BATCH_SIZE", "100")],
                1,
                BAD_BATCH,
            );
        }
    }
}

/// The same for TBLASTX (`-subject` has written no file at all: the first query's report is held back, too).
#[test]
fn tblastx_failed_batch_keeps_the_db_files_and_loses_the_held_subject_file() {
    for format in ["13", "14"] {
        for (mode, target) in [
            ("sub", ["-subject", "snucl.fna"]),
            ("db", ["-db", "db/snucl"]),
        ] {
            let mut args = vec!["tblastx", "-query", "qnucl_bad.fna"];
            args.extend(target);
            args.extend(["-evalue", "1e-20", "-max_hsps", "2"]);
            equals_frozen_files(
                &format!("badbatch_tblastx_{mode}_{format}"),
                format,
                &args,
                &[("BATCH_SIZE", "100")],
                1,
                BAD_BATCH,
            );
        }
    }
}

/// A first batch that fails to read (`firstbad` of the 15/16 tests; NCBI is the same in 13 and 14, cases run
/// against 2.17.0): no report is made, so no side file and an empty master file.
#[test]
fn a_first_batch_that_fails_to_read_leaves_an_empty_master_and_no_side_file() {
    for (program, query, target_sub, target_db) in [
        (
            "blastx",
            "qnx_firstbad.fna",
            ["-subject", "sprot.faa"],
            ["-db", "db/sprot"],
        ),
        (
            "tblastn",
            "qprot_bad.faa",
            ["-subject", "snucl.fna"],
            ["-db", "db/snucl"],
        ),
        (
            "tblastx",
            "qnucl_bad.fna",
            ["-subject", "snucl.fna"],
            ["-db", "db/snucl"],
        ),
    ] {
        for format in ["13", "14"] {
            for target in [target_sub, target_db] {
                let work = temp_dir("firstbad");
                let out = work.join("o.out");
                let mut args = vec![program, "-query", query];
                args.extend(target);
                args.extend(["-outfmt", format, "-out", out.to_str().unwrap()]);
                let output = run_in(&inputs(), &args, &[]);
                let stderr = String::from_utf8_lossy(&output.stderr);
                let label = format!("{program} {format} {target:?}");
                assert_eq!(output.status.code(), Some(1), "{label}: {stderr}");
                assert!(stderr.starts_with(BAD_BATCH), "{label}: {stderr}");
                assert!(output.stdout.is_empty(), "{label}");
                assert_eq!(names(&work), ["o.out"], "{label}");
                assert_eq!(fs::read(&out).unwrap(), b"", "{label}");
                fs::remove_dir_all(work).ok();
            }
        }
    }
}

/// Design 4.5: TBLASTN `-matrix blosum62` (lower case) with queries (probe `tblastn_sub_lcmatrix_fmt1[34]`): the
/// report data constructor throws after the formatter is made, so rc 3, an empty master and no side file.
#[test]
fn tblastn_matrix_error_leaves_an_empty_master_and_no_side_file() {
    for format in ["13", "14"] {
        for target in [["-subject", "snucl.fna"], ["-db", "db/snucl"]] {
            let work = temp_dir("lcmatrix");
            let out = work.join("o.out");
            let mut args = vec!["tblastn", "-query", "qprot_unsearched.faa"];
            args.extend(target);
            args.extend([
                "-matrix",
                "blosum62",
                "-outfmt",
                format,
                "-out",
                out.to_str().unwrap(),
            ]);
            let output = run_in(&inputs(), &args, &[]);
            let label = format!("{format} {target:?}");
            assert_eq!(output.status.code(), Some(3), "{label}");
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                "BLAST engine error: unsupported score matrix\n",
                "{label}"
            );
            assert!(output.stdout.is_empty(), "{label}");
            assert_eq!(names(&work), ["o.out"], "{label}");
            assert_eq!(fs::read(&out).unwrap(), b"", "{label}");
            fs::remove_dir_all(work).ok();
        }
    }
}

/// BLASTX's hold (DECISIONS 21/32) keeps the side files' bytes too: a `-db` search whose taxonomy database
/// has a name LOSAT would refuse for a sequence never reported is held until the epilog, and then writes what
/// NCBI writes (`holdmix_blastx_db_N`: NCBI never looks that taxid up).
#[test]
fn blastx_hold_writes_the_files_at_the_epilog_as_ncbi_does() {
    let taxdb = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/high_bytes_utf8");
    for format in ["13", "14"] {
        equals_frozen_files(
            &format!("holdmix_blastx_db_{format}"),
            format,
            &blastx_args("qnx_late.fna", ["-db", "db/sprot_mixtax"]),
            &[("BLASTDB", taxdb.to_str().unwrap())],
            0,
            "",
        );
    }
}

/// An unopenable side file under the hold: the held side files are written at the epilog, in order, and the
/// run stops at the one it cannot open with NCBI's message (the files before it stay; the master stays empty).
#[test]
fn blastx_hold_stops_at_a_side_file_it_cannot_open() {
    let taxdb = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/high_bytes_utf8");
    for (format, extension) in [("13", "json"), ("14", "xml")] {
        for (blocked, written) in [(2, 1), (1, 0)] {
            let work = temp_dir("holdblocked");
            let out = work.join("holdmix_blastx_db.out");
            fs::create_dir(work.join(format!("holdmix_blastx_db.out_{blocked}.{extension}")))
                .unwrap();
            let mut args = blastx_args("qnx_late.fna", ["-db", "db/sprot_mixtax"]);
            args.extend(["-outfmt", format, "-out", out.to_str().unwrap()]);
            let output = run_in(&inputs(), &args, &[("BLASTDB", taxdb.to_str().unwrap())]);
            let label = format!("{format} blocked {blocked}");
            assert_eq!(output.status.code(), Some(1), "{label}");
            assert!(output.stdout.is_empty(), "{label}");
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                "Command line argument error: Cannot open output file\n",
                "{label}"
            );
            assert_eq!(
                fs::read(&out).unwrap(),
                b"",
                "{label}: the master stays empty"
            );
            let mut expected = vec!["holdmix_blastx_db.out".to_string()];
            expected.push(format!("holdmix_blastx_db.out_{blocked}.{extension}"));
            for n in 1..=written {
                let name = format!("holdmix_blastx_db.out_{n}.{extension}");
                assert_eq!(
                    fs::read(work.join(&name)).unwrap(),
                    read_frozen(&format!("holdmix_blastx_db_{format}.out_{n}.{extension}")),
                    "{label}: {name}"
                );
                expected.push(name);
            }
            expected.sort();
            assert_eq!(names(&work), expected, "{label}");
            fs::remove_dir_all(work).ok();
        }
    }
}

/// DECISIONS 28/30 for 13 and 14: a taxonomy database name outside US-ASCII that the second batch reports is
/// refused after the first batch's report was made: the hold drops it, so no side file is written and the master
/// stays empty. (NCBI cuts the side file of the failing report at the byte, rc 255, in 14: a refusal is LOSAT's.)
#[test]
fn a_refusal_in_a_later_batch_leaves_no_side_file_and_an_empty_master() {
    let taxdb = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/high_bytes_utf8");
    for format in ["13", "14"] {
        for threads in ["1", "4"] {
            let work = temp_dir("taxrefused");
            let out = work.join("o.out");
            let mut args = blastx_args("qnx_late.fna", ["-db", "db/sprot_tax"]);
            args.extend([
                "-num_threads",
                threads,
                "-outfmt",
                format,
                "-out",
                out.to_str().unwrap(),
            ]);
            let output = run_in(&inputs(), &args, &[("BLASTDB", taxdb.to_str().unwrap())]);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{format} {threads}: {stderr}"
            );
            assert!(
                stderr.contains(
                    "a taxonomy database scientific name with a character outside the US-ASCII range"
                ),
                "{format} {threads}: {stderr}"
            );
            assert!(output.stdout.is_empty());
            assert_eq!(names(&work), ["o.out"], "{format} {threads}");
            assert_eq!(fs::read(&out).unwrap(), b"", "{format} {threads}");
            fs::remove_dir_all(work).ok();
        }
    }
}

fn write_fasta(path: &Path, records: &[(&str, String)]) {
    let mut text = String::new();
    for (title, sequence) in records {
        text.push_str(&format!(">{title}\n{sequence}\n"));
    }
    fs::write(path, text).unwrap();
}

/// Q2 for 13 and 14 (as `a_title_or_database_name_outside_us_ascii_is_refused_in_15_and_16`): a byte >= 0x80 in a
/// query title (also in a later batch, which the pre-scan finds), a `-subject` title or the `-db` name is refused
/// before anything is written: the master is empty and there is no side file.
#[test]
fn a_title_outside_us_ascii_is_refused_in_13_and_14_with_nothing_written() {
    let protein = "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
    let codons: std::collections::HashMap<char, &str> = [
        ('A', "GCT"),
        ('C', "TGC"),
        ('D', "GAT"),
        ('E', "GAA"),
        ('F', "TTC"),
        ('G', "GGT"),
        ('H', "CAC"),
        ('I', "ATC"),
        ('K', "AAA"),
        ('L', "CTG"),
        ('M', "ATG"),
        ('N', "AAC"),
        ('P', "CCT"),
        ('Q', "CAG"),
        ('R', "CGT"),
        ('S', "TCT"),
        ('T', "ACC"),
        ('V', "GTG"),
        ('W', "TGG"),
        ('Y', "TAC"),
    ]
    .into_iter()
    .collect();
    let nucleotide: String = protein.chars().map(|c| codons[&c]).collect();
    let mut padded = String::new();
    while padded.len() < 10_003 {
        padded.push_str("ACGTTGCAAGCTTGCA");
    }
    let refused_title = |stderr: &str, what: &str| {
        stderr.contains(&format!(
            "{what} with a character outside the US-ASCII range (which the XML output of NCBI BLAST+ cannot write) is not supported by LOSAT's BLASTX"
        ))
    };
    for format in ["13", "14"] {
        for threads in ["1", "4"] {
            let work = temp_dir("q2");
            let (query, subject) = (work.join("q.fa"), work.join("s.fa"));
            let out = work.join("o.out");
            let run = || {
                run_in(
                    &work,
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
                        "-out",
                        out.to_str().unwrap(),
                    ],
                    &[],
                )
            };
            let check = |output: &Output, what: &str, label: &str| {
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{label}: {stderr}");
                assert!(refused_title(&stderr, what), "{label}: {stderr}");
                assert!(output.stdout.is_empty(), "{label}");
                assert_eq!(names(&work), ["o.out", "q.fa", "s.fa"], "{label}");
                assert_eq!(fs::read(&out).unwrap(), b"", "{label}: the master is empty");
            };
            write_fasta(&subject, &[("s1 subject", protein.to_string())]);
            // a query title: the first query's, or the second's in the next batch
            for (first, second, label) in [
                ("caf\u{e9} q1", "q2 ok", "first"),
                ("q1 ok", "caf\u{e9} q2", "second, same batch"),
            ] {
                write_fasta(
                    &query,
                    &[(first, nucleotide.clone()), (second, nucleotide.clone())],
                );
                check(
                    &run(),
                    "a query title",
                    &format!("{format} {threads} {label}"),
                );
            }
            write_fasta(
                &query,
                &[
                    ("q1 ok", padded.clone()),
                    ("caf\u{e9} q2", nucleotide.clone()),
                ],
            );
            check(
                &run(),
                "a query title",
                &format!("{format} {threads} later batch"),
            );
            // a subject title with a hit, in a batch after one without hits
            write_fasta(
                &query,
                &[("q1 no hit", padded.clone()), ("q2 ok", nucleotide.clone())],
            );
            write_fasta(
                &subject,
                &[
                    (
                        "s1 ok",
                        "WQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNM".to_string(),
                    ),
                    ("s2 caf\u{e9}", protein.to_string()),
                ],
            );
            check(
                &run(),
                "a subject title",
                &format!("{format} {threads} subject"),
            );
            fs::remove_dir_all(work).ok();
        }
    }
}

/// DECISIONS 22 for 13 and 14: BLASTX reads a query that is not a regular file once only, so it is refused before
/// the search with nothing written (the master file is created empty before: NCBI opens `-out` with its
/// arguments).
#[cfg(unix)]
#[test]
fn a_blastx_query_that_is_not_a_regular_file_is_refused_in_13_and_14() {
    use std::io::Write;
    let work = temp_dir("pipe");
    let text = b">q1 plain\nATGGATTTCGAACAGCGGCATGAACCGTAA\n".to_vec();
    for format in ["13", "14"] {
        let out = work.join(format!("o{format}.out"));
        let mut child = command(&[])
            .args([
                "blastx",
                "-query",
                "/dev/stdin",
                "-subject",
                inputs().join("sprot.faa").to_str().unwrap(),
                "-outfmt",
                format,
                "-out",
                out.to_str().unwrap(),
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
        assert_eq!(output.status.code(), Some(1), "{format}: {stderr}");
        assert!(output.stdout.is_empty(), "{format}");
        assert!(
            stderr.contains(&format!(
                "unsupported BLASTX query input for -outfmt {format}: use a regular file"
            )),
            "{format}: {stderr}"
        );
        assert_eq!(
            fs::read(&out).unwrap(),
            b"",
            "{format}: the master is empty"
        );
    }
    assert_eq!(names(&work), ["o13.out", "o14.out"]);
    fs::remove_dir_all(work).ok();
}

/// An empty query file is `Query is Empty!` and an empty master; an empty pipe (TBLASTN, TBLASTX) a master with
/// no entry (run against NCBI BLAST+ 2.17.0 for the OF5c part 2 report). BLASTX does not read a query from
/// standard input at all (its own message for every format), so it has no pipe case.
#[test]
fn an_empty_query_gives_an_empty_master_or_a_master_with_no_entry_in_the_other_programs() {
    let master14 = "<?xml version=\"1.0\"?>\n<BlastXML2\nxmlns=\"http://www.ncbi.nlm.nih.gov\"\nxmlns:xi=\"http://www.w3.org/2003/XInclude\"\nxmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"\nxs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\">\n</BlastXML2>\n";
    let master13 = "{\n\t\"BlastJSON\": [\n\t]\n}";
    let subject = inputs().join("snucl.fna");
    let db = inputs().join("db/snucl");
    for program in ["tblastn", "tblastx", "blastx"] {
        let (subject, db) = if program == "blastx" {
            (inputs().join("sprot.faa"), inputs().join("db/sprot"))
        } else {
            (subject.clone(), db.clone())
        };
        for (format, expected) in [("14", master14), ("13", master13)] {
            for target in [
                ["-subject", subject.to_str().unwrap()],
                ["-db", db.to_str().unwrap()],
            ] {
                for threads in ["1", "4"] {
                    let root = temp_dir("emptypipe2");
                    let mut args = vec![program, "-query", "-"];
                    args.extend(target);
                    args.extend(["-num_threads", threads, "-outfmt", format, "-out", "o.out"]);
                    let output = run_with_stdin(&root, &args, b"");
                    let label = format!("{program} {format} {target:?}");
                    if program == "blastx" {
                        assert_eq!(output.status.code(), Some(1), "{label}");
                        assert!(
                            String::from_utf8_lossy(&output.stderr)
                                .contains("unsupported BLASTX stdin query/subject"),
                            "{label}"
                        );
                    } else {
                        assert_ok(&output, &label);
                        assert_eq!(names(&root), ["o.out"], "{label}");
                        assert_eq!(
                            String::from_utf8_lossy(&fs::read(root.join("o.out")).unwrap()),
                            expected,
                            "{label}"
                        );
                    }
                    fs::remove_dir_all(root).ok();
                }
            }
            // an empty query file
            let root = temp_dir("emptyfile2");
            fs::write(root.join("empty.fa"), b"").unwrap();
            let output = run_in(
                &root,
                &[
                    program,
                    "-query",
                    "empty.fa",
                    "-subject",
                    subject.to_str().unwrap(),
                    "-outfmt",
                    format,
                    "-out",
                    "o.out",
                ],
                &[],
            );
            assert_eq!(output.status.code(), Some(0), "{program} {format}");
            assert!(output.stdout.is_empty());
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                format!("Warning: [{program}] Query is Empty!\n")
            );
            assert_eq!(names(&root), ["empty.fa", "o.out"]);
            assert_eq!(fs::read(root.join("o.out")).unwrap(), b"");
            fs::remove_dir_all(root).ok();
        }
    }
}

/// 13 and 14 to standard output for BLASTX, TBLASTN and TBLASTX are NCBI's error as for the other programs
/// (`cli_v2.rs` checks the message of each program; here: nothing is written and no search starts).
#[test]
fn outfmt_13_and_14_to_standard_output_ask_for_a_file_name_in_every_program() {
    let root = temp_dir("stdout2");
    for (program, query, target) in [
        ("blastx", "qnx_late.fna", ["-subject", "sprot.faa"]),
        ("tblastn", "qprot_bad.faa", ["-subject", "snucl.fna"]),
        ("tblastx", "qnucl_valid.fna", ["-subject", "snucl.fna"]),
    ] {
        for format in ["13", "14"] {
            for out in [None, Some("-")] {
                let mut args = vec![program, "-query", query];
                args.extend(target);
                args.extend(["-outfmt", format]);
                if let Some(out) = out {
                    args.extend(["-out", out]);
                }
                let output = run_in(&inputs(), &args, &[]);
                let label = format!("{program} {format} {out:?}");
                assert_eq!(output.status.code(), Some(1), "{label}");
                assert!(output.stdout.is_empty(), "{label}");
                assert!(
                    String::from_utf8_lossy(&output.stderr)
                        .contains(&format!("Please provide a file name for outfmt {format}.")),
                    "{label}: {}",
                    String::from_utf8_lossy(&output.stderr)
                );
            }
        }
    }
    assert_eq!(names(&root).len(), 0);
    fs::remove_dir_all(root).ok();
}

/// `run_local` has no output file, so TBLASTN and TBLASTX (like BLASTP and BLASTN) stop at NCBI's XInclude error.
#[test]
fn run_local_does_not_write_outfmt_13_and_14_in_tblastn_and_tblastx() {
    use run_local_support::reader_records as records;
    use LOSAT::api::local_blast::{
        run_local_tblastn, run_local_tblastx, FormatOutput, OutputSink, ReportOutputs,
    };
    use LOSAT::blastinput::fasta_reader::ReaderConfig;
    use LOSAT::cli::{try_parse_from, Cli, Commands};
    let (pquery, nsubject) = (inputs().join("qprot_bad.faa"), inputs().join("snucl.fna"));
    let nquery = inputs().join("qnucl_valid.fna");
    let pquery = {
        // qprot_bad.faa has a bad third record: use the first two
        let valid = temp_dir("rl").join("qprot_valid.faa");
        let text = fs::read_to_string(&pquery).unwrap();
        fs::write(&valid, text.split(">q3").next().unwrap()).unwrap();
        valid
    };
    for format in ["13", "14"] {
        for program in ["tblastn", "tblastx"] {
            let mut diagnostics = Vec::new();
            let mut sink = Vec::new();
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: format,
                    sink: OutputSink::Writer(&mut sink),
                }],
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            let query = if program == "tblastn" {
                &pquery
            } else {
                &nquery
            };
            let cli: Cli = try_parse_from([
                "LOSAT",
                program,
                "-query",
                query.to_str().unwrap(),
                "-subject",
                nsubject.to_str().unwrap(),
            ])
            .unwrap();
            let error = match cli.command {
                Commands::Tblastn(args) => run_local_tblastn(
                    args,
                    &records(query, ReaderConfig::query("TBLASTN", true, true)),
                    &records(&nsubject, ReaderConfig::subject("TBLASTN", false, true)),
                    &mut outputs,
                ),
                Commands::Tblastx(args) => run_local_tblastx(
                    args,
                    &records(query, ReaderConfig::query("TBLASTX", false, true)),
                    &records(&nsubject, ReaderConfig::subject("TBLASTX", false, true)),
                    &mut outputs,
                ),
                _ => unreachable!(),
            }
            .unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains(&format!("Please provide a file name for outfmt {format}.")),
                "{program} {format}: {error}"
            );
            drop(outputs);
            assert!(sink.is_empty() && diagnostics.is_empty());
        }
    }
}

// ---- OF5c review fixes (DECISIONS 35) ----

/// The inputs of the probes `W1` of the OF5c review: three query batches (BLASTX reads batches of 10,002 bases) of
/// which the third carries an invalid residue (a reader warning), and a `-subject` file with a title outside
/// US-ASCII on a sequence without a hit, which makes BLASTX hold its output.
fn write_hold_inputs(work: &Path, third_query: &str) -> (PathBuf, PathBuf) {
    let protein = "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
    let codons: std::collections::HashMap<char, &str> = [
        ('A', "GCT"),
        ('C', "TGC"),
        ('D', "GAT"),
        ('E', "GAA"),
        ('F', "TTC"),
        ('G', "GGT"),
        ('H', "CAC"),
        ('I', "ATC"),
        ('K', "AAA"),
        ('L', "CTG"),
        ('M', "ATG"),
        ('N', "AAC"),
        ('P', "CCT"),
        ('Q', "CAG"),
        ('R', "CGT"),
        ('S', "TCT"),
        ('T', "ACC"),
        ('V', "GTG"),
        ('W', "TGG"),
        ('Y', "TAC"),
    ]
    .into_iter()
    .collect();
    let nucleotide: String = protein.chars().map(|c| codons[&c]).collect();
    let mut padded = String::new();
    while padded.len() < 10_003 {
        padded.push_str("ACGTTGCAAGCTTGCA");
    }
    let (query, subject) = (work.join("q.fa"), work.join("s.fa"));
    write_fasta(
        &query,
        &[
            ("q1 no hit", padded.clone()),
            ("q2 hit", format!("{nucleotide}{padded}")),
            ("q3 later", third_query.to_string()),
        ],
    );
    write_fasta(
        &subject,
        &[
            ("s1 subject", protein.to_string()),
            (
                "s2 caf\u{e9} no hit",
                "WQERTYIPASDFGHKLCVNMWQERTYIPASDFGHKLCVNM".to_string(),
            ),
        ],
    );
    (query, subject)
}

/// DECISIONS 35 S1 (OF5c review, probe `W1` and its later-batch-failure variant): a run that holds its side files
/// opens each one when NCBI does, so a directory in the place of side file 1 stops the run when the first report
/// is written (a `-subject` run writes it when the second arrives), with NCBI's message alone: the third batch is
/// not read, so neither its reader warning nor its failure comes first. The directory stays and the master is empty.
#[test]
fn blastx_hold_stops_at_an_unopenable_side_file_before_it_reads_the_later_batches() {
    let warned = format!("{}ACGTJJACGT{}", "ACGT".repeat(2400), "ACGT".repeat(100));
    for (format, extension) in [("13", "json"), ("14", "xml")] {
        for (label, third) in [("warning", warned.as_str()), ("bad batch", "ACGT%%%%")] {
            for threads in ["1", "4"] {
                let work = temp_dir("holdearly");
                let (query, subject) = write_hold_inputs(&work, third);
                let out = work.join("o.out");
                let side = work.join(format!("o.out_1.{extension}"));
                fs::create_dir(&side).unwrap();
                let output = run_in(
                    &work,
                    &[
                        "blastx",
                        "-query",
                        query.to_str().unwrap(),
                        "-subject",
                        subject.to_str().unwrap(),
                        "-evalue",
                        "1e-5",
                        "-num_threads",
                        threads,
                        "-outfmt",
                        format,
                        "-out",
                        out.to_str().unwrap(),
                    ],
                    &[],
                );
                let what = format!("{format} {label} {threads}");
                assert_eq!(output.status.code(), Some(1), "{what}");
                assert!(output.stdout.is_empty(), "{what}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr),
                    "Command line argument error: Cannot open output file\n",
                    "{what}"
                );
                assert_eq!(fs::read(&out).unwrap(), b"", "{what}: the master is empty");
                assert!(side.is_dir(), "{what}");
                assert_eq!(
                    names(&work),
                    [
                        "o.out",
                        format!("o.out_1.{extension}").as_str(),
                        "q.fa",
                        "s.fa"
                    ],
                    "{what}"
                );
                fs::remove_dir_all(work).ok();
            }
        }
    }
}

/// The same run with nothing in the way writes the two files and the master, the third query's warning
/// included: the early check leaves no file of its own and changes no byte.
#[test]
fn blastx_hold_with_free_side_files_writes_them_all_and_the_check_leaves_nothing_else() {
    let warned = format!("{}ACGTJJACGT{}", "ACGT".repeat(2400), "ACGT".repeat(100));
    for (format, extension) in [("13", "json"), ("14", "xml")] {
        let work = temp_dir("holdfree");
        let (query, subject) = write_hold_inputs(&work, &warned);
        let out = work.join("o.out");
        // an existing side file is replaced; a higher number stays
        fs::write(
            work.join(format!("o.out_1.{extension}")),
            b"stale stale stale",
        )
        .unwrap();
        fs::write(work.join(format!("o.out_9.{extension}")), b"stale").unwrap();
        let output = run_in(
            &work,
            &[
                "blastx",
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
                "-evalue",
                "1e-5",
                "-num_threads",
                "1",
                "-outfmt",
                format,
                "-out",
                out.to_str().unwrap(),
            ],
            &[],
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "{format}: {stderr}");
        assert!(
            stderr.contains("Ignoring invalid residues"),
            "{format}: {stderr}"
        );
        let mut expected = vec![
            "o.out".to_string(),
            format!("o.out_1.{extension}"),
            format!("o.out_2.{extension}"),
            format!("o.out_3.{extension}"),
            format!("o.out_9.{extension}"),
            "q.fa".to_string(),
            "s.fa".to_string(),
        ];
        expected.sort();
        assert_eq!(names(&work), expected, "{format}");
        assert_ne!(
            fs::read(work.join(format!("o.out_1.{extension}"))).unwrap(),
            b"stale stale stale"
        );
        assert_eq!(
            fs::read(work.join(format!("o.out_9.{extension}"))).unwrap(),
            b"stale"
        );
        fs::remove_dir_all(work).ok();
    }
}

/// DECISIONS 35 S2: a write to a side file that fails after it was opened (the file is a link to `/dev/full`)
/// ends the run with one answer in all five programs and in BLASTX with and without a hold: `BLAST failed to write
/// output`, exit status 6 (the message and status of a failed write to `-out` in outfmt 0), an empty master, and the
/// side files after the failed one not written. Deviation: NCBI reports the failed flush and goes on with exit
/// status 0 (OF5c review, probe `F`). Linux only: macOS and Windows have no `/dev/full`.
#[cfg(target_os = "linux")]
#[test]
fn a_side_file_write_error_gives_one_answer_in_every_program() {
    let query = frozen().join("qprot.faa");
    let sprot = inputs().join("sprot.faa");
    let qnucl = inputs().join("qnucl_valid.fna");
    let snucl = inputs().join("snucl.fna");
    let db_snucl = inputs().join("db/snucl");
    let db_sprot = inputs().join("db/sprot");
    let (query, sprot, qnucl, snucl, db_snucl, db_sprot) = (
        query.to_str().unwrap().to_string(),
        sprot.to_str().unwrap().to_string(),
        qnucl.to_str().unwrap().to_string(),
        snucl.to_str().unwrap().to_string(),
        db_snucl.to_str().unwrap().to_string(),
        db_sprot.to_str().unwrap().to_string(),
    );
    let strs = |args: &[&str]| args.iter().map(|a| a.to_string()).collect::<Vec<_>>();
    // (label, arguments before -outfmt, whether the run holds)
    let cases: Vec<(&str, Vec<String>, bool)> = vec![
        (
            "blastp -subject",
            strs(&["blastp", "-query", &query, "-subject", &sprot]),
            false,
        ),
        (
            "blastp -db",
            strs(&["blastp", "-query", &query, "-db", &db_sprot]),
            false,
        ),
        (
            "blastn -subject",
            strs(&["blastn", "-query", &qnucl, "-subject", &snucl]),
            false,
        ),
        (
            "blastn -db",
            strs(&["blastn", "-query", &qnucl, "-db", &db_snucl]),
            false,
        ),
        (
            "tblastn -subject",
            strs(&["tblastn", "-query", &query, "-subject", &snucl]),
            false,
        ),
        (
            "tblastn -db",
            strs(&["tblastn", "-query", &query, "-db", &db_snucl]),
            false,
        ),
        (
            "tblastx -subject",
            strs(&[
                "tblastx", "-query", &qnucl, "-subject", &snucl, "-evalue", "1e-20",
            ]),
            false,
        ),
        (
            "tblastx -db",
            strs(&[
                "tblastx", "-query", &qnucl, "-db", &db_snucl, "-evalue", "1e-20",
            ]),
            false,
        ),
        (
            "blastx -subject",
            strs(&[
                "blastx",
                "-query",
                inputs().join("qnx_late.fna").to_str().unwrap(),
                "-subject",
                &sprot,
                "-evalue",
                "1e-5",
            ]),
            false,
        ),
        (
            "blastx -db",
            strs(&[
                "blastx",
                "-query",
                inputs().join("qnx_late.fna").to_str().unwrap(),
                "-db",
                &db_sprot,
            ]),
            false,
        ),
        // the hold: the subject file has a title outside US-ASCII on a sequence without a hit
        ("blastx -subject, hold", Vec::new(), true),
    ];
    for (format, extension) in [("13", "json"), ("14", "xml")] {
        for (label, args, hold) in &cases {
            let work = temp_dir("devfull");
            let mut args = args.clone();
            if *hold {
                let (query, subject) = write_hold_inputs(&work, &"ACGT".repeat(2600));
                args = strs(&[
                    "blastx",
                    "-query",
                    query.to_str().unwrap(),
                    "-subject",
                    subject.to_str().unwrap(),
                    "-evalue",
                    "1e-5",
                ]);
            }
            let out = work.join("o.out");
            std::os::unix::fs::symlink("/dev/full", work.join(format!("o.out_1.{extension}")))
                .unwrap();
            args.extend(strs(&[
                "-num_threads",
                "1",
                "-outfmt",
                format,
                "-out",
                out.to_str().unwrap(),
            ]));
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let output = run_in(&inputs(), &args, &[]);
            let what = format!("{label} {format}");
            assert_eq!(
                String::from_utf8_lossy(&output.stderr),
                "BLAST failed to write output\n",
                "{what}"
            );
            assert_eq!(output.status.code(), Some(6), "{what}");
            assert!(output.stdout.is_empty(), "{what}");
            assert_eq!(fs::read(&out).unwrap(), b"", "{what}: the master is empty");
            assert!(
                !work.join(format!("o.out_2.{extension}")).exists(),
                "{what}: no side file after the failed one"
            );
            fs::remove_dir_all(work).ok();
        }
    }
}
