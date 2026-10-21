#![allow(warnings, clippy::all)]

//! BLAST JSON and XML2 as files (outfmt 13 and 14) of BLASTP and BLASTN: what the frozen-fixture
//! rows of the manifests (`files` column) cannot show. The names of the side files and of the entries of
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
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if (name_text.starts_with("LOSAT_") && name_text != "LOSAT_X_DEFAULT")
            || name_text.starts_with("NCBI")
        {
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
