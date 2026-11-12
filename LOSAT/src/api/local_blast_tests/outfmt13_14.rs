//! `run_local` without an output file (outfmt 13 and 14) in process: the entry's own contract
//! when the caller gives no file (an error before the search) or a receiver of the side files
//! (`OutputSink::Files`, LOSAT Web). What the command line writes for these formats is checked in
//! `LOSAT/tests/outfmt13_14_files.rs`.

use std::fs;
use std::path::{Path, PathBuf};

use crate::api::{read_input_file, read_inputs, run_to_file, validate, Inputs, Program, Role};

fn engine() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf()
}

fn inputs() -> PathBuf {
    engine().join("tests/fixtures/outfmt16")
}

fn frozen() -> PathBuf {
    engine().join("tests/fixtures/outfmt13_14")
}

fn temp_dir(stem: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "losat_outfmt1314_crate_{stem}_{}_{}",
        std::process::id(),
        nanos
    ));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// A `run_local` call without an output file and without a receiver of side files
/// (`OutputSink::Files`) cannot write 13 and 14: the shared entries of BLASTP and BLASTN stop at
/// NCBI's XInclude error.
#[test]
fn run_local_does_not_write_outfmt_13_and_14() {
    use super::support::reader_records as records;
    use crate::api::local_blast::{
        run_local_blastn, run_local_blastp, FormatOutput, OutputSink, ReportOutputs,
    };
    use crate::blastinput::fasta_reader::ReaderConfig;
    use crate::cli::{try_parse_from, Cli, Commands};
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

/// Without an output file or a receiver of side files, TBLASTN and TBLASTX (like BLASTP and BLASTN) stop at
/// NCBI's XInclude error.
#[test]
fn run_local_does_not_write_outfmt_13_and_14_in_tblastn_and_tblastx() {
    use super::support::reader_records as records;
    use crate::api::local_blast::{
        run_local_tblastn, run_local_tblastx, FormatOutput, OutputSink, ReportOutputs,
    };
    use crate::blastinput::fasta_reader::ReaderConfig;
    use crate::cli::{try_parse_from, Cli, Commands};
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

/// A receiver of side files that keeps them in memory (`SideFileReceiver`).
#[derive(Default)]
struct Received(std::sync::Mutex<Vec<(String, Vec<u8>)>>);

impl crate::report::blast_xml2::SideFileReceiver for Received {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
        self.0
            .lock()
            .unwrap()
            .push((String::from_utf8(name.to_vec()).unwrap(), bytes.to_vec()));
        Ok(())
    }
}

/// The files of `api::run_to_file` for outfmt `format` of `program` with `-out <dir>/r`: its
/// warnings, the master file and the side files `r_<n>.<ext>` in the order of their numbers.
fn reference(
    program: Program,
    query: &Path,
    subject: &Path,
    extra: &[&str],
    format: &str,
    dir: &Path,
    ext: &str,
) -> (Vec<u8>, Vec<u8>, Vec<(String, Vec<u8>)>) {
    let mut options = extra.to_vec();
    options.extend_from_slice(&["-num_threads", "1"]);
    let search = validate(program, &options).unwrap();
    let (queries, subjects) = read_inputs(
        &search,
        &read_input_file(Role::Query, query).unwrap(),
        &read_input_file(Role::Subject, subject).unwrap(),
    )
    .unwrap();
    let (query_name, subject_name) = (query.display().to_string(), subject.display().to_string());
    let mut diagnostics = Vec::new();
    run_to_file(
        &search,
        Inputs::new(&queries, &subjects).named(&query_name, &subject_name),
        format,
        &dir.join("r"),
        &mut diagnostics,
    )
    .unwrap();
    let sides = (1..)
        .map_while(|n| {
            let name = format!("r_{n}.{ext}");
            fs::read(dir.join(&name)).ok().map(|bytes| (name, bytes))
        })
        .collect();
    (diagnostics, fs::read(dir.join("r")).unwrap(), sides)
}

/// `run_local` with the caller's receiver of side files (`OutputSink::Files`, LOSAT Web) writes
/// outfmt 13 and 14 as `api::run_to_file` writes them with `-out r` in every program: the master
/// file to the writer, each side file whole to the receiver, by the names of the files, and the
/// warnings to the diagnostics. The test of the same name in `LOSAT/tests/outfmt13_14_files.rs`
/// holds `run_to_file`'s files equal to the command line's.
#[test]
fn run_local_writes_outfmt_13_and_14_to_a_receiver_as_the_command_line_writes_files() {
    use super::support::reader_records as records;
    use crate::api::local_blast::{
        run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
        FormatOutput, OutputSink, ReportOutputs,
    };
    use crate::blastinput::fasta_reader::ReaderConfig;
    use crate::cli::{try_parse_from, Cli, Commands};
    use crate::report::blast_xml2::SideFileSink;
    use std::sync::Arc;
    let work = temp_dir("receiver");
    let pquery = {
        // qprot_bad.faa has a bad third record: use the first two
        let valid = work.join("qprot_valid.faa");
        let text = fs::read_to_string(inputs().join("qprot_bad.faa")).unwrap();
        fs::write(&valid, text.split(">q3").next().unwrap()).unwrap();
        valid
    };
    let (sprot, snucl) = (inputs().join("sprot.faa"), inputs().join("snucl.fna"));
    let (nquery, xquery) = (
        inputs().join("qnucl_valid.fna"),
        inputs().join("qnx_late.fna"),
    );
    let cases: [(Program, &Path, &Path, &[&str]); 5] = [
        (Program::Blastp, &frozen().join("qprot.faa"), &sprot, &[]),
        (Program::Blastn, &nquery, &snucl, &[]),
        (Program::Blastx, &xquery, &sprot, &[]),
        (Program::Tblastn, &pquery, &snucl, &[]),
        (
            Program::Tblastx,
            &nquery,
            &snucl,
            &["-evalue", "1e-20", "-max_hsps", "2"],
        ),
    ];
    for (program, query, subject, extra) in cases {
        for (format, ext) in [("13", "json"), ("14", "xml")] {
            let what = format!("{} {format}", program.name());
            let mut argv = vec![
                program.name(),
                "-query",
                query.to_str().unwrap(),
                "-subject",
                subject.to_str().unwrap(),
            ];
            argv.extend_from_slice(extra);
            let reference_dir = work.join(format!("{}_{format}_reference", program.name()));
            fs::create_dir_all(&reference_dir).unwrap();
            let (reference_diagnostics, reference_master, reference_sides) =
                reference(program, query, subject, extra, format, &reference_dir, ext);
            assert!(
                reference_sides.len() >= 2,
                "{what}: {} side files",
                reference_sides.len()
            );

            let received = Arc::new(Received::default());
            let mut master = Vec::new();
            let mut diagnostics = Vec::new();
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: format,
                    sink: OutputSink::Files(&mut master, SideFileSink::new("r", received.clone())),
                }],
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            let mut parse_argv = vec!["LOSAT"];
            parse_argv.extend_from_slice(&argv);
            parse_argv.extend_from_slice(&["-num_threads", "1"]);
            let cli: Cli = try_parse_from(parse_argv).unwrap();
            let result = match cli.command {
                Commands::Blastp(args) => run_local_blastp(
                    args,
                    &records(query, ReaderConfig::query("BLASTP", true, true)),
                    &records(subject, ReaderConfig::subject("BLASTP", true, true)),
                    "",
                    "",
                    &mut outputs,
                ),
                Commands::Blastn(args) => run_local_blastn(
                    args,
                    &records(query, ReaderConfig::query("BLASTN", false, true)),
                    &records(subject, ReaderConfig::subject("BLASTN", false, true)),
                    &mut outputs,
                ),
                Commands::Blastx(args) => run_local_blastx(
                    args,
                    &records(query, ReaderConfig::query("BLASTX", false, true)),
                    &records(subject, ReaderConfig::subject("BLASTX", true, true)),
                    &mut outputs,
                ),
                Commands::Tblastn(args) => run_local_tblastn(
                    args,
                    &records(query, ReaderConfig::query("TBLASTN", true, true)),
                    &records(subject, ReaderConfig::subject("TBLASTN", false, true)),
                    &mut outputs,
                ),
                Commands::Tblastx(args) => run_local_tblastx(
                    args,
                    &records(query, ReaderConfig::query("TBLASTX", false, true)),
                    &records(subject, ReaderConfig::subject("TBLASTX", false, true)),
                    &mut outputs,
                ),
                _ => unreachable!(),
            };
            drop(outputs);
            result.unwrap_or_else(|error| panic!("{what}: {error:#}"));
            assert!(master == reference_master, "{what}: master");
            assert!(
                *received.0.lock().unwrap() == reference_sides,
                "{what}: side files"
            );
            assert_eq!(
                String::from_utf8_lossy(&diagnostics),
                String::from_utf8_lossy(&reference_diagnostics),
                "{what}: diagnostics"
            );
        }
    }
    fs::remove_dir_all(work).ok();
}
