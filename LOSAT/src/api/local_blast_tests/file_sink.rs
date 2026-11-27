//! No NCBI counterpart (the library's own sinks): a library `File` sink (`OutputSink::File`) gets
//! the bytes that a writer of the caller gets for the same call, in every program. Its file is
//! kept for the run (`api::local_blast::KeptFile`), so the outfmt 0 to 4 prolog that each program
//! writes before the search and the reports after it go to the one file: a call that succeeds and
//! one that stops at its first query batch write what they wrote before the file was kept (the
//! prolog with the report), and a call that stops after the prolog leaves the prolog in the file,
//! as on the command line. psiblast's rounds and a database search write one file too; sinks that
//! name one file share it; a file that cannot be created fails the call where it is first written.

use std::io::Write;

use super::support::{fixture_sequence, reader_records, temp_path, TempFasta};
use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_database_blastp,
    run_local_tblastn, run_local_tblastx, FormatOutput, FormatRefusal, KeptFile, OutputSink,
    ReportOutputs,
};
use crate::api::QueryEnd;
use crate::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
use crate::cli::{try_parse_from, Cli, Commands};

/// What one call wrote: its sink's bytes (`None`: no file was made), its error's text and its
/// diagnostics.
type Written = (Option<Vec<u8>>, Option<String>, Vec<u8>);

/// `search` with the one format `outfmt`, to a writer of the caller and to a file.
fn writer_and_file(
    outfmt: &str,
    search: impl Fn(&mut ReportOutputs<'_>) -> anyhow::Result<()>,
) -> (Written, Written) {
    let mut bytes = Vec::new();
    let mut diagnostics = Vec::new();
    let result = {
        let mut outputs =
            ReportOutputs::single(outfmt, OutputSink::Writer(&mut bytes), &mut diagnostics);
        search(&mut outputs)
    };
    let writer = (
        Some(bytes),
        result.err().map(|e| format!("{e:#}")),
        diagnostics,
    );
    let path = temp_path("file_sink.out");
    let mut diagnostics = Vec::new();
    let result = {
        let mut outputs = ReportOutputs::single(outfmt, OutputSink::File(&path), &mut diagnostics);
        search(&mut outputs)
    };
    let file = std::fs::read(&path).ok();
    let _ = std::fs::remove_file(&path);
    (
        writer,
        (file, result.err().map(|e| format!("{e:#}")), diagnostics),
    )
}

/// Checks that the file got the writer's bytes, error and diagnostics; returns the writer's.
fn same_as_writer(context: &str, (writer, file): (Written, Written)) -> (Vec<u8>, Option<String>) {
    assert_eq!(file.1, writer.1, "{context}: the error");
    assert!(file.2 == writer.2, "{context}: the diagnostics");
    let written = writer.0.expect("the writer's bytes");
    match &file.0 {
        Some(bytes) => assert!(
            *bytes == written,
            "{context}: the file has {} bytes, the writer {}",
            bytes.len(),
            written.len()
        ),
        None => assert!(
            written.is_empty(),
            "{context}: no file, {} bytes to the writer",
            written.len()
        ),
    }
    (written, writer.1)
}

/// The parsed arguments of `losat PROGRAM -query QUERY -subject SUBJECT EXTRA`.
fn command(program: &str, query: &TempFasta, subject: &TempFasta, extra: &[&str]) -> Commands {
    let (query, subject) = (
        query.0.display().to_string(),
        subject.0.display().to_string(),
    );
    let mut argv = vec!["losat", program, "-query", &query, "-subject", &subject];
    argv.extend_from_slice(extra);
    let cli: Cli = try_parse_from(argv).expect("argv");
    cli.command
}

/// The records of `fasta` as NCBI's reader reads `program`'s queries or subjects.
fn records(
    fasta: &TempFasta,
    program: &'static str,
    query: bool,
    protein: bool,
) -> Vec<FastaRecord> {
    let config = if query {
        ReaderConfig::query(program, protein, true)
    } else {
        ReaderConfig::subject(program, protein, true)
    };
    reader_records(&fasta.0, config)
}

/// `losat PROGRAM` over the records of `query` and `subject` to `outputs`.
fn run(
    program: &str,
    (query, subject): (&TempFasta, &TempFasta),
    extra: &[&str],
    outputs: &mut ReportOutputs<'_>,
) -> anyhow::Result<()> {
    let (name, query_protein, subject_protein) = match program {
        "blastn" => ("BLASTN", false, false),
        "blastp" => ("BLASTP", true, true),
        "blastx" => ("BLASTX", false, true),
        "tblastn" => ("TBLASTN", true, false),
        _ => ("TBLASTX", false, false),
    };
    let queries = records(query, name, true, query_protein);
    let subjects = records(subject, name, false, subject_protein);
    let end = QueryEnd::Input;
    match command(program, query, subject, extra) {
        Commands::Blastn(args) => run_local_blastn(args, &queries, &end, &subjects, outputs),
        Commands::Blastp(args) => {
            run_local_blastp(args, &queries, &end, &subjects, "", "", outputs)
        }
        Commands::Blastx(args) => run_local_blastx(args, &queries, &end, &subjects, outputs),
        Commands::Tblastn(args) => run_local_tblastn(args, &queries, &end, &subjects, outputs),
        Commands::Tblastx(args) => run_local_tblastx(args, &queries, &end, &subjects, outputs),
        _ => unreachable!("a program of run_local"),
    }
}

#[test]
fn a_file_sink_gets_the_bytes_of_a_writer_in_every_program() {
    let protein = fixture_sequence("CoBV.faa");
    let genome = fixture_sequence("LC738874.fasta");
    let proteins = TempFasta::new("file_sink_p.faa", &[("p1", &protein)]);
    let nucleotides = TempFasta::new("file_sink_n.fna", &[("n1", &genome[..6_000])]);
    let genomes = TempFasta::new("file_sink_g.fna", &[("g1", &genome)]);
    // A query without letters: the first query batch stops after the prolog.
    let empty = TempFasta::new("file_sink_e.fa", &[("e1", b"")]);
    let programs: [(&str, &TempFasta, &TempFasta); 5] = [
        ("blastn", &nucleotides, &nucleotides),
        ("blastp", &proteins, &proteins),
        ("blastx", &nucleotides, &proteins),
        ("tblastn", &proteins, &genomes),
        ("tblastx", &nucleotides, &nucleotides),
    ];
    for (program, query, subject) in programs {
        // outfmt 0 and 1: the prolog before the search; 7: no prolog.
        for outfmt in ["0", "1", "7"] {
            if program == "tblastx" && outfmt == "1" {
                // TBLASTX writes no query-anchored report.
                continue;
            }
            let context = format!("{program} -outfmt {outfmt}");
            let (bytes, error) = same_as_writer(
                &context,
                writer_and_file(outfmt, |outputs| {
                    run(program, (query, subject), &[], outputs)
                }),
            );
            assert_eq!(error, None, "{context}");
            assert!(
                outfmt == "7" || bytes.windows(7).any(|window| window == b"Query= "),
                "{context}: a report"
            );
            let context = format!("{context} with a query without letters");
            same_as_writer(
                &context,
                writer_and_file(outfmt, |outputs| {
                    run(program, (&empty, subject), &[], outputs)
                }),
            );
        }
    }
    // A refusal where NCBI starts the search: the file holds the prolog, as the writer
    // (`-use_sw_tback -comp_based_stats 0 -ungapped`: `redo_run::sw_ungapped_check`).
    let (bytes, error) = same_as_writer(
        "blastp -use_sw_tback -comp_based_stats 0 -ungapped",
        writer_and_file("0", |outputs| {
            run(
                "blastp",
                (&proteins, &proteins),
                &["-use_sw_tback", "-comp_based_stats", "0", "-ungapped"],
                outputs,
            )
        }),
    );
    assert!(error.is_some(), "the refusal");
    assert!(
        bytes.starts_with(b"BLASTP"),
        "the prolog before the refusal"
    );
}

/// The outputs of `formats` (a format and its sink each), every one failing the call.
fn outputs<'a>(
    formats: Vec<(&'a str, OutputSink<'a>)>,
    diagnostics: &'a mut Vec<u8>,
) -> ReportOutputs<'a> {
    ReportOutputs {
        other_outfmts: &[],
        html_posts: None,
        formats: formats
            .into_iter()
            .map(|(outfmt, sink)| FormatOutput {
                outfmt: outfmt.as_bytes(),
                refusal: FormatRefusal::FailsCall,
                sink,
            })
            .collect(),
        diagnostics,
        later_warnings: None,
        hits: None,
        observer: None,
    }
}

// The file of a `File` sink is created at the run's first write or flush to it (BLASTP's block
// search writes it as a writer, `XBlockSink`), and the sinks of a run that name one file, also
// through another spelling of its folder, write it through one handle, in the order written.
#[test]
fn a_kept_file_is_made_at_its_first_write_and_shared_by_the_sinks_naming_it() {
    let dir = temp_path("file_sink_kept");
    std::fs::create_dir_all(dir.join("sub")).expect("the folders");
    let path = dir.join("out.tsv");
    let mut kept = vec![Some(KeptFile::of(&path, &[]))];
    let again = KeptFile::of(&dir.join("sub").join("..").join("out.tsv"), &kept);
    kept.push(Some(again));
    let other = KeptFile::of(&dir.join("sub").join("out.tsv"), &kept);
    kept.push(Some(other));
    {
        let [Some(first), Some(second), Some(third)] = &mut kept[..] else {
            unreachable!("three kept files")
        };
        assert!(!path.exists());
        first.flush().expect("flush");
        assert!(path.exists());
        first.write_all(b"a\t1\n").expect("write");
        second.write_all(b"b\t2\n").expect("write");
        first.write_all(b"c\t3\n").expect("write");
        third.write_all(b"d\n").expect("write");
    }
    drop(kept);
    let written = std::fs::read(&path).expect("the file");
    let other = std::fs::read(dir.join("sub").join("out.tsv")).expect("the other file");
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!(written, b"a\t1\nb\t2\nc\t3\n");
    assert_eq!(other, b"d\n");
}

// No NCBI counterpart (one `-out` per run): two outputs of one call that name one file write it
// through one handle (`KeptFile::of`), so the file holds what each writes, in the order written,
// with nothing lost or overwritten: with outfmt 0 and 6, outfmt 0's prolog (before the search),
// then the reports in the order of the formats (BLASTP writes each output's whole report in turn).
#[test]
fn two_outputs_naming_one_file_write_it_in_turn() {
    let protein = fixture_sequence("CoBV.faa");
    let proteins = TempFasta::new("file_sink_shared.faa", &[("p1", &protein)]);
    for formats in [["6", "0"], ["0", "6"]] {
        let (mut first, mut second, mut diagnostics) = (Vec::new(), Vec::new(), Vec::new());
        {
            let sinks = vec![
                (formats[0], OutputSink::Writer(&mut first)),
                (formats[1], OutputSink::Writer(&mut second)),
            ];
            let mut outputs = outputs(sinks, &mut diagnostics);
            run("blastp", (&proteins, &proteins), &[], &mut outputs).expect("the search");
        }
        let path = temp_path("file_sink_shared.out");
        let mut diagnostics = Vec::new();
        let result = {
            let sinks = vec![
                (formats[0], OutputSink::File(&path)),
                (formats[1], OutputSink::File(&path)),
            ];
            let mut outputs = outputs(sinks, &mut diagnostics);
            run("blastp", (&proteins, &proteins), &[], &mut outputs)
        };
        let written = std::fs::read(&path).expect("the file");
        let _ = std::fs::remove_file(&path);
        result.expect("the search to one file");
        let (report, rows) = match formats[0] {
            "0" => (&first, &second),
            _ => (&second, &first),
        };
        assert!(!rows.is_empty(), "the fixture must produce hits");
        // The prolog: the report's bytes before its first query.
        let prolog = report
            .windows(9)
            .position(|window| window == b"\n\nQuery= ")
            .expect("a query's report");
        let expected = match formats[0] {
            "0" => [report.as_slice(), rows].concat(),
            _ => [&report[..prolog], rows, &report[prolog..]].concat(),
        };
        assert!(
            written == expected,
            "{formats:?}: the file has {} bytes, the outputs {}",
            written.len(),
            expected.len()
        );
    }
}

// psiblast's rounds write one file (each round's report and the round lines after the first open
// of the file append to it; before the kept file, each `open` created the file again), and a
// database search writes its `File` sink as a writer (`run_local_database`, here BLASTP's with
// psiblast's arguments): the file gets the writer's bytes (the case of the rows
// `db.ps15b.pk.mts2.*`).
#[test]
fn psiblast_rounds_of_a_database_search_write_one_file() {
    let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let query = manifest.join("tests/fasta/psiblast/ps5_pk.faa");
    let database = manifest.join("tests/db_search/paja");
    let (query_name, database_name) = (query.display().to_string(), database.display().to_string());
    let records = reader_records(&query, ReaderConfig::query("PSIBLAST", true, true));
    let search = |outputs: &mut ReportOutputs<'_>| {
        let argv = [
            "losat",
            "psiblast",
            "-query",
            &query_name,
            "-db",
            &database_name,
            "-num_iterations",
            "3",
            "-inclusion_ethresh",
            "5",
            "-max_target_seqs",
            "2",
        ];
        let cli: Cli = try_parse_from(argv).expect("argv");
        let Commands::Psiblast(args) = cli.command else {
            unreachable!("psiblast's arguments")
        };
        run_local_database_blastp(
            args.into_blastp_args()?,
            &records,
            &QueryEnd::Input,
            crate::blastinput::ncbi_environment::ApplicationSettings::default(),
            outputs,
        )
    };
    for outfmt in ["0", "7"] {
        let context = format!("psiblast -db -outfmt {outfmt}");
        let (bytes, error) = same_as_writer(&context, writer_and_file(outfmt, search));
        assert_eq!(error, None, "{context}");
        let round = match outfmt {
            "0" => &b"Results from round 2"[..],
            _ => &b"# Iteration: 2"[..],
        };
        assert!(
            bytes.windows(round.len()).any(|window| window == round),
            "{context}: a second round"
        );
    }
}

// A `File` sink whose file cannot be created (its folder does not exist) fails the call where the
// file is first written: with outfmt 0 at the prolog, before the search (as NCBI's `-out`, opened
// with the arguments), so outfmt 7 to a writer before it gets nothing; with outfmt 6 at its report
// after the search, as before the kept file, after outfmt 7's report (BLASTX opens every sink
// before the search, as before).
#[test]
fn a_file_that_cannot_be_made_fails_the_call_where_it_is_first_written() {
    let protein = fixture_sequence("CoBV.faa");
    let genome = fixture_sequence("LC738874.fasta");
    let proteins = TempFasta::new("file_sink_missing_p.faa", &[("p1", &protein)]);
    let nucleotides = TempFasta::new("file_sink_missing_n.fna", &[("n1", &genome[..6_000])]);
    let genomes = TempFasta::new("file_sink_missing_g.fna", &[("g1", &genome)]);
    let missing = temp_path("file_sink_no_such_folder").join("out.txt");
    let programs: [(&str, &TempFasta, &TempFasta); 5] = [
        ("blastn", &nucleotides, &nucleotides),
        ("blastp", &proteins, &proteins),
        ("blastx", &nucleotides, &proteins),
        ("tblastn", &proteins, &genomes),
        ("tblastx", &nucleotides, &nucleotides),
    ];
    for (program, query, subject) in programs {
        for outfmt in ["0", "6"] {
            let context = format!("{program}: -outfmt 7 to a writer, {outfmt} to {missing:?}");
            let (mut table, mut diagnostics) = (Vec::new(), Vec::new());
            let result = {
                let sinks = vec![
                    ("7", OutputSink::Writer(&mut table)),
                    (outfmt, OutputSink::File(&missing)),
                ];
                let mut outputs = outputs(sinks, &mut diagnostics);
                run(program, (query, subject), &[], &mut outputs)
            };
            assert!(result.is_err(), "{context}: the call fails");
            assert!(!missing.exists(), "{context}");
            let before_search = outfmt == "0" || program == "blastx";
            assert_eq!(
                table.is_empty(),
                before_search,
                "{context}: outfmt 7 has {} bytes",
                table.len()
            );
        }
    }
}
