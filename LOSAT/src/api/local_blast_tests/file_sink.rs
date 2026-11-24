//! No NCBI counterpart (the library's own sinks): a library `File` sink (`OutputSink::File`) gets
//! the bytes that a writer of the caller gets for the same call, in every program. Its file is
//! kept for the run (`api::local_blast::KeptFile`), so the outfmt 0 to 4 prolog that each program
//! writes before the search and the reports after it go to the one file: a call that succeeds and
//! one that stops at its first query batch write what they wrote before the file was kept (the
//! prolog with the report), and a call that stops after the prolog leaves the prolog in the file,
//! as on the command line.

use super::support::{fixture_sequence, reader_records, temp_path, TempFasta};
use crate::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_blastx, run_local_tblastn, run_local_tblastx,
    OutputSink, ReportOutputs,
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
    // A refusal where NCBI starts the search: the file holds the prolog, as the writer.
    let (bytes, error) = same_as_writer(
        "blastp -use_sw_tback -comp_based_stats 0",
        writer_and_file("0", |outputs| {
            run(
                "blastp",
                (&proteins, &proteins),
                &["-use_sw_tback", "-comp_based_stats", "0"],
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
