#![allow(warnings, clippy::all)]

//! The shared BLASTN entry `run_local` writes several output formats from one search.
//! Every format must be byte-identical to the CLI run that requests only that format,
//! the warnings must equal that run's standard error, and the observer must report the
//! exact row and pairwise section of every HSP. (Identity of outfmt 0 with NCBI is
//! checked by docs/evidence/losat_web_e2a/check_losat.py over the frozen fixtures.) (Identity with the previous release over the regression
//! manifests is checked by docs/evidence/losat_web_e1a/capture_outputs.py.)

mod run_local_support;

use std::process::Command;

use run_local_support::{
    assert_observer_ranges, fixture_sequence, reader_records, run_formats, Run, TempFasta,
};
use LOSAT::algorithm::blastn::BlastnArgs;
use LOSAT::api::local_blast::{run_local_blastn, FormatOutput, OutputSink, ReportOutputs};
use LOSAT::blastinput::fasta_reader::ReaderConfig;
use LOSAT::cli::{try_parse_from, Cli, Commands};

const FORMATS: [&str; 3] = ["0", "6", "7"];

/// Three queries against four subjects cut from two genomes. The first query overlaps
/// two subjects, so `-max_target_seqs 1` removes a subject from its hit list.
struct Inputs {
    query: TempFasta,
    subject: TempFasta,
}

impl Inputs {
    fn new() -> Self {
        let a = fixture_sequence("LC738884.fasta");
        let b = fixture_sequence("LC741431.fasta");
        Self {
            query: TempFasta::new(
                "blastn_query.fna",
                &[
                    ("q_a0", &a[0..20_000]),
                    ("q_a100k", &a[100_000..120_000]),
                    ("q_b50k", &b[50_000..60_000]),
                ],
            ),
            subject: TempFasta::new(
                "blastn_subject.fna",
                &[
                    ("s_a10k", &a[10_000..40_000]),
                    ("s_a90k", &a[90_000..130_000]),
                    ("s_a5k", &a[5_000..25_000]),
                    ("s_b45k", &b[45_000..70_000]),
                ],
            ),
        }
    }

    fn args(&self, extra: &[&str]) -> BlastnArgs {
        let query = self.query.0.display().to_string();
        let subject = self.subject.0.display().to_string();
        let mut argv = vec![
            "LOSAT", "blastn", "-query", &query, "-subject", &subject, "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv).expect("blastn argv");
        match cli.command {
            Commands::Blastn(args) => args,
            _ => unreachable!("blastn subcommand"),
        }
    }

    fn run(&self, formats: &[&str], extra: &[&str], observe: bool) -> Run {
        let queries = reader_records(&self.query.0, ReaderConfig::query("BLASTN", false, true));
        let subjects = reader_records(
            &self.subject.0,
            ReaderConfig::subject("BLASTN", false, true),
        );
        run_formats(formats, observe, |outputs| {
            run_local_blastn(self.args(extra), &queries, &subjects, outputs)
        })
    }

    /// The CLI output file and standard error of one single-format run.
    fn cli(&self, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let out = self.query.0.with_extension(format!("cli{outfmt}.out"));
        let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
            .arg("blastn")
            .arg("-query")
            .arg(&self.query.0)
            .arg("-subject")
            .arg(&self.subject.0)
            .args(extra)
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-out")
            .arg(&out)
            .env_remove("LOSAT_TIMING")
            .output()
            .expect("run LOSAT CLI");
        assert!(output.status.success(), "CLI failed: {output:?}");
        let bytes = std::fs::read(&out).expect("CLI output");
        let _ = std::fs::remove_file(&out);
        (bytes, output.stderr)
    }
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
// ```c
// if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
//    hitlist_size = args[kArgMaxTargetSequences].AsInteger();
// }
// ```
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let inputs = Inputs::new();
    let mut row_counts = Vec::new();
    for extra in [
        &[][..],
        &["-max_target_seqs", "1"][..],
        &["-task", "blastn"][..],
        &["-num_threads", "2"][..],
    ] {
        let together = inputs.run(&FORMATS, extra, true);
        assert_eq!(
            together.hits.len(),
            together.ranges[1].len(),
            "one hit record per outfmt 6 row {extra:?}"
        );
        row_counts.push(together.ranges[1].len());
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = inputs.run(&[outfmt], extra, false);
            assert!(
                together.outputs[index] == alone.outputs[0],
                "run_local outfmt {outfmt} {extra:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, extra);
            assert!(
                together.outputs[index] == cli_output,
                "CLI outfmt {outfmt} {extra:?}"
            );
            assert_eq!(
                together.diagnostics, cli_stderr,
                "warnings of two formats equal those of the single-format CLI run {extra:?}"
            );
        }
    }
    assert!(row_counts[0] > 0, "the fixture must produce hits");
    assert!(
        row_counts[1] < row_counts[0],
        "-max_target_seqs 1 must drop a subject of this fixture: {row_counts:?}"
    );
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//     // Add tab in front of field, except for the first field.
//     if (iter != m_FieldsToShow.begin())
//         m_Ostream << m_FieldDelimiter;
//     x_PrintField(*iter);
// }
// m_Ostream << "\n";
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1970-1973
// ```c++
// x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
// ```
#[test]
fn observer_ranges_are_exact_rows_and_sections_of_the_same_hsp() {
    let inputs = Inputs::new();
    for extra in [
        &[][..],
        &["-max_target_seqs", "1"][..],
        &["-task", "blastn"][..],
    ] {
        let result = inputs.run(&FORMATS, extra, true);
        assert_observer_ranges(&result, [0, 1, 2], &format!("blastn {extra:?}"));
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
// ```c
// if (args[kArgOutputFormat]) {
//     string fmt_choice =
//         NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
// ```
#[test]
fn unsupported_formats_fail_before_searching() {
    let inputs = Inputs::new();
    let queries = reader_records(&inputs.query.0, ReaderConfig::query("BLASTN", false, true));
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true),
    );
    // NCBI ignores a custom specification with outfmt 0 (blast_args.cpp:2845-2851).
    for outfmt in ["5", "6 qseqid", "abc"] {
        let (mut valid, mut invalid, mut diagnostics) = (Vec::new(), Vec::new(), Vec::new());
        let mut outputs = ReportOutputs {
            formats: vec![
                FormatOutput {
                    outfmt: "6",
                    sink: OutputSink::Writer(&mut valid),
                },
                FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(&mut invalid),
                },
            ],
            diagnostics: &mut diagnostics,
            hits: None,
            observer: None,
        };
        let result = run_local_blastn(inputs.args(&[]), &queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty() && diagnostics.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:261-274,299-300
// ```c
//         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
//     ...
//                 if (!batch_size)
//                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
// ```
// LOSAT searches NCBI's query batches, so an invalid query after the first batch is
// reported with the batch that holds it: every format succeeds, as the CLI run does (the
// bytes against NCBI are checked by docs/evidence/losat_web_e2f/check_inputs.py).
#[test]
fn an_invalid_query_after_the_first_batch_is_reported_with_its_batch() {
    let genome = fixture_sequence("LC738884.fasta");
    for n_residues in [40, 150] {
        let all_n = vec![b'N'; n_residues];
        let query = TempFasta::new(
            "blastn_batches.fna",
            &[
                ("first_batch", &genome[0..6_000]),
                ("all_n", &all_n),
                ("valid", &genome[10_000..12_000]),
            ],
        );
        let subject = TempFasta::new("blastn_batches_subject.fna", &[("s", &genome[0..20_000])]);
        let inputs = Inputs { query, subject };
        let queries = reader_records(&inputs.query.0, ReaderConfig::query("BLASTN", false, true));
        let subjects = reader_records(
            &inputs.subject.0,
            ReaderConfig::subject("BLASTN", false, true),
        );
        for outfmt in FORMATS {
            let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(&mut output),
                }],
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            let result = run_local_blastn(inputs.args(&[]), &queries, &subjects, &mut outputs);
            drop(outputs);
            assert!(
                result.is_ok(),
                "{n_residues} N, outfmt {outfmt}: {result:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, &[]);
            assert_eq!(output, cli_output, "{n_residues} N, outfmt {outfmt}");
            assert_eq!(diagnostics, cli_stderr, "{n_residues} N, outfmt {outfmt}");
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:672-674
// ```c
//     if (hsp->score == 0) {
//         return CRef<CSeq_align>();
//     }
// ```
// The preliminary search reads the N of the subject as random bases, so seeds of the query
// land in them; the traceback scores those HSPs 0, and no format or hit record shows them.
#[test]
fn hsps_of_score_zero_are_not_reported() {
    let inputs = Inputs {
        query: TempFasta::new("blastn_score0_query.fna", &[("q11", b"GTCTGTCACAA")]),
        subject: TempFasta::new("blastn_score0_subject.fna", &[("sN", &[b'N'; 100])]),
    };
    let run = inputs.run(
        &["6", "7"],
        &["-task", "blastn", "-word_size", "4", "-evalue", "1e6"],
        true,
    );
    assert!(run.outputs[0].is_empty(), "{:?}", run.outputs[0]);
    assert!(String::from_utf8_lossy(&run.outputs[1]).contains("# 0 hits found"));
    assert!(run.hits.is_empty());
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:977,989
// ```c
//     sv.SetCoding(CSeq_data::e_Ncbi4na);
//     ...
//     sv.GetStrandData(strand, buffer);
// ```
// Without -lcase_masking, the letter case of the query changes nothing: both strands come
// from the same case-free sequence data.
#[test]
fn a_lowercase_query_finds_the_same_hits_on_both_strands() {
    let genome = fixture_sequence("LC738884.fasta");
    let segment = &genome[20_000..22_000];
    let reverse: Vec<u8> = segment
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect();
    let subject = TempFasta::new(
        "blastn_case_subject.fna",
        &[("plus", segment), ("minus", &reverse)],
    );
    let lower: Vec<u8> = segment.to_ascii_lowercase();
    let mut outputs = Vec::new();
    for (name, sequence) in [("upper", segment), ("lower", &lower[..])] {
        let query = TempFasta::new(&format!("blastn_case_{name}.fna"), &[("q", sequence)]);
        // The subject file belongs to `subject`; the copy of its path in `inputs` is
        // forgotten below so that it does not delete the file.
        let inputs = Inputs {
            query,
            subject: TempFasta(subject.0.clone()),
        };
        for task in ["megablast", "blastn"] {
            let run = inputs.run(&["6"], &["-task", task], false);
            outputs.push((name, task, run.outputs[0].clone()));
        }
        std::mem::forget(inputs.subject);
    }
    for task in ["megablast", "blastn"] {
        let rows: Vec<&Vec<u8>> = outputs
            .iter()
            .filter(|(_, t, _)| *t == task)
            .map(|(_, _, rows)| rows)
            .collect();
        assert_eq!(rows[0], rows[1], "{task}: lowercase and uppercase queries");
        let text = String::from_utf8_lossy(rows[0]);
        assert!(
            text.lines().any(|row| row.contains("\tminus\t")),
            "{task}: the fixture has a minus-strand hit"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:218-225
// ```c++
//     ITERATE(vector< CConstRef<CSeq_id> >, itr, original_seqids) {
//         CRef<CSeq_id> next_seqid(new CSeq_id());
//         string id_token = NcbiEmptyString;
//
//         if (((*itr)->IsGeneral() &&
//             (*itr)->AsFastaString().find("gnl|BL_ORD_ID")
//             != string::npos) ||
// 		(*itr)->AsFastaString().find("lcl|Subject_") != string::npos) {
// ```
// The tabular subject ID of a title whose first word starts with `Subject_` is the first
// word of the decoded title (oracle RP v01: `Subject_&amp;x desc` -> `Subject_&x`; v09:
// `Subject_4.` -> `Subject_4`); the query ID stays the title's first word. A title that NCBI
// cannot decode (v06, `Subject_\x81`) stops NCBI in the tabular formats (exit 255), which
// LOSAT rejects; its outfmt 0 report names it `Unknown`.
#[test]
fn subject_ids_of_subject_titles_follow_get_seq_id_list() {
    let genome = fixture_sequence("LC738884.fasta");
    let segment = &genome[1_000..1_330];
    let mut subject_sequence = genome[5_000..5_060].to_vec();
    subject_sequence.extend_from_slice(segment);
    let write = |name: &str, title: &[u8], sequence: &[u8]| {
        let path = std::env::temp_dir().join(format!(
            "losat_run_local_subject_ids_{}_{name}",
            std::process::id()
        ));
        let mut text = b">".to_vec();
        text.extend_from_slice(title);
        text.push(b'\n');
        for line in sequence.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
        std::fs::write(&path, text).expect("write temporary FASTA");
        TempFasta(path)
    };
    for (title, query_id, subject_id) in [
        (
            &b"Subject_&amp;x desc"[..],
            &b"Subject_&amp;x"[..],
            &b"Subject_&x"[..],
        ),
        (b"Subject_4.", b"Subject_4.", b"Subject_4"),
        (b"s1 Subject_&amp;x", b"s1", b"s1"),
    ] {
        let inputs = Inputs {
            query: write("q.fna", title, segment),
            subject: write("s.fna", title, &subject_sequence),
        };
        let run = inputs.run(&["6"], &[], false);
        let row = run.outputs[0].split(|&byte| byte == b'\n').next().unwrap();
        let fields: Vec<&[u8]> = row.split(|&byte| byte == b'\t').collect();
        assert_eq!(fields[0], query_id, "{:?}", String::from_utf8_lossy(title));
        assert_eq!(
            fields[1],
            subject_id,
            "{:?}",
            String::from_utf8_lossy(title)
        );
    }
    let unknown = b"Subject_\x81 desc";
    let inputs = Inputs {
        query: write("q.fna", b"q1", segment),
        subject: write("s.fna", unknown, &subject_sequence),
    };
    let queries = reader_records(&inputs.query.0, ReaderConfig::query("BLASTN", false, true));
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true),
    );
    for outfmt in ["6", "7"] {
        let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs = ReportOutputs {
            formats: vec![FormatOutput {
                outfmt,
                sink: OutputSink::Writer(&mut output),
            }],
            diagnostics: &mut diagnostics,
            hits: None,
            observer: None,
        };
        let error = run_local_blastn(inputs.args(&[]), &queries, &subjects, &mut outputs)
            .expect_err("a Subject_ title that NCBI cannot decode is rejected");
        drop(outputs);
        let message = format!("{error:#}");
        assert!(
            message.contains("non-UTF-8") && message.contains("not supported by LOSAT's BLASTN"),
            "{message}"
        );
        assert!(output.is_empty(), "outfmt {outfmt}: nothing is written");
    }
    let run = inputs.run(&["0"], &[], false);
    let report = String::from_utf8_lossy(&run.outputs[0]);
    assert!(report.contains("\nUnknown "), "{report}");
    assert!(
        report.contains(
            "Sequence with id Subject_1 no longer exists in database...alignment skipped"
        ),
        "{report}"
    );
}

/// Runs `run_local` with one writer per format and returns the outputs, the diagnostics
/// and the result.
fn run_local_outputs(
    inputs: &Inputs,
    formats: &[&str],
    extra: &[&str],
) -> (Vec<Vec<u8>>, Vec<u8>, anyhow::Result<()>) {
    let queries = reader_records(&inputs.query.0, ReaderConfig::query("BLASTN", false, true));
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true),
    );
    let mut sinks: Vec<Vec<u8>> = vec![Vec::new(); formats.len()];
    let mut diagnostics = Vec::new();
    let result = {
        let mut outputs = ReportOutputs {
            formats: formats
                .iter()
                .zip(sinks.iter_mut())
                .map(|(outfmt, sink)| FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(sink),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            hits: None,
            observer: None,
        };
        run_local_blastn(inputs.args(extra), &queries, &subjects, &mut outputs)
    };
    (sinks, diagnostics, result)
}

/// The exit code and the message of NCBI's error that a run ended with.
fn native_error(result: anyhow::Result<()>) -> (i32, String) {
    let error = result.expect_err("the run fails");
    let native = error
        .downcast_ref::<LOSAT::cli::NativeError>()
        .unwrap_or_else(|| panic!("NCBI's error: {error:#}"));
    (native.exit, native.message.clone())
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-788
// ```c
//         catch(CBlastException & e ) {
//         	// Skip bad subject sequence
//         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
//         		seqblk_vec->push_back(subj);
//         ...
//         		warning += "Subject sequence contains no data";
//         		ERR_POST(Warning << warning);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-636
// ```c
//         } catch (const CException& e) {
//             // FIXME: is index this the right value for the 2nd arg? Also, how
//             // to determine whether the message should contain a warning or
//             // error?
//             CRef<CSearchMessage> m
// ```
// A record without residues stays in its role (oracle BI e01, e02, e06, e10): a subject
// gets its warning before the search, counts in the database statistics and never hits; a
// query gets its warning before its report, which has no hits, and is counted.
#[test]
fn records_without_residues_are_kept_with_ncbis_warnings() {
    let genome = fixture_sequence("LC738884.fasta");
    let inputs = Inputs {
        query: TempFasta::new(
            "blastn_empty_records_query.fna",
            &[
                ("q1", &genome[0..2_000]),
                ("q2 empty two", b""),
                ("q3", &genome[10_000..12_000]),
            ],
        ),
        subject: TempFasta::new(
            "blastn_empty_records_subject.fna",
            &[
                ("s1", &genome[0..5_000]),
                ("s2 no letters", b""),
                ("s3", &genome[9_000..13_000]),
            ],
        ),
    };
    let run = inputs.run(&FORMATS, &[], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [blastn] Subject_2 s2 no letters: Subject sequence contains no data\n\
         Warning: [blastn] Query_2 q2 empty two: Sequence contains no data \n"
    );
    let pairwise = String::from_utf8_lossy(&run.outputs[0]);
    assert!(
        pairwise.contains("3 sequences; 9,000 total letters"),
        "{pairwise}"
    );
    assert!(
        pairwise.contains("Query= q2 empty two\n\nLength=0\n\n\n***** No hits found *****\n"),
        "{pairwise}"
    );
    let tabular = String::from_utf8_lossy(&run.outputs[1]);
    assert!(
        tabular.lines().any(|row| row.starts_with("q1\ts1\t")),
        "{tabular}"
    );
    assert!(
        tabular.lines().any(|row| row.starts_with("q3\ts3\t")),
        "{tabular}"
    );
    assert!(
        !tabular.contains("q2") && !tabular.contains("s2"),
        "{tabular}"
    );
    let commented = String::from_utf8_lossy(&run.outputs[2]);
    assert!(
        commented.contains("# Query: q2 empty two\n# Database: ")
            && commented.contains("# BLAST processed 3 queries\n"),
        "{commented}"
    );
    for (index, outfmt) in FORMATS.into_iter().enumerate() {
        let (cli_output, cli_stderr) = inputs.cli(outfmt, &[]);
        assert_eq!(run.outputs[index], cli_output, "outfmt {outfmt}");
        assert_eq!(run.diagnostics, cli_stderr, "outfmt {outfmt}");
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:649-652
// ```c
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// A batch of queries without residues only stops NCBI after the prolog with one message
// per query (oracle BI e03).
#[test]
fn a_batch_of_queries_without_residues_stops_after_the_prolog() {
    let genome = fixture_sequence("LC738884.fasta");
    let inputs = Inputs {
        query: TempFasta::new(
            "blastn_all_empty_query.fna",
            &[("q1", b""), ("q2 two", b"")],
        ),
        subject: TempFasta::new("blastn_all_empty_subject.fna", &[("s1", &genome[0..5_000])]),
    };
    let (outputs, diagnostics, result) = run_local_outputs(&inputs, &["0", "6"], &[]);
    assert_eq!(
        native_error(result),
        (
            3,
            "BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n"
                .to_string()
        )
    );
    assert!(diagnostics.is_empty(), "{diagnostics:?}");
    let prolog = String::from_utf8_lossy(&outputs[0]);
    assert!(
        prolog.starts_with("BLASTN 2.17.0+") && !prolog.contains("Query="),
        "{prolog}"
    );
    assert!(outputs[1].is_empty());
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-980
// ```c
//    if (sbp->gbp) {
//        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
//        if (Blast_SubjectIsTranslated(program_number)) {
//            min_subject_length/=3;
//        }
//    } else {
//        min_subject_length = (Int4) (total_length/num_seqs);
//    }
//
//    if(min_subject_length <=0) {
// 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
//    }
// ```
// Subjects with fewer letters than subjects (every subject without residues, oracle BI e08)
// stop NCBI when it searches the first batch, after the subjects' warnings and the prolog.
#[test]
fn subjects_with_fewer_letters_than_subjects_stop_the_search() {
    let genome = fixture_sequence("LC738884.fasta");
    for (subjects, warnings) in [
        (
            vec![("s1", &b"A"[..]), ("s2", &b""[..])],
            "Warning: [blastn] Subject_2 s2: Subject sequence contains no data\n",
        ),
        (
            vec![("s1 one", &b""[..]), ("s2", &b""[..])],
            "Warning: [blastn] Subject_1 s1 one: Subject sequence contains no data\n\
             Warning: [blastn] Subject_2 s2: Subject sequence contains no data\n",
        ),
    ] {
        let inputs = Inputs {
            query: TempFasta::new(
                "blastn_short_subjects_query.fna",
                &[("q1", &genome[0..500])],
            ),
            subject: TempFasta::new("blastn_short_subjects.fna", &subjects),
        };
        let (outputs, diagnostics, result) = run_local_outputs(&inputs, &["0", "7"], &[]);
        assert_eq!(
            native_error(result),
            (
                3,
                "BLAST engine error: The average subject length is too short\n".to_string()
            )
        );
        assert_eq!(String::from_utf8_lossy(&diagnostics), warnings);
        assert!(String::from_utf8_lossy(&outputs[0]).starts_with("BLASTN 2.17.0+"));
        assert!(outputs[1].is_empty());
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-460
// ```c++
//     // set sequence range
//     retval->SetInt().SetFrom(from);
//     retval->SetInt().SetTo((to > 0 && to < seqlen) ? to : (seqlen-1));
// ```
// A range that starts just past a record's end gives an interval without letters, which
// NCBI sets up as a record without residues (oracle BI r3, r4).
#[test]
fn intervals_without_letters_are_records_without_residues() {
    let genome = fixture_sequence("LC738884.fasta");
    let inputs = Inputs {
        query: TempFasta::new(
            "blastn_empty_interval_query.fna",
            &[("q1", &genome[0..2_000]), ("q2 short", &genome[0..300])],
        ),
        subject: TempFasta::new(
            "blastn_empty_interval_subject.fna",
            &[("s1", &genome[0..5_000]), ("s2 short", &genome[0..300])],
        ),
    };
    let run = inputs.run(&["6"], &["-subject_loc", "301-1000"], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [blastn] Subject_2 s2 short: Subject sequence contains no data\n"
    );
    let rows = String::from_utf8_lossy(&run.outputs[0]);
    assert!(rows.lines().all(|row| row.contains("\ts1\t")), "{rows}");
    assert!(!rows.is_empty());
    let run = inputs.run(&["7"], &["-query_loc", "301-1000"], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [blastn] Query_2 q2 short: Sequence contains no data \n"
    );
    let report = String::from_utf8_lossy(&run.outputs[0]);
    assert!(
        report.contains("# Query: q2 short\n") && report.contains("# BLAST processed 2 queries\n"),
        "{report}"
    );
}
