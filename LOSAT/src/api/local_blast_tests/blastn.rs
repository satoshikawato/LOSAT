//! The shared BLASTN entry `run_local` in process: what the entry itself does, not what the
//! command line writes. The tests that compare the entry's reports with the `losat` binary go
//! through `losat::api` (`LOSAT/tests/run_local_blastn.rs`).

use super::support::{
    assert_observer_ranges, fixture_sequence, reader_records, run_formats, Run, TempFasta,
};
use crate::algorithm::blastn::BlastnArgs;
use crate::api::local_blast::{
    run_local_blastn, FormatObserver, FormatOutput, HspIndex, OutputSink, ReportOutputs,
};
use crate::api::{read_fasta, run_batches, validate, Inputs as ApiInputs, Program, Role};
use crate::blastinput::fasta_reader::ReaderConfig;
use crate::cli::{try_parse_from, Cli, Commands};

const FORMATS: [&str; 12] = [
    "0", "6", "7", "8", "9", "12", "5", "15", "16", "11", "17", "17 SQ SR",
];

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
            "losat", "blastn", "-query", &query, "-subject", &subject, "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv).expect("blastn argv");
        match cli.command {
            Commands::Blastn(args) => args,
            _ => unreachable!("blastn subcommand"),
        }
    }

    fn run(&self, formats: &[&str], extra: &[&str], observe: bool) -> Run {
        let queries = reader_records(
            &self.query.0,
            ReaderConfig::query("BLASTN", false, true).with_defline(true),
        );
        let subjects = reader_records(
            &self.subject.0,
            ReaderConfig::subject("BLASTN", false, true).with_defline(true),
        );
        run_formats(formats, observe, |outputs| {
            run_local_blastn(
                self.args(extra),
                &queries,
                &crate::api::QueryEnd::Input,
                &subjects,
                outputs,
            )
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1970-1973
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
#[test]
fn unsupported_formats_fail_before_searching() {
    let inputs = Inputs::new();
    let queries = reader_records(
        &inputs.query.0,
        ReaderConfig::query("BLASTN", false, true).with_defline(true),
    );
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true).with_defline(true),
    );
    // NCBI ignores a custom specification with outfmt 0 (c++/src/algo/blast/blastinput/blast_args.cpp:2845-2851); a delimiter
    // that a field uses between its entries is NCBI's error (blast_args.cpp:2753-2792). Outfmt 19
    // (AIRR) is igblastn's alone (BLASTN writes 1 to 4 since line RG T3b).
    for outfmt in ["19", "6 delim=; qseqid staxids", "abc"] {
        let (mut valid, mut invalid, mut diagnostics) = (Vec::new(), Vec::new(), Vec::new());
        let mut outputs = ReportOutputs {
            formats: vec![
                FormatOutput {
                    outfmt: "6",
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(&mut valid),
                },
                FormatOutput {
                    outfmt,
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(&mut invalid),
                },
            ],
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        let result = run_local_blastn(
            inputs.args(&[]),
            &queries,
            &crate::api::QueryEnd::Input,
            &subjects,
            &mut outputs,
        );
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty() && diagnostics.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}

// NCBI derives the hit list size per output format (c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973): with
// `-num_descriptions 3` it is 250 for outfmt 0 and 500 for outfmt 6. One run searches once,
// so formats that need different hit list sizes are refused before anything is written
// (`blastinput/report_counts.rs`), while formats of one hit list size still share a run.
#[test]
fn formats_needing_different_hit_list_sizes_are_refused() {
    let inputs = Inputs::new();
    let queries = reader_records(
        &inputs.query.0,
        ReaderConfig::query("BLASTN", false, true).with_defline(true),
    );
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true).with_defline(true),
    );
    let attempt = |formats: &[&str], extra: &[&str]| {
        let mut sinks: Vec<Vec<u8>> = formats.iter().map(|_| Vec::new()).collect();
        let mut diagnostics = Vec::new();
        let mut outputs = ReportOutputs {
            formats: formats
                .iter()
                .zip(sinks.iter_mut())
                .map(|(outfmt, sink)| FormatOutput {
                    outfmt,
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(sink),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        let result = run_local_blastn(
            inputs.args(extra),
            &queries,
            &crate::api::QueryEnd::Input,
            &subjects,
            &mut outputs,
        );
        drop(outputs);
        (result, sinks, diagnostics)
    };
    for (formats, extra) in [
        (&["0", "6"][..], &["-num_descriptions", "3"][..]),
        (&["6", "0"][..], &["-num_descriptions", "3"][..]),
        (&["0", "6"][..], &["-num_alignments", "3"][..]),
        (
            &["0", "6", "7"][..],
            &["-num_descriptions", "3", "-num_alignments", "2"][..],
        ),
    ] {
        let (result, sinks, _) = attempt(formats, extra);
        let error = format!("{:#}", result.expect_err("different hit list sizes"));
        assert!(
            error.contains("need different hit list sizes"),
            "{formats:?} {extra:?}: {error}"
        );
        // The option warnings of the first format may precede the refusal (BLASTP's
        // `check_options`, as the CLI writes a warning before an options error), so the
        // diagnostics are not checked; no report bytes are written.
        assert!(
            sinks.iter().all(Vec::is_empty),
            "{formats:?} {extra:?}: nothing is written"
        );
    }
    // The same formats without these options, or with one hit list size, run together.
    for extra in [&[][..], &["-max_target_seqs", "3"][..]] {
        let (result, sinks, _) = attempt(&["0", "6"], extra);
        result.unwrap_or_else(|error| panic!("{extra:?}: {error:#}"));
        assert!(
            !sinks[0].is_empty(),
            "{extra:?}: outfmt 0 writes its header"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:672-674
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:977,989
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
        let path = crate::utils::test_dir::temp_dir().join(format!(
            "losat_run_local_subject_ids_{}_{name}",
            crate::utils::test_dir::process_id()
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
    let queries = reader_records(
        &inputs.query.0,
        ReaderConfig::query("BLASTN", false, true).with_defline(true),
    );
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true).with_defline(true),
    );
    for outfmt in ["6", "7"] {
        let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs = ReportOutputs {
            formats: vec![FormatOutput {
                outfmt,
                refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                sink: OutputSink::Writer(&mut output),
            }],
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        let error = run_local_blastn(
            inputs.args(&[]),
            &queries,
            &crate::api::QueryEnd::Input,
            &subjects,
            &mut outputs,
        )
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
    let queries = reader_records(
        &inputs.query.0,
        ReaderConfig::query("BLASTN", false, true).with_defline(true),
    );
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true).with_defline(true),
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
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(sink),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        run_local_blastn(
            inputs.args(extra),
            &queries,
            &crate::api::QueryEnd::Input,
            &subjects,
            &mut outputs,
        )
    };
    (sinks, diagnostics, result)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1312-1326
// The command line's `-word_size` is at least 4 (clap's `blastn_word_size`); a caller of
// `run_local` that sets the field itself gets the error of NCBI's core (the task's lookup word
// is longer, so the search would have wrapped `word - lut + 1` and reported other hits).
#[test]
fn a_word_size_below_four_is_ncbis_core_error() {
    let inputs = Inputs::new();
    for (task, word, message) in [
        ("megablast", 0, "Word-size must be greater than zero"),
        (
            "blastn",
            3,
            "Word-size must be 4 or greater for nucleotide comparison",
        ),
        (
            "blastn-short",
            1,
            "Word-size must be 4 or greater for nucleotide comparison",
        ),
    ] {
        let queries = reader_records(
            &inputs.query.0,
            ReaderConfig::query("BLASTN", false, true).with_defline(true),
        );
        let subjects = reader_records(
            &inputs.subject.0,
            ReaderConfig::subject("BLASTN", false, true).with_defline(true),
        );
        let mut args = inputs.args(&["-task", task]);
        args.word_size = Some(word);
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        let result = {
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: "6",
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(&mut sink),
                }],
                diagnostics: &mut diagnostics,
                later_warnings: None,
                hits: None,
                observer: None,
            };
            run_local_blastn(
                args,
                &queries,
                &crate::api::QueryEnd::Input,
                &subjects,
                &mut outputs,
            )
        };
        assert_eq!(
            native_error(result),
            (
                1,
                format!("LOSAT query/options error: {message}\nPlease refer to the LOSAT user manual.\n")
            ),
            "-task {task} word size {word}"
        );
        assert!(sink.is_empty() && diagnostics.is_empty());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1263-1269 and :1422-1427
// (cited at `check_scoring_options_with`); c++/src/algo/blast/api/blast_aux.cpp:662-673
// (`ThrowIfInvalidTask`, cited at `TaskParser::check`). Values that the command line cannot give
// (a template length the parser refuses, `limit_lookup`, which has no spelling, a task the parser
// refuses) get NCBI's errors, or LOSAT's for another program's task, instead of a search.
#[test]
fn hand_set_template_lookup_and_task_are_refused() {
    let inputs = Inputs::new();
    let queries = reader_records(
        &inputs.query.0,
        ReaderConfig::query("BLASTN", false, true).with_defline(true),
    );
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true).with_defline(true),
    );
    let cases: [(&[&str], &str, fn(&mut BlastnArgs)); 5] = [
        (
            &["-task", "dc-megablast"],
            "Invalid discontiguous template parameters: template length must be 16, 18, or 21",
            |a| a.template_length = Some(17),
        ),
        // NCBI's options factory reads the task in lower case (c++/src/algo/blast/api/blast_options_handle.cpp:341):
        // DC-MEGABLAST is dc-megablast, whose template the next check refuses.
        (
            &["-task", "dc-megablast"],
            "Invalid discontiguous template parameters: template length must be 16, 18, or 21",
            |a| {
                a.task = "DC-MEGABLAST".into();
                a.template_length = Some(17)
            },
        ),
        (
            &["-task", "blastn"],
            "The limit_lookup option can only be used with word size >= 16",
            |a| a.limit_lookup = true,
        ),
        (&[], "'blastz' is not a supported task", |a| a.task = "BlastZ".into()),
        (
            &[],
            "the task tblastn is not supported by LOSAT's BLASTN (use megablast, dc-megablast, blastn, blastn-short or rmblastn)",
            |a| a.task = "tblastn".into(),
        ),
    ];
    for (extra, message, set) in cases {
        let mut args = inputs.args(extra);
        set(&mut args);
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        let result = {
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: "6",
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(&mut sink),
                }],
                diagnostics: &mut diagnostics,
                later_warnings: None,
                hits: None,
                observer: None,
            };
            run_local_blastn(
                args,
                &queries,
                &crate::api::QueryEnd::Input,
                &subjects,
                &mut outputs,
            )
        };
        let error = format!("{:#}", result.expect_err(message));
        assert!(error.contains(message), "{message}: {error}");
        assert!(sink.is_empty() && diagnostics.is_empty(), "{message}");
    }
    // The database entry checks the task first too (before the database is opened).
    let mut args = inputs.args(&[]);
    args.db = Some("no_such_database".into());
    args.task = "blastz".into();
    let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
    let result = {
        let mut outputs = ReportOutputs {
            formats: vec![FormatOutput {
                outfmt: "6",
                refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                sink: OutputSink::Writer(&mut sink),
            }],
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        crate::api::local_blast::run_local_database_blastn(
            args,
            &queries,
            &crate::api::QueryEnd::Input,
            crate::blastinput::ncbi_environment::ApplicationSettings::default(),
            &mut outputs,
        )
    };
    let error = format!("{:#}", result.expect_err("refused"));
    assert!(
        error.contains("'blastz' is not a supported task"),
        "{error}"
    );
}

/// The exit code and the message of NCBI's error that a run ended with.
fn native_error(result: anyhow::Result<()>) -> (i32, String) {
    let error = result.expect_err("the run fails");
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .unwrap_or_else(|| panic!("NCBI's error: {error:#}"));
    (
        native.exit,
        String::from_utf8(native.message.clone()).unwrap(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:649-652
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
            "LOSAT engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n"
                .to_string()
        )
    );
    assert!(diagnostics.is_empty(), "{diagnostics:?}");
    let prolog = String::from_utf8_lossy(&outputs[0]);
    assert!(
        prolog.starts_with(&format!("BLASTN (LOSAT {})", env!("CARGO_PKG_VERSION")))
            && !prolog.contains("Query="),
        "{prolog}"
    );
    assert!(outputs[1].is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:153-156
// With scores that have no Karlin-Altschul values, each query of the batch gets the error,
// and a query without residues keeps its warning after it (`scoring::karlin_error`; NCBI
// 2.17.0 with `-reward 1 -penalty -6`).
#[test]
fn a_query_without_residues_keeps_its_warning_in_a_table_error() {
    let genome = fixture_sequence("LC738884.fasta");
    let inputs = Inputs {
        query: TempFasta::new(
            "blastn_kaempty_query.fna",
            &[("q1", b""), ("q2", &genome[0..500]), ("q3 three", b"")],
        ),
        subject: TempFasta::new("blastn_kaempty_subject.fna", &[("s1", &genome[0..5_000])]),
    };
    let (outputs, diagnostics, result) =
        run_local_outputs(&inputs, &["0", "6"], &["-reward", "1", "-penalty", "-6"]);
    let error = "Error: Substitution scores 1 and -6 are not supported ";
    let warning = "Warning: Sequence contains no data ";
    assert_eq!(
        native_error(result),
        (
            3,
            format!("LOSAT engine error: {error}{warning}{error}{error}{warning}\n")
        )
    );
    assert!(diagnostics.is_empty(), "{diagnostics:?}");
    let prolog = String::from_utf8_lossy(&outputs[0]);
    assert!(
        prolog.starts_with(&format!("BLASTN (LOSAT {})", env!("CARGO_PKG_VERSION")))
            && !prolog.contains("Query="),
        "{prolog}"
    );
    assert!(outputs[1].is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
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
                "LOSAT engine error: The average subject length is too short\n".to_string()
            )
        );
        assert_eq!(String::from_utf8_lossy(&diagnostics), warnings);
        assert!(String::from_utf8_lossy(&outputs[0])
            .starts_with(&format!("BLASTN (LOSAT {})", env!("CARGO_PKG_VERSION"))));
        assert!(outputs[1].is_empty());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-460
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
// The command line checks `-strand` at parsing; `run_local` takes the arguments as a struct, so
// it checks the value itself: an unknown strand is an argument error and nothing is written,
// while `plus` and `minus` search and `both` is the plain search.
#[test]
fn run_local_rejects_an_unknown_strand() {
    let inputs = Inputs::new();
    let queries = reader_records(
        &inputs.query.0,
        ReaderConfig::query("BLASTN", false, true).with_defline(true),
    );
    let subjects = reader_records(
        &inputs.subject.0,
        ReaderConfig::subject("BLASTN", false, true).with_defline(true),
    );
    for strand in ["foo", "", "PLUS", "plus "] {
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs =
            ReportOutputs::single("6", OutputSink::Writer(&mut sink), &mut diagnostics);
        let mut args = inputs.args(&[]);
        args.strand = strand.to_string();
        let result = run_local_blastn(
            args,
            &queries,
            &crate::api::QueryEnd::Input,
            &subjects,
            &mut outputs,
        );
        drop(outputs);
        let error = result
            .expect_err("an unknown strand is an error")
            .to_string();
        assert_eq!(
            error,
            format!(
                "Argument \"strand\". Illegal value, expected `both', `minus', `plus':  `{strand}'"
            )
        );
        assert!(sink.is_empty() && diagnostics.is_empty(), "{strand:?}");
    }
    let both = inputs.run(&["6"], &["-strand", "both"], false);
    let plain = inputs.run(&["6"], &[], false);
    assert_eq!(both.outputs, plain.outputs);
    for strand in ["plus", "minus"] {
        let run = inputs.run(&["6"], &["-strand", strand], false);
        assert!(!run.outputs[0].is_empty(), "{strand}");
    }
}

/// The extensions of each batch of a search, as the engine reports them.
#[derive(Default)]
struct Extensions(Vec<i32>);

impl FormatObserver for Extensions {
    fn hsp_begin(&mut self, _format: usize, _hsp: HspIndex) {}

    fn hsp_end(&mut self, _format: usize, _hsp: HspIndex) {}

    fn batch_extensions(&mut self, extensions: i32) {
        self.0.push(extensions);
    }
}

#[test]
fn blastn_batches_are_those_of_the_search_of_all_queries() {
    // `CBatchSizeMixer` sizes BLASTN's batches by the extensions of the batch before.
    let once = std::fs::read(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/small_test.fasta"),
    )
    .expect("fixture");
    let copied = once.repeat(30);
    let queries = read_fasta(Program::Blastn, Role::Query, &copied).unwrap();
    let subjects = read_fasta(Program::Blastn, Role::Subject, &once).unwrap();
    for task in ["blastn", "megablast"] {
        let cli = crate::cli::try_parse_from::<crate::cli::Cli, _, _>([
            "losat", "blastn", "-query", "q", "-subject", "s", "-task", task,
        ])
        .unwrap();
        let crate::cli::Commands::Blastn(args) = cli.command else {
            panic!("a BLASTN command")
        };
        let mut extensions = Extensions::default();
        let mut report = Vec::new();
        let mut diagnostics = Vec::new();
        {
            let mut on_hits = |_| {};
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: "6",
                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(&mut report),
                }],
                diagnostics: &mut diagnostics,
                later_warnings: None,
                hits: Some(&mut on_hits),
                observer: Some(&mut extensions),
            };
            run_local_blastn(
                args,
                &queries,
                &crate::api::QueryEnd::Input,
                &subjects,
                &mut outputs,
            )
            .unwrap();
        }
        let search = validate(Program::Blastn, &["-task", task]).unwrap();
        let batches = run_batches(&search, ApiInputs::new(&queries, &subjects))
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(batches.len(), extensions.0.len(), "{task}");
        assert!(batches.len() > 1, "{task}");
    }
}

// Line RG item T7-W (review 3 of line RG T7): with `-parse_deflines`, LOSAT refuses outfmt 0 of a
// `-subject` record whose description id is longer than 65 bytes, where NCBI BLAST+ 2.17.0 pads
// the description line without end (`db_subject::check_pairwise_records`). The refusal is outfmt
// 0's: a call in which outfmt 0 fails alone (`FormatRefusal::FailsAlone`, a LOSAT Web format that
// the argv does not ask for) writes outfmt 6, and one in which it fails the call (the command
// line's one format) fails with it.
#[test]
fn an_outfmt_0_refusal_of_a_parse_deflines_record_fails_outfmt_0_alone_or_the_call() {
    use crate::api::local_blast::FormatRefusal;
    let a = fixture_sequence("LC738884.fasta");
    let long =
        "gnl|ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGH|ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJ long id";
    let query = TempFasta::new("t7w_parse_q.fna", &[("q1", &a[0..2000])]);
    let subject = TempFasta::new("t7w_parse_s.fna", &[(long, &a[0..2000])]);
    let queries = reader_records(
        &query.0,
        ReaderConfig::query("BLASTN", false, true).with_parse_deflines(true),
    );
    let subjects = reader_records(
        &subject.0,
        ReaderConfig::subject("BLASTN", false, true).with_parse_deflines(true),
    );
    let args = || -> BlastnArgs {
        let (query, subject) = (
            query.0.display().to_string(),
            subject.0.display().to_string(),
        );
        let argv = [
            "losat",
            "blastn",
            "-query",
            &query,
            "-subject",
            &subject,
            "-parse_deflines",
        ];
        let cli: Cli = try_parse_from(argv).expect("blastn argv");
        match cli.command {
            Commands::Blastn(args) => args,
            _ => unreachable!("blastn subcommand"),
        }
    };
    let refusal = "description id longer than 65 bytes";
    for alone in [true, false] {
        let (mut six, mut zero, mut diagnostics) = (Vec::new(), Vec::new(), Vec::new());
        let mut outputs = ReportOutputs {
            formats: vec![
                FormatOutput {
                    outfmt: "6",
                    sink: OutputSink::Writer(&mut six),
                    refusal: FormatRefusal::FailsCall,
                },
                FormatOutput {
                    outfmt: "0",
                    sink: OutputSink::Writer(&mut zero),
                    refusal: if alone {
                        FormatRefusal::FailsAlone
                    } else {
                        FormatRefusal::FailsCall
                    },
                },
            ],
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        let result = run_local_blastn(
            args(),
            &queries,
            &crate::api::QueryEnd::Input,
            &subjects,
            &mut outputs,
        );
        if alone {
            result.expect("outfmt 0 fails alone");
            let error = outputs.formats[1]
                .refusal
                .take_refused()
                .expect("outfmt 0 refused");
            assert!(format!("{error:#}").contains(refusal), "{error:#}");
            assert!(outputs.formats[0].refusal.take_refused().is_none());
            drop(outputs);
            let rows = String::from_utf8_lossy(&six).into_owned();
            assert!(
                rows.starts_with("q1\t") && rows.contains("\t100.000\t2000\t"),
                "{rows}"
            );
        } else {
            let error = result.expect_err("outfmt 0 fails the call");
            assert!(format!("{error:#}").contains(refusal), "{error:#}");
        }
    }
}
