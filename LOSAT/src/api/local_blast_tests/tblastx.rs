//! The shared TBLASTX entry `run_local` in process: what the entry itself does, not what the
//! command line writes. The tests that compare the entry's reports with the `losat` binary go
//! through `losat::api` (`LOSAT/tests/run_local_tblastx.rs`).

use super::support::{
    assert_observer_ranges, assert_tabular_ranges, fixture_sequence, reader_records, run_formats,
    Run, TempFasta,
};
use crate::algorithm::tblastx::TblastxArgs;
use crate::api::local_blast::{run_local_tblastx, FormatOutput, OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::{FastaRecord, ReaderConfig};
use crate::cli::{try_parse_from, Cli, Commands};

/// The last three are the custom tabular reports of the shared writer (line CLI C2, PR 2).
const FORMATS: [&str; 13] = [
    "0",
    "6",
    "7",
    "8",
    "9",
    "12",
    "5",
    "15",
    "16",
    "11",
    "10 delim=| std qseq sseq sframe",
    "7 std qseq sseq btop positive ppos frames qcovs qcovhsp",
    "20 qseqid sseqid evalue sstrand",
];

/// Two queries against three subjects cut from two genomes; the first query overlaps
/// two subjects.
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
                "tblastx_query.fna",
                &[("q_a0", &a[0..3_000]), ("q_a100k", &a[100_000..103_000])],
            ),
            subject: TempFasta::new(
                "tblastx_subject.fna",
                &[
                    ("s_a1k", &a[1_000..6_000]),
                    ("s_a99k", &a[99_000..105_000]),
                    ("s_a0", &a[0..2_000]),
                    ("s_b", &b[0..5_000]),
                ],
            ),
        }
    }

    fn args(&self, extra: &[&str]) -> TblastxArgs {
        let query = self.query.0.display().to_string();
        let subject = self.subject.0.display().to_string();
        let mut argv = vec![
            "losat", "tblastx", "-query", &query, "-subject", &subject, "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv).expect("tblastx argv");
        match cli.command {
            Commands::Tblastx(args) => args,
            _ => unreachable!("tblastx subcommand"),
        }
    }

    fn run(&self, formats: &[&str], extra: &[&str], observe: bool) -> Run {
        let queries = queries(&self.query.0);
        let subjects = subjects(&self.subject.0);
        run_formats(formats, observe, |outputs| {
            run_local_tblastx(self.args(extra), &queries, &subjects, outputs)
        })
    }
}

/// The query records of a FASTA file as NCBI's reader reads them for TBLASTX.
fn queries(path: &std::path::Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::query("TBLASTX", false, true))
}

/// The subject records of a FASTA file as NCBI's reader reads them for TBLASTX.
fn subjects(path: &std::path::Path) -> Vec<FastaRecord> {
    reader_records(path, ReaderConfig::subject("TBLASTX", false, true))
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
        &["-num_threads", "2"][..],
        &["-max_target_seqs", "1"][..],
    ] {
        let result = inputs.run(&FORMATS, extra, true);
        assert_observer_ranges(&result, [0, 1, 2], &format!("tblastx {extra:?}"));
        // One row per HSP in the custom tabular report too (outfmt 10, no comment lines).
        assert_eq!(
            assert_tabular_ranges(&result, 10, None, "tblastx outfmt 10"),
            result.hits.len()
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
// ```c
// arg_desc.AddDefaultKey(kArgOutputFormat, "format",
//                        kOutputFormatDescription,
//                        CArgDescriptions::eString,
//                        NStr::IntToString(dft_outfmt));
// ```
#[test]
fn unsupported_formats_fail_before_searching() {
    let inputs = Inputs::new();
    let queries = queries(&inputs.query.0);
    let subjects = subjects(&inputs.subject.0);
    // A delimiter that a field separates its entries with is NCBI's error
    // (blast_args.cpp:2753-2792, `s_ValidateCustomDelim`).
    for outfmt in ["2", "6 delim=; qseqid staxids", "abc"] {
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
        let result = run_local_tblastx(inputs.args(&[]), &queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(result.is_err(), "outfmt {outfmt:?} must be rejected");
        assert!(
            valid.is_empty() && invalid.is_empty() && diagnostics.is_empty(),
            "nothing is written when a requested format is invalid"
        );
    }
}

// NCBI derives the hit list size per output format (blast_args.cpp:2894-2973): with
// `-num_descriptions 3` it is 250 for outfmt 0 and 500 for outfmt 6. One run searches once,
// so formats that need different hit list sizes are refused before anything is written
// (`blastinput/report_counts.rs`), while formats of one hit list size still share a run.
#[test]
fn formats_needing_different_hit_list_sizes_are_refused() {
    let inputs = Inputs::new();
    let queries = queries(&inputs.query.0);
    let subjects = subjects(&inputs.subject.0);
    let attempt = |formats: &[&str], extra: &[&str]| {
        let mut sinks: Vec<Vec<u8>> = formats.iter().map(|_| Vec::new()).collect();
        let mut diagnostics = Vec::new();
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
        let result = run_local_tblastx(inputs.args(extra), &queries, &subjects, &mut outputs);
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

/// Runs `run_local` with one writer per format and returns the outputs, the diagnostics
/// and the result.
fn run_local_outputs(
    inputs: &Inputs,
    formats: &[&str],
    extra: &[&str],
) -> (Vec<Vec<u8>>, Vec<u8>, anyhow::Result<()>) {
    let queries = queries(&inputs.query.0);
    let subjects = subjects(&inputs.subject.0);
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
        run_local_tblastx(inputs.args(extra), &queries, &subjects, &mut outputs)
    };
    (sinks, diagnostics, result)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1312-1318
// ```c
//     if (options->word_size <= 0)
//     {
//         if ( !Blast_ProgramIsRpsBlast(program_number)) {
//             Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
//                                      "Word-size must be greater than zero");
//             return BLASTERR_OPTION_VALUE_INVALID;
// ```
// The command line's `-word_size` is at least 2 (clap's `protein_word_size`); a caller of
// `run_local` that sets the field itself gets the error of NCBI's core. A word size of 0 with
// `-lcase_masking` and a lowercase subject took the reference scan, whose loop over the first
// `word_length - 1` letters wrapped and read past the subject.
#[test]
fn a_word_size_of_zero_is_ncbis_core_error() {
    let a = fixture_sequence("LC738884.fasta");
    let query = TempFasta::new("tblastx_word0_query.fna", &[("q", &a[0..3_000])]);
    let lower = a[1_000..6_000].to_ascii_lowercase();
    let subject = TempFasta::new("tblastx_word0_subject.fna", &[("s", &lower[..])]);
    let (q, s) = (
        query.0.display().to_string(),
        subject.0.display().to_string(),
    );
    let cli: Cli = try_parse_from([
        "losat",
        "tblastx",
        "-query",
        &q,
        "-subject",
        &s,
        "-outfmt",
        "6",
        "-lcase_masking",
    ])
    .expect("tblastx argv");
    let Commands::Tblastx(mut args) = cli.command else {
        unreachable!("tblastx subcommand")
    };
    args.word_size = Some(0);
    let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
    let result = {
        let mut outputs = ReportOutputs {
            formats: vec![FormatOutput {
                outfmt: "6",
                sink: OutputSink::Writer(&mut sink),
            }],
            diagnostics: &mut diagnostics,
            hits: None,
            observer: None,
        };
        run_local_tblastx(
            args,
            &queries(&query.0),
            &subjects(&subject.0),
            &mut outputs,
        )
    };
    let (_, message) = native_error(result);
    assert!(
        message.contains("Word-size must be greater than zero"),
        "{message}"
    );
    assert!(sink.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1312-1325 (quoted above):
// NCBI's core refuses a word size of 0 or less, and below 4 only for a nucleotide comparison,
// which TBLASTX (translated, `Blast_ProgramIsNucleotide`, blast_program.c:82-84) is not. So a
// caller of `run_local` that sets a word size of 1 (the command line's `-word_size` is at least 2)
// gets a search: through the reference scan (`-lcase_masking` and a lowercase subject), the
// default one and at two threads. A debug build checks the unchecked reads' preconditions.
#[test]
fn a_word_size_of_one_searches() {
    let a = fixture_sequence("LC738884.fasta");
    let query = TempFasta::new("tblastx_word1_query.fna", &[("q", &a[0..1_500])]);
    let lower = a[1_000..3_000].to_ascii_lowercase();
    let subject = TempFasta::new(
        "tblastx_word1_subject.fna",
        &[("s_lower", &lower[..]), ("s", &a[2_000..4_000])],
    );
    let (q, s) = (
        query.0.display().to_string(),
        subject.0.display().to_string(),
    );
    for extra in [&[][..], &["-lcase_masking"], &["-num_threads", "2"]] {
        let mut argv = vec![
            "losat", "tblastx", "-query", &q, "-subject", &s, "-outfmt", "6",
        ];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv).expect("tblastx argv");
        let Commands::Tblastx(mut args) = cli.command else {
            unreachable!("tblastx subcommand")
        };
        args.word_size = Some(1);
        let (mut sink, mut diagnostics) = (Vec::new(), Vec::new());
        let result = {
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt: "6",
                    sink: OutputSink::Writer(&mut sink),
                }],
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            run_local_tblastx(
                args,
                &queries(&query.0),
                &subjects(&subject.0),
                &mut outputs,
            )
        };
        result.unwrap_or_else(|error| panic!("{extra:?}: {error:#}"));
        // The lowercase subject repeats the query's last 500 letters: the searches find them,
        // except where `-lcase_masking` masks that subject.
        if extra != ["-lcase_masking"] {
            assert!(!sink.is_empty(), "{extra:?}: no rows");
        }
    }
}

/// The exit code and the message of NCBI's error that a run ended with.
fn native_error(result: anyhow::Result<()>) -> (i32, String) {
    let error = result.expect_err("the run fails");
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .unwrap_or_else(|| panic!("NCBI's error: {error:#}"));
    (native.exit, native.message.clone())
}

// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312 (x_CleanAndCompress) and
// 4066 (`NStr::HtmlDecode`): NCBI reads past the end of the title `, ,` (and crashes when
// such a subject has hits) and decodes `&amp;` in outfmt 0; the tabular formats print
// only the ids. LOSAT writes the title `, ` stopped at the end of the string (approved
// exception 2 of PD-LOSAT-NCBI-DEFECTS) and the decoded title (`report/defline.rs`).
#[test]
fn outfmt0_titles_are_written_as_ncbi_makes_them() {
    let inputs = Inputs::new();
    let queries = queries(&inputs.query.0);
    let sequence = fixture_sequence("LC738884.fasta");
    for (defline, heading) in [
        (", ,", "\n> , \n"),
        ("s &amp; t", "\n> s & t\n"),
        ("R&D; x", "\n> R&D; x\n"),
        ("a&foo;b", "\n> a&foo;b\n"),
        ("a&amp;", "\n> a&amp\n"),
    ] {
        let subject = TempFasta::new(
            "tblastx_title_subject.fna",
            &[(defline, &sequence[1_000..6_000])],
        );
        let subjects = subjects(&subject.0);
        for outfmt in FORMATS {
            let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
            let mut outputs =
                ReportOutputs::single(outfmt, OutputSink::Writer(&mut report), &mut diagnostics);
            let result = run_local_tblastx(inputs.args(&[]), &queries, &subjects, &mut outputs);
            drop(outputs);
            // BLAST XML: NCBI crashes on the punctuation-only title of a subject with hits
            // (outfmt 0, 5, 15 and 16 alike), so LOSAT refuses it and writes nothing.
            if matches!(outfmt, "5" | "15" | "16") && defline == ", ," {
                let error = result.expect_err("outfmt 5, 15 and 16 refuse the punctuation title");
                assert!(error.to_string().contains("crashes on"), "{error}");
                assert!(report.is_empty(), "nothing is written");
                continue;
            }
            result.unwrap_or_else(|error| panic!("outfmt {outfmt} {defline:?}: {error}"));
            if outfmt == "0" {
                let report = String::from_utf8(report).expect("UTF-8 report");
                assert!(report.contains(heading), "{defline:?}: {report}");
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:649-652
// ```c
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// A batch whose queries all have fewer letters than a codon stops NCBI after the prolog
// when one of them has no letters, with one message per query without letters; a batch of
// queries of one or two letters only is not searched.
#[test]
fn a_batch_without_valid_contexts_and_with_queries_without_letters_stops() {
    let genome = fixture_sequence("LC738884.fasta");
    let subject = TempFasta::new("tblastx_no_data_subject.fna", &[("s1", &genome[0..5_000])]);
    for (records, message) in [
        (
            vec![("q1", &b""[..]), ("q2 two", &b""[..])],
            "BLAST engine error: Warning: Sequence contains no data Warning: Sequence contains no data \n",
        ),
        (
            vec![("q1", &b"AC"[..]), ("q2 two", &b""[..])],
            "BLAST engine error: Warning: Sequence contains no data \n",
        ),
    ] {
        let inputs = Inputs {
            query: TempFasta::new("tblastx_no_data_query.fna", &records),
            subject: TempFasta::new("tblastx_no_data_subject.fna", &[("s1", &genome[0..5_000])]),
        };
        let (outputs, diagnostics, result) = run_local_outputs(&inputs, &["0", "6"], &[]);
        assert_eq!(native_error(result), (3, message.to_string()));
        assert!(diagnostics.is_empty(), "{diagnostics:?}");
        let prolog = String::from_utf8_lossy(&outputs[0]);
        assert!(
            prolog.starts_with("TBLASTX 2.17.0+") && !prolog.contains("Query="),
            "{prolog}"
        );
        assert!(outputs[1].is_empty());
    }
    let inputs = Inputs {
        query: TempFasta::new("tblastx_short_query.fna", &[("q1", b"AC"), ("q2", b"A")]),
        subject,
    };
    let (outputs, diagnostics, result) = run_local_outputs(&inputs, &["7"], &[]);
    result.expect("a batch of short queries is reported unsearched");
    let warnings = String::from_utf8_lossy(&diagnostics);
    assert_eq!(
        warnings
            .matches("Could not calculate ungapped Karlin-Altschul parameters")
            .count(),
        2,
        "{warnings}"
    );
    let report = String::from_utf8_lossy(&outputs[0]);
    assert!(
        report.ends_with("# BLAST processed 2 queries\n"),
        "{report}"
    );
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
// Subjects with fewer letters than subjects (every subject without residues; AUTHORITY
// §E3, BI-44) stop NCBI when it searches the first batch with a valid context, after the
// subjects' warnings and the prolog.
#[test]
fn subjects_with_fewer_letters_than_subjects_stop_the_search() {
    let genome = fixture_sequence("LC738884.fasta");
    for (subjects, warnings) in [
        (
            vec![("s1", &b"A"[..]), ("s2", &b""[..])],
            "Warning: [tblastx] Subject_2 s2: Subject sequence contains no data\n",
        ),
        (
            vec![("s1 one", &b""[..]), ("s2", &b""[..])],
            "Warning: [tblastx] Subject_1 s1 one: Subject sequence contains no data\n\
             Warning: [tblastx] Subject_2 s2: Subject sequence contains no data\n",
        ),
    ] {
        let inputs = Inputs {
            query: TempFasta::new(
                "tblastx_short_subjects_query.fna",
                &[("q1", &genome[0..500])],
            ),
            subject: TempFasta::new("tblastx_short_subjects.fna", &subjects),
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
        assert!(String::from_utf8_lossy(&outputs[0]).starts_with("TBLASTX 2.17.0+"));
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
// NCBI sets up as a record without residues.
#[test]
fn intervals_without_letters_are_records_without_residues() {
    let genome = fixture_sequence("LC738884.fasta");
    let inputs = Inputs {
        query: TempFasta::new(
            "tblastx_empty_interval_query.fna",
            &[("q1", &genome[0..2_000]), ("q2 short", &genome[0..300])],
        ),
        subject: TempFasta::new(
            "tblastx_empty_interval_subject.fna",
            &[("s1", &genome[0..5_000]), ("s2 short", &genome[0..300])],
        ),
    };
    let run = inputs.run(&["6"], &["-subject_loc", "301-1000"], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [tblastx] Subject_2 s2 short: Subject sequence contains no data\n"
    );
    let rows = String::from_utf8_lossy(&run.outputs[0]);
    assert!(rows.lines().all(|row| row.contains("\ts1\t")), "{rows}");
    assert!(!rows.is_empty());
    let run = inputs.run(&["7"], &["-query_loc", "301-1000"], false);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [tblastx] Query_2 q2 short: Sequence contains no data \n"
    );
    let report = String::from_utf8_lossy(&run.outputs[0]);
    assert!(
        report.contains("# Query: q2 short\n") && report.contains("# BLAST processed 2 queries\n"),
        "{report}"
    );
}

// NCBI reference: c++/src/objtools/readers/fasta.cpp:919-935 (CFastaReader: a residue that is
// not an IUPAC nucleotide letter is removed with a warning). NCBI's check of the hit saving
// options (blast_options.c:1518-1523) comes before `Query is Empty!`.
// This is the half about the entry: a hand-set expect value of 0 (or one that underflows to 0) is
// refused with and without queries, and nothing is written. The command line's half (the
// warning for the invalid residue) is in `LOSAT/tests/run_local_tblastx.rs`.
#[test]
fn residues_are_read_as_ncbis_reader_reads_them() {
    let inputs = Inputs::new();
    let queries = queries(&inputs.query.0);
    let subjects = subjects(&inputs.subject.0);
    let run = |extra: &[&str], queries: &[FastaRecord]| {
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        let mut outputs =
            ReportOutputs::single("6", OutputSink::Writer(&mut report), &mut diagnostics);
        let result = run_local_tblastx(inputs.args(extra), queries, &subjects, &mut outputs);
        drop(outputs);
        assert!(report.is_empty() && diagnostics.is_empty());
        result.expect_err("rejected").to_string()
    };
    for evalue in ["0", "1e-400"] {
        for queries in [&queries[..], &[]] {
            assert_eq!(
                run(&["-evalue", evalue], queries),
                "BLAST query/options error: expect value or cutoff score must be greater than zero\nPlease refer to the BLAST+ user manual.\n"
            );
        }
    }
}
