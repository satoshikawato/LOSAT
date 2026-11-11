//! The BLASTN search of `LOSAT::api` (which calls the shared entry `run_local`) writes several
//! output formats from one search. Every format must be byte-identical to the CLI run that
//! requests only that format, and the warnings must equal that run's standard error. (Identity of
//! outfmt 0 with NCBI is checked by docs/evidence/losat_web_e2a/check_losat.py over the frozen
//! fixtures.) (Identity with the previous release over the regression manifests is checked by
//! docs/evidence/losat_web_e1a/capture_outputs.py.) The entry's own contract (its observer, its
//! refusals, its errors) is checked in the crate, `LOSAT/src/api/local_blast_tests/blastn.rs`.

mod run_local_support;

use std::process::Command;

use run_local_support::{cli_out_path, fixture_sequence, Files, TempFasta};
use LOSAT::api::{Outcome, Program};

const FORMATS: [&str; 12] = [
    "0", "6", "7", "8", "9", "12", "5", "15", "16", "11", "17", "17 SQ SR",
];

/// The SAM document without the CL tag of its @PG line: the bytes from the `\tCL:` that follows
/// `\n@PG\tID:0\tVN:2.17.0+` up to, not including, the first following `\tPN:blastn\n`; a
/// document without the tag stays as it is (the rule of `db_search_fixtures.py strip_sam_cl`).
fn without_sam_cl(document: &[u8]) -> Vec<u8> {
    let find = |haystack: &[u8], needle: &[u8], from: usize| {
        haystack[from..]
            .windows(needle.len())
            .position(|window| window == needle)
            .map(|at| at + from)
    };
    let head: &[u8] = b"\n@PG\tID:0\tVN:2.17.0+";
    let Some(at) = find(document, &[head, b"\tCL:"].concat(), 0) else {
        return document.to_vec();
    };
    let start = at + head.len();
    let Some(end) = find(document, b"\tPN:blastn\n", start + 4) else {
        return document.to_vec();
    };
    [&document[..start], &document[end..]].concat()
}

#[test]
fn the_sam_cl_rule_hides_the_command_line_and_nothing_else() {
    let with: &[u8] =
        b"@HD\tVN:1.2\n@SQ\tSN:Query_1\tLN:3\n@PG\tID:0\tVN:2.17.0+\tCL:/x/blastn -a \tPN:blastn\nr\t0\n";
    let without: &[u8] =
        b"@HD\tVN:1.2\n@SQ\tSN:Query_1\tLN:3\n@PG\tID:0\tVN:2.17.0+\tPN:blastn\nr\t0\n";
    let other_cl: &[u8] =
        b"@HD\tVN:1.2\n@SQ\tSN:Query_1\tLN:3\n@PG\tID:0\tVN:2.17.0+\tCL:losat blastn -b \tPN:blastn\nr\t0\n";
    assert_eq!(without_sam_cl(with), without);
    assert_eq!(without_sam_cl(other_cl), without);
    assert_eq!(without_sam_cl(without), without);
    // Any other byte is compared: a record, the header, the program name.
    let changed_record = with
        .iter()
        .map(|&b| if b == b'r' { b's' } else { b })
        .collect::<Vec<u8>>();
    assert_ne!(without_sam_cl(&changed_record), without);
    let changed_sq = String::from_utf8(with.to_vec())
        .unwrap()
        .replace("LN:3", "LN:4");
    assert_ne!(without_sam_cl(changed_sq.as_bytes()), without);
    let changed_pn = String::from_utf8(with.to_vec())
        .unwrap()
        .replace("PN:blastn", "PN:blastx");
    assert_ne!(without_sam_cl(changed_pn.as_bytes()), without);
}

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

    /// The library's search of `formats` with the options `extra`, with the files' paths as the
    /// names that the reports show, as the CLI shows them.
    fn collect(&self, formats: &[&str], extra: &[&str]) -> Outcome {
        Files::new(Program::Blastn, &self.query.0, &self.subject.0).collect(formats, extra)
    }

    /// The CLI output file and standard error of one single-format run, with a SAM document's
    /// CL tag removed (`without_sam_cl`).
    fn cli(&self, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let (bytes, stderr) = self.cli_raw(outfmt, extra);
        // A SAM document of the CLI shows its command line in the @PG line, which `run_local`
        // has not: the one rule of the SAM rows removes the tag (`without_sam_cl`).
        let bytes = if outfmt.split_whitespace().next() == Some("17") {
            without_sam_cl(&bytes)
        } else {
            bytes
        };
        (bytes, stderr)
    }

    /// The CLI output file and standard error of one single-format run, as written.
    fn cli_raw(&self, outfmt: &str, extra: &[&str]) -> (Vec<u8>, Vec<u8>) {
        let out = cli_out_path(&self.query.0, outfmt);
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
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

/// The frozen fixtures and the oracle compare a SAM document without its CL tag (both sides), so
/// they cannot see a CLI that drops the tag or a library that writes one: the CLI document has
/// `CL:` with its own command line (the program path first), a `run_local` document has none.
#[test]
fn the_cli_writes_a_cl_tag_and_the_library_none() {
    let inputs = Inputs::new();
    let head: &[u8] = b"\n@PG\tID:0\tVN:2.17.0+";
    let tag = |document: &[u8]| -> Option<String> {
        let at = document
            .windows(head.len())
            .position(|window| window == head)?
            + head.len();
        let line = document[at..].split(|&byte| byte == b'\n').next()?;
        line.starts_with(b"\tCL:")
            .then(|| String::from_utf8_lossy(line).into_owned())
    };
    for outfmt in ["17", "17 SQ SR"] {
        let (cli, _) = inputs.cli_raw(outfmt, &[]);
        let line =
            tag(&cli).unwrap_or_else(|| panic!("the CLI document of {outfmt} has no CL tag"));
        let query = inputs.query.0.display().to_string();
        assert!(
            line.starts_with(&format!(
                "\tCL:{} blastn -query {query} -subject ",
                env!("CARGO_BIN_EXE_losat")
            )) && line.ends_with(&format!(
                " -outfmt {outfmt} -out {} \tPN:blastn",
                cli_out_path(&inputs.query.0, outfmt).display()
            )),
            "{line}"
        );
        let library = inputs.collect(&[outfmt], &[]).reports.remove(0).bytes;
        assert!(
            tag(&library).is_none()
                && library
                    .windows(head.len() + 11)
                    .any(|window| window == [head, b"\tPN:blastn\n"].concat()),
            "the library document of {outfmt} has a CL tag or no @PG line"
        );
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
        // The culling writer and pipe (`-culling_limit`) act on every path of the shared entry.
        &["-culling_limit", "1"][..],
    ] {
        let together = inputs.collect(&FORMATS, extra);
        let rows = together.reports[1].spans.iter().flatten().count();
        assert_eq!(
            together.hsps.len(),
            rows,
            "one hit record per outfmt 6 row {extra:?}"
        );
        row_counts.push(rows);
        for (index, outfmt) in FORMATS.iter().enumerate() {
            let alone = inputs.collect(&[outfmt], extra);
            assert!(
                together.reports[index].bytes == alone.reports[0].bytes,
                "run_collect outfmt {outfmt} {extra:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, extra);
            assert!(
                together.reports[index].bytes == cli_output,
                "CLI outfmt {outfmt} {extra:?}"
            );
            assert_eq!(
                together.diagnostics,
                cli_stderr,
                "warnings of {} formats equal those of the single-format CLI run {extra:?}",
                FORMATS.len()
            );
        }
    }
    assert!(row_counts[0] > 0, "the fixture must produce hits");
    assert!(
        row_counts[1] < row_counts[0],
        "-max_target_seqs 1 must drop a subject of this fixture: {row_counts:?}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786
// ```c++
//     if (m_FormatType == CFormattingArgs::eTabular ||
//         m_FormatType == CFormattingArgs::eTabularWithComments ||
//         m_FormatType == CFormattingArgs::eCommaSeparatedValues ||
//         m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
// ```
// The tabular reports with a field list or a delimiter, and outfmt 10 and 20 (line CLI C2), of
// one search: each equals its run alone and the CLI run of that format, beside the default 6.
#[test]
fn custom_tabular_formats_of_one_search_match_the_cli_run_of_that_format() {
    const CUSTOM: [&str; 5] = [
        "6",
        "7 std qcovs qcovus qcovhsp sstrand staxid slen",
        "10 qseqid sseqid qseq sseq btop nident positive ppos",
        "20 delim=| qseqid sallseqid sacc stitle frames qlen",
        "6 delim=, std",
    ];
    let inputs = Inputs::new();
    for extra in [
        &[][..],
        &["-task", "blastn"][..],
        &["-num_threads", "2"][..],
    ] {
        let together = inputs.collect(&CUSTOM, extra);
        for (index, outfmt) in CUSTOM.iter().enumerate() {
            assert!(
                !together.reports[index].bytes.is_empty(),
                "outfmt {outfmt:?}"
            );
            let alone = inputs.collect(&[outfmt], extra);
            assert!(
                together.reports[index].bytes == alone.reports[0].bytes,
                "run_collect outfmt {outfmt:?} {extra:?}"
            );
            let (cli_output, cli_stderr) = inputs.cli(outfmt, extra);
            assert!(
                together.reports[index].bytes == cli_output,
                "CLI outfmt {outfmt:?} {extra:?}"
            );
            assert_eq!(together.diagnostics, cli_stderr, "{outfmt:?} {extra:?}");
        }
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
        for outfmt in FORMATS {
            let result = Files::new(Program::Blastn, &inputs.query.0, &inputs.subject.0)
                .try_collect(&[outfmt], &[]);
            let outcome = result.unwrap_or_else(|error| {
                panic!("{n_residues} N, outfmt {outfmt}: {error}");
            });
            let (cli_output, cli_stderr) = inputs.cli(outfmt, &[]);
            assert_eq!(
                outcome.reports[0].bytes, cli_output,
                "{n_residues} N, outfmt {outfmt}"
            );
            assert_eq!(
                outcome.diagnostics, cli_stderr,
                "{n_residues} N, outfmt {outfmt}"
            );
        }
    }
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
    let run = inputs.collect(&FORMATS, &[]);
    assert_eq!(
        String::from_utf8_lossy(&run.diagnostics),
        "Warning: [blastn] Subject_2 s2 no letters: Subject sequence contains no data\n\
         Warning: [blastn] Query_2 q2 empty two: Sequence contains no data \n"
    );
    let pairwise = String::from_utf8_lossy(&run.reports[0].bytes);
    assert!(
        pairwise.contains("3 sequences; 9,000 total letters"),
        "{pairwise}"
    );
    assert!(
        pairwise.contains("Query= q2 empty two\n\nLength=0\n\n\n***** No hits found *****\n"),
        "{pairwise}"
    );
    let tabular = String::from_utf8_lossy(&run.reports[1].bytes);
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
    let commented = String::from_utf8_lossy(&run.reports[2].bytes);
    assert!(
        commented.contains("# Query: q2 empty two\n# Database: ")
            && commented.contains("# BLAST processed 3 queries\n"),
        "{commented}"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eAsnText
    //       ...
    //       || m_FormatType == CFormattingArgs::eJsonSeqalign
    //       || m_FormatType == CFormattingArgs::eSAM)
    //     {
    //         x_PrintStructuredReport(results, queries);
    //         return;
    //     }
    //     ...
    //     if (results.HasWarnings()) {
    //         ERR_POST(Warning << results.GetWarningStrings());
    //     }
    // ```
    // The warnings go with the first format (outfmt 0 here); outfmt 5, 8, 9 and 12 write no query
    // result warnings, so their CLI runs equal `run_local` with that format alone (outfmt 15 and 16
    // return the same way, `x_PrintStructuredReport`).
    for (index, outfmt) in FORMATS.into_iter().enumerate() {
        let (cli_output, cli_stderr) = inputs.cli(outfmt, &[]);
        assert_eq!(run.reports[index].bytes, cli_output, "outfmt {outfmt}");
        if matches!(
            outfmt,
            "5" | "8" | "9" | "11" | "12" | "15" | "16" | "17" | "17 SQ SR"
        ) {
            let alone = inputs.collect(&[outfmt], &[]);
            assert_eq!(
                String::from_utf8_lossy(&alone.diagnostics),
                "Warning: [blastn] Subject_2 s2 no letters: Subject sequence contains no data\n"
            );
            assert_eq!(alone.diagnostics, cli_stderr, "outfmt {outfmt}");
        } else {
            assert_eq!(run.diagnostics, cli_stderr, "outfmt {outfmt}");
        }
    }
}
