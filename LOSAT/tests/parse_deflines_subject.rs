//! LOSAT's refusals of `-parse_deflines` with `-subject` input whose ids NCBI BLAST+ 2.17.0's
//! object manager would show as another sequence (`dropin-survey-20261029/DI7/DECISIONS.md` 12),
//! on a whole run of the binary of each search program.
//!
//! NCBI keeps the queries and the subjects in one scope and finds a sequence again by its id: of
//! two subjects with one id it shows the first for both, a subject with the id of a query is shown
//! with the query's title, and an accession with a version and the same accession without one are
//! one sequence (probe `DI7/ws/pd_sub.py`, 2026-11-01). A record given as both query and subject is
//! searched, except where NCBI stops (BLASTX, TBLASTN) or draws its masks at the other record's
//! places (outfmt 0; probe `DI7/ws/pd_self.py`, DECISIONS 15, 19). The reports that NCBI writes as
//! LOSAT does are the `pd.*` rows of `tests/outfmt0_manifest.tsv`.

use std::path::PathBuf;
use std::process::Command;

mod temp_support;
use temp_support::{unique_dir, TempPath};

/// The binary with no `LOSAT_*` variable that chooses a code path (`LOSAT_BLASTX_WORKER_LOG`
/// stays, as in `blastn_db_search.rs`).
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// A temporary directory of its own (the tests of this file run in parallel), removed with what is
/// in it when the value is dropped.
fn temp_dir(tag: &str) -> TempPath {
    unique_dir(&format!("losat-parse-deflines-subject-{tag}"))
}

const NUCLEOTIDES: &str =
    "ATGGCTAAACGTCTGGAAGATTTTGGTCATCCGAGCACCTGGTATCAAGCTGAATGTCGTAAACTGGGTATTCCGTAA";
const PROTEIN: &str = "MAKRLEDFGHPSTWYQAECRKLGIPMAKRLEDFGHPSTWYQAECRKLGIP";

/// The five programs with the letters of their queries and subjects.
const PROGRAMS: [(&str, &str, &str); 5] = [
    ("blastn", NUCLEOTIDES, NUCLEOTIDES),
    ("blastp", PROTEIN, PROTEIN),
    ("blastx", NUCLEOTIDES, PROTEIN),
    ("tblastn", PROTEIN, NUCLEOTIDES),
    ("tblastx", NUCLEOTIDES, NUCLEOTIDES),
];

/// Runs `program -query Q -subject S -parse_deflines -outfmt OUTFMT` on records with the deflines
/// `queries` and `subjects`; returns the exit status, stdout and stderr.
fn run(
    tag: &str,
    program: &str,
    query_letters: &str,
    subject_letters: &str,
    queries: &[&str],
    subjects: &[&str],
    outfmt: &str,
) -> (i32, Vec<u8>, String) {
    run_with(
        tag,
        program,
        query_letters,
        subject_letters,
        queries,
        subjects,
        outfmt,
        &[],
    )
}

/// `run` with the options `extra`.
#[allow(clippy::too_many_arguments)]
fn run_with(
    tag: &str,
    program: &str,
    query_letters: &str,
    subject_letters: &str,
    queries: &[&str],
    subjects: &[&str],
    outfmt: &str,
    extra: &[&str],
) -> (i32, Vec<u8>, String) {
    let dir = temp_dir(tag);
    let write = |name: &str, deflines: &[&str], letters: &str| {
        let path = dir.join(format!("{program}.{name}"));
        let text: String = deflines
            .iter()
            .map(|defline| format!(">{defline}\n{letters}\n"))
            .collect();
        std::fs::write(&path, text).unwrap();
        path
    };
    let query = write("q", queries, query_letters);
    let subject = write("s", subjects, subject_letters);
    let output = losat()
        .arg(program)
        .arg("-query")
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .args(["-parse_deflines", "-outfmt", outfmt, "-num_threads", "1"])
        .args(extra)
        .output()
        .unwrap();
    (
        output.status.code().unwrap_or(-1),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// Asserts that every program refuses the deflines `queries` and `subjects` with `expected` in its
/// message, before any byte of outfmt `outfmt`. A protein accession (`NP_`) on nucleotide letters
/// reads the same as the program's molecule (DECISIONS 7), so one name serves every program.
fn refused(tag: &str, queries: &[&str], subjects: &[&str], outfmt: &str, expected: &str) {
    for (program, query_letters, subject_letters) in PROGRAMS {
        let (status, stdout, stderr) = run(
            tag,
            program,
            query_letters,
            subject_letters,
            queries,
            subjects,
            outfmt,
        );
        assert_eq!(status, 1, "{program} {tag}: {stderr}");
        assert!(stdout.is_empty(), "{program} {tag}: stdout");
        assert!(
            stderr.contains(expected) && stderr.contains("not supported by LOSAT"),
            "{program} {tag}: {stderr}"
        );
    }
}

#[test]
fn two_subjects_with_one_id_are_refused() {
    refused(
        "dups",
        &["q1 a"],
        &["dup1 x", "dup1 y"],
        "6",
        "two subject records have the Seq-id lcl|dup1",
    );
}

#[test]
fn two_subjects_with_one_accession_with_and_without_a_version_are_refused() {
    refused(
        "dupver",
        &["q1 a"],
        &["NP_000001 x", "NP_000001.1 y"],
        "6",
        "with and without a version",
    );
}

#[test]
fn a_query_and_a_subject_with_one_id_are_refused() {
    // Other titles, or other letters: another record.
    refused(
        "title",
        &["same1 a"],
        &["same1 b"],
        "6",
        "a query and a subject record have the Seq-id lcl|same1",
    );
    refused(
        "qeqs",
        &["q1 a", "same1 b"],
        &["s1 x", "lcl|same1 y"],
        "6",
        "a query and a subject record have the Seq-id lcl|same1",
    );
    // A generated id is an id like the others (`Subject_1`, `Query_1`).
    refused(
        "gen",
        &["Subject_1 a"],
        &[""],
        "6",
        "a query and a subject record have the Seq-id lcl|Subject_1",
    );
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253
// ```c++
//         formatter.PrintProlog();
// ```
// BLASTP's query-anchored formats 1 to 4 read the queries as outfmt 0 does: the refusal of a query
// with a subject's id comes with the query's batch, after the prolog that NCBI's application writes
// before its first batch, the same bytes as outfmt 0's (line RG T3a review D3).
#[test]
fn blastp_formats_1_to_4_refuse_a_query_with_a_subject_s_id_after_outfmt_0_s_prolog() {
    // One directory: the prolog names the subject file.
    let refuse = |outfmt: &str| {
        run(
            "anchored",
            "blastp",
            PROTEIN,
            PROTEIN,
            &["same1 a"],
            &["same1 b"],
            outfmt,
        )
    };
    let (status, prolog, stderr) = refuse("0");
    assert_eq!(status, 1, "{stderr}");
    assert!(prolog.starts_with(b"BLASTP 2.17.0+\n"));
    assert!(!prolog.windows(6).any(|window| window == b"Query="));
    for outfmt in ["1", "2", "3", "4"] {
        let (status, stdout, stderr) = refuse(outfmt);
        assert_eq!(status, 1, "outfmt {outfmt}: {stderr}");
        assert_eq!(stdout, prolog, "outfmt {outfmt}");
        assert!(
            stderr.contains("a query and a subject record have the Seq-id lcl|same1"),
            "outfmt {outfmt}: {stderr}"
        );
    }
}

#[test]
fn a_query_and_a_subject_with_one_accession_with_and_without_a_version_are_refused() {
    refused(
        "qver",
        &["NP_000001 a"],
        &["NP_000001.1 x"],
        "6",
        "a query and a subject record have the Seq-id",
    );
}

#[test]
fn a_subject_named_as_a_database_ordinal_is_refused() {
    refused(
        "ordinal",
        &["q1 a"],
        &["gnl|BL_ORD_ID|7 x"],
        "6",
        "the subject Seq-id gnl|BL_ORD_ID|7 of -parse_deflines names a database ordinal",
    );
}

#[test]
fn a_subject_with_a_range_in_the_archive_is_refused() {
    refused(
        "range",
        &["q1 a"],
        &["s1:20-50 x"],
        "11",
        "a subject Seq-id with a range under -parse_deflines in the BLAST archive",
    );
}

// DECISIONS 15, 19: outfmt 6 of the programs whose NCBI reports do not change; the rows of this
// report are frozen as `pd.self.*`.
#[test]
fn one_record_as_a_query_and_a_subject_is_searched() {
    for (program, query_letters, subject_letters) in [PROGRAMS[0], PROGRAMS[1], PROGRAMS[4]] {
        let (status, stdout, stderr) = run(
            "self",
            program,
            query_letters,
            subject_letters,
            &["same1 a", "q2 b"],
            &["s2 x", "same1 a"],
            "6",
        );
        assert_eq!(status, 0, "{program}: {stderr}");
        assert!(stdout.starts_with(b"same1\tsame1\t"), "{program}");
    }
}

// DECISIONS 19 (a): NCBI stops with `eBadAlignment`, one id being DNA and protein in one scope.
#[test]
fn blastx_and_tblastn_refuse_one_record_as_a_query_and_a_subject() {
    for (program, query_letters, subject_letters) in [PROGRAMS[2], PROGRAMS[3]] {
        let (status, stdout, stderr) = run(
            "selfmol",
            program,
            query_letters,
            subject_letters,
            &["same1 a"],
            &["same1 a"],
            "6",
        );
        assert_eq!(status, 1, "{program}: {stderr}");
        assert!(stdout.is_empty(), "{program}");
        assert!(
            stderr.contains("a query and a subject record"),
            "{program}: {stderr}"
        );
    }
}

// NCBI writes the messages of both records (each title warning after the warning of its data line)
// before it stops (NCBI BLAST+ 2.17.0, these inputs, 2026-11-04); LOSAT refuses the query when it
// closes the record, as NCBI's reader would, so its messages are NCBI's first lines.
#[test]
fn blastx_and_tblastn_refuse_one_id_after_the_records_messages() {
    let hyphen = |letters: &str| format!("{}-{}", &letters[..21], &letters[21..]);
    let (nucleotides, protein) = (hyphen(NUCLEOTIDES), hyphen(PROTEIN));
    let nucleotide_title = "same1 x ACGTACGTACGTACGTACGTACGT";
    let protein_title = "same1 t ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN";
    let hyphens = "CFastaReader: Hyphens are invalid and will be ignored around line 2\n";
    let nucleotide_warning = "FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?\n";
    let protein_warning = "FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?\n";
    for (program, query, subject) in [
        (
            "blastx",
            (&nucleotides, nucleotide_title, nucleotide_warning),
            (&protein, protein_title, protein_warning),
        ),
        (
            "tblastn",
            (&protein, protein_title, protein_warning),
            (&nucleotides, nucleotide_title, nucleotide_warning),
        ),
    ] {
        let (status, stdout, stderr) = run(
            "messages",
            program,
            query.0,
            subject.0,
            &[query.1],
            &[subject.1],
            "6",
        );
        assert_eq!(status, 1, "{program}: {stderr}");
        assert!(stdout.is_empty(), "{program}");
        let messages = format!("{hyphens}{}{hyphens}{}", subject.2, query.2);
        assert!(
            stderr.starts_with(&messages)
                && stderr[messages.len()..]
                    .starts_with("Error: a query and a subject record have the Seq-id lcl|same1"),
            "{program}: {stderr}"
        );
    }
}

// DECISIONS 19 (b), (e): outfmt 0 is refused by the record's masks, not by the option: a run of one
// letter for the default DUST and SEG, lower case in either copy for `-lcase_masking`. A record
// without masks, and the other formats, are searched.
#[test]
fn outfmt_0_refuses_one_record_with_masks_only() {
    const MASKED: &str = "both a query and a subject";
    let nucleotides = format!("{NUCLEOTIDES}{}{NUCLEOTIDES}", "A".repeat(64));
    let protein = format!("{PROTEIN}{}{PROTEIN}", "Q".repeat(40));
    let lower = NUCLEOTIDES.to_ascii_lowercase();
    let lower_protein = PROTEIN.to_ascii_lowercase();
    for (program, letters, off, lcase, lower_letters) in [
        ("blastn", &nucleotides, "-dust", NUCLEOTIDES, &lower),
        ("tblastx", &nucleotides, "-seg", NUCLEOTIDES, &lower),
        ("blastp", &protein, "-seg", PROTEIN, &lower_protein),
    ] {
        let run_it = |outfmt: &str, extra: &[&str]| {
            run_with(
                "selfmask",
                program,
                letters,
                letters,
                &["same1 a"],
                &["same1 a"],
                outfmt,
                extra,
            )
        };
        let (status, _, stderr) = run_it("0", &[off, "yes"]);
        assert_eq!(status, 1, "{program}: {stderr}");
        assert!(stderr.contains(MASKED), "{program}: {stderr}");
        let (status, _, stderr) = run_it("0", &[off, "no"]);
        assert_eq!(status, 0, "{program}: {stderr}");
        let (status, _, stderr) = run_it("6", &[off, "yes"]);
        assert_eq!(status, 0, "{program}: {stderr}");
        // Lower case in one copy only: one record, whose masks are the lower case of either.
        for (query, subject, masked) in
            [(lcase, lower_letters.as_str(), true), (lcase, lcase, false)]
        {
            for option in [true, false] {
                let with_option = [off, "no", "-lcase_masking"];
                let without_option = [off, "no"];
                let extra: &[&str] = if option {
                    &with_option
                } else {
                    &without_option
                };
                let (status, _, stderr) = run_with(
                    "selflc",
                    program,
                    query,
                    subject,
                    &["same1 a"],
                    &["same1 a"],
                    "0",
                    extra,
                );
                let refused = masked && option;
                assert_eq!(status, i32::from(refused), "{program}: {stderr}");
                assert_eq!(stderr.contains(MASKED), refused, "{program}: {stderr}");
            }
        }
    }
}

// The option reads the subject's ids; without it the same file has no conflicts.
#[test]
fn the_same_files_without_the_option_are_searched() {
    for (program, query_letters, subject_letters) in PROGRAMS {
        let dir = temp_dir("plain");
        let query = dir.join(format!("{program}.q"));
        let subject = dir.join(format!("{program}.s"));
        std::fs::write(&query, format!(">same1 a\n{query_letters}\n")).unwrap();
        std::fs::write(
            &subject,
            format!(">same1 a\n{subject_letters}\n>same1 b\n{subject_letters}\n"),
        )
        .unwrap();
        let output = losat()
            .arg(program)
            .arg("-query")
            .arg(&query)
            .arg("-subject")
            .arg(&subject)
            .args(["-outfmt", "6", "-num_threads", "1"])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{program}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

/// NCBI BLAST+ 2.17.0's error archive (outfmt 11) of a `-subject` file that stops at `line`
/// (probe `ncbi-parity-20261019/SM/subnodef`, the `*.s.*.pd.f11` runs of the five programs).
fn no_defline_archive(line: u64) -> String {
    format!(
        "Blast4-archive ::= {{\n  request {{\n    body get-request-info {{\n      request-id \"Error\"\n    }}\n  }},\n  results {{\n  }},\n  messages {{\n    {{\n      code 2,\n      message \"BLAST query/options error: Query input doesn't start with a\n defline or comment, line {line}\"\n    }},\n    {{\n      code 2,\n      message \"Please refer to the BLAST+ user manual.\"\n    }}\n  }}\n}}\n"
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:209-214 and
// c++/src/app/blast/blast_app_util.hpp:172-176: the subjects' `GetAllSeqs` turns the reader's
// `eNoDefline` into a `CInputException` (the NUL of `ends` stays in the text on standard error,
// the archive drops it); the queries keep the reader's `BLAST query error` (blast_app_util.hpp:
// 181-184). The subjects are read first. Expected bytes: NCBI 2.17.0 n1 (probe `SM/subnodef`).
#[test]
fn a_subject_without_a_first_defline_ends_with_the_input_error_of_get_all_seqs() {
    for (program, query_letters, subject_letters) in PROGRAMS {
        let dir = temp_dir("nodef");
        let clean_query = dir.join(format!("{program}.q"));
        let clean_subject = dir.join(format!("{program}.s"));
        std::fs::write(&clean_query, format!(">q1 a\n{query_letters}\n")).unwrap();
        std::fs::write(&clean_subject, format!(">s1 b\n{subject_letters}\n")).unwrap();
        let bad = |name: &str, text: String| {
            let path = dir.join(format!("{program}.{name}"));
            std::fs::write(&path, text).unwrap();
            path
        };
        let subjects = [
            (
                bad(
                    "s_data",
                    format!(";c\n\n{subject_letters}\n>s2\n{subject_letters}\n"),
                ),
                3,
            ),
            (bad("s_gap", format!(">?5\n{subject_letters}\n")), 1),
        ];
        let bad_query = bad("q_data", format!("{query_letters}\n"));
        let search = |query: &std::path::Path, subject: &std::path::Path, outfmt: &str, threads| {
            let output = losat()
                .arg(program)
                .arg("-query")
                .arg(query)
                .arg("-subject")
                .arg(subject)
                .args([
                    "-parse_deflines",
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    threads,
                ])
                .output()
                .unwrap();
            (
                output.status.code().unwrap_or(-1),
                String::from_utf8_lossy(&output.stdout).into_owned(),
                String::from_utf8_lossy(&output.stderr).into_owned(),
            )
        };
        for (subject, line) in &subjects {
            let message = format!(
                "BLAST query/options error: Query input doesn't start with a defline or comment, \
                 line {line}\0\nPlease refer to the BLAST+ user manual.\n"
            );
            for threads in ["1", "4", "8", "16"] {
                for query in [&clean_query, &bad_query] {
                    let got = search(query, subject, "6", threads);
                    assert_eq!(
                        got,
                        (1, String::new(), message.clone()),
                        "{program} n{threads}"
                    );
                }
                let got = search(&clean_query, subject, "11", threads);
                assert_eq!(
                    got,
                    (1, no_defline_archive(*line), message.clone()),
                    "{program} n{threads} outfmt 11"
                );
            }
        }
        let got = search(&bad_query, &clean_subject, "6", "1");
        assert_eq!(
            got,
            (
                1,
                String::new(),
                "BLAST query error: CFastaReader: Input doesn't start with a defline or comment \
                 around line 1\n"
                    .to_string()
            ),
            "{program} query"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:327-332
// ```c++
//     } CATCH_ALL(status)
//
//     if(!m_Bah.GetMessages().empty()) {
//     	const CArgs & a = GetArgs();
//     	PrintErrorArchive(a, m_Bah.GetMessages());
//     }
// ```
// The `Run` of each of the five programs ends so. `CATCH_ALL` posts the error, which the handler
// writes to standard error at once (blast_aux.cpp:1155-1159), before `PrintErrorArchive` writes the
// error archive: in one stream of both (`2>&1`) the error comes between the archives written
// before it and the error archive (NCBI BLAST+ 2.17.0 n1, probe `ncbi-parity-20261019/SM/archorder`,
// every error of the five programs). Two errors: the reader's after the first batch's archive
// (BATCH_SIZE=1), and a `-query` file that does not open, before the search.
#[test]
fn the_error_comes_before_its_error_archive_in_one_stream() {
    const ARCHIVE: &[u8] = b"Blast4-archive ::= {";
    for (program, query_letters, subject_letters) in PROGRAMS {
        let dir = temp_dir("order");
        let subject = dir.join(format!("{program}.s"));
        std::fs::write(&subject, format!(">s1 b\n{subject_letters}\n")).unwrap();
        let bad_gi = dir.join(format!("{program}.q"));
        std::fs::write(
            &bad_gi,
            format!(">q1 a\n{query_letters}\n>gi|abc b\n{query_letters}\n"),
        )
        .unwrap();
        let missing = dir.join(format!("{program}.missing"));
        let cases: [(&PathBuf, &str, usize); 2] = [
            (
                &bad_gi,
                "BLAST query error: Could not construct seq-id from 'gi|abc'\n",
                1,
            ),
            (
                &missing,
                "Command line argument error: Argument \"query\". File is not accessible:",
                0,
            ),
        ];
        for (query, message, batches) in cases {
            for threads in ["1", "4"] {
                let search = || {
                    let mut command = losat();
                    command
                        .arg(program)
                        .arg("-query")
                        .arg(query)
                        .arg("-subject")
                        .arg(&subject)
                        .args(["-parse_deflines", "-outfmt", "11", "-num_threads", threads])
                        .env("BATCH_SIZE", "1")
                        .env_remove("ARCHIVE_FORMAT");
                    command
                };
                let apart = search().output().unwrap();
                // One file for both streams: the two descriptors share its offset, as `2>&1`.
                let both = dir.join(format!("{program}.both"));
                let file = std::fs::File::create(&both).unwrap();
                let status = search()
                    .stdout(file.try_clone().unwrap())
                    .stderr(file)
                    .status()
                    .unwrap();
                let merged = std::fs::read(&both).unwrap();
                let what = format!("{program} {} n{threads}", query.display());
                assert_eq!(apart.status.code(), Some(1), "{what}");
                assert_eq!(status.code(), Some(1), "{what}");
                assert!(apart.stderr.starts_with(message.as_bytes()), "{what}");
                let starts: Vec<usize> = (0..apart.stdout.len())
                    .filter(|&at| apart.stdout[at..].starts_with(ARCHIVE))
                    .collect();
                assert_eq!(starts.len(), batches + 1, "{what}");
                let error_archive = *starts.last().unwrap();
                let expected = [
                    &apart.stdout[..error_archive],
                    &apart.stderr[..],
                    &apart.stdout[error_archive..],
                ]
                .concat();
                assert_eq!(
                    String::from_utf8_lossy(&merged),
                    String::from_utf8_lossy(&expected),
                    "{what}"
                );
            }
        }
    }
}

/// BLASTX reads its query file once before the search with outfmt 5, 11 and 13 to 16, to refuse a
/// title its writers cannot write before any byte (`blast_xml::prescan_queries`,
/// `archive::prescan_query_titles`). Under `-parse_deflines` that read stops where the search's
/// does, so a reader error that NCBI meets before it writes a title comes first: a first line
/// that is not a defline (any format, NCBI 2.17.0, probe `SM/prescan/probe.py`), or a Seq-id NCBI
/// cannot make in the batch of the record with the title. With `BATCH_SIZE=1` that record is its
/// own batch, which NCBI writes before it reads the bad Seq-id: LOSAT refuses its title before any
/// byte (NCBI stops with an exception).
#[test]
fn blastx_reads_its_query_before_the_search_as_the_search_does() {
    let dir = temp_dir("prescan");
    let write = |name: &str, text: String| {
        let path = dir.join(format!("blastx.{name}"));
        std::fs::write(&path, text).unwrap();
        path
    };
    let subject = write("s", format!(">s1 b\n{PROTEIN}\n"));
    let later = format!(">r1 caf\u{e9} x\n{NUCLEOTIDES}\n");
    let no_defline = [
        (write("q_data", format!("{NUCLEOTIDES}\n{later}")), 1),
        (write("q_gap", format!(">?5\n{NUCLEOTIDES}\n{later}")), 1),
        (write("q_cdata", format!(";c\n\n{NUCLEOTIDES}\n{later}")), 3),
        (
            write("q_cgap", format!(";c\n\n#x\n>?5\n{NUCLEOTIDES}\n{later}")),
            4,
        ),
    ];
    let bad_gi = write(
        "q_badgi",
        format!(">r0 caf\u{e9} x\n{NUCLEOTIDES}\n>gi|abc plain\n{NUCLEOTIDES}\n"),
    );
    // A query with the id of the subject (LOSAT's refusal, which the search makes as it reads the
    // record) before a record with a title the XML cannot write: the read before the search stops
    // where the search does, so the search's refusal comes first.
    let same_id = write(
        "q_same",
        format!(">s1 a\n{NUCLEOTIDES}\n>r1 caf\u{e9} x\n{NUCLEOTIDES}\n"),
    );
    let search = |query: &std::path::Path, outfmt: &str, threads: &str, batch: Option<&str>| {
        let mut command = losat();
        command
            .arg("blastx")
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(&subject)
            .args([
                "-parse_deflines",
                "-outfmt",
                outfmt,
                "-num_threads",
                threads,
            ])
            .env_remove("BATCH_SIZE");
        if let Some(batch) = batch {
            command.env("BATCH_SIZE", batch);
        }
        let output = command.output().unwrap();
        (
            output.status.code().unwrap_or(-1),
            output.stdout,
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )
    };
    for threads in ["1", "4"] {
        for (query, line) in &no_defline {
            let message = format!(
                "BLAST query error: CFastaReader: Input doesn't start with a defline or comment \
                 around line {line}\n"
            );
            for outfmt in ["5", "11", "15", "16"] {
                let (status, stdout, stderr) = search(query, outfmt, threads, None);
                assert_eq!(
                    (status, stderr.as_str()),
                    (1, message.as_str()),
                    "line {line} outfmt {outfmt} n{threads}"
                );
                if outfmt == "5" {
                    assert!(stdout.is_empty(), "line {line} n{threads}");
                }
            }
        }
        for outfmt in ["5", "11", "16"] {
            let (status, stdout, stderr) = search(&bad_gi, outfmt, threads, None);
            assert_eq!(
                (status, stderr.as_str()),
                (
                    1,
                    "BLAST query error: Could not construct seq-id from 'gi|abc'\n"
                ),
                "gi|abc outfmt {outfmt} n{threads}"
            );
            if outfmt == "5" {
                assert!(stdout.is_empty(), "gi|abc n{threads}");
            }
        }
        let (status, stdout, stderr) = search(&same_id, "5", threads, None);
        assert_eq!(status, 1, "s1 n{threads}: {stderr}");
        assert!(stdout.is_empty(), "s1 n{threads}");
        assert!(
            stderr.contains("a query and a subject record have the Seq-id lcl|s1"),
            "s1 n{threads}: {stderr}"
        );
        let (status, stdout, stderr) = search(&bad_gi, "5", threads, Some("1"));
        assert_eq!(status, 1, "gi|abc BATCH_SIZE=1 n{threads}: {stderr}");
        assert!(stdout.is_empty(), "gi|abc BATCH_SIZE=1 n{threads}");
        assert!(
            stderr.contains("a query title with a character outside the US-ASCII range")
                && stderr.contains("not supported by LOSAT"),
            "gi|abc BATCH_SIZE=1 n{threads}: {stderr}"
        );
    }
}
