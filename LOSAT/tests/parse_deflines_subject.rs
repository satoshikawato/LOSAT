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

use std::path::{Path, PathBuf};
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
    run_in(
        &dir,
        program,
        query_letters,
        subject_letters,
        queries,
        subjects,
        outfmt,
        extra,
    )
}

/// `run_with` in the directory `dir` (the same files for every call).
#[allow(clippy::too_many_arguments)]
fn run_in(
    dir: &std::path::Path,
    program: &str,
    query_letters: &str,
    subject_letters: &str,
    queries: &[&str],
    subjects: &[&str],
    outfmt: &str,
    extra: &[&str],
) -> (i32, Vec<u8>, String) {
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
/// reads identically across programs, so one name serves every program.
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
// BLASTP's query-anchored formats 1 to 4 read the queries as outfmt 0 does: the refusal of a query
// with a subject's id comes with the query's batch, after the prolog that NCBI's application writes
// before its first batch, the same bytes as outfmt 0's (line RG T3a review D3). The reference bytes
// are LOSAT's own outfmt 0 run, which is LOSAT's refusal too, not bytes of NCBI's.
#[test]
fn blastp_formats_1_to_4_refuse_a_query_with_a_subject_s_id_after_outfmt_0_s_prolog() {
    // One directory: the prolog names the subject file.
    let dir = temp_dir("anchored");
    let refuse = |outfmt: &str| {
        run_in(
            &dir,
            "blastp",
            PROTEIN,
            PROTEIN,
            &["same1 a"],
            &["same1 b"],
            outfmt,
            &[],
        )
    };
    let (status, prolog, stderr) = refuse("0");
    assert_eq!(status, 1, "{stderr}");
    assert!(
        prolog.starts_with(concat!("BLASTP (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n").as_bytes())
    );
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
        "Blast4-archive ::= {{\n  request {{\n    body get-request-info {{\n      request-id \"Error\"\n    }}\n  }},\n  results {{\n  }},\n  messages {{\n    {{\n      code 2,\n      message \"LOSAT query/options error: Query input doesn't start with a\n defline or comment, line {line}\"\n    }},\n    {{\n      code 2,\n      message \"Please refer to the LOSAT user manual.\"\n    }}\n  }}\n}}\n"
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:209-214 and
// c++/src/app/blast/blast_app_util.hpp:172-176: the subjects' `GetAllSeqs` turns the reader's
// `eNoDefline` into a `CInputException` (the NUL of `ends` stays in the text on standard error,
// the archive drops it); the queries keep the reader's `LOSAT query error` (NCBI's `BLAST query error`) (blast_app_util.hpp:
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
                "LOSAT query/options error: Query input doesn't start with a defline or comment, \
                 line {line}\0\nPlease refer to the LOSAT user manual.\n"
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
                "LOSAT query error: CFastaReader: Input doesn't start with a defline or comment \
                 around line 1\n"
                    .to_string()
            ),
            "{program} query"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:327-332
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
                "LOSAT query error: Could not construct seq-id from 'gi|abc'\n",
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

/// BLASTX reads its query file once before the search with outfmt 5 and 13 to 16, to refuse a
/// title of punctuation alone, on which NCBI crashes, before any byte (`blast_xml::prescan_queries`;
/// outfmt 11 reads the query once, as the search does, so a reader error comes first there too).
/// Under `-parse_deflines` that read stops where the search's does, so a reader error that NCBI
/// meets before it writes a title comes first: a first line that is not a defline (any format,
/// NCBI 2.17.0, probe `SM/prescan/probe.py`), or a Seq-id NCBI cannot make in the batch of the
/// record with the title. With `BATCH_SIZE=1` that record is its own batch, which NCBI writes
/// before it reads the bad Seq-id: its title (not US-ASCII) fails the first header, so NCBI writes
/// nothing and stops with its exception (rc 255), and LOSAT writes the same.
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
                "LOSAT query error: CFastaReader: Input doesn't start with a defline or comment \
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
                    "LOSAT query error: Could not construct seq-id from 'gi|abc'\n"
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
        // The first record is its own batch, which NCBI writes before it reads `>gi|abc`: its title
        // fails the header (`BlastOutput_query-def`), nothing is written, and NCBI stops with its
        // exception (`tests/fixtures/outfmt16/prescan_badgi_blastx_5.*`, NCBI BLAST+ 2.17.0 on the
        // same two files).
        let (status, stdout, stderr) = search(&bad_gi, "5", threads, Some("1"));
        let frozen = |extension: &str| {
            // NCBI's bytes in LOSAT's words (`losat::wording`, AGENTS.md rule 8).
            let bytes = std::fs::read(
                PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/fixtures/outfmt16")
                    .join(format!("prescan_badgi_blastx_5.{extension}")),
            )
            .unwrap();
            losat::wording::losat_output_bytes(&bytes).into_owned()
        };
        let code = String::from_utf8(frozen("rc")).unwrap();
        assert_eq!(
            status,
            code.trim().parse::<i32>().unwrap(),
            "gi|abc BATCH_SIZE=1 n{threads}: {stderr}"
        );
        assert_eq!(stdout, frozen("out"), "gi|abc BATCH_SIZE=1 n{threads}");
        assert_eq!(
            stderr.as_bytes(),
            frozen("err"),
            "gi|abc BATCH_SIZE=1 n{threads}"
        );
    }
}

// BLASTP's query-anchored formats refuse the section that NCBI BLAST+ aborts on (selfseg3, the
// rejection row `db.of1.bp.selfseg3.reject`) in the search that reports the query, before its
// first byte (line RG T3a review 2 D3): a run of many queries takes the block search, which refuses
// in the query's block, after the reports of the blocks before (NCBI's own bytes before its abort
// depend on its buffers, so no prefix is NCBI's); the search of all queries refuses before any byte
// of a report, after the prolog that BLASTP writes before its first batch (line RG item T4b-f).
// `BATCH_SIZE=1` and `LOSAT_X_BPQBLOCK=1` make each query a block of its own.
#[test]
fn blastp_query_anchored_abort_is_refused_in_the_block_that_reports_it() {
    let of1 = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/of1");
    let dir = temp_dir("anchored-blocks");
    let query = dir.join("q.fa");
    let mut queries = std::fs::read(of1.join("qa.fa")).unwrap();
    queries.extend(std::fs::read(of1.join("selfseg3.fa")).unwrap());
    std::fs::write(&query, queries).unwrap();
    let run_query = |query: &Path, blocks: &str, threads: &str| {
        let output = losat()
            .arg("blastp")
            .arg("-query")
            .arg(query)
            .arg("-subject")
            .arg(of1.join("selfseg3.fa"))
            .args(["-parse_deflines", "-seg", "yes", "-outfmt", "1"])
            .args(["-num_threads", threads])
            .env("BATCH_SIZE", "1")
            .env("LOSAT_X_BPQBLOCK", blocks)
            .output()
            .unwrap();
        (
            output.status.code().unwrap_or(-1),
            String::from_utf8(output.stdout).unwrap(),
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )
    };
    let run = |blocks: &str, threads: &str| run_query(&query, blocks, threads);
    // The prolog: the bytes before the first query's report of the run of the first query alone,
    // which the search writes to its end.
    let (status, alone, stderr) = run_query(&of1.join("qa.fa"), "0", "1");
    assert_eq!(status, 0, "the first query alone: {stderr}");
    let prolog = &alone[..alone.find("\n\nQuery=").expect("the first query's report")];
    assert!(
        prolog.starts_with(concat!("BLASTP (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n")),
        "{alone}"
    );
    let refusal =
        "the record lcl|selfseg3 is both a query and a subject under -parse_deflines and \
                   NCBI BLAST+ aborts on its query-anchored alignments";
    for threads in ["1", "4"] {
        let (status, stdout, stderr) = run("0", threads);
        assert_eq!(status, 1, "n{threads}: {stderr}");
        assert_eq!(stdout, prolog, "n{threads}: exactly the prolog");
        assert!(stderr.contains(refusal), "n{threads}: {stderr}");
        let (status, stdout, stderr) = run("1", threads);
        assert_eq!(status, 1, "blocks n{threads}: {stderr}");
        assert!(stderr.contains(refusal), "blocks n{threads}: {stderr}");
        assert!(
            stdout
                .strip_prefix(prolog)
                .is_some_and(|reports| reports.starts_with("\n\nQuery= qa query protein alpha\n"))
                && !stdout.contains("Query= selfseg3"),
            "blocks n{threads}: {stdout}"
        );
    }
}

/// A run of `losat ARGS` in the fixture directory `tests/db_search/of1`: the exit status, the
/// standard output and the standard error.
fn of1_run(args: &[&str]) -> (i32, String, String) {
    let of1 = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/of1");
    let output = losat().current_dir(&of1).args(args).output().unwrap();
    (
        output.status.code().unwrap_or(-1),
        String::from_utf8(output.stdout).unwrap(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

// Line RG T3b (PD-LOSAT-NCBI-DEFECTS 1.23, 1.28): `bself.fa` is a record whose second copy of its first
// 80 letters (5 letters inserted in the middle) aligns best and whose ungapped third copy holds the
// longest segment. With the record as the query (`-query_loc 1-80`) and the subject (`-subject_loc
// 81-325`) under `-parse_deflines`, NCBI BLAST+ 2.17.0's merge of formats 1 to 4 leaves the record's
// subject rows out and aborts (`double free or corruption (out)`, exit 134; probes `bn.t3b.bself.f1` to
// `f4`); its outfmt 0 is a report. BLASTN writes its prolog before the search, as NCBI does, so the
// refusal follows the prolog.
#[test]
fn blastn_query_anchored_abort_is_refused_after_the_prolog() {
    let refusal = "the record lcl|bself is both a query and a subject under -parse_deflines and \
                   NCBI BLAST+ aborts on its query-anchored alignments";
    for outfmt in ["0", "1", "4"] {
        for threads in ["1", "4"] {
            let (status, stdout, stderr) = of1_run(&[
                "blastn",
                "-task",
                "blastn",
                "-query",
                "bself.fa",
                "-subject",
                "bself.fa",
                "-parse_deflines",
                "-query_loc",
                "1-80",
                "-subject_loc",
                "81-325",
                "-outfmt",
                outfmt,
                "-num_threads",
                threads,
            ]);
            if outfmt == "0" {
                assert_eq!(status, 0, "outfmt 0 n{threads}: {stderr}");
                assert!(stdout.contains("Query= bself"), "n{threads}: {stdout}");
                continue;
            }
            assert_eq!(status, 1, "outfmt {outfmt} n{threads}: {stderr}");
            assert!(
                stderr.contains(refusal),
                "outfmt {outfmt} n{threads}: {stderr}"
            );
            assert!(
                stdout.starts_with(concat!("BLASTN (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n"))
                    && !stdout.contains("Query="),
                "outfmt {outfmt} n{threads}: {stdout}"
            );
        }
    }
}

// Line RG T3b (PD-LOSAT-NCBI-DEFECTS 1.28): a database sequence whose Seq-align id (`bnqid`:
// `lcl|Query_1`) or display id (`bnqgi`: `gi|777|lcl|Query_1`) is the id NCBI BLAST+ gives the first
// query read without `-parse_deflines`: NCBI's formats 1 to 4 show that query's title in its
// description line (and its letters in the rows of the Seq-align id; probes `bn.t3b.d.qid.f1`,
// `bn.t3b.d.qgi.f1`). Line RG item X7-1 (PD-LOSAT-NCBI-DEFECTS 1.127): so does its outfmt 0 (NCBI
// BLAST+ 2.17.0, exit 0: `Query_1 bnq query` in the description line and the alignment's title,
// identities of 299 of 299 where the sequence's own letters give 274 of 299; probe `bn.qid.f0` of
// `x71/probe2` and `x71/fix2/of0`; rows `db.x7.bn.qid.f0.reject`, `db.x7.bn.qgi.f0.reject` hold
// NCBI's bytes before the first query's report). BLASTN refuses them all after its prolog.
#[test]
fn blastn_query_anchored_refuses_a_database_sequence_named_as_a_query() {
    for (db, refusal) in [
        (
            "bnqid",
            "database sequence 0 has the Seq-id lcl|Query_1, the id of a query (generated, or the \
             PSSM's) in the query batch that reports it",
        ),
        (
            "bnqgi",
            "database sequence 0 has the display id lcl|Query_1 (the id of its description line)",
        ),
    ] {
        for outfmt in ["0", "1", "3"] {
            let (status, stdout, stderr) = of1_run(&[
                "blastn", "-task", "blastn", "-query", "bnq.fa", "-db", db, "-outfmt", outfmt,
            ]);
            assert_eq!(status, 1, "{db} outfmt {outfmt}: {stderr}");
            assert!(stderr.contains(refusal), "{db} outfmt {outfmt}: {stderr}");
            assert!(
                stdout.starts_with(concat!("BLASTN (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n"))
                    && !stdout.contains("Query="),
                "{db} outfmt {outfmt}: {stdout}"
            );
        }
    }
}

// Line RG T3c (PD-LOSAT-NCBI-DEFECTS 1.39): the same database sequences for BLASTX (`bxqid`, `bxqgi`:
// proteins) and TBLASTX (`txqid`, `txqgi`). NCBI's merge of the Seq-align id throws (exit 255, probes
// `bx.t3c.d.qid.f1`, `tx.t3c.d.qid.f1`), and for the display id its description line shows the query's
// title (`tx.t3c.d.qgi.f1`). Line RG item X7-1 (PD-LOSAT-NCBI-DEFECTS 1.127): their outfmt 0 shows
// the query in the sequence's place (NCBI BLAST+ 2.17.0, exit 0: BLASTX `Query_1 bxq coding query
// in six frames`, `Length=1229` for the protein of 125 letters; TBLASTX the query's title; probes
// `bx.qid.f0`, `tx.qid.f0` of `x71/probe2` and `x71/fix2/of0`; rows `db.x7.bx.qid.f0.reject`,
// `db.x7.tx.qid.f0.reject` and their `qgi` rows hold NCBI's bytes before the first query's report).
// Both programs write their prolog before the search, as NCBI does, so the refusal follows it.
#[test]
fn blastx_and_tblastx_query_anchored_refuse_a_database_sequence_named_as_a_query() {
    for (program, prefix, version) in [
        (
            "blastx",
            "bx",
            concat!("BLASTX (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n"),
        ),
        (
            "tblastx",
            "tx",
            concat!("TBLASTX (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n"),
        ),
    ] {
        let query = format!("{prefix}q.fa");
        for (db, refusal) in [
            (
                format!("{prefix}qid"),
                "database sequence 0 has the Seq-id lcl|Query_1, the id of a query (generated, or \
                 the PSSM's) in the query batch that reports it",
            ),
            (
                format!("{prefix}qgi"),
                "database sequence 0 has the display id lcl|Query_1 (the id of its description \
                 line)",
            ),
        ] {
            for outfmt in ["0", "1", "3"] {
                let (status, stdout, stderr) =
                    of1_run(&[program, "-query", &query, "-db", &db, "-outfmt", outfmt]);
                assert_eq!(status, 1, "{program} {db} outfmt {outfmt}: {stderr}");
                assert!(
                    stderr.contains(refusal),
                    "{program} {db} outfmt {outfmt}: {stderr}"
                );
                assert!(
                    stdout.starts_with(version) && !stdout.contains("Query="),
                    "{program} {db} outfmt {outfmt}: {stdout}"
                );
            }
        }
    }
}

// Line RG T3c fix 1: `txself.fa` is one record read as the query and as the subject under
// `-parse_deflines`. TBLASTX's formats 1 to 4 add their Dense-segs with `fForceTranslation`, not
// `fPreserveRows`, so NCBI's `CAlnMix` takes the query's row and the record's subject rows for one
// sequence and folds the record's HSPs with itself into the query row (probes `tx.t3c.self.f1` to
// `f4`, exit 0; `report::query_anchored::aln_mix`): LOSAT writes the report, the same at 1 and 4
// threads (its bytes are NCBI's: rows `db.of1.tx.self.*` of the fixtures).
#[test]
fn tblastx_query_anchored_writes_one_record_as_a_query_and_a_subject() {
    for outfmt in ["0", "1", "4"] {
        let mut outputs = Vec::new();
        for threads in ["1", "4"] {
            let (status, stdout, stderr) = of1_run(&[
                "tblastx",
                "-query",
                "txself.fa",
                "-subject",
                "txself.fa",
                "-parse_deflines",
                "-outfmt",
                outfmt,
                "-num_threads",
                threads,
            ]);
            assert_eq!(status, 0, "outfmt {outfmt} n{threads}: {stderr}");
            assert!(
                stdout.contains("Query= txself"),
                "outfmt {outfmt} n{threads}: {stdout}"
            );
            if outfmt != "0" {
                assert!(
                    stdout.contains("For reading frame +1 of query sequence:"),
                    "outfmt {outfmt} n{threads}: {stdout}"
                );
            }
            outputs.push(stdout);
        }
        assert_eq!(outputs[0], outputs[1], "outfmt {outfmt}: n1 and n4 differ");
    }
}

/// NCBI BLAST+ 2.17.0's standard error where its merge of a query-anchored section throws at
/// `file`, `line`, in `function`, with `message` (`CATCH_ALL`'s report, exit 255), in LOSAT's
/// words (NCBI's text through `losat::wording`, AGENTS.md rule 8).
fn merge_exception(file: &str, line: u32, function: &str, message: &str) -> String {
    losat::wording::losat_text(&format!(
        "Error: NCBI C++ Exception:\n    T0 \"/home/coremake/release_build/build/\
         PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_\
         1751374548/c++/compilers/unix/../../src/objtools/alnmgr/{file}\", line {line}: Error: \
         (CAlnException::eMergeFailure) ncbi::objects::{function}() - {message}\n\n"
    ))
    .into_owned()
}

// Line RG T3c fix 2 (PD-LOSAT-NCBI-DEFECTS 1.39, rule 2): `txrep2.fa` is one record of two copies of a
// coding piece. As the query (`-query_loc 1-300`, the first copy; or after `lcl|other1` in `txrep2b.fa`)
// and the subject (`-subject_loc 301-600`, the second) under `-parse_deflines`, the first HSP of the
// record with itself that NCBI's merge takes puts subject positions on the query's own row where the
// query has none, and `x_Merge` throws `Start iterator already exists for seq2` (probes
// `f1.rep2.qs.f1`, `f1.rep2.qs.f3`, `f1.rep2b.b1.sl.f1` of `T/t3/fix1_t3c/cases_fix1.tsv`). `txmg0.fa` is
// a tandem record `A+A` and `txmg0s.fa` that record, then a subject holding `A` (`-subject_loc
// 154-303 -seg no`): the merge joins segments of the query's row and `Build` throws `Infinite loop
// detected` (review 2 of fix 1, `mg0.sx.f1`). NCBI's application catches the exception and exits with
// 255, its stream flushed: the reports of the queries before, the query's heading and description
// table, no epilog, and the exception on standard error. LOSAT writes the same, at 1 and 4 threads and
// with `BATCH_SIZE=1` (formats 1 to 4 are written batch by batch); the whole bytes are the fixture rows
// `db.of1.tx.throw.*`. Outfmt 0 of the same search is a report.
#[test]
fn tblastx_query_anchored_ends_with_ncbi_s_exception_where_its_merge_throws() {
    let start_iterator = merge_exception(
        "alnmerger.cpp",
        790,
        "CAlnMixMerger::x_Merge",
        "CAlnMixMerger::x_Merge(): Internal error: Start iterator already exists for seq2.",
    );
    let infinite_loop = merge_exception(
        "alnsegments.cpp",
        168,
        "CAlnMixSegments::Build",
        "CAlnMixSegments::Build(): Internal error: Infinite loop detected.",
    );
    let of1 = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/of1");
    let rep2: &[&str] = &["-subject", "txrep2.fa", "-subject_loc", "301-600"];
    let mg0: &[&str] = &[
        "-subject",
        "txmg0s.fa",
        "-subject_loc",
        "154-303",
        "-seg",
        "no",
    ];
    // (query, subject, `-query_loc 1-300`, `BATCH_SIZE`, outfmt, NCBI's exception)
    let cases: [(&str, &[&str], bool, Option<&str>, &str, Option<&String>); 8] = [
        ("txrep2.fa", rep2, true, None, "0", None),
        ("txrep2.fa", rep2, true, None, "1", Some(&start_iterator)),
        ("txrep2.fa", rep2, true, None, "3", Some(&start_iterator)),
        (
            "txrep2b.fa",
            rep2,
            false,
            Some("1"),
            "1",
            Some(&start_iterator),
        ),
        ("txrep2b.fa", rep2, false, None, "2", Some(&start_iterator)),
        ("txmg0.fa", mg0, false, None, "0", None),
        ("txmg0.fa", mg0, false, None, "1", Some(&infinite_loop)),
        ("txmg0.fa", mg0, false, None, "4", Some(&infinite_loop)),
    ];
    for (query, subject, query_loc, batch_size, outfmt, exception) in cases {
        let mut outputs = Vec::new();
        for threads in ["1", "4"] {
            let mut command = losat();
            command
                .current_dir(&of1)
                .args(["tblastx", "-query", query])
                .args(subject)
                .args([
                    "-parse_deflines",
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    threads,
                ]);
            if query_loc {
                command.args(["-query_loc", "1-300"]);
            }
            if let Some(size) = batch_size {
                command.env("BATCH_SIZE", size);
            }
            let output = command.output().unwrap();
            let (status, stdout, stderr) = (
                output.status.code().unwrap_or(-1),
                String::from_utf8(output.stdout).unwrap(),
                String::from_utf8_lossy(&output.stderr).into_owned(),
            );
            let what = format!("{query} {batch_size:?} outfmt {outfmt} n{threads}");
            let Some(exception) = exception else {
                assert_eq!(status, 0, "{what}: {stderr}");
                assert!(stdout.contains("Lambda"), "{what}: {stdout}");
                continue;
            };
            assert_eq!(status, 255, "{what}: {stderr}");
            assert_eq!(&stderr, exception, "{what}");
            assert!(
                stdout.starts_with(concat!("TBLASTX (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n")),
                "{what}: {stdout}"
            );
            // The report ends with the query's description table, before its alignments.
            let last = &stdout[stdout.rfind("\n\n\nQuery= ").expect("a query's report")..];
            assert!(
                last.contains("Sequences producing significant alignments:")
                    && last.ends_with("\n\n")
                    && !last.contains("Lambda")
                    && !last.contains("For reading frame"),
                "{what}: {stdout}"
            );
            // The report of `lcl|other1` before it, also in the same batch.
            assert_eq!(
                stdout.contains("Query= other1"),
                query == "txrep2b.fa",
                "{what}: {stdout}"
            );
            outputs.push(stdout);
        }
        if outputs.len() == 2 {
            assert_eq!(
                outputs[0], outputs[1],
                "{query} outfmt {outfmt}: n1 and n4 differ"
            );
        }
    }
}
