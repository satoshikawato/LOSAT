//! A first line of the query or of the `-subject` file that NCBI BLAST+ 2.17.0 tries as a Seq-id
//! when its data loaders are on (`CBlastInputReader::ReadOneSeq`, blast_fasta_input.cpp:130-162),
//! on a whole run of the binary of each search program.
//!
//! NCBI looks such a line up as a sequence (GenBank over the network, or a BLAST database), which
//! LOSAT does not do: the five programs refuse it with one text (`fasta_reader::seq_id`), the
//! refusal of `-parse_deflines`'s stop included, since NCBI tries the Seq-id first. With the data
//! loaders off (an empty `NCBI_CONFIG__BLAST__DATA_LOADERS`, or `DATA_LOADERS = none` in the
//! registry) NCBI reads the line as FASTA data, and so does LOSAT: for BLASTX the bytes below are
//! NCBI's (probe `SM/seqid/probe.py` of line NP, 2026-11-07), for the other programs the
//! `seqid0_loaders_env_empty` rows of `tests/fasta_input_fixtures.py`. NCBI is never run here with
//! such a line and its data loaders on.

use std::path::{Path, PathBuf};
use std::process::Command;

const N1: &str = "ATGCTGTGTCGCATATTTAAGGACCACTTGCATCAGGAGTGGCGCGCACCCAATATAATT";
const N2: &str = "CTGCGGTGCCACGACGTCCTCCCCTCACGGCGTTGTGCCAACGACTCGGACATGTGGACA";
const P1: &str = "MLCRIFKDHLHQEWRAPNIILRCHDVLPSRRCANDSDMWTYDMFNCDELHKGDKAHEATQ";
const P2: &str = "WKVLCSRYHPQEPLMMRHDTIDWCKYVAGLEPNNRFATFHSIQMCYHRAFNERPTWKKHR";

/// The five programs with whether their queries and subjects are protein.
const PROGRAMS: [(&str, bool, bool); 5] = [
    ("blastn", false, false),
    ("blastp", true, true),
    ("blastx", false, true),
    ("tblastn", true, false),
    ("tblastx", false, false),
];

/// A temporary directory of this test.
fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-seq-id-first-line-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Two records of the role (`q` or `s`), after the line `first` (none when empty).
fn records(first: &str, role: &str, protein: bool) -> String {
    let (a, b) = if protein { (P1, P2) } else { (N1, N2) };
    let first = if first.is_empty() {
        String::new()
    } else {
        format!("{first}\n")
    };
    format!("{first}>{role}1 one\n{a}\n{b}\n>{role}2 two\n{b}\n")
}

/// Runs `program -query Q -subject S EXTRA -num_threads 1` in `dir` (which may hold a `.ncbirc`),
/// with an empty HOME and no NCBI, BLAST or LOSAT variable but `env`; returns the exit status,
/// stdout and stderr.
fn run(
    dir: &Path,
    program: &str,
    query: &Path,
    subject: &Path,
    extra: &[&str],
    env: &[(&str, &str)],
) -> (i32, Vec<u8>, String) {
    let home = dir.join("home");
    std::fs::create_dir_all(&home).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name = name.to_string_lossy();
        if (name.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || name.starts_with("NCBI")
            || name.starts_with("BLAST")
            || name == "BATCH_SIZE"
        {
            command.env_remove(&*name);
        }
    }
    let output = command
        .current_dir(dir)
        .env("HOME", &home)
        .env("USERPROFILE", &home)
        .envs(env.iter().copied())
        .arg(program)
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(extra)
        .args(["-num_threads", "1"])
        .output()
        .unwrap();
    (
        output.status.code().unwrap_or(-1),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// LOSAT's refusal of the first line `line` of the `role` of `program`, as it ends the run.
fn refusal(program: &str, role: &str, line: &str) -> String {
    format!(
        "Error: the first line of the {role} ({line:?}) is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT's {} (start the {role} with a '>' defline, or set DATA_LOADERS=none in the [BLAST] section of .ncbirc)\n",
        program.to_uppercase()
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-151
// ```c++
//     virtual CRef<CSeq_entry> ReadOneSeq(ILineErrorListener * pMessageListener) {
//
//         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
//         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
//             try {
//                 CRef<CSeq_id> id(new CSeq_id(line, (CSeq_id::fParse_AnyRaw |
// 							CSeq_id::fParse_ValidLocal)));
// ...
//             } catch (const CSeqIdException& e) {
//                 if (NStr::Find(e.GetMsg(), "Malformatted ID") != NPOS) {
//                     // This is probably just plain fasta, so just
//                     // defer to CFastaReader
// ```
// Every program refuses the line in the query and in the subject with the same text and exit
// status, with and without `-parse_deflines`; the white space around the line is trimmed, as
// NCBI trims it. The refusal comes when the reader reads the line: a subject's before any byte, a
// query's after what the report writes before its first batch (outfmt 0: the prolog, no
// `Query=`; outfmt 6 and 11: nothing, and no error archive, since the refusal is LOSAT's, not an
// error of NCBI's).
#[test]
fn the_five_programs_refuse_a_seq_id_first_line_with_one_text() {
    let dir = temp_dir("refuse");
    for (program, query_protein, subject_protein) in PROGRAMS {
        let clean_q = dir.join(format!("{program}_clean_q.fa"));
        let clean_s = dir.join(format!("{program}_clean_s.fa"));
        std::fs::write(&clean_q, records("", "q", query_protein)).unwrap();
        std::fs::write(&clean_s, records("", "s", subject_protein)).unwrap();
        for (first, trimmed) in [
            ("lcl|abc", "lcl|abc"),
            ("12345", "12345"),
            ("  AB123456 ", "AB123456"),
        ] {
            for role in ["query", "subject"] {
                let path = dir.join(format!("{program}_{role}_{}.fa", trimmed.replace('|', "_")));
                let protein = if role == "query" {
                    query_protein
                } else {
                    subject_protein
                };
                std::fs::write(&path, records(first, &role[..1], protein)).unwrap();
                let (query, subject) = if role == "query" {
                    (&path, &clean_s)
                } else {
                    (&clean_q, &path)
                };
                for outfmt in ["6", "11", "0"] {
                    for parse_deflines in [false, true] {
                        let mut extra = vec!["-outfmt", outfmt];
                        if parse_deflines {
                            extra.push("-parse_deflines");
                        }
                        let what = format!("{program} {role} {first:?} {extra:?}");
                        let (status, stdout, stderr) =
                            run(&dir, program, query, subject, &extra, &[]);
                        assert_eq!(status, 1, "{what}: {stderr}");
                        let report = String::from_utf8_lossy(&stdout);
                        if outfmt == "0" && role == "query" {
                            let prolog = format!("{} 2.17.0+\n", program.to_uppercase());
                            assert!(
                                report.starts_with(&prolog) && !report.contains("Query="),
                                "{what}: {report}"
                            );
                        } else {
                            assert!(stdout.is_empty(), "{what}: {report}");
                        }
                        assert_eq!(stderr, refusal(program, role, trimmed), "{what}");
                    }
                }
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:76-92
// ```c++
// SDataLoaderConfig::x_LoadDataLoadersConfig(const CMetaRegistry::SEntry& sentry)
// {
//     _ASSERT(sentry.registry);
//
//     static const string kDataLoadersConfig("DATA_LOADERS");
//
//     if (sentry.registry->HasEntry("BLAST", kDataLoadersConfig)) {
//         const string& loaders = sentry.registry->Get("BLAST", kDataLoadersConfig);
//         if (NStr::FindNoCase(loaders, "blastdb") == NPOS) {
//             m_UseBlastDbs = false;
//         }
//         if (NStr::FindNoCase(loaders, "genbank") == NPOS) {
//             m_UseGenbank = false;
//         }
//         if (NStr::FindNoCase(loaders, "none") != NPOS) {
//             m_UseBlastDbs = false;
//             m_UseGenbank = false;
//         }
//     }
// ```
// With the data loaders off, by the environment or by a registry file, BLASTX reads the line as
// the letters of a first record without a defline (NCBI's bytes), and under `-parse_deflines` it
// stops there with NCBI's error of the query or of the subjects; the other programs give the
// same output for both ways of turning the loaders off.
#[test]
fn with_the_data_loaders_off_the_line_is_read_as_ncbi_reads_it() {
    let dir = temp_dir("off");
    let rc_dir = dir.join("rc");
    std::fs::create_dir_all(&rc_dir).unwrap();
    std::fs::write(rc_dir.join(".ncbirc"), "[BLAST]\nDATA_LOADERS = none\n").unwrap();
    let off: [(&Path, &[(&str, &str)]); 2] = [
        (&dir, &[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")]),
        (&rc_dir, &[]),
    ];
    for (program, query_protein, subject_protein) in PROGRAMS {
        let file = |name: &str, first: &str, role: &str, protein: bool| {
            let path = dir.join(format!("{program}_{name}.fa"));
            std::fs::write(&path, records(first, role, protein)).unwrap();
            path
        };
        let clean_q = file("clean_q", "", "q", query_protein);
        let clean_s = file("clean_s", "", "s", subject_protein);
        let lcl_q = file("lcl_q", "lcl|abc", "q", query_protein);
        let lcl_s = file("lcl_s", "lcl|abc", "s", subject_protein);
        for (role, query, subject) in [("query", &lcl_q, &clean_s), ("subject", &clean_q, &lcl_s)] {
            for parse_deflines in [false, true] {
                let mut extra = vec!["-outfmt", "6"];
                if parse_deflines {
                    extra.push("-parse_deflines");
                }
                let runs: Vec<_> = off
                    .iter()
                    .map(|(cwd, env)| run(cwd, program, query, subject, &extra, env))
                    .collect();
                let what = format!("{program} {role} {extra:?}");
                assert_eq!(runs[0], runs[1], "{what}");
                if program != "blastx" {
                    assert_eq!(
                        runs[0].0,
                        if parse_deflines { 1 } else { 0 },
                        "{what}: {}",
                        runs[0].2
                    );
                    continue;
                }
                let expected: (i32, &[u8], &str) = match (role, parse_deflines) {
                    ("query", false) => (
                        0,
                        b"q1\ts1\t100.000\t40\t0\t0\t1\t120\t1\t40\t1.65e-29\t87.0\nq2\ts1\t100.000\t20\t0\t0\t1\t60\t21\t40\t1.03e-13\t45.8\n",
                        "FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 1, 3-4\n",
                    ),
                    ("subject", false) => (
                        0,
                        b"q1\ts1\t100.000\t40\t0\t0\t1\t120\t1\t40\t1.70e-29\t87.0\nq2\ts1\t100.000\t20\t0\t0\t1\t60\t21\t40\t1.07e-13\t45.8\n",
                        "FASTA-Reader: Ignoring invalid residues at position(s): On line 1: 4\n",
                    ),
                    ("query", true) => (
                        1,
                        b"",
                        "BLAST query error: CFastaReader: Input doesn't start with a defline or comment around line 1\n",
                    ),
                    _ => (
                        1,
                        b"",
                        "BLAST query/options error: Query input doesn't start with a defline or comment, line 1\0\nPlease refer to the BLAST+ user manual.\n",
                    ),
                };
                let (status, stdout, stderr) = &runs[0];
                assert_eq!(
                    (*status, stdout.as_slice(), stderr.as_str()),
                    expected,
                    "{what}"
                );
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

// Only the input's first line is tried (a record ends before the next `>` line): a blank line or
// a comment before the Seq-id, a defline that starts with one, and a line of letters only (which
// `CSeq_id` does not take: Malformatted ID) are read with the data loaders on as with them off.
#[test]
fn only_a_first_line_that_ncbi_takes_as_a_seq_id_is_refused() {
    let dir = temp_dir("fasta");
    for (program, query_protein, subject_protein) in PROGRAMS {
        let clean_q = dir.join(format!("{program}_clean_q.fa"));
        let clean_s = dir.join(format!("{program}_clean_s.fa"));
        std::fs::write(&clean_q, records("", "q", query_protein)).unwrap();
        std::fs::write(&clean_s, records("", "s", subject_protein)).unwrap();
        for (tag, first) in [
            ("blank", "\nlcl|abc"),
            ("comment", ";c\nlcl|abc"),
            ("defline", ">lcl|abc x\nACGT"),
            ("letters", "ACGTACGTAC"),
        ] {
            for role in ["query", "subject"] {
                let path = dir.join(format!("{program}_{role}_{tag}.fa"));
                let protein = if role == "query" {
                    query_protein
                } else {
                    subject_protein
                };
                std::fs::write(&path, records(first, &role[..1], protein)).unwrap();
                let (query, subject) = if role == "query" {
                    (&path, &clean_s)
                } else {
                    (&clean_q, &path)
                };
                let extra = ["-outfmt", "6"];
                let on = run(&dir, program, query, subject, &extra, &[]);
                let off = run(
                    &dir,
                    program,
                    query,
                    subject,
                    &extra,
                    &[("NCBI_CONFIG__BLAST__DATA_LOADERS", "")],
                );
                assert_eq!(on, off, "{program} {role} {tag}");
                assert_eq!(on.0, 0, "{program} {role} {tag}: {}", on.2);
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI reference: ncbi-blast/c++/src/corelib/ncbiapp.cpp:1031-1044
// ```c
//     // Setup some debugging features from environment variables.
//     if ( !m_Environ->Get(DIAG_TRACE).empty() ) {
//         SetDiagTrace(eDT_Enable, eDT_Enable);
//     }
//     string post_level = m_Environ->Get(DIAG_POST_LEVEL);
// ```
// NCBI's application layer reads its environment and registry files (the data loaders among
// them) before the program's own code: every program refuses the settings that LOSAT does not
// reproduce with one text, before its run checks an option or opens a file (BLASTX reads them at
// the start of its run, the others in `main`; the command line's own parse comes first in all).
#[test]
fn the_settings_are_read_before_the_options_and_the_files() {
    let dir = temp_dir("settings");
    let missing = dir.join("missing.fa");
    let env = [("DIAG_POST_LEVEL", "Trace")];
    let mut texts = Vec::new();
    for (program, query_protein, _) in PROGRAMS {
        let query = dir.join(format!("{program}_q.fa"));
        std::fs::write(&query, records("", "q", query_protein)).unwrap();
        for extra in [&["-outfmt", "6"][..], &["-outfmt", "99"]] {
            let (status, stdout, stderr) = run(&dir, program, &query, &missing, extra, &env);
            assert_eq!(status, 1, "{program} {extra:?}: {stderr}");
            assert!(stdout.is_empty(), "{program} {extra:?}");
            assert!(
                stderr.starts_with("Error: the environment variable DIAG_POST_LEVEL")
                    && stderr.contains("is not supported by LOSAT"),
                "{program} {extra:?}: {stderr}"
            );
            texts.push(stderr);
        }
    }
    assert!(texts.iter().all(|text| *text == texts[0]), "{texts:?}");
    let _ = std::fs::remove_dir_all(&dir);
}
