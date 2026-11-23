//! BLASTP's refusals of a step where NCBI's search of a later query batch stops (line RG item X4,
//! fix round 2; PD-LOSAT-NCBI-DEFECTS 1.53, decision D13): NCBI prints each batch's reports before
//! it searches the next, so where its search of batch 2 stops with an exception, batch 1's reports
//! are printed (without the epilog). LOSAT refuses after them. `x4/chunk_second.faa`: a protein of
//! PajaWSV (92 letters, `x4/chunk_first.faa`) and the TK-TMK protein of psiblast's rows (398
//! letters); with `BATCH_SIZE=50` each is a batch, and with `CHUNK_SIZE=150` (the overlap 100) the
//! first is not split (92 / 50 makes one chunk) and the second is split into chunks of 143 letters
//! that NCBI would split again, which NCBI stops on with a null reference (`db.ps.it.rej.split`'s
//! case). `x4/chunk_empty.faa`: 9 and 31 letters of that protein; with `BATCH_SIZE=1 CHUNK_SIZE=100
//! OVERLAP_CHUNK_SIZE=95` the first is not split (9 / 5 makes one chunk) and the second makes 6
//! chunks of an evened size of 85, below the overlap: chunk 0 covers the batch (and would be split
//! again) and chunks 1 to 5 hold no query, on which NCBI's `CQuerySplitter::Split` stops first
//! (split_query_cxx.cpp:234,876-879). A negative pair (fix round 4): `outfmt0/xsb_ql.faa` (qa, 31
//! letters, and ql, qa ten times) with `BATCH_SIZE=1 CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200` (NCBI's
//! `size_t` sizes 100 apart): qa is one chunk (31 / 100), ql makes 3 chunks whose evened size wraps
//! past its 310 letters, so chunk 0 covers ql and chunks 1 and 2 hold no query: NCBI stops there
//! after qa's report (NCBI BLAST+ 2.17.0 n1: exit status 3, rows `db.x4.chunk.neg.*`). Where the wrapped
//! ranges leave gaps between chunks that all hold a query (`x4/chunk_second.faa` with `CHUNK_SIZE=-1
//! OVERLAP_CHUNK_SIZE=-200`: the 398 letters make chunks [0, 100) and [300, 398)), NCBI searches on
//! without the hits in the gap; the maintainer kept that setting refused (PD-LOSAT-NCBI-DEFECTS
//! version 1.1), and LOSAT refuses where the search of that batch's first chunk with a valid
//! context would read subjects in hand, after batch 1's report (review 3 of X4, row
//! `db.x4.chunk.neg.gap.6`; the step of TBLASTN's since PD-LOSAT-NCBI-DEFECTS 1.90, line RG item
//! X-bpgap).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276

use std::path::{Path, PathBuf};
use std::process::Command;

fn path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(name)
}

struct Run {
    stdout: String,
    stderr: String,
    code: Option<i32>,
}

/// `losat blastp -query QUERY -db tests/db_search/paja -outfmt OUTFMT -num_threads THREADS` with
/// `env` as the only variables that choose a code path or NCBI's batches and chunks.
fn search(query: &str, outfmt: &str, threads: usize, env: &[(&str, &str)]) -> Run {
    search_in(
        "blastp",
        query,
        &["-db", "tests/db_search/paja"],
        outfmt,
        threads,
        env,
    )
}

/// `search` of `program` against the subjects `subjects` (arguments, paths under the manifest
/// directory after `-db`, `-subject` and `-negative_seqidlist`).
fn search_in(
    program: &str,
    query: &str,
    subjects: &[&str],
    outfmt: &str,
    threads: usize,
    env: &[(&str, &str)],
) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        // `LOSAT_BLASTX_WORKER_LOG` is the log file of the build with `--all-features` (its worker
        // probe), not a code path.
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE"
                    | "CHUNK_SIZE"
                    | "OVERLAP_CHUNK_SIZE"
                    | "RAYON_NUM_THREADS"
                    | "BLASTDB"
            )
        {
            command.env_remove(name);
        }
    }
    let mut paths = false;
    let subjects: Vec<std::ffi::OsString> = subjects
        .iter()
        .map(|&arg| {
            let value = if paths {
                path(arg).into_os_string()
            } else {
                arg.into()
            };
            paths = matches!(arg, "-db" | "-subject" | "-negative_seqidlist");
            value
        })
        .collect();
    let output = command
        .envs(env.iter().copied())
        .arg(program)
        .arg("-query")
        .arg(path(query))
        .args(subjects)
        .args(["-outfmt", outfmt, "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    Run {
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
    }
}

#[test]
fn a_chunk_split_again_in_batch_2_is_refused_after_batch_1() {
    let env = [("BATCH_SIZE", "50"), ("CHUNK_SIZE", "150")];
    let alone = search("tests/fasta/x4/chunk_first.faa", "6", 1, &env);
    let alone7 = search("tests/fasta/x4/chunk_first.faa", "7", 1, &env);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert!(!alone.stdout.is_empty());
    let epilog = alone7
        .stdout
        .rfind("# BLAST processed")
        .expect("outfmt 7's last line");
    for threads in [1, 4] {
        let run = search("tests/fasta/x4/chunk_second.faa", "6", threads, &env);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        assert!(
            run.stderr.contains("would split a query chunk again"),
            "n{threads}: {}",
            run.stderr
        );
        let commented = search("tests/fasta/x4/chunk_second.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}");
        assert_eq!(commented.stdout, alone7.stdout[..epilog], "n{threads}");
        // In one batch (the default batch size) the refusal comes before the first byte.
        let one_batch = search(
            "tests/fasta/x4/chunk_second.faa",
            "6",
            threads,
            &[("CHUNK_SIZE", "150")],
        );
        assert_eq!(one_batch.code, Some(1), "n{threads}");
        assert_eq!(one_batch.stdout, "", "n{threads}");
    }
}

#[test]
fn a_chunk_without_a_query_in_batch_2_is_refused_before_a_split_again() {
    let env = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "100"),
        ("OVERLAP_CHUNK_SIZE", "95"),
    ];
    let alone7 = search("tests/fasta/x4/chunk_empty.faa", "7", 1, &[]);
    assert_eq!(alone7.code, Some(0), "{}", alone7.stderr);
    // Batch 1's comment lines and rows: the part of the run without the batches' epilog before
    // q31's `# Query:` line.
    let batch1 = alone7
        .stdout
        .find(concat!(
            "# BLASTP (LOSAT ",
            env!("CARGO_PKG_VERSION"),
            ")\n# Query: q31"
        ))
        .expect("q31's report");
    for threads in [1, 4] {
        let commented = search("tests/fasta/x4/chunk_empty.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}: {}", commented.stderr);
        assert!(
            commented
                .stderr
                .contains("a query chunk without a query is not supported by LOSAT's BLASTP"),
            "n{threads}: {}",
            commented.stderr
        );
        assert_eq!(commented.stdout, alone7.stdout[..batch1], "n{threads}");
    }
}

#[test]
fn a_negative_pair_that_leaves_chunks_without_a_query_in_batch_2_is_refused_after_batch_1() {
    let env = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-100"),
        ("OVERLAP_CHUNK_SIZE", "-200"),
    ];
    // qa's report in one batch of both queries without the variables (not split): its rows, and
    // its comment lines and rows before ql's `# Query:` line.
    let whole = search("tests/fasta/outfmt0/xsb_ql.faa", "6", 1, &[]);
    let whole7 = search("tests/fasta/outfmt0/xsb_ql.faa", "7", 1, &[]);
    assert_eq!(whole.code, Some(0), "{}", whole.stderr);
    let qa: String = whole
        .stdout
        .lines()
        .filter(|line| line.starts_with("qa\t"))
        .map(|line| format!("{line}\n"))
        .collect();
    assert!(!qa.is_empty());
    let batch1 = whole7
        .stdout
        .find(concat!(
            "# BLASTP (LOSAT ",
            env!("CARGO_PKG_VERSION"),
            ")\n# Query: ql"
        ))
        .expect("ql's report");
    for threads in [1, 4] {
        let run = search("tests/fasta/outfmt0/xsb_ql.faa", "6", threads, &env);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert!(
            run.stderr
                .contains("a query chunk without a query is not supported by LOSAT's BLASTP"),
            "n{threads}: {}",
            run.stderr
        );
        assert_eq!(run.stdout, qa, "n{threads}");
        let commented = search("tests/fasta/outfmt0/xsb_ql.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}: {}", commented.stderr);
        assert_eq!(commented.stdout, whole7.stdout[..batch1], "n{threads}");
    }
}

#[test]
fn a_negative_pair_that_leaves_gaps_between_its_chunks_in_batch_2_is_refused_after_batch_1() {
    let env = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-1"),
        ("OVERLAP_CHUNK_SIZE", "-200"),
    ];
    // Batch 1 (qsok2, 92 letters: one chunk with this pair) alone, without the variables.
    let alone = search("tests/fasta/x4/chunk_first.faa", "6", 1, &[]);
    let alone7 = search("tests/fasta/x4/chunk_first.faa", "7", 1, &[]);
    assert_eq!(alone.code, Some(0), "{}", alone.stderr);
    assert!(!alone.stdout.is_empty());
    let epilog = alone7
        .stdout
        .rfind("# BLAST processed")
        .expect("outfmt 7's last line");
    let message =
        "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this query batch";
    for threads in [1, 4] {
        let run = search("tests/fasta/x4/chunk_second.faa", "6", threads, &env);
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert!(run.stderr.contains(message), "n{threads}: {}", run.stderr);
        assert_eq!(run.stdout, alone.stdout, "n{threads}");
        let commented = search("tests/fasta/x4/chunk_second.faa", "7", threads, &env);
        assert_eq!(commented.code, Some(1), "n{threads}: {}", commented.stderr);
        assert!(
            commented.stderr.contains(message),
            "n{threads}: {}",
            commented.stderr
        );
        assert_eq!(commented.stdout, alone7.stdout[..epilog], "n{threads}");
        // In one batch (the default batch size) the refusal comes before the first byte.
        let one_batch = search(
            "tests/fasta/x4/chunk_second.faa",
            "6",
            threads,
            &[("CHUNK_SIZE", "-1"), ("OVERLAP_CHUNK_SIZE", "-200")],
        );
        assert_eq!(one_batch.code, Some(1), "n{threads}: {}", one_batch.stderr);
        assert_eq!(one_batch.stdout, "", "n{threads}");
    }
}

const GAP: [(&str, &str); 2] = [("CHUNK_SIZE", "-1"), ("OVERLAP_CHUNK_SIZE", "-200")];
const NO_SUBJECT: [&str; 4] = [
    "-db",
    "tests/db_search/taxfilter/tfp",
    "-negative_seqidlist",
    "tests/db_search/seqidlist/sl_tfp_all.txt",
];
const SHORT: [&str; 4] = [
    "-subject",
    "tests/db_search/q1c/avg_s.faa",
    "-searchsp",
    "1000",
];

// The gaps (PD-LOSAT-NCBI-DEFECTS 1.90, line RG item X-bpgap): refused only where the search of the
// batch's first chunk with a valid context would read subjects in hand. A database left without a
// subject, and a query whose letters are valid only in the gap (`x4/chunk_gapvalid.faa`: 100 X, the
// TK-TMK protein's letters 101 to 300, 98 X: neither chunk has a valid context), are searched on and
// reported as NCBI 2.17.0 reports them (exit status 0; rows `db.x4.chunk.neg.gap.nosub.*`,
// `db.x4.chunk.neg.gapvalid.*`), as is psiblast's database left without a subject.
#[test]
fn a_negative_pair_that_leaves_gaps_without_a_subject_or_a_valid_chunk_is_searched_on() {
    let second = "tests/fasta/x4/chunk_second.faa";
    for threads in [1, 4] {
        for (program, query) in [
            ("blastp", second),
            ("psiblast", "tests/fasta/x4/chunk_gapvalid.faa"),
        ] {
            let plain = search_in(program, query, &NO_SUBJECT, "7", threads, &[]);
            let gaps = search_in(program, query, &NO_SUBJECT, "7", threads, &GAP);
            assert_eq!(gaps.code, Some(0), "{program} n{threads}: {}", gaps.stderr);
            assert_eq!(gaps.stdout, plain.stdout, "{program} n{threads}");
            assert_eq!(gaps.stderr, plain.stderr, "{program} n{threads}");
        }
        let valid_in_gap = search("tests/fasta/x4/chunk_gapvalid.faa", "7", threads, &GAP);
        assert_eq!(
            valid_in_gap.code,
            Some(0),
            "n{threads}: {}",
            valid_in_gap.stderr
        );
        assert!(
            valid_in_gap.stdout.contains("# 0 hits found"),
            "n{threads}: {}",
            valid_in_gap.stdout
        );
    }
}

// A split batch whose chunks have no valid context (`x4/chunk_gapvalid.faa` with the gaps): the
// traceback's own setup of the batch fails the subject-length check of `-searchsp` against subjects
// of 2, 0 and 0 letters (NCBI 2.17.0: `Traceback failed`, exit status 3, after the warnings of the
// empty subjects), and a database whose sequence file is missing is not read (exit status 0).
#[test]
fn a_split_batch_without_a_valid_chunk_meets_the_tracebacks_setup() {
    let query = "tests/fasta/x4/chunk_gapvalid.faa";
    for threads in [1, 4] {
        let short = search_in("blastp", query, &SHORT, "7", threads, &GAP);
        assert_eq!(short.code, Some(3), "n{threads}: {}", short.stderr);
        assert!(
            short.stderr.ends_with("Traceback failed\n"),
            "n{threads}: {}",
            short.stderr
        );
        let unread = search_in(
            "blastp",
            query,
            &["-db", "tests/db_search/arch_sprot_nopsq"],
            "7",
            threads,
            &GAP,
        );
        assert_eq!(unread.code, Some(0), "n{threads}: {}", unread.stderr);
        assert!(
            unread.stdout.contains("# 0 hits found"),
            "n{threads}: {}",
            unread.stdout
        );
        // psiblast: the first query's chunks read nothing and it is reported; the next query
        // reads the missing sequence file (NCBI 2.17.0 probe `p6_gapvalid_nopsq2` of line RG item
        // X-bpgap: the database's error, exit status 2, after the first query's report). The
        // error names the file's absolute path, which no frozen row can hold.
        let next = search_in(
            "psiblast",
            "tests/fasta/x4/chunk_gapvalid_first.faa",
            &["-db", "tests/db_search/arch_sprot_nopsq"],
            "7",
            threads,
            &GAP,
        );
        assert_eq!(next.code, Some(2), "n{threads}: {}", next.stderr);
        assert_eq!(
            next.stdout.matches("# Query: ").count(),
            1,
            "n{threads}: {}",
            next.stdout
        );
        assert!(
            next.stdout.contains("# 0 hits found"),
            "n{threads}: {}",
            next.stdout
        );
        assert!(
            next.stderr
                .starts_with("LOSAT database error: Error: File (")
                && next.stderr.ends_with("arch_sprot_nopsq.psq) not found.\n"),
            "n{threads}: {}",
            next.stderr
        );
    }
}

/// A query file of the `index`-th record (from 0) of `source`, or of `letters` (`source` empty), in
/// the temporary directory (`tag` names it).
fn query_file(tag: &str, source: &str, index: usize, letters: &str) -> String {
    let text = if source.is_empty() {
        format!(">{tag}\n{letters}\n")
    } else {
        let all = std::fs::read_to_string(path(source)).expect("read the queries");
        format!(">{}", all.split('>').nth(index + 1).expect("the record"))
    };
    let file = std::env::temp_dir().join(format!("losat_bpgap_{tag}_{}.faa", std::process::id()));
    std::fs::write(&file, text).expect("write the query");
    file.to_str().expect("a path").to_string()
}

// The chunks' stops come before the errors of the setup of the batch's first chunk (NCBI 2.17.0 exits
// 3 with its null pointer, not the database's or the average subject length's error): a chunk
// without a query (`CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200`, ql of `outfmt0/xsb_ql.faa`, 310
// letters) and a chunk split again (`CHUNK_SIZE=150`, the 398 letters of `x4/chunk_second.faa`)
// against a database without sequences and against `-searchsp` with subjects too short. A batch of
// invalid queries only, which NCBI does not search, is not split and meets no stop (exit status 0,
// as without the variables).
#[test]
fn the_chunk_stops_come_before_the_setup_failures() {
    let ql = query_file("ql", "tests/fasta/outfmt0/xsb_ql.faa", 1, "");
    let tk = query_file("tk", "tests/fasta/x4/chunk_second.faa", 1, "");
    let empty = ["-db", "tests/seqdb/db/empty"];
    let without_a_query = [("CHUNK_SIZE", "-100"), ("OVERLAP_CHUNK_SIZE", "-200")];
    let split_again = [("CHUNK_SIZE", "150")];
    for threads in [1, 4] {
        for (query, env, message) in [
            (&ql, &without_a_query[..], "a query chunk without a query"),
            (&tk, &split_again[..], "would split a query chunk again"),
        ] {
            for subjects in [&empty[..], &SHORT[..]] {
                let run = search_in("blastp", query, subjects, "6", threads, env);
                assert_eq!(
                    run.code,
                    Some(1),
                    "{subjects:?} {env:?} n{threads}: {}",
                    run.stderr
                );
                assert!(
                    run.stderr.contains(message),
                    "{subjects:?} {env:?} n{threads}: {}",
                    run.stderr
                );
                assert_eq!(run.stdout, "", "{subjects:?} {env:?} n{threads}");
            }
        }
    }
    let all_x = query_file("allx", "", 0, &"X".repeat(398));
    for env in [&GAP[..], &split_again[..]] {
        let plain = search(&all_x, "7", 1, &[]);
        let run = search(&all_x, "7", 1, env);
        assert_eq!(run.code, Some(0), "{env:?}: {}", run.stderr);
        assert_eq!(
            (run.stdout, run.stderr),
            (plain.stdout, plain.stderr),
            "{env:?}"
        );
    }
    for file in [ql, tk, all_x] {
        let _ = std::fs::remove_file(file);
    }
}
