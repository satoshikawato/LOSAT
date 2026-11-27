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
//! `db.x4.chunk.neg.gap.6`; the step of TBLASTN's since PD-LOSAT-NCBI-DEFECTS 1.98, line RG item
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

// The gaps (PD-LOSAT-NCBI-DEFECTS 1.98, line RG item X-bpgap): refused only where the search of the
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

// A split batch whose chunks have no valid context reads nothing; the next batch, split too, whose
// first valid chunk reaches the read of the missing sequence file fails after the first batch's
// report (`x4/chunk_gapvalid_tk.faa`: the gapvalid record, then the TK-TMK protein, one batch each;
// NCBI 2.17.0 probe `s1_two` of line RG item X-bpgap fix round 6, n1 and n4: `Database memory map
// file error`, exit status 3; row `db.x4.chunk.neg.gapvalid.tk.nopsq.7`).
#[test]
fn a_later_split_batch_reads_the_missing_file_after_the_reports_before() {
    let nopsq = ["-db", "tests/db_search/arch_sprot_nopsq"];
    let first = search_in(
        "blastp",
        "tests/fasta/x4/chunk_gapvalid.faa",
        &nopsq,
        "7",
        1,
        &GAP,
    );
    assert_eq!(first.code, Some(0), "{}", first.stderr);
    let report = &first.stdout[..first
        .stdout
        .rfind("# BLAST processed")
        .expect("outfmt 7's last line")];
    let env = [("BATCH_SIZE", "1"), GAP[0], GAP[1]];
    for threads in [1, 4] {
        let run = search_in(
            "blastp",
            "tests/fasta/x4/chunk_gapvalid_tk.faa",
            &nopsq,
            "7",
            threads,
            &env,
        );
        assert_eq!(run.code, Some(3), "n{threads}: {}", run.stderr);
        assert!(
            run.stderr.ends_with("Database memory map file error\n"),
            "n{threads}: {}",
            run.stderr
        );
        assert_eq!(run.stdout, report, "n{threads}");
    }
}

// A refused batch of two queries with one id (`q1c/avgdup_bp.faa` under `-parse_deflines`) searched
// whole meets NCBI's setup in NCBI's order (line RG item X-bpgap fix round 6, NCBI 2.17.0 probes
// `w2`-`w4`, `u4`, `u9`, `u12`): against `aliaskeys/p_nopsq_short` (no sequence file, an average of
// 0 letters) the average subject length fails before the search's read without the Gumbel block
// (`-searchsp`, an ungapped search) and the read fails with it; an ungapped search, which NCBI never
// splits, is searched whole with sizes that would split it; a count of 0, on which NCBI crashes,
// keeps the batch's refusal.
#[test]
fn a_refused_batch_meets_the_setup_checks_before_the_search_read() {
    #[track_caller]
    fn engine_error(run: Run, error: &str) {
        assert_eq!(run.code, Some(3), "{}", run.stderr);
        assert!(run.stderr.ends_with(error), "{}", run.stderr);
        assert_eq!(run.stdout, "");
    }
    let refused = |subjects: &[&str], ungapped: bool, env: &[(&str, &str)]| {
        let mut arguments = subjects.to_vec();
        arguments.push("-parse_deflines");
        if ungapped {
            arguments.extend(["-ungapped", "-comp_based_stats", "0"]);
        }
        let query = "tests/db_search/q1c/avgdup_bp.faa";
        search_in("blastp", query, &arguments, "6", 1, env)
    };
    let (too_short, unread) = (
        "The average subject length is too short\n",
        "Database memory map file error\n",
    );
    let short = ["-db", "tests/db_search/aliaskeys/p_nopsq_short"];
    let again = [("CHUNK_SIZE", "150")];
    let searchsp = [short[0], short[1], "-searchsp", "1000"];
    engine_error(refused(&searchsp, false, &[]), too_short);
    engine_error(refused(&short, false, &[]), unread);
    engine_error(refused(&short, true, &[]), too_short);
    let subjects = ["-subject", "tests/db_search/q1c/avg_s.faa"];
    engine_error(refused(&subjects, true, &again), too_short);
    let nopsq = ["-db", "tests/db_search/arch_sprot_nopsq"];
    engine_error(refused(&nopsq, true, &again), unread);
    let count_0 = refused(&["-db", "tests/db_search/aliaskeys/p_nseq0"], true, &[]);
    assert_eq!(count_0.code, Some(1), "{}", count_0.stderr);
    assert!(
        count_0.stderr.contains("two query records have the Seq-id"),
        "{}",
        count_0.stderr
    );
}

/// A query file in the temporary directory, removed when it is dropped (also where an assertion
/// fails).
struct QueryFile(String);

impl std::ops::Deref for QueryFile {
    type Target = str;
    fn deref(&self) -> &str {
        &self.0
    }
}

impl Drop for QueryFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// A query file of the `index`-th record (from 0) of `source`, or of `letters` (`source` empty), in
/// the temporary directory (`tag` names it).
fn query_file(tag: &str, source: &str, index: usize, letters: &str) -> QueryFile {
    let text = if source.is_empty() {
        format!(">{tag}\n{letters}\n")
    } else {
        let all = std::fs::read_to_string(path(source)).expect("read the queries");
        format!(">{}", all.split('>').nth(index + 1).expect("the record"))
    };
    let file = std::env::temp_dir().join(format!("losat_bpgap_{tag}_{}.faa", std::process::id()));
    std::fs::write(&file, text).expect("write the query");
    QueryFile(file.to_str().expect("a path").to_string())
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
}

/// Line RG item bp-SF1 (NCBI BLAST+ 2.17.0 probes of the item, remaining-gaps-20261112/D/logs/bp_sf1/PROBES.md):
/// the part of the near run `near` that NCBI writes before the batch that a refusal of the same input with one more
/// batch meets: outfmt 7 up to its last line (`# BLAST processed`), outfmt 0 up to its epilog, which starts with
/// the two line ends before the database's summary (a query's report ends with its line of the search space, as in
/// NCBI's own bytes where an engine error stops it after such reports, row `db.empty.unsearched.bs1.0`), outfmt 5 up
/// to the end of its last iteration, `</Iteration>`, whose line end starts the document's end (item bp-SF2; NCBI's
/// own bytes in row `db.ub.bp.avg.b1.5`), outfmt 6 all of it.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2249-2258
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:731-732
fn reports_before(near: &Run, outfmt: &str) -> String {
    assert_eq!(near.code, Some(0), "near run: {}", near.stderr);
    let end = match outfmt {
        "7" => near
            .stdout
            .rfind("# BLAST processed")
            .expect("outfmt 7's last line"),
        "0" => near
            .stdout
            .rfind("\n\n  Database: ")
            .expect("outfmt 0's epilog"),
        "5" => near
            .stdout
            .rfind("\n</BlastOutput_iterations>")
            .expect("outfmt 5's end"),
        _ => near.stdout.len(),
    };
    near.stdout[..end].to_string()
}

/// The block search's settings of the item's runs: its default and the search of all queries.
const BLOCKS: [&[(&str, &str)]; 2] = [&[], &[("LOSAT_X_BPQBLOCK", "0")]];

// A statistics count of 0 without the Gumbel block (an alias file's `NSEQ 0`, `OLD_FSC`): NCBI
// writes the reports of the batches before the first batch with a valid query (their queries are
// X-only), then divides by the count in that batch's setup (SIGFPE). LOSAT writes them (without
// the epilog) and refuses that batch, with the block search and without it (rows
// `db.ak.bp.rej.nseq0.oldfsc.blocks*`, `db.ak.bp.rej.nseq0.oldfsc.b1*`).
#[test]
fn the_count_of_0_is_refused_after_the_unsearched_batches_before() {
    let db = ["-db", "tests/db_search/aliaskeys/p_nseq0"];
    for (query, near_query, batch_size) in [
        (
            "tests/fasta/rps/x6q.faa",
            "tests/fasta/rps/x6only.faa",
            None,
        ),
        (
            "tests/fasta/rps/allx_q649.faa",
            "tests/fasta/rps/allx.faa",
            Some("1"),
        ),
    ] {
        for outfmt in ["7", "0", "6"] {
            let mut env = vec![("OLD_FSC", "1")];
            if let Some(size) = batch_size {
                env.push(("BATCH_SIZE", size));
            }
            let near = search_in("blastp", near_query, &db, outfmt, 1, &env);
            let before = reports_before(&near, outfmt);
            for blocks in BLOCKS {
                for threads in [1, 4] {
                    let mut env = env.clone();
                    env.extend_from_slice(blocks);
                    let run = search_in("blastp", query, &db, outfmt, threads, &env);
                    let case = format!("{query} outfmt {outfmt} n{threads} {blocks:?}");
                    assert_eq!(run.code, Some(1), "{case}: {}", run.stderr);
                    assert!(
                        run.stderr
                            .contains("the statistics count 0 sequences of the database"),
                        "{case}: {}",
                        run.stderr
                    );
                    assert_eq!(run.stdout, before, "{case}");
                }
            }
        }
    }
}

// A search space that overflows NCBI's `Int8` (an alias file's `LENGTH` of 2^55, or `-dbsize` of
// it, with a query of 4,970 letters; NCBI wraps it and goes on): LOSAT writes the reports of the
// batches before, with their hits, and refuses the batch of that query; in one batch the refusal
// comes before any report (rows `db.ak.bp.rej.l55.*`).
#[test]
fn a_search_space_overflow_is_refused_after_the_batches_before() {
    let alias = ["-db", "tests/db_search/aliaskeys/p_l55"];
    let dbsize = [
        "-subject",
        "tests/db_search/taxfilter/tfp.faa",
        "-dbsize",
        "36028797018963968",
    ];
    let env = [("BATCH_SIZE", "1")];
    for (subjects, outfmts) in [
        (&alias[..], &["7", "0", "6"][..]),
        (&dbsize[..], &["7", "0", "6"]),
    ] {
        for &outfmt in outfmts {
            let near = search_in(
                "blastp",
                "tests/db_search/aliaskeys/q_hba.faa",
                subjects,
                outfmt,
                1,
                &env,
            );
            let before = reports_before(&near, outfmt);
            let hits = if outfmt == "0" {
                "Sequences producing significant alignments"
            } else {
                "hba\t"
            };
            assert!(before.contains(hits), "the near run's hits: {before}");
            for blocks in BLOCKS {
                for threads in [1, 4] {
                    let mut env = env.to_vec();
                    env.extend_from_slice(blocks);
                    let run = search_in(
                        "blastp",
                        "tests/db_search/aliaskeys/q_hba_long.faa",
                        subjects,
                        outfmt,
                        threads,
                        &env,
                    );
                    let case = format!("{subjects:?} outfmt {outfmt} n{threads} {blocks:?}");
                    assert_eq!(run.code, Some(1), "{case}: {}", run.stderr);
                    assert!(
                        run.stderr.contains(
                            "with a query of 4970 letters overflows NCBI BLAST+'s 64-bit effective search space"
                        ),
                        "{case}: {}",
                        run.stderr
                    );
                    assert_eq!(run.stdout, before, "{case}");
                }
            }
        }
        let one_batch = search_in(
            "blastp",
            "tests/db_search/aliaskeys/q_hba_long.faa",
            subjects,
            "7",
            1,
            &[],
        );
        assert_eq!(
            one_batch.code,
            Some(1),
            "{subjects:?}: {}",
            one_batch.stderr
        );
        assert_eq!(one_batch.stdout, "", "{subjects:?}");
    }
}

// Review M2 of item bp-SF1: a refusal that the call's first batch meets follows the prolog of
// outfmt 0 to 4 (the program, the references and the database), which NCBI writes before it sets
// the first batch up; outfmt 5 and 7 write nothing before it (rows `db.ak.bp.rej.l55.b1.*`,
// `db.ak.bp.rej.nseq0.oldfsc.one.0`).
#[test]
fn a_first_batch_refusal_follows_the_pairwise_prolog() {
    // The prolog ends before the first query's report, which starts with two line ends (NCBI's
    // bytes where an engine error stops the first batch, as in rows `db.empty.0` and
    // `db.t5d.blastp.sw.0`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1489-1494
    let prolog = |near: &Run| {
        assert_eq!(near.code, Some(0), "near run: {}", near.stderr);
        let end = near.stdout.find("\n\nQuery= ").expect("a query's report");
        near.stdout[..end].to_string()
    };
    let overflow =
        "with a query of 4970 letters overflows NCBI BLAST+'s 64-bit effective search space";
    let count = "the statistics count 0 sequences of the database";
    let cases: [(&str, &[(&str, &str)], &str, &str, &str); 2] = [
        (
            "tests/db_search/aliaskeys/p_l55",
            &[],
            "tests/db_search/aliaskeys/q_hba.faa",
            "tests/db_search/aliaskeys/q_hba_long.faa",
            overflow,
        ),
        (
            "tests/db_search/aliaskeys/p_nseq0",
            &[("OLD_FSC", "1")],
            "tests/fasta/rps/x6only.faa",
            "tests/db_search/aliaskeys/q_hba.faa",
            count,
        ),
    ];
    for (db, env, near_query, query, refusal) in cases {
        let db = ["-db", db];
        for outfmt in ["0", "3", "5", "7"] {
            let before = match outfmt {
                "0" | "3" => prolog(&search_in("blastp", near_query, &db, outfmt, 1, env)),
                _ => String::new(),
            };
            for blocks in BLOCKS {
                for threads in [1, 4] {
                    let mut env = env.to_vec();
                    env.extend_from_slice(blocks);
                    let run = search_in("blastp", query, &db, outfmt, threads, &env);
                    let case = format!("{db:?} outfmt {outfmt} n{threads} {blocks:?}");
                    assert_eq!(run.code, Some(1), "{case}: {}", run.stderr);
                    assert!(run.stderr.contains(refusal), "{case}: {}", run.stderr);
                    assert_eq!(run.stdout, before, "{case}");
                }
            }
        }
    }
}

// NCBI computes no search space for a query that is not valid (an X-only one): such queries over
// a length whose space for them would overflow `Int8` are reported as NCBI reports them (exit 0;
// row `db.ak.bp.l55.x6.7`).
#[test]
fn invalid_queries_meet_no_search_space_overflow() {
    let run = search_in(
        "blastp",
        "tests/fasta/rps/x6only.faa",
        &["-db", "tests/db_search/aliaskeys/p_l55"],
        "7",
        1,
        &[],
    );
    assert_eq!(run.code, Some(0), "{}", run.stderr);
    // An invalid query's outfmt 7 report is its header lines (no hit count), as NCBI's.
    assert_eq!(
        run.stdout.matches("\n# Query: x").count(),
        6,
        "{}",
        run.stdout
    );
    assert!(
        run.stdout.ends_with("# BLAST processed 6 queries\n"),
        "{}",
        run.stdout
    );
}

// An ungapped search over a subject longer than NCBI's part length (5,000,000 residues), which
// LOSAT refuses at the first batch with a valid query: the reports of the batches before (of
// X-only queries) come first.
#[test]
fn the_long_subject_refusal_of_an_ungapped_search_follows_the_unsearched_batches() {
    let dir = std::env::temp_dir().join(format!("losat-bp-sf1-long-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temporary directory");
    let subject = dir.join("long.faa");
    let mut letters = "ACDEFGHIKLMNPQRSTVWY".repeat(250_001);
    letters.truncate(5_000_001);
    std::fs::write(&subject, format!(">long\n{letters}\n")).expect("long subject");
    let subject = subject.to_str().expect("a UTF-8 path").to_string();
    // An ungapped search takes only `-comp_based_stats 0` (NCBI's option check).
    let subjects = [
        "-subject",
        subject.as_str(),
        "-ungapped",
        "-comp_based_stats",
        "0",
    ];
    let env = [("BATCH_SIZE", "1")];
    let near = search_in(
        "blastp",
        "tests/fasta/rps/allx.faa",
        &subjects,
        "7",
        1,
        &env,
    );
    let before = reports_before(&near, "7");
    for threads in [1, 4] {
        let run = search_in(
            "blastp",
            "tests/fasta/rps/allx_q649.faa",
            &subjects,
            "7",
            threads,
            &env,
        );
        assert_eq!(run.code, Some(1), "n{threads}: {}", run.stderr);
        assert!(
            run.stderr
                .contains("subject sequence 1 has 5000001 residues"),
            "n{threads}: {}",
            run.stderr
        );
        assert_eq!(run.stdout, before, "n{threads}");
    }
    std::fs::remove_dir_all(&dir).ok();
}

// Line RG item bp-SF2: NCBI's engine error of a batch's setup (`The average subject length is too
// short`: `OLD_FSC` over subjects of 2, 0 and 0 letters) after first batches that NCBI does not
// search (`BATCH_SIZE=1`: two X-only queries, then a valid one): the reports of those batches come
// first, without the epilog, then the error (exit status 3), in every format (rows `db.ub.*`). The
// same bytes with the block search's blocks of one batch (`LOSAT_X_BPQBLOCK=1`: the failing block
// starts with the valid query), of x1 and then x2 with the valid query (`=13`: the failing block
// starts with a batch that NCBI does not search), its default (one block, the search of all
// queries) and without it (`=0`).
#[test]
fn the_average_subject_length_error_follows_the_unsearched_batches_in_every_format() {
    let subjects = ["-subject", "tests/db_search/q1c/avg_s.faa"];
    let env = [("OLD_FSC", "1"), ("BATCH_SIZE", "1")];
    let blocks: [&[(&str, &str)]; 4] = [
        &[],
        &[("LOSAT_X_BPQBLOCK", "1")],
        &[("LOSAT_X_BPQBLOCK", "13")],
        &[("LOSAT_X_BPQBLOCK", "0")],
    ];
    for outfmt in ["0", "7", "6", "5", "11", "15", "10"] {
        let before = matches!(outfmt, "0" | "7" | "6" | "5").then(|| {
            let near = search_in(
                "blastp",
                "tests/fasta/rps/allx.faa",
                &subjects,
                outfmt,
                1,
                &env,
            );
            reports_before(&near, outfmt)
        });
        let mut first: Option<Run> = None;
        for blocks in blocks {
            for threads in [1, 4] {
                let mut env = env.to_vec();
                env.extend_from_slice(blocks);
                let run = search_in(
                    "blastp",
                    "tests/fasta/rps/allx_q649.faa",
                    &subjects,
                    outfmt,
                    threads,
                    &env,
                );
                let case = format!("outfmt {outfmt} n{threads} {blocks:?}");
                assert_eq!(run.code, Some(3), "{case}: {}", run.stderr);
                assert!(
                    run.stderr
                        .ends_with("engine error: The average subject length is too short\n"),
                    "{case}: {}",
                    run.stderr
                );
                if let Some(before) = &before {
                    assert_eq!(&run.stdout, before, "{case}");
                }
                match &first {
                    Some(first) => {
                        assert_eq!(run.stdout, first.stdout, "{case}");
                        assert_eq!(run.stderr, first.stderr, "{case}");
                    }
                    None => first = Some(run),
                }
            }
        }
    }
}
