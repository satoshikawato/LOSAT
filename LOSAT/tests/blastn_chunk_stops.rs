//! BLASTN's stops at query chunks that NCBI BLAST+ 2.17.0 cannot make (line RG item X-chunk). NCBI
//! searches each query batch and prints its reports before it reads the next (blastn_app.cpp:
//! 295-318), and `CBlastPrelimSearch::Run` splits a batch of two or more chunks first:
//! `CQuerySplitter::Split` makes the query factory of every chunk before any chunk is set up and
//! dereferences the missing query vector of a chunk that no query reaches (split_query_cxx.cpp:
//! 234-249,876-880): a null-pointer `CCoreException`, exit 3, at every thread count, after the
//! reports of the batches before. LOSAT writes those reports, without the epilog (for the first
//! batch, the prologs), and then refuses (`ErrorAfterInvalidBatches`), in every format, task and
//! thread count; nothing of the refused batch or of the batches after it is written. The single
//! searches are rows `db.xc.*` of `tests/db_search_manifest/blastn.xc.tsv`, NCBI's stdout frozen;
//! here what those rows cannot hold: every format, both tasks, the thread counts, `-db` and the
//! other order of the stops.
//!
//! The inputs (`tests/db_search/make_xchunk_inputs.py`): the queries are pieces of the subject
//! `pa` (4000 bases), so the first query of a pair hits it; with `BATCH_SIZE=1` the first query is
//! the first query batch. `CHUNK_SIZE=140 OVERLAP_CHUNK_SIZE=100` split the 90 bases of the second
//! query into 2 chunks of 96, below the overlap: the start of the second wraps past the query
//! (`size_t`) and that chunk holds no query. `-100/-200` split 310 bases into 3 chunks of which the
//! first covers the query (the others hold none), and 248 bases into 2 chunks of 25 with a gap of
//! 200 between them (NCBI searches on and loses the hits there; LOSAT refuses it before any
//! output). `150/100` split 1000 bases into 20 chunks of 146 that NCBI would split again (the
//! stop of exception 1 of PD-LOSAT-NCBI-DEFECTS, which LOSAT does not make).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The refusal of a query chunk without a query.
const EMPTY_CHUNK: &str = "a query chunk without a query is not supported by LOSAT's BLASTN";

/// The refusal of a negative `CHUNK_SIZE` whose wrapped chunk ranges leave gaps between them.
const NEGATIVE_GAPS: &str =
    "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this \
query batch with NCBI BLAST+'s size_t chunk ranges wrapped, which is not supported by LOSAT's \
BLASTN";

/// The start of `-outfmt 0`'s epilog (`PrintEpilog`) of a `-subject` search, after the last
/// query's report.
const EPILOG: &str = "\n\n  Database: User specified sequence set";

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta/outfmt0")
        .join(name)
}

/// A temporary directory of this test, removed when it is dropped (also when an assertion fails).
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-blastn-chunk-stops-{tag}-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The first record of the FASTA file `name`, written to `dir/<name>.first`.
fn first_record(dir: &Path, name: &str) -> PathBuf {
    let text = std::fs::read_to_string(fixture(name)).unwrap();
    let end = text[1..].find('>').expect("a second record") + 1;
    let path = dir.join(format!("{name}.first"));
    std::fs::write(&path, &text[..end]).unwrap();
    path
}

/// `losat ARGS` in `dir` with no `LOSAT_*`, `NCBI*`, `BATCH_SIZE`, `CHUNK_SIZE`,
/// `OVERLAP_CHUNK_SIZE` or `BLASTDB` variable but `LOSAT_BLASTX_WORKER_LOG` (a binary built with
/// the `blastx-worker-probe` feature needs that file name once a search has more than one thread;
/// it selects nothing), and the variables `env`.
fn losat(dir: &Path, env: &[(&str, &str)], args: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "BLASTDB" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .current_dir(dir)
        .env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"))
        .envs(env.iter().copied())
        .args(args)
        .output()
        .expect("run LOSAT")
}

/// `blastn -task TASK -query QUERY SUBJECT -evalue 1e10 -outfmt OUTFMT -num_threads THREADS`.
fn search(
    dir: &Path,
    env: &[(&str, &str)],
    task: &str,
    query: &Path,
    subject: &[&str],
    outfmt: &str,
    threads: &str,
) -> Output {
    let mut args = vec!["blastn", "-task", task, "-query", query.to_str().unwrap()];
    args.extend(subject);
    args.extend([
        "-evalue",
        "1e10",
        "-outfmt",
        outfmt,
        "-num_threads",
        threads,
    ]);
    losat(dir, env, &args)
}

/// `actual` equals `expected`, both shown as text when they differ.
fn assert_bytes(actual: &[u8], expected: &[u8], what: &str) {
    assert!(
        actual == expected,
        "{what}:\n{}\n--- expected ---\n{}",
        String::from_utf8_lossy(actual),
        String::from_utf8_lossy(expected)
    );
}

/// The output of a refused search: exit status 1, standard error ending with `message` once.
fn assert_refused(output: &Output, message: &str, case: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{case}: {stderr}");
    assert!(
        stderr.ends_with(&format!("Error: {message}\n")) && stderr.matches(message).count() == 1,
        "{case}: {stderr}"
    );
}

/// What NCBI has written when it stops in the second query batch: the report of the first query
/// as it writes it alone, without what comes after the last report (the epilog of `-outfmt 0`, the
/// end of the XML document, the count line of `-outfmt 7`); the tabular reports are whole. The XML2
/// and JSON documents hold the first batch's report back (NCBI writes it with the second), so
/// only their prolog is written.
fn written_before_the_stop(alone: &[u8], outfmt: &str) -> Vec<u8> {
    let text = String::from_utf8_lossy(alone);
    let end = match outfmt {
        "0" => text.find(EPILOG).expect("the epilog"),
        "5" => text
            .rfind("\n</BlastOutput_iterations>")
            .expect("the end of the document"),
        "7" => text.rfind("# BLAST processed").expect("the count line"),
        "15" => "{\n\"BlastOutput2\": [\n".len(),
        "16" => text.find("\n<BlastOutput2>").expect("the first report"),
        _ => text.len(),
    };
    alone[..end].to_vec()
}

/// The formats, tasks and thread counts of the stops of a split batch.
const OUTFMTS: [&str; 6] = ["0", "5", "6", "7", "15", "16"];

/// A query chunk without a query: the second query batch stops at `CQuerySplitter::Split`, after
/// the report of the first, for the three ways to leave one (an evened chunk size below the
/// overlap, a positive pair on 31 bases, a negative pair whose first range covers the batch).
#[test]
fn a_chunk_without_a_query_stops_after_the_reports_of_the_batches_before() {
    let temp = TempDir::new("empty");
    let dir = temp.0.as_path();
    let subject_file = fixture("xch_s.fna");
    let subject = ["-subject", subject_file.to_str().unwrap()];
    for (name, chunk, overlap) in [
        ("xch_under.fna", "140", "100"),
        ("xch_p95.fna", "100", "95"),
        ("xch_neg.fna", "-100", "-200"),
    ] {
        let first = first_record(dir, name);
        let query = fixture(name);
        let env = [
            ("BATCH_SIZE", "1"),
            ("CHUNK_SIZE", chunk),
            ("OVERLAP_CHUNK_SIZE", overlap),
        ];
        for task in ["blastn", "megablast"] {
            for outfmt in OUTFMTS {
                for threads in ["1", "4"] {
                    let case = format!("{name} {task} outfmt {outfmt} threads {threads}");
                    let alone = search(dir, &[], task, &first, &subject, outfmt, threads);
                    assert!(alone.status.success(), "{case}: the first query alone");
                    let refused = search(dir, &env, task, &query, &subject, outfmt, threads);
                    assert_refused(&refused, EMPTY_CHUNK, &case);
                    assert_bytes(
                        &refused.stdout,
                        &written_before_the_stop(&alone.stdout, outfmt),
                        &format!("{case}: the report of the first query"),
                    );
                }
            }
        }
    }
}

/// In the first query batch nothing is written but the prologs (`-outfmt 6` and 7 have none).
#[test]
fn a_chunk_without_a_query_in_the_first_batch_writes_the_prologs_alone() {
    let temp = TempDir::new("first");
    let dir = temp.0.as_path();
    let subject_file = fixture("xch_s.fna");
    let subject = ["-subject", subject_file.to_str().unwrap()];
    let query = fixture("xch_under.fna");
    let text = std::fs::read_to_string(&query).unwrap();
    let second = dir.join("second.fna");
    std::fs::write(&second, &text[text[1..].find('>').unwrap() + 1..]).unwrap();
    let env = [("CHUNK_SIZE", "140"), ("OVERLAP_CHUNK_SIZE", "100")];
    for outfmt in ["0", "5", "6", "7", "15", "16"] {
        for threads in ["1", "4"] {
            let case = format!("outfmt {outfmt} threads {threads}");
            let refused = search(dir, &env, "blastn", &second, &subject, outfmt, threads);
            assert_refused(&refused, EMPTY_CHUNK, &case);
            assert!(
                matches!(outfmt, "0" | "5" | "15" | "16") || refused.stdout.is_empty(),
                "{case}: {}",
                String::from_utf8_lossy(&refused.stdout)
            );
            assert!(
                !String::from_utf8_lossy(&refused.stdout).contains("Query= q90"),
                "{case}: nothing of the refused batch"
            );
        }
    }
}

/// A negative pair whose wrapped chunk ranges leave gaps between them is refused before the first
/// byte (NCBI searches on and loses the hits in the gaps): no report of the batch before it, at
/// every thread count; a negative pair that does not split the batch is searched as without the
/// variables.
#[test]
fn a_negative_pair_with_gaps_is_refused_before_any_output() {
    let temp = TempDir::new("gaps");
    let dir = temp.0.as_path();
    let subject_file = fixture("xch_s.fna");
    let subject = ["-subject", subject_file.to_str().unwrap()];
    let query = fixture("xch_gaps.fna");
    let negative = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-100"),
        ("OVERLAP_CHUNK_SIZE", "-200"),
    ];
    for task in ["blastn", "megablast"] {
        for outfmt in ["0", "6", "7", "5"] {
            for threads in ["1", "4"] {
                let case = format!("{task} outfmt {outfmt} threads {threads}");
                let refused = search(dir, &negative, task, &query, &subject, outfmt, threads);
                assert_refused(&refused, NEGATIVE_GAPS, &case);
                assert!(
                    !refused.stdout.windows(4).any(|word| word == b"q150"),
                    "{case}: the report of q150 is not written"
                );
            }
        }
    }
    // -1 over INT_MIN does not split these 1,150 letters: NCBI searches as without the variables.
    let wide = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-1"),
        ("OVERLAP_CHUNK_SIZE", "-2147483648"),
    ];
    let plain = search(
        dir,
        &[("BATCH_SIZE", "1")],
        "blastn",
        &query,
        &subject,
        "6",
        "1",
    );
    let searched = search(dir, &wide, "blastn", &query, &subject, "6", "1");
    assert!(plain.status.success() && !plain.stdout.is_empty());
    assert_bytes(&searched.stdout, &plain.stdout, "-1 over INT_MIN");
}

/// A chunk that NCBI would split again is searched once (PD-LOSAT-NCBI-DEFECTS exception 1): the
/// reports of both batches are written and the exit status is 0, at every thread count (NCBI
/// stops at one thread, and logs the exception and goes on at more).
#[test]
fn a_chunk_split_again_is_searched_once() {
    let temp = TempDir::new("again");
    let dir = temp.0.as_path();
    let subject_file = fixture("xch_s.fna");
    let subject = ["-subject", subject_file.to_str().unwrap()];
    let query = fixture("xch_again.fna");
    for threads in ["1", "4"] {
        let case = format!("150/100 threads {threads}");
        let env = [
            ("BATCH_SIZE", "1"),
            ("CHUNK_SIZE", "150"),
            ("OVERLAP_CHUNK_SIZE", "100"),
        ];
        let output = search(dir, &env, "blastn", &query, &subject, "6", threads);
        let text = String::from_utf8_lossy(&output.stdout);
        assert!(
            output.status.success() && text.contains("q70\t") && text.contains("q1000\t"),
            "{case}: {text} {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}
