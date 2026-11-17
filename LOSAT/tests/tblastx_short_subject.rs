//! TBLASTX on a subject (or its `-subject_loc` range) of fewer than 12 bases whose HSPs reach the
//! ungapped re-evaluation, where NCBI BLAST+ 2.17.0 reads the translation of a frame it never makes
//! and dies with SIGSEGV: LOSAT refuses there (PD-LOSAT-NCBI-DEFECTS 1.50, rule 3). The single
//! searches are rows `db.xs.tblastx.*` of `tests/db_search_manifest/tblastx.tsv`; here, at 1 and 4
//! threads, what those rows cannot hold: the reports of the query batches before the refusal, a
//! subject that the re-evaluation's fetch excludes (a negative seqid list), which NCBI never reads,
//! and a `-db` subject with ambiguous bases, whose letters the re-evaluation reads again.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2659-2662,2700
//! ```c
//!       	if (( BLAST_SEQSRC_EXCLUDED == BlastSeqSrcGetSequence(seq_src, &seq_arg))) {
//! ...
//!          subject_start = Blast_HSPGetTargetTranslation(target_t, hsp, NULL);
//! ```
//! The inputs (NCBI's outcomes from the probes of line RG item X-short-tbx on ws, 2026-11-17):
//! `nq` (`ATGCATTATTGGTGGTGGTTT`, frame +1 M H Y W W W F) hits `w11` (`TGGTGGTGGTT`) with
//! `-window_size 0`, and NCBI dies; so it does on `amb11` (`TGGTGGTGGNN`) with `-db`. `qa` (30 A)
//! hits the 33 A of `pa` and nothing of `w11`. With `BATCH_SIZE=1` the query `qa nq` is two
//! batches: NCBI formats the report of `qa` and dies in the second, and what reaches its standard
//! output is what its output buffer had passed on before the crash: with `-outfmt 6` as many bytes
//! as the report of `qa` alone (11,318), with `-outfmt 0` the prolog and the report of `qa` but its
//! last few bytes (50,001 bytes); in one batch, nothing with `-outfmt 6` and the first 18 bytes of
//! the prolog with `-outfmt 0`. LOSAT writes those reports and prologs whole. On a database of
//! `w12` and `w11` with parsed Seq-ids, `-negative_seqidlist` of `w11` gives NCBI's report of `w12`
//! (the same bytes as `-seqidlist` of `w12`), and of `w12` NCBI's crash.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The refusal's message on standard error.
const REFUSAL: &str = "TBLASTX: the ungapped re-evaluation of an alignment on a subject (or its \
-subject_loc range) of fewer than 12 bases, where NCBI BLAST+ reads the translation of the \
subject's frame, which it never makes for a subject this short, and crashes; LOSAT does not \
support it";

const SEARCH: [&str; 6] = ["-window_size", "0", "-seg", "no", "-evalue", "1e10"];

/// The start of `-outfmt 0`'s epilog (`PrintEpilog`) of a `-subject` search, after the last
/// query's report, whose own start is `\n\nQuery= ` after the prolog.
const EPILOG: &str = "\n\n  Database: User specified sequence set";

/// A temporary directory of this test, removed when it is dropped (also when an assertion fails).
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastx-short-subject-{tag}-{}",
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

/// Writes the FASTA `records` (id, letters) to `dir/name`.
fn fasta(dir: &Path, name: &str, records: &[(&str, &str)]) -> PathBuf {
    let path = dir.join(name);
    let text: String = records
        .iter()
        .map(|(id, letters)| format!(">{id}\n{letters}\n"))
        .collect();
    std::fs::write(&path, text).unwrap();
    path
}

/// `losat ARGS` in `dir` with no `LOSAT_*`, `NCBI*` or `BLASTDB` variable but
/// `LOSAT_BLASTX_WORKER_LOG` (a binary built with the `blastx-worker-probe` feature needs that file
/// name once a search has more than one thread; it selects nothing), `BATCH_SIZE` only as
/// `batch_size` says.
fn losat(dir: &Path, batch_size: Option<&str>, args: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || matches!(text.as_ref(), "BATCH_SIZE" | "BLASTDB")
        {
            command.env_remove(name);
        }
    }
    if let Some(size) = batch_size {
        command.env("BATCH_SIZE", size);
    }
    command
        .current_dir(dir)
        .env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"))
        .args(args)
        .output()
        .expect("run LOSAT")
}

/// The output of a refused search: exit status 1, standard error ending with the refusal (after
/// NCBI's warnings of the search, such as that of a seqid list in text).
fn assert_refused(output: &Output, case: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{case}: {stderr}");
    assert!(
        stderr.ends_with(&format!("Error: {REFUSAL}\n")) && stderr.matches(REFUSAL).count() == 1,
        "{case}: {stderr}"
    );
}

/// The report of `qa`, the first of two query batches, is written whole before the refusal of the
/// second (NCBI formats it before its crash, which loses what NCBI's output buffer still holds); in
/// the tabular report written batch by batch, in the one with a field list written after the search
/// (`report_stream`), and with `-outfmt 0`, which has the prolog and the report of `qa` and leaves
/// out the epilog. In one batch the refusal comes before any report: only `-outfmt 0`'s prolog.
#[test]
fn the_reports_of_the_query_batches_before_the_refusal_are_written() {
    let temp = TempDir::new("batches");
    let dir = temp.0.as_path();
    let qa = "A".repeat(30);
    let query = fasta(
        dir,
        "qq.fna",
        &[("qa", qa.as_str()), ("nq", "ATGCATTATTGGTGGTGGTTT")],
    );
    let first = fasta(dir, "qa.fna", &[("qa", &qa)]);
    let subject = fasta(
        dir,
        "s.fna",
        &[("pa", &"A".repeat(33)), ("w11", "TGGTGGTGGTT")],
    );
    for outfmt in ["6", "6 std", "0"] {
        for threads in ["1", "4"] {
            let search = |query: &Path, batch_size| {
                let mut args = vec!["tblastx", "-query", query.to_str().unwrap()];
                args.extend(["-subject", subject.to_str().unwrap()]);
                args.extend(SEARCH);
                args.extend(["-outfmt", outfmt, "-num_threads", threads]);
                losat(dir, batch_size, &args)
            };
            let case = format!("outfmt {outfmt} threads {threads}");
            let alone = search(&first, None);
            assert!(
                alone.status.success() && !alone.stdout.is_empty(),
                "{case}: qa alone"
            );
            let refused = search(&query, Some("1"));
            assert_refused(&refused, &case);
            let one_batch = search(&query, None);
            assert_refused(&one_batch, &case);
            if outfmt == "0" {
                // The output of `qa` alone is the prolog, the report of `qa` and the epilog.
                let text = String::from_utf8_lossy(&alone.stdout);
                let report = text.find("\n\nQuery= qa").expect("the report of qa");
                let epilog = report + text[report..].find(EPILOG).expect("the epilog");
                assert_eq!(
                    refused.stdout,
                    alone.stdout[..epilog],
                    "{case}: the prolog and the report of qa"
                );
                assert_eq!(
                    one_batch.stdout,
                    alone.stdout[..report],
                    "{case}: one batch, the prolog"
                );
            } else {
                assert_eq!(refused.stdout, alone.stdout, "{case}: the report of qa");
                assert!(one_batch.stdout.is_empty(), "{case}: one batch, no report");
            }
        }
    }
}

/// `losat makeblastdb` of the nucleotide FASTA `records` as database `name` in `dir`, with parsed
/// Seq-ids when `parse_seqids`.
fn makeblastdb(dir: &Path, name: &str, records: &[(&str, &str)], parse_seqids: bool) {
    let set = fasta(dir, &format!("{name}.fna"), records);
    let mut args = vec!["makeblastdb", "-in", set.to_str().unwrap()];
    args.extend(["-dbtype", "nucl", "-out", name]);
    if parse_seqids {
        args.push("-parse_seqids");
    }
    let made = losat(dir, None, &args);
    assert!(
        made.status.success(),
        "makeblastdb {name}: {}",
        String::from_utf8_lossy(&made.stderr)
    );
}

/// A subject of 11 bases that the negative seqid list excludes leaves the search at the
/// re-evaluation's fetch, before NCBI reads its translation: no refusal, the report of the other
/// subject; the same list without it is refused.
#[test]
fn a_subject_the_fetch_excludes_is_not_refused() {
    let temp = TempDir::new("fetch");
    let dir = temp.0.as_path();
    let set = [("w12", "TGGTGGTGGTTT"), ("w11", "TGGTGGTGGTT")];
    makeblastdb(dir, "w12w11", &set, true);
    let query = fasta(dir, "nq.fna", &[("nq", "ATGCATTATTGGTGGTGGTTT")]);
    std::fs::write(dir.join("w11.txt"), "w11\n").unwrap();
    std::fs::write(dir.join("w12.txt"), "w12\n").unwrap();
    for threads in ["1", "4"] {
        let search = |list: &str, file: &str| {
            let mut args = vec!["tblastx", "-query", query.to_str().unwrap()];
            args.extend(["-db", "w12w11", list, file]);
            args.extend(SEARCH);
            args.extend(["-outfmt", "6", "-num_threads", threads]);
            losat(dir, None, &args)
        };
        let excluded = search("-negative_seqidlist", "w11.txt");
        let kept = search("-seqidlist", "w12.txt");
        assert_eq!(excluded.status.code(), Some(0), "threads {threads}");
        assert!(
            !excluded.stdout.is_empty() && excluded.stdout == kept.stdout,
            "threads {threads}: the report of w12"
        );
        assert_refused(
            &search("-negative_seqidlist", "w12.txt"),
            &format!("threads {threads}: w11 kept"),
        );
    }
}

/// A `-db` subject of 11 bases with ambiguous bases (`amb11`, W W W in frame +1), whose letters the
/// re-evaluation's fetch reads again (NCBI's too, before it reads the translation it never made and
/// dies): refused, as without them.
#[test]
fn a_database_subject_with_ambiguous_bases_is_refused() {
    let temp = TempDir::new("amb");
    let dir = temp.0.as_path();
    makeblastdb(dir, "amb11", &[("amb11", "TGGTGGTGGNN")], false);
    let query = fasta(dir, "nq.fna", &[("nq", "ATGCATTATTGGTGGTGGTTT")]);
    for threads in ["1", "4"] {
        let mut args = vec!["tblastx", "-query", query.to_str().unwrap()];
        args.extend(["-db", "amb11"]);
        args.extend(SEARCH);
        args.extend(["-outfmt", "6", "-num_threads", threads]);
        let output = losat(dir, None, &args);
        assert_refused(&output, &format!("threads {threads}"));
        assert!(output.stdout.is_empty(), "threads {threads}: no report");
    }
}
