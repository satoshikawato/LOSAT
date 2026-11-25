//! TBLASTX on a subject (or its `-subject_loc` range) of fewer than 12 bases whose HSPs reach the
//! ungapped re-evaluation, where NCBI BLAST+ 2.17.0 reads the translation of a frame it never makes
//! and dies with SIGSEGV: LOSAT refuses there (PD-LOSAT-NCBI-DEFECTS 1.51, rule 3). The single
//! searches are rows `db.xs.tblastx.*` of `tests/db_search_manifest/tblastx.tsv`; here, at 1 and 4
//! threads, what those rows cannot hold: the reports of the query batches before the refusal, a
//! subject that the re-evaluation's fetch excludes (a negative seqid list), which NCBI never reads,
//! and a `-db` subject with ambiguous bases, whose letters the re-evaluation reads again.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2659-2662,2700
//! The inputs (NCBI's outcomes from the probes of line RG item X-short-tbx on ws, 2026-11-17):
//! `nq` (`ATGCATTATTGGTGGTGGTTT`, frame +1 M H Y W W W F) hits `w11` (`TGGTGGTGGTT`) with
//! `-window_size 0`, and NCBI dies; so it does on `amb11` (`TGGTGGTGGNN`) with `-db`. `qa` (30 A)
//! hits the 33 A of `pa` and nothing of `w11`. With `BATCH_SIZE=1` the query `qa nq` is two
//! batches: NCBI formats the report of `qa` and dies in the second, and what reaches its standard
//! output is what its output buffer had passed on before the crash: with `-outfmt 6` as many bytes
//! as the report of `qa` alone (11,318), with `-outfmt 0` the prolog and the report of `qa` but its
//! footer's last 86 bytes (50,001 of the 50,087 bytes LOSAT writes); in one batch, nothing with
//! `-outfmt 6` and the first 18 bytes of the prolog with `-outfmt 0`. LOSAT writes those reports
//! and prologs whole. On a database of `w12` and `w11` with parsed Seq-ids, `-negative_seqidlist`
//! of `w11` gives NCBI's report of `w12` (the same bytes as `-seqidlist` of `w12`), and of `w12`
//! NCBI's crash.
//!
//! The re-evaluation's fetch under a negative seqid list reads the header of each subject that
//! still has HSPs (blast_hits.c:2651-2665, seqsrc_seqdb.cpp:306-311): one that cannot be read ends
//! NCBI's search with an engine error (exit 3) after the reports of the query batches before, and
//! LOSAT refuses there, after them too (line RG item X-short-tbx review 2, probes on ws
//! 2026-11-18, whose NCBI bytes the test compares).
//!
//! The fetch is made for a list that holds only one-hit HSPs of length 0 or less too, which the
//! re-evaluation deletes after it unread (blast_engine.c:1480-1497, blast_hits.c:687-699): with
//! `-matrix PAM30 -threshold 11 -window_size 0`, the query `TGGGATTGG` (W D W) against `bb`
//! (`TGGTGGTGGTAA`, W W W) gives only an HSP of length -1, whose header NCBI reads all the same
//! (line RG item X-short-tbx-link-2, probes on ws 2026-11-25, whose NCBI bytes the test compares).

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
/// Seq-ids when `parse_seqids`. The FASTA file is named `{name}.fna` in `dir`, so the database's
/// title (the `-in` argument) is that name.
fn makeblastdb(dir: &Path, name: &str, records: &[(&str, &str)], parse_seqids: bool) {
    let file = format!("{name}.fna");
    fasta(dir, &file, records);
    let mut args = vec!["makeblastdb", "-in", file.as_str()];
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

/// The refusal of a header that the re-evaluation's fetch cannot read, sequence 1 of database `dbD`
/// (after NCBI's warning of a seqid list in text).
fn assert_fetch_refused(output: &Output, case: &str) {
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{case}: {stderr}");
    let last = stderr.lines().last().unwrap_or_default();
    assert!(
        stderr.ends_with('\n')
            && last.starts_with(
                "Error: under the negative seqid list NCBI reads the definition lines of each \
                 sequence that its search aligns again, and those of sequence 1 cannot be read ("
            )
            && last.contains("dbD.nhr")
            && last.ends_with("); this database is not supported by LOSAT's TBLASTX")
            && stderr.matches("Error: ").count() == 1,
        "{case}: {stderr}"
    );
}

/// NCBI's `-outfmt 7` of `q8` and the database `dbD` before its engine error in the next query
/// batch (61 bytes; probe of line RG item X-short-tbx review 2 on ws, 2026-11-18).
const Q8_OUTFMT_7: &str = "# TBLASTX 2.17.0+\n# Query: q8\n# Database: dbD\n# 0 hits found\n";

/// NCBI's `-outfmt 0` of `q8` and the database `dbD` before its engine error in the next query
/// batch: the prolog and the report of `q8` with its footer, no epilog (478 bytes; the same probe).
const Q8_OUTFMT_0: &str = concat!(
    "TBLASTX 2.17.0+\n\n\n",
    "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.\n",
    "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.\n",
    "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of\n",
    "protein database search programs\", Nucleic Acids Res. 25:3389-3402.\n\n\n\n",
    "Database: dbD.fna\n",
    "           3 sequences; 56 total letters\n\n\n\n",
    "Query= q8\n\n",
    "Length=8\n\n\n",
    "***** No hits found *****\n\n\n\n",
    "Lambda      K        H\n",
    "   0.318    0.134    0.401 \n\n\n",
    "Effective search space used: 26\n",
);

/// A subject whose header the re-evaluation's fetch cannot read under a negative seqid list ends
/// the search where NCBI's does, with an engine error (exit 3; LOSAT's refusal), after the reports
/// of the query batches before (the Owner's decision of 2026-11-16 on line RG's Q-X2c). Database
/// `dbD` of `w12`, `nn30` (30 N) and `z14` with parsed Seq-ids, 8 bytes 0xFF written over the start
/// of `nn30`'s header (right after `w12`'s, whose length is the header file of a database of `w12`
/// alone), and a negative seqid list of `z14`: at `-threshold 1` the query `nq` leaves preliminary
/// HSPs on `nn30`, whose fetch reads its header. In one batch the search is refused (LOSAT wrote
/// the rows of `w12` before this version); with `BATCH_SIZE=1` and `q8` (no word) before `nq`, the
/// report of `q8` is written whole first, NCBI's bytes.
#[test]
fn a_header_the_fetch_cannot_read_fails_the_batch_after_the_batches_before() {
    let temp = TempDir::new("header");
    let dir = temp.0.as_path();
    let w12 = ("w12", "TGGTGGTGGTTT");
    makeblastdb(dir, "dbW", &[w12], true);
    let n30 = "N".repeat(30);
    let set = [w12, ("nn30", n30.as_str()), ("z14", "TGGTGGTGGTTTAA")];
    makeblastdb(dir, "dbD", &set, true);
    let first = std::fs::metadata(dir.join("dbW.nhr")).unwrap().len() as usize;
    let headers = dir.join("dbD.nhr");
    let mut bytes = std::fs::read(&headers).unwrap();
    assert!(bytes.len() > first + 8, "nn30's header follows w12's");
    bytes[first..first + 8].fill(0xFF);
    std::fs::write(&headers, bytes).unwrap();
    std::fs::write(dir.join("neg_z14.txt"), "z14\n").unwrap();
    let nq = fasta(dir, "nq.fna", &[("nq", "ATGCATTATTGGTGGTGGTTT")]);
    let q8nq = fasta(
        dir,
        "q8nq.fna",
        &[("q8", "ATGCATTA"), ("nq", "ATGCATTATTGGTGGTGGTTT")],
    );
    for threads in ["1", "4"] {
        let search = |query: &Path, batch_size, outfmt| {
            let mut args = vec!["tblastx", "-query", query.to_str().unwrap()];
            args.extend(["-db", "dbD", "-negative_seqidlist", "neg_z14.txt"]);
            args.extend(SEARCH);
            args.extend(["-threshold", "1"]);
            args.extend(["-outfmt", outfmt, "-num_threads", threads]);
            losat(dir, batch_size, &args)
        };
        let one_batch = search(&nq, None, "6");
        assert_fetch_refused(&one_batch, &format!("threads {threads}: one batch"));
        assert!(one_batch.stdout.is_empty(), "threads {threads}: one batch");
        for (outfmt, before) in [("6", ""), ("7", Q8_OUTFMT_7), ("0", Q8_OUTFMT_0)] {
            let case = format!("threads {threads}: q8 then nq, outfmt {outfmt}");
            let batches = search(&q8nq, Some("1"), outfmt);
            assert_fetch_refused(&batches, &case);
            // `before` is NCBI's frozen output: LOSAT's is NCBI's in LOSAT's words.
            assert_eq!(
                String::from_utf8_lossy(&batches.stdout),
                String::from_utf8_lossy(&losat::wording::losat_output_bytes(before.as_bytes())),
                "{case}"
            );
        }
    }
}

/// NCBI's `-outfmt 0` of `q8` and the database `dbD` of `w12`, `bb`, `z14` and `pad` with PAM30,
/// before its engine error in the next query batch: the prolog and the report of `q8` with its
/// footer, no epilog (480 bytes; probe of line RG item X-short-tbx-link-2 on ws, 2026-11-25).
const Q8_PAM30_OUTFMT_0: &str = concat!(
    "TBLASTX 2.17.0+\n\n\n",
    "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.\n",
    "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.\n",
    "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of\n",
    "protein database search programs\", Nucleic Acids Res. 25:3389-3402.\n\n\n\n",
    "Database: dbD.fna\n",
    "           4 sequences; 338 total letters\n\n\n\n",
    "Query= q8\n\n",
    "Length=8\n\n\n",
    "***** No hits found *****\n\n\n\n",
    "Lambda      K        H\n",
    "   0.340    0.283     1.75 \n\n\n",
    "Effective search space used: 218\n",
);

/// A subject whose list after the reap by the preliminary E-value holds only an HSP of length 0
/// or less still reaches the re-evaluation's fetch, which reads its header under a negative seqid
/// list: one that cannot be read ends the search there, as NCBI's does (exit 3; LOSAT's refusal),
/// with sum statistics and without, after the reports of the query batches before. Database
/// `dbD` of `w12`, `bb` (`TGGTGGTGGTAA`), `z14` and `pad` (300 bases without a word) with parsed
/// Seq-ids, 8 bytes 0xFF written over the start of `bb`'s header, a negative seqid list of
/// `z14`; with PAM30, `-threshold 11` and `-window_size 0`, the query `q` (`TGGGATTGG`) leaves
/// `bb` (and `w12`) only an HSP of length -1. At `-evalue 1` the reap removes it (its linked
/// E-value is between 1 and 10, NCBI's too: exit 0 at 1, the engine error at 10), so nothing
/// reads the header; the database undamaged gives no row. LOSAT wrote no row and exited 0 before
/// this version.
#[test]
fn a_list_of_only_alignments_of_length_zero_or_less_reaches_the_fetch() {
    let temp = TempDir::new("inverted");
    let dir = temp.0.as_path();
    let w12 = ("w12", "TGGTGGTGGTTT");
    makeblastdb(dir, "dbW", &[w12], true);
    let pad = "TAA".repeat(100);
    let set = [
        w12,
        ("bb", "TGGTGGTGGTAA"),
        ("z14", "TGGTGGTGGTTTAA"),
        ("pad", pad.as_str()),
    ];
    makeblastdb(dir, "dbD", &set, true);
    makeblastdb(dir, "dbU", &set, true);
    let first = std::fs::metadata(dir.join("dbW.nhr")).unwrap().len() as usize;
    let headers = dir.join("dbD.nhr");
    let mut bytes = std::fs::read(&headers).unwrap();
    assert!(bytes.len() > first + 8, "bb's header follows w12's");
    bytes[first..first + 8].fill(0xFF);
    std::fs::write(&headers, bytes).unwrap();
    std::fs::write(dir.join("neg_z14.txt"), "z14\n").unwrap();
    let q = fasta(dir, "q.fna", &[("q", "TGGGATTGG")]);
    let q8q = fasta(dir, "q8q.fna", &[("q8", "ATGCATTA"), ("q", "TGGGATTGG")]);
    for threads in ["1", "4"] {
        let search = |db: &str, query: &Path, batch_size, extra: &[&str]| {
            let mut args = vec!["tblastx", "-query", query.to_str().unwrap()];
            args.extend(["-db", db, "-negative_seqidlist", "neg_z14.txt"]);
            args.extend(["-window_size", "0", "-seg", "no"]);
            args.extend(["-matrix", "PAM30", "-threshold", "11"]);
            args.extend(extra);
            args.extend(["-num_threads", threads]);
            losat(dir, batch_size, &args)
        };
        let deep = ["-evalue", "1e10", "-outfmt", "6"];
        for sum_stats in ["true", "false"] {
            let case = format!("threads {threads}: -sum_stats {sum_stats}");
            let mut extra = deep.to_vec();
            extra.extend(["-sum_stats", sum_stats]);
            let refused = search("dbD", &q, None, &extra);
            assert_fetch_refused(&refused, &case);
            assert!(refused.stdout.is_empty(), "{case}");
        }
        let reaped = search("dbD", &q, None, &["-evalue", "1", "-outfmt", "6"]);
        assert_eq!(
            reaped.status.code(),
            Some(0),
            "threads {threads}: -evalue 1"
        );
        assert!(reaped.stdout.is_empty(), "threads {threads}: -evalue 1");
        let undamaged = search("dbU", &q, None, &deep);
        assert_eq!(
            undamaged.status.code(),
            Some(0),
            "threads {threads}: undamaged"
        );
        assert!(undamaged.stdout.is_empty(), "threads {threads}: undamaged");
        for (outfmt, before) in [("6", ""), ("7", Q8_OUTFMT_7), ("0", Q8_PAM30_OUTFMT_0)] {
            let case = format!("threads {threads}: q8 then q, outfmt {outfmt}");
            let batches = search(
                "dbD",
                &q8q,
                Some("1"),
                &["-evalue", "1e10", "-outfmt", outfmt],
            );
            assert_fetch_refused(&batches, &case);
            // `before` is NCBI's frozen output: LOSAT's is NCBI's in LOSAT's words.
            assert_eq!(
                String::from_utf8_lossy(&batches.stdout),
                String::from_utf8_lossy(&losat::wording::losat_output_bytes(before.as_bytes())),
                "{case}"
            );
        }
    }
}
