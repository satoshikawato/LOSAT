//! BLASTN's refusals of a step where NCBI BLAST+ 2.17.0's search of a query batch stops: NCBI
//! searches each batch and prints its reports before it reads the next (blastn_app.cpp:276-312),
//! so when the search of batch k stops it has printed the reports of the batches before. LOSAT
//! writes them whole, without the epilog (for the first batch, the prologs), and then refuses
//! (`blastinput::app::ErrorAfterBatchReports`, line RG item X-short-bat-2); nothing of the refused
//! batch or of the batches after it is written, at every `-num_threads`, with the reports streamed
//! batch by batch (`LOSAT_X_RPTSTREAM`, on by default) and without. The single searches are rows
//! `db.xs.blastn.bat.*` of `tests/db_search_manifest/blastn.tsv`, which hold NCBI's bytes; here
//! every format of the rows, both report paths, `-db`, and the first batch.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:276-312
//! The inputs: `na` (93 bases) hits the subject `pa`, the same letters; `nl` is `na` ten times (930
//! bases). With `BATCH_SIZE=1` the query `na nl` is two batches, and the chunk sizes below split
//! `nl` alone (`na` is shorter than one chunk): a negative `CHUNK_SIZE` above a negative
//! `OVERLAP_CHUNK_SIZE` whose first chunk covers the batch (chunks without a query, where NCBI's
//! query splitter stops before any chunk is searched), one whose chunks leave gaps between them
//! (NCBI searches on and loses the hits in the gaps: decision D13 of PD-LOSAT-NCBI-DEFECTS), and a
//! positive pair whose evened chunk size falls below the overlap (a chunk without a query).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The refusal of a query chunk without a query (NCBI's `CQuerySplitter::Split` dereferences its
/// missing query vector, before any chunk is searched).
const EMPTY_CHUNK: &str = "a query chunk without a query is not supported by LOSAT's BLASTN";

/// The refusal of a negative `CHUNK_SIZE` whose wrapped chunk ranges leave gaps between them.
const NEGATIVE_GAPS: &str =
    "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this \
query batch with NCBI BLAST+'s chunk ranges wrapped, which is not supported by LOSAT's BLASTN";

/// `na`: the back-translation of TBLASTN's test query `qa`, one codon per residue.
const NA: &str = "ATGAAAGAACTGCTGCGTAGCGCGTTTCAGGATATTTATGAAAAATGCATGCGTCAGCTGGCGGAATTTAGCAAAGA\
TATTTATCGTCTGTGC";

/// The start of `-outfmt 0`'s epilog (`PrintEpilog`) of a `-subject` search, after the last
/// query's report.
const SUBJECT_EPILOG: &str = "\n\n  Database: User specified sequence set";

/// The start of `-outfmt 0`'s epilog of a `-db` search.
const DB_EPILOG: &str = "\n\n  Database: ";

/// The chunk sizes of the three stops, with `BATCH_SIZE=1`: (name, `CHUNK_SIZE`,
/// `OVERLAP_CHUNK_SIZE`, the refusal, whether one batch of `na nl` is refused too).
const STOPS: [(&str, &str, &str, &str, bool); 3] = [
    (
        "negative, a first chunk that covers the batch, then chunks without a query",
        "-100",
        "-200",
        EMPTY_CHUNK,
        true,
    ),
    (
        "negative, gaps between the chunks",
        "-1",
        "-200",
        NEGATIVE_GAPS,
        true,
    ),
    (
        "an evened chunk size below the overlap, a chunk without a query",
        "1400",
        "1000",
        EMPTY_CHUNK,
        // One batch of 1023 bases is split into two chunks that both hold a query.
        false,
    ),
];

/// The formats of the rows: the pairwise report, BLAST XML (never streamed), the tabular reports
/// and the archive (one Blast4 archive per batch).
const FORMATS: [&str; 5] = ["0", "5", "6", "7", "11"];

/// A temporary directory of this test, removed when it is dropped (also when an assertion fails).
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-blastn-batch-refusal-{tag}-{}",
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

/// What a search refused in the batch after `alone`'s writes of `alone` (the same search of the
/// first batch's queries alone, which NCBI runs to its end) in `outfmt`: all of it but what comes
/// after the last report (the epilog of `-outfmt 0`, from `epilog`; the end of the XML document;
/// the count line of `-outfmt 7`). The tabular reports and the archive are whole.
fn before_the_epilog<'a>(alone: &'a [u8], outfmt: &str, epilog: &str) -> &'a [u8] {
    let text = String::from_utf8_lossy(alone);
    match outfmt {
        "0" => &alone[..text.rfind(epilog).expect("the epilog")],
        "5" => {
            &alone[..text
                .rfind("\n</BlastOutput_iterations>")
                .expect("the end of the document")]
        }
        "7" => &alone[..text.rfind("# BLAST processed").expect("the count line")],
        _ => alone,
    }
}

/// What a search refused in its first batch writes in `outfmt`, from `alone` (the same search of a
/// query that is not refused): the prologs alone (`-outfmt 0`'s, up to the first report).
fn prologs<'a>(alone: &'a [u8], outfmt: &str) -> &'a [u8] {
    match outfmt {
        "0" => {
            let text = String::from_utf8_lossy(alone);
            &alone[..text.find("\n\nQuery= ").expect("the first report")]
        }
        _ => b"",
    }
}

/// The report paths: streamed batch by batch (the default) and after the loop.
const PATHS: [&[(&str, &str)]; 2] = [&[], &[("LOSAT_X_RPTSTREAM", "0")]];

/// With `BATCH_SIZE=1` the report of `na`, the first of two query batches, is written before the
/// refusal of the second, as `na` alone writes it without what comes after the last report; in one
/// batch the refusal comes after the prologs alone. `-subject`, every format of the rows, both
/// report paths, 1 and 4 threads.
#[test]
fn the_query_chunk_stops_are_refused_after_the_reports_before() {
    let temp = TempDir::new("chunk");
    let dir = temp.0.as_path();
    let subject_file = fasta(dir, "s.fna", &[("pa", NA), ("w11", "TGGTGGTGGTT")]);
    let first = fasta(dir, "na.fna", &[("na", NA)]);
    let long = NA.repeat(10);
    let queries = fasta(dir, "nl.fna", &[("na", NA), ("nl", &long)]);
    let subject = ["-subject", subject_file.to_str().unwrap()];
    for outfmt in FORMATS {
        for threads in ["1", "4"] {
            let search = |env: &[(&str, &str)], query: &Path| {
                let mut args = vec!["blastn", "-query", query.to_str().unwrap()];
                args.extend(subject);
                args.extend(["-outfmt", outfmt, "-num_threads", threads]);
                losat(dir, env, &args)
            };
            let alone = search(&[], &first);
            assert!(
                alone.status.success() && !alone.stdout.is_empty(),
                "outfmt {outfmt} threads {threads}: na alone: {}",
                String::from_utf8_lossy(&alone.stderr)
            );
            for (name, chunk, overlap, message, one_batch_refused) in STOPS {
                for path in PATHS {
                    let case = format!("{name}, outfmt {outfmt}, threads {threads}, {path:?}");
                    let sizes = [("CHUNK_SIZE", chunk), ("OVERLAP_CHUNK_SIZE", overlap)];
                    let mut env = path.to_vec();
                    env.extend(sizes);
                    let one_batch = search(&env, &queries);
                    env.push(("BATCH_SIZE", "1"));
                    let refused = search(&env, &queries);
                    assert_refused(&refused, message, &case);
                    assert_bytes(
                        &refused.stdout,
                        before_the_epilog(&alone.stdout, outfmt, SUBJECT_EPILOG),
                        &format!("{case}: the report of na"),
                    );
                    if one_batch_refused {
                        assert_refused(&one_batch, message, &case);
                        assert_bytes(
                            &one_batch.stdout,
                            prologs(&alone.stdout, outfmt),
                            &format!("{case}: one batch, the prologs"),
                        );
                    } else {
                        assert!(one_batch.status.success(), "{case}: one batch");
                    }
                }
            }
        }
    }
}

/// `losat makeblastdb` of the nucleotide FASTA `records` as database `name` in `dir`, with parsed
/// Seq-ids. The FASTA file is named `{name}.fna` in `dir`, so the database's title (the `-in`
/// argument) is that name.
fn makeblastdb(dir: &Path, name: &str, records: &[(&str, &str)]) {
    let file = format!("{name}.fna");
    fasta(dir, &file, records);
    let made = losat(
        dir,
        &[],
        &[
            "makeblastdb",
            "-in",
            &file,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            name,
        ],
    );
    assert!(
        made.status.success(),
        "makeblastdb {name}: {}",
        String::from_utf8_lossy(&made.stderr)
    );
}

/// The query chunk stops of a `-db` search (`pa w11`): the report of `na` before the refusal.
#[test]
fn a_database_search_writes_the_reports_before_the_chunk_stops() {
    let temp = TempDir::new("db");
    let dir = temp.0.as_path();
    makeblastdb(dir, "s", &[("pa", NA), ("w11", "TGGTGGTGGTT")]);
    let first = fasta(dir, "na.fna", &[("na", NA)]);
    let long = NA.repeat(10);
    let queries = fasta(dir, "nl.fna", &[("na", NA), ("nl", &long)]);
    for outfmt in ["0", "6", "7"] {
        for threads in ["1", "4"] {
            let search = |env: &[(&str, &str)], query: &Path| {
                let args = [
                    "blastn",
                    "-query",
                    query.to_str().unwrap(),
                    "-db",
                    "s",
                    "-outfmt",
                    outfmt,
                    "-num_threads",
                    threads,
                ];
                losat(dir, env, &args)
            };
            let alone = search(&[], &first);
            assert!(alone.status.success(), "outfmt {outfmt} threads {threads}");
            for (name, chunk, overlap, message, _) in STOPS {
                let case = format!("{name}, -db, outfmt {outfmt}, threads {threads}");
                let env = [
                    ("BATCH_SIZE", "1"),
                    ("CHUNK_SIZE", chunk),
                    ("OVERLAP_CHUNK_SIZE", overlap),
                ];
                let refused = search(&env, &queries);
                assert_refused(&refused, message, &case);
                assert_bytes(
                    &refused.stdout,
                    before_the_epilog(&alone.stdout, outfmt, DB_EPILOG),
                    &format!("{case}: the report of na"),
                );
            }
        }
    }
}

/// 60 bases of a fixed sequence of pseudo-random letters (a linear congruential generator from
/// `seed`), which no other of them hits at megablast's word size.
fn random_bases(seed: u64) -> String {
    let mut state = seed;
    (0..60)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[(state >> 62) as usize] as char
        })
        .collect()
}

/// The refusal of a header that the fetch under the negative seqid list cannot read, sequence 1 of
/// database `dbD` (after NCBI's warning of a seqid list in text).
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
            && last.ends_with("); this database is not supported by LOSAT's BLASTN")
            && stderr.matches("Error: ").count() == 1,
        "{case}: {stderr}"
    );
}

/// A subject whose header the fetch under a negative seqid list cannot read ends the search where
/// NCBI's does (its `BLAST engine error`, exit status 3; LOSAT's refusal), after the reports of the
/// query batches before, in the gapped search (the traceback's fetch) and with `-ungapped` (the
/// re-evaluation's): database `dbD` of `s1 s2 s3` with parsed Seq-ids, 8 bytes 0xFF written over
/// the start of `s2`'s header (right after `s1`'s, whose length is the header file of a database of
/// `s1` alone), and a negative seqid list of `s3`. `qa` (the letters of `s1`) hits `s1` alone; `qb`
/// (those of `s2`) hits `s2`, whose fetch reads its header. With `BATCH_SIZE=1` the report of `qa`
/// is written first, as `qa` alone writes it without what comes after the last report; in one batch
/// the refusal comes after the prologs.
#[test]
fn a_header_the_fetch_cannot_read_fails_the_batch_after_the_batches_before() {
    let temp = TempDir::new("header");
    let dir = temp.0.as_path();
    let (s1, s2, s3) = (random_bases(1), random_bases(2), random_bases(3));
    makeblastdb(dir, "dbW", &[("s1", &s1)]);
    makeblastdb(dir, "dbD", &[("s1", &s1), ("s2", &s2), ("s3", &s3)]);
    let first = std::fs::metadata(dir.join("dbW.nhr")).unwrap().len() as usize;
    let headers = dir.join("dbD.nhr");
    let mut bytes = std::fs::read(&headers).unwrap();
    assert!(bytes.len() > first + 8, "s2's header follows s1's");
    bytes[first..first + 8].fill(0xFF);
    std::fs::write(&headers, bytes).unwrap();
    std::fs::write(dir.join("neg_s3.txt"), "s3\n").unwrap();
    let qa = fasta(dir, "qa.fna", &[("qa", &s1)]);
    let qab = fasta(dir, "qab.fna", &[("qa", &s1), ("qb", &s2)]);
    let qb = fasta(dir, "qb.fna", &[("qb", &s2)]);
    for ungapped in [false, true] {
        for outfmt in FORMATS {
            for threads in ["1", "4"] {
                let search = |env: &[(&str, &str)], query: &Path| {
                    let mut args = vec!["blastn", "-query", query.to_str().unwrap()];
                    args.extend(["-db", "dbD", "-negative_seqidlist", "neg_s3.txt"]);
                    if ungapped {
                        args.push("-ungapped");
                    }
                    args.extend(["-outfmt", outfmt, "-num_threads", threads]);
                    losat(dir, env, &args)
                };
                let alone = search(&[], &qa);
                assert!(
                    alone.status.success() && !alone.stdout.is_empty(),
                    "ungapped {ungapped} outfmt {outfmt} threads {threads}: qa alone: {}",
                    String::from_utf8_lossy(&alone.stderr)
                );
                for path in PATHS {
                    let case = format!(
                        "ungapped {ungapped}, outfmt {outfmt}, threads {threads}, {path:?}"
                    );
                    let one_batch = search(path, &qb);
                    assert_fetch_refused(&one_batch, &format!("{case}: qb alone"));
                    assert_bytes(
                        &one_batch.stdout,
                        prologs(&alone.stdout, outfmt),
                        &format!("{case}: qb alone, the prologs"),
                    );
                    let mut env = path.to_vec();
                    env.push(("BATCH_SIZE", "1"));
                    let refused = search(&env, &qab);
                    assert_fetch_refused(&refused, &case);
                    assert_bytes(
                        &refused.stdout,
                        before_the_epilog(&alone.stdout, outfmt, DB_EPILOG),
                        &format!("{case}: the report of qa"),
                    );
                }
            }
        }
    }
}
