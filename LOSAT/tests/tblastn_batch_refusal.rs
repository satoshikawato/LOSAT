//! TBLASTN's refusals of a step where NCBI BLAST+ 2.17.0's search of a query batch crashes, does
//! not finish or stops: NCBI searches each batch and prints its reports before it reads the next
//! (tblastn_app.cpp:275-302), so when the search of batch k stops it has printed the reports of the
//! batches before. LOSAT writes them whole, without the epilog (for the first batch, the prologs),
//! and then refuses (`blastinput::app::ErrorAfterBatchReports`, line RG item X-short-bat); nothing
//! of the refused batch or of the batches after it is written, at every `-num_threads`. The single
//! searches are rows `db.xs.tblastn.bat.*` of `tests/db_search_manifest/tblastn.tsv`; here what
//! those rows cannot hold: every format, the batches searched ahead (side by side, on their own
//! pools, one after another), `-db`, the query chunk refusals, and those of `-in_pssm` (one batch:
//! after the prolog).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275,289-291,300-302
//! The inputs: `qa` (31 residues without W, G, V, N, H, T or P, so no word of it scores 13 against
//! a frame of `w11`) hits `pa`, its back-translation (93 bases), and nothing of `w11`
//! (`TGGTGGTGGTT`); `xs_q` (`MHYWWW`) hits `w11` with `-window_size 0`, where the traceback reads
//! the translation NCBI never makes for a subject of fewer than 12 bases (PD-LOSAT-NCBI-DEFECTS
//! 1.46; NCBI dies with SIGSEGV). With `BATCH_SIZE=1` the query `qa xs_q` is two batches: NCBI has
//! formatted the report of `qa` when it dies (line RG X-short review 1: NCBI wrote the first
//! query's row before its SIGSEGV; what reaches its standard output beyond the tabular reports,
//! which its formatter flushes query by query, `CBlastTabularInfo::~CBlastTabularInfo`, is what its
//! output buffer had passed on).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The X-short refusal's message on standard error.
const X_SHORT: &str = "TBLASTN: the traceback of an alignment on a subject (or its \
-subject_loc range) of fewer than 12 bases, where NCBI BLAST+ reads the translation of the \
subject's frame, which it never makes for a subject this short, and crashes; LOSAT does not \
support it";

/// The refusal of a query chunk that NCBI would split again (the chunk's search dereferences the
/// lookup table its setup did not make).
const SPLIT_AGAIN: &str = "an OVERLAP_CHUNK_SIZE so close to CHUNK_SIZE that NCBI BLAST+ would \
split a query chunk again (it then stops with a null reference) is not supported by LOSAT's \
TBLASTN";

/// The refusal of a query chunk without a query (NCBI's `CQuerySplitter::Split` dereferences its
/// missing query vector, before any chunk is searched).
const EMPTY_CHUNK: &str = "a query chunk without a query is not supported by LOSAT's TBLASTN";

/// The refusal of a negative `CHUNK_SIZE` whose wrapped chunk ranges leave gaps between them
/// (NCBI searches on with the hits in the gaps lost: PD-LOSAT-NCBI-DEFECTS D13).
const NEGATIVE_GAPS: &str =
    "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this \
query batch with NCBI BLAST+'s chunk ranges wrapped, which is not supported by LOSAT's TBLASTN";

const SEARCH: [&str; 6] = ["-window_size", "0", "-seg", "no", "-evalue", "1e10"];

/// `qa`: no word of it scores 13 against a frame of `w11` (W W W, G G G, V V V; N H H, T T T,
/// P P P).
const QA: &str = "MKELLRSAFQDIYEKCMRQLAEFSKDIYRLC";

/// The start of `-outfmt 0`'s epilog (`PrintEpilog`) of a `-subject` search, after the last
/// query's report.
const EPILOG: &str = "\n\n  Database: User specified sequence set";

/// A temporary directory of this test, removed when it is dropped (also when an assertion fails).
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-tblastn-batch-refusal-{tag}-{}",
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

/// `pa`: `QA` back-translated, one codon per residue.
fn pa() -> String {
    QA.chars()
        .map(|residue| match residue {
            'M' => "ATG",
            'K' => "AAA",
            'E' => "GAA",
            'L' => "CTG",
            'R' => "CGT",
            'S' => "AGC",
            'A' => "GCG",
            'F' => "TTT",
            'Q' => "CAG",
            'D' => "GAT",
            'I' => "ATT",
            'Y' => "TAT",
            'C' => "TGC",
            other => panic!("no codon for {other}"),
        })
        .collect()
}

/// The inputs in `dir`: the subjects `pa w11` (`s.fna`), the queries `qa xs_q` (`qq.faa`) and `qa`
/// alone (`qa.faa`).
struct Inputs {
    subject: PathBuf,
    queries: PathBuf,
    first: PathBuf,
}

fn inputs(dir: &Path) -> Inputs {
    let pa = pa();
    Inputs {
        subject: fasta(dir, "s.fna", &[("pa", &pa), ("w11", "TGGTGGTGGTT")]),
        queries: fasta(dir, "qq.faa", &[("qa", QA), ("xs_q", "MHYWWW")]),
        first: fasta(dir, "qa.faa", &[("qa", QA)]),
    }
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

/// `tblastn -query QUERY -subject SUBJECT` with `SEARCH`, `-outfmt OUTFMT -num_threads THREADS`.
fn search(
    dir: &Path,
    env: &[(&str, &str)],
    query: &Path,
    subject: &[&str],
    outfmt: &str,
    threads: &str,
) -> Output {
    let mut args = vec!["tblastn", "-query", query.to_str().unwrap()];
    args.extend(subject);
    args.extend(SEARCH);
    args.extend(["-outfmt", outfmt, "-num_threads", threads]);
    losat(dir, env, &args)
}

/// With `BATCH_SIZE=1` the report of `qa`, the first of two query batches, is written before the
/// refusal of the second, as `qa` alone writes it without what comes after the last report (the
/// epilog of `-outfmt 0`, 1 and 3, the end of the XML document, the count line of `-outfmt 7`); the
/// tabular reports, the Seq-annots and the archive of the batch are whole. The XML2 and JSON
/// documents of a `-subject` search hold the first batch's report back (NCBI writes it with the
/// second), so only their prolog is written. In one batch the refusal comes after the prologs alone
/// (`-outfmt 0`'s, the XML2 and JSON documents' first lines, `-outfmt 20`'s field names). Not here:
/// `-outfmt 13` and `14`, which need `-out` and write each query's report to a file of its own,
/// and 17, which TBLASTN does not write.
#[test]
fn the_reports_of_the_query_batches_before_the_refusal_are_written() {
    let temp = TempDir::new("formats");
    let dir = temp.0.as_path();
    let inputs = inputs(dir);
    let subject = ["-subject", inputs.subject.to_str().unwrap()];
    let batch_one = [("BATCH_SIZE", "1")];
    for outfmt in [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "6 std",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "15",
        "16",
        "18",
        "20 qseqid sseqid evalue",
    ] {
        for threads in ["1", "4"] {
            let case = format!("outfmt {outfmt} threads {threads}");
            let alone = search(dir, &[], &inputs.first, &subject, outfmt, threads);
            assert!(
                alone.status.success() && !alone.stdout.is_empty(),
                "{case}: qa alone: {}",
                String::from_utf8_lossy(&alone.stderr)
            );
            let refused = search(dir, &batch_one, &inputs.queries, &subject, outfmt, threads);
            assert_refused(&refused, X_SHORT, &case);
            let one_batch = search(dir, &[], &inputs.queries, &subject, outfmt, threads);
            assert_refused(&one_batch, X_SHORT, &case);
            let text = String::from_utf8_lossy(&alone.stdout);
            let report = text.find("\n\nQuery= qa");
            let (written, prolog): (&[u8], &[u8]) = match outfmt {
                "0" | "1" | "2" | "3" | "4" => {
                    let report = report.expect("the report of qa");
                    let epilog = report + text[report..].find(EPILOG).expect("the epilog");
                    (&alone.stdout[..epilog], &alone.stdout[..report])
                }
                "5" => {
                    let end = text
                        .rfind("\n</BlastOutput_iterations>")
                        .expect("the end of the document");
                    (&alone.stdout[..end], b"")
                }
                "7" => {
                    let count = text.rfind("# BLAST processed").expect("the count line");
                    (&alone.stdout[..count], b"")
                }
                "15" | "16" => {
                    let start = if outfmt == "15" {
                        "{\n\"BlastOutput2\": [\n".len()
                    } else {
                        text.find("\n<BlastOutput2>").expect("the first report")
                    };
                    (&alone.stdout[..start], &alone.stdout[..start])
                }
                "20 qseqid sseqid evalue" => {
                    let names = text.find('\n').expect("the field names") + 1;
                    (&alone.stdout, &alone.stdout[..names])
                }
                _ => (&alone.stdout, b""),
            };
            assert_bytes(
                &refused.stdout,
                written,
                &format!("{case}: the report of qa"),
            );
            assert_bytes(
                &one_batch.stdout,
                prolog,
                &format!("{case}: one batch, the prologs"),
            );
        }
    }
}

/// The batches searched ahead hand the earliest refused batch to the loop, which writes the
/// reports of the batches before it, as those queries alone give them, and nothing of it or of the
/// batches after it. With two subjects, `x_search_batches` searches six batches side by side at 2
/// threads; at 4 threads four side by side and the last two on pools of their own, and of five
/// batches four side by side and the fifth after them (`one_after_another`): the refused batch is
/// met in each of those places, and of two refused batches the earlier one fails.
#[test]
fn the_batches_searched_ahead_stop_at_the_earliest_refused_batch() {
    let temp = TempDir::new("ahead");
    let dir = temp.0.as_path();
    let inputs = inputs(dir);
    let subject = ["-subject", inputs.subject.to_str().unwrap()];
    const X: &str = "MHYWWW";
    // (queries, the number before the first refused one, threads with the layout at 4 threads)
    let cases: [(&[&str], usize, &[&str], &str); 4] = [
        (
            &[QA, QA, QA, X, QA, QA],
            3,
            &["1", "2", "4"],
            "side_by_side=4 one_after_another=0 own_pools=2",
        ),
        (
            &[QA, QA, QA, QA, X, QA],
            4,
            &["4"],
            "side_by_side=4 one_after_another=0 own_pools=2",
        ),
        (
            &[QA, QA, QA, QA, X],
            4,
            &["4"],
            "side_by_side=4 one_after_another=1 own_pools=0",
        ),
        (
            &[QA, QA, QA, X, QA, X],
            3,
            &["2", "4"],
            "side_by_side=4 one_after_another=0 own_pools=2",
        ),
    ];
    for (index, (letters, before, threads_list, layout)) in cases.into_iter().enumerate() {
        let ids: Vec<String> = (1..=letters.len()).map(|n| format!("q{n}")).collect();
        let records: Vec<(&str, &str)> = ids
            .iter()
            .map(String::as_str)
            .zip(letters.iter().copied())
            .collect();
        let queries = fasta(dir, &format!("ahead{index}.faa"), &records);
        let first = fasta(dir, &format!("before{index}.faa"), &records[..before]);
        for outfmt in ["0", "6"] {
            for &threads in threads_list {
                let case = format!("case {index} outfmt {outfmt} threads {threads}");
                let batch_one = [("BATCH_SIZE", "1")];
                let alone = search(dir, &batch_one, &first, &subject, outfmt, threads);
                assert!(alone.status.success(), "{case}: the queries before");
                let refused = search(
                    dir,
                    &[("BATCH_SIZE", "1"), ("LOSAT_WASI_THREADS_DEBUG", "1")],
                    &queries,
                    &subject,
                    outfmt,
                    threads,
                );
                assert_refused(&refused, X_SHORT, &case);
                if threads == "4" {
                    let stderr = String::from_utf8_lossy(&refused.stderr);
                    let stage = stderr
                        .lines()
                        .find(|line| {
                            line.starts_with(
                                "[losat-thread-stage] program=tblastn stage=query_batches ",
                            )
                        })
                        .unwrap_or_else(|| panic!("{case}: no query_batches line: {stderr}"));
                    assert!(stage.ends_with(layout), "{case}: {stage}");
                }
                let text = String::from_utf8_lossy(&alone.stdout);
                let written = match outfmt {
                    "0" => &alone.stdout[..text.rfind(EPILOG).expect("the epilog")],
                    _ => &alone.stdout[..],
                };
                assert!(!written.is_empty(), "{case}");
                assert_bytes(
                    &refused.stdout,
                    written,
                    &format!("{case}: the reports of the queries before"),
                );
            }
        }
    }
}

/// `losat makeblastdb` of the nucleotide FASTA `set` as database `name` in `dir`, with `extra`.
fn makeblastdb(dir: &Path, name: &str, set: &Path, extra: &[&str]) {
    let mut args = vec![
        "makeblastdb",
        "-in",
        set.to_str().unwrap(),
        "-dbtype",
        "nucl",
        "-out",
        name,
    ];
    args.extend(extra);
    let made = losat(dir, &[], &args);
    assert!(
        made.status.success(),
        "makeblastdb {name}: {}",
        String::from_utf8_lossy(&made.stderr)
    );
}

/// The same with `-db`: the report of `qa` before the refusal.
#[test]
fn a_database_search_writes_the_reports_before_the_refusal() {
    let temp = TempDir::new("db");
    let dir = temp.0.as_path();
    let inputs = inputs(dir);
    makeblastdb(dir, "sdb", &inputs.subject, &[]);
    let db = ["-db", "sdb"];
    for outfmt in ["0", "6"] {
        for threads in ["1", "4"] {
            let case = format!("outfmt {outfmt} threads {threads}");
            let alone = search(dir, &[], &inputs.first, &db, outfmt, threads);
            assert!(alone.status.success(), "{case}: qa alone");
            let refused = search(
                dir,
                &[("BATCH_SIZE", "1")],
                &inputs.queries,
                &db,
                outfmt,
                threads,
            );
            assert_refused(&refused, X_SHORT, &case);
            let text = String::from_utf8_lossy(&alone.stdout);
            let written = match outfmt {
                "0" => &alone.stdout[..text.rfind("\n\n  Database: ").expect("the epilog")],
                _ => &alone.stdout[..],
            };
            assert!(!written.is_empty(), "{case}");
            assert_eq!(
                String::from_utf8_lossy(&refused.stdout),
                String::from_utf8_lossy(written),
                "{case}: the report of qa"
            );
        }
    }
}

/// The refusals of a query chunk where NCBI stops in the search of the batch come after the report
/// of the batch before, as NCBI's error does (`qa`, 31 residues, and `q9`, 9, are not split), in
/// NCBI's order: `CQuerySplitter::Split` stops at a chunk without a query before any chunk is
/// searched, the search of a chunk that would be split again stops after the chunks before it.
/// - `CHUNK_SIZE=100 OVERLAP_CHUNK_SIZE=80` split a query of 310 residues into 15 chunks of 95 to
///   100, which NCBI would split again (chunk 0 stops);
/// - `CHUNK_SIZE=23 OVERLAP_CHUNK_SIZE=5` split it into 16 chunks of 22 and a last one of 38,
///   which NCBI would split again (after the 16 chunks before it are searched), also where that
///   chunk's letters are all X (no valid context: NCBI's stop comes first);
/// - `CHUNK_SIZE=-100 OVERLAP_CHUNK_SIZE=-200` (NCBI's `size_t` ranges wrap) make the same query a
///   first chunk that covers it and 2 chunks without a query, and a query of 200 residues a chunk
///   of 1 residue and a chunk without a query;
/// - `CHUNK_SIZE=100 OVERLAP_CHUNK_SIZE=95` make `qa` a first chunk that covers it (the evened
///   chunk size, 85, is below the overlap) and 5 chunks without a query.
///
/// A query of 248 residues is 2 chunks of 25 and 23 residues with 200 between them, which NCBI
/// searches on with the hits there lost: refused where the search of its first chunk with a valid
/// context would read the subjects (after that chunk's setup checks and the reads of the database
/// that fail before it), after the report of the query before (PD-LOSAT-NCBI-DEFECTS D13, the form
/// of 1.40; row `db.xs.tblastn.bat.gap.6`).
#[test]
fn the_query_chunk_stops_are_refused_after_the_reports_before() {
    let temp = TempDir::new("chunk");
    let dir = temp.0.as_path();
    let inputs = inputs(dir);
    let subject = ["-subject", inputs.subject.to_str().unwrap()];
    let long = QA.repeat(10);
    let tail_x = format!("{QA}M{}", "X".repeat(278));
    let queries = fasta(dir, "ql.faa", &[("qa", QA), ("ql", &long)]);
    let tail = fasta(dir, "qt.faa", &[("qa", QA), ("qt", &tail_x)]);
    let empty = fasta(dir, "qe.faa", &[("qa", QA), ("qe", &QA.repeat(7)[..200])]);
    let gaps = fasta(dir, "qg.faa", &[("qa", QA), ("qg", &QA.repeat(8))]);
    let short = fasta(dir, "q9.faa", &[("q9", &QA[..9])]);
    let close = fasta(dir, "q9qa.faa", &[("q9", &QA[..9]), ("qa", QA)]);
    let positive = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "100"),
        ("OVERLAP_CHUNK_SIZE", "80"),
    ];
    let last = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "23"),
        ("OVERLAP_CHUNK_SIZE", "5"),
    ];
    let negative = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-100"),
        ("OVERLAP_CHUNK_SIZE", "-200"),
    ];
    let overlap = [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "100"),
        ("OVERLAP_CHUNK_SIZE", "95"),
    ];
    for threads in ["1", "4"] {
        let alone = search(dir, &[], &inputs.first, &subject, "6", threads);
        assert!(
            alone.status.success() && !alone.stdout.is_empty(),
            "threads {threads}"
        );
        let short_alone = search(dir, &[], &short, &subject, "6", threads);
        assert!(short_alone.status.success(), "threads {threads}: q9");
        for (name, env, query, message, before) in [
            ("split again", &positive, &queries, SPLIT_AGAIN, &alone),
            (
                "the last chunk split again",
                &last,
                &queries,
                SPLIT_AGAIN,
                &alone,
            ),
            (
                "the last chunk, all X, split again",
                &last,
                &tail,
                SPLIT_AGAIN,
                &alone,
            ),
            (
                "negative, a first chunk that covers the batch, then chunks without a query",
                &negative,
                &queries,
                EMPTY_CHUNK,
                &alone,
            ),
            (
                "negative, a chunk without a query",
                &negative,
                &empty,
                EMPTY_CHUNK,
                &alone,
            ),
            (
                "an evened chunk size below the overlap, chunks without a query",
                &overlap,
                &close,
                EMPTY_CHUNK,
                &short_alone,
            ),
            (
                "negative, gaps between the chunks",
                &negative,
                &gaps,
                NEGATIVE_GAPS,
                &alone,
            ),
        ] {
            let case = format!("{name}, threads {threads}");
            let refused = search(dir, env, query, &subject, "6", threads);
            assert_refused(&refused, message, &case);
            assert_bytes(
                &refused.stdout,
                &before.stdout,
                &format!("{case}: the report of the query before"),
            );
        }
    }
}

/// The splitting values of `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` above: chunk 0 split again
/// (`100`, `80`), the last chunk split again (`23`, `5`), chunks without a query (`-100`,
/// `-200`), each with `BATCH_SIZE=1`.
const SPLITS: [[(&str, &str); 3]; 3] = [
    [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "100"),
        ("OVERLAP_CHUNK_SIZE", "80"),
    ],
    [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "23"),
        ("OVERLAP_CHUNK_SIZE", "5"),
    ],
    [
        ("BATCH_SIZE", "1"),
        ("CHUNK_SIZE", "-100"),
        ("OVERLAP_CHUNK_SIZE", "-200"),
    ],
];

/// A database that a seqid list leaves without a subject (`-negative_seqidlist` of every id,
/// `-seqidlist` of an id it does not have; NCBI posts `Filtering resulted in an empty database.`
/// with each batch's results) or whose volumes have none (`e_nA`, an alias with `NSEQ 1` and
/// `LENGTH 100`, on which NCBI reads no sequence): NCBI 2.17.0 still splits a batch and stops in
/// its search, after the report of the batch before (`qa`, no hit, and its warning, which the
/// stopped batch's stderr starts with); nothing of the stopped batch is written, its warning
/// neither (probes of line RG item X-short-bat, fix round 3). Without a stop (a negative pair
/// whose chunks leave gaps between them, which lose no hit where there is no subject) the batch
/// is reported as it is without the variables.
#[test]
fn a_split_batch_stops_on_a_database_without_a_subject() {
    let temp = TempDir::new("nosubject");
    let dir = temp.0.as_path();
    let inputs = inputs(dir);
    makeblastdb(dir, "pdb", &inputs.subject, &["-parse_seqids"]);
    std::fs::write(dir.join("all.txt"), "pa\nw11\n").unwrap();
    std::fs::write(dir.join("none.txt"), "zz\n").unwrap();
    let e_na = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/aliaskeys/e_nA");
    let queries = fasta(dir, "ql.faa", &[("qa", QA), ("ql", &QA.repeat(10))]);
    let gaps = fasta(dir, "qg.faa", &[("qa", QA), ("qg", &QA.repeat(8))]);
    let databases: [&[&str]; 3] = [
        &["-db", "pdb", "-negative_seqidlist", "all.txt"],
        &["-db", "pdb", "-seqidlist", "none.txt"],
        &["-db", e_na.to_str().unwrap()],
    ];
    for db in databases {
        for threads in ["1", "4"] {
            let alone = search(dir, &[], &inputs.first, db, "7", threads);
            let text = String::from_utf8_lossy(&alone.stdout);
            assert!(alone.status.success(), "{db:?} threads {threads}: {text}");
            let before = &alone.stdout[..text.rfind("# BLAST processed").expect("the count line")];
            for (env, message) in SPLITS.iter().zip([SPLIT_AGAIN, SPLIT_AGAIN, EMPTY_CHUNK]) {
                let case = format!("{db:?} {env:?} threads {threads}");
                let refused = search(dir, env, &queries, db, "7", threads);
                assert_refused(&refused, message, &case);
                assert_bytes(
                    &refused.stdout,
                    before,
                    &format!("{case}: the report of qa"),
                );
                let stderr = String::from_utf8_lossy(&refused.stderr);
                assert!(!stderr.contains("Query_2"), "{case}: {stderr}");
                assert!(
                    refused.stderr.starts_with(&alone.stderr),
                    "{case}: the warning of qa first: {stderr}"
                );
            }
            let case = format!("{db:?} gaps threads {threads}");
            let whole = search(dir, &SPLITS[2], &gaps, db, "7", threads);
            let unsplit = search(dir, &SPLITS[2][..1], &gaps, db, "7", threads);
            assert!(whole.status.success() && unsplit.status.success(), "{case}");
            assert_bytes(&whole.stdout, &unsplit.stdout, &case);
            assert_bytes(&whole.stderr, &unsplit.stderr, &format!("{case}: stderr"));
        }
    }
}

/// Where NCBI's search first fails to read the database, a split batch meets NCBI's chunk stops
/// in NCBI's order (2.17.0 probes of line RG item X-short-bat, fix round 3 and review 3): on
/// `n_min1` (`MIN_SEQ_LENGTH 1`) the setup of the first chunk with a valid context that is not
/// split again fails (`The average subject length is too short`; its cases never reach the
/// `break` of `ReadFails::Search`, since the chunk's own setup checks fail first), on `e_nF` (volumes without an
/// OID, `NSEQ 1`) that of the first chunk not split again, which reads the length of OID 0 before
/// its Karlin-Altschul check (`Database memory map file error`, `ReadFails::Setup`), and on
/// `arch_tbn_snucl_nonsq` (a database without its sequence file) the search of the first chunk
/// with a valid context, which reads its first sequence after that check (`Database memory map
/// file error`, `ReadFails::Search`: an all-X head goes on to the split-again stop, which `e_nF`
/// meets later, and `late` and `gaps` end at the read, which only the `break` after the first
/// valid chunk's checks gives); a chunk without a query stops NCBI before all of them, a chunk
/// split again before the chunks after it. One batch (a query of 310 residues: `ql`, `qa` ten
/// times; `xh`, 278 X then `qa` and M, whose chunks 0 to 15 of the second split are all X; `xt`,
/// the reverse; `qg`, `qa` eight times, whose negative pair leaves gaps between its chunks), so
/// nothing is written.
#[test]
fn a_split_batch_meets_a_failing_read_of_the_database_where_ncbi_does() {
    let temp = TempDir::new("readfails");
    let dir = temp.0.as_path();
    let databases = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let x = "X".repeat(278);
    let long = fasta(dir, "ql.faa", &[("ql", &QA.repeat(10))]);
    let gaps = fasta(dir, "qg.faa", &[("qg", &QA.repeat(8))]);
    let head = fasta(dir, "xh.faa", &[("xh", &format!("{x}{QA}M"))]);
    let tail = fasta(dir, "xt.faa", &[("xt", &format!("{QA}M{x}"))]);
    let min1 = "LOSAT engine error: The average subject length is too short";
    let memory = "LOSAT engine error: Database memory map file error";
    // (database, splitting values, query, the refusal or (with exit status 3) NCBI's error)
    let nosq = "arch_tbn_snucl_nonsq";
    let cases: [(&str, usize, &PathBuf, Result<&str, &str>); 15] = [
        ("aliaskeys/n_min1", 0, &long, Ok(SPLIT_AGAIN)),
        ("aliaskeys/n_min1", 1, &long, Err(min1)),
        ("aliaskeys/n_min1", 2, &long, Ok(EMPTY_CHUNK)),
        ("aliaskeys/n_min1", 1, &head, Ok(SPLIT_AGAIN)),
        ("aliaskeys/n_min1", 1, &tail, Err(min1)),
        ("aliaskeys/e_nF", 0, &long, Ok(SPLIT_AGAIN)),
        ("aliaskeys/e_nF", 1, &long, Err(memory)),
        ("aliaskeys/e_nF", 2, &long, Ok(EMPTY_CHUNK)),
        ("aliaskeys/e_nF", 1, &head, Err(memory)),
        (nosq, 0, &long, Ok(SPLIT_AGAIN)),
        (nosq, 1, &long, Err(memory)),
        (nosq, 2, &long, Ok(EMPTY_CHUNK)),
        (nosq, 1, &head, Ok(SPLIT_AGAIN)),
        (nosq, 1, &tail, Err(memory)),
        (nosq, 2, &gaps, Err(memory)),
    ];
    for (name, split, query, expected) in cases {
        let db = databases.join(name);
        for threads in ["1", "4"] {
            let case = format!("{name} {:?} {query:?} threads {threads}", SPLITS[split]);
            let run = search(
                dir,
                &SPLITS[split],
                query,
                &["-db", db.to_str().unwrap()],
                "6",
                threads,
            );
            assert!(run.stdout.is_empty(), "{case}");
            match expected {
                Ok(message) => assert_refused(&run, message, &case),
                Err(error) => {
                    let stderr = String::from_utf8_lossy(&run.stderr);
                    assert_eq!(run.status.code(), Some(3), "{case}: {stderr}");
                    assert_eq!(stderr, format!("{error}\n"), "{case}");
                }
            }
        }
    }
}

/// A batch that LOSAT refuses for two queries of one id under `-parse_deflines` (`q1`, `q2`, `q1`
/// of 300 residues each, `tests/db_search/q1c/avgdup_bp.faa`; line WQA D10) gives NCBI's error of
/// the batch's set-up where NCBI's search of the batch meets it, also where NCBI splits the batch
/// (review 4 of line RG item X-short-bat, M1; rows `db.pd.q1c.tn.avgdup.split*.subj.6` and
/// `db.pd.q1c.tn.nsm1.split.6`). NCBI splits the batch by all three records, and its chunks meet
/// their stops and the setup's checks in the order of the chunk loop:
/// - `CHUNK_SIZE=300 OVERLAP_CHUNK_SIZE=10`: 3 chunks of 307 residues, none split again; chunk 0's
///   setup fails for subjects too short for `-searchsp 1000` (`tests/db_search/q1c/avg_s.fna`),
///   and for the hit-saving parameters of `dbcounts/nsm1` (`NSEQ 4294967295`) with `-searchsp`;
/// - `403` and `202`: 4 chunks, none split again, where `q1` and `q2` alone would make a chunk
///   that NCBI splits again: the error;
/// - `453` and `227`: NCBI's chunk 0 would be split again (NCBI stops there, exit status 3), where
///   `q1` and `q2` alone would make 2 chunks that are not: LOSAT's refusal of the repeated id,
///   which stands for a stop of NCBI's chunks;
/// - unsplit: the error, as before.
#[test]
fn a_refused_batch_that_ncbi_splits_meets_the_set_up_errors_where_ncbi_does() {
    let temp = TempDir::new("refusedsplit");
    let dir = temp.0.as_path();
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let query = tests.join("q1c/avgdup_bp.faa");
    let subject = tests.join("q1c/avg_s.fna");
    let nsm1 = tests.join("dbcounts/nsm1");
    let short = "LOSAT engine error: The average subject length is too short";
    let hit_saving = "LOSAT engine error: Unknown error code 1";
    let repeated = "two query records have the Seq-id lcl|q1";
    let subjects = ["-subject", subject.to_str().unwrap(), "-searchsp", "1000"];
    let database = ["-db", nsm1.to_str().unwrap(), "-searchsp", "1000000"];
    // (subjects, CHUNK_SIZE and OVERLAP_CHUNK_SIZE, NCBI's error (exit status 3) or the refusal)
    let cases: [(&[&str], Option<(&str, &str)>, Result<&str, &str>); 5] = [
        (&subjects, Some(("300", "10")), Err(short)),
        (&subjects, Some(("403", "202")), Err(short)),
        (&subjects, Some(("453", "227")), Ok(repeated)),
        (&subjects, None, Err(short)),
        (&database, Some(("300", "10")), Err(hit_saving)),
    ];
    for (searched, split, expected) in cases {
        let env: Vec<(&str, &str)> = split
            .map(|(chunk, overlap)| vec![("CHUNK_SIZE", chunk), ("OVERLAP_CHUNK_SIZE", overlap)])
            .unwrap_or_default();
        for threads in ["1", "4"] {
            let case = format!("{searched:?} {split:?} threads {threads}");
            let mut args = vec!["tblastn", "-query", query.to_str().unwrap()];
            args.extend(searched);
            args.extend(["-parse_deflines", "-outfmt", "6", "-num_threads", threads]);
            let run = losat(dir, &env, &args);
            let stderr = String::from_utf8_lossy(&run.stderr);
            assert!(run.stdout.is_empty(), "{case}: {stderr}");
            match expected {
                Ok(refusal) => {
                    assert_eq!(run.status.code(), Some(1), "{case}: {stderr}");
                    assert!(stderr.contains(refusal), "{case}: {stderr}");
                }
                Err(error) => {
                    // After the warnings of the subjects without letters (`s2`, `s3`).
                    assert_eq!(run.status.code(), Some(3), "{case}: {stderr}");
                    assert!(
                        stderr.ends_with(&format!("{error}\n"))
                            && stderr.matches(error).count() == 1,
                        "{case}: {stderr}"
                    );
                }
            }
        }
    }
}

/// The setup of a query batch reads the length of OID 0 of a database whose statistics' total is
/// not positive (`tests/db_search/aliaskeys/e_nF`: volumes without an OID, `NSEQ 1`) before its
/// setup checks, so that read's error comes first also where the average subject length fails
/// them without the Gumbel block (`-searchsp`, `OLD_FSC`): for a batch searched whole, for one
/// that NCBI splits and for one that LOSAT refuses (`q1`, `q2`, `q1` under `-parse_deflines`).
/// NCBI 2.17.0: rows `db.ak.tbn.e0nseq.sp.6`, `db.ak.tbn.e0nseq.fsc.6`, `db.pd.q1c.tn.enf.sp.6`
/// and `db.pd.q1c.tn.enf.split.6` (line RG item X-short-bat, fix 5).
#[test]
fn the_read_of_oid_0_comes_before_the_set_up_checks() {
    let temp = TempDir::new("oid0");
    let dir = temp.0.as_path();
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search");
    let db = tests.join("aliaskeys/e_nF");
    let one = tests.join("taxfilter/q.faa");
    let repeated = tests.join("q1c/avgdup_bp.faa");
    let memory = "LOSAT engine error: Database memory map file error\n";
    let split = [("CHUNK_SIZE", "300"), ("OVERLAP_CHUNK_SIZE", "10")];
    let searchsp = ["-searchsp", "1000"];
    let deflines = ["-parse_deflines", "-searchsp", "1000"];
    // (query, environment, options)
    let cases: [(&PathBuf, &[(&str, &str)], &[&str]); 4] = [
        (&one, &[], &searchsp),
        (&one, &[("OLD_FSC", "1")], &[]),
        (&repeated, &[], &deflines),
        (&repeated, &split, &deflines),
    ];
    for (query, env, options) in cases {
        for threads in ["1", "4"] {
            let case = format!("{query:?} {env:?} {options:?} threads {threads}");
            let mut args = vec!["tblastn", "-query", query.to_str().unwrap()];
            args.extend(["-db", db.to_str().unwrap()]);
            args.extend(options);
            args.extend(["-outfmt", "6", "-num_threads", threads]);
            let run = losat(dir, env, &args);
            let stderr = String::from_utf8_lossy(&run.stderr);
            assert!(run.stdout.is_empty(), "{case}: {stderr}");
            assert_eq!(run.status.code(), Some(3), "{case}: {stderr}");
            assert_eq!(stderr, memory, "{case}");
        }
    }
}

/// The refusal of a PSSM query whose letters give no Karlin-Altschul parameters (one that `-seg`
/// masks completely), with which NCBI crashes in the setup of its preliminary search.
const MASKED_PSSM: &str = "a PSSM query whose letters give no Karlin-Altschul parameters (one \
that -seg masks completely), with which NCBI BLAST+ 2.17.0's tblastn crashes, is not supported by \
LOSAT's TBLASTN";

/// The PSSM `tests/fixtures/pssm/NAME.asn.gz` written to `dir/NAME.asn`.
fn pssm_file(dir: &Path, name: &str) -> PathBuf {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let mut pssm = Vec::new();
    std::io::Read::read_to_end(
        &mut flate2::read::GzDecoder::new(
            std::fs::File::open(tests.join(format!("fixtures/pssm/{name}.asn.gz"))).unwrap(),
        ),
        &mut pssm,
    )
    .unwrap();
    let path = dir.join(format!("{name}.asn"));
    std::fs::write(&path, pssm).unwrap();
    path
}

/// `-in_pssm` is one batch: its refusal comes after the `-outfmt 0` prolog, and with `-outfmt 6`
/// nothing is written (row `db.xs.tblastn.h11.rej`).
#[test]
fn a_pssm_query_is_refused_after_the_prolog() {
    let temp = TempDir::new("pssm");
    let dir = temp.0.as_path();
    pssm_file(dir, "x2d_d50");
    let subject = fasta(dir, "h11.fna", &[("h11", "ATGCATTATTT")]);
    for (outfmt, threads) in [("0", "1"), ("0", "4"), ("6", "1")] {
        let case = format!("outfmt {outfmt} threads {threads}");
        let mut args = vec!["tblastn", "-in_pssm", "x2d_d50.asn"];
        args.extend(["-subject", subject.to_str().unwrap()]);
        args.extend(SEARCH);
        args.extend(["-outfmt", outfmt, "-num_threads", threads]);
        let refused = losat(dir, &[], &args);
        assert_refused(&refused, X_SHORT, &case);
        let text = String::from_utf8_lossy(&refused.stdout);
        if outfmt == "0" {
            assert!(
                text.starts_with(concat!(
                    "PSITBLASTN (LOSAT ",
                    env!("CARGO_PKG_VERSION"),
                    ")\n"
                )) && text.ends_with(" total letters\n\n")
                    && !text.contains("Query="),
                "{case}: {text}"
            );
        } else {
            assert!(text.is_empty(), "{case}: {text}");
        }
    }
}

/// A PSSM query whose letters give no Karlin-Altschul parameters (`ps7_gs_masked`, which the
/// default SEG masks completely) on the database `tslice`, as the probe `pssmmasked.f0` of line RG
/// X-short-bat: NCBI crashes in the setup of its preliminary search, after its prolog, so LOSAT
/// writes the `-outfmt 0` prolog (that of the same search with `-seg no`, which NCBI runs, up to its
/// query's report) and then refuses; with `-outfmt 6` nothing is written (row
/// `db.tbn.inpssm.rej.masked`).
#[test]
fn a_masked_pssm_query_is_refused_after_the_prolog() {
    let temp = TempDir::new("pssm-masked");
    let dir = temp.0.as_path();
    let pssm = pssm_file(dir, "ps7_gs_masked");
    let db = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/tslice");
    let run = |extra: &[&str], outfmt: &str, threads: &str| {
        let mut args = vec!["tblastn", "-in_pssm", pssm.to_str().unwrap()];
        args.extend(["-db", db.to_str().unwrap(), "-evalue", "1000"]);
        args.extend(extra);
        args.extend(["-outfmt", outfmt, "-num_threads", threads]);
        losat(dir, &[], &args)
    };
    let searched = run(&["-seg", "no"], "0", "1");
    let text = String::from_utf8_lossy(&searched.stdout);
    assert!(searched.status.success(), "-seg no: {text}");
    let prolog = &searched.stdout[..text.find("\n\nQuery=").expect("the query's report")];
    assert!(
        prolog.starts_with(
            concat!("PSITBLASTN (LOSAT ", env!("CARGO_PKG_VERSION"), ")\n").as_bytes()
        ) && prolog.ends_with(b" total letters\n\n"),
        "-seg no: {text}"
    );
    for outfmt in ["0", "6"] {
        for threads in ["1", "4"] {
            let case = format!("outfmt {outfmt} threads {threads}");
            let refused = run(&[], outfmt, threads);
            assert_refused(&refused, MASKED_PSSM, &case);
            let expected: &[u8] = if outfmt == "0" { prolog } else { b"" };
            assert_bytes(&refused.stdout, expected, &case);
        }
    }
}
