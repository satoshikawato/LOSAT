//! TBLASTN's refusals of a step where NCBI BLAST+ 2.17.0's search of a query batch crashes, does
//! not finish or stops: NCBI searches each batch and prints its reports before it reads the next
//! (tblastn_app.cpp:275-302), so when the search of batch k stops it has printed the reports of the
//! batches before. LOSAT writes them whole, without the epilog (for the first batch, the prologs),
//! and then refuses (`blastinput::app::ErrorAfterBatchReports`, line RG item X-short-bat); nothing
//! of the refused batch or of the batches after it is written, at every `-num_threads`. The single
//! searches are rows `db.xs.tblastn.*` and `db.xsb.tblastn.*` of
//! `tests/db_search_manifest/tblastn.tsv`; here what those rows cannot hold: every format, the
//! batches searched ahead (side by side, on their own pools), `-db`, and the query chunk refusal.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275,289-291,300-302
//! ```c++
//!             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//! ...
//!                     results = lcl_blast.Run();
//! ...
//!                     ITERATE(CSearchResultSet, result, *results) {
//!                         formatter.PrintOneResultSet(**result, query);
//!                     }
//! ```
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

/// The refusal of a query chunk that NCBI would split again (it stops with a null reference).
const SPLIT_AGAIN: &str = "an OVERLAP_CHUNK_SIZE so close to CHUNK_SIZE that NCBI BLAST+ would \
split a query chunk again (it then stops with a null reference) is not supported by LOSAT's \
TBLASTN";

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
/// (`-outfmt 0`'s, the XML2 and JSON documents' first lines, `-outfmt 20`'s field names).
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
        "5",
        "6",
        "6 std",
        "7",
        "8",
        "10",
        "11",
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
                "0" | "1" => {
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
            assert_eq!(
                String::from_utf8_lossy(&refused.stdout),
                String::from_utf8_lossy(written),
                "{case}: the report of qa"
            );
            assert_eq!(
                String::from_utf8_lossy(&one_batch.stdout),
                String::from_utf8_lossy(prolog),
                "{case}: one batch, the prologs"
            );
        }
    }
}

/// Six batches, the fourth refused: the batches searched ahead (side by side at 2 threads, side by
/// side and on pools of their own at 4) hand the refused batch to the loop, which writes the
/// reports of the three before, as the three queries alone give them, and nothing of the fifth
/// and sixth.
#[test]
fn the_batches_searched_ahead_stop_at_the_earliest_refused_batch() {
    let temp = TempDir::new("ahead");
    let dir = temp.0.as_path();
    let inputs = inputs(dir);
    let subject = ["-subject", inputs.subject.to_str().unwrap()];
    let six = fasta(
        dir,
        "q6.faa",
        &[
            ("q1", QA),
            ("q2", QA),
            ("q3", QA),
            ("xs_q", "MHYWWW"),
            ("q5", QA),
            ("q6", QA),
        ],
    );
    let three = fasta(dir, "q3.faa", &[("q1", QA), ("q2", QA), ("q3", QA)]);
    let batch_one = [("BATCH_SIZE", "1")];
    for outfmt in ["0", "6"] {
        for threads in ["1", "2", "4"] {
            let case = format!("outfmt {outfmt} threads {threads}");
            let before = search(dir, &batch_one, &three, &subject, outfmt, threads);
            assert!(before.status.success(), "{case}: the three queries");
            let refused = search(dir, &batch_one, &six, &subject, outfmt, threads);
            assert_refused(&refused, X_SHORT, &case);
            let text = String::from_utf8_lossy(&before.stdout);
            let written = match outfmt {
                "0" => &before.stdout[..text.rfind(EPILOG).expect("the epilog")],
                _ => &before.stdout[..],
            };
            assert!(!written.is_empty(), "{case}");
            assert_eq!(
                String::from_utf8_lossy(&refused.stdout),
                String::from_utf8_lossy(written),
                "{case}: the reports of q1 to q3"
            );
        }
    }
}

/// `losat makeblastdb` of the nucleotide FASTA `set` as database `name` in `dir`.
fn makeblastdb(dir: &Path, name: &str, set: &Path) {
    let args = [
        "makeblastdb",
        "-in",
        set.to_str().unwrap(),
        "-dbtype",
        "nucl",
        "-out",
        name,
    ];
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
    makeblastdb(dir, "sdb", &inputs.subject);
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

/// The refusal of a query chunk that NCBI would split again (`CHUNK_SIZE=100
/// OVERLAP_CHUNK_SIZE=80`: `qa`, 31 residues, is one chunk; a query of 310 residues is split into
/// 15 chunks of 95, which NCBI would split again and stops with a null reference) comes after the
/// report of the batch before, as NCBI's error does. A negative `CHUNK_SIZE` that splits a batch,
/// where NCBI's wrapped chunk ranges either stop it or search with hits lost, is refused before
/// any byte.
#[test]
fn a_query_chunk_split_again_is_refused_after_the_reports_before() {
    let temp = TempDir::new("chunk");
    let dir = temp.0.as_path();
    let inputs = inputs(dir);
    let subject = ["-subject", inputs.subject.to_str().unwrap()];
    let long = QA.repeat(10);
    let queries = fasta(dir, "ql.faa", &[("qa", QA), ("ql", &long)]);
    for threads in ["1", "4"] {
        let case = format!("threads {threads}");
        let alone = search(dir, &[], &inputs.first, &subject, "6", threads);
        assert!(alone.status.success() && !alone.stdout.is_empty(), "{case}");
        let split = [
            ("BATCH_SIZE", "1"),
            ("CHUNK_SIZE", "100"),
            ("OVERLAP_CHUNK_SIZE", "80"),
        ];
        let refused = search(dir, &split, &queries, &subject, "6", threads);
        assert_refused(&refused, SPLIT_AGAIN, &case);
        assert_eq!(refused.stdout, alone.stdout, "{case}: the report of qa");
        let negative = [
            ("BATCH_SIZE", "1"),
            ("CHUNK_SIZE", "-100"),
            ("OVERLAP_CHUNK_SIZE", "-200"),
        ];
        let refused = search(dir, &negative, &queries, &subject, "6", threads);
        assert_eq!(refused.status.code(), Some(1), "{case}");
        assert!(refused.stdout.is_empty(), "{case}: nothing before");
    }
}

/// `-in_pssm` is one batch: its refusal comes after the `-outfmt 0` prolog, and with `-outfmt 6`
/// nothing is written (row `db.xs.tblastn.h11.rej`).
#[test]
fn a_pssm_query_is_refused_after_the_prolog() {
    let temp = TempDir::new("pssm");
    let dir = temp.0.as_path();
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let mut pssm = Vec::new();
    std::io::Read::read_to_end(
        &mut flate2::read::GzDecoder::new(
            std::fs::File::open(tests.join("fixtures/pssm/x2d_d50.asn.gz")).unwrap(),
        ),
        &mut pssm,
    )
    .unwrap();
    std::fs::write(dir.join("x2d_d50.asn"), pssm).unwrap();
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
                text.starts_with("PSITBLASTN 2.17.0+\n")
                    && text.ends_with(" total letters\n\n")
                    && !text.contains("Query="),
                "{case}: {text}"
            );
        } else {
            assert!(text.is_empty(), "{case}: {text}");
        }
    }
}
