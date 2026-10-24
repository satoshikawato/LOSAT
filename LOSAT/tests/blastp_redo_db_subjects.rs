//! blastp `-db` searches with `LOSAT_X_BPREDODBSUBJ` (`src/algorithm/blastp/x_redo_subjects.rs`).
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless its variable is `0`.
//! With it on, the composition redo reads the subject of each
//! match from the database when it redoes the match, instead of a copy of every subject of the
//! preliminary hit lists made before the redo, and the reports read the subjects of the final hit
//! lists again. These tests compare the reports with the switch off and on on every redo path (the
//! one-thread loop, `LOSAT_X_BPREDOPAR`, the query-parallel and the match-parallel redo of
//! `LOSAT_X_BPREDOPAR=0`) and for each set of residues the reports keep (tabular rows of the
//! default fields and Seq-aligns: the subjects with an O and those of the HSPs whose query has one;
//! outfmt 0 and other tabular fields (XML reads the same set): every final subject; the records of
//! `LOSAT_X_BPLISTEDLEAN=0`: every listed subject; the recount of the HSPs with a pyrrolysine on
//! proteins made for it), and run its shadow mode (`LOSAT_X_BPREDODBSUBJSHADOW`), which makes the
//! copy too and compares every read with it; its line at exit shows that the path ran. A
//! difference stops the run.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1381-1385
//! ```c
//!         if( program_number == eBlastTypeTblastn ) {
//!             seq_info->seq_arg.encoding = eBlastEncodingNcbi4na;
//!         } else {
//!             seq_info->seq_arg.encoding = eBlastEncodingProtein;
//!         }
//! ```

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// What one run printed: the report and standard error.
struct Run {
    stdout: Vec<u8>,
    stderr: String,
}

const SHADOW_PREFIX: &str = "[X_SHADOW] BPREDODBSUBJ ";

impl Run {
    /// The counters of the `[X_SHADOW] BPREDODBSUBJ` line, which must end with `identical`.
    fn shadow(&self) -> HashMap<String, u64> {
        let line = self
            .stderr
            .lines()
            .find(|line| line.starts_with(SHADOW_PREFIX))
            .unwrap_or_else(|| panic!("no {SHADOW_PREFIX}line in {}", self.stderr));
        assert!(line.ends_with(" identical"), "{line}");
        line[SHADOW_PREFIX.len()..]
            .split(' ')
            .filter_map(|field| field.split_once('='))
            .map(|(key, value)| (key.to_string(), value.parse().expect("a count")))
            .collect()
    }

    /// Standard error without the shadow's line.
    fn stderr_without_shadow(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with(SHADOW_PREFIX))
            .map(|line| format!("{line}\n"))
            .collect()
    }
}

/// A directory of the test `test` (the tests of this file run in parallel).
fn work_dir(test: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-bpredodbsubj-{}-{test}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A protein database of the fixture's subject proteome (`MejoMJNV.faa`, 124 sequences), made with
/// `LOSAT makedb` in `dir`; returns the database name.
fn subject_database(dir: &Path) -> PathBuf {
    protein_database(&fixture("MejoMJNV.faa"), &dir.join("MejoMJNV"))
}

/// A protein database of the records of `fasta`, made with `LOSAT makedb` as `name`; returns
/// `name`.
fn protein_database(fasta: &Path, name: &Path) -> PathBuf {
    let name = name.to_path_buf();
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .arg("makedb")
        .arg("-in")
        .arg(fasta)
        .args(["-dbtype", "prot"])
        .arg("-out")
        .arg(&name)
        .output()
        .expect("run LOSAT makedb");
    assert!(
        output.status.success(),
        "makedb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    name
}

/// `count` records of the fixture's query proteome (`MeenMJNV.faa`) from BDT62584.1 on (more
/// subjects of the database than any other: 12 in its preliminary list), written to `dir`. One is a
/// search with one query, whose matches the redo of `LOSAT_X_BPREDOPAR=0` parallelizes; two (with
/// BDT62585.1, which has subjects too) take its query-parallel redo.
fn queries(dir: &Path, count: usize) -> PathBuf {
    let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let start = text
        .find(">BDT62584.1 ")
        .expect("the record is in the fixture");
    let mut end = start;
    for _ in 0..count {
        end = text[end + 1..]
            .find("\n>")
            .map_or(text.len(), |at| end + 1 + at + 1);
    }
    let path = dir.join(format!("queries_{count}.faa"));
    std::fs::write(&path, &text[start..end]).unwrap();
    path
}

/// One blastp run of `query` against the database `db`, with `vars` as its only `LOSAT_*` variables
/// that choose a code path, at most 5 subjects per query (so that subjects of the preliminary lists
/// leave the final ones) and `-outfmt outfmt`; it must exit with 0. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in
/// CI) needs that file name as soon as a search has more than one thread, and it selects nothing.
fn blastp(query: &Path, db: &Path, threads: u32, vars: &[(&str, &str)], outfmt: &str) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || text == "BATCH_SIZE"
            || text == "CHUNK_SIZE"
            || text == "OVERLAP_CHUNK_SIZE"
        {
            command.env_remove(name);
        }
    }
    if std::env::var_os("LOSAT_BLASTX_WORKER_LOG").is_none() {
        command.env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::temp_dir().join("losat-bpredodbsubj-worker.log"),
        );
    }
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-db")
        .arg(db)
        .args([
            "-evalue",
            "0.001",
            "-max_target_seqs",
            "5",
            "-outfmt",
            outfmt,
        ])
        .args(["-num_threads", &threads.to_string()])
        .envs(vars.iter().copied())
        .output()
        .expect("run LOSAT");
    let run = Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    };
    assert!(
        output.status.success(),
        "{vars:?} threads {threads} outfmt {outfmt}: {:?} {}",
        output.status.code(),
        run.stderr
    );
    run
}

#[test]
fn reports_with_reads_per_match_are_the_copys_reports() {
    let dir = work_dir("reports");
    let db = subject_database(&dir);
    let one = queries(&dir, 1);
    let two = queries(&dir, 2);
    // Each redo path and each set of residues that the reports keep once. With
    // `LOSAT_WASI_THREADS_DEBUG=1` the `[losat-thread-stage]` lines name the paths: kappa-redo-par
    // takes precedence over kappa-match-redo and kappa-query-redo, which are printed before it. The
    // switch lists and reads the subjects search by search, the same for a block of the block
    // search (the default; every other `-db` test of blastp runs the switch in it) and for the one
    // search of all queries of `LOSAT_X_BPQBLOCK=0`, so no case needs more than two queries.
    for (query, outfmt, threads, vars) in [
        // `LOSAT_X_BPREDOPAR` (kappa-redo-par) over two queries of one batch, the default fields.
        (&two, "6", 4, &[][..]),
        // The query-parallel redo of `LOSAT_X_BPREDOPAR=0` (kappa-query-redo), outfmt 0.
        (&two, "0", 2, &[("LOSAT_X_BPREDOPAR", "0")][..]),
        // Its redo of one query's matches (kappa-match-redo), Seq-aligns.
        (&one, "8", 2, &[("LOSAT_X_BPREDOPAR", "0")][..]),
        // The one-thread loop, other tabular fields.
        (&one, "6 qseqid sseqid slen stitle sseq", 1, &[][..]),
        // `LOSAT_X_BPREDOPAR`, the records of `LOSAT_X_BPLISTEDLEAN=0`.
        (&one, "6", 2, &[("LOSAT_X_BPLISTEDLEAN", "0")][..]),
    ] {
        let with = |extra: (&'static str, &'static str)| {
            vars.iter().copied().chain([extra]).collect::<Vec<_>>()
        };
        let what = format!("{query:?} outfmt {outfmt} threads {threads} {vars:?}");
        let off = blastp(
            query,
            &db,
            threads,
            &with(("LOSAT_X_BPREDODBSUBJ", "0")),
            outfmt,
        );
        assert!(!off.stdout.is_empty(), "{what}: an empty report");
        let on = blastp(
            query,
            &db,
            threads,
            &with(("LOSAT_X_BPREDODBSUBJ", "1")),
            outfmt,
        );
        assert!(on.stdout == off.stdout, "{what}: report");
        assert_eq!(on.stderr, off.stderr, "{what}: standard error");
        let shadow = blastp(
            query,
            &db,
            threads,
            &with(("LOSAT_X_BPREDODBSUBJSHADOW", "1")),
            outfmt,
        );
        assert!(shadow.stdout == off.stdout, "{what}: shadow report");
        assert_eq!(shadow.stderr_without_shadow(), off.stderr, "{what}");
        let counts = shadow.shadow();
        assert!(counts["searches"] > 0, "{what}: {counts:?}");
        assert!(counts["lengths"] > 0, "{what}: {counts:?}");
        assert!(
            counts["reads"] > 0 && counts["residues"] > counts["reads"],
            "{what}: {counts:?}"
        );
    }
}

/// The letters before position 27 of the proteins of `pyrrolysine_recounts_are_the_copys`.
const PREFIX: &str = "MKTAYIAKQRQISFVKSHFSRQLEDN";
/// The letters after position 27.
const SUFFIX: &str = "GLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAACVRESVPALL";

// The reports recount the identities of the HSPs whose query or subject has a pyrrolysine (O), which
// the search reads as X: a query O against a subject X is not identical in NCBI's reports. With the
// tabular rows of the default fields the reports keep the residues of the subjects with an O and of
// the subjects of the HSPs whose query has one; these tests compare those rows (and those of other
// formats) with the switch off and on.
#[test]
fn pyrrolysine_recounts_are_the_copys() {
    let dir = work_dir("pyrrolysine");
    let with_o_at_128 = format!("{}O{}", &SUFFIX[..100], &SUFFIX[101..]);
    let subjects = format!(
        ">sx X at 27\n{PREFIX}X{SUFFIX}\n>so O at 128\n{PREFIX}G{with_o_at_128}\n\
         >s2 other\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA\n"
    );
    std::fs::write(dir.join("subjects.faa"), subjects).unwrap();
    let db = protein_database(&dir.join("subjects.faa"), &dir.join("subjects"));
    // qx holds an X where so holds its O: the recount of qx x so reads so's letters, which only a
    // read of the redo that found the O keeps.
    let with_x_at_128 = format!("{}X{}", &SUFFIX[..100], &SUFFIX[101..]);
    let queries = format!(
        ">qo O at 27\n{PREFIX}O{SUFFIX}\n>qg G at 27\n{PREFIX}G{SUFFIX}\n>qx X at 128\n{PREFIX}G{with_x_at_128}\n"
    );
    let query = dir.join("queries.faa");
    std::fs::write(&query, queries).unwrap();
    for (outfmt, threads) in [
        ("6", 1),
        ("6", 2),
        ("7", 1),
        ("6 qseqid sseqid pident nident mismatch positive", 2),
        ("0", 1),
    ] {
        let what = format!("outfmt {outfmt} threads {threads}");
        let off = blastp(
            &query,
            &db,
            threads,
            &[("LOSAT_X_BPREDODBSUBJ", "0")],
            outfmt,
        );
        if outfmt == "6" {
            // The recount ran: qo's O against sx's X and qx's X against so's O are mismatches.
            for pair in ["qo\tsx\t", "qx\tso\t"] {
                let row = String::from_utf8_lossy(&off.stdout)
                    .lines()
                    .find(|row| row.starts_with(pair))
                    .map(str::to_string)
                    .unwrap_or_else(|| panic!("{what}: no {pair:?} row"));
                assert!(!row.contains("\t100.000\t"), "{what}: {row}");
            }
        }
        let on = blastp(
            &query,
            &db,
            threads,
            &[("LOSAT_X_BPREDODBSUBJ", "1")],
            outfmt,
        );
        assert!(
            on.stdout == off.stdout,
            "{what}: report\n{}\n{}",
            String::from_utf8_lossy(&on.stdout),
            String::from_utf8_lossy(&off.stdout)
        );
        assert_eq!(on.stderr, off.stderr, "{what}: standard error");
        let shadow = blastp(
            &query,
            &db,
            threads,
            &[("LOSAT_X_BPREDODBSUBJSHADOW", "1")],
            outfmt,
        );
        assert!(shadow.stdout == off.stdout, "{what}: shadow report");
        assert_eq!(shadow.stderr_without_shadow(), off.stderr, "{what}");
        assert!(shadow.shadow()["reads"] > 0, "{what}");
    }
}
