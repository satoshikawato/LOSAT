//! blastp `-db` searches with `LOSAT_X_BPREDODBSUBJ` (`src/algorithm/blastp/x_redo_subjects.rs`).
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless its variable is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. With it on, the composition redo reads the subject of each
//! match from the database when it redoes the match, instead of a copy of every subject of the
//! preliminary hit lists made before the redo, and the reports read the subjects of the final hit
//! lists again. These tests compare the reports with the switch off and on on every redo path (the
//! one-thread loop, `LOSAT_X_BPREDOPAR`, the two query-parallel paths of `LOSAT_X_BPREDOPAR=0`) and
//! for each set of residues the reports keep (tabular rows of the default fields and Seq-aligns:
//! the subjects with an O and those of the HSPs whose query has one; outfmt 0, XML and other
//! tabular fields: every final subject; the records of `LOSAT_X_BPLISTEDLEAN=0`: every listed
//! subject; the recount of the HSPs with a pyrrolysine on proteins made for it), and run its shadow mode
//! (`LOSAT_X_BPREDODBSUBJSHADOW`), which makes the copy too and compares every read with it; its
//! line at exit shows that the path ran. A difference stops the run.
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

/// A directory of this test.
fn work_dir() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-bpredodbsubj-{}", std::process::id()));
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
    let output = Command::new(env!("CARGO_BIN_EXE_LOSAT"))
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

/// The record BDT62584.1 of the fixture's query proteome (`MeenMJNV.faa`; more subjects of the
/// database than any other), written to `dir`: a search with one query, which the parallel redo of
/// one query's matches takes (`LOSAT_X_BPREDOPAR=0`).
fn one_query(dir: &Path) -> PathBuf {
    let text = std::fs::read_to_string(fixture("MeenMJNV.faa")).unwrap();
    let start = text
        .find(">BDT62584.1 ")
        .expect("the record is in the fixture");
    let end = text[start + 1..]
        .find("\n>")
        .map_or(text.len(), |at| start + 1 + at + 1);
    let path = dir.join("one_query.faa");
    std::fs::write(&path, &text[start..end]).unwrap();
    path
}

/// One blastp run of `query` against the database `db`, with `vars` as its only `LOSAT_*` variables
/// that choose a code path, at most 5 subjects per query (so that subjects of the preliminary lists
/// leave the final ones) and `-outfmt outfmt`; it must exit with 0. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test --all-features`, as in
/// CI) needs that file name as soon as a search has more than one thread, and it selects nothing.
fn blastp(query: &Path, db: &Path, threads: u32, vars: &[(&str, &str)], outfmt: &str) -> Run {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
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
    let dir = work_dir();
    let db = subject_database(&dir);
    let all = fixture("MeenMJNV.faa");
    let one = one_query(&dir);
    // Each redo path and each set of residues that the reports keep once; one search of all
    // queries (`LOSAT_X_BPQBLOCK=0`) and the block search of `LOSAT_X_BPQBLOCK` (its default)
    // alternate.
    for (query, outfmt, threads, vars) in [
        (&all, "6", 1, &[("LOSAT_X_BPQBLOCK", "0")][..]),
        (&all, "6", 4, &[][..]),
        (&all, "0", 2, &[("LOSAT_X_BPREDOPAR", "0")][..]),
        (
            &one,
            "6",
            2,
            &[("LOSAT_X_BPREDOPAR", "0"), ("LOSAT_X_BPQBLOCK", "0")][..],
        ),
        (&all, "6 qseqid sseqid slen stitle sseq", 1, &[][..]),
        (&all, "5", 2, &[("LOSAT_X_BPQBLOCK", "0")][..]),
        (&all, "8", 1, &[][..]),
        (&all, "6", 2, &[("LOSAT_X_BPLISTEDLEAN", "0")][..]),
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
        assert!(!off.stdout.is_empty(), "{what}: the fixture pair has rows");
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
    let dir = work_dir().join("pyrrolysine");
    std::fs::create_dir_all(&dir).unwrap();
    let with_o_at_128 = format!("{}O{}", &SUFFIX[..100], &SUFFIX[101..]);
    let subjects = format!(
        ">sx X at 27\n{PREFIX}X{SUFFIX}\n>so O at 128\n{PREFIX}G{with_o_at_128}\n\
         >s2 other\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA\n"
    );
    std::fs::write(dir.join("subjects.faa"), subjects).unwrap();
    let db = protein_database(&dir.join("subjects.faa"), &dir.join("subjects"));
    let queries = format!(">qo O at 27\n{PREFIX}O{SUFFIX}\n>qg G at 27\n{PREFIX}G{SUFFIX}\n");
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
            // The recount ran: qo's O against sx's X is a mismatch.
            let row = String::from_utf8_lossy(&off.stdout)
                .lines()
                .find(|row| row.starts_with("qo\tsx\t"))
                .map(str::to_string)
                .unwrap_or_else(|| panic!("{what}: no qo x sx row"));
            assert!(!row.contains("\t100.000\t"), "{what}: {row}");
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
