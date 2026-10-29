//! The scan on presence bytes in BLASTX's direct seed loop (`LOSAT_X_BXSCANBYTE`, in the strict
//! set; `src/algorithm/blastx/x_seed_scan.rs` module header, tblastn's
//! `src/algorithm/tblastn/x_scan_byte.rs` for the scan) on whole searches.
//!
//! NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127
//! ```c
//! for (s = s_first; s <= s_last; s++) {
//!     index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
//!     if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
//!         numhits = bbc[index].num_used;
//!         if (numhits <= (array_size - totalhits)) {
//!             ...
//!         } else {
//!             s_range[1] = (Int4)(s - subject->sequence);
//!             return totalhits;
//!         }
//!     }
//! }
//! ```
//! The switch's scan returns NCBI's offset pairs and scan range on every call (the unit test
//! `x_scan_byte_gives_ncbi_scan_calls_on_blastx_tables` compares the calls of both forms on
//! blastx's tables), so a search with it gives the report of the scans before it
//! (`LOSAT_X_BXSCANBYTE=0`). That reports equal NCBI's is the business of the BLASTX fixtures.
//! The switch is on by default: the AVX2 forms where the CPU has AVX2, the plain loops for sparse
//! tables (at most three tenths of the backbone cells with hits) on x86-64 without AVX2 and in
//! WebAssembly, nothing on aarch64; `=2` takes the plain loops for every table and `=3` for sparse
//! ones, on any CPU. It acts in the direct loop only, which on x86-64 with AVX2 scans the subjects
//! of batches with few hits per used cell (a short query) and elsewhere every batch;
//! `LOSAT_X_BXDIAG=0` makes it scan every subject here. `LOSAT_X_REACH` shows the side each run
//! took (the scan on presence bytes, or the scan before for a table that is not sparse under
//! `=3`).
//!
//! The queries are a 300-nt piece of AP027131 (a sparse table) and 6 pieces of 1,500 to 3,900
//! letters (two query batches, of 12,290 and 3,895 letters, whose tables are not sparse); the
//! subjects are the first 150 proteins of AP027131, so there are extensions and HSPs.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const SWITCH: &str = "LOSAT_X_BXSCANBYTE";
const OFF: (&str, &str) = (SWITCH, "0");
const PLAIN: (&str, &str) = (SWITCH, "2");
const SPARSE: (&str, &str) = (SWITCH, "3");
/// The direct loop on every subject (`LOSAT_X_BXDIAG`'s fused loop off).
const DIRECT: (&str, &str) = ("LOSAT_X_BXDIAG", "0");

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_blastx_scan_byte_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The records of a FASTA file as (title, letters).
fn records(name: &str) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(fixture(name)).expect("read fixture");
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(title) = line.strip_prefix('>') {
            records.push((title.to_string(), Vec::new()));
        } else if let Some(record) = records.last_mut() {
            record.1.extend_from_slice(line.trim().as_bytes());
        }
    }
    records
}

fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) {
    let mut text = Vec::new();
    for (title, letters) in records {
        text.extend_from_slice(format!(">{title}\n").as_bytes());
        for line in letters.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(path, text).expect("write FASTA");
}

/// The inputs of the searches, removed when dropped.
struct Inputs {
    dir: PathBuf,
    /// The 300-nt query (a sparse table).
    short: PathBuf,
    /// The six pieces (two batches, neither table sparse).
    long: PathBuf,
    subject: PathBuf,
    db: PathBuf,
}

impl Inputs {
    fn new() -> Self {
        let dir = temp_path("inputs");
        std::fs::create_dir_all(&dir).unwrap();
        let genome = &records("AP027131.fasta")[0].1;
        let pieces: Vec<(String, Vec<u8>)> = (0..6)
            .map(|index| {
                let start = index * 80_000;
                let length = 1_500 + (index * 479) % 2_400;
                (
                    format!("q{index} piece {start}"),
                    genome[start..start + length].to_vec(),
                )
            })
            .collect();
        let read = vec![("r300 piece 4000".to_string(), genome[4_000..4_300].to_vec())];
        let mut proteins = records("AP027131.faa");
        proteins.truncate(150);
        let short = dir.join("short.fna");
        let long = dir.join("long.fna");
        let subject = dir.join("subject.faa");
        write_fasta(&short, &read);
        write_fasta(&long, &pieces);
        write_fasta(&subject, &proteins);
        let db = dir.join("db");
        let made = losat(&[])
            .arg("makeblastdb")
            .arg("-in")
            .arg(&subject)
            .args(["-dbtype", "prot", "-out"])
            .arg(&db)
            .output()
            .expect("run losat makeblastdb");
        assert!(
            made.status.success(),
            "makeblastdb: {}",
            String::from_utf8_lossy(&made.stderr)
        );
        Inputs {
            dir,
            short,
            long,
            subject,
            db,
        }
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The binary with `env` as its only variables that choose a code path. `LOSAT_BLASTX_WORKER_LOG`
/// names a file of its own: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one thread,
/// and it selects nothing.
fn losat(env: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .envs(env.iter().copied());
    command
}

struct Run {
    stdout: Vec<u8>,
    stderr: String,
    code: Option<i32>,
    /// The sides of the switch that the run recorded (`LOSAT_X_REACH`): (new, ref).
    reached: (bool, bool),
}

impl Run {
    /// Standard error without the diagnostic lines (`[TIMING]`, `[X_SHADOW]`).
    fn messages(&self) -> String {
        self.stderr
            .lines()
            .filter(|line| !line.starts_with("[TIMING]") && !line.starts_with("[X_SHADOW]"))
            .map(|line| format!("{line}\n"))
            .collect()
    }
}

/// One `blastx` run with `env` and `LOSAT_X_REACH` set to a directory of its own (removed here).
fn blastx(args: &[&str], env: &[(&str, &str)]) -> Run {
    let dir = temp_path("reach");
    std::fs::create_dir_all(&dir).unwrap();
    let output = losat(env)
        .env("LOSAT_X_REACH", &dir)
        .arg("blastx")
        .args(args)
        .output()
        .expect("run LOSAT");
    let (mut new, mut reference) = (false, false);
    for file in std::fs::read_dir(&dir).expect("the reach directory") {
        let text = std::fs::read_to_string(file.unwrap().path()).unwrap();
        for row in text.lines().skip(1) {
            let f: Vec<&str> = row.split('\t').collect();
            if f[0] == SWITCH {
                new |= f[1] == "1";
                reference |= f[2] == "1";
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    Run {
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        code: output.status.code(),
        reached: (new, reference),
    }
}

/// The arguments of one search of `query`: `-subject` or `-db`, outfmt, threads.
fn search(inputs: &Inputs, query: &Path, db: bool, outfmt: &str, threads: &str) -> Vec<String> {
    let target = if db {
        ["-db", inputs.db.to_str().unwrap()]
    } else {
        ["-subject", inputs.subject.to_str().unwrap()]
    };
    [
        "-query",
        query.to_str().unwrap(),
        target[0],
        target[1],
        "-outfmt",
        outfmt,
        "-num_threads",
        threads,
    ]
    .map(String::from)
    .to_vec()
}

fn avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// The sides the default records for a run whose direct loop scans a table that is `sparse` or
/// not: every table where the CPU has AVX2, sparse ones on x86-64 without AVX2 and in
/// WebAssembly, none on aarch64 (the switch does nothing there).
fn default_sides(sparse: bool) -> (bool, bool) {
    if avx2() {
        (true, false)
    } else if cfg!(any(target_arch = "x86_64", target_arch = "wasm32")) {
        (sparse, !sparse)
    } else {
        (false, false)
    }
}

/// Every arm gives the report of `=0` (the scans before), on `-subject` and `-db`, at 1 and 3
/// threads; the plain loops and the sparse rule run where their rule says (`LOSAT_X_REACH`).
#[test]
fn the_scan_on_presence_bytes_gives_the_reports_of_the_scans_before() {
    let inputs = Inputs::new();
    for (query, sparse) in [(&inputs.short, true), (&inputs.long, false)] {
        for (db, outfmt, threads) in [(false, "6", "1"), (true, "7", "3")] {
            let args = search(&inputs, query, db, outfmt, threads);
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let what = format!(
                "{} db={db} outfmt={outfmt} threads={threads}",
                query.file_name().unwrap().to_string_lossy()
            );
            let reference = blastx(&args, &[DIRECT, OFF]);
            assert_eq!(reference.code, Some(0), "{what}: {}", reference.stderr);
            assert!(!reference.stdout.is_empty(), "{what}: no report");
            assert_eq!(reference.reached, (false, false), "{what}: =0 records");
            // (arm, run, the sides it must record; `None` = not checked).
            for (arm, run, sides) in [
                (
                    "direct",
                    blastx(&args, &[DIRECT]),
                    Some(default_sides(sparse)),
                ),
                (
                    "plain",
                    blastx(&args, &[DIRECT, PLAIN]),
                    Some((true, false)),
                ),
                (
                    "sparse",
                    blastx(&args, &[DIRECT, SPARSE]),
                    Some((sparse, !sparse)),
                ),
                // The default: the fused loop where it pays (AVX2), the direct loop elsewhere.
                ("default", blastx(&args, &[]), None),
            ] {
                assert_eq!(run.code, reference.code, "{what} {arm}: {}", run.stderr);
                assert!(
                    run.stdout == reference.stdout,
                    "{what} {arm}: report differs"
                );
                assert_eq!(run.messages(), reference.messages(), "{what} {arm}");
                if let Some(sides) = sides {
                    assert_eq!(run.reached, sides, "{what} {arm}: sides (new, ref)");
                }
            }
        }
    }
}

/// The short query's batch takes the direct loop by default, also where the CPU has AVX2 (few
/// hits per used cell), so the default scans it on presence bytes wherever the switch acts.
#[test]
fn the_default_scans_a_short_query_on_presence_bytes() {
    let inputs = Inputs::new();
    let args = search(&inputs, &inputs.short, true, "6", "1");
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    let reference = blastx(&args, &[OFF]);
    assert_eq!(reference.code, Some(0), "{}", reference.stderr);
    let run = blastx(&args, &[]);
    assert_eq!(run.code, Some(0), "{}", run.stderr);
    assert!(run.stdout == reference.stdout, "report differs");
    assert_eq!(run.reached, default_sides(true), "sides (new, ref)");
}
