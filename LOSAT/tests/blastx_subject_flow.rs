//! BLASTX's subject loop searched as a flow of tasks (`LOSAT_X_BXSUBJFLOW`,
//! `src/algorithm/blastx/x_subject_flow.rs`) on whole searches.
//!
//! The switch is in the strict set (on unless `=0`). With `=0` and two or
//! more threads a query chunk's subjects are searched in waves of `threads x 64` and each wave is
//! handed to the sink after its slowest worker; with the switch the other threads search tasks of
//! `k` consecutive subjects ahead and the loop's thread hands them to the sink in subject order.
//! So the property under test is that the report does not depend on it: standard output equals
//! that of the same search with `=0`, for two, three and four threads, several `k` (64, the
//! default, and 2, 7 and 16), outfmt 0, 6 and 7, a database (`-db`) and a
//! FASTA subject set (`-subject`), and an ungapped search (whose subjects carry link cutoffs).
//! The shadow (`LOSAT_X_BXSUBJFLOWSHADOW`) searches every task again on another table and must
//! find the same subjects. The tasks are cut by residues as well (`LOSAT_X_BXTASKRES`, in the
//! strict set): a task ends after `k` subjects or before the subject that would take it past `k`
//! times the mean subject length; with `=0` they are the grid of `k` subjects. `LOSAT_TIMING`
//! shows that the flow ran (and the waves did not). The fallback to the waves where a diagonal offset could pass 2^30 (`-window_size` near 2^28, a
//! diagonal table of 1 GiB or more) is a unit test of the module. That the bytes equal NCBI's
//! is the business of the BLASTX fixtures and the campaign's comparison runs, but for the frozen
//! `-db_hard_mask` rows of `ov_m1`, whose bytes depend on what a worker's table holds when it
//! takes its next subject: they are run here with the waves and with the grid of 64.
//!
//! The query is 9,000 letters of AP027131 (one batch, one chunk, so the pool searches its
//! subjects); the subjects are the proteins of AP027131, AP027132 and AP027133 (1,711 records).

#![cfg(feature = "parallel")]

use std::path::{Path, PathBuf};
use std::process::Command;

mod temp_support;
use temp_support::InLosatWords;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with no `LOSAT_*` variable that chooses a code path. `LOSAT_BLASTX_WORKER_LOG`
/// stays: a binary built with the `blastx-worker-probe` feature (`cargo test
/// --all-features`, as in CI) needs that file name as soon as a search has more than one
/// thread, and it selects nothing.
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// The temporary directory of this test process.
fn temp_dir() -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-blastx-subject-flow-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The query (letters 100,000-109,000 of AP027131) and the subject set (its proteins and those
/// of AP027132 and AP027133) as files, and the subject set's `losat makeblastdb` database.
fn inputs(dir: &Path) -> (PathBuf, PathBuf, String) {
    inputs_of(dir, 9_000, 1)
}

/// A query of `length` letters of AP027131 from position 100,000, and the three protein sets
/// (1,711 proteins) `copies` times, each copy after the first with its names prefixed by `cN_`, as
/// a FASTA file and as a database.
fn inputs_of(dir: &Path, length: usize, copies: usize) -> (PathBuf, PathBuf, String) {
    let genome = std::fs::read_to_string(fixture("AP027131.fasta")).unwrap();
    let letters: String = genome
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(str::chars)
        .filter(|c| !c.is_ascii_whitespace())
        .skip(100_000)
        .take(length)
        .collect();
    let query = dir.join("query.fna");
    std::fs::write(
        &query,
        format!(">piece AP027131 100000-{}\n{letters}\n", 100_000 + length),
    )
    .unwrap();
    let set = dir.join("set.faa");
    let mut text = String::new();
    for copy in 0..copies {
        for name in ["AP027131.faa", "AP027132.faa", "AP027133.faa"] {
            let proteins = std::fs::read_to_string(fixture(name)).unwrap();
            if copy == 0 {
                text.push_str(&proteins);
            } else {
                text.push_str(&proteins.replace('>', &format!(">c{copy}_")));
            }
        }
    }
    std::fs::write(&set, text).unwrap();
    let db = dir.join("set").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(&set)
        .args(["-dbtype", "prot", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (query, set, db)
}

/// One search; returns standard output and standard error.
fn blastx(query: &Path, args: &[&str], threads: u32, env: &[(&str, &str)]) -> (Vec<u8>, String) {
    let output = losat()
        .arg("blastx")
        .arg("-query")
        .arg(query)
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .envs(env.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{args:?} {env:?}: {stderr}");
    (output.stdout, stderr)
}

/// Whether a search that wrote its reach files (`LOSAT_X_REACH`) into `reach` took `switch`'s own
/// path.
fn own_side_reached(reach: &Path, switch: &str) -> bool {
    std::fs::read_dir(reach).unwrap().any(|file| {
        std::fs::read_to_string(file.unwrap().path())
            .unwrap()
            .lines()
            .any(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                fields[0] == switch && fields[1] == "1"
            })
    })
}

/// The numbers of a `[TIMING] blastx_subject_flow` or `blastx_subject_waves` line: `name=value`.
fn timing(stderr: &str, line: &str, name: &str) -> u64 {
    let text = stderr
        .lines()
        .find(|l| l.starts_with(&format!("[TIMING] {line}: ")))
        .unwrap_or_else(|| panic!("no {line} line: {stderr}"));
    text.split_whitespace()
        .find_map(|word| word.strip_prefix(&format!("{name}=")))
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no {name} in {text}"))
}

#[test]
fn the_report_does_not_depend_on_the_subject_flow() {
    let dir = temp_dir();
    let (query, set, db) = inputs(&dir);
    let set = set.to_str().unwrap();
    for (subjects, outfmts) in [
        (["-db", db.as_str()], &["6", "0", "7"][..]),
        (["-subject", set], &["6"][..]),
    ] {
        for &outfmt in outfmts {
            let args = [subjects[0], subjects[1], "-outfmt", outfmt];
            let (reference, _) = blastx(&query, &args, 3, &[("LOSAT_X_BXSUBJFLOW", "0")]);
            assert!(
                reference.iter().filter(|&&b| b == b'\n').count() > 20,
                "{args:?}: too few hits to test"
            );
            for threads in [2, 3, 4] {
                for (k, cut) in [("1", "1"), ("2", "1"), ("7", "1"), ("16", "1"), ("16", "0")] {
                    let env = [("LOSAT_X_BXSUBJFLOW", k), ("LOSAT_X_BXTASKRES", cut)];
                    let (report, _) = blastx(&query, &args, threads, &env);
                    assert!(
                        report == reference,
                        "{args:?} threads {threads} {env:?}: report differs"
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn an_ungapped_search_and_the_shadow_give_the_same_report() {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-subject-flow-ungapped-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let (query, _, db) = inputs(&dir);
    for extra in [&[][..], &["-ungapped", "-comp_based_stats", "0"][..]] {
        let mut args = vec!["-db", db.as_str(), "-outfmt", "6"];
        args.extend_from_slice(extra);
        let (reference, _) = blastx(&query, &args, 4, &[("LOSAT_X_BXSUBJFLOW", "0")]);
        let (report, stderr) = blastx(
            &query,
            &args,
            4,
            &[
                ("LOSAT_X_BXSUBJFLOW", "2"),
                ("LOSAT_X_BXSUBJFLOWSHADOW", "1"),
            ],
        );
        assert!(report == reference, "{extra:?}: report differs");
        // Tasks of at most 2 subjects and twice the mean length (LOSAT_X_BXTASKRES).
        let tasks = 1_143;
        assert!(
            stderr.contains(&format!(
                "[X_SHADOW] BXSUBJFLOW tasks_compared={tasks} (subject results identical)"
            )),
            "{extra:?}: {stderr}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_flow_takes_the_subject_loop() {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-subject-flow-timing-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let (query, subjects, db) = inputs(&dir);
    let args = ["-db", db.as_str(), "-outfmt", "6"];
    let off = [("LOSAT_X_BXSUBJFLOW", "0"), ("LOSAT_TIMING", "1")];
    let (reference, stderr) = blastx(&query, &args, 2, &off);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 0);
    assert!(timing(&stderr, "blastx_subject_waves", "waves") > 0);
    // The default: the 1,711 subjects as 29 tasks of at most 64 subjects and 64 times the mean
    // length (LOSAT_X_BXTASKRES), and no wave; the grid of 64 with that switch off: 27 tasks.
    let reach = dir.join("reach-taskres");
    let reach_text = reach.to_str().unwrap().to_string();
    let (report, stderr) = blastx(
        &query,
        &args,
        2,
        &[
            ("LOSAT_TIMING", "1"),
            ("LOSAT_X_REACH", reach_text.as_str()),
        ],
    );
    assert!(report == reference);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 1);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "tasks"), 29);
    assert_eq!(timing(&stderr, "blastx_subject_waves", "waves"), 0);
    assert!(own_side_reached(&reach, "LOSAT_X_BXTASKRES"));
    let (report, stderr) = blastx(
        &query,
        &args,
        2,
        &[("LOSAT_TIMING", "1"), ("LOSAT_X_BXTASKRES", "0")],
    );
    assert!(report == reference);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "tasks"), 27);
    // A `-subject` file takes the same loop.
    let (subject_report, stderr) = blastx(
        &query,
        &["-subject", subjects.to_str().unwrap(), "-outfmt", "6"],
        2,
        &[("LOSAT_TIMING", "1")],
    );
    let (subject_reference, _) = blastx(
        &query,
        &["-subject", subjects.to_str().unwrap(), "-outfmt", "6"],
        2,
        &[("LOSAT_X_BXSUBJFLOW", "0")],
    );
    assert!(subject_report == subject_reference);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 1);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "tasks"), 29);
    assert_eq!(timing(&stderr, "blastx_subject_waves", "waves"), 0);
    // The helpers search: with 1,143 tasks of at most 2 subjects at 4 threads the loop takes tasks
    // that were searched ahead (only a helper, or the loop's thread while it waits for a helper's task,
    // searches ahead), so the helper tables are the ones under test above.
    let (report, stderr) = blastx(
        &query,
        &args,
        4,
        &[("LOSAT_X_BXSUBJFLOW", "2"), ("LOSAT_TIMING", "1")],
    );
    assert!(report == reference);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "tasks"), 1_143);
    assert!(timing(&stderr, "blastx_subject_flow", "ahead") > 0);
    assert!(timing(&stderr, "blastx_subject_flow", "taken") > 0);
    // One thread: no pool, so the subject loop of the waves runs (on this thread).
    let (report, stderr) = blastx(&query, &args, 1, &[("LOSAT_TIMING", "1")]);
    assert!(report == reference);
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 0);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_chunks_of_a_split_query_take_the_flow_one_after_another() {
    let dir = std::env::temp_dir().join(format!(
        "losat-blastx-subject-flow-chunks-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    // 30,000 letters: three query chunks (blastx/split.rs), each against 3,422 subjects (the 1,711
    // proteins twice), 57 tasks of at most 64 (LOSAT_X_BXTASKRES; the grid's 54 with it off).
    // LOSAT_X_BXCHUNKFLOW needs at least 2,048 subjects and 256 per thread: up to 8 threads the chunks run one after another, each with the flow; at 16
    // they run side by side.
    let (query, set, db) = inputs_of(&dir, 30_000, 2);
    let args = ["-db", db.as_str(), "-outfmt", "6"];
    let (reference, _) = blastx(&query, &args, 1, &[]);
    assert!(
        reference.iter().filter(|&&b| b == b'\n').count() > 20,
        "too few hits to test"
    );
    // (threads, environment, flow loops): in order at 2 and 8 threads, as one flow across the
    // three chunks (LOSAT_X_BXCHUNKCROSS), or one flow per chunk with that switch off; side by side
    // at 16 threads and with LOSAT_X_BXCHUNKFLOW off; with LOSAT_X_BXSUBJFLOW=0 the chunks run
    // side by side (2 threads) with the waves.
    for (threads, env, loops) in [
        (2, &[][..], 1),
        (8, &[][..], 1),
        (3, &[("LOSAT_X_BXBUILDPAR", "0")][..], 1),
        (2, &[("LOSAT_X_BXCHUNKCROSS", "0")][..], 3),
        (8, &[("LOSAT_X_BXCHUNKCROSS", "0")][..], 3),
        // The helpers' own diagonal tables; the batch inside the capped batch flow (K = 2 below
        // the thread count), whose helpers leave once every task is claimed.
        (2, &[("LOSAT_X_BXDIAGKEEP", "0")][..], 1),
        (4, &[("LOSAT_X_BXUNITCAP", "2")][..], 1),
        (16, &[][..], 0),
        (2, &[("LOSAT_X_BXTASKRES", "0")][..], 1),
        (2, &[("LOSAT_X_BXCHUNKFLOW", "0")][..], 0),
        (2, &[("LOSAT_X_BXSUBJFLOW", "0")][..], 0),
    ] {
        let mut env = env.to_vec();
        env.push(("LOSAT_TIMING", "1"));
        let (report, stderr) = blastx(&query, &args, threads, &env);
        assert!(
            report == reference,
            "threads {threads} {env:?}: report differs"
        );
        assert_eq!(
            timing(&stderr, "blastx_subject_flow", "loops"),
            loops,
            "threads {threads} {env:?}"
        );
        if loops > 0 {
            let grid = env.contains(&("LOSAT_X_BXTASKRES", "0"));
            let tasks = if grid { 162 } else { 171 };
            assert_eq!(timing(&stderr, "blastx_subject_flow", "tasks"), tasks);
        }
    }
    // The capped row runs inside the capped batch flow: the cap acted (`LOSAT_X_REACH`), which it
    // does not at 4 threads without `LOSAT_X_BXUNITCAP=2` (K = 8).
    for (cap, acted) in [(Some("2"), true), (None, false)] {
        let reach = dir.join(format!("reach-cap-{acted}"));
        let reach_text = reach.to_str().unwrap().to_string();
        let mut env = vec![("LOSAT_X_REACH", reach_text.as_str())];
        env.extend(cap.map(|k| ("LOSAT_X_BXUNITCAP", k)));
        let (report, _) = blastx(&query, &args, 4, &env);
        assert!(report == reference);
        assert_eq!(own_side_reached(&reach, "LOSAT_X_BXUNITCAP"), acted);
    }
    // The flow's shadow compares every task of the three chunks, across the chunks and per chunk.
    for env in [&[][..], &[("LOSAT_X_BXCHUNKCROSS", "0")][..]] {
        let mut env = env.to_vec();
        env.push(("LOSAT_X_BXSUBJFLOWSHADOW", "1"));
        let (report, stderr) = blastx(&query, &args, 3, &env);
        assert!(report == reference);
        assert!(
            stderr.contains("[X_SHADOW] BXSUBJFLOW tasks_compared=171 (subject results identical)"),
            "{stderr}"
        );
    }
    // Chunks without a valid context (all N): the loop gives them `chunk_start` only, in the flow
    // across chunks as in the serial loop. 10,000 letters, 25,000 N, 10,000 letters: four chunks
    // (blastx/split.rs), the third all N. The flow across chunks counts the tasks of all four (228,
    // the third chunk's skipped); the loop per chunk searches the other three (171).
    let letters: String = std::fs::read_to_string(&query)
        .unwrap()
        .lines()
        .skip(1)
        .collect();
    let gapped = dir.join("gapped.fna");
    std::fs::write(
        &gapped,
        format!(
            ">gapped\n{}{}{}\n",
            &letters[..10_000],
            "N".repeat(25_000),
            &letters[10_000..20_000]
        ),
    )
    .unwrap();
    let (reference, _) = blastx(&gapped, &args, 1, &[]);
    for (threads, env, loops, tasks) in [
        (2, &[][..], 1, 228),
        (8, &[][..], 1, 228),
        (2, &[("LOSAT_X_BXCHUNKCROSS", "0")][..], 3, 171),
    ] {
        let mut env = env.to_vec();
        env.push(("LOSAT_TIMING", "1"));
        let (report, stderr) = blastx(&gapped, &args, threads, &env);
        assert!(
            report == reference,
            "N stretch, threads {threads} {env:?}: report differs"
        );
        assert_eq!(
            timing(&stderr, "blastx_subject_flow", "loops"),
            loops,
            "N stretch, threads {threads} {env:?}"
        );
        assert_eq!(
            timing(&stderr, "blastx_subject_flow", "tasks"),
            tasks,
            "N stretch, threads {threads} {env:?}"
        );
    }
    // A `-subject` file takes the same rule.
    let args = ["-subject", set.to_str().unwrap(), "-outfmt", "6"];
    let (reference, _) = blastx(&query, &args, 1, &[]);
    for (threads, loops) in [(4, 1), (16, 0)] {
        let (report, stderr) = blastx(&query, &args, threads, &[("LOSAT_TIMING", "1")]);
        assert!(
            report == reference,
            "-subject threads {threads}: report differs"
        );
        assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), loops);
    }
    // Fewer than 2,048 subjects: side by side.
    std::fs::create_dir_all(dir.join("few")).unwrap();
    let (query, _, db) = inputs_of(&dir.join("few"), 30_000, 1);
    let (_, stderr) = blastx(
        &query,
        &["-db", db.as_str(), "-outfmt", "6"],
        2,
        &[("LOSAT_TIMING", "1")],
    );
    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 0);
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
// The frozen NCBI rows `db.pmask.bx.hard.ov_m1.*` (`tests/db_search_manifest/blastx.pmask.tsv`,
// two-hit and one-hit; review 4 of M4b-3b, L1). NCBI's one thread moves the diagonal offset by
// each subject's length and the window. In `ov_m1` the 64th subject, P (no masks), is one that
// NCBI's one-thread search starts on a new table below the entries that no sequence wrote, and the
// first subject of each later group of 64 is a copy of a subject with rows. With two threads or
// more the waves (`LOSAT_X_BXSUBJFLOW=0`) and the grid of 64 (`LOSAT_X_BXTASKRES=0`) give the worker
// that searched P one of those copies next, on P's table: the copy has NCBI's rows only where P's
// search leaves the offset past those entries (`DiagonalWalk::plain` and its end move), and the rows
// of one to four copies were lost before that. LOSAT_TIMING shows which loop ran.
#[test]
fn the_subject_after_a_sequence_on_a_new_table_has_ncbis_rows_in_waves_and_the_grid() {
    let manifest = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search_manifest/blastx.pmask.tsv"),
    )
    .unwrap();
    let mut lines = manifest.lines().filter(|line| !line.starts_with('#'));
    let columns: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let column = |name: &str| columns.iter().position(|c| *c == name).unwrap();
    let frozen = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/db_search");
    let mut rows = 0;
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        let row = cells[column("row_id")];
        if !row.starts_with("db.pmask.bx.hard.ov_m1.") {
            continue;
        }
        let expected_out = temp_support::expected_losat_bytes(
            std::fs::read(frozen.join(format!("{row}.out"))).unwrap(),
        );
        let expected_err = std::fs::read_to_string(frozen.join(format!("{row}.err")))
            .unwrap_or_default()
            .in_losat_words();
        let mut args: Vec<&str> = cells[column("extra_args")].split_whitespace().collect();
        assert!(
            args.iter().all(|arg| !arg.contains('=')),
            "{row}: a row without variables"
        );
        args.extend([
            "-db",
            cells[column("db")],
            "-outfmt",
            cells[column("outfmt")],
        ]);
        let query = Path::new(cells[column("query")]);
        for threads in [2, 4] {
            for switch in ["LOSAT_X_BXSUBJFLOW", "LOSAT_X_BXTASKRES"] {
                let what = format!("{row} threads {threads} {switch}=0");
                let env = [(switch, "0"), ("LOSAT_TIMING", "1")];
                let (code, stdout, stderr) = frozen_row_blastx(query, &args, threads, &env);
                assert_eq!(
                    code.map(|code| code.to_string()).as_deref(),
                    Some(cells[column("rc")]),
                    "{what}: {stderr}"
                );
                assert!(stdout == expected_out, "{what}: the bytes of NCBI differ");
                if switch == "LOSAT_X_BXSUBJFLOW" {
                    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 0, "{what}");
                    assert!(
                        timing(&stderr, "blastx_subject_waves", "waves") > 0,
                        "{what}: the waves did not run"
                    );
                } else {
                    // The grid of 64: the 600 subjects in 10 tasks.
                    assert_eq!(timing(&stderr, "blastx_subject_flow", "loops"), 1, "{what}");
                    assert_eq!(
                        timing(&stderr, "blastx_subject_flow", "tasks"),
                        10,
                        "{what}"
                    );
                }
                let stderr: String = stderr
                    .split_inclusive('\n')
                    .filter(|line| !line.starts_with("[TIMING]"))
                    .collect();
                assert_eq!(stderr, expected_err, "{what}: standard error differs");
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 2);
}

/// A frozen row's search from `LOSAT/` (the manifest's paths are relative to it) in the
/// environment of `db_search_fixtures.py check`: no `LOSAT_*` variable but the worker log (`losat`),
/// no BATCH_SIZE, NCBI* or BLAST* variable, HOME an empty directory; then `env`. Returns the exit
/// status, standard output and standard error.
fn frozen_row_blastx(
    query: &Path,
    args: &[&str],
    threads: u32,
    env: &[(&str, &str)],
) -> (Option<i32>, Vec<u8>, String) {
    let mut command = losat();
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text == "BATCH_SIZE" || text.starts_with("NCBI") || text.starts_with("BLAST") {
            command.env_remove(name);
        }
    }
    let _log = temp_support::worker_log_if_unset(&mut command, "blastx_subject_flow_frozen");
    let home = temp_support::unique_dir("losat-blastx-subject-flow-home");
    let output = command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .env("HOME", home.path())
        .arg("blastx")
        .arg("-query")
        .arg(query)
        .args(args)
        .args(["-num_threads", &threads.to_string()])
        .envs(env.iter().copied())
        .output()
        .expect("run LOSAT");
    (
        output.status.code(),
        output.stdout,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}
