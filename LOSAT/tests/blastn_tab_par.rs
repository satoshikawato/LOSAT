//! BLASTN's outfmt 6 / 7 rows formatted on the search pool (`LOSAT_X_BNTABPAR`,
//! `src/algorithm/blastn/x_tab_par.rs`) on whole searches.
//!
//! The switch is in the strict set: on unless it is `0`. It chooses between two LOSAT code paths
//! with the same result: the rows of a report written query after query on the thread of the
//! batch loop, or the rows cut into pieces of one query, each formatted into a buffer of its own
//! on the pool, and the buffers written in order. So the property under test is that the report
//! does not depend on it: for megablast and blastn with 60 queries (pieces of LC738874) against
//! LC738874, with outfmt 6 and 7 (whose header lines and last line come from the same writer) at
//! 1 and 4 threads, the default and `LOSAT_X_BNTABPAR=0` give the same bytes. `LOSAT_X_REACH`
//! shows that the default took the pool's path at 4 threads and the serial one at 1 thread (no
//! pool). That the rows equal NCBI's is the business of the blastn fixtures; the unit tests in
//! `x_tab_par.rs` compare the two writers on random hit lists, with warnings and windows of every
//! size.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

mod reach_support;

const OFF: (&str, &str) = ("LOSAT_X_BNTABPAR", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// A query file removed when dropped: 60 pieces of LC738874 of 300 to 2,700 letters, every fifth
/// reverse-complemented.
struct Queries(PathBuf);

impl Queries {
    fn new() -> Queries {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let letters: Vec<u8> = std::fs::read(fasta("LC738874.fasta"))
            .expect("read FASTA")
            .split(|&byte| byte == b'\n')
            .skip(1)
            .flat_map(|line| line.iter().copied().filter(|b| !b.is_ascii_whitespace()))
            .collect();
        let mut text = String::new();
        for i in 0..60usize {
            let start = (i * 4_513) % (letters.len() - 3_000);
            let len = 300 + (i * 977) % 2_400;
            let mut piece = letters[start..start + len].to_vec();
            if i % 5 == 4 {
                piece.reverse();
                for base in piece.iter_mut() {
                    *base = match base.to_ascii_uppercase() {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'T' => b'A',
                        other => other,
                    };
                }
            }
            text.push_str(&format!(">piece_{i} LC738874 {start}+{len}\n"));
            for line in piece.chunks(70) {
                text.push_str(std::str::from_utf8(line).unwrap());
                text.push('\n');
            }
        }
        let path = std::env::temp_dir().join(format!(
            "losat_bntabpar_{}_{}.fa",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, text).expect("write queries");
        Queries(path)
    }
}

impl Drop for Queries {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// The report of `blastn -task <task> -outfmt <outfmt>` of the queries against LC738874, with no
/// `LOSAT_*` variable of the caller but `switches` (`LOSAT_BLASTX_WORKER_LOG` stays: a binary built
/// with the `blastx-worker-probe` feature needs it at two threads or more).
fn blastn(
    queries: &Queries,
    task: &str,
    outfmt: &str,
    threads: u32,
    switches: &[(&str, &str)],
) -> Vec<u8> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG")
            || matches!(&*text, "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE")
        {
            command.env_remove(name);
        }
    }
    let output = command
        .arg("blastn")
        .args(["-task", task])
        .arg("-query")
        .arg(&queries.0)
        .arg("-subject")
        .arg(fasta("LC738874.fasta"))
        .args(["-outfmt", outfmt])
        .args(["-num_threads", &threads.to_string()])
        .envs(switches.iter().copied())
        .output()
        .expect("run LOSAT");
    assert!(
        output.status.success(),
        "{task} outfmt {outfmt} threads {threads} {switches:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

#[test]
fn tabular_rows_do_not_depend_on_the_switch() {
    let queries = Queries::new();
    for task in ["megablast", "blastn"] {
        for outfmt in ["6", "7"] {
            for threads in [1u32, 4] {
                let reach =
                    reach_support::Reach::new(&format!("bntabpar-{task}-{outfmt}-{threads}"));
                let default = blastn(&queries, task, outfmt, threads, &[reach.var()]);
                let off = blastn(&queries, task, outfmt, threads, &[OFF]);
                let rows = default
                    .split(|&b| b == b'\n')
                    .filter(|l| !l.is_empty() && l[0] != b'#')
                    .count();
                assert!(
                    rows > 100,
                    "{task} outfmt {outfmt} threads {threads}: {rows} rows"
                );
                assert!(
                    default == off,
                    "{task} outfmt {outfmt} threads {threads}: the default's report differs from LOSAT_X_BNTABPAR=0's"
                );
                let (new, reference) = reach.sides("LOSAT_X_BNTABPAR");
                assert_eq!(
                    (new, reference),
                    (threads > 1, threads == 1),
                    "{task} outfmt {outfmt} threads {threads}: sides of LOSAT_X_BNTABPAR"
                );
            }
        }
    }
}
