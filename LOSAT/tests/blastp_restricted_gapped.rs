//! BLASTP restricted gapped alignments (`s_RestrictedGappedAlign`) and their exact redo, with
//! and without `LOSAT_X_BPGAP`, against NCBI BLAST+ 2.17.0.
//!
//! `LOSAT_X_BPGAP` is in the strict set (`src/utils/xswitch.rs`): it is on by default, unless its
//! variable is `0`, or `LOSAT_X_DEFAULT=0` and its variable is not set. Every test states the
//! default it runs with: the report is NCBI's with no variable, with `LOSAT_X_BPGAP=0`, with
//! `LOSAT_X_DEFAULT=0` and with `LOSAT_X_DEFAULT=0 LOSAT_X_BPGAP=1` (the run the variable gave
//! before it joined the set). `LOSAT_X_BPGAP16` (the sixteen-bit rows, also in the set) is checked
//! the same way: with `LOSAT_X_BPGAP16=0`, and with `LOSAT_X_DEFAULT=0 LOSAT_X_BPGAP=1
//! LOSAT_X_BPGAP16=1`.
//!
//! NCBI reference (2.17.0): c++/src/algo/blast/core/blast_gapalign.c:3777-3786,3955-3957
//! ```c
//!           /* use approximate gapped alignment if the highest scoring ungapped
//!              alignment scores below the reduced cutoff */
//!           if (init_hsp->ungapped_data && init_hsp->ungapped_data->score <
//!               (Int4)(kRestrictedMult *
//!                      hit_params->cutoffs[contxt].cutoff_score)) {
//!
//!               restricted_align_array[query_index] = TRUE;
//! ...
//!            if (restricted_alignment &&
//!                gap_align->score < cutoff &&
//!                gap_align->score >= restricted_cutoff) {
//! ```
//! The blastp fixtures of `ci_fast_regressions.py` and `check_losat.py` never reach this path
//! (LOSAT_X_BPGAPSHADOW compared nothing on any of them): it needs a query whose best ungapped
//! alignment scores below 0.68 of the cutoff, and the cutoff grows with the search space. The
//! inputs here are generated:
//! - four random queries;
//! - 1,500 random subjects of 400 letters (0.6 Mb);
//! - at subjects 0, 375, 750 and 1,125, a copy of each query with 60 % of its letters replaced
//!   and an insertion or deletion about every eighth letter. These copies give gapped alignments
//!   above the cutoff whose ungapped pieces are short.
//!
//! On x86_64 with AVX2, LOSAT_X_BPGAPSHADOW printed `compared=380` on these inputs, and
//! `LOSAT_TIMING` counted five exact redos. The expected rows are NCBI BLAST+ 2.17.0's
//! (`blastp -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score"
//! -evalue 0.001 -query query.faa -subject subject.faa`), on the files this test writes.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const DEFAULT_OFF: (&str, &str) = ("LOSAT_X_DEFAULT", "0");

const OUTFMT: &str = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score";

const NCBI_ROWS: &str = concat!(
    "q0\th0\t31.389\t360\t12\t359\t13\t337\t4.40e-14\t68.2\t165\n",
    "q1\th375\t38.330\t467\t5\t450\t9\t454\t9.31e-40\t146\t368\n",
    "q2\th750\t33.947\t489\t5\t464\t4\t455\t1.08e-26\t107\t268\n",
    "q3\th1125\t36.895\t496\t5\t473\t6\t464\t3.52e-35\t132\t332\n",
);

/// xorshift64*: the inputs are the same on every platform.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    /// A letter with the background frequencies (per 10,000) of Robinson and Robinson.
    fn letter(&mut self) -> u8 {
        const LETTERS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
        const WEIGHTS: [u64; 20] = [
            825, 553, 406, 545, 137, 393, 675, 707, 227, 596, 966, 584, 242, 386, 470, 656, 534,
            108, 292, 687,
        ];
        let total: u64 = WEIGHTS.iter().sum();
        let mut pick = self.below(total);
        for (letter, weight) in LETTERS.iter().zip(WEIGHTS) {
            if pick < weight {
                return *letter;
            }
            pick -= weight;
        }
        unreachable!()
    }

    fn protein(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.letter()).collect()
    }

    /// `src` with 60 % of its letters replaced and an insertion or deletion of one or two letters
    /// about every eighth letter.
    fn distant_copy(&mut self, src: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(src.len() + 64);
        let mut i = 0;
        while i < src.len() {
            if self.below(8) == 0 {
                if self.below(2) == 0 {
                    i += 1 + self.below(2) as usize;
                    continue;
                }
                for _ in 0..1 + self.below(2) {
                    out.push(self.letter());
                }
            }
            if self.below(10) < 6 {
                out.push(self.letter());
            } else {
                out.push(src[i]);
            }
            i += 1;
        }
        out
    }
}

fn fasta(records: &[(String, Vec<u8>)]) -> String {
    let mut text = String::new();
    for (name, residues) in records {
        text.push('>');
        text.push_str(name);
        text.push('\n');
        for line in residues.chunks(60) {
            text.push_str(std::str::from_utf8(line).expect("letters"));
            text.push('\n');
        }
    }
    text
}

/// The query and subject files, written once.
fn inputs() -> &'static (PathBuf, PathBuf) {
    static FILES: OnceLock<(PathBuf, PathBuf)> = OnceLock::new();
    FILES.get_or_init(|| {
        let mut rng = Rng(0x4444);
        let queries: Vec<(String, Vec<u8>)> = (0..4)
            .map(|k| {
                let len = 300 + rng.below(200) as usize;
                (format!("q{k}"), rng.protein(len))
            })
            .collect();
        let subjects: Vec<(String, Vec<u8>)> = (0..1500)
            .map(|k| {
                if k % 375 == 0 {
                    (format!("h{k}"), rng.distant_copy(&queries[k / 375].1))
                } else {
                    (format!("r{k}"), rng.protein(400))
                }
            })
            .collect();
        // One directory per process: nextest runs each test in a process of its own, and a
        // process that wrote these files again while another one's search read them cut its input.
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
            .join("blastp_restricted_gapped")
            .join(std::process::id().to_string());
        std::fs::create_dir_all(&dir).expect("work directory");
        let files = (dir.join("query.faa"), dir.join("subject.faa"));
        std::fs::write(&files.0, fasta(&queries)).expect("query file");
        std::fs::write(&files.1, fasta(&subjects)).expect("subject file");
        files
    })
}

/// One BLASTP run with `variables` and no other `LOSAT_*` variable except
/// `LOSAT_BLASTX_WORKER_LOG` (a build with the worker probe writes to it as soon as a pool
/// starts, so the caller's value is kept, or a file of the tests is named). `LOSAT_X_DEFAULT` of
/// the test's own environment is dropped too: every test states the default it runs with.
/// Returns the report and standard error.
fn blastp(threads: usize, variables: &[(&str, &str)]) -> (String, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let (query, subject) = inputs();
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(WORKER_LOG, query.with_file_name("worker.log"));
    }
    command.envs(variables.iter().copied());
    let output = command
        .arg("blastp")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(["-outfmt", OUTFMT, "-evalue", "0.001"])
        .args(["-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{stderr}");
    (
        String::from_utf8(output.stdout).expect("UTF-8 report"),
        stderr,
    )
}

/// The CPU runs the vector restricted alignment (`x_restricted_gap.rs`: AVX2, BMI1, LZCNT).
fn vector_kernel_runs() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("avx2")
            && is_x86_feature_detected!("bmi1")
            && is_x86_feature_detected!("lzcnt")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

#[test]
fn restricted_gapped_alignments_match_ncbi_with_and_without_the_switch() {
    for threads in [1, 4] {
        for variables in [
            &[][..],
            &[("LOSAT_X_BPGAP", "0")][..],
            &[DEFAULT_OFF][..],
            &[DEFAULT_OFF, ("LOSAT_X_BPGAP", "1")][..],
            &[("LOSAT_X_BPGAP16", "0")][..],
            &[
                DEFAULT_OFF,
                ("LOSAT_X_BPGAP", "1"),
                ("LOSAT_X_BPGAP16", "1"),
            ][..],
        ] {
            let (report, _) = blastp(threads, variables);
            assert_eq!(report, NCBI_ROWS, "{threads} threads, {variables:?}");
        }
    }
}

#[test]
fn the_shadow_mode_compares_the_restricted_alignments() {
    let (report, stderr) = blastp(1, &[("LOSAT_X_BPGAPSHADOW", "1")]);
    assert_eq!(report, NCBI_ROWS);
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_SHADOW] BPGAP "))
        .unwrap_or_else(|| panic!("no shadow line in {stderr:?}"));
    let count = |key: &str| -> u64 {
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(key))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {key} in {line:?}"))
    };
    let (compared, not_taken) = (count("compared="), count("not_taken="));
    if vector_kernel_runs() {
        assert!(compared > 0 && not_taken == 0, "{line}");
        assert!(
            line.ends_with("(every restricted gapped alignment identical)"),
            "{line}"
        );
    } else {
        assert!(compared == 0 && not_taken > 0, "{line}");
    }
}
