//! BLASTP restricted gapped alignments (`s_RestrictedGappedAlign`) and their exact redo, with
//! and without `LOSAT_X_BPGAP`, against NCBI BLAST+ 2.17.0, and psiblast's first round, which
//! aligns exactly.
//!
//! `LOSAT_X_BPGAP` is in the strict set (`src/utils/xswitch.rs`): it is on by default, unless its
//! variable is `0`. Every test states the
//! default it runs with: the report is NCBI's with no variable, with `LOSAT_X_BPGAP=0`, with
//! `LOSAT_X_BPGAP16=0` (the sixteen-bit rows, also in the set) and with `LOSAT_X_BPGAPBAND=0`
//! (the rows' arrays sized for the subject letters instead of the band, also in the set).
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
//! (`[X_STATS] BPGAP` counted no restricted alignment on any of them): it needs a query whose best
//! ungapped alignment scores below 0.68 of the cutoff, and the cutoff grows with the search space.
//! The inputs here are generated:
//! - four random queries;
//! - 1,500 random subjects of 400 letters (0.6 Mb);
//! - at subjects 0, 375, 750 and 1,125, a copy of each query with 60 % of its letters replaced
//!   and an insertion or deletion about every eighth letter. These copies give gapped alignments
//!   above the cutoff whose ungapped pieces are short.
//!
//! On x86_64 with AVX2, `[X_STATS] BPGAP` counted `fast=380` on these inputs, and
//! `LOSAT_TIMING` counted five exact redos. The expected rows are NCBI BLAST+ 2.17.0's
//! (`blastp -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score"
//! -evalue 0.001 -query query.faa -subject subject.faa`), on the files this test writes.
//!
//! psiblast's first round of a FASTA query runs the same engine but aligns exactly: NCBI sets
//! restricted alignment for blastp's program only, and psiblast's is `eBlastTypePsiBlast`.
//!
//! NCBI reference (2.17.0): c++/src/algo/blast/core/blast_parameters.c:852-855
//! ```c
//!    if (program_number == eBlastTypeBlastp && gapped_calculation &&
//!       options->expect_value <= RESTRICTED_ALIGNMENT_WORST_EVALUE) {
//!       params->restricted_align = TRUE;
//!    }
//! ```
//! NCBI reference (2.17.0): c++/src/algo/blast/api/psiblast_options.cpp:50 and
//! c++/src/algo/blast/api/blast_aux.cpp:739-741
//! ```c++
//!     m_Opts->SetProgram(ePSIBlast);
//! ...
//!     case ePSIBlast:
//!     case eDeltaBlast:
//!         return eBlastTypePsiBlast;
//! ```
//! Its inputs are generated in the same way with another seed (series RG Q3d's set a2): five
//! copies of each query, with 55 % to 64 % of the letters replaced and an insertion or deletion
//! about every sixth letter, among 1,500 random subjects. With restricted alignment the hit of
//! `q3` on its 64 % copy was missing from the report. The expected rows are NCBI BLAST+ 2.17.0's
//! (`psiblast -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send
//! evalue bitscore score" -evalue 1e-3 -comp_based_stats 0|2 -query query.faa -subject
//! subject.faa -num_threads 1`).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const OUTFMT: &str = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore score";

const NCBI_ROWS: &str = concat!(
    "q0\th0\t31.389\t360\t12\t359\t13\t337\t4.40e-14\t68.2\t165\n",
    "q1\th375\t38.330\t467\t5\t450\t9\t454\t9.31e-40\t146\t368\n",
    "q2\th750\t33.947\t489\t5\t464\t4\t455\t1.08e-26\t107\t268\n",
    "q3\th1125\t36.895\t496\t5\t473\t6\t464\t3.52e-35\t132\t332\n",
);

const PSIBLAST_OUTFMT: &str =
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score";

const NCBI_PSIBLAST_ROWS_CBS0: &str = concat!(
    "q0\th0_0_55_6\t42.767\t477\t209\t29\t4\t450\t2\t444\t8.23e-52\t177\t448\n",
    "q0\th900_0_62_6\t35.345\t464\t257\t20\t5\t445\t2\t445\t3.01e-39\t143\t360\n",
    "q0\th1200_0_64_6\t34.416\t462\t244\t24\t15\t446\t13\t445\t1.79e-32\t124\t311\n",
    "q0\th300_0_58_6\t37.071\t437\t236\t22\t31\t449\t24\t439\t2.16e-32\t124\t310\n",
    "q0\th600_0_60_6\t28.270\t474\t293\t18\t3\t450\t3\t455\t7.54e-20\t87.4\t215\n",
    "q1\th675_1_60_6\t37.910\t335\t175\t16\t6\t320\t1\t322\t4.38e-40\t139\t350\n",
    "q1\th75_1_55_6\t42.683\t328\t151\t21\t3\t312\t2\t310\t6.53e-27\t103\t258\n",
    "q1\th375_1_58_6\t35.347\t331\t172\t18\t19\t320\t22\t339\t3.82e-20\t85.5\t210\n",
    "q1\th1275_1_64_6\t32.836\t335\t188\t16\t1\t315\t1\t318\t3.74e-19\t82.4\t202\n",
    "q2\th1050_2_62_6\t34.940\t498\t275\t20\t17\t480\t17\t499\t4.93e-37\t138\t348\n",
    "q2\th1350_2_64_6\t32.809\t509\t273\t24\t7\t484\t6\t476\t7.10e-30\t117\t294\n",
    "q2\th150_2_55_6\t33.733\t501\t269\t29\t15\t484\t15\t483\t6.36e-29\t115\t287\n",
    "q2\th450_2_58_6\t33.008\t512\t266\t28\t6\t474\t4\t481\t5.58e-28\t112\t280\n",
    "q2\th750_2_60_6\t35.841\t452\t232\t29\t49\t467\t53\t479\t5.92e-24\t100\t249\n",
    "q3\th825_3_60_6\t37.116\t423\t213\t24\t1\t389\t1\t404\t1.04e-25\t103\t257\n",
    "q3\th225_3_55_6\t35.904\t376\t199\t21\t44\t391\t35\t396\t1.12e-19\t85.9\t211\n",
    "q3\th525_3_58_6\t33.590\t390\t218\t21\t16\t385\t16\t384\t7.48e-16\t74.3\t181\n",
    "q3\th1125_3_62_6\t28.763\t372\t233\t16\t11\t368\t10\t363\t4.33e-10\t56.6\t135\n",
    "q3\th1425_3_64_6\t28.732\t355\t201\t21\t37\t372\t39\t360\t1.35e-09\t55.1\t131\n",
);

const NCBI_PSIBLAST_ROWS_CBS2: &str = concat!(
    "q0\th0_0_55_6\t42.558\t477\t210\t28\t4\t450\t2\t444\t2.17e-54\t184\t466\n",
    "q0\th900_0_62_6\t35.345\t464\t257\t20\t5\t445\t2\t445\t8.28e-40\t144\t364\n",
    "q0\th1200_0_64_6\t36.364\t429\t211\t27\t15\t413\t13\t409\t8.11e-34\t128\t321\n",
    "q0\th300_0_58_6\t37.443\t438\t233\t23\t31\t449\t24\t439\t1.11e-33\t127\t319\n",
    "q0\th600_0_60_6\t29.019\t479\t283\t20\t3\t450\t3\t455\t8.41e-20\t87.4\t215\n",
    "q1\th675_1_60_6\t37.798\t336\t176\t16\t6\t321\t1\t323\t2.10e-41\t142\t359\n",
    "q1\th75_1_55_6\t44.109\t331\t142\t25\t3\t312\t2\t310\t1.61e-29\t110\t276\n",
    "q1\th375_1_58_6\t36.080\t352\t180\t22\t1\t320\t1\t339\t7.35e-22\t90.5\t223\n",
    "q1\th1275_1_64_6\t33.036\t336\t186\t17\t1\t315\t1\t318\t4.38e-20\t85.1\t209\n",
    "q2\th1050_2_62_6\t35.686\t510\t277\t21\t6\t480\t6\t499\t9.48e-39\t143\t361\n",
    "q2\th150_2_55_6\t34.990\t503\t260\t32\t15\t484\t15\t483\t1.43e-30\t119\t299\n",
    "q2\th1350_2_64_6\t33.659\t511\t266\t27\t7\t484\t6\t476\t2.98e-30\t119\t297\n",
    "q2\th450_2_58_6\t33.528\t513\t262\t30\t6\t474\t4\t481\t1.21e-28\t114\t285\n",
    "q2\th750_2_60_6\t35.762\t453\t233\t29\t49\t468\t53\t480\t1.59e-24\t102\t253\n",
    "q3\th825_3_60_6\t37.825\t423\t210\t25\t1\t389\t1\t404\t8.00e-28\t109\t273\n",
    "q3\th225_3_55_6\t36.000\t375\t191\t22\t52\t394\t44\t401\t6.03e-21\t89.7\t221\n",
    "q3\th525_3_58_6\t33.418\t392\t216\t22\t16\t385\t16\t384\t3.04e-17\t78.6\t192\n",
    "q3\th1425_3_64_6\t29.213\t356\t198\t23\t37\t372\t39\t360\t1.66e-11\t60.8\t146\n",
    "q3\th1125_3_62_6\t29.708\t377\t225\t17\t10\t368\t9\t363\t2.77e-09\t54.3\t129\n",
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
        self.distant_copy_with(src, (6, 10), 8)
    }

    /// `src` with a fraction `replaced.0 / replaced.1` of its letters replaced and an insertion
    /// or deletion of one or two letters about every `indel_every`th letter.
    fn distant_copy_with(&mut self, src: &[u8], replaced: (u64, u64), indel_every: u64) -> Vec<u8> {
        let mut out = Vec::with_capacity(src.len() + 64);
        let mut i = 0;
        while i < src.len() {
            if self.below(indel_every) == 0 {
                if self.below(2) == 0 {
                    i += 1 + self.below(2) as usize;
                    continue;
                }
                for _ in 0..1 + self.below(2) {
                    out.push(self.letter());
                }
            }
            if self.below(replaced.1) < replaced.0 {
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
        write_inputs("blastp", &queries, &subjects)
    })
}

/// The query and subject files of the psiblast test (series RG Q3d's set a2), written once:
/// four random queries, and 1,500 subjects of 400 random letters except at every 75th subject
/// from 0 to 1,425, a copy of query `j % 4` (`j` the copy's rank) with 55, 58, 60, 62 or 64 % of
/// its letters replaced (the five copies of a query in that order) and an insertion or deletion
/// about every sixth letter.
fn psiblast_inputs() -> &'static (PathBuf, PathBuf) {
    static FILES: OnceLock<(PathBuf, PathBuf)> = OnceLock::new();
    FILES.get_or_init(|| {
        const REPLACED: [u64; 5] = [55, 58, 60, 62, 64];
        let mut rng = Rng(0x5102);
        let queries: Vec<(String, Vec<u8>)> = (0..4)
            .map(|k| {
                let len = 300 + rng.below(200) as usize;
                (format!("q{k}"), rng.protein(len))
            })
            .collect();
        let subjects: Vec<(String, Vec<u8>)> = (0..1500)
            .map(|k| {
                let copy = k / 75;
                if k % 75 == 0 && copy < 4 * REPLACED.len() {
                    let (query, replaced) = (copy % 4, REPLACED[copy / 4]);
                    (
                        format!("h{k}_{query}_{replaced}_6"),
                        rng.distant_copy_with(&queries[query].1, (replaced, 100), 6),
                    )
                } else {
                    (format!("r{k}"), rng.protein(400))
                }
            })
            .collect();
        write_inputs("psiblast", &queries, &subjects)
    })
}

/// Writes the query and subject files of `name`.
fn write_inputs(
    name: &str,
    queries: &[(String, Vec<u8>)],
    subjects: &[(String, Vec<u8>)],
) -> (PathBuf, PathBuf) {
    // One directory per process: nextest runs each test in a process of its own, and a
    // process that wrote these files again while another one's search read them cut its input.
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("blastp_restricted_gapped")
        .join(std::process::id().to_string())
        .join(name);
    std::fs::create_dir_all(&dir).expect("work directory");
    let files = (dir.join("query.faa"), dir.join("subject.faa"));
    std::fs::write(&files.0, fasta(queries)).expect("query file");
    std::fs::write(&files.1, fasta(subjects)).expect("subject file");
    files
}

/// One BLASTP run of the BLASTP inputs (`run`).
fn blastp(threads: usize, variables: &[(&str, &str)]) -> (String, String) {
    let args = ["-outfmt", OUTFMT, "-evalue", "0.001"];
    run("blastp", inputs(), &args, threads, variables)
}

/// One psiblast run (one round) of the psiblast inputs with `-comp_based_stats cbs` (`run`).
fn psiblast(cbs: &str, threads: usize, variables: &[(&str, &str)]) -> (String, String) {
    let args = [
        "-outfmt",
        PSIBLAST_OUTFMT,
        "-evalue",
        "1e-3",
        "-comp_based_stats",
        cbs,
    ];
    run("psiblast", psiblast_inputs(), &args, threads, variables)
}

/// One run of `program` on `files` (query, subject) with `args`, with `variables` and no other
/// `LOSAT_*` variable except `LOSAT_BLASTX_WORKER_LOG` (a build with the worker probe writes to it
/// as soon as a pool starts, so the caller's value is kept, or a file of the tests is named).
/// Returns the report and standard error.
fn run(
    program: &str,
    files: &(PathBuf, PathBuf),
    args: &[&str],
    threads: usize,
    variables: &[(&str, &str)],
) -> (String, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let (query, subject) = files;
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
        .arg(program)
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(args)
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
    for threads in [1, 4, 8, 16] {
        for variables in [
            &[][..],
            &[("LOSAT_X_BPGAP", "0")][..],
            &[("LOSAT_X_BPGAP16", "0")][..],
            &[("LOSAT_X_BPGAPBAND", "0")][..],
        ] {
            let (report, _) = blastp(threads, variables);
            assert_eq!(report, NCBI_ROWS, "{threads} threads, {variables:?}");
        }
    }
}

/// `[X_STATS] BPGAP fast= not_taken=`: the restricted alignments the vector kernel computed
/// and those it passed on to the reference function. With the kernel on this CPU every
/// alignment of the fixture is its own and the report is NCBI's.
#[test]
fn the_vector_kernel_computes_the_restricted_alignments_of_the_fixture() {
    let (report, stderr) = blastp(1, &[("LOSAT_X_STATS", "1")]);
    assert_eq!(report, NCBI_ROWS);
    let (fast, not_taken) = bpgap_counts(&stderr);
    let line = format!("fast={fast} not_taken={not_taken}");
    if vector_kernel_runs() {
        assert!(fast > 0 && not_taken == 0, "{line}");
    } else {
        assert!(fast == 0 && not_taken > 0, "{line}");
    }
}

/// The counts of `[X_STATS] BPGAP fast= not_taken=` in `stderr`.
fn bpgap_counts(stderr: &str) -> (u64, u64) {
    let line = stderr
        .lines()
        .find(|line| line.starts_with("[X_STATS] BPGAP "))
        .unwrap_or_else(|| panic!("no BPGAP line in {stderr:?}"));
    let count = |key: &str| -> u64 {
        line.split_whitespace()
            .find_map(|field| field.strip_prefix(key))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {key} in {line:?}"))
    };
    (count("fast="), count("not_taken="))
}

/// psiblast's first round of a FASTA query aligns exactly, as NCBI's (whose program is
/// `eBlastTypePsiBlast`, not `eBlastTypeBlastp`), without and with composition-based statistics.
#[test]
fn psiblast_first_round_matches_ncbi_with_exact_alignments() {
    for (cbs, rows) in [
        ("0", NCBI_PSIBLAST_ROWS_CBS0),
        ("2", NCBI_PSIBLAST_ROWS_CBS2),
    ] {
        for threads in [1, 4, 8] {
            let (report, _) = psiblast(cbs, threads, &[]);
            assert_eq!(report, rows, "-comp_based_stats {cbs}, {threads} threads");
        }
    }
}

/// `[X_STATS] BPGAP`: psiblast computes no restricted alignment on inputs where BLASTP does.
#[test]
fn psiblast_computes_no_restricted_alignment() {
    let (_, stderr) = psiblast("0", 1, &[("LOSAT_X_STATS", "1")]);
    assert_eq!(bpgap_counts(&stderr), (0, 0), "psiblast: {stderr}");
    let args = [
        "-outfmt",
        PSIBLAST_OUTFMT,
        "-evalue",
        "1e-3",
        "-comp_based_stats",
        "0",
    ];
    let (_, stderr) = run(
        "blastp",
        psiblast_inputs(),
        &args,
        1,
        &[("LOSAT_X_STATS", "1")],
    );
    let (fast, not_taken) = bpgap_counts(&stderr);
    assert!(fast + not_taken > 0, "blastp: {stderr}");
}
