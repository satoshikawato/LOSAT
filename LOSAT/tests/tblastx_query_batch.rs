//! TBLASTX's ungapped extensions past the end of a query frame (`tblastx::extension::QueryBatch`):
//! the searches that LOSAT refuses because such an extension scores higher in the next frame of the
//! query batch (PD-LOSAT-NCBI-DEFECTS, rule 3), at 1, 4, 8 and 16 threads. The searches where it
//! does not, which NCBI BLAST+ 2.17.0 reports as LOSAT does, are rows `x2c.tblastx.*` of
//! `tests/outfmt0_manifest.tsv`; one of them also runs here at 16 threads, where the unit path
//! splits the query's frames into several context groups.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:841 and :897
//! The walks count the letters of the whole query batch, so they read the NULLB between two frames
//! (BLAST_SCORE_MIN) and go on into the next frame while the X-drop allows it. Where they score
//! higher there, NCBI keeps the HSP in the seed's frame with query offsets outside it
//! (`s_AdjustInitialHSPOffsets`, blast_gapalign.c:2384-2392, whose `ASSERT(q_start >= 0)` is not in
//! release builds). The inputs (`tests/fasta/outfmt0/x2c_*`, written by line RG's `gen_x2c.py`):
//! - `d7k` (`-matrix PAM250 -threshold 999 -seg no -xdrop_ungap 20000 -evalue 1e-20`): frame 1 of
//!   the query ends with 7,000 letters of F, I and L, frame 2 starts with 40 letters of W, C and H
//!   (P); the subject's frame 1 holds the 7,000 letters mapped to Y, V and M, a G and P. P's seed
//!   extends left over the NULLB into frame 1 (+35,000 there). NCBI drops that HSP from its report
//!   (139 rows; LOSAT's dev build had 140, P/P among them), 0.19 s.
//! - `e7k` (the same options): the same from frame -1 back into frame 3. NCBI reports the HSP as
//!   `0.000 7040 7040 0 42245 21126 1 21120 0.0 11314` (`-outfmt "6 std ..."`): a query start of
//!   42,245 in a query of 21,122 bases, no identities.
//! - `g48` (`-seg no -threshold 999 -xdrop_ungap 20000 -evalue 1e-20`): frame 1 ends with 4,800
//!   letters of W C H P Y F N D G, frame 2 starts with another 4,800; the subject holds both, one G
//!   apart. The right walk of the first passes the NULLB with 35,000 and gains past it: NCBI reports
//!   `49.995 9601 4801 0 14404 43206 1 28803 0.0 17110`, a query end of 43,206 in a query of 28,803.
//! - `d7k` with `-window_size 0`: the one-hit word finder's extension, the same refusal.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The refusal's message on standard error.
const REFUSAL: &str = "Error: an ungapped extension that passes the end of a query frame into \
the frame next to it in the query batch and scores higher there";

fn manifest_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

/// `losat tblastx` from LOSAT/ on `tests/fasta/outfmt0/x2c_<case>_{q,s}.fna` with no `LOSAT_*`
/// variable but `LOSAT_BLASTX_WORKER_LOG` (a binary built with the `blastx-worker-probe` feature
/// needs that file name once a search has more than one thread; it selects nothing) and none of the
/// other variables that change a search (the list of `tests/blastn_window_masker.rs`).
fn tblastx(case: &str, args: &[&str], outfmt: &str, threads: u32) -> Output {
    tblastx_on(
        &manifest_path(&format!("tests/fasta/outfmt0/x2c_{case}_q.fna")),
        &manifest_path(&format!("tests/fasta/outfmt0/x2c_{case}_s.fna")),
        case,
        args,
        outfmt,
        threads,
    )
}

/// `tblastx` on the files `query` and `subject` (`tag` names its worker log).
fn tblastx_on(
    query: &Path,
    subject: &Path,
    tag: &str,
    args: &[&str],
    outfmt: &str,
    threads: u32,
) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE"
                    | "CHUNK_SIZE"
                    | "OVERLAP_CHUNK_SIZE"
                    | "RAYON_NUM_THREADS"
                    | "BLASTDB"
                    | "ARCHIVE_FORMAT"
                    | "WINDOW_MASKER_PATH"
                    | "BLAST_USAGE_REPORT"
                    | "BL2SEQ_LEGACY"
                    | "CTOOLKIT_COMPATIBLE"
                    | "OLD_FSC"
                    | "ADAPTIVE_CBS"
                    | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(name);
        }
    }
    let log = std::env::temp_dir().join(format!(
        "losat_tblastx_query_batch_{}_{tag}_{threads}.log",
        std::process::id()
    ));
    let output = command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .env("LOSAT_BLASTX_WORKER_LOG", &log)
        .arg("tblastx")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(args)
        .args(["-outfmt", outfmt, "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT");
    let _ = std::fs::remove_file(&log);
    output
}

/// The search of `case` with `args` is refused at 1, 4, 8 and 16 threads: exit status 1, no report,
/// the refusal's message.
fn refused(case: &str, args: &[&str]) {
    for threads in [1, 4, 8, 16] {
        let output = tblastx(case, args, "6", threads);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(
            output.status.code(),
            Some(1),
            "{case} threads {threads}: {stderr}"
        );
        assert!(
            output.stdout.is_empty(),
            "{case} threads {threads}: no report"
        );
        assert!(
            stderr.starts_with(REFUSAL)
                && stderr.ends_with("is not supported by LOSAT's TBLASTX\n"),
            "{case} threads {threads}: {stderr}"
        );
        assert_eq!(
            stderr.lines().count(),
            1,
            "{case} threads {threads}: {stderr}"
        );
    }
}

const PAM250_LEFT: [&str; 10] = [
    "-matrix",
    "PAM250",
    "-threshold",
    "999",
    "-seg",
    "no",
    "-xdrop_ungap",
    "20000",
    "-evalue",
    "1e-20",
];

/// `d7k`: a left walk from frame 2 that gains in frame 1.
#[test]
fn a_left_walk_that_gains_in_the_previous_frame_is_refused() {
    refused("d7k", &PAM250_LEFT);
}

/// `e7k`: a left walk from frame -1 that gains in frame 3.
#[test]
fn a_left_walk_from_a_minus_frame_that_gains_in_frame_3_is_refused() {
    refused("e7k", &PAM250_LEFT);
}

/// `g48`: a right walk from frame 1 that gains in frame 2.
#[test]
fn a_right_walk_that_gains_in_the_next_frame_is_refused() {
    refused(
        "g48",
        &[
            "-threshold",
            "999",
            "-seg",
            "no",
            "-xdrop_ungap",
            "20000",
            "-evalue",
            "1e-20",
        ],
    );
}

/// `d7k` with the one-hit word finder (`-window_size 0`).
#[test]
fn a_one_hit_extension_that_gains_in_the_previous_frame_is_refused() {
    let args: Vec<&str> = PAM250_LEFT
        .iter()
        .copied()
        .chain(["-window_size", "0"])
        .collect();
    refused("d7k", &args);
}

/// `f48` (row `x2c.tblastx.f48.6`): the right walk of R/R (4,800 letters, 49,995 bits) from frame
/// 1 goes on past the NULLB over P/P in frame 2 without scoring higher, and its stop covers every
/// P/P hit, so NCBI reports no P/P HSP. At 16 threads the unit path puts the two frames in
/// different context groups; the units run again in group order on one table
/// (`x_pair_par::x_units_in_order`) and give NCBI BLAST+ 2.17.0's rows, as at 1 thread.
#[test]
fn a_right_walk_past_the_nullb_covers_the_next_frames_hits_at_1_and_16_threads() {
    let expected = std::fs::read(manifest_path(
        "tests/fixtures/outfmt0/x2c.tblastx.f48.6.out",
    ))
    .expect("the frozen rows of x2c.tblastx.f48.6");
    let args = [
        "-seg",
        "no",
        "-threshold",
        "999",
        "-xdrop_ungap",
        "20000",
        "-evalue",
        "1e-40",
    ];
    for threads in [1, 16] {
        let output = tblastx(
            "f48",
            &args,
            "6 qframe sframe qstart qend sstart send evalue bitscore",
            threads,
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(0), "threads {threads}: {stderr}");
        assert!(stderr.is_empty(), "threads {threads}: {stderr}");
        assert!(
            output.stdout == expected,
            "threads {threads}: the rows differ from NCBI BLAST+ 2.17.0's"
        );
    }
}

/// Two query batches (`f48`'s query of 14,523 bases and `g48`'s of 28,803, each more than a
/// batch of 10,002 bases) against both subjects: the first batch's search reports (`f48`'s right
/// walk passes the NULLB without scoring higher), the second's is refused (`g48`). The Owner's
/// decision of 2026-11-16 (line RG, Q-X2c): the results of the query batches before are written,
/// then the error. The default tabular report is written batch by batch and the one with a field
/// list after the search (`report_stream`); both write the first batch's rows.
#[test]
fn the_results_of_the_query_batches_before_a_refusal_are_written() {
    let dir = std::env::temp_dir().join(format!("losat_x2c_two_batches_{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temporary directory");
    let join = |name: &str, parts: [&str; 2]| {
        let path = dir.join(name);
        let mut text = Vec::new();
        for part in parts {
            text.extend(
                std::fs::read(manifest_path(&format!("tests/fasta/outfmt0/{part}")))
                    .expect("x2c input"),
            );
        }
        std::fs::write(&path, text).expect("joined input");
        path
    };
    let query = join("q.fna", ["x2c_f48_q.fna", "x2c_g48_q.fna"]);
    let subject = join("s.fna", ["x2c_f48_s.fna", "x2c_g48_s.fna"]);
    let args = [
        "-seg",
        "no",
        "-threshold",
        "999",
        "-xdrop_ungap",
        "20000",
        "-evalue",
        "1e-20",
    ];
    let mut reports = Vec::new();
    for outfmt in ["6", "6 std"] {
        for threads in [1, 4] {
            let output = tblastx_on(&query, &subject, "two", &args, outfmt, threads);
            let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
            assert_eq!(
                output.status.code(),
                Some(1),
                "{outfmt} {threads}: {stderr}"
            );
            assert!(
                stderr.starts_with(REFUSAL)
                    && stderr.ends_with("is not supported by LOSAT's TBLASTX\n")
                    && stderr.lines().count() == 1,
                "{outfmt} {threads}: {stderr}"
            );
            let rows = String::from_utf8_lossy(&output.stdout).into_owned();
            assert!(
                !rows.is_empty() && rows.lines().all(|row| row.starts_with("x2cfq\t")),
                "{outfmt} {threads}: the first batch's rows only"
            );
            reports.push(rows);
        }
    }
    assert!(
        reports.iter().all(|rows| *rows == reports[0]),
        "the same rows in both formats and at 1 and 4 threads"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
