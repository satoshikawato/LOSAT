//! TBLASTN and TBLASTX `-db_soft_mask` on long database sequences with masks stored out of order
//! (series NP, M4b-T3):
//!
//! 1. Short initial alignments (at most 11 residues) starting before a 5,000,000-residue chunk
//!    whose middle lies inside the chunk are searched with NCBI BLAST+ 2.17.0's rows (`q1a1`,
//!    `q1a2`), while alignments longer than 11 residues (`q1d1`) or with their middle before the
//!    chunk (`q1c1`), where NCBI crashes with SIGSEGV, are refused.
//! 2. When NCBI's one-thread diagonal table clears (at INT4_MAX / 4) inside such a sequence (`C0`
//!    vs one-base-shorter `C2`, ~284 Mbp each, `#[ignore]`), LOSAT computes NCBI's one-thread
//!    table start from database order and reproduces NCBI's `-num_threads 1` rows (or refusal on
//!    `C2`) at every thread count.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3393-3413,3936-3953,
//! 4300-4333 and c++/src/algo/blast/core/blast_extend.c:87-107,162-173
//!
//! Setting `TBLASTN_DB_SOFT_LEAD_START_DIR=<dir>` keeps the generated FASTA, `.asnt`, and database
//! files in `<dir>` instead of removing them.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// Length of the chunked long sequence (frame +1 has 5,100,000 residues, split into two chunks).
const LONG_LEN: usize = 15_300_000;
/// Start residue of frame +1's second chunk.
const CHUNK2_START: usize = 4_999_900;
/// Control core residue inside chunk 2.
const CONTROL_AT: usize = CHUNK2_START + 20_000;

/// Probe query protein (118 residues: seven unique cores separated by `GWG` flanks).
const QUERY_AA: &str = "MGWGNDQYCGWQEWHGWGGWGYRDEECYWNWWGWGGWGWPHCNCWCGWGGWGNWCYHWMQGGPGWGGWGHNDMDQRQHEDGWGGWGEYDWYGYPPHQGWGGWGDPHWYDNYNMNRGWG";

/// Case cores used by the always-run and ignored tests.
const CORE_A1: &str = "NDQYCGWQEWH";
const CORE_A2: &str = "YRDEECYWNWW";
const CORE_E1: &str = "NWCYHWMQGGP";
const CORE_B1: &str = "HNDMDQRQHED";
const CORE_C1: &str = "EYDWYGYPPHQ";
const CORE_D1: &str = "DPHWYDNYNMNR";

/// Out-of-order dust masks `[from, to]` inclusive in bases, in stored order.
const MASKS: [(usize, usize); 3] = [
    (14_999_400, 14_999_849),
    (14_999_880, 15_179_999),
    (3_000_000, 3_299_999),
];

/// Extended tabular columns used by the rich-column frozen fixtures.
const RICH_OUTFMT: &str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore score nident positive btop";

/// Refusal phrase shared by all pre-chunk NCBI-crash refusals.
const REFUSAL_SUFFIX: &str = "an initial alignment of the chunk from residue 4999901 starts before the chunk (-db_soft_mask with masks stored out of order on a frame longer than 5,000,000 residues) and is longer than 11 residues or has its middle before the chunk";

/// Encodes a protein string using one fixed codon per amino acid.
fn protein_dna(protein: &str) -> String {
    let mut dna = String::with_capacity(protein.len() * 3);
    for ch in protein.bytes() {
        let codon = match ch {
            b'A' => "GCT",
            b'R' => "CGT",
            b'N' => "AAC",
            b'D' => "GAC",
            b'C' => "TGC",
            b'Q' => "CAG",
            b'E' => "GAG",
            b'G' => "GGT",
            b'H' => "CAC",
            b'I' => "ATC",
            b'L' => "CTG",
            b'K' => "AAG",
            b'M' => "ATG",
            b'F' => "TTC",
            b'P' => "CCG",
            b'S' => "TCT",
            b'T' => "ACC",
            b'W' => "TGG",
            b'Y' => "TAC",
            b'V' => "GTG",
            other => panic!("unexpected amino acid byte {other}"),
        };
        dna.push_str(codon);
    }
    dna
}

/// Writes FASTA records with 80-base lines, streaming poly-A fillers without huge buffers.
fn write_fasta(path: &Path, records: &[(&str, usize, &[(usize, String)])]) {
    let mut out = BufWriter::with_capacity(1 << 20, File::create(path).unwrap());
    let mut poly_a_block = Vec::with_capacity(81 * 1000);
    for _ in 0..1000 {
        poly_a_block.extend_from_slice(&[b'A'; 80]);
        poly_a_block.push(b'\n');
    }
    for &(name, length, pieces) in records {
        writeln!(out, ">{name}").unwrap();
        if pieces.is_empty() {
            let full = length / 80;
            let rest = length % 80;
            for _ in 0..(full / 1000) {
                out.write_all(&poly_a_block).unwrap();
            }
            out.write_all(&poly_a_block[..(full % 1000) * 81]).unwrap();
            if rest > 0 {
                out.write_all(&poly_a_block[..rest]).unwrap();
                writeln!(out).unwrap();
            }
            continue;
        }
        let mut seq = vec![b'A'; length];
        for (residue, dna) in pieces {
            let start = 3 * residue;
            seq[start..start + dna.len()].copy_from_slice(dna.as_bytes());
        }
        for line in seq.chunks(80) {
            out.write_all(line).unwrap();
            out.write_all(b"\n").unwrap();
        }
    }
    out.flush().unwrap();
}

/// Formats ASN.1 text masks (`algo-id 2`, `dust`, `"m4bt3"`) for `local_ids`.
fn format_asnt(local_ids: &[u32]) -> String {
    let blocks: Vec<String> = local_ids
        .iter()
        .map(|&sid| {
            let items: Vec<String> = MASKS
                .iter()
                .map(|&(from, to)| {
                    format!(
                        "        {{\n          from {from},\n          to {to},\n          id local id {sid}\n        }}"
                    )
                })
                .collect();
            format!("      packed-int {{\n{}\n      }}", items.join(",\n"))
        })
        .collect();
    let mut text = String::from(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program dust,\n  algo-options \"m4bt3\",\n  \
         masks {\n    masks {\n",
    );
    writeln!(text, "{}", blocks.join(",\n")).unwrap();
    text.push_str("    },\n    more FALSE\n  }\n}\n");
    text
}

/// Runs `losat` with `LOSAT_*` and NCBI search-changing env vars stripped, preserving
/// `LOSAT_BLASTX_WORKER_LOG` for `blastx-worker-probe` builds.
fn losat(dir: &Path) -> Command {
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
    command
        .current_dir(dir)
        .env("LOSAT_BLASTX_WORKER_LOG", dir.join("worker.log"));
    command
}

fn build_db(dir: &Path, name: &str, keep_fasta: bool) {
    let fasta_name = format!("{name}.fna");
    let asnt_name = format!("{name}.asnt");
    let output = losat(dir)
        .args([
            "makeblastdb",
            "-in",
            &fasta_name,
            "-dbtype",
            "nucl",
            "-out",
            name,
            "-mask_data",
            &asnt_name,
        ])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb {name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    if !keep_fasta {
        std::fs::remove_file(dir.join(&fasta_name)).unwrap();
    }
}

fn run_search(
    dir: &Path,
    program: &str,
    db: &str,
    extra: &[&str],
    outfmt: &str,
    threads: u32,
) -> Output {
    let query = if program == "tblastn" {
        "q.faa"
    } else {
        "qx.fna"
    };
    losat(dir)
        .arg(program)
        .args([
            "-query",
            query,
            "-db",
            db,
            "-db_soft_mask",
            "11",
            "-seg",
            "no",
        ])
        .args(extra)
        .args(["-outfmt", outfmt, "-num_threads", &threads.to_string()])
        .output()
        .expect("run losat search")
}

fn expected_fixture(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastn_db_soft_lead_start")
        .join(format!("{name}.6.out"));
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

fn check_row_case(
    dir: &Path,
    fixture: &str,
    program: &str,
    db: &str,
    extra: &[&str],
    outfmt: &str,
    failed: &mut Vec<String>,
) {
    let want = expected_fixture(fixture);
    for threads in [1, 4] {
        let output = run_search(dir, program, db, extra, outfmt, threads);
        if !output.status.success() || !output.stderr.is_empty() || output.stdout != want {
            failed.push(format!(
                "{fixture} at {threads} threads: status={} stderr={} stdout={}",
                output.status,
                String::from_utf8_lossy(&output.stderr),
                String::from_utf8_lossy(&output.stdout),
            ));
        }
    }
}

fn check_refusal_case(dir: &Path, db: &str, oid: u32, failed: &mut Vec<String>) {
    let want_msg = format!("TBLASTN: database sequence {oid}, frame 1: {REFUSAL_SUFFIX}");
    for threads in [1, 4] {
        let output = run_search(dir, "tblastn", db, &["-window_size", "0"], "6", threads);
        let stderr = String::from_utf8_lossy(&output.stderr);
        let stdout = String::from_utf8_lossy(&output.stdout);
        if output.status.code() != Some(1)
            || !output.stdout.is_empty()
            || !stderr.contains(&want_msg)
        {
            failed.push(format!(
                "refusal {db} at {threads} threads: status={} stdout={stdout} stderr={stderr}",
                output.status,
            ));
        }
    }
}

fn test_work_dir(tag: &str) -> (PathBuf, bool) {
    match std::env::var_os("TBLASTN_DB_SOFT_LEAD_START_DIR") {
        Some(dir) => (PathBuf::from(dir), true),
        None => (
            std::env::temp_dir().join(format!(
                "losat-tblastn-db-soft-lead-start-{tag}-{}",
                std::process::id()
            )),
            false,
        ),
    }
}

fn write_queries(dir: &Path) {
    std::fs::write(
        dir.join("q.faa"),
        format!(">q M4b-T3 probe query\n{QUERY_AA}\n"),
    )
    .unwrap();
    std::fs::write(
        dir.join("qx.fna"),
        format!(">qx coding DNA of q\n{}\n", protein_dna(QUERY_AA)),
    )
    .unwrap();
}

const W0_ARGS: &[&str] = &["-window_size", "0"];
const W0_CBS0_ARGS: &[&str] = &["-window_size", "0", "-comp_based_stats", "0"];
const UNG_W0_ARGS: &[&str] = &["-ungapped", "-comp_based_stats", "0", "-window_size", "0"];

fn with_work_dir(tag: &str, body: impl FnOnce(&Path, bool, &mut Vec<String>)) {
    let (dir, kept) = test_work_dir(tag);
    std::fs::create_dir_all(&dir).unwrap();
    write_queries(&dir);
    let mut failed = Vec::new();
    body(&dir, kept, &mut failed);
    if !kept {
        std::fs::remove_dir_all(&dir).unwrap();
    }
    assert!(failed.is_empty(), "{failed:#?}");
}

/// Always-run test (~31 MB peak): searches `q1a1` (11 aa at -4, gapped start +1, standard and rich
/// columns) and `q1a2` (11 aa at -5, gapped start 0) against frozen NCBI BLAST+ 2.17.0 outputs,
/// and checks refusal of `q1c1` (middle at -1) and `q1d1` (12 aa), at 1 and 4 threads.
#[test]
fn short_pre_chunk_alignments_match_ncbi_and_crashes_are_refused() {
    with_work_dir("q1", |dir, kept, failed| {
        let asnt_one = format_asnt(&[1]);
        for &(case_name, core, offset) in &[
            ("a1", CORE_A1, -4isize),
            ("a2", CORE_A2, -5isize),
            ("c1", CORE_C1, -6isize),
            ("d1", CORE_D1, -4isize),
        ] {
            let db = format!("q1{case_name}");
            let dna = protein_dna(core);
            let lead_res = (CHUNK2_START as isize + offset) as usize;
            let pieces = [(lead_res, dna.clone()), (CONTROL_AT, dna)];
            write_fasta(
                &dir.join(format!("{db}.fna")),
                &[(case_name, LONG_LEN, &pieces)],
            );
            std::fs::write(dir.join(format!("{db}.asnt")), &asnt_one).unwrap();
            build_db(dir, &db, kept);
        }

        for &(fixture, db, outfmt) in &[
            ("q1a1.tblastn_w0", "q1a1", "6"),
            ("q1a1.tblastn_w0rich", "q1a1", RICH_OUTFMT),
            ("q1a2.tblastn_w0", "q1a2", "6"),
        ] {
            check_row_case(dir, fixture, "tblastn", db, W0_ARGS, outfmt, failed);
        }
        for db in ["q1c1", "q1d1"] {
            check_refusal_case(dir, db, 0, failed);
        }
    });
}

/// Ignored test (~284 Mbp per database; ~10 s in release mode on the workstation):
/// `C0` (diagonal table clears after chunk 1 of `L`, skipping pre-chunk words so gapped
/// `-window_size 0` completes and `tblastn`/`tblastx` rows match NCBI's `-num_threads 1` output)
/// vs `C2` (last filler one base shorter: no clear after chunk 1, so `tblastn -window_size 0` is
/// refused while `-ungapped` and `tblastx` match NCBI's `-num_threads 1` output), at 1 and 4
/// threads.
#[test]
#[ignore = "writes two 284 Mbp databases; run explicitly on the workstation"]
fn one_thread_diagonal_clear_boundary_matches_ncbi_at_every_thread_count() {
    with_work_dir("q2", |dir, kept, failed| {
        let dna_a1 = protein_dna(CORE_A1);
        let dna_e1 = protein_dna(CORE_E1);
        let dna_b1 = protein_dna(CORE_B1);
        let l_pieces = [
            (CHUNK2_START - 4, dna_a1.clone()),
            (CONTROL_AT, dna_a1),
            (CHUNK2_START - 1, dna_e1.clone()),
            (CHUNK2_START + 30_000, dna_e1),
            (CHUNK2_START - 30, dna_b1),
        ];
        let filler_names: Vec<String> = (2..17).map(|k| format!("f{k}")).collect();
        let asnt_l = format_asnt(&[19]);

        for &(db, flast_len) in &[("C0", 3_434_890usize), ("C2", 3_434_889usize)] {
            let mut records: Vec<(&str, usize, &[(usize, String)])> = Vec::with_capacity(19);
            records.push(("t0", 10, &[]));
            records.push(("f1", 40_000_000, &[]));
            for name in &filler_names {
                records.push((name.as_str(), 15_000_000, &[]));
            }
            records.push(("flast", flast_len, &[]));
            records.push(("L", LONG_LEN, &l_pieces));
            write_fasta(&dir.join(format!("{db}.fna")), &records);
            std::fs::write(dir.join(format!("{db}.asnt")), &asnt_l).unwrap();
            build_db(dir, db, kept);
        }

        let cases = [
            ("C0.tblastn_ungw0", "tblastn", "C0", UNG_W0_ARGS, "6"),
            (
                "C0.tblastn_w0cbs0rich",
                "tblastn",
                "C0",
                W0_CBS0_ARGS,
                RICH_OUTFMT,
            ),
            ("C0.tblastn_w0", "tblastn", "C0", W0_ARGS, "6"),
            ("C0.tblastx_w0", "tblastx", "C0", W0_ARGS, "6"),
            ("C2.tblastn_ungw0", "tblastn", "C2", UNG_W0_ARGS, "6"),
            ("C2.tblastx_w0", "tblastx", "C2", W0_ARGS, "6"),
        ];
        for (fixture, program, db, extra, outfmt) in cases {
            check_row_case(dir, fixture, program, db, extra, outfmt, failed);
        }
        check_refusal_case(dir, "C2", 18, failed);
    });
}
