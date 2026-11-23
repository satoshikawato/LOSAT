//! TBLASTX and TBLASTN with `-db_soft_mask` where the masks of a database sequence are stored out
//! of order on a frame longer than NCBI's chunk (`MAX_DBSEQ_LEN`, 5,000,000 residues): a later
//! range of the second chunk starts before the chunk, so the scan reads the frame's letters before
//! it and finds words there (`common::translated_chunks::scan_lead`; series NP, M4b-T2).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:286-307
//! The test writes its inputs (the same letters, ids and masks as
//! `ncbi-parity-20261019/M/m4bt2/scripts/gen_lead_test.py`): 15,300,000 letters of A (their frames
//! read K and F, which give no word hit with the query) with copies of the coding DNA of a query
//! protein Q, and masks in an order that gives the second chunk of frame +1 (from residue
//! 4,999,900) a last range from residue 1,100,000 to the frame's end:
//! - `lead`: `lead1` with A, D and Q's residues 6 to 55 at residues 4,999,898 to 4,999,949: a
//!   one-hit word at the chunk's offset -1 whose first letter scores below 0, so its HSP starts at
//!   offset 0 and its seed lies at -1 (the gapped start NCBI takes from that seed, a `Uint4` wrapped
//!   below 0); `lead2` with the masks of `LAST_MASKS`, whose second chunk's first range is raised
//!   to one word on the frame's last letter (no hit there) before the range from residue
//!   1,100,000, which holds Q at residue 5,029,900.
//! - `before`: a short sequence (OID 0), then Q at residues 4,999,805 (the chunk's offsets -95 to
//!   -35, an HSP wholly before the chunk) and 4,999,875 (offsets -25 to 35, across its start).
//!
//! The rows are NCBI BLAST+ 2.17.0's (`tests/fixtures/tblastx_db_soft_lead/`, README.md there), at 1
//! and 4 threads. TBLASTN `-window_size 0` on `before` makes NCBI's gapped start read past its
//! buffer (a segmentation fault at 1 and 2 threads): LOSAT refuses it (PD-LOSAT-NCBI-DEFECTS, rule
//! 3) and names the sequence and the frame. The run writes about 60 MB into a temporary directory.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The long sequences' length: frame +1 of 5,100,000 residues, in two chunks.
const LENGTH: usize = 15_300_000;
/// The query protein Q and its coding DNA.
const Q: &str = "GMLHVWCHMFLLEMKFEKLYQIKYIRVYVIMIKRYYRYWMCIHALLKELPVCWKMYYYRD";
const QX: &str = concat!(
    "GGTATGCTTCATGTGTGGTGTCACATGTTTCTACTTGAAATGAAGTTTGAAAAGTTGTACCAAATAAAGTACATACGCGTGTACGTAATCAT",
    "GATCAAGAGATACTACCGCTATTGGATGTGCATACATGCCTTGTTAAAAGAACTGCCGGTTTGCTGGAAAATGTACTACTATCGTGAT"
);
/// The masks of `lead1` and `before1`, `[from, to]` in bases, in the order makeblastdb stores them.
const MASKS: [(usize, usize); 3] = [
    (14_999_400, 14_999_849),
    (14_999_880, 15_179_999),
    (3_000_000, 3_299_999),
];
/// The masks of `lead2`.
const LAST_MASKS: [(usize, usize); 3] = [
    (14_996_700, 15_299_996),
    (15_299_998, 15_299_998),
    (3_000_000, 3_299_999),
];

/// The searches: the fixture's name, the program, the database and the options.
const SEARCHES: [(&str, &str, &str, &[&str]); 7] = [
    ("lead.tblastx_w0", "tblastx", "lead", &["-window_size", "0"]),
    ("lead.tblastn_w0", "tblastn", "lead", &["-window_size", "0"]),
    ("lead.tblastn", "tblastn", "lead", &[]),
    (
        "before.tblastx_w0",
        "tblastx",
        "before",
        &["-window_size", "0"],
    ),
    ("before.tblastx", "tblastx", "before", &[]),
    ("before.tblastn", "tblastn", "before", &[]),
    (
        "before.tblastn_ungapped_w0",
        "tblastn",
        "before",
        &["-ungapped", "-comp_based_stats", "0", "-window_size", "0"],
    ),
];

/// A long sequence of A with the DNA of `pieces` at their frame +1 residues.
fn sequence(pieces: &[(usize, &str)]) -> Vec<u8> {
    let mut letters = vec![b'A'; LENGTH];
    for &(residue, dna) in pieces {
        letters[3 * residue..3 * residue + dna.len()].copy_from_slice(dna.as_bytes());
    }
    letters
}

fn fasta(records: &[(&str, &[u8])]) -> Vec<u8> {
    let mut text = Vec::with_capacity(2 * LENGTH);
    for (name, letters) in records {
        text.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in letters.chunks(80) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    text
}

/// The ASN.1 text of the masks of the sequences of local ids `by_id`.
fn masks(by_id: &[(u32, &[(usize, usize)])]) -> String {
    let blocks: Vec<String> = by_id
        .iter()
        .map(|(id, masks)| {
            let items: Vec<String> = masks
                .iter()
                .map(|(from, to)| {
                    format!(
                        "        {{\n          from {from},\n          to {to},\n          id local id {id}\n        }}"
                    )
                })
                .collect();
            format!("      packed-int {{\n{}\n      }}", items.join(",\n"))
        })
        .collect();
    let mut text = String::from(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program dust,\n  algo-options \"lead\",\n  \
         masks {\n    masks {\n",
    );
    writeln!(text, "{}", blocks.join(",\n")).unwrap();
    text.push_str("    },\n    more FALSE\n  }\n}\n");
    text
}

fn write_inputs(dir: &Path) {
    std::fs::write(dir.join("q.faa"), format!(">q\n{Q}\n")).unwrap();
    std::fs::write(dir.join("qx.fna"), format!(">qx\n{QX}\n")).unwrap();
    let lead_piece = format!("GCTGAT{}", &QX[15..165]);
    let lead = sequence(&[(4_999_898, lead_piece.as_str())]);
    let last = sequence(&[(5_029_900, QX)]);
    std::fs::write(
        dir.join("lead.fna"),
        fasta(&[("lead1", lead.as_slice()), ("lead2", last.as_slice())]),
    )
    .unwrap();
    std::fs::write(
        dir.join("lead_masks.asnt"),
        masks(&[(1, &MASKS[..]), (2, &LAST_MASKS[..])]),
    )
    .unwrap();
    let before = sequence(&[(4_999_805, QX), (4_999_875, QX)]);
    let decoy = vec![b'A'; 300];
    std::fs::write(
        dir.join("before.fna"),
        fasta(&[("decoy", decoy.as_slice()), ("before1", before.as_slice())]),
    )
    .unwrap();
    std::fs::write(dir.join("before_masks.asnt"), masks(&[(2, &MASKS[..])])).unwrap();
}

/// `losat` with no `LOSAT_*` variable but `LOSAT_BLASTX_WORKER_LOG` (a binary built with the
/// `blastx-worker-probe` feature needs that file name once a search has more than one thread; it
/// selects nothing; the log is in `dir`), and none of the other variables that change a search
/// (the list of `tests/tblastx_chunk_merge.rs`).
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

fn search(dir: &Path, program: &str, db: &str, args: &[&str], threads: u32) -> Output {
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
        .args(args)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

fn expected(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/tblastx_db_soft_lead")
        .join(format!("{name}.6.out"));
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

#[test]
fn chunks_whose_ranges_start_before_them_search_as_ncbi() {
    let kept = std::env::var_os("TBLASTX_DB_SOFT_LEAD_DIR");
    let dir = match &kept {
        Some(dir) => PathBuf::from(dir),
        None => {
            std::env::temp_dir().join(format!("losat-tblastx-db-soft-lead-{}", std::process::id()))
        }
    };
    std::fs::create_dir_all(&dir).unwrap();
    write_inputs(&dir);
    for db in ["lead", "before"] {
        let output = losat(&dir)
            .args([
                "makeblastdb",
                "-in",
                &format!("{db}.fna"),
                "-dbtype",
                "nucl",
            ])
            .args(["-out", db, "-mask_data", &format!("{db}_masks.asnt")])
            .output()
            .expect("run losat makeblastdb");
        assert!(
            output.status.success(),
            "makeblastdb {db}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let mut failed = Vec::new();
    for (name, program, db, args) in SEARCHES {
        let expected = expected(name);
        for threads in [1, 4] {
            let output = search(&dir, program, db, args, threads);
            if !output.status.success() || !output.stderr.is_empty() || output.stdout != expected {
                failed.push(format!(
                    "{name} at {threads} threads: {} {}",
                    output.status,
                    String::from_utf8_lossy(&output.stderr)
                ));
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3936-3953
    // The one-hit HSPs of `before` start before the chunk and are longer than HSP_MAX_WINDOW: NCBI's
    // gapped start reads the subject about 4 GiB past the frame.
    for threads in [1, 4] {
        let output = search(&dir, "tblastn", "before", &["-window_size", "0"], threads);
        let stderr = String::from_utf8_lossy(&output.stderr);
        if output.status.code() != Some(1)
            || !output.stdout.is_empty()
            || !stderr.contains("TBLASTN: database sequence 1, frame 1: an initial alignment of the chunk from residue 4999901 starts before the chunk")
        {
            failed.push(format!(
                "the refusal of before.tblastn_w0 at {threads} threads: {} {stderr}",
                output.status
            ));
        }
    }
    if kept.is_none() {
        std::fs::remove_dir_all(&dir).unwrap();
    }
    assert!(failed.is_empty(), "{failed:#?}");
}
