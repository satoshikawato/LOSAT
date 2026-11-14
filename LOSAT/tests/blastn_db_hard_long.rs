//! BLASTN's `-db_hard_mask` on database sequences longer than NCBI's chunk (`MAX_DBSEQ_LEN`,
//! 5,000,000 letters): the preliminary stage cuts each hard range into chunks
//! (`s_GetNextSubjectChunk`, blast_engine.c:221-310) and merges each chunk's list with an overlap
//! that NCBI reads one range past its list after the last chunk (blast_engine.c:579-584; series NP,
//! M4a-2b decision B-d, review T3).
//!
//! The database is generated here with a fixed PRNG (`write_inputs`): `hard1` (5,300,000 letters)
//! with the one-letter mask 1001 (review F3: its last hard range is longer than `MAX_DBSEQ_LEN`, so
//! the chunk after the mask is a full chunk and the merge across the mask reads the overlap in
//! bounds; its final chunk merges with the full chunk before it in the same range), `hard2`
//! (5,400,000 letters) with the masks 1,000-1,099 and 5,100,002-5,100,049 (its middle range is cut
//! into chunks of `MAX_DBSEQ_LEN` letters, its last chunk starts two letters before its range),
//! and `hard3` (3,000 letters) with adjacent, one-letter and three-letter masks and a mask that
//! runs to its end after other masks. The queries are copies of regions across the masks and
//! the chunk boundaries with one letter in 40 changed. The database is made by
//! `losat makeblastdb -mask_data` from ASN.1 text and searched by the four tasks at 1 and 4
//! threads. The expected rows are NCBI BLAST+ 2.17.0's (`tests/fixtures/db_search_long/`): NCBI's
//! makeblastdb and blastn `-num_threads 1` on the same inputs (`DBMASK_HARD_LONG_DIR=<dir>` keeps
//! them there; README.md of that fixture directory).
//!
//! The run writes about 14 MB into a temporary directory.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The binary with no `LOSAT_*` variable that chooses a code path (`LOSAT_BLASTX_WORKER_LOG`
/// stays, as in `blastn_db_search.rs`).
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// xorshift64*: the inputs' letters and mutations, the same on every host.
struct Prng(u64);

impl Prng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn letter(&mut self) -> u8 {
        b"ACGT"[(self.next() >> 62) as usize]
    }
}

/// The sequences' lengths.
const LENGTHS: [usize; 3] = [5_300_000, 5_400_000, 3_000];

/// The stored masks of each sequence (local ids 1 to 3), `[from, to]` inclusive, sorted.
const MASKS: [(u32, &[(usize, usize)]); 3] = [
    (1, &[(1_001, 1_001)]),
    (2, &[(1_000, 1_099), (5_100_002, 5_100_049)]),
    (
        3,
        &[
            (100, 199),
            (200, 299),
            (1_001, 1_001),
            (1_500, 1_502),
            (2_990, 2_999),
        ],
    ),
];

/// The query regions `(sequence, from, to)` (0-based, end exclusive): across the masks (one HSP
/// over `hard1`'s mask, merged across it in NCBI's run), across the chunk boundaries of the long
/// ranges (`hard1`'s final chunk starts at 5,000,900, the second chunk of `hard2`'s middle range
/// at 5,001,000), across `hard2`'s second mask, and the whole of `hard3`.
const QUERIES: [(usize, usize, usize); 6] = [
    (0, 500, 1_700),
    (0, 5_000_000, 5_002_000),
    (1, 500, 1_700),
    (1, 5_000_000, 5_002_000),
    (1, 5_099_000, 5_101_000),
    (2, 0, 3_000),
];

/// Writes `hl.fna` (the sequences), `hl_masks.asnt` (their masks) and `hlq.fna` (the queries:
/// copies of the regions with one letter in 40 changed) into `dir`.
fn write_inputs(dir: &Path) {
    let mut prng = Prng(20_261_116);
    let sequences: Vec<Vec<u8>> = LENGTHS
        .iter()
        .map(|&length| (0..length).map(|_| prng.letter()).collect())
        .collect();
    let mut fasta = String::new();
    for (i, sequence) in sequences.iter().enumerate() {
        writeln!(fasta, ">hard{} hard mask chunk test", i + 1).unwrap();
        for line in sequence.chunks(80) {
            fasta.push_str(std::str::from_utf8(line).unwrap());
            fasta.push('\n');
        }
    }
    std::fs::write(dir.join("hl.fna"), fasta).unwrap();
    let mut masks = String::from(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program dust,\n  algo-options \"hard \
         chunks\",\n  masks {\n    masks {\n",
    );
    let locations: Vec<String> = MASKS
        .iter()
        .map(|(id, ranges)| {
            let items: Vec<String> = ranges
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
    masks.push_str(&locations.join(",\n"));
    masks.push_str("\n    },\n    more FALSE\n  }\n}\n");
    std::fs::write(dir.join("hl_masks.asnt"), masks).unwrap();
    let mut queries = String::new();
    for (n, &(sequence, from, to)) in QUERIES.iter().enumerate() {
        let mut letters = sequences[sequence][from..to].to_vec();
        for letter in letters.iter_mut() {
            if prng.next().is_multiple_of(40) {
                *letter = prng.letter();
            }
        }
        writeln!(queries, ">hq{n} hard{} {from}-{to}", sequence + 1).unwrap();
        queries.push_str(std::str::from_utf8(&letters).unwrap());
        queries.push('\n');
    }
    std::fs::write(dir.join("hlq.fna"), queries).unwrap();
}

/// The frozen NCBI rows of a task.
fn expected(task: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/db_search_long")
        .join(format!("dbmask_hard_long.{task}.6.out"));
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The tabular columns: NCBI's standard twelve and the identities.
const OUTFMT: &str = "6 std nident";

#[test]
fn hard_masked_sequences_over_ncbi_s_chunk_merge_their_chunks_as_ncbi() {
    let dir = match std::env::var_os("DBMASK_HARD_LONG_DIR") {
        Some(dir) => PathBuf::from(dir),
        None => {
            std::env::temp_dir().join(format!("losat-blastn-db-hard-long-{}", std::process::id()))
        }
    };
    std::fs::create_dir_all(&dir).unwrap();
    write_inputs(&dir);
    let db = dir.join("hl").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(dir.join("hl.fna"))
        .args(["-dbtype", "nucl", "-out", &db, "-mask_data"])
        .arg(dir.join("hl_masks.asnt"))
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut failed = Vec::new();
    for task in ["megablast", "dc-megablast", "blastn", "blastn-short"] {
        let expected = expected(task);
        for threads in [1, 4] {
            let mut search = losat();
            search
                .arg("blastn")
                .arg("-query")
                .arg(dir.join("hlq.fna"))
                .args(["-db", &db, "-db_hard_mask", "11", "-task", task])
                .args(["-outfmt", OUTFMT, "-num_threads", &threads.to_string()]);
            if task == "blastn-short" {
                search.args(["-evalue", "1e-5"]);
            }
            let output = search.output().expect("run losat blastn");
            assert!(
                output.status.success() && output.stderr.is_empty(),
                "{task} at {threads} threads: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            if output.stdout != expected {
                failed.push(format!("{task} at {threads} threads"));
            }
        }
    }
    if std::env::var_os("DBMASK_HARD_LONG_DIR").is_none() {
        std::fs::remove_dir_all(&dir).unwrap();
    }
    assert!(failed.is_empty(), "rows differ from NCBI's: {failed:?}");
}
