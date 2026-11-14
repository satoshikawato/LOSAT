//! BLASTN's `-db_soft_mask` on a database sequence longer than NCBI's chunk (`MAX_DBSEQ_LEN`,
//! 5,000,000 letters): the preliminary stage cuts it into chunks (`s_GetNextSubjectChunk`,
//! blast_engine.c:221-310) and gives each the soft ranges that overlap it, shifted and clipped
//! (279-307), in the order the masks are stored (series NP, M4a-2a; review 1 T4).
//!
//! The database is generated here with a fixed PRNG (`write_inputs`: a 5.6 Mbp sequence and a
//! short one, with hand-made masks that are unsorted, overlapping, nested, and across the chunks'
//! overlap at 4,999,900-5,000,000, `MASKS`), made by `losat makeblastdb -mask_data` from ASN.1
//! text, and searched by the four tasks at 1 and 4 threads; the masks remove hits that the same
//! searches without `-db_soft_mask` find. The expected rows are NCBI BLAST+ 2.17.0's
//! (`tests/fixtures/db_search_long/`): NCBI's makeblastdb and blastn `-num_threads 1` on the same
//! inputs (`DBMASK_LONG_DIR=<dir>` keeps them there; README.md of that fixture directory).
//!
//! The run writes about 7 MB into a temporary directory and takes a few seconds (3.5 s for the
//! eight searches with the test profile on one local machine).
//!
//! The second test (review 2 T1) searches a database of two volumes of 400 sequences of about
//! 1,000 letters, most with masks (`write_volume_inputs`), at 1, 4 and 8 threads: the preliminary
//! stage reads it in many blocks, split across the workers, and reads each block's masks once
//! (`SeqDb::read_mask_block`); the rows are NCBI's (`dbmask_vols.<task>.6.out`).

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

/// The length of the long sequence: two chunks of NCBI's preliminary stage.
const LONG: usize = 5_600_000;

/// The stored masks of each sequence (local ids 1 and 2), `[from, to]` inclusive, in the order
/// makeblastdb stores them. NCBI scans the complement in that order (blast_util.c:197-218), so in
/// `long1` the sorted ones mask their letters, the overlapping pair masks their union, the nested
/// mask gives back the end of the mask around it (from 5,000,080), and the unsorted last mask
/// gives back the two before it (an inverted range, then one to the end); the mask across
/// 4,999,900-5,000,000 lies in both chunks.
const MASKS: [(u32, &[(usize, usize)]); 2] = [
    (
        1,
        &[
            (500, 699),
            (1_000, 1_199),
            (3_000_000, 3_000_399),
            (3_000_200, 3_000_599),
            (4_999_850, 5_000_199),
            (5_000_050, 5_000_079),
            (5_400_000, 5_400_099),
            (5_599_950, 5_599_999),
            (5_500_000, 5_500_099),
        ],
    ),
    (2, &[(100, 199), (1_500, 1_799)]),
];

/// The query regions `(sequence, from, to)` (0-based, end exclusive): around and inside the
/// masks, across the chunk boundary, and one without masks.
const QUERIES: [(usize, usize, usize); 10] = [
    (0, 300, 1_500),
    (0, 4_999_400, 5_000_700),
    (0, 4_999_960, 5_000_040),
    (0, 2_999_700, 3_000_900),
    (0, 3_000_250, 3_000_330),
    (0, 5_399_800, 5_400_300),
    (0, 5_599_500, 5_600_000),
    (0, 2_000_000, 2_001_000),
    (1, 0, 2_000),
    (1, 1_520, 1_600),
];

/// Writes `long.fna` (the two sequences), `long_masks.asnt` (their masks) and `lq.fna` (the
/// queries: copies of the regions with one letter in 40 changed) into `dir`.
fn write_inputs(dir: &Path) {
    let mut prng = Prng(20_261_114);
    let sequences: Vec<Vec<u8>> = [LONG, 2_000]
        .iter()
        .map(|&length| (0..length).map(|_| prng.letter()).collect())
        .collect();
    let mut fasta = String::new();
    for (i, sequence) in sequences.iter().enumerate() {
        writeln!(fasta, ">long{} soft mask chunk test", i + 1).unwrap();
        for line in sequence.chunks(80) {
            fasta.push_str(std::str::from_utf8(line).unwrap());
            fasta.push('\n');
        }
    }
    std::fs::write(dir.join("long.fna"), fasta).unwrap();
    let mut masks = String::from(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program dust,\n  algo-options \"long \
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
    std::fs::write(dir.join("long_masks.asnt"), masks).unwrap();
    let mut queries = String::new();
    for (n, &(sequence, from, to)) in QUERIES.iter().enumerate() {
        let mut letters = sequences[sequence][from..to].to_vec();
        for letter in letters.iter_mut() {
            if prng.next().is_multiple_of(40) {
                *letter = prng.letter();
            }
        }
        writeln!(queries, ">lq{n} long{} {from}-{to}", sequence + 1).unwrap();
        queries.push_str(std::str::from_utf8(&letters).unwrap());
        queries.push('\n');
    }
    std::fs::write(dir.join("lq.fna"), queries).unwrap();
}

/// The frozen NCBI rows of a task.
fn expected(task: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/db_search_long")
        .join(format!("dbmask_long.{task}.6.out"));
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

#[test]
fn chunked_database_sequences_scan_their_soft_ranges_as_ncbi() {
    let dir = match std::env::var_os("DBMASK_LONG_DIR") {
        Some(dir) => PathBuf::from(dir),
        None => {
            std::env::temp_dir().join(format!("losat-blastn-db-soft-long-{}", std::process::id()))
        }
    };
    std::fs::create_dir_all(&dir).unwrap();
    write_inputs(&dir);
    let db = dir.join("long").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(dir.join("long.fna"))
        .args(["-dbtype", "nucl", "-out", &db, "-mask_data"])
        .arg(dir.join("long_masks.asnt"))
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
                .arg(dir.join("lq.fna"))
                .args(["-db", &db, "-db_soft_mask", "11", "-task", task])
                .args(["-outfmt", "6", "-num_threads", &threads.to_string()]);
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
    if std::env::var_os("DBMASK_LONG_DIR").is_none() {
        std::fs::remove_dir_all(&dir).unwrap();
    }
    assert!(failed.is_empty(), "rows differ from NCBI's: {failed:?}");
}

/// The sequences of the two-volume database, their length range and the volume size limit
/// (`-max_file_sz`) that cuts them into two volumes.
const VOLUME_SEQUENCES: usize = 400;
const VOLUME_MAX_FILE_SZ: &str = "64000";

/// Writes `vols.fna` (400 sequences of 900 to 1,099 letters), `vols_masks.asnt` (up to three
/// masks of 20 to 139 letters per sequence, a quarter of the sequences without any; a mask may
/// overlap the one before it) and `vq.fna` (24 queries over sequences of both volumes, with one
/// letter in 40 changed: copies of a sequence's first mask with 10 letters on each side, which
/// the masks keep from seeding, and of 400 letters) into `dir`.
fn write_volume_inputs(dir: &Path) {
    let mut prng = Prng(20_261_115);
    let sequences: Vec<Vec<u8>> = (0..VOLUME_SEQUENCES)
        .map(|_| {
            let length = 900 + (prng.next() % 200) as usize;
            (0..length).map(|_| prng.letter()).collect()
        })
        .collect();
    let mut fasta = String::new();
    for (i, sequence) in sequences.iter().enumerate() {
        writeln!(fasta, ">v{} volume block test", i + 1).unwrap();
        for line in sequence.chunks(80) {
            fasta.push_str(std::str::from_utf8(line).unwrap());
            fasta.push('\n');
        }
    }
    std::fs::write(dir.join("vols.fna"), fasta).unwrap();
    let mut locations = Vec::new();
    let mut first_masks = vec![None; sequences.len()];
    for (i, sequence) in sequences.iter().enumerate() {
        let count = (prng.next() % 4) as usize;
        let mut items = Vec::new();
        let mut from = 50 + (prng.next() % 250) as usize;
        for _ in 0..count {
            let to = from + 19 + (prng.next() % 120) as usize;
            if to >= sequence.len() {
                break;
            }
            first_masks[i].get_or_insert((from, to));
            items.push(format!(
                "        {{\n          from {from},\n          to {to},\n          id local id {}\n        }}",
                i + 1
            ));
            // The next mask starts up to 20 letters inside this one, or after it.
            from = (to + (prng.next() % 260) as usize).saturating_sub(20);
        }
        if !items.is_empty() {
            locations.push(format!(
                "      packed-int {{\n{}\n      }}",
                items.join(",\n")
            ));
        }
    }
    let mut masks = String::from(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program dust,\n  algo-options \"two \
         volumes\",\n  masks {\n    masks {\n",
    );
    masks.push_str(&locations.join(",\n"));
    masks.push_str("\n    },\n    more FALSE\n  }\n}\n");
    std::fs::write(dir.join("vols_masks.asnt"), masks).unwrap();
    let mut queries = String::new();
    for n in 0..24 {
        let sequence = n * (VOLUME_SEQUENCES / 24) + 3;
        let (from, end) = match first_masks[sequence] {
            Some((from, to)) if n % 3 != 2 => (from - 10, to + 11),
            _ => {
                let from = (prng.next() % 400) as usize;
                (from, from + 400)
            }
        };
        let mut letters = sequences[sequence][from..end].to_vec();
        for letter in letters.iter_mut() {
            if prng.next().is_multiple_of(40) {
                *letter = prng.letter();
            }
        }
        writeln!(queries, ">vq{n} v{} {from}", sequence + 1).unwrap();
        queries.push_str(std::str::from_utf8(&letters).unwrap());
        queries.push('\n');
    }
    std::fs::write(dir.join("vq.fna"), queries).unwrap();
}

#[test]
fn masked_searches_of_two_volumes_in_many_blocks_equal_ncbi_at_any_thread_count() {
    let dir = match std::env::var_os("DBMASK_VOLS_DIR") {
        Some(dir) => PathBuf::from(dir),
        None => {
            std::env::temp_dir().join(format!("losat-blastn-db-soft-vols-{}", std::process::id()))
        }
    };
    std::fs::create_dir_all(&dir).unwrap();
    write_volume_inputs(&dir);
    let db = dir.join("vols").to_str().unwrap().to_string();
    let output = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(dir.join("vols.fna"))
        .args([
            "-dbtype",
            "nucl",
            "-out",
            &db,
            "-max_file_sz",
            VOLUME_MAX_FILE_SZ,
        ])
        .arg("-mask_data")
        .arg(dir.join("vols_masks.asnt"))
        .output()
        .expect("run losat makeblastdb");
    assert!(
        output.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        dir.join("vols.01.nsq").exists() && !dir.join("vols.02.nsq").exists(),
        "the database has two volumes"
    );
    let mut failed = Vec::new();
    for task in ["megablast", "dc-megablast", "blastn", "blastn-short"] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/db_search_long")
            .join(format!("dbmask_vols.{task}.6.out"));
        let expected =
            std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        for threads in [1, 4, 8] {
            let mut search = losat();
            search
                .arg("blastn")
                .arg("-query")
                .arg(dir.join("vq.fna"))
                .args(["-db", &db, "-db_soft_mask", "11", "-task", task])
                .args(["-outfmt", "6", "-num_threads", &threads.to_string()]);
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
    if std::env::var_os("DBMASK_VOLS_DIR").is_none() {
        std::fs::remove_dir_all(&dir).unwrap();
    }
    assert!(failed.is_empty(), "rows differ from NCBI's: {failed:?}");
}
