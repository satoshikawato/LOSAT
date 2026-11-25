//! BLASTP and BLASTX `-db_hard_mask` on protein database sequences with hard ranges longer than
//! NCBI's chunk (`MAX_DBSEQ_LEN`, 5,000,000 residues), out-of-order masks whose parts overlap
//! (`PartsOutOfOrder`), one-word short chunks whose HSPs merge with another chunk's, and an
//! unmasked > 5,000,000-residue sibling sequence (series NP, M4b-3d).
//!
//! The database is generated deterministically (`write_inputs`) and built with
//! `losat makeblastdb -mask_data`, then searched by `blastp` and `blastx` at 1 and 4 threads
//! across default (CBS 2), `-comp_based_stats 0`, `-ungapped -comp_based_stats 0`,
//! `-window_size 0`, and `-subject_besthit`. Expected outputs in
//! `tests/fixtures/db_search_long/dbmask_prot_hard_long.*.6.out` were frozen from official NCBI
//! BLAST+ 2.17.0 on the same inputs (`DBMASK_PROT_HARD_LONG_DIR=<dir>` keeps them there).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The binary with no `LOSAT_*` variable that chooses a code path (`LOSAT_BLASTX_WORKER_LOG`
/// stays, as in `blastn_db_hard_long.rs`).
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

const ALPHABET: &[u8] = b"WYFHCMIVKRDEPNQST";

fn make_motif(mut seed: u32, length: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(length);
    for _ in 0..length {
        seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12_345) & 0x7fff_ffff;
        out.push(ALPHABET[(seed as usize) % ALPHABET.len()]);
    }
    out
}

fn prot_to_dna(prot: &[u8]) -> String {
    let mut dna = String::with_capacity(prot.len() * 3);
    for &aa in prot {
        let codon = match aa {
            b'A' => "GCT",
            b'C' => "TGC",
            b'D' => "GAC",
            b'E' => "GAA",
            b'F' => "TTC",
            b'G' => "GGT",
            b'H' => "CAC",
            b'I' => "ATC",
            b'K' => "AAA",
            b'L' => "CTG",
            b'M' => "ATG",
            b'N' => "AAC",
            b'P' => "CCG",
            b'Q' => "CAG",
            b'R' => "CGT",
            b'S' => "TCT",
            b'T' => "ACC",
            b'V' => "GTT",
            b'W' => "TGG",
            b'Y' => "TAC",
            other => panic!("unexpected residue byte {other}"),
        };
        dna.push_str(codon);
    }
    dna
}

fn write_fasta(path: &Path, records: &[(&str, &[u8])]) {
    let mut out = String::new();
    for (name, seq) in records {
        writeln!(out, ">{name}").unwrap();
        for line in seq.chunks(80) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
    }
    std::fs::write(path, out).unwrap();
}

/// Writes `phl.faa`, `phl_masks.asnt`, `phlq.faa`, and `phlq.fna` into `dir`.
fn write_inputs(dir: &Path) {
    // Long-range subject template (5,000,300 aa):
    // - motif m0 at [500, 580)
    // - bridge motif at [4_999_840, 5_000_060) straddling 5,000,000 chunk boundaries
    // - motif m1 at [5_000_090, 5_000_170)
    let m0 = make_motif(301, 80);
    let m_bridge = make_motif(302, 220);
    let m1 = make_motif(303, 80);
    let mut long_seq = vec![b'A'; 5_000_300];
    long_seq[500..580].copy_from_slice(&m0);
    long_seq[4_999_840..5_000_060].copy_from_slice(&m_bridge);
    long_seq[5_000_090..5_000_170].copy_from_slice(&m1);

    // Out-of-order subject (400 aa): motif at [80, 280).
    let m_po = make_motif(101, 200);
    let mut po_seq = vec![b'A'; 400];
    po_seq[80..280].copy_from_slice(&m_po);

    // OneWord merge subject (240 aa): motif at [20, 120) with WWW at [70, 73).
    let mut m_ow = make_motif(202, 100);
    m_ow[50..53].copy_from_slice(b"WWW");
    let mut ow_seq = vec![b'A'; 240];
    ow_seq[20..120].copy_from_slice(&m_ow);

    write_fasta(
        &dir.join("phl.faa"),
        &[
            ("phl_lr_mid", &long_seq),
            ("phl_lr_po", &long_seq),
            ("phl_lr_plain", &long_seq),
            ("phl_po", &po_seq),
            ("phl_ow", &ow_seq),
        ],
    );

    // Masks (0-based inclusive `from..=to`):
    // 1 (phl_lr_mid):   [10, 20), [5_000_250, 5_000_280) -> LongRange in middle
    // 2 (phl_lr_po):    [620, 650), [50, 80), [5_000_250, 5_000_280) -> PartsOutOfOrder + LongRange
    // 3 (phl_lr_plain): no masks -> unmasked > 5,000,000 sibling in a -db_hard_mask database
    // 4 (phl_po):       [250, 280), [100, 130), [320, 340) -> PartsOutOfOrder overlap merge
    // 5 (phl_ow):       [40, 60), [55, 70), [72, 90), [50, 68), [200, 210) -> OneWord overlap merge
    let masks_entries: &[(u32, &[(usize, usize)])] = &[
        (1, &[(10, 19), (5_000_250, 5_000_279)]),
        (2, &[(620, 649), (50, 79), (5_000_250, 5_000_279)]),
        (4, &[(250, 279), (100, 129), (320, 339)]),
        (5, &[(40, 59), (55, 69), (72, 89), (50, 67), (200, 209)]),
    ];
    let mut asnt = String::from(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program seg,\n  algo-options \"m4b3d \
         long hard\",\n  masks {\n    masks {\n",
    );
    let blocks: Vec<String> = masks_entries
        .iter()
        .map(|(id, ranges)| {
            let items: Vec<String> = ranges
                .iter()
                .map(|(from, to)| {
                    format!(
                        "        {{\n          from {from},\n          to {to},\n          id \
                         local id {id}\n        }}"
                    )
                })
                .collect();
            format!("      packed-int {{\n{}\n      }}", items.join(",\n"))
        })
        .collect();
    asnt.push_str(&blocks.join(",\n"));
    asnt.push_str("\n    },\n    more FALSE\n  }\n}\n");
    std::fs::write(dir.join("phl_masks.asnt"), asnt).unwrap();

    let mut q_bridge_gap2 = Vec::with_capacity(218);
    q_bridge_gap2.extend_from_slice(&m_bridge[..100]);
    q_bridge_gap2.extend_from_slice(&m_bridge[102..]);

    let q_po_exact = m_po[30..180].to_vec();
    let mut q_po_gap2 = Vec::with_capacity(148);
    q_po_gap2.extend_from_slice(&m_po[30..95]);
    q_po_gap2.extend_from_slice(&m_po[97..180]);
    let mut q_po_gap12 = Vec::with_capacity(138);
    q_po_gap12.extend_from_slice(&m_po[30..90]);
    q_po_gap12.extend_from_slice(&m_po[102..180]);

    let prot_queries: Vec<(&str, &[u8])> = vec![
        ("q_bridge", &m_bridge),
        ("q_bridge_gap2", &q_bridge_gap2),
        ("q_m0", &m0),
        ("q_po_exact", &q_po_exact),
        ("q_po_gap2", &q_po_gap2),
        ("q_po_gap12", &q_po_gap12),
        ("q_ow", &m_ow),
    ];
    write_fasta(&dir.join("phlq.faa"), &prot_queries);

    let dna_owned: Vec<(&str, Vec<u8>)> = prot_queries
        .iter()
        .map(|&(name, seq)| (name, prot_to_dna(seq).into_bytes()))
        .collect();
    let dna_refs: Vec<(&str, &[u8])> = dna_owned
        .iter()
        .map(|(name, seq)| (*name, seq.as_slice()))
        .collect();
    write_fasta(&dir.join("phlq.fna"), &dna_refs);
}

fn expected(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/db_search_long")
        .join(format!("dbmask_prot_hard_long.{name}.6.out"));
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

const OUTFMT: &str = "6 std nident positive btop";

struct RunDir {
    path: PathBuf,
    keep: bool,
}

impl Drop for RunDir {
    fn drop(&mut self) {
        if !self.keep {
            let _ = std::fs::remove_dir_all(&self.path);
        }
    }
}

#[test]
fn protein_hard_masked_long_ranges_and_out_of_order_chunks_match_ncbi() {
    let run_dir = match std::env::var_os("DBMASK_PROT_HARD_LONG_DIR") {
        Some(dir) => RunDir {
            path: PathBuf::from(dir),
            keep: true,
        },
        None => RunDir {
            path: std::env::temp_dir()
                .join(format!("losat-protein-db-hard-long-{}", std::process::id())),
            keep: false,
        },
    };
    let dir = run_dir.path.as_path();
    std::fs::create_dir_all(dir).unwrap();
    write_inputs(dir);
    let db = dir.join("phl").to_str().unwrap().to_string();
    let made = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(dir.join("phl.faa"))
        .args(["-dbtype", "prot", "-out", &db, "-mask_data"])
        .arg(dir.join("phl_masks.asnt"))
        .output()
        .expect("run losat makeblastdb");
    assert!(
        made.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&made.stderr)
    );

    let modes: &[(&str, &[&str])] = &[
        ("def", &[]),
        ("cbs0", &["-comp_based_stats", "0"]),
        ("ungapped", &["-ungapped", "-comp_based_stats", "0"]),
        ("w0", &["-window_size", "0"]),
        ("besthit", &["-subject_besthit"]),
    ];

    let mut failed = Vec::new();
    for &(prog, query_file) in &[("blastp", "phlq.faa"), ("blastx", "phlq.fna")] {
        for &(mode, extra) in modes {
            let label = format!("{prog}.{mode}");
            let want = expected(&label);
            for threads in [1, 4] {
                let output = losat()
                    .arg(prog)
                    .arg("-query")
                    .arg(dir.join(query_file))
                    .args(["-db", &db, "-db_hard_mask", "21", "-outfmt", OUTFMT])
                    .args(extra)
                    .args(["-num_threads", &threads.to_string()])
                    .output()
                    .expect("run losat");
                if !output.status.success() || !output.stderr.is_empty() {
                    failed.push(format!(
                        "{label} t{threads}: status {:?}, stderr {:?}",
                        output.status.code(),
                        String::from_utf8_lossy(&output.stderr)
                    ));
                } else if output.stdout != want {
                    failed.push(format!(
                        "{label} t{threads}: {} bytes vs expected {}",
                        output.stdout.len(),
                        want.len()
                    ));
                }
            }
        }
    }
    assert!(failed.is_empty(), "mismatches:\n{}", failed.join("\n"));
}
