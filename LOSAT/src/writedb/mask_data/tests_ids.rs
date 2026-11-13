//! Which sequence a mask reaches when the database is built with `-parse_seqids` (and without):
//! the probes of NCBI BLAST+ 2.17.0 behind `tests/seqdb/mask_ids/cases.tsv` (134 cases: 98 of
//! the handle probe (`seqid::mask_handle`), 4 on the order of two handles with one packed number
//! (`n63`-`n66`, M4a-1 fix A), 28 variant Seq-ids against one FASTA file, the
//! masks the real `dustmasker` / `segmasker` write with and without `-parse_seqids`), each built
//! with the library (`process_mask_data`, `BuildDatabase::add_fasta`, `end_build`) in format
//! version 5 and 4, then read back with `SeqDb::mask_ranges` (what `blastdbcmd %m` prints) and
//! compared with the rows NCBI printed.
//!
//! A case has one of four fates, the column `loss` of the table:
//! * `match`: LOSAT gives NCBI's rows (or, for `NOMATCH`, makeblastdb's message "No sequences
//!   matched any of the masks provided." and status 0);
//! * `refuse`: NCBI's answer depends on memory addresses (two different handles with the same
//!   `m_Packed`, `giim` handles among them, write different ranges on one strand; `range_set` has
//!   the rule), and LOSAT stops with its refusal (status 1, "losat makeblastdb:
//!   ...", no destructor message) that names the Seq-ids;
//! * `env:...`: NCBI's answer needs `OBJECTS_PACK_TEXTID=0` or `OBJECTS_PACK_GENERAL=0` in the
//!   environment, which LOSAT refuses (`blastinput::ncbi_environment`); not run here.

use std::collections::BTreeMap;
use std::path::PathBuf;

use super::tests_e2e::{run_case_with, Case, DESTRUCTOR_MESSAGE};
use crate::seqdb::{MaskRanges, SeqDb, SeqType};
use crate::writedb::tests::{seqdb_fixtures, TempDir};
use crate::writedb::volume::Indices;
use crate::writedb::writer::DbVersion;

pub(super) struct Row {
    pub(super) name: String,
    pub(super) fasta: String,
    pub(super) protein: bool,
    pub(super) parse: bool,
    pub(super) files: String,
    pub(super) algorithms: Vec<i32>,
    pub(super) loss: String,
    pub(super) expected: String,
}

pub(super) fn rows() -> Vec<Row> {
    let text = std::fs::read_to_string(dir().join("cases.tsv")).unwrap();
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            Row {
                name: f[0].to_string(),
                fasta: f[1].to_string(),
                protein: f[2] == "p",
                parse: f[3] == "1",
                files: f[4].to_string(),
                algorithms: f[5].split(',').map(|a| a.parse().unwrap()).collect(),
                loss: f[6].to_string(),
                expected: f[7].to_string(),
            }
        })
        .collect()
}

pub(super) fn dir() -> PathBuf {
    seqdb_fixtures().join("mask_ids")
}

/// The ranges as `blastdbcmd -outfmt %m` prints them.
fn printed(ranges: &MaskRanges) -> String {
    match ranges {
        MaskRanges::Ranges(ranges) if ranges.is_empty() => "none".to_string(),
        MaskRanges::Ranges(ranges) => ranges.iter().map(|(a, b)| format!("{a}-{b};")).collect(),
        other => panic!("{other:?}"),
    }
}

/// The rows of the sequences that have a mask, `oid|ranges of each algorithm`, joined by a space.
pub(super) fn masked_rows(base: &str, protein: bool, algorithms: &[i32]) -> String {
    let seqtype = if protein {
        SeqType::Protein
    } else {
        SeqType::Nucleotide
    };
    let db = SeqDb::open(base, seqtype).unwrap();
    let mut rows = Vec::new();
    for oid in 0..db.num_seqs() {
        let parts: Vec<String> = algorithms
            .iter()
            .map(|&algorithm| printed(&db.mask_ranges(oid, algorithm).unwrap()))
            .collect();
        if parts.iter().any(|part| part != "none") {
            rows.push(format!("{oid}|{}", parts.join("|")));
        }
    }
    rows.join(" ")
}

#[test]
fn masks_reach_the_sequences_that_ncbi_gives_them_to() {
    let mut fates: BTreeMap<&str, usize> = BTreeMap::new();
    let (mut matched, mut no_match) = (0, 0);
    for row in rows() {
        if row.loss.starts_with("env:") {
            *fates.entry("env").or_default() += 1;
            continue;
        }
        for version in [DbVersion::V5, DbVersion::V4] {
            let case = Case {
                name: row.name.clone(),
                protein: row.protein,
                inputs: vec![dir().join("fasta").join(&row.fasta)],
                mask_data: Some(row.files.clone()),
                mask_id: None,
                mask_desc: None,
                root: dir().join("case"),
                max_file_size: 3_000_000_000,
                title: row.fasta.clone(),
                version,
            };
            let tmp = TempDir::new("mask-ids");
            let indices = Indices {
                full: row.parse,
                hash: false,
            };
            let run = run_case_with(
                &case,
                &tmp.base("db"),
                "Oct 20, 2026  1:00 PM",
                &[],
                indices,
            );
            let stderr = run.stderr_text();
            let context = format!("{} ({version:?}): {stderr}", row.name);
            match row.loss.as_str() {
                "refuse" => {
                    // The library's flow (the CLI drops the destructor message after a refusal:
                    // `tests_cli`).
                    assert_eq!(run.status, 1, "{context}");
                    assert!(stderr.starts_with("losat makeblastdb: "), "{context}");
                    assert!(stderr.contains("memory addresses"), "{context}");
                    *fates.entry("refused").or_default() += 1;
                }
                "match" if row.expected == "NOMATCH" => {
                    assert_eq!(run.status, 0, "{context}");
                    assert_eq!(stderr, DESTRUCTOR_MESSAGE, "{context}");
                    *fates.entry("no match").or_default() += 1;
                    no_match += 1;
                }
                "match" => {
                    assert_eq!(run.status, 0, "{context}");
                    assert_eq!(stderr, "", "{context}");
                    assert_eq!(
                        masked_rows(&tmp.base("db"), row.protein, &row.algorithms),
                        row.expected,
                        "{context}"
                    );
                    *fates.entry("matched").or_default() += 1;
                    matched += 1;
                }
                other => panic!("{}: unknown loss {other}", row.name),
            }
        }
    }
    println!("{fates:?}");
    // 132 cases: 3 not run (environment), 6 refused, 51 without a match, the rest matched; each
    // run twice (format versions 5 and 4).
    assert_eq!(fates["env"], 3);
    assert_eq!(fates["refused"], 2 * 7);
    assert_eq!(no_match, 2 * 51);
    assert_eq!(matched, 2 * (134 - 3 - 7) - no_match);
}
