//! Damaged copies of the database files that [`SeqDb`] opens and reads: the index (`.pin`,
//! `.nin`), sequence (`.psq`, `.nsq`) and header (`.phr`, `.nhr`) files of versions 4 and 5, two
//! alias files (`.pal`, `.nal`) and a volume index that one of them names, and the index (`.?aa`)
//! and little-endian data (`.?ac`) files of the mask columns. Each target file of a fixture under
//! `LOSAT/tests/seqdb/db` is cut at many lengths, then damaged in rounds drawn from a fixed seed:
//! bits flipped, a 4-byte word (a length, a count or an offset) rewritten, or some of its first
//! 64 bytes written. After each damage the database is opened, and its sequences, headers,
//! blocks and masks are read through the calls that the searches and `blastdbcmd` make
//! (`exercise`). Any result but a panic or a hang is accepted: an error or, where the damage
//! leaves a readable file, a value. The loops stop at `MAX_ORDINALS` sequences and `MAX_BLOCKS`
//! ranges, so a damaged count cannot keep a case going. A panic is caught and named with its
//! file and case, and the test fails listing them all. Each target must also have seeded rounds
//! that opened and seeded rounds whose outcome (opened or not, the number of errors) differs from
//! the intact copy's, so that the damage reaches the reader.
//!
//! The LMDB files (`.pdb`, `.ndb`) are damaged by `corrupt_files_give_errors_not_panics`
//! (`lmdb_read.rs`), and the index files on the `LOSAT_X_DBIDXREAD=0` path by
//! `open_gives_the_errors_of_parse` (`x_dbidxread_tests.rs`). The form follows those tests and
//! `damaged_files_never_panic` (`tax_sqlite.rs`); the generator is the xorshift of
//! `isam_tests.rs`.

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};

use super::alias::SearchPath;
use super::volume::SequenceBlock;
use super::{SeqDb, SeqDbError, SeqType};

/// Seeded rounds of damage per target file, after the cuts.
const ROUNDS: usize = 1000;
/// The most sequences read of a damaged database. The fixtures hold at most 30 (`alias_nested`).
const MAX_ORDINALS: u32 = 64;
/// The most block or header ranges read of a damaged database.
const MAX_BLOCKS: usize = 64;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db")
}

/// Marsaglia's xorshift64, as in `isam_tests.rs`.
struct Rng(u64);

impl Rng {
    /// A generator seeded from the bytes of `name`: the v4 and v5 fixtures have files of the
    /// same bytes and the same length of name, which then get different damage.
    fn seeded(name: &str) -> Self {
        let seed = name.bytes().fold(0x5eed_0000_0000_0001u64, |seed, byte| {
            seed.rotate_left(7) ^ u64::from(byte)
        });
        Self(seed | 1)
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n.max(1) as u64) as usize
    }
}

/// A directory of its own for each test (memory `parallel-test-shared-file-race`), removed at
/// the end.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-damaged-{tag}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    /// Copies every fixture file whose name starts with one of `stems` followed by a dot.
    fn copy_fixtures(&self, stems: &[&str]) {
        let mut copied = 0;
        for entry in std::fs::read_dir(fixture_dir()).unwrap() {
            let entry = entry.unwrap();
            let name = entry.file_name().into_string().unwrap();
            if stems.iter().any(|stem| {
                name.strip_prefix(stem)
                    .is_some_and(|rest| rest.starts_with('.'))
            }) {
                std::fs::copy(entry.path(), self.0.join(&name)).unwrap();
                copied += 1;
            }
        }
        assert!(copied > 0, "no fixture files for {stems:?}");
    }

    fn name(&self, stem: &str) -> String {
        self.0.join(stem).to_str().unwrap().to_string()
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// How a seeded round damages a file.
#[derive(Clone, Copy, Debug)]
enum Damage {
    /// The file cut to a length drawn below its length.
    Cut,
    /// One to four bits flipped anywhere.
    Bits,
    /// A 4-byte word at a multiple of 4 set, in either byte order, to a value that breaks a
    /// length, count or offset: 0, 1, 2^31 - 1, 2^31, 2^32 - 1, the file length minus one, the
    /// file length, plus one or doubled, or the old value plus or minus one.
    Word,
    /// One to four of the first 64 bytes set to drawn values.
    Head,
}

fn damage(orig: &[u8], rng: &mut Rng, round: usize) -> (Damage, Vec<u8>) {
    let mut bytes = orig.to_vec();
    let kind = [Damage::Cut, Damage::Bits, Damage::Word, Damage::Head][round % 4];
    if bytes.is_empty() {
        return (kind, bytes);
    }
    match kind {
        Damage::Cut => bytes.truncate(rng.below(orig.len())),
        Damage::Bits => {
            for _ in 0..1 + rng.below(4) {
                let at = rng.below(bytes.len());
                bytes[at] ^= 1 << rng.below(8);
            }
        }
        Damage::Word => {
            if bytes.len() >= 4 {
                let at = rng.below(bytes.len() / 4) * 4;
                let big_endian = rng.next().is_multiple_of(2);
                let word: [u8; 4] = bytes[at..at + 4].try_into().unwrap();
                let old = if big_endian {
                    u32::from_be_bytes(word)
                } else {
                    u32::from_le_bytes(word)
                };
                let len = orig.len() as u32;
                let values = [
                    0,
                    1,
                    i32::MAX as u32,
                    1 << 31,
                    u32::MAX,
                    len - 1,
                    len,
                    len + 1,
                    len.wrapping_mul(2),
                    old.wrapping_add(1),
                    old.wrapping_sub(1),
                ];
                let value = values[rng.below(values.len())];
                let new = if big_endian {
                    value.to_be_bytes()
                } else {
                    value.to_le_bytes()
                };
                bytes[at..at + 4].copy_from_slice(&new);
            }
        }
        Damage::Head => {
            let head = bytes.len().min(64);
            for _ in 0..1 + rng.below(4) {
                let at = rng.below(head);
                bytes[at] = rng.next() as u8;
            }
        }
    }
    (kind, bytes)
}

/// What the cases of one target file gave.
#[derive(Default)]
struct Tally {
    opened: usize,
    errors: usize,
}

impl Tally {
    /// The value of `result`, or `None` after counting its error. The text and NCBI's form of
    /// the error are made for every error, damaged file or not.
    fn check<T>(&mut self, result: Result<T, SeqDbError>) -> Option<T> {
        match result {
            Ok(value) => Some(value),
            Err(error) => {
                let _ = error.to_string();
                let _ = error.ncbi_message();
                let _ = error.ncbi_search_message();
                self.errors += 1;
                None
            }
        }
    }
}

/// Opens the database `name` (in the directory `dir`) and reads it through the calls of the
/// searches and of `blastdbcmd`.
fn exercise(dir: &Path, name: &str, seqtype: SeqType, masks: bool, tally: &mut Tally) {
    // `blastdbcmd -info` and the metadata read the alias files again, without the volumes.
    let path = SearchPath::new(dir, None, None).unwrap();
    tally.check(SeqDb::info_volumes(name, seqtype, &path));
    tally.check(SeqDb::volume_and_alias_paths(name, seqtype, &path));
    let opened = if masks {
        SeqDb::open_with_masks(name, seqtype)
    } else {
        SeqDb::open(name, seqtype)
    };
    let Some(db) = tally.check(opened) else {
        return;
    };
    tally.opened += 1;
    let _ = (db.title(), db.date(), db.total_length(), db.max_length());
    let _ = (db.min_length(), db.statistics_numbers(), db.is_version5());
    tally.check(db.has_parsed_seqids());
    let ordinals = db.num_seqs().min(MAX_ORDINALS);
    for ordinal in 0..ordinals {
        if let Some(block) = tally.check(db.read_sequence(ordinal)) {
            read_letters(&block, seqtype, tally);
        }
        match seqtype {
            SeqType::Protein => {
                tally.check(db.protein_length(ordinal));
            }
            SeqType::Nucleotide => {
                if let Some(length) = tally.check(db.nucleotide_length(ordinal)) {
                    let ranges = [
                        (0, length.min(40)),
                        (length / 2, length),
                        (length, length + 1),
                    ];
                    for (from, to) in ranges {
                        tally.check(db.nucleotide_iupac_range(ordinal, from, to));
                    }
                    if let Some(ambiguity) = tally.check(db.nucleotide_ambiguity(ordinal)) {
                        tally.check(db.nucleotide_iupac_range_with(ordinal, 0, length, &ambiguity));
                    }
                }
            }
        }
        tally.check(db.deflines(ordinal));
        tally.check(db.read_header(ordinal));
        tally.check(db.defline_set(ordinal));
        tally.check(db.defline_taxids(ordinal));
    }
    // The blocks of the sequence file, of 1 KiB each: one at a time as `next_block` gives them
    // and as `x_next_blocks` lists them, read by each of the block readers.
    let mut ranges = Vec::new();
    let mut first = 0;
    while ranges.len() < MAX_BLOCKS {
        let Some(range) = db.next_block(first, 1024) else {
            break;
        };
        ranges.push(range);
        let next = range.first_ordinal.saturating_add(range.len().max(1));
        if next <= first {
            break;
        }
        first = next;
    }
    ranges.extend(db.x_next_blocks(0, 1024, MAX_BLOCKS));
    ranges.extend(db.x_block_range(0, ordinals));
    for range in ranges {
        if let Some(block) = tally.check(db.read_block(range)) {
            read_letters(&block, seqtype, tally);
        }
        if let Some(block) = tally.check(db.x_read_block_mapped(range)) {
            read_letters(&block, seqtype, tally);
        }
        tally.check(db.x_read_block_into(range, Vec::new()));
        let end = range.first_ordinal.saturating_add(range.len());
        match seqtype {
            SeqType::Protein => {
                tally.check(db.x_protein_lengths(range.first_ordinal, end));
            }
            // `x_nucleotide_check` takes one sequence.
            SeqType::Nucleotide if range.len() == 1 => {
                tally.check(db.x_nucleotide_check(range));
            }
            SeqType::Nucleotide => {}
        }
    }
    // The header ranges that the searches visit before they report.
    for max_seqs in [2, 1024] {
        for range in db.header_ranges(max_seqs).into_iter().take(MAX_BLOCKS) {
            db.visit_deflines(range, |_, result| {
                tally.check(result);
                true
            });
            db.visit_deflines_taxids(range, |_, result| {
                tally.check(result);
                true
            });
            db.visit_defline_taxids(range, |_, result| {
                tally.check(result);
                true
            });
            db.visit_defline_set_taxids(range, |_, result| {
                tally.check(result);
                true
            });
            db.visit_sequence_headers(range, |_, result| {
                tally.check(result);
                true
            });
        }
    }
    if masks {
        if let Some(algorithms) = tally.check(db.mask_algorithms()) {
            tally.check(algorithms.descriptions());
            for id in algorithms.ids().into_iter().take(8) {
                tally.check(algorithms.details(id));
                for ordinal in 0..ordinals {
                    tally.check(db.mask_ranges(ordinal, id));
                }
            }
        }
    }
}

/// Decodes every sequence of a block in each form the searches take.
fn read_letters(block: &SequenceBlock, seqtype: SeqType, tally: &mut Tally) {
    for i in 0..block.len() {
        match seqtype {
            SeqType::Protein => {
                let _ = (
                    block.protein(i).len(),
                    block.protein_with_sentinels(i).len(),
                );
            }
            SeqType::Nucleotide => {
                let _ = (
                    block.nucleotide_packed(i).len(),
                    block.ambiguity_bytes(i).len(),
                );
                tally.check(block.nucleotide_len(i));
                tally.check(block.ambiguity_runs(i));
                tally.check(block.nucleotide_blastna(i));
                tally.check(block.x_nucleotide_blastna_mapped(i));
                tally.check(block.nucleotide_iupac(i));
                tally.check(block.nucleotide_stored(i));
            }
        }
    }
}

/// Damages each file of `targets` (in turn, the others intact) of the databases `stems`, and
/// reads the database `open` after each damage. Fails listing every case that panicked and
/// every target whose damage did not reach the reader.
fn sweep(tag: &str, stems: &[&str], open: &str, seqtype: SeqType, masks: bool, targets: &[&str]) {
    let dir = TempDir::new(tag);
    dir.copy_fixtures(stems);
    let name = dir.name(open);
    let mut intact = Tally::default();
    exercise(&dir.0, &name, seqtype, masks, &mut intact);
    assert_eq!(intact.opened, 1, "{tag}: the intact copy must open");
    let mut problems = Vec::new();
    for target in targets {
        let path = dir.0.join(target);
        let orig = std::fs::read(&path).unwrap();
        let mut rng = Rng::seeded(&format!("{tag}/{target}"));
        // Every cut length of the small files, a spread of them for the larger ones, then the
        // seeded rounds.
        let step = orig.len().div_ceil(256).max(1);
        let mut cases: Vec<(Option<(usize, Damage)>, Vec<u8>)> = (0..orig.len())
            .step_by(step)
            .map(|cut| (None, orig[..cut].to_vec()))
            .collect();
        for round in 0..ROUNDS {
            let (kind, bytes) = damage(&orig, &mut rng, round);
            cases.push((Some((round, kind)), bytes));
        }
        let mut tally = Tally::default();
        let (mut seeded_opened, mut seeded_changed) = (0, 0);
        for (case, bytes) in cases {
            std::fs::write(&path, &bytes).unwrap();
            let (opened, errors) = (tally.opened, tally.errors);
            let outcome = catch_unwind(AssertUnwindSafe(|| {
                exercise(&dir.0, &name, seqtype, masks, &mut tally)
            }));
            if case.is_some() {
                seeded_opened += usize::from(tally.opened > opened);
                if tally.opened == opened || tally.errors - errors != intact.errors {
                    seeded_changed += 1;
                }
            }
            if outcome.is_err() {
                let what = match case {
                    None => format!("cut to {} bytes", bytes.len()),
                    Some((round, kind)) => format!("round {round} ({kind:?})"),
                };
                problems.push(format!("{target}: {what} panicked"));
            }
        }
        std::fs::write(&path, &orig).unwrap();
        if seeded_opened == 0 || seeded_changed == 0 {
            problems.push(format!(
                "{target}: the seeded damage did not reach the reader (opened {seeded_opened}, \
                 changed {seeded_changed} of {ROUNDS})"
            ));
        }
    }
    assert!(
        problems.is_empty(),
        "{tag}: {} problems:\n{}",
        problems.len(),
        problems.join("\n")
    );
}

#[test]
fn damaged_protein_volume_files_never_panic() {
    sweep(
        "p4",
        &["prot_v4"],
        "prot_v4",
        SeqType::Protein,
        false,
        &["prot_v4.pin", "prot_v4.psq", "prot_v4.phr"],
    );
    sweep(
        "p5",
        &["prot_v5"],
        "prot_v5",
        SeqType::Protein,
        false,
        &["prot_v5.pin", "prot_v5.psq", "prot_v5.phr"],
    );
}

#[test]
fn damaged_nucleotide_volume_files_never_panic() {
    sweep(
        "n4",
        &["nucl_v4"],
        "nucl_v4",
        SeqType::Nucleotide,
        false,
        &["nucl_v4.nin", "nucl_v4.nsq", "nucl_v4.nhr"],
    );
    sweep(
        "n5",
        &["nucl_v5"],
        "nucl_v5",
        SeqType::Nucleotide,
        false,
        &["nucl_v5.nin", "nucl_v5.nsq", "nucl_v5.nhr"],
    );
}

#[test]
fn damaged_alias_files_never_panic() {
    sweep(
        "nal",
        &["nucl_mv"],
        "nucl_mv",
        SeqType::Nucleotide,
        false,
        &["nucl_mv.nal", "nucl_mv.03.nin"],
    );
    sweep(
        "pal",
        &["alias_nested", "prot_v4", "prot_mv", "prot_old"],
        "alias_nested",
        SeqType::Protein,
        false,
        &["alias_nested.pal", "prot_mv.pal"],
    );
}

/// The index (`.?aa`) and the little-endian data (`.?ac`) of the mask column. The big-endian
/// data file (`.?ab`) is not read when the little-endian one exists (`column::open_all`).
#[test]
fn damaged_column_files_of_masks_never_panic() {
    sweep(
        "nmask",
        &["mask_nucl_v4"],
        "mask_nucl_v4",
        SeqType::Nucleotide,
        true,
        &["mask_nucl_v4.naa", "mask_nucl_v4.nac"],
    );
    sweep(
        "pmask",
        &["mask_prot_v4"],
        "mask_prot_v4",
        SeqType::Protein,
        true,
        &["mask_prot_v4.paa", "mask_prot_v4.pac"],
    );
}
