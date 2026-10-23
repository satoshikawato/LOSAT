//! EXPERIMENT (LOSAT_X_BXDBSTREAM / LOSAT_X_BXDBSTREAMSHADOW): a `-db` search of BLASTX reads the
//! database's sequences where it searches them, a small block at a time, instead of a copy of the
//! whole database made before the first query batch.
//!
//! Without the switch, `BlastxDatabase::read_stores` copies the residues of every sequence of the
//! database, a sentinel byte after each and the start of each (`StoredProteins`) into anonymous
//! memory before the first query batch, and every stage reads the copy until the run ends: 2.71 GiB
//! for UniRef50's `uniref50.03` volume (2.76 G residues, 16.6 M sequences), 205 MiB for Swiss-Prot
//! (v4-D12 s6 probe, workstation), where NCBI maps the volume files and holds 61 to 95 MiB of
//! anonymous memory.
//!
//! With the switch on:
//!
//! - `BlastxDatabase::x_check_stores` reads the database once in the blocks of `read_stores`, with
//!   its checks in its order (the empty database, `check_protein_sequence`, the 5,000,000-residue
//!   limit), and keeps nothing but the length of the longest sequence. Its error is the error of
//!   `read_stores`, which fails the same query batch with the same text.
//! - The preliminary search reads the sequences of each pass over the database in ascending order
//!   (`XPass`): a serial pass and a query chunk (chunks searched side by side at 2+ threads) read blocks of `block_bytes`
//!   (`SeqDb::next_block`) and drop each block before the next. A task of the subject flow
//!   (`x_subject_flow.rs`) reads its own subjects with one read per volume, and a wave of the pool
//!   (where the flow is not chosen: `LOSAT_X_BXSUBJFLOW=0`, at most `k` subjects, the 2^28 bound;
//!   a wave of 64 subjects per worker) reads its subjects before its fork
//!   (`XStreamDb::read_range`). The passes that run at the same time (the query batches
//!   searched side by side at 2+ threads, the query chunks searched side by side, the tasks of the flow) each read their own
//!   blocks: nothing is shared or waited for between them.
//! - The sink reads the lengths of the subjects in its order from the index (`XCursor`, a window of
//!   `LENGTH_WINDOW` sequences per read): the gapped sink evaluates every subject with its length,
//!   with or without HSPs. The `-ungapped` sink reads the residues of the subjects with HSPs, which
//!   it counts the identities of, from a pass of its own (only the blocks that hold them).
//! - The composition redo (its speculation and its serial loop), the traceback of
//!   `-comp_based_stats 0` and `ReportedSubjects::for_batch` read each subject they need with
//!   `SeqDb::read_sequence` (`Subjects::read_encoded`), after the preliminary search of the batch,
//!   as NCBI fetches the sequences of its hit lists.
//!
//! A read that fails during the search (a database changed or damaged after the check, an error of
//! the file system) gives the error that the check gives for the same read (`XReadError` holds
//! `search_read_error`'s), which the run reports as it reports the check's (`native::search_error`:
//! NCBI's engine error, or LOSAT's refusal).
//!
//! The stored bytes are the same in every read (the sentinel before and after the residues is the
//! file's byte, which the check requires to be 0, as the copy holds), so every stage computes the
//! same values on the same subjects in the same order. The anonymous memory of the sequences is the
//! blocks of the passes that run at the same time.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:328-350, quoted
//! above `SeqDb::next_block`: NCBI's threads read the database in chunks of consecutive ordinal ids.
//!
//! `LOSAT_X_BXDBSTREAM` is a value switch of the strict set: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set. Its value is the size of a pass's block in KiB, 2 to
//! 1,048,576, 256 when it is not set or holds anything else (`1` among them). It acts for every
//! task of BLASTX with `-db`.
//!
//! `LOSAT_X_BXDBSTREAMSHADOW=1` (shadow mode, with the switch on): `read_stores` makes the copy too,
//! and every block, every window of lengths and every read of a subject is compared with it (the
//! length and every residue with the sentinel after them); the check's error is compared with the
//! copy's. The first difference stops the run (panic) after printing the line that is otherwise
//! printed at exit: `[X_SHADOW] BXDBSTREAM blocks= sequences= residues= lengths= reads= identical`.

use std::borrow::Cow;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::OnceLock;

use anyhow::{anyhow, Result};

use crate::algorithm::common::database::{search_read_error, StoredProteins};
use crate::seqdb::{SeqDb, SequenceBlock};

/// The program that the reader's errors name (`database::PROGRAM`).
const PROGRAM: &str = "BLASTX";

/// Sequences whose lengths the sink's cursor reads with one read of the index (64 KiB of offsets).
const LENGTH_WINDOW: usize = 16_384;

// No NCBI counterpart: reads the LOSAT_X_BXDBSTREAM / LOSAT_X_BXDBSTREAMSHADOW switches once; it
// does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDBSTREAM): the bytes of a pass's block, or `None` when the switch is off.
pub(crate) fn x_block_bytes() -> Option<u64> {
    static BYTES: OnceLock<Option<u64>> = OnceLock::new();
    *BYTES.get_or_init(|| {
        crate::utils::xswitch::value("LOSAT_X_BXDBSTREAM", 2..=1 << 20, 256)
            .map(|kib| kib as u64 * 1024)
    })
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAMSHADOW): whether the reads are compared with the copy.
pub(crate) fn x_shadow() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        x_block_bytes().is_some() && std::env::var_os("LOSAT_X_BXDBSTREAMSHADOW").is_some()
    })
}

/// A read of the database that failed during the search: the error of the check for the same read
/// (`search_read_error`: NCBI's `PreliminaryReadError`, or LOSAT's refusal), which the run reports
/// as the check's (`native::search_error`).
#[derive(Debug)]
pub(crate) struct XReadError(pub(crate) anyhow::Error);

impl std::fmt::Display for XReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for XReadError {}

/// The reader's error `error` of a read during the search.
fn read_error(error: crate::seqdb::SeqDbError) -> anyhow::Error {
    XReadError(search_read_error(error, PROGRAM)).into()
}

/// The database of a search that reads its sequences in blocks.
pub(crate) struct XStreamDb<'a> {
    db: &'a SeqDb,
    count: usize,
    total_length: u64,
    /// The residues of the longest sequence (`StoredProteins::max_length`, from the check).
    longest: usize,
    block_bytes: u64,
    /// Shadow mode: the copy that every read is compared with.
    shadow: Option<&'a StoredProteins>,
}

impl<'a> XStreamDb<'a> {
    /// The database `db` after its check (`BlastxDatabase::x_check_stores`, which counted
    /// `longest`).
    pub(crate) fn new(
        db: &'a SeqDb,
        longest: usize,
        block_bytes: u64,
        shadow: Option<&'a StoredProteins>,
    ) -> Self {
        let stream = Self {
            db,
            count: db.num_seqs() as usize,
            total_length: db.total_length(),
            longest,
            block_bytes,
            shadow,
        };
        if let Some(stores) = shadow {
            X_SHADOW_OPENED.store(true, Ordering::Relaxed);
            x_shadow_assert(
                (stores.len(), stores.total_length(), stores.max_length())
                    == (stream.count, stream.total_length, stream.longest),
                || {
                    format!(
                        "the copy holds {} sequences, {} residues, {} at most; the check {}, {}, {}",
                        stores.len(),
                        stores.total_length(),
                        stores.max_length(),
                        stream.count,
                        stream.total_length,
                        stream.longest
                    )
                },
            );
        }
        stream
    }

    pub(crate) fn len(&self) -> usize {
        self.count
    }

    pub(crate) fn total_length(&self) -> u64 {
        self.total_length
    }

    pub(crate) fn longest(&self) -> usize {
        self.longest
    }

    /// The block of a pass that starts at `oid`: as many sequences as fit `block_bytes`
    /// (`SeqDb::next_block`, at least one).
    fn read_from(&self, oid: usize) -> Result<XBlock> {
        let range = u32::try_from(oid)
            .ok()
            .and_then(|first| self.db.next_block(first, self.block_bytes))
            .ok_or_else(|| anyhow!("BLASTX database ordinal id {oid} is past the database"))?;
        let part = self.db.read_block(range).map_err(read_error)?;
        self.checked(XBlock {
            first: oid,
            end: range.end_ordinal() as usize,
            parts: vec![part],
        })
    }

    /// The sequences `first..end` (non-empty), one read per volume.
    pub(crate) fn read_range(&self, first: usize, end: usize) -> Result<XBlock> {
        let (first32, end32) = (u32::try_from(first)?, u32::try_from(end)?);
        let mut parts = Vec::with_capacity(1);
        let mut at = first32;
        while at < end32 {
            let range = self.db.x_block_range(at, end32).ok_or_else(|| {
                anyhow!("BLASTX database ordinal ids {first}..{end} are past the database")
            })?;
            parts.push(self.db.read_block(range).map_err(read_error)?);
            at = range.end_ordinal();
        }
        self.checked(XBlock { first, end, parts })
    }

    /// The lengths of the sequences from `oid` on, at most `LENGTH_WINDOW` of them (to the end of
    /// `oid`'s volume), from the index.
    fn read_lengths(&self, oid: usize) -> Result<Vec<u32>> {
        let first = u32::try_from(oid)?;
        let end = first.saturating_add(LENGTH_WINDOW as u32);
        let lengths = self.db.x_protein_lengths(first, end).map_err(read_error)?;
        if let Some(stores) = self.shadow {
            for (i, &length) in lengths.iter().enumerate() {
                x_shadow_assert(stores.length(oid + i) == length as usize, || {
                    format!(
                        "sequence {}: the index says {length} residues, the copy {}",
                        oid + i,
                        stores.length(oid + i)
                    )
                });
            }
            X_SHADOW_LENGTHS.fetch_add(lengths.len() as u64, Ordering::Relaxed);
        }
        Ok(lengths)
    }

    /// The residues of `oid` (`SeqDb::read_sequence`), for the stages after the preliminary
    /// search.
    fn read_one(&self, oid: usize) -> Result<Vec<u8>> {
        let block = self
            .db
            .read_sequence(u32::try_from(oid)?)
            .map_err(read_error)?;
        let residues = block.protein(0).to_vec();
        if let Some(stores) = self.shadow {
            let copy = stores.residues_with_sentinel(oid);
            x_shadow_assert(residues[..] == copy[..copy.len() - 1], || {
                format!("sequence {oid}: the read differs from the copy")
            });
            X_SHADOW_READS.fetch_add(1, Ordering::Relaxed);
        }
        Ok(residues)
    }

    /// The residues of `oid` followed by the sentinel byte, read where no block holds them: no
    /// stage does (every preliminary read has its pass's or task's block; the sink and the later
    /// stages read through `XCursor` and `read_encoded`), so it is kept out of the callers' code.
    #[cold]
    #[inline(never)]
    fn outside_block(&self, oid: usize) -> Vec<u8> {
        debug_assert!(false, "subject {oid} is read outside a block");
        let mut residues = self.read_one(oid).expect("the database read");
        residues.push(0);
        residues
    }

    /// `block` after the shadow mode's comparison with the copy.
    fn checked(&self, block: XBlock) -> Result<XBlock> {
        if let Some(stores) = self.shadow {
            let mut residues = 0u64;
            for oid in block.first..block.end {
                let read = block.residues_with_sentinel(oid);
                x_shadow_assert(read == stores.residues_with_sentinel(oid), || {
                    format!("sequence {oid}: the block differs from the copy")
                });
                residues += read.len() as u64 - 1;
            }
            X_SHADOW_BLOCKS.fetch_add(1, Ordering::Relaxed);
            X_SHADOW_SEQUENCES.fetch_add((block.end - block.first) as u64, Ordering::Relaxed);
            X_SHADOW_RESIDUES.fetch_add(residues, Ordering::Relaxed);
        }
        Ok(block)
    }
}

/// Consecutive sequences `first..end` read from the database, in one block per volume.
pub(crate) struct XBlock {
    first: usize,
    end: usize,
    parts: Vec<SequenceBlock>,
}

impl XBlock {
    fn contains(&self, oid: usize) -> bool {
        (self.first..self.end).contains(&oid)
    }

    /// The residues of `oid` (in the block) followed by the sentinel byte.
    fn residues_with_sentinel(&self, oid: usize) -> &[u8] {
        let mut i = oid - self.first;
        for part in &self.parts {
            if i < part.len() {
                return &part.protein_with_sentinels(i)[1..];
            }
            i -= part.len();
        }
        unreachable!("sequence {oid} is outside its block")
    }
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAM): the subjects of a `-db` search read in blocks
/// (`subjects::Subjects::Stream`), with the block that the current pass or task has read, which
/// the preliminary stage reads its subjects from.
#[derive(Clone, Copy)]
pub(crate) struct XStream<'a> {
    pub(crate) db: &'a XStreamDb<'a>,
    pub(crate) block: Option<&'a XBlock>,
}

impl<'a> XStream<'a> {
    /// The residues of `oid` followed by the sentinel byte: from the block, where the preliminary
    /// stage reads them; a read of the database anywhere else.
    pub(crate) fn encoded_with_sentinel(&self, oid: usize) -> Cow<'a, [u8]> {
        match self.block {
            Some(block) if block.contains(oid) => Cow::Borrowed(block.residues_with_sentinel(oid)),
            _ => Cow::Owned(self.db.outside_block(oid)),
        }
    }

    /// The residues of `oid`, as `encoded_with_sentinel` without the sentinel byte.
    pub(crate) fn encoded(&self, oid: usize) -> Cow<'a, [u8]> {
        match self.encoded_with_sentinel(oid) {
            Cow::Borrowed(bytes) => Cow::Borrowed(&bytes[..bytes.len() - 1]),
            Cow::Owned(mut bytes) => {
                bytes.pop();
                Cow::Owned(bytes)
            }
        }
    }

    /// The length of `oid`: from the block, or the index anywhere else (the sink reads it from
    /// its `XCursor`).
    pub(crate) fn length(&self, oid: usize) -> usize {
        match self.block {
            Some(block) if block.contains(oid) => block.residues_with_sentinel(oid).len() - 1,
            _ => self.db.outside_block(oid).len() - 1,
        }
    }

    /// The residues of `oid` for the stages after the preliminary search: from the block, or one
    /// read of the database.
    pub(crate) fn read_encoded(&self, oid: usize) -> Result<Cow<'a, [u8]>> {
        match self.block {
            Some(block) if block.contains(oid) => {
                let bytes = block.residues_with_sentinel(oid);
                Ok(Cow::Borrowed(&bytes[..bytes.len() - 1]))
            }
            _ => Ok(Cow::Owned(self.db.read_one(oid)?)),
        }
    }
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAM): the blocks of one pass over the subjects in ascending order,
/// read as the pass reaches them; the block before is dropped first.
pub(crate) struct XPass<'a> {
    db: &'a XStreamDb<'a>,
    block: Option<XBlock>,
}

impl<'a> XPass<'a> {
    pub(crate) fn new(db: &'a XStreamDb<'a>) -> Self {
        Self { db, block: None }
    }

    /// The subjects with the block that holds `oid`, read if the pass's block does not.
    pub(crate) fn view(&mut self, oid: usize) -> Result<XStream<'_>> {
        if !self.block.as_ref().is_some_and(|block| block.contains(oid)) {
            self.block = None;
            self.block = Some(self.db.read_from(oid)?);
        }
        Ok(XStream {
            db: self.db,
            block: self.block.as_ref(),
        })
    }
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAM): the lengths that the sink of one query batch reads of its
/// subjects in its order, from windows of the index.
pub(crate) struct XCursor<'a> {
    db: &'a XStreamDb<'a>,
    first: usize,
    lengths: Vec<u32>,
}

impl<'a> XCursor<'a> {
    pub(crate) fn new(db: &'a XStreamDb<'a>) -> Self {
        Self {
            db,
            first: 0,
            lengths: Vec::new(),
        }
    }

    /// The length of subject `oid`.
    pub(crate) fn length(&mut self, oid: usize) -> Result<usize> {
        if !(self.first..self.first + self.lengths.len()).contains(&oid) {
            self.lengths = self.db.read_lengths(oid)?;
            self.first = oid;
        }
        Ok(self.lengths[oid - self.first] as usize)
    }
}

/// Whether a database was searched in shadow mode (the line is printed only then).
static X_SHADOW_OPENED: AtomicBool = AtomicBool::new(false);
static X_SHADOW_BLOCKS: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_SEQUENCES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_RESIDUES: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_LENGTHS: AtomicU64 = AtomicU64::new(0);
static X_SHADOW_READS: AtomicU64 = AtomicU64::new(0);

/// The shadow mode's line.
fn x_shadow_line() -> String {
    format!(
        "[X_SHADOW] BXDBSTREAM blocks={} sequences={} residues={} lengths={} reads={} identical",
        X_SHADOW_BLOCKS.load(Ordering::Relaxed),
        X_SHADOW_SEQUENCES.load(Ordering::Relaxed),
        X_SHADOW_RESIDUES.load(Ordering::Relaxed),
        X_SHADOW_LENGTHS.load(Ordering::Relaxed),
        X_SHADOW_READS.load(Ordering::Relaxed),
    )
}

/// Stops the run when `same` is false (shadow mode), after the shadow line.
pub(crate) fn x_shadow_assert(same: bool, what: impl FnOnce() -> String) {
    if !same {
        eprintln!("{}", x_shadow_line());
        panic!("LOSAT_X_BXDBSTREAMSHADOW: {}", what());
    }
}

/// EXPERIMENT (LOSAT_X_BXDBSTREAMSHADOW): the shadow mode's line at exit
/// (`x_seed_scan::print_shadow_stats`), after a search of a database.
pub(crate) fn x_print_shadow_summary() {
    if X_SHADOW_OPENED.load(Ordering::Relaxed) {
        eprintln!("{}", x_shadow_line());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastx::database::BlastxDatabase;
    use crate::blastinput::ncbi_environment::ApplicationSettings;

    fn open(name: &str) -> BlastxDatabase {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/seqdb/db")
            .join(name);
        BlastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    }

    // No NCBI counterpart: LOSAT's own reads of the stored bytes (`read_stores`'s copy is the
    // reference). LOSAT/tests/seqdb/expected/prot_mv.info.tsv: 10 sequences, 801 residues, in six
    // volumes; prot_v5 holds them in one.
    #[test]
    fn every_read_gives_the_bytes_and_lengths_of_the_copy() {
        for name in ["prot_mv", "prot_v5"] {
            let database = open(name);
            let stores = database.read_stores().unwrap();
            let longest = database.x_check_stores().unwrap();
            assert_eq!(longest, stores.max_length(), "{name}");
            let n = stores.len();
            let without_sentinel = |oid: usize| {
                let bytes = stores.residues_with_sentinel(oid);
                &bytes[..bytes.len() - 1]
            };
            // One byte: a block of one sequence; 200: a few; 1 MiB: the rest of the volume.
            for block_bytes in [1, 200, 1 << 20] {
                let stream =
                    XStreamDb::new(database.x_seqdb(), longest, block_bytes, Some(&stores));
                assert_eq!(
                    (stream.len(), stream.total_length()),
                    (n, stores.total_length())
                );
                let mut pass = XPass::new(&stream);
                for oid in 0..n {
                    let view = pass.view(oid).unwrap();
                    assert_eq!(
                        &*view.encoded_with_sentinel(oid),
                        stores.residues_with_sentinel(oid),
                        "{name} {block_bytes} {oid}"
                    );
                    assert_eq!(&*view.encoded(oid), without_sentinel(oid));
                    assert_eq!(view.length(oid), stores.length(oid));
                    assert_eq!(&*view.read_encoded(oid).unwrap(), without_sentinel(oid));
                }
                for first in 0..n {
                    for end in first + 1..=n {
                        let block = stream.read_range(first, end).unwrap();
                        for oid in first..end {
                            assert_eq!(
                                block.residues_with_sentinel(oid),
                                stores.residues_with_sentinel(oid),
                                "{name} {first}..{end} {oid}"
                            );
                        }
                    }
                }
                // The sink's order restarts at 0 for the next pass of a chunk.
                let mut cursor = XCursor::new(&stream);
                for oid in (0..n).chain(0..n) {
                    assert_eq!(cursor.length(oid).unwrap(), stores.length(oid), "{oid}");
                }
                let outside = XStream {
                    db: &stream,
                    block: None,
                };
                for oid in 0..n {
                    assert_eq!(&*outside.read_encoded(oid).unwrap(), without_sentinel(oid));
                }
            }
        }
    }

    #[test]
    fn a_read_that_fails_during_the_search_gives_the_error_of_the_check() {
        // A copy of prot_v5 whose sequence file is cut short after the check.
        let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db");
        let dir =
            std::env::temp_dir().join(format!("losat-bxdbstream-read-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for extension in ["pdb", "phr", "pin", "psq"] {
            let file = format!("prot_v5.{extension}");
            std::fs::copy(fixtures.join(&file), dir.join(&file)).unwrap();
        }
        let name = dir.join("prot_v5");
        let database =
            BlastxDatabase::open(name.to_str().unwrap(), &ApplicationSettings::default()).unwrap();
        let longest = database.x_check_stores().unwrap();
        let stream = XStreamDb::new(database.x_seqdb(), longest, 1 << 20, None);
        let psq = dir.join("prot_v5.psq");
        let bytes = std::fs::read(&psq).unwrap();
        std::fs::write(&psq, &bytes[..bytes.len() / 2]).unwrap();
        let error = XPass::new(&stream).view(0).map(|_| ()).unwrap_err();
        let XReadError(kept) = error.downcast::<XReadError>().unwrap();
        assert_eq!(
            kept.to_string(),
            database.x_check_stores().unwrap_err().to_string()
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn the_check_fails_as_the_copy_does() {
        let empty = open("empty");
        assert_eq!(
            empty.x_check_stores().unwrap_err().to_string(),
            empty.read_stores().unwrap_err().to_string()
        );
    }
}
