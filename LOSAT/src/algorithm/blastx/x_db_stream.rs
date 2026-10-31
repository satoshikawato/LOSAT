//! A `-db` search of BLASTX reads the database's sequences where it searches them, a small block
//! at a time (graduated 2026-10-27; the switch was `LOSAT_X_BXDBSTREAM`), instead of a copy of the
//! whole database made before the first query batch.
//!
//! The copy (`read_stores`, deleted) held the residues of every sequence of the database, a
//! sentinel byte after each and the start of each (`StoredProteins`) in anonymous memory before the
//! first query batch, and every stage read it until the run ended: 2.71 GiB for UniRef50's
//! `uniref50.03` volume (2.76 G residues, 16.6 M sequences), 205 MiB for Swiss-Prot (v4-D12 s6
//! probe, workstation), where NCBI maps the volume files and holds 61 to 95 MiB of anonymous memory.
//!
//! The stream:
//!
//! - `BlastxDatabase::x_check_stores` reads the database once in blocks of
//!   `PROTEIN_STORE_BLOCK_BYTES`, with the checks of the copy in its order (the empty database,
//!   `check_protein_sequence`, the 5,000,000-residue limit), and keeps nothing but the length of the
//!   longest sequence. Its error is the error that the copy gave, which fails the same query batch
//!   with the same text.
//! - The preliminary search reads the sequences of each pass over the database in ascending order
//!   (`XPass`): a serial pass and a query chunk (chunks searched side by side at 2+ threads) read
//!   blocks of `block_bytes` (`SeqDb::next_block`) and drop each block before the next. A task of
//!   the subject flow (`x_subject_flow.rs`) reads its own subjects with one read per volume, and a
//!   wave of the pool (where the flow is not chosen: `LOSAT_X_BXSUBJFLOW=0`, at most `k` subjects,
//!   the 2^28 bound; a wave of 64 subjects per worker) reads its subjects before its fork
//!   (`XStreamDb::read_range`). The passes that run at the same time (the query batches searched
//!   side by side at 2+ threads, the query chunks searched side by side, the tasks of the flow)
//!   each read their own blocks: nothing is shared or waited for between them.
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
//! file's byte, which the check requires to be 0, as the copy held), so every stage computes the
//! same values on the same subjects in the same order. The anonymous memory of the sequences is the
//! blocks of the passes that run at the same time.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:328-350, quoted
//! above `SeqDb::next_block`: NCBI's threads read the database in chunks of consecutive ordinal ids.
//!
//! It acts for every task of BLASTX with `-db`. A pass reads blocks of `BLOCK_BYTES` (256 KiB; the
//! variable `LOSAT_X_BXDBSTREAM=KiB` that set it, and the shadow mode `LOSAT_X_BXDBSTREAMSHADOW=1`
//! that compared every read with the copy, went with the copy).

use std::borrow::Cow;

use anyhow::{anyhow, Result};

use crate::algorithm::common::database::search_read_error;
use crate::seqdb::{SeqDb, SequenceBlock};

/// The program that the reader's errors name (`database::PROGRAM`).
const PROGRAM: &str = "BLASTX";

/// Sequences whose lengths the sink's cursor reads with one read of the index (64 KiB of offsets).
const LENGTH_WINDOW: usize = 16_384;

/// The bytes of a pass's block (`SeqDb::next_block`). No NCBI counterpart: NCBI maps the volume
/// files. The size changes nothing but the number of reads.
const BLOCK_BYTES: u64 = 256 * 1024;

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
    /// The residues of the longest sequence (from the check).
    longest: usize,
    block_bytes: u64,
    /// EXPERIMENT (LOSAT_X_BXTASKRES): the subject flow's task bounds by tasks' `k`, made once
    /// (`x_task_starts`).
    x_task_starts: std::sync::Mutex<Vec<(usize, std::sync::Arc<[u32]>)>>,
}

impl<'a> XStreamDb<'a> {
    /// The database `db` after its check (`BlastxDatabase::x_check_stores`, which counted
    /// `longest`).
    pub(crate) fn new(db: &'a SeqDb, longest: usize) -> Self {
        Self::with_block_bytes(db, longest, BLOCK_BYTES)
    }

    /// As `new`, with blocks of `block_bytes` (the unit tests read blocks of a few sequences).
    fn with_block_bytes(db: &'a SeqDb, longest: usize, block_bytes: u64) -> Self {
        Self {
            db,
            count: db.num_seqs() as usize,
            total_length: db.statistics_numbers().0,
            longest,
            block_bytes,
            x_task_starts: std::sync::Mutex::new(Vec::new()),
        }
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
        Ok(XBlock {
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
        Ok(XBlock { first, end, parts })
    }

    /// The lengths of the sequences from `oid` on, at most `LENGTH_WINDOW` of them (to the end of
    /// `oid`'s volume), from the index.
    fn read_lengths(&self, oid: usize) -> Result<Vec<u32>> {
        let first = u32::try_from(oid)?;
        let end = first.saturating_add(LENGTH_WINDOW as u32);
        self.db.x_protein_lengths(first, end).map_err(read_error)
    }

    /// EXPERIMENT (LOSAT_X_BXTASKRES): `read_lengths`, for the bounds of the subject flow's tasks
    /// (`x_subject_flow::XTasks`).
    #[cfg(feature = "parallel")]
    pub(crate) fn x_lengths_from(&self, oid: usize) -> Result<Vec<u32>> {
        self.read_lengths(oid)
    }

    /// EXPERIMENT (LOSAT_X_BXTASKRES): the subject flow's task bounds for tasks of `k` subjects:
    /// those made before in this run, or `make`'s (kept unless it gives `None`). A flow that asks
    /// while another makes them waits for them.
    #[cfg(feature = "parallel")]
    pub(crate) fn x_task_starts(
        &self,
        k: usize,
        make: impl FnOnce() -> Option<std::sync::Arc<[u32]>>,
    ) -> Option<std::sync::Arc<[u32]>> {
        let mut held = self
            .x_task_starts
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some((_, starts)) = held.iter().find(|(held_k, _)| *held_k == k) {
            return Some(starts.clone());
        }
        let starts = make()?;
        held.push((k, starts.clone()));
        Some(starts)
    }

    /// The residues of `oid` (`SeqDb::read_sequence`), for the stages after the preliminary
    /// search.
    fn read_one(&self, oid: usize) -> Result<Vec<u8>> {
        let block = self
            .db
            .read_sequence(u32::try_from(oid)?)
            .map_err(read_error)?;
        Ok(block.protein(0).to_vec())
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

    /// The residues of `oid` (in the block) between the sentinel bytes of the volume.
    fn x_residues_with_sentinels(&self, oid: usize) -> &[u8] {
        let mut i = oid - self.first;
        for part in &self.parts {
            if i < part.len() {
                return part.protein_with_sentinels(i);
            }
            i -= part.len();
        }
        unreachable!("sequence {oid} is outside its block")
    }
}

/// The subjects of a `-db` search read in blocks (`subjects::Subjects::Stream`), with the block
/// that the current pass or task has read, which the preliminary stage reads its subjects from.
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

    /// EXPERIMENT (LOSAT_X_BXUNGAPSENT): `encoded_with_sentinel` with the sentinel byte before the
    /// residues too.
    pub(crate) fn x_encoded_with_sentinels(&self, oid: usize) -> Cow<'a, [u8]> {
        match self.block {
            Some(block) if block.contains(oid) => {
                Cow::Borrowed(block.x_residues_with_sentinels(oid))
            }
            _ => {
                let mut residues = vec![0];
                residues.extend_from_slice(&self.db.outside_block(oid));
                Cow::Owned(residues)
            }
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

/// The blocks of one pass over the subjects in ascending order, read as the pass reaches them;
/// the block before is dropped first.
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

/// The lengths that the sink of one query batch reads of its subjects in its order, from windows
/// of the index.
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

    // No NCBI counterpart: LOSAT's own reads of the stored bytes (`SeqDb::read_sequence`, one
    // sequence at a time, is the reference, as the copy of the whole database was before it was
    // deleted). LOSAT/tests/seqdb/expected/prot_mv.info.tsv: 10 sequences, 801 residues, in six
    // volumes; prot_v5 holds them in one.
    #[test]
    fn every_read_gives_the_bytes_and_lengths_of_the_files() {
        for name in ["prot_mv", "prot_v5"] {
            let database = open(name);
            let db = database.x_seqdb();
            let n = db.num_seqs() as usize;
            // The residues of every sequence, one read each.
            let residues: Vec<Vec<u8>> = (0..n)
                .map(|oid| db.read_sequence(oid as u32).unwrap().protein(0).to_vec())
                .collect();
            let total: usize = residues.iter().map(Vec::len).sum();
            if name == "prot_mv" {
                assert_eq!((n, total), (10, 801));
            }
            let longest = database.x_check_stores().unwrap();
            assert_eq!(
                longest,
                residues.iter().map(Vec::len).max().unwrap(),
                "{name}"
            );
            let with_sentinel = |oid: usize| [&residues[oid][..], &[0u8][..]].concat();
            // The same bytes with the sentinel before them.
            let with_sentinels = |oid: usize| [&[0u8][..], &with_sentinel(oid)[..]].concat();
            // One byte: a block of one sequence; 200: a few; 1 MiB: the rest of the volume.
            for block_bytes in [1, 200, 1 << 20] {
                let stream = XStreamDb::with_block_bytes(db, longest, block_bytes);
                assert_eq!(
                    (stream.len(), stream.total_length()),
                    (n, db.total_length())
                );
                let mut pass = XPass::new(&stream);
                for oid in 0..n {
                    let view = pass.view(oid).unwrap();
                    assert_eq!(
                        &*view.encoded_with_sentinel(oid),
                        &with_sentinel(oid)[..],
                        "{name} {block_bytes} {oid}"
                    );
                    assert_eq!(
                        &*view.x_encoded_with_sentinels(oid),
                        &with_sentinels(oid)[..],
                        "{name} {block_bytes} {oid}"
                    );
                    assert_eq!(&*view.encoded(oid), &residues[oid][..]);
                    assert_eq!(view.length(oid), residues[oid].len());
                    assert_eq!(&*view.read_encoded(oid).unwrap(), &residues[oid][..]);
                }
                for first in 0..n {
                    for end in first + 1..=n {
                        let block = stream.read_range(first, end).unwrap();
                        for oid in first..end {
                            assert_eq!(
                                block.residues_with_sentinel(oid),
                                &with_sentinel(oid)[..],
                                "{name} {first}..{end} {oid}"
                            );
                            assert_eq!(
                                block.x_residues_with_sentinels(oid),
                                &with_sentinels(oid)[..],
                                "{name} {first}..{end} {oid}"
                            );
                        }
                    }
                }
                // The sink's order restarts at 0 for the next pass of a chunk.
                let mut cursor = XCursor::new(&stream);
                for oid in (0..n).chain(0..n) {
                    assert_eq!(cursor.length(oid).unwrap(), residues[oid].len(), "{oid}");
                }
                let outside = XStream {
                    db: &stream,
                    block: None,
                };
                for oid in 0..n {
                    assert_eq!(&*outside.read_encoded(oid).unwrap(), &residues[oid][..]);
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
        let stream = XStreamDb::new(database.x_seqdb(), longest);
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

    // The error of the first batch with a valid context for a database without sequences is NCBI's
    // (`blastx_db_search.rs`): the reader's, not LOSAT's.
    #[test]
    fn the_check_of_an_empty_database_fails_as_ncbis_preliminary_search() {
        let error = open("empty").x_check_stores().unwrap_err();
        let error = error
            .downcast::<crate::algorithm::common::database::PreliminaryReadError>()
            .expect("NCBI's error of the preliminary search");
        assert_eq!(error.0, "Database memory map file error");
    }
}
