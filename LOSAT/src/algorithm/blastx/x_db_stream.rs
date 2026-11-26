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
//! - Before the first query batch the run reads no sequence (`BlastxDatabase::x_open_stream`): the
//!   length of the longest sequence comes from the index's offsets, and the empty database, an
//!   index that cannot be read and a longest above 5,000,000 residues give the error of the check (`x_check_stores`, which
//!   reads the database once in blocks of `PROTEIN_STORE_BLOCK_BYTES` with the checks of the copy
//!   in its order: the empty database, `check_protein_sequence`, the 5,000,000-residue limit), which
//!   fails the first query batch with a valid context. The reads of the search check each block
//!   they read the same way (`XStreamDb::x_check`) until they have checked every sequence once.
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
//! A read or a check that fails during the search (a damaged database, an error of the file system)
//! gives the error that the check gives for the same read (`XReadError` holds `search_read_error`'s
//! or the check's). Until the reads have checked every sequence, the run reports the check's first
//! error in ordinal order, whichever sequence the threads read first (`native::search_error`:
//! NCBI's engine error, or LOSAT's refusal), so the text and the query batch it fails are those of
//! the check that read the database before the first query batch. An error that the batch raises
//! before its first read but after the point where a database refused before the search fails it
//! (its setup's refusals: a search space or a diagonal length past NCBI's integers, a query split
//! past its chunks) comes after the database's error, as before (`XPastHead`). After that, a read
//! that fails (a file changed during the run) gives its own error.
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

use crate::algorithm::common::database::{check_protein_sequence, search_read_error};
use crate::seqdb::{SeqDb, SequenceBlock};

use super::database::{too_long, x_check_whole, MAX_DBSEQ_LEN};

/// The program that the reader's errors name (`database::PROGRAM`).
const PROGRAM: &str = "BLASTX";

/// Sequences whose lengths the sink's cursor reads with one read of the index (64 KiB of offsets).
const LENGTH_WINDOW: usize = 16_384;

/// The bytes of a pass's block (`SeqDb::next_block`). No NCBI counterpart: NCBI maps the volume
/// files. The size changes nothing but the number of reads.
const BLOCK_BYTES: u64 = 256 * 1024;

/// A read or a check of the database that failed during the search: the error of the check for
/// the same read (`search_read_error`: NCBI's `PreliminaryReadError`, or LOSAT's refusal of the
/// sequence, `x_check`), which the run reports as the check's (`native::search_error`).
#[derive(Debug)]
pub(crate) struct XReadError(pub(crate) anyhow::Error);

impl std::fmt::Display for XReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for XReadError {}

/// An error of a search raised after the point where a database refused before the search fails
/// it, while the reads have not checked every sequence (`Subjects::x_past_head`): the run reports
/// the database's first error instead, if the check finds one (`native::search_error`).
#[derive(Debug)]
pub(crate) struct XPastHead(pub(crate) anyhow::Error);

impl std::fmt::Display for XPastHead {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl std::error::Error for XPastHead {}

/// The reader's error `error` of a read during the search.
fn read_error(error: crate::seqdb::SeqDbError) -> anyhow::Error {
    XReadError(search_read_error(error, PROGRAM)).into()
}

/// The database of a search that reads its sequences in blocks.
pub(crate) struct XStreamDb<'a> {
    db: &'a SeqDb,
    count: usize,
    /// The residues of the longest sequence (`BlastxDatabase::x_open_stream`).
    longest: usize,
    block_bytes: u64,
    /// Whether the reads have checked every sequence (`x_check`).
    x_checked: std::sync::atomic::AtomicBool,
    /// The ordinal ranges the reads have checked, merged, until they cover the database.
    x_covered: std::sync::Mutex<Vec<(usize, usize)>>,
    /// Whether the whole check finds an error (`x_damaged`), run at most once for the stream.
    x_damaged: std::sync::OnceLock<bool>,
    /// The whole check's error, kept by `x_damaged` for the run's report (`x_first_error`).
    x_first: std::sync::Mutex<Option<anyhow::Error>>,
    /// EXPERIMENT (LOSAT_X_BXTASKRES): the subject flow's task bounds by tasks' `k`, made once
    /// (`x_task_starts`).
    #[cfg(feature = "parallel")]
    x_task_starts: std::sync::Mutex<Vec<(usize, std::sync::Arc<[u32]>)>>,
}

impl<'a> XStreamDb<'a> {
    /// The database `db` with `(longest, count)` computed by `BlastxDatabase::x_open_stream`.
    pub(crate) fn with_scan_count(db: &'a SeqDb, longest: usize, count: usize) -> Self {
        Self::with_block_bytes_and_count(db, longest, count, BLOCK_BYTES)
    }

    /// The database `db` with the longest of `BlastxDatabase::x_open_stream` (or of the check,
    /// `x_check_stores`), scanning all `db.num_seqs()` ordinals (for unit tests).
    #[cfg(test)]
    pub(crate) fn new(db: &'a SeqDb, longest: usize) -> Self {
        Self::with_block_bytes(db, longest, BLOCK_BYTES)
    }

    /// As `new`, with blocks of `block_bytes` (the unit tests read blocks of a few sequences).
    #[cfg(test)]
    fn with_block_bytes(db: &'a SeqDb, longest: usize, block_bytes: u64) -> Self {
        Self::with_block_bytes_and_count(db, longest, db.num_seqs() as usize, block_bytes)
    }

    fn with_block_bytes_and_count(
        db: &'a SeqDb,
        longest: usize,
        count: usize,
        block_bytes: u64,
    ) -> Self {
        Self {
            db,
            count,
            longest,
            block_bytes,
            x_checked: std::sync::atomic::AtomicBool::new(false),
            x_covered: std::sync::Mutex::new(Vec::new()),
            x_damaged: std::sync::OnceLock::new(),
            x_first: std::sync::Mutex::new(None),
            #[cfg(feature = "parallel")]
            x_task_starts: std::sync::Mutex::new(Vec::new()),
        }
    }

    // No NCBI counterpart: LOSAT's check of the stored bytes (NCBI searches them as they are; the
    // check is the one `x_check_stores` makes of the whole database). A read of the whole database
    // before the first query batch made it once: 0.85 s for UniRef50's uniref50.03 volume (3 GB of
    // residues) at every thread count, on one thread while the others waited (0.84 of 2.22 s for
    // one query at 16 threads), and the search then read every sequence again.
    /// Checks the sequences of `part`, just read, as `x_check_stores` checks them
    /// (`check_protein_sequence`, `MAX_DBSEQ_LEN`) until the reads have checked every sequence of
    /// the database; after that, nothing (the bytes are those already checked). A sequence that
    /// fails is an `XReadError`, which the run reports as the check's first error
    /// (`native::search_error`).
    fn x_check(&self, part: &SequenceBlock) -> Result<()> {
        if self.x_checked.load(std::sync::atomic::Ordering::Acquire) {
            return Ok(());
        }
        let first = part.first_ordinal();
        for i in 0..part.len() {
            let with_sentinels = part.protein_with_sentinels(i);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1425-1432
            // Proteins with NCBI 32-bit length <= 0 (`0` or `>= 2^31`) are skipped by NCBI and
            // represented as `[0, 0]` by `read_protein_search_block`.
            if with_sentinels.len() == 2 {
                continue;
            }
            let ordinal = first + i as u32;
            check_protein_sequence(with_sentinels, || self.db.oid(ordinal), PROGRAM)
                .map_err(XReadError)?;
            let length = with_sentinels.len() - 2;
            if length > MAX_DBSEQ_LEN {
                return Err(XReadError(too_long(
                    format_args!("database sequence {}", self.db.oid(ordinal)),
                    length,
                ))
                .into());
            }
        }
        self.x_cover(first as usize, first as usize + part.len());
        Ok(())
    }

    // No NCBI counterpart: the bookkeeping of `x_check`.
    /// Adds the checked ordinals `first..end` to `x_covered`; once the ranges cover the database,
    /// sets `x_checked` and drops them.
    fn x_cover(&self, first: usize, end: usize) {
        let mut covered = self
            .x_covered
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let mut at = covered.partition_point(|&(start, _)| start < first);
        covered.insert(at, (first, end));
        if at > 0 && covered[at - 1].1 >= first {
            covered[at - 1].1 = covered[at - 1].1.max(end);
            covered.remove(at);
            at -= 1;
        }
        while at + 1 < covered.len() && covered[at + 1].0 <= covered[at].1 {
            covered[at].1 = covered[at].1.max(covered[at + 1].1);
            covered.remove(at + 1);
        }
        if covered.len() == 1 && covered[0] == (0, self.count) {
            self.x_checked
                .store(true, std::sync::atomic::Ordering::Release);
            *covered = Vec::new();
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.count
    }

    /// The database read.
    pub(crate) fn seqdb(&self) -> &'a SeqDb {
        self.db
    }

    /// `GetTotalLength` (`SeqDb::total_length`).
    pub(crate) fn total_length(&self) -> u64 {
        self.db.total_length()
    }

    /// The length of the included sequences (`SeqDb::included_length`), which the schedules read.
    #[cfg(feature = "parallel")]
    pub(crate) fn included_length(&self) -> Result<u64, crate::seqdb::SeqDbError> {
        self.db.included_length()
    }

    /// The numbers of the statistics, with NCBI's types (`SeqDb::statistics`).
    pub(crate) fn statistics(&self) -> crate::seqdb::DbStatistics {
        self.db.statistics()
    }

    pub(crate) fn longest(&self) -> usize {
        self.longest
    }

    /// Whether the reads have checked every sequence (`x_check`).
    pub(crate) fn x_checked(&self) -> bool {
        self.x_checked.load(std::sync::atomic::Ordering::Acquire)
    }

    // No NCBI counterpart: LOSAT's check of the stored bytes, on the paths of an error only (an
    // error raised before the reads have checked every sequence).
    /// Whether the database fails the whole check (`database::x_check_whole`: the first error in
    /// ordinal order, which the check before the first query batch gave), run at most once for the
    /// stream; its error is kept for `x_first_error`.
    pub(crate) fn x_damaged(&self) -> bool {
        *self.x_damaged.get_or_init(|| match x_check_whole(self.db) {
            Ok(_) => false,
            Err(error) => {
                *self
                    .x_first
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(error);
                true
            }
        })
    }

    /// The whole check's first error (`x_damaged`), `None` where it finds none.
    pub(crate) fn x_first_error(&self) -> Option<anyhow::Error> {
        if !self.x_damaged() {
            return None;
        }
        let kept = self
            .x_first
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .take();
        kept.or_else(|| x_check_whole(self.db).err())
    }

    /// `CSeqDB::GetMinLength` (`SeqDb::min_length`).
    pub(crate) fn min_length(&self) -> i32 {
        self.db.min_length()
    }

    /// The block of a pass that starts at `ordinal`: as many sequences as fit `block_bytes`
    /// (`SeqDb::next_block`, at least one).
    fn read_from(&self, ordinal: usize) -> Result<XBlock> {
        let range = u32::try_from(ordinal)
            .ok()
            .and_then(|first| self.db.next_block(first, self.block_bytes))
            .ok_or_else(|| anyhow!("BLASTX database ordinal id {ordinal} is past the database"))?;
        let part = self
            .db
            .read_protein_search_block(range)
            .map_err(read_error)?;
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:339-344,699-707
        // A multithreaded DB-split block may end early at a negative-length OID, so `end` reflects
        // `part.len()` rather than `range.end_ordinal()`.
        let end = ordinal + part.len();
        self.x_check(&part)?;
        Ok(XBlock {
            first: ordinal,
            end,
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
            let part = self
                .db
                .read_protein_search_block(range)
                .map_err(read_error)?;
            self.x_check(&part)?;
            parts.push(part);
            at = range.end_ordinal();
        }
        Ok(XBlock { first, end, parts })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1425-1432
    // No NCBI counterpart: LOSAT's length window for parallel work scheduling; proteins whose
    // signed 32-bit `GetSeqLengthProt` length is non-positive are skipped by NCBI and counted as
    // length 0.
    /// The lengths of the sequences from `oid` on, at most `LENGTH_WINDOW` of them (to the end of
    /// `oid`'s volume), from the index (`0` for proteins skipped by NCBI with length `<= 0`).
    fn read_lengths(&self, oid: usize) -> Result<Vec<u32>> {
        let first = u32::try_from(oid)?;
        let end = first.saturating_add(LENGTH_WINDOW as u32);
        self.db
            .x_protein_lengths_i32(first, end)
            .map(|lengths| lengths.into_iter().map(|len| len.max(0) as u32).collect())
            .map_err(read_error)
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
            .read_protein_search_sequence(u32::try_from(oid)?)
            .map_err(read_error)?;
        self.x_check(&block)?;
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

    /// The length of `oid`: from the block, or a read of the database anywhere else (no stage
    /// does: the sink reads the lengths from its `XCursor`).
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
#[path = "x_db_stream_tests.rs"]
mod tests;
