//! The database of a `-db` search, read where it is searched, instead of a copy of the whole
//! database made before the first query batch: the core of TBLASTN's stream
//! (`tblastn/x_db_stream.rs` holds what only TBLASTN's stages do with it) and TBLASTX's
//! (`tblastx::subjects::Subjects::Stream`).
//!
//! A copy packs every sequence of the database into anonymous memory before the first query
//! batch, and it stays until the reports, where NCBI maps the sequence file.
//!
//! The stream:
//!
//! - `XDbStream::new` reads the database once in blocks of `STORE_BLOCK_BYTES`, with the checks
//!   of the copy in its order (`common::database::x_check_nucleotide_stores`). A block of one
//!   sequence (a long one) is checked without its 2-bit data: its last 2-bit byte and its
//!   ambiguity data only (`SeqDb::x_nucleotide_check`). Its error is the error that the copy
//!   gave, which fails the same query batch with the same text. The check keeps the number of
//!   sequences, the longest, the sum of the lengths and the first length (`XChecked`), not a
//!   length per sequence (a `u32` each was 529 MB on core_nt, held through the search; NCBI maps
//!   the index): a search reads the length of a sequence where it reads the sequence
//!   (`XDbStream::length`), a report from the index and the sequence's last byte
//!   (`XDbStream::read_length`), and TBLASTN's chunks of subjects are cut from the index's
//!   approximate lengths (`XDbStream::x_approx_lengths`).
//! - The blocks of 256 KiB of the sequence file are fixed once, from the index alone
//!   (`x_block_grid`), since `SeqDb::next_block` starts a block at the sequence asked for and
//!   threads that ask for different sequences would cut different blocks. The search of a short
//!   subject makes what it needs from the bytes of its block, without a store
//!   (`XDbStream::with_block`): each thread keeps the last block it used, and a block that another
//!   thread holds or reads is not read again. The threads that hand out the sequences one at a
//!   time read each block about once, and the error of a failed read (the text of a damaged block
//!   names the block's bytes) does not depend on the threads or their order. A thread that
//!   searches a chunk of subjects (`x_plan_chunk`) lets its block go with the chunk; the threads
//!   of the neighbouring chunks share the block that the chunk's end cuts across.
//! - A long sequence (a block of its own) is read once while anything holds it: the stages and
//!   the query batches that run at the same time share it, and a thread that asks for it while
//!   another reads it waits for that read (`XDbStream::read_shared`). With more than one thread
//!   it is read into a buffer that glibc maps (`utils::x_mapped`), taken from the pool
//!   `utils::x_mapped::X_BYTES` when it holds one, and its store is that buffer cut to the 2-bit
//!   data, not a copy; dropped, its pages go back instead of staying in the worker's arena. At one
//!   thread it is read into the heap, where a dropped store's pages serve the next read.
//! - While a stream is open (`x_frames_pooled`), a search may give its idle large buffers to
//!   `X_BYTES`. When the last stream of the process closes, the idle buffers are freed.
//! - TBLASTX makes the frames of a short subject from the bytes of its block, with the ambiguity
//!   runs that its store would keep (`XSequence::with_view`), and keeps the block for the frames
//!   of its letters when it has ambiguous bases; a long one is read alone (`read_shared`) for its
//!   frames, and again for its letters' frames when it has ambiguous bases and hits reach the
//!   re-evaluation. Each query batch reads the sequences it reports once each
//!   (`tblastx::database::read_reported`), when a report reads their letters.
//!
//! Every read packs the stored bytes as the copy did (`common::database::x_read_nucleotide_store`:
//! the packing of `PackedNucleotide::from_stored`), and a view of a block's bytes gives the letters
//! and bases of that store (`XSequence::with_view`), so every stage computes the same values on the
//! same subjects in the same order. A read that fails during the search (a database changed or
//! damaged after the check, an error of the file system) is an `XReadError` holding the error that
//! the check gives for the same read; the program decides how the batch that read it fails.
//! The check does not read the 2-bit data of a block of one sequence but its last byte, so an
//! error of the disk in the rest comes in the search's read of the sequence and not in the check:
//! TBLASTX fails the batch that reads it, if one does, as the copy failed the first batch with a
//! valid context, except that a batch without one that comes first and is reported with the batch
//! (`tblastx::report_stream`) is reported before LOSAT's refusal.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:374-376 (the sequence of
//! the search is the mapped volume's bytes, not a copy)

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock, Weak};

use anyhow::Result;

use crate::algorithm::common::database as shared;
use crate::algorithm::common::db_mask_search::{SearchMask, SubjectMasks};
use crate::algorithm::common::hard_chunk_walk::{HardStarts, NcbiStart, NcbiWalk};
use crate::seqdb::{AmbiguityRun, SeqDb, SeqDbError, SequenceBlock};
use crate::sequence::packed_subject::{x_stored_runs, PackedNucleotide, PackedRef};

/// A read of the database that failed during the search, with
/// the error that the check (`XDbStream::new`) gives for the same read.
#[derive(Debug)]
pub(crate) struct XReadError(pub(crate) anyhow::Error);

impl std::fmt::Display for XReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for XReadError {}

/// What the threads of a stream read once while anything holds it, by a number of its own: the
/// long sequences by ordinal id (`XDbStream::read_shared`), the blocks by their number in the grid
/// (`XDbStream::with_block`).
struct XSharedReads<T> {
    table: Mutex<XShared<T>>,
    /// Woken when a read ends.
    read_done: Condvar,
}

struct XShared<T> {
    /// What something holds.
    held: HashMap<usize, Weak<T>>,
    /// What a thread reads now.
    reading: HashSet<usize>,
}

/// A read of `XSharedReads::read_once`; dropped (the read made, failed or unwound), its key is no
/// longer read and the waiters wake.
struct XReading<'s, T> {
    reads: &'s XSharedReads<T>,
    key: usize,
}

impl<T> Drop for XReading<'_, T> {
    fn drop(&mut self) {
        let mut table = self.reads.table.lock().unwrap_or_else(|e| e.into_inner());
        table.reading.remove(&self.key);
        drop(table);
        self.reads.read_done.notify_all();
    }
}

impl<T> XSharedReads<T> {
    fn new() -> Self {
        Self {
            table: Mutex::new(XShared {
                held: HashMap::new(),
                reading: HashSet::new(),
            }),
            read_done: Condvar::new(),
        }
    }

    /// What `key` names: the value that something holds, or the value of the read that a thread
    /// makes now (this thread waits for it), or else `read` of this thread, whose value the next
    /// threads share. It is read again once nothing holds it, and after a read that failed. `read`
    /// runs no rayon work, so a waiting worker depends on no task that could be stolen below it.
    fn read_once(&self, key: usize, read: impl FnOnce() -> Result<T>) -> Result<Arc<T>> {
        let mut table = self.table.lock().unwrap_or_else(|e| e.into_inner());
        loop {
            if let Some(value) = table.held.get(&key).and_then(Weak::upgrade) {
                return Ok(value);
            }
            if !table.reading.contains(&key) {
                break;
            }
            table = self
                .read_done
                .wait(table)
                .unwrap_or_else(|e| e.into_inner());
        }
        table.reading.insert(key);
        drop(table);
        let reading = XReading { reads: self, key };
        let value = Arc::new(read()?);
        let mut table = self.table.lock().unwrap_or_else(|e| e.into_inner());
        // A dead `Weak` still holds its allocation.
        if table.held.len() >= 64 && table.held.len().is_power_of_two() {
            table.held.retain(|_, value| value.strong_count() > 0);
        }
        table.held.insert(key, Arc::downgrade(&value));
        drop(table);
        drop(reading);
        Ok(value)
    }
}

/// A failure of [`XDbStream::lead_start`]'s one-thread walk, kept in the stream's `OnceLock` so
/// that every caller sees the same read error (`XReadError` via `search_read_error_of`) or
/// internal error.
enum XLeadWalkError {
    Read(SeqDbError),
    Other(String),
}

impl From<SeqDbError> for XLeadWalkError {
    fn from(error: SeqDbError) -> Self {
        Self::Read(error)
    }
}

impl From<anyhow::Error> for XLeadWalkError {
    fn from(error: anyhow::Error) -> Self {
        Self::Other(error.to_string())
    }
}

/// The database of a `-db` search, read where it is searched.
pub(crate) struct XDbStream<'a> {
    db: &'a SeqDb,
    /// The program that the errors of the reads name (`common::database`).
    program: &'static str,
    /// A number of its own in the process (`XPlanBlock::stream`).
    id: u64,
    /// The number of sequences and their lengths summed up (`x_check_nucleotide_stores`).
    checked: shared::XChecked,
    /// The ordinal ids that the blocks start at, and the number of sequences (`x_block_grid`).
    grid: Vec<u32>,
    /// Whether a sequence is read into a buffer that glibc maps (`x_read_nucleotide_store`): with
    /// more than one thread. The freed stores of the workers stayed resident in their arenas; at
    /// one thread the heap's freed store serves the next read without faults (TBLASTN, C. elegans
    /// genome: 19 reads of its 7 sequences, against one read of the copy).
    mapped: bool,
    /// The long sequences that something holds or a thread reads, by ordinal id (`read_shared`).
    stores: XSharedReads<PackedNucleotide>,
    /// The blocks that something holds or a thread reads, by their number in `grid` (`with_block`).
    blocks: XSharedReads<SequenceBlock>,
    /// The database masking of the search (`-db_soft_mask`, `-db_hard_mask`; series NP, M4b PR 2),
    /// `None` without the options (`with_mask`).
    mask: Option<SearchMask>,
    /// With `-db_soft_mask`, NCBI's one-thread diagonal-table start of every lead sequence
    /// (`lead_start`, series NP, M4b-T3), computed once at the first lead sequence of the run
    /// (`word`, `window`, recorded starts or walk error). Never populated for searches without a
    /// lead sequence.
    lead_starts: OnceLock<(usize, i32, std::result::Result<HardStarts, XLeadWalkError>)>,
    /// With `-db_hard_mask` (series NP, M4b-3c part 2), the diagonal-table starts computed by
    /// `ResolvedMask::check_losat_limits` when the database contains overlapping masks (`None`
    /// when no walk ran).
    hard_starts: Option<HardStarts>,
    /// With `-db_hard_mask`, the sequences as the stages after the preliminary one read them
    /// (`read_stage`), by ordinal id, while something holds them.
    masked_stores: XSharedReads<PackedNucleotide>,
    /// The cuts made of the sequences, by their key (`x_kept_cut`).
    #[cfg(any(feature = "parallel", test))]
    cuts: Mutex<Vec<(u64, Arc<[u32]>)>>,
    /// The blocks read (`with_block`).
    #[cfg(test)]
    block_reads: AtomicUsize,
    /// The sequences whose reads fail, and whether with LOSAT's refusal (`x_failing`).
    #[cfg(test)]
    x_failing: (Vec<usize>, bool),
}

impl<'a> XDbStream<'a> {
    /// Checks `db` as the copy of the database was read (`common::database::x_check_nucleotide_stores`)
    /// and fixes the blocks (`x_block_grid`). The errors name `program`.
    pub(crate) fn new(db: &'a SeqDb, threads: usize, program: &'static str) -> Result<Self> {
        static STREAMS: AtomicU64 = AtomicU64::new(0);
        let checked = shared::x_check_nucleotide_stores(db, program)?;
        let grid = x_block_grid(db, PLAN_BLOCK_BYTES);
        X_OPEN_STREAMS.fetch_add(1, Ordering::Relaxed);
        Ok(Self {
            db,
            program,
            id: STREAMS.fetch_add(1, Ordering::Relaxed),
            checked,
            grid,
            mapped: threads > 1,
            stores: XSharedReads::new(),
            blocks: XSharedReads::new(),
            mask: None,
            lead_starts: OnceLock::new(),
            hard_starts: None,
            masked_stores: XSharedReads::new(),
            #[cfg(any(feature = "parallel", test))]
            cuts: Mutex::new(Vec::new()),
            #[cfg(test)]
            block_reads: AtomicUsize::new(0),
            #[cfg(test)]
            x_failing: (Vec::new(), false),
        })
    }

    /// The stream of a search masked with `mask` (`-db_soft_mask`, `-db_hard_mask`; `None`: not
    /// masked), whose algorithm the database's volumes know (`SeqDb::check_search_masks`, before
    /// the search).
    pub(crate) fn with_mask(mut self, mask: Option<SearchMask>) -> Self {
        self.mask = mask;
        self
    }

    /// The database masking of the search (`with_mask`).
    pub(crate) fn mask(&self) -> Option<SearchMask> {
        self.mask
    }

    /// Attaches the `-db_hard_mask` diagonal-table starts from `ResolvedMask::translated_starts`
    /// (`None` when the database has no overlapping masks).
    pub(crate) fn with_hard_starts(mut self, hard_starts: Option<HardStarts>) -> Self {
        self.hard_starts = hard_starts;
        self
    }

    /// The `-db_hard_mask` diagonal-table start for sequence `ordinal`, if this sequence runs in
    /// order on its own table (a walked sequence or a sequence that meets unwritten entries or
    /// clears inside it in NCBI's one-thread search; `None` for all other sequences).
    pub(crate) fn hard_start(&self, ordinal: usize) -> Option<NcbiStart> {
        let starts = self.hard_starts.as_ref()?;
        let ordinal = u32::try_from(ordinal).ok()?;
        starts.of(ordinal)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-483,804-855,
    // 1318-1331,1368-1369,1377,1407-1432
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107,162-173
    // Each search thread makes one diagonal table for the database (`s_BlastSetUpAuxStructures`),
    // and every query batch starts from a fresh table walking the same searched sequences in
    // database order, so a lead sequence's one-thread start (`NcbiStart`) is fixed for the run.
    /// NCBI's one-thread diagonal-table start ([`NcbiStart`]) of lead sequence `ordinal` in a
    /// `-db_soft_mask` translated search of word size `word` and two-hit window `window` (0 for
    /// one-hit). On the first call of the run, walks every searched sequence in database order
    /// with its own mask reader (`SeqDb::searched_sequences`, never `X_MASK_READERS`) and records
    /// the start of every lead sequence; later calls across threads and query batches read that
    /// table. Requesting an ordinal that the walk did not record is an internal error; read
    /// failures use the preliminary stage's errors (`search_read_error_of`, `XReadError`).
    #[cold]
    #[inline(never)]
    pub(crate) fn lead_start(&self, ordinal: usize, word: usize, window: i32) -> Result<NcbiStart> {
        let Some(mask) = self.mask.filter(|mask| !mask.hard()) else {
            anyhow::bail!(
                "LOSAT internal error: lead_start requested on a stream without database soft masking"
            );
        };
        let (walk_word, walk_window, result) = self.lead_starts.get_or_init(|| {
            (
                word,
                window,
                self.walk_lead_starts(mask.algo_id, word, window),
            )
        });
        debug_assert_eq!(
            (*walk_word, *walk_window),
            (word, window),
            "word and window are fixed for a run"
        );
        let starts = match result {
            Ok(starts) => starts,
            Err(XLeadWalkError::Read(error)) => {
                return Err(XReadError(shared::search_read_error_of(error, self.program)).into());
            }
            Err(XLeadWalkError::Other(message)) => anyhow::bail!("{message}"),
        };
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        starts.of(ordinal).ok_or_else(|| {
            anyhow::anyhow!(
                "LOSAT internal error: database sequence {} was not recorded as a soft-mask lead sequence",
                self.db.oid(ordinal)
            )
        })
    }

    /// Walks every searched sequence in database order with [`NcbiWalk`] and records the
    /// [`NcbiStart`] of each `-db_soft_mask` lead sequence (`lead_start`).
    fn walk_lead_starts(
        &self,
        algo_id: i32,
        word: usize,
        window: i32,
    ) -> std::result::Result<HardStarts, XLeadWalkError> {
        let mut walk = NcbiWalk::new(i32::try_from(word).unwrap_or(i32::MAX), window);
        let mut starts = HardStarts::default();
        let mut exits = Vec::new();
        let mut ranges = Vec::new();
        let mut split = crate::core::subject_split::HardSplitState::default();
        let mut walk_error: Option<XLeadWalkError> = None;
        let walked = self
            .db
            .searched_sequences(algo_id, |ordinal, length, masks| {
                let lead =
                    match crate::algorithm::common::translated_chunks::translated_soft_mask_exits(
                        length,
                        word,
                        &mut exits,
                        &mut split,
                        || -> std::result::Result<Cow<'_, [(i32, i32)]>, XLeadWalkError> {
                            let stored = masks.read()?;
                            crate::core::subject_split::complement_masks(
                                stored,
                                length,
                                &mut ranges,
                            );
                            Ok(Cow::Borrowed(ranges.as_slice()))
                        },
                    ) {
                        Ok(lead) => lead,
                        Err(XLeadWalkError::Read(error)) => return Err(error),
                        Err(error @ XLeadWalkError::Other(_)) => {
                            walk_error = Some(error);
                            return Ok(Some(()));
                        }
                    };
                let start = walk.translated(exits.iter().copied());
                if lead {
                    starts.push(ordinal, start);
                }
                Ok(None::<()>)
            });
        if let Some(error) = walk_error {
            return Err(error);
        }
        walked?;
        Ok(starts)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:315-320, 393-402 and
    // c++/src/algo/blast/core/blast_util.c:182-217
    // In a masked search every fetched sequence carries the complement of its masks as its
    // ranges (one range, the whole sequence, for a sequence without masks), and the masking kind.
    /// The ranges in bases of sequence `ordinal` of `length` bases in a masked search (`with_mask`):
    /// the complement of its stored masks of the search's algorithm, in the order stored
    /// (`complement_masks`; `[0, length)` for a sequence without masks); `None` in a search
    /// without the options. The preliminary stage's read of the masks: an error is an
    /// `XReadError` with the error of the search's reads (`search_read_error_of`).
    pub(crate) fn mask_ranges(
        &self,
        ordinal: usize,
        length: usize,
    ) -> Result<Option<Vec<(i32, i32)>>> {
        let Some(mask) = self.mask else {
            return Ok(None);
        };
        let mut ranges = Vec::new();
        self.with_stored_masks(ordinal, mask, true, |masks| {
            crate::core::subject_split::complement_masks(masks, length, &mut ranges);
        })?;
        Ok(Some(ranges))
    }

    /// `f` of the stored masks of sequence `ordinal` for the algorithm of `mask`, read by the
    /// thread's reader of the preliminary stage (`preliminary`: a window of blocks, with the
    /// search's read errors) or of the later stages (one sequence at a time, with the
    /// database's errors), each kept by the thread for this stream.
    fn with_stored_masks<R>(
        &self,
        ordinal: usize,
        mask: SearchMask,
        preliminary: bool,
        f: impl FnOnce(&[(u32, u32)]) -> R,
    ) -> Result<R> {
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        X_MASK_READERS.with(|cell| {
            let mut slot = cell.borrow_mut();
            if slot
                .as_ref()
                .is_none_or(|readers| readers.stream != self.id)
            {
                *slot = Some(XMaskReaders {
                    stream: self.id,
                    preliminary: SubjectMasks::default(),
                    later: SubjectMasks::with_window(1),
                });
            }
            let readers = slot.as_mut().expect("the thread's readers were set");
            let (reader, error): (_, fn(&crate::seqdb::SeqDbError, &str) -> anyhow::Error) =
                if preliminary {
                    (&mut readers.preliminary, shared::search_read_error_of)
                } else {
                    (&mut readers.later, shared::database_error_of)
                };
            let masks = reader
                .of(self.db, ordinal, mask.algo_id, error, self.program)
                .map_err(|error| anyhow::Error::from(XReadError(error)))?;
            Ok(f(masks))
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:322-370 and
    // c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1539, 1580-1589
    // With hard masks, every fetch after the preliminary stage (the traceback's, the composition
    // statistics', the ungapped re-evaluation's) is a copied ncbi4na fetch whose masked bases are
    // B, over the ambiguities; a partial fetch masks each slice it reads the same way.
    /// Sequence `ordinal` as the stages after the preliminary one read it: with `-db_hard_mask`, its
    /// store with `B` over its stored masks (`db_mask_search::mask_nucleotide_store`), shared
    /// as `read_shared` shares the store; otherwise `read_shared`.
    pub(crate) fn read_stage(&self, ordinal: usize) -> Result<Arc<PackedNucleotide>> {
        let Some(mask) = self.mask.filter(|mask| mask.hard()) else {
            return self.read_shared(ordinal);
        };
        self.masked_stores.read_once(ordinal, || {
            let store = self.read_shared(ordinal)?;
            self.with_stored_masks(ordinal, mask, false, |masks| {
                crate::algorithm::common::db_mask_search::mask_nucleotide_store(&store, masks)
            })
        })
    }

    /// The number of sequences.
    pub(crate) fn count(&self) -> usize {
        self.checked.count
    }

    /// The database read.
    pub(crate) fn seqdb(&self) -> &'a SeqDb {
        self.db
    }

    /// The number of bases of the longest sequence and of all the sequences.
    pub(crate) fn longest_and_total(&self) -> (usize, u64) {
        (self.checked.longest as usize, self.checked.total)
    }

    /// The number of bases of the first sequence (none without one), from the check.
    pub(crate) fn first_length(&self) -> usize {
        self.checked.first as usize
    }

    /// The number of bases of sequence `ordinal`, read where the search reads the sequence: a short
    /// one's from its block, which this thread keeps for the `with_block` that follows; a long
    /// one's (a block of its own) from the index and its last byte (`read_length`), without its
    /// sequence. A read that fails is an `XReadError` with the error of the sequence's read
    /// (`with_block`, `read`). Only the length is read of the block: no ambiguity data.
    pub(crate) fn length(&self, ordinal: usize) -> Result<usize> {
        self.with_grid_block(ordinal, |held| match held {
            None => self.read_length(ordinal),
            Some((block, i)) => block
                .nucleotide_len(i)
                .map_err(|error| XReadError(shared::search_read_error(error, self.program)).into()),
        })
    }

    /// The number of bases of sequence `ordinal` read from the index and the last byte of its 2-bit
    /// data (`common::database::x_nucleotide_length`, NCBI's `GetSeqLength`), without a block or
    /// the sequence: the length of the sequence's store, for the reports that read no letters.
    /// A read that fails is an `XReadError`.
    pub(crate) fn read_length(&self, ordinal: usize) -> Result<usize> {
        #[cfg(test)]
        self.x_failed(ordinal..ordinal + 1)?;
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        shared::x_nucleotide_length(self.db, ordinal, self.program)
            .map_err(|error| XReadError(error).into())
    }

    /// Calls `visit` with the length of every sequence in ordinal order as the index gives it
    /// without a read of the sequence file (`SeqDb::x_approx_lengths`: NCBI's approximate
    /// length, at most 3 bases off the stored one, the low bits of the OID for the last byte's).
    /// When the index cannot be read (it changed after the check), the sequences from there on
    /// are given by the blocks of the grid: the first of a block with the bases of a whole block
    /// and the others with none. For a cut of the sequences that is scheduling only
    /// (`x_kept_cut`).
    #[cfg(any(feature = "parallel", test))]
    pub(crate) fn x_approx_lengths(&self, mut visit: impl FnMut(u32)) {
        let count = u32::try_from(self.count()).expect("a database ordinal id fits Uint4");
        let mut visited = 0u32;
        let _ = self.db.x_approx_lengths(|length| {
            if visited < count {
                visit(length);
                visited += 1;
            }
        });
        self.x_grid_lengths(visited, visit);
    }

    /// `x_approx_lengths` of the sequences from ordinal id `from` on that the index did not give:
    /// a block of the grid as its first sequence with the bases of a whole block and the others
    /// with none.
    #[cfg(any(feature = "parallel", test))]
    fn x_grid_lengths(&self, from: u32, mut visit: impl FnMut(u32)) {
        let count = u32::try_from(self.count()).expect("a database ordinal id fits Uint4");
        let whole = u32::try_from(PLAN_BLOCK_BYTES * 4).unwrap_or(u32::MAX);
        let mut next = self.grid.partition_point(|&first| first < from);
        for ordinal in from..count {
            let starts = self.grid.get(next) == Some(&ordinal);
            next += usize::from(starts);
            visit(if starts { whole } else { 0 });
        }
    }

    /// The cut of the sequences that `make` gives for `key` (the ends of TBLASTN's chunks of
    /// subjects for their bases, from `x_approx_lengths`), made by the first caller with `key`
    /// and kept for the run: every query batch of the search and every pool (each `key`) cuts
    /// the same sequences, which on core_nt read 1 GB of the index.
    #[cfg(any(feature = "parallel", test))]
    pub(crate) fn x_kept_cut(&self, key: u64, make: impl FnOnce() -> Vec<u32>) -> Arc<[u32]> {
        let mut cuts = self
            .cuts
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some((_, cut)) = cuts.iter().find(|(kept, _)| *kept == key) {
            return Arc::clone(cut);
        }
        let cut: Arc<[u32]> = make().into();
        cuts.push((key, Arc::clone(&cut)));
        cut
    }

    /// Sequence `ordinal`, packed as the copy of the database packed it, with one read (`XReadError` if it fails).
    pub(crate) fn read(&self, ordinal: usize) -> Result<PackedNucleotide> {
        #[cfg(test)]
        self.x_failed(ordinal..ordinal + 1)?;
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        shared::x_read_nucleotide_store(self.db, ordinal, self.program, self.mapped)
            .map_err(|error| XReadError(error).into())
    }

    /// Sequence `ordinal` (`read`), shared with whatever holds it already: a stage or a query batch
    /// that runs at the same time reads it once, and a thread that asks for it while another
    /// reads it waits for that read (the threads of a stage that first align the same subject at
    /// once made one copy each). It is read again once nothing holds it.
    pub(crate) fn read_shared(&self, ordinal: usize) -> Result<Arc<PackedNucleotide>> {
        self.stores.read_once(ordinal, || self.read(ordinal))
    }

    /// `make` of sequence `ordinal` (`XSequence`), from the bytes of its block (`x_block_grid`), with
    /// the errors of its store but without the store. The thread keeps the last block it used; a
    /// block that another thread holds or reads is shared (`XSharedReads::read_once`). A block of
    /// one sequence (a long one) is read alone and shared (`read_shared`). The thread's block is
    /// out of its cell while `make` runs: a nested call (a pool thread that takes another subject
    /// while it waits) does not borrow it. A thread that searches a chunk of subjects
    /// (`x_plan_chunk`) lets its block go with the chunk.
    pub(crate) fn with_block<R>(
        &self,
        ordinal: usize,
        make: impl FnOnce(XSequence<'_>) -> Result<R>,
    ) -> Result<R> {
        self.with_grid_block(ordinal, |held| match held {
            None => make(XSequence::Store(&self.read_shared(ordinal)?)),
            Some((block, i)) => {
                let (length, runs) = self.parts_in(block, i, ordinal)?;
                make(XSequence::Block {
                    block,
                    i,
                    length,
                    runs,
                })
            }
        })
    }

    /// `f(None)` when sequence `ordinal` is a block of its own (a long one), else `f` of its block
    /// and its place in the block, with the block kept by this thread as `with_block` keeps it:
    /// out of the thread's cell while `f` runs, and back in it after.
    fn with_grid_block<R>(
        &self,
        ordinal: usize,
        f: impl FnOnce(Option<(&Arc<SequenceBlock>, usize)>) -> Result<R>,
    ) -> Result<R> {
        let ordinal = u32::try_from(ordinal).expect("a database ordinal id fits Uint4");
        let mut kept = X_PLAN_BLOCK.with(std::cell::Cell::take);
        let holds = |kept: &Option<XPlanBlock>| {
            kept.as_ref().is_some_and(|kept| {
                kept.stream == self.id
                    && (kept.first..kept.first + kept.block.len() as u32).contains(&ordinal)
            })
        };
        if !holds(&kept) {
            // Release the thread's previous block before the next one (or a long sequence) is
            // read, which keeps the peak memory down.
            drop(kept.take());
            let number = self.grid.partition_point(|&first| first <= ordinal) - 1;
            let (first, end) = (self.grid[number], self.grid[number + 1]);
            if end - first == 1 {
                return f(None);
            }
            // LOSAT_X_TBNPLANCUT: the block goes with this thread's chunk.
            crate::utils::xreach::x_reach!(
                "LOSAT_X_TBNPLANCUT",
                X_PLAN_CHUNKS.with(std::cell::Cell::get) > 0
            );
            kept = Some(XPlanBlock {
                stream: self.id,
                first,
                block: self
                    .blocks
                    .read_once(number, || self.read_grid_block(first, end))?,
            });
        }
        let block = kept.as_ref().expect("a block holds the sequence");
        let made = f(Some((&block.block, (ordinal - block.first) as usize)));
        X_PLAN_BLOCK.with(|cell| cell.set(kept));
        made
    }

    /// The stream with the blocks of `max_bytes` (`x_block_grid`), for the tests of a program:
    /// the fixtures are smaller than a block of the search.
    #[cfg(test)]
    pub(crate) fn x_regrid(mut self, max_bytes: u64) -> Self {
        self.grid = x_block_grid(self.db, max_bytes);
        self
    }

    /// The stream whose reads of the sequences `ordinals` fail (a block's read when it holds one of
    /// them) with NCBI's read error, or with LOSAT's refusal: a read that fails after the check,
    /// which the files of a test cannot give (the check keeps their handles and lengths).
    #[cfg(test)]
    pub(crate) fn x_failing(mut self, ordinals: &[usize], refusal: bool) -> Self {
        self.x_failing = (ordinals.to_vec(), refusal);
        self
    }

    /// The error of a read of the sequences `ordinals` that `x_failing` names one of.
    #[cfg(test)]
    fn x_failed(&self, ordinals: std::ops::Range<usize>) -> Result<()> {
        let (failing, refusal) = &self.x_failing;
        let error = match failing.iter().find(|ordinal| ordinals.contains(ordinal)) {
            None => return Ok(()),
            Some(&ordinal) if *refusal => shared::refusal(
                &format!(
                    "the read of database sequence {} failed",
                    self.db.oid(ordinal as u32)
                ),
                self.program,
            ),
            Some(_) => shared::PreliminaryReadError("Database memory map file error").into(),
        };
        Err(XReadError(error).into())
    }

    /// The block of the sequences `first..end` of the grid (`XReadError` if it cannot be read).
    fn read_grid_block(&self, first: u32, end: u32) -> Result<SequenceBlock> {
        #[cfg(test)]
        self.block_reads.fetch_add(1, Ordering::Relaxed);
        #[cfg(test)]
        self.x_failed(first as usize..end as usize)?;
        let range = self
            .db
            .x_block_range(first, end)
            .expect("a block of the grid is a range of the database");
        self.db
            .read_block(range)
            .map_err(|error| XReadError(shared::search_read_error(error, self.program)).into())
    }

    /// The length and the ambiguity runs of sequence `ordinal`, sequence `i` of `block`, with the
    /// errors of its store (`XReadError`).
    pub(crate) fn parts_in(
        &self,
        block: &SequenceBlock,
        i: usize,
        ordinal: usize,
    ) -> Result<(usize, Vec<AmbiguityRun>)> {
        shared::x_nucleotide_parts_in(block, i, || self.db.oid(ordinal as u32), self.program)
            .map_err(|error| XReadError(error).into())
    }
}

/// A sequence that `XDbStream::with_block` gives.
pub(crate) enum XSequence<'b> {
    /// Sequence `i` of a block of several, with its length and its ambiguity runs as stored.
    Block {
        block: &'b Arc<SequenceBlock>,
        i: usize,
        length: usize,
        runs: Vec<AmbiguityRun>,
    },
    /// A long sequence (a block of its own), read alone and shared (`XDbStream::read_shared`).
    Store(&'b Arc<PackedNucleotide>),
}

impl XSequence<'_> {
    /// The packed bytes, as the store's (`PackedNucleotide::packed`).
    pub(crate) fn packed(&self) -> &[u8] {
        match self {
            Self::Block { block, i, .. } => block.nucleotide_packed(*i),
            Self::Store(store) => store.packed(),
        }
    }

    /// The number of bases.
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Block { length, .. } => *length,
            Self::Store(store) => store.len(),
        }
    }

    /// `f` of the sequence as its store reads it (`PackedRef`): the store's view, or the block's
    /// bytes with the runs that the store keeps (`x_stored_runs`).
    pub(crate) fn with_view<R>(self, f: impl FnOnce(PackedRef<'_>) -> R) -> R {
        match self {
            Self::Block {
                block,
                i,
                length,
                mut runs,
            } => {
                let packed = block.nucleotide_packed(i);
                x_stored_runs(packed, length, &mut runs);
                f(PackedRef::x_stored(packed, length, &runs))
            }
            Self::Store(store) => f(store.view()),
        }
    }
}

impl Drop for XDbStream<'_> {
    /// The calling thread (the pool's first worker) outlives the run: its block goes with the
    /// stream. The pool's other threads end with the pool.
    fn drop(&mut self) {
        let _ = X_PLAN_BLOCK.try_with(std::cell::Cell::take);
        if X_OPEN_STREAMS.fetch_sub(1, Ordering::Relaxed) == 1 {
            crate::utils::x_mapped::X_BYTES.x_clear_idle();
        }
    }
}

/// The streams open in the process (`XDbStream::new` to its drop).
static X_OPEN_STREAMS: AtomicUsize = AtomicUsize::new(0);

/// A thread's readers of the stored masks of one stream (`XDbStream::with_stored_masks`).
#[derive(Debug)]
struct XMaskReaders {
    /// The stream's number (`XDbStream::id`).
    stream: u64,
    /// The reader of the preliminary stage (a window of sequences).
    preliminary: SubjectMasks,
    /// The reader of the later stages (one sequence).
    later: SubjectMasks,
}

thread_local! {
    /// The thread's readers of masks, for the last masked stream it read.
    static X_MASK_READERS: std::cell::RefCell<Option<XMaskReaders>> =
        const { std::cell::RefCell::new(None) };
}

/// Whether a `-db` stream is open in the process: then TBLASTN's large frames of a dropped subject
/// plan go to `utils::x_mapped::X_BYTES`, from which the next frames and reads take them
/// (`tblastn/x_db_stream.rs`). A `-subject` search that runs at the same time (library use) pools
/// its frames too; a frame taken from the pool is zero-filled before it is written, so no report
/// changes. No NCBI counterpart: only where the buffers go.
pub(crate) fn x_frames_pooled() -> bool {
    X_OPEN_STREAMS.load(Ordering::Relaxed) > 0
}

/// No NCBI counterpart: the bytes of the sequence file in a block of the grid
/// (`XDbStream::with_block`; the blocks of BLASTX's `-db` stream are 256 KiB too).
const PLAN_BLOCK_BYTES: u64 = 256 << 10;

/// No NCBI counterpart: where the reads of a stream start, and no value depends on them. The
/// ordinal ids that the blocks of at most `max_bytes` of the sequence file start at, chained from
/// the first sequence (`SeqDb::next_block`, which cuts them at the end of a volume and at the
/// sequences a list excludes), then the number of sequences. It reads the index, not the sequences.
fn x_block_grid(db: &SeqDb, max_bytes: u64) -> Vec<u32> {
    let mut grid = vec![0];
    let mut first = 0;
    while let Some(range) = db.next_block(first, max_bytes) {
        first = range.end_ordinal();
        grid.push(first);
    }
    grid.shrink_to_fit();
    grid
}

/// A block of sequences that a thread keeps (`XDbStream::with_block`), and its stream.
struct XPlanBlock {
    /// `XDbStream::id` of the stream that read it.
    stream: u64,
    first: u32,
    block: Arc<SequenceBlock>,
}

thread_local! {
    /// This thread's block.
    static X_PLAN_BLOCK: std::cell::Cell<Option<XPlanBlock>> = const { std::cell::Cell::new(None) };
    /// The chunks of subjects this thread searches now (`x_plan_chunk`; they nest).
    static X_PLAN_CHUNKS: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
}

// No NCBI counterpart (memory only; every sequence is made from the same bytes): reads the
// LOSAT_X_TBNPLANCUT variable once.
/// LOSAT_X_TBNPLANCUT, in the strict set (on unless `=0`): whether a thread's block goes with the
/// chunk of subjects the thread searches (`x_plan_chunk`).
#[cfg(any(feature = "parallel", test))]
fn x_plan_cut() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNPLANCUT"))
}

/// A thread's search of a chunk of subjects (`x_plan_chunk`), until it is dropped.
#[cfg(any(feature = "parallel", test))]
pub(crate) struct XPlanChunk(());

#[cfg(any(feature = "parallel", test))]
impl Drop for XPlanChunk {
    fn drop(&mut self) {
        X_PLAN_CHUNKS.with(|chunks| chunks.set(chunks.get() - 1));
        let _ = X_PLAN_BLOCK.try_with(std::cell::Cell::take);
    }
}

/// The calling thread searches a chunk of subjects until the value returned is dropped, and the
/// block it holds then goes with the chunk. A chunk holds about a sixteenth of the sequence file
/// over the threads and a thread's next chunk is about a sixteenth of the file further on, so past
/// 4 MiB of file the next chunk is outside the block; without this a thread kept its last block
/// until its next read, through the later stages of the search. The block is the grid's
/// (`x_block_grid`), not cut at the chunk's end: the threads of the neighbouring chunks, searched
/// at the same time, take the rest of it from the shared table (`XDbStream::with_block`), where a
/// block cut at a chunk's end would be the thread's own read. Chunks nest when a thread that waits
/// runs another flow's loop (query batches side by side): the inner one's end lets the outer block
/// go, and the outer chunk takes the block of its next sequence again. `None` with
/// `LOSAT_X_TBNPLANCUT=0`.
#[cfg(any(feature = "parallel", test))]
pub(crate) fn x_plan_chunk() -> Option<XPlanChunk> {
    x_plan_cut().then(|| {
        X_PLAN_CHUNKS.with(|chunks| chunks.set(chunks.get() + 1));
        XPlanChunk(())
    })
}

#[cfg(test)]
#[path = "x_db_stream_tests.rs"]
mod tests;
