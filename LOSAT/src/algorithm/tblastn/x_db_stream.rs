//! A `-db` search of TBLASTN reads each database sequence where it searches it (graduated
//! 2026-10-27; the switch was `LOSAT_X_TBNDBSTREAM`), instead of a copy of the whole database made
//! before the first query batch.
//!
//! The copy (`read_stores`, deleted) packed every sequence of the database into anonymous memory
//! before the first query batch, and every stage read the copy until the reports: 23.9 MiB for
//! the 100 Mb of the C. elegans genome, 32.2 MiB for a zebrafish cDNA set, alive at the peak of
//! every thread count (v4-D12 s8 probe, workstation), where NCBI maps the sequence file. The
//! reports then cloned the packed sequences they show (`database::reported_subjects`).
//!
//! The stream:
//!
//! - `XDbStream::new` reads the database once in blocks of `STORE_BLOCK_BYTES`, with the checks
//!   of the copy in its order (`common::database::x_check_nucleotide_stores`), and keeps the
//!   length of every sequence. A block of one sequence (a long one) is checked without its 2-bit
//!   data: its last 2-bit byte and its ambiguity data only (`SeqDb::x_nucleotide_check`). Its
//!   error is the error that the copy gave, which fails the same query batch with the same text.
//! - The preliminary stage of a subject reads its sequence when it makes the subject's plan
//!   (`search_seed::x_scan_plan`, `XPrelimSubject::Stream`) and drops it once the plan is made:
//!   the frames are translated from it and the stage reads nothing else of the sequence. A plan
//!   that the search keeps between query batches (`XSubjectSide`) is not made again, and its
//!   sequence is not read again. Each thread keeps the last block of 256 KiB of the sequence file
//!   it read for the plans (`XDbStream::plan_with`), from whose bytes the plans of the next short
//!   subjects are made, without a store.
//! - The composition redo and the traceback of `-comp_based_stats 0` read a subject of their
//!   lists when they first align it and hold it until the stage has done its last match
//!   (`XStageSubjects`, LOSAT_X_TBNSTAGEDROP; with `=0` until the stage ends); the reports read
//!   the sequences they show, one read each, when a report reads their letters
//!   (`database::reported_subjects`; the other reports read the lengths of the stream).
//! - A long sequence (a block of its own) is read once while anything holds it: the stages and
//!   the query batches that run at the same time share it, and a thread that asks for it while
//!   another reads it waits for that read (`XDbStream::read_shared`). With more than one thread
//!   it is read into a buffer that glibc maps (`utils::x_mapped`), taken from the pool
//!   `utils::x_mapped::X_BYTES` when it holds one, and its store is that buffer cut to the 2-bit
//!   data, not a copy; dropped, its pages go back instead of staying in the worker's arena. At one
//!   thread it is read into the heap, where a dropped store's pages serve the next read.
//! - A sequence whose plan has no residues (LOSAT_X_TBNUNITXLATE: a first frame of 1 MiB or more,
//!   in a search that keeps no plan of it) is read before its plan and held while its units
//!   translate it, as a packed subject of the copy is (`search_seed::x_lean_subject`,
//!   `search_gapped::x_preliminary_protein_hsps_by_unit`).
//! - While a stream is open (`x_frames_pooled`), the frames of 1 MiB or more of a dropped subject
//!   plan go to `X_BYTES` and the next frames and long reads take them with their pages resident
//!   (`x_subject_translate::x_frame_buffer`, `x_reuse_frame`), as BLASTN's database reads do. The
//!   frames go idle when their subject is searched, after the peak of the preliminary search, so
//!   the composition redo and the reports read their sequences into them, and at one thread each
//!   chromosome's frames are those of the one before (C. elegans genome: 6 chromosomes of 13.8 to
//!   20.9 Mb, frames of 4.6 to 7 MB). A frame takes an idle buffer that touched at most twice its
//!   bytes, since the pool does not count the pages of a lent buffer. No store is given to the
//!   pool: the stores read for the plans of 8 and 16 threads, made at once, would stay idle across
//!   the peak. When the last stream of the process closes, the idle buffers are freed.
//!
//! Every read packs the stored bytes as the copy did (`common::database::x_read_nucleotide_store`:
//! the packing of `PackedNucleotide::from_stored`), so every stage computes the same values on the
//! same subjects in the same order. A read that fails during the
//! search (a database changed or damaged after the check, an error of the file system) is an
//! `XReadError` holding the error that the check gives for the same read: the search of the
//! batches is then made again one batch at a time, so that the batch that reads it fails after the
//! reports of the batches before, as with the check's error (`tblastn_app.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:374-376 (the sequence of
//! the search is the mapped volume's bytes, not a copy)
//! ```c++
//!     else {
//!         len = seqdb.GetSequence(oid, &buf);
//!     }
//! ```
//! and c++/src/algo/blast/core/blast_kappa.c:1380-1386 (the redo fetches each match's subject)
//! ```c
//!         if( program_number == eBlastTypeTblastn ) {
//!             seq_info->seq_arg.encoding = eBlastEncodingNcbi4na;
//!         } else {
//!             seq_info->seq_arg.encoding = eBlastEncodingProtein;
//!         }
//!         if (BlastSeqSrcGetSequence(seqSrc, &seq_info->seq_arg) >= 0) {
//! ```
//!
//! It acts for every TBLASTN search with `-db`. The shadow mode (`LOSAT_X_TBNDBSTREAMSHADOW=1`,
//! which compared the check and every read with the copy) went with the copy.

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock, Weak};

use anyhow::Result;

use super::database::PROGRAM;
use crate::algorithm::common::database as shared;
use crate::seqdb::SeqDb;
use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};

/// A read of the database that failed during the search, with
/// the error that the check (`XDbStream::new`) gives for the same read.
#[derive(Debug)]
pub(super) struct XReadError(pub(super) anyhow::Error);

impl std::fmt::Display for XReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for XReadError {}

/// The sequences of `XDbStream::read_shared`.
#[derive(Default)]
struct XShared {
    /// The stores that something holds, by ordinal id.
    held: HashMap<usize, Weak<PackedNucleotide>>,
    /// The sequences that a thread reads now.
    reading: HashSet<usize>,
}

/// The read of a sequence by `XDbStream::read_shared`; dropped
/// (the read made, failed or unwound), the sequence is no longer read and the waiters wake.
struct XReading<'s, 'a> {
    stream: &'s XDbStream<'a>,
    oid: usize,
}

impl Drop for XReading<'_, '_> {
    fn drop(&mut self) {
        let mut shared = self.stream.shared.lock().unwrap_or_else(|e| e.into_inner());
        shared.reading.remove(&self.oid);
        drop(shared);
        self.stream.read_done.notify_all();
    }
}

/// The database of a `-db` search, read where it is searched.
pub(super) struct XDbStream<'a> {
    db: &'a SeqDb,
    /// A number of its own in the process (`XPlanBlock::stream`).
    id: u64,
    /// The length of every sequence, by ordinal id (`x_check_nucleotide_stores`).
    lengths: Vec<usize>,
    /// Whether a sequence is read into a buffer that glibc maps (`x_read_nucleotide_store`): with
    /// more than one thread. The freed stores of the workers stayed resident in their arenas; at
    /// one thread the heap's freed store serves the next read without faults (C. elegans genome:
    /// 19 reads of its 7 sequences, against one read of the copy).
    mapped: bool,
    /// The sequences that something holds or a thread reads (`read_shared`).
    shared: Mutex<XShared>,
    /// Woken when a read of `read_shared` ends.
    read_done: Condvar,
}

impl<'a> XDbStream<'a> {
    /// Checks `db` as the copy of the database was read (`common::database::x_check_nucleotide_stores`)
    /// and keeps the lengths.
    pub(super) fn new(db: &'a SeqDb, threads: usize) -> Result<Self> {
        static STREAMS: AtomicU64 = AtomicU64::new(0);
        let lengths = shared::x_check_nucleotide_stores(db, PROGRAM)?;
        X_OPEN_STREAMS.fetch_add(1, Ordering::Relaxed);
        Ok(Self {
            db,
            id: STREAMS.fetch_add(1, Ordering::Relaxed),
            lengths,
            mapped: threads > 1,
            shared: Mutex::new(XShared::default()),
            read_done: Condvar::new(),
        })
    }

    /// The number of sequences.
    pub(super) fn count(&self) -> usize {
        self.lengths.len()
    }

    /// The number of bases of sequence `oid`.
    pub(super) fn length(&self, oid: usize) -> usize {
        self.lengths[oid]
    }

    /// Sequence `oid`, packed as the copy of the database packed it, with one read (`XReadError` if it fails).
    pub(super) fn read(&self, oid: usize) -> Result<PackedNucleotide> {
        let ordinal = u32::try_from(oid).expect("a database ordinal id fits Uint4");
        shared::x_read_nucleotide_store(self.db, ordinal, PROGRAM, self.mapped)
            .map_err(|error| XReadError(error).into())
    }

    /// Sequence `oid` (`read`), shared with whatever holds it already: a stage or a query batch
    /// that runs at the same time reads it once, and a thread that asks for it while another
    /// reads it waits for that read (the threads of a stage that first align the same subject at
    /// once made one copy each). It is read again once nothing holds it.
    pub(super) fn read_shared(&self, oid: usize) -> Result<Arc<PackedNucleotide>> {
        let mut shared = self.shared.lock().unwrap_or_else(|e| e.into_inner());
        loop {
            if let Some(store) = shared.held.get(&oid).and_then(Weak::upgrade) {
                return Ok(store);
            }
            if !shared.reading.contains(&oid) {
                break;
            }
            shared = self
                .read_done
                .wait(shared)
                .unwrap_or_else(|e| e.into_inner());
        }
        shared.reading.insert(oid);
        drop(shared);
        let reading = XReading { stream: self, oid };
        let store = Arc::new(self.read(oid)?);
        let mut shared = self.shared.lock().unwrap_or_else(|e| e.into_inner());
        if shared.held.len() >= 1024 && shared.held.len().is_power_of_two() {
            shared.held.retain(|_, store| store.strong_count() > 0);
        }
        shared.held.insert(oid, Arc::downgrade(&store));
        drop(shared);
        drop(reading);
        Ok(store)
    }

    /// `make` of the packed bytes and the length of sequence `oid`, for the plan of its
    /// preliminary stage. The thread keeps the last block of sequences it read
    /// (`SeqDb::next_block` of `PLAN_BLOCK_BYTES` from the sequence asked for) and gives the next
    /// sequences of the block from its bytes, with the errors of their stores but without the
    /// stores: the preliminary stage asks for the subjects of a range in ascending order on each
    /// thread, so most short subjects come without a read or a copy. A block of one sequence (a
    /// long one) is read alone and shared (`read_shared`).
    pub(super) fn plan_with<R>(
        &self,
        oid: usize,
        make: impl FnOnce(&[u8], usize) -> Result<R>,
    ) -> Result<R> {
        let ordinal = u32::try_from(oid).expect("a database ordinal id fits Uint4");
        let mut kept = X_PLAN_BLOCK.with(std::cell::Cell::take);
        let held = |kept: &Option<XPlanBlock>| {
            kept.as_ref().is_some_and(|kept| {
                kept.stream == self.id
                    && (kept.first..kept.first + kept.block.len() as u32).contains(&ordinal)
            })
        };
        if !held(&kept) {
            // PROBE (LOSAT_PROBE_PLANCUT bit 2): the last block's buffer takes the next block.
            let bytes = if x_probe_plan() & 4 != 0 {
                let bytes = kept.take().map(|kept| kept.block.x_into_data());
                Some(bytes.unwrap_or_else(|| Vec::with_capacity(PLAN_BLOCK_BYTES as usize + 64)))
            } else {
                kept = None;
                None
            };
            match self.db.next_block(ordinal, PLAN_BLOCK_BYTES) {
                Some(mut range) if range.len() > 1 => {
                    // PROBE (LOSAT_PROBE_PLANCUT bit 1): the block stops at the end of the chunk.
                    if x_probe_plan() & 2 != 0 {
                        let end = X_PLAN_END.with(std::cell::Cell::get);
                        if end > ordinal && end < range.end_ordinal() {
                            range.end = range.first + (end - ordinal);
                        }
                    }
                    let block = match bytes {
                        Some(bytes) => self.db.x_read_block_into(range, bytes),
                        None => self.db.read_block(range),
                    }
                    .map_err(|error| XReadError(shared::search_read_error(error, PROGRAM)))?;
                    kept = Some(XPlanBlock {
                        stream: self.id,
                        first: range.first_ordinal,
                        block,
                    });
                }
                _ => {
                    let store = self.read_shared(oid)?;
                    return make(store.packed(), store.len());
                }
            }
        }
        let block = kept.as_ref().expect("a block holds the sequence");
        let i = (ordinal - block.first) as usize;
        let made = self
            .length_in(&block.block, i, oid)
            .and_then(|length| make(block.block.nucleotide_packed(i), length));
        X_PLAN_BLOCK.with(|cell| cell.set(kept));
        made
    }

    /// The length of sequence `oid`, sequence `i` of `block`, with the errors of its store
    /// (`XReadError`).
    fn length_in(
        &self,
        block: &crate::seqdb::SequenceBlock,
        i: usize,
        oid: usize,
    ) -> Result<usize> {
        shared::x_nucleotide_length_in(block, i, oid as u32, PROGRAM)
            .map_err(|error| XReadError(error).into())
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

/// Whether a `-db` stream is open in the process: then the large
/// frames of a dropped subject plan go to `utils::x_mapped::X_BYTES`, from which the next frames
/// and reads take them (module header). A `-subject` search that runs at the same time (library
/// use) pools its frames too; a frame taken from the pool is zero-filled before it is written, so
/// no report changes. No NCBI counterpart: only where the buffers go.
pub(super) fn x_frames_pooled() -> bool {
    X_OPEN_STREAMS.load(Ordering::Relaxed) > 0
}

/// No NCBI counterpart: the bytes of the sequence file in a block that a thread keeps for the
/// plans (`XDbStream::plan_with`; the blocks of BLASTX's `-db` stream are 256 KiB too).
const PLAN_BLOCK_BYTES: u64 = 256 << 10;

/// A block of sequences that a thread keeps for the plans, and its stream.
struct XPlanBlock {
    /// `XDbStream::id` of the stream that read it.
    stream: u64,
    first: u32,
    block: crate::seqdb::SequenceBlock,
}

thread_local! {
    /// This thread's block for the plans.
    static X_PLAN_BLOCK: std::cell::Cell<Option<XPlanBlock>> = const { std::cell::Cell::new(None) };
    /// PROBE: the ordinal id after the last subject of this thread's chunk.
    static X_PLAN_END: std::cell::Cell<u32> = const { std::cell::Cell::new(u32::MAX) };
}

/// PROBE: `LOSAT_PROBE_PLANCUT` bits: 1 = the block dropped at the end of the chunk, 2 = the block
/// cut at the end of the chunk, 4 = the last block's buffer reused.
fn x_probe_plan() -> u8 {
    static V: OnceLock<u8> = OnceLock::new();
    *V.get_or_init(|| {
        std::env::var("LOSAT_PROBE_PLANCUT")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    })
}

/// PROBE: this thread searches a chunk of subjects that ends before `end`, until it is dropped.
pub(super) struct XPlanChunk(u32);

pub(super) fn x_plan_chunk(end: usize) -> XPlanChunk {
    let end = u32::try_from(end).unwrap_or(u32::MAX);
    XPlanChunk(X_PLAN_END.with(|cell| cell.replace(end)))
}

impl Drop for XPlanChunk {
    fn drop(&mut self) {
        X_PLAN_END.with(|cell| cell.set(self.0));
        if x_probe_plan() & 1 != 0 {
            let _ = X_PLAN_BLOCK.try_with(std::cell::Cell::take);
        }
    }
}

// No NCBI counterpart: reads the LOSAT_X_TBNSTAGEDROP variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSTAGEDROP): a switch of the strict set (`utils/xswitch.rs`): a stage's
/// subject is let go after its last match (`XStageSubjects::done`) instead of at the stage's end.
pub(super) fn x_tbn_stage_drop() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNSTAGEDROP"))
}

/// The subjects that a stage after the preliminary search reads
/// (the composition redo, the traceback of `-comp_based_stats 0`): with a `-db` stream, each
/// subject of the stage's lists read when the stage first aligns it (`XDbStream::read_shared`)
/// and held until the stage has passed its last match (LOSAT_X_TBNSTAGEDROP, `done`) or, with
/// `LOSAT_X_TBNSTAGEDROP=0`, until the stage ends; otherwise the search's subjects.
pub(super) enum XStageSubjects<'a> {
    All(super::x_subject_store::SubjectLetters<'a>),
    Read {
        stream: &'a XDbStream<'a>,
        /// The subjects of the lists, ascending.
        oids: Vec<usize>,
        /// The store of each, once read, until its last match is done.
        stores: Vec<Mutex<Option<Arc<PackedNucleotide>>>>,
        /// The matches of each that are not done (LOSAT_X_TBNSTAGEDROP); empty with the switch
        /// off, where the stores stay until the stage ends.
        left: Vec<AtomicUsize>,
    },
}

/// A subject of a stage, held while a match reads it (`XStageSubjects::held`).
pub(super) enum XHeld<'a> {
    Nt(SubjectNt<'a>),
    Store(Arc<PackedNucleotide>),
}

impl XHeld<'_> {
    /// The subject for the readers of windows.
    pub(super) fn nt(&self) -> SubjectNt<'_> {
        match self {
            Self::Nt(nt) => *nt,
            Self::Store(store) => SubjectNt::Packed(store),
        }
    }
}

impl<'a> XStageSubjects<'a> {
    /// The subjects of a stage whose lists name the subjects `oids` (one item per match).
    pub(super) fn new(
        subjects: super::x_subject_store::SubjectLetters<'a>,
        oids: impl Iterator<Item = usize>,
    ) -> Self {
        let super::x_subject_store::SubjectLetters::Stream(stream) = subjects else {
            return Self::All(subjects);
        };
        let mut all: Vec<usize> = oids.collect();
        all.sort_unstable();
        let mut oids: Vec<usize> = Vec::new();
        let mut counts: Vec<usize> = Vec::new();
        for oid in all {
            if oids.last() == Some(&oid) {
                *counts.last_mut().expect("a count per subject") += 1;
            } else {
                oids.push(oid);
                counts.push(1);
            }
        }
        let stores = oids.iter().map(|_| Mutex::new(None)).collect();
        let left = if x_tbn_stage_drop() {
            counts.into_iter().map(AtomicUsize::new).collect()
        } else {
            Vec::new()
        };
        Self::Read {
            stream,
            oids,
            stores,
            left,
        }
    }

    /// Subject `oid` for the readers of windows (`SubjectLetters::nt`), read the first time.
    pub(super) fn held(&self, oid: usize) -> Result<XHeld<'_>> {
        match self {
            Self::All(subjects) => Ok(XHeld::Nt(subjects.nt(oid))),
            Self::Read {
                stream,
                oids,
                stores,
                ..
            } => {
                let at = oids
                    .binary_search(&oid)
                    .expect("every subject of the stage's lists is listed");
                let mut slot = stores[at].lock().unwrap_or_else(|e| e.into_inner());
                if let Some(store) = slot.as_ref() {
                    return Ok(XHeld::Store(Arc::clone(store)));
                }
                // Two threads may ask at once: the stream reads it once (`read_shared`).
                drop(slot);
                let store = stream.read_shared(oid)?;
                slot = stores[at].lock().unwrap_or_else(|e| e.into_inner());
                Ok(XHeld::Store(Arc::clone(slot.get_or_insert(store))))
            }
        }
    }

    /// A guard of one match of subject `oid` that calls `done` when it is dropped. Declared before
    /// the match's `held`, it runs after the match let the subject go.
    pub(super) fn x_match(&self, oid: usize) -> XMatchDone<'_, 'a> {
        XMatchDone { stage: self, oid }
    }

    /// One match of subject `oid` is done (redone, or passed by the early termination): after its
    /// last one the stage lets the subject go (LOSAT_X_TBNSTAGEDROP); a match that reads it later
    /// would read it again. No NCBI counterpart: NCBI's redo fetches each match's subject and
    /// releases it after the match (`blast_kappa.c` `s_MatchingSequenceRelease`).
    pub(super) fn done(&self, oid: usize) {
        let Self::Read {
            oids, stores, left, ..
        } = self
        else {
            return;
        };
        if left.is_empty() {
            return;
        }
        let at = oids
            .binary_search(&oid)
            .expect("every subject of the stage's lists is listed");
        let before = left[at].fetch_sub(1, Ordering::AcqRel);
        debug_assert!(
            before > 0,
            "more matches done than listed for subject {oid}"
        );
        if before == 1 {
            let store = stores[at].lock().unwrap_or_else(|e| e.into_inner()).take();
            drop(store);
        }
    }
}

/// Calls `XStageSubjects::done` for one match when it is dropped (`XStageSubjects::x_match`).
pub(super) struct XMatchDone<'s, 'a> {
    stage: &'s XStageSubjects<'a>,
    oid: usize,
}

impl Drop for XMatchDone<'_, '_> {
    fn drop(&mut self) {
        self.stage.done(self.oid);
    }
}
