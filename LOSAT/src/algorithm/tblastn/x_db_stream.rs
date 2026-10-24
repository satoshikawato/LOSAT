//! EXPERIMENT (LOSAT_X_TBNDBSTREAM / LOSAT_X_TBNDBSTREAMSHADOW): a `-db` search of TBLASTN reads
//! each database sequence where it searches it, instead of a copy of the whole database made
//! before the first query batch.
//!
//! Without the switch, `TblastnDatabase::read_stores` packs every sequence of the database
//! (`common::database::read_nucleotide_stores`) into anonymous memory before the first query
//! batch, and every stage reads the copy until the reports: 23.9 MiB for the 100 Mb of the
//! C. elegans genome, 32.2 MiB for a zebrafish cDNA set, alive at the peak of every thread count
//! (v4-D12 s8 probe, workstation), where NCBI maps the sequence file. The reports then cloned the
//! packed sequences they show (`database::reported_subjects`).
//!
//! With the switch on:
//!
//! - `XDbStream::new` reads the database once in the blocks of `read_nucleotide_stores`, with its
//!   checks in its order (`common::database::x_check_nucleotide_stores`), and keeps the length of
//!   every sequence. Its error is the error of `read_stores`, which fails the same query batch
//!   with the same text.
//! - The preliminary stage of a subject reads its sequence when it makes the subject's plan
//!   (`search_seed::x_scan_plan`, `XPrelimSubject::Stream`) and drops it once the plan is made:
//!   the frames are translated from it and the stage reads nothing else of the sequence. A plan
//!   that the search keeps between query batches (`XSubjectSide`) is not made again, and its
//!   sequence is not read again. Each thread keeps the last block of 256 KiB of the sequence file
//!   it read for the plans (`XDbStream::plan_with`), from whose bytes the plans of the next short
//!   subjects are made, without a store.
//! - The composition redo and the traceback of `-comp_based_stats 0` read a subject of their
//!   lists when they first align it and hold it until the stage ends (`XStageSubjects`); the
//!   reports read the sequences they show (`database::reported_subjects`), one read each.
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
//!   `search_gapped::x_preliminary_protein_hsps_by_unit`); so is every sequence with
//!   LOSAT_X_TBNUNITXLATESHADOW, whose units are compared with the frames.
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
//! Every read packs the stored bytes as the copy does (`common::database::x_nucleotide_store_in`,
//! `x_read_nucleotide_store`: the packing of `PackedNucleotide::from_stored`), so every stage
//! computes the same values on the same subjects in the same order. A read that fails during the
//! search (a database changed or damaged after the check, an error of the file system) is an
//! `XReadError` holding the error that the check gives for the same read: the search of the
//! batches is then made again one batch at a time, so that the batch that reads it fails after the
//! reports of the batches before, as with the check's error (`args.rs`).
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
//! `LOSAT_X_TBNDBSTREAM` is a switch of the strict set: on unless it is `0`. It acts for every
//! TBLASTN search with `-db`.
//!
//! `LOSAT_X_TBNDBSTREAMSHADOW=1` (shadow mode, with the switch on): `read_stores` makes the copy
//! too; the check's lengths and error, and every read sequence, are compared with it. The first
//! difference stops the run (panic) after printing the line that is otherwise printed at exit:
//! `[X_SHADOW] TBNDBSTREAM sequences= reads= bases= identical`.

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock, Weak};

use anyhow::Result;

use super::database::PROGRAM;
use crate::algorithm::common::database as shared;
use crate::seqdb::SeqDb;
use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};

// No NCBI counterpart: reads the LOSAT_X_TBNDBSTREAM variable once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): a switch of the strict set (`utils/xswitch.rs`).
pub(super) fn x_tbn_db_stream() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNDBSTREAM"))
}

// No NCBI counterpart: reads the LOSAT_X_TBNDBSTREAMSHADOW variable once; it does not change any
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNDBSTREAMSHADOW): shadow mode, with `LOSAT_X_TBNDBSTREAM` on.
pub(super) fn x_shadow() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        x_tbn_db_stream()
            && std::env::var("LOSAT_X_TBNDBSTREAMSHADOW").is_ok_and(|value| value == "1")
    })
}

static SHADOW_RAN: AtomicBool = AtomicBool::new(false);
static SHADOW_SEQUENCES: AtomicU64 = AtomicU64::new(0);
static SHADOW_READS: AtomicU64 = AtomicU64::new(0);
static SHADOW_BASES: AtomicU64 = AtomicU64::new(0);

/// The shadow line, printed at exit (`x_print_shadow_summary`) or before the panic of a difference.
fn shadow_line() -> String {
    format!(
        "[X_SHADOW] TBNDBSTREAM sequences={} reads={} bases={} identical",
        SHADOW_SEQUENCES.load(Ordering::Relaxed),
        SHADOW_READS.load(Ordering::Relaxed),
        SHADOW_BASES.load(Ordering::Relaxed),
    )
}

/// Stops the run at a difference of the shadow mode.
fn shadow_assert(same: bool, what: impl FnOnce() -> String) {
    if !same {
        eprintln!("{}", shadow_line());
        panic!("LOSAT_X_TBNDBSTREAMSHADOW: {}", what());
    }
}

/// Prints the shadow line once, if the shadow mode compared anything in this run.
pub(super) fn x_print_shadow_summary() {
    if SHADOW_RAN.swap(false, Ordering::Relaxed) {
        eprintln!("{}", shadow_line());
    }
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): a read of the database that failed during the search, with
/// the error that the check (`XDbStream::new`) gives for the same read.
#[derive(Debug)]
pub(super) struct XReadError(pub(super) anyhow::Error);

impl std::fmt::Display for XReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for XReadError {}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the sequences of `XDbStream::read_shared`.
#[derive(Default)]
struct XShared {
    /// The stores that something holds, by ordinal id.
    held: HashMap<usize, Weak<PackedNucleotide>>,
    /// The sequences that a thread reads now.
    reading: HashSet<usize>,
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the read of a sequence by `XDbStream::read_shared`; dropped
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

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the database of a `-db` search, read where it is searched.
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
    /// The copy of the shadow mode.
    copy: Option<Vec<PackedNucleotide>>,
}

impl<'a> XDbStream<'a> {
    /// Checks `db` as `read_nucleotide_stores` reads it and keeps the lengths; in the shadow mode
    /// also makes the copy (`read_stores`) and compares the check's result with it.
    pub(super) fn new(
        db: &'a SeqDb,
        threads: usize,
        read_stores: impl FnOnce() -> Result<Vec<PackedNucleotide>>,
    ) -> Result<Self> {
        let checked = shared::x_check_nucleotide_stores(db, PROGRAM);
        let copy = if x_shadow() {
            SHADOW_RAN.store(true, Ordering::Relaxed);
            for counter in [&SHADOW_SEQUENCES, &SHADOW_READS, &SHADOW_BASES] {
                counter.store(0, Ordering::Relaxed);
            }
            let copy = read_stores();
            match (&checked, &copy) {
                (Ok(lengths), Ok(copy)) => {
                    shadow_assert(
                        lengths.len() == copy.len()
                            && lengths
                                .iter()
                                .zip(copy)
                                .all(|(&length, store)| length == store.len()),
                        || "the lengths of the check differ from the copy's".to_string(),
                    );
                    SHADOW_SEQUENCES.store(lengths.len() as u64, Ordering::Relaxed);
                }
                (Err(check), Err(copy)) => {
                    shadow_assert(format!("{check:#}") == format!("{copy:#}"), || {
                        format!("the check's error ({check:#}) differs from the copy's ({copy:#})")
                    })
                }
                _ => shadow_assert(false, || {
                    "the check and the copy differ in success".to_string()
                }),
            }
            copy.ok()
        } else {
            None
        };
        static STREAMS: AtomicU64 = AtomicU64::new(0);
        let lengths = checked?;
        X_OPEN_STREAMS.fetch_add(1, Ordering::Relaxed);
        Ok(Self {
            db,
            id: STREAMS.fetch_add(1, Ordering::Relaxed),
            lengths,
            mapped: threads > 1,
            shared: Mutex::new(XShared::default()),
            read_done: Condvar::new(),
            copy,
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

    /// Sequence `oid`, packed as the copy packs it, with one read (`XReadError` if it fails).
    pub(super) fn read(&self, oid: usize) -> Result<PackedNucleotide> {
        let ordinal = u32::try_from(oid).expect("a database ordinal id fits Uint4");
        let store = shared::x_read_nucleotide_store(self.db, ordinal, PROGRAM, self.mapped)
            .map_err(XReadError)?;
        Ok(self.checked(oid, store))
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
            kept = None;
            match self.db.next_block(ordinal, PLAN_BLOCK_BYTES) {
                Some(range) if range.len() > 1 => {
                    let block = self
                        .db
                        .read_block(range)
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
    /// (`XReadError`); the shadow mode packs the store and compares it with the copy.
    fn length_in(
        &self,
        block: &crate::seqdb::SequenceBlock,
        i: usize,
        oid: usize,
    ) -> Result<usize> {
        let ordinal = oid as u32;
        if self.copy.is_some() {
            let store =
                shared::x_nucleotide_store_in(block, i, ordinal, PROGRAM).map_err(XReadError)?;
            return Ok(self.checked(oid, store).len());
        }
        shared::x_nucleotide_length_in(block, i, ordinal, PROGRAM)
            .map_err(|error| XReadError(error).into())
    }

    /// `store`, read for sequence `oid`, after the shadow mode's comparison with the copy.
    fn checked(&self, oid: usize, store: PackedNucleotide) -> PackedNucleotide {
        if let Some(copy) = &self.copy {
            shadow_assert(store == copy[oid], || {
                format!("database sequence {oid} differs from the copy")
            });
            SHADOW_READS.fetch_add(1, Ordering::Relaxed);
            SHADOW_BASES.fetch_add(store.len() as u64, Ordering::Relaxed);
        }
        store
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

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): whether a `-db` stream is open in the process: then the large
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
    /// EXPERIMENT (LOSAT_X_TBNDBSTREAM): this thread's block for the plans.
    static X_PLAN_BLOCK: std::cell::Cell<Option<XPlanBlock>> = const { std::cell::Cell::new(None) };
}

/// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the subjects that a stage after the preliminary search reads
/// (the composition redo, the traceback of `-comp_based_stats 0`): with a `-db` stream, each
/// subject of the stage's lists read when the stage first aligns it (`XDbStream::read_shared`)
/// and held until the stage ends; otherwise the search's subjects.
pub(super) enum XStageSubjects<'a> {
    All(super::x_subject_store::SubjectLetters<'a>),
    Read {
        stream: &'a XDbStream<'a>,
        /// The subjects of the lists, ascending.
        oids: Vec<usize>,
        /// The store of each, once read.
        stores: Vec<OnceLock<Arc<PackedNucleotide>>>,
    },
}

impl<'a> XStageSubjects<'a> {
    /// The subjects of a stage whose lists name the subjects `oids`.
    pub(super) fn new(
        subjects: super::x_subject_store::SubjectLetters<'a>,
        oids: impl Iterator<Item = usize>,
    ) -> Self {
        let super::x_subject_store::SubjectLetters::Stream(stream) = subjects else {
            return Self::All(subjects);
        };
        let mut oids: Vec<usize> = oids.collect();
        oids.sort_unstable();
        oids.dedup();
        let stores = oids.iter().map(|_| OnceLock::new()).collect();
        Self::Read {
            stream,
            oids,
            stores,
        }
    }

    /// Subject `oid` for the readers of windows (`SubjectLetters::nt`), read the first time.
    pub(super) fn nt(&self, oid: usize) -> Result<SubjectNt<'_>> {
        match self {
            Self::All(subjects) => Ok(subjects.nt(oid)),
            Self::Read {
                stream,
                oids,
                stores,
            } => {
                let at = oids
                    .binary_search(&oid)
                    .expect("every subject of the stage's lists is listed");
                if stores[at].get().is_none() {
                    // Two threads may read it at once: one store is kept.
                    let _ = stores[at].set(stream.read_shared(oid)?);
                }
                Ok(SubjectNt::Packed(
                    stores[at].get().expect("the store was set"),
                ))
            }
        }
    }
}
