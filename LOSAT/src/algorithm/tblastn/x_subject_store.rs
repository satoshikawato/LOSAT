//! The subjects of a TBLASTN search as Stage D and the report read them: the reader's letters,
//! or (EXPERIMENT, `LOSAT_X_TBNSTORE`) packed subjects whose letters are made again when a
//! subject is read (followups-20261014 DESIGN_subject_store.md).

use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};

// No NCBI counterpart: reads the LOSAT_X_TBNSTORE variable once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNSTORE): a switch of the strict set (`utils/xswitch.rs`): on unless it
/// is `0`, or `LOSAT_X_DEFAULT=0` and it is not set.
pub(super) fn x_tbn_store() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNSTORE"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
// ```c
//     virtual SBlastSequence GetCompressedPlusStrand() {
//         SBlastSequence retval(size());
//         string ncbi4na = kEmptyStr;
//         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
//         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
//         return retval;
// ```
// NCBI keeps the subjects in its scope and makes each one's search forms from it. Both variants
// give the same letters for every subject; only where they are held differs.
/// The subjects of one search, by OID.
#[derive(Clone, Copy)]
pub(super) enum SubjectLetters<'a> {
    /// The letters of the reader's records.
    Records(&'a [&'a [u8]]),
    /// EXPERIMENT (LOSAT_X_TBNSTORE): packed subjects, the letters of one made when it is read.
    Packed(&'a [PackedNucleotide]),
}

impl<'a> SubjectLetters<'a> {
    /// The number of subjects.
    pub(super) fn count(&self) -> usize {
        match self {
            Self::Records(subjects) => subjects.len(),
            Self::Packed(subjects) => subjects.len(),
        }
    }

    /// The number of letters of subject `oid`.
    pub(super) fn length(&self, oid: usize) -> usize {
        match self {
            Self::Records(subjects) => subjects[oid].len(),
            Self::Packed(subjects) => subjects[oid].len(),
        }
    }

    /// The numbers of letters of the subjects, in OID order.
    pub(super) fn lengths(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.count()).map(|oid| self.length(oid))
    }

    /// The letters of subject `oid`: the record's own, or made in `buffer` from its store.
    pub(super) fn letters<'b>(&self, oid: usize, buffer: &'b mut Vec<u8>) -> &'b [u8]
    where
        'a: 'b,
    {
        match self {
            Self::Records(subjects) => subjects[oid],
            Self::Packed(subjects) => {
                LETTERS_MADE.fetch_add(1, Ordering::Relaxed);
                subjects[oid].ascii_into(buffer);
                buffer
            }
        }
    }

    /// `f` of the letters of subject `oid`, made in a buffer that the thread keeps between
    /// subjects (a buffer above `X_KEPT_BUFFER` is freed after use, so that a thread does not keep
    /// a chromosome's letters). Re-entrant: a nested call (a pool thread that takes another
    /// subject while it waits) gets a buffer of its own.
    pub(super) fn with_letters<R>(&self, oid: usize, f: impl FnOnce(&[u8]) -> R) -> R {
        match self {
            Self::Records(subjects) => f(subjects[oid]),
            Self::Packed(subjects) => {
                LETTERS_MADE.fetch_add(1, Ordering::Relaxed);
                with_buffer(|buffer| subjects[oid].ascii_into(buffer), f)
            }
        }
    }

    /// `f` of subject `oid` as its preliminary stage reads it. With `case` (the stage reads the
    /// lower-case letters, `-lcase_masking`) it is the record's letters (`with_letters`).
    /// Without it a packed subject is given as its store: the stage translates its packed bases
    /// (the bases the letters resolve to), so a worker holds the frames alone and neither the
    /// letters nor the resolved bases (`x_with_search_letters` makes the bases where a stage
    /// reads letters).
    pub(super) fn with_search_subject<R>(
        &self,
        oid: usize,
        case: bool,
        f: impl FnOnce(SubjectNt<'_>) -> R,
    ) -> R {
        match self {
            Self::Packed(subjects) if !case => f(SubjectNt::Packed(&subjects[oid])),
            _ => self.with_letters(oid, |letters| f(SubjectNt::Letters(letters))),
        }
    }

    /// Subject `oid` for the readers of windows (the redo, the traceback's translations, the
    /// report): the record's letters, or its store, whose windows are made when read.
    pub(super) fn nt(&self, oid: usize) -> SubjectNt<'a> {
        match self {
            Self::Records(subjects) => SubjectNt::Letters(subjects[oid]),
            Self::Packed(subjects) => SubjectNt::Packed(&subjects[oid]),
        }
    }
}

/// `f` of the letters of a subject as a preliminary stage that reads letters takes them: the
/// letters, or a packed subject's bases as resolved, `A`, `C`, `G` or `T` for every base. The
/// stage resolves the letters into these same bases, which resolve to themselves.
pub(super) fn x_with_search_letters<R>(subject: SubjectNt<'_>, f: impl FnOnce(&[u8]) -> R) -> R {
    match subject {
        SubjectNt::Letters(letters) => f(letters),
        SubjectNt::Packed(store) => {
            RESOLVED_MADE.fetch_add(1, Ordering::Relaxed);
            with_buffer(|buffer| store.resolved_letters_into(buffer), f)
        }
    }
}

/// `f` of the bytes `make` writes in this thread's buffer (see `with_letters`).
fn with_buffer<R>(make: impl FnOnce(&mut Vec<u8>), f: impl FnOnce(&[u8]) -> R) -> R {
    let mut buffer = X_LETTERS.with(Cell::take);
    make(&mut buffer);
    let result = f(&buffer);
    if buffer.capacity() <= X_KEPT_BUFFER {
        X_LETTERS.with(|cell| cell.set(buffer));
    }
    result
}

/// The largest letter buffer a thread keeps between subjects.
const X_KEPT_BUFFER: usize = 8 << 20;

thread_local! {
    /// EXPERIMENT (LOSAT_X_TBNSTORE): this thread's letter buffer between subjects.
    static X_LETTERS: Cell<Vec<u8>> = const { Cell::new(Vec::new()) };
}

/// EXPERIMENT (LOSAT_X_TBNSTORE): the subject records (`read_subjects` without a range), each
/// one's letters packed as soon as it has been read, so the letters of the whole set are never
/// held at once. The records keep their deflines and no letters. If a record holds a byte without
/// an ncbi4na code (only an API caller can make one), every record keeps its letters and no store
/// is returned: the search then reads the letters as without the switch. With `threads` above one
/// a second thread packs each record while the reader reads the next one (`x_packed_while_read`).
pub(super) fn x_read_subject_stores<R: std::io::Read>(
    source: &mut crate::blastinput::fasta_reader::FastaInputSource<R>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    threads: usize,
) -> anyhow::Result<(
    Vec<crate::blastinput::fasta_reader::FastaRecord>,
    Option<Vec<PackedNucleotide>>,
)> {
    #[cfg(feature = "parallel")]
    let read = if threads > 1 {
        x_packed_while_read(source, warn)?
    } else {
        crate::blastinput::fasta_reader::read_subjects_map(source, None, warn, |record| {
            Ok(x_packed(record))
        })?
    };
    #[cfg(not(feature = "parallel"))]
    let read = {
        let _ = threads;
        crate::blastinput::fasta_reader::read_subjects_map(source, None, warn, |record| {
            Ok(x_packed(record))
        })?
    };
    if read.iter().all(|(_, store)| store.is_some()) {
        let (records, stores): (_, Vec<PackedNucleotide>) = read
            .into_iter()
            .map(|(record, store)| (record, store.expect("every record packed")))
            .unzip();
        STORED_SUBJECTS.store(stores.len() as u64, Ordering::Relaxed);
        STORED_BYTES.store(
            stores.iter().map(|store| store.heap_bytes() as u64).sum(),
            Ordering::Relaxed,
        );
        return Ok((records, Some(stores)));
    }
    Ok((
        read.into_iter()
            .map(|(mut record, store)| {
                if let Some(store) = store {
                    store.ascii_into(&mut record.sequence);
                }
                record
            })
            .collect(),
        None,
    ))
}

/// A record and its store (`None` when a byte has no ncbi4na code; the record keeps its letters
/// then, and drops them otherwise).
fn x_packed(
    mut record: crate::blastinput::fasta_reader::FastaRecord,
) -> (
    crate::blastinput::fasta_reader::FastaRecord,
    Option<PackedNucleotide>,
) {
    let store = PackedNucleotide::from_subject_letters(&record.sequence);
    if store.is_some() {
        record.sequence = Vec::new();
    }
    (record, store)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
// ```c
//     virtual SBlastSequence GetCompressedPlusStrand() {
//         SBlastSequence retval(size());
//         string ncbi4na = kEmptyStr;
//         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
//         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
// ```
// No NCBI counterpart for the thread: NCBI packs a subject when its search needs it. Here the
// reader packs the records of the whole set before the search, on the thread that reads them, so
// the packing adds to the time the search waits for its subjects (about 0.5 s for the 1.37 Gb
// zebrafish genome). With two or more search threads a second thread packs the records while the
// reader reads the next ones. The reader hands them over in batches (up to `X_PACK_BATCH`
// records or `X_PACK_BATCH_LETTERS` letters, so that a set of many short records costs few
// hand-overs) and only when that thread takes a batch (a channel without a buffer), so at most
// two batches' letters are held; the records come back in reading order. Packing is a function
// of a record alone, so the records and stores are those of the reader packing them itself.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the records read and packed on two threads.
#[cfg(feature = "parallel")]
fn x_packed_while_read<R: std::io::Read>(
    source: &mut crate::blastinput::fasta_reader::FastaInputSource<R>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> anyhow::Result<
    Vec<(
        crate::blastinput::fasta_reader::FastaRecord,
        Option<PackedNucleotide>,
    )>,
> {
    std::thread::scope(|scope| {
        let (batches, received) =
            std::sync::mpsc::sync_channel::<Vec<crate::blastinput::fasta_reader::FastaRecord>>(0);
        let packer = scope.spawn(move || {
            received
                .into_iter()
                .flatten()
                .map(x_packed)
                .collect::<Vec<_>>()
        });
        let ended = || anyhow::anyhow!("TBLASTN: the subject packing thread ended");
        let mut batch = Vec::new();
        let mut letters = 0usize;
        let read =
            crate::blastinput::fasta_reader::read_subjects_map(source, None, warn, |record| {
                letters += record.sequence.len();
                batch.push(record);
                if batch.len() >= X_PACK_BATCH || letters >= X_PACK_BATCH_LETTERS {
                    letters = 0;
                    batches
                        .send(std::mem::take(&mut batch))
                        .map_err(|_| ended())?;
                }
                Ok(())
            });
        let last = if read.is_ok() && !batch.is_empty() {
            batches.send(batch).map_err(|_| ended())
        } else {
            Ok(())
        };
        // the packer ends when the sender is dropped, also after an error of the reader
        drop(batches);
        let packed = packer
            .join()
            .map_err(|_| anyhow::anyhow!("TBLASTN: the subject packing thread panicked"))?;
        read?;
        last?;
        Ok(packed)
    })
}

/// EXPERIMENT (LOSAT_X_TBNSTORE): the records and the letters of one hand-over of
/// `x_packed_while_read`.
#[cfg(feature = "parallel")]
const X_PACK_BATCH: usize = 1024;
#[cfg(feature = "parallel")]
const X_PACK_BATCH_LETTERS: usize = 4 << 20;

// No NCBI counterpart: counters for LOSAT_X_STATS that tell whether the packed subjects were
// used; they do not change any value NCBI computes.
/// Subjects read packed, and the bytes their stores hold.
static STORED_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static STORED_BYTES: AtomicU64 = AtomicU64::new(0);
/// The times the letters of a packed subject were made.
static LETTERS_MADE: AtomicU64 = AtomicU64::new(0);
/// The times the resolved bases of a packed subject were made for its preliminary stage.
static RESOLVED_MADE: AtomicU64 = AtomicU64::new(0);
/// The subject plans translated from the packed bases.
static PACKED_PLANS: AtomicU64 = AtomicU64::new(0);

/// Counts a subject plan translated from the packed bases (`x_packed_plan`).
// No NCBI counterpart: a counter for LOSAT_X_STATS; it does not change any value NCBI computes.
pub(super) fn x_count_packed_plan() {
    PACKED_PLANS.fetch_add(1, Ordering::Relaxed);
}

/// The summary line (printed by `main` at exit) with LOSAT_X_STATS: the subjects held packed,
/// the bytes of their stores, the times their letters and resolved bases were made, and the
/// plans translated from their packed bases.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] TBNSTORE subjects={} stored_bytes={} letters_made={} resolved_made={} packed_plans={}",
            STORED_SUBJECTS.load(Ordering::Relaxed),
            STORED_BYTES.load(Ordering::Relaxed),
            LETTERS_MADE.load(Ordering::Relaxed),
            RESOLVED_MADE.load(Ordering::Relaxed),
            PACKED_PLANS.load(Ordering::Relaxed)
        );
    }
}
