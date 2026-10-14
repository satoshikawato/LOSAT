//! The subjects of a TBLASTN search as Stage D and the report read them: the reader's letters,
//! or (EXPERIMENT, `LOSAT_X_TBNSTORE`) packed subjects whose letters are made again when a
//! subject is read (followups-20261014 DESIGN_subject_store.md).

use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::sequence::packed_subject::PackedNucleotide;

// No NCBI counterpart: reads the LOSAT_X_TBNSTORE variable once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNSTORE): off unless the variable is set (it is not in the strict set).
pub(super) fn x_tbn_store() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_TBNSTORE").is_some())
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
                let mut buffer = X_LETTERS.with(Cell::take);
                LETTERS_MADE.fetch_add(1, Ordering::Relaxed);
                subjects[oid].ascii_into(&mut buffer);
                let result = f(&buffer);
                if buffer.capacity() <= X_KEPT_BUFFER {
                    X_LETTERS.with(|cell| cell.set(buffer));
                }
                result
            }
        }
    }
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
/// is returned: the search then reads the letters as without the switch.
pub(super) fn x_read_subject_stores<R: std::io::Read>(
    source: &mut crate::blastinput::fasta_reader::FastaInputSource<R>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> anyhow::Result<(
    Vec<crate::blastinput::fasta_reader::FastaRecord>,
    Option<Vec<PackedNucleotide>>,
)> {
    let read =
        crate::blastinput::fasta_reader::read_subjects_map(source, None, warn, |mut record| {
            let store = PackedNucleotide::from_subject_letters(&record.sequence);
            if store.is_some() {
                record.sequence = Vec::new();
            }
            Ok((record, store))
        })?;
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

// No NCBI counterpart: counters for LOSAT_X_STATS that tell whether the packed subjects were
// used; they do not change any value NCBI computes.
/// Subjects read packed, and the bytes their stores hold.
static STORED_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static STORED_BYTES: AtomicU64 = AtomicU64::new(0);
/// The times the letters of a packed subject were made.
static LETTERS_MADE: AtomicU64 = AtomicU64::new(0);

/// The summary line (printed by `main` at exit) with LOSAT_X_STATS: the subjects held packed,
/// the bytes of their stores, and the times their letters were made.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] TBNSTORE subjects={} stored_bytes={} letters_made={}",
            STORED_SUBJECTS.load(Ordering::Relaxed),
            STORED_BYTES.load(Ordering::Relaxed),
            LETTERS_MADE.load(Ordering::Relaxed)
        );
    }
}
