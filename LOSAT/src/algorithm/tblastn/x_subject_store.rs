//! The subjects of a TBLASTN search as Stage D and the report read them: the reader's letters,
//! or (EXPERIMENT, `LOSAT_X_TBNSTORE`) packed subjects whose letters are made again when a
//! subject is read (followups-20261014 DESIGN_subject_store.md).

use std::cell::Cell;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use crate::algorithm::common::x_db_stream::XDbStream;
use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};

// No NCBI counterpart: reads the LOSAT_X_TBNSTORE variable once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNSTORE): a switch of the strict set (`utils/xswitch.rs`): on unless it
/// is `0`.
pub(super) fn x_tbn_store() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNSTORE"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
// NCBI keeps the subjects in its scope and makes each one's search forms from it. Both variants
// give the same letters for every subject; only where they are held differs.
/// The subjects of one search, by OID.
#[derive(Clone, Copy)]
pub(super) enum SubjectLetters<'a> {
    /// The letters of the reader's records.
    Records(&'a [&'a [u8]]),
    /// EXPERIMENT (LOSAT_X_TBNSTORE): packed subjects, the letters of one made when it is read.
    Packed(&'a [PackedNucleotide]),
    /// The sequences of a `-db` database, each read where it is
    /// searched (`common/x_db_stream.rs`).
    Stream(&'a XDbStream<'a>),
    /// The numbers of bases alone, of the reported sequences of a `-db` search whose reports read
    /// no letters (`database::reported_subjects`).
    Lengths(&'a [usize]),
}

impl<'a> SubjectLetters<'a> {
    /// The number of subjects.
    pub(super) fn count(&self) -> usize {
        match self {
            Self::Records(subjects) => subjects.len(),
            Self::Packed(subjects) => subjects.len(),
            Self::Stream(stream) => stream.count(),
            Self::Lengths(lengths) => lengths.len(),
        }
    }

    /// The number of letters of subject `oid`. A `-db` stream holds no lengths: its search reads
    /// a subject's length where it reads the subject (`XPrelimSubject::Stream`), its reports
    /// from the reported sequences (`database::reported_subjects`).
    pub(super) fn length(&self, oid: usize) -> usize {
        match self {
            Self::Records(subjects) => subjects[oid].len(),
            Self::Packed(subjects) => subjects[oid].len(),
            Self::Stream(_) => unreachable!("a -db stream reads its lengths with its sequences"),
            Self::Lengths(lengths) => lengths[oid],
        }
    }

    /// The numbers of letters of the subjects in OID order, a `-db` stream's read from the index
    /// and the last byte of each sequence (`XDbStream::read_length`): for an observer of the
    /// search (a library caller's), not the search.
    pub(super) fn x_read_lengths(&self) -> anyhow::Result<Vec<usize>> {
        match self {
            Self::Stream(stream) => (0..stream.count())
                .map(|oid| stream.read_length(oid))
                .collect(),
            _ => Ok(self.lengths().collect()),
        }
    }

    /// The numbers of letters of the subjects, in OID order.
    pub(super) fn lengths(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.count()).map(|oid| self.length(oid))
    }

    /// The number of letters of the longest subject and of all the subjects (saturated): a `-db`
    /// database sums them up in its check (`XDbStream::longest_and_total`), as it holds no
    /// length per sequence.
    pub(super) fn longest_and_total(&self) -> (usize, u64) {
        if let Self::Stream(stream) = self {
            return stream.longest_and_total();
        }
        self.lengths().fold((0, 0), |(longest, total), length| {
            let length_u64 = u64::try_from(length).unwrap_or(u64::MAX);
            (
                usize::max(longest, length),
                total.saturating_add(length_u64),
            )
        })
    }

    /// The total of `longest_and_total`, saturated to `usize`.
    pub(super) fn total_bases(&self) -> usize {
        usize::try_from(self.longest_and_total().1).unwrap_or(usize::MAX)
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
            Self::Stream(_) | Self::Lengths(_) => {
                unreachable!("a database sequence has no lower-case letters")
            }
        }
    }

    /// `f` of subject `oid` as its preliminary stage reads it. With `case` (the stage reads the
    /// lower-case letters, `-lcase_masking`) it is the record's letters (`with_letters`) where it
    /// has lower-case letters. Otherwise a packed subject is given as its store: the stage
    /// translates its packed bases (the bases the letters resolve to), so a worker holds the
    /// frames alone and neither the letters nor the resolved bases. A database sequence (`-db`)
    /// has no lower-case letters, so its stage always reads its stored bases, as NCBI's (no
    /// subject mask without `-db_soft_mask`).
    ///
    /// A sequence of a `-db` stream is given unread and without its length
    /// (`XSearchSubject::Stream`): the search reads its length (`XPrelimSubject::Stream`), and the
    /// stage reads it when it makes its plan.
    pub(super) fn with_search_subject<R>(
        &self,
        oid: usize,
        case: bool,
        f: impl FnOnce(XSearchSubject<'_>) -> R,
    ) -> R {
        match self {
            Self::Packed(subjects) if !case || subjects[oid].lowercase().is_empty() => {
                f(XSearchSubject::Nt(SubjectNt::Packed(&subjects[oid])))
            }
            Self::Packed(_) => {
                // A packed subject whose lower-case letters the stage reads: the letters.
                crate::utils::xreach::x_reach!("LOSAT_X_TBNSTORE", false);
                self.with_letters(oid, |letters| {
                    f(XSearchSubject::Nt(SubjectNt::Letters(letters)))
                })
            }
            Self::Records(_) => self.with_letters(oid, |letters| {
                f(XSearchSubject::Nt(SubjectNt::Letters(letters)))
            }),
            Self::Stream(stream) => f(XSearchSubject::Stream(stream, oid)),
            Self::Lengths(_) => unreachable!("the lengths of the reports are not searched"),
        }
    }

    /// Subject `oid` for the readers of windows (the redo, the traceback's translations, the
    /// report): the record's letters, or its store, whose windows are made when read.
    /// The stages after the preliminary search read a `-db` stream's sequences through
    /// `x_db_stream::XStageSubjects`; the reports see the reported sequences as `Packed`.
    pub(super) fn nt(&self, oid: usize) -> SubjectNt<'a> {
        match self {
            Self::Records(subjects) => SubjectNt::Letters(subjects[oid]),
            Self::Packed(subjects) => SubjectNt::Packed(&subjects[oid]),
            Self::Stream(_) => unreachable!("a -db stream is read through XStageSubjects"),
            Self::Lengths(_) => unreachable!("reports that read no letters hold the lengths"),
        }
    }
}

/// A subject as `SubjectLetters::with_search_subject` gives it: its letters or its store, or a
/// sequence of a `-db` stream (its OID), not read and without its length.
#[derive(Clone, Copy)]
pub(super) enum XSearchSubject<'a> {
    Nt(SubjectNt<'a>),
    Stream(&'a XDbStream<'a>, usize),
}

/// The subject of one preliminary stage: its letters or its
/// store, or a sequence of a `-db` stream (its OID) with its length, which
/// `search_seed::x_scan_plan` reads when it makes the subject's plan and drops once the plan is
/// made.
#[derive(Clone, Copy)]
pub(super) enum XPrelimSubject<'a> {
    Nt(SubjectNt<'a>),
    Stream(&'a XDbStream<'a>, usize, usize),
}

impl XPrelimSubject<'_> {
    /// The number of bases.
    pub(super) fn len(&self) -> usize {
        match self {
            Self::Nt(subject) => subject.len(),
            Self::Stream(_, _, length) => *length,
        }
    }
}

impl<'a> From<SubjectNt<'a>> for XPrelimSubject<'a> {
    fn from(subject: SubjectNt<'a>) -> Self {
        Self::Nt(subject)
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

/// Frees this thread's letter buffer at the end of its loop over the chunks of the subject flow
/// (`x_oid_flow::x_loop_free`, LOSAT_X_TBNLOOPFREE).
#[cfg(feature = "parallel")]
pub(super) fn x_free_letters() {
    let _ = X_LETTERS.try_with(|buffer| drop(buffer.take()));
}

/// EXPERIMENT (LOSAT_X_TBNSTORE): the subject records, each one's letters packed as soon as it has
/// been read (`blastinput::x_packed_subjects::x_read_packed_subjects`: cut to `-subject_loc` first,
/// packed on a second thread at two or more threads natively, on the reader in WebAssembly). The
/// records keep their deflines and no letters; a set with a byte without an ncbi4na code keeps
/// every record's letters and gives no stores.
pub(super) fn x_read_subject_stores<R: std::io::Read>(
    source: &mut crate::blastinput::fasta_reader::FastaInputSource<R>,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    threads: usize,
) -> anyhow::Result<(
    Vec<crate::blastinput::fasta_reader::FastaRecord>,
    Option<Vec<PackedNucleotide>>,
    crate::blastinput::seq_range::Placements,
)> {
    let (records, stores, placements) =
        crate::blastinput::x_packed_subjects::x_read_packed_subjects(
            source, range, warn, threads, "TBLASTN",
        )?;
    // A byte without an ncbi4na code keeps the letters (the reference path).
    crate::utils::xreach::x_reach!("LOSAT_X_TBNSTORE", stores.is_some());
    if let Some(stores) = &stores {
        STORED_SUBJECTS.store(stores.len() as u64, Ordering::Relaxed);
        STORED_BYTES.store(
            stores.iter().map(|store| store.heap_bytes() as u64).sum(),
            Ordering::Relaxed,
        );
    }
    Ok((records, stores, placements))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
// NCBI makes the compressed subject of every search from its scope. A subject the search reads
// as letters (the records a caller gives `run_local`: LOSAT Web, the library; the records of a set
// that keeps its letters because one of them holds a byte without an ncbi4na code) was made into
// one byte per base for its plan; packed, it is a quarter of a byte per base, and the frames are
// read from the packed bases. The worker that makes the plan packs the subject, and the store is
// dropped once the plan is made.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the letters of a subject searched from its letters, packed for
/// its plan: `None` with the switch off, when the plan reads the case of its letters
/// (`mask_lowercase` and a lower-case letter), when a byte has no ncbi4na code, or when the
/// subject is longer than a store holds (`u32::MAX` bases).
pub(super) fn x_plan_store(letters: &[u8], mask_lowercase: bool) -> Option<PackedNucleotide> {
    if !x_tbn_store() {
        return None;
    }
    let store = if mask_lowercase && letters.iter().any(u8::is_ascii_lowercase) {
        None
    } else {
        PackedNucleotide::from_subject_letters(letters)
    };
    crate::utils::xreach::x_reach!("LOSAT_X_TBNSTORE", store.is_some());
    store
}

// No NCBI counterpart: counters for LOSAT_X_STATS that tell whether the packed subjects were
// used; they do not change any value NCBI computes.
/// Subjects read packed, and the bytes their stores hold.
static STORED_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static STORED_BYTES: AtomicU64 = AtomicU64::new(0);
/// The times the whole letters of a packed subject were made (not the windows the readers make).
static LETTERS_MADE: AtomicU64 = AtomicU64::new(0);
/// The subject plans translated from the packed bases.
static PACKED_PLANS: AtomicU64 = AtomicU64::new(0);

/// Counts a subject plan translated from the packed bases (`x_packed_plan`).
// No NCBI counterpart: a counter for LOSAT_X_STATS; it does not change any value NCBI computes.
pub(super) fn x_count_packed_plan() {
    PACKED_PLANS.fetch_add(1, Ordering::Relaxed);
}

/// The summary line (printed by `main` at exit) with LOSAT_X_STATS: the subjects held packed,
/// the bytes of their stores, the times their letters were made, and the plans translated from
/// their packed bases.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] TBNSTORE subjects={} stored_bytes={} letters_made={} packed_plans={}",
            STORED_SUBJECTS.load(Ordering::Relaxed),
            STORED_BYTES.load(Ordering::Relaxed),
            LETTERS_MADE.load(Ordering::Relaxed),
            PACKED_PLANS.load(Ordering::Relaxed)
        );
    }
}

#[cfg(test)]
#[path = "x_subject_store_tests.rs"]
mod tests;
