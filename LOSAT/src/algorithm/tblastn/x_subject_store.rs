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
    /// The sequences of a `-db` database, each read where it is
    /// searched (`x_db_stream.rs`).
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

    /// The number of letters of subject `oid`.
    pub(super) fn length(&self, oid: usize) -> usize {
        match self {
            Self::Records(subjects) => subjects[oid].len(),
            Self::Packed(subjects) => subjects[oid].len(),
            Self::Stream(stream) => stream.length(oid),
            Self::Lengths(lengths) => lengths[oid],
        }
    }

    /// The numbers of letters of the subjects, in OID order.
    pub(super) fn lengths(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.count()).map(|oid| self.length(oid))
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
    /// A sequence of a `-db` stream is given unread
    /// (`XPrelimSubject::Stream`): the stage reads it when it makes its plan.
    pub(super) fn with_search_subject<R>(
        &self,
        oid: usize,
        case: bool,
        f: impl FnOnce(XPrelimSubject<'_>) -> R,
    ) -> R {
        match self {
            Self::Packed(subjects) if !case || subjects[oid].lowercase().is_empty() => {
                f(XPrelimSubject::Nt(SubjectNt::Packed(&subjects[oid])))
            }
            Self::Packed(_) => {
                // A packed subject whose lower-case letters the stage reads: the letters.
                crate::utils::xreach::x_reach!("LOSAT_X_TBNSTORE", false);
                self.with_letters(oid, |letters| {
                    f(XPrelimSubject::Nt(SubjectNt::Letters(letters)))
                })
            }
            Self::Records(_) => self.with_letters(oid, |letters| {
                f(XPrelimSubject::Nt(SubjectNt::Letters(letters)))
            }),
            Self::Stream(stream) => f(XPrelimSubject::Stream(stream, oid)),
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

/// The subject of one preliminary stage: its letters or its
/// store, or a sequence of a `-db` stream, which `search_seed::x_scan_plan` reads when it makes
/// the subject's plan and drops once the plan is made.
#[derive(Clone, Copy)]
pub(super) enum XPrelimSubject<'a> {
    Nt(SubjectNt<'a>),
    Stream(&'a XDbStream<'a>, usize),
}

impl XPrelimSubject<'_> {
    /// The number of bases.
    pub(super) fn len(&self) -> usize {
        match self {
            Self::Nt(subject) => subject.len(),
            Self::Stream(stream, oid) => stream.length(*oid),
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
// ```c
//     virtual SBlastSequence GetCompressedPlusStrand() {
//         SBlastSequence retval(size());
//         string ncbi4na = kEmptyStr;
//         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
//         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
// ```
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
mod tests {
    use super::*;
    use crate::blastinput::fasta_reader::{read_subjects, FastaInputSource, ReaderConfig};
    use crate::blastinput::seq_range::{cut_subjects, SequenceRange};

    // No NCBI counterpart: when the letters of a subject are packed for its plan. The test holds
    // with the switch on (the default) and off (`LOSAT_X_TBNSTORE=0` in the test's environment).
    #[test]
    fn a_plan_packs_its_letters_unless_it_reads_their_case_or_cannot() {
        let on = x_tbn_store();
        let packed = |letters: &[u8], mask_lowercase: bool| {
            x_plan_store(letters, mask_lowercase).map(|store| {
                let mut back = Vec::new();
                store.ascii_into(&mut back);
                back
            })
        };
        // every IUPAC letter, `-` and lower case: the store gives the letters back
        let letters = b"ACGTacgtNRYKMSWBDHVnrykmswbdhv-ACGT".to_vec();
        assert_eq!(packed(&letters, false), on.then(|| letters.clone()));
        // `-lcase_masking` reads the case: letters with lower case stay letters
        assert_eq!(packed(&letters, true), None);
        let upper = b"ACGTNRYKMSWBDHV".to_vec();
        assert_eq!(packed(&upper, true), on.then(|| upper.clone()));
        // a byte without an ncbi4na code (only an API caller can give one)
        assert_eq!(packed(b"ACGU", false), None);
        assert_eq!(packed(b"", false), on.then(Vec::new));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-431
    // ```c
    //     static void
    //     s_Ncbi4naToNcbi2na(const string& ncbi4na, size_t base_length, Uint1* ncbi2na)
    //     {
    //         ...
    //         CRandom random(base_length);
    // ```
    // NCBI packs the subject it searches, the interval of `-subject_loc`, and seeds the resolution
    // of its ambiguous bases with the interval's length. With a range the reader packs the records
    // `cut_subjects` makes (at one thread and on the packing thread): their letters, deflines,
    // placements and resolved bases; the whole record packed and read in the interval resolves
    // the ambiguous bases otherwise.
    #[test]
    fn a_range_packs_the_records_cut_to_it() {
        let mut state = 5u64;
        let mut fasta = Vec::new();
        for (index, length) in [3_000usize, 2_500, 4_000].into_iter().enumerate() {
            fasta.extend_from_slice(format!(">s{index} title {index}\n").as_bytes());
            for _ in 0..length {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let draw = (state >> 33) as usize;
                fasta.push(if draw % 8 == 0 {
                    b"RYKMSWBDHVNrykmswbdhvn"[draw / 8 % 22]
                } else {
                    b"ACGTacgt"[draw / 8 % 8]
                });
            }
            fasta.push(b'\n');
        }
        let range = SequenceRange {
            from: 100,
            to: 2_800,
        };
        let config = || ReaderConfig::subject("TBLASTN", false, false);
        let whole = read_subjects(
            &mut FastaInputSource::from_bytes(&fasta, config()),
            Some(&range),
            &mut |_| Ok(()),
        )
        .unwrap();
        let (cut, expected_placements) = cut_subjects(&whole, Some(&range)).unwrap().unwrap();
        let mut resolved_otherwise = 0;
        for threads in [1, 2] {
            let (records, stores, placements) = x_read_subject_stores(
                &mut FastaInputSource::from_bytes(&fasta, config()),
                Some(&range),
                &mut |_| Ok(()),
                threads,
            )
            .unwrap();
            let stores = stores.expect("every record packed");
            assert_eq!(placements.0, expected_placements.0, "{threads} threads");
            assert_eq!(records.len(), cut.len());
            for (index, ((record, store), cut)) in records.iter().zip(&stores).zip(&cut).enumerate()
            {
                assert_eq!(
                    (&record.local_id, &record.title, record.sequence.is_empty()),
                    (&cut.local_id, &cut.title, true)
                );
                let (mut letters, mut resolved, mut expected) =
                    (Vec::new(), Vec::new(), Vec::new());
                store.ascii_into(&mut letters);
                assert_eq!(letters, cut.sequence, "{threads} threads, record {index}");
                store.resolved_letters_into(&mut resolved);
                PackedNucleotide::from_subject_letters(&cut.sequence)
                    .unwrap()
                    .resolved_letters_into(&mut expected);
                assert_eq!(resolved, expected, "{threads} threads, record {index}");
                let mut whole_resolved = Vec::new();
                PackedNucleotide::from_subject_letters(&whole[index].sequence)
                    .unwrap()
                    .resolved_letters_into(&mut whole_resolved);
                let offset = placements.offset(index);
                if whole_resolved[offset..offset + resolved.len()] != resolved[..] {
                    resolved_otherwise += 1;
                }
            }
        }
        assert!(resolved_otherwise > 0, "the seed is the interval's length");
    }
}
