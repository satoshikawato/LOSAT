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
    /// lower-case letters, `-lcase_masking`) it is the record's letters (`with_letters`) where it
    /// has lower-case letters. Otherwise a packed subject is given as its store: the stage
    /// translates its packed bases (the bases the letters resolve to), so a worker holds the
    /// frames alone and neither the letters nor the resolved bases. A database sequence (`-db`)
    /// has no lower-case letters, so its stage always reads its stored bases, as NCBI's (no
    /// subject mask without `-db_soft_mask`).
    pub(super) fn with_search_subject<R>(
        &self,
        oid: usize,
        case: bool,
        f: impl FnOnce(SubjectNt<'_>) -> R,
    ) -> R {
        match self {
            Self::Packed(subjects) if !case || subjects[oid].lowercase().is_empty() => {
                f(SubjectNt::Packed(&subjects[oid]))
            }
            Self::Packed(_) => {
                // A packed subject whose lower-case letters the stage reads: the letters.
                crate::utils::xreach::x_reach!("LOSAT_X_TBNSTORE", false);
                self.with_letters(oid, |letters| f(SubjectNt::Letters(letters)))
            }
            Self::Records(_) => self.with_letters(oid, |letters| f(SubjectNt::Letters(letters))),
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

/// EXPERIMENT (LOSAT_X_TBNSTORE): the subject records (`read_subjects`), each one's letters packed
/// as soon as it has been read, so the letters of the whole set are never held at once. With a
/// subject range (`-subject_loc`) each record is first cut to its interval, in place, into the
/// record `seq_range::cut_subjects` makes, and the placements are returned (empty without a
/// range). The records keep their deflines and no letters. If a record holds a byte without
/// an ncbi4na code (only an API caller can make one), every record keeps its letters and no store
/// is returned: the search then reads the letters as without the switch. With `threads` above one
/// a second thread packs each record while the reader reads the next one (`x_packed_while_read`),
/// except in WebAssembly, where the reader packs each record itself as at one thread: there a thread
/// started outside the search pool is one more host worker than `-num_threads` asks for.
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
    use crate::blastinput::seq_range::{record_cut, Placements, RecordPlacement};
    let mut placements = Vec::new();
    // The record `cut_subjects` makes, cut in place (`FastaRecord::cut_in_place`). An interval that
    // starts past the end has stopped the reading before (`read_subjects_map`).
    let mut cut = |mut record: crate::blastinput::fasta_reader::FastaRecord| {
        if let Some(range) = range {
            let length = record.sequence.len();
            let (from, to_exclusive) = record_cut(range, length)?;
            record.cut_in_place(from, to_exclusive);
            placements.push(RecordPlacement {
                offset: from,
                length,
            });
        }
        anyhow::Ok(record)
    };
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let read = if threads > 1 {
        x_packed_while_read(source, range, warn, &mut cut)?
    } else {
        crate::blastinput::fasta_reader::read_subjects_map(source, range, warn, |record| {
            Ok(x_packed(cut(record)?))
        })?
    };
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let read = {
        let _ = threads;
        crate::blastinput::fasta_reader::read_subjects_map(source, range, warn, |record| {
            Ok(x_packed(cut(record)?))
        })?
    };
    let placements = Placements(placements);
    let all_packed = read.iter().all(|(_, store)| store.is_some());
    // A byte without an ncbi4na code keeps the letters (the reference path).
    crate::utils::xreach::x_reach!("LOSAT_X_TBNSTORE", all_packed);
    if all_packed {
        let (records, stores): (_, Vec<PackedNucleotide>) = read
            .into_iter()
            .map(|(record, store)| (record, store.expect("every record packed")))
            .unzip();
        STORED_SUBJECTS.store(stores.len() as u64, Ordering::Relaxed);
        STORED_BYTES.store(
            stores.iter().map(|store| store.heap_bytes() as u64).sum(),
            Ordering::Relaxed,
        );
        return Ok((records, Some(stores), placements));
    }
    // The records are already cut: their letters are those of their intervals.
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
        placements,
    ))
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
// reader reads the next ones (natively; a WebAssembly build packs on the reader, see
// `x_read_subject_stores`). The reader hands them over in batches (up to `X_PACK_BATCH`
// records or `X_PACK_BATCH_LETTERS` letters, so that a set of many short records costs few
// hand-overs) and only when that thread takes a batch (a channel without a buffer), so at most
// two batches' letters are held; the records come back in reading order. Packing is a function
// of a record alone, so the records and stores are those of the reader packing them itself.
/// EXPERIMENT (LOSAT_X_TBNSTORE): the records read (each cut by `cut`) and packed on two threads.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn x_packed_while_read<R: std::io::Read>(
    source: &mut crate::blastinput::fasta_reader::FastaInputSource<R>,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    cut: &mut dyn FnMut(
        crate::blastinput::fasta_reader::FastaRecord,
    ) -> anyhow::Result<crate::blastinput::fasta_reader::FastaRecord>,
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
            crate::blastinput::fasta_reader::read_subjects_map(source, range, warn, |record| {
                // the reader cuts a record to its interval before it hands it over; a cut record
                // keeps the room of its whole letters until it is packed, so that room is counted
                let record = cut(record)?;
                letters += record.sequence.capacity();
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
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
const X_PACK_BATCH: usize = 1024;
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
const X_PACK_BATCH_LETTERS: usize = 4 << 20;

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
