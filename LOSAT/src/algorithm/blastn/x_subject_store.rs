//! The subjects' letters as a BLASTN search reads them: from the records, or (EXPERIMENT,
//! `LOSAT_X_BNSTORE`) from packed subjects made while the records were read, whose letters are
//! made again where they are read (`sequence/packed_subject.rs`).

use std::sync::atomic::{AtomicU64, Ordering};

use crate::blastinput::fasta_reader::{FastaInputSource, FastaRecord};
use crate::core::blast_encoding::encode_iupac_to_blastna;
use crate::sequence::packed_subject::PackedNucleotide;
use crate::utils::dust::MaskedInterval;

use super::coordination::collect_lowercase_masks;

// No NCBI counterpart: reads the LOSAT_X_BNSTORE variable once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BNSTORE): off unless the variable is set (it is not in the strict set).
pub(crate) fn x_bnstore() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BNSTORE").is_some())
}

/// The letters of the subjects of a search, by subject index: the records' own, or (EXPERIMENT,
/// `LOSAT_X_BNSTORE`) the packed subjects of records that keep none. Every read gives the bytes
/// the records' letters give (the comparison tests of `packed_subject.rs`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
/// ```c
///     virtual SBlastSequence GetCompressedPlusStrand() {
///         SBlastSequence retval(size());
///         string ncbi4na = kEmptyStr;
///         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
///         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
///         return retval;
/// ```
/// NCBI keeps the subjects in its scope and makes each form a stage reads from it: the packed
/// bases of the preliminary stage, the blastna letters of the traceback, the IUPAC letters of the
/// report.
#[derive(Clone, Copy)]
pub(crate) struct SubjectLetters<'a> {
    records: &'a [FastaRecord],
    /// EXPERIMENT (LOSAT_X_BNSTORE): the letters of every record, packed (the records keep none).
    stores: Option<&'a [PackedNucleotide]>,
}

impl<'a> SubjectLetters<'a> {
    /// The letters of the records.
    pub(crate) fn records(records: &'a [FastaRecord]) -> Self {
        Self {
            records,
            stores: None,
        }
    }

    /// EXPERIMENT (LOSAT_X_BNSTORE): the letters of `records` held in `stores`, one per record.
    pub(crate) fn packed(records: &'a [FastaRecord], stores: &'a [PackedNucleotide]) -> Self {
        assert_eq!(records.len(), stores.len(), "one packed subject per record");
        Self {
            records,
            stores: Some(stores),
        }
    }

    /// EXPERIMENT (LOSAT_X_BNSTORE): the packed subjects, when the records keep no letters.
    pub(crate) fn stores(&self) -> Option<&'a [PackedNucleotide]> {
        self.stores
    }

    /// The number of subjects.
    pub(crate) fn count(&self) -> usize {
        self.records.len()
    }

    /// The number of letters of subject `index`.
    pub(crate) fn len(&self, index: usize) -> usize {
        match self.stores {
            Some(stores) => stores[index].len(),
            None => self.records[index].seq().len(),
        }
    }

    /// The letters of every subject (`len` summed).
    pub(crate) fn total(&self) -> usize {
        (0..self.count()).map(|index| self.len(index)).sum()
    }

    /// The lower-case intervals of subject `index` (`collect_lowercase_masks` of its letters).
    pub(crate) fn lowercase_masks(&self, index: usize) -> Vec<MaskedInterval> {
        match self.stores {
            Some(stores) => stores[index].lowercase().to_vec(),
            None => collect_lowercase_masks(self.records[index].seq()),
        }
    }

    /// The blastna letters of subject `index` (`encode_iupac_to_blastna` of its letters).
    pub(crate) fn blastna(&self, index: usize) -> Vec<u8> {
        match self.stores {
            Some(stores) => {
                BLASTNA_MADE.fetch_add(1, Ordering::Relaxed);
                let mut letters = Vec::new();
                stores[index].blastna_into(&mut letters);
                letters
            }
            None => encode_iupac_to_blastna(self.records[index].seq()),
        }
    }

    /// The letters `[from, to)` of subject `index` in upper case (the letters of a report row).
    pub(crate) fn upper_letters(&self, index: usize, from: usize, to: usize) -> Vec<u8> {
        match self.stores {
            Some(stores) => {
                ROWS_MADE.fetch_add(1, Ordering::Relaxed);
                let mut letters = Vec::new();
                stores[index].letters_into(from, to, &mut letters);
                letters
            }
            None => self.records[index].seq()[from..to].to_ascii_uppercase(),
        }
    }
}

/// EXPERIMENT (LOSAT_X_BNSTORE): the subject records (`read_subjects` without a range), each one's
/// letters packed as soon as it has been read, so the letters of the whole set are never held at
/// once. The records keep their deflines and messages and no letters. If a record holds a byte
/// without an ncbi4na code (the CLI reader admits none), every record keeps its letters and no
/// packed subject is returned: the search then reads the letters as without the switch.
pub(crate) fn x_read_subject_stores<R: std::io::Read>(
    source: &mut FastaInputSource<R>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> anyhow::Result<(Vec<FastaRecord>, Option<Vec<PackedNucleotide>>)> {
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
/// Subjects read packed, and the bytes their packed forms hold.
static STORED_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static STORED_BYTES: AtomicU64 = AtomicU64::new(0);
/// The times the blastna letters of a packed subject were made (its traceback).
static BLASTNA_MADE: AtomicU64 = AtomicU64::new(0);
/// The report rows whose subject letters were made from a packed subject.
static ROWS_MADE: AtomicU64 = AtomicU64::new(0);

/// The summary line (printed by `main` at exit) with LOSAT_X_STATS: the subjects held packed, the
/// bytes of their packed forms, the times their blastna letters were made and the report rows
/// made from them.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] BNSTORE subjects={} stored_bytes={} blastna_made={} rows_made={}",
            STORED_SUBJECTS.load(Ordering::Relaxed),
            STORED_BYTES.load(Ordering::Relaxed),
            BLASTNA_MADE.load(Ordering::Relaxed),
            ROWS_MADE.load(Ordering::Relaxed)
        );
    }
}
