//! EXPERIMENT (LOSAT_X_TBNSTORE, LOSAT_X_TBXSTORE): the nucleotide records of `-subject`, each one's
//! letters packed (`sequence::packed_subject::PackedNucleotide`) as soon as it has been read, for
//! the programs that search them from their stores (TBLASTN `tblastn/x_subject_store.rs`, TBLASTX
//! `tblastx/subjects.rs`; followups-20261014 DESIGN_subject_store.md). Each caller records its own
//! switch's side.

use crate::sequence::packed_subject::PackedNucleotide;

/// EXPERIMENT (LOSAT_X_TBNSTORE, LOSAT_X_TBXSTORE): the subject records (`read_subjects`), each one's
/// letters packed as soon as it has been read, so the letters of the whole set are never held at
/// once. With a subject range (`-subject_loc`) each record is first cut to its interval, in place,
/// into the record `seq_range::cut_subjects` makes, and the placements are returned (empty without
/// a range). The records keep their deflines and no letters. If a record holds a byte without an
/// ncbi4na code (`iupacna_to_ncbi4na`: `U` among them), every record keeps its letters and no store
/// is returned: the search then reads the letters as without the switch. With `threads` above one
/// a second thread packs each record while the reader reads the next one (`x_packed_while_read`),
/// except in WebAssembly, where the reader packs each record itself as at one thread: there a thread
/// started outside the search pool is one more host worker than `-num_threads` asks for.
pub(crate) fn x_read_packed_subjects<R: std::io::Read>(
    source: &mut crate::blastinput::fasta_reader::FastaInputSource<R>,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    threads: usize,
    program: &'static str,
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
        x_packed_while_read(source, range, warn, &mut cut, program)?
    } else {
        crate::blastinput::fasta_reader::read_subjects_map(source, range, warn, |record| {
            Ok(x_packed(cut(record)?))
        })?
    };
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let read = {
        let _ = (threads, program);
        crate::blastinput::fasta_reader::read_subjects_map(source, range, warn, |record| {
            Ok(x_packed(cut(record)?))
        })?
    };
    let placements = Placements(placements);
    if read.iter().all(|(_, store)| store.is_some()) {
        let (records, stores): (_, Vec<PackedNucleotide>) = read
            .into_iter()
            .map(|(record, store)| (record, store.expect("every record packed")))
            .unzip();
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
/// EXPERIMENT (LOSAT_X_TBNSTORE, LOSAT_X_TBXSTORE): the records read (each cut by `cut`) and packed
/// on two threads.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn x_packed_while_read<R: std::io::Read>(
    source: &mut crate::blastinput::fasta_reader::FastaInputSource<R>,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    cut: &mut dyn FnMut(
        crate::blastinput::fasta_reader::FastaRecord,
    ) -> anyhow::Result<crate::blastinput::fasta_reader::FastaRecord>,
    program: &'static str,
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
        let ended = || anyhow::anyhow!("{program}: the subject packing thread ended");
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
            .map_err(|_| anyhow::anyhow!("{program}: the subject packing thread panicked"))?;
        read?;
        last?;
        Ok(packed)
    })
}

/// EXPERIMENT (LOSAT_X_TBNSTORE, LOSAT_X_TBXSTORE): the records and the letters of one hand-over of
/// `x_packed_while_read`.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
const X_PACK_BATCH: usize = 1024;
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
const X_PACK_BATCH_LETTERS: usize = 4 << 20;
