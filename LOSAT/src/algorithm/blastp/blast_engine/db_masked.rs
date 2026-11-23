//! The preliminary search of a database sequence in a BLASTP search with `-db_soft_mask` or
//! `-db_hard_mask` (series NP, M4b): `run_call` chooses these kernels once per call, in place of
//! the kernel that `prelim_run::run_preliminary_with` calls, so that a search without the options
//! keeps its kernel (`prelim_kernel.rs`).
//!
//! - A sequence without masks of the algorithm is searched as in a search without the options.
//! - `-db_soft_mask`: the complement of the sequence's masks is its scan ranges; the extensions
//!   read every letter (`BlastpPrelimKernel::search_subject_in` with the ranges, the scan of
//!   `-lcase_masking`'s subjects).
//! - `-db_hard_mask`: each range between the masks is searched as a sequence of its own (a hard
//!   chunk: its scan, extensions, gapped stage and purge, or with `-ungapped` its ungapped hits),
//!   its HSPs moved by the chunk's offset and merged; then the E-values of the whole sequence and
//!   the reap by the preliminary E-value (`BlastpPrelimKernel::subject_lists`, or
//!   `BlastpPrelimKernel::ungapped_lists`).
//!
//! The masks come from `SeqDb::read_mask_block`, a window of sequences at a time
//! (`algorithm::common::db_mask_search::SubjectMasks`); the first read error by ordinal id ends
//! the preliminary stage as a read error of the sequences does (`prelim_run.rs`).

use std::sync::Mutex;

use super::prelim_kernel::BlastpPrelimKernel;
use super::{
    preliminary_score_compare, BlastpPrelimValues, BlastpSubjectRef, BlastpSubjectResult,
    BlastpSubjectScratch,
};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::common::database::search_read_error_of;
use crate::algorithm::common::db_mask_search::{SearchMask, SubjectMasks, SubjectPlan};
use crate::seqdb::SeqDb;

/// What a worker keeps for the masked kernels: the reader of the stored masks, the stored masks
/// of the sequence, its plan, the letters of a hard chunk between two sentinels, and the HSPs of
/// its chunks.
#[derive(Debug, Default)]
pub(super) struct BlastpMaskScratch {
    reader: SubjectMasks,
    stored: Vec<(u32, u32)>,
    plan: SubjectPlan,
    letters: Vec<u8>,
    hits: Vec<BlastpPreliminaryHsp>,
}

/// The database masking of a call's preliminary stage: the database, the resolved algorithm and
/// kind, the program that the errors name, and the first read error of the masks by ordinal id.
pub(super) struct BlastpDbMasking<'a> {
    pub(super) db: &'a SeqDb,
    pub(super) mask: SearchMask,
    pub(super) program: &'static str,
    pub(super) error: Mutex<Option<(u32, anyhow::Error)>>,
}

impl BlastpDbMasking<'_> {
    /// The first read error of the masks by ordinal id, taken (`None`: every read succeeded).
    pub(super) fn take_error(&self) -> Option<anyhow::Error> {
        self.error
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .take()
            .map(|(_, error)| error)
    }

    fn note_error(&self, ordinal: u32, error: anyhow::Error) {
        let mut first = self
            .error
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if first.as_ref().is_none_or(|(at, _)| ordinal < *at) {
            *first = Some((ordinal, error));
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:162-176
    // The sequence's masks are its soft ranges or its hard ranges (`s_BackupSubject`).
    /// The kernel of database sequence `s_idx` (its ordinal id) with the statistics `values` (the
    /// setup's, or the sequence's own: `prelim_run::run_preliminary_with`): `kernel`'s search
    /// without masks, over the soft ranges, or over the hard chunks.
    pub(super) fn search_subject(
        &self,
        kernel: &BlastpPrelimKernel<'_>,
        scratch: &mut BlastpSubjectScratch,
        s_idx: usize,
        subject: BlastpSubjectRef<'_>,
        subject_result: &mut BlastpSubjectResult,
        values: &BlastpPrelimValues<'_>,
    ) {
        if subject.aa_len == 0 {
            kernel.search_subject(scratch, s_idx, subject, subject_result, values);
            return;
        }
        let mut masked = scratch.db_mask.take().unwrap_or_default();
        let ordinal = s_idx as u32;
        match masked.reader.of(
            self.db,
            ordinal,
            self.mask.algo_id,
            search_read_error_of,
            self.program,
        ) {
            Ok(stored) => {
                masked.stored.clear();
                masked.stored.extend_from_slice(stored);
            }
            Err(error) => {
                subject_result.reset(s_idx);
                self.note_error(ordinal, error);
                scratch.db_mask = Some(masked);
                return;
            }
        }
        if masked.stored.is_empty() {
            kernel.search_subject(scratch, s_idx, subject, subject_result, values);
        } else if !self.mask.hard() {
            let BlastpMaskScratch { stored, plan, .. } = &mut *masked;
            let ranges = plan.soft(stored, subject.aa_len);
            kernel.search_subject_in_any(
                scratch,
                s_idx,
                subject,
                subject_result,
                values,
                Some(ranges),
                false,
            );
        } else {
            hard_chunks(
                kernel,
                scratch,
                &mut masked,
                s_idx,
                subject,
                subject_result,
                values,
            );
        }
        scratch.db_mask = Some(masked);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:476-584
// Each hard chunk is searched alone: its scan, its gapped stage, the purge of the HSPs with
// common endpoints and the sort by score; its HSPs move by the chunk's offset
// (`Blast_HSPListAdjustOffsets`) and join the sequence's list (`Blast_HSPListsMerge`). NCBI's
// overlap there is not 0: it compares the chunk's offset with the left of the range after it
// (the chunker has already moved to it, blast_engine.c:254-258, 579-580) and takes the chunk
// overlap of 100. The merge only looks at the HSPs of the joined list that end past the new
// chunk's offset (blast_hits.c:2888-2899); a chunk ends where its range does and no HSP reaches
// into the mask after it, so there are none, and the merge appends the chunk's HSPs and sorts
// the list by score (`s_BlastHSPListsCombineByScore`, blast_hits.c:2749-2766; no limit on the
// HSPs of a list in the preliminary stage). TBLASTN's and TBLASTX's chunks of a long range
// overlap; this reason does not carry over. With `-subject_besthit` the filter then runs on the
// joined list after every chunk with HSPs (blast_engine.c:586-591, blast_hits.c:2537-2606).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:869-899
// After the chunks, the E-values of the sequence's HSPs with its whole length and the reap by
// the preliminary E-value.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:560-562,1482-1541
// An ungapped search (`-ungapped`) makes each chunk's list of its ungapped hits sorted by score
// (`BLAST_GetUngappedHSPList`) and merges it as above (the merge's gap flag only matters for HSPs
// in the overlap, of which there are none); after the chunks, the HSPs of a negative length go,
// then the E-values with the whole length and the reaps (`BlastpPrelimKernel::ungapped_lists`).
/// The hard chunks of a sequence with masks (`masked.stored`), searched one by one into the HSPs
/// of the sequence, and its HSP lists per query in `subject_result`. A chunk shorter than the
/// word finds nothing in NCBI's one-hit and two-hit scans alike and is passed over
/// (`SubjectPlan::next_searched_chunk`). Every chunk is searched over its own letters, the chunk
/// from the start before a mask that runs to the end included, which NCBI scans past its end (its
/// "no chunking" branch; an approved NCBI defect, PD-LOSAT-NCBI-DEFECTS exception 21, the Owner's
/// answer Q1 = B of M4b). The sequences that LOSAT refuses never come here: a mask past the
/// sequence's end (`HardMaskShape::PastEnd`) and a range longer than `MAX_DBSEQ_LEN`
/// (`HardMaskShape::LongRange`).
fn hard_chunks(
    kernel: &BlastpPrelimKernel<'_>,
    scratch: &mut BlastpSubjectScratch,
    masked: &mut BlastpMaskScratch,
    s_idx: usize,
    subject: BlastpSubjectRef<'_>,
    subject_result: &mut BlastpSubjectResult,
    values: &BlastpPrelimValues<'_>,
) {
    subject_result.reset(s_idx);
    masked.plan.start_hard(&masked.stored, subject.aa_len);
    masked.hits.clear();
    let residues = &subject.aa_seq[1..=subject.aa_len];
    let word = usize::try_from(kernel.word_size).unwrap_or(usize::MAX);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:475-486
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:807
    // With a negative window the diagonal table keeps its state from chunk to chunk, and a chunk
    // shorter than the word moves the offset (`Blast_ExtendWordExit`) though it holds no word:
    // every chunk is searched (`search_subject_in_neg`; its scan of such a chunk reads nothing).
    let every_chunk = kernel.window < 0;
    loop {
        let next = if every_chunk {
            masked.plan.next_chunk()
        } else {
            masked.plan.next_searched_chunk(word)
        };
        let Some(chunk) = next else {
            break;
        };
        // The chunk's letters between two sentinels, as a sequence of the database block.
        masked.letters.clear();
        masked.letters.push(0);
        masked
            .letters
            .extend_from_slice(&residues[chunk.offset..chunk.offset + chunk.length]);
        masked.letters.push(0);
        let chunk_subject = BlastpSubjectRef {
            aa_seq: &masked.letters,
            aa_len: chunk.length,
        };
        kernel.search_subject_in_any(
            scratch,
            s_idx,
            chunk_subject,
            subject_result,
            values,
            None,
            true,
        );
        let offset = chunk.offset as i32;
        let chunk_hits = !scratch.preliminary_hsp_list.is_empty();
        masked
            .hits
            .extend(scratch.preliminary_hsp_list.drain(..).map(|mut hsp| {
                hsp.subject_start += offset;
                hsp.subject_end += offset;
                hsp.gapped_subject_start += offset;
                hsp
            }));
        if kernel.args.subject_besthit && chunk_hits {
            // The filter of the joined list after this chunk, in the merge's order.
            masked.hits.sort_by(preliminary_score_compare);
            kernel.subject_besthit(&mut masked.hits);
        }
    }
    // The merge's order: by score, then the subject's offset and end, then the query's (a total
    // order among one query's HSPs, which are split by query below in this order).
    masked.hits.sort_by(preliminary_score_compare);
    if kernel.args.ungapped {
        kernel.ungapped_lists(
            scratch,
            s_idx,
            subject.aa_len,
            values,
            &mut masked.hits,
            subject_result,
        );
    } else {
        kernel.subject_lists(
            scratch,
            s_idx,
            subject.aa_len,
            values,
            &mut masked.hits,
            subject_result,
        );
    }
}
