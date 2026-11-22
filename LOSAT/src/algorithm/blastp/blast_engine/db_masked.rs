//! The preliminary search of a database sequence in a BLASTP search with `-db_soft_mask` or
//! `-db_hard_mask` (series NP, M4b): `run_call` chooses these kernels once per call, in place of
//! the kernel that `prelim_run::run_preliminary_with` calls, so that a search without the options
//! keeps its kernel (`prelim_kernel.rs`).
//!
//! - A sequence without masks of the algorithm is searched as in a search without the options
//!   (with `-db_hard_mask` on a new table first where its hits meet entries that no sequence wrote
//!   in NCBI's one-thread search: `algorithm::common::hard_chunk_walk::NcbiWalk`).
//! - `-db_soft_mask`: the complement of the sequence's masks is its scan ranges; the extensions
//!   read every letter (`BlastpPrelimKernel::search_subject_in` with the ranges, the scan of
//!   `-lcase_masking`'s subjects).
//! - `-db_hard_mask`: each range between the masks is searched as a sequence of its own (a hard
//!   chunk: its scan, extensions, gapped stage and purge, or with `-ungapped` its ungapped hits),
//!   its HSPs moved by the chunk's offset and merged; then the E-values of the whole sequence and
//!   the reap by the preliminary E-value (`BlastpPrelimKernel::subject_lists`, or
//!   `BlastpPrelimKernel::ungapped_lists`). Every chunk is walked, whatever its length (masks may
//!   overlap or nest, or be out of the order of their start; series NP, M4b-3b): a chunk shorter
//!   than the word gives the two-hit search its one word (`common::hard_chunk_scan::OneWord`),
//!   whose seeds go through the chunk's gapped stage or are its ungapped HSPs, and the one-hit
//!   search nothing, and the diagonal offset moves by each chunk's length plus the window.
//!
//! The masks come from `SeqDb::read_mask_block`, a window of sequences at a time
//! (`algorithm::common::db_mask_search::SubjectMasks`); the first read error by ordinal id ends
//! the preliminary stage as a read error of the sequences does (`prelim_run.rs`), and so do
//! LOSAT's refusals of a sequence whose short chunk's HSPs NCBI's interval tree splits without
//! end or its merge joins with another chunk's (`common::hard_chunk_scan::
//! short_chunk_tree_refusal`, `short_chunk_merge_refusal`).

use std::sync::Mutex;

use super::prelim_kernel::BlastpPrelimKernel;
use super::{
    blastp_context_idx_for_query_offset_with_bounds, blastp_save_init_hsp,
    preliminary_score_compare, x_soft_subject, BlastpPrelimValues, BlastpSubjectRef,
    BlastpSubjectResult, BlastpSubjectScratch, InitHSP,
};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::common::database::search_read_error_of;
use crate::algorithm::common::db_mask_search::{SearchMask, SubjectMasks, SubjectPlan};
use crate::algorithm::common::hard_chunk_scan::{
    merge_joins, short_chunk_merge_refusal, short_chunk_tree_refusal, OneWord, WordContext,
};
use crate::algorithm::common::hard_chunk_walk::{
    ChunkWork, DiagonalWalk, HardStarts, HardTable, NcbiStart,
};
use crate::core::subject_split::HardChunk;
use crate::seqdb::SeqDb;

/// What a worker keeps for the masked kernels: the reader of the stored masks, the stored masks
/// of the sequence, its plan and its hard chunks, the letters of a hard chunk between two
/// sentinels (or of a short chunk's one word), the query offsets of that word's hits, and the
/// HSPs of its chunks.
#[derive(Debug, Default)]
pub(super) struct BlastpMaskScratch {
    reader: SubjectMasks,
    stored: Vec<(u32, u32)>,
    plan: SubjectPlan,
    chunks: Vec<HardChunk>,
    letters: Vec<u8>,
    word_hits: Vec<i32>,
    hits: Vec<BlastpPreliminaryHsp>,
}

/// The database masking of a call's preliminary stage: the database, the resolved algorithm and
/// kind, with `-db_hard_mask` the sequences that the search starts on a new table
/// (`blastinput::db_mask::ResolvedMask::hard_starts`), the program that the errors name, and the
/// first read error of the masks by ordinal id.
pub(super) struct BlastpDbMasking<'a> {
    pub(super) db: &'a SeqDb,
    pub(super) mask: SearchMask,
    pub(super) starts: Option<&'a HardStarts>,
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
        let ncbi = self.starts.and_then(|starts| starts.of(ordinal));
        if masked.stored.is_empty() {
            match ncbi {
                Some(start) => plain_on_ncbi_table(
                    kernel,
                    scratch,
                    s_idx,
                    subject,
                    subject_result,
                    values,
                    start,
                ),
                None => kernel.search_subject(scratch, s_idx, subject, subject_result, values),
            }
        } else if !self.mask.hard() {
            let BlastpMaskScratch { stored, plan, .. } = &mut *masked;
            let ranges = plan.soft(stored, subject.aa_len);
            kernel.search_subject_in(
                scratch,
                s_idx,
                subject,
                subject_result,
                values,
                Some(ranges),
                false,
            );
        } else {
            let refused = hard_chunks(
                kernel,
                scratch,
                &mut masked,
                s_idx,
                ncbi,
                subject,
                subject_result,
                values,
            );
            if let Err(refused) = refused {
                subject_result.reset(s_idx);
                let oid = self.db.oid(ordinal);
                self.note_error(
                    ordinal,
                    match refused {
                        ShortChunkRefusal::Tree(length) => {
                            short_chunk_tree_refusal(oid, self.program, length)
                        }
                        ShortChunkRefusal::Merge => short_chunk_merge_refusal(oid, self.program),
                    },
                );
            }
        }
        scratch.db_mask = Some(masked);
    }
}

/// A sequence without masks of the algorithm in a `-db_hard_mask` search that starts on a new table
/// at NCBI's one-thread offset `start` (`common::hard_chunk_walk::TableStart::New`: an offset that
/// the walks before it took within the window of the entries that no sequence wrote), searched as
/// in a search without the options on that table; by the reference scan of its whole range where
/// the entries lie above the offset (the other scans take every entry to lie below it). Its
/// entries end at least the window below the offset after it, so the next sequence starts as
/// after any sequence without masks.
#[cold]
#[inline(never)]
fn plain_on_ncbi_table(
    kernel: &BlastpPrelimKernel<'_>,
    scratch: &mut BlastpSubjectScratch,
    s_idx: usize,
    subject: BlastpSubjectRef<'_>,
    subject_result: &mut BlastpSubjectResult,
    values: &BlastpPrelimValues<'_>,
    start: NcbiStart,
) {
    scratch.ncbi_table(start, kernel.window);
    if start.unwritten(kernel.window) > start.offset {
        let whole = [(0, subject.aa_len as i32)];
        kernel.search_subject_in(
            scratch,
            s_idx,
            subject,
            subject_result,
            values,
            Some(&whole),
            false,
        );
    } else {
        kernel.search_subject(scratch, s_idx, subject, subject_result, values);
    }
}

/// Why the search of a hard-masked sequence stops at a chunk shorter than the word (the caller's
/// refusal): NCBI's insertion of one of its HSPs into the interval tree of the chunk of that length
/// never ends, or NCBI's merge joins its HSPs with another chunk's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ShortChunkRefusal {
    Tree(i32),
    Merge,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:476-584
// Each hard chunk is searched alone: its scan, its gapped stage, the purge of the HSPs with common
// endpoints and the sort by score; its HSPs move by the chunk's offset
// (`Blast_HSPListAdjustOffsets`) and join the sequence's list (`Blast_HSPListsMerge`). NCBI's
// overlap there is not 0: it compares the chunk's offset with the left of the range after it (the
// chunker has already moved to it, blast_engine.c:254-258, 579-580) and takes the chunk overlap of
// 100. The merge only looks at the HSPs of the joined list that end past the new chunk's offset
// (blast_hits.c:2888-2899); a chunk ends where its range does and no HSP reaches into the mask
// after it, and no chunk that can hold HSPs starts before the end of an earlier one (else
// `HardMaskShape::PartsOutOfOrder`), so there are none, and the merge appends the chunk's HSPs and
// sorts the list by score (`s_BlastHSPListsCombineByScore`, blast_hits.c:2749-2766; no limit on the
// HSPs of a list in the preliminary stage). TBLASTN's and TBLASTX's chunks of a long range overlap;
// this reason does not carry over. With `-subject_besthit` the filter then runs on the joined list
// after every chunk with HSPs (blast_engine.c:586-591, blast_hits.c:2537-2606).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:869-899
// After the chunks, the E-values of the sequence's HSPs with its whole length and the reap by
// the preliminary E-value.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:560-562,1482-1541
// An ungapped search (`-ungapped`) makes each chunk's list of its ungapped hits sorted by score
// (`BLAST_GetUngappedHSPList`) and merges it as above (the merge's gap flag only matters for HSPs
// in the overlap, of which there are none); after the chunks, the HSPs of a negative length go,
// then the E-values with the whole length and the reaps (`BlastpPrelimKernel::ungapped_lists`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:476-500
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
// Every chunk goes through the word finder, whose exit moves the diagonal offset by the chunk's
// length (negative between overlapping masks) plus the window, or clears the table where the
// offset has reached INT4_MAX / 4, after any chunk (inside a sequence too). On a new table at
// NCBI's one-thread offset the clear falls after the same chunk as NCBI's; elsewhere the offsets
// are NCBI's plus the moves before (and a worker's start), so a clear can fall a sequence or more
// earlier than NCBI's, never inside a sequence (`DiagonalWalk::start` clears first), where it
// changes no row.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-562
// A chunk shorter than the word goes through the same stages: the seeds of its one word are its
// ungapped hits, which an ungapped search makes its HSPs, moved by the chunk's offset and merged
// as any chunk's.
/// The hard chunks of a sequence with masks (`masked.stored`) at the place `place` among the
/// searched sequences, searched one by one into the HSPs of the sequence, and its HSP lists per
/// query in `subject_result`; `Err` where the one word of a chunk shorter than the word gives
/// HSPs that NCBI's interval tree splits without end or that its merge joins with another chunk's
/// (LOSAT's refusals, the caller's: `common::hard_chunk_scan::short_chunk_tree_refusal`,
/// `short_chunk_merge_refusal`). Every chunk is walked in order
/// (`common::hard_chunk_walk::DiagonalWalk`): on a new table at NCBI's one-thread offset where its
/// hits meet the entries that no sequence wrote in NCBI's one-thread search or that search clears
/// the table inside it (`ncbi`, `common::hard_chunk_walk::NcbiWalk`), else with the diagonal offset
/// moved so that no entry of another sequence reaches its hits or theirs (`DiagonalWalk::start`,
/// `DiagonalWalk::after`; PD-LOSAT-NCBI-DEFECTS exception 23, the Owner's answer Q5 = B of M4b,
/// 2026-11-22). A
/// chunk that holds a word is searched over its own
/// letters (by the reference scan where its diagonal offset went back below an entry of the
/// table), the chunk from the start before a mask that runs to the end included, which NCBI scans
/// past its end (its "no chunking" branch; an approved NCBI defect, PD-LOSAT-NCBI-DEFECTS
/// exception 21, the Owner's answer Q1 = B of M4b); a shorter one gives the two-hit search its one
/// word, whose seeds go through the chunk's gapped stage, or are its HSPs in an ungapped search,
/// and the one-hit search nothing. The sequences that LOSAT refuses never come here: a mask that
/// ends before it starts (`HardMaskShape::Inverted`), a chunk that can hold alignments starting
/// before the end of an earlier one (`HardMaskShape::PartsOutOfOrder`), a mask past the sequence's
/// end (`HardMaskShape::PastEnd`), a range longer than `MAX_DBSEQ_LEN` (`HardMaskShape::
/// LongRange`), a mask that runs to the end after the chunk from the start with masks stored after
/// it (`HardMaskShape::TailMaskWithLaterMasks`) and a one word past the sentinel that can change a
/// result (`HardMaskShape::WordPastEnd`), so the sequence's result does not depend on the entries
/// of the sequences that the diagonal table searched before or after it.
#[allow(clippy::too_many_arguments)]
fn hard_chunks(
    kernel: &BlastpPrelimKernel<'_>,
    scratch: &mut BlastpSubjectScratch,
    masked: &mut BlastpMaskScratch,
    s_idx: usize,
    ncbi: Option<NcbiStart>,
    subject: BlastpSubjectRef<'_>,
    subject_result: &mut BlastpSubjectResult,
    values: &BlastpPrelimValues<'_>,
) -> Result<(), ShortChunkRefusal> {
    subject_result.reset(s_idx);
    masked.plan.start_hard(&masked.stored, subject.aa_len);
    masked.plan.hard_chunks(&mut masked.chunks);
    masked.hits.clear();
    let residues = &subject.aa_seq[1..=subject.aa_len];
    let (word, window) = (kernel.word_size, kernel.window);
    // The search's walk reads only what each chunk's scan reads and whether it went back: a one
    // word that covers the sentinel is taken to store entries (its lookup finds no hit where the
    // query holds no residue 0, `DiagonalWalk::new`), so `reference` is exact or more, and its
    // margins too (a larger margin moves the offset further past entries that no hit reaches).
    let (mut walk, start) = DiagonalWalk::start(
        word,
        window,
        subject.aa_len,
        &masked.chunks,
        ncbi,
        scratch.diag_offset,
    );
    let Ok(()) = start.apply(scratch, window);
    for index in 0..masked.chunks.len() {
        let chunk = masked.chunks[index];
        // The exit of this chunk's word finder clears the table (`finish_chunk`).
        let clears = scratch.diag_offset >= i32::MAX / 4;
        let searched = match walk.chunk(chunk.offset, chunk.length) {
            ChunkWork::Words { reference } => {
                // The chunk's letters between two sentinels, as a sequence of the database block.
                let letters = chunk.letters();
                masked.letters.clear();
                masked.letters.push(0);
                masked
                    .letters
                    .extend_from_slice(&residues[chunk.offset..chunk.offset + letters]);
                masked.letters.push(0);
                let chunk_subject = BlastpSubjectRef {
                    aa_seq: &masked.letters,
                    aa_len: letters,
                };
                let whole = [(0, chunk.length)];
                kernel.search_subject_in(
                    scratch,
                    s_idx,
                    chunk_subject,
                    subject_result,
                    values,
                    reference.then_some(&whole[..]),
                    true,
                );
                true
            }
            ChunkWork::OneWord => {
                let seeds = one_word(kernel, scratch, masked, residues, chunk, values);
                scratch.finish_chunk(chunk.length, window);
                if seeds.is_empty() {
                    scratch.init_hsps = seeds;
                    false
                } else if kernel.args.ungapped {
                    kernel.one_word_hsps(scratch, s_idx, subject_result, values, seeds);
                    true
                } else if kernel.one_word_gapped(
                    scratch,
                    s_idx,
                    subject_result,
                    values,
                    seeds,
                    &masked.letters,
                    chunk.length,
                ) {
                    true
                } else {
                    return Err(ShortChunkRefusal::Tree(chunk.length));
                }
            }
            ChunkWork::Nothing => {
                scratch.finish_chunk(chunk.length, window);
                false
            }
        };
        if clears {
            walk.clear();
        }
        if !searched {
            continue;
        }
        let offset = chunk.offset as i32;
        let chunk_hits = !scratch.preliminary_hsp_list.is_empty();
        // Only the HSPs of a chunk shorter than the word can meet another chunk's
        // (`HardMaskShape::PartsOutOfOrder` refuses chunks that hold words and meet).
        if chunk_hits
            && merge_joins(
                &chunk,
                masked.hits.iter().map(|hsp| hsp.subject_end),
                scratch
                    .preliminary_hsp_list
                    .iter()
                    .map(|hsp| hsp.subject_start + offset),
            )
        {
            scratch.preliminary_hsp_list.clear();
            return Err(ShortChunkRefusal::Merge);
        }
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
    let Ok(()) = walk.after(scratch.diag_offset).apply(scratch, window);
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
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-505
// The two-hit scan of a chunk shorter than the word reads the one word at the chunk's offset, from
// the sequence's stored letters (the preliminary stage's fetch is not masked).
/// The one word of a hard chunk `chunk` shorter than the word in a two-hit search
/// (`common::hard_chunk_scan::OneWord`): the stored letters from the chunk's offset (`residues`,
/// the
/// sequence's; the sentinel after the last letter, zeros past it, which only a word whose hits
/// change nothing reads: `HardMaskShape::WordPastEnd`), their lookup hits by the reference scan of
/// one word, and NCBI's two-hit test of each on the scratch's diagonal table with the contexts'
/// X-drops and cutoffs. Returns the seeds that it saved (the scratch's list of initial hits, taken;
/// the caller gives it back), with the chunk's subject offsets; `masked.letters` keeps the word's
/// letters and the sentinel after them, the subject of the chunk's gapped stage.
fn one_word(
    kernel: &BlastpPrelimKernel<'_>,
    scratch: &mut BlastpSubjectScratch,
    masked: &mut BlastpMaskScratch,
    residues: &[u8],
    chunk: HardChunk,
    values: &BlastpPrelimValues<'_>,
) -> Vec<InitHSP> {
    let word = usize::try_from(kernel.word_size).unwrap_or(0);
    let stored = &residues[chunk.offset..(chunk.offset + word).min(residues.len())];
    masked.letters.clear();
    masked.letters.extend_from_slice(stored);
    // The word and the sentinel after it, as the scan reads a subject.
    masked.letters.resize(word + 1, 0);
    let hits = &mut masked.word_hits;
    hits.clear();
    x_soft_subject::x_for_each_hit(
        kernel.lookup,
        &masked.letters,
        word,
        &[(0, kernel.word_size)],
        &mut scratch.offset_pairs,
        true,
        &mut |pair, _| hits.push(pair.q_off as i32),
    );
    let test = OneWord {
        query: kernel.concatenated_query.as_slice(),
        letters: &masked.letters[..word],
        length: chunk.length,
        word: kernel.word_size,
        window: kernel.window,
        scores: crate::utils::matrix::protein_score_rows(kernel.args.scoring.matrix),
    };
    let context = |q: i32| {
        let (index, start) = kernel.x_context_table.get_with_start(q).unwrap_or_else(|| {
            let index = blastp_context_idx_for_query_offset_with_bounds(
                kernel.contexts.as_slice(),
                kernel.context_lookup_bounds,
                q,
            );
            (index, kernel.contexts[index].frame_base)
        });
        WordContext {
            start,
            x_drop: kernel.x_dropoff_per_context[index],
            cutoff: values.word_cutoff_scores[index],
        }
    };
    let mut seeds = std::mem::take(&mut scratch.init_hsps);
    seeds.clear();
    test.test(
        &masked.word_hits,
        &mut scratch.diag_array,
        kernel.diag_mask,
        scratch.diag_offset,
        context,
        |q, seed| {
            blastp_save_init_hsp(
                &mut seeds,
                seed.q_start,
                seed.s_start,
                q,
                0,
                seed.length,
                seed.score,
            )
        },
    );
    seeds
}
