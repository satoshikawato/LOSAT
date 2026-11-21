//! The preliminary search of a database sequence in a BLASTX search with `-db_soft_mask` or
//! `-db_hard_mask` (series NP, M4b), which `compute_subject_scored` chooses with one test of the
//! search's options per subject; a search without the options keeps its path.
//!
//! - A sequence without masks of the algorithm is searched as in a search without the options.
//! - `-db_soft_mask`: the complement of the sequence's masks is its scan ranges (the reference
//!   word finder over ranges, `x_soft_subject::x_word_finder`, as for `-lcase_masking`'s
//!   subjects); the extensions read every letter.
//! - `-db_hard_mask`: each range between the masks is searched as a sequence of its own (a hard
//!   chunk: its scan, extensions and purge), its HSPs moved by the chunk's offset and merged by
//!   score. The sink then computes the E-values, the links and the reap of the merged list with
//!   the whole sequence's length (`Runtime::subject`), as NCBI does after the chunks.
//!
//! The masks come from `SeqDb::read_mask_block` a window of sequences at a time
//! (`algorithm::common::db_mask_search::SubjectMasks`, one reader per worker's diagonal table); a
//! read error is the reader's error of the search (`x_db_stream::XReadError`).

use anyhow::Result;

use super::super::args::ResolvedOptions;
use super::super::parameters::ContextParameters;
use super::super::preliminary::{sort_by_score, PreliminaryHsp};
use super::super::query_setup::PreparedQueryBatch;
use super::super::results::subject_besthit_preliminary;
use super::super::seed::{BxLookup, BxScoring, Diagonals};
use super::super::subjects::Subjects;
use super::super::x_db_stream::XReadError;
use super::{compute_subject, compute_subject_of, subject_hits_with, SubjectComputation};
use crate::algorithm::common::database::search_read_error_of;
use crate::algorithm::common::db_mask_search::{SearchMask, SubjectPlan};
use crate::seqdb::SeqDbError;

/// The program that the reader's errors name.
const PROGRAM: &str = "BLASTX";

/// A read error of the masks during the search, as a read error of the sequences
/// (`x_db_stream::read_error`).
fn mask_read_error(error: &SeqDbError, program: &str) -> anyhow::Error {
    XReadError(search_read_error_of(error, program)).into()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:162-176
// The sequence's masks are its soft ranges or its hard ranges (`s_BackupSubject`).
/// The preliminary search of database sequence `oid` (its ordinal id) in a search masked with
/// `mask`: without masks, over the soft ranges, or over the hard chunks.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(super) fn compute_masked_subject(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
    mask: SearchMask,
) -> Result<SubjectComputation> {
    let blosum62 = options.scoring == crate::config::ScoringMatrix::Blosum62;
    let plain = |diagonals: &mut Diagonals| {
        if blosum62 {
            compute_subject(
                subjects,
                batch,
                params,
                options,
                lookup,
                diagonals,
                xdrop,
                oid,
                x_containment,
            )
        } else {
            compute_subject_of(
                subjects,
                batch,
                params,
                options,
                lookup,
                diagonals,
                xdrop,
                oid,
                x_containment,
            )
        }
    };
    let Some(db) = subjects.seqdb() else {
        return plain(diagonals);
    };
    let length = subjects.length(oid);
    if length == 0 {
        return plain(diagonals);
    }
    let stored = diagonals
        .x_db_masks
        .of(db, oid as u32, mask.algo_id, mask_read_error, PROGRAM)?
        .to_vec();
    if stored.is_empty() {
        return plain(diagonals);
    }
    let masked = Masked {
        subjects,
        batch,
        params,
        options,
        lookup,
        xdrop,
        oid,
        x_containment,
        stored: &stored,
        hard: mask.hard(),
    };
    if blosum62 {
        masked.compute(super::super::seed::Blosum62Bx, diagonals)
    } else {
        masked.compute(super::super::seed::MatrixBx::of(options.scoring), diagonals)
    }
}

/// A masked sequence's search: the arguments of `compute_subject_with` and its stored masks.
struct Masked<'a> {
    subjects: Subjects<'a>,
    batch: &'a PreparedQueryBatch,
    params: &'a [ContextParameters],
    options: &'a ResolvedOptions,
    lookup: &'a BxLookup,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
    stored: &'a [(u32, u32)],
    hard: bool,
}

impl Masked<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:476-591
    // Each hard chunk is searched alone: its scan, its extensions, the purge of the HSPs with
    // common endpoints (gapped) and the sort by score; its HSPs move by the chunk's offset
    // (`Blast_HSPListAdjustOffsets`) and join the sequence's list (`Blast_HSPListsMerge`).
    // NCBI's overlap there is not 0: it compares the chunk's offset with the left of the range
    // after it (the chunker has already moved to it, blast_engine.c:254-258, 579-580) and takes
    // the chunk overlap of 100. The merge only looks at the HSPs of the joined list that end past
    // the new chunk's offset (blast_hits.c:2888-2899); a chunk ends where its range does and no
    // HSP reaches into the mask after it, so there are none and the merge appends the chunk's
    // HSPs and sorts the list by score (`s_BlastHSPListsCombineByScore`, blast_hits.c:2749-2766).
    // TBLASTN's and TBLASTX's chunks of a long range overlap; this reason does not carry over.
    // With `-subject_besthit` the filter then runs on the joined list after every chunk with HSPs
    // (blast_engine.c:586-591, blast_hits.c:2538-2580), not once on the whole sequence's list.
    /// The sequence's computation with the matrix of `scoring`: the link cutoffs of its whole
    /// length (`compute_subject_with`), and its HSPs over the soft ranges or the hard chunks. A
    /// chunk shorter than the word finds nothing in NCBI's one-hit and two-hit scans alike and is
    /// passed over (`SubjectPlan::next_searched_chunk`).
    fn compute<S: BxScoring>(
        &self,
        scoring: S,
        diagonals: &mut Diagonals,
    ) -> Result<SubjectComputation> {
        let Masked {
            subjects,
            batch,
            params,
            options,
            lookup,
            xdrop,
            oid,
            x_containment,
            stored,
            hard,
        } = *self;
        let x_lead = batch.x_ungap_sent;
        let x_buffer = if x_lead {
            subjects.x_encoded_with_sentinels(oid)
        } else {
            subjects.encoded_with_sentinel(oid)
        };
        let residues = &x_buffer[usize::from(x_lead)..x_buffer.len() - 1];
        let length = residues.len();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455
        // The link cutoffs of an ungapped search with the whole sequence's length.
        let prepared_link = if super::super::linking::computes_link_cutoffs(options) {
            Some(super::super::linking::link_cutoffs(
                batch,
                params,
                length as i32,
                subjects.total_length()?,
                None,
            ))
        } else {
            None
        };
        let mut containment = Vec::new();
        let mut plan = SubjectPlan::default();
        if !hard {
            let ranges = plan.soft(stored, length);
            let (raw, purged) = subject_hits_with(
                scoring,
                batch,
                params,
                options,
                lookup,
                diagonals,
                xdrop,
                &x_buffer,
                x_lead,
                Some(ranges),
                x_containment,
                &mut containment,
            )?;
            return Ok(SubjectComputation {
                raw,
                purged,
                containment,
                prepared_link,
            });
        }
        let word = usize::try_from(options.word_size).unwrap_or(usize::MAX);
        plan.start_hard(stored, length);
        let mut raw = Vec::new();
        let mut purged = Vec::new();
        let mut letters = Vec::new();
        while let Some(chunk) = plan.next_searched_chunk(word) {
            // The chunk's letters with the sentinel after them, and before them with `x_lead`.
            letters.clear();
            if x_lead {
                letters.push(0);
            }
            letters.extend_from_slice(&residues[chunk.offset..chunk.offset + chunk.length]);
            letters.push(0);
            let mut chunk_containment = Vec::new();
            let (chunk_raw, chunk_purged) = subject_hits_with(
                scoring,
                batch,
                params,
                options,
                lookup,
                diagonals,
                xdrop,
                &letters,
                x_lead,
                None,
                x_containment,
                &mut chunk_containment,
            )?;
            let offset = chunk.offset as i32;
            let chunk_hits = !chunk_purged.is_empty();
            raw.extend(chunk_raw.into_iter().map(|hsp| shifted(hsp, offset)));
            purged.extend(chunk_purged.into_iter().map(|hsp| shifted(hsp, offset)));
            containment.extend(
                chunk_containment
                    .into_iter()
                    .map(|(hsp, contained)| (shifted(hsp, offset), contained)),
            );
            if options.subject_besthit && chunk_hits {
                // The filter of the joined list after this chunk; the sink's own filter of the
                // whole list after the last chunk deletes nothing more (the filter keeps no pair
                // that it would delete).
                sort_by_score(&mut purged);
                subject_besthit_preliminary(&mut purged);
            }
        }
        sort_by_score(&mut purged);
        Ok(SubjectComputation {
            raw,
            purged,
            containment,
            prepared_link,
        })
    }
}

/// `hsp` moved by a chunk's offset `offset` (`Blast_HSPListAdjustOffsets`: the subject's start,
/// end and gapped start).
fn shifted(mut hsp: PreliminaryHsp, offset: i32) -> PreliminaryHsp {
    hsp.s_start += offset;
    hsp.s_end += offset;
    hsp.s_gapped_start += offset;
    hsp
}
