//! The preliminary search of a database sequence in a BLASTX search with `-db_soft_mask` or
//! `-db_hard_mask` (series NP, M4b), which `compute_subject_scored` chooses with one test of the
//! search's options per subject; a search without the options keeps its path.
//!
//! - A sequence without masks of the algorithm is searched as in a search without the options
//!   (with `-db_hard_mask` on a new table first where its hits meet entries that no sequence wrote
//!   in NCBI's one-thread search: `algorithm::common::hard_chunk_walk::NcbiWalk`).
//! - `-db_soft_mask`: the complement of the sequence's masks is its scan ranges (the reference
//!   word finder over ranges, `x_soft_subject::x_word_finder`, as for `-lcase_masking`'s
//!   subjects); the extensions read every letter.
//! - `-db_hard_mask`: each range between the masks is searched as a sequence of its own (a hard
//!   chunk: its scan, extensions and purge), its HSPs moved by the chunk's offset and merged by
//!   score. The sink then computes the E-values, the links and the reap of the merged list with
//!   the whole sequence's length (`Runtime::subject`), as NCBI does after the chunks. Every chunk
//!   is walked, whatever its length (masks may overlap or nest, or be out of the order of their
//!   start; series NP, M4b-3b): a chunk shorter than the word gives the two-hit search its one word
//!   (`common::hard_chunk_scan::OneWord`), whose seeds go through the chunk's gapped stage or are
//!   its ungapped HSPs, and the one-hit search nothing, and the diagonal offset moves by each
//!   chunk's length plus the window.
//!
//! The masks come from `SeqDb::read_mask_block` a window of sequences at a time
//! (`algorithm::common::db_mask_search::SubjectMasks`, one reader per worker's diagonal table); a
//! read error is the reader's error of the search (`x_db_stream::XReadError`). LOSAT's refusals
//! of a sequence whose short chunk's HSPs NCBI's interval tree splits without end or its merge
//! joins with another chunk's are the search's errors of that sequence
//! (`common::hard_chunk_scan::short_chunk_tree_refusal`, `short_chunk_merge_refusal`, passed on
//! as they are by `native::search_error_as_is`).

use anyhow::Result;

use super::super::args::ResolvedOptions;
use super::super::parameters::ContextParameters;
use super::super::preliminary::{sort_by_score, PreliminaryHsp};
use super::super::query_setup::PreparedQueryBatch;
use super::super::results::subject_besthit_preliminary;
use super::super::seed::{sort_initial_hits, BxLookup, Diagonals, InitHsp};
use super::super::subjects::Subjects;
use super::super::x_db_stream::XReadError;
use super::{compute_subject, compute_subject_of, one_word_hits, SubjectComputation, SubjectPiece};
use crate::algorithm::common::database::search_read_error_of;
use crate::algorithm::common::db_mask_search::{SearchMask, SubjectPlan};
use crate::algorithm::common::hard_chunk_scan::{
    merge_joins, short_chunk_merge_refusal, short_chunk_tree_refusal, OneWord, WordContext,
};
use crate::algorithm::common::hard_chunk_walk::{ChunkWork, DiagonalWalk, HardTable, NcbiStart};
use crate::core::subject_split::HardChunk;
use crate::seqdb::{SeqDb, SeqDbError};

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
    let plain = |diagonals: &mut Diagonals| {
        masked_piece(
            subjects,
            batch,
            params,
            options,
            lookup,
            diagonals,
            xdrop,
            oid,
            x_containment,
            None,
        )
    };
    let Some(db) = subjects.seqdb() else {
        return plain(diagonals);
    };
    if subjects.length(oid) == 0 {
        return plain(diagonals);
    }
    let mut reader = diagonals.x_db_masks.take().unwrap_or_default();
    let stored = reader
        .of(db, oid as u32, mask.algo_id, mask_read_error, PROGRAM)
        .map(<[(u32, u32)]>::to_vec);
    diagonals.x_db_masks = Some(reader);
    let stored = stored?;
    // With `-db_hard_mask`: NCBI's one-thread offset where the search starts the sequence on a new
    // table (`ResolvedMask::hard_starts`, the walk that `check_losat_limits` ran).
    let ncbi = if mask.hard() {
        options.db_mask.as_ref().and_then(|resolved| {
            resolved
                .hard_starts(db, options.word_size, options.window_size)
                .of(oid as u32)
        })
    } else {
        None
    };
    if stored.is_empty() {
        let Some(start) = ncbi else {
            return plain(diagonals);
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 87-107
        // A sequence without masks on a new table at NCBI's one-thread offset (`NcbiWalk`),
        // searched
        // as in a search without the options on that table; by the reference scan of its whole
        // range where the entries lie above the offset (the other scans take every entry to lie
        // below it). Its entries end at least the window below the offset after it, so the next
        // sequence starts as after any sequence without masks.
        diagonals.ncbi_table(start, options.window_size);
        if start.unwritten(options.window_size) <= start.offset {
            return plain(diagonals);
        }
        let x_buffer = if batch.x_ungap_sent {
            subjects.x_encoded_with_sentinels(oid)
        } else {
            subjects.encoded_with_sentinel(oid)
        };
        let whole = [(0, subjects.length(oid) as i32)];
        return masked_piece(
            subjects,
            batch,
            params,
            options,
            lookup,
            diagonals,
            xdrop,
            oid,
            x_containment,
            Some(SubjectPiece {
                buffer: &x_buffer,
                ranges: Some(&whole),
            }),
        );
    }
    Masked {
        db,
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
        ncbi,
    }
    .compute(diagonals)
}

/// The search of database sequence `oid`, or with `piece` of a piece of it, by the search's
/// matrix: the default path's bodies (`compute_subject`, `compute_subject_of`).
#[allow(clippy::too_many_arguments)]
fn masked_piece(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    xdrop: i32,
    oid: usize,
    x_containment: bool,
    piece: Option<SubjectPiece<'_>>,
) -> Result<SubjectComputation> {
    if options.scoring == crate::config::ScoringMatrix::Blosum62 {
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
            piece,
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
            piece,
        )
    }
}

/// A masked sequence's search: the arguments of `compute_subject_with`, its database and its
/// stored masks.
struct Masked<'a> {
    db: &'a SeqDb,
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
    /// With `-db_hard_mask`, NCBI's one-thread start of the sequence where it is searched on a new
    /// table (`common::hard_chunk_walk::NcbiWalk`).
    ncbi: Option<NcbiStart>,
}

impl Masked<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:476-591
    // Each hard chunk is searched alone: its scan, its extensions, the purge of the HSPs with
    // common endpoints (gapped) and the sort by score; its HSPs move by the chunk's offset
    // (`Blast_HSPListAdjustOffsets`) and join the sequence's list (`Blast_HSPListsMerge`). NCBI's
    // overlap there is not 0: it compares the chunk's offset with the left of the range after it
    // (the chunker has already moved to it, blast_engine.c:254-258, 579-580) and takes the chunk
    // overlap of 100. The merge only looks at the HSPs of the joined list that end past the new
    // chunk's offset (blast_hits.c:2888-2899); a chunk ends where its range does and no HSP reaches
    // into the mask after it, and no chunk that can hold HSPs starts before the end of an earlier
    // one (else `HardMaskShape::PartsOutOfOrder`), so there are none and the merge appends the
    // chunk's HSPs and sorts the list by score (`s_BlastHSPListsCombineByScore`,
    // blast_hits.c:2749-2766). TBLASTN's and TBLASTX's chunks of a long range overlap; this reason
    // does not carry over. With `-subject_besthit` the filter then runs on the joined list after
    // every chunk with HSPs (blast_engine.c:586-591, blast_hits.c:2538-2580), not once on the whole
    // sequence's list.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
    // Every chunk goes through the word finder, whose exit moves the diagonal offset by the chunk's
    // length (negative between overlapping masks) plus the window.
    /// The sequence's computation: the link cutoffs of its whole length (`compute_subject_with`),
    /// and its HSPs over the soft ranges or the hard chunks, each piece searched by the default
    /// path's bodies (`masked_piece`). Every hard chunk is walked in order
    /// (`common::hard_chunk_walk::DiagonalWalk`): a chunk that holds a word is searched over its
    /// own letters (over the scan range of the whole chunk, the reference scan, where its diagonal
    /// offset went back below an entry of the table); a shorter one gives the two-hit search its
    /// one word ([`Masked::one_word`]), whose seeds go through the chunk's gapped stage, or are its
    /// HSPs in an ungapped search (LOSAT's refusal of the sequence where NCBI's interval tree
    /// splits without end or its merge joins them with another chunk's), and the one-hit search
    /// nothing. The sequence is searched on its own: on a new table at NCBI's one-thread offset
    /// where its hits meet the entries that no sequence wrote in NCBI's one-thread search or that
    /// search clears the table inside it (`ncbi`, `common::hard_chunk_walk::NcbiWalk`), else with
    /// the diagonal offset moved so that no entry of another sequence reaches its hits or theirs
    /// (`common::hard_chunk_walk::DiagonalWalk::start`,
    /// `DiagonalWalk::after`; PD-LOSAT-NCBI-DEFECTS exception 23, the Owner's answer Q5 = B of M4b,
    /// 2026-11-22). The sequences that LOSAT refuses before the search never come here
    /// (`blastinput::db_mask::HardMaskShape`: masks out of the order of their start, past the end,
    /// a range longer than `MAX_DBSEQ_LEN`, a tail mask with masks stored after it, a one word past
    /// the sentinel that can change a result), so the sequence's result does not depend on the
    /// entries of the sequences the table searched before or after it.
    fn compute(&self, diagonals: &mut Diagonals) -> Result<SubjectComputation> {
        let Masked {
            db,
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
            ncbi,
        } = *self;
        let piece = |diagonals: &mut Diagonals, piece: SubjectPiece<'_>| {
            masked_piece(
                subjects,
                batch,
                params,
                options,
                lookup,
                diagonals,
                xdrop,
                oid,
                x_containment,
                Some(piece),
            )
        };
        // The sequence with the sentinel after its residues (and before them too for a batch with
        // `x_ungap_sent`), as `compute_subject_with` reads it.
        let x_lead = batch.x_ungap_sent;
        let x_buffer = if x_lead {
            subjects.x_encoded_with_sentinels(oid)
        } else {
            subjects.encoded_with_sentinel(oid)
        };
        let residues = &x_buffer[usize::from(x_lead)..x_buffer.len() - 1];
        let length = residues.len();
        let mut plan = SubjectPlan::default();
        if !hard {
            // The whole sequence over its soft ranges: its own link cutoffs.
            let ranges = plan.soft(stored, length);
            return piece(
                diagonals,
                SubjectPiece {
                    buffer: &x_buffer,
                    ranges: Some(ranges),
                },
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1451-1455 (cited at
        // `compute_subject_with`): the link cutoffs of the sequence's whole length.
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
        plan.start_hard(stored, length);
        let mut chunks = Vec::new();
        plan.hard_chunks(&mut chunks);
        // The search's walk takes a one word that covers the sentinel to store entries (its lookup
        // finds no hit: BLASTX's translation has no residue 0), so `reference` is exact or more,
        // and its margins too (a larger margin moves the offset further past entries that no hit
        // reaches).
        let (mut walk, start) = DiagonalWalk::start(
            options.word_size,
            options.window_size,
            length,
            &chunks,
            ncbi,
            diagonals.offset(),
        );
        start.apply(diagonals, options.window_size)?;
        let mut raw = Vec::new();
        let mut purged = Vec::new();
        let mut containment = Vec::new();
        let mut letters = Vec::new();
        for &chunk in &chunks {
            // The exit of this chunk's word finder clears the table (`Diagonals::finish_chunk`).
            let clears = diagonals.offset() >= i32::MAX / 4;
            let work = walk.chunk(chunk.offset, chunk.length);
            let found = match work {
                ChunkWork::Words { reference } => {
                    // The chunk's letters between the sentinels that the sequence has.
                    letters.clear();
                    if x_lead {
                        letters.push(0);
                    }
                    letters
                        .extend_from_slice(&residues[chunk.offset..chunk.offset + chunk.letters()]);
                    letters.push(0);
                    let whole = [(0, chunk.length)];
                    let SubjectComputation {
                        raw: chunk_raw,
                        purged: chunk_purged,
                        containment: chunk_containment,
                        ..
                    } = piece(
                        diagonals,
                        SubjectPiece {
                            buffer: &letters,
                            ranges: reference.then_some(&whole[..]),
                        },
                    )?;
                    Some((chunk_raw, chunk_purged, chunk_containment))
                }
                ChunkWork::OneWord => {
                    let seeds = self.one_word(diagonals, residues, chunk, &mut letters);
                    diagonals.finish_chunk(chunk.length)?;
                    if seeds.is_empty() {
                        None
                    } else {
                        let mut chunk_containment = Vec::new();
                        let Some((chunk_raw, chunk_purged)) = one_word_hits(
                            batch,
                            params,
                            options,
                            &letters,
                            chunk.length,
                            &seeds,
                            xdrop,
                            x_containment,
                            &mut chunk_containment,
                        )?
                        else {
                            return Err(short_chunk_tree_refusal(
                                db.oid(oid as u32),
                                PROGRAM,
                                chunk.length,
                            ));
                        };
                        Some((chunk_raw, chunk_purged, chunk_containment))
                    }
                }
                ChunkWork::Nothing => {
                    diagonals.finish_chunk(chunk.length)?;
                    None
                }
            };
            if clears {
                walk.clear();
            }
            let Some((chunk_raw, chunk_purged, chunk_containment)) = found else {
                continue;
            };
            let offset = chunk.offset as i32;
            let chunk_hits = !chunk_purged.is_empty();
            // Only the HSPs of a chunk shorter than the word can meet another chunk's
            // (`HardMaskShape::PartsOutOfOrder` refuses chunks that hold words and meet).
            if chunk_hits
                && merge_joins(
                    &chunk,
                    purged.iter().map(|hsp: &PreliminaryHsp| hsp.s_end),
                    chunk_purged.iter().map(|hsp| hsp.s_start + offset),
                )
            {
                return Err(short_chunk_merge_refusal(db.oid(oid as u32), PROGRAM));
            }
            raw.extend(chunk_raw.into_iter().map(|hsp| shifted(hsp, offset)));
            purged.extend(chunk_purged.into_iter().map(|hsp| shifted(hsp, offset)));
            containment.extend(
                chunk_containment
                    .into_iter()
                    .map(|(hsp, contained)| (shifted(hsp, offset), contained)),
            );
            if options.subject_besthit && chunk_hits {
                // The filter of the joined list after this chunk, in the merge's order.
                sort_by_score(&mut purged);
                subject_besthit_preliminary(&mut purged);
            }
        }
        walk.after(diagonals.offset())
            .apply(diagonals, options.window_size)?;
        // The merge's order: by score.
        sort_by_score(&mut purged);
        Ok(SubjectComputation {
            raw,
            purged,
            containment,
            prepared_link,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-505
    // The two-hit scan of a chunk shorter than the word reads the one word at the chunk's offset,
    // from the sequence's stored letters (the preliminary stage's fetch is not masked).
    /// The one word of hard chunk `chunk`, shorter than the word, in a two-hit search
    /// (`common::hard_chunk_scan::OneWord`): the stored letters from the chunk's offset
    /// (`residues`, the
    /// sequence's; zeros past its end, which only a word whose hits change nothing reads:
    /// `HardMaskShape::WordPastEnd`), their lookup hits by the reference scan of one word, and
    /// NCBI's two-hit test of each on the diagonal table with the contexts' X-drops and cutoffs.
    /// Returns the seeds that it saved, in the order of the gapped stage (`seed::
    /// sort_initial_hits`), with the chunk's subject offsets; `letters` keeps the word's letters
    /// and the sentinel after them, the subject of the chunk's gapped stage.
    fn one_word(
        &self,
        diagonals: &mut Diagonals,
        residues: &[u8],
        chunk: HardChunk,
        letters: &mut Vec<u8>,
    ) -> Vec<InitHsp> {
        let Masked {
            batch,
            params,
            options,
            lookup,
            ..
        } = *self;
        let word = options.word_size;
        let size = usize::try_from(word).unwrap_or(0);
        // The word and the sentinel after it, as the scan reads a subject.
        letters.clear();
        letters
            .extend_from_slice(&residues[chunk.offset..(chunk.offset + size).min(residues.len())]);
        letters.resize(size + 1, 0);
        let mut hits = Vec::new();
        let range = [(0, word)];
        let mut each = |pairs: &[(i32, i32)]| hits.extend(pairs.iter().map(|&(q, _)| q));
        match lookup {
            BxLookup::Wide(lookup) => super::super::x_soft_subject::x_scan(
                lookup,
                letters,
                &range,
                options.window_size,
                &mut each,
            ),
            BxLookup::Narrow(lookup) => {
                lookup.x_scan(letters, size, Some(&range), options.window_size, &mut each)
            }
        }
        let test = OneWord {
            query: &batch.sequence_start[1..],
            letters: &letters[..size],
            length: chunk.length,
            word,
            window: options.window_size,
            scores: crate::utils::matrix::protein_score_rows(options.scoring),
        };
        let context = |q: i32| {
            let index = batch
                .contexts
                .partition_point(|c| c.offset as i32 <= q)
                .saturating_sub(1);
            WordContext {
                start: batch.contexts[index].offset as i32,
                x_drop: params[index].word_xdrop,
                cutoff: params[index].word_cutoff,
            }
        };
        let (entries, mask, offset) = diagonals.table();
        let mut seeds = Vec::new();
        test.test(&hits, entries, mask, offset, context, |q, seed| {
            seeds.push(InitHsp {
                q_seed: q,
                s_seed: 0,
                q_start: seed.q_start,
                s_start: seed.s_start,
                length: seed.length,
                score: seed.score,
            })
        });
        sort_initial_hits(&mut seeds);
        seeds
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
