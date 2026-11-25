//! The one word that NCBI's two-hit scan reads in a hard chunk shorter than the word, in a
//! database sequence with `-db_hard_mask` whose masks overlap or nest (BLASTP and BLASTX, series
//! NP M4b-3b; TBLASTN and TBLASTX, M4b-3c part 2; the walk of the chunks' diagonal offsets is
//! `algorithm::common::hard_chunk_walk::DiagonalWalk`). Here is what those programs share:
//!
//! - [`OneWord`]: the two-hit test of that word, which saves its seeds as NCBI's word finder does
//!   (the chunk's ungapped hits, from the stored letters past the chunk), which then go through
//!   the chunk's gapped stage, or are its HSPs in an ungapped search, as any chunk's seeds;
//! - [`merge_joins`]: whether NCBI's merge of a chunk's HSPs into the sequence's list joins
//!   them, which only such a chunk's HSPs can meet;
//! - [`short_chunk_tree_refusal`] and [`short_chunk_merge_refusal`]: LOSAT's refusals of a
//!   sequence where NCBI's gapped stage of such a chunk never ends, and where its merge would
//!   join such a chunk's HSPs with another's across an overlap read past the list (or in BLASTP
//!   and BLASTX, not yet ported).

use super::diag_entry::{
    batch_context_distance, x_test_hit, DiagCell, DiagRing, DiagStruct, DiagTable,
};
use crate::utils::matrix::ScoreRows;

/// The scores and the cutoffs of the context of a hit's query offset, which the two-hit test of
/// [`OneWord`] reads where a hit needs the extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct WordContext {
    /// The context's first offset in the concatenated query.
    pub(crate) start: i32,
    /// The X-drop of the ungapped extension.
    pub(crate) x_drop: i32,
    /// The score at which the extension saves a seed.
    pub(crate) cutoff: i32,
}

/// The result of the two-hit extension of one hit of [`OneWord`] (`s_BlastAaExtendTwoHit`): the
/// seed's score, whether the extension went right, the last subject offset that its right walk
/// read, and the ungapped hit (its first query offset in the concatenated query, its first subject
/// offset in the chunk and its length), which the word finder saves where the score reaches the
/// cutoff.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OneWordSeed {
    pub(crate) score: i32,
    pub(crate) right_extend: bool,
    pub(crate) s_last_off: i32,
    pub(crate) q_start: i32,
    pub(crate) s_start: i32,
    pub(crate) length: i32,
}

/// The one word that NCBI's two-hit scan reads in a hard chunk shorter than the word
/// (`common::hard_chunk_walk::ChunkWork::OneWord`), at subject offset 0 of the chunk.
#[derive(Clone, Copy)]
pub(crate) struct OneWord<'a> {
    /// The concatenated query as the program's two-hit extension reads it: its length bounds the
    /// extension's right walk.
    pub(crate) query: &'a [u8],
    /// The stored letters from the chunk's offset on, `word` of them (past the chunk, into the
    /// mask after it; after the sequence's sentinel, zeros).
    pub(crate) letters: &'a [u8],
    /// The chunk's length (NCBI's `Int4`, below the word): the bound of the right walk.
    pub(crate) length: i32,
    pub(crate) word: i32,
    pub(crate) window: i32,
    /// The search's matrix.
    pub(crate) scores: &'a ScoreRows,
}

impl OneWord<'_> {
    /// NCBI's two-hit test of each hit of the word (`hits`: the query offsets of its lookup
    /// cell, at subject offset 0) on the diagonal table `entries` (its mask `mask`, the offset
    /// `offset`) of one query batch ([`OneWord::test_on`]).
    pub(crate) fn test(
        &self,
        hits: &[i32],
        entries: &mut [DiagStruct],
        mask: i32,
        offset: i32,
        context: impl Fn(i32) -> WordContext,
        save: impl FnMut(i32, OneWordSeed),
    ) {
        self.test_on(&mut DiagRing { entries, mask }, hits, offset, context, save);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-510, 519-606
    // The scan's query offset is unsigned and counts from the start of the hit's query batch: the
    // context test subtracts the distance to the last hit in unsigned arithmetic, so a distance
    // past the batch's start passes it. A seed whose score reaches the cutoff is saved
    // (`BlastSaveInitHsp`) and the test goes on with the next hit; the entry is updated after
    // every extension.
    /// NCBI's two-hit test of each hit of the word (`hits`: the query offsets of its lookup
    /// cell, at subject offset 0) on the diagonal table `table` (the entry of each hit's diagonal
    /// in its query batch, `common::diag_entry::DiagTable`) at the offset `offset`, with
    /// `context` giving a query offset's context: the diagonal test, the context test in the
    /// offsets of the hit's batch, the extension and the store of each hit. `save` takes each
    /// seed whose score reaches the context's cutoff, with its hit's query offset, in the order
    /// of `hits`.
    pub(crate) fn test_on(
        &self,
        table: &mut impl DiagTable,
        hits: &[i32],
        offset: i32,
        context: impl Fn(i32) -> WordContext,
        save: impl FnMut(i32, OneWordSeed),
    ) {
        self.test_on_with(
            table,
            hits,
            offset,
            self.query.len() as i32,
            |q, s| self.scores.score(self.query[q as usize], s),
            context,
            save,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-510, 519-606
    /// [`OneWord::test_on`] where `score_at` scores concatenated query offset `q`
    /// (bounded by `query_len`) against a subject residue: used where the query is held
    /// per context (TBLASTX) or scored from a PSSM's rows (TBLASTN `-in_pssm`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn test_on_with(
        &self,
        table: &mut impl DiagTable,
        hits: &[i32],
        offset: i32,
        query_len: i32,
        score_at: impl Fn(i32, u8) -> i32,
        context: impl Fn(i32) -> WordContext,
        mut save: impl FnMut(i32, OneWordSeed),
    ) {
        for &q in hits {
            let (entry, batch_start) = table.two_hit_entry(q, 0);
            let Some((last_hit, diff)) = x_test_hit(entry, 0, offset, self.window, self.word)
            else {
                continue;
            };
            let found = context(q);
            let diff = batch_context_distance(q, diff, batch_start);
            if (q.wrapping_sub(diff) as u32) < found.start as u32 {
                entry.set_last_hit(offset);
                continue;
            }
            let seed = self.extend_with(
                last_hit.wrapping_add(self.word),
                q,
                found.x_drop,
                query_len,
                &score_at,
            );
            if seed.score >= found.cutoff {
                save(q, seed);
            }
            if seed.right_extend {
                entry.set_flag(1);
                entry.set_last_hit(
                    seed.s_last_off
                        .wrapping_sub(self.word - 1)
                        .wrapping_add(offset),
                );
            } else {
                entry.set_last_hit(offset);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:886-921
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:831-868
    // The best end of the word's first walk, then the left walk back from it, bounded by the
    // chunk's start and the query's (none where the bound is negative), then the right walk only
    // where the left walk reached the first hit's end (`s_left_off`, negative where that hit lies
    // before the chunk), bounded by the chunk's length (no letters where the length is at most the
    // walk's start) and the query's end. The ungapped hit runs from the left walk's best start to
    // the right walk's best end (one past its last letter that raised the score).
    /// The two-hit extension (`s_BlastAaExtendTwoHit`) of the hit at query offset `q_right_off`
    /// and subject offset 0 whose first hit ends at `s_left_off` (signed), with the X-drop
    /// `x_drop`.
    #[cfg(test)]
    pub(crate) fn extend(&self, s_left_off: i32, q_right_off: i32, x_drop: i32) -> OneWordSeed {
        self.extend_with(
            s_left_off,
            q_right_off,
            x_drop,
            self.query.len() as i32,
            &|q, s| self.scores.score(self.query[q as usize], s),
        )
    }

    /// The two-hit extension (`s_BlastAaExtendTwoHit`) of the hit at query offset `q_right_off`
    /// and subject offset 0 whose first hit ends at `s_left_off` (signed), with the X-drop
    /// `x_drop`, the query's length `query_len` and the residue scorer `score_at`.
    fn extend_with(
        &self,
        s_left_off: i32,
        q_right_off: i32,
        x_drop: i32,
        query_len: i32,
        score_at: &impl Fn(i32, u8) -> i32,
    ) -> OneWordSeed {
        let score_at = |q: i32, s: i32| score_at(q, self.letters[s as usize]);
        let (mut left_score, mut right_d, mut score) = (0, 0, 0);
        for i in 0..self.word {
            score += score_at(q_right_off + i, i);
            if score > left_score {
                left_score = score;
                right_d = i + 1;
            }
        }
        let (q_right_off, s_right_off) = (q_right_off + right_d, right_d);
        // The left walk from one before the best end.
        let (s_off, q_off) = (s_right_off - 1, q_right_off - 1);
        let n = s_off.min(q_off);
        let (mut best_i, mut walk, mut best) = (n + 1, 0, 0);
        let mut i = n;
        while i >= 0 {
            walk += score_at(q_off - n + i, s_off - n + i);
            if walk > best {
                best = walk;
                best_i = i;
            }
            if best - walk >= x_drop {
                break;
            }
            i -= 1;
        }
        let (left_score, left_d) = (best, n - best_i + 1);
        let (mut right_score, mut s_last_off, mut right_len) = (0, s_right_off, 0);
        let right_extend = left_d >= s_right_off - s_left_off;
        if right_extend {
            let n = (self.length - s_right_off).min(query_len - q_right_off);
            let (mut walk, mut best, mut best_i) = (left_score, left_score, -1);
            let mut i = 0;
            while i < n {
                walk += score_at(q_right_off + i, s_right_off + i);
                if walk > best {
                    best = walk;
                    best_i = i;
                }
                if walk <= 0 || best - walk >= x_drop {
                    break;
                }
                i += 1;
            }
            (right_score, s_last_off, right_len) = (best, s_right_off + i, best_i + 1);
        }
        OneWordSeed {
            score: left_score.max(right_score),
            right_extend,
            s_last_off,
            q_start: q_right_off - left_d,
            s_start: s_right_off - left_d,
            length: left_d + right_len,
        }
    }
}

/// LOSAT's refusals of a sequence during its search, where the one word of a hard chunk shorter
/// than the word ([`OneWord`]) saved seeds ([`short_chunk_tree_refusal`],
/// [`short_chunk_merge_refusal`]), a type of their own so that a program that reports the other
/// errors of its search as NCBI's engine errors passes them on as they are (BLASTX's
/// `native::search_error_as_is`).
#[derive(Debug)]
pub(crate) struct ShortChunkSeedRefusal(String);

impl std::fmt::Display for ShortChunkSeedRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ShortChunkSeedRefusal {}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-545
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3820-3831
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:511-790
// The seeds of such a word go to the chunk's gapped stage, whose subject length is the chunk's:
// the interval tree of the stage's HSPs ends one past it (before its start where the length is
// negative), and the gapped extension of a seed reads the stored letters past the chunk to the
// left of its start, so an HSP can leave the tree's bounds, where NCBI's insertion can split a
// node without end (`BlastIntervalTree::try_add_hsp`): NCBI's search then never ends (its memory
// grows until it fails; series NP, M4b-3b probe 2, g3: blastp's default search and blastx's at
// -evalue 1e6). Elsewhere LOSAT runs the stage as NCBI does (`blastp::blast_engine::
// prelim_kernel::BlastpPrelimKernel::one_word_gapped`, `blastx::preliminary::gapped_observed`).
/// LOSAT's refusal (in `program`'s name) of the database sequence of OID `oid`, where the one word
/// of a part shorter than the word after overlapping masks, of `length` letters (NCBI's `Int4`:
/// negative where the masks overlap), saves seeds whose gapped alignments NCBI's interval tree
/// splits without end.
pub(crate) fn short_chunk_tree_refusal(
    oid: impl std::fmt::Display,
    program: &str,
    length: i32,
) -> anyhow::Error {
    ShortChunkSeedRefusal(format!(
        "-db_hard_mask: the sequence of OID {oid} has a part shorter than the word after \
         overlapping masks, of {length} letters to NCBI BLAST+ 2.17.0, which saves seeds from the \
         one word that it reads past that part and adds their alignments to a tree whose bounds \
         that length sets, where it splits the tree's nodes without end (its search never ends); \
         LOSAT's {program} does not search such a database with -db_hard_mask"
    ))
    .into()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:576-584
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2873-2911, 2956-2998
// After each chunk with HSPs NCBI merges them into the sequence's list with the subject split at
// the chunk's offset: the joined list's HSPs that end past the offset and the chunk's that start
// before the offset plus the overlap (0 where the chunk starts at the left of the range after its
// own, else 100; after the last range NCBI reads it past its list) go to its loop of joins, which
// LOSAT's protein programs do not port (with `HardMaskShape::LongRange` and
// `HardMaskShape::PartsOutOfOrder`, item M4b-3d of series NP). Where no HSP of one side is there,
// the merge appends the chunk's HSPs, which LOSAT does.
/// Whether NCBI's merge of the HSPs of hard chunk `chunk` into a sequence's list runs its loop of
/// joins: an HSP of the list ends past the chunk's offset (`earlier_ends`, the ends of the list's
/// HSPs in the sequence) and one of the chunk's starts before its offset plus the overlap
/// (`starts`, the starts of the chunk's HSPs in the sequence; any where NCBI reads the overlap
/// past its list of ranges).
pub(crate) fn merge_joins(
    chunk: &crate::core::subject_split::HardChunk,
    mut earlier_ends: impl Iterator<Item = i32>,
    mut starts: impl Iterator<Item = i32>,
) -> bool {
    let offset = i64::try_from(chunk.offset).unwrap_or(i64::MAX);
    let overlap = if chunk.overlap_past_list {
        i64::MAX / 2
    } else {
        i64::try_from(chunk.overlap).unwrap_or(i64::MAX / 2)
    };
    earlier_ends.any(|end| i64::from(end) > offset)
        && starts.any(|start| i64::from(start) < offset + overlap)
}

/// LOSAT's refusal (in `program`'s name) of the database sequence of OID `oid`, where the HSPs of a
/// part shorter than the word after overlapping masks meet those of another part in NCBI's merge
/// across their overlap ([`merge_joins`]), which LOSAT does not do yet.
pub(crate) fn short_chunk_merge_refusal(
    oid: impl std::fmt::Display,
    program: &str,
) -> anyhow::Error {
    ShortChunkSeedRefusal(format!(
        "-db_hard_mask: the sequence of OID {oid} has a part shorter than the word after \
         overlapping masks, whose alignments from the one word that NCBI BLAST+ 2.17.0 reads past \
         that part meet those of another part in its merge across their overlap, which LOSAT's \
         {program} does not do yet; LOSAT's {program} does not search such a database with \
         -db_hard_mask"
    ))
    .into()
}

#[cfg(test)]
#[path = "hard_chunk_scan_tests.rs"]
mod tests;
