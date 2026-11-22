//! The word finders over the hard chunks of a protein database sequence with `-db_hard_mask`
//! (BLASTP and BLASTX; series NP, M4b-3b), where the stored masks, sorted by their start, may
//! overlap or nest. NCBI searches every chunk of the sequence's hard ranges whatever its length
//! (`core::subject_split::HardChunk::length`, NCBI's `Int4`: negative for the inverted range
//! between two masks that overlap), and after each chunk moves the diagonal offset by the chunk's
//! length plus the window: back where the length is negative, so that a later chunk can meet the
//! word hits of an earlier one on its diagonals. Here is what both programs share:
//!
//! - [`DiagonalWalk`]: the diagonal offsets of a sequence's chunks, relative to the offset before
//!   its first chunk, what the word finder reads of each chunk ([`ChunkWork`]), and the margins by
//!   which the sequence's word hits reach past it ([`DiagonalMargins`]), which LOSAT's refusal
//!   before the search reads (`blastinput::db_mask::HardMaskShape::StaleDiagonals`);
//! - [`OneWord`]: the two-hit test of the one word that NCBI's two-hit scan reads in a chunk
//!   shorter than the word;
//! - [`short_chunk_seed_refusal`]: LOSAT's refusal of a sequence where that word saves a seed.

use super::x_seed_fused::x_test_hit;
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::utils::matrix::ScoreRows;

/// What the word finder reads of one hard chunk ([`DiagonalWalk::chunk`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ChunkWork {
    /// The chunk holds a word: the scan reads its words, from its offset to its length less the
    /// word. `reference`: the chunk's diagonal offset lies below an entry that an earlier chunk of
    /// the sequence may have stored (the offset went back), so the chunk is searched by the
    /// reference scan (the scan ranges of `-db_soft_mask`, one range of the whole chunk): the other
    /// scans take every entry of the table to lie below the offset of the subject they read
    /// (`blastp::x_seed_fused::XHitBf`, blastx's `x_seed_scan.rs`). Never with masks that do not
    /// overlap, where the offset only grows.
    Words { reference: bool },
    /// A chunk of a two-hit search shorter than the word (of no letters, or of a negative length,
    /// too): the scan reads the one word at the chunk's offset, past the chunk ([`OneWord`]).
    OneWord,
    /// A chunk of a one-hit search shorter than the word: its scan range is empty, nothing is read.
    Nothing,
}

/// How far a sequence's word hits reach beyond the sequence on the diagonal table
/// ([`DiagonalWalk::margins`]), in positions of the diagonal offset. Where both are 0 the table
/// entries that the sequence's chunks store read as cleared entries to every hit of a later
/// sequence, and its own hits read the entries of earlier sequences as cleared ones: the sequence's
/// result does not depend on the sequences that the same table searched before or after it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct DiagonalMargins {
    /// How far the lowest diagonal offset of a chunk whose scan reads a word lies below the
    /// sequence's offset before its first chunk: its hits can meet the entries of earlier
    /// sequences (0 where none does).
    pub(crate) below_start: i64,
    /// How far an entry that the sequence's chunks can store lies past the last value that every
    /// hit of the next sequence reads as cleared: the offset after the last chunk less the window
    /// for an entry without the flag (a later hit is at least the window past it), the offset
    /// after the last chunk for one with the flag (a later hit is not behind it). 0 where none
    /// does.
    pub(crate) past_end: i64,
}

impl DiagonalMargins {
    /// Whether the sequence's hits reach entries of another sequence, or leave entries that
    /// another sequence's hits read (either margin above 0).
    pub(crate) fn stale(self) -> bool {
        self.below_start > 0 || self.past_end > 0
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-500
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-766
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551, 588-606
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:771-802
// The chunk's length is NCBI's `Int4` and the exit of each chunk's word finder adds it and the
// window to the diagonal offset (the window is 0 in the one-hit search). The two-hit scan raises
// the end of its range to its start, so it reads one word at the chunk's offset whatever the
// chunk's length; the one-hit scan of a chunk shorter than the word reads nothing. A scanned hit
// at subject offset `s` of a chunk at diagonal offset `D` stores `s + D` without the flag, or the
// extension's last offset less the word plus one and `D` with the flag (two-hit; the one-hit
// search stores that value too, and never sets the flag). The two-hit scan reads `s` up to the
// chunk's length less the word (0 for a shorter chunk), and the extension's last offset is at most
// the chunk's length, or the end of the one word's first walk in a shorter chunk; the one-hit
// scan reads `s` up to the length less the word, and the extension's last offset is at most the
// length. Hence the bounds of the stores below.
/// The diagonal offsets of the hard chunks of one sequence in a protein program's word finder,
/// relative to the offset before the sequence's first chunk (`D_0 = 0`, and `D_{i+1} = D_i + L_i
/// + W` for a chunk of `L_i` letters and the window `W`, 0 in the one-hit search): what the word
/// finder reads of each chunk ([`DiagonalWalk::chunk`]), the bounds of the entries that the
/// chunks can store, and from them the margins of [`DiagonalWalk::margins`]. The search and
/// LOSAT's refusal before it (`blastinput::db_mask::HardMaskShape::of_scan`) both walk the
/// chunks of `core::subject_split::HardSplitState` with it, so the two cannot disagree.
#[derive(Debug, Clone)]
pub(crate) struct DiagonalWalk {
    word: i64,
    window: i64,
    /// The sequence's letters: its sentinel byte is at this offset.
    length: i64,
    /// The diagonal offset of the next chunk.
    at: i64,
    /// The highest entry without the flag that a chunk so far can store (the one-hit search's
    /// entries too).
    unflagged: Option<i64>,
    /// The highest entry with the flag that a chunk so far can store (two-hit).
    flagged: Option<i64>,
    /// The lowest diagonal offset of a chunk whose scan reads a word.
    lowest_scan: Option<i64>,
    /// The offset of the latest chunk, and the highest offset of the chunks before it.
    latest: Option<i64>,
    inner: Option<i64>,
    /// The one words so far that run past the sequence's sentinel, in walk order: the diagonal
    /// offset of each, its chunk's offset in the sequence, and whether its hits can change a
    /// result ([`DiagonalWalk::word_past_end`]; empty but in rare sequences).
    past_end_words: Vec<(i64, usize, bool)>,
}

impl DiagonalWalk {
    /// A walk of the chunks of a sequence of `length` letters in a search with the word `word` and
    /// the window `window` (0: the one-hit search).
    pub(crate) fn new(word: i32, window: i32, length: usize) -> Self {
        Self {
            word: i64::from(word),
            window: i64::from(window),
            length: i64::try_from(length).unwrap_or(i64::MAX),
            at: 0,
            unflagged: None,
            flagged: None,
            lowest_scan: None,
            latest: None,
            inner: None,
            past_end_words: Vec::new(),
        }
    }

    /// The next chunk, at `offset` in the sequence with the length `length` (NCBI's `Int4`): what
    /// its word finder reads. The diagonal offset then moves on by the length and the window.
    pub(crate) fn chunk(&mut self, offset: usize, length: i32) -> ChunkWork {
        let (d, l, w, window) = (self.at, i64::from(length), self.word, self.window);
        self.inner = self.inner.max(self.latest);
        self.latest = Some(d);
        let work = if l >= w {
            let stored = self.unflagged.max(self.flagged);
            ChunkWork::Words {
                reference: stored.is_some_and(|stored| stored > d),
            }
        } else if window > 0 {
            ChunkWork::OneWord
        } else {
            ChunkWork::Nothing
        };
        if work != ChunkWork::Nothing {
            self.lowest_scan = Some(self.lowest_scan.map_or(d, |lowest| lowest.min(d)));
            self.past_words(offset, d, work == ChunkWork::OneWord);
            let (unflagged, flagged) = if window > 0 {
                (d + (l - w).max(0), Some(d + (l - w + 1).max(1)))
            } else {
                (d + l - w + 1, None)
            };
            self.unflagged = self.unflagged.max(Some(unflagged));
            self.flagged = self.flagged.max(flagged);
        }
        self.at = d + l + window;
        work
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:323-330, 345-378
    // The preliminary stage's fetch of a protein sequence is not copied: it points into the
    // mapped file, where the next sequence's letters follow the sentinel byte.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551
    // A hit on an entry without the flag that lies at least the window behind it, or on one with
    // the flag that is not ahead of it, stores the hit's own value, as on a cleared entry.
    /// The one word of a chunk at `offset` (`one_word`, at diagonal offset `d`) runs past the
    /// sequence's sentinel when `offset + word` exceeds the length by more than one: NCBI reads its
    /// last letters from the next sequence in the database, which LOSAT does not have (it pads the
    /// word with zeros). The word's hits cannot change a result when no entry without the flag that
    /// an earlier chunk can store lies within the window behind it (so no hit of it is extended)
    /// and every later chunk that reads a word lies at least the window past it (so a later hit on
    /// a diagonal it stored reads that entry as it would read the one before it); past the
    /// sequence the margins cover it (`past_end`). Otherwise the sequence is refused
    /// (`HardMaskShape::WordPastEnd`). Sorted masks that do not overlap meet both conditions.
    /// Each word's conditions are its own (they do not depend on which other words run past the
    /// sentinel), so that the refusal of a sequence and the offset it names change only where the
    /// set of such words does (`HardMaskShape::of_scan` reads the shape at two lengths).
    fn past_words(&mut self, offset: usize, d: i64, one_word: bool) {
        for word in &mut self.past_end_words {
            if d < word.0 + self.window {
                word.2 = true;
            }
        }
        let start = i64::try_from(offset).unwrap_or(i64::MAX);
        if one_word && start.saturating_add(self.word) > self.length + 1 {
            let extends = self
                .unflagged
                .is_some_and(|stored| stored > d - self.window);
            self.past_end_words.push((d, offset, extends));
        }
    }

    /// The margins of the chunks walked so far (after the last chunk, the sequence's): how far
    /// the lowest diagonal offset of a chunk that reads a word lies below the sequence's start,
    /// and how far its highest entry lies past what every hit of the next sequence reads as
    /// cleared (the offset after the last chunk, less the window for an entry without the flag).
    /// Where either is above 0, NCBI's rows of the sequence depend on the sequences that its
    /// thread searches before and after it, and LOSAT refuses it before the search
    /// (`HardMaskShape::StaleDiagonals`). A search that started the sequence's offset higher by
    /// the first margin and the next sequence's higher by the second would isolate it instead (the
    /// alternative of question Q5 of M/m4b/OWNER_QUESTIONS.md, not taken).
    pub(crate) fn margins(&self) -> DiagonalMargins {
        let end = self.at;
        let unflagged = self
            .unflagged
            .map_or(0, |stored| stored - (end - self.window));
        let flagged = self.flagged.map_or(0, |stored| stored - end);
        DiagonalMargins {
            below_start: self.lowest_scan.map_or(0, |lowest| (-lowest).max(0)),
            past_end: unflagged.max(flagged).max(0),
        }
    }

    /// The offset in the sequence of the first chunk whose one word runs past the sequence's
    /// sentinel where its hits can change a result ([`DiagonalWalk::chunk`]'s check), if any.
    pub(crate) fn word_past_end(&self) -> Option<usize> {
        self.past_end_words
            .iter()
            .find(|&&(_, _, changes)| changes)
            .map(|&(_, offset, _)| offset)
    }

    /// The diagonal offset after the chunks walked so far, relative to the sequence's start.
    #[cfg(test)]
    pub(crate) fn end(&self) -> i64 {
        self.at
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
    // The exit of each chunk's word finder clears the table and sets the offset to the window when
    // the offset of the chunk is at least INT4_MAX / 4. Where that falls inside a sequence depends
    // on the sequences that NCBI's thread searched before it (NCBI's own dependence on its
    // threads, at most once per about 10^6 sequences of a thread); LOSAT clears the table before
    // such a sequence instead, so that no clear falls between its chunks and its result does not
    // depend on the worker's history.
    /// Whether a search whose diagonal offset is `offset` before the sequence's first chunk would
    /// clear the table after a chunk that is not the last (`Blast_ExtendWordExit`): then the table
    /// is cleared before the first chunk, as after the sequence before it.
    pub(crate) fn clears_inside(&self, offset: i32) -> bool {
        self.inner
            .is_some_and(|inner| i64::from(offset) + inner >= i64::from(i32::MAX / 4))
    }
}

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

/// The result of the two-hit extension of one hit of [`OneWord`] (`s_BlastAaExtendTwoHit`) that
/// the two-hit test reads: the seed's score, whether the extension went right, and the last subject
/// offset that its right walk read. (Where NCBI saves the seed LOSAT refuses the sequence, so the
/// seed's start and length are not kept.)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OneWordSeed {
    pub(crate) score: i32,
    pub(crate) right_extend: bool,
    pub(crate) s_last_off: i32,
}

/// The one word that NCBI's two-hit scan reads in a hard chunk shorter than the word
/// ([`ChunkWork::OneWord`]), at subject offset 0 of the chunk.
#[derive(Debug, Clone, Copy)]
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-510, 519-606
    // The scan's query offset is unsigned: the context test subtracts the distance to the last
    // hit in unsigned arithmetic, so a distance past the query's start passes it.
    /// NCBI's two-hit test of each hit of the word (`hits`: the query offsets of its lookup
    /// cell, at subject offset 0) on the diagonal table `entries` (its mask `mask`, the offset
    /// `offset`), with `context` giving a query offset's context: the diagonal test, the context
    /// test, the extension and the store of each hit. Returns whether an extension saved a seed
    /// (its score at least the context's cutoff), where the search then stops.
    pub(crate) fn test(
        &self,
        hits: &[i32],
        entries: &mut [DiagStruct],
        mask: i32,
        offset: i32,
        context: impl Fn(i32) -> WordContext,
    ) -> bool {
        for &q in hits {
            let entry = &mut entries[(q & mask) as usize];
            let Some((last_hit, diff)) = x_test_hit(entry, 0, offset, self.window, self.word)
            else {
                continue;
            };
            let found = context(q);
            if (q.wrapping_sub(diff) as u32) < found.start as u32 {
                entry.set_last_hit(offset);
                continue;
            }
            let seed = self.extend(last_hit.wrapping_add(self.word), q, found.x_drop);
            if seed.score >= found.cutoff {
                return true;
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
        false
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:886-921
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:831-868
    // The best end of the word's first walk, then the left walk back from it, bounded by the
    // chunk's start and the query's (none where the bound is negative), then the right walk only
    // where the left walk reached the first hit's end (`s_left_off`, negative where that hit lies
    // before the chunk), bounded by the chunk's length (no letters where the length is at most the
    // walk's start) and the query's end.
    /// The two-hit extension (`s_BlastAaExtendTwoHit`) of the hit at query offset `q_right_off`
    /// and subject offset 0 whose first hit ends at `s_left_off` (signed), with the X-drop
    /// `x_drop`.
    pub(crate) fn extend(&self, s_left_off: i32, q_right_off: i32, x_drop: i32) -> OneWordSeed {
        let score_at = |q: i32, s: i32| {
            self.scores
                .score(self.query[q as usize], self.letters[s as usize])
        };
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
        let (mut right_score, mut s_last_off) = (0, s_right_off);
        let right_extend = left_d >= s_right_off - s_left_off;
        if right_extend {
            let n = (self.length - s_right_off).min(self.query.len() as i32 - q_right_off);
            let (mut walk, mut best) = (left_score, left_score);
            let mut i = 0;
            while i < n {
                walk += score_at(q_right_off + i, s_right_off + i);
                best = best.max(walk);
                if walk <= 0 || best - walk >= x_drop {
                    break;
                }
                i += 1;
            }
            (right_score, s_last_off) = (best, s_right_off + i);
        }
        OneWordSeed {
            score: left_score.max(right_score),
            right_extend,
            s_last_off,
        }
    }
}

/// LOSAT's refusal of a search in which the one word of a hard chunk shorter than the word
/// ([`OneWord`]) saves a seed, a type of its own so that a program that reports the other errors
/// of its search as NCBI's engine errors passes it on as it is (BLASTX's
/// `native::search_error_as_is`).
#[derive(Debug)]
pub(crate) struct ShortChunkSeedRefusal(String);

impl std::fmt::Display for ShortChunkSeedRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ShortChunkSeedRefusal {}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-500
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158
// The seed of such a word lies past the chunk, from which NCBI's gapped search of a chunk of a
// negative length never ends (its memory grows) and whose rows show letters past the chunk
// otherwise (series NP, M4b-3b probe 2: s2, g3).
/// LOSAT's refusal (in `program`'s name) of the search of the database sequence of OID `oid`, in
/// which the one word of a part between overlapping masks shorter than the word saves a seed.
pub(crate) fn short_chunk_seed_refusal(
    oid: impl std::fmt::Display,
    program: &str,
) -> anyhow::Error {
    ShortChunkSeedRefusal(format!(
        "-db_hard_mask: the sequence of OID {oid} has a part between overlapping masks shorter \
         than the word, in which NCBI BLAST+ 2.17.0 saves a seed from the one word that it reads \
         past that part (its rows then show letters past the part, or its search never ends); \
         LOSAT's {program} does not search such a database with -db_hard_mask"
    ))
    .into()
}

#[cfg(test)]
#[path = "hard_chunk_scan_tests.rs"]
mod tests;
