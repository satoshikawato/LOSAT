//! The walk of the diagonal offsets over the hard chunks of a protein database sequence with
//! `-db_hard_mask` (BLASTP and BLASTX; series NP, M4b-3b), where the stored masks may overlap or
//! nest, or be out of the order of their start. NCBI searches every chunk of the sequence's hard
//! ranges whatever its length (`core::subject_split::HardChunk::length`, NCBI's `Int4`: negative
//! for the inverted range between two masks that overlap), and after each chunk moves the diagonal
//! offset by the chunk's length plus the window: back where the length is negative, so that a
//! later chunk can meet the word hits of an earlier one on its diagonals. [`DiagonalWalk`] is the
//! one owner of that walk: the two programs' searches (`blastp::blast_engine::db_masked`,
//! `blastx::search::db_masked`) and LOSAT's refusal before the search
//! (`blastinput::db_mask::HardMaskShape::of_scan`) all walk the chunks of
//! `core::subject_split::HardSplitState` with it, so they cannot disagree. It gives what the word
//! finder reads of each chunk ([`ChunkWork`]) and the margins by which the sequence's word hits
//! reach past it ([`DiagonalMargins`]).

/// What the word finder reads of one hard chunk ([`DiagonalWalk::chunk`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ChunkWork {
    /// The chunk holds a word: the scan reads its words, from its offset to its length less the
    /// word. `reference`: the chunk's diagonal offset lies below an entry that the table may hold
    /// (one that an earlier chunk of the sequence stored, or a cleared one after a clear inside
    /// the sequence: the offset went back), so the chunk is searched by the reference scan (the
    /// scan ranges of `-db_soft_mask`, one range of the whole chunk): the other scans take every
    /// entry of the table to lie below the offset of the subject they read
    /// (`blastp::x_seed_fused::XHitBf`, blastx's `x_seed_scan.rs`). Never with masks apart in the
    /// order of their start, where the offset only grows.
    Words { reference: bool },
    /// A chunk of a two-hit search shorter than the word (of no letters, or of a negative length,
    /// too): the scan reads the one word at the chunk's offset, past the chunk
    /// (`blastp::hard_chunk_scan::OneWord`).
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
    /// How far the lowest diagonal offset of a chunk whose scan reads a word that can hit lies
    /// below the sequence's offset before its first chunk: its hits can meet the entries of
    /// earlier sequences (0 where none does).
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1018-1080
// The chunk's length is NCBI's `Int4` and the exit of each chunk's word finder adds it and the
// window to the diagonal offset (the window is 0 in the one-hit search). The two-hit scan raises
// the end of its range to its start, so it reads one word at the chunk's offset whatever the
// chunk's length; the one-hit scan of a chunk shorter than the word reads nothing. A scanned hit
// at subject offset `s` of a chunk at diagonal offset `D` stores `s + D` without the flag, or the
// extension's last offset less the word plus one and `D` with the flag (two-hit; the one-hit
// search stores that value too, and never sets the flag). The two-hit scan reads `s` up to the
// chunk's length less the word (0 for a shorter chunk), and the extension's last offset is at most
// the chunk's length, or the end of the one word's first walk in a shorter chunk. The one-hit
// scan reads `s` up to the length less the word, and its extension's last offset is at most the
// length, or one past it: where no letter of the word scores above 0 (an exact word of the query
// with a self-score of 0, which the lookup table always holds, blast_aalookup.c:503-509) the
// right walk starts one past the word's end and reads nothing. Hence the bounds of the stores
// below.
/// The diagonal offsets of the hard chunks of one sequence in a protein program's word finder,
/// relative to the offset before the sequence's first chunk (`D_0 = 0`, and `D_{i+1} = D_i + L_i
/// + W` for a chunk of `L_i` letters and the window `W`, 0 in the one-hit search): what the word
/// finder reads of each chunk ([`DiagonalWalk::chunk`]), the bounds of the entries that the
/// chunks can store, and from them the margins of [`DiagonalWalk::margins`].
#[derive(Debug, Clone)]
pub(crate) struct DiagonalWalk {
    word: i64,
    window: i64,
    /// The sequence's letters: its sentinel byte is at this offset.
    length: i64,
    /// Whether the query's lookup table can hold a word with the residue 0 (the gap `-`), which
    /// the sentinel byte after the sequence's letters is ([`DiagonalWalk::new`]).
    zero_words: bool,
    /// The diagonal offset of the next chunk.
    at: i64,
    /// The highest entry without the flag that the table can hold so far (the one-hit search's
    /// entries too).
    unflagged: Option<i64>,
    /// The highest entry with the flag that a chunk so far can store (two-hit).
    flagged: Option<i64>,
    /// The lowest diagonal offset of a chunk whose scan reads a word that can hit.
    lowest_scan: Option<i64>,
    /// The one words so far that run past the sequence's sentinel, in walk order: the diagonal
    /// offset of each, its chunk's offset in the sequence, and whether its hits can change a
    /// result ([`DiagonalWalk::word_past_end`]; empty but in rare sequences).
    past_end_words: Vec<(i64, usize, bool)>,
}

impl DiagonalWalk {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1578
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:478-509
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:79-125
    // Every score of the residue 0 (the gap `-`) is `BLAST_SCORE_MIN`, so no neighbouring word
    // that holds it reaches the threshold; the lookup table holds such a word only as an exact
    // word of the query (a query word whose self-score is below the threshold is added as it is),
    // so only where the query holds that residue.
    /// A walk of the chunks of a sequence of `length` letters in a search with the word `word` and
    /// the window `window` (0: the one-hit search). `zero_words`: whether the query's lookup table
    /// can hold a word with the residue 0, the sentinel byte after the sequence's letters (a
    /// query that holds the gap `-`); without such a word, the one word of a chunk that covers
    /// the sentinel has no hit, stores nothing and reads no letter past it that matters.
    pub(crate) fn new(word: i32, window: i32, length: usize, zero_words: bool) -> Self {
        Self {
            word: i64::from(word),
            window: i64::from(window),
            length: i64::try_from(length).unwrap_or(i64::MAX),
            zero_words,
            at: 0,
            unflagged: None,
            flagged: None,
            lowest_scan: None,
            past_end_words: Vec::new(),
        }
    }

    /// The next chunk, at `offset` in the sequence with the length `length` (NCBI's `Int4`): what
    /// its word finder reads. The diagonal offset then moves on by the length and the window.
    pub(crate) fn chunk(&mut self, offset: usize, length: i32) -> ChunkWork {
        let (d, l, w, window) = (self.at, i64::from(length), self.word, self.window);
        let start = i64::try_from(offset).unwrap_or(i64::MAX);
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
        // A one word that covers the sentinel byte (at `length`) without a lookup word that holds
        // the residue 0 has no hit.
        let hits = match work {
            ChunkWork::Nothing => false,
            ChunkWork::OneWord => self.zero_words || start.saturating_add(w) <= self.length,
            ChunkWork::Words { .. } => true,
        };
        if hits {
            self.lowest_scan = Some(self.lowest_scan.map_or(d, |lowest| lowest.min(d)));
            self.past_words(offset, d, work == ChunkWork::OneWord);
            let (unflagged, flagged) = if window > 0 {
                (d + (l - w).max(0), Some(d + (l - w + 1).max(1)))
            } else {
                (d + l - w + 2, None)
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
    /// The one word of a chunk at `offset` (`one_word`, at diagonal offset `d`, in a search whose
    /// query can hit a word with the residue 0) runs past the sequence's sentinel when `offset +
    /// word` exceeds the length by more than one: NCBI reads its last letters from the next
    /// sequence in the database, which LOSAT does not have (it pads the word with zeros). The
    /// word's hits cannot change a result when no entry without the flag that an earlier chunk can
    /// store lies within the window behind it (so no hit of it is extended) and every later chunk
    /// that reads a word lies at least the window past it (so a later hit on a diagonal it stored
    /// reads that entry as it would read the one before it); past the sequence the margins cover
    /// it (`past_end`). Otherwise the sequence is refused (`HardMaskShape::WordPastEnd`). Masks
    /// apart in the order of their start meet both conditions. Each word's conditions are its own
    /// (they do not depend on which other words run past the sentinel), so that the refusal of a
    /// sequence and the offset it names change only where the set of such words does
    /// (`HardMaskShape::of_scan` reads the shape at two lengths).
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

    /// The letters of a word.
    pub(crate) fn word(&self) -> i32 {
        i32::try_from(self.word).unwrap_or(i32::MAX)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 110-155
    /// The sequence is the first that the table searches: NCBI's thread makes its table with every
    /// entry 0 without the flag and the offset the window (`s_BlastDiagTableNew`, the `calloc` of
    /// `BlastExtendWordNew`), so the entries lie the window below the sequence's start, where a
    /// chunk that goes back below the start meets them (a hit within the window of them is
    /// extended, one in the one-hit search below them is skipped). The search of that sequence
    /// makes its table so first (`blastp::blast_engine::BlastpSubjectScratch::new_table`,
    /// `blastx::seed::Diagonals::new_table`), and the walk of the refusal before the search starts
    /// so for the first searched sequence (`HardMaskShape::of_scan`), whose margin below the start
    /// is then no refusal.
    pub(crate) fn start_on_new_table(&mut self) {
        self.unflagged = Some(-self.window);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107, 162-173
    /// The exit of the last chunk's word finder cleared the table (its offset was at least
    /// INT4_MAX / 4: `Blast_ExtendWordExit`), as it does after any chunk where the offset falls
    /// there, inside a sequence too: the next chunk is at the window (the walk's new start), and
    /// every entry of the table is the cleared one (`last_hit` the window less than 0, 2 windows
    /// below the new start, no flag). The search's walk calls it after such a chunk, so that its
    /// [`ChunkWork::Words`] `reference` stays exact; LOSAT's refusal before the search walks the
    /// chunks without it (a clear only removes entries and moves the later chunks together).
    pub(crate) fn clear(&mut self) {
        self.at = 0;
        self.unflagged = Some(-2 * self.window);
        self.flagged = None;
    }

    /// The margins of the chunks walked so far (after the last chunk, the sequence's): how far
    /// the lowest diagonal offset of a chunk that reads a word that can hit lies below the
    /// sequence's start, and how far its highest entry lies past what every hit of the next
    /// sequence reads as cleared (the offset after the last chunk, less the window for an entry
    /// without the flag). Where either is above 0, NCBI's rows of the sequence depend on the
    /// sequences that its thread searches before and after it, and LOSAT refuses it before the
    /// search (`HardMaskShape::StaleDiagonals`). A search that started the sequence's offset higher
    /// by the first margin and the next sequence's higher by the second would isolate it instead
    /// (the alternative of question Q5 of M/m4b/OWNER_QUESTIONS.md, not taken).
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
}

#[cfg(test)]
#[path = "hard_chunk_walk_tests.rs"]
mod tests;
