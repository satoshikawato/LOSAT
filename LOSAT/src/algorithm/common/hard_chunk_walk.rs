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
//! finder reads of each chunk ([`ChunkWork`]), the margins by which the sequence's word hits
//! reach past it ([`DiagonalMargins`]), and from them how the search moves the diagonal offset
//! so that no entry of another sequence reaches a sequence's hits ([`DiagonalWalk::start`],
//! [`DiagonalWalk::after`]: PD-LOSAT-NCBI-DEFECTS exception 29, the Owner's answer Q5 = B of M4b,
//! 2026-11-22). [`NcbiWalk`] follows NCBI's one-thread offsets over the searched sequences before
//! the search and names the sequences whose hits meet entries that no sequence wrote, which the
//! search starts on such a table at NCBI's offset ([`HardStarts`], [`TableStart::New`]).

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
    /// (`common::hard_chunk_scan::OneWord`).
    OneWord,
    /// A chunk of a one-hit search shorter than the word: its scan range is empty, nothing is read.
    Nothing,
}

/// How far a sequence's word hits reach beyond the sequence on the diagonal table
/// ([`DiagonalWalk::margins`]), in positions of the diagonal offset. Where both are 0 the table
/// entries that the sequence's chunks store read as cleared entries to every hit of a later
/// sequence, and its own hits read the entries of earlier sequences as cleared ones: the sequence's
/// result does not depend on the entries of the sequences that the same table searched before or
/// after it. The search moves the offset by them ([`DiagonalWalk::start`], [`DiagonalWalk::after`]).
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

/// How the search of a sequence starts on the diagonal table ([`DiagonalWalk::start`]; the
/// search applies it to its table with [`TableStart::apply`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TableStart {
    /// A table whose every entry is the one that no sequence wrote since NCBI's thread made or
    /// cleared its table, at NCBI's one-thread offset of the sequence ([`NcbiStart`],
    /// [`HardTable::ncbi_table`]): a sequence whose hits meet such entries in NCBI's one-thread
    /// search, or inside which that search clears the table ([`NcbiWalk`]).
    New(NcbiStart),
    /// The table and its offset as they are (NCBI's offsets until the first move).
    Keep,
    /// The offset moved up by `by` before the sequence's first chunk, after the table is cleared
    /// (its offset the window) where `clear`.
    Move { clear: bool, by: i32 },
}

/// How the search moves the diagonal offset after a sequence's last chunk
/// ([`DiagonalWalk::after`]; the search applies it with [`TableEnd::apply`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TableEnd {
    /// The offset as the last chunk's exit left it.
    Keep,
    /// The offset moved up by this.
    Move(i32),
    /// The table cleared, where the move would take the offset to INT4_MAX / 4 or past it.
    Clear,
}

/// The diagonal table of a protein program's word finder as a `-db_hard_mask` search moves it
/// between its sequences (BLASTP's `blastp::blast_engine::BlastpSubjectScratch`, BLASTX's
/// `blastx::seed::Diagonals`), so that both programs' searches change their tables the one way of
/// [`TableStart::apply`] and [`TableEnd::apply`].
pub(crate) trait HardTable {
    /// The error of a move that takes the offset past Int4 (none in BLASTP, whose walk keeps the
    /// offset below INT4_MAX / 4).
    type Error;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 87-107, 131-144
    /// Every entry the one that NCBI's thread made (0 without the flag: the `calloc` of
    /// `BlastExtendWordNew`) or cleared (`s_BlastDiagClear`: the window `window` less than 0),
    /// and the offset NCBI's one-thread offset of the sequence (`start`): only this worker's table.
    fn ncbi_table(&mut self, start: NcbiStart, window: i32);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
    /// The table cleared as `s_BlastDiagClear` clears it: every entry the window `window` less
    /// than 0, the offset the window.
    fn clear_table(&mut self, window: i32);
    // No NCBI counterpart: PD-LOSAT-NCBI-DEFECTS exception 29 (the Owner's answer Q5 = B of M4b).
    /// The offset moved up by `by` between two sequences, so that no entry of one reaches the
    /// other's hits ([`DiagonalWalk::start`], [`DiagonalWalk::after`]).
    fn move_offset(&mut self, by: i32) -> Result<(), Self::Error>;
}

impl TableStart {
    /// The start applied to the search's table `table` (the search's window `window`): NCBI's
    /// table at the sequence's one-thread offset, nothing, or a clear where `clear` and the move.
    pub(crate) fn apply<T: HardTable>(self, table: &mut T, window: i32) -> Result<(), T::Error> {
        match self {
            TableStart::New(start) => table.ncbi_table(start, window),
            TableStart::Keep => {}
            TableStart::Move { clear, by } => {
                if clear {
                    table.clear_table(window);
                }
                table.move_offset(by)?;
            }
        }
        Ok(())
    }
}

impl TableEnd {
    /// The end applied to the search's table `table` (the search's window `window`).
    pub(crate) fn apply<T: HardTable>(self, table: &mut T, window: i32) -> Result<(), T::Error> {
        match self {
            TableEnd::Keep => {}
            TableEnd::Move(by) => table.move_offset(by)?,
            TableEnd::Clear => table.clear_table(window),
        }
        Ok(())
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
    /// Whether a chunk whose scan reads a word that can hit lies within reach of an entry that an
    /// earlier chunk can store (one without the flag less than the window behind it, one with the
    /// flag ahead of it): where a clear of the table falls between the two, it changes the
    /// sequence's rows.
    reads_back: bool,
    /// The diagonal offset of the last chunk walked, and the highest of the chunks before it
    /// (after whose exit a clear falls inside the sequence).
    latest: Option<i64>,
    inner: Option<i64>,
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
            reads_back: false,
            latest: None,
            inner: None,
        }
    }

    /// The next chunk, at `offset` in the sequence with the length `length` (NCBI's `Int4`): what
    /// its word finder reads. The diagonal offset then moves on by the length and the window.
    pub(crate) fn chunk(&mut self, offset: usize, length: i32) -> ChunkWork {
        let (d, l, w, window) = (self.at, i64::from(length), self.word, self.window);
        let start = i64::try_from(offset).unwrap_or(i64::MAX);
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
        // A one word that covers the sentinel byte (at `length`) without a lookup word that holds
        // the residue 0 has no hit.
        let hits = match work {
            ChunkWork::Nothing => false,
            ChunkWork::OneWord => self.zero_words || start.saturating_add(w) <= self.length,
            ChunkWork::Words { .. } => true,
        };
        if hits {
            if self.unflagged.is_some_and(|stored| stored > d - window)
                || self.flagged.is_some_and(|stored| stored > d)
            {
                self.reads_back = true;
            }
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 87-107, 110-155
    /// The sequence starts on the table of [`TableStart::New`]: every entry is the one that no
    /// sequence wrote since NCBI's thread made its table (0 without the flag: the `calloc` of
    /// `BlastExtendWordNew`) or cleared it (the window less than 0, `s_BlastDiagClear`), and the
    /// offset is NCBI's one-thread offset of the sequence (`start`), so the entries lie that offset
    /// less their value below the sequence's start, where a chunk that goes back below the start
    /// can meet them (a hit within the window of them is extended, one in the one-hit search below
    /// them is skipped). The search of the sequence makes its table so first, and so does the walk
    /// of the refusal before the search (`HardMaskShape::of_scan`).
    pub(crate) fn start_on_ncbi_table(&mut self, start: NcbiStart) {
        self.unflagged = Some(i64::from(start.unwritten(self.window())) - i64::from(start.offset));
    }

    /// The walk of the search of a sequence of `length` letters without masks of the algorithm on
    /// the table of [`TableStart::New`] at `start`, as one chunk of its letters (the search's exit
    /// moves the offset by them and the window; NCBI's chunks of a longer one move it further,
    /// which only lowers the entries below what the walk takes): its end
    /// ([`DiagonalWalk::after`]) moves the offset past the entries that no sequence wrote as well
    /// as past its own, so that the next sequence that the table searches (not always the next one
    /// in the database: a worker's next subject) meets none of them.
    pub(crate) fn plain(word: i32, window: i32, length: usize, start: NcbiStart) -> Self {
        let mut walk = Self::new(word, window, length, true);
        walk.start_on_ncbi_table(start);
        walk.chunk(0, i32::try_from(length).unwrap_or(i32::MAX));
        walk
    }

    /// The two-hit window (0: the one-hit search).
    fn window(&self) -> i32 {
        i32::try_from(self.window).unwrap_or(i32::MAX)
    }

    /// The walk of the chunks `chunks` of a sequence of `length` letters from its start, as the
    /// search's walk reads them (a one word that covers the sentinel taken to hit).
    fn ahead(
        word: i32,
        window: i32,
        length: usize,
        chunks: &[crate::core::subject_split::HardChunk],
    ) -> Self {
        let mut ahead = Self::new(word, window, length, true);
        for chunk in chunks {
            ahead.chunk(chunk.offset, chunk.length);
        }
        ahead
    }

    /// The sequence does not start on a new table, and the search moved the offset up by `before`
    /// ([`DiagonalWalk::start`]): every entry of the earlier sequences lies at least
    /// the window below the offset before the move without the flag, and not above it with the
    /// flag (each sequence's search leaves the table so, [`DiagonalWalk::after`]), so at most the
    /// window plus `before`, or `before`, below the sequence's start. No chunk that reads a word
    /// that can hit lies below `-before`, so none of its hits is within reach of them; the walk
    /// keeps them for its margin past the end, which they set where the chunks end below
    /// `-before` (a walk that goes back more than it goes on).
    pub(crate) fn start_after_others(&mut self, before: i64) {
        self.unflagged = Some(-self.window - before);
        self.flagged = Some(-before);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551, 771-802
    // A hit reads the entry of its diagonal: a two-hit hit at least the window past an entry
    // without the flag, or not behind one with it, stores itself as on a cleared entry, and so
    // does a one-hit hit not behind its entry; a hit within reach of an entry is extended with it
    // or skipped. The exit of a chunk clears the table where the offset has reached INT4_MAX / 4.
    /// How the search of a sequence of `length` letters whose hard chunks are `chunks` (in walk
    /// order) starts on a table at the diagonal offset `offset`, in a search with the word `word`
    /// and the window `window`, and the walk of its search (the margin of [`DiagonalWalk::after`]
    /// comes from it). PD-LOSAT-NCBI-DEFECTS exception 29 (the Owner's answer Q5 = B of M4b,
    /// 2026-11-22 22:34 JST): the rows differ from NCBI's one-thread rows only where an entry of
    /// another sequence reaches a hit, and are the same at every thread count.
    /// - `ncbi`, a sequence whose hits meet the entries that no sequence wrote in NCBI's one-thread
    ///   search, or inside which that search clears the table ([`NcbiWalk`], before the search):
    ///   the table of [`TableStart::New`] at NCBI's one-thread offset, so that its walk meets those
    ///   entries and those clears as NCBI's does ([`DiagonalWalk::start_on_ncbi_table`]).
    /// - Any other: the offset moves up by the margin below its start (the walk of its chunks,
    ///   [`DiagonalMargins::below_start`]), so that every entry of the earlier sequences lies out of
    ///   reach of its hits ([`DiagonalWalk::start_after_others`]): NCBI's one-thread search of it
    ///   where no entry of another sequence reaches it, whose entries that no sequence wrote it
    ///   does not meet. Where an entry of another sequence does reach it, NCBI's rows depend on
    ///   which sequences its thread searched before, which varies by run with threads (probe
    ///   `lk400` of M/m4b3b/NOTES_M4B3B.md).
    /// - It is also kept apart from the clears of its own walk: where a chunk before its last
    ///   would lie at INT4_MAX / 4 or past it after the move, the table is cleared first (a clear
    ///   after that chunk would remove the entries that the chunks after it read, and a one-hit
    ///   chunk below the chunk after the clear would meet the cleared entries; NCBI's one-thread
    ///   search has no clear inside it, else it is `ncbi`). From the cleared offset its chunks
    ///   lie within the walk's extent of the window, which reaches INT4_MAX / 4 (536,870,911
    ///   positions) only for a protein whose chunks move the offset that far: each chunk moves it
    ///   by its letters plus the window, and the chunks with words are disjoint parts of the
    ///   sequence (else `HardMaskShape::PartsOutOfOrder`), so it takes more than 536 million
    ///   letters or about 12 million masks (at the window 40 and the word 3).
    /// - A sequence with no margin below its start whose walk does not read back and has no clear
    ///   inside keeps its offset ([`TableStart::Keep`]).
    pub(crate) fn start(
        word: i32,
        window: i32,
        length: usize,
        chunks: &[crate::core::subject_split::HardChunk],
        ncbi: Option<NcbiStart>,
        offset: i32,
    ) -> (Self, TableStart) {
        let mut walk = Self::new(word, window, length, true);
        if let Some(start) = ncbi {
            walk.start_on_ncbi_table(start);
            return (walk, TableStart::New(start));
        }
        let ahead = Self::ahead(word, window, length, chunks);
        let before = ahead.margins().below_start;
        walk.start_after_others(before);
        let clear = ahead
            .inner
            .is_some_and(|inner| i64::from(offset) + before + inner >= CLEAR_AT);
        if before == 0 && !ahead.reads_back && !clear {
            return (walk, TableStart::Keep);
        }
        // Below INT4_MAX / 4 as long as the extent is (above).
        let by = i32::try_from(before).unwrap_or(i32::MAX);
        (walk, TableStart::Move { clear, by })
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:216, 257
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-573, 771-782
    // NCBI searches each query batch on a table of its own. Where every entry within reach of a
    // hit was written by a hit of the same chunk's scan, a search of several batches on one table
    // gives each batch's result (on each diagonal a batch's hits follow those of the batches
    // before it, whose entries read as out of reach); an entry that no sequence wrote, or one of
    // an earlier chunk after a walk back, reads as the batch's own hits on its diagonal and the
    // context test of the batch's offsets decide.
    /// Whether a search that holds several NCBI query batches on one table (BLASTP's) keeps the
    /// sequence of `length` letters whose hard chunks are `chunks` (none: a sequence without masks)
    /// apart, on NCBI's tables of each batch (`common::query_batch_tables::QueryBatchTables`), and
    /// at which start of those tables: a sequence that starts on a new table (`ncbi`, from
    /// [`NcbiWalk`]) at NCBI's one-thread offset, and one whose walk reads back (a chunk within
    /// reach of an entry that an earlier chunk can store) at the window plus its margin below its
    /// start, on entries of 0, which no hit reaches (every entry of another sequence out of reach,
    /// as [`DiagonalWalk::start_after_others`] takes them; no clear falls inside it there, nor in
    /// NCBI's one-thread search, else it is `ncbi`). `None`: every hit of the sequence reaches only
    /// entries of its own chunk's scan, or entries out of reach, so it shares the search's table.
    pub(crate) fn apart(
        word: i32,
        window: i32,
        length: usize,
        chunks: &[crate::core::subject_split::HardChunk],
        ncbi: Option<NcbiStart>,
    ) -> Option<NcbiStart> {
        if ncbi.is_some() {
            return ncbi;
        }
        let ahead = Self::ahead(word, window, length, chunks);
        ahead.reads_back.then(|| NcbiStart {
            offset: i32::try_from(i64::from(window) + ahead.margins().below_start)
                .unwrap_or(i32::MAX),
            cleared: false,
        })
    }

    /// How the search moves the offset after the sequence's last chunk (`self`, the walk of its
    /// search from [`DiagonalWalk::start`]), from the offset `offset` that the last chunk's exit
    /// left: up by the margin past its end, so that no entry of it or of the sequences before it
    /// reaches the next sequence's hits (every entry then lies at least the window below the
    /// offset without the flag, and not above it with the flag, as
    /// [`DiagonalWalk::start_after_others`] takes them for the next sequence), or a clear of the
    /// table where that move would take the offset to INT4_MAX / 4 or past it (a cleared table
    /// holds no entry within reach either). Where the margin is 0 the offset stays as the last
    /// chunk's exit left it.
    pub(crate) fn after(&self, offset: i32) -> TableEnd {
        let past_end = self.margins().past_end;
        if past_end == 0 {
            TableEnd::Keep
        } else if i64::from(offset) + past_end >= CLEAR_AT {
            TableEnd::Clear
        } else {
            TableEnd::Move(i32::try_from(past_end).unwrap_or(i32::MAX))
        }
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
    /// without the flag). Where either is above 0, NCBI's rows of the sequence, or of the next
    /// one, depend on the sequences that its thread searches before and after it; LOSAT's search
    /// moves the offset by them ([`DiagonalWalk::start`], [`DiagonalWalk::after`]).
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
/// The diagonal offset at which the exit of a chunk's word finder clears the table instead of
/// moving the offset on (INT4_MAX / 4, `Blast_ExtendWordExit`).
const CLEAR_AT: i64 = (i32::MAX / 4) as i64;

/// Where NCBI's one-thread search starts a sequence on its thread's diagonal table
/// ([`NcbiWalk`]): the offset before the sequence's first chunk, and whether the table was cleared
/// before it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NcbiStart {
    /// The diagonal offset before the sequence's first chunk.
    pub(crate) offset: i32,
    /// Whether the exit of an earlier chunk cleared the table (at INT4_MAX / 4).
    pub(crate) cleared: bool,
}

impl NcbiStart {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107, 131-144
    /// The value of an entry (without the flag) that no sequence wrote since NCBI's thread made
    /// its table (0: the `calloc` of `BlastExtendWordNew`) or cleared it (the window less than 0:
    /// `s_BlastDiagClear`), in a search with the window `window` (0: the one-hit search).
    pub(crate) fn unwritten(self, window: i32) -> i32 {
        if self.cleared {
            -window
        } else {
            0
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1432
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-276, 476-500
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66, 87-107, 162-173
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:230-270
// Each search (a batch of queries, or a chunk of a split query) makes its tables: a thread's
// starts with every entry 0 and the offset the window. With one thread the search reads the
// database's sequences in the order of their ids (those that the lists and the alias files'
// ranges include), skips one of no letters (the minimum length of a protein subject is 1), and
// passes every chunk of the others (one per hard range, or per `MAX_DBSEQ_LEN` letters of a
// longer range) through the word finder, whose exit clears the table where the offset has reached
// INT4_MAX / 4 (every entry the window less than 0, the offset the window), else moves the offset
// by the chunk's length (NCBI's `Int4`, negative for an inverted range) and the window. A hit at
// subject offset `s` of a chunk at offset `D` on an entry `u` without the flag stores itself as on
// a cleared entry only where `s + D - u` is at least the window (two-hit), or not below 0
// (one-hit, whose window is 0).
/// NCBI's one-thread walk of the diagonal offset over the searched sequences of a `-db_hard_mask`
/// search with the word `word` and the window `window`, in the order of the search, before the
/// search: for each sequence its [`NcbiStart`], and whether the search starts it on a table of
/// entries that no sequence wrote at that offset ([`TableStart::New`]): where a hit of its walk
/// can meet such an entry (its lowest chunk that reads a word lies less than the window past the
/// entries, or below them in the one-hit search: the first searched sequence's walk below its
/// start, a sequence searched within a few hundred letters of the start or of a clear, or one
/// after walks that went back), or where NCBI's walk clears the table after one of its chunks
/// before the last. Every other sequence meets such entries nowhere, and no clear falls inside
/// it, so that its result is NCBI's wherever no entry of another sequence reaches it, at any
/// offset ([`DiagonalWalk::start`]).
#[derive(Debug, Clone)]
pub(crate) struct NcbiWalk {
    word: i32,
    window: i32,
    offset: i64,
    cleared: bool,
}

impl NcbiWalk {
    /// The walk at the start of a search: the thread's new table, its offset the window.
    pub(crate) fn new(word: i32, window: i32) -> Self {
        Self {
            word,
            window,
            offset: i64::from(window),
            cleared: false,
        }
    }

    /// The exit of a chunk of `length` letters (NCBI's `Int4`).
    fn exit(&mut self, length: i32) {
        if self.offset >= CLEAR_AT {
            self.offset = i64::from(self.window);
            self.cleared = true;
        } else {
            self.offset += i64::from(length) + i64::from(self.window);
        }
    }

    /// The start of the next sequence.
    fn here(&self) -> NcbiStart {
        NcbiStart {
            offset: self.offset.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32,
            cleared: self.cleared,
        }
    }

    /// Whether a hit of a chunk `below` positions below the start `start` can meet an entry that
    /// no sequence wrote (a hit at the chunk's first letter is the closest).
    fn meets_unwritten(&self, start: NcbiStart, below: i64) -> bool {
        i64::from(start.offset) - below - i64::from(start.unwritten(self.window))
            < i64::from(self.window)
    }

    /// A searched sequence of `length` letters without hard ranges of its own (no masks of the
    /// algorithm): its one chunk, or the chunks of `MAX_DBSEQ_LEN` letters of a longer one, each
    /// at or past the one before it. `Some` where the search starts it on a new table
    /// ([`TableStart::New`]); a sequence of no letters is not searched.
    pub(crate) fn plain(&mut self, length: usize) -> Option<NcbiStart> {
        use crate::core::subject_split::{HardSplitState, DBSEQ_CHUNK_OVERLAP, MAX_DBSEQ_LEN};
        if length == 0 {
            return None;
        }
        let start = self.here();
        let new = self.meets_unwritten(start, 0);
        if length <= MAX_DBSEQ_LEN {
            self.exit(length as i32);
        } else {
            let mut split = HardSplitState::default();
            split.start(length, vec![(0, length as i32)], false);
            while let Some(chunk) = split.next_chunk(DBSEQ_CHUNK_OVERLAP) {
                self.exit(chunk.length);
            }
        }
        new.then_some(start)
    }

    /// A searched sequence of `length` letters whose hard chunks are `chunks` (in walk order):
    /// `Some` where the search starts it on a new table ([`TableStart::New`]: its walk, as the
    /// search's reads it, reaches an entry that no sequence wrote, or a chunk before its last
    /// exits at INT4_MAX / 4 or past it), and whether its walk reads back (a chunk within reach of
    /// an entry that an earlier chunk can store: [`DiagonalWalk::apart`] without a new table).
    pub(crate) fn masked(
        &mut self,
        length: usize,
        chunks: &[crate::core::subject_split::HardChunk],
    ) -> (Option<NcbiStart>, bool) {
        if length == 0 {
            return (None, false);
        }
        let start = self.here();
        let ahead = DiagonalWalk::ahead(self.word, self.window, length, chunks);
        let clears_inside = ahead
            .inner
            .is_some_and(|inner| self.offset + inner >= CLEAR_AT);
        let new = self.meets_unwritten(start, ahead.margins().below_start) || clears_inside;
        for chunk in chunks {
            self.exit(chunk.length);
        }
        (new.then_some(start), ahead.reads_back)
    }
}

/// The searched sequences of a `-db_hard_mask` search that start on a new table
/// ([`TableStart::New`]), by ordinal id, with their [`NcbiStart`]: what [`NcbiWalk`] gives, before
/// the search (`blastinput::db_mask::ResolvedMask::check_losat_limits`). Few: sequences near the
/// start of the database or of a clear, or after walks that went back. With them, whether the walk
/// keeps any searched sequence apart ([`DiagonalWalk::apart`]: one of these, or one whose walk
/// reads back).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct HardStarts {
    starts: Vec<(u32, NcbiStart)>,
    apart: bool,
}

impl HardStarts {
    /// Sequence `ordinal`, after every sequence already added, starts at `start`.
    pub(crate) fn push(&mut self, ordinal: u32, start: NcbiStart) {
        debug_assert!(self.starts.last().is_none_or(|&(last, _)| last < ordinal));
        self.starts.push((ordinal, start));
        self.apart = true;
    }

    /// A searched sequence whose walk reads back ([`DiagonalWalk::apart`] without a new table).
    pub(crate) fn push_reads_back(&mut self) {
        self.apart = true;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
    // A hit that meets an entry no chunk of its own scan wrote reads it through the ring of its
    // query batch's length and through the unsigned context test of the batch's query offsets, so
    // its outcome depends on where its query lies in the batch (cited at `DiagonalWalk::apart`).
    /// Whether the walk keeps a searched sequence apart ([`DiagonalWalk::apart`]): a sequence whose
    /// hits can meet entries that its chunks' own scans did not write, so that its rows depend on
    /// its query's place in NCBI's query batch (BLASTX then searches the batch as one unit,
    /// `blastx::query_setup::unit_ranges`).
    pub(crate) fn keeps_apart(&self) -> bool {
        self.apart
    }

    /// Where the search starts sequence `ordinal` on a new table (`None`: it does not).
    pub(crate) fn of(&self, ordinal: u32) -> Option<NcbiStart> {
        self.starts
            .binary_search_by_key(&ordinal, |&(at, _)| at)
            .ok()
            .map(|index| self.starts[index].1)
    }
}

#[cfg(test)]
#[path = "hard_chunk_walk_tests.rs"]
mod tests;
