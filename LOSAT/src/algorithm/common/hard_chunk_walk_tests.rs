use super::*;

/// The work, the margins, the end and the first word past the end of a walk over `chunks`
/// (`(offset, length)`) of a sequence of `length` letters, for a query whose lookup table can
/// (`zero_words`) or cannot hold a word with the residue 0.
fn walk_of(
    chunks: &[(usize, i32)],
    length: usize,
    word: i32,
    window: i32,
    zero_words: bool,
) -> (Vec<ChunkWork>, DiagonalMargins, i64, Option<usize>) {
    let mut walk = DiagonalWalk::new(word, window, length, zero_words);
    let works = chunks
        .iter()
        .map(|&(offset, length)| walk.chunk(offset, length))
        .collect();
    (works, walk.margins(), walk.end(), walk.word_past_end())
}

/// [`walk_of`] for a query without the residue 0.
fn walk(
    chunks: &[(usize, i32)],
    length: usize,
    word: i32,
    window: i32,
) -> (Vec<ChunkWork>, DiagonalMargins, i64, Option<usize>) {
    walk_of(chunks, length, word, window, false)
}

const WORDS: ChunkWork = ChunkWork::Words { reference: false };
const WORDS_BACK: ChunkWork = ChunkWork::Words { reference: true };

/// The chunks of probe 2 of M/m4b3b/NOTES_M4B3B.md (a sequence of 240 letters, word 3, window 40
/// or 0): the offsets `D_i` and what each chunk's scan reads, and the margins of the sequence.
#[test]
fn walks_of_overlapping_and_nested_masks() {
    // e3: masks [40, 60) [55, 70). Two-hit: D = 0, 80, 115; the last chunk starts past every
    // entry of the earlier ones.
    let e3 = [(0, 40), (60, -5), (70, 170)];
    let (works, margins, end, past) = walk(&e3, 240, 3, 40);
    assert_eq!(works, [WORDS, ChunkWork::OneWord, WORDS]);
    assert_eq!(margins, DiagonalMargins::default());
    assert_eq!((end, past), (325, None));
    // One-hit: D = 0, 40, 35; the short chunk reads nothing, the last goes back below the first
    // chunk's entries (39 at most) but stays above the start.
    let (works, margins, end, _) = walk(&e3, 240, 3, 0);
    assert_eq!(works, [WORDS, ChunkWork::Nothing, WORDS_BACK]);
    assert_eq!(margins, DiagonalMargins::default());
    assert_eq!(end, 205);
    // nest: masks [40, 100) [50, 60). Two-hit: D = 0, 80, 70, the last chunk below the one word's
    // entries (81 with the flag).
    let nest = [(0, 40), (100, -50), (60, 180)];
    let (works, margins, end, _) = walk(&nest, 240, 3, 40);
    assert_eq!(works, [WORDS, ChunkWork::OneWord, WORDS_BACK]);
    assert_eq!(margins, DiagonalMargins::default());
    assert_eq!(end, 290);
    // One-hit: D = 0, 40, -10: the last chunk's hits reach 10 positions below the start.
    let (_, margins, _, _) = walk(&nest, 240, 3, 0);
    assert_eq!(
        margins,
        DiagonalMargins {
            below_start: 10,
            past_end: 0
        }
    );
    // big: masks [40, 120) [60, 130), an overlap longer than the window. Two-hit: D = 0, 80, 60.
    let big = [(0, 40), (120, -60), (130, 110)];
    let (works, margins, _, _) = walk(&big, 240, 3, 40);
    assert_eq!(works, [WORDS, ChunkWork::OneWord, WORDS_BACK]);
    assert_eq!(margins, DiagonalMargins::default());
    // lk: masks [40, 200) [60, 238). Two-hit: D = 0, 80, -20 and the end 22: the one word at 80
    // leaves an entry 98 past the end less the window. The one word at 238 covers the sentinel:
    // without a query word that holds the residue 0 it has no hit, else it reaches 20 below the
    // start.
    let lk = [(0, 40), (200, -140), (238, 2)];
    let (works, margins, end, _) = walk(&lk, 240, 3, 40);
    assert_eq!(works, [WORDS, ChunkWork::OneWord, ChunkWork::OneWord]);
    assert_eq!(
        margins,
        DiagonalMargins {
            below_start: 0,
            past_end: 98
        }
    );
    assert_eq!(end, 22);
    let (_, margins, _, _) = walk_of(&lk, 240, 3, 40, true);
    assert_eq!(
        margins,
        DiagonalMargins {
            below_start: 20,
            past_end: 98
        }
    );
    // One-hit: D = 0, 40, -100 and the end -98: the first chunk's entries (39) lie 137 past it.
    let (works, margins, end, _) = walk(&lk, 240, 3, 0);
    assert_eq!(works, [WORDS, ChunkWork::Nothing, ChunkWork::Nothing]);
    assert_eq!(
        margins,
        DiagonalMargins {
            below_start: 0,
            past_end: 137
        }
    );
    assert_eq!(end, -98);
}

/// The one-hit bound counts a word whose letters score 0 at most: its extension's last offset is
/// one past the chunk (`s_BlastAaExtendOneHit`), so the entry is the chunk's length less the word
/// plus two past its offset. A chunk that ends one before the next chunk's offset would leave
/// that entry on it.
#[test]
fn one_hit_bound_covers_a_zero_score_word() {
    // Chunks of 10 letters with masks of 1 letter between them, so D = 0, 10 (one-hit).
    let (works, margins, _, _) = walk(&[(0, 10), (11, 10)], 21, 3, 0);
    // The first chunk's entries reach 0 + 10 - 3 + 2 = 9, below the second chunk's offset.
    assert_eq!(works, [WORDS, WORDS]);
    assert_eq!(margins, DiagonalMargins::default());
    // A chunk that goes back by 1 (masks [10, 12) [11, 13)): D = 0, 10, 9; 9 is not above 9.
    let (works, _, _, _) = walk(&[(0, 10), (12, -1), (13, 10)], 23, 3, 0);
    assert_eq!(works, [WORDS, ChunkWork::Nothing, WORDS]);
    // Back by 2: D = 0, 10, 8, below the first chunk's highest entry 9.
    let (works, _, _, _) = walk(&[(0, 10), (12, -2), (14, 10)], 24, 3, 0);
    assert_eq!(works, [WORDS, ChunkWork::Nothing, WORDS_BACK]);
}

/// Masks apart in the order of their start: every offset is at least the window past the
/// entries before it, so no chunk goes back, and the margins are 0, at every word (2 to 7, the
/// protein programs' words: `blast_args.cpp`) and window. A one-hit word of 1 letter that scores
/// 0 would leave an entry one past the next chunk's offset.
#[test]
fn walks_of_masks_apart_stay_within_the_sequence() {
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % below
    };
    for _ in 0..20_000 {
        let length = 1 + next(300) as usize;
        let mut chunks = Vec::new();
        let mut left = 0usize;
        while left < length {
            let right = (left + next(40) as usize).min(length);
            chunks.push((left, (right - left) as i32));
            left = right + 1 + next(30) as usize;
        }
        let word = 2 + next(6) as i32;
        let window = [0, 1, 2, 40, 1000][next(5) as usize];
        let zero_words = next(2) == 1;
        let (works, margins, _, past) = walk_of(&chunks, length, word, window, zero_words);
        assert_eq!(
            margins,
            DiagonalMargins::default(),
            "{chunks:?} {word} {window}"
        );
        assert_eq!(past, None, "{chunks:?} {word} {window}");
        assert!(
            works
                .iter()
                .all(|&work| work != ChunkWork::Words { reference: true }),
            "{chunks:?} {word} {window}"
        );
    }
}

/// The one word of a chunk at the sequence's last letter runs past its sentinel: with a query
/// that can hit a word with the residue 0, refused where an earlier entry within the window
/// behind it would extend its hits (masks [10, 60) [20, 99) [30, 40) of a sequence of 100
/// letters), not where none does; without such a query never (the word has no hit).
#[test]
fn words_past_the_end_that_can_change_a_result() {
    let chunks = [(0, 10), (60, -40), (99, -69), (40, 60)];
    let (works, margins, _, past) = walk_of(&chunks, 100, 3, 40, true);
    assert_eq!(
        works,
        [WORDS, ChunkWork::OneWord, ChunkWork::OneWord, WORDS_BACK]
    );
    assert_eq!(margins, DiagonalMargins::default());
    assert_eq!(past, Some(99));
    let (_, margins, _, past) = walk(&chunks, 100, 3, 40);
    assert_eq!(margins, DiagonalMargins::default());
    assert_eq!(past, None);
    // The same word with no entry within the window behind it and no chunk after it.
    let (_, _, _, past) = walk_of(&[(0, 10), (99, -69)], 100, 3, 40, true);
    assert_eq!(past, None);
    // A word that ends at the sentinel reads no other sequence.
    let (_, _, _, past) = walk_of(&[(0, 10), (60, -40), (98, -68), (40, 60)], 100, 3, 40, true);
    assert_eq!(past, None);
}

/// After a clear inside the sequence the next chunk is at the walk's new start and the table
/// holds only cleared entries (2 windows below it): a later chunk that goes back below them is
/// read by the reference scan, one that stays above them is not, whatever the chunks before the
/// clear stored.
#[test]
fn a_clear_inside_the_sequence_keeps_only_cleared_entries() {
    // The one word at 299 covers the sentinel and has no hit: it stores nothing.
    for (back, work) in [(-100, WORDS), (-150, WORDS_BACK)] {
        let mut walk = DiagonalWalk::new(3, 40, 300, false);
        assert_eq!(walk.chunk(0, 100), WORDS);
        walk.clear();
        assert_eq!(walk.chunk(299, back), ChunkWork::OneWord);
        // D = back + 40: -60 above the cleared entries at -80, or -110 below them (the first
        // chunk's entries, up to 98, are gone).
        assert_eq!(walk.chunk(199, 100), work, "{back}");
    }
    // Without the clear the first chunk's entries (98) are above D = 140 - 100 + 40 = 80.
    let mut walk = DiagonalWalk::new(3, 40, 300, false);
    walk.chunk(0, 100);
    walk.chunk(299, -100);
    assert_eq!(walk.chunk(199, 100), WORDS_BACK);
}

/// A sequence on a new table at NCBI's one-thread offset (`DiagonalWalk::start_on_ncbi_table`):
/// the entries that no sequence wrote lie that offset less their value below the sequence's start
/// (the window below it at the start of the search, two windows after a clear), and a chunk that
/// goes back below them is read by the reference scan, as a chunk below an earlier chunk's entries
/// is.
#[test]
fn a_new_table_holds_entries_below_the_start() {
    let at = |offset, cleared| NcbiStart { offset, cleared };
    for (ncbi, work) in [
        (None, WORDS),
        // Entries at -40: D = -60 below them.
        (Some(at(40, false)), WORDS_BACK),
        // Entries at -80 after a clear: D = -60 above them.
        (Some(at(40, true)), WORDS),
        // Entries at -5 (NCBI's offset 5): below them.
        (Some(at(5, false)), WORDS_BACK),
        // Entries at -100: above them.
        (Some(at(100, false)), WORDS),
    ] {
        let mut walk = DiagonalWalk::new(3, 40, 300, false);
        if let Some(start) = ncbi {
            walk.start_on_ncbi_table(start);
        }
        // The one word at 299 covers the sentinel and has no hit: it stores nothing.
        assert_eq!(walk.chunk(299, -100), ChunkWork::OneWord);
        // D = -60.
        assert_eq!(walk.chunk(0, 100), work, "{ncbi:?}");
    }
    // One-hit: the entries lie at the start at the start of the search; a chunk at the start
    // reads them as cleared, one below it does not.
    let mut walk = DiagonalWalk::new(3, 0, 300, false);
    walk.start_on_ncbi_table(at(0, false));
    assert_eq!(walk.chunk(0, 10), WORDS);
    let mut walk = DiagonalWalk::new(3, 0, 300, false);
    walk.start_on_ncbi_table(at(0, false));
    assert_eq!(walk.chunk(299, -10), ChunkWork::Nothing);
    assert_eq!(walk.chunk(0, 100), WORDS_BACK);
    // At NCBI's offset -98 (after a walk that went back) they lie 98 above the start.
    let mut walk = DiagonalWalk::new(3, 0, 300, false);
    walk.start_on_ncbi_table(at(-98, false));
    assert_eq!(walk.chunk(0, 230), WORDS_BACK);
}

/// The chunks `(offset, length)` as the search's hard chunks.
fn hard(chunks: &[(usize, i32)]) -> Vec<crate::core::subject_split::HardChunk> {
    chunks
        .iter()
        .map(|&(offset, length)| crate::core::subject_split::HardChunk {
            offset,
            length,
            residual: 0,
            full_range: false,
            range: 0,
            overlap: 0,
            overlap_past_list: false,
        })
        .collect()
}

/// How the search keeps each sequence on its own (PD-LOSAT-NCBI-DEFECTS exception 29, the Owner's
/// answer Q5 = B of M4b): a sequence that NCBI's one-thread walk starts within reach of the entries
/// that no sequence wrote takes NCBI's table at that offset; another moves the offset up by its
/// margin below the start, after a clear where a chunk before its last would reach INT4_MAX / 4; a
/// walk that does not read back keeps its offset; after the sequence the offset moves up by its
/// margin past the end, or the table is cleared where that would reach INT4_MAX / 4.
#[test]
fn each_sequence_is_searched_on_its_own() {
    const LIMIT: i32 = i32::MAX / 4;
    // nest, one-hit: 10 below the start (the last chunk at -10, the one before it at 40).
    let nest = hard(&[(0, 40), (100, -50), (60, 180)]);
    let first = NcbiStart {
        offset: 0,
        cleared: false,
    };
    let (_, start) = DiagonalWalk::start(3, 0, 240, &nest, Some(first), 5000);
    assert_eq!(start, TableStart::New(first));
    let (mut walk, start) = DiagonalWalk::start(3, 0, 240, &nest, None, 5000);
    assert_eq!(
        start,
        TableStart::Move {
            clear: false,
            by: 10
        }
    );
    // The chunk at 40 lies at the offset + 10 + 40 after the move.
    let (_, start) = DiagonalWalk::start(3, 0, 240, &nest, None, LIMIT - 51);
    assert_eq!(
        start,
        TableStart::Move {
            clear: false,
            by: 10
        }
    );
    let (_, start) = DiagonalWalk::start(3, 0, 240, &nest, None, LIMIT - 50);
    assert_eq!(
        start,
        TableStart::Move {
            clear: true,
            by: 10
        }
    );
    // Its search's walk: no margin past the end.
    for chunk in &nest {
        walk.chunk(chunk.offset, chunk.length);
    }
    assert_eq!(walk.after(5000), TableEnd::Keep);
    // e3, two-hit: its last chunk reads the one word's entry (80) within the window (it reads
    // back), with no margin below the start; the one word's chunk lies at 80.
    let e3 = hard(&[(0, 40), (60, -5), (70, 170)]);
    let (_, start) = DiagonalWalk::start(3, 40, 240, &e3, None, 5000);
    assert_eq!(
        start,
        TableStart::Move {
            clear: false,
            by: 0
        }
    );
    let (_, start) = DiagonalWalk::start(3, 40, 240, &e3, None, LIMIT - 81);
    assert_eq!(
        start,
        TableStart::Move {
            clear: false,
            by: 0
        }
    );
    let (_, start) = DiagonalWalk::start(3, 40, 240, &e3, None, LIMIT - 80);
    assert_eq!(start, TableStart::Move { clear: true, by: 0 });
    // Masks apart: NCBI's offsets; near INT4_MAX / 4 (the chunk at 80 would exit there) a clear
    // first, which keeps the walk's chunks together.
    let apart = hard(&[(0, 40), (50, 100), (160, 80)]);
    for (offset, want) in [
        (5000, TableStart::Keep),
        (LIMIT - 81, TableStart::Keep),
        (LIMIT - 80, TableStart::Move { clear: true, by: 0 }),
        (LIMIT - 10, TableStart::Move { clear: true, by: 0 }),
    ] {
        let (mut walk, start) = DiagonalWalk::start(3, 40, 240, &apart, None, offset);
        assert_eq!(start, want, "{offset}");
        for chunk in &apart {
            walk.chunk(chunk.offset, chunk.length);
        }
        assert_eq!(walk.after(5000), TableEnd::Keep);
    }
    // One-hit chunks 10, 2, -1, 10 (the third below the second's start): without the clear first,
    // a clear after the second chunk would put the third below the cleared entries.
    let below = hard(&[(0, 10), (12, 2), (15, -1), (20, 10)]);
    let (_, start) = DiagonalWalk::start(3, 0, 40, &below, None, 5000);
    assert_eq!(start, TableStart::Keep);
    let (_, start) = DiagonalWalk::start(3, 0, 40, &below, None, LIMIT - 12);
    assert_eq!(start, TableStart::Move { clear: true, by: 0 });
    // lk, two-hit, the search's walk (the one word at 238 covers the sentinel and is taken to
    // hit): 20 below the start and 98 past the end.
    let lk = hard(&[(0, 40), (200, -140), (238, 2)]);
    let (mut walk, start) = DiagonalWalk::start(3, 40, 240, &lk, None, 5000);
    assert_eq!(
        start,
        TableStart::Move {
            clear: false,
            by: 20
        }
    );
    for chunk in &lk {
        walk.chunk(chunk.offset, chunk.length);
    }
    assert_eq!(walk.after(5000), TableEnd::Move(98));
    assert_eq!(walk.after(LIMIT - 99), TableEnd::Move(98));
    assert_eq!(walk.after(LIMIT - 98), TableEnd::Clear);
    // A walk that goes back more than it goes on and stores nothing (one-hit, chunks shorter than
    // the word): the earlier sequences' entries would reach the next sequence's hits, 48 past the
    // end.
    let back = hard(&[(0, 2), (100, -50)]);
    let (mut walk, start) = DiagonalWalk::start(3, 0, 300, &back, None, 5000);
    assert_eq!(start, TableStart::Keep);
    for chunk in &back {
        walk.chunk(chunk.offset, chunk.length);
    }
    assert_eq!(walk.after(5000), TableEnd::Move(48));
}

/// NCBI's one-thread walk over the searched sequences (`NcbiWalk`, DECISIONS 55 of series M): the
/// sequences that the search starts on a new table at NCBI's offset, and the offsets after each.
#[test]
fn ncbi_walk_names_the_sequences_that_meet_unwritten_entries() {
    let at = |offset, cleared| Some(NcbiStart { offset, cleared });
    // nest, one-hit (10 below the start), and the probe zt shape of review 2 (masks [10,100)
    // [20,110) of 240 letters: chunks 10, -80, 130, 70 below the start).
    let nest = hard(&[(0, 40), (100, -50), (60, 180)]);
    let zt = hard(&[(0, 10), (100, -80), (110, 130)]);
    // Two-hit walks 10 and 50 below the start (chunks 40, then one word at 80, then a chunk at
    // -10 or -50).
    let shallow = hard(&[(0, 40), (100, -130), (60, 180)]);
    let deep = hard(&[(0, 40), (100, -170), (10, 180)]);
    // The first searched sequence: on the new table where its walk goes below its start, at the
    // window (two-hit) or 0 (one-hit).
    let mut walk = NcbiWalk::new(3, 0);
    assert_eq!(walk.masked(240, &nest).0, at(0, false));
    let mut walk = NcbiWalk::new(3, 40);
    assert_eq!(walk.masked(240, &shallow).0, at(40, false));
    // A first sequence without masks, or whose walk stays at or past its start, is not.
    let mut walk = NcbiWalk::new(3, 0);
    assert_eq!(walk.plain(240), None);
    let apart = hard(&[(0, 40), (50, 100), (160, 80)]);
    assert_eq!(walk.masked(240, &apart).0, None);
    // Sequences of no letters are not searched: the first sequence with letters is the first.
    let mut walk = NcbiWalk::new(3, 0);
    assert_eq!(walk.plain(0), None);
    assert_eq!(walk.masked(0, &nest).0, None);
    assert_eq!(walk.masked(240, &nest).0, at(0, false));
    // zt: a sequence of 5 letters, then the zt shape at 5 (one-hit), whose walk reaches 70 below.
    let mut walk = NcbiWalk::new(3, 0);
    assert_eq!(walk.plain(5), None);
    assert_eq!(walk.masked(240, &zt).0, at(5, false));
    // After it NCBI's offset is 5 + 10 - 80 + 130 = 65: a sequence without masks is not.
    assert_eq!(walk.plain(240), None);
    // Two-hit, after 5 letters (NCBI's offset 40 + 5 + 40 = 85): zt does not go below its start,
    // 10 below does not reach the entries at -85 within the window, 50 below does.
    for (shape, want) in [(&zt, None), (&shallow, None), (&deep, at(85, false))] {
        let mut walk = NcbiWalk::new(3, 40);
        assert_eq!(walk.plain(5), None);
        assert_eq!(walk.masked(240, shape).0, want);
    }
    // After 70 letters (150): not.
    let mut walk = NcbiWalk::new(3, 40);
    assert_eq!(walk.plain(70), None);
    assert_eq!(walk.masked(240, &deep).0, None);
    // ov_lk, one-hit: A (the first, masks [40,200) [60,238)) walks back to -98, so B, without a
    // margin below its start, meets the entries 98 above its start.
    let lk = hard(&[(0, 40), (200, -140), (238, 2)]);
    let mut walk = NcbiWalk::new(3, 0);
    assert_eq!(walk.masked(240, &lk).0, None);
    assert_eq!(
        walk.masked(240, &hard(&[(0, 230), (235, 5)])).0,
        at(-98, false)
    );
    // A sequence without masks there too.
    let mut walk = NcbiWalk::new(3, 0);
    walk.masked(240, &lk);
    assert_eq!(walk.plain(240), at(-98, false));
    // Two-hit ov_lk: A, the first, walks 20 below its start (the one word at 238 taken to hit);
    // it ends at 40 + 40 + 40 - 140 + 40 + 2 + 40 = 62, and B is not.
    let mut walk = NcbiWalk::new(3, 40);
    assert_eq!(walk.masked(240, &lk).0, at(40, false));
    assert_eq!(walk.plain(240), None);
    // Whether the walk reads back (`DiagonalWalk::apart` without a new table; line NP item
    // M4b-3b2): two-hit, ov_e3's chunk at 70 lies within the window of the entries of its one
    // word at 60; the touching masks of ov_e3c read nothing back.
    let e3 = hard(&[(0, 40), (60, -5), (70, 170)]);
    let e3c = hard(&[(0, 40), (60, 0), (70, 170)]);
    assert_eq!(NcbiWalk::new(3, 40).masked(240, &e3), (None, true));
    assert_eq!(NcbiWalk::new(3, 40).masked(240, &e3c), (None, false));
    assert_eq!(NcbiWalk::new(3, 0).masked(240, &apart), (None, false));
}

/// NCBI's clears in its one-thread walk (`NcbiWalk`): after the clear the entries lie the window
/// below 0 (two-hit), so a walk meets them only below its start by more than the offset; a clear
/// inside a sequence starts it on a new table at NCBI's offset.
#[test]
fn ncbi_walk_follows_the_clears() {
    const LIMIT: i64 = (i32::MAX / 4) as i64;
    let near = |offset: i64| NcbiWalk {
        word: 3,
        window: 40,
        offset,
        cleared: false,
    };
    let start = |offset, cleared| Some(NcbiStart { offset, cleared });
    // The exit of a sequence at INT4_MAX / 4 clears the table: the next sequence starts at the
    // window on cleared entries.
    let mut walk = near(LIMIT);
    assert_eq!(walk.plain(240), None);
    assert_eq!(
        walk.here(),
        NcbiStart {
            offset: 40,
            cleared: true
        }
    );
    // 50 below the start: 40 - 50 + 40 = 30 < 40 meets them; 10 below does not.
    let shallow = hard(&[(0, 40), (100, -130), (60, 180)]);
    let deep = hard(&[(0, 40), (100, -170), (10, 180)]);
    let mut walk = near(LIMIT);
    walk.plain(240);
    assert_eq!(walk.masked(240, &shallow).0, None);
    let mut walk = near(LIMIT);
    walk.plain(240);
    assert_eq!(walk.masked(240, &deep).0, start(40, true));
    // Plain sequences walk the offset to the clear.
    let mut walk = NcbiWalk::new(3, 40);
    while !walk.cleared {
        assert_eq!(walk.plain(1_000_000), None);
    }
    assert_eq!(walk.here().offset, 40);
    // One longer than `MAX_DBSEQ_LEN` exits once per chunk.
    let mut walk = NcbiWalk::new(3, 40);
    let long = crate::core::subject_split::MAX_DBSEQ_LEN + 1;
    walk.plain(long);
    assert!(walk.offset > (long as i64) + 80, "{}", walk.offset);
    // A clear inside a sequence whose walk does not go below its start: the chunk at 80 exits at
    // INT4_MAX / 4.
    let apart = hard(&[(0, 40), (50, 100), (160, 80)]);
    assert_eq!(near(LIMIT - 81).masked(240, &apart).0, None);
    assert_eq!(
        near(LIMIT - 80).masked(240, &apart).0,
        start((LIMIT - 80) as i32, false)
    );
    // The last chunk's exit is after the sequence: the clear falls before the next one.
    assert_eq!(near(LIMIT - 220).masked(240, &apart).0, None);
}

/// The table of the new-table starts (`HardStarts`): by ordinal id.
#[test]
fn hard_starts_by_ordinal() {
    let mut starts = HardStarts::default();
    let at = |offset| NcbiStart {
        offset,
        cleared: false,
    };
    starts.push(0, at(40));
    starts.push(7, at(-98));
    starts.push(4096, at(5));
    assert_eq!(starts.of(0), Some(at(40)));
    assert_eq!(starts.of(7), Some(at(-98)));
    assert_eq!(starts.of(4096), Some(at(5)));
    for ordinal in [1, 6, 8, 4095, 4097] {
        assert_eq!(starts.of(ordinal), None);
    }
    assert_eq!(HardStarts::default().of(0), None);
}

/// Whether the walk keeps a sequence apart (`HardStarts::keeps_apart`, line NP item M4b-3b2): a
/// sequence on a new table or one whose walk reads back; a walk with neither keeps none apart.
#[test]
fn hard_starts_keep_a_sequence_apart_on_a_new_table_or_after_a_read_back() {
    let mut starts = HardStarts::default();
    assert!(!starts.keeps_apart());
    starts.push_reads_back();
    assert!(starts.keeps_apart());
    assert_eq!(starts.of(0), None);
    let mut starts = HardStarts::default();
    starts.push(
        3,
        NcbiStart {
            offset: 10,
            cleared: false,
        },
    );
    assert!(starts.keeps_apart());
}

/// A table that records what [`TableStart::apply`] and [`TableEnd::apply`] do to it.
#[derive(Debug, Default, PartialEq, Eq)]
struct Recorded {
    ops: Vec<String>,
}

impl HardTable for Recorded {
    type Error = ();
    fn ncbi_table(&mut self, start: NcbiStart, window: i32) {
        self.ops
            .push(format!("new {} {} {window}", start.offset, start.cleared));
    }
    fn clear_table(&mut self, window: i32) {
        self.ops.push(format!("clear {window}"));
    }
    fn move_offset(&mut self, by: i32) -> Result<(), ()> {
        self.ops.push(format!("move {by}"));
        if by == i32::MAX {
            return Err(());
        }
        Ok(())
    }
}

/// The starts and ends as the searches apply them to their tables: NCBI's table, nothing, the
/// clear before the move, the move's error passed on.
#[test]
fn starts_and_ends_apply_to_the_table() {
    let ops = |start: Option<TableStart>, end: Option<TableEnd>| {
        let mut table = Recorded::default();
        if let Some(start) = start {
            start.apply(&mut table, 40).unwrap();
        }
        if let Some(end) = end {
            end.apply(&mut table, 40).unwrap();
        }
        table.ops
    };
    let start = NcbiStart {
        offset: -98,
        cleared: true,
    };
    assert_eq!(ops(Some(TableStart::New(start)), None), ["new -98 true 40"]);
    assert!(ops(Some(TableStart::Keep), Some(TableEnd::Keep)).is_empty());
    assert_eq!(
        ops(Some(TableStart::Move { clear: true, by: 7 }), None),
        ["clear 40", "move 7"]
    );
    assert_eq!(
        ops(
            Some(TableStart::Move {
                clear: false,
                by: 0
            }),
            None
        ),
        ["move 0"]
    );
    assert_eq!(ops(None, Some(TableEnd::Move(98))), ["move 98"]);
    assert_eq!(ops(None, Some(TableEnd::Clear)), ["clear 40"]);
    let mut table = Recorded::default();
    assert_eq!(
        TableStart::Move {
            clear: false,
            by: i32::MAX
        }
        .apply(&mut table, 40),
        Err(())
    );
    assert_eq!(TableEnd::Move(i32::MAX).apply(&mut table, 40), Err(()));
}

/// A sequence without masks on a new table (`DiagonalWalk::plain`, item M1 of M/m4b3b/BRIEF_FIX3.md):
/// its end moves the offset past the entries that no sequence wrote where they lie above its own,
/// so that the worker's next subject meets none of them.
#[test]
fn a_plain_sequence_on_a_new_table_ends_past_the_unwritten_entries() {
    let at = |offset, cleared| NcbiStart { offset, cleared };
    // One-hit at NCBI's offset -98 (after a walk that went back): the entries lie 98 above the
    // start, 48 past the end of 50 letters.
    let walk = DiagonalWalk::plain(3, 0, 50, at(-98, false));
    assert_eq!(walk.after(1000), TableEnd::Move(48));
    // Longer than that: its own entries reach further, past which nothing remains.
    assert_eq!(
        DiagonalWalk::plain(3, 0, 120, at(-98, false)).after(1000),
        TableEnd::Keep
    );
    // Two-hit at -60, 50 letters: the entries at 60 lie 10 past the end less the window.
    let walk = DiagonalWalk::plain(3, 40, 50, at(-60, false));
    assert_eq!(walk.after(1000), TableEnd::Move(10));
    // At a positive offset the entries lie below the start: nothing to move past.
    assert_eq!(
        DiagonalWalk::plain(3, 40, 50, at(40, false)).after(1000),
        TableEnd::Keep
    );
    assert_eq!(
        DiagonalWalk::plain(3, 0, 50, at(0, true)).after(1000),
        TableEnd::Keep
    );
    // A move that would reach INT4_MAX / 4 clears the table.
    let walk = DiagonalWalk::plain(3, 0, 50, at(-98, false));
    assert_eq!(walk.after(i32::MAX / 4 - 48), TableEnd::Clear);
}

/// Which sequences BLASTP's search of several query batches keeps apart, on NCBI's tables of each
/// batch (`DiagonalWalk::apart`): one on a new table, at NCBI's start, and one whose walk reads
/// back, at the window plus its margin below its start; any other shares the search's table.
#[test]
fn sequences_whose_hits_meet_other_entries_are_kept_apart() {
    let first = NcbiStart {
        offset: 7,
        cleared: true,
    };
    let apart = hard(&[(0, 40), (100, -50), (60, 180)]);
    assert_eq!(
        DiagonalWalk::apart(3, 0, 240, &apart, Some(first)),
        Some(first)
    );
    assert_eq!(
        DiagonalWalk::apart(3, 0, 240, &[], Some(first)),
        Some(first)
    );
    // nest, one-hit: the last chunk 10 below the start, within reach of the chunk before it.
    assert_eq!(
        DiagonalWalk::apart(3, 0, 240, &apart, None),
        Some(NcbiStart {
            offset: 10,
            cleared: false
        })
    );
    // Two-hit: the last chunk at 70 within the window of the one word's entries at 80.
    assert_eq!(
        DiagonalWalk::apart(3, 40, 240, &apart, None),
        Some(NcbiStart {
            offset: 40,
            cleared: false
        })
    );
    // Masks apart in the order of their start: no chunk reads back.
    let forward = hard(&[(0, 40), (60, 50), (130, 110)]);
    assert_eq!(DiagonalWalk::apart(3, 40, 240, &forward, None), None);
    assert_eq!(DiagonalWalk::apart(3, 40, 240, &[], None), None);
}
