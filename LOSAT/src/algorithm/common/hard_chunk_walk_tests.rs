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
    assert!(margins.stale());
    // big: masks [40, 120) [60, 130), an overlap longer than the window. Two-hit: D = 0, 80, 60.
    let big = [(0, 40), (120, -60), (130, 110)];
    let (works, margins, _, _) = walk(&big, 240, 3, 40);
    assert_eq!(works, [WORDS, ChunkWork::OneWord, WORDS_BACK]);
    assert!(!margins.stale());
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
/// entries before it, so no chunk goes back, and the margins are 0, at every word and window.
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
        let word = 1 + next(7) as i32;
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
    assert!(!margins.stale());
    assert_eq!(past, Some(99));
    let (_, margins, _, past) = walk(&chunks, 100, 3, 40);
    assert!(!margins.stale());
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

/// The first sequence that a table searches meets NCBI's new table, whose entries lie the window
/// below the sequence's start (`DiagonalWalk::start_on_new_table`): a chunk that goes back below
/// them is read by the reference scan, as a chunk below an earlier chunk's entries is.
#[test]
fn a_new_table_holds_entries_the_window_below_the_start() {
    for (new_table, work) in [(false, WORDS), (true, WORDS_BACK)] {
        let mut walk = DiagonalWalk::new(3, 40, 300, false);
        if new_table {
            walk.start_on_new_table();
        }
        // The one word at 299 covers the sentinel and has no hit: it stores nothing.
        assert_eq!(walk.chunk(299, -100), ChunkWork::OneWord);
        // D = -60, below the new table's entries at -40.
        assert_eq!(walk.chunk(0, 100), work, "{new_table}");
    }
    // One-hit: the entries lie at the start; a chunk at the start reads them as cleared, one below
    // it does not.
    let mut walk = DiagonalWalk::new(3, 0, 300, false);
    walk.start_on_new_table();
    assert_eq!(walk.chunk(0, 10), WORDS);
    let mut walk = DiagonalWalk::new(3, 0, 300, false);
    walk.start_on_new_table();
    assert_eq!(walk.chunk(299, -10), ChunkWork::Nothing);
    assert_eq!(walk.chunk(0, 100), WORDS_BACK);
}
