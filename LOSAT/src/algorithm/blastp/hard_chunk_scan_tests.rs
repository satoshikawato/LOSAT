use super::*;
use crate::algorithm::blastp::extension::extend_two_hit;
use crate::config::ScoringMatrix;
use crate::utils::matrix::{ncbistdaa, protein_score_rows};

/// The work, the margins, the end and the first word past the end of a walk over `chunks`
/// (`(offset, length)`) of a sequence of `length` letters.
fn walk(
    chunks: &[(usize, i32)],
    length: usize,
    word: i32,
    window: i32,
) -> (Vec<ChunkWork>, DiagonalMargins, i64, Option<usize>) {
    let mut walk = DiagonalWalk::new(word, window, length);
    let works = chunks
        .iter()
        .map(|&(offset, length)| walk.chunk(offset, length))
        .collect();
    (works, walk.margins(), walk.end(), walk.word_past_end())
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
    // chunk's entries (38 at most) but stays above the start.
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
    // lk: masks [40, 200) [60, 238). Two-hit: D = 0, 80, -20 and the end 22: the last chunk reaches
    // 20 below the start, and the one word at 80 leaves an entry 98 past the end less the window.
    let lk = [(0, 40), (200, -140), (238, 2)];
    let (works, margins, end, _) = walk(&lk, 240, 3, 40);
    assert_eq!(works, [WORDS, ChunkWork::OneWord, ChunkWork::OneWord]);
    assert_eq!(
        margins,
        DiagonalMargins {
            below_start: 20,
            past_end: 98
        }
    );
    assert_eq!(end, 22);
    // One-hit: D = 0, 40, -100 and the end -98: the first chunk's entries (38) lie 136 past it.
    let (works, margins, end, _) = walk(&lk, 240, 3, 0);
    assert_eq!(works, [WORDS, ChunkWork::Nothing, ChunkWork::Nothing]);
    assert_eq!(
        margins,
        DiagonalMargins {
            below_start: 0,
            past_end: 136
        }
    );
    assert_eq!(end, -98);
}

/// Masks sorted by their start that do not overlap: every offset is at least the window past the
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
            left = right + next(30) as usize;
        }
        let word = 1 + next(7) as i32;
        let window = [0, 1, 2, 40, 1000][next(5) as usize];
        let (works, margins, _, past) = walk(&chunks, length, word, window);
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

/// The one word of a chunk at the sequence's last letter runs past its sentinel: refused where an
/// earlier entry within the window behind it would extend its hits (masks [10, 60) [20, 99)
/// [30, 40) of a sequence of 100 letters), not where none does.
#[test]
fn words_past_the_end_that_can_change_a_result() {
    let chunks = [(0, 10), (60, -40), (99, -69), (40, 60)];
    let (works, margins, _, past) = walk(&chunks, 100, 3, 40);
    assert_eq!(
        works,
        [WORDS, ChunkWork::OneWord, ChunkWork::OneWord, WORDS_BACK]
    );
    assert!(!margins.stale());
    assert_eq!(past, Some(99));
    // The same word with no entry within the window behind it and no chunk after it.
    let (_, _, _, past) = walk(&[(0, 10), (99, -69)], 100, 3, 40);
    assert_eq!(past, None);
    // A word that ends at the sentinel reads no other sequence.
    let (_, _, _, past) = walk(&[(0, 10), (60, -40), (98, -68), (40, 60)], 100, 3, 40);
    assert_eq!(past, None);
}

/// The table is cleared before the sequence where the offset of a chunk that is not the last
/// reaches `INT4_MAX / 4`, never for the last chunk (whose exit clears it as after any subject).
#[test]
fn clears_before_a_sequence_whose_inner_chunk_reaches_the_limit() {
    let limit = i32::MAX / 4;
    let mut lk = DiagonalWalk::new(3, 40, 240);
    for (offset, length) in [(0, 40), (200, -140), (238, 2)] {
        lk.chunk(offset, length);
    }
    // The inner chunks are at 0 and 80.
    assert!(lk.clears_inside(limit - 80));
    assert!(!lk.clears_inside(limit - 81));
    let mut one = DiagonalWalk::new(3, 40, 240);
    one.chunk(0, 240);
    assert!(!one.clears_inside(i32::MAX));
}

/// The query of the one-word tests: `W W W` at query offset 4 among `A`s.
fn query() -> Vec<u8> {
    let (a, w) = (ncbistdaa::A, ncbistdaa::W);
    vec![a, a, a, a, w, w, w, a, a, a]
}

/// A one word `W W W` of a chunk of one letter (word 3, window 40) against [`query`].
fn one_word<'a>(query: &'a [u8], letters: &'a [u8]) -> OneWord<'a> {
    OneWord {
        query,
        letters,
        length: 1,
        word: 3,
        window: 40,
        scores: protein_score_rows(ScoringMatrix::Blosum62),
    }
}

fn context(cutoff: i32) -> impl Fn(i32) -> WordContext {
    move |_: i32| WordContext {
        start: 0,
        x_drop: 16,
        cutoff,
    }
}

/// An entry whose last hit is `diff` before subject offset 0 of a table at offset `offset`.
fn entry_behind(diff: i32, offset: i32) -> DiagStruct {
    let mut entry = DiagStruct::default();
    entry.set_last_hit(offset - diff);
    entry
}

/// The hit at query offset 4 on the diagonal table (mask 63, offset 1000) with `entry` on its
/// diagonal: whether it saved a seed, and the entry after it.
fn test_hit(entry: DiagStruct, cutoff: i32) -> (bool, DiagStruct) {
    let query = query();
    let letters = [ncbistdaa::W; 3];
    let mut entries = vec![DiagStruct::default(); 64];
    entries[4] = entry;
    let saved = one_word(&query, &letters).test(&[4], &mut entries, 63, 1000, context(cutoff));
    (saved, entries[4])
}

#[test]
fn one_word_hits_on_the_diagonal_table() {
    // A first hit a window or more behind: the hit becomes the last hit, nothing is extended.
    let (saved, entry) = test_hit(entry_behind(40, 1000), 1);
    assert!(!saved);
    assert_eq!((entry.last_hit(), entry.flag()), (1000, 0));
    // A first hit that overlaps the word: nothing changes.
    let (saved, entry) = test_hit(entry_behind(2, 1000), 1);
    assert!(!saved);
    assert_eq!(entry, entry_behind(2, 1000));
    // A first hit 10 behind: the left walk does not reach its end, so the extension stays left
    // (score 33), and the hit becomes the last hit.
    let (saved, entry) = test_hit(entry_behind(10, 1000), 34);
    assert!(!saved);
    assert_eq!((entry.last_hit(), entry.flag()), (1000, 0));
    // At the cutoff the seed is saved: LOSAT's refusal.
    assert!(test_hit(entry_behind(10, 1000), 33).0);
    // A first hit 3 behind: the left walk reaches it, the right walk reads nothing of a chunk of
    // one letter (its start is past the chunk), and the entry takes the flag at the walk's last
    // offset less the word plus one.
    let (saved, entry) = test_hit(entry_behind(3, 1000), 34);
    assert!(!saved);
    assert_eq!((entry.last_hit(), entry.flag()), (1001, 1));
    // An entry with the flag ahead of the hit: skipped; at or behind it: the hit's own value.
    let mut flagged = DiagStruct::default();
    flagged.set_flag(1);
    flagged.set_last_hit(1001);
    assert_eq!(test_hit(flagged, 1), (false, flagged));
    flagged.set_last_hit(1000);
    let (saved, entry) = test_hit(flagged, 1);
    assert!(!saved);
    assert_eq!((entry.last_hit(), entry.flag()), (1000, 0));
}

/// The context test subtracts the distance to the last hit from the unsigned query offset: a
/// last hit before the query's start passes it (and the hit is extended), one in the context
/// before the hit's does not.
#[test]
fn one_word_context_test_is_unsigned() {
    let query = query();
    let letters = [ncbistdaa::W; 3];
    let test = one_word(&query, &letters);
    // The hit at query offset 4 of a context that starts at 2.
    let in_context = |_: i32| WordContext {
        start: 2,
        x_drop: 16,
        cutoff: 1,
    };
    // 3 behind: query offset 1, before the context's start: the hit becomes the last hit.
    let mut entries = vec![DiagStruct::default(); 64];
    entries[4] = entry_behind(3, 1000);
    assert!(!test.test(&[4], &mut entries, 63, 1000, in_context));
    assert_eq!((entries[4].last_hit(), entries[4].flag()), (1000, 0));
    // 10 behind: query offset -6, a large unsigned offset: extended, and the seed is saved.
    entries[4] = entry_behind(10, 1000);
    assert!(test.test(&[4], &mut entries, 63, 1000, in_context));
}

/// The extension of the one word: NCBI's two-hit extension with the chunk's length as the bound
/// of the right walk; the first hit's end is signed.
#[test]
fn one_word_extension() {
    let query = query();
    let letters = [ncbistdaa::W; 3];
    let test = one_word(&query, &letters);
    assert_eq!(
        test.extend(-7, 4, 16),
        OneWordSeed {
            q_start: 4,
            s_start: 0,
            length: 3,
            score: 33,
            right_extend: false,
            s_last_off: 3,
        }
    );
    assert_eq!(
        test.extend(0, 4, 16),
        OneWordSeed {
            q_start: 4,
            s_start: 0,
            length: 3,
            score: 33,
            right_extend: true,
            s_last_off: 3,
        }
    );
}

/// Where the chunk holds the word, the extension is LOSAT's reference two-hit extension of the
/// chunk's letters (`extension::extend_two_hit`) at subject offset 0, for any first hit's end,
/// negative ones included.
#[test]
fn one_word_extension_is_the_reference_two_hit_extension() {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % below
    };
    let mut extended = 0;
    for _ in 0..20_000 {
        let word = 2 + next(4) as i32;
        let length = word as usize + next(8) as usize;
        let letters: Vec<u8> = (0..length).map(|_| 1 + next(20) as u8).collect();
        let query: Vec<u8> = (0..word as usize + next(12) as usize)
            .map(|_| 1 + next(20) as u8)
            .collect();
        let q = next((query.len() - word as usize + 1) as u64) as i32;
        let s_left_off = next(20) as i32 - 16;
        let x_drop = 1 + next(30) as i32;
        let test = OneWord {
            query: &query,
            letters: &letters,
            length: length as i32,
            word,
            window: 40,
            scores: protein_score_rows(ScoringMatrix::Blosum62),
        };
        let got = test.extend(s_left_off, q, x_drop);
        let want = extend_two_hit(
            ScoringMatrix::Blosum62,
            &query,
            &letters,
            s_left_off,
            0,
            q as usize,
            x_drop,
            word as usize,
        )
        .expect("the word lies within both sequences");
        assert_eq!(
            (
                got.q_start,
                got.s_start,
                got.length,
                got.score,
                got.right_extend,
                got.s_last_off
            ),
            (
                want.ungapped_data.q_start,
                want.ungapped_data.s_start,
                want.ungapped_data.length,
                want.ungapped_data.score,
                want.right_extend,
                want.s_last_off
            ),
            "{query:?} {letters:?} q={q} s_left_off={s_left_off} x_drop={x_drop}"
        );
        extended += usize::from(got.right_extend);
    }
    assert!(extended > 1_000, "only {extended} right walks");
}

#[test]
fn short_chunk_seed_refusal_names_the_oid_and_passes_as_is() {
    let error = short_chunk_seed_refusal(17, "BLASTX");
    assert!(error.is::<ShortChunkSeedRefusal>());
    let text = error.to_string();
    assert!(
        text.starts_with("-db_hard_mask: the sequence of OID 17 has a part between overlapping"),
        "{text}"
    );
    assert!(text.ends_with("LOSAT's BLASTX does not search such a database with -db_hard_mask"));
}
