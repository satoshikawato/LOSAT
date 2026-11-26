use super::*;
use crate::config::ScoringMatrix;
use crate::utils::matrix::{ncbistdaa, protein_score_rows};

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

/// The seeds that `test` saves of the hits `hits` on `entries` (mask 63, offset 1000), with the
/// query offset of each.
fn saved_seeds(
    test: &OneWord<'_>,
    hits: &[i32],
    entries: &mut [DiagStruct],
    context: impl Fn(i32) -> WordContext,
) -> Vec<(i32, OneWordSeed)> {
    let mut saved = Vec::new();
    test.test(hits, entries, 63, 1000, context, |q, seed| {
        saved.push((q, seed))
    });
    saved
}

/// The hit at query offset 4 on the diagonal table (mask 63, offset 1000) with `entry` on its
/// diagonal: whether it saved a seed, and the entry after it.
fn test_hit(entry: DiagStruct, cutoff: i32) -> (bool, DiagStruct) {
    let query = query();
    let letters = [ncbistdaa::W; 3];
    let mut entries = vec![DiagStruct::default(); 64];
    entries[4] = entry;
    let saved = saved_seeds(
        &one_word(&query, &letters),
        &[4],
        &mut entries,
        context(cutoff),
    );
    (!saved.is_empty(), entries[4])
}

/// The seed of `W W W` at query offset 4 and subject offset 0 (`one_word_extension`).
const WWW_SEED: OneWordSeed = OneWordSeed {
    score: 33,
    right_extend: false,
    s_last_off: 3,
    q_start: 4,
    s_start: 0,
    length: 3,
};

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
    // At the cutoff the seed is saved (the chunk's ungapped hit, for its gapped stage), and the
    // hit becomes the last hit as after any extension that stayed left.
    let query = query();
    let letters = [ncbistdaa::W; 3];
    let mut entries = vec![DiagStruct::default(); 64];
    entries[4] = entry_behind(10, 1000);
    let saved = saved_seeds(&one_word(&query, &letters), &[4], &mut entries, context(33));
    assert_eq!(saved, [(4, WWW_SEED)]);
    assert_eq!((entries[4].last_hit(), entries[4].flag()), (1000, 0));
    // A saved seed does not end the test of the word's hits: the next hit, on its own diagonal,
    // is tested too (here 10 behind as well).
    entries[4] = entry_behind(10, 1000);
    entries[5] = entry_behind(10, 1000);
    let saved = saved_seeds(
        &one_word(&query, &letters),
        &[4, 5],
        &mut entries,
        context(1),
    );
    assert_eq!(saved.len(), 2);
    assert_eq!(saved[0], (4, WWW_SEED));
    assert_eq!(saved[1].0, 5);
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
    assert!(saved_seeds(&test, &[4], &mut entries, in_context).is_empty());
    assert_eq!((entries[4].last_hit(), entries[4].flag()), (1000, 0));
    // 10 behind: query offset -6, a large unsigned offset: extended, and the seed is saved.
    entries[4] = entry_behind(10, 1000);
    assert_eq!(
        saved_seeds(&test, &[4], &mut entries, in_context),
        [(4, WWW_SEED)]
    );
}

/// The extension of the one word: NCBI's two-hit extension with the chunk's length as the bound
/// of the right walk; the first hit's end is signed. The ungapped hit is the word's three letters
/// (the left walk's best start, at the chunk's offset, to the first walk's best end, past the
/// chunk of one letter: the right walk reads nothing).
#[test]
fn one_word_extension() {
    let query = query();
    let letters = [ncbistdaa::W; 3];
    let test = one_word(&query, &letters);
    assert_eq!(test.extend(-7, 4, 16), WWW_SEED);
    assert_eq!(
        test.extend(0, 4, 16),
        OneWordSeed {
            right_extend: true,
            ..WWW_SEED
        }
    );
}

#[test]
fn short_chunk_refusals_name_the_oid_and_pass_as_is() {
    let error = short_chunk_tree_refusal(17, "BLASTX", -5);
    assert!(error.is::<ShortChunkSeedRefusal>());
    let text = error.to_string();
    assert!(
        text.starts_with(
            "-db_hard_mask: the sequence of OID 17 has a part shorter than the word after \
             overlapping masks, of -5 letters to NCBI BLAST+ 2.17.0"
        ),
        "{text}"
    );
    assert!(text.ends_with("LOSAT's BLASTX does not search such a database with -db_hard_mask"));
    let error = short_chunk_merge_refusal(17, "BLASTP");
    assert!(error.is::<ShortChunkSeedRefusal>());
    let text = error.to_string();
    assert!(
        text.starts_with(
            "-db_hard_mask: the sequence of OID 17 has a part shorter than the word after \
             overlapping masks, whose alignments"
        ),
        "{text}"
    );
    assert!(
        text.contains("which LOSAT's BLASTP does not do yet"),
        "{text}"
    );
    for (error, prefix, suffix) in [
        (
            chunk_merge_undefined_refusal(17, "BLASTX"),
            "-db_hard_mask: the sequence of OID 17 has alignments across a chunk overlap whose \
             merged score calculation is undefined in NCBI BLAST+",
            "LOSAT's BLASTX does not search such a database with -db_hard_mask",
        ),
        (
            last_range_overlap_merge_refusal(17, "BLASTP"),
            "-db_hard_mask: the sequence of OID 17 has alignments from an earlier range that meet \
             the final part across an overlap NCBI BLAST+ 2.17.0 reads past its list of ranges",
            "LOSAT's BLASTP does not search such a database with -db_hard_mask",
        ),
    ] {
        assert!(error.is::<ShortChunkSeedRefusal>());
        let text = error.to_string();
        assert!(text.starts_with(prefix), "{text}");
        assert!(text.ends_with(suffix), "{text}");
    }
}

/// NCBI's merge after a chunk runs its joins where an HSP of the joined list ends past the
/// chunk's offset and one of the chunk's starts before its offset plus the overlap (any start
/// where the overlap is read past the list of ranges).
#[test]
fn merges_that_join_a_chunk_with_the_list() {
    let chunk = |overlap, overlap_past_list| crate::core::subject_split::HardChunk {
        offset: 100,
        length: 2,
        residual: 0,
        full_range: false,
        range: 2,
        overlap,
        overlap_past_list,
    };
    let joins = |chunk: crate::core::subject_split::HardChunk, ends: &[i32], starts: &[i32]| {
        merge_joins(&chunk, ends.iter().copied(), starts.iter().copied())
    };
    // Earlier HSPs that end at the offset or before it: no join.
    assert!(!joins(chunk(100, false), &[40, 100], &[100]));
    // One that ends past it, with a start within the overlap: a join; past the overlap: none.
    assert!(joins(chunk(100, false), &[40, 101], &[150]));
    assert!(!joins(chunk(100, false), &[101], &[200]));
    // No overlap: the chunk's HSPs start at its offset or after it.
    assert!(!joins(chunk(0, false), &[101], &[100]));
    // The overlap read past the list: any start.
    assert!(joins(chunk(0, true), &[101], &[10_000]));
    // No HSP on one side.
    assert!(!joins(chunk(100, false), &[], &[100]));
    assert!(!joins(chunk(100, false), &[101], &[]));
}

/// `DiagCell` and `x_test_hit` on `(i32, bool)` cells give the exact same result and cell state as
/// on `DiagStruct` cells (flag set and clear, entries above and below the hit, and `OneWord` on a
/// `(i32, bool)` `DiagTable`).
#[test]
fn diag_cell_x_test_hit_matches_between_diag_struct_and_pair() {
    use crate::algorithm::common::diag_entry::{x_test_hit, DiagCell, DiagStruct, DiagTable};

    struct PairRing<'a> {
        entries: &'a mut [(i32, bool)],
        mask: i32,
    }

    impl DiagTable for PairRing<'_> {
        type Cell = (i32, bool);

        fn two_hit_entry(&mut self, q: i32, s: i32) -> (&mut (i32, bool), i32) {
            (
                &mut self.entries[(q.wrapping_sub(s) & self.mask) as usize],
                0,
            )
        }

        fn one_hit_entry(&mut self, q: i32, s: i32) -> &mut (i32, bool) {
            &mut self.entries[(s.wrapping_sub(q) & self.mask) as usize]
        }
    }

    // Test every combination of flag (0/1) and entry relative to subject_offset + diag_offset:
    // entry above hit (diff < 0), overlap (0 <= diff < word), candidate (word <= diff < window),
    // and too far (diff >= window).
    for flag in [0u32, 1u32] {
        for diff in [-50i32, -1, 0, 2, 3, 10, 39, 40, 100] {
            let offset = 1000i32;
            let s = 5i32;
            let last_hit = s + offset - diff;
            let mut diag = DiagStruct::default();
            diag.set_last_hit(last_hit);
            diag.set_flag(flag);
            let mut pair = (last_hit, flag != 0);
            let a = x_test_hit(&mut diag, s, offset, 40, 3);
            let b = x_test_hit(&mut pair, s, offset, 40, 3);
            assert_eq!(a, b, "flag={flag} diff={diff}");
            assert_eq!(
                (diag.last_hit(), diag.flag()),
                (pair.last_hit(), pair.flag()),
                "flag={flag} diff={diff}"
            );
        }
    }

    // `OneWord::test_on` on a `(i32, bool)` table matches `DiagStruct` on a candidate and a
    // flagged entry ahead of the hit.
    let query = query();
    let letters = [ncbistdaa::W; 3];
    let test = one_word(&query, &letters);
    let mut diag_entries = vec![DiagStruct::default(); 64];
    diag_entries[4] = entry_behind(10, 1000);
    let mut pair_entries = vec![(0i32, false); 64];
    pair_entries[4] = (990, false);

    let mut diag_saved = Vec::new();
    test.test(&[4], &mut diag_entries, 63, 1000, context(33), |q, seed| {
        diag_saved.push((q, seed));
    });
    let mut pair_saved = Vec::new();
    let mut ring = PairRing {
        entries: &mut pair_entries,
        mask: 63,
    };
    test.test_on(&mut ring, &[4], 1000, context(33), |q, seed| {
        pair_saved.push((q, seed));
    });
    assert_eq!(pair_saved, diag_saved);
    assert_eq!(
        (pair_entries[4].last_hit(), pair_entries[4].flag()),
        (diag_entries[4].last_hit(), diag_entries[4].flag())
    );
}
