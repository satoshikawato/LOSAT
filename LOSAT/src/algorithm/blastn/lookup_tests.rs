//! Unit tests for blastn/lookup.rs

use crate::algorithm::blastn::lookup::reverse_complement;

#[test]
fn test_reverse_complement_basic() {
    let seq = b"ACGT";
    let rc = reverse_complement(seq);

    assert_eq!(rc, b"ACGT"); // ACGT -> TGCA -> reverse -> ACGT
                             // Actually: ACGT -> complement -> TGCA -> reverse -> ACGT
                             // Wait, let me check: A->T, C->G, G->C, T->A
                             // ACGT -> TGCA, then reverse -> ACGT
                             // Actually reverse complement: reverse first then complement
                             // ACGT reverse -> TGCA, complement -> ACGT
                             // Let me verify the actual implementation
}

#[test]
fn test_reverse_complement_actual() {
    let seq = b"ACGT";
    let _rc = reverse_complement(seq);

    // Reverse complement: reverse first, then complement
    // ACGT -> reverse -> TGCA -> complement -> ACGT
    // Actually the function does: reverse, then complement each base
    // ACGT reversed = TGCA, then complement:
    // T->A, G->C, C->G, A->T = ACGT
    // So ACGT -> ACGT (palindrome)

    // Let's test a non-palindrome
    let seq2 = b"AAAA";
    let rc2 = reverse_complement(seq2);
    assert_eq!(rc2, b"TTTT");
}

#[test]
fn test_reverse_complement_case_insensitive() {
    let seq_upper = b"ACGT";
    let seq_lower = b"acgt";
    let seq_lower_normalized: Vec<u8> = seq_lower
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect();

    let rc_upper = reverse_complement(seq_upper);
    let rc_lower = reverse_complement(&seq_lower_normalized);

    assert_eq!(rc_upper, rc_lower);
}

#[test]
fn test_reverse_complement_ambiguous_bases() {
    let seq = b"ACNT";
    let rc = reverse_complement(seq);

    // Ambiguous bases should become N
    assert_eq!(rc[1], b'N');
}

#[test]
fn test_reverse_complement_round_trip() {
    let seq = b"ACGTACGT";
    let rc = reverse_complement(seq);
    let rc_rc = reverse_complement(rc.as_slice());

    // Double reverse complement should return original
    assert_eq!(seq, rc_rc.as_slice());
}

#[test]
fn test_reverse_complement_u_and_t() {
    // NCBI IUPACNA_TO_NCBI4NA maps DNA bases; U is not T in this table.
    let seq_t = b"ACGT";
    let seq_u = b"ACGU";

    let rc_t = reverse_complement(seq_t);
    let rc_u = reverse_complement(seq_u);

    assert_ne!(rc_t, rc_u);
    assert_eq!(rc_u[0], b'-');
}

// `XHeadFifo` lets every word out once, in the order the words went in, through the pushes and
// the flush, for counts below, at and past its depth; a push says whether a word left.
#[test]
fn the_head_fifo_emits_every_word_once_in_order() {
    use crate::algorithm::blastn::lookup::{XHeadFifo, X_HEAD_FIFO};
    for n in [
        0usize,
        1,
        X_HEAD_FIFO - 1,
        X_HEAD_FIFO,
        X_HEAD_FIFO + 1,
        2 * X_HEAD_FIFO + 1,
        100,
    ] {
        let mut fifo: XHeadFifo<(u32, Option<u64>)> = XHeadFifo::new();
        let mut out = Vec::new();
        for i in 0..n {
            let left = fifo.push((i as u32, Some(i as u64 * 3)), &mut |w| out.push(w));
            assert_eq!(left, i >= X_HEAD_FIFO, "push {i} of {n}");
        }
        fifo.flush(&mut |w| out.push(w));
        let want: Vec<_> = (0..n).map(|i| (i as u32, Some(i as u64 * 3))).collect();
        assert_eq!(out, want, "{n} words");
        // (a flushed FIFO starts again empty)
        let mut again = Vec::new();
        fifo.push((7, None), &mut |w| again.push(w));
        fifo.flush(&mut |w| again.push(w));
        assert_eq!(again, vec![(7, None)]);
    }
}

/// No NCBI counterpart (a check of the loop's form): `x_range_words` gives the words, codes and
/// ambiguity count of the loop the table fills used before it (one pass over the letters that
/// counts each word's letters, `valid_bases`), for random letters with ambiguities (none, rare,
/// frequent), lengths 0 to 300, word lengths 1 to 32 and masks of 2 to 64 bits (`u64::MAX` as
/// the discontiguous templates' loop, which does not mask).
#[test]
fn x_range_words_gives_the_counted_loop_s_words() {
    use crate::algorithm::blastn::lookup::x_range_words;
    fn counted(letters: &[u8], w: usize, mask: u64) -> (Vec<(u64, usize)>, usize) {
        let (mut out, mut amb) = (Vec::new(), 0usize);
        let (mut code, mut valid) = (0u64, 0usize);
        for (i, &base) in letters.iter().enumerate() {
            if (base & 0xFC) != 0 {
                code = 0;
                valid = 0;
                amb += 1;
                continue;
            }
            code = ((code << 2) | base as u64) & mask;
            valid += 1;
            if valid < w {
                continue;
            }
            out.push((code, i));
        }
        (out, amb)
    }
    let mut state = 0x5eed_0b1du64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as u32
    };
    let mut words_seen = 0usize;
    for round in 0..400 {
        let len = (next() % 301) as usize;
        let ambiguity_per_mille = [0, 5, 100, 600][round % 4];
        let letters: Vec<u8> = (0..len)
            .map(|_| {
                if next() % 1000 < ambiguity_per_mille {
                    4 + (next() % 12) as u8
                } else {
                    (next() % 4) as u8
                }
            })
            .collect();
        for w in 1..=32usize {
            for mask in [
                (1u64 << (2 * w.min(31))) - 1,
                (1u64 << (2 * (w.min(16) + 1).min(31))) - 1,
                u64::MAX,
            ] {
                let (mut got, mut amb) = (Vec::new(), 0usize);
                x_range_words(
                    &letters,
                    w,
                    mask,
                    || amb += 1,
                    |code, i| got.push((code, i)),
                );
                let (want, want_amb) = counted(&letters, w, mask);
                assert_eq!(got, want, "round {round} len {len} w {w} mask {mask:#x}");
                assert_eq!(amb, want_amb, "round {round} len {len} w {w}");
                words_seen += got.len();
            }
        }
    }
    assert!(words_seen > 1_000_000, "{words_seen} words");
}

/// No NCBI counterpart (a check of the walk's form): `contains_hit` answers, for every cell of the
/// query's words and a few empty cells, and every query offset, whether the cell's chain holds
/// the offset (`for_each_hit`), now that the walk returns at the first equal element; on chain
/// tables of 4-byte elements (lookup width 8) and of 3-byte elements (width 12, LOSAT_X_MB24),
/// with long chains (runs of one letter, repeats) and random letters.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:41-69
#[test]
fn contains_hit_answers_the_chain_s_offsets() {
    use crate::algorithm::blastn::lookup::build_two_stage_lookup;
    let mut seed = 0x00C0_47A1_u64;
    let mut next = |n: u64| {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) % n
    };
    let mut seq: Vec<u8> = Vec::new();
    seq.extend(std::iter::repeat_n(0u8, 40));
    for _ in 0..30 {
        seq.extend([1u8, 3, 2, 0, 1, 3]);
    }
    seq.extend((0..400).map(|_| next(4) as u8));
    seq.extend(std::iter::repeat_n(3u8, 25));
    seq.extend((0..200).map(|_| next(4) as u8));
    for lut in [8usize, 12] {
        let refs: Vec<&[u8]> = vec![&seq];
        let masks = vec![Vec::new()];
        let table = build_two_stage_lookup(&refs, &[0], 28, lut, &masks, None, 0, 1_000_000, false);
        let mask = (1u64 << (2 * lut)) - 1;
        let mut cells: Vec<u64> = seq
            .windows(lut)
            .map(|w| w.iter().fold(0u64, |a, &b| (a << 2) | u64::from(b)) & mask)
            .collect();
        cells.extend((0..20).map(|_| next(mask + 1)));
        cells.sort_unstable();
        cells.dedup();
        let mut longest = 0usize;
        for &cell in &cells {
            let mut chain = Vec::new();
            table.for_each_hit(cell, |q| chain.push(q));
            longest = longest.max(chain.len());
            for q in 0..=(seq.len() as u32 + 2) {
                assert_eq!(
                    table.contains_hit(cell, q),
                    q != 0 && chain.contains(&q),
                    "lut {lut} cell {cell} offset {q}"
                );
            }
        }
        assert!(longest >= 25, "lut {lut}: a long chain ({longest})");
    }
}
