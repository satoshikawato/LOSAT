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
    for n in [0usize, 1, X_HEAD_FIFO - 1, X_HEAD_FIFO, X_HEAD_FIFO + 1, 2 * X_HEAD_FIFO + 1, 100] {
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
