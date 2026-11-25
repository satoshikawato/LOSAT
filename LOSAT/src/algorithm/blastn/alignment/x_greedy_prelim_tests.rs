use super::*;
use crate::core::blast_encoding::{
    encode_subject_ncbi2na_packed, iupacna_to_ncbi4na, resolve_ncbi4na_to_ncbi2na,
};

/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_util.h:52-55
/// The one-byte copy of a packed subject (IUPAC letters resolved as the engine packs them)
/// holds the resolved base at every position, and `XPrelimSubject::chunk` returns the tail
/// that starts at a multiple of 4.
#[test]
fn x_unpack_ncbi2na_gives_the_resolved_bases() {
    const LETTERS: &[u8] = b"ACGTacgtNRYKMSWBDHVnrykmswbdhv-";
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    for len in [1usize, 3, 4, 5, 63, 64, 65, 1000] {
        let text: Vec<u8> = (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                LETTERS[(state % LETTERS.len() as u64) as usize]
            })
            .collect();
        let packed = encode_subject_ncbi2na_packed(&text);
        let codes = resolve_ncbi4na_to_ncbi2na(
            &text
                .iter()
                .map(|&letter| iupacna_to_ncbi4na(letter).unwrap_or(15))
                .collect::<Vec<u8>>(),
        );
        assert_eq!(x_unpack_ncbi2na(&packed, len), codes, "length {len}");
        let subject = XPrelimSubject::new(true, &packed, len);
        if len >= 4 {
            assert_eq!(subject.chunk(4).expect("enabled"), &codes[4..]);
        }
        assert!(XPrelimSubject::new(false, &packed, len).chunk(0).is_none());
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:516-520,1651-1684
/// A sequence of score-only calls through `XPrelimSubject::score_only` (windows kept across
/// calls, the whole copy after a fence) gives each call's result on the whole copy and leaves
/// `max_dist` where the whole-copy calls leave it: hits of a query made of mutated pieces of
/// the subject (alignments within the margins and past them), chunks at offset 0 and inside
/// the subject, hits at the ends of chunks, linear and affine scores.
#[test]
fn score_only_on_windows_gives_the_whole_copy_result() {
    const SUBJECT_LEN: usize = 30_000;
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let subject: Vec<u8> = (0..SUBJECT_LEN).map(|_| (next() % 4) as u8).collect();
    let text: Vec<u8> = subject.iter().map(|&code| b"ACGT"[code as usize]).collect();
    let packed = encode_subject_ncbi2na_packed(&text);
    // Pieces of the subject, 1-5 % substituted, some with an insertion or a deletion; 12
    // seeds spread over each piece, at the piece's offsets before its indel.
    let mut query = Vec::new();
    let mut seeds = Vec::new();
    for (piece, &(start, len)) in [
        (0usize, 300usize),
        (2_900, 6_000),
        (11_500, 900),
        (11_990, 40),
        (15_000, 2_500),
        (19_900, 200),
        (24_000, 5_999),
        (29_700, 300),
    ]
    .iter()
    .enumerate()
    {
        let q_start = query.len();
        let substitution = 100 / (piece as u64 % 5 + 1);
        let indel_at = if piece % 3 == 1 { len / 2 } else { len };
        for i in 0..len {
            if i == indel_at && piece % 2 == 0 {
                continue;
            }
            if i == indel_at {
                query.push((next() % 4) as u8);
            }
            let base = subject[start + i];
            query.push(if next() % substitution == 0 {
                (base + 1 + (next() % 3) as u8) % 4
            } else {
                base
            });
        }
        query.extend((0..50).map(|_| (next() % 4) as u8));
        let span = indel_at.min(len);
        seeds.extend((0..12).map(|k| (q_start + k * span / 12, start + k * span / 12)));
    }
    let whole = x_unpack_ncbi2na(&packed, SUBJECT_LEN);
    let (mut calls, mut fenced) = (0usize, 0usize);
    for (chunk_offset, subject_len) in [
        (0usize, SUBJECT_LEN),
        (0, 12_000),
        (8_000, 12_000),
        (11_996, 18_004),
    ] {
        for (reward, penalty, gap_open, gap_extend, x_drop) in
            [(1, -2, 0, 0, 20), (1, -3, 0, 0, 40), (2, -3, 5, 2, 30)]
        {
            let subject = XPrelimSubject::new(true, &packed, SUBJECT_LEN);
            let (mut windows, mut copy) = (GreedyAlignScratch::new(), GreedyAlignScratch::new());
            for &(q_off, s_seed) in &seeds {
                if s_seed < chunk_offset || s_seed >= chunk_offset + subject_len {
                    continue;
                }
                let s_off = s_seed - chunk_offset;
                let half = (next() % 400) as usize;
                let hit = (
                    s_off.saturating_sub(half),
                    (s_off + half).min(subject_len - 1),
                );
                let fences = x_fence::x_fence_count();
                let context = ContextQuery {
                    buf: &query,
                    start: 0,
                    len: query.len() as i32,
                };
                let got = subject.score_only(
                    chunk_offset,
                    subject_len,
                    hit,
                    context,
                    q_off as i32,
                    s_off,
                    reward,
                    penalty,
                    gap_open,
                    gap_extend,
                    x_drop,
                    &mut windows,
                );
                fenced += usize::from(x_fence::x_fence_count() != fences);
                calls += 1;
                let want = x_greedy_gapped_alignment_score_only(
                    &whole[chunk_offset..],
                    context,
                    subject_len,
                    q_off as i32,
                    s_off,
                    reward,
                    penalty,
                    gap_open,
                    gap_extend,
                    x_drop,
                    &mut copy,
                );
                assert_eq!(
                    got, want,
                    "chunk {chunk_offset}+{subject_len}, scores {reward}/{penalty} \
                         {gap_open}/{gap_extend}, seed ({q_off}, {s_off})"
                );
                assert_eq!(windows.max_dist, copy.max_dist);
            }
        }
    }
    assert!(
        fenced * 20 > calls && fenced * 2 < calls,
        "calls {calls}, fenced {fenced}"
    );
}
