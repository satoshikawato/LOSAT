use super::*;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::algorithm::tblastx::lookup::compressed::{
    build_blosum62_compressed_lookup, compressed_lookup_parameters,
    compute_encoded_compressed_index,
};

// A small deterministic generator (xorshift64).
struct XorShift(u64);

impl XorShift {
    fn below(&mut self, n: usize) -> usize {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x % n as u64) as usize
    }
}

// NCBISTDAA letters inside the 15-letter compressed alphabet and outside it (NULLB, X, U, *, O).
const LETTERS: [u8; 23] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 27,
];
const OUTSIDE: [u8; 5] = [0, 21, 24, 25, 26];

// The reference: `blastp_for_each_offset_pair`'s rounds of `scan_subject` with an array of
// `array_size` pairs; returns the hits in order and the number of rounds.
fn reference_hits(
    lookup: &BlastCompressedAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    array_size: usize,
) -> (Vec<(i32, i32)>, usize) {
    let seq_ranges = [(0, subject_length as i32)];
    let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - lookup.word_length];
    if scan_range[2] < scan_range[1] {
        scan_range[2] = scan_range[1];
    }
    let mut offset_pairs = vec![OffsetPair::default(); array_size];
    let mut hits = Vec::new();
    let mut rounds = 0;
    while scan_range[1] <= scan_range[2] {
        let count = lookup.scan_subject(
            subject,
            &seq_ranges,
            &mut offset_pairs,
            array_size as i32,
            &mut scan_range,
        );
        hits.extend(
            offset_pairs[..count as usize]
                .iter()
                .map(|pair| (pair.s_off as i32, pair.q_off as i32)),
        );
        rounds += 1;
    }
    (hits, rounds)
}

fn fused_hits(
    lookup: &BlastCompressedAaLookupTable,
    subject: &[u8],
    subject_length: usize,
) -> Vec<(i32, i32)> {
    let mut scratch = XScratch::default();
    let mut hits = Vec::new();
    x_scan_cells(lookup, subject, subject_length, &mut scratch, |s, cell| {
        hits.extend(cell.iter().map(|&q| (s, q)));
    });
    hits
}

// The reference's two-hit test (`blast_engine.rs`, the closure of `blastp_for_each_offset_pair`)
// over `hits`, with `candidate` in place of the context check and the extension.
#[allow(clippy::too_many_arguments)]
fn reference_two_hit(
    hits: &[(i32, i32)],
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    candidate: &mut dyn FnMut(i32, i32, i32, i32, &mut DiagStruct),
) {
    for &(subject_offset, query_offset) in hits {
        let diag_entry = &mut diag_array[((query_offset - subject_offset) & diag_mask) as usize];
        if diag_entry.flag() != 0 {
            let subject_plus_offset = subject_offset + diag_offset;
            if subject_plus_offset < diag_entry.last_hit() {
                continue;
            }
            diag_entry.set_last_hit(subject_plus_offset);
            diag_entry.set_flag(0);
            continue;
        }
        let last_hit = diag_entry.last_hit() - diag_offset;
        let diff = subject_offset - last_hit;
        if diff >= window {
            diag_entry.set_last_hit(subject_offset + diag_offset);
            continue;
        }
        if diff < word_size {
            continue;
        }
        candidate(query_offset, subject_offset, last_hit, diff, diag_entry);
    }
}

// A stand-in for the extension: it logs the call and writes the entry as the reference closure
// does (a flag and a later `last_hit` after a right extension, else this hit).
fn stand_in(
    log: &mut Vec<(i32, i32, i32, i32)>,
    diag_offset: i32,
    query_offset: i32,
    subject_offset: i32,
    last_hit: i32,
    diff: i32,
    entry: &mut DiagStruct,
) {
    log.push((query_offset, subject_offset, last_hit, diff));
    match (query_offset ^ subject_offset) & 3 {
        0 => {
            entry.set_flag(1);
            entry.set_last_hit(subject_offset + 7 + diag_offset);
        }
        1 => {}
        _ => entry.set_last_hit(subject_offset + diag_offset),
    }
}

fn random_subject(rng: &mut XorShift, length: usize, outside_one_in: usize) -> Vec<u8> {
    let mut subject: Vec<u8> = (0..length)
        .map(|_| {
            if outside_one_in != 0 && rng.below(outside_one_in) == 0 {
                OUTSIDE[rng.below(OUTSIDE.len())]
            } else {
                LETTERS[rng.below(LETTERS.len())]
            }
        })
        .collect();
    subject.push(0);
    subject
}

// Runs the subjects one after the other through both loops from the same diagonal table (the
// offset moves on by length + window after each subject, as `finish_subject` does) and compares
// the hits, the candidates and the table.
fn assert_same_two_hit(
    lookup: &BlastCompressedAaLookupTable,
    subjects: &[Vec<u8>],
    diag_len: usize,
    array_sizes: &[usize],
) {
    let window = 40;
    let word_size = lookup.word_length;
    let diag_mask = diag_len as i32 - 1;
    let mut fused_diag = vec![DiagStruct::default(); diag_len];
    let mut reference_diag = vec![DiagStruct::default(); diag_len];
    let mut diag_offset = window;
    let mut total_candidates = 0usize;
    for (n, subject) in subjects.iter().enumerate() {
        let length = subject.len() - 1;
        let fused = fused_hits(lookup, subject, length);
        let array_size = array_sizes[n % array_sizes.len()];
        let (reference, _) = reference_hits(lookup, subject, length, array_size);
        assert_eq!(fused, reference, "subject {n}: hits");

        let mut fused_log = Vec::new();
        x_scan_two_hit_compressed(
            lookup,
            subject,
            length,
            &mut fused_diag,
            diag_mask,
            diag_offset,
            window,
            word_size,
            |q, s, last_hit, diff, entry| {
                stand_in(&mut fused_log, diag_offset, q, s, last_hit, diff, entry)
            },
        );
        let mut reference_log = Vec::new();
        reference_two_hit(
            &reference,
            &mut reference_diag,
            diag_mask,
            diag_offset,
            window,
            word_size,
            &mut |q, s, last_hit, diff, entry| {
                stand_in(&mut reference_log, diag_offset, q, s, last_hit, diff, entry)
            },
        );
        assert_eq!(fused_log, reference_log, "subject {n}: candidates");
        assert!(
            fused_diag
                .iter()
                .zip(&reference_diag)
                .all(|(a, b)| a.raw_bits() == b.raw_bits()),
            "subject {n}: diagonal table"
        );
        total_candidates += fused_log.len();
        diag_offset += length as i32 + window;
    }
    assert!(total_candidates > 0, "the test reaches the extension");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1288-1301
// The resumed rounds of the reference scan give the rolling index only if the reciprocal gives
// the exact quotient for every index the scan forms (module header).
#[test]
fn the_reciprocal_gives_the_exact_quotient_below_every_tables_size() {
    for word_size in 5..=7usize {
        let (alphabet, recip) = compressed_lookup_parameters(word_size).expect("table");
        let bound = (alphabet as i64).pow(word_size as u32);
        for index in 0..bound {
            let index = index as i32;
            assert_eq!(
                compressed_reciprocal_divide(index, i64::from(recip)),
                index / alphabet as i32,
                "word {word_size} index {index}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:825-844
// Cells of 1 to 20 hits (every position of the newest group), each word once in the subject and
// between letters outside the alphabet: the fused loop gives each cell's hits in the order of
// NCBI's scan, with arrays that make the reference stop before many cells.
#[test]
fn cells_of_one_to_twenty_hits_come_in_ncbis_group_order() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    let mut lookup = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
    let mut subject = Vec::new();
    let mut expected = Vec::new();
    for numhits in 1..=20usize {
        let (word, index) = loop {
            let word: Vec<u8> = (0..5).map(|_| LETTERS[rng.below(LETTERS.len())]).collect();
            let encoded: Vec<u8> = word
                .iter()
                .map(|&letter| lookup.compress_table[letter as usize])
                .collect();
            let index = compute_encoded_compressed_index(5, &encoded, 15);
            if lookup.backbone[index].num_used == 0 {
                break (word, index);
            }
        };
        for k in 0..numhits {
            lookup.add_word_hit(index, (1000 * numhits + 37 * k) as i32);
        }
        let order = lookup.query_offsets_in_scan_order(index);
        let s = subject.len() as i32;
        expected.extend(order.iter().map(|&q| (s, q)));
        subject.extend_from_slice(&word);
        subject.push(OUTSIDE[rng.below(OUTSIDE.len())]);
    }
    lookup.finalize();
    assert_eq!(lookup.longest_chain, 20);
    let length = subject.len();
    subject.push(0);
    let fused = fused_hits(&lookup, &subject, length);
    assert_eq!(fused, expected);
    for array_size in [20, 21, 27, 40, 4096 + 20] {
        let (reference, rounds) = reference_hits(&lookup, &subject, length, array_size);
        assert_eq!(fused, reference, "array {array_size}");
        if array_size < 40 {
            assert!(rounds > 5, "array {array_size}: {rounds} rounds");
        }
    }
    // The 8-hit cell of NCBI's example: the first three, the newest group, the older group.
    let eight: Vec<i32> = fused
        .iter()
        .filter(|(_, q)| q / 1000 == 8)
        .map(|&(_, q)| q - 8000)
        .collect();
    assert_eq!(eight, vec![0, 37, 74, 259, 111, 148, 185, 222]);
}

// Queries of few letters and a low threshold give cells of many hits; subjects of random
// letters with letters outside the alphabet (none, rare, frequent, runs), short subjects, and
// arrays from the longest cell up (many rounds) to the reference's size. Small diagonal tables
// make hits of one subject offset share entries, so their order shows in the table.
#[test]
fn the_fused_loop_matches_the_reference_scan_and_two_hit_test() {
    let mut rng = XorShift(0xbb67_ae85_84ca_a73b);
    // Word size 5 only: the tables of 6 and 7 have 15^6 and 10^7 backbone cells.
    let word_size = 5usize;
    for (threshold, letters) in [
        (11.0, &[1u8, 11, 17, 10][..]),
        (9.0, &[1u8, 11, 17][..]),
        (13.0, &LETTERS[..]),
    ] {
        let mut query = Vec::new();
        let mut locations = Vec::new();
        for length in [120usize, 75, 200, 33] {
            let left = query.len() as i32;
            query.extend((0..length).map(|_| letters[rng.below(letters.len())]));
            locations.push((left, left + length as i32 - 1));
            query.push(0);
        }
        let lookup = build_blosum62_compressed_lookup(word_size, threshold, &query, &locations)
            .expect("compressed lookup")
            .expect("within the overflow banks");
        assert!(lookup.longest_chain > 5, "cells with overflow lists");
        let longest = lookup.longest_chain as usize;
        let mut subjects = Vec::new();
        for outside_one_in in [0usize, 50, 6, 2] {
            for _ in 0..6 {
                let length = 1 + rng.below(300);
                let mut subject = random_subject(&mut rng, length, outside_one_in);
                // Some subjects copy a stretch of the query (many hits, many extensions).
                if length > 60 && rng.below(2) == 0 {
                    let from = rng.below(query.len() - 60);
                    let to = rng.below(length - 50);
                    subject[to..to + 50].copy_from_slice(&query[from..from + 50]);
                }
                subjects.push(subject);
            }
        }
        for length in 1..=word_size + 1 {
            subjects.push(random_subject(&mut rng, length, 0));
        }
        let array_sizes = [
            longest,
            longest + 1,
            longest + 5,
            2 * longest,
            4096 + longest,
        ];
        for diag_len in [64usize, 1024, 4096] {
            assert_same_two_hit(&lookup, &subjects, diag_len, &array_sizes);
        }
    }
}
