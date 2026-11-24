use super::*;

// The bound on the appends leaves the record's doubling below `isize::MAX` bytes (a 32-bit
// target takes the `TooLarge` fallback instead of panicking), and is `i32::MAX` where
// `usize` has 64 bits.
#[test]
fn the_appends_bound_keeps_the_records_doubling_inside_isize() {
    for backbone_size in [0, 1, 15usize.pow(5) + 2, 1 << 20] {
        let record = x_max_appends(backbone_size) + 2 * backbone_size;
        assert!(
            record * 2 * core::mem::size_of::<u32>() <= isize::MAX as usize,
            "{backbone_size}"
        );
    }
    #[cfg(target_pointer_width = "64")]
    assert_eq!(x_max_appends(0), i32::MAX as usize);
}
use crate::utils::matrix::ncbistdaa;

// A small deterministic generator (xorshift64).
pub(super) struct XorShift(pub(super) u64);

impl XorShift {
    pub(super) fn below(&mut self, n: usize) -> usize {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x % n as u64) as usize
    }
}

// NCBISTDAA letters inside the 15-letter compressed alphabet and outside it (NULLB, X, U, *, O).
pub(super) const LETTERS: [u8; 23] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 27,
];
const OUTSIDE: [u8; 5] = [0, 21, 24, 25, 26];
pub(super) const FEW_LETTERS: [u8; 4] = [ncbistdaa::A, ncbistdaa::L, ncbistdaa::S, ncbistdaa::K];

// Item SC9i of line NP: `for_each_word_in` gives the cells of `for_each_word` whose last letter
// is in the range, each once (less the range's first cell), for word sizes 5 to 7, matrices of each kind, thresholds that
// make few and many neighbors (and 0: the word alone), words with letters outside the
// compressed alphabet, and every range of letters (all, none, one, any run).
#[test]
fn word_cells_of_a_range_of_last_letters_are_the_words_cells_there() {
    let mut rng = XorShift(0x6a09_e667_f3bc_c908);
    let mut seen = [0usize; 3];
    for (word_size, matrix, thresholds) in [
        (5usize, ScoringMatrix::Blosum62, &[1100, 1930, 0][..]),
        (6, ScoringMatrix::Blosum62, &[1300, 2100][..]),
        (6, ScoringMatrix::Blosum45, &[1300, 2100][..]),
        (7, ScoringMatrix::Blosum62, &[1200, 2025][..]),
        (7, ScoringMatrix::Pam30, &[1200, 2025][..]),
    ] {
        let header = XHeader::new(word_size).expect("word size");
        let scores = matrix
            .score_matrix(header.alphabet)
            .expect("compressed matrix");
        let info = CompressedNeighborInfo::new(&scores, header.alphabet);
        let alphabet = header.alphabet as u32;
        let top = alphabet.pow(word_size as u32 - 1);
        for &threshold in thresholds {
            let words = XWords::new(&header, threshold, &scores, &info);
            for round in 0..60 {
                let letters: &[u8] = match round % 3 {
                    0 => &FEW_LETTERS,
                    1 => &LETTERS,
                    _ => &[
                        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 21, 24, 25,
                    ],
                };
                let word: Vec<u8> = (0..word_size)
                    .map(|_| letters[rng.below(letters.len())])
                    .collect();
                let mut all = Vec::new();
                words.for_each_word(&word[..], &mut all);
                let mut distinct = all.clone();
                distinct.sort_unstable();
                distinct.dedup();
                assert_eq!(distinct.len(), all.len(), "an offset appends a cell twice");
                let a = rng.below(alphabet as usize + 1) as u32;
                let b = rng.below(alphabet as usize + 1) as u32;
                let mut ranges = vec![0..alphabet, a..a, a.min(b)..a.max(b)];
                ranges.extend((0..alphabet).map(|letter| letter..letter + 1));
                for range in ranges {
                    let mut want: Vec<u32> = all
                        .iter()
                        .copied()
                        .filter(|cell| range.contains(&(cell / top)))
                        .map(|cell| cell - range.start * top)
                        .collect();
                    let mut got = Vec::new();
                    words.for_each_word_in(&word[..], &range, &mut got);
                    want.sort_unstable();
                    got.sort_unstable();
                    assert_eq!(
                        got, want,
                        "{matrix:?} word {word_size} threshold {threshold} word {word:?} letters {range:?}"
                    );
                    seen[usize::from(!want.is_empty())
                        + usize::from(want.len() == all.len() && !all.is_empty())] += 1;
                }
            }
        }
    }
    assert!(seen.iter().all(|&n| n > 100), "cases seen {seen:?}");
}

// LOSAT_X_CLUTPAIRS: every matrix of NCBI's with a compressed matrix has a pair table for each
// compressed alphabet (its letters' sums fit the rows' i16); a row's merge gives the pairs in
// the order of the stable sort by sum, and `emit` the cells of the pairs of each key that a
// scan of the sorted pairs gives, at the ends of i32 too.
#[test]
fn pair_tables_hold_every_matrix() {
    for word_size in [5usize, 6, 7] {
        let header = XHeader::new(word_size).expect("word size");
        for matrix in ScoringMatrix::ALL {
            let Some(scores) = matrix.score_matrix(header.alphabet) else {
                assert_eq!(matrix, ScoringMatrix::Identity);
                continue;
            };
            let info = CompressedNeighborInfo::new(&scores, header.alphabet);
            let table = x_pair_table(&header, &info);
            let used = header.alphabet * header.alphabet;
            for a in LETTERS {
                for b in LETTERS {
                    let (a, b) = (usize::from(a), usize::from(b));
                    if usize::from(header.compress_table[a]) >= header.alphabet
                        || usize::from(header.compress_table[b]) >= header.alphabet
                    {
                        continue;
                    }
                    let row = table.row(a, b);
                    let mut want: Vec<(i32, u32)> = Vec::new();
                    for i in 0..header.alphabet {
                        for j in 0..header.alphabet {
                            want.push((
                                info.matrix_sorted[a][i] + info.matrix_sorted[b][j],
                                u32::from(info.matrix_sorted_char[a][i])
                                    + header.alphabet as u32
                                        * u32::from(info.matrix_sorted_char[b][j]),
                            ));
                        }
                    }
                    // The order of the stable sort by sum of the pairs in k order.
                    want.sort_by(|x, y| y.0.cmp(&x.0));
                    // Every group (the last makes the ones before it), then read as `emit`
                    // reads them.
                    row.make(X_PAIR_GROUPS - 1);
                    let mut got: Vec<(i32, u32)> = Vec::new();
                    for group in &row.groups {
                        let group = group.get().expect("made above");
                        for (&sum, &code) in group.sums.iter().zip(&group.codes) {
                            got.push((i32::from(sum), u32::from(code)));
                        }
                    }
                    assert_eq!(got[..used], want[..], "{matrix:?} pair ({a}, {b})");
                    assert!(got[used..]
                        .iter()
                        .all(|&(sum, code)| sum == i32::from(i16::MIN) && code == 0));
                    let mut keys = vec![
                        i32::MIN,
                        i32::MIN + 1,
                        -40_000,
                        -32_768,
                        -32_767,
                        32_767,
                        32_768,
                        40_000,
                        i32::MAX,
                    ];
                    for &(sum, _) in &want {
                        keys.extend([sum - 1, sum, sum + 1]);
                    }
                    for key in keys {
                        let scanned: Vec<u32> = want
                            .iter()
                            .filter(|&&(sum, _)| sum >= key)
                            .map(|&(_, code)| 7 + code * 3)
                            .collect();
                        let mut emitted = Vec::new();
                        row.emit(key, 7, 3, &mut emitted);
                        assert_eq!(emitted, scanned, "{matrix:?} ({a}, {b}) key {key}");
                    }
                }
            }
        }
    }
}

// LOSAT_X_CLUTPAIRS: the cells of an offset with the pair table are those of the loops of
// `add_word_hits` (`for_each_word`, and `for_each_word_in` for runs of last letters: the
// parallel builder's waves, `neighbors_below_last_paired`), each once, for word sizes 5 to 7,
// every matrix, thresholds from few neighbors to all words (i32::MIN: keys past i16), and words
// with letters outside the compressed alphabet.
#[test]
fn paired_cells_are_the_loops_cells() {
    let mut rng = XorShift(0xbb67_ae85_84ca_a73b);
    let mut seen = [0usize; 2];
    for (word_size, thresholds) in [
        (
            5usize,
            &[i32::MIN, -1, 1, 700, 1100, 1600, 2100, 3000, 7000][..],
        ),
        (6, &[1100, 1600, 2100, 3000][..]),
        (7, &[1200, 2025, 3000][..]),
    ] {
        let header = XHeader::new(word_size).expect("word size");
        let alphabet = header.alphabet as u32;
        for matrix in ScoringMatrix::ALL {
            let Some(scores) = matrix.score_matrix(header.alphabet) else {
                continue;
            };
            let info = CompressedNeighborInfo::new(&scores, header.alphabet);
            for &threshold in thresholds {
                let loops = XWords::with_pairs(&header, threshold, &scores, &info, false);
                let paired = XWords::with_pairs(&header, threshold, &scores, &info, true);
                assert!(loops.pairs.is_none() && paired.pairs.is_some());
                let rounds = if threshold <= 1 { 2 } else { 12 };
                for round in 0..rounds {
                    let letters: &[u8] = match round % 3 {
                        0 => &FEW_LETTERS,
                        1 => &LETTERS,
                        _ => &[
                            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 21, 24,
                        ],
                    };
                    let word: Vec<u8> = (0..word_size)
                        .map(|_| letters[rng.below(letters.len())])
                        .collect();
                    let mut want = Vec::new();
                    loops.for_each_word(&word[..], &mut want);
                    let mut got = Vec::new();
                    paired.for_each_word(&word[..], &mut got);
                    let mut distinct = got.clone();
                    distinct.sort_unstable();
                    distinct.dedup();
                    assert_eq!(distinct.len(), got.len(), "an offset appends a cell twice");
                    want.sort_unstable();
                    got.sort_unstable();
                    assert_eq!(
                        got, want,
                        "{matrix:?} word {word_size} threshold {threshold} word {word:?}"
                    );
                    seen[usize::from(!want.is_empty())] += 1;
                    // BLASTX's words of a length fixed by the instance (`[u8; W]`).
                    if word_size == 5 {
                        let fixed: [u8; 5] = word[..].try_into().expect("5 letters");
                        let mut fixed_got = Vec::new();
                        paired.for_each_word(&fixed, &mut fixed_got);
                        fixed_got.sort_unstable();
                        assert_eq!(fixed_got, want, "[u8; 5] word {word:?}");
                    }
                    let a = rng.below(alphabet as usize + 1) as u32;
                    let b = rng.below(alphabet as usize + 1) as u32;
                    // All words (threshold <= 1): the whole range and one letter's (keys
                    // past i16 in `neighbors_below_last_paired`).
                    let mut ranges = vec![0..alphabet, a..(a + 1).min(alphabet)];
                    if threshold > 1 {
                        ranges.extend([a..a, a.min(b)..a.max(b)]);
                        ranges.extend((0..alphabet).map(|letter| letter..letter + 1));
                    }
                    for range in ranges {
                        let mut want = Vec::new();
                        loops.for_each_word_in(&word[..], &range, &mut want);
                        let mut got = Vec::new();
                        paired.for_each_word_in(&word[..], &range, &mut got);
                        want.sort_unstable();
                        got.sort_unstable();
                        assert_eq!(
                            got, want,
                            "{matrix:?} word {word_size} threshold {threshold} word {word:?} letters {range:?}"
                        );
                    }
                }
            }
        }
    }
    assert!(seen.iter().all(|&n| n > 50), "cases seen {seen:?}");
}

fn flat_table(build: Option<XFlatBuild>) -> XFlatCompressedLookup {
    match build.expect("word size 5") {
        XFlatBuild::Table(table) => table,
        XFlatBuild::Refused(filled) => panic!("refused: {filled:?}"),
        XFlatBuild::TooLarge => panic!("too large"),
    }
}

// The reference's rounds of `scan_subject` over a subject's one range (`blastp_for_each_offset_pair`
// in two-hit mode), with an array of `array_size` pairs: the hits in order and the rounds.
fn rounds<F>(
    word_length: i32,
    subject_length: usize,
    array_size: usize,
    mut scan: F,
) -> (Vec<(u32, u32)>, usize)
where
    F: FnMut(&[(i32, i32)], &mut [OffsetPair], i32, &mut [i32; 3]) -> i32,
{
    let seq_ranges = [(0, subject_length as i32)];
    let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - word_length];
    if scan_range[2] < scan_range[1] {
        scan_range[2] = scan_range[1];
    }
    let mut offset_pairs = vec![OffsetPair::default(); array_size];
    let mut hits = Vec::new();
    let mut count = 0;
    while scan_range[1] <= scan_range[2] {
        let n = scan(
            &seq_ranges,
            &mut offset_pairs,
            array_size as i32,
            &mut scan_range,
        );
        hits.extend(
            offset_pairs[..n as usize]
                .iter()
                .map(|pair| (pair.s_off, pair.q_off)),
        );
        count += 1;
    }
    (hits, count)
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

// `x_to_scan_order` puts o_k at `x_emission_position(k, n)` for cells of 0 to 300 offsets, of
// both offset types of the tables (`i32`, `u16`).
#[test]
fn scan_order_in_place_is_the_emission_positions() {
    for n in 0..=300u32 {
        let mut expected = vec![u32::MAX; n as usize];
        for k in 0..n {
            expected[x_emission_position(k, n) as usize] = k;
        }
        let mut cell: Vec<i32> = (0..n as i32).collect();
        x_to_scan_order(&mut cell);
        assert!(
            cell.iter().map(|&v| v as u32).eq(expected.iter().copied()),
            "n {n}"
        );
        let mut narrow: Vec<u16> = (0..n as u16).collect();
        x_to_scan_order(&mut narrow);
        assert!(
            narrow
                .iter()
                .map(|&v| u32::from(v))
                .eq(expected.iter().copied()),
            "u16 n {n}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:825-841
// The permutation for cells of 1 to 40 offsets: the reference's scan order of a cell
// (`query_offsets_in_scan_order`, and the hits of the reference's `scan_subject` over a subject
// with every word once) puts the offset of rank k at `x_emission_position(k, n)`, which is a
// permutation of 0..n; the reference's overflow cells number cells(n).
#[test]
fn emission_positions_are_the_reference_scan_order_for_one_to_forty_offsets() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    let mut lookup = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
    let mut subject = Vec::new();
    let mut expected_hits = Vec::new();
    for n in 1..=40u32 {
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
        let cells_before = lookup.overflow_cells.len();
        let mut cells = 0usize;
        for k in 0..n {
            cells += usize::from(x_takes_overflow_cell(k));
            lookup.add_word_hit(index, (1000 * n + k) as i32);
        }
        assert_eq!(lookup.overflow_cells.len() - cells_before, cells, "n {n}");
        assert_eq!(cells, if n <= 5 { 0 } else { n as usize / 4 }, "n {n}");
        let mut expected = vec![-1i32; n as usize];
        for k in 0..n {
            let position = x_emission_position(k, n) as usize;
            assert_eq!(expected[position], -1, "n {n}: position {position} twice");
            expected[position] = (1000 * n + k) as i32;
        }
        assert_eq!(lookup.query_offsets_in_scan_order(index), expected, "n {n}");
        let s = subject.len() as u32;
        expected_hits.extend(expected.iter().map(|&q| (s, q as u32)));
        subject.extend_from_slice(&word);
        subject.push(OUTSIDE[rng.below(OUTSIDE.len())]);
    }
    lookup.finalize();
    let length = subject.len();
    subject.push(0);
    for array_size in [40, 41, 77, 4096 + 40] {
        let (hits, _) = rounds(5, length, array_size, |ranges, pairs, size, range| {
            lookup.scan_subject(&subject, ranges, pairs, size, range)
        });
        assert_eq!(hits, expected_hits, "array {array_size}");
    }
}

// Queries of `lengths` letters from `letters`, concatenated as BLASTP's lookup query (each query
// followed by a sentinel). Each query has one lookup segment, or `masked` gaps (SEG-masked runs)
// that cut it into several.
pub(super) fn queries(
    rng: &mut XorShift,
    lengths: &[usize],
    letters: &[u8],
    masked: usize,
) -> (Vec<u8>, Vec<(i32, i32)>, Vec<usize>) {
    let mut query = Vec::new();
    let mut locations = Vec::new();
    // The first segment of each query.
    let mut first = Vec::new();
    for &length in lengths {
        let left = query.len() as i32;
        for _ in 0..length {
            query.push(letters[rng.below(letters.len())]);
        }
        first.push(locations.len());
        let mut from = left;
        let end = left + length as i32 - 1;
        for _ in 0..masked {
            if end - from < 30 {
                break;
            }
            let gap = from + 5 + rng.below((end - from - 20) as usize) as i32;
            let gap_end = (gap + 1 + rng.below(12) as i32).min(end);
            if gap > from {
                locations.push((from, gap - 1));
            }
            from = gap_end + 1;
        }
        if from <= end {
            locations.push((from, end));
        }
        query.push(0);
    }
    first.push(locations.len());
    (query, locations, first)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
// Random queries (few letters and low thresholds, so that cells overflow; all letters with letters
// outside the compressed alphabet), SEG-masked runs (several segments per query), several NCBI
// tables and a split batch (segments outside every table), threshold 0 (the words themselves
// only): the flat table is the reference's, cell by cell in scan order, with the same PV array and
// header.
#[test]
fn flat_table_equals_reference_on_random_queries() {
    let mut rng = XorShift(0xbb67_ae85_84ca_a73b);
    let mut overflow_seen = false;
    for round in 0..10 {
        let (letters, threshold): (&[u8], i32) = match round % 5 {
            0 => (&FEW_LETTERS[..], 1100),
            1 => (&FEW_LETTERS[..], 1600),
            2 => (&LETTERS[..], 1100),
            3 => (&[1u8, 11, 17, 21, 24, 10, 3][..], 2000),
            _ => (&FEW_LETTERS[..], 0),
        };
        let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 1 + rng.below(300)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
        // NCBI's tables: consecutive queries' segments; some batches are left out (split).
        let mut tables = Vec::new();
        let mut q = 0;
        while q < lengths.len() {
            let end = (q + 1 + rng.below(2)).min(lengths.len());
            if rng.below(3) != 0 {
                tables.push(first[q]..first[end]);
            }
            q = end;
        }
        let reference =
            try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                .expect("word size 5")
                .expect("not refused");
        let flat = flat_table(x_try_build_flat_compressed_lookup(
            5, threshold, &query, &locations, &tables,
        ));
        x_assert_same_compressed_table(&flat, &reference);
        overflow_seen |= reference.longest_chain > 5;
    }
    assert!(overflow_seen, "cells with overflow lists");
    // No words at all (every segment shorter than a word).
    let (query, locations, _) = queries(&mut rng, &[3, 4, 2], &FEW_LETTERS, 0);
    let reference = try_build_blosum62_compressed_lookup(5, 1100, &query, &locations, &[])
        .expect("word size 5")
        .expect("not refused");
    let flat = flat_table(x_try_build_flat_compressed_lookup(
        5,
        1100,
        &query,
        &locations,
        &[],
    ));
    x_assert_same_compressed_table(&flat, &reference);
    assert!(flat.hits.is_empty() && flat.longest_chain == 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310
// Every search matrix but IDENTITY at word size 5 (BLOSUM45 also at 6, PAM30 at 7), random queries
// with NCBI's tables and a split batch, the word sizes' thresholds and lower ones: the bucket fill (word size 5) and the
// count-and-fill builder of the matrix's words build the reference table of those words
// (`try_build_compressed_lookup`), or refuse where it refuses, and
// `x_tbn_build_compressed_lookup_of` takes them.
#[test]
fn flat_tables_of_each_matrix_are_its_reference_tables() {
    let mut rng = XorShift(0x3c6e_f372_fe94_f82b);
    let (mut built_seen, mut refused_seen) = (0, 0);
    for matrix in ScoringMatrix::ALL {
        if matrix == ScoringMatrix::Identity {
            assert!(x_try_build_flat_of(matrix, 5, 1100, &[], &[], &[], 1, 1).is_none());
            continue;
        }
        // Word sizes 6 and 7 (backbones of 15^6 and 10^7 cells) with two matrices only.
        let larger: &[(usize, &[i32], usize)] = match matrix {
            ScoringMatrix::Blosum45 => &[(6, &[2100, 1300], NCBI_OVERFLOW_CELL_CAPACITY)],
            ScoringMatrix::Pam30 => &[(7, &[2025, 1200], NCBI_OVERFLOW_CELL_CAPACITY)],
            _ => &[],
        };
        for &(word, thresholds, capacity) in [
            (5usize, &[1930, 1100][..], NCBI_OVERFLOW_CELL_CAPACITY),
            (5, &[1600, 0][..], 4),
        ]
        .iter()
        .chain(larger)
        {
            for &threshold in thresholds {
                let lengths: Vec<usize> =
                    (0..2 + rng.below(4)).map(|_| 8 + rng.below(250)).collect();
                let letters: &[u8] = if rng.below(2) == 0 {
                    &FEW_LETTERS
                } else {
                    &LETTERS
                };
                let (query, locations, first) = queries(&mut rng, &lengths, letters, 1);
                let tables = [first[0]..first[lengths.len() - 1]];
                let what = format!("{matrix:?} word {word} threshold {threshold}");
                let reference = try_build_compressed_lookup_with_capacity(
                    matrix, word, threshold, &query, &locations, &tables, capacity,
                )
                .expect("word size");
                let counted = x_try_build_flat_of(
                    matrix,
                    word,
                    threshold,
                    &query,
                    &locations,
                    &tables,
                    capacity,
                    X_COUNT_CHUNK,
                )
                .expect("word size");
                match (&reference, counted) {
                    (Err(want), XFlatBuild::Refused(got)) if *want == got => refused_seen += 1,
                    (Ok(reference), XFlatBuild::Table(flat)) => {
                        x_assert_same_compressed_table(&flat, reference);
                        built_seen += 1;
                    }
                    (reference, _) => {
                        panic!("{what}: reference refused {}", reference.is_err())
                    }
                }
                match x_buckets_of(
                    matrix, word, threshold, &query, &locations, &tables, capacity,
                ) {
                    // Word sizes 6 and 7 and small capacities only (`x_try_build_flat_by_buckets`).
                    XBucketBuild::Declined => assert!(
                        word != 5 || capacity != NCBI_OVERFLOW_CELL_CAPACITY,
                        "{what}: declined"
                    ),
                    XBucketBuild::Built(Some(XFlatBuild::Table(flat))) => {
                        let reference = reference.as_ref().expect(&what);
                        x_assert_same_compressed_table(&flat, reference);
                    }
                    XBucketBuild::Built(_) => panic!("{what}: buckets refused"),
                }
                if let Some(reference) = reference
                    .as_ref()
                    .ok()
                    .filter(|_| capacity == NCBI_OVERFLOW_CELL_CAPACITY)
                {
                    let Some(Ok(XCompressedLookup::Flat(flat))) = x_tbn_build_compressed_lookup_of(
                        matrix, word, threshold, &query, &locations, &tables, 1,
                    ) else {
                        panic!("{what}: no flat table")
                    };
                    x_assert_same_compressed_table(&flat, reference);
                }
            }
        }
    }
    assert!(
        built_seen > 0 && refused_seen > 0,
        "{built_seen} {refused_seen}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
// The refusal at small capacities: on random queries and tables (with split batches), at
// capacities around each table's own cells and its cells of the table of all batches, the flat
// builder refuses exactly when the reference does, and otherwise builds the reference's table;
// `x_table_fills_banks` is `compressed_table_fills_banks`.
#[test]
fn flat_refusal_is_the_references_at_small_capacities() {
    let mut rng = XorShift(0x1234_5678_9abc_def1);
    let mut refused_seen = 0usize;
    let mut built_seen = 0usize;
    let header = XHeader::new(5).expect("word size 5");
    let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
    let info = CompressedNeighborInfo::new(&matrix, 15);
    for round in 0..8 {
        let threshold = [1100, 1600, 2000][round % 3];
        let words = XWords::new(&header, threshold, &matrix, &info);
        let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 5 + rng.below(350)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, &FEW_LETTERS, round % 2);
        let mut tables = Vec::new();
        let mut q = 0;
        while q < lengths.len() {
            let end = (q + 1 + rng.below(3)).min(lengths.len());
            if rng.below(4) != 0 {
                tables.push(first[q]..first[end]);
            }
            q = end;
        }
        // Each table's own cells (the reference's table of its segments).
        let mut capacities = vec![1usize];
        for table in &tables {
            let mut own = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
            own.threshold = threshold;
            own.overflow_cell_capacity = usize::MAX;
            own.add_neighboring_words(&matrix, &query, &locations[table.clone()]);
            let cells = own.overflow_cells.len();
            capacities.extend([cells.saturating_sub(1), cells, cells + 1, 2 * cells]);
            for capacity in [cells.saturating_sub(1), cells]
                .into_iter()
                .filter(|&c| c > 0)
            {
                assert_eq!(
                    x_table_fills_banks(
                        &words,
                        header.backbone_size,
                        &query,
                        &locations[table.clone()],
                        capacity
                    ),
                    compressed_table_fills_banks(
                        5,
                        threshold,
                        &matrix,
                        &query,
                        &locations[table.clone()],
                        capacity
                    ),
                    "round {round}, capacity {capacity}"
                );
            }
        }
        capacities.retain(|&capacity| capacity > 0);
        for capacity in capacities {
            let reference = try_build_compressed_lookup_with_capacity(
                Blosum62Compressed,
                5,
                threshold,
                &query,
                &locations,
                &tables,
                capacity,
            )
            .expect("word size 5");
            // The bucket fill (LOSAT_X_BPFFILL) builds the same table or declines.
            match x_try_build_flat_by_buckets(
                5, threshold, &query, &locations, &tables, capacity, 0,
            ) {
                XBucketBuild::Declined => {}
                XBucketBuild::Built(Some(XFlatBuild::Table(flat))) => {
                    let reference = reference.as_ref().unwrap_or_else(|_| {
                        panic!(
                            "round {round}, capacity {capacity}: buckets built, reference refused"
                        )
                    });
                    x_assert_same_compressed_table(&flat, reference);
                }
                XBucketBuild::Built(_) => {
                    panic!("round {round}, capacity {capacity}: buckets refused")
                }
            }
            // Chunks of one append (a test after each query offset with appends, as the
            // reference tests), of a few, and the default (a test at each table's end).
            for chunk in [1, 7, X_COUNT_CHUNK] {
                let built = x_try_build_flat_with_capacity(
                    5, threshold, &query, &locations, &tables, capacity, chunk, 0,
                )
                .expect("word size 5");
                match (&reference, built) {
                    (Err(want), XFlatBuild::Refused(got)) if *want == got => refused_seen += 1,
                    (Ok(reference), XFlatBuild::Table(flat)) => {
                        x_assert_same_compressed_table(&flat, reference);
                        built_seen += 1;
                    }
                    (reference, _) => panic!(
                        "round {round}, capacity {capacity}, chunk {chunk}: reference refused {}",
                        reference.is_err()
                    ),
                }
            }
        }
    }
    assert!(refused_seen > 0 && built_seen > 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-849 (module header)
// The bucket fill (LOSAT_X_BPFFILL): on random queries (few letters and low thresholds, so that
// cells overflow; letters outside the compressed alphabet), SEG-masked runs, NCBI's tables with
// split batches, threshold 0, and capacities from 1 to NCBI's, it builds the count-and-fill
// builder's table field by field and slot by slot, or declines where a table's appends pass the
// capacity (never at NCBI's); a query of 2^22 letters or more does not fit its packed appends.
#[test]
fn x_bucket_fill_is_the_count_and_fill_build() {
    let mut rng = XorShift(0x5be0_cd19_137e_2179);
    let (mut built, mut declined, mut overflow_seen) = (0usize, 0usize, false);
    for round in 0..24 {
        let (letters, threshold): (&[u8], i32) = match round % 5 {
            0 => (&FEW_LETTERS[..], 1100),
            1 => (&FEW_LETTERS[..], 1600),
            2 => (&LETTERS[..], 1100),
            3 => (&[1u8, 11, 17, 21, 24, 10, 3][..], 2000),
            _ => (&FEW_LETTERS[..], 0),
        };
        let lengths: Vec<usize> = (0..2 + rng.below(5)).map(|_| 1 + rng.below(400)).collect();
        let (query, locations, first) = queries(&mut rng, &lengths, letters, round % 3);
        let mut tables = Vec::new();
        let mut q = 0;
        while q < lengths.len() {
            let end = (q + 1 + rng.below(2)).min(lengths.len());
            if rng.below(3) != 0 {
                tables.push(first[q]..first[end]);
            }
            q = end;
        }
        for capacity in [1usize, 50, 400, NCBI_OVERFLOW_CELL_CAPACITY] {
            let want = x_try_build_flat_with_capacity(
                5,
                threshold,
                &query,
                &locations,
                &tables,
                capacity,
                X_COUNT_CHUNK,
                0,
            );
            let got =
                x_try_build_flat_by_buckets(5, threshold, &query, &locations, &tables, capacity, 0);
            match (got, want) {
                (XBucketBuild::Declined, _) => {
                    assert!(
                        capacity < NCBI_OVERFLOW_CELL_CAPACITY,
                        "round {round}: declined"
                    );
                    declined += 1;
                }
                (
                    XBucketBuild::Built(Some(XFlatBuild::Table(got))),
                    Some(XFlatBuild::Table(want)),
                ) => {
                    x_assert_same_flat_tables(&got, &want);
                    overflow_seen |= want.longest_chain > 5;
                    built += 1;
                }
                (XBucketBuild::Built(_), _) => {
                    panic!("round {round}, capacity {capacity}: another outcome")
                }
            }
        }
    }
    assert!(built > 0 && declined > 0 && overflow_seen);
    // Offsets up to 2^22 - 1 (bit 31 of the packed append set): a query of 2^22 - 1 letters
    // with one segment at its end; a query of 2^22 letters is declined, and so is word size 6
    // (the PV shift depends on the counts).
    let length = (1usize << 22) - 1;
    let mut query: Vec<u8> = (0..length)
        .map(|_| FEW_LETTERS[rng.below(FEW_LETTERS.len())])
        .collect();
    let locations = [(length as i32 - 400, length as i32 - 1)];
    let capacity = NCBI_OVERFLOW_CELL_CAPACITY;
    let want = x_try_build_flat_with_capacity(
        5,
        1100,
        &query,
        &locations,
        &[0..1],
        capacity,
        X_COUNT_CHUNK,
        0,
    );
    let got = x_try_build_flat_by_buckets(5, 1100, &query, &locations, &[0..1], capacity, 0);
    let (XBucketBuild::Built(Some(XFlatBuild::Table(got))), Some(XFlatBuild::Table(want))) =
        (got, want)
    else {
        panic!("2^22 - 1 letters: not two tables");
    };
    assert!(got.hits.iter().any(|&offset| offset >= 1 << 21));
    x_assert_same_flat_tables(&got, &want);
    query.push(FEW_LETTERS[0]);
    let got = x_try_build_flat_by_buckets(5, 1100, &query, &locations, &[0..1], capacity, 0);
    assert!(matches!(got, XBucketBuild::Declined), "2^22 letters");
    let got =
        x_try_build_flat_by_buckets(6, 1100, &query[..1000], &[(0, 999)], &[0..1], capacity, 0);
    assert!(matches!(got, XBucketBuild::Declined), "word size 6");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:340-404
// The scans: on random subjects (letters outside the alphabet none, rare, frequent; short
// subjects), the flat `scan_subject` gives the reference's hits and rounds for arrays from the
// longest cell up, and `for_each_cell` gives the same hits, with the branch-free list too
// (LOSAT_X_BPCLUTSCANBF).
#[test]
fn flat_scans_give_the_reference_scans_hits() {
    let mut rng = XorShift(0x3c6e_f372_fe94_f82b);
    for (threshold, letters) in [
        (1100, &FEW_LETTERS[..]),
        (900, &[1u8, 11, 17][..]),
        (1300, &LETTERS[..]),
    ] {
        let lengths = [120usize, 75, 200, 33];
        let (query, locations, _) = queries(&mut rng, &lengths, letters, 1);
        let tables = [0..locations.len()];
        let reference =
            try_build_blosum62_compressed_lookup(5, threshold, &query, &locations, &tables)
                .expect("word size 5")
                .expect("not refused");
        let flat = flat_table(x_try_build_flat_compressed_lookup(
            5, threshold, &query, &locations, &tables,
        ));
        assert!(reference.longest_chain > 5, "cells with overflow lists");
        let used = flat.starts.windows(2).filter(|w| w[1] > w[0]).count();
        assert!(used > 0, "cells with hits");
        let longest = reference.longest_chain as usize;
        let mut subjects = Vec::new();
        for outside_one_in in [0usize, 50, 6, 2] {
            for _ in 0..5 {
                let length = 1 + rng.below(300);
                let mut subject = random_subject(&mut rng, length, outside_one_in);
                if length > 60 && rng.below(2) == 0 {
                    let from = rng.below(query.len() - 60);
                    let to = rng.below(length - 50);
                    subject[to..to + 50].copy_from_slice(&query[from..from + 50]);
                }
                subjects.push(subject);
            }
        }
        for length in 1..=6 {
            subjects.push(random_subject(&mut rng, length, 0));
        }
        let mut cells = Vec::new();
        for (n, subject) in subjects.iter().enumerate() {
            let length = subject.len() - 1;
            let mut expected = None;
            for array_size in [longest, longest + 1, 2 * longest + 3, 4096 + longest] {
                let (want, want_rounds) =
                    rounds(5, length, array_size, |ranges, pairs, size, range| {
                        reference.scan_subject(subject, ranges, pairs, size, range)
                    });
                let (got, got_rounds) =
                    rounds(5, length, array_size, |ranges, pairs, size, range| {
                        flat.scan_subject::<X_BLASTP>(subject, ranges, pairs, size, range)
                    });
                assert_eq!(got, want, "subject {n}, array {array_size}");
                assert_eq!(got_rounds, want_rounds, "subject {n}, array {array_size}");
                expected.get_or_insert(want);
            }
            let mut fused = Vec::new();
            flat.for_each_cell::<false, _>(subject, length, &mut cells, |s, hits| {
                fused.extend(hits.iter().map(|&q| (s as u32, q as u32)));
            });
            assert_eq!(Some(fused), expected, "subject {n}: for_each_cell");
            // LOSAT_X_BPCLUTSCANBF: the branch-free list, on the used scratch and on a new one.
            for mut scratch in [std::mem::take(&mut cells), Vec::new()] {
                let mut fused = Vec::new();
                flat.for_each_cell::<true, _>(subject, length, &mut scratch, |s, hits| {
                    fused.extend(hits.iter().map(|&q| (s as u32, q as u32)));
                });
                assert_eq!(
                    Some(fused),
                    expected,
                    "subject {n}: branch-free for_each_cell"
                );
                cells = scratch;
            }
            // The lean loop's list (`x_list_pairs_lean`, two offsets per step) on a subject of
            // a word or more with its sentinel; a shorter one, or one without the sentinel, is
            // declined.
            assert_eq!(
                flat.x_pairs_fit(subject, length),
                length >= 5,
                "subject {n}: fits"
            );
            assert!(
                !flat.x_pairs_fit(&subject[..length], length),
                "subject {n}: no sentinel"
            );
            if length < 5 {
                continue;
            }
            let want_list = cells.clone();
            let mut lean = Vec::new();
            // SAFETY: x_pairs_fit holds (asserted above).
            unsafe { flat.x_list_pairs_lean::<X_BLASTP>(subject, length, &mut lean) };
            assert_eq!(lean, want_list, "subject {n}: x_list_pairs_lean's list");
            assert!(
                lean.capacity() >= lean.len() + X_STARTS_AHEAD,
                "subject {n}: lean reserve"
            );
            let listed: Vec<(u32, u32)> = lean
                .iter()
                .flat_map(|&(s, index)| {
                    flat.cell(index as usize)
                        .iter()
                        .map(move |&q| (s as u32, q as u32))
                })
                .collect();
            assert_eq!(Some(listed), expected, "subject {n}: x_list_pairs_lean");
        }
    }
}

// LOSAT_X_BPCLUTSCANBF's lean loop: `x * 19,088,744 >> 32` is `x / 225` below 15^6, NCBI's
// reciprocal multiply is `x / 15` there (the indices of word lengths up to 6), and the multiply
// by the inverse of 15 is `x / 15` on its multiples.
#[test]
fn the_pairs_scan_divisions_are_exact_below_fifteen_to_the_sixth() {
    for x in 0..15i64.pow(6) {
        assert_eq!((x * 19_088_744) >> 32, x / 225, "x {x}");
        assert_eq!(
            i64::from(compressed_reciprocal_divide(x as i32, 286_331_154)),
            x / 15,
            "x {x}"
        );
        if x % 15 == 0 {
            assert_eq!(
                i64::from((x as u32).wrapping_mul(0xEEEE_EEEF)),
                x / 15,
                "x {x}"
            );
        }
    }
}

// LOSAT_X_BPCLUTSCANBF's lean loop: the list two offsets per step on word length 6 too, and
// the tables and subjects it declines.
#[test]
fn the_pairs_scan_lists_word_length_six_and_declines_outside_its_range() {
    let mut rng = XorShift(0x6a09_e667_f3bc_c909);
    let (query, locations, _) = queries(&mut rng, &[200usize, 150], &LETTERS[..], 1);
    let tables = [0..locations.len()];
    let mut flat6 = flat_table(x_try_build_flat_compressed_lookup(
        6, 1100, &query, &locations, &tables,
    ));
    // A table with under 1 % of its cells used takes a coarser shift of its bits (NCBI's
    // `pv_array_bts`), which the scan declines; this one has one bit per cell.
    let probe = random_subject(&mut rng, 60, 0);
    assert_eq!(
        flat6.pv_array_bts, PV_ARRAY_BTS as i32,
        "a denser word-6 table"
    );
    assert!(flat6.x_pairs_fit(&probe, 60));
    flat6.pv_array_bts += 2;
    assert!(!flat6.x_pairs_fit(&probe, 60), "coarse bits");
    flat6.pv_array_bts -= 2;
    let mut cells = Vec::new();
    let mut found = 0;
    for outside_one_in in [0usize, 6, 2] {
        for length in [1usize, 5, 6, 7, 8, 60, 301] {
            let subject = random_subject(&mut rng, length, outside_one_in);
            assert_eq!(flat6.x_pairs_fit(&subject, length), length >= 6);
            if length < 6 {
                continue;
            }
            let mut want = Vec::new();
            flat6.for_each_cell::<true, _>(&subject, length, &mut cells, |s, hits| {
                want.extend(hits.iter().map(|&q| (s, q)));
            });
            let want_list = cells.clone();
            let mut lean = Vec::new();
            // SAFETY: x_pairs_fit holds (asserted above).
            unsafe { flat6.x_list_pairs_lean::<X_BLASTP>(&subject, length, &mut lean) };
            assert_eq!(lean, want_list, "length {length}: list");
            found += want.len();
            let got: Vec<(i32, i32)> = lean
                .iter()
                .flat_map(|&(s, index)| flat6.cell(index as usize).iter().map(move |&q| (s, q)))
                .collect();
            assert_eq!(got, want, "length {length}");
        }
    }
    assert!(found > 0, "word length 6 finds hits");
    let subject = random_subject(&mut rng, 50, 0);
    assert!(
        !flat6.x_pairs_fit(&subject, 51),
        "no byte past the last word"
    );
    assert!(!flat6.x_pairs_fit(&subject, 52), "length past the buffer");
    flat6.reciprocal_alphabet_size += 1;
    assert!(!flat6.x_pairs_fit(&subject, 50), "another reciprocal");
    flat6.reciprocal_alphabet_size -= 1;
    flat6.pv.pop();
    assert!(!flat6.x_pairs_fit(&subject, 50), "bits short of 15^6");
}

// BLASTX's bucket fill (`u16` offsets, LOSAT_X_BXCLUTPOOL on kept vectors or not) makes
// `x_buckets_on`'s table with its offsets narrowed: queries of one segment and of several in
// NCBI's tables, low-complexity ones whose cells pass 32 offsets (the scan order), and one of
// `X_U16_LETTERS` letters.
#[test]
#[cfg_attr(
    target_pointer_width = "32",
    ignore = "builds a table of over 268 million appends (threshold 1100), the most a 32-bit target takes (`x_max_appends`: the record's doubling stays below isize::MAX bytes); there the product takes the TooLarge fallback, so this table is not built"
)]
fn pooled_flat_tables_equal_the_unpooled_ones_on_kept_vectors() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    let segments = |length: usize, pieces: usize| -> Vec<(i32, i32)> {
        let step = length / pieces;
        (0..pieces)
            .map(|p| {
                let from = p * step;
                let to = if p + 1 == pieces {
                    length
                } else {
                    from + step - 3
                };
                (from as i32, to as i32 - 1)
            })
            .collect()
    };
    let wide = |query: &[u8], locations: &[(i32, i32)], tables: &[std::ops::Range<usize>]| {
        match x_buckets_on(5, 1100, query, locations, tables, usize::MAX, false, 0) {
            XBucketBuild::Built(Some(XFlatBuild::Table(table))) => x_narrow_flat(table),
            _ => panic!("the bucket fill builds this table"),
        }
    };
    let narrow = |query: &[u8],
                  locations: &[(i32, i32)],
                  tables: &[std::ops::Range<usize>],
                  pooled: bool| {
        match x_bx_buckets_on::<5>(1100, query, locations, tables, usize::MAX, pooled, 1) {
            XBucketBuild::Built(Some(XFlatBuild::Table(table))) => table,
            _ => panic!("BLASTX's bucket fill builds this table"),
        }
    };
    let same = |got: &XFlatCompressedLookup<u16>, want: &XFlatCompressedLookup<u16>, what: &str| {
        assert_eq!(got.hits, want.hits, "{what}");
        assert_eq!(got.starts, want.starts, "{what}");
        assert_eq!(got.pv, want.pv, "{what}");
        assert_eq!(got.longest_chain, want.longest_chain, "{what}");
        assert_eq!(got.pv_array_bts, want.pv_array_bts, "{what}");
        assert_eq!(got.word_length, want.word_length, "{what}");
    };
    let mut alive = None;
    let mut longest = 0;
    for (round, (length, pieces, repeat)) in [
        (3_000, 1, 0),
        (40, 1, 0),
        (9_000, 3, 0),
        (200, 1, 0),
        (20_000, 6, 0),
        (7, 1, 0),
        (1_500, 2, 40),
        (20_000, 1, 0),
        (6_000, 4, 200),
        // `random_subject` adds a letter.
        (X_U16_LETTERS - 1, 6, 0),
    ]
    .into_iter()
    .enumerate()
    {
        let mut query = random_subject(&mut rng, length, 0);
        // Low complexity: the first `repeat` letters repeat one word.
        for (i, letter) in query.iter_mut().take(repeat).enumerate() {
            *letter = [1, 5, 9, 1, 5][i % 5];
        }
        let locations = segments(length, pieces);
        let tables: Vec<std::ops::Range<usize>> = if pieces > 2 {
            vec![0..1, 2..pieces]
        } else {
            vec![0..pieces]
        };
        let want = wide(&query, &locations, &tables);
        longest = longest.max(want.longest_chain);
        let unpooled = narrow(&query, &locations, &tables, false);
        same(&unpooled, &want, &format!("round {round} unpooled"));
        let got = narrow(&query, &locations, &tables, true);
        same(&got, &want, &format!("round {round} pooled"));
        match round % 3 {
            0 => {
                if let Some(mut kept) = alive.replace(got) {
                    x_give_back_flat(&mut kept);
                }
            }
            // Values no table holds: the build must write every entry it keeps (starts, hits)
            // and clear the PV array.
            _ => {
                let mut poisoned = got;
                poisoned.starts.fill(u32::MAX);
                poisoned.hits.fill(u16::MAX);
                poisoned.pv.fill(u32::MAX);
                x_give_back_flat(&mut poisoned);
            }
        }
    }
    if let Some(mut kept) = alive {
        x_give_back_flat(&mut kept);
    }
    assert!(longest > 32, "a cell passes 32 offsets ({longest})");
}

// BLASTX's table: `u16` offsets up to `X_U16_LETTERS` letters, `i32` ones past them;
// the count-and-fill builder's table (where the bucket fill declines) narrowed is the bucket
// fill's, here on locations out of query order in two tables.
#[test]
#[cfg_attr(
    target_pointer_width = "32",
    ignore = "builds a table of over 268 million appends (threshold 1100), the most a 32-bit target takes (`x_max_appends`: the record's doubling stays below isize::MAX bytes); there the product takes the TooLarge fallback, so this table is not built"
)]
fn blastx_tables_take_u16_offsets_up_to_their_letters() {
    let mut rng = XorShift(0x2545_F491_4F6C_DD1D);
    let query = random_subject(&mut rng, X_U16_LETTERS + 1, 0);
    let all = [(0, query.len() as i32 - 1)];
    assert!(matches!(
        x_bx_build_compressed_lookup::<5>(1100, &query, &all, &[0..1], 1),
        Some(Ok(XBxCompressedLookup::Wide(_)))
    ));
    let query = &query[..X_U16_LETTERS];
    assert!(matches!(
        x_bx_build_compressed_lookup::<5>(1100, query, &[(0, query.len() as i32 - 1)], &[0..1], 1),
        Some(Ok(XBxCompressedLookup::Flat(_)))
    ));
    let swapped = [(5_000, 9_999), (0, 4_999)];
    let tables = [0..1, 1..2];
    let XBucketBuild::Built(Some(XFlatBuild::Table(got))) =
        x_bx_buckets_on::<5>(1100, query, &swapped, &tables, usize::MAX, false, 1)
    else {
        panic!("BLASTX's bucket fill builds this table");
    };
    let Some(XFlatBuild::Table(want)) = x_try_build_flat_with_capacity(
        5,
        1100,
        query,
        &swapped,
        &tables,
        NCBI_OVERFLOW_CELL_CAPACITY,
        X_COUNT_CHUNK,
        0,
    ) else {
        panic!("the count-and-fill builder builds this table");
    };
    let want = x_narrow_flat(want);
    assert_eq!(got.starts, want.starts);
    assert_eq!(got.hits, want.hits);
    assert_eq!(got.pv, want.pv);
    assert_eq!(got.longest_chain, want.longest_chain);
}

// Line NP, SC9e: BLASTX's bucket fill in parts (`x_bx_buckets_with`) makes the one-part fill's
// table, slot for slot, on a pool of 4 threads: BLOSUM62's (`u16`) and BLOSUM45's (`i32`, which
// is also the record fill's, `x_buckets_of`), on segments in one table and in two with a gap,
// low-complexity ones whose cells pass 32 offsets (the scan order), runs cut inside segments,
// more parts than offsets, and an empty query; a table past the overflow capacity is declined in
// parts as at one part. The runs (`x_offset_parts`) hold every offset once, in order.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn blastx_bucket_fill_in_parts_is_the_one_part_fill() {
    fn offsets(word: usize, runs: &[Vec<(usize, (i32, i32))>]) -> Vec<(usize, i32)> {
        let mut all = Vec::new();
        for run in runs {
            for &(index, (left, right)) in run {
                all.extend((left..=right - word as i32 + 1).map(|offset| (index, offset)));
            }
        }
        all
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("test pool");
    let mut rng = XorShift(0x2545_F491_4F6C_DD1D);
    let mut longest = 0;
    for (round, (length, segments, tables, repeat)) in [
        (
            4_400,
            vec![(3, 1_499), (1_501, 2_998), (3_001, 4_399)],
            vec![0..3],
            0,
        ),
        (
            9_000,
            vec![(0, 2_999), (3_002, 5_999), (6_004, 8_999)],
            vec![0..1, 2..3],
            0,
        ),
        (1_500, vec![(0, 1_499)], vec![0..1], 300),
        (12, vec![(0, 11)], vec![0..1], 0),
        (0, vec![], vec![0..0], 0),
    ]
    .into_iter()
    .enumerate()
    {
        let mut query = random_subject(&mut rng, length, 0);
        for (i, letter) in query.iter_mut().take(repeat).enumerate() {
            *letter = [1, 5, 9, 1, 5][i % 5];
        }
        let ranges = x_count_ranges(&tables, segments.len());
        let serial = offsets(5, &x_offset_parts(5, &ranges, &segments, 1));
        for parts in [2, 3, 7, 64] {
            let runs = x_offset_parts(5, &ranges, &segments, parts);
            assert!(runs.len() <= parts, "round {round}: at most {parts} runs");
            assert_eq!(offsets(5, &runs), serial, "round {round}: {parts} runs");
        }
        let narrow = |parts: usize, capacity: usize| {
            x_bx_buckets_on::<5>(1100, &query, &segments, &tables, capacity, false, parts)
        };
        let wide = |parts: usize, capacity: usize| {
            x_bx_buckets_with::<_, i32, [u8]>(
                ScoringMatrix::Blosum45,
                5,
                1100,
                &query,
                &segments,
                &tables,
                capacity,
                false,
                parts,
                |word| word,
            )
        };
        let XBucketBuild::Built(Some(XFlatBuild::Table(want))) = narrow(1, usize::MAX) else {
            panic!("round {round}: BLASTX's bucket fill builds this table");
        };
        let XBucketBuild::Built(Some(XFlatBuild::Table(want_wide))) = wide(1, usize::MAX) else {
            panic!("round {round}: BLASTX's bucket fill builds BLOSUM45's table");
        };
        let XBucketBuild::Built(Some(XFlatBuild::Table(record))) = x_buckets_of(
            ScoringMatrix::Blosum45,
            5,
            1100,
            &query,
            &segments,
            &tables,
            usize::MAX,
        ) else {
            panic!("round {round}: the record fill builds BLOSUM45's table");
        };
        assert_eq!(want_wide.starts, record.starts, "round {round}: B45 starts");
        assert_eq!(want_wide.hits, record.hits, "round {round}: B45 hits");
        assert_eq!(want_wide.pv, record.pv, "round {round}: B45 pv");
        assert_eq!(
            want_wide.longest_chain, record.longest_chain,
            "round {round}"
        );
        assert_eq!(want_wide.pv_array_bts, record.pv_array_bts, "round {round}");
        longest = longest.max(want.longest_chain);
        for parts in [2, 3, 5, 64] {
            let what = format!("round {round} in {parts} parts");
            let XBucketBuild::Built(Some(XFlatBuild::Table(got))) =
                pool.install(|| narrow(parts, usize::MAX))
            else {
                panic!("{what}: no table");
            };
            assert_eq!(got.starts, want.starts, "{what}: starts");
            assert_eq!(got.hits, want.hits, "{what}: hits");
            assert_eq!(got.pv, want.pv, "{what}: pv");
            assert_eq!(got.longest_chain, want.longest_chain, "{what}");
            assert_eq!(got.pv_array_bts, want.pv_array_bts, "{what}");
            let XBucketBuild::Built(Some(XFlatBuild::Table(got))) =
                pool.install(|| wide(parts, usize::MAX))
            else {
                panic!("{what}: no BLOSUM45 table");
            };
            assert_eq!(got.starts, want_wide.starts, "{what}: B45 starts");
            assert_eq!(got.hits, want_wide.hits, "{what}: B45 hits");
            assert_eq!(got.pv, want_wide.pv, "{what}: B45 pv");
            assert_eq!(got.longest_chain, want_wide.longest_chain, "{what}");
        }
        // Past the overflow capacity: declined in parts as at one part.
        let total: u32 = *want.starts.last().expect("the number of hits");
        if total > 1 {
            let capacity = total as usize / 4;
            assert!(matches!(narrow(1, capacity), XBucketBuild::Declined));
            for parts in [2, 3] {
                assert!(
                    matches!(
                        pool.install(|| narrow(parts, capacity)),
                        XBucketBuild::Declined
                    ),
                    "round {round}: declined in {parts} parts"
                );
            }
        }
    }
    assert!(longest > 32, "a cell passes 32 offsets ({longest})");
}

// Item SC9b (2) of line NP: BLASTX's table of a unit built in parts on a pool (`x_flat_in_parts`,
// BLOSUM62's with `u16` offsets written by the stripes, SC9d; word size 5's by the bucket fill in
// parts, SC9e) is the one-part
// build's, slot for slot, for each word size and for BLOSUM62 and another matrix, on segments of
// one table; also a query past `X_U16_LETTERS`. Each table is given back to the pool
// (LOSAT_X_BXCLUTPOOL) with values no table holds, so that the builds after it take vectors
// whose entries they must all write (the pool is process-wide: other tests' vectors may be taken).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn blastx_tables_built_in_parts_are_the_one_part_tables() {
    // The table's fields, and whether its offsets are `u16`.
    type Fields = ((i32, i32, i32, Vec<u32>, Vec<u32>, Vec<u32>, i32), bool);
    fn fields(built: Option<Result<XBxCompressedLookup, ()>>) -> Fields {
        fn of<H: XQueryOffset>(
            table: &mut XFlatCompressedLookup<H>,
            poison: H,
            narrow: bool,
        ) -> Fields {
            let fields = (
                table.word_length,
                table.pv_array_bts,
                table.compressed_alphabet_size,
                table.starts.clone(),
                table.hits.iter().map(|&q| q.offset()).collect(),
                table.pv.clone(),
                table.longest_chain,
            );
            table.starts.fill(u32::MAX);
            table.hits.fill(poison);
            table.pv.fill(u32::MAX);
            (fields, narrow)
        }
        match built {
            Some(Ok(XBxCompressedLookup::Flat(mut table))) => {
                let got = of(&mut table, u16::MAX, true);
                x_give_back_flat(&mut table);
                got
            }
            Some(Ok(XBxCompressedLookup::Wide(mut table))) => {
                let got = of(&mut table, i32::MIN, false);
                x_give_back_wide_flat(&mut table);
                got
            }
            _ => panic!("a flat table"),
        }
    }
    let build = |matrix: ScoringMatrix, word: usize, query: &[u8], parts: usize| {
        let scaled = x_compressed_scaled_threshold(match word {
            5 => 11.0,
            6 => 21.0,
            _ => 20.25,
        });
        let segments = [(3, 1_499), (1_501, 2_998), (3_001, 4_399)];
        let tables = [0..segments.len()];
        match (matrix, word) {
            (ScoringMatrix::Blosum62, 5) => {
                x_bx_build_compressed_lookup::<5>(scaled, query, &segments, &tables, parts)
            }
            (ScoringMatrix::Blosum62, 6) => {
                x_bx_build_compressed_lookup::<6>(scaled, query, &segments, &tables, parts)
            }
            (ScoringMatrix::Blosum62, _) => {
                x_bx_build_compressed_lookup::<7>(scaled, query, &segments, &tables, parts)
            }
            _ => x_bx_build_compressed_lookup_of(
                matrix, word, scaled, query, &segments, &tables, parts,
            )
            .map(|built| {
                built.map(|table| match table {
                    XCompressedLookup::Reference(table) => XBxCompressedLookup::Reference(table),
                    XCompressedLookup::Flat(table) => XBxCompressedLookup::Wide(table),
                })
            }),
        }
    };
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    let query = random_subject(&mut rng, 4_400, 0);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("test pool");
    for matrix in [ScoringMatrix::Blosum62, ScoringMatrix::Blosum45] {
        for word in 5..=7 {
            let want = fields(build(matrix, word, &query, 1));
            assert_eq!(want.1, matrix == ScoringMatrix::Blosum62);
            for parts in [2, 3, 5, 2] {
                let got = pool.install(|| fields(build(matrix, word, &query, parts)));
                assert!(
                    got == want,
                    "{matrix:?} word {word} in {parts} parts: not the one-part table"
                );
            }
        }
    }
    let long = random_subject(&mut rng, X_U16_LETTERS + 2_000, 0);
    let all = [(0, long.len() as i32 - 2)];
    let wide = |parts| {
        fields(x_bx_build_compressed_lookup::<6>(
            2100,
            &long,
            &all,
            &[0..1],
            parts,
        ))
    };
    let want = wide(1);
    assert!(
        pool.install(|| wide(3)) == want,
        "a query past X_U16_LETTERS in 3 parts"
    );
}

// LOSAT_X_BXCLUTPOOL, `i32` offsets (BLASTX's query batches of more than `X_U16_LETTERS`
// letters): a table built on vectors another table and another build left (smaller and
// larger tables in turn, two alive at once) equals the table built without the pool, in every
// field. The pool is process-wide, so other tests' vectors may be among those taken.
#[test]
fn pooled_wide_flat_tables_equal_the_unpooled_ones_on_kept_vectors() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    let build = |query: &[u8], pooled: bool| {
        let locations = [(0, query.len() as i32 - 1)];
        match x_buckets_on(5, 1100, query, &locations, &[0..1], usize::MAX, pooled, 0) {
            XBucketBuild::Built(Some(XFlatBuild::Table(table))) => table,
            _ => panic!("the bucket fill builds this table"),
        }
    };
    let mut alive = None;
    for (round, length) in [3_000, 40, 9_000, 200, 20_000, 7, 1_500, 20_000]
        .into_iter()
        .enumerate()
    {
        let query = random_subject(&mut rng, length, 0);
        let want = build(&query, false);
        let got = build(&query, true);
        x_assert_same_flat_tables(&got, &want);
        assert_eq!(got.pv, want.pv, "round {round}");
        assert_eq!(got.hits, want.hits, "round {round}");
        assert_eq!(got.starts, want.starts, "round {round}");
        assert_eq!(got.longest_chain, want.longest_chain, "round {round}");
        assert_eq!(got.pv_array_bts, want.pv_array_bts, "round {round}");
        match round % 3 {
            0 => {
                if let Some(mut kept) = alive.replace(got) {
                    x_give_back_wide_flat(&mut kept);
                }
            }
            // Values no table holds: the build must write every entry it keeps (starts, hits)
            // and clear the PV array.
            _ => {
                let mut poisoned = got;
                poisoned.starts.fill(u32::MAX);
                poisoned.hits.fill(-7);
                poisoned.pv.fill(u32::MAX);
                x_give_back_wide_flat(&mut poisoned);
            }
        }
    }
    if let Some(mut kept) = alive {
        x_give_back_wide_flat(&mut kept);
    }
}

// LOSAT_X_BXCLUTPOOL: a set comes back only against one lent, so tables built without the pool
// do not grow it; it holds at most as many sets as were lent at once.
#[test]
fn the_flat_pool_keeps_only_sets_it_lent() {
    let set = |n: usize| XFlatTableVectors::<u16> {
        starts: vec![0; n],
        hits: Vec::new(),
        pv: Vec::new(),
    };
    let mut pool = XFlatPool {
        kept: Vec::new(),
        lent: 0,
    };
    assert!(pool.take_back(set(1)).is_err(), "nothing lent");
    assert!(pool.kept.is_empty());
    let first = pool.lend();
    let second = pool.lend();
    assert!(first.starts.is_empty() && second.starts.is_empty());
    assert!(pool.take_back(set(2)).is_ok());
    assert!(pool.take_back(set(3)).is_ok());
    assert!(pool.take_back(set(4)).is_err(), "both sets are back");
    assert_eq!(pool.kept.len(), 2);
    assert_eq!(
        pool.lend().starts.len(),
        3,
        "the last one back goes out first"
    );
    assert_eq!((pool.kept.len(), pool.lent), (1, 1));
}
