use super::{
    build_compressed_query, exact_match_extend_left, exact_match_extend_right, SmallNaWord,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1330-1347
const NCBI_LEFT: [i64; 256] = [
    4, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
];
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1353-1358
// (64 initializers; C sets the other entries to 0.)
const NCBI_RIGHT_INITIALIZED: [i64; 64] = [
    4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
];

#[test]
fn exact_match_counts_are_ncbi_tables() {
    for xor in 0..=255u8 {
        assert_eq!(
            exact_match_extend_left(xor),
            NCBI_LEFT[xor as usize],
            "{xor}"
        );
        let right = NCBI_RIGHT_INITIALIZED
            .get(xor as usize)
            .copied()
            .unwrap_or(0);
        assert_eq!(exact_match_extend_right(xor), right, "{xor}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:473-499 (the loops,
// with `curr_letter`, an `Int4`, stored into `Uint1` bytes).
fn ncbi_compress(old: &[u8]) -> Vec<u8> {
    let len = old.len() as isize;
    let at = |index: isize| (index + 3) as usize;
    let mut new_seq = vec![0u8; old.len() + 3];
    let max_start = len.min(3);
    let mut curr_letter = 0u32;
    let mut i = 0;
    while i < len {
        curr_letter = curr_letter << 2 | u32::from(old[i as usize] & 3);
        new_seq[at(i - max_start)] = curr_letter as u8;
        i += 1;
    }
    for i in 0..max_start {
        curr_letter <<= 2;
        new_seq[at(len - (max_start - i))] = curr_letter as u8;
    }
    new_seq
}

#[test]
fn the_compressed_query_is_ncbis() {
    let mut state = 7u32;
    for len in 3..64 {
        let query: Vec<u8> = (0..len)
            .map(|_| {
                state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                ((state >> 16) % 16) as u8
            })
            .collect();
        assert_eq!(
            build_compressed_query(&query),
            ncbi_compress(&query),
            "{query:?}"
        );
    }
}

fn packed(bases: &[u8]) -> Vec<u8> {
    bases
        .chunks(4)
        .map(|chunk| {
            (0..4).fold(0u8, |byte, index| {
                byte << 2 | chunk.get(index).copied().unwrap_or(0)
            })
        })
        .collect()
}

// The fourteenth S07+ audit round: the query `K` (blastna 7) reads as `T` in the
// compressed query, so the exact-match extension of a lookup hit runs through it and the
// word starts before it; the per-base extension of the other tables stops at it.
#[test]
fn the_small_table_extension_reads_ambiguity_codes_as_bases() {
    let code = |text: &str| -> Vec<u8> {
        text.bytes()
            .map(|letter| match letter {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                b'K' => 7,
                b'D' => 11,
                _ => unreachable!(),
            })
            .collect()
    };
    let query = code("TTTCGTTGACCTAAAAAGTTCKGTTGGTATGAGAGAAGAATTTTTGGTTGDCAGAAAAAC");
    let subject =
        code("AAAGCGGCATTTCGTTGACCTAAAAAGTTCTGTTGGTATGAGAGAAGAATTTTTGGTTGACAGAAAAACGGCATTACGA");
    let compressed = build_compressed_query(&query);
    let subject_packed = packed(&subject);
    let word = SmallNaWord {
        compressed_query: &compressed,
        subject_packed: &subject_packed,
        query_length: query.len(),
        word_length: 28,
        lut_word_length: 8,
    };
    // The lookup hit of query 12 at subject 21 (scan step 21) starts its word at the
    // query start.
    assert_eq!(
        word.extend(12, 21, 0, query.len(), subject.len()),
        Some((0, 9))
    );
    // Query 33 at subject 42 runs back through the K to query 13.
    assert_eq!(
        word.extend(33, 42, 0, query.len(), subject.len()),
        Some((13, 22))
    );
    // A word that would end past the query context (here at 40) is dropped.
    assert_eq!(word.extend(33, 42, 0, 40, subject.len()), None);
}
