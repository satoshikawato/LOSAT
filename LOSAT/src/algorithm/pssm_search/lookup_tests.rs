use super::*;
use crate::algorithm::pssm_search::rows::QueryPssm;

fn pssm(rows: Vec<[i32; BLASTAA_SIZE]>) -> QueryPssm {
    QueryPssm { rows }
}

/// The hits of a word (NCBIstdaa letters) in a finalized table, through its presence bit.
fn hits(table: &BlastAaLookupTable<BackboneCell>, word: &[u8]) -> Vec<i32> {
    let mut index = 0usize;
    for &letter in word {
        index = (index << 5) | letter as usize;
    }
    if table.pv[index >> PV_ARRAY_BTS] & (1 << (index & PV_ARRAY_MASK)) == 0 {
        return Vec::new();
    }
    table.get_hits(index).to_vec()
}

// A PSSM of rows that score one letter per position high: each position's word of three
// letters reaches the threshold only with those letters, so the table holds exactly the
// query's words (the offsets in position order), and a word whose sum misses the threshold
// is not added (no exact-word rule).
#[test]
fn words_reaching_the_threshold_are_added_in_position_order() {
    let mut rows = vec![[-4i32; BLASTAA_SIZE]; 6];
    // A C D A C D (NCBIstdaa 1, 3, 4)
    for (pos, letter) in [1usize, 3, 4, 1, 3, 4].into_iter().enumerate() {
        rows[pos][letter] = 5;
    }
    let rows = pssm(rows);
    let table = build_pssm_lookup(&rows, 3, 11);
    assert_eq!(hits(&table, &[1, 3, 4]), vec![0, 3]);
    assert_eq!(hits(&table, &[3, 4, 1]), vec![1]);
    assert_eq!(hits(&table, &[4, 1, 3]), vec![2]);
    assert_eq!(table.longest_chain, 2);
    // 15 >= 11, but one other letter gives 5 + 5 - 4 = 6 < 11.
    assert!(hits(&table, &[1, 3, 3]).is_empty());
    // Threshold 16: no word reaches it (the best is 15), not even the query's own.
    let table = build_pssm_lookup(&rows, 3, 16);
    assert!(table.backbone.iter().all(|cell| cell.num_used == 0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:623-658: the words of each lookup segment, none across a masked
// position; a segment shorter than a word and an inverted one add nothing.
#[test]
fn words_come_from_the_lookup_segments_only() {
    let mut rows = vec![[-4i32; BLASTAA_SIZE]; 6];
    for (pos, letter) in [1usize, 3, 4, 1, 3, 4].into_iter().enumerate() {
        rows[pos][letter] = 5;
    }
    let rows = pssm(rows);
    // Position 3 masked: the segments [0, 2] and [4, 5]; only the word at 0 fits.
    let table = build_pssm_lookup_in(&rows, 3, 11, &[(0, 2), (4, 5)]);
    assert_eq!(hits(&table, &[1, 3, 4]), vec![0]);
    assert!(hits(&table, &[3, 4, 1]).is_empty());
    assert!(hits(&table, &[4, 1, 3]).is_empty());
    // The whole query as one segment is `build_pssm_lookup`'s table.
    let whole = build_pssm_lookup_in(&rows, 3, 11, &[(0, 5)]);
    let all = build_pssm_lookup(&rows, 3, 11);
    assert_eq!(whole.pv, all.pv);
    for word in [[1, 3, 4], [3, 4, 1], [4, 1, 3]] {
        assert_eq!(hits(&whole, &word), hits(&all, &word));
    }
    let none = build_pssm_lookup_in(&rows, 3, 11, &[(6, 5)]);
    assert!(none.backbone.iter().all(|cell| cell.num_used == 0));
}

// Chains of more than three offsets go to the overflow array in index order, with the cell's
// cursor in its first entry.
#[test]
fn long_chains_go_to_the_overflow_array() {
    let mut rows = vec![[-4i32; BLASTAA_SIZE]; 7];
    for row in rows.iter_mut() {
        row[1] = 5;
    }
    let rows = pssm(rows);
    let table = build_pssm_lookup(&rows, 3, 15);
    assert_eq!(hits(&table, &[1, 1, 1]), vec![0, 1, 2, 3, 4]);
    assert_eq!(table.overflow, vec![0, 1, 2, 3, 4]);
    assert_eq!(table.longest_chain, 5);
    assert_eq!(table.query_length, 7);
}

/// The de Bruijn sequence of order 3 over the 28 letters, unwrapped (`debruijn` of the unit
/// test below): every word of 3 letters appears once.
fn debruijn3() -> Vec<u8> {
    fn db(t: usize, p: usize, a: &mut [usize; 7], out: &mut Vec<u8>) {
        if t > 3 {
            if 3 % p == 0 {
                out.extend(a[1..=p].iter().map(|&letter| letter as u8));
            }
            return;
        }
        a[t] = a[t - p];
        db(t + 1, p, a, out);
        for letter in a[t - p] + 1..BLASTAA_SIZE {
            a[t] = letter;
            db(t + 1, t, a, out);
        }
    }
    let mut out = Vec::new();
    db(1, 1, &mut [0; 7], &mut out);
    let wrap: Vec<u8> = out[..2].to_vec();
    out.extend(wrap);
    out
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/aalookup_unit_test.cpp:310-421 testDebruijnPSSM
// NCBI compiles the case out (`#if 0`, "Needs to be fixed to actually use a PSSM"); this port runs
// it on the PSSM table: the query is the de Bruijn sequence of order 3 (every word once), the PSSM
// row of each position BLOSUM62's column of its letter (row `i` of `posMatrix`, column `j`, is the
// entry of `matrix` at row `j` and at the column of the letter of `sequence` at position one past
// `i`), threshold 11, word size 3. Every word's cell holds exactly the query words whose score with
// it reaches 11 (the brute-force count of the test): the PSSM table has no exact-word rule.
#[test]
fn debruijn_pssm_cells_hold_the_brute_force_neighbours() {
    let sequence = debruijn3();
    assert_eq!(sequence.len(), BLASTAA_SIZE.pow(3) + 2);
    let matrix = |q: usize, s: usize| crate::utils::matrix::blosum62_score(q as u8, s as u8);
    let rows = pssm(
        sequence
            .iter()
            .map(|&letter| {
                let mut row = [0i32; BLASTAA_SIZE];
                for (j, cell) in row.iter_mut().enumerate() {
                    *cell = matrix(j, letter as usize);
                }
                row
            })
            .collect(),
    );
    let table = build_pssm_lookup(&rows, 3, 11);
    // count(i, j, k) = #{(x, y, z): m[i][x] + m[j][y] + m[k][z] >= 11}, by the number of z of
    // each row k reaching each threshold.
    let sorted: Vec<Vec<i32>> = (0..BLASTAA_SIZE)
        .map(|k| {
            let mut row: Vec<i32> = (0..BLASTAA_SIZE).map(|z| matrix(k, z)).collect();
            row.sort_unstable();
            row
        })
        .collect();
    let reaching = |k: usize, t: i32| BLASTAA_SIZE - sorted[k].partition_point(|&v| v < t);
    for i in 0..BLASTAA_SIZE {
        for j in 0..BLASTAA_SIZE {
            let mut pairs = Vec::with_capacity(BLASTAA_SIZE * BLASTAA_SIZE);
            for x in 0..BLASTAA_SIZE {
                for y in 0..BLASTAA_SIZE {
                    pairs.push(matrix(i, x).saturating_add(matrix(j, y)));
                }
            }
            for k in 0..BLASTAA_SIZE {
                let count: usize = pairs
                    .iter()
                    .map(|&pair| reaching(k, 11i32.saturating_sub(pair)))
                    .sum();
                let index = (i << 10) | (j << 5) | k;
                assert_eq!(
                    table.backbone[index].num_used as usize, count,
                    "word ({i}, {j}, {k})"
                );
            }
        }
    }
}
