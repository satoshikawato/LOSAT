use super::super::x_mb_lookup_par::tests::{fill_old, Lcg};
use super::super::x_mb_lookup_par::Tables;
use super::super::{build_two_stage_lookup, compute_mb_pv_params, MaskedInterval};
use super::*;

/// Queries (BLASTNA), their offsets in the concatenation and their masks.
pub(crate) struct Case {
    pub(crate) queries: Vec<Vec<u8>>,
    pub(crate) offsets: Vec<i32>,
    pub(crate) masks: Vec<Vec<MaskedInterval>>,
}

impl Case {
    pub(crate) fn refs(&self) -> Vec<&[u8]> {
        self.queries.iter().map(Vec::as_slice).collect()
    }

    pub(crate) fn words<'a>(
        &'a self,
        refs: &'a [&'a [u8]],
        word_length: usize,
        lut_word_length: usize,
        ascending_cells: bool,
    ) -> Words<'a> {
        Words {
            queries_blastna: refs,
            query_offsets: &self.offsets,
            query_masks: &self.masks,
            word_length,
            lut_word_length,
            kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
            ascending_cells,
        }
    }
}

/// EXPERIMENT (LOSAT_X_BNCSRPACK): `case` with four more queries whose words fill a few cells
/// with more words than a byte counts: 1,200 letters A (one cell), 800 letters ACAC... (two
/// cells), and one letter repeated so that its cell gets exactly 255 (T) or 256 (G) words of a
/// table of `lut` letters.
pub(crate) fn x_with_repeats(case: &mut Case, lut: usize) {
    for seq in [
        vec![0u8; 1200],
        (0..800).map(|i| (i % 2) as u8).collect(),
        vec![3u8; 255 + lut - 1],
        vec![2u8; 256 + lut - 1],
    ] {
        x_push_query(case, seq);
    }
}

/// `case` with `seq` as one more query, unmasked, after the last query and its separator.
pub(crate) fn x_push_query(case: &mut Case, seq: Vec<u8>) {
    let at = match (case.offsets.last(), case.queries.last()) {
        (Some(&o), Some(q)) => o + q.len() as i32 + 1,
        _ => 0,
    };
    case.offsets.push(at);
    case.queries.push(seq);
    case.masks.push(Vec::new());
}

/// `n` random queries of about `len` letters: few distinct letters in half of them (so cells
/// get long chains), ambiguities (N), masked intervals, and an empty query, a query of N only
/// and a query of 5 letters (shorter than every table's words) among them.
pub(crate) fn random_case(rng: &mut Lcg, n: usize, len: usize) -> Case {
    let mut queries = Vec::new();
    let mut offsets = Vec::new();
    let mut masks = Vec::new();
    let mut at = 0i32;
    for q in 0..n + 3 {
        let seq: Vec<u8> = match q {
            0 => Vec::new(),
            1 => vec![14; 40],
            2 => vec![0, 1, 2, 3, 0],
            _ => {
                let l = len / 2 + rng.below(len as u64 + 1) as usize;
                let letters = if q % 2 == 0 { 2 } else { 4 };
                (0..l)
                    .map(|_| {
                        if rng.below(150) == 0 {
                            14
                        } else {
                            rng.below(letters) as u8
                        }
                    })
                    .collect()
            }
        };
        let mut m = Vec::new();
        if q % 3 == 1 && seq.len() > 60 {
            let s = rng.below((seq.len() - 40) as u64) as usize;
            m.push(MaskedInterval {
                start: s,
                end: s + 1 + rng.below(30) as usize,
            });
        }
        offsets.push(at);
        at += seq.len() as i32 + 1;
        queries.push(seq);
        masks.push(m);
    }
    Case {
        queries,
        offsets,
        masks,
    }
}

/// The table of `cells` cells as the loop fills it (`fill_old` with one part: the loop's
/// statements, re-threaded for `ascending_cells`): `(hashtable, next_pos, pv_array,
/// helper_array)`.
pub(crate) fn chains(
    case: &Case,
    words: &Words<'_>,
    cells: usize,
    pv_len: usize,
    pv_array_bts: usize,
) -> (Vec<u32>, Vec<u32>, Vec<PvArrayType>, Vec<u32>) {
    let query_length = case
        .queries
        .iter()
        .zip(&case.offsets)
        .map(|(seq, &off)| off.max(0) as usize + seq.len())
        .max()
        .unwrap_or(0);
    let mut hashtable = vec![0u32; cells];
    let mut next_pos = vec![0u32; query_length + 1];
    let mut pv = vec![0 as PvArrayType; pv_len];
    let mut helper = vec![0u32; (cells >> 11).max(1)];
    fill_old(
        1,
        words,
        Tables {
            hashtable: &mut hashtable,
            next_pos: &mut next_pos,
            pv_array: &mut pv,
            pv_array_bts,
            helper_array: &mut helper,
            x_h24: false,
        },
    );
    (hashtable, next_pos, pv, helper)
}

/// The presence-vector length and shift of a table of `lut_word_length` for `case`.
pub(crate) fn pv_params(case: &Case, lut_word_length: usize) -> (usize, usize) {
    let approx = case.queries.iter().map(Vec::len).sum();
    compute_mb_pv_params(1 << (2 * lut_word_length), approx, false, lut_word_length)
}

// The serial build lists the loop's chains, newest first and in reading order, for every
// `lut_word_length` of the rule, with the fill's presence vector, counters and longest chain
// (`x_assert_runs_are_the_chains`: every cell's run, `contains` on each run).
#[test]
fn the_runs_are_the_chains_in_both_orders() {
    let mut rng = Lcg(0xB4C5_0001);
    for lut in MIN_LUT..=MAX_LUT {
        for ascending in [false, true] {
            for (n, len) in [(1, 6), (3, 300), (5, 6000)] {
                let case = random_case(&mut rng, n, len);
                let refs = case.refs();
                let word_length = lut + rng.below(6) as usize;
                let words = case.words(&refs, word_length, lut, ascending);
                let cells = 1usize << (2 * lut);
                let (pv_len, bts) = pv_params(&case, lut);
                let (hashtable, next_pos, fill_pv, fill_helper) =
                    chains(&case, &words, cells, pv_len, bts);
                let mut pv = vec![0 as PvArrayType; pv_len];
                let mut helper = vec![0u32; cells >> 11];
                let csr = serial(&words, cells, &mut pv, bts, &mut helper, Vec::new());
                x_assert_runs_are_the_chains(
                    &csr,
                    &pv,
                    &helper,
                    &hashtable,
                    &next_pos,
                    &fill_pv,
                    &fill_helper,
                    &format!("lut {lut} ascending {ascending} n {n} len {len}"),
                );
            }
        }
    }
}

// LOSAT_X_BNCSRPACK: the packed head and the head of one offset per cell list the loop's chains
// (newest first), with the fill's presence vector, counters and longest chain, whatever the
// switch reads; many queries make tables with most cells non-empty.
#[test]
fn the_packed_head_lists_the_chains_too() {
    let mut rng = Lcg(0xB4C5_0002);
    for lut in MIN_LUT..=12 {
        // LOSAT_X_MBCSR's 4^12 cells: words past a byte's count, and the batch of 800,000 table
        // entries or more whose presence vector has half the words.
        let cases: &[(usize, usize, bool)] = if lut == 12 {
            &[(40, 3000, false), (2, 500, true), (1, 900_000, false)]
        } else {
            &[
                (1, 6, false),
                (3, 300, false),
                (5, 6000, false),
                (40, 3000, false),
                (2, 500, true),
            ]
        };
        for &(n, len, repeats) in cases {
            let mut case = random_case(&mut rng, n, len);
            if repeats {
                x_with_repeats(&mut case, lut);
            }
            let refs = case.refs();
            let word_length = lut + rng.below(6) as usize;
            let words = case.words(&refs, word_length, lut, false);
            let cells = 1usize << (2 * lut);
            let (pv_len, bts) = pv_params(&case, lut);
            let (hashtable, next_pos, fill_pv, fill_helper) =
                chains(&case, &words, cells, pv_len, bts);
            for packed in [false, true] {
                let mut pv = vec![0 as PvArrayType; pv_len];
                let mut helper = vec![0u32; cells >> 11];
                let build = if packed {
                    x_serial_packed
                } else {
                    x_serial_unpacked
                };
                let csr = build(&words, cells, &mut pv, bts, &mut helper, Vec::new());
                let what = format!("lut {lut} packed {packed} n {n} len {len}");
                if repeats {
                    // the cells of 255, 256 and more words
                    let mut runs: Vec<usize> = (0..cells).map(|c| csr.x_run_any(c).len()).collect();
                    runs.sort_unstable();
                    assert!(runs.contains(&255) && runs.contains(&256), "{what}");
                    assert!(runs[cells - 1] >= 1200 - lut, "{what}");
                }
                x_assert_runs_are_the_chains(
                    &csr,
                    &pv,
                    &helper,
                    &hashtable,
                    &next_pos,
                    &fill_pv,
                    &fill_helper,
                    &what,
                );
            }
        }
    }
}

#[test]
fn contains_searches_runs_of_either_order() {
    for (run, q, want) in [
        (&[][..], 1, false),
        (&[7][..], 7, true),
        (&[7][..], 8, false),
        (&[3, 9, 40][..], 9, true),
        (&[3, 9, 40][..], 3, true),
        (&[3, 9, 40][..], 40, true),
        (&[3, 9, 40][..], 10, false),
        (&[40, 9, 3][..], 9, true),
        (&[40, 9, 3][..], 40, true),
        (&[40, 9, 3][..], 3, true),
        (&[40, 9, 3][..], 4, false),
        (&[40, 9][..], 0, false),
    ] {
        assert_eq!(contains(run, q), want, "{run:?} {q}");
    }
}

#[test]
fn the_rule_takes_a_quarter_of_the_cells_and_lut_8_to_11() {
    assert!(!wanted(0, 1 << 22));
    assert!(!wanted((1 << 20) - 1, 1 << 22));
    assert!(wanted(1 << 20, 1 << 22));
    assert!(wanted(usize::MAX, 1 << 22));
    // The positions are those the loop reads: masked ranges, queries shorter than
    // `lut_word_length` and ranges shorter than `word_length` do not count.
    let case = Case {
        queries: vec![vec![0; 100], vec![1; 7], vec![2; 30]],
        offsets: vec![0, 101, 109],
        masks: vec![
            vec![MaskedInterval { start: 10, end: 20 }],
            Vec::new(),
            vec![MaskedInterval { start: 12, end: 13 }],
        ],
    };
    let refs = case.refs();
    // ranges 0..10, 20..100 (both >= 9), the 7-letter query (< lut 8), 0..12 and 13..30
    assert_eq!(
        positions(&case.words(&refs, 9, 8, false)),
        10 + 80 + 12 + 17
    );
    // word size 12: the range of 10 letters is dropped, the one of 12 is not
    assert_eq!(positions(&case.words(&refs, 12, 8, false)), 80 + 12 + 17);
    if on() {
        let big = Case {
            queries: vec![vec![0; 20_000]],
            offsets: vec![0],
            masks: vec![Vec::new()],
        };
        let refs = big.refs();
        assert!(plan(
            &big.words(&refs, 9, 8, true),
            1 << 16,
            false,
            false,
            None
        ));
        assert!(!plan(
            &big.words(&refs, 9, 8, true),
            1 << 16,
            true,
            false,
            None
        ));
        assert!(!plan(
            &big.words(&refs, 9, 8, true),
            1 << 16,
            false,
            true,
            None
        ));
        assert!(!plan(
            &big.words(&refs, 9, 7, true),
            1 << 14,
            false,
            false,
            None
        ));
        assert!(!plan(
            &big.words(&refs, 13, 12, false),
            1 << 24,
            false,
            false,
            None
        ));
        // LOSAT_X_MBCSR: 4^12 cells with a search of at least 8 letters per position and one
        // per cell.
        let positions = big.refs()[0].len() as u64;
        assert_eq!(
            plan(
                &big.words(&refs, 13, 12, false),
                1 << 24,
                false,
                false,
                Some(positions * X_MB12_LETTERS_PER_POSITION + (1 << 24))
            ),
            x_will_pack(false) && x_mb12_on()
        );
        assert!(!plan(
            &big.words(&refs, 13, 12, false),
            1 << 24,
            false,
            false,
            Some(positions * X_MB12_LETTERS_PER_POSITION + (1 << 24) - 1)
        ));
        // below the rule: chains, unless the head is packed (LOSAT_X_BNCSRPACK has no rule)
        assert_eq!(
            plan(
                &big.words(&refs, 11, 11, false),
                1 << 22,
                false,
                false,
                None
            ),
            x_will_pack(false)
        );
    }
}

// A table above the rule, built by `build_two_stage_lookup` off and on a pool, is laid out as
// CSR (with the switch on) and answers `for_each_hit`, `count_hits`, `contains_hit` and
// `longest_chain` as the loop's chains, in both orders.
#[test]
fn a_table_above_the_rule_answers_from_its_csr() {
    let mut rng = Lcg(0xB4C5_0002);
    // lut 8 in reading order (16,384 positions for 4^8 cells) and lut 9 newest first (65,536
    // for 4^9)
    for (lut, ascending, len) in [(8usize, true, 12_000usize), (9, false, 40_000)] {
        let case = random_case(&mut rng, 3, len);
        let refs = case.refs();
        let word_length = lut + 1;
        let words = case.words(&refs, word_length, lut, ascending);
        assert!(wanted(positions(&words), 1 << (2 * lut)));
        let cells = 1usize << (2 * lut);
        let (pv_len, bts) = pv_params(&case, lut);
        let (hashtable, next_pos, _, fill_helper) = chains(&case, &words, cells, pv_len, bts);
        let approx = case.queries.iter().map(Vec::len).sum();
        let check = |what: &str| {
            let lookup = build_two_stage_lookup(
                &refs,
                &case.offsets,
                word_length,
                lut,
                &case.masks,
                None,
                0,
                approx,
                ascending,
            );
            assert_eq!(lookup.x_mb_csr().is_some(), on(), "{what}");
            let longest = fill_helper.iter().fold(2usize, |m, &v| m.max(v as usize));
            assert_eq!(lookup.longest_chain(), longest, "{what}");
            let mut chain = Vec::new();
            for (cell, &head) in hashtable.iter().enumerate() {
                chain.clear();
                let mut q = head;
                while q != 0 {
                    chain.push(q);
                    q = next_pos[q as usize];
                }
                let mut seen = Vec::new();
                lookup.for_each_hit(cell as u64, |q| seen.push(q));
                assert_eq!(seen, chain, "{what} cell {cell}");
                assert_eq!(lookup.count_hits(cell as u64), chain.len(), "{what}");
                for &q in chain.iter().take(3) {
                    assert!(lookup.contains_hit(cell as u64, q), "{what} cell {cell}");
                }
                assert!(!lookup.contains_hit(cell as u64, 1 << 30), "{what}");
                if let Some(csr) = lookup.x_mb_csr() {
                    assert_eq!(csr.x_run_any(cell), chain.as_slice(), "{what} cell {cell}");
                }
            }
        };
        check(&format!("lut {lut} off a pool"));
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        crate::utils::threading::with_search_pool(4, "blastn", |_pool| {
            check(&format!("lut {lut} on a pool of 4"));
            Ok(())
        })
        .unwrap();
    }
}
