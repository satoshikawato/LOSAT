//! Outfmt 0's masks of an HSP's row (`fill_row_loc_list`) against `fill_loc_list` over the HSP's
//! `CAlnMap` (line RG, item T4a1).

use super::aln_vec::{AlnMap, MapBuffers};
use super::{fill_loc_list, fill_row_loc_list, Denseg, Mask, SeqKey};

/// A small xorshift generator: the cases are the same on every run.
struct Rng(u64);

impl Rng {
    fn below(&mut self, n: u32) -> u32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % u64::from(n)) as u32
    }
}

/// A random Dense-seg of two rows (no segment of two gaps, each row with a residue), its rows'
/// strands, widths and first positions random.
fn random_denseg(rng: &mut Rng) -> Denseg {
    loop {
        let num_segs = 1 + rng.below(6) as usize;
        let mut lens = Vec::new();
        let mut present = Vec::new();
        for _ in 0..num_segs {
            lens.push(1 + rng.below(5));
            present.push(match rng.below(4) {
                0 => [true, false],
                1 => [false, true],
                _ => [true, true],
            });
        }
        if !(0..2).all(|row| present.iter().any(|p| p[row])) {
            continue;
        }
        let widths = [1 + 2 * rng.below(2), 1 + 2 * rng.below(2)];
        let plus = [rng.below(2) == 0, rng.below(2) == 0];
        let mut starts = vec![-1i32; num_segs * 2];
        for row in 0..2 {
            let first = rng.below(20);
            let total: u32 = (0..num_segs)
                .filter(|&s| present[s][row])
                .map(|s| lens[s] * widths[row])
                .sum();
            let mut done = 0;
            for s in 0..num_segs {
                if !present[s][row] {
                    continue;
                }
                let size = lens[s] * widths[row];
                // A minus-strand row runs from its last position down.
                let start = if plus[row] {
                    first + done
                } else {
                    first + total - done - size
                };
                starts[s * 2 + row] = start as i32;
                done += size;
            }
        }
        return Denseg {
            ids: vec![SeqKey::Query, SeqKey::Subject(0)],
            lens,
            starts,
            plus: plus.to_vec(),
            widths: if widths == [1, 1] {
                Vec::new()
            } else {
                widths.to_vec()
            },
        };
    }
}

#[test]
fn row_loc_list_matches_fill_loc_list() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut buffers = MapBuffers::default();
    // One buffer of gap runs and one list for every row, as a reused `HspMasks` keeps them.
    let (mut gaps, mut list) = (Vec::new(), Vec::new());
    for case in 0..20_000 {
        let ds = random_denseg(&mut rng);
        let map = AlnMap::new_in(&ds, &mut buffers).unwrap();
        for row in 0..2 {
            let (lo, hi) = (map.seq_start(row), map.seq_stop(row));
            // Masks inside, outside and across the row's ends, adjacent, of one position, at
            // its first position (column 0 or the row's first residue), with frames.
            let mut masks: Vec<Mask> = Vec::new();
            for _ in 0..1 + rng.below(4) {
                let from = (lo + rng.below(hi - lo + 12)).saturating_sub(6);
                let to = from + rng.below(8);
                masks.push((from, to, rng.below(3) as i8 - 1));
            }
            if let Some(&(_, to, frame)) = masks.last() {
                masks.push((to + 1, to + 1 + rng.below(3), frame));
            }
            masks.push((lo, lo, 0));
            let letters: Vec<u8> = (0..ds.lens.len())
                .flat_map(|s| {
                    let letter = if ds.starts[s * 2 + row] >= 0 {
                        b'A'
                    } else {
                        b'-'
                    };
                    std::iter::repeat_n(letter, ds.lens[s] as usize)
                })
                .collect();
            list.clear();
            fill_row_loc_list(
                &letters,
                |column| letters[column] != b'-',
                (lo, hi),
                map.is_positive_strand(row),
                map.width(row),
                masks.iter().copied(),
                &mut gaps,
                &mut list,
            );
            let expected = fill_loc_list(&map, map.seq_key(row), None, &masks);
            assert_eq!(list, expected, "case {case} row {row} {masks:?}");
        }
    }
}

// No NCBI counterpart: `MaskSeek::for_each_meeting` against the scan of the whole list (the masks
// with `to >= low && from <= high`, what `fill_row_loc_list` keeps of a list), for lists in any
// order with overlaps, nested masks and repeats; sorted by `index` they are in the list's order.
#[test]
fn mask_seek_keeps_the_full_scan() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    for _ in 0..5_000 {
        let masks: Vec<(u32, u32)> = (0..rng.below(12))
            .map(|_| {
                let from = rng.below(60);
                let longest = if rng.below(4) == 0 { 40 } else { 6 };
                (from, from + rng.below(longest))
            })
            .collect();
        let seek = super::MaskSeek::new(masks.iter().copied());
        assert_eq!(seek.ranges().len(), masks.len());
        for _ in 0..8 {
            let low = rng.below(70);
            let high = low + rng.below(30);
            let mut found: Vec<(u32, (u32, u32))> = Vec::new();
            seek.for_each_meeting(low, high, |from, to, index| found.push((index, (from, to))));
            found.sort_unstable_by_key(|&(index, _)| index);
            let scan: Vec<(u32, (u32, u32))> = masks
                .iter()
                .copied()
                .zip(0..)
                .filter(|&((from, to), _)| to >= low && from <= high)
                .map(|(mask, index)| (index, mask))
                .collect();
            assert_eq!(found, scan, "masks {masks:?}, row {low}..={high}");
        }
    }
}
