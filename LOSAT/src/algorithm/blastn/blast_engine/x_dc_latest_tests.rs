//! Tests of `x_dc_latest.rs` (LOSAT_X_DCLATEST): the map against a model of each diagonal's
//! latest insert, and against the diagonal hash's chains (NCBI's window and LOSAT's).

use super::super::super::x_diag_lv::LvSeeds;
use super::super::super::{diag_hash_insert_window, DiagHashTable};
use super::{
    eq_mask, lt_mask, scalar_mask, DiagLatest, Found, MIN_BITS, MISS, SLOTS, STASH, STASHED,
};
use std::collections::HashMap;

/// A xorshift generator (the tests' inputs only).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// A value in `0..n`.
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// The windows of the tests for a word length: the band 2w + 2 < wl and around it, around the
/// word length, and the defaults (40 for dc-megablast) and above.
fn windows(wl: usize) -> Vec<usize> {
    let mut w: Vec<usize> = (1..=10).collect();
    w.extend([wl - 1, wl, wl + 1, 40, 64]);
    w.sort_unstable();
    w.dedup();
    w
}

/// What `find` returned as the retrieve's triple (`hit_len` is read by nothing, module comment
/// (i); no entry gives level 0 and no saved hit, which decide as no cell).
fn as_cell(found: Found) -> Option<(i32, i32, bool)> {
    Some((found.level, 0, found.saved))
}

/// Diagonals that share their home bucket at `bits` 6 to 10 (their products with the constant
/// share the top 10 bits), so that they fill one bucket and the stash at every size up to 1024
/// buckets.
fn crowded(n: usize) -> Vec<i32> {
    let mut out = Vec::new();
    let mut d: i32 = -5_000_000;
    let top = |d: i32| (d as u32).wrapping_mul(0x9E37_79B1) >> 22;
    let want = top(d);
    while out.len() < n {
        if top(d) == want {
            out.push(d);
        }
        d += 1;
    }
    out
}

// Every seed's decisions (rejected, and the two-hit branch: `LvSeeds::decisions`) are those of
// the latest insert of its diagonal, and `find` returns that insert itself while it is live: over
// stages of non-decreasing positions, the offset moved as NCBI moves it after each chunk and
// sometimes cleared (the model empty then, as NCBI's table), with random diagonals, diagonals that
// crowd one home bucket (the stash and the growths), levels from the seed's position to the
// chunk's end, and two seeds at one position.
#[test]
fn the_map_decides_as_the_latest_insert_of_each_diagonal() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let mut grows = 0u64;
    let mut stash_reads = 0u64;
    for wl in [16usize, 18, 21] {
        for w in windows(wl) {
            let seeds = LvSeeds::new(true, wl, w);
            let mut map = DiagLatest::default();
            let mut model: HashMap<i32, (i32, bool)> = HashMap::new();
            let crowd = crowded(80);
            let spread: Vec<i32> = (0..600)
                .map(|_| rng.below(4_000_000) as i32 - 2_000_000)
                .collect();
            let mut offset = w;
            for stage in 0..12 {
                let s_len = 200 + rng.below(4_000) as usize;
                map.begin(seeds.bound(offset));
                let mut s = offset;
                for _ in 0..2_000 {
                    s = (s + rng.below(4) as usize).min(offset + s_len - 1);
                    let diag = match rng.below(3) {
                        0 => crowd[rng.below(crowd.len() as u64) as usize],
                        _ => spread[rng.below(spread.len() as u64) as usize],
                    };
                    let bound = seeds.bound(s);
                    assert_ne!(bound, i32::MIN);
                    let found = map.find(diag, s as i32);
                    stash_reads += u64::from(found.at >= STASHED);
                    let latest = model.get(&diag).copied();
                    let expected = seeds.decisions(latest.map(|(l, v)| (l, 0, v)), s);
                    let got = seeds.decisions(as_cell(found), s);
                    assert_eq!(
                        got, expected,
                        "wl {wl} w {w} stage {stage} s {s} diag {diag}"
                    );
                    if let Some((level, saved)) = latest {
                        if level >= bound {
                            assert_eq!((found.level, found.saved), (level, saved));
                        }
                    }
                    if got.0 {
                        continue;
                    }
                    // a level at least the position, at most the chunk's end (an insert's `s_end_pos`
                    // or the end of its ungapped extension)
                    let end = offset + s_len;
                    let level = match rng.below(4) {
                        0 => s,
                        1 => s + 1,
                        _ => s + wl + rng.below(300) as usize,
                    }
                    .min(end) as i32;
                    let saved = rng.below(2) == 0;
                    map.put(found, diag, level, saved, s as i32, bound);
                    model.insert(diag, (level, saved));
                }
                if rng.below(5) == 0 {
                    // NCBI's clear (at a quarter of the largest Int4): the offset back to the window, the table empty
                    offset = w;
                    model.clear();
                } else {
                    offset += s_len + w;
                }
            }
            grows += u64::from(map.bits > MIN_BITS);
        }
    }
    assert!(grows > 0, "no map grew");
    assert!(stash_reads > 0, "no seed read the stash");
}

// `begin`: a map whose levels may reach the stage's first bound is emptied (after NCBI's clear);
// one whose levels all lie below it keeps them, dead for the whole stage.
#[test]
fn begin_empties_the_map_only_where_a_level_may_be_live() {
    let seeds = LvSeeds::new(true, 18, 40);
    let mut map = DiagLatest::default();
    map.begin(seeds.bound(40));
    let found = map.find(7, 100);
    map.put(found, 7, 500, true, 100, seeds.bound(100));
    assert_eq!(map.top, 500);
    // the next chunk starts past every level: kept, and dead
    map.begin(seeds.bound(600));
    assert_eq!(map.entries, 1);
    let found = map.find(7, 600);
    assert_eq!((found.level, found.saved), (500, true));
    assert!(found.level < seeds.bound(600));
    // after the offset's clear the old level lies above the stage's positions: emptied
    map.begin(seeds.bound(40));
    assert_eq!(map.entries, 0);
    assert_eq!(map.top, i32::MIN);
    let found = map.find(7, 40);
    assert_eq!((found.level, found.saved, found.at), (0, false, MISS));
    // no bound: emptied too
    map.put(found, 7, 60, false, 40, seeds.bound(40));
    map.begin(i32::MIN);
    assert_eq!(map.entries, 0);
}

// The stash: diagonals crowding one home bucket go there once its slots are taken by live
// entries; each is found while live, a dead one's place is taken again, and past the bucket's
// `stash_until` the stash is not read. The map grows while it has fewer than two buckets per
// entry and keeps every live entry; past that the stash takes the rest.
#[test]
fn crowded_diagonals_go_to_the_stash_and_the_map_grows_with_its_entries() {
    let seeds = LvSeeds::new(true, 16, 40);
    let crowd = crowded(200);
    let mut map = DiagLatest::default();
    map.begin(seeds.bound(40));
    let s = 1_000i32;
    let bound = seeds.bound(s as usize);
    for (k, &diag) in crowd.iter().enumerate() {
        let found = map.find(diag, s);
        assert_eq!(found.level, 0, "diagonal {k} found before its insert");
        map.put(found, diag, 2_000 + k as i32, k % 3 == 0, s, bound);
        for (j, &d) in crowd[..=k].iter().enumerate() {
            let found = map.find(d, s);
            assert_eq!(
                (found.level, found.saved),
                (2_000 + j as i32, j % 3 == 0),
                "after {k} inserts, diagonal {j}"
            );
        }
    }
    // all 200 live: 7 in the home slots at most, the rest in the stash, and the buckets doubled
    // while fewer than twice the entries
    assert!(map.stash.len() >= 200 - SLOTS as usize);
    assert!(map.buckets.len() >= 2 * STASH);
    assert!(
        map.buckets.len() <= 4 * map.entries,
        "{} buckets",
        map.buckets.len()
    );
    // later: every entry dead (levels below the bound), so a new diagonal takes a dead place and
    // the stash is not read past `stash_until`
    let s2 = 3_000i32;
    let bound2 = seeds.bound(s2 as usize);
    assert!(crowd.len() as i32 + 2_000 < bound2);
    let found = map.find(crowd[150], s2);
    assert_eq!(found.at, MISS, "a stash read past stash_until");
    let entries = map.entries;
    let fresh = crowded(201)[200];
    let found = map.find(fresh, s2);
    map.put(found, fresh, 3_500, true, s2, bound2);
    assert_eq!(map.entries, entries, "a dead place was not taken again");
    let found = map.find(fresh, s2);
    assert_eq!((found.level, found.saved), (3_500, true));
}

/// One run of seeds at non-decreasing positions on `diags` against both the map and a diagonal
/// hash whose insert takes over at `window` (`s_off - level > window`), each structure taking its
/// own decisions. Returns the seeds whose decisions differed and those where the chains returned
/// another cell than the map (an older insert, or none: a cell taken over).
fn against_chains(wl: usize, w: usize, window: i32, rng: &mut Rng) -> (u64, u64) {
    let seeds = LvSeeds::new(true, wl, w);
    let mut chains = DiagHashTable::new(w as i32);
    let mut chains_map = DiagLatest::default();
    chains_map.begin(seeds.bound(w));
    // most diagonals share NCBI's bucket (the diagonal modulo 512), so that inserts take over
    // each other's cells
    let diags: Vec<i32> = (0..40)
        .map(|k| 37 + 512 * k - 512 * 20)
        .chain((0..8).map(|k| 1_000 + k))
        .collect();
    let (mut differ, mut other) = (0u64, 0u64);
    let mut s = w;
    for _ in 0..6_000 {
        s += rng.below(2 * w as u64 + wl as u64) as usize;
        let diag = diags[rng.below(diags.len() as u64) as usize];
        let bound = seeds.bound(s);
        let cell = chains.retrieve(diag);
        let found = chains_map.find(diag, s as i32);
        let by_chains = seeds.decisions(cell, s);
        let by_map = seeds.decisions(as_cell(found), s);
        differ += u64::from(by_chains != by_map);
        let map_cell = (found.level, found.saved);
        if cell.map(|(l, _, v)| (l, v)).unwrap_or((0, false)) != map_cell {
            other += 1;
        }
        let level = (s + wl + rng.below(30) as usize) as i32;
        let hit_len = rng.below(40) as i32;
        let saved = rng.below(3) == 0;
        if !by_chains.0 {
            chains.insert_plain(diag, level, hit_len, saved, s as i32, window);
        }
        if !by_map.0 {
            chains_map.put(found, diag, level, saved, s as i32, bound);
        }
    }
    (differ, other)
}

// Against the chains of NCBI's diagonal hash: with NCBI's window in the scope (w + 1, module
// comment (iii)) the map makes every seed's decisions for every window, while the chains return
// other cells than the map (cells taken over); with LOSAT's window (`diag_hash_insert_window`,
// 2w - wl + 1 below the word length) the same outside the band 2w + 2 < wl, and inside it the
// decisions differ for some seeds (module comment (v): no seed is extended there).
#[test]
fn the_map_decides_as_the_chains_outside_the_band() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut taken_over = 0u64;
    let mut band_differ = 0u64;
    for wl in [16usize, 18, 21] {
        for w in windows(wl) {
            let ncbi = w as i32 + 1;
            let (differ, other) = against_chains(wl, w, ncbi, &mut rng);
            assert_eq!(differ, 0, "NCBI's window: wl {wl} w {w}");
            taken_over += other;
            let losat = diag_hash_insert_window(w, 0, wl);
            let (differ, _) = against_chains(wl, w, losat, &mut rng);
            if 2 * w + 2 >= wl {
                assert_eq!(differ, 0, "LOSAT's window: wl {wl} w {w}");
            } else {
                band_differ += differ;
            }
        }
    }
    assert!(taken_over > 0, "no cell was taken over");
    assert!(band_differ > 0, "the band showed no difference");
}

// The masks of the 128-bit compares (`lanes_mask` on x86_64 and wasm32 with simd128) are the lane
// by lane masks: lanes and values drawn from a few values around each other and the ends of `i32`
// (`NO_DIAG` and `NO_LEVEL` are `i32::MIN`), every lane equal, below and above.
#[test]
fn the_lane_masks_are_the_lane_by_lane_masks() {
    let mut rng = Rng(0x3A5C_0FF1);
    let values = [
        i32::MIN,
        i32::MIN + 1,
        -2,
        -1,
        0,
        1,
        2,
        3,
        i32::MAX - 1,
        i32::MAX,
    ];
    let pick = |rng: &mut Rng| values[(rng.next() % values.len() as u64) as usize];
    let (mut some_eq, mut some_lt) = (0, 0);
    for _ in 0..200_000 {
        let lanes: [i32; 8] = std::array::from_fn(|_| pick(&mut rng));
        let v = pick(&mut rng);
        let (eq, lt) = (eq_mask(&lanes, v), lt_mask(&lanes, v));
        assert_eq!(eq, scalar_mask::<false>(&lanes, v), "{lanes:?} == {v}");
        assert_eq!(lt, scalar_mask::<true>(&lanes, v), "{lanes:?} < {v}");
        some_eq += usize::from(eq != 0 && eq != 0xff);
        some_lt += usize::from(lt != 0 && lt != 0xff);
    }
    assert!(some_eq > 1000 && some_lt > 1000);
    for v in values {
        assert_eq!(eq_mask(&[v; 8], v), 0xff);
        assert_eq!(lt_mask(&[v; 8], v), 0);
    }
}
