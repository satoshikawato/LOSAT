use super::super::{advance_diag_table_offset, diag_hash_insert_window, DiagStruct};
use super::*;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

fn bucket_of(diag: i32) -> usize {
    ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize
}

/// A table without bounds, or with bounds from the start.
fn table(window: i32, levels: bool) -> DiagHashTable {
    let mut h = DiagHashTable::new(window);
    h.max_level = Vec::new();
    h.x_last = Vec::new();
    if levels {
        h.x_lv_start();
    }
    h
}

type HashState = (Vec<u32>, Vec<(i32, i32, bool, i32, u32)>, u32, u32, i32);

fn hash_state(h: &DiagHashTable) -> HashState {
    (
        h.backbone.clone(),
        h.chain[..h.occupancy as usize]
            .iter()
            .map(|c| (c.diag, c.level, c.hit_saved, c.hit_len, c.next))
            .collect(),
        h.occupancy,
        h.capacity,
        h.offset,
    )
}

/// (i): `max_level` bounds every reachable cell's level; a bucket's walk from its head.
fn check_bounds(h: &DiagHashTable, step: usize) {
    for bucket in 0..h.num_buckets {
        let mut index = h.backbone[bucket];
        while index != 0 {
            let cell = &h.chain[index as usize];
            assert_eq!(bucket_of(cell.diag), bucket, "step {step}");
            assert!(
                cell.level <= h.max_level[bucket],
                "step {step}: bucket {bucket} cell {index} level {} above its bound {}",
                cell.level,
                h.max_level[bucket]
            );
            index = cell.next;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-450
// Random inserts (diagonals d + 512 k share a bucket, so cells are matched, taken over
// and appended, the chain grows past 256 cells, levels go up and down, some negative),
// offset advances between "subjects" and clears at INT4_MAX / 4: after every step the
// bounds hold for every bucket, the table is the table without bounds, and after a clear
// every bound is `i32::MIN`.
#[test]
fn the_bounds_hold_under_inserts_advances_and_clears() {
    let mut rng = Lcg(0xD1A6_17C0);
    for case in 0..30 {
        let window = [0, 16, 40][case % 3];
        let mut on = table(window, true);
        let mut off = table(window, false);
        let base = rng.below(2000) as i32 - 1000;
        // 24 diagonals in three buckets, and (in every other case) 600 more, so that the
        // chain grows past its 256 cells between clears
        let spread = if case % 2 == 0 { 600 } else { 16 };
        let diags: Vec<i32> = (0..24 + spread)
            .map(|k| {
                if k < 24 {
                    base + 512 * (k % 8) + (k / 8) * 7
                } else {
                    rng.below(1_000_000) as i32 - 500_000
                }
            })
            .collect();
        let mut grown = false;
        let (mut hl_on, mut hl_off) = (Vec::<DiagStruct>::new(), Vec::<DiagStruct>::new());
        let (mut ln_on, mut ln_off) = (Vec::<u8>::new(), Vec::<u8>::new());
        let (mut dt_on, mut dt_off) = (window, window);
        let mut clears = 0;
        for step in 0..4000 {
            match rng.below(100) {
                0..=1 => {
                    // the next subject; sometimes past INT4_MAX / 4 (the clear)
                    let len = if rng.below(4) == 0 {
                        (i32::MAX / 4) as usize
                    } else {
                        rng.below(5000) as usize
                    };
                    let cleared = on.offset >= i32::MAX / 4;
                    advance_diag_table_offset(
                        &mut dt_on, window, len, false, &mut hl_on, &mut ln_on, &mut on,
                    );
                    advance_diag_table_offset(
                        &mut dt_off,
                        window,
                        len,
                        false,
                        &mut hl_off,
                        &mut ln_off,
                        &mut off,
                    );
                    if cleared {
                        clears += 1;
                        assert_eq!(on.occupancy, 1);
                        assert!(on.max_level.iter().all(|&m| m == i32::MIN));
                    }
                }
                _ => {
                    let diag = diags[rng.below(diags.len() as u64) as usize];
                    let s_off = on.offset.saturating_add(rng.below(3000) as i32);
                    let level = s_off + rng.below(400) as i32 - 100;
                    let hit_len = rng.below(30) as i32;
                    let saved = rng.below(2) == 0;
                    let stale = rng.below(60) as i32 - 10;
                    on.insert_memo(diag, level, hit_len, saved, s_off, stale);
                    off.insert_memo(diag, level, hit_len, saved, s_off, stale);
                }
            }
            grown |= on.capacity > 256;
            check_bounds(&on, step);
            assert_eq!(hash_state(&on), hash_state(&off), "case {case} step {step}");
            assert!(off.max_level.is_empty());
        }
        assert!(clears > 0, "case {case}: no clear");
        assert!(grown || spread < 600, "case {case}: the chain never grew");
    }
}

/// The decisions of `s_BlastnDiagHashExtendInitialHit` that read the retrieve (module
/// comment, (ii)): rejected, and the two-hit branch taken.
fn decisions(
    found: Option<(i32, i32, bool)>,
    s_off_pos: usize,
    two_hits: bool,
    word_length: usize,
    window_size: usize,
) -> (bool, bool) {
    let (last_hit, _hit_len, hit_saved) = found.unwrap_or((0, 0, false));
    let s_off_pos_i32 = s_off_pos as i32;
    let rejected = s_off_pos_i32 < last_hit;
    let s_end_pos_i32 = (s_off_pos + word_length) as i32;
    let branch = two_hits && (hit_saved || s_end_pos_i32 > last_hit + window_size as i32);
    (rejected, branch)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:833-841
// Random tables filled as a fold fills them (seed positions rising from the hash offset,
// levels at or past the seed), probed at random diagonals and positions with the bounds
// of one hit (megablast and blastn word lengths, window 0) and of two hits
// (dc-megablast's window 40 and template lengths, `scan_range` 0 and 2): the decisions
// with `retrieve_below` are those with `retrieve`, a walk that is not skipped returns the
// same cell, and walks are skipped.
#[test]
fn retrieve_below_makes_the_decisions_of_retrieve() {
    let mut rng = Lcg(0x5EED_D1A6);
    let configs: [(bool, usize, usize, usize); 6] = [
        (false, 28, 0, 0),
        (false, 11, 0, 0),
        (true, 18, 40, 0),
        (true, 21, 40, 0),
        (true, 16, 40, 2),
        (true, 11, 40, 2),
    ];
    for (case, &(two_hits, word_length, window_size, scan_range)) in configs.iter().enumerate() {
        let window = window_size as i32;
        let diag_hash_window = diag_hash_insert_window(window_size, scan_range, word_length);
        let seeds = LvSeeds::new(two_hits, word_length, window_size);
        let mut h = table(window, true);
        h.offset = window + rng.below(100_000) as i32;
        let base = h.offset as usize;
        let diags: Vec<i32> = (0..3000)
            .map(|_| rng.below(20_000) as i32 - 10_000)
            .collect();
        let (mut skipped, mut probes, mut found_below) = (0u64, 0u64, 0u64);
        let mut s = base;
        for _ in 0..20_000 {
            s += rng.below(3) as usize;
            // a probe at the current position (as the fold's seeds come in order)
            let diag = diags[rng.below(diags.len() as u64) as usize];
            let bound = seeds.bound(s);
            let full = h.retrieve(diag);
            let before = test_skips();
            let fast = h.retrieve_below::<true>(diag, &seeds, s);
            let skip = test_skips() > before;
            probes += 1;
            if skip {
                skipped += 1;
                assert_eq!(fast, None);
                if let Some((level, _, _)) = full {
                    assert!(level < bound);
                    found_below += 1;
                }
            } else {
                assert_eq!(fast, full, "case {case}: a walk returned another cell");
            }
            assert_eq!(
                decisions(fast, s, two_hits, word_length, window_size),
                decisions(full, s, two_hits, word_length, window_size),
                "case {case} diag {diag} s_off_pos {s} bound {bound}: {full:?}"
            );
            // then an insert as the fold makes one (first hit or extension end)
            let reach = if rng.below(4) == 0 { 300 } else { 3 };
            let level = (s + word_length + rng.below(reach) as usize) as i32;
            let saved = rng.below(3) == 0;
            let hit_len = if saved { 0 } else { word_length as i32 };
            h.insert_memo(diag, level, hit_len, saved, s as i32, diag_hash_window);
        }
        assert!(
            skipped * 5 > probes && found_below > 0,
            "case {case}: too few skips ({skipped} of {probes}, {found_below} with a cell)"
        );
    }
}

/// LOSAT_X_DIAGLV2 (iv) and (v): in every bucket the first cell of the latest insert's
/// diagonal holds the latest insert's values and every other cell lies at or below `other`.
fn check_memo(h: &DiagHashTable, at: &str) {
    for bucket in 0..h.num_buckets {
        let last = h.x_last[bucket];
        let mut first = true;
        let mut index = h.backbone[bucket];
        while index != 0 {
            let cell = &h.chain[index as usize];
            if first && cell.diag == last.diag {
                first = false;
                assert_eq!(
                    (cell.level, cell.hit_len, cell.hit_saved),
                    (last.level, last.hit_len, last.hit_saved),
                    "{at}: bucket {bucket} cell {index}"
                );
            } else {
                assert!(
                    cell.level <= last.other,
                    "{at}: bucket {bucket} cell {index} level {} above {last:?}",
                    cell.level
                );
            }
            index = cell.next;
        }
        assert!(
            !first || h.backbone[bucket] == 0 || last.diag == LvLast::NO_DIAG,
            "{at}: bucket {bucket} has no cell of the memo's diagonal ({last:?})"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:368-381,396-450
// LOSAT_X_DIAGLV2: tables filled as the dc-megablast fold fills them (seeds in rising
// order, an insert after every probe; diagonals d + 512 k share a bucket, so cells are
// matched, taken over and appended), with offset advances between "subjects" and clears at
// INT4_MAX / 4, at the template lengths of dc-megablast and `scan_range` 0 and 2, with runs
// of the reference fold (no memo read or kept: `retrieve_below::<false>`, `insert_plain`)
// between strips stages (each makes the memo again, `x_lv2_reset`): the memo holds (iv) and
// (v) in every stage, the table is the table without a memo, a memo answer for the latest
// insert's diagonal is the walk's cell, every answer makes the decisions of the walk, and the
// memo answers (also for the latest insert's diagonal) where the bound does not skip.
#[test]
fn the_memo_makes_the_decisions_of_retrieve() {
    let mut rng = Lcg(0x1A57_D1A6);
    let configs: [(usize, usize); 4] = [(16, 0), (18, 0), (21, 0), (11, 2)];
    for (case, &(word_length, scan_range)) in configs.iter().enumerate() {
        let window_size = 40usize;
        let window = window_size as i32;
        let diag_hash_window = diag_hash_insert_window(window_size, scan_range, word_length);
        let seeds = LvSeeds::new(true, word_length, window_size);
        let mut h = table(window, true);
        h.x_lv2_start();
        let mut plain = table(window, false);
        let base = rng.below(2000) as i32 - 1000;
        let diags: Vec<i32> = (0..2400)
            .map(|k| {
                if k < 48 {
                    base + 512 * (k % 16) + (k / 16) * 7
                } else {
                    rng.below(20_000) as i32 - 10_000
                }
            })
            .collect();
        let (mut hl_on, mut hl_off) = (Vec::<DiagStruct>::new(), Vec::<DiagStruct>::new());
        let (mut ln_on, mut ln_off) = (Vec::<u8>::new(), Vec::<u8>::new());
        let (mut dt_on, mut dt_off) = (window, window);
        let (mut memos, mut memo_found, mut clears) = (0u64, 0u64, 0u64);
        // steps left in a run of the reference fold, and the runs and stage starts
        let (mut reference_left, mut runs, mut resets) = (0usize, 0u64, 0u64);
        let mut s = h.offset as usize;
        for step in 0..30_000 {
            if reference_left == 0 && rng.below(700) == 0 {
                reference_left = 1 + rng.below(400) as usize;
                runs += 1;
            } else if reference_left == 0 && rng.below(1000) == 0 {
                // a stage after a stage
                h.x_lv2_reset();
                resets += 1;
                check_memo(&h, &format!("case {case} step {step} reset"));
            }
            if rng.below(1000) == 0 {
                // the next subject; sometimes past INT4_MAX / 4 (the clear)
                let len = if rng.below(3) == 0 {
                    (i32::MAX / 4) as usize
                } else {
                    rng.below(5000) as usize
                };
                let cleared = h.offset >= i32::MAX / 4;
                advance_diag_table_offset(
                    &mut dt_on, window, len, false, &mut hl_on, &mut ln_on, &mut h,
                );
                advance_diag_table_offset(
                    &mut dt_off,
                    window,
                    len,
                    false,
                    &mut hl_off,
                    &mut ln_off,
                    &mut plain,
                );
                if cleared {
                    clears += 1;
                    assert!(h.x_last.iter().all(|&l| l == LvLast::EMPTY));
                }
                s = h.offset as usize;
            }
            s += rng.below(3) as usize;
            let diag = diags[rng.below(diags.len() as u64) as usize];
            let full = plain.retrieve(diag);
            assert_eq!(h.retrieve(diag), full, "case {case} step {step}");
            let before = test_memos();
            let fast = if reference_left > 0 {
                h.retrieve_below::<false>(diag, &seeds, s)
            } else {
                h.retrieve_below::<true>(diag, &seeds, s)
            };
            if reference_left > 0 {
                assert_eq!(test_memos(), before, "case {case} step {step}: memo read");
            } else if test_memos() > before {
                memos += 1;
                if fast.is_some() {
                    memo_found += 1;
                    assert_eq!(fast, full, "case {case} step {step}: the memo's cell");
                }
            }
            assert_eq!(
                seeds.decisions(fast, s),
                seeds.decisions(full, s),
                "case {case} step {step} diag {diag} s_off_pos {s}: {fast:?} {full:?}"
            );
            let reach = if rng.below(4) == 0 { 300 } else { 3 };
            let level = (s + word_length + rng.below(reach) as usize) as i32;
            let saved = rng.below(3) == 0;
            let hit_len = if saved { 0 } else { word_length as i32 };
            if reference_left > 0 {
                h.insert_plain(diag, level, hit_len, saved, s as i32, diag_hash_window);
                reference_left -= 1;
                if reference_left == 0 {
                    // the next strips stage
                    h.x_lv2_reset();
                    resets += 1;
                    check_memo(&h, &format!("case {case} step {step} reset"));
                }
            } else {
                h.insert_memo(diag, level, hit_len, saved, s as i32, diag_hash_window);
            }
            plain.insert_plain(diag, level, hit_len, saved, s as i32, diag_hash_window);
            if step % 97 == 0 {
                if reference_left == 0 {
                    check_memo(&h, &format!("case {case} step {step}"));
                }
                assert_eq!(
                    hash_state(&h),
                    hash_state(&plain),
                    "case {case} step {step}"
                );
            }
        }
        assert!(clears > 0, "case {case}: no clear");
        assert!(
            runs > 10 && resets > runs,
            "case {case}: {runs} runs of the reference fold, {resets} resets"
        );
        assert!(
            memos > 1000 && memo_found > 100,
            "case {case}: too few memo answers ({memos}, {memo_found} for the latest insert)"
        );
    }
}

// LOSAT_X_DIAGLV2, debug builds: a read of the memo, and an insert that keeps it, after an
// insert that did not keep it (the reference fold's `insert_plain`) and before a stage made
// it again (`x_lv2_reset`) panic, each with its own message; after the reset both work.
#[cfg(debug_assertions)]
#[test]
#[cfg_attr(panic = "abort", ignore = "catches a panic: the build aborts on panic")]
fn a_stale_memo_panics_in_debug_builds() {
    use std::panic::{catch_unwind, AssertUnwindSafe};
    let seeds = LvSeeds::new(true, 18, 40);
    for (read, expected) in [
        (true, "a retrieve reads a stale memo"),
        (false, "an insert keeps a stale memo"),
    ] {
        let mut h = table(40, true);
        h.x_lv2_start();
        h.insert_memo(7, 500, 0, true, 480, 41);
        h.insert_plain(7 + 512, 490, 0, true, 470, 41);
        h.x_lv2_reset();
        let fresh = catch_unwind(AssertUnwindSafe(|| {
            let _ = h.retrieve_below::<true>(7, &seeds, 480);
            h.insert_memo(7 + 1024, 505, 0, true, 485, 41);
        }));
        assert!(
            fresh.is_ok(),
            "read {read}: the memo after the reset panicked"
        );
        h.insert_plain(7, 520, 0, true, 500, 41);
        let stale = catch_unwind(AssertUnwindSafe(|| {
            if read {
                let _ = h.retrieve_below::<true>(7, &seeds, 480);
            } else {
                h.insert_memo(7 + 1024, 530, 0, true, 510, 41);
            }
        }));
        let payload = stale.expect_err("a stale memo did not panic");
        let message = payload
            .downcast_ref::<String>()
            .cloned()
            .or_else(|| payload.downcast_ref::<&str>().map(|s| s.to_string()))
            .unwrap_or_default();
        assert!(message.contains(expected), "read {read}: {message}");
    }
}

#[test]
fn the_bound_is_the_smallest_limit_and_none_where_the_argument_fails() {
    // one hit: s_off_pos + 1
    assert_eq!(level_bound(1000, false, 28, 0), 1001);
    assert_eq!(level_bound(1000, false, 11, 0), 1001);
    // two hits (window 40): the window limit s_end_pos - window_size lies below
    // s_off_pos + 1 (template 18: 1000 + 18 - 40; word 11: 1000 + 11 - 40)
    assert_eq!(level_bound(1000, true, 18, 40), 978);
    assert_eq!(level_bound(1000, true, 11, 40), 971);
    // a word longer than the window (not reachable): s_off_pos + 1
    assert_eq!(level_bound(1000, true, 50, 40), 1001);
    // s_end_pos <= window_size: no bound
    assert_eq!(level_bound(10, true, 18, 40), i32::MIN);
    assert_eq!(level_bound(22, true, 18, 40), i32::MIN);
    assert_eq!(level_bound(23, true, 18, 40), 1);
    // past i32::MAX: no bound
    assert_eq!(level_bound(i32::MAX as usize, false, 28, 0), i32::MIN);
    assert_eq!(level_bound(usize::MAX, false, 28, 0), i32::MIN);
    assert_eq!(level_bound(i32::MAX as usize - 27, false, 28, 0), i32::MIN);
    assert_eq!(
        level_bound(i32::MAX as usize - 28, false, 28, 0),
        i32::MAX - 27
    );
    assert_eq!(
        level_bound(i32::MAX as usize - 28, true, 28, 40),
        i32::MAX - 40
    );
    // an empty bucket (`i32::MIN`) is never below `i32::MIN`
    let h = table(40, true);
    assert_eq!(
        h.retrieve_below::<true>(5, &LvSeeds::new(true, 18, 40), 10),
        None
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:839-841,939-941
// The bound at its edge, for every setting the command line reaches: two hits with window
// 40 and templates 16, 18 and 21 (`scan_range` 0), and one hit with words 7, 11 and 28. One
// cell on the probed diagonal at `bound - 1` (with `hit_saved` 0 and 1) is skipped
// and decides as no cell; at `bound` it is walked, and with `hit_saved` 0
// (with one hit, either) it decides otherwise than no cell, so a bound one higher would
// change a decision. The insert the fold makes at that seed then leaves the same table with
// the bounds as without, and the bounds hold.
#[test]
fn the_bound_is_tight_at_every_reachable_setting() {
    let configs: [(bool, usize, usize); 6] = [
        (true, 16, 40),
        (true, 18, 40),
        (true, 21, 40),
        (false, 7, 0),
        (false, 11, 0),
        (false, 28, 0),
    ];
    let effective = |(rejected, branch): (bool, bool)| (rejected, !rejected && branch);
    for &(two_hits, word_length, window_size) in &configs {
        let seeds = LvSeeds::new(two_hits, word_length, window_size);
        let diag_hash_window = diag_hash_insert_window(window_size, 0, word_length);
        let window = window_size as i32;
        for s_off_pos in [window_size, window_size + 1, 5_000, 1 << 20] {
            let bound = seeds.bound(s_off_pos);
            let expected = if two_hits {
                s_off_pos + word_length - window_size
            } else {
                s_off_pos + 1
            };
            assert_eq!(bound as usize, expected);
            assert_eq!(
                bound,
                level_bound(s_off_pos, two_hits, word_length, window_size)
            );
            for (level, hit_saved) in [
                (bound - 1, false),
                (bound - 1, true),
                (bound, false),
                (bound, true),
            ] {
                let at = format!(
                    "two_hits {two_hits} word {word_length} s_off_pos {s_off_pos} bound {bound} level {level} hit_saved {hit_saved}"
                );
                let diag = 4321;
                let hit_len = if hit_saved { 0 } else { word_length as i32 };
                // one cell, written by an earlier seed of this diagonal
                let earlier = (level as usize).saturating_sub(word_length) as i32;
                let mut on = table(window, true);
                let mut off = table(window, false);
                on.insert_memo(diag, level, hit_len, hit_saved, earlier, diag_hash_window);
                off.insert_memo(diag, level, hit_len, hit_saved, earlier, diag_hash_window);
                let full = off.retrieve(diag);
                assert_eq!(full, Some((level, hit_len, hit_saved)), "{at}");
                let before = test_skips();
                let fast = on.retrieve_below::<true>(diag, &seeds, s_off_pos);
                let skipped = test_skips() > before;
                let walked = decisions(full, s_off_pos, two_hits, word_length, window_size);
                let none = decisions(None, s_off_pos, two_hits, word_length, window_size);
                // `LvSeeds::decisions` gives those of the test's
                assert_eq!(seeds.decisions(full, s_off_pos), effective(walked), "{at}");
                assert_eq!(seeds.decisions(None, s_off_pos), effective(none), "{at}");
                if level < bound {
                    assert!(skipped, "{at}: not skipped");
                    assert_eq!(fast, None, "{at}");
                    assert_eq!(walked, none, "{at}");
                } else {
                    assert!(!skipped, "{at}: skipped");
                    assert_eq!(fast, full, "{at}");
                    if !hit_saved || !two_hits {
                        assert_ne!(effective(walked), effective(none), "{at}: not tight");
                    }
                }
                // the insert at the end of this seed (unless it was rejected)
                if !walked.0 {
                    let s_end_pos = (s_off_pos + word_length) as i32;
                    let ready = walked.1;
                    let len = if ready { 0 } else { word_length as i32 };
                    for h in [&mut on, &mut off] {
                        h.insert_memo(
                            diag,
                            s_end_pos,
                            len,
                            ready,
                            s_off_pos as i32,
                            diag_hash_window,
                        );
                    }
                }
                check_bounds(&on, 0);
                assert_eq!(hash_state(&on), hash_state(&off), "{at}");
            }
        }
    }
}
