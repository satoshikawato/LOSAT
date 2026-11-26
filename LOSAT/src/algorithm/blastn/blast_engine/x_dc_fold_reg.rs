//! EXPERIMENT (LOSAT_X_DCFOLDREG, a strict-set switch since 2026-11-26): the two-hit fold of
//! dc-megablast's strips stage on a map of LOSAT_X_DCLATEST (`x_dc_latest.rs`), with the map in
//! registers. No NCBI counterpart (the same statements in another layout of the code).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-841
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:882-895,923-941
//! (`s_BlastnDiagHashExtendInitialHit`, as `x_pair_seed_disc.rs` `fold_two_hit_on` makes it on a
//! map)
//!
//! `fold_two_hit_reg` makes every decision and insert of `fold_two_hit_on` on `Latest` (the fold
//! that `LOSAT_X_DCFOLDREG=0` runs), seed after seed in the same order, with the map's fields that
//! every seed reads or writes held as locals of the call (`latest::Hot`: the buckets, the shift of
//! the home bucket, `top` and `entries`) instead of being loaded and stored through the map for
//! every seed. What the scope of the maps fixes is taken once per call instead of per seed:
//! - Two hits and `scan_range == 0` (the scope of the strips stage, `in_scope`): the off-diagonal
//!   loop runs no iteration, so a seed that takes the two-hit branch with a single word
//!   (`word_type` 1) only inserts (`hit_ready` cleared), and `off_found` stays false.
//! - Where LOSAT_X_FOLDLEAN takes `type_of_word`'s early return (`TOW`: the word length equals the
//!   lookup width, always so for a discontiguous table), that return is a constant of the instance
//!   (`word_type` 1, nothing extended, the offsets unchanged), and the query context is looked up
//!   only by the extension (`extend`), which 24,191 of 188.7 M seeds reach in case Ad. With
//!   LOSAT_X_FOLDLEAN=0 (`TOW` false) every seed calls `type_of_word` with its context, as the fold
//!   of `=0` does.
//! - A seed whose diagonal the map does not hold (a miss: 95 % of the seeds of case Ad) reads level 0
//!   and no saved hit. Every position of a chunk on the maps has a level bound (`latest_bound`, and
//!   S never decreases within the stage, module comment (vi) of `x_dc_latest.rs`), so its S lies in
//!   `LvSeeds`' range: S is at least 0 (not rejected) and S plus the word length is above the window
//!   (the two-hit branch). With `TOW` such a seed inserts its S plus the word length, not saved, at
//!   once (`Hot::put_new`).
//! - The level bound of a put is S plus a constant of the call (`LvSeeds::bound_in`, no range test).

use super::super::super::x_diag_lv::LvSeeds;
use super::super::super::{
    extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi, mask_lookup_index,
    packed_kmer_at_seed_mask, type_of_word, BlastnTiming, UngappedHit,
};
use super::super::latest::{self, DiagLatest, Hot, MISS};
#[cfg(test)]
use super::super::latest::{Event, STASHED};
use super::super::{part, FoldCounts, Inputs, Seed};
use std::sync::atomic::Ordering;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-841
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:882-895,923-941
/// `fold_two_hit_on` on the map `map` (`Latest`), with the map in registers (module comment): the
/// same seeds rejected (`counts.skipped`), the same hits pushed to `ungapped_hits` (and with `PT`
/// their keys to `out`), the same inserts into the map in the same order. `TOW`: `type_of_word`
/// returns at once (LOSAT_X_FOLDLEAN with equal word lengths), else it is called for every seed.
/// Not inlined: one call per strip (or piece), and its own frame keeps the map's fields in
/// registers.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
pub(super) fn fold_two_hit_reg<const PT: bool, const TOW: bool>(
    x: &Inputs<'_>,
    seeds: &[Seed],
    map: &mut DiagLatest,
    ungapped_hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    scan_pause_ns: &mut u64,
    counts: &mut FoldCounts,
    out: &mut part::PartOut,
    pos: &[u32],
) {
    latest::count_reg_call();
    let two_stage = x.two_stage;
    let word_length = two_stage.word_length();
    // the scope of the maps (`Latest::new`, `latest_bound`)
    debug_assert!(x.window_size > 0 && x.scan_range == 0);
    debug_assert!(!TOW || word_length == two_stage.lut_word_length());
    let lv = LvSeeds::new(true, word_length, x.window_size);
    debug_assert!(
        lv.bound(x.diag_offset as usize) != i32::MIN
            && lv.bound(x.diag_offset as usize + x.s_len) != i32::MIN
    );
    // `s_off_pos`, `s_end_pos` and the window test in `i32`, as the reference fold casts them (the
    // same low 32 bits as its `usize` sums)
    let offset = x.diag_offset as u32;
    let wl = word_length as i32;
    let window = x.window_size as i32;
    counts.seeds += seeds.len() as u64;
    let mut skipped = 0u64;
    let mut hot = Hot::new(map);
    for (k, seed) in seeds.iter().enumerate() {
        // EXPERIMENT (LOSAT_X_DCLATEST): the map's bucket of a later seed (a hint only).
        if let Some(later) = seeds.get(k + latest::AHEAD) {
            hot.prefetch(later.s.wrapping_sub(later.q) as i32);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-836
        let diag = seed.s.wrapping_sub(seed.q) as i32;
        let s_pos = seed.s.wrapping_add(offset) as i32;
        let s_end_pos = s_pos.wrapping_add(wl);
        let b = hot.home(diag);
        let lanes = hot.lanes(b, diag);
        let (level, saved, at) = if lanes != 0 {
            let lane = lanes.trailing_zeros();
            let (level, saved) = hot.slot(b, lane);
            (level, saved, lane)
        } else if let Some((at, level, saved)) = hot.stashed(b, diag, s_pos) {
            (level, saved, at)
        } else if TOW {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:841-886,939-941
            // A miss (module comment): not rejected, the two-hit branch with a single word, the
            // insert.
            debug_assert!(s_pos >= 0 && s_end_pos > window);
            #[cfg(test)]
            {
                hot.note(Event::Read {
                    diag,
                    s: s_pos,
                    level: 0,
                    saved: false,
                    stash: false,
                });
                hot.note(Event::Put {
                    diag,
                    level: s_end_pos,
                    saved: false,
                    s: s_pos,
                });
            }
            hot.put_new(b, diag, s_end_pos, false, s_pos, lv.bound_in(s_pos));
            continue;
        } else {
            (0, false, MISS)
        };
        #[cfg(test)]
        hot.note(Event::Read {
            diag,
            s: s_pos,
            level,
            saved,
            stash: at >= STASHED,
        });

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:838-839
        if s_pos < level {
            skipped += 1;
            continue;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:841-851,887-895
        let two_hit = saved || s_end_pos > level.wrapping_add(window);
        let (word_type, extended, q_off, s_off) = if TOW {
            (1, 0, seed.q as usize, seed.s as usize)
        } else {
            word_of(x, seed, two_hit)
        };
        if word_type == 0 {
            continue;
        }
        let s_end_pos = s_end_pos.wrapping_add(extended as i32);

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:853-886,897,939-941
        // (`scan_range` is 0: no off-diagonal cell, `hit_ready` cleared for a single word)
        let (level, saved) = if two_hit && word_type == 1 {
            (s_end_pos, false)
        } else {
            extend::<PT>(
                x,
                seed,
                (q_off, s_off, seed.s as usize + word_length + extended),
                s_end_pos,
                ungapped_hits,
                timing,
                scan_pause_ns,
                out,
                (pos, k),
            )
        };
        #[cfg(test)]
        hot.note(Event::Put {
            diag,
            level,
            saved,
            s: s_pos,
        });
        hot.put(b, at, diag, level, saved, s_pos, lv.bound_in(s_pos));
    }
    counts.skipped += skipped;
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:841-851,887-895
/// `type_of_word` of `seed` as the fold of `=0` calls it without LOSAT_X_FOLDLEAN (`check_double`
/// in the two-hit branch): its query context's end and masks, the scan range's `s_range`, and
/// `s_MBLookup` as the masked-seed test.
#[inline(always)]
fn word_of(x: &Inputs<'_>, seed: &Seed, check_double: bool) -> (u8, usize, usize, usize) {
    let two_stage = x.two_stage;
    let lut_word_length = two_stage.lut_word_length();
    let s_len = x.s_len;
    let search_seq_packed = x.search_seq_packed;
    let ctx = &x.query_contexts[seed.context as usize];
    let q_context_end = ctx.query_offset as usize + ctx.len;
    let query_mask = ctx.masks.as_slice();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:41-69 (s_MBLookup)
    let mut is_seed_masked = |s_pos: usize, q_pos: usize| -> bool {
        if s_pos + lut_word_length > s_len {
            return true;
        }
        let kmer = mask_lookup_index(
            packed_kmer_at_seed_mask(search_seq_packed, s_pos, lut_word_length),
            lut_word_length,
        );
        if !two_stage.has_hits(kmer) {
            return true;
        }
        let q_off_1 = (q_pos + 1) as u32;
        !two_stage.contains_hit(kmer, q_off_1)
    };
    type_of_word(
        seed.q as usize,
        seed.s as usize,
        !query_mask.is_empty(),
        q_context_end,
        seed.s_range as usize,
        two_stage.word_length(),
        lut_word_length,
        check_double,
        &mut is_seed_masked,
    )
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:908-933,939-941
/// The ungapped extension of `seed` at `(q_off, s_off, s_end)` (after `type_of_word`), the cutoff
/// test (no off-diagonal cell was found: `scan_range` is 0) and, past it, the hit pushed (with `PT`
/// its key, `pos[k]` of `out.strip`); returns the level and the saved bit of the seed's insert.
/// Out of line: the seeds of the window's second hit only.
#[allow(clippy::too_many_arguments)]
#[cold]
#[inline(never)]
fn extend<const PT: bool>(
    x: &Inputs<'_>,
    seed: &Seed,
    (q_off, s_off, s_end): (usize, usize, usize),
    s_end_pos: i32,
    ungapped_hits: &mut Vec<UngappedHit>,
    timing: Option<&BlastnTiming>,
    scan_pause_ns: &mut u64,
    out: &mut part::PartOut,
    (pos, k): (&[u32], usize),
) -> (i32, bool) {
    let context_idx = seed.context as usize;
    let ctx = &x.query_contexts[context_idx];
    let cutoff_score = x.cutoff_scores[context_idx];
    let search_seq_packed = x.search_seq_packed;
    let s_len = x.s_len;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:908-920
    let x_dropoff = x.x_dropoff_scores[context_idx];
    let reduced_cutoff = x.reduced_cutoff_scores[context_idx];
    let ungapped_start = timing.is_some().then(std::time::Instant::now);
    let ungapped = if x.exact_ungapped {
        extend_hit_ungapped_exact_ncbi(
            x.encoded_query_concat_blastna,
            search_seq_packed,
            q_off,
            s_off,
            s_len,
            x_dropoff,
            x.score_matrix,
        )
    } else {
        extend_hit_ungapped_approx_ncbi(
            x.encoded_query_concat_blastna,
            search_seq_packed,
            q_off,
            s_off,
            s_end,
            s_len,
            x_dropoff,
            x.nucl_score_table,
            reduced_cutoff,
            x.score_matrix,
        )
    };
    if let Some(ungapped_start) = ungapped_start {
        let elapsed_ns = ungapped_start.elapsed().as_nanos() as u64;
        if PT {
            // a partition's own counters, added to the shared ones once per stage
            out.ungapped_ns += elapsed_ns;
            out.ungapped_calls += 1;
        } else if let Some(timing) = timing {
            timing.ungapped_ns.fetch_add(elapsed_ns, Ordering::Relaxed);
            timing.ungapped_calls.fetch_add(1, Ordering::Relaxed);
        }
        *scan_pause_ns = scan_pause_ns.saturating_add(elapsed_ns);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:923-932
    if ungapped.score < cutoff_score {
        return (s_end_pos, false);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:924-929
    let ungapped_s_end_pos = ungapped.s_start + ungapped.length + x.diag_offset as usize;
    super::super::super::x_push_ungapped(
        ungapped_hits,
        UngappedHit::saved(context_idx as u32, ctx, q_off, s_off, &ungapped),
    );
    if PT {
        out.hit_keys.push(out.key(pos, k));
    }
    (ungapped_s_end_pos as i32, true)
}
