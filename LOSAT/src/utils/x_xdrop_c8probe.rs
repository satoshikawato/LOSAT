//! EXPERIMENT (line C, s7 / C8; branch exp/c8-align8-grid-probe, not for merging). No NCBI
//! counterpart: a measuring harness and a probe kernel, compiled into the test build only (the
//! library and the LOSAT binary do not change).
//!
//! `grid8` is `align8::<false, true, false>` (score-only, the `EM` row end) with the blocks of a
//! row on a fixed grid of columns (multiples of 32 lanes, so 32-byte aligned in the 64-byte
//! aligned rows) instead of at `first + 32k` / `first + 32k + 16k'`:
//!
//! - every block is 32 lanes; block 0 starts at `first & !31`, and the lanes left of `first` (the
//!   "dead" lanes) are computed like the others. They read only lanes the row before wrote as
//!   dropped (H 0, E 0) or the masked lane 0, so their H is at most `max_m < thr` and they are
//!   dropped (stored 0); their row gap contributions are below the threshold (the argument at the
//!   head of `x_xdrop_avx2.rs`), so no surviving cell, drop decision, best cell or band changes.
//! - the H of the row before is loaded at the block's own column (`hr + b`, the place and width of
//!   one store of the row before) and moved up one lane in registers with the last lane of the
//!   previous block's load (zero for block 0, which is the masked lane 0 of `align8`). The E load
//!   is at the same place as before (`e + b`). So every load of a row reads exactly one store of
//!   the row before and can be forwarded; in `align8` the next row's E (when `first` moves) and H
//!   (when it does not move by exactly one) span two stores (line R, R1), which ties a row to the
//!   end of the row before.
//! - past the band the row stops at a block start (multiple of 32) where the row gap entering the
//!   block is below the threshold (at most 31 lanes more than `align8`'s 16-lane steps; the lanes
//!   are dropped by the same argument).
//!
//! Harness `c8_loop` (ignored): `C8_DIR` (capture directory, `LOSAT_X_DPCAPTURE`), `C8_KERN`
//! (`dev` = `align8::<false, true, false>`, `grid` = `grid8`), `C8_REPS` (0: load, check and count
//! only), `C8_LIMIT`, `C8_OUT` (append the result line). Every handled 8-bit score-only problem is
//! checked once against the recorded result with both kernels, and counted (rows, cells, blocks of
//! `grid8`); then the chosen kernel runs `C8_REPS` times over all problems. Counters: `perf stat`
//! of `C8_REPS=K` minus `C8_REPS=0`, divided by K times the rows.

use super::super::tests::{random_problem, Rng};
use super::super::x_capture::load_dir;
use super::*;
use std::time::Instant;

fn var(k: &str) -> Option<String> {
    std::env::var(k).ok().filter(|v| !v.is_empty())
}

#[derive(Default, Clone, Copy)]
struct Cnt {
    rows: u64,
    blocks: u64,
    lanes_band: u64,
}

thread_local! {
    static COUNTS: std::cell::Cell<Cnt> = const { std::cell::Cell::new(Cnt { rows: 0, blocks: 0, lanes_band: 0 }) };
}

/// `block8::<false, true>` with the H of the row before loaded at `b` and moved up one lane with
/// `prev` (the previous block's load of the same row, zero for block 0).
#[inline(always)]
unsafe fn block8g(
    c: &C8,
    rs: &mut Row,
    prev: &mut W,
    hr: *const u8,
    hw: *mut u8,
    e: *mut u8,
    b: usize,
    s: W,
    lane_mask: Option<W>,
) -> W {
    let cur = ld(hr.add(b));
    let hp = _mm256_alignr_epi8::<15>(cur, _mm256_permute2x128_si256::<0x21>(*prev, cur));
    *prev = cur;
    let ep = ld(e.add(b));
    let d = _mm256_adds_epi8(hp, s);
    let h0 = _mm256_max_epi8(d, ep);
    let v = _mm256_subs_epu8(h0, c.goe);
    let fi = _mm256_subs_epu8(pmax8(_mm256_adds_epu8(up1_8(v), c.ramp), c.kb15), c.ramp);
    let f = _mm256_max_epu8(fi, _mm256_subs_epu8(rs.carry, c.ramp));
    let mut hh = _mm256_max_epi8(h0, f);
    if let Some(m) = lane_mask {
        hh = _mm256_and_si256(hh, m);
    }
    let hg = _mm256_subs_epu8(hh, c.goe);
    let out = last8(_mm256_max_epu8(_mm256_subs_epu8(fi, c.ge1), v), c.kb15);
    rs.carry = _mm256_max_epu8(_mm256_subs_epu8(rs.carry, c.ge32), out);
    let es = _mm256_subs_epu8(ep, c.ge1);
    let en = _mm256_max_epi8(hg, es);
    let gt = _mm256_cmpgt_epi8(hh, rs.bestv);
    let thr;
    if _mm256_movemask_epi8(gt) != 0 {
        let pfx = _mm256_max_epu8(pmax8(up1_8(hh), c.kb15), rs.bestv);
        thr = _mm256_subs_epu8(pfx, c.xv);
        let rm = _mm256_movemask_epi8(_mm256_cmpgt_epi8(hh, pfx)) as u32;
        rs.rec_b = b + (31 - rm.leading_zeros()) as usize;
        rs.bestv = last8(_mm256_max_epu8(pfx, hh), c.kb15);
        rs.thrv = _mm256_subs_epu8(rs.bestv, c.xv);
        rs.slow = true;
    } else {
        thr = rs.thrv;
    }
    let drop = _mm256_cmpgt_epi8(thr, hh);
    st(hw.add(b), _mm256_andnot_si256(drop, hh));
    st(e.add(b), _mm256_andnot_si256(drop, en));
    drop
}

/// `align8::<false, true, false>` with the blocks on the 32-lane grid (see the head of this file).
/// `CNT`: count rows, blocks and band lanes into `COUNTS` (the check pass; the timed pass uses false).
#[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
pub(super) unsafe fn grid8<const CNT: bool>(
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    max_m: i32,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    check_fence: bool,
    sc: &mut XdropScratch,
) -> Option<XdropResult> {
    sc.ops.clear();
    if len1 == 0 || len2 == 0 {
        return Some(XdropResult {
            a_offset: 0,
            b_offset: 0,
            score: 0,
            fence_hit: false,
            cells: 0,
        });
    }
    let goe = gap_open.checked_add(gap_extend)?;
    let mut x = x_drop;
    if x < goe {
        x = goe;
    }
    if gap_extend < 1 || gap_open < 0 || max_m < 0 {
        return None;
    }
    let bias = x as i64 + goe as i64 + max_m as i64 + 1;
    let limit_rebase = 127 - bias - max_m as i64;
    if limit_rebase < 0 {
        return None;
    }
    if gap_extend > 4 {
        return None;
    }
    let bias = bias as i32;
    let limit_rebase = limit_rebase as i32;
    let num_extra = (x / gap_extend + 3) as usize;
    let ncols = scores.ncols();
    sc.reset();

    sc.ensure(num_extra.min(len2 + 2) + 16);
    let mut b_size = 1usize;
    {
        let h = sc.h8a.as_mut_ptr().add(PAD);
        let e = sc.e8.as_mut_ptr().add(PAD);
        *h = bias as u8;
        *e = (bias - goe) as u8;
        let mut score = -goe;
        for i in 1..=len2 {
            if score < -x {
                break;
            }
            *h.add(i) = (bias + score) as u8;
            *e.add(i) = (bias + score - goe) as u8;
            score -= gap_extend;
            b_size = i + 1;
        }
        *h.add(b_size - 1) = 0;
    }
    sc.dirty8 = b_size;

    let sat = |v: i32| v.min(255) as u8 as i8;
    let mut ramp = [0u8; 32];
    for (i, r) in ramp.iter_mut().enumerate() {
        *r = (i as i32 * gap_extend) as u8;
    }
    let c = C8 {
        goe: _mm256_set1_epi8(sat(goe)),
        ge1: _mm256_set1_epi8(sat(gap_extend)),
        ge16: _mm256_set1_epi8(sat(16 * gap_extend)),
        ge32: _mm256_set1_epi8(sat(32 * gap_extend)),
        ramp: ld(ramp.as_ptr()),
        xv: _mm256_set1_epi8(sat(x)),
        kb15: _mm256_set1_epi8(15),
        c3: _mm256_setzero_si256(),
        c10: _mm256_setzero_si256(),
        c40: _mm256_setzero_si256(),
    };
    let zero = _mm256_setzero_si256();
    let k70 = _mm256_set1_epi8(0x70);
    let k10 = _mm256_set1_epi8(0x10);
    let static_tabs: Option<*const u8> = match scores {
        Scores::Tables { t, .. } => Some(t.data.as_ptr()),
        _ => None,
    };
    let small_alphabet = ncols <= 16;
    let mut cells: u64 = 0;

    let mut base: i32 = 0;
    let mut best_rel: i32 = 0;
    let mut a_off = 0usize;
    let mut b_off = 0usize;
    let mut first = 0usize;
    let mut fence_hit = false;
    let mut wend = [b_size, 0usize];
    let mut cur = 0usize;
    let mut cnt = Cnt::default();

    let win_at = sc.x_win_at;
    let mut len2 = len2;
    for a in 1..=len1 {
        if first >= win_at {
            let d = (first - 64) & !63;
            sc.x_shift(d, false);
            first -= d;
            b_size -= d;
            len2 -= d;
            wend[0] = wend[0].saturating_sub(d);
            wend[1] = wend[1].saturating_sub(d);
        }
        sc.ensure(b_size + num_extra + 64);
        if sc.sb_filled < b_size {
            fill_cols(s, sc, (b_size + 64).min(len2 + 1), ncols, check_fence);
        }
        if sc.bad_k < b_size {
            return None;
        }
        let mut limit = len2 + 1;
        let mut bend = b_size;
        let mut fence_row = false;
        if sc.fence_k < b_size {
            limit = sc.fence_k;
            bend = sc.fence_k;
            fence_row = true;
            if bend <= first {
                fence_hit = true;
                break;
            }
        }
        let qa = q.get(a)? as usize;
        if qa >= ncols {
            return None;
        }
        cells += (b_size - first) as u64;
        let t: *const u8 = if let Some(tabs) = static_tabs {
            tabs.add(qa * TAB)
        } else {
            if sc.tab_ok & (1u32 << qa) == 0 {
                if !build_row(scores, qa, &mut sc.tabs[qa * TAB..(qa + 1) * TAB]) {
                    return None;
                }
                sc.tab_ok |= 1u32 << qa;
            }
            sc.tabs.as_ptr().add(qa * TAB)
        };
        let (hr, hw) = if cur == 0 {
            (sc.h8a.as_mut_ptr().add(PAD), sc.h8b.as_mut_ptr().add(PAD))
        } else {
            (sc.h8b.as_mut_ptr().add(PAD), sc.h8a.as_mut_ptr().add(PAD))
        };
        let e = sc.e8.as_mut_ptr().add(PAD);
        let sb = sc.sb.as_ptr().add(PAD);
        let t0 = bcast128(t.add(64));
        let t1 = bcast128(t.add(80));
        let s0 = first & !31;

        if best_rel > limit_rebase {
            let dv = _mm256_set1_epi8(best_rel as u8 as i8);
            let mut b = s0;
            while b < b_size {
                st(hr.add(b), _mm256_subs_epu8(ld(hr.add(b)), dv));
                st(e.add(b), _mm256_subs_epu8(ld(e.add(b)), dv));
                b += 32;
            }
            base += best_rel;
            best_rel = 0;
        }

        let bestv0 = _mm256_set1_epi8((bias + best_rel) as u8 as i8);
        let mut rs = Row {
            carry: zero,
            bestv: bestv0,
            thrv: _mm256_subs_epu8(bestv0, c.xv),
            rec_b: 0,
            slow: false,
        };
        let mut prev = zero;
        let mut lo_surv = usize::MAX;
        let mut hi_surv = 0usize;
        let mut b0 = s0;
        loop {
            if b0 >= bend {
                if b0 >= limit {
                    break;
                }
                if (_mm256_cvtsi256_si32(rs.carry) & 0xFF) < (_mm256_cvtsi256_si32(rs.thrv) & 0xFF)
                {
                    break;
                }
            }
            let idx = ld(sb.add(b0));
            let sv = if small_alphabet {
                _mm256_shuffle_epi8(t0, idx)
            } else {
                lut32(t0, t1, _mm256_add_epi8(idx, k70), _mm256_sub_epi8(idx, k10))
            };
            let mask = if b0 + 32 <= limit {
                None
            } else {
                Some(ld(LANE_MASK8.as_ptr().add(32 - (limit - b0))))
            };
            let d0 = block8g(&c, &mut rs, &mut prev, hr, hw, e, b0, sv, mask);
            let live = !(_mm256_movemask_epi8(d0) as u32);
            if live != 0 {
                if lo_surv == usize::MAX {
                    lo_surv = b0 + live.trailing_zeros() as usize;
                }
                hi_surv = b0 + (31 - live.leading_zeros() as usize);
            }
            b0 += 32;
            if CNT {
                cnt.blocks += 1;
            }
        }
        if CNT {
            cnt.rows += 1;
            cnt.lanes_band += (b_size - first) as u64;
        }
        let proc_end = b0;
        sc.dirty8 = sc.dirty8.max(proc_end);
        let wi = 1 - cur;
        let mut b = proc_end;
        while b < wend[wi] {
            st(hw.add(b), zero);
            b += 32;
        }
        wend[wi] = proc_end;
        cur = wi;
        if rs.slow {
            b_off = rs.rec_b + sc.x_sh;
            a_off = a;
            best_rel = (_mm256_cvtsi256_si32(rs.bestv) & 0xFF) - bias;
        }
        if fence_row {
            fence_hit = true;
            break;
        }
        if lo_surv == usize::MAX {
            break;
        }
        first = lo_surv;
        b_size = hi_surv + 1;
        if b_size <= len2 {
            b_size += 1;
        }
    }
    if CNT {
        COUNTS.with(|k| {
            let mut v = k.get();
            v.rows += cnt.rows;
            v.blocks += cnt.blocks;
            v.lanes_band += cnt.lanes_band;
            k.set(v);
        });
    }

    let best = base + best_rel;
    finish::<false, false>(sc, a_off, b_off, best, fence_hit, cells)
}

/// `grid8` against `align8::<false, true, false>` on random problems (forward and reversed, fences,
/// band windows, both alphabets, the 8-bit cost range): equal results. `LOSAT_FUZZ_CASES` sets the
/// number of problems.
#[test]
fn grid8_matches_align8_on_random_problems() {
    if !available() {
        eprintln!("no AVX2: skipped");
        return;
    }
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3000);
    let mut rng = Rng(0xC8C8_0001_9E37_79B9);
    let (mut sa, mut sg) = (XdropScratch::new(), XdropScratch::new());
    let (mut handled, mut fences, mut moved, mut slow_rows) = (0usize, 0usize, 0usize, 0u64);
    for case in 0..cases {
        let win_at = [128, 320, usize::MAX][case % 3];
        (sa.x_win_at, sg.x_win_at) = (win_at, win_at);
        let n = if rng.below(2) == 0 { 16 } else { 28 };
        let long = rng.below(5) == 0;
        let len = 1 + rng.below(if long { 2500 } else { 300 }) as usize;
        let mut p = random_problem(&mut rng, n, 1, len);
        if rng.below(4) == 0 {
            // near-identical pairs: long bands, rebases, many best updates
            p.s = p.q.clone();
            for _ in 0..rng.below(len as u64 / 10 + 1) {
                let k = rng.below(p.s.len() as u64) as usize;
                p.s[k] = 1 + rng.below((n - 1) as u64) as u8;
            }
        }
        let (go, ge, x) = match rng.below(5) {
            0 => (11, 1, 38),
            1 => (11, 1, 15 + rng.below(30) as i32),
            2 => (5, 2, 33),
            3 => (2, 1, 10 + rng.below(40) as i32),
            _ => (
                rng.below(8) as i32,
                1 + rng.below(4) as i32,
                10 + rng.below(80) as i32,
            ),
        };
        let len1 = p.q.len();
        if p.s.is_empty() {
            continue;
        }
        let mut s = p.s.clone();
        let mut len2 = s.len();
        let mut sentinel = false;
        let cf = rng.below(3) == 0;
        if cf {
            match rng.below(3) {
                0 => {
                    s.push(FENCE_SENTRY);
                    sentinel = true;
                }
                1 if len2 > 3 => {
                    let k = 1 + rng.below(len2 as u64 - 1) as usize;
                    s[k] = FENCE_SENTRY;
                }
                _ => {
                    len2 = 1 + rng.below(len2 as u64) as usize;
                    sentinel = rng.below(2) == 0;
                }
            }
        }
        let scores = Scores::Flat { data: &p.m, n };
        let Some(m) = fits_8bit(&scores, go, ge, x) else {
            continue;
        };
        let rq: Vec<u8> = p.q.iter().rev().copied().collect();
        let rs: Vec<u8> = s[..len2].iter().rev().copied().collect();
        for reverse in [false, true] {
            if reverse && sentinel {
                continue;
            }
            let (rows, cols) = if reverse {
                (
                    RowSeq::Bytes {
                        data: &rq,
                        base: len1 as isize,
                        step: -1,
                    },
                    ColSeq {
                        data: &rs,
                        base: len2 as isize - 1,
                        step: -1,
                        zero_from: len2,
                    },
                )
            } else {
                (
                    RowSeq::Bytes {
                        data: &p.q,
                        base: -1,
                        step: 1,
                    },
                    ColSeq {
                        data: &s,
                        base: 0,
                        step: 1,
                        zero_from: if sentinel { usize::MAX } else { len2 },
                    },
                )
            };
            // SAFETY: AVX2 detected above.
            let (ra, rg) = unsafe {
                (
                    align8::<false, true, false>(
                        &rows, &cols, &scores, m, len1, len2, go, ge, x, cf, &mut sa,
                    ),
                    grid8::<true>(&rows, &cols, &scores, m, len1, len2, go, ge, x, cf, &mut sg),
                )
            };
            assert_eq!(
                ra, rg,
                "case {case} rev={reverse} n={n} go={go} ge={ge} x={x} len1={len1} len2={len2} cf={cf} sentinel={sentinel} win_at={win_at}"
            );
            if let Some(r) = ra {
                handled += 1;
                // score-only results report no fence flag; count the calls with fence checking
                fences += usize::from(cf);
                moved += usize::from(sg.x_sh > 0);
                slow_rows += r.cells;
            }
        }
    }
    eprintln!("grid8 fuzz: cases={cases} handled={handled} check_fence={fences} moved={moved} cells={slow_rows}");
    if cases >= 1000 {
        assert!(handled > cases / 2 && fences > 0 && moved > 0);
    }
}

struct Prob<'a> {
    q: RowSeq<'a>,
    s: ColSeq<'a>,
    scores: Scores<'a>,
    m: i32,
    len1: usize,
    len2: usize,
    go: i32,
    ge: i32,
    x: i32,
    cf: bool,
    want: Option<XdropResult>,
}

#[inline(never)]
fn call(p: &Prob<'_>, grid: bool, sc: &mut XdropScratch) -> Option<XdropResult> {
    // SAFETY: `available()` is checked before any call.
    unsafe {
        if grid {
            grid8::<false>(
                &p.q, &p.s, &p.scores, p.m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
            )
        } else {
            align8::<false, true, false>(
                &p.q, &p.s, &p.scores, p.m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
            )
        }
    }
}

#[test]
#[ignore = "needs C8_DIR=<capture directory>"]
fn c8_loop() {
    let Some(dir) = var("C8_DIR") else {
        eprintln!("C8_DIR is not set");
        return;
    };
    assert!(available(), "the harness needs AVX2");
    let reps: usize = var("C8_REPS").and_then(|v| v.parse().ok()).unwrap_or(0);
    let limit: usize = var("C8_LIMIT")
        .and_then(|v| v.parse().ok())
        .unwrap_or(usize::MAX);
    let kern = var("C8_KERN").unwrap_or_else(|| "dev".into());
    let grid = match kern.as_str() {
        "dev" => false,
        "grid" => true,
        k => panic!("C8_KERN must be dev or grid, not {k}"),
    };
    let t_load = Instant::now();
    let probs = load_dir(std::path::Path::new(&dir)).expect("read the capture files");
    let mut list: Vec<Prob<'_>> = Vec::new();
    for p in &probs {
        if list.len() >= limit {
            break;
        }
        if p.tb || p.result.is_none() {
            continue;
        }
        let Some(m) = fits_8bit(&p.scores(), p.gap_open, p.gap_extend, p.x_drop) else {
            continue;
        };
        list.push(Prob {
            q: p.row_seq(),
            s: p.col_seq(),
            scores: p.scores(),
            m,
            len1: p.len1,
            len2: p.len2,
            go: p.gap_open,
            ge: p.gap_extend,
            x: p.x_drop,
            cf: p.check_fence,
            want: p.result,
        });
    }
    let load_s = t_load.elapsed().as_secs_f64();
    // check (both kernels) and count (grid8::<true>), once
    let mut sc = XdropScratch::new();
    let (mut ok, mut declined) = (0usize, 0usize);
    COUNTS.with(|k| k.set(Cnt::default()));
    for p in &list {
        let rd = call(p, false, &mut sc);
        // SAFETY: AVX2 checked above.
        let rg = unsafe {
            grid8::<true>(
                &p.q, &p.s, &p.scores, p.m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, &mut sc,
            )
        };
        assert_eq!(rd, rg, "grid8 vs align8");
        match rd {
            Some(r) => {
                assert_eq!(Some(r), p.want, "align8 vs recorded");
                ok += 1;
            }
            None => declined += 1,
        }
    }
    let cnt = COUNTS.with(|k| k.get());
    let t0 = Instant::now();
    for _ in 0..reps {
        for p in &list {
            std::hint::black_box(call(p, grid, &mut sc));
        }
    }
    let loop_s = t0.elapsed().as_secs_f64();
    let rows_total = cnt.rows as f64 * reps as f64;
    let line = format!(
        "C8 kern={kern} dir={dir} problems={} handled={ok} declined={declined} reps={reps} \
         load_s={load_s:.3} loop_s={loop_s:.6} ns_per_row={:.3} rows={} band_lanes_per_row={:.2} \
         grid_blocks_per_row={:.3}\n",
        list.len(),
        if rows_total > 0.0 {
            loop_s * 1e9 / rows_total
        } else {
            0.0
        },
        cnt.rows,
        cnt.lanes_band as f64 / cnt.rows.max(1) as f64,
        cnt.blocks as f64 / cnt.rows.max(1) as f64,
    );
    eprint!("{line}");
    if let Some(out) = var("C8_OUT") {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(out)
            .expect("open C8_OUT");
        f.write_all(line.as_bytes()).expect("write C8_OUT");
    }
}
