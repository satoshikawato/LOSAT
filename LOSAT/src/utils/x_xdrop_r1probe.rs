//! EXPERIMENT (line R, item R1; branch exp/r1-xdrop-counters, not for merging). No NCBI
//! counterpart: a measuring harness only; it does not change any value NCBI computes, and it is
//! compiled into the test build only (the library and the LOSAT binary do not change).
//!
//! The 256-bit kernels of `x_xdrop_avx2.rs` alone in a loop over captured (`LOSAT_X_DPCAPTURE`)
//! or synthetic problems, for hardware counters: `perf stat` of a run with `R1_REPS=K` minus a run
//! with `R1_REPS=0`, divided by K times the rows. Each run also prints the shape the counts are
//! divided by: calls, rows, cells, lanes processed, the 32-lane and 16-lane blocks per row, the
//! advance of the band start, the rows and blocks that raise the best, and a model of how each
//! row's H and E loads overlap the stores of the row before (store forwarding).
//!
//! Variables: `R1_DIR` (capture directory; the synthetic problems take the matrix and gap costs of
//! its first 8-bit score-only problem), `R1_CAT` (so8, tb8, so16, tb16), `R1_SYNTH` (identical,
//! 80pct, 60pct, 45pct; unset = the captures), `R1_REPS` (default 0), `R1_LIMIT` (most problems),
//! `R1_OUT` (append the result line to this file), `R1_LAYOUT` (`spread`: before the timed loop, the
//! kernels' rows move into one leaked arena at offsets modulo 4096 at least 512 bytes apart, see
//! `spread`; unset = the rows where the scratch put them). The line also gives each row's address
//! modulo 4096 (`lay=`).

use super::super::tests::Rng;
use super::super::x_capture::load_dir;
use super::super::{fits_8bit, ColSeq, RowSeq, Scores, XdropResult, XdropScratch, FENCE_SENTRY};
use super::{align, available};
use std::time::Instant;

const MININT: i32 = i32::MIN / 2;

fn var(k: &str) -> Option<String> {
    std::env::var(k).ok().filter(|v| !v.is_empty())
}

/// One problem: the kernel's inputs and the plain residues and matrix of the shape model.
struct Prob<'a> {
    q: RowSeq<'a>,
    s: ColSeq<'a>,
    scores: Scores<'a>,
    len1: usize,
    len2: usize,
    go: i32,
    ge: i32,
    x: i32,
    tb: bool,
    cf: bool,
    rows: &'a [u8],
    cols: &'a [u8],
    mat: &'a [i32],
    n: usize,
    want: Option<XdropResult>,
}

/// One load or store of the model, in bytes of one buffer (0, 1: the two H rows; 2: E).
#[derive(Clone, Copy)]
struct Acc {
    row: usize,
    buf: u8,
    at: i64,
    len: i64,
}

#[derive(Default)]
struct Shape {
    calls: u64,
    rows: u64,
    cells: u64,
    lanes: u64,
    /// 32-lane steps (align8: block8; align16: two block16) and 16-lane steps.
    n32: u64,
    n16: u64,
    /// 16-lane steps at or past the band end (the row-gap extension).
    n16_ext: u64,
    best_rows: u64,
    best_blocks: u64,
    rebases: u64,
    /// advance of the band start: 0, 1, 2, 3..=8, more; and rows that end the call
    adv: [u64; 5],
    zero_st: u64,
    eclear_st: u64,
    /// loads of H / E: [no older store of the last row overlaps, one store with the same start and
    /// size, inside one store (other start or size), overlaps a store but no single store holds it]
    ld_h: [u64; 4],
    ld_e: [u64; 4],
    /// problems whose rows the model could not follow (fence, a decline)
    skipped: u64,
}

impl Shape {
    fn line(&self) -> String {
        let r = self.rows.max(1) as f64;
        format!(
            "calls={} rows={} cells={} lanes={} n32={} n16={} n16_ext={} best_rows={} \
             best_blocks={} rebases={} adv0={} adv1={} adv2={} adv3_8={} adv9p={} zero_st={} \
             eclear_st={} ldh_l1={} ldh_same={} ldh_in={} ldh_fail={} lde_l1={} lde_same={} \
             lde_in={} lde_fail={} skipped={} | per row: cells {:.2} lanes {:.2} n32 {:.3} \
             n16 {:.3} fail_h {:.3} fail_e {:.3} in_h {:.3} in_e {:.3}",
            self.calls,
            self.rows,
            self.cells,
            self.lanes,
            self.n32,
            self.n16,
            self.n16_ext,
            self.best_rows,
            self.best_blocks,
            self.rebases,
            self.adv[0],
            self.adv[1],
            self.adv[2],
            self.adv[3],
            self.adv[4],
            self.zero_st,
            self.eclear_st,
            self.ld_h[0],
            self.ld_h[1],
            self.ld_h[2],
            self.ld_h[3],
            self.ld_e[0],
            self.ld_e[1],
            self.ld_e[2],
            self.ld_e[3],
            self.skipped,
            self.cells as f64 / r,
            self.lanes as f64 / r,
            self.n32 as f64 / r,
            self.n16 as f64 / r,
            self.ld_h[3] as f64 / r,
            self.ld_e[3] as f64 / r,
            self.ld_h[2] as f64 / r,
            self.ld_e[2] as f64 / r,
        )
    }
}

/// Class of a load against the stores still in the store buffer (the last row's and this row's
/// so far, youngest last): see `Shape::ld_h`.
fn classify(sb: &[Acc], ld: Acc) -> usize {
    let over = |s: &Acc| s.buf == ld.buf && s.at < ld.at + ld.len && ld.at < s.at + s.len;
    match sb.iter().rev().find(|s| over(s)) {
        None => 0,
        Some(s) if s.at == ld.at && s.len == ld.len => 1,
        Some(s) if s.at <= ld.at && ld.at + ld.len <= s.at + s.len => 2,
        Some(_) => 3,
    }
}

/// The score-only loop of `Blast_SemiGappedAlign` (as `tests::reference`, the band only), with the
/// block layout of the 256-bit kernels laid over each row, and the loads and stores of the
/// kernel's rows. `rec` = the kernel's own (row offset, first, lanes) per row from a traceback run;
/// `w8` = 8-bit lanes (align8) or 16-bit (align16).
fn shape(p: &Prob<'_>, rec: &[(u32, u32, u32)], w8: bool, limit_rebase: i32, st: &mut Shape) {
    let (len1, len2) = (p.len1, p.len2);
    let goe = p.go + p.ge;
    let x = p.x.max(goe);
    let lb: i64 = if w8 { 1 } else { 2 };
    let m = |a: u8, b: u8| p.mat[a as usize * p.n + b as usize];
    let col = |b: usize| p.cols.get(b).copied().unwrap_or(0);
    let mut best = vec![MININT; len2 + 4];
    let mut best_gap = vec![MININT; len2 + 4];
    let mut score = -goe;
    best[0] = 0;
    best_gap[0] = -goe;
    let mut b_size = 1usize;
    for i in 1..=len2 {
        if score < -x {
            break;
        }
        best[i] = score;
        best_gap[i] = score - goe;
        score -= p.ge;
        b_size = i + 1;
    }
    // row 0: byte (i16) stores of H (buffer 0) and E
    let mut sbuf: Vec<Acc> = Vec::new();
    for i in 0..b_size as i64 {
        sbuf.push(Acc { row: 0, buf: 0, at: i * lb, len: lb });
        sbuf.push(Acc { row: 0, buf: 2, at: i * lb, len: lb });
    }
    let mut wend = [b_size, 0usize];
    let mut cur = 0usize;
    let (mut best_score, mut first) = (0i32, 0usize);
    let mut base = 0i32;
    let mut sh = Shape::default();
    for a in 1..=len1 {
        let Some(&(_, kfirst, klanes)) = rec.get(a) else {
            st.skipped += 1;
            return;
        };
        if kfirst as usize != first {
            st.skipped += 1;
            return;
        }
        let proc_end = first + klanes as usize;
        let bend = b_size;
        let n32 = if bend > first + 16 {
            (bend - first - 16).div_ceil(32)
        } else {
            0
        };
        let s16 = first + 32 * n32;
        if proc_end < s16 || (proc_end - s16) % 16 != 0 {
            st.skipped += 1;
            return;
        }
        let n16 = (proc_end - s16) / 16;
        sh.rows += 1;
        sh.cells += (b_size - first) as u64;
        sh.lanes += klanes as u64;
        sh.n32 += n32 as u64;
        sh.n16 += n16 as u64;
        sh.n16_ext += (0..n16).filter(|j| s16 + 16 * j >= bend).count() as u64;
        // stores older than the last row have left the store buffer
        sbuf.retain(|s| s.row + 1 >= a);
        let (rb, wb) = (cur as u8, (1 - cur) as u8);
        // rebase at the start of the row (scalar: the same test on the best so far)
        let rebase = best_score - base > limit_rebase;
        if rebase {
            sh.rebases += 1;
            base = best_score;
            let step = if w8 { 32 } else { 16 };
            let mut b = first;
            while b < b_size {
                for buf in [rb, 2] {
                    let at = b as i64 * lb;
                    sbuf.push(Acc { row: a, buf, at, len: 32 });
                }
                b += step;
            }
        }
        // blocks: H load at b0 - 1 (buffer rb), E load at b0, then the stores at b0
        let mut blocks: Vec<(usize, usize)> = Vec::new();
        for k in 0..n32 {
            if w8 {
                blocks.push((first + 32 * k, 32));
            } else {
                blocks.push((first + 32 * k, 16));
                blocks.push((first + 32 * k + 16, 16));
            }
        }
        for j in 0..n16 {
            blocks.push((s16 + 16 * j, 16));
        }
        for &(b0, w) in &blocks {
            let len = w as i64 * lb;
            let ldh = Acc { row: a, buf: rb, at: (b0 as i64 - 1) * lb, len };
            let lde = Acc { row: a, buf: 2, at: b0 as i64 * lb, len };
            sh.ld_h[classify(&sbuf, ldh)] += 1;
            sh.ld_e[classify(&sbuf, lde)] += 1;
            sbuf.push(Acc { row: a, buf: wb, at: b0 as i64 * lb, len });
            sbuf.push(Acc { row: a, buf: 2, at: b0 as i64 * lb, len });
        }
        // the scalar row (band only), with the best-raising blocks
        let f0 = first;
        let blk = move |b: usize| -> usize {
            if b < s16 {
                if w8 {
                    (b - f0) / 32
                } else {
                    (b - f0) / 16
                }
            } else {
                (if w8 { n32 } else { 2 * n32 }) + (b - s16) / 16
            }
        };
        let orig = first;
        let qc = p.rows[a - 1];
        let mut score_val = MININT;
        let mut gap_row = MININT;
        let mut last = first;
        let before = best_score;
        let mut raised_blk = usize::MAX;
        for b in orig..b_size {
            let sc = col(b);
            if p.cf && sc == FENCE_SENTRY {
                st.skipped += 1;
                return;
            }
            let gap_col = best_gap[b];
            let next = best[b].wrapping_add(m(qc, sc));
            if score_val < gap_col {
                score_val = gap_col;
            }
            if score_val < gap_row {
                score_val = gap_row;
            }
            if best_score.wrapping_sub(score_val) > x {
                if b == first {
                    first += 1;
                } else {
                    best[b] = MININT;
                }
            } else {
                last = b;
                if score_val > best_score {
                    best_score = score_val;
                    if blk(b) != raised_blk {
                        raised_blk = blk(b);
                        sh.best_blocks += 1;
                    }
                }
                gap_row -= p.ge;
                let ext = gap_col - p.ge;
                let open = score_val - goe;
                best_gap[b] = ext.max(open);
                gap_row = gap_row.max(open);
                best[b] = score_val;
            }
            score_val = next;
        }
        if best_score > before {
            sh.best_rows += 1;
        }
        // zeroing of the written H row past this row's end, up to what it held before
        {
            let step = if w8 { 32 } else { 16 };
            let wi = 1 - cur;
            let mut b = proc_end;
            while b < wend[wi] {
                sbuf.push(Acc { row: a, buf: wb, at: b as i64 * lb, len: 32 });
                sh.zero_st += 1;
                b += step;
            }
            wend[wi] = proc_end;
            cur = wi;
        }
        if first == b_size {
            break;
        }
        // E right of the last survivor
        {
            let step = if w8 { 32 } else { 16 };
            let mut b = last + 1;
            while b < proc_end {
                sbuf.push(Acc { row: a, buf: 2, at: b as i64 * lb, len: 32 });
                sh.eclear_st += 1;
                b += step;
            }
        }
        if last < b_size - 1 {
            b_size = last + 1;
        } else {
            while gap_row >= best_score - x && b_size <= len2 {
                best[b_size] = gap_row;
                best_gap[b_size] = gap_row - goe;
                gap_row -= p.ge;
                b_size += 1;
            }
        }
        if b_size <= len2 {
            best[b_size] = MININT;
            best_gap[b_size] = MININT;
            b_size += 1;
        }
        let d = first - orig;
        sh.adv[match d {
            0 => 0,
            1 => 1,
            2 => 2,
            3..=8 => 3,
            _ => 4,
        }] += 1;
    }
    sh.calls = 1;
    st.calls += sh.calls;
    st.rows += sh.rows;
    st.cells += sh.cells;
    st.lanes += sh.lanes;
    st.n32 += sh.n32;
    st.n16 += sh.n16;
    st.n16_ext += sh.n16_ext;
    st.best_rows += sh.best_rows;
    st.best_blocks += sh.best_blocks;
    st.rebases += sh.rebases;
    for i in 0..5 {
        st.adv[i] += sh.adv[i];
    }
    st.zero_st += sh.zero_st;
    st.eclear_st += sh.eclear_st;
    for i in 0..4 {
        st.ld_h[i] += sh.ld_h[i];
        st.ld_e[i] += sh.ld_e[i];
    }
}

/// Address modulo 4096 and length in bytes of each row buffer of the scratch.
fn layout(sc: &XdropScratch) -> String {
    let f = |n: &str, p: usize, b: usize| format!("{n}:{}/{}", p % 4096, b);
    [
        f("ha", sc.ha.as_ptr() as usize, sc.ha.len() * 2),
        f("hb", sc.hb.as_ptr() as usize, sc.hb.len() * 2),
        f("e", sc.e.as_ptr() as usize, sc.e.len() * 2),
        f("h8a", sc.h8a.as_ptr() as usize, sc.h8a.len()),
        f("h8b", sc.h8b.as_ptr() as usize, sc.h8b.len()),
        f("e8", sc.e8.as_ptr() as usize, sc.e8.len()),
        f("sb", sc.sb.as_ptr() as usize, sc.sb.len()),
        f("tabs", sc.tabs.as_ptr() as usize, sc.tabs.len()),
        f("trace", sc.trace.as_ptr() as usize, sc.trace.capacity()),
    ]
    .join(",")
}

/// `R1_LAYOUT=spread`: the rows `ha`, `hb`, `e`, `h8a`, `h8b`, `e8`, `sb` copied into one leaked
/// arena, each starting at its own offset modulo 4096 (64 + 576 i), with whole pages between them;
/// lengths kept. The timed loop replays problems the scratch has already been grown for, so no
/// `grow` (which would reallocate these vectors) runs while they live in the arena; the caller
/// forgets them before the scratch drops (`unspread`).
fn spread(sc: &mut XdropScratch) {
    let lens = [
        sc.ha.len() * 2,
        sc.hb.len() * 2,
        sc.e.len() * 2,
        sc.h8a.len(),
        sc.h8b.len(),
        sc.e8.len(),
        sc.sb.len(),
    ];
    let slot = |l: usize| l.div_ceil(4096) * 4096 + 4096;
    let total: usize = lens.iter().map(|&l| slot(l)).sum::<usize>() + 8192;
    let arena: &'static mut [u8] = Vec::leak(vec![0u8; total]);
    let base = arena.as_mut_ptr() as usize;
    let mut at = base.div_ceil(4096) * 4096;
    let mut ptrs = [0usize; 7];
    for (i, &l) in lens.iter().enumerate() {
        ptrs[i] = at + 64 + 576 * i;
        at += slot(l);
    }
    assert!(at <= base + total);
    // SAFETY: each target range lies inside the leaked arena and the ranges do not overlap; the
    // vectors built here are never freed or grown (see `unspread`).
    unsafe {
        let mv16 = |v: &mut Vec<i16>, p: usize| {
            let n = v.len();
            std::ptr::copy_nonoverlapping(v.as_ptr(), p as *mut i16, n);
            let old = std::mem::replace(v, Vec::from_raw_parts(p as *mut i16, n, n));
            drop(old);
        };
        let mv8 = |v: &mut Vec<u8>, p: usize| {
            let n = v.len();
            std::ptr::copy_nonoverlapping(v.as_ptr(), p as *mut u8, n);
            let old = std::mem::replace(v, Vec::from_raw_parts(p as *mut u8, n, n));
            drop(old);
        };
        mv16(&mut sc.ha, ptrs[0]);
        mv16(&mut sc.hb, ptrs[1]);
        mv16(&mut sc.e, ptrs[2]);
        mv8(&mut sc.h8a, ptrs[3]);
        mv8(&mut sc.h8b, ptrs[4]);
        mv8(&mut sc.e8, ptrs[5]);
        mv8(&mut sc.sb, ptrs[6]);
    }
}

/// Forgets the arena-backed rows of `spread` (they must not reach the allocator).
fn unspread(sc: &mut XdropScratch) {
    std::mem::forget(std::mem::take(&mut sc.ha));
    std::mem::forget(std::mem::take(&mut sc.hb));
    std::mem::forget(std::mem::take(&mut sc.e));
    std::mem::forget(std::mem::take(&mut sc.h8a));
    std::mem::forget(std::mem::take(&mut sc.h8b));
    std::mem::forget(std::mem::take(&mut sc.e8));
    std::mem::forget(std::mem::take(&mut sc.sb));
}

/// One call as `xdrop_align` makes it with the 256-bit kernels on: `x_avx2::align`; with
/// `R1_KERN=copy` / `nodirty` (8-bit lanes only) one of the probe copies of `align8` instead.
#[inline(never)]
fn call(p: &Prob<'_>, tb: bool, sc: &mut XdropScratch) -> Option<XdropResult> {
    let kern = KERN.with(|k| k.get());
    // SAFETY: `available()` is checked before any call.
    unsafe {
        if kern != 0 && fits_8bit(&p.scores, p.go, p.ge, p.x).is_none() {
            // 16-bit problems: the align16 copy (copy, emask, emaskld)
            return match (tb, kern) {
                (false, 1) => probe::align16_probe::<false, 0>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 13) => probe::align16_probe::<false, 9>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 14) => probe::align16_probe::<false, 10>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (true, 1) => probe::align16_probe::<true, 0>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (true, 13) => probe::align16_probe::<true, 9>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (true, 14) => probe::align16_probe::<true, 10>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                _ => panic!("kernel {kern} has no 16-bit copy"),
            };
        }
        if kern != 0 {
            let m = fits_8bit(&p.scores, p.go, p.ge, p.x)?;
            return match (tb, kern) {
                (false, 13) => probe::align8_probe::<false, 9, 1000>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 14) => probe::align8_probe::<false, 10, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (true, 13) => probe::align8_probe::<true, 9, 1000>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (true, 14) => probe::align8_probe::<true, 10, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 1) => probe::align8_probe::<false, 0, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 3) => probe::align8_probe::<false, 2, 4>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 4) => probe::align8_probe::<false, 3, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 5) => probe::align8_probe::<false, 4, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 6) => probe::align8_probe::<false, 5, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 7) => probe::align8_probe::<false, 7, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 8) => probe::align8_probe::<false, 7, 4>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 11) => probe::align8_probe::<false, 0, 40>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 12) => probe::align8_probe::<false, 0, 400>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 10) => probe::align8_probe::<false, 8, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, 9) => probe::align8_probe::<false, 6, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (false, _) => probe::align8_probe::<false, 1, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (true, 1) => probe::align8_probe::<true, 0, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
                (true, _) => probe::align8_probe::<true, 1, 0>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, p.cf, sc,
                ),
            };
        }
        align(
            &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, tb, p.cf, sc,
        )
    }
}

thread_local! {
    /// 0 = `align`, 1 = probe copy, 2 = probe copy without the high-water-mark branch.
    static KERN: std::cell::Cell<u8> = const { std::cell::Cell::new(0) };
    /// The band table of the problem being run (`oracle`: recorded by kernel 9, replayed by 7 / 8).
    static TAB: std::cell::Cell<*mut Vec<(u32, u32, u32)>> =
        const { std::cell::Cell::new(std::ptr::null_mut()) };
}

#[test]
#[ignore = "needs R1_DIR=<capture directory> and R1_CAT"]
fn r1_loop() {
    let Some(dir) = var("R1_DIR") else {
        eprintln!("R1_DIR is not set");
        return;
    };
    assert!(available(), "the harness needs AVX2");
    let cat = var("R1_CAT").unwrap_or_else(|| "so8".into());
    let (tb, w8) = match cat.as_str() {
        "so8" => (false, true),
        "tb8" => (true, true),
        "so16" => (false, false),
        "tb16" => (true, false),
        _ => panic!("R1_CAT must be so8, tb8, so16 or tb16"),
    };
    let reps: usize = var("R1_REPS").and_then(|v| v.parse().ok()).unwrap_or(0);
    let limit: usize = var("R1_LIMIT")
        .and_then(|v| v.parse().ok())
        .unwrap_or(usize::MAX);
    let synth = var("R1_SYNTH");
    let kern = match var("R1_KERN").as_deref() {
        None | Some("orig") => 0u8,
        Some("copy") => 1,
        Some("nodirty") => 2,
        Some("lat4") => 3,
        Some("addr4") => 4,
        Some("nop30") => 5,
        Some("alu20") => 6,
        Some("oracle") => 7,
        Some("oraclelat4") => 8,
        Some("oraclez") => 10,
        Some("elat4") => 11,
        Some("hslat4") => 12,
        Some("emask") => 13,
        Some("emaskld") => 14,
        Some(k) => panic!(
            "R1_KERN must be orig, copy, nodirty, lat4, addr4, nop30, alu20, oracle, oraclelat4, oraclez, elat4, hslat4, emask or emaskld, not {k}"
        ),
    };
    assert!(
        kern == 0 || w8 || matches!(kern, 1 | 13 | 14),
        "the 16-bit copy runs copy, emask and emaskld only"
    );
    assert!(
        kern <= 2 || !tb || matches!(kern, 13 | 14),
        "the latency / issue probes are score-only"
    );
    KERN.with(|k| k.set(kern));
    let t_load = Instant::now();
    let probs = load_dir(std::path::Path::new(&dir)).expect("read the capture files");
    // synthetic problems: B3's construction (perf-of-20261013 wip/b3_window, b3r)
    let mut syn: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut list: Vec<Prob<'_>> = Vec::new();
    let mut p0: Option<usize> = None;
    for (i, p) in probs.iter().enumerate() {
        if !p.tb && fits_8bit(&p.scores(), p.gap_open, p.gap_extend, p.x_drop).is_some() {
            p0 = Some(i);
            break;
        }
    }
    if let Some(label) = &synth {
        let (sub_pct, indel_pct) = match label.as_str() {
            "identical" => (0u64, 0u64),
            "80pct" => (20, 1),
            "60pct" => (40, 2),
            "45pct" => (55, 3),
            _ => panic!("R1_SYNTH must be identical, 80pct, 60pct or 45pct"),
        };
        let aa: [u8; 20] = [
            1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
        ];
        let mut rng = Rng(0xB3B3_0001);
        let len = 3000usize;
        let a: Vec<u8> = (0..len).map(|_| aa[rng.below(20) as usize]).collect();
        let mut b = Vec::with_capacity(len + 100);
        for &c in &a {
            let r = rng.below(100);
            if r < indel_pct {
                if rng.below(2) == 0 {
                    continue;
                }
                b.push(aa[rng.below(20) as usize]);
                b.push(c);
            } else if r < indel_pct + sub_pct {
                b.push(aa[rng.below(20) as usize]);
            } else {
                b.push(c);
            }
        }
        syn.push((a, b));
    }
    if synth.is_some() {
        let pp = &probs[p0.expect("an 8-bit score-only problem for the synthetic matrix")];
        let (a, b) = &syn[0];
        list.push(Prob {
            q: RowSeq::Bytes {
                data: a,
                base: -1,
                step: 1,
            },
            s: ColSeq {
                data: b,
                base: 0,
                step: 1,
                zero_from: usize::MAX,
            },
            scores: pp.scores(),
            len1: a.len(),
            len2: b.len(),
            go: pp.gap_open,
            ge: pp.gap_extend,
            x: pp.x_drop,
            tb,
            cf: false,
            rows: a,
            cols: b,
            mat: &pp.matrix,
            n: pp.n,
            want: None,
        });
    } else {
        for p in &probs {
            if list.len() >= limit {
                break;
            }
            let m8 = fits_8bit(&p.scores(), p.gap_open, p.gap_extend, p.x_drop);
            if p.tb != tb || m8.is_some() != w8 || p.result.is_none() {
                continue;
            }
            list.push(Prob {
                q: p.row_seq(),
                s: p.col_seq(),
                scores: p.scores(),
                len1: p.len1,
                len2: p.len2,
                go: p.gap_open,
                ge: p.gap_extend,
                x: p.x_drop,
                tb: p.tb,
                cf: p.check_fence,
                rows: &p.rows,
                cols: &p.cols,
                mat: &p.matrix,
                n: p.n,
                want: p.result,
            });
        }
    }
    let load_s = t_load.elapsed().as_secs_f64();
    // check and shape (once)
    let (mut sc, mut sd) = (XdropScratch::new(), XdropScratch::new());
    let mut st = Shape::default();
    let (mut declined, mut ok) = (0usize, 0usize);
    let synth_want = if synth.is_some() {
        // the synthetic problem against the 128-bit kernel of the dispatch
        let p = &list[0];
        let m8 = fits_8bit(&p.scores, p.go, p.ge, p.x);
        // SAFETY: AVX2 checked above.
        Some(unsafe {
            match (tb, m8) {
                (true, Some(m)) => super::super::kernel::align8_avx::<true>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, false, &mut sd,
                ),
                (false, Some(m)) => super::super::kernel::align8_avx::<false>(
                    &p.q, &p.s, &p.scores, m, p.len1, p.len2, p.go, p.ge, p.x, false, &mut sd,
                ),
                (true, None) => super::super::kernel::align_avx::<true>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, false, &mut sd,
                ),
                (false, None) => super::super::kernel::align_avx::<false>(
                    &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, false, &mut sd,
                ),
            }
        })
    } else {
        None
    };
    // oracle: each problem's band, recorded once by kernel 9 before the checked call replays it
    let mut tabs: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); list.len()];
    for (i, p) in list.iter().enumerate() {
        if matches!(kern, 7 | 8 | 10) {
            TAB.with(|t| t.set(&mut tabs[i]));
            KERN.with(|k| k.set(9));
            // a declined problem stops at the same row start in both passes (no table read past it)
            let _ = call(p, tb, &mut sc);
            KERN.with(|k| k.set(kern));
        }
        let r = call(p, tb, &mut sc);
        let want = if synth.is_some() { synth_want.flatten() } else { p.want };
        match r {
            Some(r) => {
                assert_eq!(Some(r), want, "256-bit result vs recorded / 128-bit");
                ok += 1;
            }
            None => {
                declined += 1;
                continue;
            }
        }
        // the kernel's rows from a traceback run (the band evolves the same way)
        KERN.with(|k| k.set(0));
        let shaped = call(p, true, &mut sd).is_some();
        KERN.with(|k| k.set(kern));
        if shaped {
            let goe = p.go + p.ge;
            let x = p.x.max(goe);
            let limit_rebase = if w8 {
                let mm = fits_8bit(&p.scores, p.go, p.ge, p.x).unwrap();
                127 - (x + goe + mm + 1) - mm
            } else {
                let bias = x + goe + super::super::MAX_M + 1;
                (super::super::VMAX as i32) - bias - super::super::MAX_M
            };
            let rec: Vec<(u32, u32, u32)> = sd.rows.clone();
            shape(p, &rec, w8, limit_rebase, &mut st);
        } else {
            st.skipped += 1;
        }
    }
    // the loop the counters are taken over
    let lay_before = layout(&sc);
    let spread_on = var("R1_LAYOUT").as_deref() == Some("spread");
    if spread_on {
        spread(&mut sc);
    }
    let lay = layout(&sc);
    let t0 = Instant::now();
    for _ in 0..reps {
        for (i, p) in list.iter().enumerate() {
            if matches!(kern, 7 | 8 | 10) {
                TAB.with(|t| t.set(&mut tabs[i]));
            }
            std::hint::black_box(call(p, tb, &mut sc));
        }
    }
    let loop_s = t0.elapsed().as_secs_f64();
    if spread_on {
        unspread(&mut sc);
    }
    let rows_total = st.rows as f64 * reps as f64;
    let line = format!(
        "R1 cat={cat} synth={} dir={dir} problems={} handled={ok} declined={declined} reps={reps} \
         load_s={load_s:.3} loop_s={loop_s:.6} ns_per_row={:.3} layout={} lay0={lay_before} \
         lay={lay} kern={} | {}\n",
        synth.as_deref().unwrap_or("-"),
        list.len(),
        if rows_total > 0.0 {
            loop_s * 1e9 / rows_total
        } else {
            0.0
        },
        if spread_on { "spread" } else { "as-is" },
        ["orig", "copy", "nodirty", "lat4", "addr4", "nop30", "alu20", "oracle", "oraclelat4", "-", "oraclez", "elat4",
          "hslat4", "emask", "emaskld"]
            [kern as usize],
        st.line()
    );
    eprint!("{line}");
    if let Some(out) = var("R1_OUT") {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(out)
            .expect("open R1_OUT");
        f.write_all(line.as_bytes()).expect("write R1_OUT");
    }
}

// ---------------------------------------------------------------------------
// Probe copies of the kernel (R1_KERN): `copy` = align8 verbatim (V = 0), `nodirty` = the same with
// the per-row high-water mark `sc.dirty8` taken as a max instead of a branch (V = 1); score-only
// slope probes, same results: `lat4` (V = 2) four identity vpor on hh in every block (latency on the
// H / E / drop chain), `addr4` (V = 3) four identity adds on the band position before the next row,
// `nop30` (V = 4) 30 NOPs per row, `alu20` (V = 5) 20 independent adds per row; `oracle` (V = 7,
// recorded by V = 6) takes the next row's band from a recorded run instead of this row's survivors,
// `oraclelat4` adds lat4's latency to it, `oraclez` (V = 8) also clears E after the recorded last
// survivor instead of this row's; `elat4` / `hslat4` (L = 40 / 400) put lat4's four cycles on the
// stored E / the stored H only.
// ---------------------------------------------------------------------------
#[allow(clippy::all)]
mod probe {
    use super::super::super::kernel::finish;
    use super::super::super::{
        build_row_table, ColSeq, RowSeq, Scores, XdropResult, XdropScratch, FENCE_SENTRY, PAD,
        SCRIPT_EXTEND_GAP_A, SCRIPT_EXTEND_GAP_B, SCRIPT_GAP_IN_A, SCRIPT_SUB, TAB,
    };
    use super::super::{
        bcast128, build_row, fill_cols, last16, last8, ld, ldx, lut32, lut32x, pmax16, pmax8,
        pmax8x, st, stx, up1_16, up1_8, Row, RowX, C16, C8, KB7X2, LANE_MASK16, LANE_MASK8,
        LANE_MASK8X, W, X,
    };
    use super::super::super::{MAX_M, VMAX};
    use core::arch::x86_64::*;
    #[allow(unused_imports)]
    use super::super::super::fits_8bit;

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn block8_p<const TB: bool, const LAT: u32>(
        c: &C8,
        rs: &mut Row,
        hr: *const u8,
        hw: *mut u8,
        e: *mut u8,
        b: usize,
        s: W,
        lane_mask: Option<W>,
        hpm: W,
        emk: Option<W>,
    ) -> (W, W) {
        let hp = _mm256_andnot_si256(hpm, ld(hr.add(b).sub(1)));
        let mut ep = ld(e.add(b));
        if let Some(m) = emk {
            // probe (emaskld): lanes right of the last row's last survivor read as minus infinity
            ep = _mm256_and_si256(ep, m);
        }
        let d = _mm256_adds_epi8(hp, s);
        let h0 = _mm256_max_epi8(d, ep);
        let v = _mm256_subs_epu8(h0, c.goe);
        let fi = _mm256_subs_epu8(pmax8(_mm256_adds_epu8(up1_8(v), c.ramp), c.kb15), c.ramp);
        let f = _mm256_max_epu8(fi, _mm256_subs_epu8(rs.carry, c.ramp));
        let mut hh = _mm256_max_epi8(h0, f);
        if let Some(m) = lane_mask {
            hh = _mm256_and_si256(hh, m);
        }
        for _ in 0..(if LAT < 10 { LAT } else { 0 }) {
            // probe: one cycle of latency on the H / E / drop chain (vpor x, x, x = x)
            core::arch::asm!("vpor {0}, {0}, {0}", inout(ymm_reg) hh, options(pure, nomem, nostack));
        }
        let hg = _mm256_subs_epu8(hh, c.goe);
        let fs = _mm256_subs_epu8(f, c.ge1);
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
        let mut hs = _mm256_andnot_si256(drop, hh);
        let mut en = en;
        for _ in 0..(if LAT == 400 { 4 } else { 0 }) {
            // probe: four cycles on the stored H only (the next row's H loads)
            core::arch::asm!("vpor {0}, {0}, {0}", inout(ymm_reg) hs, options(pure, nomem, nostack));
        }
        for _ in 0..(if LAT == 40 { 4 } else { 0 }) {
            // probe: four cycles on the stored E only (the next row's E loads)
            core::arch::asm!("vpor {0}, {0}, {0}", inout(ymm_reg) en, options(pure, nomem, nostack));
        }
        if LAT == 1000 {
            // probe (emask): E of a dropped cell stored as minus infinity
            en = _mm256_andnot_si256(drop, en);
        }
        st(hw.add(b), hs);
        st(e.add(b), en);
        let scr = if TB {
            let me = _mm256_cmpgt_epi8(ep, d);
            let mf = _mm256_cmpgt_epi8(f, h0);
            let op = _mm256_andnot_si256(mf, _mm256_add_epi8(c.c3, _mm256_and_si256(me, c.c3)));
            let fb = _mm256_andnot_si256(_mm256_cmpgt_epi8(hg, es), c.c40);
            let fa = _mm256_andnot_si256(_mm256_cmpgt_epi8(hg, fs), c.c10);
            _mm256_or_si256(op, _mm256_or_si256(fb, fa))
        } else {
            _mm256_setzero_si256()
        };
        (drop, scr)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn block8x_p<const TB: bool, const LAT: u32>(
        c: &C8,
        rs: &mut RowX,
        rec_b: &mut usize,
        slow: &mut bool,
        hr: *const u8,
        hw: *mut u8,
        e: *mut u8,
        b: usize,
        s: X,
        lane_mask: Option<X>,
        hpm: X,
        emk: Option<X>,
    ) -> (X, X) {
        let goe = _mm256_castsi256_si128(c.goe);
        let ge1 = _mm256_castsi256_si128(c.ge1);
        let ramp = _mm256_castsi256_si128(c.ramp);
        let xv = _mm256_castsi256_si128(c.xv);
        let k15 = _mm256_castsi256_si128(c.kb15);
        let hp = _mm_andnot_si128(hpm, ldx(hr.add(b).sub(1)));
        let mut ep = ldx(e.add(b));
        if let Some(m) = emk {
            ep = _mm_and_si128(ep, m);
        }
        let d = _mm_adds_epi8(hp, s);
        let h0 = _mm_max_epi8(d, ep);
        let v = _mm_subs_epu8(h0, goe);
        let fi = _mm_subs_epu8(pmax8x(_mm_adds_epu8(_mm_bslli_si128::<1>(v), ramp)), ramp);
        let f = _mm_max_epu8(fi, _mm_subs_epu8(rs.carry, ramp));
        let mut hh = _mm_max_epi8(h0, f);
        if let Some(m) = lane_mask {
            hh = _mm_and_si128(hh, m);
        }
        for _ in 0..(if LAT < 10 { LAT } else { 0 }) {
            core::arch::asm!("vpor {0}, {0}, {0}", inout(xmm_reg) hh, options(pure, nomem, nostack));
        }
        let hg = _mm_subs_epu8(hh, goe);
        let fs = _mm_subs_epu8(f, ge1);
        let out = _mm_shuffle_epi8(_mm_max_epu8(_mm_subs_epu8(fi, ge1), v), k15);
        rs.carry = _mm_max_epu8(_mm_subs_epu8(rs.carry, _mm256_castsi256_si128(c.ge16)), out);
        let es = _mm_subs_epu8(ep, ge1);
        let en = _mm_max_epi8(hg, es);
        let gt = _mm_cmpgt_epi8(hh, rs.bestv);
        let thr;
        if _mm_movemask_epi8(gt) != 0 {
            let pfx = _mm_max_epu8(pmax8x(_mm_bslli_si128::<1>(hh)), rs.bestv);
            thr = _mm_subs_epu8(pfx, xv);
            let rm = _mm_movemask_epi8(_mm_cmpgt_epi8(hh, pfx)) as u32;
            *rec_b = b + (31 - rm.leading_zeros()) as usize;
            rs.bestv = _mm_shuffle_epi8(_mm_max_epu8(pfx, hh), k15);
            rs.thrv = _mm_subs_epu8(rs.bestv, xv);
            *slow = true;
        } else {
            thr = rs.thrv;
        }
        let drop = _mm_cmpgt_epi8(thr, hh);
        let mut hs = _mm_andnot_si128(drop, hh);
        let mut en = en;
        for _ in 0..(if LAT == 400 { 4 } else { 0 }) {
            core::arch::asm!("vpor {0}, {0}, {0}", inout(xmm_reg) hs, options(pure, nomem, nostack));
        }
        for _ in 0..(if LAT == 40 { 4 } else { 0 }) {
            core::arch::asm!("vpor {0}, {0}, {0}", inout(xmm_reg) en, options(pure, nomem, nostack));
        }
        if LAT == 1000 {
            en = _mm_andnot_si128(drop, en);
        }
        stx(hw.add(b), hs);
        stx(e.add(b), en);
        let scr = if TB {
            let c3 = _mm256_castsi256_si128(c.c3);
            let me = _mm_cmpgt_epi8(ep, d);
            let mf = _mm_cmpgt_epi8(f, h0);
            let op = _mm_andnot_si128(mf, _mm_add_epi8(c3, _mm_and_si128(me, c3)));
            let fb = _mm_andnot_si128(_mm_cmpgt_epi8(hg, es), _mm256_castsi256_si128(c.c40));
            let fa = _mm_andnot_si128(_mm_cmpgt_epi8(hg, fs), _mm256_castsi256_si128(c.c10));
            _mm_or_si128(op, _mm_or_si128(fb, fa))
        } else {
            _mm_setzero_si128()
        };
        (drop, scr)
    }

    #[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
    pub(super) unsafe fn align8_probe<const TB: bool, const V: u32, const L: u32>(
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
        // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 31 * gap_extend to
        // fit 8 bits (v <= 127).
        if gap_extend > 4 {
            return None;
        }
        let bias = bias as i32;
        let limit_rebase = limit_rebase as i32;
        let num_extra = (x / gap_extend + 3) as usize;
        let ncols = scores.ncols();
        sc.reset();
    
        // Row 0 goes to buffer A (as in `align_body8`).
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
        if TB {
            sc.trace.resize(b_size, SCRIPT_GAP_IN_A);
            sc.rows.push((0, 0, b_size as u32));
        }
    
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
            c3: _mm256_set1_epi8(SCRIPT_SUB as i8),
            c10: _mm256_set1_epi8(SCRIPT_EXTEND_GAP_A as i8),
            c40: _mm256_set1_epi8(SCRIPT_EXTEND_GAP_B as i8),
        };
        let zero = _mm256_setzero_si256();
        let lane0 = _mm256_setr_epi8(
            -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0,
        );
        let k70 = _mm256_set1_epi8(0x70);
        let k10 = _mm256_set1_epi8(0x10);
        let static_tabs: Option<*const u8> = match scores {
            Scores::Tables { t, .. } => Some(t.data.as_ptr()),
            _ => None,
        };
        let small_alphabet = ncols <= 16;
        let mut cells: u64 = 0;
        // V = 6 records each row's next band (first, b_size); V = 7 replays it (see the row end).
        let tabv: *mut Vec<(u32, u32, u32)> = if V >= 6 {
            super::TAB.with(|t| t.get())
        } else {
            core::ptr::null_mut()
        };
        if V == 6 {
            (*tabv).clear();
        }
        let tab: *const (u32, u32, u32) = if V >= 7 { (*tabv).as_ptr() } else { core::ptr::null() };

        let mut base: i32 = 0;
        let mut best_rel: i32 = 0;
        let mut a_off = 0usize;
        let mut b_off = 0usize;
        let mut first = 0usize;
        let mut fence_hit = false;
        let mut wend = [b_size, 0usize];
        let mut cur = 0usize;
        // V = 10: one past the last row's last survivor (row 0: its written lanes)
        let mut lastp1 = b_size;
    
        for a in 1..=len1 {
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
    
            if best_rel > limit_rebase {
                let dv = _mm256_set1_epi8(best_rel as u8 as i8);
                let mut b = first;
                while b < b_size {
                    st(hr.add(b), _mm256_subs_epu8(ld(hr.add(b)), dv));
                    st(e.add(b), _mm256_subs_epu8(ld(e.add(b)), dv));
                    b += 32;
                }
                base += best_rel;
                best_rel = 0;
            }
            let mut hpm = lane0;
            if V == 4 {
                // probe: 30 NOPs per row (issue slots and ROB entries, no scheduler entry or port)
                core::arch::asm!(".rept 30", "nop", ".endr", options(nomem, nostack, preserves_flags));
            }
            if V == 5 {
                // probe: 20 independent ALU uops per row (four chains of five adds)
                let (mut d0, mut d1, mut d2, mut d3) = (0u64, 0u64, 0u64, 0u64);
                core::arch::asm!(".rept 5", "add {0}, 1", "add {1}, 1", "add {2}, 1", "add {3}, 1", ".endr",
                    inout(reg) d0, inout(reg) d1, inout(reg) d2, inout(reg) d3, options(pure, nomem, nostack));
                core::hint::black_box((d0, d1, d2, d3));
            }
    
            let row_off = sc.trace.len();
            if TB {
                sc.trace.reserve(bend - first + num_extra + 96);
            }
            let tr = sc.trace.as_mut_ptr().add(row_off);
    
            let bestv0 = _mm256_set1_epi8((bias + best_rel) as u8 as i8);
            let mut rs = Row {
                carry: zero,
                bestv: bestv0,
                thrv: _mm256_subs_epu8(bestv0, c.xv),
                rec_b: 0,
                slow: false,
            };
            let mut lo_surv = usize::MAX;
            let mut hi_surv = 0usize;
            let mut b0 = first;
            // 32-lane blocks while more than 16 lanes of the band are left, then 16-lane blocks for the
            // rest of the band and the row-gap extension (stop test as in `align16`).
            while b0 + 16 < bend {
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
                let emk = if V == 10 {
                    let n = lastp1.saturating_sub(b0).min(32);
                    Some(ld(LANE_MASK8.as_ptr().add(32 - n)))
                } else {
                    None
                };
                let (d0, r0) = block8_p::<TB, L>(&c, &mut rs, hr, hw, e, b0, sv, mask, hpm, emk);
                hpm = zero;
                if TB {
                    st(tr.add(b0 - first), r0);
                }
                let live = !(_mm256_movemask_epi8(d0) as u32);
                if live != 0 {
                    if lo_surv == usize::MAX {
                        lo_surv = b0 + live.trailing_zeros() as usize;
                    }
                    hi_surv = b0 + (31 - live.leading_zeros() as usize);
                }
                b0 += 32;
            }
            let mut rx = RowX {
                carry: _mm256_castsi256_si128(rs.carry),
                bestv: _mm256_castsi256_si128(rs.bestv),
                thrv: _mm256_castsi256_si128(rs.thrv),
            };
            let (mut rec_b, mut slow) = (rs.rec_b, rs.slow);
            let t0x = _mm256_castsi256_si128(t0);
            let t1x = _mm256_castsi256_si128(t1);
            loop {
                if b0 >= bend {
                    if b0 >= limit {
                        break;
                    }
                    if (_mm_cvtsi128_si32(rx.carry) & 0xFF) < (_mm_cvtsi128_si32(rx.thrv) & 0xFF) {
                        break;
                    }
                }
                let idx = ldx(sb.add(b0));
                let sv = if small_alphabet {
                    _mm_shuffle_epi8(t0x, idx)
                } else {
                    lut32x(t0x, t1x, idx)
                };
                let mask = if b0 + 16 <= limit {
                    None
                } else {
                    Some(ldx(LANE_MASK8X.as_ptr().add(16 - (limit - b0))))
                };
                let emk = if V == 10 {
                    let n = lastp1.saturating_sub(b0).min(16);
                    Some(ldx(LANE_MASK8X.as_ptr().add(16 - n)))
                } else {
                    None
                };
                let (d0, r0) = block8x_p::<TB, L>(
                    &c,
                    &mut rx,
                    &mut rec_b,
                    &mut slow,
                    hr,
                    hw,
                    e,
                    b0,
                    sv,
                    mask,
                    _mm256_castsi256_si128(hpm),
                    emk,
                );
                hpm = zero;
                if TB {
                    stx(tr.add(b0 - first), r0);
                }
                let live = (!(_mm_movemask_epi8(d0) as u32)) & 0xFFFF;
                if live != 0 {
                    if lo_surv == usize::MAX {
                        lo_surv = b0 + live.trailing_zeros() as usize;
                    }
                    hi_surv = b0 + (31 - live.leading_zeros() as usize);
                }
                b0 += 16;
            }
            let proc_end = b0;
            if V == 1 {
                // probe: the high-water mark without a branch (same value)
                sc.dirty8 = sc.dirty8.max(proc_end);
            } else if proc_end > sc.dirty8 {
                sc.dirty8 = proc_end;
            }
            let wi = 1 - cur;
            let mut b = proc_end;
            while b < wend[wi] {
                st(hw.add(b), zero);
                b += 32;
            }
            wend[wi] = proc_end;
            cur = wi;
            if slow {
                b_off = rec_b;
                a_off = a;
                best_rel = (_mm_cvtsi128_si32(rx.bestv) & 0xFF) - bias;
            }
            if TB {
                sc.trace.set_len(row_off + proc_end - first);
                sc.rows
                    .push((row_off as u32, first as u32, (proc_end - first) as u32));
            }
            if fence_row {
                fence_hit = true;
                break;
            }
            if lo_surv == usize::MAX {
                break;
            }
            let mut last = hi_surv;
            if V == 8 {
                // probe: the E clearing after the last survivor at the recorded position, so that its
                // store addresses no longer wait for this row's survivors
                core::arch::asm!("/* {0} */", in(reg) last, options(nomem, nostack, preserves_flags));
                last = (*tab.add(a - 1)).2 as usize;
            }
            if V != 9 && V != 10 {
                let mut b = last + 1;
                while b < proc_end {
                    st(e.add(b), zero);
                    b += 32;
                }
            }
            lastp1 = last + 1;
            let (mut lo, mut hi) = (lo_surv, last);
            if V == 3 {
                // probe: four cycles of latency on the band position (next row's addresses, bounds)
                for _ in 0..4 {
                    core::arch::asm!("add {0}, 0", "add {1}, 0", inout(reg) lo, inout(reg) hi,
                        options(pure, nomem, nostack));
                }
            }
            first = lo;
            b_size = hi + 1;
            if b_size <= len2 {
                b_size += 1;
            }
            if V == 6 {
                (*tabv).push((first as u32, b_size as u32, last as u32));
            }
            if V == 7 || V == 8 {
                // probe: the next row's band from the recorded run of the same problem, so that it
                // no longer waits for this row's survivors (which are still computed, then dropped)
                core::arch::asm!("/* {0} {1} */", in(reg) first, in(reg) b_size,
                    options(nomem, nostack, preserves_flags));
                let (f, bs, _) = *tab.add(a - 1);
                first = f as usize;
                b_size = bs as usize;
            }
        }
    
        let best = base + best_rel;
        finish::<TB>(sc, a_off, b_off, best, fence_hit, cells)
    }

    // align16 copy (R1_KERN copy / emask / emaskld on 16-bit problems)
    #[inline(always)]
    unsafe fn em16<const V: u32>(lastp1: usize, b: usize) -> Option<W> {
        if V == 10 {
            let n = lastp1.saturating_sub(b).min(16);
            Some(ld(LANE_MASK16.as_ptr().add(16 - n)))
        } else {
            None
        }
    }

    #[inline(always)]
    unsafe fn block16_p<const TB: bool, const EM: u32>(
        c: &C16,
        rs: &mut Row,
        hr: *const i16,
        hw: *mut i16,
        e: *mut i16,
        b: usize,
        s: W,
        lane_mask: Option<W>,
        hpm: W,
        emk: Option<W>,
    ) -> (W, W) {
        // `hpm`: lanes with no diagonal (the first lane of the band), as a 0 in the H row
        let hp = _mm256_andnot_si256(hpm, ld(hr.add(b).sub(1)));
        let mut ep = ld(e.add(b));
        if let Some(m) = emk {
            // probe (emaskld): lanes right of the last row's last survivor read as minus infinity
            ep = _mm256_and_si256(ep, m);
        }
        let d = _mm256_adds_epi16(hp, s);
        let h0 = _mm256_max_epi16(d, ep);
        let v = _mm256_subs_epu16(h0, c.goe);
        // row gap opened inside the block: sat(max_j (v[j] + (j+1)g) - i*g) over j < i
        let fi = _mm256_subs_epu16(pmax16(_mm256_adds_epu16(up1_16(v), c.ramp), c.kb7), c.ramp);
        let f = _mm256_max_epi16(fi, _mm256_subs_epu16(rs.carry, c.ramp));
        let mut hh = _mm256_max_epi16(h0, f);
        if let Some(m) = lane_mask {
            hh = _mm256_and_si256(hh, m);
        }
        let hg = _mm256_subs_epu16(hh, c.goe);
        let fs = _mm256_subs_epu16(f, c.ge1);
        let out = last16(_mm256_max_epi16(_mm256_subs_epu16(fi, c.ge1), v), c.kb7);
        rs.carry = _mm256_max_epi16(_mm256_subs_epu16(rs.carry, c.ge16), out);
        let es = _mm256_subs_epu16(ep, c.ge1);
        let mut en = _mm256_max_epi16(hg, es);
        let gt = _mm256_cmpgt_epi16(hh, rs.bestv);
        let thr;
        if _mm256_movemask_epi8(gt) != 0 {
            // A new best inside this block: exact running best per lane.
            let pfx = _mm256_max_epi16(pmax16(up1_16(hh), c.kb7), rs.bestv);
            thr = _mm256_subs_epu16(pfx, c.xv);
            let rm = _mm256_movemask_epi8(_mm256_cmpgt_epi16(hh, pfx)) as u32;
            rs.rec_b = b + ((31 - rm.leading_zeros()) >> 1) as usize;
            rs.bestv = last16(_mm256_max_epi16(pfx, hh), c.kb7);
            rs.thrv = _mm256_subs_epu16(rs.bestv, c.xv);
            rs.slow = true;
        } else {
            thr = rs.thrv;
        }
        let drop = _mm256_cmpgt_epi16(thr, hh);
        st(hw.add(b), _mm256_andnot_si256(drop, hh));
        if EM == 9 {
            // probe (emask): E of a dropped cell stored as minus infinity
            en = _mm256_andnot_si256(drop, en);
        }
        st(e.add(b), en);
        let scr = if TB {
            let me = _mm256_cmpgt_epi16(ep, d);
            let mf = _mm256_cmpgt_epi16(f, h0);
            let op = _mm256_andnot_si256(mf, _mm256_add_epi16(c.c3, _mm256_and_si256(me, c.c3)));
            let fb = _mm256_andnot_si256(_mm256_cmpgt_epi16(hg, es), c.c40);
            let fa = _mm256_andnot_si256(_mm256_cmpgt_epi16(hg, fs), c.c10);
            _mm256_or_si256(op, _mm256_or_si256(fb, fa))
        } else {
            _mm256_setzero_si256()
        };
        (drop, scr)
    }

    #[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
    pub(super) unsafe fn align16_probe<const TB: bool, const V: u32>(
        q: &RowSeq<'_>,
        s: &ColSeq<'_>,
        scores: &Scores<'_>,
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
        if gap_extend < 1 || gap_open < 0 {
            return None;
        }
        let bias = x as i64 + goe as i64 + MAX_M as i64 + 1;
        let limit_rebase = VMAX - bias - MAX_M as i64;
        if limit_rebase < MAX_M as i64 || gap_extend as i64 * 8 > VMAX {
            return None;
        }
        // No NCBI counterpart: the prefix-maximum form of the row gap needs v + 15 * gap_extend to
        // fit 16 bits (v <= 32767).
        if gap_extend > 2048 {
            return None;
        }
        let bias = bias as i32;
        let limit_rebase = limit_rebase as i32;
        let num_extra = (x / gap_extend + 3) as usize;
        sc.reset();

        // Row 0 goes to buffer A (as in `align_body`).
        sc.ensure(num_extra.min(len2 + 2) + 16);
        let mut b_size = 1usize;
        {
            let h = sc.ha.as_mut_ptr().add(PAD);
            let e = sc.e.as_mut_ptr().add(PAD);
            *h = bias as i16;
            *e = (bias - goe) as i16;
            let mut score = -goe;
            for i in 1..=len2 {
                if score < -x {
                    break;
                }
                *h.add(i) = (bias + score) as i16;
                *e.add(i) = (bias + score - goe) as i16;
                score -= gap_extend;
                b_size = i + 1;
            }
            // In row 1 the scalar loop never forms the diagonal into column b_size.
            *h.add(b_size - 1) = 0;
        }
        sc.dirty = b_size;
        if TB {
            sc.trace.resize(b_size, SCRIPT_GAP_IN_A);
            sc.rows.push((0, 0, b_size as u32));
        }

        let mut ramp = [0i16; 16];
        for (i, r) in ramp.iter_mut().enumerate() {
            *r = (i as i32 * gap_extend) as i16;
        }
        let c = C16 {
            goe: _mm256_set1_epi16(goe as i16),
            ge1: _mm256_set1_epi16(gap_extend as i16),
            ge16: _mm256_set1_epi16((16 * gap_extend) as i16),
            ramp: ld(ramp.as_ptr()),
            xv: _mm256_set1_epi16(x as i16),
            kb7: ld(KB7X2.as_ptr()),
            c3: _mm256_set1_epi16(SCRIPT_SUB as i16),
            c10: _mm256_set1_epi16(SCRIPT_EXTEND_GAP_A as i16),
            c40: _mm256_set1_epi16(SCRIPT_EXTEND_GAP_B as i16),
        };
        let zero = _mm256_setzero_si256();
        let lane0 = _mm256_setr_epi16(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        let k70 = _mm256_set1_epi8(0x70);
        let k10 = _mm256_set1_epi8(0x10);
        let static_tabs: Option<*const u8> = match scores {
            Scores::Tables { t, .. } => Some(t.data.as_ptr()),
            _ => None,
        };
        let ncols = scores.ncols();
        let mut cells: u64 = 0;

        let mut base: i32 = 0;
        let mut best_rel: i32 = 0;
        let mut a_off = 0usize;
        let mut b_off = 0usize;
        let mut first = 0usize;
        let mut fence_hit = false;
        // exclusive end of the lanes last written in each H buffer
        let mut wend = [b_size, 0usize];
        let mut cur = 0usize; // buffer holding the previous row
        // V = 10: one past the last row's last survivor (row 0: its written lanes)
        let mut lastp1 = b_size;

        for a in 1..=len1 {
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
                (sc.ha.as_mut_ptr().add(PAD), sc.hb.as_mut_ptr().add(PAD))
            } else {
                (sc.hb.as_mut_ptr().add(PAD), sc.ha.as_mut_ptr().add(PAD))
            };
            let e = sc.e.as_mut_ptr().add(PAD);
            let sb = sc.sb.as_ptr().add(PAD);
            let t0lo = bcast128(t);
            let t1lo = bcast128(t.add(16));
            let t0hi = bcast128(t.add(32));
            let t1hi = bcast128(t.add(48));

            if best_rel > limit_rebase {
                let dv = _mm256_set1_epi16(best_rel as i16);
                let mut b = first;
                while b < b_size {
                    st(hr.add(b), _mm256_subs_epu16(ld(hr.add(b)), dv));
                    st(e.add(b), _mm256_subs_epu16(ld(e.add(b)), dv));
                    b += 16;
                }
                base += best_rel;
                best_rel = 0;
            }
            // no diagonal into the first lane of the band: lane 0 of the row's first block (`hpm`);
            // a store of 0 into the H row would sit right under the first block's load
            let mut hpm = lane0;

            let row_off = sc.trace.len();
            if TB {
                sc.trace.reserve(bend - first + num_extra + 96);
            }
            let tr = sc.trace.as_mut_ptr().add(row_off);

            let bestv0 = _mm256_set1_epi16((bias + best_rel) as i16);
            let mut rs = Row {
                carry: zero,
                bestv: bestv0,
                thrv: _mm256_subs_epu16(bestv0, c.xv),
                rec_b: 0,
                slow: false,
            };
            let mut lo_surv = usize::MAX;
            let mut hi_surv = 0usize;
            let mut b0 = first;
            // Two blocks (32 lanes) while more than 16 lanes of the band are left, then single blocks
            // for the rest of the band and the row-gap extension. The stop test runs on the grid
            // first + 16k once the band is passed; lanes processed past the point where the 128-bit
            // kernel stops are computed exactly and dropped (row gap below the threshold), so the
            // result is the same.
            while b0 + 16 < bend {
                // scores of lanes b0..b0+32 (lane b reads column b-1)
                let idx = _mm256_permute4x64_epi64::<0xD8>(ld(sb.add(b0)));
                let i0 = _mm256_add_epi8(idx, k70);
                let i1 = _mm256_sub_epi8(idx, k10);
                let lo = lut32(t0lo, t1lo, i0, i1);
                let hi = lut32(t0hi, t1hi, i0, i1);
                // the first block lies inside the band (b0 + 16 < bend <= limit); the second may
                // reach past the band and past `limit`
                let mask1 = if b0 + 32 <= limit {
                    None
                } else {
                    Some(ld(LANE_MASK16.as_ptr().add(32 - (limit - b0))))
                };
                let (d0, r0) = block16_p::<TB, V>(
                    &c,
                    &mut rs,
                    hr,
                    hw,
                    e,
                    b0,
                    _mm256_unpacklo_epi8(lo, hi),
                    None,
                    hpm,
                    em16::<V>(lastp1, b0),
                );
                hpm = zero;
                let (d1, r1) = block16_p::<TB, V>(
                    &c,
                    &mut rs,
                    hr,
                    hw,
                    e,
                    b0 + 16,
                    _mm256_unpackhi_epi8(lo, hi),
                    mask1,
                    zero,
                    em16::<V>(lastp1, b0 + 16),
                );
                if TB {
                    st(
                        tr.add(b0 - first),
                        _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi16(r0, r1)),
                    );
                }
                // 2 bits per lane
                let live = !((_mm256_movemask_epi8(d0) as u32 as u64)
                    | ((_mm256_movemask_epi8(d1) as u32 as u64) << 32));
                if live != 0 {
                    if lo_surv == usize::MAX {
                        lo_surv = b0 + (live.trailing_zeros() as usize >> 1);
                    }
                    hi_surv = b0 + ((63 - live.leading_zeros() as usize) >> 1);
                }
                b0 += 32;
            }
            loop {
                if b0 >= bend {
                    if b0 >= limit {
                        break;
                    }
                    if (_mm256_cvtsi256_si32(rs.carry) as i16) < (_mm256_cvtsi256_si32(rs.thrv) as i16)
                    {
                        break;
                    }
                }
                let idx = _mm256_permute4x64_epi64::<0xD8>(ld(sb.add(b0)));
                let i0 = _mm256_add_epi8(idx, k70);
                let i1 = _mm256_sub_epi8(idx, k10);
                let lo = lut32(t0lo, t1lo, i0, i1);
                let hi = lut32(t0hi, t1hi, i0, i1);
                let mask = if b0 + 16 <= limit {
                    None
                } else {
                    Some(ld(LANE_MASK16.as_ptr().add(16 - (limit - b0))))
                };
                let (d0, r0) = block16_p::<TB, V>(
                    &c,
                    &mut rs,
                    hr,
                    hw,
                    e,
                    b0,
                    _mm256_unpacklo_epi8(lo, hi),
                    mask,
                    hpm,
                    em16::<V>(lastp1, b0),
                );
                hpm = zero;
                if TB {
                    stx(
                        tr.add(b0 - first),
                        _mm256_castsi256_si128(_mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi16(
                            r0, r0,
                        ))),
                    );
                }
                let live = !(_mm256_movemask_epi8(d0) as u32);
                if live != 0 {
                    if lo_surv == usize::MAX {
                        lo_surv = b0 + (live.trailing_zeros() as usize >> 1);
                    }
                    hi_surv = b0 + ((31 - live.leading_zeros() as usize) >> 1);
                }
                b0 += 16;
            }
            let proc_end = b0;
            if proc_end > sc.dirty {
                sc.dirty = proc_end;
            }
            // lanes this buffer held beyond what the row wrote are stale
            let wi = 1 - cur;
            let mut b = proc_end;
            while b < wend[wi] {
                st(hw.add(b), zero);
                b += 16;
            }
            wend[wi] = proc_end;
            cur = wi;
            if rs.slow {
                b_off = rs.rec_b;
                a_off = a;
                best_rel = (_mm256_cvtsi256_si32(rs.bestv) as i16 as i32) - bias;
            }
            if TB {
                sc.trace.set_len(row_off + proc_end - first);
                sc.rows
                    .push((row_off as u32, first as u32, (proc_end - first) as u32));
            }
            if fence_row {
                fence_hit = true;
                break;
            }
            if lo_surv == usize::MAX {
                break;
            }
            // columns right of the last survivor leave the band: E = -inf
            let last = hi_surv;
            if V != 9 && V != 10 {
                let mut b = last + 1;
                while b < proc_end {
                    st(e.add(b), zero);
                    b += 16;
                }
            }
            lastp1 = last + 1;
            first = lo_surv;
            b_size = last + 1;
            if b_size <= len2 {
                b_size += 1;
            }
        }

        let best = base + best_rel;
        finish::<TB>(sc, a_off, b_off, best, fence_hit, cells)
    }
}
