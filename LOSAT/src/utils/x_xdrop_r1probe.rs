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

/// One call as `xdrop_align` makes it with the 256-bit kernels on: `x_avx2::align`.
#[inline(never)]
fn call(p: &Prob<'_>, tb: bool, sc: &mut XdropScratch) -> Option<XdropResult> {
    // SAFETY: `available()` is checked before any call.
    unsafe {
        align(
            &p.q, &p.s, &p.scores, p.len1, p.len2, p.go, p.ge, p.x, tb, p.cf, sc,
        )
    }
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
    for p in &list {
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
        if call(p, true, &mut sd).is_some() {
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
        for p in &list {
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
         lay={lay} | {}\n",
        synth.as_deref().unwrap_or("-"),
        list.len(),
        if rows_total > 0.0 {
            loop_s * 1e9 / rows_total
        } else {
            0.0
        },
        if spread_on { "spread" } else { "as-is" },
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
