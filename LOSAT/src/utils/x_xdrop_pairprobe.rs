//! EXPERIMENT (line C, C3 premise probe; branch exp/c3-pairlane-probe, not for merging). No NCBI
//! counterpart: a measuring harness only; it is compiled into the test build only and does not
//! change the library or the LOSAT binary.
//!
//! A pair-lane form of the `ALIGN_EX` traceback DP (`v4-C/C3_ALIGNEX_NOTES.md`): 16 independent
//! captured problems, one per 16-bit lane, in lockstep over rows a (A, the query) and columns b
//! (B, the subject). Each row processes the union of the lanes' bands; inside a lane the cells run
//! in NCBI's order (row-major), so the running best, the X-drop test and the first best cell are
//! NCBI's per lane. The values are those of the row kernels (`score - base + bias`, 0 = minus
//! infinity, per-lane base and rebase), and the argument of `xdrop_simd.rs` (whatever differs below
//! `thr` cannot change a surviving cell, a drop or a script bit the walk reads) covers the lanes
//! outside their band: a lane's H is stored as 0 wherever it drops, every column a later row reads
//! was written by the row before (the union of the next row lies inside this row's processed
//! columns plus one, and a survivor's column is inside its lane's next band), and a lane that ends
//! at its last row gets its H and E cleared, so a dead lane only ever sees values below its
//! threshold. Past a lane's subject end (b > N) the lane is masked in rows that can reach it.
//!
//! The score of lane l at (a, b) is M_l[A_l(a)][B_l(b - 1)]: two 8-lane `vpgatherdd` per step from
//! per-lane 16-bit matrices (`gather`), or the same vectors recorded in a first pass and read back
//! in order (`pre`: the cost without the gather).
//!
//! Variables: `P_DIR` (capture directory), `P_KERN` (row = `x_avx2::align`, the row kernel; tb,
//! so, tb2 = score-only pass then traceback pass to each lane's best row, batched again by that
//! row; tbpre, sopre, tb2pre = the same with recorded scores), `P_ORDER` (stream, m = len1, rows =
//! rows run (oracle), cells = cells (oracle)), `P_REPS` (default 0), `P_LIMIT`, `P_OUT`.

use super::super::x_capture::{load_dir, Captured};
use super::super::{
    fits_8bit, XdropResult, XdropScratch, FENCE_SENTRY, MAX_M, SCRIPT_EXTEND_GAP_A,
    SCRIPT_EXTEND_GAP_B, SCRIPT_GAP_IN_A, SCRIPT_GAP_IN_B, SCRIPT_OP_MASK, SCRIPT_SUB, VMAX,
};
use super::{align, available};
use core::arch::x86_64::*;
use std::time::Instant;

type W = __m256i;

const L: usize = 16;
/// Lane stride of the 16-bit matrices (28 x 28 entries and padding).
const LS: usize = 28 * 28 + 16;

fn var(k: &str) -> Option<String> {
    std::env::var(k).ok().filter(|v| !v.is_empty())
}

#[inline(always)]
unsafe fn ld<T>(p: *const T) -> W {
    _mm256_loadu_si256(p as *const W)
}
#[inline(always)]
unsafe fn st<T>(p: *mut T, v: W) {
    _mm256_storeu_si256(p as *mut W, v)
}
#[inline(always)]
unsafe fn hmin_u16(v: W) -> u16 {
    let m = _mm_min_epu16(_mm256_castsi256_si128(v), _mm256_extracti128_si256::<1>(v));
    _mm_cvtsi128_si32(_mm_minpos_epu16(m)) as u16
}
#[inline(always)]
unsafe fn hmax_u16(v: W) -> u16 {
    !hmin_u16(_mm256_xor_si256(v, _mm256_set1_epi16(-1)))
}
/// Gather slot of lane p: lanes 0, 2, .., 14 in the first index vector, 1, 3, .., 15 in the second.
#[inline(always)]
fn slot(p: usize) -> usize {
    (p & 1) * 8 + (p >> 1)
}
fn lanes_i16(v: W) -> [i16; 16] {
    let mut o = [0i16; 16];
    // SAFETY: 32-byte store into a 32-byte array.
    unsafe { st(o.as_mut_ptr(), v) };
    o
}

/// Result of one lane: the `XdropResult` fields and the ops in walk order.
#[derive(Clone, Default, Debug, PartialEq)]
struct Out {
    a_off: usize,
    b_off: usize,
    score: i32,
    cells: u64,
    ops: Vec<(u8, u32)>,
}

#[derive(Default)]
struct Count {
    batches: u64,
    steps: u64,
    rows: u64,
    kill_rows: u64,
    rebases: u64,
    clears: u64,
}

/// Scratch of one batch, reused.
struct PB {
    mats: Vec<i16>,
    cidx: Vec<i32>,
    cfill: usize,
    h: Vec<i16>,
    e: Vec<i16>,
    hw: usize,
    trace: Vec<u8>,
    /// per row a >= 1: (offset in `trace`, first column, columns)
    rinfo: Vec<(usize, usize, usize)>,
    pre: Vec<W>,
    pre_k: usize,
}

impl PB {
    fn new() -> Self {
        PB {
            mats: vec![0; L * LS + 16],
            cidx: Vec::new(),
            cfill: 0,
            h: Vec::new(),
            e: Vec::new(),
            hw: 0,
            trace: Vec::new(),
            rinfo: Vec::new(),
            pre: Vec::new(),
            pre_k: 0,
        }
    }
}

struct K {
    goe: W,
    ge: W,
    xv: W,
    one: W,
    c3: W,
    c10: W,
    c40: W,
}

/// Raw pointers of the batch scratch for one row (the vectors do not move inside a stretch of
/// steps).
#[derive(Clone, Copy)]
struct Ptrs {
    h: *mut i16,
    e: *mut i16,
    cidx: *const i32,
    mats: *const i32,
}

/// Per-lane state carried along a row.
struct St {
    f: W,
    hpp: W,
    bv: W,
    nlim: W,
    bestv: W,
    boffv: W,
    firstv: W,
    lastv: W,
    rb0: W,
    rb1: W,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:578-635 (one cell, per
/// lane; the lane operations of `x_xdrop_avx2.rs` `block16`, with the row gap sequential inside the
/// lane and the running best per lane).
#[inline(always)]
unsafe fn step<const TB: bool, const FETCH: u8, const KILL: bool>(
    c: &K,
    s_: &mut St,
    pt: Ptrs,
    cur: &mut *mut W,
    b: usize,
) -> W {
    let h = pt.h.add(b * L);
    let e = pt.e.add(b * L);
    let s = if FETCH == 1 {
        let v = ld(*cur);
        *cur = cur.add(1);
        v
    } else {
        let ci = pt.cidx.add(b * L);
        let i0 = _mm256_add_epi32(s_.rb0, ld(ci));
        let i1 = _mm256_add_epi32(s_.rb1, ld(ci.add(8)));
        let g0 = _mm256_i32gather_epi32::<2>(pt.mats, i0);
        let g1 = _mm256_i32gather_epi32::<2>(pt.mats, i1);
        let v = _mm256_blend_epi16::<0xAA>(g0, _mm256_slli_epi32::<16>(g1));
        if FETCH == 2 {
            st(*cur, v);
            *cur = cur.add(1);
        }
        v
    };
    let hp = ld(h);
    let ep = ld(e);
    let d = _mm256_adds_epi16(s_.hpp, s);
    let h0 = _mm256_max_epi16(d, ep);
    let mut hh = _mm256_max_epi16(h0, s_.f);
    if KILL {
        hh = _mm256_andnot_si256(_mm256_cmpgt_epi16(s_.bv, s_.nlim), hh);
    }
    let thr = _mm256_subs_epu16(s_.bestv, c.xv);
    let drop = _mm256_cmpgt_epi16(thr, hh);
    let gt = _mm256_cmpgt_epi16(hh, s_.bestv);
    s_.bestv = _mm256_max_epi16(s_.bestv, hh);
    s_.boffv = _mm256_blendv_epi8(s_.boffv, s_.bv, gt);
    st(h, _mm256_andnot_si256(drop, hh));
    let hg = _mm256_subs_epu16(hh, c.goe);
    let es = _mm256_subs_epu16(ep, c.ge);
    st(e, _mm256_max_epi16(hg, es));
    let fs = _mm256_subs_epu16(s_.f, c.ge);
    let scr = if TB {
        let me = _mm256_cmpgt_epi16(ep, d);
        let mf = _mm256_cmpgt_epi16(s_.f, h0);
        let op = _mm256_andnot_si256(mf, _mm256_add_epi16(c.c3, _mm256_and_si256(me, c.c3)));
        let fb = _mm256_andnot_si256(_mm256_cmpgt_epi16(hg, es), c.c40);
        let fa = _mm256_andnot_si256(_mm256_cmpgt_epi16(hg, fs), c.c10);
        _mm256_or_si256(op, _mm256_or_si256(fb, fa))
    } else {
        _mm256_setzero_si256()
    };
    s_.f = _mm256_max_epi16(fs, hg);
    let x = _mm256_or_si256(s_.bv, drop);
    s_.firstv = _mm256_min_epu16(s_.firstv, x);
    s_.lastv = _mm256_max_epi16(s_.lastv, x);
    s_.bv = _mm256_add_epi16(s_.bv, c.one);
    s_.hpp = hp;
    scr
}

/// The steps of one row, two columns at a time: the union band, then the row-gap extension while
/// a live lane's row gap reaches its threshold (NCBI :652-664; past a lane's band the full
/// recurrence gives H = F, E = F - gap_open_extend). Returns the end of the processed columns.
#[inline(always)]
unsafe fn pairs<const TB: bool, const FETCH: u8, const KILL: bool>(
    c: &K,
    s_: &mut St,
    pt: Ptrs,
    cur: &mut *mut W,
    tr: *mut u8,
    ufirst: usize,
    ubend: usize,
    alive: W,
) -> usize {
    let mut b = ufirst;
    while b < ubend {
        let r0 = step::<TB, FETCH, KILL>(c, s_, pt, cur, b);
        let r1 = step::<TB, FETCH, KILL>(c, s_, pt, cur, b + 1);
        if TB {
            st(
                tr.add((b - ufirst) * L),
                _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi16(r0, r1)),
            );
        }
        b += 2;
    }
    loop {
        let thr = _mm256_subs_epu16(s_.bestv, c.xv);
        let ok = _mm256_andnot_si256(
            _mm256_or_si256(_mm256_cmpgt_epi16(thr, s_.f), _mm256_cmpgt_epi16(s_.bv, s_.nlim)),
            alive,
        );
        if _mm256_movemask_epi8(ok) == 0 {
            return b;
        }
        let r0 = step::<TB, FETCH, KILL>(c, s_, pt, cur, b);
        let r1 = step::<TB, FETCH, KILL>(c, s_, pt, cur, b + 1);
        if TB {
            st(
                tr.add((b - ufirst) * L),
                _mm256_permute4x64_epi64::<0xD8>(_mm256_packus_epi16(r0, r1)),
            );
        }
        b += 2;
    }
}

/// Column letters of the batch's lanes into `cidx` up to column `to` (exclusive): lane p's
/// column b reads B[b - 1] (`cols[b - 1]`, 0 past the end), as the row kernels' lane b.
unsafe fn fill(pb: &mut PB, lanes: &[&Captured], to: usize) {
    if pb.cidx.len() < to * L {
        pb.cidx.resize(to * L + 64 * L, 0);
    }
    for b in pb.cfill..to {
        for p in 0..L {
            let ch = match lanes.get(p) {
                Some(q) if b >= 1 && b - 1 < q.cols.len() => q.cols[b - 1] as usize,
                _ => 0,
            };
            *pb.cidx.get_unchecked_mut(b * L + slot(p)) = (p * LS + ch) as i32;
        }
    }
    if to > pb.cfill {
        pb.cfill = to;
    }
}

/// One batch of up to 16 problems; `mlim[p]` = rows to run in lane p (len1, or the best row of the
/// score-only pass). Returns the lanes' results (cells only meaningful for full runs).
#[target_feature(enable = "avx2,bmi1,lzcnt,popcnt")]
unsafe fn batch<const TB: bool, const FETCH: u8>(
    lanes: &[&Captured],
    mlim: &[usize],
    pb: &mut PB,
    cnt: &mut Count,
    out: &mut [Out],
) {
    let np = lanes.len();
    let q0 = lanes[0];
    let goe = q0.gap_open + q0.gap_extend;
    let x = q0.x_drop.max(goe);
    let bias = x + goe + MAX_M + 1;
    let limv = _mm256_set1_epi16((VMAX as i32 - MAX_M) as i16);
    let maxn = lanes.iter().map(|q| q.len2).max().unwrap();
    // num_extra (NCBI :458) and slack: the row-gap extension and the pair step stay inside it
    let margin = (x / q0.gap_extend + 3) as usize + 8;
    let need = (maxn + margin + 8) * L;
    if pb.h.len() < need {
        pb.h.resize(need, 0);
        pb.e.resize(need, 0);
    }
    // clear what the last batch wrote
    let z = _mm256_setzero_si256();
    for b in 0..pb.hw.min(pb.h.len() / L) {
        st(pb.h.as_mut_ptr().add(b * L), z);
        st(pb.e.as_mut_ptr().add(b * L), z);
    }
    pb.hw = 0;
    pb.cfill = 0;
    pb.trace.clear();
    pb.rinfo.clear();
    pb.rinfo.push((0, 0, 0));
    // per-lane 16-bit matrices (i32 -> i16 saturating, as `build_row`)
    for (p, q) in lanes.iter().enumerate() {
        let src = q.matrix.as_ptr();
        let dst = pb.mats.as_mut_ptr().add(p * LS);
        let mut k = 0;
        while k + 16 <= 784 {
            let w = _mm256_permute4x64_epi64::<0xD8>(_mm256_packs_epi32(
                ld(src.add(k)),
                ld(src.add(k + 8)),
            ));
            st(dst.add(k), w);
            k += 16;
        }
    }
    for p in np..L {
        for k in 0..784 {
            pb.mats[p * LS + k] = 0;
        }
    }
    // row 0 per lane (NCBI :476-492)
    let mut bs0 = [0i16; 16];
    let mut nv_a = [-1i16; 16];
    let mut mv_a = [0i16; 16];
    let mut b0max = 1usize;
    for (p, q) in lanes.iter().enumerate() {
        let h = pb.h.as_mut_ptr();
        let e = pb.e.as_mut_ptr();
        *h.add(p) = bias as i16;
        *e.add(p) = (bias - goe) as i16;
        let mut score = -goe;
        let mut bs = 1usize;
        for i in 1..=q.len2 {
            if score < -x {
                break;
            }
            *h.add(i * L + p) = (bias + score) as i16;
            *e.add(i * L + p) = (bias + score - goe) as i16;
            score -= q.gap_extend;
            bs = i + 1;
        }
        *h.add((bs - 1) * L + p) = 0;
        bs0[p] = bs as i16;
        nv_a[p] = q.len2 as i16;
        mv_a[p] = mlim[p] as i16;
        b0max = b0max.max(bs);
    }
    pb.hw = b0max;
    let c = K {
        goe: _mm256_set1_epi16(goe as i16),
        ge: _mm256_set1_epi16(q0.gap_extend as i16),
        xv: _mm256_set1_epi16(x as i16),
        one: _mm256_set1_epi16(1),
        c3: _mm256_set1_epi16(SCRIPT_SUB as i16),
        c10: _mm256_set1_epi16(SCRIPT_EXTEND_GAP_A as i16),
        c40: _mm256_set1_epi16(SCRIPT_EXTEND_GAP_B as i16),
    };
    let ones = _mm256_set1_epi16(-1);
    let nv = ld(nv_a.as_ptr());
    let mv = ld(mv_a.as_ptr());
    // lanes with mlim 0 (and unused lanes) are dead from the start
    let mut nlim = _mm256_or_si256(nv, _mm256_cmpgt_epi16(_mm256_set1_epi16(1), mv));
    let mut bandfirst = z;
    let mut bsizev = ld(bs0.as_ptr());
    let mut bestv = _mm256_set1_epi16(bias as i16);
    let mut aoffv = z;
    let mut boffv = z;
    let mut cells_lo = z;
    let mut cells_hi = z;
    let mut base = [0i32; 16];
    let mut a = 1usize;
    loop {
        let alive = _mm256_cmpgt_epi16(nlim, ones);
        if _mm256_movemask_epi8(alive) == 0 {
            break;
        }
        let ufirst = hmin_u16(_mm256_or_si256(bandfirst, _mm256_xor_si256(alive, ones))) as usize;
        let ubend = hmax_u16(_mm256_and_si256(bsizev, alive)) as usize;
        let minn = hmin_u16(_mm256_or_si256(nv, _mm256_xor_si256(alive, ones))) as usize;
        let dl = _mm256_and_si256(_mm256_sub_epi16(bsizev, bandfirst), alive);
        cells_lo = _mm256_add_epi32(cells_lo, _mm256_cvtepu16_epi32(_mm256_castsi256_si128(dl)));
        cells_hi =
            _mm256_add_epi32(cells_hi, _mm256_cvtepu16_epi32(_mm256_extracti128_si256::<1>(dl)));
        // rebase (per lane, as the row kernel at the row start)
        let rbm = _mm256_and_si256(_mm256_cmpgt_epi16(bestv, limv), alive);
        if _mm256_movemask_epi8(rbm) != 0 {
            let dv = _mm256_and_si256(rbm, _mm256_subs_epu16(bestv, _mm256_set1_epi16(bias as i16)));
            for b in ufirst..pb.hw {
                let hp = pb.h.as_mut_ptr().add(b * L);
                let ep = pb.e.as_mut_ptr().add(b * L);
                st(hp, _mm256_subs_epu16(ld(hp), dv));
                st(ep, _mm256_subs_epu16(ld(ep), dv));
            }
            bestv = _mm256_subs_epu16(bestv, dv);
            let d = lanes_i16(dv);
            for p in 0..L {
                base[p] += d[p] as u16 as i32;
            }
            cnt.rebases += 1;
        }
        // row letters: lane p reads matrix row A_p(a) (letter 0 past its last row)
        let mut rb = [0i32; 16];
        for (p, q) in lanes.iter().enumerate() {
            if a <= q.len1 {
                rb[slot(p)] = q.rows[a - 1] as i32 * 28;
            }
        }
        let tr_end = ubend + margin;
        if FETCH != 1 && pb.cfill < tr_end + 2 {
            fill(pb, lanes, tr_end + 64);
        }
        if pb.h.len() < (tr_end + 2) * L {
            pb.h.resize((tr_end + 66) * L, 0);
            pb.e.resize((tr_end + 66) * L, 0);
        }
        let row_off = pb.trace.len();
        if TB {
            pb.trace.reserve((tr_end + 2 - ufirst) * L);
        }
        let mut s_ = St {
            f: z,
            hpp: z,
            bv: _mm256_set1_epi16(ufirst as i16),
            nlim,
            bestv,
            boffv,
            firstv: ones,
            lastv: ones,
            rb0: ld(rb.as_ptr()),
            rb1: ld(rb.as_ptr().add(8)),
        };
        let kill = ubend + 60 > minn;
        cnt.rows += 1;
        if kill {
            cnt.kill_rows += 1;
        }
        let pt = Ptrs {
            h: pb.h.as_mut_ptr(),
            e: pb.e.as_mut_ptr(),
            cidx: pb.cidx.as_ptr(),
            mats: pb.mats.as_ptr() as *const i32,
        };
        let tr = pb.trace.as_mut_ptr().add(row_off);
        if FETCH == 2 {
            pb.pre.reserve(tr_end + 2 - ufirst);
        }
        let mut cur: *mut W = pb.pre.as_mut_ptr().add(pb.pre_k);
        let b = if kill {
            pairs::<TB, FETCH, true>(&c, &mut s_, pt, &mut cur, tr, ufirst, ubend, alive)
        } else {
            pairs::<TB, FETCH, false>(&c, &mut s_, pt, &mut cur, tr, ufirst, ubend, alive)
        };
        // the row-gap extension stays within num_extra columns of the band end
        assert!(b + 2 <= tr_end, "row ran past the reserved columns");
        let proc_end = b;
        cnt.steps += (proc_end - ufirst) as u64;
        pb.pre_k = cur.offset_from(pb.pre.as_ptr()) as usize;
        if FETCH == 2 {
            pb.pre.set_len(pb.pre_k);
        }
        if proc_end > pb.hw {
            pb.hw = proc_end;
        }
        if TB {
            pb.trace.set_len(row_off + (proc_end - ufirst) * L);
        }
        pb.rinfo.push((row_off, ufirst, proc_end - ufirst));
        boffv = s_.boffv;
        aoffv = _mm256_blendv_epi8(aoffv, _mm256_set1_epi16(a as i16), _mm256_cmpgt_epi16(s_.bestv, bestv));
        bestv = s_.bestv;
        // band of the next row (NCBI :638-675): [first survivor, last survivor + 1 (+1 sentinel))
        let dead = _mm256_cmpeq_epi16(s_.firstv, ones);
        bandfirst = s_.firstv;
        bsizev = _mm256_sub_epi16(_mm256_add_epi16(s_.lastv, c.one), _mm256_cmpgt_epi16(nv, s_.lastv));
        // lanes that ran their last row: clear their H and E so they stay below the threshold
        let fin = _mm256_andnot_si256(dead, _mm256_and_si256(alive, _mm256_cmpgt_epi16(_mm256_set1_epi16((a + 1) as i16), mv)));
        if _mm256_movemask_epi8(fin) != 0 {
            for b in ufirst..pb.hw {
                let hp = pb.h.as_mut_ptr().add(b * L);
                let ep = pb.e.as_mut_ptr().add(b * L);
                st(hp, _mm256_andnot_si256(fin, ld(hp)));
                st(ep, _mm256_andnot_si256(fin, ld(ep)));
            }
            cnt.clears += 1;
        }
        nlim = _mm256_or_si256(nlim, _mm256_or_si256(dead, fin));
        a += 1;
    }
    cnt.batches += 1;
    // per-lane results and the walk (NCBI :682-727, `kernel::finish`)
    let best = lanes_i16(bestv);
    let ao = lanes_i16(aoffv);
    let bo = lanes_i16(boffv);
    let mut cl = [0i32; 16];
    st(cl.as_mut_ptr(), cells_lo);
    st(cl.as_mut_ptr().add(8), cells_hi);
    for p in 0..np {
        let o = &mut out[p];
        o.ops.clear();
        o.cells = cl[p] as u32 as u64;
        let score = base[p] + best[p] as u16 as i32 - bias;
        if score <= 0 {
            o.a_off = 0;
            o.b_off = 0;
            o.score = 0;
            continue;
        }
        o.a_off = ao[p] as u16 as usize;
        o.b_off = bo[p] as u16 as usize;
        o.score = score;
        if !TB {
            continue;
        }
        let (mut ai, mut bi) = (o.a_off, o.b_off);
        let mut script = SCRIPT_SUB;
        while ai > 0 || bi > 0 {
            let next = if ai == 0 {
                SCRIPT_GAP_IN_A
            } else {
                let (off, start, len) = pb.rinfo[ai];
                if bi < start || bi - start >= len {
                    break;
                }
                *pb.trace.get_unchecked(off + (bi - start) * L + p)
            };
            match script & SCRIPT_OP_MASK {
                v if v == SCRIPT_GAP_IN_A => {
                    script = next & SCRIPT_OP_MASK;
                    if next & SCRIPT_EXTEND_GAP_A != 0 {
                        script = SCRIPT_GAP_IN_A;
                    }
                }
                v if v == SCRIPT_GAP_IN_B => {
                    script = next & SCRIPT_OP_MASK;
                    if next & SCRIPT_EXTEND_GAP_B != 0 {
                        script = SCRIPT_GAP_IN_B;
                    }
                }
                _ => script = next & SCRIPT_OP_MASK,
            }
            let op = script & SCRIPT_OP_MASK;
            if op == SCRIPT_GAP_IN_A {
                bi = bi.saturating_sub(1);
            } else if op == SCRIPT_GAP_IN_B {
                ai = ai.saturating_sub(1);
            } else {
                ai = ai.saturating_sub(1);
                bi = bi.saturating_sub(1);
            }
            match o.ops.last_mut() {
                Some(l) if l.0 == op => l.1 += 1,
                _ => o.ops.push((op, 1)),
            }
        }
    }
}

/// The problems this probe takes: 16-bit traceback problems with a recorded result, 28-letter
/// matrices with entries <= MAX_M, the redo gap costs of the first one, no fence or bad letter.
fn usable(p: &Captured, first: Option<&Captured>) -> bool {
    if !p.tb || p.result.is_none() || p.n != 28 || p.len1 == 0 || p.len2 == 0 {
        return false;
    }
    if fits_8bit(&p.scores(), p.gap_open, p.gap_extend, p.x_drop).is_some() {
        return false;
    }
    if p.len1 > 30000 || p.len2 > 30000 || p.matrix.iter().any(|&v| v > MAX_M) {
        return false;
    }
    if p.rows.iter().any(|&r| r >= 28) || p.cols.iter().any(|&c| c >= 28 || c == FENCE_SENTRY) {
        return false;
    }
    if let Some(f) = first {
        if (p.gap_open, p.gap_extend, p.x_drop) != (f.gap_open, f.gap_extend, f.x_drop) {
            return false;
        }
    }
    true
}

fn to_res(o: &Out) -> XdropResult {
    XdropResult {
        a_offset: o.a_off,
        b_offset: o.b_off,
        score: o.score,
        fence_hit: false,
        cells: o.cells,
    }
}

/// All batches of one kernel over `order` (problem indices); returns per-problem results.
unsafe fn run_all<const TB: bool, const FETCH: u8>(
    list: &[&Captured],
    order: &[usize],
    two_pass: bool,
    pb: &mut PB,
    cnt: &mut Count,
    res: &mut [Out],
    cnt2: &mut Count,
) {
    let mut out: Vec<Out> = vec![Out::default(); L];
    for chunk in order.chunks(L) {
        let lanes: Vec<&Captured> = chunk.iter().map(|&i| list[i]).collect();
        let ml: Vec<usize> = lanes.iter().map(|q| q.len1).collect();
        if two_pass {
            batch::<false, FETCH>(&lanes, &ml, pb, cnt, &mut out);
        } else {
            batch::<TB, FETCH>(&lanes, &ml, pb, cnt, &mut out);
        }
        for (k, &i) in chunk.iter().enumerate() {
            res[i] = out[k].clone();
        }
    }
    if !two_pass {
        return;
    }
    // pass 2: rows 1..=a_off of every lane with a positive score, batched by a_off
    let mut o2: Vec<usize> = order.iter().copied().filter(|&i| res[i].score > 0).collect();
    o2.sort_by_key(|&i| std::cmp::Reverse(res[i].a_off));
    for chunk in o2.chunks(L) {
        let lanes: Vec<&Captured> = chunk.iter().map(|&i| list[i]).collect();
        let ml: Vec<usize> = chunk.iter().map(|&i| res[i].a_off).collect();
        batch::<true, FETCH>(&lanes, &ml, pb, cnt2, &mut out);
        for (k, &i) in chunk.iter().enumerate() {
            let r = &mut res[i];
            debug_assert_eq!((out[k].a_off, out[k].b_off, out[k].score), (r.a_off, r.b_off, r.score));
            if (out[k].a_off, out[k].b_off, out[k].score) != (r.a_off, r.b_off, r.score) {
                panic!("pass 2 differs from pass 1 at problem {i}");
            }
            r.ops = std::mem::take(&mut out[k].ops);
        }
    }
}

#[test]
#[ignore = "needs P_DIR=<capture directory>"]
fn pair_loop() {
    let Some(dir) = var("P_DIR") else {
        eprintln!("P_DIR is not set");
        return;
    };
    assert!(available(), "the probe needs AVX2");
    let kern = var("P_KERN").unwrap_or_else(|| "tb".into());
    let ord = var("P_ORDER").unwrap_or_else(|| "rows".into());
    let reps: usize = var("P_REPS").and_then(|v| v.parse().ok()).unwrap_or(0);
    let limit: usize = var("P_LIMIT").and_then(|v| v.parse().ok()).unwrap_or(usize::MAX);
    let t_load = Instant::now();
    let probs = load_dir(std::path::Path::new(&dir)).expect("read the capture files");
    let mut list: Vec<&Captured> = Vec::new();
    let mut first: Option<&Captured> = None;
    let mut skipped = 0usize;
    for p in &probs {
        if list.len() >= limit {
            break;
        }
        if !p.tb || p.result.is_none() || fits_8bit(&p.scores(), p.gap_open, p.gap_extend, p.x_drop).is_some() {
            continue;
        }
        if usable(p, first) {
            if first.is_none() {
                first = Some(p);
            }
            list.push(p);
        } else {
            skipped += 1;
        }
    }
    let load_s = t_load.elapsed().as_secs_f64();
    // the row kernel: results, ops and rows run (the oracle order)
    let mut sc = XdropScratch::new();
    let mut rows_run = vec![0usize; list.len()];
    let mut want: Vec<Out> = Vec::with_capacity(list.len());
    let mut cells_total = 0u64;
    for (i, p) in list.iter().enumerate() {
        // SAFETY: AVX2 checked above.
        let r = unsafe {
            align(
                &p.row_seq(), &p.col_seq(), &p.scores(), p.len1, p.len2, p.gap_open, p.gap_extend,
                p.x_drop, true, p.check_fence, &mut sc,
            )
        }
        .expect("row kernel declined");
        assert_eq!(Some(r), p.result, "row kernel vs recorded");
        assert_eq!(sc.ops, p.ops, "row kernel ops vs recorded");
        rows_run[i] = sc.rows.len() - 1;
        cells_total += r.cells;
        want.push(Out {
            a_off: r.a_offset,
            b_off: r.b_offset,
            score: r.score,
            cells: r.cells,
            ops: p.ops.clone(),
        });
    }
    let mut order: Vec<usize> = (0..list.len()).collect();
    match ord.as_str() {
        "stream" => {}
        "m" => order.sort_by_key(|&i| std::cmp::Reverse(list[i].len1)),
        "rows" => order.sort_by_key(|&i| std::cmp::Reverse(rows_run[i])),
        "cells" => order.sort_by_key(|&i| std::cmp::Reverse(want[i].cells)),
        _ => panic!("P_ORDER must be stream, m, rows or cells"),
    }
    let (tb, fetch_pre, two) = match kern.as_str() {
        "row" => (true, false, false),
        "tb" => (true, false, false),
        "so" => (false, false, false),
        "tb2" => (true, false, true),
        "tbpre" => (true, true, false),
        "sopre" => (false, true, false),
        "tb2pre" => (true, true, true),
        k => panic!("P_KERN must be row, tb, so, tb2, tbpre, sopre or tb2pre, not {k}"),
    };
    let mut pb = PB::new();
    let mut res: Vec<Out> = vec![Out::default(); list.len()];
    let (mut c1, mut c2) = (Count::default(), Count::default());
    let is_row = kern == "row";
    // check pass (records the scores for the `pre` kernels)
    if !is_row {
        // SAFETY: AVX2 checked above.
        unsafe {
            match (tb, fetch_pre) {
                (true, false) => run_all::<true, 0>(&list, &order, two, &mut pb, &mut c1, &mut res, &mut c2),
                (false, false) => run_all::<false, 0>(&list, &order, two, &mut pb, &mut c1, &mut res, &mut c2),
                (true, true) => run_all::<true, 2>(&list, &order, two, &mut pb, &mut c1, &mut res, &mut c2),
                (false, true) => run_all::<false, 2>(&list, &order, two, &mut pb, &mut c1, &mut res, &mut c2),
            }
        }
        let mut bad = 0usize;
        for i in 0..list.len() {
            let ok = if tb {
                res[i] == want[i]
            } else {
                to_res(&res[i]) == to_res(&want[i])
            };
            if !ok {
                if bad < 5 {
                    eprintln!(
                        "MISMATCH problem {i} len1 {} len2 {}: got {:?} ops {:?} want {:?} ops {:?}",
                        list[i].len1,
                        list[i].len2,
                        to_res(&res[i]),
                        &res[i].ops[..res[i].ops.len().min(6)],
                        to_res(&want[i]),
                        &want[i].ops[..want[i].ops.len().min(6)]
                    );
                }
                bad += 1;
            }
        }
        assert_eq!(bad, 0, "{bad} problems differ from the row kernel");
    }
    let pre_len = pb.pre.len();
    // timed loop
    let t0 = Instant::now();
    let (mut d1, mut d2) = (Count::default(), Count::default());
    for _ in 0..reps {
        if is_row {
            for &i in &order {
                let p = list[i];
                // SAFETY: AVX2 checked above.
                std::hint::black_box(unsafe {
                    align(
                        &p.row_seq(), &p.col_seq(), &p.scores(), p.len1, p.len2, p.gap_open,
                        p.gap_extend, p.x_drop, true, p.check_fence, &mut sc,
                    )
                });
            }
        } else {
            pb.pre_k = 0;
            // SAFETY: AVX2 checked above.
            unsafe {
                match (tb, fetch_pre) {
                    (true, false) => run_all::<true, 0>(&list, &order, two, &mut pb, &mut d1, &mut res, &mut d2),
                    (false, false) => run_all::<false, 0>(&list, &order, two, &mut pb, &mut d1, &mut res, &mut d2),
                    (true, true) => run_all::<true, 1>(&list, &order, two, &mut pb, &mut d1, &mut res, &mut d2),
                    (false, true) => run_all::<false, 1>(&list, &order, two, &mut pb, &mut d1, &mut res, &mut d2),
                }
            }
            std::hint::black_box(&res);
        }
    }
    let loop_s = t0.elapsed().as_secs_f64();
    let steps = c1.steps + c2.steps;
    let u = cells_total as f64 / (16.0 * c1.steps.max(1) as f64);
    let r = reps.max(1) as f64;
    let line = format!(
        "PAIR kern={kern} order={ord} dir={dir} problems={} skipped={skipped} reps={reps} load_s={load_s:.3} \
         loop_s={loop_s:.6} cells={cells_total} rows_rowkernel={} batches={} steps1={} steps2={} \
         rows1={} rows2={} kill_rows={} rebases={} clears={} pre={pre_len} u1={u:.4} \
         ns_per_step={:.3} ns_per_cell={:.4}\n",
        list.len(),
        rows_run.iter().sum::<usize>(),
        c1.batches + c2.batches,
        c1.steps,
        c2.steps,
        c1.rows,
        c2.rows,
        c1.kill_rows + c2.kill_rows,
        c1.rebases + c2.rebases,
        c1.clears + c2.clears,
        if is_row || steps == 0 { 0.0 } else { loop_s * 1e9 / (r * steps as f64) },
        loop_s * 1e9 / (r * cells_total.max(1) as f64),
    );
    eprint!("{line}");
    if let Some(out) = var("P_OUT") {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(out)
            .expect("open P_OUT");
        f.write_all(line.as_bytes()).expect("write P_OUT");
    }
}
