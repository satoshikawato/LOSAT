//! The two-register kernels against the scalar transliteration of the NCBI loops and against
//! the 128-bit NEON kernels of `xdrop_simd.rs` (random problems: forward and reversed rows and
//! columns, both lane widths, traceback and score-only, fences, scale 1 and 32, bands wider
//! than 256 cells, gap costs at both row-gap preconditions; each problem with rows of 128 / 320
//! lanes to start with that move with the band, or held, the NEON2 result equal to a scratch
//! that never moves), and the pair helpers against scalar models of their lanes.
use super::super::kernel;
use super::super::tests::{random_problem, reference, Problem, Rng};
use super::super::FENCE_SENTRY;
use super::*;

unsafe fn v(b: [u8; 16]) -> V {
    vld1q_u8(b.as_ptr())
}
unsafe fn b(x: V) -> [u8; 16] {
    let mut o = [0u8; 16];
    vst1q_u8(o.as_mut_ptr(), x);
    o
}
unsafe fn wv(x: [u8; 32]) -> W {
    ldw(x.as_ptr())
}
unsafe fn wb(x: W) -> [u8; 32] {
    let mut o = [0u8; 32];
    stw(o.as_mut_ptr(), x);
    o
}
fn u16s(x: [u8; 32]) -> [u16; 16] {
    core::array::from_fn(|i| u16::from_le_bytes([x[2 * i], x[2 * i + 1]]))
}
fn bytes(rng: &mut Rng) -> [u8; 32] {
    // mostly small values, some near the ends of the lanes
    core::array::from_fn(|_| match rng.below(4) {
        0 => 255 - rng.below(3) as u8,
        1 => rng.below(3) as u8,
        _ => rng.below(256) as u8,
    })
}

/// Each pair helper against the scalar definition of its lanes, on random inputs.
#[test]
fn neon2_pair_helpers_match_scalar_models() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    for _ in 0..20_000 {
        let x = bytes(&mut rng);
        let x16 = u16s(x);
        // SAFETY: NEON is a baseline feature of the aarch64 targets.
        unsafe {
            // one lane up
            let up = u16s(wb(up1_16(wv(x))));
            for i in 0..16 {
                assert_eq!(up[i], if i == 0 { 0 } else { x16[i - 1] });
            }
            let up = wb(up1_8(wv(x)));
            for i in 0..32 {
                assert_eq!(up[i], if i == 0 { 0 } else { x[i - 1] });
            }
            // prefix maxima
            // the signed forms are applied to `H`, which is never negative: inputs 0..=i16::MAX
            // (zeros are shifted in, so negative lanes would give a mix of 0 and the inputs)
            let xp: [u8; 32] =
                core::array::from_fn(|i| if i % 2 == 1 { x[i] & 0x7F } else { x[i] });
            let xp16 = u16s(xp);
            let pu = u16s(wb(pmax_u16w(wv(x))));
            let ps = u16s(wb(pmax_s16w(wv(xp))));
            let (mut mu, mut ms) = (0u16, 0u16);
            for i in 0..16 {
                mu = mu.max(x16[i]);
                ms = ms.max(xp16[i]);
                assert_eq!(pu[i], mu);
                assert_eq!(ps[i], ms);
            }
            let x7: [u8; 32] = core::array::from_fn(|i| x[i] & 0x7F);
            let pu = wb(pmax_u8w(wv(x)));
            let ps = wb(pmax_s8w(wv(x7)));
            let (mut mu, mut ms) = (0u8, 0u8);
            for i in 0..32 {
                mu = mu.max(x[i]);
                ms = ms.max(x7[i]);
                assert_eq!(pu[i], mu);
                assert_eq!(ps[i], ms);
            }
            // last lane to all lanes
            let l = b(last16(wv(x)));
            for i in 0..8 {
                assert_eq!(u16::from_le_bytes([l[2 * i], l[2 * i + 1]]), x16[15]);
            }
            assert_eq!(b(last8(wv(x))), [x[31]; 16]);
            // any bit
            let sparse: [u8; 32] =
                core::array::from_fn(|i| if rng.below(40) == 0 { x[i] } else { 0 });
            assert_eq!(anyw(wv(sparse)), sparse.iter().any(|&c| c != 0));
            // lane masks to 4 bits per lane
            let m8: [u8; 32] = core::array::from_fn(|i| if x[i] & 1 != 0 { 255 } else { 0 });
            let want8: u128 = (0..32).fold(0, |a, i| {
                a | if m8[i] != 0 { 0xFu128 << (4 * i) } else { 0 }
            });
            assert_eq!(nib8w(wv(m8)), want8);
            let lo: [u8; 16] = core::array::from_fn(|i| m8[i]);
            assert_eq!(nib8(v(lo)), want8 as u64);
            let m16: [u8; 32] = core::array::from_fn(|i| m8[i & !1]);
            let want16: u64 = (0..16).fold(0, |a, i| {
                a | if m16[2 * i] != 0 {
                    0xFu64 << (4 * i)
                } else {
                    0
                }
            });
            assert_eq!(nib16w(wv(m16)), want16);
            // script bytes
            let r: [u8; 32] =
                core::array::from_fn(|i| if i % 2 == 0 { (x[i] % 0x57) as u8 } else { 0 });
            let p = b(pack16(wv(r)));
            for i in 0..16 {
                assert_eq!(p[i], r[2 * i]);
            }
            // 16-bit scores from the two byte planes (an index of 32 and up gives 0)
            let tl = bytes(&mut rng);
            let th = bytes(&mut rng);
            let idx: [u8; 16] =
                core::array::from_fn(|i| if rng.below(8) == 0 { x[i] } else { x[i] % 32 });
            let s = u16s(wb(scores16(
                uint8x16x2_t(ldw(tl.as_ptr()).0, ldw(tl.as_ptr()).1),
                uint8x16x2_t(ldw(th.as_ptr()).0, ldw(th.as_ptr()).1),
                v(idx),
            )));
            for i in 0..16 {
                let k = idx[i] as usize;
                let want = if k < 32 {
                    u16::from_le_bytes([tl[k], th[k]])
                } else {
                    0
                };
                assert_eq!(s[i], want, "scores16 lane {i} index {k}");
            }
        }
    }
}

/// The 128-bit NEON kernel of `xdrop_simd.rs` or the kernels of this file; `m8` selects the
/// 8-bit lanes (as `fits_8bit` would).
unsafe fn run(
    two: bool,
    m8: Option<i32>,
    tb: bool,
    q: &RowSeq<'_>,
    s: &ColSeq<'_>,
    scores: &Scores<'_>,
    len1: usize,
    len2: usize,
    go: i32,
    ge: i32,
    x: i32,
    cf: bool,
    sc: &mut XdropScratch,
) -> Option<XdropResult> {
    match (two, tb, m8) {
        (false, true, Some(m)) => {
            kernel::align_body8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (false, false, Some(m)) => {
            kernel::align_body8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (false, true, None) => {
            kernel::align_body::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
        (false, false, None) => {
            kernel::align_body::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
        (true, true, Some(m)) if sc.x_nib => {
            align8::<true, true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (true, true, Some(m)) => {
            align8::<true, false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (true, false, Some(m)) => {
            align8::<false, false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (true, true, None) if sc.x_nib => {
            align16::<true, true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
        (true, true, None) => align16::<true, false>(q, s, scores, len1, len2, go, ge, x, cf, sc),
        (true, false, None) => align16::<false, false>(q, s, scores, len1, len2, go, ge, x, cf, sc),
    }
}

/// Repetitive sequences: every diagonal scores alike, so each row keeps the widest band
/// x_drop allows.
fn repeat_problem(rng: &mut Rng, n: usize, scale: i32, len: usize) -> Problem {
    let mut p = random_problem(rng, n, scale, 1);
    let letters = (n - 1) as u64;
    let period = 1 + rng.below(3) as usize;
    let motif: Vec<u8> = (0..period).map(|_| 1 + rng.below(letters) as u8).collect();
    p.q = (0..len).map(|i| motif[i % period]).collect();
    p.s = (0..len + rng.below(50) as usize)
        .map(|i| motif[i % period])
        .collect();
    for _ in 0..rng.below(len as u64 / 20 + 1) {
        let k = rng.below(p.s.len() as u64) as usize;
        p.s[k] = 1 + rng.below(letters) as u8;
    }
    p
}

fn junk(rng: &mut Rng, k: usize) -> Vec<u8> {
    (0..k)
        .map(|_| match rng.below(4) {
            0 => FENCE_SENTRY,
            1 => 255,
            _ => rng.below(32) as u8,
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:541-557
// Random problems run through the scalar loops (`reference`), the 128-bit NEON kernels and the
// kernels of this file, with the rows and columns given forward and reversed (as the call sites
// pass them for the left extension), on 16-bit lanes and, where the problem fits, 8-bit lanes,
// with and without traceback. Offsets, score, fence flag, edit script and cell count must be
// equal; a call these kernels decline must be one that fails a row-gap precondition. The
// rows start at 128 / 320 lanes and move with the band, or are held, and the result is equal
// either way.
// `LOSAT_FUZZ_CASES` sets the number of problems.
#[test]
fn neon2_kernels_match_scalar_loops_on_random_problems() {
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(600);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let (mut sc128, mut sc2) = (XdropScratch::new(), XdropScratch::new());
    // the same kernels with every column held, for the comparison of the rows over the band
    let mut sc2_off = XdropScratch::new();
    sc2_off.x_test_rows(None);
    // calls whose rows were moved: [width 8 / 16][so / tb][128-bit / neon2]
    let mut moved = [[[0usize; 2]; 2]; 2];
    // [width 8 / 16][tb][reverse]
    let mut handled = [[[0usize; 2]; 2]; 2];
    let (mut fences, mut ops_total, mut wide) = (0usize, 0usize, 0usize);
    // calls at the row-gap preconditions: 8-bit with 32 * gap_extend - gap_open_extend in
    // 113..=128 (handled), above 128 or 16-bit above 32768 (declined, run by the 128-bit
    // kernel); 16-bit with gap_extend >= 1000 handled
    let (mut edge_in, mut edge_out, mut big_ge) = (0usize, 0usize, 0usize);
    for case in 0..cases {
        // rows over the band (`x_room`): 128 or 320 lanes to start with, or every column
        let lanes = [Some(128), Some(320), None][case % 3];
        sc128.x_test_rows(lanes);
        sc2.x_test_rows(lanes);
        // trace layout (LOSAT_X_DPNIBBLE): four-bit codes or one byte per cell
        let nib = case % 2 == 0;
        (sc128.x_nib, sc2.x_nib) = (nib, nib);
        let n = if rng.below(2) == 0 { 16 } else { 28 };
        let mut scale = if rng.below(3) == 0 { 32 } else { 1 };
        let long = rng.below(6) == 0;
        let len = 1 + rng.below(if long { 2500 } else { 300 }) as usize;
        let kind = rng.below(10);
        if kind == 8 {
            scale = 1;
        }
        let repeats = rng.below(4) == 0;
        let p = if repeats {
            repeat_problem(&mut rng, n, scale, len)
        } else {
            random_problem(&mut rng, n, scale, len)
        };
        let (go, ge, x) = match kind {
            0 => (11 * scale, scale, 38 * scale),
            1 => (11 * scale, scale, 64 * scale + rng.below(5) as i32),
            2 => (5 * scale, 2 * scale, 33 * scale),
            3 => (5 * scale, 2 * scale, 110 * scale),
            4 => (
                rng.below(12) as i32 * scale,
                (1 + rng.below(3) as i32) * scale,
                (3 + rng.below(150) as i32) * scale,
            ),
            5 => (9 * scale, scale, 3 * scale),
            // wide bands (more than 256 cells on 16-bit lanes)
            6 => (
                rng.below(12) as i32 * scale,
                scale,
                (150 + rng.below(450) as i32) * scale,
            ),
            // small costs for the 8-bit lanes
            7 => (
                rng.below(8) as i32,
                1 + rng.below(5) as i32,
                10 + rng.below(80) as i32,
            ),
            // 8-bit lanes at the row-gap precondition (32 * ge - goe around 128)
            8 => {
                let ge = 4 + rng.below(2) as i32;
                let go = rng.below(45) as i32;
                (go, ge, go + ge + rng.below(30) as i32)
            }
            // 16-bit lanes with gap_extend up to the 128-bit kernel's limit (8 * ge <= 32767)
            _ => {
                let ge = 1000 + rng.below(3096) as i32;
                let go = rng.below(2000) as i32;
                (go, ge, go + ge + rng.below(4000) as i32)
            }
        };
        let len1 = p.q.len();
        if p.s.is_empty() {
            continue;
        }
        for tb in [false, true] {
            let mut s = p.s.clone();
            let mut len2 = s.len();
            let mut sentinel = false;
            if tb {
                match rng.below(6) {
                    0 => {
                        s.push(FENCE_SENTRY);
                        sentinel = true;
                    }
                    1 if len2 > 3 => {
                        let k = 1 + rng.below(len2 as u64 - 1) as usize;
                        s[k] = FENCE_SENTRY;
                    }
                    2 if len2 > 3 => {
                        len2 = 1 + rng.below(len2 as u64 - 1) as usize;
                        sentinel = rng.below(2) == 0;
                    }
                    _ => {}
                }
            }
            let prob = Problem {
                q: p.q.clone(),
                s: s.clone(),
                n,
                m: p.m.clone(),
            };
            let expected = reference(&prob, len1, len2, go, ge, x, tb, sentinel);
            let scores = Scores::Flat { data: &p.m, n };
            // forward forms, at an offset inside a larger buffer
            let pad_q = rng.below(5) as usize;
            let pad_s = rng.below(5) as usize;
            let mut fq = junk(&mut rng, pad_q);
            fq.extend_from_slice(&p.q);
            let mut fs = junk(&mut rng, pad_s);
            fs.extend_from_slice(&s);
            // reversed forms (`A[M - a_index]`, `B[N - b]`): the columns stop at len2
            let mut rq = junk(&mut rng, pad_q);
            rq.extend(p.q.iter().rev());
            let mut rs = junk(&mut rng, pad_s);
            rs.extend(s[..len2].iter().rev());
            let tail = junk(&mut rng, 3);
            rs.extend_from_slice(&tail);
            for reverse in [false, true] {
                if reverse && sentinel {
                    continue;
                }
                let (rows, cols) = if reverse {
                    (
                        RowSeq::Bytes {
                            data: &rq,
                            base: (pad_q + len1) as isize,
                            step: -1,
                        },
                        ColSeq {
                            data: &rs,
                            base: (pad_s + len2) as isize - 1,
                            step: -1,
                            zero_from: len2,
                        },
                    )
                } else {
                    (
                        RowSeq::Bytes {
                            data: &fq,
                            base: pad_q as isize - 1,
                            step: 1,
                        },
                        ColSeq {
                            data: &fs,
                            base: pad_s as isize,
                            step: 1,
                            zero_from: if sentinel { usize::MAX } else { len2 },
                        },
                    )
                };
                let m8 = fits_8bit(&scores, go, ge, x);
                for wide8 in [false, true] {
                    let m = if wide8 {
                        match m8 {
                            Some(m) => Some(m),
                            None => continue,
                        }
                    } else {
                        None
                    };
                    let tag = format!(
                        "case {case} tb={tb} rev={reverse} w8={wide8} n={n} scale={scale} \
                             go={go} ge={ge} x={x} len1={len1} len2={len2} rep={repeats}"
                    );
                    // SAFETY: NEON is a baseline feature of the aarch64 targets.
                    let (r128, r2, r2_off) = unsafe {
                        (
                            run(
                                false, m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                &mut sc128,
                            ),
                            run(
                                true, m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                &mut sc2,
                            ),
                            run(
                                true,
                                m,
                                tb,
                                &rows,
                                &cols,
                                &scores,
                                len1,
                                len2,
                                go,
                                ge,
                                x,
                                tb,
                                &mut sc2_off,
                            ),
                        )
                    };
                    // rows moved or held give the same offsets, score, fence flag, cells and
                    // edit script (`sc2_off` never moves; `sc2.x_sh` is 0 when its rows hold
                    // every column)
                    assert_eq!(r2, r2_off, "{tag}: neon2 rows {lanes:?} vs every column");
                    assert_eq!(sc2.ops, sc2_off.ops, "{tag}: neon2 moved rows ops");
                    assert_eq!(sc2_off.x_sh, 0, "{tag}: neon2 off moved");
                    let edge = if wide8 {
                        32 * ge - (go + ge)
                    } else {
                        16 * ge - (go + ge)
                    };
                    match (r128, r2) {
                        (Some(a), Some(b)) => {
                            assert_eq!(a, b, "{tag}: 128 vs neon2");
                            assert_eq!(sc128.ops, sc2.ops, "{tag}: 128 vs neon2 ops");
                        }
                        (None, None) => {}
                        (Some(_), None) => {
                            assert!(
                                edge > if wide8 { 128 } else { 32768 },
                                "{tag}: neon2 declined"
                            );
                            edge_out += 1;
                        }
                        (None, Some(_)) => {
                            panic!("{tag}: neon2 handled a call the 128-bit kernel declines")
                        }
                    }
                    let Some(r) = r2 else {
                        continue;
                    };
                    assert_eq!(
                        (r.a_offset, r.b_offset, r.score, r.fence_hit),
                        (expected.0, expected.1, expected.2, expected.3),
                        "{tag}: neon2 vs scalar"
                    );
                    assert_eq!(sc2.ops, expected.4, "{tag}: neon2 vs scalar ops");
                    if !r.fence_hit {
                        assert_eq!(r.cells, expected.5, "{tag}: cells");
                    }
                    handled[usize::from(!wide8)][usize::from(tb)][usize::from(reverse)] += 1;
                    moved[usize::from(!wide8)][usize::from(tb)][0] += usize::from(sc128.x_sh > 0);
                    moved[usize::from(!wide8)][usize::from(tb)][1] += usize::from(sc2.x_sh > 0);
                    fences += usize::from(r.fence_hit);
                    ops_total += sc2.ops.len();
                    if wide8 && (113..=128).contains(&edge) {
                        edge_in += 1;
                    }
                    if !wide8 && ge >= 1000 {
                        big_ge += 1;
                    }
                    if !wide8 && r.cells > 256 * len1.min(len2) as u64 && len1.min(len2) > 500 {
                        wide += 1;
                    }
                }
            }
        }
    }
    eprintln!(
        "neon2 fuzz: cases={cases} handled[w8/w16][so/tb][fwd/rev]={handled:?} fences={fences} \
             wide={wide} edge_in={edge_in} edge_out={edge_out} big_ge={big_ge} \
             moved[w8/w16][so/tb][128/neon2]={moved:?}"
    );
    if cases >= 200 {
        for (w, by_tb) in handled.iter().enumerate() {
            for (t, by_dir) in by_tb.iter().enumerate() {
                for (d, &count) in by_dir.iter().enumerate() {
                    assert!(
                        count > 0,
                        "nothing handled for width {w} tb {t} reverse {d}"
                    );
                }
            }
        }
        assert!(fences > 0 && ops_total > 0 && wide > 0);
        assert!(edge_in > 0 && edge_out > 0 && big_ge > 0);
        // the kernels of this file and the 128-bit kernels move their rows, on both lane
        // widths, with and without the traceback
        assert!(
            moved.iter().flatten().flatten().all(|&n| n > 0),
            "moved rows {moved:?}"
        );
    }
}

/// Two similar sequences of about `len` letters (a few substitutions and indels), so the
/// band follows a diagonal for the whole length and its start passes the rows' first length.
fn similar_problem(rng: &mut Rng, n: usize, len: usize) -> Problem {
    let mut p = random_problem(rng, n, 1, len);
    let letters = (n - 1) as u64;
    let mut s = Vec::with_capacity(len + 64);
    let mut i = 0;
    while i < p.q.len() {
        let u = rng.unit();
        if u < 0.01 {
            for _ in 0..=rng.below(3) {
                s.push(1 + rng.below(letters) as u8);
            }
        } else if u < 0.02 {
            i += 1 + rng.below(3) as usize;
            continue;
        }
        s.push(if rng.unit() < 0.03 {
            1 + rng.below(letters) as u8
        } else {
            p.q[i]
        });
        i += 1;
    }
    p.s = s;
    p
}

// No NCBI counterpart: the rows over the band of the NEON2 kernels (`XdropScratch::x_room`).
// Long similar problems (gap costs (5, 2, 110) on 16-bit lanes and (11, 1, 38) on 8-bit and
// 16-bit lanes, forward and reversed columns, score-only and traceback) run `align16` and
// `align8` with rows of 128 lanes to start with that move with the band and with every column
// held: offsets, score, fence flag, cell count and edit script are equal to each
// other, to the 128-bit kernel and to the scalar loops, and the moved rows stay short (they
// would hold every column up to the band's end, 1500 or more lanes, if they were not moved).
// Every call is inside the preconditions of the kernels, so none may be declined.
#[test]
fn neon2_band_window_keeps_every_result() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let (mut on, mut off, mut sc128) = (
        XdropScratch::new(),
        XdropScratch::new(),
        XdropScratch::new(),
    );
    // [width 8 / 16][so / tb]: calls handled and calls whose rows were moved
    let mut handled = [[0usize; 2]; 2];
    let mut shifted = [[0usize; 2]; 2];
    for case in 0..48 {
        let lanes = Some(128);
        on.x_test_rows(lanes);
        off.x_test_rows(None);
        sc128.x_test_rows(lanes);
        let n = if case % 2 == 0 { 16 } else { 28 };
        let len = 1500 + rng.below(2500) as usize;
        let p = similar_problem(&mut rng, n, len);
        // (5, 2, 110) needs 16-bit lanes, (11, 1, 38) fits 8-bit lanes
        let (go, ge, x) = if case % 4 < 2 {
            (5, 2, 110)
        } else {
            (11, 1, 38)
        };
        let (len1, len2) = (p.q.len(), p.s.len());
        if len2 == 0 {
            continue;
        }
        let scores = Scores::Flat { data: &p.m, n };
        let m8 = fits_8bit(&scores, go, ge, x);
        let rq: Vec<u8> = p.q.iter().rev().copied().collect();
        let rs: Vec<u8> = p.s.iter().rev().copied().collect();
        let expected = [false, true].map(|tb| reference(&p, len1, len2, go, ge, x, tb, false));
        for reverse in [false, true] {
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
                        data: &p.s,
                        base: 0,
                        step: 1,
                        zero_from: len2,
                    },
                )
            };
            for tb in [false, true] {
                // the reversed forms present the same logical sequences
                let expected = &expected[usize::from(tb)];
                // 8-bit lanes where the costs fit them, and always 16-bit lanes
                let widths: &[Option<i32>] = if m8.is_some() { &[m8, None] } else { &[None] };
                for &m in widths {
                    let tag = format!(
                        "case {case} tb={tb} rev={reverse} w8={} go={go} ge={ge} x={x} \
                             len1={len1} len2={len2}",
                        m.is_some()
                    );
                    // SAFETY: NEON is a baseline feature of the aarch64 targets.
                    let (a, b, c) = unsafe {
                        (
                            run(
                                true, m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                &mut on,
                            ),
                            run(
                                true, m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                &mut off,
                            ),
                            run(
                                false, m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb,
                                &mut sc128,
                            ),
                        )
                    };
                    let (Some(a), Some(b), Some(c)) = (a, b, c) else {
                        panic!(
                            "{tag}: declined (moved {} / every column {} / 128-bit {})",
                            a.is_some(),
                            b.is_some(),
                            c.is_some()
                        );
                    };
                    let w = usize::from(m.is_none());
                    handled[w][usize::from(tb)] += 1;
                    assert_eq!(a, b, "{tag}: moved rows vs every column");
                    assert_eq!(a, c, "{tag}: neon2 moved rows vs the 128-bit kernel");
                    assert_eq!(on.ops, off.ops, "{tag}: traceback, moved vs every column");
                    assert_eq!(on.ops, sc128.ops, "{tag}: traceback, 128-bit kernel");
                    assert_eq!(off.x_sh, 0, "{tag}: held rows moved");
                    assert_eq!(
                        (a.a_offset, a.b_offset, a.score, a.fence_hit),
                        (expected.0, expected.1, expected.2, expected.3),
                        "{tag}: against the scalar loops"
                    );
                    assert_eq!(on.ops, expected.4, "{tag}: scalar traceback");
                    if on.x_sh > 0 {
                        shifted[w][usize::from(tb)] += 1;
                        // the rows start at 128 lanes and grow only where the band needs it
                        // (`x_room`): within a few times the band, its extra lanes and the
                        // padding, not the 1500 or more columns they would hold unmoved
                        let rows_len = on.ha.len().max(on.h8a.len());
                        assert!(
                            rows_len <= 1024,
                            "{tag}: rows of {rows_len} lanes after moving {} columns",
                            on.x_sh
                        );
                    }
                }
            }
        }
    }
    eprintln!("neon2 moved rows: handled[w8/w16][so/tb]={handled:?} moved={shifted:?}");
    // 48 cases, both column directions each: 24 with (5, 2, 110) on 16-bit lanes, 24 with
    // (11, 1, 38) on 8-bit and on 16-bit lanes
    assert_eq!(
        handled,
        [[48, 48], [96, 96]],
        "handled {handled:?}, moved {shifted:?}"
    );
    assert!(
        shifted.iter().flatten().all(|&n| n > 10),
        "handled {handled:?}, moved {shifted:?}"
    );
}
