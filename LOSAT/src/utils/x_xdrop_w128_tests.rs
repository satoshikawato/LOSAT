//! These kernels against the scalar transliteration of the NCBI loops and against the 128-bit
//! kernels of `xdrop_simd.rs` (random problems: forward and reversed rows and columns, both
//! lane widths, traceback and score-only, fences, scale 1 and 32, bands wider than 256 cells,
//! gap costs at the row-gap precondition), the vector row tables and the column copy. They run
//! natively on the x86_64 build of the same code. The replay of captured problems
//! (`LOSAT_X_DPREPLAY`) is in `x_xdrop_avx2.rs`.
use super::super::kernel;
use super::super::tests::{random_problem, reference, Problem, Rng};
use super::super::x_capture::load_dir;
use super::*;

/// The 128-bit kernel of `xdrop_simd.rs` in the build its dispatch picks; `m8` selects the
/// 8-bit lanes (as `fits_8bit` would).
///
/// # Safety
/// `available()` must be true.
unsafe fn xs128(
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
    #[cfg(target_arch = "x86_64")]
    {
        let avx = super::super::cpu_level() == 2;
        match (avx, tb, m8) {
            (true, true, Some(m)) => {
                kernel::align8_avx::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (true, false, Some(m)) => {
                kernel::align8_avx::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (false, true, Some(m)) => {
                kernel::align8_sse::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (false, false, Some(m)) => {
                kernel::align8_sse::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (true, true, None) => {
                kernel::align_avx::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (true, false, None) => {
                kernel::align_avx::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (false, true, None) => {
                kernel::align_sse::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
            (false, false, None) => {
                kernel::align_sse::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        match (tb, m8) {
            (true, Some(m)) => {
                kernel::align_body8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (false, Some(m)) => {
                kernel::align_body8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
            (true, None) => kernel::align_body::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc),
            (false, None) => {
                kernel::align_body::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
        }
    }
}

/// The kernels of this file; `m8` selects the 8-bit lanes.
///
/// # Safety
/// `available()` must be true.
unsafe fn w128(
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
    match (tb, m8) {
        (true, Some(m)) => run8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc),
        (false, Some(m)) => run8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc),
        (true, None) => run16::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc),
        (false, None) => run16::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc),
    }
}

/// The vector row tables are byte for byte those of `build_row_table`, and are rejected in
/// the same cases.
#[test]
fn w128_row_tables_match_build_row_table() {
    if !available() {
        eprintln!("no SSE4.1: skipped");
        return;
    }
    let special = [
        i32::MIN,
        i32::MIN / 2,
        -40000,
        -32769,
        -32768,
        -32767,
        -129,
        -128,
        -127,
        -1,
        0,
        1,
        127,
        128,
        1999,
        2000,
        2001,
        32767,
        i32::MAX,
    ];
    let mut rng = Rng(0x5151_7A7A_0101_2324);
    let mut accepted = 0usize;
    for _ in 0..3000 {
        let rare = rng.below(8) == 0;
        let pick = |rng: &mut Rng| {
            if rare && rng.below(50) == 0 {
                special[rng.below(special.len() as u64) as usize]
            } else if rng.below(10) == 0 {
                special[rng.below(14) as usize]
            } else {
                rng.below(300) as i32 - 150
            }
        };
        let mut m = [[0i32; 28]; 28];
        for row in m.iter_mut() {
            for v in row.iter_mut() {
                *v = pick(&mut rng);
            }
        }
        let q = rng.below(30) as usize;
        let (mut a, mut b) = ([0u8; TAB], [0u8; TAB]);
        let rows = Scores::Rows28(&m);
        let ra = build_row_table(&rows, q, &mut a);
        // SAFETY: SSE4.1 checked above (or simd128).
        let rb = unsafe { build_row(&rows, q, &mut b) };
        assert_eq!(ra, rb, "Rows28 q={q}");
        if ra {
            assert_eq!(a, b, "Rows28 q={q}");
            accepted += 1;
        }
        let n = [1usize, 5, 8, 9, 16, 20, 28, 31, 32, 33][rng.below(10) as usize];
        let data: Vec<i32> = (0..n * n).map(|_| pick(&mut rng)).collect();
        let q = rng.below(n as u64 + 2) as usize;
        let flat = Scores::Flat { data: &data, n };
        let (mut a, mut b) = ([0u8; TAB], [0u8; TAB]);
        let ra = build_row_table(&flat, q, &mut a);
        // SAFETY: as above.
        let rb = unsafe { build_row(&flat, q, &mut b) };
        assert_eq!(ra, rb, "Flat n={n} q={q}");
        if ra {
            assert_eq!(a, b, "Flat n={n} q={q}");
            accepted += 1;
        }
    }
    assert!(accepted > 1000, "only {accepted} tables accepted");
}

/// The column copy fills the lane buffer, the first fence column and the first bad column as
/// the scalar loop of the 128-bit kernels does, forward and backward, near both ends of the
/// data and with `zero_from` inside the range.
#[test]
fn w128_fill_cols_matches_scalar_loop() {
    if !available() {
        eprintln!("no SSE4.1: skipped");
        return;
    }
    let mut rng = Rng(0x0BAD_C0DE_1234_5678);
    let mut vector_blocks = 0usize;
    for case in 0..4000 {
        let len = rng.below(200) as usize;
        let ncols = [4usize, 16, 25, 28, 32][rng.below(5) as usize];
        let data: Vec<u8> = (0..len)
            .map(|_| match rng.below(60) {
                0 => FENCE_SENTRY,
                1 => ncols as u8 + rng.below(20) as u8,
                2 => 255,
                _ => rng.below(ncols as u64) as u8,
            })
            .collect();
        let step: isize = if rng.below(2) == 0 { 1 } else { -1 };
        let base = rng.below(len as u64 + 10) as isize - 5;
        let zero_from = if rng.below(3) == 0 {
            rng.below(len as u64 + 20) as usize
        } else {
            usize::MAX
        };
        let s = ColSeq {
            data: &data,
            base,
            step,
            zero_from,
        };
        let check_fence = rng.below(2) == 0;
        let start = rng.below(40) as usize;
        let to = start + rng.below(150) as usize;
        let mut a = XdropScratch::new();
        let mut b = XdropScratch::new();
        a.ensure(to + 64);
        b.ensure(to + 64);
        for sc in [&mut a, &mut b] {
            sc.sb_filled = start;
            sc.fence_k = usize::MAX;
            sc.bad_k = usize::MAX;
        }
        // the scalar loop of `align_body` / `align_body8`
        for k in a.sb_filled..to {
            let ch = s.get(k);
            if ch as usize >= ncols {
                if check_fence && ch == FENCE_SENTRY {
                    if a.fence_k == usize::MAX {
                        a.fence_k = k;
                    }
                } else if a.bad_k == usize::MAX {
                    a.bad_k = k;
                }
            }
            a.sb[PAD + k + 1] = ch;
        }
        a.sb_filled = to;
        // SAFETY: SSE4.1 checked above (or simd128); the buffer holds `to + 64` lanes.
        unsafe { fill_cols(&s, &mut b, to, ncols, check_fence) };
        assert_eq!(
            (a.sb_filled, a.fence_k, a.bad_k),
            (b.sb_filled, b.fence_k, b.bad_k),
            "case {case}"
        );
        assert_eq!(a.sb, b.sb, "case {case}");
        if (start..to).filter(|k| s.get(*k) < ncols as u8).count() >= 16 {
            vector_blocks += 1;
        }
    }
    assert!(vector_blocks > 500, "{vector_blocks}");
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
// Random problems run through the scalar loops (`reference`), the 128-bit kernels of
// `xdrop_simd.rs` and the kernels of this file, with the rows and columns given forward and
// reversed (as the call sites pass them for the left extension), on 16-bit lanes and, where
// the problem fits, 8-bit lanes, with and without traceback. Offsets, score, fence flag, edit
// script and cell count must be equal; a call these kernels decline must be one that fails the
// row-gap precondition. `LOSAT_FUZZ_CASES` sets the number of problems.
#[test]
fn w128_kernels_match_scalar_loops_on_random_problems() {
    if !available() {
        eprintln!("no SSE4.1: skipped");
        return;
    }
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(600);
    let mut rng = Rng(0x7F4A_7C15_D1B5_4A33);
    let (mut sc128, mut scw) = (XdropScratch::new(), XdropScratch::new());
    // [width 8 / 16][tb][reverse]
    let mut handled = [[[0usize; 2]; 2]; 2];
    let (mut fences, mut ops_total, mut wide) = (0usize, 0usize, 0usize);
    // 8-bit calls with 16 * gap_extend - gap_open_extend in 113..=128 (handled) and above 128
    // (declined, run by the 128-bit kernel)
    let (mut edge_in, mut edge_out, mut big_ge) = (0usize, 0usize, 0usize);
    // calls whose rows were moved: [width 8 / 16][128-bit / w128]
    let mut moved = [[0usize; 2]; 2];
    for case in 0..cases {
        // rows over the band (`x_room`): 128 or 320 lanes to start with, or every column
        let lanes = [Some(128), Some(320), None][case % 3];
        sc128.x_test_rows(lanes);
        scw.x_test_rows(lanes);
        // trace layout (LOSAT_X_DPNIBBLE): four-bit codes or one byte per cell
        let nib = case % 2 == 0;
        (sc128.x_nib, scw.x_nib) = (nib, nib);
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
            // 8-bit lanes at the row-gap precondition (16 * ge - goe around 128)
            8 => {
                let ge = 8 + rng.below(4) as i32;
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
                    // SAFETY: SSE4.1 checked above (or simd128).
                    let (r128, rw) = unsafe {
                        (
                            xs128(
                                m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb, &mut sc128,
                            ),
                            w128(
                                m, tb, &rows, &cols, &scores, len1, len2, go, ge, x, tb, &mut scw,
                            ),
                        )
                    };
                    let edge = 16 * ge - (go + ge);
                    match (r128, rw) {
                        (Some(a), Some(b)) => {
                            assert_eq!(a, b, "{tag}: 128 vs w128");
                            assert_eq!(sc128.ops, scw.ops, "{tag}: 128 vs w128 ops");
                        }
                        (None, None) => {}
                        (Some(_), None) => {
                            assert!(wide8 && edge > 128, "{tag}: w128 declined");
                            edge_out += 1;
                        }
                        (None, Some(_)) => {
                            panic!("{tag}: w128 handled a call the 128-bit kernel declines")
                        }
                    }
                    let Some(r) = rw else {
                        continue;
                    };
                    assert_eq!(
                        (r.a_offset, r.b_offset, r.score, r.fence_hit),
                        (expected.0, expected.1, expected.2, expected.3),
                        "{tag}: w128 vs scalar"
                    );
                    assert_eq!(scw.ops, expected.4, "{tag}: w128 vs scalar ops");
                    if !r.fence_hit {
                        assert_eq!(r.cells, expected.5, "{tag}: cells");
                    }
                    handled[usize::from(!wide8)][usize::from(tb)][usize::from(reverse)] += 1;
                    moved[usize::from(!wide8)][0] += usize::from(sc128.x_sh > 0);
                    moved[usize::from(!wide8)][1] += usize::from(scw.x_sh > 0);
                    fences += usize::from(r.fence_hit);
                    ops_total += scw.ops.len();
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
        "w128 fuzz: cases={cases} handled[w8/w16][so/tb][fwd/rev]={handled:?} fences={fences} \
             wide={wide} edge_in={edge_in} edge_out={edge_out} big_ge={big_ge} \
             moved[w8/w16][128/w128]={moved:?}"
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
        assert!(
            moved.iter().flatten().all(|&n| n > 0),
            "moved rows {moved:?}"
        );
    }
}

/// Replay of captured problems (`LOSAT_X_DPCAPTURE` output) through the 128-bit kernel of
/// `xdrop_simd.rs` and the kernels of this file, both compared with the recorded result, then
/// both timed over all problems per lane width and traceback / score-only: "cold" = the
/// problems back to back in capture order (median of `LOSAT_X_DPREPLAY_REPS` repetitions,
/// default 7), "warm-min" = each call timed alone after its residues and matrix were read,
/// fastest of the repetitions, summed. It runs on every target with the 128-bit kernels, the
/// wasm32 simd128 and aarch64 builds included (the replay in `x_xdrop_avx2.rs` needs AVX2).
/// `LOSAT_X_DPREPLAY_OUT=<file>` appends the table to a file.
#[test]
#[ignore = "needs LOSAT_X_DPREPLAY=<capture directory>"]
fn w128_replay_captured_problems() {
    use std::time::Instant;
    let Some(dir) = std::env::var_os("LOSAT_X_DPREPLAY") else {
        eprintln!("LOSAT_X_DPREPLAY is not set: nothing to replay");
        return;
    };
    assert!(available(), "the replay needs the 128-bit kernels");
    let dir = std::path::PathBuf::from(dir);
    let probs = load_dir(&dir).expect("read the capture files");
    let reps: usize = std::env::var("LOSAT_X_DPREPLAY_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(7)
        .max(1);
    let (mut sa, mut sw) = (XdropScratch::new(), XdropScratch::new());
    // categories: [tb][8-bit] -> (problem indices, cells)
    let mut cats: [[(Vec<usize>, u64); 2]; 2] = Default::default();
    let (mut unhandled, mut declined) = (0usize, 0usize);
    let mut m8s = Vec::with_capacity(probs.len());
    for (i, p) in probs.iter().enumerate() {
        let (q, s, scores) = (p.row_seq(), p.col_seq(), p.scores());
        let m8 = fits_8bit(&scores, p.gap_open, p.gap_extend, p.x_drop);
        m8s.push(m8);
        let Some(want) = p.result else {
            unhandled += 1;
            continue;
        };
        let args = (p.len1, p.len2, p.gap_open, p.gap_extend, p.x_drop);
        // SAFETY: `available()` checked above.
        let (r128, rw) = unsafe {
            (
                xs128(
                    m8,
                    p.tb,
                    &q,
                    &s,
                    &scores,
                    args.0,
                    args.1,
                    args.2,
                    args.3,
                    args.4,
                    p.check_fence,
                    &mut sa,
                ),
                w128(
                    m8,
                    p.tb,
                    &q,
                    &s,
                    &scores,
                    args.0,
                    args.1,
                    args.2,
                    args.3,
                    args.4,
                    p.check_fence,
                    &mut sw,
                ),
            )
        };
        assert_eq!(r128, Some(want), "problem {i}: 128-bit vs recorded");
        assert_eq!(sa.ops, p.ops, "problem {i}: 128-bit ops vs recorded");
        match rw {
            Some(r) => {
                assert_eq!(r, want, "problem {i}: w128 vs recorded");
                assert_eq!(sw.ops, p.ops, "problem {i}: w128 ops vs recorded");
            }
            None => declined += 1,
        }
        let cat = &mut cats[usize::from(p.tb)][usize::from(m8.is_some())];
        cat.0.push(i);
        cat.1 += want.cells;
    }
    let checked: usize = cats.iter().flatten().map(|c| c.0.len()).sum();
    let mut report = format!(
        "w128 replay {} ({}): problems={} checked={checked} (128-bit, w128: 0 mismatches) \
             unhandled={unhandled} declined_by_w128={declined} reps={reps}\n",
        dir.display(),
        std::env::consts::ARCH,
        probs.len()
    );
    // One call as `xdrop_align` makes it: a call w128 declines goes to the 128-bit kernel.
    let call = |w: bool, i: usize, sc: &mut XdropScratch| {
        let p = &probs[i];
        let (q, s, scores) = (p.row_seq(), p.col_seq(), p.scores());
        let a = (
            p.len1,
            p.len2,
            p.gap_open,
            p.gap_extend,
            p.x_drop,
            p.check_fence,
        );
        // SAFETY: `available()` checked above.
        unsafe {
            let r = if w {
                w128(
                    m8s[i], p.tb, &q, &s, &scores, a.0, a.1, a.2, a.3, a.4, a.5, sc,
                )
            } else {
                None
            };
            if r.is_none() {
                xs128(
                    m8s[i], p.tb, &q, &s, &scores, a.0, a.1, a.2, a.3, a.4, a.5, sc,
                )
            } else {
                r
            }
        }
    };
    let touch = |i: usize| {
        let p = &probs[i];
        let mut acc = 0i64;
        for &b in p.rows.iter().chain(p.cols.iter()) {
            acc += b as i64;
        }
        for &v in &p.matrix {
            acc += v as i64;
        }
        std::hint::black_box(acc);
    };
    let timer_ns = {
        let mut v: Vec<f64> = (0..15)
            .map(|_| {
                let t0 = Instant::now();
                for _ in 0..10_000 {
                    std::hint::black_box(Instant::now().elapsed());
                }
                t0.elapsed().as_nanos() as f64 / 10_000.0
            })
            .collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    report += &format!("  timer cost {timer_ns:.1} ns/call (subtracted from warm calls)\n");
    for (tb, by_w) in cats.iter().enumerate() {
        for (w8, (list, cells)) in by_w.iter().enumerate() {
            if list.is_empty() {
                continue;
            }
            let cat = format!(
                "{}-{}",
                if tb == 1 { "traceback" } else { "score-only" },
                if w8 == 1 { 8 } else { 16 }
            );
            let mut cold = [Vec::new(), Vec::new()];
            let mut best = [vec![f64::MAX; list.len()], vec![f64::MAX; list.len()]];
            for rep in 0..reps {
                for slot in 0..2 {
                    let w = (slot + rep) % 2 == 1;
                    let sc = if w { &mut sw } else { &mut sa };
                    let t0 = Instant::now();
                    for &i in list {
                        std::hint::black_box(call(w, i, sc));
                    }
                    cold[usize::from(w)].push(t0.elapsed().as_nanos() as f64);
                    let bk = &mut best[usize::from(w)];
                    for (j, &i) in list.iter().enumerate() {
                        touch(i);
                        let t0 = Instant::now();
                        std::hint::black_box(call(w, i, sc));
                        bk[j] = bk[j].min(t0.elapsed().as_nanos() as f64 - timer_ns);
                    }
                }
            }
            let (nc, nl) = (*cells as f64, list.len() as f64);
            let med = |v: &mut Vec<f64>| {
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                v[v.len() / 2]
            };
            let (c128, cw) = (med(&mut cold[0]), med(&mut cold[1]));
            let (b128, bw): (f64, f64) = (best[0].iter().sum(), best[1].iter().sum());
            report += &format!(
                "  {cat}: problems={} cells={cells} | cold 128-bit ns/cell={:.3} ns/call={:.0} \
                     w128 ns/cell={:.3} ns/call={:.0} 128/w128={:.2} | warm-min 128-bit \
                     ns/cell={:.3} ns/call={:.0} w128 ns/cell={:.3} ns/call={:.0} 128/w128={:.2}\n",
                list.len(),
                c128 / nc,
                c128 / nl,
                cw / nc,
                cw / nl,
                c128 / cw,
                b128 / nc,
                b128 / nl,
                bw / nc,
                bw / nl,
                b128 / bw
            );
        }
    }
    eprint!("{report}");
    if let Some(out) = std::env::var_os("LOSAT_X_DPREPLAY_OUT") {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(out)
            .expect("open LOSAT_X_DPREPLAY_OUT");
        f.write_all(report.as_bytes()).expect("write report");
    }
}
