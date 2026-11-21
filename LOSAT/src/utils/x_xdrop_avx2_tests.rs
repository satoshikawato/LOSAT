//! The 256-bit kernels against the scalar transliteration of the NCBI loops and against the
//! 128-bit kernels: random problems (forward and reversed rows and columns, both lane widths,
//! traceback and score-only, fences, scale 1 and 32, bands wider than 256 cells), and the
//! replay of captured problems (`LOSAT_X_DPREPLAY`).
use super::super::kernel;
use super::super::tests::{random_problem, reference, Problem, Rng};
use super::super::x_capture::load_dir;
// EXPERIMENT: the 128-bit kernels with per-row savings (`x_xdrop_w128.rs`), compared and timed in
// the replay next to the 128-bit and the 256-bit kernels.
use super::super::x_w128 as w128;
use super::*;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kern {
    V128,
    /// the default instances (`EM` of the 8-bit lanes, LOSAT_X_DPEMASK)
    V256,
    /// LOSAT_X_DPEMASK=0: the 8-bit lanes' row end of before (16-bit lanes: the same kernel)
    V256Z,
    W128,
}

/// One call of the chosen kernel; `m8` selects the 8-bit lanes (as `fits_8bit` would).
///
/// # Safety
/// `available()` must be true.
unsafe fn run(
    k: Kern,
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
    match (k, tb, m8) {
        (Kern::V128, true, Some(m)) => {
            kernel::align8_avx::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::V128, false, Some(m)) => {
            kernel::align8_avx::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::V128, true, None) => {
            kernel::align_avx::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::V128, false, None) => {
            kernel::align_avx::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::V256, true, Some(m)) => {
            if sc.x_nib {
                align8::<true, true, true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            } else {
                align8::<true, true, false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
        }
        (Kern::V256, false, Some(m)) => {
            align8::<false, true, false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::V256 | Kern::V256Z, true, None) => {
            if sc.x_nib {
                align16::<true, true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            } else {
                align16::<true, false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
            }
        }
        (Kern::V256 | Kern::V256Z, false, None) => {
            align16::<false, false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::V256Z, true, Some(m)) => {
            if sc.x_nib {
                align8::<true, false, true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            } else {
                align8::<true, false, false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
            }
        }
        (Kern::V256Z, false, Some(m)) => {
            align8::<false, false, false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::W128, true, Some(m)) => {
            w128::run8::<true>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::W128, false, Some(m)) => {
            w128::run8::<false>(q, s, scores, m, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::W128, true, None) => {
            w128::run16::<true>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
        (Kern::W128, false, None) => {
            w128::run16::<false>(q, s, scores, len1, len2, go, ge, x, cf, sc)
        }
    }
}

/// The vector row tables are byte for byte those of `build_row_table`, and are rejected in
/// the same cases.
#[test]
fn row_tables_match_build_row_table() {
    if !available() {
        eprintln!("no AVX2: skipped");
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
    let mut rng = Rng(0x5151_7A7A_0101_2323);
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
        // SAFETY: AVX2 detected above.
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
        // SAFETY: AVX2 detected above.
        let rb = unsafe { build_row(&flat, q, &mut b) };
        assert_eq!(ra, rb, "Flat n={n} q={q}");
        if ra {
            assert_eq!(a, b, "Flat n={n} q={q}");
            accepted += 1;
        }
    }
    assert!(accepted > 1000, "only {accepted} tables accepted");
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
// Random problems run through the scalar loops (`reference`), the 128-bit kernels and the
// 256-bit kernels, with the rows and columns given forward and reversed (as the call sites
// pass them for the left extension), on 16-bit lanes and, where the problem fits, 8-bit lanes,
// with and without traceback. Offsets, score, fence flag, edit script and cell count must be
// equal. `LOSAT_FUZZ_CASES` sets the number of problems.
#[test]
fn avx2_kernels_match_scalar_loops_on_random_problems() {
    if !available() {
        eprintln!("no AVX2: skipped");
        return;
    }
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(600);
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let (mut sc128, mut sc256, mut sc256z) = (
        XdropScratch::new(),
        XdropScratch::new(),
        XdropScratch::new(),
    );
    // [width 8 / 16][tb][reverse]
    let mut handled = [[[0usize; 2]; 2]; 2];
    let (mut fences, mut ops_total, mut declined, mut wide) = (0usize, 0usize, 0usize, 0usize);
    // calls whose rows were moved: [width 8 / 16][128-bit / 256-bit / 256-bit without EM]
    let mut moved = [[0usize; 3]; 2];
    for case in 0..cases {
        // rows over the band (`x_room`): 128 or 320 lanes to start with, or every column
        let lanes = [Some(128), Some(320), None][case % 3];
        sc128.x_test_rows(lanes);
        sc256.x_test_rows(lanes);
        sc256z.x_test_rows(lanes);
        // trace layout (LOSAT_X_DPNIBBLE): four-bit codes or one byte per cell
        let nib = case % 2 == 0;
        (sc128.x_nib, sc256.x_nib, sc256z.x_nib) = (nib, nib, nib);
        let n = if rng.below(2) == 0 { 16 } else { 28 };
        let scale = if rng.below(3) == 0 { 32 } else { 1 };
        let long = rng.below(6) == 0;
        let len = 1 + rng.below(if long { 2500 } else { 300 }) as usize;
        let repeats = rng.below(4) == 0;
        let p = if repeats {
            repeat_problem(&mut rng, n, scale, len)
        } else {
            random_problem(&mut rng, n, scale, len)
        };
        let (go, ge, x) = match rng.below(8) {
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
            // small costs for the 8-bit lanes, gap_extend up to 5 (5 is left to 128 bits)
            _ => (
                rng.below(8) as i32,
                1 + rng.below(5) as i32,
                10 + rng.below(80) as i32,
            ),
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
                    // SAFETY: AVX2 detected above.
                    let r128 = unsafe {
                        run(
                            Kern::V128,
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
                            &mut sc128,
                        )
                    };
                    // SAFETY: AVX2 detected above.
                    let r256 = unsafe {
                        run(
                            Kern::V256,
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
                            &mut sc256,
                        )
                    };
                    // on 16-bit lanes V256Z is the same `align16` call
                    let r256z = if !wide8 {
                        r256
                    } else {
                        // SAFETY: AVX2 detected above.
                        unsafe {
                            run(
                                Kern::V256Z,
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
                                &mut sc256z,
                            )
                        }
                    };
                    assert_eq!(r256, r256z, "{tag}: 256 vs 256 with the E clearing");
                    if wide8 && r256.is_some() {
                        assert_eq!(sc256.ops, sc256z.ops, "{tag}: 256 vs 256z ops");
                    }
                    match (r128, r256) {
                        (Some(a), Some(b)) => {
                            assert_eq!(a, b, "{tag}: 128 vs 256");
                            assert_eq!(sc128.ops, sc256.ops, "{tag}: 128 vs 256 ops");
                        }
                        (None, None) => {}
                        (Some(_), None) => declined += 1,
                        (None, Some(_)) => {
                            panic!("{tag}: 256-bit kernel handled a call the 128-bit declines")
                        }
                    }
                    let Some(r) = r256 else {
                        continue;
                    };
                    assert_eq!(
                        (r.a_offset, r.b_offset, r.score, r.fence_hit),
                        (expected.0, expected.1, expected.2, expected.3),
                        "{tag}: 256 vs scalar"
                    );
                    assert_eq!(sc256.ops, expected.4, "{tag}: 256 vs scalar ops");
                    if !r.fence_hit {
                        assert_eq!(r.cells, expected.5, "{tag}: cells");
                    }
                    handled[usize::from(!wide8)][usize::from(tb)][usize::from(reverse)] += 1;
                    moved[usize::from(!wide8)][0] += usize::from(sc128.x_sh > 0);
                    moved[usize::from(!wide8)][1] += usize::from(sc256.x_sh > 0);
                    // on 16-bit lanes V256Z is the `align16` call of `sc256`
                    let z = if wide8 { &sc256z } else { &sc256 };
                    moved[usize::from(!wide8)][2] += usize::from(z.x_sh > 0);
                    fences += usize::from(r.fence_hit);
                    ops_total += sc256.ops.len();
                    if !wide8 && r.cells > 256 * len1.min(len2) as u64 && len1.min(len2) > 500 {
                        wide += 1;
                    }
                }
            }
        }
    }
    eprintln!(
        "avx2 fuzz: cases={cases} handled[w8/w16][so/tb][fwd/rev]={handled:?} fences={fences} \
             wide={wide} declined={declined} moved[w8/w16][128/256/256z]={moved:?}"
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
        assert!(
            moved.iter().flatten().all(|&n| n > 0),
            "moved rows {moved:?}"
        );
    }
}

/// Replay of captured problems (`LOSAT_X_DPCAPTURE` output): every problem through the scalar
/// loops, the 128-bit and the 256-bit kernels, the 128-bit kernels with per-row savings
/// (`x_xdrop_w128.rs`, "w128") and the `xdrop_align` dispatch, all compared with the recorded
/// result; then the 128-bit, w128 and 256-bit kernels timed over all problems
/// (`LOSAT_X_DPREPLAY_REPS` repetitions, default 7, median), per lane width and traceback /
/// score-only. `LOSAT_X_DPREPLAY_OUT=<file>` appends the table to a file.
#[test]
#[ignore = "needs LOSAT_X_DPREPLAY=<capture directory>"]
fn replay_captured_problems() {
    use std::time::Instant;
    let Some(dir) = std::env::var_os("LOSAT_X_DPREPLAY") else {
        eprintln!("LOSAT_X_DPREPLAY is not set: nothing to replay");
        return;
    };
    assert!(available(), "the replay needs AVX2");
    let dir = std::path::PathBuf::from(dir);
    let probs = load_dir(&dir).expect("read the capture files");
    let reps: usize = std::env::var("LOSAT_X_DPREPLAY_REPS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(7)
        .max(1);
    let (mut sa, mut sb, mut sd, mut sw) = (
        XdropScratch::new(),
        XdropScratch::new(),
        XdropScratch::new(),
        XdropScratch::new(),
    );
    // categories: [tb][8-bit] -> (problem indices, cells)
    let mut cats: [[(Vec<usize>, u64); 2]; 2] = Default::default();
    // [tb][8-bit]: rows and lanes the 256-bit kernel processes (a traceback run of the problem
    // records them; the band evolves the same way without traceback)
    let mut shape = [[(0u64, 0u64); 2]; 2];
    let (mut unhandled, mut declined, mut declined_w) = (0usize, 0usize, 0usize);
    let mut m8s = Vec::with_capacity(probs.len());
    for (i, p) in probs.iter().enumerate() {
        let (q, s, scores) = (p.row_seq(), p.col_seq(), p.scores());
        let m8 = fits_8bit(&scores, p.gap_open, p.gap_extend, p.x_drop);
        m8s.push(m8);
        let Some(want) = p.result else {
            unhandled += 1;
            continue;
        };
        assert_eq!(p.check_fence, p.tb, "problem {i}: check_fence != tb");
        let args = (p.len1, p.len2, p.gap_open, p.gap_extend, p.x_drop);
        // SAFETY: AVX2 checked above.
        let r128 = unsafe {
            run(
                Kern::V128,
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
            )
        };
        assert_eq!(r128, Some(want), "problem {i}: 128-bit vs recorded");
        assert_eq!(sa.ops, p.ops, "problem {i}: 128-bit ops vs recorded");
        // SAFETY: AVX2 checked above.
        let r256 = unsafe {
            run(
                Kern::V256,
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
                &mut sb,
            )
        };
        // SAFETY: AVX2 checked above.
        if unsafe {
            run(
                Kern::V256,
                m8,
                true,
                &q,
                &s,
                &scores,
                args.0,
                args.1,
                args.2,
                args.3,
                args.4,
                p.check_fence,
                &mut sd,
            )
        }
        .is_some()
        {
            let sh = &mut shape[usize::from(p.tb)][usize::from(m8.is_some())];
            let rows = sd.trace.x_test_layout(sd.x_nib);
            sh.0 += rows.len().saturating_sub(1) as u64;
            sh.1 += rows.iter().skip(1).map(|r| r.3 as u64).sum::<u64>();
        }
        match r256 {
            Some(r) => {
                assert_eq!(r, want, "problem {i}: 256-bit vs recorded");
                assert_eq!(sb.ops, p.ops, "problem {i}: 256-bit ops vs recorded");
            }
            None => declined += 1,
        }
        // SAFETY: AVX2 checked above.
        let r256z = unsafe {
            run(
                Kern::V256Z,
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
                &mut sb,
            )
        };
        assert_eq!(
            r256z, r256,
            "problem {i}: 256-bit with the E clearing vs without"
        );
        if r256z.is_some() {
            assert_eq!(
                sb.ops, p.ops,
                "problem {i}: 256-bit (E clearing) ops vs recorded"
            );
        }
        // SAFETY: AVX2 checked above (the 128-bit kernels need a subset of its features).
        let rw = unsafe {
            run(
                Kern::W128,
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
            )
        };
        match rw {
            Some(r) => {
                assert_eq!(r, want, "problem {i}: w128 vs recorded");
                assert_eq!(sw.ops, p.ops, "problem {i}: w128 ops vs recorded");
            }
            None => declined_w += 1,
        }
        assert_eq!(p.replay(&mut sd), Some(want), "problem {i}: dispatch");
        assert_eq!(sd.ops, p.ops, "problem {i}: dispatch ops");
        let prob = Problem {
            q: p.rows.clone(),
            s: p.cols.clone(),
            n: p.n,
            m: p.matrix.clone(),
        };
        let exp = reference(&prob, args.0, args.1, args.2, args.3, args.4, p.tb, true);
        assert_eq!(
            (want.a_offset, want.b_offset, want.score, want.fence_hit),
            (exp.0, exp.1, exp.2, exp.3),
            "problem {i}: recorded vs scalar"
        );
        assert_eq!(p.ops, exp.4, "problem {i}: recorded ops vs scalar");
        if !want.fence_hit {
            assert_eq!(want.cells, exp.5, "problem {i}: cells vs scalar");
        }
        let cat = &mut cats[usize::from(p.tb)][usize::from(m8.is_some())];
        cat.0.push(i);
        cat.1 += want.cells;
    }
    let checked: usize = cats.iter().flatten().map(|c| c.0.len()).sum();
    let mut report = format!(
        "replay {}: problems={} checked={checked} (scalar, 128-bit, 256-bit, w128, dispatch: 0 \
             mismatches) unhandled={unhandled} declined_by_256={declined} \
             declined_by_w128={declined_w} reps={reps}\n",
        dir.display(),
        probs.len()
    );
    // One call as `xdrop_align` makes it: a call the 256-bit or the w128 kernel declines goes to
    // the 128-bit kernel.
    let call = |k: Kern, i: usize, sc: &mut XdropScratch| {
        let p = &probs[i];
        let (q, s, scores) = (p.row_seq(), p.col_seq(), p.scores());
        let go = |k: Kern, sc: &mut XdropScratch| {
            // SAFETY: AVX2 checked above.
            unsafe {
                run(
                    k,
                    m8s[i],
                    p.tb,
                    &q,
                    &s,
                    &scores,
                    p.len1,
                    p.len2,
                    p.gap_open,
                    p.gap_extend,
                    p.x_drop,
                    p.check_fence,
                    sc,
                )
            }
        };
        let r = go(k, sc);
        if r.is_none() && k != Kern::V128 {
            go(Kern::V128, sc)
        } else {
            r
        }
    };
    // Reads a problem's residues and matrix, so that a "warm" call finds them in cache as the
    // search's calls do (the matrix and sequences are in use there).
    let touch = |i: usize| {
        let p = &probs[i];
        let mut acc = 0i64;
        for &b in p.rows.iter().chain(p.cols.iter()) {
            acc += b as i64;
        }
        for &v in &p.matrix {
            acc += v as i64;
        }
        if let Scores::Rows28(m) = p.scores() {
            for &v in m.iter().flatten() {
                acc += v as i64;
            }
        }
        std::hint::black_box(acc);
    };
    // cost of one pair of `Instant::now()` reads, subtracted from the warm per-call times
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
    report += &format!(
        "  cold = problems back to back in capture order; warm = each call timed alone after its \
             residues and matrix were read (timer cost {timer_ns:.1} ns/call subtracted)\n"
    );
    // kernels in report order; index into `times` / `best`
    const KERNS: [Kern; 3] = [Kern::V128, Kern::W128, Kern::V256];
    let kidx = |k: Kern| KERNS.iter().position(|&x| x == k).unwrap();
    for (tb, by_w) in cats.iter().enumerate() {
        for (w8, (list, cells)) in by_w.iter().enumerate() {
            if list.is_empty() {
                continue;
            }
            let (rows, lanes) = shape[tb][w8];
            let cat = format!(
                "{}-{}",
                if tb == 1 { "traceback" } else { "score-only" },
                if w8 == 1 { 8 } else { 16 }
            );
            report += &format!(
                "  {cat}: rows per call {:.1}, cells per row {:.1}, lanes processed per cell \
                     {:.3} (256-bit)\n",
                rows as f64 / list.len() as f64,
                *cells as f64 / rows.max(1) as f64,
                lanes as f64 / *cells as f64
            );
            for warm in [false, true] {
                let mut times = [Vec::new(), Vec::new(), Vec::new()];
                // warm: fastest time of each call over the repetitions
                let mut best = [
                    vec![f64::MAX; list.len()],
                    vec![f64::MAX; list.len()],
                    vec![f64::MAX; list.len()],
                ];
                for rep in 0..reps {
                    // the order of the kernels rotates from one repetition to the next
                    for slot in 0..KERNS.len() {
                        let k = KERNS[(slot + rep) % KERNS.len()];
                        let sc = match k {
                            Kern::V128 => &mut sa,
                            Kern::W128 => &mut sw,
                            Kern::V256 | Kern::V256Z => &mut sb,
                        };
                        let total = if warm {
                            let mut ns = 0f64;
                            let bk = &mut best[kidx(k)];
                            for (j, &i) in list.iter().enumerate() {
                                touch(i);
                                let t0 = Instant::now();
                                std::hint::black_box(call(k, i, sc));
                                let t = t0.elapsed().as_nanos() as f64 - timer_ns;
                                ns += t;
                                bk[j] = bk[j].min(t);
                            }
                            ns
                        } else {
                            let t0 = Instant::now();
                            for &i in list {
                                std::hint::black_box(call(k, i, sc));
                            }
                            t0.elapsed().as_nanos() as f64
                        };
                        times[kidx(k)].push(total);
                    }
                }
                let stat = |v: &mut Vec<f64>| {
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    (v[v.len() / 2], v[0], v[v.len() - 1])
                };
                let (nc, nl) = (*cells as f64, list.len() as f64);
                if warm {
                    let b: Vec<f64> = best.iter().map(|v| v.iter().sum()).collect();
                    report += &format!(
                        "  {cat} warm-min (fastest of {reps} per call, summed): 128-bit \
                             ns/cell={:.3} ns/call={:.0} | w128 ns/cell={:.3} ns/call={:.0} | \
                             256-bit ns/cell={:.3} ns/call={:.0} | 128/w128={:.2} 128/256={:.2}\n",
                        b[0] / nc,
                        b[0] / nl,
                        b[1] / nc,
                        b[1] / nl,
                        b[2] / nc,
                        b[2] / nl,
                        b[0] / b[1],
                        b[0] / b[2]
                    );
                }
                let st: Vec<(f64, f64, f64)> = times.iter_mut().map(stat).collect();
                let mut line = format!(
                    "  {cat} {}: problems={} cells={}",
                    if warm { "warm" } else { "cold" },
                    list.len(),
                    cells
                );
                for (name, &(m, lo, hi)) in ["128-bit", "w128", "256-bit"].iter().zip(&st) {
                    line += &format!(
                        " | {name} ns/cell={:.3} ns/call={:.0} (min {:.3} max {:.3})",
                        m / nc,
                        m / nl,
                        lo / nc,
                        hi / nc
                    );
                }
                line += &format!(
                    " | 128/w128={:.2} 128/256={:.2}\n",
                    st[0].0 / st[1].0,
                    st[0].0 / st[2].0
                );
                report += &line;
            }
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
