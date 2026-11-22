//! Differential test: the SIMD kernels against a plain scalar
//! transliteration of the NCBI loops, on random related sequences.
use super::*;

const MININT: i32 = i32::MIN / 2;

pub(super) struct Problem {
    pub(super) q: Vec<u8>,
    pub(super) s: Vec<u8>,
    pub(super) n: usize,
    pub(super) m: Vec<i32>,
}

impl Problem {
    fn q(&self, a: usize) -> u8 {
        self.q[a - 1]
    }
    fn s(&self, k: usize, len2: usize, sentinel: bool) -> u8 {
        if k >= len2 && !sentinel {
            return 0;
        }
        self.s.get(k).copied().unwrap_or(0)
    }
    fn score(&self, q: u8, s: u8) -> i32 {
        self.m[q as usize * self.n + s as usize]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:862-912,563-636
// The plain scalar transliteration of the C loops that the differential test compares the vector
// kernels with.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:835-959
// (Blast_SemiGappedAlign, score-only) and :500-727 (ALIGN_EX).
pub(super) fn reference(
    p: &Problem,
    len1: usize,
    len2: usize,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    tb: bool,
    sentinel: bool,
) -> (usize, usize, i32, bool, Vec<(u8, u32)>, u64) {
    let goe = gap_open + gap_extend;
    let x = x_drop.max(goe);
    let mut best = vec![MININT; len2 + 4];
    let mut best_gap = vec![MININT; len2 + 4];
    let mut rows: Vec<(usize, Vec<u8>)> = Vec::new();
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
        score -= gap_extend;
        b_size = i + 1;
    }
    rows.push((0, vec![SCRIPT_GAP_IN_A; b_size]));
    let (mut best_score, mut a_off, mut b_off, mut first) = (0i32, 0usize, 0usize, 0usize);
    let mut fence_hit = false;
    let mut cells = 0u64;
    for a in 1..=len1 {
        let orig = first;
        let mut row = Vec::new();
        let qc = p.q(a);
        let mut score_val = MININT;
        let mut gap_row = MININT;
        let mut last = first;
        cells += (b_size - first) as u64;
        for b in orig..b_size {
            let sc = p.s(b, len2, sentinel);
            if tb && sc == FENCE_SENTRY {
                fence_hit = true;
                break;
            }
            let gap_col = best_gap[b];
            let next = best[b].wrapping_add(p.score(qc, sc));
            let mut script = SCRIPT_SUB;
            if score_val < gap_col {
                script = SCRIPT_GAP_IN_B;
                score_val = gap_col;
            }
            if score_val < gap_row {
                script = SCRIPT_GAP_IN_A;
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
                    a_off = a;
                    b_off = b;
                }
                gap_row -= gap_extend;
                let ext = gap_col - gap_extend;
                let open = score_val - goe;
                if ext < open {
                    best_gap[b] = open;
                } else {
                    best_gap[b] = ext;
                    script += SCRIPT_EXTEND_GAP_B;
                }
                if gap_row < open {
                    gap_row = open;
                } else {
                    script += SCRIPT_EXTEND_GAP_A;
                }
                best[b] = score_val;
            }
            score_val = next;
            row.push(script);
        }
        if first == b_size || fence_hit {
            rows.push((orig, row));
            break;
        }
        if last < b_size - 1 {
            b_size = last + 1;
        } else {
            while gap_row >= best_score - x && b_size <= len2 {
                best[b_size] = gap_row;
                best_gap[b_size] = gap_row - goe;
                gap_row -= gap_extend;
                let idx = b_size - orig;
                if idx < row.len() {
                    row[idx] = SCRIPT_GAP_IN_A;
                } else {
                    row.resize(idx, SCRIPT_GAP_IN_A);
                    row.push(SCRIPT_GAP_IN_A);
                }
                b_size += 1;
            }
        }
        rows.push((orig, row));
        if b_size <= len2 {
            best[b_size] = MININT;
            best_gap[b_size] = MININT;
            b_size += 1;
        }
    }
    let mut ops = Vec::new();
    if !tb {
        return (a_off, b_off, best_score, false, ops, cells);
    }
    if best_score <= 0 {
        return (0, 0, 0, fence_hit, ops, cells);
    }
    if fence_hit {
        return (a_off, b_off, best_score, true, ops, cells);
    }
    let (mut a, mut b, mut script) = (a_off, b_off, SCRIPT_SUB);
    while a > 0 || b > 0 {
        let (start, row) = &rows[a];
        let next = row[b - start];
        match script & SCRIPT_OP_MASK {
            x if x == SCRIPT_GAP_IN_A => {
                script = next & SCRIPT_OP_MASK;
                if next & SCRIPT_EXTEND_GAP_A != 0 {
                    script = SCRIPT_GAP_IN_A;
                }
            }
            x if x == SCRIPT_GAP_IN_B => {
                script = next & SCRIPT_OP_MASK;
                if next & SCRIPT_EXTEND_GAP_B != 0 {
                    script = SCRIPT_GAP_IN_B;
                }
            }
            _ => script = next & SCRIPT_OP_MASK,
        }
        let op = script & SCRIPT_OP_MASK;
        if op == SCRIPT_GAP_IN_A {
            b = b.saturating_sub(1);
        } else if op == SCRIPT_GAP_IN_B {
            a = a.saturating_sub(1);
        } else {
            a = a.saturating_sub(1);
            b = b.saturating_sub(1);
        }
        push_ops(&mut ops, op, 1);
    }
    (a_off, b_off, best_score, false, ops, cells)
}

pub(super) struct Rng(pub(super) u64);
impl Rng {
    pub(super) fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    pub(super) fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    pub(super) fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

pub(super) fn random_problem(rng: &mut Rng, n: usize, scale: i32, len: usize) -> Problem {
    // random symmetric-ish matrix with a positive diagonal
    let mut m = vec![0i32; n * n];
    for i in 0..n {
        for j in 0..n {
            let base = if i == j {
                2 + rng.below(8) as i32
            } else {
                -(rng.below(5) as i32)
            };
            let noise = if scale > 1 {
                rng.below(scale as u64 / 2 + 1) as i32 - scale / 4
            } else {
                0
            };
            m[i * n + j] = base * scale + noise;
        }
    }
    // one "minus infinity" row/column, like the gap letter
    for i in 0..n {
        m[i] = i16::MIN as i32;
        m[i * n] = i32::MIN / 2;
    }
    let letters = (n - 1) as u64;
    let q: Vec<u8> = (0..len).map(|_| 1 + rng.below(letters) as u8).collect();
    let (ps, pi) = (rng.unit() * 0.8, rng.unit() * 0.05);
    let mut s = Vec::new();
    let mut i = 0;
    while i < q.len() {
        let u = rng.unit();
        if u < pi {
            for _ in 0..=rng.below(5) {
                s.push(1 + rng.below(letters) as u8);
            }
        } else if u < 2.0 * pi {
            i += 1 + rng.below(5) as usize;
            continue;
        }
        s.push(if rng.unit() < ps {
            1 + rng.below(letters) as u8
        } else {
            q[i]
        });
        i += 1;
    }
    if rng.below(6) == 0 && !s.is_empty() {
        let k = rng.below(s.len() as u64) as usize;
        s[k] = 0;
    }
    Problem { q, s, n, m }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:736-741,766-770
// `Rows` (memory placement only): the first element on a 64-byte boundary at every length
// the scratch takes, the contents kept when the rows grow, as `Vec::resize`.
#[test]
fn rows_start_on_64_byte_boundaries_and_keep_their_values_when_they_grow() {
    fn check<T: Copy + Default + PartialEq + std::fmt::Debug + From<u8>>() {
        let mut r = Rows::<T>::new(1024);
        assert_eq!(r.len(), 1024);
        assert_eq!(r.as_ptr() as usize % 64, 0);
        for (i, x) in r.iter_mut().enumerate() {
            *x = T::from((i % 251) as u8);
        }
        for len in [1030usize, 2100, 4096, 9001] {
            r.resize(len, T::default());
            assert_eq!(r.len(), len);
            assert_eq!(r.as_ptr() as usize % 64, 0, "len {len}");
            assert!(r[..1024]
                .iter()
                .enumerate()
                .all(|(i, &x)| x == T::from((i % 251) as u8)));
            assert!(r[1024..].iter().all(|&x| x == T::default()));
        }
        assert_eq!(Rows::<T>::new(0).len(), 0);
    }
    check::<u8>();
    check::<i16>();
}

// Random problems (matrices of 16 and 28 letters, several gap costs and X-drop values,
// some with a fence sentinel) are run through `reference` and `xdrop_align`. The
// offsets, score, fence flag, edit script and cell count are compared.
#[test]
fn simd_kernels_match_scalar_loops_on_random_problems() {
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1500);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut sc = XdropScratch::new();
    let (mut handled, mut fences, mut ops_total) = (0usize, 0usize, 0usize);
    for case in 0..cases {
        // rows over the band (`x_room`): 128 or 320 lanes to start with, or every column
        sc.x_test_rows([Some(128), Some(320), None][case % 3]);
        // trace layout (LOSAT_X_DPNIBBLE): four-bit codes or one byte per cell
        sc.x_nib = case % 2 == 0;
        let n = if rng.below(2) == 0 { 16 } else { 28 };
        let scale = if rng.below(3) == 0 { 32 } else { 1 };
        let max_len = if rng.below(5) == 0 { 2500 } else { 300 };
        let len = 1 + rng.below(max_len) as usize;
        let mut p = random_problem(&mut rng, n, scale, len);
        let (go, ge, x) = match rng.below(6) {
            0 => (11 * scale, scale, 38 * scale),
            1 => (11 * scale, scale, 64 * scale + rng.below(5) as i32),
            2 => (5 * scale, 2 * scale, 33 * scale),
            3 => (5 * scale, 2 * scale, 110 * scale),
            4 => (
                rng.below(12) as i32 * scale,
                (1 + rng.below(3) as i32) * scale,
                (3 + rng.below(150) as i32) * scale,
            ),
            _ => (9 * scale, scale, 3 * scale),
        };
        let len1 = p.q.len();
        let mut len2 = p.s.len();
        if len2 == 0 {
            continue;
        }
        for tb in [false, true] {
            let mut sentinel = false;
            if tb {
                match rng.below(6) {
                    0 => {
                        p.s.push(FENCE_SENTRY);
                        sentinel = true;
                    }
                    1 if len2 > 3 => {
                        let k = 1 + rng.below(len2 as u64 - 1) as usize;
                        p.s[k] = FENCE_SENTRY;
                    }
                    2 if len2 > 3 => {
                        len2 = 1 + rng.below(len2 as u64 - 1) as usize;
                        sentinel = rng.below(2) == 0;
                    }
                    _ => {}
                }
            }
            let expected = reference(&p, len1, len2, go, ge, x, tb, sentinel);
            let rows = RowSeq::Bytes {
                data: &p.q,
                base: -1,
                step: 1,
            };
            let cols = ColSeq {
                data: &p.s,
                base: 0,
                step: 1,
                zero_from: if sentinel { usize::MAX } else { len2 },
            };
            let scores = Scores::Flat { data: &p.m, n };
            let Some(r) = xdrop_align(
                &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut sc,
            ) else {
                continue;
            };
            handled += 1;
            fences += usize::from(r.fence_hit);
            ops_total += sc.ops.len();
            assert_eq!(
                (r.a_offset, r.b_offset, r.score, r.fence_hit),
                (expected.0, expected.1, expected.2, expected.3),
                "case {case} tb={tb} n={n} scale={scale} go={go} ge={ge} x={x} len1={len1} len2={len2}"
            );
            assert_eq!(sc.ops, expected.4, "case {case} traceback");
            if !r.fence_hit {
                assert_eq!(r.cells, expected.5, "case {case} cell count");
            }
        }
    }
    // On targets without the SIMD kernel nothing is handled.
    if cfg!(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    )) {
        assert!(
            handled > cases,
            "the SIMD kernel handled {handled} of {cases} cases"
        );
        assert!(fences > 0 && ops_total > 0);
    }
}

// No NCBI counterpart: the rows over the band of the vector kernels (`XdropScratch::x_room`).
// Long problems whose band moves right (both scores of the blastn task and a protein matrix,
// score-only and traceback, forward columns) run with rows of 128 lanes to start with that
// move with the band and with every column held: the offsets, score, fence flag, edit
// script and cell count are equal, equal to the scalar loops, and the moved rows stay short.
#[test]
fn band_window_keeps_every_result() {
    let mut rng = Rng(0x5851_F42D_4C95_7F2D);
    let (mut on, mut off) = (XdropScratch::new(), XdropScratch::new());
    let (mut shifted, mut handled) = (0usize, 0usize);
    for case in 0..40 {
        on.x_test_rows(Some(128));
        off.x_test_rows(None);
        let n = if case % 2 == 0 { 16 } else { 28 };
        let len = 1500 + rng.below(2500) as usize;
        let p = random_problem(&mut rng, n, 1, len);
        let (go, ge, x) = if case % 4 < 2 {
            (5, 2, 110)
        } else {
            (11, 1, 38)
        };
        let (len1, len2) = (p.q.len(), p.s.len());
        if len2 == 0 {
            continue;
        }
        let rows = RowSeq::Bytes {
            data: &p.q,
            base: -1,
            step: 1,
        };
        let cols = ColSeq {
            data: &p.s,
            base: 0,
            step: 1,
            zero_from: len2,
        };
        let scores = Scores::Flat { data: &p.m, n };
        for tb in [false, true] {
            let expected = reference(&p, len1, len2, go, ge, x, tb, false);
            let a = xdrop_align(
                &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut on,
            );
            let moved = on.x_sh;
            let rows_len = on.h8a.len().max(on.ha.len());
            let b = xdrop_align(
                &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut off,
            );
            let (Some(a), Some(b)) = (a, b) else {
                continue;
            };
            handled += 1;
            assert_eq!(a, b, "case {case} tb={tb} len1={len1} len2={len2}");
            assert_eq!(on.ops, off.ops, "case {case} tb={tb} traceback");
            assert_eq!(
                (a.a_offset, a.b_offset, a.score, a.fence_hit),
                (expected.0, expected.1, expected.2, expected.3),
                "case {case} tb={tb} against the scalar loops"
            );
            assert_eq!(on.ops, expected.4, "case {case} tb={tb} scalar traceback");
            if moved > 0 {
                shifted += 1;
                assert!(
                    rows_len <= 2048,
                    "case {case}: rows of {rows_len} lanes after moving {moved} columns"
                );
            }
        }
    }
    // every vector kernel moves its rows, the NEON2 ones (LOSAT_X_DPNEON2) too
    if cfg!(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    )) {
        assert!(
            handled > 40 && shifted > 20,
            "handled {handled}, moved {shifted}"
        );
    }
}

// No NCBI counterpart: the traceback rows at four bits per cell (`LOSAT_X_DPNIBBLE`). Long
// tracebacks (both scores of the blastn task and a protein matrix, 8-bit and 16-bit lanes) run
// with four-bit codes and with one byte per cell: the offsets, score, fence flag, edit script
// and cell count are equal, equal to the scalar loops; the rows have the same bands and each
// takes half its cells' bytes, rounded up.
#[test]
fn nibble_trace_keeps_every_result() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let (mut on, mut off) = (XdropScratch::new(), XdropScratch::new());
    (on.x_nib, off.x_nib) = (true, false);
    let mut handled = 0usize;
    for case in 0..40 {
        let n = if case % 2 == 0 { 16 } else { 28 };
        // scale 32 leaves the 8-bit lanes (16-bit kernels)
        let scale = if case % 3 == 2 { 32 } else { 1 };
        let len = 1500 + rng.below(2500) as usize;
        let p = random_problem(&mut rng, n, scale, len);
        let (go, ge, x) = if case % 4 < 2 {
            (5 * scale, 2 * scale, 110 * scale)
        } else {
            (11 * scale, scale, 38 * scale)
        };
        let (len1, len2) = (p.q.len(), p.s.len());
        if len2 == 0 {
            continue;
        }
        let rows = RowSeq::Bytes {
            data: &p.q,
            base: -1,
            step: 1,
        };
        let cols = ColSeq {
            data: &p.s,
            base: 0,
            step: 1,
            zero_from: len2,
        };
        let scores = Scores::Flat { data: &p.m, n };
        let expected = reference(&p, len1, len2, go, ge, x, true, false);
        let a = xdrop_align(
            &rows, &cols, &scores, len1, len2, go, ge, x, true, true, &mut on,
        );
        let b = xdrop_align(
            &rows, &cols, &scores, len1, len2, go, ge, x, true, true, &mut off,
        );
        assert_eq!(a.is_some(), b.is_some(), "case {case}: one layout declined");
        let (Some(a), Some(b)) = (a, b) else {
            continue;
        };
        handled += 1;
        assert_eq!(a, b, "case {case} len1={len1} len2={len2}");
        assert_eq!(on.ops, off.ops, "case {case} traceback");
        assert_eq!(
            (a.a_offset, a.b_offset, a.score, a.fence_hit),
            (expected.0, expected.1, expected.2, expected.3),
            "case {case} against the scalar loops"
        );
        assert_eq!(on.ops, expected.4, "case {case} scalar traceback");
        let (lon, loff) = (on.trace.x_test_layout(true), off.trace.x_test_layout(false));
        let bands =
            |r: &[(usize, usize, usize, usize)]| r.iter().map(|x| (x.2, x.3)).collect::<Vec<_>>();
        assert_eq!(bands(&lon), bands(&loff), "case {case} rows");
        // each row is followed by its record (`x_row_done`)
        let rec = |r: &[(usize, usize, usize, usize)], i: usize| {
            let above = if i == 0 { 0 } else { r[i - 1].2 };
            XTrace::x_test_rec_bytes(r[i].3, r[i].2.wrapping_sub(above) as isize)
        };
        for i in 1..lon.len() {
            if lon[i - 1].0 == lon[i].0 {
                assert_eq!(
                    lon[i].1 - lon[i - 1].1,
                    lon[i - 1].3.div_ceil(2) + rec(&lon, i - 1),
                    "case {case} row bytes"
                );
            }
        }
        assert_eq!(
            off.trace.x_test_used().iter().sum::<usize>(),
            (0..loff.len())
                .map(|i| loff[i].3 + rec(&loff, i))
                .sum::<usize>()
        );
        assert_eq!(
            on.trace.x_test_used().iter().sum::<usize>(),
            (0..lon.len())
                .map(|i| lon[i].3.div_ceil(2) + rec(&lon, i))
                .sum::<usize>(),
            "case {case} trace bytes"
        );
    }
    if cfg!(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    )) {
        assert!(handled > 30, "handled {handled}");
    }
}

// No NCBI counterpart: the rows in chunks (`XTrace`, NCBI `s_GapGetState`). With chunks of
// 512 bytes, rows (each followed by its record: cells, then the step of the first column, one
// byte each below 255, else the long form of nine bytes) are written whole into the free part
// of the chunk being filled, a row that does not fit there goes to the next chunk, a row
// longer than a chunk to a chunk of its own of `length + length/3` bytes, the scratch's own
// chunk keeps its rows and the best row's record when it grows, the walk passes from a chunk's
// first row to the last row of the chunk before, long records (600 cells, steps back, a step
// of 300 forward) are read back, and a purge keeps only the scratch's own chunk, emptied, and
// starts the steps from column 0 again.
#[test]
fn rows_that_do_not_fit_go_to_the_next_chunk() {
    X_TEST_SHIFT.with(|c| c.set(Some(9)));
    let mut t = XTrace::new();
    X_TEST_SHIFT.with(|c| c.set(None));
    assert!(t.x_walk_start().is_none());
    let row = |t: &mut XTrace, n: usize, cells: usize, v: u8| {
        let p = t.x_row(n);
        // SAFETY: `x_row` gave room for `n >= cells` bytes at `p`.
        unsafe { std::ptr::write_bytes(p, v, cells) };
        t.x_row_done(false, v as usize, cells)
    };
    let r1 = row(&mut t, 300, 200, 1);
    t.x_best(r1);
    row(&mut t, 250, 250, 2);
    row(&mut t, 100, 100, 3);
    row(&mut t, 600, 600, 4);
    row(&mut t, 10, 10, 5);
    // records: two bytes, nine for the row of 600 cells
    assert_eq!(
        t.x_test_layout(false),
        vec![
            (0, 0, 1, 200),
            (0, 202, 2, 250),
            (1, 0, 3, 100),
            (2, 0, 4, 600),
            (2, 609, 5, 10)
        ]
    );
    assert_eq!(t.x_test_used(), vec![454, 102, 621]);
    assert!(t.t[0].capacity() <= 512 && t.t[1].capacity() >= 512);
    assert!(t.t[2].capacity() >= 600 + 600 / 3);
    assert!(t.t[0][..200].iter().all(|&b| b == 1) && t.t[0][202..452].iter().all(|&b| b == 2));
    // the best row's record moved with chunk 0 when it grew for row 2
    let (c, base, p) = t.x_walk_start().expect("best row");
    assert_eq!(
        (c, base, p as usize - base as usize, t.best_start),
        (0, t.t[0].as_ptr(), 202, 1)
    );
    // SAFETY: `p` and the records below end records of this call.
    unsafe {
        let at = |c: usize, p: *const u8| p as usize - t.t[c].as_ptr() as usize;
        let (step, cells, rec) = XTrace::x_rec(p);
        assert_eq!((step, cells, at(0, rec)), (1, 200, 200));
        let (c, base, p) = t.x_walk_prev(2);
        let (step, cells, rec) = XTrace::x_rec(p);
        assert_eq!(
            (c, base, step, cells, at(1, rec)),
            (1, t.t[1].as_ptr(), 1, 100, 100)
        );
        let (c, base, p) = t.x_walk_prev(1);
        let (step, cells, rec) = XTrace::x_rec(p);
        assert_eq!(
            (c, base, step, cells, at(0, rec)),
            (0, t.t[0].as_ptr(), 1, 250, 452)
        );
        let (step, cells, rec) = XTrace::x_rec(t.t[2].as_ptr().add(609));
        assert_eq!((step, cells, at(2, rec)), (1, 600, 600));
    }
    // first columns 0, 1, ..., 255, 0, 1, ...: steps of -5 and -255 back (long records)
    for i in 0..400 {
        row(&mut t, 1, 1, i as u8);
    }
    let layout = t.x_test_layout(false);
    assert_eq!(layout.len(), 405);
    assert!(t.t.len() >= 5);
    for i in 0..400 {
        assert_eq!(
            (layout[5 + i].2, layout[5 + i].3),
            (i % 256, 1),
            "row {}",
            5 + i
        );
    }
    t.x_purge();
    assert_eq!((t.t.len(), t.tc), (1, 0));
    assert!(t.x_test_empty() && t.t[0].is_empty() && t.x_walk_start().is_none());
    row(&mut t, 300, 300, 6);
    assert_eq!(t.x_test_layout(false), vec![(0, 0, 6, 300)]);
    // a step of 255 or more forward takes the long form too
    let p = t.x_row(10);
    // SAFETY: `x_row` gave room for 10 bytes at `p`.
    unsafe { std::ptr::write_bytes(p, 7, 10) };
    t.x_row_done(false, 306, 10);
    assert_eq!(
        t.x_test_layout(false),
        vec![(0, 0, 6, 300), (0, 309, 306, 10)]
    );
    assert_eq!(t.x_test_used(), vec![328]);
}

// No NCBI counterpart: the rows in chunks (`XTrace`). Long tracebacks (both scores of the
// blastn task and a protein matrix, 8-bit and 16-bit lanes, four-bit codes and one byte per
// cell) on a scratch with chunks of 512 bytes span many chunks: the offsets, score, fence flag,
// edit script and cell count equal those of a scratch with chunks of 1 MiB (all rows in its own
// chunk) and those of the scalar loops; no row straddles two chunks; after the call
// (`x_trace_release`) the scratch keeps only its own chunks.
#[test]
fn trace_chunks_keep_every_result() {
    X_TEST_SHIFT.with(|c| c.set(Some(9)));
    let mut sc = XdropScratch::new();
    X_TEST_SHIFT.with(|c| c.set(Some(20)));
    let mut own = XdropScratch::new();
    X_TEST_SHIFT.with(|c| c.set(None));
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let (mut handled, mut spanned) = (0usize, 0usize);
    for case in 0..40 {
        let n = if case % 2 == 0 { 16 } else { 28 };
        let scale = if case % 3 == 2 { 32 } else { 1 };
        let len = 1500 + rng.below(2500) as usize;
        let p = random_problem(&mut rng, n, scale, len);
        let (go, ge, x) = if case % 4 < 2 {
            (5 * scale, 2 * scale, 110 * scale)
        } else {
            (11 * scale, scale, 38 * scale)
        };
        let (len1, len2) = (p.q.len(), p.s.len());
        if len2 == 0 {
            continue;
        }
        let rows = RowSeq::Bytes {
            data: &p.q,
            base: -1,
            step: 1,
        };
        let cols = ColSeq {
            data: &p.s,
            base: 0,
            step: 1,
            zero_from: len2,
        };
        let scores = Scores::Flat { data: &p.m, n };
        let expected = reference(&p, len1, len2, go, ge, x, true, false);
        sc.x_nib = case % 8 < 4;
        own.x_nib = sc.x_nib;
        let b = xdrop_align(
            &rows, &cols, &scores, len1, len2, go, ge, x, true, true, &mut own,
        );
        let a = xdrop_align(
            &rows, &cols, &scores, len1, len2, go, ge, x, true, true, &mut sc,
        );
        assert!(own.trace.t.len() <= 1, "case {case}: one chunk");
        let layout = sc.trace.x_test_layout(sc.x_nib);
        let used = sc.trace.x_test_used();
        for &(c, off, _, cells) in &layout {
            let end = off + trace_bytes(sc.x_nib, cells);
            assert!(
                end + 2 <= used[c],
                "case {case}: row past its chunk's used part"
            );
            assert!(end <= 512 || (c > 0 && sc.trace.t[c].capacity() > 512));
        }
        spanned += (sc.trace.t.len() >= 3) as usize;
        sc.x_trace_release();
        assert!(
            sc.trace.t.len() <= 1 && sc.trace.x_test_empty(),
            "case {case}: chunks kept after the call"
        );
        assert!(sc.trace.t.first().map_or(0, |c| c.capacity()) <= 512);
        assert_eq!(a.is_some(), b.is_some(), "case {case}: one side declined");
        let (Some(a), Some(b)) = (a, b) else {
            continue;
        };
        handled += 1;
        assert_eq!(a, b, "case {case} len1={len1} len2={len2}");
        assert_eq!(sc.ops, own.ops, "case {case} traceback");
        assert_eq!(
            (a.a_offset, a.b_offset, a.score, a.fence_hit),
            (expected.0, expected.1, expected.2, expected.3),
            "case {case} against the scalar loops"
        );
        assert_eq!(sc.ops, expected.4, "case {case} scalar traceback");
    }
    if cfg!(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128"),
        target_arch = "aarch64"
    )) {
        assert!(
            handled > 30 && spanned > 30,
            "handled {handled}, spanned {spanned}"
        );
    }
}

// No NCBI counterpart: a BLASTNA-like problem for the 32-bit lanes' tests: row letters 0..3
// (the packed subject of `s_BlastAlignPackedNucl`), column letters 0..3 with some N (14) and,
// in some problems, a gap letter (15), whose row and column hold half of `INT4_MIN` as NCBI's
// matrix does (`blast_stat.c` `BlastScoreBlkNuclMatrixCreate`); `scale` multiplies the scores.
pub(super) fn nucleotide_problem(rng: &mut Rng, scale: i32, len: usize) -> Problem {
    let n = 16;
    let (reward, penalty) = (1 + rng.below(5) as i32, 1 + rng.below(5) as i32);
    let mut m = vec![0i32; n * n];
    for i in 0..n {
        for j in 0..n {
            let v = if i < 4 && j < 4 {
                if i == j {
                    reward
                } else {
                    -penalty
                }
            } else {
                -(1 + rng.below(3) as i32)
            };
            m[i * n + j] = scale * v;
        }
    }
    for i in 0..n {
        m[i * n + 15] = i32::MIN / 2;
        m[15 * n + i] = i32::MIN / 2;
    }
    let q: Vec<u8> = (0..len).map(|_| rng.below(4) as u8).collect();
    let (ps, pi) = (rng.unit() * 0.5, rng.unit() * 0.05);
    let mut s = Vec::new();
    let mut i = 0;
    while i < q.len() {
        let u = rng.unit();
        if u < pi {
            for _ in 0..=rng.below(5) {
                s.push(rng.below(4) as u8);
            }
        } else if u < 2.0 * pi {
            i += 1 + rng.below(5) as usize;
            continue;
        }
        s.push(if rng.unit() < ps {
            rng.below(4) as u8
        } else {
            q[i]
        });
        i += 1;
    }
    if !s.is_empty() {
        for _ in 0..rng.below(3) {
            let k = rng.below(s.len() as u64) as usize;
            s[k] = 14;
        }
        if rng.below(8) == 0 {
            let k = rng.below(s.len() as u64) as usize;
            s[k] = 15;
        }
    }
    Problem { q, s, n, m }
}

// No NCBI counterpart: the row letters of a problem in one of the forms the callers pass
// (`RowSeq`): 0 bytes forward, 1 bytes reversed, 2 NCBI2NA packed forward, 3 packed reversed
// (letters 0..3), with `get(a) == q[a - 1]` (`row_seq`).
fn row_data(q: &[u8], form: usize) -> Vec<u8> {
    match form {
        0 => q.to_vec(),
        1 => q.iter().rev().copied().collect(),
        2 => {
            let mut d = vec![0u8; q.len() / 4 + 2];
            for (i, &c) in q.iter().enumerate() {
                d[i / 4] |= c << (2 * (3 - i % 4));
            }
            d
        }
        _ => {
            let len = q.len();
            let mut d = vec![0u8; len / 4 + 2];
            for a in 1..=len {
                d[(len - a) / 4] |= q[a - 1] << (2 * ((a - 1) % 4));
            }
            d
        }
    }
}

fn row_seq(data: &[u8], form: usize, len1: usize) -> RowSeq<'_> {
    match form {
        0 => RowSeq::Bytes {
            data,
            base: -1,
            step: 1,
        },
        1 => RowSeq::Bytes {
            data,
            base: len1 as isize,
            step: -1,
        },
        2 => RowSeq::PackedFwd {
            data,
            byte_offset: -1,
        },
        _ => RowSeq::PackedRev { data, len: len1 },
    }
}

// No NCBI counterpart: the column letters forward or reversed (`ColSeq`), with
// `get(k) == p.s(k, len2, sentinel)`.
fn col_seq(data: &[u8], rev: bool, zero_from: usize) -> ColSeq<'_> {
    if rev {
        ColSeq {
            data,
            base: data.len() as isize - 1,
            step: -1,
            zero_from,
        }
    } else {
        ColSeq {
            data,
            base: 0,
            step: 1,
            zero_from,
        }
    }
}

// No NCBI counterpart: `R` of the X-drop clamp (head of this file) from every score of the
// matrix, without the `MAX_M` cap of `x_drop_bound`.
fn clamp_bound(p: &Problem, len1: usize, len2: usize, go: i32, ge: i32) -> i64 {
    let max_m = p.m.iter().copied().max().unwrap_or(0).max(0) as i64;
    let goe = (go + ge) as i64;
    max_m * len1.min(len2) as i64 + 3 * goe + ge as i64 * (len1 + len2 + 1) as i64
}

// No NCBI counterpart: gap costs of the 32-bit lanes' tests: the usual ones scaled, random
// ones, and gap extensions the 16-bit lanes cannot hold at any X-drop (`gap_extend * 8 >
// 32767`).
fn wide_gap_costs(rng: &mut Rng, scale: i32) -> (i32, i32) {
    match rng.below(6) {
        0 => (11 * scale, scale),
        1 => (5 * scale, 2 * scale),
        2 => (
            rng.below(12) as i32 * scale,
            (1 + rng.below(3) as i32) * scale,
        ),
        3 => (rng.below(3000) as i32, 4096 + rng.below(3000) as i32),
        4 => (0, scale),
        _ => (2 * scale, scale),
    }
}

// 32-bit lanes (`kernel::align_body32`) and the X-drop clamp (`x_drop_for_lanes`): random
// problems (protein matrices of 28 letters and BLASTNA-like ones of 16, scores scaled up to
// 400 times; rows forward or reversed, one byte per letter or packed; columns forward or
// reversed; Flat matrices and prebuilt tables; score-only and traceback, with fences; both
// trace layouts; rows over the band or holding every column) with X-drops past the 16-bit
// lanes, up to 2^31 - 300, run through `xdrop_align` and through the scalar loops with the
// X-drop as given. The offsets, score, fence flag, edit script and cell count are compared.
// This is the unit test of both sides of `LOSAT_X_DP32` in one process: its `=0` side is the
// scalar ports, of which `reference` is the transliteration, run with the X-drop as given and
// no clamp (the blastn and blastp modules compare their own scalar ports the same way). It also
// counts the 32-bit calls that ended at a fence, that moved their rows with the band
// (`X32_SHIFTS`) and that declined for a -32768 table entry (`X32_AMB`).
#[test]
fn wide_lanes_match_scalar_loops_on_random_problems() {
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(700);
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let mut sc = XdropScratch::new();
    let (mut handled, mut wide, mut narrow, mut fences, mut ops_total) = (0, 0, 0, 0, 0);
    let mut forms = [0usize; 4];
    let mut wide_fences = 0usize;
    let (shifts0, amb0) = (X32_SHIFTS.with(|n| n.get()), X32_AMB.with(|n| n.get()));
    for case in 0..cases {
        sc.x_test_rows([Some(128), Some(320), None][case % 3]);
        sc.x_nib = case % 2 == 0;
        let nucleotide = rng.below(2) == 0;
        let scale = [1, 8, 32, 200, 400][rng.below(if nucleotide { 5 } else { 4 }) as usize];
        let max_len = if rng.below(6) == 0 { 1500 } else { 400 };
        let len = 1 + rng.below(max_len) as usize;
        let mut p = if nucleotide {
            nucleotide_problem(&mut rng, scale, len)
        } else {
            random_problem(&mut rng, 28, scale, len)
        };
        let (go, ge) = wide_gap_costs(&mut rng, scale);
        let len1 = p.q.len();
        let mut len2 = p.s.len();
        if len2 == 0 {
            continue;
        }
        let r = clamp_bound(&p, len1, len2, go, ge);
        let x: i32 = match rng.below(6) {
            0 => 26_000 + rng.below(40_000) as i32,
            1 => (26_000 + rng.below(r.max(1) as u64) as i64).min(i32::MAX as i64 - 300) as i32,
            2 => (r - 1 + rng.below(3) as i64).clamp(1, i32::MAX as i64 - 300) as i32,
            3 => (1 << 20) + rng.below(1 << 30) as i32,
            4 => i32::MAX - 300 - rng.below(1000) as i32,
            _ => 1 + rng.below(i32::MAX as u64 - 300) as i32,
        };
        let row_form = rng.below(if nucleotide { 4 } else { 2 }) as usize;
        let col_rev = rng.below(2) == 1;
        let tables = if rng.below(3) == 0 {
            let m = &p.m;
            let n = p.n;
            build_tables(&|a: u8, b: u8| m[a as usize * n + b as usize], n)
        } else {
            None
        };
        for tb in [false, true] {
            let mut sentinel = false;
            if tb {
                match rng.below(6) {
                    0 => {
                        p.s.push(FENCE_SENTRY);
                        sentinel = true;
                    }
                    1 if len2 > 3 => {
                        let k = 1 + rng.below(len2 as u64 - 1) as usize;
                        p.s[k] = FENCE_SENTRY;
                    }
                    2 if len2 > 3 => {
                        len2 = 1 + rng.below(len2 as u64 - 1) as usize;
                        sentinel = rng.below(2) == 0;
                    }
                    _ => {}
                }
            }
            let expected = reference(&p, len1, len2, go, ge, x, tb, sentinel);
            let rd = row_data(&p.q, row_form);
            let rows = row_seq(&rd, row_form, len1);
            let cd: Vec<u8> = if col_rev {
                p.s.iter().rev().copied().collect()
            } else {
                p.s.clone()
            };
            let cols = col_seq(&cd, col_rev, if sentinel { usize::MAX } else { len2 });
            for a in 1..=len1 {
                assert_eq!(rows.get(a), Some(p.q(a)), "case {case} row form {row_form}");
            }
            for k in 0..=len2 + 1 {
                assert_eq!(cols.get(k), p.s(k, len2, sentinel), "case {case} columns");
            }
            let scores = match &tables {
                Some(t) => Scores::Tables { t, n: p.n },
                None => Scores::Flat { data: &p.m, n: p.n },
            };
            let x_lanes = x_drop_for_lanes(&scores, len1, len2, go, ge, x);
            assert!(x_lanes <= x);
            let runs = X32_RUNS.with(|n| n.get());
            let got = xdrop_align(
                &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut sc,
            );
            let took_wide = X32_RUNS.with(|n| n.get()) > runs;
            if x_kernels_run() {
                assert_eq!(
                    took_wide,
                    x_dp32() && !fits_16bit(go, ge, x_lanes),
                    "case {case}: the 32-bit lanes run exactly where the 16-bit lanes end"
                );
            }
            let Some(got) = got else {
                continue;
            };
            handled += 1;
            if took_wide {
                wide += 1;
                forms[row_form] += 1;
                wide_fences += usize::from(got.fence_hit);
            } else if !fits_16bit(go, ge, x) {
                // a large X-drop that the clamp brought into the 16-bit or 8-bit lanes
                narrow += 1;
            }
            fences += usize::from(got.fence_hit);
            ops_total += sc.ops.len();
            assert_eq!(
                (got.a_offset, got.b_offset, got.score, got.fence_hit),
                (expected.0, expected.1, expected.2, expected.3),
                "case {case} tb={tb} nucleotide={nucleotide} scale={scale} go={go} ge={ge} \
                     x={x} (lanes {x_lanes}, R {r}) len1={len1} len2={len2} rows {row_form} \
                     reversed columns {col_rev} tables {}",
                tables.is_some()
            );
            assert_eq!(sc.ops, expected.4, "case {case} tb={tb} traceback");
            if !got.fence_hit {
                assert_eq!(got.cells, expected.5, "case {case} tb={tb} cell count");
            }
        }
    }
    if x_kernels_run() && x_dp32() {
        let shifts = X32_SHIFTS.with(|n| n.get()) - shifts0;
        let amb = X32_AMB.with(|n| n.get()) - amb0;
        assert!(
            wide > cases / 2 && narrow > cases / 20 && fences > 0 && ops_total > 0,
            "handled {handled}, 32-bit {wide}, clamped to 16/8-bit {narrow}, fences {fences}"
        );
        assert!(
            wide_fences > 0 && shifts > 0 && amb > 0,
            "32-bit fences {wide_fences}, row moves {shifts}, -32768 declines {amb}"
        );
        assert!(
            forms.iter().all(|&f| f > 10),
            "32-bit calls by row form {forms:?}"
        );
    }
}

// The 32-bit lanes on their own (`align_32bit`, no clamp) at X-drops up to the largest they
// hold, `i32::MAX - goe - 3 * MAX_M - 1` (every row rebases there), compared with the scalar
// loops; one past it they decline.
#[test]
fn wide_lanes_hold_x_drops_up_to_their_limit() {
    let mut rng = Rng(0x94D0_49BB_1331_11EB);
    let mut sc = XdropScratch::new();
    let (mut handled, mut top_handled) = (0usize, 0usize);
    for case in 0..80 {
        sc.x_nib = case % 2 == 0;
        let nucleotide = case % 3 == 0;
        let scale = [1, 32, 200][case % 3];
        let len = 1 + rng.below(if case % 5 == 0 { 1200 } else { 300 }) as usize;
        let p = if nucleotide {
            nucleotide_problem(&mut rng, scale, len)
        } else {
            random_problem(&mut rng, 28, scale, len)
        };
        let (go, ge) = wide_gap_costs(&mut rng, scale);
        let (len1, len2) = (p.q.len(), p.s.len());
        if len2 == 0 {
            continue;
        }
        let top = (VMAX32 - (go + ge) as i64 - 3 * MAX_M as i64 - 1) as i32;
        let x = match case % 4 {
            0 => top,
            1 => top - rng.below(1 << 20) as i32,
            2 => 26_000 + rng.below(200_000) as i32,
            _ => top + 1,
        };
        let row_form = rng.below(if nucleotide { 4 } else { 2 }) as usize;
        let rd = row_data(&p.q, row_form);
        let rows = row_seq(&rd, row_form, len1);
        let cols = col_seq(&p.s, false, len2);
        let scores = Scores::Flat { data: &p.m, n: p.n };
        for tb in [false, true] {
            let got = align_32bit(
                &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut sc,
            );
            if x > top {
                assert!(got.is_none(), "case {case}: past the 32-bit lanes");
                continue;
            }
            let Some(got) = got else {
                continue;
            };
            handled += 1;
            top_handled += usize::from(x == top);
            let expected = reference(&p, len1, len2, go, ge, x, tb, false);
            assert_eq!(
                (got.a_offset, got.b_offset, got.score, got.fence_hit),
                (expected.0, expected.1, expected.2, expected.3),
                "case {case} tb={tb} go={go} ge={ge} x={x} len1={len1} len2={len2}"
            );
            assert_eq!(sc.ops, expected.4, "case {case} tb={tb} traceback");
            assert_eq!(got.cells, expected.5, "case {case} tb={tb} cell count");
        }
    }
    if x_kernels_run() {
        assert!(
            handled > 60 && top_handled > 15,
            "handled {handled}, at the limit {top_handled}"
        );
    }
}

// The X-drop clamp (head of this file) on the scalar loops themselves: for random problems
// (protein and BLASTNA-like matrices with their minus-infinity letters, fences, cut columns),
// the offsets, score, fence flag, edit script and cells with an X-drop of `R`, `R` plus a
// little and 2^31 - 300 are the same; and `x_drop_bound` is `R`.
#[test]
fn x_drop_clamp_keeps_every_scalar_result() {
    let mut rng = Rng(0xBF58_476D_1CE4_E5B9);
    let mut compared = 0usize;
    for case in 0..400 {
        let nucleotide = case % 2 == 0;
        let scale = [1, 8, 32, 200][rng.below(4) as usize];
        let len = 1 + rng.below(if case % 7 == 0 { 800 } else { 200 }) as usize;
        let mut p = if nucleotide {
            nucleotide_problem(&mut rng, scale, len)
        } else {
            random_problem(&mut rng, 28, scale, len)
        };
        let (go, ge) = wide_gap_costs(&mut rng, scale);
        let len1 = p.q.len();
        let mut len2 = p.s.len();
        if len2 == 0 {
            continue;
        }
        for tb in [false, true] {
            let mut sentinel = false;
            if tb {
                match rng.below(4) {
                    0 => {
                        p.s.push(FENCE_SENTRY);
                        sentinel = true;
                    }
                    1 if len2 > 3 => {
                        let k = 1 + rng.below(len2 as u64 - 1) as usize;
                        p.s[k] = FENCE_SENTRY;
                    }
                    2 if len2 > 3 => {
                        len2 = 1 + rng.below(len2 as u64 - 1) as usize;
                        sentinel = rng.below(2) == 0;
                    }
                    _ => {}
                }
            }
            let r = clamp_bound(&p, len1, len2, go, ge);
            if p.m.iter().all(|&v| v <= MAX_M) {
                let scores = Scores::Flat { data: &p.m, n: p.n };
                assert_eq!(x_drop_bound(&scores, len1, len2, go, ge), r);
            }
            if r > i32::MAX as i64 - 2000 {
                continue;
            }
            let at_r = reference(&p, len1, len2, go, ge, r as i32, tb, sentinel);
            for x in [r as i32 + 1 + rng.below(1000) as i32, i32::MAX - 300] {
                let above = reference(&p, len1, len2, go, ge, x, tb, sentinel);
                assert_eq!(
                    above, at_r,
                    "case {case} tb={tb} go={go} ge={ge} R={r} x={x} len1={len1} len2={len2}"
                );
                compared += 1;
            }
        }
    }
    assert!(compared > 1000, "compared {compared}");
}

// A traceback and a score-only call with a large X-drop take the 32-bit lanes, not the scalar
// port (`xdrop_align` answers, and `X32_RUNS` counts the call), with the scalar loops'
// results; a small problem with the same X-drops takes the 16-bit or 8-bit lanes (the clamp).
#[test]
fn large_x_drop_takes_the_wide_lanes() {
    let mut rng = Rng(0x2F2E_1A3B_99C4_0D17);
    let mut sc = XdropScratch::new();
    for (len, wide) in [(900usize, true), (20, false)] {
        let mut p = random_problem(&mut rng, 28, 32, len);
        // no minus-infinity column letter (the 32-bit lanes would decline it)
        for c in p.s.iter_mut().filter(|c| **c == 0) {
            *c = 1;
        }
        let (len1, len2) = (p.q.len(), p.s.len());
        let rows = RowSeq::Bytes {
            data: &p.q,
            base: -1,
            step: 1,
        };
        let cols = col_seq(&p.s, false, len2);
        let scores = Scores::Flat { data: &p.m, n: p.n };
        for x in [110_903, 1 << 30, i32::MAX - 300] {
            for tb in [false, true] {
                let runs = X32_RUNS.with(|n| n.get());
                let got = xdrop_align(
                    &rows,
                    &cols,
                    &scores,
                    len1,
                    len2,
                    11 * 32,
                    32,
                    x,
                    tb,
                    tb,
                    &mut sc,
                );
                let expected = reference(&p, len1, len2, 11 * 32, 32, x, tb, false);
                // past the 16-bit lanes with `LOSAT_X_DP32=0` (or without vector kernels):
                // the scalar port
                if !x_kernels_run() || !x_dp32() {
                    assert!(got.is_none());
                    continue;
                }
                let got = got.expect("a vector kernel answers");
                assert_eq!(
                    X32_RUNS.with(|n| n.get()) - runs,
                    usize::from(wide),
                    "len {len} x {x} tb {tb}: the 32-bit lanes ran"
                );
                assert_eq!(
                    (got.a_offset, got.b_offset, got.score, got.fence_hit),
                    (expected.0, expected.1, expected.2, expected.3)
                );
                assert_eq!(sc.ops, expected.4);
                assert_eq!(got.cells, expected.5);
            }
        }
    }
}

// The decline for a -32768 table entry (`x32_amb_of`) is what keeps such a call exact: with a
// gap opening of 40,000 and an X-drop of 1e6, a BLASTNA gap letter (half of `INT4_MIN`, -32768 in
// the 16-bit tables) third in the columns gives the diagonal `H - 2^30` in the scalar loops, so
// the cell takes a gap (`H - 40,001`) and the 37 matches after it (+37,000) do not reach the
// best again, while with -32768 in its place they would (the scalar loops on the clamped
// matrix give another result). The 32-bit lanes take the call (`X32_RUNS`), decline it
// (`X32_AMB`), and `xdrop_align` leaves it to the scalar port (None).
#[test]
fn wide_lanes_decline_a_clamped_table_entry() {
    let n = 16;
    let mut m = vec![-1000i32; n * n];
    for i in 0..4 {
        m[i * n + i] = 1000;
    }
    for i in 0..n {
        m[i * n + 15] = i32::MIN / 2;
        m[15 * n + i] = i32::MIN / 2;
    }
    let q: Vec<u8> = (0..40usize).map(|i| ((i * 7 + i / 3) % 4) as u8).collect();
    let mut s = q.clone();
    s[2] = 15;
    let p = Problem { q, s, n, m };
    let (len1, len2) = (p.q.len(), p.s.len());
    let (go, ge, x) = (40_000, 1, 1_000_000);
    let rows = RowSeq::Bytes {
        data: &p.q,
        base: -1,
        step: 1,
    };
    let cols = col_seq(&p.s, false, len2);
    let scores = Scores::Flat { data: &p.m, n };
    let mut sc = XdropScratch::new();
    // the scalar loops: the alignment ends before the gap letter (score 2,000); with -32768 in
    // its place it would run on through it
    let expected = reference(&p, len1, len2, go, ge, x, true, false);
    let mut clamped = p.m.clone();
    for v in clamped.iter_mut() {
        *v = (*v).max(i16::MIN as i32);
    }
    let as_clamped = Problem {
        m: clamped,
        q: p.q.clone(),
        s: p.s.clone(),
        n,
    };
    assert_ne!(
        reference(&as_clamped, len1, len2, go, ge, x, true, false),
        expected,
        "the clamped entry would change the result"
    );
    for tb in [false, true] {
        let (runs, amb) = (X32_RUNS.with(|n| n.get()), X32_AMB.with(|n| n.get()));
        let direct = align_32bit(
            &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut sc,
        );
        let dispatched = xdrop_align(
            &rows, &cols, &scores, len1, len2, go, ge, x, tb, tb, &mut sc,
        );
        assert!(direct.is_none() && dispatched.is_none(), "tb {tb}");
        if x_kernels_run() {
            let runs = X32_RUNS.with(|n| n.get()) - runs;
            let amb = X32_AMB.with(|n| n.get()) - amb;
            // the direct call always; the dispatched one with `LOSAT_X_DP32` on
            assert_eq!((runs, amb), (1 + usize::from(x_dp32()), runs), "tb {tb}");
        }
    }
}

// No NCBI counterpart: the 32-bit operations of the `arch` layer against their definitions on
// random lanes, alike on x86_64, wasm32 simd128 and aarch64 NEON (`subsu32` on lanes in
// `0 ..= i32::MAX`, where the kernels use it).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128"),
    target_arch = "aarch64"
))]
#[test]
fn arch_32bit_lanes_match_their_definitions() {
    // the x86_64 `arch` layer runs inside the kernels' SSE4.1 entries (`x86_entry32!`)
    #[cfg(target_arch = "x86_64")]
    if !is_x86_feature_detected!("sse4.1") {
        return;
    }
    use super::arch::*;
    let mut rng = Rng(0x6A09_E667_F3BC_C909);
    let lane = |rng: &mut Rng| -> i32 {
        match rng.below(7) {
            0 => 0,
            1 => i32::MAX,
            2 => i32::MIN,
            3 => -1,
            4 => rng.below(70_000) as i32 - 35_000,
            5 => i16::MIN as i32 - rng.below(3) as i32 + 1,
            _ => rng.next() as i32,
        }
    };
    for _ in 0..20_000 {
        let a: [i32; 4] = core::array::from_fn(|_| lane(&mut rng));
        let b: [i32; 4] = core::array::from_fn(|_| lane(&mut rng));
        let h: [i16; 8] = core::array::from_fn(|_| lane(&mut rng) as i16);
        // SAFETY: SSE4.1 was detected above (x86_64); simd128 and NEON are compile-time
        // features of the wasm32 and aarch64 builds; each load and store covers one array of
        // 16 bytes.
        unsafe {
            let w = |v: V| {
                let mut o = [0i32; 4];
                store(o.as_mut_ptr(), v);
                o
            };
            let w16 = |v: V| {
                let mut o = [0i16; 8];
                store(o.as_mut_ptr(), v);
                o
            };
            let (va, vb, vh) = (load(a.as_ptr()), load(b.as_ptr()), load(h.as_ptr()));
            let e: [i32; 4] = core::array::from_fn(|i| a[i].wrapping_add(b[i]));
            assert_eq!(w(add32(va, vb)), e, "add32");
            let e: [i32; 4] = core::array::from_fn(|i| a[i].max(b[i]));
            assert_eq!(w(max32(va, vb)), e, "max32");
            let e: [i32; 4] = core::array::from_fn(|i| -i32::from(a[i] > b[i]));
            assert_eq!(w(gt_i32(va, vb)), e, "gt_i32");
            let (pa, pb) = (a.map(|v| v & i32::MAX), b.map(|v| v & i32::MAX));
            let e: [i32; 4] = core::array::from_fn(|i| (pa[i] - pb[i]).max(0));
            let got = subsu32(load(pa.as_ptr()), load(pb.as_ptr()));
            assert_eq!(w(got), e, "subsu32");
            assert_eq!(w(splat_i32(a[1])), [a[1]; 4], "splat_i32");
            assert_eq!(w(bcast_last32(va)), [a[3]; 4], "bcast_last32");
            assert_eq!(
                (lane0_32(va), lane3_32(va)),
                (a[0], a[3]),
                "lane0_32, lane3_32"
            );
            let e: [i32; 4] = core::array::from_fn(|i| h[i] as i32);
            assert_eq!(w(widen_lo16(vh)), e, "widen_lo16");
            let e: [i32; 4] = core::array::from_fn(|i| h[4 + i] as i32);
            assert_eq!(w(widen_hi16(vh)), e, "widen_hi16");
            let e: [i16; 8] = core::array::from_fn(|i| {
                let v = if i < 4 { a[i] } else { b[i - 4] };
                v.clamp(i16::MIN as i32, i16::MAX as i32) as i16
            });
            assert_eq!(w16(packs32(va, vb)), e, "packs32");
            for (n, got) in [(1, bsl4(va)), (2, bsl8(va))] {
                let e: [i32; 4] = core::array::from_fn(|i| if i >= n { a[i - n] } else { 0 });
                assert_eq!(w(got), e, "bsl{}", 4 * n);
            }
        }
    }
}

// No NCBI counterpart: `nibs` of the `arch` layer (`LOSAT_X_DPNIBBLE`) against its definition
// on random codes: byte k of the low half = code 2k | code 2k+1 << 4.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128"),
    target_arch = "aarch64"
))]
#[test]
fn nibs_pairs_adjacent_codes() {
    // the x86_64 `arch` layer runs inside the kernels' SSSE3 / SSE4.1 entries (`x86_entry!`)
    #[cfg(target_arch = "x86_64")]
    if !is_x86_feature_detected!("ssse3") {
        return;
    }
    let mut rng = Rng(0x0123_4567_89AB_CDEF);
    for _ in 0..1000 {
        let mut v = [0u8; 16];
        for b in v.iter_mut() {
            *b = rng.below(16) as u8;
        }
        let mut out = [0xAAu8; 16];
        // SAFETY: SSSE3 was detected above (x86_64); simd128 and NEON are compile-time features
        // of the wasm32 and aarch64 builds; the pointers cover 16 bytes.
        unsafe {
            super::arch::store64(
                out.as_mut_ptr(),
                super::arch::nibs(super::arch::load(v.as_ptr())),
            );
        }
        for k in 0..8 {
            assert_eq!(out[k], v[2 * k] | v[2 * k + 1] << 4, "byte {k} of {v:?}");
        }
        assert!(
            out[8..].iter().all(|&b| b == 0xAA),
            "store64 wrote past 8 bytes"
        );
    }
}
