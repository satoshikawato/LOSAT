use super::*;
use crate::algorithm::blastp::gapalign::GapAlignScratch;

// No NCBI counterpart: the CPU checks of the tests' rows.
/// Whether the CPU runs `rows` (the checks `x_rows` makes).
fn runs_here(rows: XRows) -> bool {
    match rows {
        XRows::Scalar => true,
        #[cfg(target_arch = "x86_64")]
        XRows::Avx2 => cpu_has_avx2(),
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        XRows::W128 => cpu_has_w128(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:1028-1038 (as
// `super::restricted_fast`)
/// `super::restricted_fast` for rows the CPU runs (asserted): the tests call it by this name.
fn restricted_fast(
    query: &[u8],
    query_base: usize,
    subject: &[u8],
    subject_base: usize,
    len1: usize,
    len2: usize,
    score_matrix: BlastpScoreMatrix<'_>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
    sc: &mut XRestrictedScratch,
    rows: XRows,
    narrow: bool,
    band: bool,
) -> Option<(usize, usize, i32)> {
    assert!(runs_here(rows), "{rows:?} on this CPU");
    // SAFETY: the CPU runs `rows` (asserted above).
    unsafe {
        super::restricted_fast(
            query,
            query_base,
            subject,
            subject_base,
            len1,
            len2,
            score_matrix,
            gap_open,
            gap_extend,
            x_drop,
            reverse,
            sc,
            rows,
            narrow,
            band,
        )
    }
}

/// xorshift64*: a fixed generator, so every failure can be replayed.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn reference(
    q: &[u8],
    qb: usize,
    s: &[u8],
    sb: usize,
    len1: usize,
    len2: usize,
    m: BlastpScoreMatrix<'_>,
    go: i32,
    ge: i32,
    x: i32,
    rev: bool,
    scratch: &mut GapAlignScratch,
) -> (usize, usize, i32) {
    restricted_gapped_align_protein(
        q,
        qb,
        s,
        sb,
        len1,
        len2,
        m,
        go,
        ge,
        x,
        rev,
        &mut scratch.dp_mem,
        &mut scratch.dp_mem_alloc,
    )
}

/// The rows this machine runs, each with its `narrow` (LOSAT_X_BPGAP16): the scalar rows; with
/// AVX2, BMI1 and LZCNT the AVX2 rows without and with the sixteen-bit rows; with SSSE3 and
/// SSE4.1 (always on wasm32 simd128) the 128-bit rows (LOSAT_X_BPGAP128), sixteen-bit rows
/// only.
fn kernels() -> Vec<(XRows, bool)> {
    #[cfg_attr(
        not(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        )),
        expect(
            unused_mut,
            reason = "the vector rows exist on x86_64 and wasm32 simd128 only"
        )
    )]
    let mut kernels = vec![(XRows::Scalar, false)];
    #[cfg(target_arch = "x86_64")]
    if cpu_has_avx2() {
        kernels.extend([(XRows::Avx2, false), (XRows::Avx2, true)]);
    }
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if cpu_has_w128() {
        kernels.push((XRows::W128, true));
    }
    kernels
}

/// LOSAT_X_BPGAP128: the 128-bit rows may leave the problem to the reference function (None)
/// only where its costs fail `fits16` or its score `best` passes `limit16` (they stop at the
/// start of a row whose best score is past it).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn w128_may_decline(m: BlastpScoreMatrix<'_>, go: i32, ge: i32, x: i32, best: i32) -> bool {
    let tables = x_score_tables(m).expect("tables");
    let p = XParams {
        x_dropoff: x.max(go + ge),
        gap_open_extend: go + ge,
        gap_extend: ge,
    };
    !fits16(&p, tables) || best > limit16(&p, tables)
}

/// `got`, the result of `rows`, is `want`; the 128-bit rows may instead give None where
/// `w128_may_decline` allows it.
fn assert_rows_result(
    rows: XRows,
    got: Option<(usize, usize, i32)>,
    want: (usize, usize, i32),
    m: BlastpScoreMatrix<'_>,
    (go, ge, x): (i32, i32, i32),
    case: &str,
) {
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if rows == XRows::W128 && got.is_none() {
        assert!(
            w128_may_decline(m, go, ge, x, want.2),
            "{case}: the 128-bit rows declined"
        );
        return;
    }
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    let _ = (rows, m, go, ge, x);
    assert_eq!(got, Some(want), "{case}");
}

/// Every kernel this machine runs (`kernels`) against the reference.
fn check(
    q: &[u8],
    qb: usize,
    s: &[u8],
    sb: usize,
    len1: usize,
    len2: usize,
    m: BlastpScoreMatrix<'_>,
    go: i32,
    ge: i32,
    x: i32,
    rev: bool,
    scratch: &mut GapAlignScratch,
) -> (usize, usize, i32) {
    let want = reference(q, qb, s, sb, len1, len2, m, go, ge, x, rev, scratch);
    for (rows, narrow) in kernels() {
        // LOSAT_X_BPGAPBAND on a scratch of its own, so that its arrays start at the first
        // rows' cells and double where the rows stop (the shared scratch holds the largest
        // `len2` so far).
        for band in [false, true] {
            let mut own = XRestrictedScratch::new();
            let sc = if band {
                &mut own
            } else {
                &mut scratch.x_restricted
            };
            let got = restricted_fast(
                q, qb, s, sb, len1, len2, m, go, ge, x, rev, sc, rows, narrow, band,
            );
            let case = format!(
                "rows={rows:?} narrow={narrow} band={band} len1={len1} len2={len2} qb={qb} sb={sb} go={go} ge={ge} x={x} rev={rev}"
            );
            assert_rows_result(rows, got, want, m, (go, ge, x), &case);
        }
    }
    want
}

/// A random protein over `alphabet` letters (codes 1..=alphabet), with a few gap letters (0)
/// and stop letters (25) when `odd` is set.
fn random_seq(rng: &mut Rng, len: usize, alphabet: u64, odd: bool) -> Vec<u8> {
    (0..len)
        .map(|_| {
            if odd && rng.below(50) == 0 {
                [0u8, 25, 21, 23, 24, 26, 27][rng.below(7) as usize]
            } else {
                1 + rng.below(alphabet) as u8
            }
        })
        .collect()
}

/// `src` with substitutions and insertions / deletions, so that the alignment runs long and
/// uses gaps.
fn mutate(rng: &mut Rng, src: &[u8], rate: u64) -> Vec<u8> {
    let mut out = Vec::with_capacity(src.len() + 16);
    for &c in src {
        match rng.below(100) {
            r if r < rate => out.push(1 + rng.below(24) as u8),
            r if r < rate + 2 => {}
            r if r < rate + 4 => {
                out.push(c);
                for _ in 0..1 + rng.below(4) {
                    out.push(1 + rng.below(24) as u8);
                }
            }
            _ => out.push(c),
        }
    }
    out
}

/// The two calls of `s_BlastProtGappedAlignment` (left, reverse; right, forward) at a seed
/// inside `q` and `s`, with the arguments `blastp_score_only_gapped_alignment_with_scratch`
/// gives them (no subject shift).
fn check_seed(
    q: &[u8],
    s: &[u8],
    q_start: usize,
    s_start: usize,
    m: BlastpScoreMatrix<'_>,
    go: i32,
    ge: i32,
    x: i32,
    scratch: &mut GapAlignScratch,
) {
    let q_length = q_start + 1;
    let s_length = s_start + 1;
    check(q, 0, s, 0, q_length, s_length, m, go, ge, x, true, scratch);
    if q_length < q.len() && s_length < s.len() {
        check(
            q,
            q_start,
            s,
            s_start,
            q.len() - q_length,
            s.len() - s_length,
            m,
            go,
            ge,
            x,
            false,
            scratch,
        );
    }
}

const MATRICES: [BlastpScoreMatrix<'static>; 4] = [
    BlastpScoreMatrix::Blosum62,
    BlastpScoreMatrix::Standard(ScoringMatrix::Pam30),
    BlastpScoreMatrix::Standard(ScoringMatrix::Blosum45),
    BlastpScoreMatrix::Standard(ScoringMatrix::Blosum80),
];

/// The sixteen-bit rows this machine runs (`narrow16`): the AVX2 rows (AVX2, BMI1 and LZCNT)
/// and the 128-bit rows (LOSAT_X_BPGAP128: SSSE3 and SSE4.1; wasm32 simd128).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn narrow_kernels() -> Vec<XRows> {
    let mut kernels = Vec::new();
    #[cfg(target_arch = "x86_64")]
    if cpu_has_avx2() {
        kernels.push(XRows::Avx2);
    }
    if cpu_has_w128() {
        kernels.push(XRows::W128);
    }
    kernels
}

/// LOSAT_X_BPGAP16: the sixteen-bit rows of `rows` (`XRows::Avx2`: `restricted_rows16_avx2`,
/// `XRows::W128`: `restricted_rows16_w128`) alone, set up as `restricted_fast` sets them up.
/// None: `fits16` does not hold (or the problem is not taken at all); Some(None): the rows
/// stopped at the 16-bit bound. Only for the rows of `narrow_kernels`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn narrow16(
    rows: XRows,
    q: &[u8],
    qb: usize,
    s: &[u8],
    sb: usize,
    len1: usize,
    len2: usize,
    m: BlastpScoreMatrix<'_>,
    go: i32,
    ge: i32,
    x: i32,
    rev: bool,
    scratch: &mut GapAlignScratch,
) -> Option<Option<(usize, usize, i32)>> {
    assert!(narrow_kernels().contains(&rows) && len1 > 0 && len2 > 0);
    // The 32-bit rows alone size `h`, `e10` and `col`.
    let sc = &mut scratch.x_restricted;
    restricted_fast(
        q,
        qb,
        s,
        sb,
        len1,
        len2,
        m,
        go,
        ge,
        x,
        rev,
        sc,
        XRows::Scalar,
        false,
        false,
    )?;
    let tables = x_score_tables(m)?;
    let p = XParams {
        x_dropoff: x.max(go + ge),
        gap_open_extend: go + ge,
        gap_extend: ge,
    };
    if !fits16(&p, tables) {
        return None;
    }
    let cells = PAD + len2 + 2 + TAIL;
    for h in sc.h16.iter_mut() {
        if h.len() < cells {
            h.resize(cells, i16::MIN);
        }
    }
    if sc.e16.len() < cells {
        sc.e16.resize(cells, i16::MIN);
    }
    Some(
        match rows {
            // SAFETY: AVX2, BMI1 and LZCNT (`narrow_kernels`); the arrays are sized for `len2`;
            // `fits16`.
            #[cfg(target_arch = "x86_64")]
            XRows::Avx2 => unsafe {
                restricted_rows16_avx2(q, qb, s, sb, len1, len2, rev, tables, &p, sc)
            },
            // SAFETY: SSSE3 and SSE4.1 on x86_64 (`narrow_kernels`); the arrays are sized for
            // `len2`; `fits16`.
            XRows::W128 => unsafe {
                restricted_rows16_w128(q, qb, s, sb, len1, len2, rev, tables, &p, sc)
            },
            XRows::Scalar => unreachable!("the scalar rows have no sixteen-bit form"),
        }
        // Arrays of `len2` cells: the rows never reach the end of `col`.
        .inspect_err(|&stop| assert_eq!(stop, XStop::Decline))
        .ok(),
    )
}

/// LOSAT_X_BPGAP16 at its bound: long alignments of W against W (BLOSUM62's largest entry,
/// 11) score just below `limit16`, where the sixteen-bit rows finish, and above `i16::MAX`,
/// where they stop and the 32-bit rows give the reference's result; forward and reverse, with
/// gap_extend 1 and 256 (the `16 * gap_extend` headroom). The AVX2 and the 128-bit
/// (LOSAT_X_BPGAP128) sixteen-bit rows finish and stop on the same problems: their rows end
/// with the same best scores.
#[test]
#[cfg_attr(
    not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )),
    expect(
        unused_variables,
        unused_mut,
        reason = "the comparison with the sixteen-bit rows below exists only on x86_64 and wasm32 simd128"
    )
)]
fn x_bpgap16_falls_back_past_the_16_bit_bound() {
    let mut scratch = GapAlignScratch::new();
    for (go, ge, x, lens) in [
        (11, 1, 38, [2_900usize, 2_990, 3_100, 4_000]),
        (3840, 256, 4096, [2_560, 2_640, 2_800, 3_100]),
    ] {
        let limit = i16::MAX as i32 - 11 - 16 * ge - 1;
        let (mut finished, mut stopped) = (0, 0);
        for len in lens {
            let mut q = vec![20u8; len];
            let mut s = vec![20u8; len];
            // A few differences so the band and the gap rows do some work.
            for k in (37..len).step_by(301) {
                q[k] = 1;
                if ge == 1 {
                    s.insert(k / 2, 5);
                }
            }
            let (ql, sl) = (q.len(), s.len());
            for rev in [false, true] {
                // As `check_seed`: reverse from a seed at the last letters, forward from one
                // at the first.
                let (l1, l2) = if rev { (ql, sl) } else { (ql - 1, sl - 1) };
                let want = check(
                    &q,
                    0,
                    &s,
                    0,
                    l1,
                    l2,
                    MATRICES[0],
                    go,
                    ge,
                    x,
                    rev,
                    &mut scratch,
                );
                #[cfg(any(
                    target_arch = "x86_64",
                    all(target_arch = "wasm32", target_feature = "simd128")
                ))]
                {
                    let mut first: Option<Option<(usize, usize, i32)>> = None;
                    for rows in narrow_kernels() {
                        let got = narrow16(
                            rows,
                            &q,
                            0,
                            &s,
                            0,
                            l1,
                            l2,
                            MATRICES[0],
                            go,
                            ge,
                            x,
                            rev,
                            &mut scratch,
                        )
                        .expect("fits16");
                        let case =
                            format!("{rows:?} go={go} ge={ge} len={len} rev={rev} want={want:?}");
                        if want.2 <= limit {
                            assert_eq!(got, Some(want), "{case}");
                            finished += 1;
                        } else if want.2 > i16::MAX as i32 {
                            assert_eq!(got, None, "{case}");
                            stopped += 1;
                        } else {
                            assert!(got.is_none() || got == Some(want), "{case}");
                        }
                        assert_eq!(*first.get_or_insert(got), got, "{case}");
                    }
                }
            }
        }
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        if !narrow_kernels().is_empty() {
            assert!(
                finished > 0 && stopped > 0,
                "go={go} ge={ge}: {finished} {stopped}"
            );
        }
    }
}

/// LOSAT_X_BPGAP16 at the edges of `fits16`: X-drop and `gap_open_extend` of 4,096 and
/// `gap_extend` of 256 (the sixteen-bit rows run), one more (they do not), with gap and stop
/// letters, whose entries (down to -32,768) saturate the 16-bit sums; for the AVX2 and the
/// 128-bit sixteen-bit rows.
#[test]
#[cfg_attr(
    not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )),
    expect(
        unused_variables,
        unused_mut,
        reason = "the comparison with the sixteen-bit rows below exists only on x86_64 and wasm32 simd128"
    )
)]
fn x_bpgap16_matches_at_the_edges_of_its_costs() {
    let mut rng = Rng(0x0BAD_5EED_1616_0001);
    let mut scratch = GapAlignScratch::new();
    let (mut taken, mut declined) = (0, 0);
    for round in 0..240 {
        let len = 20 + rng.below(300) as usize;
        let q = random_seq(&mut rng, len, 20, round % 3 == 0);
        let s = mutate(&mut rng, &q, [5, 30][round % 2]);
        if s.is_empty() {
            continue;
        }
        let m = MATRICES[round % 4];
        let (go, ge) = [(3840, 256), (3841, 256), (4095, 1), (0, 257), (11, 1)][round % 5];
        let x = [4096, 4097, 38][rng.below(3) as usize];
        let fits = x.max(go + ge) <= 4096 && go + ge <= 4096 && ge <= 256;
        let qs = rng.below(q.len() as u64) as usize;
        let ss = qs.min(s.len() - 1);
        let (ql, sl) = (qs + 1, ss + 1);
        let want = check(&q, 0, &s, 0, ql, sl, m, go, ge, x, true, &mut scratch);
        check_seed(&q, &s, qs, ss, m, go, ge, x, &mut scratch);
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        for rows in narrow_kernels() {
            let got = narrow16(rows, &q, 0, &s, 0, ql, sl, m, go, ge, x, true, &mut scratch);
            let case = format!("{rows:?} round {round} go={go} ge={ge} x={x}");
            if fits {
                assert_eq!(got, Some(Some(want)), "{case}");
                taken += 1;
            } else {
                assert_eq!(got, None, "{case}");
                declined += 1;
            }
        }
    }
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if !narrow_kernels().is_empty() {
        assert!(taken > 0 && declined > 0, "{taken} {declined}");
    }
}

/// LOSAT_X_BPGAP128: the 128-bit rows against the reference and the AVX2 sixteen-bit rows on
/// the shapes their steps meet: rows of 1 to 40 cells (every count of cells in a row's last
/// pair of steps, 1 to 16), bands that start at every column phase of `RESTRICT_SIZE`, X-drops
/// from 0 (rows whose cells are all dropped) to 4,096, gap costs from 0 up, ties of two-letter
/// alphabets, and long related sequences whose bands move and widen. The results are the
/// reference's and those of the AVX2 rows (where the CPU has them), and the 128-bit rows take
/// every problem below the 16-bit bound.
#[test]
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_bpgap128_rows_match_the_reference_and_the_avx2_rows() {
    if !cpu_has_w128() {
        return;
    }
    let mut rng = Rng(0x1281_2812_8128_1281);
    let mut scratch = GapAlignScratch::new();
    let gaps = [
        (11, 1),
        (10, 1),
        (8, 2),
        (5, 2),
        (0, 1),
        (1, 0),
        (0, 0),
        (15, 3),
        (100, 30),
    ];
    let mut taken = 0;
    for round in 0..6000 {
        let (q, s) = match round % 4 {
            0 => {
                let alphabet = [2u64, 4, 20][rng.below(3) as usize];
                let ql = 1 + rng.below(40) as usize;
                let sl = 1 + rng.below(40) as usize;
                (
                    random_seq(&mut rng, ql, alphabet, round % 3 == 0),
                    random_seq(&mut rng, sl, alphabet, round % 5 == 0),
                )
            }
            1 => {
                let len = 1 + rng.below(60) as usize;
                let rate = [0, 5, 30][rng.below(3) as usize];
                let q = random_seq(&mut rng, len, 20, round % 3 == 0);
                let s = mutate(&mut rng, &q, rate);
                (q, s)
            }
            2 => {
                let len = 20 + rng.below(if round % 40 == 2 { 2000 } else { 300 }) as usize;
                let rate = [5, 15, 30, 50][rng.below(4) as usize];
                let q = random_seq(&mut rng, len, 20, round % 3 == 0);
                let s = mutate(&mut rng, &q, rate);
                (q, s)
            }
            _ => {
                let ql = 1 + rng.below(200) as usize;
                let sl = 1 + rng.below(200) as usize;
                (
                    random_seq(&mut rng, ql, 2, false),
                    random_seq(&mut rng, sl, 2, false),
                )
            }
        };
        if s.is_empty() {
            continue;
        }
        let m = MATRICES[rng.below(4) as usize];
        let (go, ge) = gaps[rng.below(gaps.len() as u64) as usize];
        let x = [0, 1, 5, 11, 16, 38, 100, 300, 1000, 4096][rng.below(10) as usize];
        let qs = rng.below(q.len() as u64) as usize;
        let ss = rng.below(s.len() as u64) as usize;
        for (qb, sb, len1, len2, rev) in [
            (0, 0, qs + 1, ss + 1, true),
            (qs, ss, q.len() - qs - 1, s.len() - ss - 1, false),
        ] {
            if len1 == 0 || len2 == 0 {
                continue;
            }
            let want = reference(&q, qb, &s, sb, len1, len2, m, go, ge, x, rev, &mut scratch);
            let case =
                format!("round={round} len1={len1} len2={len2} go={go} ge={ge} x={x} rev={rev}");
            // Below the 16-bit bound every sixteen-bit row finishes.
            let may_stop = w128_may_decline(m, go, ge, x, want.2);
            let mut first: Option<Option<(usize, usize, i32)>> = None;
            for rows in narrow_kernels() {
                let got = narrow16(
                    rows,
                    &q,
                    qb,
                    &s,
                    sb,
                    len1,
                    len2,
                    m,
                    go,
                    ge,
                    x,
                    rev,
                    &mut scratch,
                )
                .expect("fits16");
                assert!(
                    got == Some(want) || (may_stop && got.is_none()),
                    "{rows:?} {case}: {got:?} {want:?}"
                );
                assert_eq!(*first.get_or_insert(got), got, "{rows:?} {case}");
            }
            for band in [false, true] {
                let mut own = XRestrictedScratch::new();
                let got = restricted_fast(
                    &q,
                    qb,
                    &s,
                    sb,
                    len1,
                    len2,
                    m,
                    go,
                    ge,
                    x,
                    rev,
                    &mut own,
                    XRows::W128,
                    true,
                    band,
                );
                assert_rows_result(XRows::W128, got, want, m, (go, ge, x), &case);
                taken += usize::from(got.is_some());
            }
        }
    }
    assert!(taken > 10_000, "{taken}");
}

/// The `arch` operations of LOSAT_X_BPGAP128 (`utils/xdrop_simd.rs`: `subs16`, `up16_1`,
/// `up16_2`, `up16_4`, `lut32k`) against their definitions, on random and extreme lanes and
/// every index below 32.
#[test]
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_bpgap128_layer_operations_match_their_definitions() {
    use crate::utils::xdrop_simd::arch::*;
    if !cpu_has_w128() {
        return;
    }
    let mut rng = Rng(0x5EED_0128_0128_5EED);
    let lane = |rng: &mut Rng| -> i16 {
        match rng.below(4) {
            0 => [i16::MIN, i16::MIN + 1, -1, 0, 1, i16::MAX - 1, i16::MAX][rng.below(7) as usize],
            _ => rng.next() as i16,
        }
    };
    for _ in 0..5_000 {
        let a: [i16; 8] = std::array::from_fn(|_| lane(&mut rng));
        let b: [i16; 8] = std::array::from_fn(|_| lane(&mut rng));
        let t0: [u8; 16] = std::array::from_fn(|_| rng.next() as u8);
        let t1: [u8; 16] = std::array::from_fn(|_| rng.next() as u8);
        let idx: [u8; 16] = std::array::from_fn(|_| rng.below(32) as u8);
        let mut w = [[0i16; 8]; 4];
        let mut bytes = [0u8; 16];
        // SAFETY: SSSE3 and SSE4.1 on x86_64 (checked above), simd128 on wasm32; the loads and
        // stores are of 16-byte arrays.
        unsafe {
            let (va, vb) = (load(a.as_ptr()), load(b.as_ptr()));
            store(w[0].as_mut_ptr(), subs16(va, vb));
            store(w[1].as_mut_ptr(), up16_1(va, vb));
            store(w[2].as_mut_ptr(), up16_2(va, vb));
            store(w[3].as_mut_ptr(), up16_4(va, vb));
            let tables = (load(t0.as_ptr()), load(t1.as_ptr()));
            let looked = lut32k(tables.0, tables.1, load(idx.as_ptr()), splat8(-16));
            store(bytes.as_mut_ptr(), looked);
        }
        let sat: [i16; 8] = std::array::from_fn(|i| a[i].saturating_sub(b[i]));
        assert_eq!(w[0], sat, "subs16 {a:?} {b:?}");
        for (k, got) in [(1, w[1]), (2, w[2]), (4, w[3])] {
            let want: [i16; 8] =
                std::array::from_fn(|i| if i < k { b[8 - k + i] } else { a[i - k] });
            assert_eq!(got, want, "up16_{k} {a:?} {b:?}");
        }
        let want: [u8; 16] = std::array::from_fn(|i| {
            let j = idx[i] as usize;
            if j < 16 {
                t0[j]
            } else {
                t1[j - 16]
            }
        });
        assert_eq!(bytes, want, "lut32k {idx:?}");
    }
}

#[test]
fn x_bpgap_matches_reference_on_random_problems() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut scratch = GapAlignScratch::new();
    let gaps = [
        (11, 1),
        (10, 1),
        (9, 2),
        (5, 2),
        (0, 1),
        (1, 0),
        (0, 0),
        (15, 3),
    ];
    for round in 0..3000 {
        let alphabet = [2u64, 4, 20, 24][rng.below(4) as usize];
        let ql = 1 + rng.below(200) as usize;
        let sl = 1 + rng.below(200) as usize;
        let q = random_seq(&mut rng, ql, alphabet, round % 3 == 0);
        let s = random_seq(&mut rng, sl, alphabet, round % 5 == 0);
        let m = MATRICES[rng.below(4) as usize];
        let (go, ge) = gaps[rng.below(gaps.len() as u64) as usize];
        let x = [0, 5, 12, 16, 38, 60, 150][rng.below(7) as usize];
        let qs = rng.below(ql as u64) as usize;
        let ss = rng.below(sl as u64) as usize;
        check_seed(&q, &s, qs, ss, m, go, ge, x, &mut scratch);
    }
}

#[test]
fn x_bpgap_matches_reference_on_related_sequences() {
    // Long alignments: wide bands, many rows where a gap opens, scores past 255.
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let mut scratch = GapAlignScratch::new();
    for round in 0..400 {
        let len = 20 + rng.below(if round % 10 == 0 { 3000 } else { 400 }) as usize;
        let q = random_seq(&mut rng, len, 20, round % 4 == 0);
        let s = mutate(&mut rng, &q, [5, 15, 30, 50][round % 4]);
        if s.is_empty() {
            continue;
        }
        let m = MATRICES[round % 4];
        let (go, ge) = [(11, 1), (10, 1), (8, 2)][round % 3];
        let x = [38, 16, 100, 300][rng.below(4) as usize];
        let qs = rng.below(q.len() as u64) as usize;
        let ss = qs.min(s.len() - 1);
        check_seed(&q, &s, qs, ss, m, go, ge, x, &mut scratch);
    }
}

#[test]
fn x_bpgap_matches_reference_on_ties_and_drop_boundaries() {
    // Two letters with scores that tie often; every X-drop from 0 to 40 meets cells exactly
    // at `best_score - x_dropoff` and best scores reached again later in a row.
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    let mut scratch = GapAlignScratch::new();
    for x in 0..=40 {
        for _ in 0..40 {
            let ql = 1 + rng.below(80) as usize;
            let sl = 1 + rng.below(80) as usize;
            let q = random_seq(&mut rng, ql, 2, false);
            let s = random_seq(&mut rng, sl, 2, false);
            let qs = rng.below(ql as u64) as usize;
            let ss = rng.below(sl as u64) as usize;
            for (go, ge) in [(11, 1), (3, 1), (2, 2)] {
                check_seed(&q, &s, qs, ss, MATRICES[0], go, ge, x, &mut scratch);
            }
        }
    }
    // Repeats: the same score on many diagonals.
    let unit = [1u8, 7, 9, 3, 1, 7];
    let q: Vec<u8> = unit.iter().cycle().take(300).copied().collect();
    let s: Vec<u8> = unit.iter().cycle().skip(2).take(280).copied().collect();
    for x in [0, 1, 7, 12, 38, 90] {
        check_seed(&q, &s, 150, 140, MATRICES[0], 11, 1, x, &mut scratch);
        check_seed(&q, &s, 0, 0, MATRICES[0], 11, 1, x, &mut scratch);
        check_seed(&q, &s, 299, 279, MATRICES[0], 11, 1, x, &mut scratch);
    }
}

#[test]
fn x_bpgap_matches_reference_at_the_ends() {
    // Lengths 1-3, seeds at both ends (sentinels read past the window), subject windows that
    // start inside the subject (the reverse call's `subject_base`), and the gap letter.
    let mut rng = Rng(0x0F0F_1E1E_2D2D_3C3C);
    let mut scratch = GapAlignScratch::new();
    for ql in 1..=4usize {
        for sl in 1..=4usize {
            for _ in 0..30 {
                let q = random_seq(&mut rng, ql, 20, true);
                let s = random_seq(&mut rng, sl, 20, true);
                for qs in 0..ql {
                    for ss in 0..sl {
                        check_seed(&q, &s, qs, ss, MATRICES[0], 11, 1, 38, &mut scratch);
                    }
                }
            }
        }
    }
    for _ in 0..300 {
        let q = random_seq(&mut rng, 60, 20, true);
        let s = random_seq(&mut rng, 90, 20, true);
        let sb = rng.below(30) as usize;
        let len2 = 1 + rng.below((90 - sb - 1) as u64) as usize;
        let len1 = 1 + rng.below(59) as usize;
        for rev in [false, true] {
            let qb = if rev {
                0
            } else {
                rng.below((60 - len1) as u64) as usize
            };
            let qb = qb.min(59 - len1.min(59));
            check(
                &q,
                qb,
                &s,
                sb,
                len1,
                len2,
                MATRICES[0],
                11,
                1,
                38,
                rev,
                &mut scratch,
            );
        }
    }
    // Empty problems.
    let q = [1u8, 2, 3];
    assert_eq!(
        check(
            &q,
            0,
            &q,
            0,
            0,
            2,
            MATRICES[0],
            11,
            1,
            38,
            false,
            &mut scratch
        ),
        (0, 0, 0)
    );
    assert_eq!(
        check(
            &q,
            0,
            &q,
            0,
            2,
            0,
            MATRICES[0],
            11,
            1,
            38,
            true,
            &mut scratch
        ),
        (0, 0, 0)
    );
}

#[test]
fn x_bpgap_reused_scratch_matches_fresh_scratch() {
    // One scratch over problems of very different sizes (stale cells from a longer problem
    // must not leak into a shorter one) against a fresh scratch per problem.
    let mut rng = Rng(0xAAAA_5555_3333_CCCC);
    let mut shared = GapAlignScratch::new();
    for round in 0..300 {
        let len = [5usize, 2000, 30, 700, 1][round % 5] + rng.below(20) as usize;
        let q = random_seq(&mut rng, len, 20, false);
        let s = mutate(&mut rng, &q, 10);
        if s.is_empty() {
            continue;
        }
        let qs = rng.below(q.len() as u64) as usize;
        let ss = qs.min(s.len() - 1);
        let x = [16, 38, 200][round % 3];
        check_seed(&q, &s, qs, ss, MATRICES[0], 11, 1, x, &mut shared);
        let mut fresh = GapAlignScratch::new();
        check_seed(&q, &s, qs, ss, MATRICES[0], 11, 1, x, &mut fresh);
    }
}

#[test]
fn x_bpgapband_reused_scratch_matches_reference() {
    // LOSAT_X_BPGAPBAND: one scratch whose arrays long bands grow, then reused by short bands,
    // by long subjects that an unrelated query leaves after a few cells, and by other costs.
    let mut rng = Rng(0x0BAD_5EED_1234_4321);
    let mut scratch = GapAlignScratch::new();
    let mut shared = XRestrictedScratch::new();
    let kernels = kernels();
    for round in 0..300 {
        let len = [5usize, 3000, 30, 700, 1][round % 5] + rng.below(20) as usize;
        let q = random_seq(&mut rng, len, 20, false);
        let s = if round % 2 == 0 {
            mutate(&mut rng, &q, 10)
        } else {
            random_seq(&mut rng, 3 * len + 50, 20, false)
        };
        if s.is_empty() {
            continue;
        }
        let qs = rng.below(q.len() as u64) as usize;
        let ss = qs.min(s.len() - 1);
        let (go, ge) = [(11, 1), (10, 1), (8, 2)][round % 3];
        let x = [16, 38, 200][rng.below(3) as usize];
        for (qb, sb, len1, len2, rev) in [
            (0, 0, qs + 1, ss + 1, true),
            (qs, ss, q.len() - qs - 1, s.len() - ss - 1, false),
        ] {
            if len1 == 0 || len2 == 0 {
                continue;
            }
            let m = MATRICES[0];
            let want = reference(&q, qb, &s, sb, len1, len2, m, go, ge, x, rev, &mut scratch);
            for &(rows, narrow) in &kernels {
                let got = restricted_fast(
                    &q,
                    qb,
                    &s,
                    sb,
                    len1,
                    len2,
                    m,
                    go,
                    ge,
                    x,
                    rev,
                    &mut shared,
                    rows,
                    narrow,
                    true,
                );
                let case = format!(
                    "round={round} rows={rows:?} narrow={narrow} len1={len1} len2={len2} rev={rev}"
                );
                assert_rows_result(rows, got, want, m, (go, ge, x), &case);
            }
        }
    }
}

#[test]
fn x_bpgapband_arrays_hold_the_band_not_the_subject() {
    // LOSAT_X_BPGAPBAND: an unrelated query against a subject of 50,000 letters leaves the
    // band after a few cells, so the arrays keep the first rows' cells (without the switch
    // they hold the subject's letters); a related pair of 3,000 letters grows them to its
    // band. Both give the reference's result, with every kernel of this machine (the 128-bit
    // rows of LOSAT_X_BPGAP128 use the sixteen-bit arrays only).
    let mut rng = Rng(0x5E_ED0F_BA4D_0001);
    let mut scratch = GapAlignScratch::new();
    let m = MATRICES[0];
    let q = random_seq(&mut rng, 300, 20, false);
    let s = random_seq(&mut rng, 50_000, 20, false);
    let want = reference(&q, 0, &s, 0, 300, 50_000, m, 11, 1, 38, false, &mut scratch);
    let q2 = random_seq(&mut rng, 3_000, 20, false);
    let s2 = mutate(&mut rng, &q2, 5);
    let len2 = s2.len() - 1;
    let want2 = reference(
        &q2,
        0,
        &s2,
        0,
        2_999,
        len2,
        m,
        11,
        1,
        38,
        false,
        &mut scratch,
    );
    for (rows, narrow) in kernels() {
        let held = |sc: &XRestrictedScratch| {
            #[cfg(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            ))]
            if rows == XRows::W128 {
                return sc.h16[0].len();
            }
            sc.h[0].len()
        };
        for band in [false, true] {
            let mut sc = XRestrictedScratch::new();
            let got = restricted_fast(
                &q, 0, &s, 0, 300, 50_000, m, 11, 1, 38, false, &mut sc, rows, narrow, band,
            );
            assert_eq!(got, Some(want), "{rows:?} band={band}");
            let cells = held(&sc);
            if band {
                assert!(
                    cells < 1_000,
                    "{rows:?}: {cells} cells for a band of a few cells"
                );
            } else {
                assert_eq!(cells, PAD + 50_000 + 2 + TAIL, "{rows:?}");
            }
        }
        let mut sc = XRestrictedScratch::new();
        let got = restricted_fast(
            &q2, 0, &s2, 0, 2_999, len2, m, 11, 1, 38, false, &mut sc, rows, narrow, true,
        );
        assert_eq!(got, Some(want2), "{rows:?}");
        assert!(held(&sc) > 2_000 || sc.h16[0].len() > 2_000, "{rows:?}");
    }
}

#[test]
fn x_bpgap_column_letters_are_the_reference_reads() {
    // Every window position, including bases and lengths that reach past the subject.
    let mut rng = Rng(0x7777_1111_2222_9999);
    // Room for `2 + TAIL` past the letters (`band_col`).
    let mut col = vec![0u8; 200 + 2 + TAIL];
    for _ in 0..4000 {
        let sl = rng.below(40) as usize;
        let subject = random_seq(&mut rng, sl, 27, false);
        let subject_base = rng.below(50) as usize;
        let len2 = rng.below(60) as usize;
        let from = rng.below(len2 as u64 + 30) as usize;
        let to = from + rng.below(100) as usize;
        for reverse in [false, true] {
            let ok = fill_col::<COL32>(&mut col, from, to, &subject, subject_base, len2, reverse);
            for (c, &v) in col.iter().enumerate().take(to).skip(from) {
                let want = if c == 0 || c > len2 + 1 {
                    0
                } else if reverse {
                    blastp_subject_residue::<true>(&subject, subject_base, len2, c - 1)
                } else {
                    blastp_subject_residue::<false>(&subject, subject_base, len2, c - 1)
                };
                assert_eq!(
                    v, want,
                    "c={c} base={subject_base} len2={len2} rev={reverse}"
                );
            }
            assert_eq!(ok, Ok(()));
        }
    }
}

#[test]
fn x_bpgap_leaves_unusual_problems_to_the_reference() {
    let mut scratch = GapAlignScratch::new();
    let q = [1u8, 2, 3, 4];
    let s = [1u8, 2, 30, 4];
    let kernels = kernels();
    // A letter code of 28 or more, in the subject and in the query.
    let q30 = [1u8, 30, 3, 4];
    for &(rows, narrow) in &kernels {
        for (a, b) in [(&q, &s), (&q30, &q)] {
            assert!(restricted_fast(
                a,
                0,
                b,
                0,
                3,
                3,
                MATRICES[0],
                11,
                1,
                38,
                false,
                &mut scratch.x_restricted,
                rows,
                narrow,
                x_bpgapband()
            )
            .is_none());
        }
    }
    // Gap costs and X-drop outside the checked ranges.
    for &(rows, narrow) in &kernels {
        for (go, ge, x) in [
            (-1, 1, 38),
            (11, -1, 38),
            (MAX_GAP_COST + 1, 1, 38),
            (11, 1, MAX_X_DROP + 1),
        ] {
            assert!(restricted_fast(
                &q,
                0,
                &q,
                0,
                3,
                3,
                MATRICES[0],
                go,
                ge,
                x,
                false,
                &mut scratch.x_restricted,
                rows,
                narrow,
                x_bpgapband()
            )
            .is_none());
        }
    }
}

/// Timing aid (ignored; run with `cargo test --release -- --ignored x_bpgap_time_kernels
/// --nocapture`): the reference and every kernel of this machine (`kernels`: the scalar rows,
/// the AVX2 rows, the 128-bit rows) on the same unrelated sequence pairs (the restricted
/// alignments of the preliminary stage are low-scoring ones), one after the other, three
/// times. Prints ns per call.
#[test]
#[ignore]
fn x_bpgap_time_kernels() {
    let mut rng = Rng(0x5151_7272_9393_B4B4);
    let mut problems = Vec::new();
    for _ in 0..4000 {
        let (ql, sl) = (100 + rng.below(500) as usize, 100 + rng.below(500) as usize);
        let q = random_seq(&mut rng, ql, 20, false);
        let s = random_seq(&mut rng, sl, 20, false);
        let qs = rng.below(q.len() as u64) as usize;
        let ss = rng.below(s.len() as u64) as usize;
        problems.push((q, s, qs, ss));
    }
    let mut scratch = GapAlignScratch::new();
    let m = MATRICES[0];
    let calls = |scratch: &mut GapAlignScratch, kernel: Option<(XRows, bool)>| -> (u64, i64) {
        let mut sum = 0i64;
        let mut n = 0u64;
        for (q, s, qs, ss) in &problems {
            let (qs, ss) = (*qs, *ss);
            for (qb, sb, len1, len2, rev) in [
                (0, 0, qs + 1, ss + 1, true),
                (qs, ss, q.len() - qs - 1, s.len() - ss - 1, false),
            ] {
                let r = match kernel {
                    None => reference(q, qb, s, sb, len1, len2, m, 11, 1, 38, rev, scratch),
                    Some((rows, narrow)) => restricted_fast(
                        q,
                        qb,
                        s,
                        sb,
                        len1,
                        len2,
                        m,
                        11,
                        1,
                        38,
                        rev,
                        &mut scratch.x_restricted,
                        rows,
                        narrow,
                        x_bpgapband(),
                    )
                    .unwrap(),
                };
                sum += r.2 as i64 + r.0 as i64;
                n += 1;
            }
        }
        (n, sum)
    };
    let kernels: Vec<Option<(XRows, bool)>> = std::iter::once(None)
        .chain(kernels().into_iter().map(Some))
        .collect();
    for round in 0..3 {
        for &k in &kernels {
            let t = std::time::Instant::now();
            let (n, sum) = calls(&mut scratch, k);
            let ns = t.elapsed().as_nanos() as f64;
            let name = match k {
                None => "reference".to_string(),
                Some((rows, narrow)) => format!("{rows:?} narrow={narrow}"),
            };
            eprintln!(
                "round {round} kernel {name} calls {n} {:.1} ns/call checksum {sum}",
                ns / n as f64
            );
        }
    }
}
