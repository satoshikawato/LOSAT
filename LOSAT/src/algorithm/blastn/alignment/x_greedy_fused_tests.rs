use super::super::{blast_greedy_align, GapPrelimEditBlock};
use super::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn pick<T: Copy>(&mut self, values: &[T]) -> T {
        values[self.below(values.len() as u64) as usize]
    }
}

/// One call of `BLAST_GreedyAlign` without an edit block on unpacked sequences.
struct Problem {
    seq1: Vec<u8>,
    len1: i32,
    seq2: Vec<u8>,
    len2: i32,
    reverse: bool,
    xdrop: i32,
    match_cost: i32,
    mismatch_cost: i32,
    max_dist: i32,
    /// `fence_hit` on entry (callers pass false; true must leave the ported loop running).
    fence_in: bool,
}

#[derive(Debug, PartialEq, Eq)]
struct Outcome {
    ret: i32,
    seq1_align_len: i32,
    seq2_align_len: i32,
    seed: (i32, i32, i32),
    fence_hit: bool,
}

/// `blast_greedy_align` score only (this module, in the plain build when `plain`), or with an
/// edit block (`traceback`: `x_greedy_row`, the loop this module replaces without one).
fn call(traceback: bool, plain: bool, mem: &mut GreedyNonAffineMem, p: &Problem) -> Outcome {
    X_GREEDY_FUSE_TEST_PLAIN.with(|m| m.set(plain));
    let mut edit = GapPrelimEditBlock::new();
    // values the call overwrites
    let mut seed = GreedySeed {
        start_q: -7,
        start_s: -7,
        match_length: -7,
    };
    let mut fence_hit = p.fence_in;
    let (mut seq1_align_len, mut seq2_align_len) = (-5, -5);
    let ret = blast_greedy_align(
        &p.seq1,
        p.len1,
        &p.seq2,
        p.len2,
        p.reverse,
        p.xdrop,
        p.match_cost,
        p.mismatch_cost,
        &mut seq1_align_len,
        &mut seq2_align_len,
        p.max_dist,
        mem,
        traceback.then_some(&mut edit),
        &mut fence_hit,
        &mut seed,
    );
    X_GREEDY_FUSE_TEST_PLAIN.with(|m| m.set(false));
    Outcome {
        ret,
        seq1_align_len,
        seq2_align_len,
        seed: (seed.start_q, seed.start_s, seed.match_length),
        fence_hit,
    }
}

/// `max_score`, whole (the rows of a traceback call are not the rows of a score-only call).
fn assert_same_max_score(fused: &GreedyNonAffineMem, reference: &GreedyNonAffineMem, what: &str) {
    assert!(
        fused.max_score == reference.max_score,
        "{what}: max_score differs"
    );
}

/// The memory two score-only calls leave: both rows and `max_score`, whole.
fn assert_same_memory(a: &GreedyNonAffineMem, b: &GreedyNonAffineMem, what: &str) {
    for row in 0..2 {
        assert!(
            a.last_seq2_off[row] == b.last_seq2_off[row],
            "{what}: row {row} differs"
        );
    }
    assert_same_max_score(a, b, what);
}

/// Related sequences (60-100 % identity, substitutions and indels), ambiguity bytes in one or
/// both, sometimes a fence byte in `seq2`, unrelated tails, bytes after `len1` / `len2` that
/// must not be read; lengths 0-20 often, up to about 4,000.
fn problem(rng: &mut Rng) -> Problem {
    let template_len = match rng.below(8) {
        0..=2 => rng.below(21),
        3 | 4 => 21 + rng.below(100),
        5 | 6 => 121 + rng.below(800),
        _ => 1000 + rng.below(3000),
    } as usize;
    let identity = rng.pick(&[
        1000u64, 1000, 995, 980, 950, 900, 850, 800, 750, 700, 650, 600,
    ]);
    let indel_share = rng.pick(&[0u64, 0, 100, 300, 600]); // of the differences, per mille
    let ambiguity1 = rng.pick(&[0u64, 0, 0, 5, 30]); // per mille
    let ambiguity2 = rng.pick(&[0u64, 0, 0, 5, 30]);
    let mut seq1 = Vec::new();
    let mut seq2 = Vec::new();
    for _ in 0..template_len {
        let base = rng.below(4) as u8;
        if rng.below(1000) < identity {
            seq1.push(base);
            seq2.push(base);
            continue;
        }
        let roll = rng.below(1000);
        if roll < indel_share / 2 {
            seq1.push(base);
        } else if roll < indel_share {
            seq2.push(base);
        } else {
            seq1.push(base);
            seq2.push((base + 1 + rng.below(3) as u8) & 3);
        }
    }
    if rng.below(3) == 0 {
        for _ in 0..rng.below(200) {
            seq1.push(rng.below(4) as u8);
        }
        for _ in 0..rng.below(200) {
            seq2.push(rng.below(4) as u8);
        }
    }
    for byte in seq1.iter_mut() {
        if rng.below(1000) < ambiguity1 {
            *byte = 4 + rng.below(12) as u8;
        }
    }
    for byte in seq2.iter_mut() {
        if rng.below(1000) < ambiguity2 {
            *byte = 4 + rng.below(12) as u8;
        }
    }
    if !seq2.is_empty() && rng.below(5) == 0 {
        let at = rng.below(seq2.len() as u64) as usize;
        seq2[at] = FENCE_SENTRY;
    }
    let len1 = seq1.len() as i32;
    let len2 = seq2.len() as i32;
    let reverse = rng.below(2) == 1;
    if reverse {
        // the left extension reads seq[len - 1 - i]
        seq1.reverse();
        seq2.reverse();
    }
    // bytes after the lengths: matching bases and fence bytes, which a read past the end
    // would turn into a different result
    for _ in 0..rng.below(16) {
        seq1.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
    }
    for _ in 0..rng.below(16) {
        seq2.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
    }
    // megablast's costs (reward 1 doubled), and costs with match_cost / 2 != 1
    let (match_cost, mismatch_cost) = rng.pick(&[
        (2, 4),
        (2, 4),
        (2, 4),
        (2, 6),
        (2, 2),
        (4, 6),
        (6, 8),
        (4, 4),
    ]);
    let xdrop = if template_len <= 500 {
        rng.pick(&[
            0, 1, 2, 5, 10, 13, 20, 26, 27, 40, 54, 100, 300, 1000, 100_000,
        ])
    } else {
        rng.pick(&[0, 1, 5, 13, 20, 26, 27, 40, 54, 100])
    };
    // GREEDY_MAX_COST and GREEDY_MAX_COST_FRACTION as the caller uses them, and small values
    // that stop the loop at max_dist
    let initial = 1000.min(len2 / 2 + 1);
    let max_dist = match rng.below(6) {
        0 => 1 + rng.below(4) as i32,
        1 => 1 + rng.below(12) as i32,
        2 | 3 => initial,
        4 => 2 * initial,
        _ => initial.min(30) + rng.below(40) as i32,
    };
    Problem {
        seq1,
        len1,
        seq2,
        len2,
        reverse,
        xdrop,
        match_cost,
        mismatch_cost,
        max_dist,
        fence_in: rng.below(25) == 0,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509-537,682-686
/// Sequences of calls on one scratch per path, so that the cells earlier calls left matter:
/// `x_greedy_row` (a call with an edit block runs it), this module in the AVX2 build and in
/// the plain build. After every call: the same return value, align lengths, seed and
/// `fence_hit`, and the same `max_score` as `x_greedy_row`; the two builds also leave the
/// same rows.
#[test]
fn x_fused_greedy_matches_reference_on_call_sequences() {
    let histories: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(2000);
    let mut rng = Rng(0x7F4A_7C15_9E37_79B9);
    let (mut converged, mut stopped, mut fences, mut early, mut gated) = (0, 0, 0, 0, 0);
    let mut calls = 0usize;
    let mut unit_costs = 0usize;
    let mut other_costs = 0usize;
    X_GREEDY_FUSE_TEST_REST.with(|n| n.set((0, 0)));
    for history in 0..histories {
        let mut reference = GreedyNonAffineMem::new();
        let mut fused = GreedyNonAffineMem::new();
        let mut plain = GreedyNonAffineMem::new();
        for step in 0..1 + rng.below(12) {
            let p = problem(&mut rng);
            let what = format!(
                "history {history} step {step}: len1={} len2={} reverse={} xdrop={} costs={}/{} \
                     max_dist={} fence_in={}",
                p.len1,
                p.len2,
                p.reverse,
                p.xdrop,
                p.match_cost,
                p.mismatch_cost,
                p.max_dist,
                p.fence_in
            );
            let expected = call(true, false, &mut reference, &p);
            let got = call(false, false, &mut fused, &p);
            assert_eq!(got, expected, "{what}: fused");
            assert_same_max_score(&fused, &reference, &format!("{what}: fused"));
            let got = call(false, true, &mut plain, &p);
            assert_eq!(got, expected, "{what}: fused, plain build");
            assert_same_max_score(&plain, &reference, &format!("{what}: fused, plain build"));
            assert_same_memory(
                &plain,
                &fused,
                &format!("{what}: plain build and AVX2 build"),
            );

            calls += 1;
            if p.fence_in {
                gated += 1;
            } else if expected.seq1_align_len == expected.seed.2
                && expected.ret == 0
                && (expected.seed.2 == p.len1 || expected.seed.2 == p.len2)
            {
                early += 1;
            } else if expected.fence_hit {
                fences += 1;
            } else if expected.ret < 0 {
                stopped += 1;
            } else {
                converged += 1;
            }
            if p.match_cost == 2 {
                unit_costs += 1;
            } else {
                other_costs += 1;
            }
        }
    }
    let (short, long) = X_GREEDY_FUSE_TEST_REST.with(Cell::get);
    let summary = format!(
        "calls {calls}: converged {converged}, past max_dist {stopped}, fence {fences}, early \
             return {early}, fence on entry {gated}; slides with < 8 positions left {short}, \
             slides of 8 or more {long}; match_cost 2: {unit_costs}, other: {other_costs}"
    );
    eprintln!("x_fused_greedy_matches_reference_on_call_sequences: {summary}");
    // the generator must reach every exit and every path
    assert!(converged * 3 > calls, "{summary}");
    assert!(stopped * 20 > calls, "{summary}");
    assert!(fences * 40 > calls, "{summary}");
    assert!(early > 0 && gated > 0, "{summary}");
    assert!(short > 1000 && long > 1000, "{summary}");
    assert!(other_costs * 5 > calls, "{summary}");
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:331-334,356-359,571-576
/// A fence inside the loop, reached by a slide that stops within eight positions (the
/// eight-byte test), by a slide of eight or more, and with fewer than eight positions left:
/// 0 and `fence_hit` at once, with the `max_score`, align lengths and seed `x_greedy_row`
/// leaves, and the same rows in the plain build as in the AVX2 build, on scratch with stale
/// cells.
#[test]
fn x_fused_greedy_fence_inside_the_loop() {
    // (length, substitution at, fence at, slides of the fence diagonal: (short, long))
    let cases = [
        (40usize, 5usize, 12usize, (0u64, 0u64)),
        (40, 5, 30, (0, 1)),
        (10, 3, 8, (1, 0)),
    ];
    for reverse in [false, true] {
        for &(len, substitution, fence, rest) in &cases {
            let mut seq1: Vec<u8> = (0..len).map(|i| (i * 7 % 4) as u8).collect();
            let mut seq2 = seq1.clone();
            seq2[substitution] = (seq2[substitution] + 1) & 3;
            seq2[fence] = FENCE_SENTRY;
            if reverse {
                seq1.reverse();
                seq2.reverse();
            }
            let p = Problem {
                seq1,
                len1: len as i32,
                seq2,
                len2: len as i32,
                reverse,
                xdrop: 26,
                match_cost: 2,
                mismatch_cost: 4,
                max_dist: 21,
                fence_in: false,
            };
            let mut reference = GreedyNonAffineMem::new();
            let mut fused = GreedyNonAffineMem::new();
            let mut plain = GreedyNonAffineMem::new();
            reference.ensure_base_capacity(48);
            fused.ensure_base_capacity(48);
            plain.ensure_base_capacity(48);
            for row in 0..2 {
                reference.last_seq2_off[row].fill(123);
                fused.last_seq2_off[row].fill(123);
                plain.last_seq2_off[row].fill(123);
            }
            let expected = call(true, false, &mut reference, &p);
            assert!(expected.fence_hit && expected.ret == 0, "{expected:?}");
            X_GREEDY_FUSE_TEST_REST.with(|n| n.set((0, 0)));
            let got = call(false, false, &mut fused, &p);
            let taken = X_GREEDY_FUSE_TEST_REST.with(Cell::get);
            assert_eq!(got, expected, "len {len} fence {fence} reverse {reverse}");
            assert_same_max_score(&fused, &reference, "fence");
            let got = call(false, true, &mut plain, &p);
            assert_eq!(
                got, expected,
                "len {len} fence {fence} reverse {reverse}: plain"
            );
            assert_same_memory(&plain, &fused, "fence: plain build and AVX2 build");
            // the slides of distance 1 before the fence: the shifted diagonals stop early
            assert!(
                taken.0 >= rest.0 && taken.1 == rest.1,
                "len {len} fence {fence} reverse {reverse}: slides {taken:?}"
            );
        }
    }
}

/// No NCBI counterpart: the fused loop handles a kept diagonal past the end of `seq1`
/// (`x_past_seq1_end`) in place, and this test is the independent check of that branch. On
/// random lengths, bases, directions and offsets (`seq1_index` in `len1 + 1..=len1 + 5`,
/// `seq2_index` inside `seq2`), with and without the fence sentry at the byte that
/// `x_past_seq1_end` reads: it is true exactly when the sentry is absent, and then the
/// inline update (`diag_upper = k`, the cell is `seq2_index`, `curr_extent` and `curr_diag`
/// change only when `seq1_index + seq2_index` is larger, nothing else) equals what
/// `x_cell_cold` leaves from the same inputs; with the sentry `x_cell_cold` sets the fence.
#[test]
fn past_seq1_end_branch_equals_x_cell_cold() {
    let mut rng = Rng(0x005E_ED0F_E1D5_u64);
    let (mut inline, mut fenced) = (0usize, 0usize);
    for _ in 0..4000 {
        let len1 = 8 + rng.below(41) as i32;
        let len2 = 8 + rng.below(41) as i32;
        let seq1: Vec<u8> = (0..len1).map(|_| rng.below(4) as u8).collect();
        let mut seq2: Vec<u8> = (0..len2).map(|_| rng.below(4) as u8).collect();
        let reverse = rng.below(2) == 0;
        let seq2_index = rng.below(len2 as u64) as i32;
        let seq1_index = len1 + 1 + rng.below(5) as i32;
        let sentry = rng.below(2) == 0;
        let at = if reverse {
            len2 - 1 - seq2_index
        } else {
            seq2_index
        } as usize;
        if sentry {
            seq2[at] = FENCE_SENTRY;
        }
        let k = rng.below(8) as i32;
        let hot = XHot {
            diag_lower: rng.below(19) as i32 - 9,
            diag_upper: rng.below(19) as i32 - 9,
            curr_extent: rng.below((len1 + len2 + 16) as u64) as i32 - 5,
            curr_diag: rng.below(8) as i32,
            longest_match_run: rng.below(6) as i32,
        };
        let cold_in = XColdState {
            end1_reached: rng.below(2) == 0,
            end2_reached: rng.below(2) == 0,
            seed_start_q: rng.below(9) as i32 - 1,
            seed_start_s: rng.below(9) as i32 - 1,
            fence_hit: false,
        };
        let hot_tuple = |h: &XHot| {
            (
                h.diag_lower,
                h.diag_upper,
                h.curr_extent,
                h.curr_diag,
                h.longest_match_run,
            )
        };
        let cold_tuple = |c: &XColdState| {
            (
                c.end1_reached,
                c.end2_reached,
                c.seed_start_q,
                c.seed_start_s,
                c.fence_hit,
            )
        };
        let past = unsafe {
            if reverse {
                x_past_seq1_end::<true>(seq2.as_ptr(), len1, len2, seq1_index, seq2_index)
            } else {
                x_past_seq1_end::<false>(seq2.as_ptr(), len1, len2, seq1_index, seq2_index)
            }
        };
        assert_eq!(past, !sentry, "len1 {len1} len2 {len2} reverse {reverse}");

        let before = [-9i32; 8];
        let mut current = before;
        let mut cold = XColdState { ..cold_in };
        let got = unsafe {
            let (s1, s2) = (seq1.as_ptr(), seq2.as_ptr());
            let cur = current.as_mut_ptr();
            if reverse {
                x_cell_cold::<true>(
                    s1, s2, len1, len2, cur, k, seq1_index, seq2_index, hot, &mut cold,
                )
            } else {
                x_cell_cold::<false>(
                    s1, s2, len1, len2, cur, k, seq1_index, seq2_index, hot, &mut cold,
                )
            }
        };
        if !past {
            fenced += 1;
            assert!(cold.fence_hit, "the sentry sets the fence");
            continue;
        }
        inline += 1;
        // what the inline branch does
        let mut expected = hot;
        expected.diag_upper = k;
        let mut expected_current = before;
        expected_current[k as usize] = seq2_index;
        if seq1_index + seq2_index > expected.curr_extent {
            expected.curr_extent = seq1_index + seq2_index;
            expected.curr_diag = k;
        }
        let what = format!(
            "len1 {len1} len2 {len2} reverse {reverse} seq1_index {seq1_index} \
                 seq2_index {seq2_index}"
        );
        assert_eq!(hot_tuple(&got), hot_tuple(&expected), "{what}");
        assert_eq!(cold_tuple(&cold), cold_tuple(&cold_in), "{what}");
        assert_eq!(current, expected_current, "{what}");
    }
    assert!(
        inline > 500 && fenced > 500,
        "inline {inline}, fenced {fenced}"
    );
}
