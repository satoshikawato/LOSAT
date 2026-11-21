use super::super::{
    blast_greedy_align, format_gap_prelim_edit_block_for_trace, x_greedy_row_any,
    x_greedy_row_avx2, GapPrelimEditBlock, GreedyNonAffineMem, GreedySeed, XGreedyRowKernel,
    X_GREEDY_ROW_TEST_KERNEL,
};
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
    fn range(&mut self, low: i32, high: i32) -> i32 {
        low + self.below((high - low + 1) as u64) as i32
    }
}

fn avx2() -> bool {
    std::is_x86_feature_detected!("avx2")
        && std::is_x86_feature_detected!("bmi1")
        && std::is_x86_feature_detected!("bmi2")
        && std::is_x86_feature_detected!("lzcnt")
}

/// A pair of related sequences of `len` template bases (identity, substitutions, indels,
/// ambiguity bytes in either, sometimes fence bytes in `seq2`), reversed for `reverse` (the
/// left extension reads `seq[len - 1 - i]`), with matching bases and fence bytes after the
/// lengths that a read past the end would turn into a different result.
fn related(rng: &mut Rng, len: usize, reverse: bool) -> (Vec<u8>, i32, Vec<u8>, i32) {
    let identity = rng.pick(&[1000u64, 1000, 995, 980, 950, 900, 800, 700, 600]);
    let indel_share = rng.pick(&[0u64, 0, 100, 300, 600]); // of the differences, per mille
    let ambiguity1 = rng.pick(&[0u64, 0, 0, 5, 30]); // per mille
    let ambiguity2 = rng.pick(&[0u64, 0, 0, 5, 30]);
    let (mut seq1, mut seq2) = (Vec::new(), Vec::new());
    for _ in 0..len {
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
    if !seq2.is_empty() {
        for _ in 0..rng.pick(&[0u64, 0, 0, 1, 3]) {
            let at = rng.below(seq2.len() as u64) as usize;
            seq2[at] = FENCE_SENTRY;
        }
    }
    let (len1, len2) = (seq1.len() as i32, seq2.len() as i32);
    if reverse {
        seq1.reverse();
        seq2.reverse();
    }
    for _ in 0..rng.below(16) {
        seq1.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
    }
    for _ in 0..rng.below(16) {
        seq2.push(rng.pick(&[0u8, 1, 2, 3, FENCE_SENTRY]));
    }
    (seq1, len1, seq2, len2)
}

/// One distance: the sequences, the row of distance `d - 1` and the arguments.
struct Row {
    reverse: bool,
    seq1: Vec<u8>,
    len1: i32,
    seq2: Vec<u8>,
    len2: i32,
    previous: Vec<i32>,
    /// stale cells the call may leave as they are
    current: Vec<i32>,
    tmp_diag_lower: i32,
    diag_origin: i32,
    xdrop_score: i32,
    st: XGreedyRowState,
}

/// Bands of 1-600 diagonals (one to three blocks of `X_GREEDY_BLOCK`), lengths below eight
/// and up to about 300, offsets of distance `d - 1` from rows that rise along the band and
/// from noise (`INVALID_OFFSET`, -1, past the ends), the X-drop bound from always passing to
/// dropping most diagonals, the fence hit on entry, a seed run from 0 to above eight.
fn row(rng: &mut Rng) -> Row {
    let reverse = rng.below(2) == 1;
    let len = match rng.below(6) {
        0 => rng.below(12),
        1 | 2 => 12 + rng.below(60),
        _ => 72 + rng.below(240),
    } as usize;
    let (seq1, len1, seq2, len2) = related(rng, len, reverse);
    let band_len = match rng.below(6) {
        0 => 1 + rng.below(8),
        1..=3 => 9 + rng.below(250),
        _ => 250 + rng.below(350),
    } as usize;
    // seq1_index = seq2_index + first + cell
    let first = rng.range(-(band_len as i32) - 10, 20);
    let tmp_diag_lower = rng.range(-30, 30);
    let diag_origin = tmp_diag_lower - first;
    let noise = rng.pick(&[0u64, 0, 20, 100, 400]); // per mille of the cells
    let base = rng.range(-3, len2.max(1));
    let mut previous: Vec<i32> = (0..band_len + 2)
        .map(|j| {
            if rng.below(1000) < noise {
                let any = rng.range(-2, len2 + 3);
                return rng.pick(&[INVALID_OFFSET, -1, len2, len2 + 1, any]);
            }
            // near the main diagonal of the band the rows reach furthest
            let o = first + j as i32 - 1;
            (base - o.abs() / rng.pick(&[1, 2, 4]) + rng.range(-2, 2)).max(INVALID_OFFSET)
        })
        .collect();
    // The reference loop never takes a diagonal with seq1_index < 0 (its slide would index
    // before the sequence): raise the offsets that would make one survive with it until none
    // does (raising an offset only raises the starts of its three diagonals).
    loop {
        let mut raised = false;
        for i in 0..band_len {
            let s2 = previous[i + 2]
                .max(previous[i + 1])
                .saturating_add(1)
                .max(previous[i]);
            let o = first + i as i32;
            if s2 >= 0 && s2 + o < 0 {
                previous[i + 1] = previous[i + 1].max(-o - 1);
                raised = true;
            }
        }
        if !raised {
            break;
        }
    }
    let xdrop_score = match rng.below(4) {
        0 => i32::MIN / 4,
        1 => rng.range(-20, 20),
        _ => rng.range(0, 2 * (len1 + len2) + 10),
    };
    let diag_lower = match rng.below(4) {
        0 | 1 => tmp_diag_lower,
        2 => tmp_diag_lower + rng.range(0, 3),
        _ => rng.range(tmp_diag_lower - 5, tmp_diag_lower + band_len as i32 + 5),
    };
    let st = XGreedyRowState {
        diag_lower,
        diag_upper: tmp_diag_lower + band_len as i32 - 1,
        end1_reached: rng.below(4) == 0,
        end2_reached: rng.below(4) == 0,
        curr_extent: rng.range(-5, len1 + len2 + 5),
        curr_seq2_index: rng.range(-3, len2),
        curr_diag: rng.range(-3, 3),
        longest_match_run: rng.pick(&[0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 20, 1000]),
        seed_start_q: rng.range(-1, len1),
        seed_start_s: rng.range(-1, len2),
        fence_hit: rng.below(20) == 0,
    };
    let current = (0..band_len).map(|_| rng.range(-9, 999)).collect();
    Row {
        reverse,
        seq1,
        len1,
        seq2,
        len2,
        previous,
        current,
        tmp_diag_lower,
        diag_origin,
        xdrop_score,
        st,
    }
}

type RowFn = unsafe fn(
    bool,
    *const u8,
    *const u8,
    i32,
    i32,
    &[i32],
    &mut [i32],
    i32,
    i32,
    i32,
    &mut XGreedyRowState,
) -> bool;

/// The row `run` leaves: (return value, every cell, the state).
fn run(f: RowFn, r: &Row) -> (bool, Vec<i32>, XGreedyRowState) {
    let mut current = r.current.clone();
    let mut st = r.st;
    // SAFETY: the lengths lie within the vectors (bytes follow them) and the CPU has AVX2.
    let stopped = unsafe {
        f(
            r.reverse,
            r.seq1.as_ptr(),
            r.seq2.as_ptr(),
            r.len1,
            r.len2,
            &r.previous,
            &mut current,
            r.tmp_diag_lower,
            r.diag_origin,
            r.xdrop_score,
            &mut st,
        )
    };
    (stopped, current, st)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:537-553,571-614
/// Single distances on random rows: the fused build leaves the return value, every cell of
/// the row (stale cells included) and every field of the state that `x_greedy_row_avx2`
/// (`LOSAT_X_GREEDYROWFUSE=0`) and the plain build leave. The generator must reach the
/// eight-byte test's outcomes (a slide length, a long slide, a stop at the fence byte, offsets
/// outside a sequence), rows of several blocks, rows with fewer than eight bytes or the fence
/// hit on entry, and fence stops.
#[test]
fn x_greedy_row_fused_matches_three_passes_on_random_rows() {
    if !avx2() {
        eprintln!("x_greedy_row_fused_matches_three_passes_on_random_rows: no AVX2, skipped");
        return;
    }
    let rows: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(20_000);
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    // slide length read off, longer than the run, fence byte, outside a sequence
    let mut cells = [0usize; 4];
    let (mut blocks, mut short, mut gated, mut stopped) = (0, 0, 0, 0);
    for n in 0..rows {
        let r = row(&mut rng);
        let what = format!(
            "row {n}: reverse={} len1={} len2={} band={} first={} xdrop_score={} st={:?}",
            r.reverse,
            r.len1,
            r.len2,
            r.current.len(),
            r.tmp_diag_lower - r.diag_origin,
            r.xdrop_score,
            r.st
        );
        let expected = run(x_greedy_row_avx2, &r);
        assert_eq!(run(x_greedy_row_fused_avx2, &r), expected, "{what}: fused");
        assert_eq!(run(x_greedy_row_any, &r), expected, "{what}: plain");

        blocks += usize::from(r.current.len() > X_GREEDY_BLOCK);
        short += usize::from(r.len1 < 8 || r.len2 < 8);
        gated += usize::from(r.st.fence_hit);
        stopped += usize::from(expected.0 && !r.st.fence_hit);
        if r.len1 < 8 || r.len2 < 8 || r.st.fence_hit {
            continue;
        }
        // what the eight-byte test sees on the surviving diagonals (with the run on entry)
        let first = r.tmp_diag_lower - r.diag_origin;
        for i in 0..r.current.len() {
            let p = &r.previous[i..i + 3];
            let s2 = p[2].max(p[1]).saturating_add(1).max(p[0]);
            let s1 = s2 + first + i as i32;
            if s2 < 0 || s1 + s2 < r.xdrop_score {
                continue;
            }
            if s1 > r.len1 - 8 || s2 > r.len2 - 8 {
                cells[3] += 1;
                continue;
            }
            let (at1, at2) = if r.reverse {
                ((r.len1 - 8 - s1) as usize, (r.len2 - 8 - s2) as usize)
            } else {
                (s1 as usize, s2 as usize)
            };
            let w1 = u64::from_le_bytes(r.seq1[at1..at1 + 8].try_into().unwrap());
            let w2 = u64::from_le_bytes(r.seq2[at2..at2 + 8].try_into().unwrap());
            let diff = (w1 ^ w2) | (w1 & 0xFCFC_FCFC_FCFC_FCFC);
            let index = if r.reverse {
                diff.leading_zeros()
            } else {
                diff.trailing_zeros()
            } as i32
                >> 3;
            let stop = if r.reverse { 7 - index } else { index };
            if index < 8 && r.seq2[at2 + stop as usize] == FENCE_SENTRY {
                cells[2] += 1;
            } else if index > r.st.longest_match_run.min(7) {
                cells[1] += 1;
            } else {
                cells[0] += 1;
            }
        }
    }
    let summary = format!(
        "rows {rows}: several blocks {blocks}, fewer than eight bytes {short}, fence on \
             entry {gated}, fence stops {stopped}; diagonals: slide read off {}, longer than the \
             run {}, fence byte {}, outside {}",
        cells[0], cells[1], cells[2], cells[3]
    );
    eprintln!("x_greedy_row_fused_matches_three_passes_on_random_rows: {summary}");
    assert!(cells.iter().all(|&c| c * 10 > rows), "{summary}");
    assert!(blocks * 10 > rows && short * 20 > rows, "{summary}");
    assert!(gated * 40 > rows && stopped * 40 > rows, "{summary}");
}

/// What one `blast_greedy_align` call returns and leaves.
#[derive(Debug, PartialEq, Eq)]
struct Outcome {
    ret: i32,
    seq1_align_len: i32,
    seq2_align_len: i32,
    seed: (i32, i32, i32),
    fence_hit: bool,
    edit: String,
}

/// `blast_greedy_align` with an edit block (traceback: `x_greedy_row` runs every distance), or
/// score only with the fence hit on entry (the fused loop declines it: `x_greedy_row` runs
/// distance 1 and stops at its first surviving diagonal), with the build `kernel`.
fn align(
    kernel: XGreedyRowKernel,
    mem: &mut GreedyNonAffineMem,
    (seq1, len1, seq2, len2): (&[u8], i32, &[u8], i32),
    (reverse, xdrop, match_cost, mismatch_cost, max_dist): (bool, i32, i32, i32, i32),
    traceback: bool,
) -> Outcome {
    X_GREEDY_ROW_TEST_KERNEL.with(|k| k.set(Some(kernel)));
    let mut edit = GapPrelimEditBlock::new();
    // values the call overwrites
    let mut seed = GreedySeed {
        start_q: -7,
        start_s: -7,
        match_length: -7,
    };
    let mut fence_hit = !traceback;
    let (mut seq1_align_len, mut seq2_align_len) = (-5, -5);
    let ret = blast_greedy_align(
        seq1,
        len1,
        seq2,
        len2,
        reverse,
        xdrop,
        match_cost,
        mismatch_cost,
        &mut seq1_align_len,
        &mut seq2_align_len,
        max_dist,
        mem,
        traceback.then_some(&mut edit),
        &mut fence_hit,
        &mut seed,
    );
    X_GREEDY_ROW_TEST_KERNEL.with(|k| k.set(None));
    Outcome {
        ret,
        seq1_align_len,
        seq2_align_len,
        seed: (seed.start_q, seed.start_s, seed.match_length),
        fence_hit,
        edit: format_gap_prelim_edit_block_for_trace(&edit),
    }
}

/// The scratch two calls leave, whole: both rows, `max_score` and the traceback rows.
fn assert_same_memory(a: &GreedyNonAffineMem, b: &GreedyNonAffineMem, what: &str) {
    assert!(a.last_seq2_off == b.last_seq2_off, "{what}: rows differ");
    assert!(a.max_score == b.max_score, "{what}: max_score differs");
    assert!(
        a.traceback_used == b.traceback_used && a.traceback_pool == b.traceback_pool,
        "{what}: traceback rows differ"
    );
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/greedy_align.c:509-537,682-751
/// Sequences of whole alignments on one scratch per build, so that the cells earlier calls
/// left matter: traceback calls and score-only calls with the fence hit on entry, on related
/// sequences of up to about 3,000 bases with megablast's and other costs. The fused build
/// returns the values, the edit script and the seed of `x_greedy_row_avx2`
/// (`LOSAT_X_GREEDYROWFUSE=0`) and of the plain build, and leaves the same scratch.
#[test]
fn x_greedy_row_fused_matches_three_passes_on_alignments() {
    if !avx2() {
        eprintln!("x_greedy_row_fused_matches_three_passes_on_alignments: no AVX2, skipped");
        return;
    }
    let histories: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .map_or(300, |n: usize| (n / 10).max(1));
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let (mut calls, mut converged, mut stopped, mut fences, mut long) = (0, 0, 0, 0, 0);
    for history in 0..histories {
        let mut three = GreedyNonAffineMem::new();
        let mut fused = GreedyNonAffineMem::new();
        let mut plain = GreedyNonAffineMem::new();
        for step in 0..1 + rng.below(8) {
            let len = match rng.below(4) {
                0 => rng.below(40),
                1 => 40 + rng.below(300),
                _ => 340 + rng.below(2700),
            } as usize;
            let reverse = rng.below(2) == 1;
            let (seq1, len1, seq2, len2) = related(&mut rng, len, reverse);
            // megablast's costs (reward 1 doubled), and costs with match_cost / 2 != 1
            let (match_cost, mismatch_cost) =
                rng.pick(&[(2, 4), (2, 4), (2, 4), (2, 6), (2, 2), (4, 6), (6, 8)]);
            let xdrop = rng.pick(&[0, 1, 5, 13, 20, 26, 27, 40, 54, 100, 300]);
            let max_dist = match rng.below(3) {
                0 => 1 + rng.below(12) as i32,
                _ => 1000.min(len2 / 2 + 1),
            };
            let traceback = rng.below(6) != 0;
            let input = (&seq1[..], len1, &seq2[..], len2);
            let costs = (reverse, xdrop, match_cost, mismatch_cost, max_dist);
            let what = format!(
                "history {history} step {step}: len1={len1} len2={len2} reverse={reverse} \
                     xdrop={xdrop} costs={match_cost}/{mismatch_cost} max_dist={max_dist} \
                     traceback={traceback}"
            );
            let expected = align(XGreedyRowKernel::Avx2, &mut three, input, costs, traceback);
            let got = align(XGreedyRowKernel::Fused, &mut fused, input, costs, traceback);
            assert_eq!(got, expected, "{what}: fused");
            assert_same_memory(&fused, &three, &format!("{what}: fused"));
            let got = align(XGreedyRowKernel::Plain, &mut plain, input, costs, traceback);
            assert_eq!(got, expected, "{what}: plain");
            assert_same_memory(&plain, &three, &format!("{what}: plain"));

            calls += 1;
            long += usize::from(len >= 340);
            if expected.fence_hit && traceback {
                fences += 1;
            } else if expected.ret < 0 {
                stopped += 1;
            } else {
                converged += 1;
            }
        }
    }
    let summary = format!(
        "calls {calls}: converged {converged}, past max_dist {stopped}, fence {fences}; \
             340+ template bases {long}"
    );
    eprintln!("x_greedy_row_fused_matches_three_passes_on_alignments: {summary}");
    assert!(converged * 3 > calls && long * 3 > calls, "{summary}");
    assert!(stopped * 40 > calls && fences * 40 > calls, "{summary}");
}

/// No NCBI counterpart: the row kernels handle a kept diagonal past the end of `seq1`
/// (`x_past_seq1_end`) in place, and this test is the independent check of that branch. On
/// random lengths, bases, directions and offsets (`seq1_index` in `len1 + 1..=len1 + 5`,
/// `seq2_index` inside `seq2`), with and without the fence sentry at the byte that
/// `x_past_seq1_end` reads: it is true exactly when the sentry is absent, and then the
/// inline update (the cell is `seq2_index`, `curr_extent` and `best_cell` change only when
/// `seq1_index + seq2_index` is larger, nothing else) equals what `x_greedy_cell_slow` leaves
/// from the same inputs; with the sentry `x_greedy_cell_slow` stops at the fence.
#[test]
fn past_seq1_end_branch_equals_x_greedy_cell_slow() {
    let mut rng = Rng(0x005E_ED0F_E1D5_u64);
    let (mut inline, mut fenced) = (0usize, 0usize);
    for _ in 0..4000 {
        let len1 = rng.range(8, 48);
        let len2 = rng.range(8, 48);
        let seq1: Vec<u8> = (0..len1).map(|_| rng.below(4) as u8).collect();
        let mut seq2: Vec<u8> = (0..len2).map(|_| rng.below(4) as u8).collect();
        let reverse = rng.below(2) == 0;
        let seq2_index = rng.range(0, len2 - 1);
        let seq1_index = rng.range(len1 + 1, len1 + 5);
        let sentry = rng.below(2) == 0;
        let at = if reverse {
            len2 - 1 - seq2_index
        } else {
            seq2_index
        } as usize;
        if sentry {
            seq2[at] = FENCE_SENTRY;
        }
        let curr_extent = rng.range(-5, len1 + len2 + 10);
        let longest_match_run = rng.range(0, 5);
        let cell = rng.range(0, 7) as usize;
        let k = rng.range(-20, 20);
        let (seed_q, seed_s) = (rng.range(-1, len1), rng.range(-1, len2));
        let (diag_lower, best_cell, end1_cell) = (rng.range(-9, 9), rng.below(8) as usize, 3);
        let (end1, end2) = (rng.below(2) == 0, rng.below(2) == 0);
        let make = || XGreedyRowCold {
            diag_lower,
            end1_reached: end1,
            end2_reached: end2,
            curr_extent,
            best_cell,
            end1_cell,
            longest_match_run,
            seed_start_q: seed_q,
            seed_start_s: seed_s,
            fence_hit: false,
            stopped: false,
        };
        let snapshot = |c: &XGreedyRowCold| {
            (
                c.diag_lower,
                c.end1_reached,
                c.end2_reached,
                c.curr_extent,
                c.best_cell,
                c.end1_cell,
                c.longest_match_run,
                c.seed_start_q,
                c.seed_start_s,
                c.fence_hit,
                c.stopped,
            )
        };
        let before = [-9i32; 8];
        let past = unsafe {
            if reverse {
                x_past_seq1_end::<true>(seq2.as_ptr(), len1, len2, seq1_index, seq2_index)
            } else {
                x_past_seq1_end::<false>(seq2.as_ptr(), len1, len2, seq1_index, seq2_index)
            }
        };
        assert_eq!(past, !sentry, "len1 {len1} len2 {len2} reverse {reverse}");

        let mut cold = make();
        let mut current = before;
        unsafe {
            let (s1, s2) = (seq1.as_ptr(), seq2.as_ptr());
            let cur = current.as_mut_ptr();
            if reverse {
                x_greedy_cell_slow::<true>(
                    s1, s2, len1, len2, seq1_index, seq2_index, cell, k, cur, &mut cold,
                );
            } else {
                x_greedy_cell_slow::<false>(
                    s1, s2, len1, len2, seq1_index, seq2_index, cell, k, cur, &mut cold,
                );
            }
        }
        if !past {
            fenced += 1;
            assert!(
                cold.fence_hit && cold.stopped,
                "the sentry stops the distance"
            );
            continue;
        }
        inline += 1;
        // what the inline branch does
        let mut expected = make();
        let mut expected_current = before;
        expected_current[cell] = seq2_index;
        if seq1_index + seq2_index > expected.curr_extent {
            expected.curr_extent = seq1_index + seq2_index;
            expected.best_cell = cell;
        }
        assert_eq!(
            snapshot(&cold),
            snapshot(&expected),
            "len1 {len1} len2 {len2} reverse {reverse} seq1_index {seq1_index} \
                 seq2_index {seq2_index}"
        );
        assert_eq!(current, expected_current);
    }
    assert!(
        inline > 500 && fenced > 500,
        "inline {inline}, fenced {fenced}"
    );
}
