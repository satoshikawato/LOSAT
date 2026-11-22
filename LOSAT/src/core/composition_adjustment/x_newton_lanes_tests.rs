use super::super::adjust_scores::{
    build_matrix_info, read_aa_composition, x_newton_input, x_optimize_reference_for_test,
    BlastCompositionWorkspace,
};
use super::super::redo_alignment::BlastCompoAdjustMode;
use super::*;
use crate::config::ScoringMatrix;

// No NCBI counterpart: random number generator of the test; it does not change any value NCBI computes.
fn xorshift(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1385-1399
// Test problems: realistic ones through `x_newton_input` (BLOSUM62 joint probabilities as q,
// compositions of random sequences with biased letter frequencies, pseudocounts), and harder
// ones with random q, compositions and relative entropies.
fn test_problems() -> Vec<NewtonInput> {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, 0.3176).unwrap();
    let workspace = BlastCompositionWorkspace::new_blosum62();
    let mut problems = Vec::new();
    let mut shared_q = None;
    let mut attempts = 0;
    while problems.len() < 40 && attempts < 400 {
        attempts += 1;
        let mut weights = [0.0f64; 20];
        for w in weights.iter_mut() {
            *w = 0.05 + xorshift(&mut state).powi(3) * 4.0;
        }
        let total: f64 = weights.iter().sum();
        let mut sequences = Vec::new();
        for _ in 0..2 {
            let len = 30 + (xorshift(&mut state) * 900.0) as usize;
            let mut seq = Vec::with_capacity(len);
            for _ in 0..len {
                let mut u = xorshift(&mut state) * total;
                let mut letter = 0usize;
                while letter < 19 && u >= weights[letter] {
                    u -= weights[letter];
                    letter += 1;
                }
                // NCBIstdaa codes of the 20 true amino acids (A..Y without B, X, Z, U, *).
                const TRUE_AA: [u8; 20] = [
                    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
                ];
                seq.push(TRUE_AA[letter]);
            }
            sequences.push(seq);
            // the second sequence uses other weights
            for w in weights.iter_mut() {
                *w = 0.05 + xorshift(&mut state).powi(2) * 2.0;
            }
        }
        let query = read_aa_composition(&sequences[0]);
        let subject = read_aa_composition(&sequences[1]);
        let mode = if attempts % 3 == 0 {
            BlastCompoAdjustMode::ForceFullMatrixAdjust
        } else {
            BlastCompoAdjustMode::CompositionMatrixAdjust
        };
        if let Some(input) = x_newton_input(
            &matrix_info,
            &query,
            sequences[0].len() as i32,
            &subject,
            sequences[1].len() as i32,
            mode,
            &workspace,
            &mut shared_q,
        ) {
            problems.push(input);
        }
    }
    assert!(
        problems.len() >= 30,
        "realistic problems: {}",
        problems.len()
    );
    let blosum_q = problems[0].q.clone();
    // Random joint distributions, compositions and relative entropies.
    for case in 0..24 {
        let mut input = problems[case % problems.len()].clone();
        if case % 2 == 0 {
            let mut total = 0.0;
            let q = std::rc::Rc::make_mut(&mut input.q);
            for v in q.iter_mut() {
                *v = 1.0e-4 + xorshift(&mut state).powi(2);
                total += *v;
            }
            for v in q.iter_mut() {
                *v /= total;
            }
        } else {
            input.q = blosum_q.clone();
        }
        let (mut sr, mut sc) = (0.0, 0.0);
        for k in 0..N {
            input.row[k] = 0.002 + xorshift(&mut state).powi(3);
            input.col[k] = 0.002 + xorshift(&mut state).powi(3);
            sr += input.row[k];
            sc += input.col[k];
        }
        for k in 0..N {
            input.row[k] /= sr;
            input.col[k] /= sc;
        }
        input.relative_entropy = match case % 4 {
            0 => 0.05 + xorshift(&mut state),
            1 => 0.44,
            2 => 1.5 + 2.0 * xorshift(&mut state),
            _ => 1.0e-3,
        };
        problems.push(input);
    }
    problems
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// Expected results of each test problem: the exact module and the reference port (which must
// agree with each other).
fn expected(problems: &[NewtonInput]) -> Vec<NewtonOutput> {
    problems
        .iter()
        .map(|input| {
            let mut exact = NewtonOutput::new();
            exact.status = optimize_target_frequencies_exact(
                &mut exact.x,
                &input.q[..],
                &input.row,
                &input.col,
                input.relative_entropy,
            );
            let mut reference = NewtonOutput::new();
            reference.status = x_optimize_reference_for_test(
                &mut reference.x,
                &input.q[..],
                &input.row,
                &input.col,
                input.relative_entropy,
            );
            assert_eq!(exact.status, reference.status);
            assert!(same_bits(&exact.x, &reference.x));
            reference
        })
        .collect()
}

// No NCBI counterpart: bitwise comparison of a batch result with the expected results.
fn check_batch(label: &str, batch: &[usize], outputs: &[NewtonOutput], want: &[NewtonOutput]) {
    for (position, &problem) in batch.iter().enumerate() {
        assert_eq!(
            outputs[position].status, want[problem].status,
            "{label}: problem {problem} status"
        );
        for k in 0..NN {
            assert_eq!(
                outputs[position].x[k].to_bits(),
                want[problem].x[k].to_bits(),
                "{label}: problem {problem} x[{k}]: {} vs {}",
                outputs[position].x[k],
                want[problem].x[k]
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-758,762-768
// The test solves batches of every size 1..=13 and one of 64 (mixed iteration counts, lanes
// refilled as problems end) and compares status and all 400 values by bits with the exact
// module and the reference port, then repeats batches with the literal factorisation forced on
// one lane in one round.
#[test]
fn x_newton_lanes_match_exact_and_reference_bitwise() {
    let problems = test_problems();
    let want = expected(&problems);
    let failed = want.iter().filter(|o| o.status != 0).count();
    eprintln!(
        "[x_newton_lanes test] problems={} status1={}",
        problems.len(),
        failed
    );
    let mut start = 0usize;
    let mut sizes: Vec<usize> = (1..=13).collect();
    sizes.push(64);
    for &size in &sizes {
        let batch: Vec<usize> = (0..size)
            .map(|i| (start + 7 * i) % problems.len())
            .collect();
        start += 5;
        let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
        let mut outputs = vec![NewtonOutput::new(); size];
        solve_batch(&inputs, &mut outputs);
        check_batch(&format!("batch {size}"), &batch, &outputs, &want);
    }
    // The plain (non-AVX2) compilation of the lane body.
    {
        let batch: Vec<usize> = (0..9).map(|i| (3 * i + 1) % problems.len()).collect();
        let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
        let mut outputs = vec![NewtonOutput::new(); batch.len()];
        let mut st = new_lane_state();
        solve_lanes_plain(&inputs, &mut outputs, &mut st);
        check_batch("plain", &batch, &outputs, &want);
    }
    // Literal factorisation forced on one lane in one round.
    for (round, lane) in [(0usize, 0usize), (1, 3), (2, 1), (5, 2), (9, 0)] {
        let batch: Vec<usize> = (0..10).map(|i| (round + 11 * i) % problems.len()).collect();
        let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
        let mut outputs = vec![NewtonOutput::new(); batch.len()];
        FORCE_LITERAL.with(|f| f.set(Some((round, lane))));
        solve_batch(&inputs, &mut outputs);
        FORCE_LITERAL.with(|f| f.set(None));
        check_batch(
            &format!("literal round {round} lane {lane}"),
            &batch,
            &outputs,
            &want,
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:752-758,762-768
// The eight-lane batches (aarch64): every size 1..=21 and one of 64, with the move
// to four lanes and without it, the plain compilation, and the literal factorisation forced on one
// lane (any of the eight) in one round; status and all 400 values compared by bits with the exact
// module and the reference port.
#[test]
fn x_newton_lanes8_match_exact_and_reference_bitwise() {
    let problems = test_problems();
    let want = expected(&problems);
    let mut start = 0usize;
    let mut sizes: Vec<usize> = (1..=21).collect();
    sizes.push(64);
    for hand_over_at in [HAND_OVER_LANES, 0] {
        for &size in &sizes {
            let batch: Vec<usize> = (0..size)
                .map(|i| (start + 7 * i) % problems.len())
                .collect();
            start += 5;
            let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
            let mut outputs = vec![NewtonOutput::new(); size];
            solve_batch_with(&inputs, &mut outputs, Some(hand_over_at));
            check_batch(
                &format!("lanes8 hand-over {hand_over_at} batch {size}"),
                &batch,
                &outputs,
                &want,
            );
        }
    }
    // The plain (non-AVX2) compilation of the eight-lane body; the batch ends on four lanes.
    {
        let batch: Vec<usize> = (0..17).map(|i| (3 * i + 2) % problems.len()).collect();
        let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
        let mut outputs = vec![NewtonOutput::new(); batch.len()];
        let mut st8 = new_lane_state::<F8>();
        let mut st4 = new_lane_state::<F4>();
        let [rounds8, busy8, rounds4, busy4] =
            solve_lanes8_plain(&inputs, &mut outputs, &mut st8, &mut st4, HAND_OVER_LANES);
        check_batch("lanes8 plain", &batch, &outputs, &want);
        assert!(rounds8 > 0 && rounds4 > 0, "rounds {rounds8} + {rounds4}");
        assert!(busy8 <= rounds8 * LANES8 as u64 && busy4 <= rounds4 * LANES as u64);
    }
    // Literal factorisation forced on one lane in one round.
    for (round, lane) in [
        (0usize, 0usize),
        (1, 7),
        (2, 5),
        (4, 3),
        (6, 6),
        (9, 1),
        (12, 4),
    ] {
        let batch: Vec<usize> = (0..13).map(|i| (round + 11 * i) % problems.len()).collect();
        let inputs: Vec<NewtonInput> = batch.iter().map(|&p| problems[p].clone()).collect();
        let mut outputs = vec![NewtonOutput::new(); batch.len()];
        FORCE_LITERAL.with(|f| f.set(Some((round, lane))));
        solve_batch_with(&inputs, &mut outputs, Some(HAND_OVER_LANES));
        FORCE_LITERAL.with(|f| f.set(None));
        check_batch(
            &format!("lanes8 literal round {round} lane {lane}"),
            &batch,
            &outputs,
            &want,
        );
    }
}

// No NCBI counterpart: the store returns a result only for a bitwise-equal input; it does not
// change any value NCBI computes.
#[test]
fn x_newton_lanes_store_returns_only_exact_input_matches() {
    let problems = test_problems();
    let want = expected(&problems[..5]);
    prefetch(problems[..5].to_vec());
    let mut x = [0.0f64; NN];
    // A different relative entropy (one ulp) misses.
    let p = &problems[2];
    let re = f64::from_bits(p.relative_entropy.to_bits() + 1);
    assert_eq!(take(&mut x, &p.q[..], &p.row, &p.col, re), None);
    // Exact inputs hit once each, in any order.
    for i in [3usize, 0, 4, 1, 2] {
        let p = &problems[i];
        let status = take(&mut x, &p.q[..], &p.row, &p.col, p.relative_entropy);
        assert_eq!(status, Some(want[i].status));
        assert!(same_bits(&x, &want[i].x));
        assert_eq!(
            take(&mut x, &p.q[..], &p.row, &p.col, p.relative_entropy),
            None
        );
    }
}

// No NCBI counterpart: LOSAT_X_NEWTONDROP; a prefetch drops this thread's results of earlier
// prefetches not taken (with `=0` they stay); it does not change any value NCBI computes.
#[test]
fn x_newton_drop_empties_the_store_at_the_next_prefetch() {
    let problems = test_problems();
    let want = expected(&problems[..8]);
    let mut x = [0.0f64; NN];
    let take_at = |x: &mut [f64; NN], i: usize| {
        let p = &problems[i];
        take(x, &p.q[..], &p.row, &p.col, p.relative_entropy)
    };
    // the same calls with the switch (`drop`) and with `=0`, in this process, with the results
    // kept as their batch (LOSAT_X_NEWTONBATCH, only with `drop`) or stored one by one
    for (drop, batch) in [(true, true), (true, false), (false, true), (false, false)] {
        x_prefetch(Vec::new(), true, true);
        x_prefetch(problems[..5].to_vec(), drop, batch);
        assert_eq!(
            BATCH.with(|kept| kept.borrow().is_some()),
            drop && batch,
            "a batch is kept only with drop, drop {drop} batch {batch}"
        );
        for i in [0usize, 2, 4] {
            assert_eq!(
                take_at(&mut x, i),
                Some(want[i].status),
                "drop {drop} batch {batch}"
            );
            assert!(same_bits(&x, &want[i].x));
        }
        x_prefetch(problems[5..8].to_vec(), drop, batch);
        // problems 1 and 3 were not taken before the second prefetch
        for i in [1usize, 3] {
            let kept = take_at(&mut x, i);
            if drop {
                assert_eq!(kept, None, "problem {i} dropped, batch {batch}");
            } else {
                assert_eq!(
                    kept,
                    Some(want[i].status),
                    "problem {i} kept with =0, batch {batch}"
                );
                assert!(same_bits(&x, &want[i].x));
            }
        }
        for i in 5..8 {
            assert_eq!(
                take_at(&mut x, i),
                Some(want[i].status),
                "drop {drop} batch {batch}"
            );
            assert!(same_bits(&x, &want[i].x));
        }
        // an empty prefetch drops what is left as well
        x_prefetch(problems[..2].to_vec(), drop, batch);
        assert_eq!(
            take_at(&mut x, 0),
            Some(want[0].status),
            "drop {drop} batch {batch}"
        );
        x_prefetch(Vec::new(), drop, batch);
        assert_eq!(
            take_at(&mut x, 1).is_none(),
            drop,
            "problem 1 after an empty prefetch, drop {drop} batch {batch}"
        );
    }
}

// No NCBI counterpart: LOSAT_X_NEWTONBATCH; the kept batch gives each result once, in any order,
// to the call with its input's bits (also a repeated input), searches on from its first result
// not taken, stays out of the store, and goes when its last result is taken; it does not change
// any value NCBI computes.
#[test]
fn x_newton_batch_gives_each_result_once_and_goes_when_all_are_taken() {
    let problems = test_problems();
    let mut inputs = problems[..5].to_vec();
    // problem 2 twice: both calls with its bits get its result
    inputs.push(problems[2].clone());
    let want = expected(&inputs);
    let mut x = [0.0f64; NN];
    let take_at = |x: &mut [f64; NN], i: usize| {
        let p = &inputs[i];
        take(x, &p.q[..], &p.row, &p.col, p.relative_entropy)
    };
    x_prefetch(Vec::new(), true, true);
    x_prefetch(inputs.clone(), true, true);
    assert!(
        STORE.with(|store| store.borrow().is_empty()),
        "the store stays empty"
    );
    let next = || BATCH.with(|kept| kept.borrow().as_ref().map(|batch| batch.next));
    for i in [1usize, 2, 3] {
        assert_eq!(take_at(&mut x, i), Some(want[i].status), "problem {i}");
        assert!(same_bits(&x, &want[i].x));
        assert_eq!(next(), Some(0), "problem 0 is not taken yet");
    }
    // results 0 to 3 taken: the search starts at 4
    assert_eq!(take_at(&mut x, 0), Some(want[0].status));
    assert!(same_bits(&x, &want[0].x));
    assert_eq!(next(), Some(4));
    // the second call with the bits of problem 2 gets the repeated one
    assert_eq!(take_at(&mut x, 5), Some(want[5].status));
    assert!(same_bits(&x, &want[5].x));
    assert_eq!(take_at(&mut x, 2), None, "no result of problem 2 is left");
    assert_eq!(next(), Some(4), "problem 4 is left");
    assert_eq!(take_at(&mut x, 4), Some(want[4].status));
    assert!(same_bits(&x, &want[4].x));
    assert!(
        BATCH.with(|kept| kept.borrow().is_none()),
        "the batch goes with its last result"
    );
    assert_eq!(take_at(&mut x, 1), None);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// The ignored test times `optimize_target_frequencies_exact` (one problem per call) against
// four lanes and eight lanes (with and without the move to four lanes, at every batch size) on
// the realistic test problems and checks the results by bits.
/// Timing (run with --release --ignored --nocapture).
#[test]
#[ignore]
fn x_bench_newton_lanes() {
    let problems: Vec<NewtonInput> = test_problems().into_iter().take(36).collect();
    let reps = 10;
    // rounds of the four forms in turn; the fastest round of each form is printed
    let rounds = 9;
    for batch in [5usize, 6, 7, 8, 9, 10, 11, 12, 16, 21, 32, 36] {
        let inputs: Vec<NewtonInput> = problems[..batch].to_vec();
        let mut exact = vec![NewtonOutput::new(); batch];
        let mut lanes = vec![NewtonOutput::new(); batch];
        let mut best = [f64::INFINITY; 4];
        for _ in 0..rounds {
            let t0 = std::time::Instant::now();
            for _ in 0..reps {
                for (input, output) in inputs.iter().zip(exact.iter_mut()) {
                    output.status = optimize_target_frequencies_exact(
                        &mut output.x,
                        &input.q[..],
                        &input.row,
                        &input.col,
                        input.relative_entropy,
                    );
                }
            }
            best[0] = best[0].min(t0.elapsed().as_secs_f64() / (reps * batch) as f64 * 1e6);
            for (form, t) in best[1..].iter_mut().enumerate() {
                let t0 = std::time::Instant::now();
                for _ in 0..reps {
                    match form {
                        0 => solve_batch_with(&inputs, &mut lanes, None),
                        1 => solve_batch8(&inputs, &mut lanes, HAND_OVER_LANES),
                        _ => solve_batch8(&inputs, &mut lanes, 0),
                    }
                }
                *t = t.min(t0.elapsed().as_secs_f64() / (reps * batch) as f64 * 1e6);
                for (a, b) in exact.iter().zip(&lanes) {
                    assert_eq!(a.status, b.status);
                    assert!(same_bits(&a.x, &b.x));
                }
            }
        }
        eprintln!(
            "[BENCH] batch {batch}: exact {:.1} us/problem, lanes {:.1}, lanes8 {:.1}, lanes8 without move {:.1} us/problem",
            best[0], best[1], best[2], best[3]
        );
    }
}
