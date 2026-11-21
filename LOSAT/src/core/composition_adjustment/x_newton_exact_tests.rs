use super::*;

// No NCBI counterpart: random number generator of the test; it does not change any value NCBI computes.
fn xorshift(state: &mut u64) -> f64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 11) as f64 / (1u64 << 53) as f64
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:171-177,179-183
/// Reference `Nlm_FactorLtriangPosDef` + `Nlm_SolveLtriangPosDef` on a dense copy.
fn reference_factor_solve(w: &[[f64; M]; M], b: &[f64; M]) -> ([[f64; M]; M], [f64; M]) {
    let mut a = *w;
    for i in 0..M {
        for j in 0..i {
            let mut temp = a[i][j];
            for k in 0..j {
                temp -= a[i][k] * a[j][k];
            }
            a[i][j] = temp / a[j][j];
        }
        let mut temp = a[i][i];
        for k in 0..i {
            temp -= a[i][k] * a[i][k];
        }
        a[i][i] = temp.sqrt();
    }
    let mut x = *b;
    for i in 0..M {
        let mut temp = x[i];
        for j in 0..i {
            temp -= a[i][j] * x[j];
        }
        x[i] = temp / a[i][i];
    }
    for j in (0..M).rev() {
        x[j] /= a[j][j];
        for i in 0..j {
            x[i] -= a[j][i] * x[j];
        }
    }
    (a, x)
}

// No NCBI counterpart: random test input for `dinv` and the gradient; it does not change any value NCBI computes.
fn random_parts(state: &mut u64, wide: bool) -> ([f64; NN], [f64; NN]) {
    let mut dinv = [0.0f64; NN];
    for v in dinv.iter_mut() {
        *v = if wide {
            2f64.glibc_pow(xorshift(state) * 40.0 - 30.0)
        } else {
            1.0e-4 + xorshift(state) * 0.05
        };
    }
    let mut grad = [0.0f64; NN];
    for v in grad.iter_mut() {
        *v = xorshift(state) * 4.0 - 2.0;
    }
    (dinv, grad)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:132-135,527-531
// Builds `W` the way `ScaledSymmetricProductA` and `FactorReNewtonSystem` do, for the test.
fn build_w(dinv: &[f64; NN], grad: &[f64; NN]) -> ([f64; N], [f64; N], [f64; M], [[f64; M]; M]) {
    let mut wdiag = [0.0f64; N];
    let mut wrow = [0.0f64; N];
    for i in 0..N {
        for j in 0..N {
            wdiag[j] += dinv[i * N + j];
        }
    }
    for i in 1..N {
        let mut s = 0.0;
        for j in 0..N {
            s += dinv[i * N + j];
        }
        wrow[i] = s;
    }
    let mut ws = [0.0f64; NN];
    for k in 0..NN {
        ws[k] = dinv[k] * grad[k];
    }
    let mut w39 = [0.0f64; M];
    for i in 0..N {
        for j in 0..N {
            w39[j] += 1.0 * ws[i * N + j];
        }
    }
    for i in 1..N {
        for j in 0..N {
            w39[i + N - 1] += 1.0 * ws[i * N + j];
        }
    }
    let mut s = 0.0;
    for k in 0..NN {
        s += grad[k] * ws[k];
    }
    w39[M - 1] = s;
    let mut w = [[0.0f64; M]; M];
    for j in 0..N {
        w[j][j] = wdiag[j];
    }
    for i in 1..N {
        for j in 0..N {
            w[i + N - 1][j] = 0.0 + dinv[i * N + j];
        }
        w[i + N - 1][i + N - 1] = wrow[i];
    }
    w[M - 1] = w39;
    (wdiag, wrow, w39, w)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149,151-155,171-177
// The test factors and solves 400 random Newton matrices with the structured code and with a
// copy of these loops, and compares every entry by its bits.
#[test]
fn x_structured_factor_and_solve_match_reference_bitwise_on_random_newton_matrices() {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut structured = 0usize;
    for round in 0..400 {
        let (dinv, grad) = random_parts(&mut state, round % 2 == 1);
        let (wdiag, wrow, w39, w) = build_w(&dinv, &grad);
        let parts = Wparts {
            wdiag: &wdiag,
            wrow: &wrow,
            dinv: &dinv,
            w39: &w39,
        };
        let mut b = [0.0f64; M];
        for v in b.iter_mut() {
            *v = xorshift(&mut state) * 2.0 - 1.0;
        }
        let (l_ref, x_ref) = reference_factor_solve(&w, &b);
        let mut f = Factor::new();
        let dinv_ok = dinv.iter().all(|&v| v > 0.0 && v < f64::INFINITY);
        if structure_ok(&parts, dinv_ok) {
            structured += 1;
            factor_structured(&parts, &mut f);
        } else {
            factor_literal(&parts, &mut f);
        }
        for i in 0..M {
            for j in 0..=i {
                assert_eq!(
                    f.get(i, j).to_bits(),
                    l_ref[i][j].to_bits(),
                    "round {round}: L[{i}][{j}] differs: {} vs {}",
                    f.get(i, j),
                    l_ref[i][j]
                );
            }
        }
        let mut xs = b;
        solve_structured(&mut xs, &f);
        for i in 0..M {
            assert_eq!(
                xs[i].to_bits(),
                x_ref[i].to_bits(),
                "round {round}: x[{i}] differs"
            );
        }
        // The literal path must agree too (it is the fallback).
        let mut g = Factor::new();
        factor_literal(&parts, &mut g);
        for i in 0..M {
            for j in 0..=i {
                assert_eq!(g.get(i, j).to_bits(), l_ref[i][j].to_bits());
            }
        }
    }
    assert!(
        structured > 300,
        "structured path exercised {structured} times"
    );
}
