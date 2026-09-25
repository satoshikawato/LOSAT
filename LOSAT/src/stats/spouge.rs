use crate::config::{ProteinScoringSpec, ScoringMatrix};

use super::tables::KarlinParams;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:94-112
// ```c
// typedef struct Blast_GumbelBlk {
//       double  Lambda;    /**< the unscaled Lambda value */
//       double  C;
//       double  G;         /**< G is the total penalty for extension */
//       double  a;         /**< avg(L) = a     y + b    */
//       double  Alpha;     /**< var(L) = alpha y + beta */
//       double  Sigma;     /**< cov(L) = sigma y + tau  */
//       double  a_un;      /**< Ungapped a */
//       double  Alpha_un;  /**< Ungapped alpha */
//       double  b;         /**< 2*G*(a_un - a) */
//       double  Beta;      /**< 2*G*(alpha_un - alpha) */
//       double  Tau;       /**< 2*G*(alpha_un - Sigma) */
//       Int8 db_length;    /**< total length of database */
// } Blast_GumbelBlk;
// ```
#[derive(Debug, Clone, Copy)]
pub struct BlastGumbelBlk {
    pub lambda: f64,
    pub c: f64,
    pub g: f64,
    pub a: f64,
    pub alpha: f64,
    pub sigma: f64,
    pub a_un: f64,
    pub alpha_un: f64,
    pub b: f64,
    pub beta: f64,
    pub tau: f64,
    pub db_length: i64,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:3696-3742
// ```c
// if (BLAST_Nint(values[index][0]) == gap_open &&
//    BLAST_Nint(values[index][1]) == gap_extend) {
//    gbp->Lambda = values[index][3];
//    gbp->C = values[index][8];
//    gbp->G = gap_open + gap_extend;
//    gbp->a = values[index][6];
//    gbp->Alpha = values[index][9];
//    gbp->Sigma = values[index][10];
//    gbp->a_un  = values[0][6];
//    gbp->Alpha_un = values[0][9];
//    gbp->b = 2.0 * gbp->G * (gbp->a_un - gbp->a);
//    gbp->Beta = 2.0 * gbp->G * (gbp->Alpha_un - gbp->Alpha);
//    gbp->Tau  = 2.0 * gbp->G * (gbp->Alpha_un - gbp->Sigma);
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:258-271
// ```c
// static array_of_8 blosum62_values[BLOSUM62_VALUES_MAX] = {
//     {(double) INT2_MAX, ..., 0.3176, ..., 0.7916, -3.2, 0.623757, 4.964660, 4.964660},
//     ...
//     {11, 1, (double) INT2_MAX, 0.267, 0.041, 0.14, 1.9, -30, 0.669720, 42.602800, 43.636200},
// };
// ```
pub fn lookup_protein_gumbel_params(
    spec: &ProteinScoringSpec,
    db_length: i64,
) -> Option<BlastGumbelBlk> {
    match (spec.matrix, spec.gap_open, spec.gap_extend) {
        // NCBI c++/src/algo/blast/core/blast_stat.c:183-192,3696-3742:
        // {INT2_MAX, INT2_MAX, ..., 0.2291, ..., 0.9113, ..., 9.611060, 9.611060},
        // {14, 2, ..., 0.195, ..., 1.9, ..., 0.685753, 60.736200, 61.102300};
        // gbp->b = 2*G*(a_un-a); gbp->Beta = 2*G*(Alpha_un-Alpha);
        // gbp->Tau = 2*G*(Alpha_un-Sigma);
        (ScoringMatrix::Blosum45, 14, 2) => {
            let g = (spec.gap_open + spec.gap_extend) as f64;
            let a = 1.9;
            let alpha = 60.736200;
            let sigma = 61.102300;
            let a_un = 0.9113;
            let alpha_un = 9.611060;
            Some(BlastGumbelBlk {
                lambda: 0.195,
                c: 0.685753,
                g,
                a,
                alpha,
                sigma,
                a_un,
                alpha_un,
                b: 2.0 * g * (a_un - a),
                beta: 2.0 * g * (alpha_un - alpha),
                tau: 2.0 * g * (alpha_un - sigma),
                db_length,
            })
        }
        (ScoringMatrix::Blosum62, 11, 1) => {
            let g = (spec.gap_open + spec.gap_extend) as f64;
            let a = 1.9;
            let alpha = 42.602800;
            let sigma = 43.636200;
            let a_un = 0.7916;
            let alpha_un = 4.964660;
            Some(BlastGumbelBlk {
                lambda: 0.267,
                c: 0.669720,
                g,
                a,
                alpha,
                sigma,
                a_un,
                alpha_un,
                b: 2.0 * g * (a_un - a),
                beta: 2.0 * g * (alpha_un - alpha),
                tau: 2.0 * g * (alpha_un - sigma),
                db_length,
            })
        }
        _ => None,
    }
}

// NCBI c++/src/algo/blast/core/boost_erf.c:1-35,59-260 is adapted
// from Boost. (C) Copyright John Maddock 2006.
// Boost Software License - Version 1.0 - August 17th, 2003
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated
// by a source language processor.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// NCBI boost_erf.c:76-93,252-255:
// if (z < 0) { ... }
// double ErfC(double z) { return ErfImpl(z, TRUE); }
// NCBI blast_stat.c:5216,5223 calls ErfC, not NCBI_ErfC.
#[inline]
fn erfc_ncbi(z: f64) -> f64 {
    erf_impl_ncbi(z, true)
}

// NCBI c++/src/algo/blast/core/boost_erf.c:71-249:
// static double ErfImpl(double z, int invert) {
//   if (z < 0) { if (!invert) return -ErfImpl(-z, invert);
//     else if (z < -0.5) return 2.0 - ErfImpl(-z, invert);
//     else return 1 + ErfImpl(-z, FALSE); }
//   ... if(invert) result = 1 - result; return result;
// }
#[inline]
fn erf_impl_ncbi(z: f64, mut invert: bool) -> f64 {
    if z < 0.0 {
        if !invert {
            return -erf_impl_ncbi(-z, invert);
        } else if z < -0.5 {
            return 2.0 - erf_impl_ncbi(-z, invert);
        } else {
            return 1.0 + erf_impl_ncbi(-z, false);
        }
    }
    let mut result;
    if z < 0.5 {
        // NCBI boost_erf.c:96-129:
        // if (z < 1e-10) result = z * 1.125 + z * c;
        // else result = z * (Y + p / q);
        if z < 1e-10 {
            result = if z == 0.0 {
                0.0
            } else {
                z * 1.125 + z * 0.003379167095512573896158903121545171688
            };
        } else {
            const Y: f64 = 1.044948577880859375;
            const P: [f64; 5] = [
                0.0834305892146531832907,
                -0.338165134459360935041,
                -0.0509990735146777432841,
                -0.00772758345802133288487,
                -0.000322780120964605683831,
            ];
            const Q: [f64; 5] = [
                1.0,
                0.455004033050794024546,
                0.0875222600142252549554,
                0.00858571925074406212772,
                0.000370900071787748000569,
            ];
            let zz = z * z;
            result = z * (Y + horner_ncbi(&P, zz) / horner_ncbi(&Q, zz));
        }
    } else if (invert && z < 28.0) || (!invert && z < 5.8) {
        // NCBI boost_erf.c:132-233:
        // invert = !invert; result = Y + p / q;
        // result *= expl(-z * z) / z;
        invert = !invert;
        let (y, p, q, x): (f64, &[f64], &[f64], f64);
        if z < 1.5 {
            const P: [f64; 6] = [
                -0.098090592216281240205,
                0.178114665841120341155,
                0.191003695796775433986,
                0.0888900368967884466578,
                0.0195049001251218801359,
                0.00180424538297014223957,
            ];
            const Q: [f64; 6] = [
                1.0,
                1.84759070983002217845,
                1.42628004845511324508,
                0.578052804889902404909,
                0.12385097467900864233,
                0.0113385233577001411017,
            ];
            y = 0.405935764312744140625;
            p = &P;
            q = &Q;
            x = z - 0.5;
        } else if z < 2.5 {
            const P: [f64; 6] = [
                -0.0243500476207698441272,
                0.0386540375035707201728,
                0.04394818964209516296,
                0.0175679436311802092299,
                0.00323962406290842133584,
                0.000235839115596880717416,
            ];
            const Q: [f64; 6] = [
                1.0,
                1.53991494948552447182,
                0.982403709157920235114,
                0.325732924782444448493,
                0.0563921837420478160373,
                0.00410369723978904575884,
            ];
            y = 0.50672817230224609375;
            p = &P;
            q = &Q;
            x = z - 1.5;
        } else if z < 4.5 {
            const P: [f64; 6] = [
                0.00295276716530971662634,
                0.0137384425896355332126,
                0.00840807615555585383007,
                0.00212825620914618649141,
                0.000250269961544794627958,
                0.113212406648847561139e-4,
            ];
            const Q: [f64; 6] = [
                1.0,
                1.04217814166938418171,
                0.442597659481563127003,
                0.0958492726301061423444,
                0.0105982906484876531489,
                0.000479411269521714493907,
            ];
            y = 0.5405750274658203125;
            p = &P;
            q = &Q;
            x = z - 3.5;
        } else {
            const P: [f64; 7] = [
                0.00628057170626964891937,
                0.0175389834052493308818,
                -0.212652252872804219852,
                -0.687717681153649930619,
                -2.5518551727311523996,
                -3.22729451764143718517,
                -2.8175401114513378771,
            ];
            const Q: [f64; 7] = [
                1.0,
                2.79257750980575282228,
                11.0567237927800161565,
                15.930646027911794143,
                22.9367376522880577224,
                13.5064170191802889145,
                5.48409182238641741584,
            ];
            y = 0.5579090118408203125;
            p = &P;
            q = &Q;
            x = 1.0 / z;
        }
        result = y + horner_ncbi(p, x) / horner_ncbi(q, x);
        result = exp_long_double_product_ncbi(result, z);
    } else {
        // NCBI boost_erf.c:235-240:
        // result = 0; invert = !invert;
        result = 0.0;
        invert = !invert;
    }
    if invert {
        result = 1.0 - result;
    }
    result
}

// NCBI c++/src/algo/blast/core/boost_erf.c:149,175,201,230:
// result *= expl(-z * z) / z;
// The argument is binary64, while expl and the following division/product
// evaluate in C long double before the assignment back to binary64.
// Double-double arithmetic preserves those intermediate guard bits in Rust.
#[derive(Clone, Copy)]
struct DoubleDouble {
    hi: f64,
    lo: f64,
}

impl DoubleDouble {
    #[inline]
    fn from_hi(hi: f64) -> Self {
        Self { hi, lo: 0.0 }
    }

    #[inline]
    fn add(self, other: Self) -> Self {
        let sum = self.hi + other.hi;
        let bb = sum - self.hi;
        let error = (self.hi - (sum - bb)) + (other.hi - bb) + self.lo + other.lo;
        let hi = sum + error;
        Self {
            hi,
            lo: error - (hi - sum),
        }
    }

    #[inline]
    fn mul(self, other: Self) -> Self {
        let product = self.hi * other.hi;
        let error = self.hi.mul_add(other.hi, -product)
            + self.hi * other.lo
            + self.lo * other.hi
            + self.lo * other.lo;
        let hi = product + error;
        Self {
            hi,
            lo: error - (hi - product),
        }
    }

    #[inline]
    fn div_f64(self, divisor: f64) -> Self {
        let quotient = self.hi / divisor;
        let product = Self::from_hi(quotient).mul(Self::from_hi(divisor));
        let residual = self.add(Self {
            hi: -product.hi,
            lo: -product.lo,
        });
        Self::from_hi(quotient).add(Self::from_hi(residual.hi / divisor))
    }
}

// NCBI c++/src/algo/blast/core/boost_erf.c:149,175,201,230:
// result *= expl(-z * z) / z;
// Range reduction and the exponential Taylor series reproduce the required
// long-double precision without a native C or NCBI runtime dependency.
#[inline]
fn exp_long_double_product_ncbi(result: f64, z: f64) -> f64 {
    let x = -z * z;
    let t = x / 512.0;
    let mut exponential = DoubleDouble::from_hi(1.0);
    let mut term = DoubleDouble::from_hi(1.0);
    for n in 1..=40 {
        term = term.mul(DoubleDouble::from_hi(t)).div_f64(n as f64);
        exponential = exponential.add(term);
    }
    for _ in 0..9 {
        exponential = exponential.mul(exponential);
    }
    exponential.div_f64(z).mul(DoubleDouble::from_hi(result)).hi
}

// NCBI c++/src/algo/blast/core/boost_erf.c:123-127,156-158,181-183:
// p = (((P[last] * x + P[last-1]) * x + ...) * x + P[1]) * x + P[0];
// q is evaluated in the same order.
#[inline]
fn horner_ncbi(coefficients: &[f64], x: f64) -> f64 {
    let mut result = *coefficients
        .last()
        .expect("NCBI polynomial has coefficients");
    for coefficient in coefficients[..coefficients.len() - 1].iter().rev() {
        result = result * x + coefficient;
    }
    result
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:5176-5231
// ```c
// double BLAST_SpougeStoE(Int4 y_, Blast_KarlinBlk* kbp, Blast_GumbelBlk* gbp,
//                         Int4 m_, Int4 n_)
// {
//     double scale_factor = kbp->Lambda / gbp->Lambda;
//     double db_scale_factor = (gbp->db_length) ?
//             (double)gbp->db_length/(double)n_ : 1.0;
//     ...
//     P_m_F = ErfC(-m_F / sqrt(2.0)) / 2.0;
//     ...
//     e_value = area * k_ * exp(-lambda_ * y_) * db_scale_factor;
// }
// ```
pub fn blast_spouge_stoe(
    score: i32,
    kbp: &KarlinParams,
    gbp: &BlastGumbelBlk,
    query_length: i32,
    subject_length: i32,
) -> f64 {
    let scale_factor = kbp.lambda / gbp.lambda;
    let db_scale_factor = if gbp.db_length != 0 {
        (gbp.db_length as f64) / (subject_length as f64)
    } else {
        1.0
    };

    let lambda = kbp.lambda;
    let k = kbp.k;
    let ai_hat = gbp.a * scale_factor;
    let bi_hat = gbp.b;
    let alphai_hat = gbp.alpha * scale_factor;
    let betai_hat = gbp.beta;
    let sigma_hat = gbp.sigma * scale_factor;
    let tau_hat = gbp.tau;

    let const_val = 0.398_942_280_401_432_7_f64;
    let y = score as f64;
    let m = query_length as f64;
    let n = subject_length as f64;

    let m_li_y = m - (ai_hat * y + bi_hat);
    let vi_y = (2.0 * alphai_hat / lambda).max(alphai_hat * y + betai_hat);
    let sqrt_vi_y = vi_y.sqrt();
    let m_f = m_li_y / sqrt_vi_y;
    let p_m_f = erfc_ncbi(-m_f / 2.0_f64.sqrt()) / 2.0;
    let p1 = m_li_y * p_m_f + sqrt_vi_y * const_val * (-0.5 * m_f * m_f).exp();

    let n_lj_y = n - (ai_hat * y + bi_hat);
    let vj_y = (2.0 * alphai_hat / lambda).max(alphai_hat * y + betai_hat);
    let sqrt_vj_y = vj_y.sqrt();
    let n_f = n_lj_y / sqrt_vj_y;
    let p_n_f = erfc_ncbi(-n_f / 2.0_f64.sqrt()) / 2.0;
    let p2 = n_lj_y * p_n_f + sqrt_vj_y * const_val * (-0.5 * n_f * n_f).exp();

    let c_y = (2.0 * sigma_hat / lambda).max(sigma_hat * y + tau_hat);
    let area = p1 * p2 + c_y * p_m_f * p_n_f;
    area * k * (-lambda * y).exp() * db_scale_factor
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:5236-5263
// ```c
// BLAST_SpougeEtoS(double e0,
//                  Blast_KarlinBlk* kbp,
//                  Blast_GumbelBlk* gbp,
//                  Int4 m, Int4 n)
// {
//     Int4 a=0, b, c;
//     double e;
//     double db_scale_factor = (gbp->db_length) ?
//             (double)gbp->db_length : 1.0;
//
//     b = MAX((int)(log(db_scale_factor/e0) / kbp->Lambda), 2);
//     e = BLAST_SpougeStoE(b, kbp, gbp, m, n);
//     ...
// }
// ```
pub fn blast_spouge_etos(
    evalue: f64,
    kbp: &KarlinParams,
    gbp: &BlastGumbelBlk,
    query_length: i32,
    subject_length: i32,
) -> i32 {
    let mut a = 0i32;
    let db_scale_factor = if gbp.db_length != 0 {
        gbp.db_length as f64
    } else {
        1.0
    };
    let e0 = evalue.max(f64::MIN_POSITIVE);
    let mut b = ((db_scale_factor / e0).ln() / kbp.lambda) as i32;
    if b < 2 {
        b = 2;
    }

    let mut e = blast_spouge_stoe(b, kbp, gbp, query_length, subject_length);
    if e > e0 {
        while e > e0 {
            a = b;
            b *= 2;
            e = blast_spouge_stoe(b, kbp, gbp, query_length, subject_length);
        }
    } else {
        a = 0;
    }

    while b - a > 1 {
        let c = (a + b) / 2;
        e = blast_spouge_stoe(c, kbp, gbp, query_length, subject_length);
        if e > e0 {
            a = c;
        } else {
            b = c;
        }
    }

    a
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI c++/src/algo/blast/core/boost_erf.c:252-255:
    // double ErfC(double z) { return ErfImpl(z, TRUE); }
    // NCBI c++/src/algo/blast/core/blast_stat.c:5180-5231:
    // double scale_factor = kbp->Lambda / gbp->Lambda;
    // e_value = area * k_ * exp(-lambda_ * y_) * db_scale_factor;
    // Compare every captured preliminary and scaled post-redo StoE call.
    #[test]
    fn pinned_tblastn_spouge_calls_match_ncbi_double_bits() {
        let cases = [
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/multi_query_20260924_default.spouge.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/multi_query_20260924_control.spouge.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/seg_hard_query_20260924_default.spouge.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/seg_hard_query_20260924_control.spouge.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/run_20260923_default.spouge.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/parameters_20260924/run_20260923_control.spouge.tsv")),
        ];
        for case in cases {
            for row in case.lines() {
                let f: Vec<_> = row.split('\t').collect();
                let mut karlin =
                    crate::stats::tables::lookup_protein_params_gapped(ScoringMatrix::Blosum62);
                karlin.lambda = f[4].parse().unwrap();
                karlin.k = f[5].parse().unwrap();
                let gumbel = BlastGumbelBlk {
                    lambda: f[7].parse().unwrap(),
                    c: f[8].parse().unwrap(),
                    g: f[9].parse().unwrap(),
                    a: f[10].parse().unwrap(),
                    alpha: f[11].parse().unwrap(),
                    sigma: f[12].parse().unwrap(),
                    a_un: f[13].parse().unwrap(),
                    alpha_un: f[14].parse().unwrap(),
                    b: f[15].parse().unwrap(),
                    beta: f[16].parse().unwrap(),
                    tau: f[17].parse().unwrap(),
                    db_length: f[18].parse().unwrap(),
                };
                let actual = blast_spouge_stoe(
                    f[1].parse().unwrap(),
                    &karlin,
                    &gumbel,
                    f[2].parse().unwrap(),
                    f[3].parse().unwrap(),
                );
                let expected: f64 = f[20].parse().unwrap();
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "{row} actual={actual:.17e}"
                );
            }
        }
    }

    #[test]
    fn pinned_tblastn_erfc_calls_match_ncbi_double_bits() {
        let cases = [
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/erfc_20260924/multi_query_20260924_default.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/erfc_20260924/multi_query_20260924_control.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/erfc_20260924/seg_hard_query_20260924_default.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/erfc_20260924/seg_hard_query_20260924_control.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/erfc_20260924/run_20260923_default.tsv")),
            include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../docs/evidence/tlosan_stage_d/erfc_20260924/run_20260923_control.tsv")),
        ];
        for case in cases {
            for row in case.lines() {
                let fields: Vec<_> = row.split('\t').collect();
                let z: f64 = fields[1].parse().unwrap();
                let expected: f64 = fields[2].parse().unwrap();
                let actual = erfc_ncbi(z);
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "ErfC({z:.17e}) actual {actual:.17e} expected {expected:.17e}"
                );
            }
        }
    }

    #[test]
    fn test_lookup_protein_gumbel_params_blosum62_11_1() {
        let spec = ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        };
        let gbp = lookup_protein_gumbel_params(&spec, 66738).expect("BLOSUM62 11/1 must exist");
        assert!((gbp.lambda - 0.267).abs() < 1e-12);
        assert!((gbp.c - 0.669720).abs() < 1e-12);
        assert!((gbp.a - 1.9).abs() < 1e-12);
        assert!((gbp.alpha - 42.602800).abs() < 1e-12);
        assert!((gbp.sigma - 43.636200).abs() < 1e-12);
        assert_eq!(gbp.db_length, 66738);
    }

    #[test]
    fn test_blast_spouge_stoe_matches_ncbi_formula_for_scaled_lambda() {
        let gbp = BlastGumbelBlk {
            lambda: 0.267,
            c: 0.669720,
            g: 12.0,
            a: 1.9,
            alpha: 42.602800,
            sigma: 43.636200,
            a_un: 0.7916,
            alpha_un: 4.964660,
            b: 2.0 * 12.0 * (0.7916 - 1.9),
            beta: 2.0 * 12.0 * (4.964660 - 42.602800),
            tau: 2.0 * 12.0 * (4.964660 - 43.636200),
            db_length: 5000,
        };
        let kbp = KarlinParams {
            lambda: 0.267 / 32.0,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };

        let evalue = blast_spouge_stoe(1062, &kbp, &gbp, 100, 200);
        assert!((evalue - 1.5863473416482987).abs() < 1e-12);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:5176-5231
    // ```c
    // P_m_F = ErfC(-m_F / sqrt(2.0)) / 2.0;
    // ...
    // e_value = area * k_ * exp(-lambda_ * y_) * db_scale_factor;
    // ```
    #[test]
    fn test_blast_spouge_stoe_matches_ncbi_formula_for_representative_tail_inputs() {
        let gbp = BlastGumbelBlk {
            lambda: 0.267,
            c: 0.669720,
            g: 12.0,
            a: 1.9,
            alpha: 42.602800,
            sigma: 43.636200,
            a_un: 0.7916,
            alpha_un: 4.964660,
            b: 2.0 * 12.0 * (0.7916 - 1.9),
            beta: 2.0 * 12.0 * (4.964660 - 42.602800),
            tau: 2.0 * 12.0 * (4.964660 - 43.636200),
            db_length: 5000,
        };

        let moderate_kbp = KarlinParams {
            lambda: 0.267 / 16.0,
            k: 0.082,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        let moderate = blast_spouge_stoe(128, &moderate_kbp, &gbp, 64, 96);
        assert!((moderate - 4250.4559852614275).abs() < 1e-9);

        let long_tail_kbp = KarlinParams {
            lambda: 0.267 / 64.0,
            k: 0.021,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        let long_tail = blast_spouge_stoe(2048, &long_tail_kbp, &gbp, 250, 400);
        assert!((long_tail - 4.058690669998566).abs() < 1e-12);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:5236-5263
    // ```c
    // b = MAX((int)(log(db_scale_factor/e0) / kbp->Lambda), 2);
    // e = BLAST_SpougeStoE(b, kbp, gbp, m, n);
    // ...
    // return a;
    // ```
    #[test]
    fn test_blast_spouge_etos_matches_ncbi_formula_for_representative_evalues() {
        let gbp = BlastGumbelBlk {
            lambda: 0.267,
            c: 0.669720,
            g: 12.0,
            a: 1.9,
            alpha: 42.602800,
            sigma: 43.636200,
            a_un: 0.7916,
            alpha_un: 4.964660,
            b: 2.0 * 12.0 * (0.7916 - 1.9),
            beta: 2.0 * 12.0 * (4.964660 - 42.602800),
            tau: 2.0 * 12.0 * (4.964660 - 43.636200),
            db_length: 5000,
        };

        let default_kbp = KarlinParams {
            lambda: 0.267 / 32.0,
            k: 0.041,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        assert_eq!(blast_spouge_etos(1e-6, &default_kbp, &gbp, 100, 200), 2537);

        let tight_kbp = KarlinParams {
            lambda: 0.267 / 16.0,
            k: 0.082,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        };
        assert_eq!(blast_spouge_etos(1e-30, &tight_kbp, &gbp, 64, 96), 3749);
    }
}
