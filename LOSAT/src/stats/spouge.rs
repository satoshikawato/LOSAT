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

#[link(name = "m")]
unsafe extern "C" {
    fn erfc(x: f64) -> f64;
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

#[inline]
fn erfc_ncbi(x: f64) -> f64 {
    // SAFETY: `erfc` is a pure libm function with no aliasing or lifetime
    // requirements. The input is a plain finite or infinite `f64`, which
    // matches the C ABI expected by NCBI's `ErfC` path.
    unsafe { erfc(x) }
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
}
