//! PHI-BLAST's statistics: the Karlin parameters of its table (by matrix and gap costs), the
//! search space, the effective count of the query's occurrences, the score cutoff, the e-value and
//! the bit score. NCBI computes them in `double` in the order written here.

use super::scan::Occurrence;
use crate::utils::glibc_fma_libm as libm;

/// The Karlin block of a PHI search (`s_PHIScoreBlkFill`): `H` is 1, `logK` is not set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PhiKarlin {
    pub(crate) lambda: f64,
    pub(crate) k: f64,
    pub(crate) param_c: f64,
    pub(crate) h: f64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:165-300
/// The table: matrix, `paramC`, then (gap open, gap extend, Lambda, K).
const PHI_TABLE: &[(&str, f64, &[(i32, i32, f64, f64)])] = &[
    (
        "BLOSUM62",
        0.50,
        &[
            (11, 1, 0.270, 0.047),
            (9, 2, 0.285, 0.075),
            (8, 2, 0.265, 0.046),
            (7, 2, 0.243, 0.032),
            (12, 1, 0.281, 0.057),
            (10, 1, 0.250, 0.033),
        ],
    ),
    (
        "PAM30",
        0.30,
        &[
            (9, 1, 0.295, 0.13),
            (7, 2, 0.306, 0.15),
            (6, 2, 0.292, 0.13),
            (5, 2, 0.263, 0.077),
            (10, 1, 0.309, 0.15),
            (8, 1, 0.270, 0.070),
        ],
    ),
    (
        "PAM70",
        0.35,
        &[
            (10, 1, 0.291, 0.089),
            (8, 2, 0.303, 0.13),
            (7, 2, 0.287, 0.095),
            (6, 2, 0.269, 0.079),
            (11, 1, 0.307, 0.13),
            (9, 1, 0.269, 0.058),
        ],
    ),
    (
        "BLOSUM80",
        0.40,
        &[
            (10, 1, 0.300, 0.072),
            (8, 2, 0.308, 0.089),
            (7, 2, 0.295, 0.077),
            (6, 2, 0.271, 0.051),
            (11, 1, 0.314, 0.096),
            (9, 1, 0.277, 0.046),
        ],
    ),
    (
        "BLOSUM45",
        0.60,
        &[
            (14, 2, 0.199, 0.040),
            (13, 3, 0.209, 0.057),
            (12, 3, 0.203, 0.049),
            (11, 3, 0.193, 0.037),
            (10, 3, 0.182, 0.029),
            (15, 2, 0.206, 0.049),
            (13, 2, 0.190, 0.032),
            (12, 2, 0.177, 0.023),
            (19, 1, 0.209, 0.049),
            (18, 1, 0.202, 0.041),
            (17, 1, 0.195, 0.034),
            (16, 1, 0.183, 0.024),
        ],
    ),
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:137-310
/// The table's block for `matrix` (upper case, as NCBI's options keep it) and the gap costs, or
/// NCBI's warning (with its final LF) when the table has none: the search then finds nothing.
pub(crate) fn phi_karlin(
    matrix: &str,
    gap_open: i32,
    gap_extend: i32,
) -> Result<PhiKarlin, String> {
    let Some((_, param_c, rows)) = PHI_TABLE.iter().find(|(name, _, _)| *name == matrix) else {
        return Err(format!("Matrix {matrix} not allowed in PHI-BLAST\n"));
    };
    match rows
        .iter()
        .find(|(open, extend, _, _)| *open == gap_open && *extend == gap_extend)
    {
        Some(&(_, _, lambda, k)) => Ok(PhiKarlin {
            lambda,
            k,
            param_c: *param_c,
            h: 1.0,
        }),
        None => Err(format!(
            "The combination {gap_open} for gap opening cost and {gap_extend} for gap extension is not supported in PHI-BLAST with matrix {matrix}\n"
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:752-761
/// The search space of every context: the database's length less the minimum pattern length per
/// sequence (NCBI's `Int4` product widened).
pub(crate) fn search_space(db_length: i64, db_num_seqs: i64, min_match_length: i64) -> i64 {
    db_length - db_num_seqs * min_match_length
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:360-388
/// `PhiBlastGetEffectiveNumberOfPatterns`: the occurrences that start more than half the minimum
/// pattern length after the last one counted.
pub(crate) fn effective_occurrences(occurrences: &[Occurrence], min_match_length: i64) -> i32 {
    if occurrences.len() <= 1 {
        return occurrences.len() as i32;
    }
    let mut count = 1;
    let mut last = occurrences[0].offset;
    for occurrence in &occurrences[1..] {
        if (i64::from(occurrence.offset) - i64::from(last)) * 2 > min_match_length {
            last = occurrence.offset;
            count += 1;
        }
    }
    count
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:659-711
/// `s_PhiBlastCutoffScore` for the threshold `ethresh` (5 times `-evalue`): the bisection of
/// `[1, 100]` in at most 20 steps; its lower end.
pub(crate) fn cutoff_score(
    ethresh: f64,
    search_space: i64,
    karlin: &PhiKarlin,
    effective_occurrences: i32,
    probability: f64,
) -> i32 {
    let expect = |score: i32| {
        search_space as f64
            * karlin.param_c
            * (1.0 + karlin.lambda * f64::from(score))
            * f64::from(effective_occurrences)
            * probability
            * libm::exp(-karlin.lambda * f64::from(score))
    };
    let mut low = 1;
    let mut high = 100;
    for _ in 0..20 {
        let target = (low + high) / 2;
        if expect(target) > ethresh {
            low = target;
        } else {
            high = target;
        }
        if high - low <= 1 {
            break;
        }
    }
    low
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:399-415
/// `s_HSPPHIGetEvalue`: no search space; the count of the pattern's hits in the database instead.
pub(crate) fn evalue(
    score: i32,
    karlin: &PhiKarlin,
    effective_occurrences: i32,
    num_patterns_db: i32,
) -> f64 {
    karlin.param_c
        * (1.0 + karlin.lambda * f64::from(score))
        * f64::from(effective_occurrences)
        * f64::from(num_patterns_db)
        * libm::exp(-karlin.lambda * f64::from(score))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1934-1953
/// `Blast_HSPListPHIGetBitScores`.
pub(crate) fn bit_score(score: i32, karlin: &PhiKarlin) -> f64 {
    let score_lambda = f64::from(score) * karlin.lambda;
    (score_lambda - libm::ln(karlin.param_c) - libm::ln(1.0 + score_lambda))
        / std::f64::consts::LN_2
}

#[cfg(test)]
#[path = "stats_tests.rs"]
mod tests;
