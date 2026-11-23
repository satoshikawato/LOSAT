use super::karlin_calc::blast_gcd;
#[cfg(test)]
use super::protein_options::protein_gap_existence_extend_params;
use super::protein_options::{gap_costs_row, scoring_matrix_table};
use super::protein_tables::ProteinTableRow;
use crate::config::{NuclScoringSpec, ProteinScoringSpec, ScoringMatrix};

/// Karlin-Altschul statistical parameters
#[derive(Debug, Clone, Copy)]
pub struct KarlinParams {
    /// Lambda parameter for bit score calculation
    pub lambda: f64,
    /// K parameter for E-value calculation
    pub k: f64,
    /// H parameter (entropy) for length adjustment
    pub h: f64,
    /// Alpha parameter for length correction mean
    pub alpha: f64,
    /// Beta parameter for length correction
    pub beta: f64,
}

impl Default for KarlinParams {
    fn default() -> Self {
        Self {
            lambda: 0.625,
            k: 0.041,
            h: 0.85,
            alpha: 1.5,
            beta: -2.0,
        }
    }
}

/// Entry in the statistical parameter table
/// Format: (gap_open, gap_extend, lambda, k, h, alpha, beta)
#[derive(Debug, Clone, Copy)]
struct ParamEntry {
    gap_open: i32,
    gap_extend: i32,
    lambda: f64,
    k: f64,
    h: f64,
    alpha: f64,
    beta: f64,
}

impl ParamEntry {
    const fn new(
        gap_open: i32,
        gap_extend: i32,
        lambda: f64,
        k: f64,
        h: f64,
        alpha: f64,
        beta: f64,
    ) -> Self {
        Self {
            gap_open,
            gap_extend,
            lambda,
            k,
            h,
            alpha,
            beta,
        }
    }

    fn to_karlin_params(self) -> KarlinParams {
        KarlinParams {
            lambda: self.lambda,
            k: self.k,
            h: self.h,
            alpha: self.alpha,
            beta: self.beta,
        }
    }
}

// ============================================================================
// NUCLEOTIDE STATISTICAL PARAMETERS (from NCBI blast_stat.c)
// ============================================================================
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:611-614
// Format: { gap_open, gap_extend, lambda, k, h, alpha, beta, ... }
// LOSAT uses: ParamEntry::new(gap_open, gap_extend, lambda, k, h, alpha, beta)

/// Parameters for reward=1, penalty=-5
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:611-614 (blastn_values_1_5)
const BLASTN_1_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.39, 0.747, 1.38, 1.00, 0.0),
    ParamEntry::new(3, 3, 1.39, 0.747, 1.38, 1.00, 0.0),
];

/// Parameters for reward=1, penalty=-4
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:617-623 (blastn_values_1_4)
const BLASTN_1_4: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.383, 0.738, 1.36, 1.02, 0.0),
    ParamEntry::new(1, 2, 1.36, 0.67, 1.2, 1.1, 0.0),
    ParamEntry::new(0, 2, 1.26, 0.43, 0.90, 1.4, -1.0),
    ParamEntry::new(2, 1, 1.35, 0.61, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.22, 0.35, 0.72, 1.7, -3.0),
];

/// Parameters for reward=2, penalty=-7 (even scores only)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:629-635 (blastn_values_2_7)
/// Only even scores can use these parameters; an odd score is first lowered to the even value
/// just below it, then the e-value is computed.
const BLASTN_2_7: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.69, 0.73, 1.34, 0.515, 0.0),
    ParamEntry::new(2, 4, 0.68, 0.67, 1.2, 0.55, 0.0),
    ParamEntry::new(0, 4, 0.63, 0.43, 0.90, 0.7, -1.0),
    ParamEntry::new(4, 2, 0.675, 0.62, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 2, 0.61, 0.35, 0.72, 1.7, -3.0),
];

/// Parameters for reward=1, penalty=-3
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:638-645 (blastn_values_1_3)
const BLASTN_1_3: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.374, 0.711, 1.31, 1.05, 0.0),
    ParamEntry::new(2, 2, 1.37, 0.70, 1.2, 1.1, 0.0),
    ParamEntry::new(1, 2, 1.35, 0.64, 1.1, 1.2, -1.0),
    ParamEntry::new(0, 2, 1.25, 0.42, 0.83, 1.5, -2.0),
    ParamEntry::new(2, 1, 1.34, 0.60, 1.1, 1.2, -1.0),
    ParamEntry::new(1, 1, 1.21, 0.34, 0.71, 1.7, -2.0),
];

/// Parameters for reward=2, penalty=-5 (even scores only)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:651-657 (blastn_values_2_5)
/// Only even scores can use these parameters; an odd score is first lowered to the even value
/// just below it, then the e-value is computed.
const BLASTN_2_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.675, 0.65, 1.1, 0.6, -1.0),
    ParamEntry::new(2, 4, 0.67, 0.59, 1.1, 0.6, -1.0),
    ParamEntry::new(0, 4, 0.62, 0.39, 0.78, 0.8, -2.0),
    ParamEntry::new(4, 2, 0.67, 0.61, 1.0, 0.65, -2.0),
    ParamEntry::new(2, 2, 0.56, 0.32, 0.59, 0.95, -4.0),
];

/// Parameters for reward=1, penalty=-2 (megablast task default)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:660-668 (blastn_values_1_2)
const BLASTN_1_2: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 1.28, 0.46, 0.85, 1.5, -2.0),
    ParamEntry::new(2, 2, 1.33, 0.62, 1.1, 1.2, 0.0),
    ParamEntry::new(1, 2, 1.30, 0.52, 0.93, 1.4, -2.0),
    ParamEntry::new(0, 2, 1.19, 0.34, 0.66, 1.8, -3.0),
    ParamEntry::new(3, 1, 1.32, 0.57, 1.0, 1.3, -1.0),
    ParamEntry::new(2, 1, 1.29, 0.49, 0.92, 1.4, -1.0),
    ParamEntry::new(1, 1, 1.14, 0.26, 0.52, 2.2, -5.0),
];

/// Parameters for reward=2, penalty=-3 (blastn task default)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:674-684 (blastn_values_2_3)
/// Only even scores can use these parameters; an odd score is first lowered to the even value
/// just below it, then the e-value is computed.
/// NCBI blastn task default: gap_open=5, gap_extend=2 → ParamEntry::new(5, 2, 0.625, 0.41, 0.78, 0.8, -2.0)
const BLASTN_2_3: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.55, 0.21, 0.46, 1.2, -5.0),
    ParamEntry::new(4, 4, 0.63, 0.42, 0.84, 0.75, -2.0),
    ParamEntry::new(2, 4, 0.615, 0.37, 0.72, 0.85, -3.0),
    ParamEntry::new(0, 4, 0.55, 0.21, 0.46, 1.2, -5.0),
    ParamEntry::new(3, 3, 0.615, 0.37, 0.68, 0.9, -3.0),
    ParamEntry::new(6, 2, 0.63, 0.42, 0.84, 0.75, -2.0),
    ParamEntry::new(5, 2, 0.625, 0.41, 0.78, 0.8, -2.0),
    ParamEntry::new(4, 2, 0.61, 0.35, 0.68, 0.9, -3.0),
    ParamEntry::new(2, 2, 0.515, 0.14, 0.33, 1.55, -9.0),
];

/// Parameters for reward=3, penalty=-4
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:687-694 (blastn_values_3_4)
/// Only even scores can use these parameters; an odd score is first lowered to the even value
/// just below it, then the e-value is computed.
const BLASTN_3_4: &[ParamEntry] = &[
    ParamEntry::new(6, 3, 0.389, 0.25, 0.56, 0.7, -5.0),
    ParamEntry::new(5, 3, 0.375, 0.21, 0.47, 0.8, -6.0),
    ParamEntry::new(4, 3, 0.351, 0.14, 0.35, 1.0, -9.0),
    ParamEntry::new(6, 2, 0.362, 0.16, 0.45, 0.8, -4.0),
    ParamEntry::new(5, 2, 0.330, 0.092, 0.28, 1.2, -13.0),
    ParamEntry::new(4, 2, 0.281, 0.046, 0.16, 1.8, -23.0),
];

/// Parameters for reward=4, penalty=-5
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:697-703 (blastn_values_4_5)
const BLASTN_4_5: &[ParamEntry] = &[
    ParamEntry::new(0, 0, 0.22, 0.061, 0.22, 1.0, -15.0),
    ParamEntry::new(6, 5, 0.28, 0.21, 0.47, 0.6, -7.0),
    ParamEntry::new(5, 5, 0.27, 0.17, 0.39, 0.7, -9.0),
    ParamEntry::new(4, 5, 0.25, 0.10, 0.31, 0.8, -10.0),
    ParamEntry::new(3, 5, 0.23, 0.065, 0.25, 0.9, -11.0),
];

/// Parameters for reward=1, penalty=-1
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:706-714 (blastn_values_1_1)
const BLASTN_1_1: &[ParamEntry] = &[
    ParamEntry::new(3, 2, 1.09, 0.31, 0.55, 2.0, -2.0),
    ParamEntry::new(2, 2, 1.07, 0.27, 0.49, 2.2, -3.0),
    ParamEntry::new(1, 2, 1.02, 0.21, 0.36, 2.8, -6.0),
    ParamEntry::new(0, 2, 0.80, 0.064, 0.17, 4.8, -16.0),
    ParamEntry::new(4, 1, 1.08, 0.28, 0.54, 2.0, -2.0),
    ParamEntry::new(3, 1, 1.06, 0.25, 0.46, 2.3, -4.0),
    ParamEntry::new(2, 1, 0.99, 0.17, 0.30, 3.3, -10.0),
];

/// Parameters for reward=3, penalty=-2
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:717-719 (blastn_values_3_2)
const BLASTN_3_2: &[ParamEntry] = &[ParamEntry::new(5, 5, 0.208, 0.030, 0.072, 2.9, -47.0)];

/// Parameters for reward=5, penalty=-4
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:722-724 (blastn_values_5_4)
const BLASTN_5_4: &[ParamEntry] = &[
    ParamEntry::new(10, 6, 0.163, 0.068, 0.16, 1.0, -19.0),
    ParamEntry::new(8, 6, 0.146, 0.039, 0.11, 1.3, -29.0),
];

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3955-3963
pub(crate) fn ungapped_beta(reward: i32, penalty: i32) -> f64 {
    if (reward, penalty) == (1, -1) || (reward, penalty) == (2, -3) {
        -2.0
    } else {
        0.0
    }
}

/// The Karlin-Altschul values of a BLASTN reward/penalty pair (NCBI
/// `s_GetNuclValuesArray`): the table of the pair divided by its greatest common divisor,
/// with the gap costs multiplied and Lambda and alpha divided by the divisor.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3237-3372
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3159-3170 (s_SplitArrayOf8)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3192-3215 (s_AdjustGapParametersByGcd)
/// The pair is the engine's (1/-3 under matrix-only scoring, never a reward of 0), so the divisor
/// is at least 1.
pub struct NuclValues {
    /// The affine gap costs of the table, with their values.
    normal: Vec<ParamEntry>,
    /// The values of gap costs 0 and 0 (megablast's linear gaps), when the table has them.
    linear: Option<ParamEntry>,
    gap_open_max: i32,
    gap_extend_max: i32,
    /// Scores are rounded down to even numbers for e-values and bit scores.
    pub round_down: bool,
}

impl NuclValues {
    pub fn new(reward: i32, penalty: i32) -> Result<Self, String> {
        let divisor = blast_gcd(reward, penalty).max(1);
        let (reward, penalty) = (reward / divisor, penalty / divisor);
        let (table, gap_open_max, gap_extend_max, round_down) = match (reward, penalty) {
            (1, -5) => (BLASTN_1_5, 3, 3, false),
            (1, -4) => (BLASTN_1_4, 2, 2, false),
            (2, -7) => (BLASTN_2_7, 4, 4, true),
            (1, -3) => (BLASTN_1_3, 2, 2, false),
            (2, -5) => (BLASTN_2_5, 4, 4, true),
            (1, -2) => (BLASTN_1_2, 2, 2, false),
            (2, -3) => (BLASTN_2_3, 6, 4, true),
            (3, -4) => (BLASTN_3_4, 6, 3, true),
            (1, -1) => (BLASTN_1_1, 4, 2, false),
            (3, -2) => (BLASTN_3_2, 5, 5, false),
            (4, -5) => (BLASTN_4_5, 12, 8, false),
            (5, -4) => (BLASTN_5_4, 25, 10, false),
            _ => {
                return Err(format!(
                    "Substitution scores {reward} and {penalty} are not supported"
                ))
            }
        };
        let adjust = |entry: &ParamEntry| ParamEntry {
            gap_open: entry.gap_open * divisor,
            gap_extend: entry.gap_extend * divisor,
            lambda: entry.lambda / divisor as f64,
            alpha: entry.alpha / divisor as f64,
            ..*entry
        };
        let split = table[0].gap_open == 0 && table[0].gap_extend == 0;
        Ok(Self {
            normal: table[usize::from(split)..].iter().map(adjust).collect(),
            linear: split.then(|| adjust(&table[0])),
            gap_open_max: gap_open_max * divisor,
            gap_extend_max: gap_extend_max * divisor,
            round_down,
        })
    }

    /// The table row of the gap costs: `Some` row, `None` for gap costs beyond the table
    /// (which use the ungapped block), or NCBI's message for unsupported gap costs
    /// (`Blast_KarlinBlkNuclGappedCalc`). The message names the given (not the divided)
    /// scores.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3875-3935
    fn row(&self, spec: &NuclScoringSpec) -> Result<Option<ParamEntry>, String> {
        let row = if spec.gap_open == 0 && spec.gap_extend == 0 && self.linear.is_some() {
            self.linear
        } else {
            self.normal
                .iter()
                .find(|entry| {
                    entry.gap_open == spec.gap_open && entry.gap_extend == spec.gap_extend
                })
                .copied()
        };
        if row.is_some()
            || (spec.gap_open >= self.gap_open_max && spec.gap_extend >= self.gap_extend_max)
        {
            return Ok(row);
        }
        let mut message = format!(
            "Gap existence and extension values {} and {} are not supported for substitution scores {} and {}\n",
            spec.gap_open, spec.gap_extend, spec.reward, spec.penalty
        );
        let supported = self
            .normal
            .iter()
            .map(|entry| (entry.gap_open, entry.gap_extend))
            .chain([(self.gap_open_max, self.gap_extend_max)]);
        for (gap_open, gap_extend) in supported {
            message += &format!(
                "{gap_open} and {gap_extend} are supported existence and extension values\n"
            );
        }
        message += &format!(
            "Any values more stringent than {} and {} are supported\n",
            self.gap_open_max, self.gap_extend_max
        );
        Err(message)
    }

    /// NCBI's message when the gap costs are not supported for these scores.
    pub fn check_gaps(&self, spec: &NuclScoringSpec) -> Result<(), String> {
        self.row(spec).map(|_| ())
    }

    /// The gapped Karlin block of a query context, with the alpha and beta of its length
    /// adjustment (NCBI `Blast_KarlinBlkNuclGappedCalc` and `Blast_GetNuclAlphaBeta`).
    /// `ungapped` is the context's ungapped block, which gap costs beyond the table use.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3995-4026
    pub fn gapped(
        &self,
        spec: &NuclScoringSpec,
        ungapped: &KarlinParams,
    ) -> Result<KarlinParams, String> {
        Ok(match self.row(spec)? {
            Some(row) => row.to_karlin_params(),
            None => KarlinParams {
                alpha: ungapped.lambda / ungapped.h,
                beta: ungapped_beta(spec.reward, spec.penalty),
                ..*ungapped
            },
        })
    }
}

/// Lambda, K, H, alpha and beta of a protein table row (NCBI's columns 3-7).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3617-3623, 3078-3081
fn protein_row_params(row: &ProteinTableRow) -> KarlinParams {
    KarlinParams {
        lambda: row.0[3],
        k: row.0[4],
        h: row.0[5],
        alpha: row.0[6],
        beta: row.0[7],
    }
}

/// The gapped Karlin-Altschul parameters of a protein scoring scheme: the table row of its
/// gap costs (`Blast_KarlinBlkGappedLoadFromTables` for Lambda, K and H, `BLAST_GetAlphaBeta`
/// for alpha and beta).
///
/// A gapped search checks the gap costs against the table before it starts
/// (`BLAST_ValidateOptions`,
/// `blastinput::protein_options_check::validate_protein_options`; blastx accepts 11/1 only), so
/// costs outside the table are a LOSAT defect here: it panics instead of reading
/// another row. NCBI does not check them for an ungapped search, which no LOSAT protein program
/// runs yet; a change that opens `-ungapped` must not call this with unchecked costs.
///
/// `BLAST_GetAlphaBeta` searches from row 1, so the ungapped row's costs (32767/32767, which
/// the validation accepts) get Lambda, K and H from row 0 and the alpha and beta that
/// `BLAST_CalcEffLengths` starts with, 0 and 0.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3117-3126
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:705
pub fn lookup_protein_params(spec: &ProteinScoringSpec) -> KarlinParams {
    let table = scoring_matrix_table(spec.matrix);
    let row = gap_costs_row(table, spec.gap_open, spec.gap_extend).unwrap_or_else(|| {
        panic!(
            "gap costs {}/{} are not in NCBI's {} table; they are validated before the search",
            spec.gap_open, spec.gap_extend, table.name
        )
    });
    let mut params = protein_row_params(row);
    // Gapped rows never carry the ungapped row's costs, so only row 0 misses the search from
    // row 1.
    if std::ptr::eq(row, &table.values[0]) {
        params.alpha = 0.0;
        params.beta = 0.0;
    }
    params
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:908-936
pub fn protein_scoring_supported(spec: &ProteinScoringSpec) -> bool {
    gap_costs_row(
        scoring_matrix_table(spec.matrix),
        spec.gap_open,
        spec.gap_extend,
    )
    .is_some()
}

/// Look up UNGAPPED Karlin-Altschul parameters for protein scoring scheme: the first row of
/// the matrix's table (gap costs INT2_MAX).
///
/// NCBI BLAST uses ungapped params (kbp_std) for gap_trigger calculation.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:340-345
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3128-3130
pub fn lookup_protein_params_ungapped(matrix: ScoringMatrix) -> KarlinParams {
    protein_row_params(&scoring_matrix_table(matrix).values[0])
}

/// Look up GAPPED Karlin-Altschul parameters for protein scoring scheme: the row of the
/// matrix's preferred gap costs (BLOSUM62 11/1).
///
/// NCBI BLAST uses gapped params (kbp_gap_std with gap_open=11, gap_extend=1 for BLOSUM62)
/// for length_adjustment calculation in BLAST_CalcEffLengths.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:814-819
///
/// For tblastx, `gapped_calculation` in `scoring_options` defaults to true,
/// so BLAST_GetAlphaBeta returns gapped alpha/beta values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3107-3116
#[cfg(test)]
pub fn lookup_protein_params_gapped(matrix: ScoringMatrix) -> KarlinParams {
    let (gap_open, gap_extend) = protein_gap_existence_extend_params(matrix.ncbi_name())
        .expect("every matrix has a BEST row");
    lookup_protein_params(&ProteinScoringSpec {
        matrix,
        gap_open,
        gap_extend,
    })
}

#[cfg(test)]
#[path = "tables_tests.rs"]
mod tests;
