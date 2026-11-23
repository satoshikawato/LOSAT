//! NCBI composition-adjustment helpers used by `blastp` redo post-processing.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_adjustment.h
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_mode_condition.c
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/unified_pvalues.c
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/matrix_frequency_data.c

use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::{bail, Result};
use std::sync::OnceLock;

use crate::config::ScoringMatrix;
use crate::utils::matrix::{protein_score_direct, BLASTAA_SIZE};

use super::matrix_tables::{self, CompoMatrixData, COMPO_MATRIX_DATA};
use super::redo_alignment::{BlastCompoAdjustMode, EMatrixAdjustRule};

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_constants.h:39-51
const COMPO_SCORE_MIN: i32 = i16::MIN as i32;
pub(crate) const COMPO_NUM_TRUE_AA: usize = 20;
const COMPO_LARGEST_ALPHABET: usize = BLASTAA_SIZE;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:77-80
#[allow(dead_code)]
const K_FIXED_RE_BLOSUM62: f64 = 0.44;
const K_MAXIMUM_X_SCORE: f64 = -1.0;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:64-75
const K_LAMBDA_ITERATION_LIMIT: usize = 100;
const K_LAMBDA_ERROR_TOLERANCE: f64 = 1.0e-7;
const K_LAMBDA_FUNCTION_TOLERANCE: f64 = 1.0e-5;
pub(crate) const K_COMPO_ADJUST_ERR_TOLERANCE: f64 = 1.0e-8;
pub(crate) const K_COMPO_ADJUST_ITERATION_LIMIT: usize = 2000;

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_adjustment.h:42-49
const E_BCHAR: usize = 2;
const E_CCHAR: usize = 3;
const E_DCHAR: usize = 4;
const E_ECHAR: usize = 5;
const E_ICHAR: usize = 9;
const E_LCHAR: usize = 11;
const E_NCHAR: usize = 13;
const E_QCHAR: usize = 15;
const E_XCHAR: usize = 21;
const E_ZCHAR: usize = 23;
const E_SELENOCYSTEINE: usize = 24;
const E_STOP_CHAR: usize = 25;
const E_OCHAR: usize = 26;
const E_JCHAR: usize = 27;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1097-1100
const LAMBDA_RATIO_LOWER_BOUND: f64 = 0.5;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_mode_condition.c:38-47
const HALF_CIRCLE_DEGREES: f64 = 180.0;
// NCBI's digits, not `std::f64::consts::PI` (3.14159265358979...): the rule's angle must be NCBI's
// to the last bit at the 70-degree threshold.
const PI: f64 = 3.1415926543;
const QUERY_MATCH_DISTANCE_THRESHOLD: f64 = 0.16;
const LENGTH_RATIO_THRESHOLD: f64 = 3.0;
const ANGLE_DEGREE_THRESHOLD: f64 = 70.0;
const HIGH_PAIR_THRESHOLD: f64 = 0.4;
const LENGTH_LOWER_THRESHOLD: i32 = 50;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/unified_pvalues.c:37-53
const COMPO_MIN_LAMBDA: f64 = 0.034;
const LAMBDA_BIN_SIZE: f64 = 0.001;
const LAMBDA_BIN_TOTAL: usize = 565;
const LOW_LAMBDA_BIN_CUT: usize = 35;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:44-54
const TRUE_CHAR_POSITIONS: [usize; COMPO_NUM_TRUE_AA] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];
const ALPHA_CONVERT: [i32; COMPO_LARGEST_ALPHABET] = [
    -1, 0, -1, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, -1, 18, -1, -1, -1,
    -1, -1,
];

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/matrix_frequency_data.c:383-537
// The values are in the generated `matrix_tables.rs`.
const BLOSUM62_JOINT_PROBS: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA] =
    matrix_tables::BLOSUM62_JOINT_PROBS;

const BLOSUM62_BG: [f64; COMPO_NUM_TRUE_AA] = matrix_tables::BLOSUM62_BG;

// The composition-based statistics data of `matrix`: its frequency ratios and its joint and
// background probabilities.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/matrix_frequency_data.c:1359-1369
// NCBI reference (598d8ae6): c++/src/algo/blast/core/matrix_freq_ratios.c:1734-1736
// `_PSIMatrixFrequencyRatiosNew` has the frequency ratios of the same eight matrices, and of
// BLOSUM62_20, 20A and 20B, which the 2.17.0 build has no Karlin table for (`blast_stat.c:435`,
// `BLOSUM62_20_ENABLE`), so NCBI's option validation refuses them. IDENTITY has none of the data:
// NCBI turns composition-based statistics off for it before the search.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:144-150
pub(crate) fn compo_matrix_data(matrix: ScoringMatrix) -> Option<&'static CompoMatrixData> {
    let index = match matrix {
        ScoringMatrix::Blosum62 => 0,
        ScoringMatrix::Blosum80 => 1,
        ScoringMatrix::Blosum50 => 2,
        ScoringMatrix::Blosum90 => 3,
        ScoringMatrix::Blosum45 => 4,
        ScoringMatrix::Pam30 => 5,
        ScoringMatrix::Pam70 => 6,
        ScoringMatrix::Pam250 => 7,
        ScoringMatrix::Identity => return None,
    };
    Some(&COMPO_MATRIX_DATA[index])
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/unified_pvalues.c:52-205
const P_LAMBDA_TABLE: [f64; LAMBDA_BIN_TOTAL] = [
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    1.247750e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    2.495500e-07,
    3.743250e-07,
    4.991000e-07,
    4.991000e-07,
    4.991000e-07,
    4.991000e-07,
    7.486490e-07,
    7.486490e-07,
    7.486490e-07,
    7.486490e-07,
    8.734240e-07,
    9.981990e-07,
    9.981990e-07,
    1.122974e-06,
    1.122974e-06,
    1.372523e-06,
    1.497298e-06,
    1.622073e-06,
    1.871622e-06,
    1.871622e-06,
    1.996397e-06,
    2.370721e-06,
    2.620270e-06,
    2.620270e-06,
    2.620270e-06,
    2.745045e-06,
    2.745045e-06,
    2.745045e-06,
    2.745045e-06,
    2.869820e-06,
    3.493695e-06,
    3.493695e-06,
    3.493695e-06,
    3.868019e-06,
    4.117568e-06,
    4.367117e-06,
    4.866216e-06,
    4.866216e-06,
    5.240540e-06,
    5.240540e-06,
    5.864415e-06,
    6.363514e-06,
    6.737838e-06,
    7.236937e-06,
    7.860812e-06,
    8.110361e-06,
    8.484685e-06,
    9.233334e-06,
    9.857209e-06,
    1.035631e-05,
    1.060586e-05,
    1.085541e-05,
    1.210316e-05,
    1.272703e-05,
    1.322613e-05,
    1.409955e-05,
    1.509775e-05,
    1.647027e-05,
    1.784279e-05,
    1.859144e-05,
    1.921532e-05,
    2.008874e-05,
    2.133649e-05,
    2.245946e-05,
    2.333288e-05,
    2.445585e-05,
    2.557882e-05,
    2.732567e-05,
    2.832387e-05,
    3.019549e-05,
    3.231666e-05,
    3.468738e-05,
    3.568558e-05,
    3.718288e-05,
    3.942882e-05,
    4.292251e-05,
    4.479413e-05,
    4.666575e-05,
    4.903647e-05,
    5.140719e-05,
    5.352836e-05,
    5.539998e-05,
    5.714683e-05,
    5.964232e-05,
    6.263691e-05,
    6.550673e-05,
    6.825177e-05,
    7.124636e-05,
    7.411618e-05,
    7.611258e-05,
    7.910717e-05,
    8.235131e-05,
    8.534590e-05,
    8.971302e-05,
    9.457923e-05,
    9.907112e-05,
    1.028144e-04,
    1.068072e-04,
    1.124220e-04,
    1.150423e-04,
    1.202828e-04,
    1.241509e-04,
    1.297657e-04,
    1.352558e-04,
    1.408707e-04,
    1.474838e-04,
    1.535977e-04,
    1.605851e-04,
    1.679469e-04,
    1.766811e-04,
    1.849162e-04,
    1.943991e-04,
    2.023847e-04,
    2.111190e-04,
    2.208514e-04,
    2.298352e-04,
    2.390685e-04,
    2.500487e-04,
    2.599059e-04,
    2.705118e-04,
    2.826149e-04,
    2.953419e-04,
    3.058230e-04,
    3.209207e-04,
    3.336477e-04,
    3.469986e-04,
    3.607238e-04,
    3.758215e-04,
    3.915431e-04,
    4.100098e-04,
    4.283517e-04,
    4.454458e-04,
    4.661584e-04,
    4.845003e-04,
    5.032165e-04,
    5.194372e-04,
    5.399003e-04,
    5.642314e-04,
    5.905589e-04,
    6.173855e-04,
    6.407184e-04,
    6.646751e-04,
    6.896300e-04,
    7.199503e-04,
    7.487733e-04,
    7.780954e-04,
    8.119093e-04,
    8.442260e-04,
    8.782895e-04,
    9.142246e-04,
    9.507836e-04,
    9.890894e-04,
    1.025898e-03,
    1.060710e-03,
    1.102010e-03,
    1.141065e-03,
    1.179870e-03,
    1.225413e-03,
    1.268086e-03,
    1.313254e-03,
    1.362665e-03,
    1.413199e-03,
    1.465479e-03,
    1.514516e-03,
    1.568668e-03,
    1.620325e-03,
    1.678595e-03,
    1.738861e-03,
    1.800375e-03,
    1.863511e-03,
    1.933760e-03,
    1.999641e-03,
    2.066271e-03,
    2.136269e-03,
    2.208639e-03,
    2.282256e-03,
    2.363235e-03,
    2.451451e-03,
    2.535050e-03,
    2.626385e-03,
    2.715973e-03,
    2.811676e-03,
    2.906879e-03,
    3.008695e-03,
    3.108390e-03,
    3.220687e-03,
    3.337351e-03,
    3.458507e-03,
    3.587275e-03,
    3.717165e-03,
    3.853419e-03,
    3.998407e-03,
    4.150882e-03,
    4.307349e-03,
    4.470305e-03,
    4.650105e-03,
    4.835520e-03,
    5.028047e-03,
    5.239291e-03,
    5.459394e-03,
    5.699835e-03,
    5.944019e-03,
    6.201679e-03,
    6.472939e-03,
    6.767532e-03,
    7.067490e-03,
    7.402760e-03,
    7.779704e-03,
    8.163261e-03,
    8.581007e-03,
    9.043672e-03,
    9.538778e-03,
    1.006046e-02,
    1.061571e-02,
    1.121351e-02,
    1.187918e-02,
    1.259427e-02,
    1.336213e-02,
    1.420611e-02,
    1.510486e-02,
    1.609220e-02,
    1.716439e-02,
    1.831930e-02,
    1.961309e-02,
    2.100770e-02,
    2.254767e-02,
    2.423026e-02,
    2.604985e-02,
    2.803414e-02,
    3.023304e-02,
    3.261836e-02,
    3.527531e-02,
    3.818206e-02,
    4.134298e-02,
    4.483305e-02,
    4.870456e-02,
    5.296474e-02,
    5.765964e-02,
    6.281384e-02,
    6.850944e-02,
    7.481781e-02,
    8.173895e-02,
    8.937667e-02,
    9.774769e-02,
    1.069522e-01,
    1.170782e-01,
    1.282548e-01,
    1.405027e-01,
    1.539673e-01,
    1.686191e-01,
    1.846574e-01,
    2.020771e-01,
    2.210459e-01,
    2.415374e-01,
    2.637095e-01,
    2.872893e-01,
    3.124426e-01,
    3.391945e-01,
    3.673232e-01,
    3.966229e-01,
    4.269345e-01,
    4.580085e-01,
    4.895437e-01,
    5.211873e-01,
    5.525534e-01,
    5.832051e-01,
    6.129409e-01,
    6.414552e-01,
    6.685366e-01,
    6.941538e-01,
    7.180993e-01,
    7.403272e-01,
    7.608830e-01,
    7.798100e-01,
    7.972520e-01,
    8.132446e-01,
    8.277589e-01,
    8.410031e-01,
    8.531792e-01,
    8.642644e-01,
    8.743145e-01,
    8.835475e-01,
    8.919562e-01,
    8.996102e-01,
    9.066574e-01,
    9.130520e-01,
    9.189736e-01,
    9.244128e-01,
    9.293928e-01,
    9.339722e-01,
    9.381913e-01,
    9.421451e-01,
    9.457657e-01,
    9.491230e-01,
    9.522394e-01,
    9.551589e-01,
    9.578458e-01,
    9.603525e-01,
    9.626876e-01,
    9.648771e-01,
    9.669361e-01,
    9.688645e-01,
    9.706878e-01,
    9.723798e-01,
    9.739466e-01,
    9.754310e-01,
    9.768135e-01,
    9.781311e-01,
    9.793539e-01,
    9.805091e-01,
    9.816204e-01,
    9.826387e-01,
    9.835868e-01,
    9.844858e-01,
    9.853414e-01,
    9.861397e-01,
    9.869005e-01,
    9.876234e-01,
    9.882862e-01,
    9.889344e-01,
    9.895357e-01,
    9.901063e-01,
    9.906435e-01,
    9.911373e-01,
    9.916119e-01,
    9.920576e-01,
    9.924838e-01,
    9.928730e-01,
    9.932470e-01,
    9.935931e-01,
    9.939352e-01,
    9.942430e-01,
    9.945419e-01,
    9.948185e-01,
    9.950818e-01,
    9.953442e-01,
    9.955881e-01,
    9.958246e-01,
    9.960354e-01,
    9.962523e-01,
    9.964441e-01,
    9.966252e-01,
    9.967933e-01,
    9.969595e-01,
    9.971159e-01,
    9.972685e-01,
    9.974102e-01,
    9.975482e-01,
    9.976761e-01,
    9.977899e-01,
    9.978979e-01,
    9.979996e-01,
    9.980933e-01,
    9.981938e-01,
    9.982869e-01,
    9.983711e-01,
    9.984526e-01,
    9.985349e-01,
    9.986089e-01,
    9.986793e-01,
    9.987437e-01,
    9.988027e-01,
    9.988587e-01,
    9.989131e-01,
    9.989656e-01,
    9.990185e-01,
    9.990667e-01,
    9.991120e-01,
    9.991549e-01,
    9.991993e-01,
    9.992385e-01,
    9.992718e-01,
    9.993088e-01,
    9.993393e-01,
    9.993723e-01,
    9.993982e-01,
    9.994278e-01,
    9.994549e-01,
    9.994825e-01,
    9.995078e-01,
    9.995336e-01,
    9.995536e-01,
    9.995732e-01,
    9.995946e-01,
    9.996155e-01,
    9.996358e-01,
    9.996535e-01,
    9.996697e-01,
    9.996842e-01,
    9.996983e-01,
    9.997115e-01,
    9.997263e-01,
    9.997404e-01,
    9.997520e-01,
    9.997612e-01,
    9.997721e-01,
    9.997834e-01,
    9.997917e-01,
    9.998010e-01,
    9.998097e-01,
    9.998197e-01,
    9.998272e-01,
    9.998347e-01,
    9.998427e-01,
    9.998504e-01,
    9.998564e-01,
    9.998629e-01,
    9.998688e-01,
    9.998745e-01,
    9.998795e-01,
    9.998834e-01,
    9.998882e-01,
    9.998937e-01,
    9.998978e-01,
    9.999021e-01,
    9.999067e-01,
    9.999094e-01,
    9.999134e-01,
    9.999174e-01,
    9.999209e-01,
    9.999240e-01,
    9.999275e-01,
    9.999312e-01,
    9.999336e-01,
    9.999374e-01,
    9.999390e-01,
    9.999415e-01,
    9.999435e-01,
    9.999452e-01,
    9.999480e-01,
    9.999506e-01,
    9.999517e-01,
    9.999541e-01,
    9.999555e-01,
    9.999570e-01,
    9.999581e-01,
    9.999595e-01,
    9.999608e-01,
    9.999623e-01,
    9.999636e-01,
    9.999655e-01,
    9.999665e-01,
    9.999680e-01,
    9.999691e-01,
    9.999702e-01,
    9.999712e-01,
    9.999718e-01,
    9.999733e-01,
    9.999739e-01,
    9.999751e-01,
    9.999761e-01,
    9.999771e-01,
    9.999782e-01,
    9.999792e-01,
    9.999801e-01,
    9.999807e-01,
    9.999813e-01,
    9.999816e-01,
    9.999826e-01,
    9.999831e-01,
    9.999834e-01,
    9.999836e-01,
    9.999838e-01,
    9.999841e-01,
    9.999846e-01,
    9.999854e-01,
    9.999856e-01,
    9.999859e-01,
    9.999870e-01,
    9.999874e-01,
    9.999882e-01,
    9.999885e-01,
    9.999890e-01,
    9.999892e-01,
    9.999892e-01,
    9.999893e-01,
    9.999895e-01,
    9.999900e-01,
    9.999904e-01,
    9.999907e-01,
    9.999909e-01,
    9.999914e-01,
    9.999917e-01,
    9.999923e-01,
    9.999923e-01,
    9.999925e-01,
    9.999928e-01,
    9.999932e-01,
    9.999934e-01,
    9.999937e-01,
    9.999939e-01,
    9.999942e-01,
    9.999948e-01,
    9.999950e-01,
    9.999950e-01,
    9.999952e-01,
    9.999954e-01,
    9.999955e-01,
    9.999957e-01,
    9.999957e-01,
    9.999959e-01,
    9.999959e-01,
    9.999960e-01,
    9.999960e-01,
    9.999962e-01,
    9.999963e-01,
    9.999965e-01,
    9.999965e-01,
    9.999967e-01,
    9.999967e-01,
    9.999968e-01,
    9.999973e-01,
    9.999973e-01,
    9.999974e-01,
    9.999975e-01,
    9.999977e-01,
    9.999977e-01,
    9.999977e-01,
    9.999977e-01,
    9.999979e-01,
    9.999979e-01,
    9.999979e-01,
    9.999982e-01,
    9.999983e-01,
    9.999984e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999985e-01,
    9.999987e-01,
    9.999987e-01,
    9.999988e-01,
    9.999988e-01,
    9.999988e-01,
    9.999988e-01,
    9.999989e-01,
    9.999989e-01,
    9.999989e-01,
    9.999989e-01,
    9.999992e-01,
];

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_adjustment.h:49-59
#[derive(Debug, Clone, PartialEq)]
pub struct BlastAminoAcidComposition {
    pub prob: [f64; COMPO_LARGEST_ALPHABET],
    pub num_true_amino_acids: i32,
}

impl BlastAminoAcidComposition {
    pub fn empty() -> Self {
        Self {
            prob: [0.0; COMPO_LARGEST_ALPHABET],
            num_true_amino_acids: 0,
        }
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_adjustment.h:76-90
#[derive(Debug, Clone, PartialEq)]
pub struct BlastMatrixInfo {
    pub matrix: ScoringMatrix,
    pub position_based: bool,
    pub cols: usize,
    pub ungapped_lambda: f64,
    pub start_matrix: [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    pub start_freq_ratios: [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
}

// NCBI reference (598d8ae6): c++/include/algo/blast/composition_adjustment/composition_adjustment.h:109-121
#[derive(Debug, Clone)]
pub(crate) struct BlastCompositionWorkspace {
    mat_b: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    mat_final: [[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    first_standard_freq: [f64; COMPO_NUM_TRUE_AA],
    second_standard_freq: [f64; COMPO_NUM_TRUE_AA],
}

impl BlastCompositionWorkspace {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1271-1332
    // BLOSUM62's row and column sums of the joint probabilities equal its background
    // frequencies bit for bit, so this is `new(ScoringMatrix::Blosum62)` (test).
    pub(crate) fn new_blosum62() -> Self {
        Self {
            mat_b: BLOSUM62_JOINT_PROBS,
            mat_final: [[0.0; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
            first_standard_freq: BLOSUM62_BG,
            second_standard_freq: BLOSUM62_BG,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1316-1331
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/matrix_frequency_data.c:1406-1416
    // The standard frequencies are the sums, not the background frequencies: the two differ in
    // the last bits for PAM30, PAM70 and PAM250.
    pub(crate) fn new(matrix: ScoringMatrix) -> Result<Self> {
        let Some(data) = compo_matrix_data(matrix) else {
            bail!(
                "Matrix {} not currently supported for RE based adjustment",
                matrix.ncbi_name()
            );
        };
        let mut first_standard_freq = [0.0; COMPO_NUM_TRUE_AA];
        let mut second_standard_freq = [0.0; COMPO_NUM_TRUE_AA];
        for i in 0..COMPO_NUM_TRUE_AA {
            for j in 0..COMPO_NUM_TRUE_AA {
                first_standard_freq[i] += data.joint_probs[i][j];
                second_standard_freq[j] += data.joint_probs[i][j];
            }
        }
        Ok(Self {
            mat_b: *data.joint_probs,
            mat_final: [[0.0; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
            first_standard_freq,
            second_standard_freq,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1025-1151
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdjustedProteinMatrix {
    pub scores: [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
}

impl AdjustedProteinMatrix {
    #[inline]
    pub fn score(&self, query_residue: u8, subject_residue: u8) -> i32 {
        self.scores[query_residue as usize][subject_residue as usize]
    }

    #[inline]
    pub fn is_positive(&self, query_residue: u8, subject_residue: u8) -> bool {
        self.score(query_residue, subject_residue) > 0
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1414-1530
#[derive(Debug, Clone)]
pub struct BlastAdjustScoresResult {
    pub matrix_adjust_rule: EMatrixAdjustRule,
    pub adjusted_matrix: AdjustedProteinMatrix,
    pub pvalue_for_this_pair: Option<f64>,
    pub lambda_ratio: f64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:61,1204-1248:
// Expand the composition interval outward from the alignment endpoints to the nearest stop
// residue (`eStopChar`) on each side, stepping inward from each stop by `kCompositionMargin`
// (20 residues) without shrinking inside the alignment range.
pub(crate) fn blast_get_composition_range(
    subject_data: &[u8],
    start: usize,
    finish: usize,
) -> Result<(usize, usize)> {
    if start > finish || finish > subject_data.len() {
        anyhow::bail!("invalid NCBI translated composition interval");
    }
    let mut left = 0;
    for i in (1..=start).rev() {
        if subject_data[i - 1] == 25 {
            left = if i + 20 < start { i + 20 } else { start };
            break;
        }
    }
    let mut right = subject_data.len();
    for i in finish..subject_data.len() {
        if subject_data[i] == 25 {
            right = if i.saturating_sub(20) > finish {
                i - 20
            } else {
                finish
            };
            break;
        }
    }
    Ok((left, right))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1170-1194
pub fn read_aa_composition(sequence: &[u8]) -> BlastAminoAcidComposition {
    let mut composition = BlastAminoAcidComposition::empty();
    for &residue in sequence {
        let residue_index = residue as usize;
        if residue_index < COMPO_LARGEST_ALPHABET
            && (ALPHA_CONVERT[residue_index] >= 0 || residue == 24)
        {
            composition.prob[residue_index] += 1.0;
            composition.num_true_amino_acids += 1;
        }
    }

    if composition.prob[24] > 0.0 {
        composition.prob[3] += composition.prob[24];
        composition.prob[24] = 0.0;
    }

    if composition.num_true_amino_acids > 0 {
        let denom = composition.num_true_amino_acids as f64;
        for prob in &mut composition.prob {
            *prob /= denom;
        }
    }
    composition
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2216-2232
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:657-660
pub fn build_matrix_info(matrix: ScoringMatrix, ungapped_lambda: f64) -> Result<BlastMatrixInfo> {
    let Some(data) = compo_matrix_data(matrix) else {
        bail!(
            "composition adjustment has no frequency ratios for the {} matrix",
            matrix.ncbi_name()
        );
    };

    let start_freq_ratios = *data.freq_ratios;
    let mut start_matrix = [[0i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
    int4_matrix_from_freq(
        &mut start_matrix,
        COMPO_LARGEST_ALPHABET,
        &start_freq_ratios,
        ungapped_lambda,
    );

    Ok(BlastMatrixInfo {
        matrix,
        position_based: false,
        cols: COMPO_LARGEST_ALPHABET,
        ungapped_lambda,
        start_matrix,
        start_freq_ratios,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:31-48
fn dense_matrix_new(nrows: usize, ncols: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; ncols]; nrows]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:52-74
fn ltriang_matrix_new(n: usize) -> Vec<Vec<f64>> {
    (0..n).map(|row| vec![0.0; row + 1]).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:725-727
// The stride is `m` for `alphsize == 20` and `constrain_rel_entropy` (the size of `W`).
/// Column stride of the column-major copy of the Newton matrix.
const X_WC_STRIDE: usize = 2 * COMPO_NUM_TRUE_AA;

// `Nlm_FactorLtriangPosDef` with the loops interchanged.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:120-155
// Every element (i, j) still receives exactly the same floating-point
// operations in the same order: the subtractions for k = 0, 1, ..., j-1 (each
// one `temp -= A[i][k] * A[j][k]` with the final values of column k), then the
// division by A[j][j] (or the square root on the diagonal). Only the order in
// which *different* elements are advanced changes: after column k is final it
// is subtracted from all later columns, whose elements are independent of each
// other, so the innermost loop has no loop-carried dependency.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-156
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
// The loop nest is interchanged: column k is finished and then subtracted from every later
// column. Element (i, j) still receives its subtractions for k = 0..j-1 in order, then the division
// or square root.
fn x_factor_ltriang_pos_def_cols(a: &mut [Vec<f64>], n: usize, wc: &mut [f64]) {
    let s = X_WC_STRIDE;
    debug_assert!(n <= s && wc.len() >= s * s && s.is_multiple_of(2));
    // Column-major copy. Slots above a column's diagonal and below row n are
    // scratch: they are written by the padded loops below but never read as
    // matrix elements.
    for k in 0..n {
        let col = &mut wc[k * s..(k + 1) * s];
        col[..k].fill(0.0);
        for (i, slot) in col[k..n].iter_mut().enumerate() {
            *slot = a[k + i][k];
        }
        col[n..].fill(0.0);
    }
    let end = (n + 1) & !1;
    for k in 0..n {
        let (done, rest) = wc.split_at_mut((k + 1) * s);
        let ck = &mut done[k * s..(k + 1) * s];
        let d = ck[k].sqrt();
        ck[k] = d;
        for value in &mut ck[k + 1..n] {
            *value /= d;
        }
        let ck = &*ck;
        // A multiplier that is +0.0 subtracts `y * 0.0` from its whole target
        // column.  For a finite `y` that product is +0.0, or -0.0 when `y` has
        // its sign bit set; `x - (+0.0)` is `x` for every `x` (NaN, the
        // infinities and both zeros included) and `x - (-0.0)` is `x` except
        // that it turns -0.0 into +0.0.  So such a column is left alone, apart
        // from that one case in the few rows whose `y` is negative.  The Newton
        // matrix `A D A^T` has these multipliers throughout its leading diagonal
        // block (optimize_target_freq.c:119-147 only ever adds to `W[j][j]`,
        // `W[i + alphsize - 1][j]` and the second diagonal); its negative `y`
        // are in the relative-entropy row.
        const NEG_ROWS: usize = 4;
        let mut neg_rows = [0usize; NEG_ROWS];
        let mut neg_count = 0usize;
        // Integer tests on the bit patterns (all exponent bits set = infinity
        // or NaN; the top bit = negative), without branches.
        const EXPONENT: u64 = 0x7FF0_0000_0000_0000;
        let mut not_finite = 0u32;
        let mut signs = 0u64;
        for &v in &ck[k + 1..n] {
            let bits = v.to_bits();
            not_finite += u32::from(bits & EXPONENT == EXPONENT);
            signs |= bits;
        }
        let finite = not_finite == 0;
        if finite && (signs >> 63) != 0 {
            for (offset, &v) in ck[k + 1..n].iter().enumerate() {
                if v.is_sign_negative() {
                    if neg_count < NEG_ROWS {
                        neg_rows[neg_count] = k + 1 + offset;
                    }
                    neg_count += 1;
                }
            }
        }
        let plain = finite && neg_count <= NEG_ROWS;
        // The effect of the skipped update on column `column` (see above).
        let settle = |rest: &mut [f64], column: usize| {
            for &row in &neg_rows[..neg_count] {
                if row >= column {
                    let cell = &mut rest[(column - k - 1) * s + row];
                    if cell.to_bits() == (-0.0f64).to_bits() {
                        *cell = 0.0;
                    }
                }
            }
        };
        // Trailing update, four columns at a time. Rows start at an even
        // index <= j so every column pair of rows is handled by one 2-lane
        // step; the extra rows are scratch slots of those columns.
        let mut j = k + 1;
        while j + 4 <= n {
            let (s0, s1, s2, s3) = (ck[j], ck[j + 1], ck[j + 2], ck[j + 3]);
            if plain && (s0.to_bits() | s1.to_bits() | s2.to_bits() | s3.to_bits()) == 0 {
                if neg_count != 0 {
                    for column in j..j + 4 {
                        settle(rest, column);
                    }
                }
                j += 4;
                continue;
            }
            let i0 = j & !1;
            let len = end - i0;
            let block = &mut rest[(j - k - 1) * s..(j - k + 3) * s];
            let (c0, block) = block.split_at_mut(s);
            let (c1, block) = block.split_at_mut(s);
            let (c2, c3) = block.split_at_mut(s);
            let y = &ck[i0..i0 + len];
            let c0 = &mut c0[i0..i0 + len];
            let c1 = &mut c1[i0..i0 + len];
            let c2 = &mut c2[i0..i0 + len];
            let c3 = &mut c3[i0..i0 + len];
            for i in 0..len {
                let yi = y[i];
                c0[i] -= yi * s0;
                c1[i] -= yi * s1;
                c2[i] -= yi * s2;
                c3[i] -= yi * s3;
            }
            j += 4;
        }
        while j < n {
            let sj = ck[j];
            if plain && sj.to_bits() == 0 {
                if neg_count != 0 {
                    settle(rest, j);
                }
                j += 1;
                continue;
            }
            let i0 = j & !1;
            let len = end - i0;
            let cj = &mut rest[(j - k - 1) * s + i0..(j - k - 1) * s + i0 + len];
            for (x, &y) in cj.iter_mut().zip(&ck[i0..i0 + len]) {
                *x -= y * sj;
            }
            j += 1;
        }
    }
    for (i, row) in a.iter_mut().enumerate().take(n) {
        for (k, value) in row.iter_mut().enumerate().take(i + 1) {
            *value = wc[k * s + i];
        }
    }
}

// `Nlm_SolveLtriangPosDef`, forward substitution by columns.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:147-182
// x[i] still receives `-= L[i][j] * x[j]` for j = 0, ..., i-1 in that order
// and is then divided by L[i][i]; the second loop is unchanged.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:171-184
// Forward substitution is done by columns; the back substitution is the reference loop.
fn x_solve_ltriang_pos_def_cols(x: &mut [f64], n: usize, l: &[Vec<f64>], wc: &[f64]) {
    let s = X_WC_STRIDE;
    for j in 0..n {
        let cj = &wc[j * s..j * s + n];
        let xj = x[j] / cj[j];
        x[j] = xj;
        for (xi, &lij) in x[j + 1..n].iter_mut().zip(&cj[j + 1..n]) {
            *xi -= lij * xj;
        }
    }
    for j in (0..n).rev() {
        let row = &l[j];
        x[j] /= row[j];
        let xj = x[j];
        for (xi, &lji) in x[..j].iter_mut().zip(&row[..j]) {
            *xi -= lji * xj;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:175-193
fn euclidean_norm(v: &[f64]) -> f64 {
    let mut sum = 1.0;
    let mut scale = 0.0;
    for &value in v {
        if value != 0.0 {
            let abs_value = value.abs();
            if scale < abs_value {
                sum = 1.0 + sum * (scale / abs_value) * (scale / abs_value);
                scale = abs_value;
            } else {
                sum += (abs_value / scale) * (abs_value / scale);
            }
        }
    }
    scale * sum.sqrt()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:213
fn add_vectors(y: &mut [f64], alpha: f64, x: &[f64]) {
    for (left, right) in y.iter_mut().zip(x.iter()) {
        *left += alpha * *right;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:209-223
fn step_bound(x: &[f64], step_x: &[f64], max: f64) -> f64 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:224-238
    // The quotients first, then the scan (`x_step_bound`), for the 400 values of the Newton
    // iteration (`n = alphsize * alphsize` with `alphsize == COMPO_NUM_TRUE_AA`).
    debug_assert!(x.len() <= 400 && x.len() == step_x.len());
    x_step_bound(x, step_x, max)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:230-238
// The quotients are computed first and then scanned in the original order with the same test.
// `Nlm_StepBound`. The quotients do not depend on each other;
// compute them first, then scan them in the original order.
fn x_step_bound(x: &[f64], step_x: &[f64], max: f64) -> f64 {
    let mut ratio = [0.0f64; 400];
    for ((out, &x_i), &step_i) in ratio.iter_mut().zip(x).zip(step_x) {
        *out = -x_i / step_i;
    }
    let mut alpha = max;
    for &alpha_i in &ratio[..x.len()] {
        if alpha_i >= 0.0 && alpha_i < alpha {
            alpha = alpha_i;
        }
    }
    alpha
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:401-415
fn freq_ratio_to_score_square(
    matrix: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    rows: usize,
    cols: usize,
    lambda: f64,
) {
    for row in matrix.iter_mut().take(rows) {
        for value in row.iter_mut().take(cols) {
            if *value == 0.0 {
                *value = COMPO_SCORE_MIN as f64;
            } else {
                *value = value.glibc_ln() / lambda;
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:220-334
fn calc_lambda_full_precision_from_scores(
    score: &[Vec<f64>],
    alphsize: usize,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> (f64, usize) {
    let mut f = 4.0;
    let mut left = 0.0;
    let mut right = 1.0;
    let mut x = 0.367879441171f64;
    let mut is_newton = false;
    let mut max_score = COMPO_SCORE_MIN as f64;
    let mut avg_score = 0.0;

    for i in 0..alphsize {
        if row_prob[i] == 0.0 {
            continue;
        }
        for (j, &col_p) in col_prob.iter().enumerate().take(alphsize) {
            if col_p == 0.0 {
                continue;
            }
            if max_score < score[i][j] {
                max_score = score[i][j];
            }
            avg_score += row_prob[i] * col_p * score[i][j];
        }
    }
    if max_score <= 0.0 || avg_score >= 0.0 {
        return (-1.0, K_LAMBDA_ITERATION_LIMIT);
    }

    let mut k = 0usize;
    while k < K_LAMBDA_ITERATION_LIMIT {
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:266
        // `fold` takes f before f is evaluated again, so the Newton test
        // below compares this iteration's f with the previous iteration's f.
        let fold = f;
        let lambda = -x.glibc_ln();
        let was_newton = is_newton;
        let x_pow_max_score = (-max_score * lambda).glibc_exp();
        f = -x_pow_max_score;
        let mut slope = max_score * f / x;

        for i in 0..alphsize {
            if row_prob[i] == 0.0 {
                continue;
            }
            for (j, &col_p) in col_prob.iter().enumerate().take(alphsize) {
                if col_p == 0.0 {
                    continue;
                }
                let ff = if max_score != score[i][j] {
                    let diff_score = max_score - score[i][j];
                    let ff = row_prob[i] * col_p * (-lambda * diff_score).glibc_exp();
                    slope += diff_score * ff / x;
                    ff
                } else {
                    row_prob[i] * col_p
                };
                f += ff;
            }
        }

        if f > 0.0 {
            left = x;
        } else if f < 0.0 {
            right = x;
        } else {
            break;
        }
        if right - left <= 2.0 * left * (1.0 - right) * K_LAMBDA_ERROR_TOLERANCE
            && (f / x_pow_max_score).abs() <= K_LAMBDA_FUNCTION_TOLERANCE
        {
            x = (left + right) / 2.0;
            break;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:311-330
        if (was_newton && f.abs() > 0.9 * fold.abs()) || slope >= 0.0 {
            x = (left + right) / 2.0;
        } else {
            let p = -f / slope;
            let y = x + p;
            if y <= left || y >= right {
                x = (left + right) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() <= K_LAMBDA_ERROR_TOLERANCE * x * (1.0 - x)
                    && (f / x_pow_max_score).abs() <= K_LAMBDA_FUNCTION_TOLERANCE
                {
                    break;
                }
            }
        }
        k += 1;
    }

    if k < K_LAMBDA_ITERATION_LIMIT {
        (-x.glibc_ln(), k)
    } else {
        (-1.0, k)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:136-147
pub fn blast_get_relative_entropy(
    a: &[f64; COMPO_NUM_TRUE_AA],
    b: &[f64; COMPO_NUM_TRUE_AA],
) -> f64 {
    let mut value = 0.0;
    for i in 0..COMPO_NUM_TRUE_AA {
        let temp = (a[i] + b[i]) / 2.0;
        if temp > 0.0 {
            if a[i] > 0.0 {
                value += a[i] * (a[i] / temp).glibc_ln() / 2.0;
            }
            if b[i] > 0.0 {
                value += b[i] * (b[i] / temp).glibc_ln() / 2.0;
            }
        }
    }
    if value < 0.0 {
        0.0
    } else {
        value.sqrt()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:89-107
fn apply_pseudocounts(
    probs20: &mut [f64; COMPO_NUM_TRUE_AA],
    number_of_observations: i32,
    background_probs20: &[f64; COMPO_NUM_TRUE_AA],
    pseudocounts: i32,
) {
    let mut sum = 0.0;
    for &value in probs20.iter() {
        sum += value;
    }
    if sum == 0.0 {
        sum = 1.0;
    }
    let weight = pseudocounts as f64 / (number_of_observations as f64 + pseudocounts as f64);
    for i in 0..COMPO_NUM_TRUE_AA {
        probs20[i] = (1.0 - weight) * probs20[i] / sum + weight * background_probs20[i];
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:428-441
fn gather_letter_probs(
    output_letter_probs: &mut [f64; COMPO_NUM_TRUE_AA],
    input_letter_probs: &[f64; COMPO_LARGEST_ALPHABET],
) {
    for c in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[c] >= 0 {
            output_letter_probs[ALPHA_CONVERT[c] as usize] = input_letter_probs[c];
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:444-458
fn unpack_letter_probs(
    std_probs: &mut [f64; COMPO_LARGEST_ALPHABET],
    probs: &[f64; COMPO_NUM_TRUE_AA],
) {
    for c in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[c] >= 0 {
            std_probs[c] = probs[ALPHA_CONVERT[c] as usize];
        } else {
            std_probs[c] = 0.0;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:467-478
fn set_pair_ambig_probs_to_sum(probs: &mut [f64; COMPO_LARGEST_ALPHABET]) {
    probs[E_BCHAR] = probs[E_DCHAR] + probs[E_NCHAR];
    probs[E_ZCHAR] = probs[E_ECHAR] + probs[E_QCHAR];
    probs[E_JCHAR] = probs[E_ICHAR] + probs[E_LCHAR];
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:538-604
fn true_aa_to_std_target_freqs(
    std_freq: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
) {
    let mut sum = 0.0;
    for row in freq {
        for &value in row {
            sum += value;
        }
    }
    for a in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[a] < 0 {
            for b in 0..COMPO_LARGEST_ALPHABET {
                std_freq[a][b] = 0.0;
            }
        } else {
            let small_a = ALPHA_CONVERT[a] as usize;
            for b in 0..COMPO_LARGEST_ALPHABET {
                if ALPHA_CONVERT[b] < 0 {
                    std_freq[a][b] = 0.0;
                } else {
                    let small_b = ALPHA_CONVERT[b] as usize;
                    std_freq[a][b] = freq[small_a][small_b] / sum;
                }
            }
            std_freq[a][E_BCHAR] = std_freq[a][E_DCHAR] + std_freq[a][E_NCHAR];
            std_freq[a][E_ZCHAR] = std_freq[a][E_ECHAR] + std_freq[a][E_QCHAR];
            std_freq[a][E_JCHAR] = std_freq[a][E_ICHAR] + std_freq[a][E_LCHAR];
        }
    }
    for b in 0..COMPO_LARGEST_ALPHABET {
        std_freq[E_BCHAR][b] = std_freq[E_DCHAR][b] + std_freq[E_NCHAR][b];
        std_freq[E_ZCHAR][b] = std_freq[E_ECHAR][b] + std_freq[E_QCHAR][b];
        std_freq[E_JCHAR][b] = std_freq[E_ICHAR][b] + std_freq[E_LCHAR][b];
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:607-639
fn calc_avg_score(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    row: usize,
    col_probs: &[f64; COMPO_LARGEST_ALPHABET],
) -> f64 {
    let mut score = 0.0;
    for col in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[col] >= 0 {
            score += scores[row][col] * col_probs[col];
        }
    }
    score
}

fn calc_x_score(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    col: usize,
    row_probs: &[f64; COMPO_LARGEST_ALPHABET],
) -> f64 {
    let mut score = 0.0;
    for row in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[row] >= 0 {
            score += scores[row][col] * row_probs[row];
        }
    }
    score.min(K_MAXIMUM_X_SCORE)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:708-711
fn nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:648-690
fn set_xuo_scores(
    scores: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    row_probs: &[f64; COMPO_LARGEST_ALPHABET],
    col_probs: &[f64; COMPO_LARGEST_ALPHABET],
) {
    let mut score_xx = 0.0;
    for i in 0..COMPO_LARGEST_ALPHABET {
        if ALPHA_CONVERT[i] >= 0 {
            let avg_ix = calc_avg_score(scores, i, col_probs);
            scores[i][E_XCHAR] = avg_ix.min(K_MAXIMUM_X_SCORE);
            score_xx += avg_ix * row_probs[i];
            scores[E_XCHAR][i] = calc_x_score(scores, i, row_probs);
        }
    }
    scores[E_XCHAR][E_XCHAR] = score_xx.min(K_MAXIMUM_X_SCORE);
    scores[E_BCHAR][E_XCHAR] = calc_avg_score(scores, E_BCHAR, col_probs).min(K_MAXIMUM_X_SCORE);
    scores[E_XCHAR][E_BCHAR] = calc_x_score(scores, E_BCHAR, row_probs);
    scores[E_ZCHAR][E_XCHAR] = calc_avg_score(scores, E_ZCHAR, col_probs).min(K_MAXIMUM_X_SCORE);
    scores[E_XCHAR][E_ZCHAR] = calc_x_score(scores, E_ZCHAR, row_probs);
    scores[E_JCHAR][E_XCHAR] = calc_avg_score(scores, E_JCHAR, col_probs).min(K_MAXIMUM_X_SCORE);
    scores[E_XCHAR][E_JCHAR] = calc_x_score(scores, E_JCHAR, row_probs);

    scores[E_SELENOCYSTEINE] = scores[E_CCHAR];
    for row in scores.iter_mut() {
        row[E_SELENOCYSTEINE] = row[E_CCHAR];
    }
    scores[E_OCHAR] = scores[E_XCHAR];
    for row in scores.iter_mut() {
        row[E_OCHAR] = row[E_XCHAR];
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:701-734
fn round_score_matrix(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
) -> [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET] {
    let mut matrix = [[0i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
    for row in 0..COMPO_LARGEST_ALPHABET {
        for col in 0..COMPO_LARGEST_ALPHABET {
            matrix[row][col] = if scores[row][col] < i32::MIN as f64 {
                i32::MIN
            } else {
                nint(scores[row][col])
            };
        }
    }
    matrix
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:945-958
fn int4_matrix_from_freq(
    matrix: &mut [[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    alphsize: usize,
    freq: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    lambda: f64,
) {
    let mut scores = *freq;
    freq_ratio_to_score_square(&mut scores, alphsize, alphsize, lambda);
    let rounded = round_score_matrix(&scores);
    for row in 0..alphsize {
        matrix[row][..alphsize].copy_from_slice(&rounded[row][..alphsize]);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:760-793
fn scores_std_alphabet(
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    matrix_info: &BlastMatrixInfo,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> AdjustedProteinMatrix {
    let scores = scores_std_alphabet_unrounded(target_freq, matrix_info, row_prob, col_prob);
    scores_std_alphabet_round(&scores, matrix_info)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:779-788
// This function ports lines 779-788 of `s_ScoresStdAlphabet`; `scores_std_alphabet_round` ports lines 790-796.
/// The part of `s_ScoresStdAlphabet` before the scores are rounded.
fn scores_std_alphabet_unrounded(
    target_freq: &[[f64; COMPO_NUM_TRUE_AA]; COMPO_NUM_TRUE_AA],
    matrix_info: &BlastMatrixInfo,
    row_prob: &[f64; COMPO_NUM_TRUE_AA],
    col_prob: &[f64; COMPO_NUM_TRUE_AA],
) -> [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET] {
    let mut row_prob_std = [0.0; COMPO_LARGEST_ALPHABET];
    let mut col_prob_std = [0.0; COMPO_LARGEST_ALPHABET];
    unpack_letter_probs(&mut row_prob_std, row_prob);
    set_pair_ambig_probs_to_sum(&mut row_prob_std);
    unpack_letter_probs(&mut col_prob_std, col_prob);
    set_pair_ambig_probs_to_sum(&mut col_prob_std);

    let mut scores = [[0.0; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET];
    true_aa_to_std_target_freqs(&mut scores, target_freq);
    calc_freq_ratios(&mut scores, &row_prob_std, &col_prob_std);
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:403-413
    // `Blast_FreqRatioToScore`: every non-zero entry becomes `log(entry) / Lambda`. The logs of all non-zero
    // entries are computed first in one slice call, then divided by Lambda in the same entry order.
    // ln of every non-zero ratio (through x_logclone::log_slice: libm unless
    // LOSAT_X_LOGCLONE is set), then / Lambda, as Blast_FreqRatioToScore does.
    let mut ratios = [0.0f64; COMPO_LARGEST_ALPHABET * COMPO_LARGEST_ALPHABET];
    let mut logs = [0.0f64; COMPO_LARGEST_ALPHABET * COMPO_LARGEST_ALPHABET];
    let mut count = 0usize;
    for row in &scores {
        for &value in row.iter() {
            if value != 0.0 {
                ratios[count] = value;
                count += 1;
            }
        }
    }
    crate::utils::x_logclone::log_slice(&ratios[..count], &mut logs[..count]);
    let mut next = 0usize;
    for row in &mut scores {
        for value in row.iter_mut() {
            if *value == 0.0 {
                *value = COMPO_SCORE_MIN as f64;
            } else {
                *value = logs[next] / matrix_info.ungapped_lambda;
                next += 1;
            }
        }
    }
    set_xuo_scores(&mut scores, &row_prob_std, &col_prob_std);
    scores
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:790-796
/// The rounding and the stop-codon scores of `s_ScoresStdAlphabet`.
fn scores_std_alphabet_round(
    scores: &[[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    matrix_info: &BlastMatrixInfo,
) -> AdjustedProteinMatrix {
    let mut rounded = round_score_matrix(scores);
    for index in 0..COMPO_LARGEST_ALPHABET {
        rounded[index][E_STOP_CHAR] = matrix_info.start_matrix[index][E_STOP_CHAR];
        rounded[E_STOP_CHAR][index] = matrix_info.start_matrix[E_STOP_CHAR][index];
    }
    AdjustedProteinMatrix { scores: rounded }
}

// No NCBI counterpart: reads LOSAT_X_ADJMEMO_STATS once; it does not change any value NCBI computes.
fn x_adjmemo_enabled() -> bool {
    if !crate::utils::xdiag::ON {
        return false;
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xdiag::var_os("LOSAT_X_ADJMEMO_STATS").is_some())
}
// No NCBI counterpart: call counter and set of hashed inputs for the LOSAT_X_ADJMEMO_STATS report; it does not change any value NCBI computes.
static X_ADJMEMO: std::sync::Mutex<Option<(u64, std::collections::HashSet<u64>)>> =
    std::sync::Mutex::new(None);
// No NCBI counterpart: counts calls of `Blast_AdjustScores` that reach the matrix adjustment, and distinct inputs; it does not change any value NCBI computes.
fn x_adjmemo_record(
    query_composition: &BlastAminoAcidComposition,
    subject_composition: &BlastAminoAcidComposition,
    rule: EMatrixAdjustRule,
) {
    if !x_adjmemo_enabled() {
        return;
    }
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    for v in query_composition.prob.iter() {
        v.to_bits().hash(&mut h);
    }
    query_composition.num_true_amino_acids.hash(&mut h);
    for v in subject_composition.prob.iter() {
        v.to_bits().hash(&mut h);
    }
    subject_composition.num_true_amino_acids.hash(&mut h);
    (rule as u8).hash(&mut h);
    let key = h.finish();
    let mut guard = X_ADJMEMO.lock().expect("adjmemo lock");
    let entry = guard.get_or_insert_with(|| (0, std::collections::HashSet::new()));
    entry.0 += 1;
    entry.1.insert(key);
}
// No NCBI counterpart: prints the LOSAT_X_STATS counters of the Newton lanes; it does not change any
// value NCBI computes.
/// Counter summary of the Newton lanes (printed by `main` at exit).
pub fn x_newtonexact_print_stats() {
    super::x_newton_lanes::print_stats();
}

// No NCBI counterpart: prints the LOSAT_X_ADJMEMO_STATS counters; it does not change any value NCBI computes.
pub fn x_adjmemo_print_stats() {
    if !x_adjmemo_enabled() {
        return;
    }
    let guard = X_ADJMEMO.lock().expect("adjmemo lock");
    if let Some((calls, set)) = guard.as_ref() {
        eprintln!(
            "[X_STATS] ADJMEMO matrix_adjust_calls={} distinct_inputs={}",
            calls,
            set.len()
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:95-126
fn scaled_symmetric_product_a(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
    // The register-accumulator form (`x_scaled_symmetric_product_a`) for `alphsize` up to 20
    // (`optimize_target_frequencies` asserts `alphsize == COMPO_NUM_TRUE_AA`).
    debug_assert!(alphsize <= COMPO_NUM_TRUE_AA);
    x_scaled_symmetric_product_a(w, diagonal, alphsize)
}

// `ScaledSymmetricProductA` with register accumulators.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:108-135
// W[j][j] still accumulates dd for i = 0, 1, ... in order starting from 0.0,
// W[i+alphsize-1][i+alphsize-1] accumulates for j = 0, 1, ... in order starting
// from 0.0, and every W[i+alphsize-1][j] is the single sum 0.0 + dd.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:126-137
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:112-113,121-125
// The lower triangle of the first `m` rows is zeroed first, as here.
fn x_scaled_symmetric_product_a(w: &mut [Vec<f64>], diagonal: &[f64], alphsize: usize) {
    let n = alphsize;
    let m = 2 * n - 1;
    for (row_idx, row) in w.iter_mut().enumerate().take(m) {
        for value in row.iter_mut().take(row_idx + 1) {
            *value = 0.0;
        }
    }
    let mut col_sum = [0.0f64; COMPO_NUM_TRUE_AA];
    for i in 0..n {
        let row = &diagonal[i * n..(i + 1) * n];
        for (acc, &dd) in col_sum[..n].iter_mut().zip(row) {
            *acc += dd;
        }
    }
    for j in 0..n {
        w[j][j] = col_sum[j];
    }
    for i in 1..n {
        let row = &diagonal[i * n..(i + 1) * n];
        let w_row = &mut w[i + n - 1];
        let mut acc = 0.0f64;
        for (slot, &dd) in w_row[..n].iter_mut().zip(row) {
            *slot = 0.0 + dd;
            acc += dd;
        }
        w_row[i + n - 1] = acc;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:138-165
// The same additions per element in the same order, arranged so that independent elements
// advance together.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:170-177
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:170-179
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
fn multiply_by_a(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for value in y.iter_mut() {
            *value = 0.0;
        }
    } else if beta != 1.0 {
        for value in y.iter_mut() {
            *value *= beta;
        }
    }
    let n = alphsize;
    {
        let head = &mut y[..n];
        for i in 0..n {
            let row = &x[i * n..(i + 1) * n];
            for (yj, &xv) in head.iter_mut().zip(row) {
                *yj += alpha * xv;
            }
        }
    }
    // y[i + n - 1] accumulates x[i*n + j] for j = 0..n in order; four rows at
    // a time, each with its own accumulator.
    let mut i = 1;
    while i + 4 <= n {
        let r0 = &x[i * n..(i + 1) * n];
        let r1 = &x[(i + 1) * n..(i + 2) * n];
        let r2 = &x[(i + 2) * n..(i + 3) * n];
        let r3 = &x[(i + 3) * n..(i + 4) * n];
        let (mut a0, mut a1, mut a2, mut a3) = (y[i + n - 1], y[i + n], y[i + n + 1], y[i + n + 2]);
        for j in 0..n {
            a0 += alpha * r0[j];
            a1 += alpha * r1[j];
            a2 += alpha * r2[j];
            a3 += alpha * r3[j];
        }
        y[i + n - 1] = a0;
        y[i + n] = a1;
        y[i + n + 1] = a2;
        y[i + n + 2] = a3;
        i += 4;
    }
    while i < n {
        let row = &x[i * n..(i + 1) * n];
        let mut acc = y[i + n - 1];
        for &xv in row {
            acc += alpha * xv;
        }
        y[i + n - 1] = acc;
        i += 1;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:176-205
// The same two additions per element.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:214-220
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:214-223
// (The line numbers in the older comment above are for a different revision; these are the pinned ones.)
fn multiply_by_at(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for value in y.iter_mut() {
            *value = 0.0;
        }
    } else if beta != 1.0 {
        for value in y.iter_mut() {
            *value *= beta;
        }
    }
    let n = alphsize;
    let head = &x[..n];
    for i in 0..n {
        let row = &mut y[i * n..(i + 1) * n];
        if i > 0 {
            let xi = x[i + n - 1];
            for (yk, &xj) in row.iter_mut().zip(head) {
                let t = *yk + alpha * xj;
                *yk = t + alpha * xi;
            }
        } else {
            for (yk, &xj) in row.iter_mut().zip(head) {
                *yk += alpha * xj;
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:216-231
fn residuals_linear_constraints(
    ra: &mut [f64],
    alphsize: usize,
    x: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
) {
    ra[..alphsize].copy_from_slice(&col_sums[..alphsize]);
    for i in 1..alphsize {
        ra[i + alphsize - 1] = row_sums[i];
    }
    multiply_by_a(1.0, ra, alphsize, -1.0, x);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:241-262
fn dual_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    grads: &[Vec<f64>],
    z: &[f64],
    constrain_rel_entropy: bool,
) {
    if constrain_rel_entropy {
        let eta = z[2 * alphsize - 1];
        for i in 0..resids_x.len() {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    } else {
        for i in 0..resids_x.len() {
            resids_x[i] = -grads[0][i];
        }
    }
    multiply_by_at(1.0, resids_x, alphsize, 1.0, z);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:282-323
fn calculate_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    resids_z: &mut [f64],
    values: &[f64; 2],
    grads: &[Vec<f64>],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    x: &[f64],
    z: &[f64],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> f64 {
    dual_residuals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    let norm_resids_x = euclidean_norm(resids_x);
    residuals_linear_constraints(resids_z, alphsize, x, row_sums, col_sums);
    let norm_resids_z = if constrain_rel_entropy {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];
        euclidean_norm(&resids_z[..(2 * alphsize)])
    } else {
        euclidean_norm(&resids_z[..(2 * alphsize - 1)])
    };
    (norm_resids_x * norm_resids_x + norm_resids_z * norm_resids_z).sqrt()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:366-382
struct ReNewtonSystem {
    alphsize: usize,
    constrain_rel_entropy: bool,
    // Column-major copy of the factor (see x_factor_ltriang_pos_def_cols).
    // No NCBI counterpart: a column-major copy of the lower triangle of `W` after the factorisation; it does not change any value NCBI computes.
    x_wc: Vec<f64>,
    w: Vec<Vec<f64>>,
    dinv: Vec<f64>,
    grad_re: Vec<f64>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:419
fn re_newton_system_new(alphsize: usize) -> ReNewtonSystem {
    ReNewtonSystem {
        alphsize,
        constrain_rel_entropy: true,
        // No NCBI counterpart: the column-major copy of the factor; it does not change any value NCBI computes.
        x_wc: vec![0.0; X_WC_STRIDE * X_WC_STRIDE],
        w: ltriang_matrix_new(2 * alphsize),
        dinv: vec![0.0; alphsize * alphsize],
        grad_re: vec![0.0; alphsize * alphsize],
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:408-470
fn factor_re_newton_system(
    newton_system: &mut ReNewtonSystem,
    x: &[f64],
    z: &[f64],
    grads: &[Vec<f64>],
    constrain_rel_entropy: bool,
    workspace: &mut [f64],
) {
    let alphsize = newton_system.alphsize;
    let n = alphsize * alphsize;
    let m = if constrain_rel_entropy {
        2 * alphsize
    } else {
        2 * alphsize - 1
    };

    newton_system.constrain_rel_entropy = constrain_rel_entropy;
    if constrain_rel_entropy {
        let eta = z[m - 1];
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:507-510
        // The divisor `1 - eta` is the same value for every element, so computing it once gives the
        // same quotients.
        let divisor = 1.0 - eta;
        for (out, &xi) in newton_system.dinv[..n].iter_mut().zip(&x[..n]) {
            *out = xi / divisor;
        }
    } else {
        newton_system.dinv.copy_from_slice(x);
    }

    scaled_symmetric_product_a(&mut newton_system.w, &newton_system.dinv, alphsize);
    if constrain_rel_entropy {
        newton_system.grad_re.copy_from_slice(&grads[1]);
        newton_system.w[m - 1][m - 1] = 0.0;
        for i in 0..n {
            workspace[i] = newton_system.dinv[i] * newton_system.grad_re[i];
            newton_system.w[m - 1][m - 1] += newton_system.grad_re[i] * workspace[i];
        }
        multiply_by_a(
            0.0,
            &mut newton_system.w[m - 1][..m - 1],
            alphsize,
            1.0,
            workspace,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:534-534
    // The column-major form (`x_factor_ltriang_pos_def_cols`) for `m` up to `X_WC_STRIDE`
    // (`m <= 2 * COMPO_NUM_TRUE_AA` as `alphsize == COMPO_NUM_TRUE_AA`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:143-149
    debug_assert!(m <= X_WC_STRIDE);
    x_factor_ltriang_pos_def_cols(&mut newton_system.w, m, &mut newton_system.x_wc);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:480-524
fn solve_re_newton_system(
    x: &mut [f64],
    z: &mut [f64],
    newton_system: &ReNewtonSystem,
    workspace: &mut [f64],
) {
    let alphsize = newton_system.alphsize;
    let n = alphsize * alphsize;
    let ma = 2 * alphsize - 1;
    let m = if newton_system.constrain_rel_entropy {
        ma + 1
    } else {
        ma
    };

    for i in 0..n {
        workspace[i] = x[i] * newton_system.dinv[i];
    }
    multiply_by_a(1.0, z, alphsize, -1.0, workspace);
    if newton_system.constrain_rel_entropy {
        for i in 0..n {
            z[m - 1] -= newton_system.grad_re[i] * workspace[i];
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:587-587
    // The column-major form (`x_solve_ltriang_pos_def_cols`), which reads the column-major copy
    // that the factorisation made.
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/nlm_linear_algebra.c:171-177
    debug_assert!(m <= X_WC_STRIDE);
    x_solve_ltriang_pos_def_cols(z, m, &newton_system.w, &newton_system.x_wc);
    if newton_system.constrain_rel_entropy {
        for i in 0..n {
            x[i] += newton_system.grad_re[i] * z[m - 1];
        }
    }
    multiply_by_at(1.0, x, alphsize, 1.0, z);
    for i in 0..n {
        x[i] *= newton_system.dinv[i];
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:538-565
fn evaluate_re_functions(
    values: &mut [f64; 2],
    grads: &mut [Vec<f64>],
    alphsize: usize,
    x: &[f64],
    q: &[f64],
    scores: &[f64],
    constrain_rel_entropy: bool,
) {
    values[0] = 0.0;
    values[1] = 0.0;
    for k in 0..alphsize * alphsize {
        let mut temp = (x[k] / q[k]).glibc_ln();
        values[0] += x[k] * temp;
        grads[0][k] = temp + 1.0;
        if constrain_rel_entropy {
            temp += scores[k];
            values[1] += x[k] * temp;
            grads[1][k] = temp + 1.0;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:578-598
fn compute_scores_from_probs(
    scores: &mut [f64],
    alphsize: usize,
    target_freqs: &[f64],
    row_freqs: &[f64; COMPO_NUM_TRUE_AA],
    col_freqs: &[f64; COMPO_NUM_TRUE_AA],
) {
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            scores[k] = (target_freqs[k] / (row_freqs[i] * col_freqs[j])).glibc_ln();
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:602-713
fn optimize_target_frequencies(
    x: &mut [f64],
    alphsize: usize,
    q: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1390-1392
    // No NCBI counterpart: the assert. NCBI's only caller (above) and LOSAT's
    // (`composition_matrix_adjust`) pass `COMPO_NUM_TRUE_AA`; the reordered forms (`x_step_bound`,
    // `x_scaled_symmetric_product_a`, `x_factor_ltriang_pos_def_cols`, `x_solve_ltriang_pos_def_cols`)
    // assume it (fixed-size storage for `n = 400` values and `m <= X_WC_STRIDE` rows). It does not
    // change any value NCBI computes.
    assert_eq!(alphsize, COMPO_NUM_TRUE_AA);
    // The same iteration, bit for bit, in fixed-size storage (see x_newton_exact.rs).
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
    // The port of this function is `optimize_target_frequencies_reference`. The exact module handles
    // only the problem with `constrain_rel_entropy`; the port runs the other problems.
    if constrain_rel_entropy {
        // A result computed ahead (four problems per vector, eight on aarch64, x_newton_lanes.rs)
        // for exactly this input.
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
        // No NCBI counterpart: the stored result was computed by the operations of this function on an
        // input with the same bits (q, row and column sums, relative entropy). It does not change any
        // value NCBI computes.
        if let Some(status) =
            super::x_newton_lanes::take(x, q, row_sums, col_sums, relative_entropy)
        {
            return status;
        }
        // No NCBI counterpart: call counter and timer of the xstats report; it does not change any value NCBI computes.
        let x_t0 = crate::utils::xstats::now();
        crate::utils::xstats::add(&crate::utils::xstats::COMPO_ADJ_CALLS, 1);
        let status = super::x_newton_exact::optimize_target_frequencies_exact(
            x,
            q,
            row_sums,
            col_sums,
            relative_entropy,
        );
        // No NCBI counterpart: timer of the xstats report; it does not change any value NCBI computes.
        crate::utils::xstats::add_ns(&crate::utils::xstats::COMPO_NS, x_t0);
        return status;
    }
    optimize_target_frequencies_reference(
        x,
        alphsize,
        q,
        row_sums,
        col_sums,
        constrain_rel_entropy,
        relative_entropy,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-696
// The port below with `alphsize == 20` and `constrain_rel_entropy`, for the bitwise tests of
// `x_newton_lanes`.
#[cfg(test)]
pub(super) fn x_optimize_reference_for_test(
    x: &mut [f64],
    q: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    relative_entropy: f64,
) -> i32 {
    optimize_target_frequencies_reference(
        x,
        COMPO_NUM_TRUE_AA,
        q,
        row_sums,
        col_sums,
        true,
        relative_entropy,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:686-687,693-697,745-748
// The port of this function, for the problems the exact module does not take (no relative-entropy
// constraint) and the reference of the bitwise tests of `x_newton_lanes`.
fn optimize_target_frequencies_reference(
    x: &mut [f64],
    alphsize: usize,
    q: &[f64],
    row_sums: &[f64; COMPO_NUM_TRUE_AA],
    col_sums: &[f64; COMPO_NUM_TRUE_AA],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) -> i32 {
    let n = alphsize * alphsize;
    let ma = 2 * alphsize - 1;
    let m = if constrain_rel_entropy { ma + 1 } else { ma };
    // No NCBI counterpart: call counter and timer of the xstats report; it does not change any value NCBI computes.
    let x_t0 = crate::utils::xstats::now();
    crate::utils::xstats::add(&crate::utils::xstats::COMPO_ADJ_CALLS, 1);

    let mut newton_system = re_newton_system_new(alphsize);
    let mut resids_x = vec![0.0; n];
    let mut resids_z = vec![0.0; ma + 1];
    let mut z = vec![0.0; ma + 1];
    let mut old_scores = vec![0.0; n];
    let mut workspace = vec![0.0; n];
    let mut grads = dense_matrix_new(2, n);
    let mut values = [0.0; 2];

    compute_scores_from_probs(&mut old_scores, alphsize, q, row_sums, col_sums);
    x.copy_from_slice(q);
    let mut rnorm = 0.0;
    let mut its = 0usize;
    while its <= K_COMPO_ADJUST_ITERATION_LIMIT {
        evaluate_re_functions(
            &mut values,
            &mut grads,
            alphsize,
            x,
            q,
            &old_scores,
            constrain_rel_entropy,
        );
        rnorm = calculate_residuals(
            &mut resids_x,
            alphsize,
            &mut resids_z,
            &values,
            &grads,
            row_sums,
            col_sums,
            x,
            &z,
            constrain_rel_entropy,
            relative_entropy,
        );
        #[expect(
            clippy::neg_cmp_op_on_partial_ord,
            reason = "NCBI's convergence test `!(rnorm > tol)` stops on a NaN rnorm (optimize_target_freq.c:760-762); `<=` would not"
        )]
        if !(rnorm > K_COMPO_ADJUST_ERR_TOLERANCE) {
            break;
        }
        its += 1;
        if its <= K_COMPO_ADJUST_ITERATION_LIMIT {
            // No NCBI counterpart: Newton iteration counter of the xstats report; it does not change any value NCBI computes.
            crate::utils::xstats::add(&crate::utils::xstats::NEWTON_ITERS, 1);
            factor_re_newton_system(
                &mut newton_system,
                x,
                &z,
                &grads,
                constrain_rel_entropy,
                &mut workspace,
            );
            solve_re_newton_system(
                &mut resids_x,
                &mut resids_z[..m],
                &newton_system,
                &mut workspace,
            );
            let mut alpha = step_bound(x, &resids_x, 1.0 / 0.95);
            alpha *= 0.95;
            add_vectors(x, alpha, &resids_x);
            add_vectors(&mut z[..m], alpha, &resids_z[..m]);
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:786-794
    // No NCBI counterpart: timer of the xstats report; it does not change any value NCBI computes.
    crate::utils::xstats::add_ns(&crate::utils::xstats::COMPO_NS, x_t0);
    let converged = its <= K_COMPO_ADJUST_ITERATION_LIMIT
        && rnorm <= K_COMPO_ADJUST_ERR_TOLERANCE
        && (!constrain_rel_entropy || z[m - 1] < 1.0);
    if converged {
        0
    } else {
        1
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1336-1408
fn composition_matrix_adjust(
    matrix_info: &BlastMatrixInfo,
    matrix_adjust_rule: EMatrixAdjustRule,
    query_length: i32,
    subject_length: i32,
    query_prob_std: &[f64; COMPO_LARGEST_ALPHABET],
    subject_prob_std: &[f64; COMPO_LARGEST_ALPHABET],
    pseudocounts: i32,
    specified_re: f64,
    workspace: &mut BlastCompositionWorkspace,
) -> Result<Option<AdjustedProteinMatrix>> {
    let mut row_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut col_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut row_probs, query_prob_std);
    gather_letter_probs(&mut col_probs, subject_prob_std);

    let desired_re = match matrix_adjust_rule {
        EMatrixAdjustRule::UserSpecifiedRelEntropy => specified_re,
        _ => bail!(
            "unsupported matrix-adjust rule {:?} for composition matrix adjustment",
            matrix_adjust_rule
        ),
    };

    apply_pseudocounts(
        &mut row_probs,
        query_length,
        &workspace.first_standard_freq,
        pseudocounts,
    );
    apply_pseudocounts(
        &mut col_probs,
        subject_length,
        &workspace.second_standard_freq,
        pseudocounts,
    );

    let mut x = vec![0.0; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
    let mut q = vec![0.0; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            q[i * COMPO_NUM_TRUE_AA + j] = workspace.mat_b[i][j];
        }
    }
    let status = optimize_target_frequencies(
        &mut x,
        COMPO_NUM_TRUE_AA,
        &q,
        &row_probs,
        &col_probs,
        desired_re > 0.0,
        desired_re,
    );
    if status != 0 {
        return Ok(None);
    }
    for i in 0..COMPO_NUM_TRUE_AA {
        for j in 0..COMPO_NUM_TRUE_AA {
            workspace.mat_final[i][j] = x[i * COMPO_NUM_TRUE_AA + j];
        }
    }
    Ok(Some(scores_std_alphabet(
        &workspace.mat_final,
        matrix_info,
        &row_probs,
        &col_probs,
    )))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/unified_pvalues.c:208-260
pub fn blast_composition_pvalue(lambda: f64) -> f64 {
    let initial_u_score = (lambda - COMPO_MIN_LAMBDA) / LAMBDA_BIN_SIZE;
    if initial_u_score < LOW_LAMBDA_BIN_CUT as f64 {
        return P_LAMBDA_TABLE[LOW_LAMBDA_BIN_CUT];
    }
    if initial_u_score > (LAMBDA_BIN_TOTAL - 1) as f64 {
        return 1.0;
    }
    let bin = initial_u_score as usize;
    if bin == LAMBDA_BIN_TOTAL - 1 {
        return P_LAMBDA_TABLE[LAMBDA_BIN_TOTAL - 1];
    }
    let delta = initial_u_score - bin as f64;
    (1.0 - delta) * P_LAMBDA_TABLE[bin] + delta * P_LAMBDA_TABLE[bin + 1]
}

pub fn blast_overall_p_value(p_comp: f64, p_alignment: f64) -> f64 {
    let product = p_comp * p_alignment;
    if product > 0.0 {
        product * (1.0 - product.glibc_ln())
    } else {
        product
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1414-1530
pub(crate) fn blast_adjust_scores(
    matrix_info: &BlastMatrixInfo,
    query_composition: &BlastAminoAcidComposition,
    query_length: i32,
    subject_composition: &BlastAminoAcidComposition,
    subject_length: i32,
    composition_adjust_mode: BlastCompoAdjustMode,
    composition_test_index: i32,
    workspace: &mut BlastCompositionWorkspace,
    calc_lambda: fn(&[f64], i32, i32, f64) -> Result<f64>,
) -> Result<Option<BlastAdjustScoresResult>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1433-1444
    if query_composition.num_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0
    {
        return Ok(None);
    }

    let mut permuted_query_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut permuted_match_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut permuted_query_probs, &query_composition.prob);
    gather_letter_probs(&mut permuted_match_probs, &subject_composition.prob);

    let pvalue_for_this_pair = if composition_test_index > 0 {
        Some(compute_pair_composition_pvalue(
            &permuted_query_probs,
            &permuted_match_probs,
        )?)
    } else {
        None
    };

    let mut matrix_adjust_rule = if matrix_info.position_based
        || composition_adjust_mode == BlastCompoAdjustMode::CompositionBasedStats
    {
        EMatrixAdjustRule::CompoScaleOldMatrix
    } else {
        blast_choose_matrix_adjust_rule(
            query_length,
            subject_length,
            &permuted_query_probs,
            &permuted_match_probs,
            matrix_info.matrix,
            composition_adjust_mode,
        )?
    };
    if matrix_adjust_rule != EMatrixAdjustRule::CompoScaleOldMatrix {
        // EXPERIMENT (LOSAT_X_ADJMEMO_STATS): how often the same inputs recur
        // (a memo of the adjusted matrix would be exact: the function is pure).
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1501-1507
        // No NCBI counterpart: LOSAT_X_ADJMEMO_STATS counts how many times the call that follows (the port of
        // `Blast_CompositionMatrixAdj`) is made and with how many distinct inputs; it does not change any value NCBI computes.
        x_adjmemo_record(query_composition, subject_composition, matrix_adjust_rule);
        if let Some(adjusted_matrix) = composition_matrix_adjust(
            matrix_info,
            matrix_adjust_rule,
            query_composition.num_true_amino_acids,
            subject_composition.num_true_amino_acids,
            &query_composition.prob,
            &subject_composition.prob,
            20,
            K_FIXED_RE_BLOSUM62,
            workspace,
        )? {
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1509-1517
            return Ok(Some(BlastAdjustScoresResult {
                matrix_adjust_rule,
                adjusted_matrix,
                pvalue_for_this_pair,
                lambda_ratio: 1.0,
            }));
        }
        matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
    }

    let lambda_ratio = compute_lambda_ratio(
        matrix_info,
        &query_composition.prob,
        &subject_composition.prob,
        composition_test_index > 0,
        calc_lambda,
    );
    let adjusted_matrix = build_scaled_square_matrix(
        matrix_info,
        &query_composition.prob,
        &subject_composition.prob,
        lambda_ratio,
    );

    Ok(Some(BlastAdjustScoresResult {
        matrix_adjust_rule,
        adjusted_matrix,
        pvalue_for_this_pair,
        lambda_ratio,
    }))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1485-1505
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1354-1355,1385-1399
// Newton lanes: the prefix of `blast_adjust_scores` and `composition_matrix_adjust`
// up to the `optimize_target_frequencies` call, with the same helpers, arguments and order, returning
// that call's input. `None` where the call would not reach the relative-entropy-constrained Newton
// iteration (a composition without true amino acids, position-based or `CompositionBasedStats`,
// rule `CompoScaleOldMatrix` or an unsupported rule, `desired_re <= 0`, an error). The p-value of the
// pair (`compositionTestIndex`) does not enter the input and is not computed.
// No NCBI counterpart: the input is used only to compute a result ahead of the real call, which uses
// it only when its own input has the same bits; it does not change any value NCBI computes.
/// The Newton input that `blast_adjust_scores` with these arguments would pass to
/// `optimize_target_frequencies`, if it reaches the constrained iteration. Its `q` (the
/// workspace's `mat_b`) is `shared_q`, made by the first call that returns an input: a caller
/// passes one `shared_q` per workspace, so the inputs of a batch share one copy.
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_newton_input(
    matrix_info: &BlastMatrixInfo,
    query_composition: &BlastAminoAcidComposition,
    query_length: i32,
    subject_composition: &BlastAminoAcidComposition,
    subject_length: i32,
    composition_adjust_mode: BlastCompoAdjustMode,
    workspace: &BlastCompositionWorkspace,
    shared_q: &mut Option<super::x_newton_lanes::NewtonQ>,
) -> Option<super::x_newton_lanes::NewtonInput> {
    if query_composition.num_true_amino_acids == 0 || subject_composition.num_true_amino_acids == 0
    {
        return None;
    }
    let mut permuted_query_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut permuted_match_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut permuted_query_probs, &query_composition.prob);
    gather_letter_probs(&mut permuted_match_probs, &subject_composition.prob);
    if matrix_info.position_based
        || composition_adjust_mode == BlastCompoAdjustMode::CompositionBasedStats
    {
        return None;
    }
    let matrix_adjust_rule = blast_choose_matrix_adjust_rule(
        query_length,
        subject_length,
        &permuted_query_probs,
        &permuted_match_probs,
        matrix_info.matrix,
        composition_adjust_mode,
    )
    .ok()?;
    if matrix_adjust_rule == EMatrixAdjustRule::CompoScaleOldMatrix {
        return None;
    }
    // composition_matrix_adjust(matrix_info, rule, query num_true_amino_acids, subject
    // num_true_amino_acids, query prob, subject prob, 20, K_FIXED_RE_BLOSUM62, workspace)
    let mut row_probs = [0.0; COMPO_NUM_TRUE_AA];
    let mut col_probs = [0.0; COMPO_NUM_TRUE_AA];
    gather_letter_probs(&mut row_probs, &query_composition.prob);
    gather_letter_probs(&mut col_probs, &subject_composition.prob);
    let desired_re = match matrix_adjust_rule {
        EMatrixAdjustRule::UserSpecifiedRelEntropy => K_FIXED_RE_BLOSUM62,
        _ => return None,
    };
    #[expect(
        clippy::neg_cmp_op_on_partial_ord,
        reason = "NCBI constrains the relative entropy only when desired_re > 0.0 (composition_adjustment.c:1396); `desired_re <= 0.0` would let NaN through"
    )]
    if !(desired_re > 0.0) {
        return None;
    }
    apply_pseudocounts(
        &mut row_probs,
        query_composition.num_true_amino_acids,
        &workspace.first_standard_freq,
        20,
    );
    apply_pseudocounts(
        &mut col_probs,
        subject_composition.num_true_amino_acids,
        &workspace.second_standard_freq,
        20,
    );
    let q = shared_q.get_or_insert_with(|| {
        let mut q = [0.0; COMPO_NUM_TRUE_AA * COMPO_NUM_TRUE_AA];
        for i in 0..COMPO_NUM_TRUE_AA {
            for j in 0..COMPO_NUM_TRUE_AA {
                q[i * COMPO_NUM_TRUE_AA + j] = workspace.mat_b[i][j];
            }
        }
        std::rc::Rc::new(q)
    });
    Some(super::x_newton_lanes::NewtonInput {
        q: std::rc::Rc::clone(q),
        row: row_probs,
        col: col_probs,
        relative_entropy: desired_re,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1444-1482
fn compute_pair_composition_pvalue(
    permuted_query_probs: &[f64; COMPO_NUM_TRUE_AA],
    permuted_match_probs: &[f64; COMPO_NUM_TRUE_AA],
) -> Result<f64> {
    let mut scores = dense_matrix_new(COMPO_NUM_TRUE_AA, COMPO_NUM_TRUE_AA);
    for (i, row) in scores.iter_mut().enumerate().take(COMPO_NUM_TRUE_AA) {
        for (j, value) in row.iter_mut().enumerate().take(COMPO_NUM_TRUE_AA) {
            *value = protein_score_direct(ScoringMatrix::Blosum62, i, j) as f64;
        }
    }
    let (mut lambda_for_pair, iter_count) = calc_lambda_full_precision_from_scores(
        &scores,
        COMPO_NUM_TRUE_AA,
        permuted_query_probs,
        permuted_match_probs,
    );
    if iter_count >= K_LAMBDA_ITERATION_LIMIT {
        lambda_for_pair = COMPO_MIN_LAMBDA;
    }
    Ok(blast_composition_pvalue(lambda_for_pair))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1105-1155
fn compute_lambda_ratio(
    matrix_info: &BlastMatrixInfo,
    query_prob: &[f64; COMPO_LARGEST_ALPHABET],
    subject_prob: &[f64; COMPO_LARGEST_ALPHABET],
    pvalue_adjustment: bool,
    calc_lambda: fn(&[f64], i32, i32, f64) -> Result<f64>,
) -> f64 {
    let (score_probs, obs_min, obs_max) =
        get_matrix_score_probs(&matrix_info.start_matrix, query_prob, subject_prob);
    let ratio = match calc_lambda(&score_probs, obs_min, obs_max, matrix_info.ungapped_lambda) {
        Ok(lambda) => lambda / matrix_info.ungapped_lambda,
        Err(_) => -1.0,
    };

    let mut lambda_ratio = ratio;
    if !pvalue_adjustment {
        lambda_ratio = lambda_ratio.min(1.0);
    }
    lambda_ratio.max(LAMBDA_RATIO_LOWER_BOUND)
}

fn build_scaled_square_matrix(
    matrix_info: &BlastMatrixInfo,
    row_prob: &[f64; COMPO_LARGEST_ALPHABET],
    col_prob: &[f64; COMPO_LARGEST_ALPHABET],
    lambda_ratio: f64,
) -> AdjustedProteinMatrix {
    let scaled_lambda = matrix_info.ungapped_lambda / lambda_ratio;
    let mut scores = matrix_info.start_freq_ratios;
    freq_ratio_to_score_square(
        &mut scores,
        COMPO_LARGEST_ALPHABET,
        COMPO_LARGEST_ALPHABET,
        scaled_lambda,
    );
    set_xuo_scores(&mut scores, row_prob, col_prob);
    let mut rounded = round_score_matrix(&scores);
    for index in 0..COMPO_LARGEST_ALPHABET {
        rounded[index][25] = matrix_info.start_matrix[index][25];
        rounded[25][index] = matrix_info.start_matrix[25][index];
    }
    AdjustedProteinMatrix { scores: rounded }
}

fn get_matrix_score_probs(
    matrix: &[[i32; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    query_prob: &[f64; COMPO_LARGEST_ALPHABET],
    subject_prob: &[f64; COMPO_LARGEST_ALPHABET],
) -> (Vec<f64>, i32, i32) {
    let mut obs_min = 0;
    let mut obs_max = 0;
    for row in 0..COMPO_LARGEST_ALPHABET {
        for &col in &TRUE_CHAR_POSITIONS {
            let score = matrix[row][col];
            if score > obs_max {
                obs_max = score;
            }
            if score < obs_min && score > COMPO_SCORE_MIN {
                obs_min = score;
            }
        }
    }

    let mut score_prob = vec![0.0; (obs_max - obs_min + 1) as usize];
    for row in 0..COMPO_LARGEST_ALPHABET {
        for &col in &TRUE_CHAR_POSITIONS {
            let score = matrix[row][col];
            if score >= obs_min {
                score_prob[(score - obs_min) as usize] += query_prob[row] * subject_prob[col];
            }
        }
    }
    (score_prob, obs_min, obs_max)
}

fn calc_freq_ratios(
    ratios: &mut [[f64; COMPO_LARGEST_ALPHABET]; COMPO_LARGEST_ALPHABET],
    row_prob: &[f64; COMPO_LARGEST_ALPHABET],
    col_prob: &[f64; COMPO_LARGEST_ALPHABET],
) {
    for row in 0..COMPO_LARGEST_ALPHABET {
        for col in 0..COMPO_LARGEST_ALPHABET {
            let product = row_prob[row] * col_prob[col];
            ratios[row][col] = if product > 0.0 {
                ratios[row][col] / product
            } else {
                0.0
            };
        }
    }
}

fn high_pair_frequencies(letter_probs: &[f64; COMPO_NUM_TRUE_AA], length: i32) -> bool {
    if length <= LENGTH_LOWER_THRESHOLD {
        return false;
    }
    let mut max = 0.0;
    let mut second = 0.0;
    for &prob in letter_probs {
        if prob > second {
            second = prob;
            if prob > max {
                second = max;
                max = prob;
            }
        }
    }
    (max + second) > HIGH_PAIR_THRESHOLD
}

fn high_pair_either_seq(
    query_prob: &[f64; COMPO_NUM_TRUE_AA],
    query_length: i32,
    match_prob: &[f64; COMPO_NUM_TRUE_AA],
    match_length: i32,
) -> bool {
    high_pair_frequencies(query_prob, query_length)
        || high_pair_frequencies(match_prob, match_length)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_mode_condition.c:200-204
// `angle > ANGLE_DEGREE_THRESHOLD` (compo_mode_condition.c:221) for this angle, decided from the
// cosine without a libm: glibc's `acos(cosine) * HALF_CIRCLE_DEGREES / PI` (on x86-64 with FMA,
// where NCBI wrote the tests' bytes) is above 70 exactly when the cosine is in
// [-1, ANGLE_COSINE_LIMIT). `acos` decreases; at the limit the expression is 70.0 and one double
// below it 70.00000000000001. Every double within 2^24 steps of the limit and 10^7 random cosines
// in [-1, 1] give the same decision with glibc 2.35 and 2.39 (with and without their FMA
// variants) and with musl (WASI), which also all put the limit at this double; another platform's
// `acos` could move it by a step (`utils::glibc_fma_libm`).
// No clamp and no zero guard: a cosine rounded below -1 or above 1, or a zero distance (NaN), is
// outside the range, where NCBI's `acos` gives NaN and `angle > ANGLE_DEGREE_THRESHOLD` is false
// (the relative-entropy rule).
const ANGLE_COSINE_LIMIT: f64 = f64::from_bits(0x3fd5_e3a8_7442_b4c9); // 0.34202014306613343
                                                                       // The limit is for these constants.
const _: () =
    assert!(HALF_CIRCLE_DEGREES == 180.0 && PI == 3.1415926543 && ANGLE_DEGREE_THRESHOLD == 70.0);

fn matrix_adjust_rule_angle_above_threshold(d_m_mat: f64, d_q_mat: f64, d_m_q: f64) -> bool {
    let cosine = (d_m_mat * d_m_mat + d_q_mat * d_q_mat - d_m_q * d_m_q) / 2.0 / d_m_mat / d_q_mat;
    (-1.0..ANGLE_COSINE_LIMIT).contains(&cosine)
}

fn blast_choose_matrix_adjust_rule(
    length1: i32,
    length2: i32,
    prob_array1: &[f64; COMPO_NUM_TRUE_AA],
    prob_array2: &[f64; COMPO_NUM_TRUE_AA],
    matrix: ScoringMatrix,
    composition_adjust_mode: BlastCompoAdjustMode,
) -> Result<EMatrixAdjustRule> {
    match composition_adjust_mode {
        BlastCompoAdjustMode::NoCompositionBasedStats => Ok(EMatrixAdjustRule::DontAdjustMatrix),
        BlastCompoAdjustMode::CompositionBasedStats => Ok(EMatrixAdjustRule::CompoScaleOldMatrix),
        BlastCompoAdjustMode::ForceFullMatrixAdjust => {
            Ok(EMatrixAdjustRule::UserSpecifiedRelEntropy)
        }
        BlastCompoAdjustMode::CompositionMatrixAdjust => {
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_mode_condition.c:188
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/matrix_frequency_data.c:1386-1390
            let Some(data) = compo_matrix_data(matrix) else {
                bail!(
                    "matrix {} is not supported for RE based adjustment",
                    matrix.ncbi_name()
                );
            };
            let p_matrix = data.background;

            let mut corr_factor = 0.0;
            for i in 0..COMPO_NUM_TRUE_AA {
                corr_factor += (prob_array1[i] - p_matrix[i]) * (prob_array2[i] - p_matrix[i]);
            }
            let d_m_mat = blast_get_relative_entropy(prob_array2, p_matrix);
            let d_q_mat = blast_get_relative_entropy(prob_array1, p_matrix);
            let d_m_q = blast_get_relative_entropy(prob_array2, prob_array1);

            let angle_above_threshold =
                matrix_adjust_rule_angle_above_threshold(d_m_mat, d_q_mat, d_m_q);

            let len_large = (length1.max(length2)) as f64;
            let len_small = (length1.min(length2)) as f64;
            let _ = corr_factor;
            // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_mode_condition.c:215-225
            // The length ratio is not guarded either: a zero length divides to inf or NaN as in C
            // (and `blast_adjust_scores` returns before a sequence without true amino acids).
            if high_pair_either_seq(prob_array1, length1, prob_array2, length2) {
                Ok(EMatrixAdjustRule::UserSpecifiedRelEntropy)
            } else if d_m_q > QUERY_MATCH_DISTANCE_THRESHOLD
                && len_large / len_small > LENGTH_RATIO_THRESHOLD
                && angle_above_threshold
            {
                Ok(EMatrixAdjustRule::CompoScaleOldMatrix)
            } else {
                Ok(EMatrixAdjustRule::UserSpecifiedRelEntropy)
            }
        }
    }
}

#[cfg(test)]
#[path = "adjust_scores_tests.rs"]
mod tests;
