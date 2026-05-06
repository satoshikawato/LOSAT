//! Standard protein scoring matrices - NCBI BLAST compatible
//!
//! NCBI uses two representations for protein scoring:
//! 1. NCBISTDAA sequence encoding (`blast_encoding.c`)
//! 2. Standard matrix order `ARNDCQEGHILKMFPSTWYVBJZX*` (`sm_*.c`)
//!
//! This module preserves the existing BLOSUM62 helpers and adds the other
//! standard BLAST protein matrices from NCBI's packed score-matrix tables.

use std::sync::OnceLock;

use crate::config::ScoringMatrix;

pub const BLASTAA_SIZE: usize = 28;
pub const PROTEIN_MATRIX_SIZE: usize = 25;
pub const BLOSUM62_SIZE: usize = PROTEIN_MATRIX_SIZE;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_stat.h:121-122
// ```c
// #define BLAST_SCORE_MIN INT2_MIN
// #define BLAST_SCORE_MAX INT2_MAX
// ```
pub const BLAST_SCORE_MIN: i32 = i16::MIN as i32;

/// NCBISTDAA encoding (0-27).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:115-121
/// ```c
/// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
/// '-','A','B','C','D','E','F','G','H','I','K','L','M',
/// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'};
/// ```
pub mod ncbistdaa {
    pub const GAP: u8 = 0;
    pub const A: u8 = 1;
    pub const B: u8 = 2;
    pub const C: u8 = 3;
    pub const D: u8 = 4;
    pub const E: u8 = 5;
    pub const F: u8 = 6;
    pub const G: u8 = 7;
    pub const H: u8 = 8;
    pub const I: u8 = 9;
    pub const K: u8 = 10;
    pub const L: u8 = 11;
    pub const M: u8 = 12;
    pub const N: u8 = 13;
    pub const P: u8 = 14;
    pub const Q: u8 = 15;
    pub const R: u8 = 16;
    pub const S: u8 = 17;
    pub const T: u8 = 18;
    pub const V: u8 = 19;
    pub const W: u8 = 20;
    pub const X: u8 = 21;
    pub const Y: u8 = 22;
    pub const Z: u8 = 23;
    pub const U: u8 = 24;
    pub const STOP: u8 = 25;
    pub const O: u8 = 26;
    pub const J: u8 = 27;
}

/// Standard matrix order `ARNDCQEGHILKMFPSTWYVBJZX*`.
pub mod blosum62_order {
    pub const A: u8 = 0;
    pub const R: u8 = 1;
    pub const N: u8 = 2;
    pub const D: u8 = 3;
    pub const C: u8 = 4;
    pub const Q: u8 = 5;
    pub const E: u8 = 6;
    pub const G: u8 = 7;
    pub const H: u8 = 8;
    pub const I: u8 = 9;
    pub const L: u8 = 10;
    pub const K: u8 = 11;
    pub const M: u8 = 12;
    pub const F: u8 = 13;
    pub const P: u8 = 14;
    pub const S: u8 = 15;
    pub const T: u8 = 16;
    pub const W: u8 = 17;
    pub const Y: u8 = 18;
    pub const V: u8 = 19;
    pub const B: u8 = 20;
    pub const J: u8 = 21;
    pub const Z: u8 = 22;
    pub const X: u8 = 23;
    pub const STOP: u8 = 24;
}

/// Convert NCBISTDAA index (0-27) to the packed 25x25 matrix index used by
/// NCBI's standard protein score matrices.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1566-1592
/// ```c
/// if (i == AMINOACID_TO_NCBISTDAA['U'] || i == AMINOACID_TO_NCBISTDAA['O'] ||
///     i == AMINOACID_TO_NCBISTDAA['-'] || ...)
///     continue;
/// ...
/// matrix[u_index][i] = matrix[c_index][i];
/// ...
/// matrix[o_index][i] = matrix[x_index][i];
/// ```
#[inline(always)]
pub fn ncbistdaa_to_blosum62(ncbi: u8) -> u8 {
    const TABLE: [u8; 28] = [
        23, 0, 20, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 23, 18, 22, 23,
        24, 23, 21,
    ];
    if ncbi < 28 {
        match ncbi {
            ncbistdaa::U => blosum62_order::C,
            ncbistdaa::O => blosum62_order::X,
            _ => TABLE[ncbi as usize],
        }
    } else {
        23
    }
}

/// Convert ASCII amino acid character to NCBISTDAA index (0-27).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:110-121
/// ```c
/// ... AMINOACID_TO_NCBISTDAA ...
/// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
/// '-','A','B','C','D','E','F','G','H','I','K','L','M',
/// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'};
/// ```
#[inline(always)]
pub fn aa_char_to_ncbistdaa(aa: u8) -> u8 {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:110-121
    // ```c
    // ... AMINOACID_TO_NCBISTDAA ...
    // const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
    // '-','A','B','C','D','E','F','G','H','I','K','L','M',
    // 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'};
    // ```
    const TABLE: [u8; 128] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 27, 10, 11, 12, 13, 26, 14, 15, 16, 17, 18, 24,
        19, 20, 21, 22, 23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    if aa < 128 {
        TABLE[aa as usize]
    } else {
        21
    }
}

// NCBI reference: ncbi-blast/c++/src/util/tables/sm_blosum62.c:34-92
// ```c
// static const TNCBIScore s_Blosum62PSM[25 * 25] = {
// ...
// };
// ```
pub static BLOSUM62: [i8; BLOSUM62_SIZE * BLOSUM62_SIZE] = [
    /*A*/ 4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1,
    -1, -1, -4, /*R*/ -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2,
    -3, -1, -2, 0, -1, -4, /*N*/ -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4,
    -2, -3, 4, -3, 0, -1, -4, /*D*/ -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0,
    -1, -4, -3, -3, 4, -3, 1, -1, -4, /*C*/ 0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1,
    -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4, /*Q*/ -1, 1, 0, 0, -3, 5, 2, -2, 0, -3,
    -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, -2, 4, -1, -4, /*E*/ -1, 0, 0, 2, -4, 2, 5, -2, 0,
    -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, -3, 4, -1, -4, /*G*/ 0, -2, 0, -1, -3, -2,
    -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -4, /*H*/ -2, 0, 1,
    -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, -3, 0, -1, -4,
    /*I*/ -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, 3,
    -3, -1, -4, /*L*/ -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1,
    1, -4, 3, -3, -1, -4, /*K*/ -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1,
    -3, -2, -2, 0, -3, 1, -1, -4, /*M*/ -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2,
    -1, -1, -1, -1, 1, -3, 2, -1, -1, -4, /*F*/ -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3,
    0, 6, -4, -2, -2, 1, 3, -1, -3, 0, -3, -1, -4, /*P*/ -1, -2, -2, -1, -3, -1, -1, -2, -2,
    -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4, /*S*/ 1, -1, 1, 0, -1, 0,
    0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, -2, 0, -1, -4, /*T*/ 0, -1, 0, -1,
    -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1, -1, -4, /*W*/ -3,
    -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -2, -2, -1, -4,
    /*Y*/ -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -1,
    -2, -1, -4, /*V*/ 0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1,
    4, -3, 2, -2, -1, -4, /*B*/ -2, -1, 4, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1,
    -4, -3, -3, 4, -3, 0, -1, -4, /*J*/ -1, -2, -3, -3, -1, -2, -3, -4, -3, 3, 3, -3, 2, 0,
    -3, -2, -1, -2, -1, 2, -3, 3, -3, -1, -4, /*Z*/ -1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1,
    -1, -3, -1, 0, -1, -2, -2, -2, 0, -3, 4, -1, -4, /*X*/ -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4, /***/ -4, -4, -4, -4, -4,
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1,
];

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:110-121
// ```c
// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
// '-','A','B','C','D','E','F','G','H','I','K','L','M',
// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
// 'O', 'J'};
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1559-1592
// ```c
// matrix[i][j] = BLAST_SCORE_MIN;
// ...
// matrix[u_index][i] = matrix[c_index][i];
// matrix[o_index][i] = matrix[x_index][i];
// ```
//
// NCBI reference: ncbi-blast/c++/src/util/tables/sm_blosum62.c:34-92
// ```c
// static const TNCBIScore s_Blosum62PSM[25 * 25] = {
//     /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
//              F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */
// ```
const NCBISTDAA_TO_BLOSUM62_DIRECT: [usize; BLASTAA_SIZE] = [
    blosum62_order::X as usize,
    blosum62_order::A as usize,
    blosum62_order::B as usize,
    blosum62_order::C as usize,
    blosum62_order::D as usize,
    blosum62_order::E as usize,
    blosum62_order::F as usize,
    blosum62_order::G as usize,
    blosum62_order::H as usize,
    blosum62_order::I as usize,
    blosum62_order::K as usize,
    blosum62_order::L as usize,
    blosum62_order::M as usize,
    blosum62_order::N as usize,
    blosum62_order::P as usize,
    blosum62_order::Q as usize,
    blosum62_order::R as usize,
    blosum62_order::S as usize,
    blosum62_order::T as usize,
    blosum62_order::V as usize,
    blosum62_order::W as usize,
    blosum62_order::X as usize,
    blosum62_order::Y as usize,
    blosum62_order::Z as usize,
    blosum62_order::C as usize,
    blosum62_order::STOP as usize,
    blosum62_order::X as usize,
    blosum62_order::J as usize,
];

#[inline(always)]
fn ncbistdaa_to_blosum62_direct_index(ncbi: u8) -> usize {
    if ncbi < BLASTAA_SIZE as u8 {
        NCBISTDAA_TO_BLOSUM62_DIRECT[ncbi as usize]
    } else {
        blosum62_order::X as usize
    }
}

#[inline(always)]
pub fn blosum62_score_ncbistdaa_direct(aa1_ncbi: u8, aa2_ncbi: u8) -> i32 {
    if aa1_ncbi == ncbistdaa::GAP || aa2_ncbi == ncbistdaa::GAP {
        return BLAST_SCORE_MIN;
    }
    let b1 = ncbistdaa_to_blosum62_direct_index(aa1_ncbi);
    let b2 = ncbistdaa_to_blosum62_direct_index(aa2_ncbi);
    BLOSUM62[b1 * PROTEIN_MATRIX_SIZE + b2] as i32
}

const BLOSUM45_TEXT: &str = r#"
/*A*/    5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -1, -1, -1, -1, -5,
/*R*/   -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -1, -3,  1, -1, -5,
/*N*/   -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3,  5, -3,  0, -1, -5,
/*D*/   -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3,  6, -3,  1, -1, -5,
/*C*/   -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -2, -3, -1, -5,
/*Q*/   -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3,  0, -2,  4, -1, -5,
/*E*/   -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3,  1, -3,  5, -1, -5,
/*G*/    0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -5,
/*H*/   -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3,  0, -2,  0, -1, -5,
/*I*/   -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -3,  4, -3, -1, -5,
/*L*/   -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -3,  4, -2, -1, -5,
/*K*/   -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2,  0, -3,  1, -1, -5,
/*M*/   -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -2,  2, -1, -1, -5,
/*F*/   -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -3,  1, -3, -1, -5,
/*P*/   -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -2, -3, -1, -1, -5,
/*S*/    1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1,  0, -2,  0, -1, -5,
/*T*/    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0,  0, -1, -1, -1, -5,
/*W*/   -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -4, -2, -2, -1, -5,
/*Y*/   -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -2,  0, -2, -1, -5,
/*V*/    0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -3,  2, -3, -1, -5,
/*B*/   -1, -1,  5,  6, -2,  0,  1, -1,  0, -3, -3,  0, -2, -3, -2,  0,  0, -4, -2, -3,  5, -3,  1, -1, -5,
/*J*/   -1, -3, -3, -3, -2, -2, -3, -4, -2,  4,  4, -3,  2,  1, -3, -2, -1, -2,  0,  2, -3,  4, -2, -1, -5,
/*Z*/   -1,  1,  0,  1, -3,  4,  5, -2,  0, -3, -2,  1, -1, -3, -1,  0, -1, -2, -2, -3,  1, -2,  5, -1, -5,
/*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
/***/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1
"#;

const BLOSUM50_TEXT: &str = r#"
/*A*/    5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -2, -1, -1, -5,
/*R*/   -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1, -3,  0, -1, -5,
/*N*/   -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5, -4,  0, -1, -5,
/*D*/   -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6, -4,  1, -1, -5,
/*C*/   -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -2, -3, -1, -5,
/*Q*/   -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0, -3,  4, -1, -5,
/*E*/   -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1, -3,  5, -1, -5,
/*G*/    0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -4, -2, -1, -5,
/*H*/   -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0, -3,  0, -1, -5,
/*I*/   -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4,  4, -3, -1, -5,
/*L*/   -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4,  4, -3, -1, -5,
/*K*/   -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0, -3,  1, -1, -5,
/*M*/   -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3,  2, -1, -1, -5,
/*F*/   -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4,  1, -4, -1, -5,
/*P*/   -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -3, -1, -1, -5,
/*S*/    1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0, -3,  0, -1, -5,
/*T*/    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -1, -5,
/*W*/   -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -2, -1, -5,
/*Y*/   -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -1, -2, -1, -5,
/*V*/    0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3,  2, -3, -1, -5,
/*B*/   -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6, -4,  1, -1, -5,
/*J*/   -2, -3, -4, -4, -2, -3, -3, -4, -3,  4,  4, -3,  2,  1, -3, -3, -1, -2, -1,  2, -4,  4, -3, -1, -5,
/*Z*/   -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1, -3,  5, -1, -5,
/*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
/***/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1
"#;

const BLOSUM80_TEXT: &str = r#"
/*A*/    5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -2, -1, -1, -6,
/*R*/   -2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -3,  0, -1, -6,
/*N*/   -2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -4, -3, -4,  5, -4,  0, -1, -6,
/*D*/   -2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4,  5, -5,  1, -1, -6,
/*C*/   -1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -4, -2, -4, -1, -6,
/*Q*/   -1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,  0, -4, -2,  0, -1, -3, -2, -3,  0, -3,  4, -1, -6,
/*E*/   -1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1, -2, -4, -2,  0, -1, -4, -3, -3,  1, -4,  5, -1, -6,
/*G*/    0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -5, -3, -1, -6,
/*H*/   -2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1, -2, -2, -3, -1, -2, -3,  2, -4, -1, -4,  0, -1, -6,
/*I*/   -2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -4,  3, -4, -1, -6,
/*L*/   -2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,  2,  0, -3, -3, -2, -2, -2,  1, -4,  3, -3, -1, -6,
/*K*/   -1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5, -2, -4, -1, -1, -1, -4, -3, -3, -1, -3,  1, -1, -6,
/*M*/   -1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,  6,  0, -3, -2, -1, -2, -2,  1, -3,  2, -1, -1, -6,
/*F*/   -3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,  0,  6, -4, -3, -2,  0,  3, -1, -4,  0, -4, -1, -6,
/*P*/   -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,  8, -1, -2, -5, -4, -3, -2, -4, -2, -1, -6,
/*S*/    1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2,  0, -3,  0, -1, -6,
/*T*/    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2,  1,  5, -4, -2,  0, -1, -1, -1, -1, -6,
/*W*/   -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2,  0, -5, -4, -4, 11,  2, -3, -5, -3, -3, -1, -6,
/*Y*/   -2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -3, -2, -3, -1, -6,
/*V*/    0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,  1, -1, -3, -2,  0, -3, -2,  4, -4,  2, -3, -1, -6,
/*B*/   -2, -1,  5,  5, -4,  0,  1, -1, -1, -4, -4, -1, -3, -4, -2,  0, -1, -5, -3, -4,  5, -4,  0, -1, -6,
/*J*/   -2, -3, -4, -5, -2, -3, -4, -5, -4,  3,  3, -3,  2,  0, -4, -3, -1, -3, -2,  2, -4,  3, -3, -1, -6,
/*Z*/   -1,  0,  0,  1, -4,  4,  5, -3,  0, -4, -3,  1, -1, -4, -2,  0, -1, -3, -3, -3,  0, -3,  5, -1, -6,
/*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6,
/***/   -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1
"#;

const BLOSUM90_TEXT: &str = r#"
/*A*/    5, -2, -2, -3, -1, -1, -1,  0, -2, -2, -2, -1, -2, -3, -1,  1,  0, -4, -3, -1, -2, -2, -1, -1, -6,
/*R*/   -2,  6, -1, -3, -5,  1, -1, -3,  0, -4, -3,  2, -2, -4, -3, -1, -2, -4, -3, -3, -2, -3,  0, -1, -6,
/*N*/   -2, -1,  7,  1, -4,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -5, -3, -4,  5, -4, -1, -1, -6,
/*D*/   -3, -3,  1,  7, -5, -1,  1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5,  5, -5,  1, -1, -6,
/*C*/   -1, -5, -4, -5,  9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -4, -2, -5, -1, -6,
/*Q*/   -1,  1,  0, -1, -4,  7,  2, -3,  1, -4, -3,  1,  0, -4, -2, -1, -1, -3, -3, -3, -1, -3,  5, -1, -6,
/*E*/   -1, -1, -1,  1, -6,  2,  6, -3, -1, -4, -4,  0, -3, -5, -2, -1, -1, -5, -4, -3,  1, -4,  5, -1, -6,
/*G*/    0, -3, -1, -2, -4, -3, -3,  6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -5, -3, -1, -6,
/*H*/   -2,  0,  0, -2, -5,  1, -1, -3,  8, -4, -4, -1, -3, -2, -3, -2, -2, -3,  1, -4, -1, -4,  0, -1, -6,
/*I*/   -2, -4, -4, -5, -2, -4, -4, -5, -4,  5,  1, -4,  1, -1, -4, -3, -1, -4, -2,  3, -5,  3, -4, -1, -6,
/*L*/   -2, -3, -4, -5, -2, -3, -4, -5, -4,  1,  5, -3,  2,  0, -4, -3, -2, -3, -2,  0, -5,  4, -4, -1, -6,
/*K*/   -1,  2,  0, -1, -4,  1,  0, -2, -1, -4, -3,  6, -2, -4, -2, -1, -1, -5, -3, -3, -1, -3,  1, -1, -6,
/*M*/   -2, -2, -3, -4, -2,  0, -3, -4, -3,  1,  2, -2,  7, -1, -3, -2, -1, -2, -2,  0, -4,  2, -2, -1, -6,
/*F*/   -3, -4, -4, -5, -3, -4, -5, -5, -2, -1,  0, -4, -1,  7, -4, -3, -3,  0,  3, -2, -4,  0, -4, -1, -6,
/*P*/   -1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4,  8, -2, -2, -5, -4, -3, -3, -4, -2, -1, -6,
/*S*/    1, -1,  0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2,  5,  1, -4, -3, -2,  0, -3, -1, -1, -6,
/*T*/    0, -2,  0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2,  1,  6, -4, -2, -1, -1, -2, -1, -1, -6,
/*W*/   -4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2,  0, -5, -4, -4, 11,  2, -3, -6, -3, -4, -1, -6,
/*Y*/   -3, -3, -3, -4, -4, -3, -4, -5,  1, -2, -2, -3, -2,  3, -4, -3, -2,  2,  8, -3, -4, -2, -3, -1, -6,
/*V*/   -1, -3, -4, -5, -2, -3, -3, -5, -4,  3,  0, -3,  0, -2, -3, -2, -1, -3, -3,  5, -4,  1, -3, -1, -6,
/*B*/   -2, -2,  5,  5, -4, -1,  1, -2, -1, -5, -5, -1, -4, -4, -3,  0, -1, -6, -4, -4,  5, -5,  0, -1, -6,
/*J*/   -2, -3, -4, -5, -2, -3, -4, -5, -4,  3,  4, -3,  2,  0, -4, -3, -2, -3, -2,  1, -5,  4, -4, -1, -6,
/*Z*/   -1,  0, -1,  1, -5,  5,  5, -3,  0, -4, -4,  1, -2, -4, -2, -1, -1, -4, -3, -3,  0, -4,  5, -1, -6,
/*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6,
/***/   -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1
"#;

const PAM30_TEXT: &str = r#"
/*A*/    6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2,  0, -1,-13, -8, -2, -3, -6, -3, -1,-17,
/*R*/   -7,  8, -6,-10, -8, -2, -9, -9, -2, -5, -8,  0, -4, -9, -4, -3, -6, -2,-10, -8, -7, -7, -4, -1,-17,
/*N*/   -4, -6,  8,  2,-11, -3, -2, -3,  0, -5, -7, -1, -9, -9, -6,  0, -2, -8, -4, -8,  6, -6, -3, -1,-17,
/*D*/   -3,-10,  2,  8,-14, -2,  2, -3, -4, -7,-12, -4,-11,-15, -8, -4, -5,-15,-11, -8,  6,-10,  1, -1,-17,
/*C*/   -6, -8,-11,-14, 10,-14,-14, -9, -7, -6,-15,-14,-13,-13, -8, -3, -8,-15, -4, -6,-12, -9,-14, -1,-17,
/*Q*/   -4, -2, -3, -2,-14,  8,  1, -7,  1, -8, -5, -3, -4,-13, -3, -5, -5,-13,-12, -7, -3, -5,  6, -1,-17,
/*E*/   -2, -9, -2,  2,-14,  1,  8, -4, -5, -5, -9, -4, -7,-14, -5, -4, -6,-17, -8, -6,  1, -7,  6, -1,-17,
/*G*/   -2, -9, -3, -3, -9, -7, -4,  6, -9,-11,-10, -7, -8, -9, -6, -2, -6,-15,-14, -5, -3,-10, -5, -1,-17,
/*H*/   -7, -2,  0, -4, -7,  1, -5, -9,  9, -9, -6, -6,-10, -6, -4, -6, -7, -7, -3, -6, -1, -7, -1, -1,-17,
/*I*/   -5, -5, -5, -7, -6, -8, -5,-11, -9,  8, -1, -6, -1, -2, -8, -7, -2,-14, -6,  2, -6,  5, -6, -1,-17,
/*L*/   -6, -8, -7,-12,-15, -5, -9,-10, -6, -1,  7, -8,  1, -3, -7, -8, -7, -6, -7, -2, -9,  6, -7, -1,-17,
/*K*/   -7,  0, -1, -4,-14, -3, -4, -7, -6, -6, -8,  7, -2,-14, -6, -4, -3,-12, -9, -9, -2, -7, -4, -1,-17,
/*M*/   -5, -4, -9,-11,-13, -4, -7, -8,-10, -1,  1, -2, 11, -4, -8, -5, -4,-13,-11, -1,-10,  0, -5, -1,-17,
/*F*/   -8, -9, -9,-15,-13,-13,-14, -9, -6, -2, -3,-14, -4,  9,-10, -6, -9, -4,  2, -8,-10, -2,-13, -1,-17,
/*P*/   -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -7, -6, -8,-10,  8, -2, -4,-14,-13, -6, -7, -7, -4, -1,-17,
/*S*/    0, -3,  0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5, -6, -2,  6,  0, -5, -7, -6, -1, -8, -5, -1,-17,
/*T*/   -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -7, -3, -4, -9, -4,  0,  7,-13, -6, -3, -3, -5, -6, -1,-17,
/*W*/  -13, -2, -8,-15,-15,-13,-17,-15, -7,-14, -6,-12,-13, -4,-14, -5,-13, 13, -5,-15,-10, -7,-14, -1,-17,
/*Y*/   -8,-10, -4,-11, -4,-12, -8,-14, -3, -6, -7, -9,-11,  2,-13, -7, -6, -5, 10, -7, -6, -7, -9, -1,-17,
/*V*/   -2, -8, -8, -8, -6, -7, -6, -5, -6,  2, -2, -9, -1, -8, -6, -6, -3,-15, -7,  7, -8,  0, -6, -1,-17,
/*B*/   -3, -7,  6,  6,-12, -3,  1, -3, -1, -6, -9, -2,-10,-10, -7, -1, -3,-10, -6, -8,  6, -8,  0, -1,-17,
/*J*/   -6, -7, -6,-10, -9, -5, -7,-10, -7,  5,  6, -7,  0, -2, -7, -8, -5, -7, -7,  0, -8,  6, -6, -1,-17,
/*Z*/   -3, -4, -3,  1,-14,  6,  6, -5, -1, -6, -7, -4, -5,-13, -4, -5, -6,-14, -9, -6,  0, -6,  6, -1,-17,
/*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-17,
/***/  -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,  1
"#;

const PAM70_TEXT: &str = r#"
/*A*/    5, -4, -2, -1, -4, -2, -1,  0, -4, -2, -4, -4, -3, -6,  0,  1,  1, -9, -5, -1, -1, -3, -1, -1,-11,
/*R*/   -4,  8, -3, -6, -5,  0, -5, -6,  0, -3, -6,  2, -2, -7, -2, -1, -4,  0, -7, -5, -4, -5, -2, -1,-11,
/*N*/   -2, -3,  6,  3, -7, -1,  0, -1,  1, -3, -5,  0, -5, -6, -3,  1,  0, -6, -3, -5,  5, -4, -1, -1,-11,
/*D*/   -1, -6,  3,  6, -9,  0,  3, -1, -1, -5, -8, -2, -7,-10, -4, -1, -2,-10, -7, -5,  5, -7,  2, -1,-11,
/*C*/   -4, -5, -7, -9,  9, -9, -9, -6, -5, -4,-10, -9, -9, -8, -5, -1, -5,-11, -2, -4, -8, -7, -9, -1,-11,
/*Q*/   -2,  0, -1,  0, -9,  7,  2, -4,  2, -5, -3, -1, -2, -9, -1, -3, -3, -8, -8, -4, -1, -3,  5, -1,-11,
/*E*/   -1, -5,  0,  3, -9,  2,  6, -2, -2, -4, -6, -2, -4, -9, -3, -2, -3,-11, -6, -4,  2, -5,  5, -1,-11,
/*G*/    0, -6, -1, -1, -6, -4, -2,  6, -6, -6, -7, -5, -6, -7, -3,  0, -3,-10, -9, -3, -1, -7, -3, -1,-11,
/*H*/   -4,  0,  1, -1, -5,  2, -2, -6,  8, -6, -4, -3, -6, -4, -2, -3, -4, -5, -1, -4,  0, -4,  1, -1,-11,
/*I*/   -2, -3, -3, -5, -4, -5, -4, -6, -6,  7,  1, -4,  1,  0, -5, -4, -1, -9, -4,  3, -4,  4, -4, -1,-11,
/*L*/   -4, -6, -5, -8,-10, -3, -6, -7, -4,  1,  6, -5,  2, -1, -5, -6, -4, -4, -4,  0, -6,  5, -4, -1,-11,
/*K*/   -4,  2,  0, -2, -9, -1, -2, -5, -3, -4, -5,  6,  0, -9, -4, -2, -1, -7, -7, -6, -1, -5, -2, -1,-11,
/*M*/   -3, -2, -5, -7, -9, -2, -4, -6, -6,  1,  2,  0, 10, -2, -5, -3, -2, -8, -7,  0, -6,  2, -3, -1,-11,
/*F*/   -6, -7, -6,-10, -8, -9, -9, -7, -4,  0, -1, -9, -2,  8, -7, -4, -6, -2,  4, -5, -7, -1, -9, -1,-11,
/*P*/    0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5, -7,  7,  0, -2, -9, -9, -3, -4, -5, -2, -1,-11,
/*S*/    1, -1,  1, -1, -1, -3, -2,  0, -3, -4, -6, -2, -3, -4,  0,  5,  2, -3, -5, -3,  0, -5, -2, -1,-11,
/*T*/    1, -4,  0, -2, -5, -3, -3, -3, -4, -1, -4, -1, -2, -6, -2,  2,  6, -8, -4, -1, -1, -3, -3, -1,-11,
/*W*/   -9,  0, -6,-10,-11, -8,-11,-10, -5, -9, -4, -7, -8, -2, -9, -3, -8, 13, -3,-10, -7, -5,-10, -1,-11,
/*Y*/   -5, -7, -3, -7, -2, -8, -6, -9, -1, -4, -4, -7, -7,  4, -9, -5, -4, -3,  9, -5, -4, -4, -7, -1,-11,
/*V*/   -1, -5, -5, -5, -4, -4, -4, -3, -4,  3,  0, -6,  0, -5, -3, -3, -1,-10, -5,  6, -5,  1, -4, -1,-11,
/*B*/   -1, -4,  5,  5, -8, -1,  2, -1,  0, -4, -6, -1, -6, -7, -4,  0, -1, -7, -4, -5,  5, -5,  1, -1,-11,
/*J*/   -3, -5, -4, -7, -7, -3, -5, -7, -4,  4,  5, -5,  2, -1, -5, -5, -3, -5, -4,  1, -5,  5, -4, -1,-11,
/*Z*/   -1, -2, -1,  2, -9,  5,  5, -3,  1, -4, -4, -2, -3, -9, -2, -2, -3,-10, -7, -4,  1, -4,  5, -1,-11,
/*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-11,
/***/  -11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,  1
"#;

const PAM250_TEXT: &str = r#"
/*A*/    2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0, -1,  0, -1, -8,
/*R*/   -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1, -3,  0, -1, -8,
/*N*/    0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2, -3,  1, -1, -8,
/*D*/    0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3, -3,  3, -1, -8,
/*C*/   -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -5, -1, -8,
/*Q*/    0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1, -2,  3, -1, -8,
/*E*/    0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3, -3,  3, -1, -8,
/*G*/    1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0, -4,  0, -1, -8,
/*H*/   -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1, -2,  2, -1, -8,
/*I*/   -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2,  3, -2, -1, -8,
/*L*/   -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3,  5, -3, -1, -8,
/*K*/   -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1, -3,  0, -1, -8,
/*M*/   -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2,  3, -2, -1, -8,
/*F*/   -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4,  2, -5, -1, -8,
/*P*/    1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1, -2,  0, -1, -8,
/*S*/    1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0, -2,  0, -1, -8,
/*T*/    1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1, -1, -1, -8,
/*W*/   -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -3, -6, -1, -8,
/*Y*/   -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -1, -4, -1, -8,
/*V*/    0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2,  2, -2, -1, -8,
/*B*/    0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3, -3,  2, -1, -8,
/*J*/   -1, -3, -3, -3, -5, -2, -3, -4, -2,  3,  5, -3,  3,  2, -2, -2, -1, -3, -1,  2, -3,  5, -2, -1, -8,
/*Z*/    0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2, -2,  3, -1, -8,
/*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -8,
/***/   -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1
"#;

// NCBI reference: ncbi-blast/c++/src/util/tables/raw_scoremat.c:81-109
// ```c
// for (i = 0;  i < dim;  ++i) {
//     aa1 = sym[i];
//     for (j = 0;  j < dim;  ++j) {
//         aa2 = sym[j];
//         fsm->s[aa1][aa2] = psm->scores[i * dim + j];
//     }
// }
// ```
fn parse_matrix_scores(raw: &str) -> Vec<i8> {
    let mut values = Vec::with_capacity(PROTEIN_MATRIX_SIZE * PROTEIN_MATRIX_SIZE);
    let mut token = String::new();

    for ch in raw.chars() {
        if ch == '-' || ch.is_ascii_digit() {
            token.push(ch);
        } else if !token.is_empty() {
            values.push(
                token
                    .parse::<i8>()
                    .expect("valid NCBI protein matrix score"),
            );
            token.clear();
        }
    }
    if !token.is_empty() {
        values.push(
            token
                .parse::<i8>()
                .expect("valid NCBI protein matrix score"),
        );
    }

    assert_eq!(
        values.len(),
        PROTEIN_MATRIX_SIZE * PROTEIN_MATRIX_SIZE,
        "matrix text must contain exactly 625 scores"
    );
    values
}

fn blosum45_scores() -> &'static [i8] {
    static SCORES: OnceLock<Vec<i8>> = OnceLock::new();
    SCORES
        .get_or_init(|| parse_matrix_scores(BLOSUM45_TEXT))
        .as_slice()
}

fn blosum50_scores() -> &'static [i8] {
    static SCORES: OnceLock<Vec<i8>> = OnceLock::new();
    SCORES
        .get_or_init(|| parse_matrix_scores(BLOSUM50_TEXT))
        .as_slice()
}

fn blosum80_scores() -> &'static [i8] {
    static SCORES: OnceLock<Vec<i8>> = OnceLock::new();
    SCORES
        .get_or_init(|| parse_matrix_scores(BLOSUM80_TEXT))
        .as_slice()
}

fn blosum90_scores() -> &'static [i8] {
    static SCORES: OnceLock<Vec<i8>> = OnceLock::new();
    SCORES
        .get_or_init(|| parse_matrix_scores(BLOSUM90_TEXT))
        .as_slice()
}

fn pam30_scores() -> &'static [i8] {
    static SCORES: OnceLock<Vec<i8>> = OnceLock::new();
    SCORES
        .get_or_init(|| parse_matrix_scores(PAM30_TEXT))
        .as_slice()
}

fn pam70_scores() -> &'static [i8] {
    static SCORES: OnceLock<Vec<i8>> = OnceLock::new();
    SCORES
        .get_or_init(|| parse_matrix_scores(PAM70_TEXT))
        .as_slice()
}

fn pam250_scores() -> &'static [i8] {
    static SCORES: OnceLock<Vec<i8>> = OnceLock::new();
    SCORES
        .get_or_init(|| parse_matrix_scores(PAM250_TEXT))
        .as_slice()
}

#[derive(Debug, Clone, Copy)]
pub struct ProteinScoreMatrix {
    pub scores: &'static [i8],
    pub defscore: i32,
}

// NCBI reference: ncbi-blast/c++/src/util/tables/raw_scoremat.c:131-154
// ```c
// const SNCBIPackedScoreMatrix* NCBISM_GetStandardMatrix(const char* name)
// {
//     ...
// }
// ```
#[inline(always)]
pub fn standard_protein_matrix(matrix: ScoringMatrix) -> ProteinScoreMatrix {
    match matrix {
        ScoringMatrix::Blosum45 => ProteinScoreMatrix {
            scores: blosum45_scores(),
            defscore: -5,
        },
        ScoringMatrix::Blosum50 => ProteinScoreMatrix {
            scores: blosum50_scores(),
            defscore: -5,
        },
        ScoringMatrix::Blosum62 => ProteinScoreMatrix {
            scores: &BLOSUM62,
            defscore: -4,
        },
        ScoringMatrix::Blosum80 => ProteinScoreMatrix {
            scores: blosum80_scores(),
            defscore: -6,
        },
        ScoringMatrix::Blosum90 => ProteinScoreMatrix {
            scores: blosum90_scores(),
            defscore: -6,
        },
        ScoringMatrix::Pam30 => ProteinScoreMatrix {
            scores: pam30_scores(),
            defscore: -17,
        },
        ScoringMatrix::Pam70 => ProteinScoreMatrix {
            scores: pam70_scores(),
            defscore: -11,
        },
        ScoringMatrix::Pam250 => ProteinScoreMatrix {
            scores: pam250_scores(),
            defscore: -8,
        },
    }
}

#[inline(always)]
pub fn protein_defscore(matrix: ScoringMatrix) -> i32 {
    standard_protein_matrix(matrix).defscore
}

#[inline(always)]
pub fn protein_score(matrix: ScoringMatrix, aa1_ncbi: u8, aa2_ncbi: u8) -> i32 {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1559-1579
    // ```c
    // for (i = 0; i < sbp->alphabet_size; i++) {
    //     for (j = 0; j < sbp->alphabet_size; j++) {
    //         matrix[i][j] = BLAST_SCORE_MIN;
    //     }
    // }
    // ...
    // if (i == AMINOACID_TO_NCBISTDAA['-'] ||
    //     j == AMINOACID_TO_NCBISTDAA['-']) {
    //     continue;
    // }
    // ```
    if aa1_ncbi == ncbistdaa::GAP || aa2_ncbi == ncbistdaa::GAP {
        return BLAST_SCORE_MIN;
    }
    let matrix_data = standard_protein_matrix(matrix);
    let b1 = ncbistdaa_to_blosum62(aa1_ncbi) as usize;
    let b2 = ncbistdaa_to_blosum62(aa2_ncbi) as usize;
    matrix_data.scores[b1 * PROTEIN_MATRIX_SIZE + b2] as i32
}

#[inline(always)]
pub fn protein_score_direct(matrix: ScoringMatrix, b1: usize, b2: usize) -> i32 {
    standard_protein_matrix(matrix).scores[b1 * PROTEIN_MATRIX_SIZE + b2] as i32
}

#[inline(always)]
pub fn blosum62_score(aa1_ncbi: u8, aa2_ncbi: u8) -> i32 {
    blosum62_score_ncbistdaa_direct(aa1_ncbi, aa2_ncbi)
}

#[inline(always)]
pub fn blosum62_score_direct(b1: usize, b2: usize) -> i32 {
    protein_score_direct(ScoringMatrix::Blosum62, b1, b2)
}

pub const AA_TABLE_SIZE: usize = BLOSUM62_SIZE * BLOSUM62_SIZE * BLOSUM62_SIZE;
pub static MATRIX: [i8; BLOSUM62_SIZE * BLOSUM62_SIZE] = BLOSUM62;
pub const NCBI_AA_SIZE: usize = BLASTAA_SIZE;
pub const STOP_CODON: u8 = ncbistdaa::STOP;

#[inline(always)]
pub fn aa_char_to_ncbi_index(aa: u8) -> u8 {
    ncbistdaa_to_blosum62(aa_char_to_ncbistdaa(aa))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ncbistdaa_encoding() {
        assert_eq!(aa_char_to_ncbistdaa(b'A'), 1);
        assert_eq!(aa_char_to_ncbistdaa(b'*'), 25);
        assert_eq!(aa_char_to_ncbistdaa(b'X'), 21);
        assert_eq!(aa_char_to_ncbistdaa(b'J'), 27);
    }

    #[test]
    fn test_standard_matrix_conversion() {
        assert_eq!(ncbistdaa_to_blosum62(1), 0);
        assert_eq!(ncbistdaa_to_blosum62(25), 24);
        assert_eq!(ncbistdaa_to_blosum62(21), 23);
        assert_eq!(ncbistdaa_to_blosum62(ncbistdaa::U), blosum62_order::C);
        assert_eq!(ncbistdaa_to_blosum62(ncbistdaa::O), blosum62_order::X);
    }

    #[test]
    fn test_blosum62_scores() {
        assert_eq!(blosum62_score(ncbistdaa::A, ncbistdaa::A), 4);
        assert_eq!(blosum62_score(ncbistdaa::STOP, ncbistdaa::STOP), 1);
        assert_eq!(blosum62_score(ncbistdaa::X, ncbistdaa::X), -1);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:1559-1592
    // ```c
    // matrix[i][j] = BLAST_SCORE_MIN;
    // ...
    // matrix[u_index][i] = matrix[c_index][i];
    // matrix[o_index][i] = matrix[x_index][i];
    // ```
    #[test]
    fn test_blosum62_ncbistdaa_direct_matches_generic_for_blastaa_alphabet() {
        for q in 0..BLASTAA_SIZE as u8 {
            for s in 0..BLASTAA_SIZE as u8 {
                assert_eq!(
                    blosum62_score_ncbistdaa_direct(q, s),
                    protein_score(ScoringMatrix::Blosum62, q, s),
                    "direct BLOSUM62 score mismatch for NCBISTDAA {q}/{s}"
                );
            }
        }
    }

    #[test]
    fn test_non_blosum62_standard_matrices() {
        assert_eq!(
            protein_score(ScoringMatrix::Blosum45, ncbistdaa::A, ncbistdaa::A),
            5
        );
        assert_eq!(
            protein_score(ScoringMatrix::Blosum80, ncbistdaa::W, ncbistdaa::W),
            11
        );
        assert_eq!(
            protein_score(ScoringMatrix::Pam30, ncbistdaa::A, ncbistdaa::A),
            6
        );
        assert_eq!(
            protein_score(ScoringMatrix::Pam250, ncbistdaa::A, ncbistdaa::A),
            2
        );
    }

    #[test]
    fn test_matrix_defscores_follow_ncbi() {
        assert_eq!(protein_defscore(ScoringMatrix::Blosum45), -5);
        assert_eq!(protein_defscore(ScoringMatrix::Blosum62), -4);
        assert_eq!(protein_defscore(ScoringMatrix::Blosum80), -6);
        assert_eq!(protein_defscore(ScoringMatrix::Pam30), -17);
        assert_eq!(protein_defscore(ScoringMatrix::Pam250), -8);
    }

    #[test]
    fn test_sentinel_uses_blast_score_min() {
        assert_eq!(
            protein_score(ScoringMatrix::Blosum45, ncbistdaa::GAP, ncbistdaa::A),
            BLAST_SCORE_MIN
        );
        assert_eq!(
            protein_score(ScoringMatrix::Pam30, ncbistdaa::A, ncbistdaa::GAP),
            BLAST_SCORE_MIN
        );
        assert_eq!(
            protein_score(ScoringMatrix::Blosum62, ncbistdaa::GAP, ncbistdaa::GAP),
            BLAST_SCORE_MIN
        );
    }

    #[test]
    fn test_u_and_o_follow_ncbi_matrix_aliases() {
        assert_eq!(
            protein_score(ScoringMatrix::Blosum62, ncbistdaa::U, ncbistdaa::A),
            protein_score(ScoringMatrix::Blosum62, ncbistdaa::C, ncbistdaa::A)
        );
        assert_eq!(
            protein_score(ScoringMatrix::Blosum62, ncbistdaa::O, ncbistdaa::W),
            protein_score(ScoringMatrix::Blosum62, ncbistdaa::X, ncbistdaa::W)
        );
    }
}
