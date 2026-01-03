//! BLOSUM62 scoring matrix - NCBI BLAST compatible
//!
//! NCBI uses two different encodings:
//! 1. NCBISTDAA (28 characters) - sequence encoding for lookup/scan
//! 2. BLOSUM62 matrix order (25 characters) - scoring matrix indices
//!
//! This module provides both encodings and conversion between them.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c
//!            ncbi-blast/c++/src/util/tables/sm_blosum62.c

/// NCBI BLASTAA_SIZE - size of amino acid alphabet for lookup tables
pub const BLASTAA_SIZE: usize = 28;

/// Size of BLOSUM62 matrix (25x25)
pub const BLOSUM62_SIZE: usize = 25;

/// Default score for unknown/sentinel residues (NCBI BLOSUM62 defscore)
/// Reference: ncbi-blast/c++/src/util/tables/sm_blosum62.c:95
///   const SNCBIPackedScoreMatrix NCBISM_Blosum62 = { ..., -4 };
pub const DEFSCORE: i32 = -4;

/// NCBISTDAA encoding (0-27)
/// Reference: blast_encoding.c NCBISTDAA_TO_AMINOACID
///   '-','A','B','C','D','E','F','G','H','I','K','L','M',
///   'N','P','Q','R','S','T','V','W','X','Y','Z','U','*','O','J'
pub mod ncbistdaa {
    pub const GAP: u8 = 0;   // '-'
    pub const A: u8 = 1;
    pub const B: u8 = 2;     // Asn or Asp
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
    pub const X: u8 = 21;    // Unknown
    pub const Y: u8 = 22;
    pub const Z: u8 = 23;    // Glu or Gln
    pub const U: u8 = 24;    // Selenocysteine
    pub const STOP: u8 = 25; // '*'
    pub const O: u8 = 26;    // Pyrrolysine
    pub const J: u8 = 27;    // Leu or Ile
}

/// BLOSUM62 matrix order: ARNDCQEGHILKMFPSTWYVBJZX*
/// Reference: sm_blosum62.c
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

/// Convert NCBISTDAA index (0-27) to BLOSUM62 matrix index (0-24)
/// Invalid/gap characters map to X (23)
/// Reference: NCBI uses this conversion internally for scoring
#[inline(always)]
pub fn ncbistdaa_to_blosum62(ncbi: u8) -> u8 {
    const TABLE: [u8; 28] = [
        23, // 0: '-' (gap) -> X
        0,  // 1: A -> 0
        20, // 2: B -> 20
        4,  // 3: C -> 4
        3,  // 4: D -> 3
        6,  // 5: E -> 6
        13, // 6: F -> 13
        7,  // 7: G -> 7
        8,  // 8: H -> 8
        9,  // 9: I -> 9
        11, // 10: K -> 11
        10, // 11: L -> 10
        12, // 12: M -> 12
        2,  // 13: N -> 2
        14, // 14: P -> 14
        5,  // 15: Q -> 5
        1,  // 16: R -> 1
        15, // 17: S -> 15
        16, // 18: T -> 16
        19, // 19: V -> 19
        17, // 20: W -> 17
        23, // 21: X -> 23
        18, // 22: Y -> 18
        22, // 23: Z -> 22
        23, // 24: U (selenocysteine) -> X
        24, // 25: '*' (stop) -> 24
        23, // 26: O (pyrrolysine) -> X
        21, // 27: J -> 21
    ];
    if ncbi < 28 {
        TABLE[ncbi as usize]
    } else {
        23 // Unknown -> X
    }
}

/// Convert ASCII amino acid character to NCBISTDAA index (0-27)
/// Reference: blast_encoding.c AMINOACID_TO_NCBISTDAA
#[inline(always)]
pub fn aa_char_to_ncbistdaa(aa: u8) -> u8 {
    const TABLE: [u8; 128] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 0, // '*' = 25
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 27, 10, 11, 12, 13, 26,
        // A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
        14, 15, 16, 17, 18, 24, 19, 20, 21, 22, 23, 0, 0, 0, 0, 0,
        // P   Q   R   S   T   U   V   W   X   Y   Z
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    if aa < 128 {
        TABLE[aa as usize]
    } else {
        21 // Unknown -> X
    }
}

/// BLOSUM62 matrix in NCBI packed order: ARNDCQEGHILKMFPSTWYVBJZX*
/// Source: NCBI BLAST sm_blosum62.c (verbatim copy)
/// 25 symbols (0-24): A(0), R(1), N(2), D(3), C(4), Q(5), E(6), G(7), H(8), I(9),
///                    L(10), K(11), M(12), F(13), P(14), S(15), T(16), W(17), Y(18), V(19),
///                    B(20), J(21), Z(22), X(23), *(24)
pub static BLOSUM62: [i8; BLOSUM62_SIZE * BLOSUM62_SIZE] = [
    //       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,  F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *
    /*A*/    4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4,
    /*R*/   -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4,
    /*N*/   -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4,
    /*D*/   -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4,
    /*C*/    0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4,
    /*Q*/   -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4,
    /*E*/   -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4,
    /*G*/    0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4,
    /*H*/   -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4,
    /*I*/   -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4,
    /*L*/   -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4,
    /*K*/   -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4,
    /*M*/   -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4,
    /*F*/   -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4,
    /*P*/   -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4,
    /*S*/    1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4,
    /*T*/    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4,
    /*W*/   -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4,
    /*Y*/   -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4,
    /*V*/    0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4,
    /*B*/   -2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4,
    /*J*/   -1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4,
    /*Z*/   -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4,
    /*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4,
    /***/   -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1,
];

/// Get BLOSUM62 score for two amino acids in NCBISTDAA encoding
/// Handles the conversion from NCBISTDAA (28) to BLOSUM62 order (25) internally
/// 
/// NCBI FSM (Full Score Matrix) behavior:
/// - Index 0 (gap/sentinel) returns defscore (-4)
/// - This ensures proper X-drop termination at sequence boundaries
/// 
/// Reference: ncbi-blast/c++/src/util/tables/raw_scoremat.c:91
///   fsm->s[0][i] = psm->defscore;
#[inline(always)]
pub fn blosum62_score(aa1_ncbi: u8, aa2_ncbi: u8) -> i32 {
    // NCBI FSM: index 0 (gap/sentinel) returns defscore
    if aa1_ncbi == 0 || aa2_ncbi == 0 {
        return DEFSCORE;
    }
    let b1 = ncbistdaa_to_blosum62(aa1_ncbi) as usize;
    let b2 = ncbistdaa_to_blosum62(aa2_ncbi) as usize;
    BLOSUM62[b1 * BLOSUM62_SIZE + b2] as i32
}

/// Get BLOSUM62 score for two amino acids already in BLOSUM62 order (0-24)
#[inline(always)]
pub fn blosum62_score_direct(b1: usize, b2: usize) -> i32 {
    BLOSUM62[b1 * BLOSUM62_SIZE + b2] as i32
}

// Legacy compatibility - keep old names as aliases
pub const AA_TABLE_SIZE: usize = BLOSUM62_SIZE * BLOSUM62_SIZE * BLOSUM62_SIZE;
pub static MATRIX: [i8; BLOSUM62_SIZE * BLOSUM62_SIZE] = BLOSUM62;
pub const NCBI_AA_SIZE: usize = BLASTAA_SIZE;
pub const STOP_CODON: u8 = ncbistdaa::STOP;

/// Legacy function - maps ASCII to BLOSUM62 order (for backward compatibility)
/// New code should use aa_char_to_ncbistdaa instead
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
    fn test_blosum62_conversion() {
        // A in NCBISTDAA (1) -> A in BLOSUM62 (0)
        assert_eq!(ncbistdaa_to_blosum62(1), 0);
        // * in NCBISTDAA (25) -> * in BLOSUM62 (24)
        assert_eq!(ncbistdaa_to_blosum62(25), 24);
        // X in NCBISTDAA (21) -> X in BLOSUM62 (23)
        assert_eq!(ncbistdaa_to_blosum62(21), 23);
    }

    #[test]
    fn test_blosum62_scores() {
        // A-A = 4
        assert_eq!(blosum62_score(ncbistdaa::A, ncbistdaa::A), 4);
        // *-* = 1 (NCBI BLOSUM62)
        assert_eq!(blosum62_score(ncbistdaa::STOP, ncbistdaa::STOP), 1);
        // X-X = -1
        assert_eq!(blosum62_score(ncbistdaa::X, ncbistdaa::X), -1);
    }
    
    #[test]
    fn test_sentinel_defscore() {
        // NCBI FSM behavior: index 0 (gap/sentinel) returns defscore (-4)
        // Reference: ncbi-blast/c++/src/util/tables/raw_scoremat.c:91
        //   fsm->s[0][i] = psm->defscore;
        
        // Sentinel (0) with any AA should return defscore
        assert_eq!(blosum62_score(0, ncbistdaa::A), DEFSCORE);
        assert_eq!(blosum62_score(ncbistdaa::A, 0), DEFSCORE);
        assert_eq!(blosum62_score(0, 0), DEFSCORE);
        assert_eq!(blosum62_score(0, ncbistdaa::STOP), DEFSCORE);
        assert_eq!(blosum62_score(0, ncbistdaa::X), DEFSCORE);
        
        // Verify defscore value matches NCBI
        assert_eq!(DEFSCORE, -4);
    }
}
