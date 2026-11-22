//! NCBI genetic-code tables used by translated searches.

use std::sync::OnceLock;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:588-603
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/gc.prt:105-357
// Each row below copies its pinned ncbieaa entry.
const TABLES: &[(u8, &[u8; 64])] = &[
    // NCBI gc.prt:109: id 1, ncbieaa "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        1,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:119: id 2, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
    (
        2,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:129: id 3, ncbieaa "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        3,
        b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:140: id 4, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        4,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:150: id 5, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
    (
        5,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:160: id 6, ncbieaa "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        6,
        b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:170: id 9, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    (
        9,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:180: id 10, ncbieaa "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        10,
        b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:189: id 11, ncbieaa "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        11,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:198: id 12, ncbieaa "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        12,
        b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:207: id 13, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
    (
        13,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:216: id 14, ncbieaa "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    (
        14,
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:225: id 15, ncbieaa "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        15,
        b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:234: id 16, ncbieaa "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        16,
        b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:243: id 21, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    (
        21,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:252: id 22, ncbieaa "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        22,
        b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:261: id 23, ncbieaa "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        23,
        b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:270: id 24, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
    (
        24,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:279: id 25, ncbieaa "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        25,
        b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:288: id 26, ncbieaa "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        26,
        b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:297: id 27, ncbieaa "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        27,
        b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:306: id 28, ncbieaa "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        28,
        b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:315: id 29, ncbieaa "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        29,
        b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:324: id 30, ncbieaa "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        30,
        b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:333: id 31, ncbieaa "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        31,
        b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:342: id 32, ncbieaa "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        32,
        b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:351: id 33, ncbieaa "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
    (
        33,
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    ),
];

pub struct GeneticCode {
    pub table: [u8; 64],
    /// The code's slot of `DISPLAY_TABLES`.
    display: &'static OnceLock<Box<[u8; 4096]>>,
}

// No NCBI counterpart: one slot per row of `TABLES` for the display's translation table of that
// code (`tblastx::report::display_table`; graduated `LOSAT_X_DISPTABLE`), made at its first use and
// kept for the process, so a code is tabulated once however many `GeneticCode` values name it.
static DISPLAY_TABLES: [OnceLock<Box<[u8; 4096]>>; TABLES.len()] =
    [const { OnceLock::new() }; TABLES.len()];

impl GeneticCode {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:588-600
    // A missing ID is an error at this boundary, never an implicit code-1 table.
    pub fn try_from_id(id: u8) -> Result<Self, String> {
        TABLES
            .iter()
            .position(|(code, _)| *code == id)
            .map(|row| Self {
                table: *TABLES[row].1,
                display: &DISPLAY_TABLES[row],
            })
            .ok_or_else(|| format!("unsupported genetic code {id}"))
    }

    /// No NCBI counterpart: the slot that holds this code's display translation table
    /// (`tblastx::report::display_table`).
    pub(crate) fn display_slot(&self) -> &'static OnceLock<Box<[u8; 4096]>> {
        self.display
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/gencode_singleton.c:65-69
    // Existing internal TBLASTX callers pass CLI-validated IDs; invalid direct calls fail.
    pub fn from_id(id: u8) -> Self {
        Self::try_from_id(id).unwrap_or_else(|error| panic!("{error}"))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:369-424
    pub fn get(&self, codon: &[u8]) -> u8 {
        if codon.len() != 3 {
            return b'X';
        }
        // Fast path (was LOSAT_X_CODONFAST): three unambiguous bases select exactly
        // one (i, j, k) of the loops below, so the answer is that one table entry.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:375-379,394-403
        // The path below ports this nested loop, for codons with an ambiguous base. When all three
        // bases are unambiguous, only one (i, j, k) passes the `&` tests, so the loop would end
        // with `aa` taken from `codes` at index 16 times `i`, plus 4 times `j`, plus `k`, and the table
        // is read once instead.
        if let Some(aa) = self.x_get_unambiguous(codon) {
            return aa;
        }
        let masks = [
            base_mask(codon[0]),
            base_mask(codon[1]),
            base_mask(codon[2]),
        ];
        if masks.contains(&0) {
            return b'X';
        }
        let bit = [8, 2, 1, 4];
        let mut aa = 0;
        for i in 0..4 {
            if masks[0] & bit[i] == 0 {
                continue;
            }
            for j in 0..4 {
                if masks[1] & bit[j] == 0 {
                    continue;
                }
                for k in 0..4 {
                    if masks[2] & bit[k] == 0 {
                        continue;
                    }
                    let translated = self.table[i * 16 + j * 4 + k];
                    if aa == 0 {
                        aa = translated;
                    } else if translated != aa {
                        return b'X';
                    }
                }
            }
        }
        aa
    }
}

impl GeneticCode {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:375-379,394-403
    /// This is the one (i, j, k) that `s_CodonToAA` reaches when every codon byte has a single
    /// bit set: `mapping[]` gives the index of that bit. Ambiguous bases return `None` and use
    /// the loop in `get`.
    /// Fast path of `get` (was LOSAT_X_CODONFAST): the amino acid of three unambiguous
    /// bases, `None` when a base is ambiguous or not a base.
    #[inline]
    fn x_get_unambiguous(&self, codon: &[u8]) -> Option<u8> {
        let i = X_BASE_INDEX[codon[0] as usize];
        let j = X_BASE_INDEX[codon[1] as usize];
        let k = X_BASE_INDEX[codon[2] as usize];
        ((i | j | k) < 4).then(|| self.table[i as usize * 16 + j as usize * 4 + k as usize])
    }
}

// Position of a base in NCBI's
// four-entry `mapping` array of 8, 2, 1, 4 (T, C, A, G) when `base_mask` has exactly that
// one bit, 4 otherwise.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:375-379,395-396
// The table holds i (0..3) for T/U, C, A, G in either case. Any other byte, including
// an ambiguity code, maps to 4 so that `(i | j | k) < 4` fails.
const X_BASE_INDEX: [u8; 256] = {
    let mut table = [4u8; 256];
    let bases = [(b'T', 0u8), (b'U', 0), (b'C', 1), (b'A', 2), (b'G', 3)];
    let mut n = 0;
    while n < bases.len() {
        let (base, index) = bases[n];
        table[base as usize] = index;
        table[base.to_ascii_lowercase() as usize] = index;
        n += 1;
    }
    table
};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:94-103
pub(crate) fn base_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'M' => 3,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        b'N' => 15,
        _ => 0,
    }
}

/// `NCBI4NA_TO_IUPACNA`: the letter of each ncbi4na code.
const NCBI4NA_LETTERS: &[u8; 16] = b"-ACMGRSVTWYHKDBN";

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:812-819
const NCBI4NA_COMPLEMENT: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

/// NCBI's `conversion_table` on the letters of the ncbi4na codes, in either case. Any other
/// byte stays: `U`, which the reader stores as `T` before a search translates a sequence
/// (`fasta_reader`), has no code here. LOSAT keeps the case of a letter because a lowercase
/// letter is a masked residue (`FastaRecord::sequence`).
static COMPLEMENT: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut byte = 0;
    while byte < 256 {
        table[byte] = byte as u8;
        byte += 1;
    }
    let mut code = 0;
    while code < 16 {
        let letter = NCBI4NA_LETTERS[code];
        let complement = NCBI4NA_LETTERS[NCBI4NA_COMPLEMENT[code] as usize];
        table[letter as usize] = complement;
        table[letter.to_ascii_lowercase() as usize] = complement.to_ascii_lowercase();
        code += 1;
    }
    table
};

/// The complement of one nucleotide letter (`COMPLEMENT`).
#[inline]
pub(crate) fn complement(base: u8) -> u8 {
    COMPLEMENT[base as usize]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:828-833
/// The letters of the reverse strand, last to first, without NCBI's sentinels (the callers
/// add their own) and its `FENCE_SENTRY` (a letter here is never one).
pub(crate) fn revcomp(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&base| complement(base))
        .collect()
}

#[cfg(test)]
#[path = "genetic_code_tests.rs"]
mod tests;
