//! The reading frames of a nucleotide query (`READINGFRAME`), its six translations
//! (`TranslateAll` through the toolkit's `CSeqTranslator::Translate` with the standard genetic
//! code) and its reverse complement (`GetComplementSeq`).
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:167-176
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:28-29
//! The positions are `SeqPos_t` (`int`) and the lengths `SeqLen_t` (`unsigned int`); the
//! functions below follow C++'s conversions: an expression with an `unsigned int` operand is
//! computed modulo 2^32 and stored back in an `int` (two's complement in NCBI's build). A
//! position argument is an `i64`: an `int` that NCBI's arithmetic took past its range holds its
//! mathematical value (PD-LOSAT-NCBI-DEFECTS exception 12), and its conversion to `unsigned int`
//! is that value modulo 2^32.

/// `READINGFRAME::RF_SIZE`.
pub(crate) const RF_SIZE: u32 = 3;

/// `READINGFRAME::RF_IDS`: the frame of each index (`Idx2Id`).
pub(crate) const RF_IDS: [i32; 6] = [1, 2, 3, -1, -2, -3];

/// `PlusRFIdx2MinusRFIdx`: the minus-strand frame index of a plus-strand one.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:200,203
pub(crate) fn plus_rf_idx_to_minus(prfidx: u32, na_len: u32) -> u32 {
    na_len.wrapping_sub(prfidx) % RF_SIZE + RF_SIZE
}

/// `PlusPr2PlusNA`: the plus-strand start of the codon of residue `pr` of a plus-strand frame.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:224-229
pub(crate) fn plus_pr_to_plus_na(pr: i64, rfidx: u32) -> i32 {
    (pr as u32).wrapping_mul(RF_SIZE).wrapping_add(rfidx) as i32
}

/// `MinusPr2PlusNA`: the plus-strand position of the first base of the codon of residue `pr`
/// of a minus-strand frame. `pr - 1` is `int` arithmetic, which overflows at `INT_MIN`
/// (undefined in C++): it is computed as its value (PD-LOSAT-NCBI-DEFECTS exception 12), whose
/// product with `RF_SIZE` converts it to `unsigned int` (modulo 2^32).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:230-238
pub(crate) fn minus_pr_to_plus_na(pr: i64, rfidx: u32, na_len: u32) -> i32 {
    let pr_1 = pr - 1;
    na_len
        .wrapping_sub((pr_1 as u32).wrapping_mul(RF_SIZE))
        .wrapping_sub(rfidx)
        .wrapping_sub(1) as i32
}

/// `PlusId2PlusRFIdx`: the plus-strand frame index of a plus frame (1 to 3; any positive
/// frame, `int` converted to `unsigned int`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:205-209
pub(crate) fn plus_id_to_plus_rf_idx(rf: i32) -> u32 {
    rf.wrapping_sub(1) as u32
}

/// `MinusId2PlusRFIdx`: the plus-strand frame index of the bases of a minus frame (-1 to -3;
/// any negative frame, in `unsigned int` arithmetic).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:211-215
pub(crate) fn minus_id_to_plus_rf_idx(rf: i32, na_len: u32) -> u32 {
    na_len.wrapping_add(rf as u32).wrapping_add(1) % RF_SIZE
}

/// `PlusNA2PlusPr`: the residue of a plus-strand frame whose codon holds base `na`
/// (`unsigned int` arithmetic, stored in an `int`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:251-254
pub(crate) fn plus_na_to_plus_pr(na: i64, rfidx: u32) -> i32 {
    ((na as u32).wrapping_sub(rfidx) / RF_SIZE) as i32
}

/// `PlusNA2MinusPr`: the residue of a minus-strand frame whose codon holds base `na` of the
/// plus strand (`unsigned int` arithmetic, stored in an `int`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:256-261
pub(crate) fn plus_na_to_minus_pr(na: i64, rfidx: u32, na_len: u32) -> i32 {
    (na_len
        .wrapping_sub(na as u32)
        .wrapping_sub(rfidx)
        .wrapping_add(2)
        / RF_SIZE) as i32
}

/// The `ncbieaa` string of genetic code 1 (`gc.prt`, "Standard"), in TCAG order.
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/gc.inc:40-44
const NCBIEAA_1: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// `CTrans_table::sm_BaseToIdx`: the IUPACna letters (either case, `U` as `T`, `X` as `N`) and
/// the NCBI4na codes to their 4-bit base sets; every other byte to 0 (a gap).
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:90-117
const fn base_to_idx_table() -> [u8; 256] {
    const CHAR_TO_BASE: &[u8; 16] = b"-ACMGRSVTWYHKDBN";
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 16 {
        let ch = CHAR_TO_BASE[i];
        table[ch as usize] = i as u8;
        table[ch.to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    table[b'U' as usize] = 8;
    table[b'u' as usize] = 8;
    table[b'X' as usize] = 15;
    table[b'x' as usize] = 15;
    let mut i = 0;
    while i < 16 {
        table[i] = i as u8;
        i += 1;
    }
    table
}

const BASE_TO_IDX: [u8; 256] = base_to_idx_table();

/// `CTrans_table::m_AminoAcid` of genetic code 1, indexed by the three 4-bit base sets of a
/// codon (`256 * i + 16 * j + k`; NCBI's state is that plus 1).
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:139-232
// (`X` is kept by every later merge, so stopping early changes no amino acid.)
const fn amino_acid_table() -> [u8; 4096] {
    const EXPANSIONS: [usize; 4] = [1, 2, 4, 8];
    const CODON_IDX: [usize; 9] = [0, 2, 1, 0, 3, 0, 0, 0, 0];
    let mut table = [b'X'; 4096];
    let mut i = 0;
    while i < 16 {
        let mut j = 0;
        while j < 16 {
            let mut k = 0;
            while k < 16 {
                let mut aa = 0u8;
                let mut p = 0;
                while p < 4 {
                    let x = EXPANSIONS[p];
                    if x & i != 0 {
                        let mut q = 0;
                        while q < 4 {
                            let y = EXPANSIONS[q];
                            if y & j != 0 {
                                let mut r = 0;
                                while r < 4 {
                                    let z = EXPANSIONS[r];
                                    if z & k != 0 {
                                        let cd =
                                            16 * CODON_IDX[x] + 4 * CODON_IDX[y] + CODON_IDX[z];
                                        let ch = NCBIEAA_1[cd];
                                        if aa == 0 {
                                            aa = ch;
                                        } else if aa != ch {
                                            aa = if matches!(aa, b'B' | b'D' | b'N')
                                                && matches!(ch, b'D' | b'N')
                                            {
                                                b'B'
                                            } else if matches!(aa, b'Z' | b'E' | b'Q')
                                                && matches!(ch, b'E' | b'Q')
                                            {
                                                b'Z'
                                            } else if matches!(aa, b'J' | b'I' | b'L')
                                                && matches!(ch, b'I' | b'L')
                                            {
                                                b'J'
                                            } else {
                                                b'X'
                                            };
                                        }
                                    }
                                    r += 1;
                                }
                            }
                            q += 1;
                        }
                    }
                    p += 1;
                }
                if aa != 0 {
                    table[256 * i + 16 * j + k] = aa;
                }
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }
    table
}

const AMINO_ACIDS: [u8; 4096] = amino_acid_table();

/// The amino acid of the codon `b1 b2 b3` (bytes as in the sequence).
fn codon_residue(b1: u8, b2: u8, b3: u8) -> u8 {
    let (i, j, k) = (
        usize::from(BASE_TO_IDX[usize::from(b1)]),
        usize::from(BASE_TO_IDX[usize::from(b2)]),
        usize::from(BASE_TO_IDX[usize::from(b3)]),
    );
    AMINO_ACIDS[256 * i + 16 * j + k]
}

/// `CSeqTranslator::Translate` (with `fIs5PrimePartial` and genetic code 1) with rpsbproc's `*` to
/// `X`: one residue per codon, then a last partial codon (1 or 2 bases) padded with `N` and kept
/// only when it is not `X`.
// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:4108-4120
// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3636-3730
// `fIs5PrimePartial` makes `check_start` false. The state after three bases is that of the
// three bases alone (`sm_NextState`, Genetic_code_table.cpp:119-134: state 0 reads as `NN`).
// The stop check at the end changes nothing for code 1: `m_OrfStop` is `*` only for codons
// whose every expansion is a stop in `sncbieaa`, and those are stops in `ncbieaa` too, so the
// last residue is `*` already.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:129-133
fn translate(seq: &[u8]) -> Vec<u8> {
    let mut protein = Vec::with_capacity(seq.len().div_ceil(3));
    let mut codons = seq.chunks_exact(3);
    for codon in codons.by_ref() {
        protein.push(codon_residue(codon[0], codon[1], codon[2]));
    }
    match *codons.remainder() {
        [b1] => {
            let aa = codon_residue(b1, b'N', b'N');
            if aa != b'X' {
                protein.push(aa);
            }
        }
        [b1, b2] => {
            let aa = codon_residue(b1, b2, b'N');
            if aa != b'X' {
                protein.push(aa);
            }
        }
        _ => {}
    }
    for residue in &mut protein {
        if *residue == b'*' {
            *residue = b'X';
        }
    }
    protein
}

/// `GetComplementSeq`: the reverse complement in upper case; a byte that is not an IUPACna
/// letter (in either case, `U` read as `T`) is left out.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:3061-3121
// (`toupper` of a byte above 0x7F gives no ASCII letter, so such a byte is left out too.)
pub(crate) fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .filter_map(|&byte| {
            Some(match byte.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' | b'U' => b'A',
                b'R' => b'Y',
                b'Y' => b'R',
                b'M' => b'K',
                b'K' => b'M',
                b'S' => b'W',
                b'W' => b'S',
                b'H' => b'D',
                b'D' => b'H',
                b'B' => b'V',
                b'V' => b'B',
                b'N' => b'N',
                _ => return None,
            })
        })
        .collect()
}

/// `TranslateAll` with genetic code 1: the translations of the three plus-strand frames (the
/// sequence from its first, second and third base) and of the three minus-strand frames (its
/// complement from its first, second and third base). `None` where NCBI erases the first byte
/// of an empty string (a sequence of 1 or 2 bases, or a complement of fewer than 3): libstdc++
/// then moves `SIZE_MAX` bytes (undefined, a crash).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:113-145
// The genetic code is `TDocsum::m_iGenCode`, which the archive never sets (`basedata.cpp:431`,
// which sets `m_iGenCode` to 1; the `gcode` argument of `ProcessBlastResults` is not used), so
// rpsbproc translates with the standard code whatever `-query_gencode` the search used.
pub(crate) fn translate_all(na: &[u8]) -> Option<[Vec<u8>; 6]> {
    let mut frames: [Vec<u8>; 6] = Default::default();
    if na.is_empty() {
        return Some(frames);
    }
    let minus = complement(na);
    if na.len() < 3 || minus.len() < 3 {
        return None;
    }
    for i in 0..3 {
        frames[i] = translate(&na[i..]);
        frames[i + 3] = translate(&minus[i..]);
    }
    Some(frames)
}

#[cfg(test)]
#[path = "frames_tests.rs"]
mod tests;
