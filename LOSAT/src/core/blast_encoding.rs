//! BLAST Sequence Encoding
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c
//!            ncbi-blast/c++/src/algo/blast/core/blast_util.c
//!
//! This module provides the nucleotide encodings of NCBI BLAST that the engines use:
//! - IUPAC letters to BLASTNA codes (`encode_iupac_to_blastna`,
//!   `encode_iupac_to_blastna_into`) and to ncbi4na codes (`iupacna_to_ncbi4na`)
//! - NCBI's `CRandom` (`NcbiRandom`) and the resolution of an ambiguous ncbi4na letter
//!   with it (`resolve_ncbi4na_base`)
//! - the 2-bit packed subject (ncbi2na) built in one pass (`encode_subject_ncbi2na_packed`)
//!
//! Test-only (`#[cfg(test)]`): `resolve_ncbi4na_to_ncbi2na`, the two-pass form that is
//! the oracle of `encode_subject_ncbi2na_packed`.
//!
//! # Encoding Scheme (NCBI BLAST compatible)
//! - A = 0b00 (0)
//! - C = 0b01 (1)
//! - G = 0b10 (2)
//! - T/U = 0b11 (3)
//!
//! # Packing Order
//! 4 nucleotides are packed into each byte, most significant bits first:
//! - Base 0: bits 6-7 (shift 6)
//! - Base 1: bits 4-5 (shift 4)
//! - Base 2: bits 2-3 (shift 2)
//! - Base 3: bits 0-1 (shift 0)

/// Compression ratio: 4 nucleotides per byte
pub const COMPRESSION_RATIO: usize = 4;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93
// ```c
// const Uint1 IUPACNA_TO_BLASTNA[128]={
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
// 15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};
// ```
const IUPACNA_TO_BLASTNA: [u8; 128] = [
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 10, 1, 11, 15, 15, 2,
    12, 15, 15, 7, 15, 6, 14, 15, 15, 15, 4, 9, 3, 15, 13, 8, 15, 5, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15,
];

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93
// ```c
// const Uint1 IUPACNA_TO_BLASTNA[128]={
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
// 15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};
// ```
const IUPACNA_TO_BLASTNA_FULL: [u8; 256] = {
    let mut table = [15u8; 256];
    let mut i = 0usize;
    while i < 128 {
        table[i] = IUPACNA_TO_BLASTNA[i];
        i += 1;
    }
    let mut c = b'a';
    while c <= b'z' {
        table[c as usize] = table[(c - 32) as usize];
        c += 1;
    }
    table
};

/// Encode an ASCII sequence to BLASTNA (IUPAC) codes (one base per byte).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93
/// ```c
/// const Uint1 IUPACNA_TO_BLASTNA[128]={
/// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
/// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
/// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
/// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
/// 15, 0,10, 1,11,15,15, 2,12,15,15, 7,15, 6,14,15,
/// 15,15, 4, 9, 3,15,13, 8,15, 5,15,15,15,15,15,15,
/// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
/// 15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15};
/// ```
pub fn encode_iupac_to_blastna(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &base in seq {
        out.push(IUPACNA_TO_BLASTNA_FULL[base as usize]);
    }
    out
}

/// `encode_iupac_to_blastna` into a caller-owned slice (`out.len() == seq.len()`), so the
/// search writes a context straight into its query buffer.
pub fn encode_iupac_to_blastna_into(seq: &[u8], out: &mut [u8]) {
    assert_eq!(seq.len(), out.len());
    for (slot, &base) in out.iter_mut().zip(seq) {
        *slot = IUPACNA_TO_BLASTNA_FULL[base as usize];
    }
}

/// The ncbi4na value of an IUPAC nucleotide letter, either case (bit 0 A, 1 C, 2 G, 3 T;
/// a gap is 0 and `N` 15), or `None` for another byte.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_objmgr_tools.cpp:424-425
/// ```c
/// static unsigned char ctable[16] = {0xFF, 0x00, 0x01, 0xFF, 0x02, 0xFF, 0xFF, 0xFF,
/// 		                           0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
/// ```
pub fn iupacna_to_ncbi4na(letter: u8) -> Option<u8> {
    Some(match letter.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
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
        b'-' => 0,
        _ => return None,
    })
}

// NCBI reference: c++/src/util/random_gen.cpp:98,227-230,287-308 and
// c++/include/util/random_gen.hpp:224-241
// ```c
// static const size_t kStateOffset = 12;
// m_State[0] = m_Seed = seed;
// for (int i = 1; i < kStateSize; ++i)
//     m_State[i] = 1103515245 * m_State[i-1] + 12345;
// m_RJ = kStateOffset; m_RK = kStateSize - 1;
// for (int i = 0; i < 10 * kStateSize; ++i) GetRand();
// r = m_State[m_RK] + m_State[m_RJ--];
// m_State[m_RK--] = r;
// return r >> 1;
// ```
/// NCBI's `CRandom` (the lagged Fibonacci generator), which resolves the ambiguous
/// letters of a subject (`resolve_ncbi4na_to_ncbi2na`).
pub(crate) struct NcbiRandom {
    state: [u32; 33],
    j: usize,
    k: usize,
}

impl NcbiRandom {
    pub(crate) fn new(seed: u32) -> Self {
        let mut state = [0; 33];
        state[0] = seed;
        for i in 1..state.len() {
            state[i] = state[i - 1]
                .wrapping_mul(1_103_515_245)
                .wrapping_add(12_345);
        }
        let mut random = Self {
            state,
            j: 12,
            k: 32,
        };
        for _ in 0..330 {
            random.get_rand();
        }
        random
    }

    pub(crate) fn get_rand(&mut self) -> u32 {
        let value = self.state[self.k].wrapping_add(self.state[self.j]);
        self.state[self.k] = value;
        self.j = if self.j == 0 { 32 } else { self.j - 1 };
        self.k = if self.k == 0 { 32 } else { self.k - 1 };
        value >> 1
    }
}

/// The ncbi2na codes (0 A, 1 C, 2 G, 3 T) of an ncbi4na sequence, as NCBI makes the
/// compressed subject that the preliminary search reads: each ambiguous letter becomes a
/// compatible base drawn from `CRandom` seeded with the length of the sequence.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474
/// ```c
/// void s_Ncbi4naToNcbi2na(const string & ncbi4na, int base_length,
///                         unsigned char * ncbi2na)
/// {
///     int inp_bytes   = base_length;
///     CRandom random(base_length);
///     ...
///         if (c  != 0xFF) {
///             // No ambiguities, so we can do this the easy way.
///         	ncbi2na[i] = c;
///
///         } else {
///             if (b == 0 || b == 0x0F) {
///             	//gap or N
///                 ncbi2na[i] = random.GetRand() & 0x3;
///             }
///             else {
///     ...
///             	int pick = random.GetRand() % bitcount;
/// ```
#[cfg(test)]
pub fn resolve_ncbi4na_to_ncbi2na(ncbi4na: &[u8]) -> Vec<u8> {
    let mut random = NcbiRandom::new(ncbi4na.len() as u32);
    ncbi4na
        .iter()
        .map(|&mask| resolve_ncbi4na_base(mask, &mut random))
        .collect()
}

/// The ncbi2na code of one ncbi4na code, an ambiguous one drawn from `random` (the next base
/// of `resolve_ncbi4na_to_ncbi2na`).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474 (above)
#[inline]
pub(crate) fn resolve_ncbi4na_base(mask: u8, random: &mut NcbiRandom) -> u8 {
    match mask {
        1 => 0,
        2 => 1,
        4 => 2,
        8 => 3,
        0 | 15 => (random.get_rand() & 0x3) as u8,
        _ => {
            let mut pick = random.get_rand() % mask.count_ones();
            (0..4u8)
                .filter(|bit| mask & (1 << bit) != 0)
                .find(|_| {
                    let found = pick == 0;
                    pick = pick.wrapping_sub(1);
                    found
                })
                .unwrap_or(0)
        }
    }
}

/// The value of `UNAMBIGUOUS_NCBI2NA` for a byte that is not an unambiguous letter (the
/// one bit that no ncbi2na code has).
const AMBIGUOUS_LETTER: u8 = 4;

/// The ncbi2na code of each byte that is `A`, `C`, `G` or `T` in either case, and
/// `AMBIGUOUS_LETTER` for every other byte (`encode_subject_ncbi2na_packed`).
const UNAMBIGUOUS_NCBI2NA: [u8; 256] = {
    let mut table = [AMBIGUOUS_LETTER; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// A subject in packed ncbi2na (4 bases per byte, the remainder count in the last byte),
/// with its ambiguous letters resolved as NCBI resolves them (`resolve_ncbi4na_to_ncbi2na`).
/// The letters have been checked (IUPAC); another byte is read as `N`.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
/// ```c
///     virtual SBlastSequence GetCompressedPlusStrand() {
///         SBlastSequence retval(size());
///         string ncbi4na = kEmptyStr;
///         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
///         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
///         return retval;
/// ```
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:1167-1194
/// ```c
///     // Populate the compressed sequence up to the last byte
///     for (ci = 0, i = 0; ci < retval.length-1; ci++, i+= COMPRESSION_RATIO) {
///         Uint1 a, b, c, d;
///         a = ((*source_ptr & NCBI2NA_MASK)<<6); ++source_ptr;
///     ...
///     // Set the number of bases in the last 2 bits of the last byte in the
///     // compressed sequence
///     retval.data.get()[ci] |= source.length%COMPRESSION_RATIO;
/// ```
///
/// The letters are read once: `A`, `C`, `G` and `T` through a table, the others as
/// `resolve_ncbi4na_to_ncbi2na` resolves their ncbi4na values, in order, from a `CRandom`
/// seeded with the length and made at the first of them (the generator is not drawn from
/// for unambiguous letters, so a subject without one needs none).
pub fn encode_subject_ncbi2na_packed(seq: &[u8]) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }
    let length = seq.len();
    let mut random: Option<NcbiRandom> = None;
    let mut code = |letter: u8| -> u8 {
        let code = UNAMBIGUOUS_NCBI2NA[letter as usize];
        if code != AMBIGUOUS_LETTER {
            return code;
        }
        let mask = iupacna_to_ncbi4na(letter).unwrap_or(15);
        let random = random.get_or_insert_with(|| NcbiRandom::new(length as u32));
        resolve_ncbi4na_base(mask, random)
    };
    let mut packed = Vec::with_capacity(length / COMPRESSION_RATIO + 1);
    let mut groups = seq.chunks_exact(COMPRESSION_RATIO);
    for group in &mut groups {
        let first = UNAMBIGUOUS_NCBI2NA[group[0] as usize];
        let second = UNAMBIGUOUS_NCBI2NA[group[1] as usize];
        let third = UNAMBIGUOUS_NCBI2NA[group[2] as usize];
        let fourth = UNAMBIGUOUS_NCBI2NA[group[3] as usize];
        if (first | second | third | fourth) & AMBIGUOUS_LETTER == 0 {
            packed.push((first << 6) | (second << 4) | (third << 2) | fourth);
            continue;
        }
        // At least one ambiguous letter: the four letters in order.
        let first = code(group[0]);
        let second = code(group[1]);
        let third = code(group[2]);
        let fourth = code(group[3]);
        packed.push((first << 6) | (second << 4) | (third << 2) | fourth);
    }
    let mut last = (length % COMPRESSION_RATIO) as u8;
    for (index, &letter) in groups.remainder().iter().enumerate() {
        last |= code(letter) << (6 - 2 * index);
    }
    packed.push(last);
    packed
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `encode_subject_ncbi2na_packed` as it was written before its one-pass form: every
    /// letter to ncbi4na, all of them resolved, then packed.
    fn encode_subject_ncbi2na_packed_two_pass(seq: &[u8]) -> Vec<u8> {
        if seq.is_empty() {
            return Vec::new();
        }
        let ncbi4na: Vec<u8> = seq
            .iter()
            .map(|&letter| iupacna_to_ncbi4na(letter).unwrap_or(15))
            .collect();
        let codes = resolve_ncbi4na_to_ncbi2na(&ncbi4na);
        let mut packed = vec![0u8; codes.len() / COMPRESSION_RATIO + 1];
        for (index, &code) in codes.iter().enumerate() {
            let shift = 6 - 2 * (index % COMPRESSION_RATIO);
            packed[index / COMPRESSION_RATIO] |= code << shift;
        }
        *packed.last_mut().expect("at least one byte") |= (codes.len() % COMPRESSION_RATIO) as u8;
        packed
    }

    #[test]
    fn subject_packing_equals_the_two_pass_form() {
        // Every byte value alone and between unambiguous letters.
        for byte in 0..=255u8 {
            for seq in [vec![byte], vec![b'A', byte, b'c', b'G', b't', byte, byte]] {
                assert_eq!(
                    encode_subject_ncbi2na_packed(&seq),
                    encode_subject_ncbi2na_packed_two_pass(&seq),
                    "{seq:?}"
                );
            }
        }
        // Mixed letters of every length up to 64 and a few long ones (a linear congruential
        // generator picks them; one letter in eight is ambiguous or not a letter).
        let alphabet = b"ACGTacgtACGTacgtACGTacgtACGTacgtMRSVWYHKDBNmrsvwyhkdbn-UuX*";
        let mut state = 0x2545_f491_u32;
        let mut next = || {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            state >> 8
        };
        let lengths = (0..=64).chain([1_000, 4_097, 100_003]);
        for length in lengths {
            for _ in 0..4 {
                let seq: Vec<u8> = (0..length)
                    .map(|_| {
                        let pick = next() as usize;
                        if pick.is_multiple_of(8) {
                            alphabet[32 + pick / 8 % (alphabet.len() - 32)]
                        } else {
                            alphabet[pick / 8 % 32]
                        }
                    })
                    .collect();
                assert_eq!(
                    encode_subject_ncbi2na_packed(&seq),
                    encode_subject_ncbi2na_packed_two_pass(&seq),
                    "length {length}"
                );
            }
        }
    }
}
