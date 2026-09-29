//! BLAST Sequence Encoding
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c
//!            ncbi-blast/c++/src/algo/blast/core/blast_util.c
//!
//! This module implements NCBI BLAST's 2-bit packed nucleotide format (ncbi2na)
//! for efficient sequence storage and k-mer extraction.
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

/// Bit mask for extracting a single 2-bit base
const BASE_MASK: u8 = 0x03;

/// Lookup table for encoding ASCII nucleotides to 2-bit codes
/// Returns 0xFF for invalid/ambiguous bases
const ENCODE_TABLE: [u8; 256] = {
    let mut table = [0xFFu8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table[b'U' as usize] = 3;
    table[b'u' as usize] = 3;
    table
};

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

/// Lookup table for decoding 2-bit codes to ASCII nucleotides
const DECODE_TABLE: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// A 2-bit packed nucleotide sequence
///
/// Stores 4 nucleotides per byte using NCBI BLAST's ncbi2na encoding.
/// This reduces memory usage by 4x compared to 1-byte-per-base representation
/// and improves cache efficiency for sequence scanning operations.
#[derive(Debug, Clone)]
pub struct PackedSequence {
    /// Packed sequence data (4 nucleotides per byte)
    data: Vec<u8>,
    /// Original sequence length in nucleotides
    len: usize,
    /// Positions of ambiguous bases (N, etc.) - None means no ambiguous bases
    ambiguous_positions: Option<Vec<usize>>,
}

impl PackedSequence {
    /// Create a new packed sequence from an ASCII nucleotide sequence
    ///
    /// # Arguments
    /// * `seq` - ASCII nucleotide sequence (A, C, G, T/U, case-insensitive)
    ///
    /// # Returns
    /// * `Some(PackedSequence)` if the sequence was successfully packed
    /// * `None` if the sequence is empty
    ///
    /// # Note
    /// Ambiguous bases (N, etc.) are stored as 0 and their positions are tracked
    /// separately. K-mers containing ambiguous bases will return None.
    pub fn new(seq: &[u8]) -> Option<Self> {
        if seq.is_empty() {
            return None;
        }

        let len = seq.len();
        let packed_len = (len + COMPRESSION_RATIO - 1) / COMPRESSION_RATIO;
        let mut data = vec![0u8; packed_len];
        let mut ambiguous_positions: Vec<usize> = Vec::new();

        for (i, &base) in seq.iter().enumerate() {
            let code = ENCODE_TABLE[base as usize];
            if code == 0xFF {
                // Ambiguous base - store as 0 and track position
                ambiguous_positions.push(i);
            } else {
                let byte_idx = i / COMPRESSION_RATIO;
                let bit_offset = 6 - 2 * (i % COMPRESSION_RATIO);
                data[byte_idx] |= code << bit_offset;
            }
        }

        Some(Self {
            data,
            len,
            ambiguous_positions: if ambiguous_positions.is_empty() {
                None
            } else {
                Some(ambiguous_positions)
            },
        })
    }

    /// Get the length of the sequence in nucleotides
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Check if the sequence is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the packed data as a slice
    #[inline]
    pub fn data(&self) -> &[u8] {
        &self.data
    }

    /// Extract a single base at the given position
    ///
    /// # Arguments
    /// * `pos` - Position in the sequence (0-based)
    ///
    /// # Returns
    /// * 2-bit encoded base (0=A, 1=C, 2=G, 3=T)
    ///
    /// # Panics
    /// Panics if `pos >= len`
    #[inline]
    pub fn get_base(&self, pos: usize) -> u8 {
        debug_assert!(pos < self.len, "Position out of bounds");
        let byte_idx = pos / COMPRESSION_RATIO;
        let bit_offset = 6 - 2 * (pos % COMPRESSION_RATIO);
        (self.data[byte_idx] >> bit_offset) & BASE_MASK
    }

    /// Extract a single base at the given position, returning ASCII
    ///
    /// # Arguments
    /// * `pos` - Position in the sequence (0-based)
    ///
    /// # Returns
    /// * ASCII nucleotide character (A, C, G, T)
    #[inline]
    pub fn get_base_ascii(&self, pos: usize) -> u8 {
        DECODE_TABLE[self.get_base(pos) as usize]
    }

    /// Check if a position contains an ambiguous base
    #[inline]
    pub fn is_ambiguous(&self, pos: usize) -> bool {
        if let Some(ref positions) = self.ambiguous_positions {
            positions.binary_search(&pos).is_ok()
        } else {
            false
        }
    }

    /// Check if a range contains any ambiguous bases
    ///
    /// # Arguments
    /// * `start` - Start position (inclusive)
    /// * `end` - End position (exclusive)
    #[inline]
    pub fn has_ambiguous_in_range(&self, start: usize, end: usize) -> bool {
        if let Some(ref positions) = self.ambiguous_positions {
            // Binary search for the first position >= start
            match positions.binary_search(&start) {
                Ok(_) => true, // Exact match found
                Err(idx) => {
                    // Check if any position in the range exists
                    idx < positions.len() && positions[idx] < end
                }
            }
        } else {
            false
        }
    }

    /// Extract a k-mer at the given position
    ///
    /// # Arguments
    /// * `pos` - Starting position of the k-mer (0-based)
    /// * `k` - K-mer length
    ///
    /// # Returns
    /// * `Some(u64)` - 2-bit encoded k-mer if valid
    /// * `None` - If the k-mer extends beyond the sequence or contains ambiguous bases
    #[inline]
    pub fn extract_kmer(&self, pos: usize, k: usize) -> Option<u64> {
        if pos + k > self.len {
            return None;
        }

        // Check for ambiguous bases in the k-mer range
        if self.has_ambiguous_in_range(pos, pos + k) {
            return None;
        }

        let mut kmer: u64 = 0;
        for i in 0..k {
            let base = self.get_base(pos + i);
            kmer = (kmer << 2) | (base as u64);
        }
        Some(kmer)
    }

    /// Extract a k-mer using sliding window optimization
    ///
    /// Given the previous k-mer and the new base, compute the next k-mer
    /// in O(1) time instead of O(k).
    ///
    /// # Arguments
    /// * `prev_kmer` - Previous k-mer value
    /// * `new_base` - New base to add (2-bit encoded)
    /// * `k` - K-mer length
    ///
    /// # Returns
    /// * New k-mer value
    #[inline]
    pub fn sliding_kmer(prev_kmer: u64, new_base: u8, k: usize) -> u64 {
        let mask = (1u64 << (2 * k)) - 1;
        ((prev_kmer << 2) | (new_base as u64)) & mask
    }

    /// Create an iterator over all valid k-mers in the sequence
    ///
    /// # Arguments
    /// * `k` - K-mer length
    ///
    /// # Returns
    /// Iterator yielding (position, kmer_code) pairs for valid k-mers
    pub fn iter_kmers(&self, k: usize) -> KmerIterator<'_> {
        KmerIterator::new(self, k)
    }

    /// Unpack the sequence back to ASCII format
    pub fn unpack(&self) -> Vec<u8> {
        let mut result = Vec::with_capacity(self.len);
        for i in 0..self.len {
            if self.is_ambiguous(i) {
                result.push(b'N');
            } else {
                result.push(self.get_base_ascii(i));
            }
        }
        result
    }
}

/// Iterator over k-mers in a packed sequence
pub struct KmerIterator<'a> {
    seq: &'a PackedSequence,
    k: usize,
    pos: usize,
    current_kmer: Option<u64>,
}

impl<'a> KmerIterator<'a> {
    fn new(seq: &'a PackedSequence, k: usize) -> Self {
        let mut iter = Self {
            seq,
            k,
            pos: 0,
            current_kmer: None,
        };
        // Initialize the first k-mer
        iter.initialize_kmer();
        iter
    }

    fn initialize_kmer(&mut self) {
        if self.seq.len() < self.k {
            return;
        }

        // Try to find the first valid k-mer
        while self.pos + self.k <= self.seq.len() {
            if let Some(kmer) = self.seq.extract_kmer(self.pos, self.k) {
                self.current_kmer = Some(kmer);
                return;
            }
            self.pos += 1;
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        let kmer = self.current_kmer?;
        let current_pos = self.pos;

        // Advance to next position
        self.pos += 1;

        if self.pos + self.k <= self.seq.len() {
            // Check if the new position has an ambiguous base at the end
            let new_pos = self.pos + self.k - 1;
            if self.seq.is_ambiguous(new_pos) {
                // Need to skip ahead and reinitialize
                self.pos += 1;
                self.current_kmer = None;
                while self.pos + self.k <= self.seq.len() {
                    if let Some(new_kmer) = self.seq.extract_kmer(self.pos, self.k) {
                        self.current_kmer = Some(new_kmer);
                        break;
                    }
                    self.pos += 1;
                }
            } else {
                // Use sliding window optimization
                let new_base = self.seq.get_base(new_pos);
                self.current_kmer = Some(PackedSequence::sliding_kmer(kmer, new_base, self.k));
            }
        } else {
            self.current_kmer = None;
        }

        Some((current_pos, kmer))
    }
}

/// Encode a single ASCII nucleotide to 2-bit code
///
/// # Returns
/// * `Some(code)` for valid bases (A, C, G, T/U)
/// * `None` for ambiguous or invalid bases
#[inline]
pub fn encode_base(base: u8) -> Option<u8> {
    let code = ENCODE_TABLE[base as usize];
    if code == 0xFF {
        None
    } else {
        Some(code)
    }
}

/// Encode a single ASCII nucleotide to BLASTNA (IUPAC) code.
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
#[inline]
fn encode_iupac_base_to_blastna(base: u8) -> u8 {
    IUPACNA_TO_BLASTNA_FULL[base as usize]
}

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

/// Encode an ASCII sequence to ncbi2na 2-bit codes (one base per byte).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:476-489
/// ```c
/// new_seq[i] = (Uint1)(old_seq[i] & 3);
/// ```
pub fn encode_iupac_to_ncbi2na(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &base in seq {
        out.push(IUPACNA_TO_BLASTNA_FULL[base as usize] & 0x03);
    }
    out
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
pub fn resolve_ncbi4na_to_ncbi2na(ncbi4na: &[u8]) -> Vec<u8> {
    let mut random = NcbiRandom::new(ncbi4na.len() as u32);
    ncbi4na
        .iter()
        .map(|&mask| match mask {
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
        })
        .collect()
}

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
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:1154-1187
/// ```c
/// for (i=0; i<length; i += 4) {
///     Uint1 encoded = (Uint1)(seq[i] & 3) << 6;
///     ...
///     packed[j++] = encoded;
/// }
/// packed[j] |= (Uint1)(length % 4);
/// ```
pub fn encode_subject_ncbi2na_packed(seq: &[u8]) -> Vec<u8> {
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

/// Decode a 2-bit code to ASCII nucleotide
#[inline]
pub fn decode_base(code: u8) -> u8 {
    DECODE_TABLE[(code & BASE_MASK) as usize]
}

/// Encode a k-mer from an ASCII sequence (for compatibility with existing code)
///
/// # Arguments
/// * `seq` - ASCII nucleotide sequence
/// * `start` - Starting position
/// * `k` - K-mer length
///
/// # Returns
/// * `Some(u64)` - Encoded k-mer if all bases are valid
/// * `None` - If any base is ambiguous or invalid
#[inline]
pub fn encode_kmer_from_ascii(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    if start + k > seq.len() {
        return None;
    }

    let mut kmer: u64 = 0;
    for i in 0..k {
        let base = unsafe { *seq.get_unchecked(start + i) };
        let code = ENCODE_TABLE[base as usize];
        if code == 0xFF {
            return None;
        }
        kmer = (kmer << 2) | (code as u64);
    }
    Some(kmer)
}

/// Decode a k-mer to ASCII sequence
pub fn decode_kmer(kmer: u64, k: usize) -> Vec<u8> {
    let mut result = vec![0u8; k];
    let mut code = kmer;
    for i in (0..k).rev() {
        result[i] = DECODE_TABLE[(code & 3) as usize];
        code >>= 2;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), Some(0));
        assert_eq!(encode_base(b'a'), Some(0));
        assert_eq!(encode_base(b'C'), Some(1));
        assert_eq!(encode_base(b'c'), Some(1));
        assert_eq!(encode_base(b'G'), Some(2));
        assert_eq!(encode_base(b'g'), Some(2));
        assert_eq!(encode_base(b'T'), Some(3));
        assert_eq!(encode_base(b't'), Some(3));
        assert_eq!(encode_base(b'U'), Some(3));
        assert_eq!(encode_base(b'u'), Some(3));
        assert_eq!(encode_base(b'N'), None);
        assert_eq!(encode_base(b'n'), None);
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0), b'A');
        assert_eq!(decode_base(1), b'C');
        assert_eq!(decode_base(2), b'G');
        assert_eq!(decode_base(3), b'T');
    }

    #[test]
    fn test_packed_sequence_new() {
        let seq = b"ACGT";
        let packed = PackedSequence::new(seq).unwrap();
        assert_eq!(packed.len(), 4);
        assert_eq!(packed.data().len(), 1);
    }

    #[test]
    fn test_packed_sequence_get_base() {
        let seq = b"ACGTACGT";
        let packed = PackedSequence::new(seq).unwrap();
        assert_eq!(packed.get_base(0), 0); // A
        assert_eq!(packed.get_base(1), 1); // C
        assert_eq!(packed.get_base(2), 2); // G
        assert_eq!(packed.get_base(3), 3); // T
        assert_eq!(packed.get_base(4), 0); // A
        assert_eq!(packed.get_base(5), 1); // C
        assert_eq!(packed.get_base(6), 2); // G
        assert_eq!(packed.get_base(7), 3); // T
    }

    #[test]
    fn test_packed_sequence_unpack() {
        let seq = b"ACGTACGT";
        let packed = PackedSequence::new(seq).unwrap();
        let unpacked = packed.unpack();
        assert_eq!(&unpacked, seq);
    }

    #[test]
    fn test_packed_sequence_with_ambiguous() {
        let seq = b"ACNGT";
        let packed = PackedSequence::new(seq).unwrap();
        assert!(packed.is_ambiguous(2));
        assert!(!packed.is_ambiguous(0));
        assert!(!packed.is_ambiguous(1));
        assert!(!packed.is_ambiguous(3));
        assert!(!packed.is_ambiguous(4));
    }

    #[test]
    fn test_extract_kmer() {
        let seq = b"ACGTACGT";
        let packed = PackedSequence::new(seq).unwrap();

        // ACGT = 0b00011011 = 27
        let kmer = packed.extract_kmer(0, 4).unwrap();
        assert_eq!(kmer, 0b00011011);

        // CGTA = 0b01101100 = 108
        let kmer = packed.extract_kmer(1, 4).unwrap();
        assert_eq!(kmer, 0b01101100);
    }

    #[test]
    fn test_extract_kmer_with_ambiguous() {
        let seq = b"ACNGT";
        let packed = PackedSequence::new(seq).unwrap();

        // K-mer containing N should return None
        assert!(packed.extract_kmer(0, 4).is_none());
        assert!(packed.extract_kmer(1, 4).is_none());

        // K-mer not containing N should work
        assert!(packed.extract_kmer(3, 2).is_some());
    }

    #[test]
    fn test_sliding_kmer() {
        // ACGT = 27, adding A should give CGTA = 108
        let prev = 0b00011011u64; // ACGT
        let next = PackedSequence::sliding_kmer(prev, 0, 4); // Add A
        assert_eq!(next, 0b01101100); // CGTA
    }

    #[test]
    fn test_kmer_iterator() {
        let seq = b"ACGTACGT";
        let packed = PackedSequence::new(seq).unwrap();

        let kmers: Vec<(usize, u64)> = packed.iter_kmers(4).collect();
        assert_eq!(kmers.len(), 5); // 8 - 4 + 1 = 5 k-mers

        // Verify first and last k-mers
        assert_eq!(kmers[0], (0, 0b00011011)); // ACGT
        assert_eq!(kmers[4], (4, 0b00011011)); // ACGT
    }

    #[test]
    fn test_kmer_iterator_with_ambiguous() {
        let seq = b"ACNGTACGT";
        let packed = PackedSequence::new(seq).unwrap();

        let kmers: Vec<(usize, u64)> = packed.iter_kmers(4).collect();

        // Should skip k-mers containing N
        // Valid k-mers start at positions 3, 4, 5
        assert!(kmers.iter().all(|(pos, _)| *pos >= 3));
    }

    #[test]
    fn test_encode_kmer_from_ascii() {
        let seq = b"ACGTACGT";

        let kmer = encode_kmer_from_ascii(seq, 0, 4).unwrap();
        assert_eq!(kmer, 0b00011011); // ACGT

        let kmer = encode_kmer_from_ascii(seq, 1, 4).unwrap();
        assert_eq!(kmer, 0b01101100); // CGTA
    }

    #[test]
    fn test_decode_kmer() {
        let kmer = 0b00011011u64; // ACGT
        let decoded = decode_kmer(kmer, 4);
        assert_eq!(&decoded, b"ACGT");
    }

    #[test]
    fn test_packing_matches_ncbi_blast() {
        // Test that our packing matches NCBI BLAST's format
        // NCBI BLAST packs: base0 at bits 6-7, base1 at bits 4-5, etc.
        let seq = b"ACGT";
        let packed = PackedSequence::new(seq).unwrap();

        // A=0, C=1, G=2, T=3
        // Expected byte: (0 << 6) | (1 << 4) | (2 << 2) | 3 = 0b00011011 = 27
        assert_eq!(packed.data()[0], 0b00011011);
    }

    #[test]
    fn test_partial_byte() {
        // Test sequence that doesn't fill the last byte completely
        let seq = b"ACGTA"; // 5 bases = 1 full byte + 1 partial byte
        let packed = PackedSequence::new(seq).unwrap();

        assert_eq!(packed.len(), 5);
        assert_eq!(packed.data().len(), 2);

        // Verify all bases can be extracted correctly
        assert_eq!(packed.get_base(0), 0); // A
        assert_eq!(packed.get_base(1), 1); // C
        assert_eq!(packed.get_base(2), 2); // G
        assert_eq!(packed.get_base(3), 3); // T
        assert_eq!(packed.get_base(4), 0); // A
    }
}
