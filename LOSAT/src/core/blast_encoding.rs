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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
const IUPACNA_TO_BLASTNA: [u8; 128] = [
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15,  0, 10,  1, 11, 15, 15,  2, 12, 15, 15,  7, 15,  6, 14, 15,
    15, 15,  4,  9,  3, 15, 13,  8, 15,  5, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
];

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
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
#[inline]
fn encode_iupac_base_to_blastna(base: u8) -> u8 {
    let upper = base.to_ascii_uppercase();
    let idx = upper as usize;
    if idx < IUPACNA_TO_BLASTNA.len() {
        IUPACNA_TO_BLASTNA[idx]
    } else {
        15
    }
}

/// Encode an ASCII sequence to BLASTNA (IUPAC) codes (one base per byte).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:85-93 (IUPACNA_TO_BLASTNA)
pub fn encode_iupac_to_blastna(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&base| encode_iupac_base_to_blastna(base)).collect()
}

/// Encode an ASCII sequence to ncbi2na 2-bit codes (one base per byte).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:476-489 (BlastCompressBlastnaSequence old_seq[i] & 3)
pub fn encode_iupac_to_ncbi2na(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&base| encode_iupac_base_to_blastna(base) & 0x03)
        .collect()
}

/// Encode an ASCII sequence to packed ncbi2na (4 bases per byte, remainder count in last byte).
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:1154-1187 (CompressNcbi2na)
pub fn encode_iupac_to_ncbi2na_packed(seq: &[u8]) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }

    let packed_len = seq.len() / COMPRESSION_RATIO + 1;
    let mut packed = vec![0u8; packed_len];

    let mut i = 0usize;
    let mut byte_idx = 0usize;

    while i + COMPRESSION_RATIO <= seq.len() {
        let a = (encode_iupac_base_to_blastna(seq[i]) & 0x03) << 6;
        let b = (encode_iupac_base_to_blastna(seq[i + 1]) & 0x03) << 4;
        let c = (encode_iupac_base_to_blastna(seq[i + 2]) & 0x03) << 2;
        let d = (encode_iupac_base_to_blastna(seq[i + 3]) & 0x03) << 0;
        packed[byte_idx] = a | b | c | d;
        byte_idx += 1;
        i += COMPRESSION_RATIO;
    }

    let mut last_byte = 0u8;
    while i < seq.len() {
        let bit_shift = match i % COMPRESSION_RATIO {
            0 => 6,
            1 => 4,
            2 => 2,
            _ => 0,
        };
        let code = encode_iupac_base_to_blastna(seq[i]) & 0x03;
        last_byte |= code << bit_shift;
        i += 1;
    }

    last_byte |= (seq.len() % COMPRESSION_RATIO) as u8;
    packed[byte_idx] = last_byte;

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
