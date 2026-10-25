//! A record's letters held as 2-bit codes of `A`, `C`, `G` and `T`, the runs of other bytes and the
//! intervals of lower-case bytes, from which the letters are made again exactly (any bytes).
//!
//! No NCBI counterpart: memory only. A query record's letters that wait for their batch take about a
//! quarter of a byte per base instead of one (EXPERIMENT LOSAT_X_BNQPACK,
//! `algorithm/blastn/blast_engine/x_query_drop.rs`). Unlike `packed_subject::PackedNucleotide`, the
//! codes of the other bytes are not drawn (no packed bases are read from it): only the letters are
//! made from it, so no value depends on how those bytes are held.

/// One run of equal bytes that are none of `A`, `C`, `G`, `T` in either case, upper-cased (the
/// lower-case intervals give the case back).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct OtherRun {
    offset: u32,
    length: u32,
    byte: u8,
}

/// The letters of one record (see the module doc).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedLetters {
    /// 4 bases per byte from the high bits (`A` 0, `C` 1, `G` 2, `T` 3; 0 under an other run),
    /// `length.div_ceil(4)` bytes.
    packed: Box<[u8]>,
    length: u32,
    others: Box<[OtherRun]>,
    /// `[start, end)` of every maximal run of lower-case bytes (`u8::is_ascii_lowercase`).
    lowercase: Box<[(u32, u32)]>,
}

/// Bytes `0x01` repeated in each byte of a `u64`.
const ONES: u64 = 0x0101_0101_0101_0101;

/// Packed byte -> its 4 letters.
const LETTERS: [[u8; 4]; 256] = {
    const ACGT: [u8; 4] = *b"ACGT";
    let mut letters = [[0u8; 4]; 256];
    let mut byte = 0;
    while byte < 256 {
        letters[byte] = [
            ACGT[byte >> 6],
            ACGT[(byte >> 4) & 3],
            ACGT[(byte >> 2) & 3],
            ACGT[byte & 3],
        ];
        byte += 1;
    }
    letters
};

/// The state of a packing pass: what the bytes outside the fast blocks add.
struct Packing {
    packed: Vec<u8>,
    others: Vec<OtherRun>,
    lowercase: Vec<(u32, u32)>,
    lower_start: Option<usize>,
}

impl Packing {
    /// Byte `index` (one byte at a time: the blocks that are not all `ACGT` or all `acgt` in the
    /// case of the bytes before them, and the tail).
    fn byte(&mut self, index: usize, byte: u8) {
        match (byte.is_ascii_lowercase(), self.lower_start) {
            (true, None) => self.lower_start = Some(index),
            (false, Some(start)) => {
                self.lowercase.push((start as u32, index as u32));
                self.lower_start = None;
            }
            _ => {}
        }
        let upper = byte.to_ascii_uppercase();
        let code = match upper {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                match self.others.last_mut() {
                    Some(run)
                        if run.byte == upper && (run.offset + run.length) as usize == index =>
                    {
                        run.length += 1
                    }
                    _ => self.others.push(OtherRun {
                        offset: index as u32,
                        length: 1,
                        byte: upper,
                    }),
                }
                0
            }
        };
        self.packed[index / 4] |= code << (6 - 2 * (index % 4));
    }

    /// The case bits (`0x20` per byte) that a fast block must have: those of the byte before it.
    fn case_bits(&self) -> u64 {
        if self.lower_start.is_some() {
            0x20 * ONES
        } else {
            0
        }
    }
}

/// The codes of 8 bytes that are all `A`, `C`, `G`, `T` with the case bits `case` (`0` or `0x20`
/// per byte) as 2 packed bytes, or `None`.
#[inline(always)]
fn pack8(block: u64, case: u64) -> Option<[u8; 2]> {
    let upper = block & !(0x20 * ONES);
    // A 0x41 -> 0, C 0x43 -> 1, G 0x47 -> 2, T 0x54 -> 3: bits 1-2 of the byte, bit 1 flipped by bit 2
    let codes = ((upper >> 1) & (3 * ONES)) ^ ((upper >> 2) & ONES);
    // the letter of each code: 0x41 + 2 b0 + 6 b1 + 11 b0 b1 (at most 0x54, no carry between bytes)
    let low = codes & ONES;
    let high = (codes >> 1) & ONES;
    let both = low & high;
    let letters =
        0x41 * ONES + (low << 1) + (high << 2) + (high << 1) + (both << 3) + (both << 1) + both;
    if letters != upper || block & (0x20 * ONES) != case {
        return None;
    }
    let pairs = ((codes & 0x00FF_00FF_00FF_00FF) << 2) | ((codes >> 8) & 0x00FF_00FF_00FF_00FF);
    let quads = ((pairs & 0x0000_FFFF_0000_FFFF) << 4) | ((pairs >> 16) & 0x0000_FFFF_0000_FFFF);
    Some([quads as u8, (quads >> 32) as u8])
}

impl PackedLetters {
    /// The store of `letters`, or `None` where it holds no fewer bytes than the letters (many other
    /// runs or case changes) or a record longer than a `u32` offset.
    pub fn pack(letters: &[u8]) -> Option<Self> {
        let length = u32::try_from(letters.len()).ok()?;
        let mut packing = Packing {
            packed: vec![0u8; letters.len().div_ceil(4)],
            others: Vec::new(),
            lowercase: Vec::new(),
            lower_start: None,
        };
        let done = pack_blocks(letters, &mut packing);
        for (index, &byte) in letters.iter().enumerate().skip(done) {
            packing.byte(index, byte);
        }
        if let Some(start) = packing.lower_start {
            packing.lowercase.push((start as u32, length));
        }
        let held = packing.packed.len()
            + std::mem::size_of_val(packing.others.as_slice())
            + std::mem::size_of_val(packing.lowercase.as_slice());
        if held >= letters.len() && !letters.is_empty() {
            return None;
        }
        Some(Self {
            packed: packing.packed.into_boxed_slice(),
            length,
            others: packing.others.into_boxed_slice(),
            lowercase: packing.lowercase.into_boxed_slice(),
        })
    }

    /// The number of letters.
    pub fn len(&self) -> usize {
        self.length as usize
    }

    /// True for a record without letters.
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// The bytes held (packed bases, other runs, intervals), for the memory accounting of a caller.
    pub fn heap_bytes(&self) -> usize {
        self.packed.len()
            + std::mem::size_of_val(&*self.others)
            + std::mem::size_of_val(&*self.lowercase)
    }

    /// Clears `out` and writes the letters as they were packed.
    pub fn unpack_into(&self, out: &mut Vec<u8>) {
        let length = self.len();
        out.clear();
        out.reserve_exact(length);
        let whole = length / 4;
        out.extend(
            self.packed[..whole]
                .iter()
                .flat_map(|&byte| LETTERS[usize::from(byte)]),
        );
        if length > whole * 4 {
            out.extend_from_slice(&LETTERS[usize::from(self.packed[whole])][..length - whole * 4]);
        }
        for run in self.others.iter() {
            let start = run.offset as usize;
            out[start..start + run.length as usize].fill(run.byte);
        }
        for &(start, end) in self.lowercase.iter() {
            out[start as usize..end as usize].make_ascii_lowercase();
        }
    }
}

/// Packs the leading blocks of `letters` that are all `ACGT` or all `acgt` (in the case of the byte
/// before them) and passes the others byte by byte; returns how many bytes it took (a multiple of
/// the block, the rest is the tail).
fn pack_blocks(letters: &[u8], packing: &mut Packing) -> usize {
    #[cfg(target_arch = "x86_64")]
    {
        // SAFETY: SSE2 is part of the x86_64 baseline.
        unsafe { pack_blocks_sse2(letters, packing) }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        pack_blocks_swar(letters, packing)
    }
}

/// `pack_blocks` 8 bytes at a time in a `u64` (every target; the oracle of the SSE2 form).
fn pack_blocks_swar(letters: &[u8], packing: &mut Packing) -> usize {
    let blocks = letters.len() / 8;
    for block in 0..blocks {
        let start = block * 8;
        let bytes: [u8; 8] = letters[start..start + 8].try_into().expect("8 bytes");
        match pack8(u64::from_le_bytes(bytes), packing.case_bits()) {
            Some(packed) => packing.packed[start / 4..start / 4 + 2].copy_from_slice(&packed),
            None => {
                for (index, &byte) in bytes.iter().enumerate() {
                    packing.byte(start + index, byte);
                }
            }
        }
    }
    blocks * 8
}

/// `pack_blocks` 16 bytes at a time with SSE2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn pack_blocks_sse2(letters: &[u8], packing: &mut Packing) -> usize {
    use std::arch::x86_64::*;
    let case_bit = _mm_set1_epi8(0x20);
    let [a, c, g, t] = [b'A', b'C', b'G', b'T'].map(|letter| _mm_set1_epi8(letter as i8));
    let threes = _mm_set1_epi8(3);
    let ones = _mm_set1_epi8(1);
    let low_bytes = _mm_set1_epi16(0x00FF);
    let low_halves = _mm_set1_epi32(0x0000_FFFF);
    let blocks = letters.len() / 16;
    for block in 0..blocks {
        let start = block * 16;
        let bytes = &letters[start..start + 16];
        let x = _mm_loadu_si128(bytes.as_ptr() as *const __m128i);
        let upper = _mm_andnot_si128(case_bit, x);
        let acgt = _mm_or_si128(
            _mm_or_si128(_mm_cmpeq_epi8(upper, a), _mm_cmpeq_epi8(upper, c)),
            _mm_or_si128(_mm_cmpeq_epi8(upper, g), _mm_cmpeq_epi8(upper, t)),
        );
        let case = if packing.lower_start.is_some() {
            case_bit
        } else {
            _mm_setzero_si128()
        };
        let same_case = _mm_cmpeq_epi8(_mm_and_si128(x, case_bit), case);
        if _mm_movemask_epi8(_mm_and_si128(acgt, same_case)) != 0xFFFF {
            for (index, &byte) in bytes.iter().enumerate() {
                packing.byte(start + index, byte);
            }
            continue;
        }
        // the codes of `pack8`, then each 16-bit lane's two codes and each 32-bit lane's four
        let codes = _mm_xor_si128(
            _mm_and_si128(_mm_srli_epi16(upper, 1), threes),
            _mm_and_si128(_mm_srli_epi16(upper, 2), ones),
        );
        let pairs = _mm_or_si128(
            _mm_slli_epi16(_mm_and_si128(codes, low_bytes), 2),
            _mm_srli_epi16(codes, 8),
        );
        let quads = _mm_or_si128(
            _mm_slli_epi32(_mm_and_si128(pairs, low_halves), 4),
            _mm_srli_epi32(pairs, 16),
        );
        let words = _mm_packs_epi32(quads, quads);
        let packed = _mm_cvtsi128_si32(_mm_packus_epi16(words, words)) as u32;
        packing.packed[start / 4..start / 4 + 4].copy_from_slice(&packed.to_le_bytes());
    }
    blocks * 16
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pack_with(letters: &[u8], blocks: fn(&[u8], &mut Packing) -> usize) -> Packing {
        let mut packing = Packing {
            packed: vec![0u8; letters.len().div_ceil(4)],
            others: Vec::new(),
            lowercase: Vec::new(),
            lower_start: None,
        };
        let done = blocks(letters, &mut packing);
        for (index, &byte) in letters.iter().enumerate().skip(done) {
            packing.byte(index, byte);
        }
        packing
    }

    fn bytewise(_: &[u8], _: &mut Packing) -> usize {
        0
    }

    /// A small xorshift, so the cases are the same on every run.
    struct Draws(u64);

    impl Draws {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
    }

    fn cases() -> Vec<Vec<u8>> {
        let mut draws = Draws(0x9E37_79B9_7F4A_7C15);
        let mut cases: Vec<Vec<u8>> = vec![Vec::new(), b"A".to_vec(), b"acgtACGTn".to_vec()];
        let others = b"NnRYKMSWBDHVrykmswbdhv-*Xx\x00\xff \x7f";
        for case in 0..3000 {
            let length = (draws.next() % 90) as usize;
            let letters = (0..length)
                .map(|_| match draws.next() % 8 {
                    0 => others[(draws.next() % others.len() as u64) as usize],
                    _ if case % 3 == 0 => b"acgt"[(draws.next() & 3) as usize],
                    _ if case % 3 == 1 && draws.next() % 16 == 0 => {
                        b"acgt"[(draws.next() & 3) as usize]
                    }
                    _ => b"ACGT"[(draws.next() & 3) as usize],
                })
                .collect();
            cases.push(letters);
        }
        // long runs of each case, as a reader gives them
        let mut long: Vec<u8> = (0..5000)
            .map(|_| b"ACGT"[(draws.next() & 3) as usize])
            .collect();
        long[1000..1037].fill(b'N');
        long[2001..2300].make_ascii_lowercase();
        cases.push(long);
        cases
    }

    #[test]
    fn the_letters_are_made_again_exactly() {
        for letters in cases() {
            let mut out = b"left over".to_vec();
            match PackedLetters::pack(&letters) {
                Some(packed) => {
                    assert_eq!(packed.len(), letters.len());
                    packed.unpack_into(&mut out);
                    assert_eq!(out, letters);
                }
                None => assert!(!letters.is_empty(), "an empty record packs: {letters:?}"),
            }
        }
    }

    #[test]
    fn every_block_form_packs_as_byte_by_byte() {
        for letters in cases() {
            let expected = pack_with(&letters, bytewise);
            let swar = pack_with(&letters, pack_blocks_swar);
            assert_eq!(swar.packed, expected.packed, "{letters:?}");
            assert_eq!(swar.others, expected.others, "{letters:?}");
            assert_eq!(swar.lowercase, expected.lowercase, "{letters:?}");
            #[cfg(target_arch = "x86_64")]
            {
                let sse2 = pack_with(&letters, |letters, packing| unsafe {
                    pack_blocks_sse2(letters, packing)
                });
                assert_eq!(sse2.packed, expected.packed, "{letters:?}");
                assert_eq!(sse2.others, expected.others, "{letters:?}");
                assert_eq!(sse2.lowercase, expected.lowercase, "{letters:?}");
            }
        }
    }

    #[test]
    fn a_record_of_other_bytes_keeps_its_letters() {
        let letters: Vec<u8> = (0..400).map(|index| b"NR"[index % 2]).collect();
        assert_eq!(PackedLetters::pack(&letters), None);
        let packed = PackedLetters::pack(&[b'N'; 400]).expect("one run");
        assert_eq!(packed.heap_bytes(), 100 + 12);
    }
}
