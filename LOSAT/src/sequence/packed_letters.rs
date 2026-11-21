//! A record's letters held as 2-bit codes of `A`, `C`, `G` and `T`, the runs of other bytes and the
//! intervals of lower-case bytes, from which the letters are made again exactly (any bytes).
//!
//! No NCBI counterpart: memory only. The query records that wait for their batch keep about a
//! quarter of a byte per letter instead of one (EXPERIMENT LOSAT_X_BNQPACK,
//! `algorithm/blastn/blast_engine/x_query_drop.rs`). Unlike `packed_subject::PackedNucleotide`, no
//! base is drawn for the other bytes and nothing reads the codes: only the letters are made from
//! the store, so no value of a search depends on how they are held.
//!
//! Speed: packing reads each letter once, 32 at a time with AVX2 where the CPU has it (checked at
//! run time) and 8 at a time in a `u64` elsewhere; making the letters again writes 32 from 8 packed
//! bytes with AVX2, else 4 from a table per packed byte. A block that is not all `ACGT` (or all
//! `acgt` inside a lower-case interval) goes byte by byte. Every form gives the same store
//! (`tests::every_form_packs_as_byte_by_byte`).

/// One run of equal bytes that are none of `A`, `C`, `G`, `T` in either case, upper-cased (the
/// lower-case intervals give the case back).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct OtherRun {
    offset: u32,
    length: u32,
    byte: u8,
}

/// The letters of one record (see the module doc); 32 bytes beside its allocations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedLetters {
    /// 4 bases per byte from the high bits (`A` 0, `C` 1, `G` 2, `T` 3; 0 under an other run),
    /// `length.div_ceil(4)` bytes.
    packed: Box<[u8]>,
    length: u32,
    /// The other runs and lower-case intervals, where the record has any.
    rest: Option<Box<Exceptions>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Exceptions {
    others: Box<[OtherRun]>,
    /// `[start, end)` of the lower-case bytes: the runs found while packing, then the reader's mask
    /// (`pack_masked`), which may meet or cover them (lower-casing twice changes nothing).
    lowercase: Box<[(u32, u32)]>,
}

/// The bytes an allocation of `bytes` takes from the heap: glibc's chunks, 16-byte steps with an
/// 8-byte header and 32 bytes at least (none for an empty one).
fn chunk_bytes(bytes: usize) -> usize {
    if bytes == 0 {
        0
    } else {
        ((bytes + 8 + 15) & !15).max(32)
    }
}

/// What a caller holds beside a store: its own size and an index (a list of the packed records).
const HELD_BESIDE: usize = std::mem::size_of::<PackedLetters>() + std::mem::size_of::<usize>();

/// `0x01` in each byte of a `u64`.
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

/// What a packing pass adds besides the packed bytes.
#[derive(Default)]
struct Rest {
    others: Vec<OtherRun>,
    lowercase: Vec<(u32, u32)>,
    lower_start: Option<usize>,
}

impl Rest {
    /// Letter `index` of the record, one byte at a time, into bit pair `at % 4` of `dst[at / 4]`
    /// (zero before).
    #[inline(always)]
    fn byte(&mut self, index: usize, byte: u8, dst: &mut [u8], at: usize) {
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
        dst[at / 4] |= code << (6 - 2 * (at % 4));
    }
}

/// The 2 packed bytes of 8 letters that are all `A`, `C`, `G`, `T` with the case bits `case`
/// (`0` or `0x20` in each byte), or `None`.
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
    /// The store of `letters`, or `None` where it and its place in a caller's list would take no
    /// fewer heap bytes than the letters' own allocation (no letters, a short record, many other
    /// runs or case changes; `chunk_bytes`) or the record is longer than a `u32` offset.
    pub fn pack(letters: &[u8]) -> Option<Self> {
        Self::pack_masked(letters, &[])
    }

    /// The store of `letters` with the closed ranges `masked` (in order) lower-cased after them: a
    /// reader's upper-case letters and its lowercase mask (`fasta_reader`'s records lower-case
    /// those letters), without making the lower-cased copy.
    pub fn pack_masked(letters: &[u8], masked: &[(u32, u32)]) -> Option<Self> {
        let length = u32::try_from(letters.len()).ok()?;
        let mut packed = vec![0u8; letters.len().div_ceil(4)];
        let mut rest = Rest::default();
        let done = pack_blocks(letters, &mut packed, &mut rest);
        for (index, &byte) in letters.iter().enumerate().skip(done) {
            rest.byte(index, byte, &mut packed, index);
        }
        if let Some(start) = rest.lower_start {
            rest.lowercase.push((start as u32, length));
        }
        for &(from, to) in masked {
            debug_assert!(from <= to && to < length, "a mask out of the record");
            rest.lowercase.push((from, to + 1));
        }
        let exceptions = !(rest.others.is_empty() && rest.lowercase.is_empty());
        let held = HELD_BESIDE
            + chunk_bytes(packed.len())
            + if exceptions {
                chunk_bytes(std::mem::size_of::<Exceptions>())
                    + chunk_bytes(std::mem::size_of_val(rest.others.as_slice()))
                    + chunk_bytes(std::mem::size_of_val(rest.lowercase.as_slice()))
            } else {
                0
            };
        if held >= chunk_bytes(letters.len()) {
            return None;
        }
        Some(Self {
            packed: packed.into_boxed_slice(),
            length,
            rest: exceptions.then(|| {
                Box::new(Exceptions {
                    others: rest.others.into_boxed_slice(),
                    lowercase: rest.lowercase.into_boxed_slice(),
                })
            }),
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
            + self.rest.as_ref().map_or(0, |rest| {
                std::mem::size_of_val(&*rest.others) + std::mem::size_of_val(&*rest.lowercase)
            })
    }

    /// Clears `out` and writes the letters as they were packed.
    pub fn unpack_into(&self, out: &mut Vec<u8>) {
        let length = self.len();
        out.clear();
        out.reserve_exact(length);
        let spare = &mut out.spare_capacity_mut()[..length];
        let done = unpack_blocks(&self.packed, spare);
        for (index, slot) in spare.iter_mut().enumerate().skip(done) {
            let code = (self.packed[index / 4] >> (6 - 2 * (index % 4))) & 3;
            slot.write(b"ACGT"[usize::from(code)]);
        }
        // SAFETY: `unpack_blocks` wrote the first `done` letters and the loop the others.
        unsafe { out.set_len(length) };
        let Some(rest) = &self.rest else {
            return;
        };
        for run in rest.others.iter() {
            let start = run.offset as usize;
            out[start..start + run.length as usize].fill(run.byte);
        }
        for &(start, end) in rest.lowercase.iter() {
            out[start as usize..end as usize].make_ascii_lowercase();
        }
    }

    /// EXPERIMENT (LOSAT_X_BNQPARTS): the bytes `append_to` appends.
    pub fn stored_bytes(&self) -> usize {
        self.packed.len()
            + self
                .rest
                .as_ref()
                .map_or(0, |rest| 9 * rest.others.len() + 8 * rest.lowercase.len())
    }

    /// EXPERIMENT (LOSAT_X_BNQPARTS): appends the store to `out` (the packed bases, then the other
    /// runs and the lower-case intervals), for `unpack_from`.
    pub fn append_to(&self, out: &mut Vec<u8>) -> PackedAt {
        out.extend_from_slice(&self.packed);
        let (others, lowercase) = match &self.rest {
            Some(rest) => {
                for run in rest.others.iter() {
                    out.extend_from_slice(&run.offset.to_le_bytes());
                    out.extend_from_slice(&run.length.to_le_bytes());
                    out.push(run.byte);
                }
                for &(from, to) in rest.lowercase.iter() {
                    out.extend_from_slice(&from.to_le_bytes());
                    out.extend_from_slice(&to.to_le_bytes());
                }
                (rest.others.len(), rest.lowercase.len())
            }
            None => (0, 0),
        };
        PackedAt {
            length: self.length,
            // A store's runs and intervals are fewer than its letters (`u32`).
            others: others as u32,
            lowercase: lowercase as u32,
        }
    }
}

/// EXPERIMENT (LOSAT_X_BNQPARTS): where a record's store was written into a caller's buffer
/// (`PackedLetters::append_to`): the letters' count and those of the other runs and lower-case intervals written
/// after the packed bases (9 and 8 bytes each, little-endian `u32`s and the run's byte).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PackedAt {
    length: u32,
    others: u32,
    lowercase: u32,
}

impl PackedAt {
    /// The number of letters.
    pub fn len(&self) -> usize {
        self.length as usize
    }

    /// True for a record without letters.
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// The bytes the store takes in its buffer.
    pub fn bytes(&self) -> usize {
        self.len().div_ceil(4) + 9 * self.others as usize + 8 * self.lowercase as usize
    }
}

/// EXPERIMENT (LOSAT_X_BNQPARTS): clears `out` and writes the letters of `store`, the bytes that
/// `PackedLetters::append_to` appended for `at`.
pub fn unpack_from(store: &[u8], at: PackedAt, out: &mut Vec<u8>) {
    let length = at.len();
    let (packed, mut rest) = store.split_at(length.div_ceil(4));
    out.clear();
    out.reserve_exact(length);
    let spare = &mut out.spare_capacity_mut()[..length];
    let done = unpack_blocks(packed, spare);
    for (index, slot) in spare.iter_mut().enumerate().skip(done) {
        let code = (packed[index / 4] >> (6 - 2 * (index % 4))) & 3;
        slot.write(b"ACGT"[usize::from(code)]);
    }
    // SAFETY: `unpack_blocks` wrote the first `done` letters and the loop the others.
    unsafe { out.set_len(length) };
    let word = |bytes: &[u8]| u32::from_le_bytes(bytes[..4].try_into().expect("4 bytes")) as usize;
    for _ in 0..at.others {
        let (start, run) = (word(rest), word(&rest[4..]));
        out[start..start + run].fill(rest[8]);
        rest = &rest[9..];
    }
    for _ in 0..at.lowercase {
        out[word(rest)..word(&rest[4..])].make_ascii_lowercase();
        rest = &rest[8..];
    }
}

/// Packs the leading blocks of `letters` into `packed` (blocks that are not all `ACGT` or all
/// `acgt` in the case of the letter before them go byte by byte into `rest`); returns how many
/// letters it took (the rest is the tail).
fn pack_blocks(letters: &[u8], packed: &mut [u8], rest: &mut Rest) -> usize {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU has AVX2.
        return unsafe { pack_blocks_avx2(letters, packed, rest) };
    }
    pack_blocks_swar(letters, packed, rest)
}

/// `pack_blocks` 8 letters at a time in a `u64` (every target; the oracle of the AVX2 form).
fn pack_blocks_swar(letters: &[u8], packed: &mut [u8], rest: &mut Rest) -> usize {
    let mut done = 0;
    for (block, (src, dst)) in letters
        .chunks_exact(8)
        .zip(packed.chunks_exact_mut(2))
        .enumerate()
    {
        let bytes: [u8; 8] = src.try_into().expect("8 letters");
        let case = if rest.lower_start.is_some() {
            0x20 * ONES
        } else {
            0
        };
        match pack8(u64::from_le_bytes(bytes), case) {
            Some(pair) => dst.copy_from_slice(&pair),
            None => {
                for (at, &byte) in bytes.iter().enumerate() {
                    rest.byte(block * 8 + at, byte, dst, at);
                }
            }
        }
        done += 8;
    }
    done
}

/// `pack_blocks` 32 letters at a time with AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn pack_blocks_avx2(letters: &[u8], packed: &mut [u8], rest: &mut Rest) -> usize {
    use std::arch::x86_64::*;
    let case_bit = _mm256_set1_epi8(0x20);
    let threes = _mm256_set1_epi8(3);
    let ones = _mm256_set1_epi8(1);
    // the letter of each code, in each 128-bit lane (`_mm256_shuffle_epi8` reads its lane)
    let acgt = _mm256_setr_epi8(
        65, 67, 71, 84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 65, 67, 71, 84, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0,
    );
    // byte weights [4, 1] (two codes per 16-bit lane), then 16-bit weights [16, 1] (four per 32)
    let pair_weights = _mm256_set1_epi16(0x0104);
    let quad_weights = _mm256_set1_epi32(0x0001_0010);
    let mut case = if rest.lower_start.is_some() {
        case_bit
    } else {
        _mm256_setzero_si256()
    };
    let mut done = 0;
    for (block, (src, dst)) in letters
        .chunks_exact(32)
        .zip(packed.chunks_exact_mut(8))
        .enumerate()
    {
        let x = _mm256_loadu_si256(src.as_ptr() as *const __m256i);
        let upper = _mm256_andnot_si256(case_bit, x);
        // the codes of `pack8`; a letter is valid where the code's letter is the letter
        let codes = _mm256_xor_si256(
            _mm256_and_si256(_mm256_srli_epi16(upper, 1), threes),
            _mm256_and_si256(_mm256_srli_epi16(upper, 2), ones),
        );
        let valid = _mm256_and_si256(
            _mm256_cmpeq_epi8(_mm256_shuffle_epi8(acgt, codes), upper),
            _mm256_cmpeq_epi8(_mm256_and_si256(x, case_bit), case),
        );
        if _mm256_movemask_epi8(valid) == -1 {
            let quads = _mm256_madd_epi16(_mm256_maddubs_epi16(codes, pair_weights), quad_weights);
            let words = _mm256_packus_epi32(quads, quads);
            let bytes = _mm256_packus_epi16(words, words);
            let low = _mm_cvtsi128_si32(_mm256_castsi256_si128(bytes)) as u32;
            let high = _mm_cvtsi128_si32(_mm256_extracti128_si256(bytes, 1)) as u32;
            dst.copy_from_slice(&(u64::from(low) | u64::from(high) << 32).to_le_bytes());
        } else {
            for (at, &byte) in src.iter().enumerate() {
                rest.byte(block * 32 + at, byte, dst, at);
            }
            case = if rest.lower_start.is_some() {
                case_bit
            } else {
                _mm256_setzero_si256()
            };
        }
        done += 32;
    }
    done
}

/// Writes the letters of the leading whole blocks of `packed` into `out` (as many letters as
/// `out` holds at most); returns how many it wrote.
fn unpack_blocks(packed: &[u8], out: &mut [std::mem::MaybeUninit<u8>]) -> usize {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU has AVX2.
        return unsafe { unpack_blocks_avx2(packed, out) };
    }
    let mut done = 0;
    for (dst, &byte) in out.chunks_exact_mut(4).zip(packed) {
        for (slot, &letter) in dst.iter_mut().zip(&LETTERS[usize::from(byte)]) {
            slot.write(letter);
        }
        done += 4;
    }
    done
}

/// `unpack_blocks` 32 letters from 8 packed bytes at a time with AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn unpack_blocks_avx2(packed: &[u8], out: &mut [std::mem::MaybeUninit<u8>]) -> usize {
    use std::arch::x86_64::*;
    // each packed byte to the 4 bytes of its letters (each lane holds all 8 packed bytes)
    let spread = _mm256_setr_epi8(
        0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7,
        7, 7,
    );
    let high_bits = _mm256_set1_epi32(i32::from_le_bytes([0x80, 0x20, 0x08, 0x02]));
    let low_bits = _mm256_set1_epi32(i32::from_le_bytes([0x40, 0x10, 0x04, 0x01]));
    let twos = _mm256_set1_epi8(2);
    let ones = _mm256_set1_epi8(1);
    let acgt = _mm256_setr_epi8(
        65, 67, 71, 84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 65, 67, 71, 84, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0,
    );
    let mut done = 0;
    for (src, dst) in packed.chunks_exact(8).zip(out.chunks_exact_mut(32)) {
        let word = i64::from_le_bytes(src.try_into().expect("8 packed bytes"));
        let bytes = _mm256_shuffle_epi8(_mm256_set1_epi64x(word), spread);
        let high = _mm256_and_si256(
            _mm256_cmpeq_epi8(_mm256_and_si256(bytes, high_bits), high_bits),
            twos,
        );
        let low = _mm256_and_si256(
            _mm256_cmpeq_epi8(_mm256_and_si256(bytes, low_bits), low_bits),
            ones,
        );
        let letters = _mm256_shuffle_epi8(acgt, _mm256_or_si256(high, low));
        _mm256_storeu_si256(dst.as_mut_ptr() as *mut __m256i, letters);
        done += 32;
    }
    done
}

#[cfg(test)]
#[path = "packed_letters_tests.rs"]
mod tests;
