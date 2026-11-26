//! `CSeqConvert::Pack` of IUPACna residues into the segments of a delta sequence
//! (`CDelta_ext::AddAndSplit`), and the conversions of the residues into NCBI2na,
//! NCBI4na and NCBI8na (`c++/src/util/sequtil/`).

use crate::utils::byte_search::find_first;

/// The codings that `Pack` chooses for IUPACna residues without gaps, and NCBIeaa, the one
/// literal of a protein's residues (`CSeqConvert::Pack` copies NCBIeaa when the target takes
/// no gaps).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Coding {
    Ncbi2na,
    Ncbi4na,
    Ncbieaa,
}

/// One segment of the packed residues (`IPackTarget::NewSegment`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) struct Segment {
    pub(crate) coding: Coding,
    pub(crate) start: u32,
    pub(crate) length: u32,
}

/// `kBestCodingsWithoutGaps.iupacna`: NCBI2na for `A`, `C`, `G`, `T`, `U` in either
/// case, NCBI4na for every other byte.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_best_codings.inc:34-62
/// The table without gaps is the one used: `CDelta_ext_PackTarget::GapsOK` is the
/// `gaps_ok` of `AddAndSplit`, which `CFastaReader::AssembleSeq` passes as false.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1064-1071
const fn best_coding_iupacna(residue: u8) -> Coding {
    match residue {
        b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u' => Coding::Ncbi2na,
        _ => Coding::Ncbi4na,
    }
}

/// 1 for bytes whose best coding without gaps is NCBI4na, tabulating [`best_coding_iupacna`]
/// for short-run scans in `pack_iupacna`.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1064-1071
const BEST_IS_NCBI4NA: [u8; 256] = {
    let mut table = [0; 256];
    let mut index = 0;
    while index < 256 {
        table[index] = matches!(best_coding_iupacna(index as u8), Coding::Ncbi4na) as u8;
        index += 1;
    }
    table
};

/// Whether the best coding of `residue` without gaps is NCBI2na (`A`, `C`, `G`, `T`, `U` in
/// either case), following `kBestCodingsWithoutGaps` for vectorized scans in `pack_iupacna`
/// (the unit test checks every byte). Clearing bit 5 folds lower case to upper case without
/// mapping any non-letter onto `A`–`Z`, and clearing bit 0 pairs `T` (`0x54`) with `U` (`0x55`).
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_best_codings.inc:34-62
const fn is_ncbi2na_best(residue: u8) -> bool {
    let upper = residue & 0xdf;
    (upper == b'A') | (upper == b'C') | (upper == b'G') | ((upper & 0xfe) == b'T')
}

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_shared.cpp:195-203
/// The sums are made in 64 bits. NCBI's `TSeqPos` sum wraps only for a run of
/// 4,294,967,293 letters or more, which leaves at most two other letters in the piece:
/// NCBI and LOSAT then both keep the run as one segment, which `add_and_split` refuses.
fn bytes_needed(coding: Coding, length: u32) -> u64 {
    match coding {
        Coding::Ncbi2na => u64::from(length).div_ceil(4),
        Coding::Ncbi4na => u64::from(length).div_ceil(2),
        Coding::Ncbieaa => u64::from(length),
    }
}

/// `CDelta_ext_PackTarget::GetOverhead` of a segment with data.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Delta_ext.cpp:126-129
const OVERHEAD: u64 = 128;

/// `CPacker::x_GetWideCoding` for `e_Iupacna`.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1292-1296
const WIDE_CODING: Coding = Coding::Ncbi4na;

/// A boundary's bit: its new coding (the coding of the segment it starts in the narrow
/// arrangement) is NCBI4na, not NCBI2na.
const STEP_NEW_WIDE: u8 = 1;
/// A boundary's bit: the narrow arrangement that ends here continues the wide one that
/// ended at the boundary before (`m_EndingNarrow` takes a copy of `m_EndingWide`).
const STEP_NARROW_FROM_WIDE: u8 = 2;
/// A boundary's bit: the wide arrangement that ends here continues the narrow one
/// (`m_EndingWide` takes a copy of `m_EndingNarrow`), so it continues what that one continues.
const STEP_WIDE_FROM_NARROW: u8 = 4;

/// `CSeqConvert_imp::CPacker` for IUPACna. NCBI keeps the codings of its two arrangements
/// (`SArrangement`) in blocks of 16 that they share up to their last common block; LOSAT
/// keeps, for each boundary, which arrangement each of the two that end there continues
/// (one byte), and reads the best one's codings back from the last boundary. The two
/// make the same choices; with the boundary's position LOSAT keeps 5 bytes for each
/// boundary, where a list of nodes (one coding and the index of the one before it, for
/// each arrangement) and the codings read back from it took 21.
struct Packer {
    /// `m_Boundaries`.
    boundaries: Vec<u32>,
    /// The `STEP_*` bits of each boundary but the last.
    steps: Vec<u8>,
    /// The costs of `m_EndingNarrow` and `m_EndingWide`.
    narrow_cost: u64,
    wide_cost: u64,
    /// The last coding of `m_EndingNarrow` (the new coding of the boundary before); the
    /// last coding of `m_EndingWide` is always NCBI4na.
    last_narrow: Coding,
}

impl Packer {
    /// `new_coding` is `None` for `kNoCoding` (the end of the residues).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1194-1249
    fn add_boundary(&mut self, pos: u32, new_coding: Option<Coding>) {
        let new_wide = |coding: Coding| {
            if coding == WIDE_CODING {
                STEP_NEW_WIDE
            } else {
                0
            }
        };
        let Some(&last_boundary) = self.boundaries.last() else {
            let Some(new_coding) = new_coding else {
                return;
            };
            self.boundaries.push(pos);
            self.steps.push(new_wide(new_coding));
            self.last_narrow = new_coding;
            self.wide_cost = OVERHEAD;
            self.narrow_cost = OVERHEAD;
            return;
        };
        let last_length = pos - last_boundary;
        self.boundaries.push(pos);

        self.narrow_cost += bytes_needed(self.last_narrow, last_length);
        self.wide_cost += bytes_needed(WIDE_CODING, last_length);

        let Some(new_coding) = new_coding else {
            return;
        };

        let mut step = new_wide(new_coding);
        if new_coding != WIDE_CODING && self.narrow_cost > self.wide_cost {
            self.narrow_cost = self.wide_cost;
            step |= STEP_NARROW_FROM_WIDE;
        }

        let alt_wide_cost = self.narrow_cost + OVERHEAD;
        self.narrow_cost += OVERHEAD;

        if self.wide_cost > alt_wide_cost {
            self.wide_cost = alt_wide_cost;
            step |= STEP_WIDE_FROM_NARROW;
        } else if new_coding == WIDE_CODING {
            self.narrow_cost = self.wide_cost;
            step |= STEP_NARROW_FROM_WIDE;
        }

        self.last_narrow = new_coding;
        self.steps.push(step);
    }
}

/// `CSeqConvert::Pack` of `src` and its length as `e_Iupacna` into a `CDelta_ext_PackTarget`: the
/// segments that `CDelta_ext::AddAndSplit` appends for the residues `src` (at most
/// 4,294,967,295 of them).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Delta_ext.cpp:201-205
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1058-1063,1105-1107
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1124-1192
/// A boundary is where the best coding of the residues changes (`m_SrcDensity` is 1 for
/// IUPACna). The two arrangements cost the same when the residues are one run of one or
/// two unambiguous bases (128 + 1 byte in either coding), and the comparison then takes
/// the wide one: such a sequence is NCBI4na. The segments' bytes are made by the caller
/// (`iupacna_to_ncbi2na`, `iupacna_to_ncbi4na`).
pub(super) fn pack_iupacna(src: &[u8]) -> Vec<Segment> {
    if src.is_empty() {
        return Vec::new();
    }
    let mut packer = Packer {
        boundaries: Vec::new(),
        steps: Vec::new(),
        narrow_cost: 0,
        wide_cost: 0,
        last_narrow: WIDE_CODING,
    };
    let mut pos = 0;
    while pos < src.len() {
        let wide = BEST_IS_NCBI4NA[src[pos] as usize];
        let coding = if wide == 1 {
            Coding::Ncbi4na
        } else {
            Coding::Ncbi2na
        };
        packer.add_boundary(pos as u32, Some(coding));
        let rest = &src[pos + 1..];
        let diff = |&residue: &u8| BEST_IS_NCBI4NA[residue as usize] != wide;
        if rest.len() < 20 {
            pos += 1 + rest.iter().position(diff).unwrap_or(rest.len());
            continue;
        }
        if let Some(at) = rest[..4].iter().position(diff) {
            pos += 1 + at;
            continue;
        }
        let tail = &rest[4..];
        // The run of residues of the same best coding, found with one search (LOSAT's speed).
        let run = if wide == 0 {
            find_first(tail, |residue| !is_ncbi2na_best(residue))
        } else {
            find_first(tail, is_ncbi2na_best)
        };
        pos += 5 + run.unwrap_or(tail.len());
    }
    packer.add_boundary(src.len() as u32, None);

    // The best arrangement's codings, from the last segment back: each step becomes 1
    // where the segment it starts is NCBI4na.
    let mut narrow = packer.narrow_cost < packer.wide_cost;
    for step in packer.steps.iter_mut().rev() {
        let bits = *step;
        let wide = !narrow || bits & STEP_NEW_WIDE != 0;
        narrow = if narrow {
            bits & STEP_NARROW_FROM_WIDE == 0
        } else {
            bits & STEP_WIDE_FROM_NARROW != 0 && bits & STEP_NARROW_FROM_WIDE == 0
        };
        *step = u8::from(wide);
    }
    let coding_of = |wide: u8| {
        if wide == 1 {
            Coding::Ncbi4na
        } else {
            Coding::Ncbi2na
        }
    };

    let steps = &packer.steps;
    let n = steps.len();
    let mut segments = Vec::new();
    let mut i = 0;
    while i < n {
        let wide = steps[i];
        let start = packer.boundaries[i];
        while i < n - 1 && steps[i + 1] == wide {
            i += 1;
        }
        segments.push(Segment {
            coding: coding_of(wide),
            start,
            length: packer.boundaries[i + 1] - start,
        });
        i += 1;
    }
    segments
}

/// The NCBI2na code of the unambiguous upper-case IUPACna letters (`A`: 0, `C`: 1, `G`: 2, `T`
/// and `U`: 3). LOSAT returns 0 for every other byte: NCBI's `CIupacnaTo2na` table also gives
/// the ambiguous letters a code (`B`/`M`/`S`: 1, `D`/`K`/`R`: 2, `W`/`Y`: 3, `H`/`N`/`V`: 0)
/// and maps lower case to the same code, but only upper-case `A`, `C`, `G`, `T` and `U` reach an
/// NCBI2na segment. Written without a branch so a loop over residues computes the codes of many
/// at once.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:79-351
pub(crate) const fn ncbi2na_code(residue: u8) -> u8 {
    let c = (residue == b'C') as u8;
    let g = (residue == b'G') as u8;
    let t = ((residue == b'T') | (residue == b'U')) as u8;
    c | (g << 1) | (t * 3)
}

/// The NCBI4na code of the IUPACna letters (`CIupacnaTo4na`'s rows: the table has the
/// code in the high and in the low half of a byte), which is also their NCBI8na code. The
/// table has the codes of the letters in either case; every other row is 15 (`N`), the
/// hyphen of a reader without `fParseGaps` among them.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:392-394,409-410
pub(crate) const fn ncbi4na_code(residue: u8) -> u8 {
    match residue.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' | b'U' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

const NCBI4NA_CODES: [u8; 256] = {
    let mut table = [0; 256];
    let mut index = 0;
    while index < 256 {
        table[index] = ncbi4na_code(index as u8);
        index += 1;
    }
    table
};

/// Lookup table of [`ncbi2na_code`] for short slices in [`iupacna_to_ncbi2na`], matching
/// `CIupacnaTo2na` on the unambiguous upper-case letters that reach an NCBI2na segment.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:79-351
const NCBI2NA_CODES: [u8; 256] = {
    let mut table = [0; 256];
    let mut index = 0;
    while index < 256 {
        table[index] = ncbi2na_code(index as u8);
        index += 1;
    }
    table
};

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:257-289
/// Short slices convert directly through `NCBI2NA_CODES` to avoid zeroing a stack block; longer
/// slices compute the codes of a 256-residue block without branches and pack four to a byte
/// ([`iupacna_to_ncbi2na_blocked`], [`pack_four_codes`]).
pub(crate) fn iupacna_to_ncbi2na(src: &[u8]) -> Vec<u8> {
    const BLOCK: usize = 256;
    let mut dst = Vec::with_capacity(src.len().div_ceil(4));
    if src.len() < 32 {
        let mut chunks = src.chunks_exact(4);
        dst.extend((&mut chunks).map(|four| {
            NCBI2NA_CODES[four[0] as usize] << 6
                | NCBI2NA_CODES[four[1] as usize] << 4
                | NCBI2NA_CODES[four[2] as usize] << 2
                | NCBI2NA_CODES[four[3] as usize]
        }));
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut byte = 0;
            for (i, &residue) in rest.iter().enumerate() {
                byte |= NCBI2NA_CODES[residue as usize] << (6 - 2 * i);
            }
            dst.push(byte);
        }
        return dst;
    }
    iupacna_to_ncbi2na_blocked::<BLOCK>(src, &mut dst);
    dst
}

/// Blocked path of [`iupacna_to_ncbi2na`], implementing the 4-residue packing loop and overhang
/// of `CSeqConvert_imp::x_ConvertIupacnaTo2na` via a stack block of codes.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:257-289
#[inline(never)]
fn iupacna_to_ncbi2na_blocked<const BLOCK: usize>(src: &[u8], dst: &mut Vec<u8>) {
    let mut codes = [0u8; BLOCK];
    for block in src.chunks(BLOCK) {
        let codes = &mut codes[..block.len()];
        for (code, &residue) in codes.iter_mut().zip(block) {
            *code = ncbi2na_code(residue);
        }
        let (fours, rest) = codes.as_chunks::<4>();
        dst.extend(fours.iter().map(|&four| pack_four_codes(four)));
        if !rest.is_empty() {
            let mut byte = 0;
            for (i, &code) in rest.iter().enumerate() {
                byte |= code << (6 - 2 * i);
            }
            dst.push(byte);
        }
    }
}

/// The byte of four NCBI2na codes (each 0 to 3), the first in the high bits:
/// `c0 << 6 | c1 << 4 | c2 << 2 | c3`, matching the four-table-lookup OR per output byte of
/// `CSeqConvert_imp::x_ConvertIupacnaTo2na`. With the codes as the bytes of a little-endian word
/// `w` (code `k` at bit `8 k`), the product `w * (2^30 + 2^20 + 2^10 + 1)` puts code `k` at bit
/// `30 - 2 k` (bits 24 to 31); the other terms put code `k` at bit `8 k + 30 - 10 j` for `j != k`,
/// below 22 or above 33, where two bits each never carry into bits 24 to 31.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:270-278
fn pack_four_codes(four: [u8; 4]) -> u8 {
    ((u64::from(u32::from_le_bytes(four)) * 0x4010_0401) >> 24) as u8
}

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
pub(crate) fn iupacna_to_ncbi4na(src: &[u8]) -> Vec<u8> {
    let mut dst = Vec::with_capacity(src.len().div_ceil(2));
    if let Some(last) = extend_ncbi4na_pairs(&mut dst, src) {
        dst.push(last << 4);
    }
    dst
}

/// The number of residues that [`extend_ncbi4na_pairs_blocked`] converts per block, matching the
/// pair-wise packing of `CSeqConvert_imp::x_ConvertIupacnaTo4na` and `C8naTo4na`.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:856-874
const NCBI4NA_BLOCK: usize = 64;

/// Appends the NCBI4na bytes for each pair of IUPACna residues in `src` (first in the high half)
/// and returns the code of the trailing residue when `src.len()` is odd, matching
/// `CSeqConvert_imp::x_ConvertIupacnaTo4na` and `s_FixBioseqDeltas`' NCBI8na-to-NCBI4na step.
/// Short slices map directly through `NCBI4NA_CODES` without entering the block frame; longer
/// slices delegate to [`extend_ncbi4na_pairs_blocked`].
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:343-357
fn extend_ncbi4na_pairs(dst: &mut Vec<u8>, src: &[u8]) -> Option<u8> {
    if src.len() >= NCBI4NA_BLOCK {
        return extend_ncbi4na_pairs_blocked(dst, src);
    }
    let mut chunks = src.chunks_exact(2);
    dst.extend(
        (&mut chunks)
            .map(|two| NCBI4NA_CODES[two[0] as usize] << 4 | NCBI4NA_CODES[two[1] as usize]),
    );
    chunks
        .remainder()
        .first()
        .map(|&residue| NCBI4NA_CODES[residue as usize])
}

/// Blocked path of [`extend_ncbi4na_pairs`]: computes codes without table lookups for 64-residue
/// blocks containing only `A`, `C`, `G`, `T`, `U` or `N` (either case) via
/// [`ncbi4na_code_of_base`] and packs them into 32 bytes via [`pair_byte`]; any block with
/// another IUPAC letter or non-letter, plus the remainder, maps directly through `NCBI4NA_CODES`.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:856-874
#[inline(never)]
fn extend_ncbi4na_pairs_blocked(dst: &mut Vec<u8>, src: &[u8]) -> Option<u8> {
    let (blocks, rest) = src.as_chunks::<NCBI4NA_BLOCK>();
    let mut codes = [0u8; NCBI4NA_BLOCK];
    let mut packed = [0u8; NCBI4NA_BLOCK / 2];
    for block in blocks {
        let mut any_zero = 0u8;
        for (code, &residue) in codes.iter_mut().zip(block) {
            let c = ncbi4na_code_of_base(residue);
            *code = c;
            any_zero |= u8::from(c == 0);
        }
        if any_zero == 0 {
            for (out, two) in packed.iter_mut().zip(codes.as_chunks::<2>().0) {
                *out = pair_byte(*two);
            }
            dst.extend_from_slice(&packed);
        } else {
            let (twos, _) = block.as_chunks::<2>();
            dst.extend(twos.iter().map(|&[first, second]| {
                NCBI4NA_CODES[first as usize] << 4 | NCBI4NA_CODES[second as usize]
            }));
        }
    }
    let mut chunks = rest.chunks_exact(2);
    dst.extend(
        (&mut chunks)
            .map(|two| NCBI4NA_CODES[two[0] as usize] << 4 | NCBI4NA_CODES[two[1] as usize]),
    );
    chunks
        .remainder()
        .first()
        .map(|&residue| NCBI4NA_CODES[residue as usize])
}

/// Packs two NCBI4na codes (`0..16`) into one byte with the first in the high nibble, matching
/// the pair-wise OR of `CSeqConvert_imp::x_ConvertIupacnaTo4na` and `C8naTo4na`.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:328-332
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:866-869
#[inline(always)]
fn pair_byte(two: [u8; 2]) -> u8 {
    let pair = u16::from_le_bytes(two);
    ((pair << 4) | (pair >> 8)) as u8
}

/// Computes `CIupacnaTo4na` without a table for `A`, `C`, `G`, `T`, `U` and `N` in either case
/// (`1, 2, 4, 8, 8, 15`), returning `0` for every other byte (every entry of `NCBI4NA_CODES` is
/// non-zero, so `0` signals the table fallback). Clearing bit 5 folds lower case to upper case
/// without mapping any non-letter onto `A`–`Z`, and clearing bit 0 pairs `T` (`0x54`) with `U`
/// (`0x55`).
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:392-457
#[inline(always)]
const fn ncbi4na_code_of_base(residue: u8) -> u8 {
    let upper = residue & 0xdf;
    let a = (upper == b'A') as u8;
    let c = (upper == b'C') as u8;
    let g = (upper == b'G') as u8;
    let t = ((upper & 0xfe) == b'T') as u8;
    let n = (upper == b'N') as u8;
    a | (c << 1) | (g << 2) | (t << 3) | (n * 15)
}

/// `s_FixBioseqDeltas`' `seq8na` and its NCBI4na copy in one: the NCBI4na bytes of the
/// NCBI8na codes appended, two letters to a byte, the first in the high half (`C8naTo4na`).
/// NCBI makes `seq8na` (1 byte a letter) and converts it; LOSAT writes the NCBI4na bytes
/// (half a byte a letter) as the codes come, the same bytes without the `seq8na`.
///
/// The NCBI8na codes of a literal's residues are the NCBI4na codes of its letters, as
/// `s_FixBioseqDeltas` gets them from the literal's NCBI2na or NCBI4na data: an NCBI2na
/// literal holds only `A`, `C`, `G`, `T`, `U`, whose codes 0 to 3 `C2naTo8na` turns into 1,
/// 2, 4, 8, and `C4naTo8na` splits a byte into its two codes
/// (`tests::ncbi8na_of_both_codings`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:343-357
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:856-874
pub(crate) struct Ncbi4naBuilder {
    bytes: Vec<u8>,
    /// The last byte holds one letter, in its high half.
    odd: bool,
}

impl Ncbi4naBuilder {
    /// The builder of `letters` letters, whose bytes are reserved at once. The count can be
    /// one that the input states rather than holds (a gap line), so the allocation can fail
    /// without the process being short of memory for anything else; that failure is
    /// returned (`fix_bioseq_deltas`), and no later growth can fail.
    pub(crate) fn with_letters(letters: usize) -> Result<Self, std::collections::TryReserveError> {
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(letters.div_ceil(2))?;
        Ok(Self { bytes, odd: false })
    }

    fn push_code(&mut self, code: u8) {
        match self.bytes.last_mut() {
            Some(last) if self.odd => *last |= code,
            _ => self.bytes.push(code << 4),
        }
        self.odd = !self.odd;
    }

    /// `run` letters of the NCBI8na code `code` (a gap: 15).
    pub(crate) fn push_run(&mut self, code: u8, mut run: usize) {
        if run > 0 && self.odd {
            self.push_code(code);
            run -= 1;
        }
        self.bytes
            .resize(self.bytes.len() + run / 2, code << 4 | code);
        if run % 2 == 1 {
            self.push_code(code);
        }
    }

    /// The NCBI8na codes of IUPACna residues of a literal.
    pub(crate) fn push_iupacna(&mut self, mut src: &[u8]) {
        if let (true, [first, rest @ ..]) = (self.odd, src) {
            self.push_code(NCBI4NA_CODES[*first as usize]);
            src = rest;
        }
        if let Some(last) = extend_ncbi4na_pairs(&mut self.bytes, src) {
            self.push_code(last);
        }
    }

    /// Appends up to `letters` bases from packed NCBI2na bytes (`CSeqConvert::Convert` from
    /// `e_Ncbi2na` to `e_Ncbi8na`, clamped by `AdjustLength` to `src.len() * 4`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:343-349
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:774-806
    pub(crate) fn push_ncbi2na(&mut self, src: &[u8], letters: usize) {
        let count = letters.min(src.len().saturating_mul(4));
        if count == 0 {
            return;
        }
        let start = if !self.odd {
            let full_bytes = count / 4;
            for &byte in &src[..full_bytes] {
                let b0 = 1u8 << ((byte >> 6) & 3);
                let b1 = 1u8 << ((byte >> 4) & 3);
                let b2 = 1u8 << ((byte >> 2) & 3);
                let b3 = 1u8 << (byte & 3);
                self.bytes.push((b0 << 4) | b1);
                self.bytes.push((b2 << 4) | b3);
            }
            full_bytes * 4
        } else {
            0
        };
        for i in start..count {
            let shift = 6 - 2 * (i & 3);
            let code = 1u8 << ((src[i / 4] >> shift) & 3);
            self.push_code(code);
        }
    }

    /// Appends up to `letters` bases from packed NCBI4na bytes (`CSeqConvert::Convert` from
    /// `e_Ncbi4na` to `e_Ncbi8na`, clamped by `AdjustLength` to `src.len() * 2`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:350-356
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:808-837
    pub(crate) fn push_ncbi4na(&mut self, src: &[u8], letters: usize) {
        let count = letters.min(src.len().saturating_mul(2));
        if count == 0 {
            return;
        }
        if !self.odd {
            let full_bytes = count / 2;
            self.bytes.extend_from_slice(&src[..full_bytes]);
            if count % 2 == 1 {
                self.push_code(src[full_bytes] >> 4);
            }
            return;
        }
        for i in 0..count {
            let code = if i & 1 == 0 {
                src[i / 2] >> 4
            } else {
                src[i / 2] & 0x0f
            };
            self.push_code(code);
        }
    }

    /// The count of letters appended.
    pub(crate) fn letters(&self) -> usize {
        self.bytes.len() * 2 - usize::from(self.odd)
    }

    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

#[cfg(test)]
#[path = "pack_tests.rs"]
mod tests;
