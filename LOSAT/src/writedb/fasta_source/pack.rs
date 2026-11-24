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

/// Whether the best coding of `residue` is NCBI2na ([`best_coding_iupacna`]), by arithmetic on
/// the byte, which the scan of `pack_iupacna` tests a block of residues at a time (the unit
/// test checks every byte). Clearing bit 5 makes a lower case letter upper case and moves no
/// other byte onto an upper case letter.
const fn is_ncbi2na_best(residue: u8) -> bool {
    let upper = residue & 0xdf;
    (upper == b'A') | (upper == b'C') | (upper == b'G') | (upper == b'T') | (upper == b'U')
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
        let coding = best_coding_iupacna(src[pos]);
        packer.add_boundary(pos as u32, Some(coding));
        let narrow = coding == Coding::Ncbi2na;
        // The run of residues of the same best coding, found with one search (LOSAT's speed).
        let rest = &src[pos + 1..];
        let run = if narrow {
            find_first(rest, |residue| !is_ncbi2na_best(residue))
        } else {
            find_first(rest, is_ncbi2na_best)
        };
        pos += 1 + run.unwrap_or(rest.len());
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

/// The NCBI2na code of the unambiguous IUPACna letters (`CIupacnaTo2na`'s rows `A`: 0,
/// `C`: 1, `G`: 2, `T` and `U`: 3, every other row 0; the table has the code shifted to each
/// of the four places of a byte). Only these letters reach an NCBI2na segment. Written
/// without a branch, so that a loop over residues computes the codes of many at once.
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

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:257-289
/// The codes of a block of residues are computed first (one loop without a branch) and then
/// packed four to a byte by one product each ([`pack_four_codes`]); the block on the stack keeps
/// the memory as it was (LOSAT's speed).
pub(crate) fn iupacna_to_ncbi2na(src: &[u8]) -> Vec<u8> {
    const BLOCK: usize = 256;
    let mut dst = Vec::with_capacity(src.len().div_ceil(4));
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
    dst
}

/// The byte of four NCBI2na codes (each 0 to 3), the first in the high bits:
/// `c0 << 6 | c1 << 4 | c2 << 2 | c3`. With the codes as the bytes of a little-endian word `w`
/// (code `k` at bit `8 k`), the product `w * (2^30 + 2^20 + 2^10 + 1)` puts code `k` at bit
/// `30 - 2 k` (bits 24 to 31); the other terms put code `k` at bit `8 k + 30 - 10 j` for `j != k`,
/// below 22 or above 33, where two bits each never carry into bits 24 to 31.
fn pack_four_codes(four: [u8; 4]) -> u8 {
    ((u64::from(u32::from_le_bytes(four)) * 0x4010_0401) >> 24) as u8
}

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
pub(crate) fn iupacna_to_ncbi4na(src: &[u8]) -> Vec<u8> {
    let mut dst = Vec::with_capacity(src.len().div_ceil(2));
    let mut chunks = src.chunks_exact(2);
    dst.extend(
        (&mut chunks)
            .map(|two| NCBI4NA_CODES[two[0] as usize] << 4 | NCBI4NA_CODES[two[1] as usize]),
    );
    if let [last] = chunks.remainder() {
        dst.push(NCBI4NA_CODES[*last as usize] << 4);
    }
    dst
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
pub(super) struct Ncbi4naBuilder {
    bytes: Vec<u8>,
    /// The last byte holds one letter, in its high half.
    odd: bool,
}

impl Ncbi4naBuilder {
    /// The builder of `letters` letters, whose bytes are reserved at once. The count can be
    /// one that the input states rather than holds (a gap line), so the allocation can fail
    /// without the process being short of memory for anything else; that failure is
    /// returned (`fix_bioseq_deltas`), and no later growth can fail.
    pub(super) fn with_letters(letters: usize) -> Result<Self, std::collections::TryReserveError> {
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
    pub(super) fn push_run(&mut self, code: u8, mut run: usize) {
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
    pub(super) fn push_iupacna(&mut self, mut src: &[u8]) {
        if let (true, [first, rest @ ..]) = (self.odd, src) {
            self.push_code(NCBI4NA_CODES[*first as usize]);
            src = rest;
        }
        let mut chunks = src.chunks_exact(2);
        self.bytes.extend(
            (&mut chunks)
                .map(|two| NCBI4NA_CODES[two[0] as usize] << 4 | NCBI4NA_CODES[two[1] as usize]),
        );
        if let [last] = chunks.remainder() {
            self.push_code(NCBI4NA_CODES[*last as usize]);
        }
    }

    /// The count of letters appended.
    pub(super) fn letters(&self) -> usize {
        self.bytes.len() * 2 - usize::from(self.odd)
    }

    pub(super) fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

#[cfg(test)]
#[path = "pack_tests.rs"]
mod tests;
