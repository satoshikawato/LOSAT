//! `CSeqConvert::Pack` of IUPACna or NCBI4na residues into the segments of a delta sequence
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
    /// Inlined into each scan (`pack_runs` of letters, `pack_ncbi4na` of codes), as into the one
    /// scan before there were two (LOSAT's speed: a call for each boundary).
    #[inline(always)]
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
    pack_runs(src.len(), |pos| {
        let wide = BEST_IS_NCBI4NA[src[pos] as usize];
        let rest = &src[pos + 1..];
        let diff = |&residue: &u8| BEST_IS_NCBI4NA[residue as usize] != wide;
        let end = if rest.len() < 20 {
            pos + 1 + rest.iter().position(diff).unwrap_or(rest.len())
        } else if let Some(at) = rest[..4].iter().position(diff) {
            pos + 1 + at
        } else {
            let tail = &rest[4..];
            // The run of residues of the same best coding, found with one search (LOSAT's speed).
            let run = if wide == 0 {
                find_first(tail, |residue| !is_ncbi2na_best(residue))
            } else {
                find_first(tail, is_ncbi2na_best)
            };
            pos + 5 + run.unwrap_or(tail.len())
        };
        (wide == 1, end)
    })
}

/// Whether the NCBI4na code `code` (1 to 15) is one base (1, 2, 4 or 8): the class of a
/// residue in `kBestCodingsWithoutGaps.ncbi4na`, NCBI2na for one base and NCBI4na for the other
/// codes. A letter's code is one base exactly where the letter's best IUPACna coding is NCBI2na
/// (`tests::the_class_of_a_code_is_the_class_of_its_letter`), so the residues of a letter and of
/// its code pack alike. Code 0 (a gap) is never a residue.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_best_codings.inc:34-41,74-98
#[inline(always)]
const fn is_one_base(code: u8) -> bool {
    code & code.wrapping_sub(1) == 0
}

/// Whether a byte of two NCBI4na codes (1 to 15) holds a code that is not one base, tested by
/// arithmetic so that `find_first` tests a block of bytes at once.
#[inline(always)]
const fn holds_ambiguous(byte: u8) -> bool {
    let (high, low) = (byte >> 4, byte & 0x0f);
    (high & high.wrapping_sub(1)) | (low & low.wrapping_sub(1)) != 0
}

/// Whether a byte of two NCBI4na codes (1 to 15) holds a code of one base ([`holds_ambiguous`]'s
/// arithmetic).
#[inline(always)]
const fn holds_one_base(byte: u8) -> bool {
    let (high, low) = (byte >> 4, byte & 0x0f);
    (high & high.wrapping_sub(1) == 0) | (low & low.wrapping_sub(1) == 0)
}

/// The NCBI4na code of residue `index` of `src` (two to a byte, the first in the high half).
#[inline(always)]
fn nibble(src: &[u8], index: usize) -> u8 {
    (src[index / 2] >> (4 - 4 * (index % 2))) & 0x0f
}

/// The classes of the two NCBI4na codes of a byte (the first in the high half): bit 1 is set
/// where the first code is not one base ([`is_one_base`]), bit 0 where the second is not. A
/// byte inside a run of one class has the bits 0 (codes of one base) or 3 (the other codes).
const PAIR_CLASSES: [u8; 256] = {
    let mut table = [0; 256];
    let mut byte = 0;
    while byte < 256 {
        table[byte] =
            (!is_one_base((byte >> 4) as u8) as u8) << 1 | !is_one_base((byte & 0x0f) as u8) as u8;
        byte += 1;
    }
    table
};

/// `CSeqConvert::Pack` of the NCBI4na residues `start..end` of `src` (two to a byte, the first
/// in the high half; `end - start` at most 4,294,967,295) as `e_Ncbi4na` into a
/// `CDelta_ext_PackTarget`: the segments that `CDelta_ext::AddAndSplit` appends for them, with
/// their starts counted from `start`. A residue's class is whether its code is one base
/// ([`is_one_base`]), which is the class of its letter, so the segments are those that
/// [`pack_iupacna`] gives for the residues' letters (makeblastdb's reader stores codes,
/// `reader::Ncbi4naResidues`; NCBI's packs the letters). The bytes of a segment are the caller's
/// (`ncbi4na_to_ncbi2na`, or the residues' own bytes).
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1058-1063,1078-1080,1105-1107
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1124-1192
/// NCBI's scan reads a byte of two codes at a time and adds the boundary of each half where the
/// two differ; so does LOSAT's, from `start` (which may be the low half of a byte), with the
/// classes of both codes of a byte from [`PAIR_CLASSES`]. While twelve or more bytes are left,
/// the bytes after two bytes in a row of the run's class are searched for the next byte with a
/// code of the other class ([`find_first`], LOSAT's speed: four letters, then a search when
/// twenty or more are left, as `pack_iupacna` scans letters).
pub(crate) fn pack_ncbi4na(src: &[u8], start: usize, end: usize) -> Vec<Segment> {
    // The class bits of a byte inside the current run (`inside`) are 3 (codes that are not one
    // base) or 0; `flip` starts a run of the other class at residue `pos`.
    #[inline(always)]
    fn flip(packer: &mut Packer, start: usize, pos: usize, inside: &mut u8) {
        *inside ^= 3;
        packer.add_boundary((pos - start) as u32, Some(coding_of_class(*inside == 3)));
    }
    // The boundaries in the whole byte `j` whose class bits differ from the run's by `differs`:
    // 1, the low half; 2, the high half alone (the low half is of the class before it again);
    // 3, both, a new class from the high half on. Inlined into both loops below (LOSAT's speed:
    // a call for each byte that differs).
    #[inline(always)]
    fn flip_byte(packer: &mut Packer, start: usize, j: usize, inside: &mut u8, differs: u8) {
        if differs & 2 != 0 {
            flip(packer, start, 2 * j, inside);
        }
        if differs != 3 {
            flip(packer, start, 2 * j + 1, inside);
        }
    }
    if start == end {
        return Vec::new();
    }
    let mut packer = Packer::new();
    let mut inside = if is_one_base(nibble(src, start)) {
        0
    } else {
        3
    };
    packer.add_boundary(0, Some(coding_of_class(inside == 3)));
    // The residue after `start`: the low half of its byte when `start` is even.
    let mut j = start / 2;
    if start % 2 == 0 && start + 1 < end && (PAIR_CLASSES[src[j] as usize] ^ inside) & 1 != 0 {
        flip(&mut packer, start, start + 1, &mut inside);
    }
    j += 1;
    // The whole bytes up to `end`, in each the high half, then the low half. While twelve or
    // more are left, two bytes of the run's class in a row are followed by a search for the next
    // byte with a code of the other class.
    let last = end / 2;
    let mut quiet = false;
    while j + 12 <= last {
        let differs = PAIR_CLASSES[src[j] as usize] ^ inside;
        j += 1;
        if differs != 0 {
            quiet = false;
            flip_byte(&mut packer, start, j - 1, &mut inside, differs);
        } else if quiet {
            let bytes = &src[j..last];
            let found = if inside == 3 {
                find_first(bytes, holds_one_base)
            } else {
                find_first(bytes, holds_ambiguous)
            };
            j += found.unwrap_or(bytes.len());
            quiet = false;
        } else {
            quiet = true;
        }
    }
    // The bytes left, one at a time.
    while j < last {
        let differs = PAIR_CLASSES[src[j] as usize] ^ inside;
        if differs != 0 {
            flip_byte(&mut packer, start, j, &mut inside, differs);
        }
        j += 1;
    }
    // With `end` odd, the high half of the byte `last` is the last residue.
    if end % 2 == 1 && j == last && (PAIR_CLASSES[src[last] as usize] ^ inside) & 2 != 0 {
        flip(&mut packer, start, 2 * last, &mut inside);
    }
    packer.add_boundary((end - start) as u32, None);
    packer.segments()
}

/// `CPacker::Pack` once its scan has found the runs of one best coding: `next_run(pos)` gives
/// whether the residue at `pos` is best as NCBI4na and where its run ends, for the residues
/// `0..length`. The segments of the best arrangement, each of the codings of its runs.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1124-1192
fn pack_runs(length: usize, mut next_run: impl FnMut(usize) -> (bool, usize)) -> Vec<Segment> {
    if length == 0 {
        return Vec::new();
    }
    let mut packer = Packer::new();
    let mut pos = 0;
    while pos < length {
        let (wide, end) = next_run(pos);
        packer.add_boundary(pos as u32, Some(coding_of_class(wide)));
        pos = end;
    }
    packer.add_boundary(length as u32, None);
    packer.segments()
}

/// The coding of a residue of the class `wide` (not one base: NCBI4na; one base: NCBI2na).
#[inline(always)]
fn coding_of_class(wide: bool) -> Coding {
    if wide {
        Coding::Ncbi4na
    } else {
        Coding::Ncbi2na
    }
}

impl Packer {
    fn new() -> Self {
        Self {
            boundaries: Vec::new(),
            steps: Vec::new(),
            narrow_cost: 0,
            wide_cost: 0,
            last_narrow: WIDE_CODING,
        }
    }

    /// The segments of the best arrangement, once the last boundary (the end of the residues)
    /// is added: each of the codings of its runs.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1124-1192
    /// Inlined into each scan, as `add_boundary` (LOSAT's speed: a short record's call and copy
    /// of the packer cost more than its steps).
    #[inline(always)]
    fn segments(mut self) -> Vec<Segment> {
        // The best arrangement's codings, from the last segment back: each step becomes 1
        // where the segment it starts is NCBI4na.
        let mut narrow = self.narrow_cost < self.wide_cost;
        for step in self.steps.iter_mut().rev() {
            let bits = *step;
            let wide = !narrow || bits & STEP_NEW_WIDE != 0;
            narrow = if narrow {
                bits & STEP_NARROW_FROM_WIDE == 0
            } else {
                bits & STEP_WIDE_FROM_NARROW != 0 && bits & STEP_NARROW_FROM_WIDE == 0
            };
            *step = u8::from(wide);
        }
        let steps = &self.steps;
        let n = steps.len();
        let mut segments = Vec::new();
        let mut i = 0;
        while i < n {
            let wide = steps[i];
            let start = self.boundaries[i];
            while i < n - 1 && steps[i + 1] == wide {
                i += 1;
            }
            segments.push(Segment {
                coding: coding_of_class(wide == 1),
                start,
                length: self.boundaries[i + 1] - start,
            });
            i += 1;
        }
        segments
    }
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

/// The number of letters that [`Ncbi4naBuilder::push_letters`] converts at a time
/// ([`ncbi4na_block`], LOSAT's speed).
const LETTER_BLOCK: usize = 64;

/// The NCBI4na code of `A`, `C`, `G`, `T`, `U` and `N` in either case (1, 2, 4, 8, 8, 15, their
/// rows of `NCBI4NA_CODES`); 0 for every other byte. No table and no branch, so that a loop over
/// letters computes the codes of many at once ([`ncbi4na_block`]). Clearing bit 5 folds lower case
/// to upper case, and clearing bit 0 pairs `T` (`0x54`) with `U` (`0x55`), as [`is_ncbi2na_best`]
/// does (`tests::the_codes_of_common_letters_are_the_table` checks every byte).
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:392-394,409-410
#[inline(always)]
const fn common_letter_code(letter: u8) -> u8 {
    let upper = letter & 0xdf;
    (upper == b'A') as u8
        | ((upper == b'C') as u8) << 1
        | ((upper == b'G') as u8) << 2
        | (((upper & 0xfe) == b'T') as u8) << 3
        | (((upper == b'N') as u8) * 15)
}

/// The byte of two NCBI4na codes (each 0 to 15), the first in the high half: the low byte of the
/// pair as a little-endian 16-bit number shifted left by 4, ORed with its high byte (one shift
/// and one OR, which a loop over pairs computes for many at once).
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:856-874
#[inline(always)]
pub(crate) fn pair_byte(two: [u8; 2]) -> u8 {
    let pair = u16::from_le_bytes(two);
    ((pair << 4) | (pair >> 8)) as u8
}

/// The NCBI4na bytes of [`LETTER_BLOCK`] residue letters (`x_ConvertIupacnaTo4na`'s loop): the
/// codes of a block whose letters are all `A`, `C`, `G`, `T`, `U` or `N` by comparisons
/// ([`common_letter_code`]) and then two to a byte ([`pair_byte`]), without a table and without a
/// branch a letter; a block with another letter takes the table for each pair, as
/// [`iupacna_to_ncbi4na`] does.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
#[inline(always)]
fn ncbi4na_block(block: &[u8; LETTER_BLOCK]) -> [u8; LETTER_BLOCK / 2] {
    let mut codes = [0u8; LETTER_BLOCK];
    for (code, &letter) in codes.iter_mut().zip(block) {
        *code = common_letter_code(letter);
    }
    let mut bytes = [0u8; LETTER_BLOCK / 2];
    // A fold of bytes, not of `bool`s: rustc 1.92 vectorizes only the former (`byte_search`).
    if codes
        .iter()
        .fold(0u8, |any, &code| any | u8::from(code == 0))
        == 0
    {
        for (byte, &two) in bytes.iter_mut().zip(codes.as_chunks::<2>().0) {
            *byte = pair_byte(two);
        }
    } else {
        for (byte, two) in bytes.iter_mut().zip(block.as_chunks::<2>().0) {
            *byte = NCBI4NA_CODES[two[0] as usize] << 4 | NCBI4NA_CODES[two[1] as usize];
        }
    }
    bytes
}

/// The NCBI2na codes of the four NCBI4na codes of `pair` (two bytes, the first code in the high
/// half of the first byte), each code one base (1, 2, 4, 8 for `A`, `C`, `G`, `T`: 0, 1, 2, 3) or
/// 0 (an NCBI2na code 0): four to a byte, the first in the high bits. A code's high bit is
/// whether it is 4 or 8, its low bit whether it is 2 or 8. The two bytes are one 16-bit number
/// and every step is arithmetic on it (no byte of the pair taken apart), so that a loop over pairs
/// computes many at once.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1190-1268
#[inline(always)]
const fn ncbi2na_of_pair(pair: [u8; 2]) -> u8 {
    let word = u16::from_le_bytes(pair);
    // In each half byte, bit 2: bit 2 or 3 of the code; bit 1: bit 1 or 3. Then down by one: the
    // first byte's codes in bits 4-5 and 0-1, the second byte's in bits 12-13 and 8-9.
    let bits = (((word | word >> 1) & 0x4444) | ((word | word >> 2) & 0x2222)) >> 1;
    let first = ((bits >> 2) & 0x0c) | (bits & 0x03);
    let second = ((bits >> 10) & 0x0c) | ((bits >> 8) & 0x03);
    (first << 4 | second) as u8
}

/// The NCBI2na codes of the sixteen NCBI4na codes of eight bytes, `word` (the bytes as a
/// little-endian number; the first code in the high half of the first byte), each code one base
/// (1, 2, 4, 8 for `A`, `C`, `G`, `T`: 0, 1, 2, 3) or 0 (an NCBI2na code 0): four to a byte, the
/// first in the high bits, as the four bytes of the little-endian result. A code's high bit is
/// whether it is 4 or 8, its low bit whether it is 2 or 8. Every step is arithmetic on the one
/// 64-bit number (no byte taken apart), a few operations for eight bytes.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1190-1268
#[inline(always)]
const fn ncbi2na_of_eight(word: u64) -> u32 {
    // In each half byte, bit 2: bit 2 or 3 of the code; bit 1: bit 1 or 3. Then down by one: each
    // half byte holds its NCBI2na code in its bits 0-1.
    let codes = (((word | word >> 1) & 0x4444_4444_4444_4444)
        | ((word | word >> 2) & 0x2222_2222_2222_2222))
        >> 1;
    // Each byte: its first code beside its second, in bits 0-3.
    let fours = (codes | codes >> 2) & 0x0f0f_0f0f_0f0f_0f0f;
    // Each 16 bits: the first byte's four bits above the second's, an output byte in bits 0-7.
    let bytes = ((fours << 4) | (fours >> 8)) & 0x00ff_00ff_00ff_00ff;
    // The four output bytes side by side.
    let halves = (bytes | bytes >> 8) & 0x0000_ffff_0000_ffff;
    (halves | halves >> 16) as u32
}

/// `CSeqConvert::Convert` from `e_Ncbi4na` to `e_Ncbi2na` of the first `letters` residues of
/// `src` (`x_Convert4naTo2na` at an even position), every one a base (`pack_ncbi4na` gave them
/// the NCBI2na coding): four to a byte, the first in the high bits, and 0 bits after the last
/// (the low half of `src`'s last byte is 0 for an odd count, and the bytes after the last are
/// taken as 0). These are the bytes that [`iupacna_to_ncbi2na`] makes of their letters.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:663-705
pub(crate) fn ncbi4na_to_ncbi2na(src: &[u8], letters: usize) -> Vec<u8> {
    debug_assert_eq!(src.len(), letters.div_ceil(2));
    let mut dst = Vec::with_capacity(letters.div_ceil(4));
    let (blocks, rest) = src.split_at(src.len() - src.len() % NCBI4NA_BLOCK);
    dst.extend(
        blocks
            .as_chunks::<2>()
            .0
            .iter()
            .map(|&pair| ncbi2na_of_pair(pair)),
    );
    let (eights, last) = rest.as_chunks::<8>();
    for eight in eights {
        dst.extend_from_slice(&ncbi2na_of_eight(u64::from_le_bytes(*eight)).to_le_bytes());
    }
    if !last.is_empty() {
        // The last bytes with 0 bytes after them, made byte by byte (not a call to copy them).
        let word = last
            .iter()
            .rev()
            .fold(0u64, |word, &byte| word << 8 | u64::from(byte));
        dst.extend_from_slice(&ncbi2na_of_eight(word).to_le_bytes()[..last.len().div_ceil(2)]);
    }
    dst
}

/// [`ncbi4na_to_ncbi2na`] converts the bytes of whole blocks of this many a pair at a time
/// ([`ncbi2na_of_pair`], which the loop computes as vector code), and the bytes after them eight at
/// a time ([`ncbi2na_of_eight`]), so that no byte takes the scalar loop of pairs (LOSAT's speed).
const NCBI4NA_BLOCK: usize = 64;

/// NCBI4na bytes made as letters or codes come, two letters to a byte, the first in the high
/// half; with an odd count the low half of the last byte is 0. Two owners use it:
/// - makeblastdb's nucleotide reader keeps a record's residues in one (`reader::Ncbi4naResidues`,
///   [`Self::push_letters`]), which becomes the record's NCBI4na data where that is its data
///   (`reader::SeqInst`);
/// - `s_FixBioseqDeltas`' `seq8na` and its NCBI4na copy in one (`fix_bioseq_deltas`, and
///   `asn_source` for the ASN.1 input): the NCBI4na bytes of the NCBI8na codes appended
///   (`C8naTo4na`). NCBI makes `seq8na` (1 byte a letter) and converts it; LOSAT writes the
///   NCBI4na bytes (half a byte a letter) as the codes come, the same bytes without the `seq8na`.
///
/// The NCBI8na codes of a literal's residues are the NCBI4na codes of its letters, as
/// `s_FixBioseqDeltas` gets them from the literal's NCBI2na or NCBI4na data: an NCBI2na
/// literal holds only `A`, `C`, `G`, `T`, `U`, whose codes 0 to 3 `C2naTo8na` turns into 1,
/// 2, 4, 8, and `C4naTo8na` splits a byte into its two codes
/// (`tests::ncbi8na_of_both_codings`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:343-357
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:856-874
#[derive(Debug, Default, PartialEq, Eq)]
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

    /// The NCBI4na codes of the residue letters `src` (`x_ConvertIupacnaTo4na`, the codes of
    /// `NCBI4NA_CODES`): after a pending letter's partner, whole blocks of [`LETTER_BLOCK`]
    /// letters ([`ncbi4na_block`]), then the letters after the last block two at a time through
    /// the table (a reader converts whole blocks until the end of a record,
    /// `reader::Ncbi4naResidues`). The bytes of all of `src` are reserved first, so that a
    /// reader's store is first as long as the letters it converts at once, and then grows as a
    /// `Vec` does.
    ///
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
    pub(crate) fn push_letters(&mut self, mut src: &[u8]) {
        self.bytes.reserve(src.len().div_ceil(2));
        if let (true, [first, rest @ ..]) = (self.odd, src) {
            self.push_code(NCBI4NA_CODES[*first as usize]);
            src = rest;
        }
        while let Some((block, rest)) = src.split_first_chunk::<LETTER_BLOCK>() {
            self.bytes.extend_from_slice(&ncbi4na_block(block));
            src = rest;
        }
        let (pairs, rest) = src.as_chunks::<2>();
        self.bytes.extend(
            pairs
                .iter()
                .map(|two| NCBI4NA_CODES[two[0] as usize] << 4 | NCBI4NA_CODES[two[1] as usize]),
        );
        if let [last] = rest {
            self.push_code(NCBI4NA_CODES[*last as usize]);
        }
    }

    /// The NCBI4na codes of the letters `src`, one at a time through `NCBI4NA_CODES`: the
    /// conversion of [`Self::push_letters`] without its blocks, for the tests' comparisons.
    #[cfg(test)]
    pub(crate) fn push_letters_by_table(&mut self, src: &[u8]) {
        for &letter in src {
            self.push_code(NCBI4NA_CODES[letter as usize]);
        }
    }

    /// Appends up to `letters` bases from packed NCBI2na bytes (`CSeqConvert::Convert` from
    /// `e_Ncbi2na` to `e_Ncbi8na`, clamped by `AdjustLength` to `src.len() * 4`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:343-349
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:521-531
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

    /// Appends up to `letters` bases of packed NCBI4na bytes from the base `start` of `src`
    /// (`CSeqConvert::Convert` from `e_Ncbi4na` to `e_Ncbi8na` at that position, clamped by
    /// `AdjustLength` to the bases of `src`): a copy of whole bytes where `start` and the
    /// letters appended so far are both even or both odd (after the one base that fills a
    /// pending half), else each byte from the halves of two (a shift by one base).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:350-356
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:761-767
    #[inline(always)]
    pub(crate) fn push_ncbi4na(&mut self, src: &[u8], start: usize, letters: usize) {
        let count = letters.min(src.len().saturating_mul(2).saturating_sub(start));
        if count == 0 {
            return;
        }
        let (mut from, end) = (start, start + count);
        if self.odd {
            self.push_code(nibble(src, from));
            from += 1;
        }
        let pairs = (end - from) / 2;
        if pairs > 0 {
            self.push_ncbi4na_pairs(src, from, pairs);
            from += 2 * pairs;
        }
        if from < end {
            self.push_code(nibble(src, from));
        }
    }

    /// The bytes of the `pairs` pairs of bases from the base `from` of `src`, no half pending
    /// ([`Self::push_ncbi4na`]). A function of its own, so that the loop over a delta's pieces,
    /// into which `push_ncbi4na` is inlined, keeps its registers for literals of a base or two
    /// (LOSAT's speed).
    #[inline(never)]
    fn push_ncbi4na_pairs(&mut self, src: &[u8], from: usize, pairs: usize) {
        debug_assert!(!self.odd);
        let first = from / 2;
        if from.is_multiple_of(2) {
            self.bytes.extend_from_slice(&src[first..first + pairs]);
        } else {
            self.bytes.extend(
                src[first..first + pairs + 1]
                    .windows(2)
                    .map(|two| two[0] << 4 | two[1] >> 4),
            );
        }
    }

    /// The count of letters appended.
    pub(crate) fn letters(&self) -> usize {
        self.bytes.len() * 2 - usize::from(self.odd)
    }

    /// The bytes (`letters()` of them, two to a byte).
    pub(crate) fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// No letters, the bytes' capacity kept.
    pub(crate) fn clear(&mut self) {
        self.bytes.clear();
        self.odd = false;
    }

    /// Drops the first `count` letters (the reader's residues past 2^32,
    /// `reader::Ncbi4naResidues`): the bytes after them move to the front, each made of the
    /// halves of two where `count` is odd.
    #[cold]
    pub(crate) fn drop_front(&mut self, count: usize) {
        let kept = self.letters() - count;
        let first = count / 2;
        if count.is_multiple_of(2) {
            self.bytes.drain(..first);
        } else {
            let bytes = kept.div_ceil(2);
            for j in 0..bytes {
                let next = self.bytes.get(first + j + 1).map_or(0, |&byte| byte >> 4);
                self.bytes[j] = self.bytes[first + j] << 4 | next;
            }
            self.bytes.truncate(bytes);
        }
        self.odd = kept % 2 == 1;
    }

    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

#[cfg(test)]
#[path = "pack_tests.rs"]
mod tests;
