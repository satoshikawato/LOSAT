//! `CSeqConvert::Pack` of IUPACna residues into the segments of a delta sequence
//! (`CDelta_ext::AddAndSplit`), and the conversions of the residues into NCBI2na,
//! NCBI4na and NCBI8na (`c++/src/util/sequtil/`).

use crate::utils::byte_search::find_first;

/// The codings that `Pack` chooses for IUPACna residues without gaps.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Coding {
    Ncbi2na,
    Ncbi4na,
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
/// The sums are made in 64 bits (NCBI's `TSeqPos` sum wraps for a length above
/// 4,294,967,292, which the reader never packs).
fn bytes_needed(coding: Coding, length: u32) -> u64 {
    match coding {
        Coding::Ncbi2na => u64::from(length).div_ceil(4),
        Coding::Ncbi4na => u64::from(length).div_ceil(2),
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

/// No coding before the first one of an arrangement.
const NO_NODE: u32 = u32::MAX;

/// `CPacker::SArrangement`: the codings of the segments so far and their cost. NCBI
/// keeps the codings in blocks of 16 that two arrangements share up to their last common
/// block, so that a copy costs one block; LOSAT keeps each coding as a node that points
/// to the one before it, so that a copy costs one index (the same list of codings).
#[derive(Clone, Copy)]
struct Arrangement {
    /// The last coding (`codings->GetLast()`), an index of `Packer::nodes`.
    last: u32,
    cost: u64,
}

/// `CSeqConvert_imp::CPacker` for IUPACna.
struct Packer {
    /// `m_Boundaries`.
    boundaries: Vec<u32>,
    /// The codings of both arrangements: each one with the index of the coding before it.
    nodes: Vec<(Coding, u32)>,
    /// `m_EndingNarrow` and `m_EndingWide`.
    ending_narrow: Arrangement,
    ending_wide: Arrangement,
}

impl Packer {
    /// `SArrangement::AddCoding`.
    fn add_coding(&mut self, last: u32, coding: Coding) -> u32 {
        self.nodes.push((coding, last));
        (self.nodes.len() - 1) as u32
    }

    /// `new_coding` is `None` for `kNoCoding` (the end of the residues).
    ///
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1194-1249
    fn add_boundary(&mut self, pos: u32, new_coding: Option<Coding>) {
        let Some(&last_boundary) = self.boundaries.last() else {
            let Some(new_coding) = new_coding else {
                return;
            };
            self.boundaries.push(pos);
            self.ending_narrow.last = self.add_coding(NO_NODE, new_coding);
            self.ending_wide.last = self.add_coding(NO_NODE, WIDE_CODING);
            self.ending_wide.cost = OVERHEAD;
            self.ending_narrow.cost = OVERHEAD;
            return;
        };
        let last_length = pos - last_boundary;
        self.boundaries.push(pos);

        let last_narrow = self.nodes[self.ending_narrow.last as usize].0;
        self.ending_narrow.cost += bytes_needed(last_narrow, last_length);
        self.ending_wide.cost += bytes_needed(WIDE_CODING, last_length);

        let Some(new_coding) = new_coding else {
            return;
        };

        if new_coding != WIDE_CODING && self.ending_narrow.cost > self.ending_wide.cost {
            self.ending_narrow = self.ending_wide;
        }

        let alt_wide_cost = self.ending_narrow.cost + OVERHEAD;
        self.ending_narrow.cost += OVERHEAD;

        if self.ending_wide.cost > alt_wide_cost {
            self.ending_wide = self.ending_narrow;
            self.ending_wide.cost = alt_wide_cost;
        } else if new_coding == WIDE_CODING {
            self.ending_narrow = self.ending_wide;
        }

        self.ending_narrow.last = self.add_coding(self.ending_narrow.last, new_coding);
        self.ending_wide.last = self.add_coding(self.ending_wide.last, WIDE_CODING);
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
        nodes: Vec::new(),
        ending_narrow: Arrangement {
            last: NO_NODE,
            cost: 0,
        },
        ending_wide: Arrangement {
            last: NO_NODE,
            cost: 0,
        },
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

    let best = if packer.ending_narrow.cost < packer.ending_wide.cost {
        packer.ending_narrow
    } else {
        packer.ending_wide
    };
    let n = packer.boundaries.len() - 1;
    let mut coding_vec = vec![Coding::Ncbi4na; n];
    let mut node = best.last;
    for slot in coding_vec.iter_mut().rev() {
        let (coding, previous) = packer.nodes[node as usize];
        *slot = coding;
        node = previous;
    }

    let mut segments = Vec::new();
    let mut i = 0;
    while i < n {
        let coding = coding_vec[i];
        let start = packer.boundaries[i];
        while i < n - 1 && coding_vec[i + 1] == coding {
            i += 1;
        }
        segments.push(Segment {
            coding,
            start,
            length: packer.boundaries[i + 1] - start,
        });
        i += 1;
    }
    segments
}

/// The NCBI2na code of the unambiguous IUPACna letters (`CIupacnaTo2na`'s rows `A`: 0,
/// `C`: 1, `G`: 2, `T` and `U`: 3; the table has the code shifted to each of the four
/// places of a byte). Only these letters reach an NCBI2na segment.
pub(crate) const fn ncbi2na_code(residue: u8) -> u8 {
    match residue {
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 0,
    }
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

const NCBI2NA_CODES: [u8; 256] = {
    let mut table = [0; 256];
    let mut index = 0;
    while index < 256 {
        table[index] = ncbi2na_code(index as u8);
        index += 1;
    }
    table
};

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
pub(crate) fn iupacna_to_ncbi2na(src: &[u8]) -> Vec<u8> {
    let mut dst = Vec::with_capacity(src.len().div_ceil(4));
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
    dst
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

/// Appends the NCBI8na codes of the residues of a delta literal, as `s_FixBioseqDeltas`
/// gets them from the literal's NCBI2na or NCBI4na data. Both ways give the NCBI4na code
/// of the letter: an NCBI2na literal holds only `A`, `C`, `G`, `T`, `U`, whose codes 0 to
/// 3 `C2naTo8na` turns into 1, 2, 4, 8, and `C4naTo8na` splits a byte into its two codes
/// (`tests::ncbi8na_of_both_codings`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:343-357
pub(super) fn append_ncbi8na(seq8na: &mut Vec<u8>, src: &[u8]) {
    seq8na.extend(src.iter().map(|&residue| NCBI4NA_CODES[residue as usize]));
}

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:856-874
/// The length of `src` can be one that the input states rather than holds (a gap line),
/// so the allocation of the result can fail without the process being short of memory
/// for anything else; that failure is returned (`fix_bioseq_deltas`).
pub(super) fn ncbi8na_to_ncbi4na(src: &[u8]) -> Result<Vec<u8>, std::collections::TryReserveError> {
    let mut dst = Vec::new();
    dst.try_reserve_exact(src.len().div_ceil(2))?;
    let mut chunks = src.chunks_exact(2);
    dst.extend((&mut chunks).map(|two| two[0] << 4 | two[1]));
    if let [last] = chunks.remainder() {
        dst.push(last << 4);
    }
    Ok(dst)
}

#[cfg(test)]
#[path = "pack_tests.rs"]
mod tests;
