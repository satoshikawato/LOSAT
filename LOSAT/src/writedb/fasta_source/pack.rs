//! `CSeqConvert::Pack` of IUPACna residues into the segments of a delta sequence
//! (`CDelta_ext::AddAndSplit`), and the conversions of the residues into NCBI2na,
//! NCBI4na and NCBI8na (`c++/src/util/sequtil/`).

/// The codings that `Pack` chooses for IUPACna residues without gaps.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Coding {
    Ncbi2na,
    Ncbi4na,
}

/// One segment of the packed residues (`IPackTarget::NewSegment`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) struct Segment {
    pub(super) coding: Coding,
    pub(super) start: u32,
    pub(super) length: u32,
}

/// `kBestCodingsWithoutGaps.iupacna`: NCBI2na for `A`, `C`, `G`, `T`, `U` in either
/// case, NCBI4na for every other byte.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_best_codings.inc:34-62
/// ```c++
/// #define N2 CSeqUtil::e_Ncbi2na
/// #define N4 CSeqUtil::e_Ncbi4na
/// #define NX CSeqUtil::e_Ncbi4na_expand // mix of the above
///
/// #ifdef COUNT_GAPS
/// #  define gp CSeqUtil::e_not_set
/// #else
/// #  define gp N4
/// #endif
/// ...
/// #define N_ gp
/// //////////////// iupacna
/// {
///     ...
///     N_, N2, N4, N2, N4, N4, N4, N2, N4, N4, N4, N4, N4, N4, gp, N4, // ACGN
///     N4, N4, N4, N4, N2, N2, N4, N4, N4, N4, N4, N_, N_, N_, N_, N_, // TU
///     N_, N2, N4, N2, N4, N4, N4, N2, N4, N4, N4, N4, N4, N4, gp, N4, // acgn
///     N4, N4, N4, N4, N2, N2, N4, N4, N4, N4, N4, N_, N_, N_, N_, N_, // tu
/// ```
/// The table without gaps is the one used: `CDelta_ext_PackTarget::GapsOK` is the
/// `gaps_ok` of `AddAndSplit`, which `CFastaReader::AssembleSeq` passes as false.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1064-1071
/// ```c++
///     bool gaps_ok = dst.GapsOK(CSeqUtil::GetCodingType(src_coding));
///     const SBestCodings& best_codings = (gaps_ok ? kBestCodingsWithGaps
///                                         : kBestCodingsWithoutGaps);
///     const TCoding* best_coding = NULL;
///     switch (src_coding) {
///     case CSeqUtil::e_Iupacna:
///         best_coding = best_codings.iupacna;
///         break;
/// ```
const fn best_coding_iupacna(residue: u8) -> Coding {
    match residue {
        b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u' => Coding::Ncbi2na,
        _ => Coding::Ncbi4na,
    }
}

/// 1 for the bytes whose best coding is NCBI4na (a table for the scan of `pack_iupacna`).
const BEST_IS_NCBI4NA: [u8; 256] = {
    let mut table = [0; 256];
    let mut index = 0;
    while index < 256 {
        table[index] = matches!(best_coding_iupacna(index as u8), Coding::Ncbi4na) as u8;
        index += 1;
    }
    table
};

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_shared.cpp:195-203
/// ```c++
/// SIZE_TYPE GetBytesNeeded(CSeqUtil::TCoding coding, TSeqPos length)
/// {
///     switch (coding) {
///     case CSeqUtil::e_not_set:  return 0;
///     case CSeqUtil::e_Ncbi2na:  return (length + 3) / 4;
///     case CSeqUtil::e_Ncbi4na:  return (length + 1) / 2;
///     default:                   return length;
///     }
/// }
/// ```
/// The sums are made in 64 bits (NCBI's `TSeqPos` sum wraps for a length above
/// 4,294,967,292, which the reader never packs).
fn bytes_needed(coding: Coding, length: u32) -> u64 {
    match coding {
        Coding::Ncbi2na => (u64::from(length) + 3) / 4,
        Coding::Ncbi4na => (u64::from(length) + 1) / 2,
    }
}

/// `CDelta_ext_PackTarget::GetOverhead` of a segment with data.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Delta_ext.cpp:126-129
/// ```c++
///     // Rough memory overhead per segment, in bytes.  The actual values
///     // are system specific, but this code's behavior shouldn't be.
///     SIZE_TYPE GetOverhead(TCoding coding) const
///         { return coding == CSeqUtil::e_not_set ? 96 : 128; }
/// ```
const OVERHEAD: u64 = 128;

/// `CPacker::x_GetWideCoding(e_Iupacna)`.
///
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1292-1296
/// ```c++
///     switch (coding) {
///     case CSeqUtil::e_Iupacna:
///     case CSeqUtil::e_Ncbi4na_expand:
///     case CSeqUtil::e_Ncbi8na: // ignore values >= 16 for now
///         return CSeqUtil::e_Ncbi4na;
/// ```
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
    /// ```c++
    /// void CSeqConvert_imp::CPacker::x_AddBoundary(TSeqPos pos, TCoding new_coding)
    /// {
    ///     if (m_Boundaries.empty()) {
    ///         _ASSERT(pos == 0);
    ///         m_Boundaries.push_back(pos);
    ///         m_EndingNarrow.AddCoding(new_coding);
    ///         m_EndingWide.AddCoding(m_WideCoding);
    ///         m_EndingWide.cost   = m_Target.GetOverhead(m_WideCoding);
    ///         m_EndingNarrow.cost = m_Target.GetOverhead(new_coding);
    ///         return;
    ///     }
    ///
    ///     TSeqPos last_length = pos - m_Boundaries.back();
    ///     _ASSERT(last_length > 0);
    ///     m_Boundaries.push_back(pos);
    ///
    ///     TCoding last_narrow = m_EndingNarrow.codings->GetLast();
    ///     // This always rounds up to full bytes, and as such can slightly
    ///     // overestimate the total cost of ncbi4na.
    ///     m_EndingNarrow.cost += GetBytesNeeded(last_narrow,  last_length);
    ///     m_EndingWide.cost   += GetBytesNeeded(m_WideCoding, last_length);
    ///     if (last_narrow == m_WideCoding) {
    ///         _ASSERT(m_EndingNarrow.cost == m_EndingWide.cost);
    ///     }
    ///
    ///     // _TRACE(static_cast<char>('@' + last_narrow) << last_length << ": ("
    ///     //        << m_EndingNarrow.cost << ", " << m_EndingWide.cost << ')');
    ///
    ///     _ASSERT(new_coding != last_narrow);
    ///     if (new_coding == kNoCoding) {
    ///         return;
    ///     }
    ///
    ///     if (new_coding != m_WideCoding
    ///         &&  m_EndingNarrow.cost > m_EndingWide.cost) {
    ///         m_EndingNarrow = m_EndingWide;
    ///         // _TRACE("N <- W before N'");
    ///     }
    ///
    ///     SIZE_TYPE alt_wide_cost
    ///         = m_EndingNarrow.cost + m_Target.GetOverhead(m_WideCoding);
    ///     m_EndingNarrow.cost += m_Target.GetOverhead(new_coding);
    ///
    ///     if (m_EndingWide.cost > alt_wide_cost) {
    ///         m_EndingWide = m_EndingNarrow;
    ///         m_EndingWide.cost = alt_wide_cost;
    ///         // _TRACE("W <- N");
    ///     } else if (new_coding == m_WideCoding) {
    ///         m_EndingNarrow = m_EndingWide;
    ///         // _TRACE("N <- W before W");
    ///     }
    ///     // _TRACE('(' << m_EndingNarrow.cost << ", " << m_EndingWide.cost << ')');
    ///
    ///     m_EndingNarrow.AddCoding(new_coding);
    ///     m_EndingWide.AddCoding(m_WideCoding);
    /// }
    /// ```
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

/// `CSeqConvert::Pack(src, length, e_Iupacna, dst)` into a `CDelta_ext_PackTarget`: the
/// segments that `CDelta_ext::AddAndSplit` appends for the residues `src` (at most
/// 4,294,967,295 of them).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Delta_ext.cpp:201-205
/// ```c++
///     if (allow_packing)
///     {
///         CDelta_ext_PackTarget dst(*this, gaps_ok);
///         CSeqConvert::Pack(src.data(), length, coding, dst);
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1058-1063,1105-1107
/// ```c++
/// SIZE_TYPE CSeqConvert_imp::Pack(const char* src, TSeqPos length,
///                                 TCoding src_coding, IPackTarget& dst)
/// {
///     if (length == 0) {
///         return 0;
///     }
///     if (best_coding != NULL) {
///         return CPacker(src_coding, best_coding, gaps_ok, dst)
///             .Pack(src, length);
/// ```
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1124-1192
/// ```c++
/// SIZE_TYPE CSeqConvert_imp::CPacker::Pack(const char* src, TSeqPos length)
/// {
///     const char* src_end   = src + GetBytesNeeded(m_SrcCoding, length);
///     TCoding     prev_type = kNoCoding;
///
///     for (const char* p = src;  p < src_end;  ++p) {
///         unsigned char residue;
///         TCoding       curr_type;
///         do {
///             residue = static_cast<unsigned char>(*p);
///             curr_type = m_BestCoding[residue];
///         } while (curr_type == prev_type  &&  ++p < src_end);
///
///         if (curr_type == CSeqUtil::e_Ncbi4na_expand) {
///             ...
///         } else if (p != src_end) {
///             _ASSERT(curr_type != kNoCoding);
///             x_AddBoundary(TSeqPos((p - src) * m_SrcDensity), curr_type);
///             prev_type = curr_type;
///         }
///     }
///     x_AddBoundary(length, kNoCoding);
///
///     _ASSERT(m_Boundaries.at(0) == 0);
///
///     SArrangement* best_arrangement
///         = (m_EndingNarrow.cost < m_EndingWide.cost
///            ? &m_EndingNarrow : &m_EndingWide);
///
///     // XXX - fine-tune 4na/2na boundaries to minimize wastage?
///
///     size_t n = m_Boundaries.size() - 1;
///     vector<TCoding> coding_vec(n);
///     SCodings*       codings = best_arrangement->codings;
///     for (size_t i = 0;  i < n;
///          i += codings->current_used, codings = codings->previous) {
///         _ASSERT(codings != NULL);
///         memcpy(coding_vec.data() + (n - i - codings->current_used),
///                codings->current, codings->current_used * sizeof(TCoding));
///     }
///     _ASSERT(codings == NULL);
///
///     SIZE_TYPE result = 0;
///     for (size_t i = 0;  i < n;  ++i) {
///         TCoding coding  = coding_vec[i];
///         TSeqPos start   = m_Boundaries[i];
///         while (i < n - 1  &&  coding_vec[i + 1] == coding) {
///             ++i; // merge when possible
///         }
///
///         TSeqPos len = m_Boundaries[i + 1] - start;
///         char* segment = m_Target.NewSegment(coding, len);
///         if (coding == CSeqUtil::e_not_set) { // gap
///             _ASSERT(m_GapsOK);
///             result += len;
///         } else {
///             result += CSeqConvert::Convert(src, m_SrcCoding, start, len,
///                                            segment, coding);
///         }
///     }
///
///     return result;
/// }
/// ```
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
        let wide = BEST_IS_NCBI4NA[src[pos] as usize];
        let coding = if wide == 1 {
            Coding::Ncbi4na
        } else {
            Coding::Ncbi2na
        };
        packer.add_boundary(pos as u32, Some(coding));
        pos += 1 + src[pos + 1..]
            .iter()
            .position(|&residue| BEST_IS_NCBI4NA[residue as usize] != wide)
            .unwrap_or(src.len() - pos - 1);
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
const fn ncbi2na_code(residue: u8) -> u8 {
    match residue {
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 0,
    }
}

/// The NCBI4na code of the IUPACna letters (`CIupacnaTo4na`'s rows: the table has the
/// code in the high and in the low half of a byte), which is also their NCBI8na code.
const fn ncbi4na_code(residue: u8) -> u8 {
    match residue {
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
        b'N' => 15,
        _ => 0,
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
/// ```c++
/// SIZE_TYPE CSeqConvert_imp::x_ConvertIupacnaTo2na
/// (const char* src,
///  TSeqPos pos,
///  TSeqPos length,
///  char* dst)
/// {
///     // The iupacna to ncbi2na table is constructed such that each row
///     // correspond to an iupacna letter and each column corresponds to
///     // that letter being in one of the 4 possible offsets within the
///     // ncbi2na byte
///     const Uint1* table = CIupacnaTo2na::GetTable();
///
///     const char* src_i = src + pos;
///     for ( size_t count = length / 4; count; --count ) {
///         *dst =
///             table[*src_i * 4          ] |
///             table[*(src_i + 1) * 4 + 1] |
///             table[*(src_i + 2) * 4 + 2] |
///             table[*(src_i + 3) * 4 + 3];
///         src_i += 4;
///         ++dst;
///     }
///
///     // Handle overhang
///     if ( length % 4 != 0 ) {
///         *dst = 0x0;
///         for( size_t i = 0; i < (length % 4); ++i, ++src_i ) {
///             *dst |= (char)table[static_cast<Uint1>(*src_i) * 4 + i];
///         }
///     }
///
///     return length;
/// }
/// ```
pub(super) fn iupacna_to_ncbi2na(src: &[u8]) -> Vec<u8> {
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
/// ```c++
/// SIZE_TYPE CSeqConvert_imp::x_ConvertIupacnaTo4na
/// (const char* src,
///  TSeqPos pos,
///  TSeqPos length,
///  char* dst)
/// {
///
///     // The iupacna to ncbi4na table is constructed such that each row
///     // correspond to an iupacna letter and each column corresponds to
///     // that letter being in one of the 2 possible offsets within the
///     // ncbi4na byte
///     const Uint1* table = CIupacnaTo4na::GetTable();
///
///     const char* src_i = src + pos;
///
///     for ( size_t count = length / 2; count; --count ) {
///         *dst = table[*src_i * 2] | table[*(src_i + 1) * 2 + 1];
///         src_i += 2;
///         ++dst;
///     }
///
///     // handle overhang
///     if ( length % 2 != 0 ) {
///         *dst = table[static_cast<Uint1>(*src_i) * 2];
///     }
///
///     return length;
/// }
/// ```
pub(super) fn iupacna_to_ncbi4na(src: &[u8]) -> Vec<u8> {
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
/// ```c++
///             if (L.GetSeq_data().IsNcbi2na()) {
///                 CSeqConvert::Convert(L.GetSeq_data().GetNcbi2na(),
///                                      CSeqUtil::e_Ncbi2na,
///                                      0,
///                                      L.GetLength(),
///                                      na8,
///                                      CSeqUtil::e_Ncbi8na);
///             } else if (L.GetSeq_data().IsNcbi4na()) {
///                 CSeqConvert::Convert(L.GetSeq_data().GetNcbi4na(),
///                                      CSeqUtil::e_Ncbi4na,
///                                      0,
///                                      L.GetLength(),
///                                      na8,
///                                      CSeqUtil::e_Ncbi8na);
/// ```
pub(super) fn append_ncbi8na(seq8na: &mut Vec<u8>, src: &[u8]) {
    seq8na.extend(src.iter().map(|&residue| NCBI4NA_CODES[residue as usize]));
}

/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:856-874
/// ```c++
/// SIZE_TYPE CSeqConvert_imp::x_Convert8naTo4na
/// (const char* src,
///  TSeqPos pos,
///  TSeqPos length,
///  char *dst)
/// {
///     const char* iter = src + pos;
///
///     for ( size_t i = length / 2; i; --i, ++dst ) {
///         *dst = char((*iter << 4) | (*(iter + 1)));
///         iter += 2;
///     }
///
///     if ( (length % 2) != 0 ) {
///         *dst = char((*iter << 4) & 0xf0);
///     }
///
///     return length;
/// }
/// ```
pub(super) fn ncbi8na_to_ncbi4na(src: &[u8]) -> Vec<u8> {
    let mut dst = Vec::with_capacity(src.len().div_ceil(2));
    let mut chunks = src.chunks_exact(2);
    dst.extend((&mut chunks).map(|two| two[0] << 4 | two[1]));
    if let [last] = chunks.remainder() {
        dst.push(last << 4);
    }
    dst
}

#[cfg(test)]
mod tests {
    use super::*;

    fn segments(src: &[u8]) -> Vec<(Coding, u32, u32)> {
        pack_iupacna(src)
            .into_iter()
            .map(|segment| (segment.coding, segment.start, segment.length))
            .collect()
    }

    /// One run of unambiguous bases: NCBI4na for one or two bases (equal costs, the wide
    /// arrangement wins), NCBI2na from three (makeblastdb: `G` and `GA` get an ambiguity
    /// block, `GAT` gets none).
    #[test]
    fn one_unambiguous_run() {
        assert_eq!(segments(b""), []);
        assert_eq!(segments(b"G"), [(Coding::Ncbi4na, 0, 1)]);
        assert_eq!(segments(b"GA"), [(Coding::Ncbi4na, 0, 2)]);
        assert_eq!(segments(b"UU"), [(Coding::Ncbi4na, 0, 2)]);
        for length in 3..200 {
            let src = b"GATC".repeat(50);
            assert_eq!(
                segments(&src[..length]),
                [(Coding::Ncbi2na, 0, length as u32)]
            );
        }
    }

    /// A segment costs 128 bytes, so a run of unambiguous bases is its own NCBI2na
    /// segment only where that saves more than the segments it adds.
    #[test]
    fn mixed_runs() {
        assert_eq!(segments(b"N"), [(Coding::Ncbi4na, 0, 1)]);
        assert_eq!(segments(b"ACGTN"), [(Coding::Ncbi4na, 0, 5)]);
        assert_eq!(segments(b"NACGT"), [(Coding::Ncbi4na, 0, 5)]);
        // 600 bases: 150 bytes of NCBI2na and a segment (128) against 300 of NCBI4na.
        let mut src = b"ACGT".repeat(150);
        src.push(b'N');
        assert_eq!(
            segments(&src),
            [(Coding::Ncbi2na, 0, 600), (Coding::Ncbi4na, 600, 1)]
        );
        // 400 bases: 100 + 128 against 200.
        let mut src = b"ACGT".repeat(100);
        src.push(b'N');
        assert_eq!(segments(&src), [(Coding::Ncbi4na, 0, 401)]);
        // The break-even: 128 + n/4 against n/2 for the first run.
        let mut src = b"ACGT".repeat(128);
        src.push(b'N');
        assert_eq!(segments(&src), [(Coding::Ncbi4na, 0, 513)]);
        let mut src = b"ACGT".repeat(129);
        src.push(b'N');
        assert_eq!(
            segments(&src),
            [(Coding::Ncbi2na, 0, 516), (Coding::Ncbi4na, 516, 1)]
        );
        // A run between two ambiguous ones must pay for two segments.
        let mut src = b"N".to_vec();
        src.extend(b"ACGT".repeat(200));
        src.push(b'N');
        assert_eq!(segments(&src), [(Coding::Ncbi4na, 0, 802)]);
        let mut src = b"N".to_vec();
        src.extend(b"ACGT".repeat(300));
        src.push(b'N');
        assert_eq!(
            segments(&src),
            [
                (Coding::Ncbi4na, 0, 1),
                (Coding::Ncbi2na, 1, 1200),
                (Coding::Ncbi4na, 1201, 1)
            ]
        );
    }

    /// `CPacker` against the least cost over every assignment of a coding to each run
    /// (an NCBI4na-only run must be NCBI4na; a segment costs 128 plus its bytes).
    #[test]
    fn least_cost_of_all_arrangements() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % bound
        };
        for _ in 0..300 {
            let runs = 1 + next(9) as usize;
            let mut src = Vec::new();
            let mut lengths = Vec::new();
            let mut wide = next(2) == 1;
            for _ in 0..runs {
                let length = 1 + next(if wide { 40 } else { 1400 }) as usize;
                let letters: &[u8] = if wide { b"NRYKM" } else { b"ACGTU" };
                for _ in 0..length {
                    src.push(letters[next(5) as usize]);
                }
                lengths.push((wide, length as u32));
                wide = !wide;
            }
            // The cost that `x_AddBoundary` sums: the bytes of each run in its coding, and
            // one overhead for each segment (adjacent runs of one coding are one segment).
            let cost = |codings: &[Coding]| {
                let mut cost = 0;
                for (run, &coding) in codings.iter().enumerate() {
                    if run == 0 || codings[run - 1] != coding {
                        cost += OVERHEAD;
                    }
                    cost += bytes_needed(coding, lengths[run].1);
                }
                cost
            };
            let mut least = u64::MAX;
            for mask in 0..1u32 << runs {
                let codings: Vec<Coding> = (0..runs)
                    .map(|run| {
                        if lengths[run].0 || mask >> run & 1 == 1 {
                            Coding::Ncbi4na
                        } else {
                            Coding::Ncbi2na
                        }
                    })
                    .collect();
                least = least.min(cost(&codings));
            }
            let packed = pack_iupacna(&src);
            let mut end = 0;
            for (index, segment) in packed.iter().enumerate() {
                assert_eq!(segment.start, end);
                assert!(index == 0 || packed[index - 1].coding != segment.coding);
                end += segment.length;
            }
            let mut run_start = 0;
            let packed_codings: Vec<Coding> = lengths
                .iter()
                .map(|&(_, length)| {
                    let segment = packed
                        .iter()
                        .find(|segment| run_start < segment.start + segment.length)
                        .unwrap();
                    run_start += length;
                    assert!(run_start <= segment.start + segment.length);
                    segment.coding
                })
                .collect();
            for (run, &(wide, _)) in lengths.iter().enumerate() {
                assert!(!wide || packed_codings[run] == Coding::Ncbi4na);
            }
            assert_eq!(cost(&packed_codings), least, "{lengths:?}");
            assert_eq!(end as usize, src.len());
        }
    }

    #[test]
    fn conversions() {
        assert_eq!(iupacna_to_ncbi2na(b"ACGT"), [0x1b]);
        assert_eq!(iupacna_to_ncbi2na(b"ACGTU"), [0x1b, 0xc0]);
        assert_eq!(iupacna_to_ncbi2na(b"GAT"), [0x8c]);
        assert_eq!(iupacna_to_ncbi2na(b"TTTTTTT"), [0xff, 0xfc]);
        assert_eq!(iupacna_to_ncbi4na(b"G"), [0x40]);
        assert_eq!(iupacna_to_ncbi4na(b"GA"), [0x41]);
        assert_eq!(
            iupacna_to_ncbi4na(b"ACGTUMRWSYKVHDBN"),
            [0x12, 0x48, 0x83, 0x59, 0x6a, 0xc7, 0xbd, 0xef]
        );
        assert_eq!(iupacna_to_ncbi4na(b"ACN"), [0x12, 0xf0]);
        assert_eq!(ncbi8na_to_ncbi4na(&[1, 2, 15]), [0x12, 0xf0]);
        assert_eq!(ncbi8na_to_ncbi4na(&[]), [0u8; 0]);
    }

    /// The NCBI8na codes that `s_FixBioseqDeltas` reads from a literal equal the NCBI4na
    /// codes of its letters, whether the literal is NCBI2na (`C2naTo8na`: 0, 1, 2, 3 are
    /// 1, 2, 4, 8) or NCBI4na.
    #[test]
    fn ncbi8na_of_both_codings() {
        for &base in b"ACGTU" {
            assert_eq!(
                1 << ncbi2na_code(base),
                ncbi4na_code(base),
                "{}",
                base as char
            );
        }
        let mut seq8na = Vec::new();
        append_ncbi8na(&mut seq8na, b"ACGTUN");
        assert_eq!(seq8na, [1, 2, 4, 8, 8, 15]);
    }
}
