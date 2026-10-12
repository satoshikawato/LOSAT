//! Sequence data of a database volume: the conversions of
//! `c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp` from a
//! `CSeq_data` to the bytes of the sequence file, the ambiguity data of a
//! nucleotide sequence, and the generator (`CRandom`) that picks the base
//! stored in place of an ambiguous one.

use super::bioseq::SeqData;

// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:159-165
// ```c++
//     enum {
//         kStateSize = 33  // size of state array
//     };
//     EGetRandMethod  m_RandMethod;
//     TValue          m_State[kStateSize];
//     int             m_RJ;
//     int             m_RK;
// ```
// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:98
// ```c++
// static const size_t kStateOffset = 12;
// ```
const STATE_SIZE: usize = 33;
const STATE_OFFSET: usize = 12;

/// `CRandom` with the lagged Fibonacci generator (`eGetRand_LFG`).
#[derive(Debug, Clone)]
pub(crate) struct Random {
    state: [u32; STATE_SIZE],
    rj: usize,
    rk: usize,
}

impl Random {
    // NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:227-231
    // ```c++
    // CRandom::CRandom(TValue seed) :
    //     m_RandMethod(eGetRand_LFG)
    // {
    //     SetSeed(seed);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:296-308
    // ```c++
    //     m_State[0] = m_Seed = seed;
    //
    //     // linear congruential initializer
    //     for (int i = 1;  i < kStateSize;  ++i) {
    //         m_State[i] = 1103515245 * m_State[i-1] + 12345;
    //     }
    //
    //     m_RJ = kStateOffset;
    //     m_RK = kStateSize - 1;
    //
    //     for (int i = 0;  i < 10 * kStateSize;  ++i) {
    //         GetRand();
    //     }
    // ```
    pub(crate) fn new(seed: u32) -> Self {
        let mut state = [0u32; STATE_SIZE];
        state[0] = seed;
        for i in 1..STATE_SIZE {
            state[i] = 1_103_515_245u32
                .wrapping_mul(state[i - 1])
                .wrapping_add(12345);
        }
        let mut random = Self {
            state,
            rj: STATE_OFFSET,
            rk: STATE_SIZE - 1,
        };
        for _ in 0..10 * STATE_SIZE {
            random.get_rand();
        }
        random
    }

    // NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:222-241
    // ```c++
    //     TValue r;
    //
    //     r = m_State[m_RK] + m_State[m_RJ--];
    //     m_State[m_RK--] = r;
    //
    //     if ( m_RK < 0 ) {
    //         m_RK = kStateSize - 1;
    //     }
    //     else if ( m_RJ < 0 ) {
    //         m_RJ = kStateSize - 1;
    //     }
    //
    //     return r;
    // }
    //
    //
    // inline CRandom::TValue CRandom::GetRand(void)
    // {
    //     return x_GetRand32Bits() >> 1; // discard the least-random bit
    // }
    // ```
    // `m_RJ` and `m_RK` are never negative at the same step (they start 21
    // apart and move together), so the unsigned indices wrap one at a time
    // exactly as the signed ones do.
    pub(crate) fn get_rand(&mut self) -> u32 {
        let r = self.state[self.rk].wrapping_add(self.state[self.rj]);
        self.state[self.rk] = r;
        if self.rk == 0 {
            self.rk = STATE_SIZE - 1;
            self.rj -= 1;
        } else {
            self.rk -= 1;
            self.rj = if self.rj == 0 {
                STATE_SIZE - 1
            } else {
                self.rj - 1
            };
        }
        r >> 1
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:141-147
// ```c++
//     enum {
//         eMaxLength = 0xFFF ///< Maximum length of a region.
//     };
//
//     int m_Start; ///< Starting offset (offset of first base).
//     int m_End;   ///< End offset (offset of first disincluded base.)
//     int m_Value; ///< Value of base (ambiguity letter).
// ```
const MAX_REGION_LENGTH: u32 = 0xFFF;

/// `CAmbiguousRegion`: a run of one ambiguity letter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AmbiguousRegion {
    start: u32,
    end: u32,
    value: u8,
}

impl AmbiguousRegion {
    fn length(&self) -> u32 {
        self.end - self.start
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:92-104
    // ```c++
    //     bool Append(int value, int offset)
    //     {
    //         _ASSERT(m_End && m_Value);
    //
    //         if (value == m_Value && offset == m_End) {
    //             if (Length() < eMaxLength) {
    //                 m_End ++;
    //                 return true;
    //             }
    //         }
    //
    //         return false;
    //     }
    // ```
    fn append(&mut self, value: u8, offset: u32) -> bool {
        if value == self.value && offset == self.end && self.length() < MAX_REGION_LENGTH {
            self.end += 1;
            return true;
        }
        false
    }
}

/// `CAmbigDataBuilder`: the ambiguity data of one nucleotide sequence.
struct AmbigDataBuilder {
    size: u32,
    regions: Vec<AmbiguousRegion>,
    random: Random,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:161-169
// ```c++
//         for(int i = 0; i < 16; i++) {
//             m_Log2[i] = -1;
//         }
//
//         // Only these values should be specified
//         m_Log2[1] = 0;
//         m_Log2[2] = 1;
//         m_Log2[4] = 2;
//         m_Log2[8] = 3;
// ```
const LOG2: [i8; 16] = [-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1];

impl AmbigDataBuilder {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:158-160
    // ```c++
    //     CAmbigDataBuilder(int sz)
    //         : m_Size(sz), m_Random(sz)
    //     {
    // ```
    // The generator is seeded with the length of the sequence.
    fn new(size: u32) -> Self {
        Self {
            size,
            regions: Vec::new(),
            random: Random::new(size),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:186-213
    // ```c++
    //     int Check(int data, int offset)
    //     {
    //         // If offset is past the sequence end, return zero.  This can
    //         // happen since the calling code checks for ambiguities two
    //         // letters at a time.
    //
    //         if (offset >= m_Size) {
    //             return 0;
    //         }
    //
    //         // If we recieve a non-ambiguity, return the normal
    //         // translation for this letter code.
    //
    //         _ASSERT(data != 0);
    //
    //         int rv = m_Log2[data & 0xF];
    //
    //         if (rv != -1) {
    //             return rv;
    //         }
    //
    //         // Bona-fide ambiguity; we need to make up some random junk,
    //         // and also build an ambiguous region.
    //
    //         x_AddAmbig(data, offset);
    //
    //         return x_Random(data);
    //     }
    // ```
    fn check(&mut self, data: u8, offset: u32) -> u8 {
        if offset >= self.size {
            return 0;
        }
        let rv = LOG2[usize::from(data & 0xF)];
        if rv != -1 {
            return rv as u8;
        }
        self.add_ambig(data, offset);
        self.random_base(data)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:349-359
    // ```c++
    //     void x_AddAmbig(int value, int offset)
    //     {
    //         if (m_Regions.size()) {
    //             if (m_Regions.back().Append(value, offset)) {
    //                 return;
    //             }
    //         }
    //
    //         CAmbiguousRegion r(value, offset);
    //         m_Regions.push_back(r);
    //     }
    // ```
    fn add_ambig(&mut self, value: u8, offset: u32) {
        if let Some(last) = self.regions.last_mut() {
            if last.append(value, offset) {
                return;
            }
        }
        self.regions.push(AmbiguousRegion {
            start: offset,
            end: offset + 1,
            value,
        });
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:369-410
    // ```c++
    //     int x_Random(int value)
    //     {
    //         // 0xF is the most common case of ambiguity so it's worth to
    //         // process it in fastest way, especially because it's very easy.
    //
    //         if (value == 0xF) {
    //             return m_Random.GetRand() & 0x3;
    //         }
    //
    //         if (value == 0) {
    //             cerr << "Error: '0' ambiguity code found, changing to 15." << endl;
    //             return m_Random.GetRand() & 0x3;
    //         }
    //
    //         int bitcount = ((value & 1) +
    //                         ((value >> 1) & 1) +
    //                         ((value >> 2) & 1) +
    //                         ((value >> 3) & 1));
    //
    //         // 1-bit ambiguities here, indicate an error in this class.
    //         _ASSERT(bitcount >= 2);
    //         _ASSERT(bitcount <= 3);
    //
    //         int pick = m_Random.GetRand() % bitcount;
    //
    //         for(int i = 0; i < 4; i++) {
    //             // skip 0 bits in input.
    //             if ((value & (1 << i)) == 0)
    //                 continue;
    //
    //             // If the bitcount is zero, this is the bit we want.
    //             if (! pick)
    //                 return i;
    //
    //             // Else, decrement.
    //             pick--;
    //         }
    //
    //         // This should be unreachable.
    //         _ASSERT(0);
    //         return 0;
    //     }
    // ```
    // The value 0 never reaches this function: the caller turns a gap into
    // 0xF, and `check` returns before it for the padding half byte past the
    // end of the sequence.
    fn random_base(&mut self, value: u8) -> u8 {
        if value == 0xF || value == 0 {
            return (self.random.get_rand() & 0x3) as u8;
        }
        let bitcount = u32::from(value & 1)
            + u32::from((value >> 1) & 1)
            + u32::from((value >> 2) & 1)
            + u32::from((value >> 3) & 1);
        let mut pick = self.random.get_rand() % bitcount;
        for i in 0..4u8 {
            if value & (1 << i) == 0 {
                continue;
            }
            if pick == 0 {
                return i;
            }
            pick -= 1;
        }
        0
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:243-283
    // ```c++
    //         bool new_format = false;
    //
    //         // Current technique: If the sequence length is longer than
    //         // 2^24-1, I use the new format.  Else, if any region of
    //         // ambiguities is longer than 16 bases, I use the new format.
    //         // Otherwise, I use the old format.
    //
    //         if (m_Size > 0xFFFFFF) {
    //             new_format = true;
    //         } else {
    //             for(unsigned i = 0; i<m_Regions.size(); i++) {
    //                 if (m_Regions[i].Length() > 0xF) {
    //                     new_format = true;
    //                     break;
    //                 }
    //             }
    //         }
    //
    //         int num_amb = (int) m_Regions.size();
    //
    //         // The size packed here is actually the number of *words* used
    //         // rather than the number of ambiguous regions; the new format
    //         // two words per region.
    //
    //         Uint4 amb_words = (new_format
    //                            ? (0x80000000 | (num_amb * 2))
    //                            : num_amb);
    //
    //         amb.reserve((1 + m_Regions.size())*8);
    //         s_AppendInt4(amb, amb_words);
    //
    //         for(int i = 0; i < num_amb; i++) {
    //             // Regions over 4k are split during ambiguity scanning; at
    //             // this point in the code, splitting has already happened.
    //
    //             if (new_format) {
    //                 x_PackNewAmbig(amb, m_Regions[i]);
    //             } else {
    //                 x_PackOldAmbig(amb, m_Regions[i]);
    //             }
    //         }
    // ```
    fn get_ambig(&self) -> Vec<u8> {
        let new_format =
            self.size > 0xFF_FFFF || self.regions.iter().any(|region| region.length() > 0xF);
        let num_amb = self.regions.len() as u32;
        let amb_words = if new_format {
            0x8000_0000 | num_amb.wrapping_mul(2)
        } else {
            num_amb
        };
        let mut amb = Vec::with_capacity((1 + self.regions.len()) * 8);
        amb.extend_from_slice(&amb_words.to_be_bytes());
        for region in &self.regions {
            if new_format {
                pack_new_ambig(&mut amb, region);
            } else {
                pack_old_ambig(&mut amb, region);
            }
        }
        amb
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:292-311
// ```c++
//         int length_m1 = r.Length() - 1;
//
//         _ASSERT(r.Value() <= 15);
//         _ASSERT((length_m1 >> 12) == 0);
//
//         // First word - residue value and run length.
//         char A1[4];
//
//         char ch0 = (r.Value() << 4) | (length_m1 >> 8);
//         char ch1 = length_m1 & 0xFF;
//
//         A1[0] = ch0;
//         A1[1] = ch1;
//         A1[2] = A1[3] = 0; // unused space
//
//         amb.append(A1, 4);
//
//         // Second word - starting offset.
//
//         s_AppendInt4(amb, r.Offset());
// ```
fn pack_new_ambig(amb: &mut Vec<u8>, region: &AmbiguousRegion) {
    let length_m1 = region.length() - 1;
    amb.push((region.value << 4) | (length_m1 >> 8) as u8);
    amb.push((length_m1 & 0xFF) as u8);
    amb.push(0);
    amb.push(0);
    amb.extend_from_slice(&region.start.to_be_bytes());
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:292-336
// ```c++
//         int length_m1 = r.Length() - 1;
//         int off = r.Offset();
//
//         _ASSERT(r.Value() <= 15);
//         _ASSERT((length_m1 >> 4) == 0);
//         _ASSERT(off <= 0xFFFFFF); // old form uses three byte offsets
//
//         // All in one word.
//         char A1[4];
//
//         char ch0 = (r.Value() << 4) | length_m1;
//
//         A1[0] = ch0;
//         A1[1] = (off >> 16) & 0xFF;
//         A1[2] = (off >>  8) & 0xFF;
//         A1[3] = (off      ) & 0xFF;
//
//         amb.append(A1, 4);
// ```
fn pack_old_ambig(amb: &mut Vec<u8>, region: &AmbiguousRegion) {
    let length_m1 = region.length() - 1;
    let off = region.start;
    amb.push((region.value << 4) | length_m1 as u8);
    amb.push(((off >> 16) & 0xFF) as u8);
    amb.push(((off >> 8) & 0xFF) as u8);
    amb.push((off & 0xFF) as u8);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:434-441
// ```c++
//     vector<unsigned char> ctable;
//     ctable.resize(16, 0xFF);
//
//     for(int i = 0; i<4; i++) {
//         ctable[1 << i] = i;
//     }
//
//     return ctable;
// ```
const NA4_TO_NA2: [u8; 16] = [
    0xFF, 0, 1, 0xFF, 2, 0xFF, 0xFF, 0xFF, 3, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
];

/// The bytes of one sequence in the sequence file, and its ambiguity data.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BinarySequence {
    /// `m_Sequence`: NCBIstdaa residues, or 2-bit bases whose last byte
    /// holds the number of bases in it.
    pub seq: Vec<u8>,
    /// `m_Ambig`: the ambiguity data of a nucleotide sequence; empty when
    /// the `CSeq_data` is ncbi2na, and for a protein.
    pub ambig: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:456-517
// ```c++
//     int inp_bytes   = s_DivideRoundUp(base_length, 2);
//     int blast_bytes = base_length / 4 + 1;
//     int remainder   = base_length & 3;
//     int last_byte   = blast_bytes - 1;
//
//     // Accumulator for ambiguity data.
//     CAmbigDataBuilder ambiguities(base_length);
//
//     if (!((int)inp_bytes == (int)byte_length)) {
//         cout << "ib=" << inp_bytes << ",n4sz=" << byte_length << endl;
//     }
//
//     _ASSERT((int)inp_bytes == (int)byte_length);
//
//     seq.resize(blast_bytes);
//     ...
//     for(int i = 0; i < inp_bytes; i++) {
//         // one input byte
//         uchar inp = ncbi4na[i];
//
//         // represents 2 bases
//         uchar b1 = inp >>  4;
//         uchar b2 = inp & 0xF;
//
//         // compress each to 2 bits
//         uchar c1 = ctable[b1];
//         uchar c2 = ctable[b2];
//
//         uchar half = 0;
//
//         if (((c1 | c2) & 0x80) == 0) {
//             // No ambiguities, so we can do this the easy way.
//
//             half = (c1 << 2) | c2;
//         } else {
//             // Check each element, accumulate ambiguity data.
//
//             if (! b1) {
//                 b1 = 0xF; // replace gap with 'N'
//             }
//
//             if (! b2 && (i*2+1) < base_length) {
//                 b2 = 0xF; // the last base is the sentinel
//             }
//
//             half |= ambiguities.Check(b1, i*2) << 2;
//             half |= ambiguities.Check(b2, i*2+1);
//         }
//
//         seq[i/2] |= (i & 1) ? half : (half << 4);
//     }
//     seq[last_byte] &= 255-3;
//     seq[last_byte] |= remainder;
//
//     ambiguities.GetAmbig(amb);
// ```
// `ncbi4na` holds `ceil(base_length / 2)` bytes: the Seq-data of the FASTA
// source always does (NCBI only prints a line otherwise and reads on).
/// `WriteDB_Ncbi4naToBinary`.
pub(crate) fn ncbi4na_to_binary(ncbi4na: &[u8], base_length: u32) -> BinarySequence {
    let inp_bytes = base_length.div_ceil(2) as usize;
    let blast_bytes = (base_length / 4 + 1) as usize;
    let remainder = (base_length & 3) as u8;
    let last_byte = blast_bytes - 1;
    let mut ambiguities = AmbigDataBuilder::new(base_length);
    let mut seq = vec![0u8; blast_bytes];
    for (i, &inp) in ncbi4na.iter().enumerate().take(inp_bytes) {
        let mut b1 = inp >> 4;
        let mut b2 = inp & 0xF;
        let c1 = NA4_TO_NA2[usize::from(b1)];
        let c2 = NA4_TO_NA2[usize::from(b2)];
        let mut half = 0u8;
        if (c1 | c2) & 0x80 == 0 {
            half = (c1 << 2) | c2;
        } else {
            if b1 == 0 {
                b1 = 0xF;
            }
            let second = (i as u32) * 2 + 1;
            if b2 == 0 && second < base_length {
                b2 = 0xF;
            }
            half |= ambiguities.check(b1, second - 1) << 2;
            half |= ambiguities.check(b2, second);
        }
        seq[i / 2] |= if i & 1 == 1 { half } else { half << 4 };
    }
    seq[last_byte] &= 255 - 3;
    seq[last_byte] |= remainder;
    BinarySequence {
        seq,
        ambig: ambiguities.get_ambig(),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:572-587
// ```c++
//     int base_length = si.GetLength();
//     int data_bytes  = s_DivideRoundUp(base_length, 4);
//     int blast_bytes = base_length / 4 + 1;
//     int remainder   = base_length & 3;
//     int last_byte   = blast_bytes - 1;
//
//     const vector<char> & v = si.GetSeq_data().GetNcbi2na().Get();
//
//     _ASSERT((int)data_bytes == (int)v.size());
//
//     seq.reserve(blast_bytes);
//     seq.assign(& v[0], data_bytes);
//     seq.resize(blast_bytes);
//
//     seq[last_byte] &= 255-3;
//     seq[last_byte] |= remainder;
// ```
/// `WriteDB_Ncbi2naToBinary`: no ambiguity data.
pub(crate) fn ncbi2na_to_binary(ncbi2na: &[u8], base_length: u32) -> BinarySequence {
    let data_bytes = base_length.div_ceil(4) as usize;
    let blast_bytes = (base_length / 4 + 1) as usize;
    let remainder = (base_length & 3) as u8;
    let last_byte = blast_bytes - 1;
    let mut seq = Vec::with_capacity(blast_bytes);
    seq.extend_from_slice(&ncbi2na[..data_bytes.min(ncbi2na.len())]);
    seq.resize(blast_bytes, 0);
    seq[last_byte] &= 255 - 3;
    seq[last_byte] |= remainder;
    BinarySequence {
        seq,
        ambig: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:391-394
// ```c++
// // IUPACna to NCBI4na
// const Uint1 CIupacnaTo4na::scm_Table[512] = {
//     0xf0,0xf,  0xf0,0xf,  0xf0,0xf,  0xf0,0xf,
//     0xf0,0xf,  0xf0,0xf,  0xf0,0xf,  0xf0,0xf,
// ```
// The table of 256 pairs (the letter as the first and as the second base of
// a byte) is 0xf0, 0xf (N) for every byte except the letters below, in
// either case (sequtil_tables.cpp:409-423, the rows of 'A'-'Z' and 'a'-'z'):
// A 1, B 14, C 2, D 13, G 4, H 11, K 12, M 3, R 5, S 6, T 8, U 8, V 7, W 9,
// Y 10.
const fn iupacna_to_ncbi4na_table() -> [u8; 256] {
    let mut table = [0xFu8; 256];
    let letters = b"ABCDGHKMRSTUVWY";
    let codes = [1u8, 14, 2, 13, 4, 11, 12, 3, 5, 6, 8, 8, 7, 9, 10];
    let mut i = 0;
    while i < letters.len() {
        table[letters[i] as usize] = codes[i];
        table[letters[i].to_ascii_lowercase() as usize] = codes[i];
        i += 1;
    }
    table
}
const IUPACNA_TO_NCBI4NA: [u8; 256] = iupacna_to_ncbi4na_table();

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:324-339
// ```c++
//     const Uint1* table = CIupacnaTo4na::GetTable();
//
//     const char* src_i = src + pos;
//
//     for ( size_t count = length / 2; count; --count ) {
//         *dst = table[*src_i * 2] | table[*(src_i + 1) * 2 + 1];
//         src_i += 2;
//         ++dst;
//     }
//
//     // handle overhang
//     if ( length % 2 != 0 ) {
//         *dst = table[static_cast<Uint1>(*src_i) * 2];
//     }
//
//     return length;
// ```
/// `CSeqConvert::Convert` from iupacna to ncbi4na: two bases per byte, the
/// first in the high half.
pub(crate) fn iupacna_to_ncbi4na(iupacna: &[u8]) -> Vec<u8> {
    let mut na4 = Vec::with_capacity(iupacna.len().div_ceil(2));
    let mut pairs = iupacna.chunks_exact(2);
    for pair in &mut pairs {
        na4.push(
            (IUPACNA_TO_NCBI4NA[usize::from(pair[0])] << 4)
                | IUPACNA_TO_NCBI4NA[usize::from(pair[1])],
        );
    }
    if let [last] = pairs.remainder() {
        na4.push(IUPACNA_TO_NCBI4NA[usize::from(*last)] << 4);
    }
    na4
}

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1747-1749
// ```c++
// // NCBIeaa to NCBIstdaa
// const Uint1 CEaaToStdaa::scm_Table[256] = {
//     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
// ```
// The table is 0xff except (sequtil_tables.cpp:1754-1760): '*' 25, '-' 0,
// and the capital letters A 1, B 2, C 3, D 4, E 5, F 6, G 7, H 8, I 9, J 27,
// K 10, L 11, M 12, N 13, O 26, P 14, Q 15, R 16, S 17, T 18, U 24, V 19,
// W 20, X 21, Y 22, Z 23.
const STDAA_OF_LETTER: [u8; 26] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 27, 10, 11, 12, 13, 26, 14, 15, 16, 17, 18, 24, 19, 20, 21, 22, 23,
];

const fn ncbieaa_to_ncbistdaa_table() -> [u8; 256] {
    let mut table = [0xFFu8; 256];
    table[b'*' as usize] = 25;
    table[b'-' as usize] = 0;
    let mut i = 0;
    while i < 26 {
        table[b'A' as usize + i] = STDAA_OF_LETTER[i];
        i += 1;
    }
    table
}
const NCBIEAA_TO_NCBISTDAA: [u8; 256] = ncbieaa_to_ncbistdaa_table();

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1673-1675
// ```c++
// // IUPACaa to NCBIstdaa
// const Uint1 CIupacaaToStdaa::scm_Table[256] = {
//     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
// ```
// The table is 0xff except the letters, in either case, with the codes of
// the ncbieaa table (sequtil_tables.cpp:1683-1690).
const fn iupacaa_to_ncbistdaa_table() -> [u8; 256] {
    let mut table = [0xFFu8; 256];
    let mut i = 0;
    while i < 26 {
        table[b'A' as usize + i] = STDAA_OF_LETTER[i];
        table[b'a' as usize + i] = STDAA_OF_LETTER[i];
        i += 1;
    }
    table
}
const IUPACAA_TO_NCBISTDAA: [u8; 256] = iupacaa_to_ncbistdaa_table();

fn convert_1_to_1(src: &[u8], table: &[u8; 256]) -> Vec<u8> {
    src.iter().map(|&byte| table[usize::from(byte)]).collect()
}

/// The `Seq-data` choices CWriteDB cannot write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnprocessableSequence;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1007-1036
// ```c++
//         switch(sd.Which()) {
//         case CSeq_data::e_Ncbistdaa:
//             WriteDB_StdaaToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Ncbieaa:
//             WriteDB_EaaToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Iupacaa:
//             WriteDB_IupacaaToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Ncbi2na:
//             WriteDB_Ncbi2naToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Ncbi4na:
//             WriteDB_Ncbi4naToBinary(si, m_Sequence, m_Ambig);
//             break;
//
//         case CSeq_data::e_Iupacna:
//              WriteDB_IupacnaToBinary(si, m_Sequence, m_Ambig);
//              break;
//
//         default:
//             msg = "Unable to process sequence for entry [";
//             msg += (m_Bioseq->GetId().front())->GetSeqIdString(false);
//             msg += "].";
//         }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:530-567
// ```c++
// void WriteDB_StdaaToBinary(const CSeq_inst & si, string & seq)
// {
//     // No conversion is actually done here.
//     const vector<char> & v = si.GetSeq_data().GetNcbistdaa().Get();
//
//     _ASSERT(si.GetLength() == v.size());
//     seq.assign(& v[0], v.size());
// }
//
// void WriteDB_EaaToBinary(const CSeq_inst & si, string & seq)
// {
//     const string & v = si.GetSeq_data().GetNcbieaa().Get();
//
//     _ASSERT(si.GetLength() == v.size());
//
//     // convert to string.
//     CSeqConvert::Convert(v,
//                          CSeqUtil::e_Ncbieaa,
//                          0,
//                          (int) v.size(),
//                          seq,
//                          CSeqUtil::e_Ncbistdaa);
// }
//
// void WriteDB_IupacaaToBinary(const CSeq_inst & si, string & seq)
// {
//     const string & v = si.GetSeq_data().GetIupacaa().Get();
//
//     _ASSERT(si.GetLength() == v.size());
//
//     // convert to string.
//     CSeqConvert::Convert(v,
//                          CSeqUtil::e_Iupacaa,
//                          0,
//                          (int) v.size(),
//                          seq,
//                          CSeqUtil::e_Ncbistdaa);
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:592-609
// ```c++
//     const string & v = si.GetSeq_data().GetIupacna().Get();
//
//     _ASSERT(si.GetLength() == v.size());
//
//     string tmp;
//     // convert to string.
//     CSeqConvert::Convert(v,
//                          CSeqUtil::e_Iupacna,
//                          0,
//                          (int) v.size(),
//                          tmp,
//                          CSeqUtil::e_Ncbi4na);
//
//     WriteDB_Ncbi4naToBinary(tmp.c_str(),
//                             (int) tmp.size(),
//                             (int) si.GetLength(),
//                             seq,
//                             amb);
// ```
/// `CWriteDB_Impl::x_CookSequence` for a `CBioseq` with `Seq-data`: the
/// sequence bytes and the ambiguity data of a sequence of `length` letters.
pub(crate) fn seq_data_to_binary(
    seq_data: &SeqData,
    length: u32,
) -> Result<BinarySequence, UnprocessableSequence> {
    let protein = |seq: Vec<u8>| BinarySequence {
        seq,
        ambig: Vec::new(),
    };
    Ok(match seq_data {
        SeqData::Ncbistdaa(v) => protein(v.clone()),
        SeqData::Ncbieaa(v) => protein(convert_1_to_1(v, &NCBIEAA_TO_NCBISTDAA)),
        SeqData::Iupacaa(v) => protein(convert_1_to_1(v, &IUPACAA_TO_NCBISTDAA)),
        SeqData::Ncbi2na(v) => ncbi2na_to_binary(v, length),
        SeqData::Ncbi4na(v) => ncbi4na_to_binary(v, length),
        SeqData::Iupacna(v) => ncbi4na_to_binary(&iupacna_to_ncbi4na(v), length),
        SeqData::Other => return Err(UnprocessableSequence),
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_general.cpp:80-88
// ```c++
// int WriteDB_FindSequenceLength(bool protein, const string & seq)
// {
//     if (protein) {
//         return seq.size();
//     }
//
//     int wholebytes = (int) seq.size() - 1;
//     return (wholebytes << 2) + (seq[wholebytes] & 0x3);
// }
// ```
/// `WriteDB_FindSequenceLength`: the number of letters of a sequence in the
/// form of the sequence file (`seq` is not empty).
pub(crate) fn find_sequence_length(protein: bool, seq: &[u8]) -> u32 {
    if protein {
        return seq.len() as u32;
    }
    let wholebytes = seq.len() - 1;
    ((wholebytes as u32) << 2) + u32::from(seq[wholebytes] & 0x3)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iupacna(letters: &[u8]) -> BinarySequence {
        seq_data_to_binary(&SeqData::Iupacna(letters.to_vec()), letters.len() as u32).unwrap()
    }

    fn hex(bytes: &[u8]) -> String {
        bytes.iter().map(|byte| format!("{byte:02x}")).collect()
    }

    // Bytes that makeblastdb 2.17.0 wrote for these sequences (one-sequence
    // and several-sequence databases of 2026-10-12; the generator is seeded
    // with the sequence length, so each expectation fixes its output).
    #[test]
    fn ambiguous_bases_take_makeblastdb_s_random_bases() {
        // "ACGTNACGTN"
        let b = iupacna(b"ACGTNACGTN");
        assert_eq!(hex(&b.seq), "1b86c2");
        assert_eq!(hex(&b.ambig), "00000002f0000004f0000009");
        // "N" and "R", one base each
        let b = iupacna(b"N");
        assert_eq!(
            (hex(&b.seq), hex(&b.ambig)),
            ("81".into(), "00000001f0000000".into())
        );
        let b = iupacna(b"R");
        assert_eq!(
            (hex(&b.seq), hex(&b.ambig)),
            ("01".into(), "0000000150000000".into())
        );
        // "ACGTNN"
        let b = iupacna(b"ACGTNN");
        assert_eq!(
            (hex(&b.seq), hex(&b.ambig)),
            ("1b22".into(), "00000001f1000004".into())
        );
        // every IUPAC letter
        let b = iupacna(b"ACGTRYKMSWBDHVN");
        assert_eq!(hex(&b.seq), "1b994b0b");
        assert_eq!(
            hex(&b.ambig),
            "0000000b50000004a0000005c0000006300000076000000890000009e000000ad000000bb000000c7000000df000000e"
        );
    }

    #[test]
    fn ambiguity_at_the_end_of_odd_and_even_lengths() {
        let cases: [(&[u8], &str, &str); 5] = [
            (b"ACGTN", "1b01", "00000001f0000004"),
            (b"ACGN", "1900", "00000001f0000003"),
            (b"ACN", "17", "00000001f0000002"),
            (b"AN", "22", "00000001f0000001"),
            (b"NA", "82", "00000001f0000000"),
        ];
        for (letters, seq, ambig) in cases {
            let b = iupacna(letters);
            assert_eq!(
                (hex(&b.seq), hex(&b.ambig)),
                (seq.into(), ambig.into()),
                "{letters:?}"
            );
        }
    }

    #[test]
    fn a_text_sequence_without_ambiguity_gets_an_empty_ambiguity_block() {
        // makeblastdb keeps "G" and "GA" as text Seq-data: 4 bytes of count 0.
        let b = iupacna(b"G");
        assert_eq!(
            (hex(&b.seq), hex(&b.ambig)),
            ("81".into(), "00000000".into())
        );
        let b = iupacna(b"GA");
        assert_eq!(
            (hex(&b.seq), hex(&b.ambig)),
            ("82".into(), "00000000".into())
        );
    }

    #[test]
    fn ncbi2na_has_no_ambiguity_data() {
        // "GAT", "GATC", "GATCA", "GATCGATCG" as ncbi2na
        let cases: [(&[u8], u32, &str); 4] = [
            (&[0x8c], 3, "8f"),
            (&[0x8d], 4, "8d00"),
            (&[0x8d, 0x00], 5, "8d01"),
            (&[0x8d, 0x8d, 0x80], 9, "8d8d81"),
        ];
        for (packed, length, seq) in cases {
            let b = seq_data_to_binary(&SeqData::Ncbi2na(packed.to_vec()), length).unwrap();
            assert_eq!(hex(&b.seq), seq);
            assert!(b.ambig.is_empty());
            assert_eq!(find_sequence_length(false, &b.seq), length);
        }
    }

    #[test]
    fn runs_use_the_old_format_up_to_15_and_split_at_4095() {
        let run = |n: usize| {
            let mut letters = vec![b'A'];
            letters.extend(std::iter::repeat_n(b'N', n));
            letters.push(b'C');
            iupacna(&letters).ambig
        };
        // 15 N: old format, one word: value f, length - 1 = e, offset 1.
        assert_eq!(hex(&run(15)), "00000001fe000001");
        // 16 N: new format, two words.
        assert_eq!(hex(&run(16)), "80000002f00f000000000001");
        // 4095 N: one region; 4096 N: two regions (4095 + 1).
        assert_eq!(hex(&run(4095)), "80000002fffe000000000001");
        assert_eq!(hex(&run(4096)), "80000004fffe000000000001f000000000001000");
    }

    #[test]
    fn protein_letters_become_ncbistdaa() {
        let b = seq_data_to_binary(&SeqData::Ncbieaa(b"MK*V-L*".to_vec()), 7).unwrap();
        assert_eq!(hex(&b.seq), "0c0a1913000b19");
        let b = seq_data_to_binary(
            &SeqData::Ncbieaa(b"ABCDEFGHIJKLMNOPQRSTUVWXYZ".to_vec()),
            26,
        )
        .unwrap();
        assert_eq!(
            hex(&b.seq),
            "0102030405060708091b0a0b0c0d1a0e0f101112181314151617"
        );
        let b = seq_data_to_binary(&SeqData::Iupacaa(b"mkvLa".to_vec()), 5).unwrap();
        assert_eq!(hex(&b.seq), "0c0a130b01");
        assert!(b.ambig.is_empty());
    }

    #[test]
    fn other_seq_data_is_not_processable() {
        assert_eq!(
            seq_data_to_binary(&SeqData::Other, 10),
            Err(UnprocessableSequence)
        );
    }
}
