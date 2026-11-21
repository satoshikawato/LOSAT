//! Sequence data of a database volume: the conversions of
//! `c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp` from a
//! `CSeq_data` to the bytes of the sequence file, the ambiguity data of a
//! nucleotide sequence, and the generator (`CRandom`) that picks the base
//! stored in place of an ambiguous one.

use super::bioseq::SeqData;

// NCBI reference (598d8ae6): c++/include/util/random_gen.hpp:159-165
// NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:98
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
    // NCBI reference (598d8ae6): c++/src/util/random_gen.cpp:296-308
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
const LOG2: [i8; 16] = [-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1];

impl AmbigDataBuilder {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:158-160
    // The generator is seeded with the length of the sequence.
    fn new(size: u32) -> Self {
        Self {
            size,
            regions: Vec::new(),
            random: Random::new(size),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:186-213
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
fn pack_new_ambig(amb: &mut Vec<u8>, region: &AmbiguousRegion) {
    let length_m1 = region.length() - 1;
    amb.push((region.value << 4) | (length_m1 >> 8) as u8);
    amb.push((length_m1 & 0xFF) as u8);
    amb.push(0);
    amb.push(0);
    amb.extend_from_slice(&region.start.to_be_bytes());
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:292-336
fn pack_old_ambig(amb: &mut Vec<u8>, region: &AmbiguousRegion) {
    let length_m1 = region.length() - 1;
    let off = region.start;
    amb.push((region.value << 4) | length_m1 as u8);
    amb.push(((off >> 16) & 0xFF) as u8);
    amb.push(((off >> 8) & 0xFF) as u8);
    amb.push((off & 0xFF) as u8);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:434-441
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
// The table of 256 pairs (the letter as the first and as the second base of
// a byte) is 0xf0, 0xf (N) for every byte except the letters below, in
// either case (sequtil_tables.cpp:409-423, the rows of 'A'-'Z' and 'a'-'z'):
// A 1, B 14, C 2, D 13, G 4, H 11, K 12, M 3, R 5, S 6, T 8, U 8, V 7, W 9,
// Y 10.
#[cfg(test)]
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
#[cfg(test)]
const IUPACNA_TO_NCBI4NA: [u8; 256] = iupacna_to_ncbi4na_table();

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:324-339
/// `CSeqConvert::Convert` from iupacna to ncbi4na: two bases per byte, the
/// first in the high half.
#[cfg(test)]
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

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1674-1707
// The table is 0xff except the 26 letters of either case, each the code of its capital in
// `STDAA_OF_LETTER`.
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:530-567
/// The NCBIstdaa codes of a protein's `Seq-data` as CWriteDB writes them
/// (`WriteDB_StdaaToBinary`, `WriteDB_EaaToBinary`, `WriteDB_IupacaaToBinary`: every letter of
/// the data, whatever `Seq-inst.length` says), `None` for another coding.
pub(crate) fn protein_codes(seq_data: &SeqData) -> Option<Vec<u8>> {
    match seq_data {
        SeqData::Ncbieaa(v) => Some(convert_1_to_1(v, &NCBIEAA_TO_NCBISTDAA)),
        SeqData::Iupacaa(v) => Some(convert_1_to_1(v, &IUPACAA_TO_NCBISTDAA)),
        SeqData::Ncbistdaa(v) => Some(v.clone()),
        SeqData::Ncbi2na(_) | SeqData::Ncbi4na(_) | SeqData::Other => None,
    }
}

fn convert_1_to_1(src: &[u8], table: &[u8; 256]) -> Vec<u8> {
    src.iter().map(|&byte| table[usize::from(byte)]).collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbobj.cpp:47-60
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbobj.cpp:100-111
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1709-1712
// `eCoding_Ncbi` gives a protein's NCBIstdaa codes (the bytes of the sequence file) and a
// nucleotide's NCBI4na code per base (after `s_FixBioseqDeltas`, whose gaps are 0x0f).
/// `SeqDB_SequenceHash(bioseq)`: the hash of `-hash_index` over `length` residues of
/// `seq_data`, `None` for a coding the writer cannot write (it stops before the hash is used).
pub(crate) fn sequence_hash(seq_data: &SeqData, length: u64) -> Option<u32> {
    let length = length as usize;
    Some(match seq_data {
        SeqData::Ncbieaa(v) => hash(
            v.iter()
                .take(length)
                .map(|&b| NCBIEAA_TO_NCBISTDAA[usize::from(b)]),
        ),
        SeqData::Iupacaa(v) => hash(
            v.iter()
                .take(length)
                .map(|&b| IUPACAA_TO_NCBISTDAA[usize::from(b)]),
        ),
        SeqData::Ncbistdaa(v) => hash(v.iter().take(length).copied()),
        SeqData::Ncbi2na(v) => {
            hash((0..length).map(|i| 1u8 << ((v[i / 4] >> (6 - 2 * (i % 4))) & 3)))
        }
        SeqData::Ncbi4na(v) => hash((0..length).map(|i| (v[i / 2] >> (4 - 4 * (i % 2))) & 0x0f)),
        SeqData::Other => return None,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbobj.cpp:47-60
/// `SeqDB_SequenceHash(sequence, length)` of the codes of a sequence.
fn hash(codes: impl Iterator<Item = u8>) -> u32 {
    codes.fold(0u32, |h, code| {
        h.wrapping_mul(1_103_515_245)
            .wrapping_add(u32::from(code).wrapping_add(12345))
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1690-1700
/// `x_ComputeHash(sequence, ambiguities)`: the hash of `-hash_index` of a sequence given in the
/// form of the sequence file (`AddSequence(seq, ambig)`): a protein's NCBIstdaa bytes, a
/// nucleotide's bases as NCBI4na codes (`SeqDB_UnpackAmbiguities`). The error is the message of
/// a malformed nucleotide.
pub(crate) fn raw_sequence_hash(
    protein: bool,
    seq: &[u8],
    ambig: &[u8],
) -> Result<u32, &'static str> {
    if protein {
        return Ok(hash(seq.iter().copied()));
    }
    let na8 = crate::seqdb::volume::unpack_ambiguities(seq, ambig)?;
    Ok(hash(na8.into_iter()))
}

/// The `Seq-data` choices CWriteDB cannot write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnprocessableSequence;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1007-1036
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:530-567
// LOSAT has no counterpart of `WriteDB_StdaaToBinary` and `WriteDB_IupacaaToBinary` (nor of the
// `e_Ncbistdaa` and `e_Iupacaa` cases of the switch above): LOSAT's FASTA reader never makes these
// Seq-data choices, so they are not ported (P-V3b).
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:592-609
// LOSAT has no counterpart of `WriteDB_IupacnaToBinary` (nor of the `e_Iupacna` case of the switch
// above): LOSAT's FASTA reader never makes this Seq-data choice, so it is not ported (P-V3b).
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
        SeqData::Ncbieaa(_) | SeqData::Iupacaa(_) | SeqData::Ncbistdaa(_) => {
            protein(protein_codes(seq_data).expect("a protein coding"))
        }
        SeqData::Ncbi2na(v) => ncbi2na_to_binary(v, length),
        SeqData::Ncbi4na(v) => ncbi4na_to_binary(v, length),
        SeqData::Other => return Err(UnprocessableSequence),
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_general.cpp:80-88
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
#[path = "convert_tests.rs"]
mod tests;
