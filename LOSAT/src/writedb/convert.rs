//! Sequence data of a database volume: the conversions of
//! `c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp` from a
//! `CSeq_data` to the bytes of the sequence file and the ambiguity data of a
//! nucleotide sequence (using `utils::random_gen::NcbiRandom` to pick the base
//! stored in place of an ambiguous one).

use super::bioseq::SeqData;
use super::fasta_source::pack::{iupacna_to_ncbi4na, pair_byte};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::utils::random_gen::NcbiRandom;

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
    random: NcbiRandom,
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
            random: NcbiRandom::new(size),
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

/// The number of input bytes (64 bases) that `ncbi4na_to_binary` converts at a time
/// ([`ncbi4na_blocks_to_binary`]; LOSAT's speed).
const NCBI4NA_BLOCK: usize = 32;

/// 0 for an NCBI4na base (a value below 16) that `NA4_TO_NA2` maps to a code (`A`, `C`, `G`,
/// `T`: 1, 2, 4, 8, the powers of two), non-zero for every other base, without a table: a power
/// of two has no bit in common with the number one below it, and the bit of 16 added to the base
/// keeps the result for 0 (whose number one below wraps to 255) non-zero.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:434-441
#[inline(always)]
const fn ambiguity_bits(base: u8) -> u8 {
    base.wrapping_sub(1) & (base | 0x10)
}

/// The NCBI2na code (`NA4_TO_NA2`) of an NCBI4na base for which [`ambiguity_bits`] is 0, without
/// a table: 1, 2, 4, 8 give 0, 1, 2, 3 (half the base less an eighth of it, in whole numbers).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:434-441
#[inline(always)]
const fn ncbi2na_code_of_base(base: u8) -> u8 {
    (base >> 1) - (base >> 3)
}

/// The body of `WriteDB_Ncbi4naToBinary`'s loop for the input byte `i` (`inp`): the half of its
/// two bases, by `NA4_TO_NA2` where both are `A`, `C`, `G` or `T`, otherwise by
/// [`ncbi4na_ambiguous_half`], ORed into the output byte `i / 2`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:444-518
#[inline(always)]
fn ncbi4na_byte_to_binary(
    seq: &mut [u8],
    ambiguities: &mut AmbigDataBuilder,
    i: usize,
    inp: u8,
    base_length: u32,
) {
    let c1 = NA4_TO_NA2[usize::from(inp >> 4)];
    let c2 = NA4_TO_NA2[usize::from(inp & 0xF)];
    let half = if (c1 | c2) & 0x80 == 0 {
        (c1 << 2) | c2
    } else {
        ncbi4na_ambiguous_half(ambiguities, i, inp, base_length)
    };
    or_half(seq, i, half);
}

/// The arm of `WriteDB_Ncbi4naToBinary`'s loop for an input byte `i` (`inp`) with a base that is
/// not `A`, `C`, `G` or `T`: a gap (0) becomes 15, but not the half byte past the end of an odd
/// count of bases, and both bases go through `check` (the ambiguity data, and the random base
/// of an ambiguous one), whose two codes are the half returned.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:444-518
#[inline(always)]
fn ncbi4na_ambiguous_half(
    ambiguities: &mut AmbigDataBuilder,
    i: usize,
    inp: u8,
    base_length: u32,
) -> u8 {
    let mut b1 = inp >> 4;
    let mut b2 = inp & 0xF;
    if b1 == 0 {
        b1 = 0xF;
    }
    let second = (i as u32) * 2 + 1;
    if b2 == 0 && second < base_length {
        b2 = 0xF;
    }
    (ambiguities.check(b1, second - 1) << 2) | ambiguities.check(b2, second)
}

/// ORs the half of the input byte `i` into its output byte `i / 2`: an even byte's in the high
/// four bits, an odd byte's below them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:444-518
#[inline(always)]
fn or_half(seq: &mut [u8], i: usize, half: u8) {
    seq[i / 2] |= if i & 1 == 1 { half } else { half << 4 };
}

/// `WriteDB_Ncbi4naToBinary`'s loop over the input bytes `read` ([`NCBI4NA_BLOCK`] or more), a
/// block of [`NCBI4NA_BLOCK`] bytes at a time (LOSAT's speed). Every block takes one pass of
/// arithmetic: the half of each of its bytes whose two bases are `A`, `C`, `G` or `T`
/// ([`ncbi2na_code_of_base`]), and a flag for each of its other bytes. A block without a flag
/// writes its output bytes from these halves ([`pair_byte`]) without a branch a byte. A block
/// with a flag writes each output byte once, from its two input bytes in order: a flagged byte's
/// half from [`ncbi4na_ambiguous_half`], the other's from the pass. (Writing a block's output
/// bytes by vector and then ORing the flagged halves into them made all-N input run above dev:
/// each OR read a byte just stored by vector, and the processor's memory-ordering machine
/// clears rose to several times dev's.) NCBI's loop calls `check` for exactly the flagged bytes
/// and in that order, so the ambiguity data and the random bases are NCBI's.
/// The bytes after the last block run the loop's body ([`ncbi4na_byte_to_binary`]).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:444-518
#[inline(never)]
fn ncbi4na_blocks_to_binary(
    read: &[u8],
    seq: &mut [u8],
    ambiguities: &mut AmbigDataBuilder,
    base_length: u32,
) {
    let (blocks, _) = read.as_chunks::<NCBI4NA_BLOCK>();
    let mut halves = [0u8; NCBI4NA_BLOCK];
    let mut ambiguous = [0u8; NCBI4NA_BLOCK];
    for (index, block) in blocks.iter().enumerate() {
        let start = index * NCBI4NA_BLOCK;
        let mut any = 0u8;
        for ((half, flag), &inp) in halves.iter_mut().zip(ambiguous.iter_mut()).zip(block) {
            let bits = ambiguity_bits(inp >> 4) | ambiguity_bits(inp & 0xF);
            *flag = bits;
            any |= bits;
            let codes = (ncbi2na_code_of_base(inp >> 4) << 2) | ncbi2na_code_of_base(inp & 0xF);
            *half = if bits == 0 { codes } else { 0 };
        }
        // The half of the even byte in the high four bits, that of the odd byte below it.
        let out = &mut seq[start / 2..(start + NCBI4NA_BLOCK) / 2];
        if any == 0 {
            for (byte, &two) in out.iter_mut().zip(halves.as_chunks::<2>().0) {
                *byte = pair_byte(two);
            }
        } else {
            let pairs = halves
                .as_chunks::<2>()
                .0
                .iter()
                .zip(ambiguous.as_chunks::<2>().0)
                .zip(block.as_chunks::<2>().0);
            for (k, (byte, ((&[high, low], &[high_flag, low_flag]), &[high_inp, low_inp]))) in
                out.iter_mut().zip(pairs).enumerate()
            {
                let i = start + 2 * k;
                let high = if high_flag != 0 {
                    ncbi4na_ambiguous_half(ambiguities, i, high_inp, base_length)
                } else {
                    high
                };
                let low = if low_flag != 0 {
                    ncbi4na_ambiguous_half(ambiguities, i + 1, low_inp, base_length)
                } else {
                    low
                };
                *byte = pair_byte([high, low]);
            }
        }
    }
    for (i, &inp) in read.iter().enumerate().skip(blocks.len() * NCBI4NA_BLOCK) {
        ncbi4na_byte_to_binary(seq, ambiguities, i, inp, base_length);
    }
}

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

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:444-518
// `byte_length` is the size of the data and `ncbi4na` the bytes that NCBI's pointer reaches (for
// iupacna, its NCBI4na string and the string's terminating zero). A size other than the
// `ceil(base_length / 2)` bytes that the loop reads gets a line on the standard output (`cout`,
// before the debug-only assertion); the loop then reads the bytes it reads, and fewer than that is
// a read past the data, which makeprofiledb refuses (`short_nucleotide_data`; here the loop stops
// at the end of `ncbi4na`). The Seq-data of the FASTA source always has its bytes.
/// `WriteDB_Ncbi4naToBinary`, with its line on the standard output appended to `cout`.
pub(crate) fn ncbi4na_to_binary(
    ncbi4na: &[u8],
    byte_length: usize,
    base_length: u32,
    cout: &mut Vec<u8>,
) -> BinarySequence {
    let inp_bytes = base_length.div_ceil(2) as usize;
    let blast_bytes = (base_length / 4 + 1) as usize;
    let remainder = (base_length & 3) as u8;
    let last_byte = blast_bytes - 1;
    let mut ambiguities = AmbigDataBuilder::new(base_length);
    if inp_bytes != byte_length {
        cout.extend_from_slice(format!("ib={inp_bytes},n4sz={byte_length}\n").as_bytes());
    }
    let mut seq = vec![0u8; blast_bytes];
    let read = &ncbi4na[..inp_bytes.min(ncbi4na.len())];
    // A sequence shorter than a block runs the loop's body for each byte (the blocks' loop would
    // take none of its bytes).
    if read.len() < NCBI4NA_BLOCK {
        for (i, &inp) in read.iter().enumerate() {
            ncbi4na_byte_to_binary(&mut seq, &mut ambiguities, i, inp, base_length);
        }
    } else {
        ncbi4na_blocks_to_binary(read, &mut seq, &mut ambiguities, base_length);
    }
    seq[last_byte] &= 255 - 3;
    seq[last_byte] |= remainder;
    BinarySequence {
        seq,
        ambig: ambiguities.get_ambig(),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:569-588
// The first `ceil(base_length / 4)` bytes of the data are copied (bytes past them are not read);
// fewer bytes is a read past the data, which makeprofiledb refuses (`short_nucleotide_data`).
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
        SeqData::Ncbieaa(v) => Some(ncbieaa_codes(v)),
        SeqData::Iupacaa(v) => Some(iupacaa_codes(v)),
        SeqData::Ncbistdaa(v) => Some(v.clone()),
        SeqData::Ncbi2na(_) | SeqData::Ncbi4na(_) | SeqData::Iupacna(_) | SeqData::Other => None,
    }
}

/// `CSeqConvert` from ncbieaa to ncbistdaa (`CEaaToStdaa`): 255 for a byte that is no residue.
pub(crate) fn ncbieaa_codes(letters: &[u8]) -> Vec<u8> {
    convert_1_to_1(letters, &NCBIEAA_TO_NCBISTDAA)
}

/// `CSeqConvert` from iupacaa to ncbistdaa (`CIupacaaToStdaa`): 255 for a byte that is no
/// letter.
pub(crate) fn iupacaa_codes(letters: &[u8]) -> Vec<u8> {
    convert_1_to_1(letters, &IUPACAA_TO_NCBISTDAA)
}

fn convert_1_to_1(src: &[u8], table: &[u8; 256]) -> Vec<u8> {
    src.iter().map(|&byte| table[usize::from(byte)]).collect()
}

/// `CSeqConvert`'s exception for a conversion that no function covers (`CSeqUtilException`,
/// `eInvalidCoding`): one from a nucleotide coding to NCBIstdaa (psiblast's query from an
/// alignment, `CPsiBlastInputClustalW`; makeprofiledb's query of a profile with frequency ratios
/// only, `CPssm::GetQuerySequenceData`). NCBI BLAST+ 2.17.0 reports it with `Unknown
/// conversion.` (PS6 corpus `e.len1`, `n.q`; row PS10d survey `FA_*`).
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:225
pub(crate) fn unknown_conversion() -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: "src/util/sequtil/sequtil_convert_imp.cpp",
            line: 225,
            kind: "CSeqUtilException::Invalid coding",
            function: Some("ncbi::CSeqConvert_imp::Convert"),
        },
        "Unknown conversion.",
    )
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:6123-6153,6436-6467
// `CSeqVectorTypes::sx_GetConvertTable` builds the IUPACna -> NCBI8na table from `CSeqportUtil`,
// whose `iupacna` alphabet covers ASCII 65..=89 (`A`..=`Y`, uppercase only; `U` and every other
// byte, lowercase included, stay `0xff`).
const fn seq_vector_iupacna_to_ncbi8na_table() -> [u8; 256] {
    let mut table = [0xFFu8; 256];
    table[b'A' as usize] = 1;
    table[b'B' as usize] = 14;
    table[b'C' as usize] = 2;
    table[b'D' as usize] = 13;
    table[b'G' as usize] = 4;
    table[b'H' as usize] = 11;
    table[b'K' as usize] = 12;
    table[b'M' as usize] = 3;
    table[b'N' as usize] = 15;
    table[b'R' as usize] = 5;
    table[b'S' as usize] = 6;
    table[b'T' as usize] = 8;
    table[b'V' as usize] = 7;
    table[b'W' as usize] = 9;
    table[b'Y' as usize] = 10;
    table
}
const SEQ_VECTOR_IUPACNA_TO_NCBI8NA: [u8; 256] = seq_vector_iupacna_to_ncbi8na_table();

// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:6182-6215,6622-6653
// `CSeqportUtil`'s `iupacaa` alphabet covers uppercase ASCII 65..=90 only; lowercase stays `0xff`.
const fn seq_vector_iupacaa_to_ncbistdaa_table() -> [u8; 256] {
    let mut table = [0xFFu8; 256];
    let mut i = 0;
    while i < 26 {
        table[b'A' as usize + i] = STDAA_OF_LETTER[i];
        i += 1;
    }
    table
}
const SEQ_VECTOR_IUPACAA_TO_NCBISTDAA: [u8; 256] = seq_vector_iupacaa_to_ncbistdaa_table();

/// Why `SeqDB_SequenceHash` (`CSeqVector`) cannot hash a Bioseq's `Seq-data` for `length` letters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum SequenceHashError {
    /// `CSeqVector_CI::ThrowOutOfRangeSeq_inst` (`CSeqVectorException::eOutOfRange`).
    OutOfRange(NcbiException),
    /// `copy_2bit_any` / `copy_4bit_any` compare the truncating quotient of the block end against
    /// the byte count, so 1..3 bases (2na) or 1 base (4na) past the data reads one byte past the
    /// vector (undefined behaviour in NCBI).
    PastEnd(&'static str),
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:56-62
// NCBI reference (598d8ae6): c++/include/objmgr/impl/seq_vector_cvt.hpp:47-131
/// `CSeqVector_CI`'s cache check in blocks of 1024 residues before hashing `length` residues of
/// a sequence buffer of `bytes` bytes holding `per_byte` residues per byte (1, 2 or 4).
fn check_seq_vector_bounds(
    coding: &'static str,
    bytes: usize,
    per_byte: usize,
    length: usize,
) -> Result<(), SequenceHashError> {
    const CACHE_SIZE: usize = 1024;
    let mut pos = 0usize;
    while pos < length {
        let end = (pos + CACHE_SIZE).min(length);
        if end / per_byte > bytes {
            return Err(SequenceHashError::OutOfRange(NcbiException::new(
                ThrowSite {
                    file: "src/objmgr/seq_vector_ci.cpp",
                    line: 60,
                    kind: "CSeqVectorException::eOutOfRange",
                    function: Some("ncbi::objects::ThrowOutOfRangeSeq_inst"),
                },
                format!("reference out of range of Seq-inst data: {end}"),
            )));
        }
        if end.div_ceil(per_byte) > bytes {
            return Err(SequenceHashError::PastEnd(coding));
        }
        pos = end;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbobj.cpp:47-60
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbobj.cpp:100-111
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1709-1712
// `eCoding_Ncbi` gives a protein's NCBIstdaa codes and a nucleotide's NCBI8na codes through
// `CSeqVector`'s `CSeqportUtil` tables (after `s_FixBioseqDeltas`, whose gaps are 0x0f).
/// `SeqDB_SequenceHash` of a bioseq: the hash of `-hash_index` over `length` residues of
/// `seq_data`, `Ok(None)` for a coding the writer cannot write (it stops before the hash is used),
/// or `Err` when `CSeqVector` throws `eOutOfRange` or reads past the packed bytes.
pub(crate) fn sequence_hash(
    seq_data: &SeqData,
    length: u64,
) -> Result<Option<u32>, SequenceHashError> {
    let length = length as usize;
    let value = match seq_data {
        SeqData::Ncbieaa(v) => {
            check_seq_vector_bounds("ncbieaa", v.len(), 1, length)?;
            hash(
                v.iter()
                    .take(length)
                    .map(|&b| NCBIEAA_TO_NCBISTDAA[usize::from(b)]),
            )
        }
        SeqData::Iupacaa(v) => {
            check_seq_vector_bounds("iupacaa", v.len(), 1, length)?;
            hash(
                v.iter()
                    .take(length)
                    .map(|&b| SEQ_VECTOR_IUPACAA_TO_NCBISTDAA[usize::from(b)]),
            )
        }
        SeqData::Ncbistdaa(v) => {
            check_seq_vector_bounds("ncbistdaa", v.len(), 1, length)?;
            hash(v.iter().take(length).copied())
        }
        SeqData::Iupacna(v) => {
            check_seq_vector_bounds("iupacna", v.len(), 1, length)?;
            hash(
                v.iter()
                    .take(length)
                    .map(|&b| SEQ_VECTOR_IUPACNA_TO_NCBI8NA[usize::from(b)]),
            )
        }
        SeqData::Ncbi2na(v) => {
            check_seq_vector_bounds("ncbi2na", v.len(), 4, length)?;
            hash((0..length).map(|i| 1u8 << ((v[i / 4] >> (6 - 2 * (i % 4))) & 3)))
        }
        SeqData::Ncbi4na(v) => {
            check_seq_vector_bounds("ncbi4na", v.len(), 2, length)?;
            hash((0..length).map(|i| (v[i / 2] >> (4 - 4 * (i % 2))) & 0x0f))
        }
        SeqData::Other => return Ok(None),
    };
    Ok(Some(value))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbobj.cpp:47-60
/// `SeqDB_SequenceHash` of a sequence and its length, over the codes of a sequence.
fn hash(codes: impl Iterator<Item = u8>) -> u32 {
    codes.fold(0u32, |h, code| {
        h.wrapping_mul(1_103_515_245)
            .wrapping_add(u32::from(code).wrapping_add(12345))
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1690-1700
/// `x_ComputeHash` of a sequence and its ambiguities: the hash of `-hash_index` of a sequence given
/// in the form of the sequence file (the arguments of `AddSequence`, `seq` and `ambig`): a
/// protein's NCBIstdaa bytes, a nucleotide's bases as NCBI4na codes (`SeqDB_UnpackAmbiguities`).
/// The error is the message of a malformed nucleotide.
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
// The six codings that CWriteDB converts, whatever the volume (a protein volume writes the bytes
// of a nucleotide coding, without its ambiguity data: `volume::write_sequence`); the others stop
// the publication.
/// `CWriteDB_Impl::x_CookSequence` for a `CBioseq` with `Seq-data`: the sequence bytes and the
/// ambiguity data of a sequence of `length` letters, and the line that the conversion of NCBI4na
/// writes on the standard output, appended to `cout`.
pub(crate) fn seq_data_to_binary(
    seq_data: &SeqData,
    length: u32,
    cout: &mut Vec<u8>,
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
        SeqData::Ncbi4na(v) => ncbi4na_to_binary(v, v.len(), length, cout),
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:590-610
        // Every letter converts (none for empty data: `CSeqConvert`'s template returns before
        // the conversion), and the conversion of the length reads the NCBI4na string through
        // its `c_str()`: its bytes, then the terminating zero.
        SeqData::Iupacna(v) => {
            let mut na4 = iupacna_to_ncbi4na(v);
            let byte_length = na4.len();
            na4.push(0);
            ncbi4na_to_binary(&na4, byte_length, length, cout)
        }
        SeqData::Other => return Err(UnprocessableSequence),
    })
}

/// A nucleotide `Seq-data` with fewer bytes than `x_CookSequence`'s conversion of its length
/// reads: its coding, how much it holds and the bytes that NCBI reads (of iupacna: its letters,
/// and the bytes of their NCBI4na form).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ShortNucleotideData {
    pub(crate) coding: &'static str,
    /// Whether the data is text (iupacna), whose `held` counts letters, not bytes.
    pub(crate) text: bool,
    pub(crate) held: usize,
    pub(crate) read: usize,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:569-588
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:444-528
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:590-610
// The conversions read `ceil(length / 4)` bytes of ncbi2na and `ceil(length / 2)` bytes of
// ncbi4na from the data (its first byte's address even when it is empty), and of iupacna the
// `ceil(letters / 2)` bytes of its NCBI4na string and the string's terminating zero; the assertions
// of the right size are debug-only. NCBI BLAST+ 2.17.0 reads past the data there (valgrind's
// invalid reads and uninitialised values, bytes that change between runs, SIGSEGV for empty
// ncbi2na and ncbi4na; row PS10d survey §6). One byte short for iupacna reads the zero (its two
// bases are N), which is NCBI's defined result.
/// The data of `length` letters that CWriteDB's conversion reads past; `None` for data of enough
/// bytes and for another coding.
pub(crate) fn short_nucleotide_data(
    seq_data: &SeqData,
    length: u32,
) -> Option<ShortNucleotideData> {
    let (coding, text, held, available, read) = match seq_data {
        SeqData::Ncbi2na(v) => ("ncbi2na", false, v.len(), v.len(), length.div_ceil(4)),
        SeqData::Ncbi4na(v) => ("ncbi4na", false, v.len(), v.len(), length.div_ceil(2)),
        SeqData::Iupacna(v) => (
            "iupacna",
            true,
            v.len(),
            v.len().div_ceil(2) + 1,
            length.div_ceil(2),
        ),
        _ => return None,
    };
    let read = read as usize;
    (available < read).then_some(ShortNucleotideData {
        coding,
        text,
        held,
        read,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:590-604
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:313-339
// `WriteDB_IupacnaToBinary` converts all `v.size()` bytes of the IUPACna string into NCBI4na. In
// `x_ConvertIupacnaTo4na`, the pair loop indexes `sm_IupacnaTo4na` with a signed `char` (`*src_i`),
// so any byte >= 0x80 in positions `0..(v.len() & !1)` reads before the static table (undefined
// behaviour); the single odd overhang byte casts to `Uint1` and is defined (`0xf0`).
/// The first byte >= 0x80 (with its 0-based position) in a pair position (`0..(v.len() & !1)`)
/// of an IUPACna `Seq-data`; `None` for any other coding or when every pair byte is below 0x80.
pub(crate) fn iupacna_high_byte(seq_data: &SeqData) -> Option<(usize, u8)> {
    let SeqData::Iupacna(v) = seq_data else {
        return None;
    };
    let pairs_len = v.len() & !1;
    v[..pairs_len]
        .iter()
        .enumerate()
        .find(|&(_, &byte)| byte >= 0x80)
        .map(|(index, &byte)| (index, byte))
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
