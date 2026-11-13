//! The reading frames of a nucleotide query (`READINGFRAME`), its six translations
//! (`TranslateAll` through the toolkit's `CSeqTranslator::Translate` with the standard genetic
//! code) and its reverse complement (`GetComplementSeq`).
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:167-176
//! ```c++
//! struct READINGFRAME
//! {
//! 	typedef int TFRAMEID;
//! 	typedef unsigned int TFRAMEINDEX;
//!
//! 	static constexpr const TFRAMEINDEX RF_SIZE = 3;	//reading frame size
//! 	static constexpr const TFRAMEINDEX TOTAL_RFS = RF_SIZE + RF_SIZE;	//total reading frames
//! ```
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:28-29
//! ```c++
//! const char * READINGFRAME::RF_TITLES[] = {"RF +1", "RF +2", "RF +3", "RF -1", "RF -2", "RF -3"};
//! READINGFRAME::TFRAMEID READINGFRAME::RF_IDS[] = {1, 2, 3, -1, -2, -3};
//! ```
//! The positions are `SeqPos_t` (`int`) and the lengths `SeqLen_t` (`unsigned int`); the
//! functions below follow C++'s conversions: an expression with an `unsigned int` operand is
//! computed modulo 2^32 and stored back in an `int` (two's complement in NCBI's build).

/// `READINGFRAME::RF_SIZE`.
pub(crate) const RF_SIZE: u32 = 3;

/// `READINGFRAME::RF_IDS`: the frame of each index (`Idx2Id`).
pub(crate) const RF_IDS: [i32; 6] = [1, 2, 3, -1, -2, -3];

/// `PlusRFIdx2MinusRFIdx`: the minus-strand frame index of a plus-strand one.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:200
// ```c++
// 	static TFRAMEINDEX PlusRFIdx2MinusRFIdx(TFRAMEINDEX prfidx, SeqLen_t na_len) {return (na_len - prfidx) % READINGFRAME::RF_SIZE + READINGFRAME::RF_SIZE;}
// ```
pub(crate) fn plus_rf_idx_to_minus(prfidx: u32, na_len: u32) -> u32 {
    na_len.wrapping_sub(prfidx) % RF_SIZE + RF_SIZE
}

/// `PlusPr2PlusNA`: the plus-strand start of the codon of residue `pr` of a plus-strand frame.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:224-227
// ```c++
//     static SeqPos_t PlusPr2PlusNA(SeqPos_t pr, TFRAMEINDEX rfidx, SeqLen_t na_len)
// 	{
// 		return pr * READINGFRAME::RF_SIZE + rfidx;
// 	}
// ```
pub(crate) fn plus_pr_to_plus_na(pr: i32, rfidx: u32) -> i32 {
    (pr as u32).wrapping_mul(RF_SIZE).wrapping_add(rfidx) as i32
}

/// `MinusPr2PlusNA`: the plus-strand position of the first base of the codon of residue `pr`
/// of a minus-strand frame. `None` where `pr - 1` overflows NCBI's `int` (`pr` = `INT_MIN`).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:230-237
// ```c++
// 	static SeqPos_t MinusPr2PlusNA(SeqPos_t pr, TFRAMEINDEX rfidx, SeqLen_t na_len)
// 	{
// 		// = na_len - (pr * READINGFRAME::RF_SIZE + (rfidx - READINGFRAME::RF_SIZE)) - 1
// 		// = na_len - (pr * READINGFRAME::RF_SIZE + rfidx - READINGFRAME::RF_SIZE) - 1
// 		// = na_len - (pr - 1) * READINGFRAME::RF_SIZE - rfidx - 1
// 		return na_len - (pr - 1) * READINGFRAME::RF_SIZE - rfidx - 1;
// 	}
// ```
pub(crate) fn minus_pr_to_plus_na(pr: i32, rfidx: u32, na_len: u32) -> Option<i32> {
    let pr_1 = pr.checked_sub(1)?;
    Some(
        na_len
            .wrapping_sub((pr_1 as u32).wrapping_mul(RF_SIZE))
            .wrapping_sub(rfidx)
            .wrapping_sub(1) as i32,
    )
}

/// The `ncbieaa` string of genetic code 1 (`gc.prt`, "Standard"), in TCAG order.
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/gc.inc:40-44
// ```c++
//     "  name \"Standard\" ,",
//     "  name \"SGC0\" ,",
//     "  id 1 ,",
//     "  ncbieaa  \"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\",",
//     "  sncbieaa \"---M------**--*----M---------------M----------------------------\"",
// ```
const NCBIEAA_1: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// `CTrans_table::sm_BaseToIdx`: the IUPACna letters (either case, `U` as `T`, `X` as `N`) and
/// the NCBI4na codes to their 4-bit base sets; every other byte to 0 (a gap).
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:90-117
// ```c++
// void CTrans_table::x_InitFsaTable (void)
// {
//     char         ch;
//     int          i, j, k, p, q, r, nx, st;
//     static char  charToBase [17] = "-ACMGRSVTWYHKDBN";
//     static char  baseToComp [17] = "-TGKCYSBAWRDMHVN";
//
//     // illegal characters map to 0
//     for (i = 0; i < 256; i++) {
//         sm_BaseToIdx [i] = 0;
//     }
//
//     // map iupacna alphabet to EBaseCode
//     for (i = eBase_gap; i <= eBase_N; i++) {
//         ch = charToBase [i];
//         sm_BaseToIdx [(int) ch] = i;
//         ch = (unsigned char)tolower (ch);
//         sm_BaseToIdx [(int) ch] = i;
//     }
//     sm_BaseToIdx [(int) 'U'] = eBase_T;
//     sm_BaseToIdx [(int) 'u'] = eBase_T;
//     sm_BaseToIdx [(int) 'X'] = eBase_N;
//     sm_BaseToIdx [(int) 'x'] = eBase_N;
//
//     // also map ncbi4na alphabet to EBaseCode
//     for (i = eBase_gap; i <= eBase_N; i++) {
//         sm_BaseToIdx [(int) i] = i;
//     }
// ```
const fn base_to_idx_table() -> [u8; 256] {
    const CHAR_TO_BASE: &[u8; 16] = b"-ACMGRSVTWYHKDBN";
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 16 {
        let ch = CHAR_TO_BASE[i];
        table[ch as usize] = i as u8;
        table[ch.to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    table[b'U' as usize] = 8;
    table[b'u' as usize] = 8;
    table[b'X' as usize] = 15;
    table[b'x' as usize] = 15;
    let mut i = 0;
    while i < 16 {
        table[i] = i as u8;
        i += 1;
    }
    table
}

const BASE_TO_IDX: [u8; 256] = base_to_idx_table();

/// `CTrans_table::m_AminoAcid` of genetic code 1, indexed by the three 4-bit base sets of a
/// codon (`256 * i + 16 * j + k`; NCBI's state is that plus 1).
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:139-232
// ```c++
//     static int  expansions [4] = {eBase_A, eBase_C, eBase_G, eBase_T};
//                                 // T = 0, C = 1, A = 2, G = 3
//     static int  codonIdx [9] = {0, 2, 1, 0, 3, 0, 0, 0, 0};
//     ...
//     // ambiguous codons map to unknown amino acid or not start
//     for (i = 0; i <= 4096; i++) {
//         m_AminoAcid [i] = 'X';
//         ...
//     }
//
//     // lookup amino acid for each codon in genetic code table
//     for (i = eBase_gap, st = 1; i <= eBase_N; i++) {
//         for (j = eBase_gap; j <= eBase_N; j++) {
//             for (k = eBase_gap; k <= eBase_N; k++, st++) {
//                 aa = '\0';
//                 ...
//                 // expand ambiguous IJK nucleotide symbols into component bases XYZ
//                 for (p = 0; p < 4 && go_on; p++) {
//                     x = expansions [p];
//                     if ((x & i) != 0) {
//                         for (q = 0; q < 4 && go_on; q++) {
//                             y = expansions [q];
//                             if ((y & j) != 0) {
//                                 for (r = 0; r < 4 && go_on; r++) {
//                                     z = expansions [r];
//                                     if ((z & k) != 0) {
//                                         ...
//                                         cd = 16 * codonIdx [x] + 4 * codonIdx [y] + codonIdx [z];
//
//                                         // lookup amino acid for codon XYZ
//                                         ch = (*ncbieaa) [cd];
//                                         if (aa == '\0') {
//                                             aa = ch;
//                                         } else if (aa != ch) {
//                                             // allow Asx (Asp or Asn) and Glx (Glu or Gln)
//                                             if ((aa == 'B' || aa == 'D' || aa == 'N') &&
//                                                 (ch == 'D' || ch == 'N')) {
//                                                 aa = 'B';
//                                             } else if ((aa == 'Z' || aa == 'E' || aa == 'Q') &&
//                                                        (ch == 'E' || ch == 'Q')) {
//                                                 aa = 'Z';
//                                             } else if ((aa == 'J' || aa == 'I' || aa == 'L') &&
//                                                 (ch == 'I' || ch == 'L')) {
//                                                 aa = 'J';
//                                             } else {
//                                                 aa = 'X';
//                                             }
//                                         }
//                                         ...
//                                         // drop out of loop as soon as answer is known
//                                         if (aa == 'X' && orf == 'X') {
//                                             go_on = false;
//                                         }
//                 ...
//                 // assign amino acid
//                 if (aa != '\0') {
//                     m_AminoAcid [st] = aa;
//                 }
// ```
// (`X` is kept by every later merge, so stopping early changes no amino acid.)
const fn amino_acid_table() -> [u8; 4096] {
    const EXPANSIONS: [usize; 4] = [1, 2, 4, 8];
    const CODON_IDX: [usize; 9] = [0, 2, 1, 0, 3, 0, 0, 0, 0];
    let mut table = [b'X'; 4096];
    let mut i = 0;
    while i < 16 {
        let mut j = 0;
        while j < 16 {
            let mut k = 0;
            while k < 16 {
                let mut aa = 0u8;
                let mut p = 0;
                while p < 4 {
                    let x = EXPANSIONS[p];
                    if x & i != 0 {
                        let mut q = 0;
                        while q < 4 {
                            let y = EXPANSIONS[q];
                            if y & j != 0 {
                                let mut r = 0;
                                while r < 4 {
                                    let z = EXPANSIONS[r];
                                    if z & k != 0 {
                                        let cd =
                                            16 * CODON_IDX[x] + 4 * CODON_IDX[y] + CODON_IDX[z];
                                        let ch = NCBIEAA_1[cd];
                                        if aa == 0 {
                                            aa = ch;
                                        } else if aa != ch {
                                            aa = if matches!(aa, b'B' | b'D' | b'N')
                                                && matches!(ch, b'D' | b'N')
                                            {
                                                b'B'
                                            } else if matches!(aa, b'Z' | b'E' | b'Q')
                                                && matches!(ch, b'E' | b'Q')
                                            {
                                                b'Z'
                                            } else if matches!(aa, b'J' | b'I' | b'L')
                                                && matches!(ch, b'I' | b'L')
                                            {
                                                b'J'
                                            } else {
                                                b'X'
                                            };
                                        }
                                    }
                                    r += 1;
                                }
                            }
                            q += 1;
                        }
                    }
                    p += 1;
                }
                if aa != 0 {
                    table[256 * i + 16 * j + k] = aa;
                }
                k += 1;
            }
            j += 1;
        }
        i += 1;
    }
    table
}

const AMINO_ACIDS: [u8; 4096] = amino_acid_table();

/// The amino acid of the codon `b1 b2 b3` (bytes as in the sequence).
fn codon_residue(b1: u8, b2: u8, b3: u8) -> u8 {
    let (i, j, k) = (
        usize::from(BASE_TO_IDX[usize::from(b1)]),
        usize::from(BASE_TO_IDX[usize::from(b2)]),
        usize::from(BASE_TO_IDX[usize::from(b3)]),
    );
    AMINO_ACIDS[256 * i + 16 * j + k]
}

/// `CSeqTranslator::Translate(seq, prot, fIs5PrimePartial, code 1)` with rpsbproc's `*` to `X`:
/// one residue per codon, then a last partial codon (1 or 2 bases) padded with `N` and kept
/// only when it is not `X`.
// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:4108-4120
// ```c++
// void CSeqTranslator::Translate(const string& seq,
//                                string& prot,
//                                TTranslationFlags flags,
//                                const CGenetic_code* code,
//                                bool* alt_start)
// {
//     x_Translate(seq, prot, 0, code,
//                 !(flags & fIs5PrimePartial),
//                 !(flags & fIs3PrimePartial),
//                 !(flags & fNoStop),
//                 flags & fRemoveTrailingX,
//                 alt_start);
// }
// ```
// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3636-3725
// ```c++
//     // reserve our space
//     const size_t usable_size = seq.size() > frame ? seq.size() - frame : 0;
//     const size_t mod = usable_size % 3;
//     ...
//         size_t length = usable_size / 3;
//         bool check_start = (is_5prime_complete && frame == 0);
//         bool first_time = true;
//
//         for (i = 0; i < length; ++i) {
//
//             // loop through one codon at a time
//             for (k = 0; k < 3; ++k, ++start) {
//                 state = tbl.NextCodonState(state, *start);
//             }
//             ...
//             // save translated amino acid
//             if (first_time  &&  check_start) {
//                 aa = tbl.GetStartResidue(state);
//                 prot.append(1, aa);
//             } else {
//                 aa = tbl.GetCodonResidue(state);
//                 prot.append(1, aa);
//             }
//
//             first_time = false;
//         }
//
//         if (mod) {
//             for (k = 0; k < mod; ++k, ++start) {
//                 state = tbl.NextCodonState(state, *start);
//             }
//
//             for (; k < 3; ++k) {
//                 state = tbl.NextCodonState(state, 'N');
//             }
//             ...
//             // save translated amino acid
//             char c = tbl.GetCodonResidue(state);
//             if (first_time  &&  check_start) {
//                 aa = tbl.GetStartResidue(state);
//                 prot.append(1, aa);
//             } else if (c != 'X') {
//                 // if padding was needed, trim ambiguous last residue
//                 aa = tbl.GetCodonResidue(state);
//                 prot.append(1, aa);
//             }
//         }
//     ...
//     if ( aa != '*' && include_stop && (! mod) && prot.size() > 0 && is_3prime_complete ) {
//         // check for stop codon that normally encodes an amino acid
//         aa = tbl.GetStopResidue(state);
//         if (aa == '*') {
//             prot[prot.size()-1] = aa;
//         }
//     }
// ```
// `fIs5PrimePartial` makes `check_start` false. The state after three bases is that of the
// three bases alone (`sm_NextState`, Genetic_code_table.cpp:119-134: state 0 reads as `NN`).
// The stop check at the end changes nothing for code 1: `m_OrfStop` is `*` only for codons
// whose every expansion is a stop in `sncbieaa`, and those are stops in `ncbieaa` too, so the
// last residue is `*` already.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:130-133
// ```c++
// 			CSeqTranslator::Translate(trans, dimAaData[i], CSeqTranslator::fIs5PrimePartial/* | CSeqTranslator::fIs3PrimePartial*/, &gc);
// 			for (size_t j = 0, jEnd = dimAaData[i].size(); j < jEnd; ++j)
// 				if ('*' == dimAaData[i][j]) dimAaData[i][j] = 'X';
// ```
fn translate(seq: &[u8]) -> Vec<u8> {
    let mut protein = Vec::with_capacity(seq.len().div_ceil(3));
    let mut codons = seq.chunks_exact(3);
    for codon in codons.by_ref() {
        protein.push(codon_residue(codon[0], codon[1], codon[2]));
    }
    match *codons.remainder() {
        [b1] => {
            let aa = codon_residue(b1, b'N', b'N');
            if aa != b'X' {
                protein.push(aa);
            }
        }
        [b1, b2] => {
            let aa = codon_residue(b1, b2, b'N');
            if aa != b'X' {
                protein.push(aa);
            }
        }
        _ => {}
    }
    for residue in &mut protein {
        if *residue == b'*' {
            *residue = b'X';
        }
    }
    protein
}

/// `GetComplementSeq`: the reverse complement in upper case; a byte that is not an IUPACna
/// letter (in either case, `U` read as `T`) is left out.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:3061-3121
// ```c++
// string GetComplementSeq(const string &rSeq)
// {
// 	string result(k_strEmptyString);
// 	result.reserve(rSeq.size());
// 	for (string::const_reverse_iterator riter = rSeq.rbegin(), riterEnd = rSeq.rend(); riter != riterEnd; ++riter)
// 	{
// 		char cap = toupper(*riter);
// 		switch (cap)
// 		{
// 		case 'A':
// 			result.push_back('T');
// 			break;
// 		...
// 		case 'T':
// 		case 'U':
// 			result.push_back('A');
// 			break;
// 		case 'R':	//A/G, purine
// 			result.push_back('Y');
// 		...
// 		case 'N':	//anything
// 			result.push_back('N');
// 			break;
// 		}
// 	}
// 	return result;
// }
// ```
// (`toupper` of a byte above 0x7F gives no ASCII letter, so such a byte is left out too.)
pub(crate) fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .filter_map(|&byte| {
            Some(match byte.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' | b'U' => b'A',
                b'R' => b'Y',
                b'Y' => b'R',
                b'M' => b'K',
                b'K' => b'M',
                b'S' => b'W',
                b'W' => b'S',
                b'H' => b'D',
                b'D' => b'H',
                b'B' => b'V',
                b'V' => b'B',
                b'N' => b'N',
                _ => return None,
            })
        })
        .collect()
}

/// `TranslateAll` with genetic code 1: the translations of the three plus-strand frames (the
/// sequence from its first, second and third base) and of the three minus-strand frames (its
/// complement from its first, second and third base). `None` where NCBI erases the first byte
/// of an empty string (a sequence of 1 or 2 bases, or a complement of fewer than 3): libstdc++
/// then moves `SIZE_MAX` bytes (undefined, a crash).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:113-145
// ```c++
// void TranslateAll(const string &NASeq, int gcode, string dimAaData[])
// {
// 	for (int i = 0; i < READINGFRAME::TOTAL_RFS; ++i) dimAaData[i].clear();
//
// 	size_t seqLen = NASeq.size();
// 	if (seqLen > 0)	//at least one
// 	{
// 		string trans(NASeq), cmpNa(GetComplementSeq(NASeq));
// 		CGenetic_code gc;
// 		CRef<CGenetic_code::C_E> ce(new CGenetic_code::C_E);
// 		ce->SetId(gcode);
// 		gc.Set().push_back(ce);
// 		//CmplConvert(NASeq, cmpNa);
// 		// -- positive strand, rfidx 0 -2 (rf 1 - 3)
// 		for (int i = 0; i < READINGFRAME::RF_SIZE; ++i)
// 		{
// 			CSeqTranslator::Translate(trans, dimAaData[i], CSeqTranslator::fIs5PrimePartial/* | CSeqTranslator::fIs3PrimePartial*/, &gc);
// 			for (size_t j = 0, jEnd = dimAaData[i].size(); j < jEnd; ++j)
// 				if ('*' == dimAaData[i][j]) dimAaData[i][j] = 'X';
//
// 			trans.erase(trans.begin());
// 		}
//
// 		//rfidx 3 -5 (rf -1 - -3)
// 		for (int i = READINGFRAME::RF_SIZE; i < READINGFRAME::TOTAL_RFS; ++i)
// 		{
// 			CSeqTranslator::Translate(cmpNa, dimAaData[i], CSeqTranslator::fIs5PrimePartial/* | CSeqTranslator::fIs3PrimePartial*/, &gc);
// 			for (size_t j = 0, jEnd = dimAaData[i].size(); j < jEnd; ++j)
// 				if ('*' == dimAaData[i][j]) dimAaData[i][j] = 'X';
//
// 			cmpNa.erase(cmpNa.begin());
// 		}
// 	}
// }
// ```
// The genetic code is `TDocsum::m_iGenCode`, which the archive never sets (`basedata.cpp:431`,
// `m_iGenCode = 1`; the `gcode` argument of `ProcessBlastResults` is not used), so rpsbproc
// translates with the standard code whatever `-query_gencode` the search used.
pub(crate) fn translate_all(na: &[u8]) -> Option<[Vec<u8>; 6]> {
    let mut frames: [Vec<u8>; 6] = Default::default();
    if na.is_empty() {
        return Some(frames);
    }
    let minus = complement(na);
    if na.len() < 3 || minus.len() < 3 {
        return None;
    }
    for i in 0..3 {
        frames[i] = translate(&na[i..]);
        frames[i + 3] = translate(&minus[i..]);
    }
    Some(frames)
}

#[cfg(test)]
#[path = "frames_tests.rs"]
mod tests;
