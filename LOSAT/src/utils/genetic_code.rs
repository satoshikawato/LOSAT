//! NCBI genetic-code tables used by translated searches.

// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:588-603
// ```c++
// const string kGenCode = CGen_code_table::GetNcbieaa(genetic_code);
// if (kGenCode == kEmptyStr) { return retval; }
// ```
// NCBI reference: c++/src/objects/seqfeat/gc.prt:105-357
// Each row below copies its pinned ncbieaa entry.
const TABLES: &[(u8, &[u8; 64])] = &[
    // NCBI gc.prt:109: id 1, ncbieaa "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        1,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:119: id 2, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
    (
        2,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:129: id 3, ncbieaa "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        3,
        b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:140: id 4, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        4,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:150: id 5, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"
    (
        5,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:160: id 6, ncbieaa "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        6,
        b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:170: id 9, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    (
        9,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:180: id 10, ncbieaa "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        10,
        b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:189: id 11, ncbieaa "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        11,
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:198: id 12, ncbieaa "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        12,
        b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:207: id 13, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"
    (
        13,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:216: id 14, ncbieaa "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    (
        14,
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:225: id 15, ncbieaa "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        15,
        b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:234: id 16, ncbieaa "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        16,
        b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:243: id 21, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"
    (
        21,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:252: id 22, ncbieaa "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        22,
        b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:261: id 23, ncbieaa "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        23,
        b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:270: id 24, ncbieaa "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
    (
        24,
        b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:279: id 25, ncbieaa "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        25,
        b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:288: id 26, ncbieaa "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        26,
        b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:297: id 27, ncbieaa "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        27,
        b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:306: id 28, ncbieaa "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        28,
        b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:315: id 29, ncbieaa "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        29,
        b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:324: id 30, ncbieaa "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        30,
        b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:333: id 31, ncbieaa "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        31,
        b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:342: id 32, ncbieaa "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    (
        32,
        b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    // NCBI gc.prt:351: id 33, ncbieaa "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"
    (
        33,
        b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
    ),
];

pub struct GeneticCode {
    pub table: [u8; 64],
}

impl GeneticCode {
    // NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:588-600
    // ```c++
    // const string kGenCode = CGen_code_table::GetNcbieaa(genetic_code);
    // if (kGenCode == kEmptyStr) { return retval; }
    // ```
    // A missing ID is an error at this boundary, never an implicit code-1 table.
    pub fn try_from_id(id: u8) -> Result<Self, String> {
        TABLES
            .iter()
            .find(|(code, _)| *code == id)
            .map(|(_, table)| Self { table: **table })
            .ok_or_else(|| format!("unsupported genetic code {id}"))
    }

    // NCBI reference: c++/src/algo/blast/core/gencode_singleton.c:65-69
    // ```c
    // Uint1* GenCodeSingletonFind(Uint4 gen_code_id) {
    //     return DynamicSGenCodeNodeArray_Find(g_theInstance, gen_code_id);
    // }
    // ```
    // Existing internal TBLASTX callers pass CLI-validated IDs; invalid direct calls fail.
    pub fn from_id(id: u8) -> Self {
        Self::try_from_id(id).unwrap_or_else(|error| panic!("{error}"))
    }

    // NCBI reference: c++/src/algo/blast/core/blast_util.c:369-424
    // ```c
    // static Uint1 mapping[4] = { 8, 2, 1, 4 };
    // for (i = 0; i < 4; i++) if (codon[0] & mapping[i])
    //   for (j = 0; j < 4; j++) if (codon[1] & mapping[j])
    //     for (k = 0; k < 4; k++) if (codon[2] & mapping[k]) {
    //       taa = codes[i*16+j*4+k];
    //       if (!aa) aa = taa; else if (taa != aa) aa = kXResidue;
    //     }
    // ```
    pub fn get(&self, codon: &[u8]) -> u8 {
        if codon.len() != 3 {
            return b'X';
        }
        let masks = [
            base_mask(codon[0]),
            base_mask(codon[1]),
            base_mask(codon[2]),
        ];
        if masks.contains(&0) {
            return b'X';
        }
        let bit = [8, 2, 1, 4];
        let mut aa = 0;
        for i in 0..4 {
            if masks[0] & bit[i] == 0 {
                continue;
            }
            for j in 0..4 {
                if masks[1] & bit[j] == 0 {
                    continue;
                }
                for k in 0..4 {
                    if masks[2] & bit[k] == 0 {
                        continue;
                    }
                    let translated = self.table[i * 16 + j * 4 + k];
                    if aa == 0 {
                        aa = translated;
                    } else if translated != aa {
                        return b'X';
                    }
                }
            }
        }
        aa
    }
}

// NCBI reference: c++/src/algo/blast/core/blast_encoding.c:94-103
// ```c
// const Uint1 IUPACNA_TO_NCBI4NA[128]={
//  0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
//  0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0,
// ```
fn base_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'M' => 3,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_gc_prt_all_27_codes_and_64_codons() {
        let fixture = include_str!("../../../docs/evidence/tlosan_stage_b/gc_prt_27.tsv");
        let rows: Vec<_> = fixture.lines().skip(1).collect();
        assert_eq!(rows.len(), 27);
        for row in rows {
            let fields: Vec<_> = row.split('\t').collect();
            let id: u8 = fields[0].parse().unwrap();
            let expected = fields[1].as_bytes();
            let code = GeneticCode::try_from_id(id).unwrap();
            assert_eq!(expected.len(), 64);
            for (idx, &aa) in expected.iter().enumerate() {
                let alphabet = b"TCAG";
                let codon = [alphabet[idx / 16], alphabet[idx / 4 % 4], alphabet[idx % 4]];
                assert_eq!(code.get(&codon), aa, "code {id}, codon {:?}", codon);
            }
            // NCBI gc.prt:105-357: each ncbieaa row contains exactly 64 residues.
            // ncbieaa "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
            println!("id={id} codons=64 rust_vs_pinned_ncbi=PASS");
        }
    }

    #[test]
    fn ambiguous_and_stop_codons_follow_ncbi4na_expansion() {
        for (id, _, _) in include_str!("../../../docs/evidence/tlosan_stage_b/gc_prt_27.tsv")
            .lines()
            .skip(1)
            .map(|row| {
                let mut f = row.split('\t');
                (
                    f.next().unwrap().parse::<u8>().unwrap(),
                    f.next().unwrap(),
                    f.next().unwrap(),
                )
            })
        {
            let code = GeneticCode::try_from_id(id).unwrap();
            assert_eq!(code.get(b"TTY"), b'F', "code {id}");
            assert_eq!(code.get(b"NNN"), b'X', "code {id}");
            assert_eq!(code.get(b"AT"), b'X', "code {id}");
            assert_eq!(code.get(b"A-G"), b'X', "code {id}");
            assert_eq!(code.get(b"TAA"), code.table[10], "code {id}");
            assert_eq!(code.get(b"TAG"), code.table[11], "code {id}");
            assert_eq!(code.get(b"TGA"), code.table[14], "code {id}");
        }
        let standard = GeneticCode::try_from_id(1).unwrap();
        assert_eq!(standard.get(b"TAR"), b'*');
        let code32 = GeneticCode::try_from_id(32).unwrap();
        assert_eq!(code32.get(b"TAG"), b'W');
        assert_eq!(code32.get(b"TAR"), b'X');
        // Stage A comparison-only C++ API oracle uses FindGeneticCode(32).
        // NCBI c++/src/algo/blast/api/blast_aux.cpp:588-600:
        // const string kGenCode = CGen_code_table::GetNcbieaa(genetic_code);
        let oracle = include_str!(
            "../../../docs/evidence/tlosan_stage_a/api_20260923_verified/code32_fmt6.out"
        );
        let standard_oracle = include_str!(
            "../../../docs/evidence/tlosan_stage_a/api_20260923_verified/subject32_g1_fmt6.out"
        );
        assert!(
            oracle.starts_with("q1\ts_code32\t100.000\t120\t0\t0\t1\t120\t1\t360\t3.39e-93\t251\n")
        );
        assert!(standard_oracle
            .starts_with("q1\ts_code32\t95.833\t120\t5\t0\t1\t120\t1\t360\t1.75e-85\t231\n"));
    }

    #[test]
    fn invalid_ids_never_resolve_to_standard() {
        for id in 0..=u8::MAX {
            if TABLES.iter().all(|(valid, _)| *valid != id) {
                assert_eq!(
                    GeneticCode::try_from_id(id).err(),
                    Some(format!("unsupported genetic code {id}"))
                );
            }
        }
    }
}
