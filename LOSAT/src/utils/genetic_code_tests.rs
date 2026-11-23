use super::*;

// The letters that the bases pair with, the case kept, and every other byte unchanged.
#[test]
fn complement_pairs_the_ncbi4na_letters_and_keeps_every_other_byte() {
    let (letters, complements) = (b"AGCTYRWSKMDVHBN", b"TCGARYWSMKHBDVN");
    for byte in 0..=255u8 {
        let upper = letters
            .iter()
            .position(|&letter| letter == byte)
            .map(|at| complements[at]);
        let lower = letters
            .iter()
            .position(|&letter| letter + 32 == byte)
            .map(|at| complements[at] + 32);
        assert_eq!(complement(byte), upper.or(lower).unwrap_or(byte), "{byte}");
    }
    assert_eq!(complement(b'U'), b'U');
    assert_eq!(complement(b'-'), b'-');
}

#[test]
fn revcomp_reverses_and_complements_keeping_the_case() {
    assert_eq!(revcomp(b"ACGTN"), b"NACGT");
    assert_eq!(revcomp(b"GaTtaCA"), b"TGtaAtC");
    assert_eq!(revcomp(b"AGCTYRWSKMDVHBN"), b"NVDBHKMSWYRAGCT");
    assert_eq!(revcomp(b""), b"");
}

// EXPERIMENT: the table path of `get` against NCBI's nested loops.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:375-379,394-405
// The test compares the table with `base_mask` and with a copy of this loop for all
// codes and for 25 symbols per codon position.
#[test]
fn x_base_index_matches_base_mask_for_every_byte() {
    for byte in 0..=255u8 {
        let expected = match base_mask(byte) {
            8 => 0,
            2 => 1,
            1 => 2,
            4 => 3,
            _ => 4,
        };
        assert_eq!(X_BASE_INDEX[byte as usize], expected, "byte {byte}");
    }
    let slow = |code: &GeneticCode, codon: &[u8; 3]| -> u8 {
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
            for j in 0..4 {
                for k in 0..4 {
                    if masks[0] & bit[i] == 0 || masks[1] & bit[j] == 0 || masks[2] & bit[k] == 0 {
                        continue;
                    }
                    let translated = code.table[i * 16 + j * 4 + k];
                    if aa == 0 {
                        aa = translated;
                    } else if translated != aa {
                        return b'X';
                    }
                }
            }
        }
        aa
    };
    let letters = b"ACGTUacgtuNRYKMSWBDHVn-*X";
    for &(id, _) in TABLES {
        let code = GeneticCode::try_from_id(id).unwrap();
        for &a in letters {
            for &b in letters {
                for &c in letters {
                    let codon = [a, b, c];
                    // `get` is the nested loops unless the switch is set.
                    assert_eq!(code.get(&codon), slow(&code, &codon), "code {id} {codon:?}");
                    let unambiguous = codon.iter().all(|base| b"ACGTUacgtu".contains(base));
                    match code.x_get_unambiguous(&codon) {
                        Some(aa) => {
                            assert!(unambiguous, "code {id} {codon:?}");
                            assert_eq!(aa, slow(&code, &codon), "code {id} {codon:?}");
                        }
                        None => assert!(!unambiguous, "code {id} {codon:?}"),
                    }
                }
            }
        }
    }
}

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
    // Stage A comparison-only C++ API oracle looks up code 32 with `FindGeneticCode`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:588-600:
    let oracle =
        include_str!("../../../docs/evidence/tlosan_stage_a/api_20260923_verified/code32_fmt6.out");
    let standard_oracle = include_str!(
        "../../../docs/evidence/tlosan_stage_a/api_20260923_verified/subject32_g1_fmt6.out"
    );
    assert!(oracle.starts_with("q1\ts_code32\t100.000\t120\t0\t0\t1\t120\t1\t360\t3.39e-93\t251\n"));
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
