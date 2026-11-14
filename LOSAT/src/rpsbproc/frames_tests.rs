use super::*;

// The standard code (gc.inc:43) and `x_InitFsaTransl`'s merges (Genetic_code_table.cpp:139-232,
// quoted at `amino_acid_table`): every unambiguous codon, the Asx/Glx/Leu-Ile merges, a codon
// that spans two amino acids, gaps and bytes that are not bases.
#[test]
fn codons_of_the_standard_code() {
    let bases = b"TCAG";
    let mut index = 0;
    for &b1 in bases {
        for &b2 in bases {
            for &b3 in bases {
                assert_eq!(
                    codon_residue(b1, b2, b3),
                    NCBIEAA_1[index],
                    "{b1} {b2} {b3}"
                );
                index += 1;
            }
        }
    }
    assert_eq!(codon_residue(b'a', b't', b'g'), b'M');
    assert_eq!(codon_residue(b'U', b'G', b'G'), b'W');
    // GCN = Ala, RAY = Asp or Asn, SAR = Gln or Glu, MTH = Ile or Leu.
    assert_eq!(codon_residue(b'G', b'C', b'N'), b'A');
    assert_eq!(codon_residue(b'R', b'A', b'Y'), b'B');
    assert_eq!(codon_residue(b'S', b'A', b'R'), b'Z');
    assert_eq!(codon_residue(b'M', b'T', b'H'), b'J');
    // TAR and TRA: stops only.
    assert_eq!(codon_residue(b'T', b'A', b'R'), b'*');
    assert_eq!(codon_residue(b'T', b'R', b'A'), b'*');
    // TAN: Tyr or a stop; a gap, `X` (as N) and other bytes (as a gap) are X.
    assert_eq!(codon_residue(b'T', b'A', b'N'), b'X');
    assert_eq!(codon_residue(b'-', b'A', b'A'), b'X');
    assert_eq!(codon_residue(b'G', b'G', b'X'), b'G');
    assert_eq!(codon_residue(b'G', b'G', b'Z'), b'X');
    assert_eq!(codon_residue(0xff, b'G', b'G'), b'X');
    // NCBI4na codes as bytes (code 4 = G, 2 = C, 15 = N).
    assert_eq!(codon_residue(4, 2, 15), b'A');
}

// `x_Translate` (sequence.cpp:3636-3725, quoted at `translate`): the last partial codon is kept
// only when it is not X; rpsbproc's `*` to `X` (biodata_core.cpp:130-133).
#[test]
fn translation_of_a_frame() {
    assert_eq!(translate(b"ATGTAA"), b"MX");
    assert_eq!(translate(b"ATGGC"), b"MA");
    assert_eq!(translate(b"ATGTA"), b"M");
    assert_eq!(translate(b"ATGG"), b"M");
    assert_eq!(translate(b"GC"), b"A");
    assert_eq!(translate(b""), b"");
}

// `GetComplementSeq` (basedata.cpp:3061-3121): reversed, upper case, letters that are not
// IUPACna bases left out.
#[test]
fn complement_drops_what_is_not_a_base() {
    assert_eq!(complement(b"acgtuRYMKSWHDBVN"), b"NBVHDSWMKRYAACGT");
    assert_eq!(complement(b"AC-XGT"), b"ACGT");
}

// `TranslateAll` (biodata_core.cpp:113-145): the three frames of the sequence and of its
// complement; `None` where NCBI erases from an empty string.
#[test]
fn six_frames_and_their_limits() {
    let frames = translate_all(b"ATGAAACCCGGGT").expect("long enough");
    assert_eq!(frames[0], b"MKPG");
    assert_eq!(frames[1], b"XNPG");
    assert_eq!(frames[2], b"ETRV");
    // complement: ACCCGGGTTTCAT
    assert_eq!(frames[3], b"TRVS");
    assert_eq!(frames[4], b"PGFH");
    assert_eq!(frames[5], b"PGF");
    assert!(translate_all(b"").expect("empty").iter().all(Vec::is_empty));
    assert!(translate_all(b"AC").is_none());
    assert!(translate_all(b"ACXXX").is_none());
    assert!(translate_all(b"ACG").is_some());
}

// READINGFRAME (basedata.hpp:200-237): the unsigned arithmetic wraps.
#[test]
fn frame_conversions() {
    assert_eq!(plus_rf_idx_to_minus(0, 5000), 5000 % 3 + 3);
    // (1 - 2) wraps to 2^32 - 1, a multiple of 3.
    assert_eq!(plus_rf_idx_to_minus(2, 1), 3);
    assert_eq!(plus_pr_to_plus_na(10, 2), 32);
    assert_eq!(plus_pr_to_plus_na(-1, 0), -3);
    // na_len - (pr - 1) * 3 - rfidx - 1
    assert_eq!(minus_pr_to_plus_na(0, 3, 5000), 4999);
    assert_eq!(minus_pr_to_plus_na(10, 5, 5000), 5000 - 27 - 6);
    // `pr - 1` of `INT_MIN` overflows `int` (undefined): its value, -2^31 - 1, is 2^31 - 1
    // modulo 2^32 in the product (PD-LOSAT-NCBI-DEFECTS exception 11); the residue after it
    // (`INT_MIN + 1`, no overflow) gives the base 3 below.
    let min = i64::from(i32::MIN);
    assert_eq!(minus_pr_to_plus_na(min, 3, 5000), -2_147_478_649);
    assert_eq!(minus_pr_to_plus_na(min + 1, 3, 5000), -2_147_478_649 - 3);
    // A residue past `int`'s range (a mapped position that NCBI's `int` arithmetic took there)
    // converts modulo 2^32.
    assert_eq!(minus_pr_to_plus_na(min - 1, 3, 5000), -2_147_478_649 + 3);
    assert_eq!(
        plus_pr_to_plus_na(i64::from(i32::MAX) + 1, 0),
        plus_pr_to_plus_na(min, 0)
    );
}
