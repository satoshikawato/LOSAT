use super::{
    blastna_reverse_complement_table, encode_reverse_complement_blastna_into, reverse_complement,
};
use crate::core::blast_encoding::{encode_iupac_to_blastna, encode_iupac_to_blastna_into};

// The table is the composition of the two functions for every byte.
#[test]
fn blastna_reverse_complement_table_is_the_composition() {
    let table = blastna_reverse_complement_table();
    for b in 0..=255u8 {
        assert_eq!(
            table[b as usize],
            encode_iupac_to_blastna(&reverse_complement(&[b]))[0],
            "byte {b}"
        );
    }
}

// The BLASTNA table gives an uppercase and a lowercase letter the same code, so the
// karlin counts need no uppercase copy of the query.
#[test]
fn blastna_encoding_ignores_letter_case() {
    for b in 0..=255u8 {
        assert_eq!(
            encode_iupac_to_blastna(&[b]),
            encode_iupac_to_blastna(&[b.to_ascii_uppercase()]),
            "byte {b}"
        );
    }
}

#[test]
fn strand_encodings_into_equal_the_copying_forms() {
    let mut seqs: Vec<Vec<u8>> = vec![
        Vec::new(),
        b"ACGT".to_vec(),
        b"acgt".to_vec(),
        b"ACgtRYKMSWBDHVNrykmswbdhvnUu-*0123456789".to_vec(),
        b"ACGTRYKMSWBDHVNUacgtrykmswbdhvnu-*".to_vec(),
        vec![b'A', 128, 200, 255, b'c', 0, 127],
    ];
    // Every byte alone and between letters.
    for b in 0..=255u8 {
        seqs.push(vec![b]);
        seqs.push(vec![b'A', b, b'g']);
    }
    for seq in &seqs {
        let mut plus = vec![0u8; seq.len()];
        encode_iupac_to_blastna_into(seq, &mut plus);
        assert_eq!(plus, encode_iupac_to_blastna(seq), "plus {seq:?}");
        let mut minus = vec![0u8; seq.len()];
        encode_reverse_complement_blastna_into(seq, &mut minus);
        assert_eq!(
            minus,
            encode_iupac_to_blastna(&reverse_complement(seq)),
            "minus {seq:?}"
        );
    }
}

// A lowercase (soft-masked) letter is the same base on the minus strand.
#[test]
fn lowercase_letters_are_complemented_like_uppercase() {
    assert_eq!(reverse_complement(b"ACgtRn-"), b"-NYACGT");
}
