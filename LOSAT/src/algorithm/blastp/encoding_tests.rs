use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
#[test]
fn test_encode_protein_sequence_adds_sentinels() {
    let encoded = encode_protein_sequence(b"ACD");
    assert_eq!(encoded.aa_len, 3);
    assert_eq!(encoded.aa_seq.len(), 5);
    assert_eq!(encoded.aa_seq[0], 0);
    assert_eq!(encoded.aa_seq[4], 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-121
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
// `-soft_masking true` keeps the SEG intervals of the hard query and its letters.
#[test]
fn test_soft_masked_query_keeps_letters_and_seg_intervals() {
    let seq = b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKQQQQQQQQQQQQQQQQQQQQQGWENYCWPANMSWKVLMKYEFSSD";
    let seg = SegParams::default();
    let hard = encode_protein_query_frame_with_seg(seq, Some(&seg));
    let soft = encode_protein_query_frame_soft_masked(seq, Some(&seg), false);
    assert!(!hard.seg_masks.is_empty());
    assert_eq!(soft.seg_masks, hard.seg_masks);
    assert_eq!(soft.aa_seq, encode_protein_sequence(seq).aa_seq);
    assert_eq!(hard.aa_seq_nomask.as_deref(), Some(&soft.aa_seq[..]));
    assert!(soft.aa_seq_nomask.is_none());
    assert_ne!(soft.aa_seq, hard.aa_seq);
    let none = encode_protein_query_frame_soft_masked(seq, None, false);
    assert!(none.seg_masks.is_empty() && none.aa_seq == soft.aa_seq);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1253-1254
// A lower-case run joins the SEG intervals; hard writes X and keeps the letters in the copy,
// soft keeps the letters.
#[test]
fn test_lcase_query_masks_join_seg_intervals() {
    let upper = b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKQQQQQQQQQQQQQQQQQQQQQGWENYCWPANMSWKVLMKYEFSSD";
    let mut seq = upper.to_vec();
    seq[2..6].make_ascii_lowercase();
    let seg = SegParams::default();
    let letters = encode_protein_sequence(upper).aa_seq;

    let hard = encode_protein_query_frame_lcase(&seq, Some(&seg), false, false);
    let without = encode_protein_query_frame_with_seg(upper, Some(&seg));
    let mut expected = vec![(2, 6)];
    expected.extend_from_slice(&without.seg_masks);
    assert_eq!(hard.seg_masks, expected);
    assert_eq!(hard.aa_seq_nomask.as_deref(), Some(&letters[..]));
    assert_eq!(hard.aa_seq[3..7], [ncbistdaa::X; 4]);

    let soft = encode_protein_query_frame_lcase(&seq, Some(&seg), true, false);
    assert_eq!(soft.seg_masks, expected);
    assert_eq!(soft.aa_seq, letters);
    assert!(soft.aa_seq_nomask.is_none());

    // No lower case and no SEG: no copy, as without the flag.
    let none = encode_protein_query_frame_lcase(upper, None, false, false);
    assert!(none.seg_masks.is_empty() && none.aa_seq_nomask.is_none());
    assert_eq!(none.aa_seq, letters);
}

#[test]
fn test_ncbistdaa_to_ascii_roundtrip() {
    let encoded = encode_protein_sequence(b"Acd");
    let letters: String = encoded.aa_seq[1..encoded.aa_seq.len() - 1]
        .iter()
        .map(|&b| ncbistdaa_to_ascii(b))
        .collect();
    assert_eq!(letters, "ACD");
}
