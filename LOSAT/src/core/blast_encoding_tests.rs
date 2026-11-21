use super::*;

/// `encode_subject_ncbi2na_packed` as it was written before its one-pass form: every
/// letter to ncbi4na, all of them resolved, then packed.
fn encode_subject_ncbi2na_packed_two_pass(seq: &[u8]) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }
    let ncbi4na: Vec<u8> = seq
        .iter()
        .map(|&letter| iupacna_to_ncbi4na(letter).unwrap_or(15))
        .collect();
    let codes = resolve_ncbi4na_to_ncbi2na(&ncbi4na);
    let mut packed = vec![0u8; codes.len() / COMPRESSION_RATIO + 1];
    for (index, &code) in codes.iter().enumerate() {
        let shift = 6 - 2 * (index % COMPRESSION_RATIO);
        packed[index / COMPRESSION_RATIO] |= code << shift;
    }
    *packed.last_mut().expect("at least one byte") |= (codes.len() % COMPRESSION_RATIO) as u8;
    packed
}

#[test]
fn subject_packing_equals_the_two_pass_form() {
    // Every byte value alone and between unambiguous letters.
    for byte in 0..=255u8 {
        for seq in [vec![byte], vec![b'A', byte, b'c', b'G', b't', byte, byte]] {
            assert_eq!(
                encode_subject_ncbi2na_packed(&seq),
                encode_subject_ncbi2na_packed_two_pass(&seq),
                "{seq:?}"
            );
        }
    }
    // Mixed letters of every length up to 64 and a few long ones (a linear congruential
    // generator picks them; one letter in eight is ambiguous or not a letter).
    let alphabet = b"ACGTacgtACGTacgtACGTacgtACGTacgtMRSVWYHKDBNmrsvwyhkdbn-UuX*";
    let mut state = 0x2545_f491_u32;
    let mut next = || {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        state >> 8
    };
    let lengths = (0..=64).chain([1_000, 4_097, 100_003]);
    for length in lengths {
        for _ in 0..4 {
            let seq: Vec<u8> = (0..length)
                .map(|_| {
                    let pick = next() as usize;
                    if pick.is_multiple_of(8) {
                        alphabet[32 + pick / 8 % (alphabet.len() - 32)]
                    } else {
                        alphabet[pick / 8 % 32]
                    }
                })
                .collect();
            assert_eq!(
                encode_subject_ncbi2na_packed(&seq),
                encode_subject_ncbi2na_packed_two_pass(&seq),
                "length {length}"
            );
        }
    }
}
