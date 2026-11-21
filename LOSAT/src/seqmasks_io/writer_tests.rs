use super::*;

#[test]
fn fasta_letters_follow_ncbis_one_step_mask_iterator() {
    // A mask list where one position passes two intervals: NCBI advances the iterator
    // once per position, so the second interval [3, 3] after [0, 1] is reached at 2 and
    // its letter 3 is masked; a third interval right after is reached only at its own
    // positions.
    let mut out = Vec::new();
    fasta_letters(
        b"ACGTACGT",
        &[(0, 1), (3, 3), (4, 5)],
        false,
        Vec::new(),
        &mut out,
    )
    .unwrap();
    assert_eq!(out, b"acGtacGT\n");
    let mut out = Vec::new();
    fasta_letters(&[b'A'; 61], &[(59, 60)], true, Vec::new(), &mut out).unwrap();
    assert_eq!(&out[..], [&[b'A'; 59][..], b"N\nN\n"].concat());
}

/// NCBI's per-letter loop as `fasta_letters` had it before lines were masked in place: the
/// letters and the mask walk of every position, then the lines.
fn fasta_letters_per_letter(iupac: &[u8], mask: &[(u32, u32)], hard_masking: bool) -> Vec<u8> {
    let masked = |letter: u8| {
        if hard_masking {
            b'N'
        } else {
            letter.to_ascii_lowercase()
        }
    };
    let mut out = Vec::new();
    let mut imask = 0;
    let mut line = 0;
    for (i, &letter) in iupac.iter().enumerate() {
        let i_pos = i as u64;
        let mut letter = letter;
        if let Some(&(first, second)) = mask.get(imask) {
            if i_pos >= u64::from(first) {
                if i_pos <= u64::from(second) {
                    letter = masked(letter);
                } else {
                    imask += 1;
                    if let Some(&(first, second)) = mask.get(imask) {
                        if i_pos >= u64::from(first) && i_pos <= u64::from(second) {
                            letter = masked(letter);
                        }
                    }
                }
            }
        }
        out.push(letter);
        line += 1;
        if line == 60 {
            out.push(b'\n');
            line = 0;
        }
    }
    if line > 0 {
        out.push(b'\n');
    }
    out
}

#[test]
fn fasta_letters_masked_by_line_equal_the_per_letter_walk() {
    // Lengths around the line and the 64 KiB block, and mask lists sorted, overlapping,
    // unsorted, empty, past the end, inverted and at u32::MAX: the lines masked in place give
    // the bytes of the per-letter walk.
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    for length in [0usize, 1, 59, 60, 61, 119, 120, 121, 1093, 66_000, 70_001] {
        let iupac: Vec<u8> = (0..length).map(|_| b"ACGTN"[next(5) as usize]).collect();
        for case in 0..12 {
            let mut mask: Vec<(u32, u32)> = (0..next(9))
                .map(|_| {
                    let a = next(length as u64 + 80) as u32;
                    (a, a + next(150) as u32)
                })
                .collect();
            if case % 3 != 0 {
                mask.sort_unstable();
            }
            if case == 5 {
                mask.push((u32::MAX - 1, u32::MAX));
            }
            if case == 6 {
                mask.insert(0, (u32::MAX - 1, u32::MAX));
            }
            if case == 7 && !mask.is_empty() {
                // An interval whose end is before its start.
                let (a, b) = mask[0];
                mask[0] = (b + 1, a);
            }
            for hard_masking in [false, true] {
                let mut out = b">x\n".to_vec();
                fasta_letters(&iupac, &mask, hard_masking, b">x\n".to_vec(), &mut out).unwrap();
                let want = [
                    &b">x\n"[..],
                    &b">x\n"[..],
                    &fasta_letters_per_letter(&iupac, &mask, hard_masking),
                ]
                .concat();
                assert_eq!(
                    out, want,
                    "length {length}, mask {mask:?}, hard {hard_masking}"
                );
            }
        }
    }
}
