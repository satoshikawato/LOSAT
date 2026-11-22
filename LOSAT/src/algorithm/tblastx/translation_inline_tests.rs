use super::*;

/// `translate_frame_into` (and `frame_length`) against `generate_frames`: each frame's
/// residues, appended after what `out` holds, and length, and no frame (`out` unchanged) where
/// `generate_frames` makes none (also frames 0, 4 and -4), for every length up to 40 and a long
/// sequence, with every byte the readers give (IUPAC letters, `-`, both cases), under two codes.
#[test]
fn translate_frame_into_is_the_frame_of_generate_frames() {
    let letters = b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn";
    let mut state = 7u64;
    let mut sequence = |length: usize| -> Vec<u8> {
        (0..length)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let draw = (state >> 33) as usize;
                // mostly unambiguous bases, so stop codons and ordinary residues occur
                if draw.is_multiple_of(4) {
                    letters[draw / 4 % letters.len()]
                } else {
                    b"ACGTacgt"[draw / 4 % 8]
                }
            })
            .collect()
    };
    let mut compared = 0;
    for id in [1, 4] {
        let code = GeneticCode::try_from_id(id).unwrap();
        let residues = x_codon_residues(&code);
        for length in (0..=40).chain([3_001]) {
            let seq = sequence(length);
            let frames = generate_frames(&seq, &code);
            for frame in [1, 2, 3, -1, -2, -3, 0, 4, -4] {
                let expected = frames
                    .iter()
                    .find(|candidate| candidate.frame == frame)
                    .map(|candidate| candidate.aa_seq[1..1 + candidate.aa_len].to_vec());
                let mut out = b"kept".to_vec();
                let got = translate_frame_into(&seq, frame, &residues, &mut out);
                assert_eq!(got, frame_length(length, frame), "{length} {frame}");
                assert_eq!(got, expected.as_ref().map(Vec::len), "{length} {frame}");
                assert_eq!(&out[..4], b"kept", "{length} {frame}");
                assert_eq!(
                    out[4..],
                    expected.clone().unwrap_or_default()[..],
                    "{length} {frame}"
                );
                compared += usize::from(expected.is_some());
            }
        }
    }
    // lengths 3 and 4 have two and four frames, 5 to 40 and 3,001 all six
    assert_eq!(compared, 2 * (2 + 4 + 6 * 36 + 6));
}

/// `x_frames_in_windows` against `generate_frames`: every field of every frame.
fn assert_same_frames(seq: &[u8], window: usize, code: &GeneticCode) {
    let expected = generate_frames(seq, code);
    let mut fills: Vec<(usize, usize)> = Vec::new();
    let got = x_frames_in_windows(
        seq.len(),
        window,
        &mut |from, out| {
            fills.push((from, out.len()));
            out.copy_from_slice(&seq[from..from + out.len()]);
        },
        code,
    );
    assert_eq!(
        got.len(),
        expected.len(),
        "len {} window {window}",
        seq.len()
    );
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g.frame, e.frame);
        assert_eq!(
            g.aa_seq,
            e.aa_seq,
            "frame {} len {} window {window}",
            e.frame,
            seq.len()
        );
        assert!(g.aa_seq_nomask.is_none() && e.aa_seq_nomask.is_none());
        assert_eq!(g.aa_len, e.aa_len);
        assert_eq!(g.orig_len, e.orig_len);
        assert!(g.seg_masks.is_empty() && e.seg_masks.is_empty());
    }
    // The windows cover the letters once, in order.
    let mut next = 0;
    for (from, n) in fills {
        assert_eq!(from, next);
        assert!(n > 0);
        next += n;
    }
    assert_eq!(next, seq.len());
}

#[test]
fn the_mask_letters_have_their_masks() {
    for (mask, &letter) in X_MASK_LETTERS.iter().enumerate() {
        assert_eq!(crate::utils::genetic_code::base_mask(letter) as usize, mask);
    }
}

#[test]
fn frames_in_windows_are_the_frames_of_the_letters() {
    // Letters of both cases, every IUPAC code, U, a gap and bytes that are not letters.
    let alphabet = b"ACGTACGTACGTacgtNnRYKMSWBDHVrykmswbdhvUu-*X.";
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut letter = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        alphabet[(state % alphabet.len() as u64) as usize]
    };
    let codes = [1u8, 2, 4, 11];
    for len in (0..40).chain([95, 96, 97, 98, 99, 100, 101, 1000, 3001]) {
        let seq: Vec<u8> = (0..len).map(|_| letter()).collect();
        for window in [1, 2, 3, 4, 5, 7, 32, 33, 34, 100, X_FRAME_WINDOW] {
            for &id in &codes {
                assert_same_frames(&seq, window, &GeneticCode::from_id(id));
            }
        }
    }
    let long: Vec<u8> = (0..2 * X_FRAME_WINDOW + 5).map(|_| letter()).collect();
    assert_same_frames(&long, X_FRAME_WINDOW, &GeneticCode::from_id(1));
}

#[test]
fn test_generate_frames() {
    let code = GeneticCode::from_id(1);
    let seq = b"ATGCGATCGATCGATCG";
    let frames = generate_frames(seq, &code);

    // Should generate up to 6 frames
    assert!(frames.len() <= 6);

    // Check that frames have valid frame numbers and sentinel structure
    for frame in &frames {
        assert!(frame.frame >= -3 && frame.frame <= 3 && frame.frame != 0);
        assert_eq!(frame.orig_len, seq.len());
        // Check sentinels
        assert_eq!(frame.aa_seq[0], SENTINEL_BYTE);
        assert_eq!(frame.aa_seq[frame.aa_seq.len() - 1], SENTINEL_BYTE);
        // Check aa_len matches actual amino acids
        assert_eq!(frame.aa_len, frame.aa_seq.len() - 2);
    }
}

#[test]
fn test_sentinel_positions() {
    let code = GeneticCode::from_id(1);
    let seq = b"ATGATGATGATG"; // 4 codons = 4 amino acids
    let frames = generate_frames(seq, &code);

    let frame = &frames[0]; // Frame 1
                            // aa_seq should be: [SENTINEL(0), aa0, aa1, aa2, aa3, SENTINEL(0)]
    assert_eq!(frame.aa_seq.len(), 6); // 4 AAs + 2 sentinels
    assert_eq!(frame.aa_len, 4);
    assert_eq!(frame.aa_seq[0], SENTINEL_BYTE);
    assert_eq!(frame.aa_seq[5], SENTINEL_BYTE);
    // Middle positions should be valid amino acids (not sentinels)
    // Note: NCBISTDAA encodes actual AAs as 1-27, so they are never 0
    for i in 1..5 {
        assert_ne!(frame.aa_seq[i], SENTINEL_BYTE);
    }
}

#[test]
fn test_ncbistdaa_encoding() {
    use crate::utils::matrix::ncbistdaa;

    let code = GeneticCode::from_id(1);
    // ATG = Met (M), NCBISTDAA M = 12
    let seq = b"ATG";
    let frames = generate_frames(seq, &code);
    let frame = &frames[0];
    assert_eq!(frame.aa_seq[1], ncbistdaa::M);

    // TAA = Stop (*), NCBISTDAA * = 25
    let seq_stop = b"TAA";
    let frames_stop = generate_frames(seq_stop, &code);
    let frame_stop = &frames_stop[0];
    assert_eq!(frame_stop.aa_seq[1], ncbistdaa::STOP);
}
