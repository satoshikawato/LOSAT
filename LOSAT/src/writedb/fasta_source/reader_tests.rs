use super::*;

const ALPHABETS: [Alphabet; 3] = [
    Alphabet::Protein,
    Alphabet::NucleotideGaps,
    Alphabet::Nucleotide,
];

fn xorshift64() -> impl FnMut() -> u64 {
    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    }
}

/// The residue test of each alphabet is its switch's (`actions()[byte] >= RESIDUE_MIN`) for every
/// byte, and the switch stores a residue as the byte in upper case.
#[test]
fn residue_tests_match_the_tables() {
    for alphabet in ALPHABETS {
        let actions = alphabet.actions();
        for byte in 0..=u8::MAX {
            let residue = actions[byte as usize] >= RESIDUE_MIN;
            assert_eq!(
                alphabet.is_residue(byte),
                residue,
                "{alphabet:?} {byte:#04x}"
            );
            if residue {
                assert_eq!(
                    actions[byte as usize],
                    byte.to_ascii_uppercase(),
                    "{byte:#04x}"
                );
            }
        }
    }
}

/// The run that `residue_run` finds is the one that a scan of the switch one byte at a time finds,
/// for every byte at every place of lines across the blocks of the search, and on random lines.
#[test]
fn residue_runs_match_the_byte_by_byte_scan() {
    let scan = |alphabet: Alphabet, s: &[u8]| {
        let actions = alphabet.actions();
        s.iter()
            .position(|&byte| actions[byte as usize] < RESIDUE_MIN)
            .unwrap_or(s.len())
    };
    for alphabet in ALPHABETS {
        for byte in 0..=u8::MAX {
            for len in [1, 2, 31, 32, 33, 63, 64, 65, 80, 127, 128, 129, 200] {
                for at in [0, 1, len / 2, len - 1] {
                    if at >= len {
                        continue;
                    }
                    let mut line = b"acgtACGTnN".repeat(len / 10 + 1);
                    line.truncate(len);
                    line[at] = byte;
                    assert_eq!(
                        alphabet.residue_run(&line),
                        scan(alphabet, &line),
                        "{alphabet:?} {byte:#04x} {len} {at}"
                    );
                }
            }
        }
        let mut next = xorshift64();
        const RESIDUES: &[u8] = b"ACGTNacgtnRYKMBDHVWSUEFILPQXZ*-";
        for _ in 0..3000 {
            let len = (next() % 300) as usize;
            let line: Vec<u8> = (0..len)
                .map(|_| match next() % 64 {
                    0 => next() as u8,
                    _ => RESIDUES[(next() % RESIDUES.len() as u64) as usize],
                })
                .collect();
            assert_eq!(
                alphabet.residue_run(&line),
                scan(alphabet, &line),
                "{alphabet:?} {line:?}"
            );
        }
    }
}

/// `extend_upper` appends what the byte-by-byte map appends, for every byte at every place of runs
/// of every length up to 100 and across blocks, after bytes already stored.
#[test]
fn extend_upper_is_the_byte_by_byte_map() {
    let stored = b"MKV".to_vec();
    let mut lengths: Vec<usize> = (0..=100).collect();
    lengths.extend([127, 128, 129, 300]);
    for len in lengths {
        for byte in 0..=u8::MAX {
            for at in [0, 1, 31, 32, 33, len / 2, len.saturating_sub(1)] {
                let mut run = b"acgtACGTnN*-".repeat(len / 12 + 1);
                run.truncate(len);
                if at < len {
                    run[at] = byte;
                }
                let mut expected = stored.clone();
                expected.extend(run.iter().map(u8::to_ascii_uppercase));
                let mut out = stored.clone();
                extend_upper(&mut out, &run);
                assert_eq!(out, expected, "{byte:#04x} {len} {at}");
            }
        }
    }
}

/// `data_line_counts` gives the counts of `CheckDataLine`'s loop, one byte at a time with its
/// stop at `;`, for every byte at every place and on random lines of up to 70 bytes.
#[test]
fn data_line_counts_match_the_loop() {
    let reference = |checked: &[u8], is_nuc: bool| {
        let (mut good, mut bad, mut ambig_nuc) = (0usize, 0usize, 0usize);
        for &c in checked {
            if c.is_ascii_alphabetic() || c == b'*' {
                good += 1;
                if is_nuc && is_ambig_nuc(c) {
                    ambig_nuc += 1;
                }
            } else if c == b'-' {
                good += 1;
            } else if is_space(c) || c.is_ascii_digit() {
                // White space and digits are neutral.
            } else if c == b';' {
                break;
            } else {
                bad += 1;
            }
        }
        (good, bad, ambig_nuc)
    };
    for is_nuc in [false, true] {
        for byte in 0..=u8::MAX {
            for len in [1, 2, 31, 32, 33, 60, 63, 64, 70] {
                for at in [0, len / 2, len - 1] {
                    let mut line = b"ACGTNRYacgt".repeat(7);
                    line.truncate(len);
                    line[at] = byte;
                    assert_eq!(
                        data_line_counts(&line, is_nuc),
                        reference(&line, is_nuc),
                        "{is_nuc} {byte:#04x} {len} {at}"
                    );
                }
            }
        }
        let mut next = xorshift64();
        const BYTES: &[u8] = b"ACGTNacgtnRYKMBDHVWSUEFILPQXZ*- ;0912\t";
        for _ in 0..5000 {
            let len = (next() % 71) as usize;
            let line: Vec<u8> = (0..len)
                .map(|_| match next() % 16 {
                    0 => next() as u8,
                    _ => BYTES[(next() % BYTES.len() as u64) as usize],
                })
                .collect();
            assert_eq!(
                data_line_counts(&line, is_nuc),
                reference(&line, is_nuc),
                "{is_nuc} {line:?}"
            );
        }
    }
}
