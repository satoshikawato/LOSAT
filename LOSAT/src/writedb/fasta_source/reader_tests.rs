use super::*;

const ALPHABETS: [Alphabet; 3] = [
    Alphabet::Protein,
    Alphabet::NucleotideGaps,
    Alphabet::Nucleotide,
];

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
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
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
