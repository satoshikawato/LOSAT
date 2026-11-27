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

/// The run that `residue_run` finds is the one that a scan of the switch one byte at a time finds:
/// for every byte value at boundary places (first, second, prefix edge, middle and last) across
/// search block sizes, for representative boundary bytes at every place of lines up to 64 bytes,
/// and on random lines (every place of `find_first` is tested in `byte_search_tests.rs`).
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
            for len in [
                1, 2, 10, 11, 12, 16, 27, 28, 31, 32, 33, 47, 48, 49, 63, 64, 65, 80, 127, 128,
                129, 200,
            ] {
                for at in [0, 1, 10, 11, 12, len / 2, len - 1] {
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
        // Every place of lines up to 64 bytes for representative residue and non-residue bytes.
        for &byte in b" ;-1*E\x00\xff" {
            for len in 1..=64usize {
                for at in 0..len {
                    let mut line = b"acgtACGT".repeat(8);
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

/// `extend_upper` appends what the byte-by-byte map appends, for every byte at block-boundary
/// places (`0`, `1`, `31`, `32`, `33`, `len / 2`, `len - 1`) of runs of every length up to 100
/// and across blocks, after bytes already stored.
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
/// stop at `;`, for every byte at the first, middle and last place of lines across block edges,
/// for representative bytes at every place up to 70 bytes, and on random lines of up to 70 bytes.
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
            for len in [1, 2, 7, 8, 9, 15, 16, 31, 32, 33, 48, 60, 63, 64, 70] {
                for at in [0, len / 2, len - 1] {
                    let mut line = b"ACGTacgtAC".repeat(7);
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
        for &byte in b"U; -1*N\x00\xff" {
            for len in [1, 8, 10, 16, 20, 32, 48, 60, 70] {
                for at in 0..len {
                    let mut line = b"ACGTacgtAC".repeat(7);
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

/// A record of `bases` in lines of 60 after the definition line `>name`.
fn record_of(name: &str, bases: &[u8]) -> Vec<u8> {
    let mut out = format!(">{name}\n").into_bytes();
    for line in bases.chunks(60) {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
    out
}

/// `count` bases from a xorshift sequence, in lower case when `lower` (the reader stores them in
/// upper case).
fn bases(count: usize, lower: bool) -> Vec<u8> {
    let mut next = xorshift64();
    (0..count)
        .map(|_| {
            let base = b"ACGT"[(next() % 4) as usize];
            if lower {
                base.to_ascii_lowercase()
            } else {
                base
            }
        })
        .collect()
}

/// makeblastdb's reader keeps its residue buffer from record to record, after a delta as after a
/// raw sequence: `residues` holds the record's residues (upper case, without the gaps), which the
/// delta's literals index; a record that fits is read into the same buffer, and a longer one
/// grows it from the capacity that the records before it left. Records without residues and the
/// last record are read the same way.
#[test]
fn the_residue_buffer_is_kept_after_a_delta() {
    let runs = |parts: &[&[u8]]| parts.concat();
    let long = bases(3000, false);
    let mut input = Vec::new();
    // (record, its residues, a delta).
    let mut expected: Vec<(Vec<u8>, bool)> = Vec::new();
    // A delta of NCBI2na and NCBI4na literals (a run of N between unambiguous runs).
    let d1 = runs(&[&long[..], &[b'N'; 20][..], &long[..]]);
    input.extend(record_of("d1", &d1));
    expected.push((d1, true));
    // A raw sequence, shorter.
    let r1 = bases(100, true);
    input.extend(record_of("r1", &r1));
    expected.push((r1.to_ascii_uppercase(), false));
    // A delta with a gap line and hyphens (gaps of `fParseGaps`), shorter.
    input.extend(b">d2\nACGTACGTAC\n>?10\nACGT--ACGTA-\nGGCC\n");
    expected.push((b"ACGTACGTACACGTACGTAGGCC".to_vec(), true));
    // A longer delta, which grows the buffer.
    let d3 = runs(&[
        &long[..],
        &[b'N'; 5][..],
        &long[..],
        &long[..],
        &[b'R'; 3][..],
        &long[..],
    ]);
    input.extend(record_of("d3", &d3));
    expected.push((d3, true));
    // A record without residues (virtual), then a shorter delta as the last record.
    input.extend(b">e\n");
    expected.push((Vec::new(), false));
    let d4 = runs(&[&long[..500], &[b'N'; 64][..], &long[..700]]);
    input.extend(record_of("d4", &d4));
    expected.push((d4, true));

    let mut reader = FastaReader::new(&input[..], ReaderFlags::fasta_bioseq_source(false, false));
    let mut capacity = 0;
    let mut buffer = std::ptr::null();
    for (index, (residues, delta)) in expected.iter().enumerate() {
        let record = reader.read_one_seq(&mut |_| ()).expect("a record");
        assert_eq!(
            matches!(record.inst, SeqInst::Delta { .. }),
            *delta,
            "record {index}"
        );
        assert_eq!(reader.residues(), &residues[..], "record {index}");
        let seq_data = &reader.state.seq_data;
        if seq_data.len() <= capacity {
            assert_eq!(
                (seq_data.capacity(), seq_data.as_ptr()),
                (capacity, buffer),
                "record {index} is read into the same buffer"
            );
        } else {
            assert!(seq_data.capacity() > capacity, "record {index}");
        }
        capacity = seq_data.capacity();
        buffer = seq_data.as_ptr();
    }
    assert!(reader.at_eof());
    // The records that grew the buffer were the first and the longer delta.
    assert!(capacity >= expected[3].0.len());
}

/// A masking program takes a delta's residues with their buffer, as before (`iupac_letters`), and
/// the next record grows its own; a raw sequence's residues leave a buffer of their capacity.
#[test]
fn a_masking_program_takes_the_residues_with_their_buffer() {
    let long = bases(3000, false);
    let delta = [&long[..], &[b'N'; 20][..], &long[..]].concat();
    let mut input = record_of("d", &delta);
    input.extend(record_of("r", &long));
    let mut reader = FastaReader::new(&input[..], ReaderFlags::mask_fasta_reader(false, false));
    let record = reader.read_one_seq(&mut |_| ()).expect("a record");
    assert!(matches!(record.inst, SeqInst::Delta { .. }));
    assert_eq!(reader.take_delta_residues(), delta);
    assert_eq!(reader.state.seq_data.capacity(), 0);
    let record = reader.read_one_seq(&mut |_| ()).expect("a record");
    assert!(matches!(record.inst, SeqInst::Raw { .. }));
    assert_eq!(reader.take_raw_residues(), long);
    assert_eq!(reader.state.seq_data.capacity(), long.len());
}
