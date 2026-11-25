use super::*;

fn segments(src: &[u8]) -> Vec<(Coding, u32, u32)> {
    pack_iupacna(src)
        .into_iter()
        .map(|segment| (segment.coding, segment.start, segment.length))
        .collect()
}

/// One run of unambiguous bases: NCBI4na for one or two bases (equal costs, the wide
/// arrangement wins), NCBI2na from three (makeblastdb: `G` and `GA` get an ambiguity
/// block, `GAT` gets none).
#[test]
fn one_unambiguous_run() {
    assert_eq!(segments(b""), []);
    assert_eq!(segments(b"G"), [(Coding::Ncbi4na, 0, 1)]);
    assert_eq!(segments(b"GA"), [(Coding::Ncbi4na, 0, 2)]);
    assert_eq!(segments(b"UU"), [(Coding::Ncbi4na, 0, 2)]);
    for length in 3..200 {
        let src = b"GATC".repeat(50);
        assert_eq!(
            segments(&src[..length]),
            [(Coding::Ncbi2na, 0, length as u32)]
        );
    }
}

/// A segment costs 128 bytes, so a run of unambiguous bases is its own NCBI2na
/// segment only where that saves more than the segments it adds.
#[test]
fn mixed_runs() {
    assert_eq!(segments(b"N"), [(Coding::Ncbi4na, 0, 1)]);
    assert_eq!(segments(b"ACGTN"), [(Coding::Ncbi4na, 0, 5)]);
    assert_eq!(segments(b"NACGT"), [(Coding::Ncbi4na, 0, 5)]);
    // 600 bases: 150 bytes of NCBI2na and a segment (128) against 300 of NCBI4na.
    let mut src = b"ACGT".repeat(150);
    src.push(b'N');
    assert_eq!(
        segments(&src),
        [(Coding::Ncbi2na, 0, 600), (Coding::Ncbi4na, 600, 1)]
    );
    // 400 bases: 100 + 128 against 200.
    let mut src = b"ACGT".repeat(100);
    src.push(b'N');
    assert_eq!(segments(&src), [(Coding::Ncbi4na, 0, 401)]);
    // The break-even: 128 + n/4 against n/2 for the first run.
    let mut src = b"ACGT".repeat(128);
    src.push(b'N');
    assert_eq!(segments(&src), [(Coding::Ncbi4na, 0, 513)]);
    let mut src = b"ACGT".repeat(129);
    src.push(b'N');
    assert_eq!(
        segments(&src),
        [(Coding::Ncbi2na, 0, 516), (Coding::Ncbi4na, 516, 1)]
    );
    // A run between two ambiguous ones must pay for two segments.
    let mut src = b"N".to_vec();
    src.extend(b"ACGT".repeat(200));
    src.push(b'N');
    assert_eq!(segments(&src), [(Coding::Ncbi4na, 0, 802)]);
    let mut src = b"N".to_vec();
    src.extend(b"ACGT".repeat(300));
    src.push(b'N');
    assert_eq!(
        segments(&src),
        [
            (Coding::Ncbi4na, 0, 1),
            (Coding::Ncbi2na, 1, 1200),
            (Coding::Ncbi4na, 1201, 1)
        ]
    );
}

/// `CPacker` against the least cost over every assignment of a coding to each run
/// (an NCBI4na-only run must be NCBI4na; a segment costs 128 plus its bytes).
#[test]
fn least_cost_of_all_arrangements() {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    for _ in 0..300 {
        let runs = 1 + next(9) as usize;
        let mut src = Vec::new();
        let mut lengths = Vec::new();
        let mut wide = next(2) == 1;
        for _ in 0..runs {
            let length = 1 + next(if wide { 40 } else { 1400 }) as usize;
            let letters: &[u8] = if wide { b"NRYKM" } else { b"ACGTU" };
            for _ in 0..length {
                src.push(letters[next(5) as usize]);
            }
            lengths.push((wide, length as u32));
            wide = !wide;
        }
        // The cost that `x_AddBoundary` sums: the bytes of each run in its coding, and
        // one overhead for each segment (adjacent runs of one coding are one segment).
        let cost = |codings: &[Coding]| {
            let mut cost = 0;
            for (run, &coding) in codings.iter().enumerate() {
                if run == 0 || codings[run - 1] != coding {
                    cost += OVERHEAD;
                }
                cost += bytes_needed(coding, lengths[run].1);
            }
            cost
        };
        let mut least = u64::MAX;
        for mask in 0..1u32 << runs {
            let codings: Vec<Coding> = (0..runs)
                .map(|run| {
                    if lengths[run].0 || mask >> run & 1 == 1 {
                        Coding::Ncbi4na
                    } else {
                        Coding::Ncbi2na
                    }
                })
                .collect();
            least = least.min(cost(&codings));
        }
        let packed = pack_iupacna(&src);
        let mut end = 0;
        for (index, segment) in packed.iter().enumerate() {
            assert_eq!(segment.start, end);
            assert!(index == 0 || packed[index - 1].coding != segment.coding);
            end += segment.length;
        }
        let mut run_start = 0;
        let packed_codings: Vec<Coding> = lengths
            .iter()
            .map(|&(_, length)| {
                let segment = packed
                    .iter()
                    .find(|segment| run_start < segment.start + segment.length)
                    .unwrap();
                run_start += length;
                assert!(run_start <= segment.start + segment.length);
                segment.coding
            })
            .collect();
        for (run, &(wide, _)) in lengths.iter().enumerate() {
            assert!(!wide || packed_codings[run] == Coding::Ncbi4na);
        }
        assert_eq!(cost(&packed_codings), least, "{lengths:?}");
        assert_eq!(end as usize, src.len());
    }
}

#[test]
fn conversions() {
    assert_eq!(iupacna_to_ncbi2na(b"ACGT"), [0x1b]);
    assert_eq!(iupacna_to_ncbi2na(b"ACGTU"), [0x1b, 0xc0]);
    assert_eq!(iupacna_to_ncbi2na(b"GAT"), [0x8c]);
    assert_eq!(iupacna_to_ncbi2na(b"TTTTTTT"), [0xff, 0xfc]);
    assert_eq!(iupacna_to_ncbi4na(b"G"), [0x40]);
    assert_eq!(iupacna_to_ncbi4na(b"GA"), [0x41]);
    assert_eq!(
        iupacna_to_ncbi4na(b"ACGTUMRWSYKVHDBN"),
        [0x12, 0x48, 0x83, 0x59, 0x6a, 0xc7, 0xbd, 0xef]
    );
    assert_eq!(iupacna_to_ncbi4na(b"ACN"), [0x12, 0xf0]);
}

/// The NCBI4na bytes of a delta's pieces, in either half of a byte where a piece starts.
#[test]
fn ncbi4na_of_pieces() {
    let bytes = |build: &dyn Fn(&mut Ncbi4naBuilder), letters: usize| {
        let mut builder = Ncbi4naBuilder::with_letters(letters).unwrap();
        build(&mut builder);
        assert_eq!(builder.letters(), letters);
        builder.into_bytes()
    };
    assert_eq!(bytes(&|_| {}, 0), [0u8; 0]);
    assert_eq!(
        bytes(
            &|b| {
                b.push_iupacna(b"AC");
                b.push_run(15, 1);
            },
            3
        ),
        [0x12, 0xf0]
    );
    assert_eq!(
        bytes(
            &|b| {
                b.push_iupacna(b"A");
                b.push_iupacna(b"CGT");
                b.push_run(15, 3);
                b.push_iupacna(b"ACG");
                b.push_run(15, 0);
                b.push_run(15, 2);
            },
            12
        ),
        [0x12, 0x48, 0xff, 0xf1, 0x24, 0xff]
    );
    for split in 0..=9 {
        let src = b"ACGTNRYKMS";
        let whole = bytes(&|b| b.push_iupacna(src), 10);
        assert_eq!(whole, iupacna_to_ncbi4na(src));
        let parts = bytes(
            &|b| {
                b.push_iupacna(&src[..split]);
                b.push_iupacna(&src[split..]);
            },
            10,
        );
        assert_eq!(parts, whole, "split {split}");
    }
}

/// Every row of `CIupacnaTo4na::scm_Table` (sequtil_tables.cpp:392-457, read from NCBI's
/// source on 2026-11-12): the rows of these letters in either case, 15 for every other
/// byte.
#[test]
fn ncbi4na_codes_are_the_rows_of_the_table() {
    let rows: [(u8, u8); 15] = [
        (b'A', 1),
        (b'B', 14),
        (b'C', 2),
        (b'D', 13),
        (b'G', 4),
        (b'H', 11),
        (b'K', 12),
        (b'M', 3),
        (b'R', 5),
        (b'S', 6),
        (b'T', 8),
        (b'U', 8),
        (b'V', 7),
        (b'W', 9),
        (b'Y', 10),
    ];
    for byte in 0..=255u8 {
        let expected = rows
            .iter()
            .find(|&&(letter, _)| letter == byte || letter.to_ascii_lowercase() == byte)
            .map_or(15, |&(_, code)| code);
        assert_eq!(ncbi4na_code(byte), expected, "byte {byte:#04x}");
    }
}

/// The NCBI8na codes that `s_FixBioseqDeltas` reads from a literal equal the NCBI4na
/// codes of its letters, whether the literal is NCBI2na (`C2naTo8na`: 0, 1, 2, 3 are
/// 1, 2, 4, 8) or NCBI4na.
#[test]
fn ncbi8na_of_both_codings() {
    for &base in b"ACGTU" {
        assert_eq!(
            1 << ncbi2na_code(base),
            ncbi4na_code(base),
            "{}",
            base as char
        );
    }
    let mut seq4na = Ncbi4naBuilder::with_letters(6).unwrap();
    seq4na.push_iupacna(b"ACGTUN");
    assert_eq!(seq4na.into_bytes(), [0x12, 0x48, 0x8f]);
}

/// The scan's test of the best coding is `best_coding_iupacna` for every byte (the search finds
/// the first byte a test holds for, `utils::byte_search`'s tests).
#[test]
fn the_scan_of_the_best_coding_is_the_table() {
    for byte in 0..=u8::MAX {
        assert_eq!(
            is_ncbi2na_best(byte),
            best_coding_iupacna(byte) == Coding::Ncbi2na,
            "{byte:#04x}"
        );
    }
}

/// `ncbi2na_code` returns the NCBI2na code for upper-case `A`, `C`, `G`, `T` and `U` (0 for
/// every other byte), and `iupacna_to_ncbi2na` packs the codes four to a byte, the first in the
/// high bits, for every length across its short and blocked paths and every byte at every place
/// of a 4-base chunk.
#[test]
fn ncbi2na_codes_and_their_bytes() {
    let row = |residue: u8| match residue {
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 0,
    };
    for byte in 0..=u8::MAX {
        assert_eq!(ncbi2na_code(byte), row(byte), "{byte:#04x}");
    }
    let packed = |src: &[u8]| {
        src.chunks(4)
            .map(|four| {
                four.iter().enumerate().fold(0u8, |byte, (i, &residue)| {
                    byte | row(residue) << (6 - 2 * i)
                })
            })
            .collect::<Vec<u8>>()
    };
    for byte in 0..=u8::MAX {
        for at in 0..4 {
            let mut src = b"ACGT".to_vec();
            src[at] = byte;
            assert_eq!(iupacna_to_ncbi2na(&src), packed(&src), "{byte:#04x} {at}");
        }
    }
    let src: Vec<u8> = (0..1100usize).map(|i| b"ACGTUGCAAT"[i * 7 % 10]).collect();
    for len in (0..=600).chain([1023, 1024, 1025, 1100]) {
        assert_eq!(
            iupacna_to_ncbi2na(&src[..len]),
            packed(&src[..len]),
            "{len}"
        );
    }
}

/// `pack_four_codes` packs every quadruple of 2-bit codes `(c0..=c3)` in `0..4` into
/// `(c0 << 6) | (c1 << 4) | (c2 << 2) | c3`.
#[test]
fn pack_four_codes_matches_all_256_quadruples() {
    for c0 in 0..4u8 {
        for c1 in 0..4u8 {
            for c2 in 0..4u8 {
                for c3 in 0..4u8 {
                    let expected = (c0 << 6) | (c1 << 4) | (c2 << 2) | c3;
                    assert_eq!(
                        pack_four_codes([c0, c1, c2, c3]),
                        expected,
                        "{c0} {c1} {c2} {c3}"
                    );
                }
            }
        }
    }
}
