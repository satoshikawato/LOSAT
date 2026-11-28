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
                b.push_letters(b"AC");
                b.push_run(15, 1);
            },
            3
        ),
        [0x12, 0xf0]
    );
    assert_eq!(
        bytes(
            &|b| {
                b.push_letters(b"A");
                b.push_letters(b"CGT");
                b.push_run(15, 3);
                b.push_letters(b"ACG");
                b.push_run(15, 0);
                b.push_run(15, 2);
            },
            12
        ),
        [0x12, 0x48, 0xff, 0xf1, 0x24, 0xff]
    );
    for split in 0..=9 {
        let src = b"ACGTNRYKMS";
        let whole = bytes(&|b| b.push_letters(src), 10);
        assert_eq!(whole, iupacna_to_ncbi4na(src));
        let parts = bytes(
            &|b| {
                b.push_letters(&src[..split]);
                b.push_letters(&src[split..]);
            },
            10,
        );
        assert_eq!(parts, whole, "split {split}");
    }
}

/// The codes without a table are the table's for `A`, `C`, `G`, `T`, `U` and `N` in either
/// case, and 0 for every other byte.
#[test]
fn the_codes_of_common_letters_are_the_table() {
    for byte in 0..=u8::MAX {
        let code = common_letter_code(byte);
        if b"ACGTUNacgtun".contains(&byte) {
            assert_eq!(code, NCBI4NA_CODES[usize::from(byte)], "{byte:#04x}");
        } else {
            assert_eq!(code, 0, "{byte:#04x}");
        }
    }
}

/// The residue letters of makeblastdb's nucleotide reader, in either case.
const NUCLEOTIDES: &[u8] = b"ABCDGHKMNRSTUVWYabcdghkmnrstuvwy";

/// The NCBI4na bytes of `push_letters` (blocks by comparisons where every letter is `A`, `C`,
/// `G`, `T`, `U` or `N`, else by the table; the letters after the last block by the table) are those of
/// the table one letter at a time, after an even and after an odd count of letters: every residue
/// letter of the reader (and three other bytes) at every place of the first two blocks (and the
/// letter after them), in lengths that end at that letter, at the end of the second block and
/// after it, every length up to 200, and random letters of the reader in runs of lengths across
/// blocks.
#[test]
fn letters_are_the_table_of_every_letter() {
    let built = |push: &dyn Fn(&mut Ncbi4naBuilder, &[u8]), first: &[u8], src: &[u8]| {
        let mut builder = Ncbi4naBuilder::default();
        builder.push_letters_by_table(first);
        push(&mut builder, src);
        assert_eq!(builder.letters(), first.len() + src.len());
        builder.into_bytes()
    };
    let check = |src: &[u8]| {
        for first in [&b""[..], b"N", b"AC", b"TGA"] {
            assert_eq!(
                built(&|b, s| b.push_letters(s), first, src),
                built(&|b, s| b.push_letters_by_table(s), first, src),
                "{first:?} then {src:?}"
            );
        }
    };
    let letters = b"ACGTUNacgtun";
    let bases: Vec<u8> = (0..200).map(|i| letters[(i * 7 + i / 12) % 12]).collect();
    for place in 0..2 * LETTER_BLOCK + 1 {
        for &byte in NUCLEOTIDES.iter().chain(b"-\x00\xff") {
            let mut src = bases.clone();
            src[place] = byte;
            for length in [place + 1, 2 * LETTER_BLOCK, 2 * LETTER_BLOCK + 3] {
                check(&src[..length.max(place + 1)]);
            }
        }
    }
    for length in 0..=bases.len() {
        check(&bases[..length]);
    }
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = move |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    for _ in 0..200 {
        let length = next(700) as usize;
        let ambiguous = next(4);
        let src: Vec<u8> = (0..length)
            .map(|_| {
                if next(64) < ambiguous {
                    NUCLEOTIDES[next(NUCLEOTIDES.len() as u64) as usize]
                } else {
                    letters[next(12) as usize]
                }
            })
            .collect();
        check(&src);
    }
}

/// A code is one base exactly where its letter's best coding is NCBI2na, for every byte (the
/// classes of `pack_ncbi4na` and `pack_iupacna` agree), and the tests of a byte's two codes are
/// those of each code.
#[test]
fn the_class_of_a_code_is_the_class_of_its_letter() {
    for byte in 0..=u8::MAX {
        assert_eq!(
            is_one_base(ncbi4na_code(byte)),
            is_ncbi2na_best(byte),
            "{byte:#04x}"
        );
    }
    for high in 1..16u8 {
        for low in 1..16u8 {
            let byte = high << 4 | low;
            assert_eq!(
                holds_ambiguous(byte),
                !is_one_base(high) || !is_one_base(low),
                "{byte:#04x}"
            );
            assert_eq!(
                holds_one_base(byte),
                is_one_base(high) || is_one_base(low),
                "{byte:#04x}"
            );
        }
    }
}

/// The NCBI4na codes of `letters`, two to a byte.
fn codes_of(letters: &[u8]) -> Vec<u8> {
    let mut builder = Ncbi4naBuilder::default();
    builder.push_letters_by_table(letters);
    builder.into_bytes()
}

/// The segments of `pack_ncbi4na` over the codes of any span of residues are those of
/// `pack_iupacna` over the span's letters: runs of every class of every length around the
/// search's edges and the block sizes (1 to 5, 19 to 21, 63 to 65, 4,095 to 4,097), boundaries
/// at even and odd places, spans that start and end in either half of a byte, and every length
/// up to 300.
#[test]
fn segments_of_codes_are_those_of_letters() {
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = move |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    const RUNS: [usize; 15] = [1, 2, 3, 4, 5, 19, 20, 21, 63, 64, 65, 4095, 4096, 4097, 600];
    let check = |letters: &[u8], start: usize, end: usize| {
        let codes = codes_of(letters);
        assert_eq!(
            pack_ncbi4na(&codes, start, end),
            pack_iupacna(&letters[start..end]),
            "{start}..{end} of {} letters",
            letters.len()
        );
    };
    for round in 0..400 {
        let mut letters = Vec::new();
        let runs = 1 + next(8);
        let mut wide = next(2) == 1;
        for _ in 0..runs {
            let length = if round % 2 == 0 {
                RUNS[next(RUNS.len() as u64) as usize]
            } else {
                1 + next(40) as usize
            };
            for _ in 0..length {
                letters.push(if wide {
                    b"NRYKMSWBDHVnrykmswbdhv"[next(22) as usize]
                } else {
                    b"ACGTUacgtu"[next(10) as usize]
                });
            }
            wide = !wide;
        }
        let n = letters.len();
        check(&letters, 0, n);
        for _ in 0..4 {
            let start = next(n as u64 + 1) as usize;
            let end = start + next((n - start) as u64 + 1) as usize;
            check(&letters, start, end);
        }
    }
    let src: Vec<u8> = (0..300usize)
        .map(|i| b"ACGTNACGTRACGTACGTUM"[(i * 7 + i / 20) % 20])
        .collect();
    for length in 0..=src.len() {
        check(&src[..length], 0, length);
        if length > 1 {
            check(&src[..length], 1, length);
            check(&src[..length], 1, length - 1);
        }
    }
}

/// The NCBI2na bytes of a raw sequence's codes are those of its letters
/// (`iupacna_to_ncbi2na`), for every length up to 300 and random bases in either case; and each
/// code of one base gives its NCBI2na code at each of the sixteen places of eight bytes.
#[test]
fn ncbi2na_of_codes_is_ncbi2na_of_letters() {
    for (letter, code) in [(b'A', 0u8), (b'C', 1), (b'G', 2), (b'T', 3)] {
        let nibble = ncbi4na_code(letter);
        for place in 0..16 {
            let mut eight = [0u8; 8];
            eight[place / 2] = if place % 2 == 0 { nibble << 4 } else { nibble };
            let mut four = [0u8; 4];
            four[place / 4] = code << (6 - 2 * (place % 4));
            assert_eq!(
                ncbi2na_of_eight(u64::from_le_bytes(eight)).to_le_bytes(),
                four,
                "{letter} {place}"
            );
        }
    }
    let mut state = 0x0fed_cba9_8765_4321_u64;
    let mut next = move |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    for round in 0..4 {
        let src: Vec<u8> = (0..300).map(|_| b"ACGTUacgtu"[next(10) as usize]).collect();
        for length in 0..=src.len() {
            let letters = &src[..length];
            let upper = letters.to_ascii_uppercase();
            assert_eq!(
                ncbi4na_to_ncbi2na(&codes_of(letters), length),
                iupacna_to_ncbi2na(&upper),
                "round {round} length {length}"
            );
        }
    }
}

/// `push_ncbi4na` from every base of a source, after an even and after an odd count of letters
/// (the copy of whole bytes and the shift by one base), for every count, appends the codes of
/// those bases (`push_letters_by_table` of their letters); the count is cut to the source's
/// bases.
#[test]
fn ncbi4na_from_any_base() {
    // An even count: the source has no half byte after its last base.
    let letters: Vec<u8> = (0..78usize)
        .map(|i| b"ACGTNRYKMSWBDHV"[(i * 11 + i / 15) % 15])
        .collect();
    let src = codes_of(&letters);
    for first in [&b""[..], b"N", b"AC", b"TGA"] {
        for start in 0..=letters.len() + 1 {
            for count in 0..=letters.len() + 2 - start.min(letters.len()) {
                let mut built = Ncbi4naBuilder::default();
                built.push_letters_by_table(first);
                built.push_ncbi4na(&src, start, count);
                let mut expected = Ncbi4naBuilder::default();
                expected.push_letters_by_table(first);
                let from = start.min(letters.len());
                expected.push_letters_by_table(&letters[from..(start + count).min(letters.len())]);
                assert_eq!(built, expected, "{first:?} {start} {count}");
            }
        }
    }
}

/// `drop_front` keeps the codes of the last letters, with an even and an odd count dropped and
/// kept.
#[test]
fn dropping_the_first_letters() {
    let letters: Vec<u8> = (0..40usize).map(|i| b"ACGTNRYKM"[i % 9]).collect();
    for length in 0..=letters.len() {
        for count in 0..=length {
            let mut built = Ncbi4naBuilder::default();
            built.push_letters(&letters[..length]);
            built.drop_front(count);
            let mut expected = Ncbi4naBuilder::default();
            expected.push_letters_by_table(&letters[count..length]);
            assert_eq!(built, expected, "{length} {count}");
        }
    }
}

/// The shift and OR of a pair as a 16-bit number is `first << 4 | second` for every pair of
/// bytes.
#[test]
fn pair_bytes_are_the_shift_and_or() {
    for pair in 0..=u16::MAX {
        let [first, second] = pair.to_be_bytes();
        assert_eq!(pair_byte([first, second]), first << 4 | second, "{pair}");
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
    seq4na.push_letters(b"ACGTUN");
    assert_eq!(seq4na.into_bytes(), [0x12, 0x48, 0x8f]);
}

/// The scan's test of the best coding is `best_coding_iupacna` for every byte (the search finds
/// the first byte a test holds for, `utils::byte_search`'s tests).
#[test]
fn the_scan_of_the_best_coding_is_the_table() {
    for byte in 0..=u8::MAX {
        let narrow = best_coding_iupacna(byte) == Coding::Ncbi2na;
        assert_eq!(is_ncbi2na_best(byte), narrow, "{byte:#04x}");
        assert_eq!(
            BEST_IS_NCBI4NA[byte as usize],
            u8::from(!narrow),
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
            let mut short = b"ACGT".to_vec();
            short[at] = byte;
            assert_eq!(
                iupacna_to_ncbi2na(&short),
                packed(&short),
                "short {byte:#04x} {at}"
            );
            let mut blocked = b"ACGT".repeat(10);
            blocked[at] = byte;
            assert_eq!(
                iupacna_to_ncbi2na(&blocked),
                packed(&blocked),
                "blocked {byte:#04x} {at}"
            );
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
