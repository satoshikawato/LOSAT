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
    assert_eq!(ncbi8na_to_ncbi4na(&[1, 2, 15]).unwrap(), [0x12, 0xf0]);
    assert_eq!(ncbi8na_to_ncbi4na(&[]).unwrap(), [0u8; 0]);
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
    let mut seq8na = Vec::new();
    append_ncbi8na(&mut seq8na, b"ACGTUN");
    assert_eq!(seq8na, [1, 2, 4, 8, 8, 15]);
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
