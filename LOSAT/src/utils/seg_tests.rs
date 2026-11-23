use super::*;
use crate::utils::matrix::aa_char_to_ncbistdaa;

fn fasta_sequence_by_id(contents: &str, wanted_id: &str) -> Vec<u8> {
    let mut current_id = None;
    let mut residues = Vec::new();
    for line in contents.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            let found = rest.split_whitespace().next();
            if current_id == Some(wanted_id) {
                break;
            }
            current_id = found;
            continue;
        }
        if current_id == Some(wanted_id) {
            residues.extend(
                line.bytes()
                    .map(|aa| aa_char_to_ncbistdaa(aa.to_ascii_uppercase())),
            );
        }
    }
    residues
}

// The histogram path (was LOSAT_X_SEGFAST) against the reference on
// random protein-like sequences with planted low-complexity stretches,
// ambiguity codes and stop codons.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2030-2032,2053-2060
// The test masks random sequences with both implementations and compares the intervals.
#[test]
fn x_fast_seg_matches_reference_on_random_sequences() {
    let mut state = 0x243F_6A88_85A3_08D3u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let valid: Vec<u8> = (0u8..28)
        .filter(|&c| c == 1 || (3..=20).contains(&c) || c == 22)
        .collect();
    let cases: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(4000);
    let mut masked = 0usize;
    for case in 0..cases {
        let params = match case % 5 {
            0 => SegParams::new(12, 2.2, 2.5),
            1 => SegParams::new(10, 1.8, 2.1),
            2 => SegParams::with_all(12, 2.2, 2.5, 2, 50),
            3 => SegParams::with_all(15, 2.6, 2.9, 15, 100),
            _ => SegParams::with_all(8, 1.5, 1.9, 1, 10),
        };
        let masker = SegMasker::with_params(&params);
        assert!(masker.x_fast_applies());
        let len = 1 + (next() % if case % 11 == 0 { 3000 } else { 400 }) as usize;
        let mut seq = Vec::with_capacity(len);
        while seq.len() < len {
            match next() % 10 {
                // a run over a tiny alphabet (possibly a homopolymer)
                0 | 1 => {
                    let letters = 1 + (next() % 4) as usize;
                    let pool: Vec<u8> = (0..letters)
                        .map(|_| valid[(next() % valid.len() as u64) as usize])
                        .collect();
                    let run = 1 + (next() % if case % 13 == 0 { 300 } else { 60 }) as usize;
                    for _ in 0..run {
                        seq.push(pool[(next() % pool.len() as u64) as usize]);
                    }
                }
                // letters outside the SEG alphabet (X, B, Z, *, gap, ...)
                2 => {
                    let run = 1 + (next() % 6) as usize;
                    for _ in 0..run {
                        seq.push([0u8, 2, 21, 23, 24, 25, 26, 27, 200][(next() % 9) as usize]);
                    }
                }
                _ => {
                    let run = 1 + (next() % 40) as usize;
                    for _ in 0..run {
                        seq.push(valid[(next() % valid.len() as u64) as usize]);
                    }
                }
            }
        }
        seq.truncate(len);
        let reference = if seq.len() < masker.window {
            Vec::new()
        } else {
            masker.mask_sequence_reference::<true>(&seq)
        };
        let fast = if seq.len() < masker.window {
            Vec::new()
        } else {
            masker.x_mask_sequence::<true>(&seq)
        };
        assert_eq!(fast, reference, "case {case}, length {len}");
        // segmasker's SEG, without the merge (review of line RG item K2, N1).
        if seq.len() >= masker.window {
            assert_eq!(
                masker.x_mask_sequence::<false>(&seq),
                masker.mask_sequence_reference::<false>(&seq),
                "case {case}, length {len}, unmerged"
            );
        }
        masked += reference.len();
    }
    assert!(
        masked > cases / 4,
        "only {masked} intervals in {cases} cases"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1518,1522,1540-1541
// LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26): TBLASTN's redo masks overlapping windows of one subject;
// the windows below overlap so that their segments repeat and the trim memo answers. Every
// window's intervals must equal the reference SEG's.
#[test]
fn x_trim_memo_on_overlapping_windows_matches_reference() {
    let mut state = 0x1357_9BDF_2468_ACE0u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let valid: Vec<u8> = (0u8..28)
        .filter(|&c| c == 1 || (3..=20).contains(&c) || c == 22)
        .collect();
    let mut seq = Vec::new();
    while seq.len() < 6000 {
        if next() % 4 == 0 {
            let pool: Vec<u8> = (0..1 + next() % 3)
                .map(|_| valid[(next() % valid.len() as u64) as usize])
                .collect();
            for _ in 0..8 + next() % 40 {
                seq.push(pool[(next() % pool.len() as u64) as usize]);
            }
        } else {
            for _ in 0..1 + next() % 60 {
                seq.push(valid[(next() % valid.len() as u64) as usize]);
            }
        }
    }
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1)).x_with_trim_memo();
    let hits_before = crate::utils::x_seg_trim_memo::x_hits();
    let mut masked = 0usize;
    for _ in 0..300 {
        let begin = (next() % 5000) as usize;
        let end = (begin + 50 + (next() % 900) as usize).min(seq.len());
        let window = &seq[begin..end];
        let reference = masker.mask_sequence_reference::<true>(window);
        assert_eq!(
            masker.x_mask_sequence::<true>(window),
            reference,
            "window {begin}..{end}"
        );
        assert_eq!(
            masker.x_mask_sequence::<false>(window),
            masker.mask_sequence_reference::<false>(window),
            "window {begin}..{end}, unmerged"
        );
        masked += reference.len();
    }
    assert!(masked > 300, "only {masked} intervals");
    assert!(
        crate::utils::x_seg_trim_memo::x_hits() > hits_before + 100,
        "the windows should repeat segments"
    );
}

#[test]
fn test_seg_params_default() {
    let params = SegParams::default();
    assert_eq!(params.window, 12);
    assert_eq!(params.locut, 2.2);
    assert_eq!(params.hicut, 2.5);
}

#[test]
fn test_seg_params_validation() {
    // NCBI: negative values become 0.0, not defaults
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2253-2255
    // (a `locut` above `hicut` raises `hicut` to it)
    let params = SegParams::new(0, -1.0, 1.0);
    assert_eq!(params.window, 12); // Should default to 12
    assert_eq!(params.locut, 0.0); // Negative becomes 0.0
    assert_eq!(params.hicut, 1.0); // 1.0 > 0.0, so stays at 1.0
}

#[test]
fn test_calculate_entropy() {
    let masker = SegMasker::with_defaults();

    // Low complexity: all same amino acid (entropy = 0)
    // NCBISTDAA: A=1 (not 0, which is '-')
    let low_complex = vec![1u8; 12]; // All alanine (NCBISTDAA code 1)
    let (entropy_low, _) = masker.calculate_entropy(&low_complex);
    assert!(
        entropy_low < 1.0,
        "Low complexity should have low entropy: got {}",
        entropy_low
    );

    // High complexity: all different amino acids
    // Use valid NCBISTDAA codes: 1=A, 3-20=C..W, 22=Y
    let high_complex: Vec<u8> = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13].to_vec();
    let (entropy_high, _) = masker.calculate_entropy(&high_complex);
    assert!(
        entropy_high > 2.0,
        "High complexity should have high entropy: got {}",
        entropy_high
    );
}

#[test]
fn test_mask_simple_repeat() {
    let masker = SegMasker::with_defaults();

    // Simple repeat sequence should be masked
    // NCBISTDAA: A=1 (not 0, which is '-')
    let seq: Vec<u8> = vec![1u8; 100]; // All alanine (NCBISTDAA code 1)
    let intervals = masker.mask_sequence(&seq);

    assert!(
        !intervals.is_empty(),
        "Poly-alanine sequence should be masked"
    );
    if !intervals.is_empty() {
        assert!(intervals[0].start < intervals[0].end);
    }
}

#[test]
fn test_mask_short_sequence() {
    let masker = SegMasker::with_defaults();

    // Very short sequences should return empty (shorter than window)
    // NCBISTDAA: A=1
    let seq = vec![1u8; 5];
    let intervals = masker.mask_sequence(&seq);
    assert!(intervals.is_empty());
}

#[test]
fn test_mask_complex_sequence() {
    let masker = SegMasker::with_defaults();

    // Complex sequence should not be heavily masked
    let seq: Vec<u8> = (0..100).map(|i| (i % 20) as u8).collect();
    let intervals = masker.mask_sequence(&seq);

    // Complex sequence should have few or no masked regions
    let total_masked: usize = intervals.iter().map(|i| i.end - i.start).sum();
    assert!(
        total_masked < seq.len() / 2,
        "Complex sequence should not be heavily masked"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1413-1416
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_bdt63528() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/PajaWSV.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "BDT63528.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(217, 223),
            MaskedInterval::new(568, 580),
            MaskedInterval::new(762, 774),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1413-1416
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_bdt63573() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/PajaWSV.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "BDT63573.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(146, 176),
            MaskedInterval::new(192, 209),
            MaskedInterval::new(214, 227),
            MaskedInterval::new(408, 420),
            MaskedInterval::new(587, 601),
            MaskedInterval::new(605, 619),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1413-1416
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_bdt63533() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/PajaWSV.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "BDT63533.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(293, 303),
            MaskedInterval::new(412, 427),
            MaskedInterval::new(873, 882),
            MaskedInterval::new(1077, 1087),
        ]
    );
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in /tmp/BDT63510.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1415-1438
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_bdt63510() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/PajaWSV.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "BDT63510.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(207, 228),
            MaskedInterval::new(257, 266),
            MaskedInterval::new(591, 610),
            MaskedInterval::new(789, 797),
        ]
    );
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in - -outfmt interval < tests/fasta/WSSV.faa`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1381-1382
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_yp_009220537_1() {
    let fasta = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fasta/WSSV.faa"));
    let seq = fasta_sequence_by_id(fasta, "YP_009220537.1");
    let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![MaskedInterval::new(319, 335), MaskedInterval::new(920, 931),]
    );
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in - -outfmt interval < tests/fasta/WSSV.faa`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1381-1382
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_yp_009220568_1() {
    let fasta = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fasta/WSSV.faa"));
    let seq = fasta_sequence_by_id(fasta, "YP_009220568.1");
    let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(12, 28),
            MaskedInterval::new(94, 129),
            MaskedInterval::new(280, 293),
            MaskedInterval::new(452, 536),
            MaskedInterval::new(649, 655),
            MaskedInterval::new(709, 727),
            MaskedInterval::new(1041, 1054),
            MaskedInterval::new(1139, 1149),
        ]
    );
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in - -outfmt interval < tests/fasta/WSSV.faa`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1381-1382
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_yp_009220574_1() {
    let fasta = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fasta/WSSV.faa"));
    let seq = fasta_sequence_by_id(fasta, "YP_009220574.1");
    let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(32, 39),
            MaskedInterval::new(364, 382),
            MaskedInterval::new(774, 791),
            MaskedInterval::new(827, 839),
            MaskedInterval::new(915, 923),
            MaskedInterval::new(1103, 1112),
            MaskedInterval::new(1227, 1246),
            MaskedInterval::new(1267, 1279),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1140-1157
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2125-2149
// Independent `SeqBufferSeg` API output with this actual caller's `overlaps` enabled.
// The former assertion combined adjacent locations, unlike this core caller.
#[test]
fn test_mask_sequence_matches_ncbi_core_bdv02435_1_query_defaults() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/AP027131.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "BDV02435.1");
    let masker = SegMasker::with_params(&SegParams::new(12, 2.2, 2.5));
    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);
    let expected: Vec<_> =
        include_str!("../../tests/unit/blastx_stage_e_seg_sequence_expected.tsv")
            .lines()
            .map(|line| {
                let (begin, end) = line.split_once('\t').unwrap();
                MaskedInterval::new(begin.parse().unwrap(), end.parse().unwrap())
            })
            .collect();
    assert_eq!(expected.len(), 27);
    assert_eq!(intervals, expected);
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in /tmp/BDT63134.faa -outfmt interval`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1415-1438
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_bdt63134() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/SicyWSV.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "BDT63134.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(251, 264),
            MaskedInterval::new(311, 321),
            MaskedInterval::new(348, 363),
            MaskedInterval::new(382, 402),
        ]
    );
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in /tmp/ap_subjects.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_wp_025208654_1() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/NZ_CP006932.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "WP_025208654.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(intervals, vec![MaskedInterval::new(137, 156),]);
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in /tmp/ap_subjects.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_wp_025208655_1() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/NZ_CP006932.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "WP_025208655.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(148, 167),
            MaskedInterval::new(1270, 1280),
        ]
    );
}

// Output of NCBI's segmasker 2.17.0:
// Command:
// `segmasker -in /tmp/ap_targets.faa -outfmt interval -window 10 -locut 1.8 -hicut 2.1`
//
// `segmasker` interval output is `0-based closed`; Rust stores merged
// `MaskedInterval` values as `0-based half-open`.
#[test]
fn test_mask_sequence_matches_ncbi_segmasker_bdv02435_1() {
    let fasta = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fasta/AP027131.faa"
    ));
    let seq = fasta_sequence_by_id(fasta, "BDV02435.1");
    let masker = SegMasker::with_params(&SegParams::new(10, 1.8, 2.1));

    let mut intervals = masker.mask_sequence(&seq);
    intervals.sort_by_key(|interval| interval.start);

    assert_eq!(
        intervals,
        vec![
            MaskedInterval::new(8, 19),
            MaskedInterval::new(174, 183),
            MaskedInterval::new(629, 640),
            MaskedInterval::new(697, 708),
            MaskedInterval::new(908, 918),
            MaskedInterval::new(1273, 1287),
            MaskedInterval::new(1477, 1491),
            MaskedInterval::new(1946, 1956),
            MaskedInterval::new(4091, 4105),
            MaskedInterval::new(4208, 4222),
            MaskedInterval::new(4241, 4256),
            MaskedInterval::new(4347, 4365),
            MaskedInterval::new(4568, 4578),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2125-2149
#[test]
fn test_merge_segs_inclusive_matches_ncbi_overlap_merge() {
    let mut segs = vec![(20, 30), (10, 25), (0, 5)];
    merge_segs_inclusive(64, &mut segs);
    assert_eq!(segs, vec![(10, 30), (0, 5)]);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:2140-2151
// Expected survivors come from the verbatim pinned C function, recorded in
// Session E/oracles/seg_overlap; the Rust helper never generates expected data.
#[test]
fn test_merge_segs_independent_ncbi_boundary_oracle() {
    let parse = |text: &str| -> Vec<(usize, usize)> {
        if text == "-" {
            return Vec::new();
        }
        text.split(',')
            .map(|range| {
                let (begin, end) = range.split_once(':').unwrap();
                (begin.parse().unwrap(), end.parse().unwrap())
            })
            .collect()
    };
    for row in include_str!("../../tests/unit/blastx_stage_e_seg_expected.tsv").lines() {
        let fields: Vec<_> = row.split('\t').collect();
        let mut intervals = parse(fields[2]);
        merge_segs_inclusive(fields[1].parse().unwrap(), &mut intervals);
        assert_eq!(intervals, parse(fields[3]), "{}", fields[0]);
    }
}

// No NCBI counterpart (reuse only): the two process-wide tables hold every window of at most
// 15 letters (684 count lists), and each gives the port's bits for it.
#[test]
fn x_entropy_tables_hold_every_short_window_with_the_ports_bits() {
    let mut lists = 0;
    x_each_count_list(&mut |counts| {
        lists += 1;
        let reference = entropy_from_state_vector_reference(counts).to_bits();
        let state = counts
            .iter()
            .enumerate()
            .fold(0u64, |key, (class, &count)| {
                key | ((count as u64) << (4 * class))
            });
        assert_eq!(
            x_memo_get(x_state_entropies(), state).map(f64::to_bits),
            Some(reference),
            "{counts:?}"
        );
        assert_eq!(entropy_from_state_vector(counts).to_bits(), reference);
        let histogram = counts
            .iter()
            .fold(0u64, |key, &count| key + (1u64 << (4 * (count - 1))));
        assert_eq!(
            x_memo_get(x_histogram_entropies(), histogram).map(f64::to_bits),
            Some(reference),
            "{counts:?}"
        );
    });
    assert_eq!(lists, 684);
}
