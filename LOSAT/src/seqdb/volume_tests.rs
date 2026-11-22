use super::{overlay_ambiguity_runs, AmbiguityRun, NucleotideAmbiguity};

/// A small deterministic generator (xorshift), so the test needs no fixture and no crate.
struct Rng(u64);
impl Rng {
    fn next(&mut self, below: u64) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % below
    }
}

/// The letters of the whole sequence with every run applied in stored order (the decode
/// that `SequenceBlock` does), the reference for any range.
fn whole_decode(length: usize, runs: &[AmbiguityRun]) -> Vec<u8> {
    let mut letters = vec![b'A'; length];
    overlay_ambiguity_runs(runs, false, 0, length, &mut letters);
    letters
}

fn check(length: usize, runs: Vec<AmbiguityRun>, rng: &mut Rng, ranges: usize) -> bool {
    let whole = whole_decode(length, &runs);
    let ambiguity = NucleotideAmbiguity::new(length, 0, runs);
    for _ in 0..ranges {
        let from = rng.next(length as u64 + 1) as usize;
        let to = from + rng.next((length - from) as u64 + 1) as usize;
        let mut letters = vec![b'A'; to - from];
        ambiguity.overlay(from, to, &mut letters);
        assert_eq!(letters, whole[from..to], "range {from}..{to}");
    }
    ambiguity.ordered
}

/// The ambiguity loop of `nucleotide_blastna_in` over the whole sequence.
fn whole_runs(runs: &[AmbiguityRun], seq: &mut [u8]) {
    let region_end = seq.len() as i64;
    for run in runs {
        let mut position = i64::from(run.offset as i32);
        let row_len = i64::from(run.length);
        if position + row_len <= 0 {
            continue;
        }
        if position >= region_end {
            break;
        }
        for _ in 0..row_len {
            if position >= 0 && position < region_end {
                seq[position as usize] = run.ncbi4na;
            }
            position += 1;
        }
    }
}

// LOSAT_X_BNTBRANGE: the runs over ranges give the whole decode's bytes inside the ranges,
// for ordered runs, shuffled and overlapping ones, offsets of 2^31 and more (negative as
// NCBI's Int4), and runs past the end (which stop the whole decode's loop).
#[test]
fn runs_over_ranges_equal_the_whole_decode_inside_them() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    for case in 0..400 {
        let length = 1 + rng.next(3000) as usize;
        let mut runs = Vec::new();
        for _ in 0..rng.next(40) {
            let offset = match rng.next(20) {
                0 => 0x8000_0000u32 + rng.next(1 << 20) as u32,
                1 => length as u32 + rng.next(100) as u32,
                2 => u32::MAX - rng.next(50) as u32,
                _ => rng.next(length as u64) as u32,
            };
            let longest = if rng.next(10) == 0 { 5000 } else { 20 };
            runs.push(AmbiguityRun {
                offset,
                length: 1 + rng.next(longest) as u32,
                ncbi4na: rng.next(16) as u8,
            });
        }
        if case % 2 == 0 {
            runs.sort_by_key(|run| run.offset);
        }
        let mut ranges = Vec::new();
        let mut at = 0usize;
        while at < length {
            let from = at + rng.next(200) as usize;
            let to = (from + rng.next(300) as usize).min(length);
            if from < to {
                ranges.push((from, to));
            }
            at = to + 1 + rng.next(100) as usize;
        }
        let mut whole = vec![1u8; length];
        whole_runs(&runs, &mut whole);
        let mut ranged = vec![1u8; length];
        super::x_overlay_runs_in_ranges(&runs, &ranges, &mut ranged, &super::X_NCBI4NA);
        for &(from, to) in &ranges {
            assert_eq!(
                ranged[from..to],
                whole[from..to],
                "case {case}, range {from}..{to}"
            );
        }
    }
}

#[test]
fn many_runs_and_many_ranges_equal_the_whole_decode() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    // 3000 ordered runs of 1..=20 bases with gaps of 0..=30 bases (NCBI's layout).
    let mut runs = Vec::new();
    let mut at = 5u32;
    for _ in 0..3000 {
        let length = 1 + rng.next(20) as u32;
        runs.push(AmbiguityRun {
            offset: at,
            length,
            ncbi4na: 2 + rng.next(14) as u8,
        });
        at += length + rng.next(31) as u32;
    }
    let length = at as usize + 10;
    assert!(check(length, runs.clone(), &mut rng, 5000));
    // Ranges that start inside a run, end inside one, or lie between two.
    for run in runs.iter().step_by(97) {
        let mut one = [b'A'; 2];
        let ambiguity = NucleotideAmbiguity::new(length, 0, runs.clone());
        ambiguity.overlay(run.offset as usize, run.offset as usize + 1, &mut one[..1]);
        assert_ne!(one[0], b'A');
    }
    // Stored runs that overlap or are out of order: every run applies in stored order, so
    // a later run still overwrites an earlier one.
    let mut messy = runs.clone();
    messy.swap(10, 2900);
    messy.push(AmbiguityRun {
        offset: 7,
        length: 5000,
        ncbi4na: 3,
    });
    messy.push(AmbiguityRun {
        offset: 100,
        length: 3,
        ncbi4na: 9,
    });
    assert!(!check(length, messy, &mut rng, 5000));
}

#[test]
fn a_range_outside_every_run_stays_unchanged() {
    let runs = vec![AmbiguityRun {
        offset: 10,
        length: 4,
        ncbi4na: 15,
    }];
    let ambiguity = NucleotideAmbiguity::new(40, 0, runs);
    for (from, to) in [(0, 10), (14, 40), (0, 0), (40, 40)] {
        let mut letters = vec![b'C'; to - from];
        ambiguity.overlay(from, to, &mut letters);
        assert!(letters.iter().all(|&letter| letter == b'C'), "{from}..{to}");
    }
    let mut letters = vec![b'C'; 6];
    ambiguity.overlay(8, 14, &mut letters);
    assert_eq!(letters, b"CCNNNN");
}

// LOSAT_X_BNDECODE: the by-byte decode gives the bytes of the decode `=0` runs, for the whole
// sequence and for ranges that start and end inside a stored byte, with every stored byte
// value, ordered and shuffled runs, runs of each code, offsets of 2^31 and more (negative as
// NCBI's Int4) and runs past the end.
#[test]
fn the_by_byte_decode_equals_the_by_base_decode() {
    let mut rng = Rng(0x0DDB_1A5E_5BAD_5EED);
    for case in 0..600 {
        let mut length = 1 + rng.next(if case % 3 == 0 { 9 } else { 4000 }) as usize;
        let mut packed: Vec<u8> = (0..length.div_ceil(4))
            .map(|_| rng.next(256) as u8)
            .collect();
        if case == 0 {
            // Every stored byte value once.
            packed = (0..=255u8).collect();
            length = packed.len() * 4;
        }
        let mut runs = Vec::new();
        for _ in 0..rng.next(30) {
            let offset = match rng.next(20) {
                0 => 0x8000_0000u32 + rng.next(1 << 20) as u32,
                1 => length as u32 + rng.next(100) as u32,
                2 => u32::MAX - rng.next(50) as u32,
                _ => rng.next(length as u64) as u32,
            };
            let longest = if rng.next(10) == 0 { 5000 } else { 20 };
            runs.push(AmbiguityRun {
                offset,
                length: 1 + rng.next(longest) as u32,
                ncbi4na: rng.next(16) as u8,
            });
        }
        if case % 2 == 0 {
            runs.sort_by_key(|run| run.offset);
        }
        let mut by_base = vec![0u8; length];
        super::x_blastna_by_base(&packed, &runs, &mut by_base);
        let mut by_byte = vec![0u8; length];
        super::x_blastna_by_byte(&packed, &runs, &[(0, length)], &mut by_byte);
        assert_eq!(by_byte, by_base, "case {case}, whole");
        let mut ranges = Vec::new();
        let mut at = rng.next(7) as usize;
        while at < length {
            let to = (at + 1 + rng.next(if case % 2 == 0 { 9 } else { 300 }) as usize).min(length);
            ranges.push((at, to));
            at = to + rng.next(12) as usize;
        }
        let mut by_base = vec![201u8; length];
        super::x_blastna_by_base_ranges(&packed, &runs, &ranges, &mut by_base);
        let mut by_byte = vec![201u8; length];
        super::x_blastna_by_byte(&packed, &runs, &ranges, &mut by_byte);
        assert_eq!(by_byte, by_base, "case {case}, ranges {ranges:?}");
    }
}
