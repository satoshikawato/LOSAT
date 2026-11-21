use super::*;
use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
use crate::algorithm::tblastx::translation::generate_frames;
use crate::config::ScoringMatrix;
use crate::stats::lookup_protein_params_ungapped;
use crate::utils::genetic_code::GeneticCode;

struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

// The calls of the reference loop (aa_ungapped.c:500-516) for one subject: (offset pairs,
// s_range) after each call of `scan` (pairs, array_size, s_range -> hits).
fn calls(
    word_length: i32,
    seq_ranges: &[(i32, i32)],
    array_size: usize,
    scan: &mut dyn FnMut(&mut [BlastOffsetPair], i32, &mut [i32; 3]) -> i32,
) -> Vec<(Vec<(u32, u32)>, [i32; 3])> {
    let slack = X_SCAN_BF_SLACK.max(crate::algorithm::tblastx::x_scan_byte::X_SCAN_BYTE_SLACK);
    let mut pairs = vec![BlastOffsetPair::default(); array_size + slack];
    let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - word_length];
    let mut out = Vec::new();
    while scan_range[1] <= scan_range[2] {
        let hits = scan(&mut pairs, array_size as i32, &mut scan_range);
        let got = pairs[..hits as usize]
            .iter()
            .map(|pair| (pair.q_off, pair.s_off))
            .collect();
        out.push((got, scan_range));
        assert!(out.len() < 1_000_000, "scan does not advance");
    }
    out
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:48-131 (the function compared
// against, as `s_blast_aa_scan_subject`). Not a port of an NCBI unit case: the pairs and the
// `s_range` of every call must equal the reference's, for TBLASTX tables (the six frames of
// nucleotide queries), subjects of one or several ranges (as masked subjects have, ranges of
// fewer letters than a word among them) and offset arrays that stop the scan after a few words
// (as aascan_unit_test.cpp's ScanOffsetTest scans).
#[test]
fn x_scan_subject_bf_and_compact_equal_ncbi_scan_on_random_subjects() {
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    let code = GeneticCode::from_id(1);
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let mut compared = 0usize;
    let mut stops = 0usize;
    let mut range_stops = 0usize;
    for case in 0..24 {
        let num_queries = [1, 1, 2, 6][case % 4];
        let queries: Vec<_> = (0..num_queries)
            .map(|_| {
                let len = 30 + rng.below(1_200);
                let dna: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.below(4)]).collect();
                generate_frames(&dna, &code)
            })
            .collect();
        let threshold = [13, 11, 15][case % 3];
        let (lookup, _) = build_ncbi_lookup(&queries, threshold, &karlin, true);
        let (mut compact, _) = build_ncbi_lookup(&queries, threshold, &karlin, true);
        let rank = x_compact(&mut compact);
        assert!(compact.backbone.len() < lookup.backbone.len());
        // LOSAT_X_TBXSCANBYTE: the scan on presence bytes, on any table that `fits`.
        let presence = crate::algorithm::tblastx::x_scan_byte::x_table_presence(&lookup)
            .filter(|_| x_tbx_scan_byte().on())
            .map(|(bytes, _)| bytes);
        let longest = lookup.longest_chain.max(1) as usize;
        // A translated frame in ncbistdaa: residues 1..=25 (X = 21, stop = 25 included),
        // sentinel 0 after the last one.
        let len = [3, 4, 255, 256, 257, 1_000, 20_000][case % 7];
        let mut subject: Vec<u8> = (0..len).map(|_| 1 + rng.below(25) as u8).collect();
        subject.push(0);
        // One range; three ranges, the middle one shorter than a word; ranges every ~300
        // letters with gaps between them.
        let mut cuts = vec![(0, len)];
        let mut three = vec![(0, len / 3), (len / 2, len / 2 + 2), (2 * len / 3, len)];
        three.retain(|&(l, r)| l < r);
        let mut many = Vec::new();
        let mut left = rng.below(20) as i32;
        while left < len {
            let right = (left + 1 + rng.below(600) as i32).min(len);
            many.push((left, right));
            left = right + rng.below(40) as i32;
        }
        for ranges in [&mut cuts, &mut three, &mut many] {
            if ranges.is_empty() {
                continue;
            }
            for array_size in [
                longest,
                longest + 1,
                longest + 2,
                3 * longest,
                4096 + longest,
            ] {
                let wl = lookup.word_length;
                let want = calls(wl, ranges, array_size, &mut |p, a, r| {
                    s_blast_aa_scan_subject(&lookup, &subject, ranges, p, a, r)
                });
                let got = calls(wl, ranges, array_size, &mut |p, a, r| {
                    x_scan_subject_bf(&lookup, &subject, ranges, p, a, r)
                });
                assert_eq!(
                    got, want,
                    "case {case} ranges {ranges:?} array_size {array_size}"
                );
                if let Some(presence) = &presence {
                    let byte_got = calls(wl, ranges, array_size, &mut |p, a, r| {
                        x_scan_subject_byte(&lookup, presence, &subject, ranges, p, a, r)
                    });
                    assert_eq!(
                        byte_got, want,
                        "byte: case {case} ranges {ranges:?} array_size {array_size}"
                    );
                }
                // The compacted table (`x_compact`), with the scan's dispatch and with the
                // portable rank.
                let compact_got = calls(wl, ranges, array_size, &mut |p, a, r| {
                    x_scan_subject_compact(&compact, &rank, &subject, ranges, p, a, r)
                });
                assert_eq!(
                    compact_got, want,
                    "compact: case {case} ranges {ranges:?} array_size {array_size}"
                );
                let portable_got = calls(wl, ranges, array_size, &mut |p, a, r| {
                    x_scan_subject_compact_in(&compact, &rank, &subject, ranges, p, a, r)
                });
                assert_eq!(
                    portable_got, want,
                    "portable rank: case {case} ranges {ranges:?} array_size {array_size}"
                );
                compared += want.iter().map(|(pairs, _)| pairs.len()).sum::<usize>();
                stops += want.len().saturating_sub(1);
                range_stops += want.iter().filter(|(_, r)| r[0] > 0).count();
            }
        }
    }
    // The cases above must reach the capacity stop, later ranges and many hits.
    assert!(compared > 100_000, "{compared} pairs compared");
    assert!(stops > 100, "{stops} capacity stops");
    assert!(
        range_stops > 100,
        "{range_stops} calls ending past the first range"
    );
}

// No NCBI counterpart: the choice of `x_scan_bf_for` (module text). A short query's table has
// its presence bit in fewer than half of the cells; with half of the cells set, NCBI's loop.
#[test]
fn x_table_sparse_counts_the_presence_bits() {
    let mut rng = XorShift(0x0123_4567_89ab_cdef);
    let dna: Vec<u8> = (0..2_000).map(|_| b"ACGT"[rng.below(4)]).collect();
    let code = GeneticCode::from_id(1);
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (mut lookup, _) = build_ncbi_lookup(&[generate_frames(&dna, &code)], 13, &karlin, true);
    assert!(x_table_sparse(&lookup));
    let cells = lookup.backbone.len();
    lookup.pv.fill(0);
    for index in 0..cells / 2 {
        lookup.pv[index >> 5] |= 1 << (index & 31);
    }
    assert!(!x_table_sparse(&lookup), "half of {cells} cells set");
    lookup.pv[(cells / 2 - 1) >> 5] &= !(1 << ((cells / 2 - 1) & 31));
    assert!(
        x_table_sparse(&lookup),
        "one cell fewer than half of {cells}"
    );
}
