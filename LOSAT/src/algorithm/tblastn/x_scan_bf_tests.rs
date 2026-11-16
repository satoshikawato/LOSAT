//! Tests of the branch-free scan `tblastx/x_aascan_bf.rs` (`x_scan_one_range_bf`) against NCBI's
//! `s_BlastAaScanSubject` (`tblastx/blast_aascan.rs`). Moved from `tblastn/x_scan_bf.rs` (SW6 item 3)
//! with the scan; they build their queries with tblastn's `encode_tblastn_lookup_query`, so they
//! stay under `tblastn`.

use crate::algorithm::tblastn::search_seed::encode_tblastn_lookup_query;
use crate::algorithm::tblastx::blast_aascan::s_blast_aa_scan_subject_one_range;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;
use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;
use crate::algorithm::tblastx::x_aascan_bf::{x_scan_one_range_bf, X_SCAN_BF_SLACK};
use crate::config::ScoringMatrix;
use crate::stats::lookup_protein_params_ungapped;

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

const AMINO: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

fn lookup_for(queries: &[Vec<u8>], threshold: i32, word_length: usize) -> BlastAaLookupTable {
    let frames: Vec<_> = queries
        .iter()
        .map(|query| vec![encode_tblastn_lookup_query(query, None, false)])
        .collect();
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    build_ncbi_lookup_for_profile(
        &frames,
        threshold,
        &karlin,
        false,
        ScoringMatrix::Blosum62,
        word_length,
    )
    .0
}

// The calls of `x_scan_unit` for one range: (offset pairs, s_range[1]) after each call.
fn calls(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    range: (i32, i32),
    array_size: usize,
    bf: bool,
) -> Vec<(Vec<(u32, u32)>, i32)> {
    let mut pairs = vec![BlastOffsetPair::default(); array_size + X_SCAN_BF_SLACK];
    let mut scan_range = [0, range.0, range.1];
    let mut out = Vec::new();
    while scan_range[1] <= scan_range[2] {
        let hits = if bf {
            x_scan_one_range_bf(
                lookup,
                subject,
                &mut pairs,
                array_size as i32,
                &mut scan_range,
            )
        } else {
            s_blast_aa_scan_subject_one_range(
                lookup,
                subject,
                &mut pairs,
                array_size as i32,
                &mut scan_range,
            )
        };
        let got = pairs[..hits as usize]
            .iter()
            .map(|pair| (pair.q_off, pair.s_off))
            .collect();
        out.push((got, scan_range[1]));
        assert!(out.len() < 1_000_000, "scan does not advance");
    }
    out
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127 (the function compared
// against, as `s_blast_aa_scan_subject_one_range`). Not a port of an NCBI unit case: the pairs
// and the `s_range[1]` of every call must equal the reference's, also with offset arrays that
// stop the scan after a few words (as aascan_unit_test.cpp's ScanOffsetTest scans).
#[test]
fn x_scan_bf_equals_ncbi_scan_on_random_subjects() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    let mut compared = 0usize;
    let mut stops = 0usize;
    for case in 0..24 {
        let num_queries = [1, 1, 3, 12][case % 4];
        let queries: Vec<Vec<u8>> = (0..num_queries)
            .map(|_| {
                let len = 20 + rng.below(380);
                (0..len).map(|_| AMINO[rng.below(AMINO.len())]).collect()
            })
            .collect();
        let threshold = [11, 13, 15][case % 3];
        // Word size 3 takes the loop with constant shifts, 2 and 4 the general one.
        let word_length = [3, 3, 3, 3, 2, 4][case % 6];
        let lookup = lookup_for(&queries, threshold, word_length);
        let longest = lookup.longest_chain.max(1) as usize;
        // A translated frame in ncbistdaa: residues 1..=25 (X = 21, stop = 25 included),
        // sentinel 0 after the last one.
        let len = [3, 4, 255, 256, 257, 1_000, 20_000][case % 7];
        let mut subject: Vec<u8> = (0..len).map(|_| 1 + rng.below(25) as u8).collect();
        subject.push(0);
        let wl = lookup.word_length;
        let last = len - wl;
        let ranges = [
            (0, last.max(0)),
            (1.min(last.max(0)), last.max(0)),
            (len / 3, last.max(len / 3)),
        ];
        for &range in &ranges {
            // The pairs of the range in scan order, the same at every array size.
            let mut stream: Option<Vec<(u32, u32)>> = None;
            for array_size in [
                longest,
                longest + 1,
                longest + 2,
                3 * longest,
                // NCBI's offset array (GetOffsetArraySize) and `LOSAT_X_TBNPAIRCAP`'s.
                4096 + longest,
                32768 + longest,
                longest * 1024,
            ] {
                let want = calls(&lookup, &subject, range, array_size, false);
                let got = calls(&lookup, &subject, range, array_size, true);
                assert_eq!(
                    got, want,
                    "case {case} range {range:?} array_size {array_size}"
                );
                let flat: Vec<(u32, u32)> = want
                    .iter()
                    .flat_map(|(pairs, _)| pairs.iter().copied())
                    .collect();
                match &stream {
                    None => stream = Some(flat),
                    Some(first) => assert_eq!(
                        &flat, first,
                        "case {case} range {range:?} array_size {array_size}: order"
                    ),
                }
                compared += want.iter().map(|(pairs, _)| pairs.len()).sum::<usize>();
                stops += want.len().saturating_sub(1);
            }
        }
    }
    // The cases above must reach both the capacity stop and many hits.
    assert!(compared > 100_000, "{compared} pairs compared");
    assert!(stops > 100, "{stops} capacity stops");
}
