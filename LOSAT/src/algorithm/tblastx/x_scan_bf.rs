//! EXPERIMENT (LOSAT_X_TBXSCANBF): TBLASTX's lookup scan without a branch on the presence
//! vector.
//!
//! NCBI's scan (c++/src/algo/blast/core/blast_aascan.c:48-131, `s_BlastAaScanSubject`, kept as
//! `blast_aascan.rs` `s_blast_aa_scan_subject`) walks the subject's ranges with
//! `s_DetermineScanningOffsets` (masksubj.inl:43-58) and tests the presence-vector bit of every
//! word of a range. Here each range is scanned by tblastn's branch-free loop
//! (`tblastn/x_scan_bf.rs` `x_scan_one_range_bf`, LOSAT_X_TBNSCANBF; its module text gives how
//! it keeps NCBI's pairs and its stop for room). The ranges are walked as NCBI walks them and the
//! room left after the earlier ranges of a call is passed on, so every call returns NCBI's pairs
//! and `s_range`. The offset array must have `X_SCAN_BF_SLACK` entries after `array_size`
//! (`x_offset_pairs_len`); the entries past the returned count are not part of the result.
//!
//! The branch-free loop gains where the presence bit is set at a hard-to-predict share of the
//! subject words: on the tables of short queries. A table that holds words in half of its cells
//! or more (a translated query of more than about 60 kb: every word of the 25 common letters is
//! then in it) has its bit set at nearly every subject word, the branch is predicted, and the
//! list costs more than it saves; such a table keeps NCBI's loop (`x_scan_bf_for`, decided once
//! per table). Measured on the workstation (perf-si-20261012 session SCANBF-1, one thread, the
//! first L kb of MelaMJNV x PemoMJNVA): time with the branch-free loop over NCBI's loop 0.73 at
//! 2 and 5 kb (presence density 0.16-0.29 of the cells), 0.81 at 10 kb (0.37), 0.91 at 20 kb
//! (0.44), 0.98 at 40 kb (0.49), 1.000 at 80 kb (0.51), 1.004-1.007 from 160 kb (0.53).
//!
//! The switch is read once; it is in the strict set (`utils/xswitch.rs`: on by default, `=0` or
//! `LOSAT_X_DEFAULT=0` turns it off). It acts on the scans of the unit path (`x_pair_par.rs`,
//! `LOSAT_X_PAIRPAR`) and of the reference frame loop (`blast_engine/run_impl.rs`).

use super::blast_aascan::{determine_scanning_offsets, s_blast_aa_scan_subject, BlastOffsetPair};
use super::lookup::BlastAaLookupTable;
use crate::algorithm::tblastn::x_scan_bf::{x_scan_one_range_bf, X_SCAN_BF_SLACK};

// No NCBI counterpart: reads LOSAT_X_TBXSCANBF once (a switch of the strict set).
/// EXPERIMENT (LOSAT_X_TBXSCANBF): whether TBLASTX scans with `x_scan_subject_bf`.
pub(crate) fn x_tbx_scan_bf() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBXSCANBF"))
}

// No NCBI counterpart: the length of LOSAT's offset array (NCBI allocates `array_size` entries,
// lookup_wrap.c GetOffsetArraySize).
/// EXPERIMENT (LOSAT_X_TBXSCANBF): the offset-array length for `array_size` entries, with
/// `X_SCAN_BF_SLACK` more when the switch is on.
pub(crate) fn x_offset_pairs_len(array_size: i32) -> usize {
    let slack = if x_tbx_scan_bf() { X_SCAN_BF_SLACK } else { 0 };
    array_size.max(0) as usize + slack
}

// No NCBI counterpart: chooses between two scans with the same result (module text).
/// EXPERIMENT (LOSAT_X_TBXSCANBF): whether `lookup` is scanned by `x_scan_subject_bf`: the switch is
/// on and fewer than half of the table's backbone cells have their presence bit set.
pub(crate) fn x_scan_bf_for(lookup: &BlastAaLookupTable) -> bool {
    x_tbx_scan_bf() && x_table_sparse(lookup)
}

// No NCBI counterpart (module text).
/// EXPERIMENT (LOSAT_X_TBXSCANBF): fewer than half of the table's backbone cells have their
/// presence bit set.
fn x_table_sparse(lookup: &BlastAaLookupTable) -> bool {
    let set: u64 = lookup.pv.iter().map(|word| word.count_ones() as u64).sum();
    2 * set < lookup.backbone.len() as u64
}

// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:509-511
// ```c
// hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// ```
/// EXPERIMENT (LOSAT_X_TBXSCANBF): `s_blast_aa_scan_subject`, or `x_scan_subject_bf` when
/// `branch_free` (`x_scan_bf_for` of the table; the offset array then from `x_offset_pairs_len`).
#[inline]
pub(crate) fn x_scan_subject(
    branch_free: bool,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    if branch_free {
        x_scan_subject_bf(
            lookup,
            subject,
            seq_ranges,
            offset_pairs,
            array_size,
            s_range,
        )
    } else {
        s_blast_aa_scan_subject(
            lookup,
            subject,
            seq_ranges,
            offset_pairs,
            array_size,
            s_range,
        )
    }
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-130
// ```c
// while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) {
//     s_first=subject->sequence + s_range[1];
//     s_last=subject->sequence + s_range[2];
//     ...
//     for (s = s_first; s <= s_last; s++) {
//         ...
//             if (numhits <= (array_size - totalhits)) {
//                 ...
//             } else {
//                 s_range[1] = (Int4)(s - subject->sequence);
//                 return totalhits;
//             }
//     }
//     s_range[1] = (Int4)(s - subject->sequence);
// }
// return totalhits;
// ```
/// EXPERIMENT (LOSAT_X_TBXSCANBF): `s_blast_aa_scan_subject` with each range scanned by
/// `x_scan_one_range_bf`. Same arguments, same result, same `s_range`; `offset_pairs` must hold
/// `array_size + X_SCAN_BF_SLACK` entries.
pub(crate) fn x_scan_subject_bf(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    let word_length = lookup.word_length;
    let mut totalhits: i32 = 0;
    while determine_scanning_offsets(seq_ranges, word_length, word_length, s_range) {
        let s_last = s_range[2];
        // The room left after the earlier ranges, at the first free entry (blast_aascan.c:95-97:
        // `numhits <= (array_size - totalhits)`, `offset_pairs[i + totalhits]`).
        totalhits += x_scan_one_range_bf(
            lookup,
            subject,
            &mut offset_pairs[totalhits as usize..],
            array_size - totalhits,
            s_range,
        );
        if s_range[1] <= s_last {
            // Stopped for room at word `s_range[1]` (blast_aascan.c:118-123).
            return totalhits;
        }
    }
    totalhits
}

#[cfg(test)]
mod tests {
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
    // s_range) after each call.
    fn calls(
        lookup: &BlastAaLookupTable,
        subject: &[u8],
        seq_ranges: &[(i32, i32)],
        array_size: usize,
        bf: bool,
    ) -> Vec<(Vec<(u32, u32)>, [i32; 3])> {
        let mut pairs = vec![BlastOffsetPair::default(); array_size + X_SCAN_BF_SLACK];
        let wl = lookup.word_length;
        let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - wl];
        let mut out = Vec::new();
        while scan_range[1] <= scan_range[2] {
            let hits = if bf {
                x_scan_subject_bf(
                    lookup,
                    subject,
                    seq_ranges,
                    &mut pairs,
                    array_size as i32,
                    &mut scan_range,
                )
            } else {
                s_blast_aa_scan_subject(
                    lookup,
                    subject,
                    seq_ranges,
                    &mut pairs,
                    array_size as i32,
                    &mut scan_range,
                )
            };
            let got = pairs[..hits as usize]
                .iter()
                .map(|pair| (pair.q_off, pair.s_off))
                .collect();
            out.push((got, scan_range));
            assert!(out.len() < 1_000_000, "scan does not advance");
        }
        out
    }

    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:48-131 (the function compared
    // against, as `s_blast_aa_scan_subject`). Not a port of an NCBI unit case: the pairs and the
    // `s_range` of every call must equal the reference's, for TBLASTX tables (the six frames of
    // nucleotide queries), subjects of one or several ranges (as masked subjects have, ranges of
    // fewer letters than a word among them) and offset arrays that stop the scan after a few words
    // (as aascan_unit_test.cpp's ScanOffsetTest scans).
    #[test]
    fn x_scan_subject_bf_equals_ncbi_scan_on_random_subjects() {
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
            let longest = lookup.longest_chain.max(1) as usize;
            // A translated frame in ncbistdaa: residues 1..=25 (X = 21, stop = 25 included),
            // sentinel 0 after the last one.
            let len = [3, 4, 255, 256, 257, 1_000, 20_000][case % 7];
            let mut subject: Vec<u8> = (0..len).map(|_| 1 + rng.below(25) as u8).collect();
            subject.push(0);
            let len = len as i32;
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
                    let want = calls(&lookup, &subject, ranges, array_size, false);
                    let got = calls(&lookup, &subject, ranges, array_size, true);
                    assert_eq!(
                        got, want,
                        "case {case} ranges {ranges:?} array_size {array_size}"
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
}
