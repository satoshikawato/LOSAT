//! EXPERIMENT (LOSAT_X_TBNSCANBF): tblastn's lookup scan without a branch on the presence
//! vector.
//!
//! NCBI's scan (c++/src/algo/blast/core/blast_aascan.c:75-127, `s_BlastAaScanSubject`, kept
//! as `tblastx/blast_aascan.rs` `s_blast_aa_scan_subject_one_range`) tests the presence-vector
//! bit of every subject word and branches on it. On a genome searched with one protein the bit
//! is set at a hard-to-predict share of the words, and both NCBI and LOSAT spend about a quarter
//! of their cycles on mispredicted branches (perf-si-20261012 session 19, decision 169).
//!
//! Here the scan takes the subject in blocks of `X_BLOCK` words. For each word of a block it
//! writes the word's position and table index to a candidate list and advances the list's
//! length by the presence bit, with no branch; then it copies the hits of the listed words in
//! position order, exactly as NCBI's loop does, including the stop when the offset array has no
//! room for a word's hits: `s_range[1]` is then that word's position, so the next call resumes
//! there. The words after it in the block were looked at but not used; the next call looks at
//! them again. The offset pairs and the `s_range` values of every call are NCBI's.
//!
//! The hits of a backbone cell (at most `AA_HITS_PER_CELL`) are copied as `AA_HITS_PER_CELL`
//! pairs and the count advanced by the cell's count, so the offset array must have
//! `X_SCAN_BF_SLACK` entries after `array_size`; the entries past the returned count are not
//! part of the result.
//!
//! The switch is read once, off unless set (not in the strict set). blastp keeps the scan of
//! `tblastx/blast_aascan.rs`; blastx calls this scan and the candidate loop (`x_block_candidates`)
//! under its own switch, LOSAT_X_BXSCANBF (`blastx/x_seed_scan.rs`); TBLASTX scans each subject
//! range with this loop under its own switch, LOSAT_X_TBXSCANBF (`tblastx/x_scan_bf.rs`).

use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;
use crate::algorithm::tblastx::lookup::{
    BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use std::mem::MaybeUninit;

/// EXPERIMENT (LOSAT_X_TBNSCANBF): subject words looked at before their hits are copied.
pub(crate) const X_BLOCK: usize = 256;

/// EXPERIMENT (LOSAT_X_TBNSCANBF): offset-array entries needed after `array_size` (a cell's
/// hits are written as `AA_HITS_PER_CELL` pairs, also when the array is full).
pub(crate) const X_SCAN_BF_SLACK: usize = AA_HITS_PER_CELL;

// No NCBI counterpart: reads LOSAT_X_TBNSCANBF once; off unless set.
/// EXPERIMENT (LOSAT_X_TBNSCANBF): whether tblastn scans with `x_scan_one_range_bf`.
pub(super) fn x_tbn_scan_bf() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_TBNSCANBF").is_some())
}

// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:75-127
// ```c
// index = ComputeTableIndex(word_length - 1, lookup->charsize, s_first);
// for (s = s_first; s <= s_last; s++) {
//     index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
//     if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
//         numhits = bbc[index].num_used;
//         if (numhits <= (array_size - totalhits)) {
//             ...
//             for (i = 0; i < numhits; i++) {
//                 offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
//                 offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
//             }
//             totalhits += numhits;
//         } else {
//             s_range[1] = (Int4)(s - subject->sequence);
//             return totalhits;
//         }
//     }
// }
// s_range[1] = (Int4)(s - subject->sequence);
// ```
/// EXPERIMENT (LOSAT_X_TBNSCANBF): `s_blast_aa_scan_subject_one_range` with the presence test
/// made branch-free (see the module text). Same arguments, same result, same `s_range`;
/// `offset_pairs` must hold `array_size + X_SCAN_BF_SLACK` entries.
pub(crate) fn x_scan_one_range_bf(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    if s_range[1] > s_range[2] {
        return 0;
    }
    assert!(
        array_size >= 0 && offset_pairs.len() >= array_size as usize + X_SCAN_BF_SLACK,
        "offset array of {} entries for array_size {array_size}",
        offset_pairs.len()
    );

    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let pv = lookup.pv.as_ptr();
    let backbone = lookup.backbone.as_ptr();
    let overflow = lookup.overflow.as_ptr();
    let dest = offset_pairs.as_mut_ptr();

    let s_first = s_range[1] as usize;
    let s_last = s_range[2] as usize;

    // NCBI reference: c++/include/algo/blast/core/blast_lookup.h:96-107 (ComputeTableIndex)
    let mut index: usize = 0;
    for i in 0..(word_length - 1) {
        // SAFETY: as in `s_blast_aa_scan_subject_one_range`: the scan range follows NCBI's
        // sentinel-backed subject contract.
        let ch = unsafe { *subject.get_unchecked(s_first + i) } as usize;
        index = (index << charsize) | ch;
    }

    // Each candidate is `position << 32 | index`.
    let mut cand = [MaybeUninit::<u64>::uninit(); X_BLOCK];
    let mut totalhits: i32 = 0;
    let mut block_first = s_first;
    while block_first <= s_last {
        let block_last = s_last.min(block_first + X_BLOCK - 1);
        let n = x_block_candidates(
            subject,
            pv,
            (word_length, charsize, mask),
            (block_first, block_last),
            &mut index,
            &mut cand,
        );
        for c in &cand[..n] {
            // SAFETY: the first `n` candidates were written above.
            let c = unsafe { c.assume_init() };
            let s_off = (c >> 32) as u32;
            // SAFETY: a masked index, as above.
            let cell = unsafe { &*backbone.add((c & 0xffff_ffff) as usize) };
            let numhits = cell.num_used;
            if numhits > array_size - totalhits {
                // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:118-123
                s_range[1] = s_off as i32;
                return totalhits;
            }
            // SAFETY: `totalhits + numhits <= array_size` and the array has `X_SCAN_BF_SLACK`
            // entries after `array_size`, so both copies stay inside it.
            unsafe {
                let d = dest.add(totalhits as usize);
                if numhits as usize <= AA_HITS_PER_CELL {
                    for i in 0..AA_HITS_PER_CELL {
                        *d.add(i) = BlastOffsetPair {
                            q_off: cell.entries[i] as u32,
                            s_off,
                        };
                    }
                } else {
                    let src = overflow.add(cell.entries[0] as usize);
                    for i in 0..numhits as usize {
                        *d.add(i) = BlastOffsetPair {
                            q_off: *src.add(i) as u32,
                            s_off,
                        };
                    }
                }
            }
            totalhits += numhits;
        }
        block_first = block_last + 1;
    }

    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:125-127
    s_range[1] = (s_last + 1) as i32;
    totalhits
}

/// EXPERIMENT (LOSAT_X_TBNSCANBF): writes the candidates of the words `block.0 ..= block.1`
/// (`position << 32 | index`, the list's length advanced by the presence bit) and returns how
/// many; `index` is the rolling index of the word before `block.0` and becomes that of
/// `block.1`. `table` = (word size, letter width, mask); `pv` is the table's presence vector.
/// The words read are those of NCBI's loop over the same positions (`subject[s + word size - 1]`
/// for every `s` of the block), which the caller's subject must hold.
#[inline(always)]
pub(crate) fn x_block_candidates(
    subject: &[u8],
    pv: *const u32,
    table: (usize, usize, usize),
    block: (usize, usize),
    index: &mut usize,
    cand: &mut [MaybeUninit<u64>; X_BLOCK],
) -> usize {
    // tblastn's and blastx's table (word size 3, 5-bit letters) with the shifts as constants.
    if table.0 == 3 && table.1 == 5 {
        x_candidates::<3, 5>(subject, pv, table, block, index, cand)
    } else {
        x_candidates::<0, 0>(subject, pv, table, block, index, cand)
    }
}

/// EXPERIMENT (LOSAT_X_TBNSCANBF): [`x_block_candidates`] for one table shape. `WORD` and
/// `CHARSIZE` are the table's word size and letter width when known at compile time, 0 to read
/// them from `table` (word size, letter width, mask).
#[inline(always)]
fn x_candidates<const WORD: usize, const CHARSIZE: usize>(
    subject: &[u8],
    pv: *const u32,
    table: (usize, usize, usize),
    block: (usize, usize),
    index: &mut usize,
    cand: &mut [MaybeUninit<u64>; X_BLOCK],
) -> usize {
    let word_length = if WORD != 0 { WORD } else { table.0 };
    let charsize = if CHARSIZE != 0 { CHARSIZE } else { table.1 };
    let mask = table.2;
    let mut ix = *index;
    let mut n = 0usize;
    for s in block.0..=block.1 {
        // NCBI reference: c++/include/algo/blast/core/blast_lookup.h:121-127
        // ```c
        // return ((index << charsize) | word[wordsize - 1]) & mask;
        // ```
        // SAFETY: `s <= s_last`, the positions NCBI's loop reads.
        let new_char = unsafe { *subject.get_unchecked(s + word_length - 1) } as usize;
        ix = ((ix << charsize) | new_char) & mask;
        // NCBI reference: c++/include/algo/blast/core/blast_lookup.h:55-57 (PV_TEST)
        // SAFETY: `ix` is masked to the backbone space the presence vector covers.
        let bit = (unsafe { *pv.add(ix >> PV_ARRAY_BTS) } >> (ix & PV_ARRAY_MASK)) & 1;
        // SAFETY: `n <= s - block.0 < X_BLOCK`.
        unsafe {
            cand.get_unchecked_mut(n)
                .write(((s as u64) << 32) | ix as u64);
        }
        n += bit as usize;
    }
    *index = ix;
    n
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::search_seed::encode_tblastn_lookup_query;
    use crate::algorithm::tblastx::blast_aascan::s_blast_aa_scan_subject_one_range;
    use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
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
            let last = len as i32 - wl;
            let ranges = [
                (0, last.max(0)),
                (1.min(last.max(0)), last.max(0)),
                (len as i32 / 3, last.max(len as i32 / 3)),
            ];
            for &range in &ranges {
                for array_size in [
                    longest,
                    longest + 1,
                    longest + 2,
                    3 * longest,
                    longest * 1024,
                ] {
                    let want = calls(&lookup, &subject, range, array_size, false);
                    let got = calls(&lookup, &subject, range, array_size, true);
                    assert_eq!(
                        got, want,
                        "case {case} range {range:?} array_size {array_size}"
                    );
                    compared += want.iter().map(|(pairs, _)| pairs.len()).sum::<usize>();
                    stops += want.len().saturating_sub(1);
                }
            }
        }
        // The cases above must reach both the capacity stop and many hits.
        assert!(compared > 100_000, "{compared} pairs compared");
        assert!(stops > 100, "{stops} capacity stops");
    }
}
