//! Subject scanning for TBLASTX
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:48-131
//!
//! This module contains the subject scanning functions that iterate through
//! subject sequences and find k-mer hits against the lookup table.

use super::lookup::{BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK};
use super::scan::OffsetPair;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
// ```c
// typedef union BlastOffsetPair {
//     struct {
//         Uint4 q_off;
//         Uint4 s_off;
//     } qs_offsets;
//     ...
// } BlastOffsetPair;
// ```
pub(crate) use super::scan::OffsetPair as BlastOffsetPair;

/// NCBI s_DetermineScanningOffsets - advance scan range across subject seq_ranges.
/// Reference: ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
#[inline(always)]
fn determine_scanning_offsets(
    seq_ranges: &[(i32, i32)],
    word_length: i32,
    lut_word_length: i32,
    range: &mut [i32; 3],
) -> bool {
    while range[1] > range[2] {
        range[0] += 1;
        if range[0] >= seq_ranges.len() as i32 {
            return false;
        }
        let (left, right) = seq_ranges[range[0] as usize];
        range[1] = left + word_length - lut_word_length;
        range[2] = right - lut_word_length;
    }
    true
}

/// s_BlastAaScanSubject - scan subject sequence for k-mer hits
///
/// Reference: blast_aascan.c:48-131 s_BlastAaScanSubject
///
/// This is the core scanning function that iterates through the subject
/// sequence and finds all positions where k-mers match entries in the
/// lookup table (via presence vector test).
///
/// NCBI comment from blast_aascan.c:
/// ```c
/// /** Scan a subject sequence for word hits.
///  * @param lookup Lookup table structure [in]
///  * @param subject Subject sequence [in]
///  * @param offset_pairs Array for storing query/subject offset pairs [out]
///  * @param array_size Number of elements in offset_pairs that can be filled [in]
///  * @param s_range Subject range to search (modified in place) [in/out]
///  * @return Number of hits found
///  */
/// ```
///
/// Optimizations applied:
/// - Unsafe array access to eliminate bounds checking in hot loop
/// - Pointer-style iteration matching NCBI's C implementation
/// - Inline presence vector test
#[inline]
pub fn s_blast_aa_scan_subject(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3], // [C] Int4 *s_range
) -> i32 {
    let mut totalhits: i32 = 0;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:64-72
    // ```c
    // pv = lookup->pv;
    // bbc = (AaLookupBackboneCell *) lookup->thick_backbone;
    // ovfl = (Int4 *) lookup->overflow;
    // word_length = lookup->word_length;
    // ```
    let word_length = lookup.word_length as usize;
    let word_length_i32 = lookup.word_length as i32;
    let lut_word_length_i32 = lookup.word_length as i32;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let pv = lookup.pv.as_ptr();
    let backbone = lookup.backbone.as_ptr();
    let overflow = lookup.overflow.as_ptr();

    // NCBI subject masking support (masksubj.inl). This walks subject->seq_ranges.
    // Reference: ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
    while determine_scanning_offsets(seq_ranges, word_length_i32, lut_word_length_i32, s_range) {
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        // [C] index = ComputeTableIndex(word_length - 1, lookup->charsize, s_first);
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:79-81
        // Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:96-107
        let mut index: usize = 0;
        for i in 0..(word_length - 1) {
            let ch = unsafe { *subject.get_unchecked(s_first + i) } as usize;
            index = (index << charsize) | ch;
        }

        // [C] for (s = s_first; s <= s_last; s++)
        let mut s = s_first;
        while s <= s_last {
            // [C] index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:85-87
            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:121-127
            let new_char = unsafe { *subject.get_unchecked(s + word_length - 1) } as usize;
            index = ((index << charsize) | new_char) & mask;

            // [C] if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:89-93
            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:55-57
            let pv_word = unsafe { *pv.add(index >> PV_ARRAY_BTS) };
            if (pv_word & (1u32 << (index & PV_ARRAY_MASK))) != 0 {
                // SAFETY: index is within backbone bounds (masked)
                let cell = unsafe { &*backbone.add(index) };
                let numhits = cell.num_used;

                // [C] if (numhits <= (array_size - totalhits)) { ... }
                if numhits <= array_size - totalhits {
                    // NCBI BlastOffsetPair uses Uint4 offsets.
                    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                    let s_off = s as u32;
                    let dest_base = totalhits as usize;

                    if numhits as usize <= AA_HITS_PER_CELL {
                        // NCBI: copy hits from backbone entries.
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115
                        unsafe {
                            let dest = offset_pairs.as_mut_ptr().add(dest_base);
                            for i in 0..numhits as usize {
                                (*dest.add(i)) = OffsetPair {
                                    q_off: cell.entries[i] as u32,
                                    s_off,
                                };
                            }
                        }
                    } else {
                        // NCBI: copy hits from overflow list.
                        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115
                        let cursor = cell.entries[0] as usize;
                        unsafe {
                            let src = overflow.add(cursor);
                            let dest = offset_pairs.as_mut_ptr().add(dest_base);
                            for i in 0..numhits as usize {
                                (*dest.add(i)) = OffsetPair {
                                    q_off: *src.add(i) as u32,
                                    s_off,
                                };
                            }
                        }
                    }
                    totalhits += numhits;
                } else {
                    // Not enough space in the destination array; return early
                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:118-123
                    s_range[1] = s as i32;
                    return totalhits;
                }
            }

            s += 1;
        }

        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:126-127
        s_range[1] = s as i32;
    }

    totalhits
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:495-504
// ```c
// scan_range[0] = 0;
// scan_range[1] = subject->seq_ranges[0].left;
// scan_range[2] = subject->seq_ranges[0].right - wordsize;
// ...
// hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:75-127
// ```c
// while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) {
//     s_first=subject->sequence + s_range[1];
//     s_last=subject->sequence + s_range[2];
//     ...
//     for (s = s_first; s <= s_last; s++) {
//         index = ComputeTableIndexIncremental(..., s, index);
//         if (PV_TEST(pv, index, PV_ARRAY_BTS)) { ... }
//     }
//     s_range[1] = (Int4)(s - subject->sequence);
// }
// ```
#[inline]
pub fn s_blast_aa_scan_subject_one_range(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    if s_range[1] > s_range[2] {
        return 0;
    }

    let mut totalhits: i32 = 0;

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:67-73
    // ```c
    // lookup = (BlastAaLookupTable *) lookup_wrap->lut;
    // pv = lookup->pv;
    // bbc = (AaLookupBackboneCell *) lookup->thick_backbone;
    // ovfl = (Int4 *) lookup->overflow;
    // word_length = lookup->word_length;
    // ```
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let pv = lookup.pv.as_ptr();
    let backbone = lookup.backbone.as_ptr();
    let overflow = lookup.overflow.as_ptr();

    let s_first = s_range[1] as usize;
    let s_last = s_range[2] as usize;

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:96-107
    // ```c
    // for(i = 0; i < wordsize; i++) {
    //     index = (index << charsize) | word[i];
    // }
    // ```
    let mut index: usize = 0;
    for i in 0..(word_length - 1) {
        // SAFETY: NCBI scan ranges are constructed for a sentinel-backed
        // subject buffer; the prime step reads the first `word_length - 1`
        // residues of the current scan word just like ComputeTableIndex.
        let ch = unsafe { *subject.get_unchecked(s_first + i) } as usize;
        index = (index << charsize) | ch;
    }

    let mut s = s_first;
    while s <= s_last {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:121-127
        // ```c
        // return ((index << charsize) | word[wordsize - 1]) & mask;
        // ```
        // SAFETY: The caller-provided scan range follows NCBI's subject
        // sequence contract: positions through `s + word_length - 1` are valid
        // in the sentinel-backed scan buffer.
        let new_char = unsafe { *subject.get_unchecked(s + word_length - 1) } as usize;
        index = ((index << charsize) | new_char) & mask;

        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:89-97
        // ```c
        // if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
        //     numhits = bbc[index].num_used;
        //     if (numhits <= (array_size - totalhits)) {
        // ```
        // SAFETY: `index` is masked to the lookup backbone space, and the PV
        // allocation is sized from that same space during lookup finalization.
        let pv_word = unsafe { *pv.add(index >> PV_ARRAY_BTS) };
        if (pv_word & (1u32 << (index & PV_ARRAY_MASK))) != 0 {
            // SAFETY: The same masked index addresses the finalized NCBI-style
            // thick backbone cell.
            let cell = unsafe { &*backbone.add(index) };
            let numhits = cell.num_used;

            if numhits <= array_size - totalhits {
                let s_off = s as u32;
                let dest_base = totalhits as usize;

                if numhits as usize <= AA_HITS_PER_CELL {
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-114
                    // ```c
                    // src = bbc[index].payload.entries;
                    // for (i = 0; i < numhits; i++) {
                    //     offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
                    //     offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
                    // }
                    // ```
                    unsafe {
                        // SAFETY: The NCBI space check above guarantees
                        // `numhits <= array_size - totalhits`, so this write is
                        // within the caller-provided offset-pair buffer.
                        let dest = offset_pairs.as_mut_ptr().add(dest_base);
                        for i in 0..numhits as usize {
                            (*dest.add(i)) = OffsetPair {
                                q_off: cell.entries[i] as u32,
                                s_off,
                            };
                        }
                    }
                } else {
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-114
                    // ```c
                    // src = &(ovfl[bbc[index].payload.overflow_cursor]);
                    // for (i = 0; i < numhits; i++) {
                    //     offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
                    //     offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
                    // }
                    // ```
                    let cursor = cell.entries[0] as usize;
                    unsafe {
                        // SAFETY: For overflow cells, NCBI lookup construction
                        // stores a cursor to `numhits` contiguous query offsets;
                        // the destination bound is guarded by the same
                        // `array_size - totalhits` check as the inline case.
                        let src = overflow.add(cursor);
                        let dest = offset_pairs.as_mut_ptr().add(dest_base);
                        for i in 0..numhits as usize {
                            (*dest.add(i)) = OffsetPair {
                                q_off: *src.add(i) as u32,
                                s_off,
                            };
                        }
                    }
                }
                totalhits += numhits;
            } else {
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:118-123
                // ```c
                // s_range[1] = (Int4)(s - subject->sequence);
                // return totalhits;
                // ```
                s_range[1] = s as i32;
                return totalhits;
            }
        }

        s += 1;
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:125-127
    // ```c
    // } /* end for */
    // s_range[1] = (Int4)(s - subject->sequence);
    // ```
    s_range[1] = s as i32;
    totalhits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
    use crate::algorithm::tblastx::translation::QueryFrame;
    use crate::config::ScoringMatrix;
    use crate::stats::lookup_protein_params_ungapped;
    use crate::utils::matrix::ncbistdaa;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:87-129
    // ```c
    // for (offset = from; offset <= to; offset++, seq++) {
    //     if (seq >= word_target) {
    //         BlastLookupAddWordHit(backbone, lut_word_length, charsize,
    //                               seq - lut_word_length, offset - lut_word_length);
    //     }
    // }
    // ```
    fn make_poly_a_query_frame() -> QueryFrame {
        QueryFrame {
            frame: 1,
            aa_seq: vec![
                0,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                ncbistdaa::A,
                0,
            ],
            aa_seq_nomask: None,
            aa_len: 6,
            orig_len: 6,
            seg_masks: Vec::new(),
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-456
    // ```c
    // BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
    //                                   lookup->charsize, lookup->word_length,
    //                                   query, location);
    // ```
    fn make_poly_a_lookup() -> BlastAaLookupTable {
        let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let queries = vec![vec![make_poly_a_query_frame()]];
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2746-2748
        // ```c
        // Boolean check_ideal =
        //    (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
        //     program == eBlastTypeRpsTblastn);
        // ```
        let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, true);
        lookup
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/unit_tests/api/aascan_unit_test.cpp:714-752
    // ```c
    // qsort(offset_pairs, hits, sizeof(BlastOffsetPair), compare_offsets);
    // if (s_off)
    //     BOOST_REQUIRE(offset_pairs[0].qs_offsets.s_off > s_off);
    // if (offset_pairs[i].qs_offsets.s_off == offset_pairs[i-1].qs_offsets.s_off)
    //     BOOST_REQUIRE(offset_pairs[i].qs_offsets.q_off > offset_pairs[i-1].qs_offsets.q_off);
    // ```
    #[test]
    fn test_s_blast_aa_scan_subject_small_buffer_resumes_without_splitting_subject_offset() {
        let lookup = make_poly_a_lookup();
        let subject = vec![ncbistdaa::A; 6];
        let seq_ranges = [(0, subject.len() as i32)];
        let mut offset_pairs = vec![OffsetPair::default(); 5];
        let array_size = i32::try_from(offset_pairs.len())
            .expect("NCBI BLAST offset pair array size must fit in Int4");
        let mut scan_range = [0, 0, subject.len() as i32 - lookup.word_length];

        let mut observed_subject_offsets = Vec::new();
        while scan_range[1] <= scan_range[2] {
            let hits = s_blast_aa_scan_subject(
                &lookup,
                &subject,
                &seq_ranges,
                &mut offset_pairs,
                array_size,
                &mut scan_range,
            );
            assert_eq!(hits, 4);

            let batch = &offset_pairs[..hits as usize];
            let expected_subject_offset = observed_subject_offsets.len() as u32;
            assert!(batch
                .iter()
                .all(|pair| pair.s_off == expected_subject_offset));
            assert_eq!(
                batch.iter().map(|pair| pair.q_off).collect::<Vec<_>>(),
                vec![0, 1, 2, 3]
            );
            observed_subject_offsets.push(expected_subject_offset);
        }

        assert_eq!(observed_subject_offsets, vec![0, 1, 2, 3]);
    }
}
