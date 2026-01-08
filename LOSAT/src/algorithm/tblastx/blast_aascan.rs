//! Subject scanning for TBLASTX
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:48-131
//!
//! This module contains the subject scanning functions that iterate through
//! subject sequences and find k-mer hits against the lookup table.

use super::lookup::{
    decode_kmer, encode_kmer, get_charsize, get_mask, BlastAaLookupTable, AA_HITS_PER_CELL,
    LOOKUP_ALPHABET_SIZE,
};
use super::scan::OffsetPair;
use crate::utils::matrix::blosum62_score;

#[cfg(target_arch = "x86_64")]
use super::scan::{
    compute_3mer_indices_16_avx2, compute_3mer_indices_4_scalar, compute_3mer_indices_8_sse2,
    copy_offset_pairs_overflow_avx2, copy_offset_pairs_overflow_sse2, pv_test_mask16_avx2,
    pv_test_mask4_avx2, pv_test_mask8_avx2,
};

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
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3], // [C] Int4 *s_range
) -> i32 {
    let mut totalhits: i32 = 0;

    // Precompute CPU feature flags once per scan call.
    // Use LOSAT_NO_SIMD=1 to force scalar-only processing for debugging
    #[cfg(target_arch = "x86_64")]
    let has_avx2 = is_x86_feature_detected!("avx2") && std::env::var("LOSAT_NO_SIMD").is_err();
    #[cfg(not(target_arch = "x86_64"))]
    let has_avx2 = false;

    // NCBI subject masking support (masksubj.inl). In NCBI this walks
    // `subject->seq_ranges[]`; in LOSAT each translated frame is one contiguous
    // range already baked into `s_range[1..=2]`.
    while s_range[1] <= s_range[2] {
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        // Fast fail: nothing to scan.
        if s_first > s_last || subject.len() < 5 {
            return totalhits;
        }

        let charsize = lookup.charsize as usize;
        let mask = lookup.mask as usize;

        // Cache pointers for hot loop - matches NCBI pointer-based iteration
        let pv = lookup.pv.as_ptr();
        let backbone = lookup.backbone.as_ptr();
        let overflow = lookup.overflow.as_ptr();

        // [C] index = ComputeTableIndex(word_length - 1, lookup->charsize, s_first);
        // Prime the rolling index with first (word_length - 1) residues
        // SAFETY: s_first + 1 < subject.len() is guaranteed by s_first <= s_last and subject.len() >= 5
        let mut index: usize = unsafe {
            (*subject.get_unchecked(s_first) as usize) << charsize
                | (*subject.get_unchecked(s_first + 1) as usize)
        };

        // [C] for (s = s_first; s <= s_last; s++)
        let mut s = s_first;

        // AVX2: batch PV_TEST for consecutive positions using direct 3-mer computation
        // This replaces rolling index with direct encoding, allowing better SIMD utilization
        #[cfg(target_arch = "x86_64")]
        if has_avx2 {
            unsafe {
                while s + 15 <= s_last {
                    // Compute 16 consecutive 3-mer indices directly (no rolling dependency)
                    let idxs = compute_3mer_indices_16_avx2(subject.as_ptr(), s);

                    let hit_mask = pv_test_mask16_avx2(pv, &idxs);

                    // Process hits strictly in order (NCBI-compatible)
                    for lane in 0..16usize {
                        if (hit_mask & (1u16 << lane)) == 0 {
                            continue;
                        }

                        let idx = idxs[lane];
                        let cell = &*backbone.add(idx);
                        let numhits = cell.num_used;

                        if numhits <= array_size - totalhits {
                            let s_off = (s + lane) as i32;
                            let dest_base = totalhits as usize;

                            if numhits as usize <= AA_HITS_PER_CELL {
                                let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                match numhits {
                                    1 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                    }
                                    2 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                        (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                    }
                                    3 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                        (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                        (*dest.add(2)) = OffsetPair { q_off: cell.entries[2], s_off };
                                    }
                                    _ => {}
                                }
                            } else {
                                let cursor = cell.entries[0] as usize;
                                let src = overflow.add(cursor);
                                let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                copy_offset_pairs_overflow_avx2(src, dest, numhits as usize, s_off);
                            }

                            totalhits += numhits;
                        } else {
                            // Not enough space in the destination array; return early
                            s_range[1] = (s + lane) as i32;
                            return totalhits;
                        }
                    }

                    // Update rolling index for next iteration (needed for scalar fallback)
                    // We compute the last index from the batch to maintain continuity
                    if s + 16 <= s_last {
                        let last_idx = idxs[15];
                        index = last_idx;
                    }
                    s += 16;
                }

                while s + 7 <= s_last {
                    // Compute 8 consecutive 3-mer indices directly
                    let idxs = compute_3mer_indices_8_sse2(subject.as_ptr(), s);

                    let hit_mask = pv_test_mask8_avx2(pv, &idxs);

                    // Process hits strictly in order (NCBI-compatible)
                    for lane in 0..8usize {
                        if (hit_mask & (1u8 << lane)) == 0 {
                            continue;
                        }

                        let idx = idxs[lane];
                        let cell = &*backbone.add(idx);
                        let numhits = cell.num_used;

                        if numhits <= array_size - totalhits {
                            let s_off = (s + lane) as i32;
                            let dest_base = totalhits as usize;

                            if numhits as usize <= AA_HITS_PER_CELL {
                                let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                match numhits {
                                    1 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                    }
                                    2 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                        (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                    }
                                    3 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                        (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                        (*dest.add(2)) = OffsetPair { q_off: cell.entries[2], s_off };
                                    }
                                    _ => {}
                                }
                            } else {
                                let cursor = cell.entries[0] as usize;
                                let src = overflow.add(cursor);
                                let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                copy_offset_pairs_overflow_avx2(src, dest, numhits as usize, s_off);
                            }

                            totalhits += numhits;
                        } else {
                            // Not enough space in the destination array; return early
                            s_range[1] = (s + lane) as i32;
                            return totalhits;
                        }
                    }

                    // Update rolling index for next iteration
                    if s + 8 <= s_last {
                        let last_idx = idxs[7];
                        index = last_idx;
                    }
                    s += 8;
                }

                while s + 3 <= s_last {
                    // Compute 4 consecutive 3-mer indices directly
                    let idxs = compute_3mer_indices_4_scalar(subject.as_ptr(), s);

                    let hit_mask = pv_test_mask4_avx2(pv, &idxs);

                    // Process hits strictly in order (NCBI-compatible)
                    for lane in 0..4usize {
                        if (hit_mask & (1u8 << lane)) == 0 {
                            continue;
                        }

                        let idx = idxs[lane];
                        let cell = &*backbone.add(idx);
                        let numhits = cell.num_used;

                        if numhits <= array_size - totalhits {
                            let s_off = (s + lane) as i32;
                            let dest_base = totalhits as usize;

                            if numhits as usize <= AA_HITS_PER_CELL {
                                let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                match numhits {
                                    1 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                    }
                                    2 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                        (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                    }
                                    3 => {
                                        (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                        (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                        (*dest.add(2)) = OffsetPair { q_off: cell.entries[2], s_off };
                                    }
                                    _ => {}
                                }
                            } else {
                                let cursor = cell.entries[0] as usize;
                                let src = overflow.add(cursor);
                                let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                copy_offset_pairs_overflow_avx2(src, dest, numhits as usize, s_off);
                            }

                            totalhits += numhits;
                        } else {
                            // Not enough space in the destination array; return early
                            s_range[1] = (s + lane) as i32;
                            return totalhits;
                        }
                    }

                    // Update rolling index for next iteration
                    if s + 4 <= s_last {
                        let last_idx = idxs[3];
                        index = last_idx;
                    }
                    s += 4;
                }
            }
        }

        while s <= s_last {
            // [C] index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
            // Rolling index computation: shift in the new character
            // SAFETY: s + 2 <= s_last + 2 < subject.len() is guaranteed by loop bounds
            let new_char = unsafe { *subject.get_unchecked(s + 2) } as usize;
            index = ((index << charsize) | new_char) & mask;

            // [C] if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
            // Inline presence vector test for performance
            // SAFETY: index is masked, pv array is sized to cover all possible indices
            let pv_word = unsafe { *pv.add(index >> 6) };
            if (pv_word & (1u64 << (index & 63))) != 0 {
                // SAFETY: index is within backbone bounds (masked)
                let cell = unsafe { &*backbone.add(index) };
                let numhits = cell.num_used;

                // [C] if (numhits <= (array_size - totalhits)) { ... }
                if numhits <= array_size - totalhits {
                    let s_off = s as i32;
                    let dest_base = totalhits as usize;

                    if numhits as usize <= AA_HITS_PER_CELL {
                        // Hits in backbone cell - unroll for common cases
                        // SAFETY: dest_base + numhits <= offset_pairs.len() guaranteed by array_size check
                        unsafe {
                            let dest = offset_pairs.as_mut_ptr().add(dest_base);
                            match numhits {
                                1 => {
                                    (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                }
                                2 => {
                                    (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                    (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                }
                                3 => {
                                    (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                    (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                    (*dest.add(2)) = OffsetPair { q_off: cell.entries[2], s_off };
                                }
                                _ => {}
                            }
                        }
                    } else {
                        // Hits in overflow array
                        let cursor = cell.entries[0] as usize;
                        // SAFETY: overflow bounds checked during lookup table construction
                        unsafe {
                            let src = overflow.add(cursor);
                            let dest = offset_pairs.as_mut_ptr().add(dest_base);
                            #[cfg(target_arch = "x86_64")]
                            {
                                if has_avx2 {
                                    copy_offset_pairs_overflow_avx2(src, dest, numhits as usize, s_off);
                                } else {
                                    copy_offset_pairs_overflow_sse2(src, dest, numhits as usize, s_off);
                                }
                            }
                            #[cfg(not(target_arch = "x86_64"))]
                            {
                                for i in 0..numhits as usize {
                                    (*dest.add(i)) = OffsetPair { q_off: *src.add(i), s_off };
                                }
                            }
                        }
                    }
                    totalhits += numhits;
                } else {
                    // Not enough space in the destination array; return early
                    s_range[1] = s as i32;
                    return totalhits;
                }
            }

            s += 1;
        }

        s_range[1] = s as i32;
    }

    totalhits
}

/// Lazy neighbor scan - compute neighbors dynamically during scan.
///
/// For each subject k-mer, we compute all neighboring query k-mers that would
/// produce a score >= threshold, then check if those neighbors exist in the
/// lookup table (which contains only exact query k-mers).
///
/// This trades scan-time computation for drastically reduced lookup table size.
#[inline(never)]
#[allow(dead_code)]
pub fn s_blast_aa_scan_subject_lazy(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    let mut totalhits: i32 = 0;
    let threshold = lookup.threshold;
    let row_max = &lookup.row_max;

    while s_range[1] <= s_range[2] {
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        if s_first > s_last || subject.len() < 5 {
            return totalhits;
        }

        let charsize = get_charsize();
        let mask = get_mask();
        let alphabet_size = LOOKUP_ALPHABET_SIZE;

        let pv = lookup.pv.as_ptr();
        let backbone = lookup.backbone.as_ptr();
        let overflow = lookup.overflow.as_ptr();

        // Rolling index for subject k-mer
        let mut index: usize = unsafe {
            (*subject.get_unchecked(s_first) as usize) << charsize
                | (*subject.get_unchecked(s_first + 1) as usize)
        };

        let mut s = s_first;
        while s <= s_last {
            let new_char = unsafe { *subject.get_unchecked(s + 2) } as usize;
            index = ((index << charsize) | new_char) & mask;

            // Decode subject k-mer
            let (s0, s1, s2) = decode_kmer(index);

            // Skip invalid residues
            if s0 >= alphabet_size || s1 >= alphabet_size || s2 >= alphabet_size {
                s += 1;
                continue;
            }

            // For lazy mode, we need to find all query k-mers Q such that
            // score(Q, S) >= threshold, where S is the subject k-mer.
            //
            // This is the REVERSE of the standard neighbor generation:
            // - Standard: for each query k-mer Q, find all S such that score(Q, S) >= threshold
            // - Lazy: for each subject k-mer S, find all Q such that score(Q, S) >= threshold
            //
            // We enumerate all possible Q and check if:
            // 1. score(Q, S) >= threshold
            // 2. Q exists in the lookup table (pv_test)
            //
            // Pruning with row_max:
            // score(Q, S) = blosum62(q0, s0) + blosum62(q1, s1) + blosum62(q2, s2)
            // max_score_for_s = row_max[s0] + row_max[s1] + row_max[s2]
            // If max_score_for_s < threshold, no Q can match.

            let max_possible = row_max[s0] + row_max[s1] + row_max[s2];
            if max_possible < threshold {
                s += 1;
                continue;
            }

            // Enumerate potential query k-mers with pruning
            let rm12 = row_max[s1] + row_max[s2];
            let rm2 = row_max[s2];

            for q0 in 0..alphabet_size {
                let sc0 = blosum62_score(q0 as u8, s0 as u8);
                if sc0 + rm12 < threshold {
                    continue;
                }
                for q1 in 0..alphabet_size {
                    let sc1 = sc0 + blosum62_score(q1 as u8, s1 as u8);
                    if sc1 + rm2 < threshold {
                        continue;
                    }
                    for q2 in 0..alphabet_size {
                        let total_score = sc1 + blosum62_score(q2 as u8, s2 as u8);
                        if total_score < threshold {
                            continue;
                        }

                        // This query k-mer Q = (q0, q1, q2) is a neighbor of subject k-mer S
                        let q_idx = encode_kmer(q0, q1, q2);

                        // Check if Q exists in lookup table
                        let pv_word = unsafe { *pv.add(q_idx >> 6) };
                        if (pv_word & (1u64 << (q_idx & 63))) == 0 {
                            continue;
                        }

                        // Q exists - get its hits
                        let cell = unsafe { &*backbone.add(q_idx) };
                        let numhits = cell.num_used;

                        if numhits <= 0 {
                            continue;
                        }

                        if numhits <= array_size - totalhits {
                            let s_off = s as i32;
                            let dest_base = totalhits as usize;

                            if numhits as usize <= AA_HITS_PER_CELL {
                                unsafe {
                                    let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                    for i in 0..numhits as usize {
                                        (*dest.add(i)) = OffsetPair {
                                            q_off: cell.entries[i],
                                            s_off,
                                        };
                                    }
                                }
                            } else {
                                let cursor = cell.entries[0] as usize;
                                unsafe {
                                    let src = overflow.add(cursor);
                                    let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                    for i in 0..numhits as usize {
                                        (*dest.add(i)) = OffsetPair {
                                            q_off: *src.add(i),
                                            s_off,
                                        };
                                    }
                                }
                            }
                            totalhits += numhits;
                        } else {
                            s_range[1] = s as i32;
                            return totalhits;
                        }
                    }
                }
            }

            s += 1;
        }

        s_range[1] = s as i32;
    }

    totalhits
}
