//! Subject scanning for TBLASTX
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:48-131
//!
//! This module contains the subject scanning functions that iterate through
//! subject sequences and find k-mer hits against the lookup table.

use super::lookup::{BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK};
use super::scan::OffsetPair;

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

    // NCBI subject masking support (masksubj.inl). This walks subject->seq_ranges.
    // Reference: ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
    let word_length_i32 = lookup.word_length as i32;
    let lut_word_length_i32 = lookup.word_length as i32;
    while determine_scanning_offsets(seq_ranges, word_length_i32, lut_word_length_i32, s_range) {
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        let word_length = lookup.word_length as usize;
        let charsize = lookup.charsize as usize;
        let mask = lookup.mask as usize;

        // Cache pointers for hot loop - matches NCBI pointer-based iteration
        let pv = lookup.pv.as_ptr();
        let backbone = lookup.backbone.as_ptr();
        let overflow = lookup.overflow.as_ptr();

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
