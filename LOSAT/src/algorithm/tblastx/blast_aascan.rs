//! Subject scanning for TBLASTX
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:48-131
//!
//! This module contains the subject scanning functions that iterate through
//! subject sequences and find k-mer hits against the lookup table.

use super::lookup::{
    BlastAaLookupTable, XCell, XHit, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use super::scan::OffsetPair;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
pub(crate) use super::scan::OffsetPair as BlastOffsetPair;

/// NCBI s_DetermineScanningOffsets - advance scan range across subject seq_ranges.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58
#[inline(always)]
pub(super) fn determine_scanning_offsets(
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58 (the scan moves to the next
// range of the list, never back) and c++/src/algo/blast/core/blast_aascan.c:75-123 (the words of a
// range left to right). The first range starts the scan at its left even when it holds no word
// (c++/src/algo/blast/core/aa_ungapped.c:497-500). The ranges of a database sequence with
// `-db_soft_mask` are the complement of its masks in the order stored, unsorted
// (c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2490-2506,
// c++/src/algo/blast/api/seqsrc_seqdb.cpp:398-405, c++/src/algo/blast/core/blast_util.c:211-215).
/// Where a word finder's scan of a chunk's `seq_ranges` has come, for the order check of its hits
/// (debug builds): the range it scans and the subject offset of its last hit. The scan takes the
/// ranges in list order and the words of each from left to right, so the subject offsets of the
/// hits rise within a range and fall where a range starts before the words of one scanned before
/// it: masks stored out of order (`[600, 700)` before `[300, 390)`) give the ranges `[0, 600)`,
/// `[700, 300)` (no word) and `[390, ...)`, and the scan reads `[390, 600)` again.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct ScanOrder {
    range: usize,
    last: Option<i32>,
}

impl ScanOrder {
    /// Whether the scan of `seq_ranges` (words of `word_length` letters) can give its next hit at
    /// subject offset `s_off`: in the range it scans, at or after its last hit, or in a later
    /// range that has a word there. The scan is then taken to the first such range. A range's
    /// words start from its left to its right less a word; the first range's, at its left at
    /// least.
    pub(crate) fn takes(
        &mut self,
        seq_ranges: &[(i32, i32)],
        word_length: i32,
        s_off: i32,
    ) -> bool {
        for (range, &(left, right)) in seq_ranges.iter().enumerate().skip(self.range) {
            let last_start = if range == 0 {
                (right - word_length).max(left)
            } else {
                right - word_length
            };
            let after_last = range > self.range || self.last.is_none_or(|last| s_off >= last);
            if after_last && (left..=last_start).contains(&s_off) {
                self.range = range;
                self.last = Some(s_off);
                return true;
            }
        }
        false
    }
}

/// s_BlastAaScanSubject - scan subject sequence for k-mer hits
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:48-131 s_BlastAaScanSubject
///
/// This is the core scanning function that iterates through the subject
/// sequence and finds all positions where k-mers match entries in the
/// lookup table (via presence vector test).
///
/// NCBI comment from c++/src/algo/blast/core/blast_aascan.c:
///
/// Optimizations applied:
/// - Unsafe array access to eliminate bounds checking in hot loop
/// - Pointer-style iteration matching NCBI's C implementation
/// - Inline presence vector test
#[inline]
pub(crate) fn s_blast_aa_scan_subject(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3], // [C] Int4 *s_range
) -> i32 {
    let mut totalhits: i32 = 0;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:64-73
    let word_length = lookup.word_length as usize;
    let word_length_i32 = lookup.word_length;
    let lut_word_length_i32 = lookup.word_length;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let pv = lookup.pv.as_ptr();
    let backbone = lookup.backbone.as_ptr();
    let overflow = lookup.overflow.as_ptr();

    // NCBI subject masking support (c++/src/algo/blast/core/masksubj.inl). This walks subject->seq_ranges.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58
    while determine_scanning_offsets(seq_ranges, word_length_i32, lut_word_length_i32, s_range) {
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        // [C] index = ComputeTableIndex(word_length - 1, lookup->charsize, s_first);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:79-81
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:96-107
        let mut index: usize = 0;
        for i in 0..(word_length - 1) {
            let ch = unsafe { *subject.get_unchecked(s_first + i) } as usize;
            index = (index << charsize) | ch;
        }

        // [C] for (s = s_first; s <= s_last; s++)
        let mut s = s_first;
        while s <= s_last {
            // [C] index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:85-87
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:121-127
            let new_char = unsafe { *subject.get_unchecked(s + word_length - 1) } as usize;
            index = ((index << charsize) | new_char) & mask;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:89-93
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:55-57
            let pv_word = unsafe { *pv.add(index >> PV_ARRAY_BTS) };
            if (pv_word & (1u32 << (index & PV_ARRAY_MASK))) != 0 {
                // SAFETY: `index < backbone_size`: the crate's callers pass NCBIstdaa letters (below
                // 28), so the word's index stays in the table's cells (the mask alone allows 32,768).
                let cell = unsafe { &*backbone.add(index) };
                let numhits = cell.num_used;

                // [C] if (numhits <= (array_size - totalhits)) { ... }
                if numhits <= array_size - totalhits {
                    // NCBI BlastOffsetPair uses Uint4 offsets.
                    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
                    let s_off = s as u32;
                    let dest_base = totalhits as usize;

                    if numhits as usize <= AA_HITS_PER_CELL {
                        // NCBI: copy hits from backbone entries.
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-115
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
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-115
                        let cursor = cell.cursor();
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
                    // Output buffer is full: record where the scan stopped and return
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:118-123
                    s_range[1] = s as i32;
                    return totalhits;
                }
            }

            s += 1;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:126-127
        s_range[1] = s as i32;
    }

    totalhits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:60-127 (cited at
// `s_blast_aa_scan_subject` and `s_blast_aa_scan_subject_one_range`)
// No NCBI counterpart for the form: the scan of `s_blast_aa_scan_subject` over a lookup table of
// either cell width, with checked reads. BLASTX's table of 16-bit hits (`BackboneCell16`, built
// where the fused loop scans it) is read by the reference scan only in the diagnostic stages and
// for a `-subject` with lower-case letters; the table of `Int4` hits keeps the functions above.
/// The pairs of `s_blast_aa_scan_subject` for a table of cells `C`: the same scan ranges, the same
/// pairs in the same order, the same early return when `offset_pairs` is full.
pub(crate) fn x_blast_aa_scan_subject_cells<C: XCell>(
    lookup: &BlastAaLookupTable<C>,
    subject: &[u8],
    seq_ranges: &[(i32, i32)],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    let mut totalhits: i32 = 0;
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    while determine_scanning_offsets(seq_ranges, lookup.word_length, lookup.word_length, s_range) {
        let (s_first, s_last) = (s_range[1] as usize, s_range[2] as usize);
        let mut index: usize = 0;
        for &ch in &subject[s_first..s_first + word_length - 1] {
            index = (index << charsize) | ch as usize;
        }
        let mut s = s_first;
        while s <= s_last {
            index = ((index << charsize) | subject[s + word_length - 1] as usize) & mask;
            if lookup.pv[index >> PV_ARRAY_BTS] & (1u32 << (index & PV_ARRAY_MASK)) != 0 {
                let cell = &lookup.backbone[index];
                let numhits = cell.count() as i32;
                if numhits > array_size - totalhits {
                    s_range[1] = s as i32;
                    return totalhits;
                }
                for (i, hit) in cell.hits(&lookup.overflow).iter().enumerate() {
                    offset_pairs[totalhits as usize + i] = OffsetPair {
                        q_off: hit.q() as u32,
                        s_off: s as u32,
                    };
                }
                totalhits += numhits;
            }
            s += 1;
        }
        s_range[1] = s as i32;
    }
    totalhits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-504
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127
#[inline]
pub(crate) fn s_blast_aa_scan_subject_one_range(
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:67-73
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let pv = lookup.pv.as_ptr();
    let backbone = lookup.backbone.as_ptr();
    let overflow = lookup.overflow.as_ptr();

    let s_first = s_range[1] as usize;
    let s_last = s_range[2] as usize;

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:96-107
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
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:121-127
        // SAFETY: The caller-provided scan range follows NCBI's subject
        // sequence contract: positions through `s + word_length - 1` are valid
        // in the sentinel-backed scan buffer.
        let new_char = unsafe { *subject.get_unchecked(s + word_length - 1) } as usize;
        index = ((index << charsize) | new_char) & mask;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:89-97
        // Bounds (not yet a SAFETY argument): `index < backbone_size` needs every letter of the
        // word below 28 (the mask alone allows 32,768), and every builder gives `pv`
        // `(backbone_size >> PV_ARRAY_BTS) + 1` words. The crate's callers pass NCBIstdaa letters
        // and words inside `subject` (TBLASTN's `x_scan_unit` does not scan a raised first-range
        // word past it).
        let pv_word = unsafe { *pv.add(index >> PV_ARRAY_BTS) };
        if (pv_word & (1u32 << (index & PV_ARRAY_MASK))) != 0 {
            // Bounds: `index < backbone_size` as for `pv` above addresses one of the table's cells.
            let cell = unsafe { &*backbone.add(index) };
            let numhits = cell.num_used;

            if numhits <= array_size - totalhits {
                let s_off = s as u32;
                let dest_base = totalhits as usize;

                if numhits as usize <= AA_HITS_PER_CELL {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-114
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
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:99-114
                    let cursor = cell.cursor();
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
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:118-123
                s_range[1] = s as i32;
                return totalhits;
            }
        }

        s += 1;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:125-127
    s_range[1] = s as i32;
    totalhits
}

#[cfg(test)]
#[path = "blast_aascan_tests.rs"]
mod tests;
