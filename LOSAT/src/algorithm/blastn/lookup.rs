use bio::io::fasta;
use rustc_hash::FxHashMap;
use crate::core::blast_encoding::{encode_iupac_to_ncbi2na_packed, COMPRESSION_RATIO};
use crate::utils::dust::MaskedInterval;
use super::constants::MAX_DIRECT_LOOKUP_WORD_SIZE;

/// 2-bit encoding for compact storage and hashing
pub fn encode_kmer(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    if start + k > seq.len() {
        return None;
    }
    let mut encoded: u64 = 0;
    for i in 0..k {
        let b = unsafe { *seq.get_unchecked(start + i) };
        let code = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' | b'U' | b'u' => 3,
            _ => return None,
        };
        encoded = (encoded << 2) | code;
    }
    Some(encoded)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:80-103
// ```c
// const char NCBI4NA_TO_IUPACNA[BLASTNA_SIZE] = {
//     '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
//     'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
// };
// const Uint1 IUPACNA_TO_NCBI4NA[128]={ ... };
// ```
const NCBI4NA_TO_IUPACNA: [u8; 16] = [
    b'-', b'A', b'C', b'M', b'G', b'R', b'S', b'V',
    b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:95-103
// ```c
// const Uint1 IUPACNA_TO_NCBI4NA[128]={ ... };
// ```
const IUPACNA_TO_NCBI4NA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0,
    0, 0, 5, 6, 8, 0, 7, 9, 0, 10, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:812-819
// ```c
// Uint1 conversion_table[16] = {
//   0,  8, 4, 12,
//   2, 10, 6, 14,
//   1,  9, 5, 13,
//   3, 11, 7, 15
// };
// ```
const NCBI4NA_REV_COMP: [u8; 16] = [
    0, 8, 4, 12,
    2, 10, 6, 14,
    1, 9, 5, 13,
    3, 11, 7, 15,
];

/// Generate the reverse complement of a DNA sequence in IUPACNA.
/// Uses NCBI4NA bitmask mapping to preserve ambiguous base complements.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| {
            let idx = if b < 128 { IUPACNA_TO_NCBI4NA[b as usize] } else { 0 };
            let rc = NCBI4NA_REV_COMP[idx as usize];
            NCBI4NA_TO_IUPACNA[rc as usize]
        })
        .collect()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1040-1046
// ```c
// /* Also add 1 to all indices, because lookup table indices count
//    from 1. */
// mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
// mb_lt->hashtable[ecode] = index;
// ```
pub type KmerLookup = FxHashMap<u64, Vec<u32>>;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:251-258
// ```c
// Int4* hashtable;   /**< Array of positions              */
// Int4* next_pos;    /**< Extra positions stored here     */
// PV_ARRAY_TYPE *pv_array;/**< Presence vector, used for quick presence
//                            check */
// ```
/// Direct address table for k-mer lookup - packed offsets + hits.
/// Hits store 1-based query offsets (q_off + 1), matching NCBI lookup chains.
/// Used for small word sizes (<=13) where 4^word_size fits in memory.
pub struct DirectKmerLookup {
    offsets: Vec<u32>,
    hits: Vec<u32>,
}

impl DirectKmerLookup {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
    // ```c
    // Int4 q_off = lookup->hashtable[index];
    // while (q_off) {
    //     offset_pairs[i].qs_offsets.q_off   = q_off - 1;
    //     offset_pairs[i++].qs_offsets.s_off = s_off;
    //     q_off = lookup->next_pos[q_off];
    // }
    // ```
    #[inline(always)]
    pub fn get(&self, idx: usize) -> &[u32] {
        if idx + 1 >= self.offsets.len() {
            return &[];
        }
        let start = self.offsets[idx] as usize;
        let end = self.offsets[idx + 1] as usize;
        &self.hits[start..end]
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_util.h:51-55
// ```c
// #define NCBI2NA_UNPACK_BASE(x, N) (((x)>>(2*(N))) & NCBI2NA_MASK)
// ```
#[inline(always)]
fn packed_base_at(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos / COMPRESSION_RATIO];
    let shift = 2 * (3 - (pos % COMPRESSION_RATIO));
    (byte >> shift) & 0x03
}

#[inline]
fn packed_kmer_at(packed: &[u8], start: usize, k: usize) -> u64 {
    let mut kmer = 0u64;
    for i in 0..k {
        let code = packed_base_at(packed, start + i) as u64;
        kmer = (kmer << 2) | code;
    }
    kmer
}

/// Compute database word counts for lookup filtering (limit_lookup).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1122-1177
/// ```c
/// word = (w >> shift) & mask;
/// if (!PV_TEST(pv, word, pv_array_bts)) continue;
/// if ((counts[index] & 0xf) < max_word_count) counts[index]++;
/// ```
pub fn build_db_word_counts(
    queries: &[fasta::Record],
    query_masks: &[Vec<MaskedInterval>],
    subjects: &[fasta::Record],
    lut_word_length: usize,
    max_word_count: u8,
    approx_table_entries: usize,
) -> Vec<u8> {
    if lut_word_length == 0 {
        return Vec::new();
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1250-1254
    // ```c
    // mb_lt->hashsize = 1ULL << (BITS_PER_NUC * mb_lt->lut_word_length);
    // ```
    let hashsize = 1usize << (2 * lut_word_length);
    let mut counts = vec![0u8; hashsize / 2];
    let (pv, pv_array_bts) =
        build_query_pv(queries, query_masks, lut_word_length, approx_table_entries);

    for record in subjects {
        let seq = record.seq();
        if seq.len() < lut_word_length {
            continue;
        }
        let packed = encode_iupac_to_ncbi2na_packed(seq);
        let mask = (1u64 << (2 * lut_word_length)) - 1;

        let mut pos = 0usize;
        let end = seq.len() - lut_word_length;
        let mut kmer = packed_kmer_at(&packed, 0, lut_word_length);

        loop {
            if pv_test_shift(&pv, kmer as usize, pv_array_bts) {
                db_word_count_increment(&mut counts, kmer, max_word_count);
            }

            if pos == end {
                break;
            }

            let next_base = packed_base_at(&packed, pos + lut_word_length);
            kmer = ((kmer << 2) | next_base as u64) & mask;
            pos += 1;
        }
    }

    counts
}

// ============================================================================
// Phase 2: Presence-Vector (PV) for fast k-mer filtering
// ============================================================================
// Following NCBI BLAST's approach from blast_lookup.h:
// - PV_ARRAY_TYPE is u32 (32-bit unsigned integer)
// - PV_ARRAY_BTS is 5 (bits-to-shift from lookup_index to pv_array index)
// - Each bit indicates whether a k-mer has any hits in the lookup table
// - This allows O(1) filtering before accessing the lookup table

/// Presence-Vector array type (matches NCBI BLAST's PV_ARRAY_TYPE)
type PvArrayType = u32;

/// Bits-to-shift from lookup index to PV array index (matches NCBI BLAST's PV_ARRAY_BTS)
const PV_ARRAY_BTS: usize = 5;
/// Bytes per PV array element (matches NCBI BLAST's PV_ARRAY_BYTES)
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:42-43
/// ```c
/// #define PV_ARRAY_BYTES 4
/// #define PV_ARRAY_BTS 5
/// ```
const PV_ARRAY_BYTES: usize = 4;

/// Mask for extracting bit position within a PV array element
const PV_ARRAY_MASK: usize = (1 << PV_ARRAY_BTS) - 1; // 31

/// Test if a k-mer is present in the presence vector
/// Equivalent to NCBI BLAST's PV_TEST macro
#[inline(always)]
fn pv_test(pv: &[PvArrayType], index: usize) -> bool {
    let array_idx = index >> PV_ARRAY_BTS;
    let bit_pos = index & PV_ARRAY_MASK;
    if array_idx < pv.len() {
        (pv[array_idx] & (1u32 << bit_pos)) != 0
    } else {
        false
    }
}

/// Set a bit in the presence vector
/// Equivalent to NCBI BLAST's PV_SET macro
#[inline(always)]
fn pv_set(pv: &mut [PvArrayType], index: usize) {
    let array_idx = index >> PV_ARRAY_BTS;
    let bit_pos = index & PV_ARRAY_MASK;
    if array_idx < pv.len() {
        pv[array_idx] |= 1u32 << bit_pos;
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:51-57
// ```c
// #define PV_SET(lookup, index, shift) \
//     lookup[(index) >> (shift)] |= (PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)
// #define PV_TEST(lookup, index, shift) \
//     ( lookup[(index) >> (shift)] & ((PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)) )
// ```
#[inline(always)]
fn pv_test_shift(pv: &[PvArrayType], index: usize, shift: usize) -> bool {
    let array_idx = index >> shift;
    let bit_pos = index & PV_ARRAY_MASK;
    if array_idx < pv.len() {
        (pv[array_idx] & (1u32 << bit_pos)) != 0
    } else {
        false
    }
}

#[inline(always)]
fn pv_set_shift(pv: &mut [PvArrayType], index: usize, shift: usize) {
    let array_idx = index >> shift;
    let bit_pos = index & PV_ARRAY_MASK;
    if array_idx < pv.len() {
        pv[array_idx] |= 1u32 << bit_pos;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_util.c:71-83
// ```c
// Int4 ilog2(Int8 x)
// {
//     Int4 lg = 0;
//     if (x == 0) return 0;
//     while ((x = x >> 1)) lg++;
//     return lg;
// }
// ```
#[inline]
fn ilog2(mut x: usize) -> usize {
    let mut lg = 0usize;
    if x == 0 {
        return 0;
    }
    while {
        x >>= 1;
        x != 0
    } {
        lg += 1;
    }
    lg
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1270-1306
// ```c
// if (mb_lt->lut_word_length <= 12) {
//     if (mb_lt->hashsize <= 8 * kTargetPVSize)
//         pv_size = (Int4)(mb_lt->hashsize >> PV_ARRAY_BTS);
//     else
//         pv_size = kTargetPVSize / PV_ARRAY_BYTES;
// } else {
//     pv_size = kTargetPVSize * 64 / PV_ARRAY_BYTES;
// }
// if(!lookup_options->db_filter &&
//    (approx_table_entries <= kSmallQueryCutoff ||
//     approx_table_entries >= kLargeQueryCutoff)) {
//     pv_size = pv_size / 2;
// }
// mb_lt->pv_array_bts = ilog2(mb_lt->hashsize / pv_size);
// ```
fn compute_mb_pv_params(
    hashsize: usize,
    approx_table_entries: usize,
    db_filter: bool,
    lut_word_length: usize,
) -> (usize, usize) {
    const K_TARGET_PV_SIZE: usize = 131_072;
    const K_SMALL_QUERY_CUTOFF: usize = 15_000;
    const K_LARGE_QUERY_CUTOFF: usize = 800_000;

    let mut pv_size = if lut_word_length <= 12 {
        if hashsize <= 8 * K_TARGET_PV_SIZE {
            hashsize >> PV_ARRAY_BTS
        } else {
            K_TARGET_PV_SIZE / PV_ARRAY_BYTES
        }
    } else {
        K_TARGET_PV_SIZE * 64 / PV_ARRAY_BYTES
    };

    if !db_filter
        && (approx_table_entries <= K_SMALL_QUERY_CUTOFF
            || approx_table_entries >= K_LARGE_QUERY_CUTOFF)
    {
        pv_size = pv_size / 2;
    }

    if pv_size == 0 {
        pv_size = 1;
    }
    let pv_array_bts = ilog2(hashsize / pv_size);
    (pv_size, pv_array_bts)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:829-901
// ```c
// if ((val & BLAST2NA_MASK) != 0) { ecode = 0; pos = seq + kLutWordLength; continue; }
// ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
// if (seq < pos) continue;
// PV_SET(pv_array, ecode, pv_array_bts);
// ```
fn build_query_pv(
    queries: &[fasta::Record],
    query_masks: &[Vec<MaskedInterval>],
    lut_word_length: usize,
    approx_table_entries: usize,
) -> (Vec<PvArrayType>, usize) {
    if lut_word_length == 0 {
        return (Vec::new(), PV_ARRAY_BTS);
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1250-1254
    // ```c
    // mb_lt->hashsize = 1ULL << (BITS_PER_NUC * mb_lt->lut_word_length);
    // ```
    let hashsize = 1usize << (2 * lut_word_length);
    let (pv_size, pv_array_bts) = compute_mb_pv_params(hashsize, approx_table_entries, true, lut_word_length);
    let mut pv = vec![0u32; pv_size];

    let kmer_mask: u64 = (1u64 << (2 * lut_word_length)) - 1;
    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < lut_word_length {
            continue;
        }
        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let mut current_kmer: u64 = 0;
        let mut valid_bases: usize = 0;

        for pos in 0..seq.len() {
            let base = seq[pos];
            let code = ENCODE_LUT[base as usize];
            if code == 0xFF {
                current_kmer = 0;
                valid_bases = 0;
                continue;
            }

            current_kmer = ((current_kmer << 2) | (code as u64)) & kmer_mask;
            valid_bases += 1;

            if valid_bases < lut_word_length {
                continue;
            }

            let kmer_start = pos + 1 - lut_word_length;
            if !masks.is_empty() && is_kmer_masked(masks, kmer_start, lut_word_length) {
                continue;
            }

            pv_set_shift(&mut pv, current_kmer as usize, pv_array_bts);
        }
    }

    (pv, pv_array_bts)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1047-1059
// ```c
// if (!(ecode & 1)) {
//     if ((counts[ecode / 2] >> 4) >= max_word_count) continue;
// } else {
//     if ((counts[ecode / 2] & 0xf) >= max_word_count) continue;
// }
// ```
#[inline(always)]
fn db_word_count_exceeds(counts: &[u8], word: u64, max_word_count: u8) -> bool {
    let idx = word as usize;
    let byte_idx = idx >> 1;
    if byte_idx >= counts.len() {
        return false;
    }
    let count = if (idx & 1) == 1 {
        counts[byte_idx] & 0x0f
    } else {
        counts[byte_idx] >> 4
    };
    count >= max_word_count
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1166-1177
// ```c
// index = word / 2;
// if (word & 1) {
//     if ((counts[index] & 0xf) < max_word_count) counts[index]++;
// } else {
//     if ((counts[index] >> 4) < max_word_count) counts[index] += 1 << 4;
// }
// ```
#[inline(always)]
fn db_word_count_increment(counts: &mut [u8], word: u64, max_word_count: u8) {
    let idx = word as usize;
    let byte_idx = idx >> 1;
    if byte_idx >= counts.len() {
        return;
    }
    if (idx & 1) == 1 {
        if (counts[byte_idx] & 0x0f) < max_word_count {
            counts[byte_idx] = counts[byte_idx].wrapping_add(1);
        }
    } else if (counts[byte_idx] >> 4) < max_word_count {
        counts[byte_idx] = counts[byte_idx].wrapping_add(1 << 4);
    }
}

/// Optimized lookup table with Presence-Vector for fast filtering.
/// Combines DirectKmerLookup with a bit vector for O(1) presence checking.
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:251-260
// ```c
// Int4* hashtable;   /**< Array of positions              */
// Int4* next_pos;    /**< Extra positions stored here     */
// PV_ARRAY_TYPE *pv_array;/**< Presence vector, used for quick presence
//                            check */
// Int4 pv_array_bts; /**< The exponent of 2 by which pv_array is smaller than
//                        the backbone */
// ```
pub struct PvDirectLookup {
    /// The actual lookup table storing 1-based query offsets (q_off + 1)
    lookup: DirectKmerLookup,
    /// Presence vector - bit i is set if lookup[i] is non-empty
    pv: Vec<PvArrayType>,
    /// Log2 compression factor for the PV array (pv_array_bts)
    pv_array_bts: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1088-1108
    // ```c
    // longest_chain = 2;
    // for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
    //     longest_chain = MAX(longest_chain, helper_array[index]);
    // mb_lt->longest_chain = longest_chain;
    // ```
    /// Longest chain length for any lookup bucket (used to size offset buffers).
    longest_chain: usize,
    /// Word size used for this lookup table
    #[allow(dead_code)]
    word_size: usize,
}

impl PvDirectLookup {
    /// Check if a k-mer has any hits using the presence vector (O(1))
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:51-57
    // ```c
    // #define PV_TEST(lookup, index, shift) \
    //     ( lookup[(index) >> (shift)] & ((PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)) )
    // ```
    #[inline(always)]
    pub fn has_hits(&self, kmer: u64) -> bool {
        pv_test_shift(&self.pv, kmer as usize, self.pv_array_bts)
    }

    /// Get hits for a k-mer (only call after has_hits returns true)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
    // ```c
    // Int4 q_off = lookup->hashtable[index];
    // while (q_off) {
    //     offset_pairs[i].qs_offsets.q_off   = q_off - 1;
    //     offset_pairs[i++].qs_offsets.s_off = s_off;
    //     q_off = lookup->next_pos[q_off];
    // }
    // ```
    #[inline(always)]
    pub fn get_hits(&self, kmer: u64) -> &[u32] {
        let idx = kmer as usize;
        self.lookup.get(idx)
    }

    /// Get hits for a k-mer with PV check (combined operation)
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1386-1394
    // ```c
    // if (PV_TEST(pv, index, pv_array_bts))
    //     return 1;
    // else
    //     return 0;
    // ```
    #[inline(always)]
    pub fn get_hits_checked(&self, kmer: u64) -> &[u32] {
        if self.has_hits(kmer) {
            self.get_hits(kmer)
        } else {
            &[]
        }
    }

    /// Get the maximum number of hits for any lookup bucket.
    #[inline(always)]
    pub fn longest_chain(&self) -> usize {
        self.longest_chain
    }
}

/// Two-stage lookup table (like NCBI BLAST)
/// - lut_word_length: Used for indexing (e.g., 8 for megablast)
/// - word_length: Used for extension triggering (e.g., 28 for megablast)
/// This allows O(1) direct array access even for large word_length values
pub struct TwoStageLookup {
    /// Lookup table using lut_word_length
    pv_lookup: PvDirectLookup,
    /// Lookup word length (for indexing)
    lut_word_length: usize,
    /// Extension word length (for triggering extension)
    word_length: usize,
}

impl TwoStageLookup {
    /// Check if a lut_word_length k-mer has any hits (O(1))
    #[inline(always)]
    pub fn has_hits(&self, lut_kmer: u64) -> bool {
        self.pv_lookup.has_hits(lut_kmer)
    }

    /// Get hits for a lut_word_length k-mer
    #[inline(always)]
    pub fn get_hits(&self, lut_kmer: u64) -> &[u32] {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1374-1418
        // ```c
        // if (PV_TEST(pv, index, pv_array_bts))
        //     return 1;
        // ...
        // Int4 q_off = lookup->hashtable[index];
        // while (q_off) {
        //     offset_pairs[i].qs_offsets.q_off   = q_off - 1;
        //     offset_pairs[i++].qs_offsets.s_off = s_off;
        //     q_off = lookup->next_pos[q_off];
        // }
        // ```
        self.pv_lookup.get_hits_checked(lut_kmer)
    }

    /// Get lookup word length
    #[inline(always)]
    pub fn lut_word_length(&self) -> usize {
        self.lut_word_length
    }

    /// Get extension word length
    #[inline(always)]
    pub fn word_length(&self) -> usize {
        self.word_length
    }

    /// Calculate optimal scan step
    #[inline(always)]
    pub fn scan_step(&self) -> usize {
        (self.word_length as isize - self.lut_word_length as isize + 1).max(1) as usize
    }

    /// Get the longest chain length for sizing offset buffers.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/lookup_wrap.c:255-288
    // ```c
    // switch (lookup->lut_type) {
    // case eMBLookupTable:
    //     offset_array_size = OFFSET_ARRAY_SIZE +
    //         ((BlastMBLookupTable*)lookup->lut)->longest_chain;
    //     break;
    // ...
    // }
    // ```
    #[inline(always)]
    pub fn longest_chain(&self) -> usize {
        self.pv_lookup.longest_chain()
    }
}

/// Pack (q_idx, diag) into a single u64 key for faster HashMap operations
#[inline(always)]
pub fn pack_diag_key(q_idx: u32, diag: isize) -> u64 {
    ((q_idx as u64) << 32) | ((diag as i32) as u32 as u64)
}

/// Check if a k-mer starting at position overlaps with any masked interval
/// Uses binary search for O(log n) performance instead of O(n) linear scan.
/// IMPORTANT: Intervals must be sorted by start position for binary search to work correctly.
#[inline]
pub fn is_kmer_masked(intervals: &[MaskedInterval], start: usize, kmer_len: usize) -> bool {
    if intervals.is_empty() {
        return false;
    }

    let end = start + kmer_len;

    // Binary search to find the first interval whose end > start
    // (i.e., the first interval that could potentially overlap)
    let idx = intervals.partition_point(|interval| interval.end <= start);

    // Check if this interval actually overlaps
    // An interval overlaps if: interval.start < end AND interval.end > start
    // We already know interval.end > start (from binary search), so just check start
    if idx < intervals.len() && intervals[idx].start < end {
        return true;
    }

    false
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:91-132
// ```c
// for (loc = locations; loc; loc = loc->next) {
//     Int4 from = loc->ssr->left;
//     Int4 to = loc->ssr->right;
//     if (word_length > to - from + 1) continue;
//     ...
// }
// ```
pub(crate) fn build_unmasked_ranges(seq_len: usize, masks: &[MaskedInterval]) -> Vec<(usize, usize)> {
    if seq_len == 0 {
        return Vec::new();
    }
    if masks.is_empty() {
        return vec![(0, seq_len)];
    }

    let mut sorted = masks.to_vec();
    sorted.sort_by_key(|m| m.start);
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    let mut cursor = 0usize;
    for mask in sorted {
        let start = mask.start.min(seq_len);
        let end = mask.end.min(seq_len);
        if start > cursor {
            ranges.push((cursor, start));
        }
        if end > cursor {
            cursor = end;
        }
    }
    if cursor < seq_len {
        ranges.push((cursor, seq_len));
    }
    ranges
}

/// Lookup table for ASCII to 2-bit encoding (0xFF = invalid/ambiguous)
/// Used for rolling k-mer extraction
const ENCODE_LUT: [u8; 256] = {
    let mut lut = [0xFFu8; 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'c' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'g' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b't' as usize] = 3;
    lut[b'U' as usize] = 3;
    lut[b'u' as usize] = 3;
    lut
};

/// Build optimized lookup table with Presence-Vector using rolling k-mer extraction
/// This combines:
/// 1. O(1) sliding window k-mer extraction (Phase 1 optimization)
/// 2. Presence-Vector for fast filtering (Phase 2 optimization)
pub fn build_pv_direct_lookup(
    queries: &[fasta::Record],
    query_offsets: &[i32],
    word_size: usize,
    full_word_size: usize,
    query_masks: &[Vec<MaskedInterval>],
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    approx_table_entries: usize,
    use_mb_pv: bool,
) -> PvDirectLookup {
    let safe_word_size = word_size.min(MAX_DIRECT_LOOKUP_WORD_SIZE);
    let full_word_size = full_word_size.max(safe_word_size);
    let table_size = 1usize << (2 * safe_word_size); // 4^word_size
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1270-1306
    // ```c
    // if (mb_lt->lut_word_length <= 12) {
    //     if (mb_lt->hashsize <= 8 * kTargetPVSize)
    //         pv_size = (Int4)(mb_lt->hashsize >> PV_ARRAY_BTS);
    //     else
    //         pv_size = kTargetPVSize / PV_ARRAY_BYTES;
    // } else {
    //     pv_size = kTargetPVSize * 64 / PV_ARRAY_BYTES;
    // }
    // if(!lookup_options->db_filter &&
    //    (approx_table_entries <= kSmallQueryCutoff ||
    //     approx_table_entries >= kLargeQueryCutoff)) {
    //     pv_size = pv_size / 2;
    // }
    // mb_lt->pv_array_bts = ilog2(mb_lt->hashsize / pv_size);
    // ```
    let (pv_size, pv_array_bts) = if use_mb_pv {
        compute_mb_pv_params(
            table_size,
            approx_table_entries,
            db_word_counts.is_some(),
            safe_word_size,
        )
    } else {
        let pv_size = (table_size + PV_ARRAY_MASK) >> PV_ARRAY_BTS;
        (pv_size, PV_ARRAY_BTS)
    };
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    if debug_mode {
        let offsets_bytes = (table_size + 1) * std::mem::size_of::<u32>();
        let pv_bytes = pv_size * PV_ARRAY_BYTES;
        eprintln!(
            "[DEBUG] build_pv_direct_lookup: word_size={}, table_size={} ({:.1}MB offsets), pv_size={} ({:.1}KB), pv_array_bts={}",
            safe_word_size,
            table_size,
            offsets_bytes as f64 / 1_000_000.0,
            pv_size,
            pv_bytes as f64 / 1_000.0,
            pv_array_bts,
        );
    }

    let mut counts: Vec<u32> = vec![0; table_size];

    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;
    let mut dust_skipped = 0usize;

    debug_assert_eq!(queries.len(), query_offsets.len());

    // K-mer mask for rolling window
    let kmer_mask: u64 = (1u64 << (2 * safe_word_size)) - 1;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:979-1091
    // ```c
    // mb_lt->next_pos = (Int4 *)calloc(query->length + 1, sizeof(Int4));
    // ...
    // if (mb_lt->hashtable[ecode] == 0) {
    //     PV_SET(pv_array, ecode, pv_array_bts);
    // }
    // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    // mb_lt->hashtable[ecode] = index;
    // ```
    // NCBI reference: blast_lookup.c:BlastLookupIndexQueryExactMatches (lines 79-132)
    // NCBI processes only unmasked regions (locations parameter)
    // For each location, it iterates through positions and adds k-mers
    // Reference: blast_nalookup.c:402-406, 571-575 (calls BlastLookupIndexQueryExactMatches)
    //
    // LOSAT mirrors this by iterating unmasked ranges derived from query masks.
    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        // NCBI reference: blast_lookup.c:99-100
        // if (word_length > to - from + 1) continue;
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);

        // NCBI reference: blast_lookup.c:108-121
        // Rolling window approach: word_target points to position where complete k-mer can be formed
        // Ambiguous bases reset the window: if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
        for (range_start, range_end) in ranges {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1001-1034
            // ```c
            // if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
            //     continue;
            // ```
            let range_len = range_end.saturating_sub(range_start);
            if full_word_size > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases: usize = 0;

            for pos in range_start..range_end {
                let base = seq[pos];
                let code = ENCODE_LUT[base as usize];

                // NCBI reference: blast_lookup.c:119-120
                // if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
                if code == 0xFF {
                    valid_bases = 0;
                    current_kmer = 0;
                    ambiguous_skipped += 1;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | (code as u64)) & kmer_mask;
                valid_bases += 1;
                if valid_bases < safe_word_size {
                    continue;
                }

                total_positions += 1;

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1047-1059
                // ```c
                // if (kDbFilter) {
                //    if ((counts[ecode / 2] >> 4) >= max_word_count) continue;
                //    ...
                // }
                // ```
                if let Some(counts_filter) = db_word_counts {
                    if db_word_count_exceeds(counts_filter, current_kmer, max_db_word_count) {
                        continue;
                    }
                }

                // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
                // Adds ALL hits without any frequency limit - no query-side filtering
                // if (backbone[index] == NULL) { initialize new chain }
                // else { use existing chain, realloc if full }
                // chain[chain[1] + 2] = query_offset; chain[1]++;
                let idx = current_kmer as usize;
                if idx < table_size {
                    counts[idx] = counts[idx].saturating_add(1);
                }
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
    // ```c
    // Int4 q_off = lookup->hashtable[index];
    // while (q_off) {
    //     offset_pairs[i].qs_offsets.q_off   = q_off - 1;
    //     offset_pairs[i++].qs_offsets.s_off = s_off;
    //     q_off = lookup->next_pos[q_off];
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
    // ```c
    // Int4 q_off = lookup->hashtable[index];
    // while (q_off) {
    //     offset_pairs[i].qs_offsets.q_off   = q_off - 1;
    //     offset_pairs[i++].qs_offsets.s_off = s_off;
    //     q_off = lookup->next_pos[q_off];
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1088-1108
    // ```c
    // longest_chain = 2;
    // for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
    //     longest_chain = MAX(longest_chain, helper_array[index]);
    // mb_lt->longest_chain = longest_chain;
    // ```
    let longest_chain = counts
        .iter()
        .copied()
        .max()
        .unwrap_or(0)
        .max(2) as usize;

    let mut offsets: Vec<u32> = vec![0; table_size + 1];
    let mut total_hits: u32 = 0;
    for idx in 0..table_size {
        offsets[idx] = total_hits;
        total_hits = total_hits.saturating_add(counts[idx]);
    }
    offsets[table_size] = total_hits;

    let mut hits: Vec<u32> = vec![0u32; total_hits as usize];
    let mut write_pos: Vec<u32> = offsets[..table_size].to_vec();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:33-77
    // ```c
    // if (backbone[index] == NULL) { ... }
    // ...
    // chain[chain[1] + 2] = query_offset;
    // chain[1]++;
    // ```
    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);
        let query_offset = query_offsets[q_idx] as usize;

        for (range_start, range_end) in ranges {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1001-1034
            // ```c
            // if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
            //     continue;
            // ```
            let range_len = range_end.saturating_sub(range_start);
            if full_word_size > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases: usize = 0;

            for pos in range_start..range_end {
                let base = seq[pos];
                let code = ENCODE_LUT[base as usize];
                if code == 0xFF {
                    valid_bases = 0;
                    current_kmer = 0;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | (code as u64)) & kmer_mask;
                valid_bases += 1;
                if valid_bases < safe_word_size {
                    continue;
                }

                let kmer_start = pos + 1 - safe_word_size;
                if let Some(counts_filter) = db_word_counts {
                    if db_word_count_exceeds(counts_filter, current_kmer, max_db_word_count) {
                        continue;
                    }
                }

                let idx = current_kmer as usize;
                if idx < table_size {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
                    // ```c
                    // /* Also add 1 to all indices, because lookup table indices count
                    //    from 1. */
                    // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
                    // mb_lt->hashtable[ecode] = index;
                    // ```
                    let q_off_1 = (query_offset + kmer_start + 1) as u32;
                    let pos_idx = write_pos[idx] as usize;
                    hits[pos_idx] = q_off_1;
                    write_pos[idx] = write_pos[idx].saturating_add(1);
                }
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:51-57
    // ```c
    // #define PV_SET(lookup, index, shift) \
    //     lookup[(index) >> (shift)] |= (PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)
    // ```
    let mut pv: Vec<PvArrayType> = vec![0; pv_size];
    let mut non_empty_count = 0usize;
    for idx in 0..table_size {
        if offsets[idx] != offsets[idx + 1] {
            pv_set_shift(&mut pv, idx, pv_array_bts);
            non_empty_count += 1;
        }
    }

    if debug_mode {
        let hits_bytes = hits.len() * std::mem::size_of::<u32>();
        eprintln!(
            "[DEBUG] build_pv_direct_lookup: total_positions={}, ambiguous_skipped={} ({:.1}%), dust_skipped={} ({:.1}%)",
            total_positions,
            ambiguous_skipped,
            100.0 * ambiguous_skipped as f64 / (total_positions + ambiguous_skipped).max(1) as f64,
            dust_skipped,
            100.0 * dust_skipped as f64 / total_positions.max(1) as f64
        );
        eprintln!(
            "[DEBUG] build_pv_direct_lookup: kmers_added={}, non_empty_buckets={}, hits_bytes={:.1}MB",
            total_hits,
            non_empty_count,
            hits_bytes as f64 / 1_000_000.0
        );
    }

    PvDirectLookup {
        lookup: DirectKmerLookup { offsets, hits },
        pv,
        pv_array_bts,
        longest_chain,
        word_size: safe_word_size,
    }
}

pub fn build_lookup(
    queries: &[fasta::Record],
    query_offsets: &[i32],
    word_size: usize,
    query_masks: &[Vec<MaskedInterval>],
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
) -> KmerLookup {
    let mut lookup: KmerLookup = FxHashMap::default();
    let safe_word_size = word_size.min(31);
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;
    let mut dust_skipped = 0usize;

    debug_assert_eq!(queries.len(), query_offsets.len());

    // NCBI reference: blast_lookup.c:BlastLookupIndexQueryExactMatches (lines 79-132)
    // NCBI processes only unmasked regions (locations parameter)
    // Reference: blast_nalookup.c:402-406, 571-575 (calls BlastLookupIndexQueryExactMatches)
    //
    // LOSAT processes entire sequence and filters masked regions (equivalent behavior)
    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        // NCBI reference: blast_lookup.c:99-100
        // if (word_length > to - from + 1) continue;
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);

        // NCBI reference: blast_lookup.c:108-121
        // for (offset = from; offset <= to; offset++, seq++) {
        //     if (seq >= word_target) { BlastLookupAddWordHit(...); }
        //     if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
        // }
        for i in 0..=(seq.len() - safe_word_size) {
            total_positions += 1;
            
            // NCBI reference: blast_nalookup.c:402-406
            // BlastLookupIndexQueryExactMatches processes only unmasked regions (locations)
            // LOSAT processes entire sequence and filters masked regions (equivalent)
            // Skip k-mers that overlap with DUST-masked regions
            if !masks.is_empty() && is_kmer_masked(masks, i, safe_word_size) {
                dust_skipped += 1;
                continue;
            }
            
            // NCBI reference: blast_lookup.c:119-120
            // if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
            if let Some(kmer) = encode_kmer(seq, i, safe_word_size) {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1047-1059
                // ```c
                // if (kDbFilter) {
                //    if ((counts[ecode / 2] >> 4) >= max_word_count) continue;
                //    ...
                // }
                // ```
                if let Some(counts) = db_word_counts {
                    if db_word_count_exceeds(counts, kmer, max_db_word_count) {
                        continue;
                    }
                }
                // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
                // Adds ALL hits without any frequency limit - no query-side filtering
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
                // ```c
                // /* Also add 1 to all indices, because lookup table indices count
                //    from 1. */
                // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
                // mb_lt->hashtable[ecode] = index;
                // ```
                let q_off_1 = (query_offsets[q_idx] as usize + i + 1) as u32;
                lookup.entry(kmer).or_default().push(q_off_1);
            } else {
                ambiguous_skipped += 1;
            }
        }
    }

    if debug_mode {
        eprintln!(
            "[DEBUG] build_lookup: total_positions={}, ambiguous_skipped={} ({:.1}%), dust_skipped={} ({:.1}%), unique_kmers={}",
            total_positions,
            ambiguous_skipped,
            100.0 * ambiguous_skipped as f64 / total_positions.max(1) as f64,
            dust_skipped,
            100.0 * dust_skipped as f64 / total_positions.max(1) as f64,
            lookup.len()
        );
    }

    // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
    // NCBI BLAST does NOT filter over-represented k-mers in the query
    // All k-mers are added to the lookup table regardless of frequency
    // Database word count filtering (kDbFilter) exists but is different:
    // it filters based on database counts, not query counts
    // Reference: blast_nalookup.c:1047-1060 (database word count filtering)
    //
    // REMOVED: Over-represented k-mer filtering (MAX_HITS_PER_KMER) - does not exist in NCBI BLAST

    lookup
}

/// Build two-stage lookup table (like NCBI BLAST)
/// - lut_word_length: Used for indexing (typically 8 for megablast)
/// - word_length: Used for extension triggering (typically 28 for megablast)
/// This allows O(1) direct array access even for large word_length values
pub fn build_two_stage_lookup(
    queries: &[fasta::Record],
    query_offsets: &[i32],
    word_length: usize,
    lut_word_length: usize,
    query_masks: &[Vec<MaskedInterval>],
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    approx_table_entries: usize,
) -> TwoStageLookup {
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();
    
    if debug_mode {
        eprintln!(
            "[DEBUG] build_two_stage_lookup: word_length={}, lut_word_length={}",
            word_length, lut_word_length
        );
    }

    // Build the lookup table using lut_word_length
    // We index all positions where a full word_length match is possible
    let pv_lookup = build_pv_direct_lookup(
        queries,
        query_offsets,
        lut_word_length,
        word_length,
        query_masks,
        db_word_counts,
        max_db_word_count,
        approx_table_entries,
        true,
    );
    
    // Note: The actual word_length matching is done during scanning,
    // where we check if the subject sequence has a word_length match
    // starting at the position indicated by the lut_word_length lookup
    
    TwoStageLookup {
        pv_lookup,
        lut_word_length,
        word_length,
    }
}

/// Build direct address table for k-mer lookup (O(1) access)
/// This is much faster than HashMap for small word sizes
pub fn build_direct_lookup(
    queries: &[fasta::Record],
    query_offsets: &[i32],
    word_size: usize,
    query_masks: &[Vec<MaskedInterval>],
) -> DirectKmerLookup {
    let safe_word_size = word_size.min(MAX_DIRECT_LOOKUP_WORD_SIZE);
    let table_size = 1usize << (2 * safe_word_size); // 4^word_size
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    if debug_mode {
        let offsets_bytes = (table_size + 1) * std::mem::size_of::<u32>();
        eprintln!(
            "[DEBUG] build_direct_lookup: word_size={}, table_size={} ({:.1}MB offsets)",
            safe_word_size,
            table_size,
            offsets_bytes as f64 / 1_000_000.0
        );
    }

    let mut counts: Vec<u32> = vec![0; table_size];

    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;
    let mut dust_skipped = 0usize;

    debug_assert_eq!(queries.len(), query_offsets.len());

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:979-1091
    // ```c
    // mb_lt->next_pos = (Int4 *)calloc(query->length + 1, sizeof(Int4));
    // ...
    // if (mb_lt->hashtable[ecode] == 0) {
    //     PV_SET(pv_array, ecode, pv_array_bts);
    // }
    // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    // mb_lt->hashtable[ecode] = index;
    // ```
    // NCBI reference: blast_lookup.c:BlastLookupIndexQueryExactMatches (lines 79-132)
    // NCBI processes only unmasked regions (locations parameter)
    // Reference: blast_nalookup.c:402-406, 571-575 (calls BlastLookupIndexQueryExactMatches)
    //
    // LOSAT processes entire sequence and filters masked regions (equivalent behavior)
    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        // NCBI reference: blast_lookup.c:99-100
        // if (word_length > to - from + 1) continue;
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);

        // NCBI reference: blast_lookup.c:108-121
        // for (offset = from; offset <= to; offset++, seq++) {
        //     if (seq >= word_target) { BlastLookupAddWordHit(...); }
        //     if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
        // }
        for i in 0..=(seq.len() - safe_word_size) {
            total_positions += 1;

            // NCBI reference: blast_nalookup.c:402-406
            // BlastLookupIndexQueryExactMatches processes only unmasked regions (locations)
            // LOSAT processes entire sequence and filters masked regions (equivalent)
            // Skip k-mers that overlap with DUST-masked regions
            if !masks.is_empty() && is_kmer_masked(masks, i, safe_word_size) {
                dust_skipped += 1;
                continue;
            }

            // NCBI reference: blast_lookup.c:119-120
            // if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
            if let Some(kmer) = encode_kmer(seq, i, safe_word_size) {
                // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
                // Adds ALL hits without any frequency limit - no query-side filtering
                let idx = kmer as usize;
                if idx < table_size {
                    counts[idx] = counts[idx].saturating_add(1);
                }
            } else {
                ambiguous_skipped += 1;
            }
        }
    }

    let mut offsets: Vec<u32> = vec![0; table_size + 1];
    let mut total_hits: u32 = 0;
    for idx in 0..table_size {
        offsets[idx] = total_hits;
        total_hits = total_hits.saturating_add(counts[idx]);
    }
    offsets[table_size] = total_hits;

    let mut hits: Vec<u32> = vec![0u32; total_hits as usize];
    let mut write_pos: Vec<u32> = offsets[..table_size].to_vec();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:33-77
    // ```c
    // if (backbone[index] == NULL) { ... }
    // ...
    // chain[chain[1] + 2] = query_offset;
    // chain[1]++;
    // ```
    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let query_offset = query_offsets[q_idx] as usize;
        for i in 0..=(seq.len() - safe_word_size) {
            if !masks.is_empty() && is_kmer_masked(masks, i, safe_word_size) {
                continue;
            }

            if let Some(kmer) = encode_kmer(seq, i, safe_word_size) {
                let idx = kmer as usize;
                if idx < table_size {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
                    // ```c
                    // /* Also add 1 to all indices, because lookup table indices count
                    //    from 1. */
                    // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
                    // mb_lt->hashtable[ecode] = index;
                    // ```
                    let q_off_1 = (query_offset + i + 1) as u32;
                    let pos_idx = write_pos[idx] as usize;
                    hits[pos_idx] = q_off_1;
                    write_pos[idx] = write_pos[idx].saturating_add(1);
                }
            }
        }
    }

    if debug_mode {
        let non_empty = offsets
            .windows(2)
            .filter(|w| w[0] != w[1])
            .count();
        let hits_bytes = hits.len() * std::mem::size_of::<u32>();
        eprintln!(
            "[DEBUG] build_direct_lookup: total_positions={}, ambiguous_skipped={} ({:.1}%), dust_skipped={} ({:.1}%), non_empty_buckets={}, hits_bytes={:.1}MB",
            total_positions,
            ambiguous_skipped,
            100.0 * ambiguous_skipped as f64 / total_positions.max(1) as f64,
            dust_skipped,
            100.0 * dust_skipped as f64 / total_positions.max(1) as f64,
            non_empty,
            hits_bytes as f64 / 1_000_000.0
        );
    }

    // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
    // NCBI BLAST does NOT filter over-represented k-mers in the query
    // All k-mers are added to the lookup table regardless of frequency
    // Database word count filtering (kDbFilter) exists but is different:
    // it filters based on database counts, not query counts
    // Reference: blast_nalookup.c:1047-1060 (database word count filtering)
    //
    // REMOVED: Over-represented k-mer filtering (MAX_HITS_PER_KMER) - does not exist in NCBI BLAST

    DirectKmerLookup { offsets, hits }
}


