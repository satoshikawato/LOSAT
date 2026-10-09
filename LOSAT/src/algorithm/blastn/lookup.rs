use super::constants::MAX_DIRECT_LOOKUP_WORD_SIZE;
use super::disc_lookup::{
    compute_discontiguous_index, get_disc_template_type, DiscTemplateType, DiscWordType,
};
use crate::blastinput::fasta_reader::InputRecord;
use crate::core::blast_encoding::{encode_subject_ncbi2na_packed, COMPRESSION_RATIO};
use crate::utils::dust::MaskedInterval;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:37-43
// ```c
// /** bitfield used to detect ambiguities in uncompressed
//  *  nucleotide letters
//  */
// #define BLAST2NA_MASK 0xfc
// #define BITS_PER_NUC 2
// ```
const BLAST2NA_MASK: u8 = 0xFC;

/// 2-bit encoding for compact storage and hashing from BLASTNA-encoded sequence.
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:96-105
/// ```c
/// for(i = 0; i < wordsize; i++) {
///   index = (index << charsize) | word[i];
/// }
/// ```
pub fn encode_kmer(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    if start + k > seq.len() {
        return None;
    }
    let mut encoded: u64 = 0;
    for i in 0..k {
        let base = unsafe { *seq.get_unchecked(start + i) };
        if (base & BLAST2NA_MASK) != 0 {
            return None;
        }
        encoded = (encoded << 2) | base as u64;
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
    b'-', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:95-103
// ```c
// const Uint1 IUPACNA_TO_NCBI4NA[128]={ ... };
// ```
const IUPACNA_TO_NCBI4NA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, 0, 0, 5, 6, 8, 0, 7, 9, 0, 10, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,
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
const NCBI4NA_REV_COMP: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

/// Generate the reverse complement of a DNA sequence in IUPACNA.
/// Uses NCBI4NA bitmask mapping to preserve ambiguous base complements.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:977,989
/// ```c
///     sv.SetCoding(CSeq_data::e_Ncbi4na);
///     ...
///     sv.GetStrandData(strand, buffer);
/// ```
/// NCBI takes both strands, in NCBI4NA, from sequence data that has no letter case
/// (the reader records lowercase as masks), so a lowercase input letter is the same
/// base as its uppercase form. `IUPACNA_TO_NCBI4NA` maps only uppercase letters, so
/// the letter is uppercased first; the result is uppercase.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| {
            let b = b.to_ascii_uppercase();
            let idx = if b < 128 {
                IUPACNA_TO_NCBI4NA[b as usize]
            } else {
                0
            };
            let rc = NCBI4NA_REV_COMP[idx as usize];
            NCBI4NA_TO_IUPACNA[rc as usize]
        })
        .collect()
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:109-156
// ```c
// #define NA_HITS_PER_CELL 3
// typedef struct NaLookupBackboneCell { ... } NaLookupBackboneCell;
// typedef struct BlastNaLookupTable { ... } BlastNaLookupTable;
// ```
const NA_HITS_PER_CELL: usize = 3;

#[repr(C)]
#[derive(Clone, Copy)]
union NaLookupPayload {
    entries: [u32; NA_HITS_PER_CELL],
    overflow_cursor: u32,
}

impl Default for NaLookupPayload {
    fn default() -> Self {
        NaLookupPayload { overflow_cursor: 0 }
    }
}

#[repr(C)]
#[derive(Clone, Copy)]
struct NaLookupBackboneCell {
    num_used: u32,
    payload: NaLookupPayload,
}

impl Default for NaLookupBackboneCell {
    fn default() -> Self {
        NaLookupBackboneCell {
            num_used: 0,
            payload: NaLookupPayload::default(),
        }
    }
}

/// Compact array-backed lookup table for blastn (NCBI BlastNaLookupTable).
pub struct NaLookupTable {
    backbone: Vec<NaLookupBackboneCell>,
    overflow: Vec<u32>,
    pv: Vec<PvArrayType>,
    longest_chain: usize,
    word_length: usize,
    lut_word_length: usize,
}

impl NaLookupTable {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:41-79
    // ```c
    // if (PV_TEST(pv, index, PV_ARRAY_BTS))
    //     return lookup->thick_backbone[index].num_used;
    // ...
    // if (num_hits <= NA_HITS_PER_CELL)
    //     lookup_pos = lookup->thick_backbone[index].payload.entries;
    // else
    //     lookup_pos = lookup->overflow + lookup->thick_backbone[index].payload.overflow_cursor;
    // ```
    #[inline(always)]
    pub fn get_hits_checked(&self, idx: u64) -> &[u32] {
        let idx = idx as usize;
        if idx >= self.backbone.len() {
            return &[];
        }
        if !pv_test(&self.pv, idx) {
            return &[];
        }
        let cell = &self.backbone[idx];
        let num_hits = cell.num_used as usize;
        if num_hits == 0 {
            return &[];
        }
        unsafe {
            // SAFETY: The payload field is accessed according to num_used,
            // matching the NCBI union layout (entries for small chains,
            // overflow_cursor for larger chains).
            if num_hits <= NA_HITS_PER_CELL {
                &cell.payload.entries[..num_hits]
            } else {
                let start = cell.payload.overflow_cursor as usize;
                &self.overflow[start..start + num_hits]
            }
        }
    }

    #[inline(always)]
    pub fn longest_chain(&self) -> usize {
        self.longest_chain
    }
}

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
pub fn build_db_word_counts<R: InputRecord>(
    queries_blastna: &[Vec<u8>],
    query_masks: &[Vec<MaskedInterval>],
    subjects: &[R],
    lut_word_length: usize,
    full_word_size: usize,
    max_word_count: u8,
    approx_table_entries: usize,
    subjects_packed: Option<&[Vec<u8>]>,
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:865-868
    // ```c
    // if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
    //     continue;
    // ```
    let (pv, pv_array_bts) = build_query_pv(
        queries_blastna,
        query_masks,
        lut_word_length,
        full_word_size,
        approx_table_entries,
    );

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1122-1177
    // ```c
    // word = (w >> shift) & mask;
    // if (!PV_TEST(pv, word, pv_array_bts)) continue;
    // if ((counts[index] & 0xf) < max_word_count) counts[index]++;
    // ```
    let mut scan_subject = |seq_len: usize, packed: &[u8]| {
        let mask = (1u64 << (2 * lut_word_length)) - 1;
        let mut pos = 0usize;
        let end = seq_len - lut_word_length;
        let mut kmer = packed_kmer_at(packed, 0, lut_word_length);

        loop {
            if pv_test_shift(&pv, kmer as usize, pv_array_bts) {
                db_word_count_increment(&mut counts, kmer, max_word_count);
            }

            if pos == end {
                break;
            }

            let next_base = packed_base_at(packed, pos + lut_word_length);
            kmer = ((kmer << 2) | next_base as u64) & mask;
            pos += 1;
        }
    };

    for (subject_idx, record) in subjects.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < lut_word_length {
            continue;
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // ```c
        // if (subj_is_na) {
        //     BlastSeqBlkSetSequence(subj, sequence.data.release(),
        //        ((sentinels == eSentinels) ? sequence.length - 2 :
        //         sequence.length));
        //     ...
        //     SBlastSequence compressed_seq =
        //         subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
        //                                   eNa_strand_plus, eNoSentinels);
        //     BlastSeqBlkSetCompressedSequence(subj,
        //                               compressed_seq.data.release());
        // }
        // ```
        if let Some(packed_cache) = subjects_packed {
            if let Some(packed) = packed_cache.get(subject_idx) {
                scan_subject(seq.len(), packed.as_slice());
                continue;
            }
        }

        let packed = encode_subject_ncbi2na_packed(seq);
        scan_subject(seq.len(), &packed);
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:878-893
// ```c
// val = *++seq;
// if ((val & BLAST2NA_MASK) != 0) {
//     ecode = 0;
//     pos = seq + kLutWordLength;
//     continue;
// }
// ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
// if (seq < pos) continue;
// PV_SET(pv_array, ecode, pv_array_bts);
// ```
fn build_query_pv(
    queries_blastna: &[Vec<u8>],
    query_masks: &[Vec<MaskedInterval>],
    lut_word_length: usize,
    full_word_size: usize,
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
    let (pv_size, pv_array_bts) =
        compute_mb_pv_params(hashsize, approx_table_entries, true, lut_word_length);
    let mut pv = vec![0u32; pv_size];

    let kmer_mask: u64 = (1u64 << (2 * lut_word_length)) - 1;
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
        if seq.len() < lut_word_length {
            continue;
        }
        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:856-889
        // ```c
        // for (loc = location; loc; loc = loc->next) {
        //     if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
        //         continue;
        //     ...
        //     if ((val & BLAST2NA_MASK) != 0) {
        //         ecode = 0;
        //         pos = seq + kLutWordLength;
        //         continue;
        //     }
        //     ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
        //     if (seq < pos) continue;
        //     PV_SET(pv_array, ecode, pv_array_bts);
        // }
        // ```
        for (range_start, range_end) in ranges {
            let range_len = range_end.saturating_sub(range_start);
            if full_word_size > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases: usize = 0;

            for pos in range_start..range_end {
                let base = seq[pos];
                if (base & BLAST2NA_MASK) != 0 {
                    current_kmer = 0;
                    valid_bases = 0;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | (base as u64)) & kmer_mask;
                valid_bases += 1;

                if valid_bases < lut_word_length {
                    continue;
                }

                pv_set_shift(&mut pv, current_kmer as usize, pv_array_bts);
            }
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
    /// Megablast lookup table using NCBI hashtable/next_pos chains.
    mb_lookup: MbLookupTable,
    /// Lookup word length (for indexing)
    lut_word_length: usize,
    /// Extension word length (for triggering extension)
    word_length: usize,
    /// The discontiguous templates of a discontiguous megablast table (NCBI
    /// `mb_lt->discontiguous`); `None` for a contiguous table.
    disc: Option<DiscTemplates>,
}

/// The templates of a discontiguous megablast table and the second template's chains.
///
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:238-256
/// ```c
/// typedef struct BlastMBLookupTable {
///     ...
///     Boolean discontiguous; /**< Are discontiguous words used? */
///     Int4 template_length; /**< Length of the discontiguous word template */
///     EDiscTemplateType template_type; /**< Type of the discontiguous
///                                          word template */
///     Boolean two_templates; /**< Use two templates simultaneously */
///     EDiscTemplateType second_template_type; /**< Type of the second
///                                                 discontiguous word template */
///     ...
///     Int4* hashtable2;  /**< Array of positions for second template */
///     ...
///     Int4* next_pos2;   /**< Extra positions for the second template */
/// ```
pub struct DiscTemplates {
    pub template_type: DiscTemplateType,
    /// The second template (`two_templates`), or `Contiguous` for one template.
    pub second_template_type: DiscTemplateType,
    pub two_templates: bool,
    pub template_length: usize,
    hashtable2: Vec<u32>,
    next_pos2: Vec<u32>,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:238-269
// ```c
// typedef struct BlastMBLookupTable {
//     Int4 word_length;
//     Int4 lut_word_length;
//     Int8 hashsize;
//     ...
//     Int4* hashtable;
//     Int4* next_pos;
//     PV_ARRAY_TYPE *pv_array;
//     Int4 pv_array_bts;
//     Int4 longest_chain;
// } BlastMBLookupTable;
// ```
struct MbLookupTable {
    hashtable: Vec<u32>,
    next_pos: Vec<u32>,
    pv_array: Vec<PvArrayType>,
    pv_array_bts: usize,
    longest_chain: usize,
}

impl MbLookupTable {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1386-1395
    // ```c
    // if (PV_TEST(pv, index, pv_array_bts))
    //     return 1;
    // else
    //     return 0;
    // ```
    #[inline(always)]
    fn has_hits(&self, index: u64) -> bool {
        pv_test_shift(&self.pv_array, index as usize, self.pv_array_bts)
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1406-1418
    // ```c
    // Int4 q_off = lookup->hashtable[index];
    //
    // while (q_off) {
    //     offset_pairs[i].qs_offsets.q_off   = q_off - 1;
    //     offset_pairs[i++].qs_offsets.s_off = s_off;
    //     q_off = lookup->next_pos[q_off];
    // }
    // ```
    #[inline(always)]
    fn for_each_hit(&self, index: u64, mut callback: impl FnMut(u32)) {
        if !self.has_hits(index) {
            return;
        }
        let mut q_off = self.hashtable[index as usize];
        while q_off != 0 {
            callback(q_off);
            q_off = self.next_pos[q_off as usize];
        }
    }

    #[inline(always)]
    fn contains_hit(&self, index: u64, q_off_1: u32) -> bool {
        let mut found = false;
        self.for_each_hit(index, |hit_q_off| {
            if hit_q_off == q_off_1 {
                found = true;
            }
        });
        found
    }

    #[inline(always)]
    fn count_hits(&self, index: u64) -> usize {
        let mut count = 0usize;
        self.for_each_hit(index, |_| {
            count += 1;
        });
        count
    }

    #[inline(always)]
    fn longest_chain(&self) -> usize {
        self.longest_chain
    }
}

impl TwoStageLookup {
    /// Check if a lut_word_length k-mer has any hits (O(1))
    #[inline(always)]
    pub fn has_hits(&self, lut_kmer: u64) -> bool {
        self.mb_lookup.has_hits(lut_kmer)
    }

    /// The discontiguous templates, for a discontiguous megablast table.
    #[inline(always)]
    pub fn disc(&self) -> Option<&DiscTemplates> {
        self.disc.as_ref()
    }

    /// Visits the second template's query offsets (1-based) of a discontiguous word; the
    /// presence vector is shared by both templates.
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1465-1469
    /// ```c
    ///    if (s_BlastMBLookupHasHits(mb_lt, index2)) {         \
    ///        total_hits += s_BlastMBLookupRetrieve2(mb_lt,    \
    ///                   index2, offset_pairs + total_hits,    \
    ///                   scan_range[0]);                               \
    ///    }
    /// ```
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1428-1444
    /// ```c
    /// static NCBI_INLINE Int4 s_BlastMBLookupRetrieve2(BlastMBLookupTable * lookup,
    ///                                                  Int8 index,
    ///                                                  BlastOffsetPair * offset_pairs,
    ///                                                  Int4 s_off)
    /// {
    ///     Int4 i=0;
    ///     Int4 q_off = lookup->hashtable2[index];
    ///
    ///     while (q_off) {
    ///         offset_pairs[i].qs_offsets.q_off   = q_off - 1;
    ///         offset_pairs[i++].qs_offsets.s_off = s_off;
    ///         q_off = lookup->next_pos2[q_off];
    ///     }
    ///     return i;
    /// }
    /// ```
    #[inline(always)]
    pub fn for_each_hit2(&self, index: u64, mut callback: impl FnMut(u32)) {
        let Some(disc) = self.disc.as_ref() else {
            return;
        };
        if !self.mb_lookup.has_hits(index) {
            return;
        }
        let mut q_off = disc.hashtable2[index as usize];
        while q_off != 0 {
            callback(q_off);
            q_off = disc.next_pos2[q_off as usize];
        }
    }

    /// Visit hits for a lut_word_length k-mer.
    #[inline(always)]
    pub fn for_each_hit(&self, lut_kmer: u64, callback: impl FnMut(u32)) {
        self.mb_lookup.for_each_hit(lut_kmer, callback)
    }

    #[inline(always)]
    pub fn contains_hit(&self, lut_kmer: u64, q_off_1: u32) -> bool {
        self.mb_lookup.contains_hit(lut_kmer, q_off_1)
    }

    #[inline(always)]
    pub fn count_hits(&self, lut_kmer: u64) -> usize {
        self.mb_lookup.count_hits(lut_kmer)
    }

    /// Get lookup word length. A discontiguous table scans and extends template-length
    /// words.
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:1624-1634
    /// ```c
    ///     else if (lookup_wrap->lut_type == eMBLookupTable) {
    ///         BlastMBLookupTable *lookup =
    ///                                 (BlastMBLookupTable *) lookup_wrap->lut;
    ///         if (lookup->discontiguous) {
    ///             word_length = lookup->template_length;
    ///             lut_word_length = lookup->template_length;
    ///         } else {
    ///             word_length = lookup->word_length;
    ///             lut_word_length = lookup->lut_word_length;
    ///         }
    /// ```
    #[inline(always)]
    pub fn lut_word_length(&self) -> usize {
        match &self.disc {
            Some(disc) => disc.template_length,
            None => self.lut_word_length,
        }
    }

    /// Get extension word length (the template length for a discontiguous table).
    #[inline(always)]
    pub fn word_length(&self) -> usize {
        match &self.disc {
            Some(disc) => disc.template_length,
            None => self.word_length,
        }
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
        self.mb_lookup.longest_chain()
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
pub(crate) fn build_unmasked_ranges(
    seq_len: usize,
    masks: &[MaskedInterval],
) -> Vec<(usize, usize)> {
    if seq_len == 0 {
        return Vec::new();
    }
    if masks.is_empty() {
        return vec![(0, seq_len)];
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:121-126
    // ```c
    // const int kTopFlags = CSeq_loc::fStrand_Ignore|CSeq_loc::fMerge_All|CSeq_loc::fSort;
    // ...
    // query_masks->Merge(kTopFlags, 0);
    // ```
    // Dust/lowercase masks are merged and sorted before lookup building.
    debug_assert!(masks.windows(2).all(|w| w[0].start <= w[1].start));

    let mut ranges: Vec<(usize, usize)> = Vec::new();
    let mut cursor = 0usize;
    for mask in masks {
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

/// Build optimized lookup table with Presence-Vector using rolling k-mer extraction
/// This combines:
/// 1. O(1) sliding window k-mer extraction (Phase 1 optimization)
/// 2. Presence-Vector for fast filtering (Phase 2 optimization)
pub fn build_pv_direct_lookup(
    queries_blastna: &[Vec<u8>],
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

    debug_assert_eq!(queries_blastna.len(), query_offsets.len());

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
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
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

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:117-120
                // ```c
                // if (*seq & invalid_mask)
                //     word_target = seq + lut_word_length + 1;
                // ```
                if (base & BLAST2NA_MASK) != 0 {
                    valid_bases = 0;
                    current_kmer = 0;
                    ambiguous_skipped += 1;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | (base as u64)) & kmer_mask;
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
    let longest_chain = counts.iter().copied().max().unwrap_or(0).max(2) as usize;

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
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
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
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:117-120
                // ```c
                // if (*seq & invalid_mask)
                //     word_target = seq + lut_word_length + 1;
                // ```
                if (base & BLAST2NA_MASK) != 0 {
                    valid_bases = 0;
                    current_kmer = 0;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | (base as u64)) & kmer_mask;
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

/// Build standard NCBI BlastNaLookupTable (array-backed with PV + overflow).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:548-583
/// ```c
/// lookup->word_length = opt->word_size;
/// lookup->lut_word_length = lut_width;
/// lookup->backbone_size = 1 << (BITS_PER_NUC * lookup->lut_word_length);
/// lookup->mask = lookup->backbone_size - 1;
/// lookup->scan_step = lookup->word_length - lookup->lut_word_length + 1;
/// BlastLookupIndexQueryExactMatches(...);
/// s_BlastNaLookupFinalize(thin_backbone, lookup);
/// ```
pub fn build_na_lookup(
    queries_blastna: &[Vec<u8>],
    query_offsets: &[i32],
    word_length: usize,
    lut_word_length: usize,
    query_masks: &[Vec<MaskedInterval>],
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
) -> NaLookupTable {
    debug_assert_eq!(queries_blastna.len(), query_offsets.len());
    debug_assert!(word_length >= lut_word_length);

    let table_size = 1usize << (2 * lut_word_length);
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    if debug_mode {
        let backbone_bytes = table_size * std::mem::size_of::<NaLookupBackboneCell>();
        let pv_bytes = ((table_size >> PV_ARRAY_BTS) + 1) * PV_ARRAY_BYTES;
        eprintln!(
            "[DEBUG] build_na_lookup: word_length={}, lut_word_length={}, table_size={} ({:.1}MB backbone), pv_size={} ({:.1}KB)",
            word_length,
            lut_word_length,
            table_size,
            backbone_bytes as f64 / 1_000_000.0,
            (table_size >> PV_ARRAY_BTS) + 1,
            pv_bytes as f64 / 1_000.0
        );
    }

    let mut counts: Vec<u32> = vec![0; table_size];
    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;
    let mut dust_skipped = 0usize;

    let kmer_mask: u64 = (1u64 << (2 * lut_word_length)) - 1;

    // NCBI reference: blast_lookup.c:BlastLookupIndexQueryExactMatches (lines 79-132)
    // NCBI processes only unmasked regions (locations parameter)
    // Reference: blast_nalookup.c:571-575 (calls BlastLookupIndexQueryExactMatches)
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
        if seq.len() < lut_word_length {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);

        for (range_start, range_end) in ranges {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1001-1022
            // ```c
            // if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
            //     continue;
            // ```
            let range_len = range_end.saturating_sub(range_start);
            if word_length > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases: usize = 0;

            for pos in range_start..range_end {
                let base = seq[pos];
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:117-120
                // ```c
                // if (*seq & invalid_mask)
                //     word_target = seq + lut_word_length + 1;
                // ```
                if (base & BLAST2NA_MASK) != 0 {
                    valid_bases = 0;
                    current_kmer = 0;
                    ambiguous_skipped += 1;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | (base as u64)) & kmer_mask;
                valid_bases += 1;
                if valid_bases < lut_word_length {
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

                let idx = current_kmer as usize;
                counts[idx] = counts[idx].saturating_add(1);
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:456-533
    // ```c
    // lookup->thick_backbone = (NaLookupBackboneCell *)calloc(...);
    // pv = lookup->pv = (PV_ARRAY_TYPE *)calloc(...);
    // if (num_hits > NA_HITS_PER_CELL) overflow_cells_needed += num_hits;
    // ...
    // PV_SET(pv, i, PV_ARRAY_BTS);
    // if (num_hits <= NA_HITS_PER_CELL) { ... entries ... }
    // else { ... overflow ... }
    // ```
    let mut backbone: Vec<NaLookupBackboneCell> = vec![NaLookupBackboneCell::default(); table_size];
    let mut pv: Vec<PvArrayType> = vec![0; (table_size >> PV_ARRAY_BTS) + 1];
    let mut overflow_size = 0usize;
    let mut longest_chain = 0usize;
    let mut non_empty_buckets = 0usize;

    for idx in 0..table_size {
        let num_hits = counts[idx] as usize;
        if num_hits == 0 {
            continue;
        }
        non_empty_buckets += 1;
        longest_chain = longest_chain.max(num_hits);
        backbone[idx].num_used = num_hits as u32;
        pv_set(&mut pv, idx);
        if num_hits > NA_HITS_PER_CELL {
            unsafe {
                // SAFETY: Storing overflow_cursor is valid for chains larger
                // than NA_HITS_PER_CELL, matching NCBI union usage.
                backbone[idx].payload.overflow_cursor = overflow_size as u32;
            }
            overflow_size += num_hits;
        }
    }

    let mut overflow: Vec<u32> = vec![0u32; overflow_size];
    let mut write_pos: Vec<u32> = vec![0; table_size];

    // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
    // Adds ALL hits without any frequency limit - no query-side filtering
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
        if seq.len() < lut_word_length {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);
        let query_offset = query_offsets[q_idx] as usize;

        for (range_start, range_end) in ranges {
            let range_len = range_end.saturating_sub(range_start);
            if word_length > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases: usize = 0;

            for pos in range_start..range_end {
                let base = seq[pos];
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:117-120
                // ```c
                // if (*seq & invalid_mask)
                //     word_target = seq + lut_word_length + 1;
                // ```
                if (base & BLAST2NA_MASK) != 0 {
                    valid_bases = 0;
                    current_kmer = 0;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | (base as u64)) & kmer_mask;
                valid_bases += 1;
                if valid_bases < lut_word_length {
                    continue;
                }

                let kmer_start = pos + 1 - lut_word_length;
                if let Some(counts_filter) = db_word_counts {
                    if db_word_count_exceeds(counts_filter, current_kmer, max_db_word_count) {
                        continue;
                    }
                }

                let idx = current_kmer as usize;
                let pos_idx = write_pos[idx] as usize;
                let q_off_1 = (query_offset + kmer_start + 1) as u32;
                if counts[idx] as usize > NA_HITS_PER_CELL {
                    unsafe {
                        // SAFETY: overflow_cursor is initialized for large chains.
                        let start = backbone[idx].payload.overflow_cursor as usize;
                        overflow[start + pos_idx] = q_off_1;
                    }
                } else {
                    unsafe {
                        // SAFETY: entries is valid for chains with <= NA_HITS_PER_CELL hits.
                        backbone[idx].payload.entries[pos_idx] = q_off_1;
                    }
                }
                write_pos[idx] = write_pos[idx].saturating_add(1);
            }
        }
    }

    if debug_mode {
        let overflow_bytes = overflow.len() * std::mem::size_of::<u32>();
        eprintln!(
            "[DEBUG] build_na_lookup: total_positions={}, ambiguous_skipped={} ({:.1}%), dust_skipped={} ({:.1}%), non_empty_buckets={}, overflow_bytes={:.1}MB, longest_chain={}",
            total_positions,
            ambiguous_skipped,
            100.0 * ambiguous_skipped as f64 / total_positions.max(1) as f64,
            dust_skipped,
            100.0 * dust_skipped as f64 / total_positions.max(1) as f64,
            non_empty_buckets,
            overflow_bytes as f64 / 1_000_000.0,
            longest_chain
        );
    }

    NaLookupTable {
        backbone,
        overflow,
        pv,
        longest_chain,
        word_length,
        lut_word_length,
    }
}

/// Build two-stage lookup table (like NCBI BLAST)
/// - lut_word_length: Used for indexing (typically 8 for megablast)
/// - word_length: Used for extension triggering (typically 28 for megablast)
/// This allows O(1) direct array access even for large word_length values
fn build_mb_lookup(
    queries_blastna: &[Vec<u8>],
    query_offsets: &[i32],
    word_length: usize,
    lut_word_length: usize,
    query_masks: &[Vec<MaskedInterval>],
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    approx_table_entries: usize,
    ascending_cells: bool,
) -> MbLookupTable {
    debug_assert_eq!(queries_blastna.len(), query_offsets.len());

    let table_size = 1usize << (2 * lut_word_length);
    let (pv_size, pv_array_bts) = compute_mb_pv_params(
        table_size,
        approx_table_entries,
        db_word_counts.is_some(),
        lut_word_length,
    );
    let mut hashtable = vec![0u32; table_size];
    let max_query_offset = queries_blastna
        .iter()
        .zip(query_offsets.iter())
        .map(|(seq, &query_offset)| query_offset.max(0) as usize + seq.len())
        .max()
        .unwrap_or(0);
    let mut next_pos = vec![0u32; max_query_offset + 1];
    let mut pv_array = vec![0u32; pv_size];
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();
    let kmer_mask: u64 = (1u64 << (2 * lut_word_length)) - 1;
    const K_COMPRESSION_FACTOR: usize = 2048;
    let mut helper_array =
        vec![0u32; ((table_size + K_COMPRESSION_FACTOR - 1) / K_COMPRESSION_FACTOR).max(1)];
    let mut words: Vec<(u32, u32)> = Vec::new();
    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:979-1108
    // ```c
    // mb_lt->next_pos = (Int4 *)calloc(query->length + 1, sizeof(Int4));
    // ...
    // for (loc = location; loc; loc = loc->next) {
    //     if (full_word_size > (loc->ssr->right - loc->ssr->left + 1))
    //         continue;
    //     ...
    //     if ((val & BLAST2NA_MASK) != 0) {
    //         ecode = 0;
    //         pos = seq + kLutWordLength;
    //         continue;
    //     }
    //     ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
    //     if (seq < pos)
    //         continue;
    //     ...
    //     if (mb_lt->hashtable[ecode] == 0) {
    //         PV_SET(pv_array, ecode, pv_array_bts);
    //     } else {
    //         helper_array[ecode/kCompressionFactor]++;
    //     }
    //     mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    //     mb_lt->hashtable[ecode] = index;
    // }
    // longest_chain = 2;
    // for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
    //     longest_chain = MAX(longest_chain, helper_array[index]);
    // ```
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
        if seq.len() < lut_word_length {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);
        let query_offset = query_offsets[q_idx].max(0) as usize;

        for (range_start, range_end) in ranges {
            let range_len = range_end.saturating_sub(range_start);
            if word_length > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases = 0usize;

            for pos in range_start..range_end {
                let base = seq[pos];
                if (base & BLAST2NA_MASK) != 0 {
                    current_kmer = 0;
                    valid_bases = 0;
                    ambiguous_skipped += 1;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | base as u64) & kmer_mask;
                valid_bases += 1;
                if valid_bases < lut_word_length {
                    continue;
                }

                total_positions += 1;

                if let Some(counts_filter) = db_word_counts {
                    if db_word_count_exceeds(counts_filter, current_kmer, max_db_word_count) {
                        continue;
                    }
                }

                let bucket = current_kmer as usize;
                let q_off_1 = (query_offset + (pos + 1 - lut_word_length) + 1) as u32;
                if hashtable[bucket] == 0 {
                    pv_set_shift(&mut pv_array, bucket, pv_array_bts);
                } else {
                    helper_array[bucket / K_COMPRESSION_FACTOR] =
                        helper_array[bucket / K_COMPRESSION_FACTOR].saturating_add(1);
                }
                next_pos[q_off_1 as usize] = hashtable[bucket];
                hashtable[bucket] = q_off_1;
                if ascending_cells {
                    words.push((bucket as u32, q_off_1));
                }
            }
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:74-76
    // ```c
    // /* add the hit */
    // chain[chain[1] + 2] = query_offset;
    // chain[1]++;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:289-293
    // ```c
    // lookup->final_backbone[i] = -overflow_cursor;
    // for (j = 0; j < num_hits; j++) {
    //     lookup->overflow[overflow_cursor++] =
    //         thin_backbone[i][j + 2];
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:522-526
    // ```c
    // for (j = 0; j < num_hits; j++) {
    //     lookup->overflow[overflow_cursor] =
    //         thin_backbone[i][j + 2];
    //     overflow_cursor++;
    // }
    // ```
    // The small (`eSmallNaLookupTable`) and standard (`eNaLookupTable`) tables
    // append each word to its cell and keep that order, so a cell lists query
    // offsets in the order of indexing; the megablast chain above lists them newest
    // first. Linking the words in reverse order of indexing makes the chain list
    // them in the order of indexing.
    if ascending_cells {
        for &(bucket, _) in &words {
            hashtable[bucket as usize] = 0;
        }
        for &(bucket, q_off_1) in words.iter().rev() {
            next_pos[q_off_1 as usize] = hashtable[bucket as usize];
            hashtable[bucket as usize] = q_off_1;
        }
    }

    let mut longest_chain = 2usize;
    for &value in &helper_array {
        longest_chain = longest_chain.max(value as usize);
    }

    if debug_mode {
        let hashtable_bytes = hashtable.len() * std::mem::size_of::<u32>();
        let next_pos_bytes = next_pos.len() * std::mem::size_of::<u32>();
        eprintln!(
            "[DEBUG] build_mb_lookup: lut_word_length={}, total_positions={}, ambiguous_skipped={}, hashtable={:.1}MB, next_pos={:.1}MB, longest_chain={}",
            lut_word_length,
            total_positions,
            ambiguous_skipped,
            hashtable_bytes as f64 / 1_000_000.0,
            next_pos_bytes as f64 / 1_000_000.0,
            longest_chain,
        );
    }

    MbLookupTable {
        hashtable,
        next_pos,
        pv_array,
        pv_array_bts,
        longest_chain,
    }
}

pub fn build_two_stage_lookup(
    queries_blastna: &[Vec<u8>],
    query_offsets: &[i32],
    word_length: usize,
    lut_word_length: usize,
    query_masks: &[Vec<MaskedInterval>],
    db_word_counts: Option<&[u8]>,
    max_db_word_count: u8,
    approx_table_entries: usize,
    ascending_cells: bool,
) -> TwoStageLookup {
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    if debug_mode {
        eprintln!(
            "[DEBUG] build_two_stage_lookup: word_length={}, lut_word_length={}",
            word_length, lut_word_length
        );
    }

    let mb_lookup = build_mb_lookup(
        queries_blastna,
        query_offsets,
        word_length,
        lut_word_length,
        query_masks,
        db_word_counts,
        max_db_word_count,
        approx_table_entries,
        ascending_cells,
    );

    TwoStageLookup {
        mb_lookup,
        lut_word_length,
        word_length,
        disc: None,
    }
}

/// Builds the discontiguous megablast table: the megablast table of the discontiguous
/// words of the first template (and of the second one, which shares the presence
/// vector), `lut_width` = the word size.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1253-1260
/// ```c
///    ASSERT(lut_width >= 9);
///    mb_lt->word_length = lookup_options->word_size;
/// /*   mb_lt->skip = lookup_options->skip; */
///    mb_lt->stride = lookup_options->stride > 0;
///    mb_lt->lut_word_length = lut_width;
///    mb_lt->hashsize = 1ULL << (BITS_PER_NUC * mb_lt->lut_word_length);
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1327-1331
/// ```c
///    if (lookup_options->mb_template_length > 0) {
///         /* discontiguous megablast */
///         mb_lt->scan_step = 1;
///         status = s_FillDiscMBTable(query, location, mb_lt, lookup_options);
///    }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:645-827
/// ```c
/// static Int2
/// s_FillDiscMBTable(BLAST_SequenceBlk* query, BlastSeqLoc* location,
///         BlastMBLookupTable* mb_lt,
///         const LookupTableOptions* lookup_options)
///
/// {
///    ...
///    const Int4 kCompressionFactor=2048; /* compress helper_array by this much */
///    ...
///    mb_lt->next_pos = (Int4 *)calloc(query->length + 1, sizeof(Int4));
///    ...
///    helper_array = (Uint4*) calloc(mb_lt->hashsize/kCompressionFactor,
///                                   sizeof(Uint4));
///    ...
///    template_type = s_GetDiscTemplateType(lookup_options->word_size,
///                       lookup_options->mb_template_length,
///                       (EDiscWordType)lookup_options->mb_template_type);
///    ...
///    if (kTwoTemplates) {
///       /* Use the temporaray to avoid annoying ICC warning. */
///       int temp_int = template_type + 1;
///       second_template_type =
///            mb_lt->second_template_type = (EDiscTemplateType) temp_int;
///
///       mb_lt->hashtable2 = (Int4*)calloc(mb_lt->hashsize, sizeof(Int4));
///       mb_lt->next_pos2 = (Int4*)calloc(query->length + 1, sizeof(Int4));
///       helper_array2 = (Uint4*) calloc(mb_lt->hashsize/kCompressionFactor,
///                                       sizeof(Uint4));
///       ...
///    }
///
///    mb_lt->discontiguous = TRUE;
///    mb_lt->template_length = lookup_options->mb_template_length;
///    template_length = lookup_options->mb_template_length;
///    pv_array = mb_lt->pv_array;
///    pv_array_bts = mb_lt->pv_array_bts;
///
///    for (loc = location; loc; loc = loc->next) {
///       Int4 from;
///       Int4 to;
///       Uint8 accum = 0;
///       Int4 ecode1 = 0;
///       Int4 ecode2 = 0;
///       Uint1* pos;
///       Uint1* seq;
///       Uint1 val;
///
///       /* A word is added to the table after the last base
///          in the word is read in. At that point, the start
///          offset of the word is (template_length-1) positions
///          behind. This index is also incremented, because
///          lookup table indices are 1-based (offset 0 is reserved). */
///
///       from = loc->ssr->left - (template_length - 2);
///       to = loc->ssr->right - (template_length - 2);
///       seq = query->sequence_start + loc->ssr->left;
///       pos = seq + template_length;
///
///       for (index = from; index <= to; index++) {
///          val = *++seq;
///          /* if an ambiguity is encountered, do not add
///             any words that would contain it */
///          if ((val & BLAST2NA_MASK) != 0) {
///             accum = 0;
///             pos = seq + template_length;
///             continue;
///          }
///
///          /* get next base */
///          accum = (accum << BITS_PER_NUC) | val;
///          if (seq < pos)
///             continue;
///          ...
///          ecode1 = ComputeDiscontiguousIndex(accum, template_type);
///          if (mb_lt->hashtable[ecode1] == 0) {
///             ...
///             PV_SET(pv_array, ecode1, pv_array_bts);
///          }
///          else {
///             helper_array[ecode1/kCompressionFactor]++;
///          }
///          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
///          mb_lt->hashtable[ecode1] = index;
///
///          if (!kTwoTemplates)
///             continue;
///
///          /* repeat for the second template, if applicable */
///
///          ecode2 = ComputeDiscontiguousIndex(accum, second_template_type);
///          if (mb_lt->hashtable2[ecode2] == 0) {
///             ...
///             PV_SET(pv_array, ecode2, pv_array_bts);
///          }
///          else {
///             helper_array2[ecode2/kCompressionFactor]++;
///          }
///          mb_lt->next_pos2[index] = mb_lt->hashtable2[ecode2];
///          mb_lt->hashtable2[ecode2] = index;
///       }
///    }
///
///    longest_chain = 2;
///    for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
///        longest_chain = MAX(longest_chain, helper_array[index]);
///     /* +1 because helper_array is not incremented for the first position of a
///       word */
///    mb_lt->longest_chain = longest_chain + 1;
///    sfree(helper_array);
///
///    if (kTwoTemplates) {
///       longest_chain = 2;
///       for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
///          longest_chain = MAX(longest_chain, helper_array2[index]);
///       /* +1 because helper_array2 is not incremented for the first position of a
///          word */
///       mb_lt->longest_chain += longest_chain + 1;
///       sfree(helper_array2);
///    }
///    return 0;
/// }
/// ```
/// The locations are the unmasked ranges of each context, in the order of
/// `build_mb_lookup`; `left`/`right` are offsets in the concatenated query.
#[allow(clippy::too_many_arguments)]
pub fn build_disc_mb_lookup(
    queries_blastna: &[Vec<u8>],
    query_offsets: &[i32],
    word_size: usize,
    template_length: usize,
    word_type: DiscWordType,
    query_masks: &[Vec<MaskedInterval>],
    approx_table_entries: usize,
) -> TwoStageLookup {
    debug_assert_eq!(queries_blastna.len(), query_offsets.len());
    let lut_word_length = word_size;
    let hashsize = 1usize << (2 * lut_word_length);
    let (pv_size, pv_array_bts) =
        compute_mb_pv_params(hashsize, approx_table_entries, false, lut_word_length);
    let template_type = get_disc_template_type(word_size as i32, template_length as u8, word_type);
    debug_assert!(template_type != DiscTemplateType::Contiguous);
    let two_templates = word_type == DiscWordType::TwoTemplates;
    let second_template_type = if two_templates {
        template_type.next()
    } else {
        DiscTemplateType::Contiguous
    };
    let query_length = queries_blastna
        .iter()
        .zip(query_offsets.iter())
        .map(|(seq, &query_offset)| query_offset.max(0) as usize + seq.len())
        .max()
        .unwrap_or(0);
    const K_COMPRESSION_FACTOR: usize = 2048;
    let mut hashtable = vec![0u32; hashsize];
    let mut next_pos = vec![0u32; query_length + 1];
    let mut helper_array = vec![0u32; hashsize / K_COMPRESSION_FACTOR];
    let (mut hashtable2, mut next_pos2, mut helper_array2) = if two_templates {
        (
            vec![0u32; hashsize],
            vec![0u32; query_length + 1],
            vec![0u32; hashsize / K_COMPRESSION_FACTOR],
        )
    } else {
        (Vec::new(), Vec::new(), Vec::new())
    };
    let mut pv_array = vec![0u32; pv_size];

    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);
        let query_offset = query_offsets[q_idx].max(0) as usize;
        for (range_start, range_end) in build_unmasked_ranges(seq.len(), masks) {
            let mut accum: u64 = 0;
            let mut valid_bases = 0usize;
            for pos in range_start..range_end {
                let val = seq[pos];
                if (val & BLAST2NA_MASK) != 0 {
                    accum = 0;
                    valid_bases = 0;
                    continue;
                }
                accum = (accum << 2) | val as u64;
                valid_bases += 1;
                if valid_bases < template_length {
                    continue;
                }
                // The 1-based start of the word in the concatenated query.
                let index = (query_offset + pos + 2 - template_length) as u32;
                let ecode1 = compute_discontiguous_index(accum, template_type) as usize;
                if hashtable[ecode1] == 0 {
                    pv_set_shift(&mut pv_array, ecode1, pv_array_bts);
                } else {
                    helper_array[ecode1 / K_COMPRESSION_FACTOR] =
                        helper_array[ecode1 / K_COMPRESSION_FACTOR].wrapping_add(1);
                }
                next_pos[index as usize] = hashtable[ecode1];
                hashtable[ecode1] = index;
                if !two_templates {
                    continue;
                }
                let ecode2 = compute_discontiguous_index(accum, second_template_type) as usize;
                if hashtable2[ecode2] == 0 {
                    pv_set_shift(&mut pv_array, ecode2, pv_array_bts);
                } else {
                    helper_array2[ecode2 / K_COMPRESSION_FACTOR] =
                        helper_array2[ecode2 / K_COMPRESSION_FACTOR].wrapping_add(1);
                }
                next_pos2[index as usize] = hashtable2[ecode2];
                hashtable2[ecode2] = index;
            }
        }
    }

    let mut longest_chain = 2u32;
    for &value in &helper_array {
        longest_chain = longest_chain.max(value);
    }
    let mut total_longest_chain = longest_chain as usize + 1;
    if two_templates {
        let mut longest_chain2 = 2u32;
        for &value in &helper_array2 {
            longest_chain2 = longest_chain2.max(value);
        }
        total_longest_chain += longest_chain2 as usize + 1;
    }

    TwoStageLookup {
        mb_lookup: MbLookupTable {
            hashtable,
            next_pos,
            pv_array,
            pv_array_bts,
            longest_chain: total_longest_chain,
        },
        lut_word_length,
        word_length: word_size,
        disc: Some(DiscTemplates {
            template_type,
            second_template_type,
            two_templates,
            template_length,
            hashtable2,
            next_pos2,
        }),
    }
}

/// Build direct address table for k-mer lookup (O(1) access)
/// This is much faster than HashMap for small word sizes
pub fn build_direct_lookup(
    queries_blastna: &[Vec<u8>],
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

    debug_assert_eq!(queries_blastna.len(), query_offsets.len());

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
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
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
    for (q_idx, seq_blastna) in queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
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
        let non_empty = offsets.windows(2).filter(|w| w[0] != w[1]).count();
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

#[cfg(test)]
mod reverse_complement_tests {
    use super::reverse_complement;

    // A lowercase (soft-masked) letter is the same base on the minus strand.
    #[test]
    fn lowercase_letters_are_complemented_like_uppercase() {
        assert_eq!(reverse_complement(b"ACgtRn-"), b"-NYACGT");
    }
}

#[cfg(test)]
mod chain_order_tests {
    use super::build_two_stage_lookup;

    fn hits(ascending_cells: bool) -> Vec<u32> {
        // "ACGTACGT" repeated: the 4-letter word ACGT occurs at offsets 0, 4, 8, ...
        // of the first context and of the second (offset 21 in the block).
        let unit = [0u8, 1, 2, 3];
        let context: Vec<u8> = unit.iter().cycle().take(20).copied().collect();
        let queries = vec![context.clone(), context];
        let lookup = build_two_stage_lookup(
            &queries,
            &[0, 21],
            8,
            4,
            &[Vec::new(), Vec::new()],
            None,
            0,
            40,
            ascending_cells,
        );
        let acgt = 0b00_01_10_11u64;
        let mut found = Vec::new();
        lookup.for_each_hit(acgt, |q_off_1| found.push(q_off_1 - 1));
        found
    }

    #[test]
    fn small_and_standard_cells_list_offsets_in_indexing_order() {
        let ascending = hits(true);
        assert_eq!(ascending, vec![0, 4, 8, 12, 16, 21, 25, 29, 33, 37]);
        let mut newest_first = hits(false);
        assert_eq!(newest_first.first(), Some(&37));
        newest_first.reverse();
        assert_eq!(newest_first, ascending);
    }
}

#[cfg(test)]
mod disc_lookup_tests {
    use super::build_disc_mb_lookup;
    use crate::algorithm::blastn::disc_lookup::{
        compute_discontiguous_index, DiscTemplateType, DiscWordType,
    };
    use crate::utils::dust::MaskedInterval;

    fn index_at(seq: &[u8], start: usize, len: usize, t: DiscTemplateType) -> u64 {
        let mut accum = 0u64;
        for &b in &seq[start..start + len] {
            accum = (accum << 2) | b as u64;
        }
        compute_discontiguous_index(accum, t) as u64
    }

    fn chain(lookup: &super::TwoStageLookup, index: u64, second: bool) -> Vec<u32> {
        let mut hits = Vec::new();
        if second {
            lookup.for_each_hit2(index, |q| hits.push(q));
        } else {
            lookup.for_each_hit(index, |q| hits.push(q));
        }
        hits
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:741-813
    // ```c
    //       from = loc->ssr->left - (template_length - 2);
    //       to = loc->ssr->right - (template_length - 2);
    //       seq = query->sequence_start + loc->ssr->left;
    //       pos = seq + template_length;
    //
    //       for (index = from; index <= to; index++) {
    //          val = *++seq;
    //          /* if an ambiguity is encountered, do not add
    //             any words that would contain it */
    //          if ((val & BLAST2NA_MASK) != 0) {
    //             accum = 0;
    //             pos = seq + template_length;
    //             continue;
    //          }
    // ...
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
    //          mb_lt->hashtable[ecode1] = index;
    // ```
    #[test]
    fn words_are_indexed_newest_first_with_one_based_starts() {
        // Two contexts in the concatenated query (offsets 0 and 41); the second repeats
        // the first's 18-mer at 0 and has an ambiguity code (14) inside its copy at 20.
        let word: Vec<u8> = (0..18).map(|i| (i * 7 % 4) as u8).collect();
        let mut first = word.clone();
        first.extend([1, 2, 3]);
        let mut second = word.clone();
        second.extend([0, 0]);
        second.extend(&word);
        second[20 + 5] = 14;
        let lookup = build_disc_mb_lookup(
            &[first.clone(), second.clone()],
            &[0, 41],
            11,
            18,
            DiscWordType::TwoTemplates,
            &[Vec::new(), Vec::new()],
            1000,
        );
        let disc = lookup.disc().expect("discontiguous table");
        assert_eq!(disc.template_type, DiscTemplateType::T11_18Coding);
        assert_eq!(disc.second_template_type, DiscTemplateType::T11_18Optimal);
        assert_eq!((lookup.word_length(), lookup.lut_word_length()), (18, 18));
        let coding = index_at(&word, 0, 18, DiscTemplateType::T11_18Coding);
        let optimal = index_at(&word, 0, 18, DiscTemplateType::T11_18Optimal);
        // The word at 0 of each context; the copy at 20 holds the ambiguity code.
        assert_eq!(chain(&lookup, coding, false), vec![42, 1]);
        assert_eq!(chain(&lookup, optimal, true), vec![42, 1]);
        // A masked range is not indexed.
        let masked = build_disc_mb_lookup(
            &[first, second],
            &[0, 41],
            11,
            18,
            DiscWordType::Coding,
            &[vec![MaskedInterval::new(3, 4)], Vec::new()],
            1000,
        );
        assert_eq!(chain(&masked, coding, false), vec![42]);
        assert!(masked.disc().is_some_and(|d| !d.two_templates));
        assert!(chain(&masked, optimal, true).is_empty());
    }
}
