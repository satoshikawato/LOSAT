use bio::io::fasta;
use rustc_hash::FxHashMap;
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

/// Generate the reverse complement of a DNA sequence.
/// Used for searching the minus strand of subject sequences.
/// Complement: A<->T, C<->G, ambiguous bases become N
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' | b'U' | b'u' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            _ => b'N',
        })
        .collect()
}

pub type KmerLookup = FxHashMap<u64, Vec<(u32, u32)>>;

/// Direct address table for k-mer lookup - O(1) access instead of hash lookup
/// Used for small word sizes (<=13) where 4^word_size fits in memory
/// For word_size=11: 4^11 = 4,194,304 entries (~100MB)
/// For word_size=13: 4^13 = 67,108,864 entries (~1.6GB)
pub type DirectKmerLookup = Vec<Vec<(u32, u32)>>;

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

/// Optimized lookup table with Presence-Vector for fast filtering
/// Combines DirectKmerLookup with a bit vector for O(1) presence checking
pub struct PvDirectLookup {
    /// The actual lookup table storing (query_idx, position) pairs
    lookup: DirectKmerLookup,
    /// Presence vector - bit i is set if lookup[i] is non-empty
    pv: Vec<PvArrayType>,
    /// Word size used for this lookup table
    #[allow(dead_code)]
    word_size: usize,
}

impl PvDirectLookup {
    /// Check if a k-mer has any hits using the presence vector (O(1))
    #[inline(always)]
    pub fn has_hits(&self, kmer: u64) -> bool {
        pv_test(&self.pv, kmer as usize)
    }

    /// Get hits for a k-mer (only call after has_hits returns true)
    #[inline(always)]
    pub fn get_hits(&self, kmer: u64) -> &[(u32, u32)] {
        let idx = kmer as usize;
        if idx < self.lookup.len() {
            &self.lookup[idx]
        } else {
            &[]
        }
    }

    /// Get hits for a k-mer with PV check (combined operation)
    #[inline(always)]
    pub fn get_hits_checked(&self, kmer: u64) -> &[(u32, u32)] {
        if self.has_hits(kmer) {
            self.get_hits(kmer)
        } else {
            &[]
        }
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
    pub fn get_hits(&self, lut_kmer: u64) -> &[(u32, u32)] {
        self.pv_lookup.get_hits(lut_kmer)
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
}

/// Pack (q_idx, diag) into a single u64 key for faster HashMap operations
#[inline(always)]
pub fn pack_diag_key(q_idx: u32, diag: isize) -> u64 {
    ((q_idx as u64) << 32) | ((diag as i32) as u32 as u64)
}

/// Check if a k-mer starting at position overlaps with any masked interval
#[inline]
pub fn is_kmer_masked(intervals: &[MaskedInterval], start: usize, kmer_len: usize) -> bool {
    let end = start + kmer_len;
    intervals.iter().any(|interval| {
        start < interval.end && end > interval.start
    })
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
    word_size: usize,
    query_masks: &[Vec<MaskedInterval>],
) -> PvDirectLookup {
    let safe_word_size = word_size.min(MAX_DIRECT_LOOKUP_WORD_SIZE);
    let table_size = 1usize << (2 * safe_word_size); // 4^word_size
    let pv_size = (table_size + 31) / 32; // Number of u32 elements needed for PV
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    if debug_mode {
        eprintln!(
            "[DEBUG] build_pv_direct_lookup: word_size={}, table_size={} ({:.1}MB), pv_size={} ({:.1}KB)",
            safe_word_size,
            table_size,
            (table_size * 24) as f64 / 1_000_000.0,
            pv_size,
            (pv_size * 4) as f64 / 1_000.0
        );
    }

    // Pre-allocate the table and presence vector
    let mut lookup: Vec<Vec<(u32, u32)>> = vec![Vec::new(); table_size];
    let mut pv: Vec<PvArrayType> = vec![0; pv_size];

    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;
    let mut dust_skipped = 0usize;
    let mut kmers_added = 0usize;

    // K-mer mask for rolling window
    let kmer_mask: u64 = (1u64 << (2 * safe_word_size)) - 1;

    // NCBI reference: blast_lookup.c:BlastLookupIndexQueryExactMatches (lines 79-132)
    // NCBI processes only unmasked regions (locations parameter)
    // For each location, it iterates through positions and adds k-mers
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
        // Rolling window approach: word_target points to position where complete k-mer can be formed
        // Ambiguous bases reset the window: if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
        // Rolling k-mer state
        let mut current_kmer: u64 = 0;
        let mut valid_bases: usize = 0;

        for pos in 0..seq.len() {
            let base = seq[pos];
            let code = ENCODE_LUT[base as usize];

            // NCBI reference: blast_lookup.c:119-120
            // if (*seq & invalid_mask) word_target = seq + lut_word_length + 1;
            if code == 0xFF {
                // Ambiguous base - reset the rolling window
                valid_bases = 0;
                current_kmer = 0;
                ambiguous_skipped += 1;
                continue;
            }

            // Shift in the new base
            current_kmer = ((current_kmer << 2) | (code as u64)) & kmer_mask;
            valid_bases += 1;

            // NCBI reference: blast_lookup.c:110-115
            // if (seq >= word_target) { BlastLookupAddWordHit(...); }
            // Only process if we have a complete k-mer
            if valid_bases < safe_word_size {
                continue;
            }

            // Calculate the starting position of this k-mer
            let kmer_start = pos + 1 - safe_word_size;
            total_positions += 1;

            // NCBI reference: blast_nalookup.c:402-406
            // BlastLookupIndexQueryExactMatches processes only unmasked regions (locations)
            // LOSAT processes entire sequence and filters masked regions (equivalent)
            // Skip k-mers that overlap with DUST-masked regions
            if !masks.is_empty() && is_kmer_masked(masks, kmer_start, safe_word_size) {
                dust_skipped += 1;
                continue;
            }

            // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
            // Adds ALL hits without any frequency limit - no query-side filtering
            // if (backbone[index] == NULL) { initialize new chain }
            // else { use existing chain, realloc if full }
            // chain[chain[1] + 2] = query_offset; chain[1]++;
            // Add to lookup table
            let idx = current_kmer as usize;
            if idx < table_size {
                lookup[idx].push((q_idx as u32, kmer_start as u32));
                kmers_added += 1;
            }
        }
    }

    // NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
    // NCBI BLAST does NOT filter over-represented k-mers in the query
    // All k-mers are added to the lookup table regardless of frequency
    // Database word count filtering (kDbFilter) exists but is different:
    // it filters based on database counts, not query counts
    // Reference: blast_nalookup.c:1047-1060 (database word count filtering)
    // 
    // Build presence vector for all non-empty k-mers
    let mut non_empty_count = 0usize;

    for (idx, positions) in lookup.iter().enumerate() {
        if !positions.is_empty() {
            // Set bit in presence vector
            pv_set(&mut pv, idx);
            non_empty_count += 1;
        }
    }

    if debug_mode {
        eprintln!(
            "[DEBUG] build_pv_direct_lookup: total_positions={}, ambiguous_skipped={} ({:.1}%), dust_skipped={} ({:.1}%)",
            total_positions,
            ambiguous_skipped,
            100.0 * ambiguous_skipped as f64 / (total_positions + ambiguous_skipped).max(1) as f64,
            dust_skipped,
            100.0 * dust_skipped as f64 / total_positions.max(1) as f64
        );
        eprintln!(
            "[DEBUG] build_pv_direct_lookup: kmers_added={}, non_empty_buckets={}",
            kmers_added, non_empty_count
        );
    }

    PvDirectLookup {
        lookup,
        pv,
        word_size: safe_word_size,
    }
}

pub fn build_lookup(
    queries: &[fasta::Record],
    word_size: usize,
    query_masks: &[Vec<MaskedInterval>],
) -> KmerLookup {
    let mut lookup: FxHashMap<u64, Vec<(u32, u32)>> = FxHashMap::default();
    let safe_word_size = word_size.min(31);
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;
    let mut dust_skipped = 0usize;

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
                lookup
                    .entry(kmer)
                    .or_default()
                    .push((q_idx as u32, i as u32));
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
    word_length: usize,
    lut_word_length: usize,
    query_masks: &[Vec<MaskedInterval>],
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
    let pv_lookup = build_pv_direct_lookup(queries, lut_word_length, query_masks);
    
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
    word_size: usize,
    query_masks: &[Vec<MaskedInterval>],
) -> DirectKmerLookup {
    let safe_word_size = word_size.min(MAX_DIRECT_LOOKUP_WORD_SIZE);
    let table_size = 1usize << (2 * safe_word_size); // 4^word_size
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();

    if debug_mode {
        eprintln!(
            "[DEBUG] build_direct_lookup: word_size={}, table_size={} ({:.1}MB base)",
            safe_word_size,
            table_size,
            (table_size * 24) as f64 / 1_000_000.0
        );
    }

    // Pre-allocate the table with empty vectors
    let mut lookup: Vec<Vec<(u32, u32)>> = vec![Vec::new(); table_size];

    let mut total_positions = 0usize;
    let mut ambiguous_skipped = 0usize;
    let mut dust_skipped = 0usize;

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
                    lookup[idx].push((q_idx as u32, i as u32));
                }
            } else {
                ambiguous_skipped += 1;
            }
        }
    }

    if debug_mode {
        let non_empty = lookup.iter().filter(|v| !v.is_empty()).count();
        eprintln!(
            "[DEBUG] build_direct_lookup: total_positions={}, ambiguous_skipped={} ({:.1}%), dust_skipped={} ({:.1}%), non_empty_buckets={}",
            total_positions,
            ambiguous_skipped,
            100.0 * ambiguous_skipped as f64 / total_positions.max(1) as f64,
            dust_skipped,
            100.0 * dust_skipped as f64 / total_positions.max(1) as f64,
            non_empty
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


