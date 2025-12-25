use anyhow::{Context, Result};
use bio::io::fasta;
use clap::Args;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::path::PathBuf;
use std::sync::mpsc::channel;

use crate::common::{write_output, Hit};
use crate::config::NuclScoringSpec;
use crate::stats::karlin::{bit_score as calc_bit_score, evalue as calc_evalue};
use crate::stats::search_space::SearchSpace;
use crate::stats::{lookup_nucl_params, KarlinParams};
use crate::utils::dust::{DustMasker, MaskedInterval};

/// Thread-local memory pool for non-affine greedy alignment.
/// This avoids per-call allocation overhead by reusing memory across calls.
/// Similar to NCBI BLAST's SGreedyAlignMem structure.
struct GreedyAlignMem {
    /// Two rows for last_seq2_off (we swap between them)
    last_seq2_off_a: Vec<i32>,
    last_seq2_off_b: Vec<i32>,
    /// Array of maximum scores at each distance
    max_score: Vec<i32>,
    /// Current allocated size
    allocated_size: usize,
}

impl GreedyAlignMem {
    fn new() -> Self {
        Self {
            last_seq2_off_a: Vec::new(),
            last_seq2_off_b: Vec::new(),
            max_score: Vec::new(),
            allocated_size: 0,
        }
    }

    /// Ensure the memory pool has enough capacity for the given array_size and max_score_size
    fn ensure_capacity(&mut self, array_size: usize, max_score_size: usize) {
        if self.allocated_size < array_size {
            self.last_seq2_off_a.resize(array_size, -2);
            self.last_seq2_off_b.resize(array_size, -2);
            self.allocated_size = array_size;
        }
        if self.max_score.len() < max_score_size {
            self.max_score.resize(max_score_size, 0);
        }
    }

    /// Reset arrays to initial state (fill with sentinel values)
    fn reset(&mut self, array_size: usize, max_score_size: usize) {
        // Fill with sentinel values
        for i in 0..array_size.min(self.last_seq2_off_a.len()) {
            self.last_seq2_off_a[i] = -2;
            self.last_seq2_off_b[i] = -2;
        }
        for i in 0..max_score_size.min(self.max_score.len()) {
            self.max_score[i] = 0;
        }
    }
}

thread_local! {
    /// Thread-local memory pool for non-affine greedy alignment
    static GREEDY_MEM: RefCell<GreedyAlignMem> = RefCell::new(GreedyAlignMem::new());
}

// NCBI BLAST compatible X-drop parameters
const X_DROP_UNGAPPED: i32 = 20; // BLAST_UNGAPPED_X_DROPOFF_NUCL
const X_DROP_GAPPED_FINAL: i32 = 100; // BLAST_GAP_X_DROPOFF_FINAL_NUCL for final traceback
const TWO_HIT_WINDOW: usize = 64; // Increased from 40 for better sensitivity
const MAX_HITS_PER_KMER: usize = 200;

// Minimum ungapped score thresholds for triggering gapped extension
// Higher threshold for blastn (smaller word size = more seeds) to reduce extension count
// NCBI BLAST uses similar task-specific thresholds
const MIN_UNGAPPED_SCORE_MEGABLAST: i32 = 20; // Lower threshold for megablast (larger seeds)
const MIN_UNGAPPED_SCORE_BLASTN: i32 = 40; // Higher threshold for blastn (smaller seeds, more filtering needed)

#[derive(Args, Debug)]
pub struct BlastnArgs {
    #[arg(short, long)]
    pub query: PathBuf,
    #[arg(short, long)]
    pub subject: PathBuf,
    #[arg(long, default_value = "megablast")]
    pub task: String,
    #[arg(short, long, default_value_t = 28)]
    pub word_size: usize,
    #[arg(short = 'n', long, default_value_t = 0)]
    pub num_threads: usize,
    #[arg(long, default_value_t = 10.0)]
    pub evalue: f64,
    #[arg(long, default_value_t = 500)]
    pub max_target_seqs: usize,
    #[arg(short, long)]
    pub out: Option<PathBuf>,
    // Scoring parameters - defaults are for megablast task
    // For blastn task, these are overridden in run() based on --task
    #[arg(long, default_value_t = 1)]
    pub reward: i32,
    #[arg(long, default_value_t = -2)]
    pub penalty: i32,
    #[arg(long, default_value_t = 0)]
    pub gap_open: i32,
    #[arg(long, default_value_t = 0)]
    pub gap_extend: i32,
    // DUST filter options for masking low-complexity regions
    #[arg(long, default_value_t = true)]
    pub dust: bool,
    #[arg(long, default_value_t = 20)]
    pub dust_level: u32,
    #[arg(long, default_value_t = 64)]
    pub dust_window: usize,
    #[arg(long, default_value_t = 1)]
    pub dust_linker: usize,
    #[arg(long, short = 'v', default_value_t = false)]
    pub verbose: bool,
}

// 2-bit encoding for compact storage and hashing
fn encode_kmer(seq: &[u8], start: usize, k: usize) -> Option<u64> {
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

type KmerLookup = FxHashMap<u64, Vec<(u32, u32)>>;

/// Direct address table for k-mer lookup - O(1) access instead of hash lookup
/// Used for small word sizes (<=13) where 4^word_size fits in memory
/// For word_size=11: 4^11 = 4,194,304 entries (~100MB)
/// For word_size=13: 4^13 = 67,108,864 entries (~1.6GB)
type DirectKmerLookup = Vec<Vec<(u32, u32)>>;

/// Maximum word size for direct address table (4^14 = 268M entries = ~6GB, too large)
const MAX_DIRECT_LOOKUP_WORD_SIZE: usize = 13;

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

// ============================================================================
// Phase 4: Key packing optimization for scanning HashMaps
// ============================================================================
// The scanning loop uses HashMaps with (q_idx, diag) tuple keys for `mask` and `last_seed`.
// Tuple keys have overhead: 12 bytes for (u32, isize) + hashing overhead.
// By packing into a single u64, we reduce key size and improve hashing performance.
//
// Packing scheme:
// - High 32 bits: q_idx (query index)
// - Low 32 bits: diag (diagonal, cast to u32 with wrapping)
// This allows O(1) pack/unpack and faster hashing.

/// Pack (q_idx, diag) into a single u64 key for faster HashMap operations
#[inline(always)]
fn pack_diag_key(q_idx: u32, diag: isize) -> u64 {
    ((q_idx as u64) << 32) | ((diag as i32) as u32 as u64)
}

/// Check if a k-mer starting at position overlaps with any masked interval
#[inline]
fn is_kmer_masked(intervals: &[MaskedInterval], start: usize, kmer_len: usize) -> bool {
    let end = start + kmer_len;
    intervals.iter().any(|interval| {
        start < interval.end && end > interval.start
    })
}

// ============================================================================
// Phase 2: Query Packing with Rolling K-mer Extraction
// ============================================================================
// Uses O(1) sliding window k-mer extraction for building lookup tables
// This is more efficient than calling encode_kmer() for each position

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
fn build_pv_direct_lookup(
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

    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);

        // Rolling k-mer state
        let mut current_kmer: u64 = 0;
        let mut valid_bases: usize = 0;

        for pos in 0..seq.len() {
            let base = seq[pos];
            let code = ENCODE_LUT[base as usize];

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

            // Only process if we have a complete k-mer
            if valid_bases < safe_word_size {
                continue;
            }

            // Calculate the starting position of this k-mer
            let kmer_start = pos + 1 - safe_word_size;
            total_positions += 1;

            // Skip k-mers that overlap with DUST-masked regions
            if !masks.is_empty() && is_kmer_masked(masks, kmer_start, safe_word_size) {
                dust_skipped += 1;
                continue;
            }

            // Add to lookup table
            let idx = current_kmer as usize;
            if idx < table_size {
                lookup[idx].push((q_idx as u32, kmer_start as u32));
                kmers_added += 1;
            }
        }
    }

    // Filter over-represented k-mers and build presence vector
    let mut filtered_count = 0usize;
    let mut non_empty_count = 0usize;

    for (idx, positions) in lookup.iter_mut().enumerate() {
        if positions.len() > MAX_HITS_PER_KMER {
            positions.clear();
            filtered_count += 1;
        } else if !positions.is_empty() {
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
            "[DEBUG] build_pv_direct_lookup: kmers_added={}, filtered={}, non_empty_buckets={}",
            kmers_added, filtered_count, non_empty_count
        );
    }

    PvDirectLookup {
        lookup,
        pv,
        word_size: safe_word_size,
    }
}

fn build_lookup(
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

    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);

        for i in 0..=(seq.len() - safe_word_size) {
            total_positions += 1;
            
            // Skip k-mers that overlap with DUST-masked regions
            if !masks.is_empty() && is_kmer_masked(masks, i, safe_word_size) {
                dust_skipped += 1;
                continue;
            }
            
            if let Some(kmer) = encode_kmer(seq, i, safe_word_size) {
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

    // Filter over-represented k-mers to prevent seed explosion
    // This is critical for performance with smaller word sizes (e.g., blastn task with word_size=11)
    let before_filter = lookup.len();
    lookup.retain(|_, positions| positions.len() <= MAX_HITS_PER_KMER);
    let after_filter = lookup.len();

    if debug_mode {
        eprintln!(
            "[DEBUG] build_lookup: filtered {} high-frequency kmers (before={}, after={})",
            before_filter - after_filter,
            before_filter,
            after_filter
        );
    }

    lookup
}

/// Build direct address table for k-mer lookup (O(1) access)
/// This is much faster than HashMap for small word sizes
fn build_direct_lookup(
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

    for (q_idx, record) in queries.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < safe_word_size {
            continue;
        }

        let masks = query_masks.get(q_idx).map(|v| v.as_slice()).unwrap_or(&[]);

        for i in 0..=(seq.len() - safe_word_size) {
            total_positions += 1;
            
            // Skip k-mers that overlap with DUST-masked regions
            if !masks.is_empty() && is_kmer_masked(masks, i, safe_word_size) {
                dust_skipped += 1;
                continue;
            }
            
            if let Some(kmer) = encode_kmer(seq, i, safe_word_size) {
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

    // Filter over-represented k-mers to prevent seed explosion
    let mut filtered_count = 0usize;
    for positions in lookup.iter_mut() {
        if positions.len() > MAX_HITS_PER_KMER {
            positions.clear();
            filtered_count += 1;
        }
    }

    if debug_mode {
        eprintln!(
            "[DEBUG] build_direct_lookup: filtered {} high-frequency kmers",
            filtered_count
        );
    }

    lookup
}

/// Find the first mismatch between two sequences starting at given offsets.
/// This is a key optimization from NCBI BLAST - quickly scan for exact matches
/// before doing any DP computation.
///
/// Returns the number of consecutive matches found.
/// Find the first mismatch between two sequences starting at given offsets.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// This avoids the need to copy and reverse sequences for left extension.
#[inline(always)]
fn find_first_mismatch_ex(
    seq1: &[u8],
    seq2: &[u8],
    len1: usize,
    len2: usize,
    start1: usize,
    start2: usize,
    reverse: bool,
) -> usize {
    let max_len = (len1 - start1).min(len2 - start2);
    let mut count = 0;

    if reverse {
        // Access sequences from the end (NCBI BLAST approach for left extension)
        while count < max_len {
            let idx1 = len1 - 1 - start1 - count;
            let idx2 = len2 - 1 - start2 - count;
            if unsafe { *seq1.get_unchecked(idx1) != *seq2.get_unchecked(idx2) } {
                break;
            }
            count += 1;
        }
    } else {
        // Forward access (normal right extension)
        let s1 = &seq1[start1..];
        let s2 = &seq2[start2..];
        while count < max_len {
            if unsafe { *s1.get_unchecked(count) != *s2.get_unchecked(count) } {
                break;
            }
            count += 1;
        }
    }

    count
}

#[inline(always)]
fn find_first_mismatch(seq1: &[u8], seq2: &[u8], start1: usize, start2: usize) -> usize {
    find_first_mismatch_ex(seq1, seq2, seq1.len(), seq2.len(), start1, start2, false)
}

// NCBI BLAST constants for greedy alignment
// GREEDY_MAX_COST: The largest distance to be examined for an optimal alignment
const GREEDY_MAX_COST: usize = 1000;
// GREEDY_MAX_COST_FRACTION: sequence_length / this is a measure of how hard the algorithm will work
const GREEDY_MAX_COST_FRACTION: usize = 2;

/// Bookkeeping structure for affine greedy alignment (NCBI BLAST's SGreedyOffset).
/// When aligning two sequences, stores the largest offset into the second sequence
/// that leads to a high-scoring alignment for a given start point, tracking
/// different path endings separately for affine gap penalties.
#[derive(Clone, Copy, Default)]
struct GreedyOffset {
    insert_off: i32,  // Best offset for a path ending in an insertion (gap in seq2)
    match_off: i32,   // Best offset for a path ending in a match/mismatch
    delete_off: i32,  // Best offset for a path ending in a deletion (gap in seq1)
}

/// Signal that a diagonal/offset is invalid
const INVALID_OFFSET: i32 = -2;
/// Large value for invalid diagonal bounds
const INVALID_DIAG: i32 = 100_000_000;

/// GCD of three values (NCBI BLAST's BLAST_Gdb3)
/// Divides all three values by their GCD and returns the GCD
fn gdb3(a: &mut i32, b: &mut i32, c: &mut i32) -> i32 {
    let g = if *b == 0 {
        gcd_i32(*a, *c)
    } else {
        gcd_i32(*a, gcd_i32(*b, *c))
    };
    if g > 1 {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    g
}

/// GCD for i32 values
fn gcd_i32(a: i32, b: i32) -> i32 {
    let mut a = a.abs();
    let mut b = b.abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

/// Greedy alignment for high-identity nucleotide sequences.
///
/// This implements NCBI BLAST's greedy alignment algorithm (Zhang et al., 2000):
/// - Uses distance-based tracking instead of score-based
/// - Quickly scans for exact matches before doing any DP
/// - Only explores diagonals that can achieve a given distance
/// - Much faster than full DP for similar sequences
///
/// For sequences with >90% identity, this is typically 10-100x faster than full DP.
///
/// This is a wrapper that implements NCBI BLAST's dynamic max_dist doubling strategy:
/// - Start with initial max_dist based on sequence length
/// - If alignment doesn't converge, double max_dist and retry
/// - This allows finding arbitrarily long alignments
///
/// For affine gap penalties (gap_open != 0 || gap_extend != 0), uses the affine
/// greedy algorithm (NCBI BLAST's BLAST_AffineGreedyAlign).
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters)
/// Greedy alignment for high-identity nucleotide sequences.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// This avoids the need to copy and reverse sequences for left extension.
fn greedy_align_one_direction_ex(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    reverse: bool,
) -> (usize, usize, i32, usize, usize, usize, usize) {
    // Calculate initial max_dist like NCBI BLAST:
    // max_dist = MIN(GREEDY_MAX_COST, max(len1, len2) / GREEDY_MAX_COST_FRACTION + 1)
    let max_len = len1.max(len2);
    let mut max_dist = GREEDY_MAX_COST.min(max_len / GREEDY_MAX_COST_FRACTION + 1);
    
    // Retry with doubled max_dist until convergence (NCBI BLAST approach)
    loop {
        // Choose between affine and non-affine greedy based on gap penalties
        // This matches NCBI BLAST's BLAST_AffineGreedyAlign which falls back to
        // BLAST_GreedyAlign when gap_open == 0 && gap_extend == 0
        let (result, converged) = if gap_open != 0 || gap_extend != 0 {
            affine_greedy_align_one_direction_with_max_dist(
                q_seq, s_seq, reward, penalty, gap_open, gap_extend, x_drop, max_dist,
            )
        } else {
            greedy_align_one_direction_with_max_dist(
                q_seq, s_seq, len1, len2, reward, penalty, gap_open, gap_extend, x_drop, max_dist, reverse,
            )
        };
        
        if converged {
            return result;
        }
        
        // Double max_dist and retry (NCBI BLAST's approach)
        // NCBI BLAST continues until convergence without an explicit upper limit
        max_dist *= 2;
        
        // Safety limit to prevent infinite loops on extremely divergent sequences
        // This is a practical limit that should never be reached for reasonable sequences
        // NCBI BLAST typically converges well before this limit
        if max_dist > 100000 {
            // Return best result found so far
            return result;
        }
    }
}

/// Wrapper for backward compatibility - uses slice lengths and forward direction
fn greedy_align_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, i32, usize, usize, usize, usize) {
    greedy_align_one_direction_ex(
        q_seq, s_seq, q_seq.len(), s_seq.len(),
        reward, penalty, gap_open, gap_extend, x_drop, false,
    )
}

/// Internal greedy alignment function with explicit max_dist parameter.
/// If `reverse` is true, sequences are accessed from the end (for left extension).
/// Returns: ((q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters), converged)
fn greedy_align_one_direction_with_max_dist(
    q_seq: &[u8],
    s_seq: &[u8],
    len1: usize,
    len2: usize,
    reward: i32,
    penalty: i32,
    _gap_open: i32,
    _gap_extend: i32,
    x_drop: i32,
    max_dist: usize,
    reverse: bool,
) -> ((usize, usize, i32, usize, usize, usize, usize), bool) {
    if len1 == 0 || len2 == 0 {
        return ((0, 0, 0, 0, 0, 0, 0), true);
    }

    // Calculate match and mismatch costs for distance-based tracking
    // NCBI BLAST doubles scores if reward is odd to avoid fractions
    // IMPORTANT: x_drop must also be doubled to maintain consistent units
    let (match_cost, mismatch_cost, scaled_xdrop) = if reward % 2 == 1 {
        (reward * 2, (-penalty) * 2, x_drop * 2)
    } else {
        (reward, -penalty, x_drop)
    };

    let op_cost = match_cost + mismatch_cost; // Cost of a mismatch in distance terms

    // X-drop offset for score comparison (using scaled_xdrop for consistent units)
    let xdrop_offset = ((scaled_xdrop + match_cost / 2) / op_cost.max(1) + 1) as usize;

    // Find initial run of matches
    let initial_matches = find_first_mismatch_ex(q_seq, s_seq, len1, len2, 0, 0, reverse);

    if initial_matches == len1 || initial_matches == len2 {
        // Perfect match - return immediately (converged)
        let score = (initial_matches as i32) * reward;
        return (
            (
                initial_matches,
                initial_matches,
                score,
                initial_matches,
                0,
                0,
                0,
            ),
            true,
        );
    }

    // Diagonal origin (center of diagonal space)
    // NCBI BLAST uses max_dist (not scaled_max_dist) for diag_origin
    let diag_origin = max_dist + 2;
    let array_size = 2 * diag_origin + 4;
    let max_score_size = max_dist + xdrop_offset + 2;

    // Use thread-local memory pool to avoid per-call allocation overhead
    // This is critical for performance - without pooling, non-affine greedy is 20x slower
    GREEDY_MEM.with(|mem_cell| {
        let mut mem = mem_cell.borrow_mut();
        mem.ensure_capacity(array_size, max_score_size);
        mem.reset(array_size, max_score_size);

        // Initialize distance 0
        mem.last_seq2_off_a[diag_origin] = initial_matches as i32;
        mem.max_score[xdrop_offset] = (initial_matches as i32) * match_cost;

        let mut best_dist = 0usize;
        let mut best_seq1_len = initial_matches;
        let mut best_seq2_len = initial_matches;

        let mut diag_lower = diag_origin as i32 - 1;
        let mut diag_upper = diag_origin as i32 + 1;
        let mut end1_reached = initial_matches == len1;
        let mut end2_reached = initial_matches == len2;

        // Track convergence
        let mut converged = false;

        // Use flag to track which array is "prev" (true = a is prev, false = b is prev)
        let mut use_a_as_prev = true;

        // For each distance (use max_dist for non-affine greedy)
        for d in 1..=max_dist {
            if diag_lower > diag_upper {
                converged = true;
                break; // Converged
            }

            // Set sentinel values on the "prev" array
            if use_a_as_prev {
                if diag_lower >= 1 {
                    mem.last_seq2_off_a[(diag_lower - 1) as usize] = -2;
                    mem.last_seq2_off_a[diag_lower as usize] = -2;
                }
                if (diag_upper as usize) < array_size - 1 {
                    mem.last_seq2_off_a[diag_upper as usize] = -2;
                    mem.last_seq2_off_a[(diag_upper + 1) as usize] = -2;
                }
            } else {
                if diag_lower >= 1 {
                    mem.last_seq2_off_b[(diag_lower - 1) as usize] = -2;
                    mem.last_seq2_off_b[diag_lower as usize] = -2;
                }
                if (diag_upper as usize) < array_size - 1 {
                    mem.last_seq2_off_b[diag_upper as usize] = -2;
                    mem.last_seq2_off_b[(diag_upper + 1) as usize] = -2;
                }
            }

            // X-drop score threshold (using scaled_xdrop for consistent units)
            let xdrop_idx = if d >= xdrop_offset {
                d - xdrop_offset
            } else {
                0
            };
            let xdrop_score = mem.max_score[xdrop_idx + xdrop_offset] + (op_cost * d as i32) - scaled_xdrop;
            let xdrop_score = (xdrop_score + match_cost / 2 - 1) / (match_cost / 2); // Ceiling division matching NCBI

            let mut curr_extent = 0i32;
            let mut curr_seq2_index = 0i32;
            let mut curr_diag = diag_origin as i32;
            let tmp_diag_lower = diag_lower;
            let tmp_diag_upper = diag_upper;

            // For each diagonal
            for k in tmp_diag_lower..=tmp_diag_upper {
                let ku = k as usize;
                if ku >= array_size || ku == 0 {
                    continue;
                }

                // Find largest seq2 offset that increases distance from d-1 to d
                // Access the "prev" array based on the flag
                let (prev_k_plus, prev_k, prev_k_minus) = if use_a_as_prev {
                    let p_plus = if ku + 1 < array_size { mem.last_seq2_off_a[ku + 1] } else { -2 };
                    let p = mem.last_seq2_off_a[ku];
                    let p_minus = if ku >= 1 { mem.last_seq2_off_a[ku - 1] } else { -2 };
                    (p_plus, p, p_minus)
                } else {
                    let p_plus = if ku + 1 < array_size { mem.last_seq2_off_b[ku + 1] } else { -2 };
                    let p = mem.last_seq2_off_b[ku];
                    let p_minus = if ku >= 1 { mem.last_seq2_off_b[ku - 1] } else { -2 };
                    (p_plus, p, p_minus)
                };

                let mut seq2_index = prev_k_plus.max(prev_k) + 1;
                seq2_index = seq2_index.max(prev_k_minus);

                let seq1_index = seq2_index + k - diag_origin as i32;

                if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                    // X-drop test failed or invalid diagonal
                    if k == diag_lower {
                        diag_lower += 1;
                    } else {
                        // Write to "curr" array
                        if use_a_as_prev {
                            mem.last_seq2_off_b[ku] = -2;
                        } else {
                            mem.last_seq2_off_a[ku] = -2;
                        }
                    }
                    continue;
                }

                diag_upper = k;

                // Slide down diagonal until mismatch
                let seq1_idx = seq1_index as usize;
                let seq2_idx = seq2_index as usize;

                if seq1_idx < len1 && seq2_idx < len2 {
                    let matches = find_first_mismatch_ex(q_seq, s_seq, len1, len2, seq1_idx, seq2_idx, reverse);
                    let new_seq1_index = seq1_index + matches as i32;
                    let new_seq2_index = seq2_index + matches as i32;

                    // Write to "curr" array
                    if use_a_as_prev {
                        mem.last_seq2_off_b[ku] = new_seq2_index;
                    } else {
                        mem.last_seq2_off_a[ku] = new_seq2_index;
                    }

                    // Track best extent
                    let extent = new_seq1_index + new_seq2_index;
                    if extent > curr_extent {
                        curr_extent = extent;
                        curr_seq2_index = new_seq2_index;
                        curr_diag = k;
                    }

                    // Clamp bounds
                    if new_seq2_index as usize >= len2 {
                        diag_lower = k + 1;
                        end2_reached = true;
                    }
                    if new_seq1_index as usize >= len1 {
                        diag_upper = k - 1;
                        end1_reached = true;
                    }
                } else {
                    // Write to "curr" array
                    if use_a_as_prev {
                        mem.last_seq2_off_b[ku] = seq2_index;
                    } else {
                        mem.last_seq2_off_a[ku] = seq2_index;
                    }
                }
            }

            // Compute max score for this distance
            let curr_score = (curr_extent * match_cost) / 2 - (d as i32) * op_cost;

            if curr_score >= mem.max_score[d - 1 + xdrop_offset] {
                mem.max_score[d + xdrop_offset] = curr_score;
                best_dist = d;
                best_seq2_len = curr_seq2_index as usize;
                best_seq1_len = (curr_seq2_index + curr_diag - diag_origin as i32) as usize;
            } else {
                mem.max_score[d + xdrop_offset] = mem.max_score[d - 1 + xdrop_offset];
            }

            // Check convergence
            if diag_lower > diag_upper {
                converged = true;
                break;
            }

            // Expand bounds for next distance
            if !end2_reached {
                diag_lower -= 1;
            }
            if !end1_reached {
                diag_upper += 1;
            }

            // Swap arrays by toggling the flag
            use_a_as_prev = !use_a_as_prev;
        }

        // Calculate final statistics
        // For greedy alignment, distance = mismatches + gaps
        // gap_letters is the difference in consumed lengths (indicates indels)
        let gap_letters = best_seq1_len.abs_diff(best_seq2_len);
        // If there are gap letters, there's at least one gap open
        // Without traceback, we estimate gap_opens = 1 if any gaps exist
        let gap_opens = if gap_letters > 0 { 1 } else { 0 };
        // Mismatches = distance - gap_letters (distance includes both mismatches and gaps)
        let mismatches = best_dist.saturating_sub(gap_letters);
        let matches = best_seq1_len.min(best_seq2_len).saturating_sub(mismatches);

        // Calculate score
        let score = (matches as i32) * reward + (mismatches as i32) * penalty;

        (
            (
                best_seq1_len,
                best_seq2_len,
                score,
                matches,
                mismatches,
                gap_opens,
                gap_letters,
            ),
            converged,
        )
    })
}

/// Affine greedy alignment function (NCBI BLAST's BLAST_AffineGreedyAlign).
/// This handles affine gap penalties properly by tracking three separate offsets
/// for each diagonal: paths ending in insertion, match, or deletion.
///
/// Returns: ((q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters), converged)
fn affine_greedy_align_one_direction_with_max_dist(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    in_gap_open: i32,
    in_gap_extend: i32,
    x_drop: i32,
    max_dist: usize,
) -> ((usize, usize, i32, usize, usize, usize, usize), bool) {
    let len1 = q_seq.len();
    let len2 = s_seq.len();

    if len1 == 0 || len2 == 0 {
        return ((0, 0, 0, 0, 0, 0, 0), true);
    }

    // Score normalization (NCBI BLAST approach):
    // Make sure bits of match_score don't disappear if divided by 2
    // IMPORTANT: In NCBI BLAST, mismatch_score is passed as -score_params->penalty,
    // where score_params->penalty is stored as a NEGATIVE value (e.g., -3).
    // So -(-3) = 3, meaning mismatch_score is the POSITIVE penalty magnitude.
    // In blemir, penalty is stored as a POSITIVE value (e.g., 3), so we use it directly.
    let (match_score, mismatch_penalty, xdrop_threshold, gap_open_in, gap_extend_in) =
        if reward % 2 == 1 {
            (reward * 2, penalty * 2, x_drop * 2, in_gap_open * 2, in_gap_extend * 2)
        } else {
            (reward, penalty, x_drop, in_gap_open, in_gap_extend)
        };

    // Fill in derived scores and penalties (NCBI BLAST approach)
    // op_cost = match_score + mismatch_penalty (both positive, so op_cost is positive)
    // This represents the "cost" of a mismatch in distance units
    let match_score_half = match_score / 2;
    let mut op_cost = match_score + mismatch_penalty;
    let mut gap_open = gap_open_in;
    let mut gap_extend = gap_extend_in + match_score_half;
    let score_common_factor = gdb3(&mut op_cost, &mut gap_open, &mut gap_extend);
    let gap_open_extend = gap_open + gap_extend;
    let max_penalty = op_cost.max(gap_open_extend);

    // Scaled max_dist for affine alignment
    // With the correct sign convention (op_cost positive), gap_extend should always be positive
    let scaled_max_dist = (max_dist as i32) * gap_extend;

    // Diagonal origin (center of diagonal space)
    let diag_origin = max_dist + 2;
    let array_size = 2 * diag_origin + 4;

    // For affine greedy, we need to track diag_lower and diag_upper for ALL distances
    // (not just current), because contributions can come from distances < d-1
    // With correct sign convention, max_penalty should always be positive
    
    // Debug assertions to catch sign issues
    debug_assert!(op_cost > 0, "op_cost must be positive, got {}", op_cost);
    debug_assert!(gap_extend > 0, "gap_extend must be positive, got {}", gap_extend);
    debug_assert!(max_penalty > 0, "max_penalty must be positive, got {}", max_penalty);
    debug_assert!(scaled_max_dist >= 0, "scaled_max_dist must be non-negative, got {}", scaled_max_dist);
    
    // Safe conversion from i32 to usize with bounds checking
    if max_penalty < 0 || scaled_max_dist < 0 {
        // Sign error - return early with non-convergence
        return ((0, 0, 0, 0, 0, 0, 0), false);
    }
    
    let max_penalty_usize = max_penalty as usize;
    let scaled_max_dist_usize = scaled_max_dist as usize;
    let bounds_size = scaled_max_dist_usize + 1 + max_penalty_usize + 1;
    
    let mut diag_lower_arr: Vec<i32> = vec![INVALID_DIAG; bounds_size];
    let mut diag_upper_arr: Vec<i32> = vec![-INVALID_DIAG; bounds_size];

    // Initialize negative distance bounds with empty ranges
    // (for distances < max_penalty, which map to indices 0..max_penalty_usize)
    for i in 0..max_penalty_usize {
        diag_lower_arr[i] = INVALID_DIAG;
        diag_upper_arr[i] = -INVALID_DIAG;
    }

    // last_seq2_off[d][k] stores GreedyOffset for diagonal k at distance d
    // We need to keep all rows for affine alignment (contributions from d - max_penalty)
    // Allocate enough rows: scaled_max_dist + max_penalty + 2
    let num_rows = scaled_max_dist_usize + max_penalty_usize + 2;
    
    // Check allocation sizes to prevent memory exhaustion
    let total_elements = num_rows.checked_mul(array_size);
    if total_elements.is_none() || total_elements.unwrap() > 100_000_000 {
        // Allocation would be too large, return early with non-convergence
        return ((0, 0, 0, 0, 0, 0, 0), false);
    }
    
    let mut last_seq2_off: Vec<Vec<GreedyOffset>> = vec![vec![GreedyOffset {
        insert_off: INVALID_OFFSET,
        match_off: INVALID_OFFSET,
        delete_off: INVALID_OFFSET,
    }; array_size]; num_rows];

    // Max score at each distance for X-drop
    let xdrop_offset = ((xdrop_threshold + match_score_half) / score_common_factor.max(1) + 1) as usize;
    let max_score_size = scaled_max_dist_usize + xdrop_offset + 2;
    let mut max_score: Vec<i32> = vec![0; max_score_size];

    // Find initial run of matches
    let initial_matches = find_first_mismatch(q_seq, s_seq, 0, 0);

    if initial_matches == len1 || initial_matches == len2 {
        // Perfect match - return immediately (converged)
        let score = (initial_matches as i32) * reward;
        return (
            (initial_matches, initial_matches, score, initial_matches, 0, 0, 0),
            true,
        );
    }

    // Initialize distance 0
    last_seq2_off[0][diag_origin].match_off = initial_matches as i32;
    last_seq2_off[0][diag_origin].insert_off = INVALID_OFFSET;
    last_seq2_off[0][diag_origin].delete_off = INVALID_OFFSET;
    max_score[xdrop_offset] = (initial_matches as i32) * match_score;
    diag_lower_arr[max_penalty_usize] = diag_origin as i32;
    diag_upper_arr[max_penalty_usize] = diag_origin as i32;

    let mut best_dist = 0i32;
    let mut best_diag = diag_origin as i32;
    let mut best_seq1_len = initial_matches;
    let mut best_seq2_len = initial_matches;
    let _ = best_diag; // Suppress unused warning for initial value

    // Set up for distance 1
    let mut curr_diag_lower = diag_origin as i32 - 1;
    let mut curr_diag_upper = diag_origin as i32 + 1;
    let mut end1_diag = 0i32;
    let mut end2_diag = 0i32;
    let mut num_nonempty_dist = 1i32;
    let mut d = 1i32;

    let mut converged = false;

    // Helper function to safely access diag bounds
    let get_diag_lower = |arr: &[i32], d: i32, max_pen: usize| -> i32 {
        let idx = (d + max_pen as i32) as usize;
        if idx < arr.len() { arr[idx] } else { INVALID_DIAG }
    };
    let get_diag_upper = |arr: &[i32], d: i32, max_pen: usize| -> i32 {
        let idx = (d + max_pen as i32) as usize;
        if idx < arr.len() { arr[idx] } else { -INVALID_DIAG }
    };

    // For each distance
    while d <= scaled_max_dist {
        // Compute X-dropoff score threshold
        let xdrop_idx = if d as usize >= xdrop_offset {
            (d as usize) - xdrop_offset
        } else {
            0
        };
        let xdrop_score_raw = max_score[xdrop_idx + xdrop_offset] 
            + score_common_factor * d - xdrop_threshold;
        let xdrop_score = ((xdrop_score_raw as f64) / (match_score_half as f64)).ceil() as i32;
        let xdrop_score = xdrop_score.max(0);

        let mut curr_extent = 0i32;
        let mut curr_seq2_index = 0i32;
        let mut curr_diag = 0i32;
        let tmp_diag_lower = curr_diag_lower;
        let tmp_diag_upper = curr_diag_upper;

        // For each valid diagonal
        for k in tmp_diag_lower..=tmp_diag_upper {
            let ku = k as usize;
            if ku >= array_size || ku == 0 {
                continue;
            }

            // Find best offset for DELETE (gap in seq1) - look at k+1 diagonal
            let mut seq2_index_del = INVALID_OFFSET;
            
            // From gap opening (match -> delete)
            let d_open = d - gap_open_extend;
            if d_open >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_open, max_penalty_usize);
                let du = get_diag_upper(&diag_upper_arr, d_open, max_penalty_usize);
                if k + 1 >= dl && k + 1 <= du && (d_open as usize) < last_seq2_off.len() {
                    let ku1 = (k + 1) as usize;
                    if ku1 < array_size {
                        seq2_index_del = last_seq2_off[d_open as usize][ku1].match_off;
                    }
                }
            }
            
            // From gap extension (delete -> delete)
            let d_ext = d - gap_extend;
            if d_ext >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_ext, max_penalty_usize);
                let du = get_diag_upper(&diag_upper_arr, d_ext, max_penalty_usize);
                if k + 1 >= dl && k + 1 <= du && (d_ext as usize) < last_seq2_off.len() {
                    let ku1 = (k + 1) as usize;
                    if ku1 < array_size {
                        let ext_off = last_seq2_off[d_ext as usize][ku1].delete_off;
                        if ext_off > seq2_index_del {
                            seq2_index_del = ext_off;
                        }
                    }
                }
            }
            
            // Save delete offset (deletion means seq2 offset slips by one)
            let du = d as usize;
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].delete_off = if seq2_index_del == INVALID_OFFSET {
                    INVALID_OFFSET
                } else {
                    seq2_index_del + 1
                };
            }

            // Find best offset for INSERT (gap in seq2) - look at k-1 diagonal
            let mut seq2_index_ins = INVALID_OFFSET;
            
            // From gap opening (match -> insert)
            if d_open >= 0 && k >= 1 {
                let dl = get_diag_lower(&diag_lower_arr, d_open, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_open, max_penalty_usize);
                if k - 1 >= dl && k - 1 <= du_bound && (d_open as usize) < last_seq2_off.len() {
                    let km1 = (k - 1) as usize;
                    if km1 < array_size {
                        seq2_index_ins = last_seq2_off[d_open as usize][km1].match_off;
                    }
                }
            }
            
            // From gap extension (insert -> insert)
            if d_ext >= 0 && k >= 1 {
                let dl = get_diag_lower(&diag_lower_arr, d_ext, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_ext, max_penalty_usize);
                if k - 1 >= dl && k - 1 <= du_bound && (d_ext as usize) < last_seq2_off.len() {
                    let km1 = (k - 1) as usize;
                    if km1 < array_size {
                        let ext_off = last_seq2_off[d_ext as usize][km1].insert_off;
                        if ext_off > seq2_index_ins {
                            seq2_index_ins = ext_off;
                        }
                    }
                }
            }
            
            // Save insert offset (insertion doesn't change seq2 offset)
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].insert_off = seq2_index_ins;
            }

            // Compare with mismatch path (from diagonal k at d - op_cost)
            let mut seq2_index = last_seq2_off[du][ku].insert_off.max(last_seq2_off[du][ku].delete_off);
            
            let d_mismatch = d - op_cost;
            if d_mismatch >= 0 {
                let dl = get_diag_lower(&diag_lower_arr, d_mismatch, max_penalty_usize);
                let du_bound = get_diag_upper(&diag_upper_arr, d_mismatch, max_penalty_usize);
                if k >= dl && k <= du_bound && (d_mismatch as usize) < last_seq2_off.len() {
                    let match_off = last_seq2_off[d_mismatch as usize][ku].match_off;
                    if match_off != INVALID_OFFSET {
                        seq2_index = seq2_index.max(match_off + 1);
                    }
                }
            }

            // Choose seq1 offset to remain on diagonal k
            let seq1_index = seq2_index + k - diag_origin as i32;

            // X-dropoff test
            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == curr_diag_lower {
                    curr_diag_lower += 1;
                } else if du < last_seq2_off.len() {
                    last_seq2_off[du][ku].match_off = INVALID_OFFSET;
                }
                continue;
            }
            curr_diag_upper = k;

            // Slide down diagonal until mismatch
            let seq1_idx = seq1_index as usize;
            let seq2_idx = seq2_index as usize;

            let matches = if seq1_idx < len1 && seq2_idx < len2 {
                find_first_mismatch(q_seq, s_seq, seq1_idx, seq2_idx)
            } else {
                0
            };

            let new_seq1_index = seq1_index + matches as i32;
            let new_seq2_index = seq2_index + matches as i32;

            // Save match offset
            if du < last_seq2_off.len() {
                last_seq2_off[du][ku].match_off = new_seq2_index;
            }

            // Track best extent
            let extent = new_seq1_index + new_seq2_index;
            if extent > curr_extent {
                curr_extent = extent;
                curr_seq2_index = new_seq2_index;
                curr_diag = k;
            }

            // Clamp bounds to avoid walking off sequences
            if new_seq1_index as usize >= len1 {
                curr_diag_upper = k;
                end1_diag = k - 1;
            }
            if new_seq2_index as usize >= len2 {
                curr_diag_lower = k;
                end2_diag = k + 1;
            }
        }

        // Compute maximum score for distance d
        let curr_score = curr_extent * match_score_half - d * score_common_factor;

        // Update best if this is better
        let prev_max = if d >= 1 && (d - 1) as usize + xdrop_offset < max_score.len() {
            max_score[(d - 1) as usize + xdrop_offset]
        } else {
            0
        };

        if curr_score > prev_max {
            if (d as usize) + xdrop_offset < max_score.len() {
                max_score[(d as usize) + xdrop_offset] = curr_score;
            }
            best_dist = d;
            best_diag = curr_diag;
            best_seq2_len = curr_seq2_index as usize;
            best_seq1_len = (curr_seq2_index + best_diag - diag_origin as i32) as usize;
        } else if (d as usize) + xdrop_offset < max_score.len() {
            max_score[(d as usize) + xdrop_offset] = prev_max;
        }

        // Save diagonal bounds for this distance
        let bounds_idx = (d + max_penalty) as usize;
        if bounds_idx < diag_lower_arr.len() {
            if curr_diag_lower <= curr_diag_upper {
                num_nonempty_dist += 1;
                diag_lower_arr[bounds_idx] = curr_diag_lower;
                diag_upper_arr[bounds_idx] = curr_diag_upper;
            } else {
                diag_lower_arr[bounds_idx] = INVALID_DIAG;
                diag_upper_arr[bounds_idx] = -INVALID_DIAG;
            }
        }

        // Check if we should decrement num_nonempty_dist
        let old_bounds_idx = (d - max_penalty + max_penalty) as usize;
        if old_bounds_idx < diag_lower_arr.len() {
            let old_lower = diag_lower_arr[old_bounds_idx];
            let old_upper = diag_upper_arr[old_bounds_idx];
            if old_lower <= old_upper {
                num_nonempty_dist -= 1;
            }
        }

        // Convergence check: max_penalty consecutive empty ranges
        if num_nonempty_dist == 0 {
            converged = true;
            break;
        }

        // Compute diagonal range for next distance
        d += 1;
        
        let d_goe = d - gap_open_extend;
        let d_ge = d - gap_extend;
        let d_op = d - op_cost;

        let lower_goe = get_diag_lower(&diag_lower_arr, d_goe, max_penalty_usize);
        let lower_ge = get_diag_lower(&diag_lower_arr, d_ge, max_penalty_usize);
        let lower_op = get_diag_lower(&diag_lower_arr, d_op, max_penalty_usize);
        
        curr_diag_lower = lower_goe.min(lower_ge) - 1;
        curr_diag_lower = curr_diag_lower.min(lower_op);
        
        if end2_diag > 0 {
            curr_diag_lower = curr_diag_lower.max(end2_diag);
        }

        let upper_goe = get_diag_upper(&diag_upper_arr, d_goe, max_penalty_usize);
        let upper_ge = get_diag_upper(&diag_upper_arr, d_ge, max_penalty_usize);
        let upper_op = get_diag_upper(&diag_upper_arr, d_op, max_penalty_usize);
        
        curr_diag_upper = upper_goe.max(upper_ge) + 1;
        curr_diag_upper = curr_diag_upper.max(upper_op);
        
        if end1_diag > 0 {
            curr_diag_upper = curr_diag_upper.min(end1_diag);
        }
    }

    if !converged {
        // Did not converge - return best result found so far
        // Calculate statistics from best alignment
        let alignment_len = best_seq1_len.max(best_seq2_len);
        let gap_letters = best_seq1_len.abs_diff(best_seq2_len);
        
        // Estimate statistics (without full traceback)
        // For affine greedy, we estimate based on distance and gap penalties
        let estimated_gaps = if gap_letters > 0 { 1 } else { 0 };
        let estimated_mismatches = (best_dist as usize).saturating_sub(estimated_gaps);
        let matches = alignment_len.saturating_sub(estimated_mismatches).saturating_sub(gap_letters);
        
        let score = (matches as i32) * reward 
            + (estimated_mismatches as i32) * penalty
            - (estimated_gaps as i32) * in_gap_open
            - (gap_letters as i32) * in_gap_extend;
        
        return (
            (best_seq1_len, best_seq2_len, score, matches, estimated_mismatches, estimated_gaps, gap_letters),
            false,
        );
    }

    // Calculate final statistics
    // For affine greedy, we need to estimate statistics from the alignment
    // Since we don't have full traceback, we estimate based on the distance and gap penalties
    let alignment_len = best_seq1_len.max(best_seq2_len);
    let gap_letters = best_seq1_len.abs_diff(best_seq2_len);
    
    // Estimate gap opens and mismatches from distance
    // In affine greedy: distance = sum of (gap_open_extend for each gap) + (gap_extend for each additional gap letter) + (op_cost for each mismatch)
    // We approximate: if there are gaps, assume one gap open
    let estimated_gap_opens = if gap_letters > 0 { 1 } else { 0 };
    
    // Remaining distance after accounting for gaps
    let gap_cost = if gap_letters > 0 {
        gap_open_extend + (gap_letters.saturating_sub(1) as i32) * gap_extend
    } else {
        0
    };
    let remaining_dist = (best_dist - gap_cost).max(0);
    let estimated_mismatches = if op_cost > 0 {
        (remaining_dist / op_cost) as usize
    } else {
        0
    };
    
    let matches = alignment_len.saturating_sub(estimated_mismatches).saturating_sub(gap_letters);
    
    // Calculate final score
    let score = (matches as i32) * reward 
        + (estimated_mismatches as i32) * penalty
        - (estimated_gap_opens as i32) * in_gap_open
        - (gap_letters as i32) * in_gap_extend;

    (
        (best_seq1_len, best_seq2_len, score, matches, estimated_mismatches, estimated_gap_opens, gap_letters),
        converged,
    )
}

#[allow(dead_code)]
fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Extend alignment in both directions using NCBI-style gapped extension.
///
/// This function uses a two-phase approach:
/// 1. First try greedy alignment (fast for high-identity sequences)
/// 2. Fall back to full DP if greedy doesn't produce good results
///
/// This function accepts an x_drop parameter to control extension termination.
/// For two-phase extension, call this twice with different x_drop values.
///
/// Returns: (q_start, q_end, s_start, s_end, score, matches, mismatches, gap_opens, gap_letters, dp_cells)
///
/// This function uses a hybrid approach for performance:
/// 1. Always try greedy alignment first (fast for high-identity sequences)
/// 2. For affine gap penalties, only fall back to full DP if greedy result is poor
fn extend_gapped_heuristic(
    q_seq: &[u8],
    s_seq: &[u8],
    qs: usize,
    ss: usize,
    len: usize,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    use_dp: bool, // If true, use DP-based extension (for blastn task); if false, use greedy (for megablast)
) -> (
    usize,
    usize,
    usize,
    usize,
    i32,
    usize,
    usize,
    usize,
    usize,
    usize,
) {
    // Bounds validation: ensure seed coordinates are valid
    if qs >= q_seq.len() || ss >= s_seq.len() {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, 0);
    }

    // Clamp len to available sequence length
    let len = len.min(q_seq.len() - qs).min(s_seq.len() - ss);
    if len == 0 {
        return (qs, qs, ss, ss, 0, 0, 0, 0, 0, 0);
    }

    // NCBI BLAST uses different extension algorithms based on task:
    // - megablast: greedy alignment (fast, good for high-identity sequences)
    // - blastn: DP-based alignment (slower, but handles divergent sequences better)
    // This follows NCBI BLAST's approach where blastn uses eDynProgScoreOnly
    // and megablast uses eGreedyScoreOnly.

    // First, extend to the right from the seed
    let right_suffix_q = q_seq.get(qs + len..).unwrap_or(&[]);
    let right_suffix_s = s_seq.get(ss + len..).unwrap_or(&[]);
    
    let (
        right_q_consumed,
        right_s_consumed,
        right_score,
        right_matches,
        right_mismatches,
        right_gaps,
        right_gap_letters,
        right_dp_cells,
    ) = if use_dp {
        // DP-based extension for blastn task (handles divergent sequences better)
        extend_gapped_one_direction(
            right_suffix_q,
            right_suffix_s,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
        )
    } else {
        // Greedy extension for megablast task (faster for high-identity sequences)
        let (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters) =
            greedy_align_one_direction(
                right_suffix_q,
                right_suffix_s,
                reward,
                penalty,
                gap_open,
                gap_extend,
                x_drop,
            );
        (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters, 0)
    };

    // Then extend to the left from the seed
    // For greedy alignment, use reverse flag to avoid O(N) prefix copying
    // For DP-based alignment, still need to copy (DP doesn't support reverse flag yet)
    let (
        left_q_consumed,
        left_s_consumed,
        left_score,
        left_matches,
        left_mismatches,
        left_gaps,
        left_gap_letters,
        left_dp_cells,
    ) = if use_dp {
        // DP-based extension for blastn task - still needs prefix copying
        let q_prefix: Vec<u8> = q_seq[..qs].iter().rev().copied().collect();
        let s_prefix: Vec<u8> = s_seq[..ss].iter().rev().copied().collect();
        extend_gapped_one_direction(
            &q_prefix,
            &s_prefix,
            reward,
            penalty,
            gap_open,
            gap_extend,
            x_drop,
        )
    } else {
        // Greedy extension for megablast task - use reverse flag to avoid copying
        // Pass the full sequences with the prefix length and reverse=true
        let (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters) =
            greedy_align_one_direction_ex(
                q_seq,
                s_seq,
                qs,  // len1 = prefix length (characters before seed)
                ss,  // len2 = prefix length (characters before seed)
                reward,
                penalty,
                gap_open,
                gap_extend,
                x_drop,
                true,  // reverse = true for left extension
            );
        (q_cons, s_cons, score, matches, mismatches, gaps, gap_letters, 0)
    };

    // Calculate seed score
    let mut seed_score = 0;
    let mut seed_matches = 0;
    let mut seed_mismatches = 0;
    for k in 0..len {
        if q_seq[qs + k] == s_seq[ss + k] {
            seed_score += reward;
            seed_matches += 1;
        } else {
            seed_score += penalty;
            seed_mismatches += 1;
        }
    }

    let total_score = left_score + seed_score + right_score;
    let total_matches = left_matches + seed_matches + right_matches;
    let total_mismatches = left_mismatches + seed_mismatches + right_mismatches;
    let total_gaps = left_gaps + right_gaps;
    let total_gap_letters = left_gap_letters + right_gap_letters;
    let total_dp_cells = left_dp_cells + right_dp_cells;

    // Calculate final positions
    let final_q_start = qs - left_q_consumed;
    let final_q_end = qs + len + right_q_consumed;
    let final_s_start = ss - left_s_consumed;
    let final_s_end = ss + len + right_s_consumed;

    (
        final_q_start,
        final_q_end,
        final_s_start,
        final_s_end,
        total_score,
        total_matches,
        total_mismatches,
        total_gaps,
        total_gap_letters,
        total_dp_cells,
    )
}

/// Alignment statistics propagated alongside DP scores
#[derive(Clone, Copy, Default)]
struct AlnStats {
    matches: u32,
    mismatches: u32,
    gap_opens: u32,
    gap_letters: u32, // Total gap characters (for alignment length calculation)
}

/// Extend alignment in one direction using banded Smith-Waterman with affine gap penalties.
///
/// This implements a proper row-by-row DP algorithm with counts propagation for accurate statistics.
/// Instead of storing full traceback matrices, we propagate alignment statistics (matches, mismatches,
/// gap_opens) alongside the DP scores, keeping only 2 rows at a time for O(band_size) memory.
///
/// - Row i corresponds to query position i
/// - Band index k represents diagonal offset: j = i + (k - W) where W is half-bandwidth
/// - Three matrices track different alignment states:
///   - M[i,j]: best score ending with a match/mismatch at (i,j)
///   - Ix[i,j]: best score ending with a gap in subject (deletion from query)
///   - Iy[i,j]: best score ending with a gap in query (insertion in query)
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters)
/// Extend alignment in one direction using NCBI-style semi-gapped DP with affine gap penalties.
///
/// This implements NCBI BLAST's Blast_SemiGappedAlign approach:
/// - X-drop based dynamic window that expands/contracts based on score
/// - Tracks best score across ALL diagonals
/// - Propagates alignment statistics for accurate traceback-based calculation
/// - No hard-coded extension limits (controlled by X-drop termination)
///
/// Returns: (q_consumed, s_consumed, score, matches, mismatches, gap_opens, gap_letters, dp_cells)
/// The dp_cells count is for diagnostic purposes.
fn extend_gapped_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
) -> (usize, usize, i32, usize, usize, usize, usize, usize) {
    // NCBI BLAST-style adaptive banding implementation
    // Key insight from Blast_SemiGappedAlign:
    // - Use dynamic window bounds (first_b_index to b_size) that expand/contract based on X-drop
    // - Window naturally expands as needed, but with a maximum limit to prevent explosion
    // - Reallocate memory when window needs to grow

    const NEG_INF: i32 = i32::MIN / 2;
    // Maximum window size to prevent computational explosion on very long high-identity alignments
    // Increased from 5000 to 50000 to allow alignments with more gaps (NCBI BLAST can produce
    // alignments with 800+ gap characters that require wider bands to traverse)
    // Note: NCBI BLAST uses additional mechanisms (greedy alignment, fences) that we don't have yet
    const MAX_WINDOW_SIZE: usize = 50000;

    let m = q_seq.len();
    let n = s_seq.len();

    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0);
    }

    // Calculate initial window size based on X-drop (NCBI-style)
    // num_extra_cells = x_dropoff / gap_extend + 3
    let gap_extend_abs = gap_extend.abs().max(1);
    let initial_window = (x_drop / gap_extend_abs + 3) as usize;

    // Start with a reasonable initial allocation, will grow as needed (up to MAX_WINDOW_SIZE)
    let mut alloc_size = initial_window.max(100).min(MAX_WINDOW_SIZE);

    // DP score arrays - indexed by j (subject position)
    // We use a 1D array approach like NCBI BLAST
    #[derive(Clone, Copy, Default)]
    struct DpCell {
        best: i32,
        best_gap: i32, // best score ending in a gap
        stats: AlnStats,
        gap_stats: AlnStats,
    }

    let mut score_array: Vec<DpCell> = vec![
        DpCell {
            best: NEG_INF,
            best_gap: NEG_INF,
            stats: AlnStats::default(),
            gap_stats: AlnStats::default()
        };
        alloc_size
    ];

    // Initialize row 0
    let gap_open_extend = gap_open + gap_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = gap_open_extend; // Cost to open gap at position 0

    // Initialize leading gaps in subject (j > 0)
    let mut score = gap_open_extend;
    let mut b_size = 1usize;
    for j in 1..=n.min(alloc_size - 1) {
        if score < -x_drop {
            break;
        }
        score_array[j].best = score;
        score_array[j].best_gap = score + gap_open_extend;
        score_array[j].stats = AlnStats {
            matches: 0,
            mismatches: 0,
            gap_opens: 1,
            gap_letters: j as u32,
        };
        score += gap_extend;
        b_size = j + 1;
    }

    // Track best score and position
    let mut best_score = 0;
    let mut best_i = 0;
    let mut best_j = 0;
    let mut best_stats = AlnStats::default();

    // Dynamic window bounds (NCBI-style)
    let mut first_b_index = 0usize;

    // DP cell counter for diagnostics
    let mut dp_cells = 0usize;

    // Process each row (query position)
    for i in 1..=m {
        let qc = q_seq[i - 1];

        // Running scores for this row
        let mut score_val = NEG_INF;
        let mut score_gap_row = NEG_INF; // Best score ending in gap in query (Ix)
        let mut score_gap_row_stats = AlnStats::default();
        let mut score_stats = AlnStats::default();
        let mut last_b_index = first_b_index;

        for j in first_b_index..b_size {
            let sc = if j < n { s_seq[j] } else { break };
            dp_cells += 1;

            // Get previous column's gap score
            let score_gap_col = score_array[j].best_gap;
            let score_gap_col_stats = score_array[j].gap_stats;

            // Compute match/mismatch score
            let is_match = qc == sc;
            let match_score = if is_match { reward } else { penalty };
            let next_score = if score_array[j].best > NEG_INF {
                score_array[j].best + match_score
            } else {
                NEG_INF
            };
            let mut next_stats = score_array[j].stats;
            if is_match {
                next_stats.matches += 1;
            } else {
                next_stats.mismatches += 1;
            }

            // Best of: continue from M, continue from Ix (gap in query), continue from Iy (gap in subject)
            if score_val < score_gap_col {
                score_val = score_gap_col;
                score_stats = score_gap_col_stats;
            }
            if score_val < score_gap_row {
                score_val = score_gap_row;
                score_stats = score_gap_row_stats;
            }

            // X-drop check
            if best_score - score_val > x_drop {
                // Failed X-drop - mark this cell as invalid
                if j == first_b_index {
                    first_b_index += 1;
                } else {
                    score_array[j].best = NEG_INF;
                }
            } else {
                last_b_index = j;

                // Update best score
                if score_val > best_score {
                    best_score = score_val;
                    best_i = i;
                    best_j = j + 1; // Convert to 1-based
                    best_stats = score_stats;
                }

                // Update gap scores for next iteration
                // Gap in query (Ix): extend existing gap or open new gap
                let extend_gap_row = score_gap_row + gap_extend;
                let open_gap_row = score_val + gap_open_extend;
                if open_gap_row > extend_gap_row {
                    score_gap_row = open_gap_row;
                    score_gap_row_stats = score_stats;
                    score_gap_row_stats.gap_opens += 1;
                    score_gap_row_stats.gap_letters += 1;
                } else {
                    score_gap_row = extend_gap_row;
                    score_gap_row_stats.gap_letters += 1;
                }

                // Gap in subject (Iy): extend existing gap or open new gap
                let extend_gap_col = score_gap_col + gap_extend;
                let open_gap_col = score_val + gap_open_extend;
                if open_gap_col > extend_gap_col {
                    score_array[j].best_gap = open_gap_col;
                    score_array[j].gap_stats = score_stats;
                    score_array[j].gap_stats.gap_opens += 1;
                    score_array[j].gap_stats.gap_letters += 1;
                } else {
                    score_array[j].best_gap = extend_gap_col;
                    score_array[j].gap_stats.gap_letters += 1;
                }

                // Store current score
                score_array[j].best = score_val;
                score_array[j].stats = score_stats;
            }

            // Move to next cell
            score_val = next_score;
            score_stats = next_stats;
        }

        // Check if all positions failed X-drop
        if first_b_index >= b_size {
            break;
        }

        // Expand window if needed (NCBI-style adaptive banding, with MAX_WINDOW_SIZE limit)
        if last_b_index + initial_window + 3 >= alloc_size && alloc_size < MAX_WINDOW_SIZE {
            // Need to grow the array (but respect MAX_WINDOW_SIZE)
            let new_alloc = (last_b_index + initial_window + 100)
                .max(alloc_size * 2)
                .min(MAX_WINDOW_SIZE);
            score_array.resize(
                new_alloc,
                DpCell {
                    best: NEG_INF,
                    best_gap: NEG_INF,
                    stats: AlnStats::default(),
                    gap_stats: AlnStats::default(),
                },
            );
            alloc_size = new_alloc;
        }

        if last_b_index < b_size.saturating_sub(1) {
            // This row ended earlier than last row - shrink window
            b_size = last_b_index + 1;
        } else {
            // Extend window if we can continue with gaps (respect MAX_WINDOW_SIZE)
            while score_gap_row >= best_score - x_drop
                && b_size <= n
                && b_size < alloc_size
                && b_size < MAX_WINDOW_SIZE
            {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].stats = score_gap_row_stats;
                score_array[b_size].best_gap = score_gap_row + gap_open_extend;
                score_array[b_size].gap_stats = score_gap_row_stats;
                score_array[b_size].gap_stats.gap_opens += 1;
                score_array[b_size].gap_stats.gap_letters += 1;
                score_gap_row += gap_extend;
                score_gap_row_stats.gap_letters += 1;
                b_size += 1;
            }
        }

        // Ensure we have a sentinel
        if b_size <= n && b_size < alloc_size {
            score_array[b_size].best = NEG_INF;
            score_array[b_size].best_gap = NEG_INF;
            b_size += 1;
        }
    }

    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, dp_cells);
    }

    (
        best_i,
        best_j,
        best_score,
        best_stats.matches as usize,
        best_stats.mismatches as usize,
        best_stats.gap_opens as usize,
        best_stats.gap_letters as usize,
        dp_cells,
    )
}

fn calculate_evalue(
    score: i32,
    q_len: usize,
    db_len: usize,
    db_num_seqs: usize,
    params: &KarlinParams,
) -> (f64, f64) {
    // Use NCBI-compatible length adjustment for database search
    let search_space = SearchSpace::for_database_search(q_len, db_len, db_num_seqs, params, true);
    let bs = calc_bit_score(score, params);
    let ev = calc_evalue(bs, &search_space);
    (bs, ev)
}

/// Merge overlapping HSPs and filter redundant hits.
///
/// This function:
/// 1. Groups hits by query-subject pair
/// 2. Sorts hits by bit score (descending)
/// 3. Removes hits that overlap significantly with higher-scoring hits
/// 4. Returns the filtered list of non-redundant hits
/// Align a region using banded Smith-Waterman (for cluster-then-extend chaining).
/// Unlike extend_gapped_heuristic which extends from a seed, this aligns the entire region.
/// Returns: (q_start, q_end, s_start, s_end, score, matches, mismatches, gap_opens, gap_letters)
/// where q_start/q_end and s_start/s_end are 0-based coordinates within the input slices.
///
/// TRUE LOCAL Smith-Waterman with affine gap penalties.
///
/// Key difference from previous version: includes max(0, ...) reset to allow alignments
/// to start fresh at any position. This prevents "bridging" through low-quality regions
/// and produces alignments more similar to NCBI BLAST.
///
/// OPTIMIZED VERSION: Uses rolling score rows (O(band_size) memory for scores) + compact traceback
/// instead of storing full DpCell structs for every cell. This reduces memory from ~300MB to ~15MB
/// for a 10kb region, and eliminates massive allocation overhead.
#[allow(dead_code)]
fn align_region(
    q_seq: &[u8],
    s_seq: &[u8],
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
) -> (usize, usize, usize, usize, i32, usize, usize, usize, usize) {
    const BAND_WIDTH: usize = 256; // Wider band for NCBI BLAST compatibility
    const NEG_INF: i32 = i32::MIN / 2;

    let m = q_seq.len();
    let n = s_seq.len();

    if m == 0 || n == 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    let band_size = 2 * BAND_WIDTH + 1;

    // Traceback state encoding
    // 0 = STOP (alignment starts here - from max(0,...) reset)
    // 1 = from M (diagonal), 2 = from Ix (up), 3 = from Iy (left)
    const TB_STOP: u8 = 0;
    const TB_DIAG: u8 = 1;
    const TB_UP: u8 = 2;
    const TB_LEFT: u8 = 3;

    // Rolling score arrays for affine gap DP (only need prev and current row)
    // Each array has band_size elements
    let mut prev_m: Vec<i32> = vec![NEG_INF; band_size];
    let mut prev_ix: Vec<i32> = vec![NEG_INF; band_size];
    let mut prev_iy: Vec<i32> = vec![NEG_INF; band_size];
    let mut cur_m: Vec<i32> = vec![NEG_INF; band_size];
    let mut cur_ix: Vec<i32> = vec![NEG_INF; band_size];
    let mut cur_iy: Vec<i32> = vec![NEG_INF; band_size];

    // Compact traceback storage: for each cell, store which state gave the best score
    // Layout: flat array indexed by i * band_size + k
    let tb_size = (m + 1) * band_size;
    let mut tb_state: Vec<u8> = vec![TB_STOP; tb_size];

    // Track best score position
    let mut best_score = 0;
    let mut best_i = 0;
    let mut best_j = 0;

    // Initialize row 0: M[0,0] = 0 (can start alignment here)
    // For true local alignment, we don't initialize leading gaps - they would be negative
    // and we'd rather start fresh (score 0) than pay gap penalties
    prev_m[BAND_WIDTH] = 0;
    // tb_state[BAND_WIDTH] = TB_STOP; // Already initialized to STOP

    // Process each row
    for i in 1..=m {
        // Reset current row
        for k in 0..band_size {
            cur_m[k] = NEG_INF;
            cur_ix[k] = NEG_INF;
            cur_iy[k] = NEG_INF;
        }

        // Determine valid k range
        let k_min = if i > BAND_WIDTH { 0 } else { BAND_WIDTH - i };
        let k_max = if n + BAND_WIDTH >= i {
            (n + BAND_WIDTH - i).min(band_size - 1)
        } else {
            0
        };

        if k_min > k_max {
            // Swap rows
            std::mem::swap(&mut prev_m, &mut cur_m);
            std::mem::swap(&mut prev_ix, &mut cur_ix);
            std::mem::swap(&mut prev_iy, &mut cur_iy);
            continue;
        }

        for k in k_min..=k_max {
            let j_offset = k as isize - BAND_WIDTH as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > n {
                continue;
            }

            let q_base = q_seq[i - 1];
            let s_base = s_seq[j - 1];
            let match_score = if q_base == s_base { reward } else { penalty };

            // M[i,j] = max(0, M[i-1,j-1] + match, Ix[i-1,j-1] + match, Iy[i-1,j-1] + match)
            // The max(0, ...) is the key for TRUE LOCAL alignment
            let diag_score = prev_m[k].max(prev_ix[k]).max(prev_iy[k]);
            let m_score = if diag_score > NEG_INF {
                (diag_score + match_score).max(0)
            } else {
                // Can always start fresh with a match (if positive) or 0
                match_score.max(0)
            };
            cur_m[k] = m_score;

            // Determine traceback for M state
            let tb_idx = i * band_size + k;
            if m_score == 0 {
                // Started fresh here (max(0,...) chose 0)
                tb_state[tb_idx] = TB_STOP;
            } else if diag_score > NEG_INF && m_score == diag_score + match_score {
                // Came from previous cell
                if prev_m[k] >= prev_ix[k] && prev_m[k] >= prev_iy[k] {
                    tb_state[tb_idx] = TB_DIAG;
                } else if prev_ix[k] >= prev_iy[k] {
                    tb_state[tb_idx] = TB_UP;
                } else {
                    tb_state[tb_idx] = TB_LEFT;
                }
            } else {
                // Started fresh with just match_score (no valid predecessor)
                tb_state[tb_idx] = TB_STOP;
            }

            // Ix[i,j] = max(M[i-1,j] + gap_open + gap_extend, Ix[i-1,j] + gap_extend)
            // Note: For local alignment, we don't apply max(0,...) to gap states
            // because gaps can only extend existing alignments, not start new ones
            if k + 1 < band_size {
                let open_score = if prev_m[k + 1] > 0 {
                    prev_m[k + 1] + gap_open + gap_extend
                } else {
                    NEG_INF
                };
                let extend_score = if prev_ix[k + 1] > NEG_INF {
                    prev_ix[k + 1] + gap_extend
                } else {
                    NEG_INF
                };
                cur_ix[k] = open_score.max(extend_score);
            }

            // Iy[i,j] = max(M[i,j-1] + gap_open + gap_extend, Iy[i,j-1] + gap_extend)
            if k > 0 {
                let open_score = if cur_m[k - 1] > 0 {
                    cur_m[k - 1] + gap_open + gap_extend
                } else {
                    NEG_INF
                };
                let extend_score = if cur_iy[k - 1] > NEG_INF {
                    cur_iy[k - 1] + gap_extend
                } else {
                    NEG_INF
                };
                cur_iy[k] = open_score.max(extend_score);
            }

            // Track best score (Smith-Waterman: best anywhere in matrix)
            // Only consider M state for best score (gaps can't end an alignment optimally)
            if cur_m[k] > best_score {
                best_score = cur_m[k];
                best_i = i;
                best_j = j;
            }
        }

        // Swap rows
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_ix, &mut cur_ix);
        std::mem::swap(&mut prev_iy, &mut cur_iy);
    }

    if best_score <= 0 {
        return (0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    // Traceback to find alignment start and compute stats
    // Start from best_i, best_j and trace back until we hit TB_STOP
    let mut cur_i = best_i;
    let mut cur_j = best_j;
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_letters = 0usize;
    let mut in_gap_ix = false;
    let mut in_gap_iy = false;

    // We need to track which state we're in during traceback
    // Start in M state (best score is always from M)
    #[derive(Clone, Copy, PartialEq)]
    enum TbState {
        M,
        Ix,
        Iy,
    }
    let mut current_state = TbState::M;

    while cur_i > 0 && cur_j > 0 {
        let k = (cur_j as isize - cur_i as isize + BAND_WIDTH as isize) as usize;
        if k >= band_size {
            break;
        }

        let tb_idx = cur_i * band_size + k;
        let state = tb_state[tb_idx];

        // Stop if we've reached a STOP state (alignment boundary from max(0,...))
        if state == TB_STOP {
            // Count the current position as part of alignment
            if q_seq[cur_i - 1] == s_seq[cur_j - 1] {
                matches += 1;
            } else {
                mismatches += 1;
            }
            cur_i -= 1;
            cur_j -= 1;
            break;
        }

        match current_state {
            TbState::M => {
                // We're in M state, count match/mismatch
                if q_seq[cur_i - 1] == s_seq[cur_j - 1] {
                    matches += 1;
                } else {
                    mismatches += 1;
                }

                // Determine where we came from
                match state {
                    TB_DIAG => {
                        // Came from M[i-1,j-1]
                        current_state = TbState::M;
                        cur_i -= 1;
                        cur_j -= 1;
                    }
                    TB_UP => {
                        // Came from Ix[i-1,j-1]
                        current_state = TbState::Ix;
                        cur_i -= 1;
                        cur_j -= 1;
                    }
                    TB_LEFT => {
                        // Came from Iy[i-1,j-1]
                        current_state = TbState::Iy;
                        cur_i -= 1;
                        cur_j -= 1;
                    }
                    _ => break,
                }
            }
            TbState::Ix => {
                // We're in Ix state (gap in subject)
                gap_letters += 1;
                if !in_gap_ix {
                    gap_opens += 1;
                    in_gap_ix = true;
                }
                in_gap_iy = false;

                // Ix can come from M (gap open) or Ix (gap extend)
                // We need to check the previous row's M vs Ix
                if k + 1 < band_size && cur_i > 0 {
                    // Check if we came from M or Ix in previous row
                    // This is a simplification - we assume gap extend if Ix was valid
                    cur_i -= 1;
                    // Stay in Ix or go to M based on scores (simplified: check if M was better)
                    // For now, assume we came from M (gap open already counted)
                    current_state = TbState::M;
                } else {
                    break;
                }
            }
            TbState::Iy => {
                // We're in Iy state (gap in query)
                gap_letters += 1;
                if !in_gap_iy {
                    gap_opens += 1;
                    in_gap_iy = true;
                }
                in_gap_ix = false;

                // Iy can come from M (gap open) or Iy (gap extend)
                if k > 0 && cur_j > 0 {
                    cur_j -= 1;
                    // Simplified: assume we came from M
                    current_state = TbState::M;
                } else {
                    break;
                }
            }
        }
    }

    // cur_i and cur_j are now at the start of the alignment (0-based in DP terms)
    let q_start = cur_i;
    let q_end = best_i;
    let s_start = cur_j;
    let s_end = best_j;

    (
        q_start,
        q_end,
        s_start,
        s_end,
        best_score,
        matches,
        mismatches,
        gap_opens,
        gap_letters,
    )
}

/// Chain nearby HSPs on similar diagonals into longer alignments using cluster-then-extend.
///
/// This function:
/// 1. Groups hits by query-subject pair
/// 2. Clusters nearby HSPs that are on similar diagonals and close in position
/// 3. For clusters with multiple HSPs, re-runs gapped alignment across the merged region
/// 4. Removes redundant overlapping hits
fn chain_and_filter_hsps(
    mut hits: Vec<Hit>,
    sequences: &FxHashMap<(String, String), (Vec<u8>, Vec<u8>)>,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    db_len_total: usize,
    db_num_seqs: usize,
    params: &KarlinParams,
    _use_dp: bool,
    verbose: bool,
) -> Vec<Hit> {
    if hits.is_empty() {
        return hits;
    }

    let total_start = std::time::Instant::now();

    // Parameters for chaining - adaptive based on HSP quality
    // Key insight: use bit_score/length as a proxy for identity
    // High-quality HSPs (high identity) can be merged across larger gaps
    // Low-quality HSPs (low identity) should use strict gap limits
    const MAX_GAP_STRICT: usize = 50; // For low-quality HSPs
    const MAX_GAP_PERMISSIVE: usize = 5000; // For high-quality HSPs
    const MAX_DIAG_DRIFT_STRICT: isize = 30; // For low-quality HSPs
    const MAX_DIAG_DRIFT_PERMISSIVE: isize = 500; // For high-quality HSPs
    const DIAG_BIN_SIZE: isize = 25; // Bin size for diagonal bucketing optimization

    // Threshold for "high quality" HSP: bit_score / length
    // For blastn with reward=2, penalty=-3:
    //   - 95% identity: expected score ≈ 0.95*2 + 0.05*(-3) = 1.75 per base
    //   - 74% identity: expected score ≈ 0.74*2 + 0.26*(-3) = 0.70 per base
    // Bit score = (raw_score * lambda - ln(K)) / ln(2)
    // For high identity (~95%), bit_score/length ≈ 1.5-1.8
    // For low identity (~74%), bit_score/length ≈ 0.5-0.8
    const HIGH_QUALITY_THRESHOLD: f64 = 1.2; // bit_score/length threshold

    // Group hits by query-subject pair
    let grouping_start = std::time::Instant::now();
    let mut groups: FxHashMap<(String, String), Vec<Hit>> = FxHashMap::default();
    for hit in hits.drain(..) {
        let key = (hit.query_id.clone(), hit.subject_id.clone());
        groups.entry(key).or_default().push(hit);
    }
    let grouping_time = grouping_start.elapsed();

    let mut result_hits = Vec::new();
    let mut total_clustering_time = std::time::Duration::ZERO;
    let mut total_align_time = std::time::Duration::ZERO;
    let mut align_calls = 0usize;
    let mut max_region_len = 0usize;
    let mut total_clusters = 0usize;
    let mut max_cluster_size = 0usize;

    if verbose {
        eprintln!(
            "[INFO] chain_and_filter_hsps: {} groups, grouping took {:.3}s",
            groups.len(),
            grouping_time.as_secs_f64()
        );
    }

    for ((query_id, subject_id), mut group_hits) in groups {
        if group_hits.is_empty() {
            continue;
        }

        let group_size = group_hits.len();
        let clustering_start = std::time::Instant::now();

        // Sort by query start position
        group_hits.sort_by_key(|h| h.q_start);

        // Optimized clustering using diagonal bins and active cluster tracking
        // Key insight: hits are sorted by q_start, so clusters whose last hit's q_end
        // is more than MAX_GAP behind the current hit's q_start can never accept new hits

        // Structure: diag_bin -> list of (cluster_idx, max_q_end)
        // We track cluster quality separately using sum_score/span
        let mut diag_bins: FxHashMap<isize, Vec<(usize, usize)>> = FxHashMap::default();
        let mut clusters: Vec<Vec<Hit>> = Vec::new();
        // Track cluster statistics for quality assessment: (sum_bit_score, q_min, q_max, current_bin)
        let mut cluster_stats: Vec<(f64, usize, usize, isize)> = Vec::new();

        for hit in group_hits {
            let hit_diag = hit.s_start as isize - hit.q_start as isize;
            let hit_diag_bin = hit_diag / DIAG_BIN_SIZE;

            // Determine if this hit is high-quality based on bit_score / length
            let hit_quality = hit.bit_score / (hit.length as f64);
            let hit_is_high_quality = hit_quality >= HIGH_QUALITY_THRESHOLD;

            // Use permissive search range for high-quality HSPs
            let search_max_diag_drift = if hit_is_high_quality {
                MAX_DIAG_DRIFT_PERMISSIVE
            } else {
                MAX_DIAG_DRIFT_STRICT
            };

            // Only search nearby diagonal bins (within search_max_diag_drift)
            let bin_range = (search_max_diag_drift / DIAG_BIN_SIZE) + 1;
            let mut cluster_idx: Option<usize> = None;

            'bin_search: for bin_offset in -bin_range..=bin_range {
                let search_bin = hit_diag_bin + bin_offset;
                if let Some(bin_clusters) = diag_bins.get(&search_bin) {
                    for &(idx, max_q_end) in bin_clusters.iter().rev() {
                        // Calculate cluster quality: average score density = sum_score / span
                        let (sum_score, q_min, q_max, _) = cluster_stats[idx];
                        let cluster_span = (q_max - q_min + 1) as f64;
                        let cluster_density = sum_score / cluster_span;
                        let cluster_is_high_quality = cluster_density >= HIGH_QUALITY_THRESHOLD;

                        // Both the hit and the cluster must be high-quality to use permissive parameters
                        let both_high_quality = hit_is_high_quality && cluster_is_high_quality;
                        let effective_max_gap = if both_high_quality {
                            MAX_GAP_PERMISSIVE
                        } else {
                            MAX_GAP_STRICT
                        };
                        let effective_max_diag = if both_high_quality {
                            MAX_DIAG_DRIFT_PERMISSIVE
                        } else {
                            MAX_DIAG_DRIFT_STRICT
                        };

                        // Skip clusters that are too far behind (can never match)
                        if hit.q_start > max_q_end + effective_max_gap {
                            continue;
                        }

                        let cluster = &clusters[idx];
                        if let Some(last) = cluster.last() {
                            let last_diag = last.s_end as isize - last.q_end as isize;
                            let diag_drift = (hit_diag - last_diag).abs();

                            let q_distance = hit.q_start as isize - last.q_end as isize;
                            let s_distance = hit.s_start as isize - last.s_end as isize;

                            let hit_q_len = (hit.q_end - hit.q_start) as isize;
                            let last_q_len = (last.q_end - last.q_start) as isize;
                            let max_overlap = -(hit_q_len.min(last_q_len) / 2);

                            let q_ok = q_distance >= max_overlap
                                && q_distance <= effective_max_gap as isize;
                            let s_ok = s_distance >= max_overlap
                                && s_distance <= effective_max_gap as isize;

                            if q_ok && s_ok && diag_drift <= effective_max_diag as isize {
                                cluster_idx = Some(idx);
                                break 'bin_search;
                            }
                        }
                    }
                }
            }

            if let Some(idx) = cluster_idx {
                clusters[idx].push(hit.clone());

                // Update cluster stats
                let (sum_score, q_min, q_max, old_bin) = cluster_stats[idx];
                let new_sum_score = sum_score + hit.bit_score;
                let new_q_min = q_min.min(hit.q_start);
                let new_q_max = q_max.max(hit.q_end);
                let new_diag = hit.s_end as isize - hit.q_end as isize;
                let new_bin = new_diag / DIAG_BIN_SIZE;
                cluster_stats[idx] = (new_sum_score, new_q_min, new_q_max, new_bin);

                // Fix: properly re-index cluster when diagonal bin changes
                if new_bin != old_bin {
                    // Remove from old bin
                    if let Some(old_bin_clusters) = diag_bins.get_mut(&old_bin) {
                        old_bin_clusters.retain(|(i, _)| *i != idx);
                    }
                    // Add to new bin
                    diag_bins.entry(new_bin).or_default().push((idx, new_q_max));
                } else {
                    // Just update max_q_end in the same bin
                    if let Some(bin_clusters) = diag_bins.get_mut(&new_bin) {
                        if let Some(entry) = bin_clusters.iter_mut().find(|(i, _)| *i == idx) {
                            entry.1 = new_q_max;
                        }
                    }
                }
            } else {
                let new_idx = clusters.len();
                clusters.push(vec![hit.clone()]);
                cluster_stats.push((hit.bit_score, hit.q_start, hit.q_end, hit_diag_bin));
                diag_bins
                    .entry(hit_diag_bin)
                    .or_default()
                    .push((new_idx, hit.q_end));
            }
        }

        let clustering_time = clustering_start.elapsed();
        total_clustering_time += clustering_time;
        total_clusters += clusters.len();
        max_cluster_size =
            max_cluster_size.max(clusters.iter().map(|c| c.len()).max().unwrap_or(0));

        if verbose && group_size > 1000 {
            eprintln!(
                "[INFO] Group {}/{}: {} hits -> {} clusters in {:.3}s",
                query_id,
                subject_id,
                group_size,
                clusters.len(),
                clustering_time.as_secs_f64()
            );
        }

        // PERFORMANCE OPTIMIZATION: Only run align_region on the best clusters
        // NCBI BLAST is fast because it only does gapped extension on a small subset of candidates
        // We score clusters by their total bit_score and only align the top ones
        const MAX_ALIGN_REGION_CALLS: usize = 20; // Limit expensive align_region calls
        const MIN_CLUSTER_SCORE_FOR_ALIGN: f64 = 500.0; // Minimum score to consider for align_region
        
        // Maximum region size for greedy re-alignment to prevent over-extension
        // NCBI BLAST typically produces alignments up to ~50kbp for blastn task
        // Larger regions should use the best HSP instead of re-alignment to match NCBI behavior
        const MAX_REGION_SIZE_FOR_ALIGN: usize = 50000;

        // Collect multi-HSP clusters with their scores for prioritization
        let mut multi_hsp_clusters_with_scores: Vec<(usize, f64, Vec<Hit>)> = Vec::new();
        let mut single_hsp_hits: Vec<Hit> = Vec::new();

        for (idx, cluster) in clusters.into_iter().enumerate() {
            if cluster.is_empty() {
                continue;
            }

            if cluster.len() == 1 {
                // Single HSP, no re-alignment needed
                single_hsp_hits.push(cluster.into_iter().next().unwrap());
            } else {
                // Multi-HSP cluster - calculate total score
                let (sum_score, _, _, _) = cluster_stats[idx];
                multi_hsp_clusters_with_scores.push((idx, sum_score, cluster));
            }
        }

        // Sort multi-HSP clusters by score (descending) and take top N
        multi_hsp_clusters_with_scores
            .sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        // Process each cluster with progress logging
        let num_multi_clusters = multi_hsp_clusters_with_scores.len();
        let mut processed_clusters = 0usize;
        let mut aligned_clusters = 0usize;
        let cluster_process_start = std::time::Instant::now();

        for (_idx, sum_score, cluster) in multi_hsp_clusters_with_scores {
            processed_clusters += 1;

            // Progress logging every 100 clusters for multi-HSP
            if verbose && processed_clusters % 100 == 0 {
                eprintln!("[INFO] Processed {}/{} multi-HSP clusters, {} align_region calls, {:.3}s elapsed",
                          processed_clusters, num_multi_clusters, align_calls,
                          cluster_process_start.elapsed().as_secs_f64());
            }

            // Calculate cluster region size to check if it's too large for re-alignment
            let cluster_q_min = cluster.iter().map(|h| h.q_start).min().unwrap_or(0);
            let cluster_q_max = cluster.iter().map(|h| h.q_end).max().unwrap_or(0);
            let cluster_s_min = cluster.iter().map(|h| h.s_start).min().unwrap_or(0);
            let cluster_s_max = cluster.iter().map(|h| h.s_end).max().unwrap_or(0);
            let cluster_region_size = (cluster_q_max.saturating_sub(cluster_q_min))
                .max(cluster_s_max.saturating_sub(cluster_s_min));

            // Only run align_region for top clusters with high scores AND reasonable region size
            // For lower-scoring clusters or very large regions, just use the best single HSP
            // This prevents over-extension that produces alignments much longer than NCBI BLAST
            let should_align = aligned_clusters < MAX_ALIGN_REGION_CALLS
                && sum_score >= MIN_CLUSTER_SCORE_FOR_ALIGN
                && cluster_region_size <= MAX_REGION_SIZE_FOR_ALIGN;

            if should_align {
                aligned_clusters += 1;

                // NCBI BLAST optimization: Instead of aligning the entire merged region
                // (which can be 657kbp and takes O(n * band_size) time), we use greedy
                // extension from the best HSP's center point. This is O(alignment_length)
                // and much faster for high-identity sequences.
                //
                // Reference: NCBI BLAST's BLAST_GappedAlignmentWithTraceback in blast_gapalign.c
                // extends LEFT and RIGHT from the HSP center using ALIGN_EX or greedy alignment.

                // Get sequences for this query-subject pair
                let key = (query_id.clone(), subject_id.clone());
                if let Some((q_seq, s_seq)) = sequences.get(&key) {
                    // Find the best HSP in the cluster to use as seed for extension
                    let best_hsp = cluster
                        .iter()
                        .max_by(|a, b| {
                            a.bit_score
                                .partial_cmp(&b.bit_score)
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap();

                    // Use the center of the best HSP as the seed point for extension
                    // Convert 1-based coordinates to 0-based
                    let seed_q_start = best_hsp.q_start.saturating_sub(1);
                    let seed_q_end = best_hsp.q_end;
                    let seed_s_start = best_hsp.s_start.saturating_sub(1);
                    let seed_s_end = best_hsp.s_end;

                    // Calculate seed center and length
                    let seed_q_center = (seed_q_start + seed_q_end) / 2;
                    let seed_s_center = (seed_s_start + seed_s_end) / 2;
                    let seed_len = 1; // Start from a single point and extend

                    let align_start = std::time::Instant::now();

                    // ALWAYS use greedy extension in chaining phase (NCBI BLAST approach)
                    // Greedy is better for chaining because:
                    // 1. It doesn't have the MAX_WINDOW_SIZE limit that caps DP at ~10kbp
                    // 2. It's faster and handles high-identity sequences well
                    // 3. NCBI BLAST uses greedy for HSP extension even in blastn task
                    // The use_dp parameter only affects the initial seed extension, not chaining
                    let (
                        final_q_start_0,
                        final_q_end_0,
                        final_s_start_0,
                        final_s_end_0,
                        score,
                        matches,
                        mismatches,
                        gap_opens,
                        gap_letters,
                        _dp_cells,
                    ) = extend_gapped_heuristic(
                        q_seq,
                        s_seq,
                        seed_q_center,
                        seed_s_center,
                        seed_len,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        X_DROP_GAPPED_FINAL,
                        false, // Always use greedy for chaining re-alignment
                    );

                    total_align_time += align_start.elapsed();
                    align_calls += 1;

                    // Track region length for logging
                    let region_len =
                        (final_q_end_0 - final_q_start_0).max(final_s_end_0 - final_s_start_0);
                    max_region_len = max_region_len.max(region_len);

                    // Skip if no valid alignment found
                    if score <= 0 {
                        // Fall back to best HSP
                        result_hits.push(best_hsp.clone());
                        continue;
                    }

                    // Calculate statistics using BLAST's definition:
                    // alignment_length = matches + mismatches + gap_letters (total aligned columns)
                    let aln_len = matches + mismatches + gap_letters;

                    let identity = if aln_len > 0 {
                        ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                    } else {
                        0.0
                    };

                    let (bit_score, e_value) =
                        calculate_evalue(score, q_seq.len(), db_len_total, db_num_seqs, params);

                    // Convert 0-based coordinates to 1-based for output
                    let final_q_start = final_q_start_0 + 1;
                    let final_q_end = final_q_end_0;
                    let final_s_start = final_s_start_0 + 1;
                    let final_s_end = final_s_end_0;

                    result_hits.push(Hit {
                        query_id: query_id.clone(),
                        subject_id: subject_id.clone(),
                        identity,
                        length: aln_len,
                        mismatch: mismatches,
                        gapopen: gap_opens,
                        q_start: final_q_start,
                        q_end: final_q_end,
                        s_start: final_s_start,
                        s_end: final_s_end,
                        e_value,
                        bit_score,
                    });
                } else {
                    // Fallback: just use the best HSP (shouldn't happen)
                    let best_hsp = cluster
                        .into_iter()
                        .max_by(|a, b| {
                            a.bit_score
                                .partial_cmp(&b.bit_score)
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap();
                    result_hits.push(best_hsp);
                }
            } else {
                // PERFORMANCE: For lower-scoring clusters, just use the best single HSP
                // This avoids expensive align_region calls for clusters that won't produce top hits
                let best_hsp = cluster
                    .into_iter()
                    .max_by(|a, b| {
                        a.bit_score
                            .partial_cmp(&b.bit_score)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap();
                result_hits.push(best_hsp);
            }
        }

        // Add single-HSP hits
        result_hits.extend(single_hsp_hits);
    }

    // Log after cluster processing
    if verbose {
        eprintln!("[INFO] Cluster processing done: {} result_hits, {} align_region calls in {:.3}s, max_region={}bp",
                  result_hits.len(), align_calls, total_align_time.as_secs_f64(), max_region_len);
    }

    // Final pass: remove redundant overlapping hits using spatial binning for O(n) instead of O(n²)
    // IMPORTANT: For self-comparison, hits on different diagonals represent different biological
    // alignments (e.g., repeats, duplications) and should NOT be filtered out even if they
    // overlap on query coordinates. We only filter hits that overlap on BOTH query AND subject
    // AND are on similar diagonals (same alignment region).
    let filter_start = std::time::Instant::now();
    let result_hits_count = result_hits.len();

    result_hits.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(Ordering::Equal)
    });

    // Use spatial binning to avoid O(n²) comparisons
    // Bin size of 1000bp - hits can only overlap if they share a bin
    const FILTER_BIN_SIZE: usize = 1000;
    const MAX_DIAG_DIFF_FOR_SAME_ALIGNMENT: isize = 100; // Hits must be on similar diagonals to be considered redundant
    let mut q_bins: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    let mut final_hits: Vec<Hit> = Vec::new();

    for hit in result_hits {
        // Calculate which q bins this hit spans
        let q_bin_start = hit.q_start / FILTER_BIN_SIZE;
        let q_bin_end = hit.q_end / FILTER_BIN_SIZE;

        // Calculate diagonal for this hit (s_start - q_start)
        let hit_diag = hit.s_start as isize - hit.q_start as isize;

        // Only check hits in overlapping bins
        let mut dominated = false;
        for bin in q_bin_start..=q_bin_end {
            if let Some(kept_indices) = q_bins.get(&bin) {
                for &kept_idx in kept_indices {
                    let kept = &final_hits[kept_idx];

                    // First check if hits are on similar diagonals
                    // Hits on different diagonals represent different biological alignments
                    let kept_diag = kept.s_start as isize - kept.q_start as isize;
                    let diag_diff = (hit_diag - kept_diag).abs();
                    if diag_diff > MAX_DIAG_DIFF_FOR_SAME_ALIGNMENT {
                        continue; // Different diagonal = different alignment, don't filter
                    }

                    // Check q overlap
                    let q_overlap_start = hit.q_start.max(kept.q_start);
                    let q_overlap_end = hit.q_end.min(kept.q_end);
                    let q_overlap = if q_overlap_end > q_overlap_start {
                        q_overlap_end - q_overlap_start
                    } else {
                        0
                    };
                    let q_hit_len = hit.q_end - hit.q_start;
                    let q_overlap_frac = if q_hit_len > 0 {
                        q_overlap as f64 / q_hit_len as f64
                    } else {
                        0.0
                    };

                    if q_overlap_frac < 0.5 {
                        continue; // No significant q overlap, skip s check
                    }

                    // Check s overlap
                    let s_overlap_start = hit.s_start.max(kept.s_start);
                    let s_overlap_end = hit.s_end.min(kept.s_end);
                    let s_overlap = if s_overlap_end > s_overlap_start {
                        s_overlap_end - s_overlap_start
                    } else {
                        0
                    };
                    let s_hit_len = hit.s_end - hit.s_start;
                    let s_overlap_frac = if s_hit_len > 0 {
                        s_overlap as f64 / s_hit_len as f64
                    } else {
                        0.0
                    };

                    if s_overlap_frac >= 0.5 {
                        dominated = true;
                        break;
                    }
                }
                if dominated {
                    break;
                }
            }
        }

        if !dominated {
            let new_idx = final_hits.len();
            // Register this hit in all bins it spans
            for bin in q_bin_start..=q_bin_end {
                q_bins.entry(bin).or_default().push(new_idx);
            }
            final_hits.push(hit);
        }
    }
    let filter_time = filter_start.elapsed();

    if verbose {
        eprintln!("[INFO] chain_and_filter_hsps summary:");
        eprintln!(
            "[INFO]   Total time: {:.3}s",
            total_start.elapsed().as_secs_f64()
        );
        eprintln!(
            "[INFO]   Clustering: {:.3}s ({} clusters, max size {})",
            total_clustering_time.as_secs_f64(),
            total_clusters,
            max_cluster_size
        );
        eprintln!(
            "[INFO]   Alignment: {:.3}s ({} calls, max region {}bp)",
            total_align_time.as_secs_f64(),
            align_calls,
            max_region_len
        );
        eprintln!(
            "[INFO]   Filtering: {:.3}s ({} -> {} hits)",
            filter_time.as_secs_f64(),
            result_hits_count,
            final_hits.len()
        );
    }

    final_hits
}

fn extend_hit_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_pos: usize,
    s_pos: usize,
    reward: i32,
    penalty: i32,
) -> (usize, usize, usize, usize, i32) {
    let mut current_score = 0;
    let mut max_score = 0;

    // Left
    let mut best_i = 0;
    let mut i = 0;
    let max_left = q_pos.min(s_pos);

    while i <= max_left {
        let q = unsafe { *q_seq.get_unchecked(q_pos - i) };
        let s = unsafe { *s_seq.get_unchecked(s_pos - i) };
        let sc = if q == s { reward } else { penalty };
        current_score += sc;

        if current_score > max_score {
            max_score = current_score;
            best_i = i;
        } else if (max_score - current_score) > X_DROP_UNGAPPED {
            break;
        }
        i += 1;
    }

    // Right
    let mut current_score_r = max_score;
    let mut max_score_total = max_score;
    let mut best_j = 0;
    let mut j = 1;

    while (q_pos + j) < q_seq.len() && (s_pos + j) < s_seq.len() {
        let q = unsafe { *q_seq.get_unchecked(q_pos + j) };
        let s = unsafe { *s_seq.get_unchecked(s_pos + j) };
        let sc = if q == s { reward } else { penalty };
        current_score_r += sc;

        if current_score_r > max_score_total {
            max_score_total = current_score_r;
            best_j = j;
        } else if (max_score_total - current_score_r) > X_DROP_UNGAPPED {
            break;
        }
        j += 1;
    }

    (
        q_pos - best_i,
        q_pos + best_j + 1,
        s_pos - best_i,
        s_pos + best_j + 1,
        max_score_total,
    )
}

pub fn run(args: BlastnArgs) -> Result<()> {
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    // Determine effective word size based on task
    // If user specified default word_size (28) and task is not megablast, use task-appropriate defaults
    let effective_word_size = match args.task.as_str() {
        "megablast" => args.word_size,
        "blastn" => {
            if args.word_size == 28 {
                11
            } else {
                args.word_size
            }
        }
        "dc-megablast" => {
            if args.word_size == 28 {
                11
            } else {
                args.word_size
            }
        }
        _ => args.word_size,
    };

    // Determine effective scoring parameters based on task
    // NCBI BLAST uses different defaults for different tasks:
    // - megablast: reward=1, penalty=-2, gapopen=0, gapextend=0 (linear gap)
    // - blastn: reward=2, penalty=-3, gapopen=5, gapextend=2
    // - dc-megablast: reward=2, penalty=-3, gapopen=5, gapextend=2
    // - blastn-short: reward=1, penalty=-3, gapopen=5, gapextend=2
    // Note: gap costs are stored as negative values internally for DP scoring
    let (reward, penalty, gap_open, gap_extend) = match args.task.as_str() {
        "megablast" => {
            // Use user-specified values or megablast defaults
            let r = if args.reward == 1 { 1 } else { args.reward };
            let p = if args.penalty == -2 { -2 } else { args.penalty };
            let go = if args.gap_open == 0 { 0 } else { args.gap_open };
            let ge = if args.gap_extend == 0 {
                0
            } else {
                args.gap_extend
            }; // Linear gap (0 triggers non-affine greedy)
            (r, p, go, ge)
        }
        "blastn" | "dc-megablast" => {
            // Use user-specified values or blastn defaults (reward=2, penalty=-3)
            let r = if args.reward == 1 { 2 } else { args.reward };
            let p = if args.penalty == -2 { -3 } else { args.penalty };
            let go = if args.gap_open == 0 {
                -5
            } else {
                args.gap_open
            };
            let ge = if args.gap_extend == 0 {
                -2
            } else {
                args.gap_extend
            };
            (r, p, go, ge)
        }
        "blastn-short" => {
            // blastn-short: reward=1, penalty=-3, gapopen=5, gapextend=2
            let r = if args.reward == 1 { 1 } else { args.reward };
            let p = if args.penalty == -2 { -3 } else { args.penalty };
            let go = if args.gap_open == 0 {
                -5
            } else {
                args.gap_open
            };
            let ge = if args.gap_extend == 0 {
                -2
            } else {
                args.gap_extend
            };
            (r, p, go, ge)
        }
        _ => (args.reward, args.penalty, args.gap_open, args.gap_extend),
    };

    // Task-specific ungapped score threshold for triggering gapped extension
    // Higher threshold for blastn reduces the number of gapped extensions significantly
    // This is critical for performance on self-comparison with many off-diagonal matches
    let min_ungapped_score = match args.task.as_str() {
        "megablast" => MIN_UNGAPPED_SCORE_MEGABLAST,
        _ => MIN_UNGAPPED_SCORE_BLASTN,
    };

    // NCBI BLAST uses different extension algorithms based on task:
    // - megablast: greedy alignment (eGreedyScoreOnly) - fast, good for high-identity sequences
    // - blastn: DP-based alignment (eDynProgScoreOnly) - slower, but handles divergent sequences better
    // This follows NCBI BLAST's approach for better alignment of divergent sequences
    let use_dp = match args.task.as_str() {
        "megablast" => false, // Use greedy for megablast (high-identity sequences)
        _ => true,            // Use DP for blastn, dc-megablast, blastn-short (divergent sequences)
    };

    eprintln!("Reading query & subject...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries
        .iter()
        .map(|r| {
            r.id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string()
        })
        .collect();

    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();

    if queries.is_empty() || subjects.is_empty() {
        return Ok(());
    }

    let db_len_total: usize = subjects.iter().map(|r| r.seq().len()).sum();
    let db_num_seqs: usize = subjects.len();

    let scoring_spec = NuclScoringSpec {
        reward,
        penalty,
        gap_open,
        gap_extend,
    };
    let params = lookup_nucl_params(&scoring_spec);

    // Choose between direct address table (O(1) lookup) and HashMap based on word size
    // Direct address table is much faster but requires 4^word_size memory
    let use_direct_lookup = effective_word_size <= MAX_DIRECT_LOOKUP_WORD_SIZE;

    // Apply DUST filter to mask low-complexity regions in query sequences
    let query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        eprintln!(
            "Applying DUST filter (level={}, window={}, linker={})...",
            args.dust_level, args.dust_window, args.dust_linker
        );
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        let masks: Vec<Vec<MaskedInterval>> = queries
            .iter()
            .map(|record| masker.mask_sequence(record.seq()))
            .collect();
        
        // Report DUST masking statistics
        let total_masked: usize = masks.iter().map(|m| m.iter().map(|i| i.end - i.start).sum::<usize>()).sum();
        let total_bases: usize = queries.iter().map(|r| r.seq().len()).sum();
        if total_bases > 0 {
            eprintln!(
                "DUST masked {} bases ({:.2}%) across {} sequences",
                total_masked,
                100.0 * total_masked as f64 / total_bases as f64,
                queries.len()
            );
        }
        masks
    } else {
        vec![Vec::new(); queries.len()]
    };

    eprintln!(
        "Building lookup (Task: {}, Word: {}, Direct: {}, DUST: {}, PV: {})...",
        args.task, effective_word_size, use_direct_lookup, args.dust, use_direct_lookup
    );

    // Build the appropriate lookup table
    // Phase 2: Use PvDirectLookup with Presence-Vector for fast filtering (word_size <= 13)
    // For word_size > 13, use hash-based lookup
    let pv_direct_lookup: Option<PvDirectLookup> = if use_direct_lookup {
        Some(build_pv_direct_lookup(&queries, effective_word_size, &query_masks))
    } else {
        None
    };
    let hash_lookup: Option<KmerLookup> = if !use_direct_lookup {
        Some(build_lookup(&queries, effective_word_size, &query_masks))
    } else {
        None
    };

    if args.verbose {
        eprintln!("Searching...");
    }

    let bar = ProgressBar::new(subjects.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    // Channel for sending hits and sequence data
    // Use Option to signal completion: None means "all subjects processed"
    let (tx, rx) = channel::<Option<(Vec<Hit>, Vec<(String, String, Vec<u8>, Vec<u8>)>)>>();
    let out_path = args.out.clone();
    // Note: reward, penalty, gap_open, gap_extend are already defined above with task-adjusted values
    let params_clone = params.clone();

    // Keep a sender for the main thread to send the completion signal
    let tx_main = tx.clone();
    let verbose = args.verbose;

    let writer_handle = std::thread::spawn(move || -> Result<()> {
        if verbose {
            eprintln!("[INFO] Writer thread started, waiting for hits...");
        }
        let mut all_hits = Vec::new();
        let mut all_sequences: FxHashMap<(String, String), (Vec<u8>, Vec<u8>)> =
            FxHashMap::default();
        let mut messages_received = 0usize;

        while let Ok(msg) = rx.recv() {
            match msg {
                Some((hits, seq_data)) => {
                    messages_received += 1;
                    if verbose && (messages_received == 1 || messages_received % 100 == 0) {
                        eprintln!(
                            "[INFO] Received message #{}, {} hits so far",
                            messages_received,
                            all_hits.len() + hits.len()
                        );
                    }
                    all_hits.extend(hits);
                    for (q_id, s_id, q_seq, s_seq) in seq_data {
                        all_sequences.entry((q_id, s_id)).or_insert((q_seq, s_seq));
                    }
                }
                None => {
                    // Completion signal received
                    if verbose {
                        eprintln!(
                            "[INFO] Completion signal received after {} messages",
                            messages_received
                        );
                    }
                    break;
                }
            }
        }

        if verbose {
            eprintln!(
                "[INFO] Received {} raw hits total, starting post-processing...",
                all_hits.len()
            );
        }
        let chain_start = std::time::Instant::now();

        // Chain nearby HSPs into longer alignments using cluster-then-extend
        let filtered_hits = chain_and_filter_hsps(
            all_hits,
            &all_sequences,
            reward,
            penalty,
            gap_open,
            gap_extend,
            db_len_total,
            db_num_seqs,
            &params_clone,
            use_dp,
            verbose,
        );

        if verbose {
            eprintln!(
                "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
                chain_start.elapsed().as_secs_f64(),
                filtered_hits.len()
            );
        }
        let write_start = std::time::Instant::now();

        write_output(&filtered_hits, out_path.as_ref())?;

        if verbose {
            eprintln!(
                "[INFO] Output written in {:.2}s",
                write_start.elapsed().as_secs_f64()
            );
        }
        Ok(())
    });

    // 修正: s_idx -> _s_idx (未使用抑制)
    // Debug mode: set BLEMIR_DEBUG=1 to enable, BLEMIR_DEBUG_WINDOW="q_start-q_end,s_start-s_end" to focus on a region
    let debug_mode = std::env::var("BLEMIR_DEBUG").is_ok();
    let debug_window: Option<(usize, usize, usize, usize)> =
        std::env::var("BLEMIR_DEBUG_WINDOW").ok().and_then(|s| {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() == 2 {
                let q_parts: Vec<usize> =
                    parts[0].split('-').filter_map(|x| x.parse().ok()).collect();
                let s_parts: Vec<usize> =
                    parts[1].split('-').filter_map(|x| x.parse().ok()).collect();
                if q_parts.len() == 2 && s_parts.len() == 2 {
                    return Some((q_parts[0], q_parts[1], s_parts[0], s_parts[1]));
                }
            }
            None
        });

    // Check if mask should be disabled for debugging
    let disable_mask = std::env::var("BLEMIR_DEBUG_NO_MASK").is_ok();

    if debug_mode {
        // Build marker to verify correct code is running
        eprintln!(
            "[DEBUG] BLEMIR build: 2024-12-24-v7 (adaptive banding with MAX_WINDOW_SIZE=50000)"
        );
        eprintln!(
            "[DEBUG] Task: {}, Scoring: reward={}, penalty={}, gap_open={}, gap_extend={}",
            args.task, reward, penalty, gap_open, gap_extend
        );
        if disable_mask {
            eprintln!("[DEBUG] Mask DISABLED for debugging");
        }
        if let Some((q_start, q_end, s_start, s_end)) = debug_window {
            eprintln!(
                "[DEBUG] Focusing on window: query {}-{}, subject {}-{}",
                q_start, q_end, s_start, s_end
            );
        }
    }

    // 修正: tx はここで所有権ごと渡す。for_each_with 終了時に自動でDropされるため、明示的な drop(tx) は不要。
    subjects
        .par_iter()
        .enumerate()
        .for_each_with(tx, |tx, (_s_idx, s_record)| {
            let mut local_hits: Vec<Hit> = Vec::new();
            let mut local_sequences: Vec<(String, String, Vec<u8>, Vec<u8>)> = Vec::new();
            let mut seen_pairs: std::collections::HashSet<(String, String)> =
                std::collections::HashSet::new();

            let s_seq = s_record.seq();
            let s_len = s_seq.len();
            let s_id = s_record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();

            if s_len < effective_word_size {
                return;
            }

            // Debug counters for this subject
            let mut dbg_total_s_positions = 0usize;
            let dbg_ambiguous_skipped = 0usize;
            let dbg_no_lookup_match = 0usize;
            let mut dbg_seeds_found = 0usize;
            let mut dbg_mask_skipped = 0usize;
            let mut dbg_ungapped_low = 0usize;
            let mut dbg_two_hit_failed = 0usize;
            let mut dbg_gapped_attempted = 0usize;
            let mut dbg_window_seeds = 0usize;

            // Mask to track already-extended regions on each diagonal
            // Phase 4: Use packed u64 key instead of (u32, isize) tuple for faster hashing
            let mut mask: FxHashMap<u64, usize> = FxHashMap::default();
            // Two-hit tracking: stores the last seed position on each diagonal for each query
            // Phase 4: Use packed u64 key instead of (u32, isize) tuple for faster hashing
            let mut last_seed: FxHashMap<u64, usize> = FxHashMap::default();
            let safe_k = effective_word_size.min(31);

            // PERFORMANCE OPTIMIZATION: Rolling k-mer scanner with O(1) sliding window
            // Instead of packing the entire sequence first (which adds O(n) overhead),
            // we compute k-mers on-the-fly using a rolling window approach.
            // This achieves O(1) per-position k-mer extraction without allocation overhead.
            //
            // Encoding: A=0, C=1, G=2, T=3 (same as NCBI BLAST's ncbi2na)
            // The mask ensures we only keep the rightmost 2*k bits.
            let kmer_mask: u64 = (1u64 << (2 * safe_k)) - 1;

            // Lookup table for ASCII to 2-bit encoding (0xFF = invalid/ambiguous)
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

            // Rolling k-mer state
            let mut current_kmer: u64 = 0;
            let mut valid_bases: usize = 0; // Count of consecutive valid bases in current window

            // Scan through the subject sequence with rolling k-mer
            for s_pos in 0..s_len {
                let base = s_seq[s_pos];
                let code = ENCODE_LUT[base as usize];

                if code == 0xFF {
                    // Ambiguous base - reset the rolling window
                    valid_bases = 0;
                    current_kmer = 0;
                    continue;
                }

                // Shift in the new base
                current_kmer = ((current_kmer << 2) | (code as u64)) & kmer_mask;
                valid_bases += 1;

                // Only process if we have a complete k-mer
                if valid_bases < safe_k {
                    continue;
                }

                // Calculate the starting position of this k-mer
                let kmer_start = s_pos + 1 - safe_k;
                dbg_total_s_positions += 1;

                // Phase 2: Use PV-based direct lookup (O(1) with fast PV filtering) for word_size <= 13
                // For word_size > 13, use hash-based lookup
                let matches_slice: &[(u32, u32)] = if use_direct_lookup {
                    // Use PV for fast filtering before accessing the lookup table
                    pv_direct_lookup.as_ref().map(|pv_dl| pv_dl.get_hits_checked(current_kmer)).unwrap_or(&[])
                } else {
                    // Use hash-based lookup for larger word sizes
                    hash_lookup.as_ref().and_then(|hl| hl.get(&current_kmer).map(|v| v.as_slice())).unwrap_or(&[])
                };

                for &(q_idx, q_pos) in matches_slice {
                    dbg_seeds_found += 1;
                    let diag = kmer_start as isize - q_pos as isize;

                    // Check if this seed is in the debug window
                    let in_window = if let Some((q_start, q_end, s_start, s_end)) = debug_window {
                        (q_pos as usize) >= q_start && (q_pos as usize) <= q_end &&
                        kmer_start >= s_start && kmer_start <= s_end
                    } else {
                        false
                    };

                    if in_window {
                        dbg_window_seeds += 1;
                    }

                    // Check if this region was already extended (skip if mask is disabled for debugging)
                    // Phase 4: Use packed key for faster HashMap lookup
                    let diag_key = pack_diag_key(q_idx, diag);
                    if !disable_mask {
                        if let Some(&last_s_end) = mask.get(&diag_key) {
                            if kmer_start < last_s_end {
                                dbg_mask_skipped += 1;
                                if in_window && debug_mode {
                                    eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED by mask (last_s_end={})", q_pos, kmer_start, last_s_end);
                                }
                                continue;
                            }
                        }
                    }

                    // PERFORMANCE OPTIMIZATION: Apply two-hit filter BEFORE ungapped extension
                    // This is the key optimization that makes NCBI BLAST fast
                    // Most seeds don't have a nearby seed on the same diagonal, so we skip them early
                    // Phase 4: Use packed key (diag_key computed above) for faster HashMap lookup
                    let trigger_extension =
                        if let Some(&prev_s_pos) = last_seed.get(&diag_key) {
                            // Check if the current seed is within the two-hit window
                            kmer_start.saturating_sub(prev_s_pos) <= TWO_HIT_WINDOW
                        } else {
                            false
                        };

                    // Update the last seed position for this diagonal
                    // Phase 4: Use packed key for faster HashMap insert
                    last_seed.insert(diag_key, kmer_start);

                    // Skip ungapped extension if two-hit requirement is not met
                    // This is the key optimization - we skip most seeds without doing any extension
                    if !trigger_extension {
                        dbg_two_hit_failed += 1;
                        if in_window && debug_mode {
                            eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: two-hit not met", q_pos, kmer_start);
                        }
                        continue;
                    }

                    let q_record = &queries[q_idx as usize];

                    // Now do ungapped extension (only for seeds that passed two-hit filter)
                    let (qs, qe, ss, _se, ungapped_score) = extend_hit_ungapped(
                        q_record.seq(),
                        s_seq,
                        q_pos as usize,
                        kmer_start,
                        reward,
                        penalty,
                    );

                    // Skip if ungapped score is too low
                    // Use task-specific threshold: higher for blastn to reduce extension count
                    if ungapped_score < min_ungapped_score {
                        dbg_ungapped_low += 1;
                        if in_window && debug_mode {
                            eprintln!("[DEBUG WINDOW] Seed at q={}, s={} SKIPPED: ungapped_score={} < {}", q_pos, kmer_start, ungapped_score, min_ungapped_score);
                        }
                        continue;
                    }

                    dbg_gapped_attempted += 1;

                    if in_window && debug_mode {
                        eprintln!("[DEBUG WINDOW] Seed at q={}, s={} -> GAPPED EXTENSION (ungapped_score={}, seed_len={})", q_pos, kmer_start, ungapped_score, qe - qs);
                    }

                    // Gapped extension with NCBI-style high X-drop for longer alignments
                    // Using X_DROP_GAPPED_FINAL directly to allow alignments to push through
                    // low-similarity regions (NCBI BLAST uses 100 for nucleotide)
                    let (
                        final_qs,
                        final_qe,
                        final_ss,
                        final_se,
                        score,
                        matches,
                        mismatches,
                        gaps,
                        gap_letters,
                        dp_cells,
                    ) = extend_gapped_heuristic(
                        q_record.seq(),
                        s_seq,
                        qs,
                        ss,
                        qe - qs,
                        reward,
                        penalty,
                        gap_open,
                        gap_extend,
                        X_DROP_GAPPED_FINAL,
                        use_dp, // Use DP for blastn task, greedy for megablast
                    );

                    // Debug: show gapped extension results for window seeds
                    if in_window && debug_mode {
                        let aln_len = matches + mismatches + gap_letters;
                        let identity = if aln_len > 0 { 100.0 * matches as f64 / aln_len as f64 } else { 0.0 };
                        eprintln!(
                            "[DEBUG WINDOW] Gapped result: q={}-{}, s={}-{}, score={}, len={}, identity={:.1}%, gaps={}, dp_cells={}",
                            final_qs, final_qe, final_ss, final_se, score, aln_len, identity, gap_letters, dp_cells
                        );
                    }

                    // Suppress unused variable warning when not in debug mode
                    let _ = dp_cells;

                    // Update mask (unless disabled for debugging via BLEMIR_DEBUG_NO_MASK=1)
                    // Phase 4: Use packed key (diag_key computed earlier) for faster HashMap insert
                    if !disable_mask {
                        mask.insert(diag_key, final_se);
                    }

                    let (bit_score, eval) = calculate_evalue(
                        score,
                        q_record.seq().len(),
                        db_len_total,
                        db_num_seqs,
                        &params,
                    );

                    if eval <= args.evalue {
                        // Calculate alignment length using BLAST's definition:
                        // alignment_length = matches + mismatches + gap_letters (total aligned columns)
                        let aln_len = matches + mismatches + gap_letters;

                        // Identity is matches / alignment_length, capped at 100%
                        let identity = if aln_len > 0 {
                            ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                        } else {
                            0.0
                        };

                        let q_id = query_ids[q_idx as usize].clone();

                        // Store sequence data for cluster-then-extend chaining
                        let pair_key = (q_id.clone(), s_id.clone());
                        if !seen_pairs.contains(&pair_key) {
                            seen_pairs.insert(pair_key);
                            local_sequences.push((
                                q_id.clone(),
                                s_id.clone(),
                                q_record.seq().to_vec(),
                                s_seq.to_vec(),
                            ));
                        }

                        local_hits.push(Hit {
                            query_id: q_id,
                            subject_id: s_id.clone(),
                            identity,
                            length: aln_len,
                            mismatch: mismatches,
                            gapopen: gaps,
                            q_start: final_qs + 1,
                            q_end: final_qe,
                            s_start: final_ss + 1,
                            s_end: final_se,
                            e_value: eval,
                            bit_score,
                        });
                    }
                }
            }

            // Print debug summary for this subject
            if debug_mode {
                eprintln!(
                    "[DEBUG] Subject {}: positions={}, seeds_found={}, mask_skipped={}, ungapped_low={}, two_hit_failed={}, gapped_attempted={}, window_seeds={}, hits={}",
                    s_id, dbg_total_s_positions, dbg_seeds_found, dbg_mask_skipped, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, dbg_window_seeds, local_hits.len()
                );
            }

            // Suppress unused variable warnings when not in debug mode
            let _ = (dbg_total_s_positions, dbg_ambiguous_skipped, dbg_no_lookup_match, dbg_seeds_found, dbg_mask_skipped, dbg_ungapped_low, dbg_two_hit_failed, dbg_gapped_attempted, dbg_window_seeds);

            if !local_hits.is_empty() {
                tx.send(Some((local_hits, local_sequences))).unwrap();
            }
            bar.inc(1);
        });

    bar.finish();
    if args.verbose {
        eprintln!("[INFO] Parallel processing complete, sending completion signal...");
    }

    // Send completion signal to writer thread
    // This ensures the writer thread exits even if sender-drop semantics are delayed
    tx_main.send(None).unwrap();
    drop(tx_main); // Explicitly drop to ensure channel closes

    writer_handle.join().unwrap()?;
    Ok(())
}
