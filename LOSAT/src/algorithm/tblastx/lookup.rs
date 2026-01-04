//! Lookup table construction for TBLASTX - EXACT NCBI BLAST port
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c
//!            ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h
//!
//! NCBI structure:
//! - AaLookupBackboneCell: flat array, 16 bytes per cell
//!   - num_used: i32 (4 bytes)
//!   - payload: union { entries[3], overflow_cursor } (12 bytes)
//! - overflow: flat i32 array for cells with >3 hits
//! - pv: presence vector bitfield

use crate::utils::matrix::{BLASTAA_SIZE, blosum62_score};
use crate::stats::{KarlinParams, lookup_protein_params_ungapped};
use crate::config::ScoringMatrix;
use crate::stats::karlin_calc::{
    compute_aa_composition, compute_std_aa_composition, compute_score_freq_profile,
    compute_karlin_params_ungapped, apply_check_ideal,
};
use super::translation::QueryFrame;

pub const AA_HITS_PER_CELL: usize = 3;

// Presence Vector - NCBI style
const PV_ARRAY_BTS: usize = 6;
const PV_ARRAY_MASK: usize = 63;
const PV_BUCKET_BITS: usize = 64;

/// NCBI BLAST parameters for protein lookup indexing.
///
/// We use the NCBI-style bit-shift + mask indexing scheme from:
/// - `blast_lookup.h`: `ComputeTableIndex` / `ComputeTableIndexIncremental`
/// - `blast_aalookup.c`: `BlastAaLookupTableNew`
///
/// NCBI uses BLASTAA_SIZE = 28 for all amino acid lookup operations.
pub const LOOKUP_WORD_LENGTH: usize = 3;
pub const LOOKUP_ALPHABET_SIZE: usize = BLASTAA_SIZE; // 28 (NCBI standard)

#[inline(always)]
fn ilog2(mut x: usize) -> usize {
    let mut r = 0usize;
    while x > 1 {
        x >>= 1;
        r += 1;
    }
    r
}

#[inline(always)]
fn compute_backbone_size(word_length: usize, alphabet_size: usize, charsize: usize) -> usize {
    // NCBI: for (i=0;i<word_length;i++) backbone_size |= (alphabet_size-1) << (i*charsize);
    //       backbone_size++;
    let mut backbone_size: usize = 0;
    for i in 0..word_length {
        backbone_size |= (alphabet_size - 1) << (i * charsize);
    }
    backbone_size + 1
}

#[inline(always)]
fn compute_mask(word_length: usize, charsize: usize) -> usize {
    (1usize << (word_length * charsize)) - 1
}

#[inline(always)]
fn encode_kmer_3(aa0: usize, aa1: usize, aa2: usize, charsize: usize) -> usize {
    (aa0 << (2 * charsize)) | (aa1 << charsize) | aa2
}

/// Encode a 3-mer starting at `pos` in `seq` using the NCBI-style bit-shift index.
/// Uses BLASTAA_SIZE = 28 alphabet.
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + (LOOKUP_WORD_LENGTH - 1) >= seq.len() {
        return None;
    }
    let a0 = seq[pos] as usize;
    let a1 = seq[pos + 1] as usize;
    let a2 = seq[pos + 2] as usize;
    // NCBI accepts all residues 0-27
    if a0 >= LOOKUP_ALPHABET_SIZE || a1 >= LOOKUP_ALPHABET_SIZE || a2 >= LOOKUP_ALPHABET_SIZE {
        return None;
    }
    let charsize = ilog2(LOOKUP_ALPHABET_SIZE) + 1;
    let mask = compute_mask(LOOKUP_WORD_LENGTH, charsize);
    Some(encode_kmer_3(a0, a1, a2, charsize) & mask)
}

#[inline(always)]
pub fn pv_test(pv: &[u64], index: usize) -> bool {
    unsafe { (*pv.get_unchecked(index >> PV_ARRAY_BTS) & (1u64 << (index & PV_ARRAY_MASK))) != 0 }
}

/// Get the charsize (bits per residue) for the lookup table alphabet.
#[inline(always)]
pub const fn get_charsize() -> usize {
    5 // ilog2(28) + 1 = 5
}

/// Get the mask for k-mer indices.
#[inline(always)]
pub const fn get_mask() -> usize {
    (1usize << (LOOKUP_WORD_LENGTH * get_charsize())) - 1
}

/// Encode a 3-mer from residue values.
#[inline(always)]
pub fn encode_kmer(aa0: usize, aa1: usize, aa2: usize) -> usize {
    let charsize = get_charsize();
    encode_kmer_3(aa0, aa1, aa2, charsize) & get_mask()
}

/// Decode a k-mer index to residue values.
#[inline(always)]
pub fn decode_kmer(index: usize) -> (usize, usize, usize) {
    let charsize = get_charsize();
    let residue_mask = (1usize << charsize) - 1;
    let w0 = (index >> (2 * charsize)) & residue_mask;
    let w1 = (index >> charsize) & residue_mask;
    let w2 = index & residue_mask;
    (w0, w1, w2)
}

#[inline(always)]
fn pv_set(pv: &mut [u64], index: usize) {
    pv[index >> PV_ARRAY_BTS] |= 1u64 << (index & PV_ARRAY_MASK);
}

/// NCBI AaLookupBackboneCell - EXACT port
#[repr(C)]
#[derive(Clone, Copy)]
pub struct BackboneCell {
    pub num_used: i32,
    pub entries: [i32; AA_HITS_PER_CELL],
}

impl Default for BackboneCell {
    fn default() -> Self {
        Self { num_used: 0, entries: [0; AA_HITS_PER_CELL] }
    }
}

/// NCBI-style lookup table with flat backbone
pub struct BlastAaLookupTable {
    pub backbone: Vec<BackboneCell>,
    pub overflow: Vec<i32>,
    pub pv: Vec<u64>,
    pub frame_bases: Vec<i32>,
    pub num_contexts: usize,
    pub word_length: i32,
    pub alphabet_size: i32,
    pub charsize: i32,
    pub mask: i32,
    pub longest_chain: i32,
    // Lazy neighbor generation support
    pub lazy_neighbors: bool,
    pub threshold: i32,
    pub row_max: Vec<i32>,
}

#[derive(Clone)]
pub struct QueryContext {
    pub q_idx: u32,
    pub f_idx: u8,
    pub frame: i8,
    pub aa_seq: Vec<u8>,
    pub aa_seq_nomask: Option<Vec<u8>>,
    pub aa_len: usize,
    pub orig_len: usize,
    pub frame_base: i32,
    /// NCBI: kbp[context] - Karlin parameters for this context
    /// Reference: link_hsps.c line 750-752, 866-867
    pub karlin_params: KarlinParams,
}

impl BlastAaLookupTable {
    #[inline]
    pub fn get_context_idx(&self, concat_off: i32) -> usize {
        let bases = &self.frame_bases;
        let mut lo = 0usize;
        let mut hi = self.num_contexts;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if concat_off < bases[mid] {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        lo.saturating_sub(1)
    }

    #[inline(always)]
    pub fn get_hits(&self, index: usize) -> &[i32] {
        let cell = unsafe { self.backbone.get_unchecked(index) };
        let num = cell.num_used as usize;
        if num == 0 {
            return &[];
        }
        if num <= AA_HITS_PER_CELL {
            unsafe { std::slice::from_raw_parts(cell.entries.as_ptr(), num) }
        } else {
            let cursor = cell.entries[0] as usize;
            unsafe { std::slice::from_raw_parts(self.overflow.as_ptr().add(cursor), num) }
        }
    }
}

fn compute_unmasked_intervals(seg_masks: &[(usize, usize)], aa_len: usize) -> Vec<(usize, usize)> {
    if seg_masks.is_empty() {
        return vec![(0, aa_len)];
    }
    let mut result = Vec::new();
    let mut pos = 0;
    let mut sorted = seg_masks.to_vec();
    sorted.sort_by_key(|m| m.0);
    for &(s, e) in &sorted {
        if pos < s {
            result.push((pos, s));
        }
        pos = pos.max(e);
    }
    if pos < aa_len {
        result.push((pos, aa_len));
    }
    result
}

/// Build NCBI-style lookup table with BLASTAA_SIZE = 28
/// Reference: blast_aalookup.c BlastAaLookupIndexQuery + BlastAaLookupFinalize
///
/// If `lazy_neighbors` is true, only exact matches are stored in the lookup table.
/// Neighbor matching is performed dynamically during scan, which drastically reduces
/// lookup table size but requires neighbor computation at scan time.
pub fn build_ncbi_lookup(
    queries: &[Vec<QueryFrame>],
    threshold: i32,
    _include_stop_seeds: bool, // Ignored - NCBI always uses full 28-char alphabet
    _ncbi_stop_stop_score: bool, // Ignored - always use NCBI BLOSUM62 (*-* = +1)
    lazy_neighbors: bool,
    _karlin_params: &KarlinParams, // Unused - computed per context, kept for API compatibility
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let word_length = LOOKUP_WORD_LENGTH;
    let alphabet_size = LOOKUP_ALPHABET_SIZE; // 28
    let charsize = ilog2(alphabet_size) + 1;   // 5
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);

    // Compute ideal Karlin parameters (kbp_ideal) - used for check_ideal logic
    // Reference: NCBI blast_stat.c:2754 Blast_ScoreBlkKbpIdealCalc
    let ideal_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    
    // Compute standard amino acid composition (for database/subject)
    // Reference: NCBI blast_stat.c:2759 Blast_ResFreqStdComp
    let std_comp = compute_std_aa_composition();
    
    // Build contexts and frame bases
    let mut frame_bases: Vec<i32> = Vec::new();
    let mut contexts: Vec<QueryContext> = Vec::new();
    let mut base: i32 = 0;

    for (q_idx, frames) in queries.iter().enumerate() {
        for (f_idx, frame) in frames.iter().enumerate() {
            frame_bases.push(base);
            
            // Compute context-specific Karlin parameters
            // Reference: NCBI blast_stat.c:2778-2797
            // 1. Compute amino acid composition for this context
            let ctx_comp = compute_aa_composition(&frame.aa_seq, frame.aa_len);
            
            // 2. Compute score frequency profile
            // Use standard composition for subject (database composition)
            let score_min = -4; // BLOSUM62 minimum
            let score_max = 11; // BLOSUM62 maximum
            let sfp = compute_score_freq_profile(&ctx_comp, &std_comp, score_min, score_max);
            
            // 3. Compute Karlin parameters
            let computed_params = compute_karlin_params_ungapped(&sfp)
                .unwrap_or_else(|_| {
                    // Fallback to ideal if computation fails
                    ideal_params
                });
            
            // 4. Apply check_ideal logic (tblastx uses check_ideal = TRUE)
            // Reference: NCBI blast_stat.c:2796-2797
            let final_params = apply_check_ideal(computed_params, ideal_params);
            
            contexts.push(QueryContext {
                q_idx: q_idx as u32,
                f_idx: f_idx as u8,
                frame: frame.frame,
                aa_seq: frame.aa_seq.clone(),
                aa_seq_nomask: frame.aa_seq_nomask.clone(),
                aa_len: frame.aa_len,
                orig_len: frame.orig_len,
                frame_base: base,
                karlin_params: final_params,
            });
            // NCBI concatenates translated frames by *sharing* the boundary sentinel NULLB.
            // BLAST_GetTranslation writes both leading+trailing NULLB; the next frame starts
            // at the previous frame's trailing NULLB, so the offset advances by (aa_len + 1),
            // not (aa_len + 2).
            //
            // NCBI reference (verbatim):
            //   /* Increment offset by 1 extra byte for the sentinel NULLB
            //      between frames. */
            //   offset += length + 1;
            //   frame_offsets[context+1] = offset;
            // Source: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
            base += frame.aa_seq.len() as i32 - 1;
        }
    }

    // Row max for BLOSUM62 over the lookup alphabet.
    // For NCBISTDAA residues, we compute max score against any other residue.
    let row_max: Vec<i32> = (0..alphabet_size)
        .map(|i| {
            (0..alphabet_size)
                .map(|j| blosum62_score(i as u8, j as u8))
                .max()
                .unwrap_or(-4)
        })
        .collect();

    // Phase 1a: Index exact query words
    let mut exact_offsets: Vec<Vec<i32>> = vec![Vec::new(); backbone_size];
    let mut total_exact_positions = 0usize;
    let mut skipped_invalid_residue = 0usize;
    let mut skipped_seg_mask = 0usize;

    let mut ctx_idx = 0usize;
    for (q_idx, frames) in queries.iter().enumerate() {
        for frame in frames.iter() {
            let frame_base = frame_bases[ctx_idx];
            let seq = &frame.aa_seq;

            if seq.len() >= 5 && frame.aa_len >= word_length {
                let unmasked = compute_unmasked_intervals(&frame.seg_masks, frame.aa_len);
                let total_aa_positions = frame.aa_len.saturating_sub(word_length - 1);
                let unmasked_positions: usize = unmasked.iter()
                    .map(|(s, e)| e.saturating_sub(word_length - 1).saturating_sub(*s))
                    .sum();
                skipped_seg_mask += total_aa_positions.saturating_sub(unmasked_positions);
                
                for (start, end) in unmasked {
                    for aa_pos in start..end.saturating_sub(word_length - 1) {
                        let raw_pos = aa_pos + 1;
                        let c0 = seq[raw_pos] as usize;
                        let c1 = seq[raw_pos + 1] as usize;
                        let c2 = seq[raw_pos + 2] as usize;
                        
                        // NCBI accepts all valid NCBISTDAA residues (0-27)
                        if c0 >= alphabet_size || c1 >= alphabet_size || c2 >= alphabet_size {
                            skipped_invalid_residue += 1;
                            continue;
                        }

                        let idx = encode_kmer_3(c0, c1, c2, charsize) & mask;
                        exact_offsets[idx].push(frame_base + raw_pos as i32);
                        total_exact_positions += 1;
                    }
                }
            }
            ctx_idx += 1;
        }
    }
    
    // Phase 1a diagnostics
    let unique_exact_words = exact_offsets.iter().filter(|v| !v.is_empty()).count();
    let max_exact_per_word = exact_offsets.iter().map(|v| v.len()).max().unwrap_or(0);
    eprintln!("\n=== Phase 1a: Exact Indexing Diagnostics ===");
    eprintln!("Total exact positions indexed: {}", total_exact_positions);
    eprintln!("Unique exact words: {}", unique_exact_words);
    eprintln!("Max offsets per exact word: {}", max_exact_per_word);
    eprintln!("Skipped (invalid residue): {}", skipped_invalid_residue);
    eprintln!("Skipped (SEG mask): {}", skipped_seg_mask);

    // Phase 1b: Add neighboring words (or just exact matches for lazy mode)
    let residue_mask: usize = (1usize << charsize) - 1;
    
    // Diagnostics for neighbor generation
    let mut exact_added_count = 0usize;
    let mut neighbor_added_count = 0usize;
    let mut words_with_exact_only = 0usize;
    let mut words_with_neighbors = 0usize;
    let mut max_neighbors_for_single_word = 0usize;
    let mut max_neighbor_word_idx = 0usize;
    let mut neighbor_words_generated = 0usize;

    let thin_backbone: Vec<Vec<i32>> = if lazy_neighbors {
        // LAZY MODE: Only store exact matches, neighbors computed at scan time
        eprintln!("\n=== Phase 1b: LAZY NEIGHBOR MODE ===");
        eprintln!("Neighbors will be computed dynamically during scan");
        
        let mut backbone: Vec<Vec<i32>> = vec![Vec::new(); backbone_size];
        
        for idx in 0..backbone_size {
            let offsets = &exact_offsets[idx];
            if !offsets.is_empty() {
                backbone[idx] = offsets.clone();
                exact_added_count += offsets.len();
                words_with_exact_only += 1;
            }
        }
        
        backbone
    } else {
        // STANDARD MODE: Pre-compute all neighbors at build time
        // Pass 1: Count entries for each target k-mer
        let mut entry_counts: Vec<u32> = vec![0; backbone_size];
        
        for idx in 0..backbone_size {
            let offsets = &exact_offsets[idx];
            if offsets.is_empty() {
                continue;
            }
            let num_offsets = offsets.len() as u32;

            let w0 = (idx >> (2 * charsize)) & residue_mask;
            let w1 = (idx >> charsize) & residue_mask;
            let w2 = idx & residue_mask;
            if w0 >= alphabet_size || w1 >= alphabet_size || w2 >= alphabet_size {
                continue;
            }

            let self_score = blosum62_score(w0 as u8, w0 as u8)
                + blosum62_score(w1 as u8, w1 as u8)
                + blosum62_score(w2 as u8, w2 as u8);
            
            // NCBI: if threshold==0 or self_score < threshold, add exact matches
            if threshold == 0 || self_score < threshold {
                entry_counts[idx] += num_offsets;
            }
            
            if threshold == 0 {
                continue;
            }

            // Count neighbor contributions
            let rm12 = row_max[w1] + row_max[w2];
            let rm2 = row_max[w2];
            for s0 in 0..alphabet_size {
                let sc0 = blosum62_score(w0 as u8, s0 as u8);
                if sc0 + rm12 < threshold {
                    continue;
                }
                for s1 in 0..alphabet_size {
                    let sc1 = sc0 + blosum62_score(w1 as u8, s1 as u8);
                    if sc1 + rm2 < threshold {
                        continue;
                    }
                    for s2 in 0..alphabet_size {
                        if sc1 + blosum62_score(w2 as u8, s2 as u8) >= threshold {
                            let nidx = encode_kmer_3(s0, s1, s2, charsize) & mask;
                            entry_counts[nidx] += num_offsets;
                        }
                    }
                }
            }
        }
        
        // Pass 2: Pre-allocate with exact capacity
        let mut backbone: Vec<Vec<i32>> = entry_counts.iter()
            .map(|&count| Vec::with_capacity(count as usize))
            .collect();

        // Pass 3: Add entries (no reallocation needed)
        for idx in 0..backbone_size {
            let offsets = &exact_offsets[idx];
            if offsets.is_empty() {
                continue;
            }

            let w0 = (idx >> (2 * charsize)) & residue_mask;
            let w1 = (idx >> charsize) & residue_mask;
            let w2 = idx & residue_mask;
            if w0 >= alphabet_size || w1 >= alphabet_size || w2 >= alphabet_size {
                continue;
            }

            let self_score = blosum62_score(w0 as u8, w0 as u8)
                + blosum62_score(w1 as u8, w1 as u8)
                + blosum62_score(w2 as u8, w2 as u8);
            
            // NCBI: if threshold==0 or self_score < threshold, add exact matches
            if threshold == 0 || self_score < threshold {
                backbone[idx].extend_from_slice(offsets);
                exact_added_count += offsets.len();
                words_with_exact_only += 1;
            }
            
            if threshold == 0 {
                continue;
            }

            // Add neighbors - use extend_from_slice for batch addition
            let mut local_neighbor_count = 0usize;
            let rm12 = row_max[w1] + row_max[w2];
            let rm2 = row_max[w2];
            for s0 in 0..alphabet_size {
                let sc0 = blosum62_score(w0 as u8, s0 as u8);
                if sc0 + rm12 < threshold {
                    continue;
                }
                for s1 in 0..alphabet_size {
                    let sc1 = sc0 + blosum62_score(w1 as u8, s1 as u8);
                    if sc1 + rm2 < threshold {
                        continue;
                    }
                    for s2 in 0..alphabet_size {
                        if sc1 + blosum62_score(w2 as u8, s2 as u8) >= threshold {
                            let nidx = encode_kmer_3(s0, s1, s2, charsize) & mask;
                            // Use extend_from_slice for efficient batch addition
                            backbone[nidx].extend_from_slice(offsets);
                            neighbor_added_count += offsets.len();
                            local_neighbor_count += 1;
                        }
                    }
                }
            }
            
            neighbor_words_generated += local_neighbor_count;
            if local_neighbor_count > 0 {
                words_with_neighbors += 1;
                if local_neighbor_count > max_neighbors_for_single_word {
                    max_neighbors_for_single_word = local_neighbor_count;
                    max_neighbor_word_idx = idx;
                }
            }
        }
        
        backbone
    };
    
    // Use thin_backbone directly as all_entries
    let mut counts: Vec<u32> = vec![0; backbone_size];
    let mut all_entries = thin_backbone;
    for idx in 0..backbone_size {
        counts[idx] = all_entries[idx].len() as u32;
    }
    
    // Phase 1b diagnostics
    eprintln!("\n=== Phase 1b: Neighbor Generation Diagnostics ===");
    eprintln!("Threshold: {}", threshold);
    eprintln!("Exact entries added (self_score < threshold): {}", exact_added_count);
    eprintln!("Neighbor entries added: {}", neighbor_added_count);
    eprintln!("Neighbor words generated (unique): {}", neighbor_words_generated);
    eprintln!("Words with exact-only additions: {}", words_with_exact_only);
    eprintln!("Words with neighbor generation: {}", words_with_neighbors);
    eprintln!("Max neighbors generated for single word: {}", max_neighbors_for_single_word);
    if max_neighbors_for_single_word > 0 {
        let w0 = (max_neighbor_word_idx >> (2 * charsize)) & residue_mask;
        let w1 = (max_neighbor_word_idx >> charsize) & residue_mask;
        let w2 = max_neighbor_word_idx & residue_mask;
        eprintln!("  Word with max neighbors: idx={} (residues {},{},{})", 
                  max_neighbor_word_idx, w0, w1, w2);
    }
    let expansion_factor = if total_exact_positions > 0 {
        (exact_added_count + neighbor_added_count) as f64 / total_exact_positions as f64
    } else {
        0.0
    };
    eprintln!("Expansion factor (total_entries / exact_positions): {:.2}x", expansion_factor);

    // Phase 2: Finalize
    let mut backbone: Vec<BackboneCell> = vec![BackboneCell::default(); backbone_size];
    let pv_size = (backbone_size + PV_BUCKET_BITS - 1) / PV_BUCKET_BITS;
    let mut pv: Vec<u64> = vec![0u64; pv_size];

    // Pre-finalize diagnostics: analyze high-frequency k-mers
    eprintln!("\n=== Phase 2: High-Frequency K-mer Analysis ===");
    
    // Collect (count, idx) pairs and sort by count descending
    let mut count_idx_pairs: Vec<(usize, usize)> = counts.iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(i, &c)| (c as usize, i))
        .collect();
    count_idx_pairs.sort_by(|a, b| b.0.cmp(&a.0));
    
    // Show top 10 high-frequency k-mers
    eprintln!("Top 10 high-frequency k-mers (before suppression):");
    for (rank, &(count, idx)) in count_idx_pairs.iter().take(10).enumerate() {
        let w0 = (idx >> (2 * charsize)) & residue_mask;
        let w1 = (idx >> charsize) & residue_mask;
        let w2 = idx & residue_mask;
        eprintln!("  #{}: count={}, idx={} (residues {},{},{})", 
                  rank + 1, count, idx, w0, w1, w2);
    }
    
    // Distribution histogram
    let hist_buckets = [1, 10, 100, 500, 1000, 5000, 10000, 50000, 100000, usize::MAX];
    let mut hist: Vec<usize> = vec![0; hist_buckets.len()];
    let mut hist_entries: Vec<usize> = vec![0; hist_buckets.len()];
    for &(count, _) in &count_idx_pairs {
        for (i, &bucket) in hist_buckets.iter().enumerate() {
            if count <= bucket {
                hist[i] += 1;
                hist_entries[i] += count;
                break;
            }
        }
    }
    eprintln!("Hit count distribution:");
    let mut prev = 0;
    for (i, &bucket) in hist_buckets.iter().enumerate() {
        if hist[i] > 0 {
            let label = if bucket == usize::MAX {
                format!(">{}", prev)
            } else {
                format!("{}-{}", prev + 1, bucket)
            };
            eprintln!("  {}: {} cells, {} entries", label, hist[i], hist_entries[i]);
        }
        prev = bucket;
    }

    let mut overflow_size = 0usize;
    // NCBI BLAST does not filter over-represented k-mers - all k-mers are kept
    // Reference: blast_lookup.c:33-77 (BlastLookupAddWordHit simply adds hits without filtering)
    for idx in 0..backbone_size {
        let count = counts[idx] as usize;
        if count > AA_HITS_PER_CELL {
            overflow_size += count;
        }
    }

    let mut overflow: Vec<i32> = vec![0; overflow_size];
    let mut overflow_cursor = 0usize;

    let mut longest_chain: i32 = 0;
    for idx in 0..backbone_size {
        let entries = &all_entries[idx];
        let count = entries.len();
        if count == 0 {
            continue;
        }

        if (count as i32) > longest_chain {
            longest_chain = count as i32;
        }

        pv_set(&mut pv, idx);
        backbone[idx].num_used = count as i32;

        if count <= AA_HITS_PER_CELL {
            for (i, &off) in entries.iter().enumerate() {
                backbone[idx].entries[i] = off;
            }
        } else {
            backbone[idx].entries[0] = overflow_cursor as i32;
            for &off in entries {
                overflow[overflow_cursor] = off;
                overflow_cursor += 1;
            }
        }
    }

    let nonempty = backbone.iter().filter(|c| c.num_used > 0).count();
    let total: usize = backbone.iter().map(|c| c.num_used.max(0) as usize).sum();
    eprintln!("\n=== NCBI Backbone Lookup Table (BLASTAA_SIZE={}) ===", LOOKUP_ALPHABET_SIZE);
    eprintln!("Contexts: {}", contexts.len());
    eprintln!(
        "Backbone size: {} cells ({:.1} MB)",
        backbone_size,
        (backbone_size * std::mem::size_of::<BackboneCell>()) as f64 / 1e6
    );
    eprintln!("Non-empty cells: {}", nonempty);
    eprintln!("Total entries: {}", total);
    eprintln!(
        "Overflow size: {} ({:.1} MB)",
        overflow_cursor,
        (overflow_cursor * 4) as f64 / 1e6
    );
    eprintln!("Longest chain: {}", longest_chain);
    eprintln!("==========================================\n");

    (
        BlastAaLookupTable {
            backbone,
            overflow,
            pv,
            frame_bases,
            num_contexts: contexts.len(),
            word_length: word_length as i32,
            alphabet_size: alphabet_size as i32,
            charsize: charsize as i32,
            mask: mask as i32,
            longest_chain,
            lazy_neighbors,
            threshold,
            row_max,
        },
        contexts,
    )
}

/// Build a simple direct (exact) 3-mer lookup table for tests and diagnostics.
pub fn build_direct_lookup(
    queries: &[Vec<QueryFrame>],
) -> Vec<Vec<(u32, u8, u32)>> {
    let word_length = LOOKUP_WORD_LENGTH;
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);

    let mut table: Vec<Vec<(u32, u8, u32)>> = vec![Vec::new(); backbone_size];

    for (q_idx, frames) in queries.iter().enumerate() {
        for (f_idx, frame) in frames.iter().enumerate() {
            if frame.aa_len < 3 || frame.aa_seq.len() < 5 {
                continue;
            }

            let unmasked = compute_unmasked_intervals(&frame.seg_masks, frame.aa_len);
            for (start, end) in unmasked {
                for aa_pos in start..end.saturating_sub(2) {
                    let raw_pos = aa_pos + 1;
                    let c0 = frame.aa_seq[raw_pos] as usize;
                    let c1 = frame.aa_seq[raw_pos + 1] as usize;
                    let c2 = frame.aa_seq[raw_pos + 2] as usize;
                    if c0 >= alphabet_size || c1 >= alphabet_size || c2 >= alphabet_size {
                        continue;
                    }

                    let idx = encode_kmer_3(c0, c1, c2, charsize) & mask;
                    table[idx].push((q_idx as u32, f_idx as u8, aa_pos as u32));
                }
            }
        }
    }

    table
}

/// Pre-computed neighbor map: for each k-mer, list of k-mers that are neighbors (score >= threshold)
/// This is the key optimization: compute once, use during every scan.
pub struct NeighborMap {
    /// For each k-mer index, the list of k-mer indices that are neighbors
    /// neighbor_map[subject_kmer] = [query_kmer_1, query_kmer_2, ...] such that
    /// score(query_kmer, subject_kmer) >= threshold
    pub map: Vec<Vec<u16>>,
    pub threshold: i32,
    pub backbone_size: usize,
}

impl NeighborMap {
    /// Pre-compute all neighbor relationships for a given threshold.
    /// This is O(alphabet^6) but only done once.
    pub fn new(threshold: i32) -> Self {
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        let backbone_size = compute_backbone_size(LOOKUP_WORD_LENGTH, alphabet_size, charsize);
        
        eprintln!("Pre-computing neighbor map (threshold={})...", threshold);
        let start = std::time::Instant::now();
        
        // For each possible k-mer (subject side), find all k-mers (query side) 
        // that would be neighbors
        let mut map: Vec<Vec<u16>> = vec![Vec::new(); backbone_size];
        
        // Row max for pruning
        let row_max: Vec<i32> = (0..alphabet_size)
            .map(|i| {
                (0..alphabet_size)
                    .map(|j| blosum62_score(i as u8, j as u8))
                    .max()
                    .unwrap_or(-4)
            })
            .collect();
        
        let mask = compute_mask(LOOKUP_WORD_LENGTH, charsize);
        let residue_mask = (1usize << charsize) - 1;
        
        // For each subject k-mer S
        for s_idx in 0..backbone_size {
            let s0 = (s_idx >> (2 * charsize)) & residue_mask;
            let s1 = (s_idx >> charsize) & residue_mask;
            let s2 = s_idx & residue_mask;
            
            if s0 >= alphabet_size || s1 >= alphabet_size || s2 >= alphabet_size {
                continue;
            }
            
            // Find all query k-mers Q such that score(Q, S) >= threshold
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
                        if sc1 + blosum62_score(q2 as u8, s2 as u8) >= threshold {
                            let q_idx = encode_kmer_3(q0, q1, q2, charsize) & mask;
                            map[s_idx].push(q_idx as u16);
                        }
                    }
                }
            }
        }
        
        let total_neighbors: usize = map.iter().map(|v| v.len()).sum();
        let nonempty = map.iter().filter(|v| !v.is_empty()).count();
        let elapsed = start.elapsed();
        
        eprintln!("Neighbor map computed in {:?}", elapsed);
        eprintln!("  Non-empty entries: {}", nonempty);
        eprintln!("  Total neighbor pairs: {}", total_neighbors);
        eprintln!("  Avg neighbors per k-mer: {:.1}", total_neighbors as f64 / nonempty.max(1) as f64);
        
        Self {
            map,
            threshold,
            backbone_size,
        }
    }
    
    /// Get neighbors for a subject k-mer
    #[inline]
    pub fn get_neighbors(&self, subject_kmer: usize) -> &[u16] {
        &self.map[subject_kmer]
    }
    
    /// Build expanded lookup table using neighbor map.
    /// This is more efficient than generating neighbors per query position.
    pub fn build_expanded_lookup(
        &self,
        query_lookup: &[Vec<(u32, u8, u32)>],
    ) -> Vec<Vec<(u32, u8, u32)>> {
        let start = std::time::Instant::now();
        
        // For each subject k-mer, collect all query positions that would match
        let mut expanded: Vec<Vec<(u32, u8, u32)>> = vec![Vec::new(); self.backbone_size];
        
        // Build reverse neighbor map: query_kmer -> [subject_kmer, ...]
        // This is essentially: for which subject k-mers would this query k-mer be a neighbor?
        // In our current map: subject_kmer -> [query_kmer, ...]
        // Reversed: query_kmer -> [subject_kmer, ...]
        let mut reverse_map: Vec<Vec<u16>> = vec![Vec::new(); self.backbone_size];
        for (s_kmer, neighbors) in self.map.iter().enumerate() {
            for &q_kmer in neighbors {
                reverse_map[q_kmer as usize].push(s_kmer as u16);
            }
        }
        
        // Now, for each query k-mer with positions, add those positions to all
        // subject k-mers that are its neighbors
        for (q_kmer, positions) in query_lookup.iter().enumerate() {
            if positions.is_empty() {
                continue;
            }
            // Get all subject k-mers that would match this query k-mer
            let subject_kmers = &reverse_map[q_kmer];
            for &s_kmer in subject_kmers {
                expanded[s_kmer as usize].extend_from_slice(positions);
            }
        }
        
        let total_entries: usize = expanded.iter().map(|v| v.len()).sum();
        let nonempty = expanded.iter().filter(|v| !v.is_empty()).count();
        let elapsed = start.elapsed();
        
        eprintln!("Expanded lookup built in {:?}", elapsed);
        eprintln!("  Total entries: {}", total_entries);
        eprintln!("  Non-empty buckets: {}", nonempty);
        
        expanded
    }
    
    /// Pre-compute all neighbor relationships, but only keep those where
    /// the query k-mer actually exists in the query lookup table.
    /// This dramatically reduces the number of neighbors to check during scan.
    pub fn new_filtered(threshold: i32, query_pv: &[u64]) -> Self {
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        let charsize = ilog2(alphabet_size) + 1;
        let backbone_size = compute_backbone_size(LOOKUP_WORD_LENGTH, alphabet_size, charsize);
        
        eprintln!("Pre-computing filtered neighbor map (threshold={})...", threshold);
        let start = std::time::Instant::now();
        
        // For each possible k-mer (subject side), find all k-mers (query side) 
        // that would be neighbors AND exist in the query
        let mut map: Vec<Vec<u16>> = vec![Vec::new(); backbone_size];
        
        // Row max for pruning
        let row_max: Vec<i32> = (0..alphabet_size)
            .map(|i| {
                (0..alphabet_size)
                    .map(|j| blosum62_score(i as u8, j as u8))
                    .max()
                    .unwrap_or(-4)
            })
            .collect();
        
        let mask = compute_mask(LOOKUP_WORD_LENGTH, charsize);
        let residue_mask = (1usize << charsize) - 1;
        
        // For each subject k-mer S
        for s_idx in 0..backbone_size {
            let s0 = (s_idx >> (2 * charsize)) & residue_mask;
            let s1 = (s_idx >> charsize) & residue_mask;
            let s2 = s_idx & residue_mask;
            
            if s0 >= alphabet_size || s1 >= alphabet_size || s2 >= alphabet_size {
                continue;
            }
            
            // Find all query k-mers Q such that score(Q, S) >= threshold
            // AND Q exists in the query
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
                        if sc1 + blosum62_score(q2 as u8, s2 as u8) >= threshold {
                            let q_idx = encode_kmer_3(q0, q1, q2, charsize) & mask;
                            // Only add if query k-mer exists
                            if pv_test(query_pv, q_idx) {
                                map[s_idx].push(q_idx as u16);
                            }
                        }
                    }
                }
            }
        }
        
        let total_neighbors: usize = map.iter().map(|v| v.len()).sum();
        let nonempty = map.iter().filter(|v| !v.is_empty()).count();
        let elapsed = start.elapsed();
        
        eprintln!("Filtered neighbor map computed in {:?}", elapsed);
        eprintln!("  Non-empty entries: {}", nonempty);
        eprintln!("  Total neighbor pairs: {}", total_neighbors);
        eprintln!("  Avg neighbors per k-mer: {:.1}", total_neighbors as f64 / nonempty.max(1) as f64);
        
        Self {
            map,
            threshold,
            backbone_size,
        }
    }
}

/// Lookup table with pre-computed neighbor map for fast scanning
pub struct NeighborLookup {
    /// Exact query k-mer positions: query_lookup[kmer] = [(q_idx, f_idx, aa_pos), ...]
    pub query_lookup: Vec<Vec<(u32, u8, u32)>>,
    /// Presence vector for query k-mers
    pub query_pv: Vec<u64>,
    /// Pre-computed neighbor map
    pub neighbor_map: NeighborMap,
    /// Frame bases for context lookup
    pub frame_bases: Vec<i32>,
    pub contexts: Vec<QueryContext>,
}

impl NeighborLookup {
    pub fn build(
        queries: &[Vec<QueryFrame>],
        threshold: i32,
        _karlin_params: &KarlinParams, // Unused - computed per context
    ) -> Self {
        // Build exact query lookup
        let query_lookup = build_direct_lookup(queries);
        
        // Build presence vector
        let pv_size = (query_lookup.len() + PV_BUCKET_BITS - 1) / PV_BUCKET_BITS;
        let mut query_pv = vec![0u64; pv_size];
        for (idx, entries) in query_lookup.iter().enumerate() {
            if !entries.is_empty() {
                pv_set(&mut query_pv, idx);
            }
        }
        
        // Compute ideal Karlin parameters (kbp_ideal) - used for check_ideal logic
        // Reference: NCBI blast_stat.c:2754 Blast_ScoreBlkKbpIdealCalc
        let ideal_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        
        // Compute standard amino acid composition (for database/subject)
        // Reference: NCBI blast_stat.c:2759 Blast_ResFreqStdComp
        let std_comp = compute_std_aa_composition();
        
        // Build contexts and frame bases
        let mut frame_bases = Vec::new();
        let mut contexts = Vec::new();
        let mut base: i32 = 0;
        
        for (q_idx, frames) in queries.iter().enumerate() {
            for (f_idx, frame) in frames.iter().enumerate() {
                frame_bases.push(base);
                
                // Compute context-specific Karlin parameters
                // Reference: NCBI blast_stat.c:2778-2797
                // 1. Compute amino acid composition for this context
                let ctx_comp = compute_aa_composition(&frame.aa_seq, frame.aa_len);
                
                // 2. Compute score frequency profile
                // Use standard composition for subject (database composition)
                let score_min = -4; // BLOSUM62 minimum
                let score_max = 11; // BLOSUM62 maximum
                let sfp = compute_score_freq_profile(&ctx_comp, &std_comp, score_min, score_max);
                
                // 3. Compute Karlin parameters
                let computed_params = compute_karlin_params_ungapped(&sfp)
                    .unwrap_or_else(|_| {
                        // Fallback to ideal if computation fails
                        ideal_params
                    });
                
                // 4. Apply check_ideal logic (tblastx uses check_ideal = TRUE)
                // Reference: NCBI blast_stat.c:2796-2797
                let final_params = apply_check_ideal(computed_params, ideal_params);
                
                contexts.push(QueryContext {
                    q_idx: q_idx as u32,
                    f_idx: f_idx as u8,
                    frame: frame.frame,
                    aa_seq: frame.aa_seq.clone(),
                    aa_seq_nomask: frame.aa_seq_nomask.clone(),
                    aa_len: frame.aa_len,
                    orig_len: frame.orig_len,
                    frame_base: base,
                    karlin_params: final_params,
                });
                // See build_ncbi_lookup() above: share the trailing NULLB sentinel between frames.
                base += frame.aa_seq.len() as i32 - 1;
            }
        }
        
        // Pre-compute neighbor map and filter by query presence
        let neighbor_map = NeighborMap::new_filtered(threshold, &query_pv);
        
        let total_query_entries: usize = query_lookup.iter().map(|v| v.len()).sum();
        eprintln!("Query lookup: {} entries (exact matches only)", total_query_entries);
        
        Self {
            query_lookup,
            query_pv,
            neighbor_map,
            frame_bases,
            contexts,
        }
    }
    
    /// Check if a query k-mer exists
    #[inline]
    pub fn has_query_kmer(&self, kmer: usize) -> bool {
        pv_test(&self.query_pv, kmer)
    }
    
    /// Get query positions for a k-mer
    #[inline]
    pub fn get_query_hits(&self, kmer: usize) -> &[(u32, u8, u32)] {
        &self.query_lookup[kmer]
    }
}

