//! Backbone lookup table for TBLASTX - EXACT NCBI BLAST port
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c

use crate::utils::matrix::blosum62_score;
use crate::stats::lookup_protein_params_ungapped;
use crate::config::ScoringMatrix;
use crate::stats::karlin_calc::{
    compute_aa_composition, compute_std_aa_composition, compute_score_freq_profile,
    compute_karlin_params_ungapped, apply_check_ideal,
};
use super::super::diagnostics::diagnostics_enabled;
use super::super::translation::QueryFrame;
use super::{
    QueryContext, LOOKUP_WORD_LENGTH, LOOKUP_ALPHABET_SIZE,
    PV_ARRAY_BTS, PV_ARRAY_MASK, ilog2, compute_backbone_size, compute_mask,
    encode_kmer_3, pv_set, compute_unmasked_intervals,
};

pub const AA_HITS_PER_CELL: usize = 3;

// Presence Vector bucket bits
const PV_BUCKET_BITS: usize = 64;

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

/// Build NCBI-style lookup table with presence vector and overflow.
///
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
    _karlin_params: &crate::stats::KarlinParams, // Unused - computed per context, kept for API compatibility
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    let diag_enabled = diagnostics_enabled();
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
    for (_q_idx, frames) in queries.iter().enumerate() {
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

    // Phase 1a diagnostics (guarded)
    if diag_enabled {
        let unique_exact_words = exact_offsets.iter().filter(|v| !v.is_empty()).count();
        let max_exact_per_word = exact_offsets.iter().map(|v| v.len()).max().unwrap_or(0);
        eprintln!("\n=== Phase 1a: Exact Indexing Diagnostics ===");
        eprintln!("Total exact positions indexed: {}", total_exact_positions);
        eprintln!("Unique exact words: {}", unique_exact_words);
        eprintln!("Max offsets per exact word: {}", max_exact_per_word);
        eprintln!("Skipped (invalid residue): {}", skipped_invalid_residue);
        eprintln!("Skipped (SEG mask): {}", skipped_seg_mask);
    }

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
        if diag_enabled {
            eprintln!("\n=== Phase 1b: LAZY NEIGHBOR MODE ===");
            eprintln!("Neighbors will be computed dynamically during scan");
        }

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

            // NCBI reference: blast_aalookup.c s_AddWordHits lines 492-519
            // Compute self-score of query word
            let self_score = blosum62_score(w0 as u8, w0 as u8)
                + blosum62_score(w1 as u8, w1 as u8)
                + blosum62_score(w2 as u8, w2 as u8);

            // NCBI: If threshold==0 OR self_score < threshold, add exact match explicitly
            // because neighbor search won't find it (score would be self_score < threshold)
            if threshold == 0 || self_score < threshold {
                entry_counts[idx] += num_offsets;
            }
            // NCBI: If threshold==0, skip neighbor search entirely
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

            // NCBI reference: blast_aalookup.c s_AddWordHits lines 492-519
            // Compute self-score of query word
            let self_score = blosum62_score(w0 as u8, w0 as u8)
                + blosum62_score(w1 as u8, w1 as u8)
                + blosum62_score(w2 as u8, w2 as u8);

            // NCBI: If threshold==0 OR self_score < threshold, add exact match explicitly
            // because neighbor search won't find it (score would be self_score < threshold)
            if threshold == 0 || self_score < threshold {
                backbone[idx].extend_from_slice(offsets);
                exact_added_count += offsets.len();
                if threshold == 0 {
                    words_with_exact_only += 1;
                }
            }
            // NCBI: If threshold==0, skip neighbor search entirely
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

    // Phase 1b diagnostics (guarded)
    if diag_enabled {
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
            eprintln!(
                "  Word with max neighbors: idx={} (residues {},{},{})",
                max_neighbor_word_idx, w0, w1, w2
            );
        }
        let expansion_factor = if total_exact_positions > 0 {
            (exact_added_count + neighbor_added_count) as f64 / total_exact_positions as f64
        } else {
            0.0
        };
        eprintln!("Expansion factor (total_entries / exact_positions): {:.2}x", expansion_factor);
    }

    // Phase 2: Finalize
    let mut backbone: Vec<BackboneCell> = vec![BackboneCell::default(); backbone_size];
    let pv_size = (backbone_size + PV_BUCKET_BITS - 1) / PV_BUCKET_BITS;
    let mut pv: Vec<u64> = vec![0u64; pv_size];

    // Pre-finalize diagnostics: analyze high-frequency k-mers (guarded)
    if diag_enabled {
        eprintln!("\n=== Phase 2: High-Frequency K-mer Analysis ===");

        // Collect (count, idx) pairs and sort by count descending
        let mut count_idx_pairs: Vec<(usize, usize)> = counts
            .iter()
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
            eprintln!(
                "  #{}: count={}, idx={} (residues {},{},{})",
                rank + 1,
                count,
                idx,
                w0,
                w1,
                w2
            );
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
    if diag_enabled {
        eprintln!(
            "\n=== NCBI Backbone Lookup Table (BLASTAA_SIZE={}) ===",
            LOOKUP_ALPHABET_SIZE
        );
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
    }

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
