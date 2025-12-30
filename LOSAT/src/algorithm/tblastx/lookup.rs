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

use crate::utils::dust::MaskedInterval;
use crate::utils::matrix::MATRIX;
use super::constants::MAX_HITS_PER_KMER;
use super::translation::QueryFrame;

#[inline(always)]
fn blosum62_score(aa1: usize, aa2: usize) -> i32 {
    MATRIX[aa1 * 25 + aa2] as i32
}

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
/// LOSAT currently excludes the stop-codon code (24) from seed indexing to avoid
/// runaway seed neighborhoods; therefore the effective alphabet for the lookup
/// table is 0..=23 (24 symbols).
pub const LOOKUP_WORD_LENGTH: usize = 3;
pub const LOOKUP_ALPHABET_SIZE: usize = 24; // 0..=23 (no '*')

#[inline(always)]
fn ilog2(mut x: usize) -> usize {
    // floor(log2(x)), for x >= 1
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
    // NCBI: mask = (1 << (word_length * charsize)) - 1;
    (1usize << (word_length * charsize)) - 1
}

#[inline(always)]
fn encode_kmer_3(aa0: usize, aa1: usize, aa2: usize, charsize: usize) -> usize {
    // NCBI-style index for word_length=3:
    // index = (aa0 << (2*charsize)) | (aa1 << charsize) | aa2
    (aa0 << (2 * charsize)) | (aa1 << charsize) | aa2
}

/// Encode a 3-mer starting at `pos` in `seq` using the NCBI-style bit-shift index.
///
/// Returns `None` if:
/// - the k-mer would run past the end of `seq`
/// - any residue is outside the lookup alphabet (0..=23), including stop codon (24)
pub fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + (LOOKUP_WORD_LENGTH - 1) >= seq.len() {
        return None;
    }
    let a0 = seq[pos] as usize;
    let a1 = seq[pos + 1] as usize;
    let a2 = seq[pos + 2] as usize;
    if a0 >= LOOKUP_ALPHABET_SIZE || a1 >= LOOKUP_ALPHABET_SIZE || a2 >= LOOKUP_ALPHABET_SIZE {
        return None;
    }
    let charsize = ilog2(LOOKUP_ALPHABET_SIZE) + 1;
    let mask = compute_mask(LOOKUP_WORD_LENGTH, charsize);
    Some(encode_kmer_3(a0, a1, a2, charsize) & mask)
}

/// Build a simple direct (exact) 3-mer lookup table for tests and diagnostics.
///
/// This indexes **exact** 3-mers from each translated query frame (no neighborhood generation).
/// The returned table is keyed by the same NCBI-style k-mer index used by the main lookup code.
///
/// Bucket payload: `(query_idx, frame_idx, aa_pos)` where:
/// - `query_idx` is the index into `queries`
/// - `frame_idx` is the index into `queries[query_idx]`
/// - `aa_pos` is the **logical** amino-acid position (0-based, excluding sentinels)
pub fn build_direct_lookup(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
) -> Vec<Vec<(u32, u8, u32)>> {
    let word_length = LOOKUP_WORD_LENGTH;
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);

    let mut table: Vec<Vec<(u32, u8, u32)>> = vec![Vec::new(); backbone_size];

    for (q_idx, frames) in queries.iter().enumerate() {
        let dna_masks = &query_masks[q_idx];

        for (f_idx, frame) in frames.iter().enumerate() {
            // Need at least 3 real amino acids to form a 3-mer.
            if frame.aa_len < 3 || frame.aa_seq.len() < 5 {
                continue;
            }

            let unmasked = compute_unmasked_intervals(&frame.seg_masks, frame.aa_len);
            for (start, end) in unmasked {
                // valid starts are: start <= aa_pos, aa_pos + 3 <= end  => aa_pos < end - 2
                for aa_pos in start..end.saturating_sub(2) {
                    if !dna_masks.is_empty() {
                        let (ds, de) = aa_pos_to_dna_range(aa_pos, 3, frame.frame, frame.orig_len);
                        if is_dna_pos_masked(dna_masks, ds, de) {
                            continue;
                        }
                    }

                    let raw_pos = aa_pos + 1; // +1 for sentinel at beginning
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

#[inline(always)]
pub fn pv_test(pv: &[u64], index: usize) -> bool {
    unsafe {
        (*pv.get_unchecked(index >> PV_ARRAY_BTS) & (1u64 << (index & PV_ARRAY_MASK))) != 0
    }
}

#[inline(always)]
fn pv_set(pv: &mut [u64], index: usize) {
    pv[index >> PV_ARRAY_BTS] |= 1u64 << (index & PV_ARRAY_MASK);
}

/// NCBI AaLookupBackboneCell - EXACT port
/// Reference: blast_aalookup.h:53-67
#[repr(C)]
#[derive(Clone, Copy)]
pub struct BackboneCell {
    pub num_used: i32,
    /// payload.entries[0..3] for inline storage, OR payload.entries[0] = overflow_cursor
    pub entries: [i32; AA_HITS_PER_CELL],
}

impl Default for BackboneCell {
    fn default() -> Self {
        Self { num_used: 0, entries: [0; AA_HITS_PER_CELL] }
    }
}

/// NCBI-style lookup table with flat backbone
pub struct BlastAaLookupTable {
    pub backbone: Vec<BackboneCell>,  // Flat array
    pub overflow: Vec<i32>,           // Flat overflow array
    pub pv: Vec<u64>,                 // Presence vector
    pub frame_bases: Vec<i32>,        // Context offset bases
    pub num_contexts: usize,
    /// NCBI-style indexing parameters (see `blast_aalookup.c` / `blast_lookup.h`)
    pub word_length: i32,
    pub alphabet_size: i32,
    pub charsize: i32,
    pub mask: i32,
}

#[derive(Clone)]
pub struct QueryContext {
    pub q_idx: u32,
    pub f_idx: u8,
    pub frame: i8,
    pub aa_seq: Vec<u8>,
    pub aa_len: usize,
    pub orig_len: usize,
    pub frame_base: i32,
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
    
    /// Get hits for a k-mer - returns slice of i32 offsets
    /// Reference: blast_aascan.c:100-105
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

fn aa_pos_to_dna_range(aa_pos: usize, kmer_len: usize, frame: i8, orig_len: usize) -> (usize, usize) {
    let aa_end = aa_pos + kmer_len;
    if frame > 0 {
        let offset = (frame - 1) as usize;
        (aa_pos * 3 + offset, (aa_end * 3 + offset).min(orig_len))
    } else {
        let offset = (-frame - 1) as usize;
        (orig_len.saturating_sub(aa_end * 3 + offset), orig_len - (aa_pos * 3 + offset))
    }
}

fn is_dna_pos_masked(intervals: &[MaskedInterval], start: usize, end: usize) -> bool {
    intervals.iter().any(|i| start < i.end && end > i.start)
}

fn compute_unmasked_intervals(seg_masks: &[(usize, usize)], aa_len: usize) -> Vec<(usize, usize)> {
    if seg_masks.is_empty() { return vec![(0, aa_len)]; }
    let mut result = Vec::new();
    let mut pos = 0;
    let mut sorted = seg_masks.to_vec();
    sorted.sort_by_key(|m| m.0);
    for &(s, e) in &sorted {
        if pos < s { result.push((pos, s)); }
        pos = pos.max(e);
    }
    if pos < aa_len { result.push((pos, aa_len)); }
    result
}

/// Build NCBI-style lookup table
/// Reference: blast_aalookup.c BlastAaLookupIndexQuery + BlastAaLookupFinalize
pub fn build_ncbi_lookup(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
    threshold: i32,
) -> (BlastAaLookupTable, Vec<QueryContext>) {
    // NCBI-style lookup parameters (BlastAaLookupTableNew)
    let word_length = LOOKUP_WORD_LENGTH;
    let alphabet_size = LOOKUP_ALPHABET_SIZE;
    let charsize = ilog2(alphabet_size) + 1;
    let mask = compute_mask(word_length, charsize);
    let backbone_size = compute_backbone_size(word_length, alphabet_size, charsize);

    // Build contexts and frame bases
    let mut frame_bases: Vec<i32> = Vec::new();
    let mut contexts: Vec<QueryContext> = Vec::new();
    let mut base: i32 = 0;
    
    for (q_idx, frames) in queries.iter().enumerate() {
        for (f_idx, frame) in frames.iter().enumerate() {
            frame_bases.push(base);
            contexts.push(QueryContext {
                q_idx: q_idx as u32,
                f_idx: f_idx as u8,
                frame: frame.frame,
                aa_seq: frame.aa_seq.clone(),
                aa_len: frame.aa_len,
                orig_len: frame.orig_len,
                frame_base: base,
            });
            base += frame.aa_seq.len() as i32;
        }
    }
    
    // Row max for BLOSUM62 (lookup alphabet only: 0..=23)
    let row_max: Vec<i32> = (0..alphabet_size)
        .map(|i| (0..alphabet_size).map(|j| blosum62_score(i, j)).max().unwrap_or(0))
        .collect();
    
    // Phase 1: Build thin backbone (linked lists) - count entries per cell
    let mut counts: Vec<u32> = vec![0; backbone_size];
    let mut all_entries: Vec<Vec<i32>> = vec![Vec::new(); backbone_size];
    
    let mut ctx_idx = 0usize;
    for (q_idx, frames) in queries.iter().enumerate() {
        let dna_masks = &query_masks[q_idx];
        for (_f_idx, frame) in frames.iter().enumerate() {
            let frame_base = frame_bases[ctx_idx];
            let seq = &frame.aa_seq;
            
            if seq.len() >= 5 {
                let unmasked = compute_unmasked_intervals(&frame.seg_masks, frame.aa_len);
                for (start, end) in unmasked {
                    let raw_start = (start + 1).max(1);
                    let raw_end = (end + 1).min(seq.len() - 3);
                    
                    for raw_pos in raw_start..raw_end {
                        let logical_pos = raw_pos - 1;
                        
                        if !dna_masks.is_empty() {
                            let (ds, de) = aa_pos_to_dna_range(logical_pos, 3, frame.frame, frame.orig_len);
                            if is_dna_pos_masked(dna_masks, ds, de) { continue; }
                        }
                        
                        let c0 = seq[raw_pos] as usize;
                        let c1 = seq[raw_pos + 1] as usize;
                        let c2 = seq[raw_pos + 2] as usize;
                        if c0 >= 24 || c1 >= 24 || c2 >= 24 { continue; }
                        
                        let max_score = row_max[c0] + row_max[c1] + row_max[c2];
                        if max_score < threshold { continue; }
                        
                        let concat_off = frame_base + raw_pos as i32;
                        
                        // Add neighbors
                        let rm12 = row_max[c1] + row_max[c2];
                        let rm2 = row_max[c2];
                        for s0 in 0..alphabet_size {
                            let sc0 = blosum62_score(c0, s0);
                            if sc0 + rm12 < threshold { continue; }
                            for s1 in 0..alphabet_size {
                                let sc1 = sc0 + blosum62_score(c1, s1);
                                if sc1 + rm2 < threshold { continue; }
                                for s2 in 0..alphabet_size {
                                    if sc1 + blosum62_score(c2, s2) >= threshold {
                                        // NCBI-style word index (bit-shift + mask)
                                        let idx = encode_kmer_3(s0, s1, s2, charsize) & mask;
                                        all_entries[idx].push(concat_off);
                                        counts[idx] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ctx_idx += 1;
        }
    }
    
    // Phase 2: Finalize - build thick backbone and overflow
    // Reference: blast_aalookup.c BlastAaLookupFinalize
    let mut backbone: Vec<BackboneCell> = vec![BackboneCell::default(); backbone_size];
    let pv_size = (backbone_size + PV_BUCKET_BITS - 1) / PV_BUCKET_BITS;
    let mut pv: Vec<u64> = vec![0u64; pv_size];
    
    // Calculate overflow size
    let mut overflow_size = 0usize;
    for idx in 0..backbone_size {
        let count = counts[idx] as usize;
        if count > MAX_HITS_PER_KMER {
            all_entries[idx].clear();  // Filter high-frequency
        } else if count > AA_HITS_PER_CELL {
            overflow_size += count;
        }
    }
    
    let mut overflow: Vec<i32> = vec![0; overflow_size];
    let mut overflow_cursor = 0usize;
    
    for idx in 0..backbone_size {
        let entries = &all_entries[idx];
        let count = entries.len();
        if count == 0 { continue; }
        
        pv_set(&mut pv, idx);
        backbone[idx].num_used = count as i32;
        
        if count <= AA_HITS_PER_CELL {
            // Store inline
            for (i, &off) in entries.iter().enumerate() {
                backbone[idx].entries[i] = off;
            }
        } else {
            // Store in overflow, entries[0] = cursor
            backbone[idx].entries[0] = overflow_cursor as i32;
            for &off in entries {
                overflow[overflow_cursor] = off;
                overflow_cursor += 1;
            }
        }
    }
    
    let nonempty = backbone.iter().filter(|c| c.num_used > 0).count();
    let total: usize = backbone.iter().map(|c| c.num_used.max(0) as usize).sum();
    eprintln!("\n=== NCBI Backbone Lookup Table ===");
    eprintln!("Contexts: {}", contexts.len());
    eprintln!(
        "Backbone size: {} cells ({:.1} MB)",
        backbone_size,
        (backbone_size * std::mem::size_of::<BackboneCell>()) as f64 / 1e6
    );
    eprintln!("Non-empty cells: {}", nonempty);
    eprintln!("Total entries: {}", total);
    eprintln!("Overflow size: {} ({:.1} MB)", overflow_cursor, (overflow_cursor * 4) as f64 / 1e6);
    eprintln!("==================================\n");
    
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
        },
        contexts,
    )
}
