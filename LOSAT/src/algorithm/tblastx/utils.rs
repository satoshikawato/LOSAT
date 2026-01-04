//! TBLASTX - LINE-BY-LINE NCBI BLAST PORT
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:439-619
//!            ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:48-131

use anyhow::{Context, Result};
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::sync::atomic::{AtomicI32, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::mpsc::channel;

use crate::common::{write_output_ncbi_order, Hit};
use crate::config::ScoringMatrix;
use crate::stats::{lookup_protein_params_ungapped, KarlinParams};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegMasker;

use super::args::TblastxArgs;
use super::chaining::UngappedHit;
use super::constants::{GAP_TRIGGER_BIT_SCORE, X_DROP_UNGAPPED_BITS};
use super::ncbi_cutoffs::{
    compute_eff_lengths_subject_mode_tblastx, cutoff_score_for_update_tblastx,
    cutoff_score_max_for_tblastx, gap_trigger_raw_score, x_drop_raw_score, BLAST_GAP_DECAY_RATE,
};
use super::diagnostics::{
    diagnostics_enabled, print_summary as print_diagnostics_summary, DiagnosticCounters,
};
use super::extension::{convert_coords, extend_hit_two_hit};
use super::lookup::{
    build_ncbi_lookup, decode_kmer, encode_kmer, get_charsize, get_mask,
    BlastAaLookupTable, NeighborLookup, AA_HITS_PER_CELL, LOOKUP_ALPHABET_SIZE,
};
use super::reevaluate::{
    get_num_identities_and_positives_ungapped, hsp_test,
    reevaluate_ungapped_hit_ncbi_translated,
};
use crate::utils::matrix::blosum62_score;
use super::sum_stats_linking::{
    apply_sum_stats_even_gap_linking, compute_avg_query_length_ncbi, 
    find_smallest_lambda_params, LinkingParams,
};
use super::translation::{generate_frames, QueryFrame};
use crate::stats::karlin::bit_score as calc_bit_score;

// [C] typedef struct DiagStruct { Int4 last_hit; Uint1 flag; } DiagStruct;
#[derive(Clone, Copy)]
struct DiagStruct {
    last_hit: i32,
    flag: u8,
}
impl Default for DiagStruct {
    fn default() -> Self {
        Self { last_hit: 0, flag: 0 }
    }
}

/// Query/subject offset pair (NCBI `BlastOffsetPair` equivalent for the hot path).
#[repr(C)]
#[derive(Clone, Copy, Default)]
struct OffsetPair {
    q_off: i32,
    s_off: i32,
}

// ---------------------------------------------------------------------------
// SIMD helpers (k-mer scan hot path)
// ---------------------------------------------------------------------------

// NCBI BLAST reference (c++/src/algo/blast/core/blast_aascan.c:99-115):
//   /* copy the hits. */
//   {
//       Int4 i;
//       Int4 s_off = (Int4)(s - subject->sequence);
//       for (i = 0; i < numhits; i++) {
//           offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
//           offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
//       }
//   }
//
// LOSAT performs the same copy, but can accelerate the (overflow) case where
// `numhits > AA_HITS_PER_CELL` by packing each (q_off, s_off) pair into u64 and
// storing 4 pairs at a time with AVX2. Output order is unchanged.

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn copy_offset_pairs_overflow_avx2(
    src: *const i32,
    dest: *mut OffsetPair,
    numhits: usize,
    s_off: i32,
) {
    use std::arch::x86_64::*;

    // Pack as: [low32=q_off | high32=s_off] (little-endian)
    let s_off_u64 = (s_off as u32 as u64) << 32;
    let s_vec = _mm256_set1_epi64x(s_off_u64 as i64);
    let low32_mask = _mm256_set1_epi64x(0xFFFF_FFFFu64 as i64);

    let mut i = 0usize;

    // 8 pairs per iteration (2x 4-pair stores)
    while i + 8 <= numhits {
        let q0 = _mm_loadu_si128(src.add(i) as *const __m128i);
        let q1 = _mm_loadu_si128(src.add(i + 4) as *const __m128i);

        // cvtepi32_epi64 sign-extends; mask to keep only low 32 bits.
        let q0_64 = _mm256_and_si256(_mm256_cvtepi32_epi64(q0), low32_mask);
        let q1_64 = _mm256_and_si256(_mm256_cvtepi32_epi64(q1), low32_mask);

        let p0 = _mm256_or_si256(q0_64, s_vec);
        let p1 = _mm256_or_si256(q1_64, s_vec);

        _mm256_storeu_si256(dest.add(i) as *mut __m256i, p0);
        _mm256_storeu_si256(dest.add(i + 4) as *mut __m256i, p1);

        i += 8;
    }

    // 4 pairs remainder
    if i + 4 <= numhits {
        let q = _mm_loadu_si128(src.add(i) as *const __m128i);
        let q64 = _mm256_and_si256(_mm256_cvtepi32_epi64(q), low32_mask);
        let p = _mm256_or_si256(q64, s_vec);
        _mm256_storeu_si256(dest.add(i) as *mut __m256i, p);
        i += 4;
    }

    // Tail
    while i < numhits {
        *dest.add(i) = OffsetPair { q_off: *src.add(i), s_off };
        i += 1;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn copy_offset_pairs_overflow_sse2(
    src: *const i32,
    dest: *mut OffsetPair,
    numhits: usize,
    s_off: i32,
) {
    use std::arch::x86_64::*;

    // Interleave q_off with constant s_off in 32-bit lanes:
    // q = [q3 q2 q1 q0]
    // s = [ s  s  s  s ]
    // lo = unpacklo(q,s) => [q1 s q0 s] (two OffsetPair)
    // hi = unpackhi(q,s) => [q3 s q2 s] (two OffsetPair)
    let s_vec = _mm_set1_epi32(s_off);

    let mut i = 0usize;
    while i + 4 <= numhits {
        let q = _mm_loadu_si128(src.add(i) as *const __m128i);
        let lo = _mm_unpacklo_epi32(q, s_vec);
        let hi = _mm_unpackhi_epi32(q, s_vec);

        _mm_storeu_si128(dest.add(i) as *mut __m128i, lo);
        _mm_storeu_si128(dest.add(i + 2) as *mut __m128i, hi);

        i += 4;
    }

    while i < numhits {
        *dest.add(i) = OffsetPair { q_off: *src.add(i), s_off };
        i += 1;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn pv_test_mask4_avx2(pv: *const u64, idxs: &[usize; 4]) -> u8 {
    use std::arch::x86_64::*;

    // NCBI BLAST reference (blast_aascan.c:83-92):
    //   index = ComputeTableIndexIncremental(...);
    //   if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
    //       numhits = bbc[index].num_used;
    //       ...
    //   }

    let w0 = (idxs[0] >> 6) as i32;
    let w1 = (idxs[1] >> 6) as i32;
    let w2 = (idxs[2] >> 6) as i32;
    let w3 = (idxs[3] >> 6) as i32;

    // Gather 4x u64 PV words (scale=8 bytes)
    let vindex = _mm_set_epi32(w3, w2, w1, w0);
    let pv_words = _mm256_i32gather_epi64(pv as *const i64, vindex, 8);

    // Compute bit masks: 1u64 << (idx & 63) per lane
    let b0 = (idxs[0] & 63) as i64;
    let b1 = (idxs[1] & 63) as i64;
    let b2 = (idxs[2] & 63) as i64;
    let b3 = (idxs[3] & 63) as i64;
    let shifts = _mm256_set_epi64x(b3, b2, b1, b0);
    let ones = _mm256_set1_epi64x(1);
    let bitmask = _mm256_sllv_epi64(ones, shifts);

    let hits = _mm256_and_si256(pv_words, bitmask);

    let mut m = 0u8;
    // Extract hit flags (lane order preserved: 0..3) without spilling to memory.
    // `_mm256_extract_epi64` requires a const index, which is fine here.
    let h0 = _mm256_extract_epi64::<0>(hits) as u64;
    let h1 = _mm256_extract_epi64::<1>(hits) as u64;
    let h2 = _mm256_extract_epi64::<2>(hits) as u64;
    let h3 = _mm256_extract_epi64::<3>(hits) as u64;

    if h0 != 0 { m |= 1; }
    if h1 != 0 { m |= 2; }
    if h2 != 0 { m |= 4; }
    if h3 != 0 { m |= 8; }
    m
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn pv_test_mask8_avx2(pv: *const u64, idxs: &[usize; 8]) -> u8 {
    let a0 = [idxs[0], idxs[1], idxs[2], idxs[3]];
    let a1 = [idxs[4], idxs[5], idxs[6], idxs[7]];

    let m0 = pv_test_mask4_avx2(pv, &a0);
    let m1 = pv_test_mask4_avx2(pv, &a1);
    m0 | (m1 << 4)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn pv_test_mask16_avx2(pv: *const u64, idxs: &[usize; 16]) -> u16 {
    let a0 = [idxs[0], idxs[1], idxs[2], idxs[3]];
    let a1 = [idxs[4], idxs[5], idxs[6], idxs[7]];
    let a2 = [idxs[8], idxs[9], idxs[10], idxs[11]];
    let a3 = [idxs[12], idxs[13], idxs[14], idxs[15]];

    let m0 = pv_test_mask4_avx2(pv, &a0) as u16;
    let m1 = pv_test_mask4_avx2(pv, &a1) as u16;
    let m2 = pv_test_mask4_avx2(pv, &a2) as u16;
    let m3 = pv_test_mask4_avx2(pv, &a3) as u16;

    m0 | (m1 << 4) | (m2 << 8) | (m3 << 12)
}

/// NCBI BlastInitHSP equivalent - stores initial HSP with absolute coordinates
/// Reference: blast_hits.h:BlastInitHSP, blast_extend.c:360-375 BlastSaveInitHsp
/// 
/// This structure stores HSPs after extension but before coordinate conversion.
/// Coordinates are stored as absolute positions in the concatenated buffer.
#[derive(Clone, Copy)]
struct InitHSP {
    /// Query absolute coordinate (in concatenated buffer, with sentinel)
    q_start_absolute: i32,
    /// Query end absolute coordinate (in concatenated buffer, with sentinel)
    q_end_absolute: i32,
    /// Subject coordinate (frame-relative, with sentinel)
    s_start: i32,
    /// Subject end coordinate (frame-relative, with sentinel)
    s_end: i32,
    /// Raw score from extension
    score: i32,
    /// Query context index
    ctx_idx: usize,
    /// Subject frame index
    s_f_idx: usize,
    /// Query record index
    q_idx: u32,
    /// Subject record index
    s_idx: u32,
    /// Query frame
    q_frame: i8,
    /// Subject frame
    s_frame: i8,
    /// Query original length
    q_orig_len: usize,
    /// Subject original length
    s_orig_len: usize,
}

struct WorkerState {
    tx: std::sync::mpsc::Sender<Vec<Hit>>,
    offset_pairs: Vec<OffsetPair>,
    diag_array: Vec<DiagStruct>,
    diag_offset: i32,
}

#[inline]
fn atomic_max_usize(dst: &AtomicUsize, val: usize) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
fn atomic_min_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val < cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
fn atomic_max_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

/// NCBI s_AdjustInitialHSPOffsets equivalent
/// Reference: blast_gapalign.c:2384-2392
/// 
/// NCBI code:
/// ```c
/// static NCBI_INLINE void
/// s_AdjustInitialHSPOffsets(BlastInitHSP* init_hsp, Int4 query_start)
/// {
///     init_hsp->offsets.qs_offsets.q_off -= query_start;
///     if (init_hsp->ungapped_data) {
///         init_hsp->ungapped_data->q_start -= query_start;  // ★座標変換
///     }
///     ASSERT(init_hsp->ungapped_data == NULL ||
///            init_hsp->ungapped_data->q_start >= 0);
/// }
/// ```
/// 
/// Converts absolute coordinates to context-relative coordinates.
/// This is called in BLAST_GetUngappedHSPList (blast_gapalign.c:4756-4758).
#[inline]
fn adjust_initial_hsp_offsets(
    hsp_q_absolute: i32,  // Absolute coordinate in concatenated buffer
    frame_base: i32,      // Context start position (absolute)
) -> i32 {
    // NCBI: init_hsp->ungapped_data->q_start -= query_start;
    hsp_q_absolute - frame_base
}

/// NCBI BLAST_GetUngappedHSPList equivalent
/// Reference: blast_gapalign.c:4719-4775
/// 
/// NCBI code flow:
/// 1. Get context for each InitHSP (s_GetUngappedHSPContext)
/// 2. Adjust coordinates (s_AdjustInitialHSPOffsets)
/// 3. Initialize HSP (Blast_HSPInit)
/// 4. Add to hsp_list
/// 
/// This function converts InitHSPs (with absolute coordinates) to UngappedHits
/// (with context-relative coordinates). Reevaluation is performed separately
/// by Blast_HSPListReevaluateUngapped equivalent.
/// 
/// NCBI reference (verbatim, blast_gapalign.c:4756-4768):
/// ```c
/// context = s_GetUngappedHSPContext(query_info, init_hsp);
/// s_AdjustInitialHSPOffsets(init_hsp, query_info->contexts[context].query_offset);
/// ungapped_data = init_hsp->ungapped_data;
/// Blast_HSPInit(ungapped_data->q_start,
///               ungapped_data->length+ungapped_data->q_start,
///               ungapped_data->s_start,
///               ungapped_data->length+ungapped_data->s_start,
///               init_hsp->offsets.qs_offsets.q_off,
///               init_hsp->offsets.qs_offsets.s_off,
///               context, query_info->contexts[context].frame,
///               subject->frame, ungapped_data->score, NULL, &new_hsp);
/// Blast_HSPListSaveHSP(hsp_list, new_hsp);
/// ```
fn get_ungapped_hsp_list(
    init_hsps: Vec<InitHSP>,
    contexts: &[super::lookup::QueryContext],
    s_frames: &[super::translation::QueryFrame],
) -> Vec<UngappedHit> {
    let mut ungapped_hits = Vec::new();
    
    for init_hsp in init_hsps {
        let ctx = &contexts[init_hsp.ctx_idx];
        let s_frame = &s_frames[init_hsp.s_f_idx];
        
        // NCBI: s_GetUngappedHSPContext equivalent
        // Context is already stored in init_hsp.ctx_idx
        
        // NCBI: s_AdjustInitialHSPOffsets (blast_gapalign.c:2384-2392)
        // Convert absolute coordinates to context-relative coordinates
        // NCBI: init_hsp->ungapped_data->q_start -= query_start;
        // where query_start = query_info->contexts[context].query_offset
        let q_start_relative = adjust_initial_hsp_offsets(init_hsp.q_start_absolute, ctx.frame_base);
        let q_end_relative = adjust_initial_hsp_offsets(init_hsp.q_end_absolute, ctx.frame_base);
        
        // NCBI: Blast_HSPInit (blast_hits.c:150-189)
        // query.offset = ungapped_data->q_start (context-relative, sentinel included)
        // query.end = ungapped_data->length + ungapped_data->q_start
        // 
        // In LOSAT, q_start_relative is context-relative coordinate (sentinel included).
        // We store it as array index into ctx.aa_seq (which has sentinel at index 0).
        let qs = q_start_relative as usize;
        let qe = q_end_relative as usize;
        let ss = init_hsp.s_start as usize;
        let se = init_hsp.s_end as usize;
        let len_u = qe.saturating_sub(qs);
        
        if len_u == 0 {
            continue;
        }
        
        // Convert to logical coords (subtract sentinel for reporting)
        // NCBI stores offsets as sentinel-inclusive, but reports as sentinel-exclusive
        let (qs_l, qe_l) = (qs.saturating_sub(1), qe.saturating_sub(1));
        let (ss_l, se_l) = (ss.saturating_sub(1), se.saturating_sub(1));
        
        // Create UngappedHit with context-relative coordinates (before reevaluation)
        // NCBI: Blast_HSPInit creates HSP with original extension score
        ungapped_hits.push(UngappedHit {
            q_idx: init_hsp.q_idx,
            s_idx: init_hsp.s_idx,
            ctx_idx: init_hsp.ctx_idx,
            s_f_idx: init_hsp.s_f_idx,
            q_frame: init_hsp.q_frame,
            s_frame: init_hsp.s_frame,
            q_aa_start: qs_l,
            q_aa_end: qe_l,
            s_aa_start: ss_l,
            s_aa_end: se_l,
            q_orig_len: init_hsp.q_orig_len,
            s_orig_len: init_hsp.s_orig_len,
            raw_score: init_hsp.score,  // Original extension score (before reevaluation)
            e_value: f64::INFINITY,
            num_ident: 0,  // Will be computed during reevaluation
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        });
    }
    
    // NCBI: Sort the HSP array by score
    // Reference: blast_gapalign.c:4772
    ungapped_hits.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));
    
    ungapped_hits
}

/// NCBI Blast_HSPListReevaluateUngapped equivalent
/// Reference: blast_hits.c:2609-2737, blast_engine.c:1492-1497
/// 
/// NCBI code flow:
/// 1. For each HSP in the list
/// 2. Get context and compute query_start (context start position)
/// 3. Call Blast_HSPReevaluateWithAmbiguitiesUngapped
/// 4. Call Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
/// 5. Delete HSPs that fail tests
/// 6. Sort by score (scores may have changed)
/// 
/// NCBI reference (verbatim, blast_hits.c:2694-2706):
/// ```c
/// context = hsp->context;
/// query_start = query_blk->sequence + query_info->contexts[context].query_offset;
/// if (kTranslateSubject)
///     subject_start = Blast_HSPGetTargetTranslation(target_t, hsp, NULL);
/// if (kNucleotideSubject) {
///     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesUngapped(hsp, query_start,
///         subject_start, word_params, sbp, kTranslateSubject);
/// }
/// ```
/// 
/// This function performs batch reevaluation on all HSPs in the list,
/// updating scores and trimming boundaries. HSPs that fail reevaluation
/// or subsequent tests are removed from the list.
fn reevaluate_ungapped_hsp_list(
    mut ungapped_hits: Vec<UngappedHit>,
    contexts: &[super::lookup::QueryContext],
    s_frames: &[super::translation::QueryFrame],
    cutoff_scores: &[Vec<i32>],
    args: &super::args::TblastxArgs,
) -> Vec<UngappedHit> {
    let mut kept_hits = Vec::new();
    
    for mut hit in ungapped_hits {
        let ctx = &contexts[hit.ctx_idx];
        let s_frame = &s_frames[hit.s_f_idx];
        let cutoff = cutoff_scores[hit.ctx_idx][hit.s_f_idx];
        
        // NCBI: query_start = query_blk->sequence + query_info->contexts[context].query_offset;
        // In LOSAT, ctx.aa_seq is the frame sequence (with sentinel at index 0)
        // query_start is the start of the context sequence (sentinel included)
        let query = &ctx.aa_seq;
        let subject = &s_frame.aa_seq;
        
        // NCBI: query = query_start + hsp->query.offset
        // hsp->query.offset is context-relative coordinate (sentinel included)
        // In LOSAT, q_aa_start is sentinel-exclusive, so we add 1 to get sentinel-inclusive
        let qs = (hit.q_aa_start + 1) as usize;  // +1 for sentinel
        let ss = (hit.s_aa_start + 1) as usize;  // +1 for sentinel
        let len_u = hit.q_aa_end.saturating_sub(hit.q_aa_start);
        
        if len_u == 0 {
            continue;
        }
        
        // NCBI: Blast_HSPReevaluateWithAmbiguitiesUngapped (blast_hits.c:2702-2705)
        // Reference: blast_hits.c:675-733
        let (new_qs, new_ss, new_len, new_score) = if let Some(result) =
            reevaluate_ungapped_hit_ncbi_translated(query, subject, qs, ss, len_u, cutoff)
        {
            result
        } else {
            // Deleted by NCBI reevaluation (score < cutoff_score)
            continue;
        };
        
        // NCBI: Blast_HSPGetNumIdentitiesAndPositives and Blast_HSPTest
        // Reference: blast_hits.c:2708-2720
        let query_nomask = ctx.aa_seq_nomask.as_ref().unwrap_or(&ctx.aa_seq);
        let num_ident = if let Some((num_ident, _align_length)) =
            get_num_identities_and_positives_ungapped(
                query_nomask,
                subject,
                new_qs,
                new_ss,
                new_len,
            ) {
            num_ident
        } else {
            continue;
        };
        
        // NCBI: Blast_HSPTest (blast_hits.c:2719)
        let percent_identity = args.percent_identity;
        let min_hit_length = args.min_hit_length;
        if hsp_test(num_ident, new_len, percent_identity, min_hit_length) {
            continue;
        }
        
        // Update hit with reevaluated coordinates and score
        // Convert from sentinel-inclusive to sentinel-exclusive for reporting
        let new_qe = new_qs + new_len;
        let new_se = new_ss + new_len;
        let (qs_l, qe_l) = (new_qs.saturating_sub(1), new_qe.saturating_sub(1));
        let (ss_l, se_l) = (new_ss.saturating_sub(1), new_se.saturating_sub(1));
        
        hit.q_aa_start = qs_l;
        hit.q_aa_end = qe_l;
        hit.s_aa_start = ss_l;
        hit.s_aa_end = se_l;
        hit.raw_score = new_score;
        hit.num_ident = num_ident;
        
        kept_hits.push(hit);
    }
    
    // NCBI: Sort the HSP array by score (scores may have changed!)
    // Reference: blast_hits.c:2734
    kept_hits.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));
    
    kept_hits
}

/// s_BlastAaScanSubject - NCBI BLAST style scan
/// Reference: blast_aascan.c:48-131
///
/// Optimized version: stop codons and invalid residues are handled by the
/// presence vector (pv_test) - buckets for k-mers containing invalid residues
/// are empty, so pv_test returns false and they are naturally skipped.
#[inline(never)]
/// NCBI-style subject sequence scan for k-mer hits.
///
/// This function scans the subject sequence for k-mer matches in the lookup table,
/// using a rolling index computation (NCBI `ComputeTableIndexIncremental`).
///
/// Optimizations applied:
/// - Unsafe array access to eliminate bounds checking in hot loop
/// - Pointer-style iteration matching NCBI's C implementation
/// - Inline presence vector test
#[inline]
fn s_blast_aa_scan_subject(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3], // [C] Int4 *s_range
) -> i32 {
    let mut totalhits: i32 = 0;

    // Precompute CPU feature flags once per scan call.
    #[cfg(target_arch = "x86_64")]
    let has_avx2 = is_x86_feature_detected!("avx2");
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

        // AVX2: batch PV_TEST for consecutive positions (rolling index stays scalar/ordered)
        #[cfg(target_arch = "x86_64")]
        if has_avx2 {
            unsafe {
                while s + 15 <= s_last {
                    // Compute 16 consecutive indices (must be sequential to preserve rolling dependency)
                    let mut idxs = [0usize; 16];
                    for i in 0..16 {
                        let new_char = *subject.get_unchecked(s + 2 + i) as usize;
                        index = ((index << charsize) | new_char) & mask;
                        idxs[i] = index;
                    }

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

                    s += 16;
                }

                while s + 7 <= s_last {
                    // Compute 8 consecutive indices (must be sequential to preserve rolling dependency)
                    let mut idxs = [0usize; 8];
                    for i in 0..8 {
                        let new_char = *subject.get_unchecked(s + 2 + i) as usize;
                        index = ((index << charsize) | new_char) & mask;
                        idxs[i] = index;
                    }

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

                    s += 8;
                }

                while s + 3 <= s_last {
                    // Compute 4 consecutive indices (must be sequential to preserve rolling dependency)
                    let mut idxs = [0usize; 4];
                    for i in 0..4 {
                        let new_char = *subject.get_unchecked(s + 2 + i) as usize;
                        index = ((index << charsize) | new_char) & mask;
                        idxs[i] = index;
                    }

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
fn s_blast_aa_scan_subject_lazy(
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

pub fn run(args: TblastxArgs) -> Result<()> {
    // Use neighbor_map mode for faster scanning with pre-computed neighbor relationships
    if args.neighbor_map {
        return run_with_neighbor_map(args);
    }
    
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    let only_qframe = args.only_qframe;
    let only_sframe = args.only_sframe;

    let valid_frame = |f: i8| matches!(f, 1 | 2 | 3 | -1 | -2 | -3);
    if let Some(f) = only_qframe {
        if !valid_frame(f) {
            anyhow::bail!("Invalid --only-qframe");
        }
    }
    if let Some(f) = only_sframe {
        if !valid_frame(f) {
            anyhow::bail!("Invalid --only-sframe");
        }
    }

    // [C] window = diag->window;
    let window = args.window_size as i32;
    // [C] wordsize = lookup->word_length;
    let wordsize: i32 = 3;

    // NCBI BLAST computes x_dropoff per-context using kbp[context]->Lambda:
    //   p->cutoffs[context].x_dropoff_init =
    //       (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // Reference: blast_parameters.c:219-221
    //
    // For tblastx, NCBI uses kbp_ideal->Lambda for all contexts (blast_stat.c:2796-2797):
    //   if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
    //      Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
    // So all contexts get the same x_dropoff, but we maintain per-context structure for parity.
    //
    // x_dropoff_per_context is populated after build_ncbi_lookup() creates the contexts.
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

    let diag_enabled = diagnostics_enabled();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    eprintln!("Reading queries...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries_raw
        .iter()
        .map(|r| r.id().split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();
    // NCBI tblastx low-complexity filtering uses SEG on translated protein sequences.
    // No nucleotide-level DUST masking is applied.
    //
    // NCBI reference (verbatim):
    //   else if (*ptr == 'L' || *ptr == 'T')
    //   { /* do low-complexity filtering; dust for blastn, otherwise seg.*/
    //       if (program_number == eBlastTypeBlastn
    //           || program_number == eBlastTypeMapping)
    //           SDustOptionsNew(&dustOptions);
    //       else
    //           SSegOptionsNew(&segOptions);
    //       ptr++;
    //   }
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:572-580

    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| {
            let mut frames = generate_frames(r.seq(), &query_code);
            if let Some(f) = only_qframe {
                frames.retain(|x| x.frame == f);
            }
            frames
        })
        .collect();

    if args.seg {
        let seg = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        for frames in &mut query_frames {
            for frame in frames {
                if frame.aa_seq.len() >= 3 {
                    for m in seg.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
                        frame.seg_masks.push((m.start, m.end));
                    }

                    // NCBI BLAST query masking semantics (SEG, etc.): keep an unmasked copy
                    // (`sequence_nomask`), then overwrite masked residues in the working
                    // query sequence buffer with X.
                    //
                    // NCBI reference (verbatim):
                    //   const Uint1 kProtMask = 21;     /* X in NCBISTDAA */
                    //   query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
                    //   query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
                    //   buffer[index] = kMaskingLetter;
                    //
                    // LOSAT uses NCBISTDAA encoding where X = 21.
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        const X_MASK_NCBISTDAA: u8 = 21; // NCBI: kProtMask = 21
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1); // keep last sentinel untouched
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!("Building lookup table...");
    // Note: karlin_params argument is now unused - computed per context in build_ncbi_lookup()
    // We still pass it for x_dropoff calculation (which uses ideal params for all contexts in tblastx)
    let (lookup, contexts) = build_ncbi_lookup(
        &query_frames,
        args.threshold,
        args.include_stop_seeds,
        args.ncbi_stop_stop_score,
        false, // lazy_neighbors disabled - use neighbor_map mode instead
        &ungapped_params_for_xdrop, // Used for x_dropoff calculation only
    );

    // NCBI BLAST: word_params->cutoffs[context].x_dropoff_init
    // Compute per-context x_dropoff using kbp[context]->Lambda.
    // Reference: blast_parameters.c:219-221
    //
    // For tblastx, all contexts use kbp_ideal (BLOSUM62 ungapped Lambda = 0.3176),
    // so x_dropoff = ceil(7.0 * LN2 / 0.3176) = 16 for all contexts.
    // We maintain per-context structure for NCBI parity.
    let x_dropoff_per_context: Vec<i32> = contexts
        .iter()
        .map(|_| x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ungapped_params_for_xdrop, 1.0))
        .collect();

    eprintln!("Reading subjects...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    // NCBI BLAST Karlin params for TBLASTX (ungapped-only algorithm):
    // 
    // TBLASTX is explicitly ungapped-only (blast_options.c line 869-873):
    //   "Gapped search is not allowed for tblastx"
    //
    // For bit score and E-value calculation, NCBI uses sbp->kbp (ungapped):
    //   blast_hits.c line 1833: kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    //   blast_hits.c line 1918: same pattern in Blast_HSPListGetBitScores
    //
    // For cutoff score search space calculation, NCBI uses gapped params:
    //   blast_parameters.c uses kbp_gap for eff_searchsp in cutoff calculation
    //
    // BLOSUM62 ungapped: lambda=0.3176, K=0.134
    // BLOSUM62 gapped:   lambda=0.267,  K=0.041
    let ungapped_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    // Gapped params for cutoff score search space calculation
    let gapped_params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    // Use UNGAPPED params for bit score and E-value (NCBI parity for tblastx)
    let params = ungapped_params.clone();

    // Compute NCBI-style average query length for linking cutoffs
    // Reference: blast_parameters.c:CalculateLinkHSPCutoffs (lines 1023-1026)
    // NCBI uses average over ALL contexts including zero-length (frame restriction via strand)
    let query_nucl_lengths: Vec<usize> = queries_raw.iter().map(|r| r.seq().len()).collect();
    let avg_query_length = compute_avg_query_length_ncbi(&query_nucl_lengths);

    eprintln!(
        "Searching {} queries vs {} subjects... (avg_query_length={})",
        queries_raw.len(),
        subjects_raw.len(),
        avg_query_length
    );

    let bar = ProgressBar::new(subjects_raw.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    let (tx, rx) = channel::<Vec<Hit>>();
    let out_path = args.out.clone();
    let evalue_threshold = args.evalue;

    let writer = std::thread::spawn(move || -> Result<()> {
        let mut all: Vec<Hit> = Vec::new();
        while let Ok(h) = rx.recv() {
            all.extend(h);
        }
        all.retain(|h| h.e_value <= evalue_threshold);
        // NCBI-style output ordering: query (input order) → subject (best_evalue/score/oid) → HSP (score/coords)
        // Reference: BLAST_LinkHsps() + s_EvalueCompareHSPLists() + ScoreCompareHSPs()
        write_output_ncbi_order(all, out_path.as_ref())?;
        Ok(())
    });

    // Diagonal array sizing MUST match NCBI's `s_BlastDiagTableNew`:
    // it depends only on (query_length + window_size), not on subject length.
    //
    // NCBI reference (verbatim):
    //   diag_array_length = 1;
    //   while (diag_array_length < (qlen+window_size))
    //       diag_array_length = diag_array_length << 1;
    //   diag_table->diag_array_length = diag_array_length;
    //   diag_table->diag_mask = diag_array_length-1;
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    //
    // `query_length` here is the concatenated translated query buffer length
    // (frames share a boundary sentinel). In LOSAT this equals the end offset
    // of the last query context buffer.
    let query_length: i32 = contexts
        .last()
        .map(|c| c.frame_base + (c.aa_seq.len() as i32 - 1))
        .unwrap_or(0);
    let mut diag_array_size: i32 = 1;
    while diag_array_size < (query_length + window) {
        diag_array_size <<= 1;
    }
    let diag_mask: i32 = diag_array_size - 1;

    // [C] array_size for offset_pairs
    // NCBI: GetOffsetArraySize() = OFFSET_ARRAY_SIZE (4096) + lookup->longest_chain
    // Reference: ncbi-blast/c++/include/algo/blast/core/lookup_wrap.h + lookup_wrap.c
    const OFFSET_ARRAY_SIZE: i32 = 4096;
    let offset_array_size: i32 = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);

    let lookup_ref = &lookup;
    let contexts_ref = &contexts;
    let query_ids_ref = &query_ids;
    let ungapped_params_ref = &ungapped_params;
    let _gapped_params_ref = &gapped_params;  // Unused - tblastx uses ungapped params

    subjects_raw.par_iter().enumerate().for_each_init(
        || WorkerState {
            tx: tx.clone(),
            offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
            diag_array: vec![DiagStruct::default(); diag_array_size as usize],
            // NCBI: diag_table->offset = window_size;
            // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
            diag_offset: window,
        },
        |st, (s_idx, s_rec)| {
            let mut s_frames = generate_frames(s_rec.seq(), &db_code);
            if let Some(f) = only_sframe {
                s_frames.retain(|x| x.frame == f);
            }

            let s_id = s_rec
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            let s_len = s_rec.seq().len();

            // [C] BlastOffsetPair *offset_pairs
            let offset_pairs = &mut st.offset_pairs;

            // [C] DiagStruct *diag_array = diag->hit_level_array;
            let diag_array = &mut st.diag_array;

            // [C] diag_offset = diag->offset;  (per-thread, reused across subjects)
            let mut diag_offset: i32 = st.diag_offset;

            // Precompute per-context cutoff scores using NCBI BLAST algorithm.
            // Reference: blast_parameters.c BlastInitialWordParametersUpdate
            //
            // NCBI cutoff calculation for tblastx ungapped path:
            // 1. gap_trigger from ungapped params (kbp_std)
            // 2. cutoff_score_max from BlastHitSavingParametersNew (uses user's E-value)
            // 3. Per-subject update: cutoff_score_for_update_tblastx with CUTOFF_E_TBLASTX=1e-300
            // 4. Final cutoff = MIN(update_cutoff, gap_trigger, cutoff_score_max)
            //
            // For -subject mode, subject_len_nucl is used (not per-frame AA length).
            let subject_len_nucl = s_len as i64;
            let mut cutoff_scores: Vec<Vec<i32>> = vec![vec![0; s_frames.len()]; contexts_ref.len()];
            // NCBI word_params->cutoff_score_min = min of cutoffs across all contexts
            // Reference: blast_parameters.c BlastInitialWordParametersUpdate line 401-403
            let mut cutoff_score_min = i32::MAX;
            
            // =======================================================================
            // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
            // =======================================================================
            // NCBI stores these in query_info->contexts[ctx].length_adjustment and
            // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths (blast_setup.c:700-850).
            // We precompute them here and pass to sum_stats_linking for NCBI parity.
            // Reference: blast_setup.c:846-847
            //   query_info->contexts[index].eff_searchsp = effective_search_space;
            //   query_info->contexts[index].length_adjustment = length_adjustment;
            let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
            let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
            
            // Precompute gap_trigger once (same for all contexts with same params)
            // Reference: blast_parameters.c:343-344
            let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ungapped_params_ref);
            
            for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
                // NCBI: query_length = query_info->contexts[context].query_length
                let query_len_aa = ctx.aa_len as i64;
                
                // =======================================================================
                // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                // length_adjustment and eff_searchsp from a single source of truth.
                // This mirrors BLAST_CalcEffLengths which computes and stores both values.
                // Reference: blast_setup.c:821-847
                // =======================================================================
                let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                    query_len_aa,
                    subject_len_nucl,
                    ungapped_params_ref,  // tblastx uses ungapped params (kbp_gap is NULL)
                );
                let eff_searchsp = eff_lengths.eff_searchsp;
                length_adj_per_context.push(eff_lengths.length_adjustment);
                eff_searchsp_per_context.push(eff_searchsp);
                
                // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
                // This uses the effective search space WITH length adjustment
                // Reference: blast_parameters.c:942-946
                let cutoff_score_max = cutoff_score_max_for_tblastx(
                    eff_searchsp,
                    evalue_threshold,  // User's E-value (typically 10.0)
                    ungapped_params_ref,
                );
                
                // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
                // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
                // Reference: blast_parameters.c:348-374
                let cutoff = cutoff_score_for_update_tblastx(
                    query_len_aa,
                    subject_len_nucl,  // NUCLEOTIDE length, NOT divided by 3!
                    gap_trigger,
                    cutoff_score_max,
                    BLAST_GAP_DECAY_RATE,  // 0.5
                    ungapped_params_ref,
                    1.0,  // scale_factor (standard BLOSUM62)
                );
                
                // DEBUG: Print cutoff values for first context
                static PRINTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
                if ctx_idx == 0 && !PRINTED.swap(true, std::sync::atomic::Ordering::Relaxed) {
                    eprintln!("[DEBUG CUTOFF] query_len_aa={}, subject_len_nucl={}", query_len_aa, subject_len_nucl);
                    eprintln!("[DEBUG CUTOFF] eff_searchsp={}", eff_searchsp);
                    eprintln!("[DEBUG CUTOFF] length_adjustment={}", eff_lengths.length_adjustment);
                    eprintln!("[DEBUG CUTOFF] cutoff_score_max={}", cutoff_score_max);
                    eprintln!("[DEBUG CUTOFF] gap_trigger={}", gap_trigger);
                    eprintln!("[DEBUG CUTOFF] final cutoff={}", cutoff);
                }
                
                // Track minimum cutoff for linking
                cutoff_score_min = cutoff_score_min.min(cutoff);
                
                // All subject frames use the same cutoff (NCBI behavior)
                for sf_idx in 0..s_frames.len() {
                    cutoff_scores[ctx_idx][sf_idx] = cutoff;
                }
            }
            // If no contexts, use 0 as fallback
            if cutoff_score_min == i32::MAX {
                cutoff_score_min = 0;
            }

            // NCBI: BlastSaveInitHsp equivalent - store initial HSPs with absolute coordinates
            // Reference: blast_extend.c:360-375 BlastSaveInitHsp
            let mut init_hsps: Vec<InitHSP> = Vec::new();
            
            // Statistics for HSP saving analysis (long sequences only)
            let is_long_sequence = subject_len_nucl > 600_000;
            let mut stats_hsp_saved = 0usize;
            let mut stats_hsp_filtered_by_cutoff = 0usize;
            let mut stats_hsp_filtered_by_reeval = 0usize;
            let mut stats_hsp_filtered_by_hsp_test = 0usize;
            let mut stats_score_distribution: Vec<i32> = Vec::new();

            for (s_f_idx, s_frame) in s_frames.iter().enumerate() {
                let subject = &s_frame.aa_seq;
                let s_aa_len = s_frame.aa_len;
                if subject.len() < 5 {
                    continue;
                }

                // [C] scan_range[1] = subject->seq_ranges[0].left;
                // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                let mut scan_range: [i32; 3] = [0, 1, (subject.len() - 4) as i32];

                // [C] while (scan_range[1] <= scan_range[2])
                while scan_range[1] <= scan_range[2] {
                    let prev_scan_left = scan_range[1];
                    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                    let hits = s_blast_aa_scan_subject(
                        lookup_ref,
                        subject,
                        offset_pairs,
                        offset_array_size,
                        &mut scan_range,
                    );

                    if diag_enabled && hits > 0 {
                        diagnostics
                            .base
                            .kmer_matches
                            .fetch_add(hits as usize, AtomicOrdering::Relaxed);
                    }

                    if hits == 0 && scan_range[1] == prev_scan_left {
                        // Safety guard: with correct NCBI-sized offset arrays, this should not happen.
                        // If it does, breaking avoids an infinite loop.
                        break;
                    }

                    // [C] for (i = 0; i < hits; ++i)
                    // OPTIMIZATION: Use raw pointers to eliminate bounds checking in hot loop
                    let diag_ptr = diag_array.as_mut_ptr();
                    let offset_pairs_ptr = offset_pairs.as_ptr();
                    
                    for i in 0..hits as usize {
                        // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                        let pair = unsafe { &*offset_pairs_ptr.add(i) };
                        let query_offset = pair.q_off;
                        let subject_offset = pair.s_off;

                        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                        let diag_coord = ((query_offset - subject_offset) & diag_mask) as usize;
                        
                        // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
                        let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

                        // [C] if (diag_array[diag_coord].flag)
                        if diag_entry.flag != 0 {
                            // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
                            if subject_offset + diag_offset < diag_entry.last_hit {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_masked
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                continue;
                            }
                            // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            // [C] diag_array[diag_coord].flag = 0;
                            diag_entry.last_hit = subject_offset + diag_offset;
                            diag_entry.flag = 0;
                        }
                        // [C] else
                        else {
                            // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
                            let last_hit = diag_entry.last_hit - diag_offset;
                            // [C] diff = subject_offset - last_hit;
                            let diff = subject_offset - last_hit;

                            // [C] if (diff >= window)
                            if diff >= window {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_too_far
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }

                            // [C] if (diff < wordsize)
                            if diff < wordsize {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_overlap
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                continue;
                            }

                            // [C] curr_context = BSearchContextInfo(query_offset, query_info);
                            let ctx_idx = lookup_ref.get_context_idx(query_offset);
                            let ctx = unsafe { contexts_ref.get_unchecked(ctx_idx) };
                            let q_raw = (query_offset - ctx.frame_base) as usize;
                            // NCBI uses masked sequence for extension (same as neighbor-map mode)
                            let query = &ctx.aa_seq;

                            // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                            if query_offset - diff < ctx.frame_base {
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }

                            if diag_enabled {
                                diagnostics
                                    .base
                                    .seeds_second_hit_window
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                diagnostics
                                    .base
                                    .seeds_passed
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                            }

                            // [C] cutoffs = word_params->cutoffs + curr_context;
                            let cutoff = unsafe { *cutoff_scores.get_unchecked(ctx_idx).get_unchecked(s_f_idx) };
                            // [C] cutoffs->x_dropoff (per-context x_dropoff)
                            // Reference: aa_ungapped.c:579
                            let x_dropoff = unsafe { *x_dropoff_per_context.get_unchecked(ctx_idx) };

                            // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
                            //                                   last_hit + wordsize, subject_offset, query_offset, ...)
                            // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158
                            let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
                                extend_hit_two_hit(
                                    query,
                                    subject,
                                    (last_hit + wordsize) as usize, // s_left_off (end of first hit word)
                                    subject_offset as usize, // s_right_off (second hit start)
                                    q_raw as usize,          // q_right_off (second hit start, local)
                                    x_dropoff,               // x_dropoff (per-context, NCBI parity)
                                );

                            let hsp_q: i32 = hsp_q_u as i32;
                            let hsp_s: i32 = hsp_s_u as i32;
                            let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
                            let s_last_off: i32 = s_last_off_u as i32;

                            if diag_enabled {
                                diagnostics
                                    .base
                                    .ungapped_extensions
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                if right_extend {
                                    diagnostics
                                        .base
                                        .ungapped_two_hit_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                } else {
                                    diagnostics
                                        .base
                                        .ungapped_one_hit_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                if hsp_len > 0 {
                                    diagnostics
                                        .base
                                        .extension_total_length
                                        .fetch_add(hsp_len as usize, AtomicOrdering::Relaxed);
                                    atomic_max_usize(
                                        &diagnostics.base.extension_max_length,
                                        hsp_len as usize,
                                    );
                                }
                            }

                            // [C] if (right_extend)
                            if right_extend {
                                // [C] diag_array[diag_coord].flag = 1;
                                // [C] diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                                diag_entry.flag = 1;
                                diag_entry.last_hit =
                                    s_last_off - (wordsize - 1) + diag_offset;
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .mask_updates
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                            } else {
                                // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                diag_entry.last_hit = subject_offset + diag_offset;
                            }

                            // [C] if (score >= cutoffs->cutoff_score)
                            if is_long_sequence {
                                if score >= cutoff {
                                    stats_score_distribution.push(score);
                                } else {
                                    stats_hsp_filtered_by_cutoff += 1;
                                }
                            }
                            if score >= cutoff {
                                if diag_enabled {
                                    diagnostics
                                        .ungapped_only_hits
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                // NCBI: BlastSaveInitHsp equivalent
                                // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                                // Store HSP with absolute coordinates (before coordinate conversion)
                                //
                                // hsp_q is frame-relative coordinate (with sentinel, array index)
                                // frame_base is concatenated buffer start position (with sentinel, array index)
                                // Calculate absolute coordinate in concatenated buffer
                                // NCBI: ungapped_data->q_start = q_start (absolute, sentinel included)
                                let hsp_q_absolute = ctx.frame_base + (hsp_q as i32);
                                let hsp_qe_absolute = ctx.frame_base + ((hsp_q + hsp_len) as i32);
                                
                                init_hsps.push(InitHSP {
                                    q_start_absolute: hsp_q_absolute,
                                    q_end_absolute: hsp_qe_absolute,
                                    s_start: hsp_s,
                                    s_end: hsp_s + hsp_len,
                                    score,
                                    ctx_idx,
                                    s_f_idx,
                                    q_idx: ctx.q_idx,
                                    s_idx: s_idx as u32,
                                    q_frame: ctx.frame,
                                    s_frame: s_frame.frame,
                                    q_orig_len: ctx.orig_len,
                                    s_orig_len: s_len,
                                });
                            } else if diag_enabled {
                                diagnostics
                                    .ungapped_cutoff_failed
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                atomic_min_i32(&diagnostics.ungapped_cutoff_failed_min_score, score);
                                atomic_max_i32(&diagnostics.ungapped_cutoff_failed_max_score, score);
                            }
                        }
                    }
                }

                // [C] Blast_ExtendWordExit(ewp, subject->length);
                //
                // NCBI reference (verbatim):
                //   if (ewp->diag_table->offset >= INT4_MAX / 4) {
                //       ewp->diag_table->offset = ewp->diag_table->window;
                //       s_BlastDiagClear(ewp->diag_table);
                //   } else {
                //       ewp->diag_table->offset += subject_length + ewp->diag_table->window;
                //   }
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:167-173
                if diag_offset >= i32::MAX / 4 {
                    diag_offset = window;
                    // s_BlastDiagClear(): clear all diagonal state when offset risks overflow.
                    for d in diag_array.iter_mut() {
                        *d = DiagStruct::default();
                    }
                } else {
                    diag_offset += s_aa_len as i32 + window;
                }
            }

                // NCBI: BLAST_GetUngappedHSPList equivalent
                // Reference: blast_gapalign.c:4719-4775
                // Convert InitHSPs (with absolute coordinates) to UngappedHits (with context-relative coordinates)
                let ungapped_hits = get_ungapped_hsp_list(
                    init_hsps,
                    contexts_ref,
                    &s_frames,
                );
                
                // NCBI: Blast_HSPListReevaluateUngapped equivalent
                // Reference: blast_engine.c:1492-1497, blast_hits.c:2609-2737
                // Perform batch reevaluation on all HSPs
                let ungapped_hits = reevaluate_ungapped_hsp_list(
                    ungapped_hits,
                    contexts_ref,
                    &s_frames,
                    &cutoff_scores,
                    &args,
                );
                
                if !ungapped_hits.is_empty() {
                    if diag_enabled {
                        diagnostics
                            .base
                            .hsps_before_chain
                            .fetch_add(ungapped_hits.len(), AtomicOrdering::Relaxed);
                    }
                    // NCBI CalculateLinkHSPCutoffs parameters
                    let linking_params = LinkingParams {
                        avg_query_length,
                        subject_len_nucl,
                        cutoff_score_min,
                        scale_factor: 1.0, // Standard BLOSUM62
                    gap_decay_rate: 0.5, // BLAST_GAP_DECAY_RATE
                };
                // Build NCBI-style subject frame base offsets for sum-statistics linking.
                // In NCBI, HSP coords live in a concatenated translation buffer with sentinels.
                // LOSAT uses per-frame sequences; for linking we emulate absolute offsets by
                // concatenating frames in the same order as `generate_frames()`.
                let mut subject_frame_bases: Vec<i32> = Vec::with_capacity(s_frames.len());
                let mut base: i32 = 0;
                for f in &s_frames {
                    subject_frame_bases.push(base);
                    // NCBI concatenation shares the trailing NULLB sentinel between frames:
                    //   offset += length + 1;
                    // where `length` is the number of residues (excluding sentinels).
                    // Source: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
                    base += f.aa_seq.len() as i32 - 1;
                }

                // NCBI parity: s_BlastFindSmallestLambda selects Karlin params with smallest lambda
                // Reference: blast_parameters.c:92-112, 1012-1013
                // For tblastx, all contexts use kbp_ideal (same lambda), so this is equivalent
                // to using ungapped_params directly. We maintain this structure for parity.
                let context_params: Vec<KarlinParams> = contexts_ref.iter().map(|_| params.clone()).collect();
                let linking_params_for_cutoff = find_smallest_lambda_params(&context_params)
                    .unwrap_or_else(|| params.clone());

                // Output HSP saving statistics for long sequences
                if is_long_sequence && (stats_hsp_saved > 0 || stats_hsp_filtered_by_cutoff > 0 || stats_hsp_filtered_by_reeval > 0 || stats_hsp_filtered_by_hsp_test > 0) {
                    let total_attempted = stats_hsp_saved + stats_hsp_filtered_by_cutoff + stats_hsp_filtered_by_reeval + stats_hsp_filtered_by_hsp_test;
                    eprintln!("[DEBUG HSP_SAVING] subject_len_nucl={}, cutoff={}", subject_len_nucl, cutoff_score_min);
                    eprintln!("[DEBUG HSP_SAVING] total_attempted={}, saved={}, filtered_by_cutoff={}, filtered_by_reeval={}, filtered_by_hsp_test={}", 
                        total_attempted, stats_hsp_saved, stats_hsp_filtered_by_cutoff, stats_hsp_filtered_by_reeval, stats_hsp_filtered_by_hsp_test);
                    if !stats_score_distribution.is_empty() {
                        stats_score_distribution.sort();
                        let min_score = stats_score_distribution[0];
                        let max_score = stats_score_distribution[stats_score_distribution.len() - 1];
                        let median_score = stats_score_distribution[stats_score_distribution.len() / 2];
                        let low_score_count = stats_score_distribution.iter().filter(|&&s| s < 30).count();
                        eprintln!("[DEBUG HSP_SAVING] score_distribution: min={}, max={}, median={}, low_score(<30)={}/{} ({:.2}%)", 
                            min_score, max_score, median_score, low_score_count, stats_score_distribution.len(),
                            if stats_score_distribution.len() > 0 { (low_score_count as f64 / stats_score_distribution.len() as f64) * 100.0 } else { 0.0 });
                    }
                }
                
                // NCBI Parity: Pass pre-computed length_adjustment and eff_searchsp per context
                // These values are stored in query_info->contexts[ctx] in NCBI and referenced
                // by link_hsps.c for BLAST_SmallGapSumE/BLAST_LargeGapSumE calculations.
                let linked = apply_sum_stats_even_gap_linking(
                    ungapped_hits,
                    &linking_params_for_cutoff,
                    &linking_params,
                    contexts_ref,
                    &subject_frame_bases,
                    &length_adj_per_context,
                    &eff_searchsp_per_context,
                );
                if diag_enabled {
                    diagnostics
                        .base
                        .hsps_after_chain
                        .fetch_add(linked.len(), AtomicOrdering::Relaxed);
                }
                let mut final_hits: Vec<Hit> = Vec::new();
                for h in linked {
                    if h.e_value > evalue_threshold {
                        if diag_enabled {
                            diagnostics
                                .ungapped_evalue_failed
                                .fetch_add(1, AtomicOrdering::Relaxed);
                        }
                        continue;
                    }
                    
                    if diag_enabled {
                        diagnostics
                            .ungapped_evalue_passed
                            .fetch_add(1, AtomicOrdering::Relaxed);
                    }

                    let ctx = &contexts_ref[h.ctx_idx];
                    let s_frame = &s_frames[h.s_f_idx];

                    // Compute identity/mismatch on demand (final hits only)
                    // NCBI computes identities using the *unmasked* query sequence buffer
                    // (`sequence_nomask`), even when the working query sequence was masked
                    // (e.g., SEG) to influence seeding/extension.
                    //
                    // NCBI reference (verbatim):
                    //   const Uint1* query_nomask = query_blk->sequence_nomask +
                    //       query_info->contexts[context].query_offset;
                    let len = h.q_aa_end.saturating_sub(h.q_aa_start);
                    let q0 = h.q_aa_start + 1; // +1 for sentinel
                    let s0 = h.s_aa_start + 1; // +1 for sentinel
                    let q_seq_nomask: &[u8] = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
                    let mut matches = 0usize;
                    for k in 0..len {
                        if q_seq_nomask[q0 + k] == s_frame.aa_seq[s0 + k] {
                            matches += 1;
                        }
                    }
                    let mismatch = len.saturating_sub(matches);
                    let identity = if len > 0 {
                        (matches as f64 / len as f64) * 100.0
                    } else {
                        0.0
                    };

                    let bit = calc_bit_score(h.raw_score, &params);
                    let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, ctx.frame, ctx.orig_len);
                    let (s_start, s_end) =
                        convert_coords(h.s_aa_start, h.s_aa_end, s_frame.frame, s_len);

                    final_hits.push(Hit {
                        query_id: query_ids_ref[ctx.q_idx as usize].clone(),
                        subject_id: s_id.clone(),
                        identity,
                        length: len,
                        mismatch,
                        gapopen: 0,
                        q_start,
                        q_end,
                        s_start,
                        s_end,
                        e_value: h.e_value,
                        bit_score: bit,
                        q_idx: ctx.q_idx,
                        s_idx: h.s_idx,
                        raw_score: h.raw_score,
                    });
                }

                if !final_hits.is_empty() {
                    if diag_enabled {
                        diagnostics
                            .base
                            .hsps_after_overlap_filter
                            .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                        diagnostics
                            .output_from_ungapped
                            .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                    }
                    st.tx.send(final_hits).unwrap();
                }
            }
            bar.inc(1);
            st.diag_offset = diag_offset;
        },
    );

    // Close channel so the writer can exit.
    drop(tx);

    bar.finish();
    writer.join().unwrap()?;
    if diag_enabled {
        print_diagnostics_summary(&diagnostics);
    }
    Ok(())
}

/// NCBI: Blast_HSPListPurgeHSPsWithCommonEndpoints
/// Remove HSPs that share the same start OR end coordinates.
///
/// **NCBI Reference:** `blast_hits.c` lines 2454-2535
///
/// **When to use:** NCBI only calls this in the **gapped** path
/// (`blast_engine.c` line 545, inside `if (score_options->gapped_calculation)`).
/// For ungapped tblastx, this function should NOT be called for NCBI parity.
/// It is retained here for potential future gapped implementation.
///
/// **Important:** This must be called per-subject (BlastHSPList semantics).
/// Never call on a mixed collection containing multiple subjects.
///
/// Field mapping (NCBI → LOSAT):
/// - `hsp->context` → `UngappedHit.ctx_idx`
/// - `hsp->query.offset` → `UngappedHit.q_aa_start`
/// - `hsp->query.end` → `UngappedHit.q_aa_end`
/// - `hsp->subject.offset` → `UngappedHit.s_aa_start`
/// - `hsp->subject.end` → `UngappedHit.s_aa_end`
/// - `hsp->subject.frame` → `UngappedHit.s_frame`
/// - `hsp->score` → `UngappedHit.raw_score`
///
/// **Note:** `subject.frame` is NOT in the sort comparators, but IS used
/// in duplicate detection (NCBI lines 2487, 2513).
#[allow(dead_code)]
fn purge_hsps_with_common_endpoints(hits: &mut Vec<super::chaining::UngappedHit>) {
    if hits.len() <= 1 {
        return;
    }

    // =========================================================================
    // Phase 1: Remove HSPs with common start points
    // =========================================================================
    // NCBI s_QueryOffsetCompareHSPs (blast_hits.c:2267-2321):
    // Sort order: context ASC → query.offset ASC → subject.offset ASC →
    //             score DESC → query.end ASC → subject.end ASC
    // ```c
    // if (h1->context < h2->context) return -1;
    // if (h1->context > h2->context) return 1;
    // if (h1->query.offset < h2->query.offset) return -1;
    // if (h1->query.offset > h2->query.offset) return 1;
    // if (h1->subject.offset < h2->subject.offset) return -1;
    // if (h1->subject.offset > h2->subject.offset) return 1;
    // if (h1->score < h2->score) return 1;   // DESC
    // if (h1->score > h2->score) return -1;  // DESC
    // if (h1->query.end < h2->query.end) return 1;   // shorter range first
    // if (h1->query.end > h2->query.end) return -1;
    // if (h1->subject.end < h2->subject.end) return 1;
    // if (h1->subject.end > h2->subject.end) return -1;
    // ```
    hits.sort_by(|a, b| {
        a.ctx_idx
            .cmp(&b.ctx_idx)
            .then_with(|| a.q_aa_start.cmp(&b.q_aa_start))
            .then_with(|| a.s_aa_start.cmp(&b.s_aa_start))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // DESC: higher score first
            .then_with(|| b.q_aa_end.cmp(&a.q_aa_end)) // DESC: larger end first (NCBI line 2310-2313)
            .then_with(|| b.s_aa_end.cmp(&a.s_aa_end)) // DESC: larger end first (NCBI line 2315-2318)
    });

    // NCBI duplicate detection (blast_hits.c:2482-2487):
    // ```c
    // while (i+j < hsp_count &&
    //        hsp_array[i] && hsp_array[i+j] &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.offset == hsp_array[i+j]->query.offset &&
    //        hsp_array[i]->subject.offset == hsp_array[i+j]->subject.offset &&
    //        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame)
    // ```
    let mut write_idx = 0;
    let mut i = 0;
    while i < hits.len() {
        // Keep the first HSP of each duplicate group (highest score due to sort)
        hits.swap(write_idx, i);
        let kept_ctx = hits[write_idx].ctx_idx;
        let kept_q_off = hits[write_idx].q_aa_start;
        let kept_s_off = hits[write_idx].s_aa_start;
        let kept_s_frame = hits[write_idx].s_frame;
        write_idx += 1;

        // Skip duplicates with same (context, query.offset, subject.offset, subject.frame)
        let mut j = 1;
        while i + j < hits.len() {
            let h = &hits[i + j];
            if h.ctx_idx == kept_ctx
                && h.q_aa_start == kept_q_off
                && h.s_aa_start == kept_s_off
                && h.s_frame == kept_s_frame
            {
                j += 1;
            } else {
                break;
            }
        }
        i += j;
    }
    hits.truncate(write_idx);

    // =========================================================================
    // Phase 2: Remove HSPs with common end points
    // =========================================================================
    // NCBI s_QueryEndCompareHSPs (blast_hits.c:2332-2387):
    // Sort order: context ASC → query.end ASC → subject.end ASC →
    //             score DESC → query.offset DESC → subject.offset DESC
    // ```c
    // if (h1->context < h2->context) return -1;
    // if (h1->context > h2->context) return 1;
    // if (h1->query.end < h2->query.end) return -1;
    // if (h1->query.end > h2->query.end) return 1;
    // if (h1->subject.end < h2->subject.end) return -1;
    // if (h1->subject.end > h2->subject.end) return 1;
    // if (h1->score < h2->score) return 1;   // DESC
    // if (h1->score > h2->score) return -1;  // DESC
    // if (h1->query.offset < h2->query.offset) return 1;   // shorter range first (larger offset)
    // if (h1->query.offset > h2->query.offset) return -1;
    // if (h1->subject.offset < h2->subject.offset) return 1;
    // if (h1->subject.offset > h2->subject.offset) return -1;
    // ```
    hits.sort_by(|a, b| {
        a.ctx_idx
            .cmp(&b.ctx_idx)
            .then_with(|| a.q_aa_end.cmp(&b.q_aa_end))
            .then_with(|| a.s_aa_end.cmp(&b.s_aa_end))
            .then_with(|| b.raw_score.cmp(&a.raw_score)) // DESC: higher score first
            .then_with(|| b.q_aa_start.cmp(&a.q_aa_start)) // DESC: larger offset first (NCBI line 2376-2379)
            .then_with(|| b.s_aa_start.cmp(&a.s_aa_start)) // DESC: larger offset first (NCBI line 2381-2384)
    });

    // NCBI duplicate detection (blast_hits.c:2508-2513):
    // ```c
    // while (i+j < hsp_count &&
    //        hsp_array[i] && hsp_array[i+j] &&
    //        hsp_array[i]->context == hsp_array[i+j]->context &&
    //        hsp_array[i]->query.end == hsp_array[i+j]->query.end &&
    //        hsp_array[i]->subject.end == hsp_array[i+j]->subject.end &&
    //        hsp_array[i]->subject.frame == hsp_array[i+j]->subject.frame)
    // ```
    write_idx = 0;
    i = 0;
    while i < hits.len() {
        hits.swap(write_idx, i);
        let kept_ctx = hits[write_idx].ctx_idx;
        let kept_q_end = hits[write_idx].q_aa_end;
        let kept_s_end = hits[write_idx].s_aa_end;
        let kept_s_frame = hits[write_idx].s_frame;
        write_idx += 1;

        // Skip duplicates with same (context, query.end, subject.end, subject.frame)
        let mut j = 1;
        while i + j < hits.len() {
            let h = &hits[i + j];
            if h.ctx_idx == kept_ctx
                && h.q_aa_end == kept_q_end
                && h.s_aa_end == kept_s_end
                && h.s_frame == kept_s_frame
            {
                j += 1;
            } else {
                break;
            }
        }
        i += j;
    }
    hits.truncate(write_idx);
}

/// Run TBLASTX with pre-computed neighbor map using subject-side indexing.
/// This approach:
/// 1. Index query k-mers: query_lookup[kmer] = [(q_idx, f_idx, pos), ...]
/// 2. Pre-compute neighbor relationships once
/// 3. For each subject k-mer, find all matching query positions via expanded lookup
/// 4. Apply NCBI-style sum_stats_even_gap_linking for HSP merging
///
/// Two-hit/diag logic is ported verbatim from NCBI aa_ungapped.c:s_BlastAaWordFinder_TwoHit
fn run_with_neighbor_map(args: TblastxArgs) -> Result<()> {
    use super::lookup::LOOKUP_ALPHABET_SIZE;

    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    // NCBI: window = diag->window
    let window: i32 = args.window_size as i32;
    // NCBI: wordsize = lookup->word_length (always 3 for protein/tblastx)
    let wordsize: i32 = 3;
    let threshold = args.threshold;

    // NCBI BLAST computes x_dropoff per-context using kbp[context]->Lambda:
    //   p->cutoffs[context].x_dropoff_init =
    //       (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // Reference: blast_parameters.c:219-221
    //
    // For tblastx, NCBI uses kbp_ideal->Lambda for all contexts (blast_stat.c:2796-2797).
    // x_dropoff_per_context is populated per-subject after counting total_contexts.
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    eprintln!("Reading queries (neighbor map mode)...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries_raw
        .iter()
        .map(|r| r.id().split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();

    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), &query_code))
        .collect();

    // Apply SEG masking if enabled
    if args.seg {
        let seg_masker = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        const X_MASK_NCBISTDAA: u8 = 21;
        for frames in query_frames.iter_mut() {
            for frame in frames.iter_mut() {
                if frame.aa_seq.len() >= 3 {
                    for m in seg_masker.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
                        frame.seg_masks.push((m.start, m.end));
                    }
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1);
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!("Building neighbor lookup...");
    let neighbor_lookup = NeighborLookup::build(&query_frames, threshold, &ungapped_params_for_xdrop);
    
    // Use compressed neighbor index: no expanded_lookup pre-computation
    // Instead, resolve neighbors on-the-fly during scan
    eprintln!("Using compressed neighbor index (no expanded_lookup)...");

    eprintln!("Reading subjects...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    // NCBI BLAST Karlin params for TBLASTX (ungapped-only algorithm):
    // 
    // TBLASTX is explicitly ungapped-only (blast_options.c line 869-873):
    //   "Gapped search is not allowed for tblastx"
    //
    // For bit score and E-value calculation, NCBI uses sbp->kbp (ungapped):
    //   blast_hits.c line 1833: kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    //   blast_hits.c line 1918: same pattern in Blast_HSPListGetBitScores
    //
    // For cutoff score search space calculation, NCBI uses gapped params:
    //   blast_parameters.c uses kbp_gap for eff_searchsp in cutoff calculation
    //
    // BLOSUM62 ungapped: lambda=0.3176, K=0.134
    // BLOSUM62 gapped:   lambda=0.267,  K=0.041
    let ungapped_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    // Gapped params for cutoff score search space calculation
    let gapped_params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    // Use UNGAPPED params for bit score and E-value (NCBI parity for tblastx)
    let params = ungapped_params.clone();

    // Compute NCBI-style average query length for linking cutoffs
    // Reference: blast_parameters.c:CalculateLinkHSPCutoffs (lines 1023-1026)
    let query_nucl_lengths: Vec<usize> = queries_raw.iter().map(|r| r.seq().len()).collect();
    let avg_query_length = compute_avg_query_length_ncbi(&query_nucl_lengths);

    eprintln!(
        "Searching {} queries vs {} subjects (neighbor map mode, avg_query_length={})...",
        queries_raw.len(),
        subjects_raw.len(),
        avg_query_length
    );

    let bar = ProgressBar::new(subjects_raw.len() as u64 * 6);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    let evalue_threshold = args.evalue;
    let query_frames_ref = &query_frames;
    let query_ids_ref = &query_ids;
    let neighbor_map_ref = &neighbor_lookup.neighbor_map.map;
    let query_lookup_ref = &neighbor_lookup.query_lookup;
    let query_contexts_ref = &neighbor_lookup.contexts;
    let params_ref = &params;
    let ungapped_params_ref = &ungapped_params;
    let _gapped_params_ref = &gapped_params;  // Unused - tblastx uses ungapped params

    // Collect UngappedHit for sum_stats_linking
    // Key: (q_idx, s_idx)
    let all_ungapped: std::sync::Mutex<Vec<UngappedHit>> = std::sync::Mutex::new(Vec::new());

    // Process each subject
    for (s_idx, s_rec) in subjects_raw.iter().enumerate() {
        let s_id = s_rec
            .id()
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();
        let s_len = s_rec.seq().len();
        let s_frames = generate_frames(s_rec.seq(), &db_code);
        
        // Precompute per-context cutoff scores using NCBI BLAST algorithm.
        // Reference: blast_parameters.c BlastInitialWordParametersUpdate
        //
        // NCBI cutoff calculation for tblastx ungapped path:
        // 1. gap_trigger from ungapped params (kbp_std)
        // 2. cutoff_score_max from BlastHitSavingParametersNew (uses user's E-value)
        // 3. Per-subject update: cutoff_score_for_update_tblastx with CUTOFF_E_TBLASTX=1e-300
        // 4. Final cutoff = MIN(update_cutoff, gap_trigger, cutoff_score_max)
        //
        // For -subject mode, subject_len_nucl is used (not per-frame AA length).
        let subject_len_nucl = s_len as i64;
        let total_contexts: usize = query_frames.iter().map(|f| f.len()).sum();
        let mut cutoff_scores: Vec<Vec<i32>> = vec![vec![0; s_frames.len()]; total_contexts];
        
        // =======================================================================
        // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
        // =======================================================================
        // NCBI stores these in query_info->contexts[ctx].length_adjustment and
        // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths (blast_setup.c:700-850).
        // We precompute them here and pass to sum_stats_linking for NCBI parity.
        // Reference: blast_setup.c:846-847
        //   query_info->contexts[index].eff_searchsp = effective_search_space;
        //   query_info->contexts[index].length_adjustment = length_adjustment;
        let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
        let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
        
        // NCBI BLAST: word_params->cutoffs[context].x_dropoff_init
        // Compute per-context x_dropoff using kbp[context]->Lambda.
        // Reference: blast_parameters.c:219-221
        //
        // For tblastx, all contexts use kbp_ideal (BLOSUM62 ungapped Lambda = 0.3176),
        // so x_dropoff = ceil(7.0 * LN2 / 0.3176) = 16 for all contexts.
        // We maintain per-context structure for NCBI parity.
        let x_dropoff_per_context: Vec<i32> = (0..total_contexts)
            .map(|_| x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ungapped_params_for_xdrop, 1.0))
            .collect();
        
        // Precompute gap_trigger once (same for all contexts with same params)
        // Reference: blast_parameters.c:343-344
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ungapped_params_ref);
        
        let mut ctx_idx = 0;
        for q_frames in query_frames.iter() {
            for q_frame in q_frames.iter() {
                // NCBI: query_length = query_info->contexts[context].query_length
                let query_len_aa = q_frame.aa_len as i64;
                
                // =======================================================================
                // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                // length_adjustment and eff_searchsp from a single source of truth.
                // This mirrors BLAST_CalcEffLengths which computes and stores both values.
                // Reference: blast_setup.c:821-847
                // =======================================================================
                let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                    query_len_aa,
                    subject_len_nucl,
                    ungapped_params_ref,  // tblastx uses ungapped params (kbp_gap is NULL)
                );
                let eff_searchsp = eff_lengths.eff_searchsp;
                length_adj_per_context.push(eff_lengths.length_adjustment);
                eff_searchsp_per_context.push(eff_searchsp);
                
                // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
                // This uses the effective search space WITH length adjustment
                // Reference: blast_parameters.c:942-946
                let cutoff_score_max = cutoff_score_max_for_tblastx(
                    eff_searchsp,
                    evalue_threshold,  // User's E-value (typically 10.0)
                    ungapped_params_ref,
                );
                
                // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
                // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
                // Reference: blast_parameters.c:348-374
                let cutoff = cutoff_score_for_update_tblastx(
                    query_len_aa,
                    subject_len_nucl,  // NUCLEOTIDE length, NOT divided by 3!
                    gap_trigger,
                    cutoff_score_max,
                    BLAST_GAP_DECAY_RATE,  // 0.5
                    ungapped_params_ref,
                    1.0,  // scale_factor (standard BLOSUM62)
                );
                
                // All subject frames use the same cutoff (NCBI behavior)
                for sf_idx in 0..s_frames.len() {
                    cutoff_scores[ctx_idx][sf_idx] = cutoff;
                }
                ctx_idx += 1;
            }
        }
        let cutoff_scores_ref = &cutoff_scores;
        let x_dropoff_per_context_ref = &x_dropoff_per_context;
        let length_adj_per_context_ref = &length_adj_per_context;
        let eff_searchsp_per_context_ref = &eff_searchsp_per_context;

        // Process each subject frame in parallel
        let frame_hits: Vec<UngappedHit> = s_frames
            .par_iter()
            .enumerate()
            .flat_map(|(s_f_idx, s_frame)| {
                let s_aa = &s_frame.aa_seq;
                if s_aa.len() < 5 {
                    bar.inc(1);
                    return Vec::new();
                }

                // NCBI: BlastSaveInitHsp equivalent - store initial HSPs with absolute coordinates
                // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                let s_aa_len = s_aa.len();
                
                // Count total query contexts for Vec-based diagonal tracking
                let total_q_contexts: usize = query_frames_ref.iter().map(|f| f.len()).sum();
                let max_q_aa_len = query_frames_ref.iter()
                    .flat_map(|frames| frames.iter().map(|f| f.aa_seq.len()))
                    .max()
                    .unwrap_or(1);
                
                // Compute diag_array dimensions
                let diag_array_len = s_aa_len + max_q_aa_len + 1;
                let diag_offset = max_q_aa_len as i32;
                
                // Build context index mapping: (q_idx, q_f_idx) -> flat index
                let mut ctx_base: Vec<usize> = Vec::with_capacity(query_frames_ref.len() + 1);
                ctx_base.push(0);
                for frames in query_frames_ref.iter() {
                    ctx_base.push(ctx_base.last().unwrap() + frames.len());
                }
                
                // NCBI-style DiagStruct arrays: [ctx_flat_idx][diag] = DiagStruct { last_hit, flag }
                // NCBI: last_hit stores (subject_offset + diag_offset) to avoid subtraction in hot path
                // NCBI: flag=1 means we just extended, flag=0 means ready for two-hit check
                let mut diag_array: Vec<Vec<DiagStruct>> = vec![vec![DiagStruct::default(); diag_array_len]; total_q_contexts];
                
                // NCBI aa_ungapped.c:502-610 - main scan loop
                for s_pos in 1..=(s_aa_len.saturating_sub(4)) {
                    let c0 = unsafe { *s_aa.get_unchecked(s_pos) } as usize;
                    let c1 = unsafe { *s_aa.get_unchecked(s_pos + 1) } as usize;
                    let c2 = unsafe { *s_aa.get_unchecked(s_pos + 2) } as usize;

                    if c0 >= LOOKUP_ALPHABET_SIZE || c1 >= LOOKUP_ALPHABET_SIZE || c2 >= LOOKUP_ALPHABET_SIZE {
                        continue;
                    }

                    let s_kmer = encode_kmer(c0, c1, c2);
                    
                    // Compressed neighbor index: resolve neighbors on-the-fly
                    let neighbor_kmers = &neighbor_map_ref[s_kmer];
                    if neighbor_kmers.is_empty() {
                        continue;
                    }
                    
                    // subject_offset in NCBI terms (0-based position in subject AA sequence)
                    let subject_offset: i32 = (s_pos - 1) as i32;
                    
                    // For each neighbor query k-mer, get query positions
                    for &neighbor_kmer in neighbor_kmers {
                        let query_hits = &query_lookup_ref[neighbor_kmer as usize];
                        for &(q_idx, q_f_idx, q_pos) in query_hits {
                            let ctx_flat = ctx_base[q_idx as usize] + q_f_idx as usize;
                            
                            // query_offset in NCBI terms (0-based position in query AA sequence)
                            let query_offset: i32 = q_pos as i32;
                            
                            // NCBI: diag_coord = (query_offset - subject_offset) & diag_mask
                            // We use a Vec so we compute a direct index instead of masking
                            let diag_coord = (query_offset - subject_offset + diag_offset) as usize;
                            
                            if diag_coord >= diag_array_len {
                                continue;
                            }
                            
                            let diag_entry = &mut diag_array[ctx_flat][diag_coord];
                            
                            // NCBI aa_ungapped.c:519-530 - flag==1 block
                            // "If the reset bit is set, an extension just happened."
                            if diag_entry.flag != 0 {
                                // "If we've already extended past this hit, skip it."
                                if subject_offset + diag_offset < diag_entry.last_hit {
                                    continue;
                                }
                                // "Otherwise, start a new hit." - reset and CONTINUE (don't extend this hit)
                                diag_entry.last_hit = subject_offset + diag_offset;
                                diag_entry.flag = 0;
                                continue; // NCBI: after reset, move to next hit, don't extend
                            }
                            
                            // NCBI aa_ungapped.c:533-606 - flag==0 block
                            // "If the reset bit is cleared, try to start an extension."
                            
                            // NCBI: last_hit = diag_array[diag_coord].last_hit - diag_offset
                            let last_hit = diag_entry.last_hit - diag_offset;
                            // NCBI: diff = subject_offset - last_hit
                            let diff = subject_offset - last_hit;
                            
                            // NCBI aa_ungapped.c:538-544: "if (diff >= window)"
                            if diff >= window {
                                // "We are beyond the window for this diagonal; start a new hit"
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:549-551: "if (diff < wordsize)"
                            if diff < wordsize {
                                // "If the difference is less than the wordsize (i.e. last hit and this hit overlap), give up"
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:560-573: Check if last hit is in current query context
                            // "if (query_offset - diff < query_info->contexts[curr_context].query_offset)"
                            // In neighbor-map mode, each (q_idx, q_f_idx) is its own context with frame_base=0
                            // So we check: query_offset - diff < 0
                            if query_offset - diff < 0 {
                                // "there was no last hit for this diagonal; start a new hit"
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // PASSED all two-hit checks - now extend
                            let q_frame = &query_frames_ref[q_idx as usize][q_f_idx as usize];
                            let q_aa = &q_frame.aa_seq;
                            
                            // Convert to raw positions (+1 for sentinel)
                            let q_raw = (query_offset + 1) as usize;
                            let s_raw = (subject_offset + 1) as usize;
                            
                            if q_raw + 2 >= q_aa.len() || s_raw + 2 >= s_aa.len() {
                                continue;
                            }
                            
                            // [C] cutoffs->x_dropoff (per-context x_dropoff)
                            // Reference: aa_ungapped.c:579
                            let x_dropoff = x_dropoff_per_context_ref[ctx_flat];
                            
                            // NCBI aa_ungapped.c:576-583: s_BlastAaExtendTwoHit
                            // s_left_off = last_hit + wordsize (end of first hit word)
                            // s_right_off = subject_offset (second hit start)
                            // q_right_off = query_offset (second hit start)
                            // last_hit = subject_offset which is 0-indexed, so add +1 for array index
                            // (consistent with s_raw and q_raw which also have +1)
                            let s_left_off = (last_hit + wordsize + 1) as usize;
                            let (hsp_q, hsp_qe, hsp_s, _hsp_se, score, right_extend, s_last_off) =
                                extend_hit_two_hit(
                                    q_aa,
                                    s_aa,
                                    s_left_off,      // s_left_off (end of first hit word, raw coords)
                                    s_raw,           // s_right_off (second hit start, raw coords)
                                    q_raw,           // q_right_off (second hit start, raw coords)
                                    x_dropoff,       // x_dropoff (per-context, NCBI parity)
                                );
                            
                            let hsp_len = hsp_qe.saturating_sub(hsp_q);
                            
                            // NCBI aa_ungapped.c:596-606: Update diag_array after extension
                            if right_extend {
                                // "If an extension to the right happened, reset the last hit so that
                                //  future hits to this diagonal must start over."
                                diag_entry.flag = 1;
                                diag_entry.last_hit = (s_last_off as i32) - (wordsize - 1) + diag_offset;
                            } else {
                                // "Otherwise, make the present hit into the previous hit for this diagonal"
                                diag_entry.last_hit = subject_offset + diag_offset;
                            }
                            
                            // NCBI aa_ungapped.c:588-591: Check score threshold
                            // cutoffs = word_params->cutoffs + curr_context
                            // if (score >= cutoffs->cutoff_score) ...
                            let cutoff = cutoff_scores_ref[ctx_flat][s_f_idx];
                            if score < cutoff {
                                continue;
                            }
                            
                            // NCBI: BlastSaveInitHsp equivalent
                            // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                            // Store HSP with absolute coordinates (before coordinate conversion)
                            //
                            // In neighbor-map mode, each query frame is independent with frame_base=0
                            // hsp_q is frame-relative coordinate (with sentinel)
                            // Calculate absolute coordinate in concatenated buffer
                            let frame_base = 0i32;  // Each query frame is independent in neighbor-map mode
                            let hsp_q_absolute = frame_base + (hsp_q as i32) - 1;  // -1 for sentinel
                            let hsp_qe_absolute = frame_base + ((hsp_q + hsp_len) as i32) - 1;  // -1 for sentinel
                            
                            init_hsps.push(InitHSP {
                                q_start_absolute: hsp_q_absolute,
                                q_end_absolute: hsp_qe_absolute,
                                s_start: hsp_s as i32,
                                s_end: (hsp_s + hsp_len) as i32,
                                score,
                                ctx_idx: ctx_flat,
                                s_f_idx,
                                q_idx,
                                s_idx: s_idx as u32,
                                q_frame: q_frame.frame,
                                s_frame: s_frame.frame,
                                q_orig_len: q_frame.orig_len,
                                s_orig_len: s_len,
                            });
                        } // end for query_hits
                    } // end for neighbor_kmers
                } // end for s_pos

                // NCBI: BLAST_GetUngappedHSPList equivalent
                // Reference: blast_gapalign.c:4719-4775
                // Convert InitHSPs (with absolute coordinates) to UngappedHits (with context-relative coordinates)
                // and perform batch reevaluation
                //
                // In neighbor-map mode, each query frame is independent with frame_base=0
                // Create temporary contexts for get_ungapped_hsp_list
                let mut temp_contexts = Vec::new();
                for (q_idx, frames) in query_frames_ref.iter().enumerate() {
                    for (f_idx, frame) in frames.iter().enumerate() {
                        temp_contexts.push(super::lookup::QueryContext {
                            q_idx: q_idx as u32,
                            f_idx: f_idx as u8,
                            frame: frame.frame,
                            aa_seq: frame.aa_seq.clone(),
                            aa_seq_nomask: frame.aa_seq_nomask.clone(),
                            aa_len: frame.aa_len,
                            orig_len: frame.orig_len,
                            frame_base: 0,  // Each query frame is independent
                            karlin_params: ungapped_params_ref.clone(),
                        });
                    }
                }
                // NCBI: BLAST_GetUngappedHSPList equivalent
                // Convert InitHSPs (with absolute coordinates) to UngappedHits (with context-relative coordinates)
                let local_hits = get_ungapped_hsp_list(
                    init_hsps,
                    &temp_contexts,
                    &s_frames,
                );
                
                // NCBI: Blast_HSPListReevaluateUngapped equivalent
                // Perform batch reevaluation on all HSPs
                let local_hits = reevaluate_ungapped_hsp_list(
                    local_hits,
                    &temp_contexts,
                    &s_frames,
                    cutoff_scores_ref,
                    &args,
                );

                bar.inc(1);
                local_hits
            })
            .collect();

        // Add to global collection
        all_ungapped.lock().unwrap().extend(frame_hits);
    }

    bar.finish();

    let all_ungapped = all_ungapped.into_inner().unwrap();
    eprintln!("=== Stage Counters ===");
    eprintln!("[1] Raw ungapped hits (after extension): {}", all_ungapped.len());

    // NOTE: NCBI does NOT call Blast_HSPListPurgeHSPsWithCommonEndpoints in the
    // ungapped tblastx path. The purge is only in the gapped path:
    //   blast_engine.c line 545: inside `if (score_options->gapped_calculation)`
    // Since tblastx is ungapped, we skip purge for NCBI parity.
    //
    // Additionally, if purge were needed, it must be per-subject (BlastHSPList),
    // not on the entire all_ungapped collection which mixes multiple subjects.

    // Apply sum-stats linking per subject (NCBI CalculateLinkHSPCutoffs is subject-specific)
    // Reference: blast_parameters.c:CalculateLinkHSPCutoffs
    // Step 1: Group hits by s_idx
    let mut hits_by_subject: rustc_hash::FxHashMap<u32, Vec<UngappedHit>> = rustc_hash::FxHashMap::default();
    for hit in all_ungapped {
        hits_by_subject.entry(hit.s_idx).or_default().push(hit);
    }
    
    // Step 2: Apply linking per subject with proper cutoffs
    let linked_hits: Vec<UngappedHit> = hits_by_subject
        .into_par_iter()
        .flat_map(|(s_idx, subject_hits)| {
            // Get subject length for this subject
            let subject_len_nucl = subjects_raw[s_idx as usize].seq().len() as i64;
            
            // =======================================================================
            // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
            // =======================================================================
            // NCBI stores these in query_info->contexts[ctx].length_adjustment and
            // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths (blast_setup.c:700-850).
            // We precompute them here for each subject and pass to sum_stats_linking.
            let total_contexts: usize = query_frames_ref.iter().map(|f| f.len()).sum();
            let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
            let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(total_contexts);
            
            // Compute cutoff_score_min as minimum across all query contexts
            // using the NCBI BlastInitialWordParametersUpdate logic
            // Reference: blast_parameters.c:348-374, 401-403
            let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ungapped_params_ref);
            let mut cutoff_score_min = i32::MAX;
            for frames in query_frames_ref.iter() {
                for frame in frames.iter() {
                    let query_len_aa = frame.aa_len as i64;
                    
                    // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
                    // length_adjustment and eff_searchsp from a single source of truth.
                    // Reference: blast_setup.c:821-847
                    let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                        query_len_aa,
                        subject_len_nucl,
                        ungapped_params_ref,  // tblastx uses ungapped params (kbp_gap is NULL)
                    );
                    let eff_searchsp = eff_lengths.eff_searchsp;
                    length_adj_per_context.push(eff_lengths.length_adjustment);
                    eff_searchsp_per_context.push(eff_searchsp);
                    
                    // Step 1: cutoff_score_max from BlastHitSavingParametersNew
                    let cutoff_score_max = cutoff_score_max_for_tblastx(
                        eff_searchsp,
                        evalue_threshold,
                        ungapped_params_ref,
                    );
                    
                    // Step 2: Per-subject cutoff from BlastInitialWordParametersUpdate
                    let cutoff = cutoff_score_for_update_tblastx(
                        query_len_aa,
                        subject_len_nucl,
                        gap_trigger,
                        cutoff_score_max,
                        BLAST_GAP_DECAY_RATE,
                        ungapped_params_ref,
                        1.0,
                    );
                    cutoff_score_min = cutoff_score_min.min(cutoff);
                }
            }
            if cutoff_score_min == i32::MAX {
                cutoff_score_min = 0;
            }
            
            // Create NCBI-style linking params for this subject
            let linking_params = LinkingParams {
                avg_query_length,
                subject_len_nucl,
                cutoff_score_min,
                scale_factor: 1.0,
                gap_decay_rate: 0.5,
            };

            // Subject frame bases in the concatenated translation buffer order.
            // We match `generate_frames()` ordering (forward 1,2,3 then reverse -1,-2,-3)
            // and each frame has its own leading+trailing sentinel in LOSAT's current model.
            let subject_len = subject_len_nucl as usize;
            let mut subject_frame_bases: Vec<i32> = Vec::with_capacity(6);
            let mut base: i32 = 0;
            for shift in 0..3usize {
                if shift + 3 <= subject_len {
                    let aa_len = ((subject_len - shift) / 3) as i32;
                    subject_frame_bases.push(base);
                    // NCBI concatenation shares the boundary NULLB between frames:
                    // increment is (aa_len + 1), not (aa_len + 2).
                    base += aa_len + 1;
                }
            }
            for shift in 0..3usize {
                if shift + 3 <= subject_len {
                    let aa_len = ((subject_len - shift) / 3) as i32;
                    subject_frame_bases.push(base);
                    base += aa_len + 1;
                }
            }

            // NCBI parity: s_BlastFindSmallestLambda selects Karlin params with smallest lambda
            // Reference: blast_parameters.c:92-112, 1012-1013
            // For tblastx, all contexts use kbp_ideal (same lambda), so this is equivalent
            // to using params_ref directly. We maintain this structure for parity.
            let context_params: Vec<KarlinParams> = query_contexts_ref.iter().map(|_| params_ref.clone()).collect();
            let linking_params_for_cutoff = find_smallest_lambda_params(&context_params)
                .unwrap_or_else(|| params_ref.clone());

            // NCBI Parity: Pass pre-computed length_adjustment and eff_searchsp per context
            // These values are stored in query_info->contexts[ctx] in NCBI and referenced
            // by link_hsps.c for BLAST_SmallGapSumE/BLAST_LargeGapSumE calculations.
            apply_sum_stats_even_gap_linking(
                subject_hits,
                &linking_params_for_cutoff,
                &linking_params,
                query_contexts_ref,
                &subject_frame_bases,
                &length_adj_per_context,
                &eff_searchsp_per_context,
            )
        })
        .collect();
    eprintln!("[2] After sum_stats_linking: {} hits", linked_hits.len());
    
    // Count E-value distribution before final filter
    let mut e_0 = 0usize;      // E <= 0 (0.0)
    let mut e_tiny = 0usize;   // 0 < E <= 1e-50
    let mut e_small = 0usize;  // 1e-50 < E <= 1e-10
    let mut e_med = 0usize;    // 1e-10 < E <= 10
    let mut e_large = 0usize;  // E > 10
    let mut e_inf = 0usize;    // E == INF
    for h in &linked_hits {
        if h.e_value == f64::INFINITY {
            e_inf += 1;
        } else if h.e_value <= 0.0 {
            e_0 += 1;
        } else if h.e_value <= 1e-50 {
            e_tiny += 1;
        } else if h.e_value <= 1e-10 {
            e_small += 1;
        } else if h.e_value <= 10.0 {
            e_med += 1;
        } else {
            e_large += 1;
        }
    }
    eprintln!("[3] E-value distribution (before final filter):");
    eprintln!("    E=0: {}, E<=1e-50: {}, E<=1e-10: {}, E<=10: {}, E>10: {}, E=INF: {}",
        e_0, e_tiny, e_small, e_med, e_large, e_inf);

    // Pre-generate subject frames to avoid redundant translation
    eprintln!("Pre-generating subject frames for identity calculation...");
    let subject_frames_cache: Vec<Vec<QueryFrame>> = subjects_raw
        .iter()
        .map(|s| generate_frames(s.seq(), &db_code))
        .collect();
    
    // Convert UngappedHit to Hit
    let mut final_hits: Vec<Hit> = Vec::with_capacity(linked_hits.len());
    
    for h in linked_hits {
        if h.e_value > evalue_threshold {
            continue;
        }
        
        let q_id = &query_ids_ref[h.q_idx as usize];
        let s_id = subjects_raw[h.s_idx as usize]
            .id()
            .split_whitespace()
            .next()
            .unwrap_or("unknown");
        
        // `ctx_idx` is the global query context index (NCBI `hsp->context` equivalent).
        // Use the pre-built `QueryContext` buffer to access the (possibly unmasked) query AA sequence.
        let q_ctx = &query_contexts_ref[h.ctx_idx];
        let q_aa = q_ctx.aa_seq_nomask.as_deref().unwrap_or(&q_ctx.aa_seq);
        
        // Get subject frame from cache (no redundant translation)
        let s_frame_obj = &subject_frames_cache[h.s_idx as usize][h.s_f_idx];
        let s_aa = &s_frame_obj.aa_seq;
        
        let len = h.q_aa_end.saturating_sub(h.q_aa_start);
        
        // Calculate identity using raw positions (+1 for sentinel)
        let q0 = h.q_aa_start + 1;
        let s0 = h.s_aa_start + 1;
        let mut matches = 0;
        for k in 0..len {
            if q0 + k < q_aa.len() && s0 + k < s_aa.len() && q_aa[q0 + k] == s_aa[s0 + k] {
                matches += 1;
            }
        }
        let identity = if len > 0 {
            (matches as f64 / len as f64) * 100.0
        } else {
            0.0
        };

        let bit_score = calc_bit_score(h.raw_score, &params);
        let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, h.q_frame, h.q_orig_len);
        let (s_start, s_end) = convert_coords(h.s_aa_start, h.s_aa_end, h.s_frame, h.s_orig_len);

        final_hits.push(Hit {
            query_id: q_id.clone(),
            subject_id: s_id.to_string(),
            identity,
            length: len,
            mismatch: len.saturating_sub(matches),
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: h.e_value,
            bit_score,
            q_idx: h.q_idx,
            s_idx: h.s_idx,
            raw_score: h.raw_score,
        });
    }

    eprintln!("[4] Final hits after E-value filter (threshold={}): {}", evalue_threshold, final_hits.len());
    
    // Report top hit alignment length
    if !final_hits.is_empty() {
        let max_len = final_hits.iter().map(|h| h.length).max().unwrap_or(0);
        eprintln!("[5] Top alignment length: {} AA", max_len);
    }
    eprintln!("=== End Stage Counters ===");
    
    // NCBI-style output ordering: query (input order) → subject (best_evalue/score/oid) → HSP (score/coords)
    // Reference: BLAST_LinkHsps() + s_EvalueCompareHSPLists() + ScoreCompareHSPs()
    write_output_ncbi_order(final_hits, args.out.as_ref())?;
    Ok(())
}

// =============================================================================
// Tests for purge_hsps_with_common_endpoints (NCBI parity)
// =============================================================================
#[cfg(test)]
mod purge_tests {
    use super::purge_hsps_with_common_endpoints;
    use super::super::chaining::UngappedHit;

    /// Helper to create a minimal UngappedHit for testing purge logic.
    fn make_hit(
        ctx_idx: usize,
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        s_frame: i8,
        raw_score: i32,
    ) -> UngappedHit {
        UngappedHit {
            q_idx: 0,
            s_idx: 0,
            ctx_idx,
            s_f_idx: 0,
            q_frame: 1,
            s_frame,
            q_aa_start: q_start,
            q_aa_end: q_end,
            s_aa_start: s_start,
            s_aa_end: s_end,
            q_orig_len: 1000,
            s_orig_len: 1000,
            raw_score,
            e_value: 0.0,
            num_ident: 0, // Mock value for tests
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        }
    }

    /// NCBI parity test: equivalent to `testCheckHSPCommonEndpoints` from
    /// `blasthits_unit_test.cpp` lines 1163-1221.
    ///
    /// Original NCBI test data:
    /// ```c
    /// const int kHspCountStart = 9;
    /// const int kHspCountEnd = 3;
    /// const int kScores[kHspCountStart] =
    ///     { 1044, 995, 965, 219, 160, 125, 110, 107, 103 };
    /// const int kQueryOffsets[kHspCountStart] =
    ///     { 2, 2, 2, 236, 88, 259, 278, 259, 278 };
    /// const int kQueryEnds[kHspCountStart] =
    ///     { 322, 336, 300, 322, 182, 322, 341, 341, 341 };
    /// const int kSubjectOffsets[kHspCountStart] =
    ///     { 7, 7, 7, 194, 2, 194, 197, 194, 197 };
    /// const int kSubjectEnds[kHspCountStart] =
    ///     { 292, 293, 301, 292, 96, 292, 260, 260, 266 };
    /// const int kSurvivingIndices[kHspCountEnd] = { 4, 0, 6 };
    /// ```
    ///
    /// Expected surviving HSPs after purge: indices 4, 0, 6 (in output order).
    #[test]
    fn test_purge_ncbi_parity() {
        let scores = [1044, 995, 965, 219, 160, 125, 110, 107, 103];
        let q_offsets = [2, 2, 2, 236, 88, 259, 278, 259, 278];
        let q_ends = [322, 336, 300, 322, 182, 322, 341, 341, 341];
        let s_offsets = [7, 7, 7, 194, 2, 194, 197, 194, 197];
        let s_ends = [292, 293, 301, 292, 96, 292, 260, 260, 266];

        // Create HSPs with context=0, s_frame=0 (matching NCBI test which uses defaults)
        let mut hits: Vec<UngappedHit> = (0..9)
            .map(|i| {
                make_hit(
                    0,                   // ctx_idx (context)
                    q_offsets[i],        // q_aa_start
                    q_ends[i],           // q_aa_end
                    s_offsets[i],        // s_aa_start
                    s_ends[i],           // s_aa_end
                    0,                   // s_frame
                    scores[i] as i32,    // raw_score
                )
            })
            .collect();

        // Run purge
        purge_hsps_with_common_endpoints(&mut hits);

        // Expected: 3 surviving HSPs (indices 4, 0, 6 from original)
        assert_eq!(hits.len(), 3, "Expected 3 surviving HSPs, got {}", hits.len());

        // Verify the surviving HSPs match the expected ones.
        // The order after purge depends on the phase 2 sort order.
        // We check by finding each expected HSP in the result.
        let surviving_indices = [4, 0, 6];
        for &orig_idx in &surviving_indices {
            let found = hits.iter().any(|h| {
                h.q_aa_start == q_offsets[orig_idx]
                    && h.q_aa_end == q_ends[orig_idx]
                    && h.s_aa_start == s_offsets[orig_idx]
                    && h.s_aa_end == s_ends[orig_idx]
                    && h.raw_score == scores[orig_idx] as i32
            });
            assert!(
                found,
                "Expected HSP {} (score={}, q={}-{}, s={}-{}) not found in result",
                orig_idx,
                scores[orig_idx],
                q_offsets[orig_idx],
                q_ends[orig_idx],
                s_offsets[orig_idx],
                s_ends[orig_idx]
            );
        }
    }

    /// Test that subject.frame is used in duplicate detection.
    /// HSPs with same coordinates but different subject.frame should NOT be purged.
    #[test]
    fn test_purge_different_s_frame_not_duplicate() {
        // Two HSPs with identical coordinates but different s_frame
        let mut hits = vec![
            make_hit(0, 10, 50, 20, 60, 1, 100),  // s_frame = 1
            make_hit(0, 10, 50, 20, 60, 2, 90),   // s_frame = 2
        ];

        purge_hsps_with_common_endpoints(&mut hits);

        // Both should survive because s_frame differs
        assert_eq!(hits.len(), 2, "HSPs with different s_frame should both survive");
    }

    /// Test that context (ctx_idx) is used in duplicate detection.
    /// HSPs with same coordinates but different context should NOT be purged.
    #[test]
    fn test_purge_different_context_not_duplicate() {
        // Two HSPs with identical coordinates but different context
        let mut hits = vec![
            make_hit(0, 10, 50, 20, 60, 1, 100),  // ctx_idx = 0
            make_hit(1, 10, 50, 20, 60, 1, 90),   // ctx_idx = 1
        ];

        purge_hsps_with_common_endpoints(&mut hits);

        // Both should survive because ctx_idx differs
        assert_eq!(hits.len(), 2, "HSPs with different context should both survive");
    }

    /// Regression test documenting why mixed-subject purge is WRONG.
    ///
    /// NCBI's BlastHSPList is per-subject. If we purge a mixed collection,
    /// HSPs from different subjects could incorrectly be treated as duplicates
    /// if they happen to have the same ctx_idx + coordinates + s_frame.
    ///
    /// This test shows the CORRECT behavior: purging per-subject separately
    /// preserves both HSPs, while incorrectly purging together would remove one.
    #[test]
    fn test_mixed_subject_purge_is_wrong() {
        // Two HSPs from DIFFERENT subjects but with identical coordinates
        let hit1 = UngappedHit {
            q_idx: 0,
            s_idx: 0,  // subject 0
            ctx_idx: 0,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: 10,
            q_aa_end: 50,
            s_aa_start: 20,
            s_aa_end: 60,
            q_orig_len: 100,
            s_orig_len: 100,
            raw_score: 100,
            e_value: 0.0,
            num_ident: 0, // Mock value for tests
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        };
        let hit2 = UngappedHit {
            q_idx: 0,
            s_idx: 1,  // subject 1 (DIFFERENT)
            ctx_idx: 0,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: 10,
            q_aa_end: 50,
            s_aa_start: 20,
            s_aa_end: 60,
            q_orig_len: 100,
            s_orig_len: 100,
            raw_score: 90,
            e_value: 0.0,
            num_ident: 0, // Mock value for tests
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
        };

        // CORRECT: Purge per-subject separately (NCBI BlastHSPList semantics)
        let mut subject0_hits = vec![hit1.clone()];
        let mut subject1_hits = vec![hit2.clone()];
        purge_hsps_with_common_endpoints(&mut subject0_hits);
        purge_hsps_with_common_endpoints(&mut subject1_hits);
        let correct_count = subject0_hits.len() + subject1_hits.len();
        assert_eq!(correct_count, 2, "Per-subject purge should preserve both HSPs");

        // WRONG (what we must NOT do): Purge mixed collection
        // Note: The current implementation doesn't check s_idx, so mixing would
        // incorrectly treat these as duplicates. This documents the requirement
        // that purge must be called per-subject.
        let mut mixed_hits = vec![hit1, hit2];
        purge_hsps_with_common_endpoints(&mut mixed_hits);
        // If purge incorrectly treats different subjects as duplicates:
        // mixed_hits.len() would be 1 (wrong!)
        // This test documents the CORRECT behavior requirement.
        // The implementation relies on callers to split by subject first.
    }
}
