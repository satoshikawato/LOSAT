//! NCBI-parity ungapped HSP reevaluation for nucleotide subject programs (incl. TBLASTX).
//!
//! NCBI runs `Blast_HSPListReevaluateUngapped` for **any ungapped search with a nucleotide
//! database**, including TBLASTX, before linking HSPs with sum statistics.
//! This can trim HSP boundaries and, in some cases, delete HSPs whose best rescored
//! subsegment falls below the cutoff.
//!
//! Reference:
//! - ncbi-blast/c++/src/algo/blast/core/blast_hits.c:675-733

use crate::utils::matrix::blosum62_score;
use std::sync::OnceLock;

// ---------------------------------------------------------------------------
// SIMD optimization: Reevaluation score table and gather helpers
// ---------------------------------------------------------------------------

// NCBI BLAST reference (c++/src/algo/blast/core/blast_hits.c:675-733):
//   for (index = 0; index < hsp_length; ++index) {
//       sum += matrix[*query & kResidueMask][*subject];
//       query++;
//       subject++;
//       if (sum < 0) {
//           sum = 0;
//           current_q_start = query;
//           current_s_start = subject;
//          if (score < cutoff_score) {
//             best_q_start = best_q_end = query;
//             best_s_start = best_s_end = subject;
//             score = 0;
//          }
//       } else if (sum > score) {
//          score = sum;
//          best_q_end = query;
//          best_s_end = subject;
//          best_q_start = current_q_start;
//          best_s_start = current_s_start;
//       }
//   }
//
// We batch the matrix lookup (blosum62_score) using AVX2 gather, but keep
// the exact per-residue control flow (sum updates, branches) to preserve
// 1-bit parity with NCBI.

static REEVALUATE_SCORE_TABLE_32: OnceLock<[i32; 32 * 32]> = OnceLock::new();

#[inline(always)]
fn reevaluate_score_table_32() -> &'static [i32; 32 * 32] {
    REEVALUATE_SCORE_TABLE_32.get_or_init(|| {
        let mut t = [0i32; 32 * 32];
        for a in 0..32u8 {
            for b in 0..32u8 {
                // Use blosum62_score directly to match scalar loop exactly
                // (includes sentinel/defscore handling)
                t[(a as usize) * 32 + (b as usize)] = blosum62_score(a, b);
            }
        }
        t
    })
}

/// Gather 8 BLOSUM62 scores using AVX2.
/// 
/// Loads 8 bytes from query and subject sequences, computes indices
/// (q<<5)|s, and gathers scores from the 32x32 table.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn gather_scores8_reevaluate_avx2(
    q_ptr: *const u8,
    s_ptr: *const u8,
    table: *const i32,
) -> [i32; 8] {
    use std::arch::x86_64::*;

    // Load 8 bytes from each sequence.
    let q8 = _mm_loadl_epi64(q_ptr as *const __m128i);
    let s8 = _mm_loadl_epi64(s_ptr as *const __m128i);

    // Expand to 8x u32 lanes.
    let q32 = _mm256_cvtepu8_epi32(q8);
    let s32 = _mm256_cvtepu8_epi32(s8);

    // Index into 32x32 score table: (q<<5)|s.
    let idx = _mm256_or_si256(_mm256_slli_epi32(q32, 5), s32);

    // Gather 8 i32 scores.
    let v = _mm256_i32gather_epi32(table, idx, 4);

    let mut out = [0i32; 8];
    _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, v);
    out
}

// ---------------------------------------------------------------------------
// Reevaluation functions
// ---------------------------------------------------------------------------

/// Reevaluate an ungapped HSP against the actual (possibly masked) translated sequences.
///
/// This is a direct port of NCBI:
/// - `Blast_HSPReevaluateWithAmbiguitiesUngapped` (blast_hits.c:675-733), for the
///   **translated** case (`translated == TRUE`), where `kResidueMask = 0xff`.
///
/// NCBI reference (verbatim, blast_hits.c:675-733):
/// ```c
/// Boolean
/// Blast_HSPReevaluateWithAmbiguitiesUngapped(BlastHSP* hsp, const Uint1* query_start,
///    const Uint1* subject_start, const BlastInitialWordParameters* word_params,
///    BlastScoreBlk* sbp, Boolean translated)
/// {
///    Int4 sum, score;
///    Int4** matrix;
///    const Uint1* query,* subject;
///    const Uint1* best_q_start,* best_s_start,* best_q_end,* best_s_end;
///    const Uint1* current_q_start, * current_s_start;
///    Int4 index;
///    const Uint1 kResidueMask = (translated ? 0xff : 0x0f);
///    Int4 hsp_length = hsp->query.end - hsp->query.offset;
///    Int4 cutoff_score = word_params->cutoffs[hsp->context].cutoff_score;
///
///    matrix = sbp->matrix->data;
///
///    query = query_start + hsp->query.offset;
///    subject = subject_start + hsp->subject.offset;
///    score = 0;
///    sum = 0;
///    best_q_start = best_q_end = current_q_start = query;
///    best_s_start = best_s_end = current_s_start = subject;
///
///    for (index = 0; index < hsp_length; ++index) {
///       sum += matrix[*query & kResidueMask][*subject];
///       query++;
///       subject++;
///       if (sum < 0) {
///           sum = 0;
///           current_q_start = query;
///           current_s_start = subject;
///          if (score < cutoff_score) {
///             best_q_start = best_q_end = query;
///             best_s_start = best_s_end = subject;
///             score = 0;
///          }
///       } else if (sum > score) {
///          score = sum;
///          best_q_end = query;
///          best_s_end = subject;
///          best_q_start = current_q_start;
///          best_s_start = current_s_start;
///       }
///    }
///
///    return s_UpdateReevaluatedHSPUngapped(hsp, cutoff_score, score,
///                                       query_start, subject_start, best_q_start,
///                                       best_q_end, best_s_start, best_s_end);
/// }
/// ```
///
/// Inputs are offsets into query/subject->sequence buffers (0-based, past NULLB):
/// - `query_raw` / `subject_raw` start at the first residue (sequence_start + 1)
/// - `q_off_raw` / `s_off_raw` are 0-based offsets from those starts
/// - `hsp_len` is the ungapped alignment length in residues
/// Reference: blast_util.c:112-116, blast_hits.c:675-733.
///
/// Returns updated (q_off_raw, s_off_raw, hsp_len, score) if kept, else `None`.
#[inline]
pub fn reevaluate_ungapped_hit_ncbi_translated(
    query_raw: &[u8],
    subject_raw: &[u8],
    q_off_raw: usize,
    s_off_raw: usize,
    hsp_len: usize,
    cutoff_score: i32,
) -> Option<(usize, usize, usize, i32)> {
    if hsp_len == 0 {
        return None;
    }
    // Bounds safety: NCBI assumes offsets/length are valid.
    if q_off_raw + hsp_len > query_raw.len() || s_off_raw + hsp_len > subject_raw.len() {
        return None;
    }

    // For short HSPs (< 16 residues), scalar is faster due to SIMD overhead
    // For longer HSPs, AVX2 gather can provide speedup
    #[cfg(target_arch = "x86_64")]
    {
        if hsp_len >= 16 && is_x86_feature_detected!("avx2") {
            unsafe {
                return reevaluate_ungapped_hit_ncbi_translated_avx2(
                    query_raw,
                    subject_raw,
                    q_off_raw,
                    s_off_raw,
                    hsp_len,
                    cutoff_score,
                );
            }
        }
    }

    reevaluate_ungapped_hit_ncbi_translated_scalar(
        query_raw,
        subject_raw,
        q_off_raw,
        s_off_raw,
        hsp_len,
        cutoff_score,
    )
}

/// Scalar implementation of reevaluation (fallback and reference).
#[inline(always)]
fn reevaluate_ungapped_hit_ncbi_translated_scalar(
    query_raw: &[u8],
    subject_raw: &[u8],
    q_off_raw: usize,
    s_off_raw: usize,
    hsp_len: usize,
    cutoff_score: i32,
) -> Option<(usize, usize, usize, i32)> {
    // NCBI variables (translated => residue mask 0xff, i.e. no masking needed here).
    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // Pointers expressed as offsets within the HSP segment: [0..hsp_len]
    let mut best_start: usize = 0;
    let mut best_end: usize = 0;
    let mut current_start: usize = 0;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:692-705
    // ```c
    // query = query_start + hsp->query.offset;
    // subject = subject_start + hsp->subject.offset;
    // for (index = 0; index < hsp_length; ++index) {
    //    sum += matrix[*query & kResidueMask][*subject];
    //    query++;
    //    subject++;
    //    if (sum < 0) {
    //       ...
    //    }
    // }
    // ```
    //
    // Safety: bounds were checked above, so q_ptr/s_ptr are valid for hsp_len bytes.
    let q_ptr = unsafe { query_raw.as_ptr().add(q_off_raw) };
    let s_ptr = unsafe { subject_raw.as_ptr().add(s_off_raw) };
    for idx in 0..hsp_len {
        let q = unsafe { *q_ptr.add(idx) };
        let s = unsafe { *s_ptr.add(idx) };
        // NCBI: matrix[*query & 0xff][*subject]
        sum += blosum62_score(q, s);

        if sum < 0 {
            // Start from new offset (NCBI lines 703-707)
            sum = 0;
            current_start = idx + 1;

            // If previous top score never reached cutoff, discard the front part completely.
            // (NCBI lines 708-715)
            if score < cutoff_score {
                best_start = idx + 1;
                best_end = idx + 1;
                score = 0;
            }
        } else if sum > score {
            // Remember this point as the best scoring end point (NCBI lines 716-725)
            score = sum;
            best_end = idx + 1; // query pointer is one past the current residue
            best_start = current_start;
        }
    }

    // NCBI s_UpdateReevaluatedHSPUngapped:
    // delete if score < cutoff_score
    if score < cutoff_score {
        return None;
    }

    let new_len = best_end.saturating_sub(best_start);
    if new_len == 0 {
        return None;
    }

    let new_q_off = q_off_raw + best_start;
    let new_s_off = s_off_raw + best_start;
    Some((new_q_off, new_s_off, new_len, score))
}

/// AVX2-optimized reevaluation that batches score lookup but preserves
/// per-residue control flow for NCBI parity.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn reevaluate_ungapped_hit_ncbi_translated_avx2(
    query_raw: &[u8],
    subject_raw: &[u8],
    q_off_raw: usize,
    s_off_raw: usize,
    hsp_len: usize,
    cutoff_score: i32,
) -> Option<(usize, usize, usize, i32)> {
    let table = reevaluate_score_table_32().as_ptr();

    // NCBI variables (translated => residue mask 0xff, i.e. no masking needed here).
    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // Pointers expressed as offsets within the HSP segment: [0..hsp_len]
    let mut best_start: usize = 0;
    let mut best_end: usize = 0;
    let mut current_start: usize = 0;

    let mut idx = 0usize;

    // Process blocks of 16 residues (2x8 gather calls)
    while idx + 16 <= hsp_len {
        let q_ptr = query_raw.as_ptr().add(q_off_raw + idx);
        let s_ptr = subject_raw.as_ptr().add(s_off_raw + idx);

        let scores0 = gather_scores8_reevaluate_avx2(q_ptr, s_ptr, table);
        let scores1 = gather_scores8_reevaluate_avx2(q_ptr.add(8), s_ptr.add(8), table);

        // Apply scores0[0..7] in order (strictly per-residue state machine)
        for k in 0..8usize {
            let pos = idx + k;
            sum += scores0[k];

            if sum < 0 {
                sum = 0;
                current_start = pos + 1;
                if score < cutoff_score {
                    best_start = pos + 1;
                    best_end = pos + 1;
                    score = 0;
                }
            } else if sum > score {
                score = sum;
                best_end = pos + 1;
                best_start = current_start;
            }
        }

        // Apply scores1[0..7] in order
        for k in 0..8usize {
            let pos = idx + 8 + k;
            sum += scores1[k];

            if sum < 0 {
                sum = 0;
                current_start = pos + 1;
                if score < cutoff_score {
                    best_start = pos + 1;
                    best_end = pos + 1;
                    score = 0;
                }
            } else if sum > score {
                score = sum;
                best_end = pos + 1;
                best_start = current_start;
            }
        }

        idx += 16;
    }

    // Process blocks of 8 residues
    while idx + 8 <= hsp_len {
        let q_ptr = query_raw.as_ptr().add(q_off_raw + idx);
        let s_ptr = subject_raw.as_ptr().add(s_off_raw + idx);
        let scores = gather_scores8_reevaluate_avx2(q_ptr, s_ptr, table);

        for k in 0..8usize {
            let pos = idx + k;
            sum += scores[k];

            if sum < 0 {
                sum = 0;
                current_start = pos + 1;
                if score < cutoff_score {
                    best_start = pos + 1;
                    best_end = pos + 1;
                    score = 0;
                }
            } else if sum > score {
                score = sum;
                best_end = pos + 1;
                best_start = current_start;
            }
        }

        idx += 8;
    }

    // Scalar tail
    while idx < hsp_len {
        let q = query_raw[q_off_raw + idx];
        let s = subject_raw[s_off_raw + idx];
        sum += blosum62_score(q, s);

        if sum < 0 {
            sum = 0;
            current_start = idx + 1;
            if score < cutoff_score {
                best_start = idx + 1;
                best_end = idx + 1;
                score = 0;
            }
        } else if sum > score {
            score = sum;
            best_end = idx + 1;
            best_start = current_start;
        }

        idx += 1;
    }

    // NCBI s_UpdateReevaluatedHSPUngapped:
    // delete if score < cutoff_score
    if score < cutoff_score {
        return None;
    }

    let new_len = best_end.saturating_sub(best_start);
    if new_len == 0 {
        return None;
    }

    let new_q_off = q_off_raw + best_start;
    let new_s_off = s_off_raw + best_start;
    Some((new_q_off, new_s_off, new_len, score))
}

/// Get number of identical and positive residues for an ungapped HSP.
///
/// This is a direct port of NCBI:
/// - `s_Blast_HSPGetNumIdentitiesAndPositives` (blast_hits.c:746-825), for the
///   **ungapped** case (no gap_info).
///
/// NCBI reference (verbatim, blast_hits.c:746-792):
/// ```c
/// s_Blast_HSPGetNumIdentitiesAndPositives(const Uint1* query, const Uint1* subject,
///                            			const BlastHSP* hsp, Int4* num_ident_ptr,
///                            			Int4* align_length_ptr, const BlastScoreBlk* sbp,
///                            			Int4* num_pos_ptr)
/// {
///    Int4 i, num_ident, align_length, q_off, s_off;
///    Uint1* q,* s;
///    Int4 q_length = hsp->query.end - hsp->query.offset;
///    Int4 s_length = hsp->subject.end - hsp->subject.offset;
///    Int4** matrix = NULL;
///    Int4 num_pos = 0;
///
///    q_off = hsp->query.offset;
///    s_off = hsp->subject.offset;
///
///    if ( !subject || !query || !hsp )
///       return -1;
///
///    q = (Uint1*) &query[q_off];
///    s = (Uint1*) &subject[s_off];
///
///    num_ident = 0;
///    align_length = 0;
///
///    if(NULL != sbp)
///    {
/// 	    if(sbp->protein_alphabet)
/// 		    matrix = sbp->matrix->data;
///    }
///
///    if (!hsp->gap_info) {
///       /* Ungapped case. Check that lengths are the same in query and subject,
///          then count number of matches. */
///       if (q_length != s_length)
///          return -1;
///       align_length = q_length;
///       for (i=0; i<align_length; i++) {
///          if (*q == *s)
///             num_ident++;
///          else if (NULL != matrix) {
///         	 if (matrix[*q][*s] > 0)
///         		 num_pos ++;
///             }
///          q++;
///          s++;
///       }
///   	}
/// }
/// ```
///
/// # Arguments
/// * `query_nomask` - Query sequence without masking (for identity calculation)
/// * `subject` - Subject sequence
/// * `q_off` - Query offset (0-based, into query_nomask)
/// * `s_off` - Subject offset (0-based, into subject)
/// * `hsp_len` - HSP length in residues
///
/// # Returns
/// `(num_ident, align_length)` if successful, `None` if lengths don't match
pub fn get_num_identities_and_positives_ungapped(
    query_nomask: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    hsp_len: usize,
) -> Option<(usize, usize)> {
    // NCBI: Check that lengths are the same in query and subject
    if q_off + hsp_len > query_nomask.len() || s_off + hsp_len > subject.len() {
        return None;
    }

    let align_length = hsp_len;

    // NCBI: for (i=0; i<align_length; i++)
    //   if (*q == *s) num_ident++;
    //
    // SIMD optimization: use AVX2/SSE2 to compare 32/16 bytes at a time
    let num_ident = unsafe {
        let q_ptr = query_nomask.as_ptr().add(q_off);
        let s_ptr = subject.as_ptr().add(s_off);
        
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                count_identities_avx2(q_ptr, s_ptr, align_length)
            } else if is_x86_feature_detected!("sse2") {
                count_identities_sse2(q_ptr, s_ptr, align_length)
            } else {
                count_identities_scalar(q_ptr, s_ptr, align_length)
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            count_identities_scalar(q_ptr, s_ptr, align_length)
        }
    };

    Some((num_ident, align_length))
}

/// Count identical residues using AVX2 (32 bytes at a time).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn count_identities_avx2(q_ptr: *const u8, s_ptr: *const u8, len: usize) -> usize {
    use std::arch::x86_64::*;

    let mut count = 0usize;
    let mut i = 0usize;

    // Process 32 bytes at a time
    while i + 32 <= len {
        let q_vec = _mm256_loadu_si256(q_ptr.add(i) as *const __m256i);
        let s_vec = _mm256_loadu_si256(s_ptr.add(i) as *const __m256i);
        
        // Compare bytes: 0xFF where equal, 0x00 where different
        let cmp = _mm256_cmpeq_epi8(q_vec, s_vec);
        
        // Get bitmask: each bit represents one byte comparison result
        let mask = _mm256_movemask_epi8(cmp) as u32;
        
        // Count set bits (number of equal bytes)
        count += mask.count_ones() as usize;
        
        i += 32;
    }

    // Process remaining bytes with SSE2 (16 bytes at a time)
    while i + 16 <= len {
        let q_vec = _mm_loadu_si128(q_ptr.add(i) as *const __m128i);
        let s_vec = _mm_loadu_si128(s_ptr.add(i) as *const __m128i);
        
        let cmp = _mm_cmpeq_epi8(q_vec, s_vec);
        let mask = _mm_movemask_epi8(cmp) as u16;
        
        count += mask.count_ones() as usize;
        i += 16;
    }

    // Scalar tail
    while i < len {
        if *q_ptr.add(i) == *s_ptr.add(i) {
            count += 1;
        }
        i += 1;
    }

    count
}

/// Count identical residues using SSE2 (16 bytes at a time).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn count_identities_sse2(q_ptr: *const u8, s_ptr: *const u8, len: usize) -> usize {
    use std::arch::x86_64::*;

    let mut count = 0usize;
    let mut i = 0usize;

    // Process 16 bytes at a time
    while i + 16 <= len {
        let q_vec = _mm_loadu_si128(q_ptr.add(i) as *const __m128i);
        let s_vec = _mm_loadu_si128(s_ptr.add(i) as *const __m128i);
        
        let cmp = _mm_cmpeq_epi8(q_vec, s_vec);
        let mask = _mm_movemask_epi8(cmp) as u16;
        
        count += mask.count_ones() as usize;
        i += 16;
    }

    // Scalar tail
    while i < len {
        if *q_ptr.add(i) == *s_ptr.add(i) {
            count += 1;
        }
        i += 1;
    }

    count
}

/// Count identical residues using scalar (fallback).
unsafe fn count_identities_scalar(q_ptr: *const u8, s_ptr: *const u8, len: usize) -> usize {
    let mut count = 0usize;
    for i in 0..len {
        if *q_ptr.add(i) == *s_ptr.add(i) {
            count += 1;
        }
    }
    count
}

/// Test if an HSP should be deleted based on percent identity and minimum hit length.
///
/// This is a direct port of NCBI:
/// - `s_HSPTest` (blast_hits.c:993-1001)
///
/// NCBI reference (verbatim, blast_hits.c:993-1001):
/// ```c
/// static Boolean s_HSPTest(const BlastHSP* hsp,
///                          const BlastHitSavingOptions* hit_options,
///                          Int4 align_length)
/// {
/// 	return ((hsp->num_ident * 100.0 <
/// 			align_length * hit_options->percent_identity) ||
/// 			align_length < hit_options->min_hit_length) ;
/// }
/// ```
///
/// # Arguments
/// * `num_ident` - Number of identical residues
/// * `align_length` - Alignment length
/// * `percent_identity` - Percent identity threshold (0.0 = disabled)
/// * `min_hit_length` - Minimum hit length (0 = disabled)
///
/// # Returns
/// `true` if HSP should be deleted, `false` otherwise
pub fn hsp_test(
    num_ident: usize,
    align_length: usize,
    percent_identity: f64,
    min_hit_length: usize,
) -> bool {
    // NCBI: (hsp->num_ident * 100.0 < align_length * hit_options->percent_identity)
    let identity_check = if percent_identity > 0.0 {
        (num_ident as f64 * 100.0) < (align_length as f64 * percent_identity)
    } else {
        false
    };

    // NCBI: align_length < hit_options->min_hit_length
    let length_check = if min_hit_length > 0 {
        align_length < min_hit_length
    } else {
        false
    };

    // NCBI: return (identity_check || length_check)
    identity_check || length_check
}

