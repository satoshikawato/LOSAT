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
/// Inputs are **raw offsets into the sentinel-wrapped AA sequences**:
/// - `query_raw` and `subject_raw` include sentinel bytes at index 0 and last
/// - `q_off_raw` / `s_off_raw` point at the first aligned residue (not sentinel)
/// - `hsp_len` is the ungapped alignment length in residues
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

    // NCBI variables (translated => residue mask 0xff, i.e. no masking needed here).
    let mut score: i32 = 0;
    let mut sum: i32 = 0;

    // Pointers expressed as offsets within the HSP segment: [0..hsp_len]
    let mut best_start: usize = 0;
    let mut best_end: usize = 0;
    let mut current_start: usize = 0;

    for idx in 0..hsp_len {
        let q = query_raw[q_off_raw + idx];
        let s = subject_raw[s_off_raw + idx];
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

    let mut num_ident = 0usize;
    let align_length = hsp_len;

    // NCBI: for (i=0; i<align_length; i++)
    for i in 0..align_length {
        let q = query_nomask[q_off + i];
        let s = subject[s_off + i];
        
        // NCBI: if (*q == *s) num_ident++;
        if q == s {
            num_ident += 1;
        }
        // Note: NCBI also counts positives (matrix[*q][*s] > 0), but we only need num_ident for Blast_HSPTest
    }

    Some((num_ident, align_length))
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


