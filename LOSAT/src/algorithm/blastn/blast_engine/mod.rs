//! BLASTN search coordination and execution
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c
//!            ncbi-blast/c++/src/algo/blast/core/blast_engine.c
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process:
//! - Query and subject sequence reading
//! - Lookup table construction
//! - Seed finding with two-hit extension
//! - HSP filtering and output generation
//!
//! # Module Structure
//!
//! - `run` - Main run() function for BLASTN search

mod run;

// Re-export the main run function
pub use run::run;

use rustc_hash::FxHashMap;
use crate::common::Hit;
use crate::stats::KarlinParams;
use super::filtering::{purge_hsps_with_common_endpoints, purge_hsps_with_common_endpoints_ex, reevaluate_hsp_with_ambiguities_gapped, subject_best_hit};

// Re-export for backward compatibility
pub use crate::algorithm::common::evalue::calculate_evalue_database_search as calculate_evalue;

/// Filter HSPs to remove redundant overlapping hits.
/// This function applies overlap filtering to match NCBI BLAST's behavior of outputting individual HSPs.
/// Chaining/clustering functionality has been removed to match NCBI BLAST blastn behavior.
///
/// NCBI reference: blast_traceback.c:637-669 (2-phase endpoint purge with trimming)
/// Phase 1: purge=FALSE - trim overlapping HSPs using edit scripts
/// Phase 2: Re-evaluate trimmed HSPs with Blast_HSPReevaluateWithAmbiguitiesGapped
/// Phase 3: purge=TRUE - delete remaining duplicates
pub fn filter_hsps(
    hits: Vec<Hit>,
    sequences: &FxHashMap<(String, String), (Vec<u8>, Vec<u8>)>,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    _db_len_total: usize,
    _db_num_seqs: usize,
    _params: &KarlinParams,
    _use_dp: bool,
    verbose: bool,
    query_lengths: &FxHashMap<String, usize>,  // Query ID -> length mapping
    cutoff_scores: &FxHashMap<String, i32>,    // Query ID -> cutoff score (NCBI: hit_params->cutoff_score_min)
) -> Vec<Hit> {
    if hits.is_empty() {
        return hits;
    }

    let total_start = std::time::Instant::now();

    // Always output individual HSPs (NCBI BLAST behavior)
    // NCBI BLAST's blastn does not use chaining/clustering - all HSPs are saved individually
    let mut result_hits: Vec<Hit> = hits;

    // NCBI blastn does NOT apply culling by default (culling_limit = 0)
    // Culling is only applied when -culling_limit is explicitly set
    // Reference: blast_options.c:1496 - culling_limit defaults to 0
    // For NCBI parity, we skip the domination filtering entirely
    let filter_start = std::time::Instant::now();
    let result_hits_count = result_hits.len();

    // Sort by raw_score (descending) to match NCBI BLAST's sorting order
    // NCBI reference: blast_hits.c:ScoreCompareHSPs
    // Order: score DESC → s_start ASC → s_end DESC → q_start ASC → q_end DESC
    result_hits.sort_by(|a, b| {
        b.raw_score
            .cmp(&a.raw_score)
            .then_with(|| a.s_start.cmp(&b.s_start))  // s_start ASC
            .then_with(|| b.s_end.cmp(&a.s_end))       // s_end DESC
            .then_with(|| a.q_start.cmp(&b.q_start))    // q_start ASC
            .then_with(|| b.q_end.cmp(&a.q_end))        // q_end DESC
    });

    // NOTE: NCBI does NOT have explicit dedup for identical coords+score.
    // NCBI only uses Blast_HSPListPurgeHSPsWithCommonEndpoints which handles
    // HSPs with common start OR end points (not both). Any exact duplicates
    // are handled by the endpoint purge below.
    // Reference: blast_hits.c:2455-2535

    // NCBI BLAST 2-phase endpoint purge with traceback trimming
    // NCBI reference: blast_traceback.c:637-669
    //
    // Int4 extra_start = Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
    // for (index=extra_start; index < hsp_list->hspcnt; index++) {
    //     hsp = hsp_array[index];
    //     if (!hsp) continue;
    //     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, ...);
    //     if (delete_hsp) hsp_array[index] = Blast_HSPFree(hsp_array[index]);
    // }
    // Blast_HSPListPurgeNullHSPs(hsp_list);
    // if(program_number == eBlastTypeBlastn) {
    //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // }

    let before_purge = result_hits.len();

    // Phase 1: purge=FALSE - trim overlapping HSPs using edit scripts
    // NCBI reference: blast_traceback.c:638
    let (mut result_hits, extra_start) = purge_hsps_with_common_endpoints_ex(result_hits, false);
    if verbose {
        eprintln!("[INFO]   Endpoint purge phase 1 (trim): {} -> {} hits, extra_start={}",
                  before_purge, result_hits.len(), extra_start);
    }

    // Phase 2: Re-evaluate trimmed HSPs
    // NCBI reference: blast_traceback.c:647-665
    // for (index=extra_start; index < hsp_list->hspcnt; index++) {
    //     ...
    //     delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query, query_length,
    //                      subject, subject_length, hit_params, score_params, sbp);
    //     ...
    // }
    let mut deleted_count = 0;
    for i in extra_start..result_hits.len() {
        let hit = &mut result_hits[i];

        // Get cutoff score for this query
        // NCBI reference: blast_traceback.c:654 - uses hit_params->cutoff_score_min
        // The cutoff_score is pre-computed in run() using compute_blastn_cutoff_score()
        let cutoff_score = cutoff_scores.get(&hit.query_id)
            .copied()
            .unwrap_or(22);

        // Get sequences for re-evaluation
        let key = (hit.query_id.clone(), hit.subject_id.clone());
        if let Some((q_seq, s_seq)) = sequences.get(&key) {
            // NCBI reference: blast_hits.c:479-647 Blast_HSPReevaluateWithAmbiguitiesGapped
            let delete = reevaluate_hsp_with_ambiguities_gapped(
                hit,
                q_seq,
                s_seq,
                reward,
                penalty,
                gap_open,
                gap_extend,
                cutoff_score,
            );
            if delete {
                hit.raw_score = i32::MIN; // Mark for deletion
                deleted_count += 1;
            }
        }
    }

    // Remove marked HSPs (equivalent to Blast_HSPListPurgeNullHSPs)
    // NCBI reference: blast_traceback.c:666
    result_hits.retain(|h| h.raw_score != i32::MIN);
    if verbose && deleted_count > 0 {
        eprintln!("[INFO]   Re-evaluation deleted {} trimmed HSPs", deleted_count);
    }

    // Phase 3: purge=TRUE for BLASTN
    // NCBI reference: blast_traceback.c:668
    // if(program_number == eBlastTypeBlastn) {
    //     Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    // }
    let before_phase3 = result_hits.len();
    let (mut result_hits, _) = purge_hsps_with_common_endpoints_ex(result_hits, true);
    if verbose && before_phase3 != result_hits.len() {
        eprintln!("[INFO]   Endpoint purge phase 3 (delete): {} -> {} hits", before_phase3, result_hits.len());
    }

    // NOTE: NCBI's containment filtering (BlastIntervalTreeContainsHSP) is applied
    // DURING extension, not as post-processing. The interval tree is built incrementally
    // as HSPs are found, and each new HSP is checked against existing HSPs.
    // This is different from post-processing containment which would be too aggressive
    // (e.g., removing all hits on self-comparison except the full-length diagonal).
    //
    // TODO: Implement interval tree containment checking during gapped extension
    // NCBI reference: blast_gapalign.c:3918, blast_itree.c:810-847

    // NCBI BLAST: Subject Best Hit filtering (OPTIONAL - not applied by default)
    // Reference: blast_hits.c:2537-2606 Blast_HSPListSubjectBestHit
    // Called from blast_engine.c:586-591 ONLY when subject_besthit_opts is set
    // (i.e., when user specifies -subject_besthit on command line)
    //
    // NOTE: Subject best hit filtering is disabled by default to match NCBI behavior.
    // The subject_best_hit() function is available but not called here.
    // To enable, add a command-line option and call subject_best_hit() when enabled.
    let _ = subject_best_hit; // Silence unused import warning
    let final_hits: Vec<Hit> = result_hits;

    let filter_time = filter_start.elapsed();

    if verbose {
        eprintln!("[INFO] filter_hsps summary:");
        eprintln!(
            "[INFO]   Total time: {:.3}s",
            total_start.elapsed().as_secs_f64()
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
