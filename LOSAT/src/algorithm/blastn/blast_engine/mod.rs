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
use super::filtering::purge_hsps_with_common_endpoints;

// Re-export for backward compatibility
pub use crate::algorithm::common::evalue::calculate_evalue_database_search as calculate_evalue;

/// Filter HSPs to remove redundant overlapping hits.
/// This function applies overlap filtering to match NCBI BLAST's behavior of outputting individual HSPs.
/// Chaining/clustering functionality has been removed to match NCBI BLAST blastn behavior.
pub fn filter_hsps(
    hits: Vec<Hit>,
    _sequences: &FxHashMap<(String, String), (Vec<u8>, Vec<u8>)>,
    _reward: i32,
    _penalty: i32,
    _gap_open: i32,
    _gap_extend: i32,
    _db_len_total: usize,
    _db_num_seqs: usize,
    _params: &KarlinParams,
    _use_dp: bool,
    verbose: bool,
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

    // NCBI BLAST: Purge HSPs with common start or end positions
    // NCBI reference: blast_hits.c:2455-2535 Blast_HSPListPurgeHSPsWithCommonEndpoints
    // Called from blast_engine.c:545 and blast_traceback.c:668
    let before_purge = result_hits.len();
    // TEMP: Skip endpoint purge to debug hit count
    let result_hits = if std::env::var("LOSAT_SKIP_PURGE").is_ok() {
        result_hits
    } else {
        purge_hsps_with_common_endpoints(result_hits)
    };
    if verbose && before_purge != result_hits.len() {
        eprintln!("[INFO]   Endpoint purge: {} -> {} hits", before_purge, result_hits.len());
    }

    // NCBI blastn does NOT apply culling by default (culling_limit = 0)
    // Culling is only applied when -culling_limit is explicitly set
    // Reference: blast_options.c:1496 - culling_limit defaults to 0
    // For NCBI parity, we skip the domination filtering entirely and return sorted hits
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
