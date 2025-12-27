//! Diagnostic counters for TBLASTX pipeline
//!
//! This module provides diagnostic tracking to understand where hits are lost
//! in the TBLASTX pipeline. Enabled via the LOSAT_DIAGNOSTICS environment variable.

use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use super::constants::MIN_UNGAPPED_SCORE;

/// Diagnostic counters for tracking where hits are lost in the pipeline
#[derive(Default)]
pub struct DiagnosticCounters {
    // Seed stage
    pub kmer_matches: AtomicUsize,
    pub seeds_low_score: AtomicUsize,
    pub seeds_two_hit_failed: AtomicUsize,
    pub seeds_passed: AtomicUsize,
    // Ungapped extension stage
    pub ungapped_extensions: AtomicUsize,
    pub ungapped_low_score: AtomicUsize,
    // Gapped extension stage
    pub ungapped_only_hits: AtomicUsize,
    pub ungapped_evalue_passed: AtomicUsize,  // ungapped-only hits that passed e-value filter
    pub ungapped_evalue_failed: AtomicUsize,  // ungapped-only hits that failed e-value filter
    pub gapped_extensions: AtomicUsize,
    pub gapped_evalue_passed: AtomicUsize,
    // Post-processing stage (set in chain_and_filter_hsps_protein)
    pub hsps_before_chain: AtomicUsize,
    pub hsps_after_chain: AtomicUsize,
    pub hsps_after_overlap_filter: AtomicUsize,
    // Cluster merging diagnostics
    pub clusters_single: AtomicUsize,    // Single-HSP clusters (output as-is)
    pub clusters_merged: AtomicUsize,    // Multi-HSP clusters that were merged
    pub hsps_in_merged_clusters: AtomicUsize, // Total HSPs that were part of merged clusters
    // HSP culling diagnostics
    pub output_from_ungapped: AtomicUsize,  // HSPs output from ungapped extension only
    pub output_from_gapped: AtomicUsize,    // HSPs output from gapped extension
    pub hsps_culled_dominated: AtomicUsize, // HSPs culled by domination test
}

impl DiagnosticCounters {
    /// Print a summary of all diagnostic counters
    pub fn print_summary(&self) {
        eprintln!("\n=== TBLASTX Pipeline Diagnostics ===");
        eprintln!("Seed Stage:");
        eprintln!("  K-mer matches found:        {}", self.kmer_matches.load(AtomicOrdering::Relaxed));
        eprintln!("  Seeds filtered (low score): {}", self.seeds_low_score.load(AtomicOrdering::Relaxed));
        eprintln!("  Seeds filtered (two-hit):   {}", self.seeds_two_hit_failed.load(AtomicOrdering::Relaxed));
        eprintln!("  Seeds passed to extension:  {}", self.seeds_passed.load(AtomicOrdering::Relaxed));
        eprintln!("Ungapped Extension Stage:");
        eprintln!("  Ungapped extensions run:    {}", self.ungapped_extensions.load(AtomicOrdering::Relaxed));
        eprintln!("  Filtered (low score < {}):  {}", MIN_UNGAPPED_SCORE, self.ungapped_low_score.load(AtomicOrdering::Relaxed));
        eprintln!("Ungapped-Only Hits (skipped gapped extension):");
        eprintln!("  Total ungapped-only:        {}", self.ungapped_only_hits.load(AtomicOrdering::Relaxed));
        eprintln!("  E-value passed:             {}", self.ungapped_evalue_passed.load(AtomicOrdering::Relaxed));
        eprintln!("  E-value failed:             {}", self.ungapped_evalue_failed.load(AtomicOrdering::Relaxed));
        eprintln!("Gapped Extension Stage:");
        eprintln!("  Gapped extensions run:      {}", self.gapped_extensions.load(AtomicOrdering::Relaxed));
        eprintln!("  Gapped hits (e-value ok):   {}", self.gapped_evalue_passed.load(AtomicOrdering::Relaxed));
        eprintln!("Post-Processing Stage:");
        eprintln!("  HSPs before chaining:       {}", self.hsps_before_chain.load(AtomicOrdering::Relaxed));
        eprintln!("  Clusters (single HSP):      {}", self.clusters_single.load(AtomicOrdering::Relaxed));
        eprintln!("  Clusters (merged):          {}", self.clusters_merged.load(AtomicOrdering::Relaxed));
        eprintln!("  HSPs in merged clusters:    {}", self.hsps_in_merged_clusters.load(AtomicOrdering::Relaxed));
        eprintln!("  HSPs after chaining:        {}", self.hsps_after_chain.load(AtomicOrdering::Relaxed));
        eprintln!("  HSPs culled (dominated):    {}", self.hsps_culled_dominated.load(AtomicOrdering::Relaxed));
        eprintln!("  Final hits (after filter):  {}", self.hsps_after_overlap_filter.load(AtomicOrdering::Relaxed));
        eprintln!("Output Source:");
        eprintln!("  From ungapped extension:    {}", self.output_from_ungapped.load(AtomicOrdering::Relaxed));
        eprintln!("  From gapped extension:      {}", self.output_from_gapped.load(AtomicOrdering::Relaxed));
        eprintln!("====================================\n");
    }
}

/// Check if diagnostics are enabled via environment variable
pub fn diagnostics_enabled() -> bool {
    std::env::var("LOSAT_DIAGNOSTICS").map(|v| v == "1" || v.to_lowercase() == "true").unwrap_or(false)
}

