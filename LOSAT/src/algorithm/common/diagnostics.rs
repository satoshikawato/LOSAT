//! Diagnostic counters for BLASTN and TBLASTX pipelines
//!
//! This module provides diagnostic tracking to understand where hits are lost
//! in the alignment pipelines. Enabled via the LOSAT_DIAGNOSTICS environment variable.
//!
//! The module uses a trait-based design to support both nucleotide and protein
//! alignment diagnostics while sharing common infrastructure.

use std::sync::atomic::{AtomicUsize, AtomicI32, Ordering as AtomicOrdering};

/// Check if diagnostics are enabled via environment variable
pub fn diagnostics_enabled() -> bool {
    std::env::var("LOSAT_DIAGNOSTICS")
        .map(|v| v == "1" || v.to_lowercase() == "true")
        .unwrap_or(false)
}

/// Base diagnostic counters shared by all alignment types
#[derive(Default)]
pub struct BaseDiagnosticCounters {
    // Seed stage
    pub kmer_matches: AtomicUsize,
    pub seeds_low_score: AtomicUsize,
    pub seeds_two_hit_failed: AtomicUsize,  // Seeds that don't satisfy two-hit requirement (but still extended)
    pub seeds_passed: AtomicUsize,
    // Seed distribution tracking (for two-hit analysis)
    pub seeds_first_hit: AtomicUsize,      // First seed on a diagonal (no previous seed)
    pub seeds_second_hit_window: AtomicUsize, // Second seed within window (two-hit satisfied)
    pub seeds_second_hit_too_far: AtomicUsize, // Second seed beyond window (diff > window)
    pub seeds_second_hit_overlap: AtomicUsize, // Second seed overlapping (diff < wordsize)
    // Diagonal mask/suppression statistics
    pub mask_updates: AtomicUsize,         // Number of mask updates
    pub seeds_masked: AtomicUsize,         // Seeds suppressed by diagonal mask
    // Self-comparison detection
    pub self_comparisons: AtomicUsize,     // Number of query-subject pairs that are self-comparisons
    // Ungapped extension stage
    pub ungapped_extensions: AtomicUsize,
    pub ungapped_one_hit_extensions: AtomicUsize,  // Extensions with left-only (two-hit not satisfied) - NOTE: Always 0 in NCBI BLAST-compatible mode
    pub ungapped_two_hit_extensions: AtomicUsize,  // Extensions with left+right (two-hit satisfied)
    pub ungapped_low_score: AtomicUsize,
    // Extension length statistics (for analyzing extension behavior)
    pub extension_total_length: AtomicUsize,  // Sum of all extension lengths
    pub extension_max_length: AtomicUsize,    // Maximum extension length seen
    // Post-processing stage
    pub hsps_before_chain: AtomicUsize,
    pub hsps_after_chain: AtomicUsize,
    pub hsps_after_overlap_filter: AtomicUsize,
    // Cluster merging diagnostics
    pub clusters_single: AtomicUsize,    // Single-HSP clusters (output as-is)
    pub clusters_merged: AtomicUsize,    // Multi-HSP clusters that were merged
    pub hsps_in_merged_clusters: AtomicUsize, // Total HSPs that were part of merged clusters
    // HSP culling diagnostics
    pub hsps_culled_dominated: AtomicUsize, // HSPs culled by domination test
}

/// Diagnostic counters for protein/translated alignments (TBLASTX)
pub struct ProteinDiagnosticCounters {
    pub base: BaseDiagnosticCounters,
    // Gapped extension stage (protein-specific)
    pub ungapped_only_hits: AtomicUsize,
    pub ungapped_cutoff_failed: AtomicUsize,  // ungapped hits that failed cutoff_score filter
    pub ungapped_cutoff_failed_min_score: AtomicI32,  // minimum score of hits that failed cutoff_score
    pub ungapped_cutoff_failed_max_score: AtomicI32,  // maximum score of hits that failed cutoff_score
    pub ungapped_evalue_passed: AtomicUsize,  // ungapped-only hits that passed e-value filter
    pub ungapped_evalue_failed: AtomicUsize,  // ungapped-only hits that failed e-value filter
    pub gapped_extensions: AtomicUsize,
    pub gapped_evalue_passed: AtomicUsize,
    // Output source tracking
    pub output_from_ungapped: AtomicUsize,  // HSPs output from ungapped extension only
    pub output_from_gapped: AtomicUsize,    // HSPs output from gapped extension
}

impl Default for ProteinDiagnosticCounters {
    fn default() -> Self {
        Self {
            base: BaseDiagnosticCounters::default(),
            ungapped_only_hits: AtomicUsize::new(0),
            ungapped_cutoff_failed: AtomicUsize::new(0),
            ungapped_cutoff_failed_min_score: AtomicI32::new(i32::MAX),
            ungapped_cutoff_failed_max_score: AtomicI32::new(i32::MIN),
            ungapped_evalue_passed: AtomicUsize::new(0),
            ungapped_evalue_failed: AtomicUsize::new(0),
            gapped_extensions: AtomicUsize::new(0),
            gapped_evalue_passed: AtomicUsize::new(0),
            output_from_ungapped: AtomicUsize::new(0),
            output_from_gapped: AtomicUsize::new(0),
        }
    }
}

impl ProteinDiagnosticCounters {
    /// Print a summary of all diagnostic counters for TBLASTX
    pub fn print_summary(&self, min_ungapped_score: i32) {
        eprintln!("\n=== TBLASTX Pipeline Diagnostics ===");
        
        // Self-comparison detection
        let self_comps = self.base.self_comparisons.load(AtomicOrdering::Relaxed);
        if self_comps > 0 {
            eprintln!("Self-Comparison Detection:");
            eprintln!(
                "  Self-comparison pairs:      {} (query == subject)",
                self_comps
            );
            eprintln!("  WARNING: Self-comparisons may produce very long alignments");
        }
        
        eprintln!("Seed Stage:");
        eprintln!(
            "  K-mer matches found:        {}",
            self.base.kmer_matches.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Seeds filtered (low score): {}",
            self.base.seeds_low_score.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Seeds suppressed (mask):    {}",
            self.base.seeds_masked.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Seeds filtered (two-hit):   {}",
            self.base.seeds_two_hit_failed.load(AtomicOrdering::Relaxed)
        );
        eprintln!("  Seed distribution (two-hit analysis):");
        eprintln!(
            "    First seed on diagonal:    {}",
            self.base.seeds_first_hit.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "    Second seed (in window):   {}",
            self.base.seeds_second_hit_window.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "    Second seed (too far):     {}",
            self.base.seeds_second_hit_too_far.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "    Second seed (overlapping): {}",
            self.base.seeds_second_hit_overlap.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Seeds passed to extension:  {}",
            self.base.seeds_passed.load(AtomicOrdering::Relaxed)
        );
        
        // Diagonal mask statistics
        let mask_updates = self.base.mask_updates.load(AtomicOrdering::Relaxed);
        if mask_updates > 0 {
            eprintln!("Diagonal Mask Statistics:");
            eprintln!(
                "  Mask updates:               {}",
                mask_updates
            );
        }
        eprintln!("Ungapped Extension Stage:");
        eprintln!(
            "  Total extensions:           {}",
            self.base.ungapped_extensions.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  One-hit extensions (left-only): {}",
            self.base.ungapped_one_hit_extensions.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Two-hit extensions (left+right): {}",
            self.base.ungapped_two_hit_extensions.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Filtered (low score < {}):  {}",
            min_ungapped_score,
            self.base.ungapped_low_score.load(AtomicOrdering::Relaxed)
        );
        
        // Extension length statistics
        let total_len = self.base.extension_total_length.load(AtomicOrdering::Relaxed);
        let max_len = self.base.extension_max_length.load(AtomicOrdering::Relaxed);
        let num_ext = self.base.ungapped_extensions.load(AtomicOrdering::Relaxed);
        if num_ext > 0 {
            eprintln!("  Extension lengths:");
            eprintln!(
                "    Average length:           {:.1}",
                total_len as f64 / num_ext as f64
            );
            eprintln!(
                "    Maximum length:           {}",
                max_len
            );
            if max_len > 10000 {
                eprintln!("    WARNING: Very long extensions detected (may indicate self-comparison issue)");
            }
        }
        
        eprintln!("Ungapped-Only Hits (skipped gapped extension):");
        eprintln!(
            "  Total ungapped-only:        {}",
            self.ungapped_only_hits.load(AtomicOrdering::Relaxed)
        );
        let cutoff_failed_count = self.ungapped_cutoff_failed.load(AtomicOrdering::Relaxed);
        eprintln!(
            "  Filtered (cutoff_score):     {}",
            cutoff_failed_count
        );
        if cutoff_failed_count > 0 {
            let min_score = self.ungapped_cutoff_failed_min_score.load(AtomicOrdering::Relaxed);
            let max_score = self.ungapped_cutoff_failed_max_score.load(AtomicOrdering::Relaxed);
            // Display score range if at least one value was updated
            if min_score != i32::MAX || max_score != i32::MIN {
                if min_score == i32::MAX {
                    eprintln!(
                        "    Score range:                N/A - {}",
                        max_score
                    );
                } else if max_score == i32::MIN {
                    eprintln!(
                        "    Score range:                {} - N/A",
                        min_score
                    );
                } else {
                    eprintln!(
                        "    Score range:                {} - {}",
                        min_score, max_score
                    );
                }
            } else {
                eprintln!(
                    "    Score range:                (not updated)"
                );
            }
        }
        eprintln!(
            "  E-value passed:             {}",
            self.ungapped_evalue_passed.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  E-value failed:             {}",
            self.ungapped_evalue_failed.load(AtomicOrdering::Relaxed)
        );
        eprintln!("Gapped Extension Stage:");
        eprintln!(
            "  Gapped extensions run:      {}",
            self.gapped_extensions.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Gapped hits (e-value ok):   {}",
            self.gapped_evalue_passed.load(AtomicOrdering::Relaxed)
        );
        eprintln!("Post-Processing Stage:");
        eprintln!(
            "  HSPs before chaining:       {}",
            self.base.hsps_before_chain.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Clusters (single HSP):      {}",
            self.base.clusters_single.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Clusters (merged):          {}",
            self.base.clusters_merged.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  HSPs in merged clusters:    {}",
            self.base.hsps_in_merged_clusters.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  HSPs after chaining:        {}",
            self.base.hsps_after_chain.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  HSPs culled (dominated):    {}",
            self.base.hsps_culled_dominated.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Final hits (after filter):  {}",
            self.base.hsps_after_overlap_filter.load(AtomicOrdering::Relaxed)
        );
        eprintln!("Output Source:");
        eprintln!(
            "  From ungapped extension:    {}",
            self.output_from_ungapped.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  From gapped extension:      {}",
            self.output_from_gapped.load(AtomicOrdering::Relaxed)
        );
        eprintln!("====================================\n");
    }
}

/// Diagnostic counters for nucleotide alignments (BLASTN)
#[derive(Default)]
pub struct NucleotideDiagnosticCounters {
    pub base: BaseDiagnosticCounters,
    // BLASTN-specific counters can be added here
    // For now, BLASTN uses the base counters
}

impl NucleotideDiagnosticCounters {
    /// Print a summary of all diagnostic counters for BLASTN
    pub fn print_summary(&self, min_ungapped_score: i32) {
        eprintln!("\n=== BLASTN Pipeline Diagnostics ===");
        eprintln!("Seed Stage:");
        eprintln!(
            "  K-mer matches found:        {}",
            self.base.kmer_matches.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Seeds filtered (low score): {}",
            self.base.seeds_low_score.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Seeds filtered (two-hit):   {}",
            self.base.seeds_two_hit_failed.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Seeds passed to extension:  {}",
            self.base.seeds_passed.load(AtomicOrdering::Relaxed)
        );
        eprintln!("Ungapped Extension Stage:");
        eprintln!(
            "  Ungapped extensions run:    {}",
            self.base.ungapped_extensions.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Filtered (low score < {}):  {}",
            min_ungapped_score,
            self.base.ungapped_low_score.load(AtomicOrdering::Relaxed)
        );
        eprintln!("Post-Processing Stage:");
        eprintln!(
            "  HSPs before chaining:       {}",
            self.base.hsps_before_chain.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Clusters (single HSP):      {}",
            self.base.clusters_single.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Clusters (merged):          {}",
            self.base.clusters_merged.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  HSPs in merged clusters:    {}",
            self.base.hsps_in_merged_clusters.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  HSPs after chaining:        {}",
            self.base.hsps_after_chain.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  HSPs culled (dominated):    {}",
            self.base.hsps_culled_dominated.load(AtomicOrdering::Relaxed)
        );
        eprintln!(
            "  Final hits (after filter):  {}",
            self.base.hsps_after_overlap_filter.load(AtomicOrdering::Relaxed)
        );
        eprintln!("====================================\n");
    }
}

// For backward compatibility, provide a type alias for TBLASTX
pub type DiagnosticCounters = ProteinDiagnosticCounters;


