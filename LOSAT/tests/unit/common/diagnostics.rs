//! Unit tests for common/diagnostics.rs

use LOSAT::algorithm::common::diagnostics::{
    diagnostics_enabled, BaseDiagnosticCounters, NucleotideDiagnosticCounters,
    ProteinDiagnosticCounters,
};
use std::sync::atomic::Ordering;
use std::env;

#[test]
fn test_diagnostics_enabled() {
    // Test default (disabled)
    env::remove_var("LOSAT_DIAGNOSTICS");
    assert!(!diagnostics_enabled());

    // Test enabled with "1"
    env::set_var("LOSAT_DIAGNOSTICS", "1");
    assert!(diagnostics_enabled());

    // Test enabled with "true"
    env::set_var("LOSAT_DIAGNOSTICS", "true");
    assert!(diagnostics_enabled());

    // Test enabled with "TRUE"
    env::set_var("LOSAT_DIAGNOSTICS", "TRUE");
    assert!(diagnostics_enabled());

    // Test disabled with "0"
    env::set_var("LOSAT_DIAGNOSTICS", "0");
    assert!(!diagnostics_enabled());

    // Test disabled with "false"
    env::set_var("LOSAT_DIAGNOSTICS", "false");
    assert!(!diagnostics_enabled());

    // Clean up
    env::remove_var("LOSAT_DIAGNOSTICS");
}

#[test]
fn test_base_diagnostic_counters_initialization() {
    let counters = BaseDiagnosticCounters::default();

    // All counters should start at 0
    assert_eq!(counters.kmer_matches.load(Ordering::Relaxed), 0);
    assert_eq!(counters.seeds_low_score.load(Ordering::Relaxed), 0);
    assert_eq!(counters.seeds_two_hit_failed.load(Ordering::Relaxed), 0);
    assert_eq!(counters.seeds_passed.load(Ordering::Relaxed), 0);
    assert_eq!(counters.ungapped_extensions.load(Ordering::Relaxed), 0);
    assert_eq!(counters.ungapped_low_score.load(Ordering::Relaxed), 0);
    assert_eq!(counters.hsps_before_chain.load(Ordering::Relaxed), 0);
    assert_eq!(counters.hsps_after_chain.load(Ordering::Relaxed), 0);
    assert_eq!(counters.hsps_after_overlap_filter.load(Ordering::Relaxed), 0);
    assert_eq!(counters.clusters_single.load(Ordering::Relaxed), 0);
    assert_eq!(counters.clusters_merged.load(Ordering::Relaxed), 0);
    assert_eq!(counters.hsps_in_merged_clusters.load(Ordering::Relaxed), 0);
    assert_eq!(counters.hsps_culled_dominated.load(Ordering::Relaxed), 0);
}

#[test]
fn test_base_diagnostic_counters_atomic_operations() {
    let counters = BaseDiagnosticCounters::default();

    // Test increment operations
    counters.kmer_matches.fetch_add(5, Ordering::Relaxed);
    assert_eq!(counters.kmer_matches.load(Ordering::Relaxed), 5);

    counters.seeds_passed.fetch_add(3, Ordering::Relaxed);
    assert_eq!(counters.seeds_passed.load(Ordering::Relaxed), 3);

    // Test multiple increments
    counters.ungapped_extensions.fetch_add(10, Ordering::Relaxed);
    counters.ungapped_extensions.fetch_add(5, Ordering::Relaxed);
    assert_eq!(counters.ungapped_extensions.load(Ordering::Relaxed), 15);
}

#[test]
fn test_protein_diagnostic_counters_initialization() {
    let counters = ProteinDiagnosticCounters::default();

    // Base counters should be initialized
    assert_eq!(counters.base.kmer_matches.load(Ordering::Relaxed), 0);

    // Protein-specific counters should be initialized
    assert_eq!(counters.ungapped_only_hits.load(Ordering::Relaxed), 0);
    assert_eq!(counters.ungapped_evalue_passed.load(Ordering::Relaxed), 0);
    assert_eq!(counters.ungapped_evalue_failed.load(Ordering::Relaxed), 0);
    assert_eq!(counters.gapped_extensions.load(Ordering::Relaxed), 0);
    assert_eq!(counters.gapped_evalue_passed.load(Ordering::Relaxed), 0);
    assert_eq!(counters.output_from_ungapped.load(Ordering::Relaxed), 0);
    assert_eq!(counters.output_from_gapped.load(Ordering::Relaxed), 0);
}

#[test]
fn test_protein_diagnostic_counters_atomic_operations() {
    let counters = ProteinDiagnosticCounters::default();

    // Test base counter operations
    counters.base.kmer_matches.fetch_add(7, Ordering::Relaxed);
    assert_eq!(counters.base.kmer_matches.load(Ordering::Relaxed), 7);

    // Test protein-specific counter operations
    counters.ungapped_only_hits.fetch_add(3, Ordering::Relaxed);
    assert_eq!(counters.ungapped_only_hits.load(Ordering::Relaxed), 3);

    counters.gapped_extensions.fetch_add(5, Ordering::Relaxed);
    assert_eq!(counters.gapped_extensions.load(Ordering::Relaxed), 5);
}

#[test]
fn test_nucleotide_diagnostic_counters_initialization() {
    let counters = NucleotideDiagnosticCounters::default();

    // Base counters should be initialized
    assert_eq!(counters.base.kmer_matches.load(Ordering::Relaxed), 0);
    assert_eq!(counters.base.seeds_passed.load(Ordering::Relaxed), 0);
}

#[test]
fn test_nucleotide_diagnostic_counters_atomic_operations() {
    let counters = NucleotideDiagnosticCounters::default();

    // Test base counter operations
    counters.base.kmer_matches.fetch_add(10, Ordering::Relaxed);
    assert_eq!(counters.base.kmer_matches.load(Ordering::Relaxed), 10);

    counters.base.seeds_passed.fetch_add(4, Ordering::Relaxed);
    assert_eq!(counters.base.seeds_passed.load(Ordering::Relaxed), 4);
}

#[test]
fn test_diagnostic_output_formatting() {
    let counters = ProteinDiagnosticCounters::default();

    // Set some values
    counters.base.kmer_matches.store(100, Ordering::Relaxed);
    counters.base.seeds_passed.store(50, Ordering::Relaxed);
    counters.ungapped_only_hits.store(20, Ordering::Relaxed);
    counters.gapped_extensions.store(30, Ordering::Relaxed);

    // The print_summary function should not panic
    // We can't easily test the output format without capturing stderr,
    // but we can at least verify it doesn't crash
    counters.print_summary(20);
}

#[test]
fn test_diagnostic_output_formatting_nucleotide() {
    let counters = NucleotideDiagnosticCounters::default();

    // Set some values
    counters.base.kmer_matches.store(200, Ordering::Relaxed);
    counters.base.seeds_passed.store(100, Ordering::Relaxed);

    // The print_summary function should not panic
    counters.print_summary(20);
}

