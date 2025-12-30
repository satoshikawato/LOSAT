//! Constants for TBLASTX algorithm
//!
//! This module contains NCBI BLAST compatible parameters for translated BLAST searches.

/// X-drop for ungapped protein alignments (in BITS)
/// NCBI BLAST standard value = 7 bits
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h
/// #define BLAST_UNGAPPED_X_DROPOFF_PROT 7
/// 
/// IMPORTANT: NCBI converts this to raw score at runtime:
///   x_dropoff_raw = x_dropoff_bits * ln(2) / lambda
/// For BLOSUM62 ungapped (lambda ≈ 0.3176):
///   x_dropoff_raw = 7 * 0.693 / 0.3176 ≈ 15
/// 
/// We use the raw score value directly to match NCBI behavior.
pub const X_DROP_UNGAPPED_BITS: f64 = 7.0;

/// Pre-calculated X-drop raw score for ungapped protein alignments
/// This is the converted value: 7 bits * ln(2) / 0.3176 ≈ 15.27
/// Used directly in extension functions.
pub const X_DROP_UNGAPPED: i32 = 15;

/// BLAST_GAP_X_DROPOFF_PROT for preliminary extension
pub const X_DROP_GAPPED_PRELIM: i32 = 15;

/// BLAST_GAP_X_DROPOFF_FINAL_PROT for final traceback
pub const X_DROP_GAPPED_FINAL: i32 = 25;

/// Maximum number of hits per k-mer before masking the bucket
/// NCBI BLAST doesn't have an explicit limit like this.
/// Setting high to preserve more neighbor entries.
pub const MAX_HITS_PER_KMER: usize = 50000;

/// Stop codon encoding - index 24 in NCBI BLAST matrix order (ARNDCQEGHILKMFPSTWYVBJZX*)
/// BLOSUM62 gives stop codon scores of -4 vs other AAs, +1 vs itself
pub const STOP_CODON: u8 = 24;

/// BLAST_WINDOW_SIZE_PROT from NCBI
/// Window size for two-hit requirement in protein searches
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:57
/// #define BLAST_WINDOW_SIZE_PROT 40
pub const TWO_HIT_WINDOW: usize = 40;

/// Default gap penalties for protein alignments (BLOSUM62 defaults, NCBI compatible)
/// BLAST_GAP_OPEN_PROT
pub const GAP_OPEN: i32 = -11;

/// BLAST_GAP_EXTN_PROT
pub const GAP_EXTEND: i32 = -1;

/// Threshold for triggering gapped extension without two-hit requirement
pub const HIGH_SCORE_THRESHOLD: i32 = 60;

/// Minimum ungapped score threshold for keeping hits
/// NCBI BLAST outputs down to bit score ~22.1, which corresponds to raw score ~46 for BLOSUM62.
/// However, NCBI BLAST actually uses a cutoff based on E-value (CUTOFF_E_TBLASTX), not a fixed score.
/// Setting to 14 to capture more low-scoring hits that may be improved by sum-statistics linking.
/// This matches the observed minimum score range (12-21) that was being filtered out.
pub const MIN_UNGAPPED_SCORE: i32 = 14;

/// Parameters for HSP chaining (similar to BLASTN)
/// Maximum gap in amino acids (~1000bp / 3 for amino acids)
pub const MAX_GAP_AA: usize = 333;

/// Maximum diagonal drift in amino acids (~100bp / 3 for amino acids)
pub const MAX_DIAG_DRIFT_AA: isize = 33;

/// NCBI BLAST cutoff E-value for TBLASTX ungapped extensions
/// This is used to calculate the minimum score cutoff for extensions.
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_parameters.h:80
/// #define CUTOFF_E_TBLASTX 1e-300
pub const CUTOFF_E_TBLASTX: f64 = 1e-300;

/// NCBI BLAST gap trigger bit score for protein searches
/// HSPs with bit scores above this threshold are considered significant enough
/// to potentially trigger gapped alignment (though TBLASTX doesn't use gapped mode).
/// This is also used as a ceiling for cutoff_score_max.
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:137
/// #define BLAST_GAP_TRIGGER_PROT 22.0
pub const GAP_TRIGGER_BIT_SCORE: f64 = 22.0;

/// Sentinel byte for sequence boundaries (NCBI BLAST style)
/// 
/// NCBI BLAST uses NULLB (0) as a sentinel at the beginning and end of translated
/// sequences, and between frames. When extension reaches a sentinel, the matrix
/// returns a very negative score (defscore = -4 for unknown residues), causing
/// X-drop termination.
/// 
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:120
///   const Uint1 kProtSentinel = NULLB;
/// Reference: ncbi-blast/c++/src/util/tables/sm_blosum62.c:95
///   defscore = -4 (for characters not in the matrix)
/// 
/// In LOSAT, we use 255 as the sentinel value (outside the 0-24 amino acid range).
/// The extension functions check for this value and apply a large penalty.
pub const SENTINEL_BYTE: u8 = 255;

/// Penalty applied when extension hits a sentinel byte.
/// This should be large enough to trigger X-drop termination immediately.
/// NCBI uses -4 (defscore), but we use a larger value to ensure termination.
pub const SENTINEL_PENALTY: i32 = -100;

