//! Constants for TBLASTX algorithm
//!
//! This module contains NCBI BLAST compatible parameters for translated BLAST searches.

/// NCBI BLAST compatible X-drop parameters for protein alignments
/// Using slightly higher X-drop (11 vs NCBI's 7) to allow longer extensions
/// before termination, which helps match NCBI BLAST's average hit length
pub const X_DROP_UNGAPPED: i32 = 11;

/// BLAST_GAP_X_DROPOFF_PROT for preliminary extension
pub const X_DROP_GAPPED_PRELIM: i32 = 15;

/// BLAST_GAP_X_DROPOFF_FINAL_PROT for final traceback
pub const X_DROP_GAPPED_FINAL: i32 = 25;

/// Maximum number of hits per k-mer before masking the bucket
pub const MAX_HITS_PER_KMER: usize = 200;

/// Stop codon encoding (25 = 'Z' + 1, used for non-standard amino acids)
pub const STOP_CODON: u8 = 25;

/// BLAST_WINDOW_SIZE_PROT from NCBI (was 16)
/// Window size for two-hit requirement in protein searches
pub const TWO_HIT_WINDOW: usize = 40;

/// Default gap penalties for protein alignments (BLOSUM62 defaults, NCBI compatible)
/// BLAST_GAP_OPEN_PROT
pub const GAP_OPEN: i32 = -11;

/// BLAST_GAP_EXTN_PROT
pub const GAP_EXTEND: i32 = -1;

/// Threshold for triggering gapped extension without two-hit requirement
pub const HIGH_SCORE_THRESHOLD: i32 = 60;

/// Minimum ungapped score threshold for keeping hits
/// NCBI BLAST outputs down to bit score ~22.1, which corresponds to ungapped score ~14
/// Setting to 22 to filter out short/low-scoring hits while allowing longer extensions
pub const MIN_UNGAPPED_SCORE: i32 = 22;

/// Parameters for HSP chaining (similar to BLASTN)
/// Maximum gap in amino acids (~1000bp / 3 for amino acids)
pub const MAX_GAP_AA: usize = 333;

/// Maximum diagonal drift in amino acids (~100bp / 3 for amino acids)
pub const MAX_DIAG_DRIFT_AA: isize = 33;

