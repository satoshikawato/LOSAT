//! Constants for TBLASTX algorithm
//!
//! This module contains NCBI BLAST compatible parameters for translated BLAST searches.

/// NCBI BLAST compatible X-drop parameters for protein alignments
/// NCBI BLAST uses BLAST_UNGAPPED_X_DROPOFF_PROT = 7 for TBLASTX
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:122
/// #define BLAST_UNGAPPED_X_DROPOFF_PROT 7
pub const X_DROP_UNGAPPED: i32 = 7;

/// BLAST_GAP_X_DROPOFF_TBLASTX for preliminary extension
/// NCBI BLAST sets this to 0 to disable gapped extension for TBLASTX
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:134
/// #define BLAST_GAP_X_DROPOFF_TBLASTX 0
pub const X_DROP_GAPPED_PRELIM: i32 = 0;

/// BLAST_GAP_X_DROPOFF_FINAL_TBLASTX for final traceback
/// NCBI BLAST sets this to 0 to disable gapped extension for TBLASTX
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:148
/// #define BLAST_GAP_X_DROPOFF_FINAL_TBLASTX 0
pub const X_DROP_GAPPED_FINAL: i32 = 0;

/// Maximum number of hits per k-mer before masking the bucket
pub const MAX_HITS_PER_KMER: usize = 200;

/// Stop codon encoding in NCBI BLAST order (24 = '*' position in ARNDCQEGHILKMFPSTWYVBJZX*)
/// NCBI BLAST uses 25 amino acids: ARNDCQEGHILKMFPSTWYVBJZX* (indices 0-24)
/// Reference: ncbi-blast/c++/src/util/tables/sm_blosum62.c
pub const STOP_CODON: u8 = 24;

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
/// NCBI BLAST outputs down to bit score ~22.1
/// For BLOSUM62 with gap_open=11, gap_extend=1: lambda=0.267, K=0.041
/// bit_score = (lambda * ungapped_score - ln(K)) / ln(2)
/// bit_score 22.1 corresponds to ungapped_score ~45.4
/// Setting to 22 (bit_score ~13.1) to filter out very short/low-scoring hits
/// Note: This may be too high and could be filtering out valid low-identity hits
pub const MIN_UNGAPPED_SCORE: i32 = 22;

/// Parameters for HSP chaining (similar to BLASTN)
/// Maximum gap in amino acids (~1000bp / 3 for amino acids)
pub const MAX_GAP_AA: usize = 333;

/// Maximum diagonal drift in amino acids (~100bp / 3 for amino acids)
pub const MAX_DIAG_DRIFT_AA: isize = 33;

