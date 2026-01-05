// NCBI BLAST compatible X-drop parameters
// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:122-148
pub const X_DROP_UNGAPPED: i32 = 20; // BLAST_UNGAPPED_X_DROPOFF_NUCL (blastn, megablast 共通)
pub const X_DROP_GAPPED_NUCL: i32 = 30; // BLAST_GAP_X_DROPOFF_NUCL (blastn, non-greedy)
pub const X_DROP_GAPPED_GREEDY: i32 = 25; // BLAST_GAP_X_DROPOFF_GREEDY (megablast, greedy)
pub const X_DROP_GAPPED_FINAL: i32 = 100; // BLAST_GAP_X_DROPOFF_FINAL_NUCL for final traceback (共通)

/// Two-hit window size for nucleotide searches
/// NCBI BLAST default: BLAST_WINDOW_SIZE_NUCL = 0 (one-hit mode)
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:58
/// When window_size = 0, all seeds trigger extension (one-hit mode)
/// When window_size > 0, two-hit requirement is enforced
/// NCBI reference: na_ungapped.c:656: Boolean two_hits = (window_size > 0);
pub const TWO_HIT_WINDOW: usize = 0; // NCBI BLAST default (one-hit mode)
pub const MAX_HITS_PER_KMER: usize = 200;

// Minimum ungapped score thresholds for triggering gapped extension
// Higher threshold for blastn (smaller word size = more seeds) to reduce extension count
// NCBI BLAST uses similar task-specific thresholds
pub const MIN_UNGAPPED_SCORE_MEGABLAST: i32 = 20; // Lower threshold for megablast (larger seeds)
pub const MIN_UNGAPPED_SCORE_BLASTN: i32 = 50; // Higher threshold for blastn (smaller seeds, more filtering needed)
// Increased from 40 to 50 to better handle self-comparison cases with many matches

/// Maximum word size for direct address table (4^14 = 268M entries = ~6GB, too large)
pub const MAX_DIRECT_LOOKUP_WORD_SIZE: usize = 13;

// NCBI BLAST constants for greedy alignment
// GREEDY_MAX_COST: The largest distance to be examined for an optimal alignment
pub const GREEDY_MAX_COST: usize = 1000;
// GREEDY_MAX_COST_FRACTION: sequence_length / this is a measure of how hard the algorithm will work
pub const GREEDY_MAX_COST_FRACTION: usize = 2;

/// Signal that a diagonal/offset is invalid
pub const INVALID_OFFSET: i32 = -2;
/// Large value for invalid diagonal bounds
pub const INVALID_DIAG: i32 = 100_000_000;

