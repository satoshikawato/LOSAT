// NCBI BLAST compatible X-drop parameters
// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:124-125
pub const X_DROP_UNGAPPED: i32 = 20; // BLAST_UNGAPPED_X_DROPOFF_NUCL
pub const X_DROP_GAPPED_FINAL: i32 = 100; // BLAST_GAP_X_DROPOFF_FINAL_NUCL for final traceback

/// Two-hit window size for BLASTN
/// Note: NCBI BLAST default is BLAST_WINDOW_SIZE_NUCL = 0 (two-hit requirement disabled)
/// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:58
/// However, LOSAT uses 64 for practical two-hit filtering to reduce extension count
/// This matches the behavior when window_size > 0 in na_ungapped.c:656
pub const TWO_HIT_WINDOW: usize = 64;
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


