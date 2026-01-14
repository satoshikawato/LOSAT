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
// REMOVED: MAX_HITS_PER_KMER - Over-represented k-mer filtering does not exist in NCBI BLAST
// NCBI reference: blast_lookup.c:BlastLookupAddWordHit (lines 33-77)
// NCBI BLAST adds all k-mers to lookup table regardless of frequency
// Database word count filtering (kDbFilter) exists but filters based on database counts, not query counts

/// Scan range for off-diagonal hit detection
/// NCBI reference: na_ungapped.c:658: Int4 Delta = MIN(word_params->options->scan_range, window_size - word_length);
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:63
/// ```c
/// #define BLAST_SCAN_RANGE_NUCL 0   /**< default scan range (blastn) */
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:165-168
/// ```c
/// SetWindowSize(BLAST_WINDOW_SIZE_NUCL);
/// SetOffDiagonalRange(BLAST_SCAN_RANGE_NUCL);
/// ```
pub const SCAN_RANGE_BLASTN: usize = 0;
pub const SCAN_RANGE_MEGABLAST: usize = 0;

// Minimum ungapped score thresholds for triggering gapped extension
// DEPRECATED: These fixed thresholds are replaced by dynamically calculated cutoff_score
// based on NCBI BLAST's implementation (blast_parameters.c:368-374).
// The dynamic calculation uses gap_trigger (bit score 27.0) and cutoff_score_max (from E-value).
// NCBI reference: na_ungapped.c:752: if (off_found || ungapped_data->score >= cutoffs->cutoff_score)
// NCBI reference: blast_parameters.c:343-344: gap_trigger calculation
// NCBI reference: blast_parameters.c:368-374: cutoff_score = MIN(gap_trigger * scale_factor, cutoff_score_max)
// These constants are kept for backward compatibility but are no longer used in the main code path.
pub const MIN_UNGAPPED_SCORE_MEGABLAST: i32 = 20; // DEPRECATED: Use compute_blastn_cutoff_score() instead
pub const MIN_UNGAPPED_SCORE_BLASTN: i32 = 50; // DEPRECATED: Use compute_blastn_cutoff_score() instead

/// Maximum word size for direct address table (4^14 = 268M entries = ~6GB, too large)
pub const MAX_DIRECT_LOOKUP_WORD_SIZE: usize = 13;

// NCBI BLAST threshold constants for adaptive lookup width selection
// Reference: blast_nalookup.c:130-184 (BlastChooseNaLookupTable)
// approx_table_entries thresholds for word_size=11:
// - < 12,000: use 8-bit lookup
// - < 180,000: use 10-bit lookup
// - >= 180,000: use 11-bit direct lookup
pub const LUT_WIDTH_11_THRESHOLD_8: usize = 12_000;
pub const LUT_WIDTH_11_THRESHOLD_10: usize = 180_000;

// NCBI BLAST constants for greedy alignment
// GREEDY_MAX_COST: The largest distance to be examined for an optimal alignment
pub const GREEDY_MAX_COST: usize = 1000;
// GREEDY_MAX_COST_FRACTION: sequence_length / this is a measure of how hard the algorithm will work
pub const GREEDY_MAX_COST_FRACTION: usize = 2;

/// Signal that a diagonal/offset is invalid
pub const INVALID_OFFSET: i32 = -2;
/// Large value for invalid diagonal bounds
pub const INVALID_DIAG: i32 = 100_000_000;

/// Minimum diagonal separation for HSP containment checking
/// NCBI reference: blast_nucl_options.cpp:239 (SetHitSavingOptionsDefaults)
/// NCBI reference: blast_nucl_options.cpp:259 (SetMBHitSavingOptionsDefaults)
/// Used in MB_HSP_CLOSE macro (blast_gapalign_priv.h:123-124)
/// Two HSPs are considered "close" if |diag1 - diag2| < min_diag_separation
pub const MIN_DIAG_SEPARATION_BLASTN: i32 = 50;
pub const MIN_DIAG_SEPARATION_MEGABLAST: i32 = 6;

