//! Query context information management for blastn
//!
//! This module provides functionality to manage query contexts and determine
//! which context a given query offset belongs to, matching NCBI BLAST's
//! BlastQueryInfo and BSearchContextInfo implementation.

/// Context information for a single query context
///
/// NCBI reference: blast_query_info.h:60-73
/// ```c
/// typedef struct BlastContextInfo {
///     Int4 query_offset;      /**< Offset of this query, strand or frame */
///     Int4 query_length;      /**< Length of this query, strand or frame */
///     Int8 eff_searchsp;      /**< Effective search space for this context */
///     Int4 length_adjustment; /**< Length adjustment for boundary conditions */
///     Int4 query_index;       /**< Index of query (same for all frames) */
///     Int1 frame;             /**< Frame number (-1, -2, -3, 0, 1, 2, or 3) */
///     Boolean is_valid;        /**< Determine if this context is valid */
///     Int4 segment_flags;      /**< Flags describing segments for paired reads */
/// } BlastContextInfo;
/// ```
#[derive(Clone, Debug)]
pub struct ContextInfo {
    /// Offset of this query, strand or frame in the concatenated super-query
    pub query_offset: i32,
    /// Length of this query, strand or frame
    pub query_length: i32,
    /// Effective search space for this context
    pub eff_searchsp: i64,
    /// Length adjustment for boundary conditions
    pub length_adjustment: i32,
    /// Index of query (same for all frames/strands)
    pub query_index: i32,
    /// Frame number (-1, -2, -3, 0, 1, 2, or 3)
    /// For blastn: 0 = forward strand, 1 = reverse strand
    pub frame: i8,
    /// Determine if this context is valid
    pub is_valid: bool,
}

/// Query information structure
///
/// NCBI reference: blast_query_info.h:80-91
/// ```c
/// typedef struct BlastQueryInfo {
///     Int4 first_context;  /**< Index of the first element of the context array */
///     Int4 last_context;   /**< Index of the last element of the context array */
///     int num_queries;     /**< Number of query sequences */
///     BlastContextInfo * contexts; /**< Information per context */
///     Uint4 max_length;    /**< Length of the longest among the concatenated queries */
///     Uint4 min_length;    /**< Length of the shortest among the concatenated queries */
///     struct SPHIQueryInfo* pattern_info; /**< Counts of PHI BLAST pattern occurrences */
/// } BlastQueryInfo;
/// ```
#[derive(Clone, Debug)]
pub struct QueryInfo {
    /// Index of the first element of the context array
    pub first_context: i32,
    /// Index of the last element of the context array
    pub last_context: i32,
    /// Number of query sequences
    pub num_queries: usize,
    /// Information per context
    pub contexts: Vec<ContextInfo>,
    /// Length of the longest among the concatenated queries
    pub max_length: u32,
    /// Length of the shortest among the concatenated queries
    pub min_length: u32,
}

impl QueryInfo {
    /// Create a new QueryInfo for blastn queries
    ///
    /// For blastn, each query has 2 contexts:
    /// - Context 0: forward strand
    /// - Context 1: reverse strand
    ///
    /// NCBI reference: blast_setup.c (query info initialization)
    pub fn new_blastn(queries: &[Vec<u8>]) -> Self {
        let num_queries = queries.len();
        let mut contexts = Vec::new();
        let mut query_offset = 0i32;
        let mut max_length = 0u32;
        let mut min_length = u32::MAX;

        for (q_idx, query_seq) in queries.iter().enumerate() {
            let query_len = query_seq.len() as u32;
            max_length = max_length.max(query_len);
            min_length = min_length.min(query_len);

            // Context 0: forward strand
            contexts.push(ContextInfo {
                query_offset,
                query_length: query_len as i32,
                eff_searchsp: 0, // Will be set later
                length_adjustment: 0, // Will be set later
                query_index: q_idx as i32,
                frame: 0, // Forward strand
                is_valid: true,
            });

            // Advance offset by query_length + 1 (sentry byte between queries)
            // NCBI reference: blast_util.c:1098-1101
            // Frames/queries are concatenated with sentinel bytes between them
            query_offset += query_len as i32 + 1;

            // Context 1: reverse strand
            contexts.push(ContextInfo {
                query_offset,
                query_length: query_len as i32,
                eff_searchsp: 0, // Will be set later
                length_adjustment: 0, // Will be set later
                query_index: q_idx as i32,
                frame: 1, // Reverse strand
                is_valid: true,
            });

            // Advance offset by query_length + 1 (sentry byte between queries)
            query_offset += query_len as i32 + 1;
        }

        let last_context = if contexts.is_empty() {
            0
        } else {
            (contexts.len() - 1) as i32
        };

        Self {
            first_context: 0,
            last_context,
            num_queries,
            contexts,
            max_length: if max_length == 0 { 0 } else { max_length },
            min_length: if min_length == u32::MAX { 0 } else { min_length },
        }
    }
}

/// Search BlastContextInfo structures for the specified offset
///
/// NCBI reference: blast_query_info.c:219-243
/// ```c
/// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
/// {
///     Int4 m=0, b=0, e=0, size=0;
///     size = A->last_context+1;
///     if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
///         b = MIN(n / (A->max_length + 1), size - 1);
///         e = MIN(n / (A->min_length + 1) + 1, size);
///         ASSERT(e <= size);
///     }
///     else {
///         b = 0;
///         e = size;
///     }
///     while (b < e - 1) {
///         m = (b + e) / 2;
///         if (A->contexts[m].query_offset > n)
///             e = m;
///         else
///             b = m;
///     }
///     return b;
/// }
/// ```
///
/// # Arguments
/// * `n` - Query offset to search for
/// * `query_info` - Query information structure
///
/// # Returns
/// Context index for the given query offset
pub fn bsearch_context_info(n: i32, query_info: &QueryInfo) -> usize {
    // NCBI reference: blast_query_info.c:223
    // size = A->last_context+1;
    let size = (query_info.last_context + 1) as usize;

    // NCBI reference: blast_query_info.c:225-233
    // Optimize search range if possible
    let (mut b, mut e) = if query_info.min_length > 0
        && query_info.max_length > 0
        && query_info.first_context == 0
    {
        // NCBI: b = MIN(n / (A->max_length + 1), size - 1);
        // NCBI: e = MIN(n / (A->min_length + 1) + 1, size);
        let b_val = (n / (query_info.max_length as i32 + 1)).min(size as i32 - 1);
        let e_val = (n / (query_info.min_length as i32 + 1) + 1).min(size as i32);
        (b_val.max(0) as usize, e_val.max(0) as usize)
    } else {
        (0, size)
    };

    // NCBI reference: blast_query_info.c:235-241
    // Binary search to find the context
    while b < e.saturating_sub(1) {
        let m = (b + e) / 2;
        if query_info.contexts[m].query_offset > n {
            e = m;
        } else {
            b = m;
        }
    }

    b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bsearch_context_info_single_query() {
        let queries = vec![vec![0u8; 100]]; // Single query, length 100
        let query_info = QueryInfo::new_blastn(&queries);

        // Context 0: forward strand, offset 0
        assert_eq!(bsearch_context_info(0, &query_info), 0);
        assert_eq!(bsearch_context_info(50, &query_info), 0);

        // Context 1: reverse strand, offset 101 (100 + 1 sentinel)
        assert_eq!(bsearch_context_info(101, &query_info), 1);
        assert_eq!(bsearch_context_info(150, &query_info), 1);
    }

    #[test]
    fn test_bsearch_context_info_multiple_queries() {
        let queries = vec![
            vec![0u8; 100], // Query 0, length 100
            vec![0u8; 200], // Query 1, length 200
        ];
        let query_info = QueryInfo::new_blastn(&queries);

        // Query 0, context 0: offset 0-100
        assert_eq!(bsearch_context_info(0, &query_info), 0);
        assert_eq!(bsearch_context_info(50, &query_info), 0);

        // Query 0, context 1: offset 101-201
        assert_eq!(bsearch_context_info(101, &query_info), 1);
        assert_eq!(bsearch_context_info(150, &query_info), 1);

        // Query 1, context 0: offset 202-402
        assert_eq!(bsearch_context_info(202, &query_info), 2);
        assert_eq!(bsearch_context_info(300, &query_info), 2);

        // Query 1, context 1: offset 403-603
        assert_eq!(bsearch_context_info(403, &query_info), 3);
        assert_eq!(bsearch_context_info(500, &query_info), 3);
    }
}

