//! The results of a search by query and by subject: views of the flat HSP list.
//!
//! The final hit list holds a query's HSPs together and, in them, a subject's HSPs together, as
//! NCBI's result set holds a `CSearchResults` per query and a Seq-align set per subject, so these
//! views are ranges of it and copy nothing.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_results.hpp:361-386
//! ```c++
//!     /// @param i query sequence index if result type is eDatabaseSearch,
//!     /// otherwise it's the query-subject index [in]
//!     CSearchResults & operator[](size_type i) {
//!     ...
//!     /// @note it only works for results of type eSequenceComparison
//!     CSearchResults & GetResults(size_type qi, size_type si);
//! ```
#![warn(missing_docs, unused)]

use super::run::{HspRecord, Outcome, QueryStats};

/// The results of one query of a search.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct QueryHits<'a> {
    /// The index of the query record in the inputs.
    pub q_idx: usize,
    /// The query's statistics (`Outcome::query_stats`).
    pub stats: Option<&'a QueryStats>,
    /// The query's HSPs, in the order of the hit list. Empty for a query without hits.
    pub hsps: &'a [HspRecord],
}

/// The HSPs of one query against one subject.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct SubjectHits<'a> {
    /// The index of the query record in the inputs.
    pub q_idx: usize,
    /// The index of the subject record in the inputs.
    pub s_idx: usize,
    /// The HSPs of the pair, in the order of the hit list.
    pub hsps: &'a [HspRecord],
}

/// The HSPs of query `q_idx` in a hit list, which is in the order of the queries.
fn hsps_of_query(hsps: &[HspRecord], q_idx: usize) -> &[HspRecord] {
    let start = hsps.partition_point(|hsp| hsp.q_idx < q_idx);
    let end = start + hsps[start..].partition_point(|hsp| hsp.q_idx == q_idx);
    &hsps[start..end]
}

/// The results of query `q_idx` in a hit list.
pub(super) fn query_hits<'a>(
    hsps: &'a [HspRecord],
    q_idx: usize,
    stats: Option<&'a QueryStats>,
) -> QueryHits<'a> {
    QueryHits {
        q_idx,
        stats,
        hsps: hsps_of_query(hsps, q_idx),
    }
}

impl<'a> QueryHits<'a> {
    /// The subjects with HSPs for this query, in the order of the hit list: the best first.
    pub fn subjects(&self) -> impl Iterator<Item = SubjectHits<'a>> + 'a {
        let q_idx = self.q_idx;
        self.hsps
            .chunk_by(|left, right| left.s_idx == right.s_idx)
            .map(move |hsps| SubjectHits {
                q_idx,
                s_idx: hsps[0].s_idx,
                hsps,
            })
    }

    /// The HSPs of this query against subject `s_idx`, if it has any.
    pub fn subject(&self, s_idx: usize) -> Option<SubjectHits<'a>> {
        self.subjects().find(|subject| subject.s_idx == s_idx)
    }
}

impl Outcome {
    /// The results of each query record of the inputs, in their order, also of those without
    /// hits.
    pub fn queries(&self) -> impl Iterator<Item = QueryHits<'_>> + '_ {
        self.query_stats
            .iter()
            .enumerate()
            .map(|(q_idx, stats)| query_hits(&self.hsps, q_idx, stats.as_ref()))
    }

    /// The results of query record `q_idx` of the inputs.
    pub fn query(&self, q_idx: usize) -> Option<QueryHits<'_>> {
        let stats = self.query_stats.get(q_idx)?;
        Some(query_hits(&self.hsps, q_idx, stats.as_ref()))
    }
}
