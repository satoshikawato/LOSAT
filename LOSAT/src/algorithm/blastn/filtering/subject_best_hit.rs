//! `-subject_besthit` for BLASTN's final HSP lists (`core::blast_hits::subject_besthit`, with
//! the pass over the other strand that only BLASTN runs).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2537-2606

use super::super::hsp::BlastnHsp;
use crate::core::blast_hits::{subject_besthit, BestHitExtent};

/// `Blast_HSPListSubjectBestHit` over a final HSP list of BLASTN, in its order: a query of
/// `query_len` bases, whose plus strand is context `2 * q_idx` and minus strand the next one,
/// with the HSPs' offsets in their own context (NCBI's `query.offset` and `query.end`, not the
/// report's coordinates).
pub fn subject_best_hit(hits: &mut Vec<BlastnHsp>, query_len: usize) {
    subject_besthit(hits, true, |h| BestHitExtent {
        context: i64::from(h.q_idx) * 2 + i64::from(h.query_frame < 0),
        query_frame: h.query_frame,
        offset: i64::from(h.internal_q_offset_0),
        end: i64::from(h.internal_q_end_0),
        query_length: query_len as i64,
    });
}

#[cfg(test)]
#[path = "subject_best_hit_tests.rs"]
mod tests;
