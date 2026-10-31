//! blastp's preliminary hit lists hold a 64-byte element instead of a `BlastpHsp` (168 bytes)
//! until the redo reads them.
//!
//! The preliminary stage keeps, for every query of a search, the HSP lists of up to the
//! preliminary hit list size of subjects. Each preliminary HSP is a `BlastpHsp`, most of whose
//! fields (identities, bit score, edit script, ...) the preliminary stage never sets and the redo
//! never reads: the redo turns every HSP of a list back into the `BlastpPreliminaryHsp` it was
//! made from (`preliminary_hit_from_local_match_hsp`). With many matches in a block of queries
//! these lists are the largest part of the peak at the end of the subject scan.
//!
//! The subject scan keeps an `XPrelimHsp` (the `BlastpPreliminaryHsp` and its E-value) in the
//! `BlastpHspList` / `BlastpHitList` helpers, made generic over the element (`HitListHsp`): the
//! same E-values, scores and offsets in the same comparisons, and every sort of HSPs is stable
//! (`sort_by`), so the lists end in the same order. Before the redo a list is made into the
//! `BlastpHitList` of full HSPs the redo reads (`x_expand_hit_list`, with the HSPs built by
//! `blastp_hsp_from_preliminary`): the one-thread redo does it query by query and drops the
//! expanded list after the query's redo (nothing reads it afterwards, `x_prelim_free.rs`); a
//! search with more than one thread expands every list before the redo (its parallel paths read
//! the lists of all queries), but where the redo runs in windows (LOSAT_X_BPREDOPAR),
//! the windows expand the lists as they reach them (`x_redo_par.rs`). The
//! lists of query chunks (`-query` longer than the split size) are expanded where a chunk's search
//! ends, so the chunk merge is unchanged.
//!
//! With `LOSAT_TIMING=1` a search prints `[TIMING] x_bppcompact: expanded_lists=<n>` (the hit lists
//! expanded during the search, those of its query chunks included).

use std::sync::atomic::{AtomicU64, Ordering};

use super::gapalign::BlastpPreliminaryHsp;
use super::hsp::{BlastpHitList, BlastpHspList, HitListHsp};

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:151-178
// ```c
// Blast_HSPInit(Int4 query_start, Int4 query_end, Int4 subject_start,
//               Int4 subject_end, Int4 query_gapped_start,
//               Int4 subject_gapped_start, Int4 query_context,
//               Int2 query_frame, Int2 subject_frame, Int4 score, ...)
// ```
// NCBI keeps one BlastHSP; this is the part of it the preliminary stage sets, with the E-value of
// Blast_HSPListGetEvalues (blast_engine.c:1523-1530).
/// A preliminary HSP as the subject scan keeps it.
#[derive(Debug, Clone, Copy)]
pub(crate) struct XPrelimHsp {
    pub hsp: BlastpPreliminaryHsp,
    pub e_value: f64,
}

// No NCBI counterpart: the values `BlastpHsp` holds for the same HSP
// (`blastp_hsp_from_preliminary`: starts + 1).
impl HitListHsp for XPrelimHsp {
    #[inline]
    fn e_value(&self) -> f64 {
        self.e_value
    }
    #[inline]
    fn raw_score(&self) -> i32 {
        self.hsp.raw_score
    }
    #[inline]
    fn offsets(&self) -> (usize, usize, usize, usize) {
        (
            (self.hsp.query_start + 1) as usize,
            self.hsp.query_end as usize,
            (self.hsp.subject_start + 1) as usize,
            self.hsp.subject_end as usize,
        )
    }
    #[inline]
    fn query_length(&self) -> usize {
        self.hsp.query_length
    }
}

// No NCBI counterpart: rebuilds the hit list of full HSPs from the compact one (same counts, same
// list order, the HSPs `blastp_hsp_from_preliminary` makes).
/// The `BlastpHitList` of full HSPs the redo reads for this list.
pub(crate) fn x_expand_hit_list(list: BlastpHitList<XPrelimHsp>) -> BlastpHitList {
    X_EXPANDED_LISTS.fetch_add(1, Ordering::Relaxed);
    let BlastpHitList {
        hsplist_count,
        hsplist_max,
        worst_evalue,
        low_score,
        heapified,
        hsplist_array,
        hsplist_current,
        num_hits,
    } = list;
    let mut expanded = Vec::with_capacity(hsplist_array.capacity());
    expanded.extend(hsplist_array.into_iter().map(|hsp_list| {
        BlastpHspList {
            oid: hsp_list.oid,
            query_index: hsp_list.query_index,
            hsps: hsp_list
                .hsps
                .iter()
                .map(|hsp| super::blast_engine::blastp_hsp_from_preliminary(&hsp.hsp, hsp.e_value))
                .collect(),
            best_evalue: hsp_list.best_evalue,
        }
    }));
    BlastpHitList {
        hsplist_count,
        hsplist_max,
        worst_evalue,
        low_score,
        heapified,
        hsplist_array: expanded,
        hsplist_current,
        num_hits,
    }
}

// No NCBI counterpart: hands the preliminary hit lists on as lists of full HSPs.
/// Every compact list becomes the list of full HSPs of its query (`full` holds the lists of the
/// queries whose chunks replaced them).
pub(crate) fn x_finish_compact_lists(
    full: &mut [Option<BlastpHitList>],
    compact: &mut [Option<BlastpHitList<XPrelimHsp>>],
) {
    for (full, compact) in full.iter_mut().zip(compact.iter_mut()) {
        if let Some(list) = compact.take() {
            *full = Some(x_expand_hit_list(list));
        }
    }
}

static X_EXPANDED_LISTS: AtomicU64 = AtomicU64::new(0);

/// The hit lists expanded so far (for the timing line).
pub(crate) fn x_expanded_lists() -> u64 {
    X_EXPANDED_LISTS.load(Ordering::Relaxed)
}

#[cfg(test)]
mod tests {
    use super::super::hsp::BlastpHsp;
    use super::*;

    fn same_hsp(a: &BlastpHsp, b: &BlastpHsp) -> bool {
        a.identity.to_bits() == b.identity.to_bits()
            && a.length == b.length
            && a.mismatch == b.mismatch
            && a.gapopen == b.gapopen
            && a.q_start == b.q_start
            && a.q_end == b.q_end
            && a.s_start == b.s_start
            && a.s_end == b.s_end
            && a.gapped_q_start == b.gapped_q_start
            && a.gapped_s_start == b.gapped_s_start
            && a.e_value.to_bits() == b.e_value.to_bits()
            && a.bit_score.to_bits() == b.bit_score.to_bits()
            && a.num_ident == b.num_ident
            && a.query_context == b.query_context
            && a.query_frame == b.query_frame
            && a.subject_frame == b.subject_frame
            && a.query_length == b.query_length
            && a.q_idx == b.q_idx
            && a.s_idx == b.s_idx
            && a.raw_score == b.raw_score
            && a.gap_info == b.gap_info
            && a.num_positives == b.num_positives
            && a.comp_adjustment_method == b.comp_adjustment_method
    }

    /// Panics unless the expanded compact list of query `q_idx` equals the list of full HSPs (all
    /// fields; floating-point values by their bits).
    fn assert_same_hit_lists(
        q_idx: usize,
        full: Option<&BlastpHitList>,
        compact: Option<BlastpHitList<XPrelimHsp>>,
    ) {
        let compact = compact.map(x_expand_hit_list);
        let (full, compact) = match (full, compact.as_ref()) {
            (None, None) => return,
            (Some(full), Some(compact)) => (full, compact),
            (full, compact) => panic!(
                "compact lists, query {q_idx}: full list {}, compact list {}",
                if full.is_some() { "present" } else { "absent" },
                if compact.is_some() {
                    "present"
                } else {
                    "absent"
                },
            ),
        };
        let header = (
            full.hsplist_count,
            full.hsplist_max,
            full.worst_evalue.to_bits(),
            full.low_score,
            full.heapified,
            full.hsplist_array.len(),
            full.hsplist_current,
            full.num_hits,
        );
        let compact_header = (
            compact.hsplist_count,
            compact.hsplist_max,
            compact.worst_evalue.to_bits(),
            compact.low_score,
            compact.heapified,
            compact.hsplist_array.len(),
            compact.hsplist_current,
            compact.num_hits,
        );
        assert_eq!(
            header, compact_header,
            "compact lists, query {q_idx}: hit list fields differ"
        );
        for (position, (a, b)) in full
            .hsplist_array
            .iter()
            .zip(&compact.hsplist_array)
            .enumerate()
        {
            assert!(
                a.oid == b.oid
                    && a.query_index == b.query_index
                    && a.best_evalue.to_bits() == b.best_evalue.to_bits()
                    && a.hsps.len() == b.hsps.len(),
                "compact lists, query {q_idx} list {position}: oid {} / {}, query {} / {}, best E {:e} / {:e}, HSPs {} / {}",
                a.oid,
                b.oid,
                a.query_index,
                b.query_index,
                a.best_evalue,
                b.best_evalue,
                a.hsps.len(),
                b.hsps.len(),
            );
            for (index, (x, y)) in a.hsps.iter().zip(&b.hsps).enumerate() {
                assert!(
                    same_hsp(x, y),
                    "compact lists, query {q_idx} list {position} HSP {index}: {x:?} / {y:?}"
                );
            }
        }
    }

    fn prelim(score: i32, q: (i32, i32), s: (i32, i32), e_value: f64) -> XPrelimHsp {
        XPrelimHsp {
            hsp: BlastpPreliminaryHsp {
                query_start: q.0,
                query_end: q.1,
                subject_start: s.0,
                subject_end: s.1,
                gapped_query_start: q.0 + 1,
                gapped_subject_start: s.0 + 1,
                raw_score: score,
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length: 300,
                q_idx: 0,
                s_idx: 0,
            },
            e_value,
        }
    }

    fn lists_of(oid: u32, hsps: &[XPrelimHsp]) -> (BlastpHspList, BlastpHspList<XPrelimHsp>) {
        let compact = BlastpHspList {
            oid,
            query_index: 0,
            hsps: hsps
                .iter()
                .map(|h| {
                    let mut h = *h;
                    h.hsp.s_idx = oid;
                    h
                })
                .collect::<Vec<_>>(),
            best_evalue: f64::MAX,
        };
        let full = BlastpHspList {
            oid,
            query_index: 0,
            hsps: compact
                .hsps
                .iter()
                .map(|h| super::super::blast_engine::blastp_hsp_from_preliminary(&h.hsp, h.e_value))
                .collect(),
            best_evalue: f64::MAX,
        };
        (full, compact)
    }

    #[test]
    fn offsets_match_the_expanded_hsp() {
        let h = prelim(57, (4, 90), (10, 97), 1e-12);
        let full = super::super::blast_engine::blastp_hsp_from_preliminary(&h.hsp, h.e_value);
        assert_eq!(h.offsets(), full.offsets());
        assert_eq!(h.raw_score(), full.raw_score());
        assert_eq!(h.e_value().to_bits(), full.e_value().to_bits());
    }

    // The hit list helpers on both element types: a list size small enough that the heap,
    // its replacement rule and the in-list E-value sort all run, with ties in E-value and score.
    #[test]
    fn hit_list_helpers_keep_the_same_lists_on_both_types() {
        let mut full_hits: BlastpHitList = BlastpHitList::new(3);
        let mut compact_hits: BlastpHitList<XPrelimHsp> = BlastpHitList::new(3);
        let subjects: Vec<Vec<XPrelimHsp>> = vec![
            vec![
                prelim(40, (0, 50), (0, 50), 1e-5),
                prelim(40, (60, 90), (60, 90), 1e-5),
            ],
            vec![
                prelim(40, (60, 90), (60, 90), 1e-5),
                prelim(40, (0, 50), (0, 50), 1e-5),
            ],
            vec![prelim(55, (1, 80), (3, 82), 1e-9)],
            vec![prelim(40, (0, 50), (0, 50), 1e-5)],
            vec![
                prelim(30, (5, 25), (7, 27), 2e-3),
                prelim(61, (30, 99), (31, 100), 1e-11),
            ],
            vec![
                prelim(55, (2, 81), (3, 82), 1e-9),
                prelim(55, (1, 80), (3, 82), 1e-9),
            ],
            vec![prelim(20, (0, 10), (0, 10), 9.0)],
        ];
        for (oid, hsps) in subjects.iter().enumerate() {
            let (mut full, mut compact) = lists_of(oid as u32, hsps);
            super::super::hsp::reap_hsplist_by_evalue(&mut full, 5.0);
            super::super::hsp::reap_hsplist_by_evalue(&mut compact, 5.0);
            if full.hsps.is_empty() {
                assert!(compact.hsps.is_empty());
                continue;
            }
            full_hits.update(full);
            compact_hits.update(compact);
        }
        full_hits.sort_by_evalue();
        compact_hits.sort_by_evalue();
        assert_eq!(full_hits.hsplist_count, 3);
        assert_same_hit_lists(0, Some(&full_hits), Some(compact_hits));
    }
}
