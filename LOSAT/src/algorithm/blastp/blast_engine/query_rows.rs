//! The rendered rows (`PairwiseHit`) of a call's final hits as its reports read them: a window of
//! queries at a time (`row_windows`), made by the call's producer of rows (`rows_of` in
//! `run_resolved_in_pool`) and dropped before the next window's, as BLASTN's and TBLASTX's outfmt 0
//! make theirs (`x_pairwise_par::query_windows`). Outfmt 0 takes each window's rows whole; BLAST
//! XML and the tabular fields that show the alignment read them query by query in order
//! (`QueryRows`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
//! ```c++
//! CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
//! ```
//! NCBI makes and prints the alignments of one query's result set at a time. A row is a function
//! of its HSP alone (`build_pairwise_hits`), so the rows of a window are the rows that the
//! queries of the window have among the rows of every query, in the same order.

use std::borrow::Cow;
use std::ops::Range;

use anyhow::{bail, Result};

use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::report::pairwise::PairwiseHit;
use crate::utils::threading::SearchPool;

/// The windows of the queries `0..num_queries` whose rows are made together: the windows of
/// `x_pairwise_par::query_windows` over the estimates of the queries' reports (the lengths of
/// the HSPs of their final lists, `x_pairwise_par::query_estimate_of`), made before any row,
/// with the budget of `x_pool` (`x_pairwise_par::window_bytes`). No window without a query.
pub(super) fn row_windows(
    hit_lists: &[Option<BlastpHitList>],
    num_queries: usize,
    x_pool: Option<&SearchPool<'_>>,
) -> Vec<Range<usize>> {
    use crate::report::x_pairwise_par::{query_estimate_of, query_windows, window_bytes};
    let estimates: Vec<usize> = (0..num_queries)
        .map(|q_idx| {
            query_estimate_of(
                hit_lists
                    .get(q_idx)
                    .and_then(Option::as_ref)
                    .into_iter()
                    .flat_map(|hit_list| hit_list.hsplist_array.iter().take(hit_list.hsplist_count))
                    .flat_map(|hsp_list| &hsp_list.hsps)
                    .map(|hsp| hsp.length),
            )
        })
        .collect();
    query_windows(&estimates, window_bytes(x_pool))
}

/// The rows of a report that reads its queries in order (`query`): those of the window that
/// holds the query, made when the report reaches the window (the window before it dropped
/// first), or the rows of every query made up front.
pub(crate) struct QueryRows<'r> {
    /// The windows after the current one.
    windows: std::vec::IntoIter<Range<usize>>,
    /// The queries whose rows `rows` holds.
    window: Range<usize>,
    rows: Cow<'r, [PairwiseHit]>,
    /// The first row of `rows` not handed out yet.
    next: usize,
    /// The query asked last.
    last: Option<usize>,
    rows_of: &'r dyn Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
}

impl<'r> QueryRows<'r> {
    /// The rows of the windows `windows` (`row_windows`) made by `rows_of`, or `up_front`, the
    /// rows of every query when they were made before the reports.
    pub(super) fn new(
        windows: Vec<Range<usize>>,
        up_front: Option<&'r [PairwiseHit]>,
        rows_of: &'r dyn Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
    ) -> Self {
        let (windows, window, rows) = match up_front {
            Some(rows) => (Vec::new(), 0..usize::MAX, Cow::Borrowed(rows)),
            None => (windows, 0..0, Cow::Borrowed(&[][..])),
        };
        QueryRows {
            windows: windows.into_iter(),
            window,
            rows,
            next: 0,
            last: None,
            rows_of,
        }
    }

    /// The rows of query `q_idx`'s HSPs, in the order of its final list. The queries are asked
    /// in increasing order (a query without HSPs may be skipped).
    pub(crate) fn query(&mut self, q_idx: usize) -> Result<&[PairwiseHit]> {
        if self.last.is_some_and(|last| q_idx <= last) {
            bail!("BLASTP report rows asked out of query order");
        }
        self.last = Some(q_idx);
        while q_idx >= self.window.end {
            let Some(window) = self.windows.next() else {
                bail!("BLASTP report rows do not cover the queries");
            };
            // The window before is dropped before the next one is made.
            self.rows = Cow::Borrowed(&[]);
            self.rows = Cow::Owned((self.rows_of)(window.clone())?);
            self.window = window;
            self.next = 0;
        }
        let rows = &self.rows[self.next..];
        let skipped = rows
            .iter()
            .take_while(|row| (row.hit.q_idx as usize) < q_idx)
            .count();
        let own = rows[skipped..]
            .iter()
            .take_while(|row| row.hit.q_idx as usize == q_idx)
            .count();
        let start = self.next + skipped;
        self.next = start + own;
        Ok(&self.rows[start..start + own])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Hit;

    /// A row of query `q_idx` against subject `s_idx`.
    fn row(q_idx: u32, s_idx: u32) -> PairwiseHit {
        PairwiseHit::from(Hit {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            e_value: 1e-5,
            bit_score: 25.0,
            num_ident: 10,
            query_frame: 0,
            query_length: 30,
            q_idx,
            s_idx,
            raw_score: 50,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 10,
        })
    }

    /// The (query, subject) pairs of `rows`.
    fn pairs(rows: &[PairwiseHit]) -> Vec<(u32, u32)> {
        rows.iter()
            .map(|row| (row.hit.q_idx, row.hit.s_idx))
            .collect()
    }

    /// Queries 0 and 3 without rows; 1 with two, 2 with one, 4 with three.
    fn all_rows() -> Vec<PairwiseHit> {
        vec![
            row(1, 5),
            row(1, 2),
            row(2, 7),
            row(4, 0),
            row(4, 9),
            row(4, 1),
        ]
    }

    /// The rows of the queries of `window` among `all_rows`, as the call's producer makes them.
    fn rows_of(window: Range<usize>) -> Result<Vec<PairwiseHit>> {
        Ok(all_rows()
            .into_iter()
            .filter(|row| window.contains(&(row.hit.q_idx as usize)))
            .collect())
    }

    #[test]
    fn each_query_gets_its_rows_from_windows_or_from_the_rows_made_up_front() {
        let all = all_rows();
        let expected: Vec<Vec<(u32, u32)>> = (0..5)
            .map(|q_idx| {
                pairs(
                    &all.iter()
                        .filter(|row| row.hit.q_idx as usize == q_idx)
                        .cloned()
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        for windows in [
            vec![0..5],
            vec![0..2, 2..5],
            vec![0..1, 1..2, 2..3, 3..4, 4..5],
            vec![0..3, 3..5],
        ] {
            for up_front in [None, Some(all.as_slice())] {
                // Every query asked.
                let mut rows = QueryRows::new(windows.clone(), up_front, &rows_of);
                for (q_idx, expected) in expected.iter().enumerate() {
                    assert_eq!(&pairs(rows.query(q_idx).unwrap()), expected, "{windows:?}");
                }
                // Only the queries with rows (BLAST XML, the tabular report).
                let mut rows = QueryRows::new(windows.clone(), up_front, &rows_of);
                for q_idx in [1, 2, 4] {
                    assert_eq!(pairs(rows.query(q_idx).unwrap()), expected[q_idx]);
                }
            }
        }
    }

    #[test]
    fn a_query_asked_again_or_past_the_windows_is_an_error() {
        let mut rows = QueryRows::new(vec![0..2, 2..5], None, &rows_of);
        rows.query(2).unwrap();
        assert!(rows.query(1).is_err());
        let mut rows = QueryRows::new(vec![0..2], None, &rows_of);
        rows.query(1).unwrap();
        assert!(rows.query(1).is_err());
        assert!(rows.query(4).is_err());
    }
}
