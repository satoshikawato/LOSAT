//! The rendered rows (`PairwiseHit`) of a call's final hits as its reports read them: a window of
//! queries at a time (`row_windows`), made by the call's producer of rows (`rows_of` in
//! `run_resolved_in_pool`) and dropped before the next window's, as BLASTN's and TBLASTX's outfmt 0
//! make theirs (`x_pairwise_par::query_windows`). Outfmt 0 takes each window's rows whole; the
//! tabular fields that show the alignment read them query by query in order (`QueryRows`), and so
//! does BLAST XML, whose windows are single queries.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
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
/// the HSPs of their final lists in a report of `line_length` columns a line,
/// `x_pairwise_par::query_estimate_of`), made before any row, with the budget of `x_pool`
/// (`x_pairwise_par::window_bytes`). No window without a query.
pub(super) fn row_windows(
    hit_lists: &[Option<BlastpHitList>],
    num_queries: usize,
    x_pool: Option<&SearchPool<'_>>,
    line_length: usize,
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
                line_length,
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
#[path = "query_rows_tests.rs"]
mod tests;
