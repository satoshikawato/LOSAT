//! The context of an offset of a concatenated query (NCBI's `BSearchContextInfo`) without the
//! binary search over the contexts (`XContextTable`). blastp's and blastx's fused two-hit loops
//! (`blastp/x_seed_fused.rs`, `blastx/x_seed_diag.rs`, `blastx/x_seed_scan.rs`) and TBLASTX's unit
//! path (`tblastx/x_pair_par.rs`, `LOSAT_X_TBXCTXTABLE`) take it. It lives here, below every
//! program, because it reads only the frame bases of the query's contexts.

/// The context of an offset of the concatenated query without the
/// binary search over the contexts: a coarse table gives the context at the start of the offset's
/// block, and a walk over the (small, contiguous) frame bases moves on to the last context starting
/// at or before the offset. The result is the largest context index whose `frame_base` is at most
/// the offset (0 below the first), which is what the binary search of
/// `blastp_context_idx_for_query_offset_with_bounds` returns. The coarse table has at most
/// `COARSE_ENTRIES` entries so that both arrays stay in the first cache levels: a table with one
/// entry per offset (as large as the diagonal table) cost a cache miss for every candidate hit.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:211-242
pub(crate) struct XContextTable {
    /// Frame bases in context order, then `i32::MAX`.
    starts: Vec<i32>,
    /// Context index at offset `b << shift`.
    coarse: Vec<u32>,
    shift: u32,
    span: usize,
}

impl XContextTable {
    const COARSE_ENTRIES: usize = 1 << 14;

    /// `frame_bases` in context order (non-decreasing); `span` = the offsets to cover.
    pub(crate) fn new(frame_bases: &[i32], span: usize) -> Self {
        let mut shift = 0u32;
        while (span >> shift) >= Self::COARSE_ENTRIES {
            shift += 1;
        }
        Self::with_shift(frame_bases, span, shift)
    }

    pub(crate) fn with_shift(frame_bases: &[i32], span: usize, shift: u32) -> Self {
        let mut starts = Vec::with_capacity(frame_bases.len() + 2);
        starts.extend_from_slice(frame_bases);
        starts.push(i32::MAX);
        if frame_bases.is_empty() {
            // no context: every offset gives index 0, as the binary search does
            starts.push(i32::MAX);
        }
        let blocks = (span >> shift) + 1;
        let mut coarse = Vec::with_capacity(blocks);
        let mut idx = 0usize;
        for b in 0..blocks {
            let offset = i64::try_from(b << shift).expect("query offset must fit i64");
            while idx + 1 < frame_bases.len() && i64::from(starts[idx + 1]) <= offset {
                idx += 1;
            }
            coarse.push(u32::try_from(idx).expect("context index must fit u32"));
        }
        Self {
            starts,
            coarse,
            shift,
            span,
        }
    }

    /// The context of `query_offset`, or `None` outside the table.
    #[inline(always)]
    pub(crate) fn get(&self, query_offset: i32) -> Option<usize> {
        let q = usize::try_from(query_offset)
            .ok()
            .filter(|&q| q < self.span)?;
        let mut idx = self.coarse[q >> self.shift] as usize;
        // starts[last] = i32::MAX ends the walk; the coarse entry is a context index, so idx + 1
        // is always inside `starts`.
        while self.starts[idx + 1] <= query_offset {
            idx += 1;
        }
        Some(idx)
    }

    /// [`XContextTable::get`] and the context's frame base (`starts[idx]`, the table's copy of
    /// `contexts[idx].frame_base`), so that the caller's test of the frame base reads the line the
    /// search just read instead of the context's record.
    #[inline(always)]
    pub(crate) fn get_with_start(&self, query_offset: i32) -> Option<(usize, i32)> {
        let q = usize::try_from(query_offset)
            .ok()
            .filter(|&q| q < self.span)?;
        let mut idx = self.coarse[q >> self.shift] as usize;
        while self.starts[idx + 1] <= query_offset {
            idx += 1;
        }
        Some((idx, self.starts[idx]))
    }
}
