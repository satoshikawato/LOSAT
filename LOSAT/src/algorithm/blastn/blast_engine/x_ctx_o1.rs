//! LOSAT_X_CTXO1 (graduated 2026-10-27, always on): the context of a query offset in constant time,
//! from the layout of the concatenated query, instead of a binary search over the contexts (the
//! runs of `QueryContextIndex`) or a table of four bytes per query letter (`direct_map`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-242
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-77
//!
//! The value answered is `direct_map[n]` of `QueryContextIndex` for every `n` below the
//! concatenated length: the context whose letters hold offset `n`, and 0 for a sentinel. The
//! layout above puts every context with letters right after the one before and its sentinel, and
//! gives the contexts without letters no room. `build` checks that layout (the first context with
//! letters at 0, each next one at the end of the one before plus one) and returns `None`
//! otherwise; the caller then keeps its previous form. Two forms:
//!
//! - `Div`, when every context has the same number of letters `L` (none empty): the context at
//!   `n` is `n / (L + 1)` and the sentinel is the remainder `L`, the shortcut of NCBI above (there
//!   `b = n / (len + 1)` and `e = b + 1`, so the loop does not run). The division is a
//!   multiplication by `ceil(2^64 / (L + 1))` and a shift by 64, exact for every `n < 2^32`
//!   (Lemire, Kaser and Kurz, "Faster remainder by direct computation", 2019, theorem 1 with
//!   F = 64, N = 32); `build` takes only concatenated lengths below 2^31.
//! - `Blocks`, for any other layout: one `u32` per block of `2^k` offsets. A block with at most one
//!   sentinel has the entry `rank << (k + 1) | p`: `rank` numbers (from 0, in offset order) the
//!   context with letters whose letters or sentinel hold the block's first offset, and `p` is the
//!   sentinel's place in the block, or `2^k` if the block has none. Offset `n` with
//!   `low = n mod 2^k` is the sentinel if `low == p`, in context `rank` if `low < p`, and in context
//!   `rank + 1` if `low > p` (the context after a sentinel starts right after it). A block with two
//!   sentinels or more (it holds a whole context shorter than `2^k - 1` letters) has `p = 2^(k+1) - 1`
//!   ("crowded"), and `get` answers `None` there: the caller asks its runs, which give the same
//!   value. The rank is the context's index unless some context has no letters; then a second
//!   table maps ranks to indices. `k` is the largest value (up to 30) for which the rank fits in the
//!   `32 - (k + 1)` upper bits and at most one block in `CROWDED_SHARE` is crowded; it is at least
//!   `log2` of the shortest context plus one rounded down, where no block is crowded (two sentinels
//!   are a context and one apart).
//!
//! Neither form reads more than one `u32` (the second table only with empty contexts), where the
//! binary search reads up to 16 entries of a 1.4 MB array for 30,000 queries (R8: 318 cycles a
//! lookup hit in blastn-short) and `direct_map` reads one of 4 bytes a letter. The block table
//! takes `4 * len / 2^k` bytes (375 KB for 30,000 queries of 18 to 30 letters, k = 4; 2.8 MB for both
//! strands of 31,168 cDNAs of 30 to 71,382 letters, k = 8, where the shortest one alone would give
//! k = 4 and 45 MB); a table above `MAX_TABLE_BYTES` is not built.

use super::QueryContext;

/// No NCBI counterpart: the largest block table `build` makes (bytes).
const MAX_TABLE_BYTES: usize = 64 << 20;

/// No NCBI counterpart: a block size is taken while at most one block in this many is crowded.
const CROWDED_SHARE: usize = 256;

/// No NCBI counterpart: the constant-time form of `direct_map` (see the module header).
pub(super) enum CtxO1 {
    /// Every context has `len` letters: context `n / (len + 1)`, 0 at the sentinels.
    Div {
        magic: u64,
        room: u64,
        len: u64,
        total_len: usize,
    },
    /// One entry per block of `2^shift` offsets.
    Blocks {
        shift: u32,
        entries: Vec<u32>,
        /// Index of the context of each rank, empty when every context has letters.
        rank_context: Vec<u32>,
        total_len: usize,
        /// Some block is crowded (`get` answers `None` there).
        crowded: bool,
    },
}

/// No NCBI counterpart: reads `LOSAT_X_CTXO1FORM` once (tuning, not part of the strict set):
/// `blocks` builds the block table even where every context has the same length.
fn x_force_blocks() -> bool {
    use std::sync::OnceLock;
    static FORCE: OnceLock<bool> = OnceLock::new();
    *FORCE.get_or_init(|| {
        std::env::var("LOSAT_X_CTXO1FORM")
            .map(|v| v.trim() == "blocks")
            .unwrap_or(false)
    })
}

impl CtxO1 {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-77
    /// The form for `contexts` (`total_len` = the end of the last context with letters), or
    /// `None` if they are not laid out as above, have no letters, or the table would be too big.
    pub(super) fn build(contexts: &[QueryContext], total_len: usize) -> Option<Self> {
        Self::build_with(contexts, total_len, x_force_blocks())
    }

    fn build_with(contexts: &[QueryContext], total_len: usize, force_blocks: bool) -> Option<Self> {
        if total_len == 0 || total_len >= 1usize << 31 {
            return None;
        }
        // the contexts with letters in index order, checked to follow each other with one sentinel;
        // `ends[r]` is the sentinel after context `r` (after the last one, the end of the letters)
        let mut ends: Vec<usize> = Vec::new();
        let mut rank_context: Vec<u32> = Vec::new();
        let mut any_empty = false;
        let mut min_len = usize::MAX;
        let mut max_len = 0usize;
        let mut next = 0usize;
        for (idx, ctx) in contexts.iter().enumerate() {
            if ctx.len == 0 {
                any_empty = true;
                continue;
            }
            if ctx.query_offset < 0 || ctx.query_offset as usize != next {
                return None;
            }
            ends.push(next + ctx.len);
            rank_context.push(idx as u32);
            min_len = min_len.min(ctx.len);
            max_len = max_len.max(ctx.len);
            next += ctx.len + 1;
        }
        // `total_len` reaches one past the last letters when contexts without letters follow
        // (their offset is that position); the table holds 0 there, as at a sentinel
        if ends.is_empty() || (total_len != next - 1 && total_len != next) {
            return None;
        }
        if !any_empty && min_len == max_len && !force_blocks {
            let room = (min_len + 1) as u64;
            return Some(CtxO1::Div {
                magic: u64::MAX / room + 1,
                room,
                len: min_len as u64,
                total_len,
            });
        }
        let last = ends.len() - 1;
        let rank_fits = |shift: u32| last < 1usize << (31 - shift);
        // blocks of 2^shift offsets with two sentinels or more below `total_len` (an upper bound: a
        // block with three counts twice; `ends[last]` is past the table unless empty contexts trail)
        let crowded_blocks = |shift: u32| {
            ends.windows(2)
                .filter(|w| w[1] < total_len && w[0] >> shift == w[1] >> shift)
                .count()
        };
        // 2^shift <= min_len + 1: at most one sentinel per block
        let mut shift = (usize::BITS - 1 - (min_len + 1).leading_zeros()).min(30);
        if !rank_fits(shift) {
            return None;
        }
        while shift < 30
            && rank_fits(shift + 1)
            && crowded_blocks(shift + 1)
                <= total_len.div_ceil(1usize << (shift + 1)) / CROWDED_SHARE
        {
            shift += 1;
        }
        let block = 1usize << shift;
        let blocks = total_len.div_ceil(block);
        if blocks.checked_mul(4)? > MAX_TABLE_BYTES {
            return None;
        }
        let crowded_place = (2u32 << shift) - 1;
        let mut crowded = false;
        let mut entries = Vec::with_capacity(blocks);
        let mut rank = 0usize;
        for j in 0..blocks {
            let first = j << shift;
            // the context whose letters or sentinel hold `first` (context `rank + 1` starts at
            // `ends[rank] + 1`; `first < total_len <= ends[last] + 1`)
            while rank < last && ends[rank] < first {
                rank += 1;
            }
            debug_assert!(ends[rank] >= first);
            let place = if rank < last && ends[rank + 1] < (first + block).min(total_len) {
                crowded = true;
                crowded_place
            } else {
                (ends[rank] - first).min(block) as u32
            };
            entries.push(((rank as u32) << (shift + 1)) | place);
        }
        if !any_empty {
            rank_context = Vec::new();
        }
        Some(CtxO1::Blocks {
            shift,
            entries,
            rank_context,
            total_len,
            crowded,
        })
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:225-227
    /// `direct_map[n]` for `n` below the concatenated length, `None` above it and in a crowded
    /// block.
    #[inline(always)]
    pub(super) fn get(&self, n: usize) -> Option<usize> {
        match self {
            CtxO1::Div {
                magic,
                room,
                len,
                total_len,
            } => {
                if n >= *total_len {
                    return None;
                }
                let q = ((*magic as u128 * n as u128) >> 64) as u64;
                let r = n as u64 - q * *room;
                Some(if r == *len { 0 } else { q as usize })
            }
            CtxO1::Blocks {
                shift,
                entries,
                rank_context,
                total_len,
                ..
            } => {
                if n >= *total_len {
                    return None;
                }
                let entry = entries[n >> *shift];
                let low = (n as u32) & ((1u32 << *shift) - 1);
                let crowded_place = (2u32 << *shift) - 1;
                let place = entry & crowded_place;
                if place == crowded_place {
                    return None;
                }
                let rank = (entry >> (*shift + 1)) + u32::from(low > place);
                let context = if rank_context.is_empty() {
                    rank
                } else {
                    rank_context[rank as usize]
                };
                Some(if low == place { 0 } else { context as usize })
            }
        }
    }

    /// No NCBI counterpart: whether `get` answers `None` for some offset below the length (the
    /// caller then keeps its runs).
    pub(super) fn needs_runs(&self) -> bool {
        matches!(self, CtxO1::Blocks { crowded: true, .. })
    }

    /// No NCBI counterpart: the form's name, bytes and block shift (tests).
    #[cfg(test)]
    pub(super) fn describe(&self) -> (&'static str, usize, u32) {
        match self {
            CtxO1::Div { .. } => ("div", 0, 0),
            CtxO1::Blocks {
                shift,
                entries,
                rank_context,
                ..
            } => ("blocks", 4 * (entries.len() + rank_context.len()), *shift),
        }
    }
}

#[cfg(test)]
#[path = "x_ctx_o1_tests.rs"]
mod tests;
