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

    /// No NCBI counterpart: EXPERIMENT (LOSAT_X_DCEMITHOIST) the form as a value whose `get` (the
    /// same arithmetic as `get`'s, the same answer for every offset) has its constants made once,
    /// for a loop that keeps them in registers and is compiled once per form
    /// (`x_pair_seed_disc.rs` `DcPairs`) instead of choosing the form and reading its fields for
    /// every offset.
    pub(super) fn x_view(&self) -> XView<'_> {
        match self {
            CtxO1::Div {
                magic,
                room,
                len,
                total_len,
            } => XView::Div(XDiv {
                magic: *magic,
                room: *room,
                len: *len,
                total_len: *total_len,
            }),
            CtxO1::Blocks {
                shift,
                entries,
                rank_context,
                total_len,
                ..
            } => {
                if rank_context.is_empty() {
                    XView::Blocks(XBlocks::new(*shift, entries, rank_context, *total_len))
                } else {
                    XView::BlocksRanks(XBlocks::new(*shift, entries, rank_context, *total_len))
                }
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

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_DCEMITHOIST) the answer of `CtxO1::get` for one form
/// with its constants made once (`CtxO1::x_view`), or of no form (`XNone`: always `None`, so the
/// caller asks its other forms, as `QueryContextIndex::context_for_offset` does without `x_o1`).
pub(super) trait XCtxGet: Copy {
    /// `CtxO1::get(n)` of the form the value was made from.
    fn get(&self, n: usize) -> Option<usize>;
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_DCEMITHOIST) `CtxO1` as a value (`CtxO1::x_view`), one
/// variant per loop its user compiles.
#[derive(Clone, Copy)]
pub(super) enum XView<'a> {
    Div(XDiv),
    /// `Blocks` where every context has letters (no `rank_context`).
    Blocks(XBlocks<'a, false>),
    /// `Blocks` with contexts without letters (ranks mapped through `rank_context`).
    BlocksRanks(XBlocks<'a, true>),
}

/// `CtxO1::Div` as a value.
#[derive(Clone, Copy)]
pub(super) struct XDiv {
    magic: u64,
    room: u64,
    len: u64,
    total_len: usize,
}

impl XCtxGet for XDiv {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:225-227
    /// `CtxO1::get` of `Div`, statement for statement.
    #[inline(always)]
    fn get(&self, n: usize) -> Option<usize> {
        if n >= self.total_len {
            return None;
        }
        let q = ((self.magic as u128 * n as u128) >> 64) as u64;
        let r = n as u64 - q * self.room;
        Some(if r == self.len { 0 } else { q as usize })
    }
}

/// `CtxO1::Blocks` as a value, with the masks `get` makes from `shift` made once; `RANKS` whether
/// ranks go through `rank_context` (some context has no letters).
#[derive(Clone, Copy)]
pub(super) struct XBlocks<'a, const RANKS: bool> {
    entries: &'a [u32],
    rank_context: &'a [u32],
    total_len: usize,
    shift: u32,
    low_mask: u32,
    crowded_place: u32,
}

impl<'a, const RANKS: bool> XBlocks<'a, RANKS> {
    /// The value of a `Blocks` form. It checks what makes `get`'s unchecked read of `entries`
    /// sound: every offset below `total_len` has its block's entry (`build_with` makes one entry
    /// per block of `2^shift` offsets up to `total_len`), and the masks fit in 32 bits.
    fn new(shift: u32, entries: &'a [u32], rank_context: &'a [u32], total_len: usize) -> Self {
        assert!(shift < 31 && entries.len() >= total_len.div_ceil(1usize << shift));
        debug_assert_eq!(RANKS, !rank_context.is_empty());
        XBlocks {
            entries,
            rank_context,
            total_len,
            shift,
            low_mask: (1u32 << shift) - 1,
            crowded_place: (2u32 << shift) - 1,
        }
    }
}

impl<const RANKS: bool> XCtxGet for XBlocks<'_, RANKS> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:225-227
    /// `CtxO1::get` of `Blocks`, statement for statement, with the entry read unchecked.
    #[inline(always)]
    fn get(&self, n: usize) -> Option<usize> {
        if n >= self.total_len {
            return None;
        }
        debug_assert!(n >> self.shift < self.entries.len());
        // SAFETY: `n < total_len`, and `new` checked that `entries` has the entry of every block of
        // `2^shift` offsets below `total_len`, so `n >> shift < entries.len()`.
        let entry = unsafe { *self.entries.get_unchecked(n >> self.shift) };
        let low = (n as u32) & self.low_mask;
        let place = entry & self.crowded_place;
        if place == self.crowded_place {
            return None;
        }
        let rank = (entry >> (self.shift + 1)) + u32::from(low > place);
        let context = if RANKS {
            self.rank_context[rank as usize]
        } else {
            rank
        };
        Some(if low == place { 0 } else { context as usize })
    }
}

/// No form (`QueryContextIndex` without `x_o1`): `get` answers `None` for every offset.
#[derive(Clone, Copy)]
pub(super) struct XNone;

impl XCtxGet for XNone {
    #[inline(always)]
    fn get(&self, _n: usize) -> Option<usize> {
        None
    }
}

#[cfg(test)]
#[path = "x_ctx_o1_tests.rs"]
mod tests;
