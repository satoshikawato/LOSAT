//! EXPERIMENT (LOSAT_X_CTXO1, a strict-set switch): the context of a query offset in constant time, from the layout of the
//! concatenated query, instead of a binary search over the contexts (the runs of
//! `QueryContextIndex`) or a table of four bytes per query letter (`direct_map`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-242
//! ```c
//! Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
//! {
//!     Int4 m=0, b=0, e=0, size=0;
//!
//!     size = A->last_context+1;
//!
//!     if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
//!         b = MIN(n / (A->max_length + 1), size - 1);
//!         e = MIN(n / (A->min_length + 1) + 1, size);
//!         ASSERT(e <= size);
//!     }
//!     else {
//!         b = 0;
//!         e = size;
//!     }
//!
//!     while (b < e - 1) {
//! 	m = (b + e) / 2;
//! 	if (A->contexts[m].query_offset > n)
//! 	    e = m;
//! 	else
//! 	    b = m;
//!     }
//!     return b;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-77
//! ```c
//!     if (index) {
//!         Uint4 prev_loc = qinfo->contexts[index-1].query_offset;
//!         Uint4 prev_len = qinfo->contexts[index-1].query_length;
//!
//!         Uint4 shift = prev_len ? prev_len + 1 : 0;
//!
//!         qinfo->contexts[index].query_offset = prev_loc + shift;
//! ```
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
    /// ```c
    ///         Uint4 shift = prev_len ? prev_len + 1 : 0;
    ///
    ///         qinfo->contexts[index].query_offset = prev_loc + shift;
    /// ```
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
    /// ```c
    ///     if (A->min_length > 0 && A->max_length > 0 && A->first_context == 0) {
    ///         b = MIN(n / (A->max_length + 1), size - 1);
    ///         e = MIN(n / (A->min_length + 1) + 1, size);
    /// ```
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
mod tests {
    use super::*;

    fn contexts_of(lengths: &[usize]) -> (Vec<QueryContext>, usize) {
        // the layout of run.rs: a context after the one before and its sentinel, none for no letters
        let mut out = Vec::new();
        let mut offset = 0usize;
        for (i, &len) in lengths.iter().enumerate() {
            out.push(QueryContext {
                query_idx: (i / 2) as u32,
                frame: if i % 2 == 0 { 1 } else { -1 },
                query_offset: offset as i32,
                len,
                masks: Vec::new(),
            });
            offset += if len == 0 { 0 } else { len + 1 };
        }
        let total = out
            .iter()
            .map(|c| c.query_offset as usize + c.len)
            .max()
            .unwrap_or(0);
        (out, total)
    }

    // the table of QueryContextIndex::new
    fn direct_map(contexts: &[QueryContext], total: usize) -> Vec<u32> {
        let mut map = vec![0u32; total];
        for (idx, ctx) in contexts.iter().enumerate() {
            let start = ctx.query_offset.max(0) as usize;
            let end = start + ctx.len;
            if start < end && end <= total {
                map[start..end].fill(idx as u32);
            }
        }
        map
    }

    // every offset: the table's value, or `None` in a block with two sentinels or more
    fn check(lengths: &[usize], force_blocks: bool) -> &'static str {
        let (contexts, total) = contexts_of(lengths);
        let form = CtxO1::build_with(&contexts, total, force_blocks).expect("form");
        let (name, _, shift) = form.describe();
        let map = direct_map(&contexts, total);
        let sentinels: Vec<usize> = contexts
            .iter()
            .filter(|c| c.len > 0)
            .map(|c| c.query_offset as usize + c.len)
            .collect();
        let mut unanswered = 0usize;
        for (n, &want) in map.iter().enumerate() {
            match form.get(n) {
                Some(got) => assert_eq!(
                    got, want as usize,
                    "offset {n} of lengths {lengths:?} (force_blocks {force_blocks})"
                ),
                None => {
                    let block = n >> shift;
                    let in_block = sentinels.iter().filter(|&&e| e >> shift == block).count();
                    assert!(
                        name == "blocks" && form.needs_runs() && in_block >= 2,
                        "offset {n}"
                    );
                    unanswered += 1;
                }
            }
        }
        assert_eq!(form.needs_runs(), unanswered > 0);
        assert_eq!(form.get(total), None);
        assert_eq!(form.get(total + 1000), None);
        name
    }

    fn lcg(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state >> 33
    }

    #[test]
    fn equal_lengths_divide_and_blocks_agree_with_the_table() {
        for len in [1usize, 2, 3, 7, 15, 16, 17, 24, 31, 32, 33, 150, 1000] {
            for queries in [1usize, 2, 3, 100] {
                let lengths = vec![len; 2 * queries];
                assert_eq!(check(&lengths, false), "div");
                assert_eq!(check(&lengths, true), "blocks");
            }
        }
    }

    #[test]
    fn mixed_lengths_and_empty_contexts_agree_with_the_table() {
        let mut state = 0x5eed_u64;
        for round in 0..400 {
            let queries = 1 + (lcg(&mut state) % 40) as usize;
            let (lo, span) = match round % 4 {
                0 => (18usize, 13usize),
                1 => (1, 4),
                2 => (1, 300),
                _ => (60, 2000),
            };
            let mut lengths = Vec::new();
            for _ in 0..queries {
                let len = if lcg(&mut state) % 10 == 0 {
                    0
                } else {
                    lo + (lcg(&mut state) as usize % span)
                };
                lengths.push(len);
                lengths.push(len);
            }
            if lengths.iter().all(|&l| l == 0) {
                lengths[0] = 5;
                lengths[1] = 5;
            }
            let form = check(&lengths, false);
            let equal = lengths.iter().all(|&l| l == lengths[0]);
            assert_eq!(form, if equal { "div" } else { "blocks" }, "{lengths:?}");
        }
    }

    #[test]
    fn tight_shortest_lengths_agree_with_the_table() {
        // the shortest context sets the block: lengths 2^k - 1 and 2^k - 2 are the tight cases
        for k in 1..12u32 {
            let short = (1usize << k) - 1;
            for lengths in [
                vec![short, short + 5, short, 3 * short + 1],
                vec![short + 1, short, short, short + 2, short + 9, short],
                vec![0, 0, short, short + 1, 0, 0, short, short],
            ] {
                assert_eq!(check(&lengths, false), "blocks");
            }
        }
    }

    #[test]
    fn a_short_context_among_long_ones_takes_the_long_ones_block_size() {
        let mut state = 0x0b10c_u64;
        let mut lengths = Vec::new();
        for i in 0..200 {
            let len = if i == 77 {
                5
            } else {
                300 + (lcg(&mut state) as usize % 400)
            };
            lengths.push(len);
            lengths.push(len);
        }
        assert_eq!(check(&lengths, false), "blocks");
        let (contexts, total) = contexts_of(&lengths);
        let form = CtxO1::build_with(&contexts, total, false).unwrap();
        // the shortest context alone gives blocks of 4; 256 still holds every context but two
        assert_eq!(form.describe().2, 8);
        assert!(form.needs_runs());
        // without the short one, the same blocks and none crowded
        lengths[154] = 300;
        lengths[155] = 300;
        let (contexts, total) = contexts_of(&lengths);
        let form = CtxO1::build_with(&contexts, total, false).unwrap();
        assert_eq!(form.describe().2, 8);
        assert!(!form.needs_runs());
    }

    #[test]
    fn a_short_last_context_leaves_the_last_block_uncrowded() {
        // the end of the last context is the table's length, not a sentinel inside it: the block
        // that holds the sentinel before it and that end has one sentinel
        let mut lengths = vec![400usize; 600];
        lengths.push(20);
        assert_eq!(check(&lengths, false), "blocks");
        let (contexts, total) = contexts_of(&lengths);
        let form = CtxO1::build_with(&contexts, total, false).unwrap();
        assert_eq!(form.describe().2, 8);
        assert!(!form.needs_runs());
        // with contexts without letters after it, that end is a sentinel (0) inside the table
        lengths.extend([0, 0]);
        assert_eq!(check(&lengths, false), "blocks");
    }

    #[test]
    fn contexts_without_letters_at_the_end_add_the_last_sentinel() {
        for lengths in [
            vec![7usize, 7, 0, 0],
            vec![5, 9, 0, 0, 0, 0],
            vec![0, 0, 30, 30, 0, 0],
        ] {
            let (contexts, total) = contexts_of(&lengths);
            let ends = contexts
                .iter()
                .filter(|c| c.len > 0)
                .map(|c| c.query_offset as usize + c.len)
                .max()
                .unwrap();
            assert_eq!(total, ends + 1, "{lengths:?}");
            assert_eq!(check(&lengths, false), "blocks");
        }
    }

    #[test]
    fn other_layouts_are_left_to_the_previous_form() {
        // a gap of two letters between contexts
        let (mut contexts, total) = contexts_of(&[10, 10, 12, 12]);
        contexts[2].query_offset += 1;
        contexts[3].query_offset += 1;
        assert!(CtxO1::build_with(&contexts, total + 1, false).is_none());
        // overlapping contexts
        let (mut contexts, total) = contexts_of(&[10, 10]);
        contexts[1].query_offset -= 3;
        assert!(CtxO1::build_with(&contexts, total - 3, false).is_none());
        // no letters at all
        let (contexts, _) = contexts_of(&[0, 0]);
        assert!(CtxO1::build_with(&contexts, 0, false).is_none());
        // a length that does not end at the last context (or one past it, see the next test)
        let (contexts, total) = contexts_of(&[10, 10]);
        assert!(CtxO1::build_with(&contexts, total + 2, false).is_none());
        assert!(CtxO1::build_with(&contexts, total - 1, false).is_none());
    }
}
