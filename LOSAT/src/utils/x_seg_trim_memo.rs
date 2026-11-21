//! LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26, always on): the trim of a SEG segment
//! computed once for each segment's bytes, for the maskers that ask for it.
//!
//! NCBI trims every candidate segment that `s_SegSeq` finds, by evaluating `s_GetProb` on every
//! window of every length from the segment's length down to `max(1, length - maxtrim) + 1`:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1988-1992,1998-2001,2007-2008,2013-2014
//!
//! The pair (lend, rend) is a function of the segment's letters, `maxtrim` and the alphabet;
//! no position, no entropy value and no other state of the caller enters it. The fast path
//! (`SegMasker::x_mask_sequence`) always uses the 20-letter alphabet.
//!
//! Two maskers repeat segments and use this memo (`SegMasker::x_with_trim_memo`):
//! - TBLASTN's composition-based redo masks a window of the translated subject around each
//!   alignment (blast_kappa.c:1518-1541), so a low-complexity stretch of the subject is trimmed
//!   again for every query whose window covers it (v4-C DECISIONS 31: 118,765 trims of 12,584
//!   distinct segments on AP027131 x AP027133);
//! - BLASTX masks the frames of a long query and then, again, the frames of each of its chunks
//!   (`split.rs`), whose segments are mostly those of the whole query (51 % of the trim work on
//!   AP027131 x AP027133).
//!
//! Other maskers repeat few segments (TBLASTX d04: 2 %) and do not use it, so their trims
//! are not stored.
//!
//! This module keeps, for all threads, the bytes of the segments trimmed and their pairs, in one
//! table behind one lock. A call whose bytes and `maxtrim` equal a stored entry (the whole slice
//! is compared; the hash only selects the slot) returns the stored pair; any other call evaluates
//! the trim (without the lock) and stores it. A segment one thread trimmed is thus not trimmed
//! again by another (only threads that miss it at the same moment each trim it, and a cleared
//! table trims it again), and the threads share one copy of the table, which lasts as long as the
//! process. The table grows on demand up to `X_MAX_ENTRIES` entries (or `X_MAX_ARENA` letters),
//! the bounds one thread's table had, and is then cleared as a whole; one thread fills and clears
//! it as before.
//!
//! No NCBI counterpart: the memo itself; it returns only values the trim computed for the same
//! input and does not change any value NCBI computes.

use std::sync::Mutex;

// No NCBI counterpart: the size bounds of the memo; they decide only which trims are evaluated
// again and do not change any value NCBI computes.
/// Most entries kept; the table is cleared when it would hold more. 32,768 keeps every distinct
/// segment of the measured jobs (12,584 for TBLASTN AP027131 x AP027133, 25,507 for BLASTX
/// AP027131 x AP027133; v4-C DECISIONS 32).
const X_MAX_ENTRIES: usize = 32_768;
/// Most letters kept (the stored segments' bytes); the table is cleared when it would hold more.
const X_MAX_ARENA: usize = 1 << 20;
/// Slots of a new table; it doubles while it is more than half full, up to twice
/// `X_MAX_ENTRIES` (65,536 slots of 12 bytes).
const X_FIRST_SLOTS: usize = 1_024;
/// Longest segment stored: the fast trim's own bound. Longer segments are trimmed every time.
const X_MAX_SEGMENT: usize = crate::utils::seg::X_TRIM_MAX;
const _: () = assert!(
    X_MAX_SEGMENT <= u8::MAX as usize
        && X_FIRST_SLOTS.is_power_of_two()
        && X_MAX_ARENA <= u32::MAX as usize
);

/// One stored trim. `len == 0` marks an empty slot (a stored segment has at least one letter).
#[derive(Clone, Copy, Default)]
struct Slot {
    /// The high 32 bits of the segment's hash.
    tag: u32,
    /// Offset of the segment's letters in the arena.
    start: u32,
    len: u8,
    maxtrim: u8,
    lend: u8,
    rend: u8,
}

/// The stored trims of one table.
struct Memo {
    slots: Vec<Slot>,
    arena: Vec<u8>,
    entries: usize,
}

/// No NCBI counterpart: the slot hash of a segment (64-bit multiply-rotate over 8-byte words);
/// it only places entries and does not change any value NCBI computes.
#[inline]
fn hash_of(segment: &[u8], maxtrim: u8) -> u64 {
    const K: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut h = ((segment.len() as u64) << 8 | u64::from(maxtrim)).wrapping_mul(K);
    let mut chunks = segment.chunks_exact(8);
    for chunk in &mut chunks {
        let word = u64::from_le_bytes(chunk.try_into().expect("8-byte chunk"));
        h = (h.rotate_left(23) ^ word).wrapping_mul(K);
    }
    let mut tail = 0u64;
    for (i, &byte) in chunks.remainder().iter().enumerate() {
        tail |= u64::from(byte) << (8 * i);
    }
    (h.rotate_left(23) ^ tail).wrapping_mul(K)
}

impl Memo {
    const fn new() -> Self {
        Self {
            slots: Vec::new(),
            arena: Vec::new(),
            entries: 0,
        }
    }

    /// The slot holding `segment`, or the empty slot where it would go.
    #[inline]
    fn position(&self, hash: u64, segment: &[u8], maxtrim: u8) -> usize {
        let mask = self.slots.len() - 1;
        let tag = (hash >> 32) as u32;
        let mut index = (hash as usize) & mask;
        loop {
            let slot = &self.slots[index];
            if slot.len == 0 {
                return index;
            }
            if slot.tag == tag
                && usize::from(slot.len) == segment.len()
                && slot.maxtrim == maxtrim
                && &self.arena[slot.start as usize..slot.start as usize + segment.len()] == segment
            {
                return index;
            }
            index = (index + 1) & mask;
        }
    }

    fn clear(&mut self) {
        self.slots.clear();
        self.slots.resize(X_FIRST_SLOTS, Slot::default());
        self.arena.clear();
        self.entries = 0;
    }

    /// Doubles the table and places the stored entries again.
    fn grow(&mut self) {
        let doubled = vec![Slot::default(); self.slots.len() * 2];
        let old = std::mem::replace(&mut self.slots, doubled);
        let mask = self.slots.len() - 1;
        for slot in old.into_iter().filter(|slot| slot.len != 0) {
            let start = slot.start as usize;
            let segment = &self.arena[start..start + usize::from(slot.len)];
            let mut index = (hash_of(segment, slot.maxtrim) as usize) & mask;
            while self.slots[index].len != 0 {
                index = (index + 1) & mask;
            }
            self.slots[index] = slot;
        }
    }

    /// The stored pair of `segment` with `maxtrim` (`hash = hash_of(segment, maxtrim)`), if any.
    fn stored(&mut self, hash: u64, segment: &[u8], maxtrim: u8) -> Option<(u8, u8)> {
        if self.slots.is_empty() {
            return None;
        }
        let slot = self.slots[self.position(hash, segment, maxtrim)];
        (slot.len != 0).then_some((slot.lend, slot.rend))
    }

    /// Stores the pair of `segment` with `maxtrim` unless it is stored already (another thread
    /// may have stored the same pair since `stored` was asked).
    fn store(&mut self, hash: u64, segment: &[u8], maxtrim: u8, (lend, rend): (u8, u8)) {
        if self.slots.is_empty() {
            self.clear();
        }
        let mut index = self.position(hash, segment, maxtrim);
        if self.slots[index].len != 0 {
            return;
        }
        if self.entries >= X_MAX_ENTRIES || self.arena.len() + segment.len() > X_MAX_ARENA {
            self.clear();
            index = self.position(hash, segment, maxtrim);
        } else if 2 * (self.entries + 1) > self.slots.len() {
            self.grow();
            index = self.position(hash, segment, maxtrim);
        }
        let start = self.arena.len() as u32;
        self.arena.extend_from_slice(segment);
        self.slots[index] = Slot {
            tag: (hash >> 32) as u32,
            start,
            len: segment.len() as u8,
            maxtrim,
            lend,
            rend,
        };
        self.entries += 1;
    }
}

/// The table and its lock, on a cache line of its own.
#[repr(align(64))]
struct Shared(Mutex<Memo>);

/// The segments trimmed by all threads and their (lend, rend) pairs.
static X_MEMO: Shared = Shared(Mutex::new(Memo::new()));

#[cfg(test)]
thread_local! {
    /// Memo hits of this thread so far (tests only).
    static X_HITS: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// The table, locked (a panic while it was held leaves no partial entry: an entry is written
/// whole after its letters).
fn x_table() -> std::sync::MutexGuard<'static, Memo> {
    X_MEMO
        .0
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26): the pair `trim()` returns for `segment` (offsets of the
/// trimmed segment within `segment`), from the memo when `segment` and `maxtrim` equal a stored
/// entry. `trim` must be the pure trim of exactly these letters with this `maxtrim`; it runs
/// without the memo's lock. Segments longer than `X_MAX_SEGMENT`, an empty segment and a
/// `maxtrim` above 255 are trimmed without the memo.
pub(crate) fn x_trimmed(
    segment: &[u8],
    maxtrim: usize,
    trim: impl FnOnce() -> (usize, usize),
) -> (usize, usize) {
    x_trimmed_in(x_table, segment, maxtrim, trim)
}

/// `x_trimmed` with the table given by `table` (the table, locked).
fn x_trimmed_in<T: std::ops::DerefMut<Target = Memo>>(
    table: impl Fn() -> T,
    segment: &[u8],
    maxtrim: usize,
    trim: impl FnOnce() -> (usize, usize),
) -> (usize, usize) {
    let (Ok(maxtrim), true) = (
        u8::try_from(maxtrim),
        (1..=X_MAX_SEGMENT).contains(&segment.len()),
    ) else {
        return trim();
    };
    let hash = hash_of(segment, maxtrim);
    // The guard goes at the end of the statement: `trim` runs without the lock.
    let stored = table().stored(hash, segment, maxtrim);
    if let Some((lend, rend)) = stored {
        #[cfg(test)]
        X_HITS.with(|hits| hits.set(hits.get() + 1));
        return (usize::from(lend), usize::from(rend));
    }
    let fresh = trim();
    if let (Ok(lend), Ok(rend)) = (u8::try_from(fresh.0), u8::try_from(fresh.1)) {
        table().store(hash, segment, maxtrim, (lend, rend));
    }
    fresh
}

/// Memo hits of this thread so far (tests only; a hit of the memo shared by all threads counts
/// for the thread that asked).
#[cfg(test)]
pub(crate) fn x_hits() -> u64 {
    X_HITS.with(std::cell::Cell::get)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::{Cell, RefCell};

    /// `x_trimmed` on the one table `memo`.
    impl Memo {
        fn trimmed(
            &mut self,
            segment: &[u8],
            maxtrim: usize,
            trim: impl FnOnce() -> (usize, usize),
        ) -> (usize, usize) {
            let memo = RefCell::new(std::mem::replace(self, Memo::new()));
            let pair = x_trimmed_in(|| memo.borrow_mut(), segment, maxtrim, trim);
            *self = memo.into_inner();
            pair
        }
    }

    #[test]
    fn memo_returns_the_stored_pair_for_equal_bytes_and_maxtrim_only() {
        let mut memo = Memo::new();
        let calls = Cell::new(0usize);
        let trim = |pair: (usize, usize)| {
            let calls = &calls;
            move || {
                calls.set(calls.get() + 1);
                pair
            }
        };
        let a = [1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
        let b = [1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];
        assert_eq!(memo.trimmed(&a, 50, trim((1, 9))), (1, 9));
        // A hit returns the stored pair without calling the trim.
        assert_eq!(memo.trimmed(&a, 50, trim((0, 0))), (1, 9));
        assert_eq!(calls.get(), 1);
        // One letter differs (in the tail word): evaluated.
        assert_eq!(memo.trimmed(&b, 50, trim((2, 8))), (2, 8));
        assert_eq!(calls.get(), 2);
        // Same letters, other maxtrim: evaluated.
        assert_eq!(memo.trimmed(&a, 49, trim((3, 7))), (3, 7));
        assert_eq!(calls.get(), 3);
        // A prefix is another segment.
        assert_eq!(memo.trimmed(&a[..10], 50, trim((0, 9))), (0, 9));
        assert_eq!(calls.get(), 4);
        // Too long for the memo: evaluated every time.
        let long = [7u8; X_MAX_SEGMENT + 1];
        assert_eq!(memo.trimmed(&long, 50, trim((5, 6))), (5, 6));
        assert_eq!(memo.trimmed(&long, 50, trim((5, 6))), (5, 6));
        assert_eq!(calls.get(), 6);
    }

    #[test]
    fn trim_runs_without_the_table_and_a_stored_key_is_not_stored_twice() {
        let memo = RefCell::new(Memo::new());
        let a = [3u8, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
        // The trim takes the table itself (a held table would make this borrow panic) and stores
        // the same key first, as another thread that missed it at the same moment would.
        let pair = x_trimmed_in(
            || memo.borrow_mut(),
            &a,
            50,
            || {
                memo.borrow_mut().store(hash_of(&a, 50), &a, 50, (2, 7));
                (2, 7)
            },
        );
        assert_eq!(pair, (2, 7));
        assert_eq!(memo.borrow().entries, 1);
        assert_eq!(
            x_trimmed_in(|| memo.borrow_mut(), &a, 50, || unreachable!()),
            (2, 7)
        );
    }

    #[test]
    fn memo_keeps_every_entry_through_growth_and_clears_at_its_entry_bound() {
        let mut memo = Memo::new();
        // Distinct segments: the index in base 251 over three letters, plus one more.
        let segment = |n: usize| {
            vec![
                (n % 251) as u8,
                (n / 251 % 251) as u8,
                (n / 251 / 251) as u8,
                (n % 7) as u8,
            ]
        };
        let pair = |n: usize| (n % 3, n % 5 + 3);
        let count = 3 * X_FIRST_SLOTS;
        for n in 0..count {
            assert_eq!(memo.trimmed(&segment(n), 50, || pair(n)), pair(n));
        }
        assert_eq!(memo.entries, count);
        assert!(memo.slots.len() >= 2 * count);
        // All of them are still there after the table doubled.
        for n in 0..count {
            assert_eq!(memo.trimmed(&segment(n), 50, || unreachable!()), pair(n));
        }
        // Filling past the bound clears the table once; the results stay the trim's.
        for n in count..X_MAX_ENTRIES + count + 10 {
            assert_eq!(memo.trimmed(&segment(n), 50, || pair(n)), pair(n));
        }
        assert_eq!(memo.entries, count + 10);
        assert!(memo.slots.len() <= 2 * X_MAX_ENTRIES);
        let mut evaluated = 0usize;
        for n in 0..count {
            assert_eq!(
                memo.trimmed(&segment(n), 50, || {
                    evaluated += 1;
                    pair(n)
                }),
                pair(n)
            );
        }
        assert_eq!(evaluated, count, "the first entries were cleared");
    }

    #[test]
    fn memo_clears_at_its_letter_bound() {
        let mut memo = Memo::new();
        // Distinct segments of X_MAX_SEGMENT letters: the arena fills before the entry bound.
        let segment = |n: usize| {
            let mut s = vec![(n % 251) as u8; X_MAX_SEGMENT];
            s[0] = (n / 251) as u8;
            s
        };
        let fit = X_MAX_ARENA / X_MAX_SEGMENT;
        for n in 0..fit {
            memo.trimmed(&segment(n), 50, || (0, 1));
        }
        assert_eq!(memo.entries, fit);
        assert_eq!(memo.trimmed(&segment(0), 50, || unreachable!()), (0, 1));
        // One more segment does not fit in the arena: the table is cleared first.
        memo.trimmed(&segment(fit), 50, || (0, 2));
        assert_eq!(memo.entries, 1);
        assert_eq!(memo.arena.len(), X_MAX_SEGMENT);
        assert_eq!(memo.trimmed(&segment(fit), 50, || unreachable!()), (0, 2));
    }
}
