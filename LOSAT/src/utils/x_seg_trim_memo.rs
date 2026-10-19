//! EXPERIMENT (LOSAT_X_SEGTRIMMEMO / LOSAT_X_SEGTRIMMEMOSHADOW): the trim of a SEG segment
//! computed once per thread for each segment's bytes.
//!
//! NCBI trims every candidate segment that `s_SegSeq` finds, by evaluating `s_GetProb` on every
//! window of every length from the segment's length down to `max(1, length - maxtrim) + 1`:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seg.c:1988-1992,1998-2001,2007-2008,2013-2014
//! ```c
//!    if ((seq->length-maxtrim)>minlen)
//!         minlen = seq->length-maxtrim;
//!    minprob = 1.;
//!    for (len=seq->length; len>minlen; len--)
//! ...
//!       while (shift)
//!       {
//!          prob = s_GetProb(win->state, len, win->palpha);
//!          if (prob<minprob)
//! ...
//!          shift = s_ShiftWin1(win);
//!          i++;
//! ...
//!    *leftend = *leftend + lend;
//!    *rightend = *rightend - (seq->length - rend - 1);
//! ```
//!
//! The pair (lend, rend) is a function of the segment's letters, `maxtrim` and the alphabet;
//! no position, no entropy value and no other state of the caller enters it. The callers of the
//! fast path (`SegMasker::x_mask_sequence`) always use the 20-letter alphabet. TBLASTN's
//! composition-based redo masks a window of the translated subject around each alignment
//! (blast_kappa.c:1530-1548), so a low-complexity stretch of the subject is trimmed again for
//! every query whose window covers it (v4-C DECISIONS 31: 118,765 trims of 12,584 distinct
//! segments on AP027131 x AP027133); BLASTX's query frames repeat segments too.
//!
//! This module keeps, per thread, the bytes of the segments it trimmed and their pairs. A call
//! whose bytes and `maxtrim` equal a stored entry (the whole slice is compared; the hash only
//! selects the slot) returns the stored pair; any other call evaluates the trim and stores it.
//! The table grows on demand up to `X_MAX_ENTRIES` entries (or `X_MAX_ARENA` letters) and is
//! then cleared as a whole.
//!
//! No NCBI counterpart: the memo itself; it returns only values the trim computed for the same
//! input and does not change any value NCBI computes.

use std::cell::RefCell;

/// EXPERIMENT (LOSAT_X_SEGTRIMMEMO): 0 = off, 1 = memo, 2 = shadow
/// (`LOSAT_X_SEGTRIMMEMOSHADOW`: on a hit the trim is evaluated too and compared).
pub(crate) fn x_segtrimmemo_mode() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_SEGTRIMMEMOSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_SEGTRIMMEMO") {
            1
        } else {
            0
        }
    })
}

/// EXPERIMENT (LOSAT_X_SEGTRIMMEMOSHADOW): memo hits compared with the trim so far.
pub(crate) static X_SHADOW_CALLS: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_SEGTRIMMEMOSHADOW): the shadow-mode summary (printed by `main` at exit
/// in shadow mode).
// No NCBI counterpart: prints the shadow-mode counter; it does not change any value NCBI computes.
pub fn print_shadow_stats() {
    if x_segtrimmemo_mode() == 2 {
        eprintln!(
            "[X_SHADOW] SEGTRIMMEMO memo_hits_compared={} (trimmed segments identical)",
            X_SHADOW_CALLS.load(std::sync::atomic::Ordering::Relaxed)
        );
    }
}

// No NCBI counterpart: the size bounds of the memo; they decide only which trims are evaluated
// again and do not change any value NCBI computes.
/// Most entries kept by one thread; the table is cleared when it would hold more. 32,768 keeps
/// every distinct segment of the measured jobs (12,584 for TBLASTN AP027131 x AP027133, 25,507
/// for BLASTX AP027131 x AP027133; v4-C DECISIONS 32).
const X_MAX_ENTRIES: usize = 32_768;
/// Most letters kept by one thread (the stored segments' bytes); the table is cleared when it
/// would hold more.
const X_MAX_ARENA: usize = 1 << 20;
/// Slots of a new table; it doubles while it is more than half full, up to twice
/// `X_MAX_ENTRIES`.
const X_FIRST_SLOTS: usize = 1_024;
/// Longest segment stored (the fast trim's own bound, `X_TRIM_MAX` in seg.rs). Longer segments
/// are trimmed every time.
pub(crate) const X_MAX_SEGMENT: usize = 127;

/// One stored trim. `len == 0` marks an empty slot (a stored segment has at least one letter).
#[derive(Clone, Copy, Default)]
struct Slot {
    hash: u32,
    start: u32,
    len: u8,
    maxtrim: u8,
    lend: u8,
    rend: u8,
}

#[derive(Default)]
struct Memo {
    slots: Vec<Slot>,
    arena: Vec<u8>,
    entries: usize,
    #[cfg(test)]
    hits: u64,
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
            if slot.hash == tag
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
}

thread_local! {
    /// The segments this thread trimmed and their (lend, rend) pairs.
    static X_MEMO: RefCell<Memo> = RefCell::new(Memo::default());
}

/// EXPERIMENT (LOSAT_X_SEGTRIMMEMO): the pair `trim()` returns for `segment` (offsets of the
/// trimmed segment within `segment`), from the per-thread memo when `segment` and `maxtrim`
/// equal a stored entry. `trim` must be the pure trim of exactly these letters with this
/// `maxtrim`. Segments longer than `X_MAX_SEGMENT`, an empty segment and a `maxtrim` above 255
/// are trimmed without the memo. With the shadow switch a hit evaluates `trim` as well and
/// compares.
pub(crate) fn x_trimmed(
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
    X_MEMO.with(|memo| {
        let mut memo = memo.borrow_mut();
        if memo.slots.is_empty() {
            memo.clear();
        }
        let index = memo.position(hash, segment, maxtrim);
        let slot = memo.slots[index];
        if slot.len != 0 {
            let stored = (usize::from(slot.lend), usize::from(slot.rend));
            #[cfg(test)]
            {
                memo.hits += 1;
            }
            if x_segtrimmemo_mode() == 2 {
                let fresh = trim();
                assert!(
                    fresh == stored,
                    "LOSAT_X_SEGTRIMMEMOSHADOW: memoised trim {stored:?} differs from {fresh:?} on a {}-residue segment",
                    segment.len()
                );
                X_SHADOW_CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            return stored;
        }
        let fresh = trim();
        let (Ok(lend), Ok(rend)) = (u8::try_from(fresh.0), u8::try_from(fresh.1)) else {
            return fresh;
        };
        let mut index = index;
        if memo.entries >= X_MAX_ENTRIES || memo.arena.len() + segment.len() > X_MAX_ARENA {
            memo.clear();
            index = memo.position(hash, segment, maxtrim);
        } else if 2 * (memo.entries + 1) > memo.slots.len() {
            memo.grow();
            index = memo.position(hash, segment, maxtrim);
        }
        let start = memo.arena.len() as u32;
        memo.arena.extend_from_slice(segment);
        memo.slots[index] = Slot {
            hash: (hash >> 32) as u32,
            start,
            len: segment.len() as u8,
            maxtrim,
            lend,
            rend,
        };
        memo.entries += 1;
        fresh
    })
}

/// Memo hits of this thread so far (tests only).
#[cfg(test)]
pub(crate) fn x_hits() -> u64 {
    X_MEMO.with(|memo| memo.borrow().hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn memo_returns_the_stored_pair_for_equal_bytes_and_maxtrim_only() {
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
        assert_eq!(x_trimmed(&a, 50, trim((1, 9))), (1, 9));
        assert_eq!(x_trimmed(&a, 50, trim((0, 0))), (1, 9));
        assert_eq!(calls.get(), 1);
        // One letter differs (in the tail word): evaluated.
        assert_eq!(x_trimmed(&b, 50, trim((2, 8))), (2, 8));
        assert_eq!(calls.get(), 2);
        // Same letters, other maxtrim: evaluated.
        assert_eq!(x_trimmed(&a, 49, trim((3, 7))), (3, 7));
        assert_eq!(calls.get(), 3);
        // A prefix is another segment.
        assert_eq!(x_trimmed(&a[..10], 50, trim((0, 9))), (0, 9));
        assert_eq!(calls.get(), 4);
        // Too long for the memo: evaluated every time.
        let long = [7u8; X_MAX_SEGMENT + 1];
        assert_eq!(x_trimmed(&long, 50, trim((5, 6))), (5, 6));
        assert_eq!(x_trimmed(&long, 50, trim((5, 6))), (5, 6));
        assert_eq!(calls.get(), 6);
    }

    #[test]
    fn memo_keeps_every_entry_through_growth_and_clears_at_its_bound() {
        // Distinct segments: the index in base 251 over three letters, plus a length byte.
        let segment = |n: usize| {
            vec![
                (n % 251) as u8,
                (n / 251 % 251) as u8,
                (n / 251 / 251) as u8,
                (n % 7) as u8,
            ]
        };
        let pair = |n: usize| (n % 3, n % 5 + 3);
        let calls = Cell::new(0usize);
        let count = 3 * X_FIRST_SLOTS;
        for n in 0..count {
            let _ = x_trimmed(&segment(n), 50, || {
                calls.set(calls.get() + 1);
                pair(n)
            });
        }
        assert_eq!(calls.get(), count);
        // All of them are still there after the table doubled.
        for n in 0..count {
            assert_eq!(x_trimmed(&segment(n), 50, || unreachable!()), pair(n));
        }
        // Filling past the bound clears the table; the results stay the trim's.
        for n in count..X_MAX_ENTRIES + count + 10 {
            assert_eq!(x_trimmed(&segment(n), 50, || pair(n)), pair(n));
        }
        let mut evaluated = 0usize;
        for n in 0..count {
            assert_eq!(
                x_trimmed(&segment(n), 50, || {
                    evaluated += 1;
                    pair(n)
                }),
                pair(n)
            );
        }
        assert_eq!(evaluated, count, "the first entries were cleared");
    }
}
