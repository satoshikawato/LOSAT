//! No NCBI counterpart: `x_mask_range` (the symdust loop of `CSymDustMasker::operator()`) with the
//! scalars of the triplet window in local variables and the triplet ring written twice
//! (`LOSAT_X_DUSTREG`; `=0` runs `XTripletWindow`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:52-163,229-282
//!
//! What changes. `XTripletWindow` holds the window's scalars (`head`, `len`, `start`, `stop`,
//! `l`, `r_w`, `r_v`, `num_diff`) in the same struct as its arrays, and the struct is handed by
//! reference to `find_perfect` and `shift_high`, which are not inlined. Its address is therefore
//! taken inside the loop over the letters, so every scalar is loaded from and stored to the stack
//! on every letter, and each of them becomes a store-forwarding chain (line FR s1, H1 §1.2: 98.9
//! instructions and 22.1 cycles per letter on the Sakai query). Here the scalars are a `Copy`
//! value (`XHot`) that never has its address taken: `find_perfect` receives it by value, and the
//! rare `shift_high` works on a copy that is written back. Only the three arrays (`XArrays`) stay
//! in memory, as they must (they are indexed by the data). `stop` is not kept: it is always
//! `start + len` (`XHot::stop`), one register fewer in the loop.
//!
//! The ring of triplets is also written twice (`ring[head]` and `ring[head + 64]`), so the deque
//! index `i` (0 = newest) is `ring[head + i]` for every `i < 64` without the `& 63`, and the part
//! of the window `find_perfect` reads is one contiguous slice. `find_perfect`'s loop keeps only
//! the counts, the score, the triplet and its threshold: `cnt > 0` and the score test are one
//! comparison (`x_gate_limit`, a conditional move, against the thresholds divided by ten),
//! `count < thresholds.len()` always holds (`count < len <=
//! window - 2`), and the body of the condition, which holds on few iterations, is a function of
//! its own (`x_merge_at`) that derives `pos` and `count` from the iteration's index and walks the
//! perfect list as one slice (`VecDeque::make_contiguous`), copying the intervals it passes in one
//! piece.
//!
//! A full window stays full: each letter pops one triplet and pushes one. While it is full and
//! holds two or more distinct triplets, `x_run_full` takes the letters with only the scalars that
//! change (`head`, the suffix length, `r_w`, `r_v`, `num_diff`): `start` rises with `pos`, `L`
//! follows from the suffix length, and `needs_processing` is one comparison with a table of
//! thresholds that has `u32::MAX` past its end (`x_tenth_thresholds`). The window's first letters,
//! the rare `shift_high` branch and the run after a uniform window go through `XHot`, as before.
//!
//! Why the result is the same. Every method below is the `XTripletWindow` method of the same name
//! with the same statements in the same order on the same integers; only where the values live
//! changes. The ring holds the same triplet at `head + i` as `XTripletWindow` holds at
//! `(head + i) & 63` for every `i` below `len` (each push writes both copies, `head` stays in
//! `0..64` and `len` is at most `window - 2 <= 62`, so `head + i < 128`). The combined condition
//! is the C condition: the score exceeds a threshold divided by ten (rounded down) exactly when
//! ten times the score exceeds the threshold (`x_tenth_thresholds`); where `cnt` is 0 the
//! threshold is `u32::MAX`, which the score never exceeds (at most 1,891), and otherwise it is the
//! C threshold divided by ten; the third test is always true. `pos` of iteration `i` is the start value lowered `i` times with
//! `saturating_sub(1)`, which is `saturating_sub(i)`; `count` rises by one per iteration. The loop in `x_mask_range_regs` is `DustMasker::x_mask_range` statement by
//! statement. The unit tests compare it with `mask_range_reference` and with `x_mask_range` on
//! random sequences of every kind the masker distinguishes.

use super::{DustMasker, MaskedInterval, NcbiLfgRandom, PerfectInterval, X_BASE_CODE, X_RING};
use std::collections::VecDeque;

/// No NCBI counterpart: the scalars of `XTripletWindow` (symdust.cpp:40-48 `triplets`), kept by
/// value so the loop holds them in registers.
#[derive(Clone, Copy)]
struct XHot {
    head: usize,
    len: usize,
    start: usize,
    l: usize,
    r_w: u32,
    r_v: u32,
    num_diff: u32,
}

/// No NCBI counterpart: the arrays of `XTripletWindow`, with the ring written twice.
struct XArrays {
    ring: [u8; 2 * X_RING],
    c_w: [u8; 64],
    c_v: [u8; 64],
}

/// No NCBI counterpart: the constants of one window.
#[derive(Clone, Copy)]
struct XConst<'a> {
    max_size: usize,
    low_k: u8,
    thresholds: &'a [u32],
}

/// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:275-277
#[inline(always)]
fn add_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
    let idx = (triplet & 63) as usize;
    *sum += counts[idx] as u32;
    counts[idx] += 1;
}

/// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:287-289
#[inline(always)]
fn rem_triplet(sum: &mut u32, counts: &mut [u8; 64], triplet: u8) {
    let idx = (triplet & 63) as usize;
    counts[idx] -= 1;
    *sum -= counts[idx] as u32;
}

impl XHot {
    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:43-48
    fn new() -> Self {
        Self {
            head: 0,
            len: 0,
            start: 0,
            l: 0,
            r_w: 0,
            r_v: 0,
            num_diff: 0,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:40-48
    /// `stop_` of the C class, which is always `start_` plus the deque's length: both start at 0,
    /// a push adds one to the length and to `stop_`, a pop takes one from the length and adds one
    /// to `start_` (`shift_window`, `shift_high`).
    #[inline(always)]
    fn stop(&self) -> usize {
        self.start + self.len
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:54-55
    /// Returns the oldest triplet and removes it (`XTripletWindow::pop_back`).
    #[inline(always)]
    fn pop_back(&mut self, a: &XArrays) -> u8 {
        self.len -= 1;
        a.ring[self.head + self.len]
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:60
    /// Makes `triplet` the newest element (`XTripletWindow::push_front`), in both copies.
    #[inline(always)]
    fn push_front(&mut self, a: &mut XArrays, triplet: u8) {
        self.head = (self.head + X_RING - 1) & (X_RING - 1);
        a.ring[self.head] = triplet;
        a.ring[self.head + X_RING] = triplet;
        self.len += 1;
    }

    /// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:251-256
    /// `XTripletWindow::needs_processing`.
    #[inline(always)]
    fn needs_processing(&self, k: XConst) -> bool {
        let count = self.stop() - self.l;
        let last = k.thresholds.len() - 1;
        (count < self.len)
            & (count < k.thresholds.len())
            & (10 * self.r_w > k.thresholds[count.min(last)])
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:76-114
    /// `XTripletWindow::shift_window`.
    #[inline(always)]
    fn shift_window(
        &mut self,
        a: &mut XArrays,
        k: XConst,
        triplet: u8,
        perfect_list: &mut VecDeque<PerfectInterval>,
    ) -> bool {
        if self.len >= k.max_size {
            if self.num_diff <= 1 {
                let mut hot = *self;
                let kept = shift_high(&mut hot, a, triplet, perfect_list);
                *self = hot;
                return kept;
            }

            let old = (self.pop_back(a) & 63) as usize;
            a.c_w[old] -= 1;
            self.r_w -= a.c_w[old] as u32;
            self.num_diff -= (a.c_w[old] == 0) as u32;

            // When `L` is `start_`, NCBI moves `L` up by one and removes the triplet from the
            // suffix's counts (`rem_triplet_info`).
            let at_start = (self.l == self.start) as usize;
            self.l += at_start;
            let suffix_count = a.c_v[old] - at_start as u8;
            a.c_v[old] = suffix_count;
            self.r_v -= suffix_count as u32 * at_start as u32;

            self.start += 1;
        }

        self.push_front(a, triplet);
        self.num_diff += (a.c_w[(triplet & 63) as usize] == 0) as u32;
        add_triplet(&mut self.r_w, &mut a.c_w, triplet);
        add_triplet(&mut self.r_v, &mut a.c_v, triplet);

        if a.c_v[(triplet & 63) as usize] > k.low_k {
            let mut off = self.len - (self.l - self.start) - 1;
            loop {
                let t = a.ring[self.head + off];
                rem_triplet(&mut self.r_v, &mut a.c_v, t);
                self.l += 1;
                if t == triplet {
                    break;
                }
                if off == 0 {
                    break;
                }
                off -= 1;
            }
        }

        // NCBI's `stop_` rises by one here; `stop()` already counts the push.

        if self.len >= k.max_size && self.num_diff <= 1 {
            perfect_list.clear();
            perfect_list.push_front(PerfectInterval::new(self.start, self.stop() + 1, 0, 0));
            return false;
        }

        true
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:52-69
/// `XTripletWindow::shift_high`.
#[inline(never)]
fn shift_high(
    h: &mut XHot,
    a: &mut XArrays,
    triplet: u8,
    perfect_list: &mut VecDeque<PerfectInterval>,
) -> bool {
    let old_triplet = h.pop_back(a);
    rem_triplet(&mut h.r_w, &mut a.c_w, old_triplet);
    if a.c_w[(old_triplet & 63) as usize] == 0 {
        h.num_diff -= 1;
    }
    h.start += 1;

    h.push_front(a, triplet);
    if a.c_w[(triplet & 63) as usize] == 0 {
        h.num_diff += 1;
    }
    add_triplet(&mut h.r_w, &mut a.c_w, triplet);
    // NCBI's `stop_` rises by one here; `stop()` already counts the push.

    if h.num_diff <= 1 {
        perfect_list.push_front(PerfectInterval::new(h.start, h.stop() + 1, 0, 0));
        return false;
    }

    true
}

/// No NCBI counterpart: the state of the merge `find_perfect` makes when the score test of an
/// iteration holds (the body of that `if` in `XTripletWindow::find_perfect`).
struct XMerge {
    max_perfect_score: u32,
    max_len: usize,
    read: usize,
    inserted: bool,
    old_len: usize,
    /// `pos` of the first iteration (`L - 1`, or 0); iteration `i` has `pos0.saturating_sub(i)`.
    pos0: usize,
    /// `count` of the first iteration (the suffix length); iteration `i` has `count0 + i`.
    count0: usize,
    /// One past `stop_`, the end of every interval this call adds.
    end: usize,
}

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:123-141,161-163
/// `XTripletWindow::find_perfect`, reading the window as one slice of the doubled ring. The loop
/// keeps only the counts, the score, the triplet and its threshold; the body of the score test,
/// which holds on few iterations, is `x_merge_at` (not inlined, so its state does not take the
/// loop's registers). `tenths` holds the thresholds divided by ten (`x_tenth_thresholds`), so the
/// score test compares the score itself. `count < tenths.len()` always holds: `count` runs from the
/// suffix length to `len - 1`, and `len <= window - 2 = tenths.len()`, so the thresholds of the
/// iterations are the slice `tenths[suffix_len..len]`.
#[inline(never)]
fn find_perfect(
    h: XHot,
    a: &XArrays,
    tenths: &[u32],
    merged: &mut Vec<PerfectInterval>,
    perfect_list: &mut VecDeque<PerfectInterval>,
) {
    let stop = h.stop();
    let suffix_len = stop - h.l;
    if suffix_len >= h.len {
        return;
    }
    let mut counts = a.c_v;
    let mut score = h.r_v;
    merged.clear();
    let mut m = XMerge {
        max_perfect_score: 0,
        max_len: 0,
        read: 0,
        inserted: false,
        old_len: perfect_list.len(),
        pos0: h.l.saturating_sub(1),
        count0: suffix_len,
        end: stop + 1,
    };
    // The same elements in the same order, in one slice (the list is almost always one already).
    let list: &[PerfectInterval] = perfect_list.make_contiguous();
    let triplets = &a.ring[h.head + suffix_len..h.head + h.len];
    let limits = &tenths[suffix_len..h.len];
    for (i, (&triplet, &limit)) in triplets.iter().zip(limits).enumerate() {
        let idx = (triplet & 63) as usize;
        let cnt = counts[idx];
        add_triplet(&mut score, &mut counts, triplet);
        if score > x_gate_limit(cnt, limit) {
            x_merge_at(&mut m, i, score, merged, list);
        }
    }
    if m.inserted {
        merged.extend_from_slice(&list[m.read..]);
        perfect_list.clear();
        perfect_list.extend(merged.drain(..));
    }
}

/// No NCBI counterpart: the threshold of one `find_perfect` iteration whose triplet had the count
/// `cnt` in the suffix before it was added: `limit` (a threshold divided by ten) where `cnt > 0`,
/// otherwise `u32::MAX`, which a window's score never exceeds (the score of 62 triplets is at most
/// 1,891). So `score > x_gate_limit(cnt, limit)` is NCBI's two tests (`cnt > 0` and the score above
/// the threshold) as one comparison. Whether `cnt` is 0 follows the data (the first time a
/// triplet appears in the scanned part of the window), and a branch on it is mispredicted on
/// about every fifth iteration in cDNA queries (line FR s1, H1 §1.2: 0.19 to 0.22 per
/// iteration). LLVM's x86 back end splits two tests joined by `&` into two branches again (ws,
/// dx2 listing of the Ad case), as it did for `select_unpredictable` in the Smith-Waterman scan
/// (`cmov_max`), so on x86-64 the choice is a conditional move written out. Other targets take
/// the plain `if`, which gives the same value.
#[inline(always)]
fn x_gate_limit(cnt: u8, limit: u32) -> u32 {
    #[cfg(target_arch = "x86_64")]
    {
        let mut gated = limit;
        // SAFETY: two instructions on the registers given (`test` of `cnt`, then `cmovz`:
        // `gated` takes `u32::MAX` where `cnt` is 0); no memory or stack is touched, and the
        // flags that `test` sets are clobbered as `asm!` assumes by default.
        unsafe {
            core::arch::asm!(
                "test {cnt:e}, {cnt:e}",
                "cmovz {gated:e}, {max:e}",
                cnt = in(reg) cnt as u32,
                gated = inout(reg) gated,
                max = in(reg) u32::MAX,
                options(pure, nomem, nostack),
            );
        }
        gated
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        if cnt > 0 {
            limit
        } else {
            u32::MAX
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:125-141
/// The body of the score test of `XTripletWindow::find_perfect` for iteration `i` (its `pos` and
/// `count` from `i`), in the same order. `list` is the perfect list as one slice; the intervals the
/// walk passes are copied to `merged` in one piece after the walk, in the order the C code moves
/// them one by one (it reads them only, so copying them later gives the same list).
#[inline(never)]
fn x_merge_at(
    m: &mut XMerge,
    i: usize,
    score: u32,
    merged: &mut Vec<PerfectInterval>,
    list: &[PerfectInterval],
) {
    let pos = m.pos0.saturating_sub(i);
    let count = m.count0 + i;
    let walk_start = m.read;
    while let Some(p) = list.get(m.read) {
        if pos > p.start {
            break;
        }
        if m.max_perfect_score == 0
            || m.max_len * p.score as usize > m.max_perfect_score as usize * p.len
        {
            m.max_perfect_score = p.score;
            m.max_len = p.len;
        }
        m.read += 1;
    }
    if m.inserted {
        merged.extend_from_slice(&list[walk_start..m.read]);
    }
    if m.max_perfect_score == 0
        || score as usize * m.max_len >= m.max_perfect_score as usize * count
    {
        m.max_perfect_score = score;
        m.max_len = count;
        if !m.inserted {
            m.inserted = true;
            merged.reserve(m.old_len + 8);
            merged.extend_from_slice(&list[..m.read]);
        }
        merged.push(PerfectInterval::new(pos, m.end, m.max_perfect_score, count));
    }
}

/// NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:75-84
/// The 2-bit code of `base` (`convert_iupac_to_ncbi2na`), drawing from `rng` where the table says 4.
#[inline(always)]
fn x_code(base: u8, rng: &mut NcbiLfgRandom) -> u8 {
    let c = X_BASE_CODE[base as usize];
    if c < 4 {
        c
    } else {
        (rng.next() & 0x3) as u8
    }
}

/// No NCBI counterpart: the thresholds divided by ten (rounded down), with one more entry,
/// `u32::MAX`, at the index of a suffix as long as the window. For a score `x` and a threshold
/// `t`, ten times `x` exceeds `t` exactly when `x` exceeds `t / 10` rounded down (if `x` is above
/// it, `10 * x` is at least ten times it plus ten, which is above `t`; if not, `10 * x` is at most
/// ten times it, at most `t`), so the comparisons need no multiplication. In a full window the
/// suffix length `sl` (the distance from NCBI's `L` up to its `stop_`) is at most the window's
/// length `max_size`, and `needs_processing`'s three tests (the suffix shorter than the window,
/// shorter than the table, ten times the window's score above its threshold) are one comparison
/// with this table: the C table has `max_size` entries, so for `sl < max_size` the entry is the C
/// threshold divided by ten and both length tests hold, and for `sl == max_size` the entry is
/// `u32::MAX`, which a score never exceeds (at most 1,891).
fn x_tenth_thresholds(k: XConst) -> [u32; X_RING] {
    debug_assert_eq!(k.thresholds.len(), k.max_size);
    let mut thr = [u32::MAX; X_RING];
    for (t, &c) in thr.iter_mut().zip(k.thresholds) {
        *t = c / 10;
    }
    thr
}

/// No NCBI counterpart: the loop's position, last triplet and window, handed to and from the
/// runs below.
struct XRun {
    pos: usize,
    current_triplet: u8,
    h: XHot,
}

/// No NCBI counterpart: why `x_run_full` stopped.
enum XFullExit {
    /// No letters left (`pos` is the end).
    End,
    /// The letter's `shift_window` returned false (the window became uniform): the run of
    /// `x_run_uniform` comes next, as in the loop.
    Uniform,
}

impl DustMasker {
    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:229-262, c++/src/algo/dustmask/symdust.cpp:76-114
    /// The letters of the loop of `x_mask_range` while the window is full (its length `max_size`)
    /// and holds two or more distinct triplets before each letter: the loop's save, triplet and
    /// `pos` steps, then `shift_window`'s branch for a full window (`XHot::shift_window`, the same
    /// statements on the same values), then `needs_processing` and `find_perfect`. A full window
    /// stays full (each letter pops one triplet and pushes one), so the run keeps no length, and
    /// the scalars it carries are the ones that change: `head`, the suffix length `sl` (from `L` up to
    /// `stop`),
    /// `r_w`, `r_v` and `num_diff`. `start` rises by one with `pos` on every letter, so it is
    /// `pos - delta` for the run's `delta`; `L` is `start + max_size - sl`. With `stop = start +
    /// max_size` (`XHot::stop`), `shift_window`'s updates of `L` and `stop` (`L` up by one where it
    /// equals `start`, then `start` and `stop` up by one) leave `sl` the same where `sl == max_size`
    /// and raise it by one otherwise; each step of the suffix trim raises `L`, so it lowers `sl`, and
    /// the trim reads the ring at `sl - 1` (`XHot::shift_window`'s `off`).
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn x_run_full(
        &self,
        seq: &[u8],
        run: &mut XRun,
        a: &mut XArrays,
        k: XConst,
        thr: &[u32; X_RING],
        rng: &mut NcbiLfgRandom,
        result: &mut Vec<MaskedInterval>,
        current_start: usize,
        merged: &mut Vec<PerfectInterval>,
        perfect_list: &mut VecDeque<PerfectInterval>,
    ) -> XFullExit {
        let end = seq.len();
        let max = k.max_size;
        let mut pos = run.pos;
        let mut cur = run.current_triplet;
        let delta = pos - run.h.start;
        let mut head = run.h.head;
        let mut sl = run.h.stop() - run.h.l;
        let mut r_w = run.h.r_w;
        let mut r_v = run.h.r_v;
        let mut num_diff = run.h.num_diff;
        let exit = loop {
            if pos >= end {
                break XFullExit::End;
            }
            // A letter that leaves the window with at most one distinct triplet ends the run
            // (`Uniform`), and the run starts with two or more.
            debug_assert!(num_diff > 1);
            if !perfect_list.is_empty() {
                self.save_masked_regions(result, pos - delta, current_start, perfect_list);
            }
            let triplet = ((cur << 2) & 0x3F) + x_code(seq[pos], rng);
            cur = triplet;
            pos += 1;

            // `XHot::shift_window`, the branch of a full window with two or more distinct triplets.
            let old = (a.ring[(head + max - 1) & (2 * X_RING - 1)] & 63) as usize;
            a.c_w[old] -= 1;
            r_w -= a.c_w[old] as u32;
            num_diff -= (a.c_w[old] == 0) as u32;
            // `L` is `start` (the suffix is the whole window): NCBI moves `L` up with `start` and
            // takes the triplet out of the suffix's counts, and the suffix keeps its length.
            // Otherwise the suffix gains the pushed triplet. Which case holds changes on about 3 %
            // of the letters (line FR s1, H1 §1.1: 2.44 M of 91.2 M steps in cDNA), so this is a
            // branch, not the selects of `XHot::shift_window`.
            if sl == max {
                a.c_v[old] -= 1;
                r_v -= a.c_v[old] as u32;
            } else {
                sl += 1;
            }

            head = (head + X_RING - 1) & (X_RING - 1);
            a.ring[head] = triplet;
            a.ring[head + X_RING] = triplet;
            let idx = (triplet & 63) as usize;
            num_diff += (a.c_w[idx] == 0) as u32;
            add_triplet(&mut r_w, &mut a.c_w, triplet);
            add_triplet(&mut r_v, &mut a.c_v, triplet);

            if a.c_v[idx] > k.low_k {
                loop {
                    let t = a.ring[(head + sl - 1) & (2 * X_RING - 1)];
                    rem_triplet(&mut r_v, &mut a.c_v, t);
                    sl -= 1;
                    if t == triplet || sl == 0 {
                        break;
                    }
                }
            }

            if num_diff <= 1 {
                let start = pos - delta;
                perfect_list.clear();
                perfect_list.push_front(PerfectInterval::new(start, start + max + 1, 0, 0));
                break XFullExit::Uniform;
            }

            // `XHot::needs_processing` (`x_tenth_thresholds`).
            if r_w > thr[sl & (X_RING - 1)] {
                let start = pos - delta;
                let h = XHot {
                    head,
                    len: max,
                    start,
                    l: start + max - sl,
                    r_w,
                    r_v,
                    num_diff,
                };
                find_perfect(h, a, &thr[..max], merged, perfect_list);
            }
        };
        let start = pos - delta;
        run.pos = pos;
        run.current_triplet = cur;
        run.h = XHot {
            head,
            len: max,
            start,
            l: start + max - sl,
            r_w,
            r_v,
            num_diff,
        };
        exit
    }

    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:247-262
    /// The inner loop of `x_mask_range` after a `shift_window` that returned false, statement by
    /// statement. Returns true where a `shift_window` returned true (the window is done).
    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    fn x_run_uniform(
        &self,
        seq: &[u8],
        run: &mut XRun,
        a: &mut XArrays,
        k: XConst,
        rng: &mut NcbiLfgRandom,
        result: &mut Vec<MaskedInterval>,
        current_start: usize,
        perfect_list: &mut VecDeque<PerfectInterval>,
    ) -> bool {
        let end = seq.len();
        while run.pos < end {
            if !perfect_list.is_empty() {
                self.save_masked_regions(result, run.h.start, current_start, perfect_list);
            }

            let new_triplet = ((run.current_triplet << 2) & 0x3F) + x_code(seq[run.pos], rng);
            run.current_triplet = new_triplet;

            if run.h.shift_window(a, k, new_triplet, perfect_list) {
                return true;
            }
            run.pos += 1;
        }
        false
    }
}

impl DustMasker {
    /// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:229-282
    /// `x_mask_range` statement by statement, with the window in `XHot` + `XArrays` (module doc).
    #[inline(always)]
    pub(super) fn x_mask_range_regs(
        &self,
        seq: &[u8],
        start: usize,
        stop: usize,
        rng: &mut NcbiLfgRandom,
    ) -> Vec<MaskedInterval> {
        let seq = &seq[..=stop];
        let end = seq.len();
        let mut result = Vec::new();
        let mut current_start = start;
        let k = XConst {
            max_size: self.window - 2,
            low_k: self.low_k,
            thresholds: &self.thresholds,
        };
        let thr = x_tenth_thresholds(k);
        let mut merged: Vec<PerfectInterval> = Vec::new();

        while stop > current_start + 2 {
            let mut perfect_list: VecDeque<PerfectInterval> = VecDeque::new();
            let mut h = XHot::new();
            let mut a = XArrays {
                ring: [0; 2 * X_RING],
                c_w: [0; 64],
                c_v: [0; 64],
            };

            let mut current_triplet =
                (x_code(seq[current_start], rng) << 2) + x_code(seq[current_start + 1], rng);
            let mut pos = current_start + 2;
            let mut done = false;

            while !done && pos < end {
                // Once the window is full it stays full: `x_run_full` takes the letters while the
                // window holds two or more distinct triplets before the letter (module doc).
                if h.len == k.max_size && h.num_diff > 1 {
                    let mut run = XRun {
                        pos,
                        current_triplet,
                        h,
                    };
                    let exit = self.x_run_full(
                        seq,
                        &mut run,
                        &mut a,
                        k,
                        &thr,
                        rng,
                        &mut result,
                        current_start,
                        &mut merged,
                        &mut perfect_list,
                    );
                    (pos, current_triplet, h) = (run.pos, run.current_triplet, run.h);
                    match exit {
                        XFullExit::End => break,
                        XFullExit::Uniform => {
                            let mut run = XRun {
                                pos,
                                current_triplet,
                                h,
                            };
                            done = self.x_run_uniform(
                                seq,
                                &mut run,
                                &mut a,
                                k,
                                rng,
                                &mut result,
                                current_start,
                                &mut perfect_list,
                            );
                            (pos, current_triplet, h) = (run.pos, run.current_triplet, run.h);
                            continue;
                        }
                    }
                }

                if !perfect_list.is_empty() {
                    self.save_masked_regions(
                        &mut result,
                        h.start,
                        current_start,
                        &mut perfect_list,
                    );
                }

                let new_triplet = ((current_triplet << 2) & 0x3F) + x_code(seq[pos], rng);
                current_triplet = new_triplet;
                pos += 1;

                if h.shift_window(&mut a, k, new_triplet, &mut perfect_list) {
                    if h.needs_processing(k) {
                        find_perfect(h, &a, &thr[..k.max_size], &mut merged, &mut perfect_list);
                    }
                } else {
                    let mut run = XRun {
                        pos,
                        current_triplet,
                        h,
                    };
                    done = self.x_run_uniform(
                        seq,
                        &mut run,
                        &mut a,
                        k,
                        rng,
                        &mut result,
                        current_start,
                        &mut perfect_list,
                    );
                    (pos, current_triplet, h) = (run.pos, run.current_triplet, run.h);
                }
            }

            let mut wstart = h.start;
            while !perfect_list.is_empty() {
                self.save_masked_regions(&mut result, wstart, current_start, &mut perfect_list);
                wstart += 1;
            }

            if h.start > 0 {
                current_start += h.start;
            } else {
                break;
            }
        }

        result
    }
}
