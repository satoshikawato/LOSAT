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
//! condition without short-circuit, `count < thresholds.len()` always holds (`count < len <=
//! window - 2`), and the body of the condition, which holds on few iterations, is a function of
//! its own (`x_merge_at`) that derives `pos` and `count` from the iteration's index.
//!
//! Why the result is the same. Every method below is the `XTripletWindow` method of the same name
//! with the same statements in the same order on the same integers; only where the values live
//! changes. The ring holds the same triplet at `head + i` as `XTripletWindow` holds at
//! `(head + i) & 63` for every `i` below `len` (each push writes both copies, `head` stays in
//! `0..64` and `len` is at most `window - 2 <= 62`, so `head + i < 128`). The combined condition
//! is the C condition: `&` of two `bool`s is true exactly when both are, and the third test is
//! always true. `pos` of iteration `i` is the start value lowered `i` times with
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
/// loop's registers). `count < thresholds.len()` always holds: `count` runs from the suffix length
/// to `len - 1`, and `len <= window - 2 = thresholds.len()`, so the thresholds of the iterations
/// are the slice `thresholds[suffix_len..len]`.
#[inline(never)]
fn find_perfect(
    h: XHot,
    a: &XArrays,
    thresholds: &[u32],
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
    let triplets = &a.ring[h.head + suffix_len..h.head + h.len];
    let limits = &thresholds[suffix_len..h.len];
    for (i, (&triplet, &limit)) in triplets.iter().zip(limits).enumerate() {
        let idx = (triplet & 63) as usize;
        let cnt = counts[idx];
        add_triplet(&mut score, &mut counts, triplet);
        if (cnt > 0) & (score * 10 > limit) {
            x_merge_at(&mut m, i, score, merged, perfect_list);
        }
    }
    if m.inserted {
        merged.extend(perfect_list.iter().skip(m.read).copied());
        perfect_list.clear();
        perfect_list.extend(merged.drain(..));
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/dustmask/symdust.cpp:125-141
/// The body of the score test of `XTripletWindow::find_perfect` for iteration `i` (its `pos` and
/// `count` from `i`), in the same order.
#[inline(never)]
fn x_merge_at(
    m: &mut XMerge,
    i: usize,
    score: u32,
    merged: &mut Vec<PerfectInterval>,
    perfect_list: &VecDeque<PerfectInterval>,
) {
    let pos = m.pos0.saturating_sub(i);
    let count = m.count0 + i;
    while m.read < m.old_len && pos <= perfect_list[m.read].start {
        let p = perfect_list[m.read];
        if m.max_perfect_score == 0
            || m.max_len * p.score as usize > m.max_perfect_score as usize * p.len
        {
            m.max_perfect_score = p.score;
            m.max_len = p.len;
        }
        if m.inserted {
            merged.push(p);
        }
        m.read += 1;
    }
    if m.max_perfect_score == 0
        || score as usize * m.max_len >= m.max_perfect_score as usize * count
    {
        m.max_perfect_score = score;
        m.max_len = count;
        if !m.inserted {
            m.inserted = true;
            merged.reserve(m.old_len + 8);
            merged.extend(perfect_list.iter().take(m.read).copied());
        }
        merged.push(PerfectInterval::new(pos, m.end, m.max_perfect_score, count));
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
        // NCBI reference (598d8ae6): c++/include/algo/dustmask/symdust.hpp:75-84
        let code = |base: u8, rng: &mut NcbiLfgRandom| -> u8 {
            let c = X_BASE_CODE[base as usize];
            if c < 4 {
                c
            } else {
                (rng.next() & 0x3) as u8
            }
        };
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
                (code(seq[current_start], rng) << 2) + code(seq[current_start + 1], rng);
            let mut pos = current_start + 2;
            let mut done = false;

            while !done && pos < end {
                if !perfect_list.is_empty() {
                    self.save_masked_regions(
                        &mut result,
                        h.start,
                        current_start,
                        &mut perfect_list,
                    );
                }

                let new_triplet = ((current_triplet << 2) & 0x3F) + code(seq[pos], rng);
                current_triplet = new_triplet;
                pos += 1;

                if h.shift_window(&mut a, k, new_triplet, &mut perfect_list) {
                    if h.needs_processing(k) {
                        find_perfect(h, &a, k.thresholds, &mut merged, &mut perfect_list);
                    }
                } else {
                    while pos < end {
                        if !perfect_list.is_empty() {
                            self.save_masked_regions(
                                &mut result,
                                h.start,
                                current_start,
                                &mut perfect_list,
                            );
                        }

                        let new_triplet = ((current_triplet << 2) & 0x3F) + code(seq[pos], rng);
                        current_triplet = new_triplet;

                        if h.shift_window(&mut a, k, new_triplet, &mut perfect_list) {
                            done = true;
                            break;
                        }
                        pos += 1;
                    }
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
