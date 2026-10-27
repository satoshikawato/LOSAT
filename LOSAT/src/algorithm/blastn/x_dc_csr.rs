//! EXPERIMENT (LOSAT_X_DCCSR, in the strict set: default ON since 2026-10-21, off with
//! `LOSAT_X_DCCSR=0`): the chains of a
//! one-template discontiguous megablast table (`build_disc_mb_lookup`) laid out as CSR
//! (compressed sparse rows) instead of NCBI's `hashtable` / `next_pos`: `offsets[cell] ..
//! offsets[cell + 1]` indexes `qs`, which holds the cell's 1-based query offsets in the order
//! its chain lists them. A scan reads a cell's run as one slice, so a word with hits costs a
//! miss on `offsets` and one on its run, both prefetched ahead by `Ring`, instead of one miss
//! per chain element on `next_pos` (indexed by query offset, so in no order).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:775-786
//! ```c
//!          ecode1 = ComputeDiscontiguousIndex(accum, template_type);
//!          if (mb_lt->hashtable[ecode1] == 0) {
//!             PV_SET(pv_array, ecode1, pv_array_bts);
//!          }
//!          else {
//!             helper_array[ecode1/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos[index] = mb_lt->hashtable[ecode1];
//!          mb_lt->hashtable[ecode1] = index;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
//! ```c
//!     Int4 q_off = lookup->hashtable[index];
//!
//!     while (q_off) {
//!         offset_pairs[i].qs_offsets.q_off   = q_off - 1;
//!         offset_pairs[i++].qs_offsets.s_off = s_off;
//!         q_off = lookup->next_pos[q_off];
//!     }
//! ```
//!
//! Why the runs are the chains. The fill puts each word at the head of its cell's chain, so a
//! chain lists the cell's words newest first, each word (its `index`, unique) once. The CSR is
//! built from the same words in the same order (the loop of `fill_disc_mb_tables`: the same
//! ranges, ambiguity resets, `accum`, `index` and `ecode1`): a count pass adds 1 to
//! `offsets[ecode1]` per word, a sum pass turns the counts into inclusive running sums (the end
//! of each cell's run), and a place pass, in the fill's order again, writes each word at
//! `--offsets[ecode1]`. So each run holds its cell's words newest first, and after the place
//! pass `offsets[cell]` is the start of the cell's run and `offsets[cell + 1]` its end
//! (`offsets[cells]`, no cell's end, is set to the number of words). The sum pass also makes
//! the fill's presence vector and helper counters from the counts: the fill sets the bit of a
//! cell at its first word (`hashtable[ecode1] == 0`) and adds 1 to the cell's counter at every
//! further word, i.e. bit `ecode1 & 31` of word `ecode1 >> pv_array_bts` for every non-empty
//! cell and `count - 1` on counter `ecode1 / kCompressionFactor` (wrapping, as the fill's
//! `wrapping_add`); so `pv_array`, `helper_array` and `longest_chain` are the fill's. A scan
//! that reads `qs[offsets[c] .. offsets[c + 1]]` after the same presence test visits the same
//! query offsets in the same order as the chain walk, so the offset pairs, and everything made
//! from them, are the same. `contains_hit` (the masked-word test of `s_TypeOfWord`) searches
//! the run (sorted by decreasing offset) instead of walking the chain.
//!
//! Scope (`plan`): a one-template table (the two-template tables, `coding_and_optimal`, keep
//! their chains) whose word loop reads at least a quarter as many positions as the table has
//! cells (`wanted`). A smaller table, such as one of the one-query batches of NCBI's query split,
//! keeps the chains: there the chains fit the caches, the CSR saves no miss, and its zero and
//! sum passes touch all the cells (on the workstation the CSR of such a table cost six times the
//! fill's cycles, and the scan's prefetch ring 0.25 G cycles more over 19 of them; v4-R R15).
//! The table is built on the search pool when LOSAT_X_DISCLUTPAR's fill would be
//! (`x_disc_lookup_par::csr`), otherwise here (`build`). Its `hashtable` and `next_pos` are then
//! empty: the CSR replaces them (`qs` has one element per word, `next_pos` one per query
//! position), so the table takes no more memory. Every reader of the chains of a discontiguous
//! table goes through `DiscTemplates::x_csr`: the strip scan of LOSAT_X_PAIRPARDISC and its
//! diagonal-array stage (`scan_strip_disc`), the reference scan of `search_query_batch`
//! (`run.rs`, its `x_on_word_csr`; the chains' `on_word` stays as it was), and
//! `TwoStageLookup::contains_hit` and `count_hits` where the table has no chains.
//!
//! The scans of a CSR table read their words through `Ring`: each word enters a FIFO of `RING`
//! words and leaves it, unchanged and in order, `RING` words later; a word entering prefetches
//! `offsets[cell]`, and the word `RING / 2` places back reads its two offsets and prefetches
//! the first and the last cache line of its run (x86_64 only; elsewhere the ring only delays
//! the words). The words, and so every pair and its order, are the scan's own.
//!
//! With LOSAT_TIMING `main` prints `[TIMING] disc_lookup_csr: builds=<CSR tables> parallel=<on a
//! pool> words=<their words>` (the build times are in `[TIMING] disc_lookup_build`, the CSR's in
//! its `loop` or `pool_fill`).

use super::{build_unmasked_ranges, DiscWords, PvArrayType, BLAST2NA_MASK, PV_ARRAY_MASK};
#[cfg(test)]
use super::{disc_longest_chain, fill_disc_mb_tables, DiscTables};
use crate::algorithm::blastn::disc_lookup::compute_discontiguous_index;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads LOSAT_X_DCCSR once (a switch of the strict set: on unless
/// `LOSAT_X_DCCSR=0`). It chooses between two LOSAT layouts of the same chains and changes no
/// value NCBI computes.
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DCCSR"))
}

// No NCBI counterpart: counters of the summary lines.
static BUILDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_ON_POOL: AtomicU64 = AtomicU64::new(0);
static BUILT_WORDS: AtomicU64 = AtomicU64::new(0);

/// The `[TIMING]` line (with LOSAT_TIMING, when a CSR was built), printed by `main` at exit.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if std::env::var_os("LOSAT_TIMING").is_some() && load(&BUILDS) > 0 {
        eprintln!(
            "[TIMING] disc_lookup_csr: builds={} parallel={} words={}",
            load(&BUILDS),
            load(&BUILDS_ON_POOL),
            load(&BUILT_WORDS)
        );
    }
}

/// No NCBI counterpart (layout only): the chains of a table's first template as CSR (module
/// comment): `offsets` has one element per cell and one more, `qs` one per word.
pub struct Csr {
    offsets: Vec<u32>,
    qs: Vec<u32>,
}

// EXPERIMENT (LOSAT_X_BNSHRINKFREE): the two vectors are shrunk before they are freed, so that
// glibc's mmap threshold stays put (`x_shrink_free`).
impl Drop for Csr {
    fn drop(&mut self) {
        crate::algorithm::blastn::x_shrink_free::free(std::mem::take(&mut self.offsets));
        crate::algorithm::blastn::x_shrink_free::free(std::mem::take(&mut self.qs));
    }
}

impl Csr {
    /// The 1-based query offsets of the cell's chain, in its order.
    #[inline(always)]
    pub fn run(&self, cell: usize) -> &[u32] {
        &self.qs[self.offsets[cell] as usize..self.offsets[cell + 1] as usize]
    }

    /// Whether the cell's chain holds `q_off_1` (the run is sorted by decreasing offset).
    #[inline]
    pub fn contains(&self, cell: usize, q_off_1: u32) -> bool {
        self.run(cell).binary_search_by(|q| q_off_1.cmp(q)).is_ok()
    }

    /// The number of words.
    pub fn words(&self) -> usize {
        self.qs.len()
    }

    #[inline(always)]
    fn offsets(&self) -> &[u32] {
        &self.offsets
    }

    #[inline(always)]
    fn qs(&self) -> &[u32] {
        &self.qs
    }
}

/// No NCBI counterpart: how `build_disc_mb_lookup` lays out a table.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Layout {
    /// NCBI's chains (`hashtable`, `next_pos`).
    Chains,
    /// The CSR alone.
    Csr,
}

/// No NCBI counterpart: the layout of a table with `cells` cells (module comment, scope).
pub(super) fn plan(words: &DiscWords<'_>, cells: usize) -> Layout {
    if !on() {
        return Layout::Chains;
    }
    if words.two_templates {
        crate::utils::xreach::x_reach!("LOSAT_X_DCCSR", false);
        return Layout::Chains;
    }
    let csr = wanted(positions(words), cells);
    crate::utils::xreach::x_reach!("LOSAT_X_DCCSR", csr);
    if csr {
        Layout::Csr
    } else {
        Layout::Chains
    }
}

/// No NCBI counterpart: the size rule (module comment): at least one position of the word loop
/// per four cells.
fn wanted(positions: usize, cells: usize) -> bool {
    positions.saturating_mul(4) >= cells
}

/// No NCBI counterpart: the positions the word loop reads (the unmasked ranges of all queries;
/// every word ends at one of them).
fn positions(words: &DiscWords<'_>) -> usize {
    let mut total = 0usize;
    for (q_idx, seq) in words.queries_blastna.iter().enumerate() {
        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        for (start, end) in build_unmasked_ranges(seq.len(), masks) {
            total += end.saturating_sub(start);
        }
    }
    total
}

/// The word loop of `fill_disc_mb_tables` for the first template, without the updates: calls
/// `f(ecode1, index)` for each word, in the fill's order.
#[inline(always)]
fn each_word(words: &DiscWords<'_>, mut f: impl FnMut(usize, u32)) {
    let template_length = words.template_length;
    for (q_idx, seq) in words.queries_blastna.iter().enumerate() {
        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let query_offset = words.query_offsets[q_idx].max(0) as usize;
        for (range_start, range_end) in build_unmasked_ranges(seq.len(), masks) {
            let mut accum: u64 = 0;
            let mut valid_bases = 0usize;
            for pos in range_start..range_end {
                let val = seq[pos];
                if (val & BLAST2NA_MASK) != 0 {
                    accum = 0;
                    valid_bases = 0;
                    continue;
                }
                accum = (accum << 2) | val as u64;
                valid_bases += 1;
                if valid_bases < template_length {
                    continue;
                }
                // The 1-based start of the word in the concatenated query.
                let index = (query_offset + pos + 2 - template_length) as u32;
                let ecode1 = compute_discontiguous_index(accum, words.template_type) as usize;
                f(ecode1, index);
            }
        }
    }
}

/// No NCBI counterpart (scheduling only): a hint to bring the line of `p` into the cache; no
/// architectural access, so any address is allowed (x86_64 only).
#[inline(always)]
pub(super) fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not access memory architecturally; any address is allowed.
    unsafe {
        core::arch::x86_64::_mm_prefetch::<{ core::arch::x86_64::_MM_HINT_T0 }>(p as *const i8);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}

/// The sum pass over the cells `c0 .. c1` (`c0` a multiple of 2048, `pv_array_bts >= 5`):
/// `cnt.add(c)` holds cell `c`'s count and receives the inclusive running sum from `base`
/// (module comment); returns the running sum after `c1`. Each 32 cells give 32 presence bits,
/// OR-ed into word `c >> pv_array_bts` (if it is below `pv_len`) when the word's last group or
/// the range ends; each counter `c >> 11` gets the run's growth over its 2048 cells minus their
/// non-empty cells (`sum(count - (count != 0))`, wrapping), so the loop has no per-cell work
/// beyond the running sum and one bit.
///
/// SAFETY: the counts, presence-vector words and counters of the range are the caller's alone
/// (`c0` and `c1` on the boundaries of 2048 cells and of `2^pv_array_bts` cells, or `c1` the
/// end of the table), and the pointers are valid for them.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(super) unsafe fn sum_cells(
    pv: *mut PvArrayType,
    pv_len: usize,
    helper: *mut u32,
    cnt: *mut u32,
    base: u32,
    c0: usize,
    c1: usize,
    pv_array_bts: usize,
) -> u32 {
    debug_assert!(c0 % 2048 == 0 && pv_array_bts >= 5);
    let pv_last = (1usize << pv_array_bts) - 1;
    let (mut end, mut word, mut pop, mut h0) = (base, 0u32, 0u32, base);
    let mut g = c0;
    while g < c1 {
        let n = (c1 - g).min(32);
        let mut w32 = 0u32;
        for i in 0..n {
            // SAFETY: cell `g + i < c1` is the caller's.
            unsafe {
                let x = *cnt.add(g + i);
                w32 |= ((x != 0) as u32) << i;
                end = end.wrapping_add(x);
                *cnt.add(g + i) = end;
            }
        }
        g += n;
        word |= w32;
        pop += w32.count_ones();
        if g & pv_last == 0 || g >= c1 {
            let i = (g - 1) >> pv_array_bts;
            if i < pv_len {
                // SAFETY: presence-vector word `i` holds cells of the range only (or the range
                // ends in it).
                unsafe { *pv.add(i) |= word };
            }
            word = 0;
        }
        if g & 2047 == 0 || g >= c1 {
            let h = (g - 1) >> 11;
            // SAFETY: counter `h` holds cells of the range only (or the range ends in it).
            unsafe {
                *helper.add(h) =
                    (*helper.add(h)).wrapping_add(end.wrapping_sub(h0).wrapping_sub(pop))
            };
            h0 = end;
            pop = 0;
        }
    }
    end
}

/// The CSR whose `offsets` (`cells + 1` elements) and `qs` (`words` elements) have been
/// written in full by a build.
///
/// SAFETY: every element of `offsets[..cells + 1]` and `qs[..words]` is initialized, and
/// `words <= qs.capacity()`.
pub(super) unsafe fn assemble(
    mut offsets: Vec<u32>,
    mut qs: Vec<u32>,
    cells: usize,
    words: usize,
) -> Csr {
    debug_assert!(offsets.capacity() > cells && qs.capacity() >= words);
    // SAFETY: the caller's.
    unsafe {
        offsets.set_len(cells + 1);
        qs.set_len(words);
    }
    Csr { offsets, qs }
}

/// No NCBI counterpart (counters only): a CSR build of the search.
pub(super) fn count_build(csr: &Csr, on_pool: bool) {
    BUILDS.fetch_add(1, Ordering::Relaxed);
    BUILDS_ON_POOL.fetch_add(u64::from(on_pool), Ordering::Relaxed);
    BUILT_WORDS.fetch_add(csr.words() as u64, Ordering::Relaxed);
}

/// The serial build (module comment): the count pass (each word's `offsets` element prefetched
/// and counted 16 words later), the sum pass with the presence bits and helper counters into
/// `pv_array` and `helper_array` (zeroed, of the fill's lengths), and the place pass.
pub(super) fn build(
    words: &DiscWords<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
) -> Csr {
    // Once per table: the bounds the unchecked writes below rely on.
    assert!(!words.two_templates && cells % 2048 == 0 && pv_array_bts >= 5);
    assert_eq!(helper_array.len(), cells >> 11);
    let mut offsets: Vec<u32> = Vec::with_capacity(cells + 1);
    let ob = offsets.as_mut_ptr();
    // SAFETY: `cells + 1` elements of the capacity, written here before anything reads them.
    unsafe { std::ptr::write_bytes(ob, 0, cells + 1) };
    const AHEAD: usize = 16;
    let mut ring = [0usize; AHEAD];
    let mut seen = 0usize;
    each_word(words, |cell, _| {
        assert!(cell < cells);
        prefetch(ob.wrapping_add(cell));
        let slot = seen % AHEAD;
        if seen >= AHEAD {
            // SAFETY: `ring[slot] < cells` (asserted when it entered), an initialized element.
            unsafe { *ob.add(ring[slot]) += 1 };
        }
        ring[slot] = cell;
        seen += 1;
    });
    for i in seen.saturating_sub(AHEAD)..seen {
        // SAFETY: as above.
        unsafe { *ob.add(ring[i % AHEAD]) += 1 };
    }
    let n_words = seen;
    assert!(
        n_words <= u32::MAX as usize,
        "DCCSR: more words than u32 offsets"
    );
    // SAFETY: the counts of all cells are this build's; the presence vector and the counters
    // are of the fill's lengths (`c >> 11 < helper_array.len()` for every cell).
    let total = unsafe {
        sum_cells(
            pv_array.as_mut_ptr(),
            pv_array.len(),
            helper_array.as_mut_ptr(),
            ob,
            0,
            0,
            cells,
            pv_array_bts,
        )
    };
    debug_assert_eq!(total as usize, n_words);
    // SAFETY: `cells < cells + 1`; `offsets[cells]` is no cell's end.
    unsafe { *ob.add(cells) = total };
    let mut qs: Vec<u32> = Vec::with_capacity(n_words);
    let qb = qs.as_mut_ptr();
    each_word(words, |cell, index| {
        // SAFETY: `cell < cells` (the same words as the count pass), `offsets[cell]` is the end
        // of the cell's non-empty run (>= 1), so the place `e < n_words` is in the capacity.
        unsafe {
            let e = *ob.add(cell) - 1;
            *ob.add(cell) = e;
            *qb.add(e as usize) = index;
        }
    });
    // SAFETY: `offsets` was written in full; the runs partition `0 .. n_words`, and the place
    // pass wrote each place once.
    unsafe { assemble(offsets, qs, cells, n_words) }
}

/// No NCBI counterpart (scheduling only): a spin barrier for the `n` threads of one
/// `rayon::broadcast` (`x_disc_lookup_par::csr`): the last of the `n` arrivals resets the count
/// and bumps the generation, the others wait until it changes (spinning, then yielding the CPU
/// on native targets; WebAssembly only spins). A thread that panics between barriers (`Poison`,
/// unwinding builds only: the release profile aborts) makes the waiting threads panic instead of
/// waiting for it forever.
pub(super) struct SpinBarrier {
    n: usize,
    count: AtomicUsize,
    generation: AtomicUsize,
    poisoned: std::sync::atomic::AtomicBool,
}

/// Held by each thread of a build that waits at a `SpinBarrier`: poisons the barrier when the
/// thread unwinds.
pub(super) struct Poison<'a>(pub(super) &'a SpinBarrier);

impl Drop for Poison<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.poisoned.store(true, Ordering::Release);
        }
    }
}

impl SpinBarrier {
    pub(super) fn new(n: usize) -> Self {
        Self {
            n,
            count: AtomicUsize::new(0),
            generation: AtomicUsize::new(0),
            poisoned: std::sync::atomic::AtomicBool::new(false),
        }
    }

    pub(super) fn wait(&self) {
        let g = self.generation.load(Ordering::Acquire);
        if self.count.fetch_add(1, Ordering::AcqRel) + 1 == self.n {
            self.count.store(0, Ordering::Relaxed);
            self.generation.store(g.wrapping_add(1), Ordering::Release);
            return;
        }
        let mut spins = 0u32;
        while self.generation.load(Ordering::Acquire) == g {
            if self.poisoned.load(Ordering::Acquire) {
                panic!("DCCSR: another thread of the table's build panicked");
            }
            if cfg!(target_arch = "wasm32") || spins < 1 << 16 {
                std::hint::spin_loop();
                spins += 1;
            } else {
                std::thread::yield_now();
            }
        }
    }
}

/// No NCBI counterpart (scheduling only): the depth of the scans' FIFO of words (module comment).
pub(crate) const RING: usize = 32;

/// No NCBI counterpart (scheduling only): the FIFO of a CSR scan's words (module comment). A
/// word is any `Copy` value whose cell `push` gets from `cell`.
pub(crate) struct Ring<'a, T: Copy> {
    buf: [T; RING],
    n: usize,
    offsets: &'a [u32],
    qs: &'a [u32],
}

impl<'a, T: Copy + Default> Ring<'a, T> {
    pub(crate) fn new(csr: &'a Csr) -> Self {
        Self {
            buf: [T::default(); RING],
            n: 0,
            offsets: csr.offsets(),
            qs: csr.qs(),
        }
    }

    /// The word `w` enters (its `offsets` element prefetched; the word `RING / 2` places back
    /// gets the first and the last line of its run prefetched), and the word `RING` places back
    /// leaves through `emit`.
    #[inline(always)]
    pub(crate) fn push(&mut self, w: T, cell: impl Fn(&T) -> usize, emit: &mut impl FnMut(T)) {
        const MASK: usize = RING - 1;
        let offsets = self.offsets;
        prefetch(offsets.as_ptr().wrapping_add(cell(&w)));
        if self.n >= RING / 2 {
            let c = cell(&self.buf[(self.n - RING / 2) & MASK]);
            let a = offsets[c] as usize;
            let b = offsets[c + 1] as usize;
            if b > a {
                let pa = self.qs.as_ptr().wrapping_add(a);
                let pb = self.qs.as_ptr().wrapping_add(b - 1);
                prefetch(pa);
                if (pa as usize >> 6) != (pb as usize >> 6) {
                    prefetch(pb);
                }
            }
        }
        let slot = self.n & MASK;
        if self.n >= RING {
            emit(self.buf[slot]);
        }
        self.buf[slot] = w;
        self.n += 1;
    }

    /// The words still in the FIFO leave through `emit`, in order.
    pub(crate) fn flush(&mut self, emit: &mut impl FnMut(T)) {
        for m in self.n.saturating_sub(RING)..self.n {
            emit(self.buf[m & (RING - 1)]);
        }
        self.n = 0;
    }
}

const _: () = assert!(RING.is_power_of_two() && RING >= 4);
const _: () = assert!(PV_ARRAY_MASK == 31);

/// No NCBI counterpart: the check of the unit tests of a table built with both layouts: the CSR
/// `csr` (built on the pool when `on_pool`; its presence vector and counters in `pv` and
/// `helper`) against the chains, presence vector, counters and longest chain of the fill
/// (`hashtable`, `next_pos`, `fill_pv`, `fill_helper`, `longest_chain`), and, when `on_pool`,
/// against the serial build. Panics on a difference.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn x_assert_runs_are_the_chains(
    words: &DiscWords<'_>,
    csr: &Csr,
    on_pool: bool,
    pv: &[PvArrayType],
    helper: &[u32],
    pv_array_bts: usize,
    hashtable: &[u32],
    next_pos: &[u32],
    fill_pv: &[PvArrayType],
    fill_helper: &[u32],
    longest_chain: usize,
) {
    let cells = hashtable.len();
    let mut diffs = Vec::new();
    if csr.offsets.len() != cells + 1 {
        diffs.push(format!("offsets len {} cells {cells}", csr.offsets.len()));
    } else {
        let mut chain = Vec::new();
        for (cell, &head) in hashtable.iter().enumerate() {
            chain.clear();
            let mut q = head;
            while q != 0 {
                chain.push(q);
                q = next_pos[q as usize];
            }
            if csr.run(cell) != chain.as_slice() {
                diffs.push(format!(
                    "cell {cell} run {:?} chain {:?}",
                    &csr.run(cell)[..csr.run(cell).len().min(8)],
                    &chain[..chain.len().min(8)]
                ));
                break;
            }
            if csr.run(cell).windows(2).any(|w| w[0] <= w[1]) {
                diffs.push(format!("cell {cell} run not decreasing"));
                break;
            }
            if let Some(&q) = chain.first() {
                if !csr.contains(cell, q) || csr.contains(cell, 0) {
                    diffs.push(format!("cell {cell} contains"));
                    break;
                }
            }
        }
    }
    if pv != fill_pv {
        diffs.push("presence vector".to_string());
    }
    if helper != fill_helper {
        diffs.push("helper counters".to_string());
    }
    let csr_longest = disc_longest_chain(helper, &[], false);
    if csr_longest != longest_chain {
        diffs.push(format!("longest_chain {csr_longest} fill {longest_chain}"));
    }
    if on_pool {
        let mut pv2 = vec![0 as PvArrayType; pv.len()];
        let mut helper2 = vec![0u32; helper.len()];
        let serial = build(words, cells, &mut pv2, pv_array_bts, &mut helper2);
        if serial.offsets != csr.offsets || serial.qs != csr.qs || pv2 != pv || helper2 != helper {
            diffs.push("pool build != serial build".to_string());
        }
    }
    assert!(
        diffs.is_empty(),
        "the CSR differs from the chains (cells {cells}, words {}, pool {on_pool}): {}",
        csr.words(),
        diffs.join("; ")
    );
}

/// The chains of the first template filled by the original loop, for the tests: `(hashtable,
/// next_pos, pv_array, helper_array)`.
#[cfg(test)]
fn chains(
    words: &DiscWords<'_>,
    cells: usize,
    pv_len: usize,
    pv_array_bts: usize,
) -> (Vec<u32>, Vec<u32>, Vec<PvArrayType>, Vec<u32>) {
    let query_length = words
        .queries_blastna
        .iter()
        .zip(words.query_offsets)
        .map(|(seq, &off)| off.max(0) as usize + seq.len())
        .max()
        .unwrap_or(0);
    let mut hashtable = vec![0u32; cells];
    let mut next_pos = vec![0u32; query_length + 1];
    let mut pv = vec![0 as PvArrayType; pv_len];
    let mut helper = vec![0u32; cells >> 11];
    fill_disc_mb_tables(
        words,
        DiscTables {
            hashtable: &mut hashtable,
            next_pos: &mut next_pos,
            helper_array: &mut helper,
            hashtable2: &mut [],
            next_pos2: &mut [],
            helper_array2: &mut [],
            pv_array: &mut pv,
            pv_array_bts,
        },
    );
    (hashtable, next_pos, pv, helper)
}

#[cfg(test)]
mod tests {
    use super::super::{DiscTemplateType, MaskedInterval};
    use super::*;

    /// A small deterministic generator (xorshift64*).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// Random BLASTNA queries over few distinct letters (so cells get long chains), with
    /// ambiguities and masks.
    fn queries(
        rng: &mut Rng,
        n: usize,
        len: usize,
    ) -> (Vec<Vec<u8>>, Vec<i32>, Vec<Vec<MaskedInterval>>) {
        let mut seqs = Vec::new();
        let mut offsets = Vec::new();
        let mut masks = Vec::new();
        let mut at = 0i32;
        for q in 0..n {
            let l = len / 2 + rng.below(len as u64) as usize;
            let letters = if q % 2 == 0 { 2 } else { 4 };
            let seq: Vec<u8> = (0..l)
                .map(|_| {
                    if rng.below(200) == 0 {
                        14 // N, an ambiguity
                    } else {
                        rng.below(letters) as u8
                    }
                })
                .collect();
            let mut m = Vec::new();
            if q % 3 == 1 && l > 40 {
                let s = rng.below((l - 20) as u64) as usize;
                m.push(MaskedInterval {
                    start: s,
                    end: s + 15,
                });
            }
            offsets.push(at);
            at += l as i32 + 1;
            seqs.push(seq);
            masks.push(m);
        }
        (seqs, offsets, masks)
    }

    fn words<'a>(
        seqs: &'a [&'a [u8]],
        offsets: &'a [i32],
        masks: &'a [Vec<MaskedInterval>],
        template_length: usize,
        template_type: DiscTemplateType,
    ) -> DiscWords<'a> {
        DiscWords {
            queries_blastna: seqs,
            query_offsets: offsets,
            query_masks: masks,
            template_length,
            template_type,
            second_template_type: DiscTemplateType::Contiguous,
            two_templates: false,
        }
    }

    #[test]
    fn the_runs_are_the_chains_with_the_fills_presence_vector_and_counters() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let cases = [
            (DiscTemplateType::T11_16Coding, 16, 22usize, 5usize),
            (DiscTemplateType::T11_18Optimal, 18, 22, 7),
            (DiscTemplateType::T12_21Coding, 21, 24, 8),
        ];
        for (template_type, template_length, log_cells, pv_array_bts) in cases {
            let cells = 1usize << log_cells;
            let pv_len = cells >> pv_array_bts;
            for (n, len) in [(1, 30), (7, 400), (40, 3000)] {
                let (seqs, offsets, masks) = queries(&mut rng, n, len);
                let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
                let w = words(&refs, &offsets, &masks, template_length, template_type);
                let (hashtable, next_pos, fill_pv, fill_helper) =
                    chains(&w, cells, pv_len, pv_array_bts);
                let mut pv = vec![0 as PvArrayType; pv_len];
                let mut helper = vec![0u32; cells >> 11];
                let csr = build(&w, cells, &mut pv, pv_array_bts, &mut helper);
                let longest = disc_longest_chain(&fill_helper, &[], false);
                x_assert_runs_are_the_chains(
                    &w,
                    &csr,
                    false,
                    &pv,
                    &helper,
                    pv_array_bts,
                    &hashtable,
                    &next_pos,
                    &fill_pv,
                    &fill_helper,
                    longest,
                );
                let chained: usize = (0..cells).map(|c| csr.run(c).len()).sum();
                assert_eq!(chained, csr.words());
                assert!(csr.words() > 0 || n == 1);
            }
        }
    }

    #[test]
    fn contains_finds_exactly_the_runs_offsets() {
        let csr = Csr {
            offsets: vec![0, 0, 3, 4, 4],
            qs: vec![90, 40, 7, 12],
        };
        for (cell, q, want) in [
            (0, 7, false),
            (1, 90, true),
            (1, 40, true),
            (1, 7, true),
            (1, 41, false),
            (1, 6, false),
            (2, 12, true),
            (2, 7, false),
            (3, 12, false),
        ] {
            assert_eq!(csr.contains(cell, q), want, "{cell} {q}");
        }
    }

    #[test]
    fn the_ring_emits_every_word_once_in_order() {
        let csr = Csr {
            offsets: vec![0, 2, 2, 5],
            qs: vec![9, 3, 30, 20, 10],
        };
        for n in [0usize, 1, 15, 16, 17, 31, 32, 33, 100] {
            let mut ring: Ring<'_, (u32, usize)> = Ring::new(&csr);
            let mut out = Vec::new();
            for i in 0..n {
                ring.push((i as u32, i % 3), |w| w.1, &mut |w| out.push(w.0));
            }
            ring.flush(&mut |w| out.push(w.0));
            assert_eq!(out, (0..n as u32).collect::<Vec<_>>());
        }
    }

    #[test]
    fn the_size_rule_takes_a_quarter_of_the_cells() {
        assert!(!wanted(0, 1 << 22));
        assert!(!wanted((1 << 20) - 1, 1 << 22));
        assert!(wanted(1 << 20, 1 << 22));
        assert!(wanted(usize::MAX, 1 << 22));
    }
}
