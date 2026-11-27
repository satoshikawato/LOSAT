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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
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
//! cell and `count - 1` on counter `ecode1` divided by `kCompressionFactor` (wrapping, as the fill's
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
//! diagonal-array stage (`scan_strip_disc`), the reference scan of a two-stage table
//! (`two_stage_scan.rs` `scan_two_stage_lookup`, its `x_on_word_csr`; the chains' `on_word`
//! stays as it was), and
//! `TwoStageLookup::contains_hit` and `count_hits` where the table has no chains.
//!
//! EXPERIMENT (LOSAT_X_DCQSAHEAD, in the strict set: on unless `LOSAT_X_DCQSAHEAD=0`): the place
//! pass of a table built on the pool by `x_disc_lookup_par::csr_on_pool_steps` (LOSAT_X_DCSEEDLEAN's
//! build, at two threads or more) brings in the line of `qs` each entry will be written to before
//! it writes it. The place pass reads `offsets[cell]`, decrements it and writes `qs` there; it
//! brought in the `offsets` element of the entry 16 places ahead, but not the `qs` line, and `qs`
//! has one element per word (tens of MiB for a large query set), so a store's line was often not
//! in the caches. Now the owner also reads the `offsets` element of the entry `QS_AHEAD` (8)
//! places ahead (brought in 8 entries before) and brings in the `qs` line it points to. The
//! entries are placed in the same order at the same places (the value read ahead only chooses
//! which line is brought in), so the CSR is the same; only the bringing in is new (`prefetch`,
//! nothing on targets without one). Its sides are DCSEEDLEAN's pooled build's: the serial build
//! (one thread, no pool) and `csr_on_pool` (`LOSAT_X_DCSEEDLEAN=0`) place as before.
//!
//! The scans of a CSR table read their words through `Ring`: each word enters a FIFO of `RING`
//! words and leaves it, unchanged and in order, `RING` words later. Every `GROUP` words, the
//! group that just entered gets its `offsets[cell]` brought in, and the group `RING / 2` places
//! back reads its offsets and gets the first and the last cache line of each run brought in
//! (`x_ahead`: a prefetch on x86_64; on wasm32, which has no prefetch instruction, a load whose
//! value is dropped, so the group's misses are waited for together and not one after another
//! where `emit` reads them; elsewhere nothing). The words, and so every pair and its order, are
//! the scan's own. A ring is made
//! for each chunk a scan reads, so its slots start uninitialised (no 0.8-0.9 KB zero fill per
//! chunk): a slot is read only after a word was written to it (`Ring::push`).
//!
//! With LOSAT_TIMING `main` prints `[TIMING] disc_lookup_csr: builds=<CSR tables> parallel=<on a
//! pool> words=<their words> qs_ahead=<tables whose place pass brought in qs ahead>` (the build
//! times are in `[TIMING] disc_lookup_build`, the CSR's in its `loop` or `pool_fill`).

use super::{build_unmasked_ranges, x_range_words, DiscWords, PvArrayType, PV_ARRAY_MASK};
#[cfg(test)]
use super::{disc_longest_chain, fill_disc_mb_tables, DiscTables};
use crate::algorithm::blastn::disc_lookup::compute_discontiguous_index;
use std::collections::HashMap;
use std::mem::MaybeUninit;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads LOSAT_X_DCCSR once (a switch of the strict set: on unless
/// `LOSAT_X_DCCSR=0`). It chooses between two LOSAT layouts of the same chains and changes no
/// value NCBI computes.
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DCCSR"))
}

/// No NCBI counterpart: reads LOSAT_X_DCQSAHEAD once (a switch of the strict set: on unless
/// `LOSAT_X_DCQSAHEAD=0`): the place pass of `csr_on_pool_steps` brings in `qs` ahead (module
/// comment). It changes when lines are brought in, not what is written where, so it changes no
/// value NCBI computes.
pub(super) fn qs_ahead_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_DCQSAHEAD"))
}

/// No NCBI counterpart (scheduling only, LOSAT_X_DCQSAHEAD): how many entries before its place an
/// entry's `qs` line is brought in; its `offsets` element was brought in `16 - QS_AHEAD` entries
/// before that.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) const QS_AHEAD: usize = 8;

// No NCBI counterpart: counters of the summary lines.
static BUILDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_ON_POOL: AtomicU64 = AtomicU64::new(0);
static BUILT_WORDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_QS_AHEAD: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart (counters only, LOSAT_X_DCQSAHEAD): a build whose place pass brought in
/// `qs` ahead.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn count_qs_ahead() {
    BUILDS_QS_AHEAD.fetch_add(1, Ordering::Relaxed);
}

/// The `[TIMING]` line (with LOSAT_TIMING, when a CSR was built), printed by `main` at exit.
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if std::env::var_os("LOSAT_TIMING").is_some() && load(&BUILDS) > 0 {
        eprintln!(
            "[TIMING] disc_lookup_csr: builds={} parallel={} words={} qs_ahead={}",
            load(&BUILDS),
            load(&BUILDS_ON_POOL),
            load(&BUILT_WORDS),
            load(&BUILDS_QS_AHEAD)
        );
    }
}

/// No NCBI counterpart (layout only): the chains of a table's first template as CSR (module
/// comment): `offsets` has one element per cell and one more, `qs` one per word. With
/// LOSAT_X_BNCSRPACK (`x_bn_csr.rs`) a contiguous table's head is packed: `x_bits` holds two
/// words per 64 cells (the cells' presence bits, then the number of non-empty cells before them)
/// and `offsets` one element per non-empty cell and one more, indexed by that rank.
pub struct Csr {
    offsets: Vec<u32>,
    /// EXPERIMENT (LOSAT_X_BNCSRPACK): empty for a head of one element per cell.
    x_bits: Vec<u64>,
    qs: Vec<u32>,
    /// EXPERIMENT (LOSAT_X_MBLUTKEEP): the batch loop's keep that takes `offsets` back when the CSR
    /// is dropped, with `Keep::take_buffer`'s length and written elements (`x_kept`), and
    /// whether it is offered `qs` too (LOSAT_X_CSRQSKEEP, `x_kept_with_qs`).
    x_keep: Option<(
        std::sync::Arc<super::x_mb_lut_keep::Keep>,
        usize,
        usize,
        bool,
    )>,
}

// EXPERIMENT (LOSAT_X_MBLUTKEEP): `offsets` goes back to the batch loop's keep when the build
// took its buffer from it (`x_kept`); EXPERIMENT (LOSAT_X_CSRQSKEEP) then, for a discontiguous
// build (`x_kept_with_qs`), `qs` (with the switch off, or once the keep is closed for the loop's
// last tables, it stays here and is freed with the CSR, as before).
// No NCBI counterpart (memory only): where the CSR's arrays go when its batch drops it.
impl Drop for Csr {
    fn drop(&mut self) {
        if let Some((keep, length, dirty, qs_back)) = self.x_keep.take() {
            keep.give_back_buffer(std::mem::take(&mut self.offsets), length, dirty);
            if qs_back {
                keep.give_back_qs(&mut self.qs);
            }
        }
    }
}

impl Csr {
    /// EXPERIMENT (LOSAT_X_MBLUTKEEP): the CSR whose `offsets` goes back to `keep` when it is
    /// dropped (`take` = `Keep::take_buffer`'s length and written elements). Its `qs` is freed
    /// with it (a contiguous build, `x_bn_csr`).
    pub(super) fn x_kept(
        mut self,
        keep: &std::sync::Arc<super::x_mb_lut_keep::Keep>,
        take: (usize, usize),
    ) -> Self {
        self.x_keep = Some((std::sync::Arc::clone(keep), take.0, take.1, false));
        self
    }

    /// EXPERIMENT (LOSAT_X_CSRQSKEEP): `x_kept` for a discontiguous build, which took the kept `qs`
    /// (`Keep::take_qs`): its `qs` goes back to `keep` after its `offsets` (`Keep::give_back_qs`).
    /// No NCBI counterpart (memory only).
    pub(super) fn x_kept_with_qs(
        mut self,
        keep: &std::sync::Arc<super::x_mb_lut_keep::Keep>,
        take: (usize, usize),
    ) -> Self {
        self.x_keep = Some((std::sync::Arc::clone(keep), take.0, take.1, true));
        self
    }

    /// The 1-based query offsets of the cell's chain, in its order. A packed head
    /// (LOSAT_X_BNCSRPACK) answers through `x_packed` instead.
    #[inline(always)]
    pub fn run(&self, cell: usize) -> &[u32] {
        debug_assert!(
            self.x_bits.is_empty(),
            "BNCSRPACK: a packed head reads through x_packed"
        );
        &self.qs[self.offsets[cell] as usize..self.offsets[cell + 1] as usize]
    }

    /// EXPERIMENT (LOSAT_X_BNCSRPACK): the reader of a packed head, or `None` for a head of one
    /// element per cell. The contiguous table's scans choose the form here once per table, so
    /// the scans of each form compile on their own (`XHead`), and those of a head of one element
    /// per cell, every discontiguous table's among them, as without the switch.
    #[inline(always)]
    pub(crate) fn x_packed(&self) -> Option<XPacked<'_>> {
        (!self.x_bits.is_empty()).then(|| XPacked {
            bits: &self.x_bits,
            offsets: &self.offsets,
            qs: &self.qs,
            pop: XPop::get(),
        })
    }

    /// EXPERIMENT (LOSAT_X_BNCSRPACK): `run` of a head of either form, for the readers outside
    /// the scans (the lookup's answers for one word, checks and tests).
    pub(crate) fn x_run_any(&self, cell: usize) -> &[u32] {
        match self.x_packed() {
            Some(packed) => packed.run(cell),
            None => self.run(cell),
        }
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

    /// No NCBI counterpart (tests only): whether the two CSRs hold the same arrays (`offsets`,
    /// packed head bits and `qs`, each to its length).
    #[cfg(test)]
    pub(super) fn x_same_arrays(&self, other: &Csr) -> bool {
        self.offsets == other.offsets && self.x_bits == other.x_bits && self.qs == other.qs
    }

    /// No NCBI counterpart (tests only, LOSAT_X_CSRQSKEEP): the address of `qs`'s buffer.
    #[cfg(test)]
    pub(super) fn x_qs_ptr(&self) -> *const u32 {
        self.qs.as_ptr()
    }
}

/// No NCBI counterpart (tests only, LOSAT_X_CSRQSKEEP): `build` (one CSR build, given the kept `qs`
/// it passes to `x_mb_lut_keep::qs_buffer`) on a kept `qs` of stale values with exactly the room
/// that `fresh` (the same build given `None`) allocated, which it takes, and on one an element
/// shorter, which it frees for a fresh one: both give `fresh`'s arrays (a build that allocated
/// nothing, with no words, gets the keep's empty vector only). Panics on a difference.
#[cfg(test)]
pub(super) fn x_assert_same_on_a_kept_qs(
    fresh: &Csr,
    mut build: impl FnMut(Option<Vec<u32>>) -> Csr,
    what: &str,
) {
    const STALE: u32 = 0xA5A5_A5A5;
    let room = fresh.qs.capacity();
    if room == 0 {
        assert!(
            build(Some(Vec::new())).x_same_arrays(fresh),
            "{what}: none kept"
        );
        return;
    }
    let kept = vec![STALE; room];
    let ptr = kept.as_ptr();
    let on_kept = build(Some(kept));
    assert!(on_kept.x_same_arrays(fresh), "{what}: on a kept qs");
    assert_eq!(on_kept.qs.as_ptr(), ptr, "{what}: the kept qs taken");
    let short = build(Some(vec![STALE; room - 1]));
    assert!(short.x_same_arrays(fresh), "{what}: on a kept qs too small");
    assert!(short.qs.capacity() >= room, "{what}: a fresh qs");
}

/// No NCBI counterpart (layout only): EXPERIMENT (LOSAT_X_BNCSRPACK) the head of a CSR in one of
/// its forms, chosen once per table (`Csr::x_packed`): `Csr` (one element per cell, `Csr::run`
/// and `Ring`) or `XPacked`; a scan generic over it compiles once per form, and reads the head
/// through a reference as without the switch.
pub(crate) trait XHead {
    /// The scans' FIFO of words over this form.
    type Ring<'a, T: Copy + Default>: XRing<T>
    where
        Self: 'a;
    /// The 1-based query offsets of the cell's chain, in its order.
    fn run(&self, cell: usize) -> &[u32];
    /// EXPERIMENT (LOSAT_X_BNRANKONCE): the run of the head's element `e` (`offsets[e] ..
    /// offsets[e + 1]`): a cell of a head of one element per cell, a rank of a packed head; the
    /// element the head's ring hands `emit` with each word (`XRing::push`).
    fn run_at(&self, e: usize) -> &[u32];
    fn ring<T: Copy + Default>(&self) -> Self::Ring<'_, T>;
}

/// EXPERIMENT (LOSAT_X_BNCSRPACK): `Ring`'s `push` and `flush`, for a scan generic over `XHead`.
/// EXPERIMENT (LOSAT_X_BNRANKONCE): `emit` gets each word with its run's element of the head
/// (`XHead::run_at`), so a packed head's word is ranked once, where it enters.
pub(crate) trait XRing<T> {
    /// `Ring::push`: whether a word left through `emit`.
    fn push(&mut self, w: T, cell: impl Fn(&T) -> usize, emit: &mut impl FnMut(T, usize)) -> bool;
    fn flush(&mut self, cell: impl Fn(&T) -> usize, emit: &mut impl FnMut(T, usize));
}

impl XHead for Csr {
    type Ring<'a, T: Copy + Default> = Ring<'a, T>;

    #[inline(always)]
    fn run(&self, cell: usize) -> &[u32] {
        Csr::run(self, cell)
    }

    #[inline(always)]
    fn run_at(&self, e: usize) -> &[u32] {
        Csr::run(self, e)
    }

    #[inline(always)]
    fn ring<T: Copy + Default>(&self) -> Ring<'_, T> {
        Ring::new(self)
    }
}

/// No NCBI counterpart (layout only): EXPERIMENT (LOSAT_X_BNCSRPACK) the reader of a packed head
/// (`Csr::x_packed`; the layout in `Csr`'s comment).
#[derive(Clone, Copy)]
pub(crate) struct XPacked<'a> {
    bits: &'a [u64],
    offsets: &'a [u32],
    qs: &'a [u32],
    /// EXPERIMENT (LOSAT_X_BNPOPCNT): the form of the counts of its ranks.
    pop: XPop,
}

impl<'a> XPacked<'a> {
    /// The bytes of the presence bits and their ranks, which a scan that tests each word here reads.
    #[inline(always)]
    pub(crate) fn x_bits_bytes(self) -> usize {
        std::mem::size_of_val(self.bits)
    }

    /// Whether `cell` has a run: its presence bit, `x_rank`'s test without the count. For a scan
    /// that drops a word without hits before it calls out (`two_stage_scan`'s `x_scan_csr!`).
    #[inline(always)]
    pub(crate) fn x_present(self, cell: usize) -> bool {
        self.bits[(cell >> 6) << 1] & (1u64 << (cell & 63)) != 0
    }

    /// `Csr::run` of a packed head: an empty cell has no element of it.
    #[inline(always)]
    pub(crate) fn run(self, cell: usize) -> &'a [u32] {
        match x_rank(self.bits, cell, self.pop) {
            Some(r) => self.run_at(r),
            None => &[],
        }
    }

    /// EXPERIMENT (LOSAT_X_BNRANKONCE): the run of the non-empty cell of rank `r`.
    #[inline(always)]
    pub(crate) fn run_at(self, r: usize) -> &'a [u32] {
        &self.qs[self.offsets[r] as usize..self.offsets[r + 1] as usize]
    }
}

impl XHead for XPacked<'_> {
    type Ring<'a, T: Copy + Default>
        = XRingPacked<'a, T>
    where
        Self: 'a;

    #[inline(always)]
    fn run(&self, cell: usize) -> &[u32] {
        XPacked::run(*self, cell)
    }

    #[inline(always)]
    fn run_at(&self, r: usize) -> &[u32] {
        XPacked::run_at(*self, r)
    }

    #[inline(always)]
    fn ring<T: Copy + Default>(&self) -> XRingPacked<'_, T> {
        XRingPacked {
            buf: [MaybeUninit::uninit(); RING],
            ranks: [MaybeUninit::uninit(); RING],
            n: 0,
            head: *self,
        }
    }
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNPOPCNT, in the strict set: on unless it is `0`) how
/// the packed heads count presence bits (`x_rank`, `x_presence8`, `x_sum_packed`): with x86-64's
/// `popcnt` instruction where the CPU has it, else with `count_ones`. The default x86-64 target has
/// no `popcnt`, so there `count_ones` is a dozen instructions; on aarch64 and in Wasm it is the
/// target's own count. Both forms count the same bits. Read once per table build or scan (`get`)
/// and carried by value, so the hot loops test a flag in a register.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct XPop(bool);

impl XPop {
    /// The form of this run: `popcnt` when the switch is on and the CPU has it.
    #[inline]
    pub(crate) fn get() -> Self {
        static ON: OnceLock<bool> = OnceLock::new();
        XPop(*ON.get_or_init(|| Self::x_cpu() && crate::utils::xswitch::on("LOSAT_X_BNPOPCNT")))
    }

    /// Whether this CPU has the `popcnt` instruction.
    fn x_cpu() -> bool {
        #[cfg(target_arch = "x86_64")]
        {
            std::arch::is_x86_feature_detected!("popcnt")
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            false
        }
    }

    /// Whether the counts use the `popcnt` instruction (`LOSAT_TIMING`'s `x_bn_csr` line).
    pub(crate) fn x_instruction(self) -> bool {
        self.0
    }

    /// The forms this CPU can run (`count_ones`, then `popcnt` if it has it), for the tests.
    #[cfg(test)]
    pub(super) fn x_forms() -> Vec<XPop> {
        let mut forms = vec![XPop(false)];
        if Self::x_cpu() {
            forms.push(XPop(true));
        }
        forms
    }

    /// The number of bits set in `w`.
    #[inline(always)]
    pub(crate) fn ones(self, w: u64) -> usize {
        #[cfg(target_arch = "x86_64")]
        if self.0 {
            let mut n = w;
            // SAFETY: `self.0` holds only where the CPU has `popcnt` (`get`, `x_forms`). The
            // instruction reads and writes one register (the same one: no false dependency on
            // another register's last writer), no memory and no stack.
            unsafe {
                std::arch::asm!("popcnt {n}, {n}", n = inout(reg) n, options(pure, nomem, nostack))
            };
            return n as usize;
        }
        w.count_ones() as usize
    }
}

/// No NCBI counterpart (layout only): EXPERIMENT (LOSAT_X_BNCSRPACK) the element of a packed head
/// that holds `cell`'s run start (the number of non-empty cells before it), or `None` when the cell
/// is empty. `bits`: two words per 64 cells, the presence bits and the rank of the first. The bits
/// below the cell's are counted in the form `pop` (LOSAT_X_BNPOPCNT).
#[inline(always)]
pub(super) fn x_rank(bits: &[u64], cell: usize, pop: XPop) -> Option<usize> {
    let j = (cell >> 6) << 1;
    let w = bits[j];
    let m = 1u64 << (cell & 63);
    if w & m == 0 {
        return None;
    }
    Some(bits[j + 1] as usize + pop.ones(w & (m - 1)))
}

/// No NCBI counterpart (layout only): EXPERIMENT (LOSAT_X_BNCSRPACK) the presence bits of 64 cells
/// from their counts `c` (`x_count8`): bit `k` is set when `c[k]` is not 0. Eight bytes at a time:
/// the high bit of a byte of `((x & 0x7f..) + 0x7f..) | x` is set when the byte is not 0 (no sum
/// leaves its byte), and the product of those bits (bit `8 k`) by `0x0102040810204080` puts byte
/// `k`'s at bit `56 + k` (the products' bits are distinct, so nothing carries).
#[inline(always)]
pub(super) fn x_nonzero_bits(c: &[u8; 64]) -> u64 {
    const LOW7: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    const GATHER: u64 = 0x0102_0408_1020_4080;
    let mut m = 0u64;
    for (i, b) in c.chunks_exact(8).enumerate() {
        let x = u64::from_le_bytes(b.try_into().unwrap());
        let t = (((x & LOW7) + LOW7) | x) & !LOW7;
        m |= ((t >> 7).wrapping_mul(GATHER) >> 56) << (8 * i);
    }
    m
}

/// No NCBI counterpart (layout only): EXPERIMENT (LOSAT_X_BNCSRPACK) one more word of `cell` in the
/// counts of a packed build: one byte per cell up to `u8::MAX`, and the words beyond in `more`
/// (a byte of `u8::MAX` without an entry is a count of 255).
///
/// SAFETY: `cnt.add(cell)` is a byte of the counts that no other thread uses meanwhile.
#[inline(always)]
pub(super) unsafe fn x_count8(cnt: *mut u8, more: &mut HashMap<u32, u32>, cell: usize) {
    // SAFETY: the caller's.
    let p = unsafe { &mut *cnt.add(cell) };
    if *p == u8::MAX {
        x_count_more(more, cell);
    } else {
        *p += 1;
    }
}

/// `x_count8` beyond 255 words of a cell (rare: a cell of 4^9 cells or more that many words of a
/// batch share).
#[cold]
#[inline(never)]
fn x_count_more(more: &mut HashMap<u32, u32>, cell: usize) {
    *more.entry(cell as u32).or_insert(0) += 1;
}

/// No NCBI counterpart (layout only): EXPERIMENT (LOSAT_X_BNCSRPACK) the presence bits of the
/// blocks `b0 .. b1` of `bits` (two words per 64 cells) from the counts `cnt` (`x_count8`): word
/// `2 j` gets those of the cells `64 j ..`; returns the number of non-empty cells of the blocks.
///
/// SAFETY: `cnt` holds the blocks' cells and `bits` their words, which no other thread uses
/// meanwhile. The bits are counted in the form `pop` (LOSAT_X_BNPOPCNT).
pub(super) unsafe fn x_presence8(
    bits: *mut u64,
    cnt: *const u8,
    b0: usize,
    b1: usize,
    pop: XPop,
) -> usize {
    let mut n = 0;
    for j in b0..b1 {
        // SAFETY: the caller's; the 64 counts of block `j`.
        let w = x_nonzero_bits(unsafe { &*(cnt.add(j << 6) as *const [u8; 64]) });
        // SAFETY: the caller's.
        unsafe { *bits.add(2 * j) = w };
        n += pop.ones(w);
    }
    n
}

/// No NCBI counterpart (layout only): EXPERIMENT (LOSAT_X_BNCSRPACK) the ranks of the blocks
/// `b0 .. b1` of `bits` (presence bits written, `x_presence8`): the second word of each gets the
/// number of non-empty cells before it, from `r0`, and `starts.add(r)` the count (`x_count8`) of
/// the non-empty cell of rank `r`; returns the rank after `b1`.
///
/// SAFETY: as `x_presence8`, and `starts.add(r)` for the blocks' ranks is no other thread's.
pub(super) unsafe fn x_rank_counts8(
    bits: *mut u64,
    cnt: *const u8,
    more: &HashMap<u32, u32>,
    starts: *mut u32,
    b0: usize,
    b1: usize,
    r0: usize,
) -> usize {
    let mut r = r0;
    for j in b0..b1 {
        // SAFETY: the caller's (the words of block `j`, its counts and ranks).
        unsafe {
            *bits.add(2 * j + 1) = r as u64;
            let mut w = *bits.add(2 * j);
            while w != 0 {
                let cell = (j << 6) + w.trailing_zeros() as usize;
                w &= w - 1;
                let c = *cnt.add(cell);
                *starts.add(r) = if c == u8::MAX {
                    u32::from(c) + more.get(&(cell as u32)).copied().unwrap_or(0)
                } else {
                    u32::from(c)
                };
                r += 1;
            }
        }
    }
    r
}

/// The sum pass of a packed head (LOSAT_X_BNCSRPACK) over the blocks `b0 .. b1` (64 cells each,
/// `b0` a multiple of 32): `cnt.add(r)` holds the count of the non-empty cell of rank `r` and
/// receives the inclusive running sum from `base`; returns the running sum after `b1`. The
/// presence-vector words and counters get what `sum_cells` gives them over the same cells: the
/// presence bits of each 32 cells OR-ed into word `c >> pv_array_bts` (if it is below `pv_len`),
/// and each counter `c >> 11` the run's growth over its 2048 cells minus their non-empty cells.
/// The non-empty cells of a block are counted in the form `xpop` (LOSAT_X_BNPOPCNT).
///
/// SAFETY: as `sum_cells`, for the cells of the blocks and the ranks of their non-empty cells.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(super) unsafe fn x_sum_packed(
    pv: *mut PvArrayType,
    pv_len: usize,
    helper: *mut u32,
    bits: &[u64],
    cnt: *mut u32,
    base: u32,
    b0: usize,
    b1: usize,
    pv_array_bts: usize,
    xpop: XPop,
) -> u32 {
    debug_assert!(b0.is_multiple_of(32) && pv_array_bts >= 5);
    let (mut end, mut pop, mut h0) = (base, 0u32, base);
    for j in b0..b1 {
        let w = bits[2 * j];
        let r0 = bits[2 * j + 1] as usize;
        let k = xpop.ones(w);
        for r in r0..r0 + k {
            // SAFETY: rank `r` is of a non-empty cell of the blocks.
            unsafe {
                end = end.wrapping_add(*cnt.add(r));
                *cnt.add(r) = end;
            }
        }
        pop += k as u32;
        for (half, g) in [(w as u32, j << 6), ((w >> 32) as u32, (j << 6) + 32)] {
            let i = g >> pv_array_bts;
            if half != 0 && i < pv_len {
                // SAFETY: presence-vector word `i` holds cells of the blocks only (or the range
                // ends in it).
                unsafe { *pv.add(i) |= half };
            }
        }
        if (j + 1) % 32 == 0 || j + 1 == b1 {
            let h = j >> 5;
            // SAFETY: counter `h` holds cells of the blocks only (or the range ends in it).
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
            x_range_words(
                &seq[range_start..range_end],
                template_length,
                u64::MAX,
                || {},
                |accum, i| {
                    let pos = range_start + i;
                    // The 1-based start of the word in the concatenated query.
                    let index = (query_offset + pos + 2 - template_length) as u32;
                    let ecode1 = compute_discontiguous_index(accum, words.template_type) as usize;
                    f(ecode1, index);
                },
            );
        }
    }
}

/// No NCBI counterpart (scheduling only): a hint to bring the line of `p` into the cache; no
/// architectural access, so any address is allowed (x86_64 only).
#[inline(always)]
pub(crate) fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not access memory architecturally; any address is allowed.
    unsafe {
        core::arch::x86_64::_mm_prefetch::<{ core::arch::x86_64::_MM_HINT_T0 }>(p as *const i8);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}

/// How many words a CSR build's place pass reads ahead of the word it places (`PlaceAhead`): the
/// distance of the count passes.
pub(super) const PLACE_AHEAD: usize = 16;

/// No NCBI counterpart (memory order only): the words a CSR build's place pass has read but not yet
/// placed (`build`, `x_bn_csr::serial` and `x_mb_lookup_once`'s builds on a pool). The cursor of a
/// word (the element of `offsets` or of a packed head's `starts` that its place is taken from) is
/// prefetched where the word is read, and the word is placed `PLACE_AHEAD` words later, as the
/// count passes count, so the misses on the cursors overlap instead of coming one after another.
/// The words are placed in the order they were read, so each cursor moves as it does when every
/// word is placed where it is read, and every place is the same.
pub(super) struct PlaceAhead {
    held: [(usize, u32); PLACE_AHEAD],
    seen: usize,
}

impl PlaceAhead {
    #[inline(always)]
    pub(super) fn new() -> Self {
        Self {
            held: [(0, 0); PLACE_AHEAD],
            seen: 0,
        }
    }

    /// Takes the word whose cursor is element `i`, at `cursor` (prefetched here), and whose value
    /// is `value`; returns the word read `PLACE_AHEAD` words before it, which the caller places now.
    #[inline(always)]
    pub(super) fn push<T>(
        &mut self,
        cursor: *const T,
        i: usize,
        value: u32,
    ) -> Option<(usize, u32)> {
        prefetch(cursor);
        let slot = self.seen % PLACE_AHEAD;
        let due = (self.seen >= PLACE_AHEAD).then(|| self.held[slot]);
        self.held[slot] = (i, value);
        self.seen += 1;
        due
    }

    /// The words still held, in the order they were read, for the caller to place; the ring is
    /// empty after them.
    #[inline(always)]
    pub(super) fn drain(&mut self) -> impl Iterator<Item = (usize, u32)> + '_ {
        let end = self.seen;
        self.seen = 0;
        (end.saturating_sub(PLACE_AHEAD)..end).map(move |k| self.held[k % PLACE_AHEAD])
    }
}

/// No NCBI counterpart (scheduling only): bring the line of `s[i]` toward the core ahead of its
/// use (`Ring`): a prefetch on x86_64; on wasm32, which has no prefetch instruction, a load of
/// the element whose value is dropped (nothing when `i` is out of range), so the loads of a group
/// wait for memory together instead of one after another in the scan; elsewhere nothing.
#[inline(always)]
pub(crate) fn x_ahead<T: Copy>(s: &[T], i: usize) {
    #[cfg(target_arch = "x86_64")]
    prefetch(s.as_ptr().wrapping_add(i));
    #[cfg(target_arch = "wasm32")]
    if let Some(e) = s.get(i) {
        // SAFETY: `e` refers to an element of `s`.
        let _ = unsafe { core::ptr::read_volatile(e) };
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "wasm32")))]
    let _ = (s, i);
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
    debug_assert!(c0.is_multiple_of(2048) && pv_array_bts >= 5);
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
    Csr {
        offsets,
        x_bits: Vec::new(),
        qs,
        x_keep: None,
    }
}

/// EXPERIMENT (LOSAT_X_BNCSRPACK): the CSR of a packed head whose `starts` (`present + 1`
/// elements: the run start of each non-empty cell by rank, then the number of words), `bits` (two
/// words per 64 cells, presence and rank) and `qs` (`words` elements) have been written in full.
///
/// SAFETY: every element of `starts[..present + 1]` and `qs[..words]` is initialized.
pub(super) unsafe fn x_assemble_packed(
    mut starts: Vec<u32>,
    bits: Vec<u64>,
    mut qs: Vec<u32>,
    present: usize,
    words: usize,
) -> Csr {
    debug_assert!(starts.capacity() > present && qs.capacity() >= words);
    // SAFETY: the caller's.
    unsafe {
        starts.set_len(present + 1);
        qs.set_len(words);
    }
    Csr {
        offsets: starts,
        x_bits: bits,
        qs,
        x_keep: None,
    }
}

/// EXPERIMENT (LOSAT_X_MBLUTKEEP): `buffer` emptied, with room for `len` elements (the buffer
/// of a build's `offsets`: the batch loop's kept table, `x_mb_lut_keep::Keep::take_buffer`, or an
/// empty vector). A build writes every element it uses before it reads one.
pub(super) fn x_offsets_buffer(mut buffer: Vec<u32>, len: usize) -> Vec<u32> {
    buffer.clear();
    buffer.reserve_exact(len);
    buffer
}

/// No NCBI counterpart (counters only): a CSR build of the search.
pub(super) fn count_build(csr: &Csr, on_pool: bool) {
    BUILDS.fetch_add(1, Ordering::Relaxed);
    BUILDS_ON_POOL.fetch_add(u64::from(on_pool), Ordering::Relaxed);
    BUILT_WORDS.fetch_add(csr.words() as u64, Ordering::Relaxed);
}

/// The serial build (module comment): the count pass (each word's `offsets` element prefetched
/// and counted 16 words later), the sum pass with the presence bits and helper counters into
/// `pv_array` and `helper_array` (zeroed, of the fill's lengths), and the place pass (each word's
/// cursor prefetched and the word placed 16 words later, `PlaceAhead`).
pub(super) fn build(
    words: &DiscWords<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    x_offsets: Vec<u32>,
    x_qs: Option<Vec<u32>>,
) -> Csr {
    // Once per table: the bounds the unchecked writes below rely on.
    assert!(!words.two_templates && cells.is_multiple_of(2048) && pv_array_bts >= 5);
    assert_eq!(helper_array.len(), cells >> 11);
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): `x_offsets` is the batch loop's kept table or an empty vector.
    let mut offsets = x_offsets_buffer(x_offsets, cells + 1);
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
    // EXPERIMENT (LOSAT_X_CSRQSKEEP): `x_qs` is the batch loop's kept `qs`, taken when it has room
    // for the words, or `None`. No NCBI counterpart (memory only). The place pass below writes each
    // of the places `0 .. n_words` once (the runs partition them) before anything reads one, so
    // a kept `qs` needs no zeroing.
    let mut qs: Vec<u32> = super::x_mb_lut_keep::qs_buffer(x_qs, n_words);
    let qb = qs.as_mut_ptr();
    let place = |cell: usize, index: u32| {
        // SAFETY: `cell < cells` (the same words as the count pass), `offsets[cell]` is the end
        // of the cell's non-empty run (>= 1), so the place `e < n_words` is in the capacity.
        unsafe {
            let e = *ob.add(cell) - 1;
            *ob.add(cell) = e;
            *qb.add(e as usize) = index;
        }
    };
    // Each word's cursor `offsets[cell]` prefetched, the word placed `PLACE_AHEAD` words later.
    let mut ahead = PlaceAhead::new();
    each_word(words, |cell, index| {
        if let Some((cell, index)) = ahead.push(ob.wrapping_add(cell), cell, index) {
            place(cell, index);
        }
    });
    for (cell, index) in ahead.drain() {
        place(cell, index);
    }
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
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) struct SpinBarrier {
    n: usize,
    count: AtomicUsize,
    generation: AtomicUsize,
    poisoned: std::sync::atomic::AtomicBool,
}

/// Held by each thread of a build that waits at a `SpinBarrier`: poisons the barrier when the
/// thread unwinds.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) struct Poison<'a>(pub(super) &'a SpinBarrier);

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl Drop for Poison<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.0.poisoned.store(true, Ordering::Release);
        }
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
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

/// No NCBI counterpart (scheduling only): the words of a scan's FIFO get their lines brought in
/// in groups of this many (module comment). Where `x_ahead` is a prefetch (x86_64), which does
/// not hold the core, each word as it enters (one word per group; 8 per group cost dc-megablast
/// 1.4 % more cycles there); where it is a load (wasm32), 8 words together, so their misses are
/// waited for at once.
pub(crate) const GROUP: usize = if cfg!(target_arch = "wasm32") { 8 } else { 1 };

/// No NCBI counterpart (scheduling only): the FIFO of a CSR scan's words (module comment). A
/// word is any `Copy` value whose cell `push` gets from `cell`.
pub(crate) struct Ring<'a, T: Copy> {
    /// Slot `m & (RING - 1)` holds the `m`-th word pushed while it is one of the last `RING` of
    /// the `n` words pushed; a slot no word has reached is uninitialised (never read).
    buf: [MaybeUninit<T>; RING],
    n: usize,
    offsets: &'a [u32],
    qs: &'a [u32],
}

impl<'a, T: Copy + Default> Ring<'a, T> {
    pub(crate) fn new(csr: &'a Csr) -> Self {
        Self {
            buf: [MaybeUninit::uninit(); RING],
            n: 0,
            offsets: csr.offsets(),
            qs: csr.qs(),
        }
    }

    /// The word `w` enters, and the word `RING` places back leaves through `emit`; whether one
    /// left (none while fewer than `RING` words came before). When `w` ends a group of `GROUP`
    /// words, first the group gets its `offsets` elements brought in (`x_ahead`), and the group
    /// `RING / 2` places back the first and the last line of each of their runs (module comment).
    #[inline(always)]
    pub(crate) fn push(
        &mut self,
        w: T,
        cell: impl Fn(&T) -> usize,
        emit: &mut impl FnMut(T),
    ) -> bool {
        const MASK: usize = RING - 1;
        let n = self.n;
        if (n + 1).is_multiple_of(GROUP) {
            let (offsets, qs) = (self.offsets, self.qs);
            // (constant trip counts and bounds written from `n` alone: the loops unroll, and
            // with one word per group they are the hints of that word alone, as before groups)
            #[allow(clippy::reversed_empty_ranges)] // (none in a group of one word)
            for k in 0..GROUP - 1 {
                let m = n - (GROUP - 1) + k;
                // SAFETY: `m < n` and `m > n - GROUP >= n - RING`: word `m` was pushed, and no
                // later word wrote its slot.
                x_ahead(
                    offsets,
                    cell(unsafe { self.buf[m & MASK].assume_init_ref() }),
                );
            }
            x_ahead(offsets, cell(&w));
            if n >= RING / 2 + GROUP - 1 {
                for k in 0..GROUP {
                    let m = n - (RING / 2 + GROUP - 1) + k;
                    // SAFETY: as above (`m >= n - (RING / 2 + GROUP - 1) > n - RING`).
                    let c = cell(unsafe { self.buf[m & MASK].assume_init_ref() });
                    let a = offsets[c] as usize;
                    let b = offsets[c + 1] as usize;
                    if b > a {
                        x_ahead(qs, a);
                        if (qs.as_ptr().wrapping_add(a) as usize >> 6)
                            != (qs.as_ptr().wrapping_add(b - 1) as usize >> 6)
                        {
                            x_ahead(qs, b - 1);
                        }
                    }
                }
            }
        }
        let slot = n & MASK;
        let leaves = n >= RING;
        if leaves {
            // SAFETY: `n >= RING`, so word `n - RING`, whose slot this is, was pushed and
            // written; the `RING - 1` words pushed since then wrote other slots.
            emit(unsafe { self.buf[slot].assume_init() });
        }
        self.buf[slot] = MaybeUninit::new(w);
        self.n = n + 1;
        leaves
    }

    /// The words still in the FIFO leave through `emit`, in order.
    pub(crate) fn flush(&mut self, emit: &mut impl FnMut(T)) {
        for m in self.n.saturating_sub(RING)..self.n {
            // SAFETY: `m < n` and `m >= n - RING`: word `m` was pushed, and no later word
            // wrote its slot.
            emit(unsafe { self.buf[m & (RING - 1)].assume_init() });
        }
        self.n = 0;
    }
}

impl<T: Copy + Default> XRing<T> for Ring<'_, T> {
    #[inline(always)]
    fn push(&mut self, w: T, cell: impl Fn(&T) -> usize, emit: &mut impl FnMut(T, usize)) -> bool {
        Ring::push(self, w, &cell, &mut |w| emit(w, cell(&w)))
    }

    #[inline(always)]
    fn flush(&mut self, cell: impl Fn(&T) -> usize, emit: &mut impl FnMut(T, usize)) {
        Ring::flush(self, &mut |w| emit(w, cell(&w)))
    }
}

/// No NCBI counterpart (scheduling only): EXPERIMENT (LOSAT_X_BNCSRPACK) `Ring` over a packed head
/// (`XPacked`), with the rank of each word in `buf`.
pub(crate) struct XRingPacked<'a, T: Copy> {
    /// `Ring::buf`, and the rank of each word in the same slot (uninitialised until written).
    buf: [MaybeUninit<T>; RING],
    ranks: [MaybeUninit<u32>; RING],
    n: usize,
    head: XPacked<'a>,
}

impl<T: Copy + Default> XRing<T> for XRingPacked<'_, T> {
    /// `Ring::push` on a packed head. A word of an empty cell does not enter (its run is empty:
    /// `emit` would read no offset of it); the others enter with their rank, and the groups get
    /// the head's elements of their runs and the first and the last line of their runs brought
    /// in as on a head of one element per cell.
    #[inline(always)]
    fn push(&mut self, w: T, cell: impl Fn(&T) -> usize, emit: &mut impl FnMut(T, usize)) -> bool {
        const MASK: usize = RING - 1;
        let Some(r) = x_rank(self.head.bits, cell(&w), self.head.pop) else {
            return false;
        };
        let n = self.n;
        if (n + 1).is_multiple_of(GROUP) {
            let (offsets, qs) = (self.head.offsets, self.head.qs);
            #[allow(clippy::reversed_empty_ranges)] // (none in a group of one word)
            for k in 0..GROUP - 1 {
                let m = n - (GROUP - 1) + k;
                // SAFETY: `m < n` and `m > n - GROUP >= n - RING`: word `m` was pushed with its
                // rank, and no later word wrote its slot.
                x_ahead(offsets, unsafe { self.ranks[m & MASK].assume_init() }
                    as usize);
            }
            x_ahead(offsets, r);
            if n >= RING / 2 + GROUP - 1 {
                for j in 0..GROUP {
                    let m = n - (RING / 2 + GROUP - 1) + j;
                    // SAFETY: as above (`m >= n - (RING / 2 + GROUP - 1) > n - RING`).
                    let k = unsafe { self.ranks[m & MASK].assume_init() } as usize;
                    let a = offsets[k] as usize;
                    let b = offsets[k + 1] as usize;
                    let e = b.max(a + 1) - 1;
                    x_ahead(qs, a);
                    if (qs.as_ptr().wrapping_add(a) as usize >> 6)
                        != (qs.as_ptr().wrapping_add(e) as usize >> 6)
                    {
                        x_ahead(qs, e);
                    }
                }
            }
        }
        let slot = n & MASK;
        let leaves = n >= RING;
        if leaves {
            // SAFETY: `n >= RING`, so word `n - RING`, whose slot this is, was pushed and
            // written with its rank; the `RING - 1` words pushed since then wrote other slots.
            unsafe {
                emit(
                    self.buf[slot].assume_init(),
                    self.ranks[slot].assume_init() as usize,
                )
            };
        }
        self.buf[slot] = MaybeUninit::new(w);
        self.ranks[slot] = MaybeUninit::new(r as u32);
        self.n = n + 1;
        leaves
    }

    /// `Ring::flush`.
    fn flush(&mut self, _cell: impl Fn(&T) -> usize, emit: &mut impl FnMut(T, usize)) {
        for m in self.n.saturating_sub(RING)..self.n {
            // SAFETY: `m < n` and `m >= n - RING`: word `m` was pushed with its rank, and no
            // later word wrote its slot.
            unsafe {
                emit(
                    self.buf[m & (RING - 1)].assume_init(),
                    self.ranks[m & (RING - 1)].assume_init() as usize,
                )
            };
        }
        self.n = 0;
    }
}

const _: () = assert!(
    RING.is_power_of_two() && RING >= 2 * GROUP && RING.is_multiple_of(GROUP) && GROUP >= 1
);
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
        let serial = build(
            words,
            cells,
            &mut pv2,
            pv_array_bts,
            &mut helper2,
            Vec::new(),
            None,
        );
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
#[path = "x_dc_csr_tests.rs"]
mod tests;
