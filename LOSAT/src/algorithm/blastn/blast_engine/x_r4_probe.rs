//! EXPERIMENT (exp/r4-scan-probes, v4-R R4; evidence only, never merged): probes of the lookup
//! scan of blastn (reference scan of `search_query_batch`) and dc-megablast
//! (`x_pair_seed::disc::scan_strip_disc`), chosen once from `R4_PROBE` (comma list):
//! - `pf=D` (D a power of two, at most `RING`): the scan's words are processed `D` words late
//!   through a FIFO ring; a word entering the ring prefetches its backbone cell
//!   (`hashtable[index]`). Same words, same order, same values.
//! - `st2`: with `pf`, the word `D / 2` places back in the ring reads its backbone cell and
//!   prefetches the first chain element (`next_pos[q_off]`).
//! - `st3`: `st2`, and the word `D / 4` places back reads its first chain element and prefetches
//!   the second.
//! - `nopf` (control): with `pf`, the ring only, no prefetch and no read ahead.
//! - `pv`: a cell-level presence bit vector (one bit per cell, built from `hashtable` at every
//!   call, its cost included) drops the words of empty cells before the ring; an empty cell has
//!   no chain, so the reference does nothing for it.
//! - `ctx` (dc only): the context of a seed from a table of the first run of every 256 offsets
//!   plus a forward walk, the value of `QueryContextIndex::x_direct`.
//! - `csr` (exp/r11-csr, v4-R R11): the scans read the chains from their CSR
//!   (`lookup::x_r11_csr`, built after the fill): `offsets[cell] .. offsets[cell + 1]` of `qs`,
//!   in the chain's order, after the same presence test. With `pf`, a word entering the ring
//!   prefetches `offsets[cell]`, and with `st2` the word `D / 2` places back reads its two
//!   offsets and prefetches the first and the last line of its run of `qs`.
//! - `csrb`: the CSR is built (its cost and memory), the scans walk the chains.
//! Unset or empty = the dev code paths.
use super::QueryContextIndex;
use std::sync::OnceLock;

pub(super) const RING: usize = 64;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Probe {
    pub(super) pf: usize,
    /// 0 `nopf`, 1 backbone prefetch (default with `pf`), 2 `st2`, 3 `st3`.
    pub(super) lvl: u8,
    pub(super) pv: bool,
    pub(super) ctx: bool,
    /// R11: the scans read the CSR.
    pub(super) csr: bool,
}

impl Probe {
    #[inline(always)]
    pub(super) fn on(&self) -> bool {
        self.pf > 0 || self.pv || self.ctx || self.csr
    }
}

pub(super) fn probe() -> Probe {
    static P: OnceLock<Probe> = OnceLock::new();
    *P.get_or_init(|| {
        let mut p = Probe {
            lvl: 1,
            ..Probe::default()
        };
        if let Ok(v) = std::env::var("R4_PROBE") {
            for t in v.split(',').map(str::trim) {
                match t {
                    "" => {}
                    "st2" => p.lvl = p.lvl.max(2),
                    "st3" => p.lvl = 3,
                    "nopf" => p.lvl = 0,
                    "pv" => p.pv = true,
                    "ctx" => p.ctx = true,
                    "csr" => p.csr = true,
                    // Read by `lookup::x_r11_csr` (how the CSR is built).
                    "csrb" | "csrw" | "csrchk" | "csrp" | "pvbf" => {}
                    // R16 (exp/r16-csr): the words build's arms, read there too.
                    "cpf" | "psplit" | "qpf" | "cblk" | "sum1" | "sumlean" | "sum0" | "rdx"
                    | "c4" | "c8" | "offkeep" | "cpfc" | "bar" => {}
                    _ if t.starts_with("pf=") => {
                        let d: usize = t[3..].parse().expect("R4_PROBE pf=D");
                        assert!(d.is_power_of_two() && d >= 2 && d <= RING, "R4_PROBE pf=D");
                        p.pf = d;
                    }
                    _ => panic!("R4_PROBE: unknown item {t}"),
                }
            }
            assert!(
                p.lvl == 1 || p.pf >= 4,
                "R4_PROBE st2 / st3 / nopf need pf=D >= 4"
            );
            assert!(
                !(p.csr && (p.pv || p.lvl == 3)),
                "R4_PROBE csr: no pv / st3"
            );
            eprintln!("[R4_PROBE] {p:?}");
        }
        p
    })
}

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

/// One bit per backbone cell: set iff the cell's chain is not empty.
pub(super) struct CellPv {
    bits: Vec<u64>,
}

impl CellPv {
    pub(super) fn build(hashtable: &[u32]) -> Self {
        let mut bits = vec![0u64; hashtable.len().div_ceil(64)];
        for (w, chunk) in bits.iter_mut().zip(hashtable.chunks(64)) {
            let mut b = 0u64;
            for (i, &v) in chunk.iter().enumerate() {
                b |= u64::from(v != 0) << i;
            }
            *w = b;
        }
        Self { bits }
    }

    #[inline(always)]
    pub(super) fn test(&self, cell: usize) -> bool {
        debug_assert!(cell >> 6 < self.bits.len());
        // SAFETY: `cell` is a backbone index (< hashtable.len()), so `cell >> 6 < bits.len()`.
        unsafe { (*self.bits.get_unchecked(cell >> 6) >> (cell & 63)) & 1 != 0 }
    }
}

/// The index of the first run of `x_ranges` starting after offset `b << 8`, for every block `b`.
pub(super) struct CtxTab {
    after: Vec<u32>,
}

impl CtxTab {
    pub(super) fn build(q: &QueryContextIndex) -> Self {
        let blocks = q.x_total_len.div_ceil(256) + 1;
        let mut after = Vec::with_capacity(blocks);
        let mut a = 0usize;
        for b in 0..blocks {
            let start = b << 8;
            while a < q.x_ranges.len() && q.x_ranges[a].0 <= start {
                a += 1;
            }
            after.push(a as u32);
        }
        Self { after }
    }

    /// `q.context_for_offset(n)`.
    #[inline(always)]
    pub(super) fn context(&self, q: &QueryContextIndex, n: usize) -> usize {
        if n >= q.x_total_len || !q.direct_map.is_empty() {
            return q.context_for_offset(n);
        }
        let r = &q.x_ranges;
        let mut after = self.after[n >> 8] as usize;
        while after < r.len() && r[after].0 <= n {
            after += 1;
        }
        if after == 0 {
            return 0;
        }
        let (_, end, context) = r[after - 1];
        if n < end {
            context as usize
        } else {
            0
        }
    }
}

/// One step of the FIFO of the scan's words: the word `w` enters the ring (with the prefetches
/// of level `lvl`), and the word `d` places back leaves it through `emit`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn ring_step<T: Copy>(
    buf: &mut [T; RING],
    n: &mut usize,
    d: usize,
    lvl: u8,
    hashtable: &[u32],
    next_pos: &[u32],
    cell: impl Fn(&T) -> usize,
    w: T,
    emit: &mut impl FnMut(T),
) {
    let mask = d - 1;
    if lvl >= 1 {
        prefetch(hashtable.as_ptr().wrapping_add(cell(&w)));
    }
    if lvl >= 2 && *n >= d / 2 {
        let v = hashtable[cell(&buf[(*n - d / 2) & mask])] as usize;
        if v != 0 {
            prefetch(next_pos.as_ptr().wrapping_add(v));
        }
    }
    if lvl >= 3 && *n >= d / 4 {
        let v1 = hashtable[cell(&buf[(*n - d / 4) & mask])] as usize;
        if v1 != 0 {
            let v2 = next_pos[v1] as usize;
            if v2 != 0 {
                prefetch(next_pos.as_ptr().wrapping_add(v2));
            }
        }
    }
    let slot = *n & mask;
    if *n >= d {
        emit(buf[slot]);
    }
    buf[slot] = w;
    *n += 1;
}

/// The words still in the ring, in order.
#[inline(always)]
pub(super) fn ring_flush<T: Copy>(buf: &[T; RING], n: usize, d: usize, emit: &mut impl FnMut(T)) {
    for m in n.saturating_sub(d)..n {
        emit(buf[m & (d - 1)]);
    }
}

/// The FIFO of the blastn reference scan's words `(kmer_start, s_range, kmer)`.
pub(super) struct BnRing<'a, const PV: bool> {
    buf: [(usize, usize, u64); RING],
    n: usize,
    d: usize,
    lvl: u8,
    hashtable: &'a [u32],
    next_pos: &'a [u32],
    pv: Option<&'a CellPv>,
}

impl<'a, const PV: bool> BnRing<'a, PV> {
    pub(super) fn new(
        p: Probe,
        hashtable: &'a [u32],
        next_pos: &'a [u32],
        pv: Option<&'a CellPv>,
    ) -> Self {
        Self {
            buf: [(0, 0, 0); RING],
            n: 0,
            d: p.pf,
            lvl: p.lvl,
            hashtable,
            next_pos,
            pv,
        }
    }

    #[inline(always)]
    pub(super) fn push(
        &mut self,
        w: (usize, usize, u64),
        emit: &mut impl FnMut(usize, usize, u64),
    ) {
        if PV {
            // PV instances are built with `Some`.
            if !self.pv.unwrap().test(w.2 as usize) {
                return;
            }
        }
        if self.d == 0 {
            emit(w.0, w.1, w.2);
            return;
        }
        ring_step(
            &mut self.buf,
            &mut self.n,
            self.d,
            self.lvl,
            self.hashtable,
            self.next_pos,
            |e: &(usize, usize, u64)| e.2 as usize,
            w,
            &mut |e: (usize, usize, u64)| emit(e.0, e.1, e.2),
        );
    }

    pub(super) fn flush(&mut self, emit: &mut impl FnMut(usize, usize, u64)) {
        if self.d == 0 {
            return;
        }
        ring_flush(&self.buf, self.n, self.d, &mut |e: (usize, usize, u64)| {
            emit(e.0, e.1, e.2)
        });
        self.n = 0;
    }
}

/// EXPERIMENT (exp/r11-csr, v4-R R11): one step of the FIFO of the scan's words over the CSR:
/// the word `w` enters the ring (level 1: prefetch of `offsets[cell]`), with level 2 the word
/// `d / 2` places back reads its two offsets and prefetches the first and the last line of its
/// run of `qs`, and the word `d` places back leaves through `emit`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub(super) fn ring_step_csr<T: Copy>(
    buf: &mut [T; RING],
    n: &mut usize,
    d: usize,
    lvl: u8,
    offsets: &[u32],
    qs: &[u32],
    cell: impl Fn(&T) -> usize,
    w: T,
    emit: &mut impl FnMut(T),
) {
    let mask = d - 1;
    if lvl >= 1 {
        prefetch(offsets.as_ptr().wrapping_add(cell(&w)));
    }
    if lvl >= 2 && *n >= d / 2 {
        let c = cell(&buf[(*n - d / 2) & mask]);
        let a = offsets[c] as usize;
        let b = offsets[c + 1] as usize;
        if b > a {
            let pa = qs.as_ptr().wrapping_add(a);
            let pb = qs.as_ptr().wrapping_add(b - 1);
            prefetch(pa);
            if (pa as usize >> 6) != (pb as usize >> 6) {
                prefetch(pb);
            }
        }
    }
    let slot = *n & mask;
    if *n >= d {
        emit(buf[slot]);
    }
    buf[slot] = w;
    *n += 1;
}

/// EXPERIMENT (exp/r11-csr): the FIFO of the reference scans' words `(kmer_start, s_range,
/// index)` over the CSR (`ring_step_csr`); `d == 0` = no ring.
pub(super) struct CsrRing<'a> {
    buf: [(usize, usize, u64); RING],
    n: usize,
    d: usize,
    lvl: u8,
    offsets: &'a [u32],
    qs: &'a [u32],
}

impl<'a> CsrRing<'a> {
    pub(super) fn new(p: Probe, offsets: &'a [u32], qs: &'a [u32]) -> Self {
        Self {
            buf: [(0, 0, 0); RING],
            n: 0,
            d: p.pf,
            lvl: p.lvl,
            offsets,
            qs,
        }
    }

    #[inline(always)]
    pub(super) fn push(
        &mut self,
        w: (usize, usize, u64),
        emit: &mut impl FnMut(usize, usize, u64),
    ) {
        if self.d == 0 {
            emit(w.0, w.1, w.2);
            return;
        }
        ring_step_csr(
            &mut self.buf,
            &mut self.n,
            self.d,
            self.lvl,
            self.offsets,
            self.qs,
            |e: &(usize, usize, u64)| e.2 as usize,
            w,
            &mut |e: (usize, usize, u64)| emit(e.0, e.1, e.2),
        );
    }

    pub(super) fn flush(&mut self, emit: &mut impl FnMut(usize, usize, u64)) {
        if self.d == 0 {
            return;
        }
        ring_flush(&self.buf, self.n, self.d, &mut |e: (usize, usize, u64)| {
            emit(e.0, e.1, e.2)
        });
        self.n = 0;
    }
}
