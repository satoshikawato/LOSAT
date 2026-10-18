//! EXPERIMENT (exp/r4-scan-probes, v4-R R4; evidence only, never merged): probes of the lookup
//! scan of blastn (reference scan of `search_query_batch`) and dc-megablast
//! (`x_pair_seed::disc::scan_strip_disc`), chosen once from `R4_PROBE` (comma list):
//! - `pf=D` (D a power of two, at most `RING`): the scan's words are processed `D` words late
//!   through a FIFO ring; a word entering the ring prefetches its backbone cell
//!   (`hashtable[index]`). Same words, same order, same values.
//! - `st2`: with `pf`, the word `D / 2` places back in the ring reads its backbone cell and
//!   prefetches the first chain element (`next_pos[q_off]`).
//! - `pv`: a cell-level presence bit vector (one bit per cell, built from `hashtable` at every
//!   call, its cost included) drops the words of empty cells before the ring; an empty cell has
//!   no chain, so the reference does nothing for it.
//! - `ctx` (dc only): the context of a seed from a table of the first run of every 256 offsets
//!   plus a forward walk, the value of `QueryContextIndex::x_direct`.
//! Unset or empty = the dev code paths.
use super::QueryContextIndex;
use std::sync::OnceLock;

pub(super) const RING: usize = 64;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Probe {
    pub(super) pf: usize,
    pub(super) st2: bool,
    pub(super) pv: bool,
    pub(super) ctx: bool,
}

impl Probe {
    #[inline(always)]
    pub(super) fn on(&self) -> bool {
        self.pf > 0 || self.pv || self.ctx
    }
}

pub(super) fn probe() -> Probe {
    static P: OnceLock<Probe> = OnceLock::new();
    *P.get_or_init(|| {
        let mut p = Probe::default();
        if let Ok(v) = std::env::var("R4_PROBE") {
            for t in v.split(',').map(str::trim) {
                match t {
                    "" => {}
                    "st2" => p.st2 = true,
                    "pv" => p.pv = true,
                    "ctx" => p.ctx = true,
                    _ if t.starts_with("pf=") => {
                        let d: usize = t[3..].parse().expect("R4_PROBE pf=D");
                        assert!(d.is_power_of_two() && d >= 2 && d <= RING, "R4_PROBE pf=D");
                        p.pf = d;
                    }
                    _ => panic!("R4_PROBE: unknown item {t}"),
                }
            }
            assert!(!p.st2 || p.pf > 0, "R4_PROBE st2 needs pf=D");
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

/// The FIFO of the blastn reference scan's words `(kmer_start, s_range, kmer)`.
pub(super) struct BnRing<'a, const ST2: bool, const PV: bool> {
    buf: [(usize, usize, u64); RING],
    n: usize,
    d: usize,
    hashtable: &'a [u32],
    next_pos: &'a [u32],
    pv: Option<&'a CellPv>,
}

impl<'a, const ST2: bool, const PV: bool> BnRing<'a, ST2, PV> {
    pub(super) fn new(
        d: usize,
        hashtable: &'a [u32],
        next_pos: &'a [u32],
        pv: Option<&'a CellPv>,
    ) -> Self {
        Self {
            buf: [(0, 0, 0); RING],
            n: 0,
            d,
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
            // SAFETY of unwrap: PV instances are built with `Some`.
            if !self.pv.unwrap().test(w.2 as usize) {
                return;
            }
        }
        let d = self.d;
        if d == 0 {
            emit(w.0, w.1, w.2);
            return;
        }
        let mask = d - 1;
        prefetch(self.hashtable.as_ptr().wrapping_add(w.2 as usize));
        if ST2 && self.n >= d / 2 {
            let mid = self.buf[(self.n - d / 2) & mask].2 as usize;
            let v = self.hashtable[mid] as usize;
            if v != 0 {
                prefetch(self.next_pos.as_ptr().wrapping_add(v));
            }
        }
        let slot = self.n & mask;
        if self.n >= d {
            let e = self.buf[slot];
            emit(e.0, e.1, e.2);
        }
        self.buf[slot] = w;
        self.n += 1;
    }

    pub(super) fn flush(&mut self, emit: &mut impl FnMut(usize, usize, u64)) {
        let d = self.d;
        if d == 0 {
            return;
        }
        for m in self.n.saturating_sub(d)..self.n {
            let e = self.buf[m & (d - 1)];
            emit(e.0, e.1, e.2);
        }
        self.n = 0;
    }
}
