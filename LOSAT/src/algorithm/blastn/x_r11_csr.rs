//! EXPERIMENT (exp/r11-csr, v4-R R11; evidence only, never merged): the megablast lookup's chains
//! laid out as CSR (compressed sparse rows): `offsets[cell] .. offsets[cell + 1]` indexes `qs`,
//! which holds the cell's 1-based query offsets in the order the chain lists them
//! (`hashtable[cell]`, `next_pos[...]`, ... up to 0). Chosen from `R4_PROBE` (comma list):
//! - `csr`: build the CSR after the fill and let the probed scans read it instead of the chains.
//! - `csrb`: build it only (its cost and resident memory), the scans walk the chains.
//! - `csrw` (with `csr` or `csrb`; dc tables and, since R15, the contiguous megablast / blastn
//!   tables): build it from the fill's words instead of the chains: one pass of the word loop
//!   counts the words of each cell, a prefix sum turns the counts into the ends of the runs, and
//!   a second pass of the word loop places each word at `--end[cell]`, so a run lists its words
//!   in the reverse of the fill's order, which is the chain's order. The blastn task's contiguous
//!   tables (`ascending_cells`) list each cell's words in the fill's order: the counts go to
//!   `offsets[cell + 2]`, the prefix sum makes `offsets[cell + 1]` the start of the cell, and the
//!   place pass goes forward (`p = offsets[cell + 1]; offsets[cell + 1] = p + 1; qs[p] = q`), so
//!   after it `offsets[cell]` is the start of the cell (no extra pass over the cells).
//!   (`offsets[cells + 1]` is then `qs[0]`; it is not read during the place pass.)
//!   A CSR that replaces the chains must also make the presence vector and the helper counters
//!   (`longest_chain`): the sum pass computes them from the counts into scratch arrays of the
//!   fill's lengths (the bit of every non-empty cell; `n - 1` for the cell's counter), so the
//!   cost is faithful; they are dropped at the end, and with `csrchk` compared with the fill's.
//! - `csrp` (with `csrw`; dc tables, one template, a search pool of two threads or more): the
//!   words build on the pool (`x_disc_lookup_par::x_r15_csr_par`): the fill's own readers
//!   (`Fill::read::<false>`) feed the owners of the cell blocks, once to count and once to
//!   place; otherwise the serial `attach_words`.
//! - `csrchk` (with `csrw`): also build it from the chains and compare (panics on a difference).
//! The CSR is kept in the spare capacity of `next_pos` (after `next_pos.len()`), so the table
//! structs keep their layout and the default code paths their machine code; nothing else
//! changes `next_pos`'s capacity (the tables are never cloned or grown).
//! Layout there: `offsets` (`hashtable.len() + 1` words), then `qs` (`offsets[last]` words).
use std::sync::OnceLock;
use std::time::Instant;

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Mode {
    /// Build the CSR at every megablast table build.
    pub(crate) build: bool,
    /// The probed scans read it.
    pub(crate) scan: bool,
    /// Build it from the words (dc tables; contiguous tables too since R15).
    pub(crate) words: bool,
    /// Build it on the search pool (dc tables).
    pub(crate) par: bool,
    /// Compare the words build with the chain walk.
    pub(crate) check: bool,
    /// R15 `pvbf`: the sum pass makes the PV bits and helper counters without a branch per cell.
    pub(crate) pvbf: bool,
    /// R16 `cpf`: the count and place passes prefetch the cell's `offsets` element 16 words ahead
    /// (the serial build through a ring of 16 words; `csrp` in the owners' bucket loops, as the
    /// fill's `update_all`).
    pub(crate) cpf: bool,
    /// R16 `psplit` (`csrp`): the owner's place of a bucket in two loops, the `offsets` decrements
    /// (each entry's `value` becomes its place) and then the `qs` writes.
    pub(crate) psplit: bool,
    /// R16 `qpf` (`csrp`, with `psplit`): the `qs` loop prefetches for writing (PREFETCHW) the
    /// place 16 entries ahead.
    pub(crate) qpf: bool,
    /// R16 `cblk` (`csrp`): the owners hold contiguous runs of blocks (block `b` to owner
    /// `b * parts / blocks`) instead of the fill's round robin.
    pub(crate) cblk: bool,
    /// R16 `sum1` (`csrp`, with `cblk`): one sum pass per owner from its base (the owners' word
    /// counts are the lengths of their buckets), instead of block sums, a prefix over the blocks
    /// and a second pass.
    pub(crate) sum1: bool,
    /// R16 `sumlean`: the sum pass makes each PV word from 32 cells at a time and each helper
    /// counter from the difference of the running sums minus the non-empty cells (popcount).
    pub(crate) sumlean: bool,
    /// R16 `sum0` (measurement only): the sum pass without the PV bits and helper counters (the
    /// plain prefix sum of R11); `csrchk` then skips the PV / helper check.
    pub(crate) sum0: bool,
    /// R16 `rdx` (`csrp`, measurement only): before each chunk's place read, the same read without
    /// the bucket pushes (`rdx_ns`).
    pub(crate) rdx: bool,
    /// R16 `c4` (`csrp`): the count pass's read pushes the cell codes only (4 bytes) into buckets
    /// of its own; `c8`: the same probe read pushing the fill's entries (control of `c4`).
    pub(crate) c4: bool,
    pub(crate) c8: bool,
    /// R16 `offkeep` (serial build, measurement only): the counts and running sums in a buffer the
    /// thread keeps from table to table (zeroed by a write over resident pages), copied into
    /// `next_pos`'s capacity after the place pass (`copy_ns`, not in `build_ns`).
    pub(crate) offkeep: bool,
    /// R16 `cpfc`: `cpf` in the count pass only.
    pub(crate) cpfc: bool,
    /// R16 `bar` (`csrp`, with `cblk,sum1,sumlean`): the whole build of a table as one
    /// `rayon::broadcast` whose threads pass a spin barrier between the phases, instead of one
    /// `run_jobs` per phase and chunk.
    pub(crate) bar: bool,
}

pub(crate) fn mode() -> Mode {
    static M: OnceLock<Mode> = OnceLock::new();
    *M.get_or_init(|| {
        let mut m = Mode::default();
        if let Ok(v) = std::env::var("R4_PROBE") {
            for t in v.split(',').map(str::trim) {
                match t {
                    "csr" => {
                        m.build = true;
                        m.scan = true;
                    }
                    "csrb" => m.build = true,
                    "csrw" => m.words = true,
                    "csrp" => m.par = true,
                    "csrchk" => m.check = true,
                    "pvbf" => m.pvbf = true,
                    "cpf" => m.cpf = true,
                    "psplit" => m.psplit = true,
                    "qpf" => m.qpf = true,
                    "cblk" => m.cblk = true,
                    "sum1" => m.sum1 = true,
                    "sumlean" => m.sumlean = true,
                    "sum0" => m.sum0 = true,
                    "rdx" => m.rdx = true,
                    "c4" => m.c4 = true,
                    "c8" => m.c8 = true,
                    "offkeep" => m.offkeep = true,
                    "cpfc" => m.cpfc = true,
                    "bar" => m.bar = true,
                    _ => {}
                }
            }
        }
        m
    })
}

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

/// R16 `qpf`: PREFETCHW (the line in the exclusive state, ready for a store; Broadwell and later,
/// the workstation's Ice Lake included). The default target has no `prfchw` feature, under which
/// `_MM_HINT_ET0` is a plain PREFETCHT0, hence the instruction itself.
#[inline(always)]
pub(crate) fn prefetch_w<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not access memory architecturally; any address is allowed.
    unsafe {
        core::arch::asm!("prefetchw [{0}]", in(reg) p, options(nostack, preserves_flags, readonly));
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
}

/// R16 `sumlean` (and `sum0` with `pvh == false`): cells `c0 .. c1` (`c0` a multiple of 2048, and
/// `bts >= 5`): `cnt.add(c)` holds cell `c`'s count and receives its inclusive running sum from
/// `base`; returns the running sum after `c1`. With `pvh`, each 32 cells give 32 presence bits
/// (OR-ed into word `c >> bts` when the word's last group or the range ends, if the word is below
/// `pv_len`), and each counter `c >> 11` gets the run's growth over its 2048 cells minus their
/// non-empty cells (`sum(cnt - (cnt != 0))`, R15 `pvbf`'s value, wrapping), so the loop has no
/// per-cell work beyond the running sum and one bit.
///
/// SAFETY: the counts, PV words and counters of the range are the caller's alone, and the
/// pointers are valid for them.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) unsafe fn range_lean(
    pv: *mut u32,
    pv_len: usize,
    helper: *mut u32,
    cnt: *mut u32,
    base: u32,
    c0: usize,
    c1: usize,
    bts: usize,
    pvh: bool,
) -> u32 {
    debug_assert!(c0 % 2048 == 0 && bts >= 5);
    let pv_last = (1usize << bts) - 1;
    let (mut end, mut word, mut pop, mut h0) = (base, 0u32, 0u32, base);
    let mut g = c0;
    while g < c1 {
        let n = (c1 - g).min(32);
        let mut w32 = 0u32;
        if n == 32 {
            for i in 0..32 {
                // SAFETY: cell `g + i < c1` is the caller's.
                unsafe {
                    let x = *cnt.add(g + i);
                    w32 |= ((x != 0) as u32) << i;
                    end = end.wrapping_add(x);
                    *cnt.add(g + i) = end;
                }
            }
        } else {
            for i in 0..n {
                // SAFETY: as above.
                unsafe {
                    let x = *cnt.add(g + i);
                    w32 |= ((x != 0) as u32) << i;
                    end = end.wrapping_add(x);
                    *cnt.add(g + i) = end;
                }
            }
        }
        g += n;
        if pvh {
            word |= w32;
            pop += w32.count_ones();
            if g & pv_last == 0 || g >= c1 {
                let i = (g - 1) >> bts;
                if i < pv_len {
                    // SAFETY: PV word `i` holds cells of the range only (or the range ends in it).
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
    }
    end
}

/// R16 `offkeep`: the thread's kept buffer of counts (measurement only).
thread_local! {
    static XR16_KEPT: std::cell::RefCell<Vec<u32>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Builds the CSR of the chains of `hashtable` / `next_pos` into the spare capacity of
/// `next_pos` (one pass over the cells in order; each cell's chain walked from its head, the
/// offsets written as they come; the first chain element of the cell `PF` cells ahead is
/// prefetched). Prints one line `[R11] csr ...` to stderr.
pub(crate) fn attach(hashtable: &[u32], next_pos: &mut Vec<u32>) {
    const PF: usize = 16;
    let t0 = Instant::now();
    let cells = hashtable.len();
    let n = next_pos.len();
    // At most one word per query offset (`next_pos[0]` is never a word).
    let need = cells + 1 + n;
    next_pos.reserve_exact(need);
    let t1 = Instant::now();
    let mut words = 0usize;
    {
        let (np, spare) = {
            let ptr = next_pos.as_mut_ptr();
            // SAFETY: `ptr .. ptr + n` is the initialized vector and `ptr + n .. ptr + n + need`
            // its reserved capacity; the two slices do not overlap, and the capacity part is
            // written before it is read (`offsets` in cell order, `qs` from the front).
            unsafe {
                (
                    std::slice::from_raw_parts(ptr as *const u32, n),
                    std::slice::from_raw_parts_mut(
                        ptr.add(n) as *mut std::mem::MaybeUninit<u32>,
                        need,
                    ),
                )
            }
        };
        let (offsets, qs) = spare.split_at_mut(cells + 1);
        for c in 0..cells {
            if c + PF < cells {
                let h = hashtable[c + PF] as usize;
                if h != 0 {
                    prefetch(np.as_ptr().wrapping_add(h));
                }
            }
            offsets[c].write(words as u32);
            let mut q = hashtable[c];
            while q != 0 {
                qs[words].write(q);
                words += 1;
                q = np[q as usize];
            }
        }
        offsets[cells].write(words as u32);
    }
    let t2 = Instant::now();
    // EXPERIMENT (exp/r14-wall): the chain-walk CSR build's interval (`R13_PROBE=tl`; the
    // contiguous blastn tables are built here).
    crate::utils::xtiming::x_r13("csr", Some(t0), words as u64);
    eprintln!(
        "[R11] csr cells={cells} words={words} next_pos_len={n} reserve_ns={} build_ns={} csr_bytes={}",
        (t1 - t0).as_nanos(),
        (t2 - t1).as_nanos(),
        (cells + 1 + words) * 4,
    );
}

/// The CSR `attach` left in `next_pos`'s capacity: `(offsets, qs)`.
#[inline(always)]
pub(crate) fn view<'a>(hashtable: &[u32], next_pos: &'a Vec<u32>) -> (&'a [u32], &'a [u32]) {
    let cells = hashtable.len();
    let n = next_pos.len();
    assert!(
        next_pos.capacity() >= n + cells + 1,
        "R11: no CSR in this table (R4_PROBE csr builds it)"
    );
    // SAFETY: `attach` wrote `offsets` (cells + 1 words) and then `qs` (offsets[cells] words)
    // right after the vector's `n` elements, inside its capacity; the table is not changed
    // after its build.
    unsafe {
        let base = next_pos.as_ptr().add(n);
        let offsets = std::slice::from_raw_parts(base, cells + 1);
        let words = offsets[cells] as usize;
        assert!(n + cells + 1 + words <= next_pos.capacity());
        let qs = std::slice::from_raw_parts(base.add(cells + 1), words);
        (offsets, qs)
    }
}

/// The words of a table's fill, in the fill's order: `(cell, 1-based query offset)`.
pub(crate) trait Words {
    fn each(&self, f: &mut impl FnMut(usize, u32));
    /// The chains list each cell's words in the fill's order (the blastn task's tables).
    fn ascending(&self) -> bool {
        false
    }
}

/// The fill's presence vector, its shift and its helper counters, for the scratch arrays of the
/// words builds (`pvh_scratch`) to be made and checked against.
#[derive(Clone, Copy)]
pub(crate) struct PvHelper<'a> {
    pub(crate) pv: &'a [u32],
    pub(crate) bts: usize,
    pub(crate) helper: &'a [u32],
}

/// The scratch arrays of a words build: fresh zeroed arrays of the fill's lengths.
pub(crate) struct PvScratch {
    pub(crate) pv: Vec<u32>,
    pub(crate) helper: Vec<u32>,
}

impl PvScratch {
    pub(crate) fn new(ph: &PvHelper<'_>) -> Self {
        Self {
            pv: vec![0u32; ph.pv.len()],
            helper: vec![0u32; ph.helper.len()],
        }
    }

    /// The cell `c` holds `n > 0` words: its presence bit (as `pv_set_shift`) and `n - 1` more
    /// words than the first on its counter (`c >> 11`, 2048 cells per counter).
    #[inline(always)]
    pub(crate) fn cell(&mut self, c: usize, n: u32, bts: usize) {
        let i = c >> bts;
        if i < self.pv.len() {
            self.pv[i] |= 1u32 << (c & 31);
        }
        self.helper[c >> 11] += n - 1;
    }

    /// R15 `pvbf`: cells `c0 .. c1` (`c0` a multiple of `2^bts` and of 2048): `at(c)` gives the cell's count and where to
    /// write its inclusive running sum (`base` + the counts up to it; the slot may hold the count itself). The PV word and
    /// the helper counter are kept in registers and stored when their last cell is done, so the loop has no
    /// data-dependent branch (the stores are taken once per `2^bts` / 2048 cells). Returns the running sum after `c1`.
    #[inline(always)]
    pub(crate) fn range_bf(
        &mut self,
        c0: usize,
        c1: usize,
        bts: usize,
        at: impl FnMut(usize) -> (u32, *mut u32),
    ) -> u32 {
        self.range_bf_from(0, c0, c1, bts, at)
    }

    #[inline(always)]
    pub(crate) fn range_bf_from(
        &mut self,
        base: u32,
        c0: usize,
        c1: usize,
        bts: usize,
        mut at: impl FnMut(usize) -> (u32, *mut u32),
    ) -> u32 {
        let pv_last = (1usize << bts) - 1;
        let (mut word, mut hacc, mut end) = (0u32, 0u32, base);
        for c in c0..c1 {
            let (cnt, slot) = at(c);
            let nz = (cnt != 0) as u32;
            word |= nz << (c & 31);
            hacc = hacc.wrapping_add(cnt.wrapping_sub(nz));
            end = end.wrapping_add(cnt);
            // SAFETY: `at` returns a valid slot of the counts array for cell `c`.
            unsafe { *slot = end };
            if c & pv_last == pv_last || c + 1 == c1 {
                if let Some(w) = self.pv.get_mut(c >> bts) {
                    *w |= word;
                }
                word = 0;
            }
            if c & 2047 == 2047 || c + 1 == c1 {
                self.helper[c >> 11] = self.helper[c >> 11].wrapping_add(hacc);
                hacc = 0;
            }
        }
        end
    }

    /// Panics unless the arrays are the fill's.
    pub(crate) fn check(&self, ph: &PvHelper<'_>) {
        assert!(
            self.pv == ph.pv && self.helper == ph.helper,
            "R11 csrchk: the presence vector / helper counters from the counts differ from the fill's"
        );
        eprintln!(
            "[R11] csrchk pvh equal pv_words={} helper_counters={}",
            self.pv.len(),
            self.helper.len()
        );
    }
}

/// `csrchk`: builds the CSR of the chains on a clone of `next_pos` and compares it with the one
/// in `next_pos`'s capacity (panics on a difference).
pub(crate) fn check_against_chains(hashtable: &[u32], next_pos: &Vec<u32>, words: usize) {
    let cells = hashtable.len();
    let (o1, q1) = {
        let (o, q) = view(hashtable, next_pos);
        (o.to_vec(), q.to_vec())
    };
    let mut np = next_pos.clone();
    attach(hashtable, &mut np);
    let (o2, q2) = view(hashtable, &np);
    assert!(
        o1 == o2 && q1 == q2,
        "R11 csrchk: the words build differs from the chain walk"
    );
    eprintln!("[R11] csrchk equal cells={cells} words={words}");
}

/// Builds the CSR from the words into the spare capacity of `next_pos` (see the module
/// comment); prints one line `[R11] csr ...` with the time of each pass. With `ph`, the sum pass
/// also makes the presence vector and the helper counters from the counts (scratch). With
/// `csrchk`, the chain walk of `attach` builds it again and the two are compared.
/// R16: `cpf` (the passes through a ring of 16 words that prefetches each word's `offsets`
/// element), `sumlean` / `sum0` (`range_lean`), `offkeep` (the counts in the thread's kept buffer).
pub(crate) fn attach_words<W: Words>(
    w: &W,
    hashtable: &[u32],
    next_pos: &mut Vec<u32>,
    ph: Option<PvHelper<'_>>,
) {
    let m = mode();
    let t0 = Instant::now();
    let cells = hashtable.len();
    let n = next_pos.len();
    let need = cells + 1 + n;
    let asc = w.ascending();
    next_pos.reserve_exact(need);
    // SAFETY: `sp .. sp + need` is the reserved capacity of `next_pos`, apart from the vector's
    // own elements; `offsets` is initialized here before it is read, `qs` is only written (each
    // place once) before `view` reads it.
    let sp = unsafe { next_pos.as_mut_ptr().add(n) };
    // The counts: `offsets[cell]`, or `offsets[cell + 2]` when ascending (which also needs
    // `offsets[cells + 1]`, the first element of `qs`; `n >= words + 1 >= 1` leaves room).
    let k = if asc { 2 } else { 0 };
    let z = cells + 1 + asc as usize;
    // R16 `offkeep`: the counts in the thread's kept buffer (a fresh one when it is shorter).
    let mut kept: Vec<u32> = Vec::new();
    let mut kept_fresh = 0u8;
    let head: &mut [u32] = if m.offkeep {
        kept = XR16_KEPT.with(|b| std::mem::take(&mut *b.borrow_mut()));
        if kept.len() < z {
            kept = vec![0u32; z];
            kept_fresh = 1;
        } else {
            kept[..z].fill(0);
        }
        &mut kept[..z]
    } else {
        // SAFETY: `z <= need`; every element of `sp[..z]` is initialized here.
        unsafe {
            std::ptr::write_bytes(sp, 0, z);
            std::slice::from_raw_parts_mut(sp, z)
        }
    };
    let t1 = Instant::now();
    let mut words = 0usize;
    if m.cpf || m.cpfc {
        // R16 `cpf` / `cpfc`: each word's element is prefetched, and counted 16 words later.
        let hp = head.as_ptr();
        let mut ring = [0usize; 16];
        let mut seen = 0usize;
        w.each(&mut |cell, _q| {
            prefetch(hp.wrapping_add(cell + k));
            let slot = seen & 15;
            if seen >= 16 {
                head[ring[slot] + k] += 1;
            }
            ring[slot] = cell;
            seen += 1;
        });
        for i in seen.saturating_sub(16)..seen {
            head[ring[i & 15] + k] += 1;
        }
        words = seen;
    } else {
        w.each(&mut |cell, _q| {
            head[cell + k] += 1;
            words += 1;
        });
    }
    let t2 = Instant::now();
    let mut scratch = ph.as_ref().map(PvScratch::new);
    let mut end = 0u32;
    if m.sumlean || m.sum0 {
        // R16 `sumlean` / `sum0`: 32 cells per PV group, counters by difference (or no PV/helper).
        let pvh = !m.sum0 && scratch.is_some();
        let (pvp, pvl, hlp) = match scratch.as_mut() {
            Some(s) => (s.pv.as_mut_ptr(), s.pv.len(), s.helper.as_mut_ptr()),
            None => (std::ptr::null_mut(), 0, std::ptr::null_mut()),
        };
        let bts = ph.as_ref().map_or(5, |ph| ph.bts);
        // SAFETY: `head[k .. k + cells]` are the counts; the scratch arrays are this build's and
        // of the fill's lengths (`c >> 11 < helper.len()` for every cell).
        end = unsafe {
            range_lean(
                pvp,
                pvl,
                hlp,
                head.as_mut_ptr().add(k),
                0,
                0,
                cells,
                bts,
                pvh,
            )
        };
    } else if let (true, Some(s), Some(ph)) = (m.pvbf, scratch.as_mut(), ph.as_ref()) {
        // R15 `pvbf`: the same prefix sum, PV bits and helper counters without a branch per cell.
        end = s.range_bf(0, cells, ph.bts, |c| {
            let cnt = head[c + k];
            (cnt, &mut head[c + k] as *mut u32)
        });
    } else {
        for c in 0..cells {
            let cnt = head[c + k];
            if cnt != 0 {
                if let (Some(s), Some(ph)) = (scratch.as_mut(), ph.as_ref()) {
                    s.cell(c, cnt, ph.bts);
                }
            }
            end += cnt;
            head[c + k] = end;
        }
    }
    if !asc {
        head[cells] = end;
    }
    let t3 = Instant::now();
    // `ob`: the counts / offsets (`head`), `qb`: `qs` (after `offsets` in `next_pos`'s capacity).
    let ob = head.as_mut_ptr();
    // SAFETY: `cells + 1 <= need`.
    let qb = unsafe { sp.add(cells + 1) };
    // SAFETY (both): `cell + 1 <= cells`, `p < words <= n - 1`, so `cells + 1 + p < need`; the
    // elements read and written are initialized (`offsets`, or `qs[0]` only as the written target
    // of the last write to it, never read: `cells + 1` is no cell + 1).
    #[inline(always)]
    unsafe fn place_asc(ob: *mut u32, qb: *mut u32, cell: usize, q: u32) {
        unsafe {
            let p = *ob.add(cell + 1);
            *ob.add(cell + 1) = p + 1;
            *qb.add(p as usize) = q;
        }
    }
    // SAFETY: `cell < cells`, `offsets[cell] >= 1` (it is the end of a non-empty run),
    // `e < words <= n - 1`, so `cells + 1 + e < need`.
    #[inline(always)]
    unsafe fn place_desc(ob: *mut u32, qb: *mut u32, cell: usize, q: u32) {
        unsafe {
            let e = *ob.add(cell) - 1;
            *ob.add(cell) = e;
            *qb.add(e as usize) = q;
        }
    }
    let one = if asc { place_asc } else { place_desc };
    if m.cpf {
        // R16 `cpf`: each word's element is prefetched, and the word placed 16 words later (in
        // the same order, so each run gets its words in the same order).
        let mut ring = [(0usize, 0u32); 16];
        let mut seen = 0usize;
        w.each(&mut |cell, q| {
            assert!(cell < cells);
            prefetch(ob.wrapping_add(cell + asc as usize));
            let slot = seen & 15;
            if seen >= 16 {
                let (c, x) = ring[slot];
                // SAFETY: as `place_asc` / `place_desc`.
                unsafe { one(ob, qb, c, x) };
            }
            ring[slot] = (cell, q);
            seen += 1;
        });
        for i in seen.saturating_sub(16)..seen {
            let (c, x) = ring[i & 15];
            // SAFETY: as above.
            unsafe { one(ob, qb, c, x) };
        }
    } else if asc {
        w.each(&mut |cell, q| {
            assert!(cell < cells);
            // SAFETY: as `place_asc`.
            unsafe { place_asc(ob, qb, cell, q) };
        });
    } else {
        w.each(&mut |cell, q| {
            assert!(cell < cells);
            // SAFETY: as `place_desc`.
            unsafe { place_desc(ob, qb, cell, q) };
        });
        // `offsets[cells]` was the total and is no cell's end, so it stays the total.
    }
    let t4 = Instant::now();
    // R16 `offkeep`: the offsets into `next_pos`'s capacity (not in `build_ns`), the buffer back.
    let mut copy_ns = 0u128;
    if m.offkeep {
        let tc = Instant::now();
        // SAFETY: `cells + 1 <= need` elements of the capacity; the kept buffer has `z > cells`.
        unsafe { std::ptr::copy_nonoverlapping(kept.as_ptr(), sp, cells + 1) };
        copy_ns = tc.elapsed().as_nanos();
        XR16_KEPT.with(|b| *b.borrow_mut() = std::mem::take(&mut kept));
    }
    // EXPERIMENT (exp/r13-wall): the CSR build's interval (`R13_PROBE=tl`).
    crate::utils::xtiming::x_r13("csr", Some(t0), words as u64);
    eprintln!(
        "[R11] csr words-build cells={cells} words={words} next_pos_len={n} ascending={} pvh={} zero_ns={} count_ns={} sum_ns={} place_ns={} build_ns={} csr_bytes={} r16=cpf{}cpfc{}lean{}sum0{}keep{}fresh{} copy_ns={copy_ns}",
        asc as u8,
        ph.is_some() as u8,
        (t1 - t0).as_nanos(),
        (t2 - t1).as_nanos(),
        (t3 - t2).as_nanos(),
        (t4 - t3).as_nanos(),
        (t4 - t0).as_nanos(),
        (cells + 1 + words) * 4,
        m.cpf as u8,
        m.cpfc as u8,
        m.sumlean as u8,
        m.sum0 as u8,
        m.offkeep as u8,
        kept_fresh,
    );
    if m.check {
        check_against_chains(hashtable, next_pos, words);
        if let (Some(s), Some(ph), false) = (scratch.as_ref(), ph.as_ref(), m.sum0) {
            s.check(ph);
        }
    }
}
