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
                    _ => {}
                }
            }
        }
        m
    })
}

#[inline(always)]
fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch does not access memory architecturally; any address is allowed.
    unsafe {
        core::arch::x86_64::_mm_prefetch::<{ core::arch::x86_64::_MM_HINT_T0 }>(p as *const i8);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = p;
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
pub(crate) fn attach_words<W: Words>(
    w: &W,
    hashtable: &[u32],
    next_pos: &mut Vec<u32>,
    ph: Option<PvHelper<'_>>,
) {
    let t0 = Instant::now();
    let cells = hashtable.len();
    let n = next_pos.len();
    let need = cells + 1 + n;
    let asc = w.ascending();
    next_pos.reserve_exact(need);
    let ptr = next_pos.as_mut_ptr();
    // SAFETY: `ptr + n .. ptr + n + need` is the reserved capacity of `next_pos`, apart from the
    // vector's own elements; `offsets` is initialized here before it is read, `qs` is only
    // written (each place once) before `view` reads it.
    let spare = unsafe {
        std::slice::from_raw_parts_mut(ptr.add(n) as *mut std::mem::MaybeUninit<u32>, need)
    };
    // The counts: `offsets[cell]`, or `offsets[cell + 2]` when ascending (which also needs
    // `offsets[cells + 1]`, the first element of `qs`; `n >= words + 1 >= 1` leaves room).
    let k = if asc { 2 } else { 0 };
    let z = cells + 1 + asc as usize;
    spare[..z].fill(std::mem::MaybeUninit::new(0));
    // SAFETY: every element of `spare[..z]` was just initialized.
    let head: &mut [u32] =
        unsafe { &mut *(&mut spare[..z] as *mut [std::mem::MaybeUninit<u32>] as *mut [u32]) };
    let t1 = Instant::now();
    let mut words = 0usize;
    w.each(&mut |cell, _q| {
        head[cell + k] += 1;
        words += 1;
    });
    let t2 = Instant::now();
    let mut scratch = ph.as_ref().map(PvScratch::new);
    let mut end = 0u32;
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
    if !asc {
        head[cells] = end;
    }
    let t3 = Instant::now();
    let base = spare.as_mut_ptr() as *mut u32;
    if asc {
        w.each(&mut |cell, q| {
            assert!(cell < cells);
            // SAFETY: `cell + 1 <= cells`, `p < words <= n - 1`, so `cells + 1 + p < need`; the
            // elements read and written are initialized (`offsets`, or `qs[0]` only as the
            // written target of the last write to it, never read: `cells + 1` is no cell + 1).
            unsafe {
                let p = *base.add(cell + 1);
                *base.add(cell + 1) = p + 1;
                *base.add(cells + 1 + p as usize) = q;
            }
        });
    } else {
        w.each(&mut |cell, q| {
            assert!(cell < cells);
            // SAFETY: `cell < cells`, `offsets[cell] >= 1` (it is the end of a non-empty run),
            // `e < words <= n - 1`, so `cells + 1 + e < need`.
            unsafe {
                let e = *base.add(cell) - 1;
                *base.add(cell) = e;
                *base.add(cells + 1 + e as usize) = q;
            }
        });
        // `offsets[cells]` was the total and is no cell's end, so it stays the total.
    }
    let t4 = Instant::now();
    // EXPERIMENT (exp/r13-wall): the CSR build's interval (`R13_PROBE=tl`).
    crate::utils::xtiming::x_r13("csr", Some(t0), words as u64);
    eprintln!(
        "[R11] csr words-build cells={cells} words={words} next_pos_len={n} ascending={} pvh={} zero_ns={} count_ns={} sum_ns={} place_ns={} build_ns={} csr_bytes={}",
        asc as u8,
        ph.is_some() as u8,
        (t1 - t0).as_nanos(),
        (t2 - t1).as_nanos(),
        (t3 - t2).as_nanos(),
        (t4 - t3).as_nanos(),
        (t4 - t0).as_nanos(),
        (cells + 1 + words) * 4,
    );
    if mode().check {
        check_against_chains(hashtable, next_pos, words);
        if let (Some(s), Some(ph)) = (scratch.as_ref(), ph.as_ref()) {
            s.check(ph);
        }
    }
}
