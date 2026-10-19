//! EXPERIMENT (exp/r11-csr, v4-R R11; evidence only, never merged): the megablast lookup's chains
//! laid out as CSR (compressed sparse rows): `offsets[cell] .. offsets[cell + 1]` indexes `qs`,
//! which holds the cell's 1-based query offsets in the order the chain lists them
//! (`hashtable[cell]`, `next_pos[...]`, ... up to 0). Chosen from `R4_PROBE` (comma list):
//! - `csr`: build the CSR after the fill and let the probed scans read it instead of the chains.
//! - `csrb`: build it only (its cost and resident memory), the scans walk the chains.
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
