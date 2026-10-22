//! EXPERIMENT (LOSAT_X_BNMAPPED): large per-thread buffers that glibc maps and returns when they
//! are dropped, and a pool that hands them from one thread to the next.
//!
//! No NCBI counterpart: these choose only the capacity of LOSAT buffers and when they are freed;
//! no value NCBI computes changes.
//!
//! glibc maps an allocation of at least its mmap threshold, but it raises the threshold to the
//! size of each mapped block that is freed (up to 32 MiB on 64-bit), so after a few frees the
//! large buffers of the worker threads come from their arenas, and the space they leave stays
//! resident there (`tblastn/x_subject_translate.rs` has the same case for TBLASTN's frames). A
//! buffer of at least `X_MAPPED_FROM` bytes is given a capacity above that ceiling, so that it
//! is mapped and returned when dropped (unless an arena has a free chunk that large); the pages
//! beyond its length are not resident until they are touched (a buffer lent by a pool may have
//! touched more than it holds: the pool counts that). Other allocators return large blocks on
//! their own, and 32-bit glibc's ceiling is 512 KiB, so elsewhere the capacity is what `Vec`
//! would give.
//!
//! A mapped block that is freed and mapped again is touched again: every page faults, and with
//! many threads the faults wait for each other (BLASTN of 1,500 cDNAs x the D. melanogaster
//! genome at 16 threads: 590,000 more minor faults, x1.04 wall). So a buffer that is done goes
//! to an `XPool` of its element type instead, shared by all threads, and the next buffer that
//! must grow past `X_MAPPED_FROM` takes one from there, its pages already resident. A thread's
//! arena keeps its free space for that thread only; the pool lends it to any thread, so the
//! resident buffers are those in use at the busiest moment, plus at most the pool's idle bound
//! of idle ones. Only buffers that are done after every subject use a pool: BLASTN's database
//! reads and BLASTNA letters (`X_BYTES`, idle bound `x_pool_idle_bytes`) and the scratch of the
//! ungapped hits' sort (`blastn` `X_SORT`, its own bound); a buffer that its thread keeps from
//! one subject to the next (the ungapped hits) is mapped but not pooled, since idle copies of it
//! would only add to those in use.

/// No NCBI counterpart: the buffers from which a buffer is mapped (Linux with glibc).
pub(crate) const X_MAPPED_FROM: usize = 1 << 20;

/// No NCBI counterpart: glibc's largest dynamic mmap threshold (DEFAULT_MMAP_THRESHOLD_MAX on
/// 64-bit), and one page more.
const X_MAPPED_BYTES: usize = (32 << 20) + 4096;

/// No NCBI counterpart: whether the allocator is 64-bit glibc's (the only one that keeps large
/// freed blocks this way; on 32 bits its ceiling is 512 KiB and a padded capacity would only
/// spend address space).
pub(crate) const X_MAPPED_ALLOCATOR: bool = cfg!(all(
    target_os = "linux",
    target_env = "gnu",
    target_pointer_width = "64"
));

/// No NCBI counterpart: whether a buffer of `bytes` is mapped (`x_mapped_capacity`).
pub(crate) fn x_is_mapped(bytes: usize) -> bool {
    X_MAPPED_ALLOCATOR && bytes >= X_MAPPED_FROM
}

/// No NCBI counterpart: the capacity, in elements of `size` bytes, of a buffer that must hold
/// `elements`: above glibc's ceiling once it reaches `X_MAPPED_FROM` bytes.
pub(crate) fn x_mapped_capacity(elements: usize, size: usize) -> usize {
    let size = size.max(1);
    if X_MAPPED_ALLOCATOR && elements.saturating_mul(size) >= X_MAPPED_FROM {
        elements.max(X_MAPPED_BYTES.div_ceil(size))
    } else {
        elements
    }
}

/// No NCBI counterpart: `Vec::reserve(additional)` (the same doubling), with the capacity of
/// `x_mapped_capacity` once the buffer reaches `X_MAPPED_FROM` bytes.
pub(crate) fn x_reserve_mapped<T>(v: &mut Vec<T>, additional: usize) {
    let need = v.len().saturating_add(additional);
    if v.capacity() >= need {
        return;
    }
    let grown = need.max(v.capacity().saturating_mul(2));
    let target = x_mapped_capacity(grown, std::mem::size_of::<T>());
    v.reserve_exact(target - v.len());
}

/// No NCBI counterpart: the built-in bound of `x_pool_idle_bytes`, in MiB (BLASTN of 1,500 cDNAs
/// x the D. melanogaster genome at 16 threads, max RSS: 273 MiB with 0, 317-331 with 64, 346-367
/// with 256; with 0 every reuse faults again, x1.04 wall).
const X_POOL_IDLE_MIB: usize = 64;

/// No NCBI counterpart: the idle buffers of `X_BYTES` may hold this many touched bytes
/// (`LOSAT_X_POOLIDLE` MiB, a tuning variable; built-in `X_POOL_IDLE_MIB`); the oldest ones
/// beyond it are freed.
pub(crate) fn x_pool_idle_bytes() -> usize {
    static IDLE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *IDLE.get_or_init(|| x_idle_bytes("LOSAT_X_POOLIDLE", X_POOL_IDLE_MIB))
}

/// No NCBI counterpart: the idle bound of a pool in bytes from the tuning variable `var`, in MiB
/// (0 to 65,536; unset or any other value gives `builtin_mib`).
pub(crate) fn x_idle_bytes(var: &str, builtin_mib: usize) -> usize {
    std::env::var(var)
        .ok()
        .and_then(|v| v.trim().parse::<usize>().ok())
        .filter(|v| *v <= 65_536)
        .unwrap_or(builtin_mib)
        .saturating_mul(1 << 20)
}

/// No NCBI counterpart: idle buffers of at least `X_MAPPED_FROM` bytes, shared by all threads
/// (see the module header). Each idle buffer is kept with the bytes it may have touched (the
/// most it held while lent out, as far as the pool saw it).
pub(crate) struct XPool<T> {
    state: std::sync::Mutex<XPoolState<T>>,
    // the bound of the idle buffers' touched bytes
    idle_bytes: fn() -> usize,
}

struct XPoolState<T> {
    idle: Vec<(Vec<T>, usize)>,
    idle_touched: usize,
    // address -> touched bytes of the buffers lent out
    lent: Vec<(usize, usize)>,
}

impl<T: Copy> XPool<T> {
    /// A pool whose idle buffers may hold `idle_bytes()` touched bytes.
    pub(crate) const fn new(idle_bytes: fn() -> usize) -> Self {
        Self {
            state: std::sync::Mutex::new(XPoolState {
                idle: Vec::new(),
                idle_touched: 0,
                lent: Vec::new(),
            }),
            idle_bytes,
        }
    }

    /// An empty idle buffer with room for `need` elements, if there is one: of those that have
    /// touched enough bytes already, the one that touched the fewest, else the one that touched
    /// the most (all are mapped, so capacity says nothing about the resident pages). Only a
    /// buffer that holds `need` is lent, so that it is not moved before it is given back (its
    /// touched bytes are kept by its address).
    pub(crate) fn take(&self, need: usize) -> Option<Vec<T>> {
        if !X_MAPPED_ALLOCATOR {
            return None;
        }
        let bytes = need.saturating_mul(std::mem::size_of::<T>().max(1));
        let mut st = self.state.lock().unwrap_or_else(|e| e.into_inner());
        let fits = |(_, (v, _)): &(usize, &(Vec<T>, usize))| v.capacity() >= need;
        let pick = st
            .idle
            .iter()
            .enumerate()
            .filter(fits)
            .filter(|(_, (_, touched))| *touched >= bytes)
            .min_by_key(|(_, (_, touched))| *touched)
            .or_else(|| {
                st.idle
                    .iter()
                    .enumerate()
                    .filter(fits)
                    .max_by_key(|(_, (_, touched))| *touched)
            })
            .map(|(i, _)| i)?;
        // `remove`, not `swap_remove`: `give` frees the oldest idle buffers first.
        let (v, touched) = st.idle.remove(pick);
        st.idle_touched -= touched;
        let addr = v.as_ptr() as usize;
        // A buffer lent before at this address and never given back (dropped on an error path)
        // left its entry; this one replaces it.
        st.lent.retain(|&(a, _)| a != addr);
        st.lent.push((addr, touched));
        Some(v)
    }

    /// `x_reserve_mapped(v, additional)`; when `v` then needs `X_MAPPED_FROM` bytes or more for
    /// the first time, its contents move into an idle buffer of the pool first, if there is one.
    pub(crate) fn reserve(&self, v: &mut Vec<T>, additional: usize) {
        let size = std::mem::size_of::<T>().max(1);
        let need = v.len().saturating_add(additional);
        if v.capacity() >= need {
            return;
        }
        if v.capacity().saturating_mul(size) < X_MAPPED_FROM
            && need.saturating_mul(size) >= X_MAPPED_FROM
        {
            if let Some(mut pooled) = self.take(need) {
                pooled.extend_from_slice(v);
                *v = pooled;
            }
        }
        x_reserve_mapped(v, additional);
    }

    /// Puts `v` into the pool, empty, when its capacity is at least `X_MAPPED_FROM` bytes (`v`
    /// becomes a new empty `Vec`); a smaller `v` stays as it is. Idle buffers beyond the pool's
    /// bound of touched bytes are freed, the oldest first; a buffer that alone touched
    /// more is freed at once.
    pub(crate) fn give(&self, v: &mut Vec<T>) {
        let size = std::mem::size_of::<T>().max(1);
        if !X_MAPPED_ALLOCATOR || v.capacity().saturating_mul(size) < X_MAPPED_FROM {
            return;
        }
        let mut v = std::mem::take(v);
        let used = v.len() * size;
        v.clear();
        let mut freed = Vec::new();
        {
            let mut st = self.state.lock().unwrap_or_else(|e| e.into_inner());
            let addr = v.as_ptr() as usize;
            let before = match st.lent.iter().position(|&(a, _)| a == addr) {
                Some(i) => st.lent.swap_remove(i).1,
                None => 0,
            };
            let touched = used.max(before);
            let cap = (self.idle_bytes)();
            if touched > cap {
                freed.push(v);
            } else {
                st.idle_touched += touched;
                st.idle.push((v, touched));
            }
            while st.idle_touched > cap && !st.idle.is_empty() {
                let (old, t) = st.idle.remove(0);
                st.idle_touched -= t;
                freed.push(old);
            }
        }
        drop(freed);
    }
}

/// No NCBI counterpart: the pool of byte buffers: BLASTN's database reads and BLASTNA letters
/// (LOSAT_X_BNMAPPED).
pub(crate) static X_BYTES: XPool<u8> = XPool::new(x_pool_idle_bytes);

/// No NCBI counterpart: the elements sorted at a time with `sort_by` before they are merged in
/// `x_sort_by_mapped` (its scratch stays below glibc's first mmap threshold).
const X_SORT_RUN_BYTES: usize = 64 << 10;

/// No NCBI counterpart: `v.sort_by(cmp)` (a stable sort: the same order as `sort_by` for a
/// comparator that is a total order) whose scratch, when it would need `X_MAPPED_FROM` bytes or
/// more with glibc, is a mapped buffer as long as `v` from `pool` (given back after the sort):
/// runs of `X_SORT_RUN_BYTES` sorted with `sort_by`, then merged pairwise from `v` into the
/// buffer and back (the left run first among equal elements). A scratch of half of `v`, the left
/// run copied out before each merge, held no less at the peak (Sakai x EDL933 `-subject` blastn,
/// whose largest lists hold 11.4 MiB of hits, and the D. melanogaster genome `-db`).
pub(crate) fn x_sort_by_mapped<T: Copy, F: FnMut(&T, &T) -> std::cmp::Ordering>(
    v: &mut [T],
    mut cmp: F,
    pool: &XPool<T>,
) {
    let size = std::mem::size_of::<T>().max(1);
    let n = v.len();
    if !X_MAPPED_ALLOCATOR || n.saturating_mul(size) < X_MAPPED_FROM {
        v.sort_by(cmp);
        return;
    }
    let run = (X_SORT_RUN_BYTES / size).max(1);
    for chunk in v.chunks_mut(run) {
        chunk.sort_by(&mut cmp);
    }
    x_merge_runs(v, &mut cmp, run, pool);
}

/// No NCBI counterpart: the merges of `x_sort_by_mapped`, `v` being sorted in runs of `run`
/// elements.
fn x_merge_runs<T: Copy, F: FnMut(&T, &T) -> std::cmp::Ordering>(
    v: &mut [T],
    cmp: &mut F,
    run: usize,
    pool: &XPool<T>,
) {
    let n = v.len();
    let mut buf: Vec<T> = Vec::new();
    pool.reserve(&mut buf, n);
    buf.extend_from_slice(v);
    let mut width = run;
    let mut in_v = true;
    while width < n {
        let (src, dst): (&[T], &mut [T]) = if in_v {
            (&*v, &mut buf[..])
        } else {
            (&buf[..], &mut *v)
        };
        let mut start = 0;
        while start < n {
            let mid = (start + width).min(n);
            let end = (start + 2 * width).min(n);
            let (mut i, mut j, mut k) = (start, mid, start);
            while i < mid && j < end {
                if cmp(&src[j], &src[i]) == std::cmp::Ordering::Less {
                    dst[k] = src[j];
                    j += 1;
                } else {
                    dst[k] = src[i];
                    i += 1;
                }
                k += 1;
            }
            dst[k..k + (mid - i)].copy_from_slice(&src[i..mid]);
            k += mid - i;
            dst[k..k + (end - j)].copy_from_slice(&src[j..end]);
            start = end;
        }
        width *= 2;
        in_v = !in_v;
    }
    if !in_v {
        v.copy_from_slice(&buf);
    }
    pool.give(&mut buf);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capacity_is_mapped_from_the_bound_on_glibc_only() {
        assert_eq!(x_mapped_capacity(X_MAPPED_FROM - 1, 1), X_MAPPED_FROM - 1);
        let at = x_mapped_capacity(X_MAPPED_FROM, 1);
        let rows = x_mapped_capacity(X_MAPPED_FROM / 12 + 1, 12);
        if X_MAPPED_ALLOCATOR {
            assert_eq!(at, X_MAPPED_BYTES);
            assert!(rows * 12 >= X_MAPPED_BYTES && (rows - 1) * 12 < X_MAPPED_BYTES);
        } else {
            assert_eq!(at, X_MAPPED_FROM);
            assert_eq!(rows, X_MAPPED_FROM / 12 + 1);
        }
        assert_eq!(x_mapped_capacity(40 << 20, 1), 40 << 20);
    }

    #[test]
    fn reserve_doubles_as_vec_below_the_bound() {
        let mut v: Vec<u8> = Vec::with_capacity(1000);
        v.resize(1000, 7);
        x_reserve_mapped(&mut v, 1);
        assert_eq!(v.capacity(), 2000);
        assert_eq!(v.len(), 1000);
        x_reserve_mapped(&mut v, 10);
        assert_eq!(v.capacity(), 2000);
        x_reserve_mapped(&mut v, X_MAPPED_FROM);
        assert!(v.capacity() >= 1000 + X_MAPPED_FROM);
        assert!(v.iter().all(|&b| b == 7));
    }

    #[test]
    fn sort_by_mapped_is_stable_and_equals_sort_by() {
        // (key, position): sorted by key only, so equal keys must keep their order.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        for &n in &[0usize, 1, 5, 7, 9, 70_000, 200_003] {
            let v: Vec<(u32, u32)> = (0..n as u32)
                .map(|i| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    ((state % 1000) as u32, i)
                })
                .collect();
            let mut a = v.clone();
            a.sort_by(|x, y| y.0.cmp(&x.0));
            let mut b = v.clone();
            static POOL: XPool<(u32, u32)> = XPool::new(|| 64 << 20);
            x_sort_by_mapped(&mut b, |x, y| y.0.cmp(&x.0), &POOL);
            assert_eq!(a, b, "n = {n}");
            // The merges alone, with runs of 4 (every merge path, odd tails).
            let mut c = v.clone();
            for run in c.chunks_mut(4) {
                run.sort_by(|x, y| y.0.cmp(&x.0));
            }
            x_merge_runs(
                &mut c,
                &mut |x: &(u32, u32), y: &(u32, u32)| y.0.cmp(&x.0),
                4,
                &POOL,
            );
            assert_eq!(a, c, "n = {n}, runs of 4");
        }
    }

    #[test]
    fn pool_lends_idle_buffers_and_keeps_small_ones_out() {
        let pool: XPool<u8> = XPool::new(|| 64 << 20);
        let mut small = vec![1u8; 100];
        pool.give(&mut small);
        assert_eq!(small.len(), 100);
        assert!(pool.take(1).is_none());
        let mut v: Vec<u8> = vec![5; 10];
        pool.reserve(&mut v, X_MAPPED_FROM);
        assert_eq!(&v[..], &[5u8; 10][..]);
        assert!(v.capacity() >= 10 + X_MAPPED_FROM);
        let addr = v.as_ptr();
        v.resize(X_MAPPED_FROM + 10, 0);
        pool.give(&mut v);
        // Elsewhere than glibc `give` keeps nothing and leaves `v` as it is.
        assert_eq!(v.capacity() == 0, X_MAPPED_ALLOCATOR);
        if X_MAPPED_ALLOCATOR {
            let mut w: Vec<u8> = vec![9; 3];
            pool.reserve(&mut w, X_MAPPED_FROM);
            assert_eq!(w.as_ptr(), addr);
            assert_eq!(&w[..], &[9u8; 3][..]);
        }
        // Only a buffer that holds what is asked is lent (a smaller one would be moved).
        let mut big: Vec<u8> = Vec::new();
        x_reserve_mapped(&mut big, X_MAPPED_FROM);
        let cap = big.capacity();
        pool.give(&mut big);
        assert!(pool.take(cap + 1).is_none());
        if X_MAPPED_ALLOCATOR {
            assert_eq!(pool.take(cap).map(|v| v.capacity()), Some(cap));
        }
    }

    #[test]
    fn pool_frees_beyond_its_bound() {
        // A buffer that alone touched more than the bound is not kept.
        let none: XPool<u8> = XPool::new(|| 0);
        let mut v: Vec<u8> = Vec::new();
        x_reserve_mapped(&mut v, X_MAPPED_FROM);
        v.resize(X_MAPPED_FROM, 1);
        none.give(&mut v);
        assert_eq!(v.capacity() == 0, X_MAPPED_ALLOCATOR);
        assert!(none.take(1).is_none());
        // With room for one touched buffer, the older one is freed when a second comes back.
        let one: XPool<u8> = XPool::new(|| X_MAPPED_FROM + X_MAPPED_FROM / 2);
        let (mut a, mut b): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
        x_reserve_mapped(&mut a, X_MAPPED_FROM);
        x_reserve_mapped(&mut b, X_MAPPED_FROM);
        a.resize(X_MAPPED_FROM, 1);
        b.resize(X_MAPPED_FROM, 2);
        let b_addr = b.as_ptr();
        one.give(&mut a);
        one.give(&mut b);
        if X_MAPPED_ALLOCATOR {
            let kept = one.take(1).expect("the newer buffer is kept");
            assert_eq!(kept.as_ptr(), b_addr);
            assert!(one.take(1).is_none());
        } else {
            assert!(one.take(1).is_none());
        }
    }
}
