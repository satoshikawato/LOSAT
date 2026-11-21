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
//! reads and BLASTNA letters (`X_BYTES`, idle bound `x_pool_idle_bytes`), BLASTP's database blocks
//! of `X_MAPPED_FROM` bytes or more (`X_BYTES` too), TBLASTN's `-db` reads
//! and, while a `-db` stream is open, its large frames (`X_BYTES` too, `tblastn/x_db_stream.rs`),
//! and the scratch of the ungapped hits' sort (`blastn` `X_SORT`, its own bound); a buffer that
//! its thread keeps from one subject to the next (the ungapped hits) is mapped but not pooled,
//! since idle copies of it would only add to those in use. A worker that keeps its own buffer
//! below `X_MAPPED_FROM` from one read to the next (the database blocks of BLASTN and BLASTP,
//! `x_filled_kept`) takes a buffer of the pool for a longer read only and gives it back after
//! (`x_keep`), so that it does not keep an outlier's room.

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

/// No NCBI counterpart (LOSAT_X_BPDBKEEP; LOSAT_X_TBNLUTLEAN, TBLASTN's lookup tables): the
/// capacity, in elements of `size` bytes, of a buffer
/// of `elements` that glibc maps whatever its size and whose free leaves the mmap threshold as it
/// is (above glibc's ceiling; `elements` elsewhere). Its pages beyond those touched are not
/// resident.
pub(crate) fn x_unthresholded_capacity(elements: usize, size: usize) -> usize {
    if X_MAPPED_ALLOCATOR {
        elements.max(X_MAPPED_BYTES.div_ceil(size.max(1)))
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
        self.take_touched(need, usize::MAX).map(|mut v| {
            v.clear();
            v
        })
    }

    /// `reserve` of an empty buffer for `len` elements, then the buffer filled to `len` elements:
    /// a buffer taken from the pool keeps the elements its last user left in it (`give` keeps
    /// them), so only the elements past those are set to `fill`, as `resize` sets only the
    /// elements past a buffer's length. For a reader that then writes all `len` elements (the
    /// database's block reads, `seqdb::file`), which need not have them zeroed first.
    pub(crate) fn x_filled(&self, len: usize, fill: T) -> Vec<T> {
        let mut v = Vec::new();
        if len.saturating_mul(std::mem::size_of::<T>().max(1)) >= X_MAPPED_FROM {
            if let Some(pooled) = self.take_touched(len, usize::MAX) {
                v = pooled;
                v.truncate(len);
            }
        }
        let additional = len - v.len();
        x_reserve_mapped(&mut v, additional);
        v.resize(len, fill);
        v
    }

    /// No NCBI counterpart: `x_filled(len, fill)` for a worker that keeps `kept`, its own buffer,
    /// from one read to the next (a worker's database blocks, `seqdb::file`; `x_keep` gives the
    /// buffer back). Below `X_MAPPED_FROM` bytes the buffer is `kept`, taken, with the elements its
    /// last read left in it, grown to exactly `len` if it is shorter; only the elements past its
    /// length are set to `fill`. From `X_MAPPED_FROM` bytes on, when `kept` does not hold `len`
    /// (it never does after `x_keep`), it is `x_filled`'s buffer of the pool, and `kept` stays as
    /// it is, so that the worker does not keep an outlier's room for its next reads.
    pub(crate) fn x_filled_kept(&self, kept: &mut Vec<T>, len: usize, fill: T) -> Vec<T> {
        let size = std::mem::size_of::<T>().max(1);
        if len.saturating_mul(size) >= X_MAPPED_FROM && kept.capacity() < len {
            return self.x_filled(len, fill);
        }
        let mut v = std::mem::take(kept);
        if v.capacity() < len {
            v.reserve_exact(len - v.len());
        }
        v.resize(len, fill);
        v
    }

    /// No NCBI counterpart: gives back a buffer of `x_filled_kept`: one of `X_MAPPED_FROM` bytes
    /// or more to the pool (`give`; freed where the pool keeps none), a smaller one to `kept`, so
    /// that a worker keeps only a buffer below `X_MAPPED_FROM` bytes.
    pub(crate) fn x_keep(&self, kept: &mut Vec<T>, mut v: Vec<T>) {
        if v.capacity().saturating_mul(std::mem::size_of::<T>().max(1)) >= X_MAPPED_FROM {
            self.give(&mut v);
        } else {
            *kept = v;
        }
    }

    /// `take` among the idle buffers that touched at most twice
    /// the bytes of `need`: a buffer lent for as long as it is used (a TBLASTN frame) holds its
    /// touched pages, which the pool's idle bound no longer counts.
    pub(crate) fn x_take_near(&self, need: usize) -> Option<Vec<T>> {
        let bytes = need.saturating_mul(std::mem::size_of::<T>().max(1));
        self.take_touched(need, bytes.saturating_mul(2))
            .map(|mut v| {
                v.clear();
                v
            })
    }

    /// The idle buffer `take` lends among those that touched at most `most` bytes, with the
    /// elements its last user left in it (`give`).
    fn take_touched(&self, need: usize, most: usize) -> Option<Vec<T>> {
        if !X_MAPPED_ALLOCATOR {
            return None;
        }
        let bytes = need.saturating_mul(std::mem::size_of::<T>().max(1));
        let mut st = self.state.lock().unwrap_or_else(|e| e.into_inner());
        let fits = |(_, (v, touched)): &(usize, &(Vec<T>, usize))| {
            v.capacity() >= need && *touched <= most
        };
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
        // A buffer lent before at this address and never given back (dropped on an error path,
        // or kept by its user without `x_disown`) left its entry; this one replaces it.
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

    /// Forgets the buffer at `v`, lent by `take` and kept by its
    /// user (a TBLASTN store made of a read buffer, freed and never given back), so that a buffer
    /// mapped later at the same address does not inherit its touched bytes.
    pub(crate) fn x_disown(&self, v: &[T]) {
        let addr = v.as_ptr() as usize;
        let mut st = self.state.lock().unwrap_or_else(|e| e.into_inner());
        st.lent.retain(|&(a, _)| a != addr);
    }

    /// Frees every idle buffer (the last TBLASTN `-db` stream of
    /// the process is closed, and its frames would stay resident until the next search).
    pub(crate) fn x_clear_idle(&self) {
        let freed = {
            let mut st = self.state.lock().unwrap_or_else(|e| e.into_inner());
            st.idle_touched = 0;
            std::mem::take(&mut st.idle)
        };
        drop(freed);
    }

    /// Puts `v` into the pool when its capacity is at least `X_MAPPED_FROM` bytes (`v` becomes a
    /// new empty `Vec`); a smaller `v` stays as it is. The pool keeps its elements (`Copy`
    /// values, so nothing is dropped), which `x_filled` does not write again; `take` empties it.
    /// Idle buffers beyond the pool's bound of touched bytes are freed, the oldest first; a buffer
    /// that alone touched more is freed at once.
    pub(crate) fn give(&self, v: &mut Vec<T>) {
        let size = std::mem::size_of::<T>().max(1);
        if !X_MAPPED_ALLOCATOR || v.capacity().saturating_mul(size) < X_MAPPED_FROM {
            return;
        }
        let v = std::mem::take(v);
        let used = v.len() * size;
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
/// (LOSAT_X_BNMAPPED; the preliminary stage's blocks from `X_MAPPED_FROM` bytes on), BLASTP's
/// database blocks from `X_MAPPED_FROM` bytes on (LOSAT_X_BPDBKEEP, a protein of a million
/// residues or more; the pool keeps its idle buffers to the end of the process, as BLASTN's),
/// TBLASTN's `-db` reads and large frames.
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
#[path = "x_mapped_tests.rs"]
mod tests;
