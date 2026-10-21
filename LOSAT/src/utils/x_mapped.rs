//! EXPERIMENT (LOSAT_X_XDROPTRIM, LOSAT_X_BNMAPPED): large per-thread buffers that glibc maps
//! and returns when they are dropped.
//!
//! No NCBI counterpart: these choose only the capacity of LOSAT buffers and when they are freed;
//! no value NCBI computes changes.
//!
//! glibc maps an allocation of at least its mmap threshold, but it raises the threshold to the
//! size of each mapped block that is freed (up to 32 MiB on 64-bit), so after a few frees the
//! large buffers of the worker threads come from their arenas, and the space they leave stays
//! resident there (`tblastn/x_subject_translate.rs` has the same case for TBLASTN's frames). A
//! buffer of at least `X_MAPPED_FROM` bytes is given a capacity above that ceiling, so that it
//! is always mapped and returned when dropped; the pages beyond its length are never touched and
//! are not resident. Other allocators return large blocks on their own, so elsewhere the
//! capacity is what `Vec` would give.

/// No NCBI counterpart: the buffers from which a buffer is mapped (Linux with glibc).
pub(crate) const X_MAPPED_FROM: usize = 1 << 20;

/// No NCBI counterpart: glibc's largest dynamic mmap threshold (DEFAULT_MMAP_THRESHOLD_MAX on
/// 64-bit), and one page more.
const X_MAPPED_BYTES: usize = (32 << 20) + 4096;

/// No NCBI counterpart: whether the allocator is glibc's (the only one that keeps large freed
/// blocks this way).
pub(crate) const X_MAPPED_ALLOCATOR: bool = cfg!(all(target_os = "linux", target_env = "gnu"));

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

/// No NCBI counterpart: frees `v` when its capacity is at least `X_MAPPED_FROM` bytes (with
/// glibc, then a mapped block that goes back to the system).
pub(crate) fn x_release_mapped<T>(v: &mut Vec<T>) {
    if v.capacity().saturating_mul(std::mem::size_of::<T>()) >= X_MAPPED_FROM {
        *v = Vec::new();
    }
}

/// No NCBI counterpart: the elements sorted at a time with `sort_by` before they are merged in
/// `x_sort_by_mapped` (its scratch stays below glibc's first mmap threshold).
const X_SORT_RUN_BYTES: usize = 64 << 10;

/// No NCBI counterpart: `v.sort_by(cmp)` (a stable sort: the same order as `sort_by` for a
/// comparator that is a total order) whose scratch, when it would need `X_MAPPED_FROM` bytes or
/// more with glibc, is `buf` (kept by the caller, grown with `x_reserve_mapped`): runs of
/// `X_SORT_RUN_BYTES` sorted with `sort_by`, then merged pairwise (the left run first among equal
/// elements).
pub(crate) fn x_sort_by_mapped<T: Copy, F: FnMut(&T, &T) -> std::cmp::Ordering>(
    v: &mut [T],
    mut cmp: F,
    buf: &mut Vec<T>,
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
    buf.clear();
    x_reserve_mapped(buf, n);
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
    fn reserve_doubles_as_vec_below_the_bound_and_release_keeps_small_buffers() {
        let mut v: Vec<u8> = Vec::with_capacity(1000);
        v.resize(1000, 7);
        x_reserve_mapped(&mut v, 1);
        assert_eq!(v.capacity(), 2000);
        assert_eq!(v.len(), 1000);
        x_reserve_mapped(&mut v, 10);
        assert_eq!(v.capacity(), 2000);
        x_release_mapped(&mut v);
        assert_eq!(v.capacity(), 2000);
        x_reserve_mapped(&mut v, X_MAPPED_FROM);
        assert!(v.capacity() >= 1000 + X_MAPPED_FROM);
        assert!(v.iter().all(|&b| b == 7));
        x_release_mapped(&mut v);
        assert_eq!(v.capacity(), 0);
    }

    #[test]
    fn sort_by_mapped_is_stable_and_equals_sort_by() {
        // (key, position): sorted by key only, so equal keys must keep their order.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        for &n in &[0usize, 1, 7, 70_000, 200_003] {
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
            let mut buf = vec![(7, 7); 3];
            x_sort_by_mapped(&mut b, |x, y| y.0.cmp(&x.0), &mut buf);
            assert_eq!(a, b, "n = {n}");
        }
    }
}
