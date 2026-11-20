//! `std::sort` of GCC's C++ library (libstdc++ 13, the library NCBI's rpsbproc 0.5.1 binary is
//! built with, GCC 13.2): introsort with a median-of-three pivot, a heap sort where the
//! recursion gets too deep, and a final insertion sort. It is not stable: elements that the
//! comparison does not order end in the order this algorithm leaves them, which NCBI's output
//! then shows (the domain hits of equal E-value and bit score, the lookup tables of the CDD
//! data). Every step below is the library's, so that the order is the same for every input.
//!
//! The functions take a slice and a "less" comparison (`comp(a, b)`: `a` goes before `b`), the
//! `_Iter_comp_iter` of `std::sort(first, last, comp)`. With a strict weak order libstdc++'s
//! unguarded loops stay inside the range. With a comparison that is not one (rpsbproc's of hits
//! with NaN E-values or bit scores), the algorithm is still libstdc++'s, step for step, and its
//! result is defined while its loops stay inside the vector. Two of them have no bound of their
//! own: the first scan of `__unguarded_partition` (stopped by a sentinel that the median of three
//! and the swaps leave) and the walk of `__unguarded_linear_insert` in `__final_insertion_sort`
//! after the first 16 elements (stopped by an element that a partition left before it). A
//! comparison that is not irreflexive and asymmetric, or the leftmost part when the depth limit
//! hands it to the heap sort (every other part has a partition's pivot left of it), can take them
//! past the vector's ends, where NCBI reads outside it (undefined): `try_std_sort` gives
//! `OutOfRange` before such a read. (`introsort_loop` compares positions as libstdc++'s signed
//! differences.)

/// `std::sort(first, last, comp)`.
// libstdc++ (GCC 13): bits/stl_algo.h:1940-1951
// libstdc++ (GCC 13): bits/stl_algobase.h:1513-1518
pub(crate) fn std_sort<T: Copy>(slice: &mut [T], comp: &mut impl FnMut(&T, &T) -> bool) {
    try_std_sort(slice, comp).expect("a strict weak order keeps libstdc++'s sort in its range");
}

/// A loop of libstdc++'s sort that would read outside the vector (a comparison that is not a
/// strict weak order).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct OutOfRange;

/// `std::sort(first, last, comp)` of any comparison: libstdc++'s permutation, or `OutOfRange`
/// where one of its unguarded loops would leave the slice (the module's text).
pub(crate) fn try_std_sort<T: Copy>(
    slice: &mut [T],
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> Result<(), OutOfRange> {
    if slice.is_empty() {
        return Ok(());
    }
    let lg = (usize::BITS - slice.len().leading_zeros() - 1) as usize;
    introsort_loop(slice, 0, slice.len(), lg * 2, comp)?;
    final_insertion_sort(slice, 0, slice.len(), comp)
}

/// `_S_threshold`.
// libstdc++ (GCC 13): bits/stl_algo.h:1848
const S_THRESHOLD: usize = 16;

// libstdc++ (GCC 13): bits/stl_algo.h:1916-1934
// (`__last - __first` is a signed difference: a cut past `__last`, which a comparison that is not
// a strict weak order can give, makes it negative and ends the loop.)
fn introsort_loop<T: Copy>(
    s: &mut [T],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> Result<(), OutOfRange> {
    while last > first && last - first > S_THRESHOLD {
        if depth_limit == 0 {
            partial_sort(s, first, last, last, comp);
            return Ok(());
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(s, first, last, comp)?;
        introsort_loop(s, cut, last, depth_limit, comp)?;
        last = cut;
    }
    Ok(())
}

// libstdc++ (GCC 13): bits/stl_algo.h:1891-1900
fn unguarded_partition_pivot<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> Result<usize, OutOfRange> {
    let mid = first + (last - first) / 2;
    move_median_to_first(s, first, first + 1, mid, last - 1, comp);
    unguarded_partition(s, first + 1, last, first, comp)
}

// libstdc++ (GCC 13): bits/stl_algo.h:81-102
fn move_median_to_first<T: Copy>(
    s: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    let chosen = if comp(&s[a], &s[b]) {
        if comp(&s[b], &s[c]) {
            b
        } else if comp(&s[a], &s[c]) {
            c
        } else {
            a
        }
    } else if comp(&s[a], &s[c]) {
        a
    } else if comp(&s[b], &s[c]) {
        c
    } else {
        b
    };
    s.swap(result, chosen);
}

// libstdc++ (GCC 13): bits/stl_algo.h:1867-1887
// (The first scan reads `*__first` until the comparison fails: past the slice's end it reads
// outside the vector, `OutOfRange`. The second scan stops at the pivot at the latest, since
// `comp(pivot, pivot)` is false for rpsbproc's comparison, NaN or not; past the slice's start it
// would read outside too.)
fn unguarded_partition<T: Copy>(
    s: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> Result<usize, OutOfRange> {
    loop {
        loop {
            if first >= s.len() {
                return Err(OutOfRange);
            }
            if !comp(&s[first], &s[pivot]) {
                break;
            }
            first += 1;
        }
        last = last.checked_sub(1).ok_or(OutOfRange)?;
        while comp(&s[pivot], &s[last]) {
            last = last.checked_sub(1).ok_or(OutOfRange)?;
        }
        if first >= last {
            return Ok(first);
        }
        s.swap(first, last);
        first += 1;
    }
}

// libstdc++ (GCC 13): bits/stl_algo.h:1902-1912
fn partial_sort<T: Copy>(
    s: &mut [T],
    first: usize,
    middle: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    heap_select(s, first, middle, last, comp);
    sort_heap(s, first, middle, comp);
}

// libstdc++ (GCC 13): bits/stl_algo.h:1628-1639
fn heap_select<T: Copy>(
    s: &mut [T],
    first: usize,
    middle: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    make_heap(s, first, middle, comp);
    for i in middle..last {
        if comp(&s[i], &s[first]) {
            pop_heap(s, first, middle, i, comp);
        }
    }
}

// libstdc++ (GCC 13): bits/stl_heap.h:338-361
fn make_heap<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    if last - first < 2 {
        return;
    }
    let len = last - first;
    let mut parent = (len - 2) / 2;
    loop {
        let value = s[first + parent];
        adjust_heap(s, first, parent, len, value, comp);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

// libstdc++ (GCC 13): bits/stl_heap.h:417-427
fn sort_heap<T: Copy>(
    s: &mut [T],
    first: usize,
    mut last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    while last - first > 1 {
        last -= 1;
        pop_heap(s, first, last, last, comp);
    }
}

// libstdc++ (GCC 13): bits/stl_heap.h:252-266
fn pop_heap<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    result: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    let value = s[result];
    s[result] = s[first];
    adjust_heap(s, first, 0, last - first, value, comp);
}

// libstdc++ (GCC 13): bits/stl_heap.h:221-250
fn adjust_heap<T: Copy>(
    s: &mut [T],
    first: usize,
    mut hole: usize,
    len: usize,
    value: T,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    let top = hole;
    let mut second = hole;
    // `(__len - 1) / 2` with `__len` at least 1 here (a pop of the last element has 1).
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if comp(&s[first + second], &s[first + second - 1]) {
            second -= 1;
        }
        s[first + hole] = s[first + second];
        hole = second;
    }
    if len & 1 == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        s[first + hole] = s[first + second - 1];
        hole = second - 1;
    }
    push_heap(s, first, hole, top, value, comp);
}

// libstdc++ (GCC 13): bits/stl_heap.h:132-147
// (`(__holeIndex - 1) / 2` of a hole at 0 is 0 in C++'s signed division; the loop then stops
// at `__holeIndex > __topIndex`.)
fn push_heap<T: Copy>(
    s: &mut [T],
    first: usize,
    mut hole: usize,
    top: usize,
    value: T,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    let mut parent = hole.saturating_sub(1) / 2;
    while hole > top && comp(&s[first + parent], &value) {
        s[first + hole] = s[first + parent];
        hole = parent;
        parent = hole.saturating_sub(1) / 2;
    }
    s[first + hole] = value;
}

// libstdc++ (GCC 13): bits/stl_algo.h:1852-1865
fn final_insertion_sort<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> Result<(), OutOfRange> {
    if last - first > S_THRESHOLD {
        insertion_sort(s, first, first + S_THRESHOLD, comp)?;
        // `__unguarded_insertion_sort` (stl_algo.h:1834-1842).
        for i in first + S_THRESHOLD..last {
            unguarded_linear_insert(s, i, comp)?;
        }
        Ok(())
    } else {
        insertion_sort(s, first, last, comp)
    }
}

// libstdc++ (GCC 13): bits/stl_algo.h:1810-1830
// (The walk of `__unguarded_linear_insert` stops at `__first` at the latest: `comp(val, *__first)`
// was false, and `*__first` is not moved before the walk reaches it.)
fn insertion_sort<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> Result<(), OutOfRange> {
    if first == last {
        return Ok(());
    }
    for i in first + 1..last {
        if comp(&s[i], &s[first]) {
            let value = s[i];
            s.copy_within(first..i, first + 1);
            s[first] = value;
        } else {
            unguarded_linear_insert(s, i, comp)?;
        }
    }
    Ok(())
}

// libstdc++ (GCC 13): bits/stl_algo.h:1790-1806
// (Before the slice's start `--__next` points outside the vector: `OutOfRange`, the moves made
// so far left in place as the error ends the sort.)
fn unguarded_linear_insert<T: Copy>(
    s: &mut [T],
    mut last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> Result<(), OutOfRange> {
    let value = s[last];
    let mut next = last - 1;
    while comp(&value, &s[next]) {
        s[last] = s[next];
        last = next;
        next = next.checked_sub(1).ok_or(OutOfRange)?;
    }
    s[last] = value;
    Ok(())
}

#[cfg(test)]
#[path = "std_sort_tests.rs"]
mod tests;
