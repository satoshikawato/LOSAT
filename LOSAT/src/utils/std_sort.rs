//! `std::sort` of GCC's C++ library (libstdc++ 13, the library NCBI's rpsbproc 0.5.1 binary is
//! built with, GCC 13.2): introsort with a median-of-three pivot, a heap sort where the
//! recursion gets too deep, and a final insertion sort. It is not stable: elements that the
//! comparison does not order end in the order this algorithm leaves them, which NCBI's output
//! then shows (the domain hits of equal E-value and bit score, the lookup tables of the CDD
//! data). Every step below is the library's, so that the order is the same for every input.
//!
//! The functions take a slice and a "less" comparison (`comp(a, b)`: `a` goes before `b`), the
//! `_Iter_comp_iter` of `std::sort(first, last, comp)`. The comparison must be a strict weak
//! order (libstdc++'s unguarded loops run past the ends otherwise; here the slice indexing
//! panics).

/// `std::sort(first, last, comp)`.
// libstdc++ (GCC 13): bits/stl_algo.h:1940-1951
// ```c++
//   template<typename _RandomAccessIterator, typename _Compare>
//     _GLIBCXX20_CONSTEXPR
//     inline void
//     __sort(_RandomAccessIterator __first, _RandomAccessIterator __last,
// 	   _Compare __comp)
//     {
//       if (__first != __last)
// 	{
// 	  std::__introsort_loop(__first, __last,
// 				std::__lg(__last - __first) * 2,
// 				__comp);
// 	  std::__final_insertion_sort(__first, __last, __comp);
// 	}
//     }
// ```
// libstdc++ (GCC 13): bits/stl_algobase.h:1513-1518
// ```c++
//   template<typename _Tp>
//     inline _GLIBCXX_CONSTEXPR _Tp
//     __lg(_Tp __n)
//     {
// #if __cplusplus >= 201402L
//       return std::__bit_width(make_unsigned_t<_Tp>(__n)) - 1;
// ```
pub(crate) fn std_sort<T: Copy>(slice: &mut [T], comp: &mut impl FnMut(&T, &T) -> bool) {
    if slice.is_empty() {
        return;
    }
    let lg = (usize::BITS - slice.len().leading_zeros() - 1) as usize;
    introsort_loop(slice, 0, slice.len(), lg * 2, comp);
    final_insertion_sort(slice, 0, slice.len(), comp);
}

/// `_S_threshold`.
// libstdc++ (GCC 13): bits/stl_algo.h:1848
// ```c++
//   enum { _S_threshold = 16 };
// ```
const S_THRESHOLD: usize = 16;

// libstdc++ (GCC 13): bits/stl_algo.h:1916-1934
// ```c++
//     __introsort_loop(_RandomAccessIterator __first,
// 		     _RandomAccessIterator __last,
// 		     _Size __depth_limit, _Compare __comp)
//     {
//       while (__last - __first > int(_S_threshold))
// 	{
// 	  if (__depth_limit == 0)
// 	    {
// 	      std::__partial_sort(__first, __last, __last, __comp);
// 	      return;
// 	    }
// 	  --__depth_limit;
// 	  _RandomAccessIterator __cut =
// 	    std::__unguarded_partition_pivot(__first, __last, __comp);
// 	  std::__introsort_loop(__cut, __last, __depth_limit, __comp);
// 	  __last = __cut;
// 	}
//     }
// ```
fn introsort_loop<T: Copy>(
    s: &mut [T],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    while last - first > S_THRESHOLD {
        if depth_limit == 0 {
            partial_sort(s, first, last, last, comp);
            return;
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(s, first, last, comp);
        introsort_loop(s, cut, last, depth_limit, comp);
        last = cut;
    }
}

// libstdc++ (GCC 13): bits/stl_algo.h:1891-1900
// ```c++
//     __unguarded_partition_pivot(_RandomAccessIterator __first,
// 				_RandomAccessIterator __last, _Compare __comp)
//     {
//       _RandomAccessIterator __mid = __first + (__last - __first) / 2;
//       std::__move_median_to_first(__first, __first + 1, __mid, __last - 1,
// 				  __comp);
//       return std::__unguarded_partition(__first + 1, __last, __first, __comp);
//     }
// ```
fn unguarded_partition_pivot<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(s, first, first + 1, mid, last - 1, comp);
    unguarded_partition(s, first + 1, last, first, comp)
}

// libstdc++ (GCC 13): bits/stl_algo.h:81-102
// ```c++
//     __move_median_to_first(_Iterator __result,_Iterator __a, _Iterator __b,
// 			   _Iterator __c, _Compare __comp)
//     {
//       if (__comp(__a, __b))
// 	{
// 	  if (__comp(__b, __c))
// 	    std::iter_swap(__result, __b);
// 	  else if (__comp(__a, __c))
// 	    std::iter_swap(__result, __c);
// 	  else
// 	    std::iter_swap(__result, __a);
// 	}
//       else if (__comp(__a, __c))
// 	std::iter_swap(__result, __a);
//       else if (__comp(__b, __c))
// 	std::iter_swap(__result, __c);
//       else
// 	std::iter_swap(__result, __b);
//     }
// ```
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
// ```c++
//     __unguarded_partition(_RandomAccessIterator __first,
// 			  _RandomAccessIterator __last,
// 			  _RandomAccessIterator __pivot, _Compare __comp)
//     {
//       while (true)
// 	{
// 	  while (__comp(__first, __pivot))
// 	    ++__first;
// 	  --__last;
// 	  while (__comp(__pivot, __last))
// 	    --__last;
// 	  if (!(__first < __last))
// 	    return __first;
// 	  std::iter_swap(__first, __last);
// 	  ++__first;
// 	}
//     }
// ```
fn unguarded_partition<T: Copy>(
    s: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) -> usize {
    loop {
        while comp(&s[first], &s[pivot]) {
            first += 1;
        }
        last -= 1;
        while comp(&s[pivot], &s[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        s.swap(first, last);
        first += 1;
    }
}

// libstdc++ (GCC 13): bits/stl_algo.h:1902-1912
// ```c++
//     __partial_sort(_RandomAccessIterator __first,
// 		   _RandomAccessIterator __middle,
// 		   _RandomAccessIterator __last,
// 		   _Compare __comp)
//     {
//       std::__heap_select(__first, __middle, __last, __comp);
//       std::__sort_heap(__first, __middle, __comp);
//     }
// ```
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
// ```c++
//     __heap_select(_RandomAccessIterator __first,
// 		  _RandomAccessIterator __middle,
// 		  _RandomAccessIterator __last, _Compare __comp)
//     {
//       std::__make_heap(__first, __middle, __comp);
//       for (_RandomAccessIterator __i = __middle; __i < __last; ++__i)
// 	if (__comp(__i, __first))
// 	  std::__pop_heap(__first, __middle, __i, __comp);
//     }
// ```
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
// ```c++
//     __make_heap(_RandomAccessIterator __first, _RandomAccessIterator __last,
// 		_Compare& __comp)
//     {
//       ...
//       if (__last - __first < 2)
// 	return;
//
//       const _DistanceType __len = __last - __first;
//       _DistanceType __parent = (__len - 2) / 2;
//       while (true)
// 	{
// 	  _ValueType __value = _GLIBCXX_MOVE(*(__first + __parent));
// 	  std::__adjust_heap(__first, __parent, __len, _GLIBCXX_MOVE(__value),
// 			     __comp);
// 	  if (__parent == 0)
// 	    return;
// 	  __parent--;
// 	}
//     }
// ```
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
// ```c++
//     __sort_heap(_RandomAccessIterator __first, _RandomAccessIterator __last,
// 		_Compare& __comp)
//     {
//       while (__last - __first > 1)
// 	{
// 	  --__last;
// 	  std::__pop_heap(__first, __last, __last, __comp);
// 	}
//     }
// ```
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
// ```c++
//     __pop_heap(_RandomAccessIterator __first, _RandomAccessIterator __last,
// 	       _RandomAccessIterator __result, _Compare& __comp)
//     {
//       ...
//       _ValueType __value = _GLIBCXX_MOVE(*__result);
//       *__result = _GLIBCXX_MOVE(*__first);
//       std::__adjust_heap(__first, _DistanceType(0),
// 			 _DistanceType(__last - __first),
// 			 _GLIBCXX_MOVE(__value), __comp);
//     }
// ```
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
// ```c++
//     __adjust_heap(_RandomAccessIterator __first, _Distance __holeIndex,
// 		  _Distance __len, _Tp __value, _Compare __comp)
//     {
//       const _Distance __topIndex = __holeIndex;
//       _Distance __secondChild = __holeIndex;
//       while (__secondChild < (__len - 1) / 2)
// 	{
// 	  __secondChild = 2 * (__secondChild + 1);
// 	  if (__comp(__first + __secondChild,
// 		     __first + (__secondChild - 1)))
// 	    __secondChild--;
// 	  *(__first + __holeIndex) = _GLIBCXX_MOVE(*(__first + __secondChild));
// 	  __holeIndex = __secondChild;
// 	}
//       if ((__len & 1) == 0 && __secondChild == (__len - 2) / 2)
// 	{
// 	  __secondChild = 2 * (__secondChild + 1);
// 	  *(__first + __holeIndex) = _GLIBCXX_MOVE(*(__first
// 						     + (__secondChild - 1)));
// 	  __holeIndex = __secondChild - 1;
// 	}
//       ...
//       std::__push_heap(__first, __holeIndex, __topIndex,
// 		       _GLIBCXX_MOVE(__value), __cmp);
//     }
// ```
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
// ```c++
//     __push_heap(_RandomAccessIterator __first,
// 		_Distance __holeIndex, _Distance __topIndex, _Tp __value,
// 		_Compare& __comp)
//     {
//       _Distance __parent = (__holeIndex - 1) / 2;
//       while (__holeIndex > __topIndex && __comp(__first + __parent, __value))
// 	{
// 	  *(__first + __holeIndex) = _GLIBCXX_MOVE(*(__first + __parent));
// 	  __holeIndex = __parent;
// 	  __parent = (__holeIndex - 1) / 2;
// 	}
//       *(__first + __holeIndex) = _GLIBCXX_MOVE(__value);
//     }
// ```
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
// ```c++
//     __final_insertion_sort(_RandomAccessIterator __first,
// 			   _RandomAccessIterator __last, _Compare __comp)
//     {
//       if (__last - __first > int(_S_threshold))
// 	{
// 	  std::__insertion_sort(__first, __first + int(_S_threshold), __comp);
// 	  std::__unguarded_insertion_sort(__first + int(_S_threshold), __last,
// 					  __comp);
// 	}
//       else
// 	std::__insertion_sort(__first, __last, __comp);
//     }
// ```
fn final_insertion_sort<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    if last - first > S_THRESHOLD {
        insertion_sort(s, first, first + S_THRESHOLD, comp);
        // `__unguarded_insertion_sort` (stl_algo.h:1834-1842).
        for i in first + S_THRESHOLD..last {
            unguarded_linear_insert(s, i, comp);
        }
    } else {
        insertion_sort(s, first, last, comp);
    }
}

// libstdc++ (GCC 13): bits/stl_algo.h:1810-1830
// ```c++
//     __insertion_sort(_RandomAccessIterator __first,
// 		     _RandomAccessIterator __last, _Compare __comp)
//     {
//       if (__first == __last) return;
//
//       for (_RandomAccessIterator __i = __first + 1; __i != __last; ++__i)
// 	{
// 	  if (__comp(__i, __first))
// 	    {
// 	      typename iterator_traits<_RandomAccessIterator>::value_type
// 		__val = _GLIBCXX_MOVE(*__i);
// 	      _GLIBCXX_MOVE_BACKWARD3(__first, __i, __i + 1);
// 	      *__first = _GLIBCXX_MOVE(__val);
// 	    }
// 	  else
// 	    std::__unguarded_linear_insert(__i,
// 				__gnu_cxx::__ops::__val_comp_iter(__comp));
// 	}
//     }
// ```
fn insertion_sort<T: Copy>(
    s: &mut [T],
    first: usize,
    last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        if comp(&s[i], &s[first]) {
            let value = s[i];
            s.copy_within(first..i, first + 1);
            s[first] = value;
        } else {
            unguarded_linear_insert(s, i, comp);
        }
    }
}

// libstdc++ (GCC 13): bits/stl_algo.h:1790-1806
// ```c++
//     __unguarded_linear_insert(_RandomAccessIterator __last,
// 			      _Compare __comp)
//     {
//       typename iterator_traits<_RandomAccessIterator>::value_type
// 	__val = _GLIBCXX_MOVE(*__last);
//       _RandomAccessIterator __next = __last;
//       --__next;
//       while (__comp(__val, __next))
// 	{
// 	  *__last = _GLIBCXX_MOVE(*__next);
// 	  __last = __next;
// 	  --__next;
// 	}
//       *__last = _GLIBCXX_MOVE(__val);
//     }
// ```
fn unguarded_linear_insert<T: Copy>(
    s: &mut [T],
    mut last: usize,
    comp: &mut impl FnMut(&T, &T) -> bool,
) {
    let value = s[last];
    let mut next = last - 1;
    while comp(&value, &s[next]) {
        s[last] = s[next];
        last = next;
        next -= 1;
    }
    s[last] = value;
}

#[cfg(test)]
#[path = "std_sort_tests.rs"]
mod tests;
