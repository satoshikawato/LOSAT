//! The growth of libstdc++'s `std::string` and `std::vector` in NCBI's rpsbproc's motif check:
//! the length and capacity of `GetMinimalXMap`'s string (`CxxString`) and of the vector of match
//! records (`Live::grow`, `push`, `pop`), and how NCBI ends where they cannot grow (`MapEnd`).
//! `prosite` re-exports `CxxString` and `MapEnd`.

use super::{Live, MatchError, HOST_DEPENDENT_BYTES, USER_SPACE_BYTES};

/// How NCBI's `GetMinimalXMap` (or the `push_back` of the text into its string) ends instead of
/// returning.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MapEnd {
    /// An exception of libstdc++, by its `what()`: `LENGTH_ERROR` or `BAD_ALLOC`.
    Throw(&'static str),
    /// An allocation below `USER_SPACE_BYTES` while NCBI holds this many bytes for the string
    /// (its old buffer and the new one), `HOST_DEPENDENT_BYTES` or more.
    HostDependent(u64),
}

/// `what()` of the `std::length_error` of `_M_check_length` in `append` (the length passes
/// `max_size()`).
pub(crate) const LENGTH_ERROR: &str = "basic_string::_M_replace_aux";

/// `what()` of `std::bad_alloc` (an allocation of `USER_SPACE_BYTES` or more).
pub(crate) const BAD_ALLOC: &str = "std::bad_alloc";

/// `operator new(bytes)` of NCBI's for a string while `held` bytes of its old buffer live: one
/// allocation of `USER_SPACE_BYTES` or more throws on every host (DECISIONS 37); below it, a total
/// of `HOST_DEPENDENT_BYTES` or more is the host's (DECISIONS 40).
fn allocate_string(held: u64, bytes: u64) -> Result<(), MapEnd> {
    if bytes >= USER_SPACE_BYTES {
        Err(MapEnd::Throw(BAD_ALLOC))
    } else if held.saturating_add(bytes) >= HOST_DEPENDENT_BYTES {
        Err(MapEnd::HostDependent(held.saturating_add(bytes)))
    } else {
        Ok(())
    }
}

/// The length and capacity of a `std::string` of libstdc++ (the C++11 ABI) whose bytes LOSAT
/// does not need: `GetMinimalXMap`'s `minMap` and the text that `MotifCheck` builds in it.
/// `max_size()` is 2^63 - 1: NCBI's rpsbproc 0.5.1 is C++20 code built by GCC 13 (the
/// distributed binary 13.2.0, Bioconda's 13.4.0; review 2 oracle O1), and in C++20 libstdc++ does
/// not declare `basic_string<char>` an extern template, so `append` and its `max_size()` are
/// compiled into the binary from the headers below, whatever `libstdc++.so` it loads (review 2
/// oracle O2: map lengths of 2^62 to 2^63 - 1 throw `std::bad_alloc`; R20, 2^63, throws
/// `length_error`).
// libstdc++ (GCC 13): bits/basic_string.tcc:962-969
// libstdc++ (GCC 13): bits/alloc_traits.h:571-578
// (`basic_string::max_size()` is `(_Alloc_traits::max_size(...) - 1) / 2`, below: (2^64 - 2) / 2.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CxxString {
    pub(crate) len: u64,
    /// `capacity()`: the bytes `0..=capacity` are NCBI's buffer (`_M_local_buf[16]` while it
    /// is 15).
    pub(crate) capacity: u64,
}

impl CxxString {
    const LOCAL_CAPACITY: u64 = 15;
    const MAX_SIZE: u64 = (1 << 63) - 1;

    /// An empty string (`string minMap(k_strEmptyString)`).
    pub(crate) fn new() -> Self {
        CxxString {
            len: 0,
            capacity: Self::LOCAL_CAPACITY,
        }
    }

    // libstdc++ (GCC 13): bits/basic_string.tcc:139-160
    // (The callers never ask for more than `max_size()`: `append` checks the length first.) The
    // callers (`reserve`, `_M_mutate`) free the old buffer after the copy into the new one, so
    // both live at the allocation.
    fn create(&mut self, mut capacity: u64) -> Result<(), MapEnd> {
        if capacity > self.capacity && capacity < 2 * self.capacity {
            capacity = (2 * self.capacity).min(Self::MAX_SIZE);
        }
        allocate_string(self.heap_bytes(), capacity + 1)?;
        self.capacity = capacity;
        Ok(())
    }

    /// The bytes of the string's buffer on the heap: none while it is the local buffer.
    pub(crate) fn heap_bytes(&self) -> u64 {
        if self.capacity > Self::LOCAL_CAPACITY {
            self.capacity + 1
        } else {
            0
        }
    }

    // libstdc++ (GCC 13): bits/basic_string.tcc:305-320
    // (GCC before 11 also shrinks; `minMap` is empty and local when `GetMinimalXMap` reserves,
    // so both give this.)
    pub(super) fn reserve(&mut self, res: u64) -> Result<(), MapEnd> {
        if res <= self.capacity {
            return Ok(());
        }
        self.create(res)
    }

    // libstdc++ (GCC 13): bits/basic_string.tcc:451-468
    // libstdc++ (GCC 13): bits/basic_string.h:397-401,1084-1085,1487-1488
    // libstdc++ (GCC 13): bits/basic_string.tcc:326-345
    /// `append(n, c)`.
    pub(crate) fn append(&mut self, n: u64) -> Result<(), MapEnd> {
        if Self::MAX_SIZE - self.len < n {
            return Err(MapEnd::Throw(LENGTH_ERROR));
        }
        let len = self.len + n;
        if len > self.capacity {
            self.create(len)?;
        }
        self.len = len;
        Ok(())
    }

    // libstdc++ (GCC 13): bits/basic_string.h:1561-1568
    pub(crate) fn push_back(&mut self) -> Result<(), MapEnd> {
        if self.len + 1 > self.capacity {
            self.create(self.len + 1)?;
        }
        self.len += 1;
        Ok(())
    }

    /// `clear()`: the capacity stays.
    pub(crate) fn clear(&mut self) {
        self.len = 0;
    }
}

impl Live {
    /// `push_back` of records into `curr.flex` until it holds `len`: libstdc++ makes room for
    /// twice the records it holds when it is full (one when it is empty), and frees the old
    /// buffer after the move into the new one.
    // libstdc++ (GCC 13): bits/stl_vector.h:1281-1293,1896-1903
    // libstdc++ (GCC 13): bits/vector.tcc:453-459,517-518
    // (`max_size()` is 2^60 - 1 records; 1 GiB stops the growth long before.)
    pub(super) fn grow(&mut self, len: u64) -> Result<(), MatchError> {
        while self.capacity < len {
            let new = if self.capacity == 0 {
                1
            } else {
                2 * self.capacity
            };
            self.hold(self.capacity + new + self.stacked)?;
            self.capacity = new;
        }
        Ok(())
    }

    /// `stkFlexStack.push(curr)`: the copy constructor allocates the `len` records it copies.
    // libstdc++ (GCC 13): bits/stl_vector.h:601-609
    pub(super) fn push(&mut self, len: u64) -> Result<(), MatchError> {
        self.hold(self.capacity + self.stacked + len)?;
        self.stacked += len;
        Ok(())
    }

    /// Taking the top of `stkFlexStack` into `curr` and popping it, a copy of `len` records: the
    /// assignment allocates only when they do not fit, while the old buffer and the copy on the
    /// stack live. (NCBI's `curr.flex` never shrinks and each copy was made from it, so its
    /// copies always fit; this follows libstdc++ all the same.)
    // libstdc++ (GCC 13): bits/vector.tcc:235-246
    pub(super) fn pop(&mut self, len: u64) -> Result<(), MatchError> {
        if len > self.capacity {
            self.hold(self.capacity + len + self.stacked)?;
            self.capacity = len;
        }
        self.stacked -= len;
        Ok(())
    }
}
