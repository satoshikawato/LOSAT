//! The growth of the `std::string` and the `std::vector<size_t>` of NCBI's rpsbproc motif check
//! (`CProSite::GetMinimalXMap`, `CProSite::Match`), by length and capacity, without their bytes.
//!
//! The model is written from observations of the real objects, not from the C++ library's
//! source. The tables are in `LOSAT/tests/fixtures/rpsbproc_growth/` (`README.md`: the compiler,
//! the library version, the observation programs and the dates):
//! - `string_ops.tsv`: the length, the capacity and the most heap bytes live at once after each
//!   operation of `reserve(n)`, `append(n, c)`, `push_back(c)` and `clear()` on the string, or the
//!   `what()` of the exception it threw;
//! - `vector_ops.tsv`: the same for the vector of a position's match records and the copies of it
//!   on the match's stack;
//! - `binary_allocs.tsv`: the allocations that NCBI's rpsbproc 0.5.1 made for edited motifs.
//!
//! What the tables show (the C++ standard leaves these unspecified; they are what the binary
//! does):
//! - A string holds 15 bytes in the object itself (capacity 15, no heap buffer). When a call
//!   needs a capacity above the current one, the new capacity is the larger of the capacity
//!   asked for and twice the current one, and the buffer on the heap is that capacity + 1 bytes.
//!   The old buffer lives until the bytes are copied (both are live at once). `clear` keeps the
//!   capacity; `reserve` never lowers it.
//! - `append(n, c)` throws `std::length_error` (`LENGTH_ERROR`) when the length would pass
//!   `max_size()` = 2^63 - 1 (the binary is built as C++20; as C++17 it is 2^62 - 1), `reserve`
//!   throws one with another text for a capacity above it, and an allocation that the host
//!   cannot give throws `std::bad_alloc`. A string unchanged by a throw stays as it was.
//! - A vector of 8-byte records grows by doubling from 0 (1, 2, 4, ...), the new buffer live
//!   beside the old one. `clear` and `pop_back` keep the capacity. A copy made for the stack has
//!   exactly its length as capacity (no buffer when empty). Taking a copy back (`curr = top`)
//!   keeps the capacity when the copy fits; a copy that does not fit gets a buffer of exactly
//!   its length, allocated while the old buffer and the copy on the stack are live.
//!
//! The allocation limits are the caller's (`prosite.rs`): an allocation of `USER_SPACE_BYTES` or
//! more fails on every host (`std::bad_alloc`); what NCBI holds at once of `HOST_DEPENDENT_BYTES`
//! or more depends on the host (refused).

use super::{Live, MatchError, HOST_DEPENDENT_BYTES, USER_SPACE_BYTES};

/// `what()` of the `std::length_error` that NCBI's `append(n, c)` throws for a length above
/// `max_size()` (rpsbproc prints it as `Application's execution failed: ...`).
pub(crate) const LENGTH_ERROR: &str = "basic_string::_M_replace_aux";
/// `what()` of `std::bad_alloc`.
pub(crate) const BAD_ALLOC: &str = "std::bad_alloc";
/// `what()` of the `std::length_error` that `reserve` throws for a capacity above `max_size()`, as
/// observed with the fixture's own program (`programs/strprobe.cpp`, rows of `string_ops.tsv`), not
/// with rpsbproc: rpsbproc reserves twice a pattern's length and never reaches it.
const RESERVE_LENGTH_ERROR: &str = "basic_string::_M_create";

/// `max_size()` of NCBI's string (the build is C++20): 2^63 - 1.
const MAX_SIZE: u64 = (1 << 63) - 1;
/// The characters a string holds in the object itself.
const LOCAL_CAPACITY: u64 = 15;

/// How NCBI's rpsbproc ends instead of returning a minimal map or the text of a site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MapEnd {
    /// A C++ exception by its `what()` (`LENGTH_ERROR` or `BAD_ALLOC`).
    Throw(&'static str),
    /// NCBI would hold this many bytes at once, `HOST_DEPENDENT_BYTES` or more: the result
    /// depends on the host (refused).
    HostDependent(u64),
}

/// The length and the capacity of a `std::string` of NCBI's rpsbproc 0.5.1, without its bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CxxString {
    pub(crate) len: u64,
    pub(crate) capacity: u64,
}

impl CxxString {
    /// An empty string.
    pub(crate) fn new() -> Self {
        CxxString {
            len: 0,
            capacity: LOCAL_CAPACITY,
        }
    }

    /// The bytes of its buffer on the heap: 0 while the string fits in the object.
    pub(crate) fn heap_bytes(&self) -> u64 {
        if self.capacity > LOCAL_CAPACITY {
            self.capacity.saturating_add(1)
        } else {
            0
        }
    }

    /// `reserve(res)`.
    pub(super) fn reserve(&mut self, res: u64) -> Result<(), MapEnd> {
        if res > MAX_SIZE {
            return Err(MapEnd::Throw(RESERVE_LENGTH_ERROR));
        }
        if res > self.capacity {
            self.regrow(res)?;
        }
        Ok(())
    }

    /// `append(n, c)`.
    pub(crate) fn append(&mut self, n: u64) -> Result<(), MapEnd> {
        if n > MAX_SIZE - self.len {
            return Err(MapEnd::Throw(LENGTH_ERROR));
        }
        let len = self.len + n;
        if len > self.capacity {
            self.regrow(len)?;
        }
        self.len = len;
        Ok(())
    }

    /// `push_back(c)`.
    pub(crate) fn push_back(&mut self) -> Result<(), MapEnd> {
        self.append(1)
    }

    /// `clear()`: the capacity stays.
    pub(crate) fn clear(&mut self) {
        self.len = 0;
    }

    /// A new buffer for a capacity of `wanted` (above the current one): the larger of that and
    /// twice the current capacity, with one more byte for the terminating 0.
    fn regrow(&mut self, wanted: u64) -> Result<(), MapEnd> {
        let capacity = wanted.max(self.capacity.saturating_mul(2)).min(MAX_SIZE);
        let buffer = u128::from(capacity) + 1;
        if buffer >= u128::from(USER_SPACE_BYTES) {
            return Err(MapEnd::Throw(BAD_ALLOC));
        }
        // The old buffer is released after the bytes are copied.
        let held = u128::from(self.heap_bytes()) + buffer;
        if held >= u128::from(HOST_DEPENDENT_BYTES) {
            return Err(MapEnd::HostDependent(held as u64));
        }
        #[cfg(test)]
        events::string_grows(self.heap_bytes(), buffer as u64);
        self.capacity = capacity;
        Ok(())
    }
}

/// NCBI's vector of match records `curr.flex` (`Live::capacity`) and the copies of it on the
/// stack (`Live::stacked`). Each allocation is checked against what lives at that moment
/// (`Live::hold`).
impl Live {
    /// `curr.flex` holds `len` records (pushed one at a time): its capacity doubles from 0 until
    /// it holds them.
    pub(super) fn grow(&mut self, len: u64) -> Result<(), MatchError> {
        while self.capacity < len {
            let capacity = if self.capacity == 0 {
                1
            } else {
                self.capacity.saturating_mul(2)
            };
            // The new buffer is live beside the old one, and the copies on the stack.
            self.hold(
                self.stacked
                    .saturating_add(self.capacity)
                    .saturating_add(capacity),
            )?;
            #[cfg(test)]
            events::records_grow(self.capacity, capacity);
            self.capacity = capacity;
        }
        Ok(())
    }

    /// A copy of `len` records is pushed onto the stack: a buffer of exactly `len` records (none
    /// when empty) beside the capacity of `curr.flex` and the copies already there.
    pub(super) fn push(&mut self, len: u64) -> Result<(), MatchError> {
        if len > 0 {
            self.hold(
                self.stacked
                    .saturating_add(self.capacity)
                    .saturating_add(len),
            )?;
            #[cfg(test)]
            events::records_copy(len);
        }
        self.stacked = self.stacked.saturating_add(len);
        Ok(())
    }

    /// A copy of `len` records is taken back from the top of the stack into `curr.flex`, and
    /// removed from the stack: the capacity stays when the copy fits; otherwise `curr.flex`
    /// gets a buffer of exactly `len` records while its old one and the copy on the stack live.
    pub(super) fn pop(&mut self, len: u64) -> Result<(), MatchError> {
        if len > self.capacity {
            self.hold(
                self.stacked
                    .saturating_add(self.capacity)
                    .saturating_add(len),
            )?;
            #[cfg(test)]
            events::records_copy(len);
            self.capacity = len;
        }
        self.stacked = self.stacked.saturating_sub(len);
        Ok(())
    }
}

/// The allocations that the model makes, in order, for the tests that compare them with the
/// allocations of NCBI's binary (`binary_allocs.tsv`): `m<bytes>` for a buffer, `r<bytes>` for a
/// buffer of match records.
#[cfg(test)]
#[path = "prosite_growth_events_tests.rs"]
pub(crate) mod events;

#[cfg(test)]
#[path = "prosite_growth_tests.rs"]
mod tests;
