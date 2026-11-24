//! The searches of a byte string that LOSAT's ports share: the first byte that a test holds
//! for ([`find_first`]), the first of one byte ([`find_byte`], C's `memchr`), the first
//! occurrence of a byte string ([`find`], `string::find`) and the mask of one byte in a block of
//! 16 ([`equal_bytes_16`], [`gather_high_bits`]). makeblastdb's FASTA reader, its titles
//! and its header writer search every byte of every line, so the search tests the bytes a
//! window at a time.
//!
//! No NCBI counterpart: NCBI's code calls the C and C++ libraries for these searches, or tests
//! one byte at a time (dustmasker's test of each letter for `N`); LOSAT keeps what they return
//! (the index of the first match, the bytes that match), not how they search.

/// The bytes of a window of [`find_first`]: a window's indexes fit a byte.
const WINDOW: usize = 64;

/// The bytes of a block that [`find_first`] tests at once in a haystack shorter than a window.
const BLOCK: usize = 32;

/// The index of the first byte of `window` (at most 255 bytes) that `test` holds for, or
/// `u8::MAX` where it holds for none: the least index of a match, in a loop without an early
/// exit, which the compiler's loop vectorizer turns into vector compares and a vector minimum.
///
/// The form is deliberate: the same minimum written as a fold over `enumerate()` is not
/// vectorized by rustc 1.92 (8-13 instructions a byte against about 1), so keep the `for` loop.
#[inline(always)]
fn first_in(window: &[u8], test: impl Fn(u8) -> bool) -> u8 {
    debug_assert!(window.len() < usize::from(u8::MAX));
    let mut least = u8::MAX;
    for (index, &b) in window.iter().enumerate() {
        least = least.min(if test(b) { index as u8 } else { u8::MAX });
    }
    least
}

/// Whether `test` holds for a byte of `window`, in a loop without an early exit (vectorized as
/// [`first_in`], with an OR in place of the minimum). The OR is of bytes: a fold of `bool`s is
/// not vectorized by rustc 1.92.
#[inline(always)]
fn holds_in(window: &[u8], test: impl Fn(u8) -> bool) -> bool {
    window.iter().fold(0u8, |any, &b| any | u8::from(test(b))) != 0
}

/// The index of the first byte of `haystack` that `test` holds for.
///
/// The first [`WINDOW`] bytes are searched at once ([`first_in`]); the windows after them are
/// tested first ([`holds_in`]) and searched only where `test` holds, the last window ending at
/// the end of `haystack` and overlapping the one before it (where `test` holds for no byte). A
/// haystack shorter than a window is tested as two blocks of [`BLOCK`] bytes that overlap, by a
/// fold over each, and searched only where `test` holds; one shorter than a block is searched
/// alone. The forms differ only in how many instructions they take (a short line is mostly its
/// own haystack; a line end is searched for in a whole buffer).
pub(crate) fn find_first(haystack: &[u8], test: impl Fn(u8) -> bool + Copy) -> Option<usize> {
    let len = haystack.len();
    let found = |least: u8| (least != u8::MAX).then_some(usize::from(least));
    if len < WINDOW {
        if len >= BLOCK {
            let first = &haystack.as_chunks::<BLOCK>().0[0];
            let last = &haystack[len - BLOCK..].as_chunks::<BLOCK>().0[0];
            if !holds_in(first, test) && !holds_in(last, test) {
                return None;
            }
        }
        return found(first_in(haystack, test));
    }
    if let Some(index) = found(first_in(&haystack[..WINDOW], test)) {
        return Some(index);
    }
    let mut start = WINDOW;
    while start < len {
        let base = start.min(len - WINDOW);
        let window = &haystack[base..base + WINDOW];
        if holds_in(window, test) {
            return found(first_in(window, test)).map(|index| base + index);
        }
        start = base + WINDOW;
    }
    None
}

/// The index of the first `byte` of `haystack` (C's `memchr`).
pub(crate) fn find_byte(haystack: &[u8], byte: u8) -> Option<usize> {
    find_first(haystack, |b| b == byte)
}

/// `string::find(needle)` of `haystack`: the index of the first occurrence of `needle` (0 for
/// an empty one). The first byte of `needle` is searched with [`find_byte`], and the rest is
/// compared where it is found; one search for every needle, whatever its length.
pub(crate) fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    let Some((&head, rest)) = needle.split_first() else {
        return Some(0);
    };
    let mut from = 0;
    while let Some(offset) = find_byte(&haystack[from..], head) {
        let at = from + offset;
        if haystack[at + 1..].starts_with(rest) {
            return Some(at);
        }
        from = at + 1;
    }
    None
}

/// The bytes of `block` that are `byte`, as the high bits of the bytes of two words: bit
/// `8 k + 7` of the first word for byte `k`, of the second for byte `k + 8`; every other bit is
/// clear. One test of a word tells whether its 8 bytes hold `byte` at all (`!= 0`) or are all
/// `byte`; [`gather_high_bits`] makes the mask of a word. A byte is tested without a branch: it
/// is `byte` when its XOR with `byte` is 0, and the high bit of `((x & 0x7f) + 0x7f) | x` is set
/// when the byte `x` is not 0 (0x7f + 0x7f = 0xfe, so no sum leaves its byte).
pub(crate) fn equal_bytes_16(block: &[u8; 16], byte: u8) -> (u64, u64) {
    const LOW7: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    let spread = u64::from_ne_bytes([byte; 8]);
    let equal = |bytes: &[u8]| {
        let x = u64::from_le_bytes(bytes.try_into().expect("8 bytes")) ^ spread;
        !(((x & LOW7) + LOW7) | x | LOW7)
    };
    let (low, high) = block.split_at(8);
    (equal(low), equal(high))
}

/// The high bits of the 8 bytes of `bits` (as [`equal_bytes_16`] gives them) as the 8 low bits
/// of the result, byte `k`'s at bit `k`. The high bits are moved to bits `8 k`, and the product
/// with `0x0102_0408_1020_4080` puts bit `8 k` at bit `56 + k`; its other terms land on bits that
/// are all distinct and below 56 or above 63, so no sum carries into bits 56 to 63.
pub(crate) fn gather_high_bits(bits: u64) -> u32 {
    ((bits >> 7).wrapping_mul(0x0102_0408_1020_4080) >> 56) as u32
}

#[cfg(test)]
#[path = "byte_search_tests.rs"]
pub(crate) mod tests;
