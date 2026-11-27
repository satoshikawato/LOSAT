//! The searches of a byte string that LOSAT's ports share: the first byte that a test holds
//! for ([`find_first`]), the first of one byte ([`find_byte`], C's `memchr`), the first of two
//! bytes ([`find_two_bytes`]), the first occurrence of a byte string ([`find`], `string::find`)
//! and the mask of one byte in a block of 16 ([`equal_bytes_16`], [`gather_high_bits`]).
//! makeblastdb's FASTA reader, its titles and its header writer search every byte of every line,
//! so the search tests the bytes a block or window at a time while keeping near matches
//! proportional to their distance.
//!
//! No NCBI counterpart: NCBI's code calls the C and C++ libraries for these searches, or tests
//! one byte at a time (dustmasker's test of each letter for `N`); LOSAT keeps what they return
//! (the index of the first match, the bytes that match), not how they search.

/// The bytes of a window of [`find_first`]: a window's indexes fit a byte.
const WINDOW: usize = 64;

/// Three 16-byte vector lanes, tested at once for medium haystacks of 48 to 63 bytes.
const MID: usize = 48;

/// Two 16-byte vector lanes, tested at once for medium haystacks of 32 to 47 bytes.
const BLOCK: usize = 32;

/// One 16-byte vector lane, tested at once for short and head blocks.
const SMALL: usize = 16;

/// Half of a 16-byte lane, used for 8-to-15-byte haystacks.
const HALF: usize = 8;

/// Leading bytes tested directly in a haystack of at least [`WINDOW`] bytes before the 16-byte
/// head block, so very near matches (short lines, dense hits) return without scanning a full
/// vector block.
const PREFIX: usize = 5;

/// The index of the first byte of `window` (64 bytes) that `test` holds for, or `u8::MAX` where
/// it holds for none: the least index of a match, in a loop without an early exit, which the
/// compiler's loop vectorizer turns into vector compares and a vector minimum.
///
/// The form is deliberate: the same minimum written as a fold over `enumerate()` is not
/// vectorized by rustc 1.92 (8-13 instructions a byte against about 1), so keep the `for` loop.
#[inline(always)]
fn first_in(window: &[u8; WINDOW], test: impl Fn(u8) -> bool) -> u8 {
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
fn holds_in<const N: usize>(window: &[u8; N], test: impl Fn(u8) -> bool) -> bool {
    window.iter().fold(0u8, |any, &b| any | u8::from(test(b))) != 0
}

#[inline(always)]
fn found(least: u8) -> Option<usize> {
    (least != u8::MAX).then_some(usize::from(least))
}

/// Searches `haystack` (length `N..N + SMALL`, byte 0 already ruled out) with an `N`-byte
/// leading check and, when longer than `N`, a trailing `SMALL`-byte overlapping check.
#[inline(always)]
fn find_first_split<const N: usize>(
    haystack: &[u8],
    test: impl Fn(u8) -> bool + Copy,
) -> Option<usize> {
    let first: &[u8; N] = haystack[..N].try_into().unwrap();
    if holds_in(first, test) {
        return first[1..].iter().position(|&b| test(b)).map(|i| 1 + i);
    }
    if haystack.len() == N {
        return None;
    }
    let last: &[u8; SMALL] = haystack[haystack.len() - SMALL..].try_into().unwrap();
    if holds_in(last, test) {
        return haystack[N..].iter().position(|&b| test(b)).map(|i| N + i);
    }
    None
}

/// Searches a haystack of 16 to 63 bytes whose first byte is already known not to match.
#[inline(never)]
fn find_first_medium(haystack: &[u8], test: impl Fn(u8) -> bool + Copy) -> Option<usize> {
    let len = haystack.len();
    if len >= MID {
        find_first_split::<MID>(haystack, test)
    } else if len >= BLOCK {
        find_first_split::<BLOCK>(haystack, test)
    } else {
        find_first_split::<SMALL>(haystack, test)
    }
}

/// Searches a haystack of at least [`WINDOW`] bytes after bytes `0..skip` have been ruled out.
#[inline(never)]
fn find_first_windows_from(
    haystack: &[u8],
    skip: usize,
    test: impl Fn(u8) -> bool + Copy,
) -> Option<usize> {
    let len = haystack.len();
    let first_base = skip.min(len - WINDOW);
    let first_win: &[u8; WINDOW] = haystack[first_base..first_base + WINDOW]
        .try_into()
        .unwrap();
    if let Some(index) = found(first_in(first_win, test)) {
        return Some(first_base + index);
    }
    let mut start = first_base + WINDOW;
    while start < len {
        let base = start.min(len - WINDOW);
        let window: &[u8; WINDOW] = haystack[base..base + WINDOW].try_into().unwrap();
        if holds_in(window, test) {
            return found(first_in(window, test)).map(|index| base + index);
        }
        start = base + WINDOW;
    }
    None
}

/// The bytes of an 8-byte block that are `byte`, as the high bit of each byte (all other bits
/// clear). A byte is tested without a branch: it equals `byte` when its XOR with `byte` is 0,
/// and the high bit of `((x & 0x7f) + 0x7f) | x` is set when byte `x` is non-zero (`0x7f + 0x7f
/// = 0xfe`, so no sum carries into the next byte).
#[inline(always)]
fn equal_bytes_8(bytes: &[u8; HALF], byte: u8) -> u64 {
    const LOW7: u64 = 0x7f7f_7f7f_7f7f_7f7f;
    let spread = u64::from_ne_bytes([byte; HALF]);
    let x = u64::from_le_bytes(*bytes) ^ spread;
    !(((x & LOW7) + LOW7) | x | LOW7)
}

/// The byte index (`0..8`) of the lowest set high bit in a non-zero mask from [`equal_bytes_8`].
#[inline(always)]
fn first_high_bit(mask: u64) -> usize {
    (mask.trailing_zeros() >> 3) as usize
}

/// Searches `haystack` with 8-byte word masks (`hit_at(offset)`) for short haystacks and the
/// first [`SMALL`] bytes, then [`find_first_windows_from`] for long haystacks.
#[inline(always)]
fn find_mask_8(
    haystack: &[u8],
    hit_at: impl Fn(usize) -> u64 + Copy,
    test: impl Fn(u8) -> bool + Copy,
) -> Option<usize> {
    let len = haystack.len();
    if len < HALF {
        return haystack.iter().position(|&b| test(b));
    }
    let h0 = hit_at(0);
    if h0 != 0 {
        return Some(first_high_bit(h0));
    }
    let check_tail = || {
        let tail = len - HALF;
        let ht = hit_at(tail);
        (ht != 0).then(|| tail + first_high_bit(ht))
    };
    if len < SMALL {
        return check_tail();
    }
    let h1 = hit_at(HALF);
    if h1 != 0 {
        return Some(HALF + first_high_bit(h1));
    }
    if len < WINDOW {
        let mut off = SMALL;
        while off + HALF <= len {
            let h = hit_at(off);
            if h != 0 {
                return Some(off + first_high_bit(h));
            }
            off += HALF;
        }
        return check_tail();
    }
    find_first_windows_from(haystack, SMALL, test)
}

/// The index of the first byte of `haystack` that `test` holds for.
///
/// Haystacks shorter than [`SMALL`] are scanned directly. Longer haystacks check byte 0 first so
/// dense hits return immediately, then widen progressively: 16-to-63-byte haystacks use
/// [`find_first_medium`], and haystacks of [`WINDOW`] bytes or more check bytes `1..PREFIX`, a
/// [`SMALL`]-byte head block, and then [`WINDOW`]-byte windows in [`find_first_windows_from`]
/// (the first window by [`first_in`], each later window by [`holds_in`] and [`first_in`] only
/// where `test` holds, with the last window overlapping the one before it).
#[inline(always)]
pub(crate) fn find_first(haystack: &[u8], test: impl Fn(u8) -> bool + Copy) -> Option<usize> {
    let len = haystack.len();
    if len < SMALL {
        return haystack.iter().position(|&b| test(b));
    }
    if test(haystack[0]) {
        return Some(0);
    }
    if len < WINDOW {
        return find_first_medium(haystack, test);
    }
    if test(haystack[1]) {
        return Some(1);
    }
    if test(haystack[2]) {
        return Some(2);
    }
    if test(haystack[3]) {
        return Some(3);
    }
    if test(haystack[4]) {
        return Some(4);
    }
    let head: &[u8; SMALL] = haystack[PREFIX..PREFIX + SMALL].try_into().unwrap();
    if holds_in(head, test) {
        return head.iter().position(|&b| test(b)).map(|i| PREFIX + i);
    }
    find_first_windows_from(haystack, PREFIX + SMALL, test)
}

/// The index of the first `byte` of `haystack` (C's `memchr`).
#[inline(always)]
pub(crate) fn find_byte(haystack: &[u8], byte: u8) -> Option<usize> {
    find_mask_8(
        haystack,
        |offset| equal_bytes_8(haystack[offset..offset + HALF].try_into().unwrap(), byte),
        |b| b == byte,
    )
}

/// The index of the first byte of `haystack` that equals `first` or `second`.
#[inline(always)]
pub(crate) fn find_two_bytes(haystack: &[u8], first: u8, second: u8) -> Option<usize> {
    find_mask_8(
        haystack,
        |offset| {
            let word: &[u8; HALF] = haystack[offset..offset + HALF].try_into().unwrap();
            equal_bytes_8(word, first) | equal_bytes_8(word, second)
        },
        |b| b == first || b == second,
    )
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
    let (low, high) = block.split_at(HALF);
    (
        equal_bytes_8(low.try_into().expect("8 bytes"), byte),
        equal_bytes_8(high.try_into().expect("8 bytes"), byte),
    )
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
