//! The searches of a byte string that LOSAT's ports share: the first byte that a test holds
//! for ([`find_first`]), the first of one byte ([`find_byte`], C's `memchr`) and the first
//! occurrence of a byte string ([`find`], `string::find`). makeblastdb's FASTA reader, its titles
//! and its header writer search every byte of every definition line, so the search tests the
//! bytes a block at a time.
//!
//! No NCBI counterpart: NCBI's code calls the C and C++ libraries for these searches; LOSAT keeps
//! what they return (the index of the first match), not how they search.

/// The bytes that [`find_first`] tests at once.
const BLOCK: usize = 32;

/// The index of the first byte of `haystack` that `test` holds for. The bytes are tested a block
/// of [`BLOCK`] at a time, in a loop without an early exit, so that the compiler can test a block
/// with vector compares; the block that holds one is then searched byte by byte. The last block
/// ends at the end of `haystack` and overlaps the one before it, where `test` holds for no byte;
/// only a haystack shorter than a block is searched byte by byte alone.
pub(crate) fn find_first(haystack: &[u8], test: impl Fn(u8) -> bool) -> Option<usize> {
    let holds = |block: &[u8; BLOCK]| block.iter().fold(false, |any, &b| any | test(b));
    let first = |block: &[u8], start: usize| {
        block
            .iter()
            .position(|&b| test(b))
            .map(|offset| start + offset)
    };
    let (blocks, tail) = haystack.as_chunks::<BLOCK>();
    for (index, block) in blocks.iter().enumerate() {
        if holds(block) {
            return first(block, index * BLOCK);
        }
    }
    if blocks.is_empty() || tail.is_empty() {
        return first(tail, 0);
    }
    let start = haystack.len() - BLOCK;
    let last = &haystack[start..].as_chunks::<BLOCK>().0[0];
    if holds(last) {
        first(last, start)
    } else {
        None
    }
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

#[cfg(test)]
#[path = "byte_search_tests.rs"]
pub(crate) mod tests;
