use super::*;

/// The texts of the search tests (also of `ncbi_str`'s `NStr::Replace` and `PrintableString`,
/// and of the header writer): lengths across one and two blocks of `find_first` and beyond, from
/// an alphabet of few bytes (so that needles occur, overlap and stop short at the end) with the
/// bytes that makeblastdb's titles and the report's titles search for, and texts of every byte.
pub(crate) fn search_texts() -> Vec<Vec<u8>> {
    const ALPHABET: &[u8] = b"ab\\t\x01\t. [";
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut texts = Vec::new();
    for len in (0..=100).chain([127, 128, 129, 300]) {
        for round in 0..12 {
            // Some texts draw from the first two letters only, where every needle of them
            // occurs many times and overlaps itself; some from every byte.
            let text = match round % 4 {
                0 => (0..len).map(|_| ALPHABET[(next() % 2) as usize]).collect(),
                1 => (0..len).map(|_| next() as u8).collect(),
                _ => (0..len)
                    .map(|_| ALPHABET[(next() % ALPHABET.len() as u64) as usize])
                    .collect(),
            };
            texts.push(text);
        }
    }
    texts
}

/// The needles of the search tests: every one of up to three bytes from the texts' alphabet
/// that the callers' needles share a first byte with, and the callers' own.
pub(crate) fn search_needles() -> Vec<Vec<u8>> {
    let letters = b"ab\\t\x01";
    let mut needles: Vec<Vec<u8>> = vec![b"\\t".to_vec(), b". [".to_vec(), b", [".to_vec()];
    for &x in letters {
        needles.push(vec![x]);
        for &y in letters {
            needles.push(vec![x, y]);
            for &z in b"ab" {
                needles.push(vec![x, y, z]);
            }
        }
    }
    needles
}

#[test]
fn find_first_is_the_first_index() {
    let tests: [&dyn Fn(u8) -> bool; 4] = [
        &|b| b == b'\\',
        &|b| !(b' '..=b'~').contains(&b),
        &|b| (b == b'\\') | (b == 0x01) | (b == b'\t'),
        &|_| false,
    ];
    for text in search_texts() {
        for test in tests {
            assert_eq!(
                find_first(&text, test),
                text.iter().position(|&b| test(b)),
                "{text:?}"
            );
        }
        for byte in b"ab\\t\x01\t. [z\x00\xff".iter().copied() {
            assert_eq!(
                find_byte(&text, byte),
                text.iter().position(|&b| b == byte),
                "{text:?} {byte}"
            );
        }
    }
    // A byte at each place of texts of up to two blocks and a half, alone.
    for len in [1, 31, 32, 33, 63, 64, 65, 80] {
        for at in 0..len {
            let mut text = vec![b'a'; len];
            text[at] = b'\\';
            assert_eq!(find_byte(&text, b'\\'), Some(at), "{len} {at}");
        }
        assert_eq!(find_byte(&vec![b'a'; len], b'\\'), None);
    }
}

#[test]
fn find_first_finds_a_match_at_every_place_of_every_window() {
    // One match alone, and with a second match at the end, at every place of haystacks across the
    // 32-byte folds, the first window of 64, the later windows and the overlapping last one; and
    // no match.
    let is_b = |byte: u8| byte == b'b';
    for len in 1..=260usize {
        let mut text = vec![b'a'; len];
        assert_eq!(find_first(&text, is_b), None, "{len}");
        for at in 0..len {
            text[at] = b'b';
            assert_eq!(find_first(&text, is_b), Some(at), "{len} {at}");
            if at + 1 < len {
                text[len - 1] = b'b';
                assert_eq!(find_first(&text, is_b), Some(at), "{len} {at} and the end");
                text[len - 1] = b'a';
            }
            text[at] = b'a';
        }
    }
}

#[test]
fn find_is_the_first_occurrence() {
    let texts = search_texts();
    for needle in search_needles() {
        for text in &texts {
            let expected = text
                .windows(needle.len())
                .position(|window| window == needle.as_slice());
            assert_eq!(find(text, &needle), expected, "{text:?} {needle:?}");
        }
    }
    assert_eq!(find(b"abc", b""), Some(0));
    assert_eq!(find(b"", b""), Some(0));
    assert_eq!(find(b"", b"a"), None);
    assert_eq!(find(b"ab", b"abc"), None);
    assert_eq!(find(b"xab", b"ab"), Some(1));
}

#[test]
fn the_mask_of_a_byte_is_the_bytes_equal_to_it() {
    // Every byte value at every place of the block, against every other value around it.
    for byte in 0..=u8::MAX {
        for other in 0..=u8::MAX {
            for place in 0..16 {
                let mut block = [other; 16];
                block[place] = byte;
                for wanted in [byte, other] {
                    let (low, high) = equal_bytes_16(&block, wanted);
                    let mask = gather_high_bits(low) | gather_high_bits(high) << 8;
                    let expected = (0..16)
                        .filter(|&k| block[k] == wanted)
                        .fold(0u32, |mask, k| mask | 1 << k);
                    assert_eq!(
                        mask, expected,
                        "{byte:#04x} {other:#04x} {place} {wanted:#04x}"
                    );
                    assert_eq!(low & !0x8080_8080_8080_8080, 0);
                    assert_eq!(high & !0x8080_8080_8080_8080, 0);
                }
            }
        }
    }
}
