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
