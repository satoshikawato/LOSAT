use super::*;

#[test]
fn has_mods_is_the_first_bracket_pair_with_an_equal_sign() {
    assert!(has_mods(b"a [organism=x] b"));
    assert!(has_mods(b"[a=b]"));
    assert!(has_mods(b"[[a=b]]"));
    assert!(!has_mods(b"a [organism] b [x=y]"));
    assert!(!has_mods(b"a ] [x=y]"));
    assert!(!has_mods(b"a [x=y"));
    assert!(!has_mods(b""));
}

/// `has_mods`, which starts `find_brackets` at the first bracket it searches for, is
/// `find_brackets` from the start of the title, on random titles of brackets, `=` and other bytes
/// across the blocks of the search.
#[test]
fn has_mods_from_the_first_bracket_is_find_brackets_from_the_start() {
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    const BYTES: &[u8] = b"[]=ab ";
    for len in 0..300 {
        for round in 0..8 {
            let title: Vec<u8> = (0..len)
                .map(|_| match (round, next() % 32) {
                    (0, _) => BYTES[(next() % BYTES.len() as u64) as usize],
                    (_, 0) => BYTES[(next() % 3) as usize],
                    (_, 1) => next() as u8,
                    _ => b'a',
                })
                .collect();
            assert_eq!(
                has_mods(&title),
                find_brackets(&title, &mut 0).is_some(),
                "{title:?}"
            );
        }
    }
}
