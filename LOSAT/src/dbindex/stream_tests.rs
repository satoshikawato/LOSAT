use super::*;

/// The bases of `letters` one letter at a time: `base_value` of each, its runs of 0, and the
/// values minus 1 packed four to a byte.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:103-112
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:814-876
fn letter_by_letter(letters: &[u8]) -> (Vec<u8>, Vec<(u32, u32)>) {
    let mut packed = vec![0u8; letters.len().div_ceil(CR)];
    let mut invalid: Vec<(u32, u32)> = Vec::new();
    for (pos, &r) in letters.iter().enumerate() {
        let value = match r {
            b'A' => 1,
            b'C' => 2,
            b'G' => 3,
            b'T' => 4,
            _ => 0,
        };
        if value == 0 {
            let pos = pos as u32;
            match invalid.last_mut() {
                Some(run) if run.1 == pos => run.1 = pos + 1,
                _ => invalid.push((pos, pos + 1)),
            }
        }
        let shift = 2 * (CR - 1 - pos % CR);
        packed[pos / CR] |= (value.max(1) - 1) << shift;
    }
    (packed, invalid)
}

/// Letters of every kind, in runs of both lengths around a byte, from a fixed generator.
fn letters(len: usize) -> Vec<u8> {
    const KINDS: &[u8] = b"ACGTACGTACGTNRYacgtn-*\x00\xc1";
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut out = Vec::with_capacity(len);
    while out.len() < len {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let r = KINDS[(state % KINDS.len() as u64) as usize];
        let run = 1 + ((state >> 32) % 9) as usize;
        out.extend(std::iter::repeat_n(r, run.min(len - out.len())));
    }
    out
}

#[test]
fn bases_match_letter_by_letter_for_every_block_split() {
    for len in [0, 1, 3, 4, 5, 8, 9, 1000, 4099] {
        let letters = letters(len);
        let (packed, invalid) = letter_by_letter(&letters);
        for block in [1, 2, 3, 4, 5, 7, 8, 13, 4096] {
            let mut bases = Bases::default();
            for part in letters.chunks(block) {
                bases.push(part);
            }
            bases.finish();
            assert_eq!(bases.len(), len, "len {len} block {block}");
            assert_eq!(bases.packed(), packed.as_slice(), "len {len} block {block}");
            assert_eq!(
                bases.invalid(),
                invalid.as_slice(),
                "len {len} block {block}"
            );
        }
        // Blocks of mixed sizes, empty ones among them.
        let mut bases = Bases::default();
        let mut rest = letters.as_slice();
        for size in [1, 0, 6, 2, 3, 0, 9].iter().cycle() {
            if rest.is_empty() {
                break;
            }
            let (part, after) = rest.split_at((*size).min(rest.len()));
            bases.push(part);
            rest = after;
        }
        bases.push(&[]);
        bases.finish();
        assert_eq!(bases.packed(), packed.as_slice(), "len {len} mixed blocks");
        assert_eq!(
            bases.invalid(),
            invalid.as_slice(),
            "len {len} mixed blocks"
        );
        let whole = Bases::of_letters(&letters);
        assert_eq!(whole.packed(), packed.as_slice(), "len {len}");
        assert_eq!(whole.invalid(), invalid.as_slice(), "len {len}");
    }
}
