use super::*;

/// FNV-1a (64 bits) of the elements, each as 4 bytes, little endian.
fn digest(items: &[u32]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in items.iter().flat_map(|item| item.to_le_bytes()) {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
    hash
}

fn shuffled(seed: u32, n: u32) -> (Vec<u32>, Mt19937) {
    let mut items: Vec<u32> = (0..n).collect();
    let mut engine = Mt19937::new(seed);
    shuffle(&mut items, &mut engine);
    (items, engine)
}

/// The standard's check of `mt19937` ([rand.predef]) and its first outputs after the default seed.
#[test]
fn mt19937_is_the_standard_s() {
    let mut engine = Mt19937::new(5489);
    let first: Vec<u32> = (0..3).map(|_| engine.next_u32()).collect();
    assert_eq!(first, [3_499_211_612, 581_869_302, 3_890_346_734]);
    for _ in 3..9_999 {
        engine.next_u32();
    }
    assert_eq!(engine.next_u32(), 4_123_659_995);
}

/// Permutations that GCC 13's `std::shuffle` gave a `vector<int>` of 0 to n - 1 (black-box
/// records of line RG item D3d, 2026-11-21): the elements of small vectors, and for large ones
/// the digest of the permutation and the engine's next output after it, which pins the count of
/// draws (draws again in Lemire's method included: 38,844 draws for 65,535 elements under the
/// seed 5489). 65,535 elements take two positions per draw; 65,536 and more take one.
#[test]
fn shuffle_gives_gcc_13_s_permutations() {
    assert!(shuffled(1, 0).0.is_empty());
    assert_eq!(shuffled(1, 1).0, [0]);
    assert_eq!(shuffled(5489, 6).0, [2, 3, 0, 5, 4, 1]);
    assert_eq!(
        shuffled(1_763_700_000, 10).0,
        [9, 8, 0, 6, 7, 3, 5, 2, 4, 1]
    );
    assert_eq!(shuffled(0, 7).0, [2, 4, 6, 0, 5, 1, 3]);
    assert_eq!(shuffled(u32::MAX, 5).0, [3, 0, 4, 2, 1]);
    for (seed, n, expected, next) in [
        (1_763_712_345, 2, 0x89cd_3129_1d2a_efa4, 2_426_647_230),
        (1_763_712_345, 3, 0x9d19_bb4b_d820_c026, 2_426_647_230),
        (1_763_712_345, 65_534, 0x63c6_f5aa_67c3_72a4, 767_310_204),
        (1_763_712_345, 65_535, 0x9721_6598_84a4_852b, 995_136_923),
        (1_763_712_345, 65_536, 0xe55c_e52f_84f9_d2b1, 2_339_223_694),
        (1_763_712_345, 65_537, 0x7226_6469_bfc4_44b4, 514_753_640),
        (
            1_763_712_345,
            1_000_000,
            0x796b_fbe3_4a33_9295,
            1_808_644_962,
        ),
        (5489, 65_535, 0xa4ec_20f7_5916_57ff, 904_011_123),
    ] {
        let (items, mut engine) = shuffled(seed, n);
        assert_eq!(digest(&items), expected, "seed {seed}, {n} elements");
        assert_eq!(engine.next_u32(), next, "seed {seed}, {n} elements");
    }
    // 38,844 draws: the next output is the engine's 38,845th.
    let mut fresh = Mt19937::new(5489);
    for _ in 0..38_844 {
        fresh.next_u32();
    }
    assert_eq!(fresh.next_u32(), 904_011_123);
}
