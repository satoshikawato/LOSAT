use super::*;

#[test]
fn runs_of_n_at_the_ends_always_count_and_inner_ones_only_when_longer_than_max() {
    assert_eq!(
        long_n_runs(2, b"NNACNNTNNNG"),
        vec![(0, 1), (7, 9)],
        "the run of 2 inside is not longer than 2"
    );
    assert_eq!(long_n_runs(64, b"ACGTN"), vec![(4, 4)]);
    assert_eq!(long_n_runs(u32::MAX, b"NNNN"), vec![(0, 3)]);
    assert!(long_n_runs(0, b"ACGT").is_empty());
}

#[test]
fn merges_take_the_linker_as_an_unsigned_sum() {
    let mut list = vec![(0, 9)];
    insert_merge(&mut list, (10, 20), 1);
    assert_eq!(list, vec![(0, 20)]);
    // -1 as Uint4: 20 + 4294967295 wraps to 19.
    insert_merge(&mut list, (19, 30), u32::MAX);
    assert_eq!(list, vec![(0, 30)]);
    insert_merge(&mut list, (40, 50), 1);
    assert_eq!(list, vec![(0, 30), (40, 50)]);
}

/// The letter-by-letter count of NCBI's `s_FindSegmentWithLongNs` (`Uint4` position and
/// length), which `long_n_runs` gives block by block.
fn long_n_runs_by_letter(max_ns: u32, seq: &[u8]) -> Vec<Interval> {
    let mut runs = Vec::new();
    let (mut pos, mut ns) = (0u32, 0u32);
    for &letter in seq {
        if letter == b'N' {
            ns = ns.wrapping_add(1);
        } else {
            if ns > 0 {
                if ns > max_ns || pos == 0 {
                    runs.push((pos, pos.wrapping_add(ns).wrapping_sub(1)));
                }
                pos = pos.wrapping_add(ns);
                ns = 0;
            }
            pos = pos.wrapping_add(1);
        }
    }
    if ns > 0 {
        runs.push((pos, pos.wrapping_add(ns).wrapping_sub(1)));
    }
    runs
}

#[test]
fn runs_found_block_by_block_are_those_counted_letter_by_letter() {
    // Runs and gaps of every length around the 16-letter blocks, at every offset, from a fixed
    // linear congruential sequence.
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = |bound: u64| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        (state >> 33) % bound
    };
    for _ in 0..4000 {
        let mut seq = Vec::new();
        let pieces = next(8);
        for _ in 0..pieces {
            let letter = if next(2) == 0 {
                b'N'
            } else {
                b"ACGTn"[next(5) as usize]
            };
            let length = [0, 1, 2, 3, 5, 7, 15, 16, 17, 31, 32, 33, 40][next(13) as usize];
            seq.extend(std::iter::repeat_n(letter, length));
        }
        for max_ns in [0, 1, 2, 3, 7, 8, 14, 15, 16, 17, 64, u32::MAX] {
            assert_eq!(
                long_n_runs(max_ns, &seq),
                long_n_runs_by_letter(max_ns, &seq),
                "{max_ns} {:?}",
                String::from_utf8_lossy(&seq)
            );
        }
    }
}

#[test]
fn every_block_of_16_letters_gives_the_runs_counted_letter_by_letter() {
    // Every pattern of `N` in the second block, after a block of `N` (a run that comes in) or of
    // `A`, with no letter after it or with one `N` or `A`.
    for pattern in 0..=u16::MAX {
        for (before, after) in [
            (b'N', None),
            (b'A', None),
            (b'N', Some(b'N')),
            (b'A', Some(b'A')),
        ] {
            let mut seq = vec![before; 16];
            seq.extend((0..16).map(|k| if pattern >> k & 1 == 1 { b'N' } else { b'C' }));
            seq.extend(after);
            for max_ns in [0, 1, 2, 7, 14, 15, 64] {
                assert_eq!(
                    long_n_runs(max_ns, &seq),
                    long_n_runs_by_letter(max_ns, &seq),
                    "{max_ns} {pattern:#06x} {before} {after:?}"
                );
            }
        }
    }
}
