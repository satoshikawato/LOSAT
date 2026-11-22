use super::{Arena, Lined, Links, Span, BLOCK, EMPTY, NONE, RUN};

/// A `Lined` buffer begins on a cache line and holds whole lines, for any length, keeps its
/// lines when it shrinks and grows within them, and holds the values set.
#[test]
fn lined_buffers_are_lines_of_their_own() {
    fn check<T: super::Plain + PartialEq + std::fmt::Debug>(v: T, lens: &[usize]) {
        let mut buf: Lined<T> = Lined::default();
        assert!(buf.is_empty());
        for &len in lens {
            let before = buf.lines.as_ptr();
            let had = buf.lines.len();
            buf.set_len(len);
            buf.fill(v);
            assert_eq!(buf.len(), len);
            assert_eq!(buf.as_ptr() as usize % 64, 0, "{len}");
            assert!(
                buf.lines.len() * 64 >= len * std::mem::size_of::<T>(),
                "{len}"
            );
            if had * 64 >= len * std::mem::size_of::<T>() {
                assert_eq!(buf.lines.as_ptr(), before, "{len}: lines kept");
            }
            assert!(buf.iter().all(|x| *x == v), "{len}");
        }
    }
    check(7u16, &[1, 31, 32, 33, 272, 16, 300]);
    check(EMPTY, &[1, 10, 11, 16, 64, 3]);
    check(Span { next: 3, end: 9 }, &[1, 16, 17, 65, 2]);
}

/// `Links::take` hands out distinct blocks of the arena for any cut of a piece's seeds by
/// partition (the pushes of `PieceSink::push`: a block when the partition's last one is full),
/// up to `piece_seeds` seeds: partition p's first block is block p, its later ones whole runs
/// of `RUN` blocks that begin on a run, and the runs taken are the first ones after the first
/// blocks (a piece touches the front of its arena); once the runs are given out, single blocks
/// of other partitions' runs.
#[test]
fn take_gives_each_block_once_in_runs_from_the_front() {
    let mut rng = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = |n: usize| {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        (rng % n as u64) as usize
    };
    let mut tails = 0;
    for parts in [1usize, 2, 3, 5, 16, 64, 128] {
        for piece_seeds in [1usize, 2, 63, 64, 65, 300, 4000, 4096, 16_384] {
            let blocks = Arena::blocks(parts, piece_seeds);
            assert_eq!(blocks % RUN, 0, "{parts} {piece_seeds}");
            let front = parts.next_multiple_of(RUN);
            for cut in 0..6 {
                let mut links = Links::default();
                links.reset(parts, blocks);
                let mut fill = vec![BLOCK; parts];
                let mut taken = vec![false; blocks];
                let mut mine: Vec<Vec<u16>> = vec![Vec::new(); parts];
                for k in 0..piece_seeds {
                    // all to one partition, round robin, uniform, mostly one, each
                    // partition's seeds one block and one more seed, or that and then a
                    // run and one more seed in turn (the most blocks: every run begun)
                    let p = match cut {
                        0 => 0,
                        1 => k % parts,
                        2 => next(parts),
                        3 if next(4) == 0 => next(parts),
                        3 => parts - 1,
                        4 => (k / (BLOCK + 1)) % parts,
                        _ if k < (BLOCK + 1) * parts => k / (BLOCK + 1),
                        _ => ((k - (BLOCK + 1) * parts) / (BLOCK * RUN)) % parts,
                    };
                    if fill[p] == BLOCK {
                        let b = links.take(p);
                        assert!(
                            (b as usize) < blocks && b != NONE,
                            "{parts} {piece_seeds} {cut}"
                        );
                        assert!(
                            !taken[b as usize],
                            "{parts} {piece_seeds} {cut}: block {b} twice"
                        );
                        taken[b as usize] = true;
                        mine[p].push(b);
                        fill[p] = 0;
                    }
                    fill[p] += 1;
                }
                assert_eq!(
                    usize::from(links.used),
                    taken.iter().filter(|&&t| t).count()
                );
                // (blocks of other runs only once the runs are given out)
                if usize::from(links.front) == blocks {
                    tails += 1;
                    continue;
                }
                let mut runs = 0;
                for (p, bs) in mine.iter().enumerate() {
                    if let Some(&first) = bs.first() {
                        assert_eq!(usize::from(first), p, "{parts} {piece_seeds} {cut}");
                    }
                    for (i, &b) in bs.iter().enumerate().skip(1) {
                        let b = usize::from(b);
                        if (i - 1) % RUN == 0 {
                            assert!(b >= front && b % RUN == 0, "{parts} {piece_seeds} {cut}");
                            runs += 1;
                        } else {
                            assert_eq!(
                                b,
                                usize::from(bs[i - 1]) + 1,
                                "{parts} {piece_seeds} {cut}"
                            );
                        }
                    }
                }
                // the runs taken are the first ones after the first blocks
                assert!(
                    taken[front + runs * RUN..].iter().all(|&t| !t),
                    "{parts} {piece_seeds} {cut}"
                );
                assert_eq!(
                    usize::from(links.front),
                    front + runs * RUN,
                    "{parts} {piece_seeds} {cut}"
                );
            }
        }
    }
    // (the cuts reach the arena's end)
    assert!(tails > 0);
}
