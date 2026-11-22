use super::*;
use std::collections::VecDeque;

type W = (usize, usize, u64);

/// A FIFO of `depth` words like the head's ring: a push enters the word and, once `depth`
/// words are in, hands the oldest over to `emit` (and says so); a flush hands over the rest.
struct Fifo {
    q: VecDeque<W>,
    depth: usize,
    pushes: usize,
}

impl x_dc_csr::XRing<W> for Fifo {
    fn push(&mut self, w: W, _cell: impl Fn(&W) -> usize, emit: &mut impl FnMut(W, usize)) -> bool {
        self.pushes += 1;
        self.q.push_back(w);
        if self.q.len() > self.depth {
            emit(self.q.pop_front().unwrap(), 0);
            true
        } else {
            false
        }
    }

    fn flush(&mut self, _cell: impl Fn(&W) -> usize, emit: &mut impl FnMut(W, usize)) {
        while let Some(w) = self.q.pop_front() {
            emit(w, 0);
        }
    }
}

/// `XWordBatch` hands the ring the words that pass, in the scan's order, so the ring hands
/// over the same words in the same order as when each passing word is pushed as it comes:
/// with a ring that fills at once, after a few words, at the scans' depth or never; across full
/// batches and the last part (`finish`); with none, every and every k-th word passing, and with
/// pseudo-random words and passes.
#[test]
fn the_word_batch_hands_the_ring_the_passing_words_in_order() {
    let mut x = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let lens = [0, 1, 2, 31, 32, 33, 63, 64, 65, 96, 1000];
    for depth in [0, 1, 5, x_dc_csr::RING, 10_000] {
        for len in lens {
            for rule in 0..6u64 {
                let words: Vec<(W, bool)> = (0..len)
                    .map(|i| {
                        let r = next();
                        let pass = match rule {
                            0 => false,
                            1 => true,
                            2 => i % 2 == 0,
                            3 => i % X_WORD_BATCH == X_WORD_BATCH - 1,
                            4 => r % 10 != 0,
                            _ => r % 7 == 0,
                        };
                        ((i, (r >> 40) as usize, r), pass)
                    })
                    .collect();
                let mut ring = Fifo {
                    q: VecDeque::new(),
                    depth,
                    pushes: 0,
                };
                let mut out = Vec::new();
                let mut emit = |w: W, _: usize| out.push(w);
                let mut batch = XWordBatch::new(&mut ring, &mut emit);
                for &(w, pass) in &words {
                    if pass {
                        batch.keep(w.0, w.1, w.2);
                    }
                }
                batch.finish();
                x_dc_csr::XRing::flush(&mut ring, |w: &W| w.2 as usize, &mut emit);
                let want: Vec<W> = words.iter().filter(|w| w.1).map(|w| w.0).collect();
                assert_eq!(
                    ring.pushes,
                    want.len(),
                    "depth {depth} len {len} rule {rule}"
                );
                assert_eq!(out, want, "depth {depth} len {len} rule {rule}");
            }
        }
    }
}
