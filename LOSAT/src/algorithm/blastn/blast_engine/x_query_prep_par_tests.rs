use super::*;
#[cfg(feature = "parallel")]
use crate::utils::threading::with_search_pool;

#[cfg(feature = "parallel")]
struct Rng(u64);
#[cfg(feature = "parallel")]
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// A random query with N, IUPAC letters, lower case, low-complexity runs, and now and then
/// no letters at all.
#[cfg(feature = "parallel")]
fn sequence(rng: &mut Rng) -> Vec<u8> {
    let len = match rng.below(6) {
        0 => rng.below(3) as usize,
        1 => 3 + rng.below(40) as usize,
        _ => 100 + rng.below(1500) as usize,
    };
    let mut seq = Vec::with_capacity(len);
    while seq.len() < len {
        match rng.below(8) {
            0 => {
                let unit = [b"A".as_slice(), b"AT", b"CAG", b"GGGA"][rng.below(4) as usize];
                for _ in 0..(3 + rng.below(30)) {
                    seq.extend_from_slice(unit);
                }
            }
            1 => seq.extend(std::iter::repeat_n(b'N', 1 + rng.below(12) as usize)),
            2 => seq.push(b"RYKMSWBDHVnacgtrykmswbdhv-"[rng.below(26) as usize]),
            _ => seq.push(b"ACGT"[rng.below(4) as usize]),
        }
    }
    seq.truncate(len);
    seq
}

/// The context offsets and the length of the block that `search_query_batch` gives to
/// `queries`: `context_room` for each context (its letters and one sentinel, nothing for an
/// empty query), plus context then minus context of each query, and the
/// `query_concat_length` formula of run.rs (the last context's offset and length plus 2, or
/// plus 1 when it is empty, less the two end sentinels).
#[cfg(feature = "parallel")]
fn production_layout(queries: &[Vec<u8>]) -> (Vec<i32>, usize) {
    let room = |length: usize| if length == 0 { 0 } else { length + 1 };
    let mut offsets = Vec::new();
    let mut at = 0usize;
    let mut last = (0usize, 0usize);
    for query in queries {
        let plus = at;
        let minus = plus + room(query.len());
        offsets.push(plus as i32);
        offsets.push(minus as i32);
        at = minus + room(query.len());
        last = (minus, query.len());
    }
    let total = if queries.is_empty() {
        0
    } else {
        (last.0 + last.1 + if last.1 == 0 { 1 } else { 2 }).saturating_sub(2)
    };
    (offsets, total)
}

/// A batch of `count` random queries in the production layout (`production_layout`): the
/// letters, the context offsets and the length of the block.
#[cfg(feature = "parallel")]
fn batch(rng: &mut Rng, count: usize) -> (Vec<Vec<u8>>, Vec<i32>, usize) {
    let queries: Vec<Vec<u8>> = (0..count).map(|_| sequence(rng)).collect();
    let (offsets, total) = production_layout(&queries);
    (queries, offsets, total)
}

/// A batch in a COMPACT layout, which is NOT the production one: the minus context follows
/// the plus context with no sentinel between them, and the next query follows at once (`gap`
/// is 0) or after `gap` sentinels. It tests spans that touch each other, which the production
/// layout (a sentinel after every context) never has.
#[cfg(feature = "parallel")]
fn compact_batch(rng: &mut Rng, count: usize, gap: usize) -> (Vec<Vec<u8>>, Vec<i32>, usize) {
    let queries: Vec<Vec<u8>> = (0..count).map(|_| sequence(rng)).collect();
    let mut offsets = Vec::new();
    let mut at = 0usize;
    for query in &queries {
        offsets.push(at as i32);
        offsets.push((at + query.len()) as i32);
        at += 2 * query.len() + gap;
    }
    (queries, offsets, at)
}

/// The contexts (slices of the block without its end sentinels) of a batch, in context order
/// (an empty query gives two empty slices).
#[cfg(feature = "parallel")]
fn contexts_of<'a>(block: &'a [u8], queries: &[Vec<u8>], offsets: &[i32]) -> Vec<&'a [u8]> {
    queries
        .iter()
        .zip(offsets.chunks_exact(2))
        .flat_map(|(query, pair)| {
            pair.iter().map(move |&offset| {
                if query.is_empty() {
                    &[][..]
                } else {
                    &block[offset as usize..offset as usize + query.len()]
                }
            })
        })
        .collect()
}

#[cfg(feature = "parallel")]
fn spec() -> NuclScoringSpec {
    NuclScoringSpec {
        reward: 1,
        penalty: -2,
        gap_open: 0,
        gap_extend: 0,
    }
}

/// `same_bits` compares every field by its bits: a block that differs in one field, in the sign
/// of a zero, or between `None` and `Some` is different, and a NaN equals itself.
#[test]
fn same_bits_sees_every_field() {
    let base = KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 1.5,
        beta: -2.0,
    };
    assert!(same_bits(&None, &None));
    assert!(!same_bits(&None, &Some(base)));
    assert!(same_bits(&Some(base), &Some(base)));
    let variants = [
        KarlinParams {
            lambda: 0.626,
            ..base
        },
        KarlinParams { k: 0.42, ..base },
        KarlinParams { h: 0.79, ..base },
        KarlinParams { alpha: 1.6, ..base },
        KarlinParams { beta: -2.1, ..base },
    ];
    for variant in variants {
        assert!(!same_bits(&Some(base), &Some(variant)), "{variant:?}");
    }
    let zero = KarlinParams { beta: 0.0, ..base };
    let negative_zero = KarlinParams { beta: -0.0, ..base };
    assert!(!same_bits(&Some(zero), &Some(negative_zero)));
    let nan = KarlinParams {
        k: f64::NAN,
        ..base
    };
    assert!(same_bits(&Some(nan), &Some(nan)));
}

/// Pool == serial for the buffer at 1, 2, 3 and 7 threads on random batches (N, IUPAC, lower
/// case, empty queries, one query, none) in the production layout, which never takes the
/// fallback, and in the compact layouts (spans that touch, or with sentinel gaps).
#[cfg(feature = "parallel")]
#[test]
fn pool_buffer_equals_serial_buffer() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for threads in [1usize, 2, 3, 7] {
        with_search_pool(threads, "blastn", |pool| {
            for (count, gap) in [
                (0, None),
                (1, None),
                (2, None),
                (40, None),
                (90, None),
                (0, Some(0)),
                (1, Some(0)),
                (2, Some(0)),
                (40, Some(0)),
                (40, Some(3)),
                (90, Some(0)),
            ] {
                let (queries, offsets, total) = match gap {
                    None => batch(&mut rng, count),
                    Some(gap) => compact_batch(&mut rng, count, gap),
                };
                let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
                let want = serial_buffer(&letters[..], &offsets, total);
                let got = buffer_with_mode(&letters[..], &offsets, total, pool, true);
                assert_eq!(got, want, "threads {threads} count {count} gap {gap:?}");
                // the batch loop's kept buffer, written again over other bytes
                let mut kept = got.clone();
                kept.fill(0xA5);
                fill_with_mode(&letters[..], &offsets, total, pool, true, &mut kept);
                assert_eq!(
                    kept, want,
                    "kept: threads {threads} count {count} gap {gap:?}"
                );
                if gap.is_none() && count > 1 && pool.threads() > 1 {
                    let (got, on_the_pool) =
                        parallel_buffer(&letters[..], &offsets, total, pool, &|| {});
                    assert!(on_the_pool, "the production layout took the fallback");
                    assert_eq!(got, want);
                }
            }
            Ok(())
        })
        .unwrap();
    }
}

/// Contexts that do not end inside the block, empty contexts and overlapping spans give the
/// serial buffer (the serial function stays the reference for them).
#[cfg(feature = "parallel")]
#[test]
fn irregular_layouts_equal_serial_buffer() {
    let mut rng = Rng(0x1234_5678_9abc_def1);
    with_search_pool(3, "blastn", |pool| {
        let (queries, offsets, total) = batch(&mut rng, 30);
        let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        for mutate in 0..4 {
            let mut offsets = offsets.clone();
            match mutate {
                // A context past the block (stays sentinels).
                0 => offsets[7] = total as i32 + 5,
                // A negative offset (clamped to 0, as the serial function does).
                1 => offsets[10] = -4,
                // Queries out of order and overlapping.
                2 => offsets.swap(6, 20),
                // A minus context of one query overlapping the plus context of the next.
                _ => offsets[9] = offsets[10] - 1,
            }
            let want = serial_buffer(&letters[..], &offsets, total);
            let got = buffer_with_mode(&letters[..], &offsets, total, pool, true);
            assert_eq!(got, want, "mutation {mutate}");
            // Out of order or overlapping spans are the fallback (the serial function
            // made the buffer, nothing was encoded on the pool).
            let (got, on_the_pool) = parallel_buffer(&letters[..], &offsets, total, pool, &|| {});
            assert_eq!(got, want, "mutation {mutate}");
            assert_eq!(on_the_pool, mutate < 2, "mutation {mutate}");
        }
        Ok(())
    })
    .unwrap();
}

/// Pool == serial for the Karlin blocks, bit for bit, at 1, 2, 3 and 7 threads on the
/// contexts of random batches (BLASTNA letters from the buffer, strands in both orders).
#[cfg(feature = "parallel")]
#[test]
fn pool_blocks_equal_serial_blocks() {
    let mut rng = Rng(0xdead_beef_cafe_f00d);
    let (queries, offsets, total) = batch(&mut rng, 70);
    let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let buffer = serial_buffer(&letters[..], &offsets, total);
    let contexts = contexts_of(&buffer[1..total + 1], &queries, &offsets);
    let mut compared = 0;
    for (reward, penalty) in [(1, -2), (2, -3), (1, -3), (4, -5)] {
        let spec = NuclScoringSpec {
            reward,
            penalty,
            ..spec()
        };
        let want = serial_blocks(&contexts, &spec);
        compared += want.iter().filter(|block| block.is_some()).count();
        for threads in [1usize, 2, 3, 7] {
            with_search_pool(threads, "blastn", |pool| {
                for take in [contexts.len(), 3, 2, 1, 0] {
                    let part = &contexts[..take];
                    let got = blocks_with_mode(part, &spec, pool, true);
                    assert_eq!(got.len(), take);
                    for (i, (got, want)) in got.iter().zip(&want).enumerate() {
                        assert!(
                            same_bits(got, want),
                            "threads {threads} take {take} context {i}: {got:?} / {want:?}"
                        );
                    }
                }
                Ok(())
            })
            .unwrap();
        }
    }
    assert!(compared > 200, "blocks compared: {compared}");
}

/// The parallel path acts only where a search pool of two or more threads is current.
#[cfg(feature = "parallel")]
#[test]
fn only_a_search_pool_of_two_threads_acts() {
    with_search_pool(1, "blastn", |pool| {
        assert!(!on_pool(pool, 10));
        Ok(())
    })
    .unwrap();
    with_search_pool(3, "blastn", |pool| {
        assert!(on_pool(pool, 10));
        assert!(!on_pool(pool, 1));
        assert!(!on_pool(pool, 0));
        Ok(())
    })
    .unwrap();
}

/// With the switch on, both steps run on the threads of the search pool, and on more than one
/// of them: the probe that each parallel map calls once per item sees a rayon worker of a pool
/// of three threads, and the items are spread over at least two workers. (Each probe call waits
/// until a second worker has called it, the first call that waits 10 s in vain ends all
/// waiting, so that the spread does not depend on the scheduler; a serial map has no second
/// worker and fails after those 10 s.) With the switch off, on a pool of one thread, or for one
/// item, the probe is never called. Only where the pool has threads, as `place`: wasm32 without
/// `wasm-threads` has no second worker to find.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn steps_run_on_the_pool_threads() {
    use std::collections::BTreeSet;
    use std::sync::atomic::{AtomicBool, AtomicUsize};
    use std::sync::Mutex;
    use std::time::{Duration, Instant};

    let mut rng = Rng(0x0bad_cafe_1234_5678);
    let (queries, offsets, total) = batch(&mut rng, 30);
    let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let want_buffer = serial_buffer(&letters[..], &offsets, total);
    let contexts = contexts_of(&want_buffer[1..total + 1], &queries, &offsets);
    let spec = spec();
    let want_blocks = serial_blocks(&contexts, &spec);

    let calls = AtomicUsize::new(0);
    let workers = Mutex::new(BTreeSet::new());
    let in_pool_of_three = AtomicUsize::new(0);
    let gave_up = AtomicBool::new(false);
    let probe = || {
        calls.fetch_add(1, Ordering::SeqCst);
        if let Some(index) = rayon::current_thread_index() {
            workers.lock().unwrap().insert(index);
            if rayon::current_num_threads() == 3 {
                in_pool_of_three.fetch_add(1, Ordering::SeqCst);
            }
        }
        let start = Instant::now();
        while !gave_up.load(Ordering::SeqCst) && workers.lock().unwrap().len() < 2 {
            if start.elapsed() > Duration::from_secs(10) {
                gave_up.store(true, Ordering::SeqCst);
            }
            std::thread::yield_now();
        }
    };
    let reset = || {
        calls.store(0, Ordering::SeqCst);
        in_pool_of_three.store(0, Ordering::SeqCst);
        workers.lock().unwrap().clear();
        gave_up.store(false, Ordering::SeqCst);
    };
    // The queries that have a placed context: the probe is called once for each of them.
    let placed = letters
        .iter()
        .zip(offsets.chunks_exact(2))
        .filter(|(letters, pair)| {
            pair.iter()
                .any(|&offset| place(offset, letters.len(), total).is_some())
        })
        .count();
    with_search_pool(3, "blastn", |pool| {
        reset();
        let got = buffer_with_probe(&letters[..], &offsets, total, pool, true, &probe);
        assert_eq!(got, want_buffer);
        assert_eq!(calls.load(Ordering::SeqCst), placed);
        assert_eq!(in_pool_of_three.load(Ordering::SeqCst), placed);
        {
            let workers = workers.lock().unwrap();
            assert!(
                workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                "workers {workers:?}"
            );
        }
        reset();
        let got = blocks_with_probe(&contexts, &spec, pool, true, &probe);
        assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
        assert_eq!(calls.load(Ordering::SeqCst), contexts.len());
        assert_eq!(in_pool_of_three.load(Ordering::SeqCst), contexts.len());
        {
            let workers = workers.lock().unwrap();
            assert!(
                workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                "workers {workers:?}"
            );
        }
        // Not acted on: switch off, one query, one context.
        reset();
        assert_eq!(
            buffer_with_probe(&letters[..], &offsets, total, pool, false, &probe),
            want_buffer
        );
        let got = blocks_with_probe(&contexts, &spec, pool, false, &probe);
        assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
        let one = buffer_with_probe(&letters[..1], &offsets[..2], total, pool, true, &probe);
        assert_eq!(one, serial_buffer(&letters[..1], &offsets[..2], total));
        let got = blocks_with_probe(&contexts[..1], &spec, pool, true, &probe);
        assert!(same_bits(&got[0], &want_blocks[0]));
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        Ok(())
    })
    .unwrap();
    with_search_pool(1, "blastn", |pool| {
        reset();
        assert_eq!(
            buffer_with_probe(&letters[..], &offsets, total, pool, true, &probe),
            want_buffer
        );
        let got = blocks_with_probe(&contexts, &spec, pool, true, &probe);
        assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        Ok(())
    })
    .unwrap();
}

/// A random query of 0 to `max` letters: often empty, one letter or very short, with bytes of
/// every value, IUPAC letters and lower case (the reviewer's generator, for the production
/// layout tests).
#[cfg(feature = "parallel")]
fn varied_sequence(rng: &mut Rng, max: usize) -> Vec<u8> {
    let len = match rng.below(7) {
        0 => 0,
        1 => 1,
        2 => rng.below(4) as usize,
        _ => rng.below(max as u64 + 1) as usize,
    };
    (0..len)
        .map(|_| match rng.below(10) {
            0 => rng.below(256) as u8,
            1 => b"NnRYKMSWBDHVrykmswbdhvUu-*"[rng.below(26) as usize],
            2 => b"acgt"[rng.below(4) as usize],
            _ => b"ACGT"[rng.below(4) as usize],
        })
        .collect()
}

/// The production layout (offsets of `search_query_batch`, empty queries first or last,
/// bytes 0..=255) on pools of 2, 3, 5 and 16 threads, with two score pairs:
/// the buffer byte for byte, the blocks bit for bit, and never the fallback.
#[cfg(feature = "parallel")]
#[test]
fn production_layout_buffers_and_blocks() {
    let mut rng = Rng(0x5151_7e57_0dd5_eed5);
    let mut cases = 0;
    for threads in [2usize, 3, 5, 16] {
        with_search_pool(threads, "blastn", |pool| {
            for (count, max) in [
                (2usize, 50usize),
                (3, 5),
                (7, 3000),
                (64, 400),
                (300, 60),
                (2_000, 30),
            ] {
                for _ in 0..2 {
                    let mut queries: Vec<Vec<u8>> =
                        (0..count).map(|_| varied_sequence(&mut rng, max)).collect();
                    if rng.below(3) == 0 {
                        queries.push(Vec::new());
                    }
                    if rng.below(3) == 0 {
                        queries.insert(0, Vec::new());
                    }
                    if queries.iter().all(|query| query.is_empty()) {
                        queries.push(b"ACGT".to_vec());
                    }
                    let (offsets, total) = production_layout(&queries);
                    let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
                    let want = serial_buffer(&letters[..], &offsets, total);
                    let (got, on_the_pool) =
                        parallel_buffer(&letters[..], &offsets, total, pool, &|| {});
                    assert!(on_the_pool, "the production layout took the fallback");
                    assert_eq!(got, want, "threads {threads} count {count}");
                    let got = buffer_with_mode(&letters[..], &offsets, total, pool, true);
                    assert_eq!(got, want, "threads {threads} count {count}");
                    let contexts = contexts_of(&want[1..total + 1], &queries, &offsets);
                    for (reward, penalty) in [(1, -2), (2, -3)] {
                        let spec = NuclScoringSpec {
                            reward,
                            penalty,
                            gap_open: 0,
                            gap_extend: 0,
                        };
                        let want = serial_blocks(&contexts, &spec);
                        let got = blocks_with_mode(&contexts, &spec, pool, true);
                        assert_eq!(got.len(), want.len());
                        assert!(
                            got.iter().zip(&want).all(|(a, b)| same_bits(a, b)),
                            "threads {threads} count {count} r {reward} p {penalty}"
                        );
                    }
                    cases += 1;
                }
            }
            Ok(())
        })
        .unwrap();
    }
    assert!(cases >= 40, "{cases}");
}

/// Called from a plain thread that is not a worker of the pool it is given: the work must still
/// run on that pool (this is what `SearchPool::install` adds), never on rayon's global pool.
/// Without `install` in the two parallel steps the probe sees a thread outside the pool, and
/// this test fails.
#[cfg(feature = "parallel")]
#[test]
fn caller_outside_the_pool_still_runs_on_it() {
    use std::sync::atomic::AtomicUsize;
    let mut rng = Rng(0x7777_1234_abcd_0001);
    let queries: Vec<Vec<u8>> = (0..64).map(|_| varied_sequence(&mut rng, 500)).collect();
    let (offsets, total) = production_layout(&queries);
    let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let want = serial_buffer(&letters[..], &offsets, total);
    let contexts = contexts_of(&want[1..total + 1], &queries, &offsets);
    let spec = spec();
    with_search_pool(3, "blastn", |pool| {
        std::thread::scope(|scope| {
            scope
                .spawn(|| {
                    assert!(rayon::current_thread_index().is_none());
                    let wrong = AtomicUsize::new(0);
                    let probe = || {
                        if rayon::current_num_threads() != 3
                            || rayon::current_thread_index().is_none()
                        {
                            wrong.fetch_add(1, Ordering::SeqCst);
                        }
                    };
                    let got = buffer_with_probe(&letters[..], &offsets, total, pool, true, &probe);
                    assert_eq!(got, want);
                    let got = blocks_with_probe(&contexts, &spec, pool, true, &probe);
                    assert_eq!(got.len(), contexts.len());
                    assert_eq!(
                        wrong.load(Ordering::SeqCst),
                        0,
                        "work ran off the search pool"
                    );
                })
                .join()
                .unwrap();
        });
        Ok(())
    })
    .unwrap();
}

/// A pool other than the search's: a plain thread with a private pool of 2 while the outer pool
/// of 4 exists, and the steps called from inside a parallel loop of the outer pool (nested).
#[cfg(feature = "parallel")]
#[test]
fn private_pool_and_nested_calls() {
    use rayon::prelude::*;
    use std::sync::atomic::AtomicUsize;
    let mut rng = Rng(0x0123_4567_89ab_cdef);
    let queries: Vec<Vec<u8>> = (0..200).map(|_| varied_sequence(&mut rng, 800)).collect();
    let (offsets, total) = production_layout(&queries);
    let letters: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let want = serial_buffer(&letters[..], &offsets, total);
    let contexts = contexts_of(&want[1..total + 1], &queries, &offsets);
    let spec = spec();
    let want_blocks = serial_blocks(&contexts, &spec);
    with_search_pool(4, "blastn", |outer| {
        std::thread::scope(|scope| {
            let private_thread = scope.spawn(|| {
                with_search_pool(2, "blastn", |private| {
                    let wrong = AtomicUsize::new(0);
                    let probe = || {
                        if rayon::current_num_threads() != 2 {
                            wrong.fetch_add(1, Ordering::SeqCst);
                        }
                    };
                    let got =
                        buffer_with_probe(&letters[..], &offsets, total, private, true, &probe);
                    assert_eq!(got, want);
                    let got = blocks_with_probe(&contexts, &spec, private, true, &probe);
                    assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
                    assert_eq!(wrong.load(Ordering::SeqCst), 0);
                    Ok(())
                })
                .unwrap()
            });
            outer.install(|| {
                (0..6).into_par_iter().for_each(|_| {
                    let got = buffer_with_mode(&letters[..], &offsets, total, outer, true);
                    assert_eq!(got, want);
                    let got = blocks_with_mode(&contexts, &spec, outer, true);
                    assert!(got.iter().zip(&want_blocks).all(|(a, b)| same_bits(a, b)));
                })
            });
            private_thread.join().unwrap();
        });
        Ok(())
    })
    .unwrap();
}
