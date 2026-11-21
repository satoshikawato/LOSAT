#[cfg(feature = "parallel")]
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

/// A random sequence with low-complexity runs, `N`s, lower case and other IUPAC letters.
#[cfg(feature = "parallel")]
fn sequence(rng: &mut Rng) -> Vec<u8> {
    let len = match rng.below(6) {
        0 => rng.below(4) as usize,
        1 => 3 + rng.below(40) as usize,
        _ => 100 + rng.below(3000) as usize,
    };
    let mut seq = Vec::with_capacity(len);
    while seq.len() < len {
        match rng.below(8) {
            0 => {
                let unit =
                    [b"A".as_slice(), b"AT", b"CAG", b"GGGA", b"TTTTC"][rng.below(5) as usize];
                for _ in 0..(5 + rng.below(60)) {
                    seq.extend_from_slice(unit);
                }
            }
            1 => seq.extend(std::iter::repeat_n(b'N', 1 + rng.below(12) as usize)),
            2 => seq.push(b"RYKMSWnacgt"[rng.below(11) as usize]),
            _ => seq.push(b"ACGT"[rng.below(4) as usize]),
        }
    }
    seq.truncate(len);
    seq
}

/// Parallel == serial on random sequences with low-complexity runs at 1, 2, 3 and 7 threads,
/// with empty sequences and batches of one and of two. (That the masks are made on the pool
/// is the business of `masks_are_made_on_the_pool_threads`; at 1 thread and for batches of
/// one this test compares the serial map with itself.)
#[cfg(feature = "parallel")]
#[test]
fn pool_masks_equal_serial_masks() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut seqs: Vec<Vec<u8>> = (0..90).map(|_| sequence(&mut rng)).collect();
    seqs.push(Vec::new());
    seqs.push(b"AC".to_vec());
    seqs.push(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec());
    let masked: usize = seqs
        .iter()
        .map(|s| DustMasker::new(20, 64, 1).mask_sequence(s).len())
        .sum();
    assert!(
        masked > 20,
        "the sequences must have masked regions ({masked})"
    );
    let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
    for (level, window, linker) in [(20, 64, 1), (10, 32, 2), (30, 16, 1)] {
        let expect = serial(&refs[..], level, window, linker);
        for threads in [1usize, 2, 3, 7] {
            with_search_pool(threads, "blastn", |pool| {
                for take in [refs.len(), 2, 1, 0] {
                    let part = &refs[..take];
                    let want = serial(part, level, window, linker);
                    let got = masks_with_mode(part, level, window, linker, pool, true);
                    assert_eq!(got, want, "threads {threads} take {take}");
                }
                let got = masks_with_mode(&refs[..], level, window, linker, pool, true);
                assert_eq!(got, expect, "threads {threads}");
                Ok(())
            })
            .unwrap();
        }
    }
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

/// With the switch on, the masks are made by the threads of the search pool, and by more than
/// one of them: the probe that the parallel map calls once per query sees a rayon worker of
/// a pool of three threads, and the queries are spread over at least two workers. (Each probe
/// call waits until a second worker has called it (the first call that waits 10 s in vain
/// ends all waiting), so that the spread does not depend on the scheduler; a serial map has
/// no second worker, and fails after those 10 s.) With the
/// switch off, or on a pool of one thread, or for one query, the probe is never called.
#[cfg(feature = "parallel")]
#[test]
fn masks_are_made_on_the_pool_threads() {
    use std::collections::BTreeSet;
    use std::sync::atomic::{AtomicBool, AtomicUsize};
    use std::sync::Mutex;
    use std::time::{Duration, Instant};

    let seqs: Vec<Vec<u8>> = (0..30).map(|i| vec![b"ACGT"[i % 4]; 200]).collect();
    let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
    let want = serial(&refs[..], 20, 64, 1);
    with_search_pool(3, "blastn", |pool| {
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
        calls.store(0, Ordering::SeqCst);
        in_pool_of_three.store(0, Ordering::SeqCst);
        workers.lock().unwrap().clear();
        gave_up.store(false, Ordering::SeqCst);
        let got = masks_with_probe(&refs[..], 20, 64, 1, pool, true, &probe);
        assert_eq!(got, want);
        assert_eq!(calls.load(Ordering::SeqCst), refs.len());
        assert_eq!(
            in_pool_of_three.load(Ordering::SeqCst),
            refs.len(),
            "every query is masked by a worker of the 3-thread pool"
        );
        {
            let workers = workers.lock().unwrap();
            assert!(
                workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                "workers {workers:?}"
            );
        }
        // Not acted on: switch off, one query.
        calls.store(0, Ordering::SeqCst);
        assert_eq!(
            masks_with_probe(&refs[..], 20, 64, 1, pool, false, &probe),
            want
        );
        assert_eq!(
            masks_with_probe(&refs[..1], 20, 64, 1, pool, true, &probe),
            want[..1]
        );
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        Ok(())
    })
    .unwrap();
    with_search_pool(1, "blastn", |pool| {
        let calls = AtomicUsize::new(0);
        let probe = || {
            calls.fetch_add(1, Ordering::SeqCst);
        };
        assert_eq!(
            masks_with_probe(&refs[..], 20, 64, 1, pool, true, &probe),
            want
        );
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        Ok(())
    })
    .unwrap();
}
