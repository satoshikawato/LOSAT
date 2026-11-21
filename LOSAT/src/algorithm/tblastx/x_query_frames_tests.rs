use super::*;
use crate::algorithm::tblastx::translation::generate_frames;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::threading::with_search_pool;

struct Rng(u64);
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

/// A random DNA sequence with low-complexity runs (so that SEG masks some of every frame),
/// `N`s and other IUPAC letters; short ones too (no frame, or frames of zero or one letter).
fn sequence(rng: &mut Rng) -> Vec<u8> {
    let len = match rng.below(6) {
        0 => rng.below(8) as usize,
        1 => 8 + rng.below(60) as usize,
        _ => 300 + rng.below(4000) as usize,
    };
    let mut seq = Vec::with_capacity(len);
    while seq.len() < len {
        match rng.below(8) {
            0 => {
                let unit = [b"CAG".as_slice(), b"GCA", b"AAAAAG", b"CCTCCA", b"GGT"]
                    [rng.below(5) as usize];
                for _ in 0..(5 + rng.below(40)) {
                    seq.extend_from_slice(unit);
                }
            }
            1 => seq.extend(std::iter::repeat_n(b'N', 1 + rng.below(9) as usize)),
            2 => seq.push(b"RYKMSWnacgt"[rng.below(11) as usize]),
            _ => seq.push(b"ACGT"[rng.below(4) as usize]),
        }
    }
    seq.truncate(len);
    seq
}

fn frames_of(seqs: &[Vec<u8>]) -> Vec<Vec<QueryFrame>> {
    let code = GeneticCode::from_id(1);
    seqs.iter().map(|s| generate_frames(s, &code)).collect()
}

/// `-lcase_masking`: a lower-case run reaches each frame in protein coordinates. Hard frames
/// get X there and keep their letters in the unmasked copy, soft frames keep their letters.
#[test]
fn lcase_runs_mask_each_frame() {
    let mut letters = b"ATGGCTAAACCCGGGTTTGATGAACATCTGAGCACCTGGTACCGTATC".to_vec();
    letters[9..18].make_ascii_lowercase();
    let record = FastaRecord {
        sequence: letters.clone(),
        ..Default::default()
    };
    let code = GeneticCode::from_id(1);
    let plain = generate_frames(&letters, &code);

    let mut hard = vec![plain.clone()];
    let mut soft = vec![plain.clone()];
    with_search_pool(1, "tblastx", |pool| {
        mask_lcase(&mut hard, std::slice::from_ref(&record), None, false, pool);
        mask_lcase(&mut soft, std::slice::from_ref(&record), None, true, pool);
        Ok(())
    })
    .unwrap();

    assert_eq!(hard[0][0].seg_masks, vec![(3, 6)]);
    assert_eq!(&hard[0][0].aa_seq[4..7], &[21, 21, 21]);
    assert_eq!(hard[0][0].aa_seq_nomask.as_ref(), Some(&plain[0].aa_seq));
    for (frame, (hard, plain)) in soft[0].iter().zip(hard[0].iter().zip(&plain)) {
        assert_eq!(frame.seg_masks, hard.seg_masks);
        assert_eq!(frame.aa_seq, plain.aa_seq);
        assert!(frame.aa_seq_nomask.is_none());
    }
}

/// The fields SEG writes, per frame, for comparison.
fn masked(frames: &[Vec<QueryFrame>]) -> Vec<(Vec<u8>, Option<Vec<u8>>, Vec<(usize, usize)>)> {
    frames
        .iter()
        .flatten()
        .map(|f| {
            (
                f.aa_seq.clone(),
                f.aa_seq_nomask.clone(),
                f.seg_masks.clone(),
            )
        })
        .collect()
}

/// Pool == serial on random queries with low-complexity runs at 1, 2, 3 and 7 threads, for
/// several SEG parameter sets, for batches of many queries, of one query and of none. (That the
/// frames are masked on the pool is the business of `frames_are_masked_on_the_pool_threads`.)
#[cfg(feature = "parallel")]
#[test]
fn pool_frames_equal_serial_frames() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut seqs: Vec<Vec<u8>> = (0..40).map(|_| sequence(&mut rng)).collect();
    seqs.push(Vec::new());
    seqs.push(b"ACG".to_vec());
    seqs.push(b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG".repeat(4));
    for (window, locut, hicut) in [(12, 2.2, 2.5), (10, 1.8, 2.1), (20, 3.0, 3.3)] {
        let seg = SegMasker::new(window, locut, hicut);
        let mut want_all = frames_of(&seqs);
        serial::<false>(&mut want_all, &seg);
        let intervals: usize = want_all.iter().flatten().map(|f| f.seg_masks.len()).sum();
        assert!(
            intervals > 40,
            "the frames must have masked regions ({intervals})"
        );
        for threads in [1usize, 2, 3, 7] {
            with_search_pool(threads, "tblastx", |pool| {
                for take in [seqs.len(), 2, 1, 0] {
                    let mut want = frames_of(&seqs[..take]);
                    serial::<false>(&mut want, &seg);
                    let mut got = frames_of(&seqs[..take]);
                    mask_with_mode::<false>(&mut got, &seg, pool, true);
                    assert_eq!(
                        masked(&got),
                        masked(&want),
                        "window {window} threads {threads} take {take}"
                    );
                }
                let mut got = frames_of(&seqs);
                mask_with_mode::<false>(&mut got, &seg, pool, true);
                assert_eq!(masked(&got), masked(&want_all), "threads {threads}");
                Ok(())
            })
            .unwrap();
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
/// `-soft_masking true`: every frame gets the intervals of the hard frame and keeps its letters
/// (no unmasked copy), in the serial loop and on pools of 1, 2 and 3 threads.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn soft_frames_keep_letters_and_hard_intervals() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let seqs: Vec<Vec<u8>> = (0..12).map(|_| sequence(&mut rng)).collect();
    let seg = SegMasker::new(12, 2.2, 2.5);
    let plain = frames_of(&seqs);
    let mut hard = frames_of(&seqs);
    serial::<false>(&mut hard, &seg);
    let mut soft = frames_of(&seqs);
    serial::<true>(&mut soft, &seg);
    let intervals: usize = hard.iter().flatten().map(|f| f.seg_masks.len()).sum();
    assert!(intervals > 0, "the frames must have masked regions");
    let frames = soft.iter().flatten().zip(hard.iter().flatten());
    for ((soft, hard), plain) in frames.zip(plain.iter().flatten()) {
        assert_eq!(soft.seg_masks, hard.seg_masks);
        assert_eq!(soft.aa_seq, plain.aa_seq);
        assert!(soft.aa_seq_nomask.is_none());
        assert_eq!(
            hard.aa_seq_nomask.as_ref().unwrap_or(&hard.aa_seq),
            &plain.aa_seq
        );
    }
    for threads in [1usize, 2, 3] {
        with_search_pool(threads, "tblastx", |pool| {
            let mut got = frames_of(&seqs);
            mask_with_mode::<true>(&mut got, &seg, pool, true);
            assert_eq!(masked(&got), masked(&soft), "threads {threads}");
            Ok(())
        })
        .unwrap();
    }
}

/// The parallel path acts only where a search pool of two or more threads is current and the
/// batch has two or more frames that SEG reads.
#[cfg(feature = "parallel")]
#[test]
fn only_a_search_pool_of_two_threads_acts() {
    with_search_pool(1, "tblastx", |pool| {
        assert!(!on_pool(pool, 6));
        Ok(())
    })
    .unwrap();
    with_search_pool(3, "tblastx", |pool| {
        assert!(on_pool(pool, 6));
        assert!(!on_pool(pool, 1));
        assert!(!on_pool(pool, 0));
        Ok(())
    })
    .unwrap();
    assert_eq!(masked_frame_count(&frames_of(&[b"ACGT".to_vec()])), 4);
    assert_eq!(masked_frame_count(&frames_of(&[b"ACGTAC".to_vec()])), 6);
    assert_eq!(masked_frame_count(&frames_of(&[b"AC".to_vec()])), 0);
}

/// With the switch on, the frames are masked by the threads of the search pool, and by more
/// than one of them: the probe that the parallel loop calls once per frame sees a rayon worker
/// of a pool of three threads, and the frames are spread over at least two workers. (Each probe
/// call waits until a second worker has called it (the first call that waits 10 s in vain ends
/// all waiting), so that the spread does not depend on the scheduler; a serial loop has no
/// second worker, and fails after those 10 s.) With the switch off, or on a pool of one
/// thread, the probe is never called.
#[cfg(feature = "parallel")]
#[test]
fn frames_are_masked_on_the_pool_threads() {
    use std::collections::BTreeSet;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    use std::sync::Mutex;
    use std::time::{Duration, Instant};

    let seqs: Vec<Vec<u8>> = (0..5).map(|i| b"CAGCAGTTGACCA".repeat(20 + i)).collect();
    let seg = SegMasker::new(12, 2.2, 2.5);
    let mut want = frames_of(&seqs);
    serial::<false>(&mut want, &seg);
    let frame_count = masked_frame_count(&want);
    with_search_pool(3, "tblastx", |pool| {
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
        // The block drops the `workers` guard (it holds the lock the probe takes) before the
        // next run.
        {
            let mut got = frames_of(&seqs);
            mask_with_probe::<false, _>(&mut got, &seg, pool, true, &probe);
            assert_eq!(masked(&got), masked(&want));
            assert_eq!(calls.load(Ordering::SeqCst), frame_count);
            assert_eq!(
                in_pool_of_three.load(Ordering::SeqCst),
                frame_count,
                "every frame is masked by a worker of the 3-thread pool"
            );
            let workers = workers.lock().unwrap();
            assert!(
                workers.len() >= 2 && workers.iter().all(|&index| index < 3),
                "workers {workers:?}"
            );
        }
        // Not acted on: switch off.
        calls.store(0, Ordering::SeqCst);
        let mut got = frames_of(&seqs);
        mask_with_probe::<false, _>(&mut got, &seg, pool, false, &probe);
        assert_eq!(masked(&got), masked(&want));
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        Ok(())
    })
    .unwrap();
    with_search_pool(1, "tblastx", |pool| {
        let calls = AtomicUsize::new(0);
        let probe = || {
            calls.fetch_add(1, Ordering::SeqCst);
        };
        let mut got = frames_of(&seqs);
        mask_with_probe::<false, _>(&mut got, &seg, pool, true, &probe);
        assert_eq!(masked(&got), masked(&want));
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        Ok(())
    })
    .unwrap();
}
