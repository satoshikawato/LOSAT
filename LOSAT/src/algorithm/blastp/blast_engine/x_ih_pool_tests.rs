use super::super::BlastpUngappedData;
use super::*;

fn hit(i: usize) -> InitHSP {
    let i = i as i32;
    InitHSP {
        q_seed_absolute: i,
        s_seed: i.wrapping_mul(7),
        ungapped_data: BlastpUngappedData {
            q_start: i - 3,
            s_start: i.wrapping_mul(7) - 3,
            length: 11 + i % 5,
            score: 20 + i % 13,
        },
    }
}

fn same(a: &[InitHSP], b: &[InitHSP]) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(x, y)| {
            x.q_seed_absolute == y.q_seed_absolute
                && x.s_seed == y.s_seed
                && x.ungapped_data.q_start == y.ungapped_data.q_start
                && x.ungapped_data.s_start == y.ungapped_data.s_start
                && x.ungapped_data.length == y.ungapped_data.length
                && x.ungapped_data.score == y.ungapped_data.score
        })
}

/// One subject's hits pushed as `blastp_save_init_hsp` pushes them with the switch on.
fn push_subject(pool: &Pool, init_hsps: &mut Vec<InitHSP>, count: usize, seed: usize) -> usize {
    let mut pooled = 0;
    for i in 0..count {
        if init_hsps.len() == init_hsps.capacity() && init_hsps.capacity() >= OWN_LEN {
            grow(pool, init_hsps);
            pooled += 1;
        }
        init_hsps.push(hit(seed + i));
    }
    pooled
}

/// The end of a subject as [`at_subject_end`] ends it with the switch on.
fn end_subject(pool: &Pool, init_hsps: &mut Vec<InitHSP>) {
    if init_hsps.capacity() > OWN_LEN {
        give_back(pool, init_hsps);
    } else {
        init_hsps.clear();
    }
}

// The hits of every subject equal those of a vector that grows by itself, in the same order;
// after each subject the scratch holds an array of NCBI's first size again. At one thread the
// pool then holds one array, the largest a subject needed, as the scratch held before.
#[test]
fn hits_through_the_pool_equal_the_scratch_vector() {
    let pool = Pool::new(Vec::new());
    let mut scratch = Vec::with_capacity(OWN_LEN);
    let counts = [
        0, 1, 99, 100, 101, 200, 201, 5, 1_000, 3_333, 12_800, 0, 25_601, 7, 400, 100, 12_801,
    ];
    let mut largest = 0;
    for (n, &count) in counts.iter().enumerate() {
        let pooled = push_subject(&pool, &mut scratch, count, n * 100_000);
        let reference: Vec<InitHSP> = (0..count).map(|i| hit(n * 100_000 + i)).collect();
        assert!(same(&scratch, &reference), "subject {n} ({count} hits)");
        assert_eq!(pooled > 0, count > OWN_LEN, "subject {n} ({count} hits)");
        largest = largest.max(scratch.capacity());
        end_subject(&pool, &mut scratch);
        assert_eq!(scratch.capacity(), OWN_LEN, "subject {n}");
        assert!(scratch.is_empty());
        let arrays = pool.lock().unwrap();
        assert!(arrays.iter().all(Vec::is_empty), "subject {n}");
        let held: Vec<usize> = arrays.iter().map(Vec::capacity).collect();
        let expected = if largest > OWN_LEN {
            vec![largest]
        } else {
            vec![]
        };
        assert_eq!(held, expected, "subject {n}");
    }
    // 12,800 hits fit 12,800 places (100 doubled 7 times); 25,601 need 51,200.
    assert_eq!(largest, 51_200);
    release_pool(&pool);
    assert!(pool.lock().unwrap().is_empty());
    assert_eq!(push_subject(&pool, &mut scratch, 450, 0), 3);
    assert_eq!(scratch.capacity(), 800);
}

// Threads that share one pool: every subject's hits equal the reference, and when all are done
// every scratch has taken an array of the first size back (none is left in the pool) and the
// pool holds no array twice.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn threads_sharing_the_pool_keep_their_own_hits() {
    let pool = Pool::new(Vec::new());
    std::thread::scope(|scope| {
        for t in 0..8usize {
            let pool = &pool;
            scope.spawn(move || {
                let mut scratch = Vec::with_capacity(OWN_LEN);
                let mut state = 0x9E37_79B9_7F4A_7C15u64 ^ t as u64;
                for n in 0..300usize {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    let count = match state % 10 {
                        0 => (state >> 8) as usize % 20_000,
                        1..=3 => (state >> 8) as usize % 600,
                        _ => (state >> 8) as usize % 120,
                    };
                    let seed = (t * 1_000 + n) * 100_000;
                    push_subject(pool, &mut scratch, count, seed);
                    let reference: Vec<InitHSP> = (0..count).map(|i| hit(seed + i)).collect();
                    assert!(same(&scratch, &reference), "thread {t} subject {n}");
                    end_subject(pool, &mut scratch);
                    assert_eq!(scratch.capacity(), OWN_LEN);
                }
            });
        }
    });
    let arrays = pool.lock().unwrap();
    assert_eq!(arrays.iter().filter(|a| a.capacity() <= OWN_LEN).count(), 0);
    assert!(arrays.iter().all(Vec::is_empty));
    let mut starts: Vec<*const InitHSP> = arrays.iter().map(|a| a.as_ptr()).collect();
    starts.sort();
    starts.dedup();
    assert_eq!(starts.len(), arrays.len());
}
