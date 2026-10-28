//! EXPERIMENT (LOSAT_X_BPIHPOOL): the arrays of initial hits that outgrow NCBI's first array size
//! live in one pool shared by the threads, not in each thread's subject scratch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:221-224
//! ```c
//!     init_hitlist->allocated = MIN_INIT_HITLIST_SIZE;
//!
//!     init_hitlist->init_hsp_array = (BlastInitHSP *)
//!         malloc(MIN_INIT_HITLIST_SIZE * sizeof(BlastInitHSP));
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:336-341
//! ```c
//!     if (num >= num_avail) {
//!         if (init_hitlist->do_not_reallocate)
//!             return FALSE;
//!         num_avail *= 2;
//!         match_array = (BlastInitHSP *)
//!             realloc(match_array, num_avail * sizeof(BlastInitHSP));
//! ```
//! NCBI resets a thread's list for each subject (`BlastInitHitListReset`, blast_engine.c:491) and
//! keeps its array, so each thread holds the largest array any of its subjects needed. The subject
//! scratch (`BlastpSubjectScratch::init_hsps`, cleared, never shrunk) did the same, and each
//! worker's malloc arena also kept the arrays its vector outgrew: the search held the sum over the
//! threads of each thread's largest array. On UniRef50 at 32 threads that was 53 MiB of capacity,
//! of which a subject used 0.2 MiB at most at any moment (v4-D12 `DECISIONS` 62).
//!
//! With the switch on, the scratch holds an array of NCBI's first size ([`OWN_LEN`]). A subject
//! that fills it continues in an array of twice the length from the pool (the smallest that
//! holds it, or a new one), and the scratch's array goes to the pool; a pool array that fills is
//! exchanged for a larger one from the pool in the same way, or, when the pool has none, doubled
//! in place as NCBI's `realloc` doubles it. At the subject's end ([`at_subject_end`]) its array
//! goes back to the pool and the scratch takes an array of the first size from it. The pool holds
//! the arrays the threads use at the same time, not each thread's largest. The hits, their order
//! and every value are those of the scratch's own vector; only where the arrays live changes. At
//! one thread the pool holds the thread's largest array, as the scratch did.
//!
//! `LOSAT_TIMING=1` prints `[TIMING] x_bpihpool: pooled=<n>` per search: the arrays filled at
//! [`OWN_LEN`] or more that continued through the pool.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

use super::InitHSP;

/// NCBI's first array size, `MIN_INIT_HITLIST_SIZE`.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:138-139
/// ```c
/// /** Minimal size of an array of initial word hits, allocated up front. */
/// #define MIN_INIT_HITLIST_SIZE 100
/// ```
pub(super) const OWN_LEN: usize = 100;

/// The arrays no subject uses: those of [`OWN_LEN`] (one per scratch whose subject holds a pool
/// array) and the larger ones.
type Pool = Mutex<Vec<Vec<InitHSP>>>;

static POOL: Pool = Mutex::new(Vec::new());
static POOLED: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: reads the LOSAT_X_BPIHPOOL switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPIHPOOL): strict set (on unless `0`).
pub(super) fn x_bpihpool() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPIHPOOL"))
}

// No NCBI counterpart: the room of the scratch's array; it does not change any value NCBI
// computes.
/// The array a new subject scratch holds: NCBI's first size with the switch on, else empty.
pub(super) fn own_array() -> Vec<InitHSP> {
    if x_bpihpool() {
        Vec::with_capacity(OWN_LEN)
    } else {
        Vec::new()
    }
}

/// The arrays filled at [`OWN_LEN`] or more that continued through the pool (for the timing line).
pub(super) fn x_pooled() -> u64 {
    POOLED.load(Ordering::Relaxed)
}

// No NCBI counterpart: where the array lives (blast_extend.c:336-341 above); it does not change
// any value NCBI computes.
/// Called by `blastp_save_init_hsp` when `init_hsps` is full: with the switch on and the array at
/// least [`OWN_LEN`] long, it continues in an array of twice the length; otherwise `push` grows
/// the vector as before.
#[cold]
#[inline(never)]
pub(super) fn before_full_push(init_hsps: &mut Vec<InitHSP>) {
    if init_hsps.capacity() >= OWN_LEN && x_bpihpool() {
        grow(&POOL, init_hsps);
        POOLED.fetch_add(1, Ordering::Relaxed);
    }
}

// No NCBI counterpart: gives the subject's array back (see the module header); it does not change
// any value NCBI computes.
/// At the end of a subject, with its hits no longer read: an array larger than [`OWN_LEN`] goes
/// back to the pool and the scratch takes one of [`OWN_LEN`].
#[inline]
pub(super) fn at_subject_end(init_hsps: &mut Vec<InitHSP>) {
    if init_hsps.capacity() > OWN_LEN && x_bpihpool() {
        give_back(&POOL, init_hsps);
    }
}

/// `init_hsps` is full: it continues in an array of twice its length. The scratch's array (at most
/// [`OWN_LEN`]) and an outgrown pool array go to the pool; a pool array that the pool cannot
/// replace is doubled in place.
fn grow(pool: &Pool, init_hsps: &mut Vec<InitHSP>) {
    let need = init_hsps.len() * 2;
    let fit = {
        let mut arrays = pool.lock().unwrap_or_else(PoisonError::into_inner);
        arrays
            .iter()
            .enumerate()
            .filter(|(_, array)| array.capacity() >= need)
            .min_by_key(|(_, array)| array.capacity())
            .map(|(index, _)| index)
            .map(|index| arrays.swap_remove(index))
    };
    let next = match fit {
        Some(array) => array,
        None if init_hsps.capacity() > OWN_LEN => {
            init_hsps.reserve_exact(need - init_hsps.len());
            return;
        }
        None => Vec::with_capacity(need),
    };
    let mut outgrown = std::mem::replace(init_hsps, next);
    init_hsps.extend_from_slice(&outgrown);
    outgrown.clear();
    pool.lock()
        .unwrap_or_else(PoisonError::into_inner)
        .push(outgrown);
}

/// The subject's array (larger than [`OWN_LEN`]) goes back to the pool; `init_hsps` becomes an
/// array of [`OWN_LEN`] from the pool, or a new one.
fn give_back(pool: &Pool, init_hsps: &mut Vec<InitHSP>) {
    let mut used = std::mem::take(init_hsps);
    used.clear();
    let own = {
        let mut arrays = pool.lock().unwrap_or_else(PoisonError::into_inner);
        let own = arrays
            .iter()
            .position(|array| array.capacity() <= OWN_LEN)
            .map(|index| arrays.swap_remove(index));
        arrays.push(used);
        own
    };
    *init_hsps = own.unwrap_or_else(|| Vec::with_capacity(OWN_LEN));
}

#[cfg(test)]
mod tests {
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
    }

    // Threads that share one pool: every subject's hits equal the reference, and when all are done
    // every scratch has taken an array of the first size back (none is left in the pool) and the
    // pool holds no array twice.
    #[test]
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
}
