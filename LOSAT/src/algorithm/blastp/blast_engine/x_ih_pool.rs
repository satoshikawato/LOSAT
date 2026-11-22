//! EXPERIMENT (LOSAT_X_BPIHPOOL): the arrays of initial hits that outgrow NCBI's first array size
//! live in one pool shared by the threads, not in each thread's subject scratch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:221-224
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:336-341
//! NCBI resets a thread's list for each subject (`BlastInitHitListReset`, c++/src/algo/blast/core/blast_engine.c:491) and
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
//! one thread the pool holds the thread's largest array, as the scratch did. After a search's
//! preliminary stage the pool's arrays are freed ([`release`]), as the scratches' arrays were when
//! the subject loop dropped its scratches.
//!
//! `LOSAT_TIMING=1` prints `[TIMING] x_bpihpool: pooled=<n>` per search: the arrays filled at
//! [`OWN_LEN`] or more that continued through the pool.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

use super::InitHSP;

/// NCBI's first array size, `MIN_INIT_HITLIST_SIZE`.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:138-139
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

// No NCBI counterpart: frees the arrays no subject uses; it does not change any value NCBI computes.
/// After a search's preliminary stage (no subject scan runs any more): the pool's arrays are freed.
pub(super) fn release() {
    if x_bpihpool() {
        release_pool(&POOL);
    }
}

fn release_pool(pool: &Pool) {
    let arrays = std::mem::take(&mut *pool.lock().unwrap_or_else(PoisonError::into_inner));
    drop(arrays);
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
#[path = "x_ih_pool_tests.rs"]
mod tests;
