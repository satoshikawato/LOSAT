//! The working memory of a search's parallel work, kept one value per thread.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback_mt_priv.c:218-226
//! ```c
//!     for (i = 0; i < array->num_elems; i++) {
//!         status =
//!            BLAST_GapAlignSetUp(program, seqsrc, score_options, eff_len_options,
//!               ext_options, hit_options, query_info, sbp,
//!               &array->tld[i]->score_params,
//!               &array->tld[i]->ext_params,
//!               &array->tld[i]->hit_params,
//!               &array->tld[i]->eff_len_params,
//!               &array->tld[i]->gap_align);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1533-1537
//! ```c
//! #ifdef _OPENMP
//!             tid = omp_get_thread_num();
//! #endif
//!             seqsrc = thread_data->tld[tid]->seqsrc;
//!             gap_align = thread_data->tld[tid]->gap_align;
//! ```
//! NCBI makes one set of working memory (`gap_align` and its parts) for each thread of a search
//! and finds it by the thread's number. `ScratchPool` keeps the same shape for work that Rayon
//! schedules: one slot per thread of the search's pool (`rayon::current_thread_index`), and one
//! for a caller outside it. A piece of work takes its thread's value, or makes one when the slot
//! is empty (a thread that waits in a nested parallel stage and runs another piece of the same
//! kind meanwhile), and gives it back to its slot when the slot is empty; otherwise the value is
//! dropped. So a search holds at most one value per thread, plus those of nested work in flight,
//! however Rayon cuts the work (a `map_init` makes one per piece Rayon cuts, at least twice the
//! threads and more on every steal).
//!
//! No NCBI counterpart for what a value keeps between pieces of work (`Settle`): NCBI keeps the
//! working memory of a thread at the largest size it has needed. Here a value that comes back
//! puts every buffer (or set of buffers sized together) that holds more than `STEADY_BYTES` back
//! to its steady size, so that the memory a search keeps does not grow with the largest piece of
//! work each thread has seen. Memory placement only: a value holds what a new value holds for
//! every piece of work, and what is kept between pieces is never read before it is written (the
//! work a value does depends on its arguments only, as the pieces of a `map_init` already showed,
//! which take a new value or a used one depending on the schedule).

use std::sync::{Mutex, MutexGuard, PoisonError};

/// No NCBI counterpart: the one steady size of the pooled values' buffers, in bytes.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:61,75
/// ```c
/// #define	CHUNKSIZE	2097152
/// ...
///    Int4	chunksize = MAX(CHUNKSIZE, length + length/3);
/// ```
/// It is the size from which a traceback scratch already gives its trace back after a call
/// (`xdrop_simd::X_POOL_FROM`, the four-bit rows of one NCBI chunk): the same bound now decides
/// what every buffer of a pooled value keeps between pieces of work, so one rule, not one
/// threshold per buffer, sets the memory a thread keeps.
pub(crate) const STEADY_BYTES: usize = crate::utils::xdrop_simd::X_POOL_FROM;

/// No NCBI counterpart: what a value of a `ScratchPool` keeps between pieces of work.
pub(crate) trait Settle {
    /// Called when the value comes back to its pool, after its piece of work: every buffer (or
    /// set of buffers sized together) that holds more than `STEADY_BYTES` goes back to its
    /// steady size. It changes only what is kept, never what a piece of work reads.
    fn settle(&mut self);
}

impl Settle for () {
    fn settle(&mut self) {}
}

impl<A: Settle, B: Settle> Settle for (A, B) {
    fn settle(&mut self) {
        self.0.settle();
        self.1.settle();
    }
}

/// No NCBI counterpart: whether a buffer of `bytes` passes the steady size (`STEADY_BYTES`).
#[inline]
pub(crate) fn passes(bytes: usize) -> bool {
    bytes > STEADY_BYTES
}

/// No NCBI counterpart: the bytes `v` holds (its capacity).
#[inline]
pub(crate) fn vec_bytes<T>(v: &Vec<T>) -> usize {
    v.capacity().saturating_mul(std::mem::size_of::<T>())
}

/// No NCBI counterpart: `Settle` for a plain buffer that each piece of work empties before it
/// fills it: past `STEADY_BYTES`, its contents go and its capacity becomes `STEADY_BYTES`.
pub(crate) fn settle_vec<T>(v: &mut Vec<T>) {
    if passes(vec_bytes(v)) {
        v.clear();
        v.shrink_to(STEADY_BYTES / std::mem::size_of::<T>().max(1));
    }
}

// No NCBI counterpart: a slot on its own cache line, so that two threads giving values back to
// neighbouring slots do not share one.
#[repr(align(128))]
struct Slot<T>(Mutex<Option<T>>);

/// The working memory of one kind for one search, one value per thread (module header).
pub(crate) struct ScratchPool<T: Settle> {
    /// One slot per thread of the search's pool, then one for a caller outside it.
    slots: Box<[Slot<T>]>,
}

impl<T: Settle> ScratchPool<T> {
    /// No NCBI counterpart: a pool for a search on `threads` threads (`SearchPool::threads`).
    pub(crate) fn new(threads: usize) -> Self {
        Self {
            slots: (0..=threads).map(|_| Slot(Mutex::new(None))).collect(),
        }
    }

    // No NCBI counterpart: the calling thread's slot (`tid` above): its index in the Rayon pool
    // it runs in, or the last slot outside a pool (or on a thread of a larger pool).
    fn here(&self) -> usize {
        let outside = self.slots.len() - 1;
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if let Some(index) = rayon::current_thread_index() {
            return index.min(outside);
        }
        outside
    }

    // No NCBI counterpart: a slot's lock; a slot holds plain memory, so a poisoned lock carries
    // no meaning.
    fn slot(&self, index: usize) -> MutexGuard<'_, Option<T>> {
        self.slots[index]
            .0
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
    }

    /// No NCBI counterpart: the calling thread's value, made with `make` when its slot is empty;
    /// it goes back to the slot when the lease is dropped.
    pub(crate) fn lease(&self, make: impl FnOnce() -> T) -> Lease<'_, T> {
        let slot = self.here();
        // The lock is released before a value is made.
        let held = self.slot(slot).take();
        Lease {
            pool: self,
            slot,
            value: Some(held.unwrap_or_else(make)),
        }
    }

    /// No NCBI counterpart: `work` on the calling thread's value (`lease`).
    pub(crate) fn with<R>(&self, make: impl FnOnce() -> T, work: impl FnOnce(&mut T) -> R) -> R {
        let mut lease = self.lease(make);
        work(&mut lease)
    }

    /// No NCBI counterpart: a lease taken on first use (`Deferred::get`), on the thread that uses
    /// it, for work that is handed to another thread before it starts.
    pub(crate) fn deferred(&self) -> Deferred<'_, T> {
        Deferred {
            pool: self,
            lease: None,
        }
    }

    /// Test only: how many slots hold a value.
    #[cfg(test)]
    fn kept(&self) -> usize {
        (0..self.slots.len())
            .filter(|&i| self.slot(i).is_some())
            .count()
    }
}

/// No NCBI counterpart: a value of a `ScratchPool`, held by one piece of work.
pub(crate) struct Lease<'p, T: Settle> {
    pool: &'p ScratchPool<T>,
    slot: usize,
    value: Option<T>,
}

impl<T: Settle> std::ops::Deref for Lease<'_, T> {
    type Target = T;
    fn deref(&self) -> &T {
        self.value
            .as_ref()
            .expect("a lease holds its value until it is dropped")
    }
}

impl<T: Settle> std::ops::DerefMut for Lease<'_, T> {
    fn deref_mut(&mut self) -> &mut T {
        self.value
            .as_mut()
            .expect("a lease holds its value until it is dropped")
    }
}

impl<T: Settle> Drop for Lease<'_, T> {
    // No NCBI counterpart: the value settles (`Settle`) and goes back to its slot when the slot is
    // empty (a nested lease of the same thread may have given its value back first); otherwise it
    // is dropped, outside the lock. A piece of work that unwinds drops its value.
    fn drop(&mut self) {
        let Some(mut value) = self.value.take() else {
            return;
        };
        if std::thread::panicking() {
            return;
        }
        value.settle();
        let spare = {
            let mut slot = self.pool.slot(self.slot);
            if slot.is_none() {
                *slot = Some(value);
                None
            } else {
                Some(value)
            }
        };
        drop(spare);
    }
}

/// No NCBI counterpart: a lease taken on its first use (`ScratchPool::deferred`). Work that is
/// handed to a thread of the pool (a spawned helper) takes the value of the thread that runs it,
/// and takes none when it does no work.
pub(crate) struct Deferred<'p, T: Settle> {
    pool: &'p ScratchPool<T>,
    lease: Option<Lease<'p, T>>,
}

impl<T: Settle> Deferred<'_, T> {
    /// No NCBI counterpart: the value, leased on this thread at the first call (`lease`).
    pub(crate) fn get(&mut self, make: impl FnOnce() -> T) -> &mut T {
        let pool = self.pool;
        self.lease.get_or_insert_with(|| pool.lease(make))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// The values made, and under one lock the values alive, the leases in flight and the most
    /// values alive beyond the leases in flight.
    #[derive(Default)]
    struct Counts {
        made: AtomicUsize,
        live: Mutex<[usize; 3]>,
    }

    impl Counts {
        fn live(&self) -> MutexGuard<'_, [usize; 3]> {
            self.live.lock().unwrap()
        }
        fn alive(&self) -> usize {
            self.live()[0]
        }
        fn made(&self) -> usize {
            self.made.load(Ordering::SeqCst)
        }
    }

    /// A value that counts itself (`Counts`), with a buffer to settle.
    struct Probe<'c> {
        id: usize,
        buffer: Vec<u64>,
        counts: &'c Counts,
    }

    impl<'c> Probe<'c> {
        fn make(counts: &'c Counts) -> Self {
            let mut live = counts.live();
            live[0] += 1;
            live[2] = live[2].max(live[0].saturating_sub(live[1]));
            Probe {
                id: counts.made.fetch_add(1, Ordering::SeqCst),
                buffer: Vec::new(),
                counts,
            }
        }
    }

    impl Drop for Probe<'_> {
        fn drop(&mut self) {
            self.counts.live()[0] -= 1;
        }
    }

    impl Settle for Probe<'_> {
        fn settle(&mut self) {
            settle_vec(&mut self.buffer);
        }
    }

    // A buffer past the steady size keeps the steady capacity and loses its contents; one within
    // it keeps both.
    #[test]
    fn a_buffer_settles_to_the_steady_size_only_past_it() {
        let mut small: Vec<u64> = Vec::with_capacity(1000);
        small.push(7);
        settle_vec(&mut small);
        assert_eq!((small.len(), small.capacity()), (1, 1000));
        let mut large: Vec<u64> = vec![3; STEADY_BYTES / 8 + 1];
        settle_vec(&mut large);
        assert!(large.is_empty());
        assert!(vec_bytes(&large) <= STEADY_BYTES, "{}", vec_bytes(&large));
        assert!(
            vec_bytes(&large) >= STEADY_BYTES / 2,
            "{}",
            vec_bytes(&large)
        );
    }

    // Outside a Rayon pool the caller's slot is the last one: its value comes back to it, settled,
    // and the next lease takes the same value.
    #[test]
    fn a_caller_outside_a_pool_gets_its_value_back_settled() {
        let counts = Counts::default();
        let pool: ScratchPool<Probe<'_>> = ScratchPool::new(4);
        let first = pool.with(
            || Probe::make(&counts),
            |probe| {
                probe.buffer.resize(4 * STEADY_BYTES / 8, 1);
                probe.id
            },
        );
        let (again, kept) = pool.with(
            || Probe::make(&counts),
            |probe| (probe.id, vec_bytes(&probe.buffer)),
        );
        assert_eq!(again, first);
        assert!(kept <= STEADY_BYTES, "{kept}");
        assert_eq!((counts.made(), pool.kept()), (1, 1));
        drop(pool);
        assert_eq!(counts.alive(), 0);
    }

    // Two leases at once on one thread (nested work): the inner one makes a second value; the
    // value that comes back first takes the slot, the other is dropped.
    #[test]
    fn nested_leases_make_a_value_each_and_keep_one() {
        let counts = Counts::default();
        let pool: ScratchPool<Probe<'_>> = ScratchPool::new(1);
        {
            let outer = pool.lease(|| Probe::make(&counts));
            let inner = pool.lease(|| Probe::make(&counts));
            assert_ne!(outer.id, inner.id);
            assert_eq!(counts.alive(), 2);
        }
        assert_eq!((counts.made(), counts.alive(), pool.kept()), (2, 1, 1));
    }

    // A deferred lease that is never used takes nothing; a used one leases once.
    #[test]
    fn an_unused_deferred_lease_takes_nothing() {
        let counts = Counts::default();
        let pool: ScratchPool<Probe<'_>> = ScratchPool::new(2);
        drop(pool.deferred());
        assert_eq!(counts.made(), 0);
        let mut used = pool.deferred();
        let id = used.get(|| Probe::make(&counts)).id;
        assert_eq!(used.get(|| unreachable!("leased once")).id, id);
        drop(used);
        assert_eq!((counts.made(), pool.kept()), (1, 1));
    }

    // Work that unwinds drops its value instead of giving it back.
    #[test]
    fn a_lease_that_unwinds_drops_its_value() {
        let counts = Counts::default();
        let pool: ScratchPool<Probe<'_>> = ScratchPool::new(1);
        let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            pool.with(
                || Probe::make(&counts),
                |_| panic!("a failed piece of work"),
            )
        }));
        assert!(unwound.is_err());
        assert_eq!((counts.alive(), pool.kept()), (0, 0));
    }

    // On a pool of N threads: pieces of work of uneven size (so that the threads steal from each
    // other) make at most N values, and the values kept are settled. With nested parallel work
    // in each piece, the values alive are never more than the leases in flight plus one per
    // thread, and at most N stay when the work is done.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn a_search_pool_holds_one_value_per_thread_and_nested_work_in_flight() {
        use rayon::prelude::*;
        let work = |n: u64| (0..n).fold(n, |acc, x| acc.wrapping_mul(31).wrapping_add(x));
        for threads in 1..=8usize {
            let rayon_pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let counts = Counts::default();
            let pool: ScratchPool<Probe<'_>> = ScratchPool::new(threads);
            let make = || Probe::make(&counts);
            rayon_pool.install(|| {
                (0..2_000u64).into_par_iter().for_each(|i| {
                    pool.with(make, |probe| {
                        // up to twice the steady size
                        probe
                            .buffer
                            .resize((i % 5) as usize * (STEADY_BYTES / 16), 0);
                        work(if i % 97 == 5 { 200_000 } else { 50 + i % 3_000 })
                    });
                })
            });
            assert!(
                (1..=threads).contains(&counts.made()),
                "{threads} threads: {} made",
                counts.made()
            );
            assert!(pool.kept() <= threads);
            for i in 0..=threads {
                if let Some(probe) = pool.slot(i).as_ref() {
                    assert!(vec_bytes(&probe.buffer) <= STEADY_BYTES);
                }
            }
            // `lease` counts a lease in flight from before it takes a value until after it gave
            // the value back (`live[1]`).
            let lease = |f: &(dyn Fn() + Sync)| {
                counts.live()[1] += 1;
                pool.with(make, |_| f());
                counts.live()[1] -= 1;
            };
            counts.live()[2] = 0;
            rayon_pool.install(|| {
                (0..300u64).into_par_iter().for_each(|i| {
                    lease(&|| {
                        (0..4u64).into_par_iter().for_each(|j| {
                            lease(&|| {
                                work(50 + (i * 4 + j) % 2_000);
                            })
                        })
                    })
                })
            });
            let excess = counts.live()[2];
            assert!(
                excess <= threads,
                "{threads} threads: {excess} beyond the leases"
            );
            assert!(pool.kept() <= threads);
            assert_eq!(counts.alive(), pool.kept());
            drop(pool);
            assert_eq!(counts.alive(), 0);
        }
    }
}
