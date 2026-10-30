//! EXPERIMENT (LOSAT_X_BPTREEPOOL): the interval trees of the preliminary stage live in one pool
//! shared by the threads, not in each thread's subject scratch.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3826-3827
//! ```c
//!     tree = Blast_IntervalTreeInit(0, query->length+1,
//!                                   0, subject->length+1);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4107
//! ```c
//!    tree = Blast_IntervalTreeFree(tree);
//! ```
//! NCBI makes a tree per subject and frees it when the subject's gapped stage ends. The subject
//! scratch (`BlastpSubjectScratch::preliminary_interval_tree`) kept one tree per thread, and each
//! tree kept the largest arrays any of its subjects needed: the search held the sum over the
//! threads of each thread's largest tree. On UniRef50 at 32 threads that was 34.8 MiB of capacity,
//! of which the trees in use held 1.72 MiB at most after the first 30 s (v4-D12 `DECISIONS` 73).
//!
//! With the switch on, the scratch holds an empty tree, which answers "no container" as a tree
//! without HSPs does. Before a subject's first HSP ([`before_add`]) the scratch takes the pool's
//! largest tree, or a new one, with the same ranges; at the subject's end ([`at_subject_end`]) it
//! goes back. The pool holds the trees the threads use at the same time, not each thread's
//! largest. The tree's nodes and the answers it gives are those of the scratch's own tree; only
//! where the arrays live changes. At one thread the pool holds the thread's largest tree, as the
//! scratch did. After a search's preliminary stage the pool's trees are freed ([`release`]), as
//! the scratches' trees were when the subject loop dropped its scratches.

use std::sync::{Mutex, OnceLock, PoisonError};

use crate::algorithm::blastn::interval_tree::BlastIntervalTree;

/// The trees no subject uses.
struct Pool {
    /// Trees that have listed HSPs, with their arrays.
    held: Vec<BlastIntervalTree>,
    /// Empty trees, one per scratch whose subject holds a tree of `held`.
    empty: Vec<BlastIntervalTree>,
}

static POOL: Mutex<Pool> = Mutex::new(Pool::new());

impl Pool {
    const fn new() -> Self {
        Self {
            held: Vec::new(),
            empty: Vec::new(),
        }
    }
}

// No NCBI counterpart: reads the LOSAT_X_BPTREEPOOL switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPTREEPOOL): strict set (on unless `0`).
fn x_bptreepool() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPTREEPOOL"))
}

// No NCBI counterpart: where the arrays live (see the module header); it does not change any value
// NCBI computes.
/// Called before an HSP is added to the scratch's tree: with the switch on, an empty tree is
/// exchanged for one of the pool.
#[inline]
pub(super) fn before_add(tree: &mut BlastIntervalTree) {
    if tree.x_has_no_arrays() && x_bptreepool() {
        take(&POOL, tree);
    }
}

// No NCBI counterpart: gives the subject's tree back (see the module header); it does not change
// any value NCBI computes.
/// At the end of a subject, with its tree no longer read: a tree of the pool goes back, and the
/// scratch takes an empty one.
#[inline]
pub(super) fn at_subject_end(tree: &mut BlastIntervalTree) {
    if !tree.x_has_no_arrays() && x_bptreepool() {
        give_back(&POOL, tree);
    }
}

// No NCBI counterpart: frees the trees no subject uses; it does not change any value NCBI computes.
/// After a search's preliminary stage (no subject scan runs any more): the pool's trees are freed.
pub(super) fn release() {
    if x_bptreepool() {
        release_pool(&POOL);
    }
}

/// The largest tree goes to the subject, so smaller trees stay small and the pool's arrays are
/// those of the largest subjects at the same time. `tree` is empty and holds the subject's ranges.
#[cold]
#[inline(never)]
fn take(pool: &Mutex<Pool>, tree: &mut BlastIntervalTree) {
    let held = {
        let mut trees = pool.lock().unwrap_or_else(PoisonError::into_inner);
        trees
            .held
            .iter()
            .enumerate()
            .max_by_key(|(_, held)| held.x_array_bytes())
            .map(|(index, _)| index)
            .map(|index| trees.held.swap_remove(index))
    };
    let mut held = held.unwrap_or_else(BlastIntervalTree::x_new_with_arrays);
    held.x_reset_to_ranges_of(tree);
    let empty = std::mem::replace(tree, held);
    pool.lock()
        .unwrap_or_else(PoisonError::into_inner)
        .empty
        .push(empty);
}

/// The subject's tree goes to the pool; `tree` becomes an empty one from the pool, or a new one.
#[cold]
#[inline(never)]
fn give_back(pool: &Mutex<Pool>, tree: &mut BlastIntervalTree) {
    let empty = pool
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .empty
        .pop();
    let empty = empty.unwrap_or_else(|| BlastIntervalTree::new(0, 1, 0, 1));
    let used = std::mem::replace(tree, empty);
    pool.lock()
        .unwrap_or_else(PoisonError::into_inner)
        .held
        .push(used);
}

fn release_pool(pool: &Mutex<Pool>) {
    let trees = std::mem::replace(
        &mut *pool.lock().unwrap_or_else(PoisonError::into_inner),
        Pool::new(),
    );
    drop(trees);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastn::interval_tree::{IndexMethod, TreeHsp};

    fn hsp(start: i32, end: i32, score: i32) -> TreeHsp {
        TreeHsp {
            query_offset: start,
            query_end: end,
            subject_offset: start,
            subject_end: end,
            score,
            query_frame: 0,
            query_length: 2000,
            query_context_offset: 0,
            subject_frame_sign: 1,
        }
    }

    fn add(tree: &mut BlastIntervalTree, value: TreeHsp) {
        tree.add_hsp(value, 0, IndexMethod::QueryAndSubject);
    }

    // The tree a subject borrows answers as a tree of its own, with arrays left by a longer
    // subject of another range, and the next subject starts empty.
    #[test]
    fn x_borrowed_tree_answers_as_the_scratchs_own() {
        let pool = Mutex::new(Pool::new());
        let mut scratch_tree = BlastIntervalTree::new(0, 2000, 0, 1);
        let mut own = BlastIntervalTree::new(0, 2000, 0, 1);

        // A long subject grows a pool tree.
        scratch_tree.reset_with_bounds(0, 2000, 0, 90_000);
        take(&pool, &mut scratch_tree);
        for i in 0..400 {
            add(&mut scratch_tree, hsp(i * 4, i * 4 + 60, 50));
        }
        assert!(scratch_tree.node_count() > 400);
        give_back(&pool, &mut scratch_tree);
        assert!(scratch_tree.x_has_no_arrays());
        assert_eq!(scratch_tree.node_count(), 1);

        // A short subject takes the same tree: nothing of the long subject answers.
        scratch_tree.reset_with_bounds(0, 2000, 0, 900);
        own.reset_with_bounds(0, 2000, 0, 900);
        assert!(!scratch_tree.contains_hsp(&hsp(100, 150, 10), 0, 0));
        take(&pool, &mut scratch_tree);
        assert!(!scratch_tree.x_has_no_arrays());
        assert_eq!(scratch_tree.node_count(), 1);
        assert!(!scratch_tree.contains_hsp(&hsp(100, 150, 10), 0, 0));
        for tree in [&mut scratch_tree, &mut own] {
            add(tree, hsp(100, 200, 50));
            add(tree, hsp(300, 500, 80));
        }
        for probe in [
            hsp(120, 180, 10),
            hsp(120, 180, 60),
            hsp(320, 480, 80),
            hsp(90, 180, 10),
            hsp(600, 700, 10),
        ] {
            assert_eq!(
                scratch_tree.contains_hsp(&probe, 0, 0),
                own.contains_hsp(&probe, 0, 0),
                "{probe:?}"
            );
        }
        assert_eq!(scratch_tree.node_count(), own.node_count());
    }

    // The pool keeps one tree per subject held at the same time and hands out the largest first.
    #[test]
    fn x_pool_holds_trees_in_use_at_once_largest_first() {
        let pool = Mutex::new(Pool::new());
        let mut a = BlastIntervalTree::new(0, 2000, 0, 1);
        let mut b = BlastIntervalTree::new(0, 2000, 0, 1);
        a.reset_with_bounds(0, 2000, 0, 2000);
        b.reset_with_bounds(0, 2000, 0, 2000);
        take(&pool, &mut a);
        take(&pool, &mut b);
        for i in 0..300 {
            add(&mut a, hsp(i * 4, i * 4 + 60, 50));
        }
        let large = a.x_array_bytes();
        let small = b.x_array_bytes();
        assert!(large > small);
        give_back(&pool, &mut b);
        give_back(&pool, &mut a);
        assert_eq!(pool.lock().unwrap().held.len(), 2);
        assert_eq!(pool.lock().unwrap().empty.len(), 0);

        take(&pool, &mut b);
        assert_eq!(b.x_array_bytes(), large);
        give_back(&pool, &mut b);

        release_pool(&pool);
        assert_eq!(pool.lock().unwrap().held.len(), 0);
    }
}
