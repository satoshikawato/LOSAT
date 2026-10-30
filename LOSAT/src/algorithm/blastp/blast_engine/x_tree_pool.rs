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
//! With the switch on, the scratch's tree has room for most subjects ([`OWN_NODES`], [`OWN_ENTRIES`],
//! [`OWN_CELLS`]). Before an HSP that may not fit ([`before_add`]) the subject continues in a copy
//! in the pool's largest tree, or a new one; at the subject's end ([`at_subject_end`]) that tree
//! goes back and the scratch takes a tree of its own size. The pool holds the large trees the
//! threads use at the same time, not each thread's largest. The tree's nodes and the answers it
//! gives are those of the scratch's own tree; only where the arrays live changes. At one thread the
//! pool holds the thread's largest tree, as the scratch did. After a search's preliminary stage the
//! pool's trees are freed ([`release`]), as the scratches' trees were when the subject loop dropped
//! its scratches.

use std::sync::{Mutex, OnceLock, PoisonError};

use crate::algorithm::blastn::interval_tree::BlastIntervalTree;

// No NCBI counterpart: the room of the scratch's tree (memory only); past it the tree continues in
// the pool's, and it does not change any value NCBI computes.
/// Nodes of the scratch's tree (64 bytes each). NCBI's first array has 100.
const OWN_NODES: usize = 512;
/// The most nodes [`before_add`] leaves room for per HSP; an HSP that needs more only grows the
/// array, as NCBI's `realloc` does.
const NEW_NODES_PER_HSP: usize = 64;
/// Grid listings of the scratch's tree (8 bytes each).
const OWN_ENTRIES: usize = 4096;
/// Grid cells of the scratch's tree (4 bytes each).
const OWN_CELLS: usize = 1 << 14;

/// The trees no subject uses.
struct Pool {
    /// Trees larger than a scratch's.
    large: Vec<BlastIntervalTree>,
    /// Trees of a scratch's size, one per scratch whose subject holds a tree of `large`.
    own: Vec<BlastIntervalTree>,
}

static POOL: Mutex<Pool> = Mutex::new(Pool::new());

impl Pool {
    const fn new() -> Self {
        Self {
            large: Vec::new(),
            own: Vec::new(),
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

// No NCBI counterpart: the room of the scratch's tree; it does not change any value NCBI computes.
/// The tree a new subject scratch holds: with the switch on, with the room above.
pub(super) fn own_tree(query_span: i32) -> BlastIntervalTree {
    if x_bptreepool() {
        let mut tree = new_own();
        tree.reset_with_bounds(0, query_span, 0, 1);
        tree
    } else {
        BlastIntervalTree::new(0, query_span, 0, 1)
    }
}

fn new_own() -> BlastIntervalTree {
    BlastIntervalTree::x_with_room(OWN_NODES, OWN_ENTRIES, OWN_CELLS)
}

// No NCBI counterpart: where the arrays live (see the module header); it does not change any value
// NCBI computes.
/// Called before an HSP is added to the scratch's tree: with the switch on, a tree of a scratch's
/// size that the HSP may outgrow continues in one of the pool.
#[inline]
pub(super) fn before_add(tree: &mut BlastIntervalTree) {
    if tree.x_node_capacity() <= OWN_NODES
        && tree.x_may_outgrow(OWN_NODES, NEW_NODES_PER_HSP, OWN_ENTRIES, OWN_CELLS)
        && x_bptreepool()
    {
        take(&POOL, tree);
    }
}

// No NCBI counterpart: gives the subject's tree back (see the module header); it does not change
// any value NCBI computes.
/// At the end of a subject, with its tree no longer read: a tree larger than a scratch's goes back
/// to the pool, and the scratch takes one of its own size.
#[inline]
pub(super) fn at_subject_end(tree: &mut BlastIntervalTree) {
    if tree.x_node_capacity() > OWN_NODES && x_bptreepool() {
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
/// those of the largest subjects at the same time.
#[cold]
#[inline(never)]
fn take(pool: &Mutex<Pool>, tree: &mut BlastIntervalTree) {
    let large = {
        let mut trees = pool.lock().unwrap_or_else(PoisonError::into_inner);
        trees
            .large
            .iter()
            .enumerate()
            .max_by_key(|(_, large)| large.x_array_bytes())
            .map(|(index, _)| index)
            .map(|index| trees.large.swap_remove(index))
    };
    // Twice the nodes: a pool tree is told from a scratch's by the room of its node array.
    let mut large = large
        .unwrap_or_else(|| BlastIntervalTree::x_with_room(2 * OWN_NODES, OWN_ENTRIES, OWN_CELLS));
    large.x_copy_of(tree);
    let mut own = std::mem::replace(tree, large);
    own.reset();
    pool.lock()
        .unwrap_or_else(PoisonError::into_inner)
        .own
        .push(own);
}

/// The subject's tree goes to the pool; `tree` becomes one of a scratch's size from the pool, or a
/// new one.
#[cold]
#[inline(never)]
fn give_back(pool: &Mutex<Pool>, tree: &mut BlastIntervalTree) {
    let own = pool
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .own
        .pop();
    let own = own.unwrap_or_else(new_own);
    let used = std::mem::replace(tree, own);
    pool.lock()
        .unwrap_or_else(PoisonError::into_inner)
        .large
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

    // A subject that outgrows the scratch's tree continues in the pool's with the same answers as a
    // tree of its own, in arrays left by a longer subject of another range, and the next subject
    // starts empty in a tree of a scratch's size.
    #[test]
    fn x_large_subject_answers_as_in_the_scratchs_own_tree() {
        let pool = Mutex::new(Pool::new());
        let mut scratch_tree = new_own();
        let mut plain = BlastIntervalTree::new(0, 8000, 0, 1);

        // A subject with a grid past the scratch's room takes a pool tree at once.
        scratch_tree.reset_with_bounds(0, 8000, 0, 90_000);
        assert!(scratch_tree.x_may_outgrow(OWN_NODES, NEW_NODES_PER_HSP, OWN_ENTRIES, OWN_CELLS));
        take(&pool, &mut scratch_tree);
        for i in 0..600 {
            add(&mut scratch_tree, hsp(i * 4, i * 4 + 60, 50));
        }
        give_back(&pool, &mut scratch_tree);
        assert!(scratch_tree.x_node_capacity() <= OWN_NODES);

        // A subject with a small grid stays in the scratch's tree until its nodes near the room;
        // then it continues in the tree above, whose arrays hold the other subject's old contents.
        scratch_tree.reset_with_bounds(0, 8000, 0, 8000);
        plain.reset_with_bounds(0, 8000, 0, 8000);
        assert!(!scratch_tree.x_may_outgrow(OWN_NODES, NEW_NODES_PER_HSP, OWN_ENTRIES, OWN_CELLS));
        let mut exchanged = 0;
        for i in 0..1000 {
            if scratch_tree.x_node_capacity() <= OWN_NODES
                && scratch_tree.x_may_outgrow(OWN_NODES, NEW_NODES_PER_HSP, OWN_ENTRIES, OWN_CELLS)
            {
                take(&pool, &mut scratch_tree);
                exchanged += 1;
            }
            let value = hsp(i * 7 % 7000, i * 7 % 7000 + 40 + i % 90, 30 + i % 40);
            add(&mut scratch_tree, value);
            add(&mut plain, value);
        }
        assert_eq!(exchanged, 1);
        assert_eq!(scratch_tree.node_count(), plain.node_count());
        for i in 0..3000 {
            let probe = hsp(i * 3 % 7000, i * 3 % 7000 + 10 + i % 120, 10 + i % 90);
            assert_eq!(
                scratch_tree.contains_hsp(&probe, 0, 0),
                plain.contains_hsp(&probe, 0, 0),
                "{probe:?}"
            );
        }
        // HSPs that the small tree's first arrays did not hold answer as well.
        assert!(scratch_tree.contains_hsp(&hsp(70, 80, 10), 0, 0));

        give_back(&pool, &mut scratch_tree);
        scratch_tree.reset_with_bounds(0, 8000, 0, 8000);
        assert_eq!(scratch_tree.node_count(), 1);
        assert!(!scratch_tree.contains_hsp(&hsp(70, 80, 10), 0, 0));
    }

    // The pool keeps one tree per large subject held at the same time and hands out the largest
    // first.
    #[test]
    fn x_pool_holds_trees_in_use_at_once_largest_first() {
        let pool = Mutex::new(Pool::new());
        let mut a = new_own();
        let mut b = new_own();
        a.reset_with_bounds(0, 2000, 0, 2000);
        b.reset_with_bounds(0, 2000, 0, 2000);
        take(&pool, &mut a);
        take(&pool, &mut b);
        assert_eq!(pool.lock().unwrap().own.len(), 2);
        for i in 0..900 {
            add(&mut a, hsp(i * 2, i * 2 + 60, 50));
        }
        let large = a.x_array_bytes();
        let small = b.x_array_bytes();
        assert!(large > small);
        give_back(&pool, &mut b);
        give_back(&pool, &mut a);
        assert_eq!(pool.lock().unwrap().large.len(), 2);
        assert_eq!(pool.lock().unwrap().own.len(), 0);

        take(&pool, &mut b);
        assert_eq!(b.x_array_bytes(), large);
        give_back(&pool, &mut b);

        release_pool(&pool);
        assert_eq!(pool.lock().unwrap().large.len(), 0);
    }
}
