//! NCBI BLAST Interval Tree implementation for HSP containment checking
//!
//! NCBI reference: blast_itree.c
//!
//! This module implements the interval tree data structure used by NCBI BLAST
//! to efficiently check if an HSP is contained within existing HSPs.
//!
//! The tree is organized as follows:
//! - Primary index: query offset range
//! - Secondary index (midpoint lists): subject offset range
//! - Nodes can be internal (with children) or leaves (with HSP)

#![allow(dead_code)]

// NCBI reference: blast_itree.h:83-89
// typedef enum EITreeIndexMethod {
//     eQueryOnly,
//     eQueryAndSubject,
//     eQueryOnlyStrandIndifferent
// } EITreeIndexMethod;
/// How HSPs added to an interval tree are indexed
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IndexMethod {
    /// Index by query offset only
    QueryOnly,
    /// Index by query and then by subject offset
    QueryAndSubject,
    /// Index by query offset only, strand indifferent
    QueryOnlyStrandIndifferent,
}

// NCBI reference: blast_itree.c:41-45
// enum EIntervalDirection {
//     eIntervalTreeLeft,
//     eIntervalTreeRight,
//     eIntervalTreeNeither
// };
/// Which half of a node's range for tree navigation
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum IntervalDirection {
    /// Node will handle left half of parent node
    Left,
    /// Node will handle right half of parent node
    Right,
    /// No parent node is assumed (for root allocation)
    Neither,
}

// Result of endpoint comparison
// NCBI reference: blast_itree.c:241-301 (s_HSPsHaveCommonEndpoint returns pointer to better HSP)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum EndpointResult {
    /// Input HSP is better - keep it
    KeepInput,
    /// Tree HSP is better - keep it
    KeepTree,
}

/// HSP data stored in interval tree nodes
/// NCBI reference: blast_itree.c - tree stores pointers to BlastHSP
#[derive(Clone, Copy, Debug)]
pub struct TreeHsp {
    pub query_offset: i32,
    pub query_end: i32,
    pub subject_offset: i32,
    pub subject_end: i32,
    pub score: i32,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:539-546
    // if ( index_method == eQueryOnlyStrandIndifferent &&
    //      query_info->contexts[hsp->context].frame == -1 ) {
    //     region_end = query_start - hsp->query.offset;
    //     region_start = query_start - hsp->query.end;
    //     query_start = query_start -
    //                   query_info->contexts[hsp->context].query_length - 1;
    // }
    pub query_frame: i32,
    pub query_length: i32,
    /// Query strand offset (for context comparison)
    /// NCBI reference: blast_itree.c:819 - in_q_start != tree_q_start
    pub query_context_offset: i32,
    /// Subject frame sign (positive = forward, negative = reverse)
    /// NCBI reference: blast_itree.c:823 - SIGN(subject.frame)
    pub subject_frame_sign: i32,
}

/// Interval tree node
/// NCBI reference: blast_itree.c:48-58 SIntervalNode
#[derive(Clone, Debug)]
struct IntervalNode {
    /// Left boundary of this node's range
    leftend: i32,
    /// Right boundary of this node's range
    rightend: i32,
    /// Index of left child (0 = none)
    leftptr: i32,
    /// Index of mid child (0 = none for internal nodes, or stores query_context_offset for leaves)
    midptr: i32,
    /// Index of right child (0 = none)
    rightptr: i32,
    /// HSP stored at this node (only for leaf nodes)
    hsp: Option<TreeHsp>,
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:359-362,398-404
    /// ```c
    ///             if (best_hsp == next_node->hsp)
    ///                 return TRUE;
    ///             else if (best_hsp == in_hsp)
    ///                 list_node->midptr = tmp_index;
    /// ...
    ///                 /* leaf gets removed */
    ///                 if (target_offset < midpt)
    ///                     root_node->leftptr = 0;
    /// ```
    /// NCBI unlinks such a leaf from the tree. This flag records that fact for the grid index,
    /// which still lists the leaf; a flagged leaf is never returned as a container.
    /// EXPERIMENT (LOSAT_X_ITREEFAST): this leaf was unlinked from the tree.
    x_dead: bool,
}

impl IntervalNode {
    fn new_internal(leftend: i32, rightend: i32) -> Self {
        Self {
            leftend,
            rightend,
            leftptr: 0,
            midptr: 0,
            rightptr: 0,
            hsp: None,
            x_dead: false,
        }
    }

    fn new_leaf(hsp: TreeHsp, query_context_offset: i32) -> Self {
        Self {
            leftend: 0,
            rightend: 0,
            leftptr: query_context_offset, // NCBI stores q_start offset here for leaves
            midptr: 0,
            rightptr: 0,
            hsp: Some(hsp),
            x_dead: false,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:822-831,952-994
// ```c
//     if (in_hsp->score <= tree_hsp->score &&
//         SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame) &&
//         CONTAINED_IN_HSP(tree_hsp->query.offset, tree_hsp->query.end,
//                               in_hsp->query.offset,
//                               tree_hsp->subject.offset, tree_hsp->subject.end,
//                               in_hsp->subject.offset) &&
// ...
//     while (node->hsp == NULL) {
// ...
//     return s_HSPIsContained(hsp, query_start,
//                             node->hsp, node->leftptr,
//                             min_diag_separation);
// ```
// The side indexes below answer what `BlastIntervalTreeContainsHSP` answers (does some HSP
// in the tree satisfy `s_HSPIsContained` for the input) and what
// `s_IntervalTreeHasHSPEndpoint` answers (is there an HSP with the same end point), without
// the walk. The tree stays the master copy and the predicates are the unchanged ones. The
// argument above explains why the answers are the same; LOSAT_X_ITREESHADOW checks it.
// ---------------------------------------------------------------------------
// EXPERIMENT (LOSAT_X_ITREEFAST / LOSAT_X_ITREESHADOW): two side indexes that
// answer, without walking the tree, the two questions the tree is asked most
// often and usually answers with "no".
//
// The tree is two-level: a tree over query offsets whose nodes each carry a
// tree over subject offsets. A containment query walks one path of the first
// and, at every node of it, one path of the second, so it visits a number of
// nodes proportional to the product of the two depths before it can say that
// nothing contains the HSP. The side indexes are derived from the same HSPs:
//
//   * A grid over (query offset, subject offset). Every HSP of the tree is
//     listed in each grid cell its bounding box overlaps. `s_HSPIsContained`
//     requires the start point of the input HSP to lie inside the box of the
//     tree HSP, so every tree HSP that can contain the input is listed in the
//     cell of the input's start point, and testing those with the unchanged
//     predicate decides whether a containing HSP exists.
//   * The set of (query, subject) start points and the set of end points of
//     the HSPs ever added. `s_IntervalTreeHasHSPEndpoint` only acts on tree
//     HSPs with the same start (or end) point as the input; when the set has
//     no such point the walk finds nothing, removes nothing and returns
//     FALSE, so it is skipped. The set is a bit table that can report a point
//     it does not hold, and points are never deleted from it; either at worst
//     makes the unchanged walk run when it was not needed.
//
// Why the grid answer is the tree's answer: the tree walk reaches every HSP
// that contains the input (an HSP is stored at the first node whose middle it
// straddles, or as the only leaf of a half it lies in; a containing HSP
// straddles every middle the input straddles and lies on the input's side of
// every other one, in query and in subject offsets alike), and it reaches no
// HSP that was unlinked. So both say whether some linked HSP satisfies
// `s_HSPIsContained`. Which containing HSP is returned can differ when there
// are several; callers only test for existence, except trace output, and the
// index is off whenever tracing is on.
// ---------------------------------------------------------------------------

#[cfg(test)]
thread_local! {
    /// Lets a test choose the mode of the trees it creates.
    static X_ITREE_TEST_MODE: std::cell::Cell<Option<u8>> = const { std::cell::Cell::new(None) };
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:931-935,952
/// ```c
/// BlastIntervalTreeContainsHSP(const BlastIntervalTree *tree,
///                              const BlastHSP *hsp,
///                              const BlastQueryInfo *query_info,
///                              Int4 min_diag_separation)
/// ...
///     while (node->hsp == NULL) {
/// ```
/// Reads the LOSAT_X_ITREEFAST / LOSAT_X_ITREESHADOW switches once. The mode chooses
/// whether the walk above (0), the side indexes (1), or both (2) answer a containment
/// query. It is 0 whenever tracing is on, because the index may return a different
/// container than the walk when there are several.
/// 0 = tree only, 1 = side indexes, 2 = both and compare.
fn x_itree_mode() -> u8 {
    #[cfg(test)]
    if let Some(mode) = X_ITREE_TEST_MODE.with(|mode| mode.get()) {
        return mode;
    }
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if super::tracing::enabled() {
            0
        } else if std::env::var_os("LOSAT_X_ITREESHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_ITREEFAST").is_some() {
            1
        } else {
            0
        }
    })
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:822-831
/// ```c
///     if (in_hsp->score <= tree_hsp->score &&
///         SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame) &&
///         CONTAINED_IN_HSP(tree_hsp->query.offset, tree_hsp->query.end,
///                               in_hsp->query.offset,
/// ```
/// A cell lists the tree leaves whose box (the `tree_hsp` offsets and ends in this test)
/// overlaps the cell. No value is computed here.
/// One listing of a leaf in a grid cell.
#[derive(Clone, Copy)]
struct XCellEntry {
    leaf: u32,
    /// 1 + index of the next entry of the same cell, 0 = none.
    next: u32,
}

// No NCBI counterpart: size limits of the grid index (memory and fallback thresholds);
// they do not change any value NCBI computes. Past a limit the index is dropped and the
// tree walk answers.
/// Grid cells at most (the per-cell heads take four bytes each).
const X_MAX_CELLS: usize = 1 << 16;
/// An HSP whose box overlaps more cells than this is kept in `overflow`.
const X_MAX_CELLS_PER_HSP: usize = 1024;
/// With more HSPs than this in `overflow` the index is abandoned.
const X_MAX_OVERFLOW: usize = 512;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:440-443,452-461
/// ```c
///     if (which_end == eIntervalTreeLeft)
///         target_offset = in_q_start + in_hsp->query.offset;
///     else
///         target_offset = in_q_start + in_hsp->query.end;
/// ...
///         tmp_index = root_node->midptr;
/// ```
/// `s_IntervalTreeHasHSPEndpoint` only acts on tree HSPs with the same start (or end)
/// point as the input. This set holds the start and end points of the HSPs added, as a
/// hash bit table. "May contain" is false only when no such HSP was added; a false
/// "true" only makes the unchanged walk run.
/// A set of (query, subject) points that may report points it does not hold
/// but never misses one it holds: a bit table with two bits per point.
struct XPointSet {
    bits: Vec<u64>,
    points: usize,
}

impl XPointSet {
    fn new() -> Self {
        Self {
            bits: Vec::new(),
            points: 0,
        }
    }

    fn clear(&mut self) {
        self.bits.clear();
        self.points = 0;
    }

    /// Two bit positions of a point; `salt` separates start from end points.
    #[inline]
    fn positions(&self, q: i32, s: i32, salt: u64) -> (usize, usize) {
        // splitmix64 finalizer
        let mut z = (((q as u32 as u64) << 32) | s as u32 as u64) ^ salt;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        let mask = self.bits.len() * 64 - 1;
        (z as usize & mask, (z >> 32) as usize & mask)
    }

    /// Bits must be at least 32 per point; true when the table has to grow
    /// (the caller then re-inserts every point).
    #[inline]
    fn needs_growth(&self, more: usize) -> bool {
        (self.points + more) * 32 > self.bits.len() * 64
    }

    fn grow(&mut self) {
        let words = (self.bits.len() * 2).max(128);
        self.bits = vec![0u64; words];
        self.points = 0;
    }

    #[inline]
    fn insert(&mut self, q: i32, s: i32, salt: u64) {
        let (a, b) = self.positions(q, s, salt);
        self.bits[a >> 6] |= 1u64 << (a & 63);
        self.bits[b >> 6] |= 1u64 << (b & 63);
        self.points += 1;
    }

    #[inline]
    fn may_contain(&self, q: i32, s: i32, salt: u64) -> bool {
        if self.bits.is_empty() {
            return false;
        }
        let (a, b) = self.positions(q, s, salt);
        (self.bits[a >> 6] >> (a & 63)) & (self.bits[b >> 6] >> (b & 63)) & 1 != 0
    }
}

// No NCBI counterpart: hash salts that keep start points and end points apart in
// `XPointSet`; they do not change any value NCBI computes.
const X_SALT_START: u64 = 0;
const X_SALT_END: u64 = 0x9E37_79B9_7F4A_7C15;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:164-170,195-201
// ```c
//     tree->num_alloc = size;
//     tree->num_used = 0;
//     tree->s_min = s_start;
//     tree->s_max = s_end;
//
//     /* The first structure in tree->nodes is the root */
//     s_IntervalRootNodeInit(tree, q_start, q_end, &retval);
// ...
//     tree->num_used = 1;
// ```
// The state of the side indexes, tied to the same query and subject ranges as the tree
// (`q_min`, `s_min`) and emptied whenever the tree is emptied. It is derived from the
// HSPs in the tree and holds nothing else.
struct XTreeIndex {
    mode: u8,
    /// False once the tree holds an HSP this index does not describe.
    usable: bool,
    shift: u32,
    q_min: i32,
    s_min: i32,
    nq: usize,
    ns: usize,
    /// Per cell: 1 + index of its first entry, 0 = empty. Allocated on the
    /// first HSP.
    head: Vec<u32>,
    /// The cells with a non-zero head, for `clear`.
    used_cells: Vec<u32>,
    entries: Vec<XCellEntry>,
    overflow: Vec<u32>,
    ends: XPointSet,
    /// Leaves unlinked from the tree so far (for the shadow check).
    removed: usize,
}

impl XTreeIndex {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:146-170
    /// ```c
    /// Blast_IntervalTreeInit(Int4 q_start, Int4 q_end,
    ///                        Int4 s_start, Int4 s_end)
    /// ...
    ///     tree->s_min = s_start;
    ///     tree->s_max = s_end;
    /// ...
    ///     s_IntervalRootNodeInit(tree, q_start, q_end, &retval);
    /// ```
    /// Makes an empty index for the ranges the tree is made for.
    fn new(q_min: i32, q_max: i32, s_min: i32, s_max: i32) -> Self {
        let mode = x_itree_mode();
        let mut index = Self {
            mode,
            usable: mode != 0,
            shift: 0,
            q_min,
            s_min,
            nq: 1,
            ns: 1,
            head: Vec::new(),
            used_cells: Vec::new(),
            entries: Vec::new(),
            overflow: Vec::new(),
            ends: XPointSet::new(),
            removed: 0,
        };
        index.set_bounds(q_min, q_max, s_min, s_max);
        index
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:166-170
    /// ```c
    ///     tree->s_min = s_start;
    ///     tree->s_max = s_end;
    ///
    ///     /* The first structure in tree->nodes is the root */
    ///     s_IntervalRootNodeInit(tree, q_start, q_end, &retval);
    /// ```
    /// The index follows `reset_with_bounds`: new ranges, empty index; the grid is sized
    /// from the ranges.
    fn set_bounds(&mut self, q_min: i32, q_max: i32, s_min: i32, s_max: i32) {
        let q_span = (q_max as i64 - q_min as i64).max(0) as usize;
        let s_span = (s_max as i64 - s_min as i64).max(0) as usize;
        let mut shift = 6u32;
        while ((q_span >> shift) + 1).saturating_mul((s_span >> shift) + 1) > X_MAX_CELLS {
            shift += 1;
        }
        self.shift = shift;
        self.q_min = q_min;
        self.s_min = s_min;
        self.nq = (q_span >> shift) + 1;
        self.ns = (s_span >> shift) + 1;
        self.head = Vec::new();
        self.used_cells.clear();
        self.clear();
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:195-201
    /// ```c
    ///     SIntervalNode *root = tree->nodes;
    ///
    ///     tree->num_used = 1;
    ///     root->leftptr = 0;
    ///     root->midptr = 0;
    ///     root->rightptr = 0;
    /// ```
    /// The index is emptied exactly when the tree is reset.
    fn clear(&mut self) {
        for &cell in &self.used_cells {
            self.head[cell as usize] = 0;
        }
        self.used_cells.clear();
        self.entries.clear();
        self.overflow.clear();
        self.ends.clear();
        self.removed = 0;
        self.usable = self.mode != 0;
    }

    /// No NCBI counterpart: maps an offset to a grid column; it does not change any value
    /// NCBI computes.
    /// Grid column of an absolute query offset (monotone, clamped).
    #[inline]
    fn q_cell(&self, q: i32) -> usize {
        let rel = (q as i64 - self.q_min as i64).max(0) as usize >> self.shift;
        rel.min(self.nq - 1)
    }

    /// No NCBI counterpart: maps an offset to a grid row; it does not change any value NCBI
    /// computes.
    /// Grid row of a subject offset (monotone, clamped).
    #[inline]
    fn s_cell(&self, s: i32) -> usize {
        let rel = (s as i64 - self.s_min as i64).max(0) as usize >> self.shift;
        rel.min(self.ns - 1)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:597-599,822-831
    /// ```c
    ///     nodes[new_index].leftptr = query_start;
    ///     nodes[new_index].midptr = 0;
    ///     nodes[new_index].hsp = hsp;
    /// ...
    ///         CONTAINED_IN_HSP(tree_hsp->query.offset, tree_hsp->query.end,
    /// ```
    /// Called when NCBI creates the leaf for an HSP. The leaf is listed in each cell of its
    /// box (query offsets shifted by the strand start `query_start`, as `leftptr` holds it).
    /// List `leaf` in every cell its box overlaps.
    fn register(&mut self, leaf: usize, hsp: &TreeHsp, query_start: i32) {
        let q_a = self.q_cell(query_start + hsp.query_offset.min(hsp.query_end));
        let q_b = self.q_cell(query_start + hsp.query_offset.max(hsp.query_end));
        let s_a = self.s_cell(hsp.subject_offset.min(hsp.subject_end));
        let s_b = self.s_cell(hsp.subject_offset.max(hsp.subject_end));
        let cells = (q_b - q_a + 1) * (s_b - s_a + 1);
        if cells > X_MAX_CELLS_PER_HSP {
            self.overflow.push(leaf as u32);
            if self.overflow.len() > X_MAX_OVERFLOW {
                self.usable = false;
            }
            return;
        }
        if self.head.is_empty() {
            self.head = vec![0u32; self.nq * self.ns];
        }
        for q in q_a..=q_b {
            for s in s_a..=s_b {
                let cell = q * self.ns + s;
                if self.head[cell] == 0 {
                    self.used_cells.push(cell as u32);
                }
                self.entries.push(XCellEntry {
                    leaf: leaf as u32,
                    next: self.head[cell],
                });
                self.head[cell] = self.entries.len() as u32;
            }
        }
    }
}

/// Interval tree for HSP containment checking
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.h:73-80
/// ```c
/// typedef struct BlastIntervalTree {
///     SIntervalNode *nodes;
///     Int4 num_alloc;
///     Int4 num_used;
///     Int4 s_min;
///     Int4 s_max;
/// } BlastIntervalTree;
/// ```
pub struct BlastIntervalTree {
    /// Array of nodes (index 0 is root for query-indexed tree)
    nodes: Vec<IntervalNode>,
    /// Minimum subject offset
    s_min: i32,
    /// Maximum subject offset
    s_max: i32,
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:158-167
    /// ```c
    ///     tree->nodes = (SIntervalNode *)malloc(size * sizeof(SIntervalNode));
    /// ...
    ///     tree->num_alloc = size;
    ///     tree->num_used = 0;
    ///     tree->s_min = s_start;
    ///     tree->s_max = s_end;
    /// ```
    /// The tree owns its side indexes; they are built from, and follow, the node array above.
    /// EXPERIMENT (LOSAT_X_ITREEFAST)
    x: XTreeIndex,
}

impl BlastIntervalTree {
    /// Initialize a new interval tree
    /// NCBI reference: blast_itree.c:96-166 Blast_IntervalTreeInit
    ///
    /// # Arguments
    /// * `q_min` - Minimum query offset (usually 0)
    /// * `q_max` - Maximum query offset (query_length + 1)
    /// * `s_min` - Minimum subject offset (usually 0)
    /// * `s_max` - Maximum subject offset (subject_length + 1)
    pub fn new(q_min: i32, q_max: i32, s_min: i32, s_max: i32) -> Self {
        let mut tree = Self {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:135-153
            // ```c
            // Int4 size = 100;
            // tree->nodes = (SIntervalNode *)malloc(size * sizeof(SIntervalNode));
            // tree->num_alloc = size;
            // tree->num_used = 0;
            // ```
            nodes: Vec::with_capacity(100),
            s_min,
            s_max,
            x: XTreeIndex::new(q_min, q_max, s_min, s_max),
        };

        // Create root node for query range
        let root = IntervalNode::new_internal(q_min, q_max);
        tree.nodes.push(root);

        tree
    }

    /// Reserve interval-tree node storage before a dense traceback/prune pass.
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:66-70
    /// ```c
    /// if (tree->num_used == tree->num_alloc) {
    ///     tree->num_alloc = 2 * tree->num_alloc;
    ///     tree->nodes = (SIntervalNode *)realloc(tree->nodes, tree->num_alloc *
    ///                                                  sizeof(SIntervalNode));
    /// }
    /// ```
    /// Rust keeps the same node indices and insertion order, but can reserve
    /// enough Vec storage from the observed HSP count to avoid repeated reallocs.
    pub fn reserve_nodes_for_hsps(&mut self, hsp_count: usize) {
        let target = hsp_count.saturating_mul(3).saturating_add(1).max(100);
        if target > self.nodes.capacity() {
            self.nodes.reserve(target - self.nodes.capacity());
        }
    }

    /// Reset the tree for reuse
    /// NCBI reference: blast_itree.c:168-193 Blast_IntervalTreeReset
    pub fn reset(&mut self) {
        // Keep only the root node, reset its children
        if !self.nodes.is_empty() {
            let leftend = self.nodes[0].leftend;
            let rightend = self.nodes[0].rightend;
            self.nodes.clear();
            self.nodes
                .push(IntervalNode::new_internal(leftend, rightend));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:195-201
        // ```c
        //     tree->num_used = 1;
        //     root->leftptr = 0;
        //     root->midptr = 0;
        //     root->rightptr = 0;
        // ```
        // Reset point: the side indexes are emptied with the tree.
        self.x.clear();
    }

    /// Reset the tree bounds for reuse.
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:146-170
    /// ```c
    /// Blast_IntervalTreeInit(Int4 q_start, Int4 q_end,
    ///                        Int4 s_start, Int4 s_end)
    /// ...
    /// tree->s_min = s_start;
    /// tree->s_max = s_end;
    /// s_IntervalRootNodeInit(tree, q_start, q_end, &retval);
    /// ```
    pub fn reset_with_bounds(&mut self, q_min: i32, q_max: i32, s_min: i32, s_max: i32) {
        self.s_min = s_min;
        self.s_max = s_max;
        self.nodes.clear();
        self.nodes.push(IntervalNode::new_internal(q_min, q_max));
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:166-170
        // ```c
        //     tree->s_min = s_start;
        //     tree->s_max = s_end;
        //
        //     /* The first structure in tree->nodes is the root */
        //     s_IntervalRootNodeInit(tree, q_start, q_end, &retval);
        // ```
        // Re-initialisation point: the side indexes are rebuilt for the new ranges.
        self.x.set_bounds(q_min, q_max, s_min, s_max);
    }

    /// Allocate a new internal node
    /// NCBI reference: blast_itree.c:57-112 s_IntervalNodeInit
    ///
    /// CRITICAL: Right child uses midpt + 1, not midpt!
    /// The two subregions do not overlap, may be of length one,
    /// and must completely cover the parent region.
    fn alloc_internal_node(&mut self, parent_idx: usize, which_half: IntervalDirection) -> usize {
        // NCBI reference: blast_itree.c:89-95
        // parent_node = tree->nodes + parent_index;
        // new_node = tree->nodes + new_index;
        // midpt = ((Int8) parent_node->leftend + (Int8) parent_node->rightend) / 2;
        let parent = &self.nodes[parent_idx];
        let midpt = ((parent.leftend as i64 + parent.rightend as i64) / 2) as i32;

        // NCBI reference: blast_itree.c:102-109
        // if (dir == eIntervalTreeLeft) {
        //     new_node->leftend = parent_node->leftend;
        //     new_node->rightend = (Int4) midpt;
        // } else {
        //     new_node->leftend = (Int4) midpt + 1;  // <-- +1 is CRITICAL
        //     new_node->rightend = parent_node->rightend;
        // }
        let (leftend, rightend) = match which_half {
            IntervalDirection::Left => (parent.leftend, midpt),
            IntervalDirection::Right => (midpt + 1, parent.rightend), // +1 is CRITICAL
            IntervalDirection::Neither => {
                panic!("Neither direction not valid for child allocation")
            }
        };

        let node = IntervalNode::new_internal(leftend, rightend);
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// Allocate a new root node for subject-indexed subtree
    /// NCBI reference: blast_itree.c:269-299 s_IntervalRootNodeInit
    fn alloc_subject_root_node(&mut self) -> usize {
        let node = IntervalNode::new_internal(self.s_min, self.s_max);
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// Allocate a new leaf node
    fn alloc_leaf_node(&mut self, hsp: TreeHsp, query_context_offset: i32) -> usize {
        let node = IntervalNode::new_leaf(hsp, query_context_offset);
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// Determine whether an input HSP shares a common start- or endpoint with
    /// an HSP from an interval tree.
    ///
    /// NCBI reference: blast_itree.c:241-301 s_HSPsHaveCommonEndpoint
    ///
    /// Returns: None if no common endpoint, Some(EndpointResult) indicating which HSP to keep
    fn hsps_have_common_endpoint(
        in_hsp: &TreeHsp,
        in_q_start: i32,
        tree_hsp: &TreeHsp,
        tree_q_start: i32,
        which_end: IntervalDirection,
    ) -> Option<EndpointResult> {
        // NCBI reference: blast_itree.c:250-254
        // check if alignments are from different query sequences or query strands
        if in_q_start != tree_q_start {
            return None;
        }

        // NCBI reference: blast_itree.c:256-259
        // check if alignments are from different subject strands
        if in_hsp.subject_frame_sign.signum() != tree_hsp.subject_frame_sign.signum() {
            return None;
        }

        // NCBI reference: blast_itree.c:261-268
        let match_found = match which_end {
            IntervalDirection::Left => {
                // if (which_end == eIntervalTreeLeft) {
                //     match = in_hsp->query.offset == tree_hsp->query.offset &&
                //             in_hsp->subject.offset == tree_hsp->subject.offset;
                // }
                in_hsp.query_offset == tree_hsp.query_offset
                    && in_hsp.subject_offset == tree_hsp.subject_offset
            }
            IntervalDirection::Right => {
                // else {
                //     match = in_hsp->query.end == tree_hsp->query.end &&
                //             in_hsp->subject.end == tree_hsp->subject.end;
                // }
                in_hsp.query_end == tree_hsp.query_end && in_hsp.subject_end == tree_hsp.subject_end
            }
            IntervalDirection::Neither => false,
        };

        if match_found {
            // NCBI reference: blast_itree.c:273-278
            // keep the higher scoring HSP
            if in_hsp.score > tree_hsp.score {
                return Some(EndpointResult::KeepInput);
            }
            if in_hsp.score < tree_hsp.score {
                return Some(EndpointResult::KeepTree);
            }

            // NCBI reference: blast_itree.c:280-286
            // for equal scores, pick the shorter HSP
            let in_q_length = in_hsp.query_end - in_hsp.query_offset;
            let tree_q_length = tree_hsp.query_end - tree_hsp.query_offset;
            if in_q_length > tree_q_length {
                return Some(EndpointResult::KeepTree);
            }
            if in_q_length < tree_q_length {
                return Some(EndpointResult::KeepInput);
            }

            // NCBI reference: blast_itree.c:288-293
            let in_s_length = in_hsp.subject_end - in_hsp.subject_offset;
            let tree_s_length = tree_hsp.subject_end - tree_hsp.subject_offset;
            if in_s_length > tree_s_length {
                return Some(EndpointResult::KeepTree);
            }
            if in_s_length < tree_s_length {
                return Some(EndpointResult::KeepInput);
            }

            // NCBI reference: blast_itree.c:295-297
            // HSPs are identical; favor the one already in the tree
            return Some(EndpointResult::KeepTree);
        }

        None
    }

    /// Determine whether a subtree of an interval tree contains an HSP that
    /// shares a common endpoint with the input HSP. The subtree indexes subject
    /// offsets, and represents the midpoint list of a tree node that indexes
    /// query offsets.
    ///
    /// NCBI reference: blast_itree.c:318-411 s_MidpointTreeHasHSPEndpoint
    ///
    /// Returns TRUE if the HSP should not be added to the tree because it shares
    /// an existing endpoint with a 'better' HSP already there.
    ///
    /// Side effect: Removes worse HSPs from the tree.
    fn midpoint_tree_has_hsp_endpoint(
        &mut self,
        root_index: usize,
        in_hsp: &TreeHsp,
        in_q_start: i32,
        which_end: IntervalDirection,
    ) -> bool {
        // NCBI reference: blast_itree.c:331-334
        let target_offset = match which_end {
            IntervalDirection::Left => in_hsp.subject_offset,
            IntervalDirection::Right => in_hsp.subject_end,
            IntervalDirection::Neither => return false,
        };

        let mut root_idx = root_index;

        // NCBI reference: blast_itree.c:338 - while (1)
        loop {
            // NCBI reference: blast_itree.c:343-366
            // First perform matching endpoint tests on all of the HSPs in the
            // midpoint list for the current node.
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:350-365
            // tmp_index = root_node->midptr;
            // list_node = root_node;
            // next_node = tree->nodes + tmp_index;
            // while (tmp_index != 0) {
            //     ...
            //     tmp_index = next_node->midptr;
            //     if (best_hsp == next_node->hsp)
            //         return TRUE;
            //     else if (best_hsp == in_hsp)
            //         list_node->midptr = tmp_index;
            //     list_node = next_node;
            //     next_node = tree->nodes + tmp_index;
            // }
            let mut list_idx = root_idx;
            let mut tmp_index = self.nodes[root_idx].midptr;

            while tmp_index != 0 {
                let (next_idx, result) = {
                    let tmp_node = &self.nodes[tmp_index as usize];
                    let tree_hsp = match tmp_node.hsp.as_ref() {
                        Some(hsp) => hsp,
                        None => break,
                    };
                    let tree_q_start = tmp_node.leftptr; // leftptr stores query_context_offset for leaves
                    let next_idx = tmp_node.midptr;
                    let result = Self::hsps_have_common_endpoint(
                        in_hsp,
                        in_q_start,
                        tree_hsp,
                        tree_q_start,
                        which_end,
                    );
                    (next_idx, result)
                };

                match result {
                    Some(EndpointResult::KeepTree) => return true,
                    Some(EndpointResult::KeepInput) => {
                        // Remove worse HSP from list: list_node->midptr = tmp_index
                        self.nodes[list_idx].midptr = next_idx;
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:361-362
                        // ```c
                        //             else if (best_hsp == in_hsp)
                        //                 list_node->midptr = tmp_index;
                        // ```
                        // The list unlink is the line above; the flag and the count only tell the grid index
                        // that the leaf is gone.
                        self.nodes[tmp_index as usize].x_dead = true;
                        self.x.removed += 1;
                    }
                    None => {}
                }
                list_idx = tmp_index as usize;
                tmp_index = next_idx;
            }

            // NCBI reference: blast_itree.c:368-376
            // Descend to the left or right subtree, whichever one contains the
            // endpoint from in_hsp
            let node = &self.nodes[root_idx];
            let midpt = ((node.leftend as i64 + node.rightend as i64) / 2) as i32;

            let next_child_idx = if target_offset < midpt {
                node.leftptr
            } else if target_offset > midpt {
                node.rightptr
            } else {
                0
            };

            // NCBI reference: blast_itree.c:382-383
            if next_child_idx == 0 {
                return false;
            }

            // NCBI reference: blast_itree.c:385-406
            let next_node = &self.nodes[next_child_idx as usize];
            if next_node.hsp.is_some() {
                // Reached a leaf; compare in_hsp with the alignment in the leaf
                let tree_hsp = next_node.hsp.as_ref().unwrap();
                let tree_q_start = next_node.leftptr;

                let result = Self::hsps_have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    tree_hsp,
                    tree_q_start,
                    which_end,
                );

                match result {
                    Some(EndpointResult::KeepTree) => return true,
                    Some(EndpointResult::KeepInput) => {
                        // NCBI reference: blast_itree.c:398-404
                        // leaf gets removed
                        if target_offset < midpt {
                            self.nodes[root_idx].leftptr = 0;
                        } else {
                            self.nodes[root_idx].rightptr = 0;
                        }
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:398-404
                        // ```c
                        //             else if (best_hsp == in_hsp) {
                        //                 /* leaf gets removed */
                        //                 if (target_offset < midpt)
                        //                     root_node->leftptr = 0;
                        // ```
                        // The leaf unlink is the code above; the flag and the count only tell the grid index
                        // that the leaf is gone.
                        self.nodes[next_child_idx as usize].x_dead = true;
                        self.x.removed += 1;
                        return false;
                    }
                    None => {}
                }
                return false;
            }

            // NCBI reference: blast_itree.c:408
            root_idx = next_child_idx as usize;
        }
    }

    /// Determine whether an interval tree contains one or more HSPs that share
    /// a common endpoint with the input HSP. Remove from the tree all such HSPs
    /// that are "worse" than the input.
    ///
    /// NCBI reference: blast_itree.c:428-506 s_IntervalTreeHasHSPEndpoint
    ///
    /// Returns TRUE if the HSP should not be added to the tree because it shares
    /// an existing endpoint with a 'better' HSP already there.
    fn interval_tree_has_hsp_endpoint(
        &mut self,
        in_hsp: &TreeHsp,
        in_q_start: i32,
        which_end: IntervalDirection,
    ) -> bool {
        // NCBI reference: blast_itree.c:440-443
        let target_offset = match which_end {
            IntervalDirection::Left => in_q_start + in_hsp.query_offset,
            IntervalDirection::Right => in_q_start + in_hsp.query_end,
            IntervalDirection::Neither => return false,
        };

        let mut root_idx = 0usize;

        // NCBI reference: blast_itree.c:447 - while (1)
        loop {
            // NCBI reference: blast_itree.c:455-461
            // First perform matching endpoint tests on all of the HSPs in the
            // midpoint tree for the current node
            let midptr = self.nodes[root_idx].midptr;
            if midptr != 0 {
                if self.midpoint_tree_has_hsp_endpoint(
                    midptr as usize,
                    in_hsp,
                    in_q_start,
                    which_end,
                ) {
                    return true;
                }
            }

            // NCBI reference: blast_itree.c:466-471
            // Descend to the left or right subtree
            let node = &self.nodes[root_idx];
            let midpt = ((node.leftend as i64 + node.rightend as i64) / 2) as i32;

            let next_child_idx = if target_offset < midpt {
                node.leftptr
            } else if target_offset > midpt {
                node.rightptr
            } else {
                0
            };

            // NCBI reference: blast_itree.c:477-478
            if next_child_idx == 0 {
                return false;
            }

            // NCBI reference: blast_itree.c:480-501
            let next_node = &self.nodes[next_child_idx as usize];
            if next_node.hsp.is_some() {
                // Reached a leaf
                let tree_hsp = next_node.hsp.as_ref().unwrap();
                let tree_q_start = next_node.leftptr;

                let result = Self::hsps_have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    tree_hsp,
                    tree_q_start,
                    which_end,
                );

                match result {
                    Some(EndpointResult::KeepTree) => return true,
                    Some(EndpointResult::KeepInput) => {
                        // NCBI reference: blast_itree.c:493-499
                        // leaf gets removed
                        if target_offset < midpt {
                            self.nodes[root_idx].leftptr = 0;
                        } else {
                            self.nodes[root_idx].rightptr = 0;
                        }
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:493-499
                        // ```c
                        //             else if (best_hsp == in_hsp) {
                        //                 /* leaf gets removed */
                        //                 if (target_offset < midpt)
                        //                     root_node->leftptr = 0;
                        // ```
                        // The leaf unlink is the code above; the flag and the count only tell the grid index
                        // that the leaf is gone.
                        self.nodes[next_child_idx as usize].x_dead = true;
                        self.x.removed += 1;
                        return false;
                    }
                    None => {}
                }
                return false;
            }

            // NCBI reference: blast_itree.c:503
            root_idx = next_child_idx as usize;
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:440-443,477-478
    /// ```c
    ///     if (which_end == eIntervalTreeLeft)
    ///         target_offset = in_q_start + in_hsp->query.offset;
    ///     else
    ///         target_offset = in_q_start + in_hsp->query.end;
    /// ...
    ///         if (tmp_index == 0)
    ///             return FALSE;
    /// ```
    /// `s_IntervalTreeHasHSPEndpoint` only acts on tree HSPs with the same start (or end)
    /// point as the input. If the point set has no such point the walk finds nothing, unlinks
    /// nothing and returns FALSE, so it is skipped. Otherwise the unchanged walk runs.
    /// EXPERIMENT (LOSAT_X_ITREEFAST): `interval_tree_has_hsp_endpoint`,
    /// skipped when no HSP ever added has the input's start (or end) point.
    fn x_has_hsp_endpoint(
        &mut self,
        in_hsp: &TreeHsp,
        in_q_start: i32,
        which_end: IntervalDirection,
        x_on: bool,
    ) -> bool {
        if !x_on {
            return self.interval_tree_has_hsp_endpoint(in_hsp, in_q_start, which_end);
        }
        let possible = match which_end {
            IntervalDirection::Left => self.x.ends.may_contain(
                in_q_start + in_hsp.query_offset,
                in_hsp.subject_offset,
                X_SALT_START,
            ),
            IntervalDirection::Right => self.x.ends.may_contain(
                in_q_start + in_hsp.query_end,
                in_hsp.subject_end,
                X_SALT_END,
            ),
            IntervalDirection::Neither => return false,
        };
        if possible {
            return self.interval_tree_has_hsp_endpoint(in_hsp, in_q_start, which_end);
        }
        if self.x.mode == 2 {
            // LOSAT_X_ITREESHADOW: the walk must find nothing and unlink nothing.
            let removed = self.x.removed;
            let found = self.interval_tree_has_hsp_endpoint(in_hsp, in_q_start, which_end);
            assert!(
                !found && self.x.removed == removed,
                "LOSAT_X_ITREESHADOW: the tree found a common endpoint the point set does not have"
            );
        }
        false
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:816-831,889-899
    /// ```c
    ///     if (in_q_start != tree_q_start)
    ///         return FALSE;
    ///
    ///     if (in_hsp->score <= tree_hsp->score &&
    ///         SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame) &&
    /// ...
    ///             if (s_HSPIsContained(in_hsp, in_q_start,
    ///                                  tmp_node->hsp, tmp_node->leftptr,
    ///                                  min_diag_separation)) {
    /// ```
    /// Tests the leaves listed in the grid cell of the input's start point with the same
    /// predicate (`is_hsp_contained`). Any tree HSP that contains the input contains its
    /// start point, so it is listed in that cell. Which container is returned can differ from
    /// the walk; whether there is one does not.
    /// EXPERIMENT (LOSAT_X_ITREEFAST): a linked tree HSP that contains `hsp`,
    /// from the grid cell of the start point of `hsp`.
    fn x_find_container(
        &self,
        hsp: &TreeHsp,
        query_start: i32,
        min_diag_separation: i32,
    ) -> Option<TreeHsp> {
        let x = &self.x;
        let test = |leaf: u32| -> Option<TreeHsp> {
            let node = &self.nodes[leaf as usize];
            if node.x_dead {
                return None;
            }
            let tree_hsp = node.hsp.as_ref()?;
            Self::is_hsp_contained(
                hsp,
                query_start,
                tree_hsp,
                node.leftptr,
                min_diag_separation,
            )
            .then_some(*tree_hsp)
        };
        if !x.head.is_empty() {
            let cell =
                x.q_cell(query_start + hsp.query_offset) * x.ns + x.s_cell(hsp.subject_offset);
            let mut at = x.head[cell];
            while at != 0 {
                let entry = x.entries[at as usize - 1];
                if let Some(found) = test(entry.leaf) {
                    return Some(found);
                }
                at = entry.next;
            }
        }
        for &leaf in &x.overflow {
            if let Some(found) = test(leaf) {
                return Some(found);
            }
        }
        None
    }

    /// Add an HSP to the tree
    /// NCBI reference: blast_itree.c:510-795 BlastIntervalTreeAddHSP
    ///
    /// This is the main entry point for adding HSPs to the interval tree.
    /// For eQueryAndSubject mode, it first checks for common endpoints and
    /// removes worse HSPs before adding.
    pub fn add_hsp(&mut self, hsp: TreeHsp, query_context_offset: i32, index_method: IndexMethod) {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:537-550
        // query_start = s_GetQueryStrandOffset(query_info, hsp->context);
        // if ( index_method == eQueryOnlyStrandIndifferent &&
        //      query_info->contexts[hsp->context].frame == -1 ) {
        //     region_end = query_start - hsp->query.offset;
        //     region_start = query_start - hsp->query.end;
        //     query_start = query_start -
        //                   query_info->contexts[hsp->context].query_length - 1;
        // } else {
        //     region_start = query_start + hsp->query.offset;
        //     region_end = query_start + hsp->query.end;
        // }
        let mut query_start = query_context_offset;
        let (region_start, region_end) =
            if index_method == IndexMethod::QueryOnlyStrandIndifferent && hsp.query_frame < 0 {
                let region_end = query_start - hsp.query_offset;
                let region_start = query_start - hsp.query_end;
                query_start = query_start - hsp.query_length - 1;
                (region_start, region_end)
            } else {
                (query_start + hsp.query_offset, query_start + hsp.query_end)
            };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:558-561
        // ```c
        //     if (index_method == eQueryAndSubject) {
        //
        //         ASSERT(hsp->subject.offset >= tree->s_min);
        //         ASSERT(hsp->subject.end <= tree->s_max);
        // ```
        // NCBI looks for common end points only in this mode. The side indexes are built for it
        // and are dropped (`usable = false`) for any other mode.
        // EXPERIMENT (LOSAT_X_ITREEFAST): the side indexes describe trees
        // built with eQueryAndSubject only.
        if index_method != IndexMethod::QueryAndSubject {
            self.x.usable = false;
        }
        let x_on = self.x.usable;

        // NCBI reference: blast_itree.c:558-585
        // For eQueryAndSubject, check for common endpoints before adding
        if index_method == IndexMethod::QueryAndSubject {
            // Check left endpoint
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:577-584
            // ```c
            //         if (s_IntervalTreeHasHSPEndpoint(tree, hsp, query_start,
            //                                          eIntervalTreeLeft)) {
            //             return retval;
            //         }
            //         if (s_IntervalTreeHasHSPEndpoint(tree, hsp, query_start,
            //                                          eIntervalTreeRight)) {
            //             return retval;
            // ```
            // Dispatch point of LOSAT_X_ITREEFAST (both calls below). `x_has_hsp_endpoint` runs the
            // ported `interval_tree_has_hsp_endpoint` unless the point set shows it cannot find
            // anything.
            if self.x_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Left, x_on) {
                return; // Better HSP with same endpoint already exists
            }
            // Check right endpoint
            if self.x_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Right, x_on) {
                return; // Better HSP with same endpoint already exists
            }
        }

        // NCBI reference: blast_itree.c:591-599
        // Encapsulate the input HSP in an SIntervalNode (leaf node)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:591-599
        // ```c
        //     /* encapsulate the input HSP in an SIntervalNode */
        //     root_index = 0;
        //     new_index = s_IntervalNodeInit(tree, 0, eIntervalTreeNeither, &retval);
        // ...
        //     nodes[new_index].leftptr = query_start;
        //     nodes[new_index].midptr = 0;
        //     nodes[new_index].hsp = hsp;
        // ```
        // After NCBI creates the leaf, the side indexes list it in the grid cells of its box and
        // add its start and end points to the point set. The tree insertion that follows is
        // unchanged.
        let new_index = self.alloc_leaf_node(hsp, query_start);
        if x_on {
            self.x.register(new_index, &hsp, query_start);
            if self.x.ends.needs_growth(2) {
                // a larger table, refilled from every leaf ever allocated
                self.x.ends.grow();
                while self.x.ends.needs_growth(2 * self.nodes.len()) {
                    self.x.ends.grow();
                }
                for node in &self.nodes {
                    if let Some(leaf_hsp) = node.hsp.as_ref() {
                        let start = node.leftptr;
                        self.x.ends.insert(
                            start + leaf_hsp.query_offset,
                            leaf_hsp.subject_offset,
                            X_SALT_START,
                        );
                        self.x.ends.insert(
                            start + leaf_hsp.query_end,
                            leaf_hsp.subject_end,
                            X_SALT_END,
                        );
                    }
                }
            } else {
                self.x.ends.insert(
                    query_start + hsp.query_offset,
                    hsp.subject_offset,
                    X_SALT_START,
                );
                self.x
                    .ends
                    .insert(query_start + hsp.query_end, hsp.subject_end, X_SALT_END);
            }
        }

        // Start the insertion loop
        self.add_hsp_internal(
            new_index,
            region_start,
            region_end,
            &hsp,
            false,
            index_method,
        );
    }

    /// Legacy add_hsp without index_method (defaults to QueryAndSubject)
    #[inline]
    pub fn add_hsp_compat(&mut self, hsp: TreeHsp, query_context_offset: i32) {
        self.add_hsp(hsp, query_context_offset, IndexMethod::QueryAndSubject);
    }

    /// Internal helper for adding HSP after leaf node has been allocated
    /// NCBI reference: blast_itree.c:601-793
    fn add_hsp_internal(
        &mut self,
        new_index: usize,
        region_start: i32,
        region_end: i32,
        hsp: &TreeHsp,
        mut index_subject_range: bool,
        index_method: IndexMethod,
    ) {
        let mut root_index = 0usize;
        let mut current_region_start = region_start;
        let mut current_region_end = region_end;

        // NCBI reference: blast_itree.c:603 - while (1)
        loop {
            // NCBI reference: blast_itree.c:608-609
            let middle = {
                let node = &self.nodes[root_index];
                ((node.leftend as i64 + node.rightend as i64) / 2) as i32
            };

            // NCBI reference: blast_itree.c:611-633
            if current_region_end < middle {
                // New interval belongs in left subtree
                let leftptr = self.nodes[root_index].leftptr;

                if leftptr == 0 {
                    // No left child - attach new leaf directly
                    self.nodes[root_index].leftptr = new_index as i32;
                    return;
                }

                // Check if existing node is a leaf or internal
                let old_node_is_leaf = self.nodes[leftptr as usize].hsp.is_some();
                if !old_node_is_leaf {
                    // Descend to internal node
                    root_index = leftptr as usize;
                    continue;
                }

                // NCBI reference: blast_itree.c:703-793
                // Two leaves in same subtree - need to split
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:703-746
                // NCBI continues the insertion loop after splitting; avoid recursion to match.
                let mid_index = self.split_and_insert(
                    root_index,
                    leftptr as usize,
                    IntervalDirection::Left,
                    index_subject_range,
                    index_method,
                );
                root_index = mid_index;
                continue;
            }
            // NCBI reference: blast_itree.c:635-657
            else if current_region_start > middle {
                // New interval belongs in right subtree
                let rightptr = self.nodes[root_index].rightptr;

                if rightptr == 0 {
                    // No right child - attach new leaf directly
                    self.nodes[root_index].rightptr = new_index as i32;
                    return;
                }

                // Check if existing node is a leaf or internal
                let old_node_is_leaf = self.nodes[rightptr as usize].hsp.is_some();
                if !old_node_is_leaf {
                    // Descend to internal node
                    root_index = rightptr as usize;
                    continue;
                }

                // NCBI reference: blast_itree.c:703-793
                // Two leaves in same subtree - need to split
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:703-746
                // NCBI continues the insertion loop after splitting; avoid recursion to match.
                let mid_index = self.split_and_insert(
                    root_index,
                    rightptr as usize,
                    IntervalDirection::Right,
                    index_subject_range,
                    index_method,
                );
                root_index = mid_index;
                continue;
            }
            // NCBI reference: blast_itree.c:659-701
            else {
                // The new interval crosses the center of the node
                if index_subject_range
                    || index_method == IndexMethod::QueryOnly
                    || index_method == IndexMethod::QueryOnlyStrandIndifferent
                {
                    // NCBI reference: blast_itree.c:668-675
                    // If indexing subject offsets already, or only indexing query,
                    // prepend the new node to the list of "midpoint" nodes
                    self.nodes[new_index].midptr = self.nodes[root_index].midptr;
                    self.nodes[root_index].midptr = new_index as i32;
                    return;
                } else {
                    // NCBI reference: blast_itree.c:677-700
                    // Begin another tree at root_index, that indexes the subject range
                    index_subject_range = true;

                    let midptr = self.nodes[root_index].midptr;
                    if midptr == 0 {
                        // Create new subject-indexed subtree root
                        let mid_index = self.alloc_subject_root_node();
                        self.nodes[root_index].midptr = mid_index as i32;
                    }
                    root_index = self.nodes[root_index].midptr as usize;

                    // Switch from query range to subject range
                    current_region_start = hsp.subject_offset;
                    current_region_end = hsp.subject_end;
                    continue;
                }
            }
        }
    }

    /// Handle the case where two leaves collide in the same subtree.
    /// This creates a new internal node and reattaches the old leaf.
    ///
    /// NCBI reference: blast_itree.c:703-793
    fn split_and_insert(
        &mut self,
        parent_root_index: usize,
        old_leaf_index: usize,
        which_half: IntervalDirection,
        index_subject_range: bool,
        index_method: IndexMethod,
    ) -> usize {
        // NCBI reference: blast_itree.c:709-712
        // Allocate new internal node
        let mid_index = self.alloc_internal_node(parent_root_index, which_half);

        // Get the old HSP data before we modify anything
        let old_hsp = self.nodes[old_leaf_index].hsp.unwrap();
        let old_q_start = self.nodes[old_leaf_index].leftptr; // leftptr stores query_start for leaves

        // NCBI reference: blast_itree.c:717-720
        // Attach the new internal node to parent
        match which_half {
            IntervalDirection::Left => {
                self.nodes[parent_root_index].leftptr = mid_index as i32;
            }
            IntervalDirection::Right => {
                self.nodes[parent_root_index].rightptr = mid_index as i32;
            }
            IntervalDirection::Neither => {}
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:726-743
        // if (index_subject_range) {
        //     old_region_start = old_hsp->subject.offset;
        //     old_region_end = old_hsp->subject.end;
        // } else {
        //     if ( index_method == eQueryOnlyStrandIndifferent &&
        //          query_info->contexts[old_hsp->context].frame == -1 ) {
        //         q_start = s_GetQueryStrandOffset(query_info, old_hsp->context);
        //         old_region_end = q_start - old_hsp->query.offset;
        //         old_region_start = q_start - old_hsp->query.end;
        //     } else {
        //         old_region_start = nodes[old_index].leftptr + old_hsp->query.offset;
        //         old_region_end = nodes[old_index].leftptr + old_hsp->query.end;
        //     }
        // }
        let (old_region_start, old_region_end) = if index_subject_range {
            (old_hsp.subject_offset, old_hsp.subject_end)
        } else if index_method == IndexMethod::QueryOnlyStrandIndifferent && old_hsp.query_frame < 0
        {
            let q_start = old_hsp.query_context_offset;
            let old_region_end = q_start - old_hsp.query_offset;
            let old_region_start = q_start - old_hsp.query_end;
            (old_region_start, old_region_end)
        } else {
            (
                old_q_start + old_hsp.query_offset,
                old_q_start + old_hsp.query_end,
            )
        };

        // NCBI reference: blast_itree.c:747-748
        let middle = {
            let node = &self.nodes[mid_index];
            ((node.leftend as i64 + node.rightend as i64) / 2) as i32
        };

        // NCBI reference: blast_itree.c:749-791
        // Reattach old leaf to the new internal node
        if old_region_end < middle {
            // Old leaf belongs in left subtree of new node
            self.nodes[mid_index].leftptr = old_leaf_index as i32;
        } else if old_region_start > middle {
            // Old leaf belongs in right subtree of new node
            self.nodes[mid_index].rightptr = old_leaf_index as i32;
        } else {
            // Old leaf straddles both subtrees of new node
            if index_subject_range
                || index_method == IndexMethod::QueryOnly
                || index_method == IndexMethod::QueryOnlyStrandIndifferent
            {
                // NCBI reference: blast_itree.c:771
                self.nodes[mid_index].midptr = old_leaf_index as i32;
            } else {
                // NCBI reference: blast_itree.c:773-790
                // Need to create a new subject-indexed tree for the old leaf
                let mid_index2 = self.alloc_subject_root_node();
                self.nodes[mid_index].midptr = mid_index2 as i32;

                let old_s_start = old_hsp.subject_offset;
                let old_s_end = old_hsp.subject_end;

                let middle2 = {
                    let node = &self.nodes[mid_index2];
                    ((node.leftend as i64 + node.rightend as i64) / 2) as i32
                };

                if old_s_end < middle2 {
                    self.nodes[mid_index2].leftptr = old_leaf_index as i32;
                } else if old_s_start > middle2 {
                    self.nodes[mid_index2].rightptr = old_leaf_index as i32;
                } else {
                    self.nodes[mid_index2].midptr = old_leaf_index as i32;
                }
            }
        }

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:746
        // Return the new internal node so the caller can continue the loop.
        mid_index
    }

    /// Return the tree HSP that envelops the input HSP, if any.
    /// NCBI reference: blast_itree.c:930-995 BlastIntervalTreeContainsHSP
    ///
    /// ```c
    /// while (node->hsp == NULL) {
    ///     tmp_index = node->midptr;
    ///     if (tmp_index > 0) {
    ///         if (s_MidpointTreeContainsHSP(tree, tmp_index,
    ///                                       hsp, query_start,
    ///                                       min_diag_separation)) {
    ///             return TRUE;
    ///         }
    ///     }
    ///     middle = ((Int8) node->leftend + (Int8) node->rightend) / 2;
    ///     if (region_end < middle)
    ///         tmp_index = node->leftptr;
    ///     else if (region_start > middle)
    ///         tmp_index = node->rightptr;
    ///     if (tmp_index == 0)
    ///         return FALSE;
    ///     node = tree->nodes + tmp_index;
    /// }
    /// return s_HSPIsContained(hsp, query_start, node->hsp, node->leftptr,
    ///                         min_diag_separation);
    /// ```
    ///
    /// LOSAT returns the containing HSP for tracing only; containment traversal
    /// and truth value are still the NCBI rules above.
    pub fn contains_hsp(
        &self,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> bool {
        self.containing_hsp(hsp, query_context_offset, min_diag_separation)
            .is_some()
    }

    /// Return the tree HSP that envelops the input HSP, if any.
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:930-995
    /// See the snippet on `contains_hsp`; this method preserves that traversal
    /// and exposes the matched tree HSP for Phase 6 parity diagnostics.
    pub fn containing_hsp(
        &self,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> Option<TreeHsp> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:931-935,952,992-994
        // ```c
        // BlastIntervalTreeContainsHSP(const BlastIntervalTree *tree,
        //                              const BlastHSP *hsp,
        //                              const BlastQueryInfo *query_info,
        //                              Int4 min_diag_separation)
        // ...
        //     while (node->hsp == NULL) {
        // ...
        //     return s_HSPIsContained(hsp, query_start,
        //                             node->hsp, node->leftptr,
        //                             min_diag_separation);
        // ```
        // Dispatch point of LOSAT_X_ITREEFAST / LOSAT_X_ITREESHADOW. `containing_hsp_tree` below
        // ports this function. The grid answers the same yes/no question (mode 1). In mode 2
        // the walk also runs, `assert!` compares the two answers, and the walk's result is used.
        // EXPERIMENT (LOSAT_X_ITREEFAST / LOSAT_X_ITREESHADOW)
        if self.x.usable {
            let found = self.x_find_container(hsp, query_context_offset, min_diag_separation);
            if self.x.mode == 1 {
                return found;
            }
            let tree = self.containing_hsp_tree(hsp, query_context_offset, min_diag_separation);
            assert!(
                found.is_some() == tree.is_some(),
                "LOSAT_X_ITREESHADOW: grid says contained={}, tree says contained={}",
                found.is_some(),
                tree.is_some()
            );
            return tree;
        }
        self.containing_hsp_tree(hsp, query_context_offset, min_diag_separation)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:952-988
    /// ```c
    ///     while (node->hsp == NULL) {
    /// ...
    ///         tmp_index = node->midptr;
    ///         if (tmp_index > 0) {
    ///             if (s_MidpointTreeContainsHSP(tree, tmp_index,
    /// ...
    ///         middle = ((Int8) node->leftend + (Int8) node->rightend) / 2;
    ///         if (region_end < middle)
    ///             tmp_index = node->leftptr;
    ///         else if (region_start > middle)
    ///             tmp_index = node->rightptr;
    /// ```
    /// The previous body of `containing_hsp`, moved here unchanged: the port of the
    /// walk, which runs when no switch is set.
    /// The tree walk of `containing_hsp`.
    fn containing_hsp_tree(
        &self,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> Option<TreeHsp> {
        // NCBI reference: blast_itree.c:942
        if self.nodes.is_empty() {
            return None;
        }

        // NCBI reference: blast_itree.c:944-948
        let query_start = query_context_offset;
        let region_start = query_start + hsp.query_offset;
        let region_end = query_start + hsp.query_end;

        let mut node_idx = 0usize;

        // NCBI reference: blast_itree.c:950-994
        loop {
            let node = &self.nodes[node_idx];

            // NCBI reference: blast_itree.c:990-994
            // If this is a leaf node, check containment directly
            if node.hsp.is_some() {
                let tree_hsp = node.hsp.as_ref().unwrap();
                return Self::is_hsp_contained(
                    hsp,
                    query_start,
                    tree_hsp,
                    node.leftptr, // leftptr stores query_context_offset for leaves
                    min_diag_separation,
                )
                .then_some(*tree_hsp);
            }

            // NCBI reference: blast_itree.c:960-967
            // Check midpoint tree first (contains all HSPs straddling this node)
            if node.midptr > 0 {
                if let Some(tree_hsp) = self.midpoint_tree_containing_hsp(
                    node.midptr as usize,
                    hsp,
                    query_start,
                    min_diag_separation,
                ) {
                    return Some(tree_hsp);
                }
            }

            // NCBI reference: blast_itree.c:973-985
            // Descend to appropriate subtree
            let middle = ((node.leftend as i64 + node.rightend as i64) / 2) as i32;

            let next_idx = if region_end < middle {
                node.leftptr
            } else if region_start > middle {
                node.rightptr
            } else {
                // Input straddles middle - all potential containers already checked
                0
            };

            // NCBI reference: blast_itree.c:987
            if next_idx == 0 {
                return None;
            }

            node_idx = next_idx as usize;
        }
    }

    /// Internal recursive helper for containment check (legacy - kept for reference)
    #[allow(dead_code)]
    fn contains_hsp_recursive(
        &self,
        root_index: usize,
        region_start: i32,
        region_end: i32,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> bool {
        let node = &self.nodes[root_index];

        // Check if leaf node
        if let Some(ref tree_hsp) = node.hsp {
            // Leaf node - check containment directly
            return Self::is_hsp_contained(
                hsp,
                query_context_offset,
                tree_hsp,
                node.leftptr, // leftptr stores query_context_offset for leaves
                min_diag_separation,
            );
        }

        // Internal node - first check midpoint subtree
        let midptr = node.midptr;
        if midptr > 0 {
            if self.midpoint_tree_contains_hsp(
                midptr as usize,
                hsp,
                query_context_offset,
                min_diag_separation,
            ) {
                return true;
            }
        }

        // Descend to appropriate subtree
        let middle = ((node.leftend as i64 + node.rightend as i64) / 2) as i32;

        let next_idx = if region_end < middle {
            node.leftptr
        } else if region_start > middle {
            node.rightptr
        } else {
            // Straddles middle - all potential containers already checked
            0
        };

        if next_idx == 0 {
            return false;
        }

        self.contains_hsp_recursive(
            next_idx as usize,
            region_start,
            region_end,
            hsp,
            query_context_offset,
            min_diag_separation,
        )
    }

    fn midpoint_tree_contains_hsp(
        &self,
        root_index: usize,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> bool {
        self.midpoint_tree_containing_hsp(
            root_index,
            hsp,
            query_context_offset,
            min_diag_separation,
        )
        .is_some()
    }

    /// Check midpoint subtree (subject-indexed).
    /// NCBI reference: blast_itree.c:864-927 s_MidpointTreeContainsHSP
    ///
    /// ```c
    /// while (tmp_index != 0) {
    ///     if (s_HSPIsContained(hsp, query_start, node->hsp,
    ///                          node->leftptr, min_diag_separation)) {
    ///         return TRUE;
    ///     }
    ///     tmp_index = node->midptr;
    /// }
    /// middle = ((Int8) node->leftend + (Int8) node->rightend) / 2;
    /// if (region_end < middle)
    ///     tmp_index = node->leftptr;
    /// else if (region_start > middle)
    ///     tmp_index = node->rightptr;
    /// if (tmp_index == 0)
    ///     return FALSE;
    /// ```
    ///
    /// LOSAT returns the containing midpoint-list HSP for trace output while
    /// preserving the same NCBI walk and predicate.
    fn midpoint_tree_containing_hsp(
        &self,
        root_index: usize,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> Option<TreeHsp> {
        let region_start = hsp.subject_offset;
        let region_end = hsp.subject_end;
        let mut current_idx = root_index;

        loop {
            let node = &self.nodes[current_idx];

            // Check if leaf
            if let Some(ref tree_hsp) = node.hsp {
                return Self::is_hsp_contained(
                    hsp,
                    query_context_offset,
                    tree_hsp,
                    node.leftptr,
                    min_diag_separation,
                )
                .then_some(*tree_hsp);
            }

            // Check all HSPs in midpoint list
            let mut tmp_idx = node.midptr;
            while tmp_idx != 0 {
                let tmp_node = &self.nodes[tmp_idx as usize];
                if let Some(ref tree_hsp) = tmp_node.hsp {
                    if Self::is_hsp_contained(
                        hsp,
                        query_context_offset,
                        tree_hsp,
                        tmp_node.leftptr,
                        min_diag_separation,
                    ) {
                        return Some(*tree_hsp);
                    }
                }
                tmp_idx = tmp_node.midptr;
            }

            // Descend to subtree
            let middle = ((node.leftend as i64 + node.rightend as i64) / 2) as i32;

            let next_idx = if region_end < middle {
                node.leftptr
            } else if region_start > middle {
                node.rightptr
            } else {
                0
            };

            if next_idx == 0 {
                return None;
            }

            current_idx = next_idx as usize;
        }
    }

    /// Check if one HSP is contained within another
    /// NCBI reference: blast_itree.c:809-847 s_HSPIsContained
    ///
    /// # Arguments
    /// * `in_hsp` - The HSP being checked for containment
    /// * `in_q_start` - Query context offset for in_hsp
    /// * `tree_hsp` - The HSP from the tree (potential container)
    /// * `tree_q_start` - Query context offset for tree_hsp
    /// * `min_diag_separation` - Diagonal separation threshold
    fn is_hsp_contained(
        in_hsp: &TreeHsp,
        in_q_start: i32,
        tree_hsp: &TreeHsp,
        tree_q_start: i32,
        min_diag_separation: i32,
    ) -> bool {
        // NCBI reference: blast_itree.c:819
        // Check if alignments are from different query sequences or query strands
        if in_q_start != tree_q_start {
            return false;
        }

        // NCBI reference: blast_itree.c:822-831
        // All conditions must be true:
        // 1. in_hsp->score <= tree_hsp->score
        // 2. SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame)
        // 3. CONTAINED_IN_HSP for start endpoint
        // 4. CONTAINED_IN_HSP for end endpoint
        if in_hsp.score <= tree_hsp.score
            && in_hsp.subject_frame_sign.signum() == tree_hsp.subject_frame_sign.signum()
            && Self::contained_in_hsp(
                tree_hsp.query_offset,
                tree_hsp.query_end,
                in_hsp.query_offset,
                tree_hsp.subject_offset,
                tree_hsp.subject_end,
                in_hsp.subject_offset,
            )
            && Self::contained_in_hsp(
                tree_hsp.query_offset,
                tree_hsp.query_end,
                in_hsp.query_end,
                tree_hsp.subject_offset,
                tree_hsp.subject_end,
                in_hsp.subject_end,
            )
        {
            // NCBI reference: blast_itree.c:833-834
            if min_diag_separation == 0 {
                return true;
            }

            // NCBI reference: blast_itree.c:836-843
            // MB_HSP_CLOSE at start OR end
            if Self::mb_hsp_close(
                tree_hsp.query_offset,
                tree_hsp.subject_offset,
                in_hsp.query_offset,
                in_hsp.subject_offset,
                min_diag_separation,
            ) || Self::mb_hsp_close(
                tree_hsp.query_end,
                tree_hsp.subject_end,
                in_hsp.query_end,
                in_hsp.subject_end,
                min_diag_separation,
            ) {
                return true;
            }
        }

        false
    }

    /// CONTAINED_IN_HSP macro
    /// NCBI reference: blast_gapalign_priv.h:120-121
    /// #define CONTAINED_IN_HSP(a,b,c,d,e,f) ((a <= c && b >= c) && (d <= f && e >= f))
    #[inline]
    fn contained_in_hsp(
        hsp_q_start: i32,
        hsp_q_end: i32,
        point_q: i32,
        hsp_s_start: i32,
        hsp_s_end: i32,
        point_s: i32,
    ) -> bool {
        (hsp_q_start <= point_q && hsp_q_end >= point_q)
            && (hsp_s_start <= point_s && hsp_s_end >= point_s)
    }

    /// MB_HSP_CLOSE macro
    /// NCBI reference: blast_gapalign_priv.h:123-124
    /// #define MB_HSP_CLOSE(q1, s1, q2, s2, c) (ABS(((q1)-(s1)) - ((q2)-(s2))) < c)
    #[inline]
    fn mb_hsp_close(q1: i32, s1: i32, q2: i32, s2: i32, separation: i32) -> bool {
        ((q1 - s1) - (q2 - s2)).abs() < separation
    }

    /// Get number of nodes in the tree
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tree_hsp(
        query_offset: i32,
        query_end: i32,
        subject_offset: i32,
        subject_end: i32,
        score: i32,
        subject_frame_sign: i32,
    ) -> TreeHsp {
        TreeHsp {
            query_offset,
            query_end,
            subject_offset,
            subject_end,
            score,
            query_frame: 1,
            query_length: 1000,
            query_context_offset: 0,
            subject_frame_sign,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:577-584,931-935,992-994
    /// ```c
    ///         if (s_IntervalTreeHasHSPEndpoint(tree, hsp, query_start,
    ///                                          eIntervalTreeLeft)) {
    /// ...
    /// BlastIntervalTreeContainsHSP(const BlastIntervalTree *tree,
    /// ...
    ///     return s_HSPIsContained(hsp, query_start,
    /// ```
    /// The test feeds the same random HSPs to trees that use the walk, the side indexes, and
    /// both. It compares every containment answer and the final trees.
    /// EXPERIMENT (LOSAT_X_ITREEFAST): three trees are fed the same random
    /// HSPs, one walking the tree only, one using the side indexes, and one
    /// doing both and comparing inside every call. Every containment answer
    /// must agree, and the trees must end up with the same nodes and the same
    /// number of unlinked HSPs. The HSPs are generated so that containment,
    /// shared start points, shared end points and replacements all occur.
    #[test]
    fn x_side_indexes_agree_with_tree_walk_on_random_hsps() {
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
        let trials: usize = std::env::var("LOSAT_FUZZ_CASES")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(60);
        let mut rng = Rng(0xD6E8_FEB8_6659_FD93);
        let mut contained = 0usize;
        let mut queries = 0usize;
        let mut unlinked = 0usize;
        for trial in 0..trials {
            let q_len = [300i32, 5_000, 80_000, 3_000_000][trial % 4];
            let s_len = [400i32, 9_000, 120_000, 2_000_000][(trial / 4) % 4];
            let separation = [0i32, 6, 50][trial % 3];
            let new_tree = |mode: u8| {
                X_ITREE_TEST_MODE.with(|m| m.set(Some(mode)));
                let tree = BlastIntervalTree::new(0, 2 * q_len + 1, 0, s_len + 1);
                X_ITREE_TEST_MODE.with(|m| m.set(None));
                tree
            };
            let mut walk = new_tree(0);
            let mut indexed = new_tree(1);
            let mut both = new_tree(2);
            let mut added: Vec<(TreeHsp, i32)> = Vec::new();
            let count = 200 + rng.below(1500) as usize;
            for _ in 0..count {
                let context_offset = if rng.below(4) == 0 { q_len + 1 } else { 0 };
                let kind = rng.below(10);
                let (mut hsp, context_offset) = if kind < 3 && !added.is_empty() {
                    // a piece of an earlier HSP (often contained in it)
                    let (base, offset) = added[rng.below(added.len() as u64) as usize];
                    let q_span = (base.query_end - base.query_offset).max(1);
                    let s_span = (base.subject_end - base.subject_offset).max(1);
                    let a = rng.below(q_span as u64) as i32;
                    let b = a + rng.below((q_span - a) as u64 + 1) as i32;
                    let shift = rng.below(5) as i32 - 2;
                    let mut piece = base;
                    piece.query_offset = base.query_offset + a;
                    piece.query_end = base.query_offset + b;
                    piece.subject_offset = (base.subject_offset + a.min(s_span) + shift).max(0);
                    piece.subject_end =
                        (base.subject_offset + b.min(s_span) + shift).max(piece.subject_offset);
                    piece.score = (base.score - rng.below(40) as i32 + 5).max(1);
                    (piece, offset)
                } else if kind < 5 && !added.is_empty() {
                    // shares the start or the end point of an earlier HSP
                    let (base, offset) = added[rng.below(added.len() as u64) as usize];
                    let mut other = base;
                    let change = 1 + rng.below(200) as i32;
                    if rng.below(2) == 0 {
                        other.query_end = (base.query_end + change - 100).max(base.query_offset);
                        other.subject_end =
                            (base.subject_end + change - 100).max(base.subject_offset);
                    } else {
                        other.query_offset =
                            (base.query_offset - change + 100).clamp(0, base.query_end);
                        other.subject_offset =
                            (base.subject_offset - change + 100).clamp(0, base.subject_end);
                    }
                    other.score = (base.score + rng.below(41) as i32 - 20).max(1);
                    (other, offset)
                } else {
                    let longest = [50u64, 800, 20_000][rng.below(3) as usize];
                    let length = 5 + rng.below(longest) as i32;
                    let length = length.min(q_len - 1).min(s_len - 1).max(1);
                    let query_offset = rng.below((q_len - length) as u64) as i32;
                    let subject_offset = rng.below((s_len - length) as u64) as i32;
                    let skew = rng.below(7) as i32 - 3;
                    (
                        TreeHsp {
                            query_offset,
                            query_end: query_offset + length,
                            subject_offset,
                            subject_end: (subject_offset + length + skew).max(subject_offset),
                            score: length + rng.below(30) as i32,
                            query_frame: 1,
                            query_length: q_len,
                            query_context_offset: context_offset,
                            subject_frame_sign: if rng.below(8) == 0 { -1 } else { 1 },
                        },
                        context_offset,
                    )
                };
                hsp.query_context_offset = context_offset;
                // inside the ranges the tree was built for
                hsp.query_offset = hsp.query_offset.clamp(0, q_len);
                hsp.query_end = hsp.query_end.clamp(hsp.query_offset, q_len);
                hsp.subject_offset = hsp.subject_offset.clamp(0, s_len);
                hsp.subject_end = hsp.subject_end.clamp(hsp.subject_offset, s_len);

                let a = walk.contains_hsp(&hsp, context_offset, separation);
                let b = indexed.contains_hsp(&hsp, context_offset, separation);
                let c = both.contains_hsp(&hsp, context_offset, separation);
                assert!(
                    a == b && a == c,
                    "trial {trial}: walk={a} indexed={b} both={c} for {hsp:?}"
                );
                queries += 1;
                contained += a as usize;
                if !a || rng.below(4) == 0 {
                    walk.add_hsp(hsp, context_offset, IndexMethod::QueryAndSubject);
                    indexed.add_hsp(hsp, context_offset, IndexMethod::QueryAndSubject);
                    both.add_hsp(hsp, context_offset, IndexMethod::QueryAndSubject);
                    added.push((hsp, context_offset));
                }
            }
            assert_eq!(walk.node_count(), indexed.node_count(), "trial {trial}");
            assert_eq!(walk.node_count(), both.node_count(), "trial {trial}");
            assert_eq!(walk.x.removed, indexed.x.removed, "trial {trial}");
            assert_eq!(walk.x.removed, both.x.removed, "trial {trial}");
            assert!(
                indexed.x.usable && both.x.usable,
                "trial {trial}: index abandoned"
            );
            unlinked += walk.x.removed;
            // every HSP ever offered, asked again of the final trees
            for &(hsp, context_offset) in &added {
                let a = walk.contains_hsp(&hsp, context_offset, separation);
                let b = indexed.contains_hsp(&hsp, context_offset, separation);
                let c = both.contains_hsp(&hsp, context_offset, separation);
                assert!(
                    a == b && a == c,
                    "trial {trial}: final walk={a} indexed={b} both={c}"
                );
            }
        }
        // the generator must exercise both answers and the unlinking
        assert!(
            contained * 10 > queries,
            "{contained} of {queries} contained"
        );
        assert!(
            (queries - contained) * 10 > queries,
            "{contained} of {queries} contained"
        );
        assert!(
            unlinked > trials,
            "{unlinked} HSPs unlinked in {trials} trials"
        );
    }

    #[test]
    fn test_interval_tree_basic() {
        let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);

        // Add an HSP
        let hsp1 = TreeHsp {
            query_offset: 100,
            query_end: 200,
            subject_offset: 500,
            subject_end: 600,
            score: 100,
            query_frame: 1,
            query_length: 1000,
            query_context_offset: 0,
            subject_frame_sign: 1,
        };
        tree.add_hsp(hsp1, 0, IndexMethod::QueryAndSubject);

        // Check containment of a smaller HSP within the first
        let hsp2 = TreeHsp {
            query_offset: 120,
            query_end: 180,
            subject_offset: 520,
            subject_end: 580,
            score: 50,
            query_frame: 1,
            query_length: 1000,
            query_context_offset: 0,
            subject_frame_sign: 1,
        };

        // hsp2 is contained within hsp1 (lower score, both endpoints inside)
        assert!(tree.contains_hsp(&hsp2, 0, 0));
    }

    #[test]
    fn test_interval_tree_not_contained() {
        let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);

        let hsp1 = TreeHsp {
            query_offset: 100,
            query_end: 200,
            subject_offset: 500,
            subject_end: 600,
            score: 100,
            query_frame: 1,
            query_length: 1000,
            query_context_offset: 0,
            subject_frame_sign: 1,
        };
        tree.add_hsp(hsp1, 0, IndexMethod::QueryAndSubject);

        // HSP with higher score should not be contained
        let hsp2 = TreeHsp {
            query_offset: 120,
            query_end: 180,
            subject_offset: 520,
            subject_end: 580,
            score: 150, // Higher score
            query_frame: 1,
            query_length: 1000,
            query_context_offset: 0,
            subject_frame_sign: 1,
        };

        assert!(!tree.contains_hsp(&hsp2, 0, 0));
    }

    #[test]
    fn test_interval_tree_different_strand() {
        let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);

        let hsp1 = TreeHsp {
            query_offset: 100,
            query_end: 200,
            subject_offset: 500,
            subject_end: 600,
            score: 100,
            query_frame: 1,
            query_length: 1000,
            query_context_offset: 0,
            subject_frame_sign: 1, // Forward strand
        };
        tree.add_hsp(hsp1, 0, IndexMethod::QueryAndSubject);

        // HSP on different strand should not be contained
        let hsp2 = TreeHsp {
            query_offset: 120,
            query_end: 180,
            subject_offset: 520,
            subject_end: 580,
            score: 50,
            query_frame: 1,
            query_length: 1000,
            query_context_offset: 0,
            subject_frame_sign: -1, // Reverse strand
        };

        assert!(!tree.contains_hsp(&hsp2, 0, 0));
    }

    #[test]
    fn test_mb_hsp_close() {
        // Same diagonal - should be close
        assert!(BlastIntervalTree::mb_hsp_close(100, 200, 150, 250, 50));

        // Different diagonals - should not be close
        assert!(!BlastIntervalTree::mb_hsp_close(100, 200, 150, 300, 50));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:269-301
    // ```c
    // if (in_hsp->score > tree_hsp->score)
    //     return in_hsp;
    // ...
    // /* for equal scores, pick the shorter HSP */
    // in_q_length = in_hsp->query.end - in_hsp->query.offset;
    // ```
    #[test]
    fn test_interval_tree_common_endpoint_prefers_shorter_equal_score_hsp() {
        let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);
        let longer_hsp = make_tree_hsp(100, 200, 500, 600, 100, 1);
        let shorter_hsp = make_tree_hsp(100, 180, 500, 580, 100, 1);
        let long_only_candidate = make_tree_hsp(181, 190, 581, 590, 50, 1);

        tree.add_hsp(longer_hsp, 0, IndexMethod::QueryAndSubject);
        tree.add_hsp(shorter_hsp, 0, IndexMethod::QueryAndSubject);

        assert!(!tree.contains_hsp(&long_only_candidate, 0, 0));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:819-823
    // ```c
    // if (in_q_start != tree_q_start)
    //     return FALSE;
    // ...
    // SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame)
    // ```
    #[test]
    fn test_interval_tree_containment_uses_subject_frame_sign() {
        let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);
        let tree_hsp = make_tree_hsp(100, 200, 500, 600, 100, 3);
        let contained_hsp = make_tree_hsp(120, 180, 520, 580, 50, 1);

        tree.add_hsp(tree_hsp, 0, IndexMethod::QueryAndSubject);

        assert!(tree.contains_hsp(&contained_hsp, 0, 0));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:819-847
    // ```c
    // if (in_q_start != tree_q_start)
    //     return FALSE;
    // if (in_hsp->score <= tree_hsp->score &&
    //     SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame) &&
    //     CONTAINED_IN_HSP(... in_hsp->query.offset ... in_hsp->subject.offset) &&
    //     CONTAINED_IN_HSP(... in_hsp->query.end ... in_hsp->subject.end)) {
    //     if (min_diag_separation == 0)
    //         return TRUE;
    //     if (MB_HSP_CLOSE(... start ...) || MB_HSP_CLOSE(... end ...))
    //         return TRUE;
    // }
    // ```
    #[test]
    fn test_containing_hsp_reports_ncbi_score_and_endpoint_semantics() {
        let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);
        let tree_hsp = make_tree_hsp(100, 200, 500, 600, 100, 1);
        tree.add_hsp(tree_hsp, 0, IndexMethod::QueryAndSubject);

        let lower_score_inside = make_tree_hsp(120, 180, 520, 580, 50, 1);
        let equal_score_inside = make_tree_hsp(125, 175, 525, 575, 100, 1);
        let higher_score_inside = make_tree_hsp(120, 180, 520, 580, 101, 1);
        let endpoint_outside = make_tree_hsp(120, 220, 520, 620, 50, 1);
        let opposite_subject_sign = make_tree_hsp(120, 180, 520, 580, 50, -1);

        let container = tree
            .containing_hsp(&lower_score_inside, 0, 0)
            .expect("lower-score HSP with both endpoints inside must be contained");
        assert_eq!(container.query_offset, tree_hsp.query_offset);
        assert_eq!(container.subject_offset, tree_hsp.subject_offset);
        assert!(tree.contains_hsp(&equal_score_inside, 0, 0));
        assert!(!tree.contains_hsp(&higher_score_inside, 0, 0));
        assert!(!tree.contains_hsp(&endpoint_outside, 0, 0));
        assert!(!tree.contains_hsp(&opposite_subject_sign, 0, 0));
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_itree.c:819-843
    // ```c
    // if (in_q_start != tree_q_start)
    //     return FALSE;
    // ...
    // if (min_diag_separation == 0)
    //     return TRUE;
    // if (MB_HSP_CLOSE(tree_hsp->query.offset, tree_hsp->subject.offset,
    //                  in_hsp->query.offset, in_hsp->subject.offset,
    //                  min_diag_separation) ||
    //     MB_HSP_CLOSE(tree_hsp->query.end, tree_hsp->subject.end,
    //                  in_hsp->query.end, in_hsp->subject.end,
    //                  min_diag_separation)) {
    //     return TRUE;
    // }
    // ```
    #[test]
    fn test_containment_respects_query_context_and_diag_separation() {
        let mut tree = BlastIntervalTree::new(0, 3000, 0, 2000);
        let mut tree_hsp = make_tree_hsp(100, 200, 500, 600, 100, 1);
        tree_hsp.query_context_offset = 1001;
        tree.add_hsp(tree_hsp, 1001, IndexMethod::QueryAndSubject);

        let same_context_far_diag = make_tree_hsp(150, 160, 500, 510, 50, 1);
        let same_context_close_diag = make_tree_hsp(149, 159, 500, 510, 50, 1);
        let different_context = make_tree_hsp(149, 159, 500, 510, 50, 1);

        assert!(tree.contains_hsp(&same_context_far_diag, 1001, 0));
        assert!(!tree.contains_hsp(&same_context_far_diag, 1001, 50));
        assert!(tree.contains_hsp(&same_context_close_diag, 1001, 50));
        assert!(!tree.contains_hsp(&different_context, 0, 0));
    }
}
