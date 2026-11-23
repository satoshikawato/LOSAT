//! NCBI BLAST Interval Tree implementation for HSP containment checking
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c
//!
//! This module implements the interval tree data structure used by NCBI BLAST
//! to efficiently check if an HSP is contained within existing HSPs.
//!
//! The tree is organized as follows:
//! - Primary index: query offset range
//! - Secondary index (midpoint lists): subject offset range
//! - Nodes can be internal (with children) or leaves (with HSP)

#![allow(dead_code)]

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.h:83-89
// Indexing modes: query only (`eQueryOnly`), query and subject (`eQueryAndSubject`),
// or query only ignoring strand (`eQueryOnlyStrandIndifferent`).
/// Which offsets an interval tree uses as the key for each HSP inserted into it
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum IndexMethod {
    /// Index by query offset only
    QueryOnly,
    /// Index by query and then by subject offset
    QueryAndSubject,
    /// Index by query offset only, strand indifferent (NCBI's `eQueryOnlyStrandIndifferent`;
    /// renamed so that the variants do not all start with `Query`, clippy `enum_variant_names`)
    StrandIndifferentQueryOnly,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:41-45
// Tree navigation directions: left (`eIntervalTreeLeft`), right (`eIntervalTreeRight`),
// or neither (`eIntervalTreeNeither`).
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:241-301 (s_HSPsHaveCommonEndpoint returns pointer to better HSP)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum EndpointResult {
    /// Input HSP is better - keep it
    KeepInput,
    /// Tree HSP is better - keep it
    KeepTree,
}

/// HSP data stored in interval tree nodes
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c - tree stores pointers to BlastHSP
#[derive(Clone, Copy, Debug)]
pub struct TreeHsp {
    pub query_offset: i32,
    pub query_end: i32,
    pub subject_offset: i32,
    pub subject_end: i32,
    pub score: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:539-546
    pub query_frame: i32,
    pub query_length: i32,
    /// Query strand offset (for context comparison; compares incoming and tree query starts).
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:819
    pub query_context_offset: i32,
    /// Subject frame sign (positive = forward, negative = reverse; sign of subject frame).
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:823
    pub subject_frame_sign: i32,
}

/// Interval tree node
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:48-58 SIntervalNode
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.h:56
    /// The HSP of a leaf, as 1 + its index in `BlastIntervalTree::hsps`; 0 (NULL) for internal
    /// nodes. A node is then 28 bytes, as NCBI's (with an 8-byte pointer) is 32, instead of 64
    /// with the HSP inside it.
    hsp: u32,
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:359-362,398-404
    /// NCBI unlinks such a leaf from the tree. This flag records that fact for the grid index,
    /// which still lists the leaf; a flagged leaf is never returned as a container.
    /// No NCBI counterpart: this leaf was unlinked from the tree.
    x_dead: bool,
}

// No NCBI counterpart: the node stays 28 bytes (memory only).
const _: () = assert!(std::mem::size_of::<IntervalNode>() == 28);

impl IntervalNode {
    fn new_internal(leftend: i32, rightend: i32) -> Self {
        Self {
            leftend,
            rightend,
            leftptr: 0,
            midptr: 0,
            rightptr: 0,
            hsp: 0,
            x_dead: false,
        }
    }

    /// `hsp` = 1 + the HSP's index in `BlastIntervalTree::hsps`.
    fn new_leaf(hsp: u32, query_context_offset: i32) -> Self {
        Self {
            leftend: 0,
            rightend: 0,
            leftptr: query_context_offset, // NCBI stores q_start offset here for leaves
            midptr: 0,
            rightptr: 0,
            hsp,
            x_dead: false,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:952
    /// Whether the node is a leaf (its `hsp` is not NULL).
    #[inline]
    fn is_leaf(&self) -> bool {
        self.hsp != 0
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:822-831,952-994
// The side indexes below answer what `BlastIntervalTreeContainsHSP` answers (does some HSP
// in the tree satisfy `s_HSPIsContained` for the input) and what
// `s_IntervalTreeHasHSPEndpoint` answers (is there an HSP with the same end point), without
// the walk. The tree stays the master copy and the predicates are the unchanged ones. The
// argument above explains why the answers are the same.
// ---------------------------------------------------------------------------
// No NCBI counterpart: two side indexes that
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
// index is off whenever tracing is on (the walk answers then).
// ---------------------------------------------------------------------------

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:931-935,952
/// Whether the side indexes answer a containment query in place of the walk above. They are
/// off whenever tracing is on, because the index may return a different container than the
/// walk when there are several; the walk then answers.
fn x_itree_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| !super::tracing::enabled())
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:822-831
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

// No NCBI counterpart: grid cells per HSP of a tree made for a known number of HSPs and the most
// cells such a tree gets (`x_grid_cells`); time and memory of the grid index only.
const X_GRID_CELLS_PER_HSP: usize = 64;
const X_GRID_MAX_CELLS_FOR_HSPS: usize = 1 << 12;
const _: () = assert!(X_GRID_MAX_CELLS_FOR_HSPS <= X_MAX_CELLS);

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
/// No NCBI counterpart: the grid limit (`BlastIntervalTree::reset_with_bounds_cells`) of a tree
/// that holds at most one HSP per initial HSP, as the tree NCBI makes for each subject (BLASTX) or
/// unit (TBLASTN) does: `X_GRID_CELLS_PER_HSP` cells per initial HSP, at most
/// `X_GRID_MAX_CELLS_FOR_HSPS` (four bytes per cell, zeroed on the first HSP). Such a tree is made
/// and dropped by each thread for each subject or unit; a grid of up to 256 KiB each time set how
/// far each thread's glibc arena grew. The count only sizes the grid: a smaller grid lists more
/// leaves per cell and answers the same, whatever the tree holds.
pub(crate) fn x_grid_cells(initial_hsps: usize) -> usize {
    initial_hsps
        .saturating_mul(X_GRID_CELLS_PER_HSP)
        .clamp(X_GRID_CELLS_PER_HSP, X_GRID_MAX_CELLS_FOR_HSPS)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:440-443,452-461
/// `s_IntervalTreeHasHSPEndpoint` only acts on tree HSPs with the same start (or end)
/// point as the input. This set holds the start and end points of the HSPs added, as a
/// hash bit table. "May contain" is false only when no such HSP was added; a false
/// "true" only makes the unchanged walk run.
/// A set of (query, subject) points that may report points it does not hold
/// but never misses one it holds: a bit table with two bits per point.
#[derive(Clone)]
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
        // The size follows the table's length only; the allocation is kept between trees.
        let words = (self.bits.len() * 2).max(128);
        self.bits.clear();
        self.bits.resize(words, 0);
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
// The state of the side indexes, tied to the same query and subject ranges as the tree
// (`q_min`, `s_min`) and emptied whenever the tree is emptied. It is derived from the
// HSPs in the tree and holds nothing else.
#[derive(Clone)]
struct XTreeIndex {
    /// Whether the index is used at all (off while tracing is on).
    on: bool,
    /// False once the tree holds an HSP this index does not describe.
    usable: bool,
    shift: u32,
    q_min: i32,
    s_min: i32,
    nq: usize,
    ns: usize,
    /// Per cell: 1 + index of its first entry, 0 = empty. Grown on the first HSP
    /// that needs the size, and all zero between trees: `clear` resets the cells in
    /// `used_cells`, so a new tree does not fill the grid again.
    head: Vec<u32>,
    /// The cells with a non-zero head, for `clear`.
    used_cells: Vec<u32>,
    entries: Vec<XCellEntry>,
    overflow: Vec<u32>,
    ends: XPointSet,
}

impl XTreeIndex {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:146-170
    /// Makes an empty index for the ranges the tree is made for.
    fn new(q_min: i32, q_max: i32, s_min: i32, s_max: i32) -> Self {
        let on = x_itree_on();
        let mut index = Self {
            on,
            usable: on,
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
        };
        index.set_bounds(q_min, q_max, s_min, s_max);
        index
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:166-170
    /// The index follows `reset_with_bounds`: new ranges, empty index; the grid is sized
    /// from the ranges (on the first HSP, in the allocation kept from earlier ranges).
    fn set_bounds(&mut self, q_min: i32, q_max: i32, s_min: i32, s_max: i32) {
        self.set_bounds_cells(q_min, q_max, s_min, s_max, X_MAX_CELLS);
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:166-170
    /// `set_bounds` with the grid limited to `max_cells` cells. A coarser grid still lists, in
    /// the cell of a start point, every leaf whose box holds that point (and more leaves to
    /// test), so the index answers the same.
    #[inline(always)]
    fn set_bounds_cells(
        &mut self,
        q_min: i32,
        q_max: i32,
        s_min: i32,
        s_max: i32,
        max_cells: usize,
    ) {
        let q_span = (q_max as i64 - q_min as i64).max(0) as usize;
        let s_span = (s_max as i64 - s_min as i64).max(0) as usize;
        let mut shift = 6u32;
        while ((q_span >> shift) + 1).saturating_mul((s_span >> shift) + 1) > max_cells {
            shift += 1;
        }
        self.shift = shift;
        self.q_min = q_min;
        self.s_min = s_min;
        self.nq = (q_span >> shift) + 1;
        self.ns = (s_span >> shift) + 1;
        self.clear();
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:195-201
    /// The index is emptied exactly when the tree is reset.
    fn clear(&mut self) {
        for &cell in &self.used_cells {
            self.head[cell as usize] = 0;
        }
        self.used_cells.clear();
        self.entries.clear();
        self.overflow.clear();
        self.ends.clear();
        self.usable = self.on;
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
        if self.head.len() < self.nq * self.ns {
            self.head.resize(self.nq * self.ns, 0);
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.h:73-80
#[derive(Clone)]
pub struct BlastIntervalTree {
    /// Array of nodes (index 0 is root for query-indexed tree)
    nodes: Vec<IntervalNode>,
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:597-599
    /// The HSPs the leaves point to (NCBI's leaves point to the caller's `BlastHSP`s), in the
    /// order the leaves were made; a leaf's `hsp` is 1 + its index here.
    hsps: Vec<TreeHsp>,
    /// Minimum subject offset
    s_min: i32,
    /// Maximum subject offset
    s_max: i32,
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:158-167
    /// The tree owns its side indexes; they are built from, and follow, the node array above.
    x: XTreeIndex,
}

impl BlastIntervalTree {
    /// Initialize a new interval tree
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:96-166 Blast_IntervalTreeInit
    ///
    /// # Arguments
    /// * `q_min` - Minimum query offset (usually 0)
    /// * `q_max` - Maximum query offset (query_length + 1)
    /// * `s_min` - Minimum subject offset (usually 0)
    /// * `s_max` - Maximum subject offset (subject_length + 1)
    pub fn new(q_min: i32, q_max: i32, s_min: i32, s_max: i32) -> Self {
        let mut tree = Self {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:135-165
            nodes: Vec::with_capacity(100),
            hsps: Vec::new(),
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:66-73
    /// Rust keeps the same node indices and insertion order, but can reserve
    /// enough Vec storage from the observed HSP count to avoid repeated reallocs.
    pub fn reserve_nodes_for_hsps(&mut self, hsp_count: usize) {
        let target = hsp_count.saturating_mul(3).saturating_add(1).max(100);
        if target > self.nodes.capacity() {
            self.nodes.reserve(target - self.nodes.capacity());
        }
        if hsp_count > self.hsps.capacity() {
            self.hsps.reserve(hsp_count - self.hsps.len());
        }
    }

    /// Reset the tree for reuse
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:168-193 Blast_IntervalTreeReset
    pub fn reset(&mut self) {
        // Keep only the root node, reset its children
        if !self.nodes.is_empty() {
            let leftend = self.nodes[0].leftend;
            let rightend = self.nodes[0].rightend;
            self.nodes.clear();
            self.nodes
                .push(IntervalNode::new_internal(leftend, rightend));
        }
        self.hsps.clear();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:195-201
        // Reset point: the side indexes are emptied with the tree.
        self.x.clear();
    }

    /// Reset the tree bounds for reuse.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:146-170
    pub fn reset_with_bounds(&mut self, q_min: i32, q_max: i32, s_min: i32, s_max: i32) {
        self.s_min = s_min;
        self.s_max = s_max;
        self.nodes.clear();
        self.nodes.push(IntervalNode::new_internal(q_min, q_max));
        self.hsps.clear();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:166-170
        // Re-initialisation point: the side indexes are rebuilt for the new ranges.
        self.x.set_bounds(q_min, q_max, s_min, s_max);
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:146-170
    /// `reset_with_bounds` with the grid index limited to `max_cells` cells (`X_MAX_CELLS`
    /// there). The tree is the same; the grid only answers the walk's question, so its size
    /// changes the time and memory of the answer, not the answer.
    pub fn reset_with_bounds_cells(
        &mut self,
        q_min: i32,
        q_max: i32,
        s_min: i32,
        s_max: i32,
        max_cells: usize,
    ) {
        self.s_min = s_min;
        self.s_max = s_max;
        self.nodes.clear();
        self.nodes.push(IntervalNode::new_internal(q_min, q_max));
        self.hsps.clear();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:166-170
        // Re-initialisation point: the side indexes are rebuilt for the new ranges.
        self.x
            .set_bounds_cells(q_min, q_max, s_min, s_max, max_cells.max(1));
    }

    /// Allocate a new internal node
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:57-112 s_IntervalNodeInit
    ///
    /// CRITICAL: Right child uses midpt + 1, not midpt!
    /// The two halves are disjoint, either may hold a single position, and together they span
    /// the whole parent range.
    fn alloc_internal_node(&mut self, parent_idx: usize, which_half: IntervalDirection) -> usize {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:89-95
        let parent = &self.nodes[parent_idx];
        let midpt = ((parent.leftend as i64 + parent.rightend as i64) / 2) as i32;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:102-109
        // The right child starts at one past `midpt` (the +1 offset is critical).
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:269-299 s_IntervalRootNodeInit
    fn alloc_subject_root_node(&mut self) -> usize {
        let node = IntervalNode::new_internal(self.s_min, self.s_max);
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// Allocate a new leaf node
    fn alloc_leaf_node(&mut self, hsp: TreeHsp, query_context_offset: i32) -> usize {
        self.hsps.push(hsp);
        let node = IntervalNode::new_leaf(self.hsps.len() as u32, query_context_offset);
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:952-993
    /// `node->hsp`: the HSP of a leaf, `None` (NULL) for an internal node.
    #[inline]
    fn leaf_hsp(&self, node: &IntervalNode) -> Option<&TreeHsp> {
        node.hsp
            .checked_sub(1)
            .map(|index| &self.hsps[index as usize])
    }

    /// Check whether an incoming HSP has the same start or end coordinate as an HSP stored
    /// in the tree, and say which of the two survives.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:241-301 s_HSPsHaveCommonEndpoint
    ///
    /// Returns: None if no common endpoint, Some(EndpointResult) indicating which HSP to keep
    fn hsps_have_common_endpoint(
        in_hsp: &TreeHsp,
        in_q_start: i32,
        tree_hsp: &TreeHsp,
        tree_q_start: i32,
        which_end: IntervalDirection,
    ) -> Option<EndpointResult> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:250-254
        // check if alignments are from different query sequences or query strands
        if in_q_start != tree_q_start {
            return None;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:256-259
        // HSPs on different subject strands never conflict
        if in_hsp.subject_frame_sign.signum() != tree_hsp.subject_frame_sign.signum() {
            return None;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:261-268
        let match_found = match which_end {
            IntervalDirection::Left => {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:261-263
                in_hsp.query_offset == tree_hsp.query_offset
                    && in_hsp.subject_offset == tree_hsp.subject_offset
            }
            IntervalDirection::Right => {
                in_hsp.query_end == tree_hsp.query_end && in_hsp.subject_end == tree_hsp.subject_end
            }
            IntervalDirection::Neither => false,
        };

        if match_found {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:273-278
            // the HSP with the higher score wins
            if in_hsp.score > tree_hsp.score {
                return Some(EndpointResult::KeepInput);
            }
            if in_hsp.score < tree_hsp.score {
                return Some(EndpointResult::KeepTree);
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:280-286
            // score tie: the smaller span wins (query span first, then subject span)
            let in_q_length = in_hsp.query_end - in_hsp.query_offset;
            let tree_q_length = tree_hsp.query_end - tree_hsp.query_offset;
            if in_q_length > tree_q_length {
                return Some(EndpointResult::KeepTree);
            }
            if in_q_length < tree_q_length {
                return Some(EndpointResult::KeepInput);
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:288-293
            let in_s_length = in_hsp.subject_end - in_hsp.subject_offset;
            let tree_s_length = tree_hsp.subject_end - tree_hsp.subject_offset;
            if in_s_length > tree_s_length {
                return Some(EndpointResult::KeepTree);
            }
            if in_s_length < tree_s_length {
                return Some(EndpointResult::KeepInput);
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:295-297
            // complete tie: the HSP already stored stays
            return Some(EndpointResult::KeepTree);
        }

        None
    }

    /// Search a subtree for an HSP whose endpoint matches the input HSP's. The subtree is keyed
    /// by subject offset and hangs off a query-keyed node as that node's midpoint list.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:318-411 s_MidpointTreeHasHSPEndpoint
    ///
    /// Returns true when the input HSP must not be inserted because a better HSP with the same
    /// endpoint is already stored.
    ///
    /// Side effect: Removes worse HSPs from the tree.
    fn midpoint_tree_has_hsp_endpoint(
        &mut self,
        root_index: usize,
        in_hsp: &TreeHsp,
        in_q_start: i32,
        which_end: IntervalDirection,
    ) -> bool {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:331-334
        let target_offset = match which_end {
            IntervalDirection::Left => in_hsp.subject_offset,
            IntervalDirection::Right => in_hsp.subject_end,
            IntervalDirection::Neither => return false,
        };

        let mut root_idx = root_index;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:338
        // Walk down the interval tree until reaching a leaf or match.
        loop {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:343-366
            // Compare against every HSP in this node's midpoint list first.
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:350-365
            let mut list_idx = root_idx;
            let mut tmp_index = self.nodes[root_idx].midptr;

            while tmp_index != 0 {
                let (next_idx, result) = {
                    let tmp_node = &self.nodes[tmp_index as usize];
                    let tree_hsp = match self.leaf_hsp(tmp_node) {
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
                        // Unlink the worse HSP from the list by updating `midptr` to `next_idx`.
                        self.nodes[list_idx].midptr = next_idx;
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:361-362
                        // The list unlink is the line above; the flag only tells the grid index
                        // that the leaf is gone.
                        self.nodes[tmp_index as usize].x_dead = true;
                    }
                    None => {}
                }
                list_idx = tmp_index as usize;
                tmp_index = next_idx;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:368-376
            // Continue into the child whose range holds in_hsp's endpoint
            let node = &self.nodes[root_idx];
            let midpt = ((node.leftend as i64 + node.rightend as i64) / 2) as i32;

            let next_child_idx = if target_offset < midpt {
                node.leftptr
            } else if target_offset > midpt {
                node.rightptr
            } else {
                0
            };

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:382-383
            if next_child_idx == 0 {
                return false;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:385-406
            let next_node = &self.nodes[next_child_idx as usize];
            if let Some(tree_hsp) = self.leaf_hsp(next_node) {
                // Reached a leaf; compare in_hsp with the alignment in the leaf
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
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:398-404
                        // the stored HSP lost: detach it from the tree
                        if target_offset < midpt {
                            self.nodes[root_idx].leftptr = 0;
                        } else {
                            self.nodes[root_idx].rightptr = 0;
                        }
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:398-404
                        // The leaf unlink is the code above; the flag only tells the grid index
                        // that the leaf is gone.
                        self.nodes[next_child_idx as usize].x_dead = true;
                        return false;
                    }
                    None => {}
                }
                return false;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:408
            root_idx = next_child_idx as usize;
        }
    }

    /// Look for stored HSPs whose endpoint matches the input HSP's. Stored HSPs that rank
    /// below the input are removed as the search goes.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:428-506 s_IntervalTreeHasHSPEndpoint
    ///
    /// Returns true when the input HSP must not be inserted because a better HSP with the same
    /// endpoint is already stored.
    fn interval_tree_has_hsp_endpoint(
        &mut self,
        in_hsp: &TreeHsp,
        in_q_start: i32,
        which_end: IntervalDirection,
    ) -> bool {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:440-443
        let target_offset = match which_end {
            IntervalDirection::Left => in_q_start + in_hsp.query_offset,
            IntervalDirection::Right => in_q_start + in_hsp.query_end,
            IntervalDirection::Neither => return false,
        };

        let mut root_idx = 0usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:447
        // Walk down the interval tree to test containment against midpoint lists.
        loop {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:455-461
            // Compare against every HSP in this node's midpoint tree first.
            let midptr = self.nodes[root_idx].midptr;
            if midptr != 0
                && self.midpoint_tree_has_hsp_endpoint(
                    midptr as usize,
                    in_hsp,
                    in_q_start,
                    which_end,
                )
            {
                return true;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:466-471
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

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:477-478
            if next_child_idx == 0 {
                return false;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:480-501
            let next_node = &self.nodes[next_child_idx as usize];
            if let Some(tree_hsp) = self.leaf_hsp(next_node) {
                // Reached a leaf
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
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:493-499
                        // the stored HSP lost: detach it from the tree
                        if target_offset < midpt {
                            self.nodes[root_idx].leftptr = 0;
                        } else {
                            self.nodes[root_idx].rightptr = 0;
                        }
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:493-499
                        // The leaf unlink is the code above; the flag only tells the grid index
                        // that the leaf is gone.
                        self.nodes[next_child_idx as usize].x_dead = true;
                        return false;
                    }
                    None => {}
                }
                return false;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:503
            root_idx = next_child_idx as usize;
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:440-443,477-478
    /// `s_IntervalTreeHasHSPEndpoint` only acts on tree HSPs with the same start (or end)
    /// point as the input. If the point set has no such point the walk finds nothing, unlinks
    /// nothing and returns FALSE, so it is skipped. Otherwise the unchanged walk runs.
    /// No NCBI counterpart: `interval_tree_has_hsp_endpoint`,
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
        false
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:816-831,889-899
    /// Tests the leaves listed in the grid cell of the input's start point with the same
    /// predicate (`is_hsp_contained`). Any tree HSP that contains the input contains its
    /// start point, so it is listed in that cell. Which container is returned can differ from
    /// the walk; whether there is one does not.
    /// No NCBI counterpart: a linked tree HSP that contains `hsp`,
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
            let tree_hsp = self.leaf_hsp(node)?;
            Self::is_hsp_contained(
                hsp,
                query_start,
                tree_hsp,
                node.leftptr,
                min_diag_separation,
            )
            .then_some(*tree_hsp)
        };
        if !x.used_cells.is_empty() {
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:510-795 BlastIntervalTreeAddHSP
    ///
    /// This is the main entry point for adding HSPs to the interval tree.
    /// For eQueryAndSubject mode, it first checks for common endpoints and
    /// removes worse HSPs before adding.
    pub fn add_hsp(&mut self, hsp: TreeHsp, query_context_offset: i32, index_method: IndexMethod) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:537-550
        let mut query_start = query_context_offset;
        let (region_start, region_end) =
            if index_method == IndexMethod::StrandIndifferentQueryOnly && hsp.query_frame < 0 {
                let region_end = query_start - hsp.query_offset;
                let region_start = query_start - hsp.query_end;
                query_start = query_start - hsp.query_length - 1;
                (region_start, region_end)
            } else {
                (query_start + hsp.query_offset, query_start + hsp.query_end)
            };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:558-561
        // NCBI looks for common end points only in this mode. The side indexes are built for it
        // and are dropped (`usable = false`) for any other mode.
        // No NCBI counterpart: the side indexes describe trees built with eQueryAndSubject
        // only.
        if index_method != IndexMethod::QueryAndSubject {
            self.x.usable = false;
        }
        let x_on = self.x.usable;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:558-585
        // For eQueryAndSubject, check for common endpoints before adding
        if index_method == IndexMethod::QueryAndSubject {
            // Check left endpoint
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:577-584
            // `x_has_hsp_endpoint` (both calls below) runs the ported
            // `interval_tree_has_hsp_endpoint` unless the point set shows it cannot find
            // anything.
            if self.x_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Left, x_on) {
                return; // Better HSP with same endpoint already exists
            }
            // Check right endpoint
            if self.x_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Right, x_on) {
                return; // Better HSP with same endpoint already exists
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:591-599
        // Encapsulate the input HSP in an SIntervalNode (leaf node)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:591-599
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
                    if let Some(index) = node.hsp.checked_sub(1) {
                        let leaf_hsp = &self.hsps[index as usize];
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:601-793
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:603
        // Descend the interval tree until inserting or splitting the target interval.
        loop {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:608-609
            let middle = {
                let node = &self.nodes[root_index];
                ((node.leftend as i64 + node.rightend as i64) / 2) as i32
            };

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:611-633
            if current_region_end < middle {
                // New interval belongs in left subtree
                let leftptr = self.nodes[root_index].leftptr;

                if leftptr == 0 {
                    // No left child - attach new leaf directly
                    self.nodes[root_index].leftptr = new_index as i32;
                    return;
                }

                // Check if existing node is a leaf or internal
                let old_node_is_leaf = self.nodes[leftptr as usize].is_leaf();
                if !old_node_is_leaf {
                    // Descend to internal node
                    root_index = leftptr as usize;
                    continue;
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:703-793
                // Two leaves in same subtree - need to split
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:703-746
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:635-657
            else if current_region_start > middle {
                // New interval belongs in right subtree
                let rightptr = self.nodes[root_index].rightptr;

                if rightptr == 0 {
                    // No right child - attach new leaf directly
                    self.nodes[root_index].rightptr = new_index as i32;
                    return;
                }

                // Check if existing node is a leaf or internal
                let old_node_is_leaf = self.nodes[rightptr as usize].is_leaf();
                if !old_node_is_leaf {
                    // Descend to internal node
                    root_index = rightptr as usize;
                    continue;
                }

                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:703-793
                // Two leaves in same subtree - need to split
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:703-746
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:659-701
            else {
                // The new interval crosses the center of the node
                if index_subject_range
                    || index_method == IndexMethod::QueryOnly
                    || index_method == IndexMethod::StrandIndifferentQueryOnly
                {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:668-675
                    // If indexing subject offsets already, or only indexing query,
                    // prepend the new node to the list of "midpoint" nodes
                    self.nodes[new_index].midptr = self.nodes[root_index].midptr;
                    self.nodes[root_index].midptr = new_index as i32;
                    return;
                } else {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:677-700
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:703-793
    fn split_and_insert(
        &mut self,
        parent_root_index: usize,
        old_leaf_index: usize,
        which_half: IntervalDirection,
        index_subject_range: bool,
        index_method: IndexMethod,
    ) -> usize {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:709-712
        // Allocate new internal node
        let mid_index = self.alloc_internal_node(parent_root_index, which_half);

        // Get the old HSP data before we modify anything
        let old_hsp = *self.leaf_hsp(&self.nodes[old_leaf_index]).unwrap();
        let old_q_start = self.nodes[old_leaf_index].leftptr; // leftptr stores query_start for leaves

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:717-720
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:726-743
        let (old_region_start, old_region_end) = if index_subject_range {
            (old_hsp.subject_offset, old_hsp.subject_end)
        } else if index_method == IndexMethod::StrandIndifferentQueryOnly && old_hsp.query_frame < 0
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:747-748
        let middle = {
            let node = &self.nodes[mid_index];
            ((node.leftend as i64 + node.rightend as i64) / 2) as i32
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:749-791
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
                || index_method == IndexMethod::StrandIndifferentQueryOnly
            {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:771
                self.nodes[mid_index].midptr = old_leaf_index as i32;
            } else {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:773-790
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

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:746
        // Return the new internal node so the caller can continue the loop.
        mid_index
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:102-109,603-793
    /// `add_hsp` (eQueryAndSubject) of an HSP of TBLASTN's one-hit word finder or of RPS-BLAST's
    /// gapped stage, or `false`, the tree unchanged, where NCBI's insertion would split a node
    /// `[p, p]` to the right into the inverted node `[p + 1, p]` (whose middle is `p` again:
    /// `BlastIntervalTreeAddHSP` repeats that split without end and never returns), or split a
    /// node into another inverted one (a root `[0, R]` with `R` below 0). Only ranges that leave
    /// the tree's bounds or are out of order lead there (a range in order and inside `[p, p]`
    /// straddles its middle `p`): in the TBLASTN word finder
    /// (`tblastn::search_gapped::x_add_one_hit_hsp`) one that starts past the tree's `s_max`, and
    /// one whose start is one past its end (a gapped alignment of score 0 whose right extension did
    /// not run, which a cutoff of 0 keeps: blast_gapalign.c:4290-4295, blast_stat.c:5236-5266),
    /// which never straddles a middle; in RPS-BLAST a query range that starts past the root's end,
    /// where an alias file's `LENGTH` and `NSEQ` make the concatenated length that the root covers
    /// shorter than the profiles (`rps_search::prelim`). An HSP whose query and subject ranges are
    /// in order and lie inside the tree's bounds (`x_insertion_ends`) goes to `add_hsp` at once.
    /// For any other HSP the insertion is first walked without changing the tree
    /// (`x_insert_inverts`); only if that walk makes an inverted node are NCBI's end-point checks,
    /// which run before the insertion and may drop the HSP or a leaf on its way (end points only
    /// empty a slot), run on a copy of the tree and the walk repeated there. `add_hsp` and the
    /// functions it calls keep their code.
    ///
    /// The test of the bounds is inlined into the caller (RPS-BLAST's gapped stage adds every HSP
    /// through here); the walk is out of line.
    #[inline]
    pub fn try_add_hsp(&mut self, hsp: TreeHsp, query_context_offset: i32) -> bool {
        if self.x_insertion_ends(&hsp, query_context_offset) {
            self.add_hsp(hsp, query_context_offset, IndexMethod::QueryAndSubject);
            return true;
        }
        self.x_try_add_hsp_outside(hsp, query_context_offset)
    }

    /// `try_add_hsp` of an HSP that leaves the tree's bounds or is out of order.
    #[cold]
    #[inline(never)]
    fn x_try_add_hsp_outside(&mut self, hsp: TreeHsp, query_context_offset: i32) -> bool {
        if self.x_insert_inverts(&hsp, query_context_offset) {
            let mut copy = self.clone();
            let x_on = copy.x.usable;
            let dropped =
                copy.x_has_hsp_endpoint(&hsp, query_context_offset, IntervalDirection::Left, x_on)
                    || copy.x_has_hsp_endpoint(
                        &hsp,
                        query_context_offset,
                        IntervalDirection::Right,
                        x_on,
                    );
            if !dropped && copy.x_insert_inverts(&hsp, query_context_offset) {
                return false;
            }
        }
        self.add_hsp(hsp, query_context_offset, IndexMethod::QueryAndSubject);
        true
    }

    /// Whether NCBI's insertion loop (eQueryAndSubject) surely ends for `hsp`, whatever the leaves
    /// already in the tree: its query and subject ranges are in order (start <= end) and inside the
    /// tree's bounds, so each node the insertion enters or makes holds the range it compares (a
    /// split goes to the half that holds it) and every split makes a smaller node; a node `[p, p]`
    /// holding the range sees it straddle its middle.
    #[inline(always)]
    fn x_insertion_ends(&self, hsp: &TreeHsp, query_context_offset: i32) -> bool {
        let root = &self.nodes[0];
        hsp.query_offset <= hsp.query_end
            && hsp.subject_offset <= hsp.subject_end
            && query_context_offset + hsp.query_offset >= root.leftend
            && query_context_offset + hsp.query_end <= root.rightend
            && hsp.subject_offset >= self.s_min
            && hsp.subject_end <= self.s_max
    }

    /// Whether the insertion loop of `add_hsp_internal` (eQueryAndSubject) would make an inverted
    /// node for `hsp`, walked without changing the tree: the real nodes down to the first leaf
    /// the HSP meets, then the nodes that the splits with that leaf would make (each holds only
    /// that leaf, placed as `split_and_insert` places it, and the halves of the subject root made
    /// for it when both straddle a query node's middle). A split to the left of `[p, p]` would
    /// repeat `[p, p]` without end too (a range that ends before `p` there: one whose start is one
    /// past its end, or one before the tree's lower bound, which no caller gives); it counts as
    /// such a node.
    fn x_insert_inverts(&self, hsp: &TreeHsp, query_context_offset: i32) -> bool {
        #[derive(PartialEq)]
        enum Side {
            Left,
            Right,
            Mid,
        }
        let middle = |l: i32, r: i32| ((i64::from(l) + i64::from(r)) / 2) as i32;
        let side = |start: i32, end: i32, m: i32| {
            if end < m {
                Side::Left
            } else if start > m {
                Side::Right
            } else {
                Side::Mid
            }
        };
        let (mut start, mut end) = (
            query_context_offset + hsp.query_offset,
            query_context_offset + hsp.query_end,
        );
        let mut subject_level = false;
        // The real nodes.
        let mut root = 0usize;
        let (mut l, mut r, old) = loop {
            let node = &self.nodes[root];
            let m = middle(node.leftend, node.rightend);
            let (child, bounds) = match side(start, end, m) {
                Side::Left => (node.leftptr, (node.leftend, m)),
                Side::Right => (node.rightptr, (m + 1, node.rightend)),
                Side::Mid => {
                    if subject_level || node.midptr == 0 {
                        // A midpoint list, or a new subject tree holding only this HSP.
                        return false;
                    }
                    subject_level = true;
                    (start, end) = (hsp.subject_offset, hsp.subject_end);
                    root = node.midptr as usize;
                    continue;
                }
            };
            if child == 0 {
                return false;
            }
            if !self.nodes[child as usize].is_leaf() {
                root = child as usize;
                continue;
            }
            if (bounds.0, bounds.1) == (node.leftend, node.rightend) {
                return true;
            }
            break (bounds.0, bounds.1, child as usize);
        };
        // The nodes of the splits with the leaf `old`.
        let old_node = &self.nodes[old];
        let Some(old_hsp) = self.leaf_hsp(old_node) else {
            return false;
        };
        loop {
            if l > r {
                return true;
            }
            let m = middle(l, r);
            let (old_start, old_end) = if subject_level {
                (old_hsp.subject_offset, old_hsp.subject_end)
            } else {
                (
                    old_node.leftptr + old_hsp.query_offset,
                    old_node.leftptr + old_hsp.query_end,
                )
            };
            let placed = side(old_start, old_end, m);
            let next = match side(start, end, m) {
                Side::Left if placed == Side::Left => (l, m),
                Side::Right if placed == Side::Right => (m + 1, r),
                Side::Mid if !subject_level && placed == Side::Mid => {
                    // Both in the new node's subject tree: its root `[s_min, s_max]` holds the old
                    // leaf alone (blast_itree.c:773-790) and the new HSP follows it there
                    // (:677-700); a split of that root makes one of its halves, and the next split
                    // is compared with the root, not with the query node.
                    subject_level = true;
                    (start, end) = (hsp.subject_offset, hsp.subject_end);
                    (l, r) = (self.s_min, self.s_max);
                    let sm = middle(l, r);
                    let old_side = side(old_hsp.subject_offset, old_hsp.subject_end, sm);
                    match side(start, end, sm) {
                        Side::Left if old_side == Side::Left => (l, sm),
                        Side::Right if old_side == Side::Right => (sm + 1, r),
                        _ => return false,
                    }
                }
                _ => return false,
            };
            if next == (l, r) {
                return true;
            }
            (l, r) = next;
        }
    }

    /// Return the tree HSP that envelops the input HSP, if any.
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:930-995 BlastIntervalTreeContainsHSP
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:930-995
    /// See the citation on `contains_hsp`; this method preserves that traversal
    /// and exposes the matched tree HSP for Phase 6 parity diagnostics.
    pub fn containing_hsp(
        &self,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> Option<TreeHsp> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:931-935,952,992-994
        // `containing_hsp_tree` below ports this function. The grid answers the same yes/no
        // question; the walk answers where the index is off (tracing on) or was dropped.
        if self.x.usable {
            return self.x_find_container(hsp, query_context_offset, min_diag_separation);
        }
        self.containing_hsp_tree(hsp, query_context_offset, min_diag_separation)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:952-988
    /// The port of the walk, which runs where the side indexes are off or were dropped.
    /// The tree walk of `containing_hsp`.
    fn containing_hsp_tree(
        &self,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> Option<TreeHsp> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:942
        if self.nodes.is_empty() {
            return None;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:944-948
        let query_start = query_context_offset;
        let region_start = query_start + hsp.query_offset;
        let region_end = query_start + hsp.query_end;

        let mut node_idx = 0usize;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:950-994
        loop {
            let node = &self.nodes[node_idx];

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:990-994
            // If this is a leaf node, check containment directly
            if let Some(tree_hsp) = self.leaf_hsp(node) {
                return Self::is_hsp_contained(
                    hsp,
                    query_start,
                    tree_hsp,
                    node.leftptr, // leftptr stores query_context_offset for leaves
                    min_diag_separation,
                )
                .then_some(*tree_hsp);
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:960-967
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

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:973-985
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

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:987
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
        if let Some(tree_hsp) = self.leaf_hsp(node) {
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
        if midptr > 0
            && self.midpoint_tree_contains_hsp(
                midptr as usize,
                hsp,
                query_context_offset,
                min_diag_separation,
            )
        {
            return true;
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:864-927 s_MidpointTreeContainsHSP
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
            if let Some(tree_hsp) = self.leaf_hsp(node) {
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
                if let Some(tree_hsp) = self.leaf_hsp(tmp_node) {
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:809-847 s_HSPIsContained
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:819
        // Check if alignments are from different query sequences or query strands
        if in_q_start != tree_q_start {
            return false;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:822-831
        // All conditions must be true:
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:833-834
            if min_diag_separation == 0 {
                return true;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:836-843
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign_priv.h:120-121
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign_priv.h:123-124
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
#[path = "interval_tree_tests.rs"]
mod tests;
