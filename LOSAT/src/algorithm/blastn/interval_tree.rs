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
#[derive(Clone, Debug)]
pub struct TreeHsp {
    pub query_offset: i32,
    pub query_end: i32,
    pub subject_offset: i32,
    pub subject_end: i32,
    pub score: i32,
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
        }
    }
}


/// Interval tree for HSP containment checking
/// NCBI reference: blast_itree.c:60-70 BlastIntervalTree
pub struct BlastIntervalTree {
    /// Array of nodes (index 0 is root for query-indexed tree)
    nodes: Vec<IntervalNode>,
    /// Minimum subject offset
    s_min: i32,
    /// Maximum subject offset
    s_max: i32,
    /// Flat list of all HSPs for linear containment check (temporary for debugging)
    hsps: Vec<TreeHsp>,
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
            nodes: Vec::with_capacity(64),
            s_min,
            s_max,
            hsps: Vec::with_capacity(1024),
        };

        // Create root node for query range
        let root = IntervalNode::new_internal(q_min, q_max);
        tree.nodes.push(root);

        tree
    }

    /// Reset the tree for reuse
    /// NCBI reference: blast_itree.c:168-193 Blast_IntervalTreeReset
    pub fn reset(&mut self) {
        // Keep only the root node, reset its children
        if !self.nodes.is_empty() {
            let leftend = self.nodes[0].leftend;
            let rightend = self.nodes[0].rightend;
            self.nodes.clear();
            self.nodes.push(IntervalNode::new_internal(leftend, rightend));
        }
        self.hsps.clear();
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
            IntervalDirection::Right => (midpt + 1, parent.rightend),  // +1 is CRITICAL
            IntervalDirection::Neither => panic!("Neither direction not valid for child allocation"),
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
                in_hsp.query_offset == tree_hsp.query_offset &&
                in_hsp.subject_offset == tree_hsp.subject_offset
            }
            IntervalDirection::Right => {
                // else {
                //     match = in_hsp->query.end == tree_hsp->query.end &&
                //             in_hsp->subject.end == tree_hsp->subject.end;
                // }
                in_hsp.query_end == tree_hsp.query_end &&
                in_hsp.subject_end == tree_hsp.subject_end
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
            let mut list_idx = root_idx;
            let mut tmp_index = self.nodes[root_idx].midptr;

            while tmp_index != 0 {
                let tmp_node = &self.nodes[tmp_index as usize];
                if let Some(ref tree_hsp) = tmp_node.hsp {
                    let tree_q_start = tmp_node.leftptr; // leftptr stores query_context_offset for leaves

                    let result = Self::hsps_have_common_endpoint(
                        in_hsp, in_q_start, tree_hsp, tree_q_start, which_end);

                    let next_idx = tmp_node.midptr;

                    // NCBI reference: blast_itree.c:359-362
                    match result {
                        Some(EndpointResult::KeepTree) => return true,
                        Some(EndpointResult::KeepInput) => {
                            // Remove worse HSP from list: list_node->midptr = tmp_index
                            self.nodes[list_idx].midptr = next_idx;
                        }
                        None => {
                            list_idx = tmp_index as usize;
                        }
                    }
                    tmp_index = next_idx;
                } else {
                    break;
                }
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
                    in_hsp, in_q_start, tree_hsp, tree_q_start, which_end);

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
                    midptr as usize, in_hsp, in_q_start, which_end) {
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
                    in_hsp, in_q_start, tree_hsp, tree_q_start, which_end);

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

    /// Add an HSP to the tree
    /// NCBI reference: blast_itree.c:510-795 BlastIntervalTreeAddHSP
    ///
    /// This is the main entry point for adding HSPs to the interval tree.
    /// For eQueryAndSubject mode, it first checks for common endpoints and
    /// removes worse HSPs before adding.
    pub fn add_hsp(&mut self, hsp: TreeHsp, query_context_offset: i32, index_method: IndexMethod) {
        // Also add to flat list for linear containment check (debug/verification)
        self.hsps.push(hsp.clone());

        // NCBI reference: blast_itree.c:537
        let query_start = query_context_offset;

        // NCBI reference: blast_itree.c:547-550
        let region_start = query_start + hsp.query_offset;
        let region_end = query_start + hsp.query_end;

        // NCBI reference: blast_itree.c:558-585
        // For eQueryAndSubject, check for common endpoints before adding
        if index_method == IndexMethod::QueryAndSubject {
            // Check left endpoint
            if self.interval_tree_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Left) {
                return;  // Better HSP with same endpoint already exists
            }
            // Check right endpoint
            if self.interval_tree_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Right) {
                return;  // Better HSP with same endpoint already exists
            }
        }

        // NCBI reference: blast_itree.c:591-599
        // Encapsulate the input HSP in an SIntervalNode (leaf node)
        let new_index = self.alloc_leaf_node(hsp.clone(), query_start);

        // Start the insertion loop
        self.add_hsp_internal(new_index, region_start, region_end, &hsp, query_start, false, index_method);
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
        query_start: i32,
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
                self.split_and_insert(
                    root_index,
                    new_index,
                    leftptr as usize,
                    IntervalDirection::Left,
                    index_subject_range,
                    index_method,
                    current_region_start,
                    current_region_end,
                    hsp,
                    query_start,
                );
                return;
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
                self.split_and_insert(
                    root_index,
                    new_index,
                    rightptr as usize,
                    IntervalDirection::Right,
                    index_subject_range,
                    index_method,
                    current_region_start,
                    current_region_end,
                    hsp,
                    query_start,
                );
                return;
            }
            // NCBI reference: blast_itree.c:659-701
            else {
                // The new interval crosses the center of the node
                if index_subject_range || index_method == IndexMethod::QueryOnly ||
                   index_method == IndexMethod::QueryOnlyStrandIndifferent {
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
    /// This creates a new internal node, reattaches the old leaf, then continues
    /// the insertion loop for the new leaf.
    ///
    /// NCBI reference: blast_itree.c:703-793
    fn split_and_insert(
        &mut self,
        parent_root_index: usize,
        new_leaf_index: usize,
        old_leaf_index: usize,
        which_half: IntervalDirection,
        index_subject_range: bool,
        index_method: IndexMethod,
        region_start: i32,
        region_end: i32,
        hsp: &TreeHsp,
        query_start: i32,
    ) {
        // NCBI reference: blast_itree.c:709-712
        // Allocate new internal node
        let mid_index = self.alloc_internal_node(parent_root_index, which_half);

        // Get the old HSP data before we modify anything
        let old_hsp = self.nodes[old_leaf_index].hsp.clone().unwrap();
        let old_q_start = self.nodes[old_leaf_index].leftptr;  // leftptr stores query_start for leaves

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

        // NCBI reference: blast_itree.c:726-744
        // Calculate old leaf's region
        let (old_region_start, old_region_end) = if index_subject_range {
            (old_hsp.subject_offset, old_hsp.subject_end)
        } else {
            (old_q_start + old_hsp.query_offset, old_q_start + old_hsp.query_end)
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
            if index_subject_range || index_method == IndexMethod::QueryOnly ||
               index_method == IndexMethod::QueryOnlyStrandIndifferent {
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

        // Now continue inserting the new leaf starting from the new internal node
        self.add_hsp_internal(
            new_leaf_index,
            region_start,
            region_end,
            hsp,
            query_start,
            index_subject_range,
            index_method,
        );
    }

    /// Check if the tree contains an HSP that envelops the input HSP
    /// NCBI reference: blast_itree.c:930-995 BlastIntervalTreeContainsHSP
    ///
    /// Uses proper interval tree traversal for O(log n) containment checks.
    pub fn contains_hsp(&self, hsp: &TreeHsp, query_context_offset: i32, min_diag_separation: i32) -> bool {
        // NCBI reference: blast_itree.c:942
        if self.nodes.is_empty() {
            return false;
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
                return Self::is_hsp_contained(
                    hsp,
                    query_start,
                    node.hsp.as_ref().unwrap(),
                    node.leftptr,  // leftptr stores query_context_offset for leaves
                    min_diag_separation,
                );
            }

            // NCBI reference: blast_itree.c:960-967
            // Check midpoint tree first (contains all HSPs straddling this node)
            if node.midptr > 0 {
                if self.midpoint_tree_contains_hsp(
                    node.midptr as usize,
                    hsp,
                    query_start,
                    min_diag_separation,
                ) {
                    return true;
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
                return false;
            }

            node_idx = next_idx as usize;
        }
    }

    /// Debug: Check containment with linear search (for verification)
    #[allow(dead_code)]
    pub fn contains_hsp_linear(&self, hsp: &TreeHsp, query_context_offset: i32, min_diag_separation: i32) -> bool {
        // Use linear search through all HSPs for correctness verification
        for tree_hsp in &self.hsps {
            if Self::is_hsp_contained(hsp, query_context_offset, tree_hsp, query_context_offset, min_diag_separation) {
                return true;
            }
        }
        false
    }

    /// Debug: Check containment with extra logging
    #[allow(dead_code)]
    pub fn contains_hsp_debug(&self, hsp: &TreeHsp, query_context_offset: i32, min_diag_separation: i32) -> (bool, usize) {
        let mut checked = 0usize;
        for tree_hsp in &self.hsps {
            checked += 1;
            if Self::is_hsp_contained(hsp, query_context_offset, tree_hsp, query_context_offset, min_diag_separation) {
                return (true, checked);
            }
        }
        (false, checked)
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

    /// Check midpoint subtree (subject-indexed)
    /// NCBI reference: blast_itree.c:864-927 s_MidpointTreeContainsHSP
    fn midpoint_tree_contains_hsp(
        &self,
        root_index: usize,
        hsp: &TreeHsp,
        query_context_offset: i32,
        min_diag_separation: i32,
    ) -> bool {
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
                );
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
                        return true;
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
                return false;
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
            && in_hsp.subject_frame_sign == tree_hsp.subject_frame_sign
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
}
