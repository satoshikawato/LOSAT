//! NCBI BLAST HSP culling (`-culling_limit`): the interval tree of `c++/src/algo/blast/core/hspfilter_culling.c`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c
//!
//! The interval tree of NCBI's culling writer and pipe: every HSP of a query context is
//! saved in the tree of its context unless `culling_max` HSPs already dominate it, and a
//! saved HSP takes one merit from every HSP that it dominates. NCBI's linked lists and
//! tree nodes are an arena here (indices instead of pointers), so the identity test
//! `r != y` of `s_ProcessHSPList` compares indices; the lists keep NCBI's order (a new
//! HSP goes to the front, `s_AddHSPtoList`).
//!
//! The algorithm is described in: Berman P, Zhang Z, Wolf YI, Koonin EV, Miller W.
//! Winnowing sequences from a database search. J Comput Biol. 2000 Feb-Apr;7(1-2):293-302.
//!
//! NCBI culls twice when `-culling_limit` is given: a writer in the preliminary stage and a
//! pipe after the traceback. This module is the part both share and holds no program type:
//! the tree never touches the payload `H` of an HSP; the keys it tests (the plus-strand
//! query range, the score and the subject offset; `LinkedHSP` and `hsp->score` /
//! `hsp->subject.offset` of `s_DominateTest`) are copied into the entry by the caller, who
//! also decides the context, the tree length and the strand fold
//! (`s_BlastHSPCullingRun`, program specific).

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, VecDeque};

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:52-60
/// `hsp` is the index of the payload in the forest's arena; `next` is the order of the
/// node's list (`CTreeNode::hsplist`). `score` and `subject_offset` are those of
/// `hsp->score` and `hsp->subject.offset`, which `s_DominateTest` reads. `merit` is an
/// `i64`, not NCBI's `Int4`: a culling limit near `INT4_MAX` (plus 3 in the preliminary
/// stage) overflows an `Int4`, which is undefined behaviour in C and wraps in the shipped
/// binary; LOSAT does not reproduce that overflow (see `prelim_max`).
#[derive(Clone, Copy)]
struct LinkedHsp {
    hsp: usize,
    sid: u32,
    begin: i32,
    end: i32,
    score: i32,
    subject_offset: i64,
    merit: i64,
}

/// What `s_BlastHSPCullingRun` makes of an HSP for the tree: its range on the plus strand of
/// the query (`LinkedHSP.begin` / `end`, already folded for a minus context) and the HSP's
/// `score` and `subject.offset`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CullKey {
    pub(crate) begin: i32,
    pub(crate) end: i32,
    pub(crate) score: i32,
    pub(crate) subject_offset: i64,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:200-207
/// `left` and `right` are node indices; a freed node is unlinked from its parent.
struct CTreeNode {
    begin: i32,
    end: i32,
    left: Option<usize>,
    right: Option<usize>,
    hsplist: VecDeque<usize>,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:435
const K_NUM_HSP_TO_FORK: i32 = 20;

/// The HSPs that a forest keeps for one query: one list per subject, in the order of the
/// subject's first HSP in the rip, with the HSPs of a list in their final order.
pub(crate) type CulledQuery<H> = Vec<(u32, Vec<H>)>;

/// NCBI's culling writer (or pipe) over the HSPs of one search: `s_BlastHSPCullingInit`
/// (`new`), `s_BlastHSPCullingRun` (`save`, once per HSP) and `s_BlastHSPCullingFinal`
/// (`finalize`). `trees` is `c_tree`: the root node of each tree id (NCBI's `cid`) with a tree.
///
/// NCBI frees an HSP, its `LinkedHSP` and an emptied tree node at once (`s_HSPFree`,
/// `s_CTreeNodeFree`). The arenas do the same with free lists: the slot of a freed HSP, the
/// entry of a freed `LinkedHSP` and a freed node are reused by the next one that is made, so
/// the arenas hold at most as many items as were alive together, not one for every HSP
/// offered. The index values are never compared except `r != y` between live entries, and
/// the lists keep their order, so the reuse changes no output.
pub(crate) struct CullingForest<H> {
    culling_max: i64,
    hsps: Vec<Option<H>>,
    entries: Vec<LinkedHsp>,
    nodes: Vec<CTreeNode>,
    trees: BTreeMap<usize, usize>,
    free_hsps: Vec<usize>,
    free_entries: Vec<usize>,
    free_nodes: Vec<usize>,
}

impl<H> CullingForest<H> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:487-493
    pub(crate) fn new(culling_max: i64) -> Self {
        Self {
            culling_max,
            hsps: Vec::new(),
            entries: Vec::new(),
            nodes: Vec::new(),
            trees: BTreeMap::new(),
            free_hsps: Vec::new(),
            free_entries: Vec::new(),
            free_nodes: Vec::new(),
        }
    }

    /// `s_HSPFree` of the linked HSP `entry`: its payload is dropped and its slot and entry
    /// go to the free lists.
    fn free_entry(&mut self, entry: usize) {
        let slot = self.entries[entry].hsp;
        self.hsps[slot] = None;
        self.free_hsps.push(slot);
        self.free_entries.push(entry);
    }

    /// `s_CTreeNodeFree` of an emptied node (no list, no children, already unlinked from its
    /// parent by the caller): the index goes to the free list.
    fn free_node(&mut self, node: usize) {
        debug_assert!(self.nodes[node].hsplist.is_empty());
        self.nodes[node].hsplist = VecDeque::new();
        self.free_nodes.push(node);
    }

    /// The `culling_max` of the preliminary writer: three more than the option when it is
    /// above 1. NCBI's sum is an `Int4`, which overflows for a limit of 2147483645 or more
    /// (the shipped binary wraps: 2147483645 culls nothing, 2147483646 and 2147483647 act as
    /// a limit of 1). LOSAT adds in `i64` (Owner 2026-10-22, NCBI parity DECISIONS 17: no
    /// reproduction of NCBI's undefined behaviour), so every limit above 1 is three more.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:338-340
    pub(crate) fn prelim_max(limit: i32) -> i64 {
        let limit = i64::from(limit);
        if limit > 1 {
            limit + 3
        } else {
            limit
        }
    }

    /// One HSP of `s_BlastHSPCullingRun`: the HSP of subject `sid`, with the keys `key`, goes
    /// into the tree `cid` (made with the length `qlen` on first use) unless `culling_max`
    /// HSPs already dominate it; an HSP that is not saved is freed.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:612-636
    pub(crate) fn save(&mut self, cid: usize, qlen: i32, sid: u32, key: CullKey, hsp: H) {
        let index = match self.free_hsps.pop() {
            Some(index) => {
                self.hsps[index] = Some(hsp);
                index
            }
            None => {
                self.hsps.push(Some(hsp));
                self.hsps.len() - 1
            }
        };
        let mut a = LinkedHsp {
            hsp: index,
            sid,
            begin: key.begin,
            end: key.end,
            score: key.score,
            subject_offset: key.subject_offset,
            merit: self.culling_max,
        };
        let root = match self.trees.get(&cid) {
            Some(&root) => root,
            None => {
                let root = self.ctree_new(qlen);
                self.trees.insert(cid, root);
                root
            }
        };
        if !self.save_hsp(root, &mut a) {
            self.hsps[index] = None;
            self.free_hsps.push(index);
        }
    }

    /// The HSPs that the trees keep, as `s_BlastHSPCullingFinal` puts them into the hit list
    /// of each query: the trees in order of their ids, each tree ripped into one list, a new
    /// HSP list for each subject in the order of its first HSP in the query's hit list (an
    /// index from subject to list stands for NCBI's linear search, which finds the same
    /// list), and after each tree every HSP list of the query sorted by `cmp`
    /// (`Blast_HSPListSortByScore`, glibc's stable `qsort`). The hit list is not cut to
    /// `prelim_hitlist_size`. `query_of` maps a tree id to its query
    /// (`Blast_GetQueryIndexFromContext`); the result has one entry for each query with a
    /// tree, in query order.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:516-581
    pub(crate) fn finalize(
        mut self,
        query_of: impl Fn(usize) -> usize,
        mut cmp: impl FnMut(&H, &H) -> Ordering,
    ) -> Vec<(usize, CulledQuery<H>)> {
        struct Group<H> {
            lists: CulledQuery<H>,
            position: HashMap<u32, usize>,
        }
        let mut groups: BTreeMap<usize, Group<H>> = BTreeMap::new();
        let roots: Vec<(usize, usize)> = self.trees.iter().map(|(&c, &r)| (c, r)).collect();
        for (cid, root) in roots {
            let group = groups.entry(query_of(cid)).or_insert_with(|| Group {
                lists: Vec::new(),
                position: HashMap::new(),
            });
            let mut ripped = Vec::new();
            self.rip_hsp_off_ctree(Some(root), &mut ripped);
            for entry in ripped {
                let LinkedHsp { hsp, sid, .. } = self.entries[entry];
                let hsp = self.hsps[hsp].take().expect("an HSP is ripped once");
                match group.position.get(&sid) {
                    Some(&at) => group.lists[at].1.push(hsp),
                    None => {
                        group.position.insert(sid, group.lists.len());
                        group.lists.push((sid, vec![hsp]));
                    }
                }
            }
            for (_, list) in &mut group.lists {
                list.sort_by(&mut cmp);
            }
        }
        groups
            .into_iter()
            .map(|(query, group)| (query, group.lists))
            .collect()
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120
    fn dominate_test(&self, p: &LinkedHsp, y: &LinkedHsp) -> bool {
        let (b1, b2) = (i64::from(p.begin), i64::from(y.begin));
        let (e1, e2) = (i64::from(p.end), i64::from(y.end));
        let (s1, s2) = (i64::from(p.score), i64::from(y.score));
        let l1 = e1 - b1;
        let l2 = e2 - b2;
        let overlap = e1.min(e2) - b1.max(b2);
        if 2 * overlap < l2 {
            return false;
        }
        let d = 4 * s1 * l1 + 2 * s1 * l2 - 2 * s2 * l1 - 4 * s2 * l2;
        if (s1 == s2 && b1 == b2 && l1 == l2) || d == 0 {
            if s1 != s2 {
                return s1 > s2;
            }
            if p.sid != y.sid {
                return p.sid < y.sid;
            }
            if p.subject_offset > y.subject_offset {
                return false;
            }
            return true;
        }
        d >= 0
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:123-133
    /// `merit` is an `i64` (see `LinkedHsp`).
    fn full_pass(&self, node: usize, y: &mut LinkedHsp) -> bool {
        for &p in &self.nodes[node].hsplist {
            if self.dominate_test(&self.entries[p], y) {
                y.merit -= 1;
                if y.merit <= 0 {
                    return false;
                }
            }
        }
        true
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:136-163
    /// `y` is the index of the new HSP's entry (`x` of `s_SaveHSP`).
    fn process_hsp_list(&mut self, node: usize, y: usize) -> i32 {
        // The list is unlinked in place (`retain_mut` visits from the front, as NCBI's `p`).
        let mut list = std::mem::take(&mut self.nodes[node].hsplist);
        let mut num = 0;
        let y_entry = self.entries[y];
        list.retain_mut(|r| {
            let r = *r;
            num += 1;
            if r != y && self.dominate_test(&y_entry, &self.entries[r]) {
                let merit = self.entries[r].merit - 1;
                self.entries[r].merit = merit;
                if merit <= 0 {
                    num -= 1;
                    self.free_entry(r);
                    return false;
                }
            }
            true
        });
        self.nodes[node].hsplist = list;
        num
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:166-189
    fn mark_down_hsp_list(&mut self, node: usize) -> i32 {
        let mut list = std::mem::take(&mut self.nodes[node].hsplist);
        let mut num = 0;
        list.retain_mut(|r| {
            let r = *r;
            num += 1;
            let merit = self.entries[r].merit - 1;
            self.entries[r].merit = merit;
            if merit <= 0 {
                num -= 1;
                self.free_entry(r);
                return false;
            }
            true
        });
        self.nodes[node].hsplist = list;
        num
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:228-247
    fn ctree_node_new(&mut self, parent: Option<usize>, left: bool) -> usize {
        let (begin, end) = match parent {
            None => (0, 0),
            Some(parent) => {
                let (begin, end) = (self.nodes[parent].begin, self.nodes[parent].end);
                let midpt = (begin + end) / 2;
                if left {
                    (begin, midpt)
                } else {
                    (midpt, end)
                }
            }
        };
        let fresh = CTreeNode {
            begin,
            end,
            left: None,
            right: None,
            hsplist: VecDeque::new(),
        };
        match self.free_nodes.pop() {
            Some(node) => {
                self.nodes[node] = fresh;
                node
            }
            None => {
                self.nodes.push(fresh);
                self.nodes.len() - 1
            }
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:376-381
    fn ctree_new(&mut self, qlen: i32) -> usize {
        let tree = self.ctree_node_new(None, true);
        self.nodes[tree].begin = 0;
        self.nodes[tree].end = qlen;
        tree
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:258-299
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:193-197
    fn fork_children(&mut self, node: usize) {
        let midpt = (self.nodes[node].begin + self.nodes[node].end) / 2;
        let mut list = std::mem::take(&mut self.nodes[node].hsplist);
        list.retain_mut(|r| {
            let r = *r;
            let LinkedHsp { begin, end, .. } = self.entries[r];
            let child = if end < midpt {
                Some(match self.nodes[node].left {
                    Some(left) => left,
                    None => {
                        let left = self.ctree_node_new(Some(node), true);
                        self.nodes[node].left = Some(left);
                        left
                    }
                })
            } else if begin > midpt {
                Some(match self.nodes[node].right {
                    Some(right) => right,
                    None => {
                        let right = self.ctree_node_new(Some(node), false);
                        self.nodes[node].right = Some(right);
                        right
                    }
                })
            } else {
                None
            };
            match child {
                Some(child) => {
                    self.nodes[child].hsplist.push_front(r);
                    false
                }
                None => true,
            }
        });
        self.nodes[node].hsplist = list;
    }

    /// Returns whether the node is freed (NCBI sets the pointer `*node` to null).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:319-330
    fn mark_down_ctree(&mut self, node: usize) -> bool {
        if let Some(left) = self.nodes[node].left {
            if self.mark_down_ctree(left) {
                self.nodes[node].left = None;
                self.free_node(left);
            }
        }
        if let Some(right) = self.nodes[node].right {
            if self.mark_down_ctree(right) {
                self.nodes[node].right = None;
                self.free_node(right);
            }
        }
        self.mark_down_hsp_list(node) <= 0
            && self.nodes[node].left.is_none()
            && self.nodes[node].right.is_none()
    }

    /// Returns whether the node is freed (NCBI sets the pointer `*node` to null).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:334-370
    fn process_ctree(&mut self, node: usize, x: usize) -> bool {
        let LinkedHsp { begin, end, .. } = self.entries[x];
        if begin <= self.nodes[node].begin && end >= self.nodes[node].end {
            return self.mark_down_ctree(node);
        }
        if self.nodes[node].left.is_none() && self.nodes[node].right.is_none() {
            return self.process_hsp_list(node, x) <= 0;
        }
        let midpt = (self.nodes[node].begin + self.nodes[node].end) / 2;
        let left = |writer: &mut Self| {
            if let Some(left) = writer.nodes[node].left {
                if writer.process_ctree(left, x) {
                    writer.nodes[node].left = None;
                    writer.free_node(left);
                }
            }
        };
        let right = |writer: &mut Self| {
            if let Some(right) = writer.nodes[node].right {
                if writer.process_ctree(right, x) {
                    writer.nodes[node].right = None;
                    writer.free_node(right);
                }
            }
        };
        if end < midpt {
            left(self);
            false
        } else if begin > midpt {
            right(self);
            false
        } else {
            left(self);
            right(self);
            self.process_hsp_list(node, x) <= 0
                && self.nodes[node].left.is_none()
                && self.nodes[node].right.is_none()
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:396-423
    fn rip_hsp_off_ctree(&mut self, tree: Option<usize>, out: &mut Vec<usize>) {
        let Some(tree) = tree else {
            return;
        };
        out.extend(std::mem::take(&mut self.nodes[tree].hsplist));
        let (left, right) = (self.nodes[tree].left, self.nodes[tree].right);
        self.rip_hsp_off_ctree(left, out);
        self.rip_hsp_off_ctree(right, out);
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:430-470
    fn save_hsp(&mut self, tree: usize, a: &mut LinkedHsp) -> bool {
        let mut current = Some(tree);
        let mut node = tree;
        while let Some(t) = current {
            if !self.full_pass(t, a) {
                return false;
            }
            let midpt = (self.nodes[t].begin + self.nodes[t].end) / 2;
            node = t;
            if a.end < midpt {
                current = self.nodes[t].left;
            } else if a.begin > midpt {
                current = self.nodes[t].right;
            } else {
                break;
            }
        }
        let x = match self.free_entries.pop() {
            Some(x) => {
                self.entries[x] = *a;
                x
            }
            None => {
                self.entries.push(*a);
                self.entries.len() - 1
            }
        };
        self.nodes[node].hsplist.push_front(x);
        if self.nodes[node].left.is_none() && self.nodes[node].right.is_none() {
            if self.process_hsp_list(node, x) >= K_NUM_HSP_TO_FORK {
                self.fork_children(node);
            }
            return true;
        }
        self.process_hsp_list(node, x);
        if let Some(left) = self.nodes[node].left {
            if self.process_ctree(left, x) {
                self.nodes[node].left = None;
                self.free_node(left);
            }
        }
        if let Some(right) = self.nodes[node].right {
            if self.process_ctree(right, x) {
                self.nodes[node].right = None;
                self.free_node(right);
            }
        }
        true
    }
}

#[cfg(test)]
#[path = "hspfilter_culling_tests.rs"]
mod tests;
