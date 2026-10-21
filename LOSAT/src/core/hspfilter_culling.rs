//! NCBI BLAST HSP culling (`-culling_limit`): the interval tree of `hspfilter_culling.c`.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c
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

/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:52-60
/// ```c
/// typedef struct LinkedHSP {
///     BlastHSP * hsp;
///     Int4 cid;    /* context id for hsp */
///     Int4 sid;    /* OID for hsp*/
///     Int4 begin;  /* query offset in plus strand */
///     Int4 end;    /* query end in plus strand */
///     Int4 merit;  /* how many other hsps in the tree dominates me? */
///     struct LinkedHSP *next;
/// } LinkedHSP;
/// ```
/// `hsp` is the index of the payload in the forest's arena; `next` is the order of the
/// node's list (`CTreeNode::hsplist`). `score` and `subject_offset` are those of
/// `hsp->score` and `hsp->subject.offset`, which `s_DominateTest` reads.
#[derive(Clone, Copy)]
struct LinkedHsp {
    hsp: usize,
    sid: u32,
    begin: i32,
    end: i32,
    score: i32,
    subject_offset: i64,
    merit: i32,
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

/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:200-207
/// ```c
/// typedef struct CTreeNode {
///     Int4 begin;                /* left endpoint */
///     Int4 end;                  /* right endpoint */
///     struct CTreeNode *left;    /* left child */
///     struct CTreeNode *right;   /* right child */
///     LinkedHSP *hsplist;        /* hsps belong to this node, start with low merits */
/// } CTreeNode;
/// ```
/// `left` and `right` are node indices; a freed node is unlinked from its parent.
struct CTreeNode {
    begin: i32,
    end: i32,
    left: Option<usize>,
    right: Option<usize>,
    hsplist: VecDeque<usize>,
}

/// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:435
/// ```c
///    Int4 kNumHSPtoFork = 20;  /** number of HSP to trig forking children */
/// ```
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
    culling_max: i32,
    hsps: Vec<Option<H>>,
    entries: Vec<LinkedHsp>,
    nodes: Vec<CTreeNode>,
    trees: BTreeMap<usize, usize>,
    free_hsps: Vec<usize>,
    free_entries: Vec<usize>,
    free_nodes: Vec<usize>,
}

impl<H> CullingForest<H> {
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:487-493
    /// ```c
    /// static int
    /// s_BlastHSPCullingInit(void* data, void* results)
    /// {
    ///     BlastHSPCullingData * cull_data = data;
    ///     cull_data->c_tree = calloc(cull_data->num_contexts, sizeof(CTreeNode *));
    ///     return 0;
    /// }
    /// ```
    pub(crate) fn new(culling_max: i32) -> Self {
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
    /// above 1. The sum is an `Int4`, which wraps for a limit of 2147483645 or more.
    ///
    /// NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:338-340
    /// ```c
    ///             if(params->culling_max > 1){
    ///             	params->culling_max += 3;
    ///             }
    /// ```
    pub(crate) fn prelim_max(limit: i32) -> i32 {
        if limit > 1 {
            limit.wrapping_add(3)
        } else {
            limit
        }
    }

    /// One HSP of `s_BlastHSPCullingRun`: the HSP of subject `sid`, with the keys `key`, goes
    /// into the tree `cid` (made with the length `qlen` on first use) unless `culling_max`
    /// HSPs already dominate it; an HSP that is not saved is freed.
    ///
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:612-633
    /// ```c
    ///       A.sid   = hsp_list->oid;
    ///       A.merit = params->culling_max;
    /// ...
    ///       if (! c_tree[A.cid]) {
    ///          c_tree[A.cid] = s_CTreeNew(qlen);
    ///       }
    ///
    ///       if(s_SaveHSP(c_tree[A.cid], &A)){
    ///     	 hsp_list->hsp_array[i] = NULL;
    ///       }
    /// ```
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
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:516-581
    /// ```c
    ///    for (cid=0; cid < cull_data->num_contexts; ++cid) {
    ///       if (c_tree[cid]) {
    ///         qid = Blast_GetQueryIndexFromContext(cid, params->program);
    ///         if (!results->hitlist_array[qid]) {
    ///             results->hitlist_array[qid] = Blast_HitListNew(params->prelim_hitlist_size);
    ///         }
    ///         hitlist = results->hitlist_array[qid];
    /// ...
    ///          cull_list = s_RipHSPOffCTree(c_tree[cid]);
    /// ...
    ///          while (cull_list) {
    ///             p = cull_list;
    ///             /* test to see if new hsplist has already been allocated */
    ///             allocated = FALSE;
    ///             for (sid=0; sid<hitlist->hsplist_count; ++sid) {
    ///                 list = hitlist->hsplist_array[sid];
    ///                 if (p->sid == list->oid) {
    ///                    allocated = TRUE;
    ///                    break;
    ///                 }
    ///             }
    /// ...
    ///             list->hsp_array[id] = p->hsp;
    ///             list->hspcnt++;
    /// ...
    ///          for (sid=0; sid < hitlist->hsplist_count; ++sid) {
    ///             list = hitlist->hsplist_array[sid];
    /// ...
    ///             Blast_HSPListSortByScore(list);
    /// ```
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

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:79-120
    /// ```c
    /// static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
    ///     Int8 b1 = p->begin;
    ///     Int8 b2 = y->begin;
    ///     Int8 e1 = p->end;
    ///     Int8 e2 = y->end;
    ///     Int8 s1 = p->hsp->score;
    ///     Int8 s2 = y->hsp->score;
    ///     Int8 l1 = e1 - b1;
    ///     Int8 l2 = e2 - b2;
    ///     Int8 overlap = MIN(e1,e2) - MAX(b1,b2);
    ///     Int8 d = 0;
    ///
    ///     // If not overlap by more than 50%
    ///     if(2 *overlap < l2) {
    ///     	return FALSE;
    ///     }
    /// ...
    ///     d  = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
    ///     // If identical, use oid as tie breaker
    ///     if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
    ///     	if(s1 != s2) {
    ///     		return (s1>s2);
    ///     	}
    ///     	if(p->sid != y->sid) {
    ///     		return (p->sid < y->sid);
    ///     	}
    ///
    ///     	if(p->hsp->subject.offset > y->hsp->subject.offset) {
    ///     		return FALSE;
    ///     	}
    ///     	return TRUE;
    ///     }
    ///
    ///    	if (d < 0) {
    ///    		return FALSE;
    ///     }
    ///
    ///     return TRUE;
    /// }
    /// ```
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

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:123-133
    /// ```c
    /// static Boolean s_FullPass(LinkedHSP *list, LinkedHSP *y) {
    ///     LinkedHSP *p = list;
    ///     while (p) {
    ///        if (s_DominateTest(p, y)) {
    ///           (y->merit)--;
    ///           if (y->merit <= 0) return FALSE;
    ///        }
    ///        p = p->next;
    ///     }
    ///     return TRUE;
    /// }
    /// ```
    /// `merit` is an `Int4`; the shipped NCBI binary wraps it (a culling limit near
    /// `INT4_MAX`, plus 3 in the preliminary stage, starts negative).
    fn full_pass(&self, node: usize, y: &mut LinkedHsp) -> bool {
        for &p in &self.nodes[node].hsplist {
            if self.dominate_test(&self.entries[p], y) {
                y.merit = y.merit.wrapping_sub(1);
                if y.merit <= 0 {
                    return false;
                }
            }
        }
        true
    }

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:136-163
    /// ```c
    /// static Int4 s_ProcessHSPList(LinkedHSP **list, LinkedHSP *y) {
    ///     Int4 num = 0;
    ///     LinkedHSP *p = *list, *q, *r;
    ///     q = p;
    ///     while (p) {
    ///        ++num;
    ///        r = p;
    ///        p = p->next;
    ///        if (r != y && s_DominateTest(y, r)) {
    ///           (r->merit)--;
    ///           if (r->merit <= 0) {
    /// ...
    ///              --num;
    ///              s_HSPFree(r);
    ///           } else {
    ///              q = r;
    ///           }
    ///        } else {
    ///           q = r;
    ///        }
    ///     }
    ///     return num;
    /// }
    /// ```
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
                let merit = self.entries[r].merit.wrapping_sub(1);
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

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:166-189
    /// ```c
    /// static Int4 s_MarkDownHSPList(LinkedHSP **list) {
    ///     Int4 num = 0;
    ///     LinkedHSP *p = *list, *q, *r;
    ///     q = p;
    ///     while (p) {
    ///        ++num;
    ///        r = p;
    ///        p = p->next;
    ///        (r->merit)--;
    ///        if (r->merit <= 0) {
    /// ...
    ///           --num;
    ///           s_HSPFree(r);
    ///        } else {
    ///           q = r;
    ///        }
    ///     }
    ///     return num;
    /// }
    /// ```
    fn mark_down_hsp_list(&mut self, node: usize) -> i32 {
        let mut list = std::mem::take(&mut self.nodes[node].hsplist);
        let mut num = 0;
        list.retain_mut(|r| {
            let r = *r;
            num += 1;
            let merit = self.entries[r].merit.wrapping_sub(1);
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

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:228-247
    /// ```c
    /// static CTreeNode * s_CTreeNodeNew(CTreeNode * parent, ECTreeChild dir) {
    ///     Int4 midpt;
    ///     CTreeNode * node = s_GetNode();
    ///
    ///     node->left    = NULL;
    ///     node->right   = NULL;
    ///     node->hsplist = NULL;
    ///
    ///     if (!parent) return node;
    ///
    ///     midpt = (parent->begin + parent->end) / 2;
    ///     if (dir == eLeft) {
    ///        node->begin = parent->begin;
    ///        node->end   = midpt;
    ///     } else {
    ///        node->begin = midpt;
    ///        node->end   = parent->end;
    ///     }
    ///     return node;
    /// }
    /// ```
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

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:376-381
    /// ```c
    /// static CTreeNode * s_CTreeNew(Int4 qlen) {
    ///     CTreeNode * tree = s_CTreeNodeNew(NULL, eLeft);
    ///     tree->begin = 0;
    ///     tree->end   = qlen;
    ///     return tree;
    /// }
    /// ```
    fn ctree_new(&mut self, qlen: i32) -> usize {
        let tree = self.ctree_node_new(None, true);
        self.nodes[tree].begin = 0;
        self.nodes[tree].end = qlen;
        tree
    }

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:258-299
    /// ```c
    ///     p = node->hsplist;
    ///     q = p;  /* q is predecessor of p */
    ///     midpt = (node->begin + node->end) /2;
    ///     while(p) {
    ///       child = NULL;
    ///       r = p;
    ///       if (p->end < midpt) {
    ///          if (!node->left) {
    ///             node->left = s_CTreeNodeNew(node, eLeft);
    ///          }
    ///          child = node->left;
    ///       } else if (p->begin > midpt) {
    ///          if (!node->right) {
    ///             node->right = s_CTreeNodeNew(node, eRight);
    ///          }
    ///          child = node->right;
    ///       }
    ///       p = p->next;
    ///       if (child) {
    ///          /* remove r from parent list */
    /// ...
    ///          /* and put it on the child */
    ///          s_AddHSPtoList(&(child->hsplist), r);
    ///       } else {
    ///          q = r;
    ///       }
    ///    }
    /// ```
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:193-197
    /// ```c
    /// static void s_AddHSPtoList(LinkedHSP **list, LinkedHSP *y) {
    ///     y->next = *list;
    ///     *list = y;
    ///     return;
    /// }
    /// ```
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

    /// Returns whether the node is freed (NCBI sets `*node = NULL`).
    ///
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:319-330
    /// ```c
    /// static void s_MarkDownCTree(CTreeNode ** node) {
    ///    if (! (*node)) return;
    ///
    ///    s_MarkDownCTree(&((*node)->left));
    ///    s_MarkDownCTree(&((*node)->right));
    ///    if ( s_MarkDownHSPList(&((*node)->hsplist)) <= 0
    ///      && !(*node)->left && !(*node)->right) {
    ///         s_CTreeNodeFree(*node);
    ///         *node = NULL;
    ///    }
    ///    return;
    /// }
    /// ```
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

    /// Returns whether the node is freed (NCBI sets `*node = NULL`).
    ///
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:334-370
    /// ```c
    /// static void s_ProcessCTree(CTreeNode ** node, LinkedHSP *x) {
    ///    Int4 midpt;
    ///
    ///    if (! (*node)) return;
    ///
    ///    /* first test if x includes the full range covered by node */
    ///    if (x->begin <= (*node)->begin && x->end >= (*node)->end) {
    ///       s_MarkDownCTree(node);
    ///       return;
    ///    }
    ///
    ///    /* if node reaches the leaves*/
    ///    if (!(*node)->left && !(*node)->right) {
    ///       if (s_ProcessHSPList(&((*node)->hsplist), x) <= 0) {
    ///           s_CTreeNodeFree(*node);
    ///           *node = NULL;
    ///       }
    ///       return;
    ///    }
    ///
    ///    /* recursive case */
    ///    midpt = ((*node)->begin + (*node)->end) / 2;
    ///    if (x->end < midpt) {
    ///       s_ProcessCTree(&((*node)->left), x);
    ///    } else if (x->begin > midpt) {
    ///       s_ProcessCTree(&((*node)->right), x);
    ///    } else {
    ///       s_ProcessCTree(&((*node)->left), x);
    ///       s_ProcessCTree(&((*node)->right), x);
    ///       if (s_ProcessHSPList(&((*node)->hsplist), x) <= 0
    ///        && !(*node)->left && !(*node)->right) {
    ///           s_CTreeNodeFree(*node);
    ///           *node = NULL;
    ///       }
    ///    }
    ///    return;
    /// }
    /// ```
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

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:396-423
    /// ```c
    /// static LinkedHSP * s_RipHSPOffCTree(CTreeNode *tree) {
    ///     LinkedHSP *q, *p;
    ///
    ///     if (!tree) return NULL;
    ///
    ///     q = tree->hsplist;
    ///     tree->hsplist = NULL;
    ///
    ///     /* grab left child */
    /// ...
    ///        p->next = s_RipHSPOffCTree(tree->left);
    /// ...
    ///     /* grab right child */
    /// ...
    ///        p->next = s_RipHSPOffCTree(tree->right);
    ///     }
    ///
    ///     return q;
    /// }
    /// ```
    fn rip_hsp_off_ctree(&mut self, tree: Option<usize>, out: &mut Vec<usize>) {
        let Some(tree) = tree else {
            return;
        };
        out.extend(std::mem::take(&mut self.nodes[tree].hsplist));
        let (left, right) = (self.nodes[tree].left, self.nodes[tree].right);
        self.rip_hsp_off_ctree(left, out);
        self.rip_hsp_off_ctree(right, out);
    }

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:430-470
    /// ```c
    /// static Boolean s_SaveHSP(CTreeNode *tree, LinkedHSP *A) {
    /// ...
    ///    /* Descend the tree */
    ///    while (tree) {
    ///       ASSERT(tree->begin <= A->begin);
    ///       ASSERT(tree->end   >= A->end);
    ///
    ///       if (! s_FullPass(tree->hsplist, A)) return FALSE;
    ///       midpt = (tree->begin + tree->end) /2;
    ///       node = tree;  /* record the last valid position */
    ///       if      (A->end   < midpt) tree = tree->left;
    ///       else if (A->begin > midpt) tree = tree->right;
    ///       else break;
    ///    }
    ///
    ///    /* if we get here, A is valid. copy and insert A at node */
    ///    x = s_HSPCopy(A);
    ///    s_AddHSPtoList(&(node->hsplist), x);
    ///
    ///    /* if this is the leaf, calculate update hsp number */
    ///    if (!node->left && !node->right) {
    ///        /* check for domination */
    ///        if (s_ProcessHSPList(&(node->hsplist), x) >= kNumHSPtoFork) {
    ///            /* fork this node into sub trees */
    ///            s_ForkChildren(node);
    ///        }
    ///        return TRUE;
    ///    }
    ///
    ///    /* check domination */
    ///    s_ProcessHSPList(&(node->hsplist), x);
    ///    s_ProcessCTree(&(node->left), x);
    ///    s_ProcessCTree(&(node->right), x);
    ///    return TRUE;
    /// }
    /// ```
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
mod tests {
    use super::*;

    fn key(begin: i32, end: i32, score: i32) -> CullKey {
        CullKey {
            begin,
            end,
            score,
            subject_offset: i64::from(begin),
        }
    }

    /// The payload of a test HSP is its score and its position in the input.
    type Payload = (i32, usize);

    fn forest(limit: i32) -> CullingForest<Payload> {
        CullingForest::new(limit)
    }

    fn by_score(a: &Payload, b: &Payload) -> Ordering {
        b.0.cmp(&a.0)
    }

    fn kept(forest: CullingForest<Payload>) -> Vec<(u32, Vec<i32>)> {
        forest
            .finalize(|_| 0, by_score)
            .into_iter()
            .flat_map(|(_, lists)| lists)
            .map(|(oid, list)| (oid, list.into_iter().map(|h| h.0).collect()))
            .collect()
    }

    // NCBI c++/src/algo/blast/api/setup_factory.cpp:338-340: the preliminary writer adds
    // 3 to a limit above 1, an Int4 sum that wraps.
    #[test]
    fn the_preliminary_limit_adds_three_above_one_and_wraps() {
        assert_eq!(CullingForest::<()>::prelim_max(0), 0);
        assert_eq!(CullingForest::<()>::prelim_max(1), 1);
        assert_eq!(CullingForest::<()>::prelim_max(2), 5);
        assert_eq!(CullingForest::<()>::prelim_max(2147483644), 2147483647);
        assert_eq!(CullingForest::<()>::prelim_max(2147483645), i32::MIN);
        assert_eq!(CullingForest::<()>::prelim_max(2147483647), -2147483646);
    }

    // NCBI frees a culled HSP at once (`s_HSPFree`): slots, entries and nodes are reused, so
    // a stream of ever better overlapping HSPs, each culling its predecessor, and a stream of
    // dominated HSPs both keep the arenas at a few items, and the survivors are unchanged.
    #[test]
    fn freed_slots_entries_and_nodes_are_reused() {
        let mut f = forest(1);
        for i in 0..1000 {
            f.save(0, 300, 0, key(10, 40, 10 + i), (10 + i, i as usize));
        }
        for i in 0..1000 {
            f.save(0, 300, 1, key(10, 40, 5), (5, i as usize));
        }
        assert!(f.hsps.len() <= 3 && f.entries.len() <= 3 && f.nodes.len() <= 3);
        assert_eq!(kept(f), vec![(0, vec![1009])]);
        // Disjoint HSPs fork the root; a long HSP then culls them all and frees the nodes.
        let mut f = forest(1);
        for i in 0..100 {
            f.save(0, 1000, 0, key(i * 5, i * 5 + 4, 10), (10, i as usize));
        }
        let nodes = f.nodes.len();
        f.save(0, 1000, 0, key(0, 600, 5000), (5000, 100));
        for i in 0..100 {
            f.save(0, 1000, 0, key(i * 5, i * 5 + 4, 10), (10, i as usize));
        }
        assert_eq!(f.hsps.len(), 101);
        assert!(f.nodes.len() <= nodes + 3);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:136-163: an HSP never dominates
    // itself (`r != y`), so a limit of 1 keeps the best HSP of overlapping ones, whichever
    // comes first.
    #[test]
    fn a_limit_of_one_keeps_the_dominating_hsp() {
        let mut f = forest(1);
        f.save(0, 300, 0, key(10, 40, 50), (50, 0));
        f.save(0, 300, 0, key(12, 40, 30), (30, 1));
        assert_eq!(kept(f), vec![(0, vec![50])]);
        let mut f = forest(1);
        f.save(0, 300, 0, key(12, 40, 30), (30, 0));
        f.save(0, 300, 0, key(10, 40, 50), (50, 1));
        assert_eq!(kept(f), vec![(0, vec![50])]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:79-120: with equal score and range
    // the lower OID dominates; another tree has its own HSPs.
    #[test]
    fn ties_are_broken_by_oid_and_trees_do_not_interact() {
        let mut f = forest(1);
        f.save(0, 300, 1, key(10, 40, 50), (50, 0));
        f.save(0, 300, 0, key(10, 40, 50), (50, 1));
        f.save(1, 300, 0, key(10, 40, 50), (50, 2));
        let all = f.finalize(|cid| cid, by_score);
        assert_eq!(all.len(), 2);
        assert_eq!(all[0].1, vec![(0, vec![(50, 1)])]);
        assert_eq!(all[1].1, vec![(0, vec![(50, 2)])]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:100-112: with the same OID the HSP
    // with the lower subject offset dominates, and equal keys dominate each other, so the
    // HSP saved later loses.
    #[test]
    fn same_oid_ties_follow_the_subject_offset_and_the_insertion_order() {
        let mut f = forest(1);
        let mut late = key(10, 40, 50);
        late.subject_offset = 90;
        f.save(0, 300, 0, late, (50, 0));
        f.save(0, 300, 0, key(10, 40, 50), (50, 1));
        let lists = f.finalize(|_| 0, by_score);
        assert_eq!(lists[0].1, vec![(0, vec![(50, 1)])]);
        let mut f = forest(1);
        f.save(0, 300, 0, key(10, 40, 50), (50, 0));
        f.save(0, 300, 0, key(10, 40, 50), (50, 1));
        let lists = f.finalize(|_| 0, by_score);
        assert_eq!(lists[0].1, vec![(0, vec![(50, 0)])]);
    }

    // The merit is an Int4 that NCBI's binary wraps: a negative merit loses an HSP at its
    // first dominator, and i32::MIN never reaches 0 after one decrement.
    #[test]
    fn merits_wrap_as_ncbis_int4() {
        let mut f = forest(CullingForest::<()>::prelim_max(2147483646));
        f.save(0, 300, 0, key(10, 40, 50), (50, 0));
        f.save(0, 300, 0, key(12, 40, 30), (30, 1));
        assert_eq!(kept(f), vec![(0, vec![50])]);
        let mut f = forest(CullingForest::<()>::prelim_max(2147483645));
        f.save(0, 300, 0, key(10, 40, 50), (50, 0));
        f.save(0, 300, 0, key(12, 40, 30), (30, 1));
        assert_eq!(kept(f), vec![(0, vec![50, 30])]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:258-299,453-458: a leaf with 20
    // HSPs forks; HSPs that do not overlap are all kept.
    #[test]
    fn many_disjoint_hsps_fork_and_survive() {
        let mut f = forest(1);
        for i in 0..40 {
            f.save(0, 3000, 0, key(i * 20, i * 20 + 10, 30), (30, i as usize));
        }
        assert_eq!(kept(f)[0].1.len(), 40);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:396-423,516-581: the rip takes a
    // node's list, then the left subtree, then the right subtree; the lists of a query are
    // made in the order of the first HSP of each subject.
    #[test]
    fn final_layout_follows_the_rip_order() {
        let mut f = forest(5);
        // The first HSP spans the middle (150) and stays in the root; twenty more fork
        // the root into a left and a right child.
        f.save(0, 300, 7, key(140, 160, 40), (40, 0));
        for i in 0..10 {
            f.save(
                0,
                300,
                3,
                key(i * 10, i * 10 + 5, 20),
                (20, 100 + i as usize),
            );
        }
        for i in 0..10 {
            f.save(
                0,
                300,
                5,
                key(200 + i * 10, 205 + i * 10, 21),
                (21, 200 + i as usize),
            );
        }
        let lists = f.finalize(|_| 0, by_score);
        let oids: Vec<u32> = lists[0].1.iter().map(|(oid, _)| *oid).collect();
        assert_eq!(oids, vec![7, 3, 5]);
        assert_eq!(lists[0].1[1].1.len(), 10);
        assert_eq!(lists[0].1[2].1.len(), 10);
    }

    // Queries come in order, each with its own lists and the lists of an OID found again.
    #[test]
    fn queries_come_in_order_with_their_own_lists() {
        let mut f = forest(10);
        f.save(3, 100, 4, key(0, 10, 9), (9, 0));
        f.save(1, 100, 4, key(0, 10, 8), (8, 1));
        f.save(2, 100, 2, key(50, 60, 7), (7, 2));
        f.save(2, 100, 4, key(70, 80, 6), (6, 3));
        let out = f.finalize(|cid| cid / 2, by_score);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].0, 0);
        assert_eq!(out[0].1, vec![(4, vec![(8, 1)])]);
        assert_eq!(out[1].0, 1);
        // Tree 2 rips its newest HSP first (subject 4), then tree 3 adds to subject 4's list,
        // which is sorted again.
        assert_eq!(out[1].1, vec![(4, vec![(9, 0), (6, 3)]), (2, vec![(7, 2)])]);
    }
}
