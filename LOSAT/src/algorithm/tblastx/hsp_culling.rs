//! NCBI BLAST HSP culling (`-culling_limit`) for TBLASTX.
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
//! NCBI culls twice when `-culling_limit` is given (`report::culled_hit_order`): a writer
//! in the preliminary stage and a pipe after the traceback.

use std::collections::{BTreeMap, VecDeque};

use super::blast_engine::TblastxHsp;

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
/// `hsp` is the index of the HSP in the writer's arena; `next` is the order of the node's
/// list (`CTreeNode::hsplist`).
#[derive(Clone, Copy)]
struct LinkedHsp {
    hsp: usize,
    sid: u32,
    begin: i32,
    end: i32,
    merit: i32,
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

/// The query context of a TBLASTX HSP within its query: NCBI's contexts of a translated
/// query are its frames 1, 2, 3, -1, -2, -3 (`BLAST_GetAllTranslations`).
///
/// NCBI reference: c++/src/algo/blast/core/blast_util.c:1211-1220
/// ```c
/// Int4 BLAST_FrameToContext(Int2 frame, EBlastProgramType program)
/// {
///     if (Blast_QueryIsTranslated(program) ||
///         Blast_SubjectIsTranslated(program)) {
///         ASSERT(frame >= -3 && frame <= 3 && frame != 0);
///         if (frame > 0) {
///             return frame - 1;
///         } else {
///             return 2 - frame;
///         }
/// ```
pub(crate) fn tblastx_context(frame: i32) -> usize {
    if frame > 0 {
        (frame - 1) as usize
    } else {
        (2 - frame) as usize
    }
}

/// The length in residues of each frame of a query of `length` nucleotides, as NCBI's
/// query information holds it (`contexts[context].query_length`), in context order.
///
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:193-195
/// ```c
///             for (unsigned int i = 0; i < kNumContexts; i++) {
///                 unsigned int prot_length =
///                     static_cast<unsigned int>(BLAST_GetTranslatedProteinLength(length, i));
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_util.c:922-929
/// ```c
/// size_t
/// BLAST_GetTranslatedProteinLength(size_t nucleotide_length, unsigned int context)
/// {
///     if (nucleotide_length == 0 || nucleotide_length <= context % CODON_LENGTH) {
///         return 0;
///     }
///     return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
/// }
/// ```
/// TBLASTX's queries are searched on both strands (`eNa_strand_both`), so every context
/// has this length.
pub(crate) fn tblastx_context_lengths(nucleotide_length: usize) -> [i32; 6] {
    let mut lengths = [0; 6];
    for (context, slot) in lengths.iter_mut().enumerate() {
        let shift = context % 3;
        *slot = if nucleotide_length == 0 || nucleotide_length <= shift {
            0
        } else {
            i32::try_from((nucleotide_length - shift) / 3).unwrap_or(i32::MAX)
        };
    }
    lengths
}

/// NCBI's culling writer (or pipe) over the HSPs of one query: `s_BlastHSPCullingInit`
/// (`new`), `s_BlastHSPCullingRun` (`run`) and `s_BlastHSPCullingFinal` (`finalize`).
pub(crate) struct CullingWriter {
    culling_max: i32,
    context_lengths: [i32; 6],
    hsps: Vec<Option<TblastxHsp>>,
    entries: Vec<LinkedHsp>,
    nodes: Vec<CTreeNode>,
    /// `c_tree[cid]`: the root node of each context with a tree.
    trees: BTreeMap<usize, usize>,
}

impl CullingWriter {
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
    pub(crate) fn new(culling_max: i32, context_lengths: [i32; 6]) -> Self {
        Self {
            culling_max,
            context_lengths,
            hsps: Vec::new(),
            entries: Vec::new(),
            nodes: Vec::new(),
            trees: BTreeMap::new(),
        }
    }

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:602-644
    /// ```c
    ///    for (i=0; i<hsp_list->hspcnt; ++i) {
    ///       /* wrap the hsp with a LinkedHSP structure */
    ///       A.hsp   = hsp_list->hsp_array[i];
    ///       A.cid   = isBlastn ? (A.hsp->context  - A.hsp->context % NUM_STRANDS) : A.hsp->context;
    ///       A.sid   = hsp_list->oid;
    ///       A.merit = params->culling_max;
    ///       qlen    = cull_data->query_info->contexts[A.hsp->context].query_length;
    ///       if(isBlastn && (A.hsp->context % NUM_STRANDS)) {
    ///     	  A.begin = qlen - A.hsp->query.end;
    ///     	  A.end   = qlen -  A.hsp->query.offset;
    ///       }
    ///       else {
    ///     	  A.begin = A.hsp->query.offset;
    ///     	  A.end   = A.hsp->query.end;
    ///       }
    ///       A.next  = NULL;
    ///
    ///       if (! c_tree[A.cid]) {
    ///          c_tree[A.cid] = s_CTreeNew(qlen);
    ///       }
    ///
    ///       if(s_SaveHSP(c_tree[A.cid], &A)){
    ///     	 hsp_list->hsp_array[i] = NULL;
    ///       }
    ///    }
    /// ```
    /// TBLASTX is not BLASTN: the context is the HSP's own and the offsets are the frame
    /// offsets of `Blast_HSPInit` (`sort_query_offset`, `sort_query_end`). The HSPs of
    /// `hsp_list` come in its order; an HSP that is not saved is freed.
    pub(crate) fn run(&mut self, hsp_list: Vec<TblastxHsp>, oid: u32) {
        for hsp in hsp_list {
            let cid = tblastx_context(hsp.hit.query_frame);
            let qlen = self.context_lengths[cid];
            let index = self.hsps.len();
            let mut a = LinkedHsp {
                hsp: index,
                sid: oid,
                begin: hsp.hit.sort_query_offset as i32,
                end: hsp.hit.sort_query_end as i32,
                merit: self.culling_max,
            };
            self.hsps.push(Some(hsp));
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
            }
        }
    }

    /// The HSPs that the trees keep, as `s_BlastHSPCullingFinal` puts them into the hit
    /// list of the query: the contexts in order, each tree ripped into one list, a new HSP
    /// list for each subject in the order of its first HSP, and every HSP list sorted by
    /// score (`Blast_HSPListSortByScore`, glibc's stable `qsort`).
    ///
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:516-581
    /// ```c
    ///    for (cid=0; cid < cull_data->num_contexts; ++cid) {
    ///       if (c_tree[cid]) {
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
    ///             best_evalue = (double) INT4_MAX;
    ///             for (id=0; id < list->hspcnt; ++id) {
    ///                 best_evalue = MIN(list->hsp_array[id]->evalue, best_evalue);
    ///             }
    ///             Blast_HSPListSortByScore(list);
    /// ```
    pub(crate) fn finalize(mut self) -> Vec<(u32, Vec<TblastxHsp>)> {
        let mut lists: Vec<(u32, Vec<TblastxHsp>)> = Vec::new();
        let roots: Vec<usize> = self.trees.values().copied().collect();
        for root in roots {
            let mut ripped = Vec::new();
            self.rip_hsp_off_ctree(Some(root), &mut ripped);
            for entry in ripped {
                let LinkedHsp { hsp, sid, .. } = self.entries[entry];
                let hsp = self.hsps[hsp].take().expect("an HSP is ripped once");
                match lists.iter_mut().find(|(oid, _)| *oid == sid) {
                    Some((_, list)) => list.push(hsp),
                    None => lists.push((sid, vec![hsp])),
                }
            }
            for (_, list) in &mut lists {
                list.sort_by(|a, b| crate::common::score_compare_hsps(&a.hit, &b.hit));
            }
        }
        lists
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
        let p_hsp = self.hsps[p.hsp].as_ref().expect("a saved HSP");
        let y_hsp = self.hsps[y.hsp].as_ref().expect("a saved HSP");
        let (b1, b2) = (i64::from(p.begin), i64::from(y.begin));
        let (e1, e2) = (i64::from(p.end), i64::from(y.end));
        let (s1, s2) = (
            i64::from(p_hsp.hit.raw_score),
            i64::from(y_hsp.hit.raw_score),
        );
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
            if p_hsp.hit.sort_subject_offset > y_hsp.hit.sort_subject_offset {
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
        let list = std::mem::take(&mut self.nodes[node].hsplist);
        let mut kept = VecDeque::with_capacity(list.len());
        let mut num = 0;
        let y_entry = self.entries[y];
        for r in list {
            num += 1;
            if r != y && self.dominate_test(&y_entry, &self.entries[r]) {
                let merit = self.entries[r].merit.wrapping_sub(1);
                self.entries[r].merit = merit;
                if merit <= 0 {
                    num -= 1;
                    self.hsps[self.entries[r].hsp] = None;
                    continue;
                }
            }
            kept.push_back(r);
        }
        self.nodes[node].hsplist = kept;
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
        let list = std::mem::take(&mut self.nodes[node].hsplist);
        let mut kept = VecDeque::with_capacity(list.len());
        let mut num = 0;
        for r in list {
            num += 1;
            let merit = self.entries[r].merit.wrapping_sub(1);
            self.entries[r].merit = merit;
            if merit <= 0 {
                num -= 1;
                self.hsps[self.entries[r].hsp] = None;
                continue;
            }
            kept.push_back(r);
        }
        self.nodes[node].hsplist = kept;
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
        self.nodes.push(CTreeNode {
            begin,
            end,
            left: None,
            right: None,
            hsplist: VecDeque::new(),
        });
        self.nodes.len() - 1
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
        let list = std::mem::take(&mut self.nodes[node].hsplist);
        let mut kept = VecDeque::with_capacity(list.len());
        for r in list {
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
                Some(child) => self.nodes[child].hsplist.push_front(r),
                None => kept.push_back(r),
            }
        }
        self.nodes[node].hsplist = kept;
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
            }
        }
        if let Some(right) = self.nodes[node].right {
            if self.mark_down_ctree(right) {
                self.nodes[node].right = None;
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
                }
            }
        };
        let right = |writer: &mut Self| {
            if let Some(right) = writer.nodes[node].right {
                if writer.process_ctree(right, x) {
                    writer.nodes[node].right = None;
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
        let x = self.entries.len();
        self.entries.push(*a);
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
            }
        }
        if let Some(right) = self.nodes[node].right {
            if self.process_ctree(right, x) {
                self.nodes[node].right = None;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Hit;

    fn hsp(oid: u32, frame: i32, begin: usize, end: usize, score: i32) -> TblastxHsp {
        TblastxHsp {
            hit: Hit {
                identity: 0.0,
                length: end - begin,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 3,
                s_start: 1,
                s_end: 3,
                e_value: 1e-5,
                bit_score: 0.0,
                num_ident: 0,
                query_frame: frame,
                query_length: 0,
                q_idx: 0,
                s_idx: oid,
                raw_score: score,
                sort_query_offset: begin,
                sort_query_end: end,
                sort_subject_offset: begin,
                sort_subject_end: end,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 0,
            },
            subject_frame: 1,
            num: 1,
        }
    }

    fn kept(writer: CullingWriter) -> Vec<(u32, i32)> {
        writer
            .finalize()
            .into_iter()
            .flat_map(|(oid, list)| list.into_iter().map(move |h| (oid, h.hit.raw_score)))
            .collect()
    }

    // NCBI c++/src/algo/blast/core/blast_util.c:922-929,1211-1220: frames 1, 2, 3, -1,
    // -2, -3 are the contexts 0..5, of (length - context % 3) / 3 residues.
    #[test]
    fn contexts_follow_ncbis_frames() {
        let order: Vec<usize> = [1, 2, 3, -1, -2, -3]
            .iter()
            .map(|&frame| tblastx_context(frame))
            .collect();
        assert_eq!(order, vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(tblastx_context_lengths(31), [10, 10, 9, 10, 10, 9]);
        assert_eq!(tblastx_context_lengths(0), [0; 6]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:136-163: an HSP never dominates
    // itself (`r != y`), so a limit of 1 keeps the best HSP of overlapping ones.
    #[test]
    fn a_limit_of_one_keeps_the_dominating_hsp() {
        let mut writer = CullingWriter::new(1, tblastx_context_lengths(300));
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
        assert_eq!(kept(writer), vec![(0, 50)]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:79-120: with equal score and range
    // the lower OID dominates; other contexts have their own trees.
    #[test]
    fn ties_are_broken_by_oid_and_contexts_do_not_interact() {
        let mut writer = CullingWriter::new(1, tblastx_context_lengths(300));
        writer.run(vec![hsp(1, 1, 10, 40, 50)], 1);
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, -1, 10, 40, 50)], 0);
        assert_eq!(kept(writer), vec![(0, 50), (0, 50)]);
    }

    // The merit is an Int4 that NCBI's binary wraps: a negative merit loses an HSP at its
    // first dominator, and i32::MIN never reaches 0 after one decrement.
    #[test]
    fn merits_wrap_as_ncbis_int4() {
        let lengths = tblastx_context_lengths(300);
        let mut writer = CullingWriter::new(2147483646i32.wrapping_add(3), lengths);
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
        assert_eq!(kept(writer), vec![(0, 50)]);
        let mut writer = CullingWriter::new(2147483645i32.wrapping_add(3), lengths);
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
        assert_eq!(kept(writer), vec![(0, 50), (0, 30)]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:258-299,453-458: a leaf with 20
    // HSPs forks; HSPs that do not overlap are all kept.
    #[test]
    fn many_disjoint_hsps_fork_and_survive() {
        let mut writer = CullingWriter::new(1, tblastx_context_lengths(3000));
        let list: Vec<TblastxHsp> = (0..40)
            .map(|i| hsp(0, 1, i * 20, i * 20 + 10, 30))
            .collect();
        writer.run(list, 0);
        assert_eq!(kept(writer).len(), 40);
    }
}
