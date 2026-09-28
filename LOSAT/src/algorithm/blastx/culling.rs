//! BLASTX context-local culling writer and traceback pipe.
use super::{
    preliminary::compare_score,
    results::{HitList, Hsp, HspList, ResultsObserver},
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:52-60
// ```c++
// typedef struct LinkedHSP {
//     BlastHSP * hsp;
//     Int4 cid;    /* context id for hsp */
//     Int4 sid;    /* OID for hsp*/
//     Int4 begin;  /* query offset in plus strand */
//     Int4 end;    /* query end in plus strand */
//     Int4 merit;  /* how many other hsps in the tree dominates me? */
//     struct LinkedHSP *next;
// } LinkedHSP;
// ```
#[derive(Clone)]
struct Entry {
    id: usize,
    oid: usize,
    hsp: Hsp,
    merit: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120
// ```c++
// static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
//     Int8 b1 = p->begin;
//     Int8 b2 = y->begin;
//     Int8 e1 = p->end;
//     Int8 e2 = y->end;
//     Int8 s1 = p->hsp->score;
//     Int8 s2 = y->hsp->score;
//     Int8 l1 = e1 - b1;
//     Int8 l2 = e2 - b2;
//     Int8 overlap = MIN(e1,e2) - MAX(b1,b2);
//     Int8 d = 0;
//
//     // If not overlap by more than 50%
//     if(2 *overlap < l2) {
//     	return FALSE;
//     }
//
//     /* the main criterion:
//        2 * (%diff in score) + 1 * (%diff in length) */
//     //Int8 d  = 3*s1*l1 + s1*l2 - s2*l1 - 3*s2*l2;
//     d  = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
//     // If identical, use oid as tie breaker
//     if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
//     	if(s1 != s2) {
//     		return (s1>s2);
//     	}
//     	if(p->sid != y->sid) {
//     		return (p->sid < y->sid);
//     	}
//
//     	if(p->hsp->subject.offset > y->hsp->subject.offset) {
//     		return FALSE;
//     	}
//     	return TRUE;
//     }
//
//    	if (d < 0) {
//    		return FALSE;
//     }
//
//     return TRUE;
// }
// ```
fn dominates_source(a: &Entry, b: &Entry) -> bool {
    let ah = &a.hsp.hsp;
    let bh = &b.hsp.hsp;
    let (b1, b2, e1, e2, s1, s2) = (
        ah.q_start as i64,
        bh.q_start as i64,
        ah.q_end as i64,
        bh.q_end as i64,
        ah.score as i64,
        bh.score as i64,
    );
    let (l1, l2) = (e1 - b1, e2 - b2);
    let overlap = e1.min(e2) - b1.max(b2);
    if 2 * overlap < l2 {
        return false;
    }
    let d = 4 * s1 * l1 + 2 * s1 * l2 - 2 * s2 * l1 - 4 * s2 * l2;
    if (s1 == s2 && b1 == b2 && l1 == l2) || d == 0 {
        if s1 != s2 {
            return s1 > s2;
        }
        if a.oid != b.oid {
            return a.oid < b.oid;
        }
        return ah.s_start <= bh.s_start;
    }
    d >= 0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120
// ```c++
// static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
//     Int8 b1 = p->begin;
//     Int8 b2 = y->begin;
//     Int8 e1 = p->end;
//     Int8 e2 = y->end;
//     Int8 s1 = p->hsp->score;
//     Int8 s2 = y->hsp->score;
//     Int8 l1 = e1 - b1;
//     Int8 l2 = e2 - b2;
//     Int8 overlap = MIN(e1,e2) - MAX(b1,b2);
//     Int8 d = 0;
//
//     // If not overlap by more than 50%
//     if(2 *overlap < l2) {
//     	return FALSE;
//     }
//
//     /* the main criterion:
//        2 * (%diff in score) + 1 * (%diff in length) */
//     //Int8 d  = 3*s1*l1 + s1*l2 - s2*l1 - 3*s2*l2;
//     d  = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
//     // If identical, use oid as tie breaker
//     if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
//     	if(s1 != s2) {
//     		return (s1>s2);
//     	}
//     	if(p->sid != y->sid) {
//     		return (p->sid < y->sid);
//     	}
//
//     	if(p->hsp->subject.offset > y->hsp->subject.offset) {
//     		return FALSE;
//     	}
//     	return TRUE;
//     }
//
//    	if (d < 0) {
//    		return FALSE;
//     }
//
//     return TRUE;
// }
// ```
fn dominates(a: &Entry, b: &Entry, trace: &mut dyn ResultsObserver) -> bool {
    let result = dominates_source(a, b);
    trace.cull_compare(a.oid, &a.hsp, b.oid, &b.hsp, result);
    result
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:201-249
// ```c++
// typedef struct CTreeNode {
//     Int4 begin;  /* left endpoint */
//     Int4 end;    /* right endpoint */
//     struct CTreeNode *left;    /* left child */
//     struct CTreeNode *right;   /* right child */
//     LinkedHSP *hsplist; /* hsps belong to this node, start with low merits */
// } CTreeNode;
//
// /*******Memory management layer (may be changed later)********/
// static CTreeNode * s_GetNode() {
//     return( (CTreeNode *) calloc(1,sizeof(CTreeNode)));
// }
//
// static CTreeNode * s_RetNode(CTreeNode * node) {
//     sfree(node);
//     return NULL;
// }
//
// /*************************************************************/
// /**  functions to manipulate Culling Tree (private)*/
//
// typedef enum ECTreeChild {
//     eLeft,
//     eRight,
// } ECTreeChild;
//
// /** Allocate and return a new node for use */
// static CTreeNode * s_CTreeNodeNew(CTreeNode * parent, ECTreeChild dir) {
//     Int4 midpt;
//     CTreeNode * node = s_GetNode();
//
//     node->left    = NULL;
//     node->right   = NULL;
//     node->hsplist = NULL;
//
//     if (!parent) return node;
//
//     midpt = (parent->begin + parent->end) / 2;
//     if (dir == eLeft) {
//         node->begin = parent->begin;
//         node->end   = midpt;
//     } else {
//         node->begin = midpt;
//         node->end   = parent->end;
//     }
//     return node;
// }
//
// /** Free an individual node */
// ```
struct Node {
    begin: i32,
    end: i32,
    left: Option<Box<Node>>,
    right: Option<Box<Node>>,
    list: Vec<Entry>,
}
impl Node {
    fn new(begin: i32, end: i32) -> Box<Self> {
        Box::new(Self {
            begin,
            end,
            left: None,
            right: None,
            list: Vec::new(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:136-189
    // ```c++
    // static Int4 s_ProcessHSPList(LinkedHSP **list, LinkedHSP *y) {
    //     Int4 num = 0;
    //     LinkedHSP *p = *list, *q, *r;
    //     q = p;
    //     while (p) {
    //        ++num;
    //        r = p;
    //        p = p->next;
    //        if (r != y && s_DominateTest(y, r)) {
    //           (r->merit)--;
    //           if (r->merit <= 0) {
    //              if (r == *list) {
    //                  *list = p;
    //                  q = p;
    //              } else {
    //                  q->next = p;
    //              }
    //              --num;
    //              s_HSPFree(r);
    //           } else {
    //              q = r;
    //           }
    //        } else {
    //           q = r;
    //        }
    //     }
    //     return num;
    // }
    //
    // /** decrease merit for all hsps in list; also returns the number of hsps in list */
    // static Int4 s_MarkDownHSPList(LinkedHSP **list) {
    //     Int4 num = 0;
    //     LinkedHSP *p = *list, *q, *r;
    //     q = p;
    //     while (p) {
    //        ++num;
    //        r = p;
    //        p = p->next;
    //        (r->merit)--;
    //        if (r->merit <= 0) {
    //           if (r == *list) {
    //               *list = p;
    //               q = p;
    //           } else {
    //               q->next = p;
    //           }
    //           --num;
    //           s_HSPFree(r);
    //        } else {
    //           q = r;
    //        }
    //     }
    //     return num;
    // }
    // ```
    fn process_list(&mut self, entry: Option<&Entry>, trace: &mut dyn ResultsObserver) {
        self.list.retain_mut(|h| {
            if entry.is_none_or(|x| x.id != h.id && dominates(x, h, trace)) {
                h.merit -= 1;
            }
            if h.merit <= 0 {
                trace.cull_delete(h.oid, &h.hsp, h.merit);
            }
            h.merit > 0
        });
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:258-299
    // ```c++
    // static void s_ForkChildren(CTreeNode * node) {
    //     CTreeNode * child;
    //     LinkedHSP *p, *q, *r;
    //     Int4 midpt;
    //
    //     ASSERT(node != NULL);
    //     ASSERT(node->left ==NULL);
    //     ASSERT(node->right ==NULL);
    //
    //     p = node->hsplist;
    //     q = p;  /* q is predecessor of p */
    //     midpt = (node->begin + node->end) /2;
    //     while(p) {
    //       child = NULL;
    //       r = p;
    //       if (p->end < midpt) {
    //          if (!node->left) {
    //             node->left = s_CTreeNodeNew(node, eLeft);
    //          }
    //          child = node->left;
    //       } else if (p->begin > midpt) {
    //          if (!node->right) {
    //             node->right = s_CTreeNodeNew(node, eRight);
    //          }
    //          child = node->right;
    //       }
    //       p = p->next;
    //       if (child) {
    //          /* remove r from parent list */
    //          if (r == node->hsplist) {
    //              node->hsplist = p;
    //              q = p;
    //          } else {
    //              q->next = p;
    //          }
    //          /* and put it on the child */
    //          s_AddHSPtoList(&(child->hsplist), r);
    //       } else {
    //          q = r;
    //       }
    //    }
    // }
    // ```
    fn fork(&mut self, trace: &mut dyn ResultsObserver) {
        trace.cull_fork(self.begin, self.end);
        let mid = (self.begin + self.end) / 2;
        for h in std::mem::take(&mut self.list) {
            if h.hsp.hsp.q_end < mid {
                self.left
                    .get_or_insert_with(|| Self::new(self.begin, mid))
                    .list
                    .insert(0, h);
            } else if h.hsp.hsp.q_start > mid {
                self.right
                    .get_or_insert_with(|| Self::new(mid, self.end))
                    .list
                    .insert(0, h);
            } else {
                self.list.push(h);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:319-370
    // ```c++
    // static void s_MarkDownCTree(CTreeNode ** node) {
    //    if (! (*node)) return;
    //
    //    s_MarkDownCTree(&((*node)->left));
    //    s_MarkDownCTree(&((*node)->right));
    //    if ( s_MarkDownHSPList(&((*node)->hsplist)) <= 0
    //      && !(*node)->left && !(*node)->right) {
    //         s_CTreeNodeFree(*node);
    //         *node = NULL;
    //    }
    //    return;
    // }
    //
    // /** recursively search and update merit hsps in culling tree
    //     due to addition of hsp x */
    // static void s_ProcessCTree(CTreeNode ** node, LinkedHSP *x) {
    //    Int4 midpt;
    //
    //    if (! (*node)) return;
    //
    //    /* first test if x includes the full range covered by node */
    //    if (x->begin <= (*node)->begin && x->end >= (*node)->end) {
    //       s_MarkDownCTree(node);
    //       return;
    //    }
    //
    //    /* if node reaches the leaves*/
    //    if (!(*node)->left && !(*node)->right) {
    //       if (s_ProcessHSPList(&((*node)->hsplist), x) <= 0) {
    //           s_CTreeNodeFree(*node);
    //           *node = NULL;
    //       }
    //       return;
    //    }
    //
    //    /* recursive case */
    //    midpt = ((*node)->begin + (*node)->end) / 2;
    //    if (x->end < midpt) {
    //       s_ProcessCTree(&((*node)->left), x);
    //    } else if (x->begin > midpt) {
    //       s_ProcessCTree(&((*node)->right), x);
    //    } else {
    //       s_ProcessCTree(&((*node)->left), x);
    //       s_ProcessCTree(&((*node)->right), x);
    //       if (s_ProcessHSPList(&((*node)->hsplist), x) <= 0
    //        && !(*node)->left && !(*node)->right) {
    //           s_CTreeNodeFree(*node);
    //           *node = NULL;
    //       }
    //    }
    //    return;
    // }
    // ```
    fn process_tree(
        tree: &mut Option<Box<Self>>,
        entry: Option<&Entry>,
        trace: &mut dyn ResultsObserver,
    ) {
        let Some(node) = tree.as_mut() else {
            return;
        };
        let covers =
            entry.is_none_or(|x| x.hsp.hsp.q_start <= node.begin && x.hsp.hsp.q_end >= node.end);
        if covers {
            Self::process_tree(&mut node.left, None, trace);
            Self::process_tree(&mut node.right, None, trace);
            node.process_list(None, trace);
        } else if node.left.is_none() && node.right.is_none() {
            node.process_list(entry, trace);
        } else {
            let x = entry.unwrap();
            let mid = (node.begin + node.end) / 2;
            if x.hsp.hsp.q_end < mid {
                Self::process_tree(&mut node.left, entry, trace);
                return;
            }
            if x.hsp.hsp.q_start > mid {
                Self::process_tree(&mut node.right, entry, trace);
                return;
            }
            Self::process_tree(&mut node.left, entry, trace);
            Self::process_tree(&mut node.right, entry, trace);
            node.process_list(entry, trace);
        }
        if node.list.is_empty() && node.left.is_none() && node.right.is_none() {
            *tree = None;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:430-470
    // ```c++
    // static Boolean s_SaveHSP(CTreeNode *tree, LinkedHSP *A) {
    //    Int4 midpt;
    //
    //    LinkedHSP *x;
    //    CTreeNode *node;
    //    Int4 kNumHSPtoFork = 20;  /** number of HSP to trig forking children */
    //
    //    ASSERT(tree != NULL);
    //    /* Descend the tree */
    //    while (tree) {
    //       ASSERT(tree->begin <= A->begin);
    //       ASSERT(tree->end   >= A->end);
    //
    //       if (! s_FullPass(tree->hsplist, A)) return FALSE;
    //       midpt = (tree->begin + tree->end) /2;
    //       node = tree;  /* record the last valid position */
    //       if      (A->end   < midpt) tree = tree->left;
    //       else if (A->begin > midpt) tree = tree->right;
    //       else break;
    //    }
    //
    //    /* if we get here, A is valid. copy and insert A at node */
    //    x = s_HSPCopy(A);
    //    s_AddHSPtoList(&(node->hsplist), x);
    //
    //    /* if this is the leaf, calculate update hsp number */
    //    if (!node->left && !node->right) {
    //        /* check for domination */
    //        if (s_ProcessHSPList(&(node->hsplist), x) >= kNumHSPtoFork) {
    //            /* fork this node into sub trees */
    //            s_ForkChildren(node);
    //        }
    //        return TRUE;
    //    }
    //
    //    /* check domination */
    //    s_ProcessHSPList(&(node->hsplist), x);
    //    s_ProcessCTree(&(node->left), x);
    //    s_ProcessCTree(&(node->right), x);
    //    return TRUE;
    // }
    // ```
    fn save(&mut self, mut entry: Entry, trace: &mut dyn ResultsObserver) -> bool {
        for h in &self.list {
            if dominates(h, &entry, trace) {
                entry.merit -= 1;
                if entry.merit <= 0 {
                    return false;
                }
            }
        }
        let mid = (self.begin + self.end) / 2;
        if entry.hsp.hsp.q_end < mid {
            if let Some(left) = &mut self.left {
                return left.save(entry, trace);
            }
        } else if entry.hsp.hsp.q_start > mid {
            if let Some(right) = &mut self.right {
                return right.save(entry, trace);
            }
        }
        self.list.insert(0, entry.clone());
        self.process_list(Some(&entry), trace);
        if self.left.is_none() && self.right.is_none() {
            if self.list.len() >= 20 {
                self.fork(trace);
            }
        } else {
            Self::process_tree(&mut self.left, Some(&entry), trace);
            Self::process_tree(&mut self.right, Some(&entry), trace);
        }
        true
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:396-423
    // ```c++
    // static LinkedHSP * s_RipHSPOffCTree(CTreeNode *tree) {
    //     LinkedHSP *q, *p;
    //
    //     if (!tree) return NULL;
    //
    //     q = tree->hsplist;
    //     tree->hsplist = NULL;
    //
    //     /* grab left child */
    //     if (!q) {
    //        q = s_RipHSPOffCTree(tree->left);
    //        p = q;
    //     } else {
    //        p = q;
    //        while(p->next) p=p->next;
    //        p->next = s_RipHSPOffCTree(tree->left);
    //     }
    //
    //     /* grab right child */
    //     if (!q) {
    //        q = s_RipHSPOffCTree(tree->right);
    //     } else {
    //        while(p->next) p=p->next;
    //        p->next = s_RipHSPOffCTree(tree->right);
    //     }
    //
    //     return q;
    // }
    // ```
    fn rip(self: Box<Self>, out: &mut Vec<Entry>) {
        out.extend(self.list);
        if let Some(l) = self.left {
            l.rip(out);
        }
        if let Some(r) = self.right {
            r.rip(out);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:604-647
// ```c++
// {
//    Int4 i, qlen;
//    LinkedHSP A;
//
//    BlastHSPCullingData * cull_data = data;
//    BlastHSPCullingParams* params = cull_data->params;
//    CTreeNode **c_tree = cull_data->c_tree;
//    Boolean isBlastn = (params->program == eBlastTypeBlastn);
//    if (!hsp_list) return 0;
//
//    for (i=0; i<hsp_list->hspcnt; ++i) {
//       /* wrap the hsp with a LinkedHSP structure */
//       A.hsp   = hsp_list->hsp_array[i];
//       A.cid   = isBlastn ? (A.hsp->context  - A.hsp->context % NUM_STRANDS) : A.hsp->context;
//       A.sid   = hsp_list->oid;
//       A.merit = params->culling_max;
//       qlen    = cull_data->query_info->contexts[A.hsp->context].query_length;
//       if(isBlastn && (A.hsp->context % NUM_STRANDS)) {
//     	  A.begin = qlen - A.hsp->query.end;
//     	  A.end   = qlen -  A.hsp->query.offset;
//       }
//       else {
//     	  A.begin = A.hsp->query.offset;
//     	  A.end   = A.hsp->query.end;
//       }
//       A.next  = NULL;
//
//       if (! c_tree[A.cid]) {
//          c_tree[A.cid] = s_CTreeNew(qlen);
//       }
//
//       if(s_SaveHSP(c_tree[A.cid], &A)){
//     	 hsp_list->hsp_array[i] = NULL;
//       }
//    }
//
//    /* now all good hits have moved to tree, we can remove hsp_list */
//    Blast_HSPListFree(hsp_list);
//
//    return 0;
// }
//
// /** Free the writer
//  * @param writer The writer to free [in]
// ```
pub(crate) struct Culling {
    trees: Vec<Option<Box<Node>>>,
    lengths: Vec<i32>,
    limit: i32,
    next_id: usize,
}
impl Culling {
    pub fn new(lengths: Vec<i32>, limit: i32) -> Self {
        Self {
            trees: (0..lengths.len()).map(|_| None).collect(),
            lengths,
            limit,
            next_id: 0,
        }
    }
    pub fn write(&mut self, oid: usize, hsps: Vec<Hsp>, trace: &mut dyn ResultsObserver) {
        for hsp in hsps {
            let context = hsp.hsp.context;
            let entry = Entry {
                id: self.next_id,
                oid,
                hsp,
                merit: self.limit,
            };
            self.next_id += 1;
            trace.cull_save(oid, &entry.hsp);
            let saved = self.trees[context]
                .get_or_insert_with(|| Node::new(0, self.lengths[context]))
                .save(entry, trace);
            trace.cull_saved(saved);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:517-589
    // ```c++
    //    for (cid=0; cid < cull_data->num_contexts; ++cid) {
    //       if (c_tree[cid]) {
    //          qid = Blast_GetQueryIndexFromContext(cid, params->program);
    //          if (!results->hitlist_array[qid]) {
    //             results->hitlist_array[qid] = Blast_HitListNew(params->prelim_hitlist_size);
    //          }
    //          hitlist = results->hitlist_array[qid];
    //
    //          /* collapse the linked hsps tree into one list and free the tree */
    //          cull_list = s_RipHSPOffCTree(c_tree[cid]);
    //          c_tree[cid] = s_CTreeFree(c_tree[cid]);
    //
    //          /* insert hsp list into results */
    //          while (cull_list) {
    //             p = cull_list;
    //             /* test to see if new hsplist has already been allocated */
    //             allocated = FALSE;
    //             for (sid=0; sid<hitlist->hsplist_count; ++sid) {
    //                 list = hitlist->hsplist_array[sid];
    //                 if (p->sid == list->oid) {
    //                    allocated = TRUE;
    //                    break;
    //                 }
    //             }
    //             if (!allocated) {
    //                 /* we must allocate a new hsplist*/
    //                 list = Blast_HSPListNew(0);
    //                 list->oid = p->sid;
    //                 list->query_index = qid;
    //                 if (sid >= hitlist->hsplist_current) {
    //                    /* we must increase the pool size as well */
    //                    new_allocated = MAX(kStartValue, 2*sid);
    //                    hitlist->hsplist_array = (BlastHSPList **)
    //                       realloc(hitlist->hsplist_array, new_allocated*sizeof(BlastHSPList*));
    //                    hitlist->hsplist_current = new_allocated;
    //                 }
    //                 hitlist->hsplist_array[sid] = list;
    //                 hitlist->hsplist_count++;
    //             }
    //             /* put the new hsp into the array */
    //             id = list->hspcnt;
    //             if (id >= list->allocated) {
    //                 /* we must increase the list size */
    //                 new_allocated = 2*id;
    //                 list->hsp_array = (BlastHSP**)
    //                       realloc(list->hsp_array, new_allocated*sizeof(BlastHSP*));
    //                 list->allocated = new_allocated;
    //             }
    //             p = cull_list;
    //             list->hsp_array[id] = p->hsp;
    //             list->hspcnt++;
    //             cull_list = p->next;
    //             free(p);
    //          }
    //
    //          /* sort hsplist */
    //          worst_evalue = 0.0;
    //          low_score = INT4_MAX;
    //          for (sid=0; sid < hitlist->hsplist_count; ++sid) {
    //             list = hitlist->hsplist_array[sid];
    //             best_evalue = (double) INT4_MAX;
    //             for (id=0; id < list->hspcnt; ++id) {
    //                 best_evalue = MIN(list->hsp_array[id]->evalue, best_evalue);
    //             }
    //             Blast_HSPListSortByScore(list);
    //             list->best_evalue = best_evalue;
    //             worst_evalue = MAX(worst_evalue, best_evalue);
    //             low_score = MIN(list->hsp_array[0]->score, low_score);
    //          }
    //          hitlist->worst_evalue = worst_evalue;
    //          hitlist->low_score = low_score;
    //       }
    //    }
    // ```
    pub fn finalize(&mut self, num_queries: usize, max: usize) -> Vec<HitList> {
        let mut queries: Vec<_> = (0..num_queries).map(|_| HitList::new(max)).collect();
        for (context, tree) in self.trees.iter_mut().enumerate() {
            let Some(tree) = tree.take() else {
                continue;
            };
            let mut entries = Vec::new();
            tree.rip(&mut entries);
            let query_index = context / 6;
            let q = &mut queries[query_index];
            for entry in entries {
                if let Some(list) = q.lists.iter_mut().find(|l| l.oid == entry.oid) {
                    list.hsps.push(entry.hsp);
                } else {
                    q.lists.push(HspList {
                        query_index,
                        oid: entry.oid,
                        hsps: vec![entry.hsp],
                        best_evalue: 0.0,
                    });
                }
            }
            for list in &mut q.lists {
                list.hsps.sort_by(|a, b| compare_score(&a.hsp, &b.hsp));
                list.refresh();
            }
        }
        queries
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:722-750
// ```c++
//          results->hitlist_array[qid]->hsplist_array[sid] = NULL;
//       }
//       results->hitlist_array[qid]->hsplist_count = 0;
//       Blast_HitListFree(results->hitlist_array[qid]);
//       results->hitlist_array[qid] = NULL;
//    }
//    s_BlastHSPCullingFinal(data, results);
//    return 0;
// }
//
// /** Free the pipe
//  * @param pipe The pipe to free [in]
//  * @return NULL.
//  */
// static
// BlastHSPPipe*
// s_BlastHSPCullingPipeFree(BlastHSPPipe* pipe)
// {
//    BlastHSPCullingData *data = pipe->data;
//    sfree(data->params);
//    sfree(pipe->data);
//    sfree(pipe);
//    return NULL;
// }
//
// /** create the pipe
//  * @param params Pointer to the besthit parameter [in]
//  * @param query_info BlastQueryInfo [in]
//  * @return pipe
// ```
pub(crate) fn traceback_pipe(
    queries: &mut Vec<HitList>,
    lengths: Vec<i32>,
    limit: i32,
    max: usize,
    trace: &mut dyn ResultsObserver,
) {
    let mut culling = Culling::new(lengths, limit);
    for query in queries.iter_mut() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:709-712
        // ```c++
        //         	 hsp_list = results->hitlist_array[qid]->hsplist_array[sid];
        //         	 Blast_HSPListSortByEvalue(hsp_list);
        //         	 hsp_list->best_evalue = hsp_list->hsp_array[0]->evalue;
        //          }
        // ```
        for list in &mut query.lists {
            list.sort_evalue();
            if let Some(h) = list.hsps.first() {
                list.best_evalue = h.evalue;
            }
        }
        query.lists.sort_by(super::results::compare_list);
    }
    for query in queries.iter_mut() {
        for list in std::mem::take(&mut query.lists) {
            culling.write(list.oid, list.hsps, trace);
        }
    }
    *queries = culling.finalize(queries.len(), max);
    trace.cull_final(queries.len());
    for q in queries {
        for l in &q.lists {
            trace.hsps("CULL_FINAL", l.oid, &l.hsps);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120
// ```c++
// static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
//     Int8 b1 = p->begin;
//     Int8 b2 = y->begin;
//     Int8 e1 = p->end;
//     Int8 e2 = y->end;
//     Int8 s1 = p->hsp->score;
//     Int8 s2 = y->hsp->score;
//     Int8 l1 = e1 - b1;
//     Int8 l2 = e2 - b2;
//     Int8 overlap = MIN(e1,e2) - MAX(b1,b2);
//     Int8 d = 0;
//
//     // If not overlap by more than 50%
//     if(2 *overlap < l2) {
//     	return FALSE;
//     }
//
//     /* the main criterion:
//        2 * (%diff in score) + 1 * (%diff in length) */
//     //Int8 d  = 3*s1*l1 + s1*l2 - s2*l1 - 3*s2*l2;
//     d  = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
//     // If identical, use oid as tie breaker
//     if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
//     	if(s1 != s2) {
//     		return (s1>s2);
//     	}
//     	if(p->sid != y->sid) {
//     		return (p->sid < y->sid);
//     	}
//
//     	if(p->hsp->subject.offset > y->hsp->subject.offset) {
//     		return FALSE;
//     	}
//     	return TRUE;
//     }
//
//    	if (d < 0) {
//    		return FALSE;
//     }
//
//     return TRUE;
// }
// ```
#[cfg(test)]
mod boundary_tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:53-62
    // ```c++
    //     BlastHSP * hsp;
    //     Int4 cid;    /* context id for hsp */
    //     Int4 sid;    /* OID for hsp*/
    //     Int4 begin;  /* query offset in plus strand */
    //     Int4 end;    /* query end in plus strand */
    //     Int4 merit;  /* how many other hsps in the tree dominates me? */
    //     struct LinkedHSP *next;
    // } LinkedHSP;
    //
    // /** functions to manipulate LinkedHSPs */
    // ```
    fn entry(row: &[i32]) -> Entry {
        Entry {
            id: 0,
            oid: row[3] as usize,
            merit: 0,
            hsp: Hsp {
                hsp: super::super::preliminary::PreliminaryHsp {
                    context: 0,
                    frame: 1,
                    score: row[2],
                    q_start: row[0],
                    q_end: row[1],
                    q_gapped_start: row[0],
                    s_start: row[4],
                    s_end: row[4] + 1,
                    s_gapped_start: row[4],
                },
                num: 1,
                evalue: 0.0,
                bit_score: 0.0,
                identity: 0,
                positive: 0,
                edit_script: Vec::new(),
                composition_method: 0,
            },
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120
    // ```c++
    // static Boolean s_DominateTest(LinkedHSP *p, LinkedHSP *y) {
    //     Int8 b1 = p->begin;
    //     Int8 b2 = y->begin;
    //     Int8 e1 = p->end;
    //     Int8 e2 = y->end;
    //     Int8 s1 = p->hsp->score;
    //     Int8 s2 = y->hsp->score;
    //     Int8 l1 = e1 - b1;
    //     Int8 l2 = e2 - b2;
    //     Int8 overlap = MIN(e1,e2) - MAX(b1,b2);
    //     Int8 d = 0;
    //
    //     // If not overlap by more than 50%
    //     if(2 *overlap < l2) {
    //     	return FALSE;
    //     }
    //
    //     /* the main criterion:
    //        2 * (%diff in score) + 1 * (%diff in length) */
    //     //Int8 d  = 3*s1*l1 + s1*l2 - s2*l1 - 3*s2*l2;
    //     d  = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2;
    //     // If identical, use oid as tie breaker
    //     if(((s1 == s2) && (b1==b2) && (l1 == l2)) || (d == 0)) {
    //     	if(s1 != s2) {
    //     		return (s1>s2);
    //     	}
    //     	if(p->sid != y->sid) {
    //     		return (p->sid < y->sid);
    //     	}
    //
    //     	if(p->hsp->subject.offset > y->hsp->subject.offset) {
    //     		return FALSE;
    //     	}
    //     	return TRUE;
    //     }
    //
    //    	if (d < 0) {
    //    		return FALSE;
    //     }
    //
    //     return TRUE;
    // }
    // ```
    #[test]
    fn pinned_c_dominance_boundaries() {
        let golden = include_str!("../../../tests/unit/blastx_d_boundary_expected.tsv");
        let mut count = 0;
        for line in golden.lines().filter(|l| l.starts_with("C\t")) {
            let row: Vec<i32> = line
                .split('\t')
                .skip(1)
                .map(|s| s.parse().unwrap())
                .collect();
            assert_eq!(
                dominates_source(&entry(&row[..5]), &entry(&row[5..10])),
                row[10] != 0,
                "{line}"
            );
            count += 1;
        }
        assert_eq!(count, 12);
    }
}
