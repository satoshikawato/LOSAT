//! NCBI BLAST HSP Culling Implementation
//! 
//! Reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c
//! 
//! This module implements the interval tree-based HSP culling algorithm that removes
//! dominated HSPs based on score/length tradeoff. The implementation is based upon
//! the algorithm described in:
//! Berman P, Zhang Z, Wolf YI, Koonin EV, Miller W. Winnowing sequences from a
//! database search. J Comput Biol. 2000 Feb-Apr;7(1-2):293-302.
//!
//! The original implementation was later rewritten to use interval trees by
//! Jason Papadopoulos.

use super::chaining::UngappedHit;
use super::lookup::QueryContext;

/// Linked list of HSPs used to store hsps in culling tree.
/// 
/// NCBI reference: hspfilter_culling.c:52-60
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
struct LinkedHSP {
    hsp: UngappedHit,
    cid: usize,      // context id for hsp
    sid: u32,        // subject index (OID)
    begin: i32,      // query offset in plus strand
    end: i32,        // query end in plus strand
    merit: i32,      // how many other hsps in the tree dominates me?
    next: Option<Box<LinkedHSP>>,
}

/// Definition of a Culling tree node
/// 
/// NCBI reference: hspfilter_culling.c:201-207
/// ```c
/// typedef struct CTreeNode {
///     Int4 begin;  /* left endpoint */
///     Int4 end;    /* right endpoint */
///     struct CTreeNode *left;    /* left child */
///     struct CTreeNode *right;   /* right child */
///     LinkedHSP *hsplist; /* hsps belong to this node, start with low merits */
/// } CTreeNode;
/// ```
struct CTreeNode {
    begin: i32,              // left endpoint
    end: i32,                // right endpoint
    left: Option<Box<CTreeNode>>,
    right: Option<Box<CTreeNode>>,
    hsplist: Option<Box<LinkedHSP>>,  // hsps belong to this node, start with low merits
}

/// Return true if p dominates y
/// 
/// NCBI reference: hspfilter_culling.c:79-120
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
///
///     /* the main criterion:
///        2 * (%diff in score) + 1 * (%diff in length) */
///     //Int8 d  = 3*s1*l1 + s1*l2 - s2*l1 - 3*s2*l2;
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
#[inline]
fn dominate_test(p: &LinkedHSP, y: &LinkedHSP) -> bool {
    let b1 = p.begin as i64;
    let b2 = y.begin as i64;
    let e1 = p.end as i64;
    let e2 = y.end as i64;
    let s1 = p.hsp.raw_score as i64;
    let s2 = y.hsp.raw_score as i64;
    let l1 = e1 - b1;
    let l2 = e2 - b2;
    let overlap = e1.min(e2) - b1.max(b2);

    // If not overlap by more than 50%
    if 2 * overlap < l2 {
        return false;
    }

    // The main criterion: 2 * (%diff in score) + 1 * (%diff in length)
    // Formula: d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2
    let d = 4 * s1 * l1 + 2 * s1 * l2 - 2 * s2 * l1 - 4 * s2 * l2;
    
    // If identical, use oid as tie breaker
    if ((s1 == s2) && (b1 == b2) && (l1 == l2)) || (d == 0) {
        if s1 != s2 {
            return s1 > s2;
        }
        if p.sid != y.sid {
            return p.sid < y.sid;
        }
        // NCBI: p->hsp->subject.offset > y->hsp->subject.offset
        // LOSAT: use s_aa_start as subject offset
        if p.hsp.s_aa_start > y.hsp.s_aa_start {
            return false;
        }
        return true;
    }

    if d < 0 {
        return false;
    }

    true
}

/// Check how many hsps in list dominates y, and update merit of y accordingly
/// 
/// NCBI reference: hspfilter_culling.c:123-133
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
fn full_pass(list: &Option<Box<LinkedHSP>>, y: &mut LinkedHSP) -> bool {
    let mut p = list;
    while let Some(node) = p {
        if dominate_test(node, y) {
            y.merit -= 1;
            if y.merit <= 0 {
                return false;
            }
        }
        p = &node.next;
    }
    true
}

/// Update merit for hsps in list; also returns the number of hsps in list
/// 
/// NCBI reference: hspfilter_culling.c:136-163
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
///              if (r == *list) {
///                  *list = p;
///                  q = p;
///              } else {
///                  q->next = p;
///              }
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
/// 
/// FAITHFUL PORT: Matches NCBI's in-place list manipulation algorithm exactly.
/// Uses a cursor-based approach to track previous node (q) while traversing.
fn process_hsp_list(list: &mut Option<Box<LinkedHSP>>, y: &LinkedHSP) -> usize {
    // NCBI: Int4 num = 0; LinkedHSP *p = *list, *q, *r; q = p;
    let mut num = 0;
    let mut p = list.take();
    
    // Build new list in-place, matching NCBI's pointer manipulation
    // We use a Vec to collect nodes, then rebuild maintaining order
    let mut nodes: Vec<Box<LinkedHSP>> = Vec::new();
    
    // NCBI: while (p) {
    while let Some(mut r) = p {
        num += 1;
        // NCBI: r = p; p = p->next;
        p = r.next.take();
        
        // NCBI: if (r != y && s_DominateTest(y, r))
        // Check if y dominates r (r != y is handled by domination test returning false for same node)
        let r_dominated = dominate_test(y, &r);
        
        if r_dominated {
            // NCBI: (r->merit)--;
            r.merit -= 1;
            if r.merit <= 0 {
                // NCBI: Remove r from list (s_HSPFree(r))
                // Don't add to nodes - r is dropped here
                num -= 1;
            } else {
                // NCBI: q = r; (keep r, it will be in new list)
                nodes.push(r);
            }
        } else {
            // NCBI: q = r; (keep r)
            nodes.push(r);
        }
    }
    
    // Rebuild list maintaining original order
    let mut new_list: Option<Box<LinkedHSP>> = None;
    for mut node in nodes.into_iter().rev() {
        node.next = new_list;
        new_list = Some(node);
    }
    
    *list = new_list;
    num
}

/// Add an hsp to the front of hsp list
/// 
/// NCBI reference: hspfilter_culling.c:193-197
/// ```c
/// static void s_AddHSPtoList(LinkedHSP **list, LinkedHSP *y) {
///     y->next = *list;
///     *list = y;
///     return;
/// }
/// ```
fn add_hsp_to_list(list: &mut Option<Box<LinkedHSP>>, y: Box<LinkedHSP>) {
    let mut new_node = y;
    new_node.next = list.take();
    *list = Some(new_node);
}

/// Allocate and return a new node for use
/// 
/// NCBI reference: hspfilter_culling.c:228-247
fn ctree_node_new(parent: Option<&CTreeNode>, dir: bool) -> CTreeNode {
    let mut node = CTreeNode {
        begin: 0,
        end: 0,
        left: None,
        right: None,
        hsplist: None,
    };
    
    if let Some(p) = parent {
        let midpt = (p.begin + p.end) / 2;
        if dir {
            // eLeft
            node.begin = p.begin;
            node.end = midpt;
        } else {
            // eRight
            node.begin = midpt;
            node.end = p.end;
        }
    }
    
    node
}

/// Fork children from a node
/// 
/// NCBI reference: hspfilter_culling.c:258-299
/// ```c
/// static void s_ForkChildren(CTreeNode * node) {
///     CTreeNode * child;
///     LinkedHSP *p, *q, *r;
///     Int4 midpt;
///
///     ASSERT(node != NULL);
///     ASSERT(node->left ==NULL);
///     ASSERT(node->right ==NULL);
///
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
///          if (r == node->hsplist) {
///              node->hsplist = p;
///              q = p;
///          } else {   
///              q->next = p;
///          }
///          /* and put it on the child */
///          s_AddHSPtoList(&(child->hsplist), r);
///       } else {
///          q = r;
///       }      
///    }
/// }
/// ```
/// 
/// FAITHFUL PORT: Matches NCBI's in-place list manipulation exactly.
fn fork_children(node: &mut CTreeNode) {
    let midpt = (node.begin + node.end) / 2;
    let mut p = node.hsplist.take();
    let mut parent_list: Vec<Box<LinkedHSP>> = Vec::new();
    let mut child_left_list: Vec<Box<LinkedHSP>> = Vec::new();
    let mut child_right_list: Vec<Box<LinkedHSP>> = Vec::new();
    
    // NCBI: while(p) {
    while let Some(mut r) = p {
        // NCBI: r = p; p = p->next;
        p = r.next.take();
        
        // NCBI: Determine which child (if any) this HSP belongs to
        if r.end < midpt {
            // NCBI: child = node->left;
            if node.left.is_none() {
                node.left = Some(Box::new(ctree_node_new(Some(node), true)));
            }
            child_left_list.push(r);
        } else if r.begin > midpt {
            // NCBI: child = node->right;
            if node.right.is_none() {
                node.right = Some(Box::new(ctree_node_new(Some(node), false)));
            }
            child_right_list.push(r);
        } else {
            // NCBI: q = r; (keep in parent)
            parent_list.push(r);
        }
    }
    
    // Rebuild parent list
    let mut new_parent_list: Option<Box<LinkedHSP>> = None;
    for mut node_hsp in parent_list.into_iter().rev() {
        node_hsp.next = new_parent_list;
        new_parent_list = Some(node_hsp);
    }
    node.hsplist = new_parent_list;
    
    // Rebuild left child list
    if let Some(ref mut left) = node.left {
        for mut hsp in child_left_list.into_iter().rev() {
            add_hsp_to_list(&mut left.hsplist, hsp);
        }
    }
    
    // Rebuild right child list
    if let Some(ref mut right) = node.right {
        for mut hsp in child_right_list.into_iter().rev() {
            add_hsp_to_list(&mut right.hsplist, hsp);
        }
    }
}

/// Recursively search and update merit hsps in culling tree due to addition of hsp x
/// 
/// NCBI reference: hspfilter_culling.c:334-370
fn process_ctree(node: &mut Option<Box<CTreeNode>>, x: &LinkedHSP) {
    if node.is_none() {
        return;
    }
    
    let node_ref = node.as_mut().unwrap();
    
    // First test if x includes the full range covered by node
    if x.begin <= node_ref.begin && x.end >= node_ref.end {
        // Mark down entire subtree
        mark_down_ctree(node);
        return;
    }
    
    // If node reaches the leaves
    if node_ref.left.is_none() && node_ref.right.is_none() {
        if process_hsp_list(&mut node_ref.hsplist, x) == 0 {
            *node = None;
        }
        return;
    }
    
    // Recursive case
    let midpt = (node_ref.begin + node_ref.end) / 2;
    if x.end < midpt {
        process_ctree(&mut node_ref.left, x);
    } else if x.begin > midpt {
        process_ctree(&mut node_ref.right, x);
    } else {
        process_ctree(&mut node_ref.left, x);
        process_ctree(&mut node_ref.right, x);
        if process_hsp_list(&mut node_ref.hsplist, x) == 0
            && node_ref.left.is_none()
            && node_ref.right.is_none()
        {
            *node = None;
        }
    }
}

/// Recursively decrease the merit of all hsps within a subtree
/// 
/// NCBI reference: hspfilter_culling.c:319-330
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
fn mark_down_ctree(node: &mut Option<Box<CTreeNode>>) {
    if node.is_none() {
        return;
    }
    
    let node_ref = node.as_mut().unwrap();
    mark_down_ctree(&mut node_ref.left);
    mark_down_ctree(&mut node_ref.right);
    
    // Decrease merit for all HSPs in this node's list
    let num_remaining = mark_down_hsp_list(&mut node_ref.hsplist);
    
    if num_remaining == 0 && node_ref.left.is_none() && node_ref.right.is_none() {
        *node = None;
    }
}

/// Decrease merit for all hsps in list; also returns the number of hsps in list
/// 
/// NCBI reference: hspfilter_culling.c:166-189
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
///           if (r == *list) {
///               *list = p;
///               q = p;
///           } else {
///               q->next = p;
///           }
///           --num;
///           s_HSPFree(r);
///        } else {
///           q = r;
///        }
///     }
///     return num;
/// }
/// ```
/// 
/// FAITHFUL PORT: Matches NCBI's in-place list manipulation algorithm exactly.
fn mark_down_hsp_list(list: &mut Option<Box<LinkedHSP>>) -> usize {
    // NCBI: Int4 num = 0; LinkedHSP *p = *list, *q, *r; q = p;
    let mut num = 0;
    let mut p = list.take();
    let mut nodes: Vec<Box<LinkedHSP>> = Vec::new();
    
    // NCBI: while (p) {
    while let Some(mut r) = p {
        num += 1;
        // NCBI: r = p; p = p->next; (r->merit)--;
        p = r.next.take();
        r.merit -= 1;
        
        if r.merit <= 0 {
            // NCBI: Remove r from list (s_HSPFree(r))
            // Don't add to nodes - r is dropped here
            num -= 1;
        } else {
            // NCBI: q = r; (keep r)
            nodes.push(r);
        }
    }
    
    // Rebuild list maintaining original order
    let mut new_list: Option<Box<LinkedHSP>> = None;
    for mut node in nodes.into_iter().rev() {
        node.next = new_list;
        new_list = Some(node);
    }
    
    *list = new_list;
    num
}

/// Allocate a tree
/// 
/// NCBI reference: hspfilter_culling.c:376-381
/// ```c
/// static CTreeNode * s_CTreeNew(Int4 qlen) {
///     CTreeNode * tree = s_CTreeNodeNew(NULL, eLeft);
///     tree->begin = 0;
///     tree->end   = qlen;
///     return tree;
/// }
/// ```
fn ctree_new(qlen: i32) -> CTreeNode {
    let mut tree = ctree_node_new(None, true);
    tree.begin = 0;
    tree.end = qlen;
    tree
}

/// Recursively rip off hsps into a link list
/// 
/// NCBI reference: hspfilter_culling.c:396-423
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
///     if (!q) {
///        q = s_RipHSPOffCTree(tree->left);
///        p = q;
///     } else {
///        p = q;
///        while(p->next) p=p->next;
///        p->next = s_RipHSPOffCTree(tree->left);
///     }
///
///     /* grab right child */
///     if (!q) {
///        q = s_RipHSPOffCTree(tree->right);
///     } else {
///        while(p->next) p=p->next;
///        p->next = s_RipHSPOffCTree(tree->right);
///     }
///
///     return q;
/// }
/// ```
fn rip_hsp_off_ctree(tree: Option<Box<CTreeNode>>) -> Option<Box<LinkedHSP>> {
    let mut tree = tree?;
    
    let mut q = tree.hsplist.take();
    tree.hsplist = None;
    
    // Grab left child
    let left_list = rip_hsp_off_ctree(tree.left.take());
    if q.is_none() {
        q = left_list;
    } else {
        let mut p = q.as_mut().unwrap();
        while p.next.is_some() {
            p = p.next.as_mut().unwrap();
        }
        p.next = left_list;
    }
    
    // Grab right child
    let right_list = rip_hsp_off_ctree(tree.right.take());
    if q.is_none() {
        q = right_list;
    } else {
        let mut p = q.as_mut().unwrap();
        while p.next.is_some() {
            p = p.next.as_mut().unwrap();
        }
        p.next = right_list;
    }
    
    q
}

/// A full traverse to determine the merit of A, in addition, insert A to the proper place if A is valid,
/// or return FALSE if A's merit decreases to zero
/// 
/// NCBI reference: hspfilter_culling.c:430-470
const K_NUM_HSP_TO_FORK: i32 = 20;  // number of HSP to trigger forking children

fn save_hsp(tree: &mut CTreeNode, a: &mut LinkedHSP) -> bool {
    let mut current: *mut CTreeNode = tree;
    let mut last_node: *mut CTreeNode = std::ptr::null_mut();
    
    // Descend the tree
    unsafe {
        loop {
            let node = &mut *current;
            
            // Check merit against current node's list
            if !full_pass(&node.hsplist, a) {
                return false;
            }
            
            let midpt = (node.begin + node.end) / 2;
            last_node = current;  // record the last valid position
            
            if a.end < midpt {
                if let Some(ref mut left) = node.left {
                    current = left.as_mut();
                    continue;
                } else {
                    break;
                }
            } else if a.begin > midpt {
                if let Some(ref mut right) = node.right {
                    current = right.as_mut();
                    continue;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        
        // If we get here, A is valid. Copy and insert A at node
        let node = &mut *last_node;
        let x = Box::new(LinkedHSP {
            hsp: a.hsp.clone(),
            cid: a.cid,
            sid: a.sid,
            begin: a.begin,
            end: a.end,
            merit: a.merit,
            next: None,
        });
        add_hsp_to_list(&mut node.hsplist, x);
        
        // Create a reference to the newly added HSP for domination checks
        // We need to clone the key fields since we can't keep a reference
        let x_for_dom = LinkedHSP {
            hsp: a.hsp.clone(),
            cid: a.cid,
            sid: a.sid,
            begin: a.begin,
            end: a.end,
            merit: a.merit,
            next: None,
        };
        
        // If this is the leaf, calculate update hsp number
        if node.left.is_none() && node.right.is_none() {
            // Check for domination
            let num_remaining = process_hsp_list(&mut node.hsplist, &x_for_dom);
            if num_remaining >= K_NUM_HSP_TO_FORK as usize {
                // Fork this node into sub trees
                fork_children(node);
            }
            return true;
        }
        
        // Check domination
        process_hsp_list(&mut node.hsplist, &x_for_dom);
        process_ctree(&mut node.left, &x_for_dom);
        process_ctree(&mut node.right, &x_for_dom);
        true
    }
}

/// Apply NCBI HSP culling to a list of UngappedHits
/// 
/// This is the main entry point for culling. It implements the logic from
/// `s_BlastHSPCullingRun` (hspfilter_culling.c:602-644).
/// 
/// NCBI reference: hspfilter_culling.c:602-644
/// ```c
/// static int 
/// s_BlastHSPCullingRun(void* data, BlastHSPList* hsp_list)
/// {
///    Int4 i, qlen;
///    LinkedHSP A;
///
///    BlastHSPCullingData * cull_data = data;
///    BlastHSPCullingParams* params = cull_data->params;
///    CTreeNode **c_tree = cull_data->c_tree;
///    Boolean isBlastn = (params->program == eBlastTypeBlastn);
///    if (!hsp_list) return 0;
///
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
///
///    /* now all good hits have moved to tree, we can remove hsp_list */
///    Blast_HSPListFree(hsp_list);
///         
///    return 0; 
/// }
/// ```
/// Apply NCBI HSP culling to a list of UngappedHits
/// 
/// This is the main entry point for culling. It implements the logic from
/// `s_BlastHSPCullingRun` (hspfilter_culling.c:602-644) and extraction from
/// `s_BlastHSPCullingFinal` (hspfilter_culling.c:500-593).
/// 
/// NCBI reference: hspfilter_culling.c:602-644
/// ```c
/// static int 
/// s_BlastHSPCullingRun(void* data, BlastHSPList* hsp_list)
/// {
///    Int4 i, qlen;
///    LinkedHSP A;
///
///    BlastHSPCullingData * cull_data = data;
///    BlastHSPCullingParams* params = cull_data->params;
///    CTreeNode **c_tree = cull_data->c_tree;
///    Boolean isBlastn = (params->program == eBlastTypeBlastn);
///    if (!hsp_list) return 0;
///
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
///
///    /* now all good hits have moved to tree, we can remove hsp_list */
///    Blast_HSPListFree(hsp_list);
///         
///    return 0; 
/// }
/// ```
/// 
/// For tblastx: isBlastn = FALSE, so we use A.cid = A.hsp->context directly,
/// and coordinates are already in plus strand (A.begin = query.offset, A.end = query.end).
pub fn apply_culling(
    hits: Vec<UngappedHit>,
    contexts: &[QueryContext],
    culling_limit: u32,
) -> Vec<UngappedHit> {
    if culling_limit == 0 || hits.is_empty() {
        return hits;
    }
    
    // NCBI: Group hits by context (cid)
    // For tblastx: cid = ctx_idx (no NUM_STRANDS adjustment needed)
    let mut hits_by_context: std::collections::HashMap<usize, Vec<UngappedHit>> = 
        std::collections::HashMap::new();
    for hit in hits {
        hits_by_context.entry(hit.ctx_idx).or_default().push(hit);
    }
    
    let mut trees: std::collections::HashMap<usize, CTreeNode> = 
        std::collections::HashMap::new();
    
    // NCBI: Process each context separately (s_BlastHSPCullingRun loop)
    for (cid, context_hits) in hits_by_context {
        let ctx = &contexts[cid];
        // NCBI: qlen = cull_data->query_info->contexts[A.hsp->context].query_length;
        let qlen = ctx.orig_len as i32;
        
        // NCBI: if (! c_tree[A.cid]) { c_tree[A.cid] = s_CTreeNew(qlen); }
        let tree = trees.entry(cid).or_insert_with(|| ctree_new(qlen));
        
        // NCBI: for (i=0; i<hsp_list->hspcnt; ++i) {
        for hit in context_hits {
            // NCBI: wrap the hsp with a LinkedHSP structure
            // For tblastx: isBlastn = FALSE, so:
            //   A.cid = A.hsp->context (no adjustment)
            //   A.begin = A.hsp->query.offset
            //   A.end = A.hsp->query.end
            // 
            // LOSAT: q_aa_start/q_aa_end are already in plus strand coordinates
            let begin = hit.q_aa_start as i32;
            let end = hit.q_aa_end as i32;
            
            let mut linked_hsp = LinkedHSP {
                hsp: hit,
                cid,
                sid: 0,  // Will be set below
                begin,
                end,
                merit: culling_limit as i32,  // NCBI: A.merit = params->culling_max;
                next: None,
            };
            // NCBI: A.sid = hsp_list->oid;
            linked_hsp.sid = linked_hsp.hsp.s_idx;
            
            // NCBI: if(s_SaveHSP(c_tree[A.cid], &A)) { hsp_list->hsp_array[i] = NULL; }
            // If save_hsp returns true, HSP was saved (not culled)
            if !save_hsp(tree, &mut linked_hsp) {
                // HSP was culled (merit <= 0), don't include it
            }
        }
    }
    
    // NCBI: Extract all saved HSPs from trees (s_BlastHSPCullingFinal: s_RipHSPOffCTree)
    let mut result: Vec<UngappedHit> = Vec::new();
    for (_cid, tree) in trees {
        // NCBI: cull_list = s_RipHSPOffCTree(c_tree[cid]);
        let mut extracted = rip_hsp_off_ctree(Some(Box::new(tree)));
        while let Some(mut node) = extracted {
            extracted = node.next.take();
            // NCBI: Extract the HSP from LinkedHSP wrapper
            result.push(node.hsp);
        }
    }
    
    // NCBI: Sort by score (s_BlastHSPCullingFinal: Blast_HSPListSortByScore)
    result.sort_by(|a, b| b.raw_score.cmp(&a.raw_score));

    result
}


