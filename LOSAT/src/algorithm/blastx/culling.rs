//! BLASTX context-local culling writer and traceback pipe.
use super::{
    preliminary::compare_score,
    results::{HitList, Hsp, HspList, ResultsObserver},
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:52-60
// `merit` is an `i64`, not NCBI's `Int4`, which overflows for a limit near `INT4_MAX` plus 3
// (NCBI parity DECISIONS 17: LOSAT does not reproduce the overflow).
#[derive(Clone)]
struct Entry {
    id: usize,
    oid: usize,
    hsp: Hsp,
    merit: i64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120
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
fn dominates(a: &Entry, b: &Entry, trace: &mut dyn ResultsObserver) -> bool {
    let result = dominates_source(a, b);
    trace.cull_compare(a.oid, &a.hsp, b.oid, &b.hsp, result);
    result
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:201-249
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
    #[expect(
        clippy::boxed_local,
        reason = "the recursion hands the children over as boxes; a by-value `self` would copy a whole node into every frame of a tree that can be deep"
    )]
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
pub(crate) struct Culling {
    trees: Vec<Option<Box<Node>>>,
    lengths: Vec<i32>,
    limit: i64,
    next_id: usize,
}
impl Culling {
    pub fn new(lengths: Vec<i32>, limit: i64) -> Self {
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
pub(crate) fn traceback_pipe(
    queries: &mut Vec<HitList>,
    lengths: Vec<i32>,
    limit: i32,
    max: usize,
    trace: &mut dyn ResultsObserver,
) {
    let mut culling = Culling::new(lengths, i64::from(limit));
    for query in queries.iter_mut() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:709-712
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
#[cfg(test)]
#[path = "culling_boundary_tests.rs"]
mod boundary_tests;
