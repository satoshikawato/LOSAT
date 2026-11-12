//! Tests of `hit_sort`: the port of libstdc++'s `std::list::sort` against a linked-list model of
//! the same algorithm (nodes moved by `_M_transfer`, as libstdc++ moves them), and NCBI's
//! comparators with their edges.

use super::*;

/// A small deterministic generator (xorshift64*), so the tests need no crate.
struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }
}

/// libstdc++'s `std::list` as a circular doubly linked list of nodes: node 0 of each list is its
/// head (`_M_impl._M_node`), the others hold the elements (`value`).
struct Lists {
    next: Vec<usize>,
    prev: Vec<usize>,
    value: Vec<usize>,
}

impl Lists {
    fn new() -> Self {
        Lists {
            next: Vec::new(),
            prev: Vec::new(),
            value: Vec::new(),
        }
    }

    /// A new empty list: its head node.
    fn head(&mut self) -> usize {
        let node = self.next.len();
        self.next.push(node);
        self.prev.push(node);
        self.value.push(usize::MAX);
        node
    }

    fn push_back(&mut self, head: usize, value: usize) {
        let node = self.next.len();
        self.next.push(head);
        self.prev.push(self.prev[head]);
        self.value.push(value);
        let last = self.prev[head];
        self.next[last] = node;
        self.prev[head] = node;
    }

    fn empty(&self, head: usize) -> bool {
        self.next[head] == head
    }

    /// `position->_M_transfer(first, last)`: moves the nodes `[first, last)` before `position`.
    fn transfer(&mut self, position: usize, first: usize, last: usize) {
        if position == last || first == last {
            return;
        }
        let before_first = self.prev[first];
        let before_last = self.prev[last];
        // Unlink [first, last).
        self.next[before_first] = last;
        self.prev[last] = before_first;
        // Link before `position`.
        let before_position = self.prev[position];
        self.next[before_position] = first;
        self.prev[first] = before_position;
        self.next[before_last] = position;
        self.prev[position] = before_last;
    }

    /// `_List_node_base::swap(x, y)` of two heads: the lists trade their nodes.
    fn swap(&mut self, x: usize, y: usize) {
        let nodes = |lists: &Self, head: usize| {
            (!lists.empty(head)).then(|| (lists.next[head], lists.prev[head]))
        };
        let (x_nodes, y_nodes) = (nodes(self, x), nodes(self, y));
        for (head, held) in [(x, y_nodes), (y, x_nodes)] {
            match held {
                Some((first, last)) => {
                    self.next[head] = first;
                    self.prev[head] = last;
                    self.prev[first] = head;
                    self.next[last] = head;
                }
                None => {
                    self.next[head] = head;
                    self.prev[head] = head;
                }
            }
        }
    }

    /// `_Scratch_list::merge(x, comp)` into `this`.
    fn merge(&mut self, this: usize, x: usize, comp: &mut impl FnMut(usize, usize) -> bool) {
        let mut first1 = self.next[this];
        let last1 = this;
        let mut first2 = self.next[x];
        let last2 = x;
        while first1 != last1 && first2 != last2 {
            if comp(self.value[first2], self.value[first1]) {
                let next = self.next[first2];
                self.transfer(first1, first2, next);
                first2 = next;
            } else {
                first1 = self.next[first1];
            }
        }
        if first2 != last2 {
            self.transfer(last1, first2, last2);
        }
    }

    fn values(&self, head: usize) -> Vec<usize> {
        let mut values = Vec::new();
        let mut node = self.next[head];
        while node != head {
            values.push(self.value[node]);
            node = self.next[node];
        }
        values
    }
}

/// `std::list::sort(comp)` on the linked-list model (bits/list.tcc:614-662).
fn model_sort(values: &[usize], mut comp: impl FnMut(usize, usize) -> bool) -> Vec<usize> {
    let mut lists = Lists::new();
    let list = lists.head();
    for &value in values {
        lists.push_back(list, value);
    }
    if lists.empty(list) || lists.next[lists.next[list]] == list {
        return lists.values(list);
    }
    let carry = lists.head();
    let tmp: Vec<usize> = (0..64).map(|_| lists.head()).collect();
    let mut fill = 0;
    loop {
        // __carry._M_take_one(begin()._M_node);
        let first = lists.next[list];
        let after = lists.next[first];
        lists.transfer(carry, first, after);
        let mut counter = 0;
        while counter != fill && !lists.empty(tmp[counter]) {
            lists.merge(tmp[counter], carry, &mut comp);
            lists.swap(carry, tmp[counter]);
            counter += 1;
        }
        lists.swap(carry, tmp[counter]);
        if counter == fill {
            fill += 1;
        }
        if lists.empty(list) {
            break;
        }
    }
    for counter in 1..fill {
        lists.merge(tmp[counter], tmp[counter - 1], &mut comp);
    }
    lists.values(tmp[fill - 1])
}

/// `list_sort` of the elements `0..keys.len()` by `less(keys[a], keys[b])`, with its calls.
fn port_sort(keys: &[i64], less: fn(i64, i64) -> bool) -> (Vec<usize>, Vec<(usize, usize)>) {
    let mut calls = Vec::new();
    let mut list: Vec<usize> = (0..keys.len()).collect();
    list_sort(&mut list, |&a, &b| {
        calls.push((a, b));
        less(keys[a], keys[b])
    });
    (list, calls)
}

fn model(keys: &[i64], less: fn(i64, i64) -> bool) -> (Vec<usize>, Vec<(usize, usize)>) {
    let mut calls = Vec::new();
    let values: Vec<usize> = (0..keys.len()).collect();
    let sorted = model_sort(&values, |a, b| {
        calls.push((a, b));
        less(keys[a], keys[b])
    });
    (sorted, calls)
}

#[test]
fn the_port_calls_the_comparator_as_the_linked_list_sort_does() {
    let comparators: [(&str, fn(i64, i64) -> bool); 6] = [
        ("<", |a, b| a < b),
        (">", |a, b| a > b),
        ("<=", |a, b| a <= b),
        (">=", |a, b| a >= b),
        ("always", |_, _| true),
        ("never", |_, _| false),
    ];
    let mut random = Random(0x9e37_79b9_7f4a_7c15);
    for case in 0..600 {
        let n = match case % 6 {
            0 => random.below(4) as usize,
            1 => random.below(20) as usize,
            _ => random.below(130) as usize,
        };
        // Few distinct keys: many ties.
        let distinct = 1 + random.below(if case % 2 == 0 { 3 } else { 1000 });
        let keys: Vec<i64> = (0..n).map(|_| random.below(distinct) as i64).collect();
        for (name, less) in comparators {
            let (port, port_calls) = port_sort(&keys, less);
            let (expected, expected_calls) = model(&keys, less);
            assert_eq!(port, expected, "{name} {keys:?}");
            assert_eq!(port_calls, expected_calls, "{name} {keys:?}");
        }
    }
}

#[test]
fn small_lists_sort_with_libstdcxx_s_calls() {
    let always = |_: i64, _: i64| true;
    // 0 or 1 element: no call.
    assert_eq!(port_sort(&[], always), (vec![], vec![]));
    assert_eq!(port_sort(&[7], always), (vec![0], vec![]));
    // The newer element goes first only when comp(newer, older).
    assert_eq!(port_sort(&[1, 1], always), (vec![1, 0], vec![(1, 0)]));
    assert_eq!(port_sort(&[1, 1], |a, b| a < b), (vec![0, 1], vec![(1, 0)]));
    // `>=` reverses equal elements (NCBI's total score with ties: dup9, dup1, msA give msA,
    // dup1, dup9).
    assert_eq!(
        port_sort(&[5, 5, 5], |a, b| a >= b),
        (vec![2, 1, 0], vec![(1, 0), (2, 1)])
    );
    assert_eq!(
        port_sort(&[5, 5, 5, 5], always),
        (vec![3, 2, 1, 0], vec![(1, 0), (3, 2), (3, 1), (2, 1)])
    );
    // A strict order keeps equal elements in their order (stable).
    assert_eq!(port_sort(&[2, 1, 2, 1], |a, b| a > b).0, vec![0, 2, 1, 3]);
    assert_eq!(port_sort(&[3, 1, 2], |a, b| a < b).0, vec![1, 2, 0]);
}

fn hsp(bits: f64, evalue: f64, num_ident: i64, length: i64, query: (u64, u64), subject: u64) -> HspScores {
    HspScores {
        bits,
        evalue,
        num_ident,
        length,
        query_from: query.0,
        query_to: query.1,
        subject_from: subject,
    }
}

#[test]
fn hsp_comparators_break_ties_by_e_value() {
    let a = hsp(50.0, 1e-5, 40, 50, (10, 59), 3);
    let b = hsp(50.0, 1e-6, 40, 50, (10, 59), 3);
    // Equal scores: `>` is false both ways.
    assert!(!hsp_by_score_descending(&a, &b) && !hsp_by_score_descending(&b, &a));
    // Equal starts: the lower E-value first.
    assert!(hsp_by_query_start_ascending(&b, &a) && !hsp_by_query_start_ascending(&a, &b));
    assert!(hsp_by_subject_start_ascending(&b, &a) && !hsp_by_subject_start_ascending(&a, &b));
    let c = hsp(10.0, 1.0, 1, 5, (9, 9), 2);
    assert!(hsp_by_query_start_ascending(&c, &a) && hsp_by_subject_start_ascending(&c, &a));
    // Equal identity ratios (40/50 and 4/5): the E-values.
    let d = hsp(20.0, 1e-7, 4, 5, (0, 4), 0);
    assert!(hsp_by_percent_identity_descending(&d, &a, false, false));
    assert!(!hsp_by_percent_identity_descending(&a, &d, false, false));
    // Equal ratios and E-values: false both ways.
    assert!(!hsp_by_percent_identity_descending(&a, &a, false, false));
    // Different ratios: the higher first (`>=` on unequal ratios).
    let e = hsp(20.0, 1e-9, 3, 5, (0, 4), 0);
    assert!(hsp_by_percent_identity_descending(&a, &e, false, false));
    assert!(!hsp_by_percent_identity_descending(&e, &a, false, false));
    // No identities on one side: the E-values decide, whatever the ratios.
    let none = hsp(20.0, 1e-3, 0, 5, (0, 4), 0);
    assert!(!hsp_by_percent_identity_descending(&none, &a, false, false));
    assert!(hsp_by_percent_identity_descending(&a, &none, false, false));
    let early = hsp(20.0, 1e-30, 0, 5, (0, 4), 0);
    assert!(hsp_by_percent_identity_descending(&early, &a, false, false));
}

#[test]
fn tblastx_lengths_are_nucleotides_unless_translated() {
    let a = hsp(50.0, 1e-5, 10, 30, (0, 89), 0);
    assert_eq!(a.alignment_length(false, false), 30);
    assert_eq!(a.alignment_length(false, true), 30);
    assert_eq!(a.alignment_length(true, false), 90);
    assert_eq!(a.alignment_length(true, true), 30);
    // The identity of the description table (no translation) is a third of the shown one.
    assert_eq!(percent_identity(30, a.alignment_length(true, false)) as i64, 33);
    assert_eq!(percent_identity(30, a.alignment_length(true, true)), 100.0);
    assert_eq!(percent_identity(7, 7), 100.0);
    assert_eq!(percent_identity(0, 1), 0.0);
}

#[test]
fn query_coverage_is_the_union_of_the_query_ranges() {
    let ranges = [(10, 19), (15, 30), (40, 40), (31, 35), (0, 2)];
    let hsps: Vec<HspScores> = ranges.iter().map(|&q| hsp(1.0, 1.0, 1, 1, q, 0)).collect();
    // 0-2 (3), 10-35 (26, adjacent ranges count once each), 40 (1).
    assert_eq!(query_coverage(&hsps), 30);
    assert_eq!(query_coverage(&hsps[..1]), 10);
    assert_eq!(query_coverage(&[]), 0);
    let contained = [hsp(1.0, 1.0, 1, 1, (0, 99), 0), hsp(1.0, 1.0, 1, 1, (10, 20), 0)];
    assert_eq!(query_coverage(&contained), 100);
}

fn sorted(sort: ReportSort, tblastx: bool, hits: &[Vec<HspScores>]) -> SortedHits {
    sort_query_hits(sort, tblastx, hits).expect("a sort")
}

fn by(hits: Option<u8>, hsps: Option<u8>) -> ReportSort {
    ReportSort { hits, hsps }
}

#[test]
fn nothing_is_sorted_without_a_positive_option() {
    let hits = vec![vec![hsp(1.0, 1.0, 1, 1, (0, 0), 0)]; 3];
    for sort in [by(None, None), by(Some(0), None), by(None, Some(0)), by(Some(0), Some(0))] {
        assert!(sort_query_hits(sort, false, &hits).is_none(), "{sort:?}");
        assert!(!sort.sorts());
    }
    assert!(by(Some(0), None).shows_columns() && !by(None, Some(3)).shows_columns());
    assert!(by(Some(0), Some(1)).sorts() && by(Some(4), None).sorts());
}

#[test]
fn equal_total_scores_come_in_reverse() {
    // Three copies of one subject (dup9, dup1, msA in the search's order) and a lower one.
    let copy = vec![hsp(30.0, 1e-3, 20, 25, (0, 24), 0), hsp(10.0, 1.0, 5, 8, (30, 37), 0)];
    let low = vec![hsp(35.0, 1e-4, 20, 25, (0, 24), 0)];
    let hits = vec![copy.clone(), copy.clone(), low, copy];
    let order = sorted(by(Some(2), None), false, &hits);
    // Totals 40, 40, 35, 40: `>=` puts each later 40 before the earlier ones.
    assert_eq!(order.hits, vec![3, 1, 0, 2]);
    // The HSPs keep their order.
    assert_eq!(order.hsps, vec![vec![0, 1], vec![0, 1], vec![0], vec![0, 1]]);
}

#[test]
fn the_score_sort_orders_the_hsps_of_compared_hits_only() {
    // The second HSP has the higher score.
    let two = vec![hsp(20.0, 1e-3, 10, 12, (0, 11), 0), hsp(25.0, 1e-2, 10, 12, (20, 31), 5)];
    // One hit: no comparison, the HSPs stay in their order.
    let order = sorted(by(Some(1), None), false, std::slice::from_ref(&two));
    assert_eq!(order, SortedHits { hits: vec![0], hsps: vec![vec![0, 1]] });
    // Two hits: the comparison sorts both hits' HSPs by score.
    let other = vec![hsp(22.0, 1e-4, 10, 12, (0, 11), 0)];
    let order = sorted(by(Some(1), None), false, &[other.clone(), two.clone()]);
    assert_eq!(order.hits, vec![1, 0]);
    assert_eq!(order.hsps, vec![vec![0], vec![1, 0]]);
    // Then -sorthsps 2 puts them back by query start.
    let order = sorted(by(Some(1), Some(2)), false, &[other, two]);
    assert_eq!(order.hsps, vec![vec![0], vec![0, 1]]);
}

#[test]
fn the_percent_identity_sort_reads_the_best_scoring_hsp() {
    // Hit 0: best HSP 18/20 (90 %); hit 1: best HSP 30/40 (75 %) despite an HSP of 100 %.
    let hit0 = vec![hsp(30.0, 1e-5, 18, 20, (0, 19), 0)];
    let hit1 = vec![hsp(40.0, 1e-8, 30, 40, (10, 49), 0), hsp(5.0, 1.0, 4, 4, (0, 3), 0)];
    let order = sorted(by(Some(3), None), false, &[hit1.clone(), hit0.clone()]);
    assert_eq!(order.hits, vec![1, 0]);
    // The comparison sorted hit 1's HSPs by percent identity (100 % first).
    assert_eq!(order.hsps[0], vec![1, 0]);
    // Equal percentages compare E-values of -1: false, so the order stays.
    let same = vec![hsp(30.0, 1e-9, 18, 20, (0, 19), 0)];
    let order = sorted(by(Some(3), None), false, &[hit0.clone(), same.clone()]);
    assert_eq!(order.hits, vec![0, 1]);
    let order = sorted(by(Some(3), None), false, &[same, hit0]);
    assert_eq!(order.hits, vec![0, 1]);
    // No identities: a percentage of 0 is not above 0, the order stays.
    let zero = vec![hsp(50.0, 1e-20, 0, 20, (0, 19), 0)];
    let best = vec![hsp(10.0, 1.0, 20, 20, (0, 19), 0)];
    let order = sorted(by(Some(3), None), false, &[zero.clone(), best.clone()]);
    assert_eq!(order.hits, vec![0, 1]);
}

#[test]
fn tblastx_s_identity_sort_takes_amino_acid_lengths() {
    // 9/10 against 25/30 amino acids: by amino acids 90 % against 83 %; in nucleotides (no
    // translation) 30 % against 27.8 %, the same order. The HSP sort with -sorthits 3 reads
    // the same lengths (kTranslation), without it nucleotides: here an HSP of 9 identities in 10
    // letters and one of 28 in 30 order the same either way, but the hit sort's percentages
    // differ by the factor three (GetPercentIdentity of 9 and 30 is 30 %).
    let a = vec![hsp(30.0, 1e-5, 9, 10, (0, 29), 0)];
    let b = vec![hsp(40.0, 1e-7, 25, 30, (0, 89), 0)];
    assert_eq!(sorted(by(Some(3), None), true, &[b.clone(), a.clone()]).hits, vec![1, 0]);
    assert_eq!(set_percent_identity(&a, &[0], true, true), 90.0);
    assert_eq!(set_percent_identity(&a, &[0], true, false), 30.0);
    // An HSP whose identities equal its amino-acid length is 100 % only with the translation.
    let full = vec![hsp(30.0, 1e-5, 10, 10, (0, 29), 0)];
    assert_eq!(set_percent_identity(&full, &[0], true, true), 100.0);
    assert!((set_percent_identity(&full, &[0], true, false) - 100.0 / 3.0).abs() < 1e-12);
}

#[test]
fn the_coverage_sort_breaks_ties_by_the_first_hsp_s_e_value() {
    let wide = vec![hsp(20.0, 1e-2, 10, 50, (0, 49), 0)];
    let narrow = vec![hsp(60.0, 1e-20, 10, 20, (0, 19), 0)];
    let tie_low = vec![hsp(20.0, 1e-3, 10, 50, (100, 149), 0)];
    let order = sorted(by(Some(4), None), false, &[narrow, wide, tie_low]);
    // 50, 50 (the lower E-value first), then 20.
    assert_eq!(order.hits, vec![2, 1, 0]);
}

#[test]
fn hsp_sorts_apply_to_every_hit() {
    let hit = vec![
        hsp(10.0, 1e-1, 5, 10, (50, 59), 7),
        hsp(30.0, 1e-5, 9, 10, (5, 14), 90),
        hsp(20.0, 1e-3, 8, 10, (5, 14), 3),
    ];
    let hits = vec![hit.clone(), hit];
    assert_eq!(sorted(by(None, Some(1)), false, &hits).hsps, vec![vec![1, 2, 0]; 2]);
    // Equal query starts: the lower E-value.
    assert_eq!(sorted(by(None, Some(2)), false, &hits).hsps, vec![vec![1, 2, 0]; 2]);
    assert_eq!(sorted(by(None, Some(3)), false, &hits).hsps, vec![vec![1, 2, 0]; 2]);
    assert_eq!(sorted(by(Some(0), Some(4)), false, &hits).hsps, vec![vec![2, 0, 1]; 2]);
    // The hits keep their order.
    assert_eq!(sorted(by(None, Some(4)), false, &hits).hits, vec![0, 1]);
}
