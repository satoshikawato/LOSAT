//! NCBI composition-adjusted TBLASTN result heap.
//! This remains internal while the public TBLASTN pipeline is gated.

use anyhow::{ensure, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:75-101
#[derive(Clone, Copy, Debug)]
pub(super) struct CompoHeapRecord {
    pub best_evalue: f64,
    pub best_score: i32,
    pub subject_index: i32,
}

impl CompoHeapRecord {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:89-101
    fn worse_than(self, other: Self) -> bool {
        if self.best_evalue > other.best_evalue {
            return true;
        }
        if self.best_evalue < other.best_evalue {
            return false;
        }
        if other.best_score > self.best_score {
            return true;
        }
        if other.best_score < self.best_score {
            return false;
        }
        other.subject_index > self.subject_index
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:43-54,401-422
#[allow(dead_code)] // Attached to the TBLASTN D result pipeline after all boundaries agree.
pub(super) struct CompoHeap {
    records: Vec<CompoHeapRecord>,
    threshold: usize,
    ecutoff: f64,
    worst_evalue: f64,
    capacity: usize,
    heapified: bool,
}

#[allow(dead_code)] // Source-equivalent heap is exercised by comparison fixtures.
impl CompoHeap {
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:401-419
    pub(super) fn new(threshold: usize, ecutoff: f64) -> Result<Self> {
        ensure!(
            threshold > 0,
            "NCBI composition heap threshold must be positive"
        );
        Ok(Self {
            records: Vec::new(),
            threshold,
            ecutoff,
            worst_evalue: 0.0,
            capacity: threshold.min(100),
            heapified: false,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:401-415
    pub(super) fn len(&self) -> usize {
        self.records.len()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:421
    pub(super) fn capacity(&self) -> usize {
        self.capacity
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:422
    pub(super) fn worst_evalue(&self) -> f64 {
        self.worst_evalue
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:405-410:
    pub(super) fn filled_to_cutoff(&self) -> bool {
        self.records.len() >= self.threshold && self.worst_evalue <= self.ecutoff
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:163-211,237-245
    fn heapify_down(&mut self, start: usize) {
        let mut position = start;
        loop {
            let left = position * 2 + 1;
            let right = left + 1;
            let mut largest = position;
            if left < self.records.len() && self.records[left].worse_than(self.records[position]) {
                largest = left;
            }
            if right < self.records.len() && self.records[right].worse_than(self.records[largest]) {
                largest = right;
            }
            if largest == position {
                break;
            }
            self.records.swap(position, largest);
            position = largest;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:212-245
    fn heapify_up(&mut self, mut position: usize) {
        while position > 0 {
            let parent = (position - 1) / 2;
            if !self.records[position].worse_than(self.records[parent]) {
                break;
            }
            self.records.swap(position, parent);
            position = parent;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:225-245
    fn convert_to_heap(&mut self) {
        if !self.heapified {
            self.heapified = true;
            for position in (0..self.records.len() / 2).rev() {
                self.heapify_down(position);
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    pub(super) fn would_insert(&mut self, candidate: CompoHeapRecord) -> bool {
        if self.records.len() < self.threshold
            || candidate.best_evalue <= self.ecutoff
            || candidate.best_evalue < self.worst_evalue
        {
            return true;
        }
        self.convert_to_heap();
        self.records[0].worse_than(candidate)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:292-391
    pub(super) fn insert(&mut self, candidate: CompoHeapRecord) -> Option<CompoHeapRecord> {
        if !self.heapified && self.records.len() >= self.threshold {
            self.convert_to_heap();
        }
        if !self.heapified {
            if self.records.len() >= self.capacity {
                self.capacity = (self.capacity + 100).max((1.5 * self.capacity as f64) as usize);
            }
            self.records.push(candidate);
            if self.worst_evalue < candidate.best_evalue {
                self.worst_evalue = candidate.best_evalue;
            }
            return None;
        }
        let discarded = if self.records.len() < self.threshold
            || (candidate.best_evalue <= self.ecutoff && self.worst_evalue <= self.ecutoff)
        {
            if self.records.len() >= self.capacity {
                self.capacity = (self.capacity + 100).max((1.5 * self.capacity as f64) as usize);
            }
            self.records.push(candidate);
            self.heapify_up(self.records.len() - 1);
            None
        } else if self.records[0].worse_than(candidate) {
            let previous = std::mem::replace(&mut self.records[0], candidate);
            self.heapify_down(0);
            Some(previous)
        } else {
            self.heapify_down(0);
            Some(candidate)
        };
        self.worst_evalue = self.records[0].best_evalue;
        discarded
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:439-466
    pub(super) fn pop(&mut self) -> Option<CompoHeapRecord> {
        self.convert_to_heap();
        if self.records.is_empty() {
            return None;
        }
        let result = self.records.swap_remove(0);
        if !self.records.is_empty() {
            self.heapify_down(0);
        }
        Some(result)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:64,1565-1566,1576-1579,
// 1560-1580:
//         if (evalue <= EVALUE_STRETCH * significantMatches[i].ecutoff)
#[allow(dead_code)] // Used by TBLASTN Kappa before each redo in the gated D pipeline.
pub(super) fn compo_early_termination(evalue: f64, heaps: &[CompoHeap]) -> bool {
    for heap in heaps {
        if !heap.filled_to_cutoff() || evalue <= 5.0 * heap.ecutoff {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:64,
    // 1560-1580; composition_adjustment/compo_heap.c:405-410:
    #[test]
    fn natural_112_subject_early_termination_matches_ncbi_call_state() {
        let trace = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/natural_c_d_early_20260925/early.tsv"
        ))
        .unwrap();
        let mut redone = 0;
        let mut terminated = 0;
        for row in trace
            .lines()
            .filter(|line| line.starts_with("K_EARLY_EVAL\t"))
        {
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(f[1], "0");
            let evalue = f[2].parse().unwrap();
            let n = f[3].parse().unwrap();
            let threshold = f[4].parse().unwrap();
            let worst_evalue = f[5].parse().unwrap();
            let ecutoff = f[6].parse().unwrap();
            let expected = f[7] == "1";
            let mut heap = CompoHeap::new(threshold, ecutoff).unwrap();
            heap.records = vec![
                CompoHeapRecord {
                    best_evalue: worst_evalue,
                    best_score: 0,
                    subject_index: 0,
                };
                n
            ];
            heap.worst_evalue = worst_evalue;
            assert_eq!(compo_early_termination(evalue, &[heap]), expected);
            if expected {
                terminated += 1;
            } else {
                redone += 1;
            }
        }
        assert_eq!((redone, terminated), (11, 16));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3705-3736,2500-2515;
    // composition_adjustment/compo_heap.c:252-275,330-391,439-466
    #[test]
    fn saved_ncbi_local_subject_heap_rank_and_pop_order() {
        for (path, threshold, expected) in [
            (
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../docs/evidence/tlosan_stage_d/kappa_traceback_20260924/run_20260923_default.tsv"
                ),
                500,
                (11, 11, 11, 0),
            ),
            (
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/run_20260925/ncbi.trace"
                ),
                2,
                (11, 10, 10, 1),
            ),
        ] {
            let trace = std::fs::read_to_string(path).unwrap();
            let mut heap = CompoHeap::new(threshold, 0.002).unwrap();
            let mut decisions = 0usize;
            let mut inserts = 0usize;
            let mut pops = 0usize;
            let mut rejected = 0usize;
            for line in trace.lines() {
                let f: Vec<_> = line.split('\t').collect();
                match f[0] {
                    "K_TRACE_HEAP_WOULD" => {
                        assert_eq!(heap.len(), f[4].parse().unwrap(), "{path}");
                        assert_eq!(heap.capacity(), f[6].parse().unwrap(), "{path}");
                        assert_eq!(
                            heap.worst_evalue().to_bits(),
                            f[8].parse::<f64>().unwrap().to_bits(),
                            "{path}"
                        );
                        let candidate = CompoHeapRecord {
                            best_evalue: f[2].parse().unwrap(),
                            best_score: f[3].parse().unwrap(),
                            subject_index: f[1].parse().unwrap(),
                        };
                        let accepted = heap.would_insert(candidate);
                        assert_eq!(accepted, f[9] == "1", "{path}");
                        if !accepted {
                            rejected += 1;
                        }
                        decisions += 1;
                    }
                    "K_TRACE_HEAP_INSERT" => {
                        let candidate = CompoHeapRecord {
                            best_evalue: f[2].parse().unwrap(),
                            best_score: f[3].parse().unwrap(),
                            subject_index: f[1].parse().unwrap(),
                        };
                        assert_eq!(heap.insert(candidate).map(|item| item.subject_index), None);
                        inserts += 1;
                    }
                    "K_TRACE_HEAP_INSERT_RETURN" => {
                        assert_eq!(heap.len(), f[3].parse().unwrap(), "{path}");
                        assert_eq!(
                            heap.worst_evalue().to_bits(),
                            f[4].parse::<f64>().unwrap().to_bits(),
                            "{path}"
                        );
                        assert_eq!(f[5], "-1");
                    }
                    "K_TRACE_HEAP_POP" => {
                        let actual = heap.pop();
                        assert_eq!(
                            actual.map_or(-1, |item| item.subject_index),
                            f[1].parse().unwrap(),
                            "{path}"
                        );
                        assert_eq!(heap.len(), f[2].parse().unwrap(), "{path}");
                        if let Some(item) = actual {
                            assert_eq!(
                                item.best_evalue.to_bits(),
                                f[3].parse::<f64>().unwrap().to_bits(),
                                "{path}"
                            );
                            pops += 1;
                        }
                    }
                    _ => {}
                }
            }
            assert_eq!((decisions, inserts, pops, rejected), expected, "{path}");
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:89-101,252-275,330-391
    #[test]
    fn heap_threshold_ties_and_discarded_record_follow_ncbi_comparator() {
        let mut heap = CompoHeap::new(2, 0.001).unwrap();
        let rec = |oid, e, score| CompoHeapRecord {
            subject_index: oid,
            best_evalue: e,
            best_score: score,
        };
        heap.insert(rec(4, 0.1, 30));
        heap.insert(rec(3, 0.1, 30));
        assert!(!heap.would_insert(rec(2, 0.1, 30)));
        assert!(heap.would_insert(rec(5, 0.1, 30)));
        assert_eq!(heap.insert(rec(5, 0.1, 30)).unwrap().subject_index, 3);
        assert_eq!(heap.pop().unwrap().subject_index, 4);
        assert_eq!(heap.pop().unwrap().subject_index, 5);
    }
}
