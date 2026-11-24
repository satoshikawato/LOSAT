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
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:56,417-422
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:417
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

    /// The heap as `common::database::redo_fetches` reads it at a match's turn, `before` telling
    /// whether a record's subject came before the match in the redo.
    pub(super) fn redo_heap(
        &self,
        before: impl Fn(i32) -> bool,
    ) -> crate::algorithm::common::database::RedoHeap {
        crate::algorithm::common::database::RedoHeap {
            ecutoff: self.ecutoff,
            heap_threshold: self.threshold,
            inserted_before: self
                .records
                .iter()
                .filter(|record| record.best_evalue <= self.ecutoff && before(record.subject_index))
                .count(),
            ends_early: true,
        }
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
#[path = "kappa_heap_tests.rs"]
mod tests;
