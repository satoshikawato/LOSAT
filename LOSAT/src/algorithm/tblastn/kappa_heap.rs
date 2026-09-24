//! NCBI composition-adjusted TBLASTN result heap.
//! This remains internal while the public TBLASTN pipeline is gated.

use anyhow::{ensure, Result};

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:75-101
// ```c
// if (0 == (result = CMP(place1->bestEvalue, place2->bestEvalue)) &&
//     0 == (result = CMP(place2->bestScore, place1->bestScore))) {
//     result = CMP(place2->subject_index, place1->subject_index);
// }
// return result > 0;
// ```
#[derive(Clone, Copy, Debug)]
pub(super) struct CompoHeapRecord {
    pub best_evalue: f64,
    pub best_score: i32,
    pub subject_index: i32,
}

impl CompoHeapRecord {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:89-101
    // ```c
    // if (0 == (result = CMP(place1->bestEvalue, place2->bestEvalue)) &&
    //     0 == (result = CMP(place2->bestScore, place1->bestScore)))
    //     result = CMP(place2->subject_index, place1->subject_index);
    // return result > 0;
    // ```
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:43-54,401-419
// ```c
// #define HEAP_INITIAL_CAPACITY 100
// self->heapThreshold = heapThreshold;
// self->ecutoff = ecutoff;
// self->capacity = MIN(HEAP_INITIAL_CAPACITY, heapThreshold);
// self->worstEvalue = 0;
// ```
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:401-419
    // ```c
    // self->n = 0; self->heapThreshold = heapThreshold;
    // self->ecutoff = ecutoff; self->capacity = MIN(100, heapThreshold);
    // self->worstEvalue = 0;
    // ```
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:401-415
    // ```c
    // self->n = 0; self->capacity = MIN(HEAP_INITIAL_CAPACITY, heapThreshold);
    // self->worstEvalue = 0;
    // ```
    pub(super) fn len(&self) -> usize {
        self.records.len()
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:411-415
    // ```c
    // self->capacity = MIN(HEAP_INITIAL_CAPACITY, heapThreshold);
    // ```
    pub(super) fn capacity(&self) -> usize {
        self.capacity
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:411-414
    // ```c
    // self->worstEvalue = 0;
    // ```
    pub(super) fn worst_evalue(&self) -> f64 {
        self.worst_evalue
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:163-211,237-245
    // ```c
    // if (left <= n && s_CompoHeapRecordCompare(&heapArray[left], &heapArray[i]))
    //     largest = left;
    // if (right <= n && s_CompoHeapRecordCompare(&heapArray[right], &heapArray[largest]))
    //     largest = right;
    // if (largest != i) s_CompoHeapRecordSwap(&heapArray[i], &heapArray[largest]);
    // for (i = n / 2; i >= 1; --i) s_CompoHeapifyDown(self->heapArray, i, n);
    // ```
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:212-245
    // ```c
    // while (parent >= 1 && s_CompoHeapRecordCompare(&heapArray[i], &heapArray[parent])) {
    //     s_CompoHeapRecordSwap(&heapArray[i], &heapArray[parent]);
    //     i = parent; parent /= 2;
    // }
    // ```
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:225-245
    // ```c
    // self->heapArray = self->array;
    // self->array = NULL;
    // for (i = n / 2; i >= 1; --i) s_CompoHeapifyDown(self->heapArray, i, n);
    // ```
    fn convert_to_heap(&mut self) {
        if !self.heapified {
            self.heapified = true;
            for position in (0..self.records.len() / 2).rev() {
                self.heapify_down(position);
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:252-275
    // ```c
    // if (self->n < self->heapThreshold || eValue <= self->ecutoff ||
    //     eValue < self->worstEvalue) return TRUE;
    // if (self->heapArray == NULL) s_ConvertToHeap(self);
    // return s_CompoHeapRecordCompare(&self->heapArray[1], &heapRecord);
    // ```
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:292-391
    // ```c
    // new_capacity = MAX(HEAP_MIN_RESIZE + *capacity,
    //                    (int)(HEAP_RESIZE_FACTOR * (*capacity)));
    // if (self->array && self->n >= self->heapThreshold) s_ConvertToHeap(self);
    // if (self->array != NULL) { append; worstEvalue = MAX(worstEvalue, eValue); }
    // else if (n < threshold || (eValue <= ecutoff && worstEvalue <= ecutoff))
    //     { append; s_CompoHeapifyUp(...); }
    // else { compare root with candidate; discard one; s_CompoHeapifyDown(...); }
    // ```
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:439-466
    // ```c
    // s_ConvertToHeap(self);
    // results = first->theseAlignments;
    // if (--self->n > 0) { memcpy(first, last, sizeof(BlastCompo_HeapRecord));
    //     s_CompoHeapifyDown(self->heapArray, 1, self->n); }
    // ```
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

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3705-3736,2500-2515;
    // composition_adjustment/compo_heap.c:252-275,330-391,439-466
    // ```c
    // if (BlastCompo_HeapWouldInsert(...))
    //     BlastCompo_HeapInsert(..., &discarded_aligns);
    // while (NULL != (hsp_list = BlastCompo_HeapPop(heap)))
    //     Blast_HitListUpdate(hitlist, hsp_list);
    // ```
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/composition_adjustment/compo_heap.c:89-101,252-275,330-391
    // ```c
    // if (self->n < self->heapThreshold || eValue <= self->ecutoff) return TRUE;
    // else return s_CompoHeapRecordCompare(&self->heapArray[1], &heapRecord);
    // ```
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
