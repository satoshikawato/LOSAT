use super::*;

/// A small deterministic generator (xorshift), so the tests need no crate.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> i32 {
        (self.next() % n) as i32
    }
}

fn record(rng: &mut Rng) -> XUnitHsp {
    // Few distinct values, so that ties on every key are common.
    XUnitHsp {
        q_start: rng.below(6),
        s_start: rng.below(4),
        len: 1 + rng.below(3),
        q_seed: rng.below(1000),
        s_seed: rng.below(1000),
        score: rng.below(3),
        ctx_idx: rng.below(2) as u32,
    }
}

fn init(r: &XUnitHsp) -> InitHSP {
    InitHSP {
        q_start_absolute: r.q_start,
        q_end_absolute: r.q_start + r.len,
        s_start: r.s_start,
        s_end: r.s_start + r.len,
        q_seed_absolute: r.q_seed,
        s_seed: r.s_seed,
        score: r.score,
        ctx_idx: r.ctx_idx as usize,
        s_f_idx: 0,
        q_idx: 0,
        s_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_orig_len: 100,
        s_orig_len: 100,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-304
// The record's comparator is the init comparator on the expanded records.
#[test]
fn the_record_order_is_the_init_order() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..20_000 {
        let (a, b) = (record(&mut rng), record(&mut rng));
        assert_eq!(
            <XUnitHsp as XRecord>::cmp(&a, &b),
            score_compare_init_hsps_ncbi(&init(&a), &init(&b)),
            "{a:?} {b:?}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
// Argument (a): the k-way merge of the units' sorted lists is the stable sort of their
// concatenation, record for record (the seed offsets tell equal records apart), with empty
// lists among them.
#[test]
fn the_kway_merge_is_the_stable_sort_of_the_concatenation() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    for case in 0..2_000 {
        let k = 1 + (case % 6);
        let lists: Vec<Vec<InitHSP>> = (0..k)
            .map(|_| {
                let n = rng.below(12) as usize;
                (0..n).map(|_| init(&record(&mut rng))).collect()
            })
            .collect();
        let mut expected: Vec<InitHSP> = lists.concat();
        expected.sort_by(score_compare_init_hsps_ncbi);
        let mut sorted = lists.clone();
        sorted.iter_mut().for_each(XUnitSort::x_sort);
        let unit_free = |lists: Vec<Vec<InitHSP>>| {
            // `InitHSP::init_hsp` does not read the unit; any `XExpand` of an empty subject.
            let frames: [QueryFrame; 0] = [];
            let inp = XSubject {
                s_frames_preliminary: &frames,
                s_frames: &frames,
                contexts: &[],
                cutoff_scores: &[],
                x_dropoff_per_context: &[],
                window: 40,
                wordsize: 3,
                s_idx: 0,
                s_len: 0,
                failure: &crate::algorithm::tblastx::subjects::XReadFailure::new(),
                passed: None,
                subject_besthit: false,
            };
            let unit = XExpand {
                f: 0,
                s_frame: 1,
                s_idx: 0,
                inp: &inp,
            };
            x_kway(lists, &unit)
        };
        assert_eq!(unit_free(sorted), expected, "case {case}");
    }
}

// No NCBI counterpart: the merged list's capacity (module documentation, "Memory"). A list of
// 1 MiB or more gets glibc's mapped capacity (on 64-bit glibc; elsewhere its length) and the
// same elements as the stable sort.
#[test]
fn the_kway_merge_of_a_large_frame_has_the_mapped_capacity() {
    let mut rng = Rng(0x5851_f42d_4c95_7f2d);
    let n = (1usize << 20).div_ceil(std::mem::size_of::<InitHSP>()) / 2 + 1;
    let lists: Vec<Vec<InitHSP>> = (0..2)
        .map(|_| {
            let mut list: Vec<InitHSP> = (0..n).map(|_| init(&record(&mut rng))).collect();
            list.sort_by(score_compare_init_hsps_ncbi);
            list
        })
        .collect();
    let mut expected: Vec<InitHSP> = lists.concat();
    expected.sort_by(score_compare_init_hsps_ncbi);
    let frames: [QueryFrame; 0] = [];
    let inp = XSubject {
        s_frames_preliminary: &frames,
        s_frames: &frames,
        contexts: &[],
        cutoff_scores: &[],
        x_dropoff_per_context: &[],
        window: 40,
        wordsize: 3,
        s_idx: 0,
        s_len: 0,
        failure: &crate::algorithm::tblastx::subjects::XReadFailure::new(),
        passed: None,
        subject_besthit: false,
    };
    let unit = XExpand {
        f: 0,
        s_frame: 1,
        s_idx: 0,
        inp: &inp,
    };
    let merged = x_kway(lists, &unit);
    assert_eq!(merged, expected);
    let bytes = merged.capacity() * std::mem::size_of::<InitHSP>();
    if crate::utils::x_mapped::X_MAPPED_ALLOCATOR {
        assert!(bytes > 32 << 20, "{bytes}");
    } else {
        assert_eq!(merged.capacity(), 2 * n);
    }
    assert_eq!(x_frame_capacity::<InitHSP>(3), 3);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
// The unit's sort in its task is `slice::sort_by` (stable), record for record, at every length
// around the run and pass boundaries, with a buffer of any capacity and contents.
#[test]
fn the_buffered_sort_is_the_stable_sort() {
    let mut rng = Rng(0x0bad_5eed_1234_5678);
    for n in (0..70).chain([127, 128, 129, 255, 256, 257, 1000, 4099]) {
        for _ in 0..4 {
            // seed offsets tell equal records apart
            let list: Vec<XUnitHsp> = (0..n)
                .map(|i| XUnitHsp {
                    q_seed: i,
                    ..record(&mut rng)
                })
                .collect();
            let mut expected = list.clone();
            expected.sort_by(<XUnitHsp as XRecord>::cmp);
            let mut sorted = list.clone();
            let mut buf: Vec<XUnitHsp> = (0..rng.below(8)).map(|_| record(&mut rng)).collect();
            x_stable_sort_with(&mut sorted, &mut buf, <XUnitHsp as XRecord>::cmp);
            assert_eq!(sorted, expected, "n {n}");
        }
    }
}

fn hit(rng: &mut Rng, order: u32) -> UngappedHit {
    let s = rng.below(4) as u32;
    let q = rng.below(4) as u32;
    let len = 1 + rng.below(2) as u32;
    UngappedHit {
        q_idx: 0,
        s_idx: 0,
        ctx_idx: rng.below(3) as u32,
        s_f_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_aa_start: q,
        q_aa_end: q + len,
        s_aa_start: s,
        s_aa_end: s + len,
        q_seed_off: q,
        s_seed_off: s,
        q_orig_len: 100,
        s_orig_len: 100,
        raw_score: rng.below(3),
        e_value: f64::INFINITY,
        num_ident: 0,
        hsp_list_order: order,
        ordering_method: 0,
        linked_set: false,
        start_of_chain: false,
        link_id: 0,
        chain_next_link_id: None,
        hsp_link_num: 0,
        num: 0,
    }
}

fn same_lists(a: &[UngappedHit], b: &[UngappedHit]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x_same_hit(x, y))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
// Argument (c): on a K1-sorted list the run-wise sort is the full stable sort; on a list out of
// K1 order the full sort runs. `hsp_list_order` tells equal HSPs apart.
#[test]
fn the_run_wise_sort_is_the_full_sort() {
    let mut rng = Rng(0x0123_4567_89ab_cdef);
    for case in 0..2_000 {
        let n = rng.below(40) as u32;
        let mut list: Vec<UngappedHit> = (0..n).map(|i| hit(&mut rng, i)).collect();
        if case % 2 == 0 {
            list.sort_by(x_k1);
        }
        let mut expected = list.clone();
        sort_ungapped_hits_by_score_ncbi(&mut expected);
        x_sort_k1_runs(&mut list);
        assert!(same_lists(&list, &expected), "case {case}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766
// Argument (d): the backward merge is the append and the stable sort, also for an unsorted or
// empty frame list and an empty combined list.
#[test]
fn the_frame_merge_is_the_append_and_the_sort() {
    let mut rng = Rng(0xdead_beef_cafe_f00d);
    for case in 0..2_000 {
        let mut combined: Vec<UngappedHit> = (0..rng.below(30) as u32)
            .map(|i| hit(&mut rng, i))
            .collect();
        sort_ungapped_hits_by_score_ncbi(&mut combined);
        let mut frame: Vec<UngappedHit> = (0..rng.below(30) as u32)
            .map(|i| hit(&mut rng, 100 + i))
            .collect();
        if case % 3 != 0 {
            sort_ungapped_hits_by_score_ncbi(&mut frame);
        }
        let mut expected = combined.clone();
        expected.extend(frame.iter().cloned());
        if !ungapped_hits_is_sorted_by_score_ncbi(&expected) {
            sort_ungapped_hits_by_score_ncbi(&mut expected);
        }
        x_merge_frame(&mut combined, frame);
        assert!(same_lists(&combined, &expected), "case {case}");
    }
}
