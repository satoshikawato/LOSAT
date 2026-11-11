//! `SegSet` against NCBI's list-based `CSegSetTmpl`: the results of `AddSeg`, `Merge`, `Cross`
//! and `Clip` (with `m_iFactor` 1) that the functions of RpsbProc/src/common/impl/segset_impl.hpp
//! gave, copied into a small C++ program outside the repository (`segsetgen.cpp` of the task
//! folder of item PS16a; g++ 13), on 400 random pairs of sets
//! (`LOSAT/tests/unit/rpsbproc_segset_vectors.tsv`: the operation, the two sets and the result,
//! each range `from:to:ori_from`).

use super::{Seg, SegSet};

fn parse(text: &str) -> Vec<Seg> {
    text.split(',')
        .filter(|range| !range.is_empty())
        .map(|range| {
            let mut fields = range.split(':').map(|field| field.parse::<i32>().unwrap());
            Seg {
                from: fields.next().unwrap(),
                to: fields.next().unwrap(),
                ori_from: fields.next().unwrap(),
            }
        })
        .collect()
}

#[test]
fn segset_operations_give_ncbis_ranges() {
    let table = include_str!("../../tests/unit/rpsbproc_segset_vectors.tsv");
    let mut count = 0;
    for line in table.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let (op, a, b, expected) = (fields[0], fields[1], fields[2], fields[3]);
        let mut set = SegSet::from_segs(parse(a));
        let other = SegSet::from_segs(parse(b));
        match op {
            "merge" => set.merge(&other),
            "cross" => set.cross(&other),
            "clip" => set.clip(&other),
            _ => {
                for seg in other.segs() {
                    set.add(seg.from, seg.to);
                }
            }
        }
        assert_eq!(set.segs(), parse(expected).as_slice(), "{line}");
        count += 1;
    }
    assert_eq!(count, 400);
}

// The overlap of two hits of one region (`Calculate`, cdalignproc_base.cpp:206-244): a range
// with a gap clipped out, crossed with another, and the total.
#[test]
fn clipped_range_crossed_with_another() {
    let mut this = SegSet::new();
    this.add(10, 100);
    let mut gaps = SegSet::new();
    gaps.add(40, 79);
    this.clip(&gaps);
    assert_eq!(
        this.segs(),
        &[
            Seg::new(10, 39),
            Seg {
                from: 80,
                to: 100,
                ori_from: 80
            }
        ]
    );
    assert_eq!(this.total_residues(), 51);
    let mut other = SegSet::new();
    other.add(30, 90);
    let mut overlap = this.clone();
    overlap.cross(&other);
    assert_eq!(overlap.total_residues(), 10 + 11);
    // An invalid range adds nothing.
    let mut empty = SegSet::new();
    empty.add(5, 4);
    assert!(empty.is_empty());
}
