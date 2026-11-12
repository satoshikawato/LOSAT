//! The pieces of the domain classification: `CleanAlignment`, `CalcMasterGaps`, `CountSrc`.

use super::*;
use crate::rpsbproc::archive::{Align, DenseSeg, Segs};

fn dense_seg(starts: &[i32], lens: &[u32]) -> DenseSeg {
    DenseSeg {
        ids: Vec::new(),
        starts: starts.to_vec(),
        lens: lens.to_vec(),
    }
}

fn align() -> Align {
    Align {
        segs: Segs::Other("dendiag"),
        score: Some(10),
        e_value: Some(1e-5),
        bit_score: Some(40.5),
        num_ident: None,
        score_of_other_type: None,
    }
}

// impl/basedata_impl.hpp:2-19 and biodata_blast.cpp:23-43 (quoted at `clean_alignment` and
// `parse_align`): the segments with a gap in either row are dropped; the range runs from the
// first kept segment to the end of the last one.
#[test]
fn parse_align_keeps_the_aligned_segments() {
    let seg = dense_seg(&[4, 0, -1, 70, 74, 80, 90, -1], &[70, 5, 10, 3]);
    let dom = parse_align(&seg, &align(), 42);
    assert_eq!(dom.m_starts, vec![4, 74]);
    assert_eq!(dom.s_starts, vec![0, 80]);
    assert_eq!(dom.lens, vec![70, 10]);
    assert_eq!((dom.from, dom.to), (4, 83));
    assert_eq!((dom.e_value, dom.bit_score, dom.pssm_id), (1e-5, 40.5, 42));
    // No aligned segment: NCBI's range stays `k_INT_MAX` .. 0.
    let none = parse_align(&dense_seg(&[-1, 0], &[5]), &align(), 42);
    assert_eq!((none.from, none.to), (i32::MAX, 0));
    assert!(none.s_starts.is_empty());
}

// basedata.cpp:1324-1339 (quoted at `master_gaps`): query gaps of at least 35 residues.
#[test]
fn master_gaps_of_35_residues_and_more() {
    let seg = dense_seg(&[0, 0, 44, 10, 100, 20, 134, 30], &[10, 10, 10, 10]);
    let dom = parse_align(&seg, &align(), 1);
    let gaps = master_gaps(&dom);
    // 10..43 (34 residues) is not a gap; 54..99 (46) and 110..133 (24): only the first.
    assert_eq!(gaps.total_residues(), 46);
}

// basedata.cpp:154-166 (quoted at `count_source`): one hit per source database and region;
// the PRK family (PRK, CHL, MTH, PHA, PLN, PTZ) counts as one source; unknown sources always
// count.
#[test]
fn one_hit_per_source_database() {
    let mut counts = [0; 7];
    assert!(count_source(&mut counts, b"cd00001"));
    assert!(!count_source(&mut counts, b"cd00002"));
    assert!(count_source(&mut counts, b"pfam00001"));
    assert!(count_source(&mut counts, b"PRK00001"));
    assert!(!count_source(&mut counts, b"PTZ00001"));
    assert!(count_source(&mut counts, b"NF000001"));
    assert!(count_source(&mut counts, b"NF000002"));
}
