use super::*;
use crate::common::Hit;

fn hsp(oid: u32, frame: i32, begin: usize, end: usize, score: i32) -> TblastxHsp {
    TblastxHsp {
        hit: Hit {
            identity: 0.0,
            length: end - begin,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 3,
            s_start: 1,
            s_end: 3,
            e_value: 1e-5,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: frame,
            query_length: 0,
            q_idx: 0,
            s_idx: oid,
            raw_score: score,
            sort_query_offset: begin,
            sort_query_end: end,
            sort_subject_offset: begin,
            sort_subject_end: end,
            has_sort_offsets: true,
            gap_info: None,
            num_positives: 0,
        },
        subject_frame: 1,
        num: 1,
    }
}

fn kept(writer: CullingWriter) -> Vec<(u32, i32)> {
    writer
        .finalize()
        .into_iter()
        .flat_map(|(oid, list)| list.into_iter().map(move |h| (oid, h.hit.raw_score)))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:922-929,1211-1220: frames 1, 2, 3, -1,
// -2, -3 are the contexts 0..5, of (length - context % 3) / 3 residues.
#[test]
fn contexts_follow_ncbis_frames() {
    let order: Vec<usize> = [1, 2, 3, -1, -2, -3]
        .iter()
        .map(|&frame| tblastx_context(frame))
        .collect();
    assert_eq!(order, vec![0, 1, 2, 3, 4, 5]);
    assert_eq!(tblastx_context_lengths(31), [10, 10, 9, 10, 10, 9]);
    assert_eq!(tblastx_context_lengths(0), [0; 6]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:136-163: an HSP never dominates
// itself (`r != y`), so a limit of 1 keeps the best HSP of overlapping ones.
#[test]
fn a_limit_of_one_keeps_the_dominating_hsp() {
    let mut writer = CullingWriter::new(1, tblastx_context_lengths(300));
    writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
    assert_eq!(kept(writer), vec![(0, 50)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120: with equal score and range
// the lower OID dominates; other contexts have their own trees.
#[test]
fn ties_are_broken_by_oid_and_contexts_do_not_interact() {
    let mut writer = CullingWriter::new(1, tblastx_context_lengths(300));
    writer.run(vec![hsp(1, 1, 10, 40, 50)], 1);
    writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, -1, 10, 40, 50)], 0);
    assert_eq!(kept(writer), vec![(0, 50), (0, 50)]);
}

// The merit does not wrap as NCBI's Int4 does (NCBI parity DECISIONS 17): the limits at
// the top of the range, plus 3, keep a dominated HSP.
#[test]
fn merits_near_int4_max_do_not_wrap() {
    let lengths = tblastx_context_lengths(300);
    for limit in [2147483645, 2147483646, i32::MAX] {
        let mut writer =
            CullingWriter::new(CullingForest::<TblastxHsp>::prelim_max(limit), lengths);
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
        assert_eq!(kept(writer), vec![(0, 50), (0, 30)]);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:258-299,453-458: a leaf with 20
// HSPs forks; HSPs that do not overlap are all kept.
#[test]
fn many_disjoint_hsps_fork_and_survive() {
    let mut writer = CullingWriter::new(1, tblastx_context_lengths(3000));
    let list: Vec<TblastxHsp> = (0..40)
        .map(|i| hsp(0, 1, i * 20, i * 20 + 10, 30))
        .collect();
    writer.run(list, 0);
    assert_eq!(kept(writer).len(), 40);
}
