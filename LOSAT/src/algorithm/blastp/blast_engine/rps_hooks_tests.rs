use super::*;
use crate::algorithm::blastp::hsp::{BlastpHsp, BlastpHspList};

/// A list of one subject whose first Seq-align has the given E-value and bit score.
fn list(oid: u32, e_value: f64, bit_score: f64) -> BlastpHspList {
    let mut hsp = BlastpHsp::from_hit(crate::common::Hit {
        identity: 0.0,
        length: 1,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 1,
        s_start: 1,
        s_end: 1,
        e_value,
        bit_score,
        num_ident: 0,
        query_frame: 0,
        query_length: 0,
        q_idx: 0,
        s_idx: oid,
        raw_score: 10,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 0,
    });
    hsp.e_value = e_value;
    hsp.bit_score = bit_score;
    BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![hsp],
        best_evalue: e_value,
    }
}

fn hit_list(lists: Vec<BlastpHspList>) -> BlastpHitList {
    let mut hit_list = BlastpHitList::new(100);
    hit_list.hsplist_count = lists.len();
    hit_list.hsplist_array = lists;
    hit_list
}

fn oids(hit_list: &BlastpHitList) -> Vec<u32> {
    hit_list.hsplist_array[..hit_list.hsplist_count]
        .iter()
        .map(|list| list.oid)
        .collect()
}

#[test]
fn volumes_merge_by_the_first_alignment_of_each_subject() {
    // Final: 1 (1e-5), 2 (1e-3), 3 (1e-1). Input: 10 (1e-6), 11 (1e-3, a larger bit score
    // than 2), 12 (1e-3, the same bit score as 2 or less), 13 (5).
    let mut lists = vec![Some(hit_list(vec![
        list(1, 1e-5, 50.0),
        list(2, 1e-3, 40.0),
        list(3, 1e-1, 30.0),
    ]))];
    let input = vec![Some(hit_list(vec![
        list(10, 1e-6, 55.0),
        list(11, 1e-3, 41.0),
        list(12, 1e-3, 40.0),
        list(13, 5.0, 20.0),
    ]))];
    merge_volume(&mut lists, input);
    // 10 before 1; at equal E-values 11 (larger bit score) before 2, but 12 (equal) after 2;
    // 13 goes after 3 (the final side runs out: the rest is appended).
    assert_eq!(
        oids(lists[0].as_ref().unwrap()),
        vec![10, 1, 11, 2, 12, 3, 13]
    );
}

#[test]
fn a_volume_without_hits_changes_nothing_and_an_empty_first_volume_takes_the_next() {
    let mut lists = vec![Some(hit_list(vec![list(1, 1e-5, 50.0)])), None];
    let input = vec![
        None,
        Some(hit_list(vec![list(10, 1e-6, 55.0), list(11, 1e-2, 9.0)])),
    ];
    merge_volume(&mut lists, input);
    assert_eq!(oids(lists[0].as_ref().unwrap()), vec![1]);
    assert_eq!(oids(lists[1].as_ref().unwrap()), vec![10, 11]);
}

#[test]
fn the_input_is_cut_nowhere_and_volumes_merge_in_order() {
    let mut lists = vec![Some(hit_list(vec![list(1, 1e-2, 10.0)]))];
    let second = vec![Some(hit_list(vec![
        list(10, 1e-3, 10.0),
        list(11, 1e-1, 10.0),
    ]))];
    let third = vec![Some(hit_list(vec![list(20, 1e-2, 10.0)]))];
    merge_volume(&mut lists, second);
    merge_volume(&mut lists, third);
    // 20 ties with 1 on both values: the final side's (already merged) group goes first.
    assert_eq!(oids(lists[0].as_ref().unwrap()), vec![10, 1, 20, 11]);
}

#[test]
fn volumes_are_ordered_by_file_name_then_path() {
    use crate::seqdb::SeqDb;
    use std::path::Path;
    let mut paths = ["zd/x", "ad/y", "prof_rps2v.01", "b/x", "a/x"]
        .map(Path::new)
        .to_vec();
    paths.sort_by(|a, b| SeqDb::compare_volume_paths(a, b));
    let sorted: Vec<&str> = paths.iter().map(|p| p.to_str().unwrap()).collect();
    // prof_rps2v.01, then x (a/x, b/x, zd/x by path), then y.
    assert_eq!(sorted, ["prof_rps2v.01", "a/x", "b/x", "zd/x", "ad/y"]);
}

#[test]
fn e_values_below_the_smallest_compare_equal_and_the_bit_score_decides() {
    // 1e-189 (535 bits) and 1e-190 (534 bits) are both 0 in the Seq-aligns: the larger bit
    // score first, whichever has the smaller E-value.
    let mut lists = vec![Some(hit_list(vec![
        list(1, 1e-189, 535.0),
        list(2, 1e-100, 300.0),
    ]))];
    let input = vec![Some(hit_list(vec![
        list(10, 1e-190, 534.0),
        list(11, 1e-190, 536.0),
    ]))];
    merge_volume(&mut lists, input);
    // 10 (534) after 1 (535); 11 (536) is compared with 2 (1e-100): before it.
    assert_eq!(oids(lists[0].as_ref().unwrap()), vec![1, 10, 11, 2]);
    let mut lists = vec![Some(hit_list(vec![list(1, 1e-189, 535.0)]))];
    let input = vec![Some(hit_list(vec![list(10, 1e-190, 536.0)]))];
    merge_volume(&mut lists, input);
    assert_eq!(oids(lists[0].as_ref().unwrap()), vec![10, 1]);
}
