use super::*;

fn hit(raw_score: i32, marker: usize) -> UngappedHit {
    UngappedHit {
        q_idx: marker as u32,
        s_idx: (marker + 10) as u32,
        ctx_idx: u32::try_from(marker + 20).unwrap(),
        s_f_idx: u8::try_from(marker + 30).unwrap(),
        q_frame: if marker.is_multiple_of(2) { -1 } else { 1 },
        s_frame: if marker.is_multiple_of(2) { -2 } else { 2 },
        q_aa_start: u32::try_from(marker + 40).unwrap(),
        q_aa_end: u32::try_from(marker + 50).unwrap(),
        s_aa_start: u32::try_from(marker + 60).unwrap(),
        s_aa_end: u32::try_from(marker + 70).unwrap(),
        q_seed_off: u32::try_from(marker + 80).unwrap(),
        s_seed_off: u32::try_from(marker + 90).unwrap(),
        q_orig_len: u32::try_from(marker + 100).unwrap(),
        s_orig_len: u32::try_from(marker + 110).unwrap(),
        raw_score,
        e_value: marker as f64 / 10.0,
        num_ident: u32::try_from(marker + 120).unwrap(),
        hsp_list_order: u32::try_from(marker + 130).unwrap(),
        ordering_method: (marker % 2) as u8,
        linked_set: marker.is_multiple_of(2),
        start_of_chain: !marker.is_multiple_of(2),
        link_id: u32::try_from(marker + 140).unwrap(),
        chain_next_link_id: Some(u32::try_from(marker + 150).unwrap()),
        hsp_link_num: 0,
        num: 0,
    }
}

#[test]
fn comparator_equal_hits_preserve_incoming_order_without_an_extra_key() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1355
    let mut hits = vec![hit(20, 0), hit(30, 2), hit(30, 1), hit(10, 3)];

    sort_ungapped_hits_by_index_replay(&mut hits, |a, b| b.raw_score.cmp(&a.raw_score));

    let order: Vec<usize> = hits
        .iter()
        .map(|value| value.hsp_list_order as usize)
        .collect();
    assert_eq!(order, vec![132, 131, 130, 133]);
}

#[test]
fn index_replay_preserves_each_complete_original_record() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-486
    let mut hits = vec![hit(20, 1), hit(30, 2), hit(10, 3)];
    let expected = [1, 0, 2]
        .into_iter()
        .map(|index| format!("{:?}", hits[index]))
        .collect::<Vec<_>>();

    sort_ungapped_hits_by_index_replay(&mut hits, |a, b| b.raw_score.cmp(&a.raw_score));

    assert_eq!(
        hits.iter()
            .map(|value| format!("{value:?}"))
            .collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn index_replay_preserves_link_relation_fields() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:960-974
    let mut hits = vec![hit(10, 1), hit(20, 2)];

    sort_ungapped_hits_by_index_replay(&mut hits, |a, b| b.raw_score.cmp(&a.raw_score));

    let relations = hits
        .iter()
        .map(|value| {
            (
                value.link_id,
                value.chain_next_link_id,
                value.linked_set,
                value.start_of_chain,
                value.ordering_method,
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        relations,
        vec![
            (142, Some(152), true, false, 0),
            (141, Some(151), false, true, 1)
        ]
    );
}
