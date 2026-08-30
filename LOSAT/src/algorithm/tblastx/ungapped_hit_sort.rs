//! Stable index-replay sorting for TBLASTX `UngappedHit` records.

use super::chaining::UngappedHit;
use std::cmp::Ordering;

pub(crate) type UngappedHitCompare = fn(&UngappedHit, &UngappedHit) -> Ordering;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1381
// ```c
// void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
// {
//     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:476-486
// ```c
// for (index = 0; index < total_number_of_hsps; ++index) {
//    link_hsp_array[index]->hsp = hsp_array[index];
// }
// qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//       s_RevCompareHSPsTbx);
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:990-994
// ```c
// qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//       s_RevCompareHSPsTransl);
// qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
//       s_FwdCompareHSPsTransl);
// ```
// LOSAT owns complete UngappedHit values rather than NCBI pointer arrays.
// Stable-sort their original indices with the unchanged NCBI comparator, then
// replay complete original records. Comparator-equal records retain incoming
// order without adding a LOSAT-only tie-break field.
pub(crate) fn sort_ungapped_hits_by_index_replay(
    hits: &mut [UngappedHit],
    compare: UngappedHitCompare,
) {
    if hits.len() <= 1 {
        return;
    }

    let original = hits.to_vec();
    let mut indices: Vec<usize> = (0..original.len()).collect();
    indices.sort_by(|&lhs_index, &rhs_index| compare(&original[lhs_index], &original[rhs_index]));

    for (dst, src) in indices.into_iter().enumerate() {
        hits[dst] = original[src].clone();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(raw_score: i32, marker: usize) -> UngappedHit {
        UngappedHit {
            q_idx: marker as u32,
            s_idx: (marker + 10) as u32,
            ctx_idx: marker + 20,
            s_f_idx: marker + 30,
            q_frame: if marker.is_multiple_of(2) { -1 } else { 1 },
            s_frame: if marker.is_multiple_of(2) { -2 } else { 2 },
            q_aa_start: marker + 40,
            q_aa_end: marker + 50,
            s_aa_start: marker + 60,
            s_aa_end: marker + 70,
            q_seed_off: marker + 80,
            s_seed_off: marker + 90,
            q_orig_len: marker + 100,
            s_orig_len: marker + 110,
            raw_score,
            e_value: marker as f64 / 10.0,
            num_ident: marker + 120,
            hsp_list_order: marker + 130,
            ordering_method: (marker % 2) as u8,
            linked_set: marker.is_multiple_of(2),
            start_of_chain: !marker.is_multiple_of(2),
            link_id: marker + 140,
            chain_next_link_id: Some(marker + 150),
        }
    }

    #[test]
    fn comparator_equal_hits_preserve_incoming_order_without_an_extra_key() {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1347-1355
        // ```c
        // if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
        //     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
        //     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
        //     0 == (result = BLAST_CMP(hsp1->query.offset,   hsp2->query.offset))) {
        //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
        // }
        // ```
        let mut hits = vec![hit(20, 0), hit(30, 2), hit(30, 1), hit(10, 3)];

        sort_ungapped_hits_by_index_replay(&mut hits, |a, b| b.raw_score.cmp(&a.raw_score));

        let order: Vec<usize> = hits.iter().map(|value| value.hsp_list_order).collect();
        assert_eq!(order, vec![132, 131, 130, 133]);
    }

    #[test]
    fn index_replay_preserves_each_complete_original_record() {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:476-486
        // ```c
        // for (index = 0; index < total_number_of_hsps; ++index) {
        //    link_hsp_array[index]->hsp = hsp_array[index];
        // }
        // qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
        //       s_RevCompareHSPsTbx);
        // ```
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
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:960-974
        // ```c
        // H->linked_set = linked_set;
        // H->ordering_method = ordering_method;
        // H->hsp->evalue = prob[ordering_method];
        // ```
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
}
