use super::*;
use crate::blastinput::seq_range::RangedBatch;

fn frames(lengths: &[usize]) -> Vec<Vec<QueryFrame>> {
    lengths
        .iter()
        .map(|&aa_len| {
            vec![QueryFrame {
                frame: 1,
                aa_seq: vec![0; aa_len + 2],
                aa_seq_nomask: None,
                aa_len,
                orig_len: aa_len,
                seg_masks: Vec::new(),
            }]
        })
        .collect()
}

// psiblast `-phi_pattern` with one round (PS8a fix round 1): a batch of queries without letters
// and the next query is searched as that query, under the first query's Seq-id (NCBI BLAST+
// 2.17.0 `PS/tmp/rev8a/ncbi3/er.first.0`, `er.mid.0`); a batch's report shows the searched
// query's rows with gaps only where its first query has no letters.
#[test]
fn a_phi_batch_is_searched_as_its_query_with_letters() {
    let batches = [
        RangedBatch {
            input: 0..1,
            searched: 0..1,
        },
        RangedBatch {
            input: 1..3,
            searched: 1..3,
        },
        RangedBatch {
            input: 3..3,
            searched: 3..3,
        },
    ];
    let phi = phi_batches(&batches, &frames(&[7, 0, 5]));
    assert_eq!(
        phi,
        [
            PhiBatch {
                queries: 0..1,
                searched: 0,
            },
            PhiBatch {
                queries: 1..3,
                searched: 2,
            },
        ]
    );
    assert_eq!((phi[0].first(), phi[0].shows_gaps()), (0, false));
    assert_eq!((phi[1].first(), phi[1].shows_gaps()), (1, true));
}
