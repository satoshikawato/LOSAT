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

fn phi_hsp(
    occ: usize,
    q_pat: i32,
    s_pat: i32,
    qs: i32,
    qe: i32,
    ss: i32,
    se: i32,
    score: i32,
) -> PhiHsp {
    PhiHsp {
        occurrence: occ,
        s_pattern_length: 6,
        q_pattern_offset: q_pat,
        s_pattern_offset: s_pat,
        query_start: qs,
        query_end: qe,
        subject_start: ss,
        subject_end: se,
        score,
        evalue: 0.0,
        bit_score: 0.0,
        edit_script: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:780-794
// `select_kept_prelim` keeps only the preliminary HSPs surviving culling / best-hit filtering,
// preserving multiset counts and the order of `list.hsps`.
#[test]
fn select_kept_prelim_filters_culled_hsps_and_preserves_order_and_multiplicity() {
    let a = phi_hsp(0, 10, 20, 5, 40, 15, 50, 90);
    let b_first = phi_hsp(1, 30, 40, 25, 60, 35, 70, 70);
    let b_second = phi_hsp(2, 30, 40, 25, 60, 35, 70, 70);
    let culled = phi_hsp(3, 80, 10, 75, 100, 5, 30, 50);
    let kept = vec![
        preliminary_of(&b_first, 0, 0, 200),
        preliminary_of(&a, 0, 0, 200),
    ];
    let selected = select_kept_prelim(&kept, vec![a.clone(), b_first.clone(), b_second, culled]);
    assert_eq!(selected, vec![b_first, a]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3572-3632
// Occurrence indices past `u16::MAX` (e.g. 64 and 65,600, which are congruent mod 65,536) must
// split into distinct occurrence hit lists without truncation.
#[test]
fn split_hit_list_by_occurrence_keeps_indices_past_u16_max_separate() {
    let make_hsp = |q_start: usize, score: i32| BlastpHsp {
        identity: 100.0,
        length: 40,
        mismatch: 0,
        gapopen: 0,
        q_start,
        q_end: q_start + 39,
        s_start: 1,
        s_end: 40,
        gapped_q_start: (q_start - 1) as i32,
        gapped_s_start: 0,
        e_value: 1e-10,
        bit_score: 80.0,
        num_ident: 40,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 65_720,
        q_idx: 0,
        s_idx: 0,
        raw_score: score,
        gap_info: None,
        num_positives: 40,
        comp_adjustment_method: 0,
        num: 1,
    };
    let hsp_low = make_hsp(65, 200);
    let hsp_high = make_hsp(65_601, 180);
    let mut hit_list = BlastpHitList::new(1);
    hit_list.hsplist_array.push(BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![hsp_low, hsp_high],
        best_evalue: 1e-10,
    });
    hit_list.hsplist_count = 1;
    let parts = split_hit_list_by_occurrence(Some(&hit_list), 65_689, |_oid, hsp| {
        if hsp.q_start == 65 {
            64
        } else {
            65_600
        }
    });
    assert_eq!(parts.len(), 65_689);
    let list_64 = parts[64].as_ref().expect("occurrence 64 has an HSP");
    assert_eq!(list_64.hsplist_count, 1);
    assert_eq!(list_64.hsplist_array[0].hsps.len(), 1);
    assert_eq!(list_64.hsplist_array[0].hsps[0].q_start, 65);
    let list_65600 = parts[65_600].as_ref().expect("occurrence 65600 has an HSP");
    assert_eq!(list_65600.hsplist_count, 1);
    assert_eq!(list_65600.hsplist_array[0].hsps.len(), 1);
    assert_eq!(list_65600.hsplist_array[0].hsps[0].q_start, 65_601);
    assert_eq!(parts.iter().filter(|p| p.is_some()).count(), 2);
}
