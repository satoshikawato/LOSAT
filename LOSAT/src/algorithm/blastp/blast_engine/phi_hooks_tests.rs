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
    s_pat_len: i32,
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
        s_pattern_length: s_pat_len,
        q_pattern_offset: q_pat,
        s_pattern_offset: s_pat,
        query_start: qs,
        query_end: qe,
        subject_start: ss,
        subject_end: se,
        score,
        evalue: 0.0,
        bit_score: 0.0,
        num: 1,
        edit_script: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:780-794
// `select_kept_prelim` keeps only the preliminary HSPs surviving culling / best-hit filtering
// in `list.hsps` order, distinguishing variable-length pattern hits that share pattern starts,
// preliminary bounds, and raw score via their `PhiPrelimOrdinal`.
#[test]
fn select_kept_prelim_distinguishes_same_start_variable_length_hits_in_kept_order() {
    let a = phi_hsp(0, 3, 10, 20, 5, 40, 15, 50, 90);
    let b_short = phi_hsp(1, 3, 30, 40, 25, 60, 35, 70, 70);
    let b_long = phi_hsp(2, 4, 30, 40, 25, 60, 35, 70, 70);
    let culled = phi_hsp(3, 4, 80, 10, 75, 100, 5, 30, 50);
    let all = vec![a.clone(), b_short.clone(), b_long.clone(), culled];
    let prelims: Vec<XPrelimHsp> = all
        .iter()
        .enumerate()
        .map(|(index, hsp)| {
            preliminary_of(hsp, 0, 0, 200).with_phi_ordinal(PhiPrelimOrdinal::new(index))
        })
        .collect();
    // Simulate a culling/best-hit writer that keeps `b_long` before `b_short`, then `a`,
    // and drops `culled`.
    let kept = vec![prelims[2], prelims[1], prelims[0]];
    let selected = select_kept_prelim(&kept, all);
    assert_eq!(selected, vec![b_long, b_short, a]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3572-3632
// Occurrence indices past `u16::MAX` (e.g. 64 and 65,600, congruent mod 65,536) and
// equal-coordinate HSPs from distinct occurrences split into their exact occurrence hit lists
// via `PhiOccurrenceIndex`, resetting `subject_frame` back to 0.
#[test]
fn split_hit_list_by_occurrence_keeps_indices_past_u16_max_and_reordered_ties_separate() {
    let make_hsp = |q_start: usize, score: i32, occ: usize| {
        BlastpHsp {
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
        }
        .with_phi_occurrence(PhiOccurrenceIndex::new(occ))
    };
    let hsp_low = make_hsp(65, 200, 64);
    // Identical coordinates/score to `hsp_low`, different occurrence and mismatch count, placed
    // before `hsp_high` to prove order-independent splitting.
    let mut hsp_tie = make_hsp(65, 200, 65);
    hsp_tie.mismatch = 2;
    let hsp_high = make_hsp(65_601, 180, 65_600);
    let mut hit_list = BlastpHitList::new(1);
    hit_list.hsplist_array.push(BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![hsp_tie, hsp_low, hsp_high],
        best_evalue: 1e-10,
    });
    hit_list.hsplist_count = 1;
    let parts = split_hit_list_by_occurrence(Some(&mut hit_list), 65_689, |_oid, hsp| {
        hsp.take_phi_occurrence()
            .expect("PHI occurrence attached")
            .index()
    });
    assert_eq!(parts.len(), 65_689);
    let list_64 = parts[64].as_ref().expect("occurrence 64 has an HSP");
    assert_eq!(list_64.hsplist_count, 1);
    assert_eq!(list_64.hsplist_array[0].hsps.len(), 1);
    assert_eq!(list_64.hsplist_array[0].hsps[0].q_start, 65);
    assert_eq!(list_64.hsplist_array[0].hsps[0].mismatch, 0);
    assert_eq!(list_64.hsplist_array[0].hsps[0].subject_frame, 0);
    let list_65 = parts[65].as_ref().expect("occurrence 65 has an HSP");
    assert_eq!(list_65.hsplist_count, 1);
    assert_eq!(list_65.hsplist_array[0].hsps.len(), 1);
    assert_eq!(list_65.hsplist_array[0].hsps[0].q_start, 65);
    assert_eq!(list_65.hsplist_array[0].hsps[0].mismatch, 2);
    assert_eq!(list_65.hsplist_array[0].hsps[0].subject_frame, 0);
    let list_65600 = parts[65_600].as_ref().expect("occurrence 65600 has an HSP");
    assert_eq!(list_65600.hsplist_count, 1);
    assert_eq!(list_65600.hsplist_array[0].hsps.len(), 1);
    assert_eq!(list_65600.hsplist_array[0].hsps[0].q_start, 65_601);
    assert_eq!(list_65600.hsplist_array[0].hsps[0].subject_frame, 0);
    assert_eq!(parts.iter().filter(|p| p.is_some()).count(), 3);
    assert!(hit_list.hsplist_array[0]
        .hsps
        .iter()
        .all(|hsp| hsp.subject_frame == 0));
}
