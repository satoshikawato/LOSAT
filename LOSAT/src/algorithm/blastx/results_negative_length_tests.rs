use super::*;

fn hsp(context: usize, score: i32, q_start: i32, length: i32) -> Hsp {
    Hsp {
        hsp: PreliminaryHsp {
            context,
            frame: 1,
            score,
            q_start,
            q_end: q_start + length,
            q_gapped_start: q_start,
            s_start: 5,
            s_end: 5 + length,
            s_gapped_start: 5,
        },
        num: 1,
        evalue: 1.0,
        bit_score: 0.0,
        identity: 0,
        positive: 0,
        edit_script: Vec::new(),
        composition_method: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2566-2569
/// An HSP of length -4 at the start of its context ends at -4: its end + 3 is below 0, so
/// NCBI's range ends at -4 itself and keeps a weaker HSP ending at -2 (an end + 3 of -1
/// would remove it), and removes one that ends at -4 or before.
#[test]
fn subject_besthit_keeps_the_end_below_minus_three() {
    let ends = |hsps: &[Hsp]| {
        hsps.iter()
            .map(|h| (h.hsp.q_start, h.hsp.q_end))
            .collect::<Vec<_>>()
    };
    let mut hsps = vec![hsp(0, 9, 0, -4), hsp(0, 5, 1, -3), hsp(0, 4, 0, -5)];
    subject_besthit(&mut hsps);
    assert_eq!(ends(&hsps), vec![(0, -4), (1, -2)]);
    // An end + 3 of 0 or more is the range's end, as before.
    let mut hsps = vec![hsp(0, 9, 0, -3), hsp(0, 5, 1, -2)];
    subject_besthit(&mut hsps);
    assert_eq!(ends(&hsps), vec![(0, -3)]);
}
