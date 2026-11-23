use super::*;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:998-1082
// When `db_length` exceeds `subject_length`, `y_variable` scales by the log of their ratio times
// `K` divided by `gap_decay_rate`.
// Expected values evaluated from the extracted NCBI C function.
#[test]
fn complete_database_changes_linking_cutoff() {
    let params = KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.7916,
        beta: -3.2,
    };
    for (q, subject, db, small, big, probability) in [
        (300, 900, 0, 30, 22, 0.5),
        (300, 900, 2700, 30, 30, 0.5),
        (30, 90, 90, 0, 8, 0.0),
        (30, 90, 2700, 0, 14, 0.0),
    ] {
        let c = calculate_link_hsp_cutoffs_ncbi(q, subject, db, 30, &params);
        assert_eq!(
            (c.cutoff_small_gap, c.cutoff_big_gap, c.gap_prob),
            (small, big, probability)
        );
    }
}
