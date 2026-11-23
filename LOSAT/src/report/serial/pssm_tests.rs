use super::*;

fn object(final_data: Option<PssmFinalData>) -> PssmWithParameters {
    PssmWithParameters {
        pssm: Pssm {
            is_protein: None,
            identifier: None,
            num_rows: 28,
            num_columns: 1,
            row_labels: None,
            by_row: None,
            query: None,
            intermediate_data: None,
            final_data,
        },
        params: None,
    }
}

/// psiblast's archive of a round holds the search's PSSM with the round's Karlin-Altschul
/// values written back (psiblast_impl.cpp:150-168): the six values of `finalData` change, its
/// scores and other fields stay, and an object without `finalData` stays without.
#[test]
fn the_final_karlin_values_are_written_back_into_final_data() {
    let data = PssmFinalData {
        scores: vec![1; 28],
        lambda: 0.267,
        kappa: 0.041,
        h: 0.14,
        scaling_factor: Some(100),
        lambda_ungapped: None,
        kappa_ungapped: None,
        h_ungapped: None,
        word_score_threshold: Some(11.0),
    };
    let karlin = FinalKarlin {
        lambda: 0.25,
        kappa: 0.05,
        h: 0.2,
        lambda_ungapped: 0.3,
        kappa_ungapped: 0.13,
        h_ungapped: 0.4,
    };
    let written = object(Some(data.clone())).with_final_karlin(karlin);
    assert_eq!(
        written.pssm.final_data,
        Some(PssmFinalData {
            lambda: 0.25,
            kappa: 0.05,
            h: 0.2,
            lambda_ungapped: Some(0.3),
            kappa_ungapped: Some(0.13),
            h_ungapped: Some(0.4),
            ..data
        })
    );
    assert_eq!(object(None).with_final_karlin(karlin), object(None));
}
