use super::*;

fn id(text: &str) -> SetKey {
    use crate::seqid::accguide::parse_flags::{ANY_LOCAL, PARTIAL_OK, RAW_TEXT};
    let ids = crate::seqid::parse::parse_ids(
        text.as_bytes(),
        PARTIAL_OK | ANY_LOCAL | RAW_TEXT,
        &mut |_| {},
    )
    .unwrap();
    SetKey::Id(ids.into_iter().next().unwrap().unwrap())
}

// Review 1 of line DBID pssmrows, M1: only a shown subject that a previous key can be one
// handle with is refused (blast_format.cpp:541-553, cited at `refuse_unsure_shown`).
#[test]
fn only_a_shown_unsure_subject_is_refused() {
    let found = FoundBefore {
        found: vec![true, false, false],
        unsure: [(2, (id("gb|ABC12.1|"), id("emb|ABC12.1|")))]
            .into_iter()
            .collect(),
    };
    found.refuse_unsure_shown([0, 1]).unwrap();
    let error = found.refuse_unsure_shown([0, 2]).unwrap_err().to_string();
    assert!(error.contains("gb|ABC12.1| and emb|ABC12.1|"), "{error}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:512-545 (cited at
// `seqalign_pssm_evalue`).
#[test]
fn the_pssm_reads_an_absent_evalue_as_10_and_any_negative_one_as_0() {
    use crate::report::serial::seqalign::ABSENT_EVALUE;
    // A NaN's, as the lists of a PSSM's search carry it, and a NaN itself: none in the Seq-align.
    assert_eq!(seqalign_pssm_evalue(ABSENT_EVALUE), 10.0);
    assert_eq!(seqalign_pssm_evalue(f64::NAN), 10.0);
    // A database of a negative size gives negative E-values, which are 0 in the Seq-align.
    assert_eq!(seqalign_pssm_evalue(-3.5), 0.0);
    assert_eq!(seqalign_pssm_evalue(-1.0e-200), 0.0);
    assert_eq!(seqalign_pssm_evalue(1.0e-200), 0.0);
    assert_eq!(seqalign_pssm_evalue(0.5), 0.5);
    assert_eq!(seqalign_pssm_evalue(1.0e9), 10.0);
    assert!(included(-3.5, 0.002));
    assert!(!included(ABSENT_EVALUE, 0.002));
}
