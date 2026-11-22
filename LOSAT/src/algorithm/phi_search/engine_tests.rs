use super::super::pattern::AMINOACID_TO_NCBISTDAA;
use super::*;
use crate::config::compat::ScoringMatrix;

fn encode(text: &str) -> Vec<u8> {
    text.bytes()
        .map(|byte| AMINOACID_TO_NCBISTDAA[byte as usize])
        .collect()
}

fn blosum62(q: u8, s: u8) -> i32 {
    crate::utils::matrix::protein_score_rows(ScoringMatrix::Blosum62).score(q, s)
}

#[test]
fn x_dropoffs_use_the_table_lambda() {
    // 15 and 25 bits with BLOSUM62 11/1's 0.270.
    assert_eq!(x_dropoffs(0.270, 15.0, 25.0), (38, 64));
}

#[test]
fn setup_errors_and_a_self_search() {
    let pattern = PhiPattern::compile(b"G-K-S").unwrap();
    let score = |q: u8, s: u8| blosum62(q, s);
    let query = encode("MKAGKSTLLW");
    let base = |query: &[u8], segments: &[(usize, usize)]| -> Result<i32, PhiSetupError> {
        let setup_params = PhiSetup {
            pattern: &pattern,
            query,
            segments,
            score: &score,
            karlin: stats::phi_karlin("BLOSUM62", 11, 1).unwrap(),
            db_length: 10,
            db_num_seqs: 1,
            evalue: 10.0,
            gap_open: 11,
            gap_extend: 1,
            x_drop_bits: 15.0,
            x_drop_final_bits: 25.0,
        };
        setup(&setup_params).map(|search| search.cutoff)
    };
    assert_eq!(
        base(&encode("MKAW"), &[(0, 3)]),
        Err(PhiSetupError::NotFound(b"G-K-S".to_vec()))
    );
    assert_eq!(
        PhiSetupError::WholeQuery(b"G-K-S".to_vec()).message(),
        b"The pattern (G-K-S) may not cover the entire query.".to_vec()
    );
    assert_eq!(
        base(&encode("GKS"), &[(0, 2)]),
        Err(PhiSetupError::WholeQuery(b"G-K-S".to_vec()))
    );
    let setup_params = PhiSetup {
        pattern: &pattern,
        query: &query,
        segments: &[(0, 9)],
        score: &score,
        karlin: stats::phi_karlin("BLOSUM62", 11, 1).unwrap(),
        db_length: 10,
        db_num_seqs: 1,
        evalue: 10.0,
        gap_open: 11,
        gap_extend: 1,
        x_drop_bits: 15.0,
        x_drop_final_bits: 25.0,
    };
    let search = setup(&setup_params).unwrap();
    assert_eq!(search.search_space, 7);
    assert_eq!(
        search.occurrences,
        vec![Occurrence {
            offset: 3,
            length: 3
        }]
    );
    let mut subject = query.clone();
    subject.push(0);
    let mut scratch = PssmGapScratch::new();
    let (prelim, hits) = search.prelim(&subject, 10, &mut scratch).unwrap();
    assert_eq!(hits, 1);
    assert_eq!(prelim.len(), 1);
    assert_eq!(prelim[0].score, 38);
    let kept = search
        .traceback(&subject, 10, prelim, 1, 10.0, &mut scratch)
        .unwrap();
    assert_eq!(kept.len(), 1);
    assert_eq!((kept[0].query_start, kept[0].query_end), (0, 10));
    assert_eq!(kept[0].edit_script, vec![GapEditOp::Sub(10)]);
    let split = split_by_occurrence(&[(4, kept)], 2);
    assert_eq!(split[0].len(), 1);
    assert!(split[1].is_empty());
}
