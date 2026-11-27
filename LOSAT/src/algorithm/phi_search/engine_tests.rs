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
    assert_eq!(kept[0].occurrence, 0);
}

#[test]
fn sum_statistics_link_the_preliminary_hsps_and_reap_them_by_the_preliminary_evalue() {
    let pattern = PhiPattern::compile(b"G-K-S").unwrap();
    let score = |q: u8, s: u8| blosum62(q, s);
    let query = encode("MKAGKSTLLW");
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
    let mut search = setup(&setup_params).unwrap();
    let at = |q: (i32, i32), s: (i32, i32), score: i32| PhiHsp {
        occurrence: 0,
        s_pattern_length: 3,
        q_pattern_offset: q.0,
        s_pattern_offset: s.0,
        query_start: q.0,
        query_end: q.1,
        subject_start: s.0,
        subject_end: s.1,
        score,
        evalue: 0.0,
        bit_score: 0.0,
        num: 1,
        edit_script: Vec::new(),
    };
    // In order on both sequences (one set of two), the lower score first.
    let hsps = vec![at((200, 260), (500, 560), 60), at((10, 60), (20, 70), 80)];
    // Without sum statistics: nothing changes.
    let mut unlinked = hsps.clone();
    search.link_preliminary(&mut unlinked, 1000);
    assert_eq!(unlinked, hsps);
    let link = |prelim_evalue: f64| PhiLink {
        params: crate::core::link_hsps::EvenGapLinkParams::preliminary(stats::word_cutoff(
            22.0,
            &search.karlin,
            search.cutoff,
        )),
        context: stats::link_context(&search.karlin, 300, 3, 1_000_000),
        prelim_evalue,
    };
    search.link = Some(link(f64::MAX));
    let mut linked = hsps.clone();
    search.link_preliminary(&mut linked, 1000);
    assert_eq!(
        linked
            .iter()
            .map(|hsp| (hsp.score, hsp.num))
            .collect::<Vec<_>>(),
        vec![(80, 2), (60, 2)]
    );
    let set = linked[0].evalue;
    assert!(set > 0.0 && linked[1].evalue == set, "{linked:?}");
    // The reap frees the HSPs past the preliminary e-value and keeps those at it.
    search.link = Some(link(set));
    let mut kept = hsps.clone();
    search.link_preliminary(&mut kept, 1000);
    assert_eq!(kept.len(), 2);
    search.link = Some(link(set / 2.0));
    let mut reaped = hsps;
    search.link_preliminary(&mut reaped, 1000);
    assert!(reaped.is_empty());
}
