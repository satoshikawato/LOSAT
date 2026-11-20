//! Tests of the Smith-Waterman redo of one match on a hand-made profile: a consensus with all
//! frequency ratios 0 (BLOSUM62's rows of its letters), the redo's scaled block, and RPS-BLAST's
//! roles (the profile is the redo's query, the user's query its subject).

use super::super::ProteinSubject;
use super::*;
use crate::algorithm::pssm_search::evalue::EvalueStats;
use crate::algorithm::pssm_search::gapalign::PssmGapScratch;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::compo_scale::PosMatrixInfo;
use crate::core::composition_adjustment::smith_waterman::evalue_from_score;
use crate::core::psi::KarlinBlk;
use crate::stats::KarlinParams;
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::utils::matrix::aa_char_to_ncbistdaa;

fn encode(text: &str) -> Vec<u8> {
    text.bytes().map(aa_char_to_ncbistdaa).collect()
}

/// The matrix information of a profile whose frequency ratios are all 0.
fn info_of(consensus: &[u8]) -> PosMatrixInfo {
    let sbp = crate::algorithm::pssm_search::in_pssm::score_block(
        ScoringMatrix::Blosum62,
        KarlinBlk {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
        },
    )
    .unwrap();
    let ratios = vec![[0.0; crate::core::blast_util::BLASTAA_SIZE]; consensus.len()];
    PosMatrixInfo::new(
        consensus,
        &ratios,
        &sbp,
        SCALING_FACTOR,
        KarlinBlk {
            lambda: 0.3176,
            k: 0.134,
            h: 0.4012,
        },
    )
    .unwrap()
}

/// The redo run of RPS-BLAST for `consensus`: Lambda over the scale, K, the zero PSSM.
fn run_of(consensus: &[u8]) -> RedoRun<'_> {
    let k = 1.2 * 0.041;
    RedoRun::new(
        consensus,
        consensus,
        KarlinParams {
            lambda: 0.267 / SCALING_FACTOR,
            k,
            h: 0.14,
            alpha: 1.9,
            beta: -30.0,
        },
        k.glibc_ln(),
        (11, 1),
        (25.0, 0),
        QueryPssm {
            rows: vec![[0; crate::core::blast_util::BLASTAA_SIZE]; consensus.len()],
        },
        PssmGapScratch::new(),
    )
    .with_bound_check()
}

fn redo_match(
    info: &PosMatrixInfo,
    length: usize,
    search_space: i64,
    expect: f64,
) -> RedoMatch<'_> {
    RedoMatch {
        info,
        stats: EvalueStats::new(None, search_space),
        expect_value: expect,
        query_length: length as i32,
        rps_roles: true,
        matrix: ScoringMatrix::Blosum62,
        link: None,
    }
}

/// A preliminary HSP of the profile's `q` and the subject's `s` (only the first HSP's ends and the
/// number of HSPs reach the Smith-Waterman redo).
fn hsp(q: (i32, i32), s: (i32, i32), score: i32) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: q.0,
        query_end: q.1,
        subject_start: s.0,
        subject_end: s.1,
        gapped_query_start: q.0,
        gapped_subject_start: s.0,
        raw_score: score,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
    }
}

const PROFILE: &str = "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEQVQEELLSSQ";

/// The subject: the profile's letters twice, the second copy with two changes, between letters
/// that SEG leaves (`segmasker -window 10 -locut 1.8 -hicut 2.1` masks none of the subject).
fn two_copies() -> Vec<u8> {
    let first = PROFILE;
    let second = PROFILE
        .replacen("WAAL", "WGAL", 1)
        .replacen("RQQ", "RNQ", 1);
    encode(&format!("HRDYCWTENI{first}WDKHNCPYTRFIGM{second}YCHWRDTN"))
}

// Two preliminary HSPs: the redo finds the two copies, the first copy first (the higher score);
// the window's list has the last one first, and the evaluated list is sorted by score.
#[test]
fn several_hsps_give_several_alignments() {
    let consensus = encode(PROFILE);
    let info = info_of(&consensus);
    let subject = two_copies();
    let mut run = run_of(&consensus);
    let mut sw = RpsSwScratch::new();
    let mut heap = BlastCompoHeap::new(500, 0.002);
    let n = consensus.len() as i32;
    let hsps = [
        hsp((0, n), (10, 10 + n), 300),
        hsp((0, n), (n + 24, 2 * n + 24), 250),
    ];
    let matching = redo_match(&info, consensus.len(), 20_000, 10.0);
    let alignments = matching
        .redo_one_match_smith_waterman(&mut run, &mut sw, &mut heap, 7, &hsps, &subject)
        .unwrap()
        .unwrap();
    assert_eq!(alignments.len(), 2, "{alignments:?}");
    assert!(
        alignments[0].match_start > alignments[1].match_start,
        "{alignments:?}"
    );
    assert!(alignments[1].score > alignments[0].score, "{alignments:?}");
    assert_eq!(
        (alignments[1].query_start, alignments[1].match_start),
        (0, 10)
    );
    assert_eq!(alignments[1].query_end, n);
    // The rule of a position-based adjustment.
    assert!(alignments
        .iter()
        .all(|align| align.rule == EMatrixAdjustRule::CompoScaleOldMatrix));
    let mut run = run_of(&consensus);
    let (list, best_score, best_evalue) = matching
        .redo_match_smith_waterman(&mut run, &mut sw, &mut heap, 7, &hsps, &subject)
        .unwrap()
        .unwrap();
    assert_eq!(list.oid, 7);
    assert_eq!(list.hsps.len(), 2);
    assert!(list.hsps[0].raw_score >= list.hsps[1].raw_score);
    assert_eq!(list.hsps[0].s_start, 11);
    assert_eq!(best_score, alignments[1].score);
    assert_eq!(best_evalue, list.hsps[0].e_value);
    // The heap is asked only: it stays empty.
    assert_eq!(heap.n, 0);
}

// `evaluate` is a copy of `redo_match`'s steps after the window (DECISIONS 63 of line PS): the X-drop
// redo's alignments through it give `redo_match`'s list.
#[test]
fn the_evaluation_is_the_x_drop_redos() {
    let consensus = encode(PROFILE);
    let info = info_of(&consensus);
    let subject = two_copies();
    let n = consensus.len() as i32;
    for hsps in [
        vec![
            hsp((0, n), (10, 10 + n), 300),
            hsp((0, n), (n + 24, 2 * n + 24), 250),
        ],
        vec![hsp((0, n), (n + 24, 2 * n + 24), 250)],
        vec![hsp((0, 30), (10, 40), 90), hsp((20, n), (30, 10 + n), 120)],
    ] {
        for (search_space, expect) in [(20_000, 10.0), (20_000, 1e-60), (1_000_000_000_000, 10.0)] {
            let matching = redo_match(&info, consensus.len(), search_space, expect);
            let mut run = run_of(&consensus);
            let alignments = matching.redo_one_match(&mut run, &hsps, &subject);
            let copied = matching.evaluate(&run, 7, &hsps, &subject, alignments);
            let original = matching.redo_match(&mut run_of(&consensus), 7, &hsps, &subject);
            assert_eq!(format!("{copied:?}"), format!("{original:?}"));
        }
    }
}

/// `two_copies` with a stretch of one letter between the copies, which SEG masks (the redo's
/// `subject_maybe_biased`).
fn low_complexity_between() -> Vec<u8> {
    let second = PROFILE
        .replacen("WAAL", "WGAL", 1)
        .replacen("RQQ", "RNQ", 1);
    encode(&format!(
        "HRDYCWTENI{PROFILE}QQQQQQQQQQQQQQ{second}YCHWRDTN"
    ))
}

// `redo_windows` (rpstblastn's redo) is a copy of `redo_match`'s steps (DECISIONS 63 S3): its
// instance of a protein subject, one window over the whole subject, gives `redo_match`'s list.
// The subjects reach the copy's branches: SEG masks a stretch of `low_complexity_between`
// (`subject_maybe_biased`), a subject of X alone has no true amino acid (`adjust_search_failed`),
// and the redo links its alignments where psiblast's `-sum_stats true` does (`link`).
#[test]
fn the_window_copy_is_the_x_drop_redo() {
    let consensus = encode(PROFILE);
    let info = info_of(&consensus);
    let n = consensus.len() as i32;
    let biased = low_complexity_between();
    assert!(!crate::utils::seg::SegMasker::new(10, 1.8, 2.1)
        .mask_sequence(&biased)
        .is_empty());
    let all_x = vec![aa_char_to_ncbistdaa(b'X'); biased.len()];
    let mut found = [0usize; 3];
    for (index, subject) in [two_copies(), biased, all_x].iter().enumerate() {
        for hsps in [
            vec![
                hsp((0, n), (10, 10 + n), 300),
                hsp((0, n), (n + 24, 2 * n + 24), 250),
            ],
            vec![hsp((0, n), (n + 24, 2 * n + 24), 250)],
            vec![hsp((0, 30), (10, 40), 90), hsp((20, n), (30, 10 + n), 120)],
            vec![hsp((0, 30), (10, 40), 90), hsp((0, 30), (10, 40), 90)],
            vec![
                hsp((0, 30), (10, 40), 90),
                hsp((0, 30), (n + 24, n + 54), 80),
            ],
        ] {
            for (search_space, expect) in
                [(20_000, 10.0), (20_000, 1e-60), (1_000_000_000_000, 10.0)]
            {
                for cutoff in [None, Some(1), Some(100)] {
                    let mut matching = redo_match(&info, consensus.len(), search_space, expect);
                    matching.link = cutoff.map(|cutoff| {
                        super::super::RedoLink::new(
                            cutoff,
                            crate::core::link_hsps::LinkContext {
                                lambda: 0.267,
                                log_k: 0.041f64.glibc_ln(),
                                query_length: n,
                                length_adjustment: 0,
                                eff_searchsp: search_space,
                            },
                        )
                    });
                    let copied = matching.redo_windows(
                        &mut run_of(&consensus),
                        7,
                        &hsps,
                        subject,
                        ProteinSubject,
                    );
                    let original = matching.redo_match(&mut run_of(&consensus), 7, &hsps, subject);
                    assert_eq!(format!("{copied:?}"), format!("{original:?}"));
                    found[index] += usize::from(original.is_some());
                }
            }
        }
    }
    // Alignments on both subjects of letters; none where no letter is a true amino acid.
    assert!(found[0] > 0 && found[1] > 0, "{found:?}");
    assert_eq!(found[2], 0);
}

// One preliminary HSP: the redo stops after its first alignment (`window->hspcnt > 1`).
#[test]
fn one_hsp_gives_one_alignment() {
    let consensus = encode(PROFILE);
    let info = info_of(&consensus);
    let subject = two_copies();
    let mut run = run_of(&consensus);
    let n = consensus.len() as i32;
    let alignments = redo_match(&info, consensus.len(), 20_000, 10.0)
        .redo_one_match_smith_waterman(
            &mut run,
            &mut RpsSwScratch::new(),
            &mut BlastCompoHeap::new(500, 0.002),
            7,
            &[hsp((0, n), (10, 10 + n), 300)],
            &subject,
        )
        .unwrap()
        .unwrap();
    assert_eq!(alignments.len(), 1);
    assert_eq!(alignments[0].match_start, 10, "{alignments:?}");
}

// The first alignment's significance: its e-value below `-evalue`, and the heap of redone
// matches would insert it; otherwise the match gets no alignment at all.
#[test]
fn the_first_alignment_must_be_significant_and_wanted_by_the_heap() {
    let consensus = encode(PROFILE);
    let info = info_of(&consensus);
    let subject = two_copies();
    let n = consensus.len() as i32;
    let hsps = [
        hsp((0, n), (10, 10 + n), 300),
        hsp((0, n), (n + 24, 2 * n + 24), 250),
    ];
    let alignments = |expect: f64, heap: &mut BlastCompoHeap| {
        let mut run = run_of(&consensus);
        redo_match(&info, consensus.len(), 20_000, expect)
            .redo_one_match_smith_waterman(
                &mut run,
                &mut RpsSwScratch::new(),
                heap,
                7,
                &hsps,
                &subject,
            )
            .unwrap()
            .unwrap()
            .len()
    };
    assert_eq!(alignments(10.0, &mut BlastCompoHeap::new(500, 0.002)), 2);
    // An e-value threshold below every alignment's.
    assert_eq!(alignments(1e-300, &mut BlastCompoHeap::new(500, 0.002)), 0);
    // A full heap whose worst match is better: `BlastCompo_HeapWouldInsert` says no.
    let mut full = BlastCompoHeap::new(1, 0.0);
    full.insert(
        BlastpHspList {
            oid: 1,
            query_index: 0,
            hsps: Vec::new(),
            best_evalue: 1e-300,
        },
        1e-300,
        i32::MAX,
        1,
    );
    assert_eq!(alignments(10.0, &mut full), 0);
}

// An alignment of score 0 that passes the threshold (a search space of 1): with nothing positive
// at (0, 0) the X-drop alignment has no operations and NCBI's redo of the match fails (no HSP);
// after a range forbids the positive cell at (0, 0), NCBI finds that alignment again and again.
#[test]
fn alignments_of_score_zero_fail_or_never_end() {
    let consensus = encode("W");
    let info = info_of(&consensus);
    let matching = redo_match(&info, 1, 1, 10.0);
    let subject = encode("A");
    let outcome = matching.redo_one_match_smith_waterman(
        &mut run_of(&consensus),
        &mut RpsSwScratch::new(),
        &mut BlastCompoHeap::new(500, 0.002),
        0,
        &[hsp((0, 1), (0, 1), 1)],
        &subject,
    );
    assert_eq!(outcome.unwrap().map(|list| list.len()), None);
    let subject = encode("W");
    let hsps = [hsp((0, 1), (0, 1), 11), hsp((0, 1), (0, 1), 11)];
    let outcome = matching.redo_one_match_smith_waterman(
        &mut run_of(&consensus),
        &mut RpsSwScratch::new(),
        &mut BlastCompoHeap::new(500, 0.002),
        0,
        &hsps,
        &subject,
    );
    assert_eq!(outcome.unwrap_err(), SwEndless);
    // With one HSP the loop stops after the first alignment.
    let outcome = matching.redo_one_match_smith_waterman(
        &mut run_of(&consensus),
        &mut RpsSwScratch::new(),
        &mut BlastCompoHeap::new(500, 0.002),
        0,
        &hsps[..1],
        &subject,
    );
    assert_eq!(outcome.unwrap().map(|list| list.len()), Some(1));
}

// The regime of a significant score 0 with two or more HSPs is decided at the first alignment by
// profile row 0 against the subject's first letter: positive, NCBI never ends (the refusal);
// otherwise NCBI makes every positive alignment, then fails and drops the match. With a score of 0
// not significant (a large search space) the same windows run the loop as usual.
#[test]
fn a_significant_score_zero_is_decided_by_the_first_cell() {
    let two = [hsp((0, 2), (0, 2), 11), hsp((0, 2), (0, 2), 11)];
    let outcome = |profile: &str, subject: &str, search_space: i64| {
        let consensus = encode(profile);
        let info = info_of(&consensus);
        redo_match(&info, consensus.len(), search_space, 10.0)
            .redo_one_match_smith_waterman(
                &mut run_of(&consensus),
                &mut RpsSwScratch::new(),
                &mut BlastCompoHeap::new(500, 0.002),
                0,
                &two,
                &encode(subject),
            )
            .map(|list| list.map(|list| list.len()))
    };
    // W against W first: never ends in NCBI.
    assert_eq!(outcome("WCWC", "WCWC", 1), Err(SwEndless));
    // A against K first (negative), W against W after: every positive alignment, then the
    // failed redo.
    assert_eq!(outcome("AWCH", "KWCH", 1), Ok(None));
    // The score of 0 not significant: the positive alignments.
    assert!(matches!(outcome("AWCH", "KWCH", 1_000), Ok(Some(n)) if n >= 1));
}

// `s_EvalueFromScore` reads the search space whatever the Gumbel block.
#[test]
fn the_evalue_of_a_score_is_karlin_altschuls() {
    let lambda = 0.267 / SCALING_FACTOR;
    let log_k = (1.2f64 * 0.041).glibc_ln();
    assert_eq!(
        evalue_from_score(3200, lambda, log_k, 1e6),
        1e6 * (-(lambda * 3200.0) + log_k).glibc_exp()
    );
    assert_eq!(evalue_from_score(0, lambda, log_k, 1.0), log_k.glibc_exp());
}
