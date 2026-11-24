use super::*;

fn test_redo_align_params() -> BlastRedoAlignParams {
    BlastRedoAlignParams {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2216-2231
        matrix_info:
            crate::core::composition_adjustment::adjust_scores::build_matrix_info(
                ScoringMatrix::Blosum62,
                0.3176,
            )
            .unwrap(),
        gapping_params:
            crate::core::composition_adjustment::redo_alignment::BlastCompoGappingParams {
                gap_open: 11,
                gap_extend: 1,
                x_dropoff: 50,
                context: std::cell::Cell::new(None),
            },
        compo_adjust_mode:
            crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode::CompositionMatrixAdjust,
        alphsize: crate::utils::matrix::BLASTAA_SIZE as i32,
        composition_test_index: 0,
        log_k: 0.0,
        score_divisor: 1.0,
        smith_waterman: false,
        near_identical_cutoff: 0.0,
        position_based: false,
        ccat_query_length: 80,
        query_is_translated: false,
        subject_is_translated: false,
        cutoff_score: 0,
        do_link_hsps: false,
    }
}

/// The list `reap_contained_hits` works on: each hit through `BlastpHsp::from_hit`, sorted by
/// score as `hsplist_from_distinct_alignments` does.
fn hsp_list_of(redone_hits: &[RedoneBlastpHit]) -> BlastpHspList {
    let mut hsp_list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: redone_hits
            .iter()
            .map(|redone_hit| BlastpHsp::from_hit(redone_hit.hit.clone()))
            .collect(),
        best_evalue: i32::MAX as f64,
    };
    sort_hsplist_by_score(&mut hsp_list);
    hsp_list
}

fn make_hit(score: i32, q_start: usize, q_end: usize, s_start: usize, s_end: usize) -> Hit {
    Hit {
        identity: 0.0,
        length: q_end.saturating_sub(q_start).saturating_add(1),
        mismatch: 0,
        gapopen: 0,
        q_start,
        q_end,
        s_start,
        s_end,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_frame: 0,
        query_length: 0,
        q_idx: 0,
        s_idx: 0,
        raw_score: score,
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-143
        // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: None,
        num_positives: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:224-260
#[test]
fn test_reap_contained_hits_uses_internal_offsets_ncbi_style() {
    let outer = RedoneBlastpHit {
        hit: make_hit(200, 101, 150, 201, 250),
    };
    let inner = RedoneBlastpHit {
        hit: make_hit(180, 101, 149, 201, 249),
    };

    let kept_hits = vec![outer.clone(), inner];
    let mut hsp_list = hsp_list_of(&kept_hits);
    reap_contained_hits(&mut hsp_list);
    assert_eq!(hsp_list.hsps.len(), 1);
    assert_eq!(hsp_list.hsps[0].raw_score, outer.hit.raw_score);
    assert_eq!(hsp_list.hsps[0].q_start, outer.hit.q_start);
    assert_eq!(hsp_list.hsps[0].s_start, outer.hit.s_start);
}

#[test]
fn test_hsplist_from_distinct_alignments_uses_edit_script_context() {
    let template_hit = BlastpPreliminaryHsp {
        query_start: 0,
        query_end: 20,
        subject_start: 4,
        subject_end: 24,
        gapped_query_start: 0,
        gapped_subject_start: 4,
        raw_score: 80,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx: 9,
    };
    let mut alignments = Some(blast_compo_alignment_new(
        80,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        20,
        0,
        4,
        24,
        0,
        Some(BlastCompoAlignmentContext::EditScript(vec![
            crate::common::GapEditOp::Sub(10),
            crate::common::GapEditOp::Del(1),
            crate::common::GapEditOp::Sub(9),
        ])),
    ));
    alignments.as_mut().unwrap().next = Some(blast_compo_alignment_new(
        70,
        EMatrixAdjustRule::DontAdjustMatrix,
        1,
        18,
        0,
        7,
        24,
        0,
        Some(BlastCompoAlignmentContext::EditScript(vec![
            crate::common::GapEditOp::Sub(17),
        ])),
    ));

    let hits = hsplist_from_distinct_alignments(alignments, &template_hit, 0, 0).hsps;
    assert_eq!(hits.len(), 2);
    assert_eq!(hits[0].raw_score, 80);
    assert_eq!(hits[1].raw_score, 70);
    assert_eq!(
        hits[0].gap_info,
        Some(vec![
            crate::common::GapEditOp::Sub(10),
            crate::common::GapEditOp::Del(1),
            crate::common::GapEditOp::Sub(9),
        ])
    );
    assert_eq!(hits[0].q_start, 1);
    assert_eq!(hits[0].s_start, 5);
    assert_eq!(hits[0].s_idx, 9);
}

// The list is built with room for its alignments only.
#[test]
fn x_redo_fit_builds_the_list_with_room_for_its_alignments_only() {
    let template_hit = BlastpPreliminaryHsp {
        query_start: 0,
        query_end: 20,
        subject_start: 4,
        subject_end: 24,
        gapped_query_start: 0,
        gapped_subject_start: 4,
        raw_score: 80,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx: 9,
    };
    let alignments = |count: usize| {
        let mut head: Option<Box<BlastCompoAlignment>> = None;
        for index in 0..count {
            let mut align = blast_compo_alignment_new(
                60 + 10 * ((index as i32 * 7) % 5),
                EMatrixAdjustRule::DontAdjustMatrix,
                index as i32,
                20 + index as i32,
                0,
                4 + index as i32,
                24 + index as i32,
                0,
                Some(BlastCompoAlignmentContext::EditScript(vec![
                    crate::common::GapEditOp::Sub(20),
                ])),
            );
            align.next = head;
            head = Some(align);
        }
        head
    };
    for count in [1, 2, 3, 5] {
        assert_eq!(
            super::super::x_redo_fit::x_alignment_count(alignments(count).as_deref()),
            count
        );
        let fitted = hsplist_from_distinct_alignments(alignments(count), &template_hit, 9, 0);
        assert_eq!(fitted.hsps.len(), count);
        assert_eq!(fitted.hsps.capacity(), count);
    }
}

#[test]
fn test_build_incoming_alignment_list_uses_truncating_scaled_score_and_query_index() {
    let preliminary_hits = vec![BlastpPreliminaryHsp {
        query_start: 10,
        query_end: 40,
        subject_start: 20,
        subject_end: 50,
        gapped_query_start: 25,
        gapped_subject_start: 35,
        raw_score: 7,
        query_context: 7,
        query_frame: 0,
        subject_frame: -1,
        query_length: 100,
        q_idx: 3,
        s_idx: 9,
    }];

    let alignments =
        build_incoming_alignment_list(&preliminary_hits, 3.2).expect("incoming alignment list");
    assert_eq!(alignments.score, 22);
    assert_eq!(alignments.query_index, 7);
    assert_eq!(alignments.query_start, 10);
    assert_eq!(alignments.query_end, 40);
    assert_eq!(alignments.match_start, 20);
    assert_eq!(alignments.match_end, 50);
    assert_eq!(alignments.frame, -1);
    assert_eq!(
        alignments.context,
        Some(BlastCompoAlignmentContext::PreliminaryHspIndex(0))
    );
}

/// A context of the query's own values (the last query of a batch).
fn own_context(
    query_length: i32,
    length_adjustment: i32,
    eff_searchsp: f64,
) -> EvalueAdjustContext {
    EvalueAdjustContext {
        query_length,
        length_adjustment,
        eff_searchsp,
    }
}

#[test]
fn test_adjust_evalues_for_composition_matches_ncbi_sequence_scaling() {
    // The context is not the query's own (a query of 100 letters, adjustment 10, space
    // 9000): the adjustment reads the context's three values only.
    let context = own_context(70, 20, 6000.0);
    let mut hsp_list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(80, 1, 20, 5, 24))],
        best_evalue: f64::MAX,
    };
    hsp_list.hsps[0].e_value = 2.5;

    adjust_evalues_for_composition(&mut hsp_list, 0.25, &context, 40, 1.0);

    let query_eff = 50.0;
    let subject_eff = 20.0;
    let db_to_sequence_scale = subject_eff / (6000.0 / query_eff);
    let expected = p_to_e(blast_overall_p_value(
        0.25,
        e_to_p(2.5 * db_to_sequence_scale),
    )) / db_to_sequence_scale;
    assert!((hsp_list.hsps[0].e_value - expected).abs() < 1e-12);
    assert!((hsp_list.best_evalue - expected).abs() < 1e-12);
}

/// A query context of `aa_len` letters (the other fields are not read by
/// `evalue_adjust_contexts`).
fn adjust_test_context(q_idx: u32, aa_len: usize, is_valid: bool) -> QueryContext {
    QueryContext {
        q_idx,
        frame: 0,
        aa_seq: Vec::new(),
        aa_seq_nomask: None,
        aa_len,
        orig_len: aa_len,
        frame_base: 0,
        is_valid,
        karlin_params: KarlinParams::default(),
    }
}

fn adjust_test_space(effective_space: f64, length_adjustment: i64) -> SearchSpace {
    SearchSpace {
        effective_space,
        length_adjustment,
    }
}

fn adjust_test_ids(count: usize) -> Vec<std::sync::Arc<[u8]>> {
    (0..count)
        .map(|q| std::sync::Arc::from(format!("id{q}").into_bytes()))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3543,3590-3592,3680
// (cited at `evalue_adjust_contexts`): two batches of two queries. The first query of a
// batch reads the second's values, the last query of a batch its own.
#[test]
fn evalue_adjust_contexts_are_the_next_querys_of_the_batch_or_the_querys_own() {
    let contexts = [
        adjust_test_context(0, 100, true),
        adjust_test_context(1, 200, true),
        adjust_test_context(2, 300, true),
        adjust_test_context(3, 400, true),
    ];
    let spaces = [
        adjust_test_space(1000.0, 10),
        adjust_test_space(2000.0, 20),
        adjust_test_space(3000.0, 30),
        adjust_test_space(4000.0, 40),
    ];
    let batches = [0..2, 0..2, 2..4, 2..4];
    let ids = adjust_test_ids(4);
    let blastp = UnifiedPRefusal {
        program: "BLASTP",
        query_ids: &ids,
    };
    for refuse in [None, Some(blastp)] {
        let adjust = evalue_adjust_contexts(&contexts, &spaces, &batches, 0, refuse)
            .expect("every next query has a search space");
        assert_eq!(
            adjust,
            [
                own_context(200, 20, 2000.0),
                own_context(200, 20, 2000.0),
                own_context(400, 40, 4000.0),
                own_context(400, 40, 4000.0),
            ]
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:773-784,846-847 (cited at
// `evalue_adjust_contexts`): a next query that is not valid keeps a length adjustment of 0
// and the options' search space. With `-searchsp` that is the given space; without it the
// space is 0, which unified P-values refuse (the error is the batch's, and its text names
// the query and the one after it) and a search without them never reads.
#[test]
fn evalue_adjust_contexts_of_a_next_query_that_is_not_valid() {
    let contexts = [
        adjust_test_context(0, 50, true),
        adjust_test_context(1, 100, true),
        adjust_test_context(2, 20, false),
        adjust_test_context(3, 0, true),
    ];
    // LOSAT's spaces do not read the validity: the producer must not take these values.
    let spaces = [
        adjust_test_space(500.0, 5),
        adjust_test_space(1000.0, 10),
        adjust_test_space(777.0, 7),
        adjust_test_space(888.0, 8),
    ];
    let batches = [0..1, 1..4, 1..4, 1..4];
    let ids = adjust_test_ids(4);
    let blastp = Some(UnifiedPRefusal {
        program: "BLASTP",
        query_ids: &ids,
    });
    let given = evalue_adjust_contexts(&contexts, &spaces, &batches, 123_456, blastp)
        .expect("-searchsp gives the next query a search space");
    assert_eq!(
        given,
        [
            own_context(50, 5, 500.0),
            own_context(20, 0, 123_456.0),
            own_context(0, 0, 123_456.0),
            own_context(0, 0, 123_456.0),
        ]
    );
    let plain = evalue_adjust_contexts(&contexts, &spaces, &batches, 0, None)
        .expect("no refusal without unified P-values");
    assert_eq!(plain[1], own_context(20, 0, 0.0));
    let error = evalue_adjust_contexts(&contexts, &spaces, &batches, 0, blastp)
        .expect_err("the valid query 1 is followed by a query without a search space");
    let failure = error
        .downcast::<crate::algorithm::common::database::HeaderReadError>()
        .expect("the batch's error");
    assert_eq!(failure.query, 1);
    let text = failure.to_string();
    assert!(text.contains("are not supported by LOSAT's BLASTP for query id1, which"));
    assert!(text.contains("(id2: no residue that can be scored"));
    // A query that is not valid itself has no match to adjust: queries 2 and 3 alone in a
    // batch after query 1's batch are not refused.
    let batches = [0..1, 1..2, 2..4, 2..4];
    evalue_adjust_contexts(&contexts, &spaces, &batches, 0, blastp)
        .expect("the query before the one without a search space has no match");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:153-154 (cited at
// `adjust_evalues_for_composition`): both effective lengths are at least 1.
#[test]
fn test_adjust_evalues_for_composition_keeps_effective_lengths_of_one() {
    // The context's adjustment (150) is above its length (100) and the subject's (40).
    let context = own_context(100, 150, 6000.0);
    let mut hsp_list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![
            BlastpHsp::from_hit(make_hit(80, 1, 20, 5, 24)),
            BlastpHsp::from_hit(make_hit(40, 30, 40, 26, 36)),
        ],
        best_evalue: f64::MAX,
    };
    hsp_list.hsps[0].e_value = 0.0;
    hsp_list.hsps[1].e_value = 3.0;

    adjust_evalues_for_composition(&mut hsp_list, 0.5, &context, 40, 1.0);

    let db_to_sequence_scale = 1.0 / (6000.0 / 1.0);
    let expected = p_to_e(blast_overall_p_value(
        0.5,
        e_to_p(3.0 * db_to_sequence_scale),
    )) / db_to_sequence_scale;
    // An E-value of 0 stays 0 (NCBI prints 0.0).
    assert_eq!(hsp_list.hsps[0].e_value, 0.0);
    assert_eq!(hsp_list.hsps[1].e_value, expected);
    assert_eq!(hsp_list.best_evalue, 0.0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:102-116
#[test]
fn test_normalize_scores_uses_blast_nint_and_scaled_bit_score() {
    let mut hsp_list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(9, 1, 5, 1, 5))],
        best_evalue: f64::MAX,
    };

    normalize_scores(&mut hsp_list, 0.25, 0.0, 2.0);

    assert_eq!(hsp_list.hsps[0].raw_score, 5);
    let expected = (5.0_f64 * 0.25 * 2.0) / std::f64::consts::LN_2;
    assert!((hsp_list.hsps[0].bit_score - expected).abs() < 1e-12);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1811-1893
#[test]
fn test_evaluate_and_purge_uses_spouge_for_gapped_blastp() {
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    let gumbel = BlastGumbelBlk {
        lambda: 0.267,
        c: 0.669720,
        g: 12.0,
        a: 1.9,
        alpha: 42.602800,
        sigma: 43.636200,
        a_un: 0.7916,
        alpha_un: 4.964660,
        b: 2.0 * 12.0 * (0.7916 - 1.9),
        beta: 2.0 * 12.0 * (4.964660 - 42.602800),
        tau: 2.0 * 12.0 * (4.964660 - 43.636200),
        db_length: 5000,
    };
    let mut hsp_list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(1062, 1, 10, 1, 10))],
        best_evalue: f64::MAX,
    };
    let query_workspace =
        build_query_workspace(&[1u8; 100], 100, 0, 100000.0, own_context(100, 0, 100000.0));

    let (_best_score, best_evalue) = evaluate_and_purge(
        &mut hsp_list,
        &params,
        Some(&gumbel),
        &query_workspace,
        200,
        10.0,
        None,
        None,
        1.0,
        false,
    );

    let expected = blast_spouge_stoe(1062, &params, &gumbel, 100, 200);
    assert_eq!(hsp_list.hsps.len(), 1);
    assert!((hsp_list.hsps[0].e_value - expected).abs() < 1e-12);
    assert!((best_evalue - expected).abs() < 1e-12);
}

#[test]
fn test_evaluate_and_purge_uses_scaled_lambda_for_spouge() {
    let params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    let gumbel = BlastGumbelBlk {
        lambda: 0.267,
        c: 0.669720,
        g: 12.0,
        a: 1.9,
        alpha: 42.602800,
        sigma: 43.636200,
        a_un: 0.7916,
        alpha_un: 4.964660,
        b: 2.0 * 12.0 * (0.7916 - 1.9),
        beta: 2.0 * 12.0 * (4.964660 - 42.602800),
        tau: 2.0 * 12.0 * (4.964660 - 43.636200),
        db_length: 5000,
    };
    let mut hsp_list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: vec![BlastpHsp::from_hit(make_hit(1062, 1, 10, 1, 10))],
        best_evalue: f64::MAX,
    };
    let query_workspace =
        build_query_workspace(&[1u8; 100], 100, 0, 100000.0, own_context(100, 0, 100000.0));

    let (_best_score, best_evalue) = evaluate_and_purge(
        &mut hsp_list,
        &params,
        Some(&gumbel),
        &query_workspace,
        200,
        10.0,
        None,
        None,
        32.0,
        false,
    );

    assert_eq!(hsp_list.hsps.len(), 1);
    let expected = blast_spouge_stoe(
        1062,
        &scaled_karlin_params(&params, 32.0),
        &gumbel,
        100,
        200,
    );
    assert!((hsp_list.hsps[0].e_value - expected).abs() < 1e-12);
    assert!((best_evalue - expected).abs() < 1e-12);
}

#[test]
fn test_redo_preliminary_blastp_hit_preserves_internal_offsets_with_ranges() {
    let full_query = BlastCompoSequenceData::from_ncbistdaa(&[1u8; 12]);
    let full_subject = BlastCompoSequenceData::from_ncbistdaa(&[1u8; 11]);
    let query_range = BlastCompoSequenceRange {
        begin: 4,
        end: 12,
        context: 0,
    };
    let subject_range = BlastCompoSequenceRange {
        begin: 3,
        end: 11,
        context: 0,
    };
    let query_data = full_query.slice_range(&query_range);
    let subject_data = full_subject.slice_range(&subject_range);
    let preliminary_hsp = BlastpPreliminaryHsp {
        query_start: 4,
        query_end: 12,
        subject_start: 3,
        subject_end: 11,
        gapped_query_start: 4,
        gapped_subject_start: 3,
        raw_score: 0,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 12,
        q_idx: 0,
        s_idx: 0,
    };
    let mut gap_scratch = GapAlignScratch::new();
    let redone = redo_preliminary_blastp_hit(
        &preliminary_hsp,
        &query_data,
        &query_range,
        &subject_data,
        &subject_range,
        ScoringMatrix::Blosum62,
        None,
        11,
        1,
        50,
        &mut gap_scratch,
    )
    .expect("traceback should produce a hit");

    assert_eq!(redone.query_start, 4);
    assert_eq!(redone.subject_start, 3);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760-1764: callback state retains zero-based signed offsets.
    assert_eq!(redone.query_start + 1, 5);
    assert_eq!(redone.subject_start + 1, 4);
}

#[test]
fn test_subject_seg_masks_when_enabled_and_not_near_identical() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1625-1637
    let matching_residues = vec![1u8; 80];
    let matching_seq = BlastCompoMatchingSequence::new(0, &matching_residues);
    let subject_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let params = test_redo_align_params();

    let query = BlastCompoSequenceData::from_ncbistdaa(&[2u8; 80]);
    let query_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let align = blast_compo_alignment_new(
        80,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        80,
        0,
        0,
        80,
        0,
        None,
    );
    let range = build_subject_range_data(
        &matching_seq,
        &subject_range,
        &query,
        &query_range,
        Some(&build_query_word_hashes(query.data())),
        &align,
        false,
        true,
        &params,
        None,
    );
    assert!(range.subject_maybe_biased);
    assert!(range.subject.data().contains(&ncbistdaa::X));
}

// The masked range carries the composition of its bytes,
// on the miss that masks it and on the cache hit that returns it.
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:934-936
#[test]
fn x_masked_subject_range_carries_the_composition_of_its_bytes() {
    let mut matching_residues = vec![1u8; 80];
    for (i, residue) in matching_residues.iter_mut().enumerate().skip(40) {
        *residue = 1 + (i % 19) as u8;
    }
    let matching_seq = BlastCompoMatchingSequence::new(0, &matching_residues);
    let subject_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let params = test_redo_align_params();
    let query = BlastCompoSequenceData::from_ncbistdaa(&[2u8; 80]);
    let query_range = subject_range;
    let align = blast_compo_alignment_new(
        80,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        80,
        0,
        0,
        80,
        0,
        None,
    );
    let words = build_query_word_hashes(query.data());
    let mut cache = BlastpKappaSubjectRangeCache::new(None);
    let get = |cache: &mut BlastpKappaSubjectRangeCache| {
        build_subject_range_data(
            &matching_seq,
            &subject_range,
            &query,
            &query_range,
            Some(&words),
            &align,
            false,
            true,
            &params,
            Some(cache),
        )
    };
    let miss = get(&mut cache);
    let hit = get(&mut cache);
    assert!(miss.subject.data().contains(&ncbistdaa::X));
    assert_eq!(hit.subject, miss.subject);
    for range in [&miss, &hit] {
        let expected = XSubjectComposition(read_aa_composition(range.subject.data()));
        assert_eq!(range.x_subject_composition.as_deref(), Some(&expected));
        assert_eq!(x_whole_subject_composition(range), expected.0);
    }
}

#[test]
fn test_subject_seg_skips_exact_near_identical_subject() {
    let residues = vec![1u8; 80];
    let matching_seq = BlastCompoMatchingSequence::new(0, &residues);
    let subject_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let query = BlastCompoSequenceData::from_ncbistdaa(&residues);
    let query_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let params = test_redo_align_params();
    let align = blast_compo_alignment_new(
        80,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        80,
        0,
        0,
        80,
        0,
        None,
    );

    let range = build_subject_range_data(
        &matching_seq,
        &subject_range,
        &query,
        &query_range,
        Some(&build_query_word_hashes(query.data())),
        &align,
        true,
        true,
        &params,
        None,
    );
    assert!(range.subject_maybe_biased);
    assert!(!range.subject.data().contains(&ncbistdaa::X));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1636
#[test]
fn test_subject_seg_preserves_prior_not_biased_state() {
    let residues = vec![1u8; 80];
    let matching_seq = BlastCompoMatchingSequence::new(0, &residues);
    let subject_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let query = BlastCompoSequenceData::from_ncbistdaa(&residues);
    let query_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let params = test_redo_align_params();
    let align = blast_compo_alignment_new(
        80,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        80,
        0,
        0,
        80,
        0,
        None,
    );

    let range = build_subject_range_data(
        &matching_seq,
        &subject_range,
        &query,
        &query_range,
        Some(&build_query_word_hashes(query.data())),
        &align,
        false,
        false,
        &params,
        None,
    );
    assert!(!range.subject_maybe_biased);
    assert!(!range.subject.data().contains(&ncbistdaa::X));
}

#[test]
fn test_subject_range_slice_has_fresh_leading_sentinel() {
    let matching_seq = BlastCompoMatchingSequence::new(0, &[1u8, 2, 3, 4, 5, 6]);
    let subject_range = BlastCompoSequenceRange {
        begin: 2,
        end: 5,
        context: 0,
    };
    let query = BlastCompoSequenceData::from_ncbistdaa(&[2u8, 3, 4, 5, 6, 7]);
    let query_range = BlastCompoSequenceRange {
        begin: 0,
        end: 6,
        context: 0,
    };
    let params = test_redo_align_params();
    let align = blast_compo_alignment_new(
        42,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        6,
        0,
        0,
        6,
        0,
        None,
    );

    let range = build_subject_range_data(
        &matching_seq,
        &subject_range,
        &query,
        &query_range,
        None,
        &align,
        false,
        true,
        &params,
        None,
    );

    assert_eq!(range.subject.data(), &[3u8, 4, 5]);
    assert_eq!(range.subject.buffer[range.subject.data_offset - 1], 0);
    assert_eq!(range.subject.buffer.len(), 5);
}

#[test]
fn test_subject_range_near_identical_uses_original_range_offset() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1267-1277
    let mut residues = vec![9u8; 4];
    residues.extend(std::iter::repeat_n(1u8, 80));
    residues.extend(std::iter::repeat_n(9u8, 4));
    let matching_seq = BlastCompoMatchingSequence::new(0, &residues);
    let subject_range = BlastCompoSequenceRange {
        begin: 4,
        end: 84,
        context: 0,
    };
    let query = BlastCompoSequenceData::from_ncbistdaa(&[1u8; 80]);
    let query_range = BlastCompoSequenceRange {
        begin: 0,
        end: 80,
        context: 0,
    };
    let params = test_redo_align_params();
    let align = blast_compo_alignment_new(
        80,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        80,
        0,
        4,
        84,
        0,
        None,
    );

    let range = build_subject_range_data(
        &matching_seq,
        &subject_range,
        &query,
        &query_range,
        Some(&build_query_word_hashes(query.data())),
        &align,
        true,
        true,
        &params,
        None,
    );

    assert!(range.subject_maybe_biased);
    assert!(!range.subject.data().contains(&ncbistdaa::X));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1682-1704
#[test]
#[should_panic(expected = "NCBI blast_kappa.c requires non-negative query_range.begin")]
fn test_query_range_data_rejects_negative_begin_instead_of_clamping() {
    let query = BlastCompoSequenceData::from_ncbistdaa(&[2u8, 3, 4, 5, 6, 7]);
    let query_range = BlastCompoSequenceRange {
        begin: -1,
        end: 3,
        context: 0,
    };

    let _ = build_query_range_data(&query, &query_range);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1638-1643
#[test]
#[should_panic(expected = "NCBI blast_kappa.c requires non-negative subject_range.begin")]
fn test_subject_range_data_rejects_negative_begin_instead_of_clamping() {
    let matching_seq = BlastCompoMatchingSequence::new(0, &[1u8, 2, 3, 4, 5, 6]);
    let subject_range = BlastCompoSequenceRange {
        begin: -1,
        end: 5,
        context: 0,
    };
    let query = BlastCompoSequenceData::from_ncbistdaa(&[2u8, 3, 4, 5, 6, 7]);
    let query_range = BlastCompoSequenceRange {
        begin: 0,
        end: 6,
        context: 0,
    };
    let params = test_redo_align_params();
    let align = blast_compo_alignment_new(
        42,
        EMatrixAdjustRule::DontAdjustMatrix,
        0,
        6,
        0,
        0,
        6,
        0,
        None,
    );

    let _ = build_subject_range_data(
        &matching_seq,
        &subject_range,
        &query,
        &query_range,
        None,
        &align,
        false,
        true,
        &params,
        None,
    );
}

#[test]
fn test_reap_contained_hits_discards_lower_scoring_inner_hit() {
    let hits = vec![
        RedoneBlastpHit {
            hit: make_hit(100, 11, 50, 101, 140),
        },
        RedoneBlastpHit {
            hit: make_hit(90, 21, 40, 111, 130),
        },
        RedoneBlastpHit {
            hit: make_hit(80, 61, 80, 151, 170),
        },
    ];
    let mut hsp_list = hsp_list_of(&hits);
    reap_contained_hits(&mut hsp_list);
    assert_eq!(hsp_list.hsps.len(), 2);
    assert_eq!(hsp_list.hsps[0].raw_score, 100);
    assert_eq!(hsp_list.hsps[1].raw_score, 80);
}

#[test]
fn test_recompute_num_identities_ungapped() {
    let mut hsps = vec![BlastpHsp::from_hit(make_hit(20, 1, 4, 1, 4))];
    recompute_num_identities(
        &[1u8, 2, 3, 4],
        &[1u8, 2, 5, 4],
        &mut hsps,
        ScoringMatrix::Blosum62,
    );
    assert_eq!(hsps[0].num_ident, 3);
    assert_eq!(hsps[0].length, 4);
    assert_eq!(hsps[0].mismatch, 1);
}

#[test]
fn test_recompute_num_identities_uses_query_nomask() {
    let mut hsps = vec![BlastpHsp::from_hit(make_hit(20, 1, 4, 1, 4))];
    recompute_num_identities(
        &[1u8, 2, 3, 4],
        &[1u8, 2, 3, 4],
        &mut hsps,
        ScoringMatrix::Blosum62,
    );
    assert_eq!(hsps[0].num_ident, 4);
    assert_eq!(hsps[0].num_positives, 4);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:248-259
#[test]
fn test_reap_contained_hits_requires_matching_subject_frame() {
    let mut hits = vec![
        RedoneBlastpHit {
            hit: make_hit(100, 1, 40, 101, 140),
        },
        RedoneBlastpHit {
            hit: {
                let mut hit = make_hit(90, 11, 30, 111, 130);
                hit.query_frame = 1;
                hit
            },
        },
    ];
    hits[0].hit.query_frame = 1;
    let mut hsp_list = hsp_list_of(&hits);
    hsp_list.hsps[0].subject_frame = 0;
    hsp_list.hsps[1].subject_frame = 1;

    reap_contained_hits(&mut hsp_list);

    assert_eq!(hsp_list.hsps.len(), 2);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760-1772,1934-1946
// Direct unit boundary only: this synthetic sentinel tests callback state,
// separate from natural FASTA runtime reachability. The exact pinned C
// s_NewAlignmentFromGapAlign owner provides every expected coordinate,
// score and complete script for zero and nonzero range origins.
#[test]
fn test_blastp_fence_signed_caller_matches_pinned_ncbi() {
    let expected = include_str!("../../../tests/unit/blastp_fence_signed_caller_expected.tsv");
    let mut observed = String::from("qstart\tsstart\tfence_position\tquery_origin\tsubject_origin\tscore\tquery_start\tquery_stop\tsubject_start\tsubject_stop\tedits\n");
    for line in expected.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        let qstart: usize = f[0].parse().unwrap();
        let sstart: usize = f[1].parse().unwrap();
        let position: i32 = f[2].parse().unwrap();
        let qb: i32 = f[3].parse().unwrap();
        let sb: i32 = f[4].parse().unwrap();
        let query_data = BlastCompoSequenceData::from_ncbistdaa(&[1u8; 6]);
        let mut sequence = [1u8; 6];
        if position >= 0 {
            sequence[usize::try_from(position).unwrap()] = 201;
        }
        let subject_data = BlastCompoSequenceData::from_ncbistdaa(&sequence);
        let qrange = BlastCompoSequenceRange {
            begin: qb,
            end: qb + 6,
            context: 0,
        };
        let srange = BlastCompoSequenceRange {
            begin: sb,
            end: sb + 6,
            context: 0,
        };
        let mut scratch = GapAlignScratch::new();
        let input = BlastpPreliminaryHsp {
            query_start: qb,
            query_end: qb + 6,
            subject_start: sb,
            subject_end: sb + 6,
            gapped_query_start: qb + qstart as i32,
            gapped_subject_start: sb + sstart as i32,
            raw_score: 0,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 6,
            q_idx: 0,
            s_idx: 0,
        };
        let result = redo_preliminary_blastp_hit(
            &input,
            &query_data,
            &qrange,
            &subject_data,
            &srange,
            ScoringMatrix::Blosum62,
            None,
            11,
            1,
            18,
            &mut scratch,
        )
        .expect("successful NCBI callback preserves signed payload");
        let score = result.score;
        let script = &result.edit_script;
        let (qs, qe, ss, se) = (
            result.query_start,
            result.query_end,
            result.subject_start,
            result.subject_end,
        );
        let edits = if script.is_empty() {
            "-".to_owned()
        } else {
            script
                .iter()
                .map(|op| {
                    let kind = match op {
                        GapEditOp::Sub(_) => 3,
                        GapEditOp::Del(_) => 0,
                        GapEditOp::Ins(_) => 6,
                    };
                    format!("{kind}:{}", op.num())
                })
                .collect::<Vec<_>>()
                .join(",")
        };
        let payload = format!("{score}\t{qs}\t{qe}\t{ss}\t{se}\t{edits}");
        assert_eq!(
            payload,
            f[7..].join("\t"),
            "direct callback case {}",
            f[..5].join("/")
        );
        observed.push_str(&format!("{}\t{payload}\n", f[..5].join("\t")));
    }
    if let Some(root) = std::env::var_os("LOSAT_FENCE_CALLER_DUMP") {
        std::fs::write(
            std::path::Path::new(&root).join("blastp_caller_actual.tsv"),
            observed,
        )
        .unwrap();
    }
}
