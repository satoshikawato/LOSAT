use super::*;

fn context(lambda: f64, k: f64) -> QueryContext {
    QueryContext {
        q_idx: 0,
        frame: 1,
        aa_seq: Vec::new(),
        aa_seq_nomask: None,
        aa_len: 100,
        orig_len: 300,
        frame_base: 0,
        is_valid: true,
        karlin_params: KarlinParams {
            lambda,
            k,
            h: 0.4,
            ..KarlinParams::default()
        },
    }
}

/// An HSP of context `ctx_idx` with `raw_score`, from `q_start` to the `Int4` end `q_end`
/// (at or before `q_start` for an HSP of length 0 or less).
fn hit(ctx_idx: u32, raw_score: i32, q_start: u32, q_end: i32) -> UngappedHit {
    let mut hit = UngappedHit {
        q_idx: 0,
        s_idx: 0,
        ctx_idx,
        s_f_idx: 0,
        q_frame: 1,
        s_frame: 1,
        q_aa_start: q_start,
        q_aa_end: 0,
        s_aa_start: q_start,
        s_aa_end: 0,
        q_seed_off: q_start,
        s_seed_off: q_start,
        q_orig_len: 300,
        s_orig_len: 300,
        raw_score,
        e_value: 99.0,
        num_ident: 0,
        hsp_list_order: 0,
        ordering_method: 0,
        linked_set: false,
        start_of_chain: false,
        link_id: 0,
        chain_next_link_id: None,
        hsp_link_num: 0,
        num: 0,
    };
    hit.set_int4_ends(q_end, q_end);
    hit
}

fn link_args<'a>(
    contexts: &'a [QueryContext],
    linking_params: Option<&'a LinkingParams>,
    eff_searchsp_per_context: &'a [i64],
) -> SubjectLinkArgs<'a> {
    SubjectLinkArgs {
        linking_params_for_cutoff: &contexts[0].karlin_params,
        linking_params,
        contexts_ref: contexts,
        subject_frame_bases: &[0; 6],
        length_adj_per_context: &[0, 0],
        eff_searchsp_per_context,
        x_link_on_pool: false,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:774-783 (cited at
// `subject_linking_params`): only a search with sum statistics has linking parameters.
#[test]
fn only_sum_statistics_make_linking_parameters() {
    let contexts = [context(0.3176, 0.134)];
    let made = |sum_stats| {
        subject_linking_params(sum_stats, 100, 3000, 3000, 20, &contexts[0].karlin_params)
    };
    assert!(made(false).is_none());
    let linking = made(true).unwrap();
    assert_eq!(linking.subject_len_nucl, 3000);
    assert_eq!(linking.gap_decay_rate, BLAST_GAP_DECAY_RATE);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:874-890 and
// c++/src/algo/blast/core/blast_hits.c:1887-1893 (cited at `evalues_without_linking`), with
// c++/src/algo/blast/core/blast_stat.c:4157-4171 (`BLAST_KarlinStoE_simple`).
#[test]
fn hsps_without_sum_statistics_get_their_own_evalues_in_their_order() {
    let contexts = [context(0.3176, 0.134), context(0.267, 0.041)];
    let eff_searchsp = [1_000_000, 2_500_000];
    let link = link_args(&contexts, None, &eff_searchsp);
    let hits = vec![hit(1, 40, 0, 20), hit(0, 50, 5, 30), hit(0, 30, 9, 19)];
    let evaluated = first_link_subject_hits(&link, hits, 40);
    // The list keeps its order (nothing links or sorts it), and no HSP joins a linked set.
    assert_eq!(
        evaluated.iter().map(|h| h.raw_score).collect::<Vec<_>>(),
        [40, 50, 30]
    );
    for h in &evaluated {
        let ctx = h.ctx_idx as usize;
        let expected = evalue_from_raw_score(
            h.raw_score,
            &contexts[ctx].karlin_params,
            eff_searchsp[ctx] as f64,
        );
        assert_eq!(h.e_value.to_bits(), expected.to_bits());
        assert_eq!((h.num, h.linked_set), (0, false));
    }
    // searchsp * exp(-Lambda * S + log K) of each HSP's own context and search space.
    let close = |actual: f64, expected: f64| (actual / expected - 1.0).abs() < 1e-12;
    assert!(close(evaluated[0].e_value, 2.3575385266626068));
    assert!(close(evaluated[1].e_value, 0.017002329395178577));
    assert!(first_link_subject_hits(&link, Vec::new(), 40).is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:687,699 and :452-453
// (cited at `drop_unread_int4_inverted`), with c++/src/algo/blast/core/blast_engine.c:1480-1497 and
// c++/src/algo/blast/core/blast_hits.c:2651-2665: the fetch of the subject comes first.
#[test]
fn hsps_of_length_zero_or_less_stay_up_to_the_fetch_with_or_without_sum_statistics() {
    let contexts = [context(0.3176, 0.134), context(0.267, 0.041)];
    let eff_searchsp = [1_000_000, 2_500_000];
    let hits = || vec![hit(0, 30, 10, 10), hit(0, 50, 5, 30), hit(1, 25, 0, -2)];
    assert!(hits()[0].is_int4_inverted() && hits()[2].is_int4_inverted());
    let linking =
        subject_linking_params(true, 100, 3000, 3000, 20, &contexts[0].karlin_params).unwrap();
    for linking_params in [None, Some(&linking)] {
        // The first call (the E-values, or the linking with `Int4` ends) does not drop them: the
        // fetch of the subject sees them.
        let link = link_args(&contexts, linking_params, &eff_searchsp);
        let mut kept = first_link_subject_hits(&link, hits(), 0);
        assert_eq!(kept.len(), 3, "sum statistics {}", linking_params.is_some());
        drop_unread_int4_inverted(0, &mut kept);
        assert_eq!(kept.iter().map(|h| h.raw_score).collect::<Vec<_>>(), [50]);
    }
    // The two-hit word finder saves none: its list stays as it is.
    let mut two_hit = hits();
    drop_unread_int4_inverted(40, &mut two_hit);
    assert_eq!(two_hit.len(), 3);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1353 (cited at
// `sort_and_reap_after_first_link`): `ScoreCompareHSPs` reads the query's end as `Int4`, so of two
// HSPs equal up to it, the one ending at 5 comes before the one ending at -2 (which a `u32` reads
// as 4294967294); then the reap by the preliminary E-value keeps an HSP of length 0 or less.
#[test]
fn the_sort_after_the_first_link_of_a_one_hit_list_reads_int4_ends() {
    let contexts = [context(0.3176, 0.134), context(0.267, 0.041)];
    let linking =
        subject_linking_params(true, 100, 3000, 3000, 20, &contexts[0].karlin_params).unwrap();
    let link = link_args(&contexts, Some(&linking), &[1_000_000, 2_500_000]);
    let [reeval_ns, reeval_calls, linking_ns, linking_calls] = [0; 4].map(AtomicU64::new);
    let diagnostics = DiagnosticCounters::default();
    let stage = SubjectStageArgs {
        contexts_ref: &contexts,
        cutoff_scores: &[1, 1],
        evalue_threshold: 10.0,
        timing_enabled: false,
        reeval_ns: &reeval_ns,
        reeval_calls: &reeval_calls,
        linking_ns: &linking_ns,
        linking_calls: &linking_calls,
        diag_enabled: false,
        diagnostics: &diagnostics,
        debug_output_filter: false,
    };
    let times = super::super::x_search_timing::XSearchTimes::new(false);
    let list = || {
        let mut hits = vec![hit(0, 30, 2, -2), hit(0, 30, 2, 5), hit(0, 30, 3, 1)];
        for h in &mut hits {
            h.s_aa_start = 2;
            h.s_aa_end = 7;
            h.e_value = 1.0;
        }
        // Above the threshold: the reap removes this HSP of length 0 or less.
        hits[2].e_value = 100.0;
        hits
    };
    let sorted = |window| {
        let mut lap = times.subject();
        sort_and_reap_after_first_link(list(), &link, &stage, window, &mut lap)
            .iter()
            .map(UngappedHit::q_end_int4)
            .collect::<Vec<_>>()
    };
    assert_eq!(sorted(0), [5, -2]);
    // The two-hit word finder's list keeps the `u32` comparator (its ends are never below 0).
    assert_eq!(sorted(40), [-2, 5]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1353 (cited at
// `sort_and_reap_after_first_link`): on ends of 0 or more the `Int4` comparator orders as the
// `u32` one, so a one-hit list without an HSP of length 0 or less sorts as before.
#[test]
fn the_int4_sort_is_the_u32_sort_on_lists_of_positive_lengths() {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n
    };
    for _ in 0..200 {
        let hits: Vec<UngappedHit> = (0..1 + next(40))
            .map(|order| {
                let q_start = next(6) as u32;
                let mut h = hit(next(2) as u32, 20 + next(4) as i32, q_start, 0);
                let s_start = next(6) as u32;
                h.set_int4_ends(q_start as i32 + 1 + next(3) as i32, s_start as i32 + 2);
                h.s_aa_start = s_start;
                h.hsp_list_order = order as u32;
                h
            })
            .collect();
        assert!(hits.iter().all(|h| !h.is_int4_inverted()));
        let (mut by_int4, mut by_u32) = (hits.clone(), hits);
        sort_ungapped_hits_by_score_int4(&mut by_int4);
        sort_ungapped_hits_by_score_ncbi(&mut by_u32);
        let orders =
            |list: &[UngappedHit]| list.iter().map(|h| h.hsp_list_order).collect::<Vec<_>>();
        assert_eq!(orders(&by_int4), orders(&by_u32));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:955-974 (an HSP alone in its
// chain takes the chain's E-value): an HSP of length 0 or less alone in the list gets the E-value
// of an HSP of length 1 or more alone in it, with the same score and context, so the reap reads
// NCBI's value for it.
#[test]
fn a_lone_hsp_of_length_zero_or_less_gets_the_evalue_of_a_lone_positive_one() {
    let contexts = [context(0.3176, 0.134), context(0.267, 0.041)];
    let eff_searchsp = [1_000_000, 2_500_000];
    let linking =
        subject_linking_params(true, 100, 3000, 3000, 20, &contexts[0].karlin_params).unwrap();
    for linking_params in [None, Some(&linking)] {
        let link = link_args(&contexts, linking_params, &eff_searchsp);
        let inverted = first_link_subject_hits(&link, vec![hit(1, 35, 10, 8)], 0);
        let positive = first_link_subject_hits(&link, vec![hit(1, 35, 10, 20)], 0);
        assert!(inverted[0].is_int4_inverted() && !positive[0].is_int4_inverted());
        assert_eq!(
            inverted[0].e_value.to_bits(),
            positive[0].e_value.to_bits(),
            "sum statistics {}",
            linking_params.is_some()
        );
        assert!(inverted[0].e_value.is_finite() && inverted[0].e_value > 0.0);
    }
}
