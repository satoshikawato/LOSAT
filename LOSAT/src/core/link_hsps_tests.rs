use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1023-1076 (the
// average query length, the expected length subtracted from both lengths, and the cutoffs)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441 (`BLAST_Nint`)
// An average query length of 0 (a batch of one-strand or very short queries): log(0) is -inf,
// NCBI's expected length is 0, and the lengths are not reduced.
#[test]
fn an_average_query_length_of_zero_takes_ncbis_expected_length_of_zero() {
    let kbp = KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.0,
        beta: 0.0,
    };
    // Two subjects (800 and 300 residues), linking the first: x = 0.25 * log(1100 / 800) *
    // K / 0.5 * 800, floor(log(x) / Lambda) + 1 = 9. An expected length of 15 (the one of an
    // average of 1) would give 10.
    let two = calculate_link_hsp_cutoffs(&kbp, 0, 2400, 3300, true, || 41);
    assert_eq!(
        (two.cutoff_small_gap, two.cutoff_big_gap, two.gap_prob),
        (0, 9, 0.0)
    );
    // One subject: the database is the subject, y = log(1) = 0, x = 0, and the big-gap
    // cutoff is the 32-bit signed integer conversion of floor(-inf), then plus 1, giving `INT_MIN` + 1
    // (x86-64 and Rust both give `INT_MIN` for the conversion of -inf).
    let one = calculate_link_hsp_cutoffs(&kbp, 0, 10_000, 10_000, true, || 41);
    assert_eq!(
        (one.cutoff_small_gap, one.cutoff_big_gap, one.gap_prob),
        (0, i32::MIN + 1, 0.0)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:1036-1072
// One subject (the database's length is the subject's) and a short query against a short subject
// (line RG item X-short-tbx-link): an expected length of 0 makes `y_variable` 0 and the logarithm
// of `x_variable` -inf, a negative one makes them negative and the logarithm NaN. NCBI's x86-64
// conversion gives `INT_MIN` for both, so the big-gap cutoff is one more; Rust's `as i32`
// gave 0 for NaN, a cutoff of 1, before this item.
#[test]
fn a_short_query_against_one_short_subject_takes_ncbis_big_gap_cutoff_of_int_min_plus_1() {
    let kbp = KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.0,
        beta: 0.0,
    };
    // An average of 2 letters against 12 bases (4 letters): log(0.134 * 8) / 0.4012 = 0.17, an
    // expected length of 0, y = log(4 / 4) = 0.
    let zero = calculate_link_hsp_cutoffs(&kbp, 2, 12, 12, true, || 41);
    assert_eq!(
        (zero.cutoff_small_gap, zero.cutoff_big_gap, zero.gap_prob),
        (0, i32::MIN + 1, 0.0)
    );
    // An average of 1 letter against 6 bases (2 letters): log(0.134 * 2) / 0.4012 = -3.28, an
    // expected length of -3; the subject length becomes 5, above the database's 2, so
    // y = log((5 - 3) / 5) < 0 and log(x) is NaN.
    let negative = calculate_link_hsp_cutoffs(&kbp, 1, 6, 6, true, || 41);
    assert_eq!(
        (
            negative.cutoff_small_gap,
            negative.cutoff_big_gap,
            negative.gap_prob
        ),
        (0, i32::MIN + 1, 0.0)
    );
    // A database longer than the subject (another subject of 12 bases) keeps y positive and the
    // cutoff finite: y = log(8 / 4) * K / 0.5, x = 0.25 * y * 2 * 4 = 0.3715, floor(-3.12) + 1.
    let two = calculate_link_hsp_cutoffs(&kbp, 2, 12, 24, true, || 41);
    assert_eq!(
        (two.cutoff_small_gap, two.cutoff_big_gap, two.gap_prob),
        (0, -3, 0.0)
    );
}

// NCBI blast_parameters.c:788-815 with `DEFAULT_LONGEST_INTRON` 122: the values of
// `-max_intron_length` on both sides of each branch.
#[test]
fn the_longest_intron_of_a_translated_search_follows_ncbis_branches() {
    for (option, gapped) in [
        (0, Some(40)),
        (1, None),
        (2, None),
        (4, None),
        (5, Some(1)),
        (7, Some(1)),
        (8, Some(2)),
        (121, Some(39)),
        (122, Some(40)),
        (123, Some(40)),
        (125, Some(41)),
        (i32::MAX, Some(715_827_881)),
    ] {
        assert_eq!(translated_longest_intron(option, true), gapped, "{option}");
    }
    for (option, ungapped) in [
        (0, 0),
        (1, 0),
        (4, 0),
        (5, 1),
        (125, 41),
        (i32::MAX, 715_827_881),
    ] {
        assert_eq!(
            translated_longest_intron(option, false),
            Some(ungapped),
            "{option}"
        );
    }
}

fn context() -> LinkContext {
    // BLOSUM62 11/1 gapped Lambda and K (blast_stat.c), a query of 300 letters.
    LinkContext {
        lambda: 0.267,
        log_k: 0.041f64.glibc_ln(),
        query_length: 300,
        length_adjustment: 40,
        eff_searchsp: 260 * 960,
    }
}

/// The list `even_gap_link_hsps` gives for HSPs of one context.
fn linked(
    hsps: &[LinkHsp],
    ctx: LinkContext,
    subject_length: i32,
    params: &EvenGapLinkParams,
) -> Vec<LinkedHsp> {
    let mut out = Vec::new();
    even_gap_link_hsps(hsps, &|_| ctx, subject_length, params, &mut |l| out.push(l));
    out
}

fn hsp(q: (i32, i32), s: (i32, i32), score: i32) -> LinkHsp {
    LinkHsp {
        context: 0,
        subject_frame: 0,
        query_offset: q.0,
        query_end: q.1,
        subject_offset: s.0,
        subject_end: s.1,
        score,
    }
}

#[test]
fn a_lone_hsp_gets_the_e_value_of_one_segment() {
    let hsps = [hsp((10, 60), (20, 70), 80)];
    let ctx = context();
    let out = linked(&hsps, ctx, 1000, &EvenGapLinkParams::traceback());
    assert_eq!(out.len(), 1);
    assert_eq!((out[0].index, out[0].num), (0, 1));
    let xsum = 80.0 * ctx.lambda - ctx.log_k;
    let expected = (ctx.eff_searchsp as f64) * (-xsum).glibc_exp() / 0.9;
    assert_eq!(out[0].evalue, expected);
}

#[test]
fn two_hsps_in_order_on_both_sequences_link_into_one_set() {
    // The first HSP ends before the second starts on the query and on the subject.
    let hsps = [hsp((200, 260), (500, 560), 60), hsp((10, 60), (20, 70), 80)];
    let out = linked(&hsps, context(), 1000, &EvenGapLinkParams::traceback());
    // The chain starts at the upstream HSP (the head) and follows its link downstream.
    assert_eq!(out.iter().map(|l| l.index).collect::<Vec<_>>(), vec![1, 0]);
    assert!(out.iter().all(|l| l.num == 2));
    assert_eq!(out[0].evalue, out[1].evalue);
}

#[test]
fn hsps_in_opposite_orders_on_the_sequences_stay_apart() {
    let hsps = [hsp((200, 260), (20, 70), 60), hsp((10, 60), (500, 560), 80)];
    let out = linked(&hsps, context(), 1000, &EvenGapLinkParams::traceback());
    assert!(out.iter().all(|l| l.num == 1));
    // The list is in the forward order of the query starts.
    assert_eq!(out.iter().map(|l| l.index).collect::<Vec<_>>(), vec![1, 0]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:528-531
// A new query frame of the linking starts where the context's query frame or the sign of the
// subject frame changes.
#[test]
fn hsps_on_subject_frames_of_opposite_signs_link_apart() {
    // The pair of `two_hsps_in_order_on_both_sequences_link_into_one_set`, on subject frames
    // of the same sign and of opposite signs.
    let params = EvenGapLinkParams::ungapped_translated_subject(LinkCutoffs {
        cutoff_small_gap: 0,
        cutoff_big_gap: 1,
        gap_prob: 0.0,
    });
    let pair = |f1: i8, f2: i8| {
        [
            LinkHsp {
                subject_frame: f1,
                ..hsp((200, 260), (500, 560), 60)
            },
            LinkHsp {
                subject_frame: f2,
                ..hsp((10, 60), (20, 70), 80)
            },
        ]
    };
    let same = linked(&pair(1, 3), context(), 3000, &params);
    assert!(same.iter().all(|l| l.num == 2), "{same:?}");
    let apart = linked(&pair(-1, 3), context(), 3000, &params);
    assert!(apart.iter().all(|l| l.num == 1), "{apart:?}");
    // The final list is in the order of the query starts whatever the frames.
    assert_eq!(
        apart.iter().map(|l| l.index).collect::<Vec<_>>(),
        vec![1, 0]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:563-571
// A translated subject's length and the length adjustment, both in nucleotides, are divided
// by the codon length before the adjustment is subtracted (at least 1 letter is left).
#[test]
fn a_translated_subject_divides_its_length_and_the_adjustment_by_three() {
    // A set of two takes the subject length into its e-value: 3,001 nucleotides and the
    // adjustment 40 give 1,000 - 13 = 987 letters, the length 1,027 of an untranslated subject
    // with the whole adjustment subtracted.
    let cutoffs = LinkCutoffs {
        cutoff_small_gap: 0,
        cutoff_big_gap: 1,
        gap_prob: 0.0,
    };
    let hsps = [hsp((200, 260), (500, 560), 60), hsp((10, 60), (20, 70), 80)];
    let translated = linked(
        &hsps,
        context(),
        3001,
        &EvenGapLinkParams::ungapped_translated_subject(cutoffs),
    );
    let protein = |subject_length| {
        linked(
            &hsps,
            context(),
            subject_length,
            &EvenGapLinkParams::ungapped(cutoffs, false),
        )
    };
    assert!(translated.iter().all(|l| l.num == 2), "{translated:?}");
    assert_eq!(translated, protein(1027));
    assert_ne!(translated, protein(1000));
}

#[test]
fn the_preliminary_cutoff_keeps_hsps_at_or_below_it_out_of_small_gap_sets() {
    // Close HSPs (a gap of 10 letters on both sequences) with scores at the small-gap cutoff:
    // NCBI walks for a predecessor only for a score above the cutoff.
    let hsps = [hsp((70, 120), (80, 130), 41), hsp((10, 60), (20, 70), 41)];
    let out = linked(&hsps, context(), 1000, &EvenGapLinkParams::preliminary(41));
    assert!(out.iter().all(|l| l.num == 1), "{out:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:871 (the store into the `Int2`
// `num`, as at :756) and blast_stat.c:4390-4391 (`s_BlastSumP`: it returns 0 when r is below 1)
#[test]
fn a_set_of_more_than_32767_hsps_takes_the_wrapped_int2_count() {
    // 32,770 HSPs in order on both sequences link into one large-gap set; NCBI's count wraps
    // to -32766, and the sum statistics give a count below 1 the E-value 0.
    let n = 32_770;
    let hsps: Vec<LinkHsp> = (0..n)
        .map(|i| hsp((100 * i, 100 * i + 30), (100 * i, 100 * i + 30), 60))
        .collect();
    let ctx = LinkContext {
        query_length: 100 * n + 100,
        eff_searchsp: (100 * n as i64) * (100 * n as i64),
        ..context()
    };
    let params = EvenGapLinkParams {
        gap_prob: 0.0,
        gap_decay_rate: BLAST_GAP_DECAY_RATE,
        cutoff_small_gap: 0,
        cutoff_big_gap: 1,
        query_translated: false,
        subject_translated: false,
    };
    let out = linked(&hsps, ctx, 100 * n + 100, &params);
    assert_eq!(out.len(), n as usize);
    assert!(out.iter().all(|l| l.num == -32766), "{:?}", &out[..2]);
    assert!(out.iter().all(|l| l.evalue == 0.0), "{:?}", &out[..2]);
    // One HSP fewer than the wrap keeps its count.
    let out = linked(&hsps[..32_767], ctx, 100 * n + 100, &params);
    assert!(out.iter().all(|l| l.num == 32_767), "{:?}", &out[..2]);
}
