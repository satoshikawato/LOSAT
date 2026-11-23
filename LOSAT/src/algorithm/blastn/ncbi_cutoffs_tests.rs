use super::*;

#[test]
fn test_gap_trigger_raw_score() {
    // Test with blastn default params (reward=2, penalty=-3, gap_open=5, gap_extend=2)
    // Ungapped params: lambda=0.625, k=0.41
    let ungapped_params = KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 0.8,
        beta: -2.0,
    };

    let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE_NUCL, &ungapped_params);
    // Expected: (27.0 * 0.6931471805599453 + ln(0.41)) / 0.625
    // = (18.714973875118523 - 0.8915981196417836) / 0.625
    // = 17.82337575547674 / 0.625
    // = 28.517401208762784
    // ≈ 28 or 29
    assert!((28..=29).contains(&gap_trigger));
}

#[test]
fn test_cutoff_score_max_from_evalue() {
    // Test with blastn gapped params (reward=2, penalty=-3, gap_open=5, gap_extend=2)
    // Gapped params: lambda=0.625, k=0.41 (same as ungapped for blastn)
    let gapped_params = KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 0.8,
        beta: -2.0,
    };

    let query_len = 1000;
    let subject_len = 10000;
    let eff_searchsp = query_len * subject_len;
    let evalue = 10.0;

    let cutoff_score_max = cutoff_score_max_from_evalue(evalue, eff_searchsp, &gapped_params);
    // Should be a positive value
    assert!(cutoff_score_max > 0);

    // With E at kSmallFloat, K * searchsp / E overflows for a search space of 1e12:
    // NCBI's (Int4) of infinity is INT_MIN on x86-64, and the cutoff stays 1.
    assert_eq!(
        cutoff_score_max_from_evalue(1e-300, 1_000_000_000_000, &gapped_params),
        1
    );
}

// BLAST_Cutoffs (blast_stat.c:4090-4129) with NCBI's infinite and NaN e-values: an
// infinite E gives log(0) / Lambda, which (Int4) makes INT_MIN, and a NaN fails `e > 0`;
// both keep the initial cutoff 1, as the largest finite e-value does (E2g).
#[test]
fn infinite_and_nan_evalues_keep_the_initial_cutoff() {
    let gapped_params = KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 0.8,
        beta: -2.0,
    };
    for evalue in [f64::INFINITY, f64::NAN, -f64::NAN, f64::MAX] {
        assert_eq!(
            cutoff_score_max_from_evalue(evalue, 10_000_000, &gapped_params),
            1,
            "{evalue}"
        );
    }
}

#[test]
fn test_compute_blastn_cutoff_score_ungapped() {
    // Test with blastn ungapped params (reward=2, penalty=-3)
    // Ungapped params: lambda=0.625, k=0.41
    let ungapped_params = KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 0.8,
        beta: -2.0,
    };

    let query_len = 1000i64 * 2; // Both strands
    let subject_len = 10000i64;
    let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE_NUCL, &ungapped_params);
    let scale_factor = 1.0;
    let gap_decay_rate = 0.0; // No gap decay for ungapped
    let cutoff_score_max = 10000; // Dummy value for testing

    let cutoff_score = compute_blastn_cutoff_score_ungapped(
        query_len,
        subject_len,
        gap_trigger,
        &ungapped_params,
        scale_factor,
        gap_decay_rate,
        cutoff_score_max,
    );

    // Should be a positive value
    assert!(cutoff_score > 0);
    // Should use CUTOFF_E_BLASTN = 0.05
    // Search space is min(2000, 10000) * 10000 = 2000 * 10000 = 20,000,000.
    // Expected score should be calculated from CUTOFF_E_BLASTN = 0.05
    // Final cutoff is capped at `cutoff_score_max`.
    assert!(cutoff_score <= cutoff_score_max);
}
