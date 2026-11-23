use super::*;
use crate::config::ProteinScoringSpec;
use crate::stats::lookup_protein_gumbel_params;

/// BLOSUM62 11/1's gapped block (blast_stat.c's table: Lambda 0.267, K 0.041, H 0.14).
fn blosum62_gapped() -> KarlinParams {
    KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

#[test]
fn the_gumbel_block_gives_spouges_values() {
    let gumbel = lookup_protein_gumbel_params(&ProteinScoringSpec::default(), 3_234_567).unwrap();
    let stats = EvalueStats::spouge(&gumbel);
    let kbp = blosum62_gapped();
    for (score, m, n) in [(48, 302, 250), (25, 302, 10), (120, 87, 1_400)] {
        assert_eq!(
            stats.e_value(score, &kbp, m, n).to_bits(),
            blast_spouge_stoe(score, &kbp, &gumbel, m, n).to_bits()
        );
        // RPS-BLAST's swapped roles reach Spouge's function as given.
        assert_eq!(
            stats.e_value(score, &kbp, n, m).to_bits(),
            blast_spouge_stoe(score, &kbp, &gumbel, n, m).to_bits()
        );
    }
    for (evalue, m, n) in [(50.0, 302, 10), (10.0, 87, 1_400)] {
        assert_eq!(
            stats.hit_cutoff(evalue, &kbp, m, n),
            blast_spouge_etos(evalue, &kbp, &gumbel, m, n)
        );
    }
}

// blast_stat.c:4157-4171: `searchsp * exp(-Lambda S + log K)`, whatever the lengths; -1 for
// a block with a negative value.
#[test]
fn without_the_gumbel_block_the_e_value_reads_the_search_space() {
    let kbp = blosum62_gapped();
    let stats = EvalueStats::new(None, 10_000);
    let e = stats.e_value(48, &kbp, 302, 250);
    assert_eq!(e, stats.e_value(48, &kbp, 1, 1));
    assert!((e - 0.0011139502021664525).abs() < 1e-18, "{e:e}");
    assert_eq!(
        e.to_bits(),
        (10_000f64 * (-0.267 * 48f64 + 0.041f64.glibc_ln()).glibc_exp()).to_bits()
    );
    let max = EvalueStats::new(None, i64::MAX);
    assert_eq!(
        max.e_value(48, &kbp, 302, 250),
        9_223_372_036_854_775_807f64 * (-0.267 * 48f64 + 0.041f64.glibc_ln()).glibc_exp()
    );
    for invalid in [
        KarlinParams {
            lambda: -1.0,
            ..kbp
        },
        KarlinParams { k: -1.0, ..kbp },
        KarlinParams { h: -1.0, ..kbp },
    ] {
        assert_eq!(stats.e_value(48, &invalid, 302, 250), -1.0);
    }
}

// blast_parameters.c:915-945 and blast_stat.c:4097-4129: `ceil(ln(K searchsp / E) / Lambda)`,
// at least 1 (`new_cutoff`), whatever the lengths.
#[test]
fn without_the_gumbel_block_the_hit_cutoff_reads_the_search_space() {
    let kbp = blosum62_gapped();
    let cutoff = |search_space, evalue| {
        EvalueStats::new(None, search_space)
            .hit_cutoff(evalue, &kbp, 302, 10)
            .unwrap()
    };
    // ln(0.041 * 10000 / 10) / 0.267 = 13.2...
    assert_eq!(cutoff(10_000, 10.0), 14);
    // A search space of 1 (`-searchsp 1`) gives a negative score: the cutoff stays at 1.
    assert_eq!(cutoff(1, 10.0), 1);
    // The largest `-searchsp`: ln(0.041 * 2^63 / 10) / 0.267 = 142.9...
    assert_eq!(cutoff(i64::MAX, 10.0), 143);
    assert_eq!(
        EvalueStats::new(None, 10_000).hit_cutoff(10.0, &kbp, 1, 1_000_000),
        Ok(14)
    );
    // blast_stat.c:4054-4057, 4101-4102: a block with a negative value keeps the cutoff 1.
    for invalid in [
        KarlinParams {
            lambda: -1.0,
            ..kbp
        },
        KarlinParams { k: -0.5, ..kbp },
        KarlinParams { h: -0.14, ..kbp },
    ] {
        assert_eq!(
            EvalueStats::new(None, i64::MAX).hit_cutoff(10.0, &invalid, 302, 10),
            Ok(1)
        );
    }
    // The largest search space with an E-value below kSmallFloat: the quotient is past the
    // double, the cutoff the finite ceil((ln K + ln searchsp - ln E) / Lambda)
    // (`cutoff_score_from_evalue`, line PS DECISIONS 28).
    assert_eq!(cutoff(i64::MAX, 1e-300), 2713);
}
