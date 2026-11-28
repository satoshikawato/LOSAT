use super::*;
use crate::utils::glibc_fma_libm::GlibcLibm;

/// BLOSUM62 ungapped (kbp_std).
fn blosum62_ungapped() -> KarlinParams {
    KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.7916,
        beta: -3.2,
    }
}

/// BLOSUM62 11/1 (kbp_gap).
fn blosum62_gapped() -> KarlinParams {
    KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    }
}

/// BLASTN 2/-3 (ungapped and 5/2 gapped share it).
fn blastn_2_3() -> KarlinParams {
    KarlinParams {
        lambda: 0.625,
        k: 0.41,
        h: 0.78,
        alpha: 0.0,
        beta: 0.0,
    }
}

/// The EtoS of the protein programs before this owner (`tblastx::ncbi_cutoffs::
/// cutoff_score_from_evalue` of a32dc3fe6b, the port shared since line PS DECISIONS 28), for the
/// bit comparison.
fn port_before(evalue: f64, searchsp: i64, kbp: &KarlinParams) -> i32 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return BLAST_SCORE_MIN;
    }
    let e = evalue.max(1.0e-297);
    let searchsp = searchsp as f64;
    let q = kbp.k * searchsp / e;
    let log_value = if q.is_finite() {
        q.glibc_ln()
    } else {
        kbp.k.glibc_ln() + searchsp.glibc_ln() - e.glibc_ln()
    };
    (log_value / kbp.lambda).ceil() as i32
}

/// The decay of the port before this owner (`cutoff_score_from_evalue_with_decay`), from 1.
fn port_with_decay_before(evalue: f64, searchsp: i64, rate: f64, kbp: &KarlinParams) -> i32 {
    let e = if rate > 0.0 && rate < 1.0 {
        evalue * (1.0 - rate)
    } else {
        evalue
    };
    port_before(e, searchsp, kbp).max(1)
}

/// BLASTN's hit-saving copy before this owner (`blastn::ncbi_cutoffs::cutoff_score_max_from_evalue`
/// of a32dc3fe6b).
fn blastn_before(evalue: f64, searchsp: i64, kbp: &KarlinParams) -> i32 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return 1;
    }
    let mut s = 1;
    if evalue > 0.0 {
        let e = evalue.max(1.0e-297);
        let q = kbp.k * (searchsp as f64) / e;
        let es = ncbi_int4_from_double((q.glibc_ln() / kbp.lambda).ceil());
        if es > s {
            s = es;
        }
    }
    s
}

/// BLASTN's ungapped word cutoff copy before this owner (`blastn::ncbi_cutoffs::
/// compute_blastn_cutoff_score_ungapped` of a32dc3fe6b) before its cap: its one e-value is
/// `CUTOFF_E_BLASTN` (0.05), which it raised to `kSmallFloat` before the gap decay (NCBI decays
/// first, blast_stat.c:4137-4144, then raises, :4049; the owner does the same).
fn blastn_ungapped_before(searchsp: i64, rate: f64, kbp: &KarlinParams) -> i32 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return 1;
    }
    let mut s = 1;
    let e = 0.05f64.max(1.0e-297);
    let e = if rate > 0.0 && rate < 1.0 {
        e * (1.0 - rate)
    } else {
        e
    };
    let q = kbp.k * (searchsp as f64) / e;
    let es = ncbi_int4_from_double((q.glibc_ln() / kbp.lambda).ceil());
    if es > s {
        s = es;
    }
    s
}

const EVALUES: [f64; 12] = [
    1e-300, 1e-297, 1e-200, 1e-50, 1e-10, 1e-3, 0.05, 1.0, 10.0, 1000.0, 1e6, 1e300,
];
const SPACES: [i64; 12] = [
    0,
    1,
    9,
    25,
    10_000,
    1_000_000,
    123_456_789,
    40_000_000_001,
    1_000_000_000_000,
    4_000_000_000_000_000,
    i64::MAX / 3,
    i64::MAX,
];

// The owner gives the bits of the code it replaced on valid blocks: the hit-saving and word
// cutoffs of the protein programs (the port raised to 1, with and without the decay) and BLASTN's
// two copies (x86-64 `INT_MIN`; the ungapped one at its only e-value, where its order of the
// clamp and the decay, the reverse of NCBI's, cannot change a bit).
#[test]
fn the_owner_gives_the_cutoffs_of_the_copies_it_replaced() {
    for kbp in [blosum62_ungapped(), blosum62_gapped(), blastn_2_3()] {
        for &e in &EVALUES {
            for &ss in &SPACES {
                assert_eq!(
                    cutoff_score(e, &kbp, ss, None, CutoffConversion::FiniteQuotient),
                    port_before(e, ss, &kbp).max(1),
                    "{e} {ss} {kbp:?}"
                );
                for rate in [0.0, 0.1, 0.5] {
                    assert_eq!(
                        cutoff_score(e, &kbp, ss, Some(rate), CutoffConversion::FiniteQuotient),
                        port_with_decay_before(e, ss, rate, &kbp),
                        "{e} {ss} {rate} {kbp:?}"
                    );
                }
                assert_eq!(
                    cutoff_score(e, &kbp, ss, None, CutoffConversion::X86Int4),
                    blastn_before(e, ss, &kbp),
                    "{e} {ss} {kbp:?}"
                );
            }
        }
        for &ss in &SPACES {
            for rate in [0.0, 0.1, 0.5] {
                assert_eq!(
                    cutoff_score(0.05, &kbp, ss, Some(rate), CutoffConversion::X86Int4),
                    blastn_ungapped_before(ss, rate, &kbp),
                    "{ss} {rate} {kbp:?}"
                );
            }
        }
    }
}

// blast_stat.c:4059-4061 (TBLASTX: PD-LOSAT-NCBI-DEFECTS exception 8; the other protein
// programs: line PS DECISIONS 28): a quotient
// `K searchsp / E` past the double takes `ln K + ln searchsp - ln E` under `FiniteQuotient`;
// under `X86Int4` NCBI's x86-64 `INT_MIN` keeps the cutoff of 1.
#[test]
fn a_quotient_past_the_double_follows_the_callers_rule() {
    let gapped = blosum62_gapped();
    let expected = ((0.041f64.glibc_ln() + (i64::MAX as f64).glibc_ln() - 1.0e-297f64.glibc_ln())
        / 0.267)
        .ceil();
    assert_eq!(expected, 2713.0);
    for evalue in [1e-300, 1e-297] {
        assert_eq!(
            karlin_etos_simple(evalue, &gapped, i64::MAX, CutoffConversion::FiniteQuotient),
            2713
        );
        assert_eq!(
            karlin_etos_simple(evalue, &gapped, i64::MAX, CutoffConversion::X86Int4),
            i32::MIN
        );
        assert_eq!(
            cutoff_score(evalue, &gapped, i64::MAX, None, CutoffConversion::X86Int4),
            1
        );
    }
    // The edge: the largest search space whose quotient at E = 1e-297 is finite, and the next.
    let last_finite = (f64::MAX * 1.0e-297 / 0.041) as i64;
    assert!((0.041 * last_finite as f64 / 1.0e-297).is_finite());
    let finite = karlin_etos_simple(
        1e-297,
        &gapped,
        last_finite,
        CutoffConversion::FiniteQuotient,
    );
    assert_eq!(
        finite,
        ((0.041 * last_finite as f64 / 1.0e-297).glibc_ln() / 0.267).ceil() as i32
    );
    assert_eq!(
        karlin_etos_simple(1e-297, &gapped, last_finite, CutoffConversion::X86Int4),
        finite
    );
    let past = (last_finite as f64 * (1.0 + 1e-9)) as i64;
    assert!(!(0.041 * past as f64 / 1.0e-297).is_finite());
    let overflowed = karlin_etos_simple(1e-297, &gapped, past, CutoffConversion::FiniteQuotient);
    assert!(
        (finite..=finite + 1).contains(&overflowed),
        "{finite} {overflowed}"
    );
}

// blast_stat.c:4051-4057: a negative Lambda, K or H gives `BLAST_SCORE_MIN` (H = 0 is valid);
// blast_stat.c:4101-4102: `BLAST_Cutoffs` returns at once for a block of -1s.
#[test]
fn an_invalid_block_keeps_the_cutoff() {
    let gapped = blosum62_gapped();
    for conversion in [CutoffConversion::FiniteQuotient, CutoffConversion::X86Int4] {
        assert!(karlin_etos_simple(1000.0, &gapped, 1_000_000, conversion) > 1);
        for invalid in [
            KarlinParams {
                lambda: -0.267,
                ..gapped
            },
            KarlinParams {
                k: -0.041,
                ..gapped
            },
            KarlinParams { h: -0.14, ..gapped },
            KarlinParams {
                h: f64::NEG_INFINITY,
                ..gapped
            },
        ] {
            assert_eq!(
                karlin_etos_simple(1000.0, &invalid, 1_000_000, conversion),
                BLAST_SCORE_MIN
            );
            assert_eq!(
                cutoff_score(1e-10, &invalid, 1_000_000, None, conversion),
                1
            );
            // Not -1: E is written back as the e-value of the cutoff 1, -1 for such a block.
            let (mut s, mut e) = (1, 1.0);
            blast_cutoffs(&mut s, &mut e, &invalid, 1_000_000, Some(0.1), conversion);
            assert_eq!((s, e), (1, -1.0 / 0.9));
        }
        assert_eq!(
            karlin_etos_simple(
                1000.0,
                &KarlinParams { h: 0.0, ..gapped },
                1_000_000,
                conversion
            ),
            karlin_etos_simple(1000.0, &gapped, 1_000_000, conversion)
        );
        // A block of -1s: neither the cutoff nor E changes.
        let minus_one = KarlinParams {
            lambda: -1.0,
            k: -1.0,
            h: -1.0,
            ..gapped
        };
        let (mut s, mut e) = (7, 0.5);
        blast_cutoffs(&mut s, &mut e, &minus_one, 1_000_000, Some(0.1), conversion);
        assert_eq!((s, e), (7, 0.5));
    }
}

// blast_stat.c:4108-4121: E not above 0 (0, negative, NaN: `-evalue -nan` passes NCBI's option
// check) leaves the score from E at 1.
#[test]
fn an_evalue_not_above_zero_gives_one() {
    let gapped = blosum62_gapped();
    for conversion in [CutoffConversion::FiniteQuotient, CutoffConversion::X86Int4] {
        for e in [0.0, -1.0, f64::NAN, -f64::NAN, f64::NEG_INFINITY] {
            assert_eq!(
                cutoff_score(e, &gapped, 1_000_000_000, None, conversion),
                1,
                "{e}"
            );
            assert_eq!(
                cutoff_score(e, &gapped, 1_000_000_000, Some(0.1), conversion),
                1,
                "{e}"
            );
        }
        // +inf is above 0: its score is below 1 in any search space.
        assert_eq!(
            cutoff_score(f64::INFINITY, &gapped, 1_000_000_000, None, conversion),
            1
        );
    }
}

// blast_stat.c:4131-4145: E is written back (the e-value of the cutoff, divided by the decay
// divisor) when the cutoff did not change or E was not above 0, and kept otherwise.
#[test]
fn blast_cutoffs_writes_e_back_only_when_the_cutoff_stays() {
    let gapped = blosum62_gapped();
    let conversion = CutoffConversion::FiniteQuotient;
    // A large space: the cutoff rises, E is kept.
    let (mut s, mut e) = (1, 1.0);
    blast_cutoffs(&mut s, &mut e, &gapped, 160 * 2500, Some(0.1), conversion);
    assert_eq!(
        s,
        cutoff_score(1.0, &gapped, 160 * 2500, Some(0.1), conversion)
    );
    assert!(s > 1);
    assert_eq!(e, 1.0);
    // A small space (25): es = 1 does not raise the cutoff of 1; E becomes StoE(1) / 0.9.
    let (mut s, mut e) = (1, 1.0);
    blast_cutoffs(&mut s, &mut e, &gapped, 25, Some(0.1), conversion);
    assert_eq!(s, 1);
    assert_eq!(
        e.to_bits(),
        (crate::stats::karlin::evalue_from_raw_score(1, &gapped, 25.0) / (1.0 - 0.1)).to_bits()
    );
    // E not above 0 is written back too, without the decay when there is none.
    let (mut s, mut e) = (1, 0.0);
    blast_cutoffs(&mut s, &mut e, &gapped, 160 * 2500, None, conversion);
    assert_eq!(s, 1);
    assert_eq!(
        e.to_bits(),
        crate::stats::karlin::evalue_from_raw_score(1, &gapped, (160 * 2500) as f64).to_bits()
    );
}

// blast_parameters.c:949-976: `evalue_hsp` is declared once before the loop over the contexts.
// With BLOSUM62 11/1, the gapped decay 0.1 and an average subject length of 5 (a search space of
// 25), the first context's score from E = 1 is 1 (E written back), and the next contexts compute
// from that E: 2. Without the carry every context gets 1.
#[test]
fn the_sum_statistics_carry_their_evalue_across_contexts() {
    let gapped = blosum62_gapped();
    let conversion = CutoffConversion::FiniteQuotient;
    for avg_qlen in [5, 100] {
        let mut sum = SumStatsCutoffs::new(avg_qlen, 5, 0.1, conversion);
        let cutoffs: Vec<i32> = (0..4).map(|_| sum.next_valid_context(&gapped)).collect();
        assert_eq!(cutoffs, [1, 2, 2, 2]);
        // Each context alone (E = 1) gives 1.
        assert_eq!(cutoff_score(1.0, &gapped, 25, Some(0.1), conversion), 1);
    }
    // A search space where the first context's cutoff rises keeps E = 1 for every context.
    let mut sum = SumStatsCutoffs::new(160, 2500, 0.1, conversion);
    let first = sum.next_valid_context(&gapped);
    assert!(first > 1);
    assert_eq!(
        first,
        cutoff_score(1.0, &gapped, 160 * 2500, Some(0.1), conversion)
    );
    for _ in 0..5 {
        assert_eq!(sum.next_valid_context(&gapped), first);
    }
    // A search space of 3 x 3 (TBLASTN's subject length 3 of a `-db` whose shortest length is
    // 10): the score from E is below 1, and every context's cutoff is 1 (rounding gives 1 again).
    let mut sum = SumStatsCutoffs::new(160, 3, 0.1, conversion);
    for _ in 0..4 {
        assert_eq!(sum.next_valid_context(&gapped), 1);
    }
}

/// For contexts that share one block, the state after one valid context is the state after any
/// number of them (`SumStatsCutoffs::next_valid_context`: a BLASTX unit resumes the carry of its
/// reader batch with one step), whether the first cutoff stays 1 or rises.
#[test]
fn one_valid_context_gives_the_state_of_any_number() {
    let conversion = CutoffConversion::FiniteQuotient;
    for gapped in [blosum62_gapped(), blosum62_ungapped()] {
        for (avg_qlen, avg_subject_length) in [(5, 5), (100, 5), (100, 3), (160, 2500), (1, 1)] {
            for rate in [0.1, 0.5] {
                let mut once = SumStatsCutoffs::new(avg_qlen, avg_subject_length, rate, conversion);
                once.next_valid_context(&gapped);
                let mut more = once;
                for _ in 0..6 {
                    more.next_valid_context(&gapped);
                    assert_eq!(more.evalue_hsp.to_bits(), once.evalue_hsp.to_bits());
                }
            }
        }
    }
}
