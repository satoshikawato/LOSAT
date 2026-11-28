//! The cutoff score of an e-value: NCBI's `BLAST_Cutoffs` and its `BlastKarlinEtoS_simple`, the
//! one owner for every program (the hit-saving cutoffs, the ungapped word cutoffs and the sum
//! statistics' cutoffs of BLASTN, BLASTP, PSI-BLAST, RPS-BLAST, BLASTX, TBLASTN and TBLASTX).
//! The reverse, the e-value of a score, is `stats::karlin` (`BLAST_KarlinStoE_simple`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4041-4150
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:321,360-363
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:916,943
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:949-976

use super::karlin::evalue_from_raw_score;
use super::sum_statistics::gap_decay_divisor;
use super::tables::KarlinParams;
use crate::core::blast_util::{ncbi_int4_from_double, BLAST_SCORE_MIN};
use crate::utils::glibc_fma_libm::GlibcLibm;

/// The smallest e-value that `BlastKarlinEtoS_simple` converts (`kSmallFloat`): a smaller one is
/// raised to it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4049
pub const K_SMALL_FLOAT: f64 = 1.0e-297;

/// How `BlastKarlinEtoS_simple` turns `ceil(log(K * searchsp / E) / Lambda)` into an `Int4` where C
/// does not define the conversion: a quotient `K * searchsp / E` past the range of a double
/// (+inf), NaN, or a score past `Int4`. Each caller passes its program's rule, one of the three
/// consts below; which program takes which rule is an Owner decision, not this module's.
///
/// - `TBLASTX` takes the finite value: PD-LOSAT-NCBI-DEFECTS exception 8, the Owner's answer
///   「Q1=A」 of 2026-11-12 (line PS DECISIONS 36), which names TBLASTX only.
/// - `PROTEIN` (TBLASTN, BLASTX, BLASTP, PSI-BLAST, RPS-BLAST and RPSTBLASTN) takes the finite
///   value under line PS DECISIONS 28, Owner-delegated, which applied the Owner's rule of
///   2026-10-10 (LOSAT does not reproduce NCBI's undefined behaviour). Its known difference from
///   NCBI (BLASTX `-ungapped`: NCBI 17 rows, LOSAT 16) is the open Owner question Q-ETOS8 (line
///   RG item X-etos-owner).
/// - `BLASTN` takes NCBI's x86-64 `INT_MIN` (PD-LOSAT-NCBI-DEFECTS rule 2, as decision D4
///   reproduces the conversions of `-threshold` and the X-drops), also pending Q-ETOS8.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CutoffConversion {
    /// The finite value: a quotient past the double takes the mathematically equal
    /// `ln K + ln searchsp - ln E` (NCBI's x86-64 build gives `INT_MIN` there, so its cutoff stays
    /// 1); a finite quotient keeps NCBI's order of operations, and the score converts with
    /// Rust's saturating `as` (NaN 0).
    FiniteQuotient,
    /// NCBI's x86-64 conversion (`core::blast_util::ncbi_int4_from_double`): `INT_MIN` for NaN,
    /// infinity and any value past `Int4`, so `BLAST_Cutoffs` keeps its caller's cutoff.
    X86Int4,
}

impl CutoffConversion {
    /// BLASTN's rule: NCBI's `INT_MIN` (rule 2; Q-ETOS8 open).
    pub const BLASTN: Self = Self::X86Int4;
    /// TBLASTX's rule: the finite value (exception 8, Owner 2026-11-12).
    pub const TBLASTX: Self = Self::FiniteQuotient;
    /// The rule of TBLASTN, BLASTX, BLASTP, PSI-BLAST, RPS-BLAST and RPSTBLASTN, with their PSSM
    /// searches: the finite value (line PS DECISIONS 28, Owner-delegated; Q-ETOS8 open).
    pub const PROTEIN: Self = Self::FiniteQuotient;
}

/// `BlastKarlinEtoS_simple`: the score whose e-value is `evalue` in the search space `searchsp`
/// for the block `kbp`, `BLAST_SCORE_MIN` for a block with a negative Lambda, K or H (below every
/// caller's starting cutoff of 1). `evalue` is raised to `K_SMALL_FLOAT` first. `conversion` is
/// the caller's rule for a value C does not define.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4041-4064
#[inline]
fn karlin_etos_simple(
    evalue: f64,
    kbp: &KarlinParams,
    searchsp: i64,
    conversion: CutoffConversion,
) -> i32 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return BLAST_SCORE_MIN;
    }
    // C's `MAX`, under which a NaN e-value would become `kSmallFloat` too; none arrives here (the
    // one caller, `raise_cutoff`, comes only for `e > 0`, as `BLAST_Cutoffs` does).
    let e = evalue.max(K_SMALL_FLOAT);
    let searchsp_f = searchsp as f64;
    let quotient = kbp.k * searchsp_f / e;
    match conversion {
        CutoffConversion::X86Int4 => {
            ncbi_int4_from_double((quotient.glibc_ln() / kbp.lambda).ceil())
        }
        CutoffConversion::FiniteQuotient => {
            // A quotient past the double (an e-value near `kSmallFloat` with a search space above
            // about 4e12, or the largest `-searchsp`; with TBLASTX's word cutoff of `1e-300`
            // ordinary genome pairs) is +inf in NCBI, whose x86-64 `(Int4)` of it gives `INT_MIN`
            // and a cutoff of 1 (TBLASTX: a word cutoff of 1, every ungapped HSP to the sum
            // statistics; NCBI 2.17.0 ran out of memory or time on MG1655 1.2 Mb x Sakai
            // 5.5 Mb). LOSAT does not reproduce NCBI's undefined behaviour (Owner 2026-10-10):
            // the cutoff that every HSP of e-value at most E reaches.
            let log_value = if quotient.is_finite() {
                quotient.glibc_ln()
            } else {
                kbp.k.glibc_ln() + searchsp_f.glibc_ln() - e.glibc_ln()
            };
            let score_before_ceil = log_value / kbp.lambda;
            let score = score_before_ceil.ceil();
            // NCBI computes the cutoff without diagnostics; LOSAT's debug lines (no NCBI counterpart).
            if crate::utils::xenv::debug_cutoffs_is_some() && searchsp > 40_000_000_000 {
                eprintln!("[DEBUG CUTOFF_CALC] eff_searchsp={}", searchsp);
                eprintln!("[DEBUG CUTOFF_CALC] evalue={}, e_clamped={}", evalue, e);
                eprintln!("[DEBUG CUTOFF_CALC] K={}, Lambda={}", kbp.k, kbp.lambda);
                eprintln!(
                    "[DEBUG CUTOFF_CALC] K*searchsp={:.6e}, K*searchsp/E={:.6e}",
                    kbp.k * searchsp_f,
                    quotient
                );
                eprintln!(
                    "[DEBUG CUTOFF_CALC] log(K*searchsp/E)={:.6e}, score_before_ceil={:.6e}",
                    log_value, score_before_ceil
                );
                eprintln!("[DEBUG CUTOFF_CALC] final_score={} (ceil)", score as i32);
            }
            score as i32
        }
    }
}

/// The gap decay divisor of one segment (`BLAST_GapDecayDivisor`, `gap_decay_divisor` with one
/// segment) where `BLAST_Cutoffs` applies it: with the gap decay (`dodecay`, here `Some(rate)`) and
/// a rate strictly between 0 and 1.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4113-4118,4137-4144
#[inline]
fn single_segment_divisor(gap_decay: Option<f64>) -> Option<f64> {
    gap_decay
        .filter(|&rate| rate > 0.0 && rate < 1.0)
        .map(|rate| gap_decay_divisor(rate, 1))
}

/// The score part of `BLAST_Cutoffs`: the cutoff `score` raised to the score of `evalue` when
/// that is higher, and whether NCBI then writes the e-value of the cutoff back into `E` (the
/// e-value given was not above 0, or the cutoff did not change). A block of -1s returns at once
/// (`score` kept, no write-back).
#[inline]
fn raise_cutoff(
    score: &mut i32,
    evalue: f64,
    kbp: &KarlinParams,
    searchsp: i64,
    divisor: Option<f64>,
    conversion: CutoffConversion,
) -> bool {
    if kbp.lambda == -1.0 || kbp.k == -1.0 || kbp.h == -1.0 {
        return false;
    }
    let mut es = 1;
    // `e > 0` is false for 0, a negative e-value and NaN (`-evalue -nan` passes NCBI's option
    // check, whose test is `<= 0`): the score from E stays 1.
    if evalue > 0.0 {
        let e = match divisor {
            Some(divisor) => evalue * divisor,
            None => evalue,
        };
        es = karlin_etos_simple(e, kbp, searchsp, conversion);
    }
    let changed = es > *score;
    if changed {
        *score = es;
    }
    evalue <= 0.0 || !changed
}

/// `BLAST_Cutoffs`: raises the cutoff `score` to the score of the e-value `evalue` in the search
/// space `searchsp` when that is higher (`gap_decay`: `None` for NCBI's `dodecay` false,
/// `Some(rate)` for true), and, when `evalue` was not above 0 or the cutoff did not change, sets
/// `evalue` to the e-value of the cutoff (divided by the gap decay divisor). A block of -1s
/// changes neither. `conversion` is the caller's program's rule (`CutoffConversion`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4090-4150
fn blast_cutoffs(
    score: &mut i32,
    evalue: &mut f64,
    kbp: &KarlinParams,
    searchsp: i64,
    gap_decay: Option<f64>,
    conversion: CutoffConversion,
) {
    let divisor = single_segment_divisor(gap_decay);
    if raise_cutoff(score, *evalue, kbp, searchsp, divisor, conversion) {
        let back = evalue_from_raw_score(*score, kbp, searchsp as f64);
        *evalue = match divisor {
            Some(divisor) => back / divisor,
            None => back,
        };
    }
}

/// The cutoff of a caller whose `new_cutoff` starts at 1 and who does not read `E` back:
/// `BLAST_Cutoffs` for the hit-saving cutoff (blast_parameters.c:916 and 943, no gap decay) and
/// the ungapped word cutoff (blast_parameters.c:321 and 360, with the gap decay). At least 1, so
/// a large e-value in a small search space, an e-value not above 0 or an invalid block gives 1.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:321,360-363
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:916,943
#[inline]
pub fn cutoff_score(
    evalue: f64,
    kbp: &KarlinParams,
    searchsp: i64,
    gap_decay: Option<f64>,
    conversion: CutoffConversion,
) -> i32 {
    let mut score = 1;
    raise_cutoff(
        &mut score,
        evalue,
        kbp,
        searchsp,
        single_segment_divisor(gap_decay),
        conversion,
    );
    score
}

/// The sum statistics' cutoffs of a gapped search's contexts (`BlastHitSavingParametersUpdate`
/// with link parameters): `BLAST_Cutoffs` with the gap decay `gap_decay_rate`, from a cutoff of 1,
/// in the search space of the shorter of the average query length and the average subject length
/// times the average subject length. NCBI declares the e-value (`evalue_hsp`, 1) once before the
/// loop over the contexts, so the e-value that one context's `BLAST_Cutoffs` writes back (its
/// cutoff stayed 1: an average subject length of about 5 residues or less) is the next valid
/// context's; the caller asks `next_valid_context` once per valid context, in order (an invalid
/// context is skipped, as in NCBI).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:949-976
#[derive(Clone, Copy, Debug)]
pub struct SumStatsCutoffs {
    evalue_hsp: f64,
    searchsp: i64,
    gap_decay_rate: f64,
    conversion: CutoffConversion,
}

impl SumStatsCutoffs {
    /// The cutoffs for the average context length `avg_qlen` (NCBI's `Int4` quotient of the
    /// concatenated length by the number of contexts) and the average subject length
    /// `avg_subject_length`.
    pub fn new(
        avg_qlen: i32,
        avg_subject_length: i32,
        gap_decay_rate: f64,
        conversion: CutoffConversion,
    ) -> Self {
        Self {
            evalue_hsp: 1.0,
            searchsp: i64::from(avg_qlen.min(avg_subject_length)) * i64::from(avg_subject_length),
            gap_decay_rate,
            conversion,
        }
    }

    /// The sum statistics' cutoff of the next valid context, whose block is `kbp`. For contexts
    /// that share one block, the state after the first is the state after any number of them:
    /// the e-value is written back only while the cutoff stays 1, and then it is the e-value of a
    /// score of 1 each time, so a caller that resumes after valid contexts of that block (a
    /// BLASTX unit within its reader batch) steps once for them.
    pub fn next_valid_context(&mut self, kbp: &KarlinParams) -> i32 {
        let mut score = 1;
        blast_cutoffs(
            &mut score,
            &mut self.evalue_hsp,
            kbp,
            self.searchsp,
            Some(self.gap_decay_rate),
            self.conversion,
        );
        score
    }
}

#[cfg(test)]
#[path = "cutoffs_tests.rs"]
mod tests;
