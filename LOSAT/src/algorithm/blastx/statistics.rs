//! BLASTX HSP statistics; protein-subject lengths remain in amino acids.
use super::{
    args::ResolvedOptions,
    linking::{link_even, link_parameters, link_uneven, LinkedHsp, LinkedHspList},
    parameters::ContextParameters,
    preliminary::PreliminaryHsp,
    query_setup::PreparedQueryBatch,
};
use crate::core::link_hsps::{EvenGapLinkParams, LinkCutoffs};
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::{
    config::{ProteinScoringSpec, ScoringMatrix},
    stats::{karlin::evalue_from_raw_score, lookup_protein_params, spouge::blast_spouge_stoe},
};
use anyhow::{ensure, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1811-1839
pub fn evaluate_gapped(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject_length: i32,
    db_length: i64,
    scale: f64,
) -> Result<LinkedHspList> {
    ensure!(
        options.gapped,
        "BLASTX's gapped statistics called for an ungapped search"
    );
    // Speed (graduated 2026-10-27; the switch was `LOSAT_X_BXEVALEMPTY`): a subject without HSPs
    // gets at once the list that both ways below give for no input (`link_uneven`'s first lines,
    // and the loop below: no HSP, best e-value 0), without their tables; the sink evaluates every
    // subject, and most subjects of a database have no HSP. No NCBI counterpart: it does not change
    // any value NCBI computes.
    if input.is_empty() {
        return Ok(LinkedHspList {
            hsps: Vec::new(),
            best_evalue: 0.0,
        });
    }
    let spec = ProteinScoringSpec {
        matrix: options.scoring,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    };
    let mut ka = lookup_protein_params(&spec);
    ka.lambda /= scale;
    let gumbel = super::parameters::gumbel_block(options, db_length);
    let query_lengths: Vec<_> = batch.contexts.iter().map(|c| c.length as i32).collect();
    if let Some(link) = link_parameters(options) {
        return link_uneven(
            input,
            &query_lengths,
            parameters,
            subject_length,
            &vec![ka; parameters.len()],
            gumbel.as_ref(),
            &link,
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1875-1906
    let hsps: Vec<_> = input
        .iter()
        .enumerate()
        .map(|(source_index, h)| LinkedHsp {
            hsp: h.clone(),
            num: 0,
            evalue: match &gumbel {
                Some(gumbel) => blast_spouge_stoe(
                    h.score,
                    &ka,
                    gumbel,
                    query_lengths[h.context],
                    subject_length,
                ),
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1886-1891
                None => {
                    evalue_from_raw_score(h.score, &ka, parameters[h.context].search_space as f64)
                }
            },
            source_index,
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1737-1749
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1900-1904
    let best_evalue = if hsps.is_empty() {
        0.0
    } else {
        hsps.iter()
            .fold(i32::MAX as f64, |best, h| h.evalue.min(best))
    };
    Ok(LinkedHspList { hsps, best_evalue })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:868,931-941
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1015-1017
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443
/// NCBI's `params->prelim_evalue`, the e-value of the preliminary reap (`reap_preliminary`):
/// `-evalue`, stretched 5 times with the Gumbel block (`ResolvedOptions::has_gumbel_block`) and a
/// composition-based statistics above 1, except where NCBI computes the statistics of each subject
/// (`per_subject`: `-subject`, a database length of 0), whose update of the hit-saving parameters
/// passes a composition of 0 before the subject's search (line RG item X-evalue-nan; LOSAT
/// stretched it there too).
pub fn prelim_evalue(options: &ResolvedOptions, per_subject: bool) -> f64 {
    if options.composition > 1 && options.has_gumbel_block() && !per_subject {
        5.0 * options.evalue
    } else {
        options.evalue
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:643-669
/// NCBI's preliminary reap of the HSPs whose e-value is greater than `cutoff`
/// (`prelim_evalue`); an HSP of NaN, or of any e-value with a cutoff of NaN, stays.
pub fn reap_preliminary(list: &mut LinkedHspList, cutoff: f64) -> Vec<LinkedHsp> {
    let mut deleted = Vec::new();
    list.hsps.retain(|h| {
        if h.evalue > cutoff {
            deleted.push(h.clone());
            false
        } else {
            true
        }
    });
    deleted
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-2008
pub fn reap_final(list: &mut LinkedHspList, options: &ResolvedOptions) -> Vec<LinkedHsp> {
    let mut deleted = Vec::new();
    list.hsps.retain(|h| {
        if h.evalue > options.evalue {
            deleted.push(h.clone());
            false
        } else {
            true
        }
    });
    deleted
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:102-114
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1917-1934
pub fn normalize_scores(list: &mut LinkedHspList, options: &ResolvedOptions) -> Vec<f64> {
    let ka = lookup_protein_params(&ProteinScoringSpec {
        matrix: options.scoring,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    let scale = super::kappa::redo_scale(options);
    list.hsps
        .iter_mut()
        .map(|h| {
            if scale != 1.0 {
                let value = h.hsp.score as f64 / scale;
                h.hsp.score = (value + if value >= 0.0 { 0.5 } else { -0.5 }) as i32;
            }
            (h.hsp.score as f64 * (ka.lambda / scale) * scale - ka.k.glibc_ln())
                / std::f64::consts::LN_2
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:511-532
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:767-833
pub fn identity_and_positive(
    h: &super::kappa::RedoneHsp,
    batch: &PreparedQueryBatch,
    subject: &[u8],
    matrix: ScoringMatrix,
) -> (usize, usize) {
    let p = &h.hsp;
    let c = &batch.contexts[p.context];
    let query = &batch.sequence_start_nomask[1 + c.offset..1 + c.offset + c.length];
    let mut qi = p.q_start as usize;
    let mut si = p.s_start as usize;
    let mut identity = 0;
    let mut positive = 0;
    for op in &h.edit_script {
        match *op {
            crate::common::GapEditOp::Sub(n) => {
                for _ in 0..n {
                    let q = query[qi];
                    let s = subject[si];
                    if q == s {
                        identity += 1;
                        positive += 1;
                    } else if h.identity_matrix.as_ref().map_or_else(
                        || crate::utils::matrix::protein_score(matrix, q, s) > 0,
                        |matrix| matrix.score(q, s) > 0,
                    ) {
                        positive += 1;
                    }
                    qi += 1;
                    si += 1;
                }
            }
            crate::common::GapEditOp::Del(n) => si += n as usize,
            crate::common::GapEditOp::Ins(n) => qi += n as usize,
        }
    }
    (identity, positive)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
pub(crate) fn evaluate_ungapped(
    input: &[PreliminaryHsp],
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject_length: i32,
    prepared: Option<LinkCutoffs>,
) -> Result<LinkedHspList> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1768-1769
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1830-1831
    // Both `BLAST_LinkHsps` and `Blast_HSPListGetEvalues` return immediately when the HSP list is
    // empty (keeping `best_evalue` at 0.0), so skipped subjects with no per-subject link cutoffs
    // do not reach the `prepared` unwrap below.
    if input.is_empty() {
        return Ok(LinkedHspList {
            hsps: Vec::new(),
            best_evalue: 0.0,
        });
    }
    if let Some(link) = link_parameters(options) {
        if link.longest_intron <= 0 {
            return Ok(link_even(
                input,
                batch,
                parameters,
                subject_length,
                &EvenGapLinkParams::ungapped(prepared.expect("per-subject even-gap cutoffs"), true),
            ));
        }
        let lengths: Vec<_> = batch.contexts.iter().map(|c| c.length as i32).collect();
        let ka: Vec<_> = parameters.iter().map(|p| p.ungapped).collect();
        return link_uneven(
            input,
            &lengths,
            parameters,
            subject_length,
            &ka,
            None,
            &link,
        );
    }
    let hsps: Vec<_> = input
        .iter()
        .enumerate()
        .map(|(source_index, h)| {
            let p = &parameters[h.context];
            LinkedHsp {
                hsp: h.clone(),
                num: 0,
                evalue: evalue_from_raw_score(h.score, &p.ungapped, p.search_space as f64),
                source_index,
            }
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1900-1904
    let best_evalue = if hsps.is_empty() {
        0.0
    } else {
        hsps.iter()
            .fold(i32::MAX as f64, |best, h| h.evalue.min(best))
    };
    Ok(LinkedHspList { hsps, best_evalue })
}
