//! Round 1 of a PHI search for one query, per subject: the setup (occurrences, statistics,
//! cutoff, X-dropoffs), the preliminary HSPs of every occurrence against every subject hit, the
//! traceback with e-values, reaping and bit scores. The driver (subjects, hit lists, the split
//! per occurrence, report) is blastp's (`blastp::blast_engine::phi_hooks`).

use super::align::{prelim_alignment, traceback, PhiGapCosts};
use super::pattern::{PhiError, PhiPattern};
use super::scan::{query_occurrences, subject_hits, Occurrence, Occurrences};
use super::stats::{self, PhiKarlin};
use crate::algorithm::pssm_search::gapalign::PssmGapScratch;
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::common::GapEditOp;
use crate::core::blast_util::BLASTAA_SIZE;

/// An HSP of a PHI search (`BlastHSP` with its `pat_info`).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PhiHsp {
    /// The query occurrence (`pat_info->index`).
    pub(crate) occurrence: usize,
    /// The subject hit's length (`pat_info->length`).
    pub(crate) s_pattern_length: i32,
    /// The two pattern starts (`query.gapped_start`, `subject.gapped_start`).
    pub(crate) q_pattern_offset: i32,
    pub(crate) s_pattern_offset: i32,
    pub(crate) query_start: i32,
    pub(crate) query_end: i32,
    pub(crate) subject_start: i32,
    pub(crate) subject_end: i32,
    pub(crate) score: i32,
    /// The e-value: 0 in the preliminary stage (NCBI computes none for PHI-BLAST) unless
    /// psiblast's `-sum_stats true` links the subject's HSPs (`PhiSearch::link_preliminary`),
    /// then the pattern's from the traceback (`PhiSearch::traceback`).
    pub(crate) evalue: f64,
    pub(crate) bit_score: f64,
    /// The size of the HSP's linked set (`num`): 1, or the size that the preliminary linking of
    /// `-sum_stats true` gives, which no later stage changes (the Seq-align's `sum_n`).
    pub(crate) num: i32,
    pub(crate) edit_script: Vec<GapEditOp>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1357
/// `ScoreCompareHSPs`: score down, then subject start up, subject end down, query start up, query
/// end down.
pub(crate) fn score_compare(a: &PhiHsp, b: &PhiHsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.subject_start.cmp(&b.subject_start))
        .then(b.subject_end.cmp(&a.subject_end))
        .then(a.query_start.cmp(&b.query_start))
        .then(b.query_end.cmp(&a.query_end))
}

/// Why a query is not searched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum PhiSetupError {
    /// `The pattern %s was not found in the query.` (the pattern's text as bytes).
    NotFound(Vec<u8>),
    /// `The pattern (%s) may not cover the entire query.`
    WholeQuery(Vec<u8>),
    /// LOSAT's refusal from the query's scan.
    Refusal(String),
}

impl PhiSetupError {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1065-1118
    /// NCBI's warning text, which its engine raises as an error.
    pub(crate) fn message(&self) -> Vec<u8> {
        match self {
            Self::NotFound(text) => [
                b"The pattern ".as_slice(),
                text,
                b" was not found in the query.",
            ]
            .concat(),
            Self::WholeQuery(text) => [
                b"The pattern (".as_slice(),
                text,
                b") may not cover the entire query.",
            ]
            .concat(),
            Self::Refusal(message) => message.clone().into_bytes(),
        }
    }
}

/// psiblast's `-sum_stats true` in a PHI search: the linking of each subject's preliminary HSPs
/// (`BLAST_LinkHsps`: the preliminary stage's parameters, whose small-gap cutoff is
/// `stats::word_cutoff`, and the context of `stats::link_context`) and the preliminary e-value
/// that reaps them (`-evalue`, `hit_params->prelim_evalue`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PhiLink {
    pub(crate) params: crate::core::link_hsps::EvenGapLinkParams,
    pub(crate) context: crate::core::link_hsps::LinkContext,
    pub(crate) prelim_evalue: f64,
}

/// The search of one query: its letters, matrix rows, occurrences and the stage parameters.
pub(crate) struct PhiSearch<'p> {
    pub(crate) pattern: &'p PhiPattern,
    pub(crate) query: Vec<u8>,
    pub(crate) rows: QueryPssm,
    pub(crate) occurrences: Vec<Occurrence>,
    pub(crate) karlin: PhiKarlin,
    pub(crate) effective_occurrences: i32,
    pub(crate) search_space: i64,
    pub(crate) cutoff: i32,
    pub(crate) prelim: PhiGapCosts,
    pub(crate) last: PhiGapCosts,
    /// psiblast's `-sum_stats true` (`PhiLink`), `None` without sum statistics; the caller sets
    /// it (blastp's `phi_hooks::phi_start`: the context reads the query batch's first query).
    pub(crate) link: Option<PhiLink>,
}

/// The X-dropoffs in raw score (`BlastExtensionParametersNew` with the table's Lambda).
pub(crate) fn x_dropoffs(lambda: f64, x_drop_bits: f64, x_drop_final_bits: f64) -> (i32, i32) {
    let ln2 = std::f64::consts::LN_2;
    let prelim = (x_drop_bits * ln2 / lambda) as i32;
    let last = (x_drop_final_bits * ln2 / lambda).max(f64::from(prelim)) as i32;
    (prelim, last)
}

/// The parameters of the setup: the pattern (`'p`, which the search keeps) and the query's
/// values (`'a`, which it copies).
pub(crate) struct PhiSetup<'p, 'a> {
    pub(crate) pattern: &'p PhiPattern,
    /// The query's letters (NCBIstdaa, masked as the search sees them).
    pub(crate) query: &'a [u8],
    /// The lookup segments (inclusive ranges of `query`).
    pub(crate) segments: &'a [(usize, usize)],
    /// The score of a query letter against a subject letter.
    pub(crate) score: &'a dyn Fn(u8, u8) -> i32,
    pub(crate) karlin: PhiKarlin,
    /// The database's length and sequence count (after `-dbsize`).
    pub(crate) db_length: i64,
    pub(crate) db_num_seqs: i64,
    pub(crate) evalue: f64,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    pub(crate) x_drop_bits: f64,
    pub(crate) x_drop_final_bits: f64,
}

/// The occurrences of `pattern` in the query's letters `query` within its lookup segments
/// `segments` (`Blast_SetPHIPatternInfo`), or NCBI's error of a query without one or covered by
/// one, or LOSAT's refusal of the scan. `setup` starts with it; a search whose matrix and gap
/// costs have no entry in NCBI's table calls it alone, for its errors (`blastp`'s `phi_start`).
pub(crate) fn occurrences(
    pattern: &PhiPattern,
    query: &[u8],
    segments: &[(usize, usize)],
) -> Result<Vec<Occurrence>, PhiSetupError> {
    match query_occurrences(pattern, query, segments) {
        Ok(Occurrences::Found(found)) if found.is_empty() => {
            Err(PhiSetupError::NotFound(pattern.text.clone()))
        }
        Ok(Occurrences::Found(found)) => Ok(found),
        Ok(Occurrences::WholeQuery) => Err(PhiSetupError::WholeQuery(pattern.text.clone())),
        Err(PhiError::Refusal(message)) | Err(PhiError::Ncbi(message)) => {
            Err(PhiSetupError::Refusal(message))
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:893-907
/// The query's setup: occurrences (`Blast_SetPHIPatternInfo`), search space, effective count,
/// cutoff (`5 * -evalue`), X-dropoffs and the matrix rows of the query's letters.
pub(crate) fn setup<'p>(params: &PhiSetup<'p, '_>) -> Result<PhiSearch<'p>, PhiSetupError> {
    let pattern = params.pattern;
    let occurrences = occurrences(pattern, params.query, params.segments)?;
    let search_space = stats::search_space(
        params.db_length,
        params.db_num_seqs,
        pattern.min_match_length,
    );
    let effective = stats::effective_occurrences(&occurrences, pattern.min_match_length);
    let cutoff = stats::cutoff_score(
        5.0 * params.evalue,
        search_space,
        &params.karlin,
        effective,
        pattern.probability,
    );
    let (x_prelim, x_last) = x_dropoffs(
        params.karlin.lambda,
        params.x_drop_bits,
        params.x_drop_final_bits,
    );
    let rows = QueryPssm {
        rows: params
            .query
            .iter()
            .map(|&q| {
                let mut row = [0i32; BLASTAA_SIZE];
                for (s, cell) in row.iter_mut().enumerate() {
                    *cell = (params.score)(q, s as u8);
                }
                row
            })
            .collect(),
    };
    let costs = |x_drop| PhiGapCosts {
        gap_open: params.gap_open,
        gap_extend: params.gap_extend,
        x_drop,
    };
    Ok(PhiSearch {
        pattern,
        query: params.query.to_vec(),
        rows,
        occurrences,
        karlin: params.karlin,
        effective_occurrences: effective,
        search_space,
        cutoff,
        prelim: costs(x_prelim),
        last: costs(x_last),
        link: None,
    })
}

impl PhiSearch<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:739-835
    /// `PHIGetGappedScore` of one subject (`letters`, then its sentinel at `letters.len()` in
    /// `subject`): every occurrence against every subject hit, kept from the cutoff, sorted by
    /// score; with the subject's hit count (`lookup_hits`).
    pub(crate) fn prelim(
        &self,
        subject: &[u8],
        subject_length: usize,
        scratch: &mut PssmGapScratch,
    ) -> Result<(Vec<PhiHsp>, usize), PhiError> {
        let hits = subject_hits(self.pattern, &subject[..subject_length])?;
        let mut hsps = Vec::new();
        for (index, occurrence) in self.occurrences.iter().enumerate() {
            for hit in &hits {
                let s_length = hit.end - hit.start + 1;
                let aligned = prelim_alignment(
                    &self.rows,
                    subject,
                    subject_length,
                    occurrence.offset as usize,
                    hit.start as usize,
                    occurrence.length as usize,
                    s_length as usize,
                    self.prelim,
                    scratch,
                );
                if aligned.score >= self.cutoff {
                    hsps.push(PhiHsp {
                        occurrence: index,
                        s_pattern_length: s_length,
                        q_pattern_offset: occurrence.offset,
                        s_pattern_offset: hit.start,
                        query_start: aligned.query_start,
                        query_end: aligned.query_stop,
                        subject_start: aligned.subject_start,
                        subject_end: aligned.subject_stop,
                        score: aligned.score,
                        evalue: 0.0,
                        bit_score: 0.0,
                        num: 1,
                        edit_script: Vec::new(),
                    });
                }
            }
        }
        hsps.sort_by(score_compare);
        Ok((hsps, hits.len()))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-900
    /// psiblast's `-sum_stats true` on a subject's preliminary HSPs `hsps` (`link`, sorted by
    /// score), after NCBI's search of the subject's parts: `BLAST_LinkHsps` with the PHI block,
    /// which gives each HSP its linked set's e-value and size and sorts the list by score again,
    /// then the reap of the HSPs past the preliminary e-value (`s_Blast_HSPListReapByPrelimEvalue`,
    /// the order kept). `subject_length` is the whole subject's. Without sum statistics the HSPs
    /// keep the preliminary e-value 0, which the reap keeps.
    pub(crate) fn link_preliminary(&self, hsps: &mut Vec<PhiHsp>, subject_length: usize) {
        let Some(link) = &self.link else {
            return;
        };
        crate::core::link_hsps::link_in_place(
            hsps,
            // `BlastHSP`'s 0-based offsets and exclusive ends, all in context 0.
            |hsp| crate::core::link_hsps::LinkHsp {
                context: 0,
                subject_frame: 0,
                query_offset: hsp.query_start,
                query_end: hsp.query_end,
                subject_offset: hsp.subject_start,
                subject_end: hsp.subject_end,
                score: hsp.score,
            },
            |_| link.context,
            subject_length as i32,
            &link.params,
            |hsp, linked| {
                hsp.evalue = linked.evalue;
                hsp.num = linked.num;
            },
        );
        hsps.sort_by(score_compare);
        // `hsp->evalue > cutoff` frees the HSP.
        hsps.retain(|hsp| {
            hsp.evalue.partial_cmp(&link.prelim_evalue) != Some(std::cmp::Ordering::Greater)
        });
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:740-815
    /// `s_PHITracebackFromHSPList`: each HSP's traceback (freed under the cutoff), sorted by
    /// score, e-values with the database's pattern count, reaped by `-evalue`, bit scores.
    pub(crate) fn traceback(
        &self,
        subject: &[u8],
        subject_length: usize,
        hsps: Vec<PhiHsp>,
        num_patterns_db: i32,
        evalue_threshold: f64,
        scratch: &mut PssmGapScratch,
    ) -> Result<Vec<PhiHsp>, PhiError> {
        let mut kept = Vec::with_capacity(hsps.len());
        for mut hsp in hsps {
            let occurrence = self.occurrences[hsp.occurrence];
            let aligned = traceback(
                &self.rows,
                &self.query,
                subject,
                subject_length,
                hsp.q_pattern_offset as usize,
                hsp.s_pattern_offset as usize,
                occurrence.length as usize,
                hsp.s_pattern_length as usize,
                self.pattern,
                self.last,
                scratch,
            )?;
            if aligned.score >= self.cutoff {
                hsp.score = aligned.score;
                hsp.query_start = aligned.query_start;
                hsp.query_end = aligned.query_stop;
                hsp.subject_start = aligned.subject_start;
                hsp.subject_end = aligned.subject_stop;
                hsp.edit_script = aligned.edit_script;
                kept.push(hsp);
            }
        }
        kept.sort_by(score_compare);
        for hsp in &mut kept {
            hsp.evalue = stats::evalue(
                hsp.score,
                &self.karlin,
                self.effective_occurrences,
                num_patterns_db,
            );
        }
        kept.retain(|hsp| hsp.evalue <= evalue_threshold);
        for hsp in &mut kept {
            hsp.bit_score = stats::bit_score(hsp.score, &self.karlin);
        }
        Ok(kept)
    }
}

#[cfg(test)]
#[path = "engine_tests.rs"]
mod tests;
