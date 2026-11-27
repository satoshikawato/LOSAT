//! The traceback of a PSSM search with composition-based statistics (`-comp_based_stats` 1, and
//! 2 and 3, which NCBI runs as 1 with a PSSM): NCBI's `Blast_RedoAlignmentCore_MT` for one
//! position-based query (series PS, item PS4b). Every preliminary list is redone on the PSSM
//! that NCBI's `Blast_CompositionBasedStats` scales for the subject's composition
//! (`core::composition_adjustment::compo_scale`), with the gapped alignment of `gapalign.rs` on
//! the scaled rows. The e-values and bit scores read the gapped Lambda of the PSSM file over the
//! scaling factor and the K of the factor-1 PSSM (`PosMatrixInfo::gapped_k`); the report's
//! footer keeps the file's values (`s_RestoreSearch`).
//!
//! The matches are redone in the order of the preliminary hit list, one after the other, at
//! every thread count: NCBI's one-thread loop, whose early termination test and heap insertions
//! read the matches before (DECISIONS 13 of line PS). The PSSM that the last adjustment left
//! stays for the next match, as NCBI's thread-local copy does.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3058-3066

use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::Result;

use super::evalue::EvalueStats;
use super::gapalign::{pssm_gapped_alignment_with_traceback, PssmGapScratch};
use super::rows::QueryPssm;
use super::traceback0::{final_hsp, TraceReads};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{
    fill_results_from_compo_heaps, BlastCompoHeap, BlastpHitList, BlastpHsp, BlastpHspList,
    COMPO_HEAP_EVALUE_STRETCH,
};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::blastp::x_redo_subjects::XRedoSubjects;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::common::GapEditOp;
use crate::core::composition_adjustment::adjust_scores::{
    blast_get_composition_range, read_aa_composition, BlastAminoAcidComposition,
};
use crate::core::composition_adjustment::compo_scale::{
    pssm_composition_based_stats, FrameRedo, PosMatrixInfo, RedoInstance, XDropRedo,
};
use crate::core::composition_adjustment::redo_alignment::{
    build_query_word_hashes, get_nucl_length, get_seq_frame, get_translated_length,
    test_near_identical, test_near_identical_on_copy, windows_from_translated_aligns,
    BlastCompoAlignment, BlastCompoAlignmentContext, BlastCompoSequenceData, EMatrixAdjustRule,
    NearIdenticalCopy, MINIMUM_LENGTH_NEAR_IDENTICAL, WINDOW_BORDER,
};
use crate::core::link_hsps::{link_in_place, EvenGapLinkParams, LinkContext, LinkHsp};
use crate::stats::KarlinParams;
use crate::utils::matrix::ncbistdaa;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:676
/// `SCALING_FACTOR`: the redo's scores, Lambda and gap costs are scaled by it.
pub(crate) const SCALING_FACTOR: f64 = 32.0;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2396-2399
const NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:46,63-64
// (LOCAL_LN2 and EVALUE_STRETCH)
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:955
const KAPPA_BIT_TOL: f64 = 2.0;
const LOCAL_LN2: f64 = std::f64::consts::LN_2; // the same f64 as NCBI's literal

/// What the redo of every match reads: the PSSM, its matrix information, the query, the
/// subjects, the scoring and hit-saving options.
pub(crate) struct PssmRedo1<'a> {
    /// The search's PSSM (the scores of the file): the matrix of NCBI's copy of the score block
    /// until the first adjustment overwrites it.
    pub(crate) pssm: &'a QueryPssm,
    /// The position-based `Blast_MatrixInfo` (`PosMatrixInfo::new`).
    pub(crate) info: &'a PosMatrixInfo,
    /// The query's context (its letters, its length).
    pub(crate) context: &'a QueryContext,
    pub(crate) subjects: XRedoSubjects<'a>,
    /// The gap costs of the search, before the scaling.
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    /// The matrix of the search (`sbp->matrix`), for the positives.
    pub(crate) matrix: crate::config::ScoringMatrix,
    /// `kbp_gap_psi[0]`, the PSSM file's gapped block (`sbp->kbp_gap` before the scaling).
    pub(crate) karlin: &'a KarlinParams,
    /// `extendParams->gap_x_dropoff_final`, the final gapped X-drop of the search (raw score).
    pub(crate) x_drop_final: i32,
    /// `options->gap_x_dropoff_final`, the final gapped X-drop in bits (`-xdrop_gap_final` or the
    /// protein default), which the redo converts at its own Lambda (`RedoRun::new`).
    pub(crate) x_drop_final_bits: f64,
    /// `sbp->gbp` (the redo's copies of the score block have it when the search has it:
    /// blast_stat.c:921-923 and blast_kappa.c:2865-2867) and the context's effective search
    /// space.
    pub(crate) stats: EvalueStats<'a>,
    pub(crate) expect_value: f64,
    /// `hitParams->options->hitlist_size` (psiblast's preliminary size).
    pub(crate) hitlist_size: usize,
    /// `psiOptions->inclusion_ethresh`, the heaps' e-value cutoff.
    pub(crate) inclusion_ethresh: f64,
    /// psiblast's `-sum_stats true` (`do_link_hsps`), `None` without sum statistics.
    pub(crate) link: Option<RedoLink>,
    /// Called with the heaps as the redo leaves them (the queries' batch is all of them), before
    /// they become the final hit lists: the search notes there which subjects the redo's fetches
    /// read (`blastp::blast_engine::database::ListedSubjects::note_fetched`), and refuses the
    /// search where one of them has a header that cannot be read.
    pub(crate) redone: &'a dyn Fn(&[BlastCompoHeap]) -> anyhow::Result<()>,
    /// psiblast's `-use_sw_tback` (`Blast_RedoOneMatchSmithWaterman`, CBS 0 to 3); `None` for the
    /// default X-drop redo.
    pub(crate) sw: Option<PssmSwOptions<'a>>,
}

/// Options for psiblast's Smith-Waterman position-based redo (`-use_sw_tback`, CBS 0 to 3).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039,3100-3109,3602-3622
pub(crate) struct PssmSwOptions<'a> {
    /// Whether composition-based adjustment runs (`true` for CBS 1..3, `false` for CBS 0).
    pub(crate) compo_adjust: bool,
    /// Applies fetch-excluded marks of the finished query to the batch status.
    pub(crate) on_query_done: &'a dyn Fn(
        Option<usize>,
        &[BlastCompoHeap],
        &mut crate::core::composition_adjustment::smith_waterman_match::SwBatchStatus,
    ),
    /// Builds the error/refusal when the batch's last redone match failed.
    pub(crate) on_batch_failed: &'a dyn Fn(&[BlastCompoHeap]) -> anyhow::Error,
    /// Builds the refusal when a match's SW redo does not terminate in NCBI.
    pub(crate) on_endless: &'a dyn Fn(&[BlastCompoHeap], usize) -> anyhow::Error,
}

/// psiblast's `-sum_stats true` in the redo (`do_link_hsps`): the redone alignments' cutoff and
/// the linking of `s_HitlistEvaluateAndPurge` in the place of the e-values.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2450-2461
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:410-415
#[derive(Debug, Clone, Copy)]
pub(crate) struct RedoLink {
    /// `params->cutoff_s`: the traceback's `cutoff_score_min` times the scaling factor,
    /// truncated.
    pub(crate) cutoff_s: i32,
    /// The traceback's linking parameters (`EvenGapLinkParams::traceback`).
    pub(crate) params: EvenGapLinkParams,
    /// The query's context; the linking reads the Lambda and logK of the redo's scaled block
    /// (`sbp->kbp_gap` of NCBI's copy, `RedoRun`) in the place of this one's.
    pub(crate) context: LinkContext,
}

impl RedoLink {
    /// The parameters of the redo of the search whose traceback cutoff is `cutoff_score_min`
    /// (one context) and whose query context is `context`, scaled by `SCALING_FACTOR`.
    #[cfg(test)]
    pub(crate) fn new(cutoff_score_min: i32, context: LinkContext) -> Self {
        Self::with_scale(cutoff_score_min, context, SCALING_FACTOR)
    }

    /// The parameters of the redo scaled by `scaling_factor` (`SCALING_FACTOR` for CBS 1..3,
    /// `1.0` for CBS 0 with `-use_sw_tback`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2450-2461,3100-3109
    pub(crate) fn with_scale(
        cutoff_score_min: i32,
        context: LinkContext,
        scaling_factor: f64,
    ) -> Self {
        Self {
            cutoff_s: (cutoff_score_min as f64 * scaling_factor) as i32,
            params: EvenGapLinkParams::traceback(),
            context,
        }
    }
}

/// The values of the search that the redo scales (`s_RescaleSearch`) and derives
/// (`s_GetAlignParams`), and its working memory.
pub(crate) struct RedoRun<'a> {
    /// The query (NCBIstdaa, O as X) and its unmasked letters.
    query: &'a [u8],
    query_nomask: &'a [u8],
    /// The query range's letters with U as C (`s_SequenceGetRange`), for the near-identical test.
    query_data: BlastCompoSequenceData,
    /// `query_info->words`.
    query_words: Vec<u64>,
    /// `query_composition->numTrueAminoAcids`.
    query_true_amino_acids: i32,
    /// The scaled gapped block of the redo's copy: Lambda over the scaling factor, K of the
    /// factor-1 PSSM.
    karlin: KarlinParams,
    log_k: f64,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    near_identical_cutoff: f64,
    /// The PSSM of NCBI's thread-local copy: what the last adjustment left.
    matrix: QueryPssm,
    scratch: PssmGapScratch,
    /// RPS-BLAST: whether each adjusted matrix is checked against `INT4_MAX` (the sum of its
    /// rows' best scores), and the first sum that was above it.
    bound_check: bool,
    bound_exceeded: Option<i64>,
    /// RPS-BLAST: what the redo asks of its profile where NCBI's reads leave the data
    /// (`with_profile_reads`).
    profile_reads: Option<&'a dyn TraceReads>,
}

/// `BlastCompo_Alignment`: an incoming (preliminary) alignment, scaled, or a redone one with its
/// edit script.
#[derive(Debug, Clone)]
struct Align {
    score: i32,
    rule: EMatrixAdjustRule,
    query_start: i32,
    query_end: i32,
    match_start: i32,
    match_end: i32,
    frame: i32,
    edit_script: Vec<GapEditOp>,
}

impl Align {
    /// The incoming alignment of a preliminary HSP (`s_ResultHspToDistinctAlign`): its score
    /// scaled, its subject frame.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:788-800
    /// `C`, the redo's instance, keeps each instance's own copy of the code that collects them
    /// (psiblast's code as before the frames' instance; DECISIONS 63 S3).
    fn incoming<C: RedoInstance>(hsp: &BlastpPreliminaryHsp) -> Self {
        Align {
            score: (hsp.raw_score as f64 * SCALING_FACTOR) as i32,
            rule: EMatrixAdjustRule::DontAdjustMatrix,
            query_start: hsp.query_start,
            query_end: hsp.query_end,
            match_start: hsp.subject_start,
            match_end: hsp.subject_end,
            frame: hsp.subject_frame,
            edit_script: Vec::new(),
        }
    }

    /// The alignment that a redo does not make (`RedoRun::starts_in_rows`): score 0, below every
    /// cutoff of the redo, so no list keeps it. Out of the redo's loops.
    #[cold]
    #[inline(never)]
    fn not_made(rule: EMatrixAdjustRule) -> Self {
        Align {
            score: 0,
            rule,
            query_start: 0,
            query_end: 0,
            match_start: 0,
            match_end: 0,
            frame: 0,
            edit_script: Vec::new(),
        }
    }

    /// The alignment as the composition library holds it (`BlastCompo_Alignment` of query 0),
    /// with `context`.
    #[inline(always)]
    fn compo(&self, context: Option<BlastCompoAlignmentContext>) -> BlastCompoAlignment {
        BlastCompoAlignment {
            score: self.score,
            matrix_adjust_rule: self.rule,
            query_index: 0,
            query_start: self.query_start,
            query_end: self.query_end,
            match_start: self.match_start,
            match_end: self.match_end,
            frame: self.frame,
            context,
            next: None,
        }
    }
}

/// A redone HSP of a match's list (`BlastHSP`), before its identities are counted: its e-value and
/// the size of its linked set (`num`, 1 without sum statistics).
struct Redone {
    align: Align,
    e_value: f64,
    num: i32,
}

impl<'a> RedoRun<'a> {
    /// The values of the search that the redo scales (`s_RescaleSearch` on `kbp_gap` and the gap
    /// costs, `s_GetAlignParams`) and the query's data, from the scaled gapped block `karlin`
    /// (Lambda over the scaling factor, the K of the redo's copy) with its `log_k`, the search's
    /// unscaled gap costs, the final gapped X-drop of the search in bits (`-xdrop_gap_final` or the
    /// protein default) and in raw score, and the PSSM the redo starts from (the matrix of NCBI's
    /// thread-local copy, which each adjustment overwrites). psiblast's run starts from the file's
    /// PSSM; RPS-BLAST's from zeros, one run for each profile.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2133
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2437-2450
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2374-2383
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2331 (`s_GetQueryInfo`:
    /// `Blast_ReadAaComposition` and `s_CreateWordArray` of the query)
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        query: &'a [u8],
        query_nomask: &'a [u8],
        karlin: KarlinParams,
        log_k: f64,
        gap_costs: (i32, i32),
        (x_drop_final_bits, x_drop_final): (f64, i32),
        matrix: QueryPssm,
        scratch: PssmGapScratch,
    ) -> Self {
        Self::new_with_scale(
            query,
            query_nomask,
            karlin,
            log_k,
            gap_costs,
            (x_drop_final_bits, x_drop_final),
            matrix,
            scratch,
            SCALING_FACTOR,
        )
    }

    /// `new` with an explicit `scaling_factor` (`SCALING_FACTOR` for CBS 1..3, `1.0` for CBS 0
    /// with `-use_sw_tback`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2133,3100-3109
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new_with_scale(
        query: &'a [u8],
        query_nomask: &'a [u8],
        karlin: KarlinParams,
        log_k: f64,
        gap_costs: (i32, i32),
        (x_drop_final_bits, x_drop_final): (f64, i32),
        matrix: QueryPssm,
        scratch: PssmGapScratch,
        scaling_factor: f64,
    ) -> Self {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1682-1695
        let query_u_as_c: Vec<u8> = query
            .iter()
            .map(|&residue| {
                if residue == ncbistdaa::U {
                    ncbistdaa::C
                } else {
                    residue
                }
            })
            .collect();
        let lambda = karlin.lambda;
        RedoRun {
            query,
            query_nomask,
            query_data: BlastCompoSequenceData::from_ncbistdaa(&query_u_as_c),
            query_words: build_query_word_hashes(query),
            query_true_amino_acids: read_aa_composition(query).num_true_amino_acids,
            karlin,
            log_k,
            gap_open: blast_nint(gap_costs.0 as f64 * scaling_factor),
            gap_extend: blast_nint(gap_costs.1 as f64 * scaling_factor),
            x_dropoff: crate::core::x_drop::redo_x_dropoff(x_drop_final_bits, lambda, x_drop_final),
            near_identical_cutoff: (NEAR_IDENTICAL_BITS_PER_POSITION * std::f64::consts::LN_2)
                / lambda,
            matrix,
            scratch,
            bound_check: false,
            bound_exceeded: None,
            profile_reads: None,
        }
    }

    /// The larger of `gap_x_dropoff_final` times `NCBIMATH_LN2` over `min_lambda` and
    /// `extendParams->gap_x_dropoff_final`, before NCBI's conversion to `Int4`
    /// (blast_kappa.c:2381-2383, cited at `new`): `new`'s value with the final X-drop `bits` and
    /// the redo's Lambda `lambda`.
    pub(crate) fn x_dropoff_of(bits: f64, lambda: f64, x_drop_final: i32) -> f64 {
        (bits * std::f64::consts::LN_2 / lambda).max(x_drop_final as f64)
    }

    /// Checks every matrix the composition adjustment builds: an alignment uses each row once,
    /// so the sum of the rows' best scores bounds its score, and above `INT4_MAX` NCBI's 32-bit
    /// scores overflow (undefined behavior in C).
    pub(crate) fn with_bound_check(mut self) -> Self {
        self.bound_check = true;
        self
    }

    /// The run of a redo whose query information has no valid context: the cutoff of the
    /// near-identical test is 0, NCBI's old test (`old_test_near_identical`). RPS-BLAST's redo
    /// keeps context 0 of the user query's information only, with its validity: rpstblastn's
    /// query whose frame +1 is invalid (`-strand minus`, a frame +1 without letters) gets it.
    /// A protein query has one context, valid where it has hits.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1244
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2436-2446
    /// (`Blast_GetOneQueryStructs` copies the query's contexts with `is_valid`, and
    /// `OffsetArrayToContextOffsets` sets only the offsets, lengths, frames and query indices.)
    pub(crate) fn without_valid_context(mut self) -> Self {
        self.near_identical_cutoff = 0.0;
        self
    }

    /// The first sum of best row scores above `INT4_MAX` that an adjusted matrix had.
    pub(crate) fn bound_exceeded(&self) -> Option<i64> {
        self.bound_exceeded
    }

    /// The run of an RPS-BLAST profile, whose `.psq` letters (the redo's query) can be fewer than
    /// the columns of its preliminary alignments, asking `reads` (the profile's, as the traceback
    /// without composition-based statistics does): the near-identical walk reads the query's copy
    /// as NCBI's does (`test_near_identical_on_copy`) and reports a read past it
    /// (`TraceReads::copy_read_past`), and each alignment's start is asked before it is made
    /// (`TraceReads::starts`); the profile keeps the first such read (LOSAT's refusal of the query
    /// batch). With the letters of a consistent profile neither happens.
    pub(crate) fn with_profile_reads(mut self, reads: &'a dyn TraceReads) -> Self {
        self.profile_reads = Some(reads);
        self
    }

    /// `s_TestNearIdentical` of the window's first alignment `align` on the subject `subject`
    /// (from offset 0) and the run's query copy and words: on slices (`test_near_identical`), or
    /// as NCBI reads the copy of an RPS-BLAST profile (`with_profile_reads`).
    fn near_identical(
        &self,
        subject: &BlastCompoSequenceData,
        align: &BlastCompoAlignment,
    ) -> bool {
        let Some(reads) = self.profile_reads else {
            return test_near_identical(
                subject,
                0,
                &self.query_data,
                0,
                Some(&self.query_words),
                align,
            );
        };
        let read_past = std::cell::Cell::new(false);
        let identical = test_near_identical_on_copy(
            subject,
            0,
            NearIdenticalCopy::new(&self.query_data, &self.query_words, &read_past),
            0,
            align,
        );
        if read_past.get() {
            reads.copy_read_past();
        }
        identical
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1924-1946
    /// Whether the redone alignment of the preliminary HSP `hsp` is made (`s_RedoOneAlignment`
    /// starts at its `gapped_start`, whose row NCBI's left extension reads first): with
    /// `with_profile_reads`, as the profile answers its start (`TraceReads::starts`, read D: not
    /// at or past the matrix's `L_pin` rows, nor after the profile's first fault).
    #[inline]
    fn starts_in_rows(&self, hsp: &BlastpPreliminaryHsp) -> bool {
        match self.profile_reads {
            None => true,
            Some(reads) => profile_starts(reads, hsp),
        }
    }

    /// The matrix and the alignment scratch, for the next run (RPS-BLAST redoes one profile
    /// after the other on the same memory).
    pub(crate) fn into_buffers(self) -> (QueryPssm, PssmGapScratch) {
        (self.matrix, self.scratch)
    }
}

impl PssmRedo1<'_> {
    /// The final hit lists of the queries from their preliminary lists, which it takes.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3093-3110
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3770-3789
    pub(crate) fn run(
        &self,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Vec<Option<BlastpHitList>>> {
        if let Some(sw) = &self.sw {
            return self.run_smith_waterman(sw, preliminary_hit_lists, compact_hit_lists);
        }
        let mut run = self.redo_run();
        let mut heaps: Vec<BlastCompoHeap> = preliminary_hit_lists
            .iter()
            .map(|_| BlastCompoHeap::new(self.hitlist_size, self.inclusion_ethresh))
            .collect();
        for q_idx in 0..preliminary_hit_lists.len() {
            // A list's HSPs as the redo reads them (`s_ResultHspToDistinctAlign`), with its
            // subject and best e-value, in the order of the sorted preliminary hit list.
            let matches: Vec<(u32, f64, Vec<BlastpPreliminaryHsp>)> =
                if let Some(hit_list) = compact_hit_lists[q_idx].take() {
                    let count = hit_list.hsplist_count;
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(|list| {
                            let hsps = list.hsps.iter().map(|hsp| hsp.hsp).collect();
                            (list.oid, list.best_evalue, hsps)
                        })
                        .collect()
                } else if let Some(hit_list) = preliminary_hit_lists[q_idx].take() {
                    let count = hit_list.hsplist_count;
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(|list| {
                            let hsps = list
                                .hsps
                                .iter()
                                .map(super::traceback0::preliminary_hsp)
                                .collect();
                            (list.oid, list.best_evalue, hsps)
                        })
                        .collect()
                } else {
                    Vec::new()
                };
            for (oid, best_evalue, hsps) in matches {
                if hsps.is_empty() {
                    continue;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3532
                if early_termination(best_evalue, &heaps) {
                    continue;
                }
                let Some((hsp_list, best_score, best_evalue)) =
                    self.redo_match(&mut run, oid, &hsps)?
                else {
                    continue;
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3705-3720
                let heap = &mut heaps[q_idx];
                if heap.would_insert(best_evalue, best_score, oid) {
                    heap.insert(hsp_list, best_evalue, best_score, oid);
                }
            }
        }
        (self.redone)(&heaps)?;
        Ok(fill_results_from_compo_heaps(self.hitlist_size, &mut heaps))
    }

    /// psiblast's Smith-Waterman position-based redo loop (`-use_sw_tback`, CBS 0 to 3).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3093-3110,3525-3789
    #[cold]
    #[inline(never)]
    fn run_smith_waterman(
        &self,
        sw: &PssmSwOptions<'_>,
        preliminary_hit_lists: &mut [Option<BlastpHitList>],
        compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Vec<Option<BlastpHitList>>> {
        let scaling_factor = if sw.compo_adjust { SCALING_FACTOR } else { 1.0 };
        let mut run = self.redo_run_with_scale(scaling_factor);
        let mut sw_scratch = smith_waterman::RpsSwScratch::new();
        let mut sw_status =
            crate::core::composition_adjustment::smith_waterman_match::SwBatchStatus::default();
        let mut heaps: Vec<BlastCompoHeap> = preliminary_hit_lists
            .iter()
            .map(|_| BlastCompoHeap::new(self.hitlist_size, self.inclusion_ethresh))
            .collect();
        for q_idx in 0..preliminary_hit_lists.len() {
            let matches: Vec<(u32, f64, Vec<BlastpPreliminaryHsp>)> =
                if let Some(hit_list) = compact_hit_lists[q_idx].take() {
                    let count = hit_list.hsplist_count;
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(|list| {
                            let hsps = list.hsps.iter().map(|hsp| hsp.hsp).collect();
                            (list.oid, list.best_evalue, hsps)
                        })
                        .collect()
                } else if let Some(hit_list) = preliminary_hit_lists[q_idx].take() {
                    let count = hit_list.hsplist_count;
                    hit_list
                        .hsplist_array
                        .into_iter()
                        .take(count)
                        .map(|list| {
                            let hsps = list
                                .hsps
                                .iter()
                                .map(super::traceback0::preliminary_hsp)
                                .collect();
                            (list.oid, list.best_evalue, hsps)
                        })
                        .collect()
                } else {
                    Vec::new()
                };
            let mut status_position = None;
            for (match_position, (oid, best_evalue, hsps)) in matches.into_iter().enumerate() {
                if hsps.is_empty() {
                    continue;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3532
                if early_termination(best_evalue, &heaps) {
                    continue;
                }
                let subject_read = self.subjects.read(oid)?;
                let redone = RedoMatch {
                    info: self.info,
                    stats: self.stats,
                    expect_value: self.expect_value,
                    query_length: self.context.aa_len as i32,
                    rps_roles: false,
                    matrix: self.matrix,
                    link: self.link,
                }
                .redo_match_smith_waterman_psiblast(
                    &mut run,
                    &mut sw_scratch,
                    &mut heaps[q_idx],
                    oid,
                    &hsps,
                    subject_read.residues(),
                    sw.compo_adjust,
                );
                let Some(matched) =
                    sw_status.redone(redone, || (sw.on_endless)(&heaps, match_position))?
                else {
                    status_position = Some(match_position);
                    continue;
                };
                status_position = Some(match_position);
                let Some((hsp_list, best_score, best_evalue)) = matched else {
                    continue;
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3705-3720
                let heap = &mut heaps[q_idx];
                if heap.would_insert(best_evalue, best_score, oid) {
                    heap.insert(hsp_list, best_evalue, best_score, oid);
                }
            }
            (sw.on_query_done)(status_position, &heaps, &mut sw_status);
        }
        if sw_status.failed() {
            return Err((sw.on_batch_failed)(&heaps));
        }
        (self.redone)(&heaps)?;
        Ok(fill_results_from_compo_heaps(self.hitlist_size, &mut heaps))
    }

    /// The scaled values of the search (`s_RescaleSearch` on `kbp_gap` and the gap costs, the
    /// copy's K from `_PSIUpdateLambdaK`, `s_GetAlignParams`) and the query's data.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2133
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2437-2450
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2374-2383
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2331 (`s_GetQueryInfo`:
    /// `Blast_ReadAaComposition` and `s_CreateWordArray` of the query)
    fn redo_run(&self) -> RedoRun<'_> {
        self.redo_run_with_scale(SCALING_FACTOR)
    }

    fn redo_run_with_scale(&self, scaling_factor: f64) -> RedoRun<'_> {
        let ctx = self.context;
        let query = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
        let query_nomask_full = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
        let query_nomask = &query_nomask_full[1..query_nomask_full.len() - 1];
        let lambda = self.karlin.lambda / scaling_factor;
        RedoRun::new_with_scale(
            query,
            query_nomask,
            KarlinParams {
                lambda,
                k: self.info.gapped_k,
                ..*self.karlin
            },
            self.info.gapped_k.glibc_ln(),
            (self.gap_open, self.gap_extend),
            (self.x_drop_final_bits, self.x_drop_final),
            self.pssm.clone(),
            PssmGapScratch::new(),
            scaling_factor,
        )
    }

    /// The redo of one match of the PSSM query: the subject read from the database, then
    /// [`RedoMatch::redo_match`].
    fn redo_match(
        &self,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
    ) -> Result<Option<(BlastpHspList, i32, f64)>> {
        let subject_read = self.subjects.read(oid)?;
        Ok(RedoMatch {
            info: self.info,
            stats: self.stats,
            expect_value: self.expect_value,
            query_length: self.context.aa_len as i32,
            rps_roles: false,
            matrix: self.matrix,
            link: self.link,
        }
        .redo_match(run, oid, hsps, subject_read.residues()))
    }
}

/// The subject of the X-drop redo (NCBI's matching sequence), and the letters its HSPs are
/// counted on: a protein sequence (psiblast's subject, rpsblast's query: [`ProteinSubject`]) or
/// a translated sequence's six frames (rpstblastn's query: [`FrameSubject`]), whose letters are
/// NCBI's translation buffer (`one_query->sequence`, `one_query->length` bytes: the frames in
/// the context order, one sentinel between two frames with letters).
pub(crate) trait RedoSubject<'a>: Copy {
    /// The instance of the steps the redos share.
    type Instance: RedoInstance;
    /// The letters an HSP of `frame` aligns to, as the redo's subject `subject` holds them.
    fn frame_letters(self, subject: &'a [u8], frame: i32) -> &'a [u8];
    /// The redo of one match of this subject: `RedoMatch::redo_match` (a protein subject) or
    /// `RedoMatch::redo_windows` (frames).
    fn redo(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Option<(BlastpHspList, i32, f64)>;
    /// The Smith-Waterman redo of one match of this subject (`-use_sw_tback`):
    /// `RedoMatch::redo_match_smith_waterman` (a protein subject) or
    /// `RedoMatch::redo_windows_smith_waterman` (frames). `heap` is the query's heap of redone
    /// matches.
    #[allow(clippy::too_many_arguments)]
    fn redo_smith_waterman(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        sw: &mut smith_waterman::RpsSwScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Result<Option<(BlastpHspList, i32, f64)>, smith_waterman::SwEndless>;
}

/// A protein subject ([`RedoSubject`]): the redo's subject is its letters.
#[derive(Clone, Copy)]
pub(crate) struct ProteinSubject;

impl<'a> RedoSubject<'a> for ProteinSubject {
    type Instance = XDropRedo;

    #[inline(always)]
    fn frame_letters(self, subject: &'a [u8], _frame: i32) -> &'a [u8] {
        subject
    }

    fn redo(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Option<(BlastpHspList, i32, f64)> {
        matching.redo_match(run, oid, hsps, subject)
    }

    #[inline(always)]
    fn redo_smith_waterman(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        sw: &mut smith_waterman::RpsSwScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Result<Option<(BlastpHspList, i32, f64)>, smith_waterman::SwEndless> {
        matching.redo_match_smith_waterman(run, sw, heap, oid, hsps, subject)
    }
}

/// A translated subject ([`RedoSubject`]): the letters of each frame, without sentinels, by
/// context (`BLAST_FrameToContext`).
#[derive(Clone, Copy)]
pub(crate) struct FrameSubject<'a> {
    pub(crate) frames: &'a [&'a [u8]],
}

impl<'a> RedoSubject<'a> for FrameSubject<'a> {
    type Instance = FrameRedo;

    fn frame_letters(self, _subject: &'a [u8], frame: i32) -> &'a [u8] {
        self.frames[crate::core::blast_util::frame_to_context(frame)]
    }

    fn redo(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Option<(BlastpHspList, i32, f64)> {
        matching.redo_windows(run, oid, hsps, subject, self)
    }

    #[inline(always)]
    fn redo_smith_waterman(
        self,
        matching: &RedoMatch<'_>,
        run: &mut RedoRun<'_>,
        sw: &mut smith_waterman::RpsSwScratch,
        heap: &mut BlastCompoHeap,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Result<Option<(BlastpHspList, i32, f64)>, smith_waterman::SwEndless> {
        matching.redo_windows_smith_waterman(run, sw, heap, oid, hsps, subject, self)
    }
}

/// The identities and positives that NCBI's RPS-TBLASTN redo counts for an HSP
/// (`s_ComputeNumIdentities`, which the Seq-align's `num_ident` and `num_positives` show): on the
/// redo's whole subject, the translation buffer (`one_query->sequence`), at the HSP's offsets in
/// its own frame, not on the frame's letters, so that the counts of an HSP of any frame but +1
/// are those of the letters of frame +1 (then the sentinel and the frames after it) at the
/// same offsets. The query side is the profile's consensus.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:482-531
/// (`Blast_RedoAlignmentCore` is called with the profile as `query_blk` and `one_query`, the
/// user's translated query of all frames, as `subject_blk`, blast_traceback.c:1274-1281; the
/// HSP's subject offsets are in its frame.) The offsets stay inside the buffer, which holds
/// every frame, the HSP's own among them. `C`, the redo's instance, keeps each instance's own
/// copy (the X-drop redo's code as before the Smith-Waterman redo of frames; DECISIONS 63 S3).
fn translation_buffer_counts<C: RedoInstance>(
    hsp: &BlastpHsp,
    query_nomask: &[u8],
    buffer: &[u8],
    matrix: crate::config::ScoringMatrix,
) -> (usize, usize) {
    let (mut q, mut s) = (hsp.q_start - 1, hsp.s_start.min(hsp.s_end) - 1);
    let implicit = [GapEditOp::Sub(hsp.q_end.saturating_sub(q) as u32)];
    let ops = hsp.gap_info.as_deref().unwrap_or(&implicit);
    let (mut num_ident, mut num_positives) = (0usize, 0usize);
    for op in ops {
        match *op {
            GapEditOp::Sub(n) => {
                for (&ql, &sl) in query_nomask[q..q + n as usize]
                    .iter()
                    .zip(&buffer[s..s + n as usize])
                {
                    if ql == sl {
                        num_ident += 1;
                        num_positives += 1;
                    } else if crate::utils::matrix::standard_score(matrix, ql, sl) > 0 {
                        num_positives += 1;
                    }
                }
                q += n as usize;
                s += n as usize;
            }
            GapEditOp::Del(n) => s += n as usize,
            GapEditOp::Ins(n) => q += n as usize,
        }
    }
    (num_ident, num_positives)
}

/// A window of the subject (`s_WindowInfo`'s `subject_range`): `[begin, end)` of the frame
/// `frame`'s letters. `end` may lie before `begin` (`end - begin` is NCBI's length, then below 0:
/// `rpstblastn -strand`, whose translation buffer is shorter than `GET_NUCL_LENGTH` assumes). A
/// protein subject's window (the whole sequence) is never read.
#[derive(Clone, Copy, Default)]
struct SubjectWindow {
    begin: i32,
    end: i32,
    frame: i32,
}

/// What the redo of one match reads besides its run: the matrix information, the Gumbel block
/// and the search space (`EvalueStats`), `-evalue`, the query's length and the roles of the
/// lengths in the e-values (RPS-BLAST's query is the profile and its subject the user's query:
/// `Blast_HSPListGetEvalues` swaps them).
pub(crate) struct RedoMatch<'a> {
    pub(crate) info: &'a PosMatrixInfo,
    pub(crate) stats: EvalueStats<'a>,
    pub(crate) expect_value: f64,
    pub(crate) query_length: i32,
    pub(crate) rps_roles: bool,
    /// The matrix of the search (`sbp->matrix`), for the positives.
    pub(crate) matrix: crate::config::ScoringMatrix,
    /// psiblast's `-sum_stats true` (`RedoLink`); `None` without sum statistics and for
    /// RPS-BLAST, which refuses it.
    pub(crate) link: Option<RedoLink>,
}

impl RedoMatch<'_> {
    /// Reaps redone HSPs whose E-value exceeds `-evalue` (keeping NaN E-values, because C's
    /// greater-than test is false for NaN) and verifies that the list's `best_evalue` is at most
    /// `-evalue`, returning the first surviving HSP's scaled score or `None` when none passes.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1976-2005
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:437-443
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3687
    #[inline(always)]
    fn reap_by_expect(&self, redone: &mut Vec<Redone>, best_evalue: f64) -> Option<i32> {
        redone.retain(|hsp| hsp.e_value <= self.expect_value || hsp.e_value.is_nan());
        let first = redone.first()?;
        // The list's best E-value must be at most `-evalue`: NaN E-values (a PSSM whose gapped
        // block gives them) pass the reaping above but leave it at `INT4_MAX`.
        if best_evalue > self.expect_value {
            return None;
        }
        Some(first.align.score)
    }

    /// The redo of one match (the HSPs of the query and one subject): `Blast_RedoOneMatch`, then
    /// `s_HSPListFromDistinctAlignments`, `s_HitlistReapContained`, `s_HitlistEvaluateAndPurge`,
    /// and for a list whose best e-value passes `-evalue`, `s_HSPListNormalizeScores` and
    /// `s_ComputeNumIdentities`. `None` when no HSP passes.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3670-3707
    pub(crate) fn redo_match(
        &self,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Option<(BlastpHspList, i32, f64)> {
        let alignments = self.redo_one_match(run, hsps, subject);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:321-357
        // The sort is glibc's `qsort`, a stable merge sort (see `blastp::hsp::sort_hsplist_by_score`).
        let mut redone: Vec<Redone> = alignments
            .into_iter()
            .map(|align| Redone {
                align,
                e_value: 0.0,
                num: 1,
            })
            .collect();
        redone.sort_by(|a, b| score_compare(&a.align, &b.align));
        if redone.len() > 1 {
            reap_contained::<XDropRedo>(&mut redone);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:409-446
        // psiblast's e-values are not adjusted by a composition p-value (that is blastp's and
        // blastx's, with `-comp_based_stats 2u`/`3u`).
        let query_length = self.query_length;
        let subject_length = subject.len() as i32;
        let best_evalue = match &self.link {
            Some(link) => link_redone(&mut redone, link, run, subject_length)?,
            None => {
                let mut best_evalue = i32::MAX as f64;
                for hsp in &mut redone {
                    // (RPS-BLAST's roles are swapped, `Blast_HSPListGetEvalues`: its `isRPS`.)
                    let (m, n) = if self.rps_roles {
                        (subject_length, query_length)
                    } else {
                        (query_length, subject_length)
                    };
                    hsp.e_value = self.stats.e_value(hsp.align.score, &run.karlin, m, n);
                    best_evalue = best_evalue.min(hsp.e_value);
                }
                best_evalue
            }
        };
        let best_score = self.reap_by_expect(&mut redone, best_evalue)?;
        // The list passed NCBI's test that `best_evalue` is at most
        // `hitParams->options->expect_value` (above).
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:102-116
        let template = &hsps[0];
        let mut hsp_list = BlastpHspList {
            oid,
            query_index: template.q_idx,
            hsps: Vec::with_capacity(redone.len()),
            best_evalue,
        };
        for Redone {
            align,
            e_value,
            num,
        } in redone
        {
            let score = blast_nint(align.score as f64 / SCALING_FACTOR);
            let traced = BlastpPreliminaryHsp {
                query_start: align.query_start,
                query_end: align.query_end,
                subject_start: align.match_start,
                subject_end: align.match_end,
                gapped_query_start: align.query_start,
                gapped_subject_start: align.match_start,
                raw_score: score,
                query_context: template.q_idx as i32,
                query_frame: template.query_frame,
                subject_frame: 0,
                query_length: template.query_length,
                q_idx: template.q_idx,
                s_idx: template.s_idx,
            };
            let mut hsp = final_hsp(
                &traced,
                align.edit_script,
                run.query_nomask,
                subject,
                self.matrix,
            );
            hsp.e_value = e_value;
            hsp.num = num;
            hsp.bit_score = (score as f64 * run.karlin.lambda * SCALING_FACTOR - run.log_k)
                / std::f64::consts::LN_2;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:332-342
            hsp.comp_adjustment_method = match align.rule {
                EMatrixAdjustRule::DontAdjustMatrix => 0,
                EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                _ => 2,
            };
            hsp_list.hsps.push(hsp);
        }
        Some((hsp_list, best_score, best_evalue))
    }

    /// `Blast_RedoOneMatch` for the one window of a protein query and subject (the whole query,
    /// the whole subject, `s_WindowsFromProteinAligns`): the distinct redone alignments, the
    /// list's head first.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1173-1286
    /// `params->cutoff_s` is 1 without sum statistics, else `RedoLink::cutoff_s`
    /// (blast_kappa.c:2453-2461).
    fn redo_one_match(
        &self,
        run: &mut RedoRun<'_>,
        hsps: &[BlastpPreliminaryHsp],
        subject: &[u8],
    ) -> Vec<Align> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:788-800
        let incoming: Vec<Align> = hsps
            .iter()
            .map(|hsp| Align {
                score: (hsp.raw_score as f64 * SCALING_FACTOR) as i32,
                rule: EMatrixAdjustRule::DontAdjustMatrix,
                query_start: hsp.query_start,
                query_end: hsp.query_end,
                match_start: hsp.subject_start,
                match_end: hsp.subject_end,
                frame: hsp.subject_frame,
                edit_script: Vec::new(),
            })
            .collect();
        let lambda = run.karlin.lambda;
        let query_length = run.query.len() as i32;
        let mut alignments: Vec<Align> = Vec::new();
        let mut matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
        let mut subject_data: Vec<u8> = Vec::new();
        let mut near_identical = true;
        let mut old_near_identical = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;
        for (hsp_index, in_align) in incoming.iter().enumerate() {
            if hsp_index == 0 || subject_maybe_biased {
                near_identical = preliminary_test_near_identical(
                    query_length,
                    in_align,
                    run.near_identical_cutoff,
                );
            }
            if hsp_index == 0 || (subject_maybe_biased && near_identical != old_near_identical) {
                subject_data = self.subject_range::<XDropRedo>(
                    run,
                    subject,
                    &incoming[0],
                    near_identical,
                    &mut subject_maybe_biased,
                );
            }
            if !is_contained::<XDropRedo>(in_align, &alignments, lambda) {
                let mut adjust_search_failed = false;
                if hsp_index == 0 || near_identical != old_near_identical {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1436-1444
                    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1484-1530
                    let subject_composition = read_aa_composition(&subject_data);
                    if run.query_true_amino_acids == 0
                        || subject_composition.num_true_amino_acids == 0
                    {
                        adjust_search_failed = true;
                    } else {
                        matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
                        pssm_composition_based_stats::<XDropRedo>(
                            &mut run.matrix.rows,
                            self.info,
                            &subject_composition.prob,
                        );
                        if run.bound_check && run.bound_exceeded.is_none() {
                            let bound: i64 = run
                                .matrix
                                .rows
                                .iter()
                                .map(|row| i64::from(row.iter().copied().max().unwrap_or(0).max(0)))
                                .sum();
                            if bound > i64::from(i32::MAX) {
                                run.bound_exceeded = Some(bound);
                            }
                        }
                    }
                    num_adjustments += 1;
                }
                if !adjust_search_failed {
                    let new_align = self.redo_one_alignment(
                        run,
                        &hsps[hsp_index],
                        matrix_adjust_rule,
                        &subject_data,
                    );
                    if new_align.score >= self.link.map_or(1, |link| link.cutoff_s) {
                        with_distinct_ends::<XDropRedo>(
                            new_align,
                            &mut alignments,
                            num_adjustments == 1,
                        );
                    }
                }
            }
            old_near_identical = near_identical;
        }
        alignments
    }

    /// The redo of one match of a translated subject (rpstblastn): `redo_match`'s steps (quoted
    /// there) over the windows of the frames (`redo_frame_windows`), with the frames compared in
    /// the reap, the subject's length of `s_GetSubjectLength` in the e-values and the frames kept
    /// in the HSPs. A copy (DECISIONS 63 S3): psiblast's and rpsblast's `redo_match` stays dev's
    /// code, which one generic form shared with the frames did not keep instruction for
    /// instruction; its instance of a protein subject is `redo_match`'s for the equivalence test
    /// (`the_window_copy_is_the_x_drop_redo`).
    pub(crate) fn redo_windows<'a, S: RedoSubject<'a>>(
        &self,
        run: &mut RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
        kind: S,
    ) -> Option<(BlastpHspList, i32, f64)> {
        // A protein subject is one window, the whole subject (`s_WindowsFromProteinAligns`).
        let alignments = if S::Instance::SUBJECT_IS_TRANSLATED {
            self.redo_frame_windows::<S>(run, hsps, subject)
        } else {
            self.redo_window::<S>(run, hsps, subject, SubjectWindow::default(), None)
        };
        self.evaluate_windows::<S::Instance, S>(run, oid, hsps, subject, kind, alignments)
    }

    /// The steps of the redo of a match after its windows' alignments `alignments` (the list's
    /// head first), as `redo_windows` quotes them: `s_HSPListFromDistinctAlignments`,
    /// `s_HitlistReapContained` (frames compared for a translated subject),
    /// `s_HitlistEvaluateAndPurge` (the subject's length of `s_GetSubjectLength`) and, for a list
    /// whose best e-value passes `-evalue`, `s_HSPListNormalizeScores` and
    /// `s_ComputeNumIdentities`. NCBI runs them after either redo of the windows
    /// (`Blast_RedoOneMatch` or `Blast_RedoOneMatchSmithWaterman`, blast_kappa.c:3602-3707):
    /// the X-drop redo (`redo_windows`) and the Smith-Waterman redo of frames
    /// (`redo_windows_smith_waterman`) share them, each as its own instance `C`, inlined, so that
    /// the X-drop redo's code stays as it was before the second caller (DECISIONS 63 S3).
    #[inline(always)]
    fn evaluate_windows<'a, C: RedoInstance, S: RedoSubject<'a>>(
        &self,
        run: &RedoRun<'_>,
        oid: u32,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
        kind: S,
        alignments: Vec<Align>,
    ) -> Option<(BlastpHspList, i32, f64)> {
        let mut redone: Vec<Redone> = alignments
            .into_iter()
            .map(|align| Redone {
                align,
                e_value: 0.0,
                num: 1,
            })
            .collect();
        redone.sort_by(|a, b| score_compare(&a.align, &b.align));
        if redone.len() > 1 {
            reap_contained::<C>(&mut redone);
        }
        let query_length = self.query_length;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:364-368
        let subject_length = if C::SUBJECT_IS_TRANSLATED {
            (get_nucl_length(subject.len() as i32) - 1) / 3
        } else {
            subject.len() as i32
        };
        let best_evalue = match &self.link {
            Some(link) => {
                // The RPS redo has no sum statistics (`rps_search::redo`): only psiblast links.
                assert!(
                    !C::SUBJECT_IS_TRANSLATED,
                    "the redo of frames has no sum statistics"
                );
                link_redone(&mut redone, link, run, subject_length)?
            }
            None => {
                let mut best_evalue = i32::MAX as f64;
                for hsp in &mut redone {
                    let (m, n) = if self.rps_roles {
                        (subject_length, query_length)
                    } else {
                        (query_length, subject_length)
                    };
                    hsp.e_value = self.stats.e_value(hsp.align.score, &run.karlin, m, n);
                    best_evalue = best_evalue.min(hsp.e_value);
                }
                best_evalue
            }
        };
        let best_score = self.reap_by_expect(&mut redone, best_evalue)?;
        let template = &hsps[0];
        let mut hsp_list = BlastpHspList {
            oid,
            query_index: template.q_idx,
            hsps: Vec::with_capacity(redone.len()),
            best_evalue,
        };
        for Redone {
            align,
            e_value,
            num,
        } in redone
        {
            let score = blast_nint(align.score as f64 / SCALING_FACTOR);
            let traced = BlastpPreliminaryHsp {
                query_start: align.query_start,
                query_end: align.query_end,
                subject_start: align.match_start,
                subject_end: align.match_end,
                gapped_query_start: align.query_start,
                gapped_subject_start: align.match_start,
                raw_score: score,
                query_context: template.q_idx as i32,
                query_frame: template.query_frame,
                // (`Blast_HSPInit` takes `align->frame` cast to `Int2`, blast_kappa.c:326-330.)
                subject_frame: if C::SUBJECT_IS_TRANSLATED {
                    align.frame
                } else {
                    0
                },
                query_length: template.query_length,
                q_idx: template.q_idx,
                s_idx: template.s_idx,
            };
            let mut hsp = final_hsp(
                &traced,
                align.edit_script,
                run.query_nomask,
                kind.frame_letters(subject, align.frame),
                self.matrix,
            );
            if C::SUBJECT_IS_TRANSLATED {
                // The Seq-align's counts (see `translation_buffer_counts`).
                (hsp.num_ident, hsp.num_positives) =
                    translation_buffer_counts::<C>(&hsp, run.query_nomask, subject, self.matrix);
            }
            hsp.e_value = e_value;
            hsp.num = num;
            hsp.bit_score = (score as f64 * run.karlin.lambda * SCALING_FACTOR - run.log_k)
                / std::f64::consts::LN_2;
            hsp.comp_adjustment_method = match align.rule {
                EMatrixAdjustRule::DontAdjustMatrix => 0,
                EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                _ => 2,
            };
            hsp_list.hsps.push(hsp);
        }
        Some((hsp_list, best_score, best_evalue))
    }

    /// `Blast_RedoOneMatch` of a translated subject: the distinct redone alignments, the list's
    /// head first. (A protein subject is one window, the whole subject with the alignments in
    /// their order, `s_WindowsFromProteinAligns`: [`Self::redo_window`].) The windows are the
    /// alignments' frames around them (`s_WindowsFromTranslatedAligns` with a position-based
    /// matrix, `frame_windows`), each with its alignments in `s_DistinctAlignmentsSort`'s order,
    /// and all windows add to the one list (`alignments[query_index]`). NCBI keeps
    /// `matrix_adjust_rule` from one window to the next; a translated subject's matrix is
    /// adjusted for each alignment that is redone (`redo_window`), so the rule a window starts
    /// with is never read.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1120-1143
    fn redo_frame_windows<'a, S: RedoSubject<'a>>(
        &self,
        run: &mut RedoRun<'_>,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
    ) -> Vec<Align> {
        let mut alignments: Vec<Align> = Vec::new();
        for (window, order) in
            frame_windows::<S::Instance>(hsps, run.query.len() as i32, subject.len() as i32)
        {
            let window_hsps: Vec<BlastpPreliminaryHsp> = order.iter().map(|&k| hsps[k]).collect();
            alignments =
                self.redo_window::<S>(run, &window_hsps, subject, window, Some(alignments));
        }
        alignments
    }

    /// The redo of one window of `Blast_RedoOneMatch` (`redo_one_match`'s steps, cited there;
    /// `window`, whose alignments are `hsps` in their order): the list of the windows before
    /// (`before`, none for the first) with the window's distinct redone alignments at its head. A
    /// translated subject's range is the window of its frame (`frame_range`), and the matrix is
    /// adjusted for each alignment that is not contained, to the composition around it in the
    /// window (`s_GetComposition`, `window_composition`).
    fn redo_window<'a, S: RedoSubject<'a>>(
        &self,
        run: &mut RedoRun<'_>,
        hsps: &[BlastpPreliminaryHsp],
        subject: &'a [u8],
        window: SubjectWindow,
        before: Option<Vec<Align>>,
    ) -> Vec<Align> {
        // The incoming alignments (`Align::incoming`) in the window's order.
        let incoming: Vec<Align> = hsps.iter().map(Align::incoming::<S::Instance>).collect();
        let lambda = run.karlin.lambda;
        let query_length = run.query.len() as i32;
        let mut alignments: Vec<Align> = before.unwrap_or_default();
        let mut matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
        let mut subject_data: Vec<u8> = Vec::new();
        let mut near_identical = true;
        let mut old_near_identical = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;
        for (hsp_index, in_align) in incoming.iter().enumerate() {
            if hsp_index == 0 || subject_maybe_biased {
                // Only rpstblastn's redo (the translated instance) can have the cutoff 0
                // (`RedoRun::without_valid_context`); psiblast's and rpsblast's keep the one test.
                near_identical =
                    if S::Instance::SUBJECT_IS_TRANSLATED && run.near_identical_cutoff <= 0.0 {
                        old_test_near_identical(query_length, incoming.len(), in_align)
                    } else {
                        preliminary_test_near_identical(
                            query_length,
                            in_align,
                            run.near_identical_cutoff,
                        )
                    };
            }
            if hsp_index == 0 || (subject_maybe_biased && near_identical != old_near_identical) {
                subject_data = if S::Instance::SUBJECT_IS_TRANSLATED {
                    frame_range::<S::Instance>(
                        run,
                        subject,
                        window,
                        &incoming[0],
                        near_identical,
                        &mut subject_maybe_biased,
                    )
                } else {
                    self.subject_range::<S::Instance>(
                        run,
                        subject,
                        &incoming[0],
                        near_identical,
                        &mut subject_maybe_biased,
                    )
                };
            }
            if !is_contained::<S::Instance>(in_align, &alignments, lambda) {
                let mut adjust_search_failed = false;
                if S::Instance::SUBJECT_IS_TRANSLATED
                    || hsp_index == 0
                    || near_identical != old_near_identical
                {
                    let subject_composition = if S::Instance::SUBJECT_IS_TRANSLATED {
                        window_composition::<S::Instance>(&subject_data, window, in_align)
                    } else {
                        read_aa_composition(&subject_data)
                    };
                    if run.query_true_amino_acids == 0
                        || subject_composition.num_true_amino_acids == 0
                    {
                        adjust_search_failed = true;
                    } else {
                        matrix_adjust_rule = EMatrixAdjustRule::CompoScaleOldMatrix;
                        pssm_composition_based_stats::<S::Instance>(
                            &mut run.matrix.rows,
                            self.info,
                            &subject_composition.prob,
                        );
                        if run.bound_check && run.bound_exceeded.is_none() {
                            let bound: i64 = run
                                .matrix
                                .rows
                                .iter()
                                .map(|row| i64::from(row.iter().copied().max().unwrap_or(0).max(0)))
                                .sum();
                            if bound > i64::from(i32::MAX) {
                                run.bound_exceeded = Some(bound);
                            }
                        }
                    }
                    num_adjustments += 1;
                }
                if !adjust_search_failed {
                    let new_align = if S::Instance::SUBJECT_IS_TRANSLATED {
                        redo_window_alignment(
                            run,
                            &hsps[hsp_index],
                            matrix_adjust_rule,
                            &subject_data,
                            window,
                        )
                    } else {
                        self.redo_one_alignment(
                            run,
                            &hsps[hsp_index],
                            matrix_adjust_rule,
                            &subject_data,
                        )
                    };
                    if new_align.score >= self.link.map_or(1, |link| link.cutoff_s) {
                        with_distinct_ends::<S::Instance>(
                            new_align,
                            &mut alignments,
                            num_adjustments == 1,
                        );
                    }
                }
            }
            old_near_identical = near_identical;
        }
        alignments
    }

    /// `s_SequenceGetProteinRange` for the whole subject: its letters, masked by SEG ("S 10 1.8
    /// 2.1") unless the window's first alignment is near identical; `subject_maybe_biased`
    /// becomes whether SEG masked anything.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1625-1638
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1413-1416
    fn subject_range<C: RedoInstance>(
        &self,
        run: &RedoRun<'_>,
        subject: &[u8],
        window_align: &Align,
        near_identical: bool,
        subject_maybe_biased: &mut bool,
    ) -> Vec<u8> {
        let mut data = subject.to_vec();
        if *subject_maybe_biased {
            let test_identical = near_identical && {
                let align = BlastCompoAlignment {
                    score: window_align.score,
                    matrix_adjust_rule: window_align.rule,
                    query_index: 0,
                    query_start: window_align.query_start,
                    query_end: window_align.query_end,
                    match_start: window_align.match_start,
                    match_end: window_align.match_end,
                    frame: window_align.frame,
                    context: None,
                    next: None,
                };
                run.near_identical(&BlastCompoSequenceData::from_ncbistdaa(subject), &align)
            };
            if !test_identical {
                let intervals =
                    crate::utils::seg::SegMasker::new(10, 1.8, 2.1).mask_sequence(&data);
                *subject_maybe_biased = !intervals.is_empty();
                for interval in &intervals {
                    data[interval.start..interval.end].fill(ncbistdaa::X);
                }
            }
        }
        data
    }

    /// `s_RedoOneAlignment`: the gapped alignment with traceback of the preliminary HSP `hsp`
    /// from its gapped start, on the PSSM of the last adjustment, with the scaled gap costs and
    /// X-drop (`s_NewAlignmentFromGapAlign` of the whole query and subject).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1924-1946
    fn redo_one_alignment(
        &self,
        run: &mut RedoRun<'_>,
        hsp: &BlastpPreliminaryHsp,
        rule: EMatrixAdjustRule,
        subject: &[u8],
    ) -> Align {
        if !run.starts_in_rows(hsp) {
            return Align::not_made(rule);
        }
        let traced = pssm_gapped_alignment_with_traceback(
            &run.matrix,
            subject,
            subject.len(),
            hsp.gapped_query_start as usize,
            hsp.gapped_subject_start as usize,
            run.gap_open,
            run.gap_extend,
            run.x_dropoff,
            &mut run.scratch,
        );
        Align {
            score: traced.score,
            rule,
            query_start: traced.query_start,
            query_end: traced.query_stop,
            match_start: traced.subject_start,
            match_end: traced.subject_stop,
            frame: 0,
            edit_script: traced.edit_script,
        }
    }
}

/// `s_DoSegSequenceData` of the redo's subject letters `data` ("S 10 1.8 2.1", the masked
/// letters X): `subject_maybe_biased` becomes whether SEG masked anything.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1435-1446
#[inline(always)]
fn seg_subject(data: &mut [u8], subject_maybe_biased: &mut bool) {
    let intervals = crate::utils::seg::SegMasker::new(10, 1.8, 2.1).mask_sequence(data);
    *subject_maybe_biased = !intervals.is_empty();
    for interval in &intervals {
        data[interval.start..interval.end].fill(ncbistdaa::X);
    }
}

/// `s_WindowsFromTranslatedAligns` of the incoming alignments of `hsps` with a position-based
/// matrix (`windows_from_translated_aligns`, `is_pos_based`) in a translated subject of
/// `subject_length` bytes (the translation buffer), for a query of `query_length` letters:
/// each window and the indices in `hsps` of its alignments, in the order the redo reads them.
/// `C`, the redo's instance, keeps each instance's own copy (`translation_buffer_counts`).
fn frame_windows<C: RedoInstance>(
    hsps: &[BlastpPreliminaryHsp],
    query_length: i32,
    subject_length: i32,
) -> Vec<(SubjectWindow, Vec<usize>)> {
    let mut list: Option<Box<BlastCompoAlignment>> = None;
    for (k, hsp) in hsps.iter().enumerate().rev() {
        let mut align = Box::new(
            Align::incoming::<C>(hsp)
                .compo(Some(BlastCompoAlignmentContext::PreliminaryHspIndex(k))),
        );
        align.next = list;
        list = Some(align);
    }
    let windows =
        windows_from_translated_aligns(&list, &[query_length], subject_length, WINDOW_BORDER, true)
            .expect("the redo's alignments are of query 0");
    windows
        .into_iter()
        .map(|window| {
            let mut order = Vec::with_capacity(window.hspcnt as usize);
            let mut current = window.align.as_deref();
            while let Some(align) = current {
                if let Some(BlastCompoAlignmentContext::PreliminaryHspIndex(k)) = align.context {
                    order.push(k);
                }
                current = align.next.as_deref();
            }
            let range = window.subject_range;
            (
                SubjectWindow {
                    begin: range.begin,
                    end: range.end,
                    frame: range.context,
                },
                order,
            )
        })
        .collect()
}

/// `s_SequenceGetProteinRange` of a translated subject (`self->index` below 0, `align->frame`
/// not 0) for the window `window` whose first alignment is `window_align`: NCBI's buffer of the
/// translation buffer's length plus 2 bytes (`calloc`), whose second byte on holds the frame's
/// letters as NCBI finds them in `subject` (`GET_TRANSLATED_LENGTH` letters of `GET_NUCL_LENGTH`
/// of the buffer's length, after the earlier frames' lengths plus a sentinel each), masked by SEG
/// unless the window's first alignment is near identical (`s_TestNearIdentical` on that buffer);
/// then the byte before the window is zeroed (the fit: the window's letters start at byte
/// `1 + window.begin`).
///
/// Every read is inside `subject` and the buffer, also where NCBI's frame arithmetic does not
/// match the buffer (`-strand plus|minus`, whose buffer holds three frames): with `nl` the
/// `GET_NUCL_LENGTH` of `l`, the three lengths `GET_TRANSLATED_LENGTH` of `nl` and `f` for `f` =
/// 0, 1, 2 add up to `nl - 2`, so the copy of frame index 5 ends at
/// `2 * (nl - 2) + 5 = 2 * ((l - 5) / 2) + 5 <= l`; a frame has letters only where `nl >= 3`,
/// `l >= 7`. The window's alignments end inside their frame, so before `l`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1597-1643
/// `C`, the redo's instance, keeps each instance's own copy (`translation_buffer_counts`).
fn frame_range<C: RedoInstance>(
    run: &RedoRun<'_>,
    subject: &[u8],
    window: SubjectWindow,
    window_align: &Align,
    near_identical: bool,
    subject_maybe_biased: &mut bool,
) -> Vec<u8> {
    let length = subject.len();
    let mut buffer = vec![0u8; length + 2];
    let f = get_seq_frame(window_align.frame);
    let nucl_length = get_nucl_length(length as i32);
    let frame_length = get_translated_length(nucl_length, f) as usize;
    if frame_length > 0 {
        let offset: usize = (0..f)
            .map(|i| get_translated_length(nucl_length, i) as usize + 1)
            .sum();
        buffer[1..1 + frame_length].copy_from_slice(&subject[offset..offset + frame_length]);
    }
    if *subject_maybe_biased {
        // `s_TestNearIdentical` reads the buffer (not only the frame's letters) at the window's
        // first alignment.
        let data = BlastCompoSequenceData {
            buffer,
            data_offset: 1,
            length: length as i32 + 1,
        };
        let test_identical = near_identical && run.near_identical(&data, &window_align.compo(None));
        buffer = data.buffer;
        if !test_identical {
            seg_subject(&mut buffer[1..1 + frame_length], subject_maybe_biased);
        }
    }
    buffer[window.begin as usize] = 0;
    buffer
}

/// `s_GetComposition` of a translated subject: the composition of the letters of `window`
/// (in `data`, `frame_range`'s buffer) around the alignment `align` (`Blast_GetCompositionRange`
/// of its range in the window's `end - begin` letters). Where the window holds fewer letters
/// than the alignment reaches (`-strand plus|minus`), NCBI's range ends at the alignment's end
/// (the scan to the right stops at once): the letters read here reach that end, which stays
/// inside the buffer (`frame_range`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:921-936
/// `C`, the redo's instance, keeps each instance's own copy (`translation_buffer_counts`).
fn window_composition<C: RedoInstance>(
    data: &[u8],
    window: SubjectWindow,
    align: &Align,
) -> BlastAminoAcidComposition {
    let length = (window.end - window.begin).max(0) as usize;
    let start = (align.match_start - window.begin) as usize;
    let finish = (align.match_end - window.begin) as usize;
    let letters = &data[1 + window.begin as usize..][..length.max(finish)];
    let (left, right) = blast_get_composition_range(letters, start, finish)
        .expect("an alignment's range lies in the letters read");
    read_aa_composition(&letters[left..right])
}

/// The profile's answer for the gapped start of the preliminary HSP `hsp`
/// (`RedoRun::starts_in_rows`, `TraceReads::starts`), out of the redo's loops.
#[inline(never)]
fn profile_starts(reads: &dyn TraceReads, hsp: &BlastpPreliminaryHsp) -> bool {
    reads.starts(usize::try_from(hsp.gapped_query_start).unwrap_or(usize::MAX))
}

/// `s_RedoOneAlignment` in a window of a translated subject: the gapped alignment of the
/// preliminary HSP `hsp` from its gapped start on the window's letters (in `data`,
/// `frame_range`'s buffer; `end - begin` of them, none where the window ends before it begins:
/// then only the left extension runs, NCBI's `s_start < s_length` failing), shifted back into the
/// frame (`s_NewAlignmentFromGapAlign`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1925-1926
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760-1765
fn redo_window_alignment(
    run: &mut RedoRun<'_>,
    hsp: &BlastpPreliminaryHsp,
    rule: EMatrixAdjustRule,
    data: &[u8],
    window: SubjectWindow,
) -> Align {
    if !run.starts_in_rows(hsp) {
        return Align::not_made(rule);
    }
    let traced = pssm_gapped_alignment_with_traceback(
        &run.matrix,
        &data[1 + window.begin as usize..],
        (window.end - window.begin).max(0) as usize,
        hsp.gapped_query_start as usize,
        (hsp.gapped_subject_start - window.begin) as usize,
        run.gap_open,
        run.gap_extend,
        run.x_dropoff,
        &mut run.scratch,
    );
    Align {
        score: traced.score,
        rule,
        query_start: traced.query_start,
        query_end: traced.query_stop,
        match_start: traced.subject_start + window.begin,
        match_end: traced.subject_stop + window.begin,
        frame: window.frame,
        edit_script: traced.edit_script,
    }
}

/// `BLAST_LinkHsps` on a match's redone list (psiblast's `-sum_stats true`): every HSP takes its
/// linked set's e-value and size, in the list the linking builds, sorted by score; the best
/// e-value of the list, `None` for an empty list. The scores are the redo's scaled ones and the
/// linking reads the redo's scaled Lambda and logK (`RedoLink::context`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1801-1811
fn link_redone(
    redone: &mut Vec<Redone>,
    link: &RedoLink,
    run: &RedoRun<'_>,
    subject_length: i32,
) -> Option<f64> {
    if redone.is_empty() {
        return None;
    }
    let context = LinkContext {
        lambda: run.karlin.lambda,
        log_k: run.log_k,
        ..link.context
    };
    link_in_place(
        redone,
        |hsp| LinkHsp {
            context: 0,
            subject_frame: 0,
            query_offset: hsp.align.query_start,
            query_end: hsp.align.query_end,
            subject_offset: hsp.align.match_start,
            subject_end: hsp.align.match_end,
            score: hsp.align.score,
        },
        |_| context,
        subject_length,
        &link.params,
        |hsp, linked| {
            hsp.e_value = linked.evalue;
            hsp.num = linked.num;
        },
    );
    // `Blast_HSPListSortByScore`: glibc's `qsort`, a stable merge sort, when the list is not
    // sorted (a stable sort of a sorted list keeps it).
    redone.sort_by(|a, b| score_compare(&a.align, &b.align));
    let best_evalue = redone[1..].iter().fold(redone[0].e_value, |best, hsp| {
        if hsp.e_value < best {
            hsp.e_value
        } else {
            best
        }
    });
    Some(best_evalue)
}

/// `BLAST_Nint`.
fn blast_nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

/// `BlastCompo_EarlyTermination`: every query's heap holds its hit list size and its worst
/// e-value is at most its cutoff, and `evalue` is above five times every cutoff.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1561-1579
fn early_termination(evalue: f64, heaps: &[BlastCompoHeap]) -> bool {
    heaps.iter().all(|heap| {
        heap.n >= heap.heap_threshold
            && heap.worst_evalue <= heap.ecutoff
            && evalue > COMPO_HEAP_EVALUE_STRETCH * heap.ecutoff
    })
}

/// `s_preliminaryTestNearIdentical` with a cutoff above 0 (the redo's).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1063-1076
fn preliminary_test_near_identical(query_length: i32, align: &Align, cutoff: f64) -> bool {
    if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL) {
        return false;
    }
    let align_len = (align.query_end - align.query_start).min(align.match_end - align.match_start);
    // A float division as NCBI's (0 gives an infinity or NaN, neither below the cutoff).
    if (align.score as f64) / (align_len as f64) < cutoff {
        return false;
    }
    true
}

/// `s_preliminaryTestNearIdentical` with the cutoff 0 (`RedoRun::without_valid_context`): the
/// window `hspcnt` holds one alignment, without gaps, at least as long as the profile or 50.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1078-1094
fn old_test_near_identical(query_length: i32, hspcnt: usize, align: &Align) -> bool {
    hspcnt == 1
        && align.query_end - align.query_start == align.match_end - align.match_start
        && align.match_end - align.match_start + 1
            >= query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
}

/// `KAPPA_CONTAINED_IN_HSP` (redo_alignment.c:72-73) and `CONTAINED_IN_HSP` (blast_kappa.c).
fn contained_in_hsp(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    (a <= c && b >= c) && (d <= f && e >= f)
}

/// `s_IsContained`: whether an alignment of the list holds both ends of `in_align` and scores
/// at least `in_align`'s score plus 2 bits.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:966-995
///
/// `C`, the redo's instance, keeps each instance's own copy (psiblast's code as before the frames'
/// instance; DECISIONS 63 S3).
fn is_contained<C: RedoInstance>(in_align: &Align, alignments: &[Align], lambda: f64) -> bool {
    let score_thresh = in_align.score as f64 + KAPPA_BIT_TOL * LOCAL_LN2 / lambda;
    alignments.iter().any(|align| {
        in_align.frame.signum() == align.frame.signum()
            && contained_in_hsp(
                align.query_start,
                align.query_end,
                in_align.query_start,
                align.match_start,
                align.match_end,
                in_align.match_start,
            )
            && contained_in_hsp(
                align.query_start,
                align.query_end,
                in_align.query_end,
                align.match_start,
                align.match_end,
                in_align.match_end,
            )
            && score_thresh <= align.score as f64
    })
}

/// `s_WithDistinctEnds`: `new_align` joins the list (at its head, without the alignments that
/// share an end with it) unless an alignment that shares an end (with the same adjustment) or a
/// similar end (otherwise) scores at least as much.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:400-449
///
/// `C`, the redo's instance, keeps each instance's own copy (psiblast's code as before the frames'
/// instance; DECISIONS 63 S3).
fn with_distinct_ends<C: RedoInstance>(
    new_align: Align,
    alignments: &mut Vec<Align>,
    is_same_adjustment: bool,
) {
    let excluded = alignments.iter().any(|align| {
        align.frame == new_align.frame
            && if is_same_adjustment {
                is_same_end_point(&new_align, align)
            } else {
                is_similar_end_point(&new_align, align)
            }
            && new_align.score <= align.score
    });
    if excluded {
        return;
    }
    alignments
        .retain(|align| !(align.frame == new_align.frame && is_same_end_point(&new_align, align)));
    alignments.insert(0, new_align);
}

/// `s_IsSameEndPoint`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:305-314
fn is_same_end_point(new_align: &Align, align: &Align) -> bool {
    (align.query_start == new_align.query_start && align.match_start == new_align.match_start)
        || (align.query_end == new_align.query_end && align.match_end == new_align.match_end)
}

/// `s_IsSimilarEndPoint`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:321-348
fn is_similar_end_point(new_align: &Align, align: &Align) -> bool {
    let start_contained = contained_in_hsp(
        align.query_start,
        align.query_end,
        new_align.query_start,
        align.match_start,
        align.match_end,
        new_align.match_start,
    );
    let end_contained = contained_in_hsp(
        align.query_start,
        align.query_end,
        new_align.query_end,
        align.match_start,
        align.match_end,
        new_align.match_end,
    );
    (start_contained
        && new_align.query_start - new_align.match_start == align.query_start - align.match_start)
        || (end_contained
            && new_align.query_end - new_align.match_end == align.query_end - align.match_end)
}

/// `ScoreCompareHSPs`: the score, highest first, then the subject start, the subject end (last
/// first), the query start and the query end (last first).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
fn score_compare(a: &Align, b: &Align) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.match_start.cmp(&b.match_start))
        .then_with(|| b.match_end.cmp(&a.match_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

/// `s_HitlistReapContained`: an HSP whose ends lie in an earlier HSP of the same frames and
/// whose score is at most that HSP's is deleted.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:232-262
/// Every HSP of a PSSM search has the query frame 0, and the subject frame 0 unless the subject
/// is translated (`RedoInstance::SUBJECT_IS_TRANSLATED`).
fn reap_contained<C: RedoInstance>(redone: &mut Vec<Redone>) {
    let mut keep = vec![true; redone.len()];
    for iread in 1..redone.len() {
        let hsp1 = &redone[iread].align;
        for ireadback in 0..iread {
            if !keep[ireadback] {
                continue;
            }
            let hsp2 = &redone[ireadback].align;
            if (!C::SUBJECT_IS_TRANSLATED || hsp2.frame == hsp1.frame)
                && contained_in_hsp(
                    hsp2.query_start,
                    hsp2.query_end,
                    hsp1.query_start,
                    hsp2.match_start,
                    hsp2.match_end,
                    hsp1.match_start,
                )
                && contained_in_hsp(
                    hsp2.query_start,
                    hsp2.query_end,
                    hsp1.query_end,
                    hsp2.match_start,
                    hsp2.match_end,
                    hsp1.match_end,
                )
                && hsp1.score <= hsp2.score
            {
                keep[iread] = false;
                break;
            }
        }
    }
    let mut index = 0;
    redone.retain(|_| {
        index += 1;
        keep[index - 1]
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039 (cited in the
// module): `-use_sw_tback` (`eSmithWatermanTbck`) redoes each match with
// `Blast_RedoOneMatchSmithWaterman` in the place of `Blast_RedoOneMatch` (series PR, item
// PS11c-3: RPS-BLAST's redo).
pub(crate) mod smith_waterman;

#[cfg(test)]
mod frame_tests;

#[cfg(test)]
#[path = "redo1_tests.rs"]
mod tests;
