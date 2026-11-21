//! BLASTX six-context Kappa redo adapter. No BLASTP searches or external runtime calls.
use super::{
    args::ResolvedOptions, linking::link_parameters, parameters::ContextParameters,
    preliminary::PreliminaryHsp, query_setup::PreparedQueryBatch,
};
use crate::blastinput::app::ErrorAfterBatchReports;
use crate::{
    algorithm::blastp::gapalign::{
        blast_gapped_alignment_with_traceback_with_scratch, GapAlignScratch,
    },
    config::{ProteinScoringSpec, ScoringMatrix},
    core::composition_adjustment::{
        adjust_scores::{
            build_matrix_info, AdjustedProteinMatrix, BlastAminoAcidComposition,
            BlastCompositionWorkspace,
        },
        redo_alignment::*,
        smith_waterman_match::{
            redo_one_match_smith_waterman, sw_batch_failed_error, sw_endless_error,
            unadjusted_square_matrix, SwMatchScratch,
        },
    },
    stats::{compute_blosum62_ideal_karlin_params, lookup_protein_params},
    utils::{
        glibc_fma_libm::GlibcLibm,
        matrix::BLASTAA_SIZE,
        seg::{SegMasker, SegParams},
    },
};
use anyhow::{bail, ensure, Context, Result};
use std::{cell::Cell, ffi::c_void, ptr::NonNull};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2309-2336
pub fn query_infos(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
) -> Vec<BlastCompoQueryInfo> {
    batch
        .contexts
        .iter()
        .zip(parameters)
        .map(|(c, p)| {
            let query = &batch.sequence_start[1 + c.offset..1 + c.offset + c.length];
            BlastCompoQueryInfo {
                origin: c.offset as i32,
                seq: BlastCompoSequenceData::from_ncbistdaa(query),
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
                composition: BlastAminoAcidComposition::empty(),
                eff_search_space: p.search_space as f64,
                // EXPERIMENT (LOSAT_X_BXENDMEM): no table; `get_range` reads the context's words.
                words: (!x_end_mem()).then(|| build_query_word_hashes(query)),
            }
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3109
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2438-2479
/// The composition mode of `-comp_based_stats` (`args.rs` `composition_mode`; 0 only with
/// `-use_sw_tback`, whose redo runs without composition-based statistics too).
fn compo_adjust_mode(composition: u8) -> BlastCompoAdjustMode {
    match composition {
        0 => BlastCompoAdjustMode::NoCompositionBasedStats,
        1 => BlastCompoAdjustMode::CompositionBasedStats,
        3 => BlastCompoAdjustMode::ForceFullMatrixAdjust,
        _ => BlastCompoAdjustMode::CompositionMatrixAdjust,
    }
}

/// The composition workspace of the search's matrix (`NRrecord`, `Blast_CompositionWorkspaceInit`
/// with `scoringParams->options->matrix`): BLOSUM62's as before, another matrix's made from its
/// frequency data (`BlastCompositionWorkspace::new`, cold, as TBLASTN's `use_matrix_workspace`).
fn composition_workspace(matrix: ScoringMatrix) -> Result<BlastCompositionWorkspace> {
    if matrix == ScoringMatrix::Blosum62 {
        Ok(BlastCompositionWorkspace::new_blosum62())
    } else {
        other_matrix_workspace(matrix)
    }
}

#[cold]
#[inline(never)]
fn other_matrix_workspace(matrix: ScoringMatrix) -> Result<BlastCompositionWorkspace> {
    BlastCompositionWorkspace::new(matrix)
}

/// The redo's scaling factor: 32 with composition-based statistics, 1 without (the
/// Smith-Waterman redo of `-comp_based_stats 0`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3100-3108
pub(crate) fn redo_scale(options: &ResolvedOptions) -> f64 {
    if options.composition != 0 {
        32.0
    } else {
        1.0
    }
}

pub fn redo_params(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
) -> Result<BlastRedoAlignParams> {
    ensure!(
        (options.composition != 0 || options.use_sw_tback) && options.gapped,
        "BLASTX Kappa requires a gapped composition mode or the Smith-Waterman traceback"
    );
    let scale = redo_scale(options);
    let matrix = options.scoring;
    let ka = lookup_protein_params(&ProteinScoringSpec {
        matrix,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    let scaled_lambda = ka.lambda / scale;
    // The ideal block of the matrix: BLOSUM62's as before, another matrix's as TBLASTN's redo
    // takes it (`stage_d_kappa_params.rs` `local_kappa_redo_params`).
    let ideal = if matrix == ScoringMatrix::Blosum62 {
        compute_blosum62_ideal_karlin_params().map_err(anyhow::Error::msg)?
    } else {
        crate::algorithm::tblastx::lookup::ideal_karlin_params_of(matrix)
    };
    let do_link = link_parameters(options).is_some();
    // `hit_params->cutoff_score_min` is that of the query batch: the reader batch when the unit
    // is a part of one (`reader_batch.rs`).
    let cutoff = match batch
        .reader_shape
        .as_deref()
        .and_then(|shape| shape.redo_cutoff_min)
    {
        Some(batch_min) => batch_min,
        None => parameters
            .iter()
            .filter(|p| p.valid)
            .map(|p| p.hit_cutoff)
            .min()
            .context("BLASTX Kappa requires a valid context")?,
    };
    // blast_parameters.c:454-463: the search's final X-drop (`parameters::gapped_x_drops`).
    let extension_final = super::parameters::gapped_x_drops(options, &ka).1;
    Ok(BlastRedoAlignParams {
        matrix_info: build_matrix_info(matrix, ideal.lambda / scale)?,
        gapping_params: BlastCompoGappingParams {
            gap_open: options.gap_open * scale as i32,
            gap_extend: options.gap_extend * scale as i32,
            // blast_kappa.c:2374-2382 (`core::x_drop::redo_x_dropoff`).
            x_dropoff: crate::core::x_drop::redo_x_dropoff(
                options.gap_x_dropoff_final,
                scaled_lambda,
                extension_final,
            ),
            context: Cell::new(None),
        },
        compo_adjust_mode: compo_adjust_mode(options.composition),
        alphsize: BLASTAA_SIZE as i32,
        composition_test_index: 0,
        log_k: 0.0,
        score_divisor: scale,
        // NCBI's traceback runs this redo with composition-based statistics or the Smith-Waterman
        // traceback, which redoes each match with or without them.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1487
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3038-3039
        smith_waterman: options.use_sw_tback,
        near_identical_cutoff: 1.74 * std::f64::consts::LN_2 / scaled_lambda,
        position_based: false,
        ccat_query_length: batch
            .contexts
            .iter()
            .map(|c| c.length as i32)
            .max()
            .unwrap_or(0),
        query_is_translated: true,
        subject_is_translated: false,
        cutoff_score: if do_link {
            (cutoff as f64 * scale) as i32
        } else {
            1
        },
        do_link_hsps: do_link,
    })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1896-1913
struct RedoContext<'a> {
    input: &'a [PreliminaryHsp],
    scratch: NonNull<GapAlignScratch>,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1677-1709
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1646
fn get_range<'a>(
    matching: &BlastCompoMatchingSequence<'a>,
    subject_range: &BlastCompoSequenceRange,
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    words: Option<&[u64]>,
    align: &BlastCompoAlignment,
    should_test: bool,
    biased: bool,
    params: &BlastRedoAlignParams,
) -> Result<BlastRedoRangeResult> {
    let query =
        BlastCompoSequenceData::copy_query_range_with_selenocysteine_fix(orig_query, query_range);
    let mut subject = BlastCompoSequenceData::from_ncbistdaa(matching.data);
    let mut biased = biased;
    if params.uses_composition_based_stats()
        && biased
        && (!should_test
            || !match words {
                Some(_) => {
                    test_near_identical(&subject, 0, &query, query_range.begin, words, align)
                }
                None => x_test_near_identical_on_residues(
                    &subject,
                    0,
                    &query,
                    query_range.begin,
                    orig_query.data(),
                    align,
                ),
            })
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
        // The intervals come from a per-thread memo of the last subject's SEG (byte-equal input,
        // the same masker call on a miss; x_seg_memo.rs).
        let intervals = super::x_seg_memo::x_subject_seg(subject.data(), |data| {
            SegMasker::with_params(&SegParams::new(10, 1.8, 2.1)).mask_sequence(data)
        });
        biased = !intervals.is_empty();
        for interval in intervals {
            for residue in &mut subject.buffer[1 + interval.start..1 + interval.end] {
                *residue = 21;
            }
        }
    }
    subject.fit_protein_range_in_place(subject_range);
    Ok(BlastRedoRangeResult {
        query,
        subject,
        subject_maybe_biased: biased,
        x_subject_composition: None,
    })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1896-1957
fn redo_alignment(
    incoming: &BlastCompoAlignment,
    rule: EMatrixAdjustRule,
    adjusted: Option<&AdjustedProteinMatrix>,
    query: &BlastCompoSequenceData,
    qrange: &BlastCompoSequenceRange,
    _full_query: i32,
    subject: &BlastCompoSequenceData,
    srange: &BlastCompoSequenceRange,
    _full_subject: i32,
    params: &BlastRedoAlignParams,
) -> Result<Option<Box<BlastCompoAlignment>>> {
    let Some(BlastCompoAlignmentContext::PreliminaryHspIndex(index)) = incoming.context.as_ref()
    else {
        bail!("BLASTX Kappa incoming HSP index missing")
    };
    let ptr = params
        .gapping_params
        .context
        .get()
        .context("BLASTX Kappa context missing")?;
    // Safety: redo_context installs this stack-owned state only for the synchronous redo call and restores the pointer afterward.
    let context = unsafe { &*ptr.cast::<RedoContext<'_>>().as_ptr() };
    let h = &context.input[*index];
    // Safety: the calling function retains exclusive ownership of scratch throughout the synchronous callback.
    let scratch = unsafe { &mut *context.scratch.as_ptr() };
    let mut fence = false;
    let Some(mut result) = blast_gapped_alignment_with_traceback_with_scratch(
        query.data(),
        subject.data(),
        usize::try_from(h.q_gapped_start - qrange.begin)?,
        usize::try_from(h.s_gapped_start - srange.begin)?,
        params.matrix_info.matrix,
        adjusted,
        params.gapping_params.gap_open,
        params.gapping_params.gap_extend,
        params.gapping_params.x_dropoff,
        scratch,
        Some(&mut fence),
    ) else {
        return Ok(None);
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760-1770
    Ok(Some(blast_compo_alignment_new(
        result.score,
        rule,
        result.query_start + qrange.begin,
        result.query_stop + qrange.begin,
        qrange.context,
        result.subject_start + srange.begin,
        result.subject_stop + srange.begin,
        srange.context,
        Some(BlastCompoAlignmentContext::EditScript(std::mem::take(
            &mut result.edit_script,
        ))),
    )))
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2389-2394
const CALLBACKS: BlastRedoAlignCallbacks = BlastRedoAlignCallbacks {
    calc_lambda: Some(redo_calc_lambda),
    get_range,
    redo_one_alignment: redo_alignment,
};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:770-803
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3590-3641
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3641
/// The redo's list of a match's preliminary HSPs in `context_index`, in their order, with
/// their scores scaled by `score_divisor` (each names its HSP, for the X-drop redo's callback).
fn incoming_alignments(
    input: &[PreliminaryHsp],
    context_index: usize,
    score_divisor: f64,
) -> Option<Box<BlastCompoAlignment>> {
    let mut incoming = None;
    let mut tail = &mut incoming;
    for (index, h) in input
        .iter()
        .enumerate()
        .filter(|(_, h)| h.context == context_index)
    {
        let node = blast_compo_alignment_new(
            (h.score as f64 * score_divisor) as i32,
            EMatrixAdjustRule::DontAdjustMatrix,
            h.q_start,
            h.q_end,
            h.context as i32,
            h.s_start,
            h.s_end,
            0,
            Some(BlastCompoAlignmentContext::PreliminaryHspIndex(index)),
        );
        *tail = Some(node);
        tail = &mut tail.as_mut().expect("incoming inserted").next;
    }
    incoming
}

/// The Smith-Waterman redo of the match's preliminary HSPs in `context_index` (`-use_sw_tback`,
/// item X4b-2): core's match body (`redo_one_match_smith_waterman`) on the context's windows (the
/// HSPs' query range with the border, the whole subject), with BLASTX's subject range
/// (`get_range`). `lambda` and `sw.log_k` are the scaled gapped block's, the context's
/// (`kbp_gap[context_index]`); `would_insert(context, evalue, score)` is the heap test of e-value
/// mode, which NCBI makes on the heap the window's context indexes (`kappa_batch`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3601-3622
#[allow(clippy::too_many_arguments)]
fn sw_redo_context(
    input: &[PreliminaryHsp],
    context_index: usize,
    infos: &[BlastCompoQueryInfo],
    subject: &[u8],
    params: &BlastRedoAlignParams,
    lambda: f64,
    sw: &mut SwState,
    workspace: &mut BlastCompositionWorkspace,
    matrix_state: &mut Option<AdjustedProteinMatrix>,
    would_insert: &mut dyn FnMut(usize, f64, i32) -> Result<bool>,
) -> Result<BlastRedoOneMatchResult> {
    let incoming = incoming_alignments(input, context_index, params.score_divisor);
    let matching = BlastCompoMatchingSequence::new(0, subject);
    redo_one_match_smith_waterman(
        &incoming,
        params,
        lambda,
        sw.log_k,
        &matching,
        infos,
        &CALLBACKS,
        workspace,
        matrix_state,
        sw.unadjusted.as_ref(),
        sw.cutoff_e,
        would_insert,
        &mut sw.scratch,
    )
}

/// The error of a query batch whose last redone match failed its Smith-Waterman redo: NCBI's
/// `Traceback failed` (exit 3), or LOSAT's refusal where NCBI's clean-up crashes, when the
/// batch's first query holds a redone match (`first_query_holds_match`, the heap of the batch's
/// first query is not empty; PD-LOSAT-NCBI-DEFECTS 1.53, 1.68 and 1.69, rule 3), after the reports
/// of the batches before (`ErrorAfterBatchReports`, which `native::search_error_as_is` passes).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3775-3779
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3795-3802
#[cold]
#[inline(never)]
pub(crate) fn sw_batch_failed(first_query_holds_match: bool) -> anyhow::Error {
    ErrorAfterBatchReports(sw_batch_failed_error("blastx", first_query_holds_match)).into()
}

/// LOSAT's refusal of a batch with a match whose Smith-Waterman redo does not end in NCBI (a
/// significant score of 0 in a window of two or more preliminary HSPs whose first query letter,
/// the window's first, scores the subject's first letter above 0; PD-LOSAT-NCBI-DEFECTS 1.69, in
/// the class of 1.30, rule 3), after the reports of the batches before. What makes the score of
/// 0 significant: a cutoff score of 0 with sum statistics (`do_link_hsps`), or its e-value in
/// e-value mode.
#[cold]
#[inline(never)]
pub(crate) fn sw_endless(do_link_hsps: bool) -> anyhow::Error {
    let why = if do_link_hsps {
        "a score of 0 reaches the cutoff score, a window of the match has several preliminary HSPs and the window's first query letter scores the subject's first letter above 0"
    } else {
        "the e-value of a score of 0 passes the threshold, a window of the match has several preliminary HSPs and the window's first query letter scores the subject's first letter above 0"
    };
    ErrorAfterBatchReports(sw_endless_error("blastx", why)).into()
}

/// LOSAT's refusal of `-use_sw_tback` with `-ungapped` (which only `-comp_based_stats 0` lets
/// reach a search): NCBI's traceback runs the Smith-Waterman redo for an ungapped search too, and
/// its set-up reads the gapped Karlin blocks of every context, which an ungapped search does not
/// make, and crashes (SIGSEGV) at every query batch NCBI searches, with hits or without, with
/// `-subject` and with `-db` (line RG X4b's survey, a debugger's trace; probes `u.*` of item X4b-2
/// on ws, 2026-11-21). PD-LOSAT-NCBI-DEFECTS 1.69, rule 3. The caller refuses at the traceback of
/// the first query batch that is searched, after the reports of the batches before it: batches
/// without a valid context (a two-base query alone), which NCBI reports with a warning and does
/// not search (rows `db.x4b2.rej.ungapped.ac.b1.*`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1498
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3092-3095
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2073-2074
#[cold]
#[inline(never)]
pub(crate) fn sw_ungapped_refusal() -> anyhow::Error {
    ErrorAfterBatchReports(anyhow::anyhow!(
        "-use_sw_tback with -ungapped, on which NCBI BLAST+ 2.17.0's blastx crashes at every query batch it searches, is not supported by LOSAT's BLASTX"
    ))
    .into()
}

/// The heap test of a redo that makes none: the X-drop redo, and the Smith-Waterman redo in
/// cutoff-score mode (sum statistics), whose significance reads no heap.
pub(crate) fn no_heap_test(_context: usize, _evalue: f64, _score: i32) -> Result<bool> {
    bail!("the redo's heap test was not expected here")
}

pub(crate) fn redo_context_observed(
    input: &[PreliminaryHsp],
    context_index: usize,
    infos: &[BlastCompoQueryInfo],
    subject: &[u8],
    params: &BlastRedoAlignParams,
    lambda: f64,
    scratch: &mut GapAlignScratch,
    workspace: &mut BlastCompositionWorkspace,
    matrix_state: &mut Option<AdjustedProteinMatrix>,
) -> Result<BlastRedoOneMatchResult> {
    let incoming = incoming_alignments(input, context_index, params.score_divisor);
    let matching = BlastCompoMatchingSequence::new(0, subject);
    let state = RedoContext {
        input,
        scratch: NonNull::from(scratch),
    };
    let old = params
        .gapping_params
        .context
        .replace(Some(NonNull::from(&state).cast::<c_void>()));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3641
    let result = blast_redo_one_match_with_workspace_queries_and_matrix(
        &incoming,
        params,
        &matching,
        infos,
        lambda,
        &CALLBACKS,
        workspace,
        matrix_state,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3641-3647
    params.gapping_params.context.set(old);
    result
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:305-356
#[derive(Clone, Debug)]
pub struct RedoneHsp {
    pub hsp: PreliminaryHsp,
    pub edit_script: Vec<crate::common::GapEditOp>,
    pub composition_method: i32,
    pub identity_matrix: Option<std::sync::Arc<AdjustedProteinMatrix>>,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:223-278
// BLASTX's untranslated protein subject has frame 0 throughout this owner;
// the NCBI subject-frame equality predicate is therefore always true here.
fn reap_contained_hsps(output: Vec<RedoneHsp>) -> Vec<RedoneHsp> {
    let mut slots: Vec<Option<RedoneHsp>> = output.into_iter().map(Some).collect();
    for index in 1..slots.len() {
        let h = &slots[index].as_ref().expect("live candidate HSP").hsp;
        let mut contained = false;
        for other in &slots[..index] {
            let Some(other) = other else { continue };
            let p = &other.hsp;
            let result = p.frame == h.frame
                && p.q_start <= h.q_start
                && h.q_start <= p.q_end
                && p.s_start <= h.s_start
                && h.s_start <= p.s_end
                && p.q_start <= h.q_end
                && h.q_end <= p.q_end
                && p.s_start <= h.s_end
                && h.s_end <= p.s_end
                && h.score <= p.score;
            if result {
                contained = true;
                break;
            }
        }
        if contained {
            slots[index] = None;
        }
    }
    slots.into_iter().flatten().collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3340-3367
pub struct KappaState {
    params: BlastRedoAlignParams,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3504
    // NCBI builds one query_info per thread (query_info_tld). Here the query information is never
    // written during the redo, so the per-thread states share one copy. Placement only.
    // Read-only, shared with the per-thread copies.
    infos: std::sync::Arc<Vec<BlastCompoQueryInfo>>,
    lambda: f64,
    scratch: GapAlignScratch,
    workspace: BlastCompositionWorkspace,
    matrix: Option<AdjustedProteinMatrix>,
    /// `-use_sw_tback`: what the Smith-Waterman redo reads besides the X-drop redo's state
    /// (`SwState`); `None` without it.
    sw: Option<Box<SwState>>,
}

/// What the Smith-Waterman redo of a BLASTX match (`-use_sw_tback`, item X4b-2) reads besides
/// the X-drop redo's state: the scaled block's `logK` (the significance test's e-value), `-evalue`
/// (`params->cutoff_e`), the search's matrix that the windows read without composition-based
/// statistics (`unadjusted_square_matrix`, `None` with them) and the round loop's memory.
struct SwState {
    log_k: f64,
    cutoff_e: f64,
    unadjusted: Option<AdjustedProteinMatrix>,
    scratch: SwMatchScratch,
}

impl SwState {
    /// `None` without `-use_sw_tback`. `log_k` is the logarithm of the gapped block's K, which
    /// the redo's rescaling leaves unscaled (every context of a BLASTX query has the matrix's
    /// gapped block).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2134
    fn of(options: &ResolvedOptions, k: f64) -> Option<Box<Self>> {
        options.use_sw_tback.then(|| {
            Box::new(Self {
                log_k: k.glibc_ln(),
                cutoff_e: options.evalue,
                unadjusted: (options.composition == 0)
                    .then(|| unadjusted_square_matrix(options.scoring)),
                scratch: SwMatchScratch::new(),
            })
        })
    }
}

impl KappaState {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3244-3248
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3340-3367
    pub fn new(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
    ) -> Result<Self> {
        let params = redo_params(batch, parameters, options)?;
        let infos = query_infos(batch, parameters);
        let ka = lookup_protein_params(&ProteinScoringSpec {
            matrix: options.scoring,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        });
        Ok(Self {
            params,
            infos: std::sync::Arc::new(infos),
            lambda: ka.lambda / redo_scale(options),
            scratch: GapAlignScratch::new(),
            workspace: composition_workspace(options.scoring)?,
            matrix: None,
            sw: SwState::of(options, ka.k),
        })
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3326-3338,3346-3350,3499-3502
    // Each NCBI thread owns an NRrecord (composition workspace), a redo_align_params and a matrix.
    // This builds the same per-thread objects for one more Rust thread; the query information is
    // shared because it is read-only.
    /// A second state for another thread, as NCBI c++/src/algo/blast/core/blast_kappa.c:3244-3340
    /// allocates per thread (`redo_align_params_tld`, `gap_align_tld`, `NRrecord_tld`); the query
    /// information is read-only and shared.
    #[cfg(feature = "parallel")]
    pub(crate) fn x_for_thread(
        &self,
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
    ) -> Result<Self> {
        Ok(Self {
            params: redo_params(batch, parameters, options)?,
            infos: std::sync::Arc::clone(&self.infos),
            lambda: self.lambda,
            scratch: GapAlignScratch::new(),
            workspace: composition_workspace(options.scoring)?,
            matrix: None,
            sw: self.sw.as_ref().map(|sw| {
                Box::new(SwState {
                    log_k: sw.log_k,
                    cutoff_e: sw.cutoff_e,
                    unadjusted: sw.unadjusted.clone(),
                    scratch: SwMatchScratch::new(),
                })
            }),
        })
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3502,3625-3641
    // In NCBI the matrix is thread local (matrix_tld[tid]) and is passed on from one
    // Blast_RedoOneMatch call to the next on that thread. x_take_matrix and x_set_matrix move it
    // between Rust states so that a match redone out of order sees the matrix it would have seen in
    // stream order.
    /// Whether the redo's significance is the cutoff score of sum statistics (`do_link_hsps`), not
    /// the e-value.
    pub(crate) fn do_link_hsps(&self) -> bool {
        self.params.do_link_hsps
    }

    /// The matrix a match leaves for the next one.
    #[cfg(feature = "parallel")]
    pub(crate) fn x_take_matrix(&mut self) -> Option<AdjustedProteinMatrix> {
        self.matrix.take()
    }
    // No NCBI counterpart: see x_take_matrix above (the matrix hand-over is a placement of the same
    // value); it does not change any value NCBI computes.
    /// Hand a match the matrix of its predecessor.
    pub(crate) fn x_set_matrix(&mut self, matrix: Option<AdjustedProteinMatrix>) {
        self.matrix = matrix;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3625-3647
    pub fn redo_list_observed(
        &mut self,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        subject: &[u8],
    ) -> Result<Vec<RedoneHsp>> {
        self.redo_list_tested(input, batch, subject, &mut no_heap_test)
    }

    /// `redo_list_observed` with the Smith-Waterman redo's heap test (`sw_redo_context`): each
    /// frame of the match in turn, as NCBI's redo loop calls the redo once per frame; with
    /// `-use_sw_tback` an error of a frame's redo ends the match (NCBI goes to the match's
    /// clean-up and redoes no later frame).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3644-3646
    pub(crate) fn redo_list_tested(
        &mut self,
        input: &[PreliminaryHsp],
        batch: &PreparedQueryBatch,
        subject: &[u8],
        would_insert: &mut dyn FnMut(usize, f64, i32) -> Result<bool>,
    ) -> Result<Vec<RedoneHsp>> {
        let query_index = input.first().map(|h| h.context / 6);
        ensure!(
            input.iter().all(|h| Some(h.context / 6) == query_index),
            "BLASTX Kappa requires one query/OID HSP list"
        );
        let mut output = Vec::new();

        for context_index in query_index.unwrap_or(0) * 6..query_index.unwrap_or(0) * 6 + 6 {
            if !input.iter().any(|h| h.context == context_index) {
                continue;
            }
            let result = match self.sw.as_deref_mut() {
                None => redo_context_observed(
                    input,
                    context_index,
                    &self.infos,
                    subject,
                    &self.params,
                    self.lambda,
                    &mut self.scratch,
                    &mut self.workspace,
                    &mut self.matrix,
                ),
                Some(sw) => sw_redo_context(
                    input,
                    context_index,
                    &self.infos,
                    subject,
                    &self.params,
                    self.lambda,
                    sw,
                    &mut self.workspace,
                    &mut self.matrix,
                    would_insert,
                ),
            }?;
            for (context, alignments) in result.alignments_by_query.into_iter().enumerate() {
                let mut current = alignments;
                while let Some(mut align) = current {
                    current = align.next.take();
                    let Some(BlastCompoAlignmentContext::EditScript(script)) = align.context.take()
                    else {
                        bail!("BLASTX Kappa alignment edit script missing")
                    };
                    let method = match align.matrix_adjust_rule {
                        EMatrixAdjustRule::DontAdjustMatrix => 0,
                        EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                        _ => 2,
                    };
                    output.push(RedoneHsp {
                        hsp: PreliminaryHsp {
                            context,
                            frame: batch.contexts[context].frame,
                            score: align.score,
                            q_start: align.query_start,
                            q_end: align.query_end,
                            q_gapped_start: 0,
                            s_start: align.match_start,
                            s_end: align.match_end,
                            s_gapped_start: 0,
                        },
                        edit_script: script,
                        composition_method: method,
                        identity_matrix: None,
                    });
                }
            }
            output.sort_by(|a, b| super::preliminary::compare_score(&a.hsp, &b.hsp));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:238-266
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3689-3703
        let final_matrix = self
            .matrix
            .as_ref()
            .map(|matrix| std::sync::Arc::new(matrix.clone()));
        for h in &mut output {
            h.identity_matrix = final_matrix.clone();
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3670-3673
        if output.len() > 1 {
            Ok(reap_contained_hsps(output))
        } else {
            Ok(output)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:248-260,266-276
// No NCBI counterpart: reads the switch; it changes no value NCBI computes (`x_bx_mem.rs`, whose
// module is compiled only with threads, has the switch's comment).
/// EXPERIMENT (LOSAT_X_BXENDMEM): whether a batch's end makes no word-hash table and no copy of the
/// residues.
pub(crate) fn x_end_mem() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXENDMEM"))
}

#[cfg(test)]
#[path = "kappa_contained_reap_tests.rs"]
mod contained_reap_tests;
