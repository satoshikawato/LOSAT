//! The search of a run's queries against one profile volume: every (query, volume) pair alone,
//! each valid context of each query (a protein query's one, a translated query's frames) a unit
//! of the run's pool on a worker's scratch, the query's traceback run by the worker that searches
//! its last context, the results in query order (series PR, item PS12s, D2). The outfmt 0 K of
//! each query batch is computed from the queries' marks after (`batch_gapped_k`).
//!
//! NCBI searches the queries of a batch one after another on one diagonal table, then traces the
//! lists of the whole batch; a query's results depend on the volume, the options and the query
//! only (the table's entries of earlier queries and contexts read as empty, `wordfinder.rs`), so
//! LOSAT's bytes are the same at every thread count.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1161-1188
//! ```c
//!     for (index = 0; index < query_info->num_queries; ++index) {
//!         /* Separate one query from the set: create an auxiliary query_info
//!            structure which refers to this single query. */
//!         if (Blast_GetOneQueryStructs(&one_query_info, &one_query,
//!                                     query_info, query, index) != 0)
//!            return -1;
//!
//!         /* It is OK to pass NULL for the BlastDatabaseOptions argument, because it
//!            will not be checked for RPS BLAST program types. */
//!         status = (Int4)
//!           s_BlastSearchEngineCore(program_number, &concat_db, one_query_info,
//!              one_query, lookup_wrap, gap_align, score_params,
//!              word_params, ext_params, hit_params, NULL,
//!              diagnostics, aux_struct, &hsp_list, interrupt_search,
//!              progress_info);
//!         ...
//!         if (hsp_list && hsp_list->hspcnt > 0) {
//!             hsp_list->query_index = index;
//!             /* Save the HSP list */
//!             BlastHSPStreamWrite(hsp_stream, &hsp_list);
//!         }
//!     }
//! ```

use anyhow::Result;

use super::prelim::{RpsPrelim, RpsPrelimScratch};
use super::redo::RpsRedoParams;
use super::scan::{RpsBlockPool, RpsTable};
use super::traceback::{
    filter_results, formatter_identities, RpsTraceScratch, RpsTraceback, RpsTraced,
};
use super::wordfinder::{RpsWordParams, RpsWordScratch};
use super::{RpsRows, RPS_K_MULT};
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::pssm_search::evalue::EvalueStats;
use crate::algorithm::pssm_search::redo1::{RedoRun, SCALING_FACTOR};
use crate::seqdb::rps::RpsVolumeFiles;
use crate::stats::{BlastGumbelBlk, KarlinParams};

/// The options of an RPS-BLAST search that blastp's driver does not apply itself: `-max_hsps`
/// and `-qcov_hsp_perc` (the RPS traceback filters with the length of the first query of a
/// batch, `traceback::filter_results`), and the composition-based statistics mode that the
/// command line asks for (the mode the search runs is blastp's `comp_based_stats`, which
/// `rps_hooks::check_freq_ratio_file` sets to 0 where the database has no frequency ratios).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RpsSearch {
    pub(crate) max_hsps_per_subject: usize,
    pub(crate) query_cov_hsp_perc: f64,
    pub(crate) requested_cbs: u8,
    /// `-num_threads 0`: one thread per volume of the database (NCBI's "auto").
    pub(crate) auto_threads: bool,
    /// `-xdrop_ungap`, `-xdrop_gap` and `-xdrop_gap_final` in bits as given; `None`: NCBI's RPS
    /// defaults (7, 15 and 25 bits, which `rps_hooks::rps_search` reads).
    pub(crate) x_drop_ungapped_bits: Option<f64>,
    pub(crate) x_drop_gapped_bits: Option<f64>,
    pub(crate) x_drop_final_bits: Option<f64>,
    /// `-sum_stats true`, refused where the search starts (`rps_hooks::rps_search`).
    pub(crate) sum_stats: bool,
}

/// The contexts of a translated query.
pub(crate) const NUM_FRAMES: usize = 6;

impl RpsSearch {
    /// The X-drops in bits, `-xdrop_ungap`, `-xdrop_gap` and `-xdrop_gap_final` as given or the
    /// RPS options handle's defaults (7, 15 and 25 bits), and the gap trigger (22 bits): the
    /// `RpsRun` fields of both RPS programs.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:76-89
    /// ```c++
    /// CBlastRPSOptionsHandle::SetInitialWordOptionsDefaults()
    /// {
    ///     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    ///     ...
    /// CBlastRPSOptionsHandle::SetGappedExtensionDefaults()
    /// {
    ///     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    ///     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    ///     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
    /// ```
    /// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:122,128,137,144
    /// ```c
    /// #define BLAST_UNGAPPED_X_DROPOFF_PROT 7 /**< ungapped dropoff score for all
    /// #define BLAST_GAP_X_DROPOFF_PROT 15 /**< default dropoff (all protein-
    /// #define BLAST_GAP_TRIGGER_PROT 22.0 /**< default bit score that will trigger
    /// #define BLAST_GAP_X_DROPOFF_FINAL_PROT 25 /**< default dropoff (all protein-
    /// ```
    pub(crate) fn x_drop_bits(&self) -> RpsXDropBits {
        RpsXDropBits {
            ungapped: self.x_drop_ungapped_bits.unwrap_or(7.0),
            gapped: self.x_drop_gapped_bits.unwrap_or(15.0),
            final_: self.x_drop_final_bits.unwrap_or(25.0),
            gap_trigger: 22.0,
        }
    }
}

/// `RpsSearch::x_drop_bits`.
pub(crate) struct RpsXDropBits {
    pub(crate) ungapped: f64,
    pub(crate) gapped: f64,
    pub(crate) final_: f64,
    pub(crate) gap_trigger: f64,
}

/// One user query as the search reads it: its contexts, one for a protein query (rpsblast) and
/// the six frames of a nucleotide query (rpstblastn), which the reversed search scans as its
/// subjects.
pub(crate) struct RpsQuery<'a> {
    /// The contexts in NCBI's order (frame 0; or +1, +2, +3, -1, -2, -3).
    pub(crate) frames: Vec<RpsFrame<'a>>,
    /// The ungapped block of the first valid query of the query's batch: `kbp_psi[0]` of the
    /// score block that the composition-based redo copies, which stays where the profile has no
    /// score frequencies.
    pub(crate) batch_kbp: KarlinParams,
}

/// One context of a user query: a protein query's letters (frame 0) or one frame of a
/// nucleotide query (NCBI's `subject->frame`, `subject->sequence` and `subject->length` of the
/// frame loop).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-813
/// ```c
///     for (context=first_context; context<=last_context; context++) {
///         BlastHSPList* hsp_list_for_chunks = NULL;
///
///         if (kTranslatedSubject) {
///             Uint4 i;
///             subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
///             subject->sequence = translation_buffer + frame_offsets[context] + 1;
///             subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
///             if (subject->length > 0) stat_length = subject->length;
/// ```
#[derive(Debug, Clone, Copy)]
pub(crate) struct RpsFrame<'a> {
    /// 0 for a protein query; +1, +2, +3, -1, -2 or -3.
    pub(crate) frame: i8,
    /// The letters between two sentinels as searched (SEG, lower-case masks,
    /// `Blast_MaskUnsupportedAA`); a frame without letters is two sentinels.
    pub(crate) aa_seq: &'a [u8],
    /// The letters as read (between two sentinels) where the searched ones differ: the reports'
    /// identities count on them (`traceback::formatter_identities`).
    pub(crate) aa_seq_nomask: Option<&'a [u8]>,
    /// Whether its Karlin-Altschul parameters exist (an invalid context gets no HSP: its cutoffs
    /// are `INT4_MAX`, blast_parameters.c:916-919).
    pub(crate) valid: bool,
    /// `kbp_std[context]`: the context's composition with BLOSUM62.
    pub(crate) kbp_std: KarlinParams,
    /// `query_info->contexts[context].eff_searchsp` in this volume (NCBI's `Int8`): what the
    /// e-values and the hit cutoff read without the Gumbel block (`EvalueStats`), and the
    /// Smith-Waterman redo with it (series PR, item PS11c-3). 0 for an invalid context.
    pub(crate) search_space: i64,
    /// The length of the first query of the query's batch at this context (`-qcov_hsp_perc`).
    pub(crate) first_query_length: usize,
}

impl RpsFrame<'_> {
    /// The letters without the sentinels (`subject->length`).
    pub(crate) fn length(&self) -> usize {
        self.aa_seq.len().saturating_sub(2)
    }

    /// What the context's e-values read (`EvalueStats`): the Gumbel block, if any, and the
    /// context's search space.
    pub(super) fn stats<'g>(&self, gumbel: Option<&'g BlastGumbelBlk>) -> EvalueStats<'g> {
        EvalueStats::new(gumbel, self.search_space)
    }
}

impl RpsQuery<'_> {
    /// The context of `frame` (`BLAST_FrameToContext`; 0 for frame 0).
    pub(crate) fn context(&self, frame: i32) -> &RpsFrame<'_> {
        &self.frames[crate::core::blast_util::frame_to_context(frame)]
    }

    /// The length of the query's translation buffer (NCBI's `one_query->length`; a protein
    /// query's length): the last context's offset plus its length, where each context starts
    /// one sentinel after the previous one, an empty context taking no byte.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:72-77
    /// ```c++
    ///         Uint4 prev_loc = qinfo->contexts[index-1].query_offset;
    ///         Uint4 prev_len = qinfo->contexts[index-1].query_length;
    ///
    ///         Uint4 shift = prev_len ? prev_len + 1 : 0;
    ///
    ///         qinfo->contexts[index].query_offset = prev_loc + shift;
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:379-381
    /// ```c
    ///     one_query->length =
    ///         one_query_info->contexts[num_frames-1].query_offset +
    ///         one_query_info->contexts[num_frames-1].query_length;
    /// ```
    pub(crate) fn translation_length(&self) -> usize {
        let (last, earlier) = self.frames.split_last().expect("a query has a context");
        earlier
            .iter()
            .map(|frame| match frame.length() {
                0 => 0,
                length => length + 1,
            })
            .sum::<usize>()
            + last.length()
    }

    /// `stat_length` of the preliminary e-values: the length of the last context that has
    /// letters (a protein query's length; frame -3 with both strands, +3 with `-strand plus`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:813 (quoted at
    /// `RpsFrame`)
    pub(crate) fn stat_length(&self) -> usize {
        self.frames
            .iter()
            .rev()
            .map(RpsFrame::length)
            .find(|&length| length > 0)
            .unwrap_or(0)
    }
}

/// What the search of every query of a volume reads besides the volume.
pub(crate) struct RpsRun<'a> {
    pub(crate) queries: Vec<RpsQuery<'a>>,
    /// `kbp_gap_std`: BLOSUM62 11/1's table values.
    pub(crate) gapped: KarlinParams,
    /// `sbp->gbp`: none with `-searchsp` (at any value) or `OLD_FSC`.
    pub(crate) gumbel: Option<&'a BlastGumbelBlk>,
    /// `avg_subject_length` of the hit cutoff (`BLAST_SEQSRC_MINLENGTH`).
    pub(crate) cutoff_subject_length: i32,
    /// `-evalue`.
    pub(crate) evalue: f64,
    /// `hit_params->prelim_evalue` and the e-value of the hit cutoff: `-evalue`, five times it
    /// with composition-based statistics 2 and 3 and the Gumbel block (`cbs_stretch`,
    /// blast_parameters.c:937).
    pub(crate) prelim_evalue: f64,
    /// The composition-based statistics mode the search runs (0 to 3).
    pub(crate) cbs: u8,
    /// BLOSUM62's score block for the redo's matrix information (modes 1 to 3).
    pub(crate) sbp: Option<&'a crate::core::psi::PsiScoreBlk>,
    pub(crate) hitlist_size: usize,
    pub(crate) prelim_hitlist_size: usize,
    pub(crate) window: i32,
    /// The X-drops in bits (`word_options->x_dropoff`, `ext_options->gap_x_dropoff` and
    /// `gap_x_dropoff_final`: `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final` or NCBI's RPS
    /// defaults) and the gap trigger in bits (`kOptions->gap_trigger`).
    pub(crate) x_drop_ungapped_bits: f64,
    pub(crate) x_drop_gapped_bits: f64,
    pub(crate) x_drop_final_bits: f64,
    pub(crate) gap_trigger_bits: f64,
    pub(crate) options: &'a RpsSearch,
    /// The program of LOSAT's refusals (`RPSBLAST`, `RPSTBLASTN`).
    pub(crate) program: &'static str,
}

/// The refusal of a database whose scale makes a scaled score leave NCBI's `Int4` (gap costs,
/// X-drops, cutoffs): NCBI's products overflow there.
fn scale_overflow(scale: f64, program: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "a profile database of the scale {scale} (its scaled gap costs, X-drops or score cutoffs leave 32-bit integers) is not supported by LOSAT's {program}"
    )
}

/// The refusal of an X-drop in bits whose raw score (the bits over Lambda, converted to `Int4`,
/// times the database's scale where NCBI multiplies) leaves NCBI's `Int4`, or whose quotient by
/// the scaled gap extension makes NCBI's kernels' sizes leave it (`kernel_x_drop_fits`): NCBI's
/// conversion of the double, or its arithmetic, is undefined there (DECISIONS 17 N1/N2, line PS).
/// `given` is the option's value as the command line gives it; an X-drop of NCBI's defaults gets
/// the database's refusal (`scale_overflow`), as only an extreme scale makes it leave 32 bits.
fn x_drop_overflow(
    option: &str,
    given: Option<f64>,
    bits: f64,
    scale: f64,
    program: &str,
) -> anyhow::Error {
    if given.is_none() {
        return scale_overflow(scale, program);
    }
    anyhow::anyhow!(
        "an {option} of {bits:?} bits, whose raw X-drop at the profile database's scale {scale} leaves NCBI BLAST+'s 32-bit integers (its conversion or the sizes of its alignment kernels overflow: undefined behaviour), is not supported by LOSAT's {program}"
    )
}

// The limit of NCBI's gapped kernels' sizes (`kernel_x_drop_fits`, its `Int4` arithmetic on the X-drop
// over the scaled gap extension) is the protein programs' (`core::x_drop`). Reachable here only
// where the scaled gap extension is 1 or 2 (scale 1 or 2 with BLOSUM62 11/1): the quotient is at
// most `INT4_MAX / scale` otherwise.
use crate::core::x_drop::kernel_x_drop_fits;

/// A query's ungapped X-drop (`x_dropoff_init`) from the X-drop in bits `bits` (above 0,
/// `validate_protein_options`), the query's ungapped Lambda and the database's scale: the product
/// is a double, converted to `Int4` once; `None` where that conversion is undefined (`int4_of`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:216-221 (quoted at
/// `RpsCutoffs`)
fn ungapped_x_dropoff(bits: f64, lambda: f64, scale: f64) -> Option<i32> {
    int4_of(scale * (bits * std::f64::consts::LN_2 / lambda).ceil())
}

/// NCBI's `(Int4)` conversion of a double where it is defined: the value truncated toward zero,
/// `None` where that leaves `Int4` (or the value is not a number), which is undefined behaviour in
/// C (C11 6.3.1.4).
fn int4_of(value: f64) -> Option<i32> {
    let truncated = value.trunc();
    (truncated >= f64::from(i32::MIN) && truncated <= f64::from(i32::MAX))
        .then_some(truncated as i32)
}

/// The per-query cutoffs of the preliminary stage at the database's scale.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:216-221
/// ```c
///         p->cutoffs[context].x_dropoff_init =
///             (Int4)(sbp->scale_factor *
///                    ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:340-370
/// ```c
///       if (sbp->kbp_std) {     /* this may not be set for gapped blastn */
///          kbp = sbp->kbp_std[context];
///          if (s_BlastKarlinBlkIsValid(kbp)) {
///             gap_trigger = (Int4)((kOptions->gap_trigger * NCBIMATH_LN2 +
///                                      kbp->logK) / kbp->Lambda);
///          }
///       }
///       ...
///       } else {
///          new_cutoff = gap_trigger;
///       }
///       new_cutoff *= (Int4)sbp->scale_factor;
///       new_cutoff = MIN(new_cutoff,
///                        hit_params->cutoffs[context].cutoff_score_max);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:937-982
/// ```c
///              new_cutoff = BLAST_SpougeEtoS(cbs_stretch*evalue, kbp, sbp->gbp,
///                          query_info->contexts[context].query_length,
///                          avg_subject_length);
///          ...
///          params->cutoffs[context].cutoff_score = new_cutoff;
///          params->cutoffs[context].cutoff_score_max = new_cutoff;
///       ...
///             params->cutoffs[context].cutoff_score *= (Int4) scale_factor;
///             params->cutoffs[context].cutoff_score_max *= (Int4) scale_factor;
/// ```
struct RpsCutoffs {
    x_dropoff: i32,
    word_cutoff: i32,
    hit_cutoff: i32,
}

impl RpsRun<'_> {
    /// The query's cutoffs; the refusal where the ungapped X-drop (`-xdrop_ungap`) or a scaled
    /// cutoff leaves `i32`.
    fn cutoffs(&self, query: &RpsFrame<'_>, frames: usize, scale: f64) -> Result<RpsCutoffs> {
        let int_scale = scale as i32;
        let x_dropoff = ungapped_x_dropoff(self.x_drop_ungapped_bits, query.kbp_std.lambda, scale)
            .ok_or_else(|| {
                x_drop_overflow(
                    "-xdrop_ungap",
                    self.options.x_drop_ungapped_bits,
                    self.x_drop_ungapped_bits,
                    scale,
                    self.program,
                )
            })?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:926-931
        // ```c
        //          searchsp = query_info->contexts[context].eff_searchsp;
        //
        //          /* translated RPS searches must scale the search space down */
        //          /** @todo FIXME why only scale down rpstblastn search space?  */
        //          if (program_number == eBlastTypeRpsTblastn)
        //             searchsp /= NUM_FRAMES;
        // ```
        // (`frames` is 1 for rpsblast, `NUM_FRAMES` for rpstblastn.)
        let hit_cutoff = EvalueStats::new(self.gumbel, query.search_space / frames as i64)
            .hit_cutoff(
                self.prelim_evalue,
                &self.gapped,
                (query.aa_seq.len() - 2) as i32,
                self.cutoff_subject_length,
            )
            .checked_mul(int_scale)
            .ok_or_else(|| scale_overflow(scale, self.program))?;
        let gap_trigger = crate::core::x_drop::gap_trigger_raw_score(
            self.gap_trigger_bits,
            query.kbp_std.lambda,
            query.kbp_std.k,
        );
        Ok(RpsCutoffs {
            x_dropoff,
            word_cutoff: gap_trigger
                .checked_mul(int_scale)
                .ok_or_else(|| scale_overflow(scale, self.program))?
                .min(hit_cutoff),
            hit_cutoff,
        })
    }

    /// The gapped X-drops of the preliminary stage and of the final traceback at the database's
    /// scale (`ext_params->gap_x_dropoff`, `gap_x_dropoff_final`); the refusal where one of
    /// them, or the X-drop of the composition-based redo when it runs
    /// (`gapping_params->x_dropoff`), leaves `i32` (`-xdrop_gap`, `-xdrop_gap_final`). A gapped
    /// X-drop of 0 or below is NCBI's too: the alignment kernels raise it to the gap costs
    /// (`pssm_search::gapalign`, NCBI's `ALIGN_EX`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:455-469
    /// ```c
    ///    if (sbp->kbp_gap) {
    ///       double min_lambda = s_BlastFindSmallestLambda(sbp->kbp_gap, query_info, NULL);
    ///       params->gap_x_dropoff = (Int4)
    ///           (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    ///       /* Note that this conversion from bits to raw score is done prematurely for
    ///          RPS and PSSM searches, as the Lambda used here would change after
    ///          the PSSM is rescaled */
    ///       params->gap_x_dropoff_final = (Int4)
    ///           MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda, params->gap_x_dropoff);
    ///    }
    ///
    ///    if (sbp->scale_factor > 1.0) {
    ///        ASSERT(Blast_ProgramIsRpsBlast(program_number));
    ///        params->gap_x_dropoff *= (Int4)sbp->scale_factor;
    ///        params->gap_x_dropoff_final *= (Int4)sbp->scale_factor;
    ///    }
    /// ```
    /// (`min_lambda` is the Lambda of `kbp_gap_std`, BLOSUM62 11/1's, for every query; the
    /// products are `Int4` multiplications.) The redo's X-drop is `RedoRun::x_dropoff_of` of
    /// Lambda over the scaling factor (`s_RescaleSearch`) and the scaled final X-drop
    /// (blast_kappa.c:2375-2383, quoted at `RedoRun::new`).
    ///
    /// The kernels then size their memory from the X-drop over the gap extension, in `Int4`
    /// (`kernel_x_drop_fits`); `gap_open` and `gap_extend` are the scaled gap costs.
    fn gapped_x_drops(&self, scale: f64, gap_open: i32, gap_extend: i32) -> Result<(i32, i32)> {
        let int_scale = scale as i32;
        let lambda = self.gapped.lambda;
        let scaled = |raw: i32| {
            if scale > 1.0 {
                raw.checked_mul(int_scale)
            } else {
                Some(raw)
            }
        };
        // An X-drop that the user did not give leaves 32 bits only at an extreme scale: the
        // refusal then names the database's scale.
        let gapped_overflow = || {
            x_drop_overflow(
                "-xdrop_gap",
                self.options.x_drop_gapped_bits,
                self.x_drop_gapped_bits,
                scale,
                self.program,
            )
        };
        let final_overflow = || {
            x_drop_overflow(
                "-xdrop_gap_final",
                self.options.x_drop_final_bits,
                self.x_drop_final_bits,
                scale,
                self.program,
            )
        };
        let raw_gapped = int4_of(self.x_drop_gapped_bits * std::f64::consts::LN_2 / lambda)
            .ok_or_else(gapped_overflow)?;
        let raw_final = int4_of(
            (self.x_drop_final_bits * std::f64::consts::LN_2 / lambda).max(f64::from(raw_gapped)),
        )
        .ok_or_else(final_overflow)?;
        let x_drop_gapped = scaled(raw_gapped).ok_or_else(gapped_overflow)?;
        let x_drop_final = scaled(raw_final).ok_or_else(final_overflow)?;
        // The preliminary stage's kernel runs at every mode; the CBS 0 traceback's (the final
        // X-drop) only without the composition-based redo, whose own gap extension (the costs
        // times 32) keeps its quotient below 2^26.
        let redo = self.sbp.is_some() && self.cbs > 0;
        if !kernel_x_drop_fits(x_drop_gapped, gap_open, gap_extend) {
            return Err(gapped_overflow());
        }
        if !redo && !kernel_x_drop_fits(x_drop_final, gap_open, gap_extend) {
            return Err(final_overflow());
        }
        // The redo converts its X-drop at its own Lambda (`RedoRun::new`);
        // it is checked here, before any search.
        if redo {
            int4_of(RedoRun::x_dropoff_of(
                self.x_drop_final_bits,
                lambda / SCALING_FACTOR,
                x_drop_final,
            ))
            .ok_or_else(final_overflow)?;
        }
        Ok((x_drop_gapped, x_drop_final))
    }

    /// Every query's traced results on the volume `files`, whose consensus letters are
    /// `consensus` (at the profiles' offsets), in query order, on the scratch of the run's
    /// workers (`workers`): the volume's refusals before the queries' cutoffs (`Err`), else the
    /// queries' results (`RpsVolumeResults`), which a caller reads batch by batch in the order
    /// NCBI's search of a batch meets its errors (the queries search at the same time on several
    /// threads).
    pub(crate) fn search_volume(
        &self,
        files: &RpsVolumeFiles,
        consensus: &super::volumes::Consensus<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
        workers: &RpsWorkers,
    ) -> Result<RpsVolumeResults> {
        let rows = files.rows();
        let offsets = files.offsets();
        let scale = files.aux.scale;
        // (A whole number of 1 or more, `rps_hooks::rps_search`.)
        if scale > f64::from(i32::MAX) {
            return Err(scale_overflow(scale, self.program));
        }
        let int_scale = scale as i32;
        let table = RpsTable {
            cells: files.cells().expect("checked when the volume was opened"),
            overflow: files
                .overflow()
                .expect("checked when the volume was opened"),
            pv: files.pv(),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:558-562
        // ```c
        //    params->gap_open = score_options->gap_open * (Int4)scale_factor;
        //    params->gap_extend = score_options->gap_extend * (Int4)scale_factor;
        // ```
        let (Some(gap_open), Some(gap_extend)) = (
            files.aux.gap_open.checked_mul(int_scale),
            files.aux.gap_extend.checked_mul(int_scale),
        ) else {
            return Err(scale_overflow(scale, self.program));
        };
        // The X-drops at the scale (`gapped_x_drops`, NCBI's `BlastExtensionParametersNew`).
        let (x_drop_gapped, x_drop_final) = self.gapped_x_drops(scale, gap_open, gap_extend)?;
        // Every query's cutoffs, before the search (`BlastHitSavingParametersNew` and
        // `BlastInitialWordParametersNew` set them for the whole batch before its scan); a query
        // whose cutoffs are refused keeps that refusal apart (`RpsVolumeResults::batches` reads it
        // before any search error of its batch) and is not searched.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1372-1374
        // ```c
        //     if ((status = BlastInitialWordParametersNew(program_number,
        //             word_options, hit_params, lookup_wrap, sbp, query_info,
        //             BlastSeqSrcGetAvgSeqLen(seq_src), &word_params)) != 0)
        // ```
        let mut refused: Vec<Option<anyhow::Error>> = Vec::with_capacity(self.queries.len());
        let mut cutoffs: Vec<Vec<RpsCutoffs>> = Vec::with_capacity(self.queries.len());
        for query in &self.queries {
            let frames = query
                .frames
                .iter()
                .map(|frame| {
                    if frame.valid {
                        self.cutoffs(frame, query.frames.len(), scale)
                    } else {
                        Ok(RpsCutoffs {
                            x_dropoff: 0,
                            word_cutoff: 0,
                            hit_cutoff: 0,
                        })
                    }
                })
                .collect::<Result<Vec<_>>>();
            match frames {
                Ok(frames) => {
                    cutoffs.push(frames);
                    refused.push(None);
                }
                Err(refusal) => {
                    cutoffs.push(Vec::new());
                    refused.push(Some(refusal));
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1862-1895
        // ```c
        //       kbp[kbp_context]->Lambda /= scaling_factor;
        //       ...
        //       /* Put back the unscaled value of Lambda. */
        //       kbp[kbp_context]->Lambda *= scaling_factor;
        // ```
        // NCBI divides Lambda by the scale for each HSP's e-value and multiplies it back; at a
        // scale where that does not give Lambda back, NCBI's Lambda changes from HSP to HSP. A
        // refusal of the volume after the queries' cutoffs (NCBI's e-values come after its set-up),
        // so a batch's refusal of its cutoffs comes first; nothing is searched.
        let lambda = self.gapped.lambda;
        if lambda / scale * scale != lambda {
            let searched = refused.iter().map(|_| None).collect();
            return Ok(RpsVolumeResults {
                refused,
                late: Some(anyhow::anyhow!(
                    "a profile database of the scale {scale} (NCBI BLAST+ divides Lambda by the scale for each e-value and multiplies it back, which changes Lambda at this scale) is not supported by LOSAT's {}",
                    self.program
                )),
                searched,
            });
        }
        // Each searched query's result (`None` for a query refused at its cutoffs).
        let mut results: Vec<std::sync::Mutex<Option<Result<RpsTraced>>>> = self
            .queries
            .iter()
            .map(|_| std::sync::Mutex::new(None))
            .collect();
        let traceback = RpsTraceback {
            rows,
            offsets,
            consensus,
            alphabet: files.alphabet,
            scale,
            karlin_k: &files.aux.karlin_k,
            gap_open,
            gap_extend,
            x_drop_final,
            gapped: self.gapped,
            expect_value: self.evalue,
            hitlist_size: self.hitlist_size,
            redo: self.sbp.filter(|_| self.cbs > 0).map(|sbp| RpsRedoParams {
                files,
                sbp,
                gap_open: files.aux.gap_open,
                gap_extend: files.aux.gap_extend,
                x_drop_final_bits: self.x_drop_final_bits,
            }),
        };
        let prelim_karlin = KarlinParams {
            lambda: self.gapped.lambda / scale,
            ..self.gapped
        };
        // `ADAPTIVE_CBS`: the close of the volume's stream trims each query's preliminary hit
        // list (`BlastHSPCBSStreamClose`, blast_traceback.c:1839-1846; `core::blast_hspstream`).
        let cbs_close =
            crate::core::blast_hspstream::adaptive_cbs_close(self.cbs > 0, self.hitlist_size);
        // The preliminary stage of one context of a query (NCBI's frame loop: each context
        // scanned as a subject with its own cutoffs).
        let prelim_of = |q_idx: usize, context: usize| -> RpsPrelim<'_> {
            let frame = &self.queries[q_idx].frames[context];
            let cutoffs = &cutoffs[q_idx][context];
            RpsPrelim {
                word: RpsWordParams {
                    rows: RpsRows { rows },
                    table,
                    window: self.window,
                    cutoff: cutoffs.word_cutoff,
                    dropoff: cutoffs.x_dropoff,
                },
                rows,
                offsets,
                gap_open,
                gap_extend,
                x_drop_gapped,
                hit_cutoff: cutoffs.hit_cutoff,
                karlin: prelim_karlin,
                stats: frame.stats(self.gumbel),
                prelim_evalue: self.prelim_evalue,
                prelim_hitlist_size: self.prelim_hitlist_size,
                query_index: q_idx as u32,
                program: self.program,
            }
        };
        // The units of the volume's search (D2 of PS12s): each valid context of each query whose
        // cutoffs exist, query by query in the context order, each searched on its own
        // (`RpsPrelim::search_frame`) with its HSPs kept apart. A context's search reads nothing
        // another context writes: the diagonal table's entries of earlier contexts read as empty
        // (`wordfinder.rs`) and each context's scan starts its own runs of hits.
        let mut units: Vec<(usize, usize)> = Vec::new();
        let mut first_unit = Vec::with_capacity(self.queries.len());
        let mut remaining = Vec::with_capacity(self.queries.len());
        for (q_idx, query) in self.queries.iter().enumerate() {
            first_unit.push(units.len());
            let searched = refused[q_idx].is_none();
            if searched {
                units.extend(
                    query
                        .frames
                        .iter()
                        .enumerate()
                        .filter(|(_, frame)| frame.valid)
                        .map(|(context, _)| (q_idx, context)),
                );
            }
            let count = units.len() - first_unit[q_idx];
            if searched && count == 0 {
                *results[q_idx]
                    .get_mut()
                    .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(Ok(RpsTraced {
                    hit_list: None,
                    k_mark: None,
                }));
            }
            remaining.push(std::sync::atomic::AtomicUsize::new(count));
        }
        let searched_units: Vec<std::sync::Mutex<Option<Result<Vec<BlastpPreliminaryHsp>>>>> =
            units.iter().map(|_| std::sync::Mutex::new(None)).collect();
        // The rest of a query's search on this volume, run by the worker that searched its last
        // context: its contexts' HSPs concatenated in the context order (`Blast_HSPListAppend`
        // of NCBI's frame loop), the e-values, the reap and the collector (`finish`), then the
        // traceback. A context whose search failed gives the query its error, the first in the
        // context order (where NCBI's loop, and the search before PS12s, stopped).
        let finish_query = |scratch: &mut RpsScratch, q_idx: usize| -> Result<RpsTraced> {
            let query = &self.queries[q_idx];
            let none = RpsTraced {
                hit_list: None,
                k_mark: None,
            };
            let mut combined = std::mem::take(&mut scratch.prelim.combined);
            combined.clear();
            let mut failed = None;
            let mut last = 0;
            for unit in first_unit[q_idx]..first_unit[q_idx] + units_of(&first_unit, &units, q_idx)
            {
                let searched = searched_units[unit]
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner())
                    .take()
                    .expect("each context of a query is searched before its last one finishes");
                match searched {
                    Ok(hsps) if failed.is_none() => combined.extend(hsps),
                    Ok(_) => {}
                    Err(error) => {
                        failed.get_or_insert(error);
                    }
                }
                last = units[unit].1;
            }
            let finished = match failed {
                Some(error) => Err(error),
                None => Ok(prelim_of(q_idx, last).finish(&mut combined, query.stat_length())),
            };
            combined.clear();
            scratch.prelim.combined = combined;
            let Some(mut lists) = finished? else {
                return Ok(none);
            };
            if cbs_close {
                lists.cbs_close_trim(self.hitlist_size);
            }
            let mut traced = traceback.trace_query(
                lists,
                query,
                q_idx as u32,
                self.gumbel,
                &mut scratch.trace,
            )?;
            if let Some(hit_list) = traced.hit_list.as_mut() {
                if query
                    .frames
                    .iter()
                    .any(|frame| frame.aa_seq_nomask.is_some())
                {
                    formatter_identities(hit_list, query, consensus.letters()?, offsets);
                }
                finish_volume_list(
                    hit_list,
                    self.options.max_hsps_per_subject,
                    self.options.query_cov_hsp_perc,
                    query,
                    self.hitlist_size,
                );
            }
            Ok(traced)
        };
        // One unit: its context's HSPs, kept until the query's last context is searched; the
        // worker that searches that one finishes the query at once, so a query's preliminary
        // HSPs live only while the query is in flight.
        let run_unit = |scratch: &mut RpsScratch, unit: usize| {
            let (q_idx, context) = units[unit];
            let frame = &self.queries[q_idx].frames[context];
            // The worker's table and buckets for this volume (what the worker kept stays).
            scratch.prelim.word.prepare_for(rows.len(), self.window);
            let mut hsps = Vec::new();
            let searched = prelim_of(q_idx, context)
                .search_frame(
                    &mut scratch.prelim,
                    &frame.aa_seq[1..],
                    frame.frame,
                    &mut hsps,
                )
                .map(|()| hsps);
            *searched_units[unit]
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(searched);
            if remaining[q_idx].fetch_sub(1, std::sync::atomic::Ordering::AcqRel) == 1 {
                let traced = finish_query(scratch, q_idx);
                *results[q_idx]
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner()) = Some(traced);
            }
        };
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() && pool.threads() > 1 {
            use rayon::prelude::*;
            pool.install(|| {
                (0..units.len())
                    .into_par_iter()
                    .with_max_len(1)
                    .for_each(|unit| {
                        // The worker's index in the run's pool (its scratch, `RpsWorkers`).
                        let slot = rayon::current_thread_index().unwrap_or(0);
                        workers.with(slot, |scratch| run_unit(scratch, unit));
                    });
            });
            return Ok(take_results(refused, results));
        }
        let _ = pool;
        // Without a pool the units run on the caller's thread, on the first scratch.
        for unit in 0..units.len() {
            workers.with(0, |scratch| run_unit(scratch, unit));
        }
        Ok(take_results(refused, results))
    }
}

/// The number of units of the query `q_idx` (`search_volume`).
fn units_of(first_unit: &[usize], units: &[(usize, usize)], q_idx: usize) -> usize {
    first_unit.get(q_idx + 1).copied().unwrap_or(units.len()) - first_unit[q_idx]
}

/// A volume's search of a run's queries (`RpsRun::search_volume`), read query batch by query
/// batch (`batches`) in the order NCBI's search of one batch meets its errors.
pub(crate) struct RpsVolumeResults {
    /// Each query's refusal of its cutoffs (`None` where they exist).
    refused: Vec<Option<anyhow::Error>>,
    /// The volume's refusal after the queries' cutoffs (the Lambda at the volume's scale): no
    /// query was searched.
    late: Option<anyhow::Error>,
    /// Each searched query's result: its traced lists, or the first error of its search (`None`
    /// for a query that was not searched).
    searched: Vec<Option<Result<RpsTraced>>>,
}

impl RpsVolumeResults {
    /// The results of the query batches whose numbers of queries are `spans`, in query order:
    /// each batch's traced lists, up to the first batch that fails, whose error ends the list.
    /// A batch's error is the first refusal of its queries' cutoffs (NCBI sets the cutoffs of the
    /// whole batch before its scan, `search_volume`), else the volume's refusal after them, else
    /// the first error of its queries' searches in query order (the queries search at the same
    /// time).
    pub(crate) fn batches(self, spans: &[usize]) -> Vec<Result<Vec<RpsTraced>>> {
        let Self {
            refused,
            mut late,
            searched,
        } = self;
        let mut refused = refused.into_iter();
        let mut searched = searched.into_iter();
        let mut batches = Vec::with_capacity(spans.len());
        for &count in spans {
            let first_refusal = refused.by_ref().take(count).flatten().next();
            let results: Vec<Option<Result<RpsTraced>>> = searched.by_ref().take(count).collect();
            let batch = match first_refusal.or_else(|| late.take()) {
                Some(refusal) => Err(refusal),
                None => results
                    .into_iter()
                    .map(|result| result.expect("a query whose cutoffs exist is searched"))
                    .collect(),
            };
            let failed = batch.is_err();
            batches.push(batch);
            if failed {
                break;
            }
        }
        batches
    }

    /// The results of the run's queries as one batch (`batches`): rpsblast's.
    pub(crate) fn all(self) -> Result<Vec<RpsTraced>> {
        let count = self.refused.len();
        self.batches(&[count])
            .pop()
            .expect("one batch has its result")
    }
}

/// Every query's result of a volume's search, in query order (`RpsVolumeResults`).
fn take_results(
    refused: Vec<Option<anyhow::Error>>,
    results: Vec<std::sync::Mutex<Option<Result<RpsTraced>>>>,
) -> RpsVolumeResults {
    let searched = results
        .into_iter()
        .zip(&refused)
        .map(|(result, refusal)| {
            let result = result
                .into_inner()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            match refusal {
                Some(_) => None,
                None => Some(result.expect("every searched query of the volume has its result")),
            }
        })
        .collect();
    RpsVolumeResults {
        refused,
        late: None,
        searched,
    }
}

/// The scratch of one worker of an RPS search: its preliminary stage's (the word finder's
/// diagonal table and buckets, the gapped kernels' memory) and its traceback's.
struct RpsScratch {
    prelim: RpsPrelimScratch,
    trace: RpsTraceScratch,
}

impl RpsScratch {
    /// A worker's scratch, its buckets' blocks from the run's `pool`.
    fn new(pool: &std::sync::Arc<RpsBlockPool>) -> Self {
        Self {
            prelim: RpsPrelimScratch::new(RpsWordScratch::new(pool.clone())),
            trace: RpsTraceScratch::new(),
        }
    }
}

/// The scratch of a run's workers (series PR, item PS12s, D3): one for each thread of the run's
/// pool (`rayon::current_thread_index` inside the pool; the first for a search without a pool),
/// kept from volume to volume and from query batch to query batch, where each search of a volume
/// made its own (whose buckets grew again by doubling). A worker's diagonal table stays where the
/// next volume's table has the same length (`RpsWordScratch::prepare_for`); the buckets' blocks
/// come from one pool of the run (`RpsBlockPool`, M2), to which each bucket gives its blocks back
/// as soon as its hits are read. A search's values do not depend on what its scratch held before
/// (`wordfinder.rs`), so only the memory and its reuse change.
pub(crate) struct RpsWorkers {
    slots: Vec<std::sync::Mutex<RpsScratch>>,
    pool: std::sync::Arc<RpsBlockPool>,
}

impl RpsWorkers {
    /// The scratch of a run of `threads` threads (empty until a worker searches).
    pub(crate) fn new(threads: usize) -> Self {
        let pool = std::sync::Arc::new(RpsBlockPool::default());
        Self {
            slots: (0..threads.max(1))
                .map(|_| std::sync::Mutex::new(RpsScratch::new(&pool)))
                .collect(),
            pool,
        }
    }

    /// `work` on the scratch of the worker `slot` (its index in the run's pool; 0 without a
    /// pool), or on a scratch of its own where that one is not free (no task of a search waits
    /// for the pool, so a worker never holds its scratch twice) or the slot is not one of the
    /// run's.
    fn with<T>(&self, slot: usize, work: impl FnOnce(&mut RpsScratch) -> T) -> T {
        match self.slots.get(slot).and_then(|slot| slot.try_lock().ok()) {
            Some(mut scratch) => work(&mut scratch),
            None => work(&mut RpsScratch::new(&self.pool)),
        }
    }
}

/// The gapped K that the outfmt 0 footer of each query prints: for the first query of a batch,
/// `RPS_K_MULT` times the `.aux` K of the last list the batch's traceback reads, else the
/// table's K; every other query prints the table's K. Without composition-based statistics the
/// lists of the batch are read by ordinal id and only those whose rescaled PSSM exists write
/// the K (the mark is the highest); with them the lists are read query by query in the order of
/// `s_EvalueCompareHSPLists` and each writes it (the mark is that of the last query that has
/// one).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1234-1235
/// ```c
///          sbp->kbp_gap[0]->K = RPS_K_MULT * karlin_k[hsp_list->oid];
///          sbp->kbp_gap[0]->logK = log(RPS_K_MULT * karlin_k[hsp_list->oid]);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:292-298
/// ```c++
///         } else {
///             for(size_t i = 0; i < alignments.size(); i++) {
///                 CRef<CBlastAncillaryData> s(new CBlastAncillaryData(program, static_cast<int>(i),
///                                                                     sbp,
///                                                                     qinfo));
///                 retval.push_back(s);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:98-103
/// ```c++
///     const int ctx_index = query_number * context_per_query + i;
///     ...
///     if (sbp->kbp_gap) {
///         s_InitializeKarlinBlk(sbp->kbp_gap[ctx_index], &m_GappedKarlinBlk);
///     }
/// ```
/// The traceback writes only `kbp_gap[0]`, the first query's (`PR/STATE.md` point 1, checked on
/// the survey's 39 outfmt 0 reports); an invalid query has no Karlin block in its footer.
pub(crate) fn batch_gapped_k(
    batches: &[std::ops::Range<usize>],
    k_marks: &[Option<u32>],
    karlin_k: &[f64],
    table_k: f64,
    composition_based: bool,
) -> Vec<f64> {
    let mut gapped_k = vec![table_k; k_marks.len()];
    for batch in batches {
        if batch.is_empty() {
            continue;
        }
        let marks = &k_marks[batch.clone()];
        let mark = if composition_based {
            marks.iter().flatten().next_back()
        } else {
            marks.iter().flatten().max()
        };
        if let Some(mark) = mark {
            gapped_k[batch.start] = RPS_K_MULT * karlin_k[*mark as usize];
        }
    }
    gapped_k
}

/// A volume's final hit list of a query: the filters (`-max_hsps`, `-qcov_hsp_perc`), then the
/// list sorted by E-value and cut to the hit list size, before the volumes are merged.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1777
/// ```c
///     if(hit_params->options->query_cov_hsp_perc > 0 || ... max_hsps_per_subject > 0 ||
///        ...) {
///     	s_FilterBlastResults(results, hit_params->options, query_info, program_number);
///     }
///
///     /* Re-sort the hit lists according to their best e-values, because they
///        could have changed. Only do this for a database search. */
///     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
///         Blast_HSPResultsSortByEvalue(results);
///     }
///     s_BlastPruneExtraHits(results, hit_params->options->hitlist_size);
/// ```
fn finish_volume_list(
    hit_list: &mut crate::algorithm::blastp::hsp::BlastpHitList,
    max_hsps: usize,
    query_cov_hsp_perc: f64,
    query: &RpsQuery<'_>,
    hitlist_size: usize,
) {
    filter_results(hit_list, max_hsps, query_cov_hsp_perc, |frame| {
        query.context(frame).first_query_length
    });
    hit_list.sort_by_evalue();
    hit_list.prune_by_size(hitlist_size);
}

#[cfg(test)]
mod tests {
    use super::*;

    // C11 6.3.1.4: a double converts to `int` truncated toward zero when the truncated value is
    // representable; otherwise (and for NaN) the behaviour is undefined, which `int4_of` refuses.
    #[test]
    fn int4_conversions_truncate_toward_zero_inside_32_bits_only() {
        assert_eq!(int4_of(2_147_483_647.9), Some(i32::MAX));
        assert_eq!(int4_of(2_147_483_648.0), None);
        assert_eq!(int4_of(-2_147_483_648.9), Some(i32::MIN));
        assert_eq!(int4_of(-2_147_483_649.0), None);
        assert_eq!(int4_of(-0.9), Some(0));
        assert_eq!(int4_of(-2.596), Some(-2));
        assert_eq!(int4_of(64.9), Some(64));
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e300, -1e300] {
            assert_eq!(int4_of(value), None, "{value}");
        }
    }

    // blast_parameters.c:216-221: `(Int4)(scale_factor * ceil(x_dropoff * LN2 / Lambda))`, the
    // ceiling before the scale, one conversion of the product.
    #[test]
    fn the_ungapped_x_drop_is_the_ceiling_times_the_scale() {
        let lambda = 0.3176;
        // 7 bits: ceil(15.28) = 16 raw, 1,600 at scale 100.
        assert_eq!(ungapped_x_dropoff(7.0, lambda, 1.0), Some(16));
        assert_eq!(ungapped_x_dropoff(7.0, lambda, 100.0), Some(1600));
        // A tiny X-drop is still one raw unit.
        assert_eq!(ungapped_x_dropoff(1e-300, lambda, 200.0), Some(200));
        // ceil(x) * scale just inside and just outside 32 bits.
        let bits_for = |raw: f64| raw * lambda / std::f64::consts::LN_2;
        assert_eq!(
            ungapped_x_dropoff(bits_for(21_474_836.0) - 1e-6, lambda, 100.0),
            Some(2_147_483_600)
        );
        assert_eq!(
            ungapped_x_dropoff(bits_for(21_474_837.0), lambda, 100.0),
            None
        );
        assert_eq!(ungapped_x_dropoff(1e300, lambda, 1.0), None);
    }

    /// An `RpsRun` of no query with the given X-drops in bits, BLOSUM62 11/1's gapped block and
    /// no composition-based statistics.
    fn x_drop_run<'a>(
        gumbel: &'a BlastGumbelBlk,
        options: &'a RpsSearch,
        gapped_bits: f64,
        final_bits: f64,
    ) -> RpsRun<'a> {
        RpsRun {
            queries: Vec::new(),
            gapped: KarlinParams {
                lambda: 0.267,
                k: 0.041,
                h: 0.14,
                alpha: 1.9,
                beta: -30.0,
            },
            gumbel: Some(gumbel),
            cutoff_subject_length: 0,
            evalue: 10.0,
            prelim_evalue: 10.0,
            cbs: 0,
            sbp: None,
            hitlist_size: 500,
            prelim_hitlist_size: 500,
            window: 40,
            x_drop_ungapped_bits: 7.0,
            x_drop_gapped_bits: gapped_bits,
            x_drop_final_bits: final_bits,
            gap_trigger_bits: 22.0,
            options,
            program: "RPSBLAST",
        }
    }

    /// A valid query of `len` letters (between two sentinels) with BLOSUM62's ungapped block.
    fn cutoff_query(letters: &[u8], search_space: i64) -> RpsQuery<'_> {
        let kbp = KarlinParams {
            lambda: 0.3176,
            k: 0.134,
            h: 0.4012,
            alpha: 0.0,
            beta: 0.0,
        };
        RpsQuery {
            frames: vec![RpsFrame {
                frame: 0,
                aa_seq: letters,
                aa_seq_nomask: None,
                valid: true,
                kbp_std: kbp,
                search_space,
                first_query_length: letters.len() - 2,
            }],
            batch_kbp: kbp,
        }
    }

    // blast_parameters.c:922-945 and 978-990 (quoted at `RpsCutoffs`): with the Gumbel block the
    // hit cutoff is `BLAST_SpougeEtoS` of the preliminary e-value, without it
    // `BLAST_Cutoffs(.., kbp_gap, eff_searchsp)`, both times the scale (100 here).
    #[test]
    fn the_hit_cutoff_reads_the_gumbel_block_or_the_search_space_at_the_scale() {
        use crate::config::ProteinScoringSpec;
        let letters = vec![1u8; 302];
        let options = RpsSearch::default();
        let gumbel =
            crate::stats::lookup_protein_gumbel_params(&ProteinScoringSpec::default(), 3_234_567)
                .unwrap();
        let mut run = x_drop_run(&gumbel, &options, 15.0, 25.0);
        run.cutoff_subject_length = 250;
        let query = cutoff_query(&letters, 1_000_000);
        let with_block = run.cutoffs(&query.frames[0], 1, 100.0).unwrap();
        assert_eq!(
            with_block.hit_cutoff,
            crate::stats::blast_spouge_etos(10.0, &run.gapped, &gumbel, 300, 250) * 100
        );
        // Without the block: ceil(ln(K * searchsp / E) / Lambda) = ceil(ln(4100) / 0.267) = 32.
        run.gumbel = None;
        let without = run.cutoffs(&query.frames[0], 1, 100.0).unwrap();
        assert_eq!(without.hit_cutoff, 3_200);
        // The search space is the query's, not a constant.
        let other = cutoff_query(&letters, 1_000_000_000);
        assert_eq!(
            run.cutoffs(&other.frames[0], 1, 100.0).unwrap().hit_cutoff,
            5_800
        );
        // The word cutoff is the gap trigger, at most the hit cutoff (blast_parameters.c:370).
        assert!(without.word_cutoff <= without.hit_cutoff);
        // The preliminary e-value is `-evalue` as given (the stretch is the caller's).
        run.prelim_evalue = 1e300;
        assert_eq!(
            run.cutoffs(&query.frames[0], 1, 100.0).unwrap().hit_cutoff,
            100
        );
    }

    // blast_parameters.c:455-469 (quoted at `RpsRun::gapped_x_drops`): each X-drop truncated
    // toward zero, the final one the larger of its own and the preliminary one (compared as
    // doubles), both times the scale above 1; the values NCBI's probes give (DESIGN_ps11c §2).
    #[test]
    fn gapped_x_drops_are_ncbis_conversions_at_the_scale() {
        let gumbel = BlastGumbelBlk {
            lambda: 0.267,
            c: 0.0,
            g: 0.0,
            a: 0.0,
            alpha: 0.0,
            sigma: 0.0,
            a_un: 0.0,
            alpha_un: 0.0,
            b: 0.0,
            beta: 0.0,
            tau: 0.0,
            db_length: 0,
        };
        // The X-drops as given on the command line (the refusals name them).
        let given = |gapped: f64, last: f64| RpsSearch {
            x_drop_gapped_bits: Some(gapped),
            x_drop_final_bits: Some(last),
            ..RpsSearch::default()
        };
        // BLOSUM62 11/1 at the scale.
        let drops = |gapped: f64, last: f64, scale: f64| {
            let options = given(gapped, last);
            x_drop_run(&gumbel, &options, gapped, last)
                .gapped_x_drops(scale, 11 * scale as i32, scale as i32)
                .map_err(|error| error.to_string())
        };
        // The defaults: 15 bits = 38 raw, 25 bits = 64 raw.
        assert_eq!(drops(15.0, 25.0, 1.0), Ok((38, 64)));
        assert_eq!(drops(15.0, 25.0, 100.0), Ok((3800, 6400)));
        // 0 and negative gapped X-drops are kept as they are (the kernels raise them).
        assert_eq!(drops(0.0, 25.0, 1.0), Ok((0, 64)));
        assert_eq!(drops(-1.0, -1.0, 100.0), Ok((-200, -200)));
        assert_eq!(drops(-1.0, 25.0, 100.0), Ok((-200, 6400)));
        // The final X-drop is at least the preliminary one: 40 and 10 bits give 103 and 103.
        assert_eq!(drops(40.0, 10.0, 1.0), Ok((103, 103)));
        // A final X-drop far below 32 bits takes the preliminary one: no undefined conversion.
        assert_eq!(drops(15.0, -1e300, 100.0), Ok((3800, 3800)));
        for (gapped, last, scale, option) in [
            (1e300, 25.0, 1.0, "-xdrop_gap of 1e300 bits"),
            (-1e300, 25.0, 1.0, "-xdrop_gap of -1e300 bits"),
            (15.0, 1e300, 1.0, "-xdrop_gap_final of 1e300 bits"),
            // 1e7 bits = 25,960,568 raw: inside 32 bits at scale 1, outside at scale 100.
            (1e7, 25.0, 100.0, "-xdrop_gap of 10000000.0 bits"),
            (15.0, 1e7, 100.0, "-xdrop_gap_final of 10000000.0 bits"),
        ] {
            let error = drops(gapped, last, scale).unwrap_err();
            assert!(
                error.contains(option) && error.contains(&format!("scale {scale}")),
                "{error}"
            );
        }
        assert_eq!(drops(1e7, 1e7, 1.0), Ok((25_960_568, 25_960_568)));
        // The kernels' window: 4.1e8 bits = 1,064,383,310 raw, below 2^30 - 103 over a gap
        // extension of 1; 4.2e8 bits = 1,090,343,879 raw, above it: refused for the preliminary
        // X-drop, and for the final one where the CBS 0 traceback runs (this run, no redo). At
        // scale 2 the product with the scale already leaves 32 bits.
        assert!(drops(4.1e8, 25.0, 1.0).is_ok());
        assert!(drops(15.0, 4.1e8, 1.0).is_ok());
        for scale in [1.0, 2.0] {
            let error = drops(4.2e8, 25.0, scale).unwrap_err();
            assert!(error.contains("-xdrop_gap of 420000000.0 bits"), "{error}");
            let error = drops(15.0, 4.2e8, scale).unwrap_err();
            assert!(
                error.contains("-xdrop_gap_final of 420000000.0 bits"),
                "{error}"
            );
        }
        // NCBI's defaults never leave 32 bits but at an extreme scale, which the refusal names.
        let defaults = RpsSearch::default();
        let error = x_drop_run(&gumbel, &defaults, 15.0, 25.0)
            .gapped_x_drops(1e8, 1_100_000_000, 100_000_000)
            .unwrap_err()
            .to_string();
        assert!(
            error.starts_with("a profile database of the scale 100000000 (its scaled gap costs"),
            "{error}"
        );
    }

    // blast_gapalign.c:426-429, 457-464, 641-644 (quoted at `kernel_x_drop_fits`): the X-drop
    // raised to the gap costs, over the gap extension, plus 103 is the DP window; twice a window
    // of 2^30 leaves `Int4`.
    #[test]
    fn the_kernels_window_of_an_x_drop_stays_in_32_bits() {
        assert_eq!(
            crate::core::x_drop::KERNEL_X_DROP_QUOTIENT_LIMIT,
            1_073_741_721
        );
        assert!(kernel_x_drop_fits(1_073_741_720, 11, 1));
        assert!(!kernel_x_drop_fits(1_073_741_721, 11, 1));
        // Scale 2: the quotient truncates.
        assert!(kernel_x_drop_fits(2_147_483_441, 22, 2));
        assert!(!kernel_x_drop_fits(2_147_483_442, 22, 2));
        // From scale 3 on no `Int4` X-drop reaches it.
        assert!(kernel_x_drop_fits(i32::MAX, 33, 3));
        // An X-drop below the gap costs is raised to them; negative ones too.
        assert!(kernel_x_drop_fits(0, 11, 1));
        assert!(kernel_x_drop_fits(i32::MIN, 11, 1));
    }

    // blast_traceback.c:1241-1261: with composition-based statistics only the redo runs after
    // the preliminary stage, whose gap extension is 32 times the database's: the final X-drop
    // then never reaches the kernels' window limit, the preliminary one still does.
    #[test]
    fn the_final_x_drop_reaches_the_kernels_only_without_the_redo() {
        static FREQ_RATIOS: [[f64; crate::core::blast_util::BLASTAA_SIZE];
            crate::core::blast_util::BLASTAA_SIZE] =
            [[0.0; crate::core::blast_util::BLASTAA_SIZE]; crate::core::blast_util::BLASTAA_SIZE];
        static BACKGROUND: [f64; 20] = [0.05; 20];
        let sbp = crate::core::psi::PsiScoreBlk {
            matrix: [[0; crate::core::blast_util::BLASTAA_SIZE];
                crate::core::blast_util::BLASTAA_SIZE],
            ideal_lambda: 0.3176,
            ideal_k: 0.134,
            gap_std_k: 0.041,
            gap_psi: crate::core::psi::KarlinBlk {
                lambda: 0.267,
                k: 0.041,
                h: 0.14,
            },
            freq_ratios: &FREQ_RATIOS,
            bit_scale_factor: 2,
            background: &BACKGROUND,
        };
        let gumbel = BlastGumbelBlk {
            lambda: 0.267,
            c: 0.0,
            g: 0.0,
            a: 0.0,
            alpha: 0.0,
            sigma: 0.0,
            a_un: 0.0,
            alpha_un: 0.0,
            b: 0.0,
            beta: 0.0,
            tau: 0.0,
            db_length: 0,
        };
        let options = RpsSearch {
            x_drop_gapped_bits: Some(15.0),
            x_drop_final_bits: Some(4.2e8),
            ..RpsSearch::default()
        };
        let mut run = x_drop_run(&gumbel, &options, 15.0, 4.2e8);
        run.cbs = 1;
        run.sbp = Some(&sbp);
        // A final X-drop of 4.2e8 bits leaves 32 bits at the redo's Lambda: the redo's check.
        let error = run.gapped_x_drops(1.0, 11, 1).unwrap_err().to_string();
        assert!(error.contains("an -xdrop_gap_final of"), "{error}");
        // (At the redo's Lambda a final X-drop is about 32 times its raw value, so one above the
        // kernels' window limit is always refused by this conversion.) 1.5e7 bits are
        // 1,246,107,290 at the redo's Lambda: accepted.
        run.x_drop_final_bits = 1.5e7;
        assert!(run.gapped_x_drops(1.0, 11, 1).is_ok());
        // The preliminary X-drop is checked at every mode.
        run.x_drop_gapped_bits = 4.2e8;
        let error = run.gapped_x_drops(1.0, 11, 1).unwrap_err().to_string();
        assert!(error.contains("an -xdrop_gap of"), "{error}");
    }

    // blast_kappa.c:2381-2383 (quoted at `RedoRun::new`): the redo's X-drop at Lambda over the
    // scaling factor 32; 25 bits give 2,076 raw, below the scaled final X-drop 6,400 of scale 100.
    #[test]
    fn the_redo_x_drop_is_the_larger_of_the_bits_at_the_redo_lambda_and_the_final_one() {
        let lambda = 0.267 / SCALING_FACTOR;
        assert_eq!(
            int4_of(RedoRun::x_dropoff_of(25.0, lambda, 6400)),
            Some(6400)
        );
        assert_eq!(int4_of(RedoRun::x_dropoff_of(25.0, lambda, 64)), Some(2076));
        assert_eq!(
            int4_of(RedoRun::x_dropoff_of(-1.0, lambda, -200)),
            Some(-83)
        );
        // 2.5e7 bits = 2,076,845,484 raw: inside; 1e8 bits: outside.
        assert!(int4_of(RedoRun::x_dropoff_of(2.5e7, lambda, 64)).is_some());
        assert_eq!(int4_of(RedoRun::x_dropoff_of(1e8, lambda, 64)), None);
        assert_eq!(int4_of(RedoRun::x_dropoff_of(-1e300, lambda, 64)), Some(64));
    }

    #[test]
    fn the_first_query_of_a_batch_gets_the_k_of_the_highest_mark() {
        let marks = [Some(2), None, Some(5), None, None, Some(1)];
        let karlin_k = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06];
        let k = batch_gapped_k(&[0..3, 3..5, 5..6], &marks, &karlin_k, 0.041, false);
        assert_eq!(k[0], RPS_K_MULT * 0.06);
        assert_eq!(k[1], 0.041);
        assert_eq!(k[2], 0.041);
        assert_eq!(k[3], 0.041);
        assert_eq!(k[5], RPS_K_MULT * 0.02);
        // With composition-based statistics: the mark of the last query that has one, whatever
        // its ordinal id.
        let marks = [Some(5), None, Some(2), None, None, Some(1)];
        let k = batch_gapped_k(&[0..3, 3..5, 5..6], &marks, &karlin_k, 0.041, true);
        assert_eq!(k[0], RPS_K_MULT * 0.03);
        assert_eq!(k[3], 0.041);
        assert_eq!(k[5], RPS_K_MULT * 0.02);
        let k = batch_gapped_k(&[0..3, 3..5, 5..6], &marks, &karlin_k, 0.041, false);
        assert_eq!(k[0], RPS_K_MULT * 0.06);
    }

    fn hsp(
        oid: u32,
        e_value: f64,
        q_start: usize,
        q_end: usize,
    ) -> crate::algorithm::blastp::hsp::BlastpHsp {
        let mut hsp = crate::algorithm::blastp::hsp::BlastpHsp::from_hit(crate::common::Hit {
            identity: 0.0,
            length: q_end - q_start + 1,
            mismatch: 0,
            gapopen: 0,
            q_start,
            q_end,
            s_start: 1,
            s_end: 10,
            e_value,
            bit_score: 10.0,
            num_ident: 0,
            query_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx: oid,
            raw_score: 20,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 0,
        });
        hsp.e_value = e_value;
        hsp
    }

    // tabular.cpp:973-1007 and showalign.cpp:2120-2135 (quoted at `formatter_identities`): the
    // reports count the identities on the query as read and the profile's consensus, across the
    // gaps of the edit script (after `swap_roles`: Del steps the profile, Ins the query).
    #[test]
    fn the_reports_identities_count_the_letters_as_read() {
        use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHspList};
        use crate::common::GapEditOp;
        // Profile 1 starts at 6 of the concatenated consensus (profile 0 holds 5 letters + 1).
        let enc = |text: &str| -> Vec<u8> {
            text.bytes()
                .map(crate::utils::matrix::aa_char_to_ncbistdaa)
                .collect()
        };
        let mut consensus = enc("AAAAA");
        consensus.push(0);
        consensus.extend(enc("MKTAYWWIAKQ"));
        let offsets = [0, 6, 18];
        // The query as read (its searched letters could have had X over KTAY, a lower-case run;
        // the traceback's count of 5 below is that of the masked letters).
        let letters = enc("GMKTAYIAKQ");
        // The letters as read between two sentinels (`RpsFrame::aa_seq_nomask`).
        let mut nomask = vec![0u8];
        nomask.extend(&letters);
        nomask.push(0);
        let query = RpsQuery {
            frames: vec![RpsFrame {
                frame: 0,
                aa_seq: &nomask,
                aa_seq_nomask: Some(&nomask),
                valid: true,
                kbp_std: KarlinParams::default(),
                search_space: 0,
                first_query_length: letters.len(),
            }],
            batch_kbp: KarlinParams::default(),
        };
        // Query 2..10 against profile 1..11: MKTAY / MKTAY, a gap of WW in the query, IAKQ / IAKQ.
        let mut aligned = hsp(1, 1e-5, 2, 10);
        aligned.s_start = 1;
        aligned.s_end = 11;
        aligned.length = 11;
        aligned.num_ident = 5;
        aligned.mismatch = 4;
        aligned.gap_info = Some(vec![
            GapEditOp::Sub(5),
            GapEditOp::Del(2),
            GapEditOp::Sub(4),
        ]);
        let mut hit_list = BlastpHitList::new(10);
        hit_list.update(BlastpHspList {
            oid: 1,
            query_index: 0,
            hsps: vec![aligned],
            best_evalue: 1e-5,
        });
        formatter_identities(&mut hit_list, &query, &consensus, &offsets);
        let counted = &hit_list.hsplist_array[0].hsps[0];
        assert_eq!(counted.num_ident, 9);
        assert_eq!(counted.mismatch, 0);
        assert_eq!(counted.identity, 900.0 / 11.0);
        // An Ins steps the query alone.
        let mut inserted = hsp(1, 1e-5, 1, 6);
        inserted.s_start = 1;
        inserted.s_end = 5;
        inserted.length = 6;
        inserted.gap_info = Some(vec![GapEditOp::Ins(1), GapEditOp::Sub(5)]);
        let mut hit_list = BlastpHitList::new(10);
        hit_list.update(BlastpHspList {
            oid: 1,
            query_index: 0,
            hsps: vec![inserted],
            best_evalue: 1e-5,
        });
        formatter_identities(&mut hit_list, &query, &consensus, &offsets);
        assert_eq!(hit_list.hsplist_array[0].hsps[0].num_ident, 5);
    }

    #[test]
    fn a_coverage_filter_that_removes_a_best_hsp_reorders_and_cuts_the_volumes_list() {
        // The batch's first query has 100 letters.
        let letters = [0u8; 102];
        let query = RpsQuery {
            frames: vec![RpsFrame {
                frame: 0,
                aa_seq: &letters,
                aa_seq_nomask: None,
                valid: true,
                kbp_std: KarlinParams::default(),
                search_space: 0,
                first_query_length: 100,
            }],
            batch_kbp: KarlinParams::default(),
        };
        use crate::algorithm::blastp::hsp::{BlastpHitList, BlastpHspList};
        let list = |oid: u32, hsps: Vec<_>, best: f64| BlastpHspList {
            oid,
            query_index: 0,
            hsps,
            best_evalue: best,
        };
        let build = || {
            let mut hit_list = BlastpHitList::new(10);
            // Subject 1: its best HSP covers 10 % of the 100-letter query, its second 80 %.
            hit_list.update(list(
                1,
                vec![hsp(1, 1e-10, 1, 10), hsp(1, 1e-3, 11, 90)],
                1e-10,
            ));
            hit_list.update(list(2, vec![hsp(2, 1e-5, 1, 90)], 1e-5));
            hit_list
        };
        let oids = |hit_list: &BlastpHitList| -> Vec<u32> {
            hit_list.hsplist_array[..hit_list.hsplist_count]
                .iter()
                .map(|list| list.oid)
                .collect()
        };
        let mut hit_list = build();
        hit_list.sort_by_evalue();
        assert_eq!(oids(&hit_list), vec![1, 2]);
        // 50 % coverage drops subject 1's best HSP: its best E-value becomes 1e-3, after 2's.
        let mut hit_list = build();
        finish_volume_list(&mut hit_list, 0, 50.0, &query, 10);
        assert_eq!(oids(&hit_list), vec![2, 1]);
        // And the list is cut to the hit list size after the sort.
        let mut hit_list = build();
        finish_volume_list(&mut hit_list, 0, 50.0, &query, 1);
        assert_eq!(oids(&hit_list), vec![2]);
        // Without a filter the order is the list's own.
        let mut hit_list = build();
        finish_volume_list(&mut hit_list, 0, 0.0, &query, 10);
        assert_eq!(oids(&hit_list), vec![1, 2]);
    }

    /// The files of the fixture database `rps_mini` (tests/db_search/rps/, three profiles)
    /// inflated into a directory of this test, and the database opened.
    fn rps_mini(name: &str) -> (std::path::PathBuf, crate::seqdb::SeqDb) {
        rps_mini_without(name, &[])
    }

    /// `rps_mini` without the files of the extensions `left_out`.
    fn rps_mini_without(
        name: &str,
        left_out: &[&str],
    ) -> (std::path::PathBuf, crate::seqdb::SeqDb) {
        use std::io::Read;
        let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/rps");
        let dir =
            std::env::temp_dir().join(format!("losat-rps-engine-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for extension in [
            "aux", "loo", "rps", "pin", "psq", "phr", "pdb", "pos", "pot", "ptf", "pto",
        ] {
            if left_out.contains(&extension) {
                continue;
            }
            let mut bytes = Vec::new();
            flate2::read::GzDecoder::new(
                std::fs::read(fixtures.join(format!("rps_mini.{extension}.gz")))
                    .unwrap()
                    .as_slice(),
            )
            .read_to_end(&mut bytes)
            .unwrap();
            std::fs::write(dir.join(format!("rps_mini.{extension}")), bytes).unwrap();
        }
        let db = crate::seqdb::SeqDb::open(
            dir.join("rps_mini").to_str().unwrap(),
            crate::seqdb::SeqType::Protein,
        )
        .unwrap();
        (dir, db)
    }

    /// NCBIstdaa letters of `text` between two sentinels.
    fn encoded(text: &str) -> Vec<u8> {
        const LETTERS: &[u8] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
        let mut letters = vec![0u8];
        letters.extend(text.bytes().map(|letter| {
            LETTERS
                .iter()
                .position(|&l| l == letter.to_ascii_uppercase())
                .unwrap_or(21) as u8
        }));
        letters.push(0);
        letters
    }

    // D2 of PS12s: a volume's search with every valid context of every query a unit of the pool
    // (the query finished by the worker that searches its last context) gives each query the
    // results of the search one context after another on one worker; a query whose cutoffs are
    // refused gets that refusal as its result, the others theirs.
    #[test]
    fn units_of_contexts_give_the_results_of_the_search_in_order() {
        use crate::algorithm::rps_search::volumes::{Consensus, RpsVolumeSet, RpsVolumeSetup};
        let (dir, db) = rps_mini("units");
        let volumes = RpsVolumeSet::of(&db);
        let volume = &volumes.volumes[0];
        let mut setup = RpsVolumeSetup::new("RPSTBLASTN", 40);
        assert!(setup.check(volume, false).unwrap().is_ok());
        let Ok(files) = setup.reopen(volume, false).unwrap() else {
            panic!("rps_mini opens");
        };
        // The profiles' fragments (tests/db_search/rps/build.sh) as the frames of queries: windows
        // of the WSSV proteins with the profiles' letters, and two queries of other letters.
        let proteins = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/rps/q40.faa"),
        )
        .unwrap();
        let sequences: Vec<String> = proteins
            .split('>')
            .filter(|record| !record.is_empty())
            .map(|record| record.lines().skip(1).collect())
            .collect();
        let windows = [
            (32, 0, 80),
            (32, 60, 190),
            (18, 300, 400),
            (18, 380, 480),
            (4, 1100, 1200),
            (4, 1180, 1290),
            (0, 0, 120),
            (1, 10, 90),
        ];
        let frames: Vec<Vec<u8>> = (0..18)
            .map(|i| {
                let (record, start, end) = windows[i % windows.len()];
                let sequence = &sequences[record];
                encoded(&sequence[start.min(sequence.len())..end.min(sequence.len())])
            })
            .collect();
        let kbp = KarlinParams {
            lambda: 0.3176,
            k: 0.134,
            h: 0.4012,
            alpha: 0.0,
            beta: 0.0,
        };
        let options = RpsSearch::default();
        // BLOSUM62 11/1's Gumbel block (the database's 1,000,000 letters).
        let gumbel = crate::stats::spouge::lookup_protein_gumbel_params(
            &crate::config::ProteinScoringSpec {
                matrix: crate::config::ScoringMatrix::Blosum62,
                gap_open: 11,
                gap_extend: 1,
            },
            1_000_000,
        )
        .unwrap();
        // Six queries of six frames (one of them invalid), the second with a Lambda that makes its
        // ungapped X-drop leave 32 bits (refused: NCBI's default X-drop, so the refusal names the
        // database's scale).
        let queries: Vec<RpsQuery<'_>> = (0..6)
            .map(|q| RpsQuery {
                frames: (0..NUM_FRAMES)
                    .map(|f| {
                        let aa_seq = &frames[(q * NUM_FRAMES + f) % frames.len()];
                        RpsFrame {
                            frame: [1, 2, 3, -1, -2, -3][f],
                            aa_seq,
                            aa_seq_nomask: None,
                            valid: f != 4,
                            kbp_std: KarlinParams {
                                lambda: if q == 1 { 1e-300 } else { kbp.lambda },
                                ..kbp
                            },
                            search_space: 20_000_000,
                            first_query_length: aa_seq.len() - 2,
                        }
                    })
                    .collect(),
                batch_kbp: kbp,
            })
            .collect();
        let run = RpsRun {
            queries,
            gapped: KarlinParams {
                lambda: 0.267,
                k: 0.041,
                h: 0.14,
                alpha: 1.9,
                beta: -30.0,
            },
            gumbel: Some(&gumbel),
            cutoff_subject_length: 180,
            evalue: 10.0,
            prelim_evalue: 10.0,
            cbs: 0,
            sbp: None,
            hitlist_size: 500,
            prelim_hitlist_size: 500,
            window: 40,
            x_drop_ungapped_bits: 7.0,
            x_drop_gapped_bits: 15.0,
            x_drop_final_bits: 25.0,
            gap_trigger_bits: 22.0,
            options: &options,
            program: "RPSTBLASTN",
        };
        let consensus = Consensus::new(&setup, volume, &files);
        let shown = |results: RpsVolumeResults| -> Vec<String> {
            assert!(results.late.is_none());
            results
                .refused
                .into_iter()
                .zip(results.searched)
                .map(|(refusal, result)| match (refusal, result) {
                    (Some(refusal), None) => format!("refused: {refusal}"),
                    (None, Some(Ok(traced))) => {
                        format!("{:?} {:?}", traced.hit_list, traced.k_mark)
                    }
                    (None, Some(Err(error))) => format!("failed: {error}"),
                    _ => unreachable!("a query is refused or searched"),
                })
                .collect()
        };
        let serial = shown(
            run.search_volume(
                &files,
                &consensus,
                &crate::utils::threading::SearchPool::x_serial(),
                &RpsWorkers::new(1),
            )
            .unwrap(),
        );
        assert_eq!(serial.len(), 6);
        assert!(
            serial[1].starts_with("refused: a profile database of the scale 100"),
            "{}",
            serial[1]
        );
        let with_hits = serial
            .iter()
            .filter(|shown| shown.starts_with("Some(BlastpHitList"))
            .count();
        assert!(with_hits >= 4, "{serial:?}");
        for threads in [2, 4] {
            let workers = RpsWorkers::new(threads);
            let parallel =
                crate::utils::threading::with_search_pool(threads, "RPSTBLASTN", |pool| {
                    // Twice on the same workers: the kept scratch gives the same results.
                    let first = shown(run.search_volume(&files, &consensus, pool, &workers)?);
                    let again = shown(run.search_volume(&files, &consensus, pool, &workers)?);
                    assert_eq!(first, again);
                    Ok(first)
                })
                .unwrap();
            assert_eq!(parallel, serial, "{threads} threads");
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    // S1 and N1 of the review of PS12s: a batch's refusal of a query's cutoffs comes before the
    // search errors of the batch's earlier queries, and the volume's refusal of its Lambda at the
    // scale after the cutoffs (NCBI sets the cutoffs of the whole batch before its scan,
    // blast_engine.c:1372-1374, quoted at `search_volume`; its e-values come after). Query 0 has
    // hits, so its traceback reads the consensus letters, whose file (.psq) is missing: its search
    // fails; query 1's ungapped X-drop leaves 32 bits: refused at its cutoffs.
    #[test]
    fn a_batchs_cutoff_refusal_comes_before_its_search_errors() {
        use crate::algorithm::rps_search::volumes::{Consensus, RpsVolumeSet, RpsVolumeSetup};
        let proteins = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/rps/q40.faa"),
        )
        .unwrap();
        let record: String = proteins
            .split('>')
            .filter(|record| !record.is_empty())
            .nth(32)
            .unwrap()
            .lines()
            .skip(1)
            .collect();
        let letters = encoded(&record[..80]);
        let options = RpsSearch::default();
        let gumbel = crate::stats::spouge::lookup_protein_gumbel_params(
            &crate::config::ProteinScoringSpec::default(),
            1_000_000,
        )
        .unwrap();
        let queries = || {
            let mut queries = vec![
                cutoff_query(&letters, 20_000_000),
                cutoff_query(&letters, 20_000_000),
            ];
            queries[1].frames[0].kbp_std.lambda = 1e-300;
            queries
        };
        let run = |lambda: f64| RpsRun {
            queries: queries(),
            gapped: KarlinParams {
                lambda,
                ..x_drop_run(&gumbel, &options, 15.0, 25.0).gapped
            },
            cutoff_subject_length: 180,
            ..x_drop_run(&gumbel, &options, 15.0, 25.0)
        };
        let search = |left_out: &[&str], lambda: f64| {
            let (dir, db) = rps_mini_without("cutoffs", left_out);
            let volumes = RpsVolumeSet::of(&db);
            let volume = &volumes.volumes[0];
            let mut setup = RpsVolumeSetup::new("RPSBLAST", 40);
            assert!(setup.check(volume, false).unwrap().is_ok());
            let Ok(files) = setup.reopen(volume, false).unwrap() else {
                panic!("rps_mini opens");
            };
            let consensus = Consensus::new(&setup, volume, &files);
            let results = run(lambda)
                .search_volume(
                    &files,
                    &consensus,
                    &crate::utils::threading::SearchPool::x_serial(),
                    &RpsWorkers::new(1),
                )
                .unwrap();
            let _ = std::fs::remove_dir_all(dir);
            results
        };
        let errors = |results: RpsVolumeResults, spans: &[usize]| -> Vec<String> {
            results
                .batches(spans)
                .into_iter()
                .map(|batch| match batch {
                    Ok(traced) => format!("{} lists", traced.len()),
                    Err(error) => error.to_string(),
                })
                .collect()
        };
        let lambda = 0.267;
        // With every file: query 0 has a hit list, query 1 is refused.
        let whole = search(&[], lambda);
        assert!(whole.refused[1].is_some() && whole.refused[0].is_none());
        assert!(matches!(&whole.searched[0], Some(Ok(traced)) if traced.hit_list.is_some()));
        let scale_refusal = "a profile database of the scale 100 (its scaled";
        assert!(errors(whole, &[2])[0].starts_with(scale_refusal));
        // Without the consensus letters: query 0's search fails; in one batch query 1's refusal
        // still comes first, in a batch of its own query 0's error ends the batches.
        let missing = |spans: &[usize]| errors(search(&["psq"], lambda), spans);
        assert!(
            missing(&[2])[0].starts_with(scale_refusal),
            "{:?}",
            missing(&[2])
        );
        let alone = missing(&[1, 1]);
        assert_eq!(alone.len(), 1);
        assert!(alone[0].contains("rps_mini.psq"), "{alone:?}");
        // A Lambda that the scale does not give back (0.4389734947748931 / 100 * 100 differs):
        // the volume's refusal comes after the batch's cutoff refusals, before any search.
        let late = 0.438_973_494_774_893_1;
        assert_ne!(late / 100.0 * 100.0, late);
        let lambda_refusal = "a profile database of the scale 100 (NCBI BLAST+ divides Lambda";
        let refused = search(&[], late);
        assert!(refused.searched.iter().all(Option::is_none));
        assert!(errors(refused, &[2])[0].starts_with(scale_refusal));
        let first = errors(search(&[], late), &[1, 1]);
        assert_eq!(first.len(), 1);
        assert!(first[0].starts_with(lambda_refusal), "{first:?}");
        // `all` reads the run's queries as one batch (rpsblast's).
        let Err(error) = search(&[], lambda).all() else {
            panic!("the run's queries are refused");
        };
        assert!(error.to_string().starts_with(scale_refusal), "{error}");
    }
}
