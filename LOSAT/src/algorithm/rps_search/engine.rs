//! The search of a run's queries against one profile volume: every (query, volume) pair alone,
//! the queries in parallel in the run's pool, each query's preliminary stage and traceback on a
//! worker's scratch, the results in query order. The outfmt 0 K of each query batch is computed
//! from the queries' marks after (`batch_gapped_k`).
//!
//! NCBI searches the queries of a batch one after another on one diagonal table, then traces the
//! lists of the whole batch; a query's results depend on the volume, the options and the query
//! only (the table's entries of earlier queries read as empty, `wordfinder.rs`), so LOSAT's
//! bytes are the same at every thread count.
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
use super::scan::{RpsBuckets, RpsTable};
use super::traceback::{
    filter_results, formatter_identities, RpsTraceScratch, RpsTraceback, RpsTraced,
};
use super::wordfinder::{RpsDiagTable, RpsWordParams, RpsWordScratch};
use super::{RpsRows, RPS_K_MULT};
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

/// One user query as the search reads it.
pub(crate) struct RpsQuery<'a> {
    /// The query's letters between two sentinels (`QueryContext::aa_seq`).
    pub(crate) aa_seq: &'a [u8],
    /// Whether its Karlin-Altschul parameters exist (an invalid query gets no HSP: its cutoffs
    /// are `INT4_MAX`, blast_parameters.c:916-919).
    pub(crate) valid: bool,
    /// `kbp_std`: the query's composition with BLOSUM62.
    pub(crate) kbp_std: KarlinParams,
    /// The length of the first query of the query's batch (`-qcov_hsp_perc`).
    pub(crate) first_query_length: usize,
    /// The ungapped block of the first valid query of the query's batch: `kbp_psi[0]` of the
    /// score block that the composition-based redo copies, which stays where the profile has no
    /// score frequencies.
    pub(crate) batch_kbp: KarlinParams,
    /// The query's letters as read (between two sentinels) where the searched ones differ (SEG,
    /// lower-case masks, `Blast_MaskUnsupportedAA`): the reports' identities count on them
    /// (`traceback::formatter_identities`).
    pub(crate) aa_seq_nomask: Option<&'a [u8]>,
    /// `query_info->contexts[0].eff_searchsp` of the query in this volume (NCBI's `Int8`): what
    /// the e-values and the hit cutoff read without the Gumbel block (`EvalueStats`), and the
    /// Smith-Waterman redo with it (series PR, item PS11c-3). 0 for an invalid query.
    pub(crate) search_space: i64,
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
}

/// The refusal of a database whose scale makes a scaled score leave NCBI's `Int4` (gap costs,
/// X-drops, cutoffs): NCBI's products overflow there.
fn scale_overflow(scale: f64) -> anyhow::Error {
    anyhow::anyhow!(
        "a profile database of the scale {scale} (its scaled gap costs, X-drops or score cutoffs leave 32-bit integers) is not supported by LOSAT's RPSBLAST"
    )
}

/// The refusal of an X-drop in bits whose raw score (the bits over Lambda, converted to `Int4`,
/// times the database's scale where NCBI multiplies) leaves NCBI's `Int4`, or whose quotient by
/// the scaled gap extension makes NCBI's kernels' sizes leave it (`kernel_x_drop_fits`): NCBI's
/// conversion of the double, or its arithmetic, is undefined there (DECISIONS 17 N1/N2, line PS).
/// `given` is the option's value as the command line gives it; an X-drop of NCBI's defaults gets
/// the database's refusal (`scale_overflow`), as only an extreme scale makes it leave 32 bits.
fn x_drop_overflow(option: &str, given: Option<f64>, bits: f64, scale: f64) -> anyhow::Error {
    if given.is_none() {
        return scale_overflow(scale);
    }
    anyhow::anyhow!(
        "an {option} of {bits:?} bits, whose raw X-drop at the profile database's scale {scale} leaves NCBI BLAST+'s 32-bit integers (its conversion or the sizes of its alignment kernels overflow: undefined behaviour), is not supported by LOSAT's RPSBLAST"
    )
}

/// The largest quotient of a raw X-drop by the scaled gap extension that NCBI's gapped kernels
/// keep in `Int4`, plus one: `num_extra_cells = x_dropoff / gap_extend + 3`, the DP window
/// `dp_mem_alloc = num_extra_cells + 100`, and its growth `2 * dp_mem_alloc`, which leaves `Int4`
/// once `dp_mem_alloc` is 2^30 (at the first row that passes the window's first 97 cells).
const KERNEL_X_DROP_QUOTIENT_LIMIT: i64 = (1 << 30) - 103;

/// Whether NCBI's gapped kernels (`Blast_SemiGappedAlign`, the preliminary stage's, and
/// `ALIGN_EX`, the CBS 0 traceback's) keep their sizes in `Int4` for the raw X-drop `x_dropoff`
/// at the scale with the scaled gap costs: the X-drop raised to `gap_open + gap_extend`, over
/// `gap_extend` (`Int4` division), below `KERNEL_X_DROP_QUOTIENT_LIMIT`. Reachable only where the
/// scaled gap extension is 1 or 2 (scale 1 or 2 with BLOSUM62 11/1): the quotient is at most
/// `INT4_MAX / scale` otherwise. LOSAT's kernels size their memory by the sequences, but NCBI's
/// overflow is undefined behaviour, so such an X-drop is refused.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:426-429
/// ```c
///     x_dropoff = gap_align->gap_x_dropoff;
///
///     if (x_dropoff < gap_open_extend)
///         x_dropoff = gap_open_extend;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:457-464
/// ```c
///     if (gap_extend > 0)
///         num_extra_cells = x_dropoff / gap_extend + 3;
///     else
///         num_extra_cells = N + 3;
///
///     if (num_extra_cells > gap_align->dp_mem_alloc) {
///         gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
///                                       2 * gap_align->dp_mem_alloc);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:641-644 (the same at
/// 923-926 in `Blast_SemiGappedAlign`)
/// ```c
///         if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {
///
///             gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
///                                           2 * gap_align->dp_mem_alloc);
/// ```
fn kernel_x_drop_fits(x_dropoff: i32, gap_open: i32, gap_extend: i32) -> bool {
    let gap_open_extend = i64::from(gap_open) + i64::from(gap_extend);
    let x_dropoff = i64::from(x_dropoff).max(gap_open_extend);
    gap_extend <= 0 || x_dropoff / i64::from(gap_extend) < KERNEL_X_DROP_QUOTIENT_LIMIT
}

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
    /// What the query's e-values and hit cutoff read in this volume: the Gumbel block, if any,
    /// and the query's effective search space.
    fn stats(&self, query: &RpsQuery<'_>) -> EvalueStats<'_> {
        EvalueStats::new(self.gumbel, query.search_space)
    }

    /// The query's cutoffs; the refusal where the ungapped X-drop (`-xdrop_ungap`) or a scaled
    /// cutoff leaves `i32`.
    fn cutoffs(&self, query: &RpsQuery<'_>, scale: f64) -> Result<RpsCutoffs> {
        let int_scale = scale as i32;
        let x_dropoff = ungapped_x_dropoff(self.x_drop_ungapped_bits, query.kbp_std.lambda, scale)
            .ok_or_else(|| {
                x_drop_overflow(
                    "-xdrop_ungap",
                    self.options.x_drop_ungapped_bits,
                    self.x_drop_ungapped_bits,
                    scale,
                )
            })?;
        let hit_cutoff = self
            .stats(query)
            .hit_cutoff(
                self.prelim_evalue,
                &self.gapped,
                (query.aa_seq.len() - 2) as i32,
                self.cutoff_subject_length,
            )
            .checked_mul(int_scale)
            .ok_or_else(|| scale_overflow(scale))?;
        let gap_trigger = ((self.gap_trigger_bits * std::f64::consts::LN_2 + query.kbp_std.k.ln())
            / query.kbp_std.lambda) as i32;
        Ok(RpsCutoffs {
            x_dropoff,
            word_cutoff: gap_trigger
                .checked_mul(int_scale)
                .ok_or_else(|| scale_overflow(scale))?
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
            )
        };
        let final_overflow = || {
            x_drop_overflow(
                "-xdrop_gap_final",
                self.options.x_drop_final_bits,
                self.x_drop_final_bits,
                scale,
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
        // The redo converts its X-drop at its own Lambda (`RedoRun::with_final_x_drop_bits`);
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
    /// `consensus` (at the profiles' offsets), in query order.
    pub(crate) fn search_volume(
        &self,
        files: &RpsVolumeFiles,
        consensus: &[u8],
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Result<Vec<RpsTraced>> {
        let rows = files.rows();
        let offsets = files.offsets();
        let scale = files.aux.scale;
        // (A whole number of 1 or more, `rps_hooks::rps_search`.)
        if scale > f64::from(i32::MAX) {
            return Err(scale_overflow(scale));
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
            return Err(scale_overflow(scale));
        };
        // The X-drops at the scale (`gapped_x_drops`, NCBI's `BlastExtensionParametersNew`).
        let (x_drop_gapped, x_drop_final) = self.gapped_x_drops(scale, gap_open, gap_extend)?;
        // Every searched query's cutoffs, before the search (`BlastHitSavingParametersNew` and
        // `BlastInitialWordParametersNew` set them for the whole batch before its scan).
        let cutoffs = self
            .queries
            .iter()
            .map(|query| {
                if query.valid {
                    self.cutoffs(query, scale)
                } else {
                    Ok(RpsCutoffs {
                        x_dropoff: 0,
                        word_cutoff: 0,
                        hit_cutoff: 0,
                    })
                }
            })
            .collect::<Result<Vec<_>>>()?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1862-1895
        // ```c
        //       kbp[kbp_context]->Lambda /= scaling_factor;
        //       ...
        //       /* Put back the unscaled value of Lambda. */
        //       kbp[kbp_context]->Lambda *= scaling_factor;
        // ```
        // NCBI divides Lambda by the scale for each HSP's e-value and multiplies it back; at a
        // scale where that does not give Lambda back, NCBI's Lambda changes from HSP to HSP.
        let lambda = self.gapped.lambda;
        if lambda / scale * scale != lambda {
            anyhow::bail!(
                "a profile database of the scale {scale} (NCBI BLAST+ divides Lambda by the scale for each e-value and multiplies it back, which changes Lambda at this scale) is not supported by LOSAT's RPSBLAST"
            );
        }
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
        let new_scratch = || {
            let mut diag = RpsDiagTable::new();
            diag.reset_for(rows.len() as i32, self.window);
            (
                RpsPrelimScratch::new(RpsWordScratch {
                    diag,
                    buckets: RpsBuckets::new(rows.len()),
                }),
                RpsTraceScratch::new(),
            )
        };
        let search_one = |scratch: &mut (RpsPrelimScratch, RpsTraceScratch),
                          q_idx: usize|
         -> Result<RpsTraced> {
            let query = &self.queries[q_idx];
            let none = RpsTraced {
                hit_list: None,
                k_mark: None,
            };
            if !query.valid {
                return Ok(none);
            }
            let cutoffs = &cutoffs[q_idx];
            let prelim = RpsPrelim {
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
                stats: self.stats(query),
                prelim_evalue: self.prelim_evalue,
                prelim_hitlist_size: self.prelim_hitlist_size,
                query_index: q_idx as u32,
            };
            let subject_seq = &query.aa_seq[1..];
            let Some(mut lists) = prelim.search_query(&mut scratch.0, subject_seq)? else {
                return Ok(none);
            };
            if cbs_close {
                lists.cbs_close_trim(self.hitlist_size);
            }
            let letters = &query.aa_seq[1..query.aa_seq.len() - 1];
            let mut traced = traceback.trace_query(
                lists,
                letters,
                q_idx as u32,
                &query.batch_kbp,
                self.stats(query),
                &mut scratch.1,
            )?;
            if let (Some(hit_list), Some(nomask)) = (traced.hit_list.as_mut(), query.aa_seq_nomask)
            {
                formatter_identities(hit_list, &nomask[1..nomask.len() - 1], consensus, offsets);
            }
            if let Some(hit_list) = traced.hit_list.as_mut() {
                finish_volume_list(
                    hit_list,
                    self.options.max_hsps_per_subject,
                    self.options.query_cov_hsp_perc,
                    query.first_query_length,
                    self.hitlist_size,
                );
            }
            Ok(traced)
        };
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() && pool.threads() > 1 {
            use rayon::prelude::*;
            return pool.install(|| {
                (0..self.queries.len())
                    .into_par_iter()
                    .map_init(new_scratch, search_one)
                    .collect()
            });
        }
        let _ = pool;
        let mut scratch = new_scratch();
        (0..self.queries.len())
            .map(|q_idx| search_one(&mut scratch, q_idx))
            .collect()
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
    first_query_length: usize,
    hitlist_size: usize,
) {
    filter_results(hit_list, max_hsps, query_cov_hsp_perc, first_query_length);
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
        }
    }

    /// A valid query of `len` letters (between two sentinels) with BLOSUM62's ungapped block.
    fn cutoff_query(letters: &[u8], search_space: i64) -> RpsQuery<'_> {
        RpsQuery {
            aa_seq: letters,
            valid: true,
            kbp_std: KarlinParams {
                lambda: 0.3176,
                k: 0.134,
                h: 0.4012,
                alpha: 0.0,
                beta: 0.0,
            },
            first_query_length: letters.len() - 2,
            batch_kbp: KarlinParams {
                lambda: 0.3176,
                k: 0.134,
                h: 0.4012,
                alpha: 0.0,
                beta: 0.0,
            },
            aa_seq_nomask: None,
            search_space,
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
        let with_block = run.cutoffs(&query, 100.0).unwrap();
        assert_eq!(
            with_block.hit_cutoff,
            crate::stats::blast_spouge_etos(10.0, &run.gapped, &gumbel, 300, 250) * 100
        );
        // Without the block: ceil(ln(K * searchsp / E) / Lambda) = ceil(ln(4100) / 0.267) = 32.
        run.gumbel = None;
        let without = run.cutoffs(&query, 100.0).unwrap();
        assert_eq!(without.hit_cutoff, 3_200);
        // The search space is the query's, not a constant.
        let other = cutoff_query(&letters, 1_000_000_000);
        assert_eq!(run.cutoffs(&other, 100.0).unwrap().hit_cutoff, 5_800);
        // The word cutoff is the gap trigger, at most the hit cutoff (blast_parameters.c:370).
        assert!(without.word_cutoff <= without.hit_cutoff);
        // The preliminary e-value is `-evalue` as given (the stretch is the caller's).
        run.prelim_evalue = 1e300;
        assert_eq!(run.cutoffs(&query, 100.0).unwrap().hit_cutoff, 100);
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
        assert_eq!(KERNEL_X_DROP_QUOTIENT_LIMIT, 1_073_741_721);
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
        let query = enc("GMKTAYIAKQ");
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
        finish_volume_list(&mut hit_list, 0, 50.0, 100, 10);
        assert_eq!(oids(&hit_list), vec![2, 1]);
        // And the list is cut to the hit list size after the sort.
        let mut hit_list = build();
        finish_volume_list(&mut hit_list, 0, 50.0, 100, 1);
        assert_eq!(oids(&hit_list), vec![2]);
        // Without a filter the order is the list's own.
        let mut hit_list = build();
        finish_volume_list(&mut hit_list, 0, 0.0, 100, 10);
        assert_eq!(oids(&hit_list), vec![1, 2]);
    }
}
