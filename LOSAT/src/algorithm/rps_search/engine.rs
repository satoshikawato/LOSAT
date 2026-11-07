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
use super::scan::{RpsBuckets, RpsTable};
use super::traceback::{filter_results, RpsTraceScratch, RpsTraceback, RpsTraced};
use super::wordfinder::{RpsDiagTable, RpsWordParams, RpsWordScratch};
use super::{RpsRows, RPS_K_MULT};
use crate::seqdb::rps::RpsVolumeFiles;
use crate::stats::{blast_spouge_stoe, BlastGumbelBlk, KarlinParams};

/// The options of an RPS-BLAST search that blastp's driver does not apply itself: `-max_hsps`
/// and `-qcov_hsp_perc` (the RPS traceback filters with the length of the first query of a
/// batch, `traceback::filter_results`), and the composition-based statistics mode that the
/// command line asks for (the search runs mode 0: the others are refused until item PS11b,
/// unless the database has no frequency ratios, `rps_hooks::check_freq_ratio_file`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RpsSearch {
    pub(crate) max_hsps_per_subject: usize,
    pub(crate) query_cov_hsp_perc: f64,
    pub(crate) requested_cbs: u8,
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
}

/// What the search of every query of a volume reads besides the volume.
pub(crate) struct RpsRun<'a> {
    pub(crate) queries: Vec<RpsQuery<'a>>,
    /// `kbp_gap_std`: BLOSUM62 11/1's table values.
    pub(crate) gapped: KarlinParams,
    pub(crate) gumbel: &'a BlastGumbelBlk,
    /// `avg_subject_length` of the hit cutoff (`BLAST_SEQSRC_MINLENGTH`).
    pub(crate) cutoff_subject_length: i32,
    /// `-evalue` (the preliminary e-value too, without composition-based statistics).
    pub(crate) evalue: f64,
    pub(crate) hitlist_size: usize,
    pub(crate) prelim_hitlist_size: usize,
    pub(crate) window: i32,
    /// `gap_x_dropoff` and `gap_x_dropoff_final` before the scale (from `kbp_gap_std`).
    pub(crate) x_drop_gapped: i32,
    pub(crate) x_drop_final: i32,
    /// blastp's ungapped X-drop and gap trigger in bits (`word_options->x_dropoff`,
    /// `kOptions->gap_trigger`).
    pub(crate) x_drop_ungapped_bits: f64,
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

/// `BLAST_SpougeEtoS` (the hit cutoff of an e-value), the RPS engine's own copy of
/// `stats::blast_spouge_etos` (blastp's default path inlines that one; this search adds no caller
/// to it).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:5236-5268
/// ```c
///     Int4 a=0, b, c;
///     double e;
///     double db_scale_factor = (gbp->db_length) ?
///             (double)gbp->db_length : 1.0;
///
///     b = MAX((int)(log(db_scale_factor/e0) / kbp->Lambda), 2);
///
///     e = BLAST_SpougeStoE(b, kbp, gbp, m, n);
///
///     if (e > e0) {
///         while (e> e0) {
///             a = b;
///             b *= 2;
///             e = BLAST_SpougeStoE(b, kbp, gbp, m, n);
///         }
///     } else {
///         a = 0;
///     }
///     while(b-a > 1) {
///         c = (a+b)/2;
///         e = BLAST_SpougeStoE(c, kbp, gbp, m, n);
///         if (e> e0) {
///             a = c;
///         } else {
///             b = c;
///         }
///     }
///     return a;
/// ```
fn spouge_etos(evalue: f64, kbp: &KarlinParams, gbp: &BlastGumbelBlk, m: i32, n: i32) -> i32 {
    let db_scale_factor = if gbp.db_length != 0 {
        gbp.db_length as f64
    } else {
        1.0
    };
    let e0 = evalue.max(f64::MIN_POSITIVE);
    let mut b = (((db_scale_factor / e0).ln() / kbp.lambda) as i32).max(2);
    let mut a = 0i32;
    let mut e = blast_spouge_stoe(b, kbp, gbp, m, n);
    if e > e0 {
        while e > e0 {
            a = b;
            b *= 2;
            e = blast_spouge_stoe(b, kbp, gbp, m, n);
        }
    }
    while b - a > 1 {
        let c = (a + b) / 2;
        e = blast_spouge_stoe(c, kbp, gbp, m, n);
        if e > e0 {
            a = c;
        } else {
            b = c;
        }
    }
    a
}

impl RpsRun<'_> {
    /// The query's cutoffs, `None` where a scaled one leaves `i32`.
    fn cutoffs(&self, query: &RpsQuery<'_>, scale: f64) -> Option<RpsCutoffs> {
        let int_scale = scale as i32;
        let x_dropoff = scale
            * (self.x_drop_ungapped_bits * std::f64::consts::LN_2 / query.kbp_std.lambda).ceil();
        if !x_dropoff.is_finite() || x_dropoff > f64::from(i32::MAX) {
            return None;
        }
        let hit_cutoff = spouge_etos(
            self.evalue,
            &self.gapped,
            self.gumbel,
            (query.aa_seq.len() - 2) as i32,
            self.cutoff_subject_length,
        )
        .checked_mul(int_scale)?;
        let gap_trigger = ((self.gap_trigger_bits * std::f64::consts::LN_2 + query.kbp_std.k.ln())
            / query.kbp_std.lambda) as i32;
        Some(RpsCutoffs {
            x_dropoff: x_dropoff as i32,
            word_cutoff: gap_trigger.checked_mul(int_scale)?.min(hit_cutoff),
            hit_cutoff,
        })
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:466-469
        // ```c
        //    if (sbp->scale_factor > 1.0) {
        //        ASSERT(Blast_ProgramIsRpsBlast(program_number));
        //        params->gap_x_dropoff *= (Int4)sbp->scale_factor;
        //        params->gap_x_dropoff_final *= (Int4)sbp->scale_factor;
        //    }
        // ```
        let (x_drop_gapped, x_drop_final) = if scale > 1.0 {
            (
                self.x_drop_gapped.checked_mul(int_scale),
                self.x_drop_final.checked_mul(int_scale),
            )
        } else {
            (Some(self.x_drop_gapped), Some(self.x_drop_final))
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:558-562
        // ```c
        //    params->gap_open = score_options->gap_open * (Int4)scale_factor;
        //    params->gap_extend = score_options->gap_extend * (Int4)scale_factor;
        // ```
        let (Some(x_drop_gapped), Some(x_drop_final), Some(gap_open), Some(gap_extend)) = (
            x_drop_gapped,
            x_drop_final,
            files.aux.gap_open.checked_mul(int_scale),
            files.aux.gap_extend.checked_mul(int_scale),
        ) else {
            return Err(scale_overflow(scale));
        };
        // Every searched query's cutoffs, before the search (`BlastHitSavingParametersNew` and
        // `BlastInitialWordParametersNew` set them for the whole batch before its scan).
        let cutoffs = self
            .queries
            .iter()
            .map(|query| {
                if query.valid {
                    self.cutoffs(query, scale)
                        .ok_or_else(|| scale_overflow(scale))
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
            gumbel: self.gumbel,
            expect_value: self.evalue,
            hitlist_size: self.hitlist_size,
        };
        let prelim_karlin = KarlinParams {
            lambda: self.gapped.lambda / scale,
            ..self.gapped
        };
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
                gumbel: self.gumbel,
                prelim_evalue: self.evalue,
                prelim_hitlist_size: self.prelim_hitlist_size,
                query_index: q_idx as u32,
            };
            let subject_seq = &query.aa_seq[1..];
            let Some(lists) = prelim.search_query(&mut scratch.0, subject_seq)? else {
                return Ok(none);
            };
            let letters = &query.aa_seq[1..query.aa_seq.len() - 1];
            let mut traced = traceback.trace_query(lists, letters, q_idx as u32, &mut scratch.1)?;
            if let Some(hit_list) = traced.hit_list.as_mut() {
                filter_results(
                    hit_list,
                    self.options.max_hsps_per_subject,
                    self.options.query_cov_hsp_perc,
                    query.first_query_length,
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
/// `RPS_K_MULT` times the `.aux` K of the highest profile whose rescaled PSSM any query of the
/// batch wrote (the last list the batch's traceback reads, by ordinal id), else the table's K;
/// every other query prints the table's K.
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
) -> Vec<f64> {
    let mut gapped_k = vec![table_k; k_marks.len()];
    for batch in batches {
        if batch.is_empty() {
            continue;
        }
        if let Some(mark) = k_marks[batch.clone()].iter().flatten().max() {
            gapped_k[batch.start] = RPS_K_MULT * karlin_k[*mark as usize];
        }
    }
    gapped_k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_first_query_of_a_batch_gets_the_k_of_the_highest_mark() {
        let marks = [Some(2), None, Some(5), None, None, Some(1)];
        let karlin_k = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06];
        let k = batch_gapped_k(&[0..3, 3..5, 5..6], &marks, &karlin_k, 0.041);
        assert_eq!(k[0], RPS_K_MULT * 0.06);
        assert_eq!(k[1], 0.041);
        assert_eq!(k[2], 0.041);
        assert_eq!(k[3], 0.041);
        assert_eq!(k[5], RPS_K_MULT * 0.02);
    }
}
