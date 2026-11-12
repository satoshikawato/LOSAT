//! The RPS-BLAST traceback with composition-based statistics (`-comp_based_stats` 1, 2 and 3):
//! `s_RPSComputeTraceback` hands each preliminary list to `Blast_RedoAlignmentCore` with the
//! profile's consensus sequence as the query, the user's query as the subject, the profile's
//! frequency ratios (`.freq`) and a zero PSSM that the composition adjustment fills, then
//! switches the roles back (`s_BlastHSPListRPSUpdate`). The redo is PS4b's
//! (`pssm_search::redo1::RedoMatch`), with the E-value block of this search: Lambda over the
//! scaling factor 32 and K 1.2 times the profile's `.aux` K.
//!
//! Modes 2 and 3 run mode 1's adjustment (a position-based search "can only use traditional
//! composition based stats"); their difference from mode 1 is the preliminary stage's E-value
//! threshold (`engine.rs`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1238-1260
//! ```c
//!       if (ext_params->options->compositionBasedStats > 0) {
//!           Int4 offsets[2], query_index;
//!           s_RPSFillFreqRatiosInPsiMatrix(sbp->psi_matrix, rpsblast_freqs + db_seq_start, seq_arg.seq->length);
//!           one_query_info->first_context = one_query_info->last_context = 0;
//!           one_query_info->num_queries = 1;
//!           offsets[0] = 0;
//!           offsets[1] = seq_arg.seq->length + 1;
//!           one_query_info->max_length = seq_arg.seq->length;
//!           OffsetArrayToContextOffsets(one_query_info, offsets, program_number);
//!           query_index = hsp_list->query_index;
//!           hsp_list->query_index = 0;
//!           Blast_RedoAlignmentCore(program_number, seq_arg.seq, one_query_info, sbp,
//!                                   one_query, NULL, default_db_genetic_code,
//!                                   hsp_list, NULL, score_params, ext_params,
//!                                   hit_params, psi_options, NULL);
//!           hsp_list->query_index = query_index;
//!           s_BlastHSPListRPSUpdate(program_number, hsp_list);
//! ```

use anyhow::Result;

use super::traceback::{sort_by_score, swap_roles, RpsTraceScratch, RpsTraceback, RpsTraced};
use super::RPS_K_MULT;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::pssm_search::evalue::EvalueStats;
use crate::algorithm::pssm_search::redo1::{RedoMatch, RedoRun, SCALING_FACTOR};
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::core::composition_adjustment::compo_scale::PosMatrixInfo;
use crate::core::psi::PsiScoreBlk;
use crate::seqdb::rps::RpsVolumeFiles;
use crate::stats::KarlinParams;

/// What the redo of every profile reads besides the traceback's values: the volume's frequency
/// ratios, BLOSUM62's score block (the matrix information of the position-based adjustment),
/// the `.aux` gap costs before the scale and the final X-drop in bits
/// (`options->gap_x_dropoff_final`: `-xdrop_gap_final`, or 25).
pub(crate) struct RpsRedoParams<'a> {
    pub(crate) files: &'a RpsVolumeFiles,
    pub(crate) sbp: &'a PsiScoreBlk,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    pub(crate) x_drop_final_bits: f64,
}

/// The statistics the traceback of one query reads: the batch's ungapped Karlin block and the
/// query's E-value statistics in this volume.
#[derive(Clone, Copy)]
pub(super) struct RpsQueryStats<'a> {
    pub(super) batch_kbp: &'a KarlinParams,
    pub(super) stats: EvalueStats<'a>,
}

impl RpsTraceback<'_> {
    /// The redo of the query of letters `query` (without sentinels) from its preliminary hit
    /// list. NCBI reads the lists in the order of the stream: sorted by best E-value, then the
    /// first HSP's score (down), then the ordinal id (down) (`Blast_HSPResultsReverseSort`
    /// reads from the end, `s_EvalueCompareHSPLists`); each list writes the profile's K
    /// (`RPS_K_MULT` times the `.aux` K) before its redo, so the last list read leaves the K
    /// that the query's batch may print (`k_mark`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1203-1208
    /// ```c
    ///          /* for CBS only allocate memory for PSSM matrix */
    ///          if (ext_params->options->compositionBasedStats > 0) {
    ///              sbp->psi_matrix->pssm->data = (Int4**)_PSIAllocateMatrix(
    ///                                                           seq_arg.seq->length,
    ///                                                           BLASTAA_SIZE,
    ///                                                           sizeof(Int4));
    /// ```
    /// (`_PSIAllocateMatrix` zeroes the matrix.) The Karlin block of the redo is
    /// `sbp->kbp_gap[0]` of the score block copy (`kbp_gap_std`, whose K was just written),
    /// its Lambda over the scaling factor (`s_RescaleSearch`); the first valid query of the
    /// batch gives the ungapped block that stays where the profile has no score frequencies
    /// (`kbp_psi[0]`).
    pub(super) fn trace_query_redo(
        &self,
        redo: &RpsRedoParams<'_>,
        mut prelim: BlastpHitList<XPrelimHsp>,
        query: &[u8],
        q_idx: u32,
        statistics: RpsQueryStats<'_>,
        scratch: &mut RpsTraceScratch,
    ) -> Result<RpsTraced> {
        let RpsQueryStats { batch_kbp, stats } = statistics;
        prelim.sort_by_evalue();
        let count = prelim.hsplist_count;
        let mut hit_list: Option<BlastpHitList> = None;
        let mut k_mark: Option<u32> = None;
        let kbp_psi = crate::core::psi::KarlinBlk {
            lambda: batch_kbp.lambda,
            k: batch_kbp.k,
            h: batch_kbp.h,
        };
        for list in prelim.hsplist_array.into_iter().take(count) {
            let oid = list.oid as usize;
            let start = self.offsets[oid] as usize;
            let length = (self.offsets[oid + 1] - self.offsets[oid] - 1) as usize;
            k_mark = Some(list.oid);
            if list.hsps.is_empty() {
                continue;
            }
            let consensus = &self.consensus[start..start + length];
            redo.files
                .freq_ratios(start, length, &mut scratch.freq_bytes, &mut scratch.ratios)
                .map_err(|error| {
                    anyhow::anyhow!("a profile database whose .freq file cannot be read ({error}) is not supported by LOSAT's RPSBLAST")
                })?;
            let Some(info) = PosMatrixInfo::new(
                consensus,
                &scratch.ratios,
                redo.sbp,
                SCALING_FACTOR,
                kbp_psi,
            ) else {
                anyhow::bail!(
                    "a profile whose position-based scaling for the composition-based statistics gives no Karlin-Altschul parameters (NCBI BLAST+ 2.17.0 goes on with a Lambda of -1) is not supported by LOSAT's RPSBLAST"
                );
            };
            // The best scores of the scaled start matrix's rows sum to the most an alignment can
            // score; above `INT4_MAX` NCBI's 32-bit scores overflow (undefined behavior).
            let bound: i64 = info
                .start_matrix
                .iter()
                .map(|row| i64::from(row.iter().copied().max().unwrap_or(0).max(0)))
                .sum();
            if bound > i64::from(i32::MAX) {
                anyhow::bail!(
                    "a profile whose scaled scores can add up to {bound} (above 2147483647, where NCBI BLAST+'s 32-bit scores overflow) is not supported by LOSAT's RPSBLAST"
                );
            }
            // The zero PSSM of the redo's copy.
            let mut matrix = QueryPssm {
                rows: std::mem::take(&mut scratch.pssm.rows),
            };
            matrix.rows.clear();
            matrix.rows.resize(length, [0; 28]);
            let k = RPS_K_MULT * self.karlin_k[oid];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2381-2383
            // ```c
            //     gapping_params->x_dropoff = (Int4)
            //         MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
            //             extendParams->gap_x_dropoff_final);
            // ```
            // `options->gap_x_dropoff_final` is the search's final X-drop in bits
            // (`-xdrop_gap_final`; `RpsRun::gapped_x_drops` has checked the conversion).
            let mut run = RedoRun::new(
                consensus,
                consensus,
                KarlinParams {
                    lambda: self.gapped.lambda / SCALING_FACTOR,
                    k,
                    ..self.gapped
                },
                k.ln(),
                (redo.gap_open, redo.gap_extend),
                self.x_drop_final,
                matrix,
                std::mem::take(&mut scratch.gap),
            )
            .with_final_x_drop_bits(redo.x_drop_final_bits, self.x_drop_final)
            .with_bound_check();
            let hsps: Vec<BlastpPreliminaryHsp> = list.hsps.iter().map(|hsp| hsp.hsp).collect();
            let redone = RedoMatch {
                info: &info,
                // The query's statistics in this volume: the Gumbel block, or without it
                // (`-searchsp`, `OLD_FSC`) the search space the redo copy of the score block
                // keeps (the user query's, blast_traceback.c:1238-1255).
                stats,
                expect_value: self.expect_value,
                query_length: length as i32,
                rps_roles: true,
            }
            .redo_match(&mut run, list.oid, &hsps, query);
            let exceeded = run.bound_exceeded();
            let (matrix, gap) = run.into_buffers();
            // The bound of every matrix the adjustments built (`RedoRun::with_bound_check`).
            if let Some(bound) = exceeded {
                anyhow::bail!(
                    "a profile whose composition-adjusted scores can add up to {bound} (above 2147483647, where NCBI BLAST+'s 32-bit scores overflow) is not supported by LOSAT's RPSBLAST"
                );
            }
            scratch.pssm.rows = matrix.rows;
            scratch.gap = gap;
            // The list's HSPs are switched back and sorted even if none is left, and a list
            // without HSPs is dropped (blast_traceback.c:1272-1275).
            let Some((mut hsp_list, _, _)) = redone else {
                continue;
            };
            hsp_list.query_index = q_idx;
            swap_roles(&mut hsp_list, |_| query.len());
            sort_by_score(&mut hsp_list.hsps);
            if hsp_list.hsps.is_empty() {
                continue;
            }
            hit_list
                .get_or_insert_with(|| BlastpHitList::new(self.hitlist_size))
                .update(hsp_list);
        }
        if let Some(hit_list) = hit_list.as_mut() {
            hit_list.sort_by_evalue();
        }
        Ok(RpsTraced { hit_list, k_mark })
    }
}
