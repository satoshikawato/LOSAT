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
//! rpstblastn's redo is the same with the user query's six frames as the subject (series PR,
//! item PS12b): NCBI's `one_query` is the query's translation buffer, the redo's subject is
//! translated (`subject_is_translated`: windows in the frames, `RedoMatch::redo_windows` with a
//! `FrameSubject`), its subject length for the e-values is `s_GetSubjectLength`'s and every list
//! writes `kbp_gap[0]`'s K, the one the query's batch may print. With `-use_sw_tback` both
//! programs' redo is the Smith-Waterman redo of the same windows (item PS12c for rpstblastn).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1238-1260

use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::Result;

use super::traceback::{sort_by_score, swap_roles, RpsTraceScratch, RpsTraceback, RpsTraced};
use super::RPS_K_MULT;
use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
use crate::algorithm::blastp::hsp::{BlastCompoHeap, BlastpHitList, BlastpHspList};
use crate::algorithm::blastp::x_prelim_compact::XPrelimHsp;
use crate::algorithm::pssm_search::evalue::EvalueStats;
use crate::algorithm::pssm_search::redo1::smith_waterman::SwEndless;
use crate::algorithm::pssm_search::redo1::{RedoMatch, RedoRun, RedoSubject, SCALING_FACTOR};
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::common::GapEditOp;
use crate::core::composition_adjustment::compo_scale::{PosMatrixInfo, RedoInstance};
use crate::core::psi::PsiScoreBlk;
use crate::seqdb::rps::RpsVolumeFiles;
use crate::stats::KarlinParams;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:163
/// The inclusion threshold of rpsblast's PSI options (`PSIBlastOptionsNew`), the e-value cutoff of
/// the redo's heaps (blast_kappa.c:3074-3076).
const PSI_INCLUSION_ETHRESH: f64 = 0.002;

/// The score block of the redo's matrix information (`RpsRedoParams::sbp`) of RPS-BLAST and
/// RPS-TBLASTN: the matrix's frequency ratios with the search's gapped block `gapped`
/// (`in_pssm::score_block`, psiblast's), made once per run before any batch.
pub(crate) fn redo_score_block(
    matrix: crate::config::ScoringMatrix,
    gapped: crate::core::psi::KarlinBlk,
) -> Result<PsiScoreBlk> {
    crate::algorithm::pssm_search::in_pssm::score_block(matrix, gapped)
}

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
    /// `-use_sw_tback`: NCBI's Smith-Waterman redo of each profile
    /// (`RedoSubject::redo_smith_waterman`: `RedoMatch::redo_match_smith_waterman` of rpsblast's
    /// query, `RedoMatch::redo_windows_smith_waterman` of rpstblastn's frames). The redo runs
    /// with composition-based statistics 1 to 3 only; without them NCBI's RPS traceback does not
    /// read the flag
    /// (`s_RPSComputeTraceback` takes no Smith-Waterman path and `kSmithWaterman` asks for the
    /// other mode).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:291-292
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1210-1214
    pub(crate) smith_waterman: bool,
    /// The program of LOSAT's refusals (`RPSBLAST`, `RPSTBLASTN`).
    pub(crate) program: &'static str,
}

/// The statistics the traceback of one query reads: the batch's ungapped Karlin block and the
/// query's E-value statistics in this volume.
#[derive(Clone, Copy)]
pub(super) struct RpsQueryStats<'a> {
    pub(super) batch_kbp: &'a KarlinParams,
    pub(super) stats: EvalueStats<'a>,
    /// Whether the user query's context 0 (frame +1; a protein query's only one) is valid: the
    /// only context the redo's query information keeps (`RedoRun::without_valid_context`).
    pub(super) context_0_valid: bool,
}

impl RpsTraceback<'_> {
    /// The redo of the query `query` (a protein query's letters without sentinels, or a
    /// translated query's translation buffer with its frames in `kind`) from its preliminary
    /// lists `lists`, in the order NCBI's stream reads them (`stream_lists`): a protein query's
    /// sorted by best E-value, then the first HSP's score (down), then the ordinal id (down)
    /// (`Blast_HSPResultsReverseSort` reads from the end, `s_EvalueCompareHSPLists`), a
    /// translated query's by ordinal id. Each list writes the profile's K (`RPS_K_MULT` times
    /// the `.aux` K) before its redo, so the last list read leaves the K that the query's batch
    /// may print (`k_mark`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1202-1208
    /// (`_PSIAllocateMatrix` zeroes the matrix.) The Karlin block of the redo is
    /// `sbp->kbp_gap[0]` of the score block copy (`kbp_gap_std`, whose K was just written),
    /// its Lambda over the scaling factor (`s_RescaleSearch`); the first valid query of the
    /// batch gives the ungapped block that stays where the profile has no score frequencies
    /// (`kbp_psi[0]`).
    pub(super) fn trace_query_redo<'q, S: RedoSubject<'q>>(
        &self,
        redo: &RpsRedoParams<'_>,
        lists: Vec<BlastpHspList<XPrelimHsp>>,
        query: &'q [u8],
        kind: S,
        q_idx: u32,
        statistics: RpsQueryStats<'_>,
        scratch: &mut RpsTraceScratch,
    ) -> Result<RpsTraced> {
        let RpsQueryStats {
            batch_kbp,
            stats,
            context_0_valid,
        } = statistics;
        let program = redo.program;
        let mut hit_list: Option<BlastpHitList> = None;
        let mut k_mark: Option<u32> = None;
        // The heap of redone matches that the Smith-Waterman redo's first alignment asks
        // (`BlastCompo_HeapWouldInsert`). NCBI's redo of each profile makes its own
        // (`BlastCompo_HeapInitialize` with the hit list size and the inclusion threshold) and
        // never inserts into it: RPS-BLAST's redo has no sequence source and leaves the match
        // before the insertion, so the heap stays empty and asks only whether `n` is below
        // `heapThreshold`.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3288-3292
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3702-3709
        let mut sw_heap: Option<BlastCompoHeap> = None;
        let kbp_psi = crate::core::psi::KarlinBlk {
            lambda: batch_kbp.lambda,
            k: batch_kbp.k,
            h: batch_kbp.h,
        };
        for list in lists {
            let oid = list.oid as usize;
            let start = self.offsets[oid] as usize;
            let length = (self.offsets[oid + 1] - self.offsets[oid] - 1) as usize;
            k_mark = Some(list.oid);
            if list.hsps.is_empty() {
                continue;
            }
            let consensus = &self.consensus.letters()?[start..start + length];
            redo.files
                .freq_ratios(start, length, &mut scratch.freq_bytes, &mut scratch.ratios)
                .map_err(|error| {
                    anyhow::anyhow!("a profile database whose .freq file cannot be read ({error}) is not supported by LOSAT's {program}")
                })?;
            let Some(info) = PosMatrixInfo::new(
                consensus,
                &scratch.ratios,
                redo.sbp,
                SCALING_FACTOR,
                kbp_psi,
            ) else {
                anyhow::bail!(
                    "a profile whose position-based scaling for the composition-based statistics gives no Karlin-Altschul parameters (NCBI BLAST+ 2.17.0 goes on with a Lambda of -1) is not supported by LOSAT's {program}"
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
                    "a profile whose scaled scores can add up to {bound} (above 2147483647, where NCBI BLAST+'s 32-bit scores overflow) is not supported by LOSAT's {program}"
                );
            }
            // The zero PSSM of the redo's copy.
            let mut matrix = QueryPssm {
                rows: std::mem::take(&mut scratch.pssm.rows),
            };
            matrix.rows.clear();
            matrix.rows.resize(length, [0; 28]);
            let k = RPS_K_MULT * self.karlin_k[oid];
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2380-2383
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
                k.glibc_ln(),
                (redo.gap_open, redo.gap_extend),
                (redo.x_drop_final_bits, self.x_drop_final),
                matrix,
                std::mem::take(&mut scratch.gap),
            )
            .with_bound_check();
            if !context_0_valid {
                run = run.without_valid_context();
            }
            let hsps: Vec<BlastpPreliminaryHsp> = list.hsps.iter().map(|hsp| hsp.hsp).collect();
            let redo_match = RedoMatch {
                info: &info,
                // The query's statistics in this volume: the Gumbel block, or without it
                // (`-searchsp`, `OLD_FSC`) the search space the redo copy of the score block
                // keeps (the user query's, blast_traceback.c:1238-1255). The Smith-Waterman
                // redo's significance reads that search space also with the Gumbel block
                // (`eff_search_space`, redo_alignment.c:1406).
                stats,
                expect_value: self.expect_value,
                query_length: length as i32,
                rps_roles: true,
                // rpsblast's and rpstblastn's databases are BLOSUM62's (the others are refused).
                matrix: crate::config::ScoringMatrix::Blosum62,
                // rpsblast and rpstblastn refuse `-sum_stats true`.
                link: None,
            };
            // The Smith-Waterman redo of a protein query (rpsblast) or of a translated query's
            // frames (rpstblastn, series PR, item PS12c).
            let redone = if redo.smith_waterman {
                let heap = sw_heap.get_or_insert_with(|| {
                    BlastCompoHeap::new(self.hitlist_size, PSI_INCLUSION_ETHRESH)
                });
                kind.redo_smith_waterman(
                    &redo_match,
                    &mut run,
                    &mut scratch.sw,
                    heap,
                    list.oid,
                    &hsps,
                    query,
                )
            } else {
                Ok(kind.redo(&redo_match, &mut run, list.oid, &hsps, query))
            };
            let exceeded = run.bound_exceeded();
            let (matrix, gap) = run.into_buffers();
            // The bound of every matrix the adjustments built (`RedoRun::with_bound_check`).
            if let Some(bound) = exceeded {
                anyhow::bail!(
                    "a profile whose composition-adjusted scores can add up to {bound} (above 2147483647, where NCBI BLAST+'s 32-bit scores overflow) is not supported by LOSAT's {program}"
                );
            }
            let redone = redone.map_err(|SwEndless| {
                if S::Instance::SUBJECT_IS_TRANSLATED {
                    anyhow::anyhow!(
                        "a profile whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0 (the e-value of a score of 0 passes the threshold, the profile has several preliminary HSPs in a window of a query frame and its first position scores the window's first letter above 0) is not supported by LOSAT's {program}"
                    )
                } else {
                    anyhow::anyhow!(
                        "a profile whose Smith-Waterman redo (-use_sw_tback) does not end in NCBI BLAST+ 2.17.0 (the e-value of a score of 0 passes the threshold, the profile has several preliminary HSPs and its first position scores the query's first letter above 0) is not supported by LOSAT's {program}"
                    )
                }
            })?;
            scratch.pssm.rows = matrix.rows;
            scratch.gap = gap;
            // The list's HSPs are switched back and sorted even if none is left, and a list
            // without HSPs is dropped (blast_traceback.c:1272-1275).
            let Some((mut hsp_list, _, _)) = redone else {
                continue;
            };
            hsp_list.query_index = q_idx;
            // A frame's HSPs take the frame's length (a protein query's: its length).
            swap_roles(&mut hsp_list, |frame| {
                kind.frame_letters(query, frame).len()
            });
            sort_by_score(&mut hsp_list.hsps);
            if hsp_list.hsps.is_empty() {
                continue;
            }
            // An alignment of the Smith-Waterman redo of a translated query that starts with a
            // gap in the query's row (series PR, item PS12d): no search was found that gives one
            // (`PR/ps12d/scripts/qgap_grid.sh`: every such start was in the profile's row, whose
            // reports LOSAT writes as NCBI's, `serial::seqalign::std_seg_row_range`), so nothing
            // shows what NCBI's reports write for it, and LOSAT refuses the profile.
            if S::Instance::SUBJECT_IS_TRANSLATED
                && redo.smith_waterman
                && hsp_list.hsps.iter().any(starts_with_query_gap)
            {
                anyhow::bail!(
                    "a profile whose Smith-Waterman redo (-use_sw_tback) gives an alignment that starts with a gap in the query's row, for which no NCBI BLAST+ 2.17.0 report was found to follow, is not supported by LOSAT's {program}"
                );
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

/// Whether the alignment `hsp` (switched back: its query is the search's query) starts with a gap
/// in the query's row (`serial::seqalign::gapped_row_range`).
pub(super) fn starts_with_query_gap(hsp: &crate::algorithm::blastp::hsp::BlastpHsp) -> bool {
    hsp.gap_info.as_deref().is_some_and(|ops: &[GapEditOp]| {
        crate::report::serial::seqalign::gapped_row_range(
            ops,
            crate::report::serial::seqalign::AlignRow::Query,
        )
        .is_some()
    })
}
