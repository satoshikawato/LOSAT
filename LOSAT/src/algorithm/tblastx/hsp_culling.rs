//! NCBI BLAST HSP culling (`-culling_limit`) for TBLASTX.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c
//!
//! The interval tree itself is the shared `core::hspfilter_culling::CullingForest`; this
//! module holds what is TBLASTX's own: the context of an HSP (its frame), the length of
//! each context and the keys that `s_BlastHSPCullingRun` reads from a TBLASTX HSP.
//!
//! NCBI culls twice when `-culling_limit` is given (`report::culled_hit_order`): a writer
//! in the preliminary stage and a pipe after the traceback.

use super::blast_engine::TblastxHsp;
use crate::core::hspfilter_culling::{CullKey, CullingForest};

/// The query context of a TBLASTX HSP within its query: NCBI's contexts of a translated
/// query are its frames 1, 2, 3, -1, -2, -3 (`BLAST_GetAllTranslations`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1211-1220
pub(crate) fn tblastx_context(frame: i32) -> usize {
    if frame > 0 {
        (frame - 1) as usize
    } else {
        (2 - frame) as usize
    }
}

/// The length in residues of each frame of a query of `length` nucleotides, as NCBI's
/// query information holds it (`contexts[context].query_length`), in context order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:193-195
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:922-929
/// These are the lengths of a query searched on both strands (`eNa_strand_both`). With
/// `-strand plus` or `minus` NCBI gives the contexts of the other strand length 0
/// (blast_setup_cxx.cpp:200-215), but those contexts hold no HSP, and the culling reads the
/// length of an HSP's own context only (`CullingWriter::run`), which is the same.
pub(crate) fn tblastx_context_lengths(nucleotide_length: usize) -> [i32; 6] {
    let mut lengths = [0; 6];
    for (context, slot) in lengths.iter_mut().enumerate() {
        let shift = context % 3;
        *slot = if nucleotide_length == 0 || nucleotide_length <= shift {
            0
        } else {
            i32::try_from((nucleotide_length - shift) / 3).unwrap_or(i32::MAX)
        };
    }
    lengths
}

/// NCBI's culling writer (or pipe) over the HSPs of one query: `s_BlastHSPCullingInit`
/// (`new`), `s_BlastHSPCullingRun` (`run`) and `s_BlastHSPCullingFinal` (`finalize`).
pub(crate) struct CullingWriter {
    context_lengths: [i32; 6],
    forest: CullingForest<TblastxHsp>,
}

impl CullingWriter {
    pub(crate) fn new(culling_max: i64, context_lengths: [i32; 6]) -> Self {
        Self {
            context_lengths,
            forest: CullingForest::new(culling_max),
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:602-644
    /// TBLASTX is not BLASTN: the context is the HSP's own and the offsets are the frame
    /// offsets of `Blast_HSPInit` (`sort_query_offset`, `sort_query_end`). The HSPs of
    /// `hsp_list` come in its order; an HSP that is not saved is freed.
    pub(crate) fn run(&mut self, hsp_list: Vec<TblastxHsp>, oid: u32) {
        for hsp in hsp_list {
            let cid = tblastx_context(hsp.hit.query_frame);
            let qlen = self.context_lengths[cid];
            let key = CullKey {
                begin: hsp.hit.sort_query_offset as i32,
                end: hsp.hit.sort_query_end as i32,
                score: hsp.hit.raw_score,
                subject_offset: hsp.hit.sort_subject_offset as i64,
            };
            self.forest.save(cid, qlen, oid, key, hsp);
        }
    }

    /// The HSPs that the trees keep, as `s_BlastHSPCullingFinal` puts them into the hit
    /// list of the query (`CullingForest::finalize`; all six contexts are one query).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:516-581
    pub(crate) fn finalize(self) -> Vec<(u32, Vec<TblastxHsp>)> {
        self.forest
            .finalize(
                |_| 0,
                |a, b| crate::common::score_compare_hsps(&a.hit, &b.hit),
            )
            .into_iter()
            .next()
            .map(|(_, lists)| lists)
            .unwrap_or_default()
    }
}

#[cfg(test)]
#[path = "hsp_culling_tests.rs"]
mod tests;
