//! NCBI BLAST HSP culling (`-culling_limit`) for TBLASTX.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/hspfilter_culling.c
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
/// NCBI reference: c++/src/algo/blast/core/blast_util.c:1211-1220
/// ```c
/// Int4 BLAST_FrameToContext(Int2 frame, EBlastProgramType program)
/// {
///     if (Blast_QueryIsTranslated(program) ||
///         Blast_SubjectIsTranslated(program)) {
///         ASSERT(frame >= -3 && frame <= 3 && frame != 0);
///         if (frame > 0) {
///             return frame - 1;
///         } else {
///             return 2 - frame;
///         }
/// ```
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
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:193-195
/// ```c
///             for (unsigned int i = 0; i < kNumContexts; i++) {
///                 unsigned int prot_length =
///                     static_cast<unsigned int>(BLAST_GetTranslatedProteinLength(length, i));
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_util.c:922-929
/// ```c
/// size_t
/// BLAST_GetTranslatedProteinLength(size_t nucleotide_length, unsigned int context)
/// {
///     if (nucleotide_length == 0 || nucleotide_length <= context % CODON_LENGTH) {
///         return 0;
///     }
///     return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
/// }
/// ```
/// TBLASTX's queries are searched on both strands (`eNa_strand_both`), so every context
/// has this length.
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
    pub(crate) fn new(culling_max: i32, context_lengths: [i32; 6]) -> Self {
        Self {
            context_lengths,
            forest: CullingForest::new(culling_max),
        }
    }

    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:602-644
    /// ```c
    ///    for (i=0; i<hsp_list->hspcnt; ++i) {
    ///       /* wrap the hsp with a LinkedHSP structure */
    ///       A.hsp   = hsp_list->hsp_array[i];
    ///       A.cid   = isBlastn ? (A.hsp->context  - A.hsp->context % NUM_STRANDS) : A.hsp->context;
    ///       A.sid   = hsp_list->oid;
    ///       A.merit = params->culling_max;
    ///       qlen    = cull_data->query_info->contexts[A.hsp->context].query_length;
    ///       if(isBlastn && (A.hsp->context % NUM_STRANDS)) {
    ///     	  A.begin = qlen - A.hsp->query.end;
    ///     	  A.end   = qlen -  A.hsp->query.offset;
    ///       }
    ///       else {
    ///     	  A.begin = A.hsp->query.offset;
    ///     	  A.end   = A.hsp->query.end;
    ///       }
    ///       A.next  = NULL;
    ///
    ///       if (! c_tree[A.cid]) {
    ///          c_tree[A.cid] = s_CTreeNew(qlen);
    ///       }
    ///
    ///       if(s_SaveHSP(c_tree[A.cid], &A)){
    ///     	 hsp_list->hsp_array[i] = NULL;
    ///       }
    ///    }
    /// ```
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
    /// NCBI reference: c++/src/algo/blast/core/hspfilter_culling.c:516-581
    /// ```c
    ///    for (cid=0; cid < cull_data->num_contexts; ++cid) {
    ///       if (c_tree[cid]) {
    /// ...
    ///          cull_list = s_RipHSPOffCTree(c_tree[cid]);
    /// ...
    ///          while (cull_list) {
    ///             p = cull_list;
    ///             /* test to see if new hsplist has already been allocated */
    ///             allocated = FALSE;
    ///             for (sid=0; sid<hitlist->hsplist_count; ++sid) {
    ///                 list = hitlist->hsplist_array[sid];
    ///                 if (p->sid == list->oid) {
    ///                    allocated = TRUE;
    ///                    break;
    ///                 }
    ///             }
    /// ...
    ///             list->hsp_array[id] = p->hsp;
    ///             list->hspcnt++;
    /// ...
    ///          for (sid=0; sid < hitlist->hsplist_count; ++sid) {
    ///             list = hitlist->hsplist_array[sid];
    ///             best_evalue = (double) INT4_MAX;
    ///             for (id=0; id < list->hspcnt; ++id) {
    ///                 best_evalue = MIN(list->hsp_array[id]->evalue, best_evalue);
    ///             }
    ///             Blast_HSPListSortByScore(list);
    /// ```
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
mod tests {
    use super::*;
    use crate::common::Hit;

    fn hsp(oid: u32, frame: i32, begin: usize, end: usize, score: i32) -> TblastxHsp {
        TblastxHsp {
            hit: Hit {
                identity: 0.0,
                length: end - begin,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 3,
                s_start: 1,
                s_end: 3,
                e_value: 1e-5,
                bit_score: 0.0,
                num_ident: 0,
                query_frame: frame,
                query_length: 0,
                q_idx: 0,
                s_idx: oid,
                raw_score: score,
                sort_query_offset: begin,
                sort_query_end: end,
                sort_subject_offset: begin,
                sort_subject_end: end,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 0,
            },
            subject_frame: 1,
            num: 1,
        }
    }

    fn kept(writer: CullingWriter) -> Vec<(u32, i32)> {
        writer
            .finalize()
            .into_iter()
            .flat_map(|(oid, list)| list.into_iter().map(move |h| (oid, h.hit.raw_score)))
            .collect()
    }

    // NCBI c++/src/algo/blast/core/blast_util.c:922-929,1211-1220: frames 1, 2, 3, -1,
    // -2, -3 are the contexts 0..5, of (length - context % 3) / 3 residues.
    #[test]
    fn contexts_follow_ncbis_frames() {
        let order: Vec<usize> = [1, 2, 3, -1, -2, -3]
            .iter()
            .map(|&frame| tblastx_context(frame))
            .collect();
        assert_eq!(order, vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(tblastx_context_lengths(31), [10, 10, 9, 10, 10, 9]);
        assert_eq!(tblastx_context_lengths(0), [0; 6]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:136-163: an HSP never dominates
    // itself (`r != y`), so a limit of 1 keeps the best HSP of overlapping ones.
    #[test]
    fn a_limit_of_one_keeps_the_dominating_hsp() {
        let mut writer = CullingWriter::new(1, tblastx_context_lengths(300));
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
        assert_eq!(kept(writer), vec![(0, 50)]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:79-120: with equal score and range
    // the lower OID dominates; other contexts have their own trees.
    #[test]
    fn ties_are_broken_by_oid_and_contexts_do_not_interact() {
        let mut writer = CullingWriter::new(1, tblastx_context_lengths(300));
        writer.run(vec![hsp(1, 1, 10, 40, 50)], 1);
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, -1, 10, 40, 50)], 0);
        assert_eq!(kept(writer), vec![(0, 50), (0, 50)]);
    }

    // The merit is an Int4 that NCBI's binary wraps: a negative merit loses an HSP at its
    // first dominator, and i32::MIN never reaches 0 after one decrement.
    #[test]
    fn merits_wrap_as_ncbis_int4() {
        let lengths = tblastx_context_lengths(300);
        let mut writer = CullingWriter::new(2147483646i32.wrapping_add(3), lengths);
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
        assert_eq!(kept(writer), vec![(0, 50)]);
        let mut writer = CullingWriter::new(2147483645i32.wrapping_add(3), lengths);
        writer.run(vec![hsp(0, 1, 10, 40, 50), hsp(0, 1, 12, 40, 30)], 0);
        assert_eq!(kept(writer), vec![(0, 50), (0, 30)]);
    }

    // NCBI c++/src/algo/blast/core/hspfilter_culling.c:258-299,453-458: a leaf with 20
    // HSPs forks; HSPs that do not overlap are all kept.
    #[test]
    fn many_disjoint_hsps_fork_and_survive() {
        let mut writer = CullingWriter::new(1, tblastx_context_lengths(3000));
        let list: Vec<TblastxHsp> = (0..40)
            .map(|i| hsp(0, 1, i * 20, i * 20 + 10, 30))
            .collect();
        writer.run(list, 0);
        assert_eq!(kept(writer).len(), 40);
    }
}
