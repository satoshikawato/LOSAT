//! Ungapped HSP record of TBLASTX
//!
//! `UngappedHit` is the lightweight HSP the scan, the linking and the report steps pass around.

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-101,126-128,133,158
// ```c
// typedef struct BlastSeg {
//    Int2 frame;  /**< Translation frame */
//    Int4 offset; /**< Start of hsp */
//    Int4 end;    /**< End of hsp */
//    Int4 gapped_start;/**< Where the gapped extension started. */
// } BlastSeg;
// typedef struct BlastHSP {
//    Int4 score;           /**< This HSP's raw score */
//    Int4 num_ident;       /**< Number of identical base pairs in this HSP */
//    Int4     context;     /**< Context number of query */
//    Int4 hspcnt; /**< Number of HSPs saved */
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:246
// ```c
//    Int4 length;         /**< Length of sequence. */
// ```
// No NCBI counterpart: the field widths. NCBI holds an HSP's offsets, seed offsets, context,
// identity count and its position in a list (`hspcnt`), and a sequence's length, in `Int4`, so
// every value NCBI accepts there fits `u32` (non-negative values below 2^31). The record is
// built only with checked conversions (`u32::try_from(..).expect(..)`, `u8::try_from` for the
// subject frame index 0..5), and every use widens the field back to the type its expression
// computed in before (`usize` for offsets and indexes), so values, order and output are the
// same as with `usize` fields; the record is 88 bytes instead of 152 on 64-bit targets.
/// Lightweight ungapped HSP used in the *hot path* before final output formatting.
///
/// This avoids per-hit `String` cloning and expensive identity/bit-score/e-value
/// computation while scanning/extension is running.
#[derive(Debug, Clone)]
pub struct UngappedHit {
    /// Query record index (into `query_ids`)
    pub q_idx: u32,
    /// Subject record index (into `subject_ids`)
    pub s_idx: u32,
    /// Query context index (into `tblastx::lookup::QueryContext` vector)
    pub ctx_idx: u32,
    /// Subject frame index within `generate_frames()` output for this subject (0..5)
    pub s_f_idx: u8,
    pub q_frame: i8,
    pub s_frame: i8,
    pub q_aa_start: u32,
    pub q_aa_end: u32,
    pub s_aa_start: u32,
    pub s_aa_end: u32,
    /// NCBI `BlastHSP::query.gapped_start` copied from
    /// `init_hsp->offsets.qs_offsets.q_off`.
    /// Reference: blast_gapalign.c:4756-4768, blast_extend.c:351-370.
    pub q_seed_off: u32,
    /// NCBI `BlastHSP::subject.gapped_start` copied from
    /// `init_hsp->offsets.qs_offsets.s_off`.
    /// Reference: blast_gapalign.c:4756-4768, blast_extend.c:351-370.
    pub s_seed_off: u32,
    pub q_orig_len: u32,
    pub s_orig_len: u32,
    pub raw_score: i32,
    pub e_value: f64,
    /// Number of identical residues (for Blast_HSPTest)
    /// Reference: blast_hits.c:Blast_HSPGetNumIdentitiesAndPositives
    pub num_ident: u32,
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1381
    // ```c
    // void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
    // {
    //     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
    //         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //               ScoreCompareHSPs);
    // ```
    /// Nonprinted position in the current NCBI-equivalent `BlastHSPList`
    /// before score/link qsort steps. Used only when the NCBI comparator
    /// returns equality, so no biological tie-break absent from NCBI is added.
    pub hsp_list_order: u32,
    /// NCBI ELinkOrderingMethod: 0 = small gaps, 1 = large gaps
    /// Reference: link_hsps.c line 973: H->ordering_method = ordering_method;
    pub ordering_method: u8,
    /// NCBI: H->linked_set - TRUE if this HSP is part of a multi-HSP chain
    /// Reference: link_hsps.c lines 959-962, 972
    pub linked_set: bool,
    /// NCBI: H->start_of_chain - TRUE if this HSP is the head of a chain
    /// Reference: link_hsps.c line 955
    pub start_of_chain: bool,
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:63-65
    // ```c
    // typedef struct BlastHSPLink {
    //    struct LinkHSPStruct* link[eOrderingMethods]; /**< Best
    //                                                choice of HSP to link with */
    // ```
    /// Stable wrapper id used to preserve `hsp_link.link[ordering_method]`
    /// across the final translated-query qsort/replay phase.
    pub link_id: u32,
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1022-1028
    // ```c
    // /* If the HSP has no "link" connect the "next", otherwise follow the "link"
    //    chain down, connecting them with "next" and "prev". */
    // ordering_method = H->ordering_method;
    // if (H->hsp_link.link[ordering_method] == NULL)
    // ```
    /// Stable id of the next HSP in the selected chain, if any.
    pub chain_next_link_id: Option<u32>,
    // NCBI reference: c++/src/algo/blast/core/link_hsps.c:63-66
    // ```c
    // typedef struct BlastHSPLink {
    //    struct LinkHSPStruct* link[eOrderingMethods]; /**< Best
    //                                                choice of HSP to link with */
    //    Int2 num[eOrderingMethods]; /**< number of HSP in the ordering. */
    // ```
    // NCBI reference: c++/src/algo/blast/core/link_hsps.c:1044-1045
    // ```c
    // 			/* The first one has the number of links correct. */
    // 			num_links = H->hsp_link.num[ordering_method];
    // ```
    /// `hsp_link.num[ordering_method]` of a chain head when its chain is selected (the
    /// number of HSPs in the chain). Read only for chain heads.
    pub hsp_link_num: i16,
    // NCBI reference: c++/include/algo/blast/core/blast_hits.h:135-138
    // ```c
    //    Int4 num;             /**< How many HSP's are linked together for sum
    //                               statistics evaluation? If unset (0), this HSP is
    //                               not part of a linked set, i.e. value 0 is treated
    //                               the same way as 1. */
    // ```
    /// NCBI `BlastHSP::num`. The pairwise report prints it as `Expect(n)` and in the `N`
    /// column (the Seq-align score "sum_n").
    pub num: i32,
}

#[cfg(test)]
mod tests {
    use super::UngappedHit;

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:98-100
    // ```c
    //    Int4 offset; /**< Start of hsp */
    //    Int4 end;    /**< End of hsp */
    //    Int4 gapped_start;/**< Where the gapped extension started. */
    // ```
    // No NCBI counterpart: the layout of LOSAT's record. With the `Int4` values of NCBI in `u32`
    // fields (and the subject frame index in `u8`), one ungapped HSP of the linking list takes
    // 88 bytes on 64-bit targets instead of 152. On 32-bit targets (wasm32) `usize` was 4 bytes
    // already, so nothing is pinned there.
    #[cfg(target_pointer_width = "64")]
    #[test]
    fn ungapped_hit_is_88_bytes_on_64_bit_targets() {
        assert_eq!(std::mem::size_of::<UngappedHit>(), 88);
    }
}
