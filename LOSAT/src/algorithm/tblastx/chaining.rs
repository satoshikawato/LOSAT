//! Ungapped HSP record of TBLASTX
//!
//! `UngappedHit` is the lightweight HSP the scan, the linking and the report steps pass around.

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-101,126-128,133,158
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:246
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
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768, c++/src/algo/blast/core/blast_extend.c:351-370.
    pub q_seed_off: u32,
    /// NCBI `BlastHSP::subject.gapped_start` copied from
    /// `init_hsp->offsets.qs_offsets.s_off`.
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768, c++/src/algo/blast/core/blast_extend.c:351-370.
    pub s_seed_off: u32,
    pub q_orig_len: u32,
    pub s_orig_len: u32,
    pub raw_score: i32,
    pub e_value: f64,
    /// Number of identical residues (for Blast_HSPTest)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:Blast_HSPGetNumIdentitiesAndPositives
    pub num_ident: u32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
    /// Nonprinted position in the current NCBI-equivalent `BlastHSPList`
    /// before score/link qsort steps. Used only when the NCBI comparator
    /// returns equality, so no biological tie-break absent from NCBI is added.
    pub hsp_list_order: u32,
    /// NCBI ELinkOrderingMethod: 0 = small gaps, 1 = large gaps
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c
    pub ordering_method: u8,
    /// Whether this HSP belongs to a multi-HSP chain (NCBI `H->linked_set`).
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c lines 959-962, 972
    pub linked_set: bool,
    /// Whether this HSP is the head of a chain (NCBI `H->start_of_chain`).
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c line 955
    pub start_of_chain: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:63-65
    /// Stable wrapper id used to preserve `hsp_link.link[ordering_method]`
    /// across the final translated-query qsort/replay phase.
    pub link_id: u32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1022-1028
    /// Stable id of the next HSP in the selected chain, if any.
    pub chain_next_link_id: Option<u32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:63-66
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1044-1045
    /// `hsp_link.num[ordering_method]` of a chain head when its chain is selected (the
    /// number of HSPs in the chain). Read only for chain heads.
    pub hsp_link_num: i16,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:135-138
    /// NCBI `BlastHSP::num`. The pairwise report prints it in parentheses after `Expect` and in
    /// the `N` column (the Seq-align score "sum_n").
    pub num: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:170-173 (in `Blast_HSPInit`,
// called by `BLAST_GetUngappedHSPList`, blast_gapalign.c:4760-4767)
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1054,1083 (in
// `s_BlastAaExtendOneHit`)
// No NCBI counterpart: where LOSAT keeps the `Int4` ends of an HSP whose length is 0 or negative.
// The one-hit word finder of TBLASTX (`-window_size 0`) can save such an HSP (a word whose
// running sum is reset after its best part), and NCBI's `end`, `offset` plus `length`, is then at or
// before the offset, and below 0 near a frame's start. The fields `q_aa_end` and `s_aa_end` hold
// that `Int4` in two's complement; a positive length leaves every end non-negative, where the
// `Int4` and the `u32` are the same number. These methods are the one reader and writer of the
// `Int4` ends: the one-hit path (`blast_gapalign::get_ungapped_hsp_list_int4`, the chunk offsets
// and merge of `blast_engine/run_impl.rs`, the first linking call of `subject_link.rs` and the
// sort after it) uses them from the conversion of such an HSP until the HSPs of length 0 or less
// are dropped after the re-evaluation's fetch of the subject, which NCBI makes for every list that
// still holds an HSP (`subject_link::drop_unread_int4_inverted`), where NCBI's re-evaluation
// deletes them (a length of 0 or less scores 0, below every cutoff: blast_hits.c:687-699,
// 449-453). Nothing after that sees such an HSP; the reap by the preliminary E-value before it
// reads only the E-value of one (the first linking call's, without sum statistics its own).
// In a chunk's list, from the conversion (`blast_gapalign::convert_init_hsps_unsorted`) to the move
// by the chunk's offset (`blast_engine/run_impl.rs`, `adjust_tblastx_chunk_subject_offsets`), the
// fields `s_aa_start`, `s_aa_end` and `s_seed_off` hold NCBI's chunk-local `Int4` the same way: below
// 0 for an HSP or a seed before its chunk, which only the scan of a range of masks stored out of order
// gives (`common::translated_chunks::scan_lead`); the move leaves them at 0 or more.
impl UngappedHit {
    /// NCBI's `subject.offset` (`Int4`).
    #[inline]
    pub(crate) fn s_start_int4(&self) -> i32 {
        self.s_aa_start as i32
    }

    /// NCBI's `subject.gapped_start` (`Int4`).
    #[inline]
    pub(crate) fn s_seed_int4(&self) -> i32 {
        self.s_seed_off as i32
    }

    /// NCBI's `query.end` (`Int4`).
    #[inline]
    pub(crate) fn q_end_int4(&self) -> i32 {
        self.q_aa_end as i32
    }

    /// NCBI's `subject.end` (`Int4`).
    #[inline]
    pub(crate) fn s_end_int4(&self) -> i32 {
        self.s_aa_end as i32
    }

    /// Sets NCBI's `query.end` and `subject.end` (`Int4`).
    #[inline]
    pub(crate) fn set_int4_ends(&mut self, q_end: i32, s_end: i32) {
        self.q_aa_end = q_end as u32;
        self.s_aa_end = s_end as u32;
    }

    /// Whether the HSP's length (`query.end - query.offset`, the subject's being the same) is 0
    /// or less.
    #[inline]
    pub(crate) fn is_int4_inverted(&self) -> bool {
        self.q_end_int4() <= self.q_aa_start as i32
    }
}

#[cfg(test)]
#[path = "chaining_tests.rs"]
mod tests;
