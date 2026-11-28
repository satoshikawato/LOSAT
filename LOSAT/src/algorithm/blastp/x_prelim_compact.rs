//! blastp's preliminary hit lists hold a 64-byte element instead of a `BlastpHsp` (168 bytes)
//! until the redo reads them.
//!
//! The preliminary stage keeps, for every query of a search, the HSP lists of up to the
//! preliminary hit list size of subjects. Each preliminary HSP is a `BlastpHsp`, most of whose
//! fields (identities, bit score, edit script, ...) the preliminary stage never sets and the redo
//! never reads: the redo turns every HSP of a list back into the `BlastpPreliminaryHsp` it was
//! made from (`preliminary_hit_from_local_match_hsp`). With many matches in a block of queries
//! these lists are the largest part of the peak at the end of the subject scan.
//!
//! The subject scan keeps an `XPrelimHsp` (the `BlastpPreliminaryHsp` and its E-value) in the
//! `BlastpHspList` / `BlastpHitList` helpers, made generic over the element (`HitListHsp`): the
//! same E-values, scores and offsets in the same comparisons, and every sort of HSPs is stable
//! (`sort_by`), so the lists end in the same order. Before the redo a list is made into the
//! `BlastpHitList` of full HSPs the redo reads (`x_expand_hit_list`, with the HSPs built by
//! `blastp_hsp_from_preliminary`): the one-thread redo does it query by query and drops the
//! expanded list after the query's redo (nothing reads it afterwards, `x_prelim_free.rs`); a
//! search with more than one thread expands every list before the redo (its parallel paths read
//! the lists of all queries), but where the redo runs in windows (LOSAT_X_BPREDOPAR),
//! the windows expand the lists as they reach them (`x_redo_par.rs`). The
//! lists of query chunks (`-query` longer than the split size) are expanded where a chunk's search
//! ends, so the chunk merge is unchanged.
//!
//! With `LOSAT_TIMING=1` a search prints `[TIMING] x_bppcompact: expanded_lists=<n>` (the hit lists
//! expanded during the search, those of its query chunks included).

use std::sync::atomic::{AtomicU64, Ordering};

use super::gapalign::BlastpPreliminaryHsp;
use super::hsp::{BlastpHitList, BlastpHspList, HitListHsp};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:151-178
// NCBI keeps one BlastHSP; this is the part of it the preliminary stage sets, with the E-value of
// Blast_HSPListGetEvalues (blast_engine.c:1523-1530).
/// A preliminary HSP as the subject scan keeps it.
#[derive(Debug, Clone, Copy)]
pub(crate) struct XPrelimHsp {
    pub hsp: BlastpPreliminaryHsp,
    pub e_value: f64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:780-794
/// 0-based position of a preliminary PHI-BLAST HSP in its subject's `PhiSearch::prelim` list.
///
/// Carried on `XPrelimHsp` without growing `XPrelimHsp` (64 B) by storing `index + 1` as a
/// positive `i32` in `BlastpPreliminaryHsp::query_frame` (`0` remains the unset sentinel used by
/// non-PHI protein HSPs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PhiPrelimOrdinal(usize);

impl PhiPrelimOrdinal {
    /// Wraps a 0-based index into a subject's preliminary PHI-BLAST HSP list (`Int4` range).
    #[inline]
    pub(crate) fn new(index: usize) -> Self {
        assert!(
            index < i32::MAX as usize,
            "PHI preliminary HSP index {index} exceeds Int4 range"
        );
        Self(index)
    }

    /// Returns the 0-based index into the subject's preliminary PHI-BLAST HSP list.
    #[inline]
    pub(crate) fn index(self) -> usize {
        self.0
    }

    /// The ordinal that [`XPrelimHsp::with_phi_ordinal`] attached to the preliminary HSP `hsp`, or
    /// `None` when unset (`query_frame <= 0`).
    #[inline]
    pub(crate) fn of(hsp: &BlastpPreliminaryHsp) -> Option<Self> {
        usize::try_from(hsp.query_frame)
            .ok()
            .and_then(|value| value.checked_sub(1))
            .map(Self)
    }
}

impl XPrelimHsp {
    /// Attaches the HSP's position in its subject's `PhiSearch::prelim` list.
    ///
    /// Stored in `self.hsp.query_frame` as `ordinal.index() + 1` (positive `i32`) so
    /// `XPrelimHsp` stays 64 bytes. Safe for PHI-BLAST because:
    /// - `HitListHsp`, `BlastpHitList`, `PrelimCulling`, `PrelimBestHit`, `exclude_fetched`,
    ///   `listed_ordinals`, `rename_to_ranks`, and `x_db_keep` never read `query_frame`.
    /// - Preliminary endpoint purges (`preliminary_hits_share_*_endpoint`) and
    ///   `query_split::merge_two_hsps` do not read `query_frame` (and never run for PHI-BLAST).
    /// - `ResolvedBlastpArgs::phi` forces `NoCompositionBasedStats` and disables query splitting,
    ///   so PHI-BLAST lists never pass through `x_expand_hit_list`, `x_finish_compact_lists`,
    ///   `kappa`, `x_redo_par`, or ordinary `traceback`; `phi_traceback` consumes
    ///   `compact_hit_lists` directly and builds final `BlastpHsp`s with `query_frame = 0`.
    #[inline]
    pub(crate) fn with_phi_ordinal(mut self, ordinal: PhiPrelimOrdinal) -> Self {
        self.hsp.query_frame = (ordinal.0 + 1) as i32;
        self
    }

    /// Returns the PHI preliminary ordinal previously attached by [`Self::with_phi_ordinal`],
    /// or `None` when unset (`query_frame <= 0`).
    #[inline]
    pub(crate) fn phi_ordinal(&self) -> Option<PhiPrelimOrdinal> {
        PhiPrelimOrdinal::of(&self.hsp)
    }
}

// No NCBI counterpart: the values `BlastpHsp` holds for the same HSP
// (`blastp_hsp_from_preliminary`: starts + 1).
impl HitListHsp for XPrelimHsp {
    #[inline]
    fn e_value(&self) -> f64 {
        self.e_value
    }
    #[inline]
    fn raw_score(&self) -> i32 {
        self.hsp.raw_score
    }
    #[inline]
    fn offsets(&self) -> (usize, usize, usize, usize) {
        (
            (self.hsp.query_start + 1) as usize,
            self.hsp.query_end as usize,
            (self.hsp.subject_start + 1) as usize,
            self.hsp.subject_end as usize,
        )
    }
    #[inline]
    fn query_length(&self) -> usize {
        self.hsp.query_length
    }
    #[inline]
    fn move_to_query(&mut self, q_idx: u32, offset: i32, query_length: usize) {
        self.hsp.query_context = q_idx as i32;
        self.hsp.q_idx = q_idx;
        self.hsp.query_length = query_length;
        self.hsp.query_start += offset;
        self.hsp.query_end += offset;
        self.hsp.gapped_query_start += offset;
    }
}

// No NCBI counterpart: rebuilds the hit list of full HSPs from the compact one (same counts, same
// list order, the HSPs `blastp_hsp_from_preliminary` makes).
/// The `BlastpHitList` of full HSPs the redo reads for this list.
pub(crate) fn x_expand_hit_list(list: BlastpHitList<XPrelimHsp>) -> BlastpHitList {
    X_EXPANDED_LISTS.fetch_add(1, Ordering::Relaxed);
    let BlastpHitList {
        hsplist_count,
        hsplist_max,
        worst_evalue,
        low_score,
        heapified,
        hsplist_array,
        hsplist_current,
        num_hits,
    } = list;
    let mut expanded = Vec::with_capacity(hsplist_array.capacity());
    expanded.extend(hsplist_array.into_iter().map(|hsp_list| {
        BlastpHspList {
            oid: hsp_list.oid,
            query_index: hsp_list.query_index,
            hsps: hsp_list
                .hsps
                .iter()
                .map(|hsp| super::blast_engine::blastp_hsp_from_preliminary(&hsp.hsp, hsp.e_value))
                .collect(),
            best_evalue: hsp_list.best_evalue,
        }
    }));
    BlastpHitList {
        hsplist_count,
        hsplist_max,
        worst_evalue,
        low_score,
        heapified,
        hsplist_array: expanded,
        hsplist_current,
        num_hits,
    }
}

// No NCBI counterpart: hands the preliminary hit lists on as lists of full HSPs.
/// Every compact list becomes the list of full HSPs of its query (`full` holds the lists of the
/// queries whose chunks replaced them).
pub(crate) fn x_finish_compact_lists(
    full: &mut [Option<BlastpHitList>],
    compact: &mut [Option<BlastpHitList<XPrelimHsp>>],
) {
    for (full, compact) in full.iter_mut().zip(compact.iter_mut()) {
        if let Some(list) = compact.take() {
            *full = Some(x_expand_hit_list(list));
        }
    }
}

static X_EXPANDED_LISTS: AtomicU64 = AtomicU64::new(0);

/// The hit lists expanded so far (for the timing line).
pub(crate) fn x_expanded_lists() -> u64 {
    X_EXPANDED_LISTS.load(Ordering::Relaxed)
}

#[cfg(test)]
#[path = "x_prelim_compact_tests.rs"]
mod tests;
