//! HSP extraction and coordinate adjustment for TBLASTX
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392, 4719-4775
//!
//! This module contains InitHSP (initial HSP with absolute coordinates) and
//! functions for converting to UngappedHit (context-relative coordinates).

use super::chaining::UngappedHit;
use super::extension::convert_coords;
use super::lookup::QueryContext;
use super::tracing::{trace_hsp_target, trace_match_target, trace_ungapped_hit_if_match};
use super::translation::QueryFrame;
use super::ungapped_hit_sort::sort_ungapped_hits_by_index_replay;
use std::cmp::Ordering;

/// NCBI BlastInitHSP equivalent - stores initial HSP with absolute coordinates
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:BlastInitHSP, c++/src/algo/blast/core/blast_extend.c:360-375 BlastSaveInitHsp
///
/// This structure stores HSPs after extension but before coordinate conversion.
/// Coordinates are stored as absolute positions in the concatenated buffer.
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:150-155
// The added PartialEq, Eq and Debug derives are used by the shadow modes and unit tests to
// compare two lists of these initial HSPs field by field. They do not change any value NCBI computes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct InitHSP {
    /// Query absolute coordinate (concatenated buffer, 0-based)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:311-315, c++/src/algo/blast/core/blast_util.c:112-116.
    pub q_start_absolute: i32,
    /// Query end absolute coordinate (concatenated buffer, 0-based)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:311-315, c++/src/algo/blast/core/blast_util.c:112-116.
    pub q_end_absolute: i32,
    /// Subject coordinate (frame-relative, 0-based)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:110-113.
    pub s_start: i32,
    /// Subject end coordinate (frame-relative, 0-based)
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:110-113.
    pub s_end: i32,
    /// Seed query offset stored in `offsets.qs_offsets.q_off`.
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:351-370, c++/src/algo/blast/core/blast_gapalign.c:4756-4768.
    pub q_seed_absolute: i32,
    /// Seed subject offset stored in `offsets.qs_offsets.s_off`.
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:351-370, c++/src/algo/blast/core/blast_gapalign.c:4756-4768.
    pub s_seed: i32,
    /// Raw score from extension
    pub score: i32,
    /// Query context index
    pub ctx_idx: usize,
    /// Subject frame index
    pub s_f_idx: usize,
    /// Query record index
    pub q_idx: u32,
    /// Subject record index
    pub s_idx: u32,
    /// Query frame
    pub q_frame: i8,
    /// Subject frame
    pub s_frame: i8,
    /// Query original length
    pub q_orig_len: usize,
    /// Subject original length
    pub s_orig_len: usize,
}

/// Trace an InitHSP if it matches the trace target.
///
/// This function converts absolute coordinates to nucleotide coordinates
/// for comparison with the trace target.
#[inline]
pub fn trace_init_hsp_if_match(stage: &str, init: &InitHSP, contexts: &[QueryContext]) {
    let Some(target) = trace_hsp_target() else {
        return;
    };

    let ctx = &contexts[init.ctx_idx];
    // Convert absolute query coords to context-relative (0-based).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768, c++/src/algo/blast/core/blast_query_info.c:311-315.
    let q_start_rel = adjust_initial_hsp_offsets(init.q_start_absolute, ctx.frame_base);
    let q_end_rel = adjust_initial_hsp_offsets(init.q_end_absolute, ctx.frame_base);
    if q_start_rel < 0 || q_end_rel < 0 || init.s_start < 0 || init.s_end < 0 {
        return;
    }

    // Convert to logical AA coords (0-based, half-open).
    let q_aa_start = q_start_rel as usize;
    let q_aa_end = q_end_rel as usize;
    let s_aa_start = init.s_start as usize;
    let s_aa_end = init.s_end as usize;

    let (q_start, q_end) = convert_coords(q_aa_start, q_aa_end, init.q_frame, init.q_orig_len);
    let (s_start, s_end) = convert_coords(s_aa_start, s_aa_end, init.s_frame, init.s_orig_len);

    if trace_match_target(target, q_start, q_end, s_start, s_end) {
        eprintln!(
            "[TRACE_HSP] stage={} raw_score={} ctx_idx={} s_f_idx={} q_frame={} s_frame={} q={}-{} s={}-{}",
            stage,
            init.score,
            init.ctx_idx,
            init.s_f_idx,
            init.q_frame,
            init.s_frame,
            q_start,
            q_end,
            s_start,
            s_end
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-304
// `q_start` is still the absolute offset in the concatenated query buffer here
// (s_AdjustInitialHSPOffsets runs later, in BLAST_GetUngappedHSPList), so this
// comparator orders HSPs that have identical context-relative coordinates in
// different query frames by context. ScoreCompareHSPs (blast_hits.c:1330-1353)
// cannot: it only sees context-relative offsets and ties.
pub(crate) fn score_compare_init_hsps_ncbi(a: &InitHSP, b: &InitHSP) -> Ordering {
    b.score
        .cmp(&a.score)
        .then_with(|| a.s_start.cmp(&b.s_start))
        .then_with(|| (b.s_end - b.s_start).cmp(&(a.s_end - a.s_start)))
        .then_with(|| a.q_start_absolute.cmp(&b.q_start_absolute))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234-235
// (end of BlastAaWordFinder, the tblastx word finder, once per subject chunk.)
// glibc's qsort is a stable merge sort and `slice::sort_by` is stable, so
// comparator-equal records keep their scan order in both.
pub(crate) fn sort_init_hsps_by_score_ncbi(init_hsps: &mut [InitHSP]) {
    init_hsps.sort_by(score_compare_init_hsps_ncbi);
}

/// NCBI s_AdjustInitialHSPOffsets equivalent
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392
///
/// NCBI code:
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2383-2391
///
/// Converts absolute coordinates to context-relative coordinates.
/// This is called in BLAST_GetUngappedHSPList (blast_gapalign.c:4756-4758).
#[inline]
pub fn adjust_initial_hsp_offsets(
    hsp_q_absolute: i32, // Absolute coordinate in concatenated buffer
    frame_base: i32,     // Context start position (absolute)
) -> i32 {
    // Subtract the context's query start offset from `q_start` to make it context-relative.
    hsp_q_absolute - frame_base
}

/// NCBI BLAST_GetUngappedHSPList equivalent
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4719-4775
///
/// NCBI code flow:
/// 1. Get context for each InitHSP (s_GetUngappedHSPContext)
/// 2. Adjust coordinates (s_AdjustInitialHSPOffsets)
/// 3. Initialize HSP (Blast_HSPInit)
/// 4. Add to hsp_list
///
/// This function converts InitHSPs (with absolute coordinates) to UngappedHits
/// (with context-relative coordinates). Reevaluation is performed separately
/// by Blast_HSPListReevaluateUngapped equivalent.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768
pub fn get_ungapped_hsp_list(
    init_hsps: Vec<InitHSP>,
    contexts: &[QueryContext],
    s_frames: &[QueryFrame],
) -> Vec<UngappedHit> {
    get_ungapped_hsp_list_into(init_hsps, contexts, s_frames, Vec::new())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4730-4732
// `get_ungapped_hsp_list` with the list passed in by the caller, as NCBI's caller may pass its own
// list; `get_ungapped_hsp_list` passes `Vec::new()`. EXPERIMENT (LOSAT_X_SEEDCOMPACT,
// x_pair_par.rs) passes a Vec reserved at `init_hsps.len()` (one hit per init HSP), so the list is
// not grown by doubling. The body after the assertion is the body `get_ungapped_hsp_list` had: it
// reads the list's length (`hsp_list_order`) and elements, never its capacity, so the capacity
// changes no element. The list must be empty, since `hsp_list_order` counts from its length and
// the final sort covers all of it.
/// `get_ungapped_hsp_list` into `ungapped_hits`, which must be empty (its capacity is kept).
pub fn get_ungapped_hsp_list_into(
    init_hsps: Vec<InitHSP>,
    contexts: &[QueryContext],
    s_frames: &[QueryFrame],
    ungapped_hits: Vec<UngappedHit>,
) -> Vec<UngappedHit> {
    let mut ungapped_hits =
        convert_init_hsps_unsorted::<false>(init_hsps, contexts, s_frames, ungapped_hits);
    sort_converted_hsps_by_score_ncbi(&mut ungapped_hits);
    ungapped_hits
}

// The loop of `get_ungapped_hsp_list_into` before its final sort, a function of its own so that
// EXPERIMENT (LOSAT_X_TBXKMERGE, x_kmerge.rs) can convert with the same statements and replace only
// the sort; `get_ungapped_hsp_list_into` is this function followed by
// `sort_converted_hsps_by_score_ncbi`. `INT4_ENDS` is true for the one-hit word finder
// (`get_ungapped_hsp_list_int4`), whose HSPs may have a length of 0 or less: their ends are stored
// as NCBI's `Int4` (`UngappedHit::set_int4_ends`) and the checks of a positive span are left out;
// with `INT4_ENDS` false (the two-hit word finder) every statement is the conversion's.
/// The HSPs of `init_hsps` in their order, pushed onto `ungapped_hits`, which must be empty.
pub(crate) fn convert_init_hsps_unsorted<const INT4_ENDS: bool>(
    init_hsps: Vec<InitHSP>,
    contexts: &[QueryContext],
    s_frames: &[QueryFrame],
    mut ungapped_hits: Vec<UngappedHit>,
) -> Vec<UngappedHit> {
    assert!(
        ungapped_hits.is_empty(),
        "get_ungapped_hsp_list_into: the list passed in must be empty"
    );

    for init_hsp in init_hsps {
        let ctx = &contexts[init_hsp.ctx_idx];
        let _s_frame = &s_frames[init_hsp.s_f_idx];

        // NCBI: s_GetUngappedHSPContext equivalent
        // Context is already stored in init_hsp.ctx_idx

        // NCBI: s_AdjustInitialHSPOffsets (c++/src/algo/blast/core/blast_gapalign.c:2384-2392)
        // Convert absolute coordinates to context-relative coordinates
        // by subtracting the `query_offset` of the context in `query_info`.
        let q_start_relative =
            adjust_initial_hsp_offsets(init_hsp.q_start_absolute, ctx.frame_base);
        let q_end_relative = adjust_initial_hsp_offsets(init_hsp.q_end_absolute, ctx.frame_base);
        let q_seed_relative = adjust_initial_hsp_offsets(init_hsp.q_seed_absolute, ctx.frame_base);

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392
        assert!(
            q_start_relative >= 0,
            "NCBI BLAST requires context-relative ungapped q_start >= 0"
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4760-4767
        if !INT4_ENDS {
            assert!(
                q_end_relative > q_start_relative,
                "NCBI BLAST requires positive ungapped query span"
            );
        }
        assert!(
            q_seed_relative >= 0,
            "NCBI BLAST requires non-negative ungapped query seed offset"
        );
        assert!(
            init_hsp.s_start >= 0,
            "NCBI BLAST requires non-negative ungapped subject start"
        );
        if !INT4_ENDS {
            assert!(
                init_hsp.s_end > init_hsp.s_start,
                "NCBI BLAST requires positive ungapped subject span"
            );
        }
        assert!(
            init_hsp.s_seed >= 0,
            "NCBI BLAST requires non-negative ungapped subject seed offset"
        );

        // NCBI: Blast_HSPInit (c++/src/algo/blast/core/blast_hits.c:150-189)
        // Sets `query.offset` to the 0-based context-relative `q_start` and `query.end` to `q_start` plus `length`.
        //
        // In LOSAT, q_start_relative is context-relative coordinate (0-based).
        // No NCBI counterpart: the record holds these `i32` (NCBI `Int4`) offsets in `u32`
        // (`UngappedHit`'s field widths); a non-negative `i32` always fits.
        let qs = u32::try_from(q_start_relative)
            .expect("NCBI BLAST requires non-negative ungapped query start");
        // `INT4_ENDS`: the ends are written below (`set_int4_ends`).
        let qe = if INT4_ENDS {
            0
        } else {
            u32::try_from(q_end_relative)
                .expect("NCBI BLAST requires non-negative ungapped query end")
        };
        let ss = u32::try_from(init_hsp.s_start)
            .expect("NCBI BLAST requires non-negative ungapped subject start");
        let se = if INT4_ENDS {
            0
        } else {
            u32::try_from(init_hsp.s_end)
                .expect("NCBI BLAST requires non-negative ungapped subject end")
        };
        let q_seed = u32::try_from(q_seed_relative)
            .expect("NCBI BLAST requires non-negative ungapped query seed");
        let s_seed = u32::try_from(init_hsp.s_seed)
            .expect("NCBI BLAST requires non-negative ungapped subject seed");

        // Offsets are already 0-based in query/subject->sequence buffers.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768, c++/src/algo/blast/core/blast_aascan.c:110-113.
        let (qs_l, qe_l) = (qs, qe);
        let (ss_l, se_l) = (ss, se);

        // Create UngappedHit with context-relative coordinates (before reevaluation)
        // NCBI: Blast_HSPInit creates HSP with original extension score
        let hsp_list_order = u32::try_from(ungapped_hits.len())
            .expect("hsp_list_order: an HSP list position fits u32 (NCBI's hspcnt is Int4)");
        let mut uh = UngappedHit {
            q_idx: init_hsp.q_idx,
            s_idx: init_hsp.s_idx,
            ctx_idx: u32::try_from(init_hsp.ctx_idx).expect(
                "ctx_idx: a query context index fits u32 (NCBI's BlastHSP context is Int4)",
            ),
            s_f_idx: u8::try_from(init_hsp.s_f_idx)
                .expect("s_f_idx: a subject frame index (0..5) fits u8"),
            q_frame: init_hsp.q_frame,
            s_frame: init_hsp.s_frame,
            q_aa_start: qs_l,
            q_aa_end: qe_l,
            s_aa_start: ss_l,
            s_aa_end: se_l,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4756-4768
            q_seed_off: q_seed,
            s_seed_off: s_seed,
            q_orig_len: u32::try_from(init_hsp.q_orig_len)
                .expect("q_orig_len: a query length fits u32 (NCBI's sequence length is Int4)"),
            s_orig_len: u32::try_from(init_hsp.s_orig_len)
                .expect("s_orig_len: a subject length fits u32 (NCBI's sequence length is Int4)"),
            raw_score: init_hsp.score, // Original extension score (before reevaluation)
            e_value: f64::INFINITY,
            num_ident: 0, // Will be computed during reevaluation
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4719-4724
            hsp_list_order,
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:476-480
            link_id: 0,
            chain_next_link_id: None,
            hsp_link_num: 0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:141-145,162
            num: 0,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:170-173 (in
        // `Blast_HSPInit`)
        // The one-hit word finder's ends, `start + length` as NCBI stores them (at or before the
        // start for an HSP of length 0 or less).
        if INT4_ENDS {
            uh.set_int4_ends(q_end_relative, init_hsp.s_end);
        }
        trace_ungapped_hit_if_match("after_get_ungapped_hsp_list", &uh);
        ungapped_hits.push(uh);
    }
    ungapped_hits
}

/// The final sort of `get_ungapped_hsp_list_into` (`Blast_HSPListSortByScore`).
fn sort_converted_hsps_by_score_ncbi(ungapped_hits: &mut [UngappedHit]) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4719-4775
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
    if !ungapped_hits_is_sorted_by_score_ncbi(ungapped_hits) {
        sort_ungapped_hits_by_score_ncbi(ungapped_hits);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
pub(crate) fn score_compare_ungapped_hits_ncbi(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.s_aa_start.cmp(&b.s_aa_start))
        .then_with(|| b.s_aa_end.cmp(&a.s_aa_end))
        .then_with(|| a.q_aa_start.cmp(&b.q_aa_start))
        .then_with(|| b.q_aa_end.cmp(&a.q_aa_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1355
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1366-1377,1380-1381
pub(crate) fn sort_ungapped_hits_by_score_ncbi(hits: &mut [UngappedHit]) {
    sort_ungapped_hits_by_score_with(hits, score_compare_ungapped_hits_ncbi);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1383
// The sort of `sort_ungapped_hits_by_score_ncbi` (the `u32` ends of the two-hit word finder's
// HSPs) and `sort_ungapped_hits_by_score_int4` (the one-hit word finder's `Int4` ends), with the
// comparator given: the is-sorted check, then the stable index-replay sort (glibc's qsort is a
// stable merge sort, `x_link_sort.rs`).
/// `Blast_HSPListSortByScore` with `compare` as `ScoreCompareHSPs`.
// `inline(always)`: each caller compiles to this body with its comparator in place, as each of
// the two bodies it replaces did.
#[inline(always)]
fn sort_ungapped_hits_by_score_with<C>(hits: &mut [UngappedHit], compare: C)
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering + Sync + Copy,
{
    if hits.len() <= 1 {
        return;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1381
    if ungapped_hits_is_sorted_by_score_with(hits, compare) {
        return;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
    // `ScoreCompareHSPs` is a partial comparator for TBLASTX HSPs. Do not add
    // frame or per-frame insertion-order keys when it returns equality; NCBI
    // delegates comparator-equal rows to the platform qsort implementation.
    sort_ungapped_hits_by_index_replay(hits, compare);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1347-1353
// `ScoreCompareHSPs` with the ends read as NCBI's `Int4` (`UngappedHit::q_end_int4`): the
// comparator of the one-hit word finder's lists, whose HSPs of length 0 or less may end below 0
// (`chaining.rs`, the impl of `UngappedHit`), and of the chunk merge of every list. On ends of 0
// or more it is `score_compare_ungapped_hits_ncbi`.
pub(crate) fn score_compare_ungapped_hits_int4(a: &UngappedHit, b: &UngappedHit) -> Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.s_aa_start.cmp(&b.s_aa_start))
        .then_with(|| b.s_end_int4().cmp(&a.s_end_int4()))
        .then_with(|| a.q_aa_start.cmp(&b.q_aa_start))
        .then_with(|| b.q_end_int4().cmp(&a.q_end_int4()))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1383
// `sort_ungapped_hits_by_score_ncbi` with `score_compare_ungapped_hits_int4`.
/// Sorts a list by score, its ends read as `Int4`.
pub(crate) fn sort_ungapped_hits_by_score_int4(hits: &mut [UngappedHit]) {
    sort_ungapped_hits_by_score_with(hits, score_compare_ungapped_hits_int4);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4760-4772
// `get_ungapped_hsp_list` for the one-hit word finder (`-window_size 0`), whose initial HSPs may
// have a length of 0 or less (aa_ungapped.c:1054,1083): `Blast_HSPInit` stores their ends,
// `start + length`, as they are (`convert_init_hsps_unsorted` with `INT4_ENDS`), and the list is
// sorted with the `Int4` comparator. On a list without such an HSP the result is
// `get_ungapped_hsp_list`'s.
/// The one-hit word finder's `get_ungapped_hsp_list`, which keeps NCBI's `Int4` ends.
pub(crate) fn get_ungapped_hsp_list_int4(
    init_hsps: Vec<InitHSP>,
    contexts: &[QueryContext],
    s_frames: &[QueryFrame],
) -> Vec<UngappedHit> {
    let mut hits = convert_init_hsps_unsorted::<true>(init_hsps, contexts, s_frames, Vec::new());
    sort_ungapped_hits_by_score_int4(&mut hits);
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1355-1369
pub(crate) fn ungapped_hits_is_sorted_by_score_ncbi(hits: &[UngappedHit]) -> bool {
    ungapped_hits_is_sorted_by_score_with(hits, score_compare_ungapped_hits_ncbi)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1355-1369
/// `Blast_HSPListIsSortedByScore` with `compare` as `ScoreCompareHSPs`.
fn ungapped_hits_is_sorted_by_score_with<C>(hits: &[UngappedHit], compare: C) -> bool
where
    C: Fn(&UngappedHit, &UngappedHit) -> Ordering,
{
    if hits.len() <= 1 {
        return true;
    }

    for index in 0..hits.len() - 1 {
        if compare(&hits[index], &hits[index + 1]) == Ordering::Greater {
            return false;
        }
    }
    true
}

#[cfg(test)]
#[path = "blast_gapalign_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "init_hsp_size_tests.rs"]
mod init_hsp_size_tests;
