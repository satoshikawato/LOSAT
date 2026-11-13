//! The even-gap linking of the HSPs of one subject with sum statistics, for a gapped search with an
//! untranslated query and subject: psiblast's `-sum_stats true` (NCBI core/link_hsps.c
//! `BLAST_LinkHsps` with `longest_intron <= 0`, `s_BlastEvenGapLinkHSPs`).
//!
//! NCBI's blastp and blastn have no `-sum_stats`, and the other programs that link take other
//! forms: TBLASTX links ungapped HSPs of translated sequences (`algorithm::tblastx::
//! sum_stats_linking`, its own kernels for its `UngappedHit`), and BLASTX and TBLASTN use the
//! uneven-gap linking (`s_BlastUnevenGapLinkHSPs`). So PSIBLAST's rounds (the matrix round 1 and
//! the PSSM rounds) are the users of this port.
//!
//! The port is literal, with NCBI's shortcuts (`use_current_max`, `path_changed`, the `changed`
//! flags, `next_larger`): they decide which of equal sums a chain takes, so the result depends on
//! them. The sums are `i64` where NCBI keeps `Int4` (a sum past `Int4` is undefined behaviour in
//! NCBI; Owner decision 2026-10-10). The callers give the HSPs in the list's order and read back the
//! list NCBI builds (`LinkedHsp`, in NCBI's order of `hsp_list->hsp_array` after the linking), then
//! sort it by score (`Blast_HSPListSortByScore`) as `BLAST_LinkHsps` does.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1760-1814
//! ```c
//!    /* Remove any information on number of linked HSPs from previous
//!       linking. */
//!    for (index = 0; index < hsp_list->hspcnt; ++index)
//!        hsp_list->hsp_array[index]->num = 1;
//!
//!    /* Link up the HSP's for this hsp_list. */
//!    if (link_hsp_params->longest_intron <= 0) {
//!        s_BlastEvenGapLinkHSPs(program_number, hsp_list, query_info,
//!                              subject_length, sbp, link_hsp_params,
//!                              gapped_calculation);
//! ...
//!    /* Sort the HSP array by score */
//!    Blast_HSPListSortByScore(hsp_list);
//! ```

use std::cmp::Ordering;

use crate::stats::sum_statistics::{gap_decay_divisor, ncbi_large_gap_sum_e, small_gap_sum_e};

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:66-71
// ```c
// #define BLAST_GAP_PROB 0.5         /**< Gap probability for ungapped search */
// #define BLAST_GAP_PROB_GAPPED 1.0  /**< Gap probability for gapped search */
// #define BLAST_GAP_DECAY_RATE 0.5   /**< Gap decay rate for ungapped search */
// #define BLAST_GAP_DECAY_RATE_GAPPED 0.1/**< Gap decay rate for gapped search */
// #define BLAST_GAP_SIZE 40          /**< Default gap size */
// #define BLAST_OVERLAP_SIZE 9      /**< Default overlap size */
// ```
// With a gap probability of 1, the preliminary stage (which also links with small gaps) never
// takes a large-gap set of more than one HSP (`1 - gap_prob == 0` gives it `INT4_MAX`); the
// traceback links with large gaps only and does not divide by the probability.
/// `BLAST_GAP_PROB_GAPPED`.
pub(crate) const BLAST_GAP_PROB_GAPPED: f64 = 1.0;
/// `BLAST_GAP_DECAY_RATE_GAPPED`.
pub(crate) const BLAST_GAP_DECAY_RATE_GAPPED: f64 = 0.1;
/// `BLAST_GAP_SIZE`.
pub(crate) const BLAST_GAP_SIZE: i32 = 40;
/// `BLAST_OVERLAP_SIZE`.
pub(crate) const BLAST_OVERLAP_SIZE: i32 = 9;

/// The linking parameters (`BlastLinkHSPParameters`) of a gapped protein search with sum
/// statistics.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:609-617
/// ```c
///    if (program_number == eBlastTypeBlastn || !gapped_calculation) {
///       params->gap_prob = BLAST_GAP_PROB;
///       params->gap_decay_rate = BLAST_GAP_DECAY_RATE;
///    } else {
///       params->gap_prob = BLAST_GAP_PROB_GAPPED;
///       params->gap_decay_rate = BLAST_GAP_DECAY_RATE_GAPPED;
///    }
///    params->gap_size = BLAST_GAP_SIZE;
///    params->overlap_size = BLAST_OVERLAP_SIZE;
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:633-637
/// ```c
///    if (gapped_calculation) {
///       /* FIXME, is this correct?? */
///       hit_params->link_hsp_params->cutoff_small_gap =
///                                 word_params->cutoff_score_min;
/// ```
/// A gapped search sets no `cutoff_big_gap` (`CalculateLinkHSPCutoffs` runs only for ungapped
/// searches, blast_engine.c:1448-1455), so it stays at the 0 of the `calloc`. The preliminary stage
/// sets `cutoff_small_gap` to the smallest cutoff of the initial words
/// (`BlastLinkHSPParametersUpdate`, blast_engine.c:1399, blast_setup.c:1024); the traceback and
/// the composition-based redo link with the parameters of their own setup
/// (`SThreadLocalDataArraySetup` -> `BLAST_GapAlignSetUp` -> `BlastHitSavingParametersNew`), which
/// no update reaches, so their `cutoff_small_gap` is 0 too (`EvenGapLinkParams::traceback`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EvenGapLinkParams {
    /// `gap_prob`.
    pub(crate) gap_prob: f64,
    /// `gap_decay_rate`.
    pub(crate) gap_decay_rate: f64,
    /// `gap_size`.
    pub(crate) gap_size: i32,
    /// `overlap_size`.
    pub(crate) overlap_size: i32,
    /// `cutoff_small_gap`: 0 links with large gaps only (`ignore_small_gaps`).
    pub(crate) cutoff_small_gap: i32,
    /// `cutoff_big_gap`.
    pub(crate) cutoff_big_gap: i32,
}

impl EvenGapLinkParams {
    /// The preliminary stage's parameters: `cutoff_small_gap` is the smallest cutoff score of the
    /// initial words (`word_params->cutoff_score_min`).
    pub(crate) fn preliminary(word_cutoff_score_min: i32) -> Self {
        Self {
            cutoff_small_gap: word_cutoff_score_min,
            ..Self::traceback()
        }
    }

    /// The parameters of the traceback and of the composition-based redo (no update: both cutoffs
    /// 0, large gaps only).
    pub(crate) fn traceback() -> Self {
        Self {
            gap_prob: BLAST_GAP_PROB_GAPPED,
            gap_decay_rate: BLAST_GAP_DECAY_RATE_GAPPED,
            gap_size: BLAST_GAP_SIZE,
            overlap_size: BLAST_OVERLAP_SIZE,
            cutoff_small_gap: 0,
            cutoff_big_gap: 0,
        }
    }
}

/// One HSP as the linking reads it (`BlastHSP`): its context, its query and subject ranges
/// (0-based, the end exclusive, `query.offset`/`query.end`) and its raw score. The subject is not
/// translated, so its frame is 0 for every HSP.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LinkHsp {
    pub(crate) context: usize,
    pub(crate) query_offset: i32,
    pub(crate) query_end: i32,
    pub(crate) subject_offset: i32,
    pub(crate) subject_end: i32,
    pub(crate) score: i32,
}

/// The values of a query context the linking reads: the gapped Karlin block of the context
/// (`sbp->kbp_gap[context]`, scaled in the composition-based redo) and the context's lengths
/// (`query_info->contexts[context]`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LinkContext {
    pub(crate) lambda: f64,
    pub(crate) log_k: f64,
    pub(crate) query_length: i32,
    pub(crate) length_adjustment: i32,
    pub(crate) eff_searchsp: i64,
}

/// An HSP of the linked list: its index in the input, its e-value (the linked set's) and the size
/// of its linked set (`BlastHSP::num`, 1 for an HSP that is not in a set).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LinkedHsp {
    pub(crate) index: usize,
    pub(crate) evalue: f64,
    pub(crate) num: i32,
}

/// `eOrderingMethods`: `eLinkSmallGaps` (0) and `eLinkLargeGaps` (1).
const SMALL: usize = 0;
const LARGE: usize = 1;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:67-97 (LinkHSPStruct, BlastHSPLink)
// ```c
// typedef struct BlastHSPLink {
//    struct LinkHSPStruct* link[eOrderingMethods]; /**< Best
//                                            choice of HSP to link with */
//    Int2  num[eOrderingMethods]; /**< number of HSP in the ordering. */
//    Int4  sum[eOrderingMethods]; /**< Sum-Score of HSP. */
//    double   xsum[eOrderingMethods]; /**< Sum-Score of HSP,
//                                      multiplied by the appropriate Lambda. */
//    Int4  changed; /**< Has the link been changed since previous access? */
// } BlastHSPLink;
// ...
//    Int4 linked_to;     /**< Where this HSP is linked to? */
//    Boolean linked_set;      /**< Is this HSp part of a linked set? */
//    Boolean start_of_chain;  /**< If TRUE, this HSP starts a chain along the
//                                "link" pointer. */
//    Int4 q_offset_trim;   /**< Start of trimmed hsp in query */
// ...
//    ELinkOrderingMethod ordering_method; /**< Which method (max or
//                                              no max for gaps) was
//                                              used for this HSP? */
//    struct LinkHSPStruct* prev;         /**< Previous HSP in a set, if any */
//    struct LinkHSPStruct* next;         /**< Next HSP in a set, if any */
// ```
/// `LinkHSPStruct`; the links are slots of `LinkState::nodes` (slot 0 is `hp_start`, the HSP of
/// sorted position `i` is slot `i + 1`).
#[derive(Debug, Clone, Default)]
struct Node {
    link: [Option<usize>; 2],
    num: [i32; 2],
    sum: [i64; 2],
    xsum: [f64; 2],
    changed: bool,
    linked_to: i32,
    linked_set: bool,
    start_of_chain: bool,
    q_offset_trim: i32,
    q_end_trim: i32,
    s_offset_trim: i32,
    s_end_trim: i32,
    ordering_method: usize,
    evalue: f64,
    prev: Option<usize>,
    next: Option<usize>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:100-109
// ```c
// typedef struct LinkHelpStruct {
//   LinkHSPStruct* ptr;         /**< The HSP to which the info belongs */
//   Int4 q_off_trim;            /**< query start of trimmed HSP */
//   Int4 s_off_trim;            /**< subject start of trimmed HSP */
//   Int4 sum[eOrderingMethods]; /**< raw score of linked set containing HSP(?) */
//   Int4 maxsum1;               /**< threshold for stopping link attempts (?) */
//   Int4 next_larger;           /**< offset into array of HelpStructs containing
//                                    HSP with higher score, used in bailout
//                                    calculations */
// } LinkHelpStruct;
// ```
/// `LinkHelpStruct` (`maxsum1` is read only in NCBI's disabled `if(0)` test, so it is not kept).
#[derive(Debug, Clone, Copy, Default)]
struct Helper {
    ptr: usize,
    q_off_trim: i32,
    s_off_trim: i32,
    sum: [i64; 2],
    next_larger: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:276-318
// ```c
// s_RevCompareHSPsTbn(const void *v1, const void *v2)
// ...
//    if (h1->context < h2->context)
//       return -1;
//    else if (h1->context > h2->context)
//       return 1;
//
// 	if (SIGN(h1->subject.frame) != SIGN(h2->subject.frame))
// ...
// 	if (h1->query.offset < h2->query.offset)
// 		return  1;
// 	if (h1->query.offset > h2->query.offset)
// 		return -1;
// 	if (h1->query.end < h2->query.end)
// 		return  1;
// 	if (h1->query.end > h2->query.end)
// 		return -1;
// 	if (h1->subject.offset < h2->subject.offset)
// 		return  1;
// 	if (h1->subject.offset > h2->subject.offset)
// 		return -1;
// 	if (h1->subject.end < h2->subject.end)
// 		return  1;
// 	if (h1->subject.end > h2->subject.end)
// 		return -1;
// 	return 0;
// ```
/// `s_RevCompareHSPsTbn` (the subject frames are all 0).
fn rev_compare_tbn(h1: &LinkHsp, h2: &LinkHsp) -> Ordering {
    h1.context
        .cmp(&h2.context)
        .then(h2.query_offset.cmp(&h1.query_offset))
        .then(h2.query_end.cmp(&h1.query_end))
        .then(h2.subject_offset.cmp(&h1.subject_offset))
        .then(h2.subject_end.cmp(&h1.subject_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:197-224
// ```c
// s_RevCompareHSPs(const void *v1, const void *v2)
// ...
//    if (h1->context < h2->context)
//       return -1;
//    else if (h1->context > h2->context)
//       return 1;
//
// 	if (h1->query.offset < h2->query.offset)
// 		return  1;
// 	if (h1->query.offset > h2->query.offset)
// 		return -1;
// 	/* Necessary in case both HSP's have the same query offset. */
// 	if (h1->subject.offset < h2->subject.offset)
// 		return 1;
// 	if (h1->subject.offset > h2->subject.offset)
// 		return -1;
// 	return 0;
// ```
/// `s_RevCompareHSPs`.
fn rev_compare(h1: &LinkHsp, h2: &LinkHsp) -> Ordering {
    h1.context
        .cmp(&h2.context)
        .then(h2.query_offset.cmp(&h1.query_offset))
        .then(h2.subject_offset.cmp(&h1.subject_offset))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:119-144
// ```c
// s_FwdCompareHSPs(const void* v1, const void* v2)
// ...
//     if (h1->context < h2->context)
//         return -1;
//     if (h1->context > h2->context)
//         return 1;
//
//     if (h1->query.offset < h2->query.offset)
//         return -1;
//     if (h1->query.offset > h2->query.offset)
//         return 1;
//     /* Necessary in case both HSP's have the same query offset. */
//     if (h1->subject.offset < h2->subject.offset)
//         return -1;
//     if (h1->subject.offset > h2->subject.offset)
//         return 1;
//
//     return 0;
// ```
/// `s_FwdCompareHSPs`.
fn fwd_compare(h1: &LinkHsp, h2: &LinkHsp) -> Ordering {
    h1.context
        .cmp(&h2.context)
        .then(h1.query_offset.cmp(&h2.query_offset))
        .then(h1.subject_offset.cmp(&h2.subject_offset))
}

/// Links the HSPs of one subject (`s_BlastEvenGapLinkHSPs` for an untranslated query and
/// subject) and returns the list NCBI builds, in its order: every HSP once, with its linked set's
/// e-value and size. `contexts` is indexed by `LinkHsp::context`; `subject_length` is the
/// subject's length in letters.
///
/// NCBI's `qsort` calls are stable sorts here (the comparators are partial: HSPs with the same
/// start keep the order of the previous sort), as TBLASTX's port sorts them.
pub(crate) fn even_gap_link_hsps(
    hsps: &[LinkHsp],
    contexts: &[LinkContext],
    subject_length: i32,
    params: &EvenGapLinkParams,
) -> Vec<LinkedHsp> {
    let total_number_of_hsps = hsps.len();
    if total_number_of_hsps == 0 {
        return Vec::new();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:468-493
    // ```c
    //    /* For convenience, include overlap size into the gap size */
    // 	window_size = link_hsp_params->gap_size + link_hsp_params->overlap_size + 1;
    //    trim_size = (link_hsp_params->overlap_size + 1) / 2;
    // 	gap_prob = link_hsp_params->gap_prob;
    // 	gap_decay_rate = link_hsp_params->gap_decay_rate;
    // ...
    //    /* Sort by (reverse) position. */
    //    if (kTranslatedQuery) {
    // ...
    //    } else {
    //       qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //             s_RevCompareHSPsTbn);
    //    }
    //
    //    cutoff[0] = link_hsp_params->cutoff_small_gap;
    //    cutoff[1] = link_hsp_params->cutoff_big_gap;
    //
    //    ignore_small_gaps = (cutoff[0] == 0);
    // ```
    let window_size = params.gap_size + params.overlap_size + 1;
    let trim_size = (params.overlap_size + 1) / 2;
    let gap_prob = params.gap_prob;
    let gap_decay_rate = params.gap_decay_rate;
    let cutoff = [params.cutoff_small_gap as i64, params.cutoff_big_gap as i64];
    let ignore_small_gaps = params.cutoff_small_gap == 0;

    // `link_hsp_array`: the input positions in the order of `s_RevCompareHSPsTbn`.
    let mut order: Vec<usize> = (0..total_number_of_hsps).collect();
    order.sort_by(|&a, &b| rev_compare_tbn(&hsps[a], &hsps[b]));
    // Slot 0 is `hp_start`; the HSP at sorted position `i` is slot `i + 1`.
    let mut nodes: Vec<Node> = vec![Node::default(); total_number_of_hsps + 1];
    let hsp_of = |slot: usize| -> &LinkHsp { &hsps[order[slot - 1]] };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:510-536
    // ```c
    //    /* Put entries from different strands into separate 'query_frame's. */
    //    {
    //       Int4 cur_frame=0;
    //       Int4 strand_factor = (kTranslatedQuery ? 3 : 1);
    //       for (index=0;index<number_of_hsps;index++)
    //       {
    //         H=link_hsp_array[index];
    //         H->start_of_chain = FALSE;
    //         hp_frame_number[cur_frame]++;
    //
    //         H->prev= index ? link_hsp_array[index-1] : NULL;
    //         H->next= index<(number_of_hsps-1) ? link_hsp_array[index+1] : NULL;
    //         if (H->prev != NULL &&
    //             ((H->hsp->context/strand_factor) !=
    // 	     (H->prev->hsp->context/strand_factor) ||
    //              (SIGN(H->hsp->subject.frame) != SIGN(H->prev->hsp->subject.frame))))
    //         { /* If frame switches, then start new list. */
    //            hp_frame_number[cur_frame]--;
    //            hp_frame_number[++cur_frame]++;
    //            hp_frame_start[cur_frame] = H;
    //            H->prev->next = NULL;
    //            H->prev = NULL;
    //         }
    //       }
    //       num_query_frames = cur_frame+1;
    //    }
    // ```
    // (frame start slot, number of HSPs) of each query frame (a context here).
    let mut frames: Vec<(usize, usize)> = vec![(1, 0)];
    for index in 0..total_number_of_hsps {
        let slot = index + 1;
        nodes[slot].start_of_chain = false;
        frames.last_mut().expect("one frame").1 += 1;
        nodes[slot].prev = (index > 0).then_some(slot - 1);
        nodes[slot].next = (index + 1 < total_number_of_hsps).then_some(slot + 1);
        if index > 0 && hsp_of(slot).context != hsp_of(slot - 1).context {
            frames.last_mut().expect("one frame").1 -= 1;
            frames.push((slot, 1));
            nodes[slot - 1].next = None;
            nodes[slot].prev = None;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:538-552
    // ```c
    //    for (index=0;index<number_of_hsps;index++)
    //    {
    //        Int4 q_length, s_length;
    //        H = link_hsp_array[index];
    //        hsp = H->hsp;
    //        q_length = (hsp->query.end - hsp->query.offset) / 4;
    //        s_length = (hsp->subject.end - hsp->subject.offset) / 4;
    //        H->q_offset_trim = hsp->query.offset + MIN(q_length, trim_size);
    //        H->q_end_trim = hsp->query.end - MIN(q_length, trim_size);
    //        H->s_offset_trim = hsp->subject.offset + MIN(s_length, trim_size);
    //        H->s_end_trim = hsp->subject.end - MIN(s_length, trim_size);
    //    }
    // ```
    for slot in 1..=total_number_of_hsps {
        let hsp = *hsp_of(slot);
        let q_length = (hsp.query_end - hsp.query_offset) / 4;
        let s_length = (hsp.subject_end - hsp.subject_offset) / 4;
        let node = &mut nodes[slot];
        node.q_offset_trim = hsp.query_offset + q_length.min(trim_size);
        node.q_end_trim = hsp.query_end - q_length.min(trim_size);
        node.s_offset_trim = hsp.subject_offset + s_length.min(trim_size);
        node.s_end_trim = hsp.subject_end - s_length.min(trim_size);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:451-454
    // ```c
    //    lh_helper_size = MAX(1024,hsp_list->hspcnt+5);
    //    lh_helper = (LinkHelpStruct *)
    //       calloc(lh_helper_size, sizeof(LinkHelpStruct));
    // ```
    let mut lh_helper: Vec<Helper> = vec![Helper::default(); total_number_of_hsps + 5];

    for &(frame_start, frame_count) in &frames {
        link_frame(
            &mut nodes,
            &mut lh_helper,
            frame_start,
            frame_count,
            &hsp_of,
            contexts,
            subject_length,
            LinkConstants {
                window_size,
                trim_size,
                gap_prob,
                gap_decay_rate,
                cutoff,
                ignore_small_gaps,
            },
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:988-1001
    // ```c
    //    if (kTranslatedQuery) {
    // ...
    //    } else {
    //       qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //             s_RevCompareHSPs);
    //       qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
    //             s_FwdCompareHSPs);
    //    }
    // ```
    let mut array: Vec<usize> = (1..=total_number_of_hsps).collect();
    array.sort_by(|&a, &b| rev_compare(hsp_of(a), hsp_of(b)));
    array.sort_by(|&a, &b| fwd_compare(hsp_of(a), hsp_of(b)));

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1006-1080
    // ```c
    //    /* hook up the HSP's. */
    // 	first_hsp = NULL;
    // 	for (index=0, last_hsp=NULL;index<total_number_of_hsps; index++)
    //    {
    // 		H = link_hsp_array[index];
    //
    //       /* If this is not a single piece or the start of a chain, then Skip it. */
    //       if (H->linked_set == TRUE && H->start_of_chain == FALSE)
    // 			continue;
    // ...
    // 		if (H->hsp_link.link[ordering_method] == NULL)
    // 		{
    //          /* Grab the next HSP that is not part of a chain or the start of a chain */
    // ...
    // 		}
    // 		else
    // 		{
    // 			/* The first one has the number of links correct. */
    // 			num_links = H->hsp_link.num[ordering_method];
    // 			link = H->hsp_link.link[ordering_method];
    // 			while (link)
    // 			{
    // 				H->hsp->num = num_links;
    //                 H->xsum = H->hsp_link.xsum[ordering_method];
    // 				H->next = (LinkHSPStruct*) link;
    // ...
    // 			}
    // 			/* Set these for last link in chain. */
    // 			H->hsp->num = num_links;
    // ...
    //    for (index = 0, H = first_hsp; index < hsp_list->hspcnt; index++) {
    //       hsp_list->hsp_array[index] = H->hsp;
    // ```
    // The list NCBI builds visits, in the array's order, each HSP that is not inside a chain,
    // and after a chain's head the rest of its chain.
    let mut out: Vec<LinkedHsp> = Vec::with_capacity(total_number_of_hsps);
    for &slot in &array {
        let head = &nodes[slot];
        if head.linked_set && !head.start_of_chain {
            continue;
        }
        let method = head.ordering_method;
        if head.link[method].is_none() {
            out.push(LinkedHsp {
                index: order[slot - 1],
                evalue: head.evalue,
                num: 1,
            });
            continue;
        }
        let num_links = head.num[method];
        let mut current = Some(slot);
        while let Some(h) = current {
            out.push(LinkedHsp {
                index: order[h - 1],
                evalue: nodes[h].evalue,
                num: num_links,
            });
            current = nodes[h].link[method];
        }
    }
    debug_assert_eq!(out.len(), total_number_of_hsps);
    out
}

/// The values of `s_BlastEvenGapLinkHSPs` that every query frame reads.
#[derive(Clone, Copy)]
struct LinkConstants {
    window_size: i32,
    trim_size: i32,
    gap_prob: f64,
    gap_decay_rate: f64,
    cutoff: [i64; 2],
    ignore_small_gaps: bool,
}

/// The linking of one query frame (`for (frame_index=0; ...)` of `s_BlastEvenGapLinkHSPs`).
#[allow(clippy::too_many_arguments)]
fn link_frame<'a>(
    nodes: &mut [Node],
    lh_helper: &mut [Helper],
    frame_start: usize,
    frame_count: usize,
    hsp_of: &impl Fn(usize) -> &'a LinkHsp,
    contexts: &[LinkContext],
    subject_length_orig: i32,
    k: LinkConstants,
) {
    let LinkConstants {
        window_size,
        trim_size,
        gap_prob,
        gap_decay_rate,
        cutoff,
        ignore_small_gaps,
    } = k;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:554-591
    // ```c
    //       hp_start = s_LinkHSPStructReset(hp_start);
    //       hp_start->next = hp_frame_start[frame_index];
    //       hp_frame_start[frame_index]->prev = hp_start;
    //       number_of_hsps = hp_frame_number[frame_index];
    //       query_context = hp_start->next->hsp->context;
    //       length_adjustment = query_info->contexts[query_context].length_adjustment;
    //       query_length = query_info->contexts[query_context].query_length;
    //       query_length = MAX(query_length - length_adjustment, 1);
    //       subject_length = subject_length_orig; /* in nucleotides even for tblast[nx] */
    // ...
    //       subject_length = MAX(subject_length - length_adjustment, 1);
    //
    //       lh_helper[0].ptr = hp_start;
    //       lh_helper[0].q_off_trim = 0;
    //       lh_helper[0].s_off_trim = 0;
    //       lh_helper[0].maxsum1  = -10000;
    //       lh_helper[0].next_larger  = 0;
    // ...
    //       first_pass=1;    /* do full search */
    //       path_changed=1;
    //       for (H=hp_start->next; H!=NULL; H=H->next)
    //          H->hsp_link.changed=1;
    // ```
    const HP_START: usize = 0;
    nodes[HP_START] = Node::default();
    nodes[HP_START].next = Some(frame_start);
    nodes[frame_start].prev = Some(HP_START);
    let mut number_of_hsps = frame_count;
    let query_context = hsp_of(frame_start).context;
    let ctx = &contexts[query_context];
    let length_adjustment = ctx.length_adjustment;
    let query_length = (ctx.query_length - length_adjustment).max(1);
    let subject_length = (subject_length_orig - length_adjustment).max(1);
    let eff_searchsp = ctx.eff_searchsp;

    lh_helper[0].ptr = HP_START;
    lh_helper[0].q_off_trim = 0;
    lh_helper[0].s_off_trim = 0;
    lh_helper[0].next_larger = 0;

    let mut first_pass = true;
    let mut path_changed = true;
    let mut h = nodes[HP_START].next;
    while let Some(slot) = h {
        nodes[slot].changed = true;
        h = nodes[slot].next;
    }

    let score_of = |slot: usize| hsp_of(slot).score as i64;
    let xsum_term = |slot: usize| {
        let hsp = hsp_of(slot);
        let c = &contexts[hsp.context];
        hsp.score as f64 * c.lambda - c.log_k
    };

    while number_of_hsps > 0 {
        let mut best: [Option<usize>; 2] = [None, None];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:603-656
        // ```c
        //          use_current_max=0;
        //          if (!first_pass){
        //             Int4 max0,max1;
        //             /* Find the current max sums */
        //             if(!ignore_small_gaps){
        //                max0 = -cutoff[0];
        //                max1 = -cutoff[1];
        //                for (H=hp_start->next; H!=NULL; H=H->next) {
        //                   Int4 sum0=H->hsp_link.sum[0];
        //                   Int4 sum1=H->hsp_link.sum[1];
        //                   if(sum0>=max0)
        //                   {
        //                      max0=sum0;
        //                      best[0]=H;
        //                   }
        //                   if(sum1>=max1)
        //                   {
        //                      max1=sum1;
        //                      best[1]=H;
        //                   }
        //                }
        //             } else {
        //                maxscore = -cutoff[1];
        //                for (H=hp_start->next; H!=NULL; H=H->next) {
        //                   Int4  sum=H->hsp_link.sum[1];
        //                   if(sum>=maxscore)
        //                   {
        //                      maxscore=sum;
        //                      best[1]=H;
        //                   }
        //                }
        //             }
        //             if(path_changed==0){
        //                /* No path was changed, use these max sums. */
        //                use_current_max=1;
        //             }
        //             else{
        //                /* If max path hasn't chaged, we can use it */
        //                /* Walk down best, give up if we find a removed item in path */
        //                use_current_max=1;
        //                if(!ignore_small_gaps){
        //                   for (H=best[0]; H!=NULL; H=H->hsp_link.link[0])
        //                      if (H->linked_to==-1000) {use_current_max=0; break;}
        //                }
        //                if(use_current_max)
        //                   for (H=best[1]; H!=NULL; H=H->hsp_link.link[1])
        //                      if (H->linked_to==-1000) {use_current_max=0; break;}
        //
        //             }
        //          }
        // ```
        let mut use_current_max = false;
        if !first_pass {
            if !ignore_small_gaps {
                let mut max0 = -cutoff[SMALL];
                let mut max1 = -cutoff[LARGE];
                let mut h = nodes[HP_START].next;
                while let Some(slot) = h {
                    let sum0 = nodes[slot].sum[SMALL];
                    let sum1 = nodes[slot].sum[LARGE];
                    if sum0 >= max0 {
                        max0 = sum0;
                        best[SMALL] = Some(slot);
                    }
                    if sum1 >= max1 {
                        max1 = sum1;
                        best[LARGE] = Some(slot);
                    }
                    h = nodes[slot].next;
                }
            } else {
                let mut maxscore = -cutoff[LARGE];
                let mut h = nodes[HP_START].next;
                while let Some(slot) = h {
                    let sum = nodes[slot].sum[LARGE];
                    if sum >= maxscore {
                        maxscore = sum;
                        best[LARGE] = Some(slot);
                    }
                    h = nodes[slot].next;
                }
            }
            use_current_max = true;
            if path_changed {
                let removed_in = |nodes: &[Node], start: Option<usize>, method: usize| {
                    let mut h = start;
                    while let Some(slot) = h {
                        if nodes[slot].linked_to == -1000 {
                            return true;
                        }
                        h = nodes[slot].link[method];
                    }
                    false
                };
                if !ignore_small_gaps && removed_in(nodes, best[SMALL], SMALL) {
                    use_current_max = false;
                }
                if use_current_max && removed_in(nodes, best[LARGE], LARGE) {
                    use_current_max = false;
                }
            }
        }

        if !use_current_max {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:663-692
            // ```c
            //             for (H=hp_start,H_index=1; H!=NULL; H=H->next,H_index++) {
            //                Int4 s_frame = H->hsp->subject.frame;
            //                Int4 s_off_t = H->s_offset_trim;
            //                Int4 q_off_t = H->q_offset_trim;
            //                lh_helper[H_index].ptr = H;
            //                lh_helper[H_index].q_off_trim = q_off_t;
            //                lh_helper[H_index].s_off_trim = s_off_t;
            //                for(i=0;i<eOrderingMethods;i++)
            //                   lh_helper[H_index].sum[i] = H->hsp_link.sum[i];
            // ...
            //                /* set next_larger to link back to closest entry with a sum1
            //                   larger than this */
            //                {
            //                   Int4 cur_sum=lh_helper[H_index].sum[1];
            //                   Int4 prev = H_index-1;
            //                   Int4 prev_sum = lh_helper[prev].sum[1];
            //                   while((cur_sum>=prev_sum) && (prev>0)){
            //                      prev=lh_helper[prev].next_larger;
            //                      prev_sum = lh_helper[prev].sum[1];
            //                   }
            //                   lh_helper[H_index].next_larger = prev;
            //                }
            //                H->linked_to = 0;
            //             }
            // ```
            let mut h = Some(HP_START);
            let mut h_index = 1usize;
            while let Some(slot) = h {
                let node = &mut nodes[slot];
                lh_helper[h_index].ptr = slot;
                lh_helper[h_index].q_off_trim = node.q_offset_trim;
                lh_helper[h_index].s_off_trim = node.s_offset_trim;
                lh_helper[h_index].sum = node.sum;
                lh_helper[h_index].next_larger = next_larger(lh_helper, h_index);
                node.linked_to = 0;
                h = node.next;
                h_index += 1;
            }

            if !ignore_small_gaps {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:697-768
                // ```c
                //                index=0;
                //                maxscore = -cutoff[index];
                //                H_index = 2;
                //                for (H=hp_start->next; H!=NULL; H=H->next,H_index++)
                //                {
                // ...
                //                   if (H->hsp->score > cutoff[index]) {
                //                      Int4 H_query_etrim = H->q_end_trim;
                //                      Int4 H_sub_etrim = H->s_end_trim;
                //                      Int4 H_q_et_gap = H_query_etrim+window_size;
                //                      Int4 H_s_et_gap = H_sub_etrim+window_size;
                // ...
                //                      for (H2_index=H_index-1; H2_index>1; H2_index=H2_index-1)
                //                      {
                // ...
                //                         q_off_t = lh_helper[H2_index].q_off_trim;
                //                         s_off_t = lh_helper[H2_index].s_off_trim;
                //
                //                         /* combine tests to reduce mispredicts -cfj */
                //                         b1 = q_off_t <= H_query_etrim;
                //                         b2 = s_off_t <= H_sub_etrim;
                //                         sum = lh_helper[H2_index].sum[index];
                //
                //
                //                         b4 = ( q_off_t > H_q_et_gap ) ;
                //                         b5 = ( s_off_t > H_s_et_gap ) ;
                // ...
                //                         if(q_off_t > (H_q_et_gap+trim_size))
                //                            break;
                //
                //                         if (b1|b2|b5|b4) continue;
                //
                //                         if (sum>H_hsp_sum)
                //                         {
                //                            H2=lh_helper[H2_index].ptr;
                //                            H_hsp_num=H2->hsp_link.num[index];
                //                            H_hsp_sum=H2->hsp_link.sum[index];
                //                            H_hsp_xsum=H2->hsp_link.xsum[index];
                //                            H_hsp_link=H2;
                //                         }
                //                      } /* end for H2... */
                //                   }
                //                   {
                //                      Int4 score=H->hsp->score;
                //                      double new_xsum =
                //                        H_hsp_xsum + score*kbp[H->hsp->context]->Lambda -
                //                        kbp[H->hsp->context]->logK;
                //                      Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
                //
                //                      H->hsp_link.sum[index] = new_sum;
                //                      H->hsp_link.num[index] = H_hsp_num+1;
                //                      H->hsp_link.link[index] = H_hsp_link;
                //                      lh_helper[H_index].sum[index] = new_sum;
                //                      if (new_sum >= maxscore)
                //                      {
                //                         maxscore=new_sum;
                //                         best[index]=H;
                //                      }
                //                      H->hsp_link.xsum[index] = new_xsum;
                //                      if(H_hsp_link)
                //                         ((LinkHSPStruct*)H_hsp_link)->linked_to++;
                //                   }
                //                } /* end for H=... */
                // ```
                let index = SMALL;
                let mut maxscore = -cutoff[index];
                let mut h_index = 2usize;
                let mut h = nodes[HP_START].next;
                while let Some(slot) = h {
                    let mut h_hsp_num = 0i32;
                    let mut h_hsp_sum = 0i64;
                    let mut h_hsp_xsum = 0.0f64;
                    let mut h_hsp_link: Option<usize> = None;
                    if score_of(slot) > cutoff[index] {
                        let h_query_etrim = nodes[slot].q_end_trim;
                        let h_sub_etrim = nodes[slot].s_end_trim;
                        let h_q_et_gap = h_query_etrim + window_size;
                        let h_s_et_gap = h_sub_etrim + window_size;
                        let mut h2_index = h_index - 1;
                        while h2_index > 1 {
                            let helper = lh_helper[h2_index];
                            let q_off_t = helper.q_off_trim;
                            let s_off_t = helper.s_off_trim;
                            let b1 = q_off_t <= h_query_etrim;
                            let b2 = s_off_t <= h_sub_etrim;
                            let sum = helper.sum[index];
                            let b4 = q_off_t > h_q_et_gap;
                            let b5 = s_off_t > h_s_et_gap;
                            if q_off_t > h_q_et_gap + trim_size {
                                break;
                            }
                            if !(b1 | b2 | b5 | b4) && sum > h_hsp_sum {
                                let h2 = &nodes[helper.ptr];
                                h_hsp_num = h2.num[index];
                                h_hsp_sum = h2.sum[index];
                                h_hsp_xsum = h2.xsum[index];
                                h_hsp_link = Some(helper.ptr);
                            }
                            h2_index -= 1;
                        }
                    }
                    let new_xsum = h_hsp_xsum + xsum_term(slot);
                    let new_sum = h_hsp_sum + (score_of(slot) - cutoff[index]);
                    let node = &mut nodes[slot];
                    node.sum[index] = new_sum;
                    node.num[index] = h_hsp_num + 1;
                    node.link[index] = h_hsp_link;
                    lh_helper[h_index].sum[index] = new_sum;
                    if new_sum >= maxscore {
                        maxscore = new_sum;
                        best[index] = Some(slot);
                    }
                    node.xsum[index] = new_xsum;
                    let next = node.next;
                    if let Some(link) = h_hsp_link {
                        nodes[link].linked_to += 1;
                    }
                    h = next;
                    h_index += 1;
                }
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:770-889
            // ```c
            //             index=1;
            //             maxscore = -cutoff[index];
            //             H_index = 2;
            //             for (H=hp_start->next; H!=NULL; H=H->next,H_index++)
            //             {
            // ...
            //                H->hsp_link.changed=1;
            //                H2 = H->hsp_link.link[index];
            //                if ((!first_pass) && ((H2==0) || (H2->hsp_link.changed==0)))
            //                {
            // ...
            //                   if(H2){
            //                      H_hsp_num=H2->hsp_link.num[index];
            //                      H_hsp_sum=H2->hsp_link.sum[index];
            //                      H_hsp_xsum=H2->hsp_link.xsum[index];
            //                   }
            //                   H_hsp_link=H2;
            //                   H->hsp_link.changed=0;
            //                } else if (H->hsp->score > cutoff[index]) {
            //                   Int4 H_query_etrim = H->q_end_trim;
            //                   Int4 H_sub_etrim = H->s_end_trim;
            // ...
            //                   if(!first_pass&&H2&&H2->linked_to>=0){
            //                      if(1){
            //                         /* We set this to less than the real value to keep the
            //                            original ordering in case of ties. */
            //                         H_hsp_sum=H2->hsp_link.sum[index]-1;
            // ...
            //                   for (H2_index=H_index-1; H2_index>1;)
            //                   {
            //                      Int4 b0,b1,b2;
            //                      Int4 q_off_t,s_off_t,sum,next_larger;
            //                      LinkHelpStruct * H2_helper=&lh_helper[H2_index];
            //                      sum = H2_helper->sum[index];
            //                      next_larger = H2_helper->next_larger;
            //
            //                      s_off_t = H2_helper->s_off_trim;
            //                      q_off_t = H2_helper->q_off_trim;
            //
            //                      b0 = sum <= H_hsp_sum;
            //
            //                      /* Compute the next H2_index */
            //                      H2_index--;
            //                      if(b0){	 /* If this sum is too small to beat H_hsp_sum, advance to a larger sum */
            //                         H2_index=next_larger;
            //                      }
            //
            //                      /* combine tests to reduce mispredicts -cfj */
            //                      b1 = q_off_t <= H_query_etrim;
            //                      b2 = s_off_t <= H_sub_etrim;
            // ...
            //                      if (!(b0|b1|b2) )
            //                      {
            //                         H2 = H2_helper->ptr;
            //
            //                         H_hsp_num=H2->hsp_link.num[index];
            //                         H_hsp_sum=H2->hsp_link.sum[index];
            //                         H_hsp_xsum=H2->hsp_link.xsum[index];
            //                         H_hsp_link=H2;
            //                      }
            //                   } /* end for H2_index... */
            //                } /* end if(H->score>cuttof[]) */
            //                {
            //                   Int4 score=H->hsp->score;
            //                   double new_xsum =
            //                      H_hsp_xsum + score*kbp[H->hsp->context]->Lambda -
            //                      kbp[H->hsp->context]->logK;
            //                   Int4 new_sum = H_hsp_sum + (score - cutoff[index]);
            //
            //                   H->hsp_link.sum[index] = new_sum;
            //                   H->hsp_link.num[index] = H_hsp_num+1;
            //                   H->hsp_link.link[index] = H_hsp_link;
            //                   lh_helper[H_index].sum[index] = new_sum;
            //                   lh_helper[H_index].maxsum1 = MAX(lh_helper[H_index-1].maxsum1, new_sum);
            //                   /* Update this entry's 'next_larger' field */
            //                   {
            // ...
            //                   }
            //
            //                   if (new_sum >= maxscore)
            //                   {
            //                      maxscore=new_sum;
            //                      best[index]=H;
            //                   }
            //                   H->hsp_link.xsum[index] = new_xsum;
            //                   if(H_hsp_link)
            //                      ((LinkHSPStruct*)H_hsp_link)->linked_to++;
            //                }
            //             }
            //             path_changed=0;
            //             first_pass=0;
            // ```
            let index = LARGE;
            let mut maxscore = -cutoff[index];
            let mut h_index = 2usize;
            let mut h = nodes[HP_START].next;
            while let Some(slot) = h {
                let mut h_hsp_num = 0i32;
                let mut h_hsp_sum = 0i64;
                let mut h_hsp_xsum = 0.0f64;
                let mut h_hsp_link: Option<usize> = None;
                nodes[slot].changed = true;
                let h2 = nodes[slot].link[index];
                if !first_pass && h2.is_none_or(|h2| !nodes[h2].changed) {
                    if let Some(h2) = h2 {
                        h_hsp_num = nodes[h2].num[index];
                        h_hsp_sum = nodes[h2].sum[index];
                        h_hsp_xsum = nodes[h2].xsum[index];
                    }
                    h_hsp_link = h2;
                    nodes[slot].changed = false;
                } else if score_of(slot) > cutoff[index] {
                    let h_query_etrim = nodes[slot].q_end_trim;
                    let h_sub_etrim = nodes[slot].s_end_trim;
                    if !first_pass {
                        if let Some(h2) = h2.filter(|&h2| nodes[h2].linked_to >= 0) {
                            h_hsp_sum = nodes[h2].sum[index] - 1;
                        }
                    }
                    let mut h2_index = h_index - 1;
                    while h2_index > 1 {
                        let helper = lh_helper[h2_index];
                        let sum = helper.sum[index];
                        let b0 = sum <= h_hsp_sum;
                        h2_index -= 1;
                        if b0 {
                            h2_index = helper.next_larger;
                        }
                        let b1 = helper.q_off_trim <= h_query_etrim;
                        let b2 = helper.s_off_trim <= h_sub_etrim;
                        if !(b0 | b1 | b2) {
                            let h2 = &nodes[helper.ptr];
                            h_hsp_num = h2.num[index];
                            h_hsp_sum = h2.sum[index];
                            h_hsp_xsum = h2.xsum[index];
                            h_hsp_link = Some(helper.ptr);
                        }
                    }
                }
                let new_xsum = h_hsp_xsum + xsum_term(slot);
                let new_sum = h_hsp_sum + (score_of(slot) - cutoff[index]);
                let node = &mut nodes[slot];
                node.sum[index] = new_sum;
                node.num[index] = h_hsp_num + 1;
                node.link[index] = h_hsp_link;
                lh_helper[h_index].sum[index] = new_sum;
                lh_helper[h_index].next_larger = next_larger(lh_helper, h_index);
                if new_sum >= maxscore {
                    maxscore = new_sum;
                    best[index] = Some(slot);
                }
                node.xsum[index] = new_xsum;
                let next = node.next;
                if let Some(link) = h_hsp_link {
                    nodes[link].linked_to += 1;
                }
                h = next;
                h_index += 1;
            }
            path_changed = false;
            first_pass = false;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:891-946
        // ```c
        //          if (!ignore_small_gaps)
        //          {
        // ...
        //             best[0]->hsp_link.sum[0] +=
        //                (best[0]->hsp_link.num[0])*cutoff[0];
        //
        //             prob[0] = BLAST_SmallGapSumE(window_size,
        //                          best[0]->hsp_link.num[0], best[0]->hsp_link.xsum[0],
        //                          query_length, subject_length,
        //                          query_info->contexts[query_context].eff_searchsp,
        //                          BLAST_GapDecayDivisor(gap_decay_rate,
        //                                               best[0]->hsp_link.num[0]) );
        //
        //             /* Adjust the e-value because we are performing multiple tests */
        //             if( best[0]->hsp_link.num[0] > 1 ) {
        //               if( gap_prob == 0 || (prob[0] /= gap_prob) > INT4_MAX ) {
        //                 prob[0] = INT4_MAX;
        //               }
        //             }
        //
        //             prob[1] = BLAST_LargeGapSumE(best[1]->hsp_link.num[1],
        //                          best[1]->hsp_link.xsum[1],
        //                          query_length, subject_length,
        //                          query_info->contexts[query_context].eff_searchsp,
        //                          BLAST_GapDecayDivisor(gap_decay_rate,
        //                                               best[1]->hsp_link.num[1]));
        //
        //             if( best[1]->hsp_link.num[1] > 1 ) {
        //               if( 1 - gap_prob == 0 || (prob[1] /= 1 - gap_prob) > INT4_MAX ) {
        //                 prob[1] = INT4_MAX;
        //               }
        //             }
        //             ordering_method =
        //                prob[0]<=prob[1] ? eLinkSmallGaps : eLinkLargeGaps;
        //          }
        //          else
        //          {
        //             /* We only consider the case of big gaps. */
        //             best[1]->hsp_link.sum[1] +=
        //                (best[1]->hsp_link.num[1])*cutoff[1];
        //
        //             prob[1] = BLAST_LargeGapSumE(
        // ...
        //             ordering_method = eLinkLargeGaps;
        //          }
        // ```
        let large_gap_e = |nodes: &[Node], best: usize| {
            let num = nodes[best].num[LARGE];
            ncbi_large_gap_sum_e(
                num as i16,
                nodes[best].xsum[LARGE],
                query_length,
                subject_length,
                eff_searchsp,
                gap_decay_divisor(gap_decay_rate, num as usize),
            )
        };
        let mut prob = [0.0f64; 2];
        let ordering_method;
        if !ignore_small_gaps {
            let b0 = best[SMALL].expect("a frame with HSPs has a best small-gap set");
            let b1 = best[LARGE].expect("a frame with HSPs has a best large-gap set");
            nodes[b0].sum[SMALL] += nodes[b0].num[SMALL] as i64 * cutoff[SMALL];
            let num0 = nodes[b0].num[SMALL];
            prob[SMALL] = small_gap_sum_e(
                window_size,
                num0 as i16,
                nodes[b0].xsum[SMALL],
                query_length,
                subject_length,
                eff_searchsp,
                gap_decay_divisor(gap_decay_rate, num0 as usize),
            );
            if num0 > 1 {
                if gap_prob == 0.0 {
                    prob[SMALL] = i32::MAX as f64;
                } else {
                    prob[SMALL] /= gap_prob;
                    if prob[SMALL] > i32::MAX as f64 {
                        prob[SMALL] = i32::MAX as f64;
                    }
                }
            }
            prob[LARGE] = large_gap_e(nodes, b1);
            if nodes[b1].num[LARGE] > 1 {
                if 1.0 - gap_prob == 0.0 {
                    prob[LARGE] = i32::MAX as f64;
                } else {
                    prob[LARGE] /= 1.0 - gap_prob;
                    if prob[LARGE] > i32::MAX as f64 {
                        prob[LARGE] = i32::MAX as f64;
                    }
                }
            }
            ordering_method = if prob[SMALL] <= prob[LARGE] {
                SMALL
            } else {
                LARGE
            };
        } else {
            let b1 = best[LARGE].expect("a frame with HSPs has a best large-gap set");
            nodes[b1].sum[LARGE] += nodes[b1].num[LARGE] as i64 * cutoff[LARGE];
            prob[LARGE] = large_gap_e(nodes, b1);
            ordering_method = LARGE;
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:948-980
        // ```c
        //          best[ordering_method]->start_of_chain = TRUE;
        //          best[ordering_method]->hsp->evalue    = prob[ordering_method];
        //
        //          /* remove the links that have been ordered already. */
        //          if (best[ordering_method]->hsp_link.link[ordering_method])
        //             linked_set = TRUE;
        //          else
        //             linked_set = FALSE;
        //
        //          if (best[ordering_method]->linked_to>0) path_changed=1;
        //          for (H=best[ordering_method]; H!=NULL;
        //               H=H->hsp_link.link[ordering_method])
        //          {
        //             if (H->linked_to>1) path_changed=1;
        //             H->linked_to=-1000;
        //             H->hsp_link.changed=1;
        //             /* record whether this is part of a linked set. */
        //             H->linked_set = linked_set;
        //             H->ordering_method = ordering_method;
        //             H->hsp->evalue = prob[ordering_method];
        //             if (H->next)
        //                (H->next)->prev=H->prev;
        //             if (H->prev)
        //                (H->prev)->next=H->next;
        //             number_of_hsps--;
        //          }
        // ```
        let chosen = best[ordering_method].expect("the chosen ordering has a best set");
        nodes[chosen].start_of_chain = true;
        nodes[chosen].evalue = prob[ordering_method];
        let linked_set = nodes[chosen].link[ordering_method].is_some();
        if nodes[chosen].linked_to > 0 {
            path_changed = true;
        }
        let mut h = Some(chosen);
        while let Some(slot) = h {
            if nodes[slot].linked_to > 1 {
                path_changed = true;
            }
            let node = &mut nodes[slot];
            node.linked_to = -1000;
            node.changed = true;
            node.linked_set = linked_set;
            node.ordering_method = ordering_method;
            node.evalue = prob[ordering_method];
            let (prev, next) = (node.prev, node.next);
            if let Some(next) = next {
                nodes[next].prev = prev;
            }
            if let Some(prev) = prev {
                nodes[prev].next = next;
            }
            number_of_hsps -= 1;
            h = nodes[slot].link[ordering_method];
        }
    }
}

/// `next_larger` of `lh_helper[h_index]`: the closest earlier entry whose large-gap sum is
/// larger.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:679-688 (the same at 864-873)
/// ```c
///                {
///                   Int4 cur_sum=lh_helper[H_index].sum[1];
///                   Int4 prev = H_index-1;
///                   Int4 prev_sum = lh_helper[prev].sum[1];
///                   while((cur_sum>=prev_sum) && (prev>0)){
///                      prev=lh_helper[prev].next_larger;
///                      prev_sum = lh_helper[prev].sum[1];
///                   }
///                   lh_helper[H_index].next_larger = prev;
///                }
/// ```
fn next_larger(lh_helper: &[Helper], h_index: usize) -> usize {
    let cur_sum = lh_helper[h_index].sum[LARGE];
    let mut prev = h_index - 1;
    let mut prev_sum = lh_helper[prev].sum[LARGE];
    while cur_sum >= prev_sum && prev > 0 {
        prev = lh_helper[prev].next_larger;
        prev_sum = lh_helper[prev].sum[LARGE];
    }
    prev
}

#[cfg(test)]
mod tests {
    use super::*;

    fn context() -> LinkContext {
        // BLOSUM62 11/1 gapped Lambda and K (blast_stat.c), a query of 300 letters.
        LinkContext {
            lambda: 0.267,
            log_k: 0.041f64.ln(),
            query_length: 300,
            length_adjustment: 40,
            eff_searchsp: 260 * 960,
        }
    }

    fn hsp(q: (i32, i32), s: (i32, i32), score: i32) -> LinkHsp {
        LinkHsp {
            context: 0,
            query_offset: q.0,
            query_end: q.1,
            subject_offset: s.0,
            subject_end: s.1,
            score,
        }
    }

    #[test]
    fn a_lone_hsp_gets_the_e_value_of_one_segment() {
        let hsps = [hsp((10, 60), (20, 70), 80)];
        let ctx = context();
        let out = even_gap_link_hsps(&hsps, &[ctx], 1000, &EvenGapLinkParams::traceback());
        assert_eq!(out.len(), 1);
        assert_eq!((out[0].index, out[0].num), (0, 1));
        let xsum = 80.0 * ctx.lambda - ctx.log_k;
        let expected = (ctx.eff_searchsp as f64) * (-xsum).exp() / 0.9;
        assert_eq!(out[0].evalue, expected);
    }

    #[test]
    fn two_hsps_in_order_on_both_sequences_link_into_one_set() {
        // The first HSP ends before the second starts on the query and on the subject.
        let hsps = [hsp((200, 260), (500, 560), 60), hsp((10, 60), (20, 70), 80)];
        let out = even_gap_link_hsps(&hsps, &[context()], 1000, &EvenGapLinkParams::traceback());
        // The chain starts at the upstream HSP (the head) and follows its link downstream.
        assert_eq!(out.iter().map(|l| l.index).collect::<Vec<_>>(), vec![1, 0]);
        assert!(out.iter().all(|l| l.num == 2));
        assert_eq!(out[0].evalue, out[1].evalue);
    }

    #[test]
    fn hsps_in_opposite_orders_on_the_sequences_stay_apart() {
        let hsps = [hsp((200, 260), (20, 70), 60), hsp((10, 60), (500, 560), 80)];
        let out = even_gap_link_hsps(&hsps, &[context()], 1000, &EvenGapLinkParams::traceback());
        assert!(out.iter().all(|l| l.num == 1));
        // The list is in the forward order of the query starts.
        assert_eq!(out.iter().map(|l| l.index).collect::<Vec<_>>(), vec![1, 0]);
    }

    #[test]
    fn the_preliminary_cutoff_keeps_hsps_at_or_below_it_out_of_small_gap_sets() {
        // Close HSPs (a gap of 10 letters on both sequences) with scores at the small-gap cutoff:
        // NCBI walks for a predecessor only for a score above the cutoff.
        let hsps = [hsp((70, 120), (80, 130), 41), hsp((10, 60), (20, 70), 41)];
        let out = even_gap_link_hsps(
            &hsps,
            &[context()],
            1000,
            &EvenGapLinkParams::preliminary(41),
        );
        assert!(out.iter().all(|l| l.num == 1), "{out:?}");
    }
}
